/// Full pipeline benchmark: QZ default vs QZ+new strategies vs SPRING
/// Measures compression ratio, speed, and peak memory for 500K reads.
///
/// New strategies:
///   - Sequences: syncmer template hybrid (k=31)
///   - Qualities: mean-quality reorder + fqzcomp strat=0

use std::io::BufRead;
use std::time::Instant;

use qz_lib::compression::bsc;
use qz_lib::compression::dna_utils::{kmer_to_hash, reverse_complement_hash};
use qz_lib::compression::fqzcomp;
use qz_lib::compression::template::{compress_sequences_template_hybrid, TemplateParams};
use rayon::prelude::*;

fn main() {
    let fastq_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/ERR3239334_1.500k.fastq".to_string());

    eprintln!("Reading FASTQ: {}", fastq_path);
    let file = std::fs::File::open(&fastq_path).expect("Cannot open FASTQ file");
    let reader = std::io::BufReader::new(file);

    let mut headers: Vec<String> = Vec::new();
    let mut sequences: Vec<String> = Vec::new();
    let mut qualities: Vec<String> = Vec::new();
    let mut line_num = 0u64;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        match line_num % 4 {
            0 => headers.push(line.trim_end().to_string()),
            1 => sequences.push(line.trim_end().to_string()),
            3 => qualities.push(line.trim_end().to_string()),
            _ => {}
        }
        line_num += 1;
    }

    let num_reads = sequences.len();
    let total_header_bytes: usize = headers.iter().map(|h| h.len()).sum();
    let total_seq_bytes: usize = sequences.iter().map(|s| s.len()).sum();
    let total_qual_bytes: usize = qualities.iter().map(|q| q.len()).sum();
    let read_len = sequences.first().map(|s| s.len()).unwrap_or(0);
    // Original FASTQ size: @header\nseq\n+\nqual\n per read
    let original_size = total_header_bytes + total_seq_bytes + total_qual_bytes + num_reads * 4; // 4 newlines + @ + +

    eprintln!(
        "Loaded {} reads, read_len={}, total={:.1} MB",
        num_reads,
        read_len,
        original_size as f64 / (1024.0 * 1024.0),
    );
    eprintln!(
        "  Headers: {:.1} MB, Sequences: {:.1} MB, Qualities: {:.1} MB",
        total_header_bytes as f64 / (1024.0 * 1024.0),
        total_seq_bytes as f64 / (1024.0 * 1024.0),
        total_qual_bytes as f64 / (1024.0 * 1024.0),
    );

    println!();
    println!("=== STREAM-BY-STREAM COMPARISON ===");
    println!();
    println!(
        "{:<45} {:>12} {:>8} {:>10}",
        "Stream / Method", "Size", "Ratio", "Time"
    );
    println!("{}", "-".repeat(78));

    // ── HEADERS ─────────────────────────────────────────────
    println!("HEADERS ({:.1} MB raw):", total_header_bytes as f64 / (1024.0 * 1024.0));

    let header_bytes: Vec<u8> = headers.iter().flat_map(|h| h.as_bytes()).copied().collect();
    let t = Instant::now();
    let headers_bsc = bsc::compress_parallel_adaptive(&header_bytes).unwrap();
    let headers_time = t.elapsed();
    println!(
        "  {:<43} {:>10} B {:>7.2}x {:>8.1}ms",
        "BSC adaptive (current default)",
        headers_bsc.len(),
        total_header_bytes as f64 / headers_bsc.len() as f64,
        headers_time.as_secs_f64() * 1000.0,
    );

    // ── SEQUENCES ───────────────────────────────────────────
    println!("SEQUENCES ({:.1} MB raw):", total_seq_bytes as f64 / (1024.0 * 1024.0));

    // Baseline: raw + BSC
    let seq_bytes: Vec<u8> = sequences.iter().flat_map(|s| s.as_bytes()).copied().collect();
    let t = Instant::now();
    let seq_bsc = bsc::compress_parallel_adaptive(&seq_bytes).unwrap();
    let seq_bsc_time = t.elapsed();
    println!(
        "  {:<43} {:>10} B {:>7.2}x {:>8.1}ms",
        "BSC adaptive (current default)",
        seq_bsc.len(),
        total_seq_bytes as f64 / seq_bsc.len() as f64,
        seq_bsc_time.as_secs_f64() * 1000.0,
    );

    // New: syncmer template hybrid (skip for large datasets)
    let (seq_template_len, seq_template_time) = if num_reads <= 1_000_000 {
        let params = TemplateParams {
            k: 31,
            fast: true,
            use_syncmers: true,
            ..TemplateParams::default()
        };
        let t = Instant::now();
        let seq_template = compress_sequences_template_hybrid(&sequences, &params).unwrap();
        let elapsed = t.elapsed();
        println!(
            "  {:<43} {:>10} B {:>7.2}x {:>8.1}ms",
            "Syncmer template hybrid (k=31)",
            seq_template.len(),
            total_seq_bytes as f64 / seq_template.len() as f64,
            elapsed.as_secs_f64() * 1000.0,
        );
        (seq_template.len(), elapsed)
    } else {
        println!(
            "  {:<43} {:>10}   (skipped, >1M reads)",
            "Syncmer template hybrid (k=31)", "—",
        );
        (0, std::time::Duration::ZERO)
    };

    // New: syncmer reorder (fast, no graph)
    let t = Instant::now();
    let anchor_k = 31usize;
    let s = anchor_k - 3;
    let t_end = anchor_k - s;
    let sort_keys: Vec<u64> = sequences
        .par_iter()
        .map(|seq| {
            let seq_bytes = seq.as_bytes();
            if seq_bytes.len() < anchor_k { return u64::MAX; }
            let positions = syncmers::find_syncmers_pos(anchor_k, s, &[0, t_end], seq_bytes);
            let mut min_hash = u64::MAX;
            for pos in positions {
                if pos + anchor_k > seq_bytes.len() { continue; }
                let kmer = &seq_bytes[pos..pos + anchor_k];
                if let Some(fwd) = kmer_to_hash(kmer) {
                    let rc = reverse_complement_hash(fwd, anchor_k);
                    let canon = fwd.min(rc);
                    if canon < min_hash { min_hash = canon; }
                }
            }
            min_hash
        })
        .collect();
    let mut seq_perm: Vec<u32> = (0..num_reads as u32).collect();
    seq_perm.sort_by_key(|&i| sort_keys[i as usize]);
    let unique_keys = {
        let mut sk = sort_keys.clone();
        sk.sort_unstable();
        sk.dedup();
        sk.len()
    };
    let mut reordered_seq = Vec::with_capacity(total_seq_bytes);
    for &idx in &seq_perm {
        reordered_seq.extend_from_slice(sequences[idx as usize].as_bytes());
    }
    let seq_reorder_bsc = bsc::compress_parallel_adaptive(&reordered_seq).unwrap();
    // Permutation cost (delta-zigzag-varint + BSC)
    let mut seq_perm_delta: Vec<u8> = Vec::with_capacity(num_reads * 3);
    let mut prev = 0i64;
    for &idx in &seq_perm {
        let delta = idx as i64 - prev;
        prev = idx as i64;
        let zz = ((delta << 1) ^ (delta >> 63)) as u64;
        let mut v = zz;
        while v >= 0x80 {
            seq_perm_delta.push((v as u8) | 0x80);
            v >>= 7;
        }
        seq_perm_delta.push(v as u8);
    }
    let seq_perm_bsc = bsc::compress_parallel_adaptive(&seq_perm_delta).unwrap();
    let seq_reorder_time = t.elapsed();
    let seq_reorder_total = seq_reorder_bsc.len() + seq_perm_bsc.len();
    println!(
        "  {:<43} {:>10} B {:>7.2}x {:>8.1}ms",
        "Syncmer reorder + BSC (seq+perm)",
        seq_reorder_total,
        total_seq_bytes as f64 / seq_reorder_total as f64,
        seq_reorder_time.as_secs_f64() * 1000.0,
    );
    println!(
        "    {:<41} {:>10} B  (sequences BSC'd)",
        "", seq_reorder_bsc.len(),
    );
    println!(
        "    {:<41} {:>10} B  (permutation BSC'd)",
        "", seq_perm_bsc.len(),
    );
    println!(
        "    {:<41} {:>10}    ({:.1}% unique syncmer keys)",
        "", "", 100.0 * unique_keys as f64 / num_reads as f64,
    );

    // Syncmer-reordered qualities using same permutation (free perm!)
    let syncmer_reordered_quals: Vec<&[u8]> = seq_perm.iter().map(|&i| qualities[i as usize].as_bytes()).collect();
    let t = Instant::now();
    let qual_syncmer_fqz = fqzcomp::compress(&syncmer_reordered_quals, 0).unwrap();
    let qual_syncmer_fqz_time = t.elapsed();
    // No extra perm cost — shared with sequences

    // ── QUALITIES ───────────────────────────────────────────
    println!("QUALITIES ({:.1} MB raw):", total_qual_bytes as f64 / (1024.0 * 1024.0));

    // Baseline: raw + BSC
    let qual_bytes: Vec<u8> = qualities.iter().flat_map(|q| q.as_bytes()).copied().collect();
    let t = Instant::now();
    let qual_bsc = bsc::compress_parallel_adaptive(&qual_bytes).unwrap();
    let qual_bsc_time = t.elapsed();
    println!(
        "  {:<43} {:>10} B {:>7.2}x {:>8.1}ms",
        "BSC adaptive (current default)",
        qual_bsc.len(),
        total_qual_bytes as f64 / qual_bsc.len() as f64,
        qual_bsc_time.as_secs_f64() * 1000.0,
    );

    // New: mean-quality reorder + fqzcomp
    let t = Instant::now();
    let mut perm_mean: Vec<u32> = (0..num_reads as u32).collect();
    perm_mean.sort_by_key(|&i| {
        let q = qualities[i as usize].as_bytes();
        let sum: u64 = q.iter().map(|&b| b as u64).sum();
        sum / q.len().max(1) as u64
    });
    let reordered_quals: Vec<&[u8]> = perm_mean.iter().map(|&i| qualities[i as usize].as_bytes()).collect();
    let qual_fqz = fqzcomp::compress(&reordered_quals, 0).unwrap();
    // Store mean-quality sort key (1 byte per read)
    let mean_quals: Vec<u8> = qualities.iter().map(|q| {
        let bytes = q.as_bytes();
        let sum: u64 = bytes.iter().map(|&b| b as u64).sum();
        (sum / bytes.len().max(1) as u64) as u8
    }).collect();
    let mean_key_bsc = bsc::compress_parallel_adaptive(&mean_quals).unwrap();
    let qual_new_time = t.elapsed();
    let qual_new_total = qual_fqz.len() + mean_key_bsc.len();
    println!(
        "  {:<43} {:>10} B {:>7.2}x {:>8.1}ms",
        "Mean-sort + fqzcomp (fqz+key)",
        qual_new_total,
        total_qual_bytes as f64 / qual_new_total as f64,
        qual_new_time.as_secs_f64() * 1000.0,
    );
    println!(
        "    {:<41} {:>10} B  (fqzcomp payload)",
        "", qual_fqz.len(),
    );
    println!(
        "    {:<41} {:>10} B  (mean-qual sort key)",
        "", mean_key_bsc.len(),
    );
    println!(
        "  {:<43} {:>10} B {:>7.2}x {:>8.1}ms",
        "Syncmer-sort + fqzcomp (shared perm)",
        qual_syncmer_fqz.len(),
        total_qual_bytes as f64 / qual_syncmer_fqz.len() as f64,
        qual_syncmer_fqz_time.as_secs_f64() * 1000.0,
    );

    // ── TOTALS ──────────────────────────────────────────────
    println!();
    println!("=== TOTAL COMPRESSED SIZE ===");
    println!();

    let metadata_est = 60; // archive format overhead

    let qz_default_total = headers_bsc.len() + seq_bsc.len() + qual_bsc.len() + metadata_est;
    let qz_default_time = headers_time.max(seq_bsc_time + qual_bsc_time); // headers parallel with seq+qual

    let qz_template_total = if seq_template_len > 0 {
        headers_bsc.len() + seq_template_len + qual_new_total + metadata_est
    } else { 0 };
    let qz_template_time = if seq_template_len > 0 {
        headers_time.max(seq_template_time + qual_new_time)
    } else { std::time::Duration::ZERO };

    let qz_reorder_total = headers_bsc.len() + seq_reorder_total + qual_new_total + metadata_est;
    let qz_reorder_time_seq = seq_reorder_time + qual_new_time;
    let qz_reorder_time = headers_time.max(qz_reorder_time_seq);

    // Shared perm: syncmer reorder seq+qual, ONE permutation
    let qz_shared_total = headers_bsc.len() + seq_reorder_bsc.len() + qual_syncmer_fqz.len() + seq_perm_bsc.len() + metadata_est;
    let qz_shared_time = headers_time.max(seq_reorder_time + qual_syncmer_fqz_time);

    println!(
        "{:<45} {:>12} {:>8}   {:>10}",
        "Tool", "Total", "Ratio", "Time"
    );
    println!("{}", "-".repeat(78));
    println!(
        "{:<45} {:>10} B {:>7.2}x   {:>8.1}ms",
        "QZ default (BSC all streams)",
        qz_default_total,
        original_size as f64 / qz_default_total as f64,
        qz_default_time.as_secs_f64() * 1000.0,
    );
    println!(
        "  Headers: {} B, Seq: {} B, Qual: {} B",
        headers_bsc.len(), seq_bsc.len(), qual_bsc.len(),
    );

    println!(
        "{:<45} {:>10} B {:>7.2}x   {:>8.1}ms",
        "QZ fast (reorder seq + fqzcomp qual)",
        qz_reorder_total,
        original_size as f64 / qz_reorder_total as f64,
        qz_reorder_time.as_secs_f64() * 1000.0,
    );
    println!(
        "  Headers: {} B, Seq: {} B, Qual: {} B",
        headers_bsc.len(), seq_reorder_total, qual_new_total,
    );
    let improvement_reorder = 100.0 * (1.0 - qz_reorder_total as f64 / qz_default_total as f64);
    println!(
        "  Improvement over default: {:.1}% smaller ({:+} B)",
        improvement_reorder,
        qz_default_total as i64 - qz_reorder_total as i64,
    );

    println!(
        "{:<45} {:>10} B {:>7.2}x   {:>8.1}ms",
        "QZ shared-perm (syncmer seq+qual, 1 perm)",
        qz_shared_total,
        original_size as f64 / qz_shared_total as f64,
        qz_shared_time.as_secs_f64() * 1000.0,
    );
    println!(
        "  Hdr: {} B, Seq: {} B, Qual: {} B, Perm: {} B",
        headers_bsc.len(), seq_reorder_bsc.len(), qual_syncmer_fqz.len(), seq_perm_bsc.len(),
    );
    let improvement_shared = 100.0 * (1.0 - qz_shared_total as f64 / qz_default_total as f64);
    println!(
        "  Improvement over default: {:.1}% smaller ({:+} B)",
        improvement_shared,
        qz_default_total as i64 - qz_shared_total as i64,
    );

    if qz_template_total > 0 {
        println!(
            "{:<45} {:>10} B {:>7.2}x   {:>8.1}ms",
            "QZ slow (template seq + fqzcomp qual)",
            qz_template_total,
            original_size as f64 / qz_template_total as f64,
            qz_template_time.as_secs_f64() * 1000.0,
        );
        println!(
            "  Headers: {} B, Seq: {} B, Qual: {} B",
            headers_bsc.len(), seq_template_len, qual_new_total,
        );
        let improvement_template = 100.0 * (1.0 - qz_template_total as f64 / qz_default_total as f64);
        println!(
            "  Improvement over default: {:.1}% smaller ({:+} B)",
            improvement_template,
            qz_default_total as i64 - qz_template_total as i64,
        );
    }

    // ── SPRING ──────────────────────────────────────────────
    println!();
    println!("=== SPRING COMPARISON ===");
    println!();

    let spring_bin = "/home/parlar_ai/dev/qz/.pixi/envs/default/bin/spring";
    let tmp_dir = "/tmp/qz_bench_spring";
    let spring_out = format!("{}/spring.compressed", tmp_dir);
    let spring_decomp = format!("{}/spring_decomp.fastq", tmp_dir);

    // Clean up any previous run
    let _ = std::fs::remove_dir_all(tmp_dir);
    std::fs::create_dir_all(tmp_dir).unwrap();

    // Compress with SPRING
    eprintln!("Running SPRING compress...");
    let t = Instant::now();
    let spring_compress = std::process::Command::new("/usr/bin/time")
        .args(["-v", spring_bin, "-c",
            "-i", &fastq_path,
            "-o", &spring_out,
            "-w", tmp_dir,
            "-t", "8"])
        .output()
        .expect("Failed to run SPRING");
    let spring_compress_time = t.elapsed();

    let spring_stderr = String::from_utf8_lossy(&spring_compress.stderr);
    let spring_compress_mem = parse_max_rss(&spring_stderr);

    if !spring_compress.status.success() {
        eprintln!("SPRING compress failed: {}", spring_stderr);
    }

    let spring_size = std::fs::metadata(&spring_out).map(|m| m.len()).unwrap_or(0) as usize;

    // Decompress with SPRING
    eprintln!("Running SPRING decompress...");
    let t = Instant::now();
    let spring_decompress = std::process::Command::new("/usr/bin/time")
        .args(["-v", spring_bin, "-d",
            "-i", &spring_out,
            "-o", &spring_decomp,
            "-w", tmp_dir,
            "-t", "8"])
        .output()
        .expect("Failed to run SPRING decompress");
    let spring_decompress_time = t.elapsed();

    let spring_decomp_stderr = String::from_utf8_lossy(&spring_decompress.stderr);
    let spring_decompress_mem = parse_max_rss(&spring_decomp_stderr);

    // Compress with QZ default
    eprintln!("Running QZ default compress...");
    let qz_out = format!("{}/qz_default.qz", tmp_dir);
    let qz_decomp = format!("{}/qz_default.fastq", tmp_dir);
    let qz_bin = "/home/parlar_ai/dev/qz/target/release/qz";

    let t_qz_c = Instant::now();
    let qz_compress = std::process::Command::new("/usr/bin/time")
        .args(["-v", qz_bin, "compress",
            "-i", &fastq_path,
            "-o", &qz_out,
            "-w", tmp_dir])
        .output()
        .expect("Failed to run QZ compress");
    let qz_compress_time = t_qz_c.elapsed();
    let qz_stderr = String::from_utf8_lossy(&qz_compress.stderr);
    let qz_compress_mem = parse_max_rss(&qz_stderr);
    let qz_file_size = std::fs::metadata(&qz_out).map(|m| m.len()).unwrap_or(0) as usize;

    // Decompress with QZ
    eprintln!("Running QZ default decompress...");
    let t_qz_d = Instant::now();
    let qz_decompress = std::process::Command::new("/usr/bin/time")
        .args(["-v", qz_bin, "decompress",
            "-i", &qz_out,
            "-o", &qz_decomp,
            "-w", tmp_dir])
        .output()
        .expect("Failed to run QZ decompress");
    let qz_decompress_time = t_qz_d.elapsed();
    let qz_decomp_stderr = String::from_utf8_lossy(&qz_decompress.stderr);
    let qz_decompress_mem = parse_max_rss(&qz_decomp_stderr);

    // Print comparison table
    println!(
        "{:<30} {:>10} {:>7} {:>10} {:>10} {:>10} {:>10}",
        "Tool", "Size", "Ratio", "C.Time", "D.Time", "C.Mem", "D.Mem"
    );
    println!("{}", "-".repeat(90));

    if qz_file_size > 0 {
        println!(
            "{:<30} {:>8} B {:>6.2}x {:>8.1}s {:>8.1}s {:>8} {:>8}",
            "QZ default (actual file)",
            qz_file_size,
            original_size as f64 / qz_file_size as f64,
            qz_compress_time.as_secs_f64(),
            qz_decompress_time.as_secs_f64(),
            format_mem(qz_compress_mem),
            format_mem(qz_decompress_mem),
        );
    }

    println!(
        "{:<30} {:>8} B {:>6.2}x {:>8.1}s {:>8} {:>8} {:>8}",
        "QZ shared-perm (projected)",
        qz_shared_total,
        original_size as f64 / qz_shared_total as f64,
        qz_shared_time.as_secs_f64(),
        "n/a",
        "~same", "~same",
    );

    if qz_template_total > 0 {
        println!(
            "{:<30} {:>8} B {:>6.2}x {:>8.1}s {:>8} {:>8} {:>8}",
            "QZ slow (projected)",
            qz_template_total,
            original_size as f64 / qz_template_total as f64,
            qz_template_time.as_secs_f64(),
            "n/a",
            "~same", "~same",
        );
    }

    if spring_size > 0 {
        println!(
            "{:<30} {:>8} B {:>6.2}x {:>8.1}s {:>8.1}s {:>8} {:>8}",
            "SPRING",
            spring_size,
            original_size as f64 / spring_size as f64,
            spring_compress_time.as_secs_f64(),
            spring_decompress_time.as_secs_f64(),
            format_mem(spring_compress_mem),
            format_mem(spring_decompress_mem),
        );
    }

    println!();
    println!(
        "Input: {} reads x {} bp = {} bytes ({:.1} MB)",
        num_reads, read_len, original_size,
        original_size as f64 / (1024.0 * 1024.0),
    );

    // Cleanup
    let _ = std::fs::remove_dir_all(tmp_dir);
}

/// Parse "Maximum resident set size (kbytes): 12345" from /usr/bin/time -v output
fn parse_max_rss(stderr: &str) -> Option<u64> {
    for line in stderr.lines() {
        if line.contains("Maximum resident set size") {
            if let Some(val) = line.split(':').last() {
                return val.trim().parse().ok();
            }
        }
    }
    None
}

/// Format memory in MB from kbytes
fn format_mem(kb: Option<u64>) -> String {
    match kb {
        Some(k) => format!("{} MB", k / 1024),
        None => "n/a".to_string(),
    }
}
