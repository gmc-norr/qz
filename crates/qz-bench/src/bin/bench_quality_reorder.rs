/// Benchmark: Quality score compression with sequence-derived reordering
///
/// Tests whether sorting quality strings by sequence-derived sort keys
/// (syncmer hashes) improves BSC compression.
/// The permutation is implicit (derived from sequences), so zero storage overhead.

use std::io::BufRead;
use std::time::Instant;

use qz_lib::compression::bsc;
use qz_lib::compression::dna_utils::{kmer_to_hash, reverse_complement_hash};
use qz_lib::compression::fqzcomp;
use qz_lib::compression::quality_context::{self, QualityContextConfig};
use rayon::prelude::*;

fn main() {
    let fastq_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/ERR3239334_1.500k.fastq".to_string());

    eprintln!("Reading FASTQ: {}", fastq_path);
    let file = std::fs::File::open(&fastq_path).expect("Cannot open FASTQ file");
    let reader = std::io::BufReader::new(file);

    let mut sequences: Vec<String> = Vec::new();
    let mut qualities: Vec<String> = Vec::new();
    let mut line_num = 0u64;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        match line_num % 4 {
            1 => sequences.push(line.trim_end().to_string()),
            3 => qualities.push(line.trim_end().to_string()),
            _ => {}
        }
        line_num += 1;
    }

    let num_reads = sequences.len();
    let total_qual_bytes: usize = qualities.iter().map(|q| q.len()).sum();
    let read_len = sequences.first().map(|s| s.len()).unwrap_or(0);

    eprintln!(
        "Loaded {} reads, read_len={}, quality bytes={} ({:.1} MB)\n",
        num_reads, read_len, total_qual_bytes,
        total_qual_bytes as f64 / (1024.0 * 1024.0),
    );

    println!(
        "{:<55} {:>12} {:>8} {:>10}",
        "Method", "Compressed", "Ratio", "Time"
    );
    println!("{}", "-".repeat(88));

    // === 1. Baseline: concatenated qualities + BSC (current QZ approach) ===
    let baseline_size;
    {
        let raw: Vec<u8> = qualities.iter().flat_map(|q| q.as_bytes()).copied().collect();
        let t = Instant::now();
        let compressed = bsc::compress_parallel_adaptive(&raw).unwrap();
        let elapsed = t.elapsed();
        baseline_size = compressed.len();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms",
            "1. Baseline (concat + BSC adaptive)",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
        );
    }

    // === 2. Compute syncmer sort keys from sequences ===
    let anchor_k = 31usize;
    let s = anchor_k.saturating_sub(3).max(5);
    let t_end = anchor_k - s;

    eprintln!("\nComputing syncmer sort keys (k={}, s={})...", anchor_k, s);
    let t = Instant::now();
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
    eprintln!("Sort keys computed: {:.2}s", t.elapsed().as_secs_f64());

    // Build permutation (stable sort)
    let t = Instant::now();
    let mut perm: Vec<u32> = (0..num_reads as u32).collect();
    perm.sort_by_key(|&i| sort_keys[i as usize]);
    eprintln!("Argsort: {:.2}s", t.elapsed().as_secs_f64());

    // Count unique sort keys to see grouping quality
    let mut unique_keys = sort_keys.clone();
    unique_keys.sort();
    unique_keys.dedup();
    eprintln!("Unique sort keys: {} ({:.1}% of reads)\n",
             unique_keys.len(),
             unique_keys.len() as f64 / num_reads as f64 * 100.0);

    // === 3. Reordered qualities + BSC ===
    {
        let mut reordered = Vec::with_capacity(total_qual_bytes);
        for &idx in &perm {
            reordered.extend_from_slice(qualities[idx as usize].as_bytes());
        }
        let t = Instant::now();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "2. Syncmer-reordered (k=31) + BSC adaptive",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );
    }

    // === 3. Lexicographic sort of quality strings + permutation cost ===
    {
        let mut perm_lex: Vec<u32> = (0..num_reads as u32).collect();
        perm_lex.sort_by(|&a, &b| qualities[a as usize].cmp(&qualities[b as usize]));

        let mut reordered = Vec::with_capacity(total_qual_bytes);
        for &idx in &perm_lex {
            reordered.extend_from_slice(qualities[idx as usize].as_bytes());
        }
        let t = Instant::now();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();
        let elapsed = t.elapsed();

        // Compute permutation storage cost
        let perm_raw: Vec<u8> = perm_lex.iter().flat_map(|&i| i.to_le_bytes()).collect();
        let perm_bsc = bsc::compress_parallel_adaptive(&perm_raw).unwrap();

        let total_with_perm = compressed.len() + perm_bsc.len();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let net_savings = baseline_size as i64 - total_with_perm as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  (qual only)",
            "3a. Lex-sorted qualities + BSC adaptive",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
        );
        println!(
            "{:<55} {:>10} B",
            "    Permutation (BSC'd)",
            perm_bsc.len(),
        );
        println!(
            "{:<55} {:>10} B {:>7.2}x            ({:+} B net)",
            "    TOTAL (qual + perm)",
            total_with_perm,
            total_qual_bytes as f64 / total_with_perm as f64,
            -net_savings,
        );
    }

    // === 3b. Lex sort with delta-encoded permutation ===
    {
        let mut perm_lex: Vec<u32> = (0..num_reads as u32).collect();
        perm_lex.sort_by(|&a, &b| qualities[a as usize].cmp(&qualities[b as usize]));

        // Delta-encode the permutation (consecutive indices may be close)
        let mut perm_delta: Vec<u8> = Vec::with_capacity(num_reads * 3);
        let mut prev = 0i64;
        for &idx in &perm_lex {
            let delta = idx as i64 - prev;
            prev = idx as i64;
            // Zigzag encode
            let zz = ((delta << 1) ^ (delta >> 63)) as u64;
            // Varint encode
            let mut v = zz;
            while v >= 0x80 {
                perm_delta.push((v as u8) | 0x80);
                v >>= 7;
            }
            perm_delta.push(v as u8);
        }
        let perm_delta_bsc = bsc::compress_parallel_adaptive(&perm_delta).unwrap();

        let mut reordered = Vec::with_capacity(total_qual_bytes);
        for &idx in &perm_lex {
            reordered.extend_from_slice(qualities[idx as usize].as_bytes());
        }
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();
        let total = compressed.len() + perm_delta_bsc.len();
        let net_savings = baseline_size as i64 - total as i64;
        println!(
            "{:<55} {:>10} B",
            "3b. Perm delta-zigzag-varint (BSC'd)",
            perm_delta_bsc.len(),
        );
        println!(
            "{:<55} {:>10} B {:>7.2}x            ({:+} B net)",
            "    TOTAL (qual + delta perm)",
            total,
            total_qual_bytes as f64 / total as f64,
            -net_savings,
        );
    }

    // === 5. Sort by mean quality value ===
    {
        let mut perm_mean: Vec<u32> = (0..num_reads as u32).collect();
        perm_mean.sort_by_key(|&i| {
            let q = qualities[i as usize].as_bytes();
            let sum: u64 = q.iter().map(|&b| b as u64).sum();
            sum / q.len().max(1) as u64
        });

        let mut reordered = Vec::with_capacity(total_qual_bytes);
        for &idx in &perm_mean {
            reordered.extend_from_slice(qualities[idx as usize].as_bytes());
        }
        let t = Instant::now();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  (qual only, {:+} B)",
            "4a. Sort by mean quality + BSC adaptive",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );

        // Permutation cost (raw BSC'd)
        let perm_raw: Vec<u8> = perm_mean.iter().flat_map(|&i| i.to_le_bytes()).collect();
        let perm_bsc = bsc::compress_parallel_adaptive(&perm_raw).unwrap();
        let total_with_perm = compressed.len() + perm_bsc.len();
        let net_savings = baseline_size as i64 - total_with_perm as i64;
        println!(
            "{:<55} {:>10} B",
            "    Permutation (BSC'd)",
            perm_bsc.len(),
        );
        println!(
            "{:<55} {:>10} B {:>7.2}x            ({:+} B net)",
            "    TOTAL (qual + perm)",
            total_with_perm,
            total_qual_bytes as f64 / total_with_perm as f64,
            -net_savings,
        );

        // Delta-zigzag-varint permutation
        let mut perm_delta: Vec<u8> = Vec::with_capacity(num_reads * 3);
        let mut prev = 0i64;
        for &idx in &perm_mean {
            let delta = idx as i64 - prev;
            prev = idx as i64;
            let zz = ((delta << 1) ^ (delta >> 63)) as u64;
            let mut v = zz;
            while v >= 0x80 {
                perm_delta.push((v as u8) | 0x80);
                v >>= 7;
            }
            perm_delta.push(v as u8);
        }
        let perm_delta_bsc = bsc::compress_parallel_adaptive(&perm_delta).unwrap();
        let total_delta = compressed.len() + perm_delta_bsc.len();
        let net_delta = baseline_size as i64 - total_delta as i64;
        println!(
            "{:<55} {:>10} B",
            "    Perm delta-zigzag-varint (BSC'd)",
            perm_delta_bsc.len(),
        );
        println!(
            "{:<55} {:>10} B {:>7.2}x            ({:+} B net)",
            "    TOTAL (qual + delta perm)",
            total_delta,
            total_qual_bytes as f64 / total_delta as f64,
            -net_delta,
        );

        // Mean quality as sort key: can we reconstruct from data?
        // Store just the mean quality per read (1 byte each) + BSC
        let mean_quals: Vec<u8> = qualities.iter().map(|q| {
            let bytes = q.as_bytes();
            let sum: u64 = bytes.iter().map(|&b| b as u64).sum();
            (sum / bytes.len().max(1) as u64) as u8
        }).collect();
        let mean_bsc = bsc::compress_parallel_adaptive(&mean_quals).unwrap();
        let total_mean_key = compressed.len() + mean_bsc.len();
        let net_mean = baseline_size as i64 - total_mean_key as i64;
        println!(
            "{:<55} {:>10} B",
            "    Mean-qual sort key (BSC'd, 1B/read)",
            mean_bsc.len(),
        );
        println!(
            "{:<55} {:>10} B {:>7.2}x            ({:+} B net)",
            "    TOTAL (qual + mean key)",
            total_mean_key,
            total_qual_bytes as f64 / total_mean_key as f64,
            -net_mean,
        );
    }

    // === 6. Sort by first 10 quality bytes (fast proxy) ===
    {
        let mut perm_prefix: Vec<u32> = (0..num_reads as u32).collect();
        perm_prefix.sort_by(|&a, &b| {
            let qa = qualities[a as usize].as_bytes();
            let qb = qualities[b as usize].as_bytes();
            let la = qa.len().min(10);
            let lb = qb.len().min(10);
            qa[..la].cmp(&qb[..lb])
        });

        let mut reordered = Vec::with_capacity(total_qual_bytes);
        for &idx in &perm_prefix {
            reordered.extend_from_slice(qualities[idx as usize].as_bytes());
        }
        let t = Instant::now();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "5. Sort by first-10 quality bytes + BSC adaptive",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );
    }

    // === 7. Syncmer sort key + BSC static QLFC (like SPRING uses) ===
    {
        let mut reordered = Vec::with_capacity(total_qual_bytes);
        for &idx in &perm {
            reordered.extend_from_slice(qualities[idx as usize].as_bytes());
        }
        let t = Instant::now();
        let compressed = bsc::compress_parallel(&reordered).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "6. Syncmer-reordered + BSC static",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );
    }

    // === 7. fqzcomp_qual strat=0 (fast) ===
    {
        let qual_refs: Vec<&[u8]> = qualities.iter().map(|q| q.as_bytes()).collect();
        let t = Instant::now();
        let compressed = fqzcomp::compress(&qual_refs, 0).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "7. fqzcomp_qual strat=0 (fast)",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );

        // Verify roundtrip
        let (decompressed, _lengths) = fqzcomp::decompress(&compressed, num_reads).unwrap();
        let concat: Vec<u8> = qualities.iter().flat_map(|q| q.as_bytes()).copied().collect();
        assert_eq!(decompressed, concat, "fqzcomp strat=0 roundtrip failed!");
    }

    // === 8. fqzcomp_qual strat=1 ===
    {
        let qual_refs: Vec<&[u8]> = qualities.iter().map(|q| q.as_bytes()).collect();
        let t = Instant::now();
        let compressed = fqzcomp::compress(&qual_refs, 1).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "8. fqzcomp_qual strat=1",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );
    }

    // === 9. fqzcomp_qual strat=2 ===
    {
        let qual_refs: Vec<&[u8]> = qualities.iter().map(|q| q.as_bytes()).collect();
        let t = Instant::now();
        let compressed = fqzcomp::compress(&qual_refs, 2).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "9. fqzcomp_qual strat=2",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );
    }

    // === 10. fqzcomp_qual strat=3 (max) ===
    {
        let qual_refs: Vec<&[u8]> = qualities.iter().map(|q| q.as_bytes()).collect();
        let t = Instant::now();
        let compressed = fqzcomp::compress(&qual_refs, 3).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "10. fqzcomp_qual strat=3 (max)",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );
    }

    // === 10b. Mean-quality reorder + fqzcomp strat=0 ===
    {
        let mut perm_mean: Vec<u32> = (0..num_reads as u32).collect();
        perm_mean.sort_by_key(|&i| {
            let q = qualities[i as usize].as_bytes();
            let sum: u64 = q.iter().map(|&b| b as u64).sum();
            sum / q.len().max(1) as u64
        });

        let reordered_quals: Vec<&[u8]> = perm_mean.iter().map(|&i| qualities[i as usize].as_bytes()).collect();
        let t = Instant::now();
        let compressed = fqzcomp::compress(&reordered_quals, 0).unwrap();
        let elapsed = t.elapsed();

        // Mean-qual sort key cost
        let mean_quals: Vec<u8> = qualities.iter().map(|q| {
            let bytes = q.as_bytes();
            let sum: u64 = bytes.iter().map(|&b| b as u64).sum();
            (sum / bytes.len().max(1) as u64) as u8
        }).collect();
        let mean_bsc = bsc::compress_parallel_adaptive(&mean_quals).unwrap();

        let total = compressed.len() + mean_bsc.len();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let net = baseline_size as i64 - total as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  (qual only, {:+} B)",
            "10b. Mean-sort + fqzcomp strat=0",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -(baseline_size as i64 - compressed.len() as i64),
        );
        println!(
            "{:<55} {:>10} B {:>7.2}x            ({:+} B net)",
            "     TOTAL (qual + mean key)",
            total,
            total_qual_bytes as f64 / total as f64,
            -net,
        );
    }

    // === 10c. Lex-sort + fqzcomp strat=0 ===
    {
        let mut perm_lex: Vec<u32> = (0..num_reads as u32).collect();
        perm_lex.sort_by(|&a, &b| qualities[a as usize].cmp(&qualities[b as usize]));

        let reordered_quals: Vec<&[u8]> = perm_lex.iter().map(|&i| qualities[i as usize].as_bytes()).collect();
        let t = Instant::now();
        let compressed = fqzcomp::compress(&reordered_quals, 0).unwrap();
        let elapsed = t.elapsed();

        // Mean-qual sort key cost (cheap proxy for lex)
        let mean_quals: Vec<u8> = qualities.iter().map(|q| {
            let bytes = q.as_bytes();
            let sum: u64 = bytes.iter().map(|&b| b as u64).sum();
            (sum / bytes.len().max(1) as u64) as u8
        }).collect();
        let _mean_bsc = bsc::compress_parallel_adaptive(&mean_quals).unwrap();

        // Delta perm cost
        let mut perm_delta: Vec<u8> = Vec::with_capacity(num_reads * 3);
        let mut prev = 0i64;
        for &idx in &perm_lex {
            let delta = idx as i64 - prev;
            prev = idx as i64;
            let zz = ((delta << 1) ^ (delta >> 63)) as u64;
            let mut v = zz;
            while v >= 0x80 {
                perm_delta.push((v as u8) | 0x80);
                v >>= 7;
            }
            perm_delta.push(v as u8);
        }
        let perm_delta_bsc = bsc::compress_parallel_adaptive(&perm_delta).unwrap();

        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  (qual only, {:+} B)",
            "10c. Lex-sort + fqzcomp strat=0",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -(baseline_size as i64 - compressed.len() as i64),
        );
        let total_delta = compressed.len() + perm_delta_bsc.len();
        println!(
            "{:<55} {:>10} B {:>7.2}x            ({:+} B net)",
            "     TOTAL (qual + delta perm)",
            total_delta,
            total_qual_bytes as f64 / total_delta as f64,
            -(baseline_size as i64 - total_delta as i64),
        );
    }

    // === 11. Context-modeled adaptive arithmetic (pure Rust, single-threaded) ===
    {
        let concat: Vec<u8> = qualities.iter().flat_map(|q| q.as_bytes()).copied().collect();
        let rlens: Vec<u32> = qualities.iter().map(|q| q.len() as u32).collect();
        let config = QualityContextConfig::default(); // pos_bits=6, qual_bits=4

        let t = Instant::now();
        let compressed = quality_context::compress(&concat, &rlens, &config).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "11. Context-arith (pos6+qual4, single-thread)",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );

        // Verify roundtrip
        let decompressed = quality_context::decompress(&compressed, &rlens, &config).unwrap();
        assert_eq!(decompressed, concat, "context-arith roundtrip failed!");
    }

    // === 12. Context-modeled adaptive arithmetic (parallel blocks) ===
    {
        let concat: Vec<u8> = qualities.iter().flat_map(|q| q.as_bytes()).copied().collect();
        let rlens: Vec<u32> = qualities.iter().map(|q| q.len() as u32).collect();
        let config = QualityContextConfig::default();

        let t = Instant::now();
        let compressed = quality_context::compress_parallel(&concat, &rlens, &config).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "12. Context-arith (pos6+qual4, parallel 100K blocks)",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );

        // Verify roundtrip
        let decompressed = quality_context::decompress_parallel(&compressed, &rlens, &config).unwrap();
        assert_eq!(decompressed, concat, "context-arith parallel roundtrip failed!");
    }

    // === 13. Context-modeled with more position bits ===
    {
        let concat: Vec<u8> = qualities.iter().flat_map(|q| q.as_bytes()).copied().collect();
        let rlens: Vec<u32> = qualities.iter().map(|q| q.len() as u32).collect();
        let config = QualityContextConfig { pos_bits: 7, qual_bits: 4 };

        let t = Instant::now();
        let compressed = quality_context::compress(&concat, &rlens, &config).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "13. Context-arith (pos7+qual4, single-thread)",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );
    }

    // === 14. Context-modeled with more quality bits ===
    {
        let concat: Vec<u8> = qualities.iter().flat_map(|q| q.as_bytes()).copied().collect();
        let rlens: Vec<u32> = qualities.iter().map(|q| q.len() as u32).collect();
        let config = QualityContextConfig { pos_bits: 6, qual_bits: 5 };

        let t = Instant::now();
        let compressed = quality_context::compress(&concat, &rlens, &config).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_qual_bytes as f64 / compressed.len() as f64;
        let savings = baseline_size as i64 - compressed.len() as i64;
        println!(
            "{:<55} {:>10} B {:>7.2}x {:>8.1}ms  ({:+} B)",
            "14. Context-arith (pos6+qual5, single-thread)",
            compressed.len(), ratio, elapsed.as_secs_f64() * 1000.0,
            -savings,
        );
    }

    println!("{}", "-".repeat(88));
    println!(
        "Input: {} reads x {} bp = {} quality bytes ({:.1} MB)",
        num_reads, read_len, total_qual_bytes,
        total_qual_bytes as f64 / (1024.0 * 1024.0),
    );
}
