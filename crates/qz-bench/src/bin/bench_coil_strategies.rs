/// Benchmark: test coil strategies on real FASTQ data
///
/// 1. Quality context_residuals + BSC vs raw qualities + BSC
/// 2. Lexicographic sequence sorting + BSC vs unsorted + BSC
/// 3. Consensus reference stats (to confirm it doesn't help on unsorted WGS)

use std::io::BufRead;
use std::time::Instant;

const MAX_QUALITY: usize = 46;
const NUM_POS_BINS: usize = 15;
const NUM_PREV_VALS: usize = MAX_QUALITY;
const NUM_CONTEXTS: usize = NUM_POS_BINS * NUM_PREV_VALS;
const RESCALE_THRESHOLD: u32 = 1 << 15;

fn main() {
    let fastq_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/NA12878_exome_50k.fastq".to_string());

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
    eprintln!("Loaded {} reads", num_reads);

    // === Strategy 1: Quality context residuals + BSC ===
    println!("=== QUALITY COMPRESSION ===\n");

    // Raw quality stream (current QZ approach: varint-framed, 7-bit packed)
    let raw_qual_stream: Vec<u8> = qualities
        .iter()
        .flat_map(|q| q.as_bytes().iter().copied())
        .collect();

    let t = Instant::now();
    let raw_bsc = qz_lib::compression::bsc::compress_parallel_adaptive(&raw_qual_stream).unwrap();
    let raw_time = t.elapsed();
    println!(
        "Raw qualities + BSC:          {:>10} B  ({:.2}x)  {:.1}ms",
        raw_bsc.len(),
        raw_qual_stream.len() as f64 / raw_bsc.len() as f64,
        raw_time.as_secs_f64() * 1000.0,
    );

    // Context residuals + BSC
    let t = Instant::now();
    let residuals = context_residuals(&qualities);
    let residual_time = t.elapsed();

    let t = Instant::now();
    let residual_bsc = qz_lib::compression::bsc::compress_parallel_adaptive(&residuals).unwrap();
    let bsc_time = t.elapsed();

    // Verify roundtrip
    let dec_residuals =
        qz_lib::compression::bsc::decompress_parallel(&residual_bsc).unwrap();
    let reconstructed = reconstruct_from_residuals(&dec_residuals, &qualities);
    assert_eq!(qualities, reconstructed, "Quality residual roundtrip FAILED");

    println!(
        "Context residuals + BSC:      {:>10} B  ({:.2}x)  {:.1}ms (predict: {:.1}ms + BSC: {:.1}ms)",
        residual_bsc.len(),
        raw_qual_stream.len() as f64 / residual_bsc.len() as f64,
        (residual_time + bsc_time).as_secs_f64() * 1000.0,
        residual_time.as_secs_f64() * 1000.0,
        bsc_time.as_secs_f64() * 1000.0,
    );

    let improvement = 100.0 * (1.0 - residual_bsc.len() as f64 / raw_bsc.len() as f64);
    println!("  -> {:.1}% smaller than raw+BSC\n", improvement);

    // === Strategy 2: Lexicographic sorting + BSC for sequences ===
    println!("=== SEQUENCE COMPRESSION ===\n");

    // Unsorted (current approach)
    let unsorted_stream: Vec<u8> = sequences
        .iter()
        .flat_map(|s| s.as_bytes().iter().copied())
        .collect();

    let t = Instant::now();
    let unsorted_bsc = qz_lib::compression::bsc::compress_parallel_adaptive(&unsorted_stream).unwrap();
    let unsorted_time = t.elapsed();
    println!(
        "Unsorted sequences + BSC:     {:>10} B  ({:.2}x)  {:.1}ms",
        unsorted_bsc.len(),
        unsorted_stream.len() as f64 / unsorted_bsc.len() as f64,
        unsorted_time.as_secs_f64() * 1000.0,
    );

    // Lexicographic sort
    let t = Instant::now();
    let mut sorted_seqs = sequences.clone();
    sorted_seqs.sort_unstable();
    let sort_time = t.elapsed();

    let sorted_stream: Vec<u8> = sorted_seqs
        .iter()
        .flat_map(|s| s.as_bytes().iter().copied())
        .collect();

    let t = Instant::now();
    let sorted_bsc = qz_lib::compression::bsc::compress_parallel_adaptive(&sorted_stream).unwrap();
    let sorted_bsc_time = t.elapsed();
    println!(
        "Lex-sorted sequences + BSC:   {:>10} B  ({:.2}x)  {:.1}ms (sort: {:.1}ms + BSC: {:.1}ms)",
        sorted_bsc.len(),
        sorted_stream.len() as f64 / sorted_bsc.len() as f64,
        (sort_time + sorted_bsc_time).as_secs_f64() * 1000.0,
        sort_time.as_secs_f64() * 1000.0,
        sorted_bsc_time.as_secs_f64() * 1000.0,
    );

    let seq_improvement = 100.0 * (1.0 - sorted_bsc.len() as f64 / unsorted_bsc.len() as f64);
    println!("  -> {:.1}% smaller than unsorted+BSC\n", seq_improvement);

    // === Strategy 3: Consensus reference analysis ===
    println!("=== CONSENSUS REFERENCE ANALYSIS ===\n");

    // Check if all same length
    let first_len = sequences[0].len();
    let all_same_len = sequences.iter().all(|s| s.len() == first_len);

    if all_same_len {
        let consensus = build_consensus(&sequences);
        let mut total_diffs = 0usize;

        for seq in &sequences {
            for (a, b) in seq.as_bytes().iter().zip(consensus.iter()) {
                if a != b {
                    total_diffs += 1;
                }
            }
        }

        let avg_diffs = total_diffs as f64 / num_reads as f64;
        let avg_pct = avg_diffs / first_len as f64 * 100.0;

        println!("Read length: {}", first_len);
        println!("Avg diffs from consensus: {:.1} per read ({:.1}%)", avg_diffs, avg_pct);
        println!(
            "Consensus encoding size estimate: {} bytes (vs {} raw)",
            total_diffs * 2 + num_reads,  // 2 bytes per diff + 1 varint per read
            unsorted_stream.len(),
        );

        if avg_pct > 50.0 {
            println!("  -> Consensus is WORSE than raw (>{:.0}% diffs = expansion)", 50.0);
        } else {
            println!("  -> Consensus would save {:.1}% before secondary compression",
                100.0 * (1.0 - (total_diffs * 2 + num_reads) as f64 / unsorted_stream.len() as f64));
        }
    } else {
        println!("Variable read lengths - consensus not applicable");
    }
}

fn build_consensus(sequences: &[String]) -> Vec<u8> {
    let read_length = sequences[0].len();
    let mut consensus = Vec::with_capacity(read_length);

    for pos in 0..read_length {
        let mut counts = [0u32; 5];
        for seq in sequences {
            let base = seq.as_bytes()[pos];
            let idx = match base {
                b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 4,
            };
            counts[idx] += 1;
        }
        let max_idx = counts.iter().enumerate().max_by_key(|(_, c)| *c).unwrap().0;
        consensus.push(match max_idx {
            0 => b'A', 1 => b'C', 2 => b'G', 3 => b'T', _ => b'N',
        });
    }
    consensus
}

/// Context-adaptive quality prediction residuals (from coil quality_rans.rs)
fn context_residuals(qualities: &[String]) -> Vec<u8> {
    let total_len: usize = qualities.iter().map(|q| q.len()).sum();
    let mut residuals = Vec::with_capacity(total_len);
    let mut counts = vec![[1u32; MAX_QUALITY]; NUM_CONTEXTS];

    for qual_str in qualities {
        let mut prev_quality: u8 = 0;
        for (pos, &byte) in qual_str.as_bytes().iter().enumerate() {
            let phred = byte.saturating_sub(33).min(MAX_QUALITY as u8 - 1);
            let symbol = phred as usize;

            let pos_bin = (pos / 10).min(NUM_POS_BINS - 1);
            let prev_val = (prev_quality as usize).min(NUM_PREV_VALS - 1);
            let ctx = pos_bin * NUM_PREV_VALS + prev_val;

            // Find most probable symbol as prediction
            let predicted = counts[ctx]
                .iter()
                .enumerate()
                .max_by_key(|&(_, &c)| c)
                .map(|(i, _)| i)
                .unwrap_or(0);

            // Residual: actual - predicted, shifted to unsigned
            let residual = (symbol as i32 - predicted as i32 + 128) as u8;
            residuals.push(residual);

            // Update model
            counts[ctx][symbol] += 1;
            let total: u32 = counts[ctx].iter().sum();
            if total >= RESCALE_THRESHOLD {
                for c in counts[ctx].iter_mut() {
                    *c = (*c + 1) / 2;
                }
            }

            prev_quality = phred;
        }
    }

    residuals
}

/// Reconstruct original qualities from residuals (for verification)
fn reconstruct_from_residuals(residuals: &[u8], original_qualities: &[String]) -> Vec<String> {
    let mut result = Vec::with_capacity(original_qualities.len());
    let mut counts = vec![[1u32; MAX_QUALITY]; NUM_CONTEXTS];
    let mut r_idx = 0;

    for orig_qual in original_qualities {
        let mut prev_quality: u8 = 0;
        let mut qual = String::with_capacity(orig_qual.len());

        for pos in 0..orig_qual.len() {
            let pos_bin = (pos / 10).min(NUM_POS_BINS - 1);
            let prev_val = (prev_quality as usize).min(NUM_PREV_VALS - 1);
            let ctx = pos_bin * NUM_PREV_VALS + prev_val;

            let predicted = counts[ctx]
                .iter()
                .enumerate()
                .max_by_key(|&(_, &c)| c)
                .map(|(i, _)| i)
                .unwrap_or(0);

            let residual = residuals[r_idx] as i32;
            let symbol = (residual - 128 + predicted as i32) as u8;
            r_idx += 1;

            qual.push((symbol + 33) as char);

            counts[ctx][symbol as usize] += 1;
            let total: u32 = counts[ctx].iter().sum();
            if total >= RESCALE_THRESHOLD {
                for c in counts[ctx].iter_mut() {
                    *c = (*c + 1) / 2;
                }
            }

            prev_quality = symbol;
        }

        result.push(qual);
    }

    result
}
