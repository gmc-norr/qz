/// Benchmark: Read sorting strategies for BSC compression
///
/// Tests whether sorting reads (sequences and/or qualities) by various
/// vector-space properties improves BSC compression. The goal is purely
/// to maximize byte-level redundancy for BWT — no biological meaning needed.
///
/// Strategies tested:
///   1. Baseline (original order)
///   2. Lexicographic sort
///   3. Sort by mean ASCII value (≈ GC content for sequences)
///   4. Sort by L2 norm (sum of squares — weights extremes)
///   5. Multi-key: (mean, variance)
///   6. Multi-key: (mean, min value)
///   7. Quantized profile: downsample to N bins, lex-sort that vector
///   8. Sort by first + last 16 bytes
///   9. Reverse-complement canonical + lex (sequences only)
///  10. Minimizer sort (sequences only, for reference)
///
/// For each strategy, reports:
///   - Compressed size of the sorted stream
///   - Sort key / permutation storage cost (BSC'd)
///   - Net total = sorted_compressed + perm_cost

use std::io::BufRead;
use std::time::Instant;

use qz_lib::compression::bsc;
use qz_lib::compression::dna_utils::{kmer_to_hash, reverse_complement_hash};
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
    let read_len = sequences.first().map(|s| s.len()).unwrap_or(0);
    let total_seq_bytes: usize = sequences.iter().map(|s| s.len()).sum();
    let total_qual_bytes: usize = qualities.iter().map(|q| q.len()).sum();

    eprintln!(
        "Loaded {} reads, read_len={}, seq={:.1} MB, qual={:.1} MB\n",
        num_reads,
        read_len,
        total_seq_bytes as f64 / (1024.0 * 1024.0),
        total_qual_bytes as f64 / (1024.0 * 1024.0),
    );

    // ========================================================================
    // SEQUENCE SORTING STRATEGIES
    // ========================================================================
    println!("=== SEQUENCE SORTING STRATEGIES ===\n");
    println!(
        "{:<50} {:>12} {:>8} {:>12} {:>12} {:>10}",
        "Strategy", "Seq BSC", "Ratio", "Perm cost", "Net total", "Time"
    );
    println!("{}", "-".repeat(108));

    // Pre-compute quality baseline (needed by three-level key section)
    let qual_baseline = {
        let raw: Vec<u8> = qualities.iter().flat_map(|q| q.as_bytes()).copied().collect();
        bsc::compress_parallel_adaptive(&raw).unwrap().len()
    };

    // Baseline: no sort
    let seq_baseline;
    {
        let raw: Vec<u8> = sequences.iter().flat_map(|s| s.as_bytes()).copied().collect();
        let t = Instant::now();
        let compressed = bsc::compress_parallel_adaptive(&raw).unwrap();
        let elapsed = t.elapsed();
        seq_baseline = compressed.len();
        print_result(
            "1. Baseline (original order)",
            total_seq_bytes,
            compressed.len(),
            0,
            elapsed,
        );
    }

    // --- Sort helpers ---
    // Each strategy produces a permutation. We measure:
    //   a) compressed size of reordered stream
    //   b) permutation cost (delta-zigzag-varint + BSC)

    // 2. Lexicographic sort
    bench_sort_strategy(
        "2. Lexicographic",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        0, // qual_baseline filled later
        |_i, seq, _qual| seq.as_bytes().to_vec(),
        true,
        false,
    );

    // 3. Sort by mean ASCII value
    bench_sort_strategy(
        "3. Mean ASCII value",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        0,
        |_i, seq, _qual| {
            let bytes = seq.as_bytes();
            let mean = bytes.iter().map(|&b| b as u64).sum::<u64>() / bytes.len().max(1) as u64;
            vec![mean as u8]
        },
        true,
        false,
    );

    // 4. Sort by L2 norm (sum of squares)
    bench_sort_strategy(
        "4. L2 norm (sum of squares)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        0,
        |_i, seq, _qual| {
            let bytes = seq.as_bytes();
            let sos: u64 = bytes.iter().map(|&b| (b as u64) * (b as u64)).sum();
            sos.to_le_bytes().to_vec()
        },
        true,
        false,
    );

    // 5. Multi-key: (mean, variance)
    bench_sort_strategy(
        "5. Multi-key (mean, variance)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        0,
        |_i, seq, _qual| {
            let bytes = seq.as_bytes();
            let n = bytes.len().max(1) as u64;
            let sum: u64 = bytes.iter().map(|&b| b as u64).sum();
            let mean = sum / n;
            let var: u64 = bytes.iter().map(|&b| {
                let d = (b as i64) - (mean as i64);
                (d * d) as u64
            }).sum::<u64>() / n;
            let mut key = Vec::with_capacity(3);
            key.push(mean as u8);
            key.extend_from_slice(&(var as u16).to_le_bytes());
            key
        },
        true,
        false,
    );

    // 6. Multi-key: (mean, min)
    bench_sort_strategy(
        "6. Multi-key (mean, min)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        0,
        |_i, seq, _qual| {
            let bytes = seq.as_bytes();
            let mean = bytes.iter().map(|&b| b as u64).sum::<u64>() / bytes.len().max(1) as u64;
            let min = bytes.iter().copied().min().unwrap_or(0);
            vec![mean as u8, min]
        },
        true,
        false,
    );

    // 7. Quantized profile (8 bins)
    bench_sort_strategy(
        "7. Quantized profile (8 bins)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        0,
        |_i, seq, _qual| quantize_profile(seq.as_bytes(), 8),
        true,
        false,
    );

    // 7b. Quantized profile (16 bins)
    bench_sort_strategy(
        "7b. Quantized profile (16 bins)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        0,
        |_i, seq, _qual| quantize_profile(seq.as_bytes(), 16),
        true,
        false,
    );

    // 8. Sort by first + last 16 bytes
    bench_sort_strategy(
        "8. First+last 16 bytes",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        0,
        |_i, seq, _qual| {
            let b = seq.as_bytes();
            let n = b.len();
            let mut key = Vec::with_capacity(32);
            key.extend_from_slice(&b[..16.min(n)]);
            if n > 16 {
                key.extend_from_slice(&b[n - 16..]);
            }
            key
        },
        true,
        false,
    );

    // 9. Reverse-complement canonical + lex sort
    bench_sort_strategy(
        "9. RC-canonical + lex",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        0,
        |_i, seq, _qual| {
            let b = seq.as_bytes();
            let rc = reverse_complement_bytes(b);
            if rc < b.to_vec() { rc } else { b.to_vec() }
        },
        true,
        false,
    );

    // 10. Minimizer sort (k=15, w=10)
    {
        let t = Instant::now();
        let mut perm: Vec<u32> = (0..num_reads as u32).collect();
        let min_hashes: Vec<u64> = sequences
            .par_iter()
            .map(|s| canonical_minimizer(s.as_bytes(), 15, 10))
            .collect();
        perm.sort_unstable_by_key(|&i| min_hashes[i as usize]);
        let elapsed = t.elapsed();

        let reordered: Vec<u8> = perm
            .iter()
            .flat_map(|&i| sequences[i as usize].as_bytes())
            .copied()
            .collect();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();

        // No permutation cost — decompressor can recompute from sequences
        print_result(
            "10. Minimizer (k=15,w=10) [free perm]",
            total_seq_bytes,
            compressed.len(),
            0,
            elapsed,
        );
        print_delta(compressed.len(), 0, seq_baseline);
    }

    // 11. Syncmer sort (k=31, s=28)
    {
        let anchor_k = 31usize;
        let s = anchor_k.saturating_sub(3).max(5);
        let t_end = anchor_k - s;

        let t = Instant::now();
        let sort_keys: Vec<u64> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                if sb.len() < anchor_k { return u64::MAX; }
                let positions = syncmers::find_syncmers_pos(anchor_k, s, &[0, t_end], sb);
                let mut min_hash = u64::MAX;
                for pos in positions {
                    if pos + anchor_k > sb.len() { continue; }
                    let kmer = &sb[pos..pos + anchor_k];
                    if let Some(fwd) = kmer_to_hash(kmer) {
                        let rc = reverse_complement_hash(fwd, anchor_k);
                        let canon = fwd.min(rc);
                        if canon < min_hash { min_hash = canon; }
                    }
                }
                min_hash
            })
            .collect();
        let mut perm: Vec<u32> = (0..num_reads as u32).collect();
        perm.sort_unstable_by_key(|&i| sort_keys[i as usize]);
        let elapsed = t.elapsed();

        let reordered: Vec<u8> = perm
            .iter()
            .flat_map(|&i| sequences[i as usize].as_bytes())
            .copied()
            .collect();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();

        print_result(
            "11. Syncmer (k=31,s=28) [free perm]",
            total_seq_bytes,
            compressed.len(),
            0,
            elapsed,
        );
        print_delta(compressed.len(), 0, seq_baseline);
    }

    // ========================================================================
    // IMPROVED SYNCMER / MINIMIZER STRATEGIES
    // ========================================================================
    println!("\n\n=== IMPROVED SYNCMER / MINIMIZER STRATEGIES (all free perm) ===\n");
    println!(
        "{:<55} {:>12} {:>8} {:>10}",
        "Strategy", "Seq BSC", "Ratio", "Time"
    );
    println!("{}", "-".repeat(88));
    println!(
        "{:<55} {:>10} B {:>6.2}x",
        "Baseline (no sort)", seq_baseline, total_seq_bytes as f64 / seq_baseline as f64,
    );

    // Helper: run a sequence sort by u128 key and print result
    let bench_seq_sort_u128 = |name: &str, keys: &[u128]| {
        let t = Instant::now();
        let mut perm: Vec<u32> = (0..num_reads as u32).collect();
        perm.sort_unstable_by_key(|&i| keys[i as usize]);
        let elapsed = t.elapsed();
        let reordered: Vec<u8> = perm.iter().flat_map(|&i| sequences[i as usize].as_bytes()).copied().collect();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();
        let delta = compressed.len() as i64 - seq_baseline as i64;
        println!(
            "{:<55} {:>10} B {:>6.2}x {:>8.0}ms  ({:+} B)",
            name,
            compressed.len(),
            total_seq_bytes as f64 / compressed.len() as f64,
            elapsed.as_secs_f64() * 1000.0,
            delta,
        );
        compressed.len()
    };

    // Helper: sort by Vec<u8> key
    let bench_seq_sort_bytes = |name: &str, keys: &[Vec<u8>]| {
        let t = Instant::now();
        let mut perm: Vec<u32> = (0..num_reads as u32).collect();
        perm.sort_by(|&a, &b| keys[a as usize].cmp(&keys[b as usize]));
        let elapsed = t.elapsed();
        let reordered: Vec<u8> = perm.iter().flat_map(|&i| sequences[i as usize].as_bytes()).copied().collect();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();
        let delta = compressed.len() as i64 - seq_baseline as i64;
        println!(
            "{:<55} {:>10} B {:>6.2}x {:>8.0}ms  ({:+} B)",
            name,
            compressed.len(),
            total_seq_bytes as f64 / compressed.len() as f64,
            elapsed.as_secs_f64() * 1000.0,
            delta,
        );
        compressed.len()
    };

    // --- Reference: current best syncmer (k=31, s=28) ---
    let syncmer_k31_keys: Vec<u64> = sequences
        .par_iter()
        .map(|seq| compute_min_syncmer_hash(seq.as_bytes(), 31, 28))
        .collect();
    {
        let keys128: Vec<u128> = syncmer_k31_keys.iter().map(|&h| h as u128).collect();
        bench_seq_sort_u128("Ref: Syncmer k=31 s=28 (single hash)", &keys128);
    }

    // === 1. Multi-level: (min_hash, 2nd_min_hash) ===
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                let (h1, h2) = compute_top2_syncmer_hashes(sb, 31, 28);
                ((h1 as u128) << 64) | (h2 as u128)
            })
            .collect();
        bench_seq_sort_u128("1. Syncmer (min_hash, 2nd_hash) k=31", &keys);
    }

    // === 2. (min_hash, position_of_min) — approximates genomic offset ===
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                let (hash, pos) = compute_min_syncmer_hash_with_pos(sb, 31, 28);
                ((hash as u128) << 16) | (pos as u128)
            })
            .collect();
        bench_seq_sort_u128("2. Syncmer (min_hash, position) k=31", &keys);
    }

    // === 3. (min_hash, lex suffix after min position) ===
    {
        let keys: Vec<Vec<u8>> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                let (hash, pos) = compute_min_syncmer_hash_with_pos(sb, 31, 28);
                let mut key = hash.to_be_bytes().to_vec();
                // Append the read starting from the minimizer position
                let start = (pos as usize).min(sb.len());
                key.extend_from_slice(&sb[start..]);
                key
            })
            .collect();
        bench_seq_sort_bytes("3. Syncmer (min_hash, lex_suffix) k=31", &keys);
    }

    // === 4. MinHash sketch: sorted top-4 canonical hashes ===
    {
        let keys: Vec<Vec<u8>> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                let hashes = compute_sorted_syncmer_hashes(sb, 31, 28, 4);
                hashes.iter().flat_map(|h| h.to_be_bytes()).collect()
            })
            .collect();
        bench_seq_sort_bytes("4. MinHash sketch (top-4 syncmer hashes) k=31", &keys);
    }

    // === 5. MinHash sketch top-8 ===
    {
        let keys: Vec<Vec<u8>> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                let hashes = compute_sorted_syncmer_hashes(sb, 31, 28, 8);
                hashes.iter().flat_map(|h| h.to_be_bytes()).collect()
            })
            .collect();
        bench_seq_sort_bytes("5. MinHash sketch (top-8 syncmer hashes) k=31", &keys);
    }

    // === 6. Parameter sweep: different k values with syncmers ===
    println!();
    for k in [15usize, 21, 25, 31, 41] {
        let s = k.saturating_sub(3).max(5);
        let keys: Vec<u128> = sequences
            .par_iter()
            .map(|seq| compute_min_syncmer_hash(seq.as_bytes(), k, s) as u128)
            .collect();
        bench_seq_sort_u128(
            &format!("6. Syncmer k={} s={} (single hash)", k, s),
            &keys,
        );
    }

    // === 7. Parameter sweep: minimizer with different k, w ===
    println!();
    for (k, w) in [(11, 5), (15, 10), (21, 11), (25, 15), (31, 16)] {
        let keys: Vec<u128> = sequences
            .par_iter()
            .map(|seq| canonical_minimizer(seq.as_bytes(), k, w) as u128)
            .collect();
        bench_seq_sort_u128(
            &format!("7. Minimizer k={} w={}", k, w),
            &keys,
        );
    }

    // === 8. (minimizer, 2nd minimizer) multi-level ===
    println!();
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                let (h1, h2) = compute_top2_minimizer_hashes(sb, 15, 10);
                ((h1 as u128) << 64) | (h2 as u128)
            })
            .collect();
        bench_seq_sort_u128("8. Minimizer (top-2 hashes) k=15 w=10", &keys);
    }

    // === 9. (minimizer, position of minimizer) ===
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                let (hash, pos) = canonical_minimizer_with_pos(sb, 15, 10);
                ((hash as u128) << 16) | (pos as u128)
            })
            .collect();
        bench_seq_sort_u128("9. Minimizer (hash, position) k=15 w=10", &keys);
    }

    // === 10. Open syncmers (different selection from closed) ===
    println!();
    for k in [21, 31] {
        let s = k / 2;
        let keys: Vec<u128> = sequences
            .par_iter()
            .map(|seq| compute_min_open_syncmer_hash(seq.as_bytes(), k, s) as u128)
            .collect();
        bench_seq_sort_u128(
            &format!("10. Open syncmer k={} s={}", k, s),
            &keys,
        );
    }

    // === 11. Hybrid: syncmer primary + lex secondary ===
    {
        let keys: Vec<Vec<u8>> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                let hash = compute_min_syncmer_hash(sb, 31, 28);
                let mut key = hash.to_be_bytes().to_vec();
                key.extend_from_slice(sb); // full lex as tiebreak
                key
            })
            .collect();
        bench_seq_sort_bytes("11. Syncmer k=31 + lex tiebreak", &keys);
    }

    // === 12. Hybrid: syncmer primary + RC-canonical lex secondary ===
    {
        let keys: Vec<Vec<u8>> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                let hash = compute_min_syncmer_hash(sb, 31, 28);
                let rc = reverse_complement_bytes(sb);
                let canon = if rc < sb.to_vec() { rc } else { sb.to_vec() };
                let mut key = hash.to_be_bytes().to_vec();
                key.extend_from_slice(&canon);
                key
            })
            .collect();
        bench_seq_sort_bytes("12. Syncmer k=31 + RC-canon lex tiebreak", &keys);
    }

    // === 13-17. Three-level keys: seq-derived + quality tiebreak (all free perm) ===
    println!("\n--- Three-level keys: seq hash + quality tiebreak (both streams, free perm) ---\n");
    println!(
        "{:<55} {:>12} {:>12} {:>12} {:>10}",
        "Strategy", "Seq BSC", "Qual BSC", "Net total", "vs base"
    );
    println!("{}", "-".repeat(104));

    let combined_baseline = seq_baseline + qual_baseline;
    println!(
        "{:<55} {:>10} B {:>10} B {:>10} B",
        "Baseline (no sort)",
        seq_baseline,
        qual_baseline,
        combined_baseline,
    );

    // Helper closure for combined benchmarks
    let bench_combined_u128 = |name: &str, keys: &[u128]| {
        let t = Instant::now();
        let mut perm: Vec<u32> = (0..num_reads as u32).collect();
        perm.sort_unstable_by_key(|&i| keys[i as usize]);
        let elapsed = t.elapsed();
        let seq_reord: Vec<u8> = perm.iter().flat_map(|&i| sequences[i as usize].as_bytes()).copied().collect();
        let qual_reord: Vec<u8> = perm.iter().flat_map(|&i| qualities[i as usize].as_bytes()).copied().collect();
        let seq_comp = bsc::compress_parallel_adaptive(&seq_reord).unwrap();
        let qual_comp = bsc::compress_parallel_adaptive(&qual_reord).unwrap();
        let net = seq_comp.len() + qual_comp.len();
        let delta = net as i64 - combined_baseline as i64;
        println!(
            "{:<55} {:>10} B {:>10} B {:>10} B  {:>+8} B  {:.0}ms",
            name,
            seq_comp.len(),
            qual_comp.len(),
            net,
            delta,
            elapsed.as_secs_f64() * 1000.0,
        );
    };

    // Reference: minimizer (hash, pos) k=15 — current best
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                let (hash, pos) = canonical_minimizer_with_pos(sb, 15, 10);
                ((hash as u128) << 16) | (pos as u128)
            })
            .collect();
        bench_combined_u128("Ref: Minimizer (hash, pos) k=15 w=10", &keys);
    }

    // 13. Three-level: (minimizer_hash, position, mean_quality)
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .zip(qualities.par_iter())
            .map(|(seq, qual)| {
                let sb = seq.as_bytes();
                let (hash, pos) = canonical_minimizer_with_pos(sb, 15, 10);
                let qb = qual.as_bytes();
                let mean_q = qb.iter().map(|&v| v as u64).sum::<u64>() / qb.len().max(1) as u64;
                ((hash as u128) << 24) | ((pos as u128) << 8) | (mean_q as u128)
            })
            .collect();
        bench_combined_u128("13. Minimizer (hash, pos, mean_qual) k=15", &keys);
    }

    // 14. Three-level: (minimizer_hash, mean_quality, position)
    //     Prioritizes quality grouping over position ordering
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .zip(qualities.par_iter())
            .map(|(seq, qual)| {
                let sb = seq.as_bytes();
                let (hash, pos) = canonical_minimizer_with_pos(sb, 15, 10);
                let qb = qual.as_bytes();
                let mean_q = qb.iter().map(|&v| v as u64).sum::<u64>() / qb.len().max(1) as u64;
                ((hash as u128) << 24) | ((mean_q as u128) << 16) | (pos as u128)
            })
            .collect();
        bench_combined_u128("14. Minimizer (hash, mean_qual, pos) k=15", &keys);
    }

    // 15. Three-level with syncmer: (syncmer_hash, position, mean_quality)
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .zip(qualities.par_iter())
            .map(|(seq, qual)| {
                let sb = seq.as_bytes();
                let (hash, pos) = compute_min_syncmer_hash_with_pos(sb, 31, 28);
                let qb = qual.as_bytes();
                let mean_q = qb.iter().map(|&v| v as u64).sum::<u64>() / qb.len().max(1) as u64;
                ((hash as u128) << 24) | ((pos as u128) << 8) | (mean_q as u128)
            })
            .collect();
        bench_combined_u128("15. Syncmer (hash, pos, mean_qual) k=31", &keys);
    }

    // 16. Four-level: (minimizer_hash, pos, mean_qual, qual_variance)
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .zip(qualities.par_iter())
            .map(|(seq, qual)| {
                let sb = seq.as_bytes();
                let (hash, pos) = canonical_minimizer_with_pos(sb, 15, 10);
                let qb = qual.as_bytes();
                let n = qb.len().max(1) as u64;
                let sum: u64 = qb.iter().map(|&v| v as u64).sum();
                let mean_q = sum / n;
                let var: u64 = qb.iter().map(|&v| {
                    let d = (v as i64) - (mean_q as i64);
                    (d * d) as u64
                }).sum::<u64>() / n;
                let var8 = (var.min(255)) as u128;
                ((hash as u128) << 32) | ((pos as u128) << 16) | ((mean_q as u128) << 8) | var8
            })
            .collect();
        bench_combined_u128("16. Minimizer (hash,pos,mean_q,var) k=15", &keys);
    }

    // 17. Open syncmer k=21 + position + mean quality
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .zip(qualities.par_iter())
            .map(|(seq, qual)| {
                let sb = seq.as_bytes();
                let (hash, pos) = compute_min_open_syncmer_hash_with_pos(sb, 21, 10);
                let qb = qual.as_bytes();
                let mean_q = qb.iter().map(|&v| v as u64).sum::<u64>() / qb.len().max(1) as u64;
                ((hash as u128) << 24) | ((pos as u128) << 8) | (mean_q as u128)
            })
            .collect();
        bench_combined_u128("17. Open syncmer (hash,pos,mean_q) k=21", &keys);
    }

    // 18. Minimizer hash + position + quantized quality profile (4 bins)
    {
        let keys: Vec<u128> = sequences
            .par_iter()
            .zip(qualities.par_iter())
            .map(|(seq, qual)| {
                let sb = seq.as_bytes();
                let (hash, pos) = canonical_minimizer_with_pos(sb, 15, 10);
                let qb = qual.as_bytes();
                let profile = quantize_profile(qb, 4);
                ((hash as u128) << 48)
                    | ((pos as u128) << 32)
                    | ((profile[0] as u128) << 24)
                    | ((profile[1] as u128) << 16)
                    | ((profile[2] as u128) << 8)
                    | (profile[3] as u128)
            })
            .collect();
        bench_combined_u128("18. Minimizer (hash,pos,qual_profile4) k=15", &keys);
    }

    println!("{}", "-".repeat(104));

    // ========================================================================
    // QUALITY SORTING STRATEGIES
    // ========================================================================
    println!("\n\n=== QUALITY SORTING STRATEGIES ===\n");
    println!(
        "{:<50} {:>12} {:>8} {:>12} {:>12} {:>10}",
        "Strategy", "Qual BSC", "Ratio", "Perm cost", "Net total", "Time"
    );
    println!("{}", "-".repeat(108));

    // Baseline (qual_baseline already computed at top)
    print_result(
        "1. Baseline (original order)",
        total_qual_bytes,
        qual_baseline,
        0,
        std::time::Duration::ZERO,
    );

    // 2. Lexicographic sort
    bench_sort_strategy(
        "2. Lexicographic",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        qual_baseline,
        |_i, _seq, qual| qual.as_bytes().to_vec(),
        false,
        true,
    );

    // 3. Sort by mean ASCII value
    bench_sort_strategy(
        "3. Mean ASCII value",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        qual_baseline,
        |_i, _seq, qual| {
            let bytes = qual.as_bytes();
            let mean = bytes.iter().map(|&b| b as u64).sum::<u64>() / bytes.len().max(1) as u64;
            vec![mean as u8]
        },
        false,
        true,
    );

    // 4. Sort by L2 norm (sum of squares)
    bench_sort_strategy(
        "4. L2 norm (sum of squares)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        qual_baseline,
        |_i, _seq, qual| {
            let bytes = qual.as_bytes();
            let sos: u64 = bytes.iter().map(|&b| (b as u64) * (b as u64)).sum();
            sos.to_le_bytes().to_vec()
        },
        false,
        true,
    );

    // 5. Multi-key: (mean, variance)
    bench_sort_strategy(
        "5. Multi-key (mean, variance)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        qual_baseline,
        |_i, _seq, qual| {
            let bytes = qual.as_bytes();
            let n = bytes.len().max(1) as u64;
            let sum: u64 = bytes.iter().map(|&b| b as u64).sum();
            let mean = sum / n;
            let var: u64 = bytes.iter().map(|&b| {
                let d = (b as i64) - (mean as i64);
                (d * d) as u64
            }).sum::<u64>() / n;
            let mut key = Vec::with_capacity(3);
            key.push(mean as u8);
            key.extend_from_slice(&(var as u16).to_le_bytes());
            key
        },
        false,
        true,
    );

    // 6. Multi-key: (mean, min)
    bench_sort_strategy(
        "6. Multi-key (mean, min)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        qual_baseline,
        |_i, _seq, qual| {
            let bytes = qual.as_bytes();
            let mean = bytes.iter().map(|&b| b as u64).sum::<u64>() / bytes.len().max(1) as u64;
            let min = bytes.iter().copied().min().unwrap_or(0);
            vec![mean as u8, min]
        },
        false,
        true,
    );

    // 7. Quantized profile (8 bins)
    bench_sort_strategy(
        "7. Quantized profile (8 bins)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        qual_baseline,
        |_i, _seq, qual| quantize_profile(qual.as_bytes(), 8),
        false,
        true,
    );

    // 7b. Quantized profile (16 bins)
    bench_sort_strategy(
        "7b. Quantized profile (16 bins)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        qual_baseline,
        |_i, _seq, qual| quantize_profile(qual.as_bytes(), 16),
        false,
        true,
    );

    // 8. Sort by first + last 16 bytes
    bench_sort_strategy(
        "8. First+last 16 bytes",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        seq_baseline,
        qual_baseline,
        |_i, _seq, qual| {
            let b = qual.as_bytes();
            let n = b.len();
            let mut key = Vec::with_capacity(32);
            key.extend_from_slice(&b[..16.min(n)]);
            if n > 16 {
                key.extend_from_slice(&b[n - 16..]);
            }
            key
        },
        false,
        true,
    );

    // 9. Sequence-derived sort: syncmer hash (free permutation for quality!)
    {
        let anchor_k = 31usize;
        let s = anchor_k.saturating_sub(3).max(5);
        let t_end = anchor_k - s;

        let t = Instant::now();
        let sort_keys: Vec<u64> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                if sb.len() < anchor_k { return u64::MAX; }
                let positions = syncmers::find_syncmers_pos(anchor_k, s, &[0, t_end], sb);
                let mut min_hash = u64::MAX;
                for pos in positions {
                    if pos + anchor_k > sb.len() { continue; }
                    let kmer = &sb[pos..pos + anchor_k];
                    if let Some(fwd) = kmer_to_hash(kmer) {
                        let rc = reverse_complement_hash(fwd, anchor_k);
                        let canon = fwd.min(rc);
                        if canon < min_hash { min_hash = canon; }
                    }
                }
                min_hash
            })
            .collect();
        let mut perm: Vec<u32> = (0..num_reads as u32).collect();
        perm.sort_unstable_by_key(|&i| sort_keys[i as usize]);
        let elapsed = t.elapsed();

        let reordered: Vec<u8> = perm
            .iter()
            .flat_map(|&i| qualities[i as usize].as_bytes())
            .copied()
            .collect();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();

        print_result(
            "9. Syncmer from seq [free perm]",
            total_qual_bytes,
            compressed.len(),
            0,
            elapsed,
        );
        print_delta(compressed.len(), 0, qual_baseline);
    }

    // ========================================================================
    // COMBINED: apply same sort to both streams
    // ========================================================================
    println!("\n\n=== COMBINED: SAME SORT FOR SEQ + QUAL ===\n");
    println!(
        "{:<50} {:>12} {:>12} {:>12} {:>12} {:>10}",
        "Strategy", "Seq BSC", "Qual BSC", "Perm cost", "Net total", "vs base"
    );
    println!("{}", "-".repeat(112));

    let combined_baseline = seq_baseline + qual_baseline;
    println!(
        "{:<50} {:>10} B {:>10} B {:>10}   {:>10} B {:>10}",
        "Baseline (no sort)",
        seq_baseline,
        qual_baseline,
        "-",
        combined_baseline,
        "-",
    );

    // Lex sort on sequences, apply same permutation to qualities
    bench_combined(
        "Lex sort (on sequences)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        combined_baseline,
        |_i, seq, _qual| seq.as_bytes().to_vec(),
    );

    // Mean ASCII sort on sequences
    bench_combined(
        "Mean ASCII (on sequences)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        combined_baseline,
        |_i, seq, _qual| {
            let bytes = seq.as_bytes();
            let mean = bytes.iter().map(|&b| b as u64).sum::<u64>() / bytes.len().max(1) as u64;
            vec![mean as u8]
        },
    );

    // Quantized profile (8 bins) on sequences
    bench_combined(
        "Quantized profile 8-bin (on sequences)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        combined_baseline,
        |_i, seq, _qual| quantize_profile(seq.as_bytes(), 8),
    );

    // Quantized profile (16 bins) on sequences
    bench_combined(
        "Quantized profile 16-bin (on sequences)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        combined_baseline,
        |_i, seq, _qual| quantize_profile(seq.as_bytes(), 16),
    );

    // Lex sort on qualities, apply same permutation to sequences
    bench_combined(
        "Lex sort (on qualities)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        combined_baseline,
        |_i, _seq, qual| qual.as_bytes().to_vec(),
    );

    // Mean quality sort, apply to both
    bench_combined(
        "Mean ASCII (on qualities)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        combined_baseline,
        |_i, _seq, qual| {
            let bytes = qual.as_bytes();
            let mean = bytes.iter().map(|&b| b as u64).sum::<u64>() / bytes.len().max(1) as u64;
            vec![mean as u8]
        },
    );

    // Multi-key (mean, variance) on qualities
    bench_combined(
        "Multi-key mean+var (on qualities)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        combined_baseline,
        |_i, _seq, qual| {
            let bytes = qual.as_bytes();
            let n = bytes.len().max(1) as u64;
            let sum: u64 = bytes.iter().map(|&b| b as u64).sum();
            let mean = sum / n;
            let var: u64 = bytes.iter().map(|&b| {
                let d = (b as i64) - (mean as i64);
                (d * d) as u64
            }).sum::<u64>() / n;
            let mut key = Vec::with_capacity(3);
            key.push(mean as u8);
            key.extend_from_slice(&(var as u16).to_le_bytes());
            key
        },
    );

    // RC-canonical + lex on sequences
    bench_combined(
        "RC-canonical + lex (on sequences)",
        &sequences,
        &qualities,
        total_seq_bytes,
        total_qual_bytes,
        combined_baseline,
        |_i, seq, _qual| {
            let b = seq.as_bytes();
            let rc = reverse_complement_bytes(b);
            if rc < b.to_vec() { rc } else { b.to_vec() }
        },
    );

    // Minimizer sort on sequences
    {
        let t = Instant::now();
        let min_hashes: Vec<u64> = sequences
            .par_iter()
            .map(|s| canonical_minimizer(s.as_bytes(), 15, 10))
            .collect();
        let mut perm: Vec<u32> = (0..num_reads as u32).collect();
        perm.sort_unstable_by_key(|&i| min_hashes[i as usize]);
        let elapsed = t.elapsed();

        let seq_reord: Vec<u8> = perm.iter().flat_map(|&i| sequences[i as usize].as_bytes()).copied().collect();
        let qual_reord: Vec<u8> = perm.iter().flat_map(|&i| qualities[i as usize].as_bytes()).copied().collect();

        let seq_comp = bsc::compress_parallel_adaptive(&seq_reord).unwrap();
        let qual_comp = bsc::compress_parallel_adaptive(&qual_reord).unwrap();

        let net = seq_comp.len() + qual_comp.len();
        let delta = net as i64 - combined_baseline as i64;
        println!(
            "{:<50} {:>10} B {:>10} B {:>10}   {:>10} B {:>+9} B  {:.0}ms",
            "Minimizer (k=15,w=10) [free perm]",
            seq_comp.len(),
            qual_comp.len(),
            0,
            net,
            delta,
            elapsed.as_secs_f64() * 1000.0,
        );
    }

    // Syncmer sort on sequences
    {
        let anchor_k = 31usize;
        let s = anchor_k.saturating_sub(3).max(5);
        let t_end = anchor_k - s;

        let t = Instant::now();
        let sort_keys: Vec<u64> = sequences
            .par_iter()
            .map(|seq| {
                let sb = seq.as_bytes();
                if sb.len() < anchor_k { return u64::MAX; }
                let positions = syncmers::find_syncmers_pos(anchor_k, s, &[0, t_end], sb);
                let mut min_hash = u64::MAX;
                for pos in positions {
                    if pos + anchor_k > sb.len() { continue; }
                    let kmer = &sb[pos..pos + anchor_k];
                    if let Some(fwd) = kmer_to_hash(kmer) {
                        let rc = reverse_complement_hash(fwd, anchor_k);
                        let canon = fwd.min(rc);
                        if canon < min_hash { min_hash = canon; }
                    }
                }
                min_hash
            })
            .collect();
        let mut perm: Vec<u32> = (0..num_reads as u32).collect();
        perm.sort_unstable_by_key(|&i| sort_keys[i as usize]);
        let elapsed = t.elapsed();

        let seq_reord: Vec<u8> = perm.iter().flat_map(|&i| sequences[i as usize].as_bytes()).copied().collect();
        let qual_reord: Vec<u8> = perm.iter().flat_map(|&i| qualities[i as usize].as_bytes()).copied().collect();

        let seq_comp = bsc::compress_parallel_adaptive(&seq_reord).unwrap();
        let qual_comp = bsc::compress_parallel_adaptive(&qual_reord).unwrap();

        let net = seq_comp.len() + qual_comp.len();
        let delta = net as i64 - combined_baseline as i64;
        println!(
            "{:<50} {:>10} B {:>10} B {:>10}   {:>10} B {:>+9} B  {:.0}ms",
            "Syncmer (k=31,s=28) [free perm]",
            seq_comp.len(),
            qual_comp.len(),
            0,
            net,
            delta,
            elapsed.as_secs_f64() * 1000.0,
        );
    }

    println!("{}", "-".repeat(112));

    // ========================================================================
    // PERMUTATION ENCODING COMPARISON
    // ========================================================================
    // Use lex-sort on qualities (worst permutation cost) as the test case.
    println!("\n\n=== PERMUTATION ENCODING METHODS (lex-sort on qualities) ===\n");
    {
        let mut perm_lex: Vec<u32> = (0..num_reads as u32).collect();
        perm_lex.sort_by(|&a, &b| qualities[a as usize].cmp(&qualities[b as usize]));

        // Reorder qualities (same for all encoding methods)
        let reordered: Vec<u8> = perm_lex
            .iter()
            .flat_map(|&i| qualities[i as usize].as_bytes())
            .copied()
            .collect();
        let qual_comp = bsc::compress_parallel_adaptive(&reordered).unwrap();

        println!(
            "{:<55} {:>12} {:>12} {:>12}",
            "Encoding method", "Perm size", "Net total", "vs baseline"
        );
        println!("{}", "-".repeat(95));
        println!(
            "{:<55} {:>10}   {:>10} B {:>+10} B",
            "Quality baseline (no sort)",
            "-",
            qual_baseline,
            0,
        );
        println!(
            "{:<55} {:>10}   {:>10} B",
            "Sorted quality stream (all methods share this)",
            "-",
            qual_comp.len(),
        );
        println!();

        // Method A: Forward perm, delta-zigzag-varint + BSC (current)
        let fwd_cost = compress_perm_delta_zigzag(&perm_lex);
        println!(
            "{:<55} {:>10} B {:>10} B {:>+10} B",
            "A. Forward delta-zigzag-varint + BSC",
            fwd_cost,
            qual_comp.len() + fwd_cost,
            (qual_comp.len() + fwd_cost) as i64 - qual_baseline as i64,
        );

        // Method B: Inverse perm, delta-zigzag-varint + BSC
        let inv_perm = invert_permutation(&perm_lex);
        let inv_cost = compress_perm_delta_zigzag(&inv_perm);
        println!(
            "{:<55} {:>10} B {:>10} B {:>+10} B",
            "B. Inverse delta-zigzag-varint + BSC",
            inv_cost,
            qual_comp.len() + inv_cost,
            (qual_comp.len() + inv_cost) as i64 - qual_baseline as i64,
        );

        // Method C: Raw u32 + BSC (no delta encoding)
        let raw_perm: Vec<u8> = perm_lex.iter().flat_map(|&i| i.to_le_bytes()).collect();
        let raw_cost = bsc::compress_parallel_adaptive(&raw_perm).unwrap().len();
        println!(
            "{:<55} {:>10} B {:>10} B {:>+10} B",
            "C. Raw u32 + BSC (no delta)",
            raw_cost,
            qual_comp.len() + raw_cost,
            (qual_comp.len() + raw_cost) as i64 - qual_baseline as i64,
        );

        // Method D: Raw inverse u32 + BSC
        let raw_inv: Vec<u8> = inv_perm.iter().flat_map(|&i| i.to_le_bytes()).collect();
        let raw_inv_cost = bsc::compress_parallel_adaptive(&raw_inv).unwrap().len();
        println!(
            "{:<55} {:>10} B {:>10} B {:>+10} B",
            "D. Raw inverse u32 + BSC",
            raw_inv_cost,
            qual_comp.len() + raw_inv_cost,
            (qual_comp.len() + raw_inv_cost) as i64 - qual_baseline as i64,
        );

        // Method E: Cycle encoding + BSC
        let cycle_cost = compress_perm_cycles(&perm_lex);
        println!(
            "{:<55} {:>10} B {:>10} B {:>+10} B",
            "E. Cycle encoding + BSC",
            cycle_cost,
            qual_comp.len() + cycle_cost,
            (qual_comp.len() + cycle_cost) as i64 - qual_baseline as i64,
        );

        // Method F: Sort key (1B mean qual) + BSC — recomputable!
        let mean_keys: Vec<u8> = qualities
            .iter()
            .map(|q| {
                let b = q.as_bytes();
                (b.iter().map(|&v| v as u64).sum::<u64>() / b.len().max(1) as u64) as u8
            })
            .collect();
        let key_cost = bsc::compress_parallel_adaptive(&mean_keys).unwrap().len();

        // Sort by mean key and measure quality compression
        let mut perm_mean: Vec<u32> = (0..num_reads as u32).collect();
        perm_mean.sort_by_key(|&i| mean_keys[i as usize]);
        let mean_reord: Vec<u8> = perm_mean
            .iter()
            .flat_map(|&i| qualities[i as usize].as_bytes())
            .copied()
            .collect();
        let mean_qual_comp = bsc::compress_parallel_adaptive(&mean_reord).unwrap();
        println!(
            "{:<55} {:>10} B {:>10} B {:>+10} B",
            "F. Recomputable key (1B mean) + BSC [mean sort]",
            key_cost,
            mean_qual_comp.len() + key_cost,
            (mean_qual_comp.len() + key_cost) as i64 - qual_baseline as i64,
        );

        // Method G: Sort key (3B mean+var) + BSC — recomputable!
        let mv_keys: Vec<u8> = qualities
            .iter()
            .flat_map(|q| {
                let b = q.as_bytes();
                let n = b.len().max(1) as u64;
                let sum: u64 = b.iter().map(|&v| v as u64).sum();
                let mean = sum / n;
                let var: u64 = b.iter().map(|&v| {
                    let d = (v as i64) - (mean as i64);
                    (d * d) as u64
                }).sum::<u64>() / n;
                [mean as u8, (var & 0xFF) as u8, ((var >> 8) & 0xFF) as u8]
            })
            .collect();
        let mv_cost = bsc::compress_parallel_adaptive(&mv_keys).unwrap().len();

        let mut perm_mv: Vec<u32> = (0..num_reads as u32).collect();
        perm_mv.sort_by(|&a, &b| {
            let ak = &mv_keys[a as usize * 3..a as usize * 3 + 3];
            let bk = &mv_keys[b as usize * 3..b as usize * 3 + 3];
            ak.cmp(bk)
        });
        let mv_reord: Vec<u8> = perm_mv
            .iter()
            .flat_map(|&i| qualities[i as usize].as_bytes())
            .copied()
            .collect();
        let mv_qual_comp = bsc::compress_parallel_adaptive(&mv_reord).unwrap();
        println!(
            "{:<55} {:>10} B {:>10} B {:>+10} B",
            "G. Recomputable key (3B mean+var) + BSC [mv sort]",
            mv_cost,
            mv_qual_comp.len() + mv_cost,
            (mv_qual_comp.len() + mv_cost) as i64 - qual_baseline as i64,
        );
    }

    // ========================================================================
    // BLOCK-BASED REORDERING
    // ========================================================================
    println!("\n\n=== BLOCK-BASED REORDERING (qualities) ===");
    println!("Sort blocks of reads together instead of individual reads.\n");
    {
        println!(
            "{:<50} {:>8} {:>12} {:>12} {:>12} {:>12}",
            "Strategy", "BlkSize", "Qual BSC", "Perm cost", "Net total", "vs baseline"
        );
        println!("{}", "-".repeat(110));
        println!(
            "{:<50} {:>8} {:>10} B {:>10}   {:>10} B {:>10}",
            "Baseline (no sort)", "-", qual_baseline, "-", qual_baseline, "-",
        );

        // Test block-based reordering for several strategies x block sizes
        let block_sizes = [32, 64, 128, 256, 512, 1024, 4096];

        // Strategy: Lex sort on qualities
        for &bs in &block_sizes {
            bench_block_sort(
                "Lex-sort (qualities)",
                &qualities,
                total_qual_bytes,
                qual_baseline,
                bs,
                |qual| qual.as_bytes().to_vec(),
            );
        }
        println!();

        // Strategy: Mean quality sort
        for &bs in &block_sizes {
            bench_block_sort(
                "Mean-sort (qualities)",
                &qualities,
                total_qual_bytes,
                qual_baseline,
                bs,
                |qual| {
                    let b = qual.as_bytes();
                    let mean = b.iter().map(|&v| v as u64).sum::<u64>() / b.len().max(1) as u64;
                    vec![mean as u8]
                },
            );
        }
        println!();

        // Strategy: Quantized profile 16-bin
        for &bs in &block_sizes {
            bench_block_sort(
                "Quantized 16-bin (qualities)",
                &qualities,
                total_qual_bytes,
                qual_baseline,
                bs,
                |qual| quantize_profile(qual.as_bytes(), 16),
            );
        }
        println!();

        // Strategy: Multi-key (mean, variance)
        for &bs in &block_sizes {
            bench_block_sort(
                "Mean+var (qualities)",
                &qualities,
                total_qual_bytes,
                qual_baseline,
                bs,
                |qual| {
                    let b = qual.as_bytes();
                    let n = b.len().max(1) as u64;
                    let sum: u64 = b.iter().map(|&v| v as u64).sum();
                    let mean = sum / n;
                    let var: u64 = b.iter().map(|&v| {
                        let d = (v as i64) - (mean as i64);
                        (d * d) as u64
                    }).sum::<u64>() / n;
                    let mut key = Vec::with_capacity(3);
                    key.push(mean as u8);
                    key.extend_from_slice(&(var as u16).to_le_bytes());
                    key
                },
            );
        }
    }

    // ========================================================================
    // BLOCK-BASED REORDERING (sequences)
    // ========================================================================
    println!("\n\n=== BLOCK-BASED REORDERING (sequences) ===\n");
    {
        println!(
            "{:<50} {:>8} {:>12} {:>12} {:>12} {:>12}",
            "Strategy", "BlkSize", "Seq BSC", "Perm cost", "Net total", "vs baseline"
        );
        println!("{}", "-".repeat(110));
        println!(
            "{:<50} {:>8} {:>10} B {:>10}   {:>10} B {:>10}",
            "Baseline (no sort)", "-", seq_baseline, "-", seq_baseline, "-",
        );

        let block_sizes = [32, 64, 128, 256, 512, 1024, 4096];

        // Strategy: Lex sort on sequences
        for &bs in &block_sizes {
            bench_block_sort(
                "Lex-sort (sequences)",
                &sequences,
                total_seq_bytes,
                seq_baseline,
                bs,
                |seq| seq.as_bytes().to_vec(),
            );
        }
        println!();

        // Strategy: Mean ASCII sort on sequences
        for &bs in &block_sizes {
            bench_block_sort(
                "Mean-ASCII (sequences)",
                &sequences,
                total_seq_bytes,
                seq_baseline,
                bs,
                |seq| {
                    let b = seq.as_bytes();
                    let mean = b.iter().map(|&v| v as u64).sum::<u64>() / b.len().max(1) as u64;
                    vec![mean as u8]
                },
            );
        }
        println!();

        // Strategy: Quantized profile 8-bin on sequences
        for &bs in &block_sizes {
            bench_block_sort(
                "Quantized 8-bin (sequences)",
                &sequences,
                total_seq_bytes,
                seq_baseline,
                bs,
                |seq| quantize_profile(seq.as_bytes(), 8),
            );
        }
    }

    println!("\n{}", "-".repeat(112));
    println!(
        "Input: {} reads x {}bp — seq {:.1} MB, qual {:.1} MB, combined baseline {} B",
        num_reads,
        read_len,
        total_seq_bytes as f64 / (1024.0 * 1024.0),
        total_qual_bytes as f64 / (1024.0 * 1024.0),
        combined_baseline,
    );
}

// ============================================================================
// Helpers
// ============================================================================

/// Benchmark a sorting strategy on either sequences or qualities (or both).
/// `sort_key_fn(index, seq, qual) -> Vec<u8>` produces the sort key for each read.
fn bench_sort_strategy<F>(
    name: &str,
    sequences: &[String],
    qualities: &[String],
    total_seq_bytes: usize,
    total_qual_bytes: usize,
    seq_baseline: usize,
    qual_baseline: usize,
    sort_key_fn: F,
    do_seq: bool,
    do_qual: bool,
) where
    F: Fn(usize, &String, &String) -> Vec<u8> + Sync,
{
    let num_reads = sequences.len();

    let t = Instant::now();
    // Compute sort keys
    let sort_keys: Vec<Vec<u8>> = (0..num_reads)
        .into_par_iter()
        .map(|i| sort_key_fn(i, &sequences[i], &qualities[i]))
        .collect();

    // Build permutation
    let mut perm: Vec<u32> = (0..num_reads as u32).collect();
    perm.sort_by(|&a, &b| sort_keys[a as usize].cmp(&sort_keys[b as usize]));
    let elapsed = t.elapsed();

    // Permutation cost: delta-zigzag-varint + BSC
    let perm_cost = compress_permutation(&perm);

    if do_seq {
        let reordered: Vec<u8> = perm
            .iter()
            .flat_map(|&i| sequences[i as usize].as_bytes())
            .copied()
            .collect();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();

        print_result(name, total_seq_bytes, compressed.len(), perm_cost, elapsed);
        print_delta(compressed.len(), perm_cost, seq_baseline);
    }

    if do_qual {
        let reordered: Vec<u8> = perm
            .iter()
            .flat_map(|&i| qualities[i as usize].as_bytes())
            .copied()
            .collect();
        let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();

        print_result(name, total_qual_bytes, compressed.len(), perm_cost, elapsed);
        print_delta(compressed.len(), perm_cost, qual_baseline);
    }
}

/// Benchmark a sorting strategy applied to BOTH seq and qual with the same permutation.
fn bench_combined<F>(
    name: &str,
    sequences: &[String],
    qualities: &[String],
    _total_seq_bytes: usize,
    _total_qual_bytes: usize,
    combined_baseline: usize,
    sort_key_fn: F,
) where
    F: Fn(usize, &String, &String) -> Vec<u8> + Sync,
{
    let num_reads = sequences.len();

    let t = Instant::now();
    let sort_keys: Vec<Vec<u8>> = (0..num_reads)
        .into_par_iter()
        .map(|i| sort_key_fn(i, &sequences[i], &qualities[i]))
        .collect();

    let mut perm: Vec<u32> = (0..num_reads as u32).collect();
    perm.sort_by(|&a, &b| sort_keys[a as usize].cmp(&sort_keys[b as usize]));
    let elapsed = t.elapsed();

    let perm_cost = compress_permutation(&perm);

    let seq_reord: Vec<u8> = perm
        .iter()
        .flat_map(|&i| sequences[i as usize].as_bytes())
        .copied()
        .collect();
    let qual_reord: Vec<u8> = perm
        .iter()
        .flat_map(|&i| qualities[i as usize].as_bytes())
        .copied()
        .collect();

    let seq_comp = bsc::compress_parallel_adaptive(&seq_reord).unwrap();
    let qual_comp = bsc::compress_parallel_adaptive(&qual_reord).unwrap();

    let net = seq_comp.len() + qual_comp.len() + perm_cost;
    let delta = net as i64 - combined_baseline as i64;

    println!(
        "{:<50} {:>10} B {:>10} B {:>10} B {:>10} B {:>+9} B  {:.0}ms",
        name,
        seq_comp.len(),
        qual_comp.len(),
        perm_cost,
        net,
        delta,
        elapsed.as_secs_f64() * 1000.0,
    );
}

/// Compress a permutation using delta-zigzag-varint encoding + BSC.
/// Returns compressed size in bytes.
fn compress_permutation(perm: &[u32]) -> usize {
    let mut delta_buf: Vec<u8> = Vec::with_capacity(perm.len() * 3);
    let mut prev = 0i64;
    for &idx in perm {
        let delta = idx as i64 - prev;
        prev = idx as i64;
        // Zigzag encode
        let zz = ((delta << 1) ^ (delta >> 63)) as u64;
        // Varint encode
        let mut v = zz;
        while v >= 0x80 {
            delta_buf.push((v as u8) | 0x80);
            v >>= 7;
        }
        delta_buf.push(v as u8);
    }
    let compressed = bsc::compress_parallel_adaptive(&delta_buf).unwrap();
    compressed.len()
}

/// Quantize a byte vector into `num_bins` bins of mean values.
fn quantize_profile(data: &[u8], num_bins: usize) -> Vec<u8> {
    if data.is_empty() {
        return vec![0; num_bins];
    }
    let bin_size = data.len().max(1) / num_bins.max(1);
    let mut profile = Vec::with_capacity(num_bins);
    for bin in 0..num_bins {
        let start = bin * data.len() / num_bins;
        let end = ((bin + 1) * data.len() / num_bins).min(data.len());
        if start >= end {
            profile.push(0);
        } else {
            let sum: u64 = data[start..end].iter().map(|&b| b as u64).sum();
            let _ = bin_size; // suppress unused warning
            profile.push((sum / (end - start) as u64) as u8);
        }
    }
    profile
}

/// Reverse complement a DNA byte slice.
fn reverse_complement_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}

/// Compute canonical minimizer for a sequence.
fn canonical_minimizer(seq: &[u8], k: usize, w: usize) -> u64 {
    if seq.len() < k {
        return u64::MAX;
    }

    let mut min_hash = u64::MAX;
    for start in 0..seq.len().saturating_sub(k + w - 1).max(1) {
        let end = (start + w).min(seq.len().saturating_sub(k - 1));
        for pos in start..end {
            let kmer = &seq[pos..pos + k];
            if let Some(fwd) = kmer_to_hash(kmer) {
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                if canon < min_hash {
                    min_hash = canon;
                }
            }
        }
    }
    min_hash
}

fn print_result(
    name: &str,
    original: usize,
    compressed: usize,
    perm_cost: usize,
    elapsed: std::time::Duration,
) {
    let net = compressed + perm_cost;
    let ratio = original as f64 / compressed as f64;
    if perm_cost > 0 {
        println!(
            "{:<50} {:>10} B {:>6.2}x {:>10} B {:>10} B {:>8.0}ms",
            name,
            compressed,
            ratio,
            perm_cost,
            net,
            elapsed.as_secs_f64() * 1000.0,
        );
    } else {
        println!(
            "{:<50} {:>10} B {:>6.2}x {:>10}   {:>10} B {:>8.0}ms",
            name,
            compressed,
            ratio,
            "-",
            net,
            elapsed.as_secs_f64() * 1000.0,
        );
    }
}

// ============================================================================
// Syncmer / minimizer helper functions
// ============================================================================

/// Compute minimum canonical syncmer hash for a sequence (closed syncmers).
fn compute_min_syncmer_hash(seq: &[u8], k: usize, s: usize) -> u64 {
    if seq.len() < k { return u64::MAX; }
    let t_end = k - s;
    let positions = syncmers::find_syncmers_pos(k, s, &[0, t_end], seq);
    let mut min_hash = u64::MAX;
    for pos in positions {
        if pos + k > seq.len() { continue; }
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            let canon = fwd.min(rc);
            if canon < min_hash { min_hash = canon; }
        }
    }
    min_hash
}

/// Compute minimum canonical open syncmer hash with position.
fn compute_min_open_syncmer_hash_with_pos(seq: &[u8], k: usize, s: usize) -> (u64, u16) {
    if seq.len() < k { return (u64::MAX, 0); }
    let positions = syncmers::find_syncmers_pos(k, s, &[0], seq);
    let mut min_hash = u64::MAX;
    let mut min_pos = 0u16;
    for pos in positions {
        if pos + k > seq.len() { continue; }
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            let canon = fwd.min(rc);
            if canon < min_hash {
                min_hash = canon;
                min_pos = pos as u16;
            }
        }
    }
    (min_hash, min_pos)
}

/// Compute minimum canonical open syncmer hash.
/// Open syncmer: smallest s-mer is at position 0 only (not at end).
fn compute_min_open_syncmer_hash(seq: &[u8], k: usize, s: usize) -> u64 {
    if seq.len() < k { return u64::MAX; }
    let positions = syncmers::find_syncmers_pos(k, s, &[0], seq);
    let mut min_hash = u64::MAX;
    for pos in positions {
        if pos + k > seq.len() { continue; }
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            let canon = fwd.min(rc);
            if canon < min_hash { min_hash = canon; }
        }
    }
    min_hash
}

/// Compute top-2 smallest canonical syncmer hashes.
fn compute_top2_syncmer_hashes(seq: &[u8], k: usize, s: usize) -> (u64, u64) {
    if seq.len() < k { return (u64::MAX, u64::MAX); }
    let t_end = k - s;
    let positions = syncmers::find_syncmers_pos(k, s, &[0, t_end], seq);
    let mut h1 = u64::MAX;
    let mut h2 = u64::MAX;
    for pos in positions {
        if pos + k > seq.len() { continue; }
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            let canon = fwd.min(rc);
            if canon < h1 {
                h2 = h1;
                h1 = canon;
            } else if canon < h2 && canon != h1 {
                h2 = canon;
            }
        }
    }
    (h1, h2)
}

/// Compute minimum canonical syncmer hash and its position in the read.
fn compute_min_syncmer_hash_with_pos(seq: &[u8], k: usize, s: usize) -> (u64, u16) {
    if seq.len() < k { return (u64::MAX, 0); }
    let t_end = k - s;
    let positions = syncmers::find_syncmers_pos(k, s, &[0, t_end], seq);
    let mut min_hash = u64::MAX;
    let mut min_pos = 0u16;
    for pos in positions {
        if pos + k > seq.len() { continue; }
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            let canon = fwd.min(rc);
            if canon < min_hash {
                min_hash = canon;
                min_pos = pos as u16;
            }
        }
    }
    (min_hash, min_pos)
}

/// Compute sorted top-N canonical syncmer hashes (MinHash sketch).
fn compute_sorted_syncmer_hashes(seq: &[u8], k: usize, s: usize, n: usize) -> Vec<u64> {
    if seq.len() < k { return vec![u64::MAX; n]; }
    let t_end = k - s;
    let positions = syncmers::find_syncmers_pos(k, s, &[0, t_end], seq);
    let mut hashes: Vec<u64> = Vec::with_capacity(positions.len());
    for pos in positions {
        if pos + k > seq.len() { continue; }
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            hashes.push(fwd.min(rc));
        }
    }
    hashes.sort_unstable();
    hashes.dedup();
    hashes.truncate(n);
    while hashes.len() < n {
        hashes.push(u64::MAX);
    }
    hashes
}

/// Compute top-2 smallest canonical minimizer hashes.
fn compute_top2_minimizer_hashes(seq: &[u8], k: usize, w: usize) -> (u64, u64) {
    if seq.len() < k { return (u64::MAX, u64::MAX); }
    let mut h1 = u64::MAX;
    let mut h2 = u64::MAX;

    for pos in 0..=seq.len().saturating_sub(k) {
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            let canon = fwd.min(rc);
            if canon < h1 {
                h2 = h1;
                h1 = canon;
            } else if canon < h2 && canon != h1 {
                h2 = canon;
            }
        }
    }
    let _ = w; // window parameter not used in top-2 global min
    (h1, h2)
}

/// Compute canonical minimizer with its position.
fn canonical_minimizer_with_pos(seq: &[u8], k: usize, w: usize) -> (u64, u16) {
    if seq.len() < k { return (u64::MAX, 0); }
    let mut min_hash = u64::MAX;
    let mut min_pos = 0u16;

    for start in 0..seq.len().saturating_sub(k + w - 1).max(1) {
        let end = (start + w).min(seq.len().saturating_sub(k - 1));
        for pos in start..end {
            let kmer = &seq[pos..pos + k];
            if let Some(fwd) = kmer_to_hash(kmer) {
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                if canon < min_hash {
                    min_hash = canon;
                    min_pos = pos as u16;
                }
            }
        }
    }
    (min_hash, min_pos)
}

fn print_delta(compressed: usize, perm_cost: usize, baseline: usize) {
    let net = compressed + perm_cost;
    let delta = net as i64 - baseline as i64;
    println!(
        "{:<50} {:>73}",
        "",
        format!("({:+} B vs baseline)", delta),
    );
}

/// Invert a permutation: if perm[sorted_idx] = original_idx,
/// then inv[original_idx] = sorted_idx.
fn invert_permutation(perm: &[u32]) -> Vec<u32> {
    let mut inv = vec![0u32; perm.len()];
    for (sorted_idx, &original_idx) in perm.iter().enumerate() {
        inv[original_idx as usize] = sorted_idx as u32;
    }
    inv
}

/// Compress permutation using delta-zigzag-varint + BSC (generic version).
fn compress_perm_delta_zigzag(perm: &[u32]) -> usize {
    let mut buf: Vec<u8> = Vec::with_capacity(perm.len() * 3);
    let mut prev = 0i64;
    for &idx in perm {
        let delta = idx as i64 - prev;
        prev = idx as i64;
        let zz = ((delta << 1) ^ (delta >> 63)) as u64;
        let mut v = zz;
        while v >= 0x80 {
            buf.push((v as u8) | 0x80);
            v >>= 7;
        }
        buf.push(v as u8);
    }
    bsc::compress_parallel_adaptive(&buf).unwrap().len()
}

/// Encode permutation as cycles, then compress with BSC.
/// Cycle representation: [cycle_len, elem0, elem1, ...] for each cycle.
fn compress_perm_cycles(perm: &[u32]) -> usize {
    let n = perm.len();
    let mut visited = vec![false; n];
    let mut cycle_buf: Vec<u8> = Vec::with_capacity(n * 5);

    for start in 0..n {
        if visited[start] { continue; }
        // Collect cycle
        let mut cycle = Vec::new();
        let mut cur = start;
        while !visited[cur] {
            visited[cur] = true;
            cycle.push(cur as u32);
            cur = perm[cur] as usize;
        }
        // Write cycle length as varint
        let mut clen = cycle.len() as u64;
        while clen >= 0x80 {
            cycle_buf.push((clen as u8) | 0x80);
            clen >>= 7;
        }
        cycle_buf.push(clen as u8);
        // Write cycle elements as delta-zigzag-varints
        let mut prev = 0i64;
        for &elem in &cycle {
            let delta = elem as i64 - prev;
            prev = elem as i64;
            let zz = ((delta << 1) ^ (delta >> 63)) as u64;
            let mut v = zz;
            while v >= 0x80 {
                cycle_buf.push((v as u8) | 0x80);
                v >>= 7;
            }
            cycle_buf.push(v as u8);
        }
    }
    bsc::compress_parallel_adaptive(&cycle_buf).unwrap().len()
}

/// Block-based sort benchmark.
/// Sorts blocks of `block_size` reads by the sort key of the block's first read,
/// then reorders data block-by-block (within-block order preserved).
fn bench_block_sort<F>(
    name: &str,
    data: &[String],
    total_bytes: usize,
    baseline: usize,
    block_size: usize,
    sort_key_fn: F,
) where
    F: Fn(&String) -> Vec<u8>,
{
    let num_reads = data.len();
    let num_blocks = (num_reads + block_size - 1) / block_size;

    // Compute sort key per block (use mean of sort keys within block, or just first read's key)
    let mut block_keys: Vec<(Vec<u8>, usize)> = (0..num_blocks)
        .map(|b| {
            let start = b * block_size;
            let end = (start + block_size).min(num_reads);
            // Use the median-ish key: take the sort key of the middle read in the block
            // Actually, for simplicity and speed, use the mean of all keys in the block
            // But keys are Vec<u8>... just use first read's key for sorting blocks
            let key = sort_key_fn(&data[start]);
            // Also try: aggregate key = sort key of element with median value in block
            let _ = end;
            (key, b)
        })
        .collect();
    block_keys.sort_by(|a, b| a.0.cmp(&b.0));

    // Reorder data block by block
    let mut reordered = Vec::with_capacity(total_bytes);
    for &(_, block_idx) in &block_keys {
        let start = block_idx * block_size;
        let end = (start + block_size).min(num_reads);
        for i in start..end {
            reordered.extend_from_slice(data[i].as_bytes());
        }
    }

    let compressed = bsc::compress_parallel_adaptive(&reordered).unwrap();

    // Block permutation cost: just num_blocks u16/u32 values + BSC
    let perm_cost = if num_blocks <= 65536 {
        let block_perm: Vec<u8> = block_keys
            .iter()
            .flat_map(|&(_, idx)| (idx as u16).to_le_bytes())
            .collect();
        bsc::compress_parallel_adaptive(&block_perm).unwrap().len()
    } else {
        let block_perm: Vec<u8> = block_keys
            .iter()
            .flat_map(|&(_, idx)| (idx as u32).to_le_bytes())
            .collect();
        bsc::compress_parallel_adaptive(&block_perm).unwrap().len()
    };

    let net = compressed.len() + perm_cost;
    let delta = net as i64 - baseline as i64;
    println!(
        "{:<50} {:>8} {:>10} B {:>10} B {:>10} B {:>+10} B",
        name,
        block_size,
        compressed.len(),
        perm_cost,
        net,
        delta,
    );
}
