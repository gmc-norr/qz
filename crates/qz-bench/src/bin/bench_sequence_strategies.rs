/// Benchmark: sequence compression strategies
///
/// Tests various approaches to compressing DNA sequences from FASTQ,
/// focusing on strategies from the coil codebase:
///
/// 1. Raw ASCII + BSC (current QZ default)
/// 2. 2-bit packed + BSC
/// 3. Columnar (position-transposed) + BSC
/// 4. De Bruijn graph encoding
/// 5. Arithmetic coding (order-3 Markov context)
/// 6. Columnar 2-bit packed + BSC

use std::io::BufRead;
use std::time::Instant;

fn main() {
    let fastq_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/ERR3239334_1.5m.fastq".to_string());

    eprintln!("Reading FASTQ: {}", fastq_path);
    let file = std::fs::File::open(&fastq_path).expect("Cannot open FASTQ file");
    let reader = std::io::BufReader::new(file);

    let mut sequences: Vec<String> = Vec::new();
    let mut line_num = 0u64;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        match line_num % 4 {
            1 => sequences.push(line.trim_end().to_string()),
            _ => {}
        }
        line_num += 1;
    }

    let num_reads = sequences.len();
    let total_bases: usize = sequences.iter().map(|s| s.len()).sum();
    let first_len = sequences[0].len();
    let all_same_len = sequences.iter().all(|s| s.len() == first_len);

    eprintln!(
        "Loaded {} reads, {} bases ({:.1} MB), read_len={}{}\n",
        num_reads,
        total_bases,
        total_bases as f64 / (1024.0 * 1024.0),
        first_len,
        if all_same_len { " (uniform)" } else { " (variable)" }
    );

    println!("{:<40} {:>12} {:>10} {:>12} {:>10}", "Method", "Compressed", "Ratio", "Comp Time", "Decomp");
    println!("{}", "-".repeat(86));

    // === 1. Raw ASCII + BSC (current QZ default) ===
    {
        let raw: Vec<u8> = sequences.iter().flat_map(|s| s.as_bytes().iter().copied()).collect();
        let t = Instant::now();
        let compressed = qz_lib::compression::bsc::compress_parallel_adaptive(&raw).unwrap();
        let comp_time = t.elapsed();

        let t2 = Instant::now();
        let decompressed = qz_lib::compression::bsc::decompress_parallel(&compressed).unwrap();
        let decomp_time = t2.elapsed();
        assert_eq!(raw, decompressed);

        print_row("1. Raw ASCII + BSC", total_bases, compressed.len(), comp_time, decomp_time);
    }

    // === 2. 2-bit packed + BSC ===
    {
        let seqs_bytes: Vec<Vec<u8>> = sequences.iter().map(|s| s.as_bytes().to_vec()).collect();
        let t = Instant::now();
        let packed = qz_lib::compression::dna_utils::pack_dna_2bit(&seqs_bytes);
        let pack_time = t.elapsed();

        let t = Instant::now();
        let compressed = qz_lib::compression::bsc::compress_parallel_adaptive(&packed).unwrap();
        let comp_time = t.elapsed();

        let t2 = Instant::now();
        let decompressed = qz_lib::compression::bsc::decompress_parallel(&compressed).unwrap();
        let decomp_time = t2.elapsed();

        // Verify roundtrip
        let (unpacked, _) = qz_lib::compression::dna_utils::unpack_dna_2bit(&decompressed, 0).unwrap();
        for (i, (orig, dec)) in sequences.iter().zip(unpacked.iter()).enumerate() {
            assert_eq!(orig.as_bytes(), dec.as_slice(), "2-bit roundtrip mismatch at read {}", i);
        }

        print_row_extra(
            "2. 2-bit packed + BSC",
            total_bases,
            compressed.len(),
            pack_time + comp_time,
            decomp_time,
            &format!("pack: {:.1}ms", pack_time.as_secs_f64() * 1000.0),
        );
    }

    // === 3. Columnar (position-transposed) ASCII + BSC ===
    // Instead of read-by-read, store all position-0 bases, then all position-1 bases, etc.
    // This groups positional base distributions together for BWT.
    if all_same_len {
        let t = Instant::now();
        let columnar = transpose_sequences(&sequences, first_len);
        let transpose_time = t.elapsed();

        let t = Instant::now();
        let compressed = qz_lib::compression::bsc::compress_parallel_adaptive(&columnar).unwrap();
        let comp_time = t.elapsed();

        let t2 = Instant::now();
        let decompressed = qz_lib::compression::bsc::decompress_parallel(&compressed).unwrap();
        let decomp_time = t2.elapsed();

        // Verify roundtrip
        let reconstructed = untranspose_sequences(&decompressed, num_reads, first_len);
        for (i, (orig, dec)) in sequences.iter().zip(reconstructed.iter()).enumerate() {
            assert_eq!(orig, dec, "Columnar roundtrip mismatch at read {}", i);
        }

        print_row_extra(
            "3. Columnar ASCII + BSC",
            total_bases,
            compressed.len(),
            transpose_time + comp_time,
            decomp_time,
            &format!("transpose: {:.1}ms", transpose_time.as_secs_f64() * 1000.0),
        );
    }

    // === 4. Columnar 2-bit packed + BSC ===
    // Transpose, then 2-bit encode each column (same-position bases together)
    if all_same_len {
        let t = Instant::now();
        let columnar = transpose_sequences(&sequences, first_len);
        // Pack the transposed data: each column of num_reads bases
        let packed = pack_2bit_raw(&columnar);
        let pack_time = t.elapsed();

        let t = Instant::now();
        let compressed = qz_lib::compression::bsc::compress_parallel_adaptive(&packed).unwrap();
        let comp_time = t.elapsed();

        let t2 = Instant::now();
        let decompressed = qz_lib::compression::bsc::decompress_parallel(&compressed).unwrap();
        let decomp_time = t2.elapsed();

        // Verify roundtrip (N bases map to A — lossy, N's go in separate stream in QZ)
        let unpacked = unpack_2bit_raw(&decompressed, columnar.len());
        let reconstructed = untranspose_sequences(&unpacked, num_reads, first_len);
        let mut n_mismatches = 0usize;
        for (orig, dec) in sequences.iter().zip(reconstructed.iter()) {
            for (&a, &b) in orig.as_bytes().iter().zip(dec.as_bytes().iter()) {
                if a != b {
                    assert!(a == b'N' || a == b'n', "Non-N mismatch: {} vs {}", a as char, b as char);
                    n_mismatches += 1;
                }
            }
        }
        if n_mismatches > 0 {
            eprintln!("  (columnar 2-bit: {} N bases mapped to A — handled by N-mask stream)", n_mismatches);
        }

        print_row_extra(
            "4. Columnar 2-bit + BSC",
            total_bases,
            compressed.len(),
            pack_time + comp_time,
            decomp_time,
            &format!("pack: {:.1}ms", pack_time.as_secs_f64() * 1000.0),
        );
    }

    // === 5. Greedy contig builder + BSC ===
    {
        let t = Instant::now();
        let compressed = qz_lib::compression::greedy_contig::compress_sequences_greedy(&sequences).unwrap();
        let comp_time = t.elapsed();

        // Decompress uses same DBG1 format as de Bruijn
        let t2 = Instant::now();
        let decompressed = qz_lib::compression::debruijn::decompress_sequences_debruijn(&compressed, num_reads).unwrap();
        let decomp_time = t2.elapsed();

        // Verify roundtrip
        for (i, (orig, dec)) in sequences.iter().zip(decompressed.iter()).enumerate() {
            assert_eq!(orig, dec, "Greedy contig roundtrip mismatch at read {}", i);
        }

        print_row("5. Greedy contig + BSC", total_bases, compressed.len(), comp_time, decomp_time);
    }

    // === 5b. De Bruijn graph encoding (skip if > 200K reads — too slow) ===
    if num_reads <= 200_000 {
        let t = Instant::now();
        let compressed = qz_lib::compression::debruijn::compress_sequences_debruijn(&sequences, 0).unwrap();
        let comp_time = t.elapsed();

        let t2 = Instant::now();
        let decompressed = qz_lib::compression::debruijn::decompress_sequences_debruijn(&compressed, num_reads).unwrap();
        let decomp_time = t2.elapsed();

        for (i, (orig, dec)) in sequences.iter().zip(decompressed.iter()).enumerate() {
            assert_eq!(orig, dec, "De Bruijn roundtrip mismatch at read {}", i);
        }

        print_row("5b. De Bruijn graph + BSC", total_bases, compressed.len(), comp_time, decomp_time);
    } else {
        println!("{:<40} {:>10} {:>10} {:>12} {:>10}", "5b. De Bruijn graph + BSC", "SKIPPED", "", "(too slow)", "");
    }

    // === 6. Arithmetic coding (order-3 Markov) ===
    // Only run on a subset if dataset is large (arithmetic coding is O(n) but slow per-symbol)
    {
        let max_reads_arith = 100_000; // Limit for arithmetic coding
        let subset: Vec<String> = sequences.iter().take(max_reads_arith).cloned().collect();
        let subset_bases: usize = subset.iter().map(|s| s.len()).sum();

        let t = Instant::now();
        let compressed = qz_lib::compression::arithmetic_sequence::encode_sequences_arithmetic(&subset).unwrap();
        let comp_time = t.elapsed();

        let read_lengths: Vec<usize> = subset.iter().map(|s| s.len()).collect();
        let t2 = Instant::now();
        let decompressed = qz_lib::compression::arithmetic_sequence::decode_sequences_arithmetic(
            &compressed,
            &read_lengths,
            subset.len(),
        ).unwrap();
        let decomp_time = t2.elapsed();

        for (i, (orig, dec)) in subset.iter().zip(decompressed.iter()).enumerate() {
            assert_eq!(orig, dec, "Arithmetic roundtrip mismatch at read {}", i);
        }

        let label = if subset.len() < num_reads {
            format!("6. Arithmetic ({}K subset)", subset.len() / 1000)
        } else {
            "6. Arithmetic coding".to_string()
        };
        print_row(&label, subset_bases, compressed.len(), comp_time, decomp_time);
    }

    // === 7. Lex-sorted + BSC (for reference) ===
    {
        let t = Instant::now();
        let mut sorted = sequences.clone();
        sorted.sort_unstable();
        let sorted_stream: Vec<u8> = sorted.iter().flat_map(|s| s.as_bytes().iter().copied()).collect();
        let sort_time = t.elapsed();

        let t = Instant::now();
        let compressed = qz_lib::compression::bsc::compress_parallel_adaptive(&sorted_stream).unwrap();
        let comp_time = t.elapsed();

        let t2 = Instant::now();
        let _ = qz_lib::compression::bsc::decompress_parallel(&compressed).unwrap();
        let decomp_time = t2.elapsed();

        print_row_extra(
            "7. Lex-sorted ASCII + BSC",
            total_bases,
            compressed.len(),
            sort_time + comp_time,
            decomp_time,
            &format!("sort: {:.1}ms", sort_time.as_secs_f64() * 1000.0),
        );
    }

    // === 8. Minimizer-sorted + BSC ===
    // Sort reads by their canonical minimizer (approximates genomic position)
    {
        let t = Instant::now();
        let mut indexed: Vec<(u64, usize)> = sequences
            .iter()
            .enumerate()
            .map(|(i, s)| (canonical_minimizer(s.as_bytes(), 15, 10), i))
            .collect();
        indexed.sort_unstable_by_key(|&(min_hash, _)| min_hash);

        let sorted_stream: Vec<u8> = indexed
            .iter()
            .flat_map(|&(_, idx)| sequences[idx].as_bytes().iter().copied())
            .collect();
        let sort_time = t.elapsed();

        let t = Instant::now();
        let compressed = qz_lib::compression::bsc::compress_parallel_adaptive(&sorted_stream).unwrap();
        let comp_time = t.elapsed();

        let t2 = Instant::now();
        let _ = qz_lib::compression::bsc::decompress_parallel(&compressed).unwrap();
        let decomp_time = t2.elapsed();

        print_row_extra(
            "8. Minimizer-sorted + BSC",
            total_bases,
            compressed.len(),
            sort_time + comp_time,
            decomp_time,
            &format!("sort: {:.1}ms", sort_time.as_secs_f64() * 1000.0),
        );
    }

    // === 9. Minimizer-sorted + columnar + BSC ===
    if all_same_len {
        let t = Instant::now();
        let mut indexed: Vec<(u64, usize)> = sequences
            .iter()
            .enumerate()
            .map(|(i, s)| (canonical_minimizer(s.as_bytes(), 15, 10), i))
            .collect();
        indexed.sort_unstable_by_key(|&(min_hash, _)| min_hash);

        let sorted_seqs: Vec<String> = indexed.iter().map(|&(_, idx)| sequences[idx].clone()).collect();
        let columnar = transpose_sequences(&sorted_seqs, first_len);
        let sort_time = t.elapsed();

        let t = Instant::now();
        let compressed = qz_lib::compression::bsc::compress_parallel_adaptive(&columnar).unwrap();
        let comp_time = t.elapsed();

        let t2 = Instant::now();
        let _ = qz_lib::compression::bsc::decompress_parallel(&compressed).unwrap();
        let decomp_time = t2.elapsed();

        print_row_extra(
            "9. Minimizer-sort + columnar + BSC",
            total_bases,
            compressed.len(),
            sort_time + comp_time,
            decomp_time,
            &format!("sort+transpose: {:.1}ms", sort_time.as_secs_f64() * 1000.0),
        );
    }

    println!("\n{}", "-".repeat(86));
    println!("Input: {} reads x {} bp = {} bytes ({:.1} MB)",
        num_reads, first_len, total_bases, total_bases as f64 / (1024.0 * 1024.0));
    println!("Theoretical 2-bit minimum: {} bytes ({:.2}x)",
        total_bases / 4, 4.0);
}

/// Transpose sequences: read-by-read → position-by-position
/// Output: all bases at position 0, then all at position 1, etc.
fn transpose_sequences(sequences: &[String], read_len: usize) -> Vec<u8> {
    let num_reads = sequences.len();
    let mut transposed = Vec::with_capacity(num_reads * read_len);

    for pos in 0..read_len {
        for seq in sequences {
            transposed.push(seq.as_bytes()[pos]);
        }
    }

    transposed
}

/// Reverse transpose: position-by-position → read-by-read
fn untranspose_sequences(data: &[u8], num_reads: usize, read_len: usize) -> Vec<String> {
    let mut sequences = vec![Vec::with_capacity(read_len); num_reads];

    for pos in 0..read_len {
        for (read_idx, seq) in sequences.iter_mut().enumerate() {
            seq.push(data[pos * num_reads + read_idx]);
        }
    }

    sequences
        .into_iter()
        .map(|s| String::from_utf8(s).unwrap())
        .collect()
}

/// Simple 2-bit packing of raw bytes (no varint framing, for columnar data)
fn pack_2bit_raw(data: &[u8]) -> Vec<u8> {
    let mut packed = Vec::with_capacity(data.len() / 4 + 1);
    for chunk in data.chunks(4) {
        let mut byte = 0u8;
        for (i, &base) in chunk.iter().enumerate() {
            let val = match base {
                b'A' | b'a' => 0u8,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0, // N → A (lossy for N, but N is in separate stream)
            };
            byte |= val << (2 * i);
        }
        packed.push(byte);
    }
    packed
}

/// Unpack 2-bit raw data back to ASCII bases
fn unpack_2bit_raw(packed: &[u8], expected_len: usize) -> Vec<u8> {
    let mut unpacked = Vec::with_capacity(expected_len);
    for &byte in packed {
        for i in 0..4 {
            if unpacked.len() >= expected_len {
                break;
            }
            let val = (byte >> (2 * i)) & 0b11;
            let base = match val {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => unreachable!(),
            };
            unpacked.push(base);
        }
    }
    unpacked
}

/// Compute canonical minimizer for a sequence.
/// Returns the minimum canonical k-mer hash in sliding windows of size w.
fn canonical_minimizer(seq: &[u8], k: usize, w: usize) -> u64 {
    if seq.len() < k {
        return u64::MAX;
    }

    let mut min_hash = u64::MAX;

    for start in 0..seq.len().saturating_sub(k + w - 1).max(1) {
        let end = (start + w).min(seq.len().saturating_sub(k - 1));
        for pos in start..end {
            let kmer = &seq[pos..pos + k];
            if let Some(fwd) = qz_lib::compression::dna_utils::kmer_to_hash(kmer) {
                let rc = qz_lib::compression::dna_utils::reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                if canon < min_hash {
                    min_hash = canon;
                }
            }
        }
    }

    min_hash
}

fn print_row(name: &str, original: usize, compressed: usize, comp: std::time::Duration, decomp: std::time::Duration) {
    let ratio = original as f64 / compressed as f64;
    println!(
        "{:<40} {:>10} B {:>8.2}x {:>10.1}ms {:>8.1}ms",
        name, compressed, ratio,
        comp.as_secs_f64() * 1000.0,
        decomp.as_secs_f64() * 1000.0,
    );
}

fn print_row_extra(name: &str, original: usize, compressed: usize, comp: std::time::Duration, decomp: std::time::Duration, extra: &str) {
    let ratio = original as f64 / compressed as f64;
    println!(
        "{:<40} {:>10} B {:>8.2}x {:>10.1}ms {:>8.1}ms  ({})",
        name, compressed, ratio,
        comp.as_secs_f64() * 1000.0,
        decomp.as_secs_f64() * 1000.0,
        extra,
    );
}
