/// Prototype: Transforms orthogonal to BWT that might assist BSC
///
/// BWT already captures context-dependent base distributions. These strategies
/// try to exploit structure that BWT *misses*:
///
///   A. Baseline: raw sequences + BSC
///   B. 2-bit packing: 4 bases per byte, reduces input 4x. BWT sees 4-grams.
///   C. RC canonicalization: store lexicographic min(read, revcomp(read)).
///      Doubles effective coverage per genomic region.
///   D. Mod-4 delta (sorted): (base[i] - prev_read_base[i]) mod 4.
///      Same 4-symbol alphabet, but heavily 0-biased after sorting.
///   E. Sorted + 2-bit packing
///   F. RC canon + sorted + 2-bit packing
///   G. Sorted + mod-4 delta + 2-bit packing
///   H. Small-cluster interleave: groups of N sorted reads interleaved by position.
///      Partial transposition that preserves BWT-friendly locality.
///
/// Usage: cargo run --release --bin bench_bsc_assists [fastq_path]

use std::io::BufRead;
use std::time::Instant;

use qz_lib::compression::bsc;
use rayon::prelude::*;

// ── Base encoding ────────────────────────────────────────────────────────

#[inline]
fn base_to_2bit(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

#[inline]
fn bit2_to_base(v: u8) -> u8 {
    match v & 3 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

// ── 2-bit packing: 4 bases per byte ──────────────────────────────────────

fn pack_2bit(bases: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity((bases.len() + 3) / 4);
    for chunk in bases.chunks(4) {
        let mut byte = 0u8;
        for (j, &b) in chunk.iter().enumerate() {
            byte |= base_to_2bit(b) << (6 - 2 * j);
        }
        out.push(byte);
    }
    out
}

fn unpack_2bit(packed: &[u8], num_bases: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(num_bases);
    for (i, &byte) in packed.iter().enumerate() {
        for j in 0..4 {
            if i * 4 + j >= num_bases {
                break;
            }
            out.push(bit2_to_base((byte >> (6 - 2 * j)) & 3));
        }
    }
    out
}

// ── Reverse complement ──────────────────────────────────────────────────

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Canonicalize: return (canonical_seq, was_reversed)
fn canonicalize(seq: &[u8]) -> (Vec<u8>, bool) {
    let rc = reverse_complement(seq);
    if rc < seq.to_vec() {
        (rc, true)
    } else {
        (seq.to_vec(), false)
    }
}

// ── Mod-4 delta encoding ────────────────────────────────────────────────

/// Encode read as (current - previous) mod 4 per base position.
/// Stays in {0,1,2,3} alphabet. For similar reads, mostly 0.
fn mod4_delta(current: &[u8], previous: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(current.len());
    let overlap = current.len().min(previous.len());
    for i in 0..overlap {
        let c = base_to_2bit(current[i]);
        let p = base_to_2bit(previous[i]);
        out.push((c + 4 - p) % 4); // mod-4 difference
    }
    // Remaining bases (if current longer): raw 2-bit
    for i in overlap..current.len() {
        out.push(base_to_2bit(current[i]));
    }
    out
}

/// Decode mod-4 delta back to bases.
fn mod4_undelta(delta: &[u8], previous: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(delta.len());
    let overlap = delta.len().min(previous.len());
    for i in 0..overlap {
        let p = base_to_2bit(previous[i]);
        let base = (delta[i] + p) % 4;
        out.push(bit2_to_base(base));
    }
    for i in overlap..delta.len() {
        out.push(bit2_to_base(delta[i]));
    }
    out
}

// ── Minimizer sort key ──────────────────────────────────────────────────

fn minimizer_sort_key(seq: &[u8]) -> u64 {
    let k = 21;
    if seq.len() < k {
        return u64::MAX;
    }
    let mut min_hash = u64::MAX;
    for i in 0..=seq.len() - k {
        let mut fwd: u64 = 0;
        let mut rev: u64 = 0;
        let mut valid = true;
        for j in 0..k {
            let v = match seq[i + j] {
                b'A' | b'a' => 0u64,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => {
                    valid = false;
                    break;
                }
            };
            fwd = (fwd << 2) | v;
            rev = rev | ((3 - v) << (2 * j));
        }
        if valid {
            min_hash = min_hash.min(fwd.min(rev));
        }
    }
    min_hash
}

// ── Helpers ──────────────────────────────────────────────────────────────

fn report(name: &str, compressed_size: usize, raw_size: usize, extra_cost: usize, elapsed: f64) {
    let total = compressed_size + extra_cost;
    if extra_cost > 0 {
        eprintln!(
            "  {:<42} {:>10} B  ({:.3}x, {:.4} bpb)  net {:>10} B ({:.3}x)  {:.2}s",
            name,
            compressed_size,
            raw_size as f64 / compressed_size as f64,
            8.0 * compressed_size as f64 / raw_size as f64,
            total,
            raw_size as f64 / total as f64,
            elapsed,
        );
    } else {
        eprintln!(
            "  {:<42} {:>10} B  ({:.3}x, {:.4} bpb)  {:.2}s",
            name,
            compressed_size,
            raw_size as f64 / compressed_size as f64,
            8.0 * compressed_size as f64 / raw_size as f64,
            elapsed,
        );
    }
}

// ── Main ────────────────────────────────────────────────────────────────

fn main() {
    let fastq_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/ERR3239334_1.1m.fastq".to_string());

    eprintln!("=== BSC-Assist Transforms Benchmark ===");
    eprintln!("Reading FASTQ: {}", fastq_path);

    let t0 = Instant::now();
    let file = std::fs::File::open(&fastq_path).expect("Cannot open FASTQ file");
    let reader = std::io::BufReader::new(file);

    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut line_num = 0u64;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line_num % 4 == 1 {
            sequences.push(line.trim_end().as_bytes().to_vec());
        }
        line_num += 1;
    }
    eprintln!(
        "Read {} sequences in {:.2}s",
        sequences.len(),
        t0.elapsed().as_secs_f64()
    );

    let num_reads = sequences.len();
    let read_len = sequences[0].len();
    let raw_size = num_reads * read_len;
    eprintln!("Read length: {}, raw size: {} bytes\n", read_len, raw_size);

    // ── A. Baseline ──────────────────────────────────────────────────────

    let t = Instant::now();
    let raw_concat: Vec<u8> = sequences.iter().flat_map(|s| s.iter().copied()).collect();
    let baseline_size = bsc::compress_parallel_adaptive(&raw_concat).unwrap().len();
    report("A. Baseline (raw + BSC)", baseline_size, raw_size, 0, t.elapsed().as_secs_f64());
    drop(raw_concat);

    // ── B. 2-bit packing ────────────────────────────────────────────────

    let t = Instant::now();
    let packed: Vec<u8> = sequences
        .par_iter()
        .flat_map_iter(|seq| {
            // Pack each read independently (aligned to read boundary)
            pack_2bit(seq)
        })
        .collect();
    let packed_size = bsc::compress_parallel_adaptive(&packed).unwrap().len();
    report("B. 2-bit packed + BSC", packed_size, raw_size, 0, t.elapsed().as_secs_f64());

    // Verify roundtrip
    {
        let packed_one = pack_2bit(&sequences[0]);
        let unpacked = unpack_2bit(&packed_one, read_len);
        for i in 0..read_len {
            let expected = match sequences[0][i] {
                b'A' | b'a' => b'A', b'C' | b'c' => b'C',
                b'G' | b'g' => b'G', b'T' | b't' => b'T', _ => b'A',
            };
            assert_eq!(unpacked[i], expected, "2-bit roundtrip failed at pos {}", i);
        }
    }
    drop(packed);

    // ── B2. 2-bit packing (whole-stream, not per-read) ──────────────────

    let t = Instant::now();
    let all_2bit: Vec<u8> = sequences
        .iter()
        .flat_map(|seq| seq.iter().map(|&b| base_to_2bit(b)))
        .collect();
    let packed_stream = pack_2bit_from_2bit(&all_2bit);
    let packed_stream_size = bsc::compress_parallel_adaptive(&packed_stream).unwrap().len();
    report("B2. 2-bit whole-stream + BSC", packed_stream_size, raw_size, 0, t.elapsed().as_secs_f64());
    drop(packed_stream);

    // ── C. RC canonicalization ───────────────────────────────────────────

    let t = Instant::now();
    let canon: Vec<(Vec<u8>, bool)> = sequences.par_iter().map(|s| canonicalize(s)).collect();
    let rc_flags: Vec<u8> = canon.iter().map(|(_, rev)| *rev as u8).collect();
    let canon_concat: Vec<u8> = canon.iter().flat_map(|(s, _)| s.iter().copied()).collect();
    let canon_size = bsc::compress_parallel_adaptive(&canon_concat).unwrap().len();
    let flags_size = bsc::compress_adaptive(&rc_flags).unwrap().len();
    report(
        "C. RC canonical + BSC",
        canon_size + flags_size,
        raw_size,
        0,
        t.elapsed().as_secs_f64(),
    );
    eprintln!("     (canon={} B, flags={} B, rc_frac={:.1}%)",
        canon_size, flags_size,
        100.0 * rc_flags.iter().filter(|&&f| f == 1).count() as f64 / num_reads as f64);
    drop(canon_concat);

    // ── Sort reads ──────────────────────────────────────────────────────

    eprintln!();
    let t = Instant::now();
    let mut sort_keys: Vec<(u64, usize)> = sequences
        .par_iter()
        .enumerate()
        .map(|(i, seq)| (minimizer_sort_key(seq), i))
        .collect();
    sort_keys.sort_unstable();
    let sorted_idx: Vec<usize> = sort_keys.iter().map(|&(_, i)| i).collect();
    let sort_time = t.elapsed().as_secs_f64();
    eprintln!("  Sorted by minimizer in {:.2}s", sort_time);

    // Permutation cost
    let perm_bytes: Vec<u8> = sorted_idx
        .iter()
        .flat_map(|&i| (i as u32).to_le_bytes())
        .collect();
    let perm_cost = bsc::compress_parallel_adaptive(&perm_bytes).unwrap().len();
    eprintln!("  Permutation cost: {} B\n", perm_cost);

    let sorted_seqs: Vec<&[u8]> = sorted_idx.iter().map(|&i| sequences[i].as_slice()).collect();

    // Sorted baseline
    let t = Instant::now();
    let sorted_concat: Vec<u8> = sorted_seqs.iter().flat_map(|s| s.iter().copied()).collect();
    let sorted_size = bsc::compress_parallel_adaptive(&sorted_concat).unwrap().len();
    report("D0. Sorted + BSC", sorted_size, raw_size, perm_cost, t.elapsed().as_secs_f64());
    drop(sorted_concat);

    // ── D. Sorted + mod-4 delta ─────────────────────────────────────────

    let t = Instant::now();
    // Encode: first read raw 2-bit, subsequent as mod-4 delta vs previous
    let mut delta_stream = Vec::with_capacity(raw_size);
    for (i, seq) in sorted_seqs.iter().enumerate() {
        if i == 0 {
            for &b in seq.iter() {
                delta_stream.push(base_to_2bit(b));
            }
        } else {
            let prev = sorted_seqs[i - 1];
            let delta = mod4_delta(seq, prev);
            delta_stream.extend_from_slice(&delta);
        }
    }

    // Check distribution
    let mut dist = [0u64; 4];
    for &v in &delta_stream {
        dist[v as usize] += 1;
    }
    let total_bases = delta_stream.len() as f64;
    eprintln!(
        "  Mod-4 delta distribution: 0={:.1}% 1={:.1}% 2={:.1}% 3={:.1}%",
        100.0 * dist[0] as f64 / total_bases,
        100.0 * dist[1] as f64 / total_bases,
        100.0 * dist[2] as f64 / total_bases,
        100.0 * dist[3] as f64 / total_bases,
    );

    let delta_size = bsc::compress_parallel_adaptive(&delta_stream).unwrap().len();
    report("D. Sorted + mod4-delta + BSC", delta_size, raw_size, perm_cost, t.elapsed().as_secs_f64());

    // Verify roundtrip for first 100 reads
    {
        let mut prev: Vec<u8> = Vec::new();
        let mut pos = 0;
        for i in 0..100.min(num_reads) {
            let rl = sorted_seqs[i].len();
            if i == 0 {
                for j in 0..rl {
                    let recovered = bit2_to_base(delta_stream[pos + j]);
                    let expected = match sorted_seqs[i][j] {
                        b'A' | b'a' => b'A', b'C' | b'c' => b'C',
                        b'G' | b'g' => b'G', b'T' | b't' => b'T', _ => b'A',
                    };
                    assert_eq!(recovered, expected, "Delta roundtrip failed read {} pos {}", i, j);
                }
                prev = sorted_seqs[i].to_vec();
            } else {
                let decoded = mod4_undelta(&delta_stream[pos..pos + rl], &prev);
                for j in 0..rl {
                    let expected = match sorted_seqs[i][j] {
                        b'A' | b'a' => b'A', b'C' | b'c' => b'C',
                        b'G' | b'g' => b'G', b'T' | b't' => b'T', _ => b'A',
                    };
                    assert_eq!(decoded[j], expected, "Delta roundtrip failed read {} pos {}", i, j);
                }
                prev = sorted_seqs[i].to_vec();
            }
            pos += rl;
        }
        eprintln!("  Delta roundtrip verified (first 100 reads)");
    }
    drop(delta_stream);

    // ── E. Sorted + 2-bit packing ───────────────────────────────────────

    let t = Instant::now();
    let sorted_2bit: Vec<u8> = sorted_seqs
        .iter()
        .flat_map(|seq| seq.iter().map(|&b| base_to_2bit(b)))
        .collect();
    let sorted_packed = pack_2bit_from_2bit(&sorted_2bit);
    let sorted_packed_size = bsc::compress_parallel_adaptive(&sorted_packed).unwrap().len();
    report("E. Sorted + 2-bit packed + BSC", sorted_packed_size, raw_size, perm_cost, t.elapsed().as_secs_f64());
    drop(sorted_packed);

    // ── F. RC canon + sorted ────────────────────────────────────────────

    let t = Instant::now();
    let canon_seqs: Vec<(Vec<u8>, bool)> = sequences.par_iter().map(|s| canonicalize(s)).collect();
    let mut canon_sort_keys: Vec<(u64, usize)> = canon_seqs
        .par_iter()
        .enumerate()
        .map(|(i, (seq, _))| (minimizer_sort_key(seq), i))
        .collect();
    canon_sort_keys.sort_unstable();
    let canon_sorted_idx: Vec<usize> = canon_sort_keys.iter().map(|&(_, i)| i).collect();

    let canon_sorted_concat: Vec<u8> = canon_sorted_idx
        .iter()
        .flat_map(|&i| canon_seqs[i].0.iter().copied())
        .collect();
    let canon_sorted_size = bsc::compress_parallel_adaptive(&canon_sorted_concat).unwrap().len();

    // RC flags in sorted order
    let canon_sorted_flags: Vec<u8> = canon_sorted_idx.iter().map(|&i| canon_seqs[i].1 as u8).collect();
    let canon_flags_size = bsc::compress_adaptive(&canon_sorted_flags).unwrap().len();

    // Permutation
    let canon_perm_bytes: Vec<u8> = canon_sorted_idx
        .iter()
        .flat_map(|&i| (i as u32).to_le_bytes())
        .collect();
    let canon_perm_cost = bsc::compress_parallel_adaptive(&canon_perm_bytes).unwrap().len();

    let canon_total = canon_sorted_size + canon_flags_size;
    report(
        "F. RC canon + sorted + BSC",
        canon_total,
        raw_size,
        canon_perm_cost,
        t.elapsed().as_secs_f64(),
    );
    eprintln!("     (seqs={} B, flags={} B, perm={} B)", canon_sorted_size, canon_flags_size, canon_perm_cost);
    drop(canon_sorted_concat);

    // ── G. Sorted + mod4-delta + 2-bit packed ───────────────────────────

    let t = Instant::now();
    let mut delta_2bit = Vec::with_capacity(raw_size);
    for (i, seq) in sorted_seqs.iter().enumerate() {
        if i == 0 {
            for &b in seq.iter() {
                delta_2bit.push(base_to_2bit(b));
            }
        } else {
            let prev = sorted_seqs[i - 1];
            let delta = mod4_delta(seq, prev);
            delta_2bit.extend_from_slice(&delta);
        }
    }
    let delta_packed = pack_2bit_from_2bit(&delta_2bit);
    let delta_packed_size = bsc::compress_parallel_adaptive(&delta_packed).unwrap().len();
    report("G. Sorted + mod4-delta + 2bit + BSC", delta_packed_size, raw_size, perm_cost, t.elapsed().as_secs_f64());
    drop(delta_packed);

    // ── H. Small-cluster interleave ─────────────────────────────────────
    // Group N sorted reads, interleave by position:
    // [r0[0], r1[0], ..., rN-1[0], r0[1], r1[1], ..., rN-1[1], ...]
    // Tests partial transposition that keeps BWT-friendly locality.

    for cluster_n in [2, 4, 8, 16, 32] {
        let t = Instant::now();
        let mut interleaved = Vec::with_capacity(raw_size);

        for chunk in sorted_seqs.chunks(cluster_n) {
            let cn = chunk.len();
            let rl = chunk[0].len();
            for col in 0..rl {
                for row in 0..cn {
                    if col < chunk[row].len() {
                        interleaved.push(chunk[row][col]);
                    }
                }
            }
        }

        let interleaved_size = bsc::compress_parallel_adaptive(&interleaved).unwrap().len();
        report(
            &format!("H. Sorted + interleave(N={}) + BSC", cluster_n),
            interleaved_size,
            raw_size,
            perm_cost,
            t.elapsed().as_secs_f64(),
        );
    }

    // ── I. Sorted + mod4-delta + cluster interleave ─────────────────────
    // Combine: first delta encode, then interleave within clusters.

    let best_n = 4; // try a promising cluster size
    {
        let t = Instant::now();
        // Build delta-encoded sorted sequences
        let mut delta_seqs: Vec<Vec<u8>> = Vec::with_capacity(num_reads);
        for (i, seq) in sorted_seqs.iter().enumerate() {
            if i == 0 {
                delta_seqs.push(seq.iter().map(|&b| base_to_2bit(b)).collect());
            } else {
                delta_seqs.push(mod4_delta(seq, sorted_seqs[i - 1]));
            }
        }

        // Interleave within clusters
        let mut interleaved = Vec::with_capacity(raw_size);
        for chunk in delta_seqs.chunks(best_n) {
            let cn = chunk.len();
            let rl = chunk[0].len();
            for col in 0..rl {
                for row in 0..cn {
                    if col < chunk[row].len() {
                        interleaved.push(chunk[row][col]);
                    }
                }
            }
        }

        let size = bsc::compress_parallel_adaptive(&interleaved).unwrap().len();
        report(
            &format!("I. Sorted + mod4-delta + interleave({}) + BSC", best_n),
            size,
            raw_size,
            perm_cost,
            t.elapsed().as_secs_f64(),
        );
    }

    // ── J. 2-bit stream (no packing, just 2-bit values) + BSC ───────────
    // Test whether BWT prefers {0,1,2,3} over {A,C,G,T} (ASCII)

    let t = Instant::now();
    let twobit_stream: Vec<u8> = sequences
        .iter()
        .flat_map(|seq| seq.iter().map(|&b| base_to_2bit(b)))
        .collect();
    let twobit_size = bsc::compress_parallel_adaptive(&twobit_stream).unwrap().len();
    report("J. 2-bit values (not packed) + BSC", twobit_size, raw_size, 0, t.elapsed().as_secs_f64());
    drop(twobit_stream);

    // ── K. Sorted + 2-bit values (not packed) + BSC ─────────────────────

    let t = Instant::now();
    let sorted_twobit: Vec<u8> = sorted_seqs
        .iter()
        .flat_map(|seq| seq.iter().map(|&b| base_to_2bit(b)))
        .collect();
    let sorted_twobit_size = bsc::compress_parallel_adaptive(&sorted_twobit).unwrap().len();
    report("K. Sorted + 2-bit values + BSC", sorted_twobit_size, raw_size, perm_cost, t.elapsed().as_secs_f64());

    // ── Summary ─────────────────────────────────────────────────────────

    eprintln!("\n=== Summary (sorted by compressed size) ===");
    let mut results: Vec<(&str, usize, usize)> = vec![
        ("A. Baseline", baseline_size, 0),
        ("B. 2-bit packed", packed_size, 0),
        ("B2. 2-bit whole-stream", packed_stream_size, 0),
        ("C. RC canon", canon_size + flags_size, 0),
        ("D0. Sorted", sorted_size, perm_cost),
        ("D. Sorted + mod4-delta", delta_size, perm_cost),
        ("E. Sorted + 2-bit packed", sorted_packed_size, perm_cost),
        ("F. RC canon + sorted", canon_total, canon_perm_cost),
        ("G. Sorted + delta + 2bit", delta_packed_size, perm_cost),
        ("J. 2-bit values", twobit_size, 0),
        ("K. Sorted + 2-bit values", sorted_twobit_size, perm_cost),
    ];
    results.sort_by_key(|&(_, size, _)| size);

    for (name, size, extra) in &results {
        let total = size + extra;
        if *extra > 0 {
            eprintln!(
                "  {:<32} {:>10} B  {:.3}x  (net {:.3}x with {} B perm)",
                name, size, raw_size as f64 / *size as f64,
                raw_size as f64 / total as f64, extra,
            );
        } else {
            eprintln!(
                "  {:<32} {:>10} B  {:.3}x",
                name, size, raw_size as f64 / *size as f64,
            );
        }
    }
}

/// Pack a stream of 2-bit values (0-3) into bytes, 4 values per byte.
fn pack_2bit_from_2bit(values: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity((values.len() + 3) / 4);
    for chunk in values.chunks(4) {
        let mut byte = 0u8;
        for (j, &v) in chunk.iter().enumerate() {
            byte |= (v & 3) << (6 - 2 * j);
        }
        out.push(byte);
    }
    out
}
