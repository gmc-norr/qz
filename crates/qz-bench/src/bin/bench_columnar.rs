/// Prototype: Columnar compression of reordered + shift-aligned reads
///
/// Pipeline:
///   1. Read FASTQ sequences
///   2. Sort by minimizer hash (groups similar reads)
///   3. Within sorted order, shift-align adjacent reads using shared k-mers
///   4. Build an aligned matrix with gap markers
///   5. Compute per-column consensus
///   6. Encode residual matrix (0 = matches consensus, 1-3 = alt base)
///   7. Transpose (column-major) the residual matrix
///   8. BSC compress
///
/// Compares several strategies:
///   A. Baseline: raw sequences + BSC (no reorder)
///   B. Sorted: reordered sequences + BSC (row-major, no alignment)
///   C. Sorted + delta: delta vs previous read + BSC
///   D. Sorted + columnar: transpose sorted matrix + BSC
///   E. Sorted + consensus-residual: residual matrix + BSC (column-major)
///   F. Sorted + aligned + consensus-residual: shift-align then residual + BSC
///
/// Usage: cargo run --release --bin bench_columnar [fastq_path]

use std::io::BufRead;
use std::time::Instant;

use qz_lib::compression::bsc;
use rayon::prelude::*;

// ── Base encoding ────────────────────────────────────────────────────────

const GAP: u8 = 255;

#[inline]
fn base_to_2bit(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0, // N → A (rare, handle separately if needed)
    }
}

#[inline]
fn bit2_to_base(v: u8) -> u8 {
    match v {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

// ── Minimizer sort key ───────────────────────────────────────────────────

/// Compute sort key from minimum canonical k-mer hash.
/// Uses k=21 with simple rolling hash for speed.
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
            let b = seq[i + j];
            let v = match b {
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
            let canonical = fwd.min(rev);
            min_hash = min_hash.min(canonical);
        }
    }

    min_hash
}

// ── Shift alignment via shared k-mers ────────────────────────────────────

/// Find the best shift of `query` relative to `reference` using shared k-mer positions.
/// Returns the shift offset (positive = query starts later in genomic coords).
/// Shift of 0 means the reads are already aligned.
fn find_best_shift(reference: &[u8], query: &[u8], k: usize) -> i32 {
    use rustc_hash::FxHashMap;

    if reference.len() < k || query.len() < k {
        return 0;
    }

    // Build k-mer → position map for reference
    let mut ref_kmers: FxHashMap<u64, i32> = FxHashMap::default();
    for i in 0..=reference.len() - k {
        let h = kmer_hash(&reference[i..i + k]);
        if let Some(h) = h {
            ref_kmers.entry(h).or_insert(i as i32);
        }
    }

    // For each k-mer in query, look up in reference and vote on shift
    let mut shift_votes: FxHashMap<i32, u32> = FxHashMap::default();
    for i in 0..=query.len() - k {
        let h = kmer_hash(&query[i..i + k]);
        if let Some(h) = h {
            if let Some(&ref_pos) = ref_kmers.get(&h) {
                let shift = i as i32 - ref_pos;
                *shift_votes.entry(shift).or_insert(0) += 1;
            }
        }
    }

    // Return shift with most votes
    shift_votes
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(shift, _)| shift)
        .unwrap_or(0)
}

#[inline]
fn kmer_hash(kmer: &[u8]) -> Option<u64> {
    let mut h: u64 = 0;
    for &b in kmer {
        let v = match b {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
        h = (h << 2) | v;
    }
    Some(h)
}

// ── Matrix operations ────────────────────────────────────────────────────

/// Build a column-major consensus + residual from a set of equal-length sequences.
/// Returns (consensus, residual_column_major).
///
/// Residual encoding: 0 = matches consensus, 1/2/3 = alternative base index.
/// For a consensus base C, the alternatives are the other 3 bases in a fixed order.
fn consensus_residual(sequences: &[&[u8]], read_len: usize) -> (Vec<u8>, Vec<u8>) {
    let n = sequences.len();
    let mut consensus = vec![0u8; read_len];
    // Column-major residual: column 0 for all reads, then column 1, etc.
    let mut residual = vec![0u8; n * read_len];

    for col in 0..read_len {
        // Count bases in this column
        let mut counts = [0u32; 4];
        for row in 0..n {
            if col < sequences[row].len() {
                counts[base_to_2bit(sequences[row][col]) as usize] += 1;
            }
        }

        // Consensus = most frequent base
        let consensus_base = counts
            .iter()
            .enumerate()
            .max_by_key(|&(_, &c)| c)
            .unwrap()
            .0 as u8;
        consensus[col] = consensus_base;

        // Build alternative base mapping: 1,2,3 → the other 3 bases in order
        // (This is deterministic given the consensus base, so decoder can reconstruct)
        let mut alt_idx = [0u8; 4]; // base → residual code
        alt_idx[consensus_base as usize] = 0; // match
        let mut next = 1u8;
        for b in 0..4u8 {
            if b != consensus_base {
                alt_idx[b as usize] = next;
                next += 1;
            }
        }

        // Encode residuals (column-major order)
        for row in 0..n {
            let base = if col < sequences[row].len() {
                base_to_2bit(sequences[row][col])
            } else {
                consensus_base // treat gap as consensus (shouldn't happen for equal-length)
            };
            residual[col * n + row] = alt_idx[base as usize];
        }
    }

    (consensus, residual)
}

/// Build consensus + residual for shift-aligned reads.
/// Reads are placed on a virtual coordinate system based on their shifts.
/// Returns (consensus, residual_column_major, matrix_width).
fn aligned_consensus_residual(
    sequences: &[&[u8]],
    shifts: &[i32],
) -> (Vec<u8>, Vec<u8>, usize) {
    let n = sequences.len();
    if n == 0 {
        return (Vec::new(), Vec::new(), 0);
    }

    // Find the virtual coordinate range
    let min_shift = *shifts.iter().min().unwrap();
    let max_end = shifts
        .iter()
        .zip(sequences.iter())
        .map(|(&s, seq)| s + seq.len() as i32)
        .max()
        .unwrap();
    let width = (max_end - min_shift) as usize;

    let mut consensus = vec![0u8; width];
    // Column-major residual — but only for positions covered by each read.
    // For positions not covered, we store GAP (255).
    // But this wastes space. Instead, store only covered positions per column.

    // First, count coverage per column to decide if column-major helps.
    // For simplicity, use a flat matrix with GAP markers, then strip later.
    let mut matrix = vec![GAP; n * width]; // row-major: matrix[row * width + col]

    for row in 0..n {
        let offset = (shifts[row] - min_shift) as usize;
        for j in 0..sequences[row].len() {
            matrix[row * width + offset + j] = base_to_2bit(sequences[row][j]);
        }
    }

    // Per-column consensus (ignoring GAPs)
    for col in 0..width {
        let mut counts = [0u32; 4];
        for row in 0..n {
            let v = matrix[row * width + col];
            if v != GAP {
                counts[v as usize] += 1;
            }
        }
        let total: u32 = counts.iter().sum();
        if total == 0 {
            consensus[col] = 0; // all gaps, doesn't matter
        } else {
            consensus[col] = counts
                .iter()
                .enumerate()
                .max_by_key(|&(_, &c)| c)
                .unwrap()
                .0 as u8;
        }
    }

    // Encode residuals: only for covered positions, column-major
    // Format: for each column, for each row: if covered, emit residual; else skip.
    // But decoder needs to know coverage. Store as: [residual for covered reads].
    // The decoder knows shifts and read lengths, so it knows which reads cover each column.

    // For the benchmark, we just measure the size of the residual stream.
    // Encode column-major, including GAPs (the decoder knows to skip them).
    let mut residual = Vec::with_capacity(n * width);

    for col in 0..width {
        let cons = consensus[col];
        let mut alt_idx = [0u8; 4];
        alt_idx[cons as usize] = 0;
        let mut next = 1u8;
        for b in 0..4u8 {
            if b != cons {
                alt_idx[b as usize] = next;
                next += 1;
            }
        }

        for row in 0..n {
            let v = matrix[row * width + col];
            if v == GAP {
                residual.push(GAP);
            } else {
                residual.push(alt_idx[v as usize]);
            }
        }
    }

    (consensus, residual, width)
}

/// Transpose a row-major matrix of equal-length rows to column-major.
fn transpose(sequences: &[&[u8]], read_len: usize) -> Vec<u8> {
    let n = sequences.len();
    let mut out = Vec::with_capacity(n * read_len);

    for col in 0..read_len {
        for row in 0..n {
            if col < sequences[row].len() {
                out.push(sequences[row][col]);
            } else {
                out.push(b'A'); // padding, shouldn't happen for equal-length
            }
        }
    }

    out
}

/// Delta-encode sequences: each read encoded as diff vs previous.
/// 0 = match, actual base otherwise.
fn delta_encode(sequences: &[&[u8]]) -> Vec<u8> {
    let mut out = Vec::with_capacity(sequences.len() * sequences[0].len());

    for (i, seq) in sequences.iter().enumerate() {
        if i == 0 {
            out.extend_from_slice(seq);
        } else {
            let prev = sequences[i - 1];
            let overlap = seq.len().min(prev.len());
            for j in 0..seq.len() {
                if j < overlap && seq[j] == prev[j] {
                    out.push(0);
                } else {
                    out.push(seq[j]);
                }
            }
        }
    }

    out
}

// ── Main ─────────────────────────────────────────────────────────────────

fn main() {
    let fastq_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/ERR3239334_1.1m.fastq".to_string());

    eprintln!("=== Columnar Compression Benchmark ===");
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
    let all_same_len = sequences.iter().all(|s| s.len() == read_len);
    eprintln!(
        "Read length: {} (uniform: {})",
        read_len, all_same_len
    );

    let raw_size = num_reads * read_len;
    eprintln!("Raw sequence data: {} bytes", raw_size);

    // ── Strategy A: Baseline (original order, raw BSC) ───────────────────

    eprintln!("\n--- Strategy A: Baseline (original order + BSC) ---");
    let t = Instant::now();
    let raw_concat: Vec<u8> = sequences.iter().flat_map(|s| s.iter().copied()).collect();
    let baseline_compressed = bsc::compress_parallel_adaptive(&raw_concat).unwrap();
    let baseline_size = baseline_compressed.len();
    eprintln!(
        "  Compressed: {} bytes ({:.3}x, {:.3} bits/base) in {:.2}s",
        baseline_size,
        raw_size as f64 / baseline_size as f64,
        8.0 * baseline_size as f64 / raw_size as f64,
        t.elapsed().as_secs_f64()
    );
    drop(baseline_compressed);

    // ── Sort reads by minimizer ──────────────────────────────────────────

    eprintln!("\nSorting reads by minimizer hash...");
    let t = Instant::now();
    let mut sort_keys: Vec<(u64, usize)> = sequences
        .par_iter()
        .enumerate()
        .map(|(i, seq)| (minimizer_sort_key(seq), i))
        .collect();
    sort_keys.sort_unstable();
    let sorted_indices: Vec<usize> = sort_keys.iter().map(|&(_, i)| i).collect();
    eprintln!("  Sorted in {:.2}s", t.elapsed().as_secs_f64());

    // Build sorted sequence refs
    let sorted_seqs: Vec<&[u8]> = sorted_indices.iter().map(|&i| sequences[i].as_slice()).collect();

    // ── Strategy B: Sorted + raw BSC ─────────────────────────────────────

    eprintln!("\n--- Strategy B: Sorted + BSC (row-major) ---");
    let t = Instant::now();
    let sorted_concat: Vec<u8> = sorted_seqs.iter().flat_map(|s| s.iter().copied()).collect();
    let sorted_compressed = bsc::compress_parallel_adaptive(&sorted_concat).unwrap();
    let sorted_size = sorted_compressed.len();
    eprintln!(
        "  Compressed: {} bytes ({:.3}x, {:.3} bits/base) in {:.2}s",
        sorted_size,
        raw_size as f64 / sorted_size as f64,
        8.0 * sorted_size as f64 / raw_size as f64,
        t.elapsed().as_secs_f64()
    );
    drop(sorted_compressed);

    // Permutation cost (store sorted_indices as u32, BSC compress)
    let perm_bytes: Vec<u8> = sorted_indices
        .iter()
        .flat_map(|&i| (i as u32).to_le_bytes())
        .collect();
    let perm_compressed = bsc::compress_parallel_adaptive(&perm_bytes).unwrap();
    let perm_cost = perm_compressed.len();
    eprintln!("  Permutation cost: {} bytes (BSC'd)", perm_cost);
    eprintln!(
        "  Net total: {} bytes ({:.3}x)",
        sorted_size + perm_cost,
        raw_size as f64 / (sorted_size + perm_cost) as f64
    );
    drop(perm_compressed);

    // ── Strategy C: Sorted + delta encoding ──────────────────────────────

    eprintln!("\n--- Strategy C: Sorted + delta vs previous + BSC ---");
    let t = Instant::now();
    let delta_data = delta_encode(&sorted_seqs);
    let delta_compressed = bsc::compress_parallel_adaptive(&delta_data).unwrap();
    let delta_size = delta_compressed.len();
    eprintln!(
        "  Compressed: {} bytes ({:.3}x, {:.3} bits/base) in {:.2}s",
        delta_size,
        raw_size as f64 / delta_size as f64,
        8.0 * delta_size as f64 / raw_size as f64,
        t.elapsed().as_secs_f64()
    );
    eprintln!(
        "  Net total: {} bytes ({:.3}x)",
        delta_size + perm_cost,
        raw_size as f64 / (delta_size + perm_cost) as f64
    );
    drop(delta_compressed);

    // ── Strategy D: Sorted + transpose (column-major) ────────────────────

    eprintln!("\n--- Strategy D: Sorted + transpose (column-major) + BSC ---");
    let t = Instant::now();
    let transposed = transpose(&sorted_seqs, read_len);
    let trans_compressed = bsc::compress_parallel_adaptive(&transposed).unwrap();
    let trans_size = trans_compressed.len();
    eprintln!(
        "  Compressed: {} bytes ({:.3}x, {:.3} bits/base) in {:.2}s",
        trans_size,
        raw_size as f64 / trans_size as f64,
        8.0 * trans_size as f64 / raw_size as f64,
        t.elapsed().as_secs_f64()
    );
    eprintln!(
        "  Net total: {} bytes ({:.3}x)",
        trans_size + perm_cost,
        raw_size as f64 / (trans_size + perm_cost) as f64
    );
    drop(trans_compressed);

    // ── Strategy E: Sorted + consensus-residual (column-major) ───────────

    eprintln!("\n--- Strategy E: Sorted + consensus-residual (column-major) + BSC ---");
    let t = Instant::now();
    let (consensus, residual) = consensus_residual(&sorted_seqs, read_len);
    let consensus_compressed = bsc::compress_adaptive(&consensus).unwrap();
    let residual_compressed = bsc::compress_parallel_adaptive(&residual).unwrap();
    let cons_res_size = consensus_compressed.len() + residual_compressed.len();
    eprintln!(
        "  Consensus: {} bytes compressed",
        consensus_compressed.len()
    );
    eprintln!(
        "  Residual: {} bytes compressed ({:.3} bits/base)",
        residual_compressed.len(),
        8.0 * residual_compressed.len() as f64 / raw_size as f64
    );
    eprintln!(
        "  Total: {} bytes ({:.3}x, {:.3} bits/base) in {:.2}s",
        cons_res_size,
        raw_size as f64 / cons_res_size as f64,
        8.0 * cons_res_size as f64 / raw_size as f64,
        t.elapsed().as_secs_f64()
    );
    eprintln!(
        "  Net total: {} bytes ({:.3}x)",
        cons_res_size + perm_cost,
        raw_size as f64 / (cons_res_size + perm_cost) as f64
    );

    // Verify roundtrip
    let dec_consensus = bsc::decompress(&consensus_compressed).unwrap();
    let dec_residual = bsc::decompress_parallel(&residual_compressed).unwrap();
    assert_eq!(consensus, dec_consensus);
    assert_eq!(residual, dec_residual);

    // Reconstruct sequences from consensus + residual
    for row in 0..sorted_seqs.len().min(100) {
        for col in 0..read_len {
            let cons = consensus[col];
            let res_code = residual[col * num_reads + row];
            let reconstructed = if res_code == 0 {
                bit2_to_base(cons)
            } else {
                // alt_idx maps: 1,2,3 → bases in order excluding consensus
                let mut alt_bases = Vec::new();
                for b in 0..4u8 {
                    if b != cons {
                        alt_bases.push(b);
                    }
                }
                bit2_to_base(alt_bases[(res_code - 1) as usize])
            };
            assert_eq!(
                reconstructed, sorted_seqs[row][col],
                "Mismatch at row={} col={}",
                row, col
            );
        }
    }
    eprintln!("  Roundtrip verified (first 100 reads)");
    drop(consensus_compressed);
    drop(residual_compressed);

    // ── Strategy F: Sorted + shift-aligned + consensus-residual ──────────

    eprintln!("\n--- Strategy F: Sorted + shift-aligned + consensus-residual + BSC ---");
    let t = Instant::now();

    // Process in chunks to manage memory and make alignment local
    let chunk_size = 10_000usize; // align within chunks of 10K reads
    let align_k = 15; // k-mer size for shift detection
    let mut total_residual_bytes = 0usize;
    let mut total_consensus_bytes = 0usize;
    let mut total_shift_bytes = 0usize;
    let mut total_width = 0usize;
    let mut total_covered = 0usize;
    let mut chunk_compressed_parts: Vec<Vec<u8>> = Vec::new();

    let num_chunks = (num_reads + chunk_size - 1) / chunk_size;

    for chunk_idx in 0..num_chunks {
        let start = chunk_idx * chunk_size;
        let end = (start + chunk_size).min(num_reads);
        let chunk_seqs = &sorted_seqs[start..end];
        let cn = chunk_seqs.len();

        // Compute shifts relative to first read in chunk
        let mut shifts = vec![0i32; cn];
        for i in 1..cn {
            // Align to previous read and accumulate shifts
            let delta = find_best_shift(chunk_seqs[i - 1], chunk_seqs[i], align_k);
            shifts[i] = shifts[i - 1] + delta;
        }

        // Build aligned consensus + residual
        let (cons, res, width) = aligned_consensus_residual(chunk_seqs, &shifts);

        // Count covered (non-GAP) cells
        let covered = res.iter().filter(|&&v| v != GAP).count();
        total_covered += covered;
        total_width += width;

        // Compress: strip GAPs from residual, compress separately
        // Option 1: Just compress the full residual including GAPs (simpler)
        let cons_c = bsc::compress_adaptive(&cons).unwrap();
        let res_c = bsc::compress_parallel_adaptive(&res).unwrap();

        // Shift stream: delta-encode shifts, BSC compress
        let shift_deltas: Vec<u8> = shifts
            .windows(2)
            .flat_map(|w| ((w[1] - w[0]) as i16).to_le_bytes())
            .collect();
        let shift_c = if shift_deltas.is_empty() {
            Vec::new()
        } else {
            bsc::compress_adaptive(&shift_deltas).unwrap()
        };

        total_consensus_bytes += cons_c.len();
        total_residual_bytes += res_c.len();
        total_shift_bytes += shift_c.len();

        chunk_compressed_parts.push(res_c);
    }

    let aligned_total = total_consensus_bytes + total_residual_bytes + total_shift_bytes;
    eprintln!("  Chunks: {} x {} reads", num_chunks, chunk_size);
    eprintln!("  Avg matrix width: {}", total_width / num_chunks);
    eprintln!(
        "  Coverage density: {:.1}%",
        100.0 * total_covered as f64 / (total_width * chunk_size / num_chunks) as f64
    );
    eprintln!("  Consensus: {} bytes compressed", total_consensus_bytes);
    eprintln!("  Residuals: {} bytes compressed", total_residual_bytes);
    eprintln!("  Shifts: {} bytes compressed", total_shift_bytes);
    eprintln!(
        "  Total: {} bytes ({:.3}x, {:.3} bits/base) in {:.2}s",
        aligned_total,
        raw_size as f64 / aligned_total as f64,
        8.0 * aligned_total as f64 / raw_size as f64,
        t.elapsed().as_secs_f64()
    );
    eprintln!(
        "  Net total: {} bytes ({:.3}x)",
        aligned_total + perm_cost,
        raw_size as f64 / (aligned_total + perm_cost) as f64
    );

    // ── Strategy G: Sorted + consensus-residual, GAP-free (covered only) ─

    eprintln!("\n--- Strategy G: Sorted + aligned, GAP-stripped residual + BSC ---");
    let t = Instant::now();

    let mut total_gapfree_residual = 0usize;
    let mut total_gapfree_consensus = 0usize;
    let mut total_gapfree_shifts = 0usize;

    for chunk_idx in 0..num_chunks {
        let start = chunk_idx * chunk_size;
        let end = (start + chunk_size).min(num_reads);
        let chunk_seqs = &sorted_seqs[start..end];
        let cn = chunk_seqs.len();

        let mut shifts = vec![0i32; cn];
        for i in 1..cn {
            let delta = find_best_shift(chunk_seqs[i - 1], chunk_seqs[i], align_k);
            shifts[i] = shifts[i - 1] + delta;
        }

        let (cons, res, _width) = aligned_consensus_residual(chunk_seqs, &shifts);

        // Strip GAPs: only keep covered residuals
        let gapfree: Vec<u8> = res.iter().copied().filter(|&v| v != GAP).collect();

        let cons_c = bsc::compress_adaptive(&cons).unwrap();
        let res_c = bsc::compress_parallel_adaptive(&gapfree).unwrap();

        let shift_deltas: Vec<u8> = shifts
            .windows(2)
            .flat_map(|w| ((w[1] - w[0]) as i16).to_le_bytes())
            .collect();
        let shift_c = if shift_deltas.is_empty() {
            Vec::new()
        } else {
            bsc::compress_adaptive(&shift_deltas).unwrap()
        };

        total_gapfree_consensus += cons_c.len();
        total_gapfree_residual += res_c.len();
        total_gapfree_shifts += shift_c.len();
    }

    let gapfree_total = total_gapfree_consensus + total_gapfree_residual + total_gapfree_shifts;
    eprintln!(
        "  Total: {} bytes ({:.3}x, {:.3} bits/base) in {:.2}s",
        gapfree_total,
        raw_size as f64 / gapfree_total as f64,
        8.0 * gapfree_total as f64 / raw_size as f64,
        t.elapsed().as_secs_f64()
    );
    eprintln!(
        "  Net total: {} bytes ({:.3}x)",
        gapfree_total + perm_cost,
        raw_size as f64 / (gapfree_total + perm_cost) as f64
    );

    // ── Summary ──────────────────────────────────────────────────────────

    eprintln!("\n=== Summary ===");
    eprintln!(
        "  Raw:                         {:>12} bytes",
        raw_size
    );
    eprintln!(
        "  A. Baseline (no sort):       {:>12} bytes  {:.3}x  {:.3} bpb",
        baseline_size,
        raw_size as f64 / baseline_size as f64,
        8.0 * baseline_size as f64 / raw_size as f64
    );
    eprintln!(
        "  B. Sorted row-major:         {:>12} bytes  {:.3}x  {:.3} bpb  (net {:.3}x)",
        sorted_size,
        raw_size as f64 / sorted_size as f64,
        8.0 * sorted_size as f64 / raw_size as f64,
        raw_size as f64 / (sorted_size + perm_cost) as f64
    );
    eprintln!(
        "  C. Sorted + delta:           {:>12} bytes  {:.3}x  {:.3} bpb  (net {:.3}x)",
        delta_size,
        raw_size as f64 / delta_size as f64,
        8.0 * delta_size as f64 / raw_size as f64,
        raw_size as f64 / (delta_size + perm_cost) as f64
    );
    eprintln!(
        "  D. Sorted + transpose:       {:>12} bytes  {:.3}x  {:.3} bpb  (net {:.3}x)",
        trans_size,
        raw_size as f64 / trans_size as f64,
        8.0 * trans_size as f64 / raw_size as f64,
        raw_size as f64 / (trans_size + perm_cost) as f64
    );
    eprintln!(
        "  E. Sorted + cons-residual:   {:>12} bytes  {:.3}x  {:.3} bpb  (net {:.3}x)",
        cons_res_size,
        raw_size as f64 / cons_res_size as f64,
        8.0 * cons_res_size as f64 / raw_size as f64,
        raw_size as f64 / (cons_res_size + perm_cost) as f64
    );
    eprintln!(
        "  F. Aligned + cons-residual:  {:>12} bytes  {:.3}x  {:.3} bpb  (net {:.3}x)",
        aligned_total,
        raw_size as f64 / aligned_total as f64,
        8.0 * aligned_total as f64 / raw_size as f64,
        raw_size as f64 / (aligned_total + perm_cost) as f64
    );
    eprintln!(
        "  G. Aligned + GAP-stripped:   {:>12} bytes  {:.3}x  {:.3} bpb  (net {:.3}x)",
        gapfree_total,
        raw_size as f64 / gapfree_total as f64,
        8.0 * gapfree_total as f64 / raw_size as f64,
        raw_size as f64 / (gapfree_total + perm_cost) as f64
    );
    eprintln!("  Permutation overhead:        {:>12} bytes", perm_cost);
}
