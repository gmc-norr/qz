/// Prototype: k-mer frequency prior for sequence compression
///
/// Idea: build a table of P(next_base | preceding k-mer) from all reads,
/// then use it to either:
///   - Compute theoretical conditional entropy (lower bound)
///   - Rank-recode each base (0=most likely, 3=least likely) and BSC compress
///   - Arithmetic-code with the k-mer prior directly
///
/// Strategies tested:
///   A. Baseline: raw sequences + BSC
///   B. Theoretical entropy for k-mer models (k=1..11)
///   C. Rank-recoded + BSC (recode each base by predicted rank)
///   D. Sorted + rank-recoded + BSC
///   E. Arithmetic coding with k-mer prior (on subset for validation)
///
/// Usage: cargo run --release --bin bench_kmer_prior [fastq_path]

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
        _ => 0, // N → A
    }
}

#[inline]
fn bit2_to_base(v: u8) -> u8 {
    match v {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

// ── k-mer frequency table ────────────────────────────────────────────────

/// A table mapping each k-mer to counts of the next base [A, C, G, T].
/// k-mer is encoded as 2*k bits packed into a usize.
/// For k ≤ 13 this fits in a u32 index (4^13 = 67M entries).
struct KmerTable {
    k: usize,
    /// counts[kmer_hash * 4 + base] = count of base following this k-mer
    counts: Vec<u32>,
    mask: usize,
}

impl KmerTable {
    fn new(k: usize) -> Self {
        let num_entries = 1usize << (2 * k);
        Self {
            k,
            counts: vec![0u32; num_entries * 4],
            mask: num_entries - 1,
        }
    }

    /// Memory usage in bytes
    fn memory(&self) -> usize {
        self.counts.len() * 4
    }

    /// Count all k-mer → next base transitions in a set of sequences.
    fn train(&mut self, sequences: &[Vec<u8>]) {
        for seq in sequences {
            if seq.len() <= self.k {
                continue;
            }
            // Build initial k-mer hash
            let mut hash = 0usize;
            let mut valid = true;
            for i in 0..self.k {
                let v = base_to_2bit(seq[i]) as usize;
                if seq[i] == b'N' || seq[i] == b'n' {
                    valid = false;
                }
                hash = (hash << 2) | v;
            }
            hash &= self.mask;

            // Slide through sequence
            for i in self.k..seq.len() {
                let next = base_to_2bit(seq[i]) as usize;
                let is_n = seq[i] == b'N' || seq[i] == b'n';

                if valid && !is_n {
                    self.counts[hash * 4 + next] += 1;
                }

                // Update hash: shift out oldest base, shift in new base
                hash = ((hash << 2) | next) & self.mask;
                // Check if the new k-mer window contains N
                // (approximate: we track validity by checking current base only,
                //  which is imperfect but fast. True N-tracking would need more state.)
                valid = !is_n;
            }
        }
    }

    /// Get the probability ranking for a base given its k-mer context.
    /// Returns rank 0-3 where 0 = most probable base.
    #[inline]
    fn rank_of(&self, kmer_hash: usize, base: u8) -> u8 {
        let idx = kmer_hash * 4;
        let c = [
            self.counts[idx],
            self.counts[idx + 1],
            self.counts[idx + 2],
            self.counts[idx + 3],
        ];

        // Count how many bases have strictly higher count
        let my_count = c[base as usize];
        let mut rank = 0u8;
        for b in 0..4u8 {
            if b != base && c[b as usize] > my_count {
                rank += 1;
            } else if b != base && c[b as usize] == my_count && b < base {
                // Tie-break: lower base index gets lower rank
                rank += 1;
            }
        }
        rank
    }

    /// Get the base for a given rank in a k-mer context (inverse of rank_of).
    #[inline]
    fn base_at_rank(&self, kmer_hash: usize, rank: u8) -> u8 {
        let idx = kmer_hash * 4;
        let c = [
            self.counts[idx],
            self.counts[idx + 1],
            self.counts[idx + 2],
            self.counts[idx + 3],
        ];

        // Sort bases by count descending, tie-break by base index ascending
        let mut order = [0u8, 1, 2, 3];
        order.sort_by(|&a, &b| {
            c[b as usize]
                .cmp(&c[a as usize])
                .then(a.cmp(&b))
        });

        order[rank as usize]
    }

    /// Compute conditional entropy H(next_base | k-mer) in bits/base.
    fn conditional_entropy(&self) -> f64 {
        let num_entries = 1usize << (2 * self.k);
        let mut total_bits = 0.0f64;
        let mut total_count = 0u64;

        for kmer in 0..num_entries {
            let idx = kmer * 4;
            let c = [
                self.counts[idx] as u64,
                self.counts[idx + 1] as u64,
                self.counts[idx + 2] as u64,
                self.counts[idx + 3] as u64,
            ];
            let sum: u64 = c.iter().sum();
            if sum == 0 {
                continue;
            }

            let sum_f = sum as f64;
            for &count in &c {
                if count > 0 {
                    let p = count as f64 / sum_f;
                    total_bits += -(count as f64) * p.log2();
                }
            }
            total_count += sum;
        }

        if total_count == 0 {
            2.0
        } else {
            total_bits / total_count as f64
        }
    }

    /// Compute prediction accuracy (fraction of bases where rank=0).
    fn prediction_accuracy(&self) -> f64 {
        let num_entries = 1usize << (2 * self.k);
        let mut correct = 0u64;
        let mut total = 0u64;

        for kmer in 0..num_entries {
            let idx = kmer * 4;
            let c = [
                self.counts[idx] as u64,
                self.counts[idx + 1] as u64,
                self.counts[idx + 2] as u64,
                self.counts[idx + 3] as u64,
            ];
            let sum: u64 = c.iter().sum();
            if sum == 0 {
                continue;
            }
            correct += *c.iter().max().unwrap();
            total += sum;
        }

        if total == 0 {
            0.25
        } else {
            correct as f64 / total as f64
        }
    }

    /// Serialize the table to bytes (for storage cost estimation).
    fn serialize(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.counts.len() * 4 + 8);
        out.extend_from_slice(&(self.k as u32).to_le_bytes());
        out.extend_from_slice(&(self.counts.len() as u32).to_le_bytes());
        for &c in &self.counts {
            out.extend_from_slice(&c.to_le_bytes());
        }
        out
    }
}

// ── Rank recoding ────────────────────────────────────────────────────────

/// Recode a sequence using the k-mer table: each base → its rank (0-3).
fn rank_recode_sequence(seq: &[u8], table: &KmerTable) -> Vec<u8> {
    let k = table.k;
    let mut out = Vec::with_capacity(seq.len());

    if seq.len() <= k {
        // Can't use context, store raw 2-bit
        for &b in seq {
            out.push(base_to_2bit(b));
        }
        return out;
    }

    // First k bases: no context, store raw 2-bit
    let mut hash = 0usize;
    for i in 0..k {
        let v = base_to_2bit(seq[i]);
        out.push(v);
        hash = ((hash << 2) | v as usize) & table.mask;
    }

    // Remaining bases: encode as rank
    for i in k..seq.len() {
        let base = base_to_2bit(seq[i]);
        let rank = table.rank_of(hash, base);
        out.push(rank);
        hash = ((hash << 2) | base as usize) & table.mask;
    }

    out
}

/// Decode a rank-recoded sequence back to bases.
fn rank_decode_sequence(encoded: &[u8], table: &KmerTable) -> Vec<u8> {
    let k = table.k;
    let mut out = Vec::with_capacity(encoded.len());

    if encoded.len() <= k {
        for &v in encoded {
            out.push(bit2_to_base(v));
        }
        return out;
    }

    // First k bases: raw 2-bit
    let mut hash = 0usize;
    for i in 0..k {
        let base = encoded[i];
        out.push(bit2_to_base(base));
        hash = ((hash << 2) | base as usize) & table.mask;
    }

    // Remaining: rank → base using table
    for i in k..encoded.len() {
        let rank = encoded[i];
        let base = table.base_at_rank(hash, rank);
        out.push(bit2_to_base(base));
        hash = ((hash << 2) | base as usize) & table.mask;
    }

    out
}

// ── Minimizer sort (from bench_columnar) ─────────────────────────────────

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

// ── Arithmetic coding (small-scale validation) ───────────────────────────

/// Compute the exact bit cost of encoding sequences with the k-mer prior,
/// without actually doing arithmetic coding. Uses log2(1/p) per symbol.
fn compute_bit_cost(sequences: &[Vec<u8>], table: &KmerTable) -> f64 {
    let k = table.k;
    let mut total_bits = 0.0f64;
    let mut total_bases = 0u64;

    for seq in sequences {
        if seq.len() <= k {
            total_bits += seq.len() as f64 * 2.0; // 2 bits/base for context-free
            total_bases += seq.len() as u64;
            continue;
        }

        // First k bases: 2 bits each (no context)
        let mut hash = 0usize;
        for i in 0..k {
            total_bits += 2.0;
            let v = base_to_2bit(seq[i]) as usize;
            hash = ((hash << 2) | v) & table.mask;
        }
        total_bases += k as u64;

        // Remaining bases: -log2(p) bits
        for i in k..seq.len() {
            let base = base_to_2bit(seq[i]) as usize;
            let idx = hash * 4;
            let c = [
                table.counts[idx] as f64,
                table.counts[idx + 1] as f64,
                table.counts[idx + 2] as f64,
                table.counts[idx + 3] as f64,
            ];
            let total: f64 = c.iter().sum();

            if total == 0.0 || c[base] == 0.0 {
                total_bits += 2.0; // uniform fallback
            } else {
                let p = c[base] / total;
                total_bits += -p.log2();
            }

            hash = ((hash << 2) | base) & table.mask;
            total_bases += 1;
        }
    }

    total_bits / total_bases as f64
}

// ── Main ─────────────────────────────────────────────────────────────────

fn main() {
    let fastq_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/ERR3239334_1.1m.fastq".to_string());

    eprintln!("=== k-mer Frequency Prior Benchmark ===");
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
    eprintln!("Read length: {}, raw size: {} bytes", read_len, raw_size);

    // ── Strategy A: Baseline BSC ─────────────────────────────────────────

    eprintln!("\n--- Strategy A: Baseline (original order + BSC) ---");
    let t = Instant::now();
    let raw_concat: Vec<u8> = sequences.iter().flat_map(|s| s.iter().copied()).collect();
    let baseline_compressed = bsc::compress_parallel_adaptive(&raw_concat).unwrap();
    let baseline_size = baseline_compressed.len();
    eprintln!(
        "  Compressed: {} bytes ({:.3}x, {:.4} bpb) in {:.2}s",
        baseline_size,
        raw_size as f64 / baseline_size as f64,
        8.0 * baseline_size as f64 / raw_size as f64,
        t.elapsed().as_secs_f64()
    );
    drop(baseline_compressed);
    drop(raw_concat);

    // ── Strategy B: Theoretical entropy for various k ────────────────────

    eprintln!("\n--- Strategy B: Theoretical conditional entropy H(base|k-mer) ---");
    let max_k = 11;
    let mut tables: Vec<KmerTable> = Vec::new();

    for k in 1..=max_k {
        let mem = (1usize << (2 * k)) * 4 * 4; // bytes
        if mem > 2_000_000_000 {
            eprintln!("  k={}: table too large ({} MB), skipping", k, mem / (1024 * 1024));
            break;
        }

        let t = Instant::now();
        let mut table = KmerTable::new(k);
        table.train(&sequences);
        let entropy = table.conditional_entropy();
        let accuracy = table.prediction_accuracy();
        let bit_cost = compute_bit_cost(&sequences, &table);
        let table_raw = table.serialize();
        let table_compressed = bsc::compress_adaptive(&table_raw).unwrap();

        eprintln!(
            "  k={:2}: H={:.4} bpb  pred_acc={:.1}%  actual_cost={:.4} bpb  table={} B ({} B compressed)  train={:.2}s",
            k,
            entropy,
            accuracy * 100.0,
            bit_cost,
            table.memory(),
            table_compressed.len(),
            t.elapsed().as_secs_f64()
        );

        tables.push(table);
    }

    // ── Strategy C: Rank-recoded + BSC (various k) ───────────────────────

    eprintln!("\n--- Strategy C: Rank-recoded + BSC (original order) ---");

    let mut best_k = 0;
    let mut best_rank_size = usize::MAX;

    for table in &tables {
        let k = table.k;
        let t = Instant::now();

        // Rank-recode all sequences
        let recoded: Vec<Vec<u8>> = sequences
            .par_iter()
            .map(|seq| rank_recode_sequence(seq, table))
            .collect();

        // Concatenate and BSC compress
        let recoded_concat: Vec<u8> = recoded.iter().flat_map(|s| s.iter().copied()).collect();
        let compressed = bsc::compress_parallel_adaptive(&recoded_concat).unwrap();

        // Table storage cost
        let table_raw = table.serialize();
        let table_compressed = bsc::compress_adaptive(&table_raw).unwrap();
        let total = compressed.len() + table_compressed.len();

        eprintln!(
            "  k={:2}: recoded={} B  table={} B  total={} B ({:.3}x, {:.4} bpb) in {:.2}s",
            k,
            compressed.len(),
            table_compressed.len(),
            total,
            raw_size as f64 / total as f64,
            8.0 * total as f64 / raw_size as f64,
            t.elapsed().as_secs_f64()
        );

        if total < best_rank_size {
            best_rank_size = total;
            best_k = k;
        }

        // Verify roundtrip for first 10 reads
        for i in 0..10.min(num_reads) {
            let decoded = rank_decode_sequence(&recoded[i], table);
            assert_eq!(
                decoded,
                sequences[i].iter().map(|&b| {
                    match b {
                        b'A' | b'a' => b'A',
                        b'C' | b'c' => b'C',
                        b'G' | b'g' => b'G',
                        b'T' | b't' => b'T',
                        _ => b'A',
                    }
                }).collect::<Vec<u8>>(),
                "Roundtrip failed for read {}",
                i
            );
        }
    }

    eprintln!("  Best k={} → {} bytes ({:.3}x)", best_k, best_rank_size,
        raw_size as f64 / best_rank_size as f64);

    // ── Strategy D: Sorted + rank-recoded + BSC ──────────────────────────

    eprintln!("\n--- Strategy D: Sorted + rank-recoded + BSC ---");

    // Sort reads
    let t = Instant::now();
    let mut sort_keys: Vec<(u64, usize)> = sequences
        .par_iter()
        .enumerate()
        .map(|(i, seq)| (minimizer_sort_key(seq), i))
        .collect();
    sort_keys.sort_unstable();
    let sorted_indices: Vec<usize> = sort_keys.iter().map(|&(_, i)| i).collect();
    eprintln!("  Sorted in {:.2}s", t.elapsed().as_secs_f64());

    // Permutation cost
    let perm_bytes: Vec<u8> = sorted_indices
        .iter()
        .flat_map(|&i| (i as u32).to_le_bytes())
        .collect();
    let perm_compressed = bsc::compress_parallel_adaptive(&perm_bytes).unwrap();
    let perm_cost = perm_compressed.len();
    eprintln!("  Permutation cost: {} bytes", perm_cost);
    drop(perm_compressed);

    // Sorted sequences reference
    let sorted_seqs: Vec<&[u8]> = sorted_indices.iter().map(|&i| sequences[i].as_slice()).collect();

    // Sorted baseline (no recoding)
    let sorted_concat: Vec<u8> = sorted_seqs.iter().flat_map(|s| s.iter().copied()).collect();
    let sorted_baseline = bsc::compress_parallel_adaptive(&sorted_concat).unwrap();
    eprintln!(
        "  Sorted baseline: {} B ({:.3}x, net {:.3}x)",
        sorted_baseline.len(),
        raw_size as f64 / sorted_baseline.len() as f64,
        raw_size as f64 / (sorted_baseline.len() + perm_cost) as f64
    );
    drop(sorted_baseline);

    // Train table on sorted data (same data, different order shouldn't matter for counts)
    // Use the best k from strategy C
    let _best_table = &tables[best_k - 1]; // tables[0] is k=1

    // Also try training on sorted data specifically
    let _sorted_seqs_owned: Vec<Vec<u8>> = sorted_seqs.iter().map(|s| s.to_vec()).collect();

    for table in &tables {
        let k = table.k;
        let t = Instant::now();

        let recoded: Vec<Vec<u8>> = sorted_seqs
            .par_iter()
            .map(|seq| rank_recode_sequence(seq, table))
            .collect();

        let recoded_concat: Vec<u8> = recoded.iter().flat_map(|s| s.iter().copied()).collect();
        let compressed = bsc::compress_parallel_adaptive(&recoded_concat).unwrap();

        let table_raw = table.serialize();
        let table_compressed = bsc::compress_adaptive(&table_raw).unwrap();
        let total = compressed.len() + table_compressed.len();
        let net_total = total + perm_cost;

        eprintln!(
            "  k={:2}: recoded={} B  total={} B ({:.3}x)  net={} B ({:.3}x) in {:.2}s",
            k,
            compressed.len(),
            total,
            raw_size as f64 / total as f64,
            net_total,
            raw_size as f64 / net_total as f64,
            t.elapsed().as_secs_f64()
        );
    }

    // ── Strategy E: Position-aware k-mer prior ───────────────────────────
    //
    // Illumina reads have position-dependent error profiles:
    // - Read starts: adapter/primer artifacts
    // - Read ends: higher error rate
    // Use separate k-mer tables for position bins.

    eprintln!("\n--- Strategy E: Position-binned k-mer prior + BSC ---");
    let num_pos_bins = 5; // divide read into 5 position bins
    let test_k = best_k.max(3); // use at least k=3

    {
        let t = Instant::now();

        // Train per-bin tables
        let mut bin_tables: Vec<KmerTable> = (0..num_pos_bins)
            .map(|_| KmerTable::new(test_k))
            .collect();

        for seq in &sequences {
            if seq.len() <= test_k {
                continue;
            }
            let rl = seq.len();
            let mut hash = 0usize;
            for i in 0..test_k {
                let v = base_to_2bit(seq[i]) as usize;
                hash = ((hash << 2) | v) & bin_tables[0].mask;
            }
            for i in test_k..seq.len() {
                let bin = (i * num_pos_bins) / rl;
                let bin = bin.min(num_pos_bins - 1);
                let base = base_to_2bit(seq[i]) as usize;
                bin_tables[bin].counts[hash * 4 + base] += 1;
                hash = ((hash << 2) | base) & bin_tables[0].mask;
            }
        }

        // Report per-bin entropy
        for (bin, bt) in bin_tables.iter().enumerate() {
            let entropy = bt.conditional_entropy();
            let accuracy = bt.prediction_accuracy();
            eprintln!(
                "  Bin {} (pos {}-{}): H={:.4} bpb, pred_acc={:.1}%",
                bin,
                bin * read_len / num_pos_bins,
                (bin + 1) * read_len / num_pos_bins,
                entropy,
                accuracy * 100.0
            );
        }

        // Rank-recode with position-aware tables
        let recoded: Vec<Vec<u8>> = sequences
            .par_iter()
            .map(|seq| {
                let k = test_k;
                let rl = seq.len();
                let mut out = Vec::with_capacity(rl);

                if rl <= k {
                    for &b in seq.iter() {
                        out.push(base_to_2bit(b));
                    }
                    return out;
                }

                let mut hash = 0usize;
                for i in 0..k {
                    let v = base_to_2bit(seq[i]);
                    out.push(v);
                    hash = ((hash << 2) | v as usize) & bin_tables[0].mask;
                }

                for i in k..rl {
                    let bin = (i * num_pos_bins) / rl;
                    let bin = bin.min(num_pos_bins - 1);
                    let base = base_to_2bit(seq[i]);
                    let rank = bin_tables[bin].rank_of(hash, base);
                    out.push(rank);
                    hash = ((hash << 2) | base as usize) & bin_tables[0].mask;
                }

                out
            })
            .collect();

        let recoded_concat: Vec<u8> = recoded.iter().flat_map(|s| s.iter().copied()).collect();
        let compressed = bsc::compress_parallel_adaptive(&recoded_concat).unwrap();

        // Table storage
        let mut tables_raw = Vec::new();
        for bt in &bin_tables {
            tables_raw.extend_from_slice(&bt.serialize());
        }
        let tables_compressed = bsc::compress_adaptive(&tables_raw).unwrap();

        let total = compressed.len() + tables_compressed.len();
        eprintln!(
            "  Position-aware k={}: {} B recoded + {} B tables = {} B ({:.3}x, {:.4} bpb) in {:.2}s",
            test_k,
            compressed.len(),
            tables_compressed.len(),
            total,
            raw_size as f64 / total as f64,
            8.0 * total as f64 / raw_size as f64,
            t.elapsed().as_secs_f64()
        );
    }

    // ── Strategy F: Mixed-order context (blend k=1..best_k) ──────────────

    eprintln!("\n--- Strategy F: Mixed-order rank recoding + BSC ---");
    {
        let t = Instant::now();

        // For each position, use the highest available context order.
        // Positions 0..k use lower orders. Position k+ use order k.
        // Use the table for the actual context length available.
        let recoded: Vec<Vec<u8>> = sequences
            .par_iter()
            .map(|seq| {
                let rl = seq.len();
                let mut out = Vec::with_capacity(rl);

                // Position 0: no context
                out.push(base_to_2bit(seq[0]));

                // Position 1..max_k: use available context order
                for pos in 1..rl {
                    let available_k = pos.min(tables.len()); // max order we can use
                    let table = &tables[available_k - 1]; // 0-indexed
                    let k = table.k;

                    // Build hash from preceding k bases
                    let mut hash = 0usize;
                    for j in (pos - k)..pos {
                        let v = base_to_2bit(seq[j]) as usize;
                        hash = ((hash << 2) | v) & table.mask;
                    }

                    let base = base_to_2bit(seq[pos]);
                    let rank = table.rank_of(hash, base);
                    out.push(rank);
                }

                out
            })
            .collect();

        let recoded_concat: Vec<u8> = recoded.iter().flat_map(|s| s.iter().copied()).collect();
        let compressed = bsc::compress_parallel_adaptive(&recoded_concat).unwrap();

        // Table storage: need all tables up to best_k
        let mut all_tables_raw = Vec::new();
        for table in &tables {
            all_tables_raw.extend_from_slice(&table.serialize());
        }
        let all_tables_compressed = bsc::compress_adaptive(&all_tables_raw).unwrap();

        let total = compressed.len() + all_tables_compressed.len();
        eprintln!(
            "  Mixed k=1..{}: {} B recoded + {} B tables = {} B ({:.3}x, {:.4} bpb) in {:.2}s",
            tables.len(),
            compressed.len(),
            all_tables_compressed.len(),
            total,
            raw_size as f64 / total as f64,
            8.0 * total as f64 / raw_size as f64,
            t.elapsed().as_secs_f64()
        );
    }

    // ── Summary ──────────────────────────────────────────────────────────

    eprintln!("\n=== Summary ===");
    eprintln!(
        "  Raw:                         {:>12} bytes  {:.4} bpb",
        raw_size,
        8.0
    );
    eprintln!(
        "  A. Baseline BSC:             {:>12} bytes  {:.4} bpb  ({:.3}x)",
        baseline_size,
        8.0 * baseline_size as f64 / raw_size as f64,
        raw_size as f64 / baseline_size as f64
    );
    eprintln!(
        "  C. Best rank-recode+BSC:     {:>12} bytes  {:.4} bpb  ({:.3}x)  [k={}]",
        best_rank_size,
        8.0 * best_rank_size as f64 / raw_size as f64,
        raw_size as f64 / best_rank_size as f64,
        best_k
    );
    eprintln!("  B. Theoretical entropy: see per-k results above");
    eprintln!("  D-F: see detailed results above");
}
