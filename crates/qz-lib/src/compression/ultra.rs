/// Read reordering for better BSC compression, with two encoding modes:
///
/// 1. **Reorder-local** (encoding_type=9): Group similar reads together,
///    then BSC-compress raw sequences.
///
/// 2. **HARC delta** (encoding_type=8): Group similar reads, then encode
///    each read as a delta against the previous read (0=match, 1-5=actual base).
///
/// Both modes store the permutation so decompression restores original read order.
///
/// Grouping uses an inverted index of minimizer hashes (k=21, window=10).
/// Each read has ~13 minimizers. Reads sharing ANY minimizer are grouped
/// at the first occurrence. Singletons stay in input order, preserving
/// BSC's BWT-friendly natural locality.
///
/// Chunked processing: 5M reads per chunk, temp files, pipelined I/O.

use anyhow::Result;
use tracing::info;
use std::time::Instant;

#[allow(unused_imports)]
use super::*;

// ── Parameters ────────────────────────────────────────────────────────────

const CHUNK_SIZE: usize = 5_000_000;

/// Dictionary window size for center hash (HARC/SPRING style).
const DICT_WINDOW: usize = 32;

// ── Ultra compression levels ──────────────────────────────────────────────

#[derive(Debug, Clone, Copy)]
pub(super) struct UltraLevel {
    pub level: u8,
    pub chunk_size: usize,
    pub max_parallel_chunks: usize,
    pub quality_sub_block: usize,
}

const ULTRA_LEVELS: [UltraLevel; 5] = [
    UltraLevel { level: 1, chunk_size: 1_000_000, max_parallel_chunks: 4, quality_sub_block: 250_000 },
    UltraLevel { level: 2, chunk_size: 2_000_000, max_parallel_chunks: 3, quality_sub_block: 500_000 },
    UltraLevel { level: 3, chunk_size: 5_000_000, max_parallel_chunks: 2, quality_sub_block: 500_000 },
    UltraLevel { level: 4, chunk_size: 8_000_000, max_parallel_chunks: 2, quality_sub_block: 1_000_000 },
    UltraLevel { level: 5, chunk_size: 10_000_000, max_parallel_chunks: 1, quality_sub_block: 1_000_000 },
];

fn available_memory_bytes() -> Option<u64> {
    let contents = std::fs::read_to_string("/proc/meminfo").ok()?;
    for line in contents.lines() {
        if line.starts_with("MemAvailable:") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                return parts[1].parse::<u64>().ok().map(|kb| kb * 1024);
            }
        }
    }
    None
}

fn estimate_peak_memory(level: &UltraLevel) -> u64 {
    let per_read_bytes: u64 = 1500;
    let io_overhead: u64 = 2 * 1024 * 1024 * 1024;
    level.chunk_size as u64 * per_read_bytes * level.max_parallel_chunks as u64 + io_overhead
}

fn auto_select_level() -> u8 {
    let available = match available_memory_bytes() {
        Some(bytes) => bytes,
        None => {
            info!("Could not detect available memory, defaulting to ultra level 2");
            return 2;
        }
    };

    let budget = available * 80 / 100;
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8);

    let max_level_for_cores: u8 = if cores < 4 { 1 }
        else if cores < 8 { 2 }
        else if cores < 16 { 3 }
        else { 5 };

    let mut selected = 1u8;
    for level in &ULTRA_LEVELS {
        if estimate_peak_memory(level) <= budget && level.level <= max_level_for_cores {
            selected = level.level;
        }
    }

    info!(
        "Auto-selected ultra level {} (available: {:.1} GB, budget: {:.1} GB, cores: {})",
        selected,
        available as f64 / 1e9,
        budget as f64 / 1e9,
        cores,
    );

    selected
}

pub(super) fn resolve_ultra_level(requested: u8) -> UltraLevel {
    let level = if requested == 0 {
        auto_select_level()
    } else {
        requested.clamp(1, 5)
    };
    ULTRA_LEVELS[(level - 1) as usize]
}

// ── Center hash ─────────────────────────────────────────────────────────

/// Lookup table for base → 2-bit encoding (branchless).
static BASE_2BIT_LUT: [u64; 256] = {
    let mut t = [0u64; 256];
    t[b'C' as usize] = 1; t[b'c' as usize] = 1;
    t[b'G' as usize] = 2; t[b'g' as usize] = 2;
    t[b'T' as usize] = 3; t[b't' as usize] = 3;
    t
};

/// Compute a 2-bit hash of a k-mer at a fixed position in the sequence.
#[inline]
fn hash_kmer_at(seq: &[u8], start: usize, k: usize) -> u64 {
    let mut h: u64 = 0;
    for i in start..start + k {
        h = (h << 2) | BASE_2BIT_LUT[seq[i] as usize];
    }
    h
}

/// Compute two dictionary-window hashes from the center of a read.
#[inline]
fn center_hashes(seq: &[u8]) -> (u64, u64) {
    let len = seq.len();
    let center = len / 2;
    let w = if len > 100 { DICT_WINDOW } else { len * DICT_WINDOW / 100 };
    let w = w.min(center);
    if len < 2 * w {
        return (0, 0);
    }
    (hash_kmer_at(seq, center - w, w), hash_kmer_at(seq, center, w))
}

// ── Reordering (group matching reads, preserve input order) ──────────────

struct ReorderResult {
    /// Reorder position → original read index.
    order: Vec<u32>,
}

/// Reorder reads by grouping those sharing center-hash windows.
/// Reads in groups of 2+ are pulled together at the first occurrence.
/// Singletons stay in input order, preserving BSC's BWT-friendly locality.
fn reorder_chunk(records: &[crate::io::FastqRecord]) -> ReorderResult {
    use rayon::prelude::*;
    use rustc_hash::FxHashMap;

    let n = records.len();
    if n == 0 {
        return ReorderResult { order: Vec::new() };
    }

    // Compute center hashes in parallel
    let hashes: Vec<(u64, u64)> = records.par_iter().map(|r| {
        center_hashes(&r.sequence)
    }).collect();

    // Build hash → read indices map (h1 only)
    let mut groups: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
    for (i, &(h1, _)) in hashes.iter().enumerate() {
        groups.entry(h1).or_default().push(i as u32);
    }

    // Count grouped reads for logging
    let grouped_count: usize = groups.values().filter(|g| g.len() > 1).map(|g| g.len()).sum();
    let group_count = groups.values().filter(|g| g.len() > 1).count();
    info!("  Grouped {} reads into {} groups ({:.1}%)",
        grouped_count, group_count, 100.0 * grouped_count as f64 / n as f64);

    // Build output order: group matching reads at first occurrence
    let mut used = vec![false; n];
    let mut order = Vec::with_capacity(n);

    for i in 0..n {
        if used[i] { continue; }
        let (h1, _) = hashes[i];
        if let Some(group) = groups.get(&h1) {
            if group.len() > 1 {
                for &idx in group {
                    if !used[idx as usize] {
                        used[idx as usize] = true;
                        order.push(idx);
                    }
                }
            } else {
                used[i] = true;
                order.push(i as u32);
            }
        } else {
            used[i] = true;
            order.push(i as u32);
        }
    }

    ReorderResult { order }
}

// ── Minimizer Union-Find reordering ─────────────────────────────────────

/// Union-Find (disjoint set) with path compression and union by rank.
struct UnionFind {
    parent: Vec<u32>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            self.parent[x as usize] = self.parent[self.parent[x as usize] as usize];
            x = self.parent[x as usize];
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb { return; }
        if self.rank[ra as usize] < self.rank[rb as usize] {
            self.parent[ra as usize] = rb;
        } else if self.rank[ra as usize] > self.rank[rb as usize] {
            self.parent[rb as usize] = ra;
        } else {
            self.parent[rb as usize] = ra;
            self.rank[ra as usize] += 1;
        }
    }
}

/// Minimizer parameters.
const UF_MINIMIZER_K: usize = 21;
const UF_MINIMIZER_W: u16 = 11; // must be odd for minimizer-iter
/// Skip minimizers that appear in more reads than this (repetitive k-mers).
const UF_MAX_BUCKET: usize = 50;

/// Compute canonical minimizer (hash, read_idx) pairs for all reads in parallel.
fn compute_minimizer_pairs(records: &[crate::io::FastqRecord]) -> Vec<(u64, u32)> {
    use rayon::prelude::*;
    use minimizer_iter::MinimizerBuilder;

    let min_seq_len = UF_MINIMIZER_K + UF_MINIMIZER_W as usize;
    let pairs_per_read: Vec<Vec<(u64, u32)>> = records.par_iter().enumerate().map(|(i, r)| {
        let seq = &r.sequence;
        if seq.len() < min_seq_len {
            return Vec::new();
        }
        let positions: Vec<usize> = MinimizerBuilder::<u64>::new()
            .minimizer_size(UF_MINIMIZER_K)
            .width(UF_MINIMIZER_W)
            .canonical()
            .iter_pos(seq)
            .map(|(pos, _is_rc)| pos)
            .collect();

        let mut pairs = Vec::with_capacity(positions.len());
        let mut prev_hash = u64::MAX;
        for pos in positions {
            if pos + UF_MINIMIZER_K > seq.len() { continue; }
            let kmer = &seq[pos..pos + UF_MINIMIZER_K];
            if let Some(fwd) = dna_utils::kmer_to_hash(kmer) {
                let rc = dna_utils::reverse_complement_hash(fwd, UF_MINIMIZER_K);
                let canon = fwd.min(rc);
                if canon != prev_hash {
                    pairs.push((canon, i as u32));
                    prev_hash = canon;
                }
            }
        }
        pairs
    }).collect();

    let total: usize = pairs_per_read.iter().map(|p| p.len()).sum();
    let mut all_pairs = Vec::with_capacity(total);
    for pairs in pairs_per_read {
        all_pairs.extend(pairs);
    }
    all_pairs
}

/// Emit groups in input order: when hitting an unplaced read, pull its entire group.
fn emit_groups_in_order(n: usize, roots: &[u32]) -> (ReorderResult, usize, usize) {
    use rustc_hash::FxHashMap;

    let mut components: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for (i, &root) in roots.iter().enumerate() {
        components.entry(root).or_default().push(i as u32);
    }

    let grouped_count: usize = components.values().filter(|g| g.len() > 1).map(|g| g.len()).sum();
    let group_count = components.values().filter(|g| g.len() > 1).count();

    let mut used = vec![false; n];
    let mut order = Vec::with_capacity(n);

    for i in 0..n {
        if used[i] { continue; }
        let root = roots[i];
        if let Some(members) = components.get(&root) {
            if members.len() > 1 {
                for &idx in members {
                    if !used[idx as usize] {
                        used[idx as usize] = true;
                        order.push(idx);
                    }
                }
            } else {
                used[i] = true;
                order.push(i as u32);
            }
        } else {
            used[i] = true;
            order.push(i as u32);
        }
    }

    (ReorderResult { order }, grouped_count, group_count)
}

/// Strategy A: No reorder — identity permutation. Tests whether the
/// compression gain is from block size alone.
fn reorder_identity(records: &[crate::io::FastqRecord]) -> ReorderResult {
    let n = records.len();
    info!("  Identity order (no reorder): {} reads", n);
    ReorderResult { order: (0..n as u32).collect() }
}

/// Strategy B: Selective UF — only union reads sharing ≥ 2 minimizers.
///
/// 1. Compute (hash, read_idx) pairs, sort by hash
/// 2. For each hash run, generate (read_a, read_b) edges (a < b)
/// 3. Sort edges, count duplicates → pairs with count ≥ 2 share 2+ minimizers
/// 4. Union only those pairs
fn reorder_selective_uf(records: &[crate::io::FastqRecord]) -> ReorderResult {
    use rayon::prelude::*;

    let n = records.len();
    if n == 0 {
        return ReorderResult { order: Vec::new() };
    }

    let t = Instant::now();

    // Step 1: Compute and sort minimizer pairs
    let mut all_pairs = compute_minimizer_pairs(records);
    let t1 = t.elapsed().as_secs_f64();

    all_pairs.par_sort_unstable_by_key(|&(hash, _)| hash);
    let t2 = t.elapsed().as_secs_f64();

    // Step 2: Generate edges from hash runs
    let mut edges: Vec<(u32, u32)> = Vec::new();
    let mut skipped_repetitive = 0usize;
    let mut i = 0;
    while i < all_pairs.len() {
        let hash = all_pairs[i].0;
        let run_start = i;
        while i < all_pairs.len() && all_pairs[i].0 == hash {
            i += 1;
        }
        let run_len = i - run_start;
        if run_len < 2 { continue; }
        if run_len > UF_MAX_BUCKET {
            skipped_repetitive += 1;
            continue;
        }
        // Generate all edges (a < b) within run
        for j in run_start..i {
            for k in (j + 1)..i {
                let a = all_pairs[j].1;
                let b = all_pairs[k].1;
                let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                edges.push((lo, hi));
            }
        }
    }
    drop(all_pairs);
    let t3 = t.elapsed().as_secs_f64();

    // Step 3: Sort edges and count duplicates
    edges.par_sort_unstable();
    let t4 = t.elapsed().as_secs_f64();

    // Step 4: Union pairs with count ≥ 2
    let mut uf = UnionFind::new(n);
    let mut unions = 0usize;
    let mut ei = 0;
    while ei < edges.len() {
        let edge = edges[ei];
        let edge_start = ei;
        while ei < edges.len() && edges[ei] == edge {
            ei += 1;
        }
        if ei - edge_start >= 2 {
            uf.union(edge.0, edge.1);
            unions += 1;
        }
    }
    drop(edges);

    let mut roots = Vec::with_capacity(n);
    for i in 0..n as u32 { roots.push(uf.find(i)); }
    drop(uf);

    let (result, grouped, groups) = emit_groups_in_order(n, &roots);

    info!("  Selective UF (≥2 shared): {} in {} groups ({:.1}%), {} unions, {} rep-skipped  \
           [min={:.1}s sort1={:.1}s edges={:.1}s sort2={:.1}s total={:.2}s]",
        grouped, groups, 100.0 * grouped as f64 / n as f64,
        unions, skipped_repetitive,
        t1, t2 - t1, t3 - t2, t4 - t3, t.elapsed().as_secs_f64());

    result
}

/// Strategy C: UF with all minimizers (original, too aggressive).
fn reorder_chunk_minimizer_uf(records: &[crate::io::FastqRecord]) -> ReorderResult {
    use rayon::prelude::*;

    let n = records.len();
    if n == 0 {
        return ReorderResult { order: Vec::new() };
    }

    let t = Instant::now();
    let mut all_pairs = compute_minimizer_pairs(records);
    let t1 = t.elapsed().as_secs_f64();

    all_pairs.par_sort_unstable_by_key(|&(hash, _)| hash);
    let t2 = t.elapsed().as_secs_f64();

    let mut uf = UnionFind::new(n);
    let mut unions = 0usize;
    let mut skipped = 0usize;
    let mut i = 0;
    while i < all_pairs.len() {
        let hash = all_pairs[i].0;
        let run_start = i;
        while i < all_pairs.len() && all_pairs[i].0 == hash { i += 1; }
        let run_len = i - run_start;
        if run_len < 2 { continue; }
        if run_len > UF_MAX_BUCKET { skipped += 1; continue; }
        let first = all_pairs[run_start].1;
        for j in (run_start + 1)..i {
            uf.union(first, all_pairs[j].1);
            unions += 1;
        }
    }
    drop(all_pairs);
    let t3 = t.elapsed().as_secs_f64();

    let mut roots = Vec::with_capacity(n);
    for i in 0..n as u32 { roots.push(uf.find(i)); }
    drop(uf);

    let (result, grouped, groups) = emit_groups_in_order(n, &roots);

    info!("  UF-all minimizer: {} in {} groups ({:.1}%), {} unions, {} rep-skipped  \
           [min={:.1}s sort={:.1}s uf={:.1}s total={:.2}s]",
        grouped, groups, 100.0 * grouped as f64 / n as f64,
        unions, skipped, t1, t2 - t1, t3 - t2, t.elapsed().as_secs_f64());

    result
}

// ── Delta encoding helpers ──────────────────────────────────────────────

/// Lookup table for base → delta code (branchless): 0=match, 1=A, 2=C, 3=G, 4=T, 5=N.
static BASE_DELTA_LUT: [u8; 256] = {
    let mut t = [5u8; 256];
    t[b'A' as usize] = 1; t[b'a' as usize] = 1;
    t[b'C' as usize] = 2; t[b'c' as usize] = 2;
    t[b'G' as usize] = 3; t[b'g' as usize] = 3;
    t[b'T' as usize] = 4; t[b't' as usize] = 4;
    t
};

/// Encode base as delta code: 0=match, 1=A, 2=C, 3=G, 4=T, 5=N.
#[inline]
fn base_to_delta_code(b: u8) -> u8 {
    BASE_DELTA_LUT[b as usize]
}

#[inline]
fn delta_code_to_base(c: u8) -> u8 {
    match c { 1 => b'A', 2 => b'C', 3 => b'G', 4 => b'T', _ => b'N' }
}

// ── Per-column arithmetic coding for sequences ──────────────────────────
//
// Predicts base[read_i][pos_j] from base[read_{i-1}][pos_j] (cross-read context).
// For grouped reads sharing center hashes, consecutive bases at the same position
// match ~90% of the time → very few bits per base. Singletons get ~2 bits/base.

const SEQ_N_SYMBOLS: usize = 5; // A=0, C=1, G=2, T=3, N=4
const SEQ_N_PREV: usize = 6;    // prev base codes: A,C,G,T,N,START
const SEQ_MAX_POS: usize = 512;  // max positions supported
const SEQ_ARITH_RESCALE: u32 = 1 << 15;

/// Marker byte at start of chunk blob to distinguish arithmetic from BSC format.
const ARITH_MARKER: u8 = 0xFE;

#[inline]
fn seq_code_to_base(c: u8) -> u8 {
    match c { 0 => b'A', 1 => b'C', 2 => b'G', 3 => b'T', _ => b'N' }
}

// ── Range coder (LZMA-style, from quality_ctx.rs pattern) ────────────────

const SEQ_RC_TOP: u32 = 1 << 24;

struct SeqRangeDecoder<'a> {
    range: u32,
    code: u32,
    input: &'a [u8],
    pos: usize,
}

impl<'a> SeqRangeDecoder<'a> {
    fn new(input: &'a [u8]) -> Self {
        let mut dec = Self { range: 0xFFFF_FFFF, code: 0, input, pos: 0 };
        if !input.is_empty() { dec.pos = 1; }
        for _ in 0..4 { dec.code = (dec.code << 8) | dec.next_byte() as u32; }
        dec
    }

    #[inline(always)]
    fn next_byte(&mut self) -> u8 {
        if self.pos < self.input.len() {
            let b = self.input[self.pos];
            self.pos += 1;
            b
        } else { 0 }
    }

    #[inline(always)]
    fn normalize(&mut self) {
        while self.range < SEQ_RC_TOP {
            self.code = (self.code << 8) | self.next_byte() as u32;
            self.range <<= 8;
        }
    }

    #[inline(always)]
    fn decode(&mut self, cum_freqs: &[u32; SEQ_N_SYMBOLS + 1], total: u32) -> usize {
        let r = self.range / total;
        let offset = (self.code / r).min(total - 1);
        let mut sym = 0;
        while sym + 1 < SEQ_N_SYMBOLS && cum_freqs[sym + 1] <= offset { sym += 1; }
        let cum = cum_freqs[sym];
        let freq = cum_freqs[sym + 1] - cum;
        self.code -= cum * r;
        if cum + freq < total { self.range = r * freq; }
        else { self.range -= r * cum; }
        self.normalize();
        sym
    }
}

struct SeqModel {
    cum_freqs: [u32; SEQ_N_SYMBOLS + 1],
    total: u32,
}

impl SeqModel {
    fn new() -> Self {
        let mut cum = [0u32; SEQ_N_SYMBOLS + 1];
        for i in 0..=SEQ_N_SYMBOLS { cum[i] = (i * 4) as u32; }
        Self { cum_freqs: cum, total: (SEQ_N_SYMBOLS * 4) as u32 }
    }

    #[inline(always)]
    fn update(&mut self, sym: usize) {
        for i in (sym + 1)..=SEQ_N_SYMBOLS {
            self.cum_freqs[i] += 1;
        }
        self.total += 1;
        if self.total >= SEQ_ARITH_RESCALE {
            self.rescale();
        }
    }

    fn rescale(&mut self) {
        let mut cum = 0u32;
        self.cum_freqs[0] = 0;
        for i in 0..SEQ_N_SYMBOLS {
            let freq = self.cum_freqs[i + 1] - self.cum_freqs[i];
            let new_freq = (freq >> 1).max(1);
            cum += new_freq;
            self.cum_freqs[i + 1] = cum;
        }
        self.total = cum;
    }
}

/// Decode arithmetic-coded sequences.
fn decode_sequences_arith(
    arith_data: &[u8],
    read_lengths: &[usize],
) -> Result<Vec<Vec<u8>>> {
    let num_reads = read_lengths.len();
    let max_len = read_lengths.iter().copied().max().unwrap_or(0);
    let max_pos = max_len.min(SEQ_MAX_POS);
    let num_contexts = max_pos * SEQ_N_PREV;

    let mut prev_col = vec![5u8; max_len];
    let mut models: Vec<SeqModel> = (0..num_contexts).map(|_| SeqModel::new()).collect();
    let mut decoder = SeqRangeDecoder::new(arith_data);

    let mut sequences = Vec::with_capacity(num_reads);

    for i in 0..num_reads {
        let read_len = read_lengths[i];
        let mut seq = Vec::with_capacity(read_len);

        for j in 0..read_len {
            let pos = j.min(max_pos - 1);
            let prev = prev_col[j] as usize;
            let ctx_idx = pos * SEQ_N_PREV + prev;

            let base = decoder.decode(&models[ctx_idx].cum_freqs, models[ctx_idx].total);
            models[ctx_idx].update(base);

            seq.push(seq_code_to_base(base as u8));
            prev_col[j] = base as u8;
        }

        sequences.push(seq);
    }

    Ok(sequences)
}

/// Decompress BSC blocks from blob format [num_blocks: u32, [block_len: u32, data]* ]
fn bsc_decompress_from_blob(data: &[u8], off: &mut usize) -> Result<Vec<u8>> {
    if *off + 4 > data.len() {
        anyhow::bail!("truncated BSC blob");
    }
    let num_blocks = super::read_le_u32(data, *off)? as usize;
    *off += 4;

    let mut result = Vec::new();
    for _ in 0..num_blocks {
        if *off + 4 > data.len() {
            anyhow::bail!("truncated BSC block length");
        }
        let block_len = super::read_le_u32(data, *off)? as usize;
        *off += 4;
        if *off + block_len > data.len() {
            anyhow::bail!("truncated BSC block data");
        }
        let decompressed = bsc::decompress_mt(&data[*off..*off + block_len])?;
        result.extend_from_slice(&decompressed);
        *off += block_len;
    }

    Ok(result)
}

// ── Compression ─────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
enum HarcMode {
    /// encoding_type=9: sort + raw BSC
    ReorderLocal,
    /// encoding_type=8: sort + previous-read delta
    Delta,
}

impl HarcMode {
    fn encoding_type(self) -> u8 {
        match self { HarcMode::ReorderLocal => 9, HarcMode::Delta => 8 }
    }

    fn name(self) -> &'static str {
        match self {
            HarcMode::ReorderLocal => "reorder-local",
            HarcMode::Delta => "HARC delta",
        }
    }
}

struct ChunkResult {
    h_blocks: Vec<Vec<u8>>,
    q_blocks: Vec<Vec<u8>>,
    seq_data: Vec<u8>,
    num_reads: usize,
}

/// Quality input: packed stream for BSC quality path.
enum QualInput {
    None,
    Bsc(Vec<u8>),
}

/// Compress quality_ctx blocks from borrowed references (avoids cloning record data).
fn compress_quality_ctx_refs(qual_refs: &[&[u8]], seq_refs: &[&[u8]], sub_block_reads: usize) -> Result<Vec<Vec<u8>>> {
    use rayon::prelude::*;
    let n = qual_refs.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    if n <= sub_block_reads {
        let ctx_blob = quality_ctx::compress_qualities_ctx(qual_refs, seq_refs)?;
        return Ok(vec![ctx_blob]);
    }
    let num_blocks = (n + sub_block_reads - 1) / sub_block_reads;
    (0..num_blocks)
        .into_par_iter()
        .map(|i| {
            let start = i * sub_block_reads;
            let end = (start + sub_block_reads).min(n);
            quality_ctx::compress_qualities_ctx(
                &qual_refs[start..end], &seq_refs[start..end],
            )
        })
        .collect::<Result<Vec<_>>>()
}

fn compress_chunk(
    records: Vec<crate::io::FastqRecord>,
    mode: HarcMode,
    quality_mode: QualityMode,
    quality_binning: QualityBinning,
    no_quality: bool,
    bsc_static: bool,
    use_quality_ctx: bool,
    quality_sub_block: usize,
) -> Result<ChunkResult> {
    let n = records.len();
    let compress_fn: fn(&[u8]) -> Result<Vec<u8>> = if bsc_static {
        bsc::compress
    } else {
        bsc::compress_adaptive_no_lzp
    };

    // Reorder strategy: identity (no reorder) is best — natural input order
    // has better BWT locality than any hash-based grouping we've tried.
    // QZ_REORDER env var allows benchmarking alternatives.
    let reorder = match mode {
        HarcMode::ReorderLocal => {
            let strategy = std::env::var("QZ_REORDER").unwrap_or_default();
            match strategy.as_str() {
                "center" => reorder_chunk(&records),
                "selective" => reorder_selective_uf(&records),
                "uf" => reorder_chunk_minimizer_uf(&records),
                _ => reorder_identity(&records),
            }
        }
        HarcMode::Delta => reorder_chunk(&records),
    };

    // Build order stream (u32 LE)
    let mut order_stream = Vec::with_capacity(n * 4);
    for &idx in &reorder.order {
        order_stream.extend_from_slice(&idx.to_le_bytes());
    }

    // Build read lengths stream
    let mut rl_stream = Vec::new();
    for &idx in &reorder.order {
        dna_utils::write_varint(&mut rl_stream, records[idx as usize].sequence.len() as u64);
    }

    // Build all raw streams from records before dropping them
    let mut header_stream = Vec::with_capacity(n * 64);
    for &idx in &reorder.order {
        let record = &records[idx as usize];
        dna_utils::write_varint(&mut header_stream, record.id.len() as u64);
        header_stream.extend_from_slice(&record.id);
    }

    // Collect quality inputs for BSC path (quality_ctx borrows directly from records)
    let qual_input = if no_quality || use_quality_ctx {
        QualInput::None
    } else {
        let mut qual_stream = Vec::with_capacity(n * 136);
        for &idx in &reorder.order {
            if let Some(ref qual) = records[idx as usize].quality {
                let quantized = quality::quantize_quality(qual, quality_mode);
                dna_utils::write_varint(&mut qual_stream, quantized.len() as u64);
                let packed = columnar::pack_qualities(&quantized, quality_binning);
                qual_stream.extend_from_slice(&packed);
            }
        }
        QualInput::Bsc(qual_stream)
    };

    // Build seq_data based on mode
    match mode {
        HarcMode::ReorderLocal => {
            let mut seq_stream = Vec::with_capacity(n * 150);
            for &idx in &reorder.order {
                seq_stream.extend_from_slice(&records[idx as usize].sequence);
            }

            // Parallel: quality compression || (BSC headers + BSC sequences)
            // Records kept alive for quality_ctx to borrow (no clone needed)
            let (q_result, hs_result) = rayon::join(
                || -> Result<Vec<Vec<u8>>> {
                    if no_quality {
                        Ok(Vec::new())
                    } else if use_quality_ctx {
                        // Borrow directly from records — avoids ~3 GB clone for 5M reads
                        let qual_refs: Vec<&[u8]> = reorder.order.iter()
                            .map(|&idx| records[idx as usize].quality.as_deref().unwrap_or(&[]))
                            .collect();
                        let seq_refs: Vec<&[u8]> = reorder.order.iter()
                            .map(|&idx| records[idx as usize].sequence.as_slice())
                            .collect();
                        compress_quality_ctx_refs(&qual_refs, &seq_refs, quality_sub_block)
                    } else {
                        match qual_input {
                            QualInput::Bsc(ref qual_stream) => compress_stream_to_bsc_blocks(qual_stream, bsc_static),
                            QualInput::None => Ok(Vec::new()),
                        }
                    }
                },
                || -> Result<(Vec<Vec<u8>>, Vec<u8>)> {
                    let h_blocks = compress_stream_to_bsc_blocks(&header_stream, bsc_static)?;
                    let seq_streams = vec![order_stream, seq_stream, rl_stream];
                    let num_seq_streams = seq_streams.len() as u8;

                    // BSC compress each seq stream. Split large streams into
                    // <=750MB blocks to stay within libsais's working limits.
                    // LZP enabled: reduces data before BWT even for DNA (~25% hit
                    // rate), making the expensive BWT sort faster overall.
                    const MAX_BSC_BLOCK: usize = 750 * 1024 * 1024;
                    use rayon::prelude::*;
                    let compressed_blocks: Vec<Vec<Vec<u8>>> = seq_streams.into_par_iter().map(|data| {
                        if data.is_empty() {
                            Ok(Vec::new())
                        } else if data.len() <= MAX_BSC_BLOCK {
                            let compressed = bsc::compress_adaptive_mt(&data)?;
                            Ok(vec![compressed])
                        } else {
                            let blocks: Vec<Vec<u8>> = data.par_chunks(MAX_BSC_BLOCK)
                                .map(|chunk| bsc::compress_adaptive_mt(chunk))
                                .collect::<Result<Vec<_>>>()?;
                            Ok(blocks)
                        }
                    }).collect::<Result<Vec<_>>>()?;

                    let mut blob = Vec::new();
                    blob.push(num_seq_streams);
                    for blocks in &compressed_blocks {
                        blob.extend_from_slice(&(blocks.len() as u32).to_le_bytes());
                        for block in blocks {
                            blob.extend_from_slice(&(block.len() as u32).to_le_bytes());
                            blob.extend_from_slice(block);
                        }
                    }
                    Ok((h_blocks, blob))
                },
            );

            let (h_blocks, seq_data) = hs_result?;
            return Ok(ChunkResult {
                h_blocks,
                q_blocks: q_result?,
                seq_data,
                num_reads: n,
            });
        }
        HarcMode::Delta => {
            // Delta encode every read against the previous.
            let mut delta_stream = Vec::with_capacity(n * 150);
            let mut prev_seq: Vec<u8> = Vec::new();

            for i in 0..n {
                let orig_idx = reorder.order[i] as usize;
                let cur_seq = &records[orig_idx].sequence;
                let read_len = cur_seq.len();

                if prev_seq.is_empty() {
                    for &b in cur_seq {
                        delta_stream.push(base_to_delta_code(b));
                    }
                } else {
                    let overlap = read_len.min(prev_seq.len());
                    for j in 0..read_len {
                        if j < overlap && cur_seq[j] == prev_seq[j] {
                            delta_stream.push(0);
                        } else {
                            delta_stream.push(base_to_delta_code(cur_seq[j]));
                        }
                    }
                }

                prev_seq.clear();
                prev_seq.extend_from_slice(cur_seq);
            }

            let seq_streams = vec![order_stream, delta_stream, rl_stream];
            let num_seq_streams = seq_streams.len() as u8;

            // BSC compress with large blocks (250MB)
            const HARC_BSC_BLOCK_SIZE: usize = 250 * 1024 * 1024;

            // Parallel: quality compression || (BSC headers + BSC sequences)
            // Records kept alive for quality_ctx to borrow (no clone needed)
            let (q_result, hs_result) = rayon::join(
                || -> Result<Vec<Vec<u8>>> {
                    if no_quality {
                        Ok(Vec::new())
                    } else if use_quality_ctx {
                        let qual_refs: Vec<&[u8]> = reorder.order.iter()
                            .map(|&idx| records[idx as usize].quality.as_deref().unwrap_or(&[]))
                            .collect();
                        let seq_refs: Vec<&[u8]> = reorder.order.iter()
                            .map(|&idx| records[idx as usize].sequence.as_slice())
                            .collect();
                        compress_quality_ctx_refs(&qual_refs, &seq_refs, quality_sub_block)
                    } else {
                        match qual_input {
                            QualInput::Bsc(ref qual_stream) => compress_stream_to_bsc_blocks(qual_stream, bsc_static),
                            QualInput::None => Ok(Vec::new()),
                        }
                    }
                },
                || -> Result<(Vec<Vec<u8>>, Vec<u8>)> {
                    let h_blocks = compress_stream_to_bsc_blocks(&header_stream, bsc_static)?;
                    let compressed_blocks: Vec<Vec<Vec<u8>>> = {
                        use rayon::prelude::*;
                        seq_streams.into_iter().map(|data| {
                            if data.is_empty() {
                                Ok(Vec::new())
                            } else if data.len() <= HARC_BSC_BLOCK_SIZE {
                                let compressed = compress_fn(&data)?;
                                Ok(vec![compressed])
                            } else {
                                let mut blocks: Vec<Vec<u8>> = data.par_chunks(HARC_BSC_BLOCK_SIZE)
                                    .map(|chunk| compress_fn(chunk))
                                    .collect::<Result<Vec<_>>>()?;
                                for b in &mut blocks { b.shrink_to_fit(); }
                                Ok(blocks)
                            }
                        }).collect::<Result<Vec<_>>>()?
                    };

                    let mut blob = Vec::new();
                    blob.push(num_seq_streams);
                    for blocks in &compressed_blocks {
                        blob.extend_from_slice(&(blocks.len() as u32).to_le_bytes());
                        for block in blocks {
                            blob.extend_from_slice(&(block.len() as u32).to_le_bytes());
                            blob.extend_from_slice(block);
                        }
                    }
                    Ok((h_blocks, blob))
                },
            );

            let (h_blocks, seq_data) = hs_result?;
            return Ok(ChunkResult {
                h_blocks,
                q_blocks: q_result?,
                seq_data,
                num_reads: n,
            });
        }
    }
}

// ── Public API: compress ─────────────────────────────────────────────────

fn compress_inner(args: &CompressConfig, mode: HarcMode) -> Result<()> {
    compress_inner_with(args, mode, CHUNK_SIZE, 2, 500_000)
}

fn compress_inner_with(args: &CompressConfig, mode: HarcMode, chunk_size: usize, max_parallel_chunks: usize, quality_sub_block: usize) -> Result<()> {
    use std::io::{Write, BufWriter};
    use crate::cli::QualityCompressor;

    let start_time = Instant::now();
    let bsc_static = args.advanced.bsc_static;
    let no_quality = args.no_quality;
    let quality_mode = args.quality_mode;
    let quality_binning = if no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(quality_mode)
    };

    // Auto-select quality_ctx for lossless mode with fixed-length reads
    // (same logic as streaming path: explicit request or large dataset)
    let use_quality_ctx = !no_quality
        && quality_mode == QualityMode::Lossless
        && (args.advanced.quality_compressor == QualityCompressor::QualityCtx
            || true); // always auto-select for harc path (chunks are 5M reads >> 100K threshold)
    if use_quality_ctx {
        info!("Using context-adaptive quality compression (quality_ctx)");
    }

    info!("{} compression: {} encoding",
        mode.name(),
        match mode { HarcMode::ReorderLocal => "raw BSC", HarcMode::Delta => "delta" });

    let working_dir = &args.working_dir;
    let h_tmp_path = working_dir.join(".qz_harc_h.tmp");
    let s_tmp_path = working_dir.join(".qz_harc_s.tmp");
    let q_tmp_path = working_dir.join(".qz_harc_q.tmp");

    struct TmpCleanup(Vec<std::path::PathBuf>);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            for p in &self.0 {
                let _ = std::fs::remove_file(p);
            }
        }
    }
    let _cleanup = TmpCleanup(vec![
        h_tmp_path.clone(), s_tmp_path.clone(), q_tmp_path.clone(),
    ]);

    let mut h_tmp = BufWriter::new(std::fs::File::create(&h_tmp_path)?);
    let mut s_tmp = BufWriter::new(std::fs::File::create(&s_tmp_path)?);
    let mut q_tmp = BufWriter::new(std::fs::File::create(&q_tmp_path)?);

    let mut h_num_blocks: u32 = 0;
    let mut q_num_blocks: u32 = 0;
    let mut s_chunks: u32 = 0;
    let mut num_reads: usize = 0;
    let mut original_size: usize = 0;

    let mut reader = crate::io::FastqReader::from_path_or_stdin(&args.input[0], args.fasta)?;

    // Read first chunk to check variable-length reads
    let (first_records, _, first_orig) = read_chunk_records(&mut reader, chunk_size)?;
    if first_records.is_empty() {
        // Empty file — write minimal archive and return
        h_tmp.flush()?; s_tmp.flush()?; q_tmp.flush()?;
        drop(h_tmp); drop(s_tmp); drop(q_tmp);
        use std::io::Write as _;
        let mut out: Box<dyn Write> = if crate::cli::is_stdio_path(&args.output) {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(std::fs::File::create(&args.output)?))
        };
        // v2 prefix + empty 52-byte body
        out.write_all(&ARCHIVE_MAGIC)?;
        out.write_all(&[ARCHIVE_VERSION, 0])?;
        out.write_all(&60u32.to_le_bytes())?; // header_size = 8 + 52
        out.write_all(&[mode.encoding_type()])?;
        out.write_all(&[0u8; 51])?;
        out.flush()?;
        return Ok(());
    }

    let quality_compressor_used = if use_quality_ctx {
        QualityCompressor::QualityCtx
    } else {
        QualityCompressor::Bsc
    };

    // Process up to max_parallel_chunks simultaneously to overlap BWT work across chunks,
    // while keeping memory bounded.

    let mut chunk_idx = 0usize;
    let mut pending: Vec<(Vec<crate::io::FastqRecord>, usize)> = vec![(first_records, first_orig)];

    loop {
        // Fill batch: read up to max_parallel_chunks
        while pending.len() < max_parallel_chunks {
            let (records, _, orig) = read_chunk_records(&mut reader, chunk_size)?;
            if records.is_empty() { break; }
            pending.push((records, orig));
        }
        if pending.is_empty() { break; }

        let batch = std::mem::take(&mut pending);
        let batch_len = batch.len();
        let start_idx = chunk_idx;

        if batch_len > 1 {
            info!("Compressing chunks {}..{} in parallel ({} chunks)",
                start_idx, start_idx + batch_len - 1, batch_len);
        } else {
            info!("Chunk {}: {} reads", start_idx, batch[0].0.len());
        }

        // Compress batch in parallel (read next batch on main thread)
        let (next_batch, results) = std::thread::scope(|scope| {
            let compress_handle = scope.spawn(|| -> Vec<Result<(ChunkResult, usize)>> {
                use rayon::prelude::*;
                batch.into_par_iter().enumerate().map(|(i, (records, orig))| {
                    let t = Instant::now();
                    let n = records.len();
                    let result = compress_chunk(
                        records, mode, quality_mode, quality_binning, no_quality,
                        bsc_static, use_quality_ctx, quality_sub_block,
                    )?;
                    info!("  Chunk {}: {} reads, {:.2}s",
                        start_idx + i, n, t.elapsed().as_secs_f64());
                    Ok((result, orig))
                }).collect()
            });

            // Read next batch while compressing
            let mut next = Vec::new();
            for _ in 0..max_parallel_chunks {
                match read_chunk_records(&mut reader, chunk_size) {
                    Ok((records, _, orig)) => {
                        if records.is_empty() { break; }
                        next.push((records, orig));
                    }
                    Err(e) => { next.push((Vec::new(), 0)); eprintln!("Read error: {e}"); break; }
                }
            }

            let compressed = match compress_handle.join() {
                Ok(v) => v,
                Err(e) => std::panic::resume_unwind(e),
            };
            (next, compressed)
        });

        // Write results in order
        for result in results {
            let (chunk, orig) = result?;
            h_num_blocks += write_blocks_to_tmp(chunk.h_blocks, &mut h_tmp)?;
            if !chunk.q_blocks.is_empty() {
                q_num_blocks += write_blocks_to_tmp(chunk.q_blocks, &mut q_tmp)?;
            }
            {
                use std::io::Write as _;
                s_tmp.write_all(&(chunk.seq_data.len() as u32).to_le_bytes())?;
                s_tmp.write_all(&chunk.seq_data)?;
                s_chunks += 1;
            }
            num_reads += chunk.num_reads;
            original_size += orig;
        }

        chunk_idx += batch_len;
        pending = next_batch;
    }

    h_tmp.flush()?;
    s_tmp.flush()?;
    q_tmp.flush()?;
    drop(h_tmp);
    drop(s_tmp);
    drop(q_tmp);

    info!("Processed {} reads in {} chunks", num_reads, chunk_idx);

    let h_tmp_size = std::fs::metadata(&h_tmp_path)?.len() as usize;
    let s_tmp_size = std::fs::metadata(&s_tmp_path)?.len() as usize;
    let q_tmp_size = std::fs::metadata(&q_tmp_path)?.len() as usize;

    let h_size = if h_num_blocks > 0 { 4 + h_tmp_size } else { 0 };
    let s_size = 4 + s_tmp_size;
    let q_size = if q_num_blocks > 0 { 4 + q_tmp_size } else { 0 };

    info!("Writing output file...");
    let mut out: Box<dyn std::io::Write> = if crate::cli::is_stdio_path(&args.output) {
        Box::new(BufWriter::new(std::io::stdout().lock()))
    } else {
        Box::new(BufWriter::new(std::fs::File::create(&args.output)?))
    };

    // v2 prefix: magic + version + reserved + header_size
    let header_size: u32 = (V2_PREFIX_SIZE + 52) as u32; // 8 + 52 = 60
    out.write_all(&ARCHIVE_MAGIC)?;
    out.write_all(&[ARCHIVE_VERSION, 0])?;
    out.write_all(&header_size.to_le_bytes())?;

    // Header body
    out.write_all(&[mode.encoding_type()])?;
    out.write_all(&[0u8])?; // flags
    out.write_all(&[binning_to_code(quality_binning)])?;
    out.write_all(&[compressor_to_code(quality_compressor_used)])?;
    out.write_all(&[seq_compressor_to_code(SequenceCompressor::Bsc)])?;
    out.write_all(&[header_compressor_to_code(HeaderCompressor::Bsc)])?;
    out.write_all(&[0u8; 3])?;
    out.write_all(&0u16.to_le_bytes())?;
    out.write_all(&[0u8])?;

    out.write_all(&(num_reads as u64).to_le_bytes())?;
    out.write_all(&(h_size as u64).to_le_bytes())?;
    out.write_all(&(s_size as u64).to_le_bytes())?;
    out.write_all(&0u64.to_le_bytes())?; // nmasks_len
    out.write_all(&(q_size as u64).to_le_bytes())?;

    if h_num_blocks > 0 {
        out.write_all(&h_num_blocks.to_le_bytes())?;
        let mut f = std::fs::File::open(&h_tmp_path)?;
        std::io::copy(&mut f, &mut out)?;
    }

    out.write_all(&s_chunks.to_le_bytes())?;
    {
        let mut f = std::fs::File::open(&s_tmp_path)?;
        std::io::copy(&mut f, &mut out)?;
    }

    if q_num_blocks > 0 {
        out.write_all(&q_num_blocks.to_le_bytes())?;
        let mut f = std::fs::File::open(&q_tmp_path)?;
        std::io::copy(&mut f, &mut out)?;
    }

    out.flush()?;

    let total = h_size + s_size + q_size + header_size as usize;
    let elapsed = start_time.elapsed();
    info!("Compression completed in {:.2}s", elapsed.as_secs_f64());
    info!("Original size: {} bytes", original_size);
    info!("Compressed size: {} bytes", total);
    info!("Compression ratio: {:.2}x", original_size as f64 / total as f64);
    info!("Stream breakdown:");
    info!("  Headers:    {} bytes ({:.1}%)", h_size, 100.0 * h_size as f64 / total as f64);
    info!("  Sequences:  {} bytes ({:.1}%)", s_size, 100.0 * s_size as f64 / total as f64);
    info!("  Qualities:  {} bytes ({:.1}%)", q_size, 100.0 * q_size as f64 / total as f64);

    Ok(())
}

pub(super) fn compress_harc(args: &CompressConfig) -> Result<()> {
    compress_inner(args, HarcMode::Delta)
}

pub(super) fn compress_reorder_local_with_level(args: &CompressConfig, level: UltraLevel) -> Result<()> {
    info!("Ultra level {}: chunk_size={}M, parallel_chunks={}, quality_sub_block={}K",
        level.level, level.chunk_size / 1_000_000, level.max_parallel_chunks, level.quality_sub_block / 1000);
    compress_inner_with(args, HarcMode::ReorderLocal, level.chunk_size, level.max_parallel_chunks, level.quality_sub_block)
}

// ── Decoding ────────────────────────────────────────────────────────────

/// Decode result: sequences in reordered order + permutation for restoring original order.
pub(super) struct HarcDecoded {
    pub sequences: Vec<Vec<u8>>,
    /// order[i] = original read index of the i-th read in reordered order.
    pub order: Vec<u32>,
}

pub(super) fn decode_harc_sequences(seq_region: &[u8], num_reads: usize) -> Result<HarcDecoded> {
    decode_chunked(seq_region, num_reads, true)
}

pub(super) fn decode_reorder_local(seq_region: &[u8], num_reads: usize) -> Result<HarcDecoded> {
    decode_chunked(seq_region, num_reads, false)
}

fn decode_chunked(seq_region: &[u8], num_reads: usize, is_delta: bool) -> Result<HarcDecoded> {
    use anyhow::Context;

    if seq_region.len() < 4 {
        anyhow::bail!("HARC: sequences region too small");
    }

    let num_chunks = super::read_le_u32(seq_region, 0)? as usize;
    let mut off = 4usize;

    // Phase 1: parse chunk boundaries (fast, sequential pointer math)
    let mut chunk_slices: Vec<&[u8]> = Vec::with_capacity(num_chunks);
    for chunk_idx in 0..num_chunks {
        if off + 4 > seq_region.len() {
            anyhow::bail!("HARC: truncated chunk {} length", chunk_idx);
        }
        let chunk_len = super::read_le_u32(seq_region, off)? as usize;
        off += 4;
        if off + chunk_len > seq_region.len() {
            anyhow::bail!("HARC: truncated chunk {} data", chunk_idx);
        }
        chunk_slices.push(&seq_region[off..off + chunk_len]);
        off += chunk_len;
    }

    // Phase 2: decode all chunks in parallel (inverse BWT is the bottleneck)
    use rayon::prelude::*;
    let decoded_chunks: Vec<Result<(Vec<Vec<u8>>, Vec<u32>)>> = chunk_slices
        .into_par_iter()
        .enumerate()
        .map(|(chunk_idx, chunk_data)| {
            if is_delta {
                decode_chunk_delta(chunk_data)
                    .with_context(|| format!("Failed to decode HARC delta chunk {}", chunk_idx))
            } else {
                decode_chunk_reorder_local(chunk_data)
                    .with_context(|| format!("Failed to decode reorder-local chunk {}", chunk_idx))
            }
        })
        .collect();

    // Phase 3: concatenate results with offset order indices
    let mut all_sequences = Vec::with_capacity(num_reads);
    let mut all_order = Vec::with_capacity(num_reads);
    for result in decoded_chunks {
        let (seqs, order) = result?;
        let chunk_offset = all_sequences.len() as u32;
        all_sequences.extend(seqs);
        all_order.extend(order.into_iter().map(|idx| idx + chunk_offset));
    }

    if all_sequences.len() != num_reads {
        anyhow::bail!("HARC: decoded {} sequences but expected {}",
            all_sequences.len(), num_reads);
    }

    Ok(HarcDecoded {
        sequences: all_sequences,
        order: all_order,
    })
}

/// Decompress sub-streams from blob format.
/// Format: [num_streams: u8] [for each: num_blocks: u32, [block_len: u32, block_data]* ]
fn decompress_streams(data: &[u8]) -> Result<(Vec<Vec<u8>>, usize)> {
    if data.is_empty() {
        return Ok((Vec::new(), 0));
    }

    let num_streams = data[0] as usize;
    let mut off = 1usize;

    // Parse all compressed blocks for all streams
    let mut stream_blocks: Vec<Vec<&[u8]>> = Vec::with_capacity(num_streams);
    for i in 0..num_streams {
        if off + 4 > data.len() {
            anyhow::bail!("HARC: truncated stream {} num_blocks", i);
        }
        let num_blocks = super::read_le_u32(data, off)? as usize;
        off += 4;

        let mut blocks = Vec::with_capacity(num_blocks);
        for b in 0..num_blocks {
            if off + 4 > data.len() {
                anyhow::bail!("HARC: truncated stream {} block {} length", i, b);
            }
            let block_len = super::read_le_u32(data, off)? as usize;
            off += 4;
            if off + block_len > data.len() {
                anyhow::bail!("HARC: truncated stream {} block {} data", i, b);
            }
            blocks.push(&data[off..off + block_len]);
            off += block_len;
        }
        stream_blocks.push(blocks);
    }

    // Decompress all blocks in parallel, then concatenate per stream
    use rayon::prelude::*;
    let decompressed: Vec<Vec<u8>> = stream_blocks.into_iter().map(|blocks| {
        if blocks.is_empty() {
            return Ok(Vec::new());
        }
        let dec_blocks: Vec<Vec<u8>> = if blocks.len() == 1 {
            // Single block: use BSC-internal MT (parallel inverse BWT)
            vec![bsc::decompress_mt(blocks[0])?]
        } else {
            // Multiple blocks: rayon parallelism across blocks
            blocks.par_iter().map(|block| {
                bsc::decompress(block)
            }).collect::<Result<Vec<_>>>()?
        };

        let total_len: usize = dec_blocks.iter().map(|b| b.len()).sum();
        let mut combined = Vec::with_capacity(total_len);
        for block in dec_blocks {
            combined.extend_from_slice(&block);
        }
        Ok(combined)
    }).collect::<Result<Vec<_>>>()?;

    Ok((decompressed, num_streams))
}

fn parse_order(order_data: &[u8]) -> Vec<u32> {
    let n = order_data.len() / 4;
    (0..n).map(|i| {
        super::read_le_u32(order_data, i * 4).unwrap_or(0)
    }).collect()
}

fn decode_chunk_reorder_local(data: &[u8]) -> Result<(Vec<Vec<u8>>, Vec<u32>)> {
    if data.is_empty() {
        anyhow::bail!("reorder-local: empty chunk data");
    }

    if data[0] == ARITH_MARKER {
        // Arithmetic mode: [ARITH_MARKER] [order BSC blob] [arith_len + data] [rl BSC blob]
        let mut off = 1;

        // Decompress order
        let order_data = bsc_decompress_from_blob(data, &mut off)?;
        let order = parse_order(&order_data);
        let num_reads = order.len();

        // Read arithmetic data
        if off + 4 > data.len() {
            anyhow::bail!("reorder-local: truncated arith length");
        }
        let arith_len = super::read_le_u32(data, off)? as usize;
        off += 4;
        if off + arith_len > data.len() {
            anyhow::bail!("reorder-local: truncated arith data");
        }
        let arith_data = &data[off..off + arith_len];
        off += arith_len;

        // Decompress read lengths
        let rl_data = bsc_decompress_from_blob(data, &mut off)?;
        let mut rl_off = 0;
        let mut read_lengths = Vec::with_capacity(num_reads);
        for _ in 0..num_reads {
            read_lengths.push(
                dna_utils::read_varint(&rl_data, &mut rl_off).unwrap_or(0) as usize
            );
        }

        // Decode sequences from arithmetic data
        let sequences = decode_sequences_arith(arith_data, &read_lengths)?;

        Ok((sequences, order))
    } else {
        // BSC mode (legacy): [num_streams=3] [BSC blocks...]
        let (streams, num_streams) = decompress_streams(data)?;
        if num_streams != 3 {
            anyhow::bail!("reorder-local: expected 3 streams, got {}", num_streams);
        }

        let order = parse_order(&streams[0]);
        let seq_data = &streams[1];
        let rl_data = &streams[2];
        let num_reads = order.len();

        let mut rl_off = 0;
        let mut read_lengths = Vec::with_capacity(num_reads);
        for _ in 0..num_reads {
            read_lengths.push(
                dna_utils::read_varint(rl_data, &mut rl_off).unwrap_or(0) as usize
            );
        }

        let mut seq_off = 0;
        let mut sequences = Vec::with_capacity(num_reads);
        for i in 0..num_reads {
            let len = read_lengths[i];
            if seq_off + len > seq_data.len() {
                anyhow::bail!("reorder-local: truncated sequence data at read {}", i);
            }
            sequences.push(seq_data[seq_off..seq_off + len].to_vec());
            seq_off += len;
        }

        Ok((sequences, order))
    }
}

fn decode_chunk_delta(data: &[u8]) -> Result<(Vec<Vec<u8>>, Vec<u32>)> {
    let (streams, num_streams) = decompress_streams(data)?;
    if num_streams != 3 {
        anyhow::bail!("HARC delta: expected 3 streams, got {}", num_streams);
    }

    let order = parse_order(&streams[0]);
    let delta_data = &streams[1];
    let rl_data = &streams[2];
    let num_reads = order.len();

    // Parse read lengths
    let mut rl_off = 0;
    let mut read_lengths = Vec::with_capacity(num_reads);
    for _ in 0..num_reads {
        read_lengths.push(
            dna_utils::read_varint(rl_data, &mut rl_off).unwrap_or(0) as usize
        );
    }

    // Decode: every read is delta-encoded against the previous
    let mut sequences = Vec::with_capacity(num_reads);
    let mut prev_seq: Vec<u8> = Vec::new();
    let mut delta_off = 0;

    for i in 0..num_reads {
        let read_len = read_lengths[i];
        let mut seq = Vec::with_capacity(read_len);

        if prev_seq.is_empty() {
            for _ in 0..read_len {
                if delta_off < delta_data.len() {
                    seq.push(delta_code_to_base(delta_data[delta_off]));
                    delta_off += 1;
                } else {
                    seq.push(b'N');
                }
            }
        } else {
            let overlap = read_len.min(prev_seq.len());
            for j in 0..read_len {
                if delta_off >= delta_data.len() {
                    seq.push(b'N');
                    continue;
                }
                let db = delta_data[delta_off];
                delta_off += 1;
                if db == 0 && j < overlap {
                    seq.push(prev_seq[j]);
                } else {
                    seq.push(delta_code_to_base(db));
                }
            }
        }

        prev_seq.clear();
        prev_seq.extend_from_slice(&seq);
        sequences.push(seq);
    }

    Ok((sequences, order))
}
