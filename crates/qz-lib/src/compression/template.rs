/// Loose-template sequence compression
///
/// Key insight: references don't need to be true genomic assemblies.
/// Even "imperfect" templates (up to ~20% mismatch) reduce encoding cost
/// vs raw bases, because encoding mismatches is cheaper than encoding
/// every base from scratch.
///
/// Approach:
/// 1. Build de Bruijn graph with SMALL k (11-15) → dense graph even at low coverage
/// 2. Extract unitigs as "loose templates" (not true contigs)
/// 3. Map ALL reads with HIGH mismatch tolerance (≤30 mismatches out of 150bp = 20%)
/// 4. Delta-encode each read against its best-matching template
/// 5. BSC-compress all streams
///
/// This works even at 0.25x coverage because small k-mers appear in many
/// reads regardless of true genomic overlap.

use super::bsc;
use super::dna_utils::*;
use anyhow::Result;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

// ── Configurable parameters ─────────────────────────────────────────────────

/// Parameters for template-based compression
pub struct TemplateParams {
    pub k: usize,
    pub mismatch_threshold: usize,
    pub tip_trim_threshold: u16,
    /// Use parallel template extraction + reduced mapping threshold for speed
    pub fast: bool,
    /// Use syncmer-based k-mer selection instead of fixed-stride for template indexing
    pub use_syncmers: bool,
}

impl Default for TemplateParams {
    fn default() -> Self {
        Self {
            k: 13,
            mismatch_threshold: 30,    // ~20% of 150bp
            tip_trim_threshold: 1,
            fast: false,
            use_syncmers: false,
        }
    }
}

// ── Data structures ─────────────────────────────────────────────────────────

struct TemplateGraph {
    kmers: FxHashMap<u64, u16>,
    k: usize,
}

struct Template {
    sequence: Vec<u8>,
}

#[derive(Clone)]
struct ReadMapping {
    template_id: usize,
    position: u32,
    is_rc: bool,
}

struct TemplateKmerHit {
    template_id: u32,
    position: u32,
    is_forward: bool,
}

// ── Graph construction ──────────────────────────────────────────────────────

fn canonical_hash(seq: &[u8], k: usize) -> Option<u64> {
    let hash = kmer_to_hash(seq)?;
    let rc = reverse_complement_hash(hash, k);
    Some(hash.min(rc))
}

impl TemplateGraph {
    fn from_sequences(sequences: &[Vec<u8>], k: usize) -> Self {
        let chunk_size = (sequences.len() / rayon::current_num_threads().max(1)).max(1000);

        let local_maps: Vec<FxHashMap<u64, u16>> = sequences
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut local: FxHashMap<u64, u16> = FxHashMap::default();
                local.reserve(chunk.len() * (150 - k + 1));
                for seq in chunk {
                    if seq.len() < k { continue; }
                    for i in 0..=seq.len() - k {
                        if let Some(canon) = canonical_hash(&seq[i..i + k], k) {
                            let count = local.entry(canon).or_insert(0);
                            *count = count.saturating_add(1);
                        }
                    }
                }
                local
            })
            .collect();

        let estimated_size = local_maps.iter().map(|m| m.len()).max().unwrap_or(0);
        let mut kmers: FxHashMap<u64, u16> = FxHashMap::default();
        kmers.reserve(estimated_size);
        for local in local_maps {
            for (hash, count) in local {
                let entry = kmers.entry(hash).or_insert(0);
                *entry = entry.saturating_add(count);
            }
        }

        eprintln!("  Template graph: {} unique canonical {}-mers from {} reads",
                 kmers.len(), k, sequences.len());

        TemplateGraph { kmers, k }
    }

    fn tip_trim(&mut self, min_count: u16) {
        let before = self.kmers.len();
        self.kmers.retain(|_, count| *count > min_count);
        eprintln!("  Tip trimming (count <= {}): {} → {} k-mers",
                 min_count, before, self.kmers.len());
    }

    fn hash_to_sequence(hash: u64, k: usize) -> Vec<u8> {
        let mut seq = vec![0u8; k];
        let mut h = hash;
        for i in (0..k).rev() {
            seq[i] = idx_to_base(h as usize & 3);
            h >>= 2;
        }
        seq
    }

}

// ── Template (greedy walk) extraction ───────────────────────────────────────

#[inline(always)]
fn base_to_2bit(b: u8) -> u64 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// Single greedy walk using rolling fwd+rc hash: O(1) per step instead of O(k).
fn walk_from_seed(
    graph: &TemplateGraph,
    start_hash: u64,
    k: usize,
    mask: u64,
    rc_shift: usize, // = 2 * (k - 1)
    max_walk_length: usize,
    min_length: usize,
) -> Option<Template> {
    let start_seq = TemplateGraph::hash_to_sequence(start_hash, k);
    let start_seq = if kmer_to_hash(&start_seq) == Some(start_hash) {
        start_seq
    } else {
        reverse_complement(&TemplateGraph::hash_to_sequence(start_hash, k))
    };
    let start_fwd = kmer_to_hash(&start_seq)?;
    let start_rc = reverse_complement_hash(start_fwd, k);

    let mut local_visited: FxHashSet<u64> = FxHashSet::default();
    local_visited.insert(start_fwd.min(start_rc));

    // Walk RIGHT: rolling fwd+rc hash
    let mut right_bases: Vec<u8> = Vec::with_capacity(1000);
    let mut cur_fwd = start_fwd;
    let mut cur_rc = start_rc;
    while right_bases.len() < max_walk_length {
        let mut best_canon = 0u64;
        let mut best_base = 0u8;
        let mut best_count = 0u16;
        let mut best_fwd = 0u64;
        let mut best_rc = 0u64;

        for bval in 0u64..4 {
            let cand_fwd = ((cur_fwd << 2) | bval) & mask;
            let cand_rc = (cur_rc >> 2) | ((3 - bval) << rc_shift);
            let canon = cand_fwd.min(cand_rc);

            if let Some(&count) = graph.kmers.get(&canon) {
                if count > best_count {
                    best_count = count;
                    best_canon = canon;
                    best_base = idx_to_base(bval as usize);
                    best_fwd = cand_fwd;
                    best_rc = cand_rc;
                }
            }
        }

        if best_count == 0 || local_visited.contains(&best_canon) { break; }
        local_visited.insert(best_canon);
        right_bases.push(best_base);
        cur_fwd = best_fwd;
        cur_rc = best_rc;
    }

    // Walk LEFT: rolling fwd+rc hash
    let mut left_bases: Vec<u8> = Vec::with_capacity(1000);
    cur_fwd = start_fwd;
    cur_rc = start_rc;
    while left_bases.len() < max_walk_length {
        let mut best_canon = 0u64;
        let mut best_base = 0u8;
        let mut best_count = 0u16;
        let mut best_fwd = 0u64;
        let mut best_rc = 0u64;

        for bval in 0u64..4 {
            let cand_fwd = (bval << rc_shift) | (cur_fwd >> 2);
            let cand_rc = ((cur_rc << 2) | (3 - bval)) & mask;
            let canon = cand_fwd.min(cand_rc);

            if let Some(&count) = graph.kmers.get(&canon) {
                if count > best_count {
                    best_count = count;
                    best_canon = canon;
                    best_base = idx_to_base(bval as usize);
                    best_fwd = cand_fwd;
                    best_rc = cand_rc;
                }
            }
        }

        if best_count == 0 || local_visited.contains(&best_canon) { break; }
        local_visited.insert(best_canon);
        left_bases.push(best_base);
        cur_fwd = best_fwd;
        cur_rc = best_rc;
    }

    left_bases.reverse();
    let total_len = left_bases.len() + k + right_bases.len();

    if total_len >= min_length {
        let mut seq = Vec::with_capacity(total_len);
        seq.extend_from_slice(&left_bases);
        seq.extend_from_slice(&start_seq);
        seq.extend_from_slice(&right_bases);
        Some(Template { sequence: seq })
    } else {
        None
    }
}

/// Extract templates via sequential greedy walks with diversity tracking.
/// Uses rolling hash walks for speed + used_seeds set to ensure templates
/// cover distinct regions of k-mer space.
fn extract_templates(graph: &TemplateGraph, min_length: usize) -> Vec<Template> {
    let k = graph.k;
    let max_walk_length = 5_000;
    let max_templates = 50_000;
    let min_seed_count = 3u16;
    let mask: u64 = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

    // Pre-filter: only collect seeds with count ≥ 3 before sorting (was 52M → 985K)
    let mut seeds: Vec<(u64, u16)> = graph.kmers.iter()
        .filter(|&(_, c)| *c >= min_seed_count)
        .map(|(&h, &c)| (h, c))
        .collect();
    seeds.sort_unstable_by(|a, b| b.1.cmp(&a.1));

    eprintln!("  {} seeds with count ≥ {} (out of {} total k-mers)",
             seeds.len(), min_seed_count, graph.kmers.len());

    let rc_shift = 2 * (k - 1);
    let mut used_seeds: FxHashSet<u64> = FxHashSet::default();
    let mut templates: Vec<Template> = Vec::new();

    for &(start_hash, start_count) in &seeds {
        if templates.len() >= max_templates { break; }
        if start_count < min_seed_count { break; }
        if used_seeds.contains(&start_hash) { continue; }

        if let Some(template) = walk_from_seed(graph, start_hash, k, mask, rc_shift, max_walk_length, min_length) {
            // Mark template k-mers as used seeds using rolling fwd+rc hash
            let seq = &template.sequence;
            if seq.len() >= k {
                let mut fwd = 0u64;
                for i in 0..k {
                    fwd = (fwd << 2) | base_to_2bit(seq[i]);
                }
                fwd &= mask;
                let mut rc = reverse_complement_hash(fwd, k); // O(k) once
                used_seeds.insert(fwd.min(rc));

                for i in k..seq.len() {
                    let bval = base_to_2bit(seq[i]);
                    fwd = ((fwd << 2) | bval) & mask;
                    rc = (rc >> 2) | ((3 - bval) << rc_shift);
                    used_seeds.insert(fwd.min(rc));
                }
            }
            templates.push(template);
        }
    }

    templates.sort_unstable_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));

    let total_bases: usize = templates.iter().map(|t| t.sequence.len()).sum();
    let max_len = templates.first().map(|t| t.sequence.len()).unwrap_or(0);
    let long_templates = templates.iter().filter(|t| t.sequence.len() >= 150).count();
    eprintln!("  Extracted {} templates ({} total bases, longest={}, {} ≥ 150bp)",
             templates.len(), total_bases, max_len, long_templates);

    templates
}

/// Extend template tips through count=1 k-mers in the full graph.
/// After extracting from count≥2 filtered graph, templates are truncated at
/// singleton boundaries. This function extends them using the original graph.
fn extend_tips(templates: &mut [Template], graph: &TemplateGraph) -> (usize, usize) {
    let k = graph.k;
    let mask: u64 = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let rc_shift = 2 * (k - 1);
    let max_extension = 500;

    let mut extended_count = 0usize;
    let mut total_bases_added = 0usize;

    for template in templates.iter_mut() {
        let seq = &template.sequence;
        if seq.len() < k { continue; }

        let orig_len = seq.len();

        // Build visited set from existing template k-mers (avoid re-entering)
        let mut visited: FxHashSet<u64> = FxHashSet::default();
        let mut fwd = 0u64;
        for i in 0..k {
            fwd = (fwd << 2) | base_to_2bit(seq[i]);
        }
        fwd &= mask;
        let mut rc = reverse_complement_hash(fwd, k);
        visited.insert(fwd.min(rc));
        for i in k..seq.len() {
            let bval = base_to_2bit(seq[i]);
            fwd = ((fwd << 2) | bval) & mask;
            rc = (rc >> 2) | ((3 - bval) << rc_shift);
            visited.insert(fwd.min(rc));
        }

        // Extend RIGHT from last k-mer
        // fwd/rc are already set to the last k-mer from the loop above
        let mut right_ext: Vec<u8> = Vec::new();
        while right_ext.len() < max_extension {
            let mut best_count = 0u16;
            let mut best_base = 0u8;
            let mut best_fwd = 0u64;
            let mut best_rc = 0u64;
            let mut best_canon = 0u64;

            for bval in 0u64..4 {
                let cand_fwd = ((fwd << 2) | bval) & mask;
                let cand_rc = (rc >> 2) | ((3 - bval) << rc_shift);
                let canon = cand_fwd.min(cand_rc);
                if let Some(&count) = graph.kmers.get(&canon) {
                    if count > best_count {
                        best_count = count;
                        best_base = idx_to_base(bval as usize);
                        best_fwd = cand_fwd;
                        best_rc = cand_rc;
                        best_canon = canon;
                    }
                }
            }

            if best_count == 0 || visited.contains(&best_canon) { break; }
            visited.insert(best_canon);
            right_ext.push(best_base);
            fwd = best_fwd;
            rc = best_rc;
        }

        // Extend LEFT from first k-mer
        fwd = 0u64;
        for i in 0..k {
            fwd = (fwd << 2) | base_to_2bit(seq[i]);
        }
        fwd &= mask;
        rc = reverse_complement_hash(fwd, k);

        let mut left_ext: Vec<u8> = Vec::new();
        while left_ext.len() < max_extension {
            let mut best_count = 0u16;
            let mut best_base = 0u8;
            let mut best_fwd = 0u64;
            let mut best_rc = 0u64;
            let mut best_canon = 0u64;

            for bval in 0u64..4 {
                let cand_fwd = (bval << rc_shift) | (fwd >> 2);
                let cand_rc = ((rc << 2) | (3 - bval)) & mask;
                let canon = cand_fwd.min(cand_rc);
                if let Some(&count) = graph.kmers.get(&canon) {
                    if count > best_count {
                        best_count = count;
                        best_base = idx_to_base(bval as usize);
                        best_fwd = cand_fwd;
                        best_rc = cand_rc;
                        best_canon = canon;
                    }
                }
            }

            if best_count == 0 || visited.contains(&best_canon) { break; }
            visited.insert(best_canon);
            left_ext.push(best_base);
            fwd = best_fwd;
            rc = best_rc;
        }
        left_ext.reverse();

        if !left_ext.is_empty() || !right_ext.is_empty() {
            let mut new_seq = Vec::with_capacity(left_ext.len() + seq.len() + right_ext.len());
            new_seq.extend_from_slice(&left_ext);
            new_seq.extend_from_slice(seq);
            new_seq.extend_from_slice(&right_ext);
            total_bases_added += new_seq.len() - orig_len;
            template.sequence = new_seq;
            extended_count += 1;
        }
    }

    (extended_count, total_bases_added)
}

/// Error-correct templates using k-mer voting from the graph.
/// Finds positions where a k-mer has count=1 but a neighbor has count≥3,
/// tries all 3 alternative bases and picks the one with highest k-mer count.
fn error_correct_templates(templates: &mut [Template], graph: &TemplateGraph) -> usize {
    let k = graph.k;
    let mask: u64 = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let rc_shift = 2 * (k - 1);
    let mut total_corrections = 0usize;

    for template in templates.iter_mut() {
        if template.sequence.len() < k + 2 { continue; }

        // Phase 1: compute rolling hashes and counts (immutable borrow of sequence)
        let n_kmers = template.sequence.len() - k + 1;
        let mut hashes_fwd: Vec<u64> = Vec::with_capacity(n_kmers);
        let mut hashes_rc: Vec<u64> = Vec::with_capacity(n_kmers);
        let mut counts: Vec<u16> = Vec::with_capacity(n_kmers);

        let mut fwd = 0u64;
        for i in 0..k {
            fwd = (fwd << 2) | base_to_2bit(template.sequence[i]);
        }
        fwd &= mask;
        let mut rc = reverse_complement_hash(fwd, k);
        let canon = fwd.min(rc);
        hashes_fwd.push(fwd);
        hashes_rc.push(rc);
        counts.push(*graph.kmers.get(&canon).unwrap_or(&0));

        for i in k..template.sequence.len() {
            let bval = base_to_2bit(template.sequence[i]);
            fwd = ((fwd << 2) | bval) & mask;
            rc = (rc >> 2) | ((3 - bval) << rc_shift);
            let canon = fwd.min(rc);
            hashes_fwd.push(fwd);
            hashes_rc.push(rc);
            counts.push(*graph.kmers.get(&canon).unwrap_or(&0));
        }

        // Phase 2: collect corrections (no borrow of sequence)
        let mut corrections: Vec<(usize, u8)> = Vec::new();
        let mut corrected_pos: FxHashSet<usize> = FxHashSet::default();

        for pos in 0..n_kmers {
            if counts[pos] >= 2 { continue; }

            // Check if left neighbor has high count → error at pos+k-1 (trailing edge)
            if pos > 0 && counts[pos - 1] >= 3 {
                let error_idx = pos + k - 1;
                if error_idx >= template.sequence.len() || corrected_pos.contains(&error_idx) { continue; }

                let cur_fwd = hashes_fwd[pos];
                let cur_rc = hashes_rc[pos];
                let orig_bval = base_to_2bit(template.sequence[error_idx]);
                let mut best_count = counts[pos];
                let mut best_bval = orig_bval;

                for bval in 0u64..4 {
                    if bval == orig_bval { continue; }
                    let new_fwd = (cur_fwd & !3u64) | bval;
                    let new_rc = (cur_rc & !(3u64 << rc_shift)) | ((3 - bval) << rc_shift);
                    let new_canon = new_fwd.min(new_rc);
                    if let Some(&count) = graph.kmers.get(&new_canon) {
                        if count > best_count {
                            best_count = count;
                            best_bval = bval;
                        }
                    }
                }

                if best_bval != orig_bval && best_count >= 3 {
                    corrections.push((error_idx, idx_to_base(best_bval as usize)));
                    corrected_pos.insert(error_idx);
                }
            }

            // Check if right neighbor has high count → error at pos (leading edge)
            if pos + 1 < n_kmers && counts[pos + 1] >= 3 {
                let error_idx = pos;
                if corrected_pos.contains(&error_idx) { continue; }

                let cur_fwd = hashes_fwd[pos];
                let cur_rc = hashes_rc[pos];
                let orig_bval = base_to_2bit(template.sequence[error_idx]);
                let mut best_count = counts[pos];
                let mut best_bval = orig_bval;

                for bval in 0u64..4 {
                    if bval == orig_bval { continue; }
                    let new_fwd = (cur_fwd & !(3u64 << rc_shift)) | (bval << rc_shift);
                    let new_rc = (cur_rc & !3u64) | (3 - bval);
                    let new_canon = new_fwd.min(new_rc);
                    if let Some(&count) = graph.kmers.get(&new_canon) {
                        if count > best_count {
                            best_count = count;
                            best_bval = bval;
                        }
                    }
                }

                if best_bval != orig_bval && best_count >= 3 {
                    corrections.push((error_idx, idx_to_base(best_bval as usize)));
                    corrected_pos.insert(error_idx);
                }
            }
        }

        // Phase 3: apply corrections (mutable borrow of sequence)
        for (idx, base) in &corrections {
            template.sequence[*idx] = *base;
        }
        total_corrections += corrections.len();
    }

    total_corrections
}

/// Fast template extraction: filter graph to count≥2 (fits L3 cache), then
/// walk sequentially with used_seeds tracking for diversity.
/// Much faster than full extraction because 52M→3M k-mers fits in cache.
/// After extraction, extends tips through count=1 k-mers and error-corrects.
fn extract_templates_fast(graph: &mut TemplateGraph, min_length: usize) -> Vec<Template> {
    let before = graph.kmers.len();

    // Save full graph, replace with filtered for cache-friendly extraction
    let full_kmers = std::mem::take(&mut graph.kmers);
    graph.kmers = full_kmers.iter()
        .filter(|&(_, c)| *c >= 2)
        .map(|(&h, &c)| (h, c))
        .collect();

    let count2 = graph.kmers.len();
    eprintln!("  Fast filter: {} → {} k-mers (count ≥ 2)", before, count2);

    let mut templates = extract_templates(graph, min_length);

    // Only do tip extension + error correction if full graph fits in L3 cache (~30MB)
    // Each k-mer entry: ~20 bytes in FxHashMap → 50M entries ≈ 1GB
    // At >60M k-mers, the graph is too large for cache-friendly walking
    let max_graph_for_tips = 60_000_000;

    if before <= max_graph_for_tips {
        // Restore full graph for tip extension + error correction
        graph.kmers = full_kmers;
        let t = std::time::Instant::now();

        let (ext_count, ext_bases) = extend_tips(&mut templates, graph);
        eprintln!("  Tip extension: {} templates extended, {} bases added ({:.1}s)",
                 ext_count, ext_bases, t.elapsed().as_secs_f64());

        let corrections = error_correct_templates(&mut templates, graph);
        eprintln!("  Error correction: {} bases corrected", corrections);
    } else {
        eprintln!("  Skipping tip extension (graph too large: {}M k-mers > {}M threshold)",
                 before / 1_000_000, max_graph_for_tips / 1_000_000);
        // Error-correct using only count≥2 graph (still in memory, much smaller)
        let corrections = error_correct_templates(&mut templates, graph);
        eprintln!("  Error correction (count≥2 only): {} bases corrected", corrections);
        // Restore full kmers for potential future use (callers may need it)
        graph.kmers = full_kmers;
    }

    templates
}

// ── Read mapping with high mismatch tolerance ───────────────────────────────

fn build_template_index(
    templates: &[Template],
    k: usize,
) -> FxHashMap<u64, Vec<TemplateKmerHit>> {
    let total_bases: usize = templates.iter().map(|t| t.sequence.len()).sum();
    // Sample every `step` positions — for small k, index is very large
    let step = if total_bases > 100_000_000 { 8 }
               else if total_bases > 10_000_000 { 4 }
               else { 2 };

    let mut index: FxHashMap<u64, Vec<TemplateKmerHit>> = FxHashMap::default();
    index.reserve(total_bases / step);

    for (tid, template) in templates.iter().enumerate() {
        if template.sequence.len() < k { continue; }
        // Only index templates ≥ read length (short ones aren't useful)
        // Actually, index all — even shorter templates can match partial reads
        let mut pos = 0;
        while pos + k <= template.sequence.len() {
            let kmer = &template.sequence[pos..pos + k];
            if let Some(fwd) = kmer_to_hash(kmer) {
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                let is_forward = fwd <= rc;
                let bucket = index.entry(canon).or_default();
                // Cap entries per k-mer to avoid slow lookups on repetitive k-mers
                if bucket.len() < 100 {
                    bucket.push(TemplateKmerHit {
                        template_id: tid as u32,
                        position: pos as u32,
                        is_forward,
                    });
                }
            }
            pos += step;
        }
    }

    eprintln!("  Template index: {} entries (step={})",
             index.values().map(|v| v.len()).sum::<usize>(), step);

    index
}

fn map_reads_to_templates(
    sequences: &[Vec<u8>],
    templates: &[Template],
    index: &FxHashMap<u64, Vec<TemplateKmerHit>>,
    k: usize,
    mismatch_threshold: usize,
) -> Vec<Option<ReadMapping>> {
    sequences
        .par_iter()
        .map(|read| {
            if read.len() < k { return None; }

            // Use more anchor positions for small k (more hits per k-mer)
            let num_anchors = 8;
            let step = read.len().saturating_sub(k) / (num_anchors - 1).max(1);

            let mut best: Option<ReadMapping> = None;
            let mut best_dist = usize::MAX;

            for anchor_idx in 0..num_anchors {
                let anchor_pos = (anchor_idx * step).min(read.len().saturating_sub(k));
                if anchor_pos + k > read.len() { continue; }
                let kmer = &read[anchor_pos..anchor_pos + k];

                let fwd = match kmer_to_hash(kmer) {
                    Some(h) => h,
                    None => continue,
                };
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                let read_is_fwd = fwd <= rc;

                let hits = match index.get(&canon) {
                    Some(h) => h,
                    None => continue,
                };

                // Limit per-anchor work
                let max_hits = 50;
                for (hit_idx, hit) in hits.iter().enumerate() {
                    if hit_idx >= max_hits { break; }

                    let template = &templates[hit.template_id as usize];
                    let maps_fwd = read_is_fwd == hit.is_forward;

                    let read_start: i64 = if maps_fwd {
                        hit.position as i64 - anchor_pos as i64
                    } else {
                        hit.position as i64 - (read.len() as i64 - anchor_pos as i64 - k as i64)
                    };

                    if read_start < 0 { continue; }
                    let read_start = read_start as usize;
                    if read_start + read.len() > template.sequence.len() { continue; }

                    let region = &template.sequence[read_start..read_start + read.len()];

                    let dist = if maps_fwd {
                        hamming_distance_within(region, read, mismatch_threshold)
                    } else {
                        let rc_region = reverse_complement(region);
                        hamming_distance_within(&rc_region, read, mismatch_threshold)
                    };

                    if let Some(d) = dist {
                        if d < best_dist {
                            best_dist = d;
                            best = Some(ReadMapping {
                                template_id: hit.template_id as usize,
                                position: read_start as u32,
                                is_rc: !maps_fwd,
                            });
                            if d == 0 { return best; }
                        }
                    }
                }

                if best_dist == 0 { break; }
            }

            best
        })
        .collect()
}

// ── Syncmer-based read mapping ──────────────────────────────────────────────

fn build_template_index_syncmer(
    templates: &[Template],
    k: usize,
) -> FxHashMap<u64, Vec<TemplateKmerHit>> {
    // Closed syncmers with s close to k for high density
    // density ≈ 2/(k-s+1) for closed syncmers
    // s = k-3 → density ≈ 2/4 = 50%, matching stride step=2
    let s = k.saturating_sub(3).max(5);
    let t_end = k - s; // target positions: [0, k-s] for closed syncmers

    let mut index: FxHashMap<u64, Vec<TemplateKmerHit>> = FxHashMap::default();

    let mut total_indexed = 0usize;
    for (tid, template) in templates.iter().enumerate() {
        if template.sequence.len() < k { continue; }

        // Closed syncmers: smallest s-mer at position 0 or k-s
        let positions = syncmers::find_syncmers_pos(k, s, &[0, t_end], &template.sequence);

        for pos in positions {
            if pos + k > template.sequence.len() { continue; }
            let kmer = &template.sequence[pos..pos + k];
            if let Some(fwd) = kmer_to_hash(kmer) {
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                let is_forward = fwd <= rc;
                let bucket = index.entry(canon).or_default();
                if bucket.len() < 100 {
                    bucket.push(TemplateKmerHit {
                        template_id: tid as u32,
                        position: pos as u32,
                        is_forward,
                    });
                    total_indexed += 1;
                }
            }
        }
    }

    eprintln!("  Template index (closed syncmer): {} entries, {} unique k-mers (s={}, targets=[0,{}])",
             total_indexed, index.len(), s, t_end);

    index
}

fn map_reads_to_templates_syncmer(
    sequences: &[Vec<u8>],
    templates: &[Template],
    index: &FxHashMap<u64, Vec<TemplateKmerHit>>,
    k: usize,
    mismatch_threshold: usize,
) -> Vec<Option<ReadMapping>> {
    let s = k.saturating_sub(3).max(5);
    let t_end = k - s;

    sequences
        .par_iter()
        .map(|read| {
            if read.len() < k { return None; }

            // Get syncmer positions from the read (content-addressed selection)
            let positions = syncmers::find_syncmers_pos(k, s, &[0, t_end], read);

            if positions.is_empty() { return None; }

            let mut best: Option<ReadMapping> = None;
            let mut best_dist = usize::MAX;

            // Evenly space anchors across all syncmer positions
            let max_anchors = 16;
            let n_pos = positions.len();
            let step_pos = if n_pos > max_anchors { n_pos / max_anchors } else { 1 };
            for idx in (0..n_pos).step_by(step_pos).take(max_anchors) {
                let anchor_pos = positions[idx];
                if anchor_pos + k > read.len() { continue; }
                let kmer = &read[anchor_pos..anchor_pos + k];

                let fwd = match kmer_to_hash(kmer) {
                    Some(h) => h,
                    None => continue,
                };
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                let read_is_fwd = fwd <= rc;

                let hits = match index.get(&canon) {
                    Some(h) => h,
                    None => continue,
                };

                let max_hits = 50;
                for (hit_idx, hit) in hits.iter().enumerate() {
                    if hit_idx >= max_hits { break; }

                    let template = &templates[hit.template_id as usize];
                    let maps_fwd = read_is_fwd == hit.is_forward;

                    let read_start: i64 = if maps_fwd {
                        hit.position as i64 - anchor_pos as i64
                    } else {
                        hit.position as i64 - (read.len() as i64 - anchor_pos as i64 - k as i64)
                    };

                    if read_start < 0 { continue; }
                    let read_start = read_start as usize;
                    if read_start + read.len() > template.sequence.len() { continue; }

                    let region = &template.sequence[read_start..read_start + read.len()];

                    let dist = if maps_fwd {
                        hamming_distance_within(region, read, mismatch_threshold)
                    } else {
                        let rc_region = reverse_complement(region);
                        hamming_distance_within(&rc_region, read, mismatch_threshold)
                    };

                    if let Some(d) = dist {
                        if d < best_dist {
                            best_dist = d;
                            best = Some(ReadMapping {
                                template_id: hit.template_id as usize,
                                position: read_start as u32,
                                is_rc: !maps_fwd,
                            });
                            if d == 0 { return best; }
                        }
                    }
                }

                if best_dist == 0 { break; }
            }

            best
        })
        .collect()
}

// ── Hybrid encoding (SPRING-style) ──────────────────────────────────────────

/// Compress sequences using hybrid encoding (SPRING-style).
///
/// Like compress_sequences_template_subst, uses the de Bruijn graph pipeline
/// (steps 1-3) to build templates and map reads. But instead of substituting
/// matching bases with \0 in a single stream, it:
///
/// - **Singletons**: kept as raw ASCII in one large BSC stream (preserves BWT context)
/// - **Mapped reads**: REMOVED from BSC stream entirely, encoded as compact metadata:
///   - Flag stream: 1 byte per read (0=singleton, 1=mapped)
///   - Mapping metadata: (template_id, position, is_rc, mismatch_count) per mapped read
///   - Mismatch positions: delta-encoded per mapped read
///   - Mismatch bases: one base per mismatch
///   - Consensus: 2-bit packed templates
///
/// This avoids inflating the BSC stream with \0 bytes for mapped reads.
pub fn compress_sequences_template_hybrid(
    sequences: &[String],
    params: &TemplateParams,
) -> Result<Vec<u8>> {
    let seqs: Vec<Vec<u8>> = sequences.iter().map(|s| s.as_bytes().to_vec()).collect();
    if seqs.is_empty() { return Ok(Vec::new()); }

    let k = params.k;
    let mismatch_threshold = params.mismatch_threshold;

    eprintln!("Template-hybrid: {} reads, k={}, max_mismatches={}",
             seqs.len(), k, mismatch_threshold);

    // Steps 1-3: Build graph, extract templates, map reads
    eprintln!("  Step 1: Building graph (k={})...", k);
    let t = std::time::Instant::now();
    let mut graph = TemplateGraph::from_sequences(&seqs, k);

    if !params.fast {
        let total_kmers = graph.kmers.len();
        let singletons_count = graph.kmers.values().filter(|&&c| c == 1).count();
        let singleton_frac = singletons_count as f64 / total_kmers.max(1) as f64;
        if params.tip_trim_threshold > 0 && singleton_frac < 0.95 {
            graph.tip_trim(params.tip_trim_threshold);
        }
    }
    eprintln!("  Graph built: {:.1}s", t.elapsed().as_secs_f64());

    let min_length = seqs.iter().map(|s| s.len()).max().unwrap_or(150);
    eprintln!("  Step 2: Extracting templates (min_length={})...", min_length);
    let t = std::time::Instant::now();
    let templates = if params.fast {
        extract_templates_fast(&mut graph, min_length)
    } else {
        extract_templates(&graph, min_length)
    };
    drop(graph);
    eprintln!("  Templates extracted: {:.1}s", t.elapsed().as_secs_f64());

    let map_threshold = mismatch_threshold;
    eprintln!("  Step 3: Mapping reads (mismatch_threshold={}, syncmers={})...",
             map_threshold, params.use_syncmers);
    let t = std::time::Instant::now();
    let mappings = if params.use_syncmers {
        let index = build_template_index_syncmer(&templates, k);
        let m = map_reads_to_templates_syncmer(&seqs, &templates, &index, k, map_threshold);
        drop(index);
        m
    } else {
        let index = build_template_index(&templates, k);
        let m = map_reads_to_templates(&seqs, &templates, &index, k, map_threshold);
        drop(index);
        m
    };

    let mapped = mappings.iter().filter(|m| m.is_some()).count();
    eprintln!("  Mapping: {}/{} mapped ({:.1}%) ({:.1}s)",
             mapped, seqs.len(), mapped as f64 / seqs.len() as f64 * 100.0,
             t.elapsed().as_secs_f64());

    // Step 4: Hybrid encoding — separate singletons from mapped reads
    eprintln!("  Step 4: Hybrid encoding...");

    // Determine which templates are used and remap IDs
    let mut used_ids: FxHashSet<usize> = FxHashSet::default();
    for m in mappings.iter().flatten() {
        used_ids.insert(m.template_id);
    }
    let mut id_map: FxHashMap<usize, usize> = FxHashMap::default();
    let mut used_templates: Vec<&Template> = Vec::new();
    let mut sorted_ids: Vec<usize> = used_ids.into_iter().collect();
    sorted_ids.sort();
    for &uid in &sorted_ids {
        id_map.insert(uid, used_templates.len());
        used_templates.push(&templates[uid]);
    }

    // Stream 1: Flag stream — 1 byte per read (0=singleton, 1=mapped)
    let mut flag_stream: Vec<u8> = Vec::with_capacity(seqs.len());

    // Stream 2: Singleton bases — raw ASCII, one big BSC block
    let mut singleton_stream: Vec<u8> = Vec::new();

    // Stream 3: Mapping metadata — (template_id, position, is_rc) per mapped read
    let mut map_meta_stream: Vec<u8> = Vec::new();

    // Stream 4: Mismatch positions — delta-encoded per mapped read
    let mut mismatch_pos_stream: Vec<u8> = Vec::new();

    // Stream 5: Mismatch bases — one base per mismatch
    let mut mismatch_base_stream: Vec<u8> = Vec::new();

    let mut total_mismatches = 0u64;

    for (read_idx, mapping) in mappings.iter().enumerate() {
        let read = &seqs[read_idx];
        match mapping {
            Some(m) => {
                flag_stream.push(1);

                let mapped_id = id_map[&m.template_id];
                write_varint(&mut map_meta_stream, mapped_id as u64);
                write_varint(&mut map_meta_stream, m.position as u64);
                map_meta_stream.push(if m.is_rc { 1 } else { 0 });

                // Compute mismatches
                let template = &templates[m.template_id];
                let region = &template.sequence[m.position as usize..m.position as usize + read.len()];

                let compare_ref = if m.is_rc {
                    reverse_complement(region)
                } else {
                    region.to_vec()
                };

                let mut mm_count = 0u32;
                let mut prev_pos = 0u16;
                // Collect mismatches first to write count before positions
                let mut mm_positions: Vec<u16> = Vec::new();
                let mut mm_bases: Vec<u8> = Vec::new();

                for (i, (&r, &c)) in read.iter().zip(compare_ref.iter()).enumerate() {
                    if r != c {
                        mm_positions.push(i as u16);
                        mm_bases.push(r);
                        mm_count += 1;
                    }
                }

                total_mismatches += mm_count as u64;

                // Write mismatch count
                write_varint(&mut mismatch_pos_stream, mm_count as u64);

                // Delta-encoded positions
                for &pos in &mm_positions {
                    write_varint(&mut mismatch_pos_stream, (pos - prev_pos) as u64);
                    prev_pos = pos;
                }

                // Mismatch bases (raw, not delta-encoded — only 4 possible values)
                mismatch_base_stream.extend_from_slice(&mm_bases);
            }
            None => {
                flag_stream.push(0);
                singleton_stream.extend_from_slice(read);
            }
        }
    }

    let avg_mm = if mapped > 0 { total_mismatches as f64 / mapped as f64 } else { 0.0 };
    let num_singletons = seqs.len() - mapped;
    eprintln!("  Hybrid: {} mapped (avg {:.1} mm), {} singletons",
             mapped, avg_mm, num_singletons);
    eprintln!("  Raw streams: flags={}, singletons={}, map_meta={}, mm_pos={}, mm_bases={}",
             flag_stream.len(), singleton_stream.len(), map_meta_stream.len(),
             mismatch_pos_stream.len(), mismatch_base_stream.len());

    // Template consensus for decompression
    let consensus_seqs: Vec<Vec<u8>> = used_templates.iter().map(|t| t.sequence.clone()).collect();
    let consensus_bytes = pack_dna_2bit(&consensus_seqs);

    // Step 5: BSC-compress all streams
    eprintln!("  Step 5: BSC-compressing...");
    let t = std::time::Instant::now();

    // Use adaptive BSC for the big singleton stream (like baseline)
    let bsc_singletons = bsc::compress_parallel_adaptive(&singleton_stream)?;
    let bsc_flags = bsc::compress_parallel(&flag_stream)?;
    let bsc_map_meta = bsc::compress_parallel(&map_meta_stream)?;
    let bsc_mm_pos = bsc::compress_parallel(&mismatch_pos_stream)?;
    let bsc_mm_bases = bsc::compress_parallel(&mismatch_base_stream)?;
    let bsc_consensus = bsc::compress_parallel(&consensus_bytes)?;

    eprintln!("  BSC done: {:.1}s", t.elapsed().as_secs_f64());
    eprintln!("  BSC: singletons={}, flags={}, map_meta={}, mm_pos={}, mm_bases={}, consensus={}",
             bsc_singletons.len(), bsc_flags.len(), bsc_map_meta.len(),
             bsc_mm_pos.len(), bsc_mm_bases.len(), bsc_consensus.len());

    let total_compressed = bsc_singletons.len() + bsc_flags.len() + bsc_map_meta.len()
        + bsc_mm_pos.len() + bsc_mm_bases.len() + bsc_consensus.len();
    let total_bases: usize = seqs.iter().map(|s| s.len()).sum();
    eprintln!("  Total: {} bytes ({:.2}x)", total_compressed,
             total_bases as f64 / total_compressed as f64);

    // Serialize (THB1 = Template HyBrid v1)
    let mut output = Vec::new();
    output.extend_from_slice(b"THB1");
    write_u64(&mut output, seqs.len() as u64);
    write_u64(&mut output, bsc_singletons.len() as u64);
    write_u64(&mut output, bsc_flags.len() as u64);
    write_u64(&mut output, bsc_map_meta.len() as u64);
    write_u64(&mut output, bsc_mm_pos.len() as u64);
    write_u64(&mut output, bsc_mm_bases.len() as u64);
    write_u64(&mut output, bsc_consensus.len() as u64);
    output.extend_from_slice(&bsc_singletons);
    output.extend_from_slice(&bsc_flags);
    output.extend_from_slice(&bsc_map_meta);
    output.extend_from_slice(&bsc_mm_pos);
    output.extend_from_slice(&bsc_mm_bases);
    output.extend_from_slice(&bsc_consensus);

    let total = output.len();
    eprintln!("  Template-hybrid: {} → {} bytes ({:.2}x)",
             total_bases, total, total_bases as f64 / total as f64);

    Ok(output)
}

