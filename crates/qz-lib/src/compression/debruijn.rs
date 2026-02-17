/// De Bruijn graph-based sequence compression
///
/// Two-pass approach:
/// 1. Build canonical de Bruijn graph from reads → extract unitigs (maximal unbranched paths)
/// 2. Map each read to a unitig + position + orientation + mismatches
///
/// RC-awareness: canonical k-mers (min of k-mer and its RC) mean forward and
/// reverse-complement reads map to the same graph edges automatically.

use super::dna_utils::*;
use super::bsc;
use anyhow::Result;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

const MISMATCH_THRESHOLD: usize = 8;

/// Auto-select k-mer size based on dataset characteristics.
///
/// Key insight: for the de Bruijn graph to produce useful unitigs, k-mers must
/// appear multiple times. With N reads of length L, total k-mers ≈ N*(L-k+1).
/// We want enough repetition that after tip-trimming, substantial graph structure
/// remains.
///
/// Heuristic: pick k so that average k-mer multiplicity >= 3.
/// multiplicity ≈ total_kmers / unique_kmers
/// For DNA, unique k-mers ≈ min(4^k, genome_size). Since we don't know genome size,
/// we use the empirical rule: start from k=31 and reduce until multiplicity is adequate.
fn auto_select_k(num_reads: usize, avg_read_len: usize) -> usize {
    let total_bases = num_reads as f64 * avg_read_len as f64;

    // Simple tier-based selection calibrated to real data:
    // - Very shallow (<5M bases): k=15 — lots of repetition, short unitigs but high mapping
    // - Shallow (<50M bases):     k=19 — balanced
    // - Moderate (<500M bases):   k=23 — good specificity
    // - Deep (<2B bases):         k=27 — high specificity
    // - Very deep (>=2B bases):   k=31 — maximum specificity
    let k = if total_bases < 5_000_000.0 {
        15
    } else if total_bases < 50_000_000.0 {
        19
    } else if total_bases < 500_000_000.0 {
        23
    } else if total_bases < 2_000_000_000.0 {
        27
    } else {
        31
    };

    // k must be odd (canonical k-mer hashing works best with odd k to avoid palindromes)
    k
}

// ── Data structures ──────────────────────────────────────────────────────────

struct DeBruijnGraph {
    /// Canonical k-mer hash → occurrence count
    kmers: FxHashMap<u64, u16>,
    k: usize,
}

struct Unitig {
    sequence: Vec<u8>,
}

#[derive(Clone)]
struct ReadMapping {
    unitig_id: usize,
    position: u32,
    is_rc: bool,
}

// ── Pass 1: Build graph + extract unitigs ────────────────────────────────────

/// Compute canonical k-mer hash: min(hash, rc_hash)
fn canonical_hash(seq: &[u8], k: usize) -> Option<u64> {
    let hash = kmer_to_hash(seq)?;
    let rc = reverse_complement_hash(hash, k);
    Some(hash.min(rc))
}

impl DeBruijnGraph {
    /// Build de Bruijn graph from sequences by counting canonical k-mers (parallel)
    fn from_sequences(sequences: &[Vec<u8>], k: usize) -> Self {
        // Parallel k-mer counting: each thread builds a local map, then merge
        let chunk_size = (sequences.len() / rayon::current_num_threads().max(1)).max(1000);

        let local_maps: Vec<FxHashMap<u64, u16>> = sequences
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut local: FxHashMap<u64, u16> = FxHashMap::default();
                local.reserve(chunk.len() * 128); // ~128 k-mers per 150bp read
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

        // Merge local maps into a single global map
        let estimated_size = local_maps.iter().map(|m| m.len()).max().unwrap_or(0);
        let mut kmers: FxHashMap<u64, u16> = FxHashMap::default();
        kmers.reserve(estimated_size);
        for local in local_maps {
            for (hash, count) in local {
                let entry = kmers.entry(hash).or_insert(0);
                *entry = entry.saturating_add(count);
            }
        }

        eprintln!("  De Bruijn graph: {} unique canonical {}-mers from {} reads",
                 kmers.len(), k, sequences.len());

        DeBruijnGraph { kmers, k }
    }

    /// Remove k-mers with count <= threshold (tip trimming for error removal)
    fn tip_trim(&mut self, min_count: u16) {
        let before = self.kmers.len();
        self.kmers.retain(|_, count| *count > min_count);
        eprintln!("  Tip trimming (count <= {}): {} → {} k-mers",
                 min_count, before, self.kmers.len());
    }

    /// Check if a canonical k-mer exists in the graph
    fn contains(&self, canon_hash: u64) -> bool {
        self.kmers.contains_key(&canon_hash)
    }

    /// Get the forward sequence of a k-mer from its hash
    fn hash_to_sequence(hash: u64, k: usize) -> Vec<u8> {
        let mut seq = vec![0u8; k];
        let mut h = hash;
        for i in (0..k).rev() {
            seq[i] = idx_to_base(h as usize & 3);
            h >>= 2;
        }
        seq
    }

    /// Find successor k-mers: try extending the (k-1)-suffix with each base.
    /// Returns list of (canonical_hash, appended_base) for existing successors.
    fn successors(&self, kmer_seq: &[u8]) -> Vec<(u64, u8)> {
        let k = self.k;
        let suffix = &kmer_seq[1..]; // (k-1)-suffix
        let mut result = Vec::new();
        for &base in &[b'A', b'C', b'G', b'T'] {
            let mut candidate = Vec::with_capacity(k);
            candidate.extend_from_slice(suffix);
            candidate.push(base);
            if let Some(canon) = canonical_hash(&candidate, k) {
                if self.contains(canon) {
                    result.push((canon, base));
                }
            }
        }
        result
    }

    /// Find predecessor k-mers: try prepending each base to the (k-1)-prefix.
    fn predecessors(&self, kmer_seq: &[u8]) -> Vec<(u64, u8)> {
        let k = self.k;
        let prefix = &kmer_seq[..k - 1]; // (k-1)-prefix
        let mut result = Vec::new();
        for &base in &[b'A', b'C', b'G', b'T'] {
            let mut candidate = Vec::with_capacity(k);
            candidate.push(base);
            candidate.extend_from_slice(prefix);
            if let Some(canon) = canonical_hash(&candidate, k) {
                if self.contains(canon) {
                    result.push((canon, base));
                }
            }
        }
        result
    }
}

/// Extract unitigs (maximal unbranched paths) from the de Bruijn graph.
fn extract_unitigs(graph: &DeBruijnGraph) -> Vec<Unitig> {
    let k = graph.k;
    let mut visited: FxHashSet<u64> = FxHashSet::default();
    visited.reserve(graph.kmers.len());
    let mut unitigs: Vec<Unitig> = Vec::new();

    // Collect all canonical hashes for iteration
    let all_hashes: Vec<u64> = graph.kmers.keys().copied().collect();

    for &start_hash in &all_hashes {
        if visited.contains(&start_hash) { continue; }

        // Recover the forward sequence for this canonical hash
        // We need to pick a consistent orientation. Use the forward one (hash == canonical).
        let start_seq = DeBruijnGraph::hash_to_sequence(start_hash, k);
        // Verify this is actually the canonical orientation
        let start_seq = if kmer_to_hash(&start_seq) == Some(start_hash) {
            start_seq
        } else {
            reverse_complement(&DeBruijnGraph::hash_to_sequence(start_hash, k))
        };

        // Walk RIGHT from this k-mer
        let mut right_bases: Vec<u8> = Vec::new();
        let mut current_seq = start_seq.clone();
        loop {
            let succs = graph.successors(&current_seq);
            if succs.len() != 1 { break; }

            let (succ_canon, base) = succs[0];
            if succ_canon == start_hash { break; } // cycle
            if visited.contains(&succ_canon) { break; }

            // Check that the successor also has exactly 1 predecessor
            // (to ensure we're on an unbranched path)
            let mut next_seq = Vec::with_capacity(k);
            next_seq.extend_from_slice(&current_seq[1..]);
            next_seq.push(base);
            let preds = graph.predecessors(&next_seq);
            if preds.len() != 1 { break; }

            visited.insert(succ_canon);
            right_bases.push(base);
            current_seq = next_seq;
        }

        // Walk LEFT from the start k-mer
        let mut left_bases: Vec<u8> = Vec::new();
        current_seq = start_seq.clone();
        loop {
            let preds = graph.predecessors(&current_seq);
            if preds.len() != 1 { break; }

            let (pred_canon, base) = preds[0];
            if pred_canon == start_hash { break; } // cycle
            if visited.contains(&pred_canon) { break; }

            // Check successor uniqueness from predecessor
            let mut prev_seq = Vec::with_capacity(k);
            prev_seq.push(base);
            prev_seq.extend_from_slice(&current_seq[..k - 1]);
            let succs = graph.successors(&prev_seq);
            if succs.len() != 1 { break; }

            visited.insert(pred_canon);
            left_bases.push(base);
            current_seq = prev_seq;
        }

        visited.insert(start_hash);

        // Assemble unitig: left_bases (reversed) + start_kmer + right_bases
        left_bases.reverse();
        let mut unitig_seq = Vec::with_capacity(left_bases.len() + k + right_bases.len());
        unitig_seq.extend_from_slice(&left_bases);
        unitig_seq.extend_from_slice(&start_seq);
        unitig_seq.extend_from_slice(&right_bases);

        unitigs.push(Unitig { sequence: unitig_seq });
    }

    // Sort by length descending (longest unitigs first)
    unitigs.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));

    let total_bases: usize = unitigs.iter().map(|u| u.sequence.len()).sum();
    let max_len = unitigs.first().map(|u| u.sequence.len()).unwrap_or(0);
    eprintln!("  Extracted {} unitigs ({} total bases, longest={})",
             unitigs.len(), total_bases, max_len);

    unitigs
}

// ── Pass 2: Map reads to unitigs ─────────────────────────────────────────────

/// Index entry: position in a unitig where a k-mer occurs
struct UnitigKmerHit {
    unitig_id: usize,
    position: u32,       // position of the k-mer in the unitig
    is_forward: bool,    // true if the k-mer is in forward orientation at this position
}

/// Build an index of all k-mers in all unitigs for fast read mapping.
/// Samples every `step` positions to reduce index size for large datasets.
fn build_unitig_index(unitigs: &[Unitig], k: usize) -> FxHashMap<u64, Vec<UnitigKmerHit>> {
    let total_unitig_bases: usize = unitigs.iter().map(|u| u.sequence.len()).sum();
    // Sample every `step` positions to keep index manageable
    let step = if total_unitig_bases > 100_000_000 { 4 } else { 1 };

    let mut index: FxHashMap<u64, Vec<UnitigKmerHit>> = FxHashMap::default();
    index.reserve(total_unitig_bases / step.max(1));

    for (uid, unitig) in unitigs.iter().enumerate() {
        if unitig.sequence.len() < k { continue; }
        let mut pos = 0;
        while pos + k <= unitig.sequence.len() {
            let kmer = &unitig.sequence[pos..pos + k];
            if let Some(fwd_hash) = kmer_to_hash(kmer) {
                let rc_hash = reverse_complement_hash(fwd_hash, k);
                let canon = fwd_hash.min(rc_hash);
                let is_forward = fwd_hash <= rc_hash;

                index.entry(canon).or_default().push(UnitigKmerHit {
                    unitig_id: uid,
                    position: pos as u32,
                    is_forward,
                });
            }
            pos += step;
        }
    }

    index
}

/// Map reads to unitigs using the k-mer index (parallel)
fn map_reads_to_unitigs(
    sequences: &[Vec<u8>],
    unitigs: &[Unitig],
    index: &FxHashMap<u64, Vec<UnitigKmerHit>>,
    k: usize,
) -> Vec<Option<ReadMapping>> {
    let mappings: Vec<Option<ReadMapping>> = sequences
        .par_iter()
        .map(|read| {
            if read.len() < k { return None; }

            let anchor_positions = [
                0,
                read.len() / 3,
                2 * read.len() / 3,
                read.len().saturating_sub(k),
            ];

            let mut best_mapping: Option<ReadMapping> = None;
            let mut best_mismatches = usize::MAX;

            for &anchor_pos in &anchor_positions {
                if anchor_pos + k > read.len() { continue; }
                let anchor_kmer = &read[anchor_pos..anchor_pos + k];

                let fwd_hash = match kmer_to_hash(anchor_kmer) {
                    Some(h) => h,
                    None => continue,
                };
                let rc_hash = reverse_complement_hash(fwd_hash, k);
                let canon = fwd_hash.min(rc_hash);
                let read_is_forward = fwd_hash <= rc_hash;

                let hits = match index.get(&canon) {
                    Some(h) => h,
                    None => continue,
                };

                for hit in hits {
                    let unitig = &unitigs[hit.unitig_id];
                    let read_maps_forward = read_is_forward == hit.is_forward;

                    let read_start: i64 = if read_maps_forward {
                        hit.position as i64 - anchor_pos as i64
                    } else {
                        hit.position as i64 - (read.len() as i64 - anchor_pos as i64 - k as i64)
                    };

                    if read_start < 0 { continue; }
                    let read_start = read_start as usize;
                    if read_start + read.len() > unitig.sequence.len() { continue; }

                    let unitig_region = &unitig.sequence[read_start..read_start + read.len()];

                    let (compare_seq, is_rc) = if read_maps_forward {
                        (read.as_slice(), false)
                    } else {
                        (read.as_slice(), true)
                    };

                    let dist = if is_rc {
                        let rc_region = reverse_complement(unitig_region);
                        match hamming_distance_within(&rc_region, compare_seq, MISMATCH_THRESHOLD) {
                            Some(d) => d,
                            None => continue,
                        }
                    } else {
                        match hamming_distance_within(unitig_region, compare_seq, MISMATCH_THRESHOLD) {
                            Some(d) => d,
                            None => continue,
                        }
                    };

                    if dist < best_mismatches {
                        best_mismatches = dist;
                        best_mapping = Some(ReadMapping {
                            unitig_id: hit.unitig_id,
                            position: read_start as u32,
                            is_rc,
                        });
                        if dist == 0 { break; }
                    }
                }

                if best_mismatches == 0 { break; }
            }

            best_mapping
        })
        .collect();

    let mapped = mappings.iter().filter(|m| m.is_some()).count();
    let singleton = mappings.len() - mapped;
    eprintln!("  Read mapping: {} / {} mapped ({:.1}%), {} singletons",
             mapped, sequences.len(),
             mapped as f64 / sequences.len() as f64 * 100.0,
             singleton);

    mappings
}

// ── Encoding / Serialization ─────────────────────────────────────────────────

/// Encode reads against their unitig assignments into the 4-stream format.
/// Returns the raw (pre-BSC) bytes for the sequence archive.
fn encode_reads(
    sequences: &[Vec<u8>],
    unitigs: &[Unitig],
    mappings: &[Option<ReadMapping>],
) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    // Determine which unitigs are actually used
    let mut used_unitig_ids: FxHashSet<usize> = FxHashSet::default();
    for m in mappings.iter().flatten() {
        used_unitig_ids.insert(m.unitig_id);
    }

    // Remap unitig IDs to compact range
    let mut unitig_id_map: FxHashMap<usize, usize> = FxHashMap::default();
    let mut used_unitigs: Vec<&Unitig> = Vec::new();
    for &uid in &used_unitig_ids {
        unitig_id_map.insert(uid, used_unitigs.len());
        used_unitigs.push(&unitigs[uid]);
    }

    // Stream 0: consensus sequences (2-bit packed)
    let consensus_seqs: Vec<Vec<u8>> = used_unitigs.iter().map(|u| u.sequence.clone()).collect();
    let consensus_bytes = pack_dna_2bit(&consensus_seqs);

    // Stream 1: metadata + Stream 2: mismatches
    let mut meta_bytes: Vec<u8> = Vec::new();
    let mut mismatch_bytes: Vec<u8> = Vec::new();

    write_varint(&mut meta_bytes, sequences.len() as u64);

    // Stream 3: singletons
    let mut singleton_seqs: Vec<Vec<u8>> = Vec::new();

    for (read_idx, mapping) in mappings.iter().enumerate() {
        match mapping {
            Some(m) => {
                let compact_id = unitig_id_map[&m.unitig_id];
                meta_bytes.push(1); // has_contig
                write_varint(&mut meta_bytes, compact_id as u64);
                write_varint(&mut meta_bytes, m.position as u64);
                meta_bytes.push(if m.is_rc { 1 } else { 0 });
                write_varint(&mut meta_bytes, sequences[read_idx].len() as u64);

                // Compute mismatches
                let read = &sequences[read_idx];
                let unitig = &unitigs[m.unitig_id];
                let region = &unitig.sequence[m.position as usize..m.position as usize + read.len()];

                let (compare_read, compare_ref) = if m.is_rc {
                    (reverse_complement(read), region.to_vec())
                } else {
                    (read.clone(), region.to_vec())
                };

                let mut mismatches: Vec<(u16, u8)> = Vec::new();
                for (i, (&r, &c)) in compare_read.iter().zip(compare_ref.iter()).enumerate() {
                    if r != c {
                        mismatches.push((i as u16, r)); // store the READ's base
                    }
                }

                write_varint(&mut mismatch_bytes, mismatches.len() as u64);
                let mut prev_pos: u16 = 0;
                for (pos, base) in &mismatches {
                    write_varint(&mut mismatch_bytes, (*pos - prev_pos) as u64);
                    mismatch_bytes.push(*base);
                    prev_pos = *pos;
                }
            }
            None => {
                meta_bytes.push(0); // singleton
                singleton_seqs.push(sequences[read_idx].clone());
            }
        }
    }

    let singleton_bytes = pack_dna_2bit(&singleton_seqs);

    eprintln!("  Encoding: {} unitigs used, {} singletons packed",
             used_unitigs.len(), singleton_seqs.len());

    (consensus_bytes, meta_bytes, mismatch_bytes, singleton_bytes)
}

/// Compress sequences using de Bruijn graph unitig extraction.
///
/// `k_override`: k-mer size. 0 = auto-select based on dataset size.
pub fn compress_sequences_debruijn(sequences: &[String], k_override: usize) -> Result<Vec<u8>> {
    let seqs: Vec<Vec<u8>> = sequences.iter().map(|s| s.as_bytes().to_vec()).collect();

    if seqs.is_empty() {
        return Ok(Vec::new());
    }

    let max_readlen = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    let avg_readlen = seqs.iter().map(|s| s.len()).sum::<usize>() / seqs.len();

    // Select k: use override if provided, otherwise auto-select
    let k = if k_override > 0 {
        // Clamp to valid range and ensure odd
        let mut k = k_override.clamp(9, 31);
        if k % 2 == 0 { k -= 1; } // make odd
        k
    } else {
        auto_select_k(seqs.len(), avg_readlen)
    };

    eprintln!("De Bruijn compression: {} reads, max_readlen={}, k={} ({})",
             seqs.len(), max_readlen, k,
             if k_override > 0 { "user-specified" } else { "auto-selected" });

    // If reads are shorter than k, fall back to packing them all as singletons
    if max_readlen < k {
        eprintln!("  Reads shorter than k={}, encoding all as singletons", k);
        let packed = pack_dna_2bit(&seqs);
        let compressed = bsc::compress_parallel(&packed)?;
        let mut output = Vec::new();
        output.extend_from_slice(b"DBG1");
        write_u64(&mut output, seqs.len() as u64);
        write_u64(&mut output, 0); // 0 consensus bytes
        write_u64(&mut output, 0); // 0 meta bytes
        write_u64(&mut output, 0); // 0 mismatch bytes
        write_u64(&mut output, compressed.len() as u64);
        output.extend_from_slice(&compressed);
        return Ok(output);
    }

    // Pass 1: Build graph + extract unitigs
    eprintln!("  Pass 1: Building de Bruijn graph...");
    let mut graph = DeBruijnGraph::from_sequences(&seqs, k);

    // Adaptive tip trimming: if most k-mers are singletons, trimming kills
    // the graph. Check singleton fraction and adjust threshold.
    let total_kmers = graph.kmers.len();
    let singletons = graph.kmers.values().filter(|&&c| c == 1).count();
    let singleton_frac = singletons as f64 / total_kmers.max(1) as f64;

    let trim_threshold = if singleton_frac > 0.90 {
        // >90% singletons: skip trimming entirely (shallow data, most k-mers are real)
        eprintln!("  Singleton fraction: {:.1}% — skipping tip trimming (shallow coverage)",
                 singleton_frac * 100.0);
        0
    } else if singleton_frac > 0.70 {
        // 70-90%: gentle trim
        eprintln!("  Singleton fraction: {:.1}% — gentle tip trimming (threshold=1)",
                 singleton_frac * 100.0);
        1
    } else {
        // <70%: deep coverage, aggressive trim removes errors effectively
        eprintln!("  Singleton fraction: {:.1}% — standard tip trimming (threshold=1)",
                 singleton_frac * 100.0);
        1
    };

    if trim_threshold > 0 {
        graph.tip_trim(trim_threshold);
    }

    let unitigs = extract_unitigs(&graph);

    // Pass 2: Map reads to unitigs
    eprintln!("  Pass 2: Mapping reads to unitigs...");
    let index = build_unitig_index(&unitigs, k);
    let mappings = map_reads_to_unitigs(&seqs, &unitigs, &index, k);

    // Encode into 4 streams
    eprintln!("  Encoding streams...");
    let (consensus_raw, meta_raw, mismatch_raw, singleton_raw) =
        encode_reads(&seqs, &unitigs, &mappings);

    // BSC-compress each stream
    eprintln!("  BSC-compressing 4 streams...");
    let bsc_consensus = bsc::compress_parallel(&consensus_raw)?;
    let bsc_meta = bsc::compress_parallel(&meta_raw)?;
    let bsc_mismatches = bsc::compress_parallel(&mismatch_raw)?;
    let bsc_singletons = bsc::compress_parallel(&singleton_raw)?;

    eprintln!("  Raw: consensus={}, meta={}, mismatches={}, singletons={}",
             consensus_raw.len(), meta_raw.len(), mismatch_raw.len(), singleton_raw.len());
    eprintln!("  BSC: consensus={}, meta={}, mismatches={}, singletons={}",
             bsc_consensus.len(), bsc_meta.len(), bsc_mismatches.len(), bsc_singletons.len());

    // Serialize: magic + num_reads + 4 stream lengths + 4 streams
    let mut output = Vec::new();
    output.extend_from_slice(b"DBG1");
    write_u64(&mut output, seqs.len() as u64);
    write_u64(&mut output, bsc_consensus.len() as u64);
    write_u64(&mut output, bsc_meta.len() as u64);
    write_u64(&mut output, bsc_mismatches.len() as u64);
    write_u64(&mut output, bsc_singletons.len() as u64);
    output.extend_from_slice(&bsc_consensus);
    output.extend_from_slice(&bsc_meta);
    output.extend_from_slice(&bsc_mismatches);
    output.extend_from_slice(&bsc_singletons);

    let total = output.len();
    let original = sequences.iter().map(|s| s.len()).sum::<usize>();
    eprintln!("  De Bruijn sequence compression: {} → {} bytes ({:.2}x)",
             original, total, original as f64 / total as f64);

    Ok(output)
}

// ── Decompression ────────────────────────────────────────────────────────────

/// Decompress de Bruijn-encoded sequences
pub fn decompress_sequences_debruijn(data: &[u8], num_reads: usize) -> Result<Vec<String>> {
    // Read and verify magic
    if data.len() < 4 || &data[0..4] != b"DBG1" {
        anyhow::bail!("Invalid de Bruijn magic bytes");
    }
    let mut offset = 4;

    let stored_num_reads = read_u64(data, &mut offset)
        .ok_or_else(|| anyhow::anyhow!("Failed to read num_reads"))? as usize;
    if stored_num_reads != num_reads {
        anyhow::bail!("Read count mismatch: expected {}, got {}", num_reads, stored_num_reads);
    }

    let bsc_consensus_len = read_u64(data, &mut offset)
        .ok_or_else(|| anyhow::anyhow!("Failed to read consensus length"))? as usize;
    let bsc_meta_len = read_u64(data, &mut offset)
        .ok_or_else(|| anyhow::anyhow!("Failed to read meta length"))? as usize;
    let bsc_mismatches_len = read_u64(data, &mut offset)
        .ok_or_else(|| anyhow::anyhow!("Failed to read mismatches length"))? as usize;
    let bsc_singletons_len = read_u64(data, &mut offset)
        .ok_or_else(|| anyhow::anyhow!("Failed to read singletons length"))? as usize;

    // Handle short-read fallback (all singletons)
    if bsc_consensus_len == 0 && bsc_meta_len == 0 && bsc_mismatches_len == 0 {
        let singleton_compressed = &data[offset..offset + bsc_singletons_len];
        let singleton_raw = bsc::decompress_parallel(singleton_compressed)?;
        let (seqs, _) = unpack_dna_2bit(&singleton_raw, 0)?;
        return Ok(seqs.into_iter().map(|s| String::from_utf8_lossy(&s).to_string()).collect());
    }

    // BSC-decompress all 4 streams
    let consensus_compressed = &data[offset..offset + bsc_consensus_len];
    offset += bsc_consensus_len;
    let meta_compressed = &data[offset..offset + bsc_meta_len];
    offset += bsc_meta_len;
    let mismatches_compressed = &data[offset..offset + bsc_mismatches_len];
    offset += bsc_mismatches_len;
    let singletons_compressed = &data[offset..offset + bsc_singletons_len];

    let consensus_raw = bsc::decompress_parallel(consensus_compressed)?;
    let meta_raw = bsc::decompress_parallel(meta_compressed)?;
    let mismatch_raw = bsc::decompress_parallel(mismatches_compressed)?;
    let singleton_raw = bsc::decompress_parallel(singletons_compressed)?;

    // Unpack consensus unitigs
    let (unitig_seqs, _) = unpack_dna_2bit(&consensus_raw, 0)?;

    // Unpack singletons
    let (singleton_seqs, _) = unpack_dna_2bit(&singleton_raw, 0)?;

    // Parse metadata and reconstruct reads
    let mut meta_offset = 0;
    let meta_num_reads = read_varint(&meta_raw, &mut meta_offset)
        .ok_or_else(|| anyhow::anyhow!("Failed to read metadata num_reads"))? as usize;

    if meta_num_reads != num_reads {
        anyhow::bail!("Metadata read count mismatch");
    }

    let mut mismatch_offset = 0;
    let mut singleton_idx = 0;
    let mut sequences: Vec<String> = Vec::with_capacity(num_reads);

    for _ in 0..num_reads {
        if meta_offset >= meta_raw.len() {
            anyhow::bail!("Truncated metadata");
        }

        let has_contig = meta_raw[meta_offset];
        meta_offset += 1;

        if has_contig == 1 {
            let contig_id = read_varint(&meta_raw, &mut meta_offset)
                .ok_or_else(|| anyhow::anyhow!("Failed to read contig_id"))? as usize;
            let position = read_varint(&meta_raw, &mut meta_offset)
                .ok_or_else(|| anyhow::anyhow!("Failed to read position"))? as usize;
            let is_rc = meta_raw[meta_offset] != 0;
            meta_offset += 1;
            let read_len = read_varint(&meta_raw, &mut meta_offset)
                .ok_or_else(|| anyhow::anyhow!("Failed to read read_length"))? as usize;

            // Extract region from unitig
            if contig_id >= unitig_seqs.len() {
                anyhow::bail!("Invalid contig_id {}", contig_id);
            }
            let unitig = &unitig_seqs[contig_id];
            if position + read_len > unitig.len() {
                anyhow::bail!("Read extends beyond unitig");
            }
            let mut read = unitig[position..position + read_len].to_vec();

            // Apply mismatches
            let num_mismatches = read_varint(&mismatch_raw, &mut mismatch_offset)
                .ok_or_else(|| anyhow::anyhow!("Failed to read mismatch count"))? as usize;

            let mut prev_pos: u16 = 0;
            for _ in 0..num_mismatches {
                let delta = read_varint(&mismatch_raw, &mut mismatch_offset)
                    .ok_or_else(|| anyhow::anyhow!("Failed to read mismatch delta"))? as u16;
                let base = mismatch_raw[mismatch_offset];
                mismatch_offset += 1;

                let pos = prev_pos + delta;
                if (pos as usize) < read.len() {
                    read[pos as usize] = base;
                }
                prev_pos = pos;
            }

            // Reverse complement if needed
            if is_rc {
                read = reverse_complement(&read);
            }

            sequences.push(String::from_utf8_lossy(&read).to_string());
        } else {
            // Singleton
            if singleton_idx >= singleton_seqs.len() {
                anyhow::bail!("Ran out of singletons");
            }
            let seq = &singleton_seqs[singleton_idx];
            singleton_idx += 1;
            sequences.push(String::from_utf8_lossy(seq).to_string());
        }
    }

    Ok(sequences)
}
