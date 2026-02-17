/// Greedy contig-based sequence compression (patent-safe, order-preserving)
///
/// Approach inspired by SSAKE-style greedy assembly with minimizer-based dictionary:
/// 1. Build minimizer dictionary: guaranteed overlap detection for ≥ (w+k-1) bp overlaps
/// 2. Pre-filter isolated reads (no shared minimizers → instant singletons)
/// 3. Greedy contig extension: chain overlapping reads (Hamming ≤ threshold)
/// 4. Map ALL reads to contigs (parallel, preserving original order)
/// 5. Encode as (contig_id, position, is_rc, mismatches) + singletons
///
/// Patent-safe: reads are NOT reordered. Original order is preserved.
///
/// Performance optimizations (for 5M+ reads):
/// - minimizer-iter crate: O(n) minimizer computation via monotone queue
/// - Capped dictionary entries per minimizer (skip repetitive k-mers)
/// - Pre-filter isolated reads in parallel (~97% at 0.25x coverage)
/// - Progress reporting every 100K reads

use super::bsc;
use super::dna_utils::*;
use anyhow::Result;
use minimizer_iter::MinimizerBuilder;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

// ── Parameters ──────────────────────────────────────────────────────────────

const DICT_K: usize = 23;                 // k-mer size for dictionary
const MINIMIZER_W: u16 = 9;               // window size (must be odd for canonical mode)
                                           // overlap guarantee ≥ w + k - 1 = 31bp
const MISMATCH_THRESH_BUILD: usize = 4;   // Hamming threshold during contig building
const MISMATCH_THRESH_MAP: usize = 8;     // Hamming threshold during read mapping
const MAX_SEARCH: usize = 200;            // Max candidates per dictionary lookup
const MIN_OVERLAP: usize = 31;            // Minimum overlap for extension (= w + k - 1)
const MAX_ENTRIES_PER_MINIMIZER: usize = 50; // Cap dictionary entries (skip repetitive k-mers)

// ── Data structures ─────────────────────────────────────────────────────────

struct Contig {
    sequence: Vec<u8>,
}

#[derive(Clone)]
struct ReadMapping {
    contig_id: u32,
    position: u32,
    is_rc: bool,
}

/// Dictionary entry: (read_index, position_of_kmer_in_read, is_forward)
#[derive(Clone, Copy)]
struct DictEntry {
    read_idx: u32,
    kmer_pos: u16,     // position of minimizer k-mer within the read
    is_forward: bool,  // true if fwd_hash <= rc_hash for this k-mer in the read
}

// ── Minimizer computation (via minimizer-iter crate) ────────────────────────

/// Compute minimizer positions for a sequence using minimizer-iter's canonical mode.
/// Returns Vec<(position, canonical_hash, is_forward)> using our own hash function
/// for dictionary consistency.
///
/// Uses minimizer-iter's monotone queue algorithm: O(n) amortized instead of O(n*w).
/// Canonical mode ensures the same minimizer is selected regardless of strand.
fn compute_minimizers(seq: &[u8], k: usize, w: u16) -> Vec<(usize, u64, bool)> {
    if seq.len() < k + w as usize - 1 {
        return Vec::new();
    }

    // Use minimizer-iter to find canonical minimizer positions (O(n) via monotone queue)
    let positions: Vec<usize> = MinimizerBuilder::<u64>::new()
        .minimizer_size(k)
        .width(w)
        .canonical()
        .iter_pos(seq)
        .map(|(pos, _is_rc)| pos)
        .collect();

    // At each minimizer position, compute our canonical hash for dictionary key
    let mut result = Vec::with_capacity(positions.len());
    for pos in positions {
        if pos + k > seq.len() { continue; }
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            let canon = fwd.min(rc);
            let is_fwd = fwd <= rc;
            result.push((pos, canon, is_fwd));
        }
    }

    result
}

// ── Dictionary ──────────────────────────────────────────────────────────────

/// Build dictionary using minimizers from each read.
/// Caps entries per minimizer to skip repetitive k-mers.
fn build_dictionary(
    sequences: &[Vec<u8>],
    k: usize,
    w: u16,
) -> FxHashMap<u64, Vec<DictEntry>> {
    let mut dict: FxHashMap<u64, Vec<DictEntry>> = FxHashMap::default();
    let est_per_read = sequences.first()
        .map(|s| 2 * s.len().saturating_sub(k) / (w as usize + 1) + 1)
        .unwrap_or(0);
    dict.reserve(sequences.len() * est_per_read / 2); // many will be capped

    for (idx, seq) in sequences.iter().enumerate() {
        for (pos, canon, is_fwd) in compute_minimizers(seq, k, w) {
            let bucket = dict.entry(canon).or_default();
            if bucket.len() < MAX_ENTRIES_PER_MINIMIZER {
                bucket.push(DictEntry {
                    read_idx: idx as u32,
                    kmer_pos: pos as u16,
                    is_forward: is_fwd,
                });
            }
        }

        if idx % 500_000 == 0 && idx > 0 {
            eprintln!("    ... dictionary: {}/{} reads processed, {} unique minimizers",
                     idx, sequences.len(), dict.len());
        }
    }

    // Report capping stats
    let capped = dict.values().filter(|v| v.len() >= MAX_ENTRIES_PER_MINIMIZER).count();
    if capped > 0 {
        eprintln!("    Dictionary capping: {} minimizers hit cap of {} entries",
                 capped, MAX_ENTRIES_PER_MINIMIZER);
    }

    dict
}

// ── Pre-filtering ───────────────────────────────────────────────────────────

/// Pre-filter: identify reads that share at least one minimizer with another read.
/// Reads with no shared minimizers are guaranteed to be singletons (no contig building needed).
/// This is the key optimization for low-coverage data (~97% of reads are isolated at 0.25x).
fn identify_reads_with_neighbors(
    sequences: &[Vec<u8>],
    dict: &FxHashMap<u64, Vec<DictEntry>>,
    k: usize,
    w: u16,
) -> Vec<bool> {
    sequences
        .par_iter()
        .enumerate()
        .map(|(idx, seq)| {
            let minimizers = compute_minimizers(seq, k, w);
            for (_, canon, _) in &minimizers {
                if let Some(entries) = dict.get(canon) {
                    for entry in entries {
                        if entry.read_idx as usize != idx {
                            return true;
                        }
                    }
                }
            }
            false
        })
        .collect()
}

// ── Contig building ─────────────────────────────────────────────────────────

/// Evaluate a single candidate for overlap extension.
/// Returns (hamming_distance, oriented_read, read_start_in_contig) if valid.
fn evaluate_candidate(
    contig: &[u8],
    read: &[u8],
    probe_pos: usize,       // position of shared k-mer in contig
    entry_kmer_pos: usize,  // position of shared k-mer in read (original orientation)
    same_strand: bool,
    k: usize,
    threshold: usize,
    extend_right: bool,     // true = right extension, false = left extension
) -> Option<(usize, Vec<u8>, usize)> {
    let clen = contig.len();

    let read_kmer_pos = if same_strand {
        entry_kmer_pos
    } else {
        read.len() - entry_kmer_pos - k
    };

    if extend_right {
        // Right extension: read must extend past contig end
        if read_kmer_pos > probe_pos {
            return None;
        }
        let read_start = probe_pos - read_kmer_pos;
        let overlap = clen - read_start;
        if overlap < MIN_OVERLAP || overlap > read.len() {
            return None;
        }
        let extension = read.len() - overlap;
        if extension == 0 { return None; }

        let oriented_read = if same_strand {
            read.to_vec()
        } else {
            reverse_complement(read)
        };

        let contig_overlap = &contig[read_start..clen];
        let read_overlap = &oriented_read[..overlap];
        let d = hamming_distance_within(contig_overlap, read_overlap, threshold)?;
        Some((d, oriented_read, read_start))
    } else {
        // Left extension: read must extend before contig start
        if read_kmer_pos <= probe_pos {
            return None;
        }
        let extension = read_kmer_pos - probe_pos;
        if extension == 0 { return None; }
        let overlap = read.len() - extension;
        if overlap < MIN_OVERLAP || overlap > clen {
            return None;
        }

        let oriented_read = if same_strand {
            read.to_vec()
        } else {
            reverse_complement(read)
        };

        let contig_overlap = &contig[..overlap];
        let read_overlap = &oriented_read[extension..extension + overlap];
        let d = hamming_distance_within(contig_overlap, read_overlap, threshold)?;
        Some((d, oriented_read, extension))
    }
}

/// Try to find a read that overlaps the right end of the contig using minimizers.
fn find_right_extension(
    contig: &[u8],
    k: usize,
    w: u16,
    dict: &FxHashMap<u64, Vec<DictEntry>>,
    sequences: &[Vec<u8>],
    used: &[bool],
    threshold: usize,
) -> Option<(u32, Vec<u8>, usize)> {
    let clen = contig.len();
    if clen < k { return None; }

    // Compute minimizers from the tail of the contig (last read_len bases)
    let probe_start = clen.saturating_sub(sequences[0].len());
    let tail = &contig[probe_start..];
    let minimizers = compute_minimizers(tail, k, w);

    let mut best: Option<(u32, Vec<u8>, usize)> = None;
    let mut best_dist = usize::MAX;

    for (local_pos, canon, contig_is_fwd) in &minimizers {
        let probe_pos = probe_start + local_pos; // position in full contig

        let candidates = match dict.get(canon) {
            Some(c) => c,
            None => continue,
        };

        let mut checked = 0;
        for entry in candidates.iter().rev() {
            if used[entry.read_idx as usize] { continue; }
            checked += 1;
            if checked > MAX_SEARCH { break; }

            let read = &sequences[entry.read_idx as usize];
            let same_strand = *contig_is_fwd == entry.is_forward;

            if let Some((d, oriented_read, read_start)) = evaluate_candidate(
                contig, read, probe_pos, entry.kmer_pos as usize,
                same_strand, k, threshold, true,
            ) {
                if d < best_dist {
                    best_dist = d;
                    best = Some((entry.read_idx, oriented_read, read_start));
                    if d == 0 { return best; }
                }
            }
        }

        if best_dist == 0 { break; }
    }

    best
}

/// Try to find a read that overlaps the left end of the contig using minimizers.
fn find_left_extension(
    contig: &[u8],
    k: usize,
    w: u16,
    dict: &FxHashMap<u64, Vec<DictEntry>>,
    sequences: &[Vec<u8>],
    used: &[bool],
    threshold: usize,
) -> Option<(u32, Vec<u8>, usize)> {
    if contig.len() < k { return None; }

    // Compute minimizers from the head of the contig (first read_len bases)
    let probe_end = contig.len().min(sequences[0].len());
    let head = &contig[..probe_end];
    let minimizers = compute_minimizers(head, k, w);

    let mut best: Option<(u32, Vec<u8>, usize)> = None;
    let mut best_dist = usize::MAX;

    for (probe_pos, canon, contig_is_fwd) in &minimizers {
        let candidates = match dict.get(canon) {
            Some(c) => c,
            None => continue,
        };

        let mut checked = 0;
        for entry in candidates.iter().rev() {
            if used[entry.read_idx as usize] { continue; }
            checked += 1;
            if checked > MAX_SEARCH { break; }

            let read = &sequences[entry.read_idx as usize];
            let same_strand = *contig_is_fwd == entry.is_forward;

            if let Some((d, oriented_read, extension)) = evaluate_candidate(
                contig, read, *probe_pos, entry.kmer_pos as usize,
                same_strand, k, threshold, false,
            ) {
                if d < best_dist {
                    best_dist = d;
                    best = Some((entry.read_idx, oriented_read, extension));
                    if d == 0 { return best; }
                }
            }
        }

        if best_dist == 0 { break; }
    }

    best
}

/// Build contigs greedily — only processes reads that have shared minimizers.
fn build_contigs(
    sequences: &[Vec<u8>],
    dict: &FxHashMap<u64, Vec<DictEntry>>,
    has_neighbor: &[bool],
    k: usize,
    w: u16,
) -> Vec<Contig> {
    let n = sequences.len();
    let mut used = vec![false; n];
    let mut contigs: Vec<Contig> = Vec::new();
    let read_len = sequences[0].len();
    let mut reads_in_contigs = 0usize;
    let mut singletons_skipped = 0usize;

    for seed_idx in 0..n {
        if used[seed_idx] { continue; }
        used[seed_idx] = true;

        // Fast path: isolated reads become instant singletons
        if !has_neighbor[seed_idx] {
            singletons_skipped += 1;
            contigs.push(Contig { sequence: sequences[seed_idx].clone() });
            continue;
        }

        reads_in_contigs += 1;
        let mut contig_seq = sequences[seed_idx].clone();

        // === Extend RIGHT ===
        for _ in 0..10000 {
            match find_right_extension(&contig_seq, k, w, dict, sequences, &used, MISMATCH_THRESH_BUILD) {
                Some((read_idx, oriented_read, read_start)) => {
                    let overlap = contig_seq.len() - read_start;
                    let new_bases = &oriented_read[overlap..];
                    contig_seq.extend_from_slice(new_bases);
                    used[read_idx as usize] = true;
                    reads_in_contigs += 1;
                }
                None => break,
            }
        }

        // === Extend LEFT ===
        for _ in 0..10000 {
            match find_left_extension(&contig_seq, k, w, dict, sequences, &used, MISMATCH_THRESH_BUILD) {
                Some((read_idx, oriented_read, extension_len)) => {
                    let new_bases = &oriented_read[..extension_len];
                    let mut new_contig = new_bases.to_vec();
                    new_contig.extend_from_slice(&contig_seq);
                    contig_seq = new_contig;
                    used[read_idx as usize] = true;
                    reads_in_contigs += 1;
                }
                None => break,
            }
        }

        contigs.push(Contig { sequence: contig_seq });

        // Progress reporting
        let total_processed = singletons_skipped + reads_in_contigs;
        if total_processed % 100_000 == 0 && total_processed > 0 {
            eprintln!("    ... contigs: {}/{} processed, {} contigs, {} chained, {} singletons",
                     seed_idx, n, contigs.len(), reads_in_contigs, singletons_skipped);
        }
    }

    let total_bases: usize = contigs.iter().map(|c| c.sequence.len()).sum();
    let max_len = contigs.iter().map(|c| c.sequence.len()).max().unwrap_or(0);
    let long_contigs = contigs.iter().filter(|c| c.sequence.len() > read_len * 2).count();
    eprintln!("  Built {} contigs, {} reads chained ({:.1}%), {} instant singletons",
             contigs.len(), reads_in_contigs,
             reads_in_contigs as f64 / n as f64 * 100.0,
             singletons_skipped);
    eprintln!("  Total={}bp, max={}bp, {} long contigs (>{}bp)",
             total_bases, max_len, long_contigs, read_len * 2);

    contigs
}

// ── Read mapping ────────────────────────────────────────────────────────────

/// Build index of contigs for fast read mapping
fn build_contig_index(
    contigs: &[Contig],
    k: usize,
) -> FxHashMap<u64, Vec<(u32, u32, bool)>> {
    let mut index: FxHashMap<u64, Vec<(u32, u32, bool)>> = FxHashMap::default();

    let total_bases: usize = contigs.iter().map(|c| c.sequence.len()).sum();
    let step = if total_bases > 50_000_000 { 4 } else { 1 };

    for (cid, contig) in contigs.iter().enumerate() {
        if contig.sequence.len() < k { continue; }
        let mut pos = 0;
        while pos + k <= contig.sequence.len() {
            let kmer = &contig.sequence[pos..pos + k];
            if let Some(fwd) = kmer_to_hash(kmer) {
                let rc = reverse_complement_hash(fwd, k);
                let canon = fwd.min(rc);
                let is_forward = fwd <= rc;
                index.entry(canon).or_default().push((cid as u32, pos as u32, is_forward));
            }
            pos += step;
        }
    }

    index
}

/// Map reads to contigs in parallel
fn map_reads_to_contigs(
    sequences: &[Vec<u8>],
    contigs: &[Contig],
    index: &FxHashMap<u64, Vec<(u32, u32, bool)>>,
    k: usize,
) -> Vec<Option<ReadMapping>> {
    sequences
        .par_iter()
        .map(|read| {
            if read.len() < k { return None; }

            let anchors = [
                0,
                read.len() / 4,
                read.len() / 2,
                3 * read.len() / 4,
                read.len().saturating_sub(k),
            ];

            let mut best: Option<ReadMapping> = None;
            let mut best_dist = usize::MAX;

            for &anchor_pos in &anchors {
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

                for &(cid, cpos, hit_is_fwd) in hits {
                    let contig = &contigs[cid as usize];
                    let maps_fwd = read_is_fwd == hit_is_fwd;

                    let read_start: i64 = if maps_fwd {
                        cpos as i64 - anchor_pos as i64
                    } else {
                        cpos as i64 - (read.len() as i64 - anchor_pos as i64 - k as i64)
                    };

                    if read_start < 0 { continue; }
                    let read_start = read_start as usize;
                    if read_start + read.len() > contig.sequence.len() { continue; }

                    let region = &contig.sequence[read_start..read_start + read.len()];

                    let dist = if maps_fwd {
                        hamming_distance_within(region, read, MISMATCH_THRESH_MAP)
                    } else {
                        let rc_region = reverse_complement(region);
                        hamming_distance_within(&rc_region, read, MISMATCH_THRESH_MAP)
                    };

                    if let Some(d) = dist {
                        if d < best_dist {
                            best_dist = d;
                            best = Some(ReadMapping {
                                contig_id: cid,
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

// ── Encoding ────────────────────────────────────────────────────────────────

fn encode_reads(
    sequences: &[Vec<u8>],
    contigs: &[Contig],
    mappings: &[Option<ReadMapping>],
) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    let contig_seqs: Vec<Vec<u8>> = contigs.iter().map(|c| c.sequence.clone()).collect();
    let consensus_bytes = pack_dna_2bit(&contig_seqs);

    let mut meta_bytes: Vec<u8> = Vec::new();
    let mut mismatch_bytes: Vec<u8> = Vec::new();
    write_varint(&mut meta_bytes, sequences.len() as u64);

    let mut singleton_seqs: Vec<Vec<u8>> = Vec::new();

    for (read_idx, mapping) in mappings.iter().enumerate() {
        match mapping {
            Some(m) => {
                meta_bytes.push(1);
                write_varint(&mut meta_bytes, m.contig_id as u64);
                write_varint(&mut meta_bytes, m.position as u64);
                meta_bytes.push(if m.is_rc { 1 } else { 0 });
                write_varint(&mut meta_bytes, sequences[read_idx].len() as u64);

                let read = &sequences[read_idx];
                let contig = &contigs[m.contig_id as usize];
                let region = &contig.sequence[m.position as usize..m.position as usize + read.len()];

                let (compare_read, compare_ref) = if m.is_rc {
                    (reverse_complement(read), region.to_vec())
                } else {
                    (read.clone(), region.to_vec())
                };

                let mut mismatches: Vec<(u16, u8)> = Vec::new();
                for (i, (&r, &c)) in compare_read.iter().zip(compare_ref.iter()).enumerate() {
                    if r != c {
                        mismatches.push((i as u16, r));
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
                meta_bytes.push(0);
                singleton_seqs.push(sequences[read_idx].clone());
            }
        }
    }

    let singleton_bytes = pack_dna_2bit(&singleton_seqs);

    let mapped = mappings.iter().filter(|m| m.is_some()).count();
    eprintln!("  Encoded: {} mapped, {} singletons",
             mapped, singleton_seqs.len());

    (consensus_bytes, meta_bytes, mismatch_bytes, singleton_bytes)
}

// ── Public API ──────────────────────────────────────────────────────────────

pub fn compress_sequences_greedy(sequences: &[String]) -> Result<Vec<u8>> {
    let seqs: Vec<Vec<u8>> = sequences.iter().map(|s| s.as_bytes().to_vec()).collect();

    if seqs.is_empty() {
        return Ok(Vec::new());
    }

    let k = DICT_K;
    let w = MINIMIZER_W;
    eprintln!("Greedy contig compression: {} reads, k={}, w={} (overlap guarantee ≥ {}bp)",
             seqs.len(), k, w, w as usize + k - 1);

    // Step 1: Build minimizer dictionary (with capping)
    eprintln!("  Step 1: Building minimizer dictionary (cap={} entries/minimizer)...",
             MAX_ENTRIES_PER_MINIMIZER);
    let t = std::time::Instant::now();
    let dict = build_dictionary(&seqs, k, w);
    eprintln!("  Dictionary: {} unique minimizers, {} entries ({:.1}s)",
             dict.len(),
             dict.values().map(|v| v.len()).sum::<usize>(),
             t.elapsed().as_secs_f64());

    // Step 2: Pre-filter isolated reads
    eprintln!("  Step 2: Identifying reads with shared minimizers...");
    let t = std::time::Instant::now();
    let has_neighbor = identify_reads_with_neighbors(&seqs, &dict, k, w);
    let num_with_neighbors = has_neighbor.iter().filter(|&&b| b).count();
    eprintln!("  {} / {} reads have shared minimizers ({:.1}%), {} isolated ({:.1}s)",
             num_with_neighbors, seqs.len(),
             num_with_neighbors as f64 / seqs.len() as f64 * 100.0,
             seqs.len() - num_with_neighbors,
             t.elapsed().as_secs_f64());

    // Step 3: Build contigs greedily (only from reads with neighbors)
    eprintln!("  Step 3: Building contigs...");
    let t = std::time::Instant::now();
    let contigs = build_contigs(&seqs, &dict, &has_neighbor, k, w);
    eprintln!("  Contig building: {:.1}s", t.elapsed().as_secs_f64());

    // Step 4: Build contig index and map ALL reads
    eprintln!("  Step 4: Mapping reads to contigs...");
    let t = std::time::Instant::now();
    let contig_index = build_contig_index(&contigs, k);
    let mappings = map_reads_to_contigs(&seqs, &contigs, &contig_index, k);

    let mapped = mappings.iter().filter(|m| m.is_some()).count();
    eprintln!("  Mapping: {}/{} reads mapped ({:.1}%) ({:.1}s)",
             mapped, seqs.len(), mapped as f64 / seqs.len() as f64 * 100.0,
             t.elapsed().as_secs_f64());

    // Step 5: Encode
    eprintln!("  Step 5: Encoding...");
    let (consensus_raw, meta_raw, mismatch_raw, singleton_raw) =
        encode_reads(&seqs, &contigs, &mappings);

    // Step 6: BSC-compress
    eprintln!("  Step 6: BSC-compressing...");
    let bsc_consensus = bsc::compress_parallel(&consensus_raw)?;
    let bsc_meta = bsc::compress_parallel(&meta_raw)?;
    let bsc_mismatches = bsc::compress_parallel(&mismatch_raw)?;
    let bsc_singletons = bsc::compress_parallel(&singleton_raw)?;

    eprintln!("  Raw streams: consensus={}, meta={}, mismatches={}, singletons={}",
             consensus_raw.len(), meta_raw.len(), mismatch_raw.len(), singleton_raw.len());
    eprintln!("  BSC streams: consensus={}, meta={}, mismatches={}, singletons={}",
             bsc_consensus.len(), bsc_meta.len(), bsc_mismatches.len(), bsc_singletons.len());

    // Serialize (DBG1 format — compatible with debruijn.rs decompressor)
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
    let original: usize = sequences.iter().map(|s| s.len()).sum();
    eprintln!("  Greedy contig: {} → {} bytes ({:.2}x)",
             original, total, original as f64 / total as f64);

    Ok(output)
}
