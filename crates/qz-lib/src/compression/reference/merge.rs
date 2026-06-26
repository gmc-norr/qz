//! Reference-aware v5 archive merge (paired reference `archive_type 2` + single-end
//! reference `archive_type 4`).
//!
//! Used by the NUMA-compress track for reference mode (shard input → compress each
//! shard against the SAME reference on a pinned process → merge), and exposed via the
//! reference `qz merge` path. Unlike single-end/paired merge, reference archives carry
//! 4 whole-archive coverage globals (PackedBacking / NBitmap / IntervalMap /
//! ReferenceMeta) that are NOT verbatim-copyable: each shard's IntervalMap covers only
//! the reference intervals ITS reads touched, so the merged archive needs the globals
//! re-derived over the UNION of all shards' covered intervals.
//!
//! Why that re-derivation is sound and bounded:
//! * **Per-chunk columns are shard-invariant.** Positions are stored ABSOLUTE
//!   `(ref_id, ref_pos)` (see `positions.rs`); edits/strands/read-len/mapped-flags and
//!   the inline FallbackPool are per-read functions of `(reference, reads)`, never of
//!   the coverage map. So every per-(chunk,mate) frame relocates VERBATIM (its
//!   per-block CRC is position-independent), exactly like single-end/paired merge.
//! * **The 3 coverage globals are reference-direct.** PackedBacking/NBitmap are the
//!   raw normalized reference bases over the covered intervals (`stream_packed_backing`
//!   pulls `reference[ref_id][ref_pos]`, NOT a read-consensus), and IntervalMap is the
//!   covered-interval set itself. Re-deriving them from `(reference, union-imap)` is
//!   deterministic. ReferenceMeta depends only on the reference, so it is shard-invariant.
//! * **Union correctness.** Each read lies in exactly one shard; that shard's coverage
//!   marks the read's whole contiguous span, so the span is contained in one of that
//!   shard's intervals. Unioning + coalescing the shard intervals therefore keeps every
//!   read's span inside a single merged interval — decode's `span_in_single_interval`
//!   check holds, and `ref_to_base` resolves against the merged map (which decode
//!   recomputes from the merged IntervalMap, so the flat base indices stay consistent).
//!
//! Inputs are fully validated up front with the SAME contracts decode runs
//! (`validate_reference_directory{,_single}` + the footer span/CRC checks in
//! `read_v5_footer`), and the per-mate `ChunkDecodedSizes` table is rebuilt so the
//! merged archive stays NUMA direct-write decodable.

use std::io::{Read, Seek, SeekFrom, Write};
use std::path::PathBuf;

use anyhow::{Result, bail};
use strobealign::io::fasta::RefSequence;

use super::backing::{Interval, IntervalMap, read_intervalmap_blockstream};
use super::encode::stream_globals_parallel;
use super::format::{validate_reference_directory, validate_reference_directory_single};
use super::refmeta::reference_digest;
use crate::compression::chunk_directory::{
    ChunkDirEntry, ChunkDirectory, GLOBAL_SENTINEL, StreamRole, encode_decoded_sizes,
    read_v5_footer, validate_directory_spans, write_footer_and_locator,
};
use crate::compression::compress_impl::write_role_blocks;
use crate::compression::decompress_impl::read_fixed_header;
use crate::compression::read_chunk_layout;

/// One validated reference shard ready to relocate.
struct RShard {
    path: PathBuf,
    dir: ChunkDirectory,
    /// This shard's covered IntervalMap (parsed from its IntervalMap global).
    imap: IntervalMap,
    /// `decoded_sizes_per_mate[mate][chunk]`, length == n_mates (1 single / 2 paired).
    per_mate: Vec<Vec<u64>>,
}

/// Merge N reference archives of the SAME `archive_type` (2 or 4) into one, streamed
/// to `out`, re-deriving the coverage globals from `refs` (the normalized reference
/// the shards were compressed against). `refs` MUST be loaded via
/// `mapping::load_reference_normalized` so its `(ref_id-order, normalized-bytes)` view
/// matches the per-shard encoders' — otherwise the re-derived globals would not line
/// up with the per-chunk absolute positions.
pub(crate) fn merge_reference_core(
    inputs: &[PathBuf],
    refs: &[RefSequence],
    out: &mut dyn Write,
) -> Result<()> {
    if inputs.is_empty() {
        bail!("qz merge: need at least one input archive");
    }

    // The archive_type is fixed by the first input; require the rest to agree. 2 =
    // paired reference (2 mates), 4 = single-end reference (1 mate).
    let (hdr0, header_end0) = read_fixed_header(&inputs[0])?;
    let archive_type = hdr0.archive_type;
    let n_mates: u8 = match archive_type {
        2 => 2,
        4 => 1,
        other => bail!("qz merge: not a reference archive (archive_type {other})"),
    };
    let validate: fn(&ChunkDirectory) -> Result<()> = if archive_type == 2 {
        validate_reference_directory
    } else {
        validate_reference_directory_single
    };

    // The coverage globals are re-derived from `refs`, and decode reconstructs every
    // read's bases from that re-derived backing (it never re-reads a FASTA). So if `refs`
    // is not the reference the shards were compressed against, the merged archive is
    // silently corrupt AND passes verify (verify only checks internal CRC/structure).
    // Each shard stores an FNV digest of its reference in the ReferenceMeta global; gate
    // on it. Requiring every shard to equal `reference_digest(refs)` also forces the
    // shards to agree with each other (no mixed-reference merge).
    let want_digest = reference_digest(refs);

    // Pass 1: validate + load every shard (front-header identity, mode directory
    // contract, IntervalMap global, reference digest, per-mate decoded-size table).
    let mut header_bytes: Option<Vec<u8>> = None;
    let mut shards: Vec<RShard> = Vec::with_capacity(inputs.len());

    for path in inputs {
        let (hdr, header_end) = read_fixed_header(path)?;
        if hdr.archive_type != archive_type {
            bail!(
                "qz merge: cannot mix archive types ({} has type {}, expected {})",
                path.display(),
                hdr.archive_type,
                archive_type
            );
        }
        // Front-header byte-identity gate (same reference config across shards).
        let mut f = std::fs::File::open(path)?;
        let mut hbuf = vec![0u8; header_end as usize];
        f.read_exact(&mut hbuf)?;
        match &header_bytes {
            None => header_bytes = Some(hbuf),
            Some(h0) => {
                if header_end != header_end0 || &hbuf != h0 {
                    bail!(
                        "qz merge: incompatible archive configs (front headers differ) — cannot merge {}",
                        path.display()
                    );
                }
            }
        }

        let dir = read_v5_footer(path, header_end)?;
        // Full reference structural validation (the same §9 contract decode runs).
        validate(&dir)?;

        // Parse this shard's IntervalMap global (same role/read the decoder uses).
        let e_imap = dir
            .entries
            .iter()
            .find(|e| e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::IntervalMap)
            .ok_or_else(|| {
                anyhow::anyhow!("qz merge: {} missing IntervalMap global", path.display())
            })?;
        let imap_bytes = read_frame(&mut f, e_imap.offset, e_imap.length)?;
        let imap = read_intervalmap_blockstream(&imap_bytes)?;

        // Gate the shard's reference against `refs` via the ReferenceMeta digest. The
        // entry payload begins `[meta_version:1B][digest:8B LE]` (mirrors decode's read).
        let e_meta = dir
            .entries
            .iter()
            .find(|e| e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::ReferenceMeta)
            .ok_or_else(|| {
                anyhow::anyhow!("qz merge: {} missing ReferenceMeta global", path.display())
            })?;
        let meta_bytes = read_frame(&mut f, e_meta.offset, e_meta.length)?;
        if meta_bytes.len() < 9 {
            bail!("qz merge: {} ReferenceMeta shorter than 9-byte prefix", path.display());
        }
        let shard_digest = u64::from_le_bytes(meta_bytes[1..9].try_into().unwrap());
        if shard_digest != want_digest {
            bail!(
                "qz merge: {} was compressed against a different reference \
                 (stored digest {shard_digest:#018x} != supplied-FASTA digest {want_digest:#018x}) \
                 — merging would silently corrupt the output",
                path.display()
            );
        }

        // Validated, CRC-checked per-mate decoded sizes (rejects a corrupt table).
        let layout = read_chunk_layout(path)?;
        if layout.decoded_sizes_per_mate.len() != n_mates as usize
            || layout
                .decoded_sizes_per_mate
                .iter()
                .any(|m| m.len() != dir.num_chunks as usize)
        {
            bail!(
                "qz merge: {} has no ChunkDecodedSizes(n_mates={n_mates}) table; every reference input must have one",
                path.display()
            );
        }

        shards.push(RShard { path: path.clone(), dir, imap, per_mate: layout.decoded_sizes_per_mate });
    }
    let header_bytes = header_bytes.expect("inputs non-empty");

    // Union the shards' covered intervals into one coalesced merged IntervalMap.
    let merged_imap = union_intervalmaps(shards.iter().map(|s| &s.imap))?;

    // Emit the shared front header.
    out.write_all(&header_bytes)?;
    let mut pos: u64 = header_end0;

    // Pass 2: copy per-chunk frames verbatim in FOOTER order (globals excluded — they
    // are re-derived below); renumber chunk_index; mate/role/codec/record_count carried
    // verbatim. Concatenate the per-mate size tables row-major (`sizes[chunk*n_mates + mate]`).
    let mut merged: Vec<ChunkDirEntry> = Vec::new();
    let mut merged_sizes: Vec<u64> = Vec::new();
    let mut chunk_base: u32 = 0;
    let mut num_reads: u64 = 0;
    let overflow = || anyhow::anyhow!("qz merge: archive metadata overflow (forged input?)");

    for shard in &shards {
        let mut f = std::fs::File::open(&shard.path)?;
        for e in shard.dir.entries.iter().filter(|e| e.chunk_index != GLOBAL_SENTINEL) {
            let frame = read_frame(&mut f, e.offset, e.length)?;
            out.write_all(&frame)?;
            merged.push(ChunkDirEntry {
                chunk_index: chunk_base.checked_add(e.chunk_index).ok_or_else(overflow)?,
                role: e.role,
                mate: e.mate,
                codec: e.codec,
                offset: pos,
                length: e.length,
                record_count: e.record_count,
            });
            pos = pos.checked_add(e.length).ok_or_else(overflow)?;
        }
        // Append this shard's chunks row-major: [m0_c0, m1_c0, m0_c1, m1_c1, …].
        for c in 0..shard.dir.num_chunks as usize {
            for mate in 0..n_mates as usize {
                merged_sizes.push(shard.per_mate[mate][c]);
            }
        }
        chunk_base = chunk_base.checked_add(shard.dir.num_chunks).ok_or_else(overflow)?;
        num_reads = num_reads.checked_add(shard.dir.num_reads).ok_or_else(overflow)?;
    }

    // Re-derive the 4 whole-archive globals over the UNION imap (PackedBacking → NBitmap
    // → IntervalMap → ReferenceMeta). The PARALLEL variant: the chr20-sized backing
    // re-derivation is the dominant merge tax, so it BSC-compresses its blocks across
    // all (now-idle) cores rather than single-threaded. Decodes identically to the
    // encoder's globals (the block partition is internal to each global).
    stream_globals_parallel(out, &merged_imap, refs, &mut pos, &mut merged)?;

    // Rebuild the single merged ChunkDecodedSizes(n_mates) global (every input had one).
    debug_assert_eq!(merged_sizes.len(), chunk_base as usize * n_mates as usize);
    let payload = encode_decoded_sizes(n_mates, chunk_base, &merged_sizes);
    write_role_blocks(
        &[(chunk_base, payload)],
        StreamRole::ChunkDecodedSizes,
        0, // mate (the global is always mate 0)
        0, // codec (raw metadata)
        GLOBAL_SENTINEL,
        out,
        &mut pos,
        &mut merged,
    )?;

    // Self-check the rebuilt directory with decode's own validators, then write tail.
    let dir = ChunkDirectory { num_reads, num_chunks: chunk_base, entries: merged };
    validate_directory_spans(&dir.entries, header_end0, pos)?;
    validate(&dir)?;
    write_footer_and_locator(out, &dir)?;
    out.flush()?;
    Ok(())
}

/// Coalesce N covered IntervalMaps into one (union of covered positions), merging
/// overlapping AND adjacent intervals so the result is a strictly-ordered,
/// non-overlapping set `IntervalMap::from_intervals` accepts. Gap-merge filler is
/// intentionally NOT re-applied across shards: every read's span is already fully
/// covered within one shard's interval, so no real coverage is lost (and the merged
/// archive is marginally smaller than re-filling cross-shard gaps would make it).
fn union_intervalmaps<'a>(maps: impl Iterator<Item = &'a IntervalMap>) -> Result<IntervalMap> {
    let mut all: Vec<Interval> = Vec::new();
    for m in maps {
        all.extend_from_slice(m.intervals());
    }
    all.sort_by_key(|iv| (iv.ref_id, iv.ref_start));

    let mut merged: Vec<Interval> = Vec::with_capacity(all.len());
    for iv in all {
        if let Some(last) = merged.last_mut()
            && last.ref_id == iv.ref_id
        {
            // Inputs come from validated IntervalMaps, so ref_start+length never overflows.
            let last_end = last.ref_start + last.length;
            if iv.ref_start <= last_end {
                let iv_end = iv.ref_start + iv.length;
                if iv_end > last_end {
                    last.length = iv_end - last.ref_start;
                }
                continue;
            }
        }
        merged.push(iv);
    }
    IntervalMap::from_intervals(merged)
}

/// Read a `[offset, offset+length)` block frame from `f` into a fresh buffer.
fn read_frame(f: &mut std::fs::File, offset: u64, length: u64) -> Result<Vec<u8>> {
    f.seek(SeekFrom::Start(offset))?;
    let mut buf = vec![0u8; length as usize];
    f.read_exact(&mut buf)?;
    Ok(buf)
}
