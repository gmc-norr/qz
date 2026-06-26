//! v5 archive merge: type-dispatched merge engine for single-end, paired, and
//! reference. Single-end/paired merge through `merge_archives` (verbatim frame copy);
//! reference goes through `merge_reference_archives_to_path` (it needs the FASTA to
//! re-derive the coverage globals — see `reference::merge`).
//!
//! Used by the NUMA-compress track (shard input → compress each shard on a pinned
//! process → merge) and exposed as `qz merge`. Approach: directory-driven
//! block-granular copy — the encoder is never re-run; each shard's block frames
//! are relocated verbatim (per-block CRCs are position-independent), the directory
//! is rebuilt with shifted offsets + renumbered `chunk_index`, and the per-shard
//! `ChunkDecodedSizes` tables are merged into one so the result stays NUMA
//! direct-write decodable. Seek-free: streams to `out`.
//!
//! Inputs are fully validated up front with the SAME contract decode runs
//! (`validate_single_end_directory`), and per-chunk decoded sizes come from
//! `read_chunk_layout` (which CRC-checks the source table) — so a malformed or
//! corrupt input is rejected, never laundered into a fresh, valid-looking output.

use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

use super::chunk_directory::{
    ChunkDirEntry, ChunkDirectory, GLOBAL_SENTINEL, StreamRole, encode_decoded_sizes,
    read_v5_footer, validate_directory_spans, write_footer_and_locator,
};
use super::compress_impl::write_role_blocks;
use super::decompress_impl::{
    decode_block_bounds, read_fixed_header, validate_single_end_directory,
    validate_single_end_directory_structure,
};
use super::read_chunk_layout;

struct Shard {
    path: PathBuf,
    dir: ChunkDirectory,
    /// Validated, CRC-checked per-chunk decoded output byte lengths (from
    /// `read_chunk_layout`). Length == `dir.num_chunks`.
    decoded_sizes: Vec<u64>,
}

/// Read the `archive_type` byte (header byte 3) of every input; require they all
/// agree. Returns the common type, or errors on disagreement / unreadable header.
fn common_archive_type(inputs: &[PathBuf]) -> Result<u8> {
    if inputs.is_empty() {
        bail!("qz merge: need at least one input archive");
    }
    let mut ty: Option<u8> = None;
    for p in inputs {
        let (hdr, _) = read_fixed_header(p)?;
        match ty {
            None => ty = Some(hdr.archive_type),
            Some(t) if t != hdr.archive_type => bail!(
                "qz merge: cannot mix archive types ({} has type {}, expected {})",
                p.display(),
                hdr.archive_type,
                t
            ),
            Some(_) => {}
        }
    }
    Ok(ty.expect("inputs non-empty"))
}

/// Merge N v5 archives of the same `archive_type` into one, streamed to `out`.
/// Dispatches on the type: single-end (0) and paired (1) merge here; reference
/// (2/4) is handled separately by `merge_reference_archives_to_path` (it needs the
/// FASTA); cluster (3) cannot merge.
pub fn merge_archives(inputs: &[PathBuf], out: &mut dyn Write) -> Result<()> {
    match common_archive_type(inputs)? {
        0 => merge_single_core(inputs, out),
        1 => merge_paired_core(inputs, out),
        2 | 4 => bail!(
            "qz merge: reference archives need `qz merge --reference <fasta>` (the per-shard coverage globals are re-derived from the reference, not concatenated)"
        ),
        3 => bail!(
            "qz merge: cluster archives cannot be merged (order is dropped per-archive; a merged set would not be reconstructible)"
        ),
        other => bail!("qz merge: unknown archive_type {other}"),
    }
}

/// File-output variant (atomic O_EXCL temp + rename), mirroring `compress()`.
pub fn merge_archives_to_path(inputs: &[PathBuf], output: &Path, force: bool) -> Result<()> {
    merge_to_path_with(inputs, output, force, merge_archives)
}

/// Merge N reference archives (`archive_type` 2 = paired reference / 4 = single-end
/// reference) into the file at `output`, re-deriving the coverage globals from
/// `reference` (the SAME FASTA the shards were compressed against). Atomic O_EXCL
/// temp + rename, like `compress()`.
///
/// Reference merge is split from the generic `merge_archives` (which has no reference
/// arm) because it needs the reference bytes to rebuild PackedBacking/NBitmap/
/// IntervalMap/ReferenceMeta over the union of the shards' covered intervals — see
/// [`super::reference::merge_reference_core`]. The reference is loaded ONCE here and
/// shared with the streaming core.
pub fn merge_reference_archives_to_path(
    inputs: &[PathBuf],
    reference: &Path,
    output: &Path,
    force: bool,
) -> Result<()> {
    let refs = super::reference::load_reference_normalized(reference)?;
    merge_to_path_with(inputs, output, force, move |inputs, out| {
        super::reference::merge_reference_core(inputs, &refs, out)
    })
}

/// Paired (`archive_type 1`) v5 merge: stitch N complete paired archives into one,
/// streamed to `out`. Mirrors `merge_single_core` — directory-driven verbatim frame copy,
/// rebuilt directory with renumbered `chunk_index`, `mate` carried verbatim, and ONE
/// rebuilt `ChunkDecodedSizes` global so NUMA direct-write decode still applies — with two
/// paired specifics: (1) the table is `n_mates=2` (per-mate R1/R2 decoded sizes, row-major
/// `sizes[chunk*2 + mate]`); (2) structural validation uses `validate_paired_directory`
/// (the 6-segment, mate-tagged contract). Like single-end, it does NOT re-CRC every block:
/// `read_v5_footer` already validates per-entry spans, and decode CRC-checks each block
/// frame at read time — so the eager whole-archive CRC walk is unnecessary.
///
/// Every input must carry the `ChunkDecodedSizes(n_mates=2)` table (the streaming paired
/// encoder always emits it, and the NUMA-compress shards go through that same engine), so
/// requiring it here is not a back-compat hazard — it keeps the merged paired archive
/// direct-write-decodable instead of silently dropping the table and forcing the slower
/// part-file decode path.
fn merge_paired_core(inputs: &[PathBuf], out: &mut dyn Write) -> Result<()> {
    if inputs.is_empty() {
        bail!("qz merge: need at least one input archive");
    }

    // Pass 1: validate + load each paired input (front-header identity + paired directory
    // contract + the per-mate decoded-size table). No per-block CRC walk.
    let mut header_bytes: Option<Vec<u8>> = None;
    let mut header_end0: u64 = 0;
    struct PShard {
        path: PathBuf,
        dir: ChunkDirectory,
        // decoded_sizes_per_mate[mate][chunk] for this shard (always 2 mates).
        per_mate: Vec<Vec<u64>>,
    }
    let mut shards: Vec<PShard> = Vec::with_capacity(inputs.len());

    for path in inputs {
        let (hdr, header_end) = read_fixed_header(path)?;
        if hdr.archive_type != 1 {
            bail!(
                "qz merge: {} is not a paired archive (archive_type {})",
                path.display(),
                hdr.archive_type
            );
        }
        // Front-header byte-identity gate (same config across shards).
        let mut f = std::fs::File::open(path)?;
        let mut hbuf = vec![0u8; header_end as usize];
        f.read_exact(&mut hbuf)?;
        match &header_bytes {
            None => {
                header_bytes = Some(hbuf);
                header_end0 = header_end;
            }
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
        // Full paired structural validation (the same §-contract decode runs).
        super::paired::format::validate_paired_directory(&dir)?;
        // Validated, CRC-checked per-mate decoded sizes (rejects a corrupt table).
        let layout = read_chunk_layout(path)?;
        if layout.decoded_sizes_per_mate.len() != 2
            || layout.decoded_sizes_per_mate.iter().any(|m| m.len() != dir.num_chunks as usize)
        {
            bail!(
                "qz merge: {} has no ChunkDecodedSizes(n_mates=2) table; every paired input must have one",
                path.display()
            );
        }
        shards.push(PShard { path: path.clone(), dir, per_mate: layout.decoded_sizes_per_mate });
    }
    let header_bytes = header_bytes.expect("inputs non-empty");

    // Emit the shared front header.
    out.write_all(&header_bytes)?;
    let mut pos: u64 = header_end0;

    // Pass 2: copy per-chunk frames verbatim in FOOTER order; renumber chunk_index; mate,
    // role, codec, record_count carried verbatim. Concatenate the per-mate size tables into
    // one row-major array (sizes[chunk*2 + mate]) for the rebuilt global.
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
        // Append this shard's chunks row-major: [r1_c0, r2_c0, r1_c1, r2_c1, …].
        for c in 0..shard.dir.num_chunks as usize {
            merged_sizes.push(shard.per_mate[0][c]);
            merged_sizes.push(shard.per_mate[1][c]);
        }
        chunk_base = chunk_base.checked_add(shard.dir.num_chunks).ok_or_else(overflow)?;
        num_reads = num_reads.checked_add(shard.dir.num_reads).ok_or_else(overflow)?;
    }

    // Rebuild the single merged ChunkDecodedSizes(n_mates=2) global (every input had one →
    // unconditional). Mirrors the paired streaming encoder exactly.
    debug_assert_eq!(merged_sizes.len(), chunk_base as usize * 2);
    let payload = encode_decoded_sizes(2, chunk_base, &merged_sizes);
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
    super::paired::format::validate_paired_directory(&dir)?;
    write_footer_and_locator(out, &dir)?;
    out.flush()?;
    Ok(())
}

/// Single-end (`archive_type 0`) v5 archive merge: stitch N complete single-end
/// archives (`inputs`, ≥1) into one, streamed to `out`.
///
/// Every input must be a complete, well-formed single-end (`archive_type 0`) v5
/// archive carrying a `ChunkDecodedSizes` table, and all must share a
/// **byte-identical front header** (same compression config). The output decodes
/// losslessly to the concatenation of the inputs' reads and keeps a rebuilt
/// `ChunkDecodedSizes` global so NUMA direct-write decode still applies.
fn merge_single_core(inputs: &[PathBuf], out: &mut dyn Write) -> Result<()> {
    if inputs.is_empty() {
        bail!("qz merge: need at least one input archive");
    }

    // Pass 1: validate + load every input (full single-end contract + config gate).
    let mut header_bytes: Option<Vec<u8>> = None;
    let mut header_end0: u64 = 0;
    let mut has_rc_canon0 = false;
    let mut shards: Vec<Shard> = Vec::with_capacity(inputs.len());

    for path in inputs {
        let (hdr, header_end) = read_fixed_header(path)?;
        if hdr.archive_type != 0 {
            bail!(
                "qz merge: {} is not a single-end archive (archive_type {})",
                path.display(),
                hdr.archive_type
            );
        }

        // Front-header byte-identity gate (same config across shards).
        let mut f = std::fs::File::open(path)?;
        let mut hbuf = vec![0u8; header_end as usize];
        f.read_exact(&mut hbuf)?;
        match &header_bytes {
            None => {
                header_bytes = Some(hbuf);
                header_end0 = header_end;
                has_rc_canon0 = hdr.has_rc_canon();
            }
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
        // FULL single-end structural validation — the same contract decode runs
        // (inline frame cross-checks, codec uniformity, role/chunk pairing,
        // per-chunk record-count agreement, RcFlags rule, global contract).
        let max_block_len = decode_block_bounds(hdr.encoding_type).0;
        validate_single_end_directory(&dir, &mut f, max_block_len, hdr.has_rc_canon())?;

        // Validated, CRC-checked per-chunk decoded sizes (rejects a corrupt table).
        let layout = read_chunk_layout(path)?;
        if layout.per_chunk_decoded_bytes.len() != dir.num_chunks as usize {
            bail!(
                "qz merge: {} has no ChunkDecodedSizes table; every input must have one",
                path.display()
            );
        }
        shards.push(Shard {
            path: path.clone(),
            dir,
            decoded_sizes: layout.per_chunk_decoded_bytes,
        });
    }
    let header_bytes = header_bytes.expect("inputs is non-empty");

    // Emit the shared front header.
    out.write_all(&header_bytes)?;
    let mut pos: u64 = header_end0;

    // Pass 2: copy per-chunk block frames; rebuild directory + concat size tables.
    let mut merged: Vec<ChunkDirEntry> = Vec::new();
    let mut merged_sizes: Vec<u64> = Vec::new();
    let mut chunk_base: u32 = 0;
    let mut num_reads: u64 = 0;

    let overflow = || anyhow::anyhow!("qz merge: archive metadata overflow (forged input?)");

    for shard in &shards {
        let mut f = std::fs::File::open(&shard.path)?;
        merged_sizes.extend_from_slice(&shard.decoded_sizes);

        // Copy per-chunk frames in FOOTER ORDER (the order entries appear in the
        // directory — exactly what single-end decode consumes; `build_role_index`
        // preserves footer order and pairs records by position). Do NOT sort by
        // offset: a valid v5 archive need not have offset order == footer order, so
        // reordering could silently mis-pair records. For qz-produced archives the
        // two coincide, so reads stay sequential. Globals are excluded — the one
        // merged table is rebuilt below.
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

        chunk_base = chunk_base.checked_add(shard.dir.num_chunks).ok_or_else(overflow)?;
        num_reads = num_reads.checked_add(shard.dir.num_reads).ok_or_else(overflow)?;
    }

    // Rebuild the single merged ChunkDecodedSizes global (every input had one →
    // unconditional). Mirrors the single-end encoder exactly.
    debug_assert_eq!(merged_sizes.len(), chunk_base as usize);
    let payload = encode_decoded_sizes(1, chunk_base, &merged_sizes);
    write_role_blocks(
        &[(chunk_base, payload)],
        StreamRole::ChunkDecodedSizes,
        0, // mate
        0, // codec (raw metadata)
        GLOBAL_SENTINEL,
        out,
        &mut pos,
        &mut merged,
    )?;

    // Self-check the rebuilt directory with decode's own validators, then write tail.
    // The file-based inline-frame half is intentionally NOT re-run: merged frames are
    // verbatim copies of inputs already validated in pass 1, and `out` may be stdout.
    let dir = ChunkDirectory {
        num_reads,
        num_chunks: chunk_base,
        entries: merged,
    };
    validate_directory_spans(&dir.entries, header_end0, pos)?;
    validate_single_end_directory_structure(&dir, has_rc_canon0)?;
    write_footer_and_locator(out, &dir)?;
    out.flush()?;
    Ok(())
}

/// Shared atomic file-output wrapper for the merge entries: writes the result of
/// `f(inputs, out)` through an O_EXCL sibling temp and renames on success (so a
/// failed merge never leaves a partial archive), with the same
/// reject-existing-unless-`force` semantics as `compress()`. For streaming to an
/// arbitrary writer (e.g. stdout) call the writer-based entry (`merge_archives`)
/// directly.
fn merge_to_path_with(
    inputs: &[PathBuf],
    output: &Path,
    force: bool,
    f: impl FnOnce(&[PathBuf], &mut dyn Write) -> Result<()>,
) -> Result<()> {
    if output.exists() && !force {
        bail!(
            "Output file already exists: {}\nUse --force to overwrite",
            output.display()
        );
    }
    // Unpredictable O_EXCL sibling temp; rename atomically on success.
    let tmp = {
        let mut name = output
            .file_name()
            .ok_or_else(|| anyhow::anyhow!("output path has no file name: {}", output.display()))?
            .to_os_string();
        name.push(super::atomic_tmp_suffix());
        output.with_file_name(name)
    };
    struct TmpCleanup(PathBuf);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    let file = super::create_new_for_write(&tmp)
        .map_err(|e| anyhow::anyhow!("Failed to create temp output {}: {}", tmp.display(), e))?;
    let guard = TmpCleanup(tmp.clone());
    let mut out = BufWriter::new(file);
    f(inputs, &mut out)?;
    out.flush()?;
    drop(out);
    std::mem::forget(guard); // success: keep the renamed output
    if let Err(e) = std::fs::rename(&tmp, output) {
        let _ = std::fs::remove_file(&tmp);
        bail!("Failed to rename temp file to {}: {}", output.display(), e);
    }
    Ok(())
}

/// Thin wrapper over [`merge_archives`] retained for the `qz merge` CLI and the
/// existing single-end test suite. Errors if any input is not single-end, with the
/// same per-input "is not a single-end archive" message the single-end core uses.
pub fn merge_single_end_archives(inputs: &[PathBuf], out: &mut dyn Write) -> Result<()> {
    if inputs.is_empty() {
        bail!("qz merge: need at least one input archive");
    }
    for p in inputs {
        let (hdr, _) = read_fixed_header(p)?;
        if hdr.archive_type != 0 {
            bail!(
                "qz merge: {} is not a single-end archive (archive_type {})",
                p.display(),
                hdr.archive_type
            );
        }
    }
    merge_archives(inputs, out)
}

/// Merge N single-end v5 archives into the file at `output`, with the same
/// atomic-write + reject-existing-unless-`force` semantics as `compress()`:
/// writes through an O_EXCL sibling temp and renames on success, so a failed merge
/// never leaves a partial archive. For streaming to an arbitrary writer (e.g.
/// stdout) call [`merge_single_end_archives`] directly.
pub fn merge_single_end_archives_to_path(
    inputs: &[PathBuf],
    output: &Path,
    force: bool,
) -> Result<()> {
    merge_to_path_with(inputs, output, force, merge_single_end_archives)
}

/// Read a `[offset, offset+length)` block frame from `f` into a fresh buffer.
fn read_frame(f: &mut std::fs::File, offset: u64, length: u64) -> Result<Vec<u8>> {
    f.seek(SeekFrom::Start(offset))?;
    let mut buf = vec![0u8; length as usize];
    f.read_exact(&mut buf)?;
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{AdvancedOptions, CompressConfig};
    use crate::compression::chunk_directory::{GLOBAL_SENTINEL, StreamRole, read_v5_footer};
    use crate::compression::decompress_impl::read_fixed_header;

    /// A CRC-corrupt source ChunkDecodedSizes table must be rejected by the merge
    /// (via read_chunk_layout's frame-CRC check), never laundered into the output.
    #[test]
    fn merge_rejects_corrupt_source_decoded_sizes_table() {
        let td = tempfile::TempDir::new().unwrap();
        let d = td.path();

        // Build a small single-end archive (chunk_records = 4 → a few chunks).
        let input = d.join("in.fastq");
        let mut s = String::new();
        for i in 0..8 {
            s.push_str(&format!("@r{i}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"));
        }
        std::fs::write(&input, &s).unwrap();
        let arc = d.join("a.qz");
        let cfg = CompressConfig {
            input: vec![input],
            output: arc.clone(),
            working_dir: d.to_path_buf(),
            threads: 1,
            advanced: AdvancedOptions { chunk_records: 4, ..AdvancedOptions::default() },
            ..CompressConfig::default()
        };
        crate::compression::compress(&cfg).unwrap();

        // Locate the ChunkDecodedSizes table frame via the footer (pub(crate) here).
        let (_, header_end) = read_fixed_header(&arc).unwrap();
        let dir = read_v5_footer(&arc, header_end).unwrap();
        let g = dir
            .entries
            .iter()
            .find(|e| e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::ChunkDecodedSizes)
            .expect("single-end archive has a ChunkDecodedSizes global");

        // Flip a byte INSIDE the payload (the 12-byte frame header precedes it).
        let mut bytes = std::fs::read(&arc).unwrap();
        bytes[g.offset as usize + 12] ^= 0xFF;
        let corrupt = d.join("corrupt.qz");
        std::fs::write(&corrupt, &bytes).unwrap();

        // Merge must reject it.
        let out = d.join("out.qz");
        let mut w = std::io::BufWriter::new(std::fs::File::create(&out).unwrap());
        let r = merge_single_end_archives(&[corrupt], &mut w);
        // Pin the failure to the table CRC path (read_chunk_layout's frame-CRC
        // check emits "block frame CRC32 mismatch"), not an unrelated error.
        let err = r.unwrap_err().to_string();
        assert!(err.contains("CRC"), "expected a table CRC failure, got: {err}");
    }

    /// Build a tiny single-end fqz archive (`Fqzcomp` quality, small fqz sub-block
    /// size) so a chunk spans MULTIPLE `Qual` sub-blocks — the production NUMA-shard
    /// shape. `chunk_records = 8`, `quality_ctx_block_size = 2` → many Qual frames.
    fn build_fqz_archive(d: &Path, name: &str, n: usize) -> PathBuf {
        use crate::cli::QualityCompressor;
        let input = d.join(format!("{name}.fastq"));
        let mut s = String::new();
        for i in 0..n {
            // Vary quality so fqz has something to model.
            let q: String = (0..12).map(|j| (b'I' - ((i + j) % 5) as u8) as char).collect();
            s.push_str(&format!("@r{i}\nACGTACGTACGT\n+\n{q}\n"));
        }
        std::fs::write(&input, &s).unwrap();
        let arc = d.join(format!("{name}.qz"));
        let cfg = CompressConfig {
            input: vec![input],
            output: arc.clone(),
            working_dir: d.to_path_buf(),
            threads: 1,
            advanced: AdvancedOptions {
                chunk_records: 8,
                quality_compressor: QualityCompressor::Fqzcomp,
                quality_ctx_block_size: 2,
                ..AdvancedOptions::default()
            },
            ..CompressConfig::default()
        };
        crate::compression::compress(&cfg).unwrap();
        arc
    }

    /// Count `Qual` directory entries belonging to `chunk_index`.
    fn qual_entries_in_chunk(dir: &ChunkDirectory, chunk_index: u32) -> usize {
        dir.entries
            .iter()
            .filter(|e| e.role == StreamRole::Qual && e.chunk_index == chunk_index)
            .count()
    }

    /// The merge must (a) only ever be tested against the production fqz shape where
    /// a chunk spans >1 Qual sub-block, and (b) preserve that shape — every shard's
    /// Qual frames are relocated verbatim, so the merged chunk 0 still has >1 Qual
    /// frame and chunk indices are renumbered contiguously across shards.
    #[test]
    fn merge_preserves_multi_qual_subblocks_per_chunk() {
        let td = tempfile::TempDir::new().unwrap();
        let d = td.path();

        // Shard A: 8 reads → 1 chunk, multiple Qual sub-blocks (qblock = 2).
        let a = build_fqz_archive(d, "fa", 8);
        let (_, header_end_a) = read_fixed_header(&a).unwrap();
        let dir_a = read_v5_footer(&a, header_end_a).unwrap();
        assert!(
            qual_entries_in_chunk(&dir_a, 0) > 1,
            "shard must have >1 Qual sub-block in chunk 0 (production fqz shape); got {}",
            qual_entries_in_chunk(&dir_a, 0)
        );

        // Shard B: 4 reads → 1 chunk, also multiple Qual sub-blocks.
        let b = build_fqz_archive(d, "fb", 4);
        let (_, header_end_b) = read_fixed_header(&b).unwrap();
        let dir_b = read_v5_footer(&b, header_end_b).unwrap();
        assert!(
            qual_entries_in_chunk(&dir_b, 0) > 1,
            "second shard must also have >1 Qual sub-block in chunk 0"
        );

        // Merge the two shards.
        let merged = d.join("merged_fqz.qz");
        let mut w = std::io::BufWriter::new(std::fs::File::create(&merged).unwrap());
        merge_single_end_archives(&[a.clone(), b.clone()], &mut w).unwrap();
        use std::io::Write;
        w.flush().unwrap();
        drop(w);

        // The merged archive must preserve the multi-Qual-per-chunk shape AND
        // renumber chunk indices contiguously (shard B's chunk → merged chunk 1).
        let (_, header_end_m) = read_fixed_header(&merged).unwrap();
        let dir_m = read_v5_footer(&merged, header_end_m).unwrap();
        assert!(
            qual_entries_in_chunk(&dir_m, 0) > 1,
            "merged chunk 0 must retain >1 Qual sub-block; got {}",
            qual_entries_in_chunk(&dir_m, 0)
        );
        assert!(
            qual_entries_in_chunk(&dir_m, 1) > 1,
            "merged chunk 1 (from shard B) must be present with >1 Qual sub-block; got {}",
            qual_entries_in_chunk(&dir_m, 1)
        );
        assert_eq!(
            dir_m.num_chunks,
            dir_a.num_chunks + dir_b.num_chunks,
            "merged num_chunks must equal the sum of the shards'"
        );
        assert_eq!(
            dir_m.num_reads,
            dir_a.num_reads + dir_b.num_reads,
            "merged num_reads must equal the sum of the shards'"
        );

        // Chunk indices present in the merged directory must be exactly 0..num_chunks.
        let mut chunks: Vec<u32> = dir_m
            .entries
            .iter()
            .filter(|e| e.chunk_index != GLOBAL_SENTINEL)
            .map(|e| e.chunk_index)
            .collect();
        chunks.sort_unstable();
        chunks.dedup();
        let expected: Vec<u32> = (0..dir_m.num_chunks).collect();
        assert_eq!(chunks, expected, "merged chunk indices must be contiguous 0..num_chunks");
    }
}
