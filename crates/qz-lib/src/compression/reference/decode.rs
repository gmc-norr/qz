//! Reference-direct paired-end sequence DECODE + verify (`--reference`, spec §4).
//!
//! The decompress half: reconstruct every read from the backing +
//! per-(chunk,mate) fallback pools alone — the reference is never required.
//! Bounded-RAM + seekable (`RefReader`): globals load once, per-chunk entries
//! read on demand, chunks reconstruct bounded-parallel into an in-order writer.
//! The encode half lives in [`super::encode`].
use crate::cli::DecompressConfig;
use anyhow::{Result, bail};

use super::backing::{IntervalMap, PackedBacking};
use super::edits::{decode_edits, reconstruct_read_from_backing};
use super::format::{validate_reference_directory, validate_reference_directory_single};
use super::{backing, fallback, positions, refmeta};
use crate::compression::chunk_directory::{
    ChunkDirEntry, ChunkDirectory, GLOBAL_SENTINEL as DIR_GLOBAL_SENTINEL, StreamRole,
};
use crate::compression::dna_utils::read_varint;
use crate::compression::{FixedHeader, QualityBinning, code_to_binning};
// The shared single-end block-streaming engine primitives reference decode is
// converged onto: flat per-mate header/qual stream producers + cursors, and the
// batch-parallel FASTQ assembler. Sequences remain reference-specific (the
// chunk-granular edit-reconstruction in `reconstruct_mate`), but headers,
// qualities, and emission now run through the SAME engine as single-end/paired.
use crate::compression::decompress_impl::{
    FORMAT_BATCH_SIZE, StreamCodec, StreamCursor, StreamIndex, WRITE_BATCH, build_quality_inputs,
    build_role_index_mate_range, decode_block_bounds, drain_one_qual, flush_remaining,
    format_batch_parallel, push_and_flush, spawn_stream_producer,
};

/// Decompress a paired-reference archive to two FASTQ files. Opens the archive
/// and reads the footer + globals, then streams each chunk's entries on demand
/// (seekable, bounded-RAM) and reconstructs via `decode_reference_streaming`,
/// which streams each mate straight into its file sink (the proven reference-free
/// decode) without materialising a full-FASTQ Vec. Writes
/// `<prefix>_R1.fastq` / `<prefix>_R2.fastq` atomically (mirrors
/// `decompress_paired_v5`).
pub fn decompress_reference(a: &DecompressConfig) -> Result<()> {
    use crate::compression::paired::{TMP_NONCE, TmpPair, publish_pair, with_suffix};
    use std::io::Write;
    use std::sync::atomic::Ordering;

    if a.output.is_empty() {
        bail!("reference decompress needs an output prefix");
    }
    let prefix = &a.output[0];
    if crate::cli::is_stdio_path(prefix) {
        bail!("reference decompress needs a named output prefix (no stdout)");
    }
    if a.gzipped {
        bail!("gzipped output not supported in reference mode (writes plain _R1.fastq/_R2.fastq)");
    }
    let out1 = with_suffix(prefix, "_R1.fastq");
    let out2 = with_suffix(prefix, "_R2.fastq");

    if !a.force && (out1.exists() || out2.exists()) {
        bail!(
            "Output file already exists: {} or {}\nUse --force to overwrite",
            out1.display(),
            out2.display()
        );
    }

    // Seekable decode: open the file + read/validate the footer, then load the
    // small whole-archive globals once and read each per-chunk entry on demand via
    // positioned reads — the archive is never fully materialised. The
    // reconstructed FASTQ output is likewise streamed below (no full r1/r2 Vecs).
    let (reader, dir, hdr) = open_ref(&a.input)?;

    // Atomic two-file publish (mirrors decompress_paired_v5). One nonce per call
    // shared by tmp1/tmp2 + publish_pair's .bak names.
    let pid = std::process::id();
    let nonce = TMP_NONCE.fetch_add(1, Ordering::Relaxed);
    let tmp1 = with_suffix(&out1, &format!(".{pid}.{nonce}.tmp"));
    let tmp2 = with_suffix(&out2, &format!(".{pid}.{nonce}.tmp"));
    let mut guard = TmpPair {
        a: tmp1.clone(),
        b: tmp2.clone(),
        armed: true,
    };
    // O_EXCL parity with the paired path: names are already pid+nonce unique, so
    // create_new also rejects a pre-planted symlink/file. Stream each mate's
    // reconstructed FASTQ straight into its BufWriter — no full-FASTQ Vec.
    let mut w1 = std::io::BufWriter::new(crate::compression::create_new_for_write(&tmp1)?);
    let mut w2 = std::io::BufWriter::new(crate::compression::create_new_for_write(&tmp2)?);
    // Engine-aligned decode: headers + qualities stream through the shared
    // block-parallel producers/cursors while sequences are reconstructed
    // chunk-by-chunk (par_iter within a chunk) and spliced via
    // `format_batch_parallel`. Bounded-RAM: globals (resident) + one chunk's
    // reconstructed seqs + the producers' bounded channels — constant in read count.
    decode_reference_streaming(
        &a.input,
        reader,
        &dir,
        &hdr,
        &mut RefSink::Two(&mut w1, &mut w2),
        None,
    )?;
    w1.flush()?;
    w2.flush()?;
    drop(w1);
    drop(w2);

    publish_pair(&tmp1, &tmp2, &out1, &out2, nonce)?;
    guard.disarm();
    Ok(())
}

/// Decode a paired-reference (`archive_type` 2) v5 archive as a single
/// **interleaved** FASTQ stream — R1/R2 records alternating — to `a.output[0]`
/// (a file, `-`/stdout, or gzipped). Reuses `decode_reference_streaming` via
/// `RefSink::Interleaved`; the single-sink dispatch owns the atomic-temp / stdout
/// / gzip output machinery (`force`/exists is checked by the caller).
pub fn decompress_reference_interleaved(a: &DecompressConfig) -> Result<()> {
    use std::time::Instant;

    if a.output.is_empty() {
        bail!("interleaved reference decompress needs one output (-o <file>|-)");
    }
    let (reader, dir, hdr) = open_ref(&a.input)?;
    // R1 + R2 records (cosmetic, for the progress log).
    let num_records = dir.num_reads.saturating_mul(2) as usize;
    let input = a.input.clone();
    crate::compression::decompress_impl::decompress_streaming_bsc_dispatch(
        a,
        Instant::now(),
        num_records,
        move |w| {
            decode_reference_streaming(
                &input,
                reader,
                &dir,
                &hdr,
                &mut RefSink::Interleaved(w),
                None,
            )
        },
    )
}

/// Decompress a single-end reference archive (`archive_type=4`) to ONE FASTQ.
/// Opens the archive, reads the footer + globals, streams each chunk's mate-1
/// entries on demand, reconstructs reference-free, and writes the single output
/// atomically (sibling O_EXCL `.tmp` → rename). The output path is the user's
/// `output[0]` verbatim (no `_R1`/`_R2` suffix — single-end produces one file).
pub fn decompress_reference_single(a: &DecompressConfig) -> Result<()> {
    use std::io::Write;
    if a.output.is_empty() {
        bail!("reference decompress needs an output path");
    }
    let out = &a.output[0];
    if crate::cli::is_stdio_path(out) {
        bail!("reference decompress needs a named output path (no stdout)");
    }
    if a.gzipped {
        bail!("gzipped output not supported in reference mode (writes plain .fastq)");
    }
    if !a.force && out.exists() {
        bail!(
            "Output file already exists: {}\nUse --force to overwrite",
            out.display()
        );
    }

    let (reader, dir, hdr) = open_ref_single(&a.input)?;

    // Atomic single-file publish: stream to a sibling O_EXCL `.tmp`, rename on
    // success, drop-guard removes the temp on any abort (mirrors single-end
    // default's atomic commit, reusing the shared `atomic_tmp_suffix`).
    let mut tmp_os = out.clone().into_os_string();
    tmp_os.push(crate::compression::atomic_tmp_suffix());
    let tmp: std::path::PathBuf = tmp_os.into();
    struct TmpCleanup(std::path::PathBuf, bool);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            if self.1 {
                let _ = std::fs::remove_file(&self.0);
            }
        }
    }
    let mut guard = TmpCleanup(tmp.clone(), true);

    let mut w = std::io::BufWriter::new(crate::compression::create_new_for_write(&tmp)?);
    decode_reference_single_streaming(&a.input, reader, &dir, &hdr, &mut w, None)?;
    w.flush()?;
    drop(w);

    std::fs::rename(&tmp, out)
        .map_err(|e| anyhow::anyhow!("Failed to rename temp file to {}: {e}", out.display()))?;
    guard.1 = false; // published — disarm the cleanup
    Ok(())
}

/// Decode chunks `[a, b)` of a reference v5 archive into `out_parts[0]` (R1) and
/// `out_parts[1]` (R2) for the range-restricted NUMA-shard decode (`decode_chunk_range`).
/// Mirrors `decompress_reference`'s setup (`open_ref` → footer + globals + reference
/// validator), but writes the reconstructed FASTQ straight into `out_parts` (no
/// temp/publish) and threads `Some((a, b))` into the shared streaming core so only the
/// requested chunks are walked. Concatenating the parts in chunk order byte-equals a
/// full decode: the range-filtered header/quality cursors start at chunk `a` and the
/// chunk-granular sequence columns are read positionally, so per-record alignment is
/// preserved. The reference path has its own validation (`open_ref` /
/// `validate_reference_directory`), so this does NOT call `require_streamable_v5`.
pub(crate) fn decode_reference_range(
    archive: &std::path::Path,
    a: u32,
    b: u32,
    out_parts: &[std::path::PathBuf],
    regions: &[crate::compression::DirectWriteRegion],
) -> Result<()> {
    use std::io::Write;
    if out_parts.len() != 2 {
        bail!(
            "reference decode_chunk_range expects 2 output parts, got {}",
            out_parts.len()
        );
    }
    // Seekable open: file handle + length, validated footer, reference semantic
    // validator (same as `decompress_reference`). Returns `(reader, dir, hdr)`.
    let (reader, dir, hdr) = open_ref(archive)?;

    if regions.is_empty() {
        // Part-file path: each mate written from byte 0; the driver concatenates.
        let f1 = crate::compression::create_new_for_write(&out_parts[0])?;
        let f2 = crate::compression::create_new_for_write(&out_parts[1])?;
        let mut w1 = std::io::BufWriter::new(f1);
        let mut w2 = std::io::BufWriter::new(f2);
        decode_reference_streaming(
            archive,
            reader,
            &dir,
            &hdr,
            &mut RefSink::Two(&mut w1, &mut w2),
            Some((a, b)),
        )?;
        w1.flush()?;
        w2.flush()?;
        return Ok(());
    }

    // Direct-write path: verify each mate's region against its per-mate decoded-size
    // slice, seek into the two pre-sized shared outputs, decode through bounded sinks
    // (identical reconstruction — only the sink differs). Any mismatch/overrun/short
    // write is a DirectWriteIntegrityError so `auto` falls back to a correct in-process
    // decode. Mirrors `decode_chunk_range_paired`'s direct branch.
    let ierr = |m: String| -> anyhow::Error { crate::compression::DirectWriteIntegrityError(m).into() };
    if regions.len() != 2 {
        return Err(ierr(format!("paired reference direct-write expects 2 regions, got {}", regions.len())));
    }
    let sizes = ref_decoded_sizes_per_mate(archive, 2)?;
    let (mut w1, exp1) = open_direct_sink(&out_parts[0], &sizes[0], a, b, regions[0])?;
    let (mut w2, exp2) = open_direct_sink(&out_parts[1], &sizes[1], a, b, regions[1])?;
    let res = decode_reference_streaming(
        archive,
        reader,
        &dir,
        &hdr,
        &mut RefSink::Two(&mut w1, &mut w2),
        Some((a, b)),
    );
    if w1.exceeded {
        return Err(ierr(format!("R1 decoded output exceeded region {exp1} (ChunkDecodedSizes underestimate)")));
    }
    if w2.exceeded {
        return Err(ierr(format!("R2 decoded output exceeded region {exp2} (ChunkDecodedSizes underestimate)")));
    }
    res?;
    w1.flush()?;
    w2.flush()?;
    if w1.written != exp1 {
        return Err(ierr(format!("R1 wrote {} != region {}", w1.written, exp1)));
    }
    if w2.written != exp2 {
        return Err(ierr(format!("R2 wrote {} != region {}", w2.written, exp2)));
    }
    Ok(())
}

/// Type alias for a reference direct-write sink: the shared `BoundedCountingWriter`
/// over a `BufWriter<File>` (caps the decode write at the table-derived region length).
type DirectSink = crate::compression::decompress_impl::BoundedCountingWriter<std::io::BufWriter<std::fs::File>>;

/// Load + validate a reference archive's per-mate `ChunkDecodedSizes` table (classifying
/// a missing/corrupt table as a `DirectWriteIntegrityError`, so `auto` falls back). The
/// returned outer length is exactly `n_mates`.
fn ref_decoded_sizes_per_mate(archive: &std::path::Path, n_mates: usize) -> Result<Vec<Vec<u64>>> {
    let ierr = |m: String| -> anyhow::Error { crate::compression::DirectWriteIntegrityError(m).into() };
    let layout = crate::compression::decompress_impl::read_chunk_layout_impl(archive)
        .map_err(|e| ierr(format!("layout read failed: {e}")))?;
    if layout.decoded_sizes_per_mate.len() != n_mates {
        return Err(ierr(format!(
            "reference archive has no {n_mates}-mate ChunkDecodedSizes table (mates={})",
            layout.decoded_sizes_per_mate.len()
        )));
    }
    Ok(layout.decoded_sizes_per_mate)
}

/// Verify one mate's direct-write `region` against its per-chunk `sizes`, open the
/// pre-sized output, seek to the region base, and wrap it in a bounded sink. Returns
/// the sink + the table-derived expected length (a region mismatch is a
/// `DirectWriteIntegrityError` via `verify_direct_region`). Shared by the paired and
/// single-end reference direct paths.
fn open_direct_sink(
    out_path: &std::path::Path,
    sizes: &[u64],
    a: u32,
    b: u32,
    region: crate::compression::DirectWriteRegion,
) -> Result<(DirectSink, u64)> {
    use std::io::{Seek, SeekFrom};
    let exp = crate::compression::decompress_impl::verify_direct_region(sizes, a, b, region)?;
    let mut f = std::fs::OpenOptions::new().write(true).open(out_path)?;
    f.seek(SeekFrom::Start(region.offset))?;
    let sink = crate::compression::decompress_impl::BoundedCountingWriter::new(
        std::io::BufWriter::new(f),
        exp,
    );
    Ok((sink, exp))
}

/// Decode chunks `[a, b)` of a single-end reference (archive_type 4) archive into
/// `out_parts[0]` for the range-restricted NUMA-shard decode. The single-end analogue
/// of [`decode_reference_range`]: same setup via `open_ref_single` (footer + globals +
/// the single-mate reference validator `validate_reference_directory_single`, asserting
/// `archive_type == 4` — `open_ref` itself hard-bails on any type ≠ 2, so type-4 MUST
/// open through `open_ref_single`), but writes ONE reconstructed FASTQ (mate 1 only) and
/// threads `Some((a, b))` into the streaming core so only the requested chunks are walked.
/// Concatenating the parts in chunk order byte-equals a full decode (the range-filtered
/// header/quality cursors start at chunk `a`; the chunk-granular sequence columns are read
/// positionally). The reference path has its own validation, so this does NOT call
/// `require_streamable_v5`.
pub(crate) fn decode_reference_range_single(
    archive: &std::path::Path,
    a: u32,
    b: u32,
    out_parts: &[std::path::PathBuf],
    regions: &[crate::compression::DirectWriteRegion],
) -> Result<()> {
    use std::io::Write;
    if out_parts.len() != 1 {
        bail!(
            "single-end reference decode_chunk_range expects 1 output part, got {}",
            out_parts.len()
        );
    }
    // `open_ref_single` asserts archive_type == 4 and runs the single-mate validator.
    // (Calling `open_ref` here would wrongly bail — it requires archive_type == 2.)
    let (reader, dir, hdr) = open_ref_single(archive)?;

    if regions.is_empty() {
        // Part-file path: written from byte 0; the driver concatenates.
        let f1 = crate::compression::create_new_for_write(&out_parts[0])?;
        let mut w1 = std::io::BufWriter::new(f1);
        decode_reference_single_streaming(archive, reader, &dir, &hdr, &mut w1, Some((a, b)))?;
        w1.flush()?;
        return Ok(());
    }

    // Direct-write path: ONE mate, one pre-sized output (mirrors the paired branch with
    // a single sink).
    let ierr = |m: String| -> anyhow::Error { crate::compression::DirectWriteIntegrityError(m).into() };
    if regions.len() != 1 {
        return Err(ierr(format!("single-end reference direct-write expects 1 region, got {}", regions.len())));
    }
    let sizes = ref_decoded_sizes_per_mate(archive, 1)?;
    let (mut w1, exp1) = open_direct_sink(&out_parts[0], &sizes[0], a, b, regions[0])?;
    let res = decode_reference_single_streaming(archive, reader, &dir, &hdr, &mut w1, Some((a, b)));
    if w1.exceeded {
        return Err(ierr(format!("decoded output exceeded region {exp1} (ChunkDecodedSizes underestimate)")));
    }
    res?;
    w1.flush()?;
    if w1.written != exp1 {
        return Err(ierr(format!("wrote {} != region {}", w1.written, exp1)));
    }
    Ok(())
}

/// Find a global entry by role (exactly one exists post-validate). A global
/// entry carries `chunk_index == GLOBAL_SENTINEL`; reference globals map 1:1 to
/// their same-named `StreamRole`, so the role byte alone disambiguates.
fn find_global(footer: &ChunkDirectory, role: StreamRole) -> Result<&ChunkDirEntry> {
    footer
        .entries
        .iter()
        .find(|e| e.chunk_index == DIR_GLOBAL_SENTINEL && e.role == role)
        .ok_or_else(|| anyhow::anyhow!("missing global role {role:?}"))
}

/// A reference archive open for **seekable** decode: the file handle plus its
/// length, captured ONCE at open time. Bundling the two means the per-entry
/// bounds check reuses the hoisted `len` instead of an `fstat` per entry, and
/// Phase 4 can share ONE thing across rayon workers (`Arc<RefReader>`) — the
/// positioned reads carry no shared seek cursor, so concurrent `read_entry`
/// calls are safe.
struct RefReader {
    file: std::fs::File,
    /// File length captured once in `open_ref`. The footer validator
    /// (`validate_and_parse_footer`) already proved every entry's
    /// `[offset, offset+length)` lies within `[header_end, footer_start)` (and
    /// `footer_start <= len`), so this is defense-in-depth, not the primary guard.
    len: u64,
}

impl RefReader {
    /// Read exactly `length` bytes at absolute `offset` into a fresh `Vec`. The
    /// seekable-decode reader: each v5 directory entry is read on demand here
    /// instead of slicing a whole-archive `Vec` into memory.
    ///
    /// **Positioned reads only** (`read_exact_at` on unix, no shared seek cursor)
    /// so the same `&RefReader` can be shared across rayon workers (Phase 4). Do
    /// NOT introduce a `File::seek` on this path.
    ///
    /// **Untrusted footer:** the footer CRC is recomputable, not a MAC, so
    /// `offset` and `length` are attacker-controlled. Before allocating `length`
    /// bytes we bound `offset.checked_add(length) <= self.len` (matching the
    /// byte-slice bounds guarantee the in-memory reader gave; `self.len` is the
    /// file length captured once in `open_ref`, NOT a per-call `fstat`) and
    /// `usize::try_from(length)` so a hostile huge `length` errors cleanly instead
    /// of triggering a multi-TB allocation.
    fn read_entry(&self, offset: u64, length: u64) -> Result<Vec<u8>> {
        let file_len = self.len;
        let end = offset
            .checked_add(length)
            .ok_or_else(|| anyhow::anyhow!("entry offset+length overflow"))?;
        if end > file_len {
            bail!("entry payload [{offset}, {end}) out of bounds (file len {file_len})");
        }
        let len = usize::try_from(length)
            .map_err(|_| anyhow::anyhow!("entry length {length} exceeds usize"))?;
        let mut buf = vec![0u8; len];
        #[cfg(unix)]
        {
            use std::os::unix::fs::FileExt;
            // Fills the whole buffer or errors on early EOF; positioned (no shared
            // cursor), so it's concurrency-safe.
            self.file.read_exact_at(&mut buf, offset)?;
        }
        #[cfg(not(unix))]
        {
            // Portable fallback: positioned reads where available, else seek+read.
            // Loops until the buffer is full (a single positioned read may short-read).
            #[cfg(windows)]
            {
                use std::os::windows::fs::FileExt;
                let mut filled = 0usize;
                while filled < len {
                    let n = self
                        .file
                        .seek_read(&mut buf[filled..], offset + filled as u64)?;
                    if n == 0 {
                        bail!("entry payload early EOF at {}", offset + filled as u64);
                    }
                    filled += n;
                }
            }
            #[cfg(not(windows))]
            {
                use std::io::{Read, Seek, SeekFrom};
                // No positioned-read API: clone the handle so we don't perturb a
                // shared cursor, then seek+read_exact. (Non-unix, non-windows is not
                // a production target — unix is.)
                let mut f = self.file.try_clone()?;
                f.seek(SeekFrom::Start(offset))?;
                f.read_exact(&mut buf)?;
            }
        }
        Ok(buf)
    }
}

/// Read one BSC role entry's payload on demand (positioned read) and decode the
/// single v4 BSC block back to its raw bytes.
fn decode_bsc_role_file(reader: &RefReader, e: &ChunkDirEntry) -> Result<Vec<u8>> {
    // `record_count` is an untrusted footer field; clamp to usize before it sizes
    // the per-stream decode cap.
    let cap = crate::compression::codecs::stream_decode_cap(usize::try_from(e.record_count)?);
    let payload = reader.read_entry(e.offset, e.length)?;
    crate::compression::paired::streams::decode_bsc_stream(&payload, cap)
}

/// The reference archive's whole-archive globals, decoded ONCE and held resident
/// for the lifetime of a decode. Each field is fully OWNED (none borrow the
/// archive), so the decoder can read the small global entries up front and then
/// stream chunks against this struct with peak RAM bounded by one chunk's working
/// set + these globals (constant in read count).
struct Globals {
    backing: PackedBacking,
    imap: IntervalMap,
    meta: refmeta::ReferenceMeta,
}

/// Decode the 4 global entries (PackedBacking, NBitmap, IntervalMap,
/// ReferenceMeta) ONCE from `reader` and assemble a resident `Globals`. Shared by
/// the file-based decode (`decode_reference_streaming`) + deep-verify paths.
/// Reads each global on demand via `RefReader::read_entry` (positioned reads;
/// same roles, same order, same validation).
fn decode_globals(reader: &RefReader, footer: &ChunkDirectory) -> Result<Globals> {
    let e_back = find_global(footer, StreamRole::PackedBacking)?;
    // Footer counts are untrusted (the CRC is recomputable, not a MAC). The
    // backing base count drives a `vec![0u8; total_bases/4]` allocation inside
    // `read_blocks_into`; clamp to usize so a hostile u64 errors cleanly rather
    // than truncating on a 32-bit target.
    let total_bases = usize::try_from(e_back.record_count).map_err(|_| {
        anyhow::anyhow!(
            "reference backing base count {} exceeds usize",
            e_back.record_count
        )
    })?;
    let e_nbm = find_global(footer, StreamRole::NBitmap)?;
    let backing = backing::read_packed_backing_with_bitmap(
        &reader.read_entry(e_back.offset, e_back.length)?,
        &reader.read_entry(e_nbm.offset, e_nbm.length)?,
        total_bases,
    )?;

    let e_imap = find_global(footer, StreamRole::IntervalMap)?;
    let imap =
        backing::read_intervalmap_blockstream(&reader.read_entry(e_imap.offset, e_imap.length)?)?;
    // The interval map must cover exactly the bases the backing claims.
    if imap.total_bases() != total_bases as u64 {
        bail!(
            "interval map covers {} bases but backing has {}",
            imap.total_bases(),
            total_bases
        );
    }

    let e_meta = find_global(footer, StreamRole::ReferenceMeta)?;
    let meta =
        refmeta::read_refmeta_blockstream(&reader.read_entry(e_meta.offset, e_meta.length)?)?;

    Ok(Globals {
        backing,
        imap,
        meta,
    })
}

/// Per-mate flat streamable indices for the engine-aligned reference decode. Unlike
/// paired's `PairedMateIndices` there is NO sequence stream (sequences are
/// reconstructed from the chunk-granular columns by `reconstruct_mate`) and never
/// `RcFlags` (strand lives in the per-chunk `Strands` bitplane, applied inside
/// `reconstruct_read`). The `headers` index is the **columnar** stream; it is empty
/// (and unused) when the archive uses R2 `HeaderDelta`, in which case both mates'
/// headers are reconstructed per chunk via `decode_chunk_headers_delta` instead.
struct RefMateStreams {
    headers: StreamIndex,
    qualities: Option<StreamIndex>,
    /// `Fqz` (CODEC_FQZCOMP → decoded ASCII) or `Bsc` (CODEC_BSC → bit-packed).
    /// `None` for a no-quality archive.
    quality_codec: Option<StreamCodec>,
}

/// Build one mate's flat header + optional quality stream indices by concatenating
/// its per-chunk directory entries in chunk order (`build_role_index_mate_range` sorts
/// by chunk_index and expands the shared `write_block_stream` segment framing into the
/// per-block entries the engine producers consume). Reference quality is fqzcomp
/// (production) or BSC (sub-threshold / lossy); the removed quality_ctx codec
/// (code 4) is rejected (same as paired/single-end).
///
/// `chunk_range = Some((a, b))` restricts the index to chunks `[a, b)` (the
/// range-restricted NUMA-shard decode path); `None` is the full archive
/// (byte-identical to a full-archive `build_role_index_mate_range(.., None)`). The `has_qual`
/// presence probe scans the WHOLE directory (a Qual entry for ANY chunk of this mate
/// means quality is present for the range too — reference quality is uniform per
/// archive), so the range only filters which blocks the producers consume.
fn build_ref_mate_streams(
    archive_path: &std::path::Path,
    dir: &ChunkDirectory,
    mate: u8,
    chunk_range: Option<(u32, u32)>,
) -> Result<RefMateStreams> {
    let headers =
        build_role_index_mate_range(archive_path, dir, StreamRole::Headers, mate, chunk_range)?;
    let has_qual = dir
        .entries
        .iter()
        .any(|e| e.role == StreamRole::Qual && e.mate == mate);
    let (qualities, quality_codec) = if has_qual {
        use crate::compression::codec_ids::{CODEC_BSC, CODEC_FQZCOMP, CODEC_QUALITY_CTX};
        let codec_byte = dir
            .entries
            .iter()
            .find(|e| e.role == StreamRole::Qual && e.mate == mate)
            .map(|e| e.codec)
            .expect("has_qual implies a Qual entry");
        let codec = match codec_byte {
            c if c == CODEC_FQZCOMP => StreamCodec::Fqz,
            c if c == CODEC_BSC => StreamCodec::Bsc,
            c if c == CODEC_QUALITY_CTX => bail!(
                "reference mate {mate}: quality_ctx quality backend has been removed (code 4) — recompress"
            ),
            c => bail!("reference mate {mate}: unknown quality codec byte {c}"),
        };
        (
            Some(build_role_index_mate_range(
                archive_path,
                dir,
                StreamRole::Qual,
                mate,
                chunk_range,
            )?),
            Some(codec),
        )
    } else {
        (None, None)
    };
    Ok(RefMateStreams {
        headers,
        qualities,
        quality_codec,
    })
}

/// Emit one chunk's already-reconstructed `seqs` (RC already applied by
/// `reconstruct_read`) as FASTQ by pulling matching headers + qualities from the
/// flat per-mate engine cursors and formatting in `FORMAT_BATCH_SIZE` sub-batches
/// via `format_batch_parallel`. Mirrors the per-mate body of
/// `bounded_write_records_paired`, minus the sequence cursor (sequences come from
/// `seqs`) and minus RC (already applied). `record_base` is the global read index of
/// `seqs[0]`, used only for error messages.
/// Headers come from EITHER the streaming columnar cursor (`h_cur`, the common path)
/// OR a pre-decoded per-chunk id slice (`decoded_headers`, the R2-HeaderDelta path
/// where R2 ids are reconstructed against R1 before formatting — see
/// `decode_chunk_headers_delta`). Exactly one is `Some`.
#[allow(clippy::too_many_arguments)]
fn emit_chunk_mate(
    seqs: &[Vec<u8>],
    mut h_cur: Option<&mut StreamCursor>,
    decoded_headers: Option<&[Vec<u8>]>,
    mut q_cur: Option<&mut StreamCursor>,
    is_fqz: bool,
    const_qual_len: usize,
    bits_per_qual: usize,
    quality_binning: QualityBinning,
    is_fasta: bool,
    out: &mut dyn std::io::Write,
    buf: &mut Vec<u8>,
    record_base: u64,
) -> Result<()> {
    debug_assert!(h_cur.is_some() ^ decoded_headers.is_some(), "exactly one header source");
    let no_rc: Vec<u8> = Vec::new(); // reference applies RC inside reconstruct_read
    let est_per_record = 80 + 150 + 4 + 150 + 1;
    let cap = FORMAT_BATCH_SIZE.min(seqs.len().max(1));
    let mut cursor_headers: Vec<Vec<u8>> = Vec::with_capacity(cap);
    let mut pq: Vec<Vec<u8>> = Vec::new();
    let mut ql: Vec<usize> = Vec::new();
    let mut done = 0usize;
    while done < seqs.len() {
        let batch_n = (seqs.len() - done).min(FORMAT_BATCH_SIZE);
        // Header batch: a slice of the pre-decoded ids, or a fresh batch pulled from
        // the streaming cursor.
        let header_batch: &[Vec<u8>] = if let Some(ids) = decoded_headers {
            &ids[done..done + batch_n]
        } else {
            let c = h_cur.as_deref_mut().expect("cursor when not pre-decoded");
            cursor_headers.clear();
            for _ in 0..batch_n {
                cursor_headers.push(c.read_one_record_varint()?);
            }
            &cursor_headers
        };
        pq.clear();
        ql.clear();
        if let Some(q) = q_cur.as_deref_mut() {
            for k in 0..batch_n {
                let seq_len = seqs[done + k].len();
                drain_one_qual(
                    q,
                    const_qual_len,
                    bits_per_qual,
                    seq_len,
                    record_base + (done + k) as u64,
                    &mut pq,
                    &mut ql,
                )?;
            }
        }
        let qi = build_quality_inputs(q_cur.is_some(), is_fqz, batch_n, &pq, &ql, quality_binning);
        let assembled = format_batch_parallel(
            header_batch,
            &seqs[done..done + batch_n],
            &no_rc,
            &qi,
            is_fasta,
            est_per_record,
        )?;
        for chunk_out in &assembled {
            if push_and_flush(out, buf, chunk_out)? {
                return Ok(()); // downstream pipe closed
            }
        }
        done += batch_n;
    }
    Ok(())
}

/// Output target for `decode_reference_streaming`: `Two` writes R1→`out1`,
/// R2→`out2` (the default two-file / two-region paired-reference decode);
/// `Interleaved` merges both mates into ONE sink (R1/R2 records alternating) for
/// `qz decompress --interleaved`.
enum RefSink<'a> {
    Two(&'a mut dyn std::io::Write, &'a mut dyn std::io::Write),
    Interleaved(&'a mut dyn std::io::Write),
}

/// Interleaved twin of [`emit_chunk_mate`]: emit one chunk's R1[k] then R2[k]
/// alternately to ONE sink. Both mates' sequences are already reconstructed; this
/// pulls each mate's headers + qualities from its own cursor (or pre-decoded delta
/// ids), builds a 2·batch_n interleaved batch, and formats once via
/// `format_batch_parallel`. Bounded-RAM (one chunk's two mates resident).
#[allow(clippy::too_many_arguments)]
fn emit_chunk_pair_interleaved(
    r1_seqs: &[Vec<u8>],
    mut h1_cur: Option<&mut StreamCursor>,
    r1_decoded: Option<&[Vec<u8>]>,
    mut q1_cur: Option<&mut StreamCursor>,
    is_fqz_q1: bool,
    r2_seqs: &[Vec<u8>],
    mut h2_cur: Option<&mut StreamCursor>,
    r2_decoded: Option<&[Vec<u8>]>,
    mut q2_cur: Option<&mut StreamCursor>,
    is_fqz_q2: bool,
    const_qual_len: usize,
    bits_per_qual: usize,
    quality_binning: QualityBinning,
    is_fasta: bool,
    out: &mut dyn std::io::Write,
    buf: &mut Vec<u8>,
    record_base: u64,
) -> Result<()> {
    if r1_seqs.len() != r2_seqs.len() {
        bail!(
            "reference interleaved: chunk R1/R2 read count mismatch ({} vs {})",
            r1_seqs.len(),
            r2_seqs.len()
        );
    }
    let no_rc: Vec<u8> = Vec::new(); // reference applies RC inside reconstruct_read
    let est_per_record = 80 + 150 + 4 + 150 + 1;
    let n = r1_seqs.len();
    let mut h1b: Vec<Vec<u8>> = Vec::new();
    let mut h2b: Vec<Vec<u8>> = Vec::new();
    let mut pq1: Vec<Vec<u8>> = Vec::new();
    let mut ql1: Vec<usize> = Vec::new();
    let mut pq2: Vec<Vec<u8>> = Vec::new();
    let mut ql2: Vec<usize> = Vec::new();
    let mut done = 0usize;
    while done < n {
        let batch_n = (n - done).min(FORMAT_BATCH_SIZE);
        // Header batches (pre-decoded delta ids, or pulled from the streaming cursor).
        let h1_batch: &[Vec<u8>] = if let Some(ids) = r1_decoded {
            &ids[done..done + batch_n]
        } else {
            let c = h1_cur.as_deref_mut().expect("r1 cursor when not pre-decoded");
            h1b.clear();
            for _ in 0..batch_n {
                h1b.push(c.read_one_record_varint()?);
            }
            &h1b
        };
        let h2_batch: &[Vec<u8>] = if let Some(ids) = r2_decoded {
            &ids[done..done + batch_n]
        } else {
            let c = h2_cur.as_deref_mut().expect("r2 cursor when not pre-decoded");
            h2b.clear();
            for _ in 0..batch_n {
                h2b.push(c.read_one_record_varint()?);
            }
            &h2b
        };
        // Quality batches.
        pq1.clear();
        ql1.clear();
        if let Some(q) = q1_cur.as_deref_mut() {
            for k in 0..batch_n {
                drain_one_qual(
                    q, const_qual_len, bits_per_qual, r1_seqs[done + k].len(),
                    record_base + (done + k) as u64, &mut pq1, &mut ql1,
                )?;
            }
        }
        let qi1 = build_quality_inputs(q1_cur.is_some(), is_fqz_q1, batch_n, &pq1, &ql1, quality_binning);
        pq2.clear();
        ql2.clear();
        if let Some(q) = q2_cur.as_deref_mut() {
            for k in 0..batch_n {
                drain_one_qual(
                    q, const_qual_len, bits_per_qual, r2_seqs[done + k].len(),
                    record_base + (done + k) as u64, &mut pq2, &mut ql2,
                )?;
            }
        }
        let qi2 = build_quality_inputs(q2_cur.is_some(), is_fqz_q2, batch_n, &pq2, &ql2, quality_binning);

        // Interleave R1[k], R2[k] → one 2·batch_n batch, format once.
        let mut hi: Vec<Vec<u8>> = Vec::with_capacity(2 * batch_n);
        let mut si: Vec<Vec<u8>> = Vec::with_capacity(2 * batch_n);
        let mut qii: Vec<_> = Vec::with_capacity(2 * batch_n);
        let mut q1it = qi1.into_iter();
        let mut q2it = qi2.into_iter();
        for k in 0..batch_n {
            hi.push(h1_batch[k].clone());
            si.push(r1_seqs[done + k].clone());
            qii.push(q1it.next().unwrap());
            hi.push(h2_batch[k].clone());
            si.push(r2_seqs[done + k].clone());
            qii.push(q2it.next().unwrap());
        }
        let assembled = format_batch_parallel(&hi, &si, &no_rc, &qii, is_fasta, est_per_record)?;
        for chunk_out in &assembled {
            if push_and_flush(out, buf, chunk_out)? {
                return Ok(()); // downstream pipe closed
            }
        }
        done += batch_n;
    }
    Ok(())
}

/// Decode one chunk's R1 (always columnar) and R2 (columnar OR delta-vs-R1) header
/// ids, for archives that use the R2 `HeaderDelta` representation. Mirrors the paired
/// per-chunk reconstruction (`paired::mod`): bounded-RAM (one chunk's ids resident),
/// positioned on-demand reads. Returns `(r1_ids, r2_ids)`, both `record_count`-long.
fn decode_chunk_headers_delta(
    reader: &RefReader,
    dir: &ChunkDirectory,
    chunk: u32,
) -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>)> {
    let find = |mate: u8, role: StreamRole| {
        dir.entries
            .iter()
            .find(|e| e.chunk_index == chunk && e.mate == mate && e.role == role)
    };
    let e_r1 = find(1, StreamRole::Headers)
        .ok_or_else(|| anyhow::anyhow!("chunk {chunk}: no R1 headers"))?;
    let r1_ids = crate::compression::paired::decode_columnar_headers(
        &reader.read_entry(e_r1.offset, e_r1.length)?,
        e_r1.record_count,
    )?;
    // R2: independent columnar (role Headers) OR delta-vs-R1 (role HeaderDelta).
    let r2_ids = if let Some(e) = find(2, StreamRole::Headers) {
        crate::compression::paired::decode_columnar_headers(
            &reader.read_entry(e.offset, e.length)?,
            e.record_count,
        )?
    } else if let Some(e) = find(2, StreamRole::HeaderDelta) {
        let ops = crate::compression::paired::streams::decode_bsc_stream(
            &reader.read_entry(e.offset, e.length)?,
            crate::compression::codecs::stream_decode_cap(e.record_count as usize),
        )?;
        let ids = crate::compression::paired::header_delta::decode(
            &ops,
            &r1_ids,
            e.record_count as usize,
        )?;
        if ids.len() as u64 != e.record_count {
            bail!("chunk {chunk}: R2 delta decoded {} ids != {}", ids.len(), e.record_count);
        }
        ids
    } else {
        bail!("chunk {chunk}: no R2 header role");
    };
    Ok((r1_ids, r2_ids))
}

/// Engine-aligned reference decode (spec §4): converge headers + qualities onto the
/// shared single-end block-streaming engine while keeping the chunk-granular
/// sequence reconstruction as the one justified reference-specific divergence.
///
/// Both mates' header + quality streams are decoded by the shared block-parallel
/// producers (`spawn_stream_producer`) feeding `StreamCursor`s — continuous across
/// the whole archive. Sequences are reconstructed chunk-by-chunk (`reconstruct_mate`,
/// par_iter within a chunk: reads are independent — mate-1 positions are absolute,
/// mate-2 deltas reference the same pair's R1 anchor, so there is no cross-read
/// prefix-sum) and spliced with the streamed headers/qualities via
/// `format_batch_parallel`.
///
/// Byte-identical to a sequential decode: chunks are walked `0..num_chunks` and the
/// flat header/qual cursors yield records in the same chunk order (the per-mate
/// builders sort their entries by chunk_index), so header/qual record i lines up with
/// reconstructed sequence i.
///
/// Bounded-RAM (constant in read count): the resident globals + ONE chunk's
/// reconstructed seqs + the producer pipelines' bounded channels. `reader` is used
/// only on this (the orchestrating) thread for the globals + per-chunk column reads;
/// the producers open their own file handles via `archive_path`, so nothing shares
/// `reader` and no `Arc` is needed. `w1`/`w2` need not be `Send`.
fn decode_reference_streaming(
    archive_path: &std::path::Path,
    reader: RefReader,
    dir: &ChunkDirectory,
    hdr: &FixedHeader,
    sink: &mut RefSink<'_>,
    chunk_range: Option<(u32, u32)>,
) -> Result<()> {
    use std::sync::mpsc::sync_channel;
    // Mirrors the engine's BOUNDED_CHANNEL_CAP (one block in the cursor + up to this
    // many buffered) — backpressure tuning, not correctness.
    const CHAN_CAP: usize = 2;

    let (max_block_len, max_inflight) = decode_block_bounds(hdr.encoding_type);
    let quality_binning = code_to_binning(hdr.quality_binning_code)?;
    // Range-filtered per-mate header/quality indices: the producers start at the
    // FIRST chunk of the range, so the cursors align with the range's
    // chunk-by-chunk reconstructed sequences (each chunk's seq columns are read
    // positionally by `reconstruct_mate`; no cross-chunk pre-skip is needed).
    // `None` = full archive, byte-identical to the original calls.
    let r1 = build_ref_mate_streams(archive_path, dir, 1, chunk_range)?;
    let r2 = build_ref_mate_streams(archive_path, dir, 2, chunk_range)?;

    // R2 headers may be delta-vs-R1 (role HeaderDelta) instead of independent columnar.
    // When so, R2's header stream isn't self-contained, so we skip the streaming header
    // producers and reconstruct BOTH mates' ids per chunk (R1 columnar → ids; R2 delta
    // applied against them — `decode_chunk_headers_delta`), bounded by one chunk's ids.
    // Columnar-only archives (incl. all pre-existing ones) keep the streaming path.
    let r2_uses_delta = dir
        .entries
        .iter()
        .any(|e| e.mate == 2 && e.role == StreamRole::HeaderDelta);

    // Fqzcomp quality decodes to FINAL raw 8-bit ASCII (varint-framed [len][bytes]);
    // its per-record reader must use bits_per_qual = 8 so packed_len == len. Only
    // BSC-packed quality uses the binning's bits-per-quality. Quality codec is uniform
    // per archive (both mates the same), so one bits_per_qual covers both.
    let is_fqz_q1 = matches!(r1.quality_codec, Some(StreamCodec::Fqz));
    let is_fqz_q2 = matches!(r2.quality_codec, Some(StreamCodec::Fqz));
    let bits_per_qual = if is_fqz_q1 || is_fqz_q2 {
        8
    } else {
        quality_binning.bits_per_quality()
    };
    let const_qual_len = hdr.const_qual_len as usize; // reference reads vary ⇒ 0
    let is_fasta = hdr.is_fasta();
    // Walk only the requested chunk range; `None` = the whole archive (byte-identical).
    let (chunk_lo, chunk_hi) = chunk_range.unwrap_or((0, dir.num_chunks));

    // Globals (backing + interval map + reference meta) loaded once, resident.
    let globals = decode_globals(&reader, dir)?;

    let (h1_tx, h1_rx) = sync_channel::<Result<(u32, Vec<u8>), String>>(CHAN_CAP);
    let (h2_tx, h2_rx) = sync_channel::<Result<(u32, Vec<u8>), String>>(CHAN_CAP);
    let q1_pair = r1
        .qualities
        .as_ref()
        .map(|_| sync_channel::<Result<(u32, Vec<u8>), String>>(CHAN_CAP));
    let q2_pair = r2
        .qualities
        .as_ref()
        .map(|_| sync_channel::<Result<(u32, Vec<u8>), String>>(CHAN_CAP));

    std::thread::scope(|scope| -> Result<()> {
        // Header producers (both mates independent columnar streams). Skipped for
        // R2-HeaderDelta archives, which reconstruct headers per chunk instead.
        if !r2_uses_delta {
        spawn_stream_producer(
            scope,
            archive_path,
            &r1.headers,
            StreamCodec::Columnar,
            h1_tx,
            max_inflight,
            max_block_len,
        );
        spawn_stream_producer(
            scope,
            archive_path,
            &r2.headers,
            StreamCodec::Columnar,
            h2_tx,
            max_inflight,
            max_block_len,
        );
        } // end if !r2_uses_delta (header producers)
        // Quality producers (fqzcomp → decoded ASCII; BSC → bit-packed).
        if let (Some(idx), Some((tx, _))) = (r1.qualities.as_ref(), q1_pair.as_ref()) {
            spawn_stream_producer(
                scope,
                archive_path,
                idx,
                if is_fqz_q1 {
                    StreamCodec::Fqz
                } else {
                    StreamCodec::Bsc
                },
                tx.clone(),
                max_inflight,
                max_block_len,
            );
        }
        if let (Some(idx), Some((tx, _))) = (r2.qualities.as_ref(), q2_pair.as_ref()) {
            spawn_stream_producer(
                scope,
                archive_path,
                idx,
                if is_fqz_q2 {
                    StreamCodec::Fqz
                } else {
                    StreamCodec::Bsc
                },
                tx.clone(),
                max_inflight,
                max_block_len,
            );
        }

        // Header cursors exist only for the streaming (columnar) path; the delta path
        // reconstructs ids per chunk instead (so h1_rx/h2_rx go unused and disconnect).
        let mut h1_cur = (!r2_uses_delta).then(move || StreamCursor::new(h1_rx, "r1_headers"));
        let mut h2_cur = (!r2_uses_delta).then(move || StreamCursor::new(h2_rx, "r2_headers"));
        let mut q1_cur = q1_pair.map(|(_, rx)| StreamCursor::new(rx, "r1_qualities"));
        let mut q2_cur = q2_pair.map(|(_, rx)| StreamCursor::new(rx, "r2_qualities"));

        let mut buf1: Vec<u8> = Vec::with_capacity(WRITE_BATCH + 512);
        let mut buf2: Vec<u8> = Vec::with_capacity(WRITE_BATCH + 512);
        // Self-contained part: numbering starts at 0 (the parent concatenates parts
        // in chunk order, so part-local numbering matches a full decode's chunk-0
        // start). `record_base` is only an error-message/progress counter passed to
        // `emit_chunk_mate`/`drain_one_qual` — never an index into a whole-archive
        // structure — so a range starting mid-archive is safe (verified).
        let mut record_base = 0u64;

        for c in chunk_lo..chunk_hi {
            // Headers: streaming columnar cursors, OR (R2-HeaderDelta archives) both
            // mates' ids reconstructed per chunk (R1 columnar → ids; R2 delta vs them).
            let delta_hdrs = if r2_uses_delta {
                Some(decode_chunk_headers_delta(&reader, dir, c)?)
            } else {
                None
            };
            let (r1_decoded, r2_decoded) = match &delta_hdrs {
                Some((a, b)) => (Some(a.as_slice()), Some(b.as_slice())),
                None => (None, None),
            };

            // Reconstruct BOTH mates before emitting: Mate 1 (absolute positions,
            // capturing per-read anchors) then Mate 2 (positions delta/escape vs the
            // SAME pair's R1 anchor). The two-sink path emits R1→w1 then R2→w2 (bytes
            // unchanged from before); the interleaved path needs both mates resident to
            // alternate records.
            let (r1_seqs, r1_anchors) = reconstruct_mate(&reader, dir, c, 1, &globals, None)?;
            let (r2_seqs, _r2_anchors) =
                reconstruct_mate(&reader, dir, c, 2, &globals, Some(&r1_anchors))?;
            if r2_seqs.len() != r1_seqs.len() {
                bail!(
                    "chunk {c}: R1/R2 read count mismatch ({} vs {})",
                    r1_seqs.len(),
                    r2_seqs.len()
                );
            }
            match sink {
                RefSink::Two(w1, w2) => {
                    emit_chunk_mate(
                        &r1_seqs, h1_cur.as_mut(), r1_decoded, q1_cur.as_mut(), is_fqz_q1,
                        const_qual_len, bits_per_qual, quality_binning, is_fasta,
                        &mut **w1, &mut buf1, record_base,
                    )?;
                    emit_chunk_mate(
                        &r2_seqs, h2_cur.as_mut(), r2_decoded, q2_cur.as_mut(), is_fqz_q2,
                        const_qual_len, bits_per_qual, quality_binning, is_fasta,
                        &mut **w2, &mut buf2, record_base,
                    )?;
                }
                RefSink::Interleaved(w) => {
                    emit_chunk_pair_interleaved(
                        &r1_seqs, h1_cur.as_mut(), r1_decoded, q1_cur.as_mut(), is_fqz_q1,
                        &r2_seqs, h2_cur.as_mut(), r2_decoded, q2_cur.as_mut(), is_fqz_q2,
                        const_qual_len, bits_per_qual, quality_binning, is_fasta,
                        &mut **w, &mut buf1, record_base,
                    )?;
                }
            }
            record_base += r1_seqs.len() as u64;
        }

        match sink {
            RefSink::Two(w1, w2) => {
                flush_remaining(&mut **w1, &buf1)?;
                flush_remaining(&mut **w2, &buf2)
            }
            RefSink::Interleaved(w) => flush_remaining(&mut **w, &buf1),
        }
    })
}

/// Single-end engine-aligned reference decode: one mate (mate 1, absolute
/// positions, never anchors). Headers + qualities stream through the shared
/// block-parallel producers; sequences are reconstructed chunk-by-chunk via
/// `reconstruct_mate(.., 1, .., None)`. Byte-identical to a sequential decode
/// (chunks walked `chunk_lo..chunk_hi`, the flat header/qual cursor record i
/// lines up with reconstructed sequence i). Bounded-RAM: globals + one chunk's
/// reconstructed seqs + the producers' bounded channels.
fn decode_reference_single_streaming(
    archive_path: &std::path::Path,
    reader: RefReader,
    dir: &ChunkDirectory,
    hdr: &FixedHeader,
    w: &mut dyn std::io::Write,
    chunk_range: Option<(u32, u32)>,
) -> Result<()> {
    use std::sync::mpsc::sync_channel;
    const CHAN_CAP: usize = 2;

    let (max_block_len, max_inflight) = decode_block_bounds(hdr.encoding_type);
    let quality_binning = code_to_binning(hdr.quality_binning_code)?;
    // Single-end reference is always columnar R1Headers (no HeaderDelta).
    let r1 = build_ref_mate_streams(archive_path, dir, 1, chunk_range)?;

    let is_fqz_q = matches!(r1.quality_codec, Some(StreamCodec::Fqz));
    let bits_per_qual = if is_fqz_q {
        8
    } else {
        quality_binning.bits_per_quality()
    };
    let const_qual_len = hdr.const_qual_len as usize; // reference reads vary ⇒ 0
    let is_fasta = hdr.is_fasta();
    let (chunk_lo, chunk_hi) = chunk_range.unwrap_or((0, dir.num_chunks));

    let globals = decode_globals(&reader, dir)?;

    let (h_tx, h_rx) = sync_channel::<Result<(u32, Vec<u8>), String>>(CHAN_CAP);
    let q_pair = r1
        .qualities
        .as_ref()
        .map(|_| sync_channel::<Result<(u32, Vec<u8>), String>>(CHAN_CAP));

    std::thread::scope(|scope| -> Result<()> {
        spawn_stream_producer(
            scope,
            archive_path,
            &r1.headers,
            StreamCodec::Columnar,
            h_tx,
            max_inflight,
            max_block_len,
        );
        if let (Some(idx), Some((tx, _))) = (r1.qualities.as_ref(), q_pair.as_ref()) {
            spawn_stream_producer(
                scope,
                archive_path,
                idx,
                if is_fqz_q {
                    StreamCodec::Fqz
                } else {
                    StreamCodec::Bsc
                },
                tx.clone(),
                max_inflight,
                max_block_len,
            );
        }

        let mut h_cur = StreamCursor::new(h_rx, "headers");
        let mut q_cur = q_pair.map(|(_, rx)| StreamCursor::new(rx, "qualities"));

        let mut buf: Vec<u8> = Vec::with_capacity(WRITE_BATCH + 512);
        let mut record_base = 0u64;

        for c in chunk_lo..chunk_hi {
            // Mate 1: absolute positions, no anchors (r1_anchor_per_read = None).
            let (seqs, _anchors) = reconstruct_mate(&reader, dir, c, 1, &globals, None)?;
            emit_chunk_mate(
                &seqs,
                Some(&mut h_cur),
                None,
                q_cur.as_mut(),
                is_fqz_q,
                const_qual_len,
                bits_per_qual,
                quality_binning,
                is_fasta,
                w,
                &mut buf,
                record_base,
            )?;
            record_base += seqs.len() as u64;
        }
        flush_remaining(w, &buf)
    })
}

/// Open a reference archive for **seekable** decode: open the file, read only the
/// 64-byte front-header prefix (asserting v5 `archive_type == expect_type`: 2 for
/// paired, 4 for single-end), read + validate the directory footer FROM THE FILE
/// via `read_v5_footer`, run the matching reference semantic validator (single
/// vs paired), and return `(RefReader, directory, header)`. The whole archive is
/// NOT read into memory — callers (`decompress_reference{,_single}`, both verify
/// paths) read the small globals up front and each per-chunk entry on demand via
/// `RefReader::read_entry`. The file length is captured ONCE here (into
/// `RefReader::len`) so per-entry reads bounds-check against it without an `fstat`
/// per entry. Mirrors the paired file-based open (`decompress_paired_v5`), reusing
/// the shared `read_v5_footer`.
fn open_ref_typed(
    input: &std::path::Path,
    expect_type: u8,
) -> Result<(RefReader, ChunkDirectory, FixedHeader)> {
    use std::io::Read;
    // Front header → header_end (footer scan anchor) + archive_type check. Read
    // only the fixed 64-byte prefix; the footer reader seeks to EOF for the rest.
    let mut hbuf = [0u8; 64];
    let mut file = std::fs::File::open(input)
        .map_err(|e| anyhow::anyhow!("open reference archive {}: {e}", input.display()))?;
    let n = file.read(&mut hbuf)?;
    if n < crate::compression::V2_PREFIX_SIZE {
        bail!("reference v5 header truncated");
    }
    let hdr = crate::compression::FixedHeader::parse_v5(&hbuf[..n])?;
    if hdr.archive_type != expect_type {
        bail!(
            "not a reference v5 archive of expected type {expect_type} (archive_type {})",
            hdr.archive_type
        );
    }
    let header_end = u32::from_le_bytes(hbuf[4..8].try_into().unwrap()) as u64;
    let dir = crate::compression::chunk_directory::read_v5_footer(input, header_end)?;
    if expect_type == 4 {
        validate_reference_directory_single(&dir)?;
    } else {
        validate_reference_directory(&dir)?;
    }
    // Capture the file length ONCE. The footer validator already proved every
    // entry's span lies within `[header_end, footer_start) <= len`, so per-entry
    // reads need no fresh `fstat` — they reuse this `len` for a defense-in-depth
    // bounds check. The header read above advanced the cursor, but `read_entry`
    // only does positioned reads, so the cursor position is irrelevant.
    let len = file
        .metadata()
        .map_err(|e| anyhow::anyhow!("stat reference archive {}: {e}", input.display()))?
        .len();
    Ok((RefReader { file, len }, dir, hdr))
}

/// Open a PAIRED reference archive (archive_type 2) for seekable decode.
fn open_ref(input: &std::path::Path) -> Result<(RefReader, ChunkDirectory, FixedHeader)> {
    open_ref_typed(input, 2)
}

/// Open a SINGLE-END reference archive (archive_type 4) for seekable decode.
fn open_ref_single(input: &std::path::Path) -> Result<(RefReader, ChunkDirectory, FixedHeader)> {
    open_ref_typed(input, 4)
}

/// Verify a reference archive.
///
/// **Deep** (`fast == false`): full reconstruct via
/// `decode_reference_streaming` (which already enforces every content
/// invariant — edit.pos < read_len, span ⊂ interval, count relations, R2 delta
/// bounds — and errors otherwise), driven through two per-mate `CrcWriter` sinks.
/// R1 and R2 each get their own CRC32 over their reconstructed FASTQ bytes
/// (reported as `crc32` and `r2_crc32`).
///
/// **Fast** (`fast == true`): per directory entry, CRC-verify every v4 block
/// WITHOUT bsc-decoding, dispatching by block-framing style: Consensus / NBitmap
/// use the `[num_blocks u32][v4...]` framing (verified via
/// `backing::verify_block_stream_crc`); every other entry uses the
/// `streams::write_block_stream` framing (verified via
/// `streams::decode_block_payloads`). In both cases Σ per-block record_count is
/// cross-checked against the entry's declared `record_count`.
/// Cosmetic `encoding_type` reported in a reference archive's `VerifyResult`
/// (reference archives are identified on the wire by v5 `archive_type == 2`, not
/// by this byte; this is purely a display tag for `qz verify`).
const REFERENCE_REPORTED_ENCODING: u8 = 0xF2;

pub(crate) fn verify_reference(
    input: &std::path::Path,
    fast: bool,
) -> Result<crate::compression::VerifyResult> {
    use crate::cli::{HeaderCompressor, QualityCompressor};
    use crate::compression::paired::CrcWriter;
    use crate::compression::paired::streams;
    use crate::compression::{VerifyMode, VerifyResult};

    let started = std::time::Instant::now();
    let (reader, dir, hdr) = open_ref(input)?;
    let num_reads = dir.num_reads.saturating_mul(2) as usize; // 2 mates per pair

    if fast {
        let mut blocks_verified: u32 = 0;
        let mut total_bytes: u64 = 0;
        // Cross-check that NBitmap and Consensus declare the same base count
        // (already enforced by the footer validator, but re-asserted here).
        let mut consensus_bases: Option<u64> = None;
        let mut nbitmap_bases: Option<u64> = None;

        for e in &dir.entries {
            // Read each entry's payload on demand (positioned read) — no
            // whole-archive buffer. `e_bytes` (owned Vec) derefs to `&[u8]` at
            // every consumer below.
            let e_bytes = reader.read_entry(e.offset, e.length)?;
            let (blocks, rc_sum) = match e.role {
                StreamRole::PackedBacking | StreamRole::NBitmap => {
                    backing::verify_block_stream_crc(&e_bytes)?
                }
                StreamRole::ChunkDecodedSizes => {
                    // Cross-cutting global written via the shared `write_role_blocks` as ONE
                    // v5 block frame ([block_len][record_count][crc32][payload]) — NOT the
                    // reference block-stream the `_` arm reads. Verify that single frame's CRC
                    // (mirrors `read_chunk_layout`'s parse). Counts as 1 block, no records.
                    if e_bytes.len() < 12 {
                        anyhow::bail!("ChunkDecodedSizes entry shorter than 12-byte frame header");
                    }
                    let block_len = u32::from_le_bytes(e_bytes[0..4].try_into().unwrap()) as usize;
                    let record_count = u32::from_le_bytes(e_bytes[4..8].try_into().unwrap());
                    let crc = u32::from_le_bytes(e_bytes[8..12].try_into().unwrap());
                    let payload = &e_bytes[12..];
                    if payload.len() != block_len {
                        anyhow::bail!(
                            "ChunkDecodedSizes frame: block_len {block_len} != payload {}",
                            payload.len()
                        );
                    }
                    crate::compression::bsc::verify_block_frame_crc(crc, record_count, payload)?;
                    // The frame's record_count == the directory entry's record_count
                    // (both = num_chunks, set together by `write_role_blocks`); return it
                    // so the `rc_sum == e.record_count` cross-check below holds.
                    (1u32, record_count as u64)
                }
                StreamRole::ReferenceMeta => {
                    // ReferenceMeta = [meta_version:1B][digest:8B] then a record
                    // block-stream; CRC-walk only the block-stream tail.
                    let stream = e_bytes.get(9..).ok_or_else(|| {
                        anyhow::anyhow!("ReferenceMeta entry shorter than 9-byte prefix")
                    })?;
                    let pls = streams::decode_block_payloads(stream)?;
                    let rc_sum: u64 = pls.iter().map(|(rc, _)| *rc as u64).sum();
                    (pls.len() as u32, rc_sum)
                }
                _ => {
                    let pls = streams::decode_block_payloads(&e_bytes)?; // CRC-verifies every block
                    let rc_sum: u64 = pls.iter().map(|(rc, _)| *rc as u64).sum();
                    (pls.len() as u32, rc_sum)
                }
            };
            if rc_sum != e.record_count {
                bail!(
                    "entry (chunk {}, mate {}, role {:?}) block record_count sum {} != {}",
                    e.chunk_index,
                    e.mate,
                    e.role,
                    rc_sum,
                    e.record_count
                );
            }
            match e.role {
                StreamRole::PackedBacking => consensus_bases = Some(rc_sum),
                StreamRole::NBitmap => nbitmap_bases = Some(rc_sum),
                _ => {}
            }
            blocks_verified += blocks;
            total_bytes += e_bytes.len() as u64;
        }

        if consensus_bases != nbitmap_bases {
            bail!(
                "N-bitmap base count {:?} != consensus base count {:?}",
                nbitmap_bases,
                consensus_bases
            );
        }

        return Ok(VerifyResult {
            num_reads,
            encoding_type: REFERENCE_REPORTED_ENCODING,
            header_compressor: HeaderCompressor::Columnar,
            quality_compressor: QualityCompressor::Auto,
            headers_compressed_len: 0,
            sequences_compressed_len: 0,
            qualities_compressed_len: 0,
            crc32: 0,
            // Fast verify is a per-block CRC walk that does not separate the R1/R2
            // reconstruction, so there is no per-mate CRC to report here. Per-mate
            // CRCs are a deep-verify feature.
            r2_crc32: None,
            total_bytes,
            blocks_verified,
            mode: VerifyMode::Fast,
            elapsed_secs: started.elapsed().as_secs_f64(),
        });
    }

    // ---- Deep: full reconstruct + per-mate CRC over the FASTQ bytes ----
    // Two CrcWriter sinks (one per mate) drive the reconstruction without ever
    // materialising a full-FASTQ Vec — peak verify RAM is bounded by one chunk's
    // working set plus the globals, constant in decoded output size. R1 and R2
    // get separate CRC32s (no crc32_combine exists, so we report both): `crc32`
    // carries the R1 CRC and `r2_crc32` the R2 CRC. `total_bytes` is the sum.
    let mut s1 = CrcWriter {
        crc: flate2::Crc::new(),
        bytes: 0,
    };
    let mut s2 = CrcWriter {
        crc: flate2::Crc::new(),
        bytes: 0,
    };
    // Same engine-aligned decode the production path uses, driving the two per-mate
    // CrcWriter sinks. The streamed byte sequence is exactly what a file decode would
    // write, so the reported per-mate CRCs are the production output's CRCs.
    decode_reference_streaming(
        input,
        reader,
        &dir,
        &hdr,
        &mut RefSink::Two(&mut s1, &mut s2),
        None,
    )?;
    let r1_crc = s1.crc.sum();
    let r2_crc = s2.crc.sum();
    let total_bytes = s1.bytes + s2.bytes;

    Ok(VerifyResult {
        num_reads,
        encoding_type: REFERENCE_REPORTED_ENCODING,
        header_compressor: HeaderCompressor::Columnar,
        quality_compressor: QualityCompressor::Auto,
        headers_compressed_len: 0,
        sequences_compressed_len: 0,
        qualities_compressed_len: 0,
        crc32: r1_crc,
        r2_crc32: Some(r2_crc),
        total_bytes,
        blocks_verified: 0,
        mode: VerifyMode::Deep,
        elapsed_secs: started.elapsed().as_secs_f64(),
    })
}

/// Verify a **single-end** reference archive (`archive_type == 4`).
///
/// **Deep** (`fast == false`): full reconstruct via
/// `decode_reference_single_streaming` into ONE `CrcWriter` sink — peak RAM is
/// bounded by one chunk's working set plus the globals, constant in output size.
/// All FASTQ bytes fold into the single `crc32`; `r2_crc32` is `None` (no mate 2).
///
/// **Fast** (`fast == true`): a per-block CRC32 walk with NO reconstruction
/// (mirrors `verify_reference`'s fast branch) — reports `blocks_verified` with
/// `crc32 == 0` and `r2_crc32 == None`.
pub(crate) fn verify_reference_single(
    input: &std::path::Path,
    fast: bool,
) -> Result<crate::compression::VerifyResult> {
    use crate::cli::{HeaderCompressor, QualityCompressor};
    use crate::compression::paired::streams;
    use crate::compression::paired::CrcWriter;
    use crate::compression::{VerifyMode, VerifyResult};

    let started = std::time::Instant::now();
    let (reader, dir, hdr) = open_ref_single(input)?;
    // Single-end: one read per record, NOT doubled (paired uses `* 2`).
    let num_reads = dir.num_reads as usize;

    if fast {
        let mut blocks_verified: u32 = 0;
        let mut total_bytes: u64 = 0;
        // Cross-check that NBitmap and Consensus declare the same base count
        // (already enforced by the footer validator, but re-asserted here).
        let mut consensus_bases: Option<u64> = None;
        let mut nbitmap_bases: Option<u64> = None;

        for e in &dir.entries {
            // Positioned read of each entry's payload — no whole-archive buffer.
            let e_bytes = reader.read_entry(e.offset, e.length)?;
            let (blocks, rc_sum) = match e.role {
                StreamRole::PackedBacking | StreamRole::NBitmap => {
                    backing::verify_block_stream_crc(&e_bytes)?
                }
                StreamRole::ChunkDecodedSizes => {
                    // Cross-cutting global written via the shared `write_role_blocks` as ONE
                    // v5 block frame ([block_len][record_count][crc32][payload]) — NOT the
                    // reference block-stream the `_` arm reads. Verify that single frame's CRC
                    // (mirrors `read_chunk_layout`'s parse). Counts as 1 block, no records.
                    if e_bytes.len() < 12 {
                        anyhow::bail!("ChunkDecodedSizes entry shorter than 12-byte frame header");
                    }
                    let block_len = u32::from_le_bytes(e_bytes[0..4].try_into().unwrap()) as usize;
                    let record_count = u32::from_le_bytes(e_bytes[4..8].try_into().unwrap());
                    let crc = u32::from_le_bytes(e_bytes[8..12].try_into().unwrap());
                    let payload = &e_bytes[12..];
                    if payload.len() != block_len {
                        anyhow::bail!(
                            "ChunkDecodedSizes frame: block_len {block_len} != payload {}",
                            payload.len()
                        );
                    }
                    crate::compression::bsc::verify_block_frame_crc(crc, record_count, payload)?;
                    // The frame's record_count == the directory entry's record_count
                    // (both = num_chunks, set together by `write_role_blocks`); return it
                    // so the `rc_sum == e.record_count` cross-check below holds.
                    (1u32, record_count as u64)
                }
                StreamRole::ReferenceMeta => {
                    // [meta_version:1B][digest:8B] then a record block-stream;
                    // CRC-walk only the block-stream tail.
                    let stream = e_bytes.get(9..).ok_or_else(|| {
                        anyhow::anyhow!("ReferenceMeta entry shorter than 9-byte prefix")
                    })?;
                    let pls = streams::decode_block_payloads(stream)?;
                    let rc_sum: u64 = pls.iter().map(|(rc, _)| *rc as u64).sum();
                    (pls.len() as u32, rc_sum)
                }
                _ => {
                    let pls = streams::decode_block_payloads(&e_bytes)?; // CRC-verifies every block
                    let rc_sum: u64 = pls.iter().map(|(rc, _)| *rc as u64).sum();
                    (pls.len() as u32, rc_sum)
                }
            };
            if rc_sum != e.record_count {
                bail!(
                    "entry (chunk {}, mate {}, role {:?}) block record_count sum {} != {}",
                    e.chunk_index,
                    e.mate,
                    e.role,
                    rc_sum,
                    e.record_count
                );
            }
            match e.role {
                StreamRole::PackedBacking => consensus_bases = Some(rc_sum),
                StreamRole::NBitmap => nbitmap_bases = Some(rc_sum),
                _ => {}
            }
            blocks_verified += blocks;
            total_bytes += e_bytes.len() as u64;
        }

        if consensus_bases != nbitmap_bases {
            bail!(
                "N-bitmap base count {:?} != consensus base count {:?}",
                nbitmap_bases,
                consensus_bases
            );
        }

        return Ok(VerifyResult {
            num_reads,
            encoding_type: REFERENCE_REPORTED_ENCODING,
            header_compressor: HeaderCompressor::Columnar,
            quality_compressor: QualityCompressor::Auto,
            headers_compressed_len: 0,
            sequences_compressed_len: 0,
            qualities_compressed_len: 0,
            crc32: 0,
            r2_crc32: None,
            total_bytes,
            blocks_verified,
            mode: VerifyMode::Fast,
            elapsed_secs: started.elapsed().as_secs_f64(),
        });
    }

    // ---- Deep: full reconstruct + single CRC over the one FASTQ ----
    let mut sink = CrcWriter { crc: flate2::Crc::new(), bytes: 0 };
    decode_reference_single_streaming(input, reader, &dir, &hdr, &mut sink, None)?;
    let crc = sink.crc.sum();
    let total_bytes = sink.bytes;

    Ok(VerifyResult {
        num_reads,
        encoding_type: REFERENCE_REPORTED_ENCODING,
        header_compressor: HeaderCompressor::Columnar,
        quality_compressor: QualityCompressor::Auto,
        headers_compressed_len: 0,
        sequences_compressed_len: 0,
        qualities_compressed_len: 0,
        crc32: crc,
        r2_crc32: None, // single-end reference: one stream, no per-mate split
        total_bytes,
        blocks_verified: 0,
        mode: VerifyMode::Deep,
        elapsed_secs: started.elapsed().as_secs_f64(),
    })
}

/// Reconstruct one mate's per-read sequences for a chunk (spec §4.4) and return,
/// per read in original order, the read's mapped anchor `Some((ref_id, ref_pos))`
/// (`None` for fallback reads). RC is applied here (inside `reconstruct_read`) so the
/// returned sequences are final FASTQ bytes. For mate 2, `r1_anchor_per_read` is the
/// mate-1 group's per-read anchor vector (length N); the per-mapped-R2 anchor passed
/// to `decode_positions_mate2` is built by selecting the anchor at each mapped R2's
/// read index.
///
/// The per-read reconstruction runs in `par_iter` — reads are independent (mate-1
/// positions are absolute; mate-2 deltas reference the SAME pair's R1 anchor, so
/// there is no cross-read prefix-sum). The two sequential streams (the `ReadLen`
/// varint run and the `FallbackPool` reader) are drained into `Vec`s in serial
/// pre-passes first, so the parallel loop only does independent indexed work. All the
/// untrusted-input validation the old serial loop performed (popcount == M, ref_id <
/// num_refs, span ⊂ one interval, span within backing, edit.pos < read_len, exact
/// pool drain, exact ReadLen consumption, mapped/fallback counts) is preserved.
fn reconstruct_mate(
    reader: &RefReader,
    footer: &ChunkDirectory,
    chunk_index: u32,
    mate: u8,
    globals: &Globals,
    r1_anchor_per_read: Option<&[Option<(u32, u64)>]>,
) -> Result<(Vec<Vec<u8>>, Vec<Option<(u32, u64)>>)> {
    let Globals {
        backing,
        imap,
        meta,
    } = globals;
    let find = |role: StreamRole| -> Result<&ChunkDirEntry> {
        footer
            .entries
            .iter()
            .find(|e| e.chunk_index == chunk_index && e.mate == mate && e.role == role)
            .ok_or_else(|| {
                anyhow::anyhow!("chunk {chunk_index} mate {mate}: missing role {role:?}")
            })
    };

    // This (chunk,mate)'s fallback literals live in its own FallbackPool entry
    // (absent ⇒ this mate has no fallbacks; the validator guarantees that ⇒ no
    // read here is Fallback). Read its bytes on demand and decode block-lazily.
    // `pool_bytes` is kept in scope so it outlives the borrowing `pool` reader.
    let pool_entry = footer.entries.iter().find(|e| {
        e.chunk_index == chunk_index && e.mate == mate && e.role == StreamRole::FallbackPool
    });
    let pool_bytes = match pool_entry {
        Some(e) => Some(reader.read_entry(e.offset, e.length)?),
        None => None,
    };
    let mut pool = match pool_bytes.as_ref() {
        Some(b) => Some(fallback::FallbackPoolReader::new(b)?),
        None => None,
    };

    let e_flags = find(StreamRole::MappedFlags)?;
    // `record_count` is untrusted (footer); clamp to usize before it sizes
    // length-N flag/position vectors and the reconstruction loop bound.
    let n_chunk = usize::try_from(e_flags.record_count)?;
    let flags_raw = decode_bsc_role_file(reader, e_flags)?;
    if flags_raw.len() != n_chunk.div_ceil(8) {
        bail!("chunk {chunk_index} mate {mate}: mapped_flags wrong length");
    }
    // Validate unused high bits in the last flags byte are zero.
    if !n_chunk.is_multiple_of(8) {
        let mask = !((1u8 << (n_chunk % 8)) - 1);
        if flags_raw[flags_raw.len() - 1] & mask != 0 {
            bail!("chunk {chunk_index} mate {mate}: non-zero unused mapped_flags bits");
        }
    }
    let is_mapped = |i: usize| (flags_raw[i / 8] >> (i % 8)) & 1 == 1;

    let e_pos = find(StreamRole::Positions)?;
    // `record_count` is untrusted (footer); clamp to usize (M sizes the
    // length-M mapped-position / strand / edit vectors).
    let m = usize::try_from(e_pos.record_count)?;
    // M ≤ N and F = N − M (re-derived here; the footer validator also enforces it).
    let f = n_chunk.checked_sub(m).ok_or_else(|| {
        anyhow::anyhow!("chunk {chunk_index} mate {mate}: M ({m}) > N ({n_chunk})")
    })?;
    // The set-bit count of mapped_flags MUST equal M. The footer validator defers
    // this content invariant (§9.3(c)) to decode time, so enforce it HERE — before
    // the reconstruction below indexes the length-M mapped_positions / strand_raw /
    // edit_recs vectors by mapped slot. Without this, a crafted archive whose flags
    // popcount exceeds M (e.g. M=0 with a bit set) would drive an out-of-bounds index
    // panic (a reachable untrusted-input DoS); a popcount below M would silently
    // mis-reconstruct reads.
    let mapped_count = (0..n_chunk).filter(|&i| is_mapped(i)).count();
    if mapped_count != m {
        bail!("chunk {chunk_index} mate {mate}: mapped_flags popcount {mapped_count} != M {m}");
    }
    let pos_raw = decode_bsc_role_file(reader, e_pos)?;

    let e_strand = find(StreamRole::Strands)?;
    let strand_raw = decode_bsc_role_file(reader, e_strand)?;
    if strand_raw.len() != m.div_ceil(8) {
        bail!("chunk {chunk_index} mate {mate}: strands wrong length");
    }

    let e_rl = find(StreamRole::ReadLen)?;
    let rl_raw = decode_bsc_role_file(reader, e_rl)?;

    let e_ec = find(StreamRole::EditCount)?;
    let e_ep = find(StreamRole::EditPos)?;
    let e_eb = find(StreamRole::EditBase)?;
    let edit_recs = decode_edits(
        &decode_bsc_role_file(reader, e_ec)?,
        &decode_bsc_role_file(reader, e_ep)?,
        &decode_bsc_role_file(reader, e_eb)?,
        m,
    )?;

    // --- Decode the M mapped positions to (ref_id, ref_pos) genome coords. ---
    let mapped_positions: Vec<(u32, u64)> = if mate == 2 {
        // Build the per-mapped-R2 anchor list: for each mapped R2 (in read order),
        // the SAME pair's R1 anchor.
        let r1_anchor = r1_anchor_per_read
            .ok_or_else(|| anyhow::anyhow!("mate 2 decode without R1 anchors"))?;
        if r1_anchor.len() != n_chunk {
            bail!(
                "chunk {chunk_index}: R1 anchor count {} != R2 N_chunk {n_chunk}",
                r1_anchor.len()
            );
        }
        let mut per_mapped_anchor: Vec<Option<(u32, u64)>> = Vec::with_capacity(m);
        for i in 0..n_chunk {
            if is_mapped(i) {
                per_mapped_anchor.push(r1_anchor[i]);
            }
        }
        positions::decode_positions_mate2(&pos_raw, &per_mapped_anchor, m)?
    } else {
        positions::decode_positions_mate1(&pos_raw, m)?
    };

    // --- Pre-pass 1: parse all M read lengths up front. The ReadLen stream is a
    //     sequential varint run, so it cannot be indexed in parallel — drain it once
    //     into a Vec and re-check exact consumption (the old loop's `rl_o` check). ---
    let mut read_lens: Vec<usize> = Vec::with_capacity(m);
    let mut rl_o = 0usize;
    for _ in 0..m {
        let rl = read_varint(&rl_raw, &mut rl_o)
            .ok_or_else(|| anyhow::anyhow!("read_len: truncated varint"))?;
        let rl =
            u32::try_from(rl).map_err(|_| anyhow::anyhow!("read_len overflows u32"))? as usize;
        read_lens.push(rl);
    }
    if rl_o != rl_raw.len() {
        bail!("chunk {chunk_index} mate {mate}: read_len has trailing bytes");
    }

    // --- Pre-pass 2: pull all F fallback literals up front (the FallbackPool reader
    //     is sequential), then assert the pool is fully drained (the old "trailing
    //     fallback literals" check). ---
    let mut fallback_lits: Vec<Vec<u8>> = Vec::with_capacity(f);
    for _ in 0..f {
        let lit = pool
            .as_mut()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "chunk {chunk_index} mate {mate}: fallback read but no FallbackPool entry"
                )
            })?
            .next_literal()?
            .ok_or_else(|| {
                anyhow::anyhow!("chunk {chunk_index} mate {mate}: fallback pool exhausted early")
            })?;
        fallback_lits.push(lit);
    }
    if let Some(p) = pool.as_mut()
        && p.next_literal()?.is_some()
    {
        bail!("chunk {chunk_index} mate {mate}: trailing fallback literals");
    }

    // --- Pre-pass 3: per-read plan from the flag prefix-popcount, so each read
    //     indexes its own column slot independently (no running counters inside the
    //     parallel loop). ---
    enum ReadPlan {
        Mapped(usize),
        Fallback(usize),
    }
    let mut plan: Vec<ReadPlan> = Vec::with_capacity(n_chunk);
    let (mut next_mapped, mut next_fallback) = (0usize, 0usize);
    for i in 0..n_chunk {
        if is_mapped(i) {
            plan.push(ReadPlan::Mapped(next_mapped));
            next_mapped += 1;
        } else {
            plan.push(ReadPlan::Fallback(next_fallback));
            next_fallback += 1;
        }
    }
    if next_mapped != m || next_fallback != f {
        bail!("chunk {chunk_index} mate {mate}: mapped/fallback count mismatch");
    }

    // --- Reconstruct reads in PARALLEL. Each mapped read validates its own bounds
    //     (ref_id < num_refs, span ⊂ one interval, span within backing, edit.pos <
    //     read_len) then `reconstruct_read` (applies edits + RC); each fallback read
    //     clones its pre-pulled literal. `backing`/`imap`/`meta` are shared read-only
    //     refs; `backing.slice` is `&self` returning an owned Vec. ---
    use rayon::prelude::*;
    let reconstructed: Vec<(Vec<u8>, Option<(u32, u64)>)> = plan
        .par_iter()
        .enumerate()
        .map(|(i, p)| -> Result<(Vec<u8>, Option<(u32, u64)>)> {
            match *p {
                ReadPlan::Fallback(fidx) => Ok((fallback_lits[fidx].clone(), None)),
                ReadPlan::Mapped(midx) => {
                    let (ref_id, ref_pos) = mapped_positions[midx];
                    let read_len = read_lens[midx];
                    let is_rc = (strand_raw[midx / 8] >> (midx % 8)) & 1 == 1;
                    let edits = &edit_recs[midx];

                    // Bounds: ref_id < num_refs; ref_pos + read_len ≤ contig_len.
                    if (ref_id as u64) >= meta.num_refs {
                        bail!(
                            "chunk {chunk_index} mate {mate} read {i}: ref_id {ref_id} >= num_refs {}",
                            meta.num_refs
                        );
                    }
                    let span_end = ref_pos.checked_add(read_len as u64).ok_or_else(|| {
                        anyhow::anyhow!(
                            "chunk {chunk_index} mate {mate} read {i}: ref_pos+read_len overflow"
                        )
                    })?;
                    if span_end > meta.contig_len(ref_id) {
                        bail!(
                            "chunk {chunk_index} mate {mate} read {i}: ref_pos {ref_pos} + read_len {read_len} > contig_len {}",
                            meta.contig_len(ref_id)
                        );
                    }
                    // Resolve to a flat backing base and validate the span ⊂ one interval.
                    let flat = imap.ref_to_base(ref_id, ref_pos).ok_or_else(|| {
                        anyhow::anyhow!(
                            "chunk {chunk_index} mate {mate} read {i}: ({ref_id},{ref_pos}) not covered"
                        )
                    })?;
                    if !imap.span_in_single_interval(flat, read_len as u64) {
                        bail!("chunk {chunk_index} mate {mate} read {i}: span not in one interval");
                    }
                    let b = flat as usize;
                    if b + read_len > backing.len() {
                        bail!("chunk {chunk_index} mate {mate} read {i}: span past backing end");
                    }
                    for e in edits {
                        if (e.pos as usize) >= read_len {
                            bail!(
                                "chunk {chunk_index} mate {mate} read {i}: edit pos {} >= read_len {read_len}",
                                e.pos
                            );
                        }
                    }
                    Ok((
                        reconstruct_read_from_backing(backing, b, read_len, edits, is_rc),
                        Some((ref_id, ref_pos)),
                    ))
                }
            }
        })
        .collect::<Result<Vec<_>>>()?;

    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(n_chunk);
    let mut anchors: Vec<Option<(u32, u64)>> = Vec::with_capacity(n_chunk);
    for (s, a) in reconstructed {
        seqs.push(s);
        anchors.push(a);
    }
    Ok((seqs, anchors))
}

#[cfg(test)]
mod quality_codec_tests {
    #[test]
    fn reference_quality_entries_are_fqzcomp() {
        // Scope: this asserts the codec BYTE on an all-mapped dataset. End-to-end
        // byte-identical roundtrip over rev-comp + unmappable (fallback) + multi-chunk
        // inputs is covered by `reference_roundtrip_byte_identical` and the chr20 e2e.
        use super::*; // brings `open_ref` + `StreamRole` into scope
        use crate::cli::{AdvancedOptions, CompressConfig, QualityMode, ReferenceOptions};

        // Relies on the DEFAULT chunk size (200 pairs ⇒ 1 chunk). No test mutates
        // `QZ_REF_CHUNK` (multi-chunk tests set `chunk_records` directly), so there is
        // nothing to guard against.

        let d = tempfile::tempdir().unwrap();
        let refseq: String =
            std::iter::repeat_n("ACGTTGCAACGTACGTTTGACCGTACGTACGTACGTACGT", 50).collect();
        let refp = d.path().join("ref.fa");
        std::fs::write(&refp, format!(">c\n{refseq}\n")).unwrap();
        let (mut r1, mut r2) = (String::new(), String::new());
        for i in 0..200usize {
            let s = (i * 3) % (refseq.len() - 60);
            let win = &refseq[s..s + 60];
            let q: String = std::iter::repeat_n('I', 60).collect();
            r1.push_str(&format!("@r{i}/1\n{win}\n+\n{q}\n"));
            r2.push_str(&format!("@r{i}/2\n{win}\n+\n{q}\n"));
        }
        let r1p = d.path().join("R1.fastq");
        let r2p = d.path().join("R2.fastq");
        std::fs::write(&r1p, &r1).unwrap();
        std::fs::write(&r2p, &r2).unwrap();
        let out = d.path().join("out.qz");
        let c = CompressConfig {
            input: vec![r1p, r2p],
            output: out.clone(),
            working_dir: d.path().to_path_buf(),
            threads: 1,
            no_quality: false,
            fasta: false,
            quality_mode: QualityMode::Lossless,
            ultra: None,
            force: true,
            advanced: AdvancedOptions::default(),
            reference: Some(ReferenceOptions {
                reference: refp,
                reference_index: None,
                reference_fast: false,
                reference_window: 2,
            }),
            cluster: None,
            require_prebuilt_index: false,
        };
        crate::compression::compress(&c).unwrap();

        let (_reader, dir, _hdr) = open_ref(&out).unwrap();
        // Both R1Qual and R2Qual fold to StreamRole::Qual in the unified directory.
        let quals: Vec<u8> = dir
            .entries
            .iter()
            .filter(|e| e.role == StreamRole::Qual)
            .map(|e| e.codec)
            .collect();
        assert!(!quals.is_empty(), "no quality entries");
        assert!(
            quals.iter().all(|&x| x == super::super::format::CODEC_FQZCOMP),
            "quality roles must be fqzcomp(6): {quals:?}"
        );
    }

    #[test]
    fn reference_honors_quality_ctx_block_size() {
        // Follow-up B: the user's --quality-ctx-block-size must control the fqzcomp
        // sub-block granularity in reference mode (it was previously hardcoded to the
        // default). With 200 mapped pairs and block size 64, each mate's quality entry
        // must frame ceil(200/64) = 4 sub-blocks.
        use super::*; // brings `open_ref` + `StreamRole` into scope
        use crate::cli::{AdvancedOptions, CompressConfig, QualityMode, ReferenceOptions};

        // Relies on the DEFAULT chunk size (200 pairs ⇒ 1 chunk so the chunk-0 quality
        // entry frames all 200 reads). No test mutates `QZ_REF_CHUNK`, so there is
        // nothing to guard against.

        let d = tempfile::tempdir().unwrap();
        let refseq: String =
            std::iter::repeat_n("ACGTTGCAACGTACGTTTGACCGTACGTACGTACGTACGT", 50).collect();
        let refp = d.path().join("ref.fa");
        std::fs::write(&refp, format!(">c\n{refseq}\n")).unwrap();
        let (mut r1, mut r2) = (String::new(), String::new());
        for i in 0..200usize {
            let s = (i * 3) % (refseq.len() - 60);
            let win = &refseq[s..s + 60];
            let q: String = std::iter::repeat_n('I', 60).collect();
            r1.push_str(&format!("@r{i}/1\n{win}\n+\n{q}\n"));
            r2.push_str(&format!("@r{i}/2\n{win}\n+\n{q}\n"));
        }
        let r1p = d.path().join("R1.fastq");
        let r2p = d.path().join("R2.fastq");
        std::fs::write(&r1p, &r1).unwrap();
        std::fs::write(&r2p, &r2).unwrap();
        let out = d.path().join("out.qz");
        // small block size forces multiple fqzcomp sub-blocks
        let adv = AdvancedOptions {
            quality_ctx_block_size: 64,
            ..Default::default()
        };
        let c = CompressConfig {
            input: vec![r1p, r2p],
            output: out.clone(),
            working_dir: d.path().to_path_buf(),
            threads: 1,
            no_quality: false,
            fasta: false,
            quality_mode: QualityMode::Lossless,
            ultra: None,
            force: true,
            advanced: adv,
            reference: Some(ReferenceOptions {
                reference: refp,
                reference_index: None,
                reference_fast: false,
                reference_window: 2,
            }),
            cluster: None,
            require_prebuilt_index: false,
        };
        crate::compression::compress(&c).unwrap();

        let (reader, dir, _hdr) = open_ref(&out).unwrap();
        // R1Qual folds to StreamRole::Qual, mate 1, chunk 0.
        let e = dir
            .entries
            .iter()
            .find(|e| e.chunk_index == 0 && e.mate == 1 && e.role == StreamRole::Qual)
            .expect("R1Qual entry for chunk 0");
        let payload = reader.read_entry(e.offset, e.length).unwrap();
        let blocks = crate::compression::paired::streams::decode_block_payloads(&payload).unwrap();
        assert_eq!(
            blocks.len(),
            4,
            "block size 64 over 200 reads must frame 4 fqzcomp sub-blocks"
        );
        assert_eq!(
            blocks.iter().map(|(c, _)| *c).sum::<u32>(),
            200,
            "sub-block counts must sum to 200"
        );
    }
}


#[cfg(test)]
mod streaming_decode_oracle {
    //! Lossless oracle for the engine-aligned reference decode
    //! (`decode_reference_streaming`): compress a synthetic paired dataset mixing
    //! mapped / reverse-complement / planted-edit / unmappable (fallback) reads
    //! across MULTIPLE chunks, decode through the production streaming path, and
    //! assert the output is byte-identical to the original R1/R2 input bytes. The
    //! multi-chunk shape pins that the flat header/qual engine cursors stay
    //! record-aligned with the per-chunk reconstructed sequences across chunk
    //! boundaries (header/qual record i must line up with reconstructed sequence i).

    use super::*; // open_ref, decode_reference_streaming
    use crate::cli::{AdvancedOptions, CompressConfig, QualityMode, ReferenceOptions};
    use std::path::Path;

    /// Deterministic low-repeat ACGT sequence (matches the integration generators
    /// so substrings seed reliably against the reference).
    fn make_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            v.push(b"ACGT"[((x >> 33) & 3) as usize]);
        }
        v
    }

    fn revcomp(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => other,
            })
            .collect()
    }

    /// Build a synthetic paired dataset against a deterministic reference: a MIX of
    /// mapped reads (R1 forward substrings, R2 reverse-complement downstream slices),
    /// a handful with planted substitutions (edit path), and two unmappable random
    /// pairs (literal fallback path). Writes the FASTAs/FASTQs under `dir` and returns
    /// `(ref_path, r1_path, r2_path, r1_bytes, r2_bytes)` — the strings are the EXACT
    /// bytes a byte-identical decode must reproduce.
    fn make_dataset(
        dir: &Path,
    ) -> (
        std::path::PathBuf,
        std::path::PathBuf,
        std::path::PathBuf,
        String,
        String,
    ) {
        let refseq = make_seq(2000, 7);
        let rf = dir.join("ref.fa");
        std::fs::write(
            &rf,
            format!(">chr0\n{}\n", std::str::from_utf8(&refseq).unwrap()),
        )
        .unwrap();

        let n_pairs = 200usize;
        let mut s1 = String::new();
        let mut s2 = String::new();
        let mut rng = 0x1234_5678_9abc_def0u64;
        let mut next = || {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (rng >> 33) as usize
        };
        for i in 0..n_pairs {
            let unmappable = i >= n_pairs - 2; // last two: random ⇒ fallback
            let rlen = 100 + (next() % 21); // 100..=120
            let (r1bytes, r2bytes) = if unmappable {
                (
                    make_seq(rlen, 9000 + i as u64),
                    make_seq(rlen, 7000 + i as u64),
                )
            } else {
                let max_start = refseq.len() - rlen;
                let st1 = next() % (max_start + 1);
                let mut r1 = refseq[st1..st1 + rlen].to_vec();
                let st2 = (st1 + 150).min(max_start);
                let mut r2 = revcomp(&refseq[st2..st2 + rlen]);
                if i % 17 == 0 {
                    let p = next() % rlen;
                    r1[p] = if r1[p] == b'A' { b'C' } else { b'A' };
                }
                if i % 23 == 0 {
                    let p = next() % rlen;
                    r2[p] = if r2[p] == b'G' { b'T' } else { b'G' };
                }
                (r1, r2)
            };
            let q: String = "I".repeat(r1bytes.len());
            s1.push_str(&format!(
                "@read_{i}/1\n{}\n+\n{q}\n",
                std::str::from_utf8(&r1bytes).unwrap()
            ));
            let q2: String = "I".repeat(r2bytes.len());
            s2.push_str(&format!(
                "@read_{i}/2\n{}\n+\n{q2}\n",
                std::str::from_utf8(&r2bytes).unwrap()
            ));
        }
        let r1p = dir.join("r1.fastq");
        let r2p = dir.join("r2.fastq");
        std::fs::write(&r1p, &s1).unwrap();
        std::fs::write(&r2p, &s2).unwrap();
        (rf, r1p, r2p, s1, s2)
    }

    /// Decode `archive` into two in-memory sinks via the production engine-aligned
    /// streaming path. A fresh `open_ref` per call (the `RefReader` is consumed by
    /// value).
    fn decode_at(archive: &Path) -> (Vec<u8>, Vec<u8>) {
        let (reader, dir, hdr) = open_ref(archive).unwrap();
        let mut r1 = Vec::new();
        let mut r2 = Vec::new();
        decode_reference_streaming(
            archive,
            reader,
            &dir,
            &hdr,
            &mut RefSink::Two(&mut r1, &mut r2),
            None,
        )
        .unwrap();
        (r1, r2)
    }

    #[test]
    fn reference_decode_streaming_is_lossless_multichunk() {
        // Force ≥3 chunks (200 pairs / 64 ⇒ 4 chunks) by driving `chunk_records`
        // directly, so the engine-aligned driver walks multiple chunks and the flat
        // header/qual cursors must stay record-aligned with the per-chunk
        // reconstructed sequences across chunk boundaries.
        let d = tempfile::tempdir().unwrap();
        let (refp, r1p, r2p, s1, s2) = make_dataset(d.path());
        let out = d.path().join("oracle.qz");
        let c = CompressConfig {
            input: vec![r1p, r2p],
            output: out.clone(),
            working_dir: d.path().to_path_buf(),
            threads: 1,
            no_quality: false,
            fasta: false,
            quality_mode: QualityMode::Lossless,
            ultra: None,
            force: true,
            advanced: AdvancedOptions {
                chunk_records: 64,
                ..AdvancedOptions::default()
            },
            reference: Some(ReferenceOptions {
                reference: refp,
                reference_index: None,
                reference_fast: false,
                reference_window: 2,
            }),
            cluster: None,
            require_prebuilt_index: false,
        };
        crate::compression::compress(&c).unwrap();

        // Confirm the archive genuinely spans ≥3 chunks (otherwise the multi-chunk
        // cursor alignment is never exercised — a single chunk is trivially aligned).
        let (_r, dir, _hdr) = open_ref(&out).unwrap();
        assert!(
            dir.num_chunks >= 3,
            "oracle archive must span ≥3 chunks to exercise the multi-chunk driver (got {})",
            dir.num_chunks
        );

        // Lossless: the engine-aligned streaming decode reproduces the input bytes.
        let (got_r1, got_r2) = decode_at(&out);
        assert_eq!(got_r1, s1.as_bytes(), "R1 streaming decode not lossless");
        assert_eq!(got_r2, s2.as_bytes(), "R2 streaming decode not lossless");
    }
}

#[cfg(test)]
mod single_streaming_decode_oracle {
    //! Byte-identity oracle for the single-end reference decode
    //! (`decompress_reference_single`): compress a synthetic single-end dataset
    //! mixing mapped / reverse-complement / planted-edit / unmappable (fallback)
    //! reads across MULTIPLE chunks, decode to ONE FASTQ, and assert the output
    //! is byte-identical to the original input bytes.
    use crate::cli::{
        AdvancedOptions, CompressConfig, DecompressConfig, QualityMode, ReferenceOptions,
    };

    fn make_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            v.push(b"ACGT"[((x >> 33) & 3) as usize]);
        }
        v
    }
    fn revcomp(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                o => o,
            })
            .collect()
    }

    #[test]
    fn reference_single_decode_is_lossless_multichunk() {
        let d = tempfile::tempdir().unwrap();
        let refseq = make_seq(2000, 7);
        let refp = d.path().join("ref.fa");
        std::fs::write(
            &refp,
            format!(">chr0\n{}\n", std::str::from_utf8(&refseq).unwrap()),
        )
        .unwrap();

        let n = 200usize;
        let mut s = String::new();
        let mut rng = 0x1234_5678_9abc_def0u64;
        let mut next = || {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (rng >> 33) as usize
        };
        for i in 0..n {
            let unmappable = i >= n - 2; // last two: random ⇒ fallback
            let rlen = 100 + (next() % 21);
            let bytes = if unmappable {
                make_seq(rlen, 9000 + i as u64)
            } else {
                let max_start = refseq.len() - rlen;
                let st = next() % (max_start + 1);
                let mut r = refseq[st..st + rlen].to_vec();
                if i % 5 == 0 {
                    // some reverse-complement reads (RC path)
                    r = revcomp(&r);
                }
                if i % 17 == 0 {
                    let p = next() % rlen;
                    r[p] = if r[p] == b'A' { b'C' } else { b'A' }; // planted edit
                }
                r
            };
            let q: String = "I".repeat(bytes.len());
            s.push_str(&format!(
                "@read_{i}\n{}\n+\n{q}\n",
                std::str::from_utf8(&bytes).unwrap()
            ));
        }
        let rp = d.path().join("reads.fastq");
        std::fs::write(&rp, &s).unwrap();
        let out = d.path().join("single.qz");

        let c = CompressConfig {
            input: vec![rp],
            output: out.clone(),
            working_dir: d.path().to_path_buf(),
            threads: 1,
            force: true,
            advanced: AdvancedOptions {
                chunk_records: 64, // ≥3 chunks
                ..AdvancedOptions::default()
            },
            reference: Some(ReferenceOptions {
                reference: refp,
                reference_index: None,
                reference_fast: false,
                reference_window: 2,
            }),
            quality_mode: QualityMode::Lossless,
            ..CompressConfig::default()
        };
        crate::compression::compress(&c).unwrap();

        // Genuinely multi-chunk + actually type 4.
        let (_r, dir, hdr) = super::open_ref_single(&out).unwrap();
        assert_eq!(hdr.archive_type, 4, "must be single-end reference");
        assert!(dir.num_chunks >= 3, "want ≥3 chunks, got {}", dir.num_chunks);

        let decp = d.path().join("dec.fastq");
        let dc = DecompressConfig {
            input: out,
            output: vec![decp.clone()],
            working_dir: d.path().to_path_buf(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 6,
            force: true,
        };
        super::decompress_reference_single(&dc).unwrap();
        assert_eq!(
            std::fs::read_to_string(&decp).unwrap(),
            s,
            "single-end reference decode not byte-identical"
        );
    }
}
