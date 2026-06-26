//! Compression orchestrators: chunked pipelined compression with disk-backed block storage.

use super::codecs;
use super::*;
use crate::cli::{CompressConfig, QualityCompressor, QualityMode};
use anyhow::Result;
use std::time::Instant;
use tracing::info;

// ---------------------------------------------------------------------------
// Per-stream BSC / fqz block helpers (shared by the streaming engine + reference)
// ---------------------------------------------------------------------------

/// Compress quality strings into one fqzcomp blob; if it exceeds `cap`,
/// recursively bisect by record count until every emitted blob is `<= cap`.
/// Hard-errors if a single record's blob still exceeds the cap (pathological;
/// real reads are bounded). Order-preserving. Returns framed `(record_count, blob)`.
///
/// This is the **byte-capped** quality splitter used by reference mode (which
/// bisects by output size, not record count). Every element must be a real
/// quality string. It operates on borrowed quality strings so reference mode can
/// encode quality without building throwaway records.
pub(crate) fn fqz_blocks_capped_quals(quals: &[&[u8]], cap: usize) -> Result<Vec<(u32, Vec<u8>)>> {
    let blob = codecs::compress_qualities_fqzcomp_quals(quals)?;
    if blob.len() <= cap {
        return Ok(vec![(quals.len() as u32, blob)]);
    }
    if quals.len() <= 1 {
        anyhow::bail!(
            "fqzcomp: single-record quality blob {} B exceeds block cap {} B",
            blob.len(),
            cap
        );
    }
    let mid = quals.len() / 2;
    let mut out = fqz_blocks_capped_quals(&quals[..mid], cap)?;
    out.extend(fqz_blocks_capped_quals(&quals[mid..], cap)?);
    Ok(out)
}

/// Columnar-encode `headers` into framed `(record_count, payload)` blocks where
/// each payload is `[count u32][columnar_blob]` and each block stays `<= cap` bytes
/// on disk (the block-stream framing rejects any block over `streams::MAX_BLOCK`).
///
/// The common case — the chunk's columnar blob fits under `cap` — returns a SINGLE
/// block byte-identical to the pre-cap encoder, so every existing archive is
/// unchanged and there is no ratio cost. Only when a chunk's columnar blob would
/// exceed the 64 MiB block cap (reachable at large `chunk_records` with long,
/// high-entropy read IDs) does this bisect by record count and recurse, so every
/// emitted block is independently decodable — without the cap such a chunk produced
/// an archive that failed its own decode/verify. Bails if a single header's payload
/// still exceeds `cap` (pathological; real read IDs are bounded). Order-preserving.
///
/// Mirrors [`fqz_blocks_capped_quals`]: same cap-triggered bisection, same framing.
/// All three columnar-header consumers handle the resulting N blocks: the single-end
/// producer (`bounded_columnar_header_producer`) and the paired/reference decoder
/// (`paired::decode_columnar_headers`) decode each block independently and
/// concatenate in block order.
/// Records per columnar header sub-block. The header path was historically the only
/// stream still operating at whole-chunk granularity (sequence/quality are split into
/// many block-sized jobs), which made each 2.5M-record chunk's columnar header one
/// coarse ~5s unit that straggled at the end of compress while most workers sat idle.
/// Splitting it into independent sub-range blocks (compressed in parallel, the decoders
/// already concatenate N blocks in order) collapses that tail. 250K is the measured
/// sweet spot: header-job ~4.4× faster for ~+0.15% archive size (each sub-range loses a
/// little cross-boundary columnar context + rebuilds its own combo dict). `QZ_HEADER_SUB`
/// overrides for tuning (0 = whole-chunk, the legacy single-unit behavior).
pub(crate) const HEADER_SUB_BLOCK_RECORDS: usize = 250_000;

pub(crate) fn columnar_blocks_capped(headers: &[&[u8]], cap: usize) -> Result<Vec<(u32, Vec<u8>)>> {
    let sub = std::env::var("QZ_HEADER_SUB")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(HEADER_SUB_BLOCK_RECORDS);
    if sub == 0 || headers.len() <= sub {
        return columnar_blocks_capped_one(headers, cap);
    }
    use rayon::prelude::*;
    let ranges: Vec<&[&[u8]]> = headers.chunks(sub).collect();
    let parts: Vec<Result<Vec<(u32, Vec<u8>)>>> = ranges
        .par_iter()
        .map(|r| columnar_blocks_capped_one(r, cap))
        .collect();
    let mut out = Vec::new();
    for p in parts {
        out.extend(p?);
    }
    Ok(out)
}

/// Compress one (sub-)range of headers as columnar, bisecting only if the compressed
/// payload exceeds `cap` (the original `columnar_blocks_capped` behavior).
fn columnar_blocks_capped_one(headers: &[&[u8]], cap: usize) -> Result<Vec<(u32, Vec<u8>)>> {
    let blob = header_col::compress_headers_columnar_bytes(headers)?;
    let mut payload = Vec::with_capacity(4 + blob.len());
    payload.extend_from_slice(&(headers.len() as u32).to_le_bytes());
    payload.extend_from_slice(&blob);
    if payload.len() <= cap {
        return Ok(vec![(headers.len() as u32, payload)]);
    }
    if headers.len() <= 1 {
        anyhow::bail!(
            "columnar headers: single-record block {} B exceeds block cap {} B",
            payload.len(),
            cap
        );
    }
    let mid = headers.len() / 2;
    let mut out = columnar_blocks_capped_one(&headers[..mid], cap)?;
    out.extend(columnar_blocks_capped_one(&headers[mid..], cap)?);
    Ok(out)
}

/// Unified chunked compression pipeline.
///
/// Replaces the former compress_chunked_bsc, compress_chunked_bsc_reorder_local,
/// and compress_chunked_fqzcomp with a single function parameterised by `mode`.
///
/// Quality strategy is auto-detected:
/// - Fqzcomp: lossless + ≥100K reads (explicit Fqzcomp, or Auto)
/// - BSC: otherwise
///
/// Features: pipelined I/O (read next chunk while compressing current), parallel
/// h||s||q compression, constant-length detection, cross-chunk validation.
/// How `compress_chunked` is parameterised.
///
/// - `Default`: 2.5M-record chunks, ≤64 MiB blocks, RawWithHints/Raw encoding.
/// - `Ultra`: big blocks (≤750 MiB), `UltraBigBlock` encoding, raw BSC, no reorder.
#[derive(Clone, Copy, PartialEq)]
pub(super) enum ChunkedMode {
    Default,
    Ultra {
        bsc_block_mb: usize,
        chunk_records: usize,
        quality_sub_block: usize,
        compress_window: usize,
    },
}

/// Write one stream's blocks for a chunk to `output_file`, appending a directory
/// entry per block. `offset`/`length` in each entry describe the **full v5 block
/// frame** (12-byte prefix + payload), and `byte_cursor` tracks the absolute output
/// position so the next frame's offset is correct. Used by the v5 chunk-major
/// non-sort pipeline. (Top-level `fn` rather than a closure so the borrow checker
/// does not fight capturing `output_file`/`byte_cursor`/`directory` while these
/// are also live across the FIFO-drain loop.)
// Generic over the writer (`W: Write + ?Sized`) so it serves BOTH the legacy path's
// `Box<dyn Write>` and the streaming writer's `Box<dyn Write + Send>` (the latter
// moved into a spawned thread). `Box<dyn Write>` and `Box<dyn Write + Send>` both
// implement `Write`, so legacy call sites (`&mut output_file`) compile unchanged and
// the framing/CRC/offset math stays byte-identical across both engines.
#[allow(clippy::too_many_arguments)]
pub(super) fn write_role_blocks<W: std::io::Write + ?Sized>(
    blocks: &[(u32, Vec<u8>)],
    role: crate::compression::chunk_directory::StreamRole,
    mate: u8,
    codec: u8,
    chunk_index: u32,
    output_file: &mut W,
    byte_cursor: &mut u64,
    directory: &mut Vec<crate::compression::chunk_directory::ChunkDirEntry>,
) -> Result<()> {
    use crate::compression::BLOCK_PREFIX_SIZE;
    let mut framing = Vec::with_capacity(BLOCK_PREFIX_SIZE);
    for (record_count, block) in blocks {
        framing.clear();
        super::bsc::write_block_frame_header(&mut framing, block.len() as u32, *record_count, block);
        let frame_offset = *byte_cursor;
        let frame_len = (framing.len() + block.len()) as u64;
        output_file.write_all(&framing)?;
        output_file.write_all(block)?;
        directory.push(crate::compression::chunk_directory::ChunkDirEntry {
            chunk_index,
            role,
            mate,
            codec,
            offset: frame_offset,
            length: frame_len,
            record_count: *record_count as u64,
        });
        *byte_cursor += frame_len;
    }
    Ok(())
}

/// Paired physical contract: one ChunkDirEntry whose payload is a `write_block_stream`
/// blob (`[num_blocks][framed blocks…]`) for the whole (mate,role) segment. `record_count`
/// is the chunk's pair count (NOT per-block). Seek-free.
#[allow(clippy::too_many_arguments)]
fn write_segment_blob<W: std::io::Write + ?Sized>(
    blocks: &[(u32, Vec<u8>)],
    role: crate::compression::chunk_directory::StreamRole,
    mate: u8,
    codec: u8,
    record_count: u64,
    chunk_index: u32,
    output_file: &mut W,
    byte_cursor: &mut u64,
    directory: &mut Vec<crate::compression::chunk_directory::ChunkDirEntry>,
) -> Result<()> {
    let mut blob = Vec::new();
    crate::compression::paired::streams::write_block_stream(&mut blob, blocks);
    output_file.write_all(&blob)?;
    directory.push(crate::compression::chunk_directory::ChunkDirEntry {
        chunk_index,
        role,
        mate,
        codec,
        offset: *byte_cursor,
        length: blob.len() as u64,
        record_count,
    });
    *byte_cursor += blob.len() as u64;
    Ok(())
}

pub(super) fn compress_chunked(args: &CompressConfig, mode: ChunkedMode) -> Result<()> {
    compress_chunked_streaming(args, mode)
}

// ===========================================================================
// Block-granular streaming single-end compress engine
// ===========================================================================
//
// Producer (this thread) frames records, stages per-role raw blocks, and dispatches
// them to a bounded worker pool. Workers compress and forward results to a single
// writer thread, which orders blocks by chunk and writes the per-chunk role sequence
// + ChunkDecodedSizes global + footer/locator.

use crate::compression::chunk_directory::StreamRole;
use crate::compression::stream_compress::{
    ArchivePlan, BlockResult, ChunkManifest, CompletedSegment, CompressJob, HeaderColStager,
    JobOutput, OrderedBlockBuffer, QualStager, SeqQualHeaderStager, compress_job,
    compress_worker_count, resolve_archive_plan,
};
use std::sync::mpsc::SyncSender;

/// Selects the per-mode physical write contract in `write_ordered` (spec §6).
#[derive(Clone, Copy, PartialEq)]
pub(crate) enum CompressMode {
    /// One ChunkDirEntry PER BLOCK (H/S/Q/RC), mate 0, + ChunkDecodedSizes global.
    Single,
    /// One ChunkDirEntry per (mate,role) SEGMENT (write_block_stream blob), no
    /// ChunkDecodedSizes; validate_paired_directory before the footer.
    Paired,
}

/// Messages on the single results channel the writer drains. `Manifest` is sent by
/// the producer (per chunk close); `Block`/`Segment`/`Err` by workers.
pub(crate) enum WorkerMsg {
    Manifest(ChunkManifest),
    Block(BlockResult),
    Segment(CompletedSegment),
    Err(String),
}

/// What the writer thread returns on a clean finish — the per-role compressed byte
/// totals + original size + metadata size for the final stats log line.
pub(crate) struct WriterDone {
    pub original_size: usize,
    pub headers_len: usize,
    pub sequences_len: usize,
    pub qualities_len: usize,
    pub rc_flags_len: usize,
    pub metadata_size: usize,
}

/// Validate one record's lengths against the chunk-0 global const-length profile.
/// Mirrors legacy `validate_const_lengths` semantics for a single record: when a
/// global const length is set (non-zero), every subsequent record MUST match it, or
/// const-length framing (no per-record varint) would silently mis-slice on decode.
fn validate_const_one(
    rec: &crate::io::FastqRecord,
    const_seq_len: usize,
    const_qual_len: usize,
    no_quality: bool,
    chunk_index: u32,
) -> Result<()> {
    if const_seq_len > 0 && rec.sequence.len() != const_seq_len {
        anyhow::bail!(
            "Constant sequence length mismatch: chunk 0 had {} but chunk {} has {} (or variable).",
            const_seq_len,
            chunk_index,
            rec.sequence.len()
        );
    }
    if const_qual_len > 0 && !no_quality {
        let qlen = rec.quality.as_ref().map_or(0, |q| q.len());
        if qlen != const_qual_len {
            anyhow::bail!(
                "Constant quality length mismatch: chunk 0 had {} but chunk {} has {} (or variable).",
                const_qual_len,
                chunk_index,
                qlen
            );
        }
    }
    Ok(())
}

/// Lightweight discriminant for `QualStager`, copied out before staging so we don't
/// hold a `match &mut self.qual` borrow while calling its mutating methods.
#[derive(Clone, Copy)]
enum QualKind {
    Fqz,
    Bsc,
    None,
}

/// Per-chunk staging + dispatch. Owns the three (+rc) stagers and the per-chunk
/// manifest accumulators. Makes the same framing/role decisions as the
/// chunk-at-a-time path used to, but flushes blocks as they fill rather than
/// per-chunk.
struct ProducerState<'a> {
    plan: &'a ArchivePlan,
    seq: SeqQualHeaderStager,
    // Exactly one of these two is active, per `use_columnar_headers`.
    hdr_bsc: Option<SeqQualHeaderStager>,
    hdr_col: Option<HeaderColStager>,
    qual: QualStager,
    rc: Option<SeqQualHeaderStager>,
    // Scratch framing buffers, cleared per record.
    fh: Vec<u8>,
    fs: Vec<u8>,
    fq: Vec<u8>,
    frc: Vec<u8>,
    // Per-chunk accumulators (all reset at chunk close). The writer sums the
    // per-chunk records/total_bases/original_size across chunks, mirroring legacy's
    // `num_reads += reads` / `original_size += orig` per FIFO drain.
    chunk_index: u32,
    sequence_blocks: u32,
    header_blocks: u32,
    qual_blocks: u32,
    rc_blocks: u32,
    decoded_output_bytes: u64,
    chunk_records: u64,       // records staged in the current (open) chunk
    chunk_total_bases: u64,   // Σ sequence bases in the current chunk
    chunk_original_size: u64, // raw FASTQ bytes in the current chunk (log-only)
    /// Running total of records across ALL chunks — used only as a stable record
    /// index in error messages (mirrors legacy `records_to_streams`' diagnostics).
    records_seen: u64,
}

impl<'a> ProducerState<'a> {
    fn new(plan: &'a ArchivePlan, target: usize) -> Self {
        let (hdr_bsc, hdr_col) = if plan.use_columnar_headers {
            (None, Some(HeaderColStager::new(0)))
        } else {
            (Some(SeqQualHeaderStager::new(StreamRole::Headers, 0, target)), None)
        };
        let qual = if plan.no_quality {
            QualStager::None
        } else if plan.use_fqzcomp {
            QualStager::fqz(plan.quality_ctx_block_size, 0)
        } else {
            QualStager::bsc(0, target)
        };
        let rc = if plan.rc_canon {
            Some(SeqQualHeaderStager::new(StreamRole::RcFlags, 0, target))
        } else {
            None
        };
        Self {
            plan,
            seq: SeqQualHeaderStager::new(StreamRole::Sequence, 0, target),
            hdr_bsc,
            hdr_col,
            qual,
            rc,
            fh: Vec::new(),
            fs: Vec::new(),
            fq: Vec::new(),
            frc: Vec::new(),
            chunk_index: 0,
            sequence_blocks: 0,
            header_blocks: 0,
            qual_blocks: 0,
            rc_blocks: 0,
            decoded_output_bytes: 0,
            chunk_records: 0,
            chunk_total_bases: 0,
            chunk_original_size: 0,
            records_seen: 0,
        }
    }

    fn dispatch_job(job_tx: &SyncSender<CompressJob>, job: CompressJob) -> Result<()> {
        // A send error means the pool aborted; the real cause is already queued to the
        // writer (first-error-wins), so this secondary error is harmless.
        job_tx
            .send(job)
            .map_err(|_| anyhow::anyhow!("compress workers exited"))
    }

    fn stage_one(
        &mut self,
        rec: &crate::io::FastqRecord,
        job_tx: &SyncSender<CompressJob>,
    ) -> Result<()> {
        let p = self.plan;
        self.fh.clear();
        self.fs.clear();
        self.fq.clear();
        self.frc.clear();
        // Legacy passes `skip_quality_in_stream` (= no_quality || use_fqzcomp) as
        // `records_to_streams`' `no_quality` arg, so fqz/no-quality records emit NO
        // quality into the framed stream (fqz compresses quality separately). Mirror
        // that exactly — `p.no_quality` alone would leave quality in `fq` for fqz.
        crate::compression::frame_record(
            rec,
            p.quality_mode,
            p.stream_quality_binning,
            p.skip_quality_in_stream,
            p.sequence_hints,
            p.rc_canon,
            p.const_seq_len,
            p.const_qual_len,
            !p.use_columnar_headers, // build_header_stream only for the BSC header path
            &mut self.fh,
            &mut self.fs,
            &mut self.fq,
            &mut self.frc,
        )
        .map_err(|e| anyhow::anyhow!("FASTQ record {}: {e}", self.records_seen))?;

        // Sequence: always BSC-staged.
        if let Some(b) = self.seq.stage_record_bytes(self.chunk_index, &self.fs) {
            Self::dispatch_job(job_tx, CompressJob::Bsc(b))?;
            self.sequence_blocks += 1;
        }

        // Headers: columnar buffers the id; BSC routes framed header bytes.
        if let Some(col) = self.hdr_col.as_mut() {
            col.stage(&rec.id);
        } else if let Some(hs) = self.hdr_bsc.as_mut()
            && let Some(b) = hs.stage_record_bytes(self.chunk_index, &self.fh)
        {
            Self::dispatch_job(job_tx, CompressJob::Bsc(b))?;
            self.header_blocks += 1;
        }

        // Quality: fqz stages raw ASCII; BSC routes framed (binned) quality bytes;
        // no_quality stages nothing.
        let qual_kind = match self.qual {
            QualStager::Fqz { .. } => QualKind::Fqz,
            QualStager::Bsc(_) => QualKind::Bsc,
            QualStager::None => QualKind::None,
        };
        match qual_kind {
            QualKind::Fqz => {
                // Mirrors the legacy fqz guard (compress_impl.rs:202): every record must
                // carry quality, else the framed record_count desyncs the decoder.
                let q = rec.quality.as_deref().ok_or_else(|| {
                    anyhow::anyhow!(
                        "fqzcomp quality encode: every record must carry quality \
                         (record_count framing requires it)"
                    )
                })?;
                if let Some(b) = self.qual.stage_fqz(self.chunk_index, q) {
                    Self::dispatch_job(job_tx, CompressJob::Fqz(b))?;
                    self.qual_blocks += 1;
                }
            }
            QualKind::Bsc => {
                // `fq` is empty when the record has no quality — staging it appends nothing.
                let chunk_index = self.chunk_index;
                let block = self
                    .qual
                    .bsc_stager()
                    .and_then(|s| s.stage_record_bytes(chunk_index, &self.fq));
                if let Some(b) = block {
                    Self::dispatch_job(job_tx, CompressJob::Bsc(b))?;
                    self.qual_blocks += 1;
                }
            }
            QualKind::None => {}
        }

        // rc_flags: one byte per record (from frame_record's rc buffer).
        if let Some(rc) = self.rc.as_mut()
            && let Some(b) = rc.stage_record_bytes(self.chunk_index, &self.frc)
        {
            Self::dispatch_job(job_tx, CompressJob::Bsc(b))?;
            self.rc_blocks += 1;
        }

        // Accumulators (mirror read_chunk_records).
        // Two different quality lengths, deliberately:
        //  - decoded_qual_len drives decoded_output_bytes (the ChunkDecodedSizes global,
        //    an ARCHIVE byte) and is 0 under no_quality because the reconstructed FASTQ
        //    carries no quality — this MUST stay guarded for byte-identity vs legacy.
        //  - actual_qual_len drives chunk_original_size (the raw-input-size stat, log-only)
        //    and counts the quality bytes present in the INPUT even under no_quality,
        //    matching legacy read_chunk_records (mod.rs).
        let actual_qual_len = rec.quality.as_ref().map_or(0, |q| q.len());
        let decoded_qual_len = if p.no_quality { 0 } else { actual_qual_len };
        let add = crate::compression::decompress_impl::fastq_output_len(
            p.fasta,
            rec.id.len(),
            rec.sequence.len(),
            decoded_qual_len,
        );
        self.decoded_output_bytes = self
            .decoded_output_bytes
            .checked_add(add)
            .ok_or_else(|| anyhow::anyhow!("decoded size overflow in chunk {}", self.chunk_index))?;
        self.chunk_records += 1;
        self.records_seen += 1;
        self.chunk_total_bases += rec.sequence.len() as u64;
        self.chunk_original_size += (rec.id.len() + rec.sequence.len() + actual_qual_len + 3) as u64;
        Ok(())
    }

    fn maybe_close_chunk(
        &mut self,
        chunk_size: usize,
        job_tx: &SyncSender<CompressJob>,
        manifest_tx: &SyncSender<WorkerMsg>,
    ) -> Result<()> {
        if self.chunk_records as usize >= chunk_size {
            self.close_chunk(job_tx, manifest_tx)?;
        }
        Ok(())
    }

    /// Flush all stagers for the current chunk, emit the manifest, and reset for the
    /// next chunk. A chunk with zero staged records emits NO manifest and NO blocks
    /// (matches legacy's empty-trailing-chunk behaviour).
    fn close_chunk(
        &mut self,
        job_tx: &SyncSender<CompressJob>,
        manifest_tx: &SyncSender<WorkerMsg>,
    ) -> Result<()> {
        if self.chunk_records == 0 {
            return Ok(());
        }

        // Build the segment list in physical write order: Headers → Sequence → Qual → RcFlags.
        // The DISPATCH order of jobs is irrelevant (the writer reorders by segment order);
        // only `segments` order drives the writer. The single-end goldens are the gate.
        use crate::compression::codec_ids::{
            CODEC_COLUMNAR, header_compressor_to_codec_id, quality_compressor_to_codec_id,
            sequence_compressor_to_codec_id,
        };
        use crate::compression::stream_compress::{SegSource, SegmentSpec};
        let seq_codec = sequence_compressor_to_codec_id(SequenceCompressor::Bsc);
        let qual_codec = quality_compressor_to_codec_id(self.plan.quality_compressor_used);
        let mut segments: Vec<SegmentSpec> = Vec::with_capacity(4);

        // Headers flush / columnar dispatch.
        if let Some(col) = self.hdr_col.as_mut() {
            if let Some(job) = col.take_chunk_job(self.chunk_index) {
                Self::dispatch_job(job_tx, CompressJob::Columnar(job))?;
            }
            segments.push(SegmentSpec {
                mate: 0,
                role: StreamRole::Headers,
                codec: CODEC_COLUMNAR,
                source: SegSource::HeaderSegment,
            });
        } else {
            let hs = self.hdr_bsc.as_mut().expect("bsc header stager present");
            if let Some(b) = hs.flush(self.chunk_index) {
                Self::dispatch_job(job_tx, CompressJob::Bsc(b))?;
                self.header_blocks += 1;
            }
            let hdr_codec = header_compressor_to_codec_id(self.plan.header_compressor);
            segments.push(SegmentSpec {
                mate: 0,
                role: StreamRole::Headers,
                codec: hdr_codec,
                source: SegSource::Blocks(self.header_blocks),
            });
        }

        // Sequence flush.
        if let Some(b) = self.seq.flush(self.chunk_index) {
            Self::dispatch_job(job_tx, CompressJob::Bsc(b))?;
            self.sequence_blocks += 1;
        }
        segments.push(SegmentSpec {
            mate: 0,
            role: StreamRole::Sequence,
            codec: seq_codec,
            source: SegSource::Blocks(self.sequence_blocks),
        });

        // Quality flush.
        let qual_kind = match self.qual {
            QualStager::Fqz { .. } => QualKind::Fqz,
            QualStager::Bsc(_) => QualKind::Bsc,
            QualStager::None => QualKind::None,
        };
        match qual_kind {
            QualKind::Fqz => {
                if let Some(b) = self.qual.flush_fqz(self.chunk_index) {
                    Self::dispatch_job(job_tx, CompressJob::Fqz(b))?;
                    self.qual_blocks += 1;
                }
            }
            QualKind::Bsc => {
                let chunk_index = self.chunk_index;
                let block = self.qual.bsc_stager().and_then(|s| s.flush(chunk_index));
                if let Some(b) = block {
                    Self::dispatch_job(job_tx, CompressJob::Bsc(b))?;
                    self.qual_blocks += 1;
                }
            }
            QualKind::None => {}
        }
        segments.push(SegmentSpec {
            mate: 0,
            role: StreamRole::Qual,
            codec: qual_codec,
            source: SegSource::Blocks(self.qual_blocks),
        });

        // rc flush.
        if let Some(rc) = self.rc.as_mut()
            && let Some(b) = rc.flush(self.chunk_index)
        {
            Self::dispatch_job(job_tx, CompressJob::Bsc(b))?;
            self.rc_blocks += 1;
        }
        if self.plan.rc_canon {
            segments.push(SegmentSpec {
                mate: 0,
                role: StreamRole::RcFlags,
                codec: seq_codec,
                source: SegSource::Blocks(self.rc_blocks),
            });
        }

        let manifest = ChunkManifest {
            chunk_index: self.chunk_index,
            segments,
            // Single-end: mate 0 only; mate 1 is unused (0).
            decoded_output_bytes: [self.decoded_output_bytes, 0],
            records: self.chunk_records,
            total_bases: self.chunk_total_bases,
            original_size: self.chunk_original_size,
        };
        manifest_tx
            .send(WorkerMsg::Manifest(manifest))
            .map_err(|_| anyhow::anyhow!("writer exited"))?;

        // Reset stagers + per-chunk accumulators/counters.
        self.seq.reset_for_next_chunk();
        if let Some(hs) = self.hdr_bsc.as_mut() {
            hs.reset_for_next_chunk();
        }
        self.qual.reset_for_next_chunk();
        if let Some(rc) = self.rc.as_mut() {
            rc.reset_for_next_chunk();
        }
        self.sequence_blocks = 0;
        self.header_blocks = 0;
        self.qual_blocks = 0;
        self.rc_blocks = 0;
        self.decoded_output_bytes = 0;
        self.chunk_records = 0;
        self.chunk_total_bases = 0;
        self.chunk_original_size = 0;
        self.chunk_index += 1;
        Ok(())
    }

    fn close_final_chunk(
        &mut self,
        job_tx: &SyncSender<CompressJob>,
        manifest_tx: &SyncSender<WorkerMsg>,
    ) -> Result<()> {
        self.close_chunk(job_tx, manifest_tx)
    }
}

/// Generic bounded worker pool. Spawns `n_workers` rayon-style threads into `scope`, each
/// pulling jobs of type `J` under `job_rx`'s mutex, running `work`, and forwarding the
/// resulting `R`s on `res_tx`. A semantic `Err` or a caught panic sets `abort` (first-error
/// -wins) and routes a single `R` built via `on_error`. The work/error closures are SHARED
/// across every scoped thread via an `Arc` (refcount clone), so they outlive this frame
/// (which returns before the caller's `scope` joins the threads) — they're `Send + Sync`.
///
/// Worker loop (identical to the old `spawn_compress_workers`): check `abort`; lock + `recv`
/// (break on disconnect); `catch_unwind(work(job))`:
///   - `Ok(Ok(results))` → send each `r`; if a send fails (consumer gone) return.
///   - `Ok(Err(e))`      → set `abort`, send `on_error(e.to_string())`, return.
///   - `Err(panic)`      → set `abort`, send `on_error("compress worker panicked: <msg>")`, return.
///
/// `on_error` receives the RAW message on the semantic-`Err` path; on the PANIC path the pool
/// itself prepends the fixed `"compress worker panicked: "` prefix (shared by every consumer, so
/// the historic single/paired panic text is byte-for-byte preserved) before calling `on_error`.
/// `on_error` just wraps whatever string it receives into the mode's error result.
/// The standard scoped-thread lifetime shape: `scope` is `&'scope Scope<'scope, 'env>`,
/// `abort` is borrowed `'env` (caller-owned), and the `Arc`-wrapped closures are `'env`.
pub(crate) fn spawn_worker_pool<'scope, 'env: 'scope, J: Send + 'env, R: Send + 'env>(
    scope: &'scope std::thread::Scope<'scope, 'env>,
    n_workers: usize,
    job_rx: &std::sync::Arc<std::sync::Mutex<std::sync::mpsc::Receiver<J>>>,
    res_tx: &std::sync::mpsc::SyncSender<R>,
    abort: &'env std::sync::atomic::AtomicBool,
    work: impl Fn(J) -> anyhow::Result<Vec<R>> + Send + Sync + 'env,
    on_error: impl Fn(String) -> R + Send + Sync + 'env,
) {
    use std::sync::atomic::Ordering;
    use std::sync::Arc;
    // `work`/`on_error` are SHARED by every scoped thread via an `Arc` (each thread owns a
    // refcount, so the closures live as long as the threads — they cannot be borrowed from
    // this frame, which returns before the caller's `scope` joins the threads). `job_rx`/
    // `res_tx` are CLONED per worker (Arc / SyncSender) — NOT borrowed for `'env` — so the
    // caller can `drop` its own clones right after spawning (the producer-unblock-on-
    // -disconnect shutdown contract). `abort` is borrowed `'env` (caller-owned, not dropped).
    let work = Arc::new(work);
    let on_error = Arc::new(on_error);
    for _ in 0..n_workers {
        let job_rx = job_rx.clone();
        let res_tx = res_tx.clone();
        let work = work.clone();
        let on_error = on_error.clone();
        scope.spawn(move || {
            loop {
                if abort.load(Ordering::Relaxed) {
                    break;
                }
                let job = {
                    let g = job_rx.lock().unwrap();
                    g.recv()
                };
                let Ok(job) = job else { break }; // job_tx dropped → done
                let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| work(job)));
                match outcome {
                    Ok(Ok(results)) => {
                        for r in results {
                            if res_tx.send(r).is_err() {
                                return;
                            }
                        }
                    }
                    Ok(Err(e)) => {
                        abort.store(true, Ordering::Relaxed);
                        let _ = res_tx.send(on_error(e.to_string()));
                        return;
                    }
                    Err(p) => {
                        // Same panic-message downcast as the old inline loop, so worker-panic
                        // text is unchanged: `compress worker panicked: <payload>`.
                        let msg = p
                            .downcast_ref::<&str>()
                            .map(|s| s.to_string())
                            .or_else(|| p.downcast_ref::<String>().cloned())
                            .unwrap_or_else(|| "<opaque panic>".to_string());
                        abort.store(true, Ordering::Relaxed);
                        let _ = res_tx.send(on_error(format!("compress worker panicked: {msg}")));
                        return;
                    }
                }
            }
        });
    }
}

/// Map a `compress_job` `JobOutput` into the writer-channel `WorkerMsg`s, in the SAME order
/// the old inline worker loop emitted them: `Blocks` → one `WorkerMsg::Block` per element,
/// `Segments` → one `WorkerMsg::Segment` per element.
fn job_output_to_worker_msgs(out: JobOutput) -> Vec<WorkerMsg> {
    match out {
        JobOutput::Blocks(results) => results.into_iter().map(WorkerMsg::Block).collect(),
        JobOutput::Segments(segs) => segs.into_iter().map(WorkerMsg::Segment).collect(),
    }
}

/// Spawn `n_workers` rayon-style threads into `scope`, each pulling `CompressJob`s under
/// `job_rx`'s mutex, compressing via `compress_job`, and forwarding `WorkerMsg`s. A worker
/// panic or semantic error is caught and routed as a first-error-wins `WorkerMsg::Err`.
/// Shared by single & paired orchestrators. A thin wrapper over the generic
/// `spawn_worker_pool`: `compress_job(_, plan)` is the work; its `JobOutput` is fanned out
/// via `job_output_to_worker_msgs`; `WorkerMsg::Err` builds the error message.
pub(crate) fn spawn_compress_workers<'scope, 'env: 'scope>(
    scope: &'scope std::thread::Scope<'scope, 'env>,
    n_workers: usize,
    job_rx: &std::sync::Arc<std::sync::Mutex<std::sync::mpsc::Receiver<CompressJob>>>,
    res_tx: &std::sync::mpsc::SyncSender<WorkerMsg>,
    plan: &'env ArchivePlan,
    abort: &'env std::sync::atomic::AtomicBool,
) {
    spawn_worker_pool(
        scope,
        n_workers,
        job_rx,
        res_tx,
        abort,
        move |job| Ok(job_output_to_worker_msgs(compress_job(job, plan)?)),
        WorkerMsg::Err,
    );
}

fn compress_chunked_streaming(args: &CompressConfig, mode: ChunkedMode) -> Result<()> {
    // Clamp to >= 1 so a library caller passing `chunk_records == 0` reads one record per
    // chunk instead of producing an empty first chunk and a misleading "Empty input file"
    // error (the paired producer already clamps; this matches it). The CLI sets a sane default.
    let chunk_size: usize = match mode {
        ChunkedMode::Ultra { chunk_records, .. } => chunk_records,
        ChunkedMode::Default => args.advanced.chunk_records,
    }
    .max(1);
    let mut reader = crate::io::FastqReader::from_path_or_stdin(&args.input[0], args.fasta)?;
    let (first_chunk, _b, _o) = read_chunk_records(&mut reader, chunk_size)?;
    if first_chunk.is_empty() {
        anyhow::bail!("Empty input file");
    }
    let plan = resolve_archive_plan(args, mode, &first_chunk)?;
    run_streaming_compress(args, mode, plan, first_chunk, reader, false)
}

/// Shared streaming-compress body used by BOTH the whole-file path
/// (`compress_chunked_streaming`) and the NUMA byte-range path
/// (`compress_byte_range`). Takes an already-resolved `plan` and the already-read
/// `first_chunk` + the live `reader` (positioned just after it).
///
/// `validate_first_chunk` is the spec §5.1.3 switch: the whole-file path passes
/// `false` (its first chunk DEFINED the const lengths), the NUMA path passes `true`
/// (its first chunk came from the parent's head — a different record set — so its
/// lengths must be checked against the pinned const lengths or const-length framing
/// would silently mis-slice on decode).
pub(super) fn run_streaming_compress<R: std::io::BufRead>(
    args: &CompressConfig,
    mode: ChunkedMode,
    plan: ArchivePlan,
    first_chunk: Vec<crate::io::FastqRecord>,
    mut reader: crate::io::FastqReader<R>,
    validate_first_chunk: bool,
) -> Result<()> {
    use std::io::Write;
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::mpsc::sync_channel;
    use std::sync::{Arc, Mutex};

    let start_time = Instant::now();
    let chunk_size: usize = match mode {
        ChunkedMode::Ultra { chunk_records, .. } => chunk_records,
        ChunkedMode::Default => args.advanced.chunk_records,
    };
    let bsc_block_mb = match mode {
        ChunkedMode::Ultra { bsc_block_mb, .. } => bsc_block_mb,
        ChunkedMode::Default => args.advanced.bsc_block_size_mb,
    };
    super::bsc::check_cuda_vram(bsc_block_mb * 1024 * 1024);
    info!("Chunked compression (streaming): {} records per chunk", chunk_size);

    if plan.const_seq_len > 0 {
        info!("Constant sequence length detected: {} bp", plan.const_seq_len);
    }
    if plan.const_qual_len > 0 {
        info!("Constant quality length detected: {} bp", plan.const_qual_len);
    }

    // Atomic output: stream to a sibling O_EXCL `.tmp` and rename on success, so an
    // aborted compress (producer/worker error, panic) never leaves a usable partial
    // archive at the final path — the drop guard removes the temp on any early return.
    // Stdout bypasses this. Mirrors the decode side (decompress_impl).
    let is_stdout = crate::cli::is_stdio_path(&args.output);
    struct TmpCleanup(std::path::PathBuf);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    let tmp_output = if is_stdout {
        None
    } else {
        let mut name = args
            .output
            .file_name()
            .ok_or_else(|| anyhow::anyhow!("output path has no file name: {:?}", args.output))?
            .to_os_string();
        name.push(super::atomic_tmp_suffix());
        Some(args.output.with_file_name(name))
    };
    let tmp_guard = tmp_output.as_ref().map(|p| TmpCleanup(p.clone()));
    let write_path = tmp_output.as_deref().unwrap_or(args.output.as_path());

    // Front header (byte-identical to legacy). `Send` writer handle so it can move into
    // the writer thread (use the Send `Stdout`, NOT the !Send `stdout().lock()` guard).
    let mut output_file: Box<dyn Write + Send> = if is_stdout {
        Box::new(std::io::BufWriter::with_capacity(4 << 20, std::io::stdout()))
    } else {
        Box::new(std::io::BufWriter::new(super::create_new_for_write(write_path)?))
    };
    // const_seq_len/const_qual_len are usize — write_v5_archive_header takes usize (NO `as u32`).
    let header_size = super::write_v5_archive_header(
        &mut output_file,
        plan.encoding_type,
        plan.stream_quality_binning,
        plan.quality_compressor_used,
        plan.header_compressor,
        plan.fasta,
        plan.const_seq_len,
        plan.const_qual_len,
        0, // archive_type: single-end
    )?;

    // Codecs now travel in the manifest segments (derived in `ProducerState::close_chunk`),
    // so the orchestrator no longer needs per-codec locals for the writer.
    let target = plan.bsc_block_mb * 1024 * 1024;

    // ── Worker pool follows the thread setting (`-t`) — the user's parallelism / RAM
    // dial — capped by a memory safety budget so a huge `-t` can't OOM. Peak RSS scales
    // with the pool's working set (≈ workers × block × ~7 libsais workspace) and is
    // constant in read count. Ultra runs ONE block at a time (big-block internal MT-BWT;
    // concurrent forward-BWT blows up libsais memory). fqz quality blocks are larger than
    // a BSC block but FEW (one per sub-chunk, ~1:5 with BSC blocks), so their raw hold is
    // bounded by `queue_depth` and does not dominate RSS. `QZ_COMPRESS_WORKERS` overrides
    // the count entirely (tuning/benchmarking). ──
    let is_ultra = matches!(mode, ChunkedMode::Ultra { .. });
    let n_workers = std::env::var("QZ_COMPRESS_WORKERS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n >= 1)
        .unwrap_or_else(|| {
            compress_worker_count(rayon::current_num_threads().max(1), target, is_ultra)
        });
    // Ultra MUST run one block at a time even when QZ_COMPRESS_WORKERS is set: ultra's
    // big blocks parallelise internally via libbsc's MT-BWT, and running concurrent
    // forward-BWTs blows up libsais memory (see finding_ultra_mt_bwt_big_blocks). The env
    // override is a tuning dial for the default 25 MB-block path only.
    let n_workers = if is_ultra { 1 } else { n_workers };
    // queue_depth == n_workers (≤ queue_depth queued + n_workers in-flight raw blocks).
    let queue_depth = n_workers;
    let (job_tx, job_rx) = sync_channel::<CompressJob>(queue_depth);
    let (res_tx, res_rx) = sync_channel::<WorkerMsg>(queue_depth.max(1) * 2);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let abort = Arc::new(AtomicBool::new(false));

    std::thread::scope(|scope| -> Result<()> {
        // ── Workers: pull → compress → emit. On error: set abort + send Err + exit. ──
        // `&*abort` derefs the `Arc<AtomicBool>` to the `&AtomicBool` the pool borrows for
        // `'env`; the Arc itself is still owned here for the producer's own `abort` checks.
        spawn_compress_workers(scope, n_workers, &job_rx, &res_tx, &plan, &abort);

        // CORRECTION 1: drop the main thread's `job_rx` clone NOW so the ONLY receivers
        // live inside the workers. Otherwise, on a worker error the surviving workers see
        // `abort` and exit WITHOUT draining the queue; the producer, blocked on a full
        // `job_tx.send()`, would never reach its next abort check, and because the main
        // thread still held a `job_rx` the channel would never disconnect → permanent
        // hang. With this drop, once all workers exit the last receiver drops, the channel
        // disconnects, and the producer's blocked send returns Err (→ "compress workers
        // exited", a harmless secondary error — the writer keeps the FIRST error).
        drop(job_rx);

        // Producer's OWN results-channel sender for manifests. Dropping the original
        // `res_tx` here means `res_rx` closes only once every worker AND the producer have
        // dropped their clone — no premature close, no hang.
        let manifest_tx = res_tx.clone();
        drop(res_tx);

        // ── Writer thread: sole error aggregator. Returns Err if it saw any
        // WorkerMsg::Err (producer or worker root cause), else finalizes Ok(WriterDone).
        // Codecs now travel in the manifest segments, so the writer needs none of the
        // per-codec locals — only the mode + header size. ──
        let writer = scope.spawn(move || -> Result<WriterDone> {
            write_ordered(output_file, CompressMode::Single, header_size, res_rx)
        });

        // ── Producer body (this thread). On its OWN error: funnel to the writer as a
        // WorkerMsg::Err (the writer is the single decision point), set abort, stop.
        // Check `abort` each iteration so a worker error halts production promptly. ──
        let prod_res: Result<()> = (|| {
            let mut state = ProducerState::new(&plan, target);
            // Whole-file path (`validate_first_chunk == false`): the first_chunk records
            // DEFINED the const lengths → NOT re-validated. NUMA byte-range path
            // (`validate_first_chunk == true`): the first chunk came from THIS shard's
            // head but the const lengths were pinned by the PARENT (a different record
            // set), so every record must be checked or const-length framing would
            // silently mis-slice on decode (spec §5.1.3).
            for rec in &first_chunk {
                if abort.load(Ordering::Relaxed) {
                    return Ok(());
                }
                if validate_first_chunk {
                    validate_const_one(
                        rec,
                        plan.const_seq_len,
                        plan.const_qual_len,
                        plan.no_quality,
                        state.chunk_index,
                    )?;
                }
                state.stage_one(rec, &job_tx)?;
                state.maybe_close_chunk(chunk_size, &job_tx, &manifest_tx)?;
            }
            drop(first_chunk);
            while let Some(rec) = reader.next()? {
                if abort.load(Ordering::Relaxed) {
                    return Ok(());
                }
                validate_const_one(
                    &rec,
                    plan.const_seq_len,
                    plan.const_qual_len,
                    plan.no_quality,
                    state.chunk_index,
                )?;
                state.stage_one(&rec, &job_tx)?;
                state.maybe_close_chunk(chunk_size, &job_tx, &manifest_tx)?;
            }
            if !abort.load(Ordering::Relaxed) {
                state.close_final_chunk(&job_tx, &manifest_tx)?;
            }
            Ok(())
        })();
        if let Err(e) = &prod_res {
            abort.store(true, Ordering::Relaxed);
            let _ = manifest_tx.send(WorkerMsg::Err(e.to_string())); // writer records the first error
        }
        drop(job_tx); // workers drain remaining jobs and exit
        drop(manifest_tx); // last non-worker sender gone → res_rx closes after workers finish

        // The writer's Result is authoritative (it carries producer/worker root causes
        // via WorkerMsg::Err). Don't `prod_res?` — that would mask the real error with a
        // secondary "workers exited" send error.
        let done = writer
            .join()
            .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;
        log_compression_stats(
            done.original_size,
            done.headers_len,
            done.sequences_len,
            done.qualities_len,
            done.rc_flags_len,
            done.metadata_size,
            start_time.elapsed(),
        );
        Ok(())
    })?; // any producer/worker/writer error: tmp_guard drops here → temp removed

    // Success: disarm the cleanup guard and atomically publish the temp file.
    if let (Some(tmp), Some(guard)) = (tmp_output.as_ref(), tmp_guard) {
        std::mem::forget(guard);
        if let Err(e) = std::fs::rename(tmp, &args.output) {
            let _ = std::fs::remove_file(tmp);
            anyhow::bail!("Failed to rename temp file to {:?}: {}", args.output, e);
        }
    }
    Ok(())
}

/// The writer thread: the SOLE error aggregator. Drains `res_rx`, orders blocks by
/// chunk via `OrderedBlockBuffer`, and writes the per-chunk segments per the active
/// `CompressMode` (§6): Single = per-block entries (Headers → Sequence → Qual → RcFlags)
/// plus the ChunkDecodedSizes global; Paired = per-(mate,role) `write_block_stream` blobs,
/// no ChunkDecodedSizes, `validate_paired_directory` before the footer. First
/// WorkerMsg::Err wins.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_ordered(
    mut output_file: Box<dyn std::io::Write + Send>,
    mode: CompressMode,
    header_size: usize,
    res_rx: std::sync::mpsc::Receiver<WorkerMsg>,
) -> Result<WriterDone> {
    use crate::compression::chunk_directory::{
        ChunkDirEntry, ChunkDirectory, GLOBAL_SENTINEL, StreamRole, encode_decoded_sizes,
        write_footer_and_locator,
    };
    use std::io::Write;

    let mut cursor: u64 = header_size as u64;
    let mut dir: Vec<ChunkDirEntry> = Vec::new();
    let mut decoded_sizes: Vec<u64> = Vec::new();
    let mut buf = OrderedBlockBuffer::new();
    let mut num_chunks: u32 = 0;
    let mut num_reads: u64 = 0;
    let mut original_size: u64 = 0;
    let mut total_bases: u64 = 0; // log-only; does NOT affect archive bytes

    for msg in res_rx {
        match msg {
            WorkerMsg::Manifest(m) => {
                num_chunks = num_chunks.max(m.chunk_index + 1);
                num_reads += m.records;
                original_size += m.original_size;
                total_bases += m.total_bases;
                buf.set_manifest(m);
            }
            WorkerMsg::Block(r) => {
                buf.add_block(r)?;
            }
            WorkerMsg::Segment(s) => {
                buf.add_header_segment(s.chunk_index, s.mate, s.role, s.codec, s.blocks)?;
            }
            WorkerMsg::Err(e) => {
                // First-error-wins. Dropping `res_rx` (by returning) makes any blocked
                // worker/producer send error out so they exit; the scope still joins them.
                return Err(anyhow::anyhow!(e));
            }
        }
        for cc in buf.drain_ready() {
            match mode {
                CompressMode::Single => {
                    for seg in &cc.segments {
                        write_role_blocks(
                            &seg.blocks,
                            seg.role,
                            seg.mate,
                            seg.codec,
                            cc.chunk_index,
                            &mut output_file,
                            &mut cursor,
                            &mut dir,
                        )?;
                    }
                    decoded_sizes.push(cc.decoded_output_bytes[0]);
                }
                CompressMode::Paired => {
                    for seg in &cc.segments {
                        write_segment_blob(
                            &seg.blocks,
                            seg.role,
                            seg.mate,
                            seg.codec,
                            cc.records,
                            cc.chunk_index,
                            &mut output_file,
                            &mut cursor,
                            &mut dir,
                        )?;
                    }
                    // Per-mate decoded sizes, row-major (sizes[chunk*2 + mate]), for the
                    // ChunkDecodedSizes(n_mates=2) global that feeds paired direct-write.
                    decoded_sizes.push(cc.decoded_output_bytes[0]);
                    decoded_sizes.push(cc.decoded_output_bytes[1]);
                }
            }
        }
    }

    // Clean close (all senders dropped, no Err seen). Assert the buffer fully drained,
    // converting a count-mismatch hang/truncation into a clean error.
    if !buf.all_emitted(num_chunks) {
        anyhow::bail!(
            "streaming writer: ordered buffer did not fully drain ({} chunks expected) — \
             a block was lost or mis-routed",
            num_chunks
        );
    }
    info!(
        "Read {} records in {} chunks ({} bases)",
        num_reads, num_chunks, total_bases
    );

    // Per-role compressed totals (frame bytes) for the stats line, computed BEFORE
    // moving `dir` into the footer struct.
    let role_bytes = |want: StreamRole| -> usize {
        dir.iter()
            .filter(|e| e.role == want)
            .map(|e| e.length as usize)
            .sum()
    };
    let headers_len = role_bytes(StreamRole::Headers);
    let sequences_len = role_bytes(StreamRole::Sequence);
    let qualities_len = role_bytes(StreamRole::Qual);
    let rc_flags_len = role_bytes(StreamRole::RcFlags);

    // ChunkDecodedSizes global: single-end emits n_mates=1 (one size/chunk), paired emits
    // n_mates=2 (R1,R2 row-major). Both feed the NUMA direct-write decode pre-sizing. The
    // block's record_count is num_chunks (NOT the size count) for both — `read_chunk_layout`
    // cross-checks it against the footer chunk count. Single's bytes are unchanged: for
    // single decoded_sizes.len() == num_chunks, so this is byte-identical to the old emit.
    {
        let n_mates: u8 = if mode == CompressMode::Paired { 2 } else { 1 };
        debug_assert_eq!(decoded_sizes.len(), num_chunks as usize * n_mates as usize);
        // Debug-only test hook: undersize the last entry so direct-write decode exercises
        // the overrun → integrity → auto-fallback path (single and paired alike).
        #[cfg(debug_assertions)]
        if std::env::var_os("QZ_NUMA_FAKE_UNDERSIZE").is_some()
            && let Some(last) = decoded_sizes.last_mut()
        {
            *last = last.saturating_sub(1);
        }
        let dsz = encode_decoded_sizes(n_mates, num_chunks, &decoded_sizes);
        write_role_blocks(
            &[(num_chunks, dsz)],
            StreamRole::ChunkDecodedSizes,
            0, // mate (global)
            0, // codec (raw metadata)
            GLOBAL_SENTINEL,
            &mut output_file,
            &mut cursor,
            &mut dir,
        )?;
    }

    if mode == CompressMode::Paired {
        let directory = ChunkDirectory {
            num_reads,
            num_chunks,
            entries: dir,
        };
        crate::compression::paired::format::validate_paired_directory(&directory)?;
        let footer_len = write_footer_and_locator(output_file.as_mut(), &directory)?;
        output_file.flush()?;
        let metadata_size = header_size + footer_len + 20;
        return Ok(WriterDone {
            original_size: original_size as usize,
            headers_len,
            sequences_len,
            qualities_len,
            rc_flags_len,
            metadata_size,
        });
    }

    let directory = ChunkDirectory {
        num_reads,
        num_chunks,
        entries: dir,
    };
    let footer_len = write_footer_and_locator(output_file.as_mut(), &directory)?;
    output_file.flush()?;

    let metadata_size = header_size + footer_len + 20;
    Ok(WriterDone {
        original_size: original_size as usize,
        headers_len,
        sequences_len,
        qualities_len,
        rc_flags_len,
        metadata_size,
    })
}

/// Reject configurations that would silently produce an unreadable/corrupt
/// archive, or silently ignore a requested option. This is the single home for
/// every cross-option compression guard — adding a new constraint means adding
/// it here (and a test), rather than threading another `bail!` into the middle
/// of `compress()`. Pure: depends only on `args`, runs before any work.
fn validate_compress_config(args: &CompressConfig) -> Result<()> {
    // Cap bsc_block_size_mb. The default streaming decompressor refuses any
    // compressed block > 64 MiB (corruption guard), and BSC output never
    // exceeds input + header — so capping input at 64 MiB guarantees the
    // archive is always decompressible. The ultra (UltraBigBlock) decode path
    // lifts that to 768 MiB, so the ultra *encoder* input is capped at 750 MiB
    // — strictly below the decode cap, leaving margin for an incompressible
    // block whose BSC output exceeds its input.
    let is_ultra_path = args.ultra.is_some();
    let max_block_mb = if is_ultra_path { 750 } else { 64 };
    if args.advanced.bsc_block_size_mb == 0 || args.advanced.bsc_block_size_mb > max_block_mb {
        anyhow::bail!(
            "bsc_block_size_mb must be in 1..={} (got {})",
            max_block_mb,
            args.advanced.bsc_block_size_mb
        );
    }

    // Validate rc_canon compatibility.
    //
    // The encoding-type selection is a priority chain (ultra > rc_canon >
    // hints), so rc_canon wins and the archive is tagged `RcCanon` (type 6).
    // But `records_to_streams` emits the hint byte whenever `sequence_hints` is
    // set, INDEPENDENT of the chosen type — and the type-6 decoder reports
    // `has_sequence_hints()==false`, so it never strips that extra byte and reads
    // it as sequence data. The combo therefore silently corrupts every sequence.
    // Reject it.
    if args.advanced.rc_canon && args.advanced.sequence_hints {
        anyhow::bail!(
            "--rc-canon cannot be combined with --sequence-hints \
             (it would inject bytes the rc-canon decoder does not strip)"
        );
    }

    // Validate --quality-compressor fqzcomp.
    // The fqz encode path feeds raw quality ASCII to fqzcomp and forces
    // `stream_quality_binning = None` (it has no lossy/binned path), so an explicit
    // Fqzcomp + a lossy quality mode would SILENTLY store full-resolution qualities
    // — ignoring the requested reduction — and an explicit Fqzcomp + --no-quality
    // ignores the codec choice entirely (no quality stream is written). Reject both
    // rather than silently honor only part of the request. (Unreachable in
    // production: the CLI and Python binding leave quality_compressor = Auto, and
    // Auto resolves lossy / no-quality to BSC. This guards direct library callers.)
    if args.advanced.quality_compressor == QualityCompressor::Fqzcomp {
        if args.no_quality {
            anyhow::bail!("--quality-compressor fqzcomp cannot be used with --no-quality");
        }
        if args.quality_mode != QualityMode::Lossless {
            anyhow::bail!(
                "--quality-compressor fqzcomp requires lossless quality mode (default); \
                 it has no lossy/binned path and would silently store full-resolution qualities"
            );
        }
    }

    // Validate --ultra options.
    let has_ultra = args.ultra.is_some();
    if has_ultra {
        let mode_name = "--ultra";
        if args.advanced.sequence_hints || args.advanced.rc_canon {
            anyhow::bail!("{mode_name} cannot be combined with other encoding options");
        }
    }

    if args.input.len() != 1 {
        anyhow::bail!("Expected exactly 1 input file, got {}", args.input.len());
    }

    Ok(())
}

pub(super) fn compress(args: &CompressConfig) -> Result<()> {
    // Reject invalid / silently-lossy config combinations before any work.
    validate_compress_config(args)?;

    // Honor -t. Idempotent: the top-level `compress()` already built the global
    // pool before dispatch; this also covers direct callers (e.g. the stdout
    // path) and keeps `actual_threads` for logging. See
    // `super::ensure_global_thread_pool` for why the global pool (not an owned
    // pool + install()) is the right mechanism here.
    let actual_threads = super::ensure_global_thread_pool(args.threads);

    info!("Using {} threads for compression", actual_threads);
    info!("Input files: {:?}", args.input);
    info!("Output: {:?}", args.output);

    // Warn loudly when the user opted into lossy quality compression. Silent
    // data loss is the worst kind — make it visible in the log so a downstream
    // pipeline can spot the choice in archived runs.
    match args.quality_mode {
        crate::cli::QualityMode::Lossless => {}
        crate::cli::QualityMode::Discard => {
            tracing::warn!(
                "Lossy quality mode 'discard': all quality scores will be replaced with a single \
                 placeholder value. Original qualities will NOT be recoverable from the archive."
            );
        }
        other => {
            tracing::warn!(
                "Lossy quality mode {:?}: quality scores will be re-binned. Original qualities \
                 will NOT be recoverable from the archive.",
                other,
            );
        }
    }
    if args.no_quality {
        tracing::warn!(
            "Quality scores disabled (--no-quality): the archive will not contain quality data."
        );
    }

    // Ultra mode: level-based compression with auto-tuning
    if let Some(requested_level) = args.ultra {
        let level = ultra::resolve_ultra_level(requested_level)?;
        return ultra::compress_ultra(args, level);
    }

    // Single archive architecture: the unified chunked v5 pipeline handles every
    // supported configuration (BSC/Columnar headers, BSC sequences, BSC/fqzcomp
    // quality). `validate_compress_config` has already rejected any unsupported
    // option, so there is no in-memory fallback.
    compress_chunked(args, ChunkedMode::Default)
}

#[cfg(test)]
mod config_guard_tests {
    use super::*;
    use crate::cli::{CompressConfig, QualityCompressor, QualityMode};
    use std::path::PathBuf;

    fn base() -> CompressConfig {
        CompressConfig {
            input: vec![PathBuf::from("in.fastq")],
            ..CompressConfig::default()
        }
    }

    #[test]
    fn explicit_fqzcomp_lossless_is_accepted() {
        let mut a = base();
        a.advanced.quality_compressor = QualityCompressor::Fqzcomp;
        a.quality_mode = QualityMode::Lossless;
        assert!(validate_compress_config(&a).is_ok());
    }

    #[test]
    fn explicit_fqzcomp_lossy_is_rejected() {
        // The fqz encode path forces stream_quality_binning = None, so a lossy mode
        // would silently store full-resolution qualities — reject, never ignore.
        for mode in [
            QualityMode::IlluminaBin,
            QualityMode::Illumina4,
            QualityMode::Binary,
        ] {
            let mut a = base();
            a.advanced.quality_compressor = QualityCompressor::Fqzcomp;
            a.quality_mode = mode;
            assert!(
                validate_compress_config(&a).is_err(),
                "explicit fqzcomp + {mode:?} must be rejected"
            );
        }
    }

    #[test]
    fn explicit_fqzcomp_no_quality_is_rejected() {
        let mut a = base();
        a.advanced.quality_compressor = QualityCompressor::Fqzcomp;
        a.no_quality = true;
        assert!(validate_compress_config(&a).is_err());
    }
}

#[cfg(test)]
mod fqzcomp_subblock_tests {
    use super::*;

    /// Distinct per-record quality strings (the byte-capped splitter used by
    /// reference mode operates on borrowed quality slices).
    fn sample_quals(n: u8) -> Vec<Vec<u8>> {
        (0..n)
            .map(|i| vec![b'I', b'5', b'#', b'0' + (i % 9)])
            .collect()
    }

    #[test]
    fn fqz_blocks_capped_quals_enforces_cap() {
        let quals = sample_quals(64);
        let refs: Vec<&[u8]> = quals.iter().map(|q| q.as_slice()).collect();
        // Largest single-record blob: cap == this guarantees recursion terminates
        // (each leaf is one record and fits), and forces splitting iff full > cap.
        let max_single = refs
            .iter()
            .map(|q| {
                codecs::compress_qualities_fqzcomp_quals(std::slice::from_ref(q))
                    .unwrap()
                    .len()
            })
            .max()
            .unwrap();
        let full = codecs::compress_qualities_fqzcomp_quals(&refs).unwrap().len();
        assert!(
            full > max_single,
            "64-record blob must exceed a single-record blob"
        );
        let cap = max_single;

        let blocks = fqz_blocks_capped_quals(&refs, cap).unwrap();
        assert!(
            blocks.len() > 1,
            "cap == single-record size must force splitting"
        );
        assert!(
            blocks.iter().all(|(_, b)| b.len() <= cap),
            "every block within cap"
        );
        assert_eq!(
            blocks.iter().map(|(c, _)| *c).sum::<u32>(),
            64,
            "counts sum to input"
        );
    }

    #[test]
    fn fqz_blocks_capped_quals_errors_below_single_record() {
        let quals = sample_quals(8);
        let refs: Vec<&[u8]> = quals.iter().map(|q| q.as_slice()).collect();
        let min_single = refs
            .iter()
            .map(|q| {
                codecs::compress_qualities_fqzcomp_quals(std::slice::from_ref(q))
                    .unwrap()
                    .len()
            })
            .min()
            .unwrap();
        // A cap below the smallest single-record blob cannot be satisfied → Err, not panic.
        let res = fqz_blocks_capped_quals(&refs, min_single.saturating_sub(1));
        assert!(res.is_err(), "sub-single-record cap must return Err");
    }
}

#[cfg(test)]
mod columnar_cap_tests {
    use super::*;

    /// Distinct, mildly-compressible read IDs (grow the columnar blob with count).
    fn sample_headers(n: usize) -> Vec<Vec<u8>> {
        (0..n)
            .map(|i| format!("@SRR{:08}.{} read_{}/1", i, i * 7 + 3, i).into_bytes())
            .collect()
    }

    #[test]
    fn columnar_blocks_capped_single_block_under_cap() {
        let hdrs = sample_headers(64);
        let refs: Vec<&[u8]> = hdrs.iter().map(|h| h.as_slice()).collect();
        // Generous cap → exactly one block, byte-identical to the pre-cap encoder
        // (`[count u32][columnar_blob]`). This guarantees existing archives are unchanged.
        let blocks = columnar_blocks_capped(&refs, usize::MAX).unwrap();
        assert_eq!(blocks.len(), 1, "fits under cap → single block");
        let blob = header_col::compress_headers_columnar_bytes(&refs).unwrap();
        let mut expect = (refs.len() as u32).to_le_bytes().to_vec();
        expect.extend_from_slice(&blob);
        assert_eq!(blocks[0].0, 64, "framed record_count == header count");
        assert_eq!(
            blocks[0].1, expect,
            "single-block payload byte-identical to old path"
        );
    }

    #[test]
    fn columnar_blocks_capped_splits_over_cap_and_reassembles() {
        let hdrs = sample_headers(256);
        let refs: Vec<&[u8]> = hdrs.iter().map(|h| h.as_slice()).collect();
        let full_len = columnar_blocks_capped(&refs, usize::MAX).unwrap()[0]
            .1
            .len();
        // Cap forces multiple blocks but stays well above any single-header payload,
        // so recursion terminates without hitting the single-record bail.
        let cap = full_len / 3;
        let blocks = columnar_blocks_capped(&refs, cap).unwrap();
        assert!(blocks.len() > 1, "cap below full size must split");
        assert!(
            blocks.iter().all(|(_, p)| p.len() <= cap),
            "every emitted block within cap"
        );
        assert_eq!(
            blocks.iter().map(|(c, _)| *c as usize).sum::<usize>(),
            256,
            "framed counts sum to input"
        );
        // Each block is a self-contained columnar encoding; decoding each and
        // concatenating in block order reconstructs every id in input order.
        let mut got: Vec<Vec<u8>> = Vec::new();
        for (rc, payload) in &blocks {
            let count = u32::from_le_bytes(payload[0..4].try_into().unwrap());
            assert_eq!(count, *rc, "payload count prefix == framed record_count");
            let ids =
                header_col::decompress_headers_columnar(&payload[4..], count as usize).unwrap();
            got.extend(ids);
        }
        assert_eq!(got, hdrs, "reassembled ids match input");
    }

    #[test]
    fn columnar_blocks_capped_errors_on_unsatisfiable_cap() {
        let hdrs = sample_headers(4);
        let refs: Vec<&[u8]> = hdrs.iter().map(|h| h.as_slice()).collect();
        // A cap below even a single-header payload cannot be satisfied → Err, not
        // a panic or unbounded recursion.
        assert!(
            columnar_blocks_capped(&refs, 1).is_err(),
            "sub-single-header cap must Err"
        );
    }
}

#[cfg(test)]
mod decoded_sizes_tests {
    #[test]
    fn single_end_archive_carries_decoded_sizes_global() {
        use crate::compression::chunk_directory::{
            GLOBAL_SENTINEL, StreamRole, parse_decoded_sizes, read_v5_footer,
        };
        let d = tempfile::tempdir().unwrap();
        let fq = d.path().join("in.fastq");
        let mut s = String::new();
        for i in 0..5 {
            s.push_str(&format!("@read{i}\nACGTACGT\n+\nIIIIIIII\n"));
        }
        std::fs::write(&fq, &s).unwrap();
        let arc = d.path().join("out.qz");
        let mut cfg = crate::cli::CompressConfig {
            input: vec![fq],
            output: arc.clone(),
            ..Default::default()
        };
        cfg.advanced.chunk_records = 2; // usize on AdvancedOptions -> 3 chunks
        crate::compression::compress(&cfg).unwrap();

        let header_end = crate::compression::test_front_header_end(&arc).unwrap(); // returns stored header_size
        let footer = read_v5_footer(&arc, header_end).unwrap();
        let g = footer
            .entries
            .iter()
            .find(|e| e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::ChunkDecodedSizes)
            .expect("ChunkDecodedSizes global present");
        assert_eq!(g.record_count, footer.num_chunks as u64);
        assert_eq!(g.mate, 0);
        let bytes = std::fs::read(&arc).unwrap();
        // +12 skips the per-block frame header: [block_len u32 | record_count u32 | crc32 u32].
        let payload = &bytes[(g.offset as usize + 12)..(g.offset as usize + g.length as usize)];
        let (n_mates, num_chunks, sizes) = parse_decoded_sizes(payload).unwrap();
        assert_eq!((n_mates, num_chunks), (1, footer.num_chunks));
        assert_eq!(sizes.iter().sum::<u64>(), 5 * (6 + 8 + 8 + 5)); // "@readN"=6, seq8, qual8, +5
    }

    /// Flipping a byte inside the ChunkDecodedSizes payload must cause
    /// `read_chunk_layout` to return `Err` (the frame CRC check fires).
    #[test]
    fn corrupt_decoded_sizes_table_is_rejected() {
        use crate::compression::chunk_directory::{
            GLOBAL_SENTINEL, StreamRole, read_v5_footer,
        };
        let d = tempfile::tempdir().unwrap();
        let fq = d.path().join("in.fastq");
        let mut s = String::new();
        for i in 0..5 {
            s.push_str(&format!("@read{i}\nACGTACGT\n+\nIIIIIIII\n"));
        }
        std::fs::write(&fq, &s).unwrap();
        let arc = d.path().join("out.qz");
        let mut cfg = crate::cli::CompressConfig {
            input: vec![fq],
            output: arc.clone(),
            ..Default::default()
        };
        cfg.advanced.chunk_records = 2; // 3 chunks
        crate::compression::compress(&cfg).unwrap();

        // Locate the ChunkDecodedSizes global entry via the footer.
        let header_end = crate::compression::test_front_header_end(&arc).unwrap();
        let footer = read_v5_footer(&arc, header_end).unwrap();
        let g = footer
            .entries
            .iter()
            .find(|e| e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::ChunkDecodedSizes)
            .expect("ChunkDecodedSizes global present");

        // Flip one byte inside the payload region (after the 12-byte frame header).
        let payload_start = g.offset as usize + 12;
        let mut bytes = std::fs::read(&arc).unwrap();
        assert!(
            payload_start < g.offset as usize + g.length as usize,
            "payload region must be non-empty"
        );
        bytes[payload_start] ^= 0xFF;
        let corrupted = d.path().join("corrupted.qz");
        std::fs::write(&corrupted, &bytes).unwrap();

        // read_chunk_layout must reject the corrupted archive (CRC mismatch).
        let result = crate::compression::read_chunk_layout(&corrupted);
        assert!(result.is_err(), "corrupted ChunkDecodedSizes table must be rejected");
        let err_msg = result.unwrap_err().to_string().to_lowercase();
        assert!(
            err_msg.contains("crc"),
            "error should mention CRC, got: {err_msg}"
        );
    }

    /// `verify --fast` must CRC-walk the ChunkDecodedSizes global block: flipping a
    /// byte inside its payload must make fast verify return Err, while the un-corrupted
    /// control verifies clean. (Without the global-block walk, the corruption slips
    /// through fast verify — it only ever affected the decode-time table read.)
    #[test]
    fn verify_fast_crc_walks_decoded_sizes_global() {
        use crate::compression::chunk_directory::{
            GLOBAL_SENTINEL, StreamRole, read_v5_footer,
        };
        let d = tempfile::tempdir().unwrap();
        let fq = d.path().join("in.fastq");
        let mut s = String::new();
        for i in 0..6 {
            s.push_str(&format!("@read{i}\nACGTACGT\n+\nIIIIIIII\n"));
        }
        std::fs::write(&fq, &s).unwrap();
        let arc = d.path().join("out.qz");
        let mut cfg = crate::cli::CompressConfig {
            input: vec![fq],
            output: arc.clone(),
            ..Default::default()
        };
        cfg.advanced.chunk_records = 2; // 3 chunks → multi-chunk single-end
        crate::compression::compress(&cfg).unwrap();

        let verify_fast = |path: &std::path::Path| {
            crate::compression::verify(&crate::cli::VerifyConfig {
                input: path.to_path_buf(),
                working_dir: d.path().to_path_buf(),
                num_threads: 1,
                fast: true,
            })
        };

        // Control: the un-corrupted archive must pass fast verify.
        verify_fast(&arc).expect("control archive must pass fast verify");

        // Locate the ChunkDecodedSizes global entry via the footer and flip one byte
        // inside its payload region (after the 12-byte frame header).
        let header_end = crate::compression::test_front_header_end(&arc).unwrap();
        let footer = read_v5_footer(&arc, header_end).unwrap();
        let g = footer
            .entries
            .iter()
            .find(|e| e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::ChunkDecodedSizes)
            .expect("ChunkDecodedSizes global present");
        let payload_start = g.offset as usize + 12;
        assert!(
            payload_start < g.offset as usize + g.length as usize,
            "payload region must be non-empty"
        );
        let mut bytes = std::fs::read(&arc).unwrap();
        bytes[payload_start] ^= 0xFF;
        let corrupted = d.path().join("corrupted.qz");
        std::fs::write(&corrupted, &bytes).unwrap();

        // Corrupted: fast verify must now fail on the ChunkDecodedSizes CRC.
        let result = verify_fast(&corrupted);
        assert!(
            result.is_err(),
            "corrupted ChunkDecodedSizes table must fail fast verify"
        );
        let err_msg = result.unwrap_err().to_string().to_lowercase();
        assert!(
            err_msg.contains("crc"),
            "error should mention CRC, got: {err_msg}"
        );
    }
}
