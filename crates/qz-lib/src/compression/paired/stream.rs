//! Block-granular streaming paired-end compress (spec Increment 2). Two mate readers are
//! zipped per pair; each mate frames into its own stagers (per-mate ArchivePlan). Headers
//! for BOTH mates are compressed by ONE `PairedHeaders` job per chunk (which makes the R2
//! smaller-of decision). Seq/quality stream block-granular per mate. The writer emits the
//! 6-segment per-chunk layout (CompressMode::Paired) — byte-identical to compress_paired_v5.

use anyhow::{Result, bail};
use std::io::Write;
use std::sync::mpsc::SyncSender;

use crate::cli::CompressConfig;
use crate::compression::chunk_directory::StreamRole;
use crate::compression::codec_ids::{CODEC_BSC, CODEC_COLUMNAR};
use crate::compression::compress_impl::{CompressMode, WorkerMsg, spawn_compress_workers, write_ordered};
use crate::compression::stream_compress::{
    ArchivePlan, ChunkManifest, CompressJob, HeaderColStager, PairedHeaderJob, QualStager,
    SeqQualHeaderStager, SegSource, SegmentSpec, compress_worker_count,
};
use crate::io::FastqRecord;
use super::resolve_mate_quality;

/// Per-mate sequence + quality stagers (+ the chunk's header ids). Headers are NOT
/// dispatched per mate — the producer emits a combined `PairedHeaders` job at chunk close.
struct MateStager {
    mate: u8,
    plan: ArchivePlan,
    qual_codec: u8,
    seq: SeqQualHeaderStager,
    qual: QualStager,
    hdr: HeaderColStager,
    fh: Vec<u8>, fs: Vec<u8>, fq: Vec<u8>, frc: Vec<u8>, // frc unused (rc_canon off for paired)
    seq_blocks: u32,
    qual_blocks: u32,
}

impl MateStager {
    fn new(mate: u8, plan: ArchivePlan, qual_codec: u8, target: usize) -> Self {
        let qual = if plan.no_quality {
            QualStager::None
        } else if plan.use_fqzcomp {
            QualStager::fqz(plan.quality_ctx_block_size, mate)
        } else {
            QualStager::bsc(mate, target)
        };
        Self {
            mate,
            seq: SeqQualHeaderStager::new(StreamRole::Sequence, mate, target),
            qual,
            hdr: HeaderColStager::new(mate),
            plan, qual_codec,
            fh: Vec::new(), fs: Vec::new(), fq: Vec::new(), frc: Vec::new(),
            seq_blocks: 0, qual_blocks: 0,
        }
    }

    /// Frame + stage one record. Emits seq/quality jobs as blocks fill; collects the
    /// header id (compressed later by the combined PairedHeaders job).
    fn stage_one(&mut self, chunk_index: u32, rec: &FastqRecord, job_tx: &SyncSender<CompressJob>) -> Result<()> {
        let p = &self.plan;
        self.fh.clear(); self.fs.clear(); self.fq.clear(); self.frc.clear();
        crate::compression::frame_record(
            rec, p.quality_mode, p.stream_quality_binning, p.skip_quality_in_stream,
            p.sequence_hints, p.rc_canon, p.const_seq_len, p.const_qual_len,
            false /* build_header_stream: paired headers always columnar, never BSC */,
            &mut self.fh, &mut self.fs, &mut self.fq, &mut self.frc,
        )?;
        if let Some(b) = self.seq.stage_record_bytes(chunk_index, &self.fs) {
            dispatch(job_tx, CompressJob::Bsc(b))?; self.seq_blocks += 1;
        }
        self.hdr.stage(&rec.id);
        match &mut self.qual {
            QualStager::Fqz { .. } => {
                let q = rec.quality.as_deref().ok_or_else(||
                    anyhow::anyhow!("fqzcomp quality encode: every record must carry quality"))?;
                if let Some(b) = self.qual.stage_fqz(chunk_index, q) {
                    dispatch(job_tx, CompressJob::Fqz(b))?; self.qual_blocks += 1;
                }
            }
            QualStager::Bsc(_) => {
                if let Some(b) = self.qual.bsc_stager().and_then(|s| s.stage_record_bytes(chunk_index, &self.fq)) {
                    dispatch(job_tx, CompressJob::Bsc(b))?; self.qual_blocks += 1;
                }
            }
            QualStager::None => {}
        }
        Ok(())
    }

    /// Flush seq + quality at chunk close, returning this mate's seq+qual SegmentSpecs and
    /// the chunk's header ids (consumed by the combined PairedHeaders job). Resets counters.
    fn close(&mut self, chunk_index: u32, job_tx: &SyncSender<CompressJob>) -> Result<(Vec<SegmentSpec>, Vec<Vec<u8>>)> {
        if let Some(b) = self.seq.flush(chunk_index) { dispatch(job_tx, CompressJob::Bsc(b))?; self.seq_blocks += 1; }
        match &mut self.qual {
            QualStager::Fqz { .. } => if let Some(b) = self.qual.flush_fqz(chunk_index) { dispatch(job_tx, CompressJob::Fqz(b))?; self.qual_blocks += 1; },
            QualStager::Bsc(_) => if let Some(b) = self.qual.bsc_stager().and_then(|s| s.flush(chunk_index)) { dispatch(job_tx, CompressJob::Bsc(b))?; self.qual_blocks += 1; },
            QualStager::None => {}
        }
        let segs = vec![
            // header spec is pushed by the orchestrator (HeaderSegment, filled by PairedHeaders);
            SegmentSpec { mate: self.mate, role: StreamRole::Sequence, codec: CODEC_BSC, source: SegSource::Blocks(self.seq_blocks) },
            SegmentSpec { mate: self.mate, role: StreamRole::Qual, codec: self.qual_codec, source: SegSource::Blocks(self.qual_blocks) },
        ];
        let ids = self.hdr.take_chunk_job(chunk_index).map(|j| j.headers).unwrap_or_default();
        self.seq.reset_for_next_chunk();
        self.qual.reset_for_next_chunk();
        self.seq_blocks = 0; self.qual_blocks = 0;
        Ok((segs, ids))
    }
}

fn dispatch(job_tx: &SyncSender<CompressJob>, job: CompressJob) -> Result<()> {
    job_tx.send(job).map_err(|_| anyhow::anyhow!("compress workers exited"))
}

/// Read up to `n` records from `reader` into an owned Vec (the planning-prelude's
/// first-chunk materialization — the ONLY whole-chunk hold; spec §3). Generic over the
/// inner reader so both whole-file (`FileReader`) and NUMA byte-range (`Take<File>`)
/// callers share it.
pub(crate) fn read_up_to<R: std::io::BufRead>(reader: &mut crate::io::FastqReader<R>, n: usize) -> Result<Vec<FastqRecord>> {
    let mut v = Vec::with_capacity(n.min(1 << 16));
    while v.len() < n {
        match reader.next()? { Some(r) => v.push(r), None => break }
    }
    Ok(v)
}

/// Per-chunk accumulators threaded through `stage_pair` → `close_paired_chunk`.
#[derive(Default)]
struct ChunkStats {
    total_bases: u64,
    original_size: u64,
    /// Per-mate decoded FASTQ output bytes (feed the ChunkDecodedSizes(2) global).
    decoded: [u64; 2],
}
impl ChunkStats {
    fn reset(&mut self) { *self = ChunkStats::default(); }
}

/// Frame + stage one pair into the two mate stagers and accumulate per-chunk stats.
fn stage_pair(
    m1: &mut MateStager, m2: &mut MateStager, chunk_index: u32,
    a: &FastqRecord, b: &FastqRecord, stats: &mut ChunkStats, no_quality: bool,
    job_tx: &SyncSender<CompressJob>,
) -> Result<()> {
    use crate::compression::decompress_impl::fastq_output_len;
    m1.stage_one(chunk_index, a, job_tx)?;
    m2.stage_one(chunk_index, b, job_tx)?;
    let qa = a.quality.as_ref().map_or(0, |q| q.len());
    let qb = b.quality.as_ref().map_or(0, |q| q.len());
    stats.total_bases += (a.sequence.len() + b.sequence.len()) as u64;
    stats.original_size += (a.id.len() + a.sequence.len() + qa + 3) as u64;
    stats.original_size += (b.id.len() + b.sequence.len() + qb + 3) as u64;
    // Per-mate decoded bytes: quality length is 0 under no_quality (the decoder emits an
    // empty quality line), matching single-end's `decoded_qual_len`. Paired is never FASTA.
    let (da, db) = if no_quality { (0, 0) } else { (qa, qb) };
    stats.decoded[0] += fastq_output_len(false, a.id.len(), a.sequence.len(), da);
    stats.decoded[1] += fastq_output_len(false, b.id.len(), b.sequence.len(), db);
    Ok(())
}

/// Flush both mates, emit the combined PairedHeaders job, and send the 6-segment manifest.
#[allow(clippy::too_many_arguments)]
fn close_paired_chunk(
    m1: &mut MateStager, m2: &mut MateStager, chunk_index: u32, pairs: u64,
    stats: &ChunkStats, target: usize,
    job_tx: &SyncSender<CompressJob>, manifest_tx: &SyncSender<WorkerMsg>,
) -> Result<()> {
    let (segs1, r1_ids) = m1.close(chunk_index, job_tx)?;
    let (segs2, r2_ids) = m2.close(chunk_index, job_tx)?;
    dispatch(job_tx, CompressJob::PairedHeaders(PairedHeaderJob {
        chunk_index, r1_ids, r2_ids, bsc_block_mb: target / (1024 * 1024),
    }))?;
    // 6-segment manifest in fixed order: R1 H/S/Q, R2 H/S/Q.
    let mut segments = Vec::with_capacity(6);
    segments.push(SegmentSpec { mate: 1, role: StreamRole::Headers, codec: CODEC_COLUMNAR, source: SegSource::HeaderSegment });
    segments.extend(segs1);
    segments.push(SegmentSpec { mate: 2, role: StreamRole::Headers, codec: CODEC_COLUMNAR, source: SegSource::HeaderSegment });
    segments.extend(segs2);
    let manifest = ChunkManifest { chunk_index, segments, decoded_output_bytes: stats.decoded,
        records: pairs, total_bases: stats.total_bases, original_size: stats.original_size };
    manifest_tx.send(WorkerMsg::Manifest(manifest)).map_err(|_| anyhow::anyhow!("writer exited"))?;
    Ok(())
}

/// Run the paired producer on THIS thread: stage the prelude's first chunk, then STREAM
/// record-by-record from both readers (NO whole-chunk buffering — spec §2 "no AoS"), closing
/// a chunk every `chunk_records` pairs. `first1`/`first2` are the already-read first chunk.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_paired_producer<R: std::io::BufRead>(
    plan1: ArchivePlan, plan2: ArchivePlan,
    qcodec1: u8, qcodec2: u8,
    first1: Vec<FastqRecord>, first2: Vec<FastqRecord>,
    reader1: &mut crate::io::FastqReader<R>,
    reader2: &mut crate::io::FastqReader<R>,
    chunk_records: usize, no_quality: bool,
    job_tx: &SyncSender<CompressJob>, manifest_tx: &SyncSender<WorkerMsg>,
    abort: &std::sync::atomic::AtomicBool,
    target: usize,
) -> Result<()> {
    use std::sync::atomic::Ordering;
    let mut m1 = MateStager::new(1, plan1, qcodec1, target);
    let mut m2 = MateStager::new(2, plan2, qcodec2, target);
    let mut chunk_index: u32 = 0;
    let mut in_chunk: u64 = 0;
    let mut stats = ChunkStats::default();

    // Phase A: the prelude's first chunk (already materialized for quality resolution).
    if first1.len() != first2.len() { bail!("unequal mate counts"); }
    for (a, b) in first1.iter().zip(first2.iter()) {
        if abort.load(Ordering::Relaxed) { return Ok(()); }
        stage_pair(&mut m1, &mut m2, chunk_index, a, b, &mut stats, no_quality, job_tx)?;
        in_chunk += 1;
        if in_chunk as usize == chunk_records {
            close_paired_chunk(&mut m1, &mut m2, chunk_index, in_chunk, &stats, target, job_tx, manifest_tx)?;
            chunk_index += 1; in_chunk = 0; stats.reset();
        }
    }
    drop(first1); drop(first2);

    // Phase B: stream record-by-record from both readers — bounded RAM, no chunk Vecs.
    loop {
        if abort.load(Ordering::Relaxed) { return Ok(()); }
        let a = reader1.next()?;
        let b = reader2.next()?;
        match (a, b) {
            (Some(x), Some(y)) => {
                stage_pair(&mut m1, &mut m2, chunk_index, &x, &y, &mut stats, no_quality, job_tx)?;
                in_chunk += 1;
                if in_chunk as usize == chunk_records {
                    close_paired_chunk(&mut m1, &mut m2, chunk_index, in_chunk, &stats, target, job_tx, manifest_tx)?;
                    chunk_index += 1; in_chunk = 0; stats.reset();
                }
            }
            (None, None) => break,
            _ => bail!("unequal mate counts"),
        }
    }
    // Final partial chunk.
    if in_chunk > 0 {
        close_paired_chunk(&mut m1, &mut m2, chunk_index, in_chunk, &stats, target, job_tx, manifest_tx)?;
    }
    Ok(())
}

pub(crate) fn compress_paired_streaming(a: &CompressConfig) -> Result<()> {
    super::reject_unsupported(a)?;
    let chunk_records = a.advanced.chunk_records.max(1);

    // Two RECORD-streaming readers (NOT chunk-buffering). Paired inputs are always files.
    let mut reader1 = crate::io::FastqReader::from_path(&a.input[0], false)?;
    let mut reader2 = crate::io::FastqReader::from_path(&a.input[1], false)?;

    // First-chunk planning prelude: materialize ONLY the first chunk from each mate (the one
    // whole-chunk hold; spec §3) to resolve per-mate quality (≥100K reads → fqz). Steady state
    // then streams record-by-record, so peak RSS is constant in read count.
    let first1 = read_up_to(&mut reader1, chunk_records)?;
    let first2 = read_up_to(&mut reader2, chunk_records)?;
    match (first1.is_empty(), first2.is_empty()) {
        (true, true) => bail!("empty paired input"),
        (a, b) if a != b => bail!("unequal mate counts"),
        _ => {}
    }
    let (q1, qcodec1) = resolve_mate_quality(a, first1.len());
    let (q2, qcodec2) = resolve_mate_quality(a, first2.len());
    let plan1 = crate::compression::stream_compress::paired_mate_plan(a, q1);
    let plan2 = crate::compression::stream_compress::paired_mate_plan(a, q2);

    compress_paired_streaming_from_readers(a, reader1, reader2, first1, first2, plan1, plan2, qcodec1, qcodec2)
}

/// Shared paired streaming core (generic over the inner reader, so whole-file `FileReader`
/// and NUMA byte-range `Take<File>` callers share it — the paired analog of single-end's
/// `run_streaming_compress`). Given the two readers, their already-read first chunks, and
/// the resolved per-mate plans, set up the worker pool + writer and run the producer.
/// `compress_paired_streaming` (whole-file) resolves quality from the prelude; the NUMA
/// `compress_byte_range_paired` passes a plan PINNED by the driver. Const lengths are always
/// 0 for paired (variable framing), so there is no first-chunk length validation to do.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compress_paired_streaming_from_readers<R: std::io::BufRead>(
    a: &CompressConfig,
    mut reader1: crate::io::FastqReader<R>,
    mut reader2: crate::io::FastqReader<R>,
    first1: Vec<FastqRecord>, first2: Vec<FastqRecord>,
    plan1: ArchivePlan, plan2: ArchivePlan,
    qcodec1: u8, qcodec2: u8,
) -> Result<()> {
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::mpsc::sync_channel;
    use std::sync::{Arc, Mutex};

    let chunk_records = a.advanced.chunk_records.max(1);
    let target = plan1.bsc_block_mb * 1024 * 1024;

    // Output: write DIRECTLY to a.output (no inner atomic temp) — the public compress()
    // dispatch already wraps non-stdout in an O_EXCL temp + rename (mod.rs:1012–1061), and
    // for stdout it routes here with a.output = the stdout sentinel. Matches legacy
    // compress_paired_v5_inner (paired/mod.rs:264–273). The writer thread needs `Send`.
    let mut output_file: Box<dyn Write + Send> = if crate::cli::is_stdio_path(&a.output) {
        Box::new(std::io::BufWriter::new(std::io::stdout()))
    } else {
        let file = std::fs::OpenOptions::new().create(true).write(true).truncate(true).open(&a.output)?;
        Box::new(std::io::BufWriter::new(file))
    };
    let header_size = crate::compression::write_v5_archive_header(
        &mut output_file, 0 /* Raw */, crate::compression::columnar::QualityBinning::None,
        crate::cli::QualityCompressor::Bsc, crate::cli::HeaderCompressor::Columnar,
        false, 0, 0, 1 /* archive_type = paired */)?;

    let n_workers = std::env::var("QZ_COMPRESS_WORKERS").ok().and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n >= 1)
        .unwrap_or_else(|| compress_worker_count(rayon::current_num_threads().max(1), target, false));
    let queue_depth = n_workers;
    let (job_tx, job_rx) = sync_channel::<CompressJob>(queue_depth);
    let (res_tx, res_rx) = sync_channel::<WorkerMsg>(queue_depth.max(1) * 2);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let abort = Arc::new(AtomicBool::new(false));
    // The worker pool's compress_job only reads block-size/static (same for both mates) +
    // each job's own fields — so ONE shared plan suffices. plan1 is moved into the producer,
    // so clone it for the workers (ArchivePlan: Clone, derived in Task 4 Step 2).
    let plan_for_workers = plan1.clone();

    std::thread::scope(|scope| -> Result<()> {
        spawn_compress_workers(scope, n_workers, &job_rx, &res_tx, &plan_for_workers, &abort);
        drop(job_rx);
        let manifest_tx = res_tx.clone();
        drop(res_tx);
        let writer = scope.spawn(move || -> Result<_> {
            write_ordered(output_file, CompressMode::Paired, header_size, res_rx)
        });
        // Producer runs on THIS (main) thread, reading both files in lockstep. The mutable
        // reader borrows never cross a thread boundary, so they need not be Send/Sync.
        let prod = run_paired_producer(plan1, plan2, qcodec1, qcodec2, first1, first2,
            &mut reader1, &mut reader2, chunk_records, a.no_quality, &job_tx, &manifest_tx, &abort, target);
        if let Err(e) = &prod {
            abort.store(true, Ordering::Relaxed);
            let _ = manifest_tx.send(WorkerMsg::Err(e.to_string()));
        }
        drop(job_tx); drop(manifest_tx);
        let _done = writer.join().map_err(|_| anyhow::anyhow!("writer thread panicked"))??;
        Ok(())
    })?;
    Ok(())
}
