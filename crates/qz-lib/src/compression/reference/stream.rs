//! Block-granular streaming reference compress (spec Increment 3). Reference's
//! alignment front end (mapping → diff → columns → consensus) stays OUTSIDE the
//! shared block-streaming core; only the *concurrency primitive* is shared: the
//! per-chunk `encode_chunk` work is driven by the generic `spawn_worker_pool`
//! (the same pool single/paired use), and a thin chunk-ordered writer reorders
//! the encoded chunks + streams the whole-archive globals in its tail.
//!
//! The dividing line (vs single/paired): reference's per-read data is *chunk-bound*
//! by alignment, so `encode_chunk` (unchanged) returns one self-contained
//! `EncodedChunk { bytes, entries }` per chunk rather than block-granular jobs. The
//! producer folds coverage SERIALLY in chunk order (IntervalMap byte-identity
//! depends on it), finalizes the `IntervalMap` after the last chunk, and hands it to
//! the writer to STREAM the globals (PackedBacking is never materialized). Output is
//! byte-identical to the legacy bespoke producer/token-semaphore consumer.

use anyhow::{Result, anyhow, bail};
use std::collections::BTreeMap;
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use crate::cli::CompressConfig;
use crate::compression::chunk_directory::{
    self, ChunkDirEntry, ChunkDirectory, StreamRole,
};
use crate::compression::compress_impl::spawn_worker_pool;
use crate::io::FastqRecord;

use super::backing::IntervalMap;
use super::encode::{EncodedChunk, encode_chunk};
use super::mapping::Mapper;
use super::{
    coverage_to_intervalmap, encode, map_and_diff, mark_cov, reference_prelude,
    resolve_ref_chunk_records_cfg,
};

/// Reference's OWN pool job type — never enters the shared block-streaming core.
/// One job = one chunk's per-mate column/header/quality encode (`encode_chunk`).
enum RefJob {
    Chunk {
        chunk_index: u32,
        records: Vec<(FastqRecord, FastqRecord)>,
        placed: Vec<(super::Placed, super::Placed)>,
        qctx_block_size: usize,
    },
}

/// Reference's OWN pool/producer result type. `Chunk` is a worker result; `Globals`
/// is sent by the PRODUCER (not the pool — globals are inherently serial + streamed,
/// not a job); `Err` carries a worker/panic message (built by the pool's `on_error`).
pub(super) enum RefResult {
    Chunk { chunk_index: u32, ec: EncodedChunk },
    Globals { imap: IntervalMap },
    Err(String),
}

/// Which reference-producer phase a timing sample belongs to.
pub(super) enum RefPhase {
    /// `map_and_diff` — strobealign mapping + Hamming diff. The measured gate (~60–72% of wall).
    Map,
    /// The serial in-order coverage fold (`mark_cov`). Measured <1% of wall.
    Cov,
    /// Blocked sending a chunk to the worker pool (pool starvation; ~0 when mapping is the gate).
    Send,
    /// `coverage_to_intervalmap` finalize (one O(genome-length) pass; serial producer tail).
    Finalize,
}

/// Per-phase wall-time accumulator for the reference producer, printed to stderr when
/// `QZ_REF_TIMING` is set. Pure diagnostic: timing is observational (the encoded bytes are
/// unchanged) and `Instant::now()` runs only per-chunk (chunk = `chunk_records` reads), so the
/// cost is negligible — and `tick()` returns `None` (no clock read) when the env var is unset.
/// Shared by the paired (`stream.rs`) and single-end (`stream_single.rs`) producers so the
/// phase split is identical across both. Confirmed `2026-06-25`: the gate is `Map`, not `Cov`.
pub(super) struct RefProdTiming {
    enabled: bool,
    map: Duration,
    cov: Duration,
    send: Duration,
    finalize: Duration,
}

impl RefProdTiming {
    pub(super) fn new() -> Self {
        Self {
            enabled: std::env::var_os("QZ_REF_TIMING").is_some(),
            map: Duration::ZERO,
            cov: Duration::ZERO,
            send: Duration::ZERO,
            finalize: Duration::ZERO,
        }
    }

    /// Start a phase timer. `None` when disabled → zero clock reads on the shipped path.
    pub(super) fn tick(&self) -> Option<Instant> {
        self.enabled.then(Instant::now)
    }

    /// Accumulate the time elapsed since `t` into `phase` (no-op when `t` is `None`).
    pub(super) fn record(&mut self, phase: RefPhase, t: Option<Instant>) {
        if let Some(t) = t {
            let d = t.elapsed();
            match phase {
                RefPhase::Map => self.map += d,
                RefPhase::Cov => self.cov += d,
                RefPhase::Send => self.send += d,
                RefPhase::Finalize => self.finalize += d,
            }
        }
    }

    /// Emit the phase breakdown to stderr (no-op when disabled).
    pub(super) fn report(&self, label: &str, chunks: u32) {
        if !self.enabled {
            return;
        }
        let s = |d: Duration| d.as_secs_f64();
        eprintln!(
            "[QZ_REF_TIMING] {label}: chunks={chunks} map_and_diff={:.3}s cov_fold={:.3}s \
             send_blocked={:.3}s finalize={:.3}s",
            s(self.map),
            s(self.cov),
            s(self.send),
            s(self.finalize),
        );
    }
}

/// Debug-only fault injection for the abort/shutdown contract regression tests.
/// Reads `QZ_REF_FAIL` once (per call site); shipped builds compile this out
/// entirely (zero overhead, no env read) — mirrors the `QZ_NUMA_FAKE_UNDERSIZE`
/// pattern. Recognized values:
///   - `worker:N`   → `run_ref_job` returns Err on chunk index N
///   - `producer:N` → `run_reference_producer` returns Err after dispatching N chunks
///
/// The two trigger strings are distinct so the two tests can run in parallel
/// without a serial guard (each ignores the other's value).
#[cfg(debug_assertions)]
fn ref_fail_target(kind: &str) -> Option<u32> {
    let v = std::env::var("QZ_REF_FAIL").ok()?;
    let n = v.strip_prefix(kind)?.strip_prefix(':')?;
    n.parse::<u32>().ok()
}

/// The pool work fn: encode ONE chunk via the unchanged `encode_chunk` (with its
/// `chunk_index`). Returns exactly one result. Globals are NOT a pool job.
fn run_ref_job(job: RefJob) -> Result<Vec<RefResult>> {
    let RefJob::Chunk { chunk_index, records, placed, qctx_block_size } = job;
    // Debug-only failpoint: bail with a recognizable message on the target chunk so
    // the abort-contract test can assert the WORKER error is surfaced (and that the
    // call still returns + leaves no partial archive).
    #[cfg(debug_assertions)]
    if ref_fail_target("worker") == Some(chunk_index) {
        bail!("injected worker failure at chunk {chunk_index}");
    }
    let ec = encode_chunk(chunk_index, &records, &placed, qctx_block_size)?;
    Ok(vec![RefResult::Chunk { chunk_index, ec }])
}

/// Run the reference producer on THIS (main) thread: pull pairs prefix-then-readers,
/// chunk at `chunk_records` (mirrors the legacy producer's `prefix`-then-readers
/// drain), fold coverage SERIALLY in chunk order, and dispatch one `RefJob::Chunk`
/// per chunk. After the last chunk, finalize the `IntervalMap` and send it to the
/// writer to STREAM as the whole-archive globals.
#[allow(clippy::too_many_arguments)]
fn run_reference_producer(
    mapper: &Mapper,
    read_len: usize,
    qctx_block_size: usize,
    prefix: Vec<(FastqRecord, FastqRecord)>,
    reader1: &mut crate::io::FastqReader<crate::io::fastq::FileReader>,
    reader2: &mut crate::io::FastqReader<crate::io::fastq::FileReader>,
    chunk_records: usize,
    job_tx: &SyncSender<RefJob>,
    globals_tx: &SyncSender<RefResult>,
    abort: &AtomicBool,
) -> Result<()> {
    // Coverage bitmap over the normalized reference (1 bit/base; genome-bounded).
    let contig_lens: Vec<u64> = mapper
        .references()
        .iter()
        .map(|r| r.sequence.len() as u64)
        .collect();
    let mut coverage = super::coverage::CoverageMap::new(&contig_lens);

    let mut prefix = prefix.into_iter();
    let mut chunk_index: u32 = 0u32;
    let mut buf: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(chunk_records);

    // Pull one pair: drain the WHOLE sampled prefix first (so no sampled pair is
    // dropped or rechunked), then stream both readers in lockstep.
    let mut next_pair = |reader1: &mut crate::io::FastqReader<crate::io::fastq::FileReader>,
                         reader2: &mut crate::io::FastqReader<crate::io::fastq::FileReader>|
     -> Result<Option<(FastqRecord, FastqRecord)>> {
        if let Some(p) = prefix.next() {
            return Ok(Some(p));
        }
        match (reader1.next()?, reader2.next()?) {
            (Some(a), Some(b)) => Ok(Some((a, b))),
            (None, None) => Ok(None),
            _ => bail!("reference mode: R1/R2 have unequal read counts"),
        }
    };

    let mut timing = RefProdTiming::new();

    // Map+diff one chunk, fold its coverage in order, and dispatch its encode job.
    let emit = |buf: &mut Vec<(FastqRecord, FastqRecord)>,
                ci: &mut u32,
                cov: &mut super::coverage::CoverageMap,
                tm: &mut RefProdTiming|
     -> Result<()> {
        let records = std::mem::take(buf);
        let t = tm.tick();
        let placed = map_and_diff(mapper, &records)?;
        tm.record(RefPhase::Map, t);
        // Serial in-order coverage fold (IntervalMap byte-identity depends on it).
        let t = tm.tick();
        for (p1, p2) in &placed {
            mark_cov(cov, p1);
            mark_cov(cov, p2);
        }
        tm.record(RefPhase::Cov, t);
        let t = tm.tick();
        job_tx
            .send(RefJob::Chunk { chunk_index: *ci, records, placed, qctx_block_size })
            .map_err(|_| anyhow!("reference workers exited"))?;
        tm.record(RefPhase::Send, t);
        *ci += 1;
        Ok(())
    };

    while let Some(p) = next_pair(reader1, reader2)? {
        if abort.load(Ordering::Relaxed) {
            return Ok(());
        }
        buf.push(p);
        if buf.len() == chunk_records {
            emit(&mut buf, &mut chunk_index, &mut coverage, &mut timing)?;
            // Debug-only failpoint: after dispatching the target number of chunks,
            // bail so the abort-contract test can assert the PRODUCER error is
            // surfaced (writer drains, sees `abort`, returns Ok; `w.and(prod)` wins).
            #[cfg(debug_assertions)]
            if ref_fail_target("producer") == Some(chunk_index) {
                bail!("injected producer failure after {chunk_index} chunks");
            }
        }
    }
    // Final partial chunk.
    if !buf.is_empty() {
        emit(&mut buf, &mut chunk_index, &mut coverage, &mut timing)?;
    }

    // Globals: finalize coverage (one O(genome-length) pass over the bitmap — a
    // serial producer tail, modest at chr20 scale but ~genome-sized at WGS) and
    // hand the IntervalMap to the writer to STREAM (not a pool job). `refs` are
    // borrowed by the writer from `&mapper`.
    let t = timing.tick();
    let imap = coverage_to_intervalmap(coverage, read_len)?;
    timing.record(RefPhase::Finalize, t);
    timing.report("paired", chunk_index);
    globals_tx
        .send(RefResult::Globals { imap })
        .map_err(|_| anyhow!("reference writer exited"))?;
    Ok(())
}

/// The thin chunk-ordered writer SHARED by paired (type 2) and single-end (type 4)
/// reference: reorder encoded chunks by `chunk_index`, flush each in order (shifting
/// per-entry offsets by the running base), then in the TAIL stream the whole-archive
/// globals + write the v5 footer + locator. `validate` is the mode's directory
/// self-check (`validate_reference_directory` for paired, `…_single` for single-end);
/// it only inspects the finished directory and never emits bytes, so both modes'
/// output stays byte-identical to a hard-coded call. Never panics on abort: if the
/// producer aborted, it returns `Ok(())` WITHOUT finalizing so the producer's real
/// error wins via `w.and(prod)`.
pub(super) fn write_reference(
    mut out: Box<dyn Write + Send>,
    header_size: u64,
    mapper: &Mapper,
    res_rx: Receiver<RefResult>,
    abort: &AtomicBool,
    validate: fn(&ChunkDirectory) -> Result<()>,
    n_mates: u8,
) -> Result<()> {
    let mut cur: u64 = header_size;
    let mut next: u32 = 0;
    let mut pending: BTreeMap<u32, EncodedChunk> = BTreeMap::new();
    let mut entries: Vec<ChunkDirEntry> = Vec::new();
    let mut num_reads: u64 = 0;
    let mut num_chunks: u32 = 0;
    let mut imap: Option<IntervalMap> = None;
    // Per-chunk [R1, R2] decoded output bytes, collected in chunk order as chunks flush
    // → one `ChunkDecodedSizes(n_mates=2)` global in the tail (for NUMA direct-write).
    let mut decoded_sizes: Vec<[u64; 2]> = Vec::new();

    for msg in res_rx {
        match msg {
            // A worker error is AUTHORITATIVE: return it (the producer may
            // secondarily see "workers exited"; `w.and(prod)` discards that).
            RefResult::Err(e) => return Err(anyhow!(e)),
            RefResult::Globals { imap: m } => imap = Some(m),
            RefResult::Chunk { chunk_index, ec } => {
                num_chunks = num_chunks.max(chunk_index + 1);
                pending.insert(chunk_index, ec);
                while let Some(ec) = pending.remove(&next) {
                    let base = cur;
                    decoded_sizes.push(ec.decoded_output_bytes);
                    out.write_all(&ec.bytes)?;
                    for mut e in ec.entries {
                        // num_reads = Σ over chunks of the mate-1 Headers entry's
                        // record_count (= the chunk's pair count) — exactly the
                        // legacy `num_pairs` (one ++ per pair, summed per chunk).
                        accumulate_num_reads(&mut num_reads, &e);
                        e.offset += base;
                        entries.push(e);
                    }
                    cur += ec.bytes.len() as u64;
                    next += 1;
                }
            }
        }
    }

    // Producer aborted (set `abort` + dropped its senders without sending Globals):
    // do NOT finalize a partial archive — its real error surfaces via `w.and(prod)`.
    if abort.load(Ordering::Relaxed) {
        return Ok(());
    }

    let imap = imap.ok_or_else(|| anyhow!("reference: missing globals result"))?;
    let refs = mapper.references();
    // Stream the 4 whole-archive globals block-by-block (PackedBacking → NBitmap →
    // IntervalMap → ReferenceMeta), NOT materialized.
    encode::stream_globals(&mut out, &imap, refs, &mut cur, &mut entries)?;

    // One ChunkDecodedSizes(n_mates=2) global (row-major `[r1_c0, r2_c0, r1_c1, …]`) so
    // reference NUMA decode can direct-write into pre-sized R1/R2 outputs (mirrors the
    // single/paired encoder). The reference directory validator treats it as a
    // cross-cutting global (skipped, not a RefRole).
    encode::emit_decoded_sizes_global(&mut out, n_mates, num_chunks, &decoded_sizes, &mut cur, &mut entries)?;

    // v5 directory footer + trailing 20-byte locator (seek-free).
    let dir = ChunkDirectory {
        num_reads,
        num_chunks,
        entries,
    };
    // Encode-time self-check (mode-specific, passed in by the caller): catch an
    // encoder bug now rather than emit a structurally-invalid archive that only
    // fails on decode. The validator never emits bytes, so output is unaffected.
    validate(&dir)?;
    chunk_directory::write_footer_and_locator(&mut out, &dir)?;
    out.flush()?;
    Ok(())
}

/// Accumulate `num_reads` from a chunk's directory entry: each chunk contributes its
/// mate-1 Headers entry's `record_count` (= the chunk's pair count). Summed across
/// chunks this equals the legacy `num_pairs` (which incremented once per pair).
fn accumulate_num_reads(num_reads: &mut u64, e: &ChunkDirEntry) {
    if e.role == StreamRole::Headers && e.mate == 1 {
        *num_reads += e.record_count;
    }
}

/// Streaming reference compress entry point. Mirrors `compress_paired_streaming`'s
/// scope/shutdown shape: prelude → spawn the shared worker pool → `drop(job_rx)` →
/// writer thread (owns `res_rx`) → run the producer on the main thread (holding a
/// `globals_tx = res_tx.clone()`) → drop the senders → join. Writer error is
/// authoritative (`w.and(prod)`).
/// `ranges` is `Some([r1, r2])` for a NUMA byte-range compress worker (each mate
/// bounded to its record-aligned shard) or `None` for a whole-file compress; it only
/// affects the prelude's readers — everything downstream (header, pool, writer,
/// producer) is identical, and the fixed reference front header means every shard's
/// header is byte-identical (so the reference merge's front-header gate passes).
pub(crate) fn compress_reference_streaming(
    a: &CompressConfig,
    ranges: Option<&[crate::compression::ByteRange]>,
) -> Result<()> {
    // Pass-0 planning prelude: profile read_len + insert μ/σ, build the mapper, and
    // buffer the WHOLE sampled prefix (drained prefix-then-readers by the producer).
    let (mapper, read_len, mut r1, mut r2, prefix, qctx_block_size) = reference_prelude(a, ranges)?;
    let chunk_records = resolve_ref_chunk_records_cfg(a);

    // Reference rejects stdout (gated in reject_unsupported); the public compress()
    // dispatch hands us the sibling O_EXCL `.tmp` path (with a drop guard removing it
    // on any error) — write directly to it. `File::create` truncates the reserved
    // temp; the writer thread needs `Send`.
    let file = std::fs::File::create(&a.output)
        .map_err(|e| anyhow!("create reference archive {}: {e}", a.output.display()))?;
    let mut out: Box<dyn Write + Send> = Box::new(std::io::BufWriter::new(file));
    // v5 front header (reference = archive_type 2). The codec/quality fields here are
    // metadata only — the directory's per-entry codec is authoritative.
    let header_size = crate::compression::write_v5_archive_header(
        &mut out,
        0, // encoding_type (Raw)
        crate::compression::columnar::QualityBinning::None,
        crate::cli::QualityCompressor::Bsc,
        crate::cli::HeaderCompressor::Columnar,
        false, // fasta
        0,     // const_seq_len
        0,     // const_qual_len
        2,     // archive_type = reference
    )? as u64;

    // Worker count = reference's window (chunk-bound RAM): each chunk job is heavy
    // (fqzcomp quality + BSC columns), so it must NOT run ~46-wide like the block
    // default. `QZ_COMPRESS_WORKERS` overrides for benchmarking.
    let n_workers = std::env::var("QZ_COMPRESS_WORKERS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n >= 1)
        .unwrap_or_else(|| {
            a.reference
                .as_ref()
                .map(|r| r.reference_window)
                .unwrap_or(4)
                .max(1)
        });

    let (job_tx, job_rx) = sync_channel::<RefJob>(n_workers);
    let (res_tx, res_rx) = sync_channel::<RefResult>(n_workers.max(1) * 2);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let abort = Arc::new(AtomicBool::new(false));

    std::thread::scope(|scope| -> Result<()> {
        spawn_worker_pool(
            scope,
            n_workers,
            &job_rx,
            &res_tx,
            &abort,
            run_ref_job,
            RefResult::Err,
        );
        drop(job_rx);

        // The producer sends its final globals message through this clone.
        let globals_tx = res_tx.clone();
        let writer = {
            let abort = abort.clone();
            let mapper = &mapper;
            scope.spawn(move || {
                write_reference(
                    out,
                    header_size,
                    mapper,
                    res_rx,
                    &abort,
                    super::validate_reference_directory,
                    2, // paired reference: R1 + R2
                )
            })
        };
        // Only `globals_tx` (producer) + the workers hold senders now; the writer's
        // `res_rx` ends once they all drop.
        drop(res_tx);

        // Producer runs on THIS (main) thread, reading both files in lockstep — the
        // mutable reader borrows never cross a thread boundary.
        let prod = run_reference_producer(
            &mapper,
            read_len,
            qctx_block_size,
            prefix,
            &mut r1,
            &mut r2,
            chunk_records,
            &job_tx,
            &globals_tx,
            &abort,
        );
        if prod.is_err() {
            // Deliberately just SET abort (don't route the error to the writer the way
            // paired/stream.rs does): the writer's abort-check returns Ok without
            // finalizing, and `w.and(prod)` below surfaces this producer error. Simpler
            // than the paired shape and equivalent — do not "align" it by sending a
            // WorkerMsg::Err here (reference has no such message and doesn't need one).
            abort.store(true, Ordering::Relaxed);
        }
        // Close senders → the writer's `res_rx` ends and it finalizes (or, on abort,
        // returns Ok without finalizing).
        drop(job_tx);
        drop(globals_tx);

        let w = writer
            .join()
            .map_err(|_| anyhow!("reference writer panicked"))?;
        // Writer error authoritative; the producer's error only surfaces when the
        // writer closed clean (abort path).
        w.and(prod)
    })
}
