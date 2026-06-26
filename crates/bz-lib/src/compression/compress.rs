use crate::cli::{AdvancedOptions, CompressConfig};
use crate::compression::archive::{
    ArchiveHeader, CHUNK_FLAG_FQZ_QUALITY, ChunkHeader, FLAG_CONSENSUS_DELTA, FLAG_QUALITY_LOSSY,
};
use crate::compression::streams::{self, BamStreams, ChunkQualityData, NUM_STREAMS};
use crate::io::bam::{RawBamReader, RawBamRecord, SortOrder};
use anyhow::Result;
use flate2::Crc;
use qz_lib::compression::bsc;
use qz_lib::compression::fqzcomp;
use std::collections::VecDeque;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use std::time::Instant;
use tracing::info;

/// Working memory one in-flight chunk holds during compression, per record:
/// raw columnar streams + compressed outputs + BSC/libsais workspace. Measured
/// ~1.3-1.4 KB/record on HG002 chr20; 1400 is a slightly conservative estimate
/// used only to bound the auto-selected pipeline window.
const PER_RECORD_CHUNK_BYTES: u64 = 1400;
/// Fixed compress overhead outside the per-chunk pipeline (BGZF reader buffers,
/// header payload, allocator pool retention).
const COMPRESS_BASE_OVERHEAD_BYTES: u64 = 4 * 1024 * 1024 * 1024;
/// Safety ceiling on the auto-selected window. Returns diminish past this on a
/// 72-core box (window 32 still only reached ~31 busy cores) while RSS keeps
/// climbing, so cap regardless of how much RAM is available.
const AUTO_WINDOW_CAP: usize = 16;
/// Fallback window when available RAM cannot be read (non-Linux, no /proc).
const AUTO_WINDOW_FALLBACK: usize = 4;

/// Available system memory in bytes (Linux `/proc/meminfo` `MemAvailable`).
pub(super) fn available_memory_bytes() -> Option<u64> {
    let contents = std::fs::read_to_string("/proc/meminfo").ok()?;
    for line in contents.lines() {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            return rest
                .split_whitespace()
                .next()?
                .parse::<u64>()
                .ok()
                .map(|kb| kb * 1024);
        }
    }
    None
}

/// Resolve the effective compress pipeline window. A non-zero `requested` value
/// is used as-is (floored at 1). `requested == 0` means auto: pick the largest
/// window whose estimated peak RSS stays under ~70% of available RAM, bounded by
/// the core count and [`AUTO_WINDOW_CAP`]. The window only caps concurrency, so
/// over-provisioning is harmless on small files (fewer chunks simply never fill
/// it) — and it never changes the archive bytes.
fn resolve_pipeline_window(requested: usize, chunk_size: usize) -> usize {
    if requested != 0 {
        return requested.max(1);
    }
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8);
    let per_chunk = (chunk_size as u64)
        .saturating_mul(PER_RECORD_CHUNK_BYTES)
        .max(1);
    let mem_window = match available_memory_bytes() {
        Some(avail) => {
            let budget = avail / 10 * 7; // ~70% of available RAM
            let usable = budget.saturating_sub(COMPRESS_BASE_OVERHEAD_BYTES);
            (usable / per_chunk).max(1) as usize
        }
        None => AUTO_WINDOW_FALLBACK,
    };
    mem_window.min(cores).clamp(1, AUTO_WINDOW_CAP)
}

/// A compression level preset: chunk size + pipeline window.
///
/// NOTE: bz's ratio is **~flat across levels** — every stream is split into fixed
/// 25 MB BSC blocks regardless of chunk size, and the per-position consensus
/// saturates at a few hundred K reads, so chunk size is effectively a
/// **memory/speed** knob, not a ratio knob (unlike qz's `ULTRA_LEVELS`, where
/// bigger BSC blocks genuinely help). Measured on HG002 chr20 (lossless): L1
/// 5.2 GB / 1.998×, L2 9.9 GB / 2.000× (fastest), L3 18 GB / 2.001×. Lower levels
/// are nearly free on ratio; L3 only costs memory here but stays available for
/// data where bigger chunks might help (e.g. low coverage).
struct BzLevel {
    level: u8,
    chunk_size: usize,
    window: usize,
}

const BZ_LEVELS: [BzLevel; 3] = [
    // L1: lowest memory — small chunks. The downshift target under RAM pressure.
    BzLevel {
        level: 1,
        chunk_size: 500_000,
        window: 6,
    },
    // L2: balanced / fastest — the default. Window widened (vs the original 3) to
    // use the memory freed by the sequences_ascii + lazy-qual removals for more
    // pipeline parallelism — profiling showed compress plateaus at ~13 cores when
    // window=3 starves the pipeline.
    BzLevel {
        level: 2,
        chunk_size: 1_500_000,
        window: 6,
    },
    // L3: large chunks — most memory, ~same ratio here; explicit-only.
    BzLevel {
        level: 3,
        chunk_size: 5_000_000,
        window: 3,
    },
];

/// Estimated peak RSS for a level: chunks-in-flight × per-chunk + base overhead.
fn estimate_level_peak(l: &BzLevel) -> u64 {
    (l.chunk_size as u64)
        .saturating_mul(PER_RECORD_CHUNK_BYTES)
        .saturating_mul(l.window as u64)
        .saturating_add(COMPRESS_BASE_OVERHEAD_BYTES)
}

/// Auto-select a level. Because ratio is ~flat (see [`BzLevel`]), levels are a
/// memory/speed trade: default to the balanced/fastest level (L2) and downshift
/// to the low-memory level (L1) only when L2's estimated peak won't fit ~70% of
/// available RAM. L3 is never auto-selected (more memory, no ratio gain here).
fn auto_select_level() -> u8 {
    let budget = available_memory_bytes().map(|a| a / 10 * 7);
    match budget {
        Some(b) if estimate_level_peak(&BZ_LEVELS[1]) > b => 1, // tight RAM → L1
        _ => 2,                                                 // balanced default
    }
}

/// Resolve `opts.level` into concrete `chunk_size` + `compress_window`. `level 0`
/// auto-selects by RAM; levels 1–3 use the preset (overriding any explicit
/// chunk_size/window). An out-of-range level is treated as auto.
fn apply_level(mut opts: AdvancedOptions) -> AdvancedOptions {
    let lvl = if opts.level == 0 || opts.level > 3 {
        auto_select_level()
    } else {
        opts.level
    };
    if let Some(l) = BZ_LEVELS.iter().find(|l| l.level == lvl) {
        opts.chunk_size = l.chunk_size;
        opts.compress_window = l.window;
        opts.level = lvl; // record the resolved level
        info!(
            "Compression level {} (chunk_size={}, window={}, est. peak ~{:.1} GB)",
            lvl,
            l.chunk_size,
            l.window,
            estimate_level_peak(l) as f64 / 1e9
        );
    }
    opts
}

/// Resolve `config`'s compression level to the concrete `1`/`2`/`3` it will use
/// (level `0`/out-of-range = auto, picked by RAM via [`auto_select_level`]).
///
/// The NUMA-compress driver calls this ONCE up front to **freeze** the level
/// before sharding: `apply_level` re-runs `auto_select_level` (which reads live
/// `/proc/meminfo`) on every call, so the prescan and each re-exec'd worker —
/// separate processes, running at different times while siblings allocate GBs —
/// could otherwise resolve DIFFERENT levels (hence different `chunk_size`) and the
/// worker's chunking would diverge from the prescan's virtual-offset partition
/// (silent record loss or a hard fail). Pinning the resolved level into the config
/// makes every downstream `apply_level` an idempotent no-op. This mirrors
/// `apply_level`'s level choice with NO side effects (no chunk_size/window change,
/// no logging).
pub fn resolve_compress_level(config: &CompressConfig) -> u8 {
    let lvl = config.advanced.level;
    if lvl == 0 || lvl > 3 {
        auto_select_level()
    } else {
        lvl
    }
}

/// Compress a single stream with BSC (the only backend).
fn compress_one(data: &[u8], opts: &AdvancedOptions) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }

    // Select BSC variant based on config.
    // use_lzp only affects the adaptive path; the static path has no LZP knob.
    let mut compressed = match (opts.bsc_adaptive, opts.use_lzp) {
        (true, true) => bsc::compress_parallel_adaptive(data)?,
        (true, false) => bsc::compress_parallel_adaptive_no_lzp(data)?,
        (false, _) => bsc::compress_parallel(data)?,
    };
    compressed.shrink_to_fit();
    Ok(compressed)
}

/// Result returned by a chunk compression thread.
struct ChunkResult {
    /// 15 compressed streams.
    compressed: Vec<Vec<u8>>,
    /// Number of BAM records in this chunk.
    num_records: usize,
    /// Per-chunk archive flags (e.g. CHUNK_FLAG_FQZ_QUALITY).
    chunk_flags: u8,
    /// Raw (uncompressed) sizes per stream, for logging.
    raw_sizes: [usize; NUM_STREAMS],
    /// CRC32 over the original BAM record bytes of this chunk (fidelity check).
    content_crc32: u32,
}

/// Compress one chunk of BAM records.
///
/// Runs on the calling thread (or a spawned thread) and uses the global rayon
/// pool for stream-level parallelism.  Returns a [`ChunkResult`] containing the
/// 15 compressed streams plus metadata needed for archive assembly.
fn compress_one_chunk(
    records: Vec<RawBamRecord>,
    opts: AdvancedOptions,
    use_fqz_quality: bool,
    chunk_idx: u32,
) -> Result<ChunkResult> {
    let n = records.len();
    info!("Chunk {}: {} records (compressing)", chunk_idx, n);

    let mut result = streams::records_to_streams(&records, opts.quality_reduction.as_ref());

    // Fidelity hash. Decompress/verify recompute it over the reconstructed records
    // to prove the chunk round-trips. In lossless mode it's a plain CRC over the
    // original record bytes. In lossy quality-reduction mode (only when the chunk
    // actually used the reduction — i.e. quality was available), the quality
    // region is *intentionally* altered, so the hash must cover the REDUCED
    // quality (what decode reconstructs): splice [pre-qual] + [reduced raw qual]
    // + [post-qual] per record. The reduced ascii lives in quality_data.
    let lossy = opts.quality_reduction.is_some() && !result.quality_data.has_unavailable_quality;
    let content_crc32 = {
        let mut content_crc = Crc::new();
        if lossy {
            for (rec, qual_ascii) in records.iter().zip(&result.quality_data.qualities_ascii) {
                let data = rec.as_bytes();
                let qual = rec.qual_bytes();
                let qstart = qual.as_ptr() as usize - data.as_ptr() as usize;
                let qend = qstart + qual.len();
                content_crc.update(&data[..qstart]);
                // reduced quality, back to raw BAM Phred (ascii - 33). Use
                // wrapping_sub to match decode's ascii_qual_to_bam exactly, so
                // the content CRC can never diverge from what decode reconstructs
                // (every reduced ascii byte is ≥ 33 today, but keep them in lockstep).
                let reduced_raw: Vec<u8> =
                    qual_ascii.iter().map(|&a| a.wrapping_sub(33)).collect();
                content_crc.update(&reduced_raw);
                content_crc.update(&data[qend..]);
            }
        } else {
            for rec in &records {
                content_crc.update(rec.as_bytes());
            }
        }
        content_crc.sum()
    };

    // Per-chunk decision: use the fqzcomp quality path only if enabled AND no
    // 0xFF (unavailable) quality. The BSC fallback path instead consumes the raw
    // quality stream, which is built lazily here from the records (still alive) —
    // the 0xFF bytes it needs aren't recoverable from quality_data.qualities_ascii
    // (masked there). Normal chunks skip this entirely.
    let chunk_uses_fqz_quality = use_fqz_quality && !result.quality_data.has_unavailable_quality;
    if !(chunk_uses_fqz_quality && n > 0) {
        streams::populate_qual_stream(&mut result.streams, &records);
    }

    drop(records);

    let streams::StreamsResult {
        streams: bam_streams,
        quality_data,
    } = result;

    // Capture raw sizes before consuming bam_streams
    let raw_slices = bam_streams.as_slices();
    let mut raw_sizes = [0usize; NUM_STREAMS];
    for (i, s) in raw_slices.iter().enumerate() {
        raw_sizes[i] = s.len();
    }

    let compressed = if chunk_uses_fqz_quality && n > 0 {
        let fqz_block_size = opts.fqz_block_size;
        let (bsc_result, fqz_result) = rayon::join(
            || compress_streams(&bam_streams, &opts, true),
            || compress_quality_fqz_parallel(&quality_data, fqz_block_size),
        );
        let mut compressed = bsc_result?;
        compressed[13] = fqz_result?;
        compressed
    } else {
        compress_streams(&bam_streams, &opts, false)?
    };

    let mut chunk_flags = 0u8;
    if chunk_uses_fqz_quality {
        chunk_flags |= CHUNK_FLAG_FQZ_QUALITY;
    }

    Ok(ChunkResult {
        compressed,
        num_records: n,
        chunk_flags,
        raw_sizes,
        content_crc32,
    })
}

/// Compress a BAM file to a BZ archive.
pub fn compress(config: &CompressConfig) -> Result<()> {
    let start_time = Instant::now();
    // Validate options at the library boundary. The CLI validates parsed args, but
    // the public `bz_lib::compress` entry can be called programmatically with raw
    // options — an unchecked `fqz_block_size == 0` would reach `.step_by(0)`
    // and panic instead of returning an error.
    config.advanced.validate()?;
    // Resolve the level preset (chunk_size + window) before anything reads them.
    let opts = apply_level(config.advanced.clone());
    let opts = &opts;

    info!(
        "BZ compression config: chunk_size={}, fqz_block={}, \
         lzp={}, adaptive={}, qual_comp={}, window={}",
        opts.chunk_size,
        opts.fqz_block_size,
        opts.use_lzp,
        opts.bsc_adaptive,
        opts.quality_compressor,
        opts.compress_window,
    );

    // Use multi-threaded BGZF reader for parallel inflate. The main thread reads
    // each chunk synchronously between dispatches, so a too-small reader pool
    // leaves serial troughs in the compress pipeline (off-CPU sampling showed the
    // running-thread count dropping to ~9 during reads). min(cores, 16) keeps the
    // reader ahead of the encoders: ~5% faster compress on HG002 chr20,
    // byte-identical archive (reader count never affects output).
    let bgzf_workers = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::new(1).unwrap())
        .min(NonZeroUsize::new(16).unwrap());
    let mut reader = RawBamReader::from_path_mt(&config.input, bgzf_workers)?;

    // Verify the BAM is coordinate-sorted — required for consensus-delta encoding.
    match crate::io::bam::parse_sort_order(&reader.header_raw) {
        SortOrder::Coordinate => {}
        SortOrder::Other(so) => {
            anyhow::bail!(
                "BZ requires a coordinate-sorted BAM file, but this file declares SO:{}. \
                 Sort it first with: samtools sort -o sorted.bam {}",
                so,
                config.input.display(),
            );
        }
        SortOrder::Unknown => {
            tracing::warn!(
                "BAM file does not declare a sort order (no @HD SO: tag). \
                 BZ requires coordinate-sorted input — compression may produce \
                 incorrect results if the file is not sorted. \
                 Sort with: samtools sort -o sorted.bam {}",
                config.input.display(),
            );
        }
    }

    // Compress SAM header + ref dict together
    let mut header_payload = Vec::new();
    let header_raw_len = u32::try_from(reader.header_raw.len())
        .map_err(|_| anyhow::anyhow!("SAM header too large ({} bytes)", reader.header_raw.len()))?;
    header_payload.extend_from_slice(&header_raw_len.to_le_bytes());
    header_payload.extend_from_slice(&reader.header_raw);
    header_payload.extend_from_slice(&reader.ref_dict_bytes);
    let sam_header_compressed = bsc::compress_parallel(&header_payload)?;

    // Write to a randomized sibling temp file created with O_EXCL, then rename
    // atomically on success. The random name avoids collisions between concurrent
    // jobs (a deterministic `foo.bz.tmp` is shared by outputs `foo` and `foo.bz`),
    // and O_EXCL refuses to follow a pre-planted symlink (which would truncate an
    // unrelated writable file). `tempfile`'s `TempPath` removes the file on any
    // early return (error or panic); `into_parts` hands us the owned `File` so the
    // writer type (and the later `seek`-to-patch-header step) are unchanged.
    let (output_file, temp_path) = tempfile::Builder::new()
        .prefix(".bz-")
        .suffix(".tmp")
        .tempfile_in(crate::compression::archive::output_parent(&config.output))?
        .into_parts();
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, output_file);

    let mut flags = FLAG_CONSENSUS_DELTA;
    if opts.quality_reduction.is_some() {
        flags |= FLAG_QUALITY_LOSSY;
    }
    let header = ArchiveHeader {
        flags,
        num_records: 0, // placeholder — patched after all chunks
        num_chunks: 0,  // placeholder — patched after all chunks
        sam_header_compressed,
    };
    header.write_to(&mut writer)?;

    let use_fqz_quality = opts.quality_compressor == 0;
    let pipeline_window = resolve_pipeline_window(opts.compress_window, opts.chunk_size);
    if opts.compress_window == 0 {
        info!(
            "Auto-selected compress_window={} ({} cores, {})",
            pipeline_window,
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(8),
            match available_memory_bytes() {
                Some(b) => format!("{} GB available", b / (1024 * 1024 * 1024)),
                None => "RAM unknown".to_string(),
            },
        );
    }

    if pipeline_window > 1 {
        info!(
            "Parallel chunk pipeline: window={} chunks compressing simultaneously",
            pipeline_window
        );
    }

    // Compress every chunk through the shared windowed pipeline. `None` = no chunk
    // budget (read to EOF). It writes each ChunkHeader + streams to `writer` in
    // chunk order and returns the running totals (patched into the header below).
    // The NUMA range worker (`compress_chunk_range`) drives the SAME pipeline with
    // a pre-seeked reader and a bounded chunk budget — so a shard's chunk bytes are
    // identical to the single-process run's.
    let (total_records, num_chunks) = run_compress_pipeline(
        &mut reader,
        opts,
        use_fqz_quality,
        pipeline_window,
        None,
        &mut writer,
    )?;

    info!("Read {} records in {} chunks", total_records, num_chunks);

    // Patch num_records and num_chunks in the archive header.
    // Layout: [2B magic][1B version][1B reserved][1B flags][8B num_records][4B num_chunks]
    const NUM_RECORDS_OFFSET: u64 = 2 + 1 + 1 + 1; // = 5
    const NUM_CHUNKS_OFFSET: u64 = NUM_RECORDS_OFFSET + 8; // = 13

    writer.flush()?;
    let file = writer.get_mut();
    file.seek(SeekFrom::Start(NUM_RECORDS_OFFSET))?;
    file.write_all(&total_records.to_le_bytes())?;
    file.seek(SeekFrom::Start(NUM_CHUNKS_OFFSET))?;
    file.write_all(&num_chunks.to_le_bytes())?;
    file.flush()?;
    drop(writer);

    // Atomically rename the temp into place. `persist` consumes the `TempPath`,
    // so its Drop no longer removes the file; on failure the error carries the
    // `TempPath` and cleanup still happens when it drops.
    temp_path
        .persist(&config.output)
        .map_err(|e| anyhow::anyhow!("Failed to finalize output {:?}: {}", config.output, e))?;

    let output_size = std::fs::metadata(&config.output)?.len();
    let elapsed = start_time.elapsed().as_secs_f64();
    info!(
        "Compressed {} records to {} bytes ({:.2} MB) in {:.2}s",
        total_records,
        output_size,
        output_size as f64 / (1024.0 * 1024.0),
        elapsed
    );

    Ok(())
}

/// Names of the first 15 columnar streams, for per-chunk compression logging.
/// (The remaining `NUM_STREAMS - 15` streams are v8–v12 derivation-metadata
/// side-channels, not individually logged. Logging never affects archive bytes.)
const STREAM_NAMES: [&str; 15] = [
    "consensus",
    "ref_id",
    "pos",
    "mapq",
    "bin",
    "flag",
    "next_ref_id",
    "next_pos",
    "tlen",
    "read_name",
    "cigar",
    "seq_diff",
    "seq_extra",
    "qual",
    "aux",
];

/// Spawn one chunk-compression thread. `idx` is log-only — it never enters the
/// encoding — so a chunk's bytes depend solely on its records, independent of
/// where in the file a worker starts (the property NUMA sharding relies on).
fn spawn_chunk_thread(
    records: Vec<RawBamRecord>,
    opts: &AdvancedOptions,
    use_fqz_quality: bool,
    idx: u32,
) -> std::thread::JoinHandle<Result<ChunkResult>> {
    let opts_clone = opts.clone();
    std::thread::spawn(move || compress_one_chunk(records, opts_clone, use_fqz_quality, idx))
}

/// Write one completed chunk (its `ChunkHeader` + 15 compressed streams) to
/// `writer`. Returns the chunk's record count.
fn write_chunk_result<W: Write>(result: ChunkResult, writer: &mut W) -> Result<u32> {
    let ChunkResult {
        compressed,
        num_records: n,
        chunk_flags,
        raw_sizes,
        content_crc32,
    } = result;

    for (i, name) in STREAM_NAMES.iter().enumerate() {
        let raw = raw_sizes[i];
        let comp = compressed[i].len();
        let ratio = if comp > 0 { raw as f64 / comp as f64 } else { 0.0 };
        info!(
            "  {:>12}: {:>10} raw -> {:>10} compressed ({:.2}x)",
            name, raw, comp, ratio
        );
    }

    let mut stream_sizes = [0u32; NUM_STREAMS];
    let mut crc = Crc::new();
    for (i, cs) in compressed.iter().enumerate() {
        if cs.len() > u32::MAX as usize {
            anyhow::bail!("Compressed stream {} exceeds 4GB ({} bytes)", i, cs.len());
        }
        stream_sizes[i] = cs.len() as u32;
        crc.update(cs);
    }
    let chunk_header = ChunkHeader {
        num_records: n as u32,
        chunk_flags,
        crc32: crc.sum(),
        content_crc32,
        stream_sizes,
    };
    chunk_header.write_to(writer)?;
    for cs in &compressed {
        writer.write_all(cs)?;
    }
    Ok(n as u32)
}

/// The shared windowed chunk-compression pipeline, driving up to `pipeline_window`
/// `compress_one_chunk` threads concurrently and writing their `ChunkHeader` +
/// streams to `writer` in chunk order. Reads from `reader` until EOF or until
/// `max_chunks` chunks have been dispatched (`None` = no budget). Returns
/// `(total_records, num_chunks)`.
///
/// `pending` is always drained (joined) on any early `?`, so a mid-pipeline error
/// never leaves BSC worker threads running detached. Single-process `compress()`
/// drives this with `max_chunks = None`; the NUMA range worker
/// [`compress_chunk_range`] drives it with a pre-seeked reader and a bounded
/// budget — producing chunks byte-identical to the full run.
fn run_compress_pipeline<R: Read, W: Write>(
    reader: &mut RawBamReader<R>,
    opts: &AdvancedOptions,
    use_fqz_quality: bool,
    pipeline_window: usize,
    max_chunks: Option<u32>,
    writer: &mut W,
) -> Result<(u64, u32)> {
    type Handle = std::thread::JoinHandle<Result<ChunkResult>>;
    let budget = max_chunks.unwrap_or(u32::MAX);
    let chunk_size = opts.chunk_size;
    let mut pending: VecDeque<Handle> = VecDeque::new();
    let mut dispatch_idx: u32 = 0; // chunks dispatched so far (also the log index)
    let mut total_records: u64 = 0;
    let mut num_chunks: u32 = 0;

    let pipeline_result: Result<()> = (|| {
        // Prime: fill the window, bounded by the chunk budget and input EOF.
        while dispatch_idx < budget && pending.len() < pipeline_window {
            let records = reader.read_chunk(chunk_size)?;
            if records.is_empty() {
                break;
            }
            pending.push_back(spawn_chunk_thread(records, opts, use_fqz_quality, dispatch_idx));
            dispatch_idx += 1;
        }
        // Steady state: write the front completed chunk, then refill one if the
        // budget and input still allow. The `while let` naturally drains the tail.
        while let Some(handle) = pending.pop_front() {
            let result = handle
                .join()
                .unwrap_or_else(|_| Err(anyhow::anyhow!("chunk compression thread panicked")))?;
            total_records += write_chunk_result(result, writer)? as u64;
            num_chunks += 1;

            if dispatch_idx < budget {
                let records = reader.read_chunk(chunk_size)?;
                if !records.is_empty() {
                    pending.push_back(spawn_chunk_thread(
                        records,
                        opts,
                        use_fqz_quality,
                        dispatch_idx,
                    ));
                    dispatch_idx += 1;
                }
            }
        }
        Ok(())
    })();

    // Always join any remaining in-flight threads before propagating an error.
    for h in pending.drain(..) {
        let _ = h.join();
    }
    pipeline_result?;
    Ok((total_records, num_chunks))
}

/// Per-chunk layout of a BAM under a given compress config, for NUMA-compress
/// sharding. Built by [`read_bam_compress_layout`].
pub struct BzCompressLayout {
    /// Total BAM records.
    pub num_records: u64,
    /// Total chunks under this config's `chunk_size`.
    pub num_chunks: u32,
    /// Resolved records-per-chunk (after level presets).
    pub chunk_size: usize,
    /// BGZF virtual position at the start of each chunk (`len == num_chunks`).
    /// `chunk_start_vpos[k]` is where a worker seeks to begin compressing chunk k.
    pub chunk_start_vpos: Vec<u64>,
    /// Record count of each chunk (`len == num_chunks`; all == `chunk_size` except
    /// possibly the last). Feeds the gate's balanced range partitioner.
    pub per_chunk_reads: Vec<u64>,
}

/// Prescan a BAM once (multithreaded inflate, NO compression), counting records
/// and capturing the BGZF virtual position at every `chunk_size`-record boundary.
/// The NUMA-compress driver uses this to split the input into contiguous chunk
/// ranges that each worker seeks to and compresses independently. Cost is one
/// inflate+parse pass over the input — the serial part of a compress shard; peak
/// RSS is one chunk's records (each chunk is dropped before the next is read).
pub fn read_bam_compress_layout(input: &Path, config: &CompressConfig) -> Result<BzCompressLayout> {
    config.advanced.validate()?;
    let opts = apply_level(config.advanced.clone());
    let chunk_size = opts.chunk_size;

    let bgzf_workers = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::new(1).unwrap())
        .min(NonZeroUsize::new(16).unwrap());
    let mut reader = RawBamReader::from_path_mt(input, bgzf_workers)?;

    // Cluster/consensus encoding requires coordinate sort; fail fast here (the
    // workers skip this re-check, trusting the driver's prescan gate).
    if let SortOrder::Other(so) = crate::io::bam::parse_sort_order(&reader.header_raw) {
        anyhow::bail!(
            "BZ requires a coordinate-sorted BAM file, but this file declares SO:{so}. \
             Sort it first with: samtools sort -o sorted.bam {}",
            input.display(),
        );
    }

    let mut chunk_start_vpos: Vec<u64> = Vec::new();
    let mut per_chunk_reads: Vec<u64> = Vec::new();
    let mut num_records: u64 = 0;
    loop {
        let vpos = reader.virtual_position();
        let records = reader.read_chunk(chunk_size)?;
        if records.is_empty() {
            break;
        }
        chunk_start_vpos.push(vpos);
        per_chunk_reads.push(records.len() as u64);
        num_records += records.len() as u64;
    }
    let num_chunks = u32::try_from(chunk_start_vpos.len())
        .map_err(|_| anyhow::anyhow!("too many chunks ({})", chunk_start_vpos.len()))?;
    Ok(BzCompressLayout {
        num_records,
        num_chunks,
        chunk_size,
        chunk_start_vpos,
        per_chunk_reads,
    })
}

/// Conservative per-worker peak compress RSS for a NUMA shard, used to filter
/// eligible nodes in the gate (`mem_filtered_nodes`). This is the SOLE OOM guard
/// the gate applies for compress: a re-exec'd worker resolves a FIXED level preset,
/// and `apply_level` always pins a non-zero `compress_window`, so
/// `resolve_pipeline_window` never reads `/proc/meminfo` — the worker does NOT
/// self-limit its window by RAM, so under-budgeting here can genuinely
/// over-subscribe a node. RSS is bounded by the window (concurrent in-flight
/// chunks), not the chunk count, so this stays **constant in total read count**.
pub fn compress_peak_rss_bound(num_records: u64, num_chunks: u32) -> u64 {
    // The per-process base (BGZF reader buffers + mimalloc pool retention) is paid
    // IN FULL by every worker process — do NOT halve it. Only the window term is a
    // per-worker share (a 2-way worker holds ~half the chunks; RSS tracks the
    // concurrent window, not the chunk count).
    let per_worker_chunks = (num_chunks as u64).div_ceil(2).max(1);
    let window = (AUTO_WINDOW_CAP as u64).min(per_worker_chunks);
    let records_per_chunk = num_records / (num_chunks.max(1) as u64);
    COMPRESS_BASE_OVERHEAD_BYTES + window * records_per_chunk.max(1) * PER_RECORD_CHUNK_BYTES
}

/// Compress chunks `[chunk_start, chunk_end)` of `input` into `out_part`, seeking
/// to `start_vpos` (the `chunk_start` boundary from [`read_bam_compress_layout`]).
///
/// The part is self-contained: if `chunk_start == 0` it writes the global
/// `ArchiveHeader` carrying the WHOLE-archive `total_records`/`total_chunks` (so
/// no later header patching is needed) followed by its chunks; otherwise it writes
/// only the chunk stream. The driver concatenates parts in chunk order to form the
/// final archive, byte-identical to a single-process `compress`. Output is atomic
/// (temp file + rename).
///
/// `expected_records` is the prescan's record count for this range
/// (`sum(per_chunk_reads[chunk_start..chunk_end])`); the worker bails if it
/// compresses a different number. With a frozen level (see [`resolve_compress_level`])
/// the counts always match by construction, so this is a defensive cross-check that
/// turns any future chunk-size divergence — or an input changing under us — into a
/// hard error with NO output, rather than silently dropping or duplicating records.
#[allow(clippy::too_many_arguments)]
pub fn compress_chunk_range(
    input: &Path,
    chunk_start: u32,
    chunk_end: u32,
    start_vpos: u64,
    total_records: u64,
    total_chunks: u32,
    expected_records: u64,
    config: &CompressConfig,
    out_part: &Path,
) -> Result<u64> {
    if chunk_start >= chunk_end {
        anyhow::bail!("compress_chunk_range: empty range [{chunk_start}, {chunk_end})");
    }
    config.advanced.validate()?;
    let opts = apply_level(config.advanced.clone());
    let use_fqz_quality = opts.quality_compressor == 0;
    let pipeline_window = resolve_pipeline_window(opts.compress_window, opts.chunk_size);

    let bgzf_workers = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::new(1).unwrap())
        .min(NonZeroUsize::new(16).unwrap());
    // Chunk 0 starts right after the header (no seek); later chunks seek to their
    // captured virtual position while keeping multithreaded inflate.
    let mut reader = if chunk_start == 0 {
        RawBamReader::from_path_mt(input, bgzf_workers)?
    } else {
        RawBamReader::from_path_mt_seek(input, bgzf_workers, start_vpos)?
    };

    let (out_file, temp_path) = tempfile::Builder::new()
        .prefix(".bz-part-")
        .suffix(".tmp")
        .tempfile_in(crate::compression::archive::output_parent(out_part))?
        .into_parts();
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, out_file);

    // Part 0 carries the global archive header with the whole-archive totals.
    if chunk_start == 0 {
        let mut header_payload = Vec::new();
        let header_raw_len = u32::try_from(reader.header_raw.len()).map_err(|_| {
            anyhow::anyhow!("SAM header too large ({} bytes)", reader.header_raw.len())
        })?;
        header_payload.extend_from_slice(&header_raw_len.to_le_bytes());
        header_payload.extend_from_slice(&reader.header_raw);
        header_payload.extend_from_slice(&reader.ref_dict_bytes);
        let sam_header_compressed = bsc::compress_parallel(&header_payload)?;

        let mut flags = FLAG_CONSENSUS_DELTA;
        if opts.quality_reduction.is_some() {
            flags |= FLAG_QUALITY_LOSSY;
        }
        let header = ArchiveHeader {
            flags,
            num_records: total_records,
            num_chunks: total_chunks,
            sam_header_compressed,
        };
        header.write_to(&mut writer)?;
    }

    let max_chunks = chunk_end - chunk_start;
    let (records_written, chunks_written) = run_compress_pipeline(
        &mut reader,
        &opts,
        use_fqz_quality,
        pipeline_window,
        Some(max_chunks),
        &mut writer,
    )?;
    if chunks_written != max_chunks {
        anyhow::bail!(
            "compress_chunk_range: expected {max_chunks} chunks in [{chunk_start}, {chunk_end}), \
             got {chunks_written} (input shorter than the layout, or a stale layout)"
        );
    }
    if records_written != expected_records {
        anyhow::bail!(
            "compress_chunk_range: range [{chunk_start}, {chunk_end}) compressed {records_written} \
             records but the prescan layout expected {expected_records} — chunk_size divergence \
             between the prescan and this worker, or the input changed under us. Aborting with no \
             output to avoid silently dropping/duplicating records."
        );
    }

    writer.flush()?;
    drop(writer);
    temp_path
        .persist(out_part)
        .map_err(|e| anyhow::anyhow!("failed to finalize part {:?}: {}", out_part, e))?;
    Ok(records_written)
}

/// Compress quality scores with htscodecs fqzcomp, parallel per block (v12).
///
/// Two transforms vs the old fqz v2 path: (1) feed **raw Phred** (the
/// ASCII Phred+33 stored in `quality_data` minus 33), and (2) **reorient
/// reverse-strand reads** — reverse the quality of `FLAG&0x10` reads so all reads
/// run in sequencing-cycle order. fqzcomp's position/history model exploits the
/// now-consistent cycle trajectory; a probe measured −3.9% on full-resolution
/// WGS quality vs fqz v2 (lossless). Reorientation is reversed on decode
/// using the FLAG stream, so nothing extra is stored. The framing matches the v2
/// multiblock layout: [num_blocks][block_len, record_count, crc32, fqz_blob]...
fn compress_quality_fqz_parallel(
    quality_data: &ChunkQualityData,
    block_size: usize,
) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    let n = quality_data.qualities_ascii.len();
    if n == 0 {
        return Ok(Vec::new());
    }

    let block_ranges: Vec<(usize, usize)> = (0..n)
        .step_by(block_size)
        .map(|start| (start, (start + block_size).min(n)))
        .collect();

    // Compress blocks in parallel. Each block produces (record_count, blob).
    let compressed_blocks: Vec<Result<(u32, Vec<u8>)>> = block_ranges
        .par_iter()
        .map(|&(start, end)| {
            // Raw Phred, reverse-strand reads reoriented to cycle order.
            let reoriented: Vec<Vec<u8>> = (start..end)
                .map(|i| {
                    let mut raw: Vec<u8> = quality_data.qualities_ascii[i]
                        .iter()
                        .map(|&a| a.wrapping_sub(33))
                        .collect();
                    if quality_data.strands[i] {
                        raw.reverse();
                    }
                    raw
                })
                .collect();
            let refs: Vec<&[u8]> = reoriented.iter().map(|v| v.as_slice()).collect();
            let blob = fqzcomp::compress(&refs, 0)?;
            Ok(((end - start) as u32, blob))
        })
        .collect();

    // Multiblock framing: [num_blocks][block_len, record_count, crc32, fqz_blob].
    // The per-block CRC32 covers (record_count || fqz_blob); redundant with bz's
    // outer per-chunk CRC32 but kept for layout parity / defense in depth.
    use flate2::Crc;
    let mut result = Vec::new();
    result.extend_from_slice(&(compressed_blocks.len() as u32).to_le_bytes());
    for block_result in compressed_blocks {
        let (record_count, block) = block_result?;
        let crc = {
            let mut h = Crc::new();
            h.update(&record_count.to_le_bytes());
            h.update(&block);
            h.sum()
        };
        result.extend_from_slice(&(block.len() as u32).to_le_bytes());
        result.extend_from_slice(&record_count.to_le_bytes());
        result.extend_from_slice(&crc.to_le_bytes());
        result.extend_from_slice(&block);
    }

    Ok(result)
}

/// Compress the columnar streams in four parallel groups. When `skip_quality` is
/// set (the fqzcomp path), stream 13 (quality) is left as an empty placeholder
/// for the caller to replace with the fqzcomp blob; otherwise (the BSC fallback)
/// it's BSC-compressed inline with the aux stream.
fn compress_streams(
    streams: &BamStreams,
    opts: &AdvancedOptions,
    skip_quality: bool,
) -> Result<Vec<Vec<u8>>> {
    let slices = streams.as_slices();

    // Group A: consensus + ref_id + pos + mapq + bin + flag (0-5)
    // Group B: next_ref_id + next_pos + tlen (6-8)
    // Group C: read_name + cigar + seq_diff + seq_extra (9-12)
    // Group D: aux (14), or quality+aux (13-14) when not skipping quality.
    let group_d_start = if skip_quality { 14 } else { 13 };
    let ((group_a, group_b), (group_c, group_d)) = rayon::join(
        || {
            rayon::join(
                || compress_group_indexed(&slices[0..6], 0, opts),
                || compress_group_indexed(&slices[6..9], 6, opts),
            )
        },
        || {
            rayon::join(
                || compress_group_indexed(&slices[9..13], 9, opts),
                || compress_group_indexed(&slices[group_d_start..15], group_d_start, opts),
            )
        },
    );

    let mut result = Vec::with_capacity(NUM_STREAMS);
    for v in group_a? {
        result.push(v);
    }
    for v in group_b? {
        result.push(v);
    }
    for v in group_c? {
        result.push(v);
    }
    if skip_quality {
        result.push(Vec::new()); // slot 13 placeholder (replaced by the fqzcomp blob)
    }
    for v in group_d? {
        result.push(v);
    }
    // Derivation streams: MD (15-17), MC (18-19), tlen_derivable (20),
    // next_derivable (21), NM (nm_present 22, nm_index 23, nm_type 24).
    for v in compress_group_indexed(&slices[15..25], 15, opts)? {
        result.push(v);
    }

    Ok(result)
}

/// Compress a group of streams sequentially (within a rayon task). All streams
/// use BSC; `start_index` is retained for the call-site grouping structure.
fn compress_group_indexed(
    slices: &[&[u8]],
    _start_index: usize,
    opts: &AdvancedOptions,
) -> Result<Vec<Vec<u8>>> {
    let mut results = Vec::with_capacity(slices.len());
    for &data in slices.iter() {
        results.push(compress_one(data, opts)?);
    }
    Ok(results)
}

#[cfg(test)]
mod rss_bound_tests {
    use super::*;

    /// Regression guard: the per-process base must be budgeted IN FULL (it was once
    /// halved, under-budgeting each worker by ~2 GiB — the gate's sole compress OOM
    /// guard, since apply_level pins a fixed window so the worker never self-limits).
    #[test]
    fn compress_peak_rss_bound_budgets_full_base() {
        // Single chunk -> window collapses to 1, isolating the base term.
        let n_records = 1_000_000u64;
        let bound = compress_peak_rss_bound(n_records, 1);
        // window factor collapses to 1 for a single chunk.
        let window_term = n_records * PER_RECORD_CHUNK_BYTES;
        assert_eq!(
            bound,
            COMPRESS_BASE_OVERHEAD_BYTES + window_term,
            "bound must be FULL base + window term (no /2 on the base)"
        );
    }

    /// The window term is the per-worker chunk share (ceil(n/2)) capped by the
    /// window cap; the base stays full.
    #[test]
    fn compress_peak_rss_bound_window_is_per_worker_share() {
        let n_records = 15_000_000u64;
        let n_chunks = 10u32;
        let per_worker_chunks = (n_chunks as u64).div_ceil(2); // 5
        let window = (AUTO_WINDOW_CAP as u64).min(per_worker_chunks);
        let records_per_chunk = n_records / n_chunks as u64;
        let expect =
            COMPRESS_BASE_OVERHEAD_BYTES + window * records_per_chunk * PER_RECORD_CHUNK_BYTES;
        assert_eq!(compress_peak_rss_bound(n_records, n_chunks), expect);
    }

    /// resolve_compress_level is idempotent on an already-concrete level (1/2/3) —
    /// the property the NUMA driver relies on to freeze the level once.
    #[test]
    fn resolve_compress_level_is_idempotent_on_explicit_levels() {
        for lvl in 1u8..=3 {
            let cfg = CompressConfig {
                input: std::path::PathBuf::from("x.bam"),
                output: std::path::PathBuf::from("x.bz"),
                working_dir: std::path::PathBuf::from("."),
                advanced: AdvancedOptions {
                    level: lvl,
                    ..Default::default()
                },
            };
            assert_eq!(resolve_compress_level(&cfg), lvl);
        }
    }

    /// Auto (level 0) resolves to a concrete in-range level (1 or 2; never 0/3+).
    #[test]
    fn resolve_compress_level_auto_picks_concrete() {
        let cfg = CompressConfig {
            input: std::path::PathBuf::from("x.bam"),
            output: std::path::PathBuf::from("x.bz"),
            working_dir: std::path::PathBuf::from("."),
            advanced: AdvancedOptions {
                level: 0,
                ..Default::default()
            },
        };
        let lvl = resolve_compress_level(&cfg);
        assert!(lvl == 1 || lvl == 2, "auto must resolve to L1 or L2, got {lvl}");
    }
}
