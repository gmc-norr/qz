use crate::cli::{AdvancedOptions, CompressConfig};
use crate::compression::archive::{ArchiveHeader, ChunkHeader, CHUNK_FLAG_QUALITY_CTX, FLAG_CONSENSUS_DELTA};
use crate::compression::streams::{self, BamStreams, ChunkQualityData, NUM_STREAMS};
use crate::io::bam::{RawBamReader, RawBamRecord};
use anyhow::Result;
use flate2::Crc;
use qz_lib::compression::bsc;
use qz_lib::compression::quality_ctx;
use std::collections::VecDeque;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::time::Instant;
use tracing::info;

/// Stream categories for compressor dispatch.
#[derive(Clone, Copy)]
enum StreamKind {
    /// Consensus data (index 0)
    Consensus,
    /// Alignment metadata: ref_id, pos, mapq, bin, flag, next_ref_id, next_pos, tlen (1-8)
    Alignment,
    /// Variable-length read data: read_name, cigar, seq_diff, seq_extra (9-12)
    ReadData,
    /// Quality scores (index 13)
    Quality,
    /// Auxiliary tags (index 14)
    Aux,
}

/// Map stream index to its kind.
fn stream_kind(index: usize) -> StreamKind {
    match index {
        0 => StreamKind::Consensus,
        1..=8 => StreamKind::Alignment,
        9..=12 => StreamKind::ReadData,
        13 => StreamKind::Quality,
        14 => StreamKind::Aux,
        _ => StreamKind::ReadData, // shouldn't happen
    }
}

/// Compress a single stream using the configured compressor for its kind.
fn compress_one(data: &[u8], kind: StreamKind, opts: &AdvancedOptions) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }

    // Determine which compressor to use based on stream kind and config
    let use_zstd = match kind {
        StreamKind::Alignment => opts.alignment_compressor == 1,
        StreamKind::Aux => opts.aux_compressor == 1,
        // Consensus, ReadData, Quality always use BSC
        _ => false,
    };

    let mut compressed = if use_zstd {
        zstd::bulk::compress(data, 3)?
    } else {
        // Select BSC variant based on config.
        // use_lzp only affects the adaptive path; the static path has no LZP knob.
        match (opts.bsc_adaptive, opts.use_lzp) {
            (true, true)  => bsc::compress_parallel_adaptive(data)?,
            (true, false) => bsc::compress_parallel_adaptive_no_lzp(data)?,
            (false, _)    => bsc::compress_parallel(data)?,
        }
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
    /// Per-chunk archive flags (e.g. CHUNK_FLAG_QUALITY_CTX).
    chunk_flags: u8,
    /// Raw (uncompressed) sizes per stream, for logging.
    raw_sizes: [usize; NUM_STREAMS],
}

/// Compress one chunk of BAM records.
///
/// Runs on the calling thread (or a spawned thread) and uses the global rayon
/// pool for stream-level parallelism.  Returns a [`ChunkResult`] containing the
/// 15 compressed streams plus metadata needed for archive assembly.
fn compress_one_chunk(
    records: Vec<RawBamRecord>,
    opts: AdvancedOptions,
    use_quality_ctx: bool,
    chunk_idx: u32,
) -> Result<ChunkResult> {
    let n = records.len();
    info!("Chunk {}: {} records (compressing)", chunk_idx, n);

    let result = streams::records_to_streams(&records);
    drop(records);

    let streams::StreamsResult { streams: bam_streams, quality_data } = result;

    // Capture raw sizes before consuming bam_streams
    let raw_slices = bam_streams.as_slices();
    let mut raw_sizes = [0usize; NUM_STREAMS];
    for (i, s) in raw_slices.iter().enumerate() {
        raw_sizes[i] = s.len();
    }

    // Per-chunk decision: use quality_ctx only if enabled AND no 0xFF quality
    let chunk_uses_quality_ctx = use_quality_ctx && !quality_data.has_unavailable_quality;

    let compressed = if chunk_uses_quality_ctx && n > 0 {
        let qctx_block_size = opts.quality_ctx_block_size;
        let (bsc_result, qctx_result) = rayon::join(
            || compress_streams_skip_qual(&bam_streams, &opts),
            || compress_quality_ctx_parallel(&quality_data, qctx_block_size),
        );
        let mut compressed = bsc_result?;
        compressed[13] = qctx_result?;
        compressed
    } else {
        compress_all_streams(&bam_streams, &opts)?
    };

    let mut chunk_flags = 0u8;
    if chunk_uses_quality_ctx {
        chunk_flags |= CHUNK_FLAG_QUALITY_CTX;
    }

    Ok(ChunkResult { compressed, num_records: n, chunk_flags, raw_sizes })
}

/// Compress a BAM file to a BZ archive.
pub fn compress(config: &CompressConfig) -> Result<()> {
    let start_time = Instant::now();
    let opts = &config.advanced;

    info!(
        "BZ compression config: chunk_size={}, quality_ctx_block={}, \
         lzp={}, adaptive={}, qual_comp={}, align_comp={}, aux_comp={}, window={}",
        opts.chunk_size, opts.quality_ctx_block_size,
        opts.use_lzp, opts.bsc_adaptive, opts.quality_compressor,
        opts.alignment_compressor, opts.aux_compressor, opts.compress_window,
    );

    // Use multi-threaded BGZF reader for parallel inflate
    let bgzf_workers = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::new(1).unwrap())
        .min(NonZeroUsize::new(4).unwrap());
    let mut reader = RawBamReader::from_path_mt(&config.input, bgzf_workers)?;

    // Compress SAM header + ref dict together
    let mut header_payload = Vec::new();
    let header_raw_len = u32::try_from(reader.header_raw.len())
        .map_err(|_| anyhow::anyhow!("SAM header too large ({} bytes)", reader.header_raw.len()))?;
    header_payload.extend_from_slice(&header_raw_len.to_le_bytes());
    header_payload.extend_from_slice(&reader.header_raw);
    header_payload.extend_from_slice(&reader.ref_dict_bytes);
    let sam_header_compressed = bsc::compress_parallel(&header_payload)?;

    // Write to a temp file, then rename atomically on success.
    // The drop guard removes the temp file on any early return (error or panic).
    let temp_output = config.output.with_extension("bz.tmp");
    struct TmpCleanup(std::path::PathBuf);
    impl Drop for TmpCleanup {
        fn drop(&mut self) { let _ = std::fs::remove_file(&self.0); }
    }
    let tmp_guard = TmpCleanup(temp_output.clone());

    let output_file = std::fs::File::create(&temp_output)?;
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, output_file);

    let flags = FLAG_CONSENSUS_DELTA;
    let header = ArchiveHeader {
        flags,
        num_records: 0,       // placeholder — patched after all chunks
        num_chunks: 0,        // placeholder — patched after all chunks
        alignment_compressor: opts.alignment_compressor,
        aux_compressor: opts.aux_compressor,
        sam_header_compressed,
    };
    header.write_to(&mut writer)?;

    let mut total_records: u64 = 0;
    let mut num_chunks: u32 = 0;  // written chunk count (patched into archive header)
    let use_quality_ctx = opts.quality_compressor == 0;
    let pipeline_window = opts.compress_window.max(1);

    if pipeline_window > 1 {
        info!(
            "Parallel chunk pipeline: window={} chunks compressing simultaneously",
            pipeline_window
        );
    }

    // ---------------------------------------------------------------------------
    // Windowed parallel pipeline
    //
    // Each chunk is dispatched as a std::thread::spawn'd task that compresses all
    // 15 streams using the rayon pool.  Up to `pipeline_window` chunks are in
    // flight simultaneously, keeping more rayon threads busy.
    //
    // Results are drained from the front of the deque and written in order, so the
    // archive byte layout is deterministic regardless of window size.
    // ---------------------------------------------------------------------------

    type Handle = std::thread::JoinHandle<Result<ChunkResult>>;
    let mut pending: VecDeque<Handle> = VecDeque::new();
    let mut dispatch_idx: u32 = 0;  // logical index for log messages only

    // Spawn one chunk compression thread, consuming `records`.
    let mut spawn_chunk = |records: Vec<RawBamRecord>, pending: &mut VecDeque<Handle>| {
        let opts_clone = opts.clone();
        let idx = dispatch_idx;
        dispatch_idx += 1;
        let handle: Handle = std::thread::spawn(move || {
            compress_one_chunk(records, opts_clone, use_quality_ctx, idx)
        });
        pending.push_back(handle);
    };

    // Write one completed ChunkResult to the output file.
    // Updates total_records and num_chunks.
    let stream_names = [
        "consensus", "ref_id", "pos", "mapq", "bin", "flag",
        "next_ref_id", "next_pos", "tlen",
        "read_name", "cigar", "seq_diff", "seq_extra", "qual", "aux",
    ];
    macro_rules! write_chunk {
        ($result:expr) => {{
            let ChunkResult { compressed, num_records: n, chunk_flags, raw_sizes } = $result;

            for (i, name) in stream_names.iter().enumerate() {
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
                    anyhow::bail!(
                        "Compressed stream {} exceeds 4GB ({} bytes)", i, cs.len()
                    );
                }
                stream_sizes[i] = cs.len() as u32;
                crc.update(cs);
            }
            let chunk_header = ChunkHeader {
                num_records: n as u32,
                chunk_flags,
                crc32: crc.sum(),
                stream_sizes,
            };
            chunk_header.write_to(&mut writer)?;
            for cs in &compressed {
                writer.write_all(cs)?;
            }

            total_records += n as u64;
            num_chunks += 1;
        }};
    }

    // ── Run the pipeline in a closure so `?` propagates from the closure, not
    // from `compress()`.  This guarantees that `pending` is always drained after
    // the closure exits — whether it succeeded or failed — so no BSC threads are
    // left running detached when an error short-circuits the pipeline.
    let pipeline_result: Result<()> = (|| {
        // ── Prime ──────────────────────────────────────────────────────────
        let first_records = reader.read_chunk(opts.chunk_size)?;
        if first_records.is_empty() {
            return Ok(());
        }
        spawn_chunk(first_records, &mut pending);

        for _ in 1..pipeline_window {
            let records = reader.read_chunk(opts.chunk_size)?;
            if records.is_empty() {
                break;
            }
            spawn_chunk(records, &mut pending);
        }

        // ── Steady-state ───────────────────────────────────────────────────
        loop {
            let result = pending.pop_front().expect("pending non-empty")
                .join().unwrap_or_else(|_| Err(anyhow::anyhow!("chunk compression thread panicked")))?;
            write_chunk!(result);

            let records = reader.read_chunk(opts.chunk_size)?;
            if records.is_empty() {
                break;
            }
            spawn_chunk(records, &mut pending);
        }

        // ── Drain remaining ────────────────────────────────────────────────
        while !pending.is_empty() {
            let result = pending.pop_front().expect("pending non-empty")
                .join().unwrap_or_else(|_| Err(anyhow::anyhow!("chunk compression thread panicked")))?;
            write_chunk!(result);
        }

        Ok(())
    })();

    // Always join any remaining in-flight threads before propagating an error,
    // so BSC threads don't continue burning CPU after a failure.
    for h in pending.drain(..) {
        let _ = h.join();
    }

    pipeline_result?;

    info!(
        "Read {} records in {} chunks",
        total_records, num_chunks
    );

    // Patch num_records and num_chunks in the archive header.
    // Layout: [2B magic][1B version][1B reserved][1B flags][1B align_comp][1B aux_comp][8B num_records][4B num_chunks]
    const NUM_RECORDS_OFFSET: u64 = 2 + 1 + 1 + 1 + 1 + 1; // = 7
    const NUM_CHUNKS_OFFSET: u64 = NUM_RECORDS_OFFSET + 8;   // = 15

    writer.flush()?;
    let file = writer.get_mut();
    file.seek(SeekFrom::Start(NUM_RECORDS_OFFSET))?;
    file.write_all(&total_records.to_le_bytes())?;
    file.seek(SeekFrom::Start(NUM_CHUNKS_OFFSET))?;
    file.write_all(&num_chunks.to_le_bytes())?;
    file.flush()?;
    drop(writer);

    // Disarm the drop guard — rename takes ownership of the file.
    // If rename fails we still want cleanup, so re-arm by forgetting the disarm.
    std::mem::forget(tmp_guard);
    if let Err(e) = std::fs::rename(&temp_output, &config.output) {
        let _ = std::fs::remove_file(&temp_output);
        return Err(anyhow::anyhow!("Failed to rename temp file to {:?}: {}", config.output, e));
    }

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

/// Compress quality scores using quality_ctx with parallel block compression.
fn compress_quality_ctx_parallel(quality_data: &ChunkQualityData, block_size: usize) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    let n = quality_data.qualities_ascii.len();
    if n == 0 {
        return Ok(Vec::new());
    }

    // Build block ranges
    let block_ranges: Vec<(usize, usize)> = (0..n)
        .step_by(block_size)
        .map(|start| (start, (start + block_size).min(n)))
        .collect();

    // Compress blocks in parallel
    let compressed_blocks: Vec<Result<Vec<u8>>> = block_ranges
        .par_iter()
        .map(|&(start, end)| {
            let qual_refs: Vec<&[u8]> = quality_data.qualities_ascii[start..end]
                .iter()
                .map(|q| q.as_slice())
                .collect();
            let seq_refs: Vec<&[u8]> = quality_data.sequences_ascii[start..end]
                .iter()
                .map(|s| s.as_slice())
                .collect();
            quality_ctx::compress_qualities_ctx(&qual_refs, &seq_refs)
        })
        .collect();

    // Combine into multiblock format
    let mut result = Vec::new();
    result.extend_from_slice(&(compressed_blocks.len() as u32).to_le_bytes());
    for block_result in compressed_blocks {
        let block = block_result?;
        result.extend_from_slice(&(block.len() as u32).to_le_bytes());
        result.extend_from_slice(&block);
    }

    Ok(result)
}

/// Compress 14 streams (skip quality at index 13, handled by quality_ctx).
fn compress_streams_skip_qual(streams: &BamStreams, opts: &AdvancedOptions) -> Result<Vec<Vec<u8>>> {
    let slices = streams.as_slices();

    // Group A: consensus + ref_id + pos + mapq + bin + flag (indices 0-5)
    // Group B: next_ref_id + next_pos + tlen (indices 6-8)
    // Group C: read_name + cigar + seq_diff + seq_extra (indices 9-12)
    // Group D: aux only (index 14) — quality skipped
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
                || compress_group_indexed(&slices[14..15], 14, opts),
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
    // Slot 13 = quality placeholder (will be replaced by quality_ctx)
    result.push(Vec::new());
    for v in group_d? {
        result.push(v);
    }

    Ok(result)
}

/// Compress all 15 streams (BSC-only fallback path when quality_ctx is disabled).
fn compress_all_streams(streams: &BamStreams, opts: &AdvancedOptions) -> Result<Vec<Vec<u8>>> {
    let slices = streams.as_slices();

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
                || compress_group_indexed(&slices[13..15], 13, opts),
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
    for v in group_d? {
        result.push(v);
    }

    Ok(result)
}

/// Compress a group of streams sequentially (within a rayon task).
/// `start_index` is the global stream index of the first slice.
fn compress_group_indexed(slices: &[&[u8]], start_index: usize, opts: &AdvancedOptions) -> Result<Vec<Vec<u8>>> {
    let mut results = Vec::with_capacity(slices.len());
    for (i, &data) in slices.iter().enumerate() {
        let kind = stream_kind(start_index + i);
        let compressed = compress_one(data, kind, opts)?;
        results.push(compressed);
    }
    Ok(results)
}
