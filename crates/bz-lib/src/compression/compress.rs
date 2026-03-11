use crate::cli::{AdvancedOptions, CompressConfig};
use crate::compression::archive::{ArchiveHeader, ChunkHeader, CHUNK_FLAG_QUALITY_CTX, FLAG_CONSENSUS_DELTA};
use crate::compression::streams::{self, BamStreams, ChunkQualityData, NUM_STREAMS};
use crate::io::bam::RawBamReader;
use anyhow::Result;
use qz_lib::compression::bsc;
use qz_lib::compression::quality_ctx;
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
        // Select BSC variant based on config
        match (opts.bsc_adaptive, opts.use_lzp) {
            (true, true) => bsc::compress_parallel_adaptive(data)?,
            (true, false) => bsc::compress_parallel_adaptive_no_lzp(data)?,
            (false, true) => bsc::compress_parallel(data)?,
            (false, false) => bsc::compress_parallel(data)?, // static, no LZP
        }
    };
    compressed.shrink_to_fit();
    Ok(compressed)
}

/// Compress a BAM file to a BZ archive.
pub fn compress(config: &CompressConfig) -> Result<()> {
    let start_time = Instant::now();
    let opts = &config.advanced;

    info!(
        "BZ compression config: chunk_size={}, quality_ctx_block={}, \
         lzp={}, adaptive={}, qual_comp={}, align_comp={}, aux_comp={}",
        opts.chunk_size, opts.quality_ctx_block_size,
        opts.use_lzp, opts.bsc_adaptive, opts.quality_compressor,
        opts.alignment_compressor, opts.aux_compressor,
    );

    if config.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build_global()
            .ok(); // ignore if already initialized
    }

    // Use multi-threaded BGZF reader for parallel inflate
    let bgzf_workers = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::new(1).unwrap())
        .min(NonZeroUsize::new(4).unwrap());
    let mut reader = RawBamReader::from_path_mt(&config.input, bgzf_workers)?;

    // Compress SAM header + ref dict together
    let mut header_payload = Vec::new();
    let header_raw_len = reader.header_raw.len() as u32;
    header_payload.extend_from_slice(&header_raw_len.to_le_bytes());
    header_payload.extend_from_slice(&reader.header_raw);
    header_payload.extend_from_slice(&reader.ref_dict_bytes);
    let sam_header_compressed = bsc::compress_parallel(&header_payload)?;

    // Write to a temp file, then rename atomically on success
    let temp_output = config.output.with_extension("bz.tmp");
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

    // Stream chunks directly to disk
    let mut total_records: u64 = 0;
    let mut num_chunks: u32 = 0;
    let use_quality_ctx = opts.quality_compressor == 0;

    loop {
        let records = reader.read_chunk(opts.chunk_size)?;
        if records.is_empty() {
            break;
        }

        let n = records.len();
        info!("Chunk {}: {} records", num_chunks, n);

        // Build columnar streams with consensus-delta encoding
        let result = streams::records_to_streams(&records);
        drop(records);

        let streams::StreamsResult { streams: bam_streams, quality_data } = result;

        // Per-chunk decision: use quality_ctx only if enabled AND this chunk has no 0xFF quality
        let chunk_uses_quality_ctx = use_quality_ctx && !quality_data.has_unavailable_quality;

        let compressed = if chunk_uses_quality_ctx && n > 0 {
            // Run BSC (skip quality) and quality_ctx in parallel
            let qctx_block_size = opts.quality_ctx_block_size;
            let (bsc_result, qctx_result) = rayon::join(
                || compress_streams_skip_qual(&bam_streams, opts),
                || compress_quality_ctx_parallel(&quality_data, qctx_block_size),
            );
            let mut compressed = bsc_result?;
            compressed[13] = qctx_result?;
            compressed
        } else {
            // BSC-only path (includes quality)
            compress_all_streams(&bam_streams, opts)?
        };

        let mut chunk_flags = 0u8;
        if chunk_uses_quality_ctx {
            chunk_flags |= CHUNK_FLAG_QUALITY_CTX;
        }

        // Log per-stream sizes
        let stream_names = [
            "consensus", "ref_id", "pos", "mapq", "bin", "flag",
            "next_ref_id", "next_pos", "tlen",
            "read_name", "cigar", "seq_diff", "seq_extra", "qual", "aux",
        ];
        let raw_slices = bam_streams.as_slices();
        for (i, name) in stream_names.iter().enumerate() {
            let raw = raw_slices[i].len();
            let comp = compressed[i].len();
            let ratio = if comp > 0 { raw as f64 / comp as f64 } else { 0.0 };
            info!(
                "  {:>12}: {:>10} raw -> {:>10} compressed ({:.2}x)",
                name, raw, comp, ratio
            );
        }

        // Write chunk directly to output
        let mut stream_sizes = [0u32; NUM_STREAMS];
        for (i, cs) in compressed.iter().enumerate() {
            if cs.len() > u32::MAX as usize {
                anyhow::bail!(
                    "Compressed stream {} exceeds 4GB ({} bytes)", i, cs.len()
                );
            }
            stream_sizes[i] = cs.len() as u32;
        }
        let chunk_header = ChunkHeader {
            num_records: n as u32,
            chunk_flags,
            stream_sizes,
        };
        chunk_header.write_to(&mut writer)?;
        for cs in &compressed {
            writer.write_all(cs)?;
        }

        total_records += n as u64;
        num_chunks += 1;
    }

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

    // Atomic rename: temp file -> final output
    std::fs::rename(&temp_output, &config.output).map_err(|e| {
        // Clean up temp file on rename failure
        let _ = std::fs::remove_file(&temp_output);
        anyhow::anyhow!("Failed to rename temp file to {:?}: {}", config.output, e)
    })?;

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
