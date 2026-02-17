use crate::cli::CompressConfig;
use crate::compression::archive::{ArchiveHeader, ChunkHeader, FLAG_CONSENSUS_DELTA, FLAG_QUALITY_CTX};
use crate::compression::streams::{self, BamStreams, ChunkQualityData, NUM_STREAMS};
use crate::io::bam::RawBamReader;
use anyhow::Result;
use qz_lib::compression::bsc;
use qz_lib::compression::quality_ctx;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::time::Instant;
use tracing::info;

/// Number of BAM records per compression chunk.
const CHUNK_SIZE: usize = 2_500_000;

/// Block size for quality_ctx (number of reads per block).
const QUALITY_CTX_BLOCK_SIZE: usize = 500_000;

/// Compress a BAM file to a BZ archive.
pub fn compress(config: &CompressConfig) -> Result<()> {
    let start_time = Instant::now();

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

    // Process chunks
    let mut total_records: u64 = 0;
    let mut chunk_data: Vec<(u32, Vec<Vec<u8>>)> = Vec::new();
    let mut use_quality_ctx = true;

    loop {
        let records = reader.read_chunk(CHUNK_SIZE)?;
        if records.is_empty() {
            break;
        }

        let n = records.len();
        info!("Chunk {}: {} records", chunk_data.len(), n);

        // Build columnar streams with consensus-delta encoding
        let result = streams::records_to_streams(&records);
        drop(records);

        let streams::StreamsResult { streams: bam_streams, quality_data } = result;

        // Check if quality_ctx is usable (no 0xFF quality bytes)
        if quality_data.has_unavailable_quality {
            use_quality_ctx = false;
        }

        let compressed = if use_quality_ctx && n > 0 {
            // Run BSC (skip quality) and quality_ctx in parallel
            let (bsc_result, qctx_result) = rayon::join(
                || compress_streams_skip_qual(&bam_streams),
                || compress_quality_ctx_parallel(&quality_data),
            );
            let mut compressed = bsc_result?;
            compressed[13] = qctx_result?;
            compressed
        } else {
            // BSC-only path
            compress_streams(&bam_streams)?
        };

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

        total_records += n as u64;
        chunk_data.push((n as u32, compressed));
    }

    info!(
        "Read {} records in {} chunks",
        total_records,
        chunk_data.len()
    );

    // Write archive
    let output_file = std::fs::File::create(&config.output)?;
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, output_file);

    let mut flags = FLAG_CONSENSUS_DELTA;
    if use_quality_ctx {
        flags |= FLAG_QUALITY_CTX;
    }

    let header = ArchiveHeader {
        flags,
        num_records: total_records,
        num_chunks: chunk_data.len() as u32,
        sam_header_compressed,
    };
    header.write_to(&mut writer)?;

    for (num_records, compressed_streams) in &chunk_data {
        let mut stream_sizes = [0u32; NUM_STREAMS];
        for (i, cs) in compressed_streams.iter().enumerate() {
            stream_sizes[i] = cs.len() as u32;
        }

        let chunk_header = ChunkHeader {
            num_records: *num_records,
            stream_sizes,
        };
        chunk_header.write_to(&mut writer)?;

        for cs in compressed_streams {
            writer.write_all(cs)?;
        }
    }

    writer.flush()?;

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
fn compress_quality_ctx_parallel(quality_data: &ChunkQualityData) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    let n = quality_data.qualities_ascii.len();
    if n == 0 {
        return Ok(Vec::new());
    }

    // Build block ranges
    let block_ranges: Vec<(usize, usize)> = (0..n)
        .step_by(QUALITY_CTX_BLOCK_SIZE)
        .map(|start| (start, (start + QUALITY_CTX_BLOCK_SIZE).min(n)))
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

/// Compress 14 streams with BSC (skip quality at index 13).
/// Quality is handled separately by quality_ctx.
fn compress_streams_skip_qual(streams: &BamStreams) -> Result<Vec<Vec<u8>>> {
    let slices = streams.as_slices();

    // Group A: consensus + ref_id + pos + mapq + bin + flag (indices 0-5)
    // Group B: next_ref_id + next_pos + tlen (indices 6-8)
    // Group C: read_name + cigar + seq_diff + seq_extra (indices 9-12)
    // Group D: aux only (index 14) — quality skipped
    let ((group_a, group_b), (group_c, group_d)) = rayon::join(
        || {
            rayon::join(
                || compress_group(&slices[0..6]),
                || compress_group(&slices[6..9]),
            )
        },
        || {
            rayon::join(
                || compress_group(&slices[9..13]),
                || compress_group(&slices[14..15]),  // aux only
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

/// Compress all 15 streams with BSC (BSC-only fallback path).
fn compress_streams(streams: &BamStreams) -> Result<Vec<Vec<u8>>> {
    let slices = streams.as_slices();

    let ((group_a, group_b), (group_c, group_d)) = rayon::join(
        || {
            rayon::join(
                || compress_group(&slices[0..6]),
                || compress_group(&slices[6..9]),
            )
        },
        || {
            rayon::join(
                || compress_group(&slices[9..13]),
                || compress_group(&slices[13..15]),
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
fn compress_group(slices: &[&[u8]]) -> Result<Vec<Vec<u8>>> {
    let mut results = Vec::with_capacity(slices.len());
    for &data in slices {
        let mut compressed = if data.is_empty() {
            Vec::new()
        } else {
            bsc::compress_parallel_adaptive_no_lzp(data)?
        };
        compressed.shrink_to_fit();
        results.push(compressed);
    }
    Ok(results)
}
