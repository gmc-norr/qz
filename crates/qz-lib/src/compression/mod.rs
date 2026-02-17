mod columnar;
mod n_mask;
mod quality;
mod quality_delta;
mod quality_model;
mod read_id;
mod zstd_dict;
pub mod arithmetic_quality;
pub mod arithmetic_sequence;
pub mod bsc;
pub mod debruijn;
pub mod dna_utils;
pub mod fqzcomp;
pub mod greedy_contig;
pub mod quality_context;
pub mod openzl;
pub mod template;
mod ultra;
pub mod quality_ctx;
pub mod header_col;
mod codecs;
mod compress_impl;
mod decompress_impl;

use crate::cli::{CompressConfig, DecompressConfig, QualityMode, QualityCompressor, SequenceCompressor, HeaderCompressor};
use anyhow::{Context, Result};
use std::time::Instant;
use tracing::info;

/// Number of BSC blocks to decompress in each parallel batch
const DECOMPRESS_BATCH_SIZE: usize = 8;

/// Minimum reads for quality_ctx to outperform BSC (below this, models don't converge)
const MIN_READS_QUALITY_CTX: usize = 100_000;

/// QZ archive magic bytes (file identification)
const ARCHIVE_MAGIC: [u8; 2] = *b"QZ";
/// Current archive format version
const ARCHIVE_VERSION: u8 = 2;
/// Size of the v2 prefix: magic(2) + version(1) + reserved(1) + header_size(4)
const V2_PREFIX_SIZE: usize = 8;
/// I/O buffer size for reading archive files during decompression
const IO_BUFFER_SIZE: usize = 8 * 1024 * 1024; // 8 MB

/// Read a little-endian u16 from `data` at `offset`, returning an error on truncation.
fn read_le_u16(data: &[u8], offset: usize) -> Result<u16> {
    data.get(offset..offset + 2)
        .and_then(|s| <[u8; 2]>::try_from(s).ok())
        .map(u16::from_le_bytes)
        .ok_or_else(|| anyhow::anyhow!("truncated archive at offset {offset}"))
}

/// Read a little-endian u32 from `data` at `offset`, returning an error on truncation.
fn read_le_u32(data: &[u8], offset: usize) -> Result<u32> {
    data.get(offset..offset + 4)
        .and_then(|s| <[u8; 4]>::try_from(s).ok())
        .map(u32::from_le_bytes)
        .ok_or_else(|| anyhow::anyhow!("truncated archive at offset {offset}"))
}

/// Read a little-endian u64 from `data` at `offset`, returning an error on truncation.
fn read_le_u64(data: &[u8], offset: usize) -> Result<u64> {
    data.get(offset..offset + 8)
        .and_then(|s| <[u8; 8]>::try_from(s).ok())
        .map(u64::from_le_bytes)
        .ok_or_else(|| anyhow::anyhow!("truncated archive at offset {offset}"))
}

/// Convert QualityMode to QualityBinning
fn quality_mode_to_binning(mode: QualityMode) -> QualityBinning {
    match mode {
        QualityMode::Lossless => QualityBinning::None, // None = lossless (keep original)
        QualityMode::IlluminaBin => QualityBinning::Illumina8,
        QualityMode::Illumina4 => QualityBinning::FourLevel,
        QualityMode::Binary => QualityBinning::Binary { threshold: 20 },
        QualityMode::Qvz => QualityBinning::None, // TODO: Implement QVZ
        QualityMode::Discard => QualityBinning::Binary { threshold: 255 }, // All qualities -> same value
    }
}

/// Convert QualityBinning to a byte code for storage
fn binning_to_code(binning: QualityBinning) -> u8 {
    match binning {
        QualityBinning::None => 0, // Lossless
        QualityBinning::Illumina8 => 1,
        QualityBinning::Binary { .. } => 2,
        QualityBinning::FourLevel => 3,
    }
}

/// Convert byte code back to QualityBinning
fn code_to_binning(code: u8) -> Result<QualityBinning> {
    match code {
        0 => Ok(QualityBinning::None), // Lossless
        1 => Ok(QualityBinning::Illumina8),
        2 => Ok(QualityBinning::Binary { threshold: 20 }),
        3 => Ok(QualityBinning::FourLevel),
        _ => anyhow::bail!("Invalid quality binning code: {}", code),
    }
}

/// Convert QualityCompressor to a byte code for storage
fn compressor_to_code(compressor: QualityCompressor) -> u8 {
    match compressor {
        QualityCompressor::Zstd => 0,
        QualityCompressor::Bsc => 1,
        QualityCompressor::OpenZl => 2,
        QualityCompressor::Fqzcomp => 3,
        QualityCompressor::QualityCtx => 4,
    }
}

/// Convert byte code back to QualityCompressor
fn code_to_compressor(code: u8) -> Result<QualityCompressor> {
    match code {
        0 => Ok(QualityCompressor::Zstd),
        1 => Ok(QualityCompressor::Bsc),
        2 => Ok(QualityCompressor::OpenZl),
        3 => Ok(QualityCompressor::Fqzcomp),
        4 => Ok(QualityCompressor::QualityCtx),
        _ => anyhow::bail!("Invalid quality compressor code: {}", code),
    }
}

/// Convert SequenceCompressor to a byte code for storage
pub(crate) fn seq_compressor_to_code(compressor: SequenceCompressor) -> u8 {
    match compressor {
        SequenceCompressor::Zstd => 0,
        SequenceCompressor::Bsc => 1,
        SequenceCompressor::OpenZl => 2,
    }
}

/// Convert byte code back to SequenceCompressor
fn code_to_seq_compressor(code: u8) -> Result<SequenceCompressor> {
    match code {
        0 => Ok(SequenceCompressor::Zstd),
        1 => Ok(SequenceCompressor::Bsc),
        2 => Ok(SequenceCompressor::OpenZl),
        _ => anyhow::bail!("Invalid sequence compressor code: {}", code),
    }
}

/// Convert HeaderCompressor to a byte code for storage
pub(crate) fn header_compressor_to_code(compressor: HeaderCompressor) -> u8 {
    match compressor {
        HeaderCompressor::Zstd => 0,
        HeaderCompressor::Bsc => 1,
        HeaderCompressor::OpenZl => 2,
    }
}

/// Convert byte code back to HeaderCompressor
fn code_to_header_compressor(code: u8) -> Result<HeaderCompressor> {
    match code {
        0 => Ok(HeaderCompressor::Zstd),
        1 => Ok(HeaderCompressor::Bsc),
        2 => Ok(HeaderCompressor::OpenZl),
        _ => anyhow::bail!("Invalid header compressor code: {}", code),
    }
}

pub use columnar::{
    compress_columnar,
    QualityBinning,
};

// Advanced and tier2 features not currently exposed in public API

/// Dispatch BSC parallel compression: static (--bsc-static) or adaptive-no-LZP (default).
/// Default uses adaptive QLFC without LZP — LZP doesn't help DNA/quality data
/// (BWT already captures patterns; SPRING also disables LZP).
fn bsc_compress_parallel(data: &[u8], use_static: bool) -> Result<Vec<u8>> {
    if use_static {
        bsc::compress_parallel(data)
    } else {
        bsc::compress_parallel_adaptive_no_lzp(data)
    }
}

/// Encode a sequence into seq_stream using inline delta against a cache of previously seen reads.
/// If a similar read exists in the cache (same min syncmer hash, >60% positional identity,
/// same length), encodes as delta (0x00=match, actual base=mismatch). Otherwise emits raw.
/// Stream format: [varint(len)] [flag: 0x00=raw, 0x01=delta] [ref_idx varint if delta] [bytes]
fn encode_sequence_delta(
    seq: &[u8],
    seq_stream: &mut Vec<u8>,
    cache: &mut std::collections::HashMap<u64, usize>,
    reads: &mut Vec<Vec<u8>>,
) -> Result<bool> {
    let hash = dna_utils::compute_min_syncmer_hash(seq);
    let current_idx = reads.len();

    let mut used_delta = false;

    if hash != u64::MAX {
        if let Some(&ref_idx) = cache.get(&hash) {
            let ref_seq = &reads[ref_idx];
            if ref_seq.len() == seq.len() {
                // Count matching bases
                let matches = seq.iter().zip(ref_seq.iter()).filter(|(a, b)| a == b).count();
                let identity = matches as f64 / seq.len() as f64;

                if identity > 0.6 {
                    // Delta encode: [varint(len)] [0x01] [ref_idx varint] [delta bytes]
                    write_varint(seq_stream, seq.len())?;
                    seq_stream.push(0x01);
                    write_varint(seq_stream, ref_idx)?;
                    for (i, &base) in seq.iter().enumerate() {
                        if base == ref_seq[i] {
                            seq_stream.push(0x00); // match
                        } else {
                            seq_stream.push(base); // mismatch: actual base
                        }
                    }
                    used_delta = true;
                }
            }
        }
    }

    if !used_delta {
        // Raw encode: [varint(len)] [0x00] [raw bases]
        write_varint(seq_stream, seq.len())?;
        seq_stream.push(0x00);
        seq_stream.extend_from_slice(seq);
    }

    // Always update cache and store the read
    if hash != u64::MAX {
        cache.insert(hash, current_idx);
    }
    reads.push(seq.to_vec());

    Ok(used_delta)
}

/// Read a chunk of FASTQ records into a Vec.
/// Returns (records, total_bases, original_size).
fn read_chunk_records<R: std::io::BufRead>(
    reader: &mut crate::io::FastqReader<R>,
    chunk_size: usize,
) -> Result<(Vec<crate::io::FastqRecord>, usize, usize)> {
    let mut records = Vec::with_capacity(chunk_size.min(1_000_000));
    let mut chunk_bases = 0usize;
    let mut chunk_orig = 0usize;

    for _ in 0..chunk_size {
        match reader.next()? {
            Some(record) => {
                chunk_bases += record.sequence.len();
                chunk_orig += record.id.len() + record.sequence.len()
                    + record.quality.as_ref().map(|q| q.len()).unwrap_or(0) + 3;
                records.push(record);
            }
            None => break,
        }
    }

    Ok((records, chunk_bases, chunk_orig))
}

/// Convert a slice of FastqRecords into raw compression streams.
/// Returns (header_stream, seq_stream, qual_stream).
/// Convert records to compression streams.
/// Returns (header_stream, seq_stream, qual_stream, rc_flags).
/// rc_flags is empty unless rc_canon=true.
fn records_to_streams(
    records: &[crate::io::FastqRecord],
    quality_mode: QualityMode,
    quality_binning: QualityBinning,
    no_quality: bool,
    sequence_hints: bool,
    sequence_delta: bool,
    rc_canon: bool,
    const_seq_len: usize,
    const_qual_len: usize,
) -> Result<(Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>)> {
    let n = records.len();
    let mut header_stream = Vec::with_capacity(n * 64);
    let mut seq_stream = Vec::with_capacity(n * 152);
    let mut qual_stream = Vec::with_capacity(n * 136);
    let mut rc_flags: Vec<u8> = if rc_canon { Vec::with_capacity(n) } else { Vec::new() };

    let mut delta_cache = std::collections::HashMap::new();
    let mut delta_reads = Vec::new();

    for record in records {
        write_varint(&mut header_stream, record.id.len())?;
        header_stream.extend_from_slice(&record.id);

        if sequence_delta {
            encode_sequence_delta(&record.sequence, &mut seq_stream, &mut delta_cache, &mut delta_reads)?;
        } else {
            if const_seq_len == 0 {
                write_varint(&mut seq_stream, record.sequence.len())?;
            }
            if sequence_hints {
                seq_stream.push(dna_utils::compute_sequence_hint(&record.sequence));
            }
            if rc_canon {
                let (canon_seq, was_reversed) = dna_utils::canonicalize_sequence(&record.sequence);
                seq_stream.extend_from_slice(&canon_seq);
                rc_flags.push(was_reversed as u8);
            } else {
                seq_stream.extend_from_slice(&record.sequence);
            }
        }

        if !no_quality {
            if let Some(ref qual) = record.quality {
                let quantized = quality::quantize_quality(qual, quality_mode);
                if const_qual_len == 0 {
                    write_varint(&mut qual_stream, quantized.len())?;
                }
                let packed = columnar::pack_qualities(&quantized, quality_binning);
                qual_stream.extend_from_slice(&packed);
            }
        }
    }

    Ok((header_stream, seq_stream, qual_stream, rc_flags))
}

/// Detect constant sequence/quality lengths from a set of records.
/// Returns (const_seq_len, const_qual_len) where 0 means variable lengths.
fn detect_const_lengths(records: &[crate::io::FastqRecord], no_quality: bool) -> (usize, usize) {
    if records.is_empty() {
        return (0, 0);
    }
    let first_seq_len = records[0].sequence.len();
    let const_seq_len = if records.iter().all(|r| r.sequence.len() == first_seq_len) {
        first_seq_len
    } else {
        0
    };
    let const_qual_len = if !no_quality {
        if let Some(ref q) = records[0].quality {
            let first_qual_len = q.len();
            if records.iter().all(|r| r.quality.as_ref().is_some_and(|q| q.len() == first_qual_len)) {
                first_qual_len
            } else {
                0
            }
        } else {
            0
        }
    } else {
        0
    };
    (const_seq_len, const_qual_len)
}

/// Sort records by reorder sort key, returning a new sorted Vec.
fn sort_records_by_key(records: Vec<crate::io::FastqRecord>) -> Vec<crate::io::FastqRecord> {
    use rayon::prelude::*;
    let keys: Vec<u128> = records
        .par_iter()
        .map(|r| {
            let qual_bytes = r.quality.as_ref().map(|q| q.as_slice()).unwrap_or(&[]);
            dna_utils::reorder_sort_key(&r.sequence, qual_bytes)
        })
        .collect();
    let mut keyed: Vec<(u128, crate::io::FastqRecord)> = keys.into_iter()
        .zip(records.into_iter())
        .collect();
    keyed.par_sort_unstable_by_key(|(k, _)| *k);
    keyed.into_iter().map(|(_, r)| r).collect()
}

/// Format bytes as human-readable string (e.g., "1.23 GiB").
fn humanize_bytes(bytes: usize) -> String {
    const KIB: f64 = 1024.0;
    const MIB: f64 = KIB * 1024.0;
    const GIB: f64 = MIB * 1024.0;
    let b = bytes as f64;
    if b >= GIB { format!("{:.2} GiB", b / GIB) }
    else if b >= MIB { format!("{:.1} MiB", b / MIB) }
    else if b >= KIB { format!("{:.0} KiB", b / KIB) }
    else { format!("{} B", bytes) }
}

/// Compress a single stream into BSC blocks using rayon par_iter.
/// Blocks are shrunk to actual compressed size to avoid excess capacity
/// (BSC allocates output buffers = input + header, but compressed data is typically ~40%).
fn compress_stream_to_bsc_blocks(data: &[u8], bsc_static: bool) -> Result<Vec<Vec<u8>>> {
    use rayon::prelude::*;
    const BSC_BLOCK_SIZE: usize = 25 * 1024 * 1024; // 25 MB: matches SPRING

    if data.is_empty() {
        return Ok(Vec::new());
    }

    let compress_fn: fn(&[u8]) -> Result<Vec<u8>> = if bsc_static {
        bsc::compress
    } else {
        bsc::compress_adaptive_no_lzp
    };

    let mut blocks: Vec<Vec<u8>> = data.par_chunks(BSC_BLOCK_SIZE)
        .map(|chunk| compress_fn(chunk))
        .collect::<Result<Vec<_>>>()?;

    // Shrink each block to actual compressed size. Without this, 1400
    // accumulated blocks x 25 MB excess capacity = ~35 GB wasted.
    for b in &mut blocks {
        b.shrink_to_fit();
    }

    Ok(blocks)
}

/// Write compressed blocks to a temp file in BSC multi-block body format:
/// [block_len: u32, block_data]... (the num_blocks header is written separately).
fn write_blocks_to_tmp(blocks: Vec<Vec<u8>>, tmp: &mut std::io::BufWriter<std::fs::File>) -> Result<u32> {
    use std::io::Write;
    let mut count = 0u32;
    for block in blocks {
        tmp.write_all(&(block.len() as u32).to_le_bytes())?;
        tmp.write_all(&block)?;
        count += 1;
    }
    Ok(count)
}

/// Write the chunked archive header (shared by all chunked compression paths).
/// Returns the metadata size in bytes (including v2 prefix).
fn write_archive_header(
    out: &mut impl std::io::Write,
    encoding_type: u8,
    quality_binning: QualityBinning,
    quality_compressor: QualityCompressor,
    no_quality: bool,
    num_reads: usize,
    headers_len: usize,
    sequences_len: usize,
    qualities_len: usize,
    const_seq_len: usize,
    const_qual_len: usize,
    rc_flags_len: usize,
) -> Result<usize> {
    let has_const_lengths = const_seq_len > 0 || const_qual_len > 0;
    let flags: u8 = if has_const_lengths { 0x02 } else { 0x00 };
    let const_len_size = if has_const_lengths { 8 } else { 0 };
    let body_size: usize = 9 + 2 + 1 + 40 + const_len_size;
    let header_size = V2_PREFIX_SIZE + body_size;

    // v2 prefix: magic + version + reserved + header_size
    out.write_all(&ARCHIVE_MAGIC)?;
    out.write_all(&[ARCHIVE_VERSION, 0])?;
    out.write_all(&(header_size as u32).to_le_bytes())?;

    // Header body (unchanged layout)
    out.write_all(&[encoding_type])?;
    out.write_all(&[flags])?;
    out.write_all(&[binning_to_code(quality_binning)])?;
    out.write_all(&[compressor_to_code(quality_compressor)])?;
    out.write_all(&[seq_compressor_to_code(SequenceCompressor::Bsc)])?;
    out.write_all(&[header_compressor_to_code(HeaderCompressor::Bsc)])?;
    out.write_all(&[0u8])?;                          // quality_model = disabled
    out.write_all(&[0u8])?;                          // quality_delta = disabled
    out.write_all(&[0u8])?;                          // quality_dict = disabled
    out.write_all(&0u16.to_le_bytes())?;             // template_prefix_len = 0
    out.write_all(&[0u8])?;                          // has_comment = false

    out.write_all(&(num_reads as u64).to_le_bytes())?;
    out.write_all(&(headers_len as u64).to_le_bytes())?;
    out.write_all(&(sequences_len as u64).to_le_bytes())?;
    out.write_all(&0u64.to_le_bytes())?;             // nmasks_len = 0
    out.write_all(&(if no_quality { 0 } else { qualities_len } as u64).to_le_bytes())?;

    if has_const_lengths {
        out.write_all(&(const_seq_len as u32).to_le_bytes())?;
        out.write_all(&(const_qual_len as u32).to_le_bytes())?;
    }

    let has_rc = rc_flags_len > 0;
    Ok(header_size + if has_rc { 8 } else { 0 })
}

/// Log compression stats.
fn log_compression_stats(
    original_size: usize,
    headers_len: usize,
    sequences_len: usize,
    qualities_len: usize,
    rc_flags_len: usize,
    metadata_size: usize,
    elapsed: std::time::Duration,
) {
    let total_compressed = headers_len + sequences_len + qualities_len + rc_flags_len + metadata_size;
    info!("Compression completed in {:.2}s", elapsed.as_secs_f64());
    info!("Original size: {} bytes", original_size);
    info!("Compressed size: {} bytes", total_compressed);
    info!("Compression ratio: {:.2}x", original_size as f64 / total_compressed as f64);
    info!("Stream breakdown:");
    info!("  Headers:   {} bytes ({:.1}%)", headers_len, 100.0 * headers_len as f64 / total_compressed as f64);
    info!("  Sequences: {} bytes ({:.1}%)", sequences_len, 100.0 * sequences_len as f64 / total_compressed as f64);
    info!("  Qualities: {} bytes ({:.1}%)", qualities_len, 100.0 * qualities_len as f64 / total_compressed as f64);
    if rc_flags_len > 0 {
        info!("  RC flags:  {} bytes ({:.1}%)", rc_flags_len, 100.0 * rc_flags_len as f64 / total_compressed as f64);
    }
    info!("  Metadata:  {} bytes ({:.1}%)", metadata_size, 100.0 * metadata_size as f64 / total_compressed as f64);
}

/// Optional RC (reverse-complement) flags stream appended after qualities.
struct RcStreamParams<'a> {
    flags_len: usize,
    num_blocks: u32,
    tmp_path: &'a std::path::Path,
}

/// Write the final archive from temp files (shared by reorder and global-reorder paths).
fn write_chunked_archive(
    output_path: &std::path::Path,
    quality_binning: QualityBinning,
    quality_compressor: QualityCompressor,
    no_quality: bool,
    encoding_type: u8,
    num_reads: usize,
    headers_len: usize,
    sequences_len: usize,
    qualities_len: usize,
    h_num_blocks: u32,
    s_num_blocks: u32,
    q_num_blocks: u32,
    h_tmp_path: &std::path::Path,
    s_tmp_path: &std::path::Path,
    q_tmp_path: &std::path::Path,
    original_size: usize,
    start_time: Instant,
    rc_stream: Option<RcStreamParams<'_>>,
    const_seq_len: usize,
    const_qual_len: usize,
) -> Result<()> {
    use std::io::{Write, BufReader, BufWriter};

    info!("Writing output file...");
    let mut output_file: Box<dyn Write> = if crate::cli::is_stdio_path(output_path) {
        Box::new(BufWriter::with_capacity(4 * 1024 * 1024, std::io::stdout().lock()))
    } else {
        Box::new(std::fs::File::create(output_path)?)
    };

    let rc_flags_len = rc_stream.as_ref().map(|rc| rc.flags_len).unwrap_or(0);
    let metadata_size = write_archive_header(
        &mut output_file, encoding_type, quality_binning, quality_compressor,
        no_quality, num_reads, headers_len, sequences_len, qualities_len,
        const_seq_len, const_qual_len, rc_flags_len,
    )?;

    let copy_stream = |num_blocks: u32, tmp_path: &std::path::Path, out: &mut dyn Write| -> Result<()> {
        if num_blocks == 0 { return Ok(()); }
        out.write_all(&num_blocks.to_le_bytes())?;
        let mut tmp_reader = BufReader::new(std::fs::File::open(tmp_path)?);
        std::io::copy(&mut tmp_reader, out)?;
        Ok(())
    };

    copy_stream(h_num_blocks, h_tmp_path, &mut *output_file)?;
    copy_stream(s_num_blocks, s_tmp_path, &mut *output_file)?;
    copy_stream(q_num_blocks, q_tmp_path, &mut *output_file)?;

    // Append RC flags stream after qualities (encoding_type=6 signals its presence)
    if let Some(ref rc) = rc_stream {
        if rc.num_blocks > 0 {
            output_file.write_all(&(rc.flags_len as u64).to_le_bytes())?;
            copy_stream(rc.num_blocks, rc.tmp_path, &mut *output_file)?;
        }
    }

    output_file.flush()?;

    log_compression_stats(original_size, headers_len, sequences_len, qualities_len, rc_flags_len, metadata_size, start_time.elapsed());

    Ok(())
}

/// Decompress headers using the appropriate compressor (BSC, Zstd, or OpenZL).
/// Handles template-encoded headers vs raw headers transparently.
fn decompress_headers_dispatch(
    header_compressor: HeaderCompressor,
    headers: &[u8],
    template_prefix_len: usize,
    read_id_template: &read_id::ReadIdTemplate,
    num_reads: usize,
) -> Result<Vec<Vec<u8>>> {
    match header_compressor {
        HeaderCompressor::Bsc => {
            if template_prefix_len > 0 {
                let header_data = bsc::decompress_parallel(headers)
                    .context("Failed to decompress headers (BSC)")?;
                let strings = read_id::decode_read_ids(&header_data, read_id_template, num_reads)
                    .context("Failed to decode template-encoded read IDs")?;
                Ok(strings.into_iter().map(String::into_bytes).collect())
            } else {
                codecs::decompress_headers_bsc(headers, num_reads)
                    .context("Failed to decompress headers (BSC)")
            }
        }
        HeaderCompressor::Zstd => {
            let header_data = decompress_zstd(headers)
                .context("Failed to decompress headers (zstd)")?;
            let strings = read_id::decode_read_ids(&header_data, read_id_template, num_reads)
                .context("Failed to decode read IDs")?;
            Ok(strings.into_iter().map(String::into_bytes).collect())
        }
        HeaderCompressor::OpenZl => {
            if template_prefix_len > 0 {
                let header_data = openzl::decompress_parallel(headers)
                    .context("Failed to decompress headers (OpenZL)")?;
                let strings = read_id::decode_read_ids(&header_data, read_id_template, num_reads)
                    .context("Failed to decode template-encoded read IDs")?;
                Ok(strings.into_iter().map(String::into_bytes).collect())
            } else {
                codecs::decompress_headers_openzl(headers, num_reads)
                    .context("Failed to decompress headers (OpenZL)")
            }
        }
    }
}

pub fn compress(args: &CompressConfig) -> Result<()> {
    compress_impl::compress(args)
}


/// Decompress a QZ archive back to FASTQ.
pub fn decompress(args: &DecompressConfig) -> Result<()> {
    decompress_impl::decompress(args)
}

/// Write variable-length integer
pub(crate) fn write_varint<W: std::io::Write>(writer: &mut W, mut value: usize) -> std::io::Result<()> {
    while value >= 0x80 {
        writer.write_all(&[((value & 0x7F) | 0x80) as u8])?;
        value >>= 7;
    }
    writer.write_all(&[value as u8])
}

/// Maximum varint size in bytes (10 bytes = up to 70 bits, enough for usize on 64-bit)
const MAX_VARINT_BYTES: usize = 10;

/// Read variable-length integer
fn read_varint(data: &[u8], offset: &mut usize) -> Option<usize> {
    let mut value = 0usize;
    let mut shift = 0;
    let mut bytes_read = 0;

    loop {
        if *offset >= data.len() {
            return None;
        }

        let byte = data[*offset];
        *offset += 1;
        bytes_read += 1;

        value |= ((byte & 0x7F) as usize) << shift;

        if byte & 0x80 == 0 {
            return Some(value);
        }

        shift += 7;
        if bytes_read >= MAX_VARINT_BYTES {
            return None; // Malformed varint: too many continuation bytes
        }
    }
}

/// Compress data with Zstd (alternative to gzip)
fn compress_zstd(data: &[u8], level: i32) -> Result<Vec<u8>> {
    zstd::bulk::compress(data, level)
        .map_err(|e| anyhow::anyhow!("Zstd compression failed: {}", e))
}

/// Decompress data with Zstd
fn decompress_zstd(compressed: &[u8]) -> Result<Vec<u8>> {
    zstd::bulk::decompress(compressed, 100_000_000) // 100MB limit
        .map_err(|e| anyhow::anyhow!("Zstd decompression failed: {}", e))
}
