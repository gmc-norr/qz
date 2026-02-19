//! Stream codec helpers: compress/decompress functions for headers, sequences, and qualities.
//!
//! Each function handles a single stream (headers OR sequences OR qualities)
//! using a specific compressor (BSC, Zstd, OpenZL, fqzcomp, arithmetic).
//! Called by compress() and decompress() in the parent module.

use anyhow::Result;
use tracing::info;
use crate::cli::{QualityCompressor, SequenceCompressor, HeaderCompressor};
use super::*;

// ── ByteBackend: trait-based dispatch for BSC/Zstd/OpenZL ─────────────────

/// Backend byte-level compressor (BSC, Zstd, OpenZL).
///
/// Abstracts the `&[u8] -> Result<Vec<u8>>` compression/decompression
/// shared by BSC, Zstd, and OpenZL. Special-case codecs (Fqzcomp, QualityCtx)
/// operate on structured data and are excluded.
#[derive(Clone, Copy, Debug)]
pub(crate) enum ByteBackend {
    Bsc { use_static: bool },
    Zstd { level: i32 },
    OpenZl,
}

pub(crate) trait ByteCompressor {
    fn compress(&self, data: &[u8]) -> Result<Vec<u8>>;
    fn decompress(&self, data: &[u8]) -> Result<Vec<u8>>;
}

impl ByteCompressor for ByteBackend {
    fn compress(&self, data: &[u8]) -> Result<Vec<u8>> {
        match self {
            ByteBackend::Bsc { use_static } => bsc_compress_parallel(data, *use_static),
            ByteBackend::Zstd { level } => compress_zstd(data, *level),
            ByteBackend::OpenZl => openzl::compress_parallel(data),
        }
    }

    fn decompress(&self, data: &[u8]) -> Result<Vec<u8>> {
        match self {
            ByteBackend::Bsc { .. } => bsc::decompress_parallel(data),
            ByteBackend::Zstd { .. } => decompress_zstd(data),
            ByteBackend::OpenZl => openzl::decompress_parallel(data),
        }
    }
}

impl ByteBackend {
    /// Create from QualityCompressor, if it is a byte-level backend.
    /// Returns None for Fqzcomp and QualityCtx (they need structured data).
    #[allow(dead_code)]
    pub(crate) fn from_quality_compressor(
        qc: QualityCompressor,
        bsc_static: bool,
        level: i32,
    ) -> Option<Self> {
        match qc {
            QualityCompressor::Bsc => Some(ByteBackend::Bsc { use_static: bsc_static }),
            QualityCompressor::Zstd => Some(ByteBackend::Zstd { level }),
            QualityCompressor::OpenZl => Some(ByteBackend::OpenZl),
            QualityCompressor::Fqzcomp | QualityCompressor::QualityCtx => None,
        }
    }

    /// Create from SequenceCompressor.
    #[allow(dead_code)]
    pub(crate) fn from_sequence_compressor(
        sc: SequenceCompressor,
        bsc_static: bool,
        level: i32,
    ) -> Self {
        match sc {
            SequenceCompressor::Bsc => ByteBackend::Bsc { use_static: bsc_static },
            SequenceCompressor::Zstd => ByteBackend::Zstd { level },
            SequenceCompressor::OpenZl => ByteBackend::OpenZl,
        }
    }

    /// Create from HeaderCompressor.
    #[allow(dead_code)]
    pub(crate) fn from_header_compressor(
        hc: HeaderCompressor,
        bsc_static: bool,
        level: i32,
    ) -> Self {
        match hc {
            HeaderCompressor::Bsc | HeaderCompressor::Columnar => ByteBackend::Bsc { use_static: bsc_static },
            HeaderCompressor::Zstd => ByteBackend::Zstd { level },
            HeaderCompressor::OpenZl => ByteBackend::OpenZl,
        }
    }
}

// ── Shared stream-building helpers ────────────────────────────────────────

/// Build varint-prefixed header byte stream from records.
fn build_header_stream(records: &[crate::io::FastqRecord]) -> Result<Vec<u8>> {
    let mut stream = Vec::new();
    for record in records {
        write_varint(&mut stream, record.id.len())?;
        stream.extend_from_slice(&record.id);
    }
    Ok(stream)
}

/// Decode varint-prefixed header strings from decompressed data.
fn decode_header_stream(decompressed: &[u8], num_reads: usize) -> Result<Vec<Vec<u8>>> {
    let mut headers = Vec::with_capacity(num_reads);
    let mut offset = 0;
    for _ in 0..num_reads {
        let hdr_len = read_varint(decompressed, &mut offset)
            .ok_or_else(|| anyhow::anyhow!("Failed to read header length varint"))?;
        if offset + hdr_len > decompressed.len() {
            anyhow::bail!("Truncated header data at offset {}", offset);
        }
        headers.push(decompressed[offset..offset + hdr_len].to_vec());
        offset += hdr_len;
    }
    Ok(headers)
}

/// Build varint-prefixed packed quality byte stream from records.
fn build_quality_stream(records: &[crate::io::FastqRecord], binning: columnar::QualityBinning) -> Result<Vec<u8>> {
    use std::io::Write;
    let mut stream = Vec::new();
    for record in records {
        if let Some(qual) = &record.quality {
            write_varint(&mut stream, qual.len())?;
            let packed = columnar::pack_qualities(qual, binning);
            stream.write_all(&packed)?;
        }
    }
    Ok(stream)
}

/// Build varint-prefixed quality model delta byte stream from records.
fn build_quality_model_stream(records: &[crate::io::FastqRecord], model: &quality_model::QualityModel) -> Result<Vec<u8>> {
    use std::io::Write;
    let mut stream = Vec::new();
    for record in records {
        if let Some(qual) = &record.quality {
            write_varint(&mut stream, qual.len())?;
            let deltas = quality_model::encode_with_model(qual, model);
            let packed = quality_model::pack_deltas(&deltas);
            stream.write_all(&packed)?;
        }
    }
    Ok(stream)
}

// ── Headers ────────────────────────────────────────────────────────────────

/// Compress headers stream with template encoding + Zstd.
pub(crate) fn compress_headers(records: &[crate::io::FastqRecord], level: i32) -> Result<(Vec<u8>, read_id::ReadIdTemplate)> {
    let read_ids: Vec<String> = records.iter().map(|r| String::from_utf8_lossy(&r.id).into_owned()).collect();
    let encoded = read_id::compress_read_ids(&read_ids)?;
    let compressed = compress_zstd(&encoded.encoded_data, level)?;
    Ok((compressed, encoded.template))
}

/// Compress headers as raw ASCII + BSC (no template encoding).
/// Benchmark showed Raw + BSC (6.90x) dramatically beats Template + Zstd (3.64x).
pub(crate) fn compress_headers_bsc_with(records: &[crate::io::FastqRecord], bsc_static: bool) -> Result<Vec<u8>> {
    let header_stream = build_header_stream(records)?;
    bsc_compress_parallel(&header_stream, bsc_static)
}

/// Decompress raw BSC-compressed headers.
pub(super) fn decompress_headers_bsc(compressed: &[u8], num_reads: usize) -> Result<Vec<Vec<u8>>> {
    let decompressed = bsc::decompress_parallel(compressed)?;
    decode_header_stream(&decompressed, num_reads)
}

/// Compress headers with OpenZL.
pub(super) fn compress_headers_openzl(records: &[crate::io::FastqRecord]) -> Result<Vec<u8>> {
    let header_stream = build_header_stream(records)?;
    openzl::compress_parallel(&header_stream)
}

/// Decompress raw OpenZL-compressed headers.
pub(super) fn decompress_headers_openzl(compressed: &[u8], num_reads: usize) -> Result<Vec<Vec<u8>>> {
    let decompressed = openzl::decompress_parallel(compressed)?;
    decode_header_stream(&decompressed, num_reads)
}

/// Compress headers using columnar encoding (parse Illumina fields into typed columns + parallel BSC).
/// Falls back to raw BSC for non-Illumina headers.
pub(crate) fn compress_headers_columnar(records: &[crate::io::FastqRecord]) -> Result<Vec<u8>> {
    let header_strs: Vec<String> = records.iter()
        .map(|r| String::from_utf8_lossy(&r.id).into_owned())
        .collect();
    let header_refs: Vec<&str> = header_strs.iter().map(|s| s.as_str()).collect();
    header_col::compress_headers_columnar(&header_refs)
}

/// Decompress columnar-encoded headers.
pub(super) fn decompress_headers_columnar(compressed: &[u8], num_reads: usize) -> Result<Vec<Vec<u8>>> {
    header_col::decompress_headers_columnar(compressed, num_reads)
}

// ── Sequences ──────────────────────────────────────────────────────────────

/// Compress sequences as raw ASCII + BSC (no 2-bit encoding, no N-mask).
pub(crate) fn compress_sequences_raw_bsc_with(records: &[crate::io::FastqRecord], bsc_static: bool, sequence_hints: bool, sequence_delta: bool) -> Result<(Vec<u8>, Vec<u8>)> {
    let mut seq_stream = Vec::new();

    if sequence_delta {
        let mut delta_cache = std::collections::HashMap::new();
        let mut delta_reads = Vec::new();
        let mut delta_count = 0usize;
        for record in records {
            if encode_sequence_delta(&record.sequence, &mut seq_stream, &mut delta_cache, &mut delta_reads)? {
                delta_count += 1;
            }
        }
        info!("Inline delta: {}/{} reads delta-encoded ({:.1}%)",
            delta_count, records.len(), delta_count as f64 / records.len().max(1) as f64 * 100.0);
    } else {
        for record in records {
            write_varint(&mut seq_stream, record.sequence.len())?;
            if sequence_hints {
                seq_stream.push(dna_utils::compute_sequence_hint(&record.sequence));
            }
            seq_stream.extend_from_slice(&record.sequence);
        }
    }

    let compressed = bsc_compress_parallel(&seq_stream, bsc_static)?;
    Ok((compressed, Vec::new())) // empty nmasks
}

/// Decompress raw ASCII BSC-compressed sequences.
pub(super) fn decompress_sequences_raw_bsc(compressed: &[u8], num_reads: usize, encoding_type: u8, const_seq_len: usize) -> Result<Vec<Vec<u8>>> {
    let decompressed = bsc::decompress_parallel(compressed)?;
    let mut sequences = Vec::with_capacity(num_reads);
    let mut offset = 0;
    let skip_hints = encoding_type == 4;
    let decode_delta = encoding_type == 5;

    // Delta decoding needs a cache of all previously decoded sequences
    let mut delta_reads: Vec<Vec<u8>> = if decode_delta { Vec::with_capacity(num_reads) } else { Vec::new() };

    for _ in 0..num_reads {
        let seq_len = if const_seq_len > 0 {
            const_seq_len
        } else {
            read_varint(&decompressed, &mut offset)
                .ok_or_else(|| anyhow::anyhow!("Failed to read sequence length varint"))?
        };

        if decode_delta {
            // Read flag byte
            if offset >= decompressed.len() {
                anyhow::bail!("Truncated delta flag at offset {}", offset);
            }
            let flag = decompressed[offset];
            offset += 1;

            if flag == 0x01 {
                // Delta-encoded: read ref_idx varint, then delta bytes
                let ref_idx = read_varint(&decompressed, &mut offset)
                    .ok_or_else(|| anyhow::anyhow!("Failed to read delta ref_idx varint"))?;
                if ref_idx >= delta_reads.len() {
                    anyhow::bail!("Invalid delta ref_idx {} (only {} reads decoded so far)", ref_idx, delta_reads.len());
                }
                if offset + seq_len > decompressed.len() {
                    anyhow::bail!("Truncated delta sequence data at offset {}", offset);
                }
                let ref_seq = &delta_reads[ref_idx];
                let mut seq = Vec::with_capacity(seq_len);
                for i in 0..seq_len {
                    let byte = decompressed[offset + i];
                    if byte == 0x00 {
                        // Match: copy from reference
                        if i < ref_seq.len() {
                            seq.push(ref_seq[i]);
                        } else {
                            seq.push(b'N');
                        }
                    } else {
                        seq.push(byte); // Mismatch: actual base
                    }
                }
                offset += seq_len;
                delta_reads.push(seq.clone());
                sequences.push(seq);
            } else {
                // Raw: read seq_len bytes directly
                if offset + seq_len > decompressed.len() {
                    anyhow::bail!("Truncated sequence data at offset {}", offset);
                }
                let seq_bytes = decompressed[offset..offset + seq_len].to_vec();
                offset += seq_len;
                delta_reads.push(seq_bytes.clone());
                sequences.push(seq_bytes);
            }
        } else {
            if skip_hints {
                offset += 1; // skip hint byte
            }
            if offset + seq_len > decompressed.len() {
                anyhow::bail!("Truncated sequence data at offset {}", offset);
            }
            sequences.push(decompressed[offset..offset + seq_len].to_vec());
            offset += seq_len;
        }
    }

    Ok(sequences)
}

/// Compress sequences as raw ASCII + OpenZL.
pub(super) fn compress_sequences_raw_openzl(records: &[crate::io::FastqRecord]) -> Result<(Vec<u8>, Vec<u8>)> {
    let mut seq_stream = Vec::new();
    for record in records {
        write_varint(&mut seq_stream, record.sequence.len())?;
        seq_stream.extend_from_slice(&record.sequence);
    }

    let compressed = openzl::compress_parallel(&seq_stream)?;
    Ok((compressed, Vec::new())) // empty nmasks
}

/// Decompress raw ASCII OpenZL-compressed sequences.
pub(super) fn decompress_sequences_raw_openzl(compressed: &[u8], num_reads: usize, const_seq_len: usize) -> Result<Vec<Vec<u8>>> {
    let decompressed = openzl::decompress_parallel(compressed)?;
    let mut sequences = Vec::with_capacity(num_reads);
    let mut offset = 0;

    for _ in 0..num_reads {
        let seq_len = if const_seq_len > 0 {
            const_seq_len
        } else {
            read_varint(&decompressed, &mut offset)
                .ok_or_else(|| anyhow::anyhow!("Failed to read sequence length varint"))?
        };
        if offset + seq_len > decompressed.len() {
            anyhow::bail!("Truncated sequence data at offset {}", offset);
        }
        sequences.push(decompressed[offset..offset + seq_len].to_vec());
        offset += seq_len;
    }

    Ok(sequences)
}

/// Compress sequences as 2-bit encoded + N-mask + BSC.
pub(super) fn compress_sequences_2bit_bsc_with(records: &[crate::io::FastqRecord], bsc_static: bool) -> Result<(Vec<u8>, Vec<u8>)> {
    use std::io::Write;

    let mut seq_stream = Vec::new();
    let mut nmask_stream = Vec::new();
    for record in records {
        let encoding = n_mask::encode_with_n_mask(&record.sequence);
        write_varint(&mut seq_stream, encoding.length)?;
        seq_stream.write_all(&encoding.sequence_2bit)?;
        nmask_stream.write_all(&encoding.n_mask)?;
    }

    let (seq_result, nmask_result) = rayon::join(
        || bsc_compress_parallel(&seq_stream, bsc_static),
        || bsc_compress_parallel(&nmask_stream, bsc_static),
    );
    Ok((seq_result?, nmask_result?))
}

/// Decompress 2-bit BSC-compressed sequences with N-mask.
pub(super) fn decompress_sequences_2bit_bsc(sequences: &[u8], nmasks: &[u8], num_reads: usize, const_seq_len: usize) -> Result<Vec<Vec<u8>>> {
    let (seq_data, nmask_data) = rayon::join(
        || bsc::decompress_parallel(sequences),
        || bsc::decompress_parallel(nmasks),
    );
    let seq_data = seq_data?;
    let nmask_data = nmask_data?;

    let mut decoded = Vec::with_capacity(num_reads);
    let mut seq_offset = 0;
    let mut nmask_offset = 0;

    for _ in 0..num_reads {
        let seq_len = if const_seq_len > 0 {
            const_seq_len
        } else {
            read_varint(&seq_data, &mut seq_offset)
                .ok_or_else(|| anyhow::anyhow!("Failed to read sequence length varint"))?
        };
        let seq_2bit_len = (seq_len + 3) / 4;
        if seq_offset + seq_2bit_len > seq_data.len() {
            anyhow::bail!("Truncated 2-bit sequence data at offset {}", seq_offset);
        }
        let sequence_2bit = &seq_data[seq_offset..seq_offset + seq_2bit_len];
        seq_offset += seq_2bit_len;

        let nmask_len = (seq_len + 7) / 8;
        let n_mask = if nmask_offset + nmask_len <= nmask_data.len() {
            let mask = &nmask_data[nmask_offset..nmask_offset + nmask_len];
            nmask_offset += nmask_len;
            mask
        } else {
            &[]
        };

        let encoding = n_mask::NMaskEncoding {
            sequence_2bit: sequence_2bit.to_vec(),
            n_mask: n_mask.to_vec(),
            length: seq_len,
        };
        decoded.push(n_mask::decode_with_n_mask(&encoding));
    }

    Ok(decoded)
}

// ── Qualities ──────────────────────────────────────────────────────────────

/// Compress qualities stream (standard mode, dispatches by compressor).
pub(crate) fn compress_qualities_with(records: &[crate::io::FastqRecord], binning: columnar::QualityBinning, level: i32, compressor: QualityCompressor, bsc_static: bool) -> Result<Vec<u8>> {
    if compressor == QualityCompressor::Fqzcomp {
        // Wrap in multi-block format (1 block) for consistency with chunked path
        let blob = compress_qualities_fqzcomp(records)?;
        let mut out = Vec::with_capacity(8 + blob.len());
        out.extend_from_slice(&1u32.to_le_bytes()); // num_blocks = 1
        out.extend_from_slice(&(blob.len() as u32).to_le_bytes()); // block_len
        out.extend_from_slice(&blob);
        return Ok(out);
    }

    let quality_stream = build_quality_stream(records, binning)?;
    let backend = ByteBackend::from_quality_compressor(compressor, bsc_static, level)
        .expect("Fqzcomp handled above, QualityCtx handled separately");
    backend.compress(&quality_stream)
}

/// Compress qualities stream with dictionary.
pub(super) fn compress_qualities_with_dict(
    records: &[crate::io::FastqRecord],
    binning: columnar::QualityBinning,
    dictionary: &[u8],
    level: i32,
    compressor: QualityCompressor,
) -> Result<Vec<u8>> {
    let quality_stream = build_quality_stream(records, binning)?;

    match compressor {
        QualityCompressor::Zstd => zstd_dict::compress_with_dict(&quality_stream, dictionary, level),
        QualityCompressor::Bsc => {
            // BSC doesn't use dictionary - just use regular BSC compression
            info!("Note: BSC doesn't support dictionary mode, using standard BSC");
            bsc::compress_parallel_adaptive(&quality_stream) // dict path always adaptive
        }
        QualityCompressor::OpenZl => {
            info!("Note: OpenZL doesn't support dictionary mode, using standard OpenZL");
            openzl::compress_parallel(&quality_stream)
        }
        QualityCompressor::Fqzcomp => {
            info!("Note: fqzcomp doesn't support dictionary mode, using standard fqzcomp");
            compress_qualities_fqzcomp(records)
        }
        QualityCompressor::QualityCtx => unreachable!("quality_ctx handled separately"),
    }
}

/// Compress qualities stream with positional modeling.
pub(super) fn compress_qualities_with_model(
    records: &[crate::io::FastqRecord],
    level: i32,
    quality_compressor: QualityCompressor,
    bsc_static: bool,
) -> Result<(Vec<u8>, quality_model::QualityModel)> {
    let quality_strings: Vec<Vec<u8>> = records
        .iter()
        .filter_map(|r| r.quality.clone())
        .collect();

    if quality_strings.is_empty() {
        anyhow::bail!("No quality scores to model");
    }

    let model = quality_model::build_quality_model(&quality_strings)?;
    info!("Built quality model: {} positions, median quality {:.1}",
        model.positional_medians.len(),
        model.positional_medians.iter().map(|&q| q as f64).sum::<f64>() / model.positional_medians.len() as f64
    );

    let quality_stream = build_quality_model_stream(records, &model)?;

    let backend = match quality_compressor {
        QualityCompressor::Fqzcomp => {
            info!("Note: fqzcomp doesn't support quality modeling, using BSC for model deltas");
            ByteBackend::Bsc { use_static: bsc_static }
        }
        QualityCompressor::QualityCtx => unreachable!("quality_ctx handled separately"),
        _ => ByteBackend::from_quality_compressor(quality_compressor, bsc_static, level)
            .expect("Fqzcomp/QualityCtx handled above"),
    };
    let compressed = backend.compress(&quality_stream)?;

    Ok((compressed, model))
}

/// Compress qualities stream with delta encoding between adjacent reads.
pub(super) fn compress_qualities_with_delta(records: &[crate::io::FastqRecord], level: i32) -> Result<Vec<u8>> {
    use std::io::Write;

    let quality_strings: Vec<Vec<u8>> = records
        .iter()
        .filter_map(|r| r.quality.clone())
        .collect();

    if quality_strings.is_empty() {
        anyhow::bail!("No quality scores for delta encoding");
    }

    let encoded_deltas = quality_delta::encode_quality_deltas(&quality_strings);

    let stats = quality_delta::analyze_deltas(&encoded_deltas);
    info!("Quality delta encoding: {:.1}% zeros, {:.1}% small (|Δ|≤3), max delta: {}",
        stats.zero_percent, stats.small_percent, stats.max_delta);

    let mut quality_stream = Vec::new();
    for deltas in &encoded_deltas {
        write_varint(&mut quality_stream, deltas.len())?;
        let packed = quality_delta::pack_deltas(deltas);
        quality_stream.write_all(&packed)?;
    }

    compress_zstd(&quality_stream, level)
}

/// Compress qualities with model and dictionary.
pub(super) fn compress_qualities_with_model_and_dict(
    records: &[crate::io::FastqRecord],
    dictionary: &[u8],
    level: i32,
) -> Result<(Vec<u8>, quality_model::QualityModel)> {
    let quality_strings: Vec<Vec<u8>> = records
        .iter()
        .filter_map(|r| r.quality.clone())
        .collect();

    if quality_strings.is_empty() {
        anyhow::bail!("No quality scores to model");
    }

    let model = quality_model::build_quality_model(&quality_strings)?;
    info!("Built quality model: {} positions, median quality {:.1}",
        model.positional_medians.len(),
        model.positional_medians.iter().map(|&q| q as f64).sum::<f64>() / model.positional_medians.len() as f64
    );

    let quality_stream = build_quality_model_stream(records, &model)?;

    let compressed = zstd_dict::compress_with_dict(&quality_stream, dictionary, level)?;

    Ok((compressed, model))
}

/// Compress qualities with delta and dictionary.
pub(super) fn compress_qualities_with_delta_and_dict(
    records: &[crate::io::FastqRecord],
    dictionary: &[u8],
    level: i32,
) -> Result<Vec<u8>> {
    use std::io::Write;

    let quality_strings: Vec<Vec<u8>> = records
        .iter()
        .filter_map(|r| r.quality.clone())
        .collect();

    if quality_strings.is_empty() {
        anyhow::bail!("No quality scores for delta encoding");
    }

    let encoded_deltas = quality_delta::encode_quality_deltas(&quality_strings);

    let stats = quality_delta::analyze_deltas(&encoded_deltas);
    info!("Quality delta encoding: {:.1}% zeros, {:.1}% small (|Δ|≤3), max delta: {}",
        stats.zero_percent, stats.small_percent, stats.max_delta);

    let mut quality_stream = Vec::new();
    for deltas in &encoded_deltas {
        write_varint(&mut quality_stream, deltas.len())?;
        let packed = quality_delta::pack_deltas(deltas);
        quality_stream.write_all(&packed)?;
    }

    zstd_dict::compress_with_dict(&quality_stream, dictionary, level)
}

// ── Quality decompression ──────────────────────────────────────────────────

/// Dispatch quality decompression by compressor type.
pub(super) fn decompress_qualities_data(compressed: &[u8], compressor: QualityCompressor) -> Result<Vec<u8>> {
    match compressor {
        QualityCompressor::Fqzcomp => decompress_qualities_fqzcomp_multiblock(compressed),
        QualityCompressor::QualityCtx => unreachable!("quality_ctx decompression handled separately"),
        _ => {
            let backend = ByteBackend::from_quality_compressor(compressor, false, 0)
                .expect("Fqzcomp/QualityCtx handled above");
            backend.decompress(compressed)
        }
    }
}

/// Decompress multi-block fqzcomp quality data.
///
/// Format: [num_blocks: u32][block_len: u32, block_data]...
fn decompress_qualities_fqzcomp_multiblock(compressed: &[u8]) -> Result<Vec<u8>> {
    if compressed.len() < 4 {
        anyhow::bail!("fqzcomp multi-block: data too short");
    }

    let num_blocks = read_le_u32(compressed, 0)? as usize;

    if num_blocks == 0 {
        return Ok(Vec::new());
    }

    decompress_fqzcomp_blocks(compressed, num_blocks)
}

fn decompress_fqzcomp_blocks(compressed: &[u8], num_blocks: usize) -> Result<Vec<u8>> {
    let mut offset = 4; // skip num_blocks
    let mut output = Vec::new();

    for i in 0..num_blocks {
        if offset + 4 > compressed.len() {
            anyhow::bail!("fqzcomp multi-block: truncated block {} header", i);
        }
        let block_len = read_le_u32(compressed, offset)? as usize;
        offset += 4;
        if offset + block_len > compressed.len() {
            anyhow::bail!("fqzcomp multi-block: truncated block {} data ({} bytes, {} available)", i, block_len, compressed.len() - offset);
        }
        let block_data = &compressed[offset..offset + block_len];
        offset += block_len;

        let decompressed = decompress_qualities_fqzcomp(block_data)?;
        output.extend_from_slice(&decompressed);
    }

    Ok(output)
}

/// Compress quality scores using mean-quality reordering + fqzcomp.
///
/// Strategy: sort reads by mean quality score (1 byte key per read), then
/// compress reordered qualities with fqzcomp strat=0.
pub(super) fn compress_qualities_fqzcomp(records: &[crate::io::FastqRecord]) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    const FQZCOMP_SUB_CHUNK: usize = 500_000;

    let qual_strs: Vec<&[u8]> = records.iter()
        .filter_map(|r| r.quality.as_deref())
        .collect();

    if qual_strs.is_empty() {
        return Ok(Vec::new());
    }

    // Compute mean quality per read (1 byte each)
    let sort_keys: Vec<u8> = qual_strs.iter()
        .map(|q| {
            let sum: u64 = q.iter().map(|&b| b as u64).sum();
            (sum / q.len().max(1) as u64) as u8
        })
        .collect();

    // Sort by mean quality (stable sort)
    let mut indices: Vec<usize> = (0..qual_strs.len()).collect();
    indices.sort_by_key(|&i| sort_keys[i]);

    // Reorder quality strings
    let sorted_quals: Vec<&[u8]> = indices.iter().map(|&i| qual_strs[i]).collect();

    let num_sub_chunks = (sorted_quals.len() + FQZCOMP_SUB_CHUNK - 1) / FQZCOMP_SUB_CHUNK;
    info!("fqzcomp: {} reads, sorting by mean quality, {} sub-chunks, compressing with strat=0",
        qual_strs.len(), num_sub_chunks);

    // Compress sub-chunks in parallel with rayon, and BSC-compress sort keys concurrently
    let (sub_chunk_blobs, sort_keys_bsc) = rayon::join(
        || -> Result<Vec<Vec<u8>>> {
            sorted_quals.par_chunks(FQZCOMP_SUB_CHUNK)
                .map(|chunk| fqzcomp::compress(chunk, 0))
                .collect::<Result<Vec<_>>>()
        },
        || bsc::compress_parallel_adaptive(&sort_keys),
    );
    let sub_chunk_blobs = sub_chunk_blobs?;
    let sort_keys_bsc = sort_keys_bsc?;

    let fqzcomp_total: usize = sub_chunk_blobs.iter().map(|b| b.len()).sum();
    info!("fqzcomp: sort keys {} B (BSC'd), fqzcomp {} B in {} sub-chunks",
        sort_keys_bsc.len(), fqzcomp_total, sub_chunk_blobs.len());

    // Pack: [num_reads][sort_keys_bsc_len][sort_keys_bsc][num_sub_chunks][sub_chunk_len, data]...
    let mut output = Vec::with_capacity(12 + sort_keys_bsc.len() + fqzcomp_total + sub_chunk_blobs.len() * 4);
    output.extend_from_slice(&(qual_strs.len() as u32).to_le_bytes());
    output.extend_from_slice(&(sort_keys_bsc.len() as u32).to_le_bytes());
    output.extend_from_slice(&sort_keys_bsc);
    output.extend_from_slice(&(sub_chunk_blobs.len() as u32).to_le_bytes());
    for blob in &sub_chunk_blobs {
        output.extend_from_slice(&(blob.len() as u32).to_le_bytes());
        output.extend_from_slice(blob);
    }

    Ok(output)
}

/// Decompress fqzcomp quality data with mean-quality unsort.
///
/// Format: [num_reads: u32][sort_keys_bsc_len: u32][sort_keys_bsc]
///         [num_sub_chunks: u32][sub_chunk_len: u32, sub_chunk_data]...
fn decompress_qualities_fqzcomp(compressed: &[u8]) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    const FQZCOMP_SUB_CHUNK: usize = 500_000;

    if compressed.is_empty() {
        return Ok(Vec::new());
    }

    if compressed.len() < 12 {
        anyhow::bail!("fqzcomp quality data too short");
    }

    let mut offset = 0;

    // Parse header
    let num_reads = read_le_u32(compressed, offset)? as usize;
    offset += 4;
    let sort_keys_len = read_le_u32(compressed, offset)? as usize;
    offset += 4;

    info!("fqzcomp decompress: {} reads, sort_keys_bsc {} bytes, total blob {} bytes",
        num_reads, sort_keys_len, compressed.len());

    if offset + sort_keys_len > compressed.len() {
        anyhow::bail!("fqzcomp: truncated sort keys");
    }

    // BSC-decompress sort keys
    let sort_keys_data = &compressed[offset..offset + sort_keys_len];
    offset += sort_keys_len;

    // Parse sub-chunk table
    if offset + 4 > compressed.len() {
        anyhow::bail!("fqzcomp: truncated sub-chunk header");
    }
    let num_sub_chunks = read_le_u32(compressed, offset)? as usize;
    offset += 4;

    // Collect sub-chunk slices
    let mut sub_chunk_slices: Vec<&[u8]> = Vec::with_capacity(num_sub_chunks);
    for i in 0..num_sub_chunks {
        if offset + 4 > compressed.len() {
            anyhow::bail!("fqzcomp: truncated sub-chunk {} header", i);
        }
        let chunk_len = read_le_u32(compressed, offset)? as usize;
        offset += 4;
        if offset + chunk_len > compressed.len() {
            anyhow::bail!("fqzcomp: truncated sub-chunk {} data", i);
        }
        sub_chunk_slices.push(&compressed[offset..offset + chunk_len]);
        offset += chunk_len;
    }

    // Compute reads per sub-chunk (must match compression side)
    let reads_per_sub: Vec<usize> = (0..num_sub_chunks).map(|i| {
        let start = i * FQZCOMP_SUB_CHUNK;
        FQZCOMP_SUB_CHUNK.min(num_reads - start)
    }).collect();

    info!("fqzcomp decompress: {} sub-chunks, reads per sub: {:?}", num_sub_chunks, &reads_per_sub);

    // Decompress sort keys and fqzcomp sub-chunks in parallel
    let (sort_keys_result, sub_results): (Result<Vec<u8>>, Vec<Result<(Vec<u8>, Vec<i32>)>>) = rayon::join(
        || bsc::decompress_parallel(sort_keys_data),
        || sub_chunk_slices.par_iter().zip(reads_per_sub.par_iter())
            .map(|(data, &n_reads)| fqzcomp::decompress(data, n_reads))
            .collect(),
    );
    let sort_keys = sort_keys_result?;

    // Concatenate decoded qualities and lengths from all sub-chunks
    let mut decoded_quals = Vec::new();
    let mut lengths: Vec<i32> = Vec::with_capacity(num_reads);
    for res in sub_results {
        let (quals, lens) = res?;
        decoded_quals.extend_from_slice(&quals);
        lengths.extend_from_slice(&lens);
    }

    // Reconstruct permutation (same stable sort as compression)
    let mut indices: Vec<usize> = (0..num_reads).collect();
    indices.sort_by_key(|&i| sort_keys[i]);

    // Split decoded qualities into individual reads
    let mut sorted_quals: Vec<&[u8]> = Vec::with_capacity(num_reads);
    let mut q_offset = 0;
    for i in 0..num_reads {
        let len = lengths[i] as usize;
        if q_offset + len > decoded_quals.len() {
            anyhow::bail!("fqzcomp: truncated quality data at read {}", i);
        }
        sorted_quals.push(&decoded_quals[q_offset..q_offset + len]);
        q_offset += len;
    }

    // Inverse permute: unsorted[indices[i]] = sorted[i]
    let mut unsorted_quals: Vec<&[u8]> = vec![&[]; num_reads];
    for (sorted_idx, &original_idx) in indices.iter().enumerate() {
        unsorted_quals[original_idx] = sorted_quals[sorted_idx];
    }

    // Build varint-prefixed stream (raw ASCII, compatible with QualityBinning::None)
    let total_est: usize = unsorted_quals.iter().map(|q| q.len() + 2).sum();
    let mut output = Vec::with_capacity(total_est);
    for qual in &unsorted_quals {
        write_varint(&mut output, qual.len())?;
        output.extend_from_slice(qual);
    }

    Ok(output)
}
