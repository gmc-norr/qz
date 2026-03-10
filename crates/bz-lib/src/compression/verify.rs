use crate::cli::{VerifyConfig, VerifyResult};
use crate::compression::archive::{ArchiveHeader, ChunkHeader, CHUNK_FLAG_QUALITY_CTX};
use crate::compression::decompress::{
    decompress_streams, read_chunk_streams,
};
use crate::compression::streams::{self, ChunkConsensus};
use anyhow::{Context, Result};
use flate2::Crc;
use qz_lib::compression::bsc;
use std::io::BufReader;
use std::time::Instant;
use tracing::info;

/// Verify a BZ archive by decompressing all streams, validating their
/// structural integrity, and computing a CRC32 over the decompressed data.
///
/// Unlike the old approach that reconstructed every BAM record (nearly as
/// expensive as full decompression), this version:
/// - CRCs the decompressed stream bytes directly
/// - Validates stream structure (varints, lengths) without reconstructing records
/// - Skips quality_ctx decode entirely (quality_ctx compressed data is CRC'd as-is)
///
/// This makes verify significantly faster than decompress while still catching
/// any corruption in the compressed data or stream structure.
///
/// Note: the CRC32 produced here is over the decompressed stream data, NOT over
/// the reconstructed BAM output. It serves as an integrity fingerprint for the
/// archive contents.
pub fn verify(config: &VerifyConfig) -> Result<VerifyResult> {
    let start_time = Instant::now();

    if config.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build_global()
            .ok();
    }

    let file = std::fs::File::open(&config.input)?;
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);

    // Read and validate archive header
    let header = ArchiveHeader::read_from(&mut reader)?;
    info!(
        "Verifying archive: {} records in {} chunks (align_comp={}, aux_comp={})",
        header.num_records, header.num_chunks,
        header.alignment_compressor, header.aux_compressor,
    );

    // Decompress SAM header (validates it's not corrupt)
    let header_payload = bsc::decompress_parallel(&header.sam_header_compressed)
        .context("Failed to decompress SAM header")?;
    if header_payload.len() < 4 {
        anyhow::bail!("SAM header payload too short ({} bytes)", header_payload.len());
    }
    let header_raw_len =
        u32::from_le_bytes(header_payload[0..4].try_into().unwrap()) as usize;
    if 4 + header_raw_len > header_payload.len() {
        anyhow::bail!(
            "SAM header length {} exceeds payload size {}",
            header_raw_len, header_payload.len() - 4
        );
    }

    let mut crc = Crc::new();
    let mut total_bytes: u64 = 0;
    let mut total_verified: u64 = 0;

    // Hash the SAM header
    crc.update(&header_payload[4..4 + header_raw_len]);
    total_bytes += header_raw_len as u64;

    for chunk_idx in 0..header.num_chunks {
        let ch = ChunkHeader::read_from(&mut reader)?;
        let compressed_streams = read_chunk_streams(&mut reader, &ch)?;
        let uses_qctx = (ch.chunk_flags & CHUNK_FLAG_QUALITY_CTX) != 0;

        // Decompress all streams (this is the main integrity check —
        // BSC/zstd decompression will fail on corrupt data)
        let decompressed = decompress_streams(
            &compressed_streams, uses_qctx,
            header.alignment_compressor, header.aux_compressor,
        )?;

        let num_records = ch.num_records as usize;

        info!(
            "Verifying chunk {}/{}: {} records",
            chunk_idx + 1, header.num_chunks, num_records
        );

        // Validate consensus deserializes correctly
        let _consensus = ChunkConsensus::deserialize(&decompressed[0])
            .with_context(|| format!("chunk {}", chunk_idx))?;

        // Lightweight structural validation: walk all variable-length streams
        // to ensure varints and lengths are consistent
        validate_stream_structure(&decompressed, num_records, uses_qctx)
            .with_context(|| format!("chunk {} structural validation", chunk_idx))?;

        // CRC the decompressed streams directly
        for stream in &decompressed {
            crc.update(stream);
            total_bytes += stream.len() as u64;
        }

        // For quality_ctx, also CRC the compressed quality data
        // (since decompress_streams passes it through as empty)
        if uses_qctx {
            crc.update(&compressed_streams[13]);
            total_bytes += compressed_streams[13].len() as u64;
        }

        total_verified += num_records as u64;
    }

    let elapsed = start_time.elapsed().as_secs_f64();
    info!(
        "Verified {} records ({} bytes) in {:.2}s — CRC32: {:08x}",
        total_verified, total_bytes, elapsed, crc.sum()
    );

    Ok(VerifyResult {
        num_records: total_verified,
        num_chunks: header.num_chunks,
        alignment_compressor: header.alignment_compressor,
        aux_compressor: header.aux_compressor,
        crc32: crc.sum(),
        total_bytes,
        elapsed_secs: elapsed,
    })
}

/// Walk the variable-length streams to validate structural integrity without
/// reconstructing records. This catches corruption in varint encoding, stream
/// truncation, or inconsistent record counts.
fn validate_stream_structure(
    decompressed: &[Vec<u8>],
    num_records: usize,
    uses_quality_ctx: bool,
) -> Result<()> {
    // Fixed-width streams: just check total size
    // Stream 1 (ref_id): 4 bytes per record
    // Stream 2 (pos): 4 bytes per record
    // Stream 3 (mapq): 1 byte per record
    // Stream 4 (bin): 2 bytes per record
    // Stream 5 (flag): 2 bytes per record
    // Stream 6 (next_ref_id): 4 bytes per record
    // Stream 7 (next_pos): 4 bytes per record
    // Stream 8 (tlen): 4 bytes per record
    let expected_sizes: &[(usize, usize, &str)] = &[
        (1, 4, "ref_id"),
        (2, 4, "pos"),
        (3, 1, "mapq"),
        (4, 2, "bin"),
        (5, 2, "flag"),
        (6, 4, "next_ref_id"),
        (7, 4, "next_pos"),
        (8, 4, "tlen"),
    ];

    for &(stream_idx, bytes_per_record, name) in expected_sizes {
        let expected = num_records.checked_mul(bytes_per_record)
            .ok_or_else(|| anyhow::anyhow!(
                "stream {} ({}) size overflow: {} records × {} bytes",
                stream_idx, name, num_records, bytes_per_record
            ))?;
        let actual = decompressed[stream_idx].len();
        if actual != expected {
            anyhow::bail!(
                "stream {} ({}) size mismatch: expected {} bytes ({} records × {}), got {}",
                stream_idx, name, expected, num_records, bytes_per_record, actual
            );
        }
    }

    // Variable-length streams: walk varints to check they're well-formed
    // Stream 9: read_name (varint length + data)
    let mut off = 0;
    for rec_i in 0..num_records {
        let len = streams::read_varint(&decompressed[9], &mut off)
            .with_context(|| format!("read_name varint at record {rec_i}"))?;
        off = off.checked_add(len)
            .filter(|&e| e <= decompressed[9].len())
            .ok_or_else(|| anyhow::anyhow!("read_name overflow at record {rec_i}"))?;
    }
    if off != decompressed[9].len() {
        anyhow::bail!("stream 9 (read_name): {} trailing bytes", decompressed[9].len() - off);
    }

    // Stream 10: cigar (varint n_ops + n_ops*4 bytes)
    off = 0;
    for rec_i in 0..num_records {
        let n_ops = streams::read_varint(&decompressed[10], &mut off)
            .with_context(|| format!("cigar varint at record {rec_i}"))?;
        let byte_len = n_ops.checked_mul(4)
            .ok_or_else(|| anyhow::anyhow!("cigar overflow at record {rec_i}"))?;
        off = off.checked_add(byte_len)
            .filter(|&e| e <= decompressed[10].len())
            .ok_or_else(|| anyhow::anyhow!("cigar stream overflow at record {rec_i}"))?;
    }
    if off != decompressed[10].len() {
        anyhow::bail!("stream 10 (cigar): {} trailing bytes", decompressed[10].len() - off);
    }

    // Stream 11: seq_diff (varint count + packed nibbles)
    off = 0;
    for rec_i in 0..num_records {
        let count = streams::read_varint(&decompressed[11], &mut off)
            .with_context(|| format!("seq_diff varint at record {rec_i}"))?;
        let packed_len = (count + 1) / 2;
        off = off.checked_add(packed_len)
            .filter(|&e| e <= decompressed[11].len())
            .ok_or_else(|| anyhow::anyhow!("seq_diff overflow at record {rec_i}"))?;
    }
    if off != decompressed[11].len() {
        anyhow::bail!("stream 11 (seq_diff): {} trailing bytes", decompressed[11].len() - off);
    }

    // Stream 12: seq_extra (varint count + packed nibbles)
    off = 0;
    for rec_i in 0..num_records {
        let count = streams::read_varint(&decompressed[12], &mut off)
            .with_context(|| format!("seq_extra varint at record {rec_i}"))?;
        let packed_len = (count + 1) / 2;
        off = off.checked_add(packed_len)
            .filter(|&e| e <= decompressed[12].len())
            .ok_or_else(|| anyhow::anyhow!("seq_extra overflow at record {rec_i}"))?;
    }
    if off != decompressed[12].len() {
        anyhow::bail!("stream 12 (seq_extra): {} trailing bytes", decompressed[12].len() - off);
    }

    // Stream 13: quality (varint l_seq + l_seq bytes) — only if not quality_ctx
    if !uses_quality_ctx {
        off = 0;
        for rec_i in 0..num_records {
            let l_seq = streams::read_varint(&decompressed[13], &mut off)
                .with_context(|| format!("quality varint at record {rec_i}"))?;
            off = off.checked_add(l_seq)
                .filter(|&e| e <= decompressed[13].len())
                .ok_or_else(|| anyhow::anyhow!("quality overflow at record {rec_i}"))?;
        }
        if off != decompressed[13].len() {
            anyhow::bail!("stream 13 (quality): {} trailing bytes", decompressed[13].len() - off);
        }
    }

    // Stream 14: aux (varint length + data)
    off = 0;
    for rec_i in 0..num_records {
        let len = streams::read_varint(&decompressed[14], &mut off)
            .with_context(|| format!("aux varint at record {rec_i}"))?;
        off = off.checked_add(len)
            .filter(|&e| e <= decompressed[14].len())
            .ok_or_else(|| anyhow::anyhow!("aux overflow at record {rec_i}"))?;
    }
    if off != decompressed[14].len() {
        anyhow::bail!("stream 14 (aux): {} trailing bytes", decompressed[14].len() - off);
    }

    Ok(())
}
