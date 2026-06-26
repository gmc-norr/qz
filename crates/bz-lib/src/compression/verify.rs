use crate::cli::{VerifyConfig, VerifyResult};
use crate::compression::archive::{ArchiveHeader, CHUNK_FLAG_FQZ_QUALITY, ChunkHeader};
use crate::compression::decompress::{
    decompress_sam_header, decompress_streams, read_chunk_streams, validate_fqz_multiblock_framing,
};
use crate::compression::streams::{self, ChunkConsensus};
use anyhow::{Context, Result};
use flate2::Crc;
use std::io::BufReader;
use std::time::Instant;
use tracing::info;

/// Verify a BZ archive by decompressing all streams, validating their
/// structural integrity, and computing a CRC32 over the decompressed data.
///
/// This is a fast *integrity* check, not a full *fidelity* proof. It:
/// - CRCs the decompressed stream bytes directly
/// - Validates stream structure (varints, lengths) without reconstructing records
/// - Skips fqz decode entirely (fqz compressed data is CRC'd as-is)
///
/// This makes verify significantly faster than decompress while still catching
/// corruption in the compressed data or stream structure.
///
/// Note: the CRC32 produced here is over the decompressed stream data, NOT over
/// the reconstructed BAM output, so verify alone does not prove the archive
/// round-trips. The round-trip guarantee is enforced at *decompress* time: each
/// chunk stores a `content_crc32` over the original records, which `decompress`
/// recomputes over the reconstructed records and rejects on mismatch. A
/// future `--deep` verify could re-run that reconstruction here.
pub fn verify(config: &VerifyConfig) -> Result<VerifyResult> {
    let start_time = Instant::now();

    let file = std::fs::File::open(&config.input)?;
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);

    // Read and validate archive header
    let header = ArchiveHeader::read_from(&mut reader)?;
    info!(
        "Verifying archive: {} records in {} chunks",
        header.num_records, header.num_chunks,
    );

    // Decompress SAM header (validates it's not corrupt)
    let header_payload = decompress_sam_header(&header.sam_header_compressed)?;
    if header_payload.len() < 4 {
        anyhow::bail!(
            "SAM header payload too short ({} bytes)",
            header_payload.len()
        );
    }
    let header_raw_len = u32::from_le_bytes(header_payload[0..4].try_into().unwrap()) as usize;
    if 4 + header_raw_len > header_payload.len() {
        anyhow::bail!(
            "SAM header length {} exceeds payload size {}",
            header_raw_len,
            header_payload.len() - 4
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
        let uses_fqz_quality = (ch.chunk_flags & CHUNK_FLAG_FQZ_QUALITY) != 0;
        let num_records = ch.num_records as usize;

        // Decompress all streams (this is the main integrity check —
        // BSC decompression will fail on corrupt data). `num_records` bounds each
        // stream's decompressed output against memory-exhausting untrusted input.
        let decompressed = decompress_streams(&compressed_streams, uses_fqz_quality, num_records)?;

        info!(
            "Verifying chunk {}/{}: {} records",
            chunk_idx + 1,
            header.num_chunks,
            num_records
        );

        // Validate consensus deserializes correctly
        let _consensus = ChunkConsensus::deserialize(&decompressed[0])
            .with_context(|| format!("chunk {}", chunk_idx))?;

        // Lightweight structural validation: walk all variable-length streams
        // to ensure varints and lengths are consistent
        validate_stream_structure(&decompressed, num_records, uses_fqz_quality)
            .with_context(|| format!("chunk {} structural validation", chunk_idx))?;

        // CRC the decompressed streams directly
        for stream in &decompressed {
            crc.update(stream);
            total_bytes += stream.len() as u64;
        }

        // For fqz, stream 13 is the raw fqzcomp multiblock blob (passed
        // through empty by decompress_streams), so the structural walk above never
        // inspects it. Validate its framing and per-block CRC here — otherwise a
        // corrupted fqz payload (with the outer chunk CRC patched to match) would
        // pass verification yet fail to decompress. This checks the inner CRCs
        // without running the expensive fqz entropy decode.
        if uses_fqz_quality {
            validate_fqz_multiblock_framing(&compressed_streams[13], num_records)
                .with_context(|| format!("chunk {} fqz quality framing", chunk_idx))?;
            crc.update(&compressed_streams[13]);
            total_bytes += compressed_streams[13].len() as u64;
        }

        total_verified += num_records as u64;
    }

    // The per-chunk record counts must add up to the header's declared total.
    if total_verified != header.num_records {
        anyhow::bail!(
            "record count mismatch: archive header declares {} records but chunks total {} — archive is corrupted",
            header.num_records,
            total_verified
        );
    }

    let elapsed = start_time.elapsed().as_secs_f64();
    info!(
        "Verified {} records ({} bytes) in {:.2}s — CRC32: {:08x}",
        total_verified,
        total_bytes,
        elapsed,
        crc.sum()
    );

    Ok(VerifyResult {
        num_records: total_verified,
        num_chunks: header.num_chunks,
        crc32: crc.sum(),
        total_bytes,
        elapsed_secs: elapsed,
    })
}

/// Count set bits among the first `n` bits of a little-endian bit-packed bitmap
/// (`bit i` = `byte[i/8] >> (i%8) & 1`). Tolerates a short/truncated bitmap by
/// treating missing bits as unset, so a corrupt archive yields a size-mismatch
/// error downstream rather than an index panic here.
fn count_set_bits(bitmap: &[u8], n: usize) -> usize {
    (0..n)
        .filter(|&i| bitmap.get(i / 8).is_some_and(|b| (b >> (i % 8)) & 1 != 0))
        .count()
}

/// Walk the variable-length streams to validate structural integrity without
/// reconstructing records. This catches corruption in varint encoding, stream
/// truncation, or inconsistent record counts.
fn validate_stream_structure(
    decompressed: &[Vec<u8>],
    num_records: usize,
    uses_fqz_quality: bool,
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
        // streams 6 (next_ref_id) and 7 (next_pos) are checked below — they are
        // not fixed-width (see comment there). stream 8 (tlen) is likewise not
        // fixed-width: v9 stores only non-derivable tlen values.
    ];

    for &(stream_idx, bytes_per_record, name) in expected_sizes {
        let expected = num_records.checked_mul(bytes_per_record).ok_or_else(|| {
            anyhow::anyhow!(
                "stream {} ({}) size overflow: {} records × {} bytes",
                stream_idx,
                name,
                num_records,
                bytes_per_record
            )
        })?;
        let actual = decompressed[stream_idx].len();
        if actual != expected {
            anyhow::bail!(
                "stream {} ({}) size mismatch: expected {} bytes ({} records × {}), got {}",
                stream_idx,
                name,
                expected,
                num_records,
                bytes_per_record,
                actual
            );
        }
    }

    // Streams 6 (next_ref_id) and 7 (next_pos): v10 PNEXT/RNEXT derivation omits
    // BOTH for any record whose primary mate is in this chunk (a single bit in the
    // next_derivable bitmap, stream 21, governs both). So each holds 4 bytes per
    // NON-derivable record, not per record. Derive the expected size from the
    // bitmap's popcount instead of assuming fixed width.
    let derivable = count_set_bits(&decompressed[21], num_records);
    let non_derivable = num_records - derivable; // derivable <= num_records by construction
    let expected_mate = non_derivable.checked_mul(4).ok_or_else(|| {
        anyhow::anyhow!("mate stream size overflow: {non_derivable} records × 4 bytes")
    })?;
    for (stream_idx, name) in [(6usize, "next_ref_id"), (7usize, "next_pos")] {
        let actual = decompressed[stream_idx].len();
        if actual != expected_mate {
            anyhow::bail!(
                "stream {} ({}) size mismatch: expected {} bytes ({} non-derivable records × 4), got {}",
                stream_idx,
                name,
                expected_mate,
                non_derivable,
                actual
            );
        }
    }

    // Variable-length streams: walk varints to check they're well-formed
    // Stream 9: read_name (varint length + data)
    let mut off = 0;
    for rec_i in 0..num_records {
        let len = streams::read_varint(&decompressed[9], &mut off)
            .with_context(|| format!("read_name varint at record {rec_i}"))?;
        off = off
            .checked_add(len)
            .filter(|&e| e <= decompressed[9].len())
            .ok_or_else(|| anyhow::anyhow!("read_name overflow at record {rec_i}"))?;
    }
    if off != decompressed[9].len() {
        anyhow::bail!(
            "stream 9 (read_name): {} trailing bytes",
            decompressed[9].len() - off
        );
    }

    // Stream 10: cigar (varint n_ops + n_ops*4 bytes)
    off = 0;
    for rec_i in 0..num_records {
        let n_ops = streams::read_varint(&decompressed[10], &mut off)
            .with_context(|| format!("cigar varint at record {rec_i}"))?;
        let byte_len = n_ops
            .checked_mul(4)
            .ok_or_else(|| anyhow::anyhow!("cigar overflow at record {rec_i}"))?;
        off = off
            .checked_add(byte_len)
            .filter(|&e| e <= decompressed[10].len())
            .ok_or_else(|| anyhow::anyhow!("cigar stream overflow at record {rec_i}"))?;
    }
    if off != decompressed[10].len() {
        anyhow::bail!(
            "stream 10 (cigar): {} trailing bytes",
            decompressed[10].len() - off
        );
    }

    // Stream 11: seq_diff (varint count + packed nibbles)
    off = 0;
    for rec_i in 0..num_records {
        let count = streams::read_varint(&decompressed[11], &mut off)
            .with_context(|| format!("seq_diff varint at record {rec_i}"))?;
        let packed_len = count.div_ceil(2);
        off = off
            .checked_add(packed_len)
            .filter(|&e| e <= decompressed[11].len())
            .ok_or_else(|| anyhow::anyhow!("seq_diff overflow at record {rec_i}"))?;
    }
    if off != decompressed[11].len() {
        anyhow::bail!(
            "stream 11 (seq_diff): {} trailing bytes",
            decompressed[11].len() - off
        );
    }

    // Stream 12: seq_extra (varint count + packed nibbles)
    off = 0;
    for rec_i in 0..num_records {
        let count = streams::read_varint(&decompressed[12], &mut off)
            .with_context(|| format!("seq_extra varint at record {rec_i}"))?;
        let packed_len = count.div_ceil(2);
        off = off
            .checked_add(packed_len)
            .filter(|&e| e <= decompressed[12].len())
            .ok_or_else(|| anyhow::anyhow!("seq_extra overflow at record {rec_i}"))?;
    }
    if off != decompressed[12].len() {
        anyhow::bail!(
            "stream 12 (seq_extra): {} trailing bytes",
            decompressed[12].len() - off
        );
    }

    // Stream 13: quality (varint l_seq + l_seq bytes) — only if not fqz
    if !uses_fqz_quality {
        off = 0;
        for rec_i in 0..num_records {
            let l_seq = streams::read_varint(&decompressed[13], &mut off)
                .with_context(|| format!("quality varint at record {rec_i}"))?;
            off = off
                .checked_add(l_seq)
                .filter(|&e| e <= decompressed[13].len())
                .ok_or_else(|| anyhow::anyhow!("quality overflow at record {rec_i}"))?;
        }
        if off != decompressed[13].len() {
            anyhow::bail!(
                "stream 13 (quality): {} trailing bytes",
                decompressed[13].len() - off
            );
        }
    }

    // Stream 14: aux (varint length + data)
    off = 0;
    for rec_i in 0..num_records {
        let len = streams::read_varint(&decompressed[14], &mut off)
            .with_context(|| format!("aux varint at record {rec_i}"))?;
        off = off
            .checked_add(len)
            .filter(|&e| e <= decompressed[14].len())
            .ok_or_else(|| anyhow::anyhow!("aux overflow at record {rec_i}"))?;
    }
    if off != decompressed[14].len() {
        anyhow::bail!(
            "stream 14 (aux): {} trailing bytes",
            decompressed[14].len() - off
        );
    }

    Ok(())
}
