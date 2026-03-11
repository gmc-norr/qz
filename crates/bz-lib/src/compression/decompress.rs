use crate::cli::DecompressConfig;
use crate::compression::archive::{ArchiveHeader, ChunkHeader, CHUNK_FLAG_QUALITY_CTX};
use crate::compression::streams::{self, ChunkConsensus, NUM_STREAMS};
use crate::io::bam::BamWriter;
use anyhow::{Context, Result};
use qz_lib::compression::bsc;
use qz_lib::compression::quality_ctx;
use flate2::Crc;
use std::io::{BufReader, Read, Write};
use std::num::NonZeroUsize;
use std::time::Instant;
use tracing::info;

/// Read a little-endian i32 from a stream at the given offset, advancing it by 4.
/// Returns an error if the stream is too short.
#[inline]
fn read_i32(stream: &[u8], off: &mut usize, name: &str) -> Result<i32> {
    let end = off.checked_add(4)
        .filter(|&e| e <= stream.len())
        .ok_or_else(|| anyhow::anyhow!("Truncated {name} stream at offset {}", *off))?;
    let val = i32::from_le_bytes(stream[*off..end].try_into().unwrap());
    *off = end;
    Ok(val)
}

/// Read a little-endian u16 from a stream at the given offset, advancing it by 2.
#[inline]
fn read_u16(stream: &[u8], off: &mut usize, name: &str) -> Result<u16> {
    let end = off.checked_add(2)
        .filter(|&e| e <= stream.len())
        .ok_or_else(|| anyhow::anyhow!("Truncated {name} stream at offset {}", *off))?;
    let val = u16::from_le_bytes(stream[*off..end].try_into().unwrap());
    *off = end;
    Ok(val)
}

/// Read a little-endian u32 from a stream at the given offset, advancing it by 4.
#[inline]
fn read_u32(stream: &[u8], off: &mut usize, name: &str) -> Result<u32> {
    let end = off.checked_add(4)
        .filter(|&e| e <= stream.len())
        .ok_or_else(|| anyhow::anyhow!("Truncated {name} stream at offset {}", *off))?;
    let val = u32::from_le_bytes(stream[*off..end].try_into().unwrap());
    *off = end;
    Ok(val)
}

/// Read a single byte from a stream at the given offset, advancing it by 1.
#[inline]
fn read_byte(stream: &[u8], off: &mut usize, name: &str) -> Result<u8> {
    if *off >= stream.len() {
        anyhow::bail!("Truncated {name} stream at offset {}", *off);
    }
    let val = stream[*off];
    *off += 1;
    Ok(val)
}

/// Get a sub-slice of a stream, advancing the offset.
#[inline]
fn read_slice<'a>(stream: &'a [u8], off: &mut usize, len: usize, name: &str) -> Result<&'a [u8]> {
    let end = off.checked_add(len)
        .filter(|&e| e <= stream.len())
        .ok_or_else(|| anyhow::anyhow!("Truncated {name} stream: need {} bytes at offset {}, have {}", len, *off, stream.len()))?;
    let slice = &stream[*off..end];
    *off = end;
    Ok(slice)
}

/// Read all compressed stream data for a chunk from the reader and verify CRC32.
///
/// The CRC32 is computed over the concatenated compressed stream payloads and
/// compared to `chunk_header.crc32`.  A mismatch means the archive has been
/// corrupted on disk; we bail before attempting (potentially slow) BSC decompression.
pub(super) fn read_chunk_streams<R: Read>(reader: &mut R, chunk_header: &ChunkHeader) -> Result<Vec<Vec<u8>>> {
    let mut compressed_streams: Vec<Vec<u8>> = Vec::with_capacity(NUM_STREAMS);
    let mut crc = Crc::new();
    for i in 0..NUM_STREAMS {
        let size = chunk_header.stream_sizes[i] as usize;
        let mut data = vec![0u8; size];
        if size > 0 {
            reader.read_exact(&mut data)?;
            crc.update(&data);
        }
        compressed_streams.push(data);
    }
    let computed = crc.sum();
    if computed != chunk_header.crc32 {
        anyhow::bail!(
            "Chunk CRC32 mismatch: stored {:08x}, computed {:08x} — archive is corrupted",
            chunk_header.crc32, computed
        );
    }
    Ok(compressed_streams)
}

/// Decompress a chunk using quality_ctx — parallel sequence reconstruction,
/// then quality_ctx decompression, then record assembly.
fn decompress_chunk_quality_ctx<W: Write>(
    decompressed: &[Vec<u8>],
    quality_compressed: &[u8],
    consensus: &ChunkConsensus,
    num_records: usize,
    bam_writer: &mut BamWriter<W>,
) -> Result<()> {
    use rayon::prelude::*;

    // Phase 1 (sequential): parse stream offsets and field values.
    let mut records: Vec<RecordInfo> = Vec::with_capacity(num_records);

    let mut ref_id_off = 0usize;
    let mut pos_off = 0usize;
    let mut mapq_off = 0usize;
    let mut bin_off = 0usize;
    let mut flag_off = 0usize;
    let mut next_ref_id_off = 0usize;
    let mut next_pos_off = 0usize;
    let mut tlen_off = 0usize;
    let mut read_name_off = 0usize;
    let mut cigar_off = 0usize;
    let mut seq_diff_off = 0usize;
    let mut seq_extra_off = 0usize;
    let mut aux_off = 0usize;

    let mut prev_ref_id: i32 = 0;
    let mut prev_pos: i32 = 0;

    for rec_i in 0..num_records {
        let delta_ref = read_i32(&decompressed[1], &mut ref_id_off, "ref_id")
            .with_context(|| format!("record {rec_i}"))?;
        let ref_id = prev_ref_id.wrapping_add(delta_ref);

        let delta_pos = read_i32(&decompressed[2], &mut pos_off, "pos")
            .with_context(|| format!("record {rec_i}"))?;
        let pos = if delta_ref == 0 {
            prev_pos.wrapping_add(delta_pos)
        } else {
            delta_pos
        };
        prev_ref_id = ref_id;
        prev_pos = pos;

        let mapq = read_byte(&decompressed[3], &mut mapq_off, "mapq")
            .with_context(|| format!("record {rec_i}"))?;

        let bin_val = read_u16(&decompressed[4], &mut bin_off, "bin")
            .with_context(|| format!("record {rec_i}"))?;

        let flag = read_u16(&decompressed[5], &mut flag_off, "flag")
            .with_context(|| format!("record {rec_i}"))?;

        let delta_next_ref = read_i32(&decompressed[6], &mut next_ref_id_off, "next_ref_id")
            .with_context(|| format!("record {rec_i}"))?;
        let next_ref_id = ref_id.wrapping_add(delta_next_ref);

        let delta_next_pos = read_i32(&decompressed[7], &mut next_pos_off, "next_pos")
            .with_context(|| format!("record {rec_i}"))?;
        let next_pos = pos.wrapping_add(delta_next_pos);

        let tlen_raw = read_u32(&decompressed[8], &mut tlen_off, "tlen")
            .with_context(|| format!("record {rec_i}"))?;
        let tlen = streams::zigzag_decode(tlen_raw);

        let rn_len = streams::read_varint(&decompressed[9], &mut read_name_off)?;
        let rn_start = read_name_off;
        read_name_off = read_name_off.checked_add(rn_len)
            .filter(|&e| e <= decompressed[9].len())
            .ok_or_else(|| anyhow::anyhow!("read_name overflow at record {rec_i}"))?;

        let n_cigar_op = streams::read_varint(&decompressed[10], &mut cigar_off)?;
        let cigar_start = cigar_off;
        let cigar_byte_len = n_cigar_op.checked_mul(4)
            .ok_or_else(|| anyhow::anyhow!("cigar overflow at record {rec_i}"))?;
        cigar_off = cigar_off.checked_add(cigar_byte_len)
            .filter(|&e| e <= decompressed[10].len())
            .ok_or_else(|| anyhow::anyhow!("cigar stream overflow at record {rec_i}"))?;

        let diff_count = streams::read_varint(&decompressed[11], &mut seq_diff_off)?;
        let diff_start = seq_diff_off;
        let diff_packed_len = diff_count.checked_add(1)
            .map(|n| n / 2)
            .ok_or_else(|| anyhow::anyhow!("diff_count overflow at record {rec_i}: {diff_count}"))?;
        seq_diff_off = seq_diff_off.checked_add(diff_packed_len)
            .filter(|&e| e <= decompressed[11].len())
            .ok_or_else(|| anyhow::anyhow!("seq_diff overflow at record {rec_i}"))?;

        let extra_count = streams::read_varint(&decompressed[12], &mut seq_extra_off)?;
        let extra_start = seq_extra_off;
        let extra_packed_len = extra_count.checked_add(1)
            .map(|n| n / 2)
            .ok_or_else(|| anyhow::anyhow!("extra_count overflow at record {rec_i}: {extra_count}"))?;
        seq_extra_off = seq_extra_off.checked_add(extra_packed_len)
            .filter(|&e| e <= decompressed[12].len())
            .ok_or_else(|| anyhow::anyhow!("seq_extra overflow at record {rec_i}"))?;

        let aux_len = streams::read_varint(&decompressed[14], &mut aux_off)?;
        let aux_start = aux_off;
        aux_off = aux_off.checked_add(aux_len)
            .filter(|&e| e <= decompressed[14].len())
            .ok_or_else(|| anyhow::anyhow!("aux stream overflow at record {rec_i}"))?;

        let l_seq = diff_count.checked_add(extra_count)
            .ok_or_else(|| anyhow::anyhow!("l_seq overflow at record {rec_i}: diff={diff_count} extra={extra_count}"))?;

        records.push(RecordInfo {
            ref_id, pos, mapq, bin_val, flag, next_ref_id, next_pos, tlen,
            rn_start, rn_len, cigar_start, n_cigar_op,
            diff_start, diff_packed_len, diff_count,
            extra_start, extra_packed_len, extra_count,
            aux_start, aux_len, l_seq,
        });
    }

    // Phase 2 (parallel): reconstruct sequences and convert to ASCII.
    let ascii_sequences: Vec<Vec<u8>> = records
        .par_iter()
        .map(|rec| {
            let diff_nibbles = streams::unpack_nibble_pairs(
                &decompressed[11][rec.diff_start..rec.diff_start + rec.diff_packed_len],
                rec.diff_count,
            );
            let extra_nibbles = streams::unpack_nibble_pairs(
                &decompressed[12][rec.extra_start..rec.extra_start + rec.extra_packed_len],
                rec.extra_count,
            );
            let seq_nibbles = streams::reconstruct_sequence_nibbles(
                consensus, rec.ref_id, rec.pos,
                &decompressed[10][rec.cigar_start..rec.cigar_start + rec.n_cigar_op * 4],
                rec.n_cigar_op, &diff_nibbles, &extra_nibbles, rec.l_seq,
            );
            streams::nibbles_to_ascii(&seq_nibbles)
        })
        .collect();

    // Phase 3: decompress quality with quality_ctx.
    let quality_arrays = quality_ctx::decompress_quality_ctx_multiblock(
        quality_compressed,
        &ascii_sequences,
    )?;

    // Phase 4 (sequential): assemble records and write.
    let mut data_buf = Vec::with_capacity(512);
    for (rec_idx, rec) in records.iter().enumerate() {
        let seq_packed = ascii_to_packed_seq(&ascii_sequences[rec_idx]);

        let read_name = &decompressed[9][rec.rn_start..rec.rn_start + rec.rn_len];
        let cigar_bytes = &decompressed[10][rec.cigar_start..rec.cigar_start + rec.n_cigar_op * 4];
        let aux_bytes = &decompressed[14][rec.aux_start..rec.aux_start + rec.aux_len];

        let record_len = 32 + read_name.len() + cigar_bytes.len() + seq_packed.len() + rec.l_seq + aux_bytes.len();
        data_buf.clear();
        data_buf.reserve(record_len);

        // BAM format limits: QNAME ≤ 254 bytes, n_cigar_op fits in u16, l_seq fits in i32
        if read_name.len() > 254 {
            anyhow::bail!("record {rec_idx}: read name too long ({} bytes, BAM max 254)", read_name.len());
        }
        if rec.n_cigar_op > u16::MAX as usize {
            anyhow::bail!("record {rec_idx}: n_cigar_op {} exceeds BAM u16 limit", rec.n_cigar_op);
        }
        if rec.l_seq > i32::MAX as usize {
            anyhow::bail!("record {rec_idx}: l_seq {} exceeds BAM i32 limit", rec.l_seq);
        }
        data_buf.extend_from_slice(&rec.ref_id.to_le_bytes());
        data_buf.extend_from_slice(&rec.pos.to_le_bytes());
        data_buf.push(read_name.len() as u8);
        data_buf.push(rec.mapq);
        data_buf.extend_from_slice(&rec.bin_val.to_le_bytes());
        data_buf.extend_from_slice(&(rec.n_cigar_op as u16).to_le_bytes());
        data_buf.extend_from_slice(&rec.flag.to_le_bytes());
        data_buf.extend_from_slice(&(rec.l_seq as i32).to_le_bytes());
        data_buf.extend_from_slice(&rec.next_ref_id.to_le_bytes());
        data_buf.extend_from_slice(&rec.next_pos.to_le_bytes());
        data_buf.extend_from_slice(&rec.tlen.to_le_bytes());
        data_buf.extend_from_slice(read_name);
        data_buf.extend_from_slice(cigar_bytes);
        data_buf.extend_from_slice(&seq_packed);

        // Quality: convert ASCII Phred+33 back to BAM Phred
        let qual = &quality_arrays[rec_idx];
        let qual_start = data_buf.len();
        data_buf.resize(qual_start + qual.len(), 0);
        streams::ascii_qual_to_bam(qual, &mut data_buf[qual_start..qual_start + qual.len()]);

        data_buf.extend_from_slice(aux_bytes);

        bam_writer.write_record(&data_buf)?;
    }

    Ok(())
}

/// Decompress a BZ archive back to a BAM file.
pub fn decompress(config: &DecompressConfig) -> Result<()> {
    let start_time = Instant::now();

    let file = std::fs::File::open(&config.input)?;
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);

    // Read archive header
    let header = ArchiveHeader::read_from(&mut reader)?;
    info!(
        "Archive: {} records in {} chunks (align_comp={}, aux_comp={})",
        header.num_records, header.num_chunks,
        header.alignment_compressor, header.aux_compressor,
    );

    // Decompress SAM header
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
    let header_raw = &header_payload[4..4 + header_raw_len];
    let ref_dict_bytes = &header_payload[4 + header_raw_len..];

    // Open BAM output with multi-threaded BGZF writer
    let bgzf_workers = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::new(1).unwrap())
        .min(NonZeroUsize::new(4).unwrap());
    let mut bam_writer = BamWriter::from_path_mt(&config.output, bgzf_workers)?;
    bam_writer.write_header(header_raw, ref_dict_bytes)?;

    let mut total_written: u64 = 0;

    for chunk_idx in 0..header.num_chunks {
        let chunk_header = ChunkHeader::read_from(&mut reader)?;
        let num_records = chunk_header.num_records as usize;

        info!(
            "Decompressing chunk {}/{}: {} records",
            chunk_idx + 1,
            header.num_chunks,
            num_records
        );

        let compressed_streams = read_chunk_streams(&mut reader, &chunk_header)?;
        let chunk_uses_quality_ctx = (chunk_header.chunk_flags & CHUNK_FLAG_QUALITY_CTX) != 0;

        let decompressed = decompress_streams(
            &compressed_streams,
            chunk_uses_quality_ctx,
            header.alignment_compressor,
            header.aux_compressor,
        )?;

        let consensus = ChunkConsensus::deserialize(&decompressed[0])
            .with_context(|| format!("chunk {}", chunk_idx))?;

        if chunk_uses_quality_ctx && num_records > 0 {
            decompress_chunk_quality_ctx(
                &decompressed,
                &compressed_streams[13],
                &consensus,
                num_records,
                &mut bam_writer,
            )?;
        } else {
            decompress_chunk_bsc(
                &decompressed, &consensus, num_records, &mut bam_writer,
            )?;
        }

        total_written += num_records as u64;
    }

    bam_writer.finish()?;

    let elapsed = start_time.elapsed().as_secs_f64();
    info!(
        "Decompressed {} records in {:.2}s",
        total_written, elapsed
    );

    Ok(())
}

/// Per-record intermediate data extracted in single pass.
struct RecordInfo {
    ref_id: i32,
    pos: i32,
    mapq: u8,
    bin_val: u16,
    flag: u16,
    next_ref_id: i32,
    next_pos: i32,
    tlen: i32,
    rn_start: usize,
    rn_len: usize,
    cigar_start: usize,
    n_cigar_op: usize,
    diff_start: usize,
    diff_packed_len: usize,
    diff_count: usize,
    extra_start: usize,
    extra_packed_len: usize,
    extra_count: usize,
    aux_start: usize,
    aux_len: usize,
    l_seq: usize,
}


/// Convert ASCII sequence (ACGTN) to BAM packed nibble format.
fn ascii_to_packed_seq(ascii: &[u8]) -> Vec<u8> {
    static ASCII_TO_NIBBLE: [u8; 256] = {
        let mut t = [15u8; 256]; // default N
        t[b'A' as usize] = 1;
        t[b'C' as usize] = 2;
        t[b'G' as usize] = 4;
        t[b'T' as usize] = 8;
        t[b'N' as usize] = 15;
        t[b'a' as usize] = 1;
        t[b'c' as usize] = 2;
        t[b'g' as usize] = 4;
        t[b't' as usize] = 8;
        t[b'n' as usize] = 15;
        t
    };

    let out_len = (ascii.len() + 1) / 2;
    let mut packed = vec![0u8; out_len];
    let pairs = ascii.len() / 2;

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("ssse3") {
            unsafe {
                ascii_to_packed_seq_ssse3(ascii, &mut packed, pairs, &ASCII_TO_NIBBLE);
            }
            return packed;
        }
    }

    for i in 0..pairs {
        let hi = ASCII_TO_NIBBLE[ascii[i * 2] as usize];
        let lo = ASCII_TO_NIBBLE[ascii[i * 2 + 1] as usize];
        packed[i] = (hi << 4) | lo;
    }
    if ascii.len() % 2 != 0 {
        packed[pairs] = ASCII_TO_NIBBLE[ascii[pairs * 2] as usize] << 4;
    }
    packed
}

/// SSSE3: hash ASCII bases to nibbles via pshufb, then pack pairs with maddubs.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
unsafe fn ascii_to_packed_seq_ssse3(
    ascii: &[u8], packed: &mut [u8], pairs: usize, scalar_lut: &[u8; 256],
) {
    use std::arch::x86_64::*;
    // Hash: (byte >> 1) & 0x0F maps A→0, C→1, G→3, T→10, N→7
    // LUT[idx] gives the BAM nibble value.
    static NIBBLE_LUT: [u8; 16] = [1, 2, 15, 4, 15, 15, 15, 15, 15, 15, 8, 15, 15, 15, 15, 15];
    let lut = unsafe { _mm_loadu_si128(NIBBLE_LUT.as_ptr() as *const __m128i) };
    let low4_mask = _mm_set1_epi8(0x0F);
    // Multiplier: even byte ×16, odd byte ×1 → packs nibble pair into one byte
    let pack_mult = _mm_set_epi8(1, 16, 1, 16, 1, 16, 1, 16, 1, 16, 1, 16, 1, 16, 1, 16);

    // Each SIMD iteration: 16 ASCII bytes → 8 packed bytes
    let simd_chunks = pairs / 8;
    for i in 0..simd_chunks {
        unsafe {
            let input = _mm_loadu_si128(ascii.as_ptr().add(i * 16) as *const __m128i);
            // (byte >> 1) & 0x0F — _mm_srli_epi16 shifts 16-bit lanes; the leaked
            // cross-byte bit is cleared by the & 0x0F mask.
            let hashed = _mm_and_si128(_mm_srli_epi16(input, 1), low4_mask);
            let nibbles = _mm_shuffle_epi8(lut, hashed);
            // Pack pairs: maddubs gives 16*even + odd in each 16-bit lane
            let packed16 = _mm_maddubs_epi16(nibbles, pack_mult);
            // Narrow 16-bit → 8-bit (values fit in [0,255])
            let packed8 = _mm_packus_epi16(packed16, packed16);
            // Store low 8 bytes
            std::ptr::copy_nonoverlapping(
                &packed8 as *const __m128i as *const u8,
                packed.as_mut_ptr().add(i * 8),
                8,
            );
        }
    }

    // Scalar tail
    let processed = simd_chunks * 8;
    for i in processed..pairs {
        let hi = scalar_lut[ascii[i * 2] as usize];
        let lo = scalar_lut[ascii[i * 2 + 1] as usize];
        packed[i] = (hi << 4) | lo;
    }
    // Odd trailing base
    if ascii.len() % 2 != 0 {
        packed[pairs] = scalar_lut[ascii[pairs * 2] as usize] << 4;
    }
}

/// Decompress a chunk using BSC-only quality path (no quality_ctx).
fn decompress_chunk_bsc<W: Write>(
    decompressed: &[Vec<u8>],
    consensus: &ChunkConsensus,
    num_records: usize,
    bam_writer: &mut BamWriter<W>,
) -> Result<()> {
    let mut ref_id_off = 0usize;
    let mut pos_off = 0usize;
    let mut mapq_off = 0usize;
    let mut bin_off = 0usize;
    let mut flag_off = 0usize;
    let mut next_ref_id_off = 0usize;
    let mut next_pos_off = 0usize;
    let mut tlen_off = 0usize;
    let mut read_name_off = 0usize;
    let mut cigar_off = 0usize;
    let mut seq_diff_off = 0usize;
    let mut seq_extra_off = 0usize;
    let mut qual_off = 0usize;
    let mut aux_off = 0usize;

    let mut prev_ref_id: i32 = 0;
    let mut prev_pos: i32 = 0;

    for rec_i in 0..num_records {
        let delta_ref = read_i32(&decompressed[1], &mut ref_id_off, "ref_id")
            .with_context(|| format!("record {rec_i}"))?;
        let ref_id = prev_ref_id.wrapping_add(delta_ref);

        let delta_pos = read_i32(&decompressed[2], &mut pos_off, "pos")
            .with_context(|| format!("record {rec_i}"))?;
        let pos = if delta_ref == 0 {
            prev_pos.wrapping_add(delta_pos)
        } else {
            delta_pos
        };
        prev_ref_id = ref_id;
        prev_pos = pos;

        let mapq = read_byte(&decompressed[3], &mut mapq_off, "mapq")
            .with_context(|| format!("record {rec_i}"))?;

        let bin_val = read_u16(&decompressed[4], &mut bin_off, "bin")
            .with_context(|| format!("record {rec_i}"))?;

        let flag = read_u16(&decompressed[5], &mut flag_off, "flag")
            .with_context(|| format!("record {rec_i}"))?;

        let delta_next_ref = read_i32(&decompressed[6], &mut next_ref_id_off, "next_ref_id")
            .with_context(|| format!("record {rec_i}"))?;
        let next_ref_id_val = ref_id.wrapping_add(delta_next_ref);

        let delta_next_pos = read_i32(&decompressed[7], &mut next_pos_off, "next_pos")
            .with_context(|| format!("record {rec_i}"))?;
        let next_pos_val = pos.wrapping_add(delta_next_pos);

        let tlen_raw = read_u32(&decompressed[8], &mut tlen_off, "tlen")
            .with_context(|| format!("record {rec_i}"))?;
        let tlen_val = streams::zigzag_decode(tlen_raw);

        let rn_len = streams::read_varint(&decompressed[9], &mut read_name_off)?;
        if rn_len > 254 {
            anyhow::bail!("record {rec_i}: read name length {rn_len} exceeds BAM limit of 254");
        }
        let read_name = read_slice(&decompressed[9], &mut read_name_off, rn_len, "read_name")
            .with_context(|| format!("record {rec_i}"))?;

        let n_cigar_op = streams::read_varint(&decompressed[10], &mut cigar_off)?;
        let cigar_byte_len = n_cigar_op.checked_mul(4)
            .ok_or_else(|| anyhow::anyhow!("cigar overflow at record {rec_i}"))?;
        let cigar_bytes = read_slice(&decompressed[10], &mut cigar_off, cigar_byte_len, "cigar")
            .with_context(|| format!("record {rec_i}"))?;

        let diff_count = streams::read_varint(&decompressed[11], &mut seq_diff_off)?;
        let diff_packed_len = diff_count.checked_add(1)
            .map(|n| n / 2)
            .ok_or_else(|| anyhow::anyhow!("diff_count overflow at record {rec_i}: {diff_count}"))?;
        let diff_data = read_slice(&decompressed[11], &mut seq_diff_off, diff_packed_len, "seq_diff")
            .with_context(|| format!("record {rec_i}"))?;
        let diff_nibbles = streams::unpack_nibble_pairs(diff_data, diff_count);

        let extra_count = streams::read_varint(&decompressed[12], &mut seq_extra_off)?;
        let extra_packed_len = extra_count.checked_add(1)
            .map(|n| n / 2)
            .ok_or_else(|| anyhow::anyhow!("extra_count overflow at record {rec_i}: {extra_count}"))?;
        let extra_data = read_slice(&decompressed[12], &mut seq_extra_off, extra_packed_len, "seq_extra")
            .with_context(|| format!("record {rec_i}"))?;
        let extra_nibbles = streams::unpack_nibble_pairs(extra_data, extra_count);

        let l_seq = streams::read_varint(&decompressed[13], &mut qual_off)?;
        let qual_bytes = read_slice(&decompressed[13], &mut qual_off, l_seq, "quality")
            .with_context(|| format!("record {rec_i}"))?;

        let aux_len = streams::read_varint(&decompressed[14], &mut aux_off)?;
        let aux_bytes = read_slice(&decompressed[14], &mut aux_off, aux_len, "aux")
            .with_context(|| format!("record {rec_i}"))?;

        let record_data = streams::reconstruct_record(
            consensus, ref_id, pos, mapq, bin_val, flag,
            next_ref_id_val, next_pos_val, tlen_val,
            read_name, cigar_bytes, n_cigar_op as u16,
            &diff_nibbles, &extra_nibbles, qual_bytes, aux_bytes, l_seq,
        );

        bam_writer.write_record(&record_data)?;
    }

    Ok(())
}

/// Decompress a single stream using the appropriate compressor.
fn decompress_one(data: &[u8], use_zstd: bool) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }
    if use_zstd {
        // Use streaming decode instead of bulk::decompress to avoid hardcoded output size cap.
        // zstd frames contain the decompressed size in the header; decode_all handles this.
        Ok(zstd::stream::decode_all(data)?)
    } else {
        bsc::decompress_parallel(data)
    }
}

/// Decompress all 15 streams in parallel groups, dispatching by compressor type.
pub(super) fn decompress_streams(
    compressed: &[Vec<u8>],
    use_quality_ctx: bool,
    alignment_compressor: u8,
    aux_compressor: u8,
) -> Result<Vec<Vec<u8>>> {
    let align_zstd = alignment_compressor == 1;
    let aux_zstd = aux_compressor == 1;

    let ((group_a, group_b), (group_c, group_d)) = rayon::join(
        || {
            rayon::join(
                || -> Result<Vec<Vec<u8>>> {
                    // consensus (index 0) always BSC, alignment (1-5) depends on config
                    let mut results = Vec::with_capacity(6);
                    results.push(decompress_one(&compressed[0], false)?);
                    for data in &compressed[1..6] {
                        results.push(decompress_one(data, align_zstd)?);
                    }
                    Ok(results)
                },
                || -> Result<Vec<Vec<u8>>> {
                    // alignment streams 6-8
                    let mut results = Vec::with_capacity(3);
                    for data in &compressed[6..9] {
                        results.push(decompress_one(data, align_zstd)?);
                    }
                    Ok(results)
                },
            )
        },
        || {
            rayon::join(
                || -> Result<Vec<Vec<u8>>> {
                    // read data streams 9-12: always BSC
                    let mut results = Vec::with_capacity(4);
                    for data in &compressed[9..13] {
                        results.push(decompress_one(data, false)?);
                    }
                    Ok(results)
                },
                || -> Result<Vec<Vec<u8>>> {
                    if use_quality_ctx {
                        // quality_ctx: pass through raw, decompress aux
                        let qual_passthrough = Vec::new();
                        let aux = decompress_one(&compressed[14], aux_zstd)?;
                        Ok(vec![qual_passthrough, aux])
                    } else {
                        // quality (13) always BSC, aux (14) depends on config
                        let qual = decompress_one(&compressed[13], false)?;
                        let aux = decompress_one(&compressed[14], aux_zstd)?;
                        Ok(vec![qual, aux])
                    }
                },
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
