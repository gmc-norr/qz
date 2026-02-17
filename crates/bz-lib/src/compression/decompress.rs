use crate::cli::DecompressConfig;
use crate::compression::archive::{ArchiveHeader, ChunkHeader, FLAG_QUALITY_CTX};
use crate::compression::streams::{self, ChunkConsensus, NUM_STREAMS};
use crate::io::bam::BamWriter;
use anyhow::Result;
use qz_lib::compression::bsc;
use qz_lib::compression::quality_ctx;
use std::io::{BufReader, Read, Write};
use std::num::NonZeroUsize;
use std::time::Instant;
use tracing::info;

/// Decompress a BZ archive back to a BAM file.
pub fn decompress(config: &DecompressConfig) -> Result<()> {
    let start_time = Instant::now();

    if config.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build_global()
            .ok();
    }

    let file = std::fs::File::open(&config.input)?;
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);

    // Read archive header
    let header = ArchiveHeader::read_from(&mut reader)?;
    let use_quality_ctx = (header.flags & FLAG_QUALITY_CTX) != 0;
    info!(
        "Archive: {} records in {} chunks (quality_ctx: {})",
        header.num_records, header.num_chunks, use_quality_ctx
    );

    // Decompress SAM header
    let header_payload = bsc::decompress_parallel(&header.sam_header_compressed)?;
    let header_raw_len =
        u32::from_le_bytes(header_payload[0..4].try_into().unwrap()) as usize;
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

        // Read all compressed stream data for this chunk
        let mut compressed_streams: Vec<Vec<u8>> = Vec::with_capacity(NUM_STREAMS);
        for i in 0..NUM_STREAMS {
            let size = chunk_header.stream_sizes[i] as usize;
            let mut data = vec![0u8; size];
            if size > 0 {
                reader.read_exact(&mut data)?;
            }
            compressed_streams.push(data);
        }

        // Decompress non-quality streams with BSC, quality handled separately
        let decompressed = decompress_streams(&compressed_streams, use_quality_ctx)?;

        // Parse consensus
        let consensus = ChunkConsensus::deserialize(&decompressed[0]);

        if use_quality_ctx && num_records > 0 {
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
    // This is fast cursor arithmetic — no heavy compute.
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

    for _ in 0..num_records {
        let delta_ref = i32::from_le_bytes(
            decompressed[1][ref_id_off..ref_id_off + 4].try_into().unwrap(),
        );
        ref_id_off += 4;
        let ref_id = prev_ref_id.wrapping_add(delta_ref);

        let delta_pos = i32::from_le_bytes(
            decompressed[2][pos_off..pos_off + 4].try_into().unwrap(),
        );
        pos_off += 4;
        let pos = if delta_ref == 0 {
            prev_pos.wrapping_add(delta_pos)
        } else {
            delta_pos
        };
        prev_ref_id = ref_id;
        prev_pos = pos;

        let mapq = decompressed[3][mapq_off];
        mapq_off += 1;

        let bin_val = u16::from_le_bytes(
            decompressed[4][bin_off..bin_off + 2].try_into().unwrap(),
        );
        bin_off += 2;

        let flag = u16::from_le_bytes(
            decompressed[5][flag_off..flag_off + 2].try_into().unwrap(),
        );
        flag_off += 2;

        let delta_next_ref = i32::from_le_bytes(
            decompressed[6][next_ref_id_off..next_ref_id_off + 4].try_into().unwrap(),
        );
        next_ref_id_off += 4;
        let next_ref_id = ref_id.wrapping_add(delta_next_ref);

        let delta_next_pos = i32::from_le_bytes(
            decompressed[7][next_pos_off..next_pos_off + 4].try_into().unwrap(),
        );
        next_pos_off += 4;
        let next_pos = pos.wrapping_add(delta_next_pos);

        let tlen_raw = u32::from_le_bytes(
            decompressed[8][tlen_off..tlen_off + 4].try_into().unwrap(),
        );
        tlen_off += 4;
        let tlen = streams::zigzag_decode(tlen_raw);

        let rn_len = streams::read_varint(&decompressed[9], &mut read_name_off);
        let rn_start = read_name_off;
        read_name_off += rn_len;

        let n_cigar_op = streams::read_varint(&decompressed[10], &mut cigar_off);
        let cigar_start = cigar_off;
        cigar_off += n_cigar_op * 4;

        let diff_count = streams::read_varint(&decompressed[11], &mut seq_diff_off);
        let diff_start = seq_diff_off;
        let diff_packed_len = (diff_count + 1) / 2;
        seq_diff_off += diff_packed_len;

        let extra_count = streams::read_varint(&decompressed[12], &mut seq_extra_off);
        let extra_start = seq_extra_off;
        let extra_packed_len = (extra_count + 1) / 2;
        seq_extra_off += extra_packed_len;

        let aux_len = streams::read_varint(&decompressed[14], &mut aux_off);
        let aux_start = aux_off;
        aux_off += aux_len;

        let l_seq = diff_count + extra_count;

        records.push(RecordInfo {
            ref_id, pos, mapq, bin_val, flag, next_ref_id, next_pos, tlen,
            rn_start, rn_len, cigar_start, n_cigar_op,
            diff_start, diff_packed_len, diff_count,
            extra_start, extra_packed_len, extra_count,
            aux_start, aux_len, l_seq,
        });
    }

    // Phase 2 (parallel): reconstruct sequences and convert to ASCII.
    // Each record is independent — consensus/decompressed are read-only.
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
            seq_nibbles.iter().map(|&n| streams::nibble_to_ascii(n)).collect()
        })
        .collect();

    // Phase 3: decompress quality with quality_ctx.
    let quality_arrays = quality_ctx::decompress_quality_ctx_multiblock(
        quality_compressed,
        &ascii_sequences,
    )?;

    // Phase 4 (sequential): assemble records and write.
    // Multi-threaded BGZF writer handles deflate parallelism.
    let mut data_buf = Vec::with_capacity(512);
    for (rec_idx, rec) in records.iter().enumerate() {
        let seq_packed = ascii_to_packed_seq(&ascii_sequences[rec_idx]);

        let read_name = &decompressed[9][rec.rn_start..rec.rn_start + rec.rn_len];
        let cigar_bytes = &decompressed[10][rec.cigar_start..rec.cigar_start + rec.n_cigar_op * 4];
        let aux_bytes = &decompressed[14][rec.aux_start..rec.aux_start + rec.aux_len];

        let record_len = 32 + read_name.len() + cigar_bytes.len() + seq_packed.len() + rec.l_seq + aux_bytes.len();
        data_buf.clear();
        data_buf.reserve(record_len);

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

        // Quality: convert ASCII Phred+33 back to BAM Phred inline
        let qual = &quality_arrays[rec_idx];
        for &q in qual.iter() {
            data_buf.push(q.saturating_sub(33));
        }

        data_buf.extend_from_slice(aux_bytes);

        bam_writer.write_record(&data_buf)?;
    }

    Ok(())
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

    let mut packed = Vec::with_capacity((ascii.len() + 1) / 2);
    for chunk in ascii.chunks(2) {
        let hi = ASCII_TO_NIBBLE[chunk[0] as usize];
        let lo = if chunk.len() > 1 { ASCII_TO_NIBBLE[chunk[1] as usize] } else { 0 };
        packed.push((hi << 4) | lo);
    }
    packed
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

    for _ in 0..num_records {
        let delta_ref = i32::from_le_bytes(
            decompressed[1][ref_id_off..ref_id_off + 4].try_into().unwrap(),
        );
        ref_id_off += 4;
        let ref_id = prev_ref_id.wrapping_add(delta_ref);

        let delta_pos = i32::from_le_bytes(
            decompressed[2][pos_off..pos_off + 4].try_into().unwrap(),
        );
        pos_off += 4;
        let pos = if delta_ref == 0 {
            prev_pos.wrapping_add(delta_pos)
        } else {
            delta_pos
        };
        prev_ref_id = ref_id;
        prev_pos = pos;

        let mapq = decompressed[3][mapq_off];
        mapq_off += 1;

        let bin_val = u16::from_le_bytes(
            decompressed[4][bin_off..bin_off + 2].try_into().unwrap(),
        );
        bin_off += 2;

        let flag = u16::from_le_bytes(
            decompressed[5][flag_off..flag_off + 2].try_into().unwrap(),
        );
        flag_off += 2;

        let delta_next_ref = i32::from_le_bytes(
            decompressed[6][next_ref_id_off..next_ref_id_off + 4].try_into().unwrap(),
        );
        next_ref_id_off += 4;
        let next_ref_id_val = ref_id.wrapping_add(delta_next_ref);

        let delta_next_pos = i32::from_le_bytes(
            decompressed[7][next_pos_off..next_pos_off + 4].try_into().unwrap(),
        );
        next_pos_off += 4;
        let next_pos_val = pos.wrapping_add(delta_next_pos);

        let tlen_raw = u32::from_le_bytes(
            decompressed[8][tlen_off..tlen_off + 4].try_into().unwrap(),
        );
        tlen_off += 4;
        let tlen_val = streams::zigzag_decode(tlen_raw);

        let rn_len = streams::read_varint(&decompressed[9], &mut read_name_off);
        let read_name = &decompressed[9][read_name_off..read_name_off + rn_len];
        read_name_off += rn_len;

        let n_cigar_op = streams::read_varint(&decompressed[10], &mut cigar_off);
        let cigar_byte_len = n_cigar_op * 4;
        let cigar_bytes = &decompressed[10][cigar_off..cigar_off + cigar_byte_len];
        cigar_off += cigar_byte_len;

        let diff_count = streams::read_varint(&decompressed[11], &mut seq_diff_off);
        let diff_packed_len = (diff_count + 1) / 2;
        let diff_nibbles = streams::unpack_nibble_pairs(
            &decompressed[11][seq_diff_off..seq_diff_off + diff_packed_len],
            diff_count,
        );
        seq_diff_off += diff_packed_len;

        let extra_count = streams::read_varint(&decompressed[12], &mut seq_extra_off);
        let extra_packed_len = (extra_count + 1) / 2;
        let extra_nibbles = streams::unpack_nibble_pairs(
            &decompressed[12][seq_extra_off..seq_extra_off + extra_packed_len],
            extra_count,
        );
        seq_extra_off += extra_packed_len;

        let l_seq = streams::read_varint(&decompressed[13], &mut qual_off);
        let qual_bytes = &decompressed[13][qual_off..qual_off + l_seq];
        qual_off += l_seq;

        let aux_len = streams::read_varint(&decompressed[14], &mut aux_off);
        let aux_bytes = &decompressed[14][aux_off..aux_off + aux_len];
        aux_off += aux_len;

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

/// Decompress all 15 streams in parallel groups.
fn decompress_streams(compressed: &[Vec<u8>], use_quality_ctx: bool) -> Result<Vec<Vec<u8>>> {
    let ((group_a, group_b), (group_c, group_d)) = rayon::join(
        || {
            rayon::join(
                || decompress_group(&compressed[0..6]),
                || decompress_group(&compressed[6..9]),
            )
        },
        || {
            rayon::join(
                || decompress_group(&compressed[9..13]),
                || {
                    if use_quality_ctx {
                        let qual_passthrough = Vec::new();
                        let aux = if compressed[14].is_empty() {
                            Vec::new()
                        } else {
                            bsc::decompress_parallel(&compressed[14])?
                        };
                        Ok(vec![qual_passthrough, aux])
                    } else {
                        decompress_group(&compressed[13..15])
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

fn decompress_group(compressed: &[Vec<u8>]) -> Result<Vec<Vec<u8>>> {
    let mut results = Vec::with_capacity(compressed.len());
    for data in compressed {
        let decompressed = if data.is_empty() {
            Vec::new()
        } else {
            bsc::decompress_parallel(data)?
        };
        results.push(decompressed);
    }
    Ok(results)
}
