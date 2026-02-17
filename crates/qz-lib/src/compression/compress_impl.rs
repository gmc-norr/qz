//! Compression orchestrators: chunked pipelined compression with disk-backed block storage.

use anyhow::Result;
use std::time::Instant;
use tracing::info;
use crate::cli::{CompressConfig, QualityCompressor, QualityMode};
use super::*;
use super::codecs;

/// Unified chunked compression pipeline.
///
/// Replaces the former compress_chunked_bsc, compress_chunked_bsc_reorder_local,
/// and compress_chunked_fqzcomp with a single function parameterised by `sort_chunks`.
///
/// - `sort_chunks=false`: 2.5M chunk size, in-memory block accumulation (default path)
/// - `sort_chunks=true`: 5M chunk size, temp-file block storage (local reorder path)
///
/// Quality strategy is auto-detected:
/// - Fqzcomp: if quality_compressor == Fqzcomp
/// - QualityCtx: if explicit or auto-detected (lossless + ≥100K reads)
/// - BSC: otherwise
///
/// Features: pipelined I/O (read next chunk while compressing current), parallel
/// h||s||q compression, constant-length detection, cross-chunk validation.
pub(super) fn compress_chunked(args: &CompressConfig, sort_chunks: bool) -> Result<()> {
    use std::io::Write;

    let chunk_size: usize = if sort_chunks { 5_000_000 } else { 2_500_000 };
    let start_time = Instant::now();
    let input_path = &args.input[0];
    let bsc_static = args.advanced.bsc_static;
    let no_quality = args.no_quality;
    let quality_mode = args.quality_mode;
    let sequence_hints = args.advanced.sequence_hints;
    let sequence_delta = args.advanced.sequence_delta;
    let rc_canon = args.advanced.rc_canon;
    let quality_binning = if no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(quality_mode)
    };
    let use_fqzcomp = args.advanced.quality_compressor == QualityCompressor::Fqzcomp;

    if sort_chunks {
        info!("Chunked compression (sorted): {} records per chunk", chunk_size);
    } else {
        info!("Chunked compression: {} records per chunk", chunk_size);
    }

    let mut reader = crate::io::FastqReader::from_path_or_stdin(input_path, args.fasta)?;

    // Read first chunk
    let (mut cur_records, mut cur_bases, mut cur_orig) =
        read_chunk_records(&mut reader, chunk_size)?;

    if cur_records.is_empty() {
        anyhow::bail!("Empty input file");
    }

    // Decide quality strategy from first chunk (must be consistent across all chunks)
    let collect_for_ctx = !no_quality && quality_mode == QualityMode::Lossless && !use_fqzcomp;
    let use_quality_ctx = collect_for_ctx
        && (args.advanced.quality_compressor == QualityCompressor::QualityCtx
            || cur_records.len() >= MIN_READS_QUALITY_CTX);

    // Detect constant lengths from first chunk
    // const_qual_len only applies to BSC quality path (fqzcomp/quality_ctx have their own framing)
    let uses_packed_quality = !use_fqzcomp && !use_quality_ctx;
    let (global_const_seq_len, global_const_qual_len) = if !sort_chunks && !sequence_delta {
        let (seq, qual) = detect_const_lengths(&cur_records, no_quality);
        (seq, if uses_packed_quality { qual } else { 0 })
    } else {
        (0, 0) // reorder path doesn't use const-length (records get reordered across chunks)
    };
    if global_const_seq_len > 0 {
        info!("Constant sequence length detected: {} bp", global_const_seq_len);
    }
    if global_const_qual_len > 0 {
        info!("Constant quality length detected: {} bp", global_const_qual_len);
    }
    let quality_compressor_used = if use_fqzcomp {
        QualityCompressor::Fqzcomp
    } else if use_quality_ctx {
        info!("Using context-adaptive quality compression (quality_ctx)");
        QualityCompressor::QualityCtx
    } else {
        QualityCompressor::Bsc
    };

    // For fqzcomp/quality_ctx, skip quality in records_to_streams (compressed separately)
    let skip_quality_in_stream = no_quality || use_fqzcomp || use_quality_ctx;
    // Quality binning for stream building (fqzcomp uses raw ASCII, quality_ctx uses raw bytes)
    let stream_quality_binning = if use_fqzcomp { QualityBinning::None } else { quality_binning };

    // Block accumulation: in-memory for non-sorted, temp files for sorted
    let mut all_h_blocks: Vec<Vec<u8>> = Vec::new();
    let mut all_s_blocks: Vec<Vec<u8>> = Vec::new();
    let mut all_q_blocks: Vec<Vec<u8>> = Vec::new();
    let mut all_rc_blocks: Vec<Vec<u8>> = Vec::new();

    // Temp file infrastructure (only used when sort_chunks=true)
    struct TmpCleanup(Vec<std::path::PathBuf>);
    impl Drop for TmpCleanup {
        fn drop(&mut self) { for p in &self.0 { let _ = std::fs::remove_file(p); } }
    }
    let working_dir = &args.working_dir;
    let h_tmp_path = working_dir.join(".qz_chunked_h.tmp");
    let s_tmp_path = working_dir.join(".qz_chunked_s.tmp");
    let q_tmp_path = working_dir.join(".qz_chunked_q.tmp");
    let rc_tmp_path = working_dir.join(".qz_chunked_rc.tmp");

    let (mut h_tmp, mut s_tmp, mut q_tmp, mut rc_tmp, _cleanup);
    let (mut h_num_blocks, mut s_num_blocks, mut q_num_blocks, mut rc_num_blocks) = (0u32, 0u32, 0u32, 0u32);

    if sort_chunks {
        let mut tmp_paths = vec![h_tmp_path.clone(), s_tmp_path.clone(), q_tmp_path.clone()];
        if rc_canon { tmp_paths.push(rc_tmp_path.clone()); }
        _cleanup = Some(TmpCleanup(tmp_paths));
        h_tmp = Some(std::io::BufWriter::new(std::fs::File::create(&h_tmp_path)?));
        s_tmp = Some(std::io::BufWriter::new(std::fs::File::create(&s_tmp_path)?));
        q_tmp = Some(std::io::BufWriter::new(std::fs::File::create(&q_tmp_path)?));
        rc_tmp = if rc_canon { Some(std::io::BufWriter::new(std::fs::File::create(&rc_tmp_path)?)) } else { None };
    } else {
        _cleanup = None;
        h_tmp = None;
        s_tmp = None;
        q_tmp = None;
        rc_tmp = None;
    }

    let mut num_reads: usize = 0;
    let mut total_bases: usize = 0;
    let mut original_size: usize = 0;
    let mut chunk_idx: usize = 0;

    while !cur_records.is_empty() {
        let cur_reads = cur_records.len();

        // Optional sort
        if sort_chunks {
            info!("Chunk {}: {} reads, sorting...", chunk_idx, cur_reads);
            let sort_start = Instant::now();
            cur_records = sort_records_by_key(cur_records);
            info!("  Sorted in {:.2}s", sort_start.elapsed().as_secs_f64());
        } else {
            info!("Chunk {}: {} reads", chunk_idx, cur_reads);
        }

        // Pipeline: compress current chunk + read next chunk
        let (next_result, compress_result) = std::thread::scope(|scope| {
            let records = std::mem::take(&mut cur_records);

            let compress_handle = scope.spawn(move || -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>)> {
                let chunk_t0 = Instant::now();

                // Build streams from records
                let (h_data, s_data, q_data, rc_data) = records_to_streams(
                    &records, quality_mode, stream_quality_binning, skip_quality_in_stream,
                    sequence_hints, sequence_delta, rc_canon,
                    global_const_seq_len, global_const_qual_len,
                )?;

                // Parallel compress: headers ‖ (sequences ‖ quality)
                let (h_result, (s_result, q_result)) = rayon::join(
                    || -> Result<Vec<Vec<u8>>> {
                        let t = Instant::now();
                        let r = compress_stream_to_bsc_blocks(&h_data, bsc_static);
                        info!("  BSC headers: {:.2}s", t.elapsed().as_secs_f64());
                        r
                    },
                    || rayon::join(
                    || -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>)> {
                        let t = Instant::now();
                        let s_blocks = compress_stream_to_bsc_blocks(&s_data, bsc_static)?;
                        info!("  BSC sequences: {:.2}s", t.elapsed().as_secs_f64());
                        let rc_blocks = if rc_data.is_empty() {
                            Vec::new()
                        } else {
                            compress_stream_to_bsc_blocks(&rc_data, bsc_static)?
                        };
                        Ok((s_blocks, rc_blocks))
                    },
                    || -> Result<Vec<Vec<u8>>> {
                        let qt = Instant::now();
                        let qr = if no_quality || (q_data.is_empty() && !use_fqzcomp && !use_quality_ctx) {
                            Ok(Vec::new())
                        } else if use_fqzcomp {
                            let q_blob = codecs::compress_qualities_fqzcomp(&records)?;
                            Ok(vec![q_blob])
                        } else if use_quality_ctx {
                            use rayon::prelude::*;
                            let qual_refs: Vec<&[u8]> = records.iter()
                                .filter_map(|r| r.quality.as_deref())
                                .collect();
                            let seq_refs: Vec<&[u8]> = records.iter()
                                .map(|r| r.sequence.as_slice())
                                .collect();
                            let sub_block_reads: usize = args.advanced.quality_ctx_block_size;
                            let n = qual_refs.len();
                            if n <= sub_block_reads {
                                let blob = quality_ctx::compress_qualities_ctx(&qual_refs, &seq_refs)?;
                                Ok(vec![blob])
                            } else {
                                let num_sub = (n + sub_block_reads - 1) / sub_block_reads;
                                let blobs: Vec<Vec<u8>> = (0..num_sub)
                                    .into_par_iter()
                                    .map(|i| {
                                        let start = i * sub_block_reads;
                                        let end = (start + sub_block_reads).min(n);
                                        quality_ctx::compress_qualities_ctx(
                                            &qual_refs[start..end], &seq_refs[start..end],
                                        )
                                    })
                                    .collect::<Result<Vec<_>>>()?;
                                Ok(blobs)
                            }
                        } else {
                            compress_stream_to_bsc_blocks(&q_data, bsc_static)
                        };
                        info!("  Quality: {:.2}s", qt.elapsed().as_secs_f64());
                        qr
                    },
                    ),
                );

                let h_blocks = h_result?;
                let (s_blocks, rc_blocks) = s_result?;
                let q_blocks = q_result?;

                info!("  Chunk compress total: {:.2}s", chunk_t0.elapsed().as_secs_f64());
                Ok((h_blocks, s_blocks, q_blocks, rc_blocks))
            });

            // Read next chunk on main thread (overlaps with compression)
            let io_t = Instant::now();
            let next = read_chunk_records(&mut reader, chunk_size);
            info!("  I/O read: {:.2}s", io_t.elapsed().as_secs_f64());

            let compressed = match compress_handle.join() {
                Ok(v) => v,
                Err(e) => std::panic::resume_unwind(e),
            };
            (next, compressed)
        });

        let (h_blk, s_blk, q_blk, rc_blk) = compress_result?;

        // Accumulate blocks
        if sort_chunks {
            h_num_blocks += write_blocks_to_tmp(h_blk, h_tmp.as_mut().unwrap())?;
            s_num_blocks += write_blocks_to_tmp(s_blk, s_tmp.as_mut().unwrap())?;
            if use_fqzcomp {
                for q_blob in q_blk {
                    if !q_blob.is_empty() {
                        let qt = q_tmp.as_mut().unwrap();
                        qt.write_all(&(q_blob.len() as u32).to_le_bytes())?;
                        qt.write_all(&q_blob)?;
                        q_num_blocks += 1;
                    }
                }
            } else {
                q_num_blocks += write_blocks_to_tmp(q_blk, q_tmp.as_mut().unwrap())?;
            }
            if let Some(ref mut rc_file) = rc_tmp {
                rc_num_blocks += write_blocks_to_tmp(rc_blk, rc_file)?;
            }
        } else {
            all_h_blocks.extend(h_blk);
            all_s_blocks.extend(s_blk);
            all_q_blocks.extend(q_blk);
            if rc_canon { all_rc_blocks.extend(rc_blk); }
        }

        num_reads += cur_reads;
        total_bases += cur_bases;
        original_size += cur_orig;
        chunk_idx += 1;

        let (nr, nb, no) = next_result?;
        cur_records = nr;
        cur_bases = nb;
        cur_orig = no;

        // Validate constant-length consistency across chunks
        if !cur_records.is_empty() && !sort_chunks {
            let (chunk_const_seq, chunk_const_qual) = detect_const_lengths(&cur_records, no_quality);
            if global_const_seq_len > 0 && chunk_const_seq != global_const_seq_len {
                anyhow::bail!("Constant sequence length mismatch: chunk 0 had {} but chunk {} has {} (or variable). \
                    This shouldn't happen with well-formed Illumina data.", global_const_seq_len, chunk_idx, chunk_const_seq);
            }
            if global_const_qual_len > 0 && chunk_const_qual != global_const_qual_len {
                anyhow::bail!("Constant quality length mismatch: chunk 0 had {} but chunk {} has {} (or variable). \
                    This shouldn't happen with well-formed Illumina data.", global_const_qual_len, chunk_idx, chunk_const_qual);
            }
        }
    }

    let encoding_type: u8 = if rc_canon { 6 } else if sequence_delta { 5 } else if sequence_hints { 4 } else { 0 };

    info!("Read {} records in {} chunks ({} bases)", num_reads, chunk_idx, total_bases);

    if sort_chunks {
        // Flush and close temp files, then write archive
        for tmp in [h_tmp.as_mut(), s_tmp.as_mut(), q_tmp.as_mut()] {
            if let Some(t) = tmp { t.flush()?; }
        }
        if let Some(ref mut rc_file) = rc_tmp {
            rc_file.flush()?;
        }
        drop(h_tmp); drop(s_tmp); drop(q_tmp); drop(rc_tmp);

        let h_tmp_size = std::fs::metadata(&h_tmp_path)?.len() as usize;
        let s_tmp_size = std::fs::metadata(&s_tmp_path)?.len() as usize;
        let q_tmp_size = std::fs::metadata(&q_tmp_path)?.len() as usize;
        let rc_tmp_size = if rc_canon { std::fs::metadata(&rc_tmp_path)?.len() as usize } else { 0 };

        let headers_len = if h_num_blocks > 0 { 4 + h_tmp_size } else { 0 };
        let sequences_len = if s_num_blocks > 0 { 4 + s_tmp_size } else { 0 };
        let qualities_len = if q_num_blocks > 0 { 4 + q_tmp_size } else { 0 };
        let rc_flags_len = if rc_num_blocks > 0 { 4 + rc_tmp_size } else { 0 };

        let rc_stream = if rc_canon {
            Some(RcStreamParams { flags_len: rc_flags_len, num_blocks: rc_num_blocks, tmp_path: &rc_tmp_path })
        } else {
            None
        };
        write_chunked_archive(
            &args.output, stream_quality_binning, quality_compressor_used, no_quality, encoding_type,
            num_reads, headers_len, sequences_len, qualities_len,
            h_num_blocks, s_num_blocks, q_num_blocks,
            &h_tmp_path, &s_tmp_path, &q_tmp_path,
            original_size, start_time, rc_stream,
            global_const_seq_len, global_const_qual_len,
        )
    } else {
        // Write archive directly from in-memory blocks
        let h_data_size: usize = all_h_blocks.iter().map(|b| 4 + b.len()).sum();
        let s_data_size: usize = all_s_blocks.iter().map(|b| 4 + b.len()).sum();
        let q_data_size: usize = all_q_blocks.iter().map(|b| 4 + b.len()).sum();
        let rc_data_size: usize = all_rc_blocks.iter().map(|b| 4 + b.len()).sum();

        let headers_len = if !all_h_blocks.is_empty() { 4 + h_data_size } else { 0 };
        let sequences_len = if !all_s_blocks.is_empty() { 4 + s_data_size } else { 0 };
        let qualities_len = if !all_q_blocks.is_empty() { 4 + q_data_size } else { 0 };
        let rc_flags_len = if !all_rc_blocks.is_empty() { 4 + rc_data_size } else { 0 };

        info!("Compressed blocks: headers={} ({}) seq={} ({}) qual={} ({})",
            all_h_blocks.len(), humanize_bytes(h_data_size),
            all_s_blocks.len(), humanize_bytes(s_data_size),
            all_q_blocks.len(), humanize_bytes(q_data_size));
        if rc_canon {
            info!("  RC flags: {} ({})", all_rc_blocks.len(), humanize_bytes(rc_data_size));
        }

        info!("Writing output file...");
        let mut output_file: Box<dyn Write> = if crate::cli::is_stdio_path(&args.output) {
            Box::new(std::io::BufWriter::with_capacity(4 * 1024 * 1024, std::io::stdout().lock()))
        } else {
            Box::new(std::io::BufWriter::new(std::fs::File::create(&args.output)?))
        };

        let metadata_size = write_archive_header(
            &mut output_file, encoding_type, stream_quality_binning, quality_compressor_used,
            no_quality, num_reads, headers_len, sequences_len, qualities_len,
            global_const_seq_len, global_const_qual_len, rc_flags_len,
        )?;

        // Write streams in multi-block format: [num_blocks: u32][block_len: u32, block_data]...
        let write_blocks = |blocks: &[Vec<u8>], out: &mut dyn Write| -> Result<()> {
            if blocks.is_empty() { return Ok(()); }
            out.write_all(&(blocks.len() as u32).to_le_bytes())?;
            for block in blocks {
                out.write_all(&(block.len() as u32).to_le_bytes())?;
                out.write_all(block)?;
            }
            Ok(())
        };

        write_blocks(&all_h_blocks, &mut *output_file)?;
        write_blocks(&all_s_blocks, &mut *output_file)?;
        write_blocks(&all_q_blocks, &mut *output_file)?;

        if !all_rc_blocks.is_empty() {
            output_file.write_all(&(rc_flags_len as u64).to_le_bytes())?;
            write_blocks(&all_rc_blocks, &mut *output_file)?;
        }
        output_file.flush()?;

        log_compression_stats(original_size, headers_len, sequences_len, qualities_len, rc_flags_len, metadata_size, start_time.elapsed());

        Ok(())
    }
}

/// Global reorder: two-pass bucket sort with bounded memory per bucket.
///
/// Pass 1: Stream through FASTQ, compute sort key for each record, write to one of 256
/// bucket temp files based on top bits of the sort key.
/// Pass 2: For each bucket in order, read all records, sort by full key, build streams,
/// compress to BSC blocks, and write to output temp files.
///
/// Memory is bounded by the largest bucket (~input_size/256 + BSC working memory).
pub(super) fn compress_global_reorder_bsc(args: &CompressConfig) -> Result<()> {
    use std::io::{Write, Read as IoRead, BufWriter, BufReader};
    const NUM_BUCKETS: usize = 256;

    let start_time = Instant::now();
    let input_path = &args.input[0];
    let bsc_static = args.advanced.bsc_static;
    let no_quality = args.no_quality;
    let quality_mode = args.quality_mode;
    let quality_binning = if no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(quality_mode)
    };

    info!("Global reorder mode: two-pass bucket sort with {} buckets", NUM_BUCKETS);

    let working_dir = &args.working_dir;

    // --- Temp file setup ---
    // Bucket files for pass 1
    let bucket_paths: Vec<std::path::PathBuf> = (0..NUM_BUCKETS)
        .map(|i| working_dir.join(format!(".qz_bucket_{:03}.tmp", i)))
        .collect();
    // Output stream temp files (same as chunked mode)
    let h_tmp_path = working_dir.join(".qz_reorder_h.tmp");
    let s_tmp_path = working_dir.join(".qz_reorder_s.tmp");
    let q_tmp_path = working_dir.join(".qz_reorder_q.tmp");

    // Cleanup guard for ALL temp files
    struct TmpCleanup(Vec<std::path::PathBuf>);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            for p in &self.0 {
                let _ = std::fs::remove_file(p);
            }
        }
    }
    let mut cleanup_paths: Vec<std::path::PathBuf> = bucket_paths.clone();
    cleanup_paths.push(h_tmp_path.clone());
    cleanup_paths.push(s_tmp_path.clone());
    cleanup_paths.push(q_tmp_path.clone());
    let _cleanup = TmpCleanup(cleanup_paths);

    // --- Pass 1: bucket records by sort key ---
    info!("Pass 1: bucketing records by sort key...");
    let pass1_start = Instant::now();

    let mut bucket_writers: Vec<BufWriter<std::fs::File>> = bucket_paths.iter()
        .map(|p| Ok(BufWriter::new(std::fs::File::create(p)?)))
        .collect::<Result<Vec<_>>>()?;
    let mut bucket_counts = vec![0usize; NUM_BUCKETS];
    let mut num_reads: usize = 0;
    let mut total_bases: usize = 0;
    let mut original_size: usize = 0;

    {
        let mut reader = crate::io::FastqReader::from_path_or_stdin(input_path, args.fasta)?;
        while let Some(record) = reader.next()? {
            let qual_bytes = record.quality.as_ref().map(|q| q.as_slice()).unwrap_or(&[]);
            let key = dna_utils::reorder_sort_key(&record.sequence, qual_bytes);
            let bucket = (key >> 56) as usize % NUM_BUCKETS;

            // Write record to bucket file: [key:16][id_len:u32][id][seq_len:u32][seq][qual_len:u32][qual]
            let w = &mut bucket_writers[bucket];
            w.write_all(&key.to_le_bytes())?;

            w.write_all(&(record.id.len() as u32).to_le_bytes())?;
            w.write_all(&record.id)?;

            w.write_all(&(record.sequence.len() as u32).to_le_bytes())?;
            w.write_all(&record.sequence)?;

            let qual_len = record.quality.as_ref().map(|q| q.len()).unwrap_or(0);
            w.write_all(&(qual_len as u32).to_le_bytes())?;
            if let Some(ref q) = record.quality {
                w.write_all(q)?;
            }

            total_bases += record.sequence.len();
            original_size += record.id.len() + record.sequence.len() + qual_len + 3;
            bucket_counts[bucket] += 1;
            num_reads += 1;
        }
    }

    // Flush and close bucket writers
    for w in &mut bucket_writers {
        w.flush()?;
    }
    drop(bucket_writers);

    let non_empty_buckets = bucket_counts.iter().filter(|&&c| c > 0).count();
    let max_bucket = bucket_counts.iter().max().copied().unwrap_or(0);
    info!("Pass 1 done in {:.2}s: {} reads into {} non-empty buckets (max bucket: {} reads)",
        pass1_start.elapsed().as_secs_f64(), num_reads, non_empty_buckets, max_bucket);

    // --- Pass 2: sort each bucket, compress ---
    info!("Pass 2: sorting buckets and compressing...");
    let pass2_start = Instant::now();

    let mut h_tmp = BufWriter::new(std::fs::File::create(&h_tmp_path)?);
    let mut s_tmp = BufWriter::new(std::fs::File::create(&s_tmp_path)?);
    let mut q_tmp = BufWriter::new(std::fs::File::create(&q_tmp_path)?);

    let mut h_num_blocks: u32 = 0;
    let mut s_num_blocks: u32 = 0;
    let mut q_num_blocks: u32 = 0;

    for bucket_idx in 0..NUM_BUCKETS {
        let count = bucket_counts[bucket_idx];
        if count == 0 {
            continue;
        }

        // Read all records from this bucket file
        let mut keyed_records: Vec<(u128, crate::io::FastqRecord)> = Vec::with_capacity(count);
        {
            let file = std::fs::File::open(&bucket_paths[bucket_idx])?;
            let mut br = BufReader::new(file);

            for _ in 0..count {
                let mut key_buf = [0u8; 16];
                br.read_exact(&mut key_buf)?;
                let key = u128::from_le_bytes(key_buf);

                let mut len_buf = [0u8; 4];

                br.read_exact(&mut len_buf)?;
                let id_len = u32::from_le_bytes(len_buf) as usize;
                let mut id_bytes = vec![0u8; id_len];
                br.read_exact(&mut id_bytes)?;

                br.read_exact(&mut len_buf)?;
                let seq_len = u32::from_le_bytes(len_buf) as usize;
                let mut seq_bytes = vec![0u8; seq_len];
                br.read_exact(&mut seq_bytes)?;

                br.read_exact(&mut len_buf)?;
                let qual_len = u32::from_le_bytes(len_buf) as usize;
                let quality = if qual_len > 0 {
                    let mut qual_bytes = vec![0u8; qual_len];
                    br.read_exact(&mut qual_bytes)?;
                    Some(qual_bytes)
                } else {
                    None
                };

                keyed_records.push((key, crate::io::FastqRecord {
                    id: id_bytes,
                    sequence: seq_bytes,
                    quality,
                }));
            }
        }

        // Sort by full key within bucket
        keyed_records.sort_unstable_by_key(|(k, _)| *k);

        // Build streams and compress (sequential: h -> s -> q)
        let records: Vec<crate::io::FastqRecord> = keyed_records.into_iter().map(|(_, r)| r).collect();
        let (h_data, s_data, q_data, _rc_unused) = records_to_streams(&records, quality_mode, quality_binning, no_quality, args.advanced.sequence_hints, args.advanced.sequence_delta, false, 0, 0)?;
        drop(records);

        let h_blocks = compress_stream_to_bsc_blocks(&h_data, bsc_static)?;
        drop(h_data);
        h_num_blocks += write_blocks_to_tmp(h_blocks, &mut h_tmp)?;

        let s_blocks = compress_stream_to_bsc_blocks(&s_data, bsc_static)?;
        drop(s_data);
        s_num_blocks += write_blocks_to_tmp(s_blocks, &mut s_tmp)?;

        if !no_quality && !q_data.is_empty() {
            let q_blocks = compress_stream_to_bsc_blocks(&q_data, bsc_static)?;
            q_num_blocks += write_blocks_to_tmp(q_blocks, &mut q_tmp)?;
        }

        // Delete bucket file immediately to free disk space
        let _ = std::fs::remove_file(&bucket_paths[bucket_idx]);
    }

    info!("Pass 2 done in {:.2}s", pass2_start.elapsed().as_secs_f64());

    // Flush and close output temp files
    h_tmp.flush()?;
    s_tmp.flush()?;
    q_tmp.flush()?;
    drop(h_tmp);
    drop(s_tmp);
    drop(q_tmp);

    let h_tmp_size = std::fs::metadata(&h_tmp_path)?.len() as usize;
    let s_tmp_size = std::fs::metadata(&s_tmp_path)?.len() as usize;
    let q_tmp_size = std::fs::metadata(&q_tmp_path)?.len() as usize;

    let headers_len = if h_num_blocks > 0 { 4 + h_tmp_size } else { 0 };
    let sequences_len = if s_num_blocks > 0 { 4 + s_tmp_size } else { 0 };
    let qualities_len = if q_num_blocks > 0 { 4 + q_tmp_size } else { 0 };

    info!("Read {} records in {} non-empty buckets ({} bases)", num_reads, non_empty_buckets, total_bases);

    // Write final output (identical archive format)
    let encoding_type: u8 = if args.advanced.sequence_delta { 5 } else if args.advanced.sequence_hints { 4 } else { 0 };
    write_chunked_archive(
        &args.output, quality_binning, QualityCompressor::Bsc, no_quality, encoding_type,
        num_reads, headers_len, sequences_len, qualities_len,
        h_num_blocks, s_num_blocks, q_num_blocks,
        &h_tmp_path, &s_tmp_path, &q_tmp_path,
        original_size, start_time, None,
        0, 0, // const_seq_len, const_qual_len: bucket reorder path uses records_to_streams with varints
    )
}
pub(super) fn compress(args: &CompressConfig) -> Result<()> {
    let start_time = Instant::now();

    // Set up thread pool
    let num_threads = if args.threads == 0 {
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8)
    } else {
        args.threads
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .ok(); // Ignore error if already initialized

    info!("Using {} threads for compression", num_threads);
    info!("Compression level: {}", args.advanced.compression_level);
    info!("Input files: {:?}", args.input);
    info!("Output: {:?}", args.output);

    // Validate --sequence-hints compatibility
    if args.advanced.sequence_hints && args.advanced.sequence_compressor != SequenceCompressor::Bsc {
        anyhow::bail!("--sequence-hints requires BSC sequence compressor (default)");
    }

    // Validate --sequence-delta compatibility
    if args.advanced.sequence_delta {
        if args.advanced.sequence_hints {
            anyhow::bail!("--sequence-delta cannot be combined with --sequence-hints");
        }
        if args.advanced.sequence_compressor != SequenceCompressor::Bsc {
            anyhow::bail!("--sequence-delta requires BSC sequence compressor (default)");
        }
    }

    // Validate --quality-compressor quality-ctx
    if args.advanced.quality_compressor == QualityCompressor::QualityCtx {
        if args.no_quality {
            anyhow::bail!("--quality-compressor quality-ctx cannot be used with --no-quality");
        }
        if args.quality_mode != QualityMode::Lossless {
            anyhow::bail!("--quality-compressor quality-ctx requires lossless quality mode (default)");
        }
    }

    // Validate --local-reorder / --ultra / --fast-ultra compatibility
    let has_ultra = args.ultra.is_some();
    if args.advanced.local_reorder || has_ultra || args.advanced.fast_ultra {
        let mode_name = if args.advanced.local_reorder { "--local-reorder" } else if args.advanced.fast_ultra { "--fast-ultra" } else { "--ultra" };
        if (args.advanced.local_reorder as u8 + has_ultra as u8 + args.advanced.fast_ultra as u8) > 1 {
            anyhow::bail!("--local-reorder, --ultra, and --fast-ultra cannot be combined");
        }
        if args.advanced.sequence_hints || args.advanced.sequence_delta || args.advanced.rc_canon {
            anyhow::bail!("{mode_name} cannot be combined with other encoding options");
        }
        if args.advanced.sequence_compressor != SequenceCompressor::Bsc {
            anyhow::bail!("{mode_name} requires BSC sequence compressor (default)");
        }
        if args.advanced.header_compressor != HeaderCompressor::Bsc {
            anyhow::bail!("{mode_name} requires BSC header compressor (default)");
        }
    }

    if args.input.len() != 1 {
        anyhow::bail!("Expected exactly 1 input file, got {}", args.input.len());
    }

    // Chunked pipeline: BSC for headers/sequences, no advanced quality encoding
    let can_stream = args.advanced.sequence_compressor == SequenceCompressor::Bsc
        && args.advanced.header_compressor == HeaderCompressor::Bsc
        && !args.advanced.quality_modeling
        && !args.advanced.quality_delta
        && !args.advanced.dict_training
        && !args.advanced.twobit
        && !args.advanced.header_template;

    // Reorder modes require the chunked pipeline
    if let Some(mode) = args.advanced.reorder {
        if !can_stream {
            anyhow::bail!("--reorder requires BSC compressors for headers/sequences and no advanced encoding options");
        }
        return match mode {
            crate::cli::ReorderMode::Local => compress_chunked(args, true),
            crate::cli::ReorderMode::Global => compress_global_reorder_bsc(args),
        };
    }

    // Local reorder mode: center-hash grouping + delta encoding
    if args.advanced.local_reorder {
        return ultra::compress_harc(args);
    }

    // Ultra mode: level-based compression with auto-tuning
    if let Some(requested_level) = args.ultra {
        let level = ultra::resolve_ultra_level(requested_level);
        return ultra::compress_reorder_local_with_level(args, level);
    }

    // Backwards compat: --fast-ultra maps to ultra level 2
    if args.advanced.fast_ultra {
        let level = ultra::resolve_ultra_level(2);
        return ultra::compress_reorder_local_with_level(args, level);
    }

    // Unified chunked pipeline handles BSC, fqzcomp, and quality_ctx quality
    if can_stream {
        return compress_chunked(args, false);
    }

    // Non-streaming mode: read all records into memory first
    // (required for advanced quality encoding, non-BSC compressors)
    compress_in_memory(args, start_time)
}

/// In-memory compression path for advanced quality encodings and non-BSC compressors.
///
/// Reads all records into memory, applies optional delta/RLE/debruijn/arithmetic encoding,
/// then compresses all three streams (headers, sequences, qualities) in parallel.
fn compress_in_memory(args: &CompressConfig, start_time: Instant) -> Result<()> {
    use crate::io::{FastqReader, FastqRecord};
    use std::io::Write;

    info!("Reading FASTQ file...");
    let input_path = &args.input[0];
    let mut reader = FastqReader::from_path_or_stdin(input_path, false)?; // false = FASTQ format
    let mut records = Vec::new();

    while let Some(record) = reader.next()? {
        records.push(record);
    }

    info!("Read {} records", records.len());

    // Apply quality quantization if needed
    let processed_records: Vec<FastqRecord> = if !args.no_quality {
        records
            .into_iter()
            .map(|mut record| {
                if let Some(qual) = &record.quality {
                    record.quality = Some(quality::quantize_quality(qual, args.quality_mode).into_owned());
                }
                record
            })
            .collect()
    } else {
        records
    };

    // Determine encoding type and quality binning
    let encoding_type: u8 = if args.advanced.sequence_hints {
        4 // Sequence hints for BWT clustering
    } else if args.advanced.sequence_delta {
        5 // Inline delta encoding
    } else {
        0 // No special encoding
    };

    let quality_binning = if args.no_quality {
        QualityBinning::None
    } else {
        quality_mode_to_binning(args.quality_mode)
    };

    // Quality pre-computation (must complete before parallel compression block)
    // Train quality dictionary if requested
    let quality_dict_opt = if args.advanced.dict_training && !args.no_quality {
        info!("Training zstd dictionary for quality scores...");
        let quality_strings: Vec<Vec<u8>> = processed_records
            .iter()
            .filter_map(|r| r.quality.clone())
            .collect();

        if !quality_strings.is_empty() {
            let dict_size_bytes = args.advanced.dict_size * 1024; // Convert KB to bytes
            let dict = zstd_dict::train_from_quality_scores(&quality_strings, dict_size_bytes, 0.1)?;
            Some(dict)
        } else {
            None
        }
    } else {
        None
    };

    // Determine quality compression strategy
    // Priority: delta > modeling > standard
    let (quality_delta_opt, quality_model_opt) = if args.advanced.quality_delta && !args.no_quality {
        info!("Applying quality delta encoding...");
        let qual_compressed = if let Some(ref dict) = quality_dict_opt {
            codecs::compress_qualities_with_delta_and_dict(&processed_records, dict, args.advanced.compression_level)?
        } else {
            codecs::compress_qualities_with_delta(&processed_records, args.advanced.compression_level)?
        };
        (Some(qual_compressed), None)
    } else if args.advanced.quality_modeling && !args.no_quality {
        info!("Building positional quality model...");
        let (qual_compressed, model) = if let Some(ref dict) = quality_dict_opt {
            codecs::compress_qualities_with_model_and_dict(&processed_records, dict, args.advanced.compression_level)?
        } else {
            codecs::compress_qualities_with_model(&processed_records, args.advanced.compression_level, args.advanced.quality_compressor, args.advanced.bsc_static)?
        };
        (None, Some((qual_compressed, model)))
    } else {
        (None, None)
    };

    // Helper closure: resolve quality compression (priority: delta > model > dict > default)
    let bsc_static = args.advanced.bsc_static;
    let resolve_quality = || -> Result<Vec<u8>> {
        if let Some(ref qual) = quality_delta_opt {
            Ok(qual.clone())
        } else if let Some((ref qual, _)) = quality_model_opt {
            Ok(qual.clone())
        } else if args.no_quality {
            Ok(Vec::new())
        } else if let Some(ref dict) = quality_dict_opt {
            codecs::compress_qualities_with_dict(&processed_records, quality_binning, dict, args.advanced.compression_level, args.advanced.quality_compressor)
        } else {
            codecs::compress_qualities_with(&processed_records, quality_binning, args.advanced.compression_level, args.advanced.quality_compressor, bsc_static)
        }
    };

    // Compress all three streams in parallel:
    //   Left:  headers
    //   Right: sequences + qualities (with nested rayon::join)
    let (header_result, seq_qual_result) = rayon::join(
        || -> Result<(Vec<u8>, read_id::ReadIdTemplate)> {
            match args.advanced.header_compressor {
                HeaderCompressor::Bsc => {
                    if args.advanced.header_template {
                        let read_ids: Vec<String> = processed_records.iter().map(|r| String::from_utf8_lossy(&r.id).into_owned()).collect();
                        let encoded = read_id::compress_read_ids(&read_ids)?;
                        if encoded.template.prefix.is_empty() {
                            // Template analysis found no common Illumina structure, fall back to raw BSC
                            info!("No template structure found in headers, falling back to raw BSC{}...", if args.advanced.bsc_static { " (static)" } else { " (adaptive)" });
                            let compressed = codecs::compress_headers_bsc_with(&processed_records, args.advanced.bsc_static)?;
                            let dummy_template = read_id::ReadIdTemplate {
                                prefix: String::new(),
                                has_comment: false,
                                common_comment: None,
                            };
                            Ok((compressed, dummy_template))
                        } else {
                            info!("Compressing read IDs with template encoding + BSC{}...", if args.advanced.bsc_static { " (static)" } else { " (adaptive)" });
                            let compressed = bsc_compress_parallel(&encoded.encoded_data, args.advanced.bsc_static)?;
                            Ok((compressed, encoded.template))
                        }
                    } else {
                        info!("Compressing read IDs with raw BSC{}...", if args.advanced.bsc_static { " (static)" } else { " (adaptive)" });
                        let compressed = codecs::compress_headers_bsc_with(&processed_records, args.advanced.bsc_static)?;
                        let dummy_template = read_id::ReadIdTemplate {
                            prefix: String::new(),
                            has_comment: false,
                            common_comment: None,
                        };
                        Ok((compressed, dummy_template))
                    }
                }
                HeaderCompressor::Zstd => {
                    info!("Compressing read IDs with template encoding...");
                    codecs::compress_headers(&processed_records, args.advanced.compression_level)
                }
                HeaderCompressor::OpenZl => {
                    info!("Compressing read IDs with raw OpenZL...");
                    let compressed = codecs::compress_headers_openzl(&processed_records)?;
                    let dummy_template = read_id::ReadIdTemplate {
                        prefix: String::new(),
                        has_comment: false,
                        common_comment: None,
                    };
                    Ok((compressed, dummy_template))
                }
            }
        },
        || -> Result<(Vec<u8>, Vec<u8>, Vec<u8>)> {
            match args.advanced.sequence_compressor {
                    SequenceCompressor::Bsc => {
                        let twobit = args.advanced.twobit;
                        info!("Compressing sequences{} and qualities in parallel with BSC{}...",
                            if twobit { " (2-bit)" } else { "" },
                            if args.advanced.bsc_static { " (static)" } else { " (adaptive)" });
                        let bsc_static = args.advanced.bsc_static;
                        let (seq_result, qual_result) = rayon::join(
                            || if twobit {
                                codecs::compress_sequences_2bit_bsc_with(&processed_records, bsc_static)
                            } else {
                                codecs::compress_sequences_raw_bsc_with(&processed_records, bsc_static, args.advanced.sequence_hints, args.advanced.sequence_delta)
                            },
                            &resolve_quality,
                        );
                        let (sequences, nmasks) = seq_result?;
                        let qualities = qual_result?;
                        Ok((sequences, nmasks, qualities))
                    }
                    SequenceCompressor::Zstd => {
                        info!("Compressing with columnar format (N-mask lossless)...");
                        let (_h, sequences, nmasks, qualities_columnar, _stats) = compress_columnar(&processed_records, quality_binning, args.advanced.compression_level)?;
                        let qualities = if let Some(ref qual) = quality_delta_opt {
                            qual.clone()
                        } else if let Some((ref qual, _)) = quality_model_opt {
                            qual.clone()
                        } else {
                            qualities_columnar
                        };
                        Ok((sequences, nmasks, qualities))
                    }
                    SequenceCompressor::OpenZl => {
                        info!("Compressing sequences and qualities in parallel with OpenZL...");
                        let (seq_result, qual_result) = rayon::join(
                            || codecs::compress_sequences_raw_openzl(&processed_records),
                            &resolve_quality,
                        );
                        let (sequences, nmasks) = seq_result?;
                        let qualities = qual_result?;
                        Ok((sequences, nmasks, qualities))
                    }
                }
        },
    );
    let (headers, read_id_template) = header_result?;
    let (sequences, nmasks, qualities) = seq_qual_result?;

    // Compute stats after all streams are compressed
    let stats = columnar::ColumnarStats {
        num_reads: processed_records.len(),
        total_bases: processed_records.iter().map(|r| r.sequence.len()).sum(),
        original_size: processed_records.iter().map(|r| r.id.len() + r.sequence.len() + r.quality.as_ref().map(|q| q.len()).unwrap_or(0) + 3).sum(),
        headers_size: headers.len(),
        sequences_size: sequences.len(),
        qualities_size: qualities.len(),
        total_compressed: headers.len() + sequences.len() + qualities.len(),
    };

    // Build archive header body into a buffer, then write v2 prefix + body + data
    info!("Writing output file...");
    let mut hdr = Vec::with_capacity(256);

    // Header body fields
    hdr.push(encoding_type);
    hdr.push(0u8); // flags: 0x00 = no const-length fields

    hdr.push(binning_to_code(quality_binning));
    hdr.push(compressor_to_code(args.advanced.quality_compressor));
    hdr.push(seq_compressor_to_code(args.advanced.sequence_compressor));
    hdr.push(header_compressor_to_code(args.advanced.header_compressor));

    // Quality modeling flag and model
    if let Some((_, ref model)) = quality_model_opt {
        hdr.push(1); // modeling enabled
        let model_bytes = quality_model::serialize_model(model);
        hdr.extend_from_slice(&(model_bytes.len() as u16).to_le_bytes());
        hdr.extend_from_slice(&model_bytes);
    } else {
        hdr.push(0); // modeling disabled
    }

    // Quality delta flag
    hdr.push(if quality_delta_opt.is_some() { 1 } else { 0 });

    // Quality dictionary
    if let Some(ref dict) = quality_dict_opt {
        hdr.push(1); // dictionary enabled
        hdr.extend_from_slice(&(dict.len() as u32).to_le_bytes());
        hdr.extend_from_slice(dict);
    } else {
        hdr.push(0); // dictionary disabled
    }

    // Read ID template metadata
    let template_prefix_bytes = read_id_template.prefix.as_bytes();
    hdr.extend_from_slice(&(template_prefix_bytes.len() as u16).to_le_bytes());
    hdr.extend_from_slice(template_prefix_bytes);
    hdr.push(if read_id_template.has_comment { 1 } else { 0 });

    // Common comment (only when template is in use and has_comment is true)
    if read_id_template.has_comment {
        if let Some(ref cc) = read_id_template.common_comment {
            let bytes = cc.as_bytes();
            hdr.extend_from_slice(&(bytes.len() as u16).to_le_bytes());
            hdr.extend_from_slice(bytes);
        } else {
            hdr.extend_from_slice(&0u16.to_le_bytes());
        }
    }

    // Stream lengths
    hdr.extend_from_slice(&(stats.num_reads as u64).to_le_bytes());
    hdr.extend_from_slice(&(headers.len() as u64).to_le_bytes());
    hdr.extend_from_slice(&(sequences.len() as u64).to_le_bytes());
    hdr.extend_from_slice(&(nmasks.len() as u64).to_le_bytes());
    hdr.extend_from_slice(&(qualities.len() as u64).to_le_bytes());

    // Write v2 prefix + header body + data streams
    let header_size = V2_PREFIX_SIZE + hdr.len();
    let mut output_file = std::fs::File::create(&args.output)?;
    output_file.write_all(&ARCHIVE_MAGIC)?;
    output_file.write_all(&[ARCHIVE_VERSION, 0])?;
    output_file.write_all(&(header_size as u32).to_le_bytes())?;
    output_file.write_all(&hdr)?;

    output_file.write_all(&headers)?;
    output_file.write_all(&sequences)?;
    output_file.write_all(&nmasks)?;
    output_file.write_all(&qualities)?;

    let metadata_size = header_size;
    let total_compressed = headers.len() + sequences.len() + nmasks.len() + qualities.len() + metadata_size;

    let elapsed = start_time.elapsed();
    info!("Compression completed in {:.2}s", elapsed.as_secs_f64());

    // Show compression stats
    info!("Original size: {} bytes", stats.original_size);
    info!("Compressed size: {} bytes", total_compressed);
    info!("Compression ratio: {:.2}x", stats.original_size as f64 / total_compressed as f64);

    // Show detailed stream breakdown
    info!("Stream breakdown:");
    info!("  Headers:   {} bytes ({:.1}%)", headers.len(), 100.0 * headers.len() as f64 / total_compressed as f64);
    info!("  Sequences: {} bytes ({:.1}%)", sequences.len(), 100.0 * sequences.len() as f64 / total_compressed as f64);
    info!("  N-masks:   {} bytes ({:.1}%)", nmasks.len(), 100.0 * nmasks.len() as f64 / total_compressed as f64);
    info!("  Qualities: {} bytes ({:.1}%)", qualities.len(), 100.0 * qualities.len() as f64 / total_compressed as f64);
    info!("  Metadata:  {} bytes ({:.1}%)", metadata_size, 100.0 * metadata_size as f64 / total_compressed as f64);

    Ok(())
}
