//! Streaming BSC decompression: parallel block decompression through bounded channels.

use anyhow::{Context, Result};
use std::time::Instant;
use tracing::info;
use crate::cli::DecompressConfig;
use super::*;
use super::codecs;

/// Block-by-block BSC stream decoder for memory-efficient decompression.
///
/// Reads compressed BSC blocks from a stream region of the archive file,
/// decompresses them in parallel batches (using rayon), and sends
/// decompressed blocks through a bounded channel to the consumer thread.
fn stream_decompressor(
    path: &std::path::Path,
    offset: u64,
    stream_len: usize,
    tx: std::sync::mpsc::SyncSender<std::result::Result<Vec<u8>, String>>,
) {
    if let Err(e) = stream_decompressor_inner(path, offset, stream_len, &tx) {
        let _ = tx.send(Err(e.to_string()));
    }
}

fn stream_decompressor_inner(
    path: &std::path::Path,
    offset: u64,
    stream_len: usize,
    tx: &std::sync::mpsc::SyncSender<std::result::Result<Vec<u8>, String>>,
) -> Result<()> {
    use std::io::{Read, Seek, SeekFrom};
    use rayon::prelude::*;

    if stream_len == 0 {
        return Ok(());
    }

    let mut file = std::io::BufReader::with_capacity(
        4 * 1024 * 1024,
        std::fs::File::open(path)?,
    );
    file.seek(SeekFrom::Start(offset))?;

    let mut buf4 = [0u8; 4];
    file.read_exact(&mut buf4)?;
    let num_blocks = u32::from_le_bytes(buf4) as usize;

    let mut blocks_done = 0;

    while blocks_done < num_blocks {
        let batch_size = DECOMPRESS_BATCH_SIZE.min(num_blocks - blocks_done);

        // Read compressed blocks from file (sequential I/O)
        let mut compressed_blocks = Vec::with_capacity(batch_size);
        for _ in 0..batch_size {
            file.read_exact(&mut buf4)?;
            let block_len = u32::from_le_bytes(buf4) as usize;
            let mut data = vec![0u8; block_len];
            file.read_exact(&mut data)?;
            compressed_blocks.push(data);
        }

        // Decompress batch in parallel using rayon
        let decompressed: Result<Vec<Vec<u8>>> = compressed_blocks
            .into_par_iter()
            .map(|block| bsc::decompress(&block))
            .collect();
        let decompressed = decompressed?;

        // Send decompressed blocks through channel (may block if full)
        for block in decompressed {
            if tx.send(Ok(block)).is_err() {
                anyhow::bail!("Decompressor channel closed: receiver dropped before all blocks were consumed");
            }
        }

        blocks_done += batch_size;
    }

    Ok(())
}

/// Buffer that receives decompressed blocks from a background decompressor
/// thread via a bounded channel. Provides varint/bytes reading interface
/// for streaming record reconstruction.
struct ChannelStreamBuffer {
    rx: std::sync::mpsc::Receiver<std::result::Result<Vec<u8>, String>>,
    buf: Vec<u8>,
    pos: usize,
}

impl ChannelStreamBuffer {
    fn new(rx: std::sync::mpsc::Receiver<std::result::Result<Vec<u8>, String>>) -> Self {
        Self {
            rx,
            buf: Vec::new(),
            pos: 0,
        }
    }

    /// Ensure at least `min` unread bytes are available.
    /// Receives decompressed blocks from the channel as needed.
    fn fill(&mut self, min: usize) -> Result<bool> {
        while self.buf.len() - self.pos < min {
            // Compact consumed data: copy_within + truncate avoids drain's iterator overhead
            if self.pos > 0 {
                let remaining = self.buf.len() - self.pos;
                if remaining > 0 {
                    self.buf.copy_within(self.pos.., 0);
                }
                self.buf.truncate(remaining);
                self.pos = 0;
            }
            match self.rx.recv() {
                Ok(Ok(block)) => self.buf.extend_from_slice(&block),
                Ok(Err(e)) => anyhow::bail!("Stream decompressor error: {}", e),
                Err(_) => return Ok(self.buf.len() - self.pos >= min),
            }
        }
        Ok(true)
    }

    fn read_varint(&mut self) -> Result<usize> {
        let mut value = 0usize;
        let mut shift = 0;
        let mut bytes_read = 0;
        loop {
            if !self.fill(1)? {
                anyhow::bail!("Unexpected end of stream reading varint");
            }
            let byte = self.buf[self.pos];
            self.pos += 1;
            bytes_read += 1;
            value |= ((byte & 0x7F) as usize) << shift;
            if byte & 0x80 == 0 {
                return Ok(value);
            }
            shift += 7;
            if bytes_read >= MAX_VARINT_BYTES {
                anyhow::bail!("Malformed varint: too many continuation bytes");
            }
        }
    }

    fn read_bytes(&mut self, n: usize) -> Result<&[u8]> {
        if !self.fill(n)? {
            anyhow::bail!(
                "Unexpected end of stream: needed {} bytes, have {}",
                n,
                self.buf.len() - self.pos
            );
        }
        let start = self.pos;
        self.pos += n;
        Ok(&self.buf[start..start + n])
    }
}

/// Check if archive can use the streaming BSC decompression path.
/// Reads the v2 header prefix + body to check all flags.
pub(super) fn can_stream_decompress(args: &DecompressConfig) -> Result<bool> {
    use std::io::Read;

    let mut file = std::fs::File::open(&args.input)?;
    let mut header = [0u8; 60]; // v2 prefix (8) + base body (52)
    if file.read_exact(&mut header).is_err() {
        return Ok(false);
    }

    // Validate v2 magic
    if header[0..2] != ARCHIVE_MAGIC {
        anyhow::bail!("Not a QZ archive (missing magic bytes)");
    }

    let b: usize = V2_PREFIX_SIZE; // body starts after v2 prefix
    let encoding_type = header[b];
    let flags = header[b + 1];
    let arithmetic_enabled = flags & 0x01 != 0;
    let quality_compressor = header[b + 3];
    let sequence_compressor = header[b + 4];
    let header_compressor = header[b + 5];
    let quality_model_enabled = header[b + 6] != 0;
    let quality_delta_enabled = header[b + 7] != 0;
    let quality_dict_present = header[b + 8] != 0;
    let template_prefix_len = u16::from_le_bytes([header[b + 9], header[b + 10]]) as usize;
    let nmasks_len = read_le_u64(&header, b + 36)?;

    Ok((encoding_type == 0 || encoding_type == 4 || encoding_type == 6) // 0=raw, 4=raw+hints, 6=rc_canon
        && !arithmetic_enabled
        && !quality_model_enabled
        && !quality_delta_enabled
        && !quality_dict_present
        && template_prefix_len == 0
        && nmasks_len == 0
        && quality_compressor == 1  // BSC
        && sequence_compressor == 1 // BSC
        && header_compressor == 1)  // BSC
}

/// Memory-efficient parallel streaming decompression for BSC archives.
///
/// Spawns 3 decompressor threads (one per stream: headers, sequences, qualities)
/// that read and decompress BSC blocks in parallel batches using rayon.
/// Decompressed blocks are fed through bounded channels to the main thread,
/// which reconstructs FASTQ records and writes them to the output file.
///
/// Peak memory: ~300 MB (3 streams x ~4 decompressed blocks x 25 MB + output buffer)
/// regardless of input size. Uses all available CPU cores via rayon.
pub(super) fn decompress_streaming_bsc(args: &DecompressConfig) -> Result<()> {
    use std::io::{Read, Write};

    let start_time = Instant::now();

    info!("Input file: {:?}", args.input);
    info!("Output files: {:?}", args.output);
    info!("Streaming decompression mode (parallel BSC)");

    // Read v2 archive header: prefix (8) + body (52 base, +8 if const-length fields)
    let mut file = std::fs::File::open(&args.input)
        .with_context(|| format!("Failed to open archive: {:?}", args.input))?;
    let mut header = [0u8; 68]; // v2 prefix (8) + body (52) + optional const-length (8)
    file.read_exact(&mut header[..60])
        .context("Failed to read archive header")?;

    // Validate v2 magic
    if header[0..2] != ARCHIVE_MAGIC {
        anyhow::bail!("Not a QZ archive (missing magic bytes)");
    }
    let stored_header_size = read_le_u32(&header, 4)? as usize;

    let b: usize = V2_PREFIX_SIZE; // body starts after prefix
    let encoding_type = header[b];
    let flags = header[b + 1];
    let has_sequence_hints = encoding_type == 4;
    let has_rc_canon = encoding_type == 6;
    let quality_binning = code_to_binning(header[b + 2])?;
    let _has_comment = header[b + 11] != 0;

    let num_reads = read_le_u64(&header, b + 12)? as usize;
    let headers_len = read_le_u64(&header, b + 20)? as usize;
    let sequences_len = read_le_u64(&header, b + 28)? as usize;
    let nmasks_len = read_le_u64(&header, b + 36)? as usize;
    let qualities_len = read_le_u64(&header, b + 44)? as usize;

    // Read constant-length fields if present
    let const_length_present = flags & 0x02 != 0;
    let (const_seq_len, const_qual_len) = if const_length_present {
        file.read_exact(&mut header[60..68])
            .context("Failed to read constant-length header fields")?;
        let sl = read_le_u32(&header, 60)? as usize;
        let ql = read_le_u32(&header, 64)? as usize;
        if sl > 0 { info!("Constant sequence length: {} bp", sl); }
        if ql > 0 { info!("Constant quality length: {} bp", ql); }
        (sl, ql)
    } else {
        (0, 0)
    };
    let data_offset = stored_header_size;
    drop(file);

    // For RC canon archives, read the rc_flags_len from after the standard streams
    let (rc_flags_len, rc_flags_offset) = if has_rc_canon {
        let q_end = data_offset as u64 + headers_len as u64 + sequences_len as u64 + nmasks_len as u64 + qualities_len as u64;
        let mut file = std::fs::File::open(&args.input)?;
        use std::io::Seek;
        file.seek(std::io::SeekFrom::Start(q_end))?;
        let mut len_buf = [0u8; 8];
        file.read_exact(&mut len_buf)?;
        let rc_len = u64::from_le_bytes(len_buf) as usize;
        let rc_off = q_end + 8;
        (rc_len, rc_off)
    } else {
        (0, 0u64)
    };

    info!("Archive: {} reads, headers={} seq={} qual={}{}",
        num_reads, humanize_bytes(headers_len), humanize_bytes(sequences_len), humanize_bytes(qualities_len),
        if has_rc_canon { format!(" rc_flags={}", humanize_bytes(rc_flags_len)) } else { String::new() });

    // Compute stream positions in file
    let h_offset = data_offset as u64;
    let s_offset = h_offset + headers_len as u64;
    let q_offset = s_offset + sequences_len as u64 + nmasks_len as u64;
    let has_quality = qualities_len > 0;
    let bits_per_qual = quality_binning.bits_per_quality();

    // Open output file (or stdout if path is "-")
    if args.output.is_empty() {
        anyhow::bail!("No output file specified");
    }
    let output_path = &args.output[0];
    let is_stdout = crate::cli::is_stdio_path(output_path);

    // Spawn parallel decompressor threads and stream records to output
    info!("Decompressing {} records with parallel BSC...", num_reads);
    let archive_path = &args.input;

    // Use parallel gzip when --gzipped, otherwise plain buffered writer.
    // ParCompress owns internal threads, so we create it outside the scope
    // and call finish() after.
    if args.gzipped {
        use gzp::deflate::Gzip;
        use gzp::par::compress::{ParCompress, ParCompressBuilder};
        use gzp::Compression;
        use gzp::ZWriter;

        let num_gz_threads = (args.num_threads / 2).max(2);
        // ParCompress::from_writer needs Send + 'static; Stdout qualifies, StdoutLock does not
        let mut output: ParCompress<Gzip> = if is_stdout {
            ParCompressBuilder::new()
                .num_threads(num_gz_threads)
                .map_err(|e| anyhow::anyhow!("gzp error: {e}"))?
                .compression_level(Compression::new(args.gzip_level))
                .from_writer(std::io::stdout())
        } else {
            let file = std::fs::File::create(output_path)
                .with_context(|| format!("Failed to create output file: {:?}", output_path))?;
            ParCompressBuilder::new()
                .num_threads(num_gz_threads)
                .map_err(|e| anyhow::anyhow!("gzp error: {e}"))?
                .compression_level(Compression::new(args.gzip_level))
                .from_writer(file)
        };

        stream_records_to_writer(&mut output, num_reads, archive_path, h_offset, headers_len, s_offset, sequences_len, q_offset, qualities_len, has_quality, const_seq_len, const_qual_len, has_sequence_hints, bits_per_qual, quality_binning, has_rc_canon, rc_flags_len, rc_flags_offset)?;

        output.finish().map_err(|e| anyhow::anyhow!("gzp finish error: {e}"))?;
    } else if is_stdout {
        let mut output = std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, std::io::stdout().lock());

        stream_records_to_writer(&mut output, num_reads, archive_path, h_offset, headers_len, s_offset, sequences_len, q_offset, qualities_len, has_quality, const_seq_len, const_qual_len, has_sequence_hints, bits_per_qual, quality_binning, has_rc_canon, rc_flags_len, rc_flags_offset)?;

        output.flush()?;
    } else {
        let file = std::fs::File::create(output_path)
            .with_context(|| format!("Failed to create output file: {:?}", output_path))?;
        let mut output = std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, file);

        stream_records_to_writer(&mut output, num_reads, archive_path, h_offset, headers_len, s_offset, sequences_len, q_offset, qualities_len, has_quality, const_seq_len, const_qual_len, has_sequence_hints, bits_per_qual, quality_binning, has_rc_canon, rc_flags_len, rc_flags_offset)?;

        output.flush()?;
    }

    info!("Decompressed {} records", num_reads);
    let elapsed = start_time.elapsed();
    info!("Decompression completed in {:.2}s", elapsed.as_secs_f64());

    Ok(())
}

/// Stream decompressed records to an output writer.
///
/// Extracted from decompress_streaming_bsc to allow different writer types
/// (plain buffered vs parallel gzip) without duplicating the decompression logic.
fn stream_records_to_writer(
    output: &mut dyn std::io::Write,
    num_reads: usize,
    archive_path: &std::path::Path,
    h_offset: u64,
    headers_len: usize,
    s_offset: u64,
    sequences_len: usize,
    q_offset: u64,
    qualities_len: usize,
    has_quality: bool,
    const_seq_len: usize,
    const_qual_len: usize,
    has_sequence_hints: bool,
    bits_per_qual: usize,
    quality_binning: QualityBinning,
    has_rc_canon: bool,
    rc_flags_len: usize,
    rc_flags_offset: u64,
) -> Result<()> {
    std::thread::scope(|scope| -> Result<()> {
        // Bounded channels: capacity 2 = max 2 pre-decompressed blocks buffered per stream
        let (h_tx, h_rx) = std::sync::mpsc::sync_channel(2);
        let (s_tx, s_rx) = std::sync::mpsc::sync_channel(2);
        let h_path = archive_path.to_path_buf();
        let s_path = archive_path.to_path_buf();

        scope.spawn(move || stream_decompressor(&h_path, h_offset, headers_len, h_tx));
        scope.spawn(move || stream_decompressor(&s_path, s_offset, sequences_len, s_tx));

        let q_rx = if has_quality {
            let (q_tx, q_rx) = std::sync::mpsc::sync_channel(2);
            let q_path = archive_path.to_path_buf();
            scope.spawn(move || stream_decompressor(&q_path, q_offset, qualities_len, q_tx));
            Some(q_rx)
        } else {
            None
        };

        // RC flags stream (encoding_type=6)
        let rc_rx = if has_rc_canon && rc_flags_len > 0 {
            let (rc_tx, rc_rx) = std::sync::mpsc::sync_channel(2);
            let rc_path = archive_path.to_path_buf();
            scope.spawn(move || stream_decompressor(&rc_path, rc_flags_offset, rc_flags_len, rc_tx));
            Some(rc_rx)
        } else {
            None
        };

        let mut h_buf = ChannelStreamBuffer::new(h_rx);
        let mut s_buf = ChannelStreamBuffer::new(s_rx);
        let mut q_buf = q_rx.map(ChannelStreamBuffer::new);
        let mut rc_buf = rc_rx.map(ChannelStreamBuffer::new);

        for i in 0..num_reads {
            // Header: varint(len) + raw bytes
            let h_len = h_buf.read_varint()?;
            output.write_all(h_buf.read_bytes(h_len)?)?;
            output.write_all(b"\n")?;

            // Sequence: [varint(len)] [+ hint byte] + raw bytes
            let s_len = if const_seq_len > 0 { const_seq_len } else { s_buf.read_varint()? };
            if has_sequence_hints {
                s_buf.read_bytes(1)?; // skip hint byte
            }
            let seq_bytes = s_buf.read_bytes(s_len)?;

            // If RC-canonicalized, check flag and reverse-complement if needed
            if let Some(ref mut rc) = rc_buf {
                let flag = rc.read_bytes(1)?[0];
                if flag != 0 {
                    let revcomp = dna_utils::reverse_complement(seq_bytes);
                    output.write_all(&revcomp)?;
                } else {
                    output.write_all(seq_bytes)?;
                }
            } else {
                output.write_all(seq_bytes)?;
            }
            output.write_all(b"\n+\n")?;

            // Quality: [varint(orig_len)] + packed bytes
            if let Some(ref mut q) = q_buf {
                let q_len = if const_qual_len > 0 { const_qual_len } else { q.read_varint()? };
                let packed_len = (q_len * bits_per_qual + 7) / 8;
                let packed = q.read_bytes(packed_len)?;
                columnar::unpack_qualities_to_writer(packed, q_len, quality_binning, output)?;
            }
            output.write_all(b"\n")?;

            if (i + 1) % 10_000_000 == 0 {
                info!("  {} million records written...", (i + 1) / 1_000_000);
            }
        }

        output.flush()?;
        Ok(())
    })
}

/// Write FASTQ records to a writer in 2MB batches.
fn write_records_batched(output: &mut dyn std::io::Write, records: &[crate::io::fastq::FastqRecord]) -> Result<()> {
    const WRITE_BATCH: usize = 2 * 1024 * 1024;
    let mut buf = Vec::with_capacity(WRITE_BATCH + 1024);
    for record in records {
        buf.extend_from_slice(&record.id);
        buf.push(b'\n');
        buf.extend_from_slice(&record.sequence);
        buf.extend_from_slice(b"\n+\n");
        if let Some(qual) = &record.quality {
            buf.extend_from_slice(qual);
        }
        buf.push(b'\n');
        if buf.len() >= WRITE_BATCH {
            output.write_all(&buf)?;
            buf.clear();
        }
    }
    if !buf.is_empty() {
        output.write_all(&buf)?;
    }
    Ok(())
}

/// Parsed archive header metadata.
struct ArchiveHeader {
    encoding_type: u8,
    arithmetic_enabled: bool,
    quality_binning: QualityBinning,
    quality_compressor: QualityCompressor,
    sequence_compressor: SequenceCompressor,
    header_compressor: HeaderCompressor,
    quality_model_opt: Option<quality_model::QualityModel>,
    quality_delta_enabled: bool,
    quality_dict_opt: Option<Vec<u8>>,
    read_id_template: read_id::ReadIdTemplate,
    num_reads: usize,
    /// Byte offset where stream data begins.
    data_offset: usize,
    headers_len: usize,
    sequences_len: usize,
    nmasks_len: usize,
    qualities_len: usize,
    /// Constant sequence length (0 = variable, >0 = all sequences this length).
    const_seq_len: usize,
    /// Constant quality length (0 = variable, >0 = all qualities this length).
    const_qual_len: usize,
}

/// Parse the archive header from raw archive bytes.
///
/// Returns parsed metadata and the byte offset where stream data starts.
fn parse_archive_header(data: &[u8]) -> Result<ArchiveHeader> {
    if data.len() < V2_PREFIX_SIZE + 50 {
        anyhow::bail!("Invalid archive: too small");
    }

    // Validate v2 magic
    if data[0..2] != ARCHIVE_MAGIC {
        anyhow::bail!("Not a QZ archive (missing magic bytes)");
    }
    let _version = data[2];
    let _stored_header_size = read_le_u32(data, 4)? as usize;

    // Body starts after v2 prefix
    let mut offset = V2_PREFIX_SIZE;

    let encoding_type = data[offset];
    offset += 1;

    let flags = data[offset];
    offset += 1;
    let arithmetic_enabled = flags & 0x01 != 0;
    let const_length_present = flags & 0x02 != 0;
    // Skip read lengths if arithmetic mode was used
    if arithmetic_enabled {
        let num_lengths = read_le_u32(data, offset)? as usize;
        offset += 4;
        offset += num_lengths * 4;
    }

    let quality_binning = code_to_binning(data[offset])?;
    offset += 1;

    let quality_compressor = if offset < data.len() && data[offset] <= 4 {
        let compressor = code_to_compressor(data[offset])?;
        offset += 1;
        compressor
    } else {
        QualityCompressor::Zstd
    };

    // Fqzcomp stores raw ASCII quality bytes (not bit-packed), so override binning
    let quality_binning = if quality_compressor == QualityCompressor::Fqzcomp {
        QualityBinning::None
    } else {
        quality_binning
    };

    let sequence_compressor = code_to_seq_compressor(data[offset])?;
    offset += 1;

    let header_compressor = code_to_header_compressor(data[offset])?;
    offset += 1;

    let quality_model_enabled = data[offset] != 0;
    offset += 1;
    let quality_model_opt = if quality_model_enabled {
        let model_size = read_le_u16(data, offset)? as usize;
        offset += 2;
        let model_bytes = &data[offset..offset + model_size];
        offset += model_size;
        Some(quality_model::deserialize_model(model_bytes)?)
    } else {
        None
    };

    let quality_delta_enabled = data[offset] != 0;
    offset += 1;

    let quality_dict_opt = if data[offset] != 0 {
        offset += 1;
        let dict_size = read_le_u32(data, offset)? as usize;
        offset += 4;
        let dict = data[offset..offset + dict_size].to_vec();
        offset += dict_size;
        Some(dict)
    } else {
        offset += 1;
        None
    };

    let template_prefix_len = read_le_u16(data, offset)? as usize;
    offset += 2;
    let template_prefix = String::from_utf8_lossy(&data[offset..offset + template_prefix_len]).to_string();
    offset += template_prefix_len;
    let template_has_comment = data[offset] != 0;
    offset += 1;

    let common_comment = if template_prefix_len > 0 && template_has_comment {
        let cc_len = read_le_u16(data, offset)? as usize;
        offset += 2;
        if cc_len > 0 {
            let cc = String::from_utf8_lossy(&data[offset..offset + cc_len]).to_string();
            offset += cc_len;
            Some(cc)
        } else {
            None
        }
    } else {
        None
    };

    let read_id_template = read_id::ReadIdTemplate {
        prefix: template_prefix,
        has_comment: template_has_comment,
        common_comment,
    };

    let num_reads = read_le_u64(data, offset)? as usize;
    offset += 8;
    let headers_len = read_le_u64(data, offset)? as usize;
    offset += 8;
    let sequences_len = read_le_u64(data, offset)? as usize;
    offset += 8;
    let nmasks_len = read_le_u64(data, offset)? as usize;
    offset += 8;
    let qualities_len = read_le_u64(data, offset)? as usize;
    offset += 8;

    // Read constant-length fields if flags bit 1 is set
    let (const_seq_len, const_qual_len) = if const_length_present {
        let sl = read_le_u32(data, offset)? as usize;
        offset += 4;
        let ql = read_le_u32(data, offset)? as usize;
        offset += 4;
        (sl, ql)
    } else {
        (0, 0)
    };

    if data.len() < offset + headers_len + sequences_len + nmasks_len + qualities_len {
        anyhow::bail!("Invalid archive: data truncated");
    }

    Ok(ArchiveHeader {
        encoding_type,
        arithmetic_enabled,
        quality_binning,
        quality_compressor,
        sequence_compressor,
        header_compressor,
        quality_model_opt,
        quality_delta_enabled,
        quality_dict_opt,
        read_id_template,
        num_reads,
        data_offset: offset,
        headers_len,
        sequences_len,
        nmasks_len,
        qualities_len,
        const_seq_len,
        const_qual_len,
    })
}

pub(super) fn decompress(args: &DecompressConfig) -> Result<()> {
    use std::io::Write;

    // If input is stdin, spool to a temp file first (decompression needs seeking)
    let (_stdin_tmp, args) = if crate::cli::is_stdio_path(&args.input) {
        info!("Reading archive from stdin...");
        let mut tmp = tempfile::NamedTempFile::new_in(&args.working_dir)
            .context("Failed to create temp file for stdin")?;
        std::io::copy(&mut std::io::stdin().lock(), &mut tmp)
            .context("Failed to spool stdin to temp file")?;
        let path = tmp.path().to_path_buf();
        let mut new_args = args.clone();
        new_args.input = path;
        (Some(tmp), std::borrow::Cow::Owned(new_args))
    } else {
        (None, std::borrow::Cow::Borrowed(args))
    };
    let args = &*args;

    let start_time = Instant::now();

    // Try streaming decompression for standard BSC archives (>100x less memory)
    if can_stream_decompress(args)? {
        return decompress_streaming_bsc(args);
    }

    info!("Input file: {:?}", args.input);
    info!("Output files: {:?}", args.output);

    // Memory-map the compressed archive (no allocation, OS pages in on demand)
    info!("Memory-mapping compressed file...");
    let input_file = std::fs::File::open(&args.input)
        .with_context(|| format!("Failed to open archive: {:?}", args.input))?;
    let archive_data = unsafe { memmap2::Mmap::map(&input_file) }
        .context("Failed to mmap archive")?;

    // Parse archive header
    let hdr = parse_archive_header(&archive_data)?;
    let mut offset = hdr.data_offset;
    let headers = &archive_data[offset..offset + hdr.headers_len];
    offset += hdr.headers_len;
    let sequences = &archive_data[offset..offset + hdr.sequences_len];
    offset += hdr.sequences_len;
    let nmasks = &archive_data[offset..offset + hdr.nmasks_len];
    offset += hdr.nmasks_len;
    let qualities = &archive_data[offset..offset + hdr.qualities_len];
    let num_reads = hdr.num_reads;
    let template_prefix_len = hdr.read_id_template.prefix.len();
    let encoding_type = hdr.encoding_type;
    let arithmetic_enabled = hdr.arithmetic_enabled;
    let quality_binning = hdr.quality_binning;
    let quality_compressor = hdr.quality_compressor;
    let sequence_compressor = hdr.sequence_compressor;
    let header_compressor = hdr.header_compressor;
    let quality_delta_enabled = hdr.quality_delta_enabled;
    let qualities_len = hdr.qualities_len;
    let nmasks_len = hdr.nmasks_len;
    let read_id_template = &hdr.read_id_template;
    let quality_model_opt = &hdr.quality_model_opt;
    let quality_dict_opt = &hdr.quality_dict_opt;
    let const_seq_len = hdr.const_seq_len;
    let const_qual_len = hdr.const_qual_len;
    if const_seq_len > 0 { info!("Constant sequence length: {} bp", const_seq_len); }
    if const_qual_len > 0 { info!("Constant quality length: {} bp", const_qual_len); }

    // Decompress based on encoding mode
    info!("Decompressing...");
    let records = if arithmetic_enabled {
        anyhow::bail!("Arithmetic-coded archives are no longer supported (encoding removed)")
    } else if encoding_type == 3 {
        anyhow::bail!("De Bruijn graph archives are no longer supported (encoding_type=3 removed)")
    } else if encoding_type == 7 {
        anyhow::bail!("Factorized archives are no longer supported (encoding_type=7 removed)")
    } else if encoding_type == 8 || encoding_type == 9 {
        // Local reorder + delta (8) or ultra/big-block (9): both store permutation for original order restore
        let mode_name = if encoding_type == 8 { "local-reorder" } else { "ultra" };
        let use_quality_ctx = quality_compressor == QualityCompressor::QualityCtx;
        info!("Decompressing {} sequences (encoding_type={}, quality={})...",
            mode_name, encoding_type, if use_quality_ctx { "quality_ctx" } else { "BSC" });

        // Decode headers and sequences in parallel
        // (quality_ctx needs sequences as context, so can't decompress qualities in parallel)
        let (header_result, seq_result) = rayon::join(
            || codecs::decompress_headers_bsc(headers, num_reads)
                .context("Failed to decompress headers"),
            || -> Result<ultra::HarcDecoded> {
                if encoding_type == 8 {
                    ultra::decode_harc_sequences(sequences, num_reads)
                } else {
                    ultra::decode_reorder_local(sequences, num_reads)
                }
            },
        );

        let read_ids = header_result?;
        let decoded = seq_result?;

        // Decode qualities (may need sequences for quality_ctx)
        let mut qual_strings: Vec<Option<Vec<u8>>> = vec![None; num_reads];
        if qualities_len > 0 {
            if use_quality_ctx {
                // Quality_ctx: decompress using sequences as context (in reordered order)
                let ctx_qualities = quality_ctx::decompress_quality_ctx_multiblock(
                    qualities, &decoded.sequences,
                ).context("Failed to decompress qualities (quality_ctx)")?;
                for (i, q) in ctx_qualities.into_iter().enumerate() {
                    qual_strings[i] = Some(q);
                }
            } else {
                // BSC: decompress and parse packed quality bytes
                let qualities_data = codecs::decompress_qualities_data(qualities, quality_compressor)
                    .context("Failed to decompress qualities")?;
                let mut qual_offset = 0;
                for i in 0..num_reads {
                    if qual_offset < qualities_data.len() {
                        let q_len = if const_qual_len > 0 {
                            const_qual_len
                        } else {
                            read_varint(&qualities_data, &mut qual_offset)
                                .ok_or_else(|| anyhow::anyhow!("Failed to read quality length at read {}", i))?
                        };
                        let bits_per_qual = quality_binning.bits_per_quality();
                        let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
                        if qual_offset + q_encoded_len <= qualities_data.len() {
                            let quality_str = columnar::unpack_qualities(
                                &qualities_data[qual_offset..qual_offset + q_encoded_len],
                                q_len, quality_binning,
                            );
                            qual_offset += q_encoded_len;
                            qual_strings[i] = Some(quality_str);
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        // Build records in reordered order, then un-permute to original order
        let mut records: Vec<Option<crate::io::FastqRecord>> = (0..num_reads).map(|_| None).collect();
        for (reorder_pos, (id, seq)) in read_ids.into_iter()
            .zip(decoded.sequences.into_iter())
            .enumerate()
        {
            let orig_idx = decoded.order[reorder_pos] as usize;
            let qual = std::mem::take(&mut qual_strings[reorder_pos]);
            records[orig_idx] = Some(crate::io::FastqRecord { id, sequence: seq, quality: qual });
        }

        records.into_iter().enumerate()
            .map(|(i, r)| r.ok_or_else(|| anyhow::anyhow!("missing record at index {i} after un-permutation")))
            .collect::<Result<Vec<_>>>()?
    } else if encoding_type == 1 || encoding_type == 2 {
        anyhow::bail!("Delta/RLE archives are no longer supported (encoding_type={} removed)", encoding_type)
    } else {
        // Normal mode: decompress all streams in parallel
        info!("Decompressing streams in parallel...");
        let use_quality_ctx = quality_compressor == QualityCompressor::QualityCtx;

        // Step 1: Parallel decompression of headers, sequences, qualities
        // For quality_ctx: skip quality decompression here (needs sequences first)
        let (header_result, (seq_result, qual_result)) = rayon::join(
            || decompress_headers_dispatch(
                header_compressor, headers, template_prefix_len, &read_id_template, num_reads,
            ),
            || rayon::join(
                || -> Result<Vec<Vec<u8>>> {
                    match sequence_compressor {
                        SequenceCompressor::Bsc => {
                            if nmasks_len > 0 {
                                codecs::decompress_sequences_2bit_bsc(sequences, nmasks, num_reads, const_seq_len)
                                    .context("Failed to decompress 2-bit sequences (BSC)")
                            } else {
                                codecs::decompress_sequences_raw_bsc(sequences, num_reads, encoding_type, const_seq_len)
                                    .context("Failed to decompress sequences (BSC)")
                            }
                        }
                        SequenceCompressor::OpenZl => {
                            codecs::decompress_sequences_raw_openzl(sequences, num_reads, const_seq_len)
                                .context("Failed to decompress sequences (OpenZL)")
                        }
                        SequenceCompressor::Zstd => {
                            let sequences_data = decompress_zstd(sequences)
                                .context("Failed to decompress sequences (zstd)")?;
                            let nmasks_data = decompress_zstd(nmasks)
                                .context("Failed to decompress N-masks (zstd)")?;
                            let mut decoded = Vec::with_capacity(num_reads);
                            let mut seq_offset = 0;
                            let mut nmask_offset = 0;
                            for _ in 0..num_reads {
                                let seq_len = read_varint(&sequences_data, &mut seq_offset)
                                    .ok_or_else(|| anyhow::anyhow!("Failed to read sequence length"))?;
                                let seq_2bit_len = (seq_len + 3) / 4;
                                if seq_offset + seq_2bit_len > sequences_data.len() {
                                    anyhow::bail!("Truncated sequence data");
                                }
                                let sequence_2bit = &sequences_data[seq_offset..seq_offset + seq_2bit_len];
                                seq_offset += seq_2bit_len;
                                let nmask_len = (seq_len + 7) / 8;
                                let n_mask = if nmask_offset + nmask_len <= nmasks_data.len() {
                                    let mask = &nmasks_data[nmask_offset..nmask_offset + nmask_len];
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
                    }
                },
                || -> Result<Vec<u8>> {
                    if use_quality_ctx || qualities_len == 0 {
                        Ok(Vec::new())
                    } else if let Some(dict) = quality_dict_opt {
                        zstd_dict::decompress_with_dict(qualities, dict)
                            .context("Failed to decompress quality scores (dictionary mode)")
                    } else {
                        codecs::decompress_qualities_data(qualities, quality_compressor)
                            .context("Failed to decompress quality scores")
                    }
                },
            ),
        );

        let read_ids = header_result?;
        let decoded_sequences = seq_result?;
        let qualities_data = qual_result?;

        // For quality_ctx: decompress qualities using sequences as context
        // Must happen after sequences are available (breaks seq||qual parallelism)
        let quality_ctx_strings = if use_quality_ctx {
            info!("Decompressing quality scores (quality_ctx)...");
            Some(quality_ctx::decompress_quality_ctx_multiblock(qualities, &decoded_sequences)
                .context("Failed to decompress quality scores (quality_ctx)")?)
        } else {
            None
        };

        // Step 2: Construct records from decoded headers and sequences
        let mut records: Vec<_> = read_ids.into_iter().zip(decoded_sequences).map(|(id, seq)| {
            crate::io::FastqRecord::new(id, seq, None)
        }).collect();

        // Step 3: Decode quality strings and attach to records
        if let Some(quality_strings) = quality_ctx_strings {
            for (i, qual) in quality_strings.into_iter().enumerate() {
                if i < records.len() {
                    records[i].quality = Some(qual);
                }
            }
        } else if !qualities_data.is_empty() {
            if quality_delta_enabled {
                let mut encoded_deltas = Vec::new();
                let mut qual_offset = 0;
                while qual_offset < qualities_data.len() {
                    let q_len = if const_qual_len > 0 {
                        const_qual_len
                    } else {
                        read_varint(&qualities_data, &mut qual_offset).ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?
                    };
                    if qual_offset + q_len > qualities_data.len() {
                        break;
                    }
                    let deltas = quality_delta::unpack_deltas(&qualities_data[qual_offset..qual_offset + q_len]);
                    qual_offset += q_len;
                    encoded_deltas.push(deltas);
                }
                let decoded_qualities = quality_delta::decode_quality_deltas(&encoded_deltas)?;
                for (i, quality) in decoded_qualities.into_iter().enumerate() {
                    if i < records.len() {
                        records[i].quality = Some(quality);
                    }
                }
            } else {
                // Two-pass parallel quality decoding:
                // Pass 1: Sequential scan to locate record boundaries (fast varint parsing)
                use rayon::prelude::*;
                let bits_per_qual = quality_binning.bits_per_quality();
                let mut qual_entries: Vec<(usize, usize)> = Vec::with_capacity(records.len());
                let mut qual_offset = 0;
                for _ in 0..records.len() {
                    if qual_offset >= qualities_data.len() { break; }
                    let q_len = if const_qual_len > 0 {
                        const_qual_len
                    } else {
                        read_varint(&qualities_data, &mut qual_offset)
                            .ok_or_else(|| anyhow::anyhow!("Failed to read quality length"))?
                    };
                    let data_len = if quality_compressor == QualityCompressor::Fqzcomp || quality_model_opt.is_some() {
                        q_len
                    } else {
                        (q_len * bits_per_qual + 7) / 8
                    };
                    if qual_offset + data_len > qualities_data.len() { break; }
                    qual_entries.push((qual_offset, q_len));
                    qual_offset += data_len;
                }

                // Pass 2: Parallel decode (each record is independent)
                let decoded_quals: Vec<Vec<u8>> = qual_entries.par_iter().map(|&(data_off, q_len)| {
                    if quality_compressor == QualityCompressor::Fqzcomp {
                        let raw = &qualities_data[data_off..data_off + q_len];
                        raw.to_vec()
                    } else if let Some(model) = quality_model_opt {
                        let deltas = quality_model::unpack_deltas(&qualities_data[data_off..data_off + q_len]);
                        quality_model::decode_with_model(&deltas, model)
                    } else {
                        let q_encoded_len = (q_len * bits_per_qual + 7) / 8;
                        columnar::unpack_qualities(&qualities_data[data_off..data_off + q_encoded_len], q_len, quality_binning)
                    }
                }).collect();

                for (i, qual) in decoded_quals.into_iter().enumerate() {
                    if i < records.len() {
                        records[i].quality = Some(qual);
                    }
                }
            }
        }

        records
    };

    info!("Decompressed {} records", records.len());

    // Write output
    info!("Writing output file...");
    if args.output.is_empty() {
        anyhow::bail!("No output file specified");
    }

    let output_path = &args.output[0];
    let is_stdout = crate::cli::is_stdio_path(output_path);

    if args.gzipped {
        use gzp::deflate::Gzip;
        use gzp::par::compress::{ParCompress, ParCompressBuilder};
        use gzp::Compression;
        use gzp::ZWriter;

        let num_gz_threads = (args.num_threads / 2).max(2);
        let mut output: ParCompress<Gzip> = if is_stdout {
            ParCompressBuilder::new()
                .num_threads(num_gz_threads)
                .map_err(|e| anyhow::anyhow!("gzp error: {e}"))?
                .compression_level(Compression::new(args.gzip_level))
                .from_writer(std::io::stdout())
        } else {
            let file = std::fs::File::create(output_path)
                .with_context(|| format!("Failed to create output file: {:?}", output_path))?;
            ParCompressBuilder::new()
                .num_threads(num_gz_threads)
                .map_err(|e| anyhow::anyhow!("gzp error: {e}"))?
                .compression_level(Compression::new(args.gzip_level))
                .from_writer(file)
        };

        write_records_batched(&mut output, &records)?;
        output.finish().map_err(|e| anyhow::anyhow!("gzp finish error: {e}"))?;
    } else if is_stdout {
        let mut output = std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, std::io::stdout().lock());
        write_records_batched(&mut output, &records)?;
        output.flush()?;
    } else {
        let file = std::fs::File::create(output_path)
            .with_context(|| format!("Failed to create output file: {:?}", output_path))?;
        let mut output = std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, file);
        write_records_batched(&mut output, &records)?;
        output.flush()?;
    }

    let elapsed = start_time.elapsed();
    info!("Decompression completed in {:.2}s", elapsed.as_secs_f64());

    Ok(())
}

