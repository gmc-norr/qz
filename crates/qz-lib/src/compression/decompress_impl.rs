//! Bounded streaming BSC decompression: per-stream per-block-frame decompression.
//! Each stream is consumed block-by-block via bounded mpsc channels, keeping
//! at most ~3 decompressed blocks per active stream in memory at peak.

use super::*;
use crate::cli::{DecompressConfig, HeaderCompressor, QualityCompressor, VerifyConfig};
use anyhow::{Context, Result};
use std::time::Instant;
use tracing::info;

/// Metadata describing a decoded v5 archive, returned by the verify path and
/// consumed when building a [`VerifyResult`]. Only the fields the v5 verify path
/// actually reads are carried (the legacy v4 stream-major header carried more).
pub(super) struct ArchiveHeader {
    pub(super) encoding_type: u8,
    pub(super) quality_compressor: QualityCompressor,
    pub(super) header_compressor: HeaderCompressor,
    pub(super) num_reads: usize,
    pub(super) headers_len: usize,
    pub(super) sequences_len: usize,
    pub(super) qualities_len: usize,
}

/// Per-block compressed-size ceiling enforced by the streaming BSC
/// decompressor and the columnar-header decompressor.
///
/// The compress side caps BSC input blocks at 64 MiB
/// (`compress_impl.rs::compress()` rejects `bsc_block_size_mb > 64`). libbsc's
/// worst-case output is `input + LIBBSC_HEADER_SIZE` (28 bytes; see
/// `third_party/libbsc/libbsc/libbsc/libbsc.cpp:80`), and the rust wrapper
/// pads its output buffer by another 1024 bytes. A 4 KiB margin on top of the
/// 64 MiB input cap covers both. Without this margin a 64 MiB block of
/// incompressible data could compress to slightly over 64 MiB and become
/// undecompressible — see CHANGELOG round-5.
pub(crate) const BSC_MAX_BLOCK_LEN: usize = 64 * 1024 * 1024 + 4096;

/// Maximum decompressed BSC block size for the `--ultra` big-block path.
/// Unlike the default/bounded path (25 MB blocks, bounded by [`BSC_MAX_BLOCK_LEN`]),
/// ultra deliberately uses large blocks for better BWT context: sequence blocks
/// up to `max_seq_block_mb` (750 MB) and header/quality blocks up to 100 MB
/// (large ultra blocks, max level 3). 768 MiB gives headroom over the 750 MB
/// sequence block while still bounding a corrupt-header allocation to a single,
/// finite amount.
pub(crate) const BSC_MAX_BLOCK_LEN_ULTRA: usize = 768 * 1024 * 1024;

/// Budget used to derive the per-stream `max_inflight` block count. The
/// per-stream in-flight cap is `DECODE_MEM_BUDGET / max_block_len`, clamped
/// to `MAX_PARALLEL_DECOMPRESS`.
///
/// **Important: this is a PER-STREAM cap.** Up to ~3 streams (headers,
/// sequences, qualities) decode concurrently, so worst-case in-flight memory
/// is `max_inflight × ~3 × max_block_len`. For ultra (768 MiB cap,
/// max_inflight ≈ 5) that is roughly 5 × 3 × 768 MiB ≈ 11 GiB in the
/// worst case; steady-state is lower because writer backpressure drains
/// blocks as soon as the consumer catches up. The empirical RSS bound is
/// verified by the RSS tests added in later tasks.
///
/// Default archives (64 MiB cap) yield max_inflight ≈ 64, clamped to
/// `MAX_PARALLEL_DECOMPRESS` (16) — i.e. unchanged from the old constant.
///
/// `pub(crate)` so the reference decode (`reference::reference_decode_max_inflight`)
/// derives its bounded-parallel chunk concurrency from the same budget.
pub(crate) const DECODE_MEM_BUDGET: usize = 4 * 1024 * 1024 * 1024;

/// Derive the per-archive `(max_block_len, max_inflight)` decode bounds from the
/// archive's `encoding_type`. `UltraBigBlock` (10) uses the 768 MiB block cap
/// and a smaller in-flight count; everything else keeps the 64 MiB cap and the
/// historical concurrency.
pub(super) fn decode_block_bounds(encoding_type: u8) -> (usize, usize) {
    let max_block_len = if encoding_type == EncodingType::UltraBigBlock as u8 {
        BSC_MAX_BLOCK_LEN_ULTRA
    } else {
        BSC_MAX_BLOCK_LEN
    };
    let max_inflight = (DECODE_MEM_BUDGET / max_block_len).clamp(1, MAX_PARALLEL_DECOMPRESS);
    (max_block_len, max_inflight)
}

// ── Verify support ──────────────────────────────────────────────────────

/// Writer that computes CRC32 of all bytes written, without storing them.
struct HashWriter {
    crc: flate2::Crc,
}

impl HashWriter {
    fn new() -> Self {
        Self {
            crc: flate2::Crc::new(),
        }
    }

    fn checksum(&self) -> u32 {
        self.crc.sum()
    }

    fn total_bytes(&self) -> u64 {
        self.crc.amount() as u64
    }
}

impl std::io::Write for HashWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.crc.update(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Result of archive verification.
#[derive(Debug, Clone)]
pub struct VerifyResult {
    pub num_reads: usize,
    pub encoding_type: u8,
    pub header_compressor: HeaderCompressor,
    pub quality_compressor: QualityCompressor,
    pub headers_compressed_len: usize,
    pub sequences_compressed_len: usize,
    pub qualities_compressed_len: usize,
    /// CRC32 (IEEE) over the reconstructed FASTQ bytes. Only populated by deep
    /// verify; set to 0 in fast mode (use `mode` to distinguish). For paired and
    /// reference deep verify (two mates), this carries the R1-only CRC and
    /// `r2_crc32` carries R2; single-end folds all bytes into this one value.
    pub crc32: u32,
    /// In a two-mate (reference/paired) deep verify, the CRC32 over the R2 FASTQ
    /// bytes (with `crc32` holding R1's). `None` for single-end and for fast
    /// verify, which does not separate the per-mate reconstruction.
    pub r2_crc32: Option<u32>,
    /// In deep mode: total bytes of reconstructed FASTQ. In fast mode: total
    /// compressed bytes whose CRC32s were verified.
    pub total_bytes: u64,
    /// Number of per-block frames whose per-block CRC32 was verified. Only set in
    /// fast mode (zero in deep mode, since deep verify catches corruption via
    /// the codecs themselves while reconstructing records).
    pub blocks_verified: u32,
    pub mode: VerifyMode,
    pub elapsed_secs: f64,
}

/// Verification mode reported back to the caller.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VerifyMode {
    /// Full decompress + reconstruct + CRC32 over FASTQ output bytes.
    Deep,
    /// Per-block CRC32 verification only — does not invoke BSC/codec decoders.
    Fast,
}

/// Validate that `input` is a v5 chunk-major archive the bounded-streaming core
/// can decode. Returns `Ok(())` for a supported archive; errors otherwise.
///
/// qz now only ever writes v5; legacy v4 (stream-major) archives are rejected
/// with a recompress hint (there are no field archives to read back).
fn require_streamable_v5(input: &std::path::Path) -> Result<()> {
    use std::io::Read;

    // v5 chunk-major: read prefix + 12-byte body + optional 8-byte const tail.
    let mut probe = [0u8; super::V2_PREFIX_SIZE + super::V5_FIXED_BODY_SIZE + 8];
    let mut f = std::fs::File::open(input)?;
    let n = f.read(&mut probe)?;
    if n < super::V2_PREFIX_SIZE
        || super::read_archive_version(&probe[..n])? != super::ARCHIVE_VERSION_CHUNK_MAJOR
    {
        anyhow::bail!(
            "v4/legacy archives are no longer supported — recompress with current qz \
             (this build only reads v5 chunk-major archives)."
        );
    }
    // Strict parse: validates header_size, rejects unknown flags, requires
    // the const tail when flagged.
    let h = super::FixedHeader::parse_v5(&probe[..n])?;
    let streamable = crate::compression::EncodingType::from_u8(h.encoding_type)
        .is_some_and(|e| e.is_bounded_streamable());
    // Codecs the streaming core can actually decode (unified codec-id namespace).
    use crate::compression::codec_ids::{CODEC_BSC, CODEC_COLUMNAR, CODEC_FQZCOMP};
    let header_ok = matches!(h.header_compressor_code, CODEC_BSC | CODEC_COLUMNAR);
    let seq_ok = h.sequence_compressor_code == CODEC_BSC;
    let qual_ok = matches!(h.quality_compressor_code, CODEC_BSC | CODEC_FQZCOMP);
    if !(streamable && header_ok && seq_ok && qual_ok) {
        anyhow::bail!(
            "unsupported v5 archive (encoding_type={}, codecs h={} s={} q={})",
            h.encoding_type,
            h.header_compressor_code,
            h.sequence_compressor_code,
            h.quality_compressor_code
        );
    }
    Ok(())
}

/// Build the FASTQ output sink (parallel-gzip, stdout, or atomic-temp file),
/// run `decode` against it, then commit (flush/finish + atomic rename) and emit
/// the progress logging.
///
/// Memory: bounded writers keep at most ~3 decompressed blocks per active
/// stream resident at peak (BOUNDED_CHANNEL_CAP=2 + one block being consumed).
/// Fqzcomp quality is bounded-streamable (record-capped multi-block, decoded via
/// `StreamCodec::Fqz`).
///
/// This owns only the writer/atomic-commit machinery; the actual decode is the
/// supplied closure (the v5 chunk-major decode body).
pub(crate) fn decompress_streaming_bsc_dispatch(
    args: &DecompressConfig,
    start_time: Instant,
    num_reads: usize,
    decode: impl FnOnce(&mut dyn std::io::Write) -> Result<()>,
) -> Result<()> {
    use std::io::Write;
    if args.output.is_empty() {
        anyhow::bail!("No output file specified");
    }
    let output_path = &args.output[0];
    let is_stdout = crate::cli::is_stdio_path(output_path);

    // Atomic output: stream to a sibling `.tmp` and rename on success (see the
    // identical comment in `decompress_streaming_bsc`). Stdout bypasses this.
    struct TmpCleanup(std::path::PathBuf);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    let tmp_output = if is_stdout {
        None
    } else {
        let mut name = output_path
            .file_name()
            .ok_or_else(|| anyhow::anyhow!("output path has no file name: {:?}", output_path))?
            .to_os_string();
        name.push(atomic_tmp_suffix());
        Some(output_path.with_file_name(name))
    };
    let tmp_guard = tmp_output.as_ref().map(|p| TmpCleanup(p.clone()));
    let write_path = tmp_output.as_deref().unwrap_or(output_path);

    info!(
        "Decompressing {} records with bounded streaming BSC...",
        num_reads
    );

    if args.gzipped {
        use gzp::Compression;
        use gzp::ZWriter;
        use gzp::deflate::Gzip;
        use gzp::par::compress::{ParCompress, ParCompressBuilder};

        let num_gz_threads = (args.num_threads / 2).max(2);
        let mut output: ParCompress<Gzip> = if is_stdout {
            ParCompressBuilder::new()
                .num_threads(num_gz_threads)
                .map_err(|e| anyhow::anyhow!("gzp error: {e}"))?
                .compression_level(Compression::new(args.gzip_level))
                .from_writer(std::io::stdout())
        } else {
            let file = create_new_for_write(write_path)
                .with_context(|| format!("Failed to create output file: {:?}", write_path))?;
            ParCompressBuilder::new()
                .num_threads(num_gz_threads)
                .map_err(|e| anyhow::anyhow!("gzp error: {e}"))?
                .compression_level(Compression::new(args.gzip_level))
                .from_writer(file)
        };
        decode(&mut output)?;
        output
            .finish()
            .map_err(|e| anyhow::anyhow!("gzp finish error: {e}"))?;
    } else if is_stdout {
        let mut output =
            std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, std::io::stdout().lock());
        decode(&mut output)?;
        output.flush()?;
    } else {
        let file = create_new_for_write(write_path)
            .with_context(|| format!("Failed to create output file: {:?}", write_path))?;
        let mut output = std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, file);
        decode(&mut output)?;
        output.flush()?;
    }

    if let (Some(tmp), Some(guard)) = (tmp_output.as_ref(), tmp_guard) {
        std::mem::forget(guard);
        if let Err(e) = std::fs::rename(tmp, output_path) {
            let _ = std::fs::remove_file(tmp);
            anyhow::bail!("Failed to rename temp file to {:?}: {}", output_path, e);
        }
    }

    info!("Decompressed {} records", num_reads);
    let elapsed = start_time.elapsed();
    info!("Decompression completed in {:.2}s", elapsed.as_secs_f64());
    Ok(())
}

pub(super) fn decompress_streaming_bsc(args: &DecompressConfig) -> Result<()> {
    let start_time = Instant::now();

    info!("Input file: {:?}", args.input);
    info!("Output files: {:?}", args.output);
    info!("Streaming decompression mode (parallel BSC)");

    // qz only writes v5 chunk-major archives. v4/legacy archives are rejected
    // upstream (`require_streamable_v5`); guard here too so this entry point is
    // safe if called directly.
    if !is_chunk_major_v5(&args.input)? {
        anyhow::bail!(
            "v4/legacy archives are no longer supported — recompress with current qz \
             (this build only reads v5 chunk-major archives)."
        );
    }

    // v5 carries no num_reads in the front header (it lives in the footer);
    // peek it from the footer for the progress log only. The authoritative
    // count + validation happen during decode in `build_indices_from_footer`.
    let num_reads = peek_v5_num_reads(&args.input) as usize;
    decompress_streaming_bsc_dispatch(args, start_time, num_reads, |output| {
        decompress_v5_to_writer(&args.input, output)
    })
}

/// Glue between `decompress_streaming_bsc`'s raw-u8 codec codes and the bounded
/// writer (`bounded_write_records_independent`). Builds per-stream indices,
/// validates them, picks `StreamCodec` for each stream, and dispatches.
///
/// `quality_compressor_code` / `header_compressor_code` are the raw bytes from
/// the archive header body. `require_streamable_v5` has already narrowed
/// these to a small set:
///
/// - header_compressor ∈ {1=BSC, 3=Columnar}
/// - sequence_compressor = 1 (BSC)
/// - quality_compressor ∈ {1=BSC, 6=Fqzcomp} when qualities_len > 0
///
/// Streaming-decompress core: consumes prepared per-stream indices.
///
/// This is the layout-agnostic engine behind bounded streaming decompress. It
/// takes already-built `StreamIndex` values (the producer/cursor machinery
/// reads each block by its absolute `block.byte_offset`). Index construction +
/// cross-stream `validate_indices` live in the v5 footer-dispatch caller.
#[allow(clippy::too_many_arguments)]
fn decompress_streaming_with_indices(
    output: &mut dyn std::io::Write,
    archive_path: &std::path::Path,
    headers_idx: StreamIndex,
    sequences_idx: StreamIndex,
    qualities_idx: Option<StreamIndex>,
    rc_idx: Option<StreamIndex>,
    num_reads: u64,
    has_quality: bool,
    const_seq_len: usize,
    const_qual_len: usize,
    has_sequence_hints: bool,
    bits_per_qual: usize,
    quality_binning: QualityBinning,
    has_rc_canon: bool,
    header_compressor_code: u8,
    quality_compressor_code: u8,
    is_fasta: bool,
    max_block_len: usize,
    max_inflight: usize,
) -> Result<()> {
    let _ = has_quality;
    let _ = has_rc_canon;
    // Map raw codec codes to StreamCodec variants. `require_streamable_v5`
    // already filtered these down to BSC + Columnar for headers, BSC-only for
    // sequence, and {BSC, Fqzcomp} for quality. Any unrecognized code here is a
    // programmer error (gate-vs-dispatch desync) — bail with a diagnostic.
    let header_codec = match header_compressor_code {
        crate::compression::codec_ids::CODEC_BSC => StreamCodec::Bsc,
        crate::compression::codec_ids::CODEC_COLUMNAR => StreamCodec::Columnar,
        c => anyhow::bail!(
            "bounded writer: unexpected header_compressor code {c} \
             (gate should restrict to BSC=1 or Columnar=3)"
        ),
    };
    // Sequence is BSC-only under the current gate.
    let sequence_codec = StreamCodec::Bsc;

    // Quality codec: BSC (code 1) or Fqzcomp (code 6) reach this independent
    // writer. Fqzcomp decodes to FINAL raw-ASCII quality bytes framed as
    // `[varint(len)][len bytes]` per record. Force `bits_per_qual = 8` so the
    // cursor's varint read pulls all `len` bytes (the header's logical binning
    // may be < 8 bits); the writer then emits them via `QualityInput::Decoded`
    // (no unpacking) — see `is_fqz_quality` in the writer. `quality_binning` is
    // unused on the fqz path.
    let is_fqz = qualities_idx.is_some()
        && quality_compressor_code == crate::compression::codec_ids::CODEC_FQZCOMP;
    let quality_codec = if qualities_idx.is_some() {
        Some(if is_fqz {
            StreamCodec::Fqz
        } else {
            StreamCodec::Bsc
        })
    } else {
        None
    };
    let eff_bits_per_qual = if is_fqz { 8 } else { bits_per_qual };

    bounded_write_records_independent(
        output,
        num_reads,
        archive_path,
        &headers_idx,
        &sequences_idx,
        qualities_idx.as_ref(),
        rc_idx.as_ref(),
        header_codec,
        sequence_codec,
        quality_codec,
        const_seq_len,
        const_qual_len,
        has_sequence_hints,
        eff_bits_per_qual,
        quality_binning,
        is_fasta,
        max_inflight,
        max_block_len,
    )
}

/// Cap the `Vec::with_capacity` hint for offset scans against the actual
/// decompressed stream size. Each record consumes at least 1 byte (its varint
/// length prefix), so `count > data.len()` is impossible for legitimate
/// streams. Without this clamp a malformed `num_reads = 11 billion` (the
/// upper bound left by the parse-time cap) would trigger a 176 GB Vec
/// pre-allocation here.
///
/// Exposed `pub(super)` so codec entry points (`codecs.rs`, `header_col.rs`)
/// can apply the same clamp without duplicating the invariant.
#[inline]
pub(super) fn scan_capacity(count: usize, data_len: usize) -> usize {
    count.min(data_len)
}

pub(super) fn decompress(args: &DecompressConfig) -> Result<()> {
    // Honor -t. Idempotent: top-level `decompress()` already built the global
    // pool before dispatch; this also covers direct callers. See
    // `super::ensure_global_thread_pool`.
    let actual_threads = super::ensure_global_thread_pool(args.num_threads);
    info!("Using {} threads for decompression", actual_threads);

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

    // qz only writes v5 chunk-major archives, all of which are bounded-streamable.
    // Reject anything else (legacy v4 / corrupt) with a recompress hint.
    require_streamable_v5(&args.input)?;
    decompress_streaming_bsc(args)
}

/// The single-end front-header metadata the bounded streaming writer consumes,
/// resolved from the v5 fixed header. Factored out so the front-header parse
/// (via [`FixedHeader::parse_v5`]) lives in ONE place: both the full decode
/// (`decompress_v5_to_writer`) and the NUMA range decode
/// (`decode_chunk_range_single`) read it from [`read_single_end_stream_meta`].
///
/// `header_codec`/`sequence_codec`/`quality_codec`/`bits_per_qual` are resolved
/// EXACTLY as `decompress_streaming_with_indices` does (raw codec code →
/// `StreamCodec`; fqz forces 8 bits) so a meta-driven decode is byte-identical
/// to the raw-code-driven decode. `quality_codec` is the "as-if-present"
/// resolution; the writer gates it on `qualities_idx.is_some()`, so a
/// quality-less archive still decodes correctly (the field is ignored).
///
/// `max_block_len`/`max_inflight` are NOT stored — call sites derive them via
/// `decode_block_bounds(meta.encoding_type)`.
pub(super) struct SingleEndStreamMeta {
    pub(super) header_end: u64,
    pub(super) encoding_type: u8,
    pub(super) has_rc_canon: bool,
    pub(super) header_codec: StreamCodec,
    pub(super) sequence_codec: StreamCodec,
    pub(super) quality_codec: Option<StreamCodec>,
    pub(super) const_seq_len: usize,
    pub(super) const_qual_len: usize,
    pub(super) has_sequence_hints: bool,
    pub(super) bits_per_qual: usize,
    pub(super) quality_binning: QualityBinning,
    pub(super) is_fasta: bool,
}

/// Parse the v5 single-end front header exactly as `decompress_v5_to_writer`
/// does (via [`FixedHeader::parse_v5`]; `header_end` from prefix bytes [4..8]),
/// resolving the codec/binning fields the bounded writer needs. This is the
/// single front-header parse site for the single-end decode paths.
pub(super) fn read_single_end_stream_meta(
    archive: &std::path::Path,
) -> Result<SingleEndStreamMeta> {
    use std::io::Read;
    // Read the full v5 header (prefix + 12 body + optional 8 tail).
    let mut hbuf = [0u8; super::V2_PREFIX_SIZE + super::V5_FIXED_BODY_SIZE + 8];
    let mut f = std::fs::File::open(archive)?;
    let n = f.read(&mut hbuf)?;
    let h = super::FixedHeader::parse_v5(&hbuf[..n])?;
    if h.archive_type != 0 {
        anyhow::bail!(
            "not a single-end v5 archive (archive_type {})",
            h.archive_type
        );
    }
    let header_end = read_le_u32(&hbuf, 4)? as u64;
    let quality_binning = code_to_binning(h.quality_binning_code)?;

    // Resolve codecs the same way `decompress_streaming_with_indices` does:
    // header ∈ {BSC, Columnar}; sequence is BSC-only under the current gate;
    // quality is Fqz (code 6) or BSC, with fqz forcing 8 bits/qual.
    let header_codec = match h.header_compressor_code {
        crate::compression::codec_ids::CODEC_BSC => StreamCodec::Bsc,
        crate::compression::codec_ids::CODEC_COLUMNAR => StreamCodec::Columnar,
        c => anyhow::bail!(
            "bounded writer: unexpected header_compressor code {c} \
             (gate should restrict to BSC=1 or Columnar=3)"
        ),
    };
    let sequence_codec = StreamCodec::Bsc;
    let is_fqz = h.quality_compressor_code == crate::compression::codec_ids::CODEC_FQZCOMP;
    let quality_codec = Some(if is_fqz {
        StreamCodec::Fqz
    } else {
        StreamCodec::Bsc
    });
    let bits_per_qual = if is_fqz {
        8
    } else {
        quality_binning.bits_per_quality()
    };

    Ok(SingleEndStreamMeta {
        header_end,
        encoding_type: h.encoding_type,
        has_rc_canon: h.has_rc_canon(),
        header_codec,
        sequence_codec,
        quality_codec,
        const_seq_len: h.const_seq_len as usize,
        const_qual_len: h.const_qual_len as usize,
        has_sequence_hints: h.has_sequence_hints(),
        bits_per_qual,
        quality_binning,
        is_fasta: h.is_fasta(),
    })
}

/// Read a v5 archive's front header and the footer anchor (`header_end`), mode-
/// agnostically. Parses the v5 front header (prefix + 12-byte body + optional
/// 8-byte const tail) via [`FixedHeader::parse_v5`] and derives `header_end` from
/// prefix bytes [4..8] — exactly as the single-end/paired/reference front-header
/// reads do. Returns the raw [`FixedHeader`] (NO archive_type assertion; callers
/// that need a specific mode dispatch on `hdr.archive_type` themselves) plus
/// `header_end` (the byte offset the footer reader scans from). This is the paired
/// analogue of [`read_single_end_stream_meta`], kept general so the NUMA paired
/// branch can reuse the unchanged paired streaming tail.
pub(super) fn read_fixed_header(archive: &std::path::Path) -> Result<(super::FixedHeader, u64)> {
    use std::io::Read;
    let mut hbuf = [0u8; super::V2_PREFIX_SIZE + super::V5_FIXED_BODY_SIZE + 8];
    let mut f = std::fs::File::open(archive)?;
    let n = f.read(&mut hbuf)?;
    if n < super::V2_PREFIX_SIZE {
        anyhow::bail!("v5 header truncated ({n} bytes)");
    }
    let hdr = super::FixedHeader::parse_v5(&hbuf[..n])?;
    let header_end = read_le_u32(&hbuf, 4)? as u64;
    Ok((hdr, header_end))
}

/// Decode a v5 chunk-major archive into `output`.
///
/// This is the single implementation of the v5 decode body, shared by both
/// `decompress_to_writer_impl` (writes into a caller-provided sink) and
/// `decompress_streaming_bsc` (writes into the file/stdout writer it creates).
/// The caller is responsible only for opening the output sink; everything v5 —
/// re-reading the front header, footer-derived index construction
/// (`build_indices_from_footer`), cross-stream validation, and the streaming
/// decode — lives here so there is exactly one copy.
///
/// Precondition: `archive_path` is a version-5 archive whose codecs/encoding
/// have already been accepted by `require_streamable_v5`. The front header is
/// re-parsed here (cheap, strict) via `read_single_end_stream_meta` so the body
/// is self-contained.
fn decompress_v5_to_writer(
    archive_path: &std::path::Path,
    output: &mut dyn std::io::Write,
) -> Result<()> {
    let meta = read_single_end_stream_meta(archive_path)?;
    let (max_block_len, max_inflight) = decode_block_bounds(meta.encoding_type);
    let plan =
        build_indices_from_footer(archive_path, meta.header_end, max_block_len, meta.has_rc_canon)?;
    let has_quality = plan.qualities.is_some();
    validate_indices(
        plan.num_reads,
        &[
            ("headers", Some(&plan.headers)),
            ("sequences", Some(&plan.sequences)),
            ("qualities", plan.qualities.as_ref()),
            ("rc_flags", plan.rc_flags.as_ref()),
        ],
    )?;
    decompress_streaming_with_indices(
        output,
        archive_path,
        plan.headers,
        plan.sequences,
        plan.qualities,
        plan.rc_flags,
        plan.num_reads,
        has_quality,
        meta.const_seq_len,
        meta.const_qual_len,
        meta.has_sequence_hints,
        meta.quality_binning.bits_per_quality(),
        meta.quality_binning,
        meta.has_rc_canon,
        // `decompress_streaming_with_indices` re-resolves codecs from these raw
        // codes (and the plan's quality presence), so feed it the raw header
        // codes — identical to the pre-refactor call.
        match meta.header_codec {
            StreamCodec::Bsc => crate::compression::codec_ids::CODEC_BSC,
            StreamCodec::Columnar => crate::compression::codec_ids::CODEC_COLUMNAR,
            StreamCodec::Fqz => unreachable!("header codec is never Fqz"),
        },
        match meta.quality_codec {
            Some(StreamCodec::Fqz) => crate::compression::codec_ids::CODEC_FQZCOMP,
            _ => crate::compression::codec_ids::CODEC_BSC,
        },
        meta.is_fasta,
        max_block_len,
        max_inflight,
    )
}

/// True if `archive_path` is a version-5 chunk-major archive. Reads only the
/// 8-byte v2 prefix. Used to route between the v5 decode and the v4 path.
fn is_chunk_major_v5(archive_path: &std::path::Path) -> Result<bool> {
    use std::io::Read;
    let mut f0 = std::fs::File::open(archive_path)?;
    let mut prefix = [0u8; super::V2_PREFIX_SIZE];
    f0.read_exact(&mut prefix)?;
    Ok(super::read_archive_version(&prefix)? == super::ARCHIVE_VERSION_CHUNK_MAJOR)
}

/// Best-effort peek of a v5 archive's read count from the first footer field,
/// for progress logging only. Returns 0 on any inconsistency — the
/// authoritative count and full validation happen during decode in
/// `build_indices_from_footer`, so this never affects correctness.
fn peek_v5_num_reads(archive_path: &std::path::Path) -> u64 {
    use crate::compression::chunk_directory::LOCATOR_LEN;
    use std::io::{Read, Seek, SeekFrom};
    (|| -> Result<u64> {
        let file_size = std::fs::metadata(archive_path)?.len();
        if file_size < LOCATOR_LEN as u64 {
            return Ok(0);
        }
        let mut f = std::fs::File::open(archive_path)?;
        f.seek(SeekFrom::Start(file_size - LOCATOR_LEN as u64))?;
        let mut loc = [0u8; LOCATOR_LEN];
        f.read_exact(&mut loc)?;
        let footer_len = u64::from_le_bytes(loc[0..8].try_into().unwrap());
        // footer body starts with num_reads (u64). Bounds-check before seeking.
        let Some(footer_start) = (file_size - LOCATOR_LEN as u64).checked_sub(footer_len) else {
            return Ok(0);
        };
        if footer_start
            .checked_add(8)
            .is_none_or(|end| end > file_size)
        {
            return Ok(0);
        }
        f.seek(SeekFrom::Start(footer_start))?;
        let mut nb = [0u8; 8];
        f.read_exact(&mut nb)?;
        Ok(u64::from_le_bytes(nb))
    })()
    .unwrap_or(0)
}

/// Decompress an archive directly into `output`, routing the bounded streaming
/// path (for eligible v4 archives) or the mmap/fallback path.
///
/// This is the internal implementation for the public
/// `compression::decompress_to_writer` API. It must not create its own output
/// file — it writes everything to the provided writer.
///
/// **BrokenPipe**: the bounded writers already convert `ErrorKind::BrokenPipe`
/// to `Ok(())` via `write_or_broken_pipe`. For the mmap/fallback path this
/// function also catches `BrokenPipe` at the top level and returns `Ok(())`.
pub(super) fn decompress_to_writer_impl(
    archive_path: &std::path::Path,
    output: &mut dyn std::io::Write,
) -> Result<()> {
    // qz only writes v5 chunk-major archives, decoded by the single shared v5
    // body. v4/legacy or otherwise-unsupported archives are rejected with a
    // recompress hint.
    require_streamable_v5(archive_path)?;
    decompress_v5_to_writer(archive_path, output)
}

// ── Verify ──────────────────────────────────────────────────────────────

/// Parse the streaming header and write records to `hasher`.
///
/// Routes through `decompress_streaming_with_indices` so the bounded code path
/// (bounded channels, per-block cursor, robustness fixes) is fully exercised
/// during deep verify. Only v5 chunk-major archives are supported; legacy v4
/// archives are rejected with a recompress hint.
fn verify_streaming(input: &std::path::Path, hasher: &mut HashWriter) -> Result<ArchiveHeader> {
    use std::io::Read;

    let mut file = std::fs::File::open(input)
        .with_context(|| format!("Failed to open archive: {:?}", input))?;
    let mut header = [0u8; 68];
    file.read_exact(&mut header[..60])
        .context("Failed to read archive header")?;

    // ── v5 chunk-major: reconstruct FASTQ into the hasher ────────────────
    // Deep verify reconstructs the full FASTQ into the HashWriter (the sink),
    // exercising codec decode + record assembly — NOT just block CRCs. This is
    // the same body as `decompress_v5_to_writer`, but the output is the hasher.
    if super::read_archive_version(&header[..super::V2_PREFIX_SIZE])?
        == super::ARCHIVE_VERSION_CHUNK_MAJOR
    {
        let h = super::FixedHeader::parse_v5(&header[..60])?;
        // No `archive_type` guard needed here: the public `verify` dispatch routes
        // a paired (type 1) archive to `verify_paired_v5` before reaching this
        // single-end path, and `build_indices_from_footer` rejects any `mate != 0`
        // entry as a backstop, so only single-end (type 0) frames arrive here.
        let header_end = read_le_u32(&header, 4)? as u64;
        let (max_block_len, max_inflight) = decode_block_bounds(h.encoding_type);
        let plan = build_indices_from_footer(input, header_end, max_block_len, h.has_rc_canon())?;
        let has_quality = plan.qualities.is_some();
        validate_indices(
            plan.num_reads,
            &[
                ("headers", Some(&plan.headers)),
                ("sequences", Some(&plan.sequences)),
                ("qualities", plan.qualities.as_ref()),
                ("rc_flags", plan.rc_flags.as_ref()),
            ],
        )?;
        let quality_binning = code_to_binning(h.quality_binning_code)?;
        // Compressed per-role byte totals for the VerifyResult metadata (v5 has no
        // contiguous stream regions — sum each role's block payload bytes).
        let headers_bytes: u64 = plan.headers.blocks.iter().map(|b| b.byte_len as u64).sum();
        let sequences_bytes: u64 = plan
            .sequences
            .blocks
            .iter()
            .map(|b| b.byte_len as u64)
            .sum();
        let qualities_bytes: u64 = plan
            .qualities
            .as_ref()
            .map(|idx| idx.blocks.iter().map(|b| b.byte_len as u64).sum())
            .unwrap_or(0);
        let num_reads = plan.num_reads as usize;
        let header_compressor = codec_ids::codec_id_to_header_compressor(h.header_compressor_code)?;
        let quality_compressor = codec_ids::codec_id_to_quality_compressor(h.quality_compressor_code)?;
        let is_fasta = h.is_fasta();
        info!(
            "Verifying {} reads (v5 chunk-major streaming)...",
            num_reads
        );
        decompress_streaming_with_indices(
            hasher,
            input,
            plan.headers,
            plan.sequences,
            plan.qualities,
            plan.rc_flags,
            plan.num_reads,
            has_quality,
            h.const_seq_len as usize,
            h.const_qual_len as usize,
            h.has_sequence_hints(),
            quality_binning.bits_per_quality(),
            quality_binning,
            h.has_rc_canon(),
            h.header_compressor_code,
            h.quality_compressor_code,
            is_fasta,
            max_block_len,
            max_inflight,
        )?;
        return Ok(ArchiveHeader {
            encoding_type: h.encoding_type,
            quality_compressor,
            header_compressor,
            num_reads,
            headers_len: headers_bytes as usize,
            sequences_len: sequences_bytes as usize,
            qualities_len: qualities_bytes as usize,
        });
    }

    anyhow::bail!(
        "v4/legacy archives are no longer supported — recompress with current qz \
         (this build only reads v5 chunk-major archives)."
    )
}

/// Resolve a CLI thread count: 0 means "auto-detect", otherwise honor it.
fn resolve_num_threads(requested: usize) -> usize {
    if requested == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(8)
    } else {
        requested
    }
}

/// Verify a QZ archive: fully decompress all streams, compute CRC32, report metadata.
pub(super) fn verify(config: &VerifyConfig) -> Result<VerifyResult> {
    let num_threads = resolve_num_threads(config.num_threads);
    // See compress() and decompress() for the rayon global-pool caveat.
    let actual_threads = match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
    {
        Ok(()) => num_threads,
        Err(_) => {
            let existing = rayon::current_num_threads();
            if existing != num_threads {
                tracing::warn!(
                    "Requested {} threads but the rayon pool was already \
                     initialized to {} threads (process-wide setting). \
                     Using {} threads.",
                    num_threads,
                    existing,
                    existing
                );
            }
            existing
        }
    };
    info!("Using {} threads for verification", actual_threads);

    let start_time = Instant::now();

    // Handle stdin: spool to temp file
    let (_stdin_tmp, input_path) = if crate::cli::is_stdio_path(&config.input) {
        info!("Reading archive from stdin...");
        let mut tmp = tempfile::NamedTempFile::new_in(&config.working_dir)
            .context("Failed to create temp file for stdin")?;
        std::io::copy(&mut std::io::stdin().lock(), &mut tmp)
            .context("Failed to spool stdin to temp file")?;
        let path = tmp.path().to_path_buf();
        (Some(tmp), std::borrow::Cow::Owned(path))
    } else {
        (None, std::borrow::Cow::Borrowed(&config.input))
    };
    let input_path: &std::path::Path = &input_path;

    if config.fast {
        return verify_fast(input_path, start_time);
    }

    let mut hasher = HashWriter::new();

    // qz only writes v5 chunk-major archives; reject anything else up front.
    require_streamable_v5(input_path)?;
    info!("Verifying (streaming mode)...");
    let hdr = verify_streaming(input_path, &mut hasher)?;

    let elapsed = start_time.elapsed();
    info!("Verification completed in {:.2}s", elapsed.as_secs_f64());

    Ok(VerifyResult {
        num_reads: hdr.num_reads,
        encoding_type: hdr.encoding_type,
        header_compressor: hdr.header_compressor,
        quality_compressor: hdr.quality_compressor,
        headers_compressed_len: hdr.headers_len,
        sequences_compressed_len: hdr.sequences_len,
        qualities_compressed_len: hdr.qualities_len,
        crc32: hasher.checksum(),
        r2_crc32: None, // single-end: one stream, no per-mate split
        total_bytes: hasher.total_bytes(),
        blocks_verified: 0,
        mode: VerifyMode::Deep,
        elapsed_secs: elapsed.as_secs_f64(),
    })
}

/// Fast verify: walk every block frame in each stream and verify its
/// CRC32 without invoking the inner codec. Catches bit-rot in O(IO + CRC)
/// time. Does NOT reconstruct FASTQ output — the deep CRC32 over decoded
/// bytes is not computed (caller can run plain `verify` for that).
///
/// Only v5 chunk-major archives are supported; legacy v4 archives are rejected
/// with a recompress hint. Archives with encoding_type 10 (`UltraBigBlock`) are
/// fully supported: they use the same per-block CRC framing as default archives.
fn verify_fast(input_path: &std::path::Path, start_time: Instant) -> Result<VerifyResult> {
    use crate::compression::bsc;

    info!("Verifying (fast mode: per-block CRC32 only)...");

    // ── v5 chunk-major: per-block inline CRC walk (no reconstruction) ────
    // v5 has no contiguous stream regions; fast verify deliberately does NOT
    // reconstruct — it builds the decode plan from the footer and CRC-checks
    // every role index's blocks.
    {
        use std::io::Read;
        let mut probe = std::fs::File::open(input_path)?;
        let mut hbuf = [0u8; super::V2_PREFIX_SIZE + super::V5_FIXED_BODY_SIZE + 8];
        let n = probe.read(&mut hbuf)?;
        if n >= super::V2_PREFIX_SIZE
            && super::read_archive_version(&hbuf[..super::V2_PREFIX_SIZE])?
                == super::ARCHIVE_VERSION_CHUNK_MAJOR
        {
            let h = super::FixedHeader::parse_v5(&hbuf[..n])?;
            // No `archive_type` guard needed here: the public `verify` dispatch routes
            // a paired (type 1) archive to `verify_paired_v5` before reaching this
            // single-end path, and `build_indices_from_footer` rejects any `mate != 0`
            // entry as a backstop, so only single-end (type 0) frames arrive here.
            let header_end = read_le_u32(&hbuf, 4)? as u64;
            let (max_block_len, _max_inflight) = decode_block_bounds(h.encoding_type);
            let plan =
                build_indices_from_footer(input_path, header_end, max_block_len, h.has_rc_canon())?;
            validate_indices(
                plan.num_reads,
                &[
                    ("headers", Some(&plan.headers)),
                    ("sequences", Some(&plan.sequences)),
                    ("qualities", plan.qualities.as_ref()),
                    ("rc_flags", plan.rc_flags.as_ref()),
                ],
            )?;

            fn crc_walk(
                file: &mut std::fs::File,
                idx: &StreamIndex,
                blocks: &mut usize,
                bytes: &mut u64,
            ) -> Result<()> {
                use std::io::{Read, Seek, SeekFrom};
                for b in &idx.blocks {
                    file.seek(SeekFrom::Start(b.byte_offset))?;
                    let mut payload = vec![0u8; b.byte_len as usize];
                    file.read_exact(&mut payload)?;
                    bsc::verify_block_frame_crc(b.expected_crc, b.record_count, &payload)?;
                    *blocks += 1;
                    *bytes += b.byte_len as u64;
                }
                Ok(())
            }

            let mut file = std::fs::File::open(input_path)?;
            let mut total_blocks: usize = 0;
            let mut total_compressed: u64 = 0;
            crc_walk(
                &mut file,
                &plan.headers,
                &mut total_blocks,
                &mut total_compressed,
            )?;
            crc_walk(
                &mut file,
                &plan.sequences,
                &mut total_blocks,
                &mut total_compressed,
            )?;
            if let Some(idx) = plan.qualities.as_ref() {
                crc_walk(&mut file, idx, &mut total_blocks, &mut total_compressed)?;
            }
            if let Some(idx) = plan.rc_flags.as_ref() {
                crc_walk(&mut file, idx, &mut total_blocks, &mut total_compressed)?;
            }

            // The whole-archive `ChunkDecodedSizes` global block (chunk_index ==
            // GLOBAL_SENTINEL, role 18) lives OUTSIDE the DecodePlan's per-role indices,
            // so the role-index walk above never touches it. Walk its frame too — without
            // this, a bit-flip in the table payload slips through `verify --fast` (it only
            // ever fails the decode-time table read in `read_chunk_layout`). The block is
            // OPTIONAL (archives predating the table omit it), so only walk it if present.
            {
                use crate::compression::chunk_directory::{GLOBAL_SENTINEL, StreamRole, read_v5_footer};
                use std::io::{Read, Seek, SeekFrom};
                let dir = read_v5_footer(input_path, header_end)?;
                if let Some(e) = dir.entries.iter().find(|e| {
                    e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::ChunkDecodedSizes
                }) {
                    if e.length < 12 {
                        anyhow::bail!(
                            "ChunkDecodedSizes frame too short ({} bytes, need >= 12)",
                            e.length
                        );
                    }
                    file.seek(SeekFrom::Start(e.offset))?;
                    let mut frame = vec![0u8; e.length as usize];
                    file.read_exact(&mut frame)?;
                    let block_len = u32::from_le_bytes(frame[0..4].try_into().unwrap()) as usize;
                    let record_count = u32::from_le_bytes(frame[4..8].try_into().unwrap());
                    let crc = u32::from_le_bytes(frame[8..12].try_into().unwrap());
                    let payload = &frame[12..];
                    if payload.len() != block_len {
                        anyhow::bail!(
                            "ChunkDecodedSizes frame: block_len field {block_len} != actual payload bytes {}",
                            payload.len()
                        );
                    }
                    bsc::verify_block_frame_crc(crc, record_count, payload)?;
                    total_blocks += 1;
                    total_compressed += payload.len() as u64;
                }
            }

            let elapsed = start_time.elapsed();
            info!(
                "Fast verification completed in {:.2}s: {} blocks across {} bytes (v5 chunk-major)",
                elapsed.as_secs_f64(),
                total_blocks,
                total_compressed,
            );

            return Ok(VerifyResult {
                num_reads: plan.num_reads as usize,
                encoding_type: h.encoding_type,
                header_compressor: codec_ids::codec_id_to_header_compressor(h.header_compressor_code)?,
                quality_compressor: codec_ids::codec_id_to_quality_compressor(h.quality_compressor_code)?,
                headers_compressed_len: plan
                    .headers
                    .blocks
                    .iter()
                    .map(|b| b.byte_len as usize)
                    .sum(),
                sequences_compressed_len: plan
                    .sequences
                    .blocks
                    .iter()
                    .map(|b| b.byte_len as usize)
                    .sum(),
                qualities_compressed_len: plan
                    .qualities
                    .as_ref()
                    .map(|idx| idx.blocks.iter().map(|b| b.byte_len as usize).sum())
                    .unwrap_or(0),
                crc32: 0,       // not computed in fast mode
                r2_crc32: None, // single-end
                total_bytes: total_compressed,
                blocks_verified: total_blocks as u32,
                mode: VerifyMode::Fast,
                elapsed_secs: elapsed.as_secs_f64(),
            });
        }
    }

    anyhow::bail!(
        "v4/legacy archives are no longer supported — recompress with current qz \
         (this build only reads v5 chunk-major archives)."
    )
}

/// One block's location and metadata, from the per-block frame.
#[derive(Debug, Clone)]
pub(super) struct BlockEntry {
    pub record_count: u32,
    pub byte_offset: u64, // absolute offset in the archive file/data where the payload starts
    pub byte_len: u32,    // payload length (does NOT include the 12-byte framing prefix)
    pub expected_crc: u32,
}

#[derive(Debug, Clone)]
pub(super) struct StreamIndex {
    pub blocks: Vec<BlockEntry>,
    pub total_records: u64,
}

/// Cross-check that every *active* stream's `total_records` equals
/// `archive_num_reads`. Inactive streams pass `None`. The stream name is
/// preserved in error messages for diagnostic clarity.
///
/// Used by the bounded streaming decoder to catch malformed archives at
/// index-build time, before any payload decompression.
pub(super) fn validate_indices(
    archive_num_reads: u64,
    streams: &[(&str, Option<&StreamIndex>)],
) -> Result<()> {
    for (name, idx) in streams {
        if let Some(idx) = idx
            && idx.total_records != archive_num_reads
        {
            anyhow::bail!(
                "stream length mismatch: {name} total records {} ≠ archive num_reads {}",
                idx.total_records,
                archive_num_reads
            );
        }
    }
    Ok(())
}

/// Cursor over a single stream's bounded mpsc-delivered decompressed blocks.
/// Reads one record at a time, refilling from the channel when the current
/// block is drained. Used by the bounded streaming writer to emit FASTQ
/// records without buffering full streams in memory.
pub(super) struct StreamCursor {
    name: &'static str,
    rx: std::sync::mpsc::Receiver<Result<(u32, Vec<u8>), String>>,
    current: Option<CurrentBlock>,
}

struct CurrentBlock {
    records_left: u32,
    data: Vec<u8>,
    pos: usize,
}

impl StreamCursor {
    pub(super) fn new(
        rx: std::sync::mpsc::Receiver<Result<(u32, Vec<u8>), String>>,
        name: &'static str,
    ) -> Self {
        Self {
            name,
            rx,
            current: None,
        }
    }

    /// Pull the next block from the channel, replacing `current`. Returns
    /// `Ok(true)` if a block arrived, `Ok(false)` if the channel closed
    /// cleanly (end-of-stream), or `Err` if the producer sent an error.
    fn refill(&mut self) -> Result<bool> {
        match self.rx.recv() {
            Ok(Ok((rc, data))) => {
                self.current = Some(CurrentBlock {
                    records_left: rc,
                    data,
                    pos: 0,
                });
                Ok(true)
            }
            Ok(Err(e)) => anyhow::bail!("stream {}: {}", self.name, e),
            Err(_) => Ok(false), // channel closed by producer
        }
    }

    /// Read one varint-length-prefixed record. Returns the payload bytes.
    // The `loop` is a refill scaffold: on a drained block it resets `current`
    // and refills before reading. It returns on the first successful read.
    #[allow(clippy::never_loop)]
    pub(super) fn read_one_record_varint(&mut self) -> Result<Vec<u8>> {
        loop {
            if self.current.as_ref().is_some_and(|c| c.records_left == 0) {
                self.current = None; // drained — refill below
            }
            if self.current.is_none() && !self.refill()? {
                anyhow::bail!(
                    "stream {}: unexpected end of stream while reading record",
                    self.name
                );
            }
            let block = self.current.as_mut().unwrap();
            let len = read_varint(&block.data, &mut block.pos).ok_or_else(|| {
                anyhow::anyhow!(
                    "stream {}: truncated varint in block at pos {}",
                    self.name,
                    block.pos,
                )
            })?;
            // checked_add: `len` is an untrusted varint (up to usize::MAX), so a
            // plain `block.pos + len` can wrap and slip past the bounds check.
            let end = block
                .pos
                .checked_add(len)
                .filter(|&e| e <= block.data.len())
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "stream {}: record claims {} bytes but block only has {} left",
                        self.name,
                        len,
                        block.data.len() - block.pos
                    )
                })?;
            let out = block.data[block.pos..end].to_vec();
            block.pos = end;
            block.records_left -= 1;
            return Ok(out);
        }
    }

    /// Read exactly `len` bytes from the stream as one record. Used when the
    /// archive header indicates const_seq_len / const_qual_len != 0 (no
    /// per-record varint length prefix).
    // The `loop` is a refill scaffold: on a drained block it resets `current`
    // and refills before reading. It returns on the first successful read.
    #[allow(clippy::never_loop)]
    pub(super) fn read_one_record_const(&mut self, len: usize) -> Result<Vec<u8>> {
        loop {
            if self.current.as_ref().is_some_and(|c| c.records_left == 0) {
                self.current = None;
            }
            if self.current.is_none() && !self.refill()? {
                anyhow::bail!("stream {}: EOF before const-length record", self.name);
            }
            let block = self.current.as_mut().unwrap();
            // checked_add: `len` comes from the archive header (const_*_len), so
            // guard against `block.pos + len` wrapping past the bounds check.
            let end = block
                .pos
                .checked_add(len)
                .filter(|&e| e <= block.data.len())
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "stream {}: const-length record needs {} bytes, only {} remain in block",
                        self.name,
                        len,
                        block.data.len() - block.pos
                    )
                })?;
            let out = block.data[block.pos..end].to_vec();
            block.pos = end;
            block.records_left -= 1;
            return Ok(out);
        }
    }

    /// Read one byte from the stream (used for RC flags: 1 byte per record).
    /// Same block refill semantics as `read_one_record_varint`.
    // The `loop` is a refill scaffold: on a drained block it resets `current`
    // and refills before reading. It returns on the first successful read.
    #[allow(clippy::never_loop)]
    pub(super) fn read_one_byte(&mut self) -> Result<u8> {
        loop {
            if self.current.as_ref().is_some_and(|c| c.records_left == 0) {
                self.current = None;
            }
            if self.current.is_none() && !self.refill()? {
                anyhow::bail!("stream {}: unexpected EOF reading byte", self.name);
            }
            let block = self.current.as_mut().unwrap();
            if block.pos >= block.data.len() {
                anyhow::bail!("stream {}: ran out of bytes in block", self.name);
            }
            let b = block.data[block.pos];
            block.pos += 1;
            block.records_left -= 1;
            return Ok(b);
        }
    }

    /// Read one variable-length sequence record laid out as
    /// `[varint(L)] [hint_byte?] [L bytes of seq]`. The hint byte is consumed
    /// and discarded — sequence_hints are encoder-side prefetch metadata that
    /// the writer does not need. Matches `scan_seq_offsets` in the legacy path.
    // The `loop` is a refill scaffold: on a drained block it resets `current`
    // and refills before reading. It returns on the first successful read.
    #[allow(clippy::never_loop)]
    pub(super) fn read_one_record_seq_with_hint(&mut self, has_hints: bool) -> Result<Vec<u8>> {
        loop {
            if self.current.as_ref().is_some_and(|c| c.records_left == 0) {
                self.current = None;
            }
            if self.current.is_none() && !self.refill()? {
                anyhow::bail!(
                    "stream {}: unexpected end of stream while reading seq record",
                    self.name
                );
            }
            let block = self.current.as_mut().unwrap();
            let len = read_varint(&block.data, &mut block.pos).ok_or_else(|| {
                anyhow::anyhow!(
                    "stream {}: truncated varint in seq block at pos {}",
                    self.name,
                    block.pos,
                )
            })?;
            let hint_extra = if has_hints { 1 } else { 0 };
            // checked_add throughout: `len` is an untrusted varint. Compute the
            // post-record end (incl. the optional hint byte) without wrapping.
            let end = block
                .pos
                .checked_add(hint_extra)
                .and_then(|p| p.checked_add(len))
                .filter(|&e| e <= block.data.len())
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "stream {}: seq record claims {} bytes (+{} hint) but block only has {} left",
                        self.name,
                        len,
                        hint_extra,
                        block.data.len() - block.pos,
                    )
                })?;
            if has_hints {
                // Skip the hint byte (encoder-side prefetch metadata, unused on decode).
                block.pos += 1;
            }
            let out = block.data[block.pos..end].to_vec();
            block.pos = end;
            block.records_left -= 1;
            return Ok(out);
        }
    }

    /// Read one variable-length quality record laid out as
    /// `[varint(l_seq)] [packed_bytes]` where the packed byte length is
    /// `(l_seq * bits_per_qual + 7) / 8`. Returns the *packed* bytes plus
    /// `l_seq` so the writer can call `columnar::unpack_qualities_to_writer`.
    /// Mirrors `scan_qual_offsets`.
    // The `loop` is a refill scaffold (see read_one_record_varint).
    #[allow(clippy::never_loop)]
    pub(super) fn read_one_record_qual_varint(
        &mut self,
        bits_per_qual: usize,
    ) -> Result<(usize, Vec<u8>)> {
        loop {
            if self.current.as_ref().is_some_and(|c| c.records_left == 0) {
                self.current = None;
            }
            if self.current.is_none() && !self.refill()? {
                anyhow::bail!(
                    "stream {}: unexpected end of stream while reading qual record",
                    self.name
                );
            }
            let block = self.current.as_mut().unwrap();
            let l_seq = read_varint(&block.data, &mut block.pos).ok_or_else(|| {
                anyhow::anyhow!(
                    "stream {}: truncated varint in qual block at pos {}",
                    self.name,
                    block.pos,
                )
            })?;
            let packed_len = l_seq
                .checked_mul(bits_per_qual)
                .and_then(|n| n.checked_add(7))
                .map(|n| n / 8)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "stream {}: qual length overflow l_seq={} bits_per_qual={}",
                        self.name,
                        l_seq,
                        bits_per_qual,
                    )
                })?;
            // checked_add: packed_len is bounded by the checked_mul above, but
            // `block.pos + packed_len` is still guarded against wrapping.
            let end = block
                .pos
                .checked_add(packed_len)
                .filter(|&e| e <= block.data.len())
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "stream {}: qual record needs {} packed bytes (l_seq={}), only {} remain in block",
                        self.name,
                        packed_len,
                        l_seq,
                        block.data.len() - block.pos,
                    )
                })?;
            let out = block.data[block.pos..end].to_vec();
            block.pos = end;
            block.records_left -= 1;
            return Ok((l_seq, out));
        }
    }
}

/// Maximum number of decompress tasks dispatched to rayon per producer at any
/// one time. Combined with `BOUNDED_CHANNEL_CAP`, this caps the in-flight
/// memory per stream at `(BOUNDED_CHANNEL_CAP + MAX_PARALLEL_DECOMPRESS) *
/// max_block_size`.
///
/// Tuned on a 10M-read benchmark (Phase 2): the parallel-format writer
/// absorbs producer output much faster than the original serial loop, so the
/// previous cap of 4 left the producers idle. Sweeping 4 → 8 → 16 → 32
/// showed 16 hits the sweet spot — 17.3 s wall vs 26.3 s at 4 and 17.6 s at
/// 32, with peak RSS ≈ 5.9 GB (well under the legacy in-memory path's
/// 7.6 GB). Picked empirically; bump cautiously if memory headroom changes.
/// Upper bound on the per-stream in-flight decode count. The actual cap is
/// `decode_block_bounds`-derived (`DECODE_MEM_BUDGET / max_block_len`), clamped
/// to this value so small-block (default) archives keep the historical 16.
///
/// `pub(crate)` so the reference decode reuses the same concurrency ceiling for
/// its bounded-parallel per-chunk reconstruction.
pub(crate) const MAX_PARALLEL_DECOMPRESS: usize = 16;

/// Run a pipelined parallel-decompress producer loop.
///
/// The producer thread reads compressed payloads sequentially from disk
/// (cheap I/O — 5–25 MB per block) and dispatches the per-block CRC check
/// plus `decompress` closure to a rayon thread. Results are collected via a
/// helper `mpsc` channel and forwarded to the bounded output `tx` strictly
/// in original block order, so the writer sees the same byte sequence it
/// would have seen from the old sequential producer.
///
/// Bounded-memory invariant:
///   in-flight decompress tasks ≤ MAX_PARALLEL_DECOMPRESS
///   buffered output blocks    ≤ tx channel cap (BOUNDED_CHANNEL_CAP)
/// so peak per-stream memory is `(BOUNDED_CHANNEL_CAP +
/// MAX_PARALLEL_DECOMPRESS) * max_block_size`.
///
/// Backpressure: when the writer is slow, `tx.send` blocks. While blocked
/// the producer thread stops draining `done_rx`, the in-flight slot count
/// stops decreasing, the dispatch loop stalls on the inner `done_rx.recv`,
/// and no further blocks are pulled from disk or queued onto rayon.
///
/// Errors: per-block CRC and decompress failures are reported through `tx`
/// as `Err(String)` (matching the previous behavior). I/O errors while
/// reading compressed payloads bail out of the producer via `?`. BrokenPipe
/// on `tx.send` (writer dropped) causes a graceful return.
fn run_parallel_producer<F>(
    archive_path: &std::path::Path,
    index: &StreamIndex,
    tx: std::sync::mpsc::SyncSender<Result<(u32, Vec<u8>), String>>,
    label: &'static str,
    max_inflight: usize,
    decompress: F,
) -> Result<()>
where
    F: Fn(usize, u32, Vec<u8>) -> Result<(u32, Vec<u8>), String> + Send + Sync + 'static,
{
    use std::collections::BTreeMap;
    use std::io::{Read, Seek, SeekFrom};
    use std::sync::Arc;

    let mut file = std::fs::File::open(archive_path)?;
    let decompress = Arc::new(decompress);

    // The "done" buffer mirrors the in-flight cap: it never needs to hold more
    // than `max_inflight` dispatched-but-unforwarded results.
    let done_cap = max_inflight.max(1);
    let (done_tx, done_rx) =
        std::sync::mpsc::sync_channel::<(usize, Result<(u32, Vec<u8>), String>)>(done_cap);

    let mut next_send_idx: usize = 0;
    let mut dispatched: usize = 0;
    let mut completed: BTreeMap<usize, Result<(u32, Vec<u8>), String>> = BTreeMap::new();

    // Helper: forward in-order results to `tx`. Returns `Ok(true)` if the
    // writer hung up (caller should bail early), `Ok(false)` otherwise.
    let flush_ordered = |completed: &mut BTreeMap<usize, Result<(u32, Vec<u8>), String>>,
                         next_send_idx: &mut usize|
     -> Result<bool> {
        while let Some(result) = completed.remove(next_send_idx) {
            if tx.send(result).is_err() {
                return Ok(true);
            }
            *next_send_idx += 1;
        }
        Ok(false)
    };

    for (idx, block) in index.blocks.iter().enumerate() {
        // Backpressure: drain done_rx until the in-flight count is below the
        // cap. This blocks naturally if rayon is saturated, and (transitively
        // via tx.send inside flush_ordered) if the writer is slow.
        while dispatched - next_send_idx >= max_inflight {
            let (done_idx, result) = done_rx
                .recv()
                .map_err(|e| anyhow::anyhow!("{label}: done_rx disconnected: {e}"))?;
            completed.insert(done_idx, result);
            if flush_ordered(&mut completed, &mut next_send_idx)? {
                return Ok(());
            }
        }

        // Sequential I/O: read this block's payload off disk.
        file.seek(SeekFrom::Start(block.byte_offset))?;
        let mut payload = vec![0u8; block.byte_len as usize];
        file.read_exact(&mut payload)?;
        let expected_crc = block.expected_crc;
        let record_count = block.record_count;
        let done_tx_clone = done_tx.clone();
        let decompress_cloned = Arc::clone(&decompress);
        let label_owned = label;

        rayon::spawn(move || {
            // rayon::spawn ABORTS the process if a panic escapes the closure
            // (rayon's AbortIfPanic → process::abort, which defeats panic="unwind"
            // and PyO3's panic→exception conversion). This decode runs on untrusted
            // input, so convert any panic into a normal decode Err — mirrors the
            // reference decode path (reference/decode.rs ~line 351).
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                match crate::compression::bsc::verify_block_frame_crc(
                    expected_crc,
                    record_count,
                    &payload,
                ) {
                    Ok(()) => decompress_cloned(idx, record_count, payload),
                    Err(e) => Err(format!("{label_owned} block {idx}: {e}")),
                }
            }))
            .unwrap_or_else(|_| Err(format!("{label_owned} block {idx}: decode panicked")));
            // Receiver lives until this task is reaped; ignore send error
            // (only possible if producer exited early on BrokenPipe).
            let _ = done_tx_clone.send((idx, result));
        });
        dispatched += 1;
    }

    // Dispose of our handle so a BrokenPipe early-return doesn't leave
    // straggler tasks blocked forever on send (the producer's local handle
    // is the only one we'd hold past the loop).
    drop(done_tx);

    // Drain remaining in-flight tasks in order.
    while next_send_idx < dispatched {
        let (done_idx, result) = done_rx
            .recv()
            .map_err(|e| anyhow::anyhow!("{label}: done_rx disconnected (drain): {e}"))?;
        completed.insert(done_idx, result);
        if flush_ordered(&mut completed, &mut next_send_idx)? {
            return Ok(());
        }
    }
    Ok(())
}

/// Producer thread: walk a `StreamIndex` and emit `(record_count,
/// decompressed_block)` tuples through a bounded mpsc channel.
///
/// Block decompression is dispatched to rayon (up to
/// `MAX_PARALLEL_DECOMPRESS` in flight) and results are forwarded in
/// original block order — preserving the writer's input contract while
/// recovering the within-stream parallelism the old non-bounded path
/// exploited via `bsc::decompress_parallel`.
///
/// Errors are serialized into `String`s through the channel (closest analogue
/// of `Result<_, anyhow::Error>` that's `Send + 'static` for a mpsc message).
/// If the receiver hangs up (BrokenPipe path), the producer exits cleanly via
/// `tx.send` returning `Err`.
pub(super) fn bounded_bsc_producer(
    archive_path: &std::path::Path,
    index: &StreamIndex,
    tx: std::sync::mpsc::SyncSender<Result<(u32, Vec<u8>), String>>,
    max_inflight: usize,
    max_block_len: usize,
) -> Result<()> {
    // Ultra archives (the only ones whose cap exceeds the default) carry large
    // blocks where the per-block inverse-BWT dominates and single-threaded decode
    // is the bottleneck. Use libbsc's OpenMP MT inverse-BWT for those. Inverse
    // BWT is memory-light (unlike the forward BWT), so running a few concurrently
    // under the bounded `max_inflight` pipeline is safe. Default archives keep the
    // single-threaded per-block path (many small blocks → rayon inter-block
    // parallelism already saturates the cores). Output is identical either way.
    let use_mt = max_block_len > BSC_MAX_BLOCK_LEN;
    run_parallel_producer(
        archive_path,
        index,
        tx,
        "bsc",
        max_inflight,
        move |idx, record_count, payload| {
            let decoded = if use_mt {
                crate::compression::bsc::decompress_mt_max(&payload, max_block_len)
            } else {
                crate::compression::bsc::decompress_one_block(&payload, max_block_len)
            };
            match decoded {
                Ok(d) => Ok((record_count, d)),
                Err(e) => Err(format!("block {idx} decompress: {e}")),
            }
        },
    )
}

/// Producer for record-capped multi-block fqzcomp quality streams.
///
/// Each v4 block carries one record-capped fqzcomp blob (≤ `quality_ctx_block_size`
/// records, locally mean-quality-sorted at encode). `decompress_qualities_fqzcomp`
/// decodes a block to the `[varint(len)][raw ASCII bytes]` per-record stream —
/// byte-for-byte the same layout the quality `StreamCursor` consumes for BSC
/// quality (with identity binning / `bits_per_qual = 8`), so no re-framing is
/// needed here.
///
/// DoS: `decompress_qualities_fqzcomp` already bounds its allocations
/// (`FQZ_MAX_READS_PER_BLOB`, per-read length / truncation / negative-length
/// checks). On top of that, fqz block decode transiently holds ~2× a block's
/// decoded quality (the concat + the varint stream), and a block can expand well
/// past its compressed size, so we clamp in-flight fqz block decodes below the
/// BSC ceiling to keep peak resident memory bounded.
fn bounded_fqz_producer(
    archive_path: &std::path::Path,
    index: &StreamIndex,
    tx: std::sync::mpsc::SyncSender<Result<(u32, Vec<u8>), String>>,
    max_inflight: usize,
) -> Result<()> {
    const MAX_INFLIGHT_FQZ: usize = 4;
    let mi = max_inflight.clamp(1, MAX_INFLIGHT_FQZ);
    run_parallel_producer(
        archive_path,
        index,
        tx,
        "fqz",
        mi,
        move |idx, record_count, payload| {
            match crate::compression::codecs::decompress_qualities_fqzcomp(&payload) {
                Ok(stream) => Ok((record_count, stream)),
                Err(e) => Err(format!("fqz quality block {idx} decompress: {e}")),
            }
        },
    )
}

/// Producer for v4-framed columnar header streams.
///
/// Each payload starts with `[chunk_reads: u32][columnar_blob]`. The producer
/// decodes the columnar blob into individual header strings, then re-encodes
/// them as a varint-prefixed byte stream `[varint(len)][bytes]*` so the
/// downstream `StreamCursor::read_one_record_varint` works identically to
/// BSC-emitted header streams. Cost: ~2 small allocations per block.
pub(super) fn bounded_columnar_header_producer(
    archive_path: &std::path::Path,
    index: &StreamIndex,
    tx: std::sync::mpsc::SyncSender<Result<(u32, Vec<u8>), String>>,
    max_inflight: usize,
) -> Result<()> {
    run_parallel_producer(
        archive_path,
        index,
        tx,
        "columnar",
        max_inflight,
        |idx, record_count, payload| {
            // Payload layout: [chunk_reads: u32][columnar_blob]
            if payload.len() < 4 {
                return Err(format!(
                    "columnar block {idx}: payload too short ({} bytes) — missing chunk_reads",
                    payload.len()
                ));
            }
            let chunk_reads =
                u32::from_le_bytes([payload[0], payload[1], payload[2], payload[3]]) as usize;
            if chunk_reads as u32 != record_count {
                return Err(format!(
                    "columnar block {idx}: framing record_count={record_count} disagrees with \
                 payload chunk_reads={chunk_reads}"
                ));
            }
            let blob = &payload[4..];
            let varint_stream =
                match crate::compression::header_col::decompress_headers_columnar_to_varint_stream(
                    blob,
                    chunk_reads,
                ) {
                    Ok(v) => v,
                    Err(e) => {
                        return Err(format!("columnar block {idx} decompress: {e}"));
                    }
                };
            Ok((record_count, varint_stream))
        },
    )
}

/// Codec dispatch for the bounded-writer: which producer to spawn for a given
/// stream. `Columnar` is headers-only; the others carry raw byte payloads.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(super) enum StreamCodec {
    Bsc,
    Columnar, // headers only
    Fqz,      // qualities only (record-capped multi-block fqzcomp)
}

/// Channel capacity for bounded producer→writer hand-off. Two-block lookahead
/// keeps the writer fed without ballooning resident memory: peak per-stream
/// buffer is ~3 decompressed blocks (one being produced + 2 in the channel).
const BOUNDED_CHANNEL_CAP: usize = 2;

/// Spawn the right producer thread for the requested codec, attached to the
/// provided `scope`. The producer's `Result` is dropped on join — errors are
/// surfaced through the channel via the producer's existing `tx.send(Err(...))`
/// path. A future polish pass could keep the JoinHandle to also propagate
/// producer panics to the main thread.
pub(super) fn spawn_stream_producer<'scope>(
    scope: &'scope std::thread::Scope<'scope, '_>,
    archive_path: &std::path::Path,
    index: &StreamIndex,
    codec: StreamCodec,
    tx: std::sync::mpsc::SyncSender<Result<(u32, Vec<u8>), String>>,
    max_inflight: usize,
    max_block_len: usize,
) {
    let path = archive_path.to_path_buf();
    let index = index.clone();
    scope.spawn(move || {
        let res = match codec {
            StreamCodec::Bsc => {
                bounded_bsc_producer(&path, &index, tx.clone(), max_inflight, max_block_len)
            }
            StreamCodec::Columnar => {
                bounded_columnar_header_producer(&path, &index, tx.clone(), max_inflight)
            }
            StreamCodec::Fqz => bounded_fqz_producer(&path, &index, tx.clone(), max_inflight),
        };
        if let Err(e) = res {
            // Surface I/O / framing errors through the channel — the writer
            // sees them on the next `recv`. Ignore send failure (writer hung up).
            let _ = tx.send(Err(format!("{codec:?} producer error: {e}")));
        }
    });
}

/// Write `buf` to `output`, treating `BrokenPipe` (downstream consumer closed)
/// as a clean early-exit signal.
///
/// Returns:
/// - `Ok(false)` on a successful write.
/// - `Ok(true)`  when the downstream pipe was closed (`ErrorKind::BrokenPipe`);
///   the caller should stop writing and return `Ok(())`.
/// - `Err(_)`    for all other I/O errors (still propagated to the caller).
#[inline]
fn write_or_broken_pipe(output: &mut dyn std::io::Write, buf: &[u8]) -> Result<bool> {
    match output.write_all(buf) {
        Ok(()) => Ok(false),
        Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => Ok(true),
        Err(e) => Err(e.into()),
    }
}

/// Output write-buffer flush threshold (~2 MiB): amortize syscalls without
/// holding the whole FASTQ in memory. Shared by both bounded writers.
pub(super) const WRITE_BATCH: usize = 2 * 1024 * 1024;

/// Drain one sequence record from the sequence cursor, handling the
/// const-length vs varint framing and stripping the leading sequence-hint byte
/// when present. Single source of truth shared by both bounded writers (their
/// const/varint/hint branching must stay identical or sequences corrupt).
#[inline]
fn drain_one_sequence(
    s_cur: &mut StreamCursor,
    const_seq_len: usize,
    has_sequence_hints: bool,
) -> Result<Vec<u8>> {
    if const_seq_len > 0 {
        let rec_len = const_seq_len + if has_sequence_hints { 1 } else { 0 };
        let bytes = s_cur.read_one_record_const(rec_len)?;
        Ok(if has_sequence_hints {
            bytes[1..].to_vec()
        } else {
            bytes
        })
    } else {
        s_cur.read_one_record_seq_with_hint(has_sequence_hints)
    }
}

/// Append one formatted sub-chunk to `buf`, flushing to `output` once it crosses
/// `WRITE_BATCH`. Returns `Ok(true)` when the downstream pipe closed (caller
/// should stop). Shared inner write loop of both bounded writers.
#[inline]
pub(super) fn push_and_flush(
    output: &mut dyn std::io::Write,
    buf: &mut Vec<u8>,
    chunk_out: &[u8],
) -> Result<bool> {
    buf.extend_from_slice(chunk_out);
    if buf.len() >= WRITE_BATCH {
        if write_or_broken_pipe(output, buf)? {
            return Ok(true);
        }
        buf.clear();
    }
    Ok(false)
}

/// Final flush for a bounded writer: write any buffered bytes, then flush,
/// treating a downstream `BrokenPipe` (consumer closed mid-write/flush) as
/// clean success. Shared tail of both bounded writers.
pub(super) fn flush_remaining(output: &mut dyn std::io::Write, buf: &[u8]) -> Result<()> {
    if !buf.is_empty() && write_or_broken_pipe(output, buf)? {
        return Ok(());
    }
    match output.flush() {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => Ok(()),
        Err(e) => Err(e.into()),
    }
}

/// Quality input variant for `assemble_fastq_record`.
///
/// The bounded writers see qualities in two shapes depending on the codec:
/// * `Packed` — variable- or const-length records with bit-packed Phred bytes
///   that need `columnar::unpack_qualities_to_writer` to decode.
/// * `Decoded` — already-unpacked ASCII Phred+33 bytes (the fqzcomp path, where
///   the producer decoded each block to final per-record quality bytes).
/// * `None` — `--discard-quality` archives (no quality stream).
#[derive(Clone, Copy)]
pub(super) enum QualityInput<'a> {
    None,
    Packed {
        l_seq: usize,
        packed: &'a [u8],
        binning: QualityBinning,
    },
    Decoded(&'a [u8]),
}

/// Per-record batch buffer size used by the parallel-format bounded writers.
///
/// Trade-off: smaller batches reduce transient memory but amortize dispatch
/// overhead worse; larger batches do the opposite. 25 K records × ~350 B/record
/// ≈ 9 MiB transient input buffer + ~9 MiB output buffer per batch — cheap.
pub(super) const FORMAT_BATCH_SIZE: usize = 25_000;

/// Sub-chunk size for rayon par_chunks inside a batch.
///
/// Splits the batch so each rayon worker formats ~SUB_CHUNK records into one
/// `Vec<u8>` output. With 25 K records and 512 per worker we get ~50 chunks —
/// enough parallelism to saturate 16-72 cores without per-record allocator
/// pressure.
const FORMAT_SUB_CHUNK: usize = 512;

/// Exact decoded output byte length for one record, matching `assemble_fastq_record`
/// byte-for-byte. `header_len` INCLUDES the leading '@'/'>' (FastqRecord.id stores
/// it — io::fastq tests confirm). `qual_len` is 0 for FASTA and for no-quality FASTQ
/// (which still emits the empty `\n+\n\n`). RC preserves length. Single source of
/// truth for the per-chunk size table; cross-checked by the test
/// `fastq_output_len_matches_formatter` below.
pub(crate) fn fastq_output_len(is_fasta: bool, header_len: usize, seq_len: usize, qual_len: usize) -> u64 {
    if is_fasta { (header_len + seq_len + 2) as u64 } else { (header_len + seq_len + qual_len + 5) as u64 }
}

/// Assemble one FASTQ record into `out`, applying RC and quality unpacking
/// as needed. Pure function — no I/O, no shared state.
///
/// Output bytes: `<header>\n<sequence>\n+\n<quality>\n` (or
/// `<header>\n<sequence>\n` when `quality == QualityInput::None`).
///
/// Invariants matching the legacy per-record loop:
/// * Header bytes already start with `@` (Revision C); we never prepend it.
/// * `rc_flag != 0` ⇒ emit `reverse_complement_canonical(sequence)` instead of
///   `sequence` (the lossless involution the encoder canonicalized with).
///   Quality bytes are NOT reversed: the encoder stored qualities aligned to
///   the canonical sequence (or to the original when rc_canon is off), so the
///   raw order matches whichever direction the writer emits.
/// * `QualityInput::Packed` uses `columnar::unpack_qualities_to_writer` —
///   identical to the serial path.
#[inline]
fn assemble_fastq_record(
    header: &[u8],
    sequence: &[u8],
    rc_flag: u8,
    quality: QualityInput<'_>,
    is_fasta: bool,
    out: &mut Vec<u8>,
) -> Result<()> {
    out.extend_from_slice(header);
    out.push(b'\n');

    if rc_flag != 0 {
        let rc_seq = dna_utils::reverse_complement_canonical(sequence);
        out.extend_from_slice(&rc_seq);
    } else {
        out.extend_from_slice(sequence);
    }

    if is_fasta {
        // FASTA: header ('>' prefix) + sequence only — no '+'/quality lines.
        out.push(b'\n');
        return Ok(());
    }

    out.extend_from_slice(b"\n+\n");
    match quality {
        QualityInput::None => {
            // Discarded-quality archives still emit the 4-line FASTQ shape
            // with an empty quality line: `\n+\n\n`. Matches the serial
            // legacy writer.
        }
        QualityInput::Packed {
            l_seq,
            packed,
            binning,
        } => {
            columnar::unpack_qualities_to_writer(packed, l_seq, binning, out)
                .map_err(|e| anyhow::anyhow!("quality unpack: {e}"))?;
        }
        QualityInput::Decoded(qual) => {
            out.extend_from_slice(qual);
        }
    }
    out.push(b'\n');
    Ok(())
}

/// Parallel FASTQ assembly over a batch of pre-pulled record inputs.
///
/// Splits the batch into `FORMAT_SUB_CHUNK`-sized slices, formats each slice
/// into its own `Vec<u8>` via rayon, and returns the per-slice outputs in
/// original order. The caller writes them sequentially to preserve FASTQ
/// record order.
///
/// Why slices, not individual records: per-record `Vec` allocation × 25 K
/// records would add millions of allocator calls per batch. Sub-chunks
/// amortize allocation while still keeping `num_threads × est_chunk_bytes`
/// well under 1 GiB of transient memory.
pub(super) fn format_batch_parallel(
    headers: &[Vec<u8>],
    sequences: &[Vec<u8>],
    rc_flags: &[u8], // empty slice = no RC
    quality_inputs: &[QualityInput<'_>],
    is_fasta: bool,
    est_per_record: usize,
) -> Result<Vec<Vec<u8>>> {
    use rayon::prelude::*;

    let n = headers.len();
    debug_assert_eq!(sequences.len(), n);
    debug_assert_eq!(quality_inputs.len(), n);
    debug_assert!(rc_flags.is_empty() || rc_flags.len() == n);

    let has_rc = !rc_flags.is_empty();

    (0..n)
        .into_par_iter()
        .step_by(FORMAT_SUB_CHUNK)
        .map(|start| -> Result<Vec<u8>> {
            let end = (start + FORMAT_SUB_CHUNK).min(n);
            let mut out = Vec::with_capacity((end - start) * est_per_record + 64);
            for i in start..end {
                let rc = if has_rc { rc_flags[i] } else { 0 };
                assemble_fastq_record(
                    &headers[i],
                    &sequences[i],
                    rc,
                    quality_inputs[i],
                    is_fasta,
                    &mut out,
                )?;
            }
            Ok(out)
        })
        .collect()
}

/// Bounded streaming writer: emit records in parallel-formatted batches from
/// per-stream cursors fed by codec-specific producers. Used when
/// `quality_compressor` is BSC, Fqzcomp, or no-quality.
///
/// Memory: `(BOUNDED_CHANNEL_CAP + MAX_PARALLEL_DECOMPRESS) × max_block_size`
/// per active stream (one in the cursor + 2 in the mpsc channel + up to
/// `MAX_PARALLEL_DECOMPRESS` in flight on rayon) — bounded by the producer
/// pipeline, see `run_parallel_producer`. Plus one batch of pulled record
/// bytes (~`FORMAT_BATCH_SIZE × avg_record_size`, ≈ 9-15 MiB at 150 bp /
/// 25 K records) and the per-batch parallel output (~same magnitude).
///
/// FASTQ emit shape:  `<header>\n<sequence>\n+\n<quality>\n`
/// Header bytes already include the leading `@` per Revision C of the spec.
pub(super) fn bounded_write_records_independent(
    output: &mut dyn std::io::Write,
    num_reads: u64,
    archive_path: &std::path::Path,
    headers_idx: &StreamIndex,
    sequences_idx: &StreamIndex,
    qualities_idx: Option<&StreamIndex>,
    rc_idx: Option<&StreamIndex>,
    header_codec: StreamCodec,
    sequence_codec: StreamCodec,
    quality_codec: Option<StreamCodec>,
    const_seq_len: usize,
    const_qual_len: usize,
    has_sequence_hints: bool,
    bits_per_qual: usize,
    quality_binning: QualityBinning,
    is_fasta: bool,
    max_inflight: usize,
    max_block_len: usize,
) -> Result<()> {
    use std::sync::mpsc::sync_channel;

    if num_reads == 0 {
        return output.flush().map_err(Into::into);
    }

    // Allocate the per-stream channels OUTSIDE the scope so the receivers can
    // move into the cursors (which live in the parent thread) while the senders
    // move into producer threads spawned inside `scope`.
    let (h_tx, h_rx) = sync_channel::<Result<(u32, Vec<u8>), String>>(BOUNDED_CHANNEL_CAP);
    let (s_tx, s_rx) = sync_channel::<Result<(u32, Vec<u8>), String>>(BOUNDED_CHANNEL_CAP);
    // Fqzcomp quality decodes to final raw-ASCII bytes (emitted via
    // `QualityInput::Decoded`); BSC quality stays bit-packed (`Packed`).
    let is_fqz_quality = matches!(quality_codec, Some(StreamCodec::Fqz));
    let q_pair = if qualities_idx.is_some() && quality_codec.is_some() {
        Some(sync_channel::<Result<(u32, Vec<u8>), String>>(
            BOUNDED_CHANNEL_CAP,
        ))
    } else {
        None
    };
    let rc_pair = if rc_idx.is_some() {
        Some(sync_channel::<Result<(u32, Vec<u8>), String>>(
            BOUNDED_CHANNEL_CAP,
        ))
    } else {
        None
    };

    std::thread::scope(|scope| -> Result<()> {
        // Spawn one producer per active stream. Headers + sequences are always
        // active; qualities and rc-flags are conditional.
        spawn_stream_producer(
            scope,
            archive_path,
            headers_idx,
            header_codec,
            h_tx,
            max_inflight,
            max_block_len,
        );
        spawn_stream_producer(
            scope,
            archive_path,
            sequences_idx,
            sequence_codec,
            s_tx,
            max_inflight,
            max_block_len,
        );
        if let (Some(idx), Some(codec), Some((tx, _))) =
            (qualities_idx, quality_codec, q_pair.as_ref())
        {
            spawn_stream_producer(
                scope,
                archive_path,
                idx,
                codec,
                tx.clone(),
                max_inflight,
                max_block_len,
            );
        }
        if let (Some(idx), Some((tx, _))) = (rc_idx, rc_pair.as_ref()) {
            // RC flags are always v4-framed BSC blocks (n_masks/rc_flags share
            // the same codec: BSC raw bytes).
            spawn_stream_producer(
                scope,
                archive_path,
                idx,
                StreamCodec::Bsc,
                tx.clone(),
                max_inflight,
                max_block_len,
            );
        }

        let mut h_cur = StreamCursor::new(h_rx, "headers");
        let mut s_cur = StreamCursor::new(s_rx, "sequences");
        let mut q_cur = q_pair.map(|(_, rx)| StreamCursor::new(rx, "qualities"));
        let mut rc_cur = rc_pair.map(|(_, rx)| StreamCursor::new(rx, "rc_flags"));

        let est_seq = if const_seq_len > 0 {
            const_seq_len
        } else {
            150
        };
        let est_per_record = 80 + est_seq + 4 + est_seq + 1;
        let mut buf: Vec<u8> = Vec::with_capacity(WRITE_BATCH + est_per_record + 256);

        // Process records in `FORMAT_BATCH_SIZE` chunks: drain inputs serially
        // from the cursors (cursors are stateful and not `Sync`), then format
        // the batch into FASTQ bytes in parallel via `format_batch_parallel`.
        //
        // Memory note: peak transient memory = one batch of pulled record bytes
        // (~`FORMAT_BATCH_SIZE × avg_record_size`, ~9 MiB at 150 bp) plus the
        // parallel output (~same magnitude).
        let mut emitted: u64 = 0;
        let mut headers: Vec<Vec<u8>> = Vec::with_capacity(FORMAT_BATCH_SIZE);
        let mut sequences: Vec<Vec<u8>> = Vec::with_capacity(FORMAT_BATCH_SIZE);
        let mut rc_flags: Vec<u8> = if rc_cur.is_some() {
            Vec::with_capacity(FORMAT_BATCH_SIZE)
        } else {
            Vec::new()
        };
        // Quality is stored as owned packed bytes per record (or empty when
        // no quality stream); `QualityInput` borrows into these on assembly.
        let mut packed_quals: Vec<Vec<u8>> = if q_cur.is_some() {
            Vec::with_capacity(FORMAT_BATCH_SIZE)
        } else {
            Vec::new()
        };
        let mut qual_lengths: Vec<usize> = if q_cur.is_some() {
            Vec::with_capacity(FORMAT_BATCH_SIZE)
        } else {
            Vec::new()
        };

        while emitted < num_reads {
            headers.clear();
            sequences.clear();
            rc_flags.clear();
            packed_quals.clear();
            qual_lengths.clear();

            let batch_n = (num_reads - emitted).min(FORMAT_BATCH_SIZE as u64) as usize;

            // Drain the batch from the cursors (serial — cursors hold &mut state).
            for _ in 0..batch_n {
                headers.push(h_cur.read_one_record_varint()?);

                sequences.push(drain_one_sequence(
                    &mut s_cur,
                    const_seq_len,
                    has_sequence_hints,
                )?);

                if let Some(rc) = rc_cur.as_mut() {
                    rc_flags.push(rc.read_one_byte()?);
                }

                if let Some(q) = q_cur.as_mut() {
                    let l_seq = if const_qual_len > 0 {
                        let packed_len = (const_qual_len * bits_per_qual).div_ceil(8);
                        packed_quals.push(q.read_one_record_const(packed_len)?);
                        const_qual_len
                    } else {
                        let (l_seq, packed) = q.read_one_record_qual_varint(bits_per_qual)?;
                        packed_quals.push(packed);
                        l_seq
                    };
                    // FASTQ requires the quality line to match the sequence length.
                    // Sequence and quality lengths come from independent streams, so
                    // a corrupt/crafted archive could disagree — reject rather than
                    // emit a structurally invalid record.
                    let seq_len = sequences.last().map_or(0, |s| s.len());
                    if l_seq != seq_len {
                        anyhow::bail!(
                            "record {}: quality length {} != sequence length {} (corrupt archive)",
                            emitted + qual_lengths.len() as u64,
                            l_seq,
                            seq_len
                        );
                    }
                    qual_lengths.push(l_seq);
                }
            }

            // Build per-record QualityInput borrows for this batch. Fqzcomp
            // blocks decode to FINAL raw-ASCII quality bytes (the producer ran
            // `decompress_qualities_fqzcomp`, and `bits_per_qual == 8` made the
            // varint read pull all `l_seq` bytes), so emit them directly via
            // `Decoded`. The `Packed` path would re-interpret them as bit-packed
            // Phred and corrupt the output. BSC quality stays `Packed`.
            let quality_inputs: Vec<QualityInput<'_>> = if q_cur.is_some() {
                if is_fqz_quality {
                    (0..batch_n)
                        .map(|i| QualityInput::Decoded(&packed_quals[i]))
                        .collect()
                } else {
                    (0..batch_n)
                        .map(|i| QualityInput::Packed {
                            l_seq: qual_lengths[i],
                            packed: &packed_quals[i],
                            binning: quality_binning,
                        })
                        .collect()
                }
            } else {
                vec![QualityInput::None; batch_n]
            };

            // Format the batch in parallel.
            let assembled = format_batch_parallel(
                &headers,
                &sequences,
                &rc_flags,
                &quality_inputs,
                is_fasta,
                est_per_record,
            )?;

            // Write each sub-chunk's bytes sequentially to preserve order.
            for chunk_out in &assembled {
                if push_and_flush(output, &mut buf, chunk_out)? {
                    return Ok(());
                }
            }

            emitted += batch_n as u64;
        }

        flush_remaining(output, &buf)
    })
}

/// Cursor over the R2 `HeaderDelta` op stream for the unified paired writer.
///
/// The delta op stream is tiny (~1-3 B/record — typically one substituted mate
/// field byte), so on first use we drain the whole BSC-decoded stream into one
/// `ops` buffer, then apply ops **record-by-record** against the lockstep R1 id via
/// `header_delta::apply_one_op`. Peak extra RAM = the whole ops stream, a small
/// fraction of the FASTQ (and far less than materialising a chunk).
struct R2DeltaCursor {
    rx: std::sync::mpsc::Receiver<Result<(u32, Vec<u8>), String>>,
    ops: Vec<u8>,
    o: usize,
    loaded: bool,
    /// Cumulative decompressed-size cap for the concatenated op stream (a
    /// content-derived bound). A forged HeaderDelta stream could otherwise chain
    /// many valid blocks into an unbounded allocation, so refuse past the cap —
    /// mirrors the `max_total` guard `decode_bsc_stream` applies on the legacy path.
    cap: usize,
}
impl R2DeltaCursor {
    fn new(rx: std::sync::mpsc::Receiver<Result<(u32, Vec<u8>), String>>, cap: usize) -> Self {
        Self {
            rx,
            ops: Vec::new(),
            o: 0,
            loaded: false,
            cap,
        }
    }
    fn ensure_loaded(&mut self) -> Result<()> {
        if self.loaded {
            return Ok(());
        }
        // Drain every decoded block (in order) until the producer drops its sender.
        while let Ok(msg) = self.rx.recv() {
            let (_rc, block) =
                msg.map_err(|e| anyhow::anyhow!("R2 header-delta producer: {e}"))?;
            if self.ops.len().saturating_add(block.len()) > self.cap {
                anyhow::bail!(
                    "R2 header-delta op stream exceeds cap of {} bytes (decompression bomb?)",
                    self.cap
                );
            }
            self.ops.extend_from_slice(&block);
        }
        self.loaded = true;
        Ok(())
    }
    fn next_id(&mut self, r1: &[u8]) -> Result<Vec<u8>> {
        self.ensure_loaded()?;
        crate::compression::paired::header_delta::apply_one_op(&self.ops, &mut self.o, r1)
    }
}

/// Drain one quality record from a mate's quality cursor into `packed_quals` /
/// `qual_lengths`, mirroring the single-end const/varint branching, and reject a
/// quality length that disagrees with the paired sequence length. Shared by both
/// mates of `bounded_write_records_paired`.
#[inline]
pub(super) fn drain_one_qual(
    q: &mut StreamCursor,
    const_qual_len: usize,
    bits_per_qual: usize,
    seq_len: usize,
    record_global_idx: u64,
    packed_quals: &mut Vec<Vec<u8>>,
    qual_lengths: &mut Vec<usize>,
) -> Result<()> {
    let l_seq = if const_qual_len > 0 {
        let packed_len = (const_qual_len * bits_per_qual).div_ceil(8);
        packed_quals.push(q.read_one_record_const(packed_len)?);
        const_qual_len
    } else {
        let (l_seq, packed) = q.read_one_record_qual_varint(bits_per_qual)?;
        packed_quals.push(packed);
        l_seq
    };
    if l_seq != seq_len {
        anyhow::bail!(
            "record {record_global_idx}: quality length {l_seq} != sequence length {seq_len} (corrupt archive)"
        );
    }
    qual_lengths.push(l_seq);
    Ok(())
}

/// Output target for the paired bounded streaming writer. `Two` writes R1→`out1`
/// and R2→`out2` (the default two-file / two-region paired decode). `Interleaved`
/// merges both mates into ONE sink, emitting record-aligned pairs alternately
/// (`r0/1, r0/2, r1/1, …`) — the form `qz decompress --interleaved` produces for
/// aligners like `bwa mem -p` / `strobealign --interleaved`.
pub(super) enum PairedSink<'a> {
    Two(&'a mut dyn std::io::Write, &'a mut dyn std::io::Write),
    Interleaved(&'a mut dyn std::io::Write),
}

impl PairedSink<'_> {
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            PairedSink::Two(a, b) => {
                a.flush()?;
                b.flush()
            }
            PairedSink::Interleaved(w) => w.flush(),
        }
    }
}

/// Paired bounded streaming writer — the paired adapter over the same
/// per-stream block-parallel producer + cursor + batch-assembly primitives the
/// single-end `bounded_write_records_independent` uses. Runs R1's and R2's stream
/// sets (headers/sequence/quality) concurrently in one `std::thread::scope` and
/// reconstructs them record-aligned (record i of R1 pairs with record i of R2);
/// the `PairedSink` then either writes the two mates to separate sinks (`Two`) or
/// interleaves them into one (`Interleaved`). R2 headers are either an independent
/// columnar stream or a `HeaderDelta` op stream applied per-record against the
/// lockstep R1 id (`r2.header_is_delta`).
///
/// Replaces the old serial chunk-materialising `decompress_paired_v5` core. Memory
/// is block-bounded per stream (no chunk materialisation), like single-end.
#[allow(clippy::too_many_arguments)]
pub(super) fn bounded_write_records_paired(
    sink: &mut PairedSink<'_>,
    num_pairs: u64,
    archive_path: &std::path::Path,
    r1: &PairedMateIndices,
    r2: &PairedMateIndices,
    const_seq_len: usize,
    const_qual_len: usize,
    has_sequence_hints: bool,
    bits_per_qual: usize,
    quality_binning: QualityBinning,
    is_fasta: bool,
    max_inflight: usize,
    max_block_len: usize,
) -> Result<()> {
    use std::sync::mpsc::sync_channel;

    if num_pairs == 0 {
        return sink.flush().map_err(Into::into);
    }

    let chan = || sync_channel::<Result<(u32, Vec<u8>), String>>(BOUNDED_CHANNEL_CAP);
    // R1 streams: headers (columnar), sequence, optional quality.
    let (h1_tx, h1_rx) = chan();
    let (s1_tx, s1_rx) = chan();
    let q1_pair = r1.qualities.as_ref().map(|_| chan());
    // R2 streams: headers (columnar) OR header-delta ops, sequence, optional quality.
    let (h2_tx, h2_rx) = chan();
    let (s2_tx, s2_rx) = chan();
    let q2_pair = r2.qualities.as_ref().map(|_| chan());

    let is_fqz_q1 = matches!(r1.quality_codec, Some(StreamCodec::Fqz));
    let is_fqz_q2 = matches!(r2.quality_codec, Some(StreamCodec::Fqz));

    std::thread::scope(|scope| -> Result<()> {
        // R1 producers.
        spawn_stream_producer(
            scope,
            archive_path,
            &r1.headers,
            StreamCodec::Columnar,
            h1_tx,
            max_inflight,
            max_block_len,
        );
        spawn_stream_producer(
            scope,
            archive_path,
            &r1.sequences,
            StreamCodec::Bsc,
            s1_tx,
            max_inflight,
            max_block_len,
        );
        if let (Some(idx), Some((tx, _))) = (r1.qualities.as_ref(), q1_pair.as_ref()) {
            spawn_stream_producer(
                scope,
                archive_path,
                idx,
                if is_fqz_q1 {
                    StreamCodec::Fqz
                } else {
                    StreamCodec::Bsc
                },
                tx.clone(),
                max_inflight,
                max_block_len,
            );
        }
        // R2 producers. Header stream is columnar OR a BSC op stream (delta).
        spawn_stream_producer(
            scope,
            archive_path,
            &r2.headers,
            if r2.header_is_delta {
                StreamCodec::Bsc
            } else {
                StreamCodec::Columnar
            },
            h2_tx,
            max_inflight,
            max_block_len,
        );
        spawn_stream_producer(
            scope,
            archive_path,
            &r2.sequences,
            StreamCodec::Bsc,
            s2_tx,
            max_inflight,
            max_block_len,
        );
        if let (Some(idx), Some((tx, _))) = (r2.qualities.as_ref(), q2_pair.as_ref()) {
            spawn_stream_producer(
                scope,
                archive_path,
                idx,
                if is_fqz_q2 {
                    StreamCodec::Fqz
                } else {
                    StreamCodec::Bsc
                },
                tx.clone(),
                max_inflight,
                max_block_len,
            );
        }

        // Cursors. R2 headers: columnar cursor OR delta op-cursor.
        let mut h1_cur = StreamCursor::new(h1_rx, "r1_headers");
        let mut s1_cur = StreamCursor::new(s1_rx, "r1_sequences");
        let mut q1_cur = q1_pair.map(|(_, rx)| StreamCursor::new(rx, "r1_qualities"));
        let mut s2_cur = StreamCursor::new(s2_rx, "r2_sequences");
        let mut q2_cur = q2_pair.map(|(_, rx)| StreamCursor::new(rx, "r2_qualities"));
        // Content-derived cap for the concatenated R2 delta op stream (per-record
        // bound from the pair count), same family of bound the legacy BSC stream
        // decode uses.
        let delta_cap = crate::compression::codecs::stream_decode_cap(num_pairs as usize);
        let (mut h2_cur, mut r2_delta) = if r2.header_is_delta {
            (None, Some(R2DeltaCursor::new(h2_rx, delta_cap)))
        } else {
            (Some(StreamCursor::new(h2_rx, "r2_headers")), None)
        };

        let est_seq = if const_seq_len > 0 { const_seq_len } else { 150 };
        let est_per_record = 80 + est_seq + 4 + est_seq + 1;
        let mut buf1: Vec<u8> = Vec::with_capacity(WRITE_BATCH + est_per_record + 256);
        let mut buf2: Vec<u8> = Vec::with_capacity(WRITE_BATCH + est_per_record + 256);

        // Per-mate reusable batch buffers.
        let mut h1: Vec<Vec<u8>> = Vec::with_capacity(FORMAT_BATCH_SIZE);
        let mut s1: Vec<Vec<u8>> = Vec::with_capacity(FORMAT_BATCH_SIZE);
        let mut pq1: Vec<Vec<u8>> = Vec::new();
        let mut ql1: Vec<usize> = Vec::new();
        let mut h2: Vec<Vec<u8>> = Vec::with_capacity(FORMAT_BATCH_SIZE);
        let mut s2: Vec<Vec<u8>> = Vec::with_capacity(FORMAT_BATCH_SIZE);
        let mut pq2: Vec<Vec<u8>> = Vec::new();
        let mut ql2: Vec<usize> = Vec::new();
        let no_rc: Vec<u8> = Vec::new(); // paired never emits RcFlags

        let mut emitted: u64 = 0;
        while emitted < num_pairs {
            h1.clear();
            s1.clear();
            pq1.clear();
            ql1.clear();
            h2.clear();
            s2.clear();
            pq2.clear();
            ql2.clear();

            let batch_n = (num_pairs - emitted).min(FORMAT_BATCH_SIZE as u64) as usize;

            for k in 0..batch_n {
                // R1 record.
                h1.push(h1_cur.read_one_record_varint()?);
                s1.push(drain_one_sequence(&mut s1_cur, const_seq_len, has_sequence_hints)?);
                if let Some(q) = q1_cur.as_mut() {
                    let seq_len = s1.last().map_or(0, |s| s.len());
                    drain_one_qual(
                        q,
                        const_qual_len,
                        bits_per_qual,
                        seq_len,
                        emitted + k as u64,
                        &mut pq1,
                        &mut ql1,
                    )?;
                }

                // R2 record. Header from columnar cursor or delta vs the R1 id we
                // just drained (record-aligned).
                let r2_id = if let Some(dc) = r2_delta.as_mut() {
                    dc.next_id(h1.last().unwrap())?
                } else {
                    h2_cur.as_mut().unwrap().read_one_record_varint()?
                };
                h2.push(r2_id);
                s2.push(drain_one_sequence(&mut s2_cur, const_seq_len, has_sequence_hints)?);
                if let Some(q) = q2_cur.as_mut() {
                    let seq_len = s2.last().map_or(0, |s| s.len());
                    drain_one_qual(
                        q,
                        const_qual_len,
                        bits_per_qual,
                        seq_len,
                        emitted + k as u64,
                        &mut pq2,
                        &mut ql2,
                    )?;
                }
            }

            let qi1 = build_quality_inputs(q1_cur.is_some(), is_fqz_q1, batch_n, &pq1, &ql1, quality_binning);
            let qi2 = build_quality_inputs(q2_cur.is_some(), is_fqz_q2, batch_n, &pq2, &ql2, quality_binning);

            match sink {
                PairedSink::Two(out1, out2) => {
                    let a1 = format_batch_parallel(&h1, &s1, &no_rc, &qi1, is_fasta, est_per_record)?;
                    let a2 = format_batch_parallel(&h2, &s2, &no_rc, &qi2, is_fasta, est_per_record)?;
                    for chunk_out in &a1 {
                        if push_and_flush(&mut **out1, &mut buf1, chunk_out)? {
                            return Ok(());
                        }
                    }
                    for chunk_out in &a2 {
                        if push_and_flush(&mut **out2, &mut buf2, chunk_out)? {
                            return Ok(());
                        }
                    }
                }
                PairedSink::Interleaved(out) => {
                    // Merge R1[k], R2[k] into one batch of 2·batch_n records, then format
                    // once → a single sink. Per-record `QualityInput`s borrow pq1/pq2,
                    // which outlive the batch, so the moved borrows stay valid.
                    let mut hi: Vec<Vec<u8>> = Vec::with_capacity(2 * batch_n);
                    let mut si: Vec<Vec<u8>> = Vec::with_capacity(2 * batch_n);
                    let mut qii: Vec<_> = Vec::with_capacity(2 * batch_n);
                    let mut h1i = h1.drain(..);
                    let mut s1i = s1.drain(..);
                    let mut h2i = h2.drain(..);
                    let mut s2i = s2.drain(..);
                    let mut q1i = qi1.into_iter();
                    let mut q2i = qi2.into_iter();
                    for _ in 0..batch_n {
                        hi.push(h1i.next().unwrap());
                        si.push(s1i.next().unwrap());
                        qii.push(q1i.next().unwrap());
                        hi.push(h2i.next().unwrap());
                        si.push(s2i.next().unwrap());
                        qii.push(q2i.next().unwrap());
                    }
                    drop((h1i, s1i, h2i, s2i));
                    let ai = format_batch_parallel(&hi, &si, &no_rc, &qii, is_fasta, est_per_record)?;
                    for chunk_out in &ai {
                        if push_and_flush(&mut **out, &mut buf1, chunk_out)? {
                            return Ok(());
                        }
                    }
                }
            }

            emitted += batch_n as u64;
        }

        // Full-consumption invariant (the legacy path's per-chunk `o != ops.len()`
        // check): a malformed/forged HeaderDelta stream with trailing bytes after
        // the last record must error, not be silently ignored.
        if let Some(dc) = r2_delta.as_ref()
            && dc.o != dc.ops.len()
        {
            anyhow::bail!(
                "R2 header-delta op stream has {} trailing bytes (corrupt archive)",
                dc.ops.len() - dc.o
            );
        }

        match sink {
            PairedSink::Two(out1, out2) => {
                flush_remaining(&mut **out1, &buf1)?;
                flush_remaining(&mut **out2, &buf2)
            }
            PairedSink::Interleaved(out) => flush_remaining(&mut **out, &buf1),
        }
    })
}

/// Build the per-record `QualityInput` borrows for one mate's batch (mirrors the
/// single-end logic): fqzcomp → already-decoded ASCII (`Decoded`); BSC → bit-packed
/// (`Packed`); no quality stream → `None`.
pub(super) fn build_quality_inputs<'a>(
    has_quality: bool,
    is_fqz: bool,
    batch_n: usize,
    packed_quals: &'a [Vec<u8>],
    qual_lengths: &'a [usize],
    quality_binning: QualityBinning,
) -> Vec<QualityInput<'a>> {
    if !has_quality {
        return vec![QualityInput::None; batch_n];
    }
    if is_fqz {
        (0..batch_n)
            .map(|i| QualityInput::Decoded(&packed_quals[i]))
            .collect()
    } else {
        (0..batch_n)
            .map(|i| QualityInput::Packed {
                l_seq: qual_lengths[i],
                packed: &packed_quals[i],
                binning: quality_binning,
            })
            .collect()
    }
}

/// Streaming paired decode driver: derive writer params from the v5 front header,
/// build per-mate stream indices, and run the two-sink block-streaming writer.
///
/// The caller (`paired::decompress_paired_v5`) owns the front-header parse, footer
/// read + `validate_paired_directory`, and the output tmp-file/publish guard.
/// `w1`/`w2` receive R1/R2 FASTQ, record-aligned.
pub(super) fn decompress_paired_streaming_v5(
    archive_path: &std::path::Path,
    hdr: &super::FixedHeader,
    dir: &crate::compression::chunk_directory::ChunkDirectory,
    sink: &mut PairedSink<'_>,
) -> Result<()> {
    let r1 = build_paired_mate_indices(archive_path, dir, 1)?;
    let r2 = build_paired_mate_indices(archive_path, dir, 2)?;
    bounded_write_records_paired_from_indices(archive_path, hdr, &r1, &r2, sink)
}

/// Tail of [`decompress_paired_streaming_v5`], factored so the full paired path
/// and the NUMA chunk-range path share it byte-for-byte: given prebuilt per-mate
/// indices, derive the block bounds / quality binning / per-archive `bits_per_qual`,
/// cross-check the bulk streams agree on the pair count, then drive
/// [`bounded_write_records_paired`]. Pure factoring — the full path's behavior is
/// unchanged (it now passes whole-archive `r1`/`r2`; the NUMA path passes
/// chunk-range-sliced indices). The pair count is taken from the supplied indices
/// (so a slice writes only that range's pairs).
pub(super) fn bounded_write_records_paired_from_indices(
    archive_path: &std::path::Path,
    hdr: &super::FixedHeader,
    r1: &PairedMateIndices,
    r2: &PairedMateIndices,
    sink: &mut PairedSink<'_>,
) -> Result<()> {
    let (max_block_len, max_inflight) = decode_block_bounds(hdr.encoding_type);
    let quality_binning = code_to_binning(hdr.quality_binning_code)?;

    // Fqzcomp quality decodes to FINAL raw 8-bit ASCII (varint-framed: [len][len
    // bytes]); the per-record varint reader must use bits_per_qual = 8 so packed_len
    // == len. Only BSC-packed quality uses the binning's bits-per-quality. (The
    // legacy split_varint_records path ignored bits, masking this.) Quality codec is
    // uniform per archive, so one value covers both mates — enforce that here. A crafted
    // archive mixing fqz on one mate and BSC on the other would mis-decode under the
    // single `bits_per_qual` below (not encoder-reachable, but reject it explicitly).
    if r1.quality_codec != r2.quality_codec {
        anyhow::bail!(
            "paired decode: R1 and R2 quality codecs differ ({:?} vs {:?}) — corrupt archive",
            r1.quality_codec,
            r2.quality_codec
        );
    }
    let any_fqz = matches!(r1.quality_codec, Some(StreamCodec::Fqz))
        || matches!(r2.quality_codec, Some(StreamCodec::Fqz));
    let bits_per_qual = if any_fqz {
        8
    } else {
        quality_binning.bits_per_quality()
    };

    // num pairs = mate-1 header record total; cross-check the other bulk streams
    // agree (mispaired mates would otherwise corrupt silently).
    let num_pairs = r1.headers.total_records;
    for (name, idx) in [
        ("r1.sequences", &r1.sequences),
        ("r2.headers", &r2.headers),
        ("r2.sequences", &r2.sequences),
    ] {
        if idx.total_records != num_pairs {
            anyhow::bail!(
                "paired decode: {name} has {} records != {num_pairs} pairs",
                idx.total_records
            );
        }
    }

    bounded_write_records_paired(
        sink,
        num_pairs,
        archive_path,
        r1,
        r2,
        hdr.const_seq_len as usize,
        hdr.const_qual_len as usize,
        hdr.has_sequence_hints(),
        bits_per_qual,
        quality_binning,
        hdr.is_fasta(),
        max_inflight,
        max_block_len,
    )
}

/// Per-stream indices for a v5 (chunk-major) archive, plus the metadata the
/// streaming core needs. Built from the footer directory via `build_role_index`
/// and fed to the producer/cursor core.
pub(super) struct DecodePlan {
    pub headers: StreamIndex,
    pub sequences: StreamIndex,
    pub qualities: Option<StreamIndex>,
    pub rc_flags: Option<StreamIndex>,
    pub num_reads: u64,
}

/// Build one role's `StreamIndex`, optionally restricted to chunks `[start, end)`.
/// `None` = the whole archive (preserves the original `build_role_index` behavior).
/// Frame `offset`/`length` cover the 12-byte v4 header, so the payload starts
/// at `offset + 12` and spans `length - 12`.
pub(super) fn build_role_index_range(
    dir: &crate::compression::chunk_directory::ChunkDirectory,
    role: crate::compression::chunk_directory::StreamRole,
    chunk_range: Option<(u32, u32)>,
) -> Result<StreamIndex> {
    let mut blocks = Vec::new();
    let mut total_records = 0u64;
    for e in dir.entries.iter().filter(|e| {
        e.role == role
            && match chunk_range {
                Some((a, b)) => e.chunk_index >= a && e.chunk_index < b,
                None => true,
            }
    }) {
        let record_count = u32::try_from(e.record_count).map_err(|_| {
            anyhow::anyhow!(
                "v5 single-end entry record_count {} exceeds u32",
                e.record_count
            )
        })?;
        blocks.push(BlockEntry {
            record_count,
            byte_offset: e.offset + 12,
            byte_len: (e.length - 12) as u32,
            expected_crc: 0, // adopted from the inline header during validation
        });
        total_records += e.record_count;
    }
    Ok(StreamIndex {
        blocks,
        total_records,
    })
}

/// Build one role's `StreamIndex` from the directory, preserving entry order.
pub(super) fn build_role_index(
    dir: &crate::compression::chunk_directory::ChunkDirectory,
    role: crate::compression::chunk_directory::StreamRole,
) -> Result<StreamIndex> {
    build_role_index_range(dir, role, None)
}

/// Mate-filtered variant of [`build_role_index`] for the paired/reference decoders.
///
/// Two differences from the single-end builder:
/// 1. Filters frames by `(role, mate)` instead of role alone.
/// 2. Populates each block's `expected_crc` directly from its 12-byte inline frame
///    header. The bounded producers (`run_parallel_producer`) read only the payload
///    at `byte_offset` and verify against `BlockEntry.expected_crc`, so the CRC must
///    be adopted here (single-end does it in a separate pass over its single mate;
///    the per-mate paired builder folds it in). The inline `block_len`/`record_count`
///    are cross-checked against the footer (cheap corruption guard) while we're here.
///
/// Frames are taken in footer order, which the paired/reference footer validators
/// (`validate_paired_directory` / `validate_reference_directory`) constrain to chunk
/// order — so cursors pair record N across roles by position.
pub(super) fn build_role_index_mate_range(
    archive_path: &std::path::Path,
    dir: &crate::compression::chunk_directory::ChunkDirectory,
    role: crate::compression::chunk_directory::StreamRole,
    mate: u8,
    chunk_range: Option<(u32, u32)>,
) -> Result<StreamIndex> {
    use std::io::{Read, Seek, SeekFrom};
    let mut file = std::fs::File::open(archive_path)?;
    let mut blocks = Vec::new();
    let mut total_records = 0u64;
    // Consume this (role, mate)'s entries in CHUNK order. The footer is validated
    // for chunk contiguity but NOT for cross-entry ordering, so a forged/reordered
    // footer (CRC recomputed) could otherwise make the per-stream cursors yield
    // blocks out of chunk order → silent record MISPAIRING (the single-end engine
    // explicitly guards this in build_indices_from_footer). Sorting by chunk_index
    // gives the paired streaming path the same chunk-order guarantee and makes it
    // byte-identical to the legacy decode_chunk_v5, which processes chunks 0..n.
    let mut role_entries: Vec<_> = dir
        .entries
        .iter()
        .filter(|e| {
            e.role == role
                && e.mate == mate
                && chunk_range
                    .map(|(a, b)| e.chunk_index >= a && e.chunk_index < b)
                    .unwrap_or(true)
        })
        .collect();
    role_entries.sort_by_key(|e| e.chunk_index);
    for e in role_entries {
        // A paired/reference role entry is a `write_block_stream` SEGMENT (one per
        // chunk-mate-role), NOT a single bare frame like single-end. Layout:
        //   [num_blocks: u32 LE] then num_blocks × ([block_len u32][record_count u32][crc32 u32][payload])
        // Expand each inner frame into a BlockEntry pointing at its payload so the
        // single-end bounded producers (seek byte_offset, read byte_len, verify
        // expected_crc) consume it unchanged.
        if e.length < 4 {
            anyhow::bail!("paired/ref entry: segment length {} < 4", e.length);
        }
        file.seek(SeekFrom::Start(e.offset))?;
        let mut cnt = [0u8; 4];
        file.read_exact(&mut cnt)
            .with_context(|| format!("paired/ref entry (role {role:?}, mate {mate}): num_blocks"))?;
        let num_blocks = u32::from_le_bytes(cnt);
        let seg_end = e.offset + e.length;
        let mut pos = e.offset + 4; // cursor sits here (just past num_blocks)
        let mut entry_records = 0u64;
        for _ in 0..num_blocks {
            if pos + 12 > seg_end {
                anyhow::bail!("paired/ref entry: frame header exceeds segment");
            }
            let inline = crate::compression::bsc::read_block_frame_header(&mut file)
                .with_context(|| format!("paired/ref entry (role {role:?}, mate {mate}): frame header"))?;
            let payload_start = pos + 12;
            let payload_end = payload_start
                .checked_add(inline.block_len as u64)
                .ok_or_else(|| anyhow::anyhow!("paired/ref entry: frame length overflow"))?;
            if payload_end > seg_end {
                anyhow::bail!("paired/ref entry: frame payload exceeds segment");
            }
            blocks.push(BlockEntry {
                record_count: inline.record_count,
                byte_offset: payload_start,
                byte_len: inline.block_len,
                expected_crc: inline.expected_crc,
            });
            entry_records += u64::from(inline.record_count);
            file.seek(SeekFrom::Start(payload_end))?; // skip payload to next frame
            pos = payload_end;
        }
        if pos != seg_end {
            anyhow::bail!(
                "paired/ref entry: {} trailing bytes in segment",
                seg_end - pos
            );
        }
        // Footer/inline consistency: the segment's summed inner record_counts must
        // equal the footer's claimed record_count. Rejects a forged footer
        // record_count (e.g. an inflated value aimed at a giant allocation) — the
        // single-end builder makes the equivalent cross-check per bare frame.
        if entry_records != e.record_count {
            anyhow::bail!(
                "paired/ref entry: inline record_count sum {entry_records} != footer {} (role {role:?}, mate {mate})",
                e.record_count
            );
        }
        total_records += entry_records;
    }
    Ok(StreamIndex {
        blocks,
        total_records,
    })
}

/// Per-mate stream indices for the unified paired decode writer. One of these is
/// built per mate (R1 = mate 1, R2 = mate 2). `header_is_delta` selects how the
/// writer turns the header stream into ids: `false` ⇒ a normal columnar header
/// stream; `true` ⇒ `headers` is the R2 `HeaderDelta` BSC op stream, applied
/// per-record against the lockstep R1 id.
pub(super) struct PairedMateIndices {
    pub headers: StreamIndex,
    pub header_is_delta: bool,
    pub sequences: StreamIndex,
    pub qualities: Option<StreamIndex>,
    /// Quality producer codec for this mate: `Fqz` (CODEC_FQZCOMP → decoded ASCII)
    /// or `Bsc` (CODEC_BSC → bit-packed). `None` when there is no quality stream.
    pub quality_codec: Option<StreamCodec>,
}

/// Build the per-mate stream indices for one mate of a paired v5 archive.
/// Quality / rc-flags presence is derived from the directory (an absent role ⇒
/// `None`). R2 may carry either independent columnar `Headers` or a `HeaderDelta`
/// op stream; R1 always carries columnar `Headers`.
pub(super) fn build_paired_mate_indices(
    archive_path: &std::path::Path,
    dir: &crate::compression::chunk_directory::ChunkDirectory,
    mate: u8,
) -> Result<PairedMateIndices> {
    build_paired_mate_indices_range(archive_path, dir, mate, None)
}

/// Chunk-range variant of [`build_paired_mate_indices`]: each role's index is
/// restricted to chunks `[a, b)` via [`build_role_index_mate_range`]. `None` =
/// the whole archive (the `build_paired_mate_indices` delegate, byte-identical to
/// the pre-range builder). The presence/codec derivation is unchanged — it reads
/// the directory (which spans all chunks), so a mate that has a Qual role
/// anywhere still resolves the codec correctly even when slicing a subrange.
pub(super) fn build_paired_mate_indices_range(
    archive_path: &std::path::Path,
    dir: &crate::compression::chunk_directory::ChunkDirectory,
    mate: u8,
    chunk_range: Option<(u32, u32)>,
) -> Result<PairedMateIndices> {
    use crate::compression::chunk_directory::StreamRole;
    let has = |role: StreamRole| dir.entries.iter().any(|e| e.role == role && e.mate == mate);

    let (headers, header_is_delta) = if has(StreamRole::Headers) {
        (
            build_role_index_mate_range(archive_path, dir, StreamRole::Headers, mate, chunk_range)?,
            false,
        )
    } else if has(StreamRole::HeaderDelta) {
        (
            build_role_index_mate_range(
                archive_path,
                dir,
                StreamRole::HeaderDelta,
                mate,
                chunk_range,
            )?,
            true,
        )
    } else {
        anyhow::bail!("paired archive: mate {mate} has no Headers or HeaderDelta role");
    };

    let sequences =
        build_role_index_mate_range(archive_path, dir, StreamRole::Sequence, mate, chunk_range)?;

    let (qualities, quality_codec) = if has(StreamRole::Qual) {
        use crate::compression::codec_ids::{CODEC_BSC, CODEC_FQZCOMP, CODEC_QUALITY_CTX};
        let codec_byte = dir
            .entries
            .iter()
            .find(|e| e.role == StreamRole::Qual && e.mate == mate)
            .map(|e| e.codec)
            .expect("has(Qual) implies a Qual entry");
        let codec = match codec_byte {
            c if c == CODEC_FQZCOMP => StreamCodec::Fqz,
            c if c == CODEC_BSC => StreamCodec::Bsc,
            c if c == CODEC_QUALITY_CTX => anyhow::bail!(
                "paired mate {mate}: quality_ctx quality backend has been removed (code 4) — recompress"
            ),
            c => anyhow::bail!("paired mate {mate}: unknown quality codec byte {c}"),
        };
        (
            Some(build_role_index_mate_range(
                archive_path,
                dir,
                StreamRole::Qual,
                mate,
                chunk_range,
            )?),
            Some(codec),
        )
    } else {
        (None, None)
    };

    Ok(PairedMateIndices {
        headers,
        header_is_delta,
        sequences,
        qualities,
        quality_codec,
    })
}

/// Strict contract for single-end whole-archive globals: AT MOST ONE global entry,
/// and it must be the `ChunkDecodedSizes` table (mate 0, raw codec 0,
/// `record_count == num_chunks`). Any other `GLOBAL_SENTINEL` role/shape is rejected.
/// Per-chunk entries are validated by the role/chunk pairing loop in
/// [`validate_single_end_directory`].
fn validate_single_end_global_entries(dir: &crate::compression::chunk_directory::ChunkDirectory) -> Result<()> {
    use crate::compression::chunk_directory::{GLOBAL_SENTINEL, StreamRole};
    let mut globals = 0usize;
    for (i, e) in dir.entries.iter().enumerate() {
        if e.chunk_index != GLOBAL_SENTINEL { continue; }
        if e.role != StreamRole::ChunkDecodedSizes {
            anyhow::bail!("v5 single-end entry {i}: unexpected global role {:?}", e.role);
        }
        if e.mate != 0 {
            anyhow::bail!("v5 single-end entry {i}: ChunkDecodedSizes mate {} != 0", e.mate);
        }
        if e.codec != 0 {
            anyhow::bail!("v5 single-end entry {i}: ChunkDecodedSizes codec {} != 0 (raw)", e.codec);
        }
        if e.record_count != dir.num_chunks as u64 {
            anyhow::bail!(
                "v5 single-end entry {i}: ChunkDecodedSizes record_count {} != num_chunks {}",
                e.record_count, dir.num_chunks
            );
        }
        globals += 1;
        if globals > 1 {
            anyhow::bail!("v5 single-end entry {i}: more than one ChunkDecodedSizes global");
        }
    }
    Ok(())
}

/// Read a v5 single-end archive's footer and build the per-role `DecodePlan`.
/// Mode-agnostic structural validation (locator magic, footer CRC, allocation
/// bound, per-entry span/overflow + global non-overlap) is delegated to the shared
/// `chunk_directory::read_v5_footer` — the same validator the paired/reference
/// decoders use. On top of that this adds the single-end-only checks: the per-block
/// payload cap, the inline 12-byte v4 block-header cross-check, and the role/chunk
/// pairing contract; then adopts each block's CRC from its inline header.
/// Single-end-specific structural validation of a v5 directory.
///
/// Factored out of [`build_indices_from_footer`] so the WHOLE-directory contract
/// is enforced identically whether decoding the full archive or a NUMA chunk
/// range ([`build_indices_from_footer_range`]). Covers every check the
/// mode-agnostic `read_v5_footer` can't: per-entry inline-frame cross-checks +
/// payload cap, Nmask rejection, single-end mate/HeaderDelta rejection, role/chunk
/// ordering + codec uniformity, per-chunk cross-role record-count agreement, the
/// RcFlags-iff-`has_rc_canon` rule, and the non-empty headers/sequences rule.
///
/// Operates on the WHOLE directory (range builders still validate the entire
/// directory, then filter). The RcFlags/non-empty rules are re-derived directly
/// from `dir` (the directory is the source of truth, identical to the plan-based
/// checks the original computed).
pub(super) fn validate_single_end_directory(
    dir: &crate::compression::chunk_directory::ChunkDirectory,
    file: &mut std::fs::File,
    max_block_len: usize,
    has_rc_canon: bool,
) -> Result<()> {
    use std::io::{Seek, SeekFrom};

    // Single-end-specific per-entry checks the mode-agnostic validator can't do:
    // every single-end entry is a per-BLOCK frame carrying a 12-byte inline v4
    // header, so cross-check its record_count/length against that header and
    // enforce the per-block payload cap. (read_v5_footer already proved every
    // frame span lies within the payload region, so seeking to `e.offset` and
    // reading the 12-byte header is in-bounds.)
    for (i, e) in dir.entries.iter().enumerate() {
        if e.length < 12 {
            anyhow::bail!("v5 entry {i}: frame length {} < 12", e.length);
        }
        if e.record_count == 0 {
            anyhow::bail!("v5 entry {i}: record_count 0 (empty frames are never emitted)");
        }
        let payload_len = (e.length - 12) as usize;
        if payload_len > max_block_len {
            anyhow::bail!(
                "v5 entry {i}: payload {payload_len} exceeds max_block_len {max_block_len}"
            );
        }
        file.seek(SeekFrom::Start(e.offset))?;
        let inline = crate::compression::bsc::read_block_frame_header(file)
            .with_context(|| format!("v5 entry {i}: inline header"))?;
        if inline.block_len as u64 != e.length - 12 {
            anyhow::bail!(
                "v5 entry {i}: inline block_len {} != footer {}",
                inline.block_len,
                e.length - 12
            );
        }
        if u64::from(inline.record_count) != e.record_count {
            anyhow::bail!(
                "v5 entry {i}: inline record_count {} != footer {}",
                inline.record_count,
                e.record_count
            );
        }
    }

    // ── Directory-only structural contract ─────────────────────────────────
    validate_single_end_directory_structure(dir, has_rc_canon)
}

/// Directory-only single-end structural contract: role/chunk pairing, codec
/// uniformity, per-chunk record-count agreement, Headers+Sequence presence, the
/// RcFlags-iff-`has_rc_canon` rule, and the `ChunkDecodedSizes` global contract
/// (via `validate_single_end_global_entries`). Reads NOTHING from a file, so the
/// archive-merge self-check can run it on a rebuilt directory whose backing output
/// is a non-seekable stream.
pub(super) fn validate_single_end_directory_structure(
    dir: &crate::compression::chunk_directory::ChunkDirectory,
    has_rc_canon: bool,
) -> Result<()> {
    use crate::compression::chunk_directory::StreamRole;

    // ── Role/chunk pairing contract ────────────────────────────────────────
    // build_role_index preserves footer order, so the cursors pair record N's
    // header/sequence/quality by position. Without this check, an independently
    // reordered role (or a chunk with mismatched per-role counts) would silently
    // mispair records while total-record-count validation still passed.
    use std::collections::BTreeMap;
    // Single-end v5 never emits Nmask frames, and no nmask cursor is wired.
    if dir.entries.iter().any(|e| e.role == StreamRole::Nmask) {
        anyhow::bail!("v5 directory: unexpected Nmask frame");
    }
    // Enforce the strict global-entry contract: at most one ChunkDecodedSizes global,
    // with the right shape. Any other GLOBAL_SENTINEL entry is rejected here.
    validate_single_end_global_entries(dir)?;
    let nchunks = dir.num_chunks;
    let mut last_chunk: BTreeMap<u8, i64> = BTreeMap::new(); // role → last chunk_index seen
    let mut role_codec: BTreeMap<u8, u8> = BTreeMap::new(); // role → codec (must be uniform)
    let mut per_chunk: BTreeMap<(u8, u32), u64> = BTreeMap::new(); // (role,chunk) → Σ record_count
    for (i, e) in dir.entries.iter().enumerate() {
        // Whole-archive global metadata (chunk_index == GLOBAL_SENTINEL): the strict
        // contract (≤1 global, role == ChunkDecodedSizes, mate 0, codec 0,
        // record_count == num_chunks) is enforced by validate_single_end_global_entries
        // above. The continue here simply excludes globals from the per-chunk pairing
        // logic below (mirrors read_chunk_layout_impl — no decode cursor consumes globals).
        if e.chunk_index == crate::compression::chunk_directory::GLOBAL_SENTINEL {
            continue;
        }
        if e.mate != 0 {
            anyhow::bail!("v5 single-end entry {i}: unexpected mate {}", e.mate);
        }
        if e.role == StreamRole::HeaderDelta {
            anyhow::bail!("v5 single-end entry {i}: unexpected HeaderDelta role");
        }
        if e.chunk_index >= nchunks {
            anyhow::bail!(
                "v5 entry {i}: chunk_index {} >= num_chunks {}",
                e.chunk_index,
                nchunks
            );
        }
        let r = e.role as u8;
        // Footer order must equal chunk order within each role (the order the
        // cursors consume records).
        let prev = last_chunk.entry(r).or_insert(-1);
        if (e.chunk_index as i64) < *prev {
            anyhow::bail!("v5 directory: role {r} frames out of chunk order at entry {i}");
        }
        *prev = e.chunk_index as i64;
        // A role's frames must all use one codec.
        match role_codec.get(&r) {
            None => {
                role_codec.insert(r, e.codec);
            }
            Some(&c) if c == e.codec => {}
            Some(&c) => anyhow::bail!(
                "v5 entry {i}: role {r} codec {} != {c} (mixed codecs)",
                e.codec
            ),
        }
        *per_chunk.entry((r, e.chunk_index)).or_insert(0) += e.record_count;
    }
    // Mandatory streams must exist BEFORE the per-chunk pairing loop below: the loop
    // indexes `counts[0]`, and a malformed footer using only non-core roles (e.g. a forged
    // Nmask-only directory) would leave `active_roles` empty → an out-of-bounds panic
    // (a DoS on hostile input). Requiring Headers + Sequence here guarantees `active_roles`
    // is non-empty and turns the malformed case into a clean error.
    if !role_codec.contains_key(&(StreamRole::Headers as u8))
        || !role_codec.contains_key(&(StreamRole::Sequence as u8))
    {
        anyhow::bail!("v5 directory: missing headers or sequences stream");
    }

    // Active roles (present in the directory) must agree on per-chunk record
    // counts for every chunk — i.e. each chunk contributes the same number of
    // records to headers, sequences, qualities, and rc-flags.
    let active_roles: Vec<u8> = [
        StreamRole::Headers,
        StreamRole::Sequence,
        StreamRole::Qual,
        StreamRole::RcFlags,
    ]
    .iter()
    .map(|r| *r as u8)
    .filter(|&r| role_codec.contains_key(&r))
    .collect();
    for c in 0..nchunks {
        let counts: Vec<u64> = active_roles
            .iter()
            .map(|&r| per_chunk.get(&(r, c)).copied().unwrap_or(0))
            .collect();
        let first = counts[0];
        if first == 0 || counts.iter().any(|&n| n != first) {
            anyhow::bail!(
                "v5 directory: per-chunk record counts disagree across roles at chunk {c}"
            );
        }
    }

    // Role activeness: an RcFlags stream must be present iff the encoding is
    // RC-canonical (spec footer rule 6). Headers + Sequence presence was checked above.
    // Derived directly from `dir` (the directory is the source of truth — the
    // original computed these from the WHOLE-archive plan, which is equivalent).
    let has_rc_flags = role_codec.contains_key(&(StreamRole::RcFlags as u8));
    if has_rc_flags != has_rc_canon {
        anyhow::bail!(
            "v5 directory: RcFlags stream presence ({}) disagrees with encoding has_rc_canon ({})",
            has_rc_flags,
            has_rc_canon
        );
    }
    Ok(())
}

pub(super) fn build_indices_from_footer(
    archive_path: &std::path::Path,
    header_end: u64,
    max_block_len: usize,
    has_rc_canon: bool,
) -> Result<DecodePlan> {
    use crate::compression::chunk_directory::{StreamRole, read_v5_footer};

    // Shared, mode-agnostic structural validation — the SINGLE source of truth,
    // also used by the paired and reference decoders. Covers: locator magic,
    // footer_len cap + underflow, header-region overlap, footer CRC32, directory
    // parse (truncation / entry-count bound / exact consumption / role-byte
    // validity), per-entry offset+length overflow + frame-span bounds within
    // `[header_end, footer_start)`, and global frame non-overlap. Keeping these in
    // one place means a future hardening patch can't land in only one copy.
    let dir = read_v5_footer(archive_path, header_end)?;

    // Single-end-specific whole-directory structural validation.
    let mut file = std::fs::File::open(archive_path)?;
    validate_single_end_directory(&dir, &mut file, max_block_len, has_rc_canon)?;

    // Build per-role indices, adopting each block's CRC from the inline header.
    let mut plan = DecodePlan {
        headers: build_role_index(&dir, StreamRole::Headers)?,
        sequences: build_role_index(&dir, StreamRole::Sequence)?,
        qualities: {
            let q = build_role_index(&dir, StreamRole::Qual)?;
            if q.blocks.is_empty() { None } else { Some(q) }
        },
        rc_flags: {
            let r = build_role_index(&dir, StreamRole::RcFlags)?;
            if r.blocks.is_empty() { None } else { Some(r) }
        },
        num_reads: dir.num_reads,
    };
    adopt_inline_crcs(&mut file, &mut plan.headers)?;
    adopt_inline_crcs(&mut file, &mut plan.sequences)?;
    if let Some(q) = plan.qualities.as_mut() {
        adopt_inline_crcs(&mut file, q)?;
    }
    if let Some(r) = plan.rc_flags.as_mut() {
        adopt_inline_crcs(&mut file, r)?;
    }
    Ok(plan)
}

/// Like [`build_indices_from_footer`] but restricted to chunks `chunk_range`
/// (`None` = whole archive). Validates the WHOLE directory (so structural
/// corruption is caught regardless of which chunks a shard decodes), then builds
/// per-role indices filtered to the range. Used by the NUMA range decoder.
pub(super) fn build_indices_from_footer_range(
    archive_path: &std::path::Path,
    header_end: u64,
    max_block_len: usize,
    has_rc_canon: bool,
    chunk_range: Option<(u32, u32)>,
) -> Result<DecodePlan> {
    use crate::compression::chunk_directory::{StreamRole, read_v5_footer};
    let dir = read_v5_footer(archive_path, header_end)?;
    let mut file = std::fs::File::open(archive_path)?;
    validate_single_end_directory(&dir, &mut file, max_block_len, has_rc_canon)?; // whole-dir
    let mut plan = DecodePlan {
        headers: build_role_index_range(&dir, StreamRole::Headers, chunk_range)?,
        sequences: build_role_index_range(&dir, StreamRole::Sequence, chunk_range)?,
        qualities: { let q = build_role_index_range(&dir, StreamRole::Qual, chunk_range)?; if q.blocks.is_empty() { None } else { Some(q) } },
        rc_flags: { let r = build_role_index_range(&dir, StreamRole::RcFlags, chunk_range)?; if r.blocks.is_empty() { None } else { Some(r) } },
        num_reads: 0,
    };
    plan.num_reads = plan.headers.total_records;
    if plan.headers.blocks.is_empty() || plan.sequences.blocks.is_empty() {
        anyhow::bail!("v5 chunk range: empty headers or sequences");
    }
    adopt_inline_crcs(&mut file, &mut plan.headers)?;
    adopt_inline_crcs(&mut file, &mut plan.sequences)?;
    if let Some(q) = plan.qualities.as_mut() { adopt_inline_crcs(&mut file, q)?; }
    if let Some(r) = plan.rc_flags.as_mut() { adopt_inline_crcs(&mut file, r)?; }
    Ok(plan)
}

/// Fill each block's `expected_crc` from its inline v4 header (the footer does
/// not store CRCs; the inline header is the source of truth).
fn adopt_inline_crcs(file: &mut std::fs::File, idx: &mut StreamIndex) -> Result<()> {
    use std::io::{Seek, SeekFrom};
    for b in idx.blocks.iter_mut() {
        file.seek(SeekFrom::Start(b.byte_offset - 12))?;
        let hdr = crate::compression::bsc::read_block_frame_header(file)?;
        b.expected_crc = hdr.expected_crc;
    }
    Ok(())
}

// ── NUMA chunk layout ─────────────────────────────────────────────────────────

/// Per-archive chunk layout for NUMA orchestration.
///
/// Returned by [`read_chunk_layout_impl`] / `compression::read_chunk_layout`.
/// Carries only the footer metadata needed to decide how to shard decode across
/// NUMA domains — no payload bytes are read.
///
/// `shardable` is true when the archive_type is 0 (single-end), 1 (paired-end),
/// 2 (reference), or 4 (single-end reference) — chunks are independent enough to
/// assign to separate workers. Reference archives (types 2 and 4) carry
/// whole-archive globals (PackedBacking/IntervalMap/NBitmap/ReferenceMeta) that
/// must be loaded once per worker; their resident size is folded into
/// `reference_resident_bytes`. The NUMA driver gates reference sharding on its
/// `INCREMENT_B_READY`, not on this flag (only the mixed-codec legacy paired path
/// is `shardable == false`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChunkLayout {
    /// On-disk `archive_type` byte: 0=single-end, 1=paired-end, 2=reference,
    /// 4=single-end reference (one output, no ChunkDecodedSizes table).
    pub archive_type: u8,
    /// On-disk `encoding_type` byte (e.g. 0=Raw, 10=UltraBigBlock).
    pub encoding_type: u8,
    /// Number of chunks in the archive (from the footer).
    pub num_chunks: u32,
    /// Total reads across all chunks (from the footer).
    pub num_reads: u64,
    /// Per-chunk read count, indexed by chunk_index (ascending).
    /// For single-end archives each entry is reads in that chunk.
    /// For paired-end archives each entry counts PAIRS (R1 == R2 record count).
    pub per_chunk_reads: Vec<u64>,
    /// True when the archive can be range-decoded by independent workers. Only the
    /// mixed-codec legacy paired path (mate-2 mixing columnar `Headers` and
    /// `HeaderDelta` across chunks) is unshardable; single-end, paired, reference,
    /// and single-end reference (type 4) are all shardable. (Reference sharding is
    /// separately gated by the driver's `INCREMENT_B_READY`, not by this flag.)
    pub shardable: bool,
    /// Estimated RESIDENT (decoded) reference-globals bytes that every worker loads
    /// before decoding — derived from the `PackedBacking` global entry's
    /// `record_count` (≈ 3/8 B/base for the 2-bit consensus + N-bitmap, plus
    /// IntervalMap/meta headroom), NOT from compressed lengths. Zero for single-end
    /// and paired-end; non-zero for reference (type 2) and single-end reference
    /// (type 4). Folded into the gate's per-node memory check.
    pub reference_resident_bytes: u64,

    /// Per-chunk decoded OUTPUT byte length, indexed by chunk_index, from the
    /// `ChunkDecodedSizes` global block. EMPTY when the archive predates the table
    /// (the NUMA driver then falls back to part-file assembly). For paired this is
    /// the per-chunk SUM across mates.
    pub per_chunk_decoded_bytes: Vec<u64>,

    /// Per-MATE, per-chunk decoded OUTPUT byte length: `decoded_sizes_per_mate[mate]`
    /// is the per-chunk vector for that mate (outer len = n_mates: 1 for single-end /
    /// single-end reference, 2 for paired R1/R2). This is what direct-write needs —
    /// each mate has its OWN pre-sized output file, so the per-mate prefix sums give
    /// the byte offset where each worker's chunk-range begins. EMPTY (outer len 0)
    /// when the archive predates the `ChunkDecodedSizes` table, exactly mirroring
    /// `per_chunk_decoded_bytes`. The combined `per_chunk_decoded_bytes[c]` equals
    /// `Σ_mate decoded_sizes_per_mate[mate][c]`.
    pub decoded_sizes_per_mate: Vec<Vec<u64>>,
}

/// Footer-only chunk layout query — no payload decompression.
///
/// Reads the v5 archive header and footer (no block payloads), runs the
/// appropriate structural validator (`validate_paired_directory` or
/// `validate_reference_directory`), and returns a [`ChunkLayout`].
pub(super) fn read_chunk_layout_impl(archive: &std::path::Path) -> Result<ChunkLayout> {
    use chunk_directory::{GLOBAL_SENTINEL, StreamRole, read_v5_footer};

    // Read enough of the file to parse the fixed front header.
    // The header is at most V2_PREFIX_SIZE + 12 + 8 = 28 bytes; read 64 to be safe.
    let mut file = std::fs::File::open(archive)
        .with_context(|| format!("read_chunk_layout: open {:?}", archive))?;
    let mut hbuf = [0u8; 64];
    let n = {
        use std::io::Read;
        file.read(&mut hbuf)?
    };
    let hbuf = &hbuf[..n];

    let h = FixedHeader::parse_v5(hbuf)?;
    let archive_type = h.archive_type;

    // `header_size` is stored in the v2 prefix at bytes 4..8 (little-endian u32).
    // Same as the pattern used in the streaming decoders (`read_le_u32(&hbuf, 4)`).
    let header_end = read_le_u32(hbuf, 4)? as u64;

    // Parse the footer — runs the generic span/overlap/bounds validation.
    let dir = read_v5_footer(archive, header_end)?;

    // Run the per-mode structural validator.
    match archive_type {
        1 => {
            paired::format::validate_paired_directory(&dir)?;
        }
        2 => {
            reference::validate_reference_directory(&dir)?;
        }
        4 => {
            // Single-end reference: per-chunk mate-1-only contract (no mate-2 entries).
            reference::validate_reference_directory_single(&dir)?;
        }
        // Single-end: run the cheap STRUCTURAL validator (role/chunk-pairing contract,
        // no file I/O) for parity with types 1/2/4, so a forged footer is rejected before
        // the NUMA planner allocates worker ranges from it. The full single-end contract
        // — `validate_single_end_directory` + `validate_single_end_global_entries`, which
        // read the file to CRC-check globals — still runs on the actual decode path
        // (build_indices_from_footer); the one global read below (ChunkDecodedSizes) is
        // CRC- and shape-validated inline.
        0 => {
            validate_single_end_directory_structure(&dir, h.has_rc_canon())?;
        }
        _ => {}
    }

    // DoS guard: num_chunks is an unbounded u32 in the footer; every chunk has >=1
    // entry, so this holds for any real archive and bounds the vec alloc below.
    if dir.num_chunks as usize > dir.entries.len() {
        anyhow::bail!(
            "v5 footer: num_chunks {} > entries {} (forged footer)",
            dir.num_chunks, dir.entries.len()
        );
    }

    // Cluster archives (archive_type 3) decode through a bespoke, order-dropping path
    // (`cluster::decompress_cluster`) whose per-read data is bucket-bound, not
    // chunk-range addressable — they are NOT shardable. Report a valid non-shardable
    // layout (rather than bailing on the unrecognized type below) so the NUMA driver
    // routes them to in-process decode instead of crashing under `--numa N`. Per-chunk
    // tables are left empty: a non-shardable layout never feeds the partitioner.
    if archive_type == 3 {
        return Ok(ChunkLayout {
            archive_type,
            encoding_type: h.encoding_type,
            num_chunks: dir.num_chunks,
            num_reads: dir.num_reads,
            per_chunk_reads: Vec::new(),
            shardable: false,
            reference_resident_bytes: 0,
            per_chunk_decoded_bytes: Vec::new(),
            decoded_sizes_per_mate: Vec::new(),
        });
    }

    // Derive per-chunk read counts from the footer entries.
    // For single-end: use the Headers role (one entry per chunk, mate=0).
    // For paired-end: use mate=1 Headers (R1; record count == pair count).
    // For reference: use mate=1 MappedFlags (R1); global-sentinel entries are skipped.
    //
    // The footer `num_reads` is authoritative for the total; per_chunk_reads is
    // assembled from the per-chunk entries for the canonical role so callers can
    // assign chunk ranges to workers.
    let canonical = match archive_type {
        0 | 1 => StreamRole::Headers,
        2 | 4 => StreamRole::MappedFlags,
        other => anyhow::bail!("unknown archive_type {other}"),
    };
    let canonical_mate: u8 = if archive_type == 0 { 0 } else { 1 };

    let n = dir.num_chunks as usize;
    let mut per_chunk_reads = vec![0u64; n];
    let mut reference_resident_bytes = 0u64;
    let mut decoded_sizes_span: Option<(u64, u64)> = None;
    for e in &dir.entries {
        if e.chunk_index == GLOBAL_SENTINEL {
            if super::is_v5_reference(archive_type) && e.role == StreamRole::PackedBacking {
                reference_resident_bytes = reference_resident_bytes
                    .saturating_add(e.record_count.saturating_mul(3) / 8)
                    .saturating_add(e.record_count / 4);
            }
            if e.role == StreamRole::ChunkDecodedSizes {
                decoded_sizes_span = Some((e.offset, e.length));
            }
            continue;
        }
        if e.role == canonical && e.mate == canonical_mate && (e.chunk_index as usize) < n {
            per_chunk_reads[e.chunk_index as usize] = per_chunk_reads[e.chunk_index as usize]
                .checked_add(e.record_count)
                .ok_or_else(|| anyhow::anyhow!("per-chunk read count overflow at chunk {}", e.chunk_index))?;
        }
    }

    let total: u64 = per_chunk_reads.iter().try_fold(0u64, |a, &c| {
        if c == 0 { anyhow::bail!("v5 footer: a chunk has 0 canonical reads"); }
        a.checked_add(c).ok_or_else(|| anyhow::anyhow!("total read count overflow"))
    })?;
    if total != dir.num_reads {
        anyhow::bail!("v5 footer: Σ per-chunk canonical reads {total} != num_reads {}", dir.num_reads);
    }

    let (per_chunk_decoded_bytes, decoded_sizes_per_mate) = match decoded_sizes_span {
        Some((off, len)) if len > 12 => {
            use std::io::{Read, Seek, SeekFrom};
            let mut f = std::fs::File::open(archive)?;
            f.seek(SeekFrom::Start(off))?;
            let mut frame = vec![0u8; len as usize];
            f.read_exact(&mut frame)?;
            let block_len = u32::from_le_bytes(frame[0..4].try_into().unwrap()) as usize;
            let record_count = u32::from_le_bytes(frame[4..8].try_into().unwrap());
            let crc = u32::from_le_bytes(frame[8..12].try_into().unwrap());
            let payload = &frame[12..];
            if payload.len() != block_len {
                anyhow::bail!(
                    "ChunkDecodedSizes frame: block_len field {block_len} != actual payload bytes {}",
                    payload.len()
                );
            }
            // Canonical frame-CRC check (shared with every other v5 block) — preserves
            // stored/computed hex in the error message for corruption triage.
            super::bsc::verify_block_frame_crc(crc, record_count, payload)?;
            let (n_mates, num_chunks, sizes) =
                crate::compression::chunk_directory::parse_decoded_sizes(payload)?;
            if num_chunks != dir.num_chunks {
                anyhow::bail!("ChunkDecodedSizes num_chunks {num_chunks} != footer {}", dir.num_chunks);
            }
            // The frame's record_count and the payload's num_chunks are both written as
            // the chunk count; require agreement (the CRC covers both, so this only trips
            // on an inconsistent re-frame — a cheap explicit invariant guard).
            if record_count != num_chunks {
                anyhow::bail!(
                    "ChunkDecodedSizes: frame record_count {record_count} != payload num_chunks {num_chunks}"
                );
            }
            let mut per_chunk = vec![0u64; num_chunks as usize];
            // Per-mate transpose: row-major sizes[c*n_mates + m] -> per_mate[m][c].
            let mut per_mate = vec![vec![0u64; num_chunks as usize]; n_mates as usize];
            for c in 0..num_chunks as usize {
                for m in 0..n_mates as usize {
                    let v = sizes[c * n_mates as usize + m];
                    per_mate[m][c] = v;
                    per_chunk[c] = per_chunk[c]
                        .checked_add(v)
                        .ok_or_else(|| anyhow::anyhow!("ChunkDecodedSizes per-chunk overflow"))?;
                }
            }
            (per_chunk, per_mate)
        }
        _ => (Vec::new(), Vec::new()),
    };

    let shardable = !(archive_type == 1 && paired_uses_mixed_r2_header(&dir));

    Ok(ChunkLayout {
        archive_type: h.archive_type,
        encoding_type: h.encoding_type,
        num_chunks: dir.num_chunks,
        num_reads: dir.num_reads,
        per_chunk_reads,
        shardable,
        reference_resident_bytes,
        per_chunk_decoded_bytes,
        decoded_sizes_per_mate,
    })
}

/// True iff mate-2 mixes columnar `Headers` and `HeaderDelta` representations
/// across chunks (the legacy `decode_chunk_v5` serial path).
fn paired_uses_mixed_r2_header(dir: &chunk_directory::ChunkDirectory) -> bool {
    use chunk_directory::StreamRole;
    let has_cols = dir.entries.iter().any(|e| e.mate == 2 && e.role == StreamRole::Headers);
    let has_delta = dir.entries.iter().any(|e| e.mate == 2 && e.role == StreamRole::HeaderDelta);
    has_cols && has_delta
}

/// Conservative per-worker peak-RSS estimate (bytes) for the NUMA gate's
/// memory check.
///
/// Returns an upper bound on how much memory one decode worker will use for
/// the largest block it could encounter, based on the archive's encoding and
/// type. This is intentionally conservative: it covers the decompressed block
/// + the in-flight channel buffer, not the steady-state average.
///
/// - Ultra archives (encoding_type 10) use 768 MiB blocks → ~3 × 768 MiB at peak.
/// - Default archives use 64 MiB blocks → ~3 × 64 MiB at peak.
/// - Reference archives add `reference_resident_bytes` on top (the whole-archive
///   globals must be resident for every worker).
///
/// The multiplier of 3 comes from `(BOUNDED_CHANNEL_CAP + 1)` blocks per stream
/// (channel holds up to 2 + the one currently being consumed).
pub(super) fn decode_peak_rss_bound_impl(
    encoding_type: u8,
    archive_type: u8,
    reference_resident_bytes: u64,
) -> u64 {
    let (max_block_len, max_inflight) = decode_block_bounds(encoding_type);
    let per_stream = (BOUNDED_CHANNEL_CAP as u64 + max_inflight as u64) * max_block_len as u64;
    let streams: u64 = match archive_type { 0 => 4, 1 | 2 | 4 => 6, _ => 6 };
    per_stream.saturating_mul(streams).saturating_add(reference_resident_bytes)
}

// ── NUMA chunk-range decode ─────────────────────────────────────────────────

/// A `Write` adaptor that counts bytes and hard-caps the total at `limit`. Used by
/// the direct-write path so a `ChunkDecodedSizes` UNDERESTIMATE (the decoded range
/// is bigger than the reserved region) is caught as a CLEAN error before it can
/// overwrite the next worker's region. On the first write that would exceed `limit`
/// it sets `exceeded` and returns an error (the caller converts that into a
/// `DirectWriteIntegrityError`), writing nothing past the limit.
pub(crate) struct BoundedCountingWriter<W: std::io::Write> {
    inner: W,
    pub written: u64,
    limit: u64,
    pub exceeded: bool,
}

impl<W: std::io::Write> BoundedCountingWriter<W> {
    pub(crate) fn new(inner: W, limit: u64) -> Self {
        Self { inner, written: 0, limit, exceeded: false }
    }
}

impl<W: std::io::Write> std::io::Write for BoundedCountingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let would = match self.written.checked_add(buf.len() as u64) {
            Some(w) => w,
            None => {
                self.exceeded = true;
                return Err(std::io::Error::other("direct-write count overflow"));
            }
        };
        if would > self.limit {
            self.exceeded = true; // caller converts this to DirectWriteIntegrityError
            return Err(std::io::Error::other(format!(
                "direct-write exceeds region limit {} (would reach {would})",
                self.limit
            )));
        }
        let n = self.inner.write(buf)?;
        self.written += n as u64;
        Ok(n)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// Decode only chunks `[chunk_start, chunk_end)` of a v5 archive into
/// `out_parts`. Routes by `archive_type`: single-end (1 part), paired-end (2 parts,
/// R1/R2), reference (2 parts, R1/R2 — via `reference::decode_reference_range`), and
/// single-end reference (type 4: 1 part — via `reference::decode_reference_range_single`).
///
/// A non-empty `regions` requests the verified direct-write path — one
/// [`DirectWriteRegion`] per output (single-end: 1; paired: 2). Direct-write is
/// rejected for gzip, reference, and single-end reference; empty `regions` is the
/// part-file path used by every mode.
#[allow(clippy::too_many_arguments)]
pub(super) fn decode_chunk_range_impl(
    archive: &std::path::Path,
    chunk_start: u32,
    chunk_end: u32,
    threads: usize,
    gzipped: bool,
    gzip_level: u32,
    out_parts: &[std::path::PathBuf],
    regions: &[crate::compression::DirectWriteRegion],
) -> Result<()> {
    if !regions.is_empty() && gzipped {
        return Err(crate::compression::DirectWriteIntegrityError(
            "gzip output is incompatible with direct-write".into(),
        )
        .into());
    }
    let _ = super::ensure_global_thread_pool(threads);
    let layout = read_chunk_layout_impl(archive)?;
    if chunk_start > chunk_end || chunk_end > layout.num_chunks {
        anyhow::bail!("chunk range [{chunk_start},{chunk_end}) out of bounds (num_chunks={})", layout.num_chunks);
    }
    if gzipped && layout.archive_type != 0 {
        anyhow::bail!("--gzipped is only supported for single-end archives");
    }
    match layout.archive_type {
        0 => decode_chunk_range_single(archive, chunk_start, chunk_end, threads, gzipped, gzip_level, out_parts, regions),
        1 => decode_chunk_range_paired(archive, chunk_start, chunk_end, out_parts, regions),
        2 => {
            if out_parts.len() != 2 {
                anyhow::bail!("reference decode_chunk_range expects 2 output parts, got {}", out_parts.len());
            }
            crate::compression::reference::decode_reference_range(archive, chunk_start, chunk_end, out_parts, regions)
        }
        4 => {
            if out_parts.len() != 1 {
                anyhow::bail!("single-end reference decode_chunk_range expects 1 output part, got {}", out_parts.len());
            }
            crate::compression::reference::decode_reference_range_single(archive, chunk_start, chunk_end, out_parts, regions)
        }
        other => anyhow::bail!("unknown archive_type {other}"),
    }
}

#[allow(clippy::too_many_arguments)]
fn decode_chunk_range_single(
    archive: &std::path::Path,
    a: u32,
    b: u32,
    threads: usize,
    gzipped: bool,
    gzip_level: u32,
    out_parts: &[std::path::PathBuf],
    regions: &[crate::compression::DirectWriteRegion],
) -> Result<()> {
    if out_parts.len() != 1 {
        anyhow::bail!("single-end decode_chunk_range expects 1 output part, got {}", out_parts.len());
    }
    if !regions.is_empty() {
        if regions.len() != 1 {
            return Err(crate::compression::DirectWriteIntegrityError(format!(
                "single-end direct-write expects 1 region, got {}",
                regions.len()
            ))
            .into());
        }
        return decode_chunk_range_single_direct(archive, a, b, gzipped, gzip_level, &out_parts[0], regions[0]);
    }
    // Defense-in-depth: the full single-end decode entry gates on this (rejects
    // legacy/unstreamable encodings + unsupported codecs); the NUMA path must too,
    // so a forged or out-of-contract single-end archive can't slip through the
    // chunk-range decoder. Single-end only — paired has validate_paired_directory.
    require_streamable_v5(archive)?;
    let meta = read_single_end_stream_meta(archive)?;
    let (max_block_len, max_inflight) = decode_block_bounds(meta.encoding_type);
    let plan = build_indices_from_footer_range(archive, meta.header_end, max_block_len, meta.has_rc_canon, Some((a, b)))?;
    let num_reads = plan.num_reads;

    let file = create_new_for_write(&out_parts[0])?;
    if gzipped {
        use gzp::{Compression, ZWriter};
        use gzp::deflate::Gzip;
        use gzp::par::compress::ParCompressBuilder;
        let num_gz_threads = (threads / 2).max(2);
        let mut out: gzp::par::compress::ParCompress<Gzip> = ParCompressBuilder::new()
            .num_threads(num_gz_threads)
            .map_err(|e| anyhow::anyhow!("gzp error: {e}"))?
            .compression_level(Compression::new(gzip_level))
            .from_writer(file);
        run_single_writer(&mut out, num_reads, archive, &plan, &meta, max_inflight, max_block_len)?;
        out.finish().map_err(|e| anyhow::anyhow!("gzp finish error: {e}"))?;
    } else {
        let mut out = std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, file);
        run_single_writer(&mut out, num_reads, archive, &plan, &meta, max_inflight, max_block_len)?;
        use std::io::Write;
        out.flush()?;
    }
    Ok(())
}

/// Verify a worker-supplied direct-write `region` against the authoritative
/// per-mate decoded-size vector `sizes` (one entry per chunk, ascending), returning
/// the table-derived expected length for chunks `[a, b)`. Shared by the single-end
/// and paired direct-write paths so the "never silent corruption, always a clean
/// integrity error" contract is identical for every mate/output:
///
/// * `region.offset` must equal `Σ sizes[0..a]` (the byte where this shard begins),
/// * `region.len` must equal `Σ sizes[a..b]` (this shard's exact decoded length).
///
/// Any out-of-bounds range, arithmetic overflow, or mismatch is a
/// [`crate::compression::DirectWriteIntegrityError`] — so `auto` falls back to a
/// correct in-process decode (which ignores the table entirely).
pub(crate) fn verify_direct_region(
    sizes: &[u64],
    a: u32,
    b: u32,
    region: crate::compression::DirectWriteRegion,
) -> Result<u64> {
    let ierr = |m: String| -> anyhow::Error { crate::compression::DirectWriteIntegrityError(m).into() };
    let n = sizes.len() as u32;
    if a > b || b > n {
        return Err(ierr(format!("chunk range [{a},{b}) out of bounds (n={n})")));
    }
    let mut expected_offset: u64 = 0;
    for &s in &sizes[..a as usize] {
        expected_offset = expected_offset.checked_add(s).ok_or_else(|| ierr("offset overflow".into()))?;
    }
    let mut expected_len: u64 = 0;
    for &s in &sizes[a as usize..b as usize] {
        expected_len = expected_len.checked_add(s).ok_or_else(|| ierr("len overflow".into()))?;
    }
    if region.offset != expected_offset {
        return Err(ierr(format!("base {} != table {}", region.offset, expected_offset)));
    }
    if region.len != expected_len {
        return Err(ierr(format!("len {} != table {}", region.len, expected_len)));
    }
    Ok(expected_len)
}

/// Verified direct-write path for single-end chunk-range decode.
///
/// Writes chunks `[a, b)` straight into the pre-sized shared file `out_path` at
/// `region.offset`, never creating a part file. The integrity contract is the whole
/// point: the supplied `offset`/`len` are checked against the authoritative
/// `ChunkDecodedSizes` table (a wrong size table is a CLEAN failure, never silent
/// corruption), the decode write is hard-capped at `region.len` via
/// `BoundedCountingWriter` (catches a table UNDERestimate as an overrun), and a
/// short write (table OVERestimate) is caught at the end. Every failure is a
/// [`crate::compression::DirectWriteIntegrityError`] so `auto` can clean up and fall
/// back in-process. Gzip is rejected upstream (incompatible with offset writes).
fn decode_chunk_range_single_direct(
    archive: &std::path::Path,
    a: u32,
    b: u32,
    gzipped: bool,
    gzip_level: u32,
    out_path: &std::path::Path,
    region: crate::compression::DirectWriteRegion,
) -> Result<()> {
    use std::io::{Seek, SeekFrom, Write};
    let _ = (gzipped, gzip_level); // gzip+region rejected upstream
    let ierr = |m: String| -> anyhow::Error { crate::compression::DirectWriteIntegrityError(m).into() };

    // 1) Verify the supplied region against the authoritative table (checked).
    // A corrupt/CRC-failed ChunkDecodedSizes frame is itself a table problem, so
    // classify a layout-read failure as an integrity error (not a hard worker
    // failure) — `auto` then falls back to in-process decode, which ignores the
    // table and produces correct output. Upholds the "never silent corruption,
    // always clean fallback" contract for the present-but-corrupt-table case.
    let layout = read_chunk_layout_impl(archive).map_err(|e| ierr(format!("layout read failed: {e}")))?;
    if layout.per_chunk_decoded_bytes.is_empty() {
        return Err(ierr("archive has no ChunkDecodedSizes table".into()));
    }
    // Single-end: the lone mate's per-chunk sizes are the combined per-chunk sizes.
    let expected_len = verify_direct_region(&layout.per_chunk_decoded_bytes, a, b, region)?;

    // 2) Same plan/meta/range setup as the part-file (None) path — including the
    //    defense-in-depth streamable-v5 gate.
    require_streamable_v5(archive)?;
    let meta = read_single_end_stream_meta(archive)?;
    let (max_block_len, max_inflight) = decode_block_bounds(meta.encoding_type);
    let plan = build_indices_from_footer_range(archive, meta.header_end, max_block_len, meta.has_rc_canon, Some((a, b)))?;
    let num_reads = plan.num_reads;

    // 3) Seek to base, wrap in a bounded sink, decode; classify a confirmed overrun
    //    (table underestimate) and a short write (table overestimate) as integrity.
    let mut f = std::fs::OpenOptions::new().write(true).open(out_path)?;
    f.seek(SeekFrom::Start(region.offset))?;
    let mut out = BoundedCountingWriter::new(std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, f), expected_len);
    let res = run_single_writer(&mut out, num_reads, archive, &plan, &meta, max_inflight, max_block_len);
    if out.exceeded {
        return Err(ierr(format!(
            "decoded output exceeded region {expected_len} (ChunkDecodedSizes underestimate)"
        )));
    }
    res?;
    out.flush()?;
    if out.written != expected_len {
        return Err(ierr(format!("wrote {} != region {}", out.written, expected_len)));
    }
    Ok(())
}

/// Decode chunks `[a, b)` of a paired v5 archive into `out_parts[0]` (R1) and
/// `out_parts[1]` (R2). Mirrors `decompress_paired_v5`'s front-header read +
/// footer + `validate_paired_directory` gate, slices the per-mate indices to the
/// chunk range, and reuses the same `bounded_write_records_paired_from_indices`
/// tail the full paired path uses — so concatenating the parts byte-equals a full
/// decode. The mixed-codec legacy path (per-chunk R2-header codec) can't be
/// streamed, so a mixed archive is rejected (the NUMA gate's `shardable` flag
/// already reports this; this is the defense-in-depth bail).
///
/// `regions` selects the output discipline (the paired analogue of single-end's
/// `Option<DirectWriteRegion>`):
/// * empty ⇒ **part-file** path: each mate's output is created from byte 0 and the
///   driver concatenates the parts.
/// * `[r1, r2]` ⇒ **direct-write** path: each mate's decoded bytes are written
///   straight into its pre-sized shared output at `r.offset`, after verifying
///   `r.offset`/`r.len` against that mate's `decoded_sizes_per_mate` slice and
///   bounding the write to `r.len`. Any mismatch/overrun/short write is a
///   [`crate::compression::DirectWriteIntegrityError`] (so `auto` falls back to a
///   correct in-process decode). Both writers share the same setup as the part-file
///   path — only the sink differs — so a direct decode is byte-identical to a
///   part-file decode plus concat.
fn decode_chunk_range_paired(
    archive: &std::path::Path,
    a: u32,
    b: u32,
    out_parts: &[std::path::PathBuf],
    regions: &[crate::compression::DirectWriteRegion],
) -> Result<()> {
    use crate::compression::chunk_directory::read_v5_footer;
    use std::io::Write;
    if out_parts.len() != 2 {
        anyhow::bail!("paired decode_chunk_range expects 2 output parts, got {}", out_parts.len());
    }
    let direct = !regions.is_empty();
    if direct && regions.len() != 2 {
        // A parent-side region-assembly bug — classify as integrity so `auto` recovers.
        return Err(crate::compression::DirectWriteIntegrityError(format!(
            "paired direct-write expects 2 regions, got {}",
            regions.len()
        ))
        .into());
    }

    // ── Shared setup (identical for part-file and direct): header, footer, validate,
    //    mixed-codec reject, per-mate range indices. ──
    let (hdr, header_end) = read_fixed_header(archive)?;
    if hdr.archive_type != 1 {
        anyhow::bail!("not a paired v5 archive (archive_type {})", hdr.archive_type);
    }
    if hdr.encoding_type != 0 {
        anyhow::bail!(
            "paired v5 archive has unexpected encoding_type {} (paired is always 0/Raw)",
            hdr.encoding_type
        );
    }
    let dir = read_v5_footer(archive, header_end)?;
    crate::compression::paired::format::validate_paired_directory(&dir)?;
    if paired_uses_mixed_r2_header(&dir) {
        anyhow::bail!("paired archive uses the mixed-codec legacy path; --numa unsupported (use --numa off)");
    }
    let r1 = build_paired_mate_indices_range(archive, &dir, 1, Some((a, b)))?;
    let r2 = build_paired_mate_indices_range(archive, &dir, 2, Some((a, b)))?;

    if !direct {
        // Part-file path: each mate written from byte 0; the driver concatenates.
        let f1 = create_new_for_write(&out_parts[0])?;
        let f2 = create_new_for_write(&out_parts[1])?;
        let mut w1 = std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, f1);
        let mut w2 = std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, f2);
        bounded_write_records_paired_from_indices(
            archive,
            &hdr,
            &r1,
            &r2,
            &mut PairedSink::Two(&mut w1, &mut w2),
        )?;
        w1.flush()?;
        w2.flush()?;
        return Ok(());
    }

    // Direct-write path: verify each mate's region against its per-mate decoded-size
    // slice, seek into the two pre-sized shared outputs, decode through bounded sinks.
    use std::io::{Seek, SeekFrom};
    let ierr = |m: String| -> anyhow::Error { crate::compression::DirectWriteIntegrityError(m).into() };
    let layout = read_chunk_layout_impl(archive).map_err(|e| ierr(format!("layout read failed: {e}")))?;
    if layout.decoded_sizes_per_mate.len() != 2 {
        return Err(ierr(format!(
            "paired archive has no 2-mate ChunkDecodedSizes table (mates={})",
            layout.decoded_sizes_per_mate.len()
        )));
    }
    let exp1 = verify_direct_region(&layout.decoded_sizes_per_mate[0], a, b, regions[0])?;
    let exp2 = verify_direct_region(&layout.decoded_sizes_per_mate[1], a, b, regions[1])?;

    let mut f1 = std::fs::OpenOptions::new().write(true).open(&out_parts[0])?;
    let mut f2 = std::fs::OpenOptions::new().write(true).open(&out_parts[1])?;
    f1.seek(SeekFrom::Start(regions[0].offset))?;
    f2.seek(SeekFrom::Start(regions[1].offset))?;
    let mut w1 = BoundedCountingWriter::new(std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, f1), exp1);
    let mut w2 = BoundedCountingWriter::new(std::io::BufWriter::with_capacity(IO_BUFFER_SIZE, f2), exp2);
    let res = bounded_write_records_paired_from_indices(
        archive,
        &hdr,
        &r1,
        &r2,
        &mut PairedSink::Two(&mut w1, &mut w2),
    );
    // Classify a confirmed overrun (table underestimate) before the raw io error.
    if w1.exceeded {
        return Err(ierr(format!(
            "R1 decoded output exceeded region {exp1} (ChunkDecodedSizes underestimate)"
        )));
    }
    if w2.exceeded {
        return Err(ierr(format!(
            "R2 decoded output exceeded region {exp2} (ChunkDecodedSizes underestimate)"
        )));
    }
    res?;
    w1.flush()?;
    w2.flush()?;
    // A short write is a table overestimate.
    if w1.written != exp1 {
        return Err(ierr(format!("R1 wrote {} != region {}", w1.written, exp1)));
    }
    if w2.written != exp2 {
        return Err(ierr(format!("R2 wrote {} != region {}", w2.written, exp2)));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_single_writer(
    out: &mut dyn std::io::Write,
    num_reads: u64,
    archive: &std::path::Path,
    plan: &DecodePlan,
    meta: &SingleEndStreamMeta,
    max_inflight: usize,
    max_block_len: usize,
) -> Result<()> {
    bounded_write_records_independent(
        out, num_reads, archive,
        &plan.headers, &plan.sequences, plan.qualities.as_ref(), plan.rc_flags.as_ref(),
        meta.header_codec, meta.sequence_codec, meta.quality_codec,
        meta.const_seq_len, meta.const_qual_len, meta.has_sequence_hints,
        meta.bits_per_qual, meta.quality_binning, meta.is_fasta,
        max_inflight, max_block_len,
    )
}

#[cfg(test)]
mod tests_bounded {
    use super::*;

    // ── BoundedCountingWriter overrun hardening ──────────────────────────

    #[test]
    fn bounded_counting_writer_rejects_overrun() {
        use std::io::Write;
        let mut buf = Vec::new();
        let mut w = BoundedCountingWriter::new(&mut buf, 4);
        assert!(w.write_all(b"ABC").is_ok());
        assert!(w.write_all(b"DE").is_err()); // would reach 5 > limit 4
        assert!(w.exceeded);
    }

    // ── StreamCursor malformed-length hardening ──────────────────────────

    /// LEB128 encoding of `usize::MAX` (10 bytes: nine 0xFF + 0x01). Used to
    /// drive the `block.pos + len` overflow in the streaming record readers.
    fn varint_usize_max() -> Vec<u8> {
        let mut v = vec![0xFFu8; 9];
        v.push(0x01);
        v
    }

    #[test]
    fn read_one_record_varint_rejects_overflowing_length() {
        // A record whose varint length is usize::MAX must produce a clean error,
        // not an arithmetic overflow / slice-bounds panic (block.pos + len wraps).
        let (tx, rx) = std::sync::mpsc::channel();
        tx.send(Ok((1u32, varint_usize_max()))).unwrap();
        drop(tx);
        let mut cur = StreamCursor::new(rx, "test");
        assert!(cur.read_one_record_varint().is_err());
    }

    #[test]
    fn read_one_record_seq_with_hint_rejects_overflowing_length() {
        let (tx, rx) = std::sync::mpsc::channel();
        tx.send(Ok((1u32, varint_usize_max()))).unwrap();
        drop(tx);
        let mut cur = StreamCursor::new(rx, "test");
        assert!(cur.read_one_record_seq_with_hint(true).is_err());
    }

    // ── assemble_fastq_record unit tests ─────────────────────────────────

    #[test]
    fn test_assemble_no_quality() {
        let mut out = Vec::new();
        assemble_fastq_record(b"@read1", b"ACGT", 0, QualityInput::None, false, &mut out).unwrap();
        // Discarded-quality archives still emit the 4-line FASTQ shape
        // with empty quality.
        assert_eq!(&out, b"@read1\nACGT\n+\n\n");
    }

    #[test]
    fn test_assemble_fasta() {
        // With is_fasta, emit FASTA: header ('>' prefix) + sequence, no +/quality.
        let mut out = Vec::new();
        assemble_fastq_record(b">seq1", b"ACGT", 0, QualityInput::None, true, &mut out).unwrap();
        assert_eq!(&out, b">seq1\nACGT\n");
    }

    #[test]
    fn test_assemble_decoded_quality_no_rc() {
        let mut out = Vec::new();
        assemble_fastq_record(
            b"@read1",
            b"ACGT",
            0,
            QualityInput::Decoded(b"IIII"),
            false,
            &mut out,
        )
        .unwrap();
        assert_eq!(&out, b"@read1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn test_assemble_decoded_quality_with_rc() {
        let mut out = Vec::new();
        assemble_fastq_record(
            b"@read2",
            b"ACGT",
            1, // RC flag set → emit reverse-complement
            QualityInput::Decoded(b"IIII"),
            false,
            &mut out,
        )
        .unwrap();
        // reverse_complement("ACGT") = "ACGT" — pick an asymmetric seq for clarity.
        // Quality is NOT reversed (stored aligned to written orientation).
        assert_eq!(&out, b"@read2\nACGT\n+\nIIII\n");

        // Now a sequence whose RC differs.
        let mut out2 = Vec::new();
        assemble_fastq_record(
            b"@read3",
            b"AAGG",
            1,
            QualityInput::Decoded(b"abcd"),
            false,
            &mut out2,
        )
        .unwrap();
        // rc(AAGG) = CCTT.
        assert_eq!(&out2, b"@read3\nCCTT\n+\nabcd\n");
    }

    #[test]
    fn test_assemble_packed_quality_roundtrip() {
        // Pack a known quality string with QualityBinning::None (7 bits) and
        // verify assemble_fastq_record reproduces the original ASCII.
        let qual_ascii: &[u8] = b"!ABCDIJ"; // Phred+33: 0, 32, 33, 34, 35, 40, 41
        let binning = QualityBinning::None;
        let packed = columnar::pack_qualities(qual_ascii, binning).unwrap();

        let mut out = Vec::new();
        assemble_fastq_record(
            b"@r",
            b"ACGTACG",
            0,
            QualityInput::Packed {
                l_seq: qual_ascii.len(),
                packed: &packed,
                binning,
            },
            false,
            &mut out,
        )
        .unwrap();
        // QualityBinning::None preserves the 7 low bits, so the unpack
        // produces the original ASCII bytes back.
        assert_eq!(&out, b"@r\nACGTACG\n+\n!ABCDIJ\n");
    }

    #[test]
    fn test_assemble_packed_quality_const_length() {
        // Mimic the const-length path: writer pulls `(const_qual_len * bits + 7)/8`
        // bytes; binning is Illumina8 (3 bits). l_seq = const_qual_len.
        let qual_ascii: &[u8] = b"IIIIIIII"; // 8 bases at Phred 40
        let binning = QualityBinning::Illumina8;
        let packed = columnar::pack_qualities(qual_ascii, binning).unwrap();

        let mut out = Vec::new();
        assemble_fastq_record(
            b"@r",
            b"AAAAAAAA",
            0,
            QualityInput::Packed {
                l_seq: 8,
                packed: &packed,
                binning,
            },
            false,
            &mut out,
        )
        .unwrap();
        // Illumina8 decodes Phred 40 -> bin 7 -> Phred 40 -> ASCII 'I'.
        assert_eq!(&out, b"@r\nAAAAAAAA\n+\nIIIIIIII\n");
    }

    #[test]
    fn test_format_batch_parallel_preserves_order() {
        // 5 records, no quality, no RC. Verifies par_chunks order is preserved.
        let headers: Vec<Vec<u8>> = (0..5).map(|i| format!("@r{i}").into_bytes()).collect();
        let seqs: Vec<Vec<u8>> = (0..5).map(|i| format!("SEQ{i}").into_bytes()).collect();
        let rc: Vec<u8> = Vec::new();
        let quals: Vec<QualityInput<'_>> = vec![QualityInput::None; 5];

        let out = format_batch_parallel(&headers, &seqs, &rc, &quals, false, 32).unwrap();
        let combined: Vec<u8> = out.into_iter().flatten().collect();
        let expected =
            b"@r0\nSEQ0\n+\n\n@r1\nSEQ1\n+\n\n@r2\nSEQ2\n+\n\n@r3\nSEQ3\n+\n\n@r4\nSEQ4\n+\n\n";
        assert_eq!(&combined, expected);
    }

    #[test]
    fn test_cross_stream_index_validation_active_streams_only() {
        let headers_idx = StreamIndex {
            blocks: vec![BlockEntry {
                record_count: 50,
                byte_offset: 0,
                byte_len: 0,
                expected_crc: 0,
            }],
            total_records: 50,
        };
        let sequences_idx = StreamIndex {
            blocks: vec![BlockEntry {
                record_count: 50,
                byte_offset: 0,
                byte_len: 0,
                expected_crc: 0,
            }],
            total_records: 50,
        };
        validate_indices(
            50,
            &[
                ("headers", Some(&headers_idx)),
                ("sequences", Some(&sequences_idx)),
                ("qualities", None), // inactive — None passes through
            ],
        )
        .expect("active streams agree, inactive=None must validate");
    }

    #[test]
    fn test_cross_stream_index_validation_total_mismatch_bails() {
        let h = StreamIndex {
            blocks: vec![],
            total_records: 50,
        };
        let s = StreamIndex {
            blocks: vec![],
            total_records: 49,
        };
        let err =
            validate_indices(50, &[("headers", Some(&h)), ("sequences", Some(&s))]).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("stream length mismatch") || msg.contains("≠ archive num_reads"),
            "expected mismatch error, got: {msg}"
        );
    }

    #[test]
    fn test_stream_cursor_reads_varint_records_across_block_boundary() {
        use super::super::write_varint;
        use std::sync::mpsc;
        // Build two blocks: block 0 = [record 0 (5 bytes), record 1 (10 bytes)]
        //                   block 1 = [record 2 (8 bytes)]
        let mut b0 = Vec::new();
        write_varint(&mut b0, 5).unwrap();
        b0.extend_from_slice(b"AAAAA");
        write_varint(&mut b0, 10).unwrap();
        b0.extend_from_slice(b"BBBBBBBBBB");

        let mut b1 = Vec::new();
        write_varint(&mut b1, 8).unwrap();
        b1.extend_from_slice(b"CCCCCCCC");

        let (tx, rx) = mpsc::sync_channel::<Result<(u32, Vec<u8>), String>>(2);
        tx.send(Ok((2, b0))).unwrap();
        tx.send(Ok((1, b1))).unwrap();
        drop(tx);

        let mut cursor = StreamCursor::new(rx, "test");
        let r0 = cursor.read_one_record_varint().expect("r0");
        assert_eq!(r0, b"AAAAA");
        let r1 = cursor.read_one_record_varint().expect("r1");
        assert_eq!(r1, b"BBBBBBBBBB");
        let r2 = cursor.read_one_record_varint().expect("r2");
        assert_eq!(r2, b"CCCCCCCC");

        // Drained cursor — next call should error cleanly.
        let err = cursor.read_one_record_varint();
        assert!(err.is_err());
    }

    #[test]
    fn test_stream_cursor_read_one_byte() {
        use std::sync::mpsc;
        // Block of 5 bytes, 5 records (1 byte per record — RC flags pattern).
        let block = vec![0u8, 1, 0, 1, 1];
        let (tx, rx) = mpsc::sync_channel::<Result<(u32, Vec<u8>), String>>(1);
        tx.send(Ok((5, block))).unwrap();
        drop(tx);

        let mut cursor = StreamCursor::new(rx, "rc_test");
        assert_eq!(cursor.read_one_byte().unwrap(), 0);
        assert_eq!(cursor.read_one_byte().unwrap(), 1);
        assert_eq!(cursor.read_one_byte().unwrap(), 0);
        assert_eq!(cursor.read_one_byte().unwrap(), 1);
        assert_eq!(cursor.read_one_byte().unwrap(), 1);
    }

    #[test]
    fn test_stream_cursor_read_one_record_const() {
        use std::sync::mpsc;
        // Const-length sequences: each record is 4 bytes, no varint length prefix.
        let block: Vec<u8> = b"ACGT".repeat(3); // 3 records, 4 bytes each
        let (tx, rx) = mpsc::sync_channel::<Result<(u32, Vec<u8>), String>>(1);
        tx.send(Ok((3, block))).unwrap();
        drop(tx);

        let mut cursor = StreamCursor::new(rx, "const_test");
        assert_eq!(cursor.read_one_record_const(4).unwrap(), b"ACGT");
        assert_eq!(cursor.read_one_record_const(4).unwrap(), b"ACGT");
        assert_eq!(cursor.read_one_record_const(4).unwrap(), b"ACGT");
    }

    #[test]
    fn test_per_stream_decompressor_emits_blocks_via_bounded_channel() {
        use std::sync::mpsc::sync_channel;
        use std::thread;

        // Build a synthetic per-block-framed BSC stream with 2 blocks, 5 records each.
        let block_a = crate::compression::bsc::compress_adaptive_no_lzp(b"ABCDEFGHIJ").unwrap();
        let block_b = crate::compression::bsc::compress_adaptive_no_lzp(b"KLMNOPQRST").unwrap();
        let mut stream = Vec::new();
        stream.extend_from_slice(&2u32.to_le_bytes()); // num_blocks
        // Block A payload starts at 4 (num_blocks) + 12 (frame header).
        let off_a = 4 + 12;
        crate::compression::bsc::write_block_frame_header(
            &mut stream,
            block_a.len() as u32,
            5,
            &block_a,
        );
        stream.extend_from_slice(&block_a);
        let off_b = off_a + block_a.len() + 12;
        crate::compression::bsc::write_block_frame_header(
            &mut stream,
            block_b.len() as u32,
            5,
            &block_b,
        );
        stream.extend_from_slice(&block_b);

        // Write to a temp file so the decompressor thread can seek+read.
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), &stream).unwrap();
        let path = tmp.path().to_path_buf();

        // Construct the per-stream index directly (the live decode path builds it
        // from the v5 footer via `build_role_index`; here we know the offsets).
        let index = StreamIndex {
            blocks: vec![
                BlockEntry {
                    record_count: 5,
                    byte_offset: off_a as u64,
                    byte_len: block_a.len() as u32,
                    expected_crc: crate::compression::bsc::compute_block_frame_crc(5, &block_a),
                },
                BlockEntry {
                    record_count: 5,
                    byte_offset: off_b as u64,
                    byte_len: block_b.len() as u32,
                    expected_crc: crate::compression::bsc::compute_block_frame_crc(5, &block_b),
                },
            ],
            total_records: 10,
        };
        assert_eq!(index.blocks.len(), 2);

        let (tx, rx) = sync_channel::<Result<(u32, Vec<u8>), String>>(3);
        let path_clone = path.clone();
        let index_clone = index.clone();
        let handle = thread::spawn(move || {
            bounded_bsc_producer(
                &path_clone,
                &index_clone,
                tx,
                MAX_PARALLEL_DECOMPRESS,
                BSC_MAX_BLOCK_LEN,
            )
            .unwrap();
        });

        let (rc0, payload0) = rx.recv().unwrap().unwrap();
        let (rc1, payload1) = rx.recv().unwrap().unwrap();
        assert_eq!(rc0, 5);
        assert_eq!(payload0, b"ABCDEFGHIJ");
        assert_eq!(rc1, 5);
        assert_eq!(payload1, b"KLMNOPQRST");
        handle.join().unwrap();
    }

    #[test]
    fn test_bounded_columnar_header_producer_emits_varint_records() {
        use std::sync::mpsc::sync_channel;
        use std::thread;

        // Build a single-block columnar header stream:
        //   [num_blocks=1][frame header][chunk_reads:u32][columnar_blob]
        let headers: Vec<String> = (0..3).map(|i| format!("@r{i}/1")).collect();
        let header_refs: Vec<&str> = headers.iter().map(|s| s.as_str()).collect();
        let blob = crate::compression::header_col::compress_headers_columnar(&header_refs).unwrap();
        let mut payload = Vec::new();
        payload.extend_from_slice(&(headers.len() as u32).to_le_bytes());
        payload.extend_from_slice(&blob);

        let mut stream = Vec::new();
        stream.extend_from_slice(&1u32.to_le_bytes()); // num_blocks
        crate::compression::bsc::write_block_frame_header(
            &mut stream,
            payload.len() as u32,
            headers.len() as u32,
            &payload,
        );
        stream.extend_from_slice(&payload);

        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), &stream).unwrap();
        let path = tmp.path().to_path_buf();

        // Construct the per-stream index directly (the live decode path builds it
        // from the v5 footer via `build_role_index`; the single block payload
        // starts at 4 (num_blocks) + 12 (frame header)).
        let index = StreamIndex {
            blocks: vec![BlockEntry {
                record_count: headers.len() as u32,
                byte_offset: (4 + 12) as u64,
                byte_len: payload.len() as u32,
                expected_crc: crate::compression::bsc::compute_block_frame_crc(
                    headers.len() as u32,
                    &payload,
                ),
            }],
            total_records: headers.len() as u64,
        };

        let (tx, rx) = sync_channel::<Result<(u32, Vec<u8>), String>>(2);
        let path_clone = path.clone();
        let index_clone = index.clone();
        let handle = thread::spawn(move || {
            bounded_columnar_header_producer(
                &path_clone,
                &index_clone,
                tx,
                MAX_PARALLEL_DECOMPRESS,
            )
            .unwrap();
        });

        let (rc, varint_stream) = rx.recv().unwrap().unwrap();
        assert_eq!(rc, 3);

        // Decode the varint-prefixed stream and verify the round-trip.
        let mut pos = 0usize;
        let mut decoded: Vec<String> = Vec::new();
        while pos < varint_stream.len() {
            let len = super::read_varint(&varint_stream, &mut pos).unwrap();
            decoded.push(String::from_utf8(varint_stream[pos..pos + len].to_vec()).unwrap());
            pos += len;
        }
        assert_eq!(decoded, headers);
        handle.join().unwrap();
    }

    #[test]
    fn footer_indices_match_directory() {
        use crate::compression::chunk_directory::{ChunkDirEntry, ChunkDirectory, StreamRole};
        // Two header frames and one sequence frame; verify per-role StreamIndex.
        let dir = ChunkDirectory {
            num_reads: 7,
            num_chunks: 2,
            entries: vec![
                ChunkDirEntry {
                    chunk_index: 0,
                    role: StreamRole::Headers,
                    mate: 0,
                    codec: 1,
                    offset: 100,
                    length: 12 + 30,
                    record_count: 4,
                },
                ChunkDirEntry {
                    chunk_index: 0,
                    role: StreamRole::Sequence,
                    mate: 0,
                    codec: 1,
                    offset: 142,
                    length: 12 + 50,
                    record_count: 4,
                },
                ChunkDirEntry {
                    chunk_index: 1,
                    role: StreamRole::Headers,
                    mate: 0,
                    codec: 1,
                    offset: 204,
                    length: 12 + 20,
                    record_count: 3,
                },
            ],
        };
        let idx = build_role_index(&dir, StreamRole::Headers).unwrap();
        assert_eq!(idx.total_records, 7);
        assert_eq!(idx.blocks.len(), 2);
        assert_eq!(idx.blocks[0].byte_offset, 112); // 100 + 12
        assert_eq!(idx.blocks[0].byte_len, 30);
        assert_eq!(idx.blocks[1].byte_offset, 216); // 204 + 12
        assert_eq!(idx.blocks[1].record_count, 3);
    }

    #[test]
    fn build_role_index_range_filters_chunks() {
        use crate::compression::chunk_directory::{ChunkDirEntry, ChunkDirectory, StreamRole};
        let mk = |chunk: u32, rc: u64, off: u64| ChunkDirEntry {
            chunk_index: chunk,
            role: StreamRole::Headers,
            mate: 0,
            codec: 1,
            offset: off,
            length: 100,
            record_count: rc,
        };
        let dir = ChunkDirectory {
            num_reads: 30,
            num_chunks: 3,
            entries: vec![mk(0, 10, 1000), mk(1, 10, 2000), mk(2, 10, 3000)],
        };
        // whole archive
        let all = build_role_index_range(&dir, StreamRole::Headers, None).unwrap();
        assert_eq!(all.total_records, 30);
        assert_eq!(all.blocks.len(), 3);
        // chunks [1,3)
        let tail = build_role_index_range(&dir, StreamRole::Headers, Some((1, 3))).unwrap();
        assert_eq!(tail.total_records, 20);
        assert_eq!(tail.blocks.len(), 2);
        assert_eq!(tail.blocks[0].byte_offset, 2000 + 12);
        // existing builder == whole-archive range
        let legacy = build_role_index(&dir, StreamRole::Headers).unwrap();
        assert_eq!(legacy.total_records, all.total_records);
        assert_eq!(legacy.blocks.len(), all.blocks.len());
    }

    #[test]
    fn mate_range_predicate_keeps_chunk_order() {
        // The retain+sort invariant: filtering to [a,b) then sorting by chunk_index
        // yields a contiguous, ascending subset.
        let mut chunks: Vec<u32> = vec![3, 0, 2, 1, 4];
        let (a, b) = (1u32, 4u32);
        chunks.retain(|&c| c >= a && c < b);
        chunks.sort_unstable();
        assert_eq!(chunks, vec![1, 2, 3]);
    }

    #[test]
    fn fastq_output_len_matches_formatter() {
        let mut out = Vec::new();
        assemble_fastq_record(b"@r1 x", b"ACGT", 0, QualityInput::Decoded(b"IIII"), false, &mut out).unwrap();
        assert_eq!(out.len() as u64, fastq_output_len(false, b"@r1 x".len(), 4, 4));
        out.clear();
        assemble_fastq_record(b"@r1", b"ACGT", 1, QualityInput::Decoded(b"IIII"), false, &mut out).unwrap();
        assert_eq!(out.len() as u64, fastq_output_len(false, 3, 4, 4)); // RC keeps length
        out.clear();
        assemble_fastq_record(b"@r2", b"ACGT", 0, QualityInput::None, false, &mut out).unwrap();
        assert_eq!(out.len() as u64, fastq_output_len(false, 3, 4, 0)); // no-quality FASTQ
        out.clear();
        assemble_fastq_record(b">r3", b"ACGTAC", 0, QualityInput::None, true, &mut out).unwrap();
        assert_eq!(out.len() as u64, fastq_output_len(true, 3, 6, 0)); // FASTA
    }

    #[test]
    fn single_end_validator_global_contract() {
        use crate::compression::chunk_directory::{ChunkDirEntry, ChunkDirectory, StreamRole, GLOBAL_SENTINEL};
        let good = ChunkDirEntry { chunk_index: GLOBAL_SENTINEL, role: StreamRole::ChunkDecodedSizes, mate: 0, codec: 0, offset: 120, length: 25, record_count: 1 };
        let seq = ChunkDirEntry { chunk_index: 0, role: StreamRole::Sequence, mate: 0, codec: 1, offset: 100, length: 20, record_count: 1 };
        let ok = ChunkDirectory { num_reads: 1, num_chunks: 1, entries: vec![seq.clone(), good.clone()] };
        assert!(validate_single_end_global_entries(&ok).is_ok());
        // no globals at all is fine
        let none = ChunkDirectory { num_reads: 1, num_chunks: 1, entries: vec![seq.clone()] };
        assert!(validate_single_end_global_entries(&none).is_ok());
        // foreign global role rejected
        let mut foreign = ok.clone(); foreign.entries[1].role = StreamRole::Sequence;
        assert!(validate_single_end_global_entries(&foreign).is_err());
        // duplicate ChunkDecodedSizes global rejected
        let mut dup = ok.clone(); dup.entries.push(good.clone());
        assert!(validate_single_end_global_entries(&dup).is_err());
        // wrong mate / codec / record_count rejected
        let mut wm = ok.clone(); wm.entries[1].mate = 1; assert!(validate_single_end_global_entries(&wm).is_err());
        let mut wc = ok.clone(); wc.entries[1].codec = 1; assert!(validate_single_end_global_entries(&wc).is_err());
        let mut wr = ok.clone(); wr.entries[1].record_count = 2; assert!(validate_single_end_global_entries(&wr).is_err());
    }
}
