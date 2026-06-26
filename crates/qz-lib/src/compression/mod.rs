pub mod bsc;
mod columnar;
pub mod dna_utils;
pub mod fqzcomp;
mod header_token;
mod quality;
pub(crate) mod chunk_directory;
pub(crate) mod codec_ids;
mod codecs;
mod compress_impl;
mod decompress_impl;
mod merge;
pub(crate) mod stream_compress;
pub mod header_col;
pub(crate) mod paired;
pub mod reference;
mod ultra;
mod range_compress;
pub(crate) mod cluster;

use crate::cli::{
    CompressConfig, DecompressConfig, HeaderCompressor, QualityCompressor, QualityMode,
    SequenceCompressor,
};
use anyhow::Result;
use tracing::info;

/// Minimum reads for fqzcomp's context model to outperform BSC (below this, the
/// model doesn't converge). Named for the legacy `--quality-ctx-block-size` flag.
const MIN_READS_QUALITY_CTX: usize = 100_000;

/// QZ archive magic bytes (file identification)
const ARCHIVE_MAGIC: [u8; 2] = *b"QZ";
/// Chunk-major physical layout — the only on-disk format qz writes and reads.
/// v5 archives carry no num_reads / stream-length fields in the front header
/// (those live in the directory footer). `encoding_type` carries sequence
/// semantics. Legacy v4 (stream-major) archives are no longer decodable.
pub(crate) const ARCHIVE_VERSION_CHUNK_MAJOR: u8 = 5;
/// Size of the v5 fixed body: 12 bytes (6 live codec/flag bytes at offsets 0–5,
/// then 6 reserved-zero legacy-metadata bytes at offsets 6–11).
pub(crate) const V5_FIXED_BODY_SIZE: usize = 12;

/// Size of the v2+ prefix: magic(2) + version(1) + archive_type(1) + header_size(4).
/// Byte 3 is `archive_type` (0=single-end, 1=paired, 2=reference, 3=cluster-reorder, 4=single-end-reference).
pub(crate) const V2_PREFIX_SIZE: usize = 8;
// ── Archive flags byte (offset `hdr_off::FLAGS`) ─────────────────────────────
/// bit1: constant per-record sequence/quality lengths (8 trailing body bytes).
pub(crate) const FLAG_CONST_LENGTHS: u8 = 0x02;
/// bit2: input was FASTA (no quality, '>' headers), so the decoder emits FASTA
/// (`>id\nseq\n`) rather than FASTQ.
pub(crate) const FLAG_FASTA: u8 = 0x04;

/// Byte offsets of every field in the fixed v5 archive body, measured from the
/// start of the body (i.e. after the 8-byte [`V2_PREFIX_SIZE`] prefix). This is
/// the single source of truth for the layout that [`FixedHeader::write_v5`],
/// [`FixedHeader::parse_v5`], and the streaming decoders in `decompress_impl`
/// share — previously each reader hard-coded `b + 12` / `b + 36` style magic
/// numbers that had to be kept in sync by hand.
pub(crate) mod hdr_off {
    pub const ENCODING_TYPE: usize = 0;
    pub const FLAGS: usize = 1;
    pub const QUALITY_BINNING: usize = 2;
    pub const QUALITY_COMPRESSOR: usize = 3;
    pub const SEQUENCE_COMPRESSOR: usize = 4;
    pub const HEADER_COMPRESSOR: usize = 5;
    // Body bytes 6–11 are reserved padding (legacy quality-model/delta/dict,
    // template-prefix-len, has-comment — all removed). Always zero; never read.
}
/// I/O buffer size for reading archive files during decompression
const IO_BUFFER_SIZE: usize = 8 * 1024 * 1024; // 8 MB

/// Read a little-endian u16 from `data` at `offset`, returning an error on truncation.
fn read_le_u16(data: &[u8], offset: usize) -> Result<u16> {
    offset
        .checked_add(2)
        .and_then(|end| data.get(offset..end))
        .and_then(|s| <[u8; 2]>::try_from(s).ok())
        .map(u16::from_le_bytes)
        .ok_or_else(|| anyhow::anyhow!("truncated archive at offset {offset}"))
}

/// Read a little-endian u32 from `data` at `offset`, returning an error on truncation.
fn read_le_u32(data: &[u8], offset: usize) -> Result<u32> {
    offset
        .checked_add(4)
        .and_then(|end| data.get(offset..end))
        .and_then(|s| <[u8; 4]>::try_from(s).ok())
        .map(u32::from_le_bytes)
        .ok_or_else(|| anyhow::anyhow!("truncated archive at offset {offset}"))
}

/// Convert QualityMode to QualityBinning
fn quality_mode_to_binning(mode: QualityMode) -> QualityBinning {
    match mode {
        QualityMode::Lossless => QualityBinning::None, // None = lossless (keep original)
        QualityMode::IlluminaBin => QualityBinning::Illumina8,
        QualityMode::Illumina4 => QualityBinning::FourLevel,
        QualityMode::Binary => QualityBinning::Binary { threshold: 20 },
        QualityMode::Discard => QualityBinning::Binary { threshold: 255 }, // All qualities -> same value
    }
}

/// Return the human-readable display name for a QZ encoding-type byte code.
///
/// This is the single source of truth used by both the CLI's `verify` output
/// and any other caller that needs to display an encoding type. Keeping it here
/// prevents the CLI from drifting out of sync with new encoding types.
pub fn encoding_type_display_name(code: u8) -> &'static str {
    match code {
        0 => "raw",
        4 => "raw+hints",
        6 => "rc-canon",
        10 => "ultra big-block",
        11 => "cluster-reorder",
        _ => "unknown",
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
pub(crate) fn code_to_binning(code: u8) -> Result<QualityBinning> {
    match code {
        0 => Ok(QualityBinning::None), // Lossless
        1 => Ok(QualityBinning::Illumina8),
        2 => Ok(QualityBinning::Binary { threshold: 20 }),
        3 => Ok(QualityBinning::FourLevel),
        _ => anyhow::bail!("Invalid quality binning code: {}", code),
    }
}

/// Resolve `Auto` to a concrete quality compressor based on input characteristics.
/// Non-Auto values are returned unchanged.
pub(crate) fn resolve_quality_compressor(
    requested: QualityCompressor,
    record_count: usize,
    quality_mode: QualityMode,
    no_quality: bool,
    modeling_or_delta: bool,
) -> QualityCompressor {
    if requested != QualityCompressor::Auto {
        return requested;
    }
    // Context-modeled fqzcomp only helps for lossless quality with enough reads
    // to amortize its model.
    if no_quality || quality_mode != QualityMode::Lossless || record_count < MIN_READS_QUALITY_CTX {
        return QualityCompressor::Bsc;
    }
    // Quality modeling / delta encode model-deltas through a byte backend; fqzcomp
    // has no model-delta path. Resolving to it would make the encoder substitute
    // BSC for the deltas while the header still advertised fqzcomp — an unreadable
    // archive. Resolve to BSC so the stored code always matches the bytes written.
    // (See the explicit-Fqzcomp + modeling reject in compress_impl, which covers
    // the non-Auto case.)
    if modeling_or_delta {
        return QualityCompressor::Bsc;
    }
    // Single-end, paired, and reference: fqzcomp is the only context-modeled
    // quality backend and is bounded-streaming (record-capped multi-block).
    QualityCompressor::Fqzcomp
}

#[cfg(test)]
mod streaming_equiv_tests;

#[cfg(test)]
mod resolve_quality_compressor_tests {
    use super::*;

    const SCALE: usize = MIN_READS_QUALITY_CTX + 1;

    #[test]
    fn explicit_choice_passes_through() {
        // A non-Auto request is honored verbatim regardless of the other inputs.
        for req in [QualityCompressor::Bsc, QualityCompressor::Fqzcomp] {
            assert_eq!(
                resolve_quality_compressor(req, SCALE, QualityMode::Lossless, false, false),
                req
            );
        }
    }

    #[test]
    fn lossless_at_scale_picks_fqzcomp() {
        // Auto + lossless + enough reads resolves to fqzcomp for every mode now
        // (single-end, paired, reference all share the one context-modeled backend).
        assert_eq!(
            resolve_quality_compressor(
                QualityCompressor::Auto,
                SCALE,
                QualityMode::Lossless,
                false,
                false,
            ),
            QualityCompressor::Fqzcomp
        );
    }

    #[test]
    fn modeling_or_delta_resolves_to_bsc() {
        // F1: model deltas need a byte backend; never resolve Auto to fqzcomp when
        // modeling/delta is on (would store a code that mismatches the bytes).
        assert_eq!(
            resolve_quality_compressor(
                QualityCompressor::Auto,
                SCALE,
                QualityMode::Lossless,
                false,
                true, // modeling_or_delta
            ),
            QualityCompressor::Bsc
        );
    }

    #[test]
    fn below_threshold_or_no_quality_is_bsc() {
        // Too few reads → BSC.
        assert_eq!(
            resolve_quality_compressor(
                QualityCompressor::Auto,
                MIN_READS_QUALITY_CTX - 1,
                QualityMode::Lossless,
                false,
                false,
            ),
            QualityCompressor::Bsc
        );
        // no_quality → BSC.
        assert_eq!(
            resolve_quality_compressor(
                QualityCompressor::Auto,
                SCALE,
                QualityMode::Lossless,
                true,
                false,
            ),
            QualityCompressor::Bsc
        );
    }
}

/// The record/sequence framing recorded in the archive's `encoding_type` byte.
///
/// Centralizes the magic numbers the encoder writes (`compress_impl` /
/// `ultra`) and the decoder dispatches on, so "type 4 means sequence
/// hints", "type 6 appends an RC stream", and "which framings the streaming
/// path handles" each live in exactly one place instead of being re-derived as
/// `encoding_type == 4` / `== 6` in every reader. Legacy values 1/2/3/5/7/8/9
/// are removed encodings and map to `None`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum EncodingType {
    Raw = 0,
    RawWithHints = 4,
    RcCanon = 6,
    /// Bounded-streaming big-block ultra encoding: raw BSC streams with
    /// per-block record-count metadata, decoded via the streaming path.
    UltraBigBlock = 10,
    /// Order-drop cluster encoding: sequences reordered by minimizer-cluster,
    /// compressed with long-range LZ (zstd), decoded via the streaming path.
    ClusterReorder = 11,
}

impl EncodingType {
    pub fn to_u8(self) -> u8 {
        self as u8
    }

    /// Decode the byte. Removed encodings (1/2/3/5/7/8/9) and unknown values → `None`.
    pub fn from_u8(code: u8) -> Option<EncodingType> {
        Some(match code {
            0 => EncodingType::Raw,
            4 => EncodingType::RawWithHints,
            6 => EncodingType::RcCanon,
            10 => EncodingType::UltraBigBlock,
            11 => EncodingType::ClusterReorder,
            _ => return None,
        })
    }

    /// Sequence hints (type 4) reorder reads for BWT clustering on decode.
    pub fn has_sequence_hints(self) -> bool {
        matches!(self, EncodingType::RawWithHints)
    }

    /// RC-canonical (type 6) appends a reverse-complement flags stream.
    pub fn has_rc_canon(self) -> bool {
        matches!(self, EncodingType::RcCanon)
    }

    /// The bounded-streaming decoder only handles the order-preserving raw
    /// framings; every remaining encoding qualifies.
    pub fn is_bounded_streamable(self) -> bool {
        matches!(
            self,
            EncodingType::Raw
                | EncodingType::RawWithHints
                | EncodingType::RcCanon
                | EncodingType::UltraBigBlock
                | EncodingType::ClusterReorder
        )
    }
}

pub use columnar::QualityBinning;
pub use reference::reference_index_sidecar_path;
pub use reference::{build_reference_index, ReferenceIndexBuildStats};
pub use reference::{ensure_reference_index, reference_index_status, ReferenceIndexStatus};

/// Read just the first record's sequence length from a FASTQ file. Used to
/// resolve the seeding profile for reference-index preflight and `qz index --like`.
pub fn peek_first_read_length(path: &std::path::Path) -> anyhow::Result<usize> {
    let mut reader = crate::io::fastq::FastqReader::from_path(path, false)?;
    let rec = reader
        .next()?
        .ok_or_else(|| anyhow::anyhow!("{} has no reads", path.display()))?;
    Ok(rec.sequence.len())
}

/// Read the first record's sequence length within a record-aligned byte range of a
/// PLAIN (non-gzip) FASTQ. Used by the NUMA driver to resolve each shard's profile.
/// Returns `Ok(None)` if the range contains no record.
pub fn peek_read_length_in_range(
    path: &std::path::Path,
    start: u64,
    end: u64,
) -> anyhow::Result<Option<usize>> {
    let mut reader = crate::io::fastq::FastqReader::from_range(path, start, end, false)?;
    Ok(reader.next()?.map(|rec| rec.sequence.len()))
}

// Advanced and tier2 features not currently exposed in public API

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
                chunk_orig += record.id.len()
                    + record.sequence.len()
                    + record.quality.as_ref().map(|q| q.len()).unwrap_or(0)
                    + 3;
                records.push(record);
            }
            None => break,
        }
    }

    Ok((records, chunk_bases, chunk_orig))
}

/// Append one record's framed bytes to the per-stream buffers. Offsets are
/// managed by callers (push `buf.len()` before calling for the record's start
/// offset). Used by the streaming stagers so framing can never drift. (Cap-check
/// errors carry no record index — it isn't available at this layer; callers may
/// prefix one if needed.)
#[allow(clippy::too_many_arguments)]
pub(crate) fn frame_record(
    record: &crate::io::FastqRecord,
    quality_mode: QualityMode,
    quality_binning: QualityBinning,
    no_quality: bool,
    sequence_hints: bool,
    rc_canon: bool,
    const_seq_len: usize,
    const_qual_len: usize,
    build_header_stream: bool,
    header_stream: &mut Vec<u8>,
    seq_stream: &mut Vec<u8>,
    qual_stream: &mut Vec<u8>,
    rc_flags: &mut Vec<u8>,
) -> Result<()> {
    const PER_RECORD_MAX_BYTES: usize = codecs::PER_RECORD_MAX_BYTES;
    if record.id.len() > PER_RECORD_MAX_BYTES {
        anyhow::bail!(
            "FASTQ record header is {} bytes, exceeds {}-byte per-record cap for bounded streaming",
            record.id.len(),
            PER_RECORD_MAX_BYTES
        );
    }
    if record.sequence.len() > PER_RECORD_MAX_BYTES {
        anyhow::bail!(
            "FASTQ record sequence is {} bytes, exceeds {}-byte per-record cap for bounded streaming",
            record.sequence.len(),
            PER_RECORD_MAX_BYTES
        );
    }
    if let Some(ref qual) = record.quality
        && qual.len() > PER_RECORD_MAX_BYTES
    {
        anyhow::bail!(
            "FASTQ record quality is {} bytes, exceeds {}-byte per-record cap for bounded streaming",
            qual.len(),
            PER_RECORD_MAX_BYTES
        );
    }

    if build_header_stream {
        write_varint(header_stream, record.id.len())?;
        header_stream.extend_from_slice(&record.id);
    }

    if const_seq_len == 0 {
        write_varint(seq_stream, record.sequence.len())?;
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

    if !no_quality && let Some(ref qual) = record.quality {
        let quantized = quality::quantize_quality(qual, quality_mode);
        if const_qual_len == 0 {
            write_varint(qual_stream, quantized.len())?;
        }
        let packed = columnar::pack_qualities(&quantized, quality_binning)?;
        qual_stream.extend_from_slice(&packed);
    }
    Ok(())
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
            if records.iter().all(|r| {
                r.quality
                    .as_ref()
                    .is_some_and(|q| q.len() == first_qual_len)
            }) {
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

/// Per-block prefix size in the qz archive v4 wire format:
/// 4-byte length + 4-byte record_count + 4-byte CRC32 = 12 bytes.
pub(crate) const BLOCK_PREFIX_SIZE: usize = 12;

/// The v5 fixed front header (v2 prefix + 12-byte body + optional 8-byte
/// const-length tail) — the header shape every chunk-major compress path emits
/// and every streaming decoder reads. Centralizing it here means the on-disk
/// layout lives in exactly one place: a new field is added once, to
/// [`FixedHeader::write_v5`]/[`FixedHeader::parse_v5`] and the [`hdr_off`]
/// constants, instead of being threaded by hand through several hand-rolled
/// writers/readers (the drift that silently dropped the FASTA flag on the ultra
/// path).
///
/// Codes are stored raw (`u8`) rather than as enums so `parse_v5` is infallible
/// on codec codes and matches the streaming readers, which compare raw codes and
/// validate (or not) per call site. num_reads and per-stream byte lengths are
/// NOT carried here — they live in the v5 directory footer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct FixedHeader {
    pub encoding_type: u8,
    pub flags: u8,
    /// archive_type byte (byte 3 of the on-disk prefix).
    /// 0 = single-end, 1 = paired-end, 2 = paired reference, 3 = cluster,
    /// 4 = single-end reference. Values > 4 are rejected by `parse_v5`.
    pub archive_type: u8,
    pub quality_binning_code: u8,
    pub quality_compressor_code: u8,
    pub sequence_compressor_code: u8,
    pub header_compressor_code: u8,
    // Body bytes 6–11 (6 bytes: a u8 trio + a u16 + a u8) are reserved padding.
    // They held the legacy quality-modeling/delta/dict/template-prefix/has-comment
    // metadata (all removed with the v4 format). The encoder always writes them as
    // 0 and the decoder skips them; they are kept in the 12-byte body for forward
    // layout stability, so they are NOT modeled as struct fields.
    pub const_seq_len: u32,
    pub const_qual_len: u32,
}

impl FixedHeader {
    pub fn is_fasta(&self) -> bool {
        self.flags & FLAG_FASTA != 0
    }
    pub fn has_const_lengths(&self) -> bool {
        self.flags & FLAG_CONST_LENGTHS != 0
    }
    /// Whether the encoding type carries per-record sequence hints (derived from
    /// `encoding_type`). Unknown codes map to `false`.
    pub fn has_sequence_hints(&self) -> bool {
        EncodingType::from_u8(self.encoding_type).is_some_and(|e| e.has_sequence_hints())
    }
    /// Whether the encoding type carries reverse-complement canonicalization
    /// flags (derived from `encoding_type`). Unknown codes map to `false`.
    pub fn has_rc_canon(&self) -> bool {
        EncodingType::from_u8(self.encoding_type).is_some_and(|e| e.has_rc_canon())
    }

    /// On-disk size of the v5 front header (prefix + 12-byte body + optional tail).
    pub fn v5_byte_len(&self) -> usize {
        V2_PREFIX_SIZE + V5_FIXED_BODY_SIZE + if self.has_const_lengths() { 8 } else { 0 }
    }

    /// Serialize the v5 front header: v2 prefix (version 5) + the 12-byte body
    /// (6 codec/flag bytes + 6 reserved-zero bytes) + optional const-length tail.
    /// Omits num_reads and the four stream lengths (they live in the footer).
    pub fn write_v5(&self, out: &mut impl std::io::Write) -> Result<usize> {
        let header_size = self.v5_byte_len();
        out.write_all(&ARCHIVE_MAGIC)?;
        out.write_all(&[ARCHIVE_VERSION_CHUNK_MAJOR, self.archive_type])?;
        out.write_all(&(header_size as u32).to_le_bytes())?;

        out.write_all(&[self.encoding_type])?;
        out.write_all(&[self.flags])?;
        out.write_all(&[self.quality_binning_code])?;
        out.write_all(&[self.quality_compressor_code])?;
        out.write_all(&[self.sequence_compressor_code])?;
        out.write_all(&[self.header_compressor_code])?;
        // Reserved body bytes 6–11 (legacy quality-model/delta/dict/template/comment
        // metadata, removed): always zero. Kept for the fixed 12-byte body layout.
        out.write_all(&[0u8; 6])?;

        if self.has_const_lengths() {
            out.write_all(&self.const_seq_len.to_le_bytes())?;
            out.write_all(&self.const_qual_len.to_le_bytes())?;
        }
        Ok(header_size)
    }

    /// Parse the v5 front header (STRICT). Validates magic, version, the stored
    /// `header_size`, rejects unknown flag bits, and requires the const-length
    /// tail to be present exactly when `FLAG_CONST_LENGTHS` is set. num_reads /
    /// stream lengths are set to 0 (supplied from the footer by the decoder).
    pub fn parse_v5(data: &[u8]) -> Result<FixedHeader> {
        if data.len() < V2_PREFIX_SIZE || data[0..2] != ARCHIVE_MAGIC {
            anyhow::bail!("Not a QZ archive (missing magic bytes)");
        }
        if data[2] != ARCHIVE_VERSION_CHUNK_MAJOR {
            anyhow::bail!(
                "parse_v5: expected chunk-major version {ARCHIVE_VERSION_CHUNK_MAJOR}, got {}",
                data[2]
            );
        }
        let archive_type = data[3];
        if archive_type > 4 {
            anyhow::bail!(
                "parse_v5: unknown archive_type {archive_type} (valid: 0=single-end, 1=paired, 2=reference, 3=cluster-reorder, 4=single-end-reference)"
            );
        }
        let stored_header_size = read_le_u32(data, 4)? as usize;
        let b = V2_PREFIX_SIZE;
        if data.len() < b + V5_FIXED_BODY_SIZE {
            anyhow::bail!("Invalid v5 archive: header truncated");
        }
        let flags = data[b + hdr_off::FLAGS];
        // Forward-safety: v5 single-end/ultra only ever sets const-lengths/fasta.
        // Reject FLAG_ARITHMETIC (v5 carries no read-length table) and any unknown bit.
        const KNOWN_V5_FLAGS: u8 = FLAG_CONST_LENGTHS | FLAG_FASTA;
        if flags & !KNOWN_V5_FLAGS != 0 {
            anyhow::bail!(
                "v5 header: unsupported flag bits 0x{:02x}",
                flags & !KNOWN_V5_FLAGS
            );
        }
        let has_const = flags & FLAG_CONST_LENGTHS != 0;
        let expected_size = b + V5_FIXED_BODY_SIZE + if has_const { 8 } else { 0 };
        if stored_header_size != expected_size {
            anyhow::bail!(
                "v5 header_size {stored_header_size} disagrees with flags 0x{flags:02x} (expected {expected_size})"
            );
        }
        if data.len() < expected_size {
            anyhow::bail!(
                "v5 header truncated: need {expected_size} bytes, have {}",
                data.len()
            );
        }
        let (const_seq_len, const_qual_len) = if has_const {
            (
                read_le_u32(data, b + V5_FIXED_BODY_SIZE)?,
                read_le_u32(data, b + V5_FIXED_BODY_SIZE + 4)?,
            )
        } else {
            (0, 0)
        };
        Ok(FixedHeader {
            encoding_type: data[b + hdr_off::ENCODING_TYPE],
            flags,
            archive_type,
            quality_binning_code: data[b + hdr_off::QUALITY_BINNING],
            quality_compressor_code: data[b + hdr_off::QUALITY_COMPRESSOR],
            sequence_compressor_code: data[b + hdr_off::SEQUENCE_COMPRESSOR],
            header_compressor_code: data[b + hdr_off::HEADER_COMPRESSOR],
            // Body bytes 6–11 (reserved legacy metadata) are skipped here; the
            // `data.len() < b + V5_FIXED_BODY_SIZE` check above guarantees they exist.
            const_seq_len,
            const_qual_len,
        })
    }
}

/// Read the archive version byte after validating the magic. Accepts any
/// version; callers dispatch (5 = chunk-major; legacy versions are rejected).
pub(crate) fn read_archive_version(data: &[u8]) -> Result<u8> {
    if data.len() < V2_PREFIX_SIZE || data[0..2] != ARCHIVE_MAGIC {
        anyhow::bail!("Not a QZ archive (missing magic bytes)");
    }
    Ok(data[2])
}

/// Write the v5 chunk-major front header (no num_reads / stream lengths — those
/// live in the directory footer). Returns the header byte size. Mirrors
/// `write_archive_header`'s argument shape minus the length/count fields.
pub(crate) fn write_v5_archive_header(
    out: &mut impl std::io::Write,
    encoding_type: u8,
    quality_binning: QualityBinning,
    quality_compressor: QualityCompressor,
    header_compressor: HeaderCompressor,
    fasta: bool,
    const_seq_len: usize,
    const_qual_len: usize,
    archive_type: u8,
) -> Result<usize> {
    // NB: no write-side archive_type guard on purpose — `parse_v5` is the single
    // validation point, and `v5_header_rejects_unknown_archive_type` writes an
    // out-of-range value here to exercise that parser rejection.
    let has_const_lengths = const_seq_len > 0 || const_qual_len > 0;
    let mut flags: u8 = if has_const_lengths {
        FLAG_CONST_LENGTHS
    } else {
        0
    };
    if fasta {
        flags |= FLAG_FASTA;
    }
    let hdr = FixedHeader {
        encoding_type,
        flags,
        archive_type,
        quality_binning_code: binning_to_code(quality_binning),
        quality_compressor_code: codec_ids::quality_compressor_to_codec_id(quality_compressor),
        sequence_compressor_code: codec_ids::sequence_compressor_to_codec_id(SequenceCompressor::Bsc),
        header_compressor_code: codec_ids::header_compressor_to_codec_id(header_compressor),
        const_seq_len: const_seq_len as u32,
        const_qual_len: const_qual_len as u32,
    };
    hdr.write_v5(out)
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
    let total_compressed =
        headers_len + sequences_len + qualities_len + rc_flags_len + metadata_size;
    info!("Compression completed in {:.2}s", elapsed.as_secs_f64());
    info!("Original size: {} bytes", original_size);
    info!("Compressed size: {} bytes", total_compressed);
    info!(
        "Compression ratio: {:.2}x",
        original_size as f64 / total_compressed as f64
    );
    info!("Stream breakdown:");
    info!(
        "  Headers:   {} bytes ({:.1}%)",
        headers_len,
        100.0 * headers_len as f64 / total_compressed as f64
    );
    info!(
        "  Sequences: {} bytes ({:.1}%)",
        sequences_len,
        100.0 * sequences_len as f64 / total_compressed as f64
    );
    info!(
        "  Qualities: {} bytes ({:.1}%)",
        qualities_len,
        100.0 * qualities_len as f64 / total_compressed as f64
    );
    if rc_flags_len > 0 {
        info!(
            "  RC flags:  {} bytes ({:.1}%)",
            rc_flags_len,
            100.0 * rc_flags_len as f64 / total_compressed as f64
        );
    }
    info!(
        "  Metadata:  {} bytes ({:.1}%)",
        metadata_size,
        100.0 * metadata_size as f64 / total_compressed as f64
    );
}

/// A process-unique suffix for atomic-output temp files. Combined with
/// `create_new` (O_EXCL), the sibling temp path is unpredictable — so it can't
/// collide with a concurrent job and can't be pre-planted as a symlink that a
/// later `File::create` would follow and truncate. Uniqueness: pid + a
/// per-process counter + a wall-clock nanosecond stamp (no extra dependency).
pub(crate) fn atomic_tmp_suffix() -> std::ffi::OsString {
    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::ffi::OsString::from(format!(".{}.{}.{}.tmp", std::process::id(), n, nanos))
}

/// Open a fresh file for writing with O_EXCL semantics (fails if the path already
/// exists, and never follows a symlink). Used to reserve atomic-output temp files.
pub(crate) fn create_new_for_write(path: &std::path::Path) -> std::io::Result<std::fs::File> {
    std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
}

/// Initialize rayon's process-global thread pool to honor `-t` / `num_threads`,
/// at every top-level entry. Previously only the single-end paths
/// (`compress_impl::compress`, `decompress_impl::decompress`) built the pool, so
/// a **paired or reference** invocation — which never routes through those —
/// silently fell back to rayon's default (all cores), ignoring `-t`. Calling
/// this from the top-level `compress`/`decompress` dispatch fixes that for every
/// mode.
///
/// rayon's global pool can be built only once per process; the first qz-lib call
/// wins (subsequent differing requests are ignored, with a warning). This is the
/// correct mechanism here despite the global state: the chunked pipeline spawns
/// `std::thread::scope` workers that each run rayon work, and a non-global pool's
/// `install()` thread-local does NOT propagate to `std::thread` children — only
/// the global pool does. Idempotent, so per-mode callers may call it again.
pub(crate) fn ensure_global_thread_pool(requested: usize) -> usize {
    let num_threads = if requested == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(8)
    } else {
        requested
    };
    match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
    {
        Ok(()) => num_threads,
        Err(_) => {
            let existing = rayon::current_num_threads();
            if existing != num_threads {
                tracing::warn!(
                    "Requested {} threads but the rayon global pool was already \
                     initialized to {} threads (process-wide; first qz-lib call \
                     wins). Using {}.",
                    num_threads,
                    existing,
                    existing
                );
            }
            existing
        }
    }
}

/// Per-mode capability profile for [`reject_unsupported_common`]. Paired and
/// reference share nearly all of their compress-config gating; this captures the
/// few genuine deltas so the common checks live in exactly one place (adding a new
/// unsupported advanced flag means editing one function, not two).
pub(crate) struct ModeCaps {
    /// Mode name interpolated into error messages (`"paired"` / `"reference"`).
    pub mode: &'static str,
    /// In addition to `Auto`, does this mode accept an explicit `Fqzcomp` quality
    /// compressor? (paired: yes; reference: `Auto` only.)
    pub allows_explicit_fqzcomp: bool,
    /// Is stdout (`"-"`) output permitted? (paired: yes — the v5 encoder is
    /// seek-free; reference: no — it needs a named output file.)
    pub allows_stdout: bool,
}

/// Shared compress-config rejection for the multi-stream modes (paired, reference).
/// Both reject the same set of unsupported options (FASTA, no-quality, lossy, ultra,
/// non-default block size / compressors / advanced flags, stdin inputs); the only
/// real deltas are the allowed quality compressor and whether stdout is permitted,
/// carried in `caps`. Mode-specific extras (e.g. reference's `--reference-index`)
/// are layered at the call site.
pub(crate) fn reject_unsupported_common(a: &CompressConfig, caps: &ModeCaps) -> Result<()> {
    let mode = caps.mode;
    if a.fasta {
        anyhow::bail!("FASTA not supported in {mode} mode");
    }
    if a.no_quality {
        anyhow::bail!("--no-quality not supported in {mode} mode");
    }
    if !matches!(a.quality_mode, QualityMode::Lossless) {
        anyhow::bail!("lossy quality not supported in {mode} mode");
    }
    if a.ultra.is_some() {
        anyhow::bail!("--ultra not supported in {mode} mode");
    }
    let x = &a.advanced;
    // Blocks must stay at the single-end 25 MB default: the decode helper caps each
    // block at `MAX_BLOCK` (64 MiB + 4096); a larger block could become undecodable.
    if x.bsc_block_size_mb != 25 {
        anyhow::bail!("custom --bsc-block-size not supported in {mode} mode");
    }
    let quality_ok = matches!(x.quality_compressor, QualityCompressor::Auto)
        || (caps.allows_explicit_fqzcomp
            && matches!(x.quality_compressor, QualityCompressor::Fqzcomp));
    if !matches!(x.header_compressor, HeaderCompressor::Columnar)
        || !matches!(x.sequence_compressor, SequenceCompressor::Bsc)
        || !quality_ok
    {
        anyhow::bail!("non-default compressors not supported in {mode} mode");
    }
    if x.bsc_static || x.sequence_hints || x.rc_canon {
        anyhow::bail!("advanced options not supported in {mode} mode");
    }
    for p in &a.input {
        if crate::cli::is_stdio_path(p) {
            anyhow::bail!("{mode} mode needs named input files");
        }
    }
    if !caps.allows_stdout && crate::cli::is_stdio_path(&a.output) {
        anyhow::bail!("{mode} mode needs a named output file");
    }
    Ok(())
}

/// Compress FASTQ input(s) into a QZ archive per `args`.
///
/// Threading note: `args.threads` initializes rayon's **process-global** pool on the
/// first qz-lib call in the process; later calls with a different count are ignored
/// (a warning is logged). Long-lived hosts (services, Python sessions) are therefore
/// pinned to the thread count of the first compress/decompress call.
pub fn compress(args: &CompressConfig) -> Result<()> {
    // Honor -t for EVERY mode (single-end, paired, reference) before dispatch —
    // paired/reference don't route through compress_impl::compress, which used
    // to be the only place the pool was built. Idempotent.
    let _ = ensure_global_thread_pool(args.threads);

    // FASTA inputs carry no quality data. Force quality compression off so the
    // Auto quality compressor never tries to encode an empty quality alphabet
    // (which errors on large inputs). Normalizing here protects every entry
    // point — CLI, Python, and direct library use — in one place.
    let normalized;
    let args = if args.fasta && !args.no_quality {
        normalized = CompressConfig {
            no_quality: true,
            ..args.clone()
        };
        &normalized
    } else {
        args
    };

    // Validate input count: single-end (1) or paired-end (2) only.
    match args.input.len() {
        0 => anyhow::bail!("no input file given"),
        1 | 2 => {}
        n => anyhow::bail!("compress accepts one or two input files (got {n})"),
    }

    // Reference mode is opt-in. Run the full reject gate so every unsupported
    // option errors before any compression work (or any output file) happens.
    if args.reference.is_some() {
        // Reference mode accepts single-end (1) or paired (2) inputs. The input
        // count was already validated to be 1 or 2 above; nothing more to gate
        // here besides the shared reject pass.
        reference::reject_unsupported(args)?;
    }
    cluster::reject_unsupported(args)?;

    // Stdio output bypasses the atomic-write path — pipe writes can't be renamed.
    // Paired v5 is seek-free, so paired→stdout is supported here. Reference+stdout
    // is already rejected upstream by `reference::reject_unsupported` (run above
    // when `--reference` is set), so it cannot reach this branch. Cluster+stdout
    // is also rejected upstream by `cluster::reject_unsupported` (allows_stdout:
    // false, run at line 880 above), so a --cluster job can never reach this path.
    if crate::cli::is_stdio_path(&args.output) {
        if args.input.len() == 2 {
            return paired::compress_paired_v5(args);
        }
        return compress_impl::compress(args);
    }

    if args.output.exists() && !args.force {
        anyhow::bail!(
            "Output file already exists: {}\nUse --force to overwrite",
            args.output.display()
        );
    }

    // Compress to a randomized sibling temp file created with O_EXCL, then rename
    // atomically on success. The random name avoids colliding with a concurrent
    // job and, with O_EXCL, refuses to follow a pre-planted symlink (which would
    // truncate an unrelated writable file). A drop guard removes the temp on any
    // early return (error or panic). The per-encoder writers below reopen this
    // path; since we just created it as a regular file at an unpredictable name,
    // that reopen targets our file, not an attacker's symlink.
    let tmp_output = {
        let mut name = args
            .output
            .file_name()
            .ok_or_else(|| {
                anyhow::anyhow!("output path has no file name: {}", args.output.display())
            })?
            .to_os_string();
        name.push(atomic_tmp_suffix());
        args.output.with_file_name(name)
    };

    struct TmpCleanup(std::path::PathBuf);
    impl Drop for TmpCleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    // Reserve the temp with O_EXCL before any encoder opens it by path.
    create_new_for_write(&tmp_output).map_err(|e| {
        anyhow::anyhow!(
            "Failed to create temp output {}: {}",
            tmp_output.display(),
            e
        )
    })?;
    let tmp_guard = TmpCleanup(tmp_output.clone());

    let mut tmp_args = args.clone();
    tmp_args.output = tmp_output.clone();
    if tmp_args.reference.is_some() {
        if tmp_args.input.len() == 1 {
            reference::compress_reference_single(&tmp_args)?; // single-end → v5 archive_type 4
        } else {
            reference::compress_reference(&tmp_args)?; // paired → v5 archive_type 2
        }
    } else if tmp_args.cluster.is_some() {
        if tmp_args.input.len() == 2 {
            cluster::compress_cluster_paired(&tmp_args)?; // paired cluster: R1+R2 → mate-tagged v5
        } else {
            cluster::compress_cluster(&tmp_args)?; // single-end cluster: 1 input → v5 archive_type 3
        }
    } else if tmp_args.input.len() == 2 {
        paired::compress_paired_v5(&tmp_args)?; // writes directly to tmp_args.output (the .tmp path)
    } else {
        compress_impl::compress(&tmp_args)?;
    }

    // Success path: rename tmp → final. Disarm the drop guard first so a rename
    // failure (only e.g. cross-device) re-arms it ourselves and reports clearly.
    std::mem::forget(tmp_guard);
    if let Err(e) = std::fs::rename(&tmp_output, &args.output) {
        let _ = std::fs::remove_file(&tmp_output);
        anyhow::bail!(
            "Failed to rename temp file to {}: {}",
            args.output.display(),
            e
        );
    }
    Ok(())
}

/// Decompress a QZ archive back to FASTQ.
///
/// Threading note: `args.num_threads` initializes rayon's **process-global** pool on
/// the first qz-lib call in the process; later calls with a different count are ignored
/// (a warning is logged). Long-lived hosts (services, Python sessions) are therefore
/// pinned to the thread count of the first compress/decompress call.
pub fn decompress(args: &DecompressConfig) -> Result<()> {
    // Honor -t for EVERY mode before dispatch — paired/reference decode does not
    // route through decompress_impl::decompress (the only place the pool used to
    // be built). Idempotent.
    let _ = ensure_global_thread_pool(args.num_threads);

    // Dispatch: reference and paired archives decode through their own bounded
    // paths (two output files); cluster falls through the force-check then routes;
    // single-end falls through entirely. A malformed-paired archive ERRORS here
    // rather than falling into the single-end parser.
    let atype = v5_archive_type(&args.input)?;
    // Output-count contract: every non-cluster mode takes exactly ONE `-o` — single-end
    // and single-end-reference write it directly; paired and (non-cluster) reference use
    // it as a PREFIX (→ `_R1.fastq`/`_R2.fastq`). Reject extra outputs LOUDLY instead of
    // silently dropping `output[1..]` (the old behavior, which produced misnamed files).
    // Only paired-cluster takes two exact paths, enforced inside `cluster::decompress_cluster`.
    if matches!(atype, Some(0 | 1 | 2 | 4)) && args.output.len() != 1 {
        anyhow::bail!(
            "this archive decompresses to a single -o argument (paired/reference modes use \
             it as a prefix → _R1.fastq/_R2.fastq); got {} -o values",
            args.output.len()
        );
    }
    match atype {
        Some(0) => {} // single-end v5: fall through to the force-check + single-end decode
        Some(1) => return paired::decompress_paired_v5(args),
        Some(2) => return reference::decompress_reference(args),
        Some(4) => return reference::decompress_reference_single(args),
        Some(3) => {} // cluster: fall through so the shared force check runs first
        Some(other) => anyhow::bail!("unknown v5 archive_type {other}"),
        None => {
            if legacy_marker(&args.input)?.is_some() {
                anyhow::bail!("legacy archive no longer supported; recompress with the current qz");
            }
            // Not a legacy marker → legitimate single-end v4: fall through.
        }
    }

    // Refuse to clobber an existing output file unless --force is given,
    // mirroring `compress`. Stdout ('-') is exempt.
    if !args.force {
        for out in &args.output {
            if !crate::cli::is_stdio_path(out) && out.exists() {
                anyhow::bail!(
                    "Output file already exists: {}\nUse --force to overwrite",
                    out.display()
                );
            }
        }
    }
    if atype == Some(3) {
        return cluster::decompress_cluster(args);
    }
    decompress_impl::decompress(args)
}

/// Decompress a **two-mate** archive (paired `archive_type` 1 or paired-reference
/// `archive_type` 2) as a single **interleaved** FASTQ stream — R1/R2 records
/// alternating (`r0/1, r0/2, r1/1, r1/2, …`) — to one `-o` target (a file,
/// `-`/stdout, or gzipped). This is the streamable paired form consumed by
/// aligners like `bwa mem -p`, `bowtie2 --interleaved`, and `strobealign
/// --interleaved`. Single-output archive types (single-end 0, single-reference 4,
/// single-cluster 3) have nothing to interleave and are rejected.
pub fn decompress_interleaved(args: &DecompressConfig) -> Result<()> {
    let _ = ensure_global_thread_pool(args.num_threads);

    if args.output.len() != 1 {
        anyhow::bail!(
            "interleaved decode writes a single stream; pass exactly one -o (file or '-'), got {}",
            args.output.len()
        );
    }
    // Decode reads the trailing footer locator, so the archive INPUT must be seekable.
    // A piped-in archive ('-i -') can't be classified, so reject it with an actionable
    // message rather than the misleading "not a paired archive" the classifier would give.
    if crate::cli::is_stdio_path(&args.input) {
        anyhow::bail!(
            "interleaved decode needs a seekable archive file as input (not stdin); \
             pass the archive path with -i <file.qz>"
        );
    }

    let atype = v5_archive_type(&args.input)?;
    match atype {
        Some(1) | Some(2) => {} // paired / paired-reference: the two-mate modes
        Some(0) | Some(4) => anyhow::bail!(
            "--interleaved applies only to paired archives; this is a single-end archive — decode it normally"
        ),
        Some(3) => anyhow::bail!(
            "--interleaved applies only to paired archives; cluster archives are decoded without --interleaved"
        ),
        Some(other) => anyhow::bail!("unknown v5 archive_type {other}"),
        None => {
            if legacy_marker(&args.input)?.is_some() {
                anyhow::bail!("legacy archive no longer supported; recompress with the current qz");
            }
            anyhow::bail!("--interleaved applies only to paired archives (this is not a paired archive)");
        }
    }

    // The single-sink dispatch renames the atomic temp over the target, so the
    // existing-output guard (mirrored from `decompress`) must run here. Stdout exempt.
    if !args.force {
        let out = &args.output[0];
        if !crate::cli::is_stdio_path(out) && out.exists() {
            anyhow::bail!(
                "Output file already exists: {}\nUse --force to overwrite",
                out.display()
            );
        }
    }

    match atype {
        Some(1) => paired::decompress_paired_interleaved_v5(args),
        Some(2) => reference::decompress_reference_interleaved(args),
        _ => unreachable!("atype validated to 1 or 2 above"),
    }
}

/// Decompress an archive directly to a writer.
///
/// Uses the bounded streaming path when the archive is eligible (per
/// `require_streamable_v5`); otherwise falls back to the mmap path
/// which materialises records and writes them sequentially.
///
/// **BrokenPipe**: if the writer returns `io::ErrorKind::BrokenPipe`
/// (downstream pipe consumer closed), this function returns `Ok(())`.
/// This matches Unix semantics for `qz decompress | head -n 4`. Other
/// I/O errors still propagate.
pub fn decompress_to_writer(
    archive_path: &std::path::Path,
    writer: &mut dyn std::io::Write,
) -> Result<()> {
    // Paired and reference archives produce two output streams; there is no
    // single-stream representation. Reject here rather than silently mis-decoding
    // as single-end. (Legacy paired/reference markers also error, via `legacy_marker`.)
    match v5_archive_type(archive_path)? {
        Some(0) => {} // single-end v5: fall through
        Some(1) => anyhow::bail!(
            "paired archive: use `qz decompress -o <prefix>`; single-stream output unsupported"
        ),
        Some(2) => anyhow::bail!(
            "reference archive: use `qz decompress -o <prefix>`; single-stream output unsupported"
        ),
        Some(3) => anyhow::bail!(
            "cluster archive: decompresses to a named file (use `decompress`), not a streaming writer"
        ),
        Some(4) => anyhow::bail!(
            "single-end reference archives decompress to a named file (use `decompress`), not a writer"
        ),
        Some(other) => anyhow::bail!("unknown v5 archive_type {other}"),
        None => {
            if legacy_marker(archive_path)?.is_some() {
                anyhow::bail!("legacy archive no longer supported; recompress with the current qz");
            }
            // Not a legacy marker → legitimate single-end v4: fall through.
        }
    }
    decompress_impl::decompress_to_writer_impl(archive_path, writer)
}

pub use decompress_impl::{ChunkLayout, VerifyMode, VerifyResult};
pub use merge::{
    merge_archives, merge_archives_to_path, merge_reference_archives_to_path,
    merge_single_end_archives, merge_single_end_archives_to_path,
};
pub use range_compress::{ByteRange, PlanOverride, compress_byte_range, resolve_plan_override};

/// Footer-only chunk layout for NUMA orchestration. Does not decode payloads.
///
/// Reads the v5 front header and directory footer (no block payloads), validates
/// the structure, and returns a [`ChunkLayout`] describing per-chunk read counts,
/// archive type, and reference-global resident size. Used by NUMA-aware callers
/// to assign chunk ranges to per-socket workers.
pub fn read_chunk_layout(archive: &std::path::Path) -> anyhow::Result<ChunkLayout> {
    decompress_impl::read_chunk_layout_impl(archive)
}

/// Test-only helper: parse the v5 front header from `path` and return the
/// byte length of the header (== the footer's `header_end` start). Mirrors how
/// `read_chunk_layout_impl` derives `header_end`: the stored u32 at bytes 4..8
/// is already the complete header size (== `v5_byte_len()`), so we return that
/// directly. Reads only the first 64 bytes (the fixed header is far smaller).
#[cfg(test)]
pub(crate) fn test_front_header_end(path: &std::path::Path) -> Result<u64> {
    use std::io::Read;
    let mut file = std::fs::File::open(path)?;
    let mut hbuf = [0u8; 64];
    let n = file.read(&mut hbuf)?;
    Ok(FixedHeader::parse_v5(&hbuf[..n])?.v5_byte_len() as u64)
}

/// Where a worker must write its chunk-range in the shared final file. The worker
/// independently recomputes both from the ChunkDecodedSizes table and requires a match.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DirectWriteRegion {
    pub offset: u64,
    pub len: u64,
}

/// Distinguishes a direct-write integrity/table failure (table mismatch, overrun,
/// short write) from ordinary I/O errors, so `auto` can clean up and fall back
/// in-process (which decodes correctly without the table) instead of failing.
#[derive(Debug)]
pub struct DirectWriteIntegrityError(pub String);

impl std::fmt::Display for DirectWriteIntegrityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "direct-write integrity: {}", self.0)
    }
}

impl std::error::Error for DirectWriteIntegrityError {}

/// Decode only chunks `[chunk_start, chunk_end)` of a v5 archive into `out_parts`
/// (1 path for single-end, 2 for paired/reference). For NUMA multi-process decode.
///
/// With `regions` empty, each part file is created and written from byte 0 (the
/// part-file assembly path). With `regions` non-empty (single-end OR paired) the
/// worker writes its chunk-range directly into pre-sized shared output file(s) —
/// one [`DirectWriteRegion`] per output (single-end: 1; paired: 2, R1/R2) — at each
/// region's `offset`, after verifying `offset`/`len` against the authoritative
/// per-mate `ChunkDecodedSizes` table and bounding each write to its `len`.
#[allow(clippy::too_many_arguments)]
pub fn decode_chunk_range(
    archive: &std::path::Path,
    chunk_start: u32,
    chunk_end: u32,
    threads: usize,
    gzipped: bool,
    gzip_level: u32,
    out_parts: &[std::path::PathBuf],
    regions: &[DirectWriteRegion],
) -> Result<()> {
    decompress_impl::decode_chunk_range_impl(archive, chunk_start, chunk_end, threads, gzipped, gzip_level, out_parts, regions)
}

/// Conservative worker peak-RSS estimate for the NUMA gate's memory check.
///
/// Returns an upper bound on how much memory one decode worker will use, in
/// bytes, based on `encoding_type` (from [`ChunkLayout::encoding_type`]),
/// `archive_type` (from [`ChunkLayout::archive_type`]), and the reference
/// resident byte total (from [`ChunkLayout::reference_resident_bytes`]).
pub fn decode_peak_rss_bound(
    encoding_type: u8,
    archive_type: u8,
    reference_resident_bytes: u64,
) -> u64 {
    decompress_impl::decode_peak_rss_bound_impl(encoding_type, archive_type, reference_resident_bytes)
}

/// The v5 `archive_type` byte (0=single-end, 1=paired, 2=reference,
/// 3=cluster-reorder, 4=single-end-reference) if `path` is
/// a v5 chunk-major archive, else `None` (v4 / legacy — caller uses the
/// `legacy_marker` probe). Reads only the 4-byte prefix; never errors on a
/// short/non-QZ file (returns None so downstream parsers produce the precise
/// error).
///
/// Stdin (`-`) is not seekable and cannot be re-read by path, so it returns
/// `None` here — the caller then falls through to the single-end stdin-spool
/// path (paired/reference archives are never read from stdin). Without this
/// short-circuit, `File::open("-")` would abort stdin decode/verify.
fn v5_archive_type(path: &std::path::Path) -> Result<Option<u8>> {
    use std::io::Read;
    if crate::cli::is_stdio_path(path) {
        return Ok(None);
    }
    let mut f = std::fs::File::open(path)?;
    let mut p = [0u8; 4];
    let mut got = 0;
    while got < 4 {
        match f.read(&mut p[got..])? {
            0 => break,
            k => got += k,
        }
    }
    if got < 4 || p[0..2] != ARCHIVE_MAGIC || p[2] != ARCHIVE_VERSION_CHUNK_MAJOR {
        return Ok(None);
    }
    Ok(Some(p[3]))
}

/// The legacy paired/reference marker byte (`0xF0`/`0xF1`/`0xF2`) if `path` is a
/// non-v5 archive whose byte 8 carries one — i.e. a legacy v4 paired/reference
/// archive that is no longer decodable and must be rejected with a recompress
/// hint. Returns `None` otherwise.
///
/// This is only ever consulted in the `v5_archive_type == None` branch (the file
/// is not a v5 chunk-major archive). A legitimate single-end **v4** archive lands
/// its `encoding_type` at byte 8 — a small value (0/4/5/6/10), never `0xF0+` — so
/// it correctly returns `None` here and falls through to the single-end decode.
/// Tolerant like `v5_archive_type`: a too-short / unreadable file → `Ok(None)`.
/// stdin (`-`) is not seekable; it returns `None` (single-end stdin spool path).
fn legacy_marker(path: &std::path::Path) -> Result<Option<u8>> {
    use std::io::Read;
    if crate::cli::is_stdio_path(path) {
        return Ok(None);
    }
    let mut f = std::fs::File::open(path)?;
    let mut p = [0u8; 9];
    let mut got = 0;
    while got < 9 {
        match f.read(&mut p[got..])? {
            0 => break,
            k => got += k,
        }
    }
    if got < 9 || p[0..2] != ARCHIVE_MAGIC {
        return Ok(None);
    }
    match p[8] {
        0xF0..=0xF2 => Ok(Some(p[8])),
        _ => Ok(None),
    }
}

/// Whether `path` is a paired QZ archive. Errors on malformed-paired bytes.
/// Public helper so single-output API consumers (e.g. the Python binding) can
/// reject paired archives without reaching into the private `paired` module.
pub fn is_paired_archive(path: &std::path::Path) -> Result<bool> {
    // Answers "is this a v5 paired archive?". A legacy (non-v5) archive is not a
    // v5 paired archive, so it classifies as `Ok(false)` here; the recompress-hint
    // rejection lives in the decompress/verify dispatch (via `legacy_marker`).
    Ok(v5_archive_type(path)? == Some(1))
}

/// Whether `path` is a reference-mode QZ archive (v5 `archive_type == 2`).
/// Public helper so single-output API consumers (e.g. the Python binding) can
/// reject reference archives — which reconstruct two FASTQ files — without
/// reaching into the public `reference` module.
pub fn is_reference_archive(path: &std::path::Path) -> Result<bool> {
    // Answers "is this a v5 reference archive?". Reference archives are now v5
    // (archive_type 2). A legacy (non-v5) reference archive is not a v5 reference
    // archive, so it classifies as `Ok(false)` here; the recompress-hint rejection
    // lives in the decompress/verify dispatch (via `legacy_marker`).
    Ok(v5_archive_type(path)? == Some(2))
}

/// Whether `path` is a PAIRED cluster archive — a v5 `archive_type == 3` archive
/// whose directory carries mate-2 entries, so it reconstructs TWO FASTQ files.
/// Single-end cluster (one output) returns `false`. Public helper so single-output
/// API consumers (the Python binding) can reject the two-output case with a clear
/// message instead of failing deep in `decompress_cluster_paired`.
pub fn is_cluster_paired_archive(path: &std::path::Path) -> Result<bool> {
    if v5_archive_type(path)? != Some(3) {
        return Ok(false);
    }
    // Cluster single-vs-paired is a property of the directory (mate-2 presence),
    // not the front header, so read the footer to classify.
    use std::io::Read;
    let mut f = std::fs::File::open(path)?;
    let mut hb = [0u8; 8];
    f.read_exact(&mut hb)?;
    let header_end = u32::from_le_bytes(hb[4..8].try_into().unwrap()) as u64;
    let dir = chunk_directory::read_v5_footer(path, header_end)?;
    Ok(cluster::decode::is_paired(&dir))
}

/// Internal: a v5 reference-mode archive (paired type 2 OR single-end type 4). Used by
/// the NUMA driver to decide whether to load whole-archive globals before sharding.
/// NOTE: distinct from the PUBLIC `is_reference_archive`, which stays type-2-only
/// (two-output reference) so qz-python does not reject the single-output type-4 archive.
fn is_v5_reference(at: u8) -> bool {
    matches!(at, 2 | 4)
}

/// Whether `path` is a single-end reference QZ archive (v5 `archive_type == 4`).
/// Public helper so single-output API consumers classify it like single-end
/// default (one output) but still load the reference-globals memory bound.
pub fn is_reference_single_archive(path: &std::path::Path) -> Result<bool> {
    Ok(v5_archive_type(path)? == Some(4))
}

/// Verify a QZ archive: decompress all streams and compute CRC32.
pub fn verify(config: &crate::cli::VerifyConfig) -> Result<VerifyResult> {
    // Tri-state classification routes paired-reference / paired / single-end to
    // their respective (deep/fast) verifiers.
    match v5_archive_type(&config.input)? {
        Some(0) => decompress_impl::verify(config), // single-end v5
        Some(1) => paired::verify_paired_v5(&config.input, config.fast),
        Some(2) => reference::verify_reference(&config.input, config.fast),
        Some(3) => cluster::verify_cluster(&config.input, config.fast),
        Some(4) => reference::verify_reference_single(&config.input, config.fast),
        Some(other) => anyhow::bail!("unknown v5 archive_type {other}"),
        None => {
            if legacy_marker(&config.input)?.is_some() {
                anyhow::bail!("legacy archive no longer supported; recompress with the current qz");
            }
            // Not a legacy marker → legitimate single-end v4.
            decompress_impl::verify(config)
        }
    }
}

/// Write variable-length integer
pub(crate) fn write_varint<W: std::io::Write>(
    writer: &mut W,
    mut value: usize,
) -> std::io::Result<()> {
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

        // Reject values that don't fit in usize rather than silently dropping
        // the high bits (the 10th byte, at shift 63, only has room for 1 bit).
        let chunk = (byte & 0x7F) as usize;
        if shift >= usize::BITS as usize || (chunk << shift) >> shift != chunk {
            return None;
        }
        value |= chunk << shift;

        if byte & 0x80 == 0 {
            return Some(value);
        }

        shift += 7;
        if bytes_read >= MAX_VARINT_BYTES {
            return None; // Malformed varint: too many continuation bytes
        }
    }
}

#[cfg(test)]
mod fixed_header_tests {
    use super::*;

    #[test]
    fn encoding_type_roundtrip() {
        for e in [
            EncodingType::Raw,
            EncodingType::RawWithHints,
            EncodingType::RcCanon,
            EncodingType::UltraBigBlock,
            EncodingType::ClusterReorder,
        ] {
            assert_eq!(EncodingType::from_u8(e.to_u8()), Some(e));
        }
        // Removed/unknown encodings decode to None (5 = old delta, 8/9 = old
        // reorder-local).
        for code in [1u8, 2, 3, 5, 7, 8, 9, 255] {
            assert_eq!(
                EncodingType::from_u8(code),
                None,
                "code {code} must be None"
            );
        }
    }

    #[test]
    fn encoding_type_capability_magic_gate() {
        // The named EncodingType capabilities must reproduce the exact magic-number
        // gate the streaming decoder dispatches on, so this stays behavior-preserving.
        for code in 0u8..=12 {
            let e = EncodingType::from_u8(code);
            assert_eq!(
                e.is_some_and(|e| e.is_bounded_streamable()),
                code == 0 || code == 4 || code == 6 || code == 10 || code == 11,
                "stream parity for encoding_type {code}"
            );
            assert_eq!(
                e.is_some_and(|e| e.has_sequence_hints()),
                code == 4,
                "hints {code}"
            );
            assert_eq!(e.is_some_and(|e| e.has_rc_canon()), code == 6, "rc {code}");
        }
    }
}

#[cfg(test)]
mod encoding_type_tests {
    use super::*;

    #[test]
    fn ultra_big_block_encoding_properties() {
        let e = EncodingType::from_u8(10).expect("10 maps to UltraBigBlock");
        assert!(matches!(e, EncodingType::UltraBigBlock));
        // Bounded-streamable: ultra decodes through the streaming path.
        assert!(e.is_bounded_streamable());
        // Raw BSC, no hint byte on the wire — must NOT claim sequence hints.
        assert!(!e.has_sequence_hints());
        // No reverse-complement canonicalization.
        assert!(!e.has_rc_canon());
    }

    #[test]
    fn cluster_reorder_encoding_roundtrips() {
        assert_eq!(EncodingType::from_u8(11), Some(EncodingType::ClusterReorder));
        assert_eq!(EncodingType::ClusterReorder.to_u8(), 11);
        assert!(EncodingType::ClusterReorder.is_bounded_streamable());
        assert!(!EncodingType::ClusterReorder.has_sequence_hints());
        assert!(!EncodingType::ClusterReorder.has_rc_canon());
        assert_eq!(encoding_type_display_name(11), "cluster-reorder");
    }
}

#[cfg(test)]
mod v5_header_tests {
    use super::*;

    #[test]
    fn v5_header_roundtrips_without_lengths() {
        let hdr = FixedHeader {
            encoding_type: EncodingType::Raw.to_u8(),
            flags: FLAG_CONST_LENGTHS,
            archive_type: 0,
            quality_binning_code: 0,
            quality_compressor_code: 1,
            sequence_compressor_code: 1,
            header_compressor_code: 1,
            const_seq_len: 150,
            const_qual_len: 150,
        };
        let mut buf = Vec::new();
        let n = hdr.write_v5(&mut buf).unwrap();
        assert_eq!(n, buf.len());
        assert_eq!(buf[2], ARCHIVE_VERSION_CHUNK_MAJOR);
        assert_eq!(
            read_archive_version(&buf).unwrap(),
            ARCHIVE_VERSION_CHUNK_MAJOR
        );
        let parsed = FixedHeader::parse_v5(&buf).unwrap();
        // num_reads / per-stream lengths are NOT carried by v5 — they live in
        // the directory footer, so only the 12-byte body fields round-trip here.
        assert_eq!(parsed.encoding_type, hdr.encoding_type);
        assert_eq!(parsed.flags, hdr.flags);
        assert_eq!(parsed.const_seq_len, 150);
        assert_eq!(parsed.const_qual_len, 150);
    }

    #[test]
    fn v5_header_archive_type_roundtrips() {
        let mut buf = Vec::new();
        write_v5_archive_header(
            &mut buf,
            0,
            QualityBinning::None,
            QualityCompressor::Bsc,
            HeaderCompressor::Bsc,
            false,
            0,
            0,
            /*archive_type=*/ 1,
        )
        .unwrap();
        assert_eq!(buf[3], 1, "archive_type must be byte 3");
        let h = FixedHeader::parse_v5(&buf).unwrap();
        assert_eq!(h.archive_type, 1);
    }

    #[test]
    fn v5_header_rejects_unknown_archive_type() {
        let mut buf = Vec::new();
        write_v5_archive_header(
            &mut buf,
            0,
            QualityBinning::None,
            QualityCompressor::Bsc,
            HeaderCompressor::Bsc,
            false,
            0,
            0,
            /*archive_type=*/ 9,
        )
        .unwrap();
        assert!(FixedHeader::parse_v5(&buf).is_err());
    }

    #[test]
    fn v5_header_accepts_reference_archive_type() {
        let mut buf = Vec::new();
        write_v5_archive_header(
            &mut buf,
            0,
            QualityBinning::None,
            QualityCompressor::Bsc,
            HeaderCompressor::Bsc,
            false,
            0,
            0,
            /*archive_type=*/ 2,
        )
        .unwrap();
        assert_eq!(buf[3], 2);
        assert_eq!(FixedHeader::parse_v5(&buf).unwrap().archive_type, 2);
    }

    #[test]
    fn v5_header_accepts_4_rejects_5() {
        // archive_type 4 (single-end reference) is now valid.
        let mut ok = Vec::new();
        write_v5_archive_header(
            &mut ok,
            /*encoding_type=*/ 0,
            QualityBinning::None,
            QualityCompressor::Bsc,
            HeaderCompressor::Columnar,
            false,
            /*const_seq_len=*/ 0,
            /*const_qual_len=*/ 0,
            /*archive_type=*/ 4,
        )
        .unwrap();
        assert_eq!(ok[3], 4, "archive_type must be byte 3");
        assert_eq!(FixedHeader::parse_v5(&ok).unwrap().archive_type, 4);

        // archive_type 5 is still out of range.
        let mut bad = Vec::new();
        write_v5_archive_header(
            &mut bad,
            0,
            QualityBinning::None,
            QualityCompressor::Bsc,
            HeaderCompressor::Columnar,
            false,
            0,
            0,
            /*archive_type=*/ 5,
        )
        .unwrap();
        assert!(FixedHeader::parse_v5(&bad).is_err());
    }
}

#[cfg(test)]
mod varint_tests {
    use super::*;

    #[test]
    fn test_read_varint_rejects_overflowing_10_byte_value() {
        // 9 continuation bytes (0x80 = value 0, continue) then a 10th byte whose
        // value (2) does not fit in the single bit remaining at shift 63. This
        // overflows usize and must be rejected, not silently truncated.
        let data = [0x80u8, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x02];
        let mut off = 0;
        assert_eq!(read_varint(&data, &mut off), None);
    }

    #[test]
    fn test_read_varint_accepts_valid_value() {
        let data = [0xACu8, 0x02]; // 44 | (2 << 7) = 300
        let mut off = 0;
        assert_eq!(read_varint(&data, &mut off), Some(300));
    }
}
