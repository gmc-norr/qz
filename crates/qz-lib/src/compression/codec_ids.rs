//! The single codec-id namespace for a v5 archive — used by BOTH the per-entry
//! footer-directory `codec` byte and the single-end front-header codec bytes.
//!
//! These are the role-local `codec` field values stored in each directory entry
//! (`ChunkDirEntry`) of a v5 archive footer, and the SAME values written into the
//! front header's `quality_compressor_code` / `sequence_compressor_code` /
//! `header_compressor_code` fields. One byte value means one codec everywhere —
//! across `archive_type` (single-end, paired, reference) AND across the
//! front-header vs directory fields (e.g. fqzcomp quality is `6` in both). Encode
//! maps compressor enums in via [`quality_compressor_to_codec_id`] /
//! [`sequence_compressor_to_codec_id`] / [`header_compressor_to_codec_id`]; decode
//! maps back via [`codec_id_to_quality_compressor`] /
//! [`codec_id_to_header_compressor`]. `paired::format` and `reference::format`
//! validate directory entries against these constants.
//!
//! Single-end decode reads its codecs from the front-header copy; the directory
//! `codec` byte is advisory for single-end (uniformity-checked, never dispatched
//! on). Paired/reference dispatch on the directory copy and ignore the
//! front-header copy. Both copies carry identical values.

use crate::cli::{HeaderCompressor, QualityCompressor, SequenceCompressor};
use anyhow::Result;

/// BSC (BWT + adaptive QLFC) — the default backend for most roles.
pub(crate) const CODEC_BSC: u8 = 1;
/// Reserved: was the OpenZL backend (removed). Kept reserved so the codec
/// namespace never re-assigns value 2 to a different codec; the reverse mappers
/// surface a dedicated "removed" error for it.
pub(crate) const CODEC_OPENZL: u8 = 2;
/// Columnar header codec (used for the R1/R2 header roles).
pub(crate) const CODEC_COLUMNAR: u8 = 3;
/// Reserved: was the context-adaptive `quality_ctx` quality codec (removed —
/// fqzcomp is now the only context-modeled quality backend). Kept reserved so the
/// codec namespace never re-assigns value 4; the reverse mapper surfaces a
/// dedicated "removed" error for it.
pub(crate) const CODEC_QUALITY_CTX: u8 = 4;
/// Packed 2-bit consensus backing (reference mode only).
pub(crate) const CODEC_PACKED_CONSENSUS: u8 = 5;
/// fqzcomp quality codec (htscodecs) — quality roles in all modes.
pub(crate) const CODEC_FQZCOMP: u8 = 6;
/// zstd long-distance-matching codec — cluster-mode sequence role ONLY (the removed
/// quality/header zstd codes 0/2 stay rejected). Level + windowLog live in the block
/// payload prelude, not the directory.
pub(crate) const CODEC_ZSTD_LONG: u8 = 7;

/// Map a **resolved** quality compressor to its directory `codec` id (`Auto` must be
/// resolved before storage). Fqzcomp is `CODEC_FQZCOMP` = 6 in BOTH the front-header
/// quality-compressor byte and the directory entry — one shared codec-id namespace.
pub(crate) fn quality_compressor_to_codec_id(qc: QualityCompressor) -> u8 {
    match qc {
        QualityCompressor::Bsc => CODEC_BSC,
        QualityCompressor::Fqzcomp => CODEC_FQZCOMP,
        QualityCompressor::Auto => unreachable!("Auto must be resolved before storage"),
    }
}

/// Map a sequence compressor to its directory `codec` id.
pub(crate) fn sequence_compressor_to_codec_id(sc: SequenceCompressor) -> u8 {
    match sc {
        SequenceCompressor::Bsc => CODEC_BSC,
    }
}

/// Map a header compressor to its directory `codec` id.
pub(crate) fn header_compressor_to_codec_id(hc: HeaderCompressor) -> u8 {
    match hc {
        HeaderCompressor::Bsc => CODEC_BSC,
        HeaderCompressor::Columnar => CODEC_COLUMNAR,
    }
}

/// Map a codec id back to its quality compressor (single-end front-header decode).
/// The removed zstd/OpenZL/quality_ctx backends keep dedicated error messages.
pub(crate) fn codec_id_to_quality_compressor(code: u8) -> Result<QualityCompressor> {
    match code {
        0 => anyhow::bail!("zstd quality backend has been removed (code 0)"),
        CODEC_BSC => Ok(QualityCompressor::Bsc),
        CODEC_OPENZL => anyhow::bail!("OpenZL quality backend has been removed (code 2)"),
        CODEC_QUALITY_CTX => {
            anyhow::bail!("quality_ctx quality backend has been removed (code 4) — recompress")
        }
        CODEC_FQZCOMP => Ok(QualityCompressor::Fqzcomp),
        _ => anyhow::bail!("invalid quality compressor codec id: {code}"),
    }
}

/// Map a codec id back to its header compressor (single-end front-header decode).
pub(crate) fn codec_id_to_header_compressor(code: u8) -> Result<HeaderCompressor> {
    match code {
        0 => anyhow::bail!("zstd header backend has been removed (code 0)"),
        CODEC_BSC => Ok(HeaderCompressor::Bsc),
        CODEC_OPENZL => anyhow::bail!("OpenZL header backend has been removed (code 2)"),
        CODEC_COLUMNAR => Ok(HeaderCompressor::Columnar),
        _ => anyhow::bail!("invalid header compressor codec id: {code}"),
    }
}

#[cfg(test)]
mod cluster_codec_tests {
    use super::*;

    #[test]
    fn zstd_long_codec_id_is_seven() {
        assert_eq!(CODEC_ZSTD_LONG, 7);
    }
}
