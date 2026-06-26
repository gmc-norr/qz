use std::path::PathBuf;

use serde::{Deserialize, Serialize};
use serde_repr::{Deserialize_repr, Serialize_repr};

/// Core compression configuration.
///
/// Production fields are directly on this struct. Experimental/advanced options
/// (compressor selection, encoding variants) live in `advanced`.
#[derive(Clone)]
pub struct CompressConfig {
    /// Input FASTQ file(s) (one for single-end, two for paired-end)
    pub input: Vec<PathBuf>,
    /// Output QZ archive file
    pub output: PathBuf,
    /// Working directory for temporary files
    pub working_dir: PathBuf,
    /// Number of threads (0 = auto-detect)
    pub threads: usize,
    /// Do not preserve quality values
    pub no_quality: bool,
    /// Input is FASTA format (no quality scores)
    pub fasta: bool,
    /// Quality compression mode
    pub quality_mode: QualityMode,
    /// Ultra compression with optional level (1-3, 0=auto)
    pub ultra: Option<u8>,
    /// Overwrite output if it already exists (default false — refuse to clobber)
    pub force: bool,
    /// Advanced/experimental options (compressor selection, encoding variants, etc.)
    pub advanced: AdvancedOptions,
    /// Opt-in reference-based paired-end mode (compress-time only; decode is
    /// reference-free). `None` = ordinary single/paired-end compression.
    pub reference: Option<ReferenceOptions>,
    /// Opt-in order-drop cluster mode (reorders reads by minimizer cluster; output is a
    /// permutation of input). `None` = ordinary order-preserving compression.
    pub cluster: Option<ClusterOptions>,
    /// Reference mode only: require a usable prebuilt strobealign index. When `true`,
    /// a missing / stale / corrupt sidecar is a hard error instead of an auto-rebuild
    /// (the CLI sets this when `--build-index` is NOT given). Default `false` keeps the
    /// library's transparent auto-build for API consumers (qz-python) and tests.
    pub require_prebuilt_index: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ClusterLevel { Fast, #[default] Balanced, Max }
impl ClusterLevel {
    pub fn zstd_level(self) -> i32 { match self { Self::Fast => 12, Self::Balanced => 16, Self::Max => 19 } }
}

#[derive(Clone, Debug, Default)]
pub struct ClusterOptions { pub level: ClusterLevel }

/// Options for reference-based paired-end compression.
///
/// The reference FASTA is used only at compress time to derive a consensus and
/// encode reads as edits against it; it is never stored in the archive and is
/// never needed to decode.
#[derive(Clone, Debug)]
pub struct ReferenceOptions {
    /// Reference FASTA used at compress time only (never stored, never needed to decode).
    pub reference: std::path::PathBuf,
    /// Optional prebuilt strobealign index (validated against `reference`). Reserved; rejected in v1.
    pub reference_index: Option<std::path::PathBuf>,
    /// Faster mapping via sparser syncmer seeds (s=12 vs the read-length default):
    /// ~5–10% faster compress at any core count, lossless, ~ratio-neutral on
    /// high-coverage data; may add a few literal fallbacks on low-coverage/divergent
    /// inputs. The sparser index is cached separately from the default one.
    pub reference_fast: bool,
    /// Encode-pipeline depth: how many chunks are encoded CONCURRENTLY before
    /// being written in order. Higher = faster compress but more peak memory
    /// (~`window` × chunk working set). The CLI default is 4 (the measured sweet
    /// spot); lower it to save memory, raise it for a little more speed. Output
    /// is byte-identical regardless of the value.
    pub reference_window: usize,
}

impl Default for CompressConfig {
    fn default() -> Self {
        Self {
            input: Vec::new(),
            output: PathBuf::new(),
            working_dir: PathBuf::from("."),
            threads: 0,
            no_quality: false,
            fasta: false,
            quality_mode: QualityMode::Lossless,
            ultra: None,
            force: false,
            advanced: AdvancedOptions::default(),
            reference: None,
            cluster: None,
            require_prebuilt_index: false,
        }
    }
}

/// Advanced/experimental compression options.
///
/// These control compressor selection, encoding variants, and other
/// features not exposed in the production CLI. Used by integration tests and
/// benchmark code.
#[derive(Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct AdvancedOptions {
    /// Quality score compressor
    pub quality_compressor: QualityCompressor,
    /// Sequence compressor
    pub sequence_compressor: SequenceCompressor,
    /// Header compressor
    pub header_compressor: HeaderCompressor,
    /// Use BSC static coder instead of adaptive
    pub bsc_static: bool,
    /// Reverse-complement canonicalization
    pub rc_canon: bool,
    /// Prepend a syncmer-derived hint byte before each read's sequence
    pub sequence_hints: bool,
    /// Number of reads per fqzcomp quality sub-block (default 500K). Kept equal to
    /// `codecs::FQZCOMP_SUB_CHUNK` so each default single-end fqz block is exactly
    /// one sub-chunk (see that const for the parallelism implication). The field
    /// name is retained from the legacy `--quality-ctx-block-size` flag.
    pub quality_ctx_block_size: usize,
    /// BSC block size in MB for parallel compression (default 25)
    pub bsc_block_size_mb: usize,
    /// Records per chunk for chunked compression (default 2_500_000)
    pub chunk_records: usize,
    /// Number of chunks to compress simultaneously (default 4).
    /// Each concurrent chunk submits its BSC blocks to the shared rayon pool,
    /// keeping more threads busy. A window ≥ 2 raises 72-core utilisation from
    /// ~49% to ~97% (35 blocks/chunk × window concurrent tasks). The default
    /// of 4 trades extra RAM for stable saturation when chunk sizes are uneven.
    pub compress_window: usize,
}

impl Default for AdvancedOptions {
    fn default() -> Self {
        Self {
            quality_compressor: QualityCompressor::Auto,
            sequence_compressor: SequenceCompressor::Bsc,
            header_compressor: HeaderCompressor::Columnar,
            bsc_static: false,
            rc_canon: false,
            sequence_hints: false,
            quality_ctx_block_size: 500_000,
            bsc_block_size_mb: 25,
            chunk_records: 2_500_000,
            compress_window: 4,
        }
    }
}

#[derive(Clone)]
pub struct DecompressConfig {
    /// Input QZ archive
    pub input: PathBuf,
    /// Output FASTQ file(s)
    pub output: Vec<PathBuf>,
    /// Working directory for temporary files
    pub working_dir: PathBuf,
    /// Number of threads
    pub num_threads: usize,
    /// Output gzipped FASTQ
    pub gzipped: bool,
    /// Gzip compression level (0-9)
    pub gzip_level: u32,
    /// Overwrite the output file if it already exists.
    pub force: bool,
}

#[derive(Clone)]
pub struct VerifyConfig {
    /// Input QZ archive
    pub input: PathBuf,
    /// Working directory for temporary files
    pub working_dir: PathBuf,
    /// Number of threads
    pub num_threads: usize,
    /// Fast mode: only verify per-block CRC32 checksums without invoking BSC
    /// decompression. Catches bit-rot in O(IO + CRC) time instead of O(BWT),
    /// but does not reconstruct the FASTQ output (so the deep-verify CRC32
    /// over the decompressed FASTQ bytes is not computed).
    pub fast: bool,
}

impl Default for VerifyConfig {
    fn default() -> Self {
        Self {
            input: PathBuf::new(),
            working_dir: PathBuf::from("."),
            num_threads: 0,
            fast: false,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize_repr, Deserialize_repr)]
#[repr(u8)]
pub enum QualityMode {
    /// Lossless quality preservation
    Lossless = 0,
    /// Illumina 8-level binning
    IlluminaBin = 1,
    /// Illumina 4-level binning (more aggressive)
    Illumina4 = 2,
    /// Binary thresholding
    Binary = 3,
    /// Discard quality scores entirely
    Discard = 5,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize_repr, Deserialize_repr)]
#[repr(u8)]
pub enum QualityCompressor {
    /// BSC/BWT (best raw byte-stream compression)
    Bsc = 1,
    /// fqzcomp context-modeled compression
    Fqzcomp = 3,
    /// Pick automatically: Fqzcomp for lossless inputs >=100k reads
    /// (single-end/paired), else Bsc. Never stored in the archive — resolved to a
    /// concrete choice at compress time.
    Auto = 5,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize_repr, Deserialize_repr)]
#[repr(u8)]
pub enum SequenceCompressor {
    /// BSC on raw ASCII sequences (best compression, default)
    Bsc = 1,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize_repr, Deserialize_repr)]
#[repr(u8)]
pub enum HeaderCompressor {
    /// BSC on raw headers
    Bsc = 1,
    /// Columnar encoding: parse Illumina fields into typed columns + parallel BSC (default)
    Columnar = 3,
}

pub fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8)
}

/// Check if a path represents stdin/stdout (the `-` convention).
pub fn is_stdio_path(p: &std::path::Path) -> bool {
    p.as_os_str() == "-"
}

#[cfg(test)]
mod cluster_cfg_tests {
    use super::*;
    #[test]
    fn cluster_options_default_is_balanced() {
        let c = CompressConfig::default();
        assert!(c.cluster.is_none());
        assert_eq!(ClusterLevel::default(), ClusterLevel::Balanced);
        assert_eq!(ClusterLevel::Balanced.zstd_level(), 16);
        assert_eq!(ClusterLevel::Fast.zstd_level(), 12);
        assert_eq!(ClusterLevel::Max.zstd_level(), 19);
    }
}
