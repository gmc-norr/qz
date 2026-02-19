use std::path::PathBuf;

use serde::{Serialize, Deserialize};
use serde_repr::{Serialize_repr, Deserialize_repr};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize_repr, Deserialize_repr)]
#[repr(u8)]
pub enum ReorderMode {
    /// Sort within each 5M-record chunk (fast, bounded memory, local ordering only)
    Local = 0,
    /// Two-pass bucket sort across entire file (better ordering, bounded memory)
    Global = 1,
}

/// Core compression configuration.
///
/// Production fields are directly on this struct. Experimental/advanced options
/// (compressor selection, encoding variants, reordering) live in `advanced`.
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
    /// Ultra compression with optional level (1-5, 0=auto)
    pub ultra: Option<u8>,
    /// Advanced/experimental options (compressor selection, encoding variants, etc.)
    pub advanced: AdvancedOptions,
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
            advanced: AdvancedOptions::default(),
        }
    }
}

/// Advanced/experimental compression options.
///
/// These control compressor selection, encoding variants, reordering, and other
/// features not exposed in the production CLI. Used by integration tests and
/// benchmark code.
#[derive(Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct AdvancedOptions {
    /// Enable positional quality modeling (experimental)
    pub quality_modeling: bool,
    /// Enable quality delta encoding between adjacent reads (experimental)
    pub quality_delta: bool,
    /// Enable zstd dictionary training
    pub dict_training: bool,
    /// Dictionary size in KB
    pub dict_size: usize,
    /// Compression level for legacy zstd mode (1-22)
    pub compression_level: i32,
    /// Quality score compressor
    pub quality_compressor: QualityCompressor,
    /// Sequence compressor
    pub sequence_compressor: SequenceCompressor,
    /// Header compressor
    pub header_compressor: HeaderCompressor,
    /// Use BSC static coder instead of adaptive
    pub bsc_static: bool,
    /// Use 2-bit sequence encoding (4 bases per byte + N-mask bitmap)
    pub twobit: bool,
    /// Use template-based header encoding
    pub header_template: bool,
    /// Reorder reads by content similarity
    pub reorder: Option<ReorderMode>,
    /// Reverse-complement canonicalization
    pub rc_canon: bool,
    /// Prepend a syncmer-derived hint byte before each read's sequence
    pub sequence_hints: bool,
    /// Inline delta encoding against cached similar reads
    pub sequence_delta: bool,
    /// Local reordering with delta encoding
    pub local_reorder: bool,
    /// Deprecated: use ultra level 2 instead
    pub fast_ultra: bool,
    /// Number of reads per quality_ctx sub-block (default 500K)
    pub quality_ctx_block_size: usize,
    /// BSC block size in MB for parallel compression (default 25)
    pub bsc_block_size_mb: usize,
    /// Records per chunk for chunked compression (default 2_500_000)
    pub chunk_records: usize,
}

impl Default for AdvancedOptions {
    fn default() -> Self {
        Self {
            quality_modeling: false,
            quality_delta: false,
            dict_training: false,
            dict_size: 64,
            compression_level: 3,
            quality_compressor: QualityCompressor::Bsc,
            sequence_compressor: SequenceCompressor::Bsc,
            header_compressor: HeaderCompressor::Columnar,
            bsc_static: false,
            twobit: false,
            header_template: false,
            rc_canon: false,
            reorder: None,
            sequence_hints: false,
            sequence_delta: false,
            local_reorder: false,
            fast_ultra: false,
            quality_ctx_block_size: 500_000,
            bsc_block_size_mb: 25,
            chunk_records: 2_500_000,
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
}

#[derive(Clone)]
pub struct VerifyConfig {
    /// Input QZ archive
    pub input: PathBuf,
    /// Working directory for temporary files
    pub working_dir: PathBuf,
    /// Number of threads
    pub num_threads: usize,
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
    /// QVZ lossy compression
    Qvz = 4,
    /// Discard quality scores entirely
    Discard = 5,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize_repr, Deserialize_repr)]
#[repr(u8)]
pub enum QualityCompressor {
    /// Legacy zstd (faster, ~20% larger)
    Zstd = 0,
    /// BSC/BWT (best compression, default)
    Bsc = 1,
    /// OpenZL format-aware compression
    OpenZl = 2,
    /// fqzcomp context-modeled compression
    Fqzcomp = 3,
    /// Context-adaptive range coding
    QualityCtx = 4,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize_repr, Deserialize_repr)]
#[repr(u8)]
pub enum SequenceCompressor {
    /// Legacy zstd with 2-bit encoding + N-mask
    Zstd = 0,
    /// BSC on raw ASCII sequences (best compression, default)
    Bsc = 1,
    /// OpenZL format-aware compression
    OpenZl = 2,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize_repr, Deserialize_repr)]
#[repr(u8)]
pub enum HeaderCompressor {
    /// Legacy template + zstd
    Zstd = 0,
    /// BSC on raw headers
    Bsc = 1,
    /// OpenZL format-aware compression
    OpenZl = 2,
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
