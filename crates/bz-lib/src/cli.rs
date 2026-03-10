use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

pub struct CompressConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub working_dir: PathBuf,
    pub threads: usize,
    pub advanced: AdvancedOptions,
}

pub struct DecompressConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub working_dir: PathBuf,
    pub threads: usize,
}

pub struct ExtractConfig {
    pub input: PathBuf,
    pub output_prefix: String,
    pub working_dir: PathBuf,
    pub threads: usize,
}

pub struct VerifyConfig {
    pub input: PathBuf,
    pub threads: usize,
}

/// Result of archive verification.
pub struct VerifyResult {
    /// Total number of BAM records verified.
    pub num_records: u64,
    /// Number of chunks in the archive.
    pub num_chunks: u32,
    /// Alignment stream compressor (0=bsc, 1=zstd).
    pub alignment_compressor: u8,
    /// Aux stream compressor (0=bsc, 1=zstd).
    pub aux_compressor: u8,
    /// CRC32 over decompressed stream data (integrity fingerprint).
    pub crc32: u32,
    /// Total bytes of data hashed for the CRC.
    pub total_bytes: u64,
    /// Time taken in seconds.
    pub elapsed_secs: f64,
}

/// Tunable compression parameters for BZ BAM compression.
///
/// Passed via `--config config.json` on the CLI.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct AdvancedOptions {
    /// Number of BAM records per compression chunk (500K–10M).
    pub chunk_size: usize,

    /// Sub-block size for quality_ctx compression (100K–2M records).
    pub quality_ctx_block_size: usize,

    /// Enable LZP preprocessing in BSC (can help some streams).
    pub use_lzp: bool,

    /// Use adaptive QLFC coder (true) or static (false).
    pub bsc_adaptive: bool,

    /// Quality compressor: 0=quality_ctx (default), 1=bsc.
    pub quality_compressor: u8,

    /// Compressor for alignment-metadata streams (ref_id, pos, mapq, bin, flag,
    /// next_ref_id, next_pos, tlen): 0=bsc (default), 1=zstd.
    pub alignment_compressor: u8,

    /// Compressor for aux tags stream: 0=bsc (default), 1=zstd.
    pub aux_compressor: u8,
}

impl AdvancedOptions {
    /// Validate that all option values are in acceptable ranges.
    pub fn validate(&self) -> Result<()> {
        if self.chunk_size == 0 {
            bail!("chunk_size must be > 0");
        }
        if self.quality_ctx_block_size == 0 {
            bail!("quality_ctx_block_size must be > 0");
        }
        if self.quality_compressor > 1 {
            bail!("quality_compressor must be 0 (quality_ctx) or 1 (bsc)");
        }
        if self.alignment_compressor > 1 {
            bail!("alignment_compressor must be 0 (bsc) or 1 (zstd)");
        }
        if self.aux_compressor > 1 {
            bail!("aux_compressor must be 0 (bsc) or 1 (zstd)");
        }
        Ok(())
    }
}

impl Default for AdvancedOptions {
    fn default() -> Self {
        Self {
            chunk_size: 2_500_000,
            quality_ctx_block_size: 500_000,
            use_lzp: false,
            bsc_adaptive: true,
            quality_compressor: 0,
            alignment_compressor: 0,
            aux_compressor: 0,
        }
    }
}
