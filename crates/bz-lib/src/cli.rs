use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

pub struct CompressConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub working_dir: PathBuf,
    pub advanced: AdvancedOptions,
}

pub struct DecompressConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub working_dir: PathBuf,
}

pub struct ExtractConfig {
    pub input: PathBuf,
    pub output_prefix: String,
    pub working_dir: PathBuf,
    /// Overwrite existing `_R1.qz`/`_R2.qz`/`_SE.qz` outputs. Threaded into the
    /// downstream qz `CompressConfig.force`; without it a re-run bails even when
    /// the CLI was given `--force`.
    pub force: bool,
}

pub struct VerifyConfig {
    pub input: PathBuf,
}

/// Result of archive verification.
pub struct VerifyResult {
    /// Total number of BAM records verified.
    pub num_records: u64,
    /// Number of chunks in the archive.
    pub num_chunks: u32,
    /// CRC32 over decompressed stream data (integrity fingerprint).
    pub crc32: u32,
    /// Total bytes of data hashed for the CRC.
    pub total_bytes: u64,
    /// Time taken in seconds.
    pub elapsed_secs: f64,
}

/// How to flatten the quality of confident (non-variant) bases in the lossy
/// variant-aware reduction. Kept bases always retain full resolution.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum FlattenScheme {
    /// All flattened bases get this single Phred value.
    Single(u8),
    /// Phred < `thresh` → `low`, else `high` (preserves a "bad base" signal).
    TwoBin { thresh: u8, low: u8, high: u8 },
    /// Map to 8 Illumina-style representative bins.
    Coarse8,
}

impl FlattenScheme {
    /// Flatten one raw-Phred value per the scheme.
    pub fn flatten(self, q: u8) -> u8 {
        match self {
            FlattenScheme::Single(v) => v,
            FlattenScheme::TwoBin { thresh, low, high } => {
                if q < thresh {
                    low
                } else {
                    high
                }
            }
            FlattenScheme::Coarse8 => match q {
                0..=2 => 2,
                3..=9 => 6,
                10..=19 => 15,
                20..=24 => 22,
                25..=29 => 27,
                30..=34 => 33,
                35..=39 => 37,
                _ => 40,
            },
        }
    }
}

/// Lossy variant-aware quality reduction (Crumble/Quartz-style). Off by default.
///
/// Quality at confident reference-matching positions is flattened per `scheme`;
/// quality at candidate-variant columns (+ a `window`), per-read mismatches, and
/// inserted bases is kept at full resolution. The classification uses bz's
/// per-position consensus base counts (depth + allele distribution).
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct QualityReduction {
    /// A column is a variant site when a non-consensus allele has ≥ this many reads.
    pub min_alt_count: u16,
    /// …and that allele is ≥ this fraction of the column depth.
    pub min_alt_frac: f32,
    /// Keep full quality within this many reference positions of a variant site.
    pub window: i32,
    /// How to flatten confident bases.
    pub scheme: FlattenScheme,
}

impl QualityReduction {
    /// Preset levels 1–3 (higher = more flattening: stricter variant detection +
    /// narrower window). `scheme` is chosen independently.
    pub fn level(level: u8, scheme: FlattenScheme) -> Self {
        let (min_alt_count, min_alt_frac, window) = match level {
            1 => (2u16, 0.05f32, 20i32), // conservative — keep more
            2 => (3, 0.10, 8),
            _ => (4, 0.20, 3), // aggressive — flatten more
        };
        Self {
            min_alt_count,
            min_alt_frac,
            window,
            scheme,
        }
    }
}

/// Tunable compression parameters for BZ BAM compression.
///
/// Passed via `--config config.json` on the CLI.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct AdvancedOptions {
    /// Number of BAM records per compression chunk (500K–10M).
    pub chunk_size: usize,

    /// Sub-block size for fqz compression (100K–2M records).
    pub fqz_block_size: usize,

    /// Enable LZP preprocessing in BSC (can help some streams).
    pub use_lzp: bool,

    /// Use adaptive QLFC coder (true) or static (false).
    pub bsc_adaptive: bool,

    /// Quality compressor: 0=fqz (default), 1=bsc.
    pub quality_compressor: u8,

    /// Number of chunks to compress simultaneously. **0 = auto** (the default):
    /// the effective window is derived at runtime from core count and available
    /// RAM (see `resolve_pipeline_window` in compress.rs). A non-zero value pins
    /// the window explicitly.
    ///
    /// This is the dominant compress-speed lever. A file has `num_records /
    /// chunk_size` chunks (e.g. ~8 for HG002 chr20 at 2.5M); the window caps how
    /// many compress concurrently, so a window below the chunk count starves the
    /// rayon pool. Measured on HG002 chr20 (18.3M reads, 72 cores): window 4 →
    /// ~10 cores busy / 57.6s; window 8 → ~19 cores / 36s (−38%). The archive
    /// bytes are independent of the window — it is a pure speed/parallelism knob.
    ///
    /// Peak RSS scales with `chunk_size × window`: at chunk_size=2.5M each
    /// in-flight chunk holds ~3-3.7 GB, so window 8 peaks ~30 GB. Auto-selection
    /// keeps the estimate under ~70% of available RAM.
    pub compress_window: usize,

    /// Compression level: a memory/speed preset (all lossless). bz's ratio is
    /// ~flat across levels (fixed 25 MB BSC blocks), so this trades peak memory
    /// and speed, not ratio. `0` = auto (default): balanced level, downshifted
    /// under low RAM. 1 = lowest memory, 2 = balanced/fastest, 3 = largest chunks
    /// (most memory). A non-zero level overrides explicit `chunk_size`/`window`.
    #[serde(default)]
    pub level: u8,

    /// Lossy variant-aware quality reduction. `None` (default) = lossless.
    /// When set, quality at confident reference-matching positions is flattened
    /// while candidate-variant sites keep full resolution. Marks the archive lossy.
    #[serde(default)]
    pub quality_reduction: Option<QualityReduction>,
}

impl AdvancedOptions {
    /// Validate that all option values are in acceptable ranges.
    pub fn validate(&self) -> Result<()> {
        if self.chunk_size == 0 {
            bail!("chunk_size must be > 0");
        }
        if self.fqz_block_size == 0 {
            bail!("fqz_block_size must be > 0");
        }
        if self.quality_compressor > 1 {
            bail!("quality_compressor must be 0 (fqz) or 1 (bsc)");
        }
        // compress_window == 0 is valid: it means "auto" (resolved at runtime
        // from core count + available RAM by resolve_pipeline_window).
        Ok(())
    }
}

impl Default for AdvancedOptions {
    fn default() -> Self {
        Self {
            chunk_size: 2_500_000,
            fqz_block_size: 500_000,
            use_lzp: false,
            bsc_adaptive: true,
            quality_compressor: 0,
            compress_window: 0, // 0 = auto (resolve_pipeline_window picks by cores + RAM)
            level: 0,           // 0 = auto-select a memory/ratio level by available RAM
            quality_reduction: None, // lossless by default
        }
    }
}
