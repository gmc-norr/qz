/// Ultra compression: routes through the unified `compress_chunked` pipeline,
/// emitting `UltraBigBlock` (encoding_type=10) archives that are fully
/// bounded-streaming decodable. All reorder machinery (encoding_type 8/9) has
/// been removed; those archives are no longer supported.
use anyhow::Result;
use tracing::info;

#[allow(unused_imports)]
use super::*;

// ── Ultra compression levels ──────────────────────────────────────────────

/// All parameters for one ultra compression level. This is the single source
/// of truth shared by `compress_ultra` (which builds `ChunkedMode::Ultra` from
/// these fields) and `estimate_peak_memory` (which models RAM from them).
#[derive(Debug, Clone, Copy)]
pub(super) struct UltraLevel {
    pub level: u8,
    /// Number of FASTQ records per streaming chunk sent to BSC.
    pub chunk_records: usize,
    /// BSC input block size in mebibytes (controls BWT context window).
    pub bsc_block_mb: usize,
    /// Quality stream sub-block size in records (inner BSC block for qualities).
    pub quality_sub_block: usize,
    /// Number of chunks compressed concurrently (pipeline window). This is the
    /// real concurrency that bounds compress-time peak memory.
    pub compress_window: usize,
}

const ULTRA_LEVELS: [UltraLevel; 3] = [
    // Level 1: Fast — smaller chunks, 2-chunk window
    UltraLevel {
        level: 1,
        chunk_records: 2_500_000,
        bsc_block_mb: 188,
        quality_sub_block: 250_000,
        compress_window: 2,
    },
    // Level 2: Balanced — medium chunks, 2-chunk window
    UltraLevel {
        level: 2,
        chunk_records: 3_500_000,
        bsc_block_mb: 375,
        quality_sub_block: 500_000,
        compress_window: 2,
    },
    // Level 3: Good compression — large chunks, sequential (window=1)
    UltraLevel {
        level: 3,
        chunk_records: 5_000_000,
        bsc_block_mb: 750,
        quality_sub_block: 500_000,
        compress_window: 1,
    },
];

fn available_memory_bytes() -> Option<u64> {
    let contents = std::fs::read_to_string("/proc/meminfo").ok()?;
    for line in contents.lines() {
        if line.starts_with("MemAvailable:") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                return parts[1].parse::<u64>().ok().map(|kb| kb * 1024);
            }
        }
    }
    None
}

fn estimate_peak_memory(level: &UltraLevel) -> u64 {
    // Conservative heuristic upper bound for `--ultra 0` auto-selection only.
    // NOTE: the streaming engine compresses ultra ONE big block at a time
    // (`compress_worker_count(is_ultra=true) == 1`), so real peak RSS is bounded by a
    // single big block's libsais workspace — NOT `chunk_records × compress_window`. This
    // product intentionally over-estimates; `compress_window` here is a sizing factor for
    // the gate, not the pipeline's actual concurrency.
    let per_read_bytes: u64 = 1500;
    let io_overhead: u64 = 2 * 1024 * 1024 * 1024;
    level.chunk_records as u64 * per_read_bytes * level.compress_window as u64 + io_overhead
}

fn auto_select_level() -> u8 {
    let available = match available_memory_bytes() {
        Some(bytes) => bytes,
        None => {
            info!("Could not detect available memory, defaulting to ultra level 2");
            return 2;
        }
    };

    let budget = available * 80 / 100;
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8);

    // Cap auto-select at 3. Empirically on 10M-read Illumina WGS, levels 4/5
    // produce identical archives to level 3 but decompress 2-3x slower (larger
    // BSC blocks → fewer blocks → less decompress parallelism). Users who want
    // higher levels can opt in explicitly; auto picks the empirical sweet spot.
    let max_level_for_cores: u8 = if cores < 4 {
        1
    } else if cores < 8 {
        2
    } else {
        3
    };

    let mut selected = 1u8;
    for level in &ULTRA_LEVELS {
        if estimate_peak_memory(level) <= budget && level.level <= max_level_for_cores {
            selected = level.level;
        }
    }

    info!(
        "Auto-selected ultra level {} (available: {:.1} GB, budget: {:.1} GB, cores: {})",
        selected,
        available as f64 / 1e9,
        budget as f64 / 1e9,
        cores,
    );

    selected
}

/// Resolve a requested ultra level to a concrete level in `1..=3`.
///
/// `0` means auto (memory/core-based selection, already capped at 3). Anything
/// outside `0..=3` is rejected — we do NOT clamp, so an explicit out-of-range
/// level surfaces as an error rather than silently coercing.
pub(super) fn resolve_ultra_level(requested: u8) -> Result<u8> {
    match requested {
        0 => Ok(auto_select_level()), // already caps at <=3
        1..=3 => Ok(requested),
        other => anyhow::bail!("--ultra level must be 1, 2, 3, or omitted for auto (got {other})"),
    }
}

/// Ultra compression: looks up the resolved level (1..=3) in `ULTRA_LEVELS`
/// (the single source of truth for all per-level parameters) and routes through
/// the unified `compress_chunked` pipeline, emitting an `UltraBigBlock`
/// (bounded-streamable) archive — no reorder.
pub(super) fn compress_ultra(args: &CompressConfig, level: u8) -> Result<()> {
    let params = &ULTRA_LEVELS[(level - 1) as usize]; // level validated to 1..=3 by resolve_ultra_level
    info!(
        "Ultra level {}: chunk_records={}M, bsc_block={}MB, quality_sub_block={}K, compress_window={}",
        params.level,
        params.chunk_records / 1_000_000,
        params.bsc_block_mb,
        params.quality_sub_block / 1000,
        params.compress_window,
    );
    // check_cuda_vram is called inside compress_chunked once bsc_block_mb is
    // resolved from the ChunkedMode::Ultra variant — no need to duplicate it here.
    super::compress_impl::compress_chunked(args, ultra_mode(level))
}

/// Build the `ChunkedMode::Ultra` for a resolved level (1..=3) from `ULTRA_LEVELS`
/// (the single source of truth).
fn ultra_mode(level: u8) -> super::compress_impl::ChunkedMode {
    let params = &ULTRA_LEVELS[(level - 1) as usize];
    super::compress_impl::ChunkedMode::Ultra {
        bsc_block_mb: params.bsc_block_mb,
        chunk_records: params.chunk_records,
        quality_sub_block: params.quality_sub_block,
        compress_window: params.compress_window,
    }
}
