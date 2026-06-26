/// FFI bindings to libbsc (Block-Sorting Compressor)
///
/// libbsc is the compression library used by official SPRING.
/// This module provides safe Rust wrappers around the C API.
///
/// Source: https://github.com/IlyaGrebnov/libbsc
use anyhow::{Context, Result};
use flate2::Crc;
use libc::{c_int, c_uchar};
use std::sync::Once;

// OpenMP thread control
unsafe extern "C" {
    fn omp_set_num_threads(num_threads: c_int);
}

// Error codes from libbsc
const LIBBSC_NO_ERROR: c_int = 0;
const _LIBBSC_BAD_PARAMETER: c_int = -1;
const _LIBBSC_NOT_ENOUGH_MEMORY: c_int = -2;
const _LIBBSC_NOT_COMPRESSIBLE: c_int = -3;
const _LIBBSC_NOT_SUPPORTED: c_int = -4;
const _LIBBSC_UNEXPECTED_EOB: c_int = -5;
const _LIBBSC_DATA_CORRUPT: c_int = -6;

// Compression parameters (matching official SPRING defaults)
const LIBBSC_BLOCKSORTER_BWT: c_int = 1; // Use BWT for compatibility
const _LIBBSC_BLOCKSORTER_ST7: c_int = 7; // ST7 is what SPRING uses
const LIBBSC_CODER_QLFC_STATIC: c_int = 1;
const LIBBSC_CODER_QLFC_ADAPTIVE: c_int = 2;

// Features (bitmask)
const LIBBSC_FEATURE_FASTMODE: c_int = 1;
const LIBBSC_FEATURE_MULTITHREADING: c_int = 2;
#[cfg(feature = "cuda")]
const LIBBSC_FEATURE_CUDA: c_int = 8;

/// Runtime GPU gate (cuda builds only). A `--features cuda` binary uses the GPU BWT by
/// default; `QZ_GPU=off` (set directly, or by the `--gpu off` CLI flag — which also reaches
/// NUMA worker processes through the inherited environment) forces the CPU BWT path instead.
/// Recognized falsey values: `off`, `0`, `false`. Read once and cached; the CLI sets the env
/// before any compression, so the cache always observes the final value.
#[cfg(feature = "cuda")]
pub(crate) fn gpu_enabled() -> bool {
    use std::sync::OnceLock;
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        !matches!(
            std::env::var("QZ_GPU").ok().as_deref(),
            Some("off") | Some("0") | Some("false")
        )
    })
}

/// Returns LIBBSC_FEATURE_CUDA only when the `cuda` feature is built AND the GPU is enabled at
/// runtime (`QZ_GPU` / `--gpu`); 0 otherwise. Output is byte-identical either way — the gate
/// only selects the BWT backend, never the archive bytes.
#[inline]
fn cuda_feature() -> c_int {
    #[cfg(feature = "cuda")]
    {
        if gpu_enabled() {
            LIBBSC_FEATURE_CUDA
        } else {
            0
        }
    }
    #[cfg(not(feature = "cuda"))]
    {
        0
    }
}

/// True iff a CUDA device is present and its memory is queryable. Used by `--gpu require` to
/// fail loud at startup instead of silently falling back to the CPU BWT. cuda builds only.
#[cfg(feature = "cuda")]
pub fn cuda_device_available() -> bool {
    unsafe extern "C" {
        fn cudaMemGetInfo(free: *mut usize, total: *mut usize) -> c_int;
    }
    let (mut free, mut total) = (0usize, 0usize);
    unsafe { cudaMemGetInfo(&mut free, &mut total) == 0 && free > 0 }
}

/// Query GPU VRAM and warn if it's too small for BWT at the given block size.
/// `max_block_bytes` is the largest single BSC block that will be passed to BWT.
/// libcubwt needs ~20.5× the input length in device memory for forward BWT.
#[cfg(feature = "cuda")]
pub fn check_cuda_vram(max_block_bytes: usize) {
    use tracing::warn;

    // GPU disabled at runtime (`QZ_GPU=off` / `--gpu off`) → CPU BWT, no VRAM concern.
    if !gpu_enabled() {
        return;
    }

    unsafe extern "C" {
        fn cudaMemGetInfo(free: *mut usize, total: *mut usize) -> c_int;
    }

    let mut free: usize = 0;
    let mut total: usize = 0;
    let status = unsafe { cudaMemGetInfo(&mut free, &mut total) };
    if status != 0 {
        warn!(
            "Could not query GPU VRAM (cudaMemGetInfo returned {}). GPU acceleration may fail.",
            status
        );
        return;
    }

    // libcubwt allocates with max_length = n + n/32, then ~20.5× that
    let max_length = max_block_bytes + max_block_bytes / 32;
    let estimated_vram = (max_length as f64 * 20.5) as usize;

    if free < estimated_vram {
        warn!(
            "GPU VRAM may be insufficient: {:.0} MB free, ~{:.0} MB needed for {:.0} MB blocks. \
             BWT will fall back to CPU for large blocks.",
            free as f64 / 1_048_576.0,
            estimated_vram as f64 / 1_048_576.0,
            max_block_bytes as f64 / 1_048_576.0,
        );
    }
}

#[cfg(not(feature = "cuda"))]
pub fn check_cuda_vram(_max_block_bytes: usize) {}

// Header size
const LIBBSC_HEADER_SIZE: usize = 28;

// Block size (matching official SPRING: 25 MB blocks for better compression)
const BSC_BLOCK_SIZE: usize = 25 * 1024 * 1024; // 25 MB

#[link(name = "libbsc", kind = "static")]
unsafe extern "C" {
    /// Initialize libbsc library
    fn bsc_init(features: c_int) -> c_int;

    /// Compress a memory block
    fn bsc_compress(
        input: *const c_uchar,
        output: *mut c_uchar,
        n: c_int,
        lzp_hash_size: c_int,
        lzp_min_len: c_int,
        block_sorter: c_int,
        coder: c_int,
        features: c_int,
    ) -> c_int;

    /// Get information about compressed block
    fn bsc_block_info(
        block_header: *const c_uchar,
        header_size: c_int,
        p_block_size: *mut c_int,
        p_data_size: *mut c_int,
        features: c_int,
    ) -> c_int;

    /// Decompress a memory block
    fn bsc_decompress(
        input: *const c_uchar,
        input_size: c_int,
        output: *mut c_uchar,
        output_size: c_int,
        features: c_int,
    ) -> c_int;
}

/// Ensure bsc_init is called exactly once.
///
/// Init does NOT include `LIBBSC_FEATURE_MULTITHREADING` — that would pin process-wide
/// OpenMP state and break the project's "rayon owns parallelism" invariant
/// (CLAUDE.md). Per-call enablement on the few MT entry points works without it.
///
/// `cuda_feature()` evaluates to `LIBBSC_FEATURE_CUDA` only when the `cuda`
/// Cargo feature is enabled at compile time (zero otherwise). Users opt in to
/// GPU acceleration via the build flag, so initializing the CUDA path at
/// startup is the intended behavior for CUDA builds. Non-CUDA builds pay no
/// cost — the flag literal is `0`.
static BSC_INIT: Once = Once::new();
static BSC_INIT_RESULT: std::sync::atomic::AtomicI32 = std::sync::atomic::AtomicI32::new(i32::MIN);

fn ensure_initialized() -> Result<()> {
    BSC_INIT.call_once(|| {
        let features = LIBBSC_FEATURE_FASTMODE | cuda_feature();
        let result = unsafe { bsc_init(features) };
        BSC_INIT_RESULT.store(result, std::sync::atomic::Ordering::SeqCst);
    });
    let result = BSC_INIT_RESULT.load(std::sync::atomic::Ordering::SeqCst);
    if result != LIBBSC_NO_ERROR {
        anyhow::bail!("BSC initialization failed with code {}", result);
    }
    Ok(())
}

/// Max OpenMP threads for BSC-internal multithreading.
///
/// Used only by the `*_mt` entry points below. Those entry points MUST NOT be
/// called concurrently from rayon worker threads — they enable OpenMP per-call
/// and call `omp_set_num_threads(BSC_MT_THREADS)`, which mutates a
/// process-wide OpenMP team size. The combined effective parallelism is
/// `(callers in flight) * BSC_MT_THREADS`; oversubscription past the physical
/// CPU count is wasteful and can cause memory blowups inside libsais.
const BSC_MT_THREADS: c_int = 12;

/// Compress data using BSC (matching official SPRING settings)
/// Official SPRING uses -p flag which disables LZP preprocessing
pub fn compress(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_params(data, 0, 0, LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_STATIC)
}

/// Compress data with custom parameters
pub fn compress_with_params(
    data: &[u8],
    lzp_hash_size: i32,
    lzp_min_len: i32,
    block_sorter: i32,
    coder: i32,
) -> Result<Vec<u8>> {
    ensure_initialized()?;

    if data.is_empty() {
        return Ok(vec![]);
    }
    if data.len() > c_int::MAX as usize {
        anyhow::bail!(
            "BSC input too large: {} bytes (max {})",
            data.len(),
            c_int::MAX
        );
    }

    // Allocate output buffer (worst case: input size + header)
    let mut output = vec![0u8; data.len() + LIBBSC_HEADER_SIZE + 1024];

    let compressed_size = unsafe {
        bsc_compress(
            data.as_ptr(),
            output.as_mut_ptr(),
            data.len() as c_int,
            lzp_hash_size as c_int,
            lzp_min_len as c_int,
            block_sorter,
            coder,
            LIBBSC_FEATURE_FASTMODE | cuda_feature(),
        )
    };

    if compressed_size < 0 {
        anyhow::bail!(
            "BSC compression failed with error code: {}",
            compressed_size
        );
    }

    output.truncate(compressed_size as usize);
    output.shrink_to_fit();
    Ok(output)
}

/// Compress using adaptive coder (better compression, slightly slower)
pub fn compress_adaptive(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_params(
        data,
        16,
        128,
        LIBBSC_BLOCKSORTER_BWT,
        LIBBSC_CODER_QLFC_ADAPTIVE,
    )
}

/// Compress using adaptive coder without LZP preprocessing.
/// Faster than full adaptive (skips LZP hash table), same QLFC quality.
/// Best for DNA/quality data where LZP doesn't help much (BWT captures patterns).
pub fn compress_adaptive_no_lzp(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_params(
        data,
        0,
        0,
        LIBBSC_BLOCKSORTER_BWT,
        LIBBSC_CODER_QLFC_ADAPTIVE,
    )
}

/// Block-parallel adaptive compression without LZP.
pub fn compress_parallel_adaptive_no_lzp(data: &[u8]) -> Result<Vec<u8>> {
    compress_parallel_with(data, compress_adaptive_no_lzp)
}

/// Compress a flat byte stream into frame-ready record-aligned blocks.
///
/// `record_byte_offsets[i]` is the byte offset where record `i` begins in
/// `data`. The caller MUST append a trailing sentinel `data.len()` entry,
/// so the slice's length is `num_records + 1`. The sentinel makes the
/// "find largest offset ≤ target" search uniform and prevents an off-by-one
/// where the final record always lands in its own extra block.
///
/// The encoder walks targets at `target_block_bytes` strides and, for each
/// target boundary, finds the largest record offset ≤ that target. The block
/// runs `[prev_boundary_record, this_boundary_record)`.
///
/// **Forward-progress rule.** If no record offset > `prev_boundary` is also
/// ≤ the target boundary (i.e., a single record exceeds `target_block_bytes`),
/// the encoder takes the *next* record boundary — yielding an oversize block
/// containing exactly one record — never an empty interior block.
///
/// Returns one tuple per block: `(record_count, compressed_block_bytes)`,
/// where the compressed bytes are the output of `compress_adaptive_no_lzp`
/// on that block's slice of `data` (no per-block frame prefix — callers add the
/// framing).
/// `static_coder` selects libbsc's QLFC entropy coder per block: `true` uses
/// the static coder ([`compress`], no LZP) — faster to encode AND decode at a
/// small ratio cost; `false` uses the adaptive coder ([`compress_adaptive_no_lzp`])
/// — best ratio. The coder is recorded in each block header, so the decoder
/// auto-detects it; archives from either setting decode identically.
pub fn compress_parallel_with_breakpoints(
    data: &[u8],
    record_byte_offsets: &[usize],
    target_block_bytes: usize,
    static_coder: bool,
) -> Result<Vec<(u32, Vec<u8>)>> {
    assert!(target_block_bytes > 0, "target_block_bytes must be > 0");
    if data.is_empty() || record_byte_offsets.len() < 2 {
        return Ok(Vec::new());
    }
    assert_eq!(
        *record_byte_offsets.last().unwrap(),
        data.len(),
        "record_byte_offsets must end with sentinel == data.len()"
    );

    let num_records = record_byte_offsets.len() - 1; // exclude sentinel
    let mut block_record_ends: Vec<usize> = vec![0];
    let mut cur_record = 0usize;

    while cur_record < num_records {
        let cur_byte = record_byte_offsets[cur_record];
        let target_byte = cur_byte + target_block_bytes;
        // Number of offsets in [cur_record..=num_records] that are ≤ target_byte.
        let partition_idx = record_byte_offsets[cur_record..=num_records]
            .partition_point(|&off| off <= target_byte);
        let mut next_record = cur_record + partition_idx.saturating_sub(1);
        if next_record <= cur_record {
            // Forward-progress: single record > target. Emit a one-record
            // oversize block instead of an empty interior one.
            next_record = cur_record + 1;
        }
        next_record = next_record.min(num_records);
        block_record_ends.push(next_record);
        cur_record = next_record;
    }

    use rayon::prelude::*;
    let plans: Vec<(usize, usize, u32)> = block_record_ends
        .windows(2)
        .map(|w| {
            let start_byte = record_byte_offsets[w[0]];
            let end_byte = record_byte_offsets[w[1]];
            (start_byte, end_byte, (w[1] - w[0]) as u32)
        })
        .collect();

    // Ultra emits a few large blocks (target ≥128 MiB), so rayon inter-block
    // parallelism is starved and a single-threaded forward-BWT on each big block
    // dominates. Use libbsc's OpenMP MT path for those. Each `*_mt` call spawns
    // ~BSC_MT_THREADS OpenMP threads and runs a memory-heavy libsais forward-BWT,
    // so we BOUND how many run concurrently to ≈cores/BSC_MT_THREADS: total
    // OpenMP threads stay near the core count and concurrent libsais memory stays
    // bounded — avoiding the oversubscription + memory blowup the `*_mt` docs warn
    // about (an unbounded `par_iter` + MT would do exactly that). Ultra's handful
    // of big blocks still run in parallel; the cap only protects the pathological
    // many-big-blocks case. (The static `--fast` coder keeps the rayon path —
    // there is no static MT entry point.)
    //
    // CORRECTNESS: libbsc's OpenMP MT-BWT intermittently produces a CORRUPT,
    // non-decodable stream on SMALL, highly repetitive (degenerate) blocks — a
    // parallel suffix-sort race on near-equal suffixes when the per-thread
    // partitions are tiny (reproduced at <=640 KiB; >=1 MiB and large blocks are
    // unaffected; single-thread BSC is always correct). The MT regime is selected
    // by the *target* block size, but the *actual* block can be far smaller (e.g.
    // a small input whose whole stream is one tiny block). So we gate MT on the
    // ACTUAL block length and route blocks below `MIN_MT_BLOCK_BYTES` through the
    // reliable single-thread adaptive coder (which is fast for small blocks
    // anyway). See finding_ultra_mtbsc_degenerate_corruption.
    const BIG_BLOCK_MT_THRESHOLD: usize = 128 * 1024 * 1024;
    // Actual-block floor for the MT path. 16 MiB is ~25x the largest observed
    // failure (~640 KiB) and well below real full-chunk ultra blocks, so MT speed
    // is preserved where it matters and tiny blocks never reach the buggy path.
    const MIN_MT_BLOCK_BYTES: usize = 16 * 1024 * 1024;
    if !static_coder && target_block_bytes >= BIG_BLOCK_MT_THRESHOLD {
        let cores = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(8);
        let max_concurrent = (cores / BSC_MT_THREADS as usize).max(1);
        let mut out: Vec<(u32, Vec<u8>)> = Vec::with_capacity(plans.len());
        for group in plans.chunks(max_concurrent) {
            // Compress up to `max_concurrent` big blocks at once (each MT-parallel
            // internally), then move to the next group. `chunks` + ordered collect
            // preserve block order. Blocks below MIN_MT_BLOCK_BYTES use the
            // single-thread adaptive coder (correct on degenerate input).
            let group_out: Vec<(u32, Vec<u8>)> = group
                .par_iter()
                .map(|(start, end, record_count)| {
                    let block = &data[*start..*end];
                    let compressed = if block.len() >= MIN_MT_BLOCK_BYTES {
                        compress_adaptive_mt_no_lzp(block)?
                    } else {
                        compress_adaptive_no_lzp(block)?
                    };
                    Ok::<_, anyhow::Error>((*record_count, compressed))
                })
                .collect::<Result<Vec<_>>>()?;
            out.extend(group_out);
        }
        return Ok(out);
    }

    let block_compress: fn(&[u8]) -> Result<Vec<u8>> = if static_coder {
        compress
    } else {
        compress_adaptive_no_lzp
    };

    plans
        .par_iter()
        .map(|(start, end, record_count)| {
            let block = &data[*start..*end];
            let compressed = block_compress(block)?;
            Ok((*record_count, compressed))
        })
        .collect()
}

/// Compress using adaptive coder with BSC-internal multithreading, no LZP.
/// Best for raw DNA streams where LZP adds overhead without helping.
pub fn compress_adaptive_mt_no_lzp(data: &[u8]) -> Result<Vec<u8>> {
    compress_mt_inner(data, 0, 0)
}

fn compress_mt_inner(data: &[u8], lzp_hash_size: i32, lzp_min_len: i32) -> Result<Vec<u8>> {
    ensure_initialized()?;

    if data.is_empty() {
        return Ok(vec![]);
    }
    if data.len() > c_int::MAX as usize {
        anyhow::bail!(
            "BSC input too large: {} bytes (max {})",
            data.len(),
            c_int::MAX
        );
    }

    let mut output = vec![0u8; data.len() + LIBBSC_HEADER_SIZE + 1024];

    let compressed_size = unsafe {
        // Limit OpenMP threads per BSC call (only affects MT path, not default rayon path).
        // Side effect: mutates process-wide OMP team size — see BSC_MT_THREADS doc.
        omp_set_num_threads(BSC_MT_THREADS);
        bsc_compress(
            data.as_ptr(),
            output.as_mut_ptr(),
            data.len() as c_int,
            lzp_hash_size as c_int,
            lzp_min_len as c_int,
            LIBBSC_BLOCKSORTER_BWT,
            LIBBSC_CODER_QLFC_ADAPTIVE,
            LIBBSC_FEATURE_FASTMODE | LIBBSC_FEATURE_MULTITHREADING | cuda_feature(),
        )
    };

    if compressed_size < 0 {
        anyhow::bail!(
            "BSC compression failed with error code: {}",
            compressed_size
        );
    }

    output.truncate(compressed_size as usize);
    output.shrink_to_fit();
    Ok(output)
}


/// Internal helper: compress data by splitting into blocks and compressing each in parallel.
/// `block_compressor` is called on each block.
///
/// Wire layout (per-block frame):
///   `[num_blocks: u32 LE]`
///   then, per block: `[block_len: u32 LE][record_count: u32 LE][crc32: u32 LE][compressed BSC payload]`
///
/// `record_count` is set to 0 for this entry point — archives produced
/// here are decoded via mmap or buffered-streaming routes that do not consume
/// `record_count` semantically. The frame CRC covers `record_count || payload`
/// (see `write_block_frame_header`), which is consistent with the decoder.
fn compress_parallel_with<F>(data: &[u8], block_compressor: F) -> Result<Vec<u8>>
where
    F: Fn(&[u8]) -> Result<Vec<u8>> + Sync,
{
    use rayon::prelude::*;

    ensure_initialized()?;

    if data.is_empty() {
        return Ok(vec![]);
    }

    // For small data, use single-block compression (no overhead)
    if data.len() <= BSC_BLOCK_SIZE {
        let compressed = block_compressor(data)?;
        // record_count=0: this entry point's archives are decoded via
        // mmap/buffered paths that don't consume record_count.
        let mut output = Vec::with_capacity(4 + 12 + compressed.len());
        output.extend_from_slice(&1u32.to_le_bytes()); // num_blocks = 1
        write_block_frame_header(&mut output, compressed.len() as u32, 0, &compressed);
        output.extend_from_slice(&compressed);
        return Ok(output);
    }

    // Split into blocks and compress in parallel
    let blocks: Vec<&[u8]> = data.chunks(BSC_BLOCK_SIZE).collect();
    let num_blocks = blocks.len();

    let compressed_blocks: Vec<Vec<u8>> = blocks
        .par_iter()
        .map(|block| block_compressor(block))
        .collect::<Result<Vec<_>>>()?;

    let mut output = Vec::new();
    output.extend_from_slice(&(num_blocks as u32).to_le_bytes());

    for block in compressed_blocks {
        // record_count=0: placeholder for legacy in-memory path (see doc above)
        write_block_frame_header(&mut output, block.len() as u32, 0, &block);
        output.extend_from_slice(&block);
    }

    Ok(output)
}

/// Compress data by splitting into blocks and compressing each in parallel.
///
/// Each block is compressed independently with BSC using rayon's thread pool.
/// Format: [num_blocks: u32][block_compressed_len: u32, block_data]...
/// Falls back to single-block for data smaller than BSC_BLOCK_SIZE.
pub fn compress_parallel(data: &[u8]) -> Result<Vec<u8>> {
    compress_parallel_with(data, compress)
}

/// Compress data by splitting into blocks and compressing each in parallel
/// using the adaptive coder (LZP + adaptive QLFC for better compression).
///
/// Same multi-block format as compress_parallel — decompression is identical.
pub fn compress_parallel_adaptive(data: &[u8]) -> Result<Vec<u8>> {
    compress_parallel_with(data, compress_adaptive)
}

/// Decompress exactly one BSC block — the payload of a single block frame.
///
/// Unlike `decompress_parallel`, this function takes no inner batch loop
/// and allocates no extra blocks beyond what libbsc needs for the single
/// block. It is the building block of the bounded streaming decompressor:
/// the per-stream decompressor thread calls `decompress_one_block` once
/// per block from its `StreamIndex`, then sends the result through a
/// bounded channel to the writer.
///
/// Callers must NOT pass a multi-block stream here — the input must be a
/// single block's payload bytes (no `[num_blocks][block_len]...` framing).
/// Use `decompress_parallel` for buffered multi-block decompression.
pub fn decompress_one_block(data: &[u8], max_block_len: usize) -> Result<Vec<u8>> {
    // Re-uses the existing single-block `decompress_max` helper; the wrapper
    // exists to give the bounded streaming path a stable, intent-revealing
    // entry point separate from `decompress_parallel`'s batched flavor.
    // `max_block_len` is threaded through so callers (e.g. the ultra path)
    // can accept blocks larger than the default 64 MiB cap.
    decompress_max(data, max_block_len)
}

/// Decompress multi-block BSC-compressed data (from `compress_parallel*`).
///
/// Blocks are decompressed in parallel using rayon, then concatenated in order.
/// Wire layout (per-block frame):
///   `[num_blocks: u32 LE]`
///   then, per block: `[block_len: u32 LE][record_count: u32 LE][crc32: u32 LE][BSC payload: block_len bytes]`
///
/// The CRC32 covers `record_count || payload` and is verified **before**
/// invoking libbsc, so disk corruption is reported with a clear "CRC32
/// mismatch" error rather than bubbling through libbsc as a cryptic
/// `block_info failed` or producing silently wrong output.
pub fn decompress_parallel(data: &[u8]) -> Result<Vec<u8>> {
    // Unbounded by default (backward-compatible). Callers decoding UNTRUSTED input
    // should prefer `decompress_parallel_max` with a known output bound.
    decompress_parallel_max(data, usize::MAX)
}

/// Like [`decompress_parallel`], but bails if the cumulative decompressed output
/// would exceed `max_total` bytes. This bounds memory on untrusted input where a
/// BSC "bomb" (tiny compressed payload → huge decompressed output, stacked across
/// many blocks) could otherwise force an unbounded allocation and abort the process.
/// Each block is additionally capped at `min(BSC_MAX_BLOCK_LEN, max_total)`, so no
/// single block can exceed the bound either. The transient peak is the legitimate
/// `max_total` plus at most `threads × per-block cap` of in-flight decode.
pub fn decompress_parallel_max(data: &[u8], max_total: usize) -> Result<Vec<u8>> {
    use rayon::prelude::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    if data.is_empty() {
        return Ok(vec![]);
    }

    if data.len() < 4 {
        anyhow::bail!("BSC parallel: data too small for header");
    }

    let num_blocks = super::read_le_u32(data, 0)? as usize;

    // Each block frame contributes at least 12 bytes of prefix (4 len + 4 record_count + 4 crc)
    // + 1 byte of payload. Reject obviously malformed inputs before we Vec::with_capacity.
    if num_blocks > 0 && num_blocks.saturating_mul(13) > data.len().saturating_sub(4) {
        anyhow::bail!(
            "BSC parallel: num_blocks={} exceeds remaining payload ({} bytes)",
            num_blocks,
            data.len() - 4,
        );
    }

    // First pass: collect (block_slice, record_count, expected_crc) tuples
    let mut offset = 4;
    let mut blocks: Vec<(&[u8], u32, u32)> = Vec::with_capacity(num_blocks);
    for idx in 0..num_blocks {
        if offset + 12 > data.len() {
            anyhow::bail!("BSC parallel: truncated block header for block {idx}");
        }
        let block_len = super::read_le_u32(data, offset)? as usize;
        offset += 4;
        let record_count = super::read_le_u32(data, offset)?;
        offset += 4;
        let expected_crc = super::read_le_u32(data, offset)?;
        offset += 4;

        if offset + block_len > data.len() {
            anyhow::bail!("BSC parallel: truncated block data for block {idx}");
        }
        blocks.push((
            &data[offset..offset + block_len],
            record_count,
            expected_crc,
        ));
        offset += block_len;
    }

    // Second pass: verify CRC + decompress each block in parallel, short-circuit on first
    // error. Each block is capped at `min(BSC_MAX_BLOCK_LEN, max_total)` and a shared
    // running total bails as soon as the cumulative decompressed output exceeds
    // `max_total` — bounding a BSC bomb on untrusted input.
    let per_block_cap = super::decompress_impl::BSC_MAX_BLOCK_LEN.min(max_total.max(1));
    let running = AtomicUsize::new(0);
    let decompressed_blocks: Vec<Vec<u8>> = blocks
        .par_iter()
        .enumerate()
        .map(|(idx, (block, record_count, expected_crc))| {
            verify_block_frame_crc(*expected_crc, *record_count, block)
                .with_context(|| format!("BSC block {idx}"))?;
            let out = decompress_max(block, per_block_cap)?;
            let total_after = running
                .fetch_add(out.len(), Ordering::Relaxed)
                .saturating_add(out.len());
            if total_after > max_total {
                anyhow::bail!(
                    "BSC parallel: decompressed output exceeds cap of {max_total} bytes (bomb?)"
                );
            }
            Ok(out)
        })
        .collect::<Result<Vec<_>>>()?;

    let total_size: usize = decompressed_blocks.iter().map(|b| b.len()).sum();
    let mut output = Vec::with_capacity(total_size);
    for block in decompressed_blocks {
        output.extend_from_slice(&block);
    }

    Ok(output)
}

/// Parsed per-block frame header. Fields match the on-disk layout
/// `[block_len: u32 LE][record_count: u32 LE][crc32: u32 LE]`.
pub struct BlockFrameHeader {
    pub block_len: u32,
    pub record_count: u32,
    /// CRC32 over `record_count || payload`, as written by `write_block_frame_header`.
    pub expected_crc: u32,
}

/// Read a per-block frame header from `reader`. Does NOT read the payload —
/// caller pulls `block_len` bytes after this returns.
pub fn read_block_frame_header<R: std::io::Read>(reader: &mut R) -> Result<BlockFrameHeader> {
    let mut buf = [0u8; 12];
    reader.read_exact(&mut buf)?;
    Ok(BlockFrameHeader {
        block_len: u32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]),
        record_count: u32::from_le_bytes([buf[4], buf[5], buf[6], buf[7]]),
        expected_crc: u32::from_le_bytes([buf[8], buf[9], buf[10], buf[11]]),
    })
}

/// Compute the per-block-frame CRC32 over `record_count || payload`.
/// The single source of truth for the frame CRC scope; if the primitive
/// changes, both reader and writer see the change at once.
pub(crate) fn compute_block_frame_crc(record_count: u32, payload: &[u8]) -> u32 {
    let mut crc = Crc::new();
    crc.update(&record_count.to_le_bytes());
    crc.update(payload);
    crc.sum()
}

/// Write a per-block frame header to `out`. The CRC32 is computed over
/// `record_count.to_le_bytes() || payload` — `record_count` is included
/// in the CRC scope so that record-count-only corruption (which would
/// misalign the bounded streaming cursor mid-stream) is caught at the
/// block layer with a clear diagnostic rather than silently producing
/// wrong output. This must match what `verify_block_frame_crc` reads on
/// decode.
pub fn write_block_frame_header(out: &mut Vec<u8>, block_len: u32, record_count: u32, payload: &[u8]) {
    let checksum = compute_block_frame_crc(record_count, payload);
    out.extend_from_slice(&block_len.to_le_bytes());
    out.extend_from_slice(&record_count.to_le_bytes());
    out.extend_from_slice(&checksum.to_le_bytes());
}

/// Verify a block frame's CRC32. Returns Err with a descriptive message on
/// mismatch so callers can bubble it up with stream/block context.
pub fn verify_block_frame_crc(expected: u32, record_count: u32, payload: &[u8]) -> Result<()> {
    let actual = compute_block_frame_crc(record_count, payload);
    if actual != expected {
        anyhow::bail!(
            "block frame CRC32 mismatch: stored {:08x}, computed {:08x} — archive is corrupted",
            expected,
            actual
        );
    }
    Ok(())
}

/// Decompress BSC-compressed data
pub fn decompress(data: &[u8]) -> Result<Vec<u8>> {
    decompress_with_features(
        data,
        LIBBSC_FEATURE_FASTMODE | cuda_feature(),
        super::decompress_impl::BSC_MAX_BLOCK_LEN,
    )
}

/// Like [`decompress`], but permits blocks up to `max_block_len` decompressed
/// bytes. Used by the `--ultra` path, whose blocks legitimately exceed the
/// default [`BSC_MAX_BLOCK_LEN`] cap (see [`BSC_MAX_BLOCK_LEN_ULTRA`]).
///
/// [`BSC_MAX_BLOCK_LEN`]: super::decompress_impl::BSC_MAX_BLOCK_LEN
/// [`BSC_MAX_BLOCK_LEN_ULTRA`]: super::decompress_impl::BSC_MAX_BLOCK_LEN_ULTRA
pub fn decompress_max(data: &[u8], max_block_len: usize) -> Result<Vec<u8>> {
    decompress_with_features(
        data,
        LIBBSC_FEATURE_FASTMODE | cuda_feature(),
        max_block_len,
    )
}

/// Decompress with BSC-internal multithreading (OpenMP parallel inverse BWT),
/// permitting blocks up to `max_block_len` decompressed bytes (for the `--ultra`
/// path's large blocks).
pub fn decompress_mt_max(data: &[u8], max_block_len: usize) -> Result<Vec<u8>> {
    unsafe {
        omp_set_num_threads(BSC_MT_THREADS);
    }
    decompress_with_features(
        data,
        LIBBSC_FEATURE_FASTMODE | LIBBSC_FEATURE_MULTITHREADING | cuda_feature(),
        max_block_len,
    )
}

fn decompress_with_features(data: &[u8], features: c_int, max_block_len: usize) -> Result<Vec<u8>> {
    ensure_initialized()?;

    if data.is_empty() {
        return Ok(vec![]);
    }

    if data.len() < LIBBSC_HEADER_SIZE {
        anyhow::bail!(
            "BSC: compressed data too small (need at least {} bytes for header)",
            LIBBSC_HEADER_SIZE
        );
    }
    if data.len() > c_int::MAX as usize {
        anyhow::bail!(
            "BSC compressed block too large: {} bytes (max {})",
            data.len(),
            c_int::MAX
        );
    }

    let mut block_size: c_int = 0;
    let mut data_size: c_int = 0;

    let info_result = unsafe {
        bsc_block_info(
            data.as_ptr(),
            data.len() as c_int,
            &mut block_size,
            &mut data_size,
            features,
        )
    };

    if info_result != LIBBSC_NO_ERROR {
        anyhow::bail!("BSC block_info failed with error code: {}", info_result);
    }

    if data_size <= 0 {
        anyhow::bail!("BSC: invalid decompressed size: {}", data_size);
    }

    // Bound the allocation against the caller's per-block maximum so a corrupt
    // block header reporting a huge `data_size` can't drive a multi-GB zero-fill
    // (× rayon workers) — an OOM-amplification DoS. The default/bounded path
    // passes `BSC_MAX_BLOCK_LEN` (25 MB blocks); the `--ultra` path passes the
    // larger `BSC_MAX_BLOCK_LEN_ULTRA` for its big-block `UltraBigBlock` archives.
    if data_size as usize > max_block_len {
        anyhow::bail!(
            "BSC: decompressed block size {} exceeds maximum {} (corrupt header)",
            data_size,
            max_block_len
        );
    }

    let mut output = vec![0u8; data_size as usize];

    let decompress_result = unsafe {
        bsc_decompress(
            data.as_ptr(),
            data.len() as c_int,
            output.as_mut_ptr(),
            data_size,
            features,
        )
    };

    if decompress_result != LIBBSC_NO_ERROR {
        anyhow::bail!(
            "BSC decompression failed with error code: {}",
            decompress_result
        );
    }

    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip_simple() {
        let data = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let compressed = compress(data).unwrap();
        let decompressed = decompress(&compressed).unwrap();
        assert_eq!(data.as_slice(), decompressed.as_slice());
    }

    #[test]
    fn test_decompress_parallel_max_caps_output() {
        // A highly compressible stream: tiny compressed, larger decompressed — the
        // shape an untrusted "bomb" exploits. The cap must reject when it is below the
        // real output size, and roundtrip cleanly when it is at/above.
        let data = vec![b'A'; 200_000];
        let compressed = compress_parallel(&data).unwrap();
        assert!(
            compressed.len() < data.len(),
            "expected real compression for this test"
        );

        // Cap >= output: full roundtrip.
        let ok = decompress_parallel_max(&compressed, data.len()).unwrap();
        assert_eq!(ok, data, "cap >= output must decode losslessly");

        // Cap below the real output: must bail (no unbounded allocation), not panic.
        // Either the per-block cap (`decompress_max`) or the cumulative cap fires,
        // depending on block layout; both bound the output. The property under test is
        // that it returns Err cleanly.
        let err = decompress_parallel_max(&compressed, data.len() - 1).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("exceeds cap") || msg.contains("exceeds maximum"),
            "expected a cap rejection, got: {err}"
        );

        // The unbounded wrapper still decodes (backward-compatible).
        assert_eq!(decompress_parallel(&compressed).unwrap(), data);
    }

    #[test]
    fn test_roundtrip_genomic() {
        let data = b"ACGTACGTNNNACGTACGTACGTACGTNNACGTACGTACGTACGTACGT\
                      ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let compressed = compress(data).unwrap();
        let decompressed = decompress(&compressed).unwrap();
        assert_eq!(data.as_slice(), decompressed.as_slice());
    }

    #[test]
    fn test_roundtrip_quality() {
        let data = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\
                      HHHHHHHHHHGGGGGGGGFFFFFFFFFFEEEEEEEEDDDDDDDDDCCCCCCC";
        let compressed = compress(data).unwrap();
        let decompressed = decompress(&compressed).unwrap();
        assert_eq!(data.as_slice(), decompressed.as_slice());
    }

    #[test]
    fn test_compression_ratio() {
        // Repetitive data should compress well
        let data = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let compressed = compress(data).unwrap();
        println!(
            "Original: {} bytes, Compressed: {} bytes, Ratio: {:.2}x",
            data.len(),
            compressed.len(),
            data.len() as f64 / compressed.len() as f64
        );
        assert!(compressed.len() < data.len());
    }

    #[test]
    fn test_empty_data() {
        let data = b"";
        let compressed = compress(data).unwrap();
        assert_eq!(compressed.len(), 0);
        let decompressed = decompress(&compressed).unwrap();
        assert_eq!(decompressed.len(), 0);
    }

    #[test]
    fn test_adaptive_coder() {
        let data = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let compressed_static = compress(data).unwrap();
        let compressed_adaptive = compress_adaptive(data).unwrap();

        // Both should decompress correctly
        let decompressed = decompress(&compressed_adaptive).unwrap();
        assert_eq!(data.as_slice(), decompressed.as_slice());

        println!(
            "Static coder: {} bytes, Adaptive coder: {} bytes",
            compressed_static.len(),
            compressed_adaptive.len()
        );
    }

    #[test]
    fn test_block_framing_roundtrip() {
        // Build a synthetic single-block stream and round-trip it through
        // the frame helpers. The CRC must cover `record_count || payload`, so
        // a flipped bit in record_count is caught at the block layer.
        let payload = b"hello world".to_vec();
        let record_count: u32 = 7;
        let mut buf = Vec::new();
        super::write_block_frame_header(&mut buf, payload.len() as u32, record_count, &payload);
        buf.extend_from_slice(&payload);

        // Re-read using the parser; assert exact fields.
        let mut cursor = std::io::Cursor::new(&buf);
        let hdr = super::read_block_frame_header(&mut cursor).expect("read frame header");
        assert_eq!(hdr.block_len, payload.len() as u32);
        assert_eq!(hdr.record_count, record_count);

        let mut read_payload = vec![0u8; hdr.block_len as usize];
        use std::io::Read;
        cursor.read_exact(&mut read_payload).expect("read payload");
        assert_eq!(read_payload, payload);

        // CRC verification on the same (record_count, payload) tuple returns Ok.
        super::verify_block_frame_crc(hdr.expected_crc, record_count, &payload).expect("crc ok");

        // Flip a byte in record_count and verify the CRC fails.
        let mut tampered = buf.clone();
        // Buffer layout: [block_len:4][record_count:4][crc:4][payload].
        // record_count byte 0 is at offset 4.
        // Sanity-check the layout comment: block_len occupies bytes 0..4,
        // so byte 4 is the low byte of record_count.
        assert_eq!(tampered[4], (record_count & 0xff) as u8);
        tampered[4] ^= 0x01;
        let mut cursor2 = std::io::Cursor::new(&tampered);
        let hdr2 = super::read_block_frame_header(&mut cursor2).unwrap();
        let mut p2 = vec![0u8; hdr2.block_len as usize];
        cursor2.read_exact(&mut p2).unwrap();
        let err =
            super::verify_block_frame_crc(hdr2.expected_crc, hdr2.record_count, &p2).unwrap_err();
        assert!(format!("{err}").contains("CRC32 mismatch"));
    }

    #[test]
    fn test_decompress_one_block_no_batching() {
        // Compress a small input through compress_adaptive_no_lzp (single-block path),
        // then assert decompress_one_block returns the original.
        let original = b"ACGTACGTACGTACGTACGTACGTACGTACGT".repeat(1000);
        let compressed = super::compress_adaptive_no_lzp(&original).expect("compress");
        let decompressed = super::decompress_one_block(
            &compressed,
            super::super::decompress_impl::BSC_MAX_BLOCK_LEN,
        )
        .expect("decompress_one_block");
        assert_eq!(decompressed, original);
    }

    #[test]
    fn test_compress_with_breakpoints_record_aligned_blocks() {
        // Build a stream of 1000 varint-prefixed records, each ~100 bytes.
        let mut data: Vec<u8> = Vec::new();
        let mut record_offsets: Vec<usize> = Vec::new();
        for i in 0..1000 {
            record_offsets.push(data.len());
            let payload = format!("record_{i:04}_").repeat(7); // ~100 bytes
            data.push(payload.len() as u8); // varint length (single byte, len < 128)
            data.extend_from_slice(payload.as_bytes());
        }
        record_offsets.push(data.len()); // sentinel
        let target = data.len() / 10; // force ~10 blocks
        let blocks =
            super::compress_parallel_with_breakpoints(&data, &record_offsets, target, false)
                .expect("compress_with_breakpoints");

        // Sum of record_counts must equal total records (1000).
        let total_records: u32 = blocks.iter().map(|(rc, _)| rc).sum();
        assert_eq!(total_records, 1000, "sum of record_counts mismatches");

        // Round-trip: decompress each block and concatenate; bytes match the input.
        let mut reconstructed = Vec::new();
        for (_, compressed) in &blocks {
            reconstructed.extend_from_slice(
                &super::decompress_one_block(
                    compressed,
                    super::super::decompress_impl::BSC_MAX_BLOCK_LEN,
                )
                .unwrap(),
            );
        }
        assert_eq!(
            reconstructed, data,
            "compressed/decompressed round-trip differs"
        );

        // No block (interior or trailing) should have record_count == 0.
        for (i, (rc, _)) in blocks.iter().enumerate() {
            assert!(
                *rc > 0,
                "block {i} has record_count=0 (interior or trailing)"
            );
        }
    }

    #[test]
    fn test_compress_with_breakpoints_oversized_record_takes_next_boundary() {
        // One record exceeds target_block_bytes. Forward-progress rule says it
        // lands in its own block (size > target), not in an empty interior block.
        let huge_record: Vec<u8> = vec![b'X'; 1_000_000];
        let mut data = Vec::new();
        let mut record_offsets = Vec::new();

        // Record 0: small (50 bytes).
        record_offsets.push(data.len());
        data.push(50);
        data.extend_from_slice(&[b'a'; 50]);

        // Record 1: 1 MB.
        record_offsets.push(data.len());
        let mut v = 1_000_000usize;
        while v >= 0x80 {
            data.push(((v & 0x7F) | 0x80) as u8);
            v >>= 7;
        }
        data.push(v as u8);
        data.extend_from_slice(&huge_record);

        // Record 2: small (50 bytes).
        record_offsets.push(data.len());
        data.push(50);
        data.extend_from_slice(&[b'b'; 50]);
        record_offsets.push(data.len()); // sentinel

        let blocks =
            super::compress_parallel_with_breakpoints(&data, &record_offsets, 1024, false).unwrap();
        let total_records: u32 = blocks.iter().map(|(rc, _)| rc).sum();
        assert_eq!(
            total_records, 3,
            "expected 3 records total, got {total_records}"
        );
    }

    #[test]
    fn test_compress_with_breakpoints_sentinel_only_returns_empty() {
        // num_records == 0 (sentinel only). Function must return Ok(Vec::new()),
        // not crash, not produce a phantom empty block.
        let blocks = super::compress_parallel_with_breakpoints(&[], &[0], 1024, false).unwrap();
        assert!(
            blocks.is_empty(),
            "sentinel-only input must produce zero blocks"
        );
    }

    #[test]
    fn test_compress_with_breakpoints_static_coder_roundtrips() {
        // The static-QLFC coder (static_coder=true) must produce blocks that
        // decode byte-exact, just like the adaptive coder. The block header
        // records the coder so the decoder auto-detects it.
        let mut data: Vec<u8> = Vec::new();
        let mut record_offsets: Vec<usize> = Vec::new();
        for i in 0..1000 {
            record_offsets.push(data.len());
            let payload = format!("record_{i:04}_").repeat(7);
            data.push(payload.len() as u8);
            data.extend_from_slice(payload.as_bytes());
        }
        record_offsets.push(data.len());
        let target = data.len() / 10;

        let blocks =
            super::compress_parallel_with_breakpoints(&data, &record_offsets, target, true)
                .unwrap();
        let mut reconstructed = Vec::new();
        for (_, compressed) in &blocks {
            reconstructed.extend_from_slice(
                &super::decompress_one_block(
                    compressed,
                    super::super::decompress_impl::BSC_MAX_BLOCK_LEN,
                )
                .unwrap(),
            );
        }
        assert_eq!(reconstructed, data, "static-coder round-trip differs");
    }

    #[test]
    fn decompress_one_block_respects_ultra_cap() {
        // ~70 MiB of compressible data → one block whose decompressed size > 64 MiB.
        let raw = vec![b'A'; 70 * 1024 * 1024];
        let compressed = compress_adaptive_mt_no_lzp(&raw).unwrap();
        // Default cap (64 MiB + 4 KiB margin) rejects it.
        assert!(
            decompress_one_block(
                &compressed,
                super::super::decompress_impl::BSC_MAX_BLOCK_LEN
            )
            .is_err()
        );
        // Ultra cap (768 MiB) accepts it and round-trips.
        let out = decompress_one_block(
            &compressed,
            super::super::decompress_impl::BSC_MAX_BLOCK_LEN_ULTRA,
        )
        .unwrap();
        assert_eq!(out, raw);
    }
}
