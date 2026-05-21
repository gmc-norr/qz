/// FFI bindings to libbsc (Block-Sorting Compressor)
///
/// libbsc is the compression library used by official SPRING.
/// This module provides safe Rust wrappers around the C API.
///
/// Source: https://github.com/IlyaGrebnov/libbsc

use anyhow::Result;
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
const LIBBSC_BLOCKSORTER_BWT: c_int = 1;  // Use BWT for compatibility
const _LIBBSC_BLOCKSORTER_ST7: c_int = 7;  // ST7 is what SPRING uses
const LIBBSC_CODER_QLFC_STATIC: c_int = 1;
const LIBBSC_CODER_QLFC_ADAPTIVE: c_int = 2;

// Features (bitmask)
const LIBBSC_FEATURE_FASTMODE: c_int = 1;
const LIBBSC_FEATURE_MULTITHREADING: c_int = 2;
#[cfg(feature = "cuda")]
const LIBBSC_FEATURE_CUDA: c_int = 8;

/// Returns LIBBSC_FEATURE_CUDA when cuda feature is enabled, 0 otherwise.
#[inline]
fn cuda_feature() -> c_int {
    #[cfg(feature = "cuda")]
    { LIBBSC_FEATURE_CUDA }
    #[cfg(not(feature = "cuda"))]
    { 0 }
}

/// Query GPU VRAM and warn if it's too small for BWT at the given block size.
/// `max_block_bytes` is the largest single BSC block that will be passed to BWT.
/// libcubwt needs ~20.5× the input length in device memory for forward BWT.
#[cfg(feature = "cuda")]
pub fn check_cuda_vram(max_block_bytes: usize) {
    use tracing::warn;

    unsafe extern "C" {
        fn cudaMemGetInfo(free: *mut usize, total: *mut usize) -> c_int;
    }

    let mut free: usize = 0;
    let mut total: usize = 0;
    let status = unsafe { cudaMemGetInfo(&mut free, &mut total) };
    if status != 0 {
        warn!("Could not query GPU VRAM (cudaMemGetInfo returned {}). GPU acceleration may fail.", status);
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
const BSC_BLOCK_SIZE: usize = 25 * 1024 * 1024;  // 25 MB

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
static BSC_INIT_RESULT: std::sync::atomic::AtomicI32 =
    std::sync::atomic::AtomicI32::new(i32::MIN);

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
        anyhow::bail!("BSC input too large: {} bytes (max {})", data.len(), c_int::MAX);
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
        anyhow::bail!("BSC compression failed with error code: {}", compressed_size);
    }

    output.truncate(compressed_size as usize);
    output.shrink_to_fit();
    Ok(output)
}

/// Compress using adaptive coder (better compression, slightly slower)
pub fn compress_adaptive(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_params(data, 16, 128, LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_ADAPTIVE)
}

/// Compress using adaptive coder without LZP preprocessing.
/// Faster than full adaptive (skips LZP hash table), same QLFC quality.
/// Best for DNA/quality data where LZP doesn't help much (BWT captures patterns).
pub fn compress_adaptive_no_lzp(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_params(data, 0, 0, LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_ADAPTIVE)
}

/// Block-parallel adaptive compression without LZP.
pub fn compress_parallel_adaptive_no_lzp(data: &[u8]) -> Result<Vec<u8>> {
    compress_parallel_with(data, compress_adaptive_no_lzp)
}

/// Compress using adaptive coder with BSC-internal multithreading (OpenMP) and LZP.
/// Best for structured non-DNA streams (order indices, read lengths) where LZP helps.
pub fn compress_adaptive_mt(data: &[u8]) -> Result<Vec<u8>> {
    compress_mt_inner(data, 16, 128)
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
        anyhow::bail!("BSC input too large: {} bytes (max {})", data.len(), c_int::MAX);
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
        anyhow::bail!("BSC compression failed with error code: {}", compressed_size);
    }

    output.truncate(compressed_size as usize);
    output.shrink_to_fit();
    Ok(output)
}

/// Compute CRC32 (IEEE) of a slice using flate2's bundled CRC. Matches the
/// algorithm bz-lib uses, so qz and bz can share verification tooling.
///
/// Single source of truth for the qz-v3 per-block CRC algorithm; called by
/// every writer/reader that materializes the `[block_len][crc32][payload]`
/// framing (bsc, openzl, columnar headers, fqzcomp, quality_ctx).
#[inline]
pub(crate) fn block_crc32(data: &[u8]) -> u32 {
    let mut crc = Crc::new();
    crc.update(data);
    crc.sum()
}

/// Internal helper: compress data by splitting into blocks and compressing each in parallel.
/// `block_compressor` is called on each block.
///
/// Wire layout (qz archive v3+):
///   `[num_blocks: u32 LE]`
///   then, per block: `[block_len: u32 LE][crc32: u32 LE][compressed BSC payload]`
///
/// The CRC32 is computed over the compressed BSC payload (post-BSC bytes) and
/// verified on read before invoking libbsc. This catches bit-rot and disk
/// corruption with a clear error instead of bubbling through libbsc as a
/// cryptic "block_info failed" or producing silently wrong output.
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
        let crc = block_crc32(&compressed);
        let mut output = Vec::with_capacity(4 + 8 + compressed.len());
        output.extend_from_slice(&1u32.to_le_bytes()); // 1 block
        output.extend_from_slice(&(compressed.len() as u32).to_le_bytes());
        output.extend_from_slice(&crc.to_le_bytes());
        output.extend_from_slice(&compressed);
        return Ok(output);
    }

    // Split into blocks and compress + CRC in parallel
    let blocks: Vec<&[u8]> = data.chunks(BSC_BLOCK_SIZE).collect();
    let num_blocks = blocks.len();

    let compressed_blocks: Vec<(Vec<u8>, u32)> = blocks
        .par_iter()
        .map(|block| {
            let compressed = block_compressor(block)?;
            let crc = block_crc32(&compressed);
            Ok::<_, anyhow::Error>((compressed, crc))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut output = Vec::new();
    output.extend_from_slice(&(num_blocks as u32).to_le_bytes());

    for (block, crc) in compressed_blocks {
        output.extend_from_slice(&(block.len() as u32).to_le_bytes());
        output.extend_from_slice(&crc.to_le_bytes());
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

/// Decompress multi-block BSC-compressed data (from `compress_parallel*`).
///
/// Blocks are decompressed in parallel using rayon, then concatenated in order.
/// Wire layout (qz archive v3+):
///   `[num_blocks: u32 LE]`
///   then, per block: `[block_len: u32 LE][crc32: u32 LE][BSC payload: block_len bytes]`
///
/// The CRC32 over each block's compressed payload is verified **before**
/// invoking libbsc, so disk corruption is reported with a clear "Chunk CRC32
/// mismatch" error rather than bubbling through libbsc as a cryptic
/// `block_info failed` or producing silently wrong output.
pub fn decompress_parallel(data: &[u8]) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    if data.is_empty() {
        return Ok(vec![]);
    }

    if data.len() < 4 {
        anyhow::bail!("BSC parallel: data too small for header");
    }

    let num_blocks = super::read_le_u32(data, 0)? as usize;

    // Each block contributes at least 8 bytes of prefix (4 len + 4 crc) + 1 byte of payload.
    // Reject obviously malformed inputs before we Vec::with_capacity(num_blocks).
    if num_blocks > 0 && num_blocks.saturating_mul(9) > data.len().saturating_sub(4) {
        anyhow::bail!(
            "BSC parallel: num_blocks={} exceeds remaining payload ({} bytes)",
            num_blocks,
            data.len() - 4,
        );
    }

    // First pass: collect (block_slice, expected_crc) pairs (sequential, just pointer math)
    let mut offset = 4;
    let mut blocks: Vec<(&[u8], u32)> = Vec::with_capacity(num_blocks);
    for idx in 0..num_blocks {
        if offset + 8 > data.len() {
            anyhow::bail!("BSC parallel: truncated block header for block {idx}");
        }
        let block_len = super::read_le_u32(data, offset)? as usize;
        offset += 4;
        let expected_crc = super::read_le_u32(data, offset)?;
        offset += 4;

        if offset + block_len > data.len() {
            anyhow::bail!("BSC parallel: truncated block data for block {idx}");
        }
        blocks.push((&data[offset..offset + block_len], expected_crc));
        offset += block_len;
    }

    // Second pass: verify CRC + decompress each block in parallel, short-circuit on first error.
    let decompressed_blocks: Vec<Vec<u8>> = blocks
        .par_iter()
        .enumerate()
        .map(|(idx, (block, expected_crc))| {
            let actual_crc = block_crc32(block);
            if actual_crc != *expected_crc {
                anyhow::bail!(
                    "BSC block {idx} CRC32 mismatch: stored {:08x}, computed {:08x} — archive is corrupted",
                    expected_crc, actual_crc
                );
            }
            decompress(block)
        })
        .collect::<Result<Vec<_>>>()?;

    let total_size: usize = decompressed_blocks.iter().map(|b| b.len()).sum();
    let mut output = Vec::with_capacity(total_size);
    for block in decompressed_blocks {
        output.extend_from_slice(&block);
    }

    Ok(output)
}

/// Verify all per-block CRC32s in a qz-v3 multi-block stream **whose outer
/// framing is `[num_blocks: u32 LE]` then per block
/// `[block_len: u32 LE][crc32: u32 LE][payload: block_len bytes]`** —
/// without invoking the inner codec.
///
/// This framing is shared by qz-v3 BSC, OpenZL, columnar-header, fqzcomp, and
/// quality_ctx multi-block streams, so a single walker covers all of them.
/// It does **not** apply to:
/// - Raw-zstd streams (no outer length/num_blocks header — just a zstd frame).
/// - The outer per-chunk wrapper used by encoding_type 8/9 (ultra/local-reorder),
///   which is `[num_chunks][chunk_len, chunk_data]...` with no per-block CRC.
///
/// Callers (e.g. `verify_fast`) must dispatch by archive compressor/encoding
/// before calling this helper. Returns the number of blocks verified.
/// Used by `qz verify --fast` to catch bit-rot in O(IO + CRC) time.
pub fn verify_parallel_crcs(data: &[u8]) -> Result<usize> {
    if data.is_empty() {
        return Ok(0);
    }
    if data.len() < 4 {
        anyhow::bail!("BSC parallel: data too small for header");
    }
    let num_blocks = super::read_le_u32(data, 0)? as usize;
    if num_blocks > 0 && num_blocks.saturating_mul(9) > data.len().saturating_sub(4) {
        anyhow::bail!(
            "BSC parallel: num_blocks={} exceeds remaining payload ({} bytes)",
            num_blocks,
            data.len() - 4,
        );
    }
    let mut offset = 4;
    for idx in 0..num_blocks {
        if offset + 8 > data.len() {
            anyhow::bail!("BSC parallel: truncated block header for block {idx}");
        }
        let block_len = super::read_le_u32(data, offset)? as usize;
        offset += 4;
        let expected_crc = super::read_le_u32(data, offset)?;
        offset += 4;
        if offset + block_len > data.len() {
            anyhow::bail!("BSC parallel: truncated block data for block {idx}");
        }
        let actual_crc = block_crc32(&data[offset..offset + block_len]);
        if actual_crc != expected_crc {
            anyhow::bail!(
                "BSC block {idx} CRC32 mismatch: stored {:08x}, computed {:08x} — archive is corrupted",
                expected_crc, actual_crc
            );
        }
        offset += block_len;
    }
    Ok(num_blocks)
}

/// Decompress BSC-compressed data
pub fn decompress(data: &[u8]) -> Result<Vec<u8>> {
    decompress_with_features(data, LIBBSC_FEATURE_FASTMODE | cuda_feature())
}

/// Decompress with BSC-internal multithreading (OpenMP parallel inverse BWT).
pub fn decompress_mt(data: &[u8]) -> Result<Vec<u8>> {
    unsafe { omp_set_num_threads(BSC_MT_THREADS); }
    decompress_with_features(data, LIBBSC_FEATURE_FASTMODE | LIBBSC_FEATURE_MULTITHREADING | cuda_feature())
}

fn decompress_with_features(data: &[u8], features: c_int) -> Result<Vec<u8>> {
    ensure_initialized()?;

    if data.is_empty() {
        return Ok(vec![]);
    }

    if data.len() < LIBBSC_HEADER_SIZE {
        anyhow::bail!("BSC: compressed data too small (need at least {} bytes for header)", LIBBSC_HEADER_SIZE);
    }
    if data.len() > c_int::MAX as usize {
        anyhow::bail!("BSC compressed block too large: {} bytes (max {})", data.len(), c_int::MAX);
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
        anyhow::bail!("BSC decompression failed with error code: {}", decompress_result);
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
        println!("Original: {} bytes, Compressed: {} bytes, Ratio: {:.2}x",
                 data.len(), compressed.len(), data.len() as f64 / compressed.len() as f64);
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

        println!("Static coder: {} bytes, Adaptive coder: {} bytes",
                 compressed_static.len(), compressed_adaptive.len());
    }
}
