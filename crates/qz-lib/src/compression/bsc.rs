/// FFI bindings to libbsc (Block-Sorting Compressor)
///
/// libbsc is the compression library used by official SPRING.
/// This module provides safe Rust wrappers around the C API.
///
/// Source: https://github.com/IlyaGrebnov/libbsc

use anyhow::Result;
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

/// Ensure bsc_init is called exactly once
static BSC_INIT: Once = Once::new();

fn ensure_initialized() {
    BSC_INIT.call_once(|| {
        unsafe {
            // Init with both features so we can use either per-call
            let features = LIBBSC_FEATURE_FASTMODE | LIBBSC_FEATURE_MULTITHREADING | cuda_feature();
            let result = bsc_init(features);
            if result != LIBBSC_NO_ERROR {
                eprintln!("Warning: BSC initialization failed with code {}", result);
            }
        }
    });
}

/// Max OpenMP threads for BSC-internal multithreading (compress_adaptive_mt only).
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
    ensure_initialized();

    if data.is_empty() {
        return Ok(vec![]);
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
    ensure_initialized();

    if data.is_empty() {
        return Ok(vec![]);
    }

    let mut output = vec![0u8; data.len() + LIBBSC_HEADER_SIZE + 1024];

    let compressed_size = unsafe {
        // Limit OpenMP threads per BSC call (only affects MT path, not default rayon path)
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
    Ok(output)
}

/// Internal helper: compress data by splitting into blocks and compressing each in parallel.
/// `block_compressor` is called on each block.
fn compress_parallel_with<F>(data: &[u8], block_compressor: F) -> Result<Vec<u8>>
where
    F: Fn(&[u8]) -> Result<Vec<u8>> + Sync,
{
    use rayon::prelude::*;

    ensure_initialized();

    if data.is_empty() {
        return Ok(vec![]);
    }

    // For small data, use single-block compression (no overhead)
    if data.len() <= BSC_BLOCK_SIZE {
        let compressed = block_compressor(data)?;
        let mut output = Vec::with_capacity(4 + 4 + compressed.len());
        output.extend_from_slice(&1u32.to_le_bytes()); // 1 block
        output.extend_from_slice(&(compressed.len() as u32).to_le_bytes());
        output.extend_from_slice(&compressed);
        return Ok(output);
    }

    // Split into blocks and compress in parallel
    let blocks: Vec<&[u8]> = data.chunks(BSC_BLOCK_SIZE).collect();
    let num_blocks = blocks.len();

    let compressed_blocks: Vec<Result<Vec<u8>>> = blocks
        .par_iter()
        .map(|block| block_compressor(block))
        .collect();

    // Check for errors
    let mut output = Vec::new();
    output.extend_from_slice(&(num_blocks as u32).to_le_bytes());

    for result in compressed_blocks {
        let block = result?;
        output.extend_from_slice(&(block.len() as u32).to_le_bytes());
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

/// Decompress multi-block BSC-compressed data (from compress_parallel).
///
/// Blocks are decompressed in parallel using rayon, then concatenated in order.
/// Format: [num_blocks: u32][block_compressed_len: u32, block_data]...
pub fn decompress_parallel(data: &[u8]) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    if data.is_empty() {
        return Ok(vec![]);
    }

    if data.len() < 4 {
        anyhow::bail!("BSC parallel: data too small for header");
    }

    let num_blocks = super::read_le_u32(data, 0)? as usize;

    // First pass: collect block slices (sequential, just pointer math)
    let mut offset = 4;
    let mut block_slices = Vec::with_capacity(num_blocks);
    for _ in 0..num_blocks {
        if offset + 4 > data.len() {
            anyhow::bail!("BSC parallel: truncated block length");
        }
        let block_len = super::read_le_u32(data, offset)? as usize;
        offset += 4;

        if offset + block_len > data.len() {
            anyhow::bail!("BSC parallel: truncated block data");
        }
        block_slices.push(&data[offset..offset + block_len]);
        offset += block_len;
    }

    // Second pass: decompress all blocks in parallel
    let decompressed_blocks: Vec<Result<Vec<u8>>> = block_slices
        .par_iter()
        .map(|block| decompress(block))
        .collect();

    // Concatenate in order
    let total_size: usize = decompressed_blocks.iter()
        .filter_map(|r| r.as_ref().ok())
        .map(|b| b.len())
        .sum();
    let mut output = Vec::with_capacity(total_size);
    for result in decompressed_blocks {
        output.extend_from_slice(&result?);
    }

    Ok(output)
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
    ensure_initialized();

    if data.is_empty() {
        return Ok(vec![]);
    }

    if data.len() < LIBBSC_HEADER_SIZE {
        anyhow::bail!("BSC: compressed data too small (need at least {} bytes for header)", LIBBSC_HEADER_SIZE);
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
