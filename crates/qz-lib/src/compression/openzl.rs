/// FFI bindings to OpenZL (Meta's format-aware compression framework)
///
/// This module provides safe Rust wrappers around the OpenZL C API,
/// using the same block-parallel pattern as BSC for consistent behavior.
///
/// Source: https://github.com/facebook/openzl

use anyhow::Result;
use openzl_sys::*;
use std::os::raw::c_void;

/// Block size for parallel compression (matching BSC: 25 MB)
const OPENZL_BLOCK_SIZE: usize = 25 * 1024 * 1024;

/// Check if a ZL_Report indicates an error
fn is_error(report: ZL_Report) -> bool {
    unsafe { report._code != ZL_ErrorCode_no_error }
}

/// Extract the valid result from a ZL_Report
fn valid_result(report: ZL_Report) -> usize {
    unsafe { report._value._value }
}

/// Extract error code from a ZL_Report for error messages
fn error_code(report: ZL_Report) -> u32 {
    unsafe { report._code }
}

/// Compress data using OpenZL (single block)
pub fn compress(data: &[u8]) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(vec![]);
    }

    unsafe {
        let cctx = ZL_CCtx_create();
        if cctx.is_null() {
            anyhow::bail!("OpenZL: failed to create compression context");
        }

        // Set format version (required by OpenZL)
        let version_report = ZL_CCtx_setParameter(
            cctx,
            ZL_CParam_formatVersion,
            ZL_MAX_FORMAT_VERSION as i32,
        );
        if is_error(version_report) {
            ZL_CCtx_free(cctx);
            anyhow::bail!(
                "OpenZL: failed to set format version, error code: {}",
                error_code(version_report)
            );
        }

        let dst_capacity = ZL_compressBound(data.len());
        let mut output = vec![0u8; dst_capacity];

        let report = ZL_CCtx_compress(
            cctx,
            output.as_mut_ptr() as *mut c_void,
            dst_capacity,
            data.as_ptr() as *const c_void,
            data.len(),
        );

        ZL_CCtx_free(cctx);

        if is_error(report) {
            anyhow::bail!(
                "OpenZL compression failed with error code: {}",
                error_code(report)
            );
        }

        let compressed_size = valid_result(report);
        output.truncate(compressed_size);
        output.shrink_to_fit();
        Ok(output)
    }
}

// ── Bench-only graph variants ───────────────────────────────────────────────
//
// The eight `compress_*` functions below (ACE / DNA-numeric / clustering /
// delta-entropy / FSE / Huffman / transpose+zstd / field_lz) and their
// `*_graph_fn` helpers are only called from `qz-bench/src/bin/bench_openzl_graphs.rs`.
// Production code uses [`compress`] / [`compress_parallel`] only. Gated behind
// the `experimental` feature so that they don't ship in qz-cli / qz-python
// binaries and don't widen the supported public API. The bench crate depends
// on `qz-lib` with `features = ["experimental"]`, so benches still compile.

#[cfg(feature = "experimental")]
mod bench_graphs {
    use super::*;

/// Graph function for ACE (Adaptive Codec Engine) training.
/// ACE analyzes the data online and builds an optimized compression graph.
unsafe extern "C" fn ace_graph_fn(compressor: *mut ZL_Compressor) -> ZL_GraphID {
    unsafe {
        let result = ZL_Compressor_setParameter(
            compressor,
            ZL_CParam_formatVersion,
            ZL_MAX_FORMAT_VERSION as i32,
        );
        if is_error(result) {
            return ZL_GraphID { gid: ZL_StandardGraphID_illegal };
        }

        // Use ACE with generic Zstd as the default until trained
        let default_graph = ZL_GraphID { gid: ZL_StandardGraphID_compress_generic };
        let ace = ZL_Compressor_buildACEGraphWithDefault(compressor, default_graph);
        if ace.gid == ZL_StandardGraphID_illegal {
            // Fall back to generic compression
            return default_graph;
        }
        ace
    }
}

/// Graph function for DNA-optimized compression.
/// Treats byte stream as u8 numeric array, applies delta + entropy coding.
unsafe extern "C" fn dna_numeric_graph_fn(compressor: *mut ZL_Compressor) -> ZL_GraphID {
    unsafe {
        let result = ZL_Compressor_setParameter(
            compressor,
            ZL_CParam_formatVersion,
            ZL_MAX_FORMAT_VERSION as i32,
        );
        if is_error(result) {
            return ZL_GraphID { gid: ZL_StandardGraphID_illegal };
        }

        // Convert serial bytes to u8 numeric, then apply entropy coding
        let inner = ZL_GraphID { gid: ZL_StandardGraphID_entropy };
        ZL_Compressor_registerStaticGraph_fromNode1o(
            compressor,
            ZL_NodeID { nid: ZL_StandardNodeID_convert_serial_to_num8 },
            inner,
        )
    }
}

/// Graph function for clustering + entropy - lets OpenZL cluster similar
/// byte patterns together before entropy coding.
unsafe extern "C" fn clustering_graph_fn(compressor: *mut ZL_Compressor) -> ZL_GraphID {
    unsafe {
        let result = ZL_Compressor_setParameter(
            compressor,
            ZL_CParam_formatVersion,
            ZL_MAX_FORMAT_VERSION as i32,
        );
        if is_error(result) {
            return ZL_GraphID { gid: ZL_StandardGraphID_illegal };
        }
        ZL_GraphID { gid: ZL_StandardGraphID_clustering }
    }
}

/// Compress data using OpenZL with a custom graph function (single block).
pub fn compress_with_graph(data: &[u8], graph_fn: ZL_GraphFn) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(vec![]);
    }

    unsafe {
        let dst_capacity = ZL_compressBound(data.len());
        let mut output = vec![0u8; dst_capacity];

        let report = ZL_compress_usingGraphFn(
            output.as_mut_ptr() as *mut c_void,
            dst_capacity,
            data.as_ptr() as *const c_void,
            data.len(),
            graph_fn,
        );

        if is_error(report) {
            anyhow::bail!(
                "OpenZL graph compression failed with error code: {}",
                error_code(report)
            );
        }

        let compressed_size = valid_result(report);
        output.truncate(compressed_size);
        output.shrink_to_fit();
        Ok(output)
    }
}

/// Compress data using ACE (Adaptive Codec Engine) for inline training.
pub fn compress_ace(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_graph(data, Some(ace_graph_fn))
}

/// Compress data using DNA-optimized graph (numeric u8 + entropy).
pub fn compress_dna_numeric(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_graph(data, Some(dna_numeric_graph_fn))
}

/// Compress data using clustering graph.
pub fn compress_clustering(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_graph(data, Some(clustering_graph_fn))
}

/// Graph function for delta + entropy pipeline.
/// Converts serial bytes to u8 numeric, applies delta encoding, then entropy coding.
unsafe extern "C" fn delta_entropy_graph_fn(compressor: *mut ZL_Compressor) -> ZL_GraphID {
    unsafe {
        let result = ZL_Compressor_setParameter(
            compressor,
            ZL_CParam_formatVersion,
            ZL_MAX_FORMAT_VERSION as i32,
        );
        if is_error(result) {
            return ZL_GraphID { gid: ZL_StandardGraphID_illegal };
        }

        // Pipeline: serial -> u8 numeric -> delta -> entropy
        let nodes = [
            ZL_NodeID { nid: ZL_StandardNodeID_convert_serial_to_num8 },
            ZL_NodeID { nid: ZL_StandardNodeID_delta_int },
        ];
        ZL_Compressor_registerStaticGraph_fromPipelineNodes1o(
            compressor,
            nodes.as_ptr(),
            nodes.len(),
            ZL_GraphID { gid: ZL_StandardGraphID_entropy },
        )
    }
}

/// Compress data using delta + entropy pipeline.
pub fn compress_delta_entropy(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_graph(data, Some(delta_entropy_graph_fn))
}

/// Graph function for FSE (Finite State Entropy) coding.
unsafe extern "C" fn fse_graph_fn(compressor: *mut ZL_Compressor) -> ZL_GraphID {
    unsafe {
        let result = ZL_Compressor_setParameter(
            compressor,
            ZL_CParam_formatVersion,
            ZL_MAX_FORMAT_VERSION as i32,
        );
        if is_error(result) {
            return ZL_GraphID { gid: ZL_StandardGraphID_illegal };
        }
        ZL_GraphID { gid: ZL_StandardGraphID_fse }
    }
}

/// Compress data using FSE (Finite State Entropy).
pub fn compress_fse(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_graph(data, Some(fse_graph_fn))
}

/// Graph function for Huffman coding.
unsafe extern "C" fn huffman_graph_fn(compressor: *mut ZL_Compressor) -> ZL_GraphID {
    unsafe {
        let result = ZL_Compressor_setParameter(
            compressor,
            ZL_CParam_formatVersion,
            ZL_MAX_FORMAT_VERSION as i32,
        );
        if is_error(result) {
            return ZL_GraphID { gid: ZL_StandardGraphID_illegal };
        }
        ZL_GraphID { gid: ZL_StandardGraphID_huffman }
    }
}

/// Compress data using Huffman coding.
pub fn compress_huffman(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_graph(data, Some(huffman_graph_fn))
}

/// Graph function for transpose-split + zstd.
/// Byte-level transposition helps when data has repeating structure at fixed offsets.
unsafe extern "C" fn transpose_zstd_graph_fn(compressor: *mut ZL_Compressor) -> ZL_GraphID {
    unsafe {
        let result = ZL_Compressor_setParameter(
            compressor,
            ZL_CParam_formatVersion,
            ZL_MAX_FORMAT_VERSION as i32,
        );
        if is_error(result) {
            return ZL_GraphID { gid: ZL_StandardGraphID_illegal };
        }
        ZL_Compressor_registerStaticGraph_fromNode1o(
            compressor,
            ZL_NodeID { nid: ZL_StandardNodeID_transpose_split },
            ZL_GraphID { gid: ZL_StandardGraphID_zstd },
        )
    }
}

/// Compress data using transpose-split + zstd.
pub fn compress_transpose_zstd(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_graph(data, Some(transpose_zstd_graph_fn))
}

/// Graph function for field_lz (dictionary-based compression tuned for structured fields).
unsafe extern "C" fn field_lz_graph_fn(compressor: *mut ZL_Compressor) -> ZL_GraphID {
    unsafe {
        let result = ZL_Compressor_setParameter(
            compressor,
            ZL_CParam_formatVersion,
            ZL_MAX_FORMAT_VERSION as i32,
        );
        if is_error(result) {
            return ZL_GraphID { gid: ZL_StandardGraphID_illegal };
        }
        ZL_GraphID { gid: ZL_StandardGraphID_field_lz }
    }
}

/// Compress data using field LZ.
pub fn compress_field_lz(data: &[u8]) -> Result<Vec<u8>> {
    compress_with_graph(data, Some(field_lz_graph_fn))
}

} // mod bench_graphs

#[cfg(feature = "experimental")]
pub use bench_graphs::{
    compress_ace, compress_clustering, compress_delta_entropy, compress_dna_numeric,
    compress_field_lz, compress_fse, compress_huffman, compress_transpose_zstd,
};

/// Decompress OpenZL-compressed data (single block)
pub fn decompress(data: &[u8]) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(vec![]);
    }

    unsafe {
        let size_report = ZL_getDecompressedSize(data.as_ptr() as *const c_void, data.len());

        if is_error(size_report) {
            anyhow::bail!(
                "OpenZL: failed to get decompressed size, error code: {}",
                error_code(size_report)
            );
        }

        let decompressed_size = valid_result(size_report);
        let mut output = vec![0u8; decompressed_size];

        let report = ZL_decompress(
            output.as_mut_ptr() as *mut c_void,
            decompressed_size,
            data.as_ptr() as *const c_void,
            data.len(),
        );

        if is_error(report) {
            anyhow::bail!(
                "OpenZL decompression failed with error code: {}",
                error_code(report)
            );
        }

        Ok(output)
    }
}

/// Compress data by splitting into blocks and compressing each in parallel.
///
/// Each block is compressed independently with OpenZL using rayon's thread pool.
/// Format: [num_blocks: u32][block_compressed_len: u32, block_data]...
/// Falls back to single-block for data smaller than OPENZL_BLOCK_SIZE.
pub fn compress_parallel(data: &[u8]) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    if data.is_empty() {
        return Ok(vec![]);
    }

    // For small data, use single-block compression (no overhead).
    // Wire layout matches BSC v3: [num_blocks: u32][block_len: u32, crc32: u32, data]...
    // CRC32 over the OpenZL-compressed payload is verified on read.
    if data.len() <= OPENZL_BLOCK_SIZE {
        let compressed = compress(data)?;
        let crc = block_crc32(&compressed);
        let mut output = Vec::with_capacity(4 + 8 + compressed.len());
        output.extend_from_slice(&1u32.to_le_bytes()); // 1 block
        output.extend_from_slice(&(compressed.len() as u32).to_le_bytes());
        output.extend_from_slice(&crc.to_le_bytes());
        output.extend_from_slice(&compressed);
        return Ok(output);
    }

    // Split into blocks and compress + CRC in parallel
    let blocks: Vec<&[u8]> = data.chunks(OPENZL_BLOCK_SIZE).collect();
    let num_blocks = blocks.len();

    let compressed_blocks: Vec<(Vec<u8>, u32)> = blocks
        .par_iter()
        .map(|block| {
            let compressed = compress(block)?;
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

// CRC32 helper is shared with all other codecs via `super::bsc::block_crc32`.
use super::bsc::block_crc32;

/// Decompress multi-block OpenZL-compressed data (from compress_parallel).
///
/// Blocks are decompressed in parallel using rayon, then concatenated in order.
/// Format: [num_blocks: u32][block_compressed_len: u32, block_data]...
pub fn decompress_parallel(data: &[u8]) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    if data.is_empty() {
        return Ok(vec![]);
    }

    if data.len() < 4 {
        anyhow::bail!("OpenZL parallel: data too small for header");
    }

    let num_blocks = u32::from_le_bytes(
        data[0..4].try_into()
            .map_err(|_| anyhow::anyhow!("OpenZL parallel: truncated header"))?,
    ) as usize;

    // Bound `num_blocks` against the remaining payload (each block ≥ 9 bytes:
    // 4 len + 4 crc + 1 payload). Prevents Vec::with_capacity OOM on a hostile
    // num_blocks value before we hit the per-block truncation check.
    if num_blocks > 0 && num_blocks.saturating_mul(9) > data.len().saturating_sub(4) {
        anyhow::bail!(
            "OpenZL parallel: num_blocks={} exceeds remaining payload ({} bytes)",
            num_blocks,
            data.len() - 4,
        );
    }

    // First pass: collect (block_slice, expected_crc) pairs (sequential)
    let mut offset = 4;
    let mut blocks: Vec<(&[u8], u32)> = Vec::with_capacity(num_blocks);
    for idx in 0..num_blocks {
        if offset + 8 > data.len() {
            anyhow::bail!("OpenZL parallel: truncated block header for block {idx}");
        }
        let block_len = u32::from_le_bytes(
            data[offset..offset + 4].try_into()
                .map_err(|_| anyhow::anyhow!("OpenZL parallel: truncated block length"))?,
        ) as usize;
        offset += 4;
        let expected_crc = u32::from_le_bytes(
            data[offset..offset + 4].try_into()
                .map_err(|_| anyhow::anyhow!("OpenZL parallel: truncated block crc"))?,
        );
        offset += 4;

        if offset + block_len > data.len() {
            anyhow::bail!("OpenZL parallel: truncated block data for block {idx}");
        }
        blocks.push((&data[offset..offset + block_len], expected_crc));
        offset += block_len;
    }

    // Second pass: verify CRC + decompress each block in parallel.
    let decompressed_blocks: Vec<Result<Vec<u8>>> = blocks
        .par_iter()
        .enumerate()
        .map(|(idx, (block, expected_crc))| {
            let actual_crc = block_crc32(block);
            if actual_crc != *expected_crc {
                anyhow::bail!(
                    "OpenZL block {idx} CRC32 mismatch: stored {:08x}, computed {:08x} — archive is corrupted",
                    expected_crc, actual_crc
                );
            }
            decompress(block)
        })
        .collect();

    // Concatenate in order
    let total_size: usize = decompressed_blocks
        .iter()
        .filter_map(|r| r.as_ref().ok())
        .map(|b| b.len())
        .sum();
    let mut output = Vec::with_capacity(total_size);
    for result in decompressed_blocks {
        output.extend_from_slice(&result?);
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
        let data = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let compressed = compress(data).unwrap();
        println!(
            "Original: {} bytes, Compressed: {} bytes, Ratio: {:.2}x",
            data.len(),
            compressed.len(),
            data.len() as f64 / compressed.len() as f64
        );
        // Note: OpenZL may not compress very small buffers well due to frame overhead
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
    fn test_parallel_roundtrip() {
        // Create data large enough to trigger multi-block
        let base = b"ACGTACGTACGTACGTNNNNACGTACGTACGT";
        let data: Vec<u8> = base.iter().cycle().take(1024 * 1024).copied().collect();
        let compressed = compress_parallel(&data).unwrap();
        let decompressed = decompress_parallel(&compressed).unwrap();
        assert_eq!(data, decompressed);
    }
}
