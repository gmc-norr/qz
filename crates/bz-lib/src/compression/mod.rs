pub mod archive;
mod compress;
mod decompress;
mod extract;
pub mod md_codec;
pub mod streams;
mod verify;

use crate::cli::{CompressConfig, DecompressConfig, ExtractConfig, VerifyConfig, VerifyResult};
use anyhow::Result;

/// NUMA decode-sharding surface: pre-scan a BZ archive's chunk layout, decode a
/// chunk range to a part file, and bound a worker's peak RSS. Used by `bz-cli`'s
/// NUMA driver (mirrors qz-lib's `read_chunk_layout` / `decode_chunk_range`).
pub use decompress::{
    decode_chunk_range, decode_peak_rss_bound, read_bz_chunk_layout, BzChunkLayout,
};

/// NUMA compress-sharding surface (mirror of the decode side): prescan a BAM's
/// per-chunk virtual-offset layout, then compress a disjoint chunk range to a
/// self-contained part the driver concatenates. Used by `bz-cli`'s NUMA driver.
pub use compress::{
    compress_chunk_range, compress_peak_rss_bound, read_bam_compress_layout,
    resolve_compress_level, BzCompressLayout,
};

pub fn compress(config: &CompressConfig) -> Result<()> {
    compress::compress(config)
}

pub fn decompress(config: &DecompressConfig) -> Result<()> {
    decompress::decompress(config)
}

pub fn extract(config: &ExtractConfig) -> Result<()> {
    extract::extract(config)
}

pub fn verify(config: &VerifyConfig) -> Result<VerifyResult> {
    verify::verify(config)
}
