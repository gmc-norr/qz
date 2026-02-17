pub mod archive;
mod compress;
mod decompress;
pub mod streams;

use crate::cli::{CompressConfig, DecompressConfig};
use anyhow::Result;

pub fn compress(config: &CompressConfig) -> Result<()> {
    compress::compress(config)
}

pub fn decompress(config: &DecompressConfig) -> Result<()> {
    decompress::decompress(config)
}
