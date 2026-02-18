pub mod archive;
mod compress;
mod decompress;
mod extract;
pub mod streams;

use crate::cli::{CompressConfig, DecompressConfig, ExtractConfig};
use anyhow::Result;

pub fn compress(config: &CompressConfig) -> Result<()> {
    compress::compress(config)
}

pub fn decompress(config: &DecompressConfig) -> Result<()> {
    decompress::decompress(config)
}

pub fn extract(config: &ExtractConfig) -> Result<()> {
    extract::extract(config)
}
