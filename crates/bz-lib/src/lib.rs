pub mod cli;
pub mod compression;
pub mod io;

pub use cli::{CompressConfig, DecompressConfig};
pub use compression::{compress, decompress};
