pub mod cli;
pub mod compression;
pub mod io;

pub use cli::{
    AdvancedOptions, CompressConfig, DecompressConfig, ExtractConfig, FlattenScheme,
    QualityReduction, VerifyConfig, VerifyResult,
};
pub use compression::{compress, decompress, extract, verify};

/// True iff a CUDA device is present and queryable — re-exported from qz-lib (bz's BSC BWT
/// runs through qz-lib) so `bz --gpu require` can fail loud. cuda builds only.
#[cfg(feature = "cuda")]
pub fn cuda_device_available() -> bool {
    qz_lib::compression::bsc::cuda_device_available()
}
