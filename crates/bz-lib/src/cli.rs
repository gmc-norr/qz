use std::path::PathBuf;

pub struct CompressConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub working_dir: PathBuf,
    pub threads: usize,
}

pub struct DecompressConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub working_dir: PathBuf,
    pub threads: usize,
}

pub struct ExtractConfig {
    pub input: PathBuf,
    pub output_prefix: String,
    pub working_dir: PathBuf,
    pub threads: usize,
}
