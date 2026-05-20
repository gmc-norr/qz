use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use std::path::PathBuf;
use std::sync::Once;

use qz_lib::cli::{CompressConfig, DecompressConfig, QualityMode};
use qz_lib::compression;

/// Install a tracing subscriber once per process so `info!`/`warn!` from
/// qz-lib reach stderr instead of being silently dropped. Honors RUST_LOG.
fn ensure_tracing_initialized() {
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        use tracing_subscriber::EnvFilter;
        let filter = EnvFilter::try_from_default_env()
            .unwrap_or_else(|_| EnvFilter::new("info"));
        let _ = tracing_subscriber::fmt()
            .with_env_filter(filter)
            .with_writer(std::io::stderr)
            .try_init();
    });
}

/// Compress a FASTQ file to QZ format.
///
/// Gzipped input is auto-detected from file contents.
///
/// Args:
///     input: Path to input FASTQ file
///     output: Path to output QZ archive
///     quality_mode: Quality mode - "lossless", "illumina-bin", or "discard" (default: "lossless")
///     ultra: Ultra compression level 1-5, or 0 for auto (default: None = disabled)
///     fasta: Input is FASTA format (default: False)
///     no_quality: Discard quality scores (default: False)
///     threads: Number of threads (0 = auto, default: 0)
///     working_dir: Working directory for temporary files (default: ".")
///     force: Overwrite output if it exists (default: False)
///
/// Example:
///     >>> import qz
///     >>> qz.compress("reads.fastq", "compressed.qz")
///     >>> qz.compress("reads.fastq", "out.qz", ultra=3)
#[pyfunction]
#[pyo3(signature = (
    input,
    output,
    quality_mode = "lossless",
    ultra = None,
    fasta = false,
    no_quality = false,
    threads = 0,
    working_dir = ".",
    force = false,
))]
#[allow(clippy::too_many_arguments)]
fn compress(
    py: Python<'_>,
    input: String,
    output: String,
    quality_mode: &str,
    ultra: Option<u8>,
    fasta: bool,
    no_quality: bool,
    threads: usize,
    working_dir: &str,
    force: bool,
) -> PyResult<()> {
    ensure_tracing_initialized();
    let qmode = parse_quality_mode(quality_mode)?;

    let config = CompressConfig {
        input: vec![PathBuf::from(input)],
        output: PathBuf::from(output),
        working_dir: PathBuf::from(working_dir),
        threads,
        no_quality: no_quality || qmode == QualityMode::Discard,
        fasta,
        quality_mode: qmode,
        ultra,
        force,
        ..CompressConfig::default()
    };

    // Release the GIL during compression. Without this, a multi-minute
    // compress call would block the entire Python interpreter — no other
    // thread runs, signals (Ctrl-C) are not delivered, and any concurrent
    // Python work effectively halts.
    py.allow_threads(|| compression::compress(&config))
        .map_err(|e| PyRuntimeError::new_err(format!("{:?}", e)))
}

/// Decompress a QZ archive to FASTQ format.
///
/// Args:
///     input: Path to input QZ archive
///     output: Path to output FASTQ file
///     working_dir: Working directory for temporary files (default: ".")
///     threads: Number of threads (0 = auto, default: 0)
///     gzipped: Output gzipped FASTQ (default: False)
///     gzip_level: Gzip compression level 0-9 (default: 6)
///
/// Example:
///     >>> import qz
///     >>> qz.decompress("compressed.qz", "reads.fastq")
///     >>> qz.decompress("compressed.qz", "reads.fastq.gz", gzipped=True)
#[pyfunction]
#[pyo3(signature = (
    input,
    output,
    working_dir = ".",
    threads = 0,
    gzipped = false,
    gzip_level = 6,
))]
fn decompress(
    py: Python<'_>,
    input: String,
    output: String,
    working_dir: &str,
    threads: usize,
    gzipped: bool,
    gzip_level: u32,
) -> PyResult<()> {
    ensure_tracing_initialized();
    let num_threads = if threads == 0 {
        qz_lib::cli::num_cpus()
    } else {
        threads
    };

    let config = DecompressConfig {
        input: PathBuf::from(input),
        output: vec![PathBuf::from(output)],
        working_dir: PathBuf::from(working_dir),
        num_threads,
        gzipped,
        gzip_level,
    };

    // See compress() above for the GIL rationale.
    py.allow_threads(|| compression::decompress(&config))
        .map_err(|e| PyRuntimeError::new_err(format!("{:?}", e)))
}

/// Get QZ version string.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

fn parse_quality_mode(mode: &str) -> PyResult<QualityMode> {
    match mode {
        "lossless" => Ok(QualityMode::Lossless),
        "illumina-bin" => Ok(QualityMode::IlluminaBin),
        "discard" => Ok(QualityMode::Discard),
        _ => Err(PyRuntimeError::new_err(format!(
            "Invalid quality mode: '{}'. Must be one of: lossless, illumina-bin, discard",
            mode
        ))),
    }
}

/// QZ - High-performance FASTQ compression
#[pymodule]
fn qz(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compress, m)?)?;
    m.add_function(wrap_pyfunction!(decompress, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
