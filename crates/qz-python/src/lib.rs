use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::path::PathBuf;
use std::sync::Once;

use qz_lib::cli::{CompressConfig, DecompressConfig, QualityMode, VerifyConfig};
use qz_lib::compression;

// Pooling global allocator for the Python extension, matching the `qz` CLI
// (see qz-cli/src/main.rs). Compress allocates a fresh ~25 MB output buffer per
// BSC block plus libsais workspace and frees it; under glibc malloc those large
// allocations are mmap'd/munmap'd to the kernel every block, so each one
// re-faults and re-zeros its pages (~44% of compress CPU was in the kernel).
// mimalloc recycles them in userspace → ~18% faster compress from Python too.
//
// Setting a Rust `#[global_allocator]` in a cdylib only governs Rust-side
// allocations (qz-lib's buffers); CPython keeps its own allocator for Python
// objects. PyO3 copies data across the boundary (e.g. building `bytes` from a
// Rust `Vec`), so no buffer is ever allocated by one allocator and freed by the
// other — the swap is safe.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Run `f` (which releases the GIL internally and returns a `Result`),
/// converting both error returns and panics into a catchable `PyRuntimeError`.
/// The `Ok` value is forwarded to the caller (e.g. a `VerifyResult`); use `()`
/// for the side-effecting compress/decompress paths.
///
/// Without the panic guard, a panic in qz-lib would propagate to the C/Python
/// boundary as PyO3's `PanicException` — a `BaseException` that a normal
/// `except Exception:` won't catch — or, under `panic="abort"`, abort the whole
/// interpreter. Catching it here yields an ordinary, catchable exception.
fn run_guarded<T, F, E>(f: F) -> PyResult<T>
where
    F: FnOnce() -> Result<T, E>,
    E: std::fmt::Display,
{
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)) {
        // `{e:#}` is anyhow's clean single-line cause chain — Debug (`{:?}`) would leak the
        // internal representation and (when enabled) a backtrace into the Python exception.
        Ok(result) => result.map_err(|e| PyRuntimeError::new_err(format!("{e:#}"))),
        Err(panic) => {
            let msg = panic
                .downcast_ref::<&str>()
                .map(|s| s.to_string())
                .or_else(|| panic.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "unknown panic".to_string());
            Err(PyRuntimeError::new_err(format!("qz panicked: {msg}")))
        }
    }
}

/// Install a tracing subscriber once per process so `info!`/`warn!` from
/// qz-lib reach stderr instead of being silently dropped. Honors RUST_LOG.
fn ensure_tracing_initialized() {
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        use tracing_subscriber::EnvFilter;
        let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
        let _ = tracing_subscriber::fmt()
            .with_env_filter(filter)
            .with_writer(std::io::stderr)
            .try_init();
    });
}

/// `input` to `compress` accepts either a single FASTQ path (single-end) or a
/// list of paths (paired-end: `[R1, R2]`). PyO3 tries the variants in order, so a
/// Python `str` extracts as `Single` and a `list`/tuple as `Multi`; the input
/// count (1 single-end, 2 paired) is validated by the library.
#[derive(FromPyObject)]
enum InputArg {
    Single(String),
    Multi(Vec<String>),
}

impl InputArg {
    fn into_paths(self) -> Vec<PathBuf> {
        match self {
            InputArg::Single(s) => vec![PathBuf::from(s)],
            InputArg::Multi(v) => v.into_iter().map(PathBuf::from).collect(),
        }
    }
}

/// Compress a FASTQ file to QZ format.
///
/// Gzipped input is auto-detected from file contents.
///
/// Args:
///     input: Path to the input FASTQ file (single-end), or a list of two paths
///         ``[R1, R2]`` for paired-end. A one-element list is single-end.
///     output: Path to output QZ archive
///     quality_mode: Quality mode - "lossless", "illumina-bin", or "discard" (default: "lossless")
///     ultra: Ultra compression level 1-3, or 0 for auto (default: None = disabled)
///     fasta: Input is FASTA format (default: False)
///     no_quality: Discard quality scores (default: False)
///     threads: Number of threads (0 = auto, default: 0). NOTE: this initializes
///         rayon's process-global pool on the first qz call; later calls with a
///         different count are ignored (the first call in the process wins).
///     working_dir: Working directory for temporary files (default: ".")
///     force: Overwrite output if it exists (default: False)
///     fast: Faster compression and decompression at a small ratio cost
///         (static QLFC entropy coder for header/sequence streams) (default: False)
///
/// Example:
///     >>> import qz
///     >>> qz.compress("reads.fastq", "compressed.qz")
///     >>> qz.compress(["R1.fastq", "R2.fastq"], "sample.qz")   # paired-end
///     >>> qz.compress("reads.fastq", "out.qz", ultra=3)
///     >>> qz.compress("reads.fastq", "out.qz", fast=True)
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
    fast = false,
))]
#[allow(clippy::too_many_arguments)]
fn compress(
    py: Python<'_>,
    input: InputArg,
    output: String,
    quality_mode: &str,
    ultra: Option<u8>,
    fasta: bool,
    no_quality: bool,
    threads: usize,
    working_dir: &str,
    force: bool,
    fast: bool,
) -> PyResult<()> {
    ensure_tracing_initialized();
    let qmode = parse_quality_mode(quality_mode)?;

    // The ultra level is 1-3, or 0 for auto; reject out-of-range values
    // instead of clamping.
    if let Some(level) = ultra
        && level != 0
        && !(1..=3).contains(&level)
    {
        return Err(PyRuntimeError::new_err(format!(
            "ultra level must be 1, 2, 3, or omitted for auto (got {level})"
        )));
    }

    // One path = single-end, two = paired-end. The library validates the count
    // (it bails on 0 or >2) and dispatches paired by `input.len()`.
    let input_paths = input.into_paths();

    let mut config = CompressConfig {
        input: input_paths,
        output: PathBuf::from(output),
        working_dir: PathBuf::from(working_dir),
        threads,
        // FASTA inputs carry no quality data; force quality off so the Auto
        // quality compressor doesn't try (and fail) to encode an empty
        // quality alphabet on large inputs.
        no_quality: no_quality || fasta || qmode == QualityMode::Discard,
        fasta,
        quality_mode: qmode,
        ultra,
        force,
        ..CompressConfig::default()
    };
    config.advanced.bsc_static |= fast;

    // Release the GIL during compression. Without this, a multi-minute
    // compress call would block the entire Python interpreter — no other
    // thread runs, signals (Ctrl-C) are not delivered, and any concurrent
    // Python work effectively halts.
    run_guarded(|| py.allow_threads(|| compression::compress(&config)))
}

/// Decompress a QZ archive to FASTQ format.
///
/// Args:
///     input: Path to input QZ archive
///     output: Path to output FASTQ file
///     working_dir: Working directory for temporary files (default: ".")
///     threads: Number of threads (0 = auto, default: 0). NOTE: this initializes
///         rayon's process-global pool on the first qz call; later calls with a
///         different count are ignored (the first call in the process wins).
///     gzipped: Output gzipped FASTQ (default: False)
///     gzip_level: Gzip compression level 0-9 (default: 6)
///     force: Overwrite the output file if it exists (default: False)
///     interleaved: For a paired or paired-reference archive, write ONE interleaved
///         FASTQ (R1/R2 records alternating) to `output` instead of rejecting it.
///         Mutually exclusive with `split`. (default: False)
///     split: For a paired or paired-reference archive, write the two mates to
///         SEPARATE files, using `output` as a prefix: `<output>_R1.fastq` and
///         `<output>_R2.fastq` (the same convention as `qz decompress -o <prefix>`).
///         Only valid for two-mate archives. Mutually exclusive with `interleaved`.
///         (default: False)
///
/// A two-mate archive (paired or paired-reference) needs either `split=True`
/// (separate `_R1`/`_R2` files) or `interleaved=True` (one stream); with neither,
/// it errors, since this API is otherwise single-output.
///
/// Example:
///     >>> import qz
///     >>> qz.decompress("compressed.qz", "reads.fastq")
///     >>> qz.decompress("compressed.qz", "reads.fastq.gz", gzipped=True)
///     >>> qz.decompress("paired.qz", "sample", split=True)            # sample_R1.fastq + sample_R2.fastq
///     >>> qz.decompress("paired.qz", "interleaved.fastq", interleaved=True)
#[pyfunction]
#[pyo3(signature = (
    input,
    output,
    working_dir = ".",
    threads = 0,
    gzipped = false,
    gzip_level = 6,
    force = false,
    interleaved = false,
    split = false,
))]
#[allow(clippy::too_many_arguments)]
fn decompress(
    py: Python<'_>,
    input: String,
    output: String,
    working_dir: &str,
    threads: usize,
    gzipped: bool,
    gzip_level: u32,
    force: bool,
    interleaved: bool,
    split: bool,
) -> PyResult<()> {
    ensure_tracing_initialized();
    if gzip_level > 9 {
        return Err(PyRuntimeError::new_err(format!(
            "gzip_level must be 0-9, got {gzip_level}"
        )));
    }
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
        force,
    };

    if interleaved && split {
        return Err(PyRuntimeError::new_err(
            "interleaved=True and split=True are mutually exclusive: interleaved writes ONE \
             stream, split writes separate _R1/_R2 files",
        ));
    }

    // Interleaved decode merges a two-mate archive (paired type 1 / paired-reference
    // type 2) into ONE output stream, so it IS single-output: route it here instead of
    // the two-file rejection below. The library validates the archive shape (single-end
    // / cluster → clear error).
    if interleaved {
        return run_guarded(|| py.allow_threads(|| compression::decompress_interleaved(&config)));
    }

    // Split decode writes the two mates to SEPARATE files, using `output` as a prefix
    // (`<output>_R1.fastq` / `<output>_R2.fastq`) exactly like `qz decompress -o <prefix>`.
    // Only paired (type 1) and paired-reference (type 2) produce two files; `compression::
    // decompress` already expands the single output[0] prefix into the two mate paths, so
    // we just skip the two-file rejection below and let it run. Single-output archives
    // (single-end, single-end-reference, single-end cluster) and paired cluster (which
    // needs two exact paths, not a prefix) error with a clear hint.
    if split {
        if qz_lib::compression::is_paired_archive(&config.input)
            .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?
            || qz_lib::compression::is_reference_archive(&config.input)
                .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?
        {
            return run_guarded(|| py.allow_threads(|| compression::decompress(&config)));
        }
        if qz_lib::compression::is_cluster_paired_archive(&config.input)
            .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?
        {
            return Err(PyRuntimeError::new_err(
                "paired cluster archives must be decompressed with the CLI \
                 (`qz decompress -o <R1> <R2>`), not split",
            ));
        }
        return Err(PyRuntimeError::new_err(
            "split=True is only valid for paired (type 1) or paired-reference (type 2) \
             archives; this archive is single-output, so omit split",
        ));
    }

    // The non-interleaved Python decode API is single-output; archives that reconstruct
    // two FASTQ files — paired (type 1), reference (type 2), and PAIRED cluster (type 3
    // with mate-2 entries) — must go through the CLI (or `interleaved=True` above for
    // paired/reference). Reject them explicitly rather than letting
    // `compression::decompress` write two files or fail deep in decode with a confusing
    // message. Single-end cluster (type 3, one output) and single-end reference (type 4,
    // one output) decode fine through this single-output API and are intentionally NOT
    // rejected. A legacy (non-v5) archive returns `Ok(false)` from these predicates and
    // instead hits the recompress-hint bail inside `compression::decompress`.
    if qz_lib::compression::is_reference_archive(&config.input)
        .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?
    {
        return Err(PyRuntimeError::new_err(
            "reference archives must be decompressed with the CLI (`qz decompress -o <prefix>`), \
             or pass interleaved=True for a single interleaved FASTQ",
        ));
    }
    if qz_lib::compression::is_paired_archive(&config.input)
        .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?
    {
        return Err(PyRuntimeError::new_err(
            "paired archives must be decompressed with the CLI (`qz decompress -o <prefix>`), \
             or pass interleaved=True for a single interleaved FASTQ",
        ));
    }
    if qz_lib::compression::is_cluster_paired_archive(&config.input)
        .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?
    {
        return Err(PyRuntimeError::new_err(
            "paired cluster archives must be decompressed with the CLI (`qz decompress -o <prefix>`)",
        ));
    }

    // See compress() above for the GIL rationale.
    run_guarded(|| py.allow_threads(|| compression::decompress(&config)))
}

/// Verify a QZ archive's integrity, without writing any output file.
///
/// Unlike `decompress` (which is single-output only), verify works for **every**
/// archive type — single-end, paired, reference, and cluster — because it never
/// reconstructs the output to disk.
///
/// Returns a dict of verification stats on success. Raises `RuntimeError` if the
/// archive is corrupt, truncated, tampered with, an unsupported/legacy format,
/// or otherwise fails its checks — including a `--cluster` archive whose
/// order-independent multiset checksum no longer matches.
///
/// Args:
///     input: Path to the QZ archive to verify
///     fast: Fast mode — walk every per-block CRC32 without a full decode
///         (catches bit-rot/truncation quickly). The default deep mode fully
///         reconstructs the FASTQ through a CRC32. (default: False)
///     threads: Number of threads (0 = auto, default: 0). NOTE: this initializes
///         rayon's process-global pool on the first qz call; later calls with a
///         different count are ignored (the first call in the process wins).
///     working_dir: Working directory for temporary files (default: ".")
///
/// Returns:
///     dict with keys: ``ok`` (always True when it returns), ``mode``
///     ("deep" | "fast"), ``num_reads``, ``encoding_type``, ``encoding``,
///     ``header_compressor``, ``quality_compressor``, ``headers_compressed_len``,
///     ``sequences_compressed_len``, ``qualities_compressed_len``, ``crc32``
///     (deep only; 0 in fast mode), ``r2_crc32`` (int for two-mate deep verify,
///     else None), ``blocks_verified`` (fast only; 0 in deep mode),
///     ``total_bytes``, ``elapsed_secs``.
///
/// Example:
///     >>> import qz
///     >>> info = qz.verify("reads.qz")
///     >>> info["num_reads"], info["mode"]
///     (1000000, 'deep')
///     >>> qz.verify("reads.qz", fast=True)["blocks_verified"]
///     42
#[pyfunction]
#[pyo3(signature = (
    input,
    fast = false,
    threads = 0,
    working_dir = ".",
))]
fn verify(
    py: Python<'_>,
    input: String,
    fast: bool,
    threads: usize,
    working_dir: &str,
) -> PyResult<PyObject> {
    ensure_tracing_initialized();

    let num_threads = if threads == 0 {
        qz_lib::cli::num_cpus()
    } else {
        threads
    };

    let config = VerifyConfig {
        input: PathBuf::from(input),
        working_dir: PathBuf::from(working_dir),
        num_threads,
        fast,
    };

    // Release the GIL for the (potentially long) deep verify; see compress().
    let result = run_guarded(|| py.allow_threads(|| compression::verify(&config)))?;

    let mode = match result.mode {
        compression::VerifyMode::Deep => "deep",
        compression::VerifyMode::Fast => "fast",
    };

    let dict = PyDict::new(py);
    dict.set_item("ok", true)?;
    dict.set_item("mode", mode)?;
    dict.set_item("num_reads", result.num_reads)?;
    dict.set_item("encoding_type", result.encoding_type)?;
    dict.set_item(
        "encoding",
        compression::encoding_type_display_name(result.encoding_type),
    )?;
    dict.set_item("header_compressor", format!("{:?}", result.header_compressor))?;
    dict.set_item(
        "quality_compressor",
        format!("{:?}", result.quality_compressor),
    )?;
    dict.set_item("headers_compressed_len", result.headers_compressed_len)?;
    dict.set_item("sequences_compressed_len", result.sequences_compressed_len)?;
    dict.set_item("qualities_compressed_len", result.qualities_compressed_len)?;
    dict.set_item("crc32", result.crc32)?;
    dict.set_item("r2_crc32", result.r2_crc32)?;
    dict.set_item("blocks_verified", result.blocks_verified)?;
    dict.set_item("total_bytes", result.total_bytes)?;
    dict.set_item("elapsed_secs", result.elapsed_secs)?;

    Ok(dict.into_any().unbind())
}

/// Merge several QZ archives into one (the `qz merge` primitive).
///
/// Each input's compressed block frames are relocated verbatim — the encoder is
/// never re-run — and the directory is rebuilt with shifted offsets and
/// renumbered chunks, so the result is a single valid archive. The reads appear
/// in the order the inputs are listed.
///
/// All inputs must share one archive type and config: single-end (type 0) and
/// paired (type 1) merge directly; reference archives (type 2/4) need
/// ``reference`` (the same FASTA they were compressed against), because their
/// per-shard coverage globals are re-derived rather than concatenated; cluster
/// archives (type 3) cannot be merged. These constraints are enforced by the
/// library and surface as ``RuntimeError``.
///
/// Args:
///     inputs: List of QZ archive paths to merge (>=1), in output order.
///     output: Path to the merged QZ archive.
///     reference: FASTA path required to merge reference archives (type 2/4);
///         omit for single-end/paired (default: None).
///     force: Overwrite the output file if it already exists (default: False).
///
/// Example:
///     >>> import qz
///     >>> qz.merge(["a.qz", "b.qz"], "merged.qz")          # single-end or paired
///     >>> qz.merge(shards, "merged.qz", reference="ref.fa")  # reference (type 2/4)
#[pyfunction]
#[pyo3(signature = (
    inputs,
    output,
    reference = None,
    force = false,
))]
fn merge(
    py: Python<'_>,
    inputs: Vec<String>,
    output: String,
    reference: Option<String>,
    force: bool,
) -> PyResult<()> {
    ensure_tracing_initialized();

    let input_paths: Vec<PathBuf> = inputs.into_iter().map(PathBuf::from).collect();
    let output_path = PathBuf::from(output);

    // Release the GIL: merge is I/O heavy (verbatim frame copy of whole archives)
    // and reference merge re-derives coverage globals on the rayon pool. See
    // compress() for the full rationale.
    run_guarded(|| {
        py.allow_threads(|| match reference {
            Some(fasta) => compression::merge_reference_archives_to_path(
                &input_paths,
                &PathBuf::from(fasta),
                &output_path,
                force,
            ),
            None => compression::merge_archives_to_path(&input_paths, &output_path, force),
        })
    })
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
    m.add_function(wrap_pyfunction!(verify, m)?)?;
    m.add_function(wrap_pyfunction!(merge, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
