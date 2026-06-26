use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;
use tracing::info;

use qz_lib::cli::{CompressConfig, DecompressConfig, QualityMode as LibQualityMode};

// Pooling global allocator. Compress allocates a fresh ~25 MB output buffer per
// BSC block (bsc.rs) plus libsais workspace, then frees it; under glibc malloc
// these large allocations are mmap'd and munmap'd back to the kernel every block,
// so each one re-faults and re-zeros its pages. Profiling showed ~44% of compress
// CPU time was in the kernel page-fault path. mimalloc recycles freed buffers in
// userspace, eliminating the recurring faults: ~18% faster compress, byte-identical
// output, no ratio change, portable (no NUMA/machine assumptions). Decompress is
// compute-bound (QLFC decode) and unaffected.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod numa;

#[derive(Parser)]
#[command(name = "qz")]
#[command(author = "QZ Contributors")]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(about = "High-performance FASTQ compression", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Enable debug mode: verbose logging, full backtraces, and system diagnostics on crash
    #[arg(long, global = true)]
    debug: bool,
}

#[derive(Subcommand)]
enum Commands {
    /// Compress FASTQ files
    Compress(CompressArgs),
    /// Decompress QZ archives
    Decompress(DecompressArgs),
    /// Verify archive integrity
    Verify(VerifyArgs),
    /// Merge compatible QZ archives into one (single-end, paired, or reference)
    Merge(MergeArgs),
    /// Build a strobealign reference index for reference-based compression
    Index(IndexArgs),
    /// Internal: one NUMA decode worker (hidden — not a stable interface).
    #[command(name = "_numa-decode-worker", hide = true)]
    NumaDecodeWorker(NumaWorkerArgs),
    /// Internal: one NUMA compress worker (hidden — not a stable interface).
    #[command(name = "_numa-compress-worker", hide = true)]
    NumaCompressWorker(NumaCompressWorkerArgs),
}

/// Quality mode (subset exposed in CLI)
#[derive(Clone, Copy, Debug, ValueEnum, PartialEq, Eq)]
enum CliQualityMode {
    /// Lossless quality preservation
    Lossless,
    /// Illumina 8-level binning
    IlluminaBin,
    /// Discard quality scores entirely
    Discard,
}

/// GPU acceleration mode (cuda builds only — the flag is absent from CPU builds).
#[cfg(feature = "cuda")]
#[derive(Clone, Copy, Debug, ValueEnum, PartialEq, Eq)]
enum GpuMode {
    /// Use the GPU BWT, falling back to the CPU per-block on failure (default).
    Auto,
    /// Force the CPU BWT even though this binary was built with GPU support.
    Off,
    /// Error at startup if no usable CUDA device is present (no silent CPU fallback).
    Require,
}

/// Apply `--gpu` by setting `QZ_GPU` for this process (and, via the inherited environment,
/// any NUMA worker processes it spawns), and enforcing `require`. Must run before any
/// compression and before the NUMA driver forks workers — the call sites place it right
/// after argument parsing, while still single-threaded, so the `set_var` is sound.
#[cfg(feature = "cuda")]
fn apply_gpu_mode(mode: GpuMode) -> anyhow::Result<()> {
    match mode {
        // Respect an inherited/explicit QZ_GPU (default: GPU on).
        GpuMode::Auto => {}
        GpuMode::Off => unsafe { std::env::set_var("QZ_GPU", "off") },
        GpuMode::Require => {
            unsafe { std::env::set_var("QZ_GPU", "on") };
            if !qz_lib::compression::bsc::cuda_device_available() {
                anyhow::bail!(
                    "--gpu require: no usable CUDA device detected (cudaMemGetInfo failed). \
                     Omit --gpu require to allow the CPU-BWT fallback."
                );
            }
        }
    }
    Ok(())
}

#[derive(Parser)]
struct CompressArgs {
    /// Input FASTQ file(s). Specify once for single-end, twice for paired-end.
    #[arg(short, long, value_name = "FILE", required = true, num_args = 1.., action = clap::ArgAction::Append)]
    input: Vec<PathBuf>,

    /// Output QZ archive file
    #[arg(short, long, value_name = "FILE", required = true)]
    output: PathBuf,

    /// Working directory for temporary files
    #[arg(short, long, default_value = ".")]
    working_dir: PathBuf,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = qz_lib::cli::num_cpus())]
    threads: usize,

    /// Input is FASTA format (no quality scores)
    #[arg(long)]
    fasta: bool,

    /// Do not preserve quality values
    #[arg(long)]
    no_quality: bool,

    /// Quality compression mode
    #[arg(short, long, value_enum, default_value = "lossless")]
    quality_mode: CliQualityMode,

    /// Ultra compression with optional level (1-3, default: auto).
    /// Level 1: fast (~8 GB RAM)
    /// Level 3: high (~17 GB RAM)
    /// Auto mode selects the highest level that fits available RAM.
    #[arg(long, value_name = "LEVEL", default_missing_value = "0", num_args = 0..=1)]
    ultra: Option<u8>,

    /// JSON config file for advanced compression options (overrides defaults)
    #[arg(long, value_name = "FILE")]
    config: Option<PathBuf>,

    /// Faster compression and decompression at a small ratio cost: uses
    /// libbsc's static QLFC entropy coder for the header/sequence streams
    /// instead of the adaptive coder.
    #[arg(long)]
    fast: bool,

    /// Overwrite the output file if it already exists
    #[arg(short = 'f', long)]
    force: bool,

    /// Reference FASTA for reference-based compression (single-end or paired-end).
    /// Used only at compress time; never stored, never needed to decode. Accepts
    /// one -i input (single-end) or two (paired-end).
    #[arg(long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Prebuilt strobealign index for --reference (reserved; rejected in this release)
    #[arg(long, value_name = "INDEX", requires = "reference")]
    reference_index: Option<PathBuf>,

    /// Faster reference mapping via sparser seeds: ~5-10% faster compress at any
    /// core count, lossless, ~ratio-neutral on high-coverage data. May add a few
    /// literal fallbacks (slightly larger) on low-coverage or divergent inputs.
    #[arg(long, requires = "reference")]
    reference_fast: bool,

    /// Build the strobealign reference index inline if it's missing or stale
    /// (default: error and tell you to run `qz index` first). Only with --reference.
    #[arg(long, requires = "reference")]
    build_index: bool,

    /// Reference-mode encode pipeline depth (chunks encoded concurrently).
    /// Higher = faster compress but more memory (~window × chunk working set).
    /// 4 is the measured sweet spot with the continuous-pipeline encoder (peak
    /// speed at modest RAM); lower it to save memory, raise it for a little more
    /// speed. Output is byte-identical regardless.
    #[arg(long, value_name = "N", default_value_t = 4, requires = "reference",
          value_parser = clap::value_parser!(u32).range(1..=8))]
    reference_window: u32,

    /// NUMA multi-process compress: 'auto' (default; shards only on real NUMA
    /// hardware, else in-process), 'off', or an explicit worker count N.
    #[arg(long, value_name = "auto|off|N", default_value = "auto",
          value_parser = numa::parse_numa_mode)]
    numa: numa::NumaMode,

    /// Order-drop cluster mode: reorder reads by minimizer cluster for better sequence
    /// ratio (output is a permutation of input). Optional level: fast|balanced|max.
    #[arg(long, value_name = "LEVEL", num_args = 0..=1, default_missing_value = "balanced",
          conflicts_with_all = ["ultra", "reference"])]
    cluster: Option<String>,

    /// GPU acceleration mode (only present in `--features cuda` builds): auto = use the GPU
    /// BWT, falling back to CPU per-block (default); off = force the CPU BWT; require = error
    /// at startup if no usable CUDA device. Output is byte-identical regardless.
    #[cfg(feature = "cuda")]
    #[arg(long, value_enum, default_value_t = GpuMode::Auto)]
    gpu: GpuMode,
}

#[derive(Parser)]
struct DecompressArgs {
    /// Input QZ archive
    #[arg(short, long, value_name = "FILE", required = true)]
    input: PathBuf,

    /// Output FASTQ destination. Single-end, single-end-reference, and cluster
    /// single-end archives take ONE path. Paired-end and (non-cluster) reference
    /// archives take ONE path used as a PREFIX (-o sample → sample_R1.fastq /
    /// sample_R2.fastq). Only paired CLUSTER archives take TWO exact paths
    /// (-o R1.fastq -o R2.fastq). Passing the wrong count is rejected, not ignored.
    #[arg(short, long, value_name = "FILE", required = true, num_args = 1..)]
    output: Vec<PathBuf>,

    /// Working directory for temporary files
    #[arg(short, long, default_value = ".")]
    working_dir: PathBuf,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = qz_lib::cli::num_cpus())]
    threads: usize,

    /// Output gzipped FASTQ
    #[arg(short, long)]
    gzipped: bool,

    /// Gzip compression level (0-9)
    #[arg(long, default_value = "6", value_parser = clap::value_parser!(u32).range(0..=9))]
    gzip_level: u32,

    /// Overwrite the output file if it already exists
    #[arg(short, long)]
    force: bool,

    /// Emit a paired or paired-reference archive as a single INTERLEAVED FASTQ stream
    /// (R1/R2 records alternating) to ONE -o target (a file, '-' for stdout, or with
    /// --gzipped) — the form aligners accept directly (bwa mem -p, bowtie2 --interleaved,
    /// strobealign --interleaved). Only valid for two-mate archives; runs in-process
    /// (ignores --numa).
    #[arg(long)]
    interleaved: bool,

    /// NUMA multi-process decode: 'auto' (default; shards only on real NUMA hardware,
    /// else in-process), 'off', or an explicit worker count N.
    #[arg(long, value_name = "auto|off|N", default_value = "auto",
          value_parser = numa::parse_numa_mode)]
    numa: numa::NumaMode,

    /// GPU acceleration mode (only present in `--features cuda` builds): auto (default) | off |
    /// require. Selects the BWT backend used during decode; output is byte-identical.
    #[cfg(feature = "cuda")]
    #[arg(long, value_enum, default_value_t = GpuMode::Auto)]
    gpu: GpuMode,
}

#[derive(Parser)]
struct VerifyArgs {
    /// Input QZ archive
    #[arg(short, long, value_name = "FILE", required = true)]
    input: PathBuf,

    /// Working directory for temporary files
    #[arg(short, long, default_value = ".")]
    working_dir: PathBuf,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = qz_lib::cli::num_cpus())]
    threads: usize,

    /// Fast mode: verify per-block CRC32s only, without invoking BSC
    /// decompression. Catches bit-rot in O(IO + CRC) instead of O(BWT); does
    /// not compute the deep CRC32 over the reconstructed FASTQ bytes.
    #[arg(long)]
    fast: bool,
}

#[derive(Parser)]
struct MergeArgs {
    /// Input QZ archives to merge (≥1), in output order. Single-end, paired, or — with
    /// `--reference` — reference archives; all inputs must share one archive type + config.
    #[arg(value_name = "ARCHIVE", required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Output merged QZ archive ('-' for stdout)
    #[arg(short, long, value_name = "FILE", required = true)]
    output: PathBuf,

    /// Overwrite the output file if it already exists
    #[arg(short = 'f', long)]
    force: bool,

    /// Reference FASTA for merging REFERENCE archives (`archive_type` 2/4). Required to
    /// merge reference shards: the per-shard coverage globals are not concatenable, so
    /// they are re-derived over the union of covered intervals from this FASTA (the SAME
    /// one the shards were compressed against). Reference output is to a file, not stdout.
    #[arg(long, value_name = "FASTA")]
    reference: Option<PathBuf>,
}

#[derive(Parser)]
struct IndexArgs {
    /// Reference FASTA to index.
    #[arg(value_name = "FASTA")]
    reference: PathBuf,

    /// Read length to tune the index for (selects the seeding profile). Default 150.
    /// Any length in a profile bucket works (149/150/151 all build the 150 index).
    #[arg(short = 'r', long, value_name = "N", conflicts_with = "like")]
    read_length: Option<usize>,

    /// Infer the read length from the first read of this FASTQ instead of -r.
    #[arg(long, value_name = "FASTQ")]
    like: Option<PathBuf>,

    /// Build the sparser '--reference-fast' index variant.
    #[arg(long)]
    reference_fast: bool,

    /// Number of threads.
    #[arg(short = 't', long, default_value_t = qz_lib::cli::num_cpus())]
    threads: usize,
}

#[derive(Parser)]
struct NumaWorkerArgs {
    #[arg(long)] input: PathBuf,
    #[arg(long)] chunk_start: u32,
    #[arg(long)] chunk_end: u32,
    #[arg(long)] node: u32,
    #[arg(long)] threads: usize,
    /// Boolean SWITCH (present iff gzipped) — never `--gzipped <value>`.
    #[arg(long)] gzipped: bool,
    #[arg(long)] gzip_level: u32,
    /// This worker's index in the shard (for test-hook failure targeting).
    #[arg(long)] worker_id: u32,
    /// 1 part path (single-end) or 2 (paired/reference R1,R2).
    #[arg(long = "out-part")] out_parts: Vec<PathBuf>,
    /// Byte offset(s) in the shared output file(s) where this shard's region begins —
    /// one per `--out-part` for direct-write (single-end: 1; paired: 2, R1/R2),
    /// positionally paired with `--out-region-len`. Empty ⇒ part-file path.
    #[arg(long = "out-base-offset")] out_base_offsets: Vec<u64>,
    /// Byte length(s) of this shard's direct-write region(s), one per `--out-part`
    /// (positionally paired with `--out-base-offset`).
    #[arg(long = "out-region-len")] out_region_lens: Vec<u64>,
}

#[derive(Parser)]
struct NumaCompressWorkerArgs {
    #[arg(long)] input: PathBuf,            // R1 (or single)
    #[arg(long)] input2: Option<PathBuf>,   // R2 (paired; Task 13)
    #[arg(long)] byte_start: u64,
    #[arg(long)] byte_end: u64,
    #[arg(long)] byte_start2: Option<u64>,
    #[arg(long)] byte_end2: Option<u64>,
    #[arg(long)] node: u32,
    #[arg(long)] threads: usize,
    #[arg(long)] output: PathBuf,           // this shard's part archive
    #[arg(long)] worker_id: u32,
    #[arg(long)] plan: String,              // serialized PlanOverride (JSON)
    #[arg(long)] working_dir: PathBuf,
    // Compress flags needed to rebuild CompressConfig (single-end Phase 1):
    #[arg(long)] fasta: bool,
    #[arg(long)] no_quality: bool,
    #[arg(long)] fast: bool,
    /// Optional advanced-options JSON (serialized AdvancedOptions); empty = default.
    #[arg(long)] advanced: Option<String>,
    // Reference-mode shard (reference compress sharding): the FASTA + window/fast so the
    // worker rebuilds config.reference and compress_byte_range takes the reference arm.
    #[arg(long)] reference: Option<PathBuf>,
    #[arg(long, default_value_t = 4)] reference_window: u32,
    #[arg(long)] reference_fast: bool,
}

impl VerifyArgs {
    fn into_config(self) -> qz_lib::cli::VerifyConfig {
        qz_lib::cli::VerifyConfig {
            input: self.input,
            working_dir: self.working_dir,
            num_threads: self.threads,
            fast: self.fast,
        }
    }
}

impl CompressArgs {
    fn into_config(self) -> anyhow::Result<CompressConfig> {
        // Reject flag combinations whose intent is ambiguous. Previously these
        // were silently coerced — e.g. `--fasta --quality-mode lossless` would
        // strip quality despite the user explicitly asking for lossless.
        if self.fasta && self.quality_mode != CliQualityMode::Lossless {
            anyhow::bail!(
                "--fasta is incompatible with --quality-mode {:?} (FASTA has no quality data)",
                self.quality_mode
            );
        }
        if self.no_quality && self.quality_mode != CliQualityMode::Lossless {
            anyhow::bail!(
                "--no-quality is incompatible with --quality-mode {:?} (choose one)",
                self.quality_mode
            );
        }

        let quality_mode = match self.quality_mode {
            CliQualityMode::Lossless => LibQualityMode::Lossless,
            CliQualityMode::IlluminaBin => LibQualityMode::IlluminaBin,
            CliQualityMode::Discard => LibQualityMode::Discard,
        };

        let mut advanced = if let Some(ref config_path) = self.config {
            let json_str = std::fs::read_to_string(config_path)
                .map_err(|e| anyhow::anyhow!("Failed to read config file {config_path:?}: {e}"))?;
            serde_json::from_str(&json_str)
                .map_err(|e| anyhow::anyhow!("Failed to parse config file {config_path:?}: {e}"))?
        } else {
            qz_lib::cli::AdvancedOptions::default()
        };
        // `--fast` forces the static QLFC coder (faster encode+decode, small
        // ratio cost); a JSON config may also enable it.
        advanced.bsc_static |= self.fast;

        // The ultra level is 1-3, or 0 for auto (see --help). Reject
        // out-of-range values loudly instead of silently clamping them.
        if let Some(level) = self.ultra
            && level != 0
            && !(1..=3).contains(&level)
        {
            anyhow::bail!("--ultra level must be 1, 2, 3, or omitted for auto (got {level})");
        }

        let cluster = match self.cluster.as_deref() {
            None => None,
            Some("fast") => Some(qz_lib::cli::ClusterOptions { level: qz_lib::cli::ClusterLevel::Fast }),
            Some("balanced") => Some(qz_lib::cli::ClusterOptions { level: qz_lib::cli::ClusterLevel::Balanced }),
            Some("max") => Some(qz_lib::cli::ClusterOptions { level: qz_lib::cli::ClusterLevel::Max }),
            Some(other) => anyhow::bail!("--cluster level must be fast|balanced|max, got {other}"),
        };

        Ok(CompressConfig {
            input: self.input,
            output: self.output,
            working_dir: self.working_dir,
            threads: self.threads,
            fasta: self.fasta,
            // FASTA inputs carry no quality data; force quality off so the Auto
            // quality compressor doesn't try (and fail) to encode an empty
            // quality alphabet on large inputs.
            no_quality: self.no_quality
                || self.fasta
                || self.quality_mode == CliQualityMode::Discard,
            quality_mode,
            ultra: self.ultra,
            force: self.force,
            advanced,
            // Reference mode: require a prebuilt index unless --build-index was given.
            // (build_index has `requires = "reference"`, so it's only true with a reference.)
            // Evaluated BEFORE `reference` below so it can read `self.reference` before the move.
            require_prebuilt_index: self.reference.is_some() && !self.build_index,
            reference: self
                .reference
                .map(|reference| qz_lib::cli::ReferenceOptions {
                    reference,
                    reference_index: self.reference_index,
                    reference_fast: self.reference_fast,
                    reference_window: self.reference_window as usize,
                }),
            cluster,
        })
    }
}

impl DecompressArgs {
    fn into_config(self) -> DecompressConfig {
        DecompressConfig {
            input: self.input,
            output: self.output,
            working_dir: self.working_dir,
            num_threads: self.threads,
            gzipped: self.gzipped,
            gzip_level: self.gzip_level,
            force: self.force,
        }
    }
}

fn install_panic_hook(debug: bool) {
    std::panic::set_hook(Box::new(move |info| {
        let msg = match info.payload().downcast_ref::<&str>() {
            Some(s) => *s,
            None => match info.payload().downcast_ref::<String>() {
                Some(s) => s.as_str(),
                None => "(unknown panic payload)",
            },
        };
        let location = info
            .location()
            .map(|l| format!("{}:{}:{}", l.file(), l.line(), l.column()))
            .unwrap_or_else(|| "(unknown location)".to_string());

        eprintln!();
        eprintln!("QZ crashed with an internal error. This is a bug.");
        eprintln!("  Panic: {msg}");
        eprintln!("  At:    {location}");
        if debug {
            eprintln!();
            eprintln!("--- System diagnostics ---");
            print_system_info();
        } else {
            eprintln!();
            eprintln!("Tip: re-run with --debug for system diagnostics and a full backtrace.");
        }
        eprintln!();
        eprintln!("Please report this at: https://github.com/your-org/qz/issues");
        eprintln!("Include the above information and the command you ran.");
    }));
}

fn print_system_info() {
    eprintln!("  Version:  qz {}", env!("CARGO_PKG_VERSION"));
    eprintln!("  OS:       {}", std::env::consts::OS);
    eprintln!("  Arch:     {}", std::env::consts::ARCH);
    let cpus = std::thread::available_parallelism()
        .map(|n| n.to_string())
        .unwrap_or_else(|_| "unknown".to_string());
    eprintln!("  CPUs:     {}", cpus);
    // Read available memory from /proc/meminfo (Linux)
    if let Ok(meminfo) = std::fs::read_to_string("/proc/meminfo") {
        for line in meminfo.lines().take(3) {
            eprintln!("  {}", line);
        }
    }
    let args: Vec<String> = std::env::args().collect();
    eprintln!("  Command:  {}", args.join(" "));
}

// Worker exit codes (stable cross-process ABI — the parent in `numa::driver`
// classifies these to decide auto-fallback vs hard failure):
//   1  forced test failure (QZ_NUMA_FAIL) / any general Err propagated from run()
//   2  arg invariant violated: --out-base-offset / --out-region-len counts disagree,
//      or the region count != --out-part count (parent-side assembly bug)
//   3  NUMA pin failure (parent auto-mode falls back to in-process)
//   4  DirectWriteIntegrityError: bad/corrupt size table, region mismatch, or overrun
//      (parent auto-mode cleans the temp and falls back to in-process decode)
// Anything that newly exits with one of these codes for a DIFFERENT reason would be
// misclassified by the parent — keep this table and driver.rs in sync.
fn run_numa_worker(a: NumaWorkerArgs) -> Result<()> {
    numa::pin::worker_test_hooks(a.worker_id, "worker");
    numa::pin::pin_worker_to_node(a.node);
    // Direct-write regions: one (offset,len) per `--out-part`, paired positionally.
    // Empty ⇒ part-file path. A length mismatch is a parent-side assembly bug.
    if a.out_base_offsets.len() != a.out_region_lens.len() {
        eprintln!("numa worker: --out-base-offset / --out-region-len count mismatch");
        std::process::exit(2);
    }
    if !a.out_base_offsets.is_empty() && a.out_base_offsets.len() != a.out_parts.len() {
        eprintln!("numa worker: direct-write region count != --out-part count");
        std::process::exit(2);
    }
    let regions: Vec<qz_lib::compression::DirectWriteRegion> = a
        .out_base_offsets
        .iter()
        .zip(a.out_region_lens.iter())
        .map(|(&offset, &len)| {
            // Debug-only test hook: perturb the supplied region len so the worker's
            // base/len verification against the table fails -> integrity exit 4.
            // Shadow (not `mut`) so release builds — where the hook compiles out —
            // don't trip `unused_mut`.
            #[cfg(debug_assertions)]
            let len = if std::env::var_os("QZ_NUMA_FAKE_REGION").is_some() {
                len + 1
            } else {
                len
            };
            qz_lib::compression::DirectWriteRegion { offset, len }
        })
        .collect();
    match qz_lib::compression::decode_chunk_range(
        &a.input,
        a.chunk_start,
        a.chunk_end,
        a.threads,
        a.gzipped,
        a.gzip_level,
        &a.out_parts,
        &regions,
    ) {
        Ok(()) => Ok(()),
        Err(e) => {
            // A direct-write integrity failure (bad table/region/overrun) is
            // distinguishable so the parent's `auto` mode can clean up and fall
            // back to in-process decode. Any other error is a genuine failure.
            if e.downcast_ref::<qz_lib::compression::DirectWriteIntegrityError>().is_some() {
                eprintln!("numa worker: {e}");
                std::process::exit(4);
            }
            Err(e)
        }
    }
}

// Compress-worker exit codes (parent classifies): 0 ok; 1 general Err / forced test
// failure; 3 NUMA pin failure (parent auto-mode falls back to in-process compress).
fn run_numa_compress_worker(a: NumaCompressWorkerArgs) -> Result<()> {
    numa::pin::worker_test_hooks(a.worker_id, "compress worker");
    numa::pin::pin_worker_to_node(a.node);

    let advanced: qz_lib::cli::AdvancedOptions = match a.advanced.as_deref() {
        Some(s) if !s.is_empty() => serde_json::from_str(s)
            .map_err(|e| anyhow::anyhow!("worker: bad --advanced JSON: {e}"))?,
        _ => qz_lib::cli::AdvancedOptions::default(),
    };
    let plan: qz_lib::compression::PlanOverride = serde_json::from_str(&a.plan)
        .map_err(|e| anyhow::anyhow!("worker: bad --plan JSON: {e}"))?;

    let mut inputs = vec![a.input.clone()];
    let mut ranges = vec![qz_lib::compression::ByteRange { start: a.byte_start, end: a.byte_end }];
    if let Some(input2) = a.input2.clone() {
        inputs.push(input2);
        let s2 = a.byte_start2.ok_or_else(|| anyhow::anyhow!(
            "worker: --input2 given but --byte-start2 missing"))?;
        let e2 = a.byte_end2.ok_or_else(|| anyhow::anyhow!(
            "worker: --input2 given but --byte-end2 missing"))?;
        ranges.push(qz_lib::compression::ByteRange { start: s2, end: e2 });
    }
    let config = CompressConfig {
        input: inputs,
        output: a.output.clone(),
        working_dir: a.working_dir.clone(),
        threads: a.threads,
        fasta: a.fasta,
        no_quality: a.no_quality,
        // non-lossless quality modes gate in-process via the driver precheck (Task 7)
        quality_mode: LibQualityMode::Lossless,
        ultra: None,
        force: true, // part file lives in a private TempDir
        advanced: { let mut adv = advanced; adv.bsc_static |= a.fast; adv },
        // Workers are load-only: the parent driver ensured the index before spawning,
        // so a worker must never silently rebuild — require it when in reference mode.
        require_prebuilt_index: a.reference.is_some(),
        reference: a.reference.clone().map(|reference| qz_lib::cli::ReferenceOptions {
            reference,
            reference_index: None,
            reference_fast: a.reference_fast,
            reference_window: a.reference_window as usize,
        }),
        cluster: None,
    };
    qz_lib::compression::compress_byte_range(&config, &ranges, &plan)
}

fn run(cli: Cli) -> Result<()> {
    if std::env::var("QZ_NO_BANNER").is_err() {
        eprintln!(
            "QZ v{} - Columnar FASTQ compression",
            env!("CARGO_PKG_VERSION")
        );
        eprintln!("Compression: BSC/BWT-based columnar encoding");
        eprintln!();
    }

    if cli.debug {
        eprintln!("--- Debug mode ---");
        print_system_info();
        eprintln!();
    }

    match cli.command {
        Commands::Compress(args) => {
            info!("Starting compression...");
            let mode = args.numa; // not copied into CompressConfig
            let build_index = args.build_index; // CLI-only policy; not in CompressConfig
            #[cfg(feature = "cuda")]
            apply_gpu_mode(args.gpu)?; // CLI-only; sets QZ_GPU before the NUMA driver forks
            let config = args.into_config()?;
            numa::compress_driver::compress_with_numa(mode, &config, build_index)?;
            info!("Compression complete!");
        }
        Commands::Decompress(args) => {
            info!("Starting decompression...");
            let mode = args.numa;
            let interleaved = args.interleaved; // CLI-only; not copied into the config
            #[cfg(feature = "cuda")]
            apply_gpu_mode(args.gpu)?; // CLI-only; sets QZ_GPU before the NUMA driver forks
            let config = args.into_config(); // numa/interleaved are NOT copied into the config
            if interleaved {
                // Interleaved merges both mates into one stream → in-process single-sink
                // decode (NUMA shards into per-mate regions, which interleaving can't use).
                qz_lib::compression::decompress_interleaved(&config)?;
            } else {
                numa::driver::decompress_with_numa(mode, &config)?;
            }
            info!("Decompression complete!");
        }
        Commands::NumaDecodeWorker(a) => {
            run_numa_worker(a)?;
        }
        Commands::NumaCompressWorker(a) => {
            run_numa_compress_worker(a)?;
        }
        Commands::Verify(args) => {
            info!("Verifying archive...");
            let input_display = args.input.display().to_string();
            let config = args.into_config();
            let result = qz_lib::compression::verify(&config)?;

            // Encoding name is derived from the library's single source of truth
            // so it can't drift out of sync when new encoding types are added.
            let encoding_name = format!(
                "{} ({})",
                qz_lib::compression::encoding_type_display_name(result.encoding_type),
                result.encoding_type,
            );

            eprintln!("Archive:     {}", input_display);
            // v5 fast verify walks every block frame in every stream, so it always
            // covers the whole archive — there is no partial-coverage case.
            eprintln!("Status:      OK");
            eprintln!(
                "Mode:        {}",
                match result.mode {
                    qz_lib::compression::VerifyMode::Deep => "deep (full decompress + FASTQ CRC32)",
                    qz_lib::compression::VerifyMode::Fast => "fast (per-block CRC32 only)",
                }
            );
            eprintln!("Reads:       {}", result.num_reads);
            eprintln!("Encoding:    {}", encoding_name);
            eprintln!(
                "Headers:     {} bytes ({:?})",
                result.headers_compressed_len, result.header_compressor
            );
            eprintln!("Sequences:   {} bytes", result.sequences_compressed_len);
            eprintln!(
                "Qualities:   {} bytes ({:?})",
                result.qualities_compressed_len, result.quality_compressor
            );
            match result.mode {
                qz_lib::compression::VerifyMode::Deep => {
                    match result.r2_crc32 {
                        Some(r2) => {
                            // Two-mate (reference) deep verify reports a per-mate CRC.
                            eprintln!(
                                "R1 CRC32:    {:08x} (over reconstructed R1 FASTQ bytes)",
                                result.crc32
                            );
                            eprintln!(
                                "R2 CRC32:    {:08x} (over reconstructed R2 FASTQ bytes)",
                                r2
                            );
                        }
                        None => {
                            eprintln!(
                                "CRC32:       {:08x} (over reconstructed FASTQ bytes)",
                                result.crc32
                            );
                        }
                    }
                    eprintln!("FASTQ size:  {} bytes", result.total_bytes);
                }
                qz_lib::compression::VerifyMode::Fast => {
                    eprintln!(
                        "Blocks:      {} verified (per-block CRC32 only)",
                        result.blocks_verified
                    );
                    eprintln!("Comp size:   {} bytes", result.total_bytes);
                }
            }
            eprintln!("Verified in: {:.2}s", result.elapsed_secs);

            info!("Verification complete!");
        }
        Commands::Index(args) => {
            let read_len = match (args.read_length, &args.like) {
                (Some(r), _) => r,
                (None, Some(fq)) => qz_lib::compression::peek_first_read_length(fq)?,
                (None, None) => {
                    eprintln!(
                        "qz index: no -r/--read-length or --like given; using the 150 bp profile"
                    );
                    150
                }
            };
            let stats = qz_lib::compression::build_reference_index(
                &args.reference,
                read_len,
                args.reference_fast,
                args.threads,
            )?;
            eprintln!(
                "qz index: wrote {} ({} reference(s), {} randstrobes, {:.2}s)",
                stats.path.display(),
                stats.references,
                stats.randstrobes,
                stats.elapsed_secs
            );
        }
        Commands::Merge(args) => {
            use std::io::Write;
            info!("Merging {} archive(s)...", args.inputs.len());
            if let Some(reference) = &args.reference {
                // Reference merge re-derives the coverage globals from the FASTA, so it
                // needs a seekable file output (not stdout).
                if args.output.as_os_str() == "-" {
                    anyhow::bail!("qz merge --reference cannot stream to stdout (re-derives globals)");
                }
                qz_lib::compression::merge_reference_archives_to_path(
                    &args.inputs,
                    reference,
                    &args.output,
                    args.force,
                )?;
            } else if args.output.as_os_str() == "-" {
                let stdout = std::io::stdout();
                let mut out = std::io::BufWriter::new(stdout.lock());
                qz_lib::compression::merge_archives(&args.inputs, &mut out)?;
                out.flush()?;
            } else {
                qz_lib::compression::merge_archives_to_path(
                    &args.inputs,
                    &args.output,
                    args.force,
                )?;
            }
            info!("Merge complete!");
        }
    }

    Ok(())
}

fn main() {
    let cli = Cli::parse();
    let debug = cli.debug;

    // Enable full backtraces when --debug is set (anyhow captures these at error creation)
    if debug {
        // SAFETY: only called at startup before any threads are spawned
        unsafe { std::env::set_var("RUST_BACKTRACE", "full") };
    }

    let log_filter = if debug {
        tracing_subscriber::EnvFilter::new("debug")
    } else {
        tracing_subscriber::EnvFilter::try_from_default_env()
            .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info"))
    };

    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .with_env_filter(log_filter)
        .init();

    install_panic_hook(debug);

    if let Err(e) = run(cli) {
        eprintln!();
        eprintln!("Error: {:#}", e);
        if debug {
            eprintln!();
            eprintln!("--- System diagnostics ---");
            print_system_info();
        } else {
            eprintln!();
            eprintln!("Tip: re-run with --debug for verbose logging and system diagnostics.");
        }
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    /// `--reference-fast` (and `--reference-index`) carry `requires = "reference"`,
    /// so Clap must REJECT them when `--reference` is absent, with a message that
    /// names `reference`.
    #[test]
    fn reference_fast_requires_reference() {
        let res = Cli::try_parse_from([
            "qz",
            "compress",
            "-i",
            "r1.fastq",
            "-i",
            "r2.fastq",
            "-o",
            "out.qz",
            "--reference-fast",
        ]);
        let err = res
            .err()
            .expect("--reference-fast without --reference must error");
        let msg = err.to_string();
        assert!(
            msg.contains("reference"),
            "error should mention reference, got: {msg}"
        );
    }

    #[test]
    fn reference_index_requires_reference() {
        let res = Cli::try_parse_from([
            "qz",
            "compress",
            "-i",
            "r1.fastq",
            "-i",
            "r2.fastq",
            "-o",
            "out.qz",
            "--reference-index",
            "idx",
        ]);
        let err = res
            .err()
            .expect("--reference-index without --reference must error");
        assert!(err.to_string().contains("reference"));
    }

    /// With `--reference` present, the reserved flags parse (they are rejected
    /// later by the library gate, not by Clap).
    #[test]
    fn reference_flags_parse_with_reference() {
        Cli::try_parse_from([
            "qz",
            "compress",
            "-i",
            "r1.fastq",
            "-i",
            "r2.fastq",
            "-o",
            "out.qz",
            "--reference",
            "ref.fa",
            "--reference-fast",
        ])
        .expect("flags should parse when --reference is given");
    }

    #[test]
    fn index_r_and_like_conflict() {
        let res = Cli::try_parse_from(["qz", "index", "ref.fa", "-r", "150", "--like", "reads.fq"]);
        assert!(res.is_err(), "-r and --like must be mutually exclusive");
    }

    #[test]
    fn index_minimal_parses() {
        let res = Cli::try_parse_from(["qz", "index", "ref.fa"]);
        assert!(res.is_ok(), "bare `qz index ref.fa` must parse");
    }

    #[test]
    fn numa_compress_worker_subcommand_parses() {
        let cli = Cli::try_parse_from([
            "qz", "_numa-compress-worker",
            "--input", "in.fastq",
            "--byte-start", "0", "--byte-end", "100",
            "--node", "1", "--threads", "8",
            "--output", "w0.qz", "--worker-id", "0",
            "--plan", "{}",
            "--working-dir", ".",
        ]).expect("worker subcommand must parse");
        match cli.command {
            Commands::NumaCompressWorker(a) => {
                assert_eq!(a.byte_start, 0);
                assert_eq!(a.byte_end, 100);
                assert_eq!(a.node, 1);
                assert_eq!(a.worker_id, 0);
            }
            _ => panic!("expected NumaCompressWorker"),
        }
    }
}
