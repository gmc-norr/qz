use anyhow::Result;
use clap::{Parser, Subcommand};
#[cfg(feature = "cuda")]
use clap::ValueEnum;
use std::path::PathBuf;
use tracing::info;

mod numa;

// Pooling global allocator. Compress allocates a fresh large output buffer per
// BSC block plus a libsais workspace, then frees it; under glibc malloc these
// large allocations are mmap'd and munmap'd back to the kernel every block, so
// each one re-faults and re-zeros its pages. mimalloc recycles freed buffers in
// userspace, eliminating the recurring faults: ~8% faster compress with *lower*
// peak RSS at scale (HG002 chr20, 8-chunk pipeline), byte-identical output, no
// ratio change. (Mirrors qz-cli; see its note for the full rationale.)
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Parser)]
#[command(name = "bz")]
#[command(author = "BZ Contributors")]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(about = "High-performance BAM compression using columnar encoding and BSC/BWT", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Enable debug mode: verbose logging, full backtraces, and system diagnostics on crash
    #[arg(long, global = true)]
    debug: bool,
}

#[derive(Subcommand)]
enum Commands {
    /// Compress BAM files
    Compress(CompressArgs),
    /// Decompress BZ archives to BAM
    Decompress(DecompressArgs),
    /// Extract FASTQ from BAM and compress to QZ format
    Extract(ExtractArgs),
    /// Verify archive integrity without decompressing to disk
    Verify(VerifyArgs),
    /// Internal: one NUMA decode worker (hidden — not a stable interface).
    #[command(name = "_numa-decode-worker", hide = true)]
    NumaDecodeWorker(numa::NumaWorkerArgs),
    /// Internal: one NUMA compress worker (hidden — not a stable interface).
    #[command(name = "_numa-compress-worker", hide = true)]
    NumaCompressWorker(numa::NumaCompressWorkerArgs),
}

/// GPU acceleration mode (cuda builds only — the flag is absent from CPU builds).
#[cfg(feature = "cuda")]
#[derive(Clone, Copy, Debug, ValueEnum, PartialEq, Eq)]
enum GpuMode {
    /// Use the GPU BWT, falling back to the CPU per-block on failure.
    Auto,
    /// Force the CPU BWT (bz default — GPU regresses bz).
    Off,
    /// Error at startup if no usable CUDA device is present.
    Require,
}

/// Apply `--gpu` by setting `QZ_GPU` for this process (inherited by any NUMA worker
/// processes), and enforcing `require`. MUST run while still single-threaded — before
/// `build_global()` and before the NUMA driver forks — so the `set_var` is sound.
#[cfg(feature = "cuda")]
fn apply_gpu_mode(mode: GpuMode) -> Result<()> {
    match mode {
        // Respect an inherited/explicit QZ_GPU (process default: GPU on).
        GpuMode::Auto => {}
        GpuMode::Off => unsafe { std::env::set_var("QZ_GPU", "off") },
        GpuMode::Require => {
            unsafe { std::env::set_var("QZ_GPU", "on") };
            if !bz_lib::cuda_device_available() {
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
    /// Input BAM file
    #[arg(short, long, value_name = "FILE", required = true)]
    input: PathBuf,
    /// Output BZ archive file
    #[arg(short, long, value_name = "FILE", required = true)]
    output: PathBuf,
    /// Working directory for temporary files
    #[arg(short, long, default_value = ".")]
    working_dir: PathBuf,
    /// Number of threads (0 = auto-detect)
    #[arg(short = 't', long, default_value = "0")]
    threads: usize,
    /// Advanced options JSON config file
    #[arg(long, value_name = "FILE")]
    config: Option<PathBuf>,
    /// Compression level 0-3 (lossless) — a memory/speed preset; ratio is ~flat
    /// across levels. 0 = auto (default: balanced, downshifts under low RAM);
    /// 1 = lowest memory, 2 = balanced/fastest, 3 = largest chunks (most memory).
    #[arg(short = 'l', long, value_name = "LEVEL", value_parser = clap::value_parser!(u8).range(0..=3))]
    level: Option<u8>,
    /// LOSSY: variant-aware quality reduction level 1-3 (higher = more
    /// flattening). Flattens quality at confident reference-matching positions
    /// while keeping full resolution at candidate-variant sites. Off by default.
    #[arg(long, value_name = "LEVEL", value_parser = clap::value_parser!(u8).range(1..=3))]
    reduce_quality: Option<u8>,
    /// Flatten scheme for --reduce-quality: "twobin" (default), "coarse8", or
    /// "single" (one fixed Q40).
    #[arg(long, value_name = "SCHEME", default_value = "twobin")]
    quality_scheme: String,
    /// NUMA multi-process compress: 'auto' (default; shards only on real NUMA
    /// hardware, ≥2 sockets, above a size threshold), 'off', or N worker processes.
    /// Each worker compresses a disjoint chunk range; the archive is byte-identical
    /// to a single-process compress.
    #[arg(long, default_value = "auto", value_name = "MODE",
          value_parser = numa_core::parse_numa_mode)]
    numa: numa_core::NumaMode,
    /// GPU acceleration (cuda builds only): auto | off (bz default) | require. GPU regresses
    /// bz (many small BSC blocks serialize through one GPU), so bz defaults to off; auto or
    /// require opts back in. Output is byte-identical.
    #[cfg(feature = "cuda")]
    #[arg(long, value_enum, default_value_t = GpuMode::Off)]
    gpu: GpuMode,
    /// Overwrite the output file if it already exists
    #[arg(short = 'f', long)]
    force: bool,
}

#[derive(Parser)]
struct DecompressArgs {
    /// Input BZ archive
    #[arg(short, long, value_name = "FILE", required = true)]
    input: PathBuf,
    /// Output BAM file
    #[arg(short, long, value_name = "FILE", required = true)]
    output: PathBuf,
    /// Working directory for temporary files
    #[arg(short, long, default_value = ".")]
    working_dir: PathBuf,
    /// Number of threads (0 = auto-detect)
    #[arg(short = 't', long, default_value = "0")]
    threads: usize,
    /// NUMA multi-process decode: 'auto' (default; shards only on real NUMA
    /// hardware, ≥2 sockets, above a size threshold), 'off', or N worker processes.
    #[arg(long, default_value = "auto", value_name = "MODE",
          value_parser = numa_core::parse_numa_mode)]
    numa: numa_core::NumaMode,
    /// GPU acceleration (cuda builds only): auto | off (bz default) | require. Selects the
    /// BWT backend used during decode; output is byte-identical.
    #[cfg(feature = "cuda")]
    #[arg(long, value_enum, default_value_t = GpuMode::Off)]
    gpu: GpuMode,
    /// Overwrite the output file if it already exists
    #[arg(short = 'f', long)]
    force: bool,
}

#[derive(Parser)]
struct ExtractArgs {
    /// Input BAM file
    #[arg(short, long, value_name = "FILE", required = true)]
    input: PathBuf,
    /// Output prefix (creates {prefix}_R1.qz and {prefix}_R2.qz for paired-end, or {prefix}_SE.qz for single-end)
    #[arg(short, long, value_name = "PREFIX", required = true)]
    output: String,
    /// Working directory for temporary files
    #[arg(short, long, default_value = ".")]
    working_dir: PathBuf,
    /// Number of threads (0 = auto-detect)
    #[arg(short = 't', long, default_value = "0")]
    threads: usize,
    /// GPU acceleration (cuda builds only): auto | off (bz default) | require. Output is
    /// byte-identical.
    #[cfg(feature = "cuda")]
    #[arg(long, value_enum, default_value_t = GpuMode::Off)]
    gpu: GpuMode,
    /// Overwrite output files if they already exist
    #[arg(short = 'f', long)]
    force: bool,
}

#[derive(Parser)]
struct VerifyArgs {
    /// Input BZ archive
    #[arg(short, long, value_name = "FILE", required = true)]
    input: PathBuf,
    /// Number of threads (0 = auto-detect)
    #[arg(short = 't', long, default_value = "0")]
    threads: usize,
}

/// Refuse to clobber an existing output file unless `force` is set.
fn check_overwrite(path: &std::path::Path, force: bool) -> Result<()> {
    if path.exists() && !force {
        anyhow::bail!(
            "Output file already exists: {}\nUse --force to overwrite",
            path.display()
        );
    }
    Ok(())
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
        eprintln!("BZ crashed with an internal error. This is a bug.");
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
    eprintln!("  Version:  bz {}", env!("CARGO_PKG_VERSION"));
    eprintln!("  OS:       {}", std::env::consts::OS);
    eprintln!("  Arch:     {}", std::env::consts::ARCH);
    let cpus = std::thread::available_parallelism()
        .map(|n| n.to_string())
        .unwrap_or_else(|_| "unknown".to_string());
    eprintln!("  CPUs:     {}", cpus);
    if let Ok(meminfo) = std::fs::read_to_string("/proc/meminfo") {
        for line in meminfo.lines().take(3) {
            eprintln!("  {}", line);
        }
    }
    let args: Vec<String> = std::env::args().collect();
    eprintln!("  Command:  {}", args.join(" "));
}

fn run(cli: Cli) -> Result<()> {
    // A NUMA decode worker is a hidden re-exec subprocess: pin to its node and build
    // its own (node-local) rayon pool, BEFORE the generic banner + global pool init
    // below. Handle it first and return.
    if let Commands::NumaDecodeWorker(args) = &cli.command {
        return numa::run_decode_worker(args);
    }
    if let Commands::NumaCompressWorker(args) = &cli.command {
        return numa::run_compress_worker(args);
    }

    eprintln!(
        "BZ v{} - Columnar BAM compression with consensus-delta encoding",
        env!("CARGO_PKG_VERSION")
    );
    eprintln!();

    if cli.debug {
        eprintln!("--- Debug mode ---");
        print_system_info();
        eprintln!();
    }

    // Initialise the rayon global thread pool once here, in the binary.
    // Library code must not call build_global() — that mutates global state
    // and silently fails when called more than once, which is wrong for a lib.
    let threads: usize = match &cli.command {
        Commands::Compress(a) => a.threads,
        Commands::Decompress(a) => a.threads,
        Commands::Extract(a) => a.threads,
        Commands::Verify(a) => a.threads,
        // Handled above (their own pinned pools); never reach here.
        Commands::NumaDecodeWorker(_) => 0,
        Commands::NumaCompressWorker(_) => 0,
    };
    // Apply --gpu (sets QZ_GPU, inherited by NUMA workers) while still single-threaded,
    // BEFORE build_global() spawns the rayon pool.
    #[cfg(feature = "cuda")]
    {
        // bz defaults GPU off everywhere. The NUMA worker subcommands already returned
        // above (they inherit QZ_GPU); the `_` arm here is `verify`, which has no --gpu flag
        // and so stays off like the rest of bz.
        let gpu = match &cli.command {
            Commands::Compress(a) => a.gpu,
            Commands::Decompress(a) => a.gpu,
            Commands::Extract(a) => a.gpu,
            _ => GpuMode::Off,
        };
        apply_gpu_mode(gpu)?;
    }
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .map_err(|e| anyhow::anyhow!("Failed to initialise thread pool: {e}"))?;
    }

    match cli.command {
        Commands::Compress(args) => {
            check_overwrite(&args.output, args.force)?;
            info!("Starting compression...");
            let advanced = if let Some(ref config_path) = args.config {
                let json_str = std::fs::read_to_string(config_path).map_err(|e| {
                    anyhow::anyhow!("Failed to read config file {config_path:?}: {e}")
                })?;
                let opts: bz_lib::AdvancedOptions =
                    serde_json::from_str(&json_str).map_err(|e| {
                        anyhow::anyhow!("Failed to parse config file {config_path:?}: {e}")
                    })?;
                opts.validate()?;
                opts
            } else {
                bz_lib::AdvancedOptions::default()
            };
            // CLI --level / --reduce-quality override any config value.
            let mut advanced = advanced;
            if let Some(level) = args.level {
                advanced.level = level;
            }
            if let Some(level) = args.reduce_quality {
                let scheme = match args.quality_scheme.as_str() {
                    "twobin" => bz_lib::FlattenScheme::TwoBin {
                        thresh: 20,
                        low: 6,
                        high: 40,
                    },
                    "coarse8" => bz_lib::FlattenScheme::Coarse8,
                    "single" => bz_lib::FlattenScheme::Single(40),
                    other => anyhow::bail!(
                        "unknown --quality-scheme {other:?} (expected twobin|coarse8|single)"
                    ),
                };
                advanced.quality_reduction = Some(bz_lib::QualityReduction::level(level, scheme));
            }
            let numa_mode = args.numa;
            let force = args.force;
            // Resolve the thread budget for the NUMA gate (the global rayon pool was
            // already sized from `args.threads` above; 0 = all detected cores).
            let total_threads = if args.threads > 0 {
                args.threads
            } else {
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(1)
            };
            let config = bz_lib::CompressConfig {
                input: args.input,
                output: args.output,
                working_dir: args.working_dir,
                advanced,
            };
            numa::compress_with_numa(numa_mode, &config, total_threads, force)?;
            info!("Compression complete!");
        }
        Commands::Decompress(args) => {
            check_overwrite(&args.output, args.force)?;
            info!("Starting decompression...");
            let numa_mode = args.numa;
            let force = args.force;
            // Resolve the thread budget for the NUMA gate (the global rayon pool was
            // already sized from `args.threads` above; 0 = all detected cores).
            let total_threads = if args.threads > 0 {
                args.threads
            } else {
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(1)
            };
            let config = bz_lib::DecompressConfig {
                input: args.input,
                output: args.output,
                working_dir: args.working_dir,
            };
            numa::decompress_with_numa(numa_mode, &config, total_threads, force)?;
            info!("Decompression complete!");
        }
        Commands::Extract(args) => {
            // Extract writes up to three outputs derived from --output prefix.
            // Refuse to clobber any of them unless --force is set.
            for suffix in ["_R1.qz", "_R2.qz", "_SE.qz"] {
                let path = PathBuf::from(format!("{}{}", args.output, suffix));
                check_overwrite(&path, args.force)?;
            }
            info!("Starting BAM to QZ extraction...");
            let config = bz_lib::ExtractConfig {
                input: args.input,
                output_prefix: args.output,
                working_dir: args.working_dir,
                force: args.force,
            };
            bz_lib::extract(&config)?;
            info!("Extraction complete!");
        }
        Commands::Verify(args) => {
            let config = bz_lib::VerifyConfig {
                input: args.input.clone(),
            };
            match bz_lib::verify(&config) {
                Ok(result) => {
                    eprintln!("Archive:     {:?}", args.input);
                    eprintln!("Status:      OK");
                    eprintln!("Records:     {}", result.num_records);
                    eprintln!("Chunks:      {}", result.num_chunks);
                    eprintln!("CRC32:       {:08x}", result.crc32);
                    eprintln!("Data size:   {} bytes", result.total_bytes);
                    eprintln!("Verified in: {:.2}s", result.elapsed_secs);
                }
                Err(e) => {
                    eprintln!("Archive:     {:?}", args.input);
                    eprintln!("Status:      FAILED");
                    eprintln!("Error:       {:#}", e);
                    std::process::exit(1);
                }
            }
        }
        // Handled at the top of `run` (before the banner / global pool init).
        Commands::NumaDecodeWorker(_) => unreachable!("numa worker dispatched before run body"),
        Commands::NumaCompressWorker(_) => {
            unreachable!("numa compress worker dispatched before run body")
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
