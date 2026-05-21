use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;
use tracing::info;

use qz_lib::cli::{
    CompressConfig, DecompressConfig, QualityMode as LibQualityMode,
};

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

#[derive(Parser)]
struct CompressArgs {
    /// Input FASTQ file
    #[arg(short, long, value_name = "FILE", required = true)]
    input: PathBuf,

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

    /// Ultra compression with optional level (1-5, default: auto).
    /// Level 1: fast (~8 GB RAM)
    /// Level 3: high (~17 GB RAM)
    /// Level 5: extreme (~17 GB RAM)
    /// Auto mode selects the highest level that fits available RAM.
    #[arg(long, value_name = "LEVEL", default_missing_value = "0", num_args = 0..=1)]
    ultra: Option<u8>,

    /// JSON config file for advanced compression options (overrides defaults)
    #[arg(long, value_name = "FILE")]
    config: Option<PathBuf>,

    /// Overwrite the output file if it already exists
    #[arg(short = 'f', long)]
    force: bool,
}

#[derive(Parser)]
struct DecompressArgs {
    /// Input QZ archive
    #[arg(short, long, value_name = "FILE", required = true)]
    input: PathBuf,

    /// Output FASTQ file
    #[arg(short, long, value_name = "FILE", required = true)]
    output: PathBuf,

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
    #[arg(long, default_value = "6")]
    gzip_level: u32,
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

        let advanced = if let Some(ref config_path) = self.config {
            let json_str = std::fs::read_to_string(config_path)
                .map_err(|e| anyhow::anyhow!("Failed to read config file {config_path:?}: {e}"))?;
            serde_json::from_str(&json_str)
                .map_err(|e| anyhow::anyhow!("Failed to parse config file {config_path:?}: {e}"))?
        } else {
            qz_lib::cli::AdvancedOptions::default()
        };

        Ok(CompressConfig {
            input: vec![self.input],
            output: self.output,
            working_dir: self.working_dir,
            threads: self.threads,
            fasta: self.fasta,
            no_quality: self.no_quality || self.quality_mode == CliQualityMode::Discard,
            quality_mode,
            ultra: self.ultra,
            force: self.force,
            advanced,
        })
    }
}

impl DecompressArgs {
    fn into_config(self) -> DecompressConfig {
        DecompressConfig {
            input: self.input,
            output: vec![self.output],
            working_dir: self.working_dir,
            num_threads: self.threads,
            gzipped: self.gzipped,
            gzip_level: self.gzip_level,
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

fn run(cli: Cli) -> Result<()> {
    if std::env::var("QZ_NO_BANNER").is_err() {
        eprintln!("QZ v{} - Columnar FASTQ compression", env!("CARGO_PKG_VERSION"));
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
            let config = args.into_config()?;
            qz_lib::compression::compress(&config)?;
            info!("Compression complete!");
        }
        Commands::Decompress(args) => {
            info!("Starting decompression...");
            let config = args.into_config();
            qz_lib::compression::decompress(&config)?;
            info!("Decompression complete!");
        }
        Commands::Verify(args) => {
            info!("Verifying archive...");
            let input_display = args.input.display().to_string();
            let config = args.into_config();
            let result = qz_lib::compression::verify(&config)?;

            // String, not &str, so the borrow lives until after the eprintln chain.
            let encoding_name: String = match result.encoding_type {
                0 => "default (0)".to_string(),
                4 => "raw+hints (4)".to_string(),
                6 => "rc-canon (6)".to_string(),
                8 => "local-reorder (8)".to_string(),
                9 => "ultra (9)".to_string(),
                n => format!("unknown ({n})"),
            };

            eprintln!("Archive:     {}", input_display);
            // Distinguish "OK" (all data checked) from "PARTIAL" (fast verify
            // skipped one or more non-v3 streams). A user grepping for "OK"
            // shouldn't mistake partial coverage for full integrity.
            let status = if result.streams_skipped > 0 {
                format!("PARTIAL ({} stream(s) skipped)", result.streams_skipped)
            } else {
                "OK".to_string()
            };
            eprintln!("Status:      {}", status);
            eprintln!("Mode:        {}", match result.mode {
                qz_lib::compression::VerifyMode::Deep => "deep (full decompress + FASTQ CRC32)",
                qz_lib::compression::VerifyMode::Fast => "fast (per-block CRC32 only)",
            });
            eprintln!("Reads:       {}", result.num_reads);
            eprintln!("Encoding:    {}", encoding_name);
            eprintln!("Headers:     {} bytes ({:?})", result.headers_compressed_len, result.header_compressor);
            eprintln!("Sequences:   {} bytes", result.sequences_compressed_len);
            eprintln!("Qualities:   {} bytes ({:?})", result.qualities_compressed_len, result.quality_compressor);
            match result.mode {
                qz_lib::compression::VerifyMode::Deep => {
                    eprintln!("CRC32:       {:08x} (over reconstructed FASTQ bytes)", result.crc32);
                    eprintln!("FASTQ size:  {} bytes", result.total_bytes);
                }
                qz_lib::compression::VerifyMode::Fast => {
                    eprintln!("Blocks:      {} verified (per-block CRC32 only)", result.blocks_verified);
                    eprintln!("Comp size:   {} bytes", result.total_bytes);
                    if result.streams_skipped > 0 {
                        eprintln!(
                            "Skipped:     {} stream(s) not v3-framed (e.g. raw zstd) \
                             — rerun without --fast for full coverage",
                            result.streams_skipped
                        );
                    }
                }
            }
            eprintln!("Verified in: {:.2}s", result.elapsed_secs);

            info!("Verification complete!");
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
