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
}

#[derive(Subcommand)]
enum Commands {
    /// Compress FASTQ files
    Compress(CompressArgs),
    /// Decompress QZ archives
    Decompress(DecompressArgs),
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

    /// Number of threads (0 = auto-detect)
    #[arg(short = 't', long, default_value = "0")]
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

impl CompressArgs {
    fn into_config(self) -> CompressConfig {
        let quality_mode = match self.quality_mode {
            CliQualityMode::Lossless => LibQualityMode::Lossless,
            CliQualityMode::IlluminaBin => LibQualityMode::IlluminaBin,
            CliQualityMode::Discard => LibQualityMode::Discard,
        };
        CompressConfig {
            input: vec![self.input],
            output: self.output,
            working_dir: self.working_dir,
            threads: self.threads,
            fasta: self.fasta,
            no_quality: self.no_quality || self.quality_mode == CliQualityMode::Discard,
            quality_mode,
            ultra: self.ultra,
            ..CompressConfig::default()
        }
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

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let cli = Cli::parse();

    if std::env::var("QZ_NO_BANNER").is_err() {
        eprintln!("QZ v{} - Columnar FASTQ compression", env!("CARGO_PKG_VERSION"));
        eprintln!("Compression: BSC/BWT-based columnar encoding");
        eprintln!();
    }

    match cli.command {
        Commands::Compress(args) => {
            info!("Starting compression...");
            let config = args.into_config();
            qz_lib::compression::compress(&config)?;
            info!("Compression complete!");
        }
        Commands::Decompress(args) => {
            info!("Starting decompression...");
            let config = args.into_config();
            qz_lib::compression::decompress(&config)?;
            info!("Decompression complete!");
        }
    }

    Ok(())
}
