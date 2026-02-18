use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;
use tracing::info;

#[derive(Parser)]
#[command(name = "bz")]
#[command(author = "BZ Contributors")]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(about = "High-performance BAM compression using columnar encoding and BSC/BWT", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compress BAM files
    Compress(CompressArgs),
    /// Decompress BZ archives to BAM
    Decompress(DecompressArgs),
    /// Extract FASTQ from BAM and compress to QZ format
    Extract(ExtractArgs),
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
}

#[derive(Parser)]
struct ExtractArgs {
    /// Input BAM file
    #[arg(short, long, value_name = "FILE", required = true)]
    input: PathBuf,
    /// Output prefix (creates {prefix}_R1.qz and {prefix}_R2.qz for paired-end, or {prefix}.qz for single-end)
    #[arg(short, long, value_name = "PREFIX", required = true)]
    output: String,
    /// Working directory for temporary files
    #[arg(short, long, default_value = ".")]
    working_dir: PathBuf,
    /// Number of threads (0 = auto-detect)
    #[arg(short = 't', long, default_value = "0")]
    threads: usize,
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

    eprintln!(
        "BZ v{} - Columnar BAM compression with consensus-delta encoding",
        env!("CARGO_PKG_VERSION")
    );
    eprintln!();

    match cli.command {
        Commands::Compress(args) => {
            info!("Starting compression...");
            let config = bz_lib::CompressConfig {
                input: args.input,
                output: args.output,
                working_dir: args.working_dir,
                threads: args.threads,
            };
            bz_lib::compress(&config)?;
            info!("Compression complete!");
        }
        Commands::Decompress(args) => {
            info!("Starting decompression...");
            let config = bz_lib::DecompressConfig {
                input: args.input,
                output: args.output,
                working_dir: args.working_dir,
                threads: args.threads,
            };
            bz_lib::decompress(&config)?;
            info!("Decompression complete!");
        }
        Commands::Extract(args) => {
            info!("Starting BAM to QZ extraction...");
            let config = bz_lib::ExtractConfig {
                input: args.input,
                output_prefix: args.output,
                working_dir: args.working_dir,
                threads: args.threads,
            };
            bz_lib::extract(&config)?;
            info!("Extraction complete!");
        }
    }

    Ok(())
}
