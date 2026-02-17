/// Benchmark quality_ctx block size impact on compression ratio and speed.
/// Run: cargo run --release -p qz-bench --bin bench_qctx_blocksize
use std::path::PathBuf;
use std::time::Instant;

fn main() -> anyhow::Result<()> {
    let input_file = std::env::args().nth(1).unwrap_or_else(|| "real_data/ERR3239334_1.1m.fastq".to_string());
    let input = PathBuf::from(&input_file);
    if !input.exists() {
        eprintln!("Need real_data/ERR3239334_1.1m.fastq");
        std::process::exit(1);
    }

    let block_sizes: Vec<usize> = std::env::args().nth(2)
        .map(|s| s.split(',').filter_map(|x| x.trim().parse().ok()).collect())
        .unwrap_or_else(|| vec![50_000, 100_000, 150_000, 250_000, 500_000, 1_000_000]);

    println!("{:<12} {:>12} {:>10} {:>12} {:>12}", "BlockSize", "QualBytes", "Ratio", "CompTime", "DecTime");
    println!("{}", "-".repeat(62));

    for &bs in &block_sizes {
        let out_path = PathBuf::from(format!("/tmp/qctx_bs_{}.qz", bs));
        let dec_path = PathBuf::from(format!("/tmp/qctx_bs_{}.fastq", bs));

        let config = qz_lib::cli::CompressConfig {
            input: vec![input.clone()],
            output: out_path.clone(),
            advanced: qz_lib::cli::AdvancedOptions {
                quality_ctx_block_size: bs,
                ..Default::default()
            },
            ..qz_lib::cli::CompressConfig::default()
        };

        // Compress
        let t0 = Instant::now();
        qz_lib::compression::compress(&config)?;
        let comp_time = t0.elapsed().as_secs_f64();

        // Get compressed size
        let compressed_size = std::fs::metadata(&out_path)?.len();
        let original_size = std::fs::metadata(&input)?.len();
        let ratio = original_size as f64 / compressed_size as f64;

        // Read quality size from archive header (offset 44..52)
        let mut f = std::fs::File::open(&out_path)?;
        use std::io::{Read, Seek, SeekFrom};
        f.seek(SeekFrom::Start(44))?;
        let mut buf = [0u8; 8];
        f.read_exact(&mut buf)?;
        let qual_size = u64::from_le_bytes(buf);

        // Decompress
        let dec_config = qz_lib::cli::DecompressConfig {
            input: out_path.clone(),
            output: vec![dec_path.clone()],
            working_dir: PathBuf::from("."),
            num_threads: 0,
            gzipped: false,
            gzip_level: 6,
        };
        let t1 = Instant::now();
        qz_lib::compression::decompress(&dec_config)?;
        let dec_time = t1.elapsed().as_secs_f64();

        println!("{:<12} {:>12} {:>10.3}x {:>11.2}s {:>11.2}s",
            format!("{}K", bs / 1000), qual_size, ratio, comp_time, dec_time);

        // Cleanup
        let _ = std::fs::remove_file(&out_path);
        let _ = std::fs::remove_file(&dec_path);
    }

    Ok(())
}
