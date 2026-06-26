//! Bench bounded streaming decompression: wall-clock + peak RSS.
//!
//! Usage: cargo run --release -p qz-bench --bin bench_bounded_streaming -- <archive.qz>

use std::process::exit;
use std::time::Instant;

#[cfg(target_os = "linux")]
fn read_vm_hwm_kb() -> u64 {
    let s = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    for line in s.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            return rest
                .split_whitespace()
                .next()
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
        }
    }
    0
}

#[cfg(not(target_os = "linux"))]
fn read_vm_hwm_kb() -> u64 {
    0
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <archive.qz>", args[0]);
        exit(1);
    }
    let archive = std::path::PathBuf::from(&args[1]);
    if !archive.exists() {
        eprintln!("Archive not found: {}", archive.display());
        exit(1);
    }
    let temp = tempfile::TempDir::new()?;
    let decoded = temp.path().join("decoded.fastq");

    let t = Instant::now();
    qz_lib::compression::decompress(&qz_lib::cli::DecompressConfig {
        input: archive,
        output: vec![decoded],
        working_dir: temp.path().to_path_buf(),
        num_threads: 0,
        gzipped: false,
        gzip_level: 6,
        force: false,
    })?;
    let wall = t.elapsed();
    let peak_rss_mb = read_vm_hwm_kb() / 1024;

    println!("wall_clock_ms = {}", wall.as_millis());
    println!("peak_rss_mb   = {}", peak_rss_mb);
    Ok(())
}
