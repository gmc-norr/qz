//! Time BSC (qz production 25 MB block-parallel) compress + decompress on a RAW
//! file blob, so it can be compared apples-to-apples against `zstd --long` on the
//! same blob (the clustered sequence stream). Reports sizes + wall + throughput,
//! and asserts a lossless roundtrip.
//!
//! Usage: bench_bsc_raw <blob>

use qz_lib::compression::bsc;
use std::fs;
use std::time::Instant;

fn main() {
    let path = std::env::args().nth(1).expect("usage: bench_bsc_raw <blob>");
    let data = fs::read(&path).expect("read blob");
    let n = data.len() as f64;
    eprintln!("[in] {} bytes ({:.2} GB)", data.len(), n / 1e9);

    let tc = Instant::now();
    let comp = bsc::compress_parallel(&data).expect("bsc compress");
    let ct = tc.elapsed().as_secs_f64();

    let td = Instant::now();
    let back = bsc::decompress_parallel(&comp).expect("bsc decompress");
    let dt = td.elapsed().as_secs_f64();

    assert_eq!(back, data, "BSC roundtrip MISMATCH (not lossless)");

    let mbps = |bytes: f64, secs: f64| bytes / 1e6 / secs;
    println!(
        "BSC(25MB)  comp_bytes={}  bpb={:.3}  comp={:.1}s ({:.0} MB/s)  decomp={:.1}s ({:.0} MB/s)  ROUNDTRIP_OK",
        comp.len(),
        comp.len() as f64 * 8.0 / n, // bits per input byte (≈ bits/base for a pure-ACGT blob)
        ct,
        mbps(n, ct),
        dt,
        mbps(n, dt),
    );
}
