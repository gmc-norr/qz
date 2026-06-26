//! Lever-A decision gate: does scheduling MANY 256 MB BSC flushes concurrently
//! recover the throughput the fused pipeline loses to small-flush work-starvation?
//!
//! The fused reorder pipeline flushes a 256 MB buffer (~11 blocks ≈ 11/72 cores)
//! synchronously in the gather thread → BSC runs at ~149 MB/s (measured
//! in-pipeline 143 GB / 961 s = 148.8 MB/s, and a standalone single 256 MB blob
//! @72t = 149 MB/s — two independent paths, exact match). Lever A decouples
//! `flush_bsc` onto a bounded worker pool so N flushes compress concurrently,
//! putting N×11 blocks on the rayon pool at once.
//!
//! This bench reproduces that: for N in a sweep, spawn N OS threads each calling
//! `bsc::compress_parallel` on a distinct warm 256 MB slice, time the batch, and
//! report AGGREGATE MB/s. Reading per the spec gate (§8k):
//!   N=1  → should reproduce ~149 MB/s (the in-pipeline number)
//!   N≥3  → ~330+  ⇒ join-barriers don't re-serialize, BUILD lever A (~1.3×)
//!          ~200   ⇒ barrier caps recovery at ~1.15×, reconsider
//!
//! Usage: bench_bsc_concurrent <blob> [slice_mb=256] [maxN=4]

use qz_lib::compression::bsc;
use std::fs;
use std::time::Instant;

fn main() {
    let path = std::env::args().nth(1).expect("usage: bench_bsc_concurrent <blob> [slice_mb] [maxN]");
    let slice_mb: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(256);
    let max_n: usize = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(4);
    let slice_bytes = slice_mb * 1024 * 1024;

    let data = fs::read(&path).expect("read blob");
    let avail = data.len() / slice_bytes;
    assert!(avail >= 1, "blob too small for one {slice_mb} MB slice");
    eprintln!(
        "[in] {} bytes ({:.2} GB) → {} distinct {} MB slices available",
        data.len(),
        data.len() as f64 / 1e9,
        avail,
        slice_mb,
    );

    // Carve up to max_n distinct slices (wrap if the blob is smaller than
    // max_n * slice). Distinct windows avoid any shared read-cache super-effect.
    let make_slice = |i: usize| -> &[u8] {
        let start = (i % avail) * slice_bytes;
        &data[start..start + slice_bytes]
    };

    // Warm: touch every byte we'll use (the blob is already read() into RAM, but
    // make the page residency explicit so N=1 and N=4 see identical conditions).
    {
        let mut acc = 0u8;
        for i in 0..max_n {
            for &b in make_slice(i).iter().step_by(4096) {
                acc ^= b;
            }
        }
        std::hint::black_box(acc);
    }
    eprintln!("warmed.\n");

    let single_mbps = (slice_bytes as f64) / 1e6; // numerator for N=1 reference

    println!(
        "{:>3}  {:>10}  {:>12}  {:>10}  {:>12}",
        "N", "wall_s", "agg_MB/s", "per_MB/s", "speedup_vs_N1"
    );
    let mut n1_agg = 0.0f64;
    for n in 1..=max_n {
        let slices: Vec<&[u8]> = (0..n).map(make_slice).collect();
        let t = Instant::now();
        let comps: Vec<usize> = std::thread::scope(|s| {
            let handles: Vec<_> = slices
                .iter()
                .map(|sl| s.spawn(move || bsc::compress_parallel(sl).expect("bsc").len()))
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        let wall = t.elapsed().as_secs_f64();
        let total_in = (n * slice_bytes) as f64;
        let agg = total_in / 1e6 / wall;
        let per = single_mbps / wall; // throughput one flush would see if alone in this wall
        if n == 1 {
            n1_agg = agg;
        }
        let _ = comps;
        println!(
            "{:>3}  {:>10.2}  {:>12.0}  {:>10.0}  {:>12.2}",
            n, wall, agg, per, agg / n1_agg
        );
    }
    println!("\nGate (§8k): N=1 ~149 confirms in-pipeline number; N≥3 ≥330 ⇒ BUILD lever A.");
}
