//! Measure the reference-aware archive **merge tax** (the serial tail the NUMA
//! reference-compress shard path would pay): re-deriving the coverage globals over
//! the union interval map + relocating every per-chunk frame.
//!
//! Usage: bench_reference_merge <ref.fa> <out.qz> <shard0.qz> <shard1.qz> [...]
//!
//! The shard archives must already exist (build them with `qz compress --reference
//! <ref.fa>` over disjoint read slices). Reports wall time + merged size so the
//! end-to-end win (`T_inproc / (T_worker + T_merge)`) can be computed.

use mimalloc::MiMalloc;
use std::path::PathBuf;
use std::time::Instant;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!("Usage: bench_reference_merge <ref.fa> <out.qz> <shard0.qz> <shard1.qz> [...]");
        std::process::exit(1);
    }
    let reference = PathBuf::from(&args[1]);
    let out = PathBuf::from(&args[2]);
    let shards: Vec<PathBuf> = args[3..].iter().map(PathBuf::from).collect();

    let t = Instant::now();
    qz_lib::compression::merge_reference_archives_to_path(&shards, &reference, &out, true)
        .expect("reference merge failed");
    let el = t.elapsed();

    let sz = std::fs::metadata(&out).map(|m| m.len()).unwrap_or(0);
    eprintln!(
        "merge_reference: {} shards -> {} in {:.3}s (merged {} bytes)",
        shards.len(),
        out.display(),
        el.as_secs_f64(),
        sz
    );
}
