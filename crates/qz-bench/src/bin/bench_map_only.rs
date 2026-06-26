//! Standalone profiling of the strobealign map-only path qz uses in
//! reference-based compression. Isolates JUST the mapper (index load + seeding +
//! find_hits + chaining + pairing) — no qz diff/encode — so we can see where the
//! mapper's own time goes and whether it's worth optimizing.
//!
//! Usage:
//!   bench_map_only <ref.fa> <R1.fastq> <R2.fastq> [n_pairs] [--fast]
//!
//! Mirrors `reference::mapping::Mapper::build` + `map_pair` exactly: same
//! normalization, same SeedingParameters, same McsStrategy::Always, same
//! mu/sigma=300/100, loading the prebuilt `.sti` sidecar if present.

use std::io::{BufRead, BufReader};
use std::time::Instant;

// Match the production `qz` binary (qz-cli sets this) so allocator cost in the
// profile is production-faithful, not glibc-inflated.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use rayon::prelude::*;
use strobealign::chainer::{Chainer, ChainingParameters};
use strobealign::index::StrobemerIndex;
use strobealign::io::fasta::{RefSequence, read_ref};
use strobealign::maponly::{DEFAULT_RESCUE_DISTANCE, map_pair_deterministic};
use strobealign::mcsstrategy::McsStrategy;
use strobealign::nam::get_nams_by_chaining;
use strobealign::seeding::SeedingParameters;

fn normalize_ref(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        *b = match b.to_ascii_uppercase() {
            c @ (b'A' | b'C' | b'G' | b'T') => c,
            _ => b'N',
        };
    }
}

/// Read up to `n` sequence lines (the 2nd of every 4) from a plain FASTQ file.
fn read_seqs(path: &str, n: usize) -> Vec<Vec<u8>> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {path}: {e}"));
    let mut r = BufReader::with_capacity(1 << 20, f);
    let mut out = Vec::with_capacity(n);
    let mut line = Vec::new();
    let mut which = 0u8;
    while out.len() < n {
        line.clear();
        let read = r.read_until(b'\n', &mut line).unwrap();
        if read == 0 {
            break;
        }
        if which == 1 {
            // strip trailing \n (and \r)
            while matches!(line.last(), Some(b'\n' | b'\r')) {
                line.pop();
            }
            out.push(line.clone());
        }
        which = (which + 1) % 4;
    }
    out
}

fn sidecar_path(reference: &str, read_len: usize, tag: &str) -> std::path::PathBuf {
    std::path::PathBuf::from(format!("{reference}.qz-r{read_len}{tag}.sti"))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: bench_map_only <ref.fa> <R1.fastq> <R2.fastq> [n_pairs] [--fast]");
        std::process::exit(2);
    }
    let ref_path = &args[1];
    let r1_path = &args[2];
    let r2_path = &args[3];
    let n_pairs: usize = args
        .get(4)
        .and_then(|s| s.parse().ok())
        .unwrap_or(1_000_000);
    let fast = args.iter().any(|a| a == "--fast");

    let threads = rayon::current_num_threads();
    eprintln!("bench_map_only: n_pairs={n_pairs} fast={fast} rayon_threads={threads}");

    // ---- Load reads ----
    let t = Instant::now();
    let (r1, r2) = rayon::join(
        || read_seqs(r1_path, n_pairs),
        || read_seqs(r2_path, n_pairs),
    );
    let n = r1.len().min(r2.len());
    let pairs: Vec<(Vec<u8>, Vec<u8>)> = r1.into_iter().zip(r2).take(n).collect();
    let read_len = pairs.first().map(|(a, _)| a.len()).unwrap_or(0);
    eprintln!(
        "loaded {} pairs (read_len={read_len}) in {:.2}s",
        pairs.len(),
        t.elapsed().as_secs_f64()
    );

    // ---- Build mapper config (mirrors Mapper::build) ----
    let t = Instant::now();
    let mut references: Vec<RefSequence> = read_ref(ref_path).unwrap();
    for r in references.iter_mut() {
        normalize_ref(&mut r.sequence);
    }
    let mut params = SeedingParameters::new(read_len);
    let mut tag = "";
    if fast {
        params = params.with_k_s(None, Some(12)).unwrap();
        tag = "-fast";
    }
    let bits = params.syncmer.pick_bits(&references);
    let sidecar = sidecar_path(ref_path, read_len, tag);
    let index: StrobemerIndex = if sidecar.exists() {
        match strobealign::index::read_index(&sidecar, params.clone(), bits) {
            Ok(ix) => {
                eprintln!("loaded cached index {}", sidecar.display());
                ix
            }
            Err(e) => panic!("cached index unusable: {e}"),
        }
    } else {
        eprintln!("building index (no sidecar at {})", sidecar.display());
        strobealign::indexer::make_index(&references, params.clone(), bits, 0.0002, threads).0
    };
    let chainer = Chainer::new(index.k(), ChainingParameters::default());
    let mcs = McsStrategy::default();
    let (mu, sigma) = (300.0f32, 100.0f32);
    eprintln!(
        "mapper ready in {:.2}s (k={}, randstrobes={})",
        t.elapsed().as_secs_f64(),
        index.k(),
        index.len()
    );

    // ---- A) Multi-threaded throughput (production-shaped: par_iter map_pair) ----
    // Collect mappings in order (par_iter collect is order-preserving) so we can
    // checksum them: ANY change to the mapper that alters output flips the digest.
    // The tuple mirrors `to_mapping` (ref_id, projected_ref_start, is_revcomp).
    let t = Instant::now();
    let maps: Vec<(Option<(u32, u64, bool)>, Option<(u32, u64, bool)>)> = pairs
        .par_iter()
        .map(|(a, b)| {
            let (n1, n2) = map_pair_deterministic(
                a,
                b,
                &index,
                &chainer,
                DEFAULT_RESCUE_DISTANCE,
                mcs,
                mu,
                sigma,
            );
            let f = |n: Option<strobealign::nam::Nam>| {
                n.map(|n| (n.ref_id as u32, n.projected_ref_start() as u64, n.is_revcomp))
            };
            (f(n1), f(n2))
        })
        .collect();
    let dt = t.elapsed().as_secs_f64();
    let mut mapped = 0usize;
    let mut digest: u64 = 0xcbf29ce484222325; // FNV-1a 64 offset basis
    let mix = |x: u64, d: &mut u64| {
        *d ^= x;
        *d = d.wrapping_mul(0x100000001b3);
    };
    for (m1, m2) in &maps {
        for m in [m1, m2] {
            match m {
                Some((rid, rs, rc)) => {
                    mapped += 1;
                    mix(1, &mut digest);
                    mix(*rid as u64, &mut digest);
                    mix(*rs, &mut digest);
                    mix(*rc as u64, &mut digest);
                }
                None => mix(0, &mut digest),
            }
        }
    }
    eprintln!(
        "\n[MT map_pair] {} pairs in {:.3}s = {:.0} pairs/s ({:.0} reads/s); mapped_mates={}",
        pairs.len(),
        dt,
        pairs.len() as f64 / dt,
        2.0 * pairs.len() as f64 / dt,
        mapped
    );
    eprintln!("[CHECKSUM] mapping digest = {digest:#018x}  (must be invariant across optimizations)");

    // ---- B) Single-threaded throughput on a subset (clean per-pair cost) ----
    let sub = pairs.len().min(100_000);
    let t = Instant::now();
    let mut acc = 0usize;
    for (a, b) in &pairs[..sub] {
        let (n1, n2) = map_pair_deterministic(
            a,
            b,
            &index,
            &chainer,
            DEFAULT_RESCUE_DISTANCE,
            mcs,
            mu,
            sigma,
        );
        acc += n1.is_some() as usize + n2.is_some() as usize;
    }
    let dt = t.elapsed().as_secs_f64();
    eprintln!(
        "[ST map_pair] {} pairs in {:.3}s = {:.0} pairs/s (1 thread); mapped_mates={}",
        sub,
        dt,
        sub as f64 / dt,
        acc
    );

    // ---- C) Per-stage breakdown via the mapper's own NamDetails timers ----
    // get_nams_by_chaining records time_randstrobes (seeding), time_find_hits
    // (index lookup), time_chaining. Sum over a single-threaded subset of mates.
    let mut t_seed = 0.0;
    let mut t_hits = 0.0;
    let mut t_chain = 0.0;
    let mut n_rs = 0usize;
    let mut n_anchors = 0usize;
    let mut n_nams = 0usize;
    let t = Instant::now();
    for (a, b) in &pairs[..sub] {
        for seq in [a, b] {
            let (d, nams) =
                get_nams_by_chaining(seq, &index, &chainer, DEFAULT_RESCUE_DISTANCE, mcs);
            t_seed += d.time_randstrobes;
            t_hits += d.time_find_hits;
            t_chain += d.time_chaining;
            n_rs += d.n_randstrobes;
            n_anchors += d.n_anchors;
            n_nams += nams.len();
        }
    }
    let wall = t.elapsed().as_secs_f64();
    let sum = t_seed + t_hits + t_chain;
    eprintln!(
        "\n[per-stage, ST, {} mates, wall {:.3}s] instrumented-sum {:.3}s:",
        2 * sub,
        wall,
        sum
    );
    eprintln!(
        "  seeding (randstrobes_query) : {:.3}s  {:5.1}%",
        t_seed,
        100.0 * t_seed / sum
    );
    eprintln!(
        "  find_hits (index lookups)   : {:.3}s  {:5.1}%",
        t_hits,
        100.0 * t_hits / sum
    );
    eprintln!(
        "  chaining                    : {:.3}s  {:5.1}%",
        t_chain,
        100.0 * t_chain / sum
    );
    eprintln!(
        "  avg per mate: {:.1} randstrobes, {:.1} anchors, {:.2} NAMs",
        n_rs as f64 / (2 * sub) as f64,
        n_anchors as f64 / (2 * sub) as f64,
        n_nams as f64 / (2 * sub) as f64,
    );
}
