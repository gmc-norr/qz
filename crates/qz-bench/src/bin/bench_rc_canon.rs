//! FTO-clear "help zstd" probe: reverse-complement canonicalization of clustered
//! reads. zstd (and BSC) match raw bytes and cannot see that a read and its
//! reverse-complement are the same locus, so at high coverage ~half the donors
//! at each locus are invisible. This flips each read to a canonical orientation
//! (so its minimum canonical syncmer — the SAME minimizer used to cluster — reads
//! forward) and stores 1 strand bit/read. The clustering ORDER is preserved
//! (canonical minimizer is orientation-invariant), so this isolates the pure
//! RC effect from any re-sort.
//!
//! Input: a CLUSTERED, fixed-length, separator-free seq blob (e.g. clust.seq,
//! 1 byte/base, all reads `--readlen` long). Output: the canonicalized blob (zstd
//! it externally and compare to the un-canonicalized baseline) + the BSC cost of
//! the strand bits (must be subtracted from any win).
//!
//! Usage: bench_rc_canon <clust.seq> <readlen> <out.canon.seq>

use rayon::prelude::*;
use std::fs;
use std::time::Instant;

use qz_lib::compression::bsc;
use qz_lib::compression::dna_utils::{canonical_minimizer_orientation, reverse_complement_canonical};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: bench_rc_canon <clust.seq> <readlen> <out.canon.seq>");
        std::process::exit(2);
    }
    let inpath = &args[1];
    let readlen: usize = args[2].parse().expect("readlen");
    let outpath = &args[3];

    let t0 = Instant::now();
    let data = fs::read(inpath).expect("read input");
    assert_eq!(
        data.len() % readlen,
        0,
        "input {} not a multiple of readlen {}",
        data.len(),
        readlen
    );
    let nreads = data.len() / readlen;
    eprintln!(
        "[load] {} reads x {}bp = {:.2} Gbase ({:.0}s)",
        nreads,
        readlen,
        data.len() as f64 / 1e9,
        t0.elapsed().as_secs_f64()
    );

    // Parallel: per read, decide orientation and emit canonical bytes + strand bit.
    let t1 = Instant::now();
    let results: Vec<(Vec<u8>, bool)> = data
        .par_chunks(readlen)
        .map(|read| {
            if canonical_minimizer_orientation(read) {
                (reverse_complement_canonical(read), true)
            } else {
                (read.to_vec(), false)
            }
        })
        .collect();
    let flipped = results.iter().filter(|(_, b)| *b).count();
    eprintln!(
        "[canon] {} reads flipped ({:.1}%)  {:.0}s",
        flipped,
        100.0 * flipped as f64 / nreads as f64,
        t1.elapsed().as_secs_f64()
    );

    // Assemble canonicalized blob in original (clustered) order.
    let mut out = Vec::with_capacity(data.len());
    let mut bits = vec![0u8; nreads.div_ceil(8)];
    for (i, (seq, flip)) in results.iter().enumerate() {
        out.extend_from_slice(seq);
        if *flip {
            bits[i >> 3] |= 1 << (i & 7);
        }
    }
    fs::write(outpath, &out).expect("write canon blob");
    eprintln!("[write] {} bytes -> {}", out.len(), outpath);

    // Cost of the strand-bit side stream (BSC, what qz would use).
    let bits_bsc = bsc::compress_parallel(&bits).expect("bsc bits");
    eprintln!(
        "[strand-bits] raw {} B ({:.2} MB) -> BSC {} B ({:.2} MB)  = the tax to subtract",
        bits.len(),
        bits.len() as f64 / 1e6,
        bits_bsc.len(),
        bits_bsc.len() as f64 / 1e6
    );

    eprintln!("[done] total {:.0}s", t0.elapsed().as_secs_f64());
    println!(
        "RC_CANON nreads={} flipped={} strand_bits_bsc={}",
        nreads,
        flipped,
        bits_bsc.len()
    );
}
