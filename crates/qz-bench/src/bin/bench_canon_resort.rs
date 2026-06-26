//! Follow-up to bench_rc_canon: RC-canonicalize clustered reads, THEN re-sort by a
//! selectable key to bring co-oriented near-duplicates adjacent (feeds zstd's
//! repeat-offset path). Read order is DROPPED in this mode, so re-sorting costs
//! NOTHING (no permutation stored) — we are free to pick the most
//! compression-friendly order.
//!
//! QZ_SORT (env):
//!   none     keep the existing clustered order (== RC-canon only; reproduces baseline)
//!   seq      global lexicographic sort by the canonicalized read bytes
//!   hashseq  sort by (min canonical syncmer hash, canonicalized read bytes)
//!            — keeps the locus grouping, refines within-locus by sequence
//!
//! Input: clustered fixed-length separator-free seq blob (clust.seq). Output: the
//! canonicalized+resorted blob (zstd it externally, compare to 122.45 MB cluster
//! baseline and 117.88 MB RC-canon-only). Strand bits (1/read, in the NEW order)
//! BSC'd and reported as the tax to subtract.
//!
//! Usage: QZ_SORT=hashseq bench_canon_resort <clust.seq> <readlen> <out.seq>

use rayon::prelude::*;
use std::fs;
use std::time::Instant;

use qz_lib::compression::bsc;
use qz_lib::compression::dna_utils::{
    canonical_minimizer_orientation, compute_min_syncmer_hash_pos, reverse_complement_canonical,
};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: QZ_SORT=none|seq|hashseq bench_canon_resort <clust.seq> <readlen> <out.seq>");
        std::process::exit(2);
    }
    let inpath = &args[1];
    let readlen: usize = args[2].parse().expect("readlen");
    let outpath = &args[3];
    let mode = std::env::var("QZ_SORT").unwrap_or_else(|_| "hashseq".into());
    // QZ_NOCANON=1 -> skip RC-canon (raw orientation); RC-invariant hash still clusters.
    // Measures how much co-orientation buys UNDER each sort mode (kill-gate #2).
    let nocanon = std::env::var("QZ_NOCANON").map(|v| v == "1").unwrap_or(false);

    let t0 = Instant::now();
    let data = fs::read(inpath).expect("read input");
    assert_eq!(data.len() % readlen, 0, "input not a multiple of readlen");
    let nreads = data.len() / readlen;
    eprintln!(
        "[load] {} reads x {}bp = {:.2} Gbase ({:.0}s)",
        nreads, readlen, data.len() as f64 / 1e9, t0.elapsed().as_secs_f64()
    );

    // RC-canonicalize in place into a contiguous buffer; record (flip, min_hash).
    let t1 = Instant::now();
    let mut canon = vec![0u8; data.len()];
    let metas: Vec<(bool, u64, u16)> = canon
        .par_chunks_mut(readlen)
        .zip(data.par_chunks(readlen))
        .map(|(out, read)| {
            let flip = if nocanon { false } else { canonical_minimizer_orientation(read) };
            if flip {
                out.copy_from_slice(&reverse_complement_canonical(read));
            } else {
                out.copy_from_slice(read);
            }
            let (h, pos) = compute_min_syncmer_hash_pos(out); // RC-invariant hash either way
            (flip, h, pos)
        })
        .collect();
    let flipped = metas.iter().filter(|(f, _, _)| *f).count();
    eprintln!(
        "[canon] {} flipped ({:.1}%)  {:.0}s",
        flipped, 100.0 * flipped as f64 / nreads as f64, t1.elapsed().as_secs_f64()
    );

    // Build the read order.
    let t2 = Instant::now();
    let mut idx: Vec<u32> = (0..nreads as u32).collect();
    let rd = |i: u32| -> &[u8] {
        let s = i as usize * readlen;
        &canon[s..s + readlen]
    };
    // suffix from the minimizer anchor (aligns shifted overlaps at the shared seed)
    let anchored = |i: u32| -> &[u8] {
        let s = i as usize * readlen;
        let p = (metas[i as usize].2 as usize).min(readlen);
        &canon[s + p..s + readlen]
    };
    match mode.as_str() {
        "none" => {}
        "seq" => idx.par_sort_unstable_by(|&a, &b| rd(a).cmp(rd(b))),
        "hashseq" => idx.par_sort_unstable_by(|&a, &b| {
            metas[a as usize].1.cmp(&metas[b as usize].1).then_with(|| rd(a).cmp(rd(b)))
        }),
        "anchorseq" => idx.par_sort_unstable_by(|&a, &b| {
            metas[a as usize].1.cmp(&metas[b as usize].1).then_with(|| anchored(a).cmp(anchored(b)))
        }),
        // bidirectional anchor: align at the shared minimizer, order by downstream suffix,
        // then by the upstream prefix read OUTWARD from the anchor (both flanks).
        "anchorbi" => idx.par_sort_unstable_by(|&a, &b| {
            metas[a as usize].1.cmp(&metas[b as usize].1)
                .then_with(|| anchored(a).cmp(anchored(b)))
                .then_with(|| {
                    let sa = a as usize * readlen;
                    let sb = b as usize * readlen;
                    let pa = (metas[a as usize].2 as usize).min(readlen);
                    let pb = (metas[b as usize].2 as usize).min(readlen);
                    canon[sa..sa + pa].iter().rev().cmp(canon[sb..sb + pb].iter().rev())
                })
        }),
        other => panic!("unknown QZ_SORT={other}"),
    }
    eprintln!("[sort:{}] {:.0}s", mode, t2.elapsed().as_secs_f64());

    // Emit reordered canonical blob + reordered strand bits.
    let mut out = Vec::with_capacity(data.len());
    let mut bits = vec![0u8; nreads.div_ceil(8)];
    for (newi, &orig) in idx.iter().enumerate() {
        out.extend_from_slice(rd(orig));
        if metas[orig as usize].0 {
            bits[newi >> 3] |= 1 << (newi & 7);
        }
    }
    fs::write(outpath, &out).expect("write blob");
    let bits_bsc = bsc::compress_parallel(&bits).expect("bsc bits");
    eprintln!(
        "[write] {} bytes -> {} | strand-bits BSC {} B ({:.2} MB)",
        out.len(), outpath, bits_bsc.len(), bits_bsc.len() as f64 / 1e6
    );
    eprintln!("[done] total {:.0}s", t0.elapsed().as_secs_f64());
    println!("CANON_RESORT mode={} nreads={} flipped={} strand_bits_bsc={}", mode, nreads, flipped, bits_bsc.len());
}
