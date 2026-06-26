//! Probe: does a shared-mate fqzcomp blob with FQZ_FREAD2 on mate-2 reads beat
//! production's current separate-mate, mean-quality-sorted blobs?
//!
//! See docs/superpowers/specs/2026-06-09-fqzcomp-read-flags-design.md (rev 2).
//!
//! Both arms reuse the EXACT production framing of
//! `compress_qualities_fqzcomp_quals` (codecs.rs): per-read mean-quality stable sort,
//! 500K-record sub-chunks, fqzcomp per sub-chunk, plus the BSC'd sort-key side
//! stream that lets decode recompute the permutation. The only difference:
//!   Arm A — R1 and R2 encoded as SEPARATE sorted blob sets, flags = 0 (== prod).
//!   Arm B — R1+R2 in ONE sorted blob set, FQZ_FREAD2 on R2, + a BSC'd mate map.
//!
//! Usage: bench_fqz_flags <R1.fastq> <R2.fastq> [num_pairs] [reps]

use anyhow::{Context, Result, bail};
use qz_lib::compression::{bsc, fqzcomp};
use rayon::prelude::*;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const FQZCOMP_SUB_CHUNK: usize = 500_000;
const FQZ_FREAD2: u32 = 128; // fqzcomp_qual.h
const DEFAULT_PAIRS: usize = 2_000_000;
const DEFAULT_REPS: usize = 3;

/// Production's per-read sort key: mean quality, one byte.
fn mean_key(q: &[u8]) -> u8 {
    let sum: u64 = q.iter().map(|&b| b as u64).sum();
    (sum / q.len().max(1) as u64) as u8
}

/// Read the first `n` quality strings from a FASTQ file (every 4th line).
fn read_qualities(path: &str, n: usize) -> Result<Vec<Vec<u8>>> {
    let f = std::fs::File::open(path).with_context(|| format!("open {path}"))?;
    let mut r = BufReader::with_capacity(1 << 20, f);
    let mut out = Vec::with_capacity(n);
    let mut line = Vec::new();
    let mut line_no = 0usize;
    while out.len() < n {
        line.clear();
        let read = r.read_until(b'\n', &mut line)?;
        if read == 0 {
            break;
        }
        line_no += 1;
        if line_no.is_multiple_of(4) {
            // quality line
            while matches!(line.last(), Some(b'\n' | b'\r')) {
                line.pop();
            }
            out.push(std::mem::take(&mut line));
        }
    }
    Ok(out)
}

/// One production-faithful fqzcomp-with-sort encoding of a read set.
struct Encoded {
    blobs: Vec<Vec<u8>>, // fqz blob per 500K sub-chunk, in sorted order
    keys_bsc: Vec<u8>,   // BSC'd per-read mean keys (original order) = permutation side channel
    n: usize,
}

impl Encoded {
    fn total(&self) -> usize {
        self.blobs.iter().map(|b| b.len()).sum::<usize>() + self.keys_bsc.len()
    }
}

/// Encode a read set with production framing. If `flags` is given, the per-read
/// flags are permuted into sorted order and passed to `compress_with_flags`.
fn encode_set(quals: &[&[u8]], flags: Option<&[u32]>) -> Result<Encoded> {
    let n = quals.len();
    let keys: Vec<u8> = quals.par_iter().map(|q| mean_key(q)).collect();

    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by_key(|&i| keys[i]); // stable — matches production

    let sorted: Vec<&[u8]> = idx.iter().map(|&i| quals[i]).collect();
    let sorted_flags: Option<Vec<u32>> = flags.map(|f| idx.iter().map(|&i| f[i]).collect());

    // Sub-chunk ranges, fqz each in parallel (mirrors production par_chunks).
    let starts: Vec<usize> = (0..n).step_by(FQZCOMP_SUB_CHUNK).collect();
    let (blobs, keys_bsc) = rayon::join(
        || -> Result<Vec<Vec<u8>>> {
            starts
                .par_iter()
                .map(|&s| {
                    let e = (s + FQZCOMP_SUB_CHUNK).min(n);
                    match &sorted_flags {
                        Some(sf) => fqzcomp::compress_with_flags(&sorted[s..e], &sf[s..e], 0),
                        None => fqzcomp::compress(&sorted[s..e], 0),
                    }
                })
                .collect()
        },
        || bsc::compress_parallel_adaptive(&keys),
    );

    Ok(Encoded {
        blobs: blobs?,
        keys_bsc: keys_bsc?,
        n,
    })
}

/// Recompute the sort permutation from the (decoded) keys and invert the encode,
/// returning the read set in ORIGINAL order. Proves the side channel suffices.
fn decode_set(enc: &Encoded) -> Result<Vec<Vec<u8>>> {
    let keys = bsc::decompress_parallel(&enc.keys_bsc)?;
    if keys.len() != enc.n {
        bail!("decoded keys len {} != n {}", keys.len(), enc.n);
    }
    let mut idx: Vec<usize> = (0..enc.n).collect();
    idx.sort_by_key(|&i| keys[i]); // same stable sort the encoder used

    // Decode all sub-chunks → sorted-order quality strings.
    let mut sorted: Vec<Vec<u8>> = Vec::with_capacity(enc.n);
    let mut produced = 0usize;
    for blob in &enc.blobs {
        let remaining = enc.n - produced;
        let this = remaining.min(FQZCOMP_SUB_CHUNK);
        let (concat, lengths) = fqzcomp::decompress(blob, this)?;
        let mut off = 0usize;
        for &l in &lengths {
            let l = l as usize;
            sorted.push(concat[off..off + l].to_vec());
            off += l;
        }
        produced += this;
    }
    if sorted.len() != enc.n {
        bail!("decoded {} reads, expected {}", sorted.len(), enc.n);
    }

    // Invert: original[idx[k]] = sorted[k]
    let mut original = vec![Vec::new(); enc.n];
    for (k, &orig_i) in idx.iter().enumerate() {
        original[orig_i] = std::mem::take(&mut sorted[k]);
    }
    Ok(original)
}

/// Bit-pack mate flags (1 = R2) and BSC them — the demux side channel a real
/// shared-mate blob must store to split mates on decode.
fn mate_map_bsc(n_total: usize, n_r1: usize) -> Result<Vec<u8>> {
    let mut bits = vec![0u8; n_total.div_ceil(8)];
    for i in n_r1..n_total {
        bits[i / 8] |= 1u8 << (i % 8);
    }
    bsc::compress_parallel_adaptive(&bits)
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: {} <R1.fastq> <R2.fastq> [num_pairs] [reps]",
            args[0]
        );
        std::process::exit(2);
    }
    let r1_path = &args[1];
    let r2_path = &args[2];
    let pairs = args
        .get(3)
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_PAIRS);
    let reps = args
        .get(4)
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_REPS);

    eprintln!("Reading {pairs} quality strings from each of R1/R2 ...");
    let (r1, r2) = rayon::join(
        || read_qualities(r1_path, pairs),
        || read_qualities(r2_path, pairs),
    );
    let r1 = r1?;
    let r2 = r2?;
    if r1.len() != r2.len() {
        bail!(
            "R1 has {} reads but R2 has {} — not paired by index",
            r1.len(),
            r2.len()
        );
    }
    let n_r1 = r1.len();
    let n_total = n_r1 * 2;
    let raw_bytes: usize =
        r1.iter().map(|q| q.len()).sum::<usize>() + r2.iter().map(|q| q.len()).sum::<usize>();
    eprintln!("Loaded {n_r1} pairs ({n_total} reads, {raw_bytes} raw quality bytes)\n");

    // Borrowed views.
    let r1_ref: Vec<&[u8]> = r1.iter().map(|v| v.as_slice()).collect();
    let r2_ref: Vec<&[u8]> = r2.iter().map(|v| v.as_slice()).collect();

    // Combined set (R1 ++ R2) + flags for Arm B.
    let combined: Vec<&[u8]> = r1_ref.iter().chain(r2_ref.iter()).copied().collect();
    let mut flags = vec![0u32; n_total];
    for f in flags.iter_mut().take(n_total).skip(n_r1) {
        *f = FQZ_FREAD2;
    }

    let mut a_best = usize::MAX;
    let mut b_best = usize::MAX;
    let mut a_keys = 0usize;
    let mut b_keys = 0usize;
    let mut b_matemap = 0usize;
    let mut a_time = f64::MAX;
    let mut b_time = f64::MAX;

    for rep in 0..reps {
        // ---- Arm A: separate mates, flags=0 (production) ----
        let t = Instant::now();
        let a1 = encode_set(&r1_ref, None)?;
        let a2 = encode_set(&r2_ref, None)?;
        let a_dt = t.elapsed().as_secs_f64();
        let a_total = a1.total() + a2.total();

        // ---- Arm B: shared mate blob + FQZ_FREAD2 + mate map ----
        let t = Instant::now();
        let b = encode_set(&combined, Some(&flags))?;
        let mm = mate_map_bsc(n_total, n_r1)?;
        let b_dt = t.elapsed().as_secs_f64();
        let b_total = b.total() + mm.len();

        // ---- Lossless roundtrip (rep 0 only; deterministic) ----
        if rep == 0 {
            let d1 = decode_set(&a1)?;
            let d2 = decode_set(&a2)?;
            if d1 != r1 || d2 != r2 {
                bail!("Arm A roundtrip MISMATCH");
            }
            let db = decode_set(&b)?;
            if db[..n_r1] != r1[..] || db[n_r1..] != r2[..] {
                bail!("Arm B roundtrip MISMATCH");
            }
            eprintln!("roundtrip: both arms lossless ✓\n");
        }

        eprintln!(
            "rep {}: A {} B {} ({:+.3}%)  [A {:.2}s, B {:.2}s]",
            rep,
            a_total,
            b_total,
            (a_total as f64 - b_total as f64) / a_total as f64 * 100.0,
            a_dt,
            b_dt
        );

        if a_total < a_best {
            a_best = a_total;
            a_keys = a1.keys_bsc.len() + a2.keys_bsc.len();
        }
        if b_total < b_best {
            b_best = b_total;
            b_keys = b.keys_bsc.len();
            b_matemap = mm.len();
        }
        a_time = a_time.min(a_dt);
        b_time = b_time.min(b_dt);
    }

    let saving = (a_best as f64 - b_best as f64) / a_best as f64 * 100.0;
    println!("\n================ fqz flags probe ================");
    println!("pairs                : {n_r1}  ({n_total} reads)");
    println!("raw quality bytes    : {raw_bytes}");
    println!("--- Arm A (production: separate mates, flags=0) ---");
    println!("  total bytes        : {a_best}   (incl. sort-keys {a_keys})");
    println!("  best wall          : {a_time:.2}s");
    println!("--- Arm B (shared mate blob + FQZ_FREAD2) ---");
    println!("  total bytes        : {b_best}   (incl. sort-keys {b_keys}, mate-map {b_matemap})");
    println!("  best wall          : {b_time:.2}s");
    println!("------------------------------------------------");
    println!("saving (A-B)/A       : {saving:+.3}%");
    let gate = if saving >= 1.0 {
        "PASS (>=1.0% — proceed to reference-mode integration)"
    } else if saving > 0.0 {
        "MARGINAL (0-1% — ask before building)"
    } else {
        "FAIL (<=0 — stop, record negative result)"
    };
    println!("gate                 : {gate}");
    println!("================================================");
    Ok(())
}
