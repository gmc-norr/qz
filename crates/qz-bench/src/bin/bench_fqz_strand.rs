//! Probe: does feeding fqzcomp quality in SEQUENCING orientation (what qz's
//! FASTQ reference mode already does) beat REFERENCE orientation (what a BAM
//! stores, and what bz had to un-reverse to gain -3.3%)?
//!
//! If sequencing-order wins, qz is already at the good end and the "strand
//! reorient" lever is moot for qz. If reference-order wins, qz is leaving
//! something on the table.
//!
//! Inputs are two files, one quality string per line:
//!   <ref_order.txt>  — QUAL exactly as stored in the BAM (reference orientation)
//!   <seq_order.txt>  — QUAL un-reversed for reverse-strand reads (sequencing order)
//! Both pass through the EXACT production framing (mean-quality sort + 500K
//! sub-chunks + fqzcomp + BSC'd sort-key side stream).
//!
//! Usage: bench_fqz_strand <ref_order.txt> <seq_order.txt> [num_reads] [reps]

use anyhow::{Context, Result, bail};
use qz_lib::compression::{bsc, fqzcomp};
use rayon::prelude::*;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const FQZCOMP_SUB_CHUNK: usize = 500_000;
const DEFAULT_READS: usize = 1_000_000;
const DEFAULT_REPS: usize = 3;

fn mean_key(q: &[u8]) -> u8 {
    let sum: u64 = q.iter().map(|&b| b as u64).sum();
    (sum / q.len().max(1) as u64) as u8
}

fn read_lines(path: &str, n: usize) -> Result<Vec<Vec<u8>>> {
    let f = std::fs::File::open(path).with_context(|| format!("open {path}"))?;
    let mut r = BufReader::with_capacity(1 << 20, f);
    let mut out = Vec::with_capacity(n);
    let mut line = Vec::new();
    while out.len() < n {
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }
        out.push(std::mem::take(&mut line));
    }
    Ok(out)
}

struct Encoded {
    blobs: Vec<Vec<u8>>,
    keys_bsc: Vec<u8>,
    n: usize,
}
impl Encoded {
    fn total(&self) -> usize {
        self.blobs.iter().map(|b| b.len()).sum::<usize>() + self.keys_bsc.len()
    }
}

/// Production-faithful: mean-quality stable sort, 500K sub-chunks, fqzcomp each,
/// BSC the per-read keys (the permutation side channel).
fn encode_set(quals: &[&[u8]]) -> Result<Encoded> {
    let n = quals.len();
    let keys: Vec<u8> = quals.par_iter().map(|q| mean_key(q)).collect();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by_key(|&i| keys[i]);
    let sorted: Vec<&[u8]> = idx.iter().map(|&i| quals[i]).collect();
    let starts: Vec<usize> = (0..n).step_by(FQZCOMP_SUB_CHUNK).collect();
    let (blobs, keys_bsc) = rayon::join(
        || -> Result<Vec<Vec<u8>>> {
            starts
                .par_iter()
                .map(|&s| {
                    let e = (s + FQZCOMP_SUB_CHUNK).min(n);
                    fqzcomp::compress(&sorted[s..e], 0)
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

fn decode_set(enc: &Encoded, original: &[Vec<u8>]) -> Result<()> {
    let keys = bsc::decompress_parallel(&enc.keys_bsc)?;
    let mut idx: Vec<usize> = (0..enc.n).collect();
    idx.sort_by_key(|&i| keys[i]);
    let mut sorted: Vec<Vec<u8>> = Vec::with_capacity(enc.n);
    let mut produced = 0usize;
    for blob in &enc.blobs {
        let this = (enc.n - produced).min(FQZCOMP_SUB_CHUNK);
        let (concat, lengths) = fqzcomp::decompress(blob, this)?;
        let mut off = 0usize;
        for &l in &lengths {
            let l = l as usize;
            sorted.push(concat[off..off + l].to_vec());
            off += l;
        }
        produced += this;
    }
    let mut restored = vec![Vec::new(); enc.n];
    for (k, &orig_i) in idx.iter().enumerate() {
        restored[orig_i] = std::mem::take(&mut sorted[k]);
    }
    if restored != original {
        bail!("roundtrip MISMATCH");
    }
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: {} <ref_order.txt> <seq_order.txt> [num_reads] [reps]",
            args[0]
        );
        std::process::exit(2);
    }
    let n = args
        .get(3)
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_READS);
    let reps = args
        .get(4)
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_REPS);

    eprintln!("Reading up to {n} quality lines from each file ...");
    let (refo, seqo) = rayon::join(|| read_lines(&args[1], n), || read_lines(&args[2], n));
    let refo = refo?;
    let seqo = seqo?;
    if refo.len() != seqo.len() {
        bail!("file lengths differ: {} vs {}", refo.len(), seqo.len());
    }
    eprintln!("Loaded {} reads each\n", refo.len());

    let refo_ref: Vec<&[u8]> = refo.iter().map(|v| v.as_slice()).collect();
    let seqo_ref: Vec<&[u8]> = seqo.iter().map(|v| v.as_slice()).collect();

    let mut ref_best = usize::MAX;
    let mut seq_best = usize::MAX;
    for rep in 0..reps {
        let t = Instant::now();
        let er = encode_set(&refo_ref)?;
        let dr = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let es = encode_set(&seqo_ref)?;
        let ds = t.elapsed().as_secs_f64();

        if rep == 0 {
            decode_set(&er, &refo).context("ref-order roundtrip")?;
            decode_set(&es, &seqo).context("seq-order roundtrip")?;
            eprintln!("roundtrip: both orientations lossless ✓\n");
        }
        let (rt, st) = (er.total(), es.total());
        eprintln!(
            "rep {rep}: ref-order {rt}  seq-order {st}  (seq saves {:+.3}%)  [{dr:.2}s / {ds:.2}s]",
            (rt as f64 - st as f64) / rt as f64 * 100.0
        );
        ref_best = ref_best.min(rt);
        seq_best = seq_best.min(st);
    }

    let saving = (ref_best as f64 - seq_best as f64) / ref_best as f64 * 100.0;
    println!("\n============ fqz strand-orientation probe ============");
    println!("reads                 : {}", refo.len());
    println!("reference-order bytes : {ref_best}   (BAM convention; bz's pre-reorient input)");
    println!("sequencing-order bytes: {seq_best}   (qz FASTQ convention)");
    println!("seq-order saving      : {saving:+.3}%");
    println!("-----------------------------------------------------");
    if saving > 0.5 {
        println!("=> sequencing-order WINS by {saving:.2}%. qz already feeds this");
        println!("   orientation, so the strand-reorient lever is ALREADY captured");
        println!("   (the bz -3.3% recovered exactly this; qz never lost it).");
    } else if saving.abs() <= 0.5 {
        println!("=> orientation is ~neutral on this data ({saving:+.2}%).");
    } else {
        println!("=> reference-order wins ({saving:+.2}%) — unexpected; investigate.");
    }
    println!("=====================================================");
    Ok(())
}
