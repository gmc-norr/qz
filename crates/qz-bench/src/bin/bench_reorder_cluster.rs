//! Reorder-cluster decision-gate: minimizer-sort a FASTQ so coverage-overlapping
//! reads become adjacent, WITHOUT diff/assembly and WITHOUT a permutation stream
//! (order is DROPPED). Output is a pure byte-range permutation of the input records
//! (identical record bytes, reordered) so an original-vs-clustered qz comparison
//! differs ONLY in order. Clustering key = qz's own `compute_min_syncmer_hash`
//! (k=21,s=10 open syncmer, canonical) — faithful to what a real mode would use.
//!
//! Usage:
//!   bench_reorder_cluster <in1.fastq> <out1.fastq>                      (single-end)
//!   bench_reorder_cluster <in1.fastq> <out1.fastq> <in2> <out2>         (paired, key=R1, pairs-as-units)
//!
//! Inputs must be PLAIN (uncompressed) FASTQ. Holds the file(s) in RAM.

use rayon::prelude::*;
use std::fs;
use std::io::{BufWriter, Write};
use std::time::Instant;

use qz_lib::compression::dna_utils::compute_min_syncmer_hash;

#[derive(Clone, Copy)]
struct Rec {
    start: usize,
    len: usize,
    seq_start: usize,
    seq_len: usize,
}

/// Index every 4-line FASTQ record as a byte range + the seq-line range.
fn parse(buf: &[u8]) -> Vec<Rec> {
    let n = buf.len();
    let find = |from: usize| -> Option<usize> {
        buf[from..].iter().position(|&b| b == b'\n').map(|i| from + i)
    };
    let mut recs = Vec::with_capacity(n / 300 + 16);
    let mut p = 0usize;
    while p < n {
        let start = p;
        let nl1 = match find(p) { Some(x) => x, None => break };
        let seq_start = nl1 + 1;
        let nl2 = match find(seq_start) { Some(x) => x, None => break };
        let seq_len = nl2 - seq_start;
        let nl3 = match find(nl2 + 1) { Some(x) => x, None => break };
        let nl4 = match find(nl3 + 1) { Some(x) => x, None => break };
        let end = nl4 + 1;
        recs.push(Rec { start, len: end - start, seq_start, seq_len });
        p = end;
    }
    recs
}

fn keys_for(buf: &[u8], recs: &[Rec]) -> Vec<u64> {
    recs.par_iter()
        .map(|r| compute_min_syncmer_hash(&buf[r.seq_start..r.seq_start + r.seq_len]))
        .collect()
}

/// Sorted record order. Primary key = minimizer (buckets co-locus reads). Tiebreak:
/// original index (default) OR, when `seqsort` (QZ_SEQSORT set, experiment #4), the
/// full read SEQUENCE — so near-identical reads within a minimizer bucket become
/// exactly adjacent, maximizing long-range LZ matches / BWT runs.
fn cluster_order(keys: &[u64], buf: &[u8], recs: &[Rec], seqsort: bool) -> Vec<u32> {
    let mut order: Vec<u32> = (0..keys.len() as u32).collect();
    order.par_sort_unstable_by(|&a, &b| {
        let c = keys[a as usize].cmp(&keys[b as usize]);
        if c != std::cmp::Ordering::Equal {
            return c;
        }
        if seqsort {
            let ra = &recs[a as usize];
            let rb = &recs[b as usize];
            buf[ra.seq_start..ra.seq_start + ra.seq_len]
                .cmp(&buf[rb.seq_start..rb.seq_start + rb.seq_len])
        } else {
            a.cmp(&b)
        }
    });
    order
}

fn write_reordered(path: &str, buf: &[u8], recs: &[Rec], order: &[u32]) -> std::io::Result<usize> {
    let f = fs::File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 22, f);
    let mut written = 0usize;
    for &i in order {
        let r = &recs[i as usize];
        w.write_all(&buf[r.start..r.start + r.len])?;
        written += r.len;
    }
    w.flush()?;
    Ok(written)
}

fn bucket_stats(keys: &[u64], order: &[u32]) -> (usize, usize, usize) {
    // distinct keys, max bucket size, count of no-key reads — order is key-sorted.
    let mut distinct = 0usize;
    let mut max_bucket = 0usize;
    let mut cur = 0usize;
    let mut prev: Option<u64> = None;
    let mut nokey = 0usize;
    for &i in order {
        let k = keys[i as usize];
        if k == u64::MAX {
            nokey += 1;
        }
        match prev {
            Some(p) if p == k => cur += 1,
            _ => {
                if cur > max_bucket {
                    max_bucket = cur;
                }
                distinct += 1;
                cur = 1;
                prev = Some(k);
            }
        }
    }
    if cur > max_bucket {
        max_bucket = cur;
    }
    (distinct, max_bucket, nokey)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let paired = args.len() == 5;
    if !(args.len() == 3 || paired) {
        eprintln!("usage: bench_reorder_cluster <in1> <out1> [<in2> <out2>]");
        std::process::exit(2);
    }

    let t0 = Instant::now();
    let buf1 = fs::read(&args[1]).expect("read in1");
    let recs1 = parse(&buf1);
    eprintln!(
        "[read] in1 {} reads, {:.2} GB, {:.1}s",
        recs1.len(),
        buf1.len() as f64 / 1e9,
        t0.elapsed().as_secs_f64()
    );

    let seqsort = std::env::var("QZ_SEQSORT").is_ok();
    let tk = Instant::now();
    let keys = keys_for(&buf1, &recs1);
    let order = cluster_order(&keys, &buf1, &recs1, seqsort);
    let (distinct, max_bucket, nokey) = bucket_stats(&keys, &order);
    eprintln!(
        "[cluster] key=R1 min-syncmer seqsort={}  distinct_keys={} max_bucket={} no_key={} ({:.1}s)",
        seqsort, distinct, max_bucket, nokey,
        tk.elapsed().as_secs_f64()
    );

    let tw = Instant::now();
    let w1 = write_reordered(&args[2], &buf1, &recs1, &order).expect("write out1");
    let in1_bytes: usize = recs1.iter().map(|r| r.len).sum();
    assert_eq!(w1, in1_bytes, "out1 byte count != sum of input record bytes (not a permutation!)");
    eprintln!("[write] out1 {} reads, {} bytes ({:.1}s)", recs1.len(), w1, tw.elapsed().as_secs_f64());

    if paired {
        let buf2 = fs::read(&args[3]).expect("read in2");
        let recs2 = parse(&buf2);
        assert_eq!(
            recs1.len(),
            recs2.len(),
            "paired: R1 and R2 record counts differ ({} vs {})",
            recs1.len(),
            recs2.len()
        );
        // pairs-as-units: R2 follows the SAME order computed from R1 keys.
        let w2 = write_reordered(&args[4], &buf2, &recs2, &order).expect("write out2");
        let in2_bytes: usize = recs2.iter().map(|r| r.len).sum();
        assert_eq!(w2, in2_bytes, "out2 byte count != sum of input record bytes");
        eprintln!("[write] out2 {} reads, {} bytes (pairs-as-units, R1-keyed)", recs2.len(), w2);
    }

    eprintln!("[done] total {:.1}s", t0.elapsed().as_secs_f64());
}
