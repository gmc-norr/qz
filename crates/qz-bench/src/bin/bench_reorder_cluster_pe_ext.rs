//! Memory-bounded PAIRED-END reorder-cluster for whole-WGS inputs (the paired
//! sibling of bench_reorder_cluster_ext). Streams gzipped R1+R2 in lockstep, keys
//! each PAIR by R1's minimizer, range-partitions into temp buckets storing BOTH
//! mate records in one entry, then sorts each bucket by (key, pair-index) and emits
//! the two mates under the SAME permutation — pairs-as-units, pairing preserved with
//! NO order stream. Output order is byte-identical to the in-RAM bench_reorder_cluster
//! paired mode (same (key, original-index) comparator; range-partition is just a
//! bucketed external sort). Order is DROPPED. Peak RAM = the largest single bucket
//! (~2 GB at WGS), not the dataset.
//!
//! Key = qz's own `compute_min_syncmer_hash` (k=21,s=10).
//!
//! Usage: bench_reorder_cluster_pe_ext <R1.fastq[.gz]> <R2.fastq[.gz]> <outR1> <outR2> <tmpdir>

use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::time::Instant;

use qz_lib::compression::dna_utils::compute_min_syncmer_hash;

const NBUCKETS: usize = 4096;
const SAMPLE_READS: usize = 8_000_000;

fn open_maybe_gz(path: &str) -> Box<dyn BufRead> {
    let mut f = File::open(path).expect("open input");
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic).unwrap_or(0);
    use std::io::Seek;
    f.seek(std::io::SeekFrom::Start(0)).unwrap();
    if n == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        Box::new(BufReader::with_capacity(1 << 22, MultiGzDecoder::new(BufReader::with_capacity(1 << 22, f))))
    } else {
        Box::new(BufReader::with_capacity(1 << 22, f))
    }
}

/// Read one 4-line FASTQ record into `rec` (cleared first); return the seq-line byte
/// range within `rec` (excluding newline), or None at EOF.
fn read_record(r: &mut dyn BufRead, rec: &mut Vec<u8>) -> Option<(usize, usize)> {
    rec.clear();
    if r.read_until(b'\n', rec).ok()? == 0 {
        return None;
    }
    let seq_start = rec.len();
    if r.read_until(b'\n', rec).ok()? == 0 {
        return None;
    }
    let seq_end = rec.len() - if rec.last() == Some(&b'\n') { 1 } else { 0 };
    if r.read_until(b'\n', rec).ok()? == 0 {
        return None;
    }
    if r.read_until(b'\n', rec).ok()? == 0 {
        return None;
    }
    Some((seq_start, seq_end))
}

/// Quantile splitters from a prefix sample of R1 keys (balanced, key-ordered buckets).
fn compute_splitters(path: &str) -> Vec<u64> {
    let mut reader = open_maybe_gz(path);
    let mut rec = Vec::with_capacity(512);
    let mut keys: Vec<u64> = Vec::with_capacity(SAMPLE_READS.min(1 << 20));
    while keys.len() < SAMPLE_READS {
        match read_record(&mut reader, &mut rec) {
            Some((ss, se)) => keys.push(compute_min_syncmer_hash(&rec[ss..se])),
            None => break,
        }
    }
    keys.par_sort_unstable();
    let n = keys.len().max(1);
    (1..NBUCKETS).map(|i| keys[(i * n / NBUCKETS).min(n - 1)]).collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 6 {
        eprintln!("usage: bench_reorder_cluster_pe_ext <R1[.gz]> <R2[.gz]> <outR1> <outR2> <tmpdir>");
        std::process::exit(2);
    }
    let (in1, in2, out1p, out2p, tmp) = (&args[1], &args[2], &args[3], &args[4], &args[5]);
    let _ = fs::remove_dir_all(tmp);
    fs::create_dir_all(tmp).expect("mkdir tmp");
    let t0 = Instant::now();

    let splitters = compute_splitters(in1);
    eprintln!("[sample] {} splitters ({:.0}s)", splitters.len(), t0.elapsed().as_secs_f64());

    // ---- Pass 1: lockstep stream R1+R2 -> buckets keyed by R1 minimizer ----
    // entry layout: [key u64][idx u64][len1 u32][rec1][len2 u32][rec2]
    let mut writers: Vec<BufWriter<File>> = (0..NBUCKETS)
        .map(|i| BufWriter::with_capacity(1 << 18, File::create(format!("{tmp}/b{i:04}")).unwrap()))
        .collect();
    let mut r1 = open_maybe_gz(in1);
    let mut r2 = open_maybe_gz(in2);
    let mut rec1 = Vec::with_capacity(512);
    let mut rec2 = Vec::with_capacity(512);
    let mut idx: u64 = 0;
    let mut nokey: u64 = 0;
    loop {
        let a = read_record(&mut r1, &mut rec1);
        let b = read_record(&mut r2, &mut rec2);
        match (a, b) {
            (Some((ss, se)), Some(_)) => {
                let key = compute_min_syncmer_hash(&rec1[ss..se]);
                if key == u64::MAX {
                    nokey += 1;
                }
                let bucket = splitters.partition_point(|&s| s <= key);
                let w = &mut writers[bucket];
                w.write_all(&key.to_le_bytes()).unwrap();
                w.write_all(&idx.to_le_bytes()).unwrap();
                w.write_all(&(rec1.len() as u32).to_le_bytes()).unwrap();
                w.write_all(&rec1).unwrap();
                w.write_all(&(rec2.len() as u32).to_le_bytes()).unwrap();
                w.write_all(&rec2).unwrap();
                idx += 1;
            }
            (None, None) => break,
            _ => panic!("R1/R2 record-count mismatch at pair {idx}"),
        }
    }
    for w in &mut writers {
        w.flush().unwrap();
    }
    drop(writers);
    let npairs = idx;
    eprintln!("[pass1] {} pairs -> {} buckets (no_key={}) {:.0}s", npairs, NBUCKETS, nokey, t0.elapsed().as_secs_f64());

    // ---- Pass 2: per bucket, sort by (key, idx), emit both mates; delete bucket ----
    let t2 = Instant::now();
    let mut o1 = BufWriter::with_capacity(1 << 22, File::create(out1p).expect("create outR1"));
    let mut o2 = BufWriter::with_capacity(1 << 22, File::create(out2p).expect("create outR2"));
    let mut written: u64 = 0;
    let mut max_bucket = 0u64;
    let (mut b1, mut b2): (u64, u64) = (0, 0);
    for i in 0..NBUCKETS {
        let fp = format!("{tmp}/b{i:04}");
        let data = fs::read(&fp).unwrap();
        // parse entries -> (key, idx, off1, len1, off2, len2)
        let mut ents: Vec<(u64, u64, usize, usize, usize, usize)> = Vec::new();
        let mut p = 0usize;
        while p + 20 <= data.len() {
            let key = u64::from_le_bytes(data[p..p + 8].try_into().unwrap());
            let pidx = u64::from_le_bytes(data[p + 8..p + 16].try_into().unwrap());
            let len1 = u32::from_le_bytes(data[p + 16..p + 20].try_into().unwrap()) as usize;
            let off1 = p + 20;
            let lp2 = off1 + len1;
            let len2 = u32::from_le_bytes(data[lp2..lp2 + 4].try_into().unwrap()) as usize;
            let off2 = lp2 + 4;
            ents.push((key, pidx, off1, len1, off2, len2));
            p = off2 + len2;
        }
        if ents.len() as u64 > max_bucket {
            max_bucket = ents.len() as u64;
        }
        ents.par_sort_unstable_by_key(|e| (e.0, e.1));
        for (_, _, off1, len1, off2, len2) in &ents {
            o1.write_all(&data[*off1..*off1 + *len1]).unwrap();
            o2.write_all(&data[*off2..*off2 + *len2]).unwrap();
            b1 += *len1 as u64;
            b2 += *len2 as u64;
        }
        written += ents.len() as u64;
        let _ = fs::remove_file(&fp); // cap peak disk: free bucket as soon as it's drained
    }
    o1.flush().unwrap();
    o2.flush().unwrap();
    let _ = fs::remove_dir_all(tmp);
    assert_eq!(written, npairs, "pass2 wrote {written} pairs, pass1 read {npairs}");
    eprintln!(
        "[pass2] sorted+wrote {} pairs (R1 {} B, R2 {} B), max_bucket={} ({:.0}s)",
        written, b1, b2, max_bucket, t2.elapsed().as_secs_f64()
    );
    eprintln!("[done] total {:.0}s", t0.elapsed().as_secs_f64());
}
