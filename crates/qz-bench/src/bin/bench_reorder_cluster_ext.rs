//! Memory-bounded streaming version of the reorder-cluster gate, for whole-WGS
//! inputs that don't fit the in-RAM sort (see bench_reorder_cluster for the simple
//! in-RAM one). Streams a (possibly gzipped) FASTQ, buckets each record by the
//! HIGH bits of its minimizer key into temp files, then sorts each bucket in RAM
//! and concatenates. Because the bucket index is monotonic in the key and each
//! bucket is internally key-sorted, the output is the SAME global minimizer sort
//! as the in-RAM version — at ~1 GB peak instead of ~188 GB. Order is DROPPED.
//!
//! Key = qz's own `compute_min_syncmer_hash` (k=21,s=10). Output is a byte-range
//! permutation of the input records (raw 4-line records copied verbatim).
//!
//! Usage: bench_reorder_cluster_ext <in.fastq[.gz]> <out.fastq> <tmpdir>

use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::time::Instant;

use qz_lib::compression::dna_utils::compute_min_syncmer_hash;

const NBUCKETS: usize = 4096;
const SAMPLE_READS: usize = 8_000_000;

/// Estimate NBUCKETS-1 quantile splitters from a prefix sample of the input so each
/// bucket is a BALANCED key RANGE (not raw high bits, which pile up because the key is
/// a MIN syncmer hash biased small; not a hash, which kills key-order locality). This
/// gives balanced buckets AND globally key-sorted output -> bounded RAM AND best ratio.
/// Prefix sampling is representative for random/flowcell-order input (our case).
fn compute_splitters(path: &str, parallel: bool) -> Vec<u64> {
    let mut reader = open_maybe_gz(path);
    let mut rec = Vec::with_capacity(512);
    let mut keys: Vec<u64>;
    if parallel {
        // A1: collect the sample seqs serially (read is cheap), hash them in parallel.
        let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(SAMPLE_READS.min(1 << 20));
        while seqs.len() < SAMPLE_READS {
            match read_record(&mut reader, &mut rec) {
                Some((ss, se)) => seqs.push(rec[ss..se].to_vec()),
                None => break,
            }
        }
        keys = seqs.par_iter().map(|s| compute_min_syncmer_hash(s)).collect();
    } else {
        keys = Vec::with_capacity(SAMPLE_READS.min(1 << 20));
        while keys.len() < SAMPLE_READS {
            match read_record(&mut reader, &mut rec) {
                Some((ss, se)) => keys.push(compute_min_syncmer_hash(&rec[ss..se])),
                None => break,
            }
        }
    }
    keys.par_sort_unstable();
    let n = keys.len().max(1);
    // splitter[i-1] = key at rank i/NBUCKETS, for i in 1..NBUCKETS
    (1..NBUCKETS).map(|i| keys[(i * n / NBUCKETS).min(n - 1)]).collect()
}

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

/// Read one 4-line FASTQ record into `rec` (cleared first), returning the seq-line
/// byte range within `rec` (excluding its newline), or None at EOF.
fn read_record(r: &mut dyn BufRead, rec: &mut Vec<u8>) -> Option<(usize, usize)> {
    rec.clear();
    let l1 = r.read_until(b'\n', rec).ok()?;
    if l1 == 0 {
        return None;
    }
    let seq_start = rec.len();
    let l2 = r.read_until(b'\n', rec).ok()?;
    if l2 == 0 {
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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: bench_reorder_cluster_ext <in.fastq[.gz]> <out.fastq> <tmpdir>");
        std::process::exit(2);
    }
    let (inp, outp, tmp) = (&args[1], &args[2], &args[3]);
    let _ = fs::remove_dir_all(tmp);
    fs::create_dir_all(tmp).expect("mkdir tmp");
    // A1: parallel per-read minimizer hashing (default). QZ_SERIAL_HASH=1 reverts to
    // the old single-threaded hash for the byte-identical A/B.
    let parallel = std::env::var("QZ_SERIAL_HASH").map(|v| v != "1").unwrap_or(true);
    let t0 = Instant::now();

    // ---- Sample pass: balanced quantile splitters (key-ordered) ----
    let splitters = compute_splitters(inp, parallel);
    eprintln!(
        "[sample] {} splitters from prefix sample ({:.0}s, {})",
        splitters.len(), t0.elapsed().as_secs_f64(), if parallel { "par-hash" } else { "serial-hash" }
    );

    // ---- Pass 1: stream -> buckets by key range (balanced + key-ordered) ----
    let mut writers: Vec<BufWriter<File>> = (0..NBUCKETS)
        .map(|i| BufWriter::with_capacity(1 << 18, File::create(format!("{tmp}/b{i:04}")).unwrap()))
        .collect();
    let mut reader = open_maybe_gz(inp);
    let mut rec = Vec::with_capacity(512);
    let mut nreads: u64 = 0;
    let mut nokey: u64 = 0;
    // Batch records into an arena, hash the batch in parallel, then route serially to
    // bucket files. Bounded RAM (one BATCH-record arena resident); byte-identical output
    // (same key -> same bucket -> same pass2 sort) regardless of hash parallelism.
    const BATCH: usize = 1_000_000;
    let mut arena: Vec<u8> = Vec::with_capacity(BATCH * 320);
    let mut meta: Vec<(usize, usize, usize, usize)> = Vec::with_capacity(BATCH); // (start, rec_len, ss, se)
    loop {
        arena.clear();
        meta.clear();
        while meta.len() < BATCH {
            match read_record(&mut reader, &mut rec) {
                Some((ss, se)) => {
                    let start = arena.len();
                    arena.extend_from_slice(&rec);
                    meta.push((start, rec.len(), ss, se));
                }
                None => break,
            }
        }
        if meta.is_empty() {
            break;
        }
        let keys: Vec<u64> = if parallel {
            meta.par_iter().map(|&(st, _ln, ss, se)| compute_min_syncmer_hash(&arena[st + ss..st + se])).collect()
        } else {
            meta.iter().map(|&(st, _ln, ss, se)| compute_min_syncmer_hash(&arena[st + ss..st + se])).collect()
        };
        for (i, &(st, ln, _ss, _se)) in meta.iter().enumerate() {
            let key = keys[i];
            if key == u64::MAX {
                nokey += 1;
            }
            let bucket = splitters.partition_point(|&s| s <= key);
            let w = &mut writers[bucket];
            w.write_all(&key.to_le_bytes()).unwrap();
            w.write_all(&(ln as u32).to_le_bytes()).unwrap();
            w.write_all(&arena[st..st + ln]).unwrap();
            nreads += 1;
        }
    }
    for w in &mut writers {
        w.flush().unwrap();
    }
    drop(writers);
    eprintln!(
        "[pass1] {} reads streamed to {} buckets (no_key={}) {:.0}s",
        nreads, NBUCKETS, nokey, t0.elapsed().as_secs_f64()
    );

    // ---- Pass 2: sort each bucket by full key, concatenate in bucket order ----
    let t2 = Instant::now();
    let out = File::create(outp).expect("create out");
    let mut out = BufWriter::with_capacity(1 << 22, out);
    let mut written: u64 = 0;
    let mut max_bucket = 0u64;
    let mut total_out: u64 = 0;
    for i in 0..NBUCKETS {
        let data = fs::read(format!("{tmp}/b{i:04}")).unwrap();
        // parse [key u64][len u32][record] entries
        let mut ents: Vec<(u64, usize, usize)> = Vec::new();
        let mut p = 0usize;
        while p + 12 <= data.len() {
            let key = u64::from_le_bytes(data[p..p + 8].try_into().unwrap());
            let len = u32::from_le_bytes(data[p + 8..p + 12].try_into().unwrap()) as usize;
            let off = p + 12;
            ents.push((key, off, len));
            p = off + len;
        }
        if ents.len() as u64 > max_bucket {
            max_bucket = ents.len() as u64;
        }
        ents.par_sort_unstable_by_key(|e| e.0);
        for (_, off, len) in &ents {
            out.write_all(&data[*off..*off + *len]).unwrap();
            total_out += *len as u64;
        }
        written += ents.len() as u64;
    }
    out.flush().unwrap();
    let _ = fs::remove_dir_all(tmp);
    assert_eq!(written, nreads, "pass2 wrote {written} records, pass1 read {nreads}");
    eprintln!(
        "[pass2] sorted+wrote {} reads, {} out bytes, max_bucket={} ({:.0}s)",
        written, total_out, max_bucket, t2.elapsed().as_secs_f64()
    );
    eprintln!("[done] total {:.0}s", t0.elapsed().as_secs_f64());
}
