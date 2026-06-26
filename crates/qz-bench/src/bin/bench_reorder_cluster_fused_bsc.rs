//! FUSED reorder->BSC prototype (task 2): cluster paired reads (pairs-as-units) AND
//! BSC-compress the SEQUENCE stream INLINE as buckets drain — NO clustered-plain
//! intermediate on disk, bounded RAM. Proves the production fused path: the BSC
//! compress overlaps/ hides under the (already-required) bucket-sort I/O, peak RAM is
//! the worst single bucket + the flush buffers (NOT the dataset), and the 2× 328 GB
//! materialize+re-read of the naive (reorder-to-plain → qz-reads-plain) path is gone.
//!
//! Seq-only: the sequence stream is the bulk (~143 GB of the 328 GB) and the ratio
//! lever; headers/quals fuse identically through the SAME permutation (extra streams).
//! Key = qz's `compute_min_syncmer_hash`. Output order byte-identical to pe_ext.
//!
//! Usage: bench_reorder_cluster_fused_bsc <R1[.gz]> <R2[.gz]> <tmpdir>

use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::time::Instant;

use qz_lib::compression::bsc;
use qz_lib::compression::dna_utils::compute_min_syncmer_hash;

const NBUCKETS: usize = 4096;
const SAMPLE_READS: usize = 8_000_000;
const FLUSH: usize = 256 * 1024 * 1024; // accumulate ~256 MB of clustered seq, then BSC it

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

/// Read one 4-line record into `rec`; return seq-line byte range (no newline) or None.
fn read_record(r: &mut dyn BufRead, rec: &mut Vec<u8>) -> Option<(usize, usize)> {
    rec.clear();
    if r.read_until(b'\n', rec).ok()? == 0 {
        return None;
    }
    let ss = rec.len();
    if r.read_until(b'\n', rec).ok()? == 0 {
        return None;
    }
    let se = rec.len() - if rec.last() == Some(&b'\n') { 1 } else { 0 };
    if r.read_until(b'\n', rec).ok()? == 0 {
        return None;
    }
    if r.read_until(b'\n', rec).ok()? == 0 {
        return None;
    }
    Some((ss, se))
}

/// Append ONLY the seq line of one 4-line record into `buf` (header/+/qual discarded via
/// `scratch`); return the seq byte range within `buf` (newline stripped), or None at EOF.
/// Used by the parallel-hash pass-1 batch reader (seq-only arena).
fn read_seq_into(r: &mut dyn BufRead, buf: &mut Vec<u8>, scratch: &mut Vec<u8>) -> Option<(usize, usize)> {
    scratch.clear();
    if r.read_until(b'\n', scratch).ok()? == 0 {
        return None; // header
    }
    let s = buf.len();
    if r.read_until(b'\n', buf).ok()? == 0 {
        return None;
    }
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    let e = buf.len();
    scratch.clear();
    if r.read_until(b'\n', scratch).ok()? == 0 {
        return None; // +
    }
    scratch.clear();
    if r.read_until(b'\n', scratch).ok()? == 0 {
        return None; // qual
    }
    Some((s, e))
}

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

/// Flush a clustered-seq accumulator through BSC (parallel 25 MB blocks), return bytes,
/// reset the buffer. Mirrors what a fused producer streams into the block compressor.
fn flush_bsc(buf: &mut Vec<u8>, total: &mut u64, secs: &mut f64) {
    if buf.is_empty() {
        return;
    }
    let t = Instant::now();
    let c = bsc::compress_parallel(buf).expect("bsc");
    *secs += t.elapsed().as_secs_f64();
    *total += c.len() as u64;
    buf.clear();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: bench_reorder_cluster_fused_bsc <R1[.gz]> <R2[.gz]> <tmpdir>");
        std::process::exit(2);
    }
    let (in1, in2, tmp) = (&args[1], &args[2], &args[3]);
    let _ = fs::remove_dir_all(tmp);
    fs::create_dir_all(tmp).expect("mkdir tmp");
    let t0 = Instant::now();

    let splitters = compute_splitters(in1);
    eprintln!("[sample] {:.0}s", t0.elapsed().as_secs_f64());

    // ---- Pass 1: lockstep -> buckets keyed by R1 minimizer, SEQ ONLY ----
    // entry: [key u64][idx u64][len1 u32][seq1][len2 u32][seq2]
    let mut writers: Vec<BufWriter<File>> = (0..NBUCKETS)
        .map(|i| BufWriter::with_capacity(1 << 18, File::create(format!("{tmp}/b{i:04}")).unwrap()))
        .collect();
    // PARALLEL-HASH pass 1: serial-read a BATCH of seqs into two arenas, hash the batch
    // across all cores (rayon), then serial-write to buckets. idx is assigned in global
    // read order so the (key,idx) sort — and thus the output — is unchanged vs serial.
    let mut r1 = open_maybe_gz(in1);
    let mut r2 = open_maybe_gz(in2);
    const BATCH: usize = 400_000;
    let mut buf1: Vec<u8> = Vec::with_capacity(BATCH * 160);
    let mut buf2: Vec<u8> = Vec::with_capacity(BATCH * 160);
    let mut off1: Vec<(usize, usize)> = Vec::with_capacity(BATCH);
    let mut off2: Vec<(usize, usize)> = Vec::with_capacity(BATCH);
    let mut scratch = Vec::with_capacity(512);
    let mut idx: u64 = 0;
    // pass-1 stage profile: read (serial parse into arenas) vs hash (parallel) vs write
    // (serial route+bucket). hash_s is cache/disk-independent -> it answers "would SIMD help?"
    let (mut read_s, mut hash_s, mut write_s) = (0.0f64, 0.0f64, 0.0f64);
    loop {
        let tr = Instant::now();
        buf1.clear();
        buf2.clear();
        off1.clear();
        off2.clear();
        let mut got = 0usize;
        for _ in 0..BATCH {
            match read_seq_into(&mut r1, &mut buf1, &mut scratch) {
                Some(r) => off1.push(r),
                None => break,
            }
            match read_seq_into(&mut r2, &mut buf2, &mut scratch) {
                Some(r) => off2.push(r),
                None => panic!("R1/R2 mismatch at pair {}", idx + got as u64),
            }
            got += 1;
        }
        read_s += tr.elapsed().as_secs_f64();
        if got == 0 {
            break;
        }
        // parallel hash R1 seqs (the ~35-min serial bottleneck, now across all cores)
        let th = Instant::now();
        let keys: Vec<u64> = off1[..got]
            .par_iter()
            .map(|&(s, e)| compute_min_syncmer_hash(&buf1[s..e]))
            .collect();
        hash_s += th.elapsed().as_secs_f64();
        // serial write to buckets
        let tw = Instant::now();
        for i in 0..got {
            let key = keys[i];
            let bucket = splitters.partition_point(|&s| s <= key);
            let (s1, e1) = off1[i];
            let (s2, e2) = off2[i];
            let w = &mut writers[bucket];
            w.write_all(&key.to_le_bytes()).unwrap();
            w.write_all(&idx.to_le_bytes()).unwrap();
            w.write_all(&((e1 - s1) as u32).to_le_bytes()).unwrap();
            w.write_all(&buf1[s1..e1]).unwrap();
            w.write_all(&((e2 - s2) as u32).to_le_bytes()).unwrap();
            w.write_all(&buf2[s2..e2]).unwrap();
            idx += 1;
        }
        write_s += tw.elapsed().as_secs_f64();
    }
    for w in &mut writers {
        w.flush().unwrap();
    }
    drop(writers);
    let npairs = idx;
    let pass1 = t0.elapsed().as_secs_f64();
    eprintln!(
        "[pass1] {npairs} pairs (PARALLEL hash) {pass1:.0}s | read(parse) {read_s:.0}s  hash {hash_s:.0}s  write(route+bucket) {write_s:.0}s",
    );
    if std::env::var("QZ_PASS1_ONLY").is_ok() {
        let _ = fs::remove_dir_all(tmp);
        println!("PASS1_ONLY npairs={npairs} pass1_s={pass1:.0} read_s={read_s:.0} hash_s={hash_s:.0} write_s={write_s:.0}");
        return;
    }

    // ---- Pass 2: per bucket sort by (key,idx), append seq to flush buffers, BSC inline ----
    let t2 = Instant::now();
    let mut buf1: Vec<u8> = Vec::with_capacity(FLUSH + (1 << 20));
    let mut buf2: Vec<u8> = Vec::with_capacity(FLUSH + (1 << 20));
    let (mut bsc1, mut bsc2): (u64, u64) = (0, 0);
    let (mut bsc1s, mut bsc2s): (f64, f64) = (0.0, 0.0);
    let mut written: u64 = 0;
    let mut max_bucket = 0u64;
    let mut sort_s = 0.0f64;
    for i in 0..NBUCKETS {
        let fp = format!("{tmp}/b{i:04}");
        let data = fs::read(&fp).unwrap();
        let mut ents: Vec<(u64, u64, usize, usize, usize, usize)> = Vec::new();
        let mut p = 0usize;
        while p + 20 <= data.len() {
            let key = u64::from_le_bytes(data[p..p + 8].try_into().unwrap());
            let pidx = u64::from_le_bytes(data[p + 8..p + 16].try_into().unwrap());
            let l1 = u32::from_le_bytes(data[p + 16..p + 20].try_into().unwrap()) as usize;
            let o1 = p + 20;
            let lp2 = o1 + l1;
            let l2 = u32::from_le_bytes(data[lp2..lp2 + 4].try_into().unwrap()) as usize;
            let o2 = lp2 + 4;
            ents.push((key, pidx, o1, l1, o2, l2));
            p = o2 + l2;
        }
        if ents.len() as u64 > max_bucket {
            max_bucket = ents.len() as u64;
        }
        let ts = Instant::now();
        ents.par_sort_unstable_by_key(|e| (e.0, e.1));
        sort_s += ts.elapsed().as_secs_f64();
        for (_, _, o1, l1, o2, l2) in &ents {
            buf1.extend_from_slice(&data[*o1..*o1 + *l1]);
            buf2.extend_from_slice(&data[*o2..*o2 + *l2]);
            if buf1.len() >= FLUSH {
                flush_bsc(&mut buf1, &mut bsc1, &mut bsc1s);
            }
            if buf2.len() >= FLUSH {
                flush_bsc(&mut buf2, &mut bsc2, &mut bsc2s);
            }
        }
        written += ents.len() as u64;
        let _ = fs::remove_file(&fp);
    }
    flush_bsc(&mut buf1, &mut bsc1, &mut bsc1s);
    flush_bsc(&mut buf2, &mut bsc2, &mut bsc2s);
    let _ = fs::remove_dir_all(tmp);
    assert_eq!(written, npairs);
    let pass2 = t2.elapsed().as_secs_f64();
    let total = t0.elapsed().as_secs_f64();
    eprintln!("[pass2] sort {sort_s:.0}s, BSC inline R1 {bsc1s:.0}s R2 {bsc2s:.0}s, {pass2:.0}s");
    eprintln!("[done] max_bucket={max_bucket}  total {total:.0}s");
    println!(
        "FUSED npairs={npairs} bsc_seq_R1={bsc1} bsc_seq_R2={bsc2} bsc_seq_total={} pass1_s={pass1:.0} pass2_s={pass2:.0} total_s={total:.0} sort_s={sort_s:.0} bsc_s={:.0}",
        bsc1 + bsc2,
        bsc1s + bsc2s
    );
}
