//! LEVER A: fused reorder->BSC with the `flush_bsc` step DECOUPLED onto a bounded
//! worker pool. Identical to `bench_reorder_cluster_fused_bsc` except that pass-2
//! flushes don't compress inline (blocking the bucket-drain gather loop ~1.7s per
//! 256 MB); instead the full buffer is swapped out and handed to a bounded channel
//! feeding N BSC workers, so 3-4 flushes compress CONCURRENTLY while the gather
//! loop keeps draining buckets. §8k/§8m: the in-pipeline BSC stage was small-flush
//! work-starved at ~149 MB/s; the concurrency gate measured N=4 -> 335 MB/s.
//!
//! Buffer boundaries are byte-identical to the sync bench (same FLUSH trigger after
//! each append), so each 256 MB buffer compresses to the SAME bytes -> bsc_seq_total
//! is byte-identical; only the WALL (pass2_s / total_s) drops. A/B against the sync
//! bench on the same input proves the speedup AND the byte-identicality.
//!
//! Env: QZ_POOL_WORKERS (default 4), QZ_POOL_DEPTH (channel slack, default 2).
//! Usage: bench_reorder_cluster_fused_pool <R1[.gz]> <R2[.gz]> <tmpdir>

use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::sync_channel;
use std::sync::Arc;
use std::time::Instant;

use qz_lib::compression::bsc;
use qz_lib::compression::dna_utils::compute_min_syncmer_hash;

const NBUCKETS: usize = 4096;
const SAMPLE_READS: usize = 8_000_000;
const FLUSH: usize = 256 * 1024 * 1024;

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

fn read_seq_into(r: &mut dyn BufRead, buf: &mut Vec<u8>, scratch: &mut Vec<u8>) -> Option<(usize, usize)> {
    scratch.clear();
    if r.read_until(b'\n', scratch).ok()? == 0 {
        return None;
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
        return None;
    }
    scratch.clear();
    if r.read_until(b'\n', scratch).ok()? == 0 {
        return None;
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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: bench_reorder_cluster_fused_pool <R1[.gz]> <R2[.gz]> <tmpdir>");
        std::process::exit(2);
    }
    let (in1, in2, tmp) = (&args[1], &args[2], &args[3]);
    let workers: usize = std::env::var("QZ_POOL_WORKERS").ok().and_then(|s| s.parse().ok()).unwrap_or(4);
    let depth: usize = std::env::var("QZ_POOL_DEPTH").ok().and_then(|s| s.parse().ok()).unwrap_or(2);
    let _ = fs::remove_dir_all(tmp);
    fs::create_dir_all(tmp).expect("mkdir tmp");
    let t0 = Instant::now();

    let splitters = compute_splitters(in1);
    eprintln!("[sample] {:.0}s", t0.elapsed().as_secs_f64());

    // ---- Pass 1: identical to sync bench ----
    let mut writers: Vec<BufWriter<File>> = (0..NBUCKETS)
        .map(|i| BufWriter::with_capacity(1 << 18, File::create(format!("{tmp}/b{i:04}")).unwrap()))
        .collect();
    let mut r1 = open_maybe_gz(in1);
    let mut r2 = open_maybe_gz(in2);
    const BATCH: usize = 400_000;
    let mut buf1: Vec<u8> = Vec::with_capacity(BATCH * 160);
    let mut buf2: Vec<u8> = Vec::with_capacity(BATCH * 160);
    let mut off1: Vec<(usize, usize)> = Vec::with_capacity(BATCH);
    let mut off2: Vec<(usize, usize)> = Vec::with_capacity(BATCH);
    let mut scratch = Vec::with_capacity(512);
    let mut idx: u64 = 0;
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
        let th = Instant::now();
        let keys: Vec<u64> = off1[..got]
            .par_iter()
            .map(|&(s, e)| compute_min_syncmer_hash(&buf1[s..e]))
            .collect();
        hash_s += th.elapsed().as_secs_f64();
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

    // ---- Pass 2: bucket sort + append, BSC via bounded WORKER POOL ----
    // Channel carries (stream_id, full 256 MB buffer). Bounded at `workers + depth`
    // so peak in-flight buffers (and thus RAM) is capped; the gather loop blocks on
    // send when full -> backpressure. Workers compress concurrently on the rayon pool.
    let t2 = Instant::now();
    let bsc1 = Arc::new(AtomicU64::new(0));
    let bsc2 = Arc::new(AtomicU64::new(0));
    let busy_ns = Arc::new(AtomicU64::new(0)); // sum of per-flush compress time (diag)
    let (tx, rx) = sync_channel::<(u8, Vec<u8>)>(workers + depth);
    let rx = Arc::new(std::sync::Mutex::new(rx));
    let mut handles = Vec::with_capacity(workers);
    for _ in 0..workers {
        let rx = Arc::clone(&rx);
        let bsc1 = Arc::clone(&bsc1);
        let bsc2 = Arc::clone(&bsc2);
        let busy_ns = Arc::clone(&busy_ns);
        handles.push(std::thread::spawn(move || loop {
            let msg = {
                let lock = rx.lock().unwrap();
                lock.recv()
            };
            let (sid, buf) = match msg {
                Ok(m) => m,
                Err(_) => break, // all senders dropped
            };
            let t = Instant::now();
            let c = bsc::compress_parallel(&buf).expect("bsc");
            busy_ns.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
            if sid == 1 {
                bsc1.fetch_add(c.len() as u64, Ordering::Relaxed);
            } else {
                bsc2.fetch_add(c.len() as u64, Ordering::Relaxed);
            }
        }));
    }

    let mut buf1: Vec<u8> = Vec::with_capacity(FLUSH + (1 << 20));
    let mut buf2: Vec<u8> = Vec::with_capacity(FLUSH + (1 << 20));
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
                let full = std::mem::replace(&mut buf1, Vec::with_capacity(FLUSH + (1 << 20)));
                tx.send((1, full)).unwrap();
            }
            if buf2.len() >= FLUSH {
                let full = std::mem::replace(&mut buf2, Vec::with_capacity(FLUSH + (1 << 20)));
                tx.send((2, full)).unwrap();
            }
        }
        written += ents.len() as u64;
        let _ = fs::remove_file(&fp);
    }
    // final partial flushes
    if !buf1.is_empty() {
        tx.send((1, std::mem::take(&mut buf1))).unwrap();
    }
    if !buf2.is_empty() {
        tx.send((2, std::mem::take(&mut buf2))).unwrap();
    }
    drop(tx);
    for h in handles {
        h.join().unwrap();
    }
    let _ = fs::remove_dir_all(tmp);
    assert_eq!(written, npairs);
    let bsc1 = bsc1.load(Ordering::Relaxed);
    let bsc2 = bsc2.load(Ordering::Relaxed);
    let bsc_busy = busy_ns.load(Ordering::Relaxed) as f64 / 1e9;
    let pass2 = t2.elapsed().as_secs_f64();
    let total = t0.elapsed().as_secs_f64();
    eprintln!(
        "[pass2] sort {sort_s:.0}s, BSC pool(w={workers}) busy-sum {bsc_busy:.0}s overlapped, {pass2:.0}s"
    );
    eprintln!("[done] max_bucket={max_bucket}  total {total:.0}s");
    println!(
        "FUSEDPOOL npairs={npairs} workers={workers} bsc_seq_R1={bsc1} bsc_seq_R2={bsc2} bsc_seq_total={} pass1_s={pass1:.0} pass2_s={pass2:.0} total_s={total:.0} sort_s={sort_s:.0} bsc_busy_s={bsc_busy:.0}",
        bsc1 + bsc2
    );
}
