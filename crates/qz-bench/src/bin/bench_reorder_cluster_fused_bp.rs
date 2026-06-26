//! Block-parse variant of the fused reorder->BSC prototype (task: kill the 732 s
//! serial `read_until` parse, profiled as 75% of pass1). Same clustering + parallel
//! hash + inline BSC as `bench_reorder_cluster_fused_bsc`, but the FASTQ parse is a
//! BLOCK SCANNER: read large blocks with `File::read`, find every '\n' in one tight
//! `memchr_iter` pass, and serve records from precomputed line offsets — no per-line
//! function-call overhead, no BufReader. Output order is byte-identical to the
//! read_until path (validate on chr20). Plain (uncompressed) input only.
//!
//! Usage: bench_reorder_cluster_fused_bp <R1.fastq> <R2.fastq> <tmpdir>
//!   QZ_PASS1_ONLY=1  -> stop after pass1, print read/hash/write profile
//!   QZ_BP_BLOCK_MB=N -> block size in MiB (default 64; small values stress carry logic)

use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufWriter, Read, Write};
use std::time::Instant;

use qz_lib::compression::bsc;
use qz_lib::compression::dna_utils::compute_min_syncmer_hash;

const NBUCKETS: usize = 4096;
const SAMPLE_READS: usize = 8_000_000;
const FLUSH: usize = 256 * 1024 * 1024;

/// Block-buffered FASTQ reader. Reads ~block bytes at a time, finds all newlines once
/// per block, and serves the seq line of each 4-line record as a (start,end) range into
/// its internal buffer (valid until the next `next()` that triggers a refill — callers
/// copy immediately). Carries the partial trailing record across block boundaries.
struct BlockReader {
    f: File,
    buf: Vec<u8>,
    filled: usize,    // valid bytes buf[0..filled]
    nls: Vec<usize>,  // '\n' offsets within buf[0..filled]
    line: usize,      // next record consumes nls[line..line+4]
    rec_start: usize, // byte offset where the current (unconsumed) record region begins
    block: usize,
    eof: bool,
}

impl BlockReader {
    fn new(path: &str, block: usize) -> BlockReader {
        let f = File::open(path).expect("open input");
        let mut br = BlockReader {
            f,
            buf: vec![0u8; block + (1 << 20)],
            filled: 0,
            nls: Vec::new(),
            line: 0,
            rec_start: 0,
            block,
            eof: false,
        };
        br.refill();
        br
    }

    fn refill(&mut self) {
        // carry the incomplete trailing record [rec_start..filled] to the front
        let rem = self.filled - self.rec_start;
        if self.rec_start > 0 {
            self.buf.copy_within(self.rec_start..self.filled, 0);
        }
        self.filled = rem;
        self.rec_start = 0;
        self.line = 0;
        if self.buf.len() < self.filled + self.block {
            self.buf.resize(self.filled + self.block, 0);
        }
        // read until we have ~block fresh bytes or EOF
        let target = rem + self.block;
        while self.filled < target && !self.eof {
            match self.f.read(&mut self.buf[self.filled..]) {
                Ok(0) => self.eof = true,
                Ok(n) => self.filled += n,
                Err(e) => panic!("read error: {e}"),
            }
        }
        self.nls.clear();
        self.nls
            .extend(memchr::memchr_iter(b'\n', &self.buf[..self.filled]));
    }

    /// Seq byte range of the next record within `self.buf`, or None at EOF.
    fn next(&mut self) -> Option<(usize, usize)> {
        if self.line + 4 > self.nls.len() {
            if self.eof {
                return None;
            }
            self.refill();
            if self.line + 4 > self.nls.len() {
                return None;
            }
        }
        let seq_start = self.nls[self.line] + 1;
        let seq_end = self.nls[self.line + 1];
        self.line += 4;
        self.rec_start = self.nls[self.line - 1] + 1;
        Some((seq_start, seq_end))
    }
}

/// Pull the next record's seq into `arena`, returning its range there (mirrors the
/// read_until path's interface; copies immediately so the reader buffer can recycle).
fn next_seq(r: &mut BlockReader, arena: &mut Vec<u8>) -> Option<(usize, usize)> {
    let (s, e) = r.next()?;
    let a = arena.len();
    arena.extend_from_slice(&r.buf[s..e]);
    Some((a, arena.len()))
}

fn compute_splitters(path: &str, block: usize) -> Vec<u64> {
    let mut r = BlockReader::new(path, block);
    let mut arena = Vec::with_capacity(256);
    let mut keys: Vec<u64> = Vec::with_capacity(SAMPLE_READS.min(1 << 20));
    while keys.len() < SAMPLE_READS {
        arena.clear();
        match next_seq(&mut r, &mut arena) {
            Some((s, e)) => keys.push(compute_min_syncmer_hash(&arena[s..e])),
            None => break,
        }
    }
    keys.par_sort_unstable();
    let n = keys.len().max(1);
    (1..NBUCKETS).map(|i| keys[(i * n / NBUCKETS).min(n - 1)]).collect()
}

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
        eprintln!("usage: bench_reorder_cluster_fused_bp <R1.fastq> <R2.fastq> <tmpdir>");
        std::process::exit(2);
    }
    let (in1, in2, tmp) = (&args[1], &args[2], &args[3]);
    let block_mb: usize = std::env::var("QZ_BP_BLOCK_MB").ok().and_then(|s| s.parse().ok()).unwrap_or(64);
    let block = block_mb * 1024 * 1024;
    let _ = fs::remove_dir_all(tmp);
    fs::create_dir_all(tmp).expect("mkdir tmp");
    let t0 = Instant::now();

    let splitters = compute_splitters(in1, block);
    eprintln!("[sample] block={block_mb}MiB {:.0}s", t0.elapsed().as_secs_f64());

    let mut writers: Vec<BufWriter<File>> = (0..NBUCKETS)
        .map(|i| BufWriter::with_capacity(1 << 18, File::create(format!("{tmp}/b{i:04}")).unwrap()))
        .collect();
    let mut r1 = BlockReader::new(in1, block);
    let mut r2 = BlockReader::new(in2, block);
    const BATCH: usize = 400_000;
    let mut buf1: Vec<u8> = Vec::with_capacity(BATCH * 160);
    let mut buf2: Vec<u8> = Vec::with_capacity(BATCH * 160);
    let mut off1: Vec<(usize, usize)> = Vec::with_capacity(BATCH);
    let mut off2: Vec<(usize, usize)> = Vec::with_capacity(BATCH);
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
            match next_seq(&mut r1, &mut buf1) {
                Some(r) => off1.push(r),
                None => break,
            }
            match next_seq(&mut r2, &mut buf2) {
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
    eprintln!("[pass1] {npairs} pairs (BLOCK-PARSE) {pass1:.0}s | read(parse) {read_s:.0}s  hash {hash_s:.0}s  write(route+bucket) {write_s:.0}s");
    if std::env::var("QZ_PASS1_ONLY").is_ok() {
        let _ = fs::remove_dir_all(tmp);
        println!("PASS1_ONLY npairs={npairs} pass1_s={pass1:.0} read_s={read_s:.0} hash_s={hash_s:.0} write_s={write_s:.0}");
        return;
    }

    // Pass 2: per bucket sort by (key,idx), inline BSC the clustered seq (same as fused).
    let t2 = Instant::now();
    let mut b1: Vec<u8> = Vec::with_capacity(FLUSH + (1 << 20));
    let mut b2: Vec<u8> = Vec::with_capacity(FLUSH + (1 << 20));
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
            b1.extend_from_slice(&data[*o1..*o1 + *l1]);
            b2.extend_from_slice(&data[*o2..*o2 + *l2]);
            if b1.len() >= FLUSH {
                flush_bsc(&mut b1, &mut bsc1, &mut bsc1s);
            }
            if b2.len() >= FLUSH {
                flush_bsc(&mut b2, &mut bsc2, &mut bsc2s);
            }
        }
        written += ents.len() as u64;
        let _ = fs::remove_file(&fp);
    }
    flush_bsc(&mut b1, &mut bsc1, &mut bsc1s);
    flush_bsc(&mut b2, &mut bsc2, &mut bsc2s);
    let _ = fs::remove_dir_all(tmp);
    assert_eq!(written, npairs);
    let pass2 = t2.elapsed().as_secs_f64();
    let total = t0.elapsed().as_secs_f64();
    eprintln!("[pass2] sort {sort_s:.0}s, BSC inline R1 {bsc1s:.0}s R2 {bsc2s:.0}s, {pass2:.0}s");
    eprintln!("[done] max_bucket={max_bucket}  total {total:.0}s");
    println!(
        "FUSED_BP npairs={npairs} bsc_seq_R1={bsc1} bsc_seq_R2={bsc2} bsc_seq_total={} pass1_s={pass1:.0} pass2_s={pass2:.0} total_s={total:.0} sort_s={sort_s:.0} bsc_s={:.0}",
        bsc1 + bsc2,
        bsc1s + bsc2s
    );
}
