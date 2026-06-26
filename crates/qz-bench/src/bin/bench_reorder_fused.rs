//! FUSED order-drop seq pipeline (measurement of plan (1)): one process that
//!   parse -> (syncmer key + RC-canon orientation + minimizer-anchor pos) per read
//!   -> external range-partition bucket sort (bounded RAM)
//!   -> within-bucket sort by (key, anchor-suffix)
//!   -> stream the canonicalized seq (1 byte/base, no separators) to STDOUT
//! Pipe stdout into `zstd -16 --long=31 -T16` for the compressed sequence stream.
//! This removes the disk round-trips of the 3-step (reorder-fastq | canon-resort | zstd)
//! chain: no full-FASTQ materialize, no seq re-extraction, no canon blob on disk.
//! Strand bits (1/read, in emit order) go to <bits_out>. Parallel hash (A1). Bounded RAM.
//!
//! Output seq order == bench_canon_resort QZ_SORT=anchorseq, so the piped zstd size must
//! match the ~105 MB anchorseq result (fusion preserves ratio).
//!
//! Usage: bench_reorder_fused <in.fastq[.gz]> <tmpdir> <bits_out>   (seq stream -> stdout)

use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::time::Instant;

use qz_lib::compression::dna_utils::{
    canonical_minimizer_orientation, compute_min_syncmer_hash, compute_min_syncmer_hash_pos,
    reverse_complement_canonical,
};

const NBUCKETS: usize = 4096;
const SAMPLE_READS: usize = 8_000_000;
const BATCH: usize = 1_000_000;

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

/// Read one 4-line record; return seq-line byte range within `rec` or None at EOF.
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

/// Balanced quantile splitters on the canonical syncmer key (RC-invariant, so computed
/// on raw seqs). Parallel hash.
fn compute_splitters(path: &str) -> Vec<u64> {
    let mut reader = open_maybe_gz(path);
    let mut rec = Vec::with_capacity(512);
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(SAMPLE_READS.min(1 << 20));
    while seqs.len() < SAMPLE_READS {
        match read_record(&mut reader, &mut rec) {
            Some((ss, se)) => seqs.push(rec[ss..se].to_vec()),
            None => break,
        }
    }
    let mut keys: Vec<u64> = seqs.par_iter().map(|s| compute_min_syncmer_hash(s)).collect();
    keys.par_sort_unstable();
    let n = keys.len().max(1);
    (1..NBUCKETS).map(|i| keys[(i * n / NBUCKETS).min(n - 1)]).collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: bench_reorder_fused <in.fastq[.gz]> <tmpdir> <bits_out>   (seq -> stdout)");
        std::process::exit(2);
    }
    let (inp, tmp, bits_out) = (&args[1], &args[2], &args[3]);
    // QZ_PLAIN=1 -> plain order-drop baseline: SAME RC-invariant clustering, but NO RC-canon
    // (flip always false) and NO anchor-suffix sort (within-bucket by key only). Isolates the
    // ratio contribution of the RC-canon + anchorseq stack at deep WGS.
    let plain = std::env::var("QZ_PLAIN").map(|v| v == "1").unwrap_or(false);
    eprintln!("[mode] {}", if plain { "PLAIN order-drop (no canon, no anchor)" } else { "FULL stack (RC-canon + anchorseq)" });
    let _ = fs::remove_dir_all(tmp);
    fs::create_dir_all(tmp).expect("mkdir tmp");
    let t0 = Instant::now();

    let splitters = compute_splitters(inp);
    eprintln!("[sample] {} splitters ({:.0}s)", splitters.len(), t0.elapsed().as_secs_f64());

    // ---- Pass 1: parse + canon + key + anchorpos -> buckets ----
    // bucket entry layout: [key u64][apos u16][flip u8][len u16][canon seq len bytes]
    let mut writers: Vec<BufWriter<File>> = (0..NBUCKETS)
        .map(|i| BufWriter::with_capacity(1 << 18, File::create(format!("{tmp}/b{i:04}")).unwrap()))
        .collect();
    let mut reader = open_maybe_gz(inp);
    let mut rec = Vec::with_capacity(512);
    let mut nreads: u64 = 0;
    let mut arena: Vec<u8> = Vec::with_capacity(BATCH * 320);
    let mut meta: Vec<(usize, usize, usize)> = Vec::with_capacity(BATCH); // (start, ss, se)
    loop {
        arena.clear();
        meta.clear();
        while meta.len() < BATCH {
            match read_record(&mut reader, &mut rec) {
                Some((ss, se)) => {
                    let start = arena.len();
                    arena.extend_from_slice(&rec[ss..se]); // store ONLY the seq
                    meta.push((start, se - ss, 0));
                }
                None => break,
            }
        }
        if meta.is_empty() {
            break;
        }
        // parallel: canon orientation + key + anchor pos on the canonical seq
        let computed: Vec<(u64, u16, bool, Vec<u8>)> = meta
            .par_iter()
            .map(|&(st, len, _)| {
                let seq = &arena[st..st + len];
                if plain {
                    // raw orientation, RC-invariant key still buckets identically to FULL
                    let (key, apos) = compute_min_syncmer_hash_pos(seq);
                    (key, apos, false, seq.to_vec())
                } else {
                    let flip = canonical_minimizer_orientation(seq);
                    let canon = if flip { reverse_complement_canonical(seq) } else { seq.to_vec() };
                    let (key, apos) = compute_min_syncmer_hash_pos(&canon);
                    (key, apos, flip, canon)
                }
            })
            .collect();
        for (key, apos, flip, canon) in &computed {
            let bucket = splitters.partition_point(|&s| s <= *key);
            let w = &mut writers[bucket];
            w.write_all(&key.to_le_bytes()).unwrap();
            w.write_all(&apos.to_le_bytes()).unwrap();
            w.write_all(&[*flip as u8]).unwrap();
            w.write_all(&(canon.len() as u16).to_le_bytes()).unwrap();
            w.write_all(canon).unwrap();
            nreads += 1;
        }
    }
    for w in &mut writers {
        w.flush().unwrap();
    }
    drop(writers);
    eprintln!("[pass1] {} reads -> {} buckets ({:.0}s)", nreads, NBUCKETS, t0.elapsed().as_secs_f64());

    // ---- Pass 2: within-bucket sort by (key, anchor-suffix); stream canon seq to stdout ----
    let t2 = Instant::now();
    let stdout = std::io::stdout();
    let mut out = BufWriter::with_capacity(1 << 22, stdout.lock());
    let mut bits = vec![0u8; (nreads as usize).div_ceil(8)];
    let mut emit: usize = 0;
    let mut written: u64 = 0;
    for i in 0..NBUCKETS {
        let data = fs::read(format!("{tmp}/b{i:04}")).unwrap();
        // parse entries -> (key, apos, flip, seq_off, seq_len)
        let mut ents: Vec<(u64, u16, bool, usize, usize)> = Vec::new();
        let mut p = 0usize;
        while p + 13 <= data.len() {
            let key = u64::from_le_bytes(data[p..p + 8].try_into().unwrap());
            let apos = u16::from_le_bytes(data[p + 8..p + 10].try_into().unwrap());
            let flip = data[p + 10] != 0;
            let len = u16::from_le_bytes(data[p + 11..p + 13].try_into().unwrap()) as usize;
            let off = p + 13;
            ents.push((key, apos, flip, off, len));
            p = off + len;
        }
        // sort by key, then by the anchor-suffix canon[apos..] (FULL) or key only (PLAIN)
        if plain {
            ents.par_sort_unstable_by_key(|a| a.0);
        } else {
            ents.par_sort_unstable_by(|a, b| {
                a.0.cmp(&b.0).then_with(|| {
                    let sa = &data[a.3 + (a.1 as usize).min(a.4)..a.3 + a.4];
                    let sb = &data[b.3 + (b.1 as usize).min(b.4)..b.3 + b.4];
                    sa.cmp(sb)
                })
            });
        }
        for (_k, _ap, flip, off, len) in &ents {
            out.write_all(&data[*off..*off + *len]).unwrap();
            if *flip {
                bits[emit >> 3] |= 1 << (emit & 7);
            }
            emit += 1;
            written += *len as u64;
        }
    }
    out.flush().unwrap();
    fs::write(bits_out, &bits).expect("write strand bits");
    let _ = fs::remove_dir_all(tmp);
    assert_eq!(emit as u64, nreads);
    eprintln!(
        "[pass2] sorted+streamed {} reads, {} seq bytes, strand-bits {} B ({:.0}s)",
        nreads, written, bits.len(), t2.elapsed().as_secs_f64()
    );
    eprintln!("[done] pipeline (excl. zstd) {:.0}s", t0.elapsed().as_secs_f64());
}
