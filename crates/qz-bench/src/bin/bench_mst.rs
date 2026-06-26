//! BOUNDED Mstcom-style donor-diff prototype (§8r): does encoding each read as a diff
//! vs its single best REAL-READ donor beat cluster+zstd-long? The donor is found online
//! (earlier reads only) via a minimizer bucket; the SHIFT comes from the shared
//! minimizer's position (Mstcom's Hamming-shifting, cheap minimizers-only "speedy").
//! NO consensus / NO assembly (FTO-clear): the diff target is always an actual input read.
//!
//! The donor reference is an EXPLICIT back-pointer (i - donor_idx): the minimizer bucket
//! is an encode-only device to find a good recent donor; the decoder just follows the
//! stored pointer (a bucket rank is NOT decodable — the decoder can't rebuild read i's
//! bucket before read i exists). In clustered input order the donor is usually recent, so
//! the back-pointer stays small + BSC-friendly. Residuals → BSC (NOT a new coder, §8p dead).
//!
//! Read i = VERBATIM (no good donor) or DIFF(back, shift, exception-list); an "exception"
//! is any read position whose base != the shift-aligned donor base, or that is uncovered
//! (outside the donor). Streams → BSC; bpb = sum/total_bases. Replay decoder asserts lossless.
//!
//! A/B target: cluster+zstd-16/19 --long on the SAME reads (0.600/0.575 bpb). KILL: if this
//! doesn't clearly beat cluster+zstd-long, the reference-free ratio is at its IP-safe ceiling.
//!
//! Env: QZ_MST_K (13), QZ_MST_DEPTH (32), QZ_MST_TBITS (24), QZ_MST_MAXEXC (def read_len/3),
//!      QZ_MST_MAXREADS (subset). Usage: bench_mst <clustered_r1.fastq>

use qz_lib::compression::bsc;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

fn envu(name: &str, def: usize) -> usize {
    std::env::var(name).ok().and_then(|s| s.parse().ok()).unwrap_or(def)
}

fn minimizer(read: &[u8], kk: usize) -> (u64, u16) {
    let n = read.len();
    if n < kk {
        let mut h = 0u64;
        for &b in read {
            h = h.wrapping_mul(131).wrapping_add(b as u64);
        }
        return (h, 0);
    }
    let mut pow = 1u64;
    for _ in 0..kk - 1 {
        pow = pow.wrapping_mul(131);
    }
    let mut h = 0u64;
    for &b in &read[..kk] {
        h = h.wrapping_mul(131).wrapping_add(b as u64);
    }
    let (mut best, mut bestpos) = (h, 0usize);
    for j in 1..=n - kk {
        h = h
            .wrapping_sub((read[j - 1] as u64).wrapping_mul(pow))
            .wrapping_mul(131)
            .wrapping_add(read[j + kk - 1] as u64);
        if h < best {
            best = h;
            bestpos = j;
        }
    }
    (best, bestpos as u16)
}

/// Fixed bucket table: each bucket holds the last DEPTH (read_idx, min_pos). Encode-only.
struct Buckets {
    idx: Vec<u32>,
    mpos: Vec<u16>,
    cnt: Vec<u32>,
    depth: usize,
    mask: usize,
}
impl Buckets {
    fn new(tbits: usize, depth: usize) -> Self {
        let nb = 1usize << tbits;
        Buckets { idx: vec![0u32; nb * depth], mpos: vec![0u16; nb * depth], cnt: vec![0u32; nb], depth, mask: nb - 1 }
    }
    #[inline]
    fn bucket(&self, h: u64) -> usize {
        (h as usize) & self.mask
    }
    #[inline]
    fn push(&mut self, b: usize, ri: u32, mp: u16) {
        let c = self.cnt[b] as usize;
        let slot = b * self.depth + (c % self.depth);
        self.idx[slot] = ri;
        self.mpos[slot] = mp;
        self.cnt[b] += 1;
    }
    #[inline]
    fn live(&self, b: usize) -> usize {
        (self.cnt[b] as usize).min(self.depth)
    }
    #[inline]
    fn at(&self, b: usize, rank: usize) -> (u32, u16) {
        let c = self.cnt[b] as usize;
        let slot = b * self.depth + ((c - 1 - rank) % self.depth);
        (self.idx[slot], self.mpos[slot])
    }
}

#[inline]
fn count_exc(read: &[u8], donor: &[u8], shift: i32, cap: usize) -> usize {
    let dl = donor.len() as i32;
    let mut e = 0usize;
    for (p, &rb) in read.iter().enumerate() {
        let dp = p as i32 + shift;
        if dp < 0 || dp >= dl || donor[dp as usize] != rb {
            e += 1;
            if e > cap {
                return e;
            }
        }
    }
    e
}

fn put_varint(v: &mut Vec<u8>, mut x: usize) {
    loop {
        let b = (x & 0x7f) as u8;
        x >>= 7;
        if x == 0 {
            v.push(b);
            break;
        }
        v.push(b | 0x80);
    }
}
fn get_varint(v: &[u8], i: &mut usize) -> usize {
    let (mut x, mut s) = (0usize, 0u32);
    loop {
        let b = v[*i];
        *i += 1;
        x |= ((b & 0x7f) as usize) << s;
        if b & 0x80 == 0 {
            break;
        }
        s += 7;
    }
    x
}
fn zigzag(s: i32) -> usize {
    ((s << 1) ^ (s >> 31)) as u32 as usize
}
fn unzig(z: usize) -> i32 {
    let z = z as u32;
    ((z >> 1) as i32) ^ -((z & 1) as i32)
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: bench_mst <clustered_r1.fastq>");
    let kk = envu("QZ_MST_K", 13);
    let depth = envu("QZ_MST_DEPTH", 32);
    let tbits = envu("QZ_MST_TBITS", 24);
    let maxexc_env = std::env::var("QZ_MST_MAXEXC").ok().and_then(|s| s.parse::<usize>().ok());
    let maxreads = envu("QZ_MST_MAXREADS", usize::MAX);

    // ---- load reads ----
    let t0 = Instant::now();
    let mut r = BufReader::with_capacity(1 << 22, File::open(&path).expect("open"));
    let mut buf: Vec<u8> = Vec::new();
    let mut offs: Vec<u32> = vec![0];
    let mut line = Vec::with_capacity(256);
    let mut ln = 0u64;
    while r.read_until(b'\n', &mut line).unwrap() != 0 {
        if ln % 4 == 1 {
            let s = if line.last() == Some(&b'\n') { &line[..line.len() - 1] } else { &line[..] };
            buf.extend_from_slice(s);
            offs.push(buf.len() as u32);
            if offs.len() > maxreads {
                break;
            }
        }
        ln += 1;
        line.clear();
    }
    let nreads = offs.len() - 1;
    let total_bases = buf.len();
    let read = |i: usize| -> &[u8] { &buf[offs[i] as usize..offs[i + 1] as usize] };
    eprintln!(
        "[load] {nreads} reads, {:.2} Gbase, {:.0}s  (K={kk} DEPTH={depth} TBITS={tbits})",
        total_bases as f64 / 1e9,
        t0.elapsed().as_secs_f64()
    );

    let tm = Instant::now();
    let mins: Vec<(u64, u16)> = (0..nreads).into_par_iter().map(|i| minimizer(read(i), kk)).collect();
    eprintln!("[minimizers] {:.0}s", tm.elapsed().as_secs_f64());

    // ---- ENCODE ----
    let te = Instant::now();
    let mut flags: Vec<u8> = Vec::new();
    let (mut facc, mut fbits) = (0u8, 0u8);
    let mut lens: Vec<u8> = Vec::new();
    let mut donor_back: Vec<u8> = Vec::new();
    let mut shifts: Vec<u8> = Vec::new();
    let mut nexcs: Vec<u8> = Vec::new();
    let mut exc_pos: Vec<u8> = Vec::new();
    let mut exc_base: Vec<u8> = Vec::new();
    let mut verb: Vec<u8> = Vec::new();

    let mut bk = Buckets::new(tbits, depth);
    let mut n_diff = 0usize;
    let mut exc_total = 0u64;
    let mut back_sum = 0u64;
    for i in 0..nreads {
        let ri = read(i);
        let li = ri.len();
        put_varint(&mut lens, li);
        let (mh, mp) = mins[i];
        let b = bk.bucket(mh);
        let maxexc = maxexc_env.unwrap_or(li / 3);
        let mut best = usize::MAX;
        let mut best_donor = u32::MAX;
        let mut best_shift = 0i32;
        let live = bk.live(b);
        for rank in 0..live {
            let (dj, dmp) = bk.at(b, rank);
            let shift = dmp as i32 - mp as i32;
            let cap = best.min(maxexc);
            let e = count_exc(ri, read(dj as usize), shift, cap);
            if e < best {
                best = e;
                best_donor = dj;
                best_shift = shift;
                if best == 0 {
                    break;
                }
            }
        }
        let bit = if best <= maxexc && best_donor != u32::MAX { 1u8 } else { 0u8 };
        facc |= bit << fbits;
        fbits += 1;
        if fbits == 8 {
            flags.push(facc);
            facc = 0;
            fbits = 0;
        }
        if bit == 1 {
            let back = i as u32 - best_donor; // ≥1, donor earlier
            put_varint(&mut donor_back, back as usize);
            back_sum += back as u64;
            put_varint(&mut shifts, zigzag(best_shift));
            let donor = read(best_donor as usize);
            let dl = donor.len() as i32;
            let mut prev = 0usize;
            let mut ne = 0usize;
            for (p, &rb) in ri.iter().enumerate() {
                let dp = p as i32 + best_shift;
                if dp < 0 || dp >= dl || donor[dp as usize] != rb {
                    put_varint(&mut exc_pos, p - prev);
                    prev = p;
                    exc_base.push(rb);
                    ne += 1;
                }
            }
            put_varint(&mut nexcs, ne);
            exc_total += ne as u64;
            n_diff += 1;
        } else {
            verb.extend_from_slice(ri);
        }
        bk.push(b, i as u32, mp);
    }
    if fbits > 0 {
        flags.push(facc);
    }
    let enc_s = te.elapsed().as_secs_f64();

    let cz = |v: &[u8]| if v.is_empty() { 0 } else { bsc::compress_parallel(v).expect("bsc").len() };
    let parts = [
        ("flags", cz(&flags)),
        ("lens", cz(&lens)),
        ("donor", cz(&donor_back)),
        ("shift", cz(&shifts)),
        ("nexc", cz(&nexcs)),
        ("excpos", cz(&exc_pos)),
        ("excbase", cz(&exc_base)),
        ("verb", cz(&verb)),
    ];
    let total: usize = parts.iter().map(|p| p.1).sum();
    let bpb = total as f64 * 8.0 / total_bases as f64;

    // ---- DECODE (replay, prove lossless) ----
    let td = Instant::now();
    let mut dbuf: Vec<u8> = Vec::with_capacity(total_bases);
    let mut doffs: Vec<u32> = Vec::with_capacity(nreads + 1);
    doffs.push(0);
    let (mut fi, mut li_p, mut db_p, mut sh_p, mut ne_p, mut ep_p, mut eb_p, mut vb_p) =
        (0usize, 0usize, 0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
    for i in 0..nreads {
        let bit = (flags[fi / 8] >> (fi % 8)) & 1;
        fi += 1;
        let li = get_varint(&lens, &mut li_p);
        let start = dbuf.len();
        dbuf.resize(start + li, 0);
        if bit == 1 {
            let back = get_varint(&donor_back, &mut db_p);
            let donor_idx = i - back;
            let shift = unzig(get_varint(&shifts, &mut sh_p));
            let ne = get_varint(&nexcs, &mut ne_p);
            let ds = doffs[donor_idx] as usize;
            let de = doffs[donor_idx + 1] as usize;
            let dl = (de - ds) as i32;
            // copy shift-aligned donor bases for covered positions
            for p in 0..li {
                let dp = p as i32 + shift;
                if dp >= 0 && dp < dl {
                    dbuf[start + p] = dbuf[ds + dp as usize];
                }
            }
            // overwrite exceptions (covered mismatches + uncovered)
            let mut prev = 0usize;
            for k in 0..ne {
                let d = get_varint(&exc_pos, &mut ep_p);
                let p = if k == 0 { d } else { prev + d };
                prev = p;
                dbuf[start + p] = exc_base[eb_p];
                eb_p += 1;
            }
        } else {
            dbuf[start..start + li].copy_from_slice(&verb[vb_p..vb_p + li]);
            vb_p += li;
        }
        doffs.push(dbuf.len() as u32);
    }
    let dec_s = td.elapsed().as_secs_f64();
    let lossless = dbuf == buf && doffs == offs;

    // ---- BSC-raw anchor ----
    let raw_bsc = bsc::compress_parallel(&buf).expect("bsc").len();

    let diff_pct = 100.0 * n_diff as f64 / nreads as f64;
    let avg_exc = if n_diff > 0 { exc_total as f64 / n_diff as f64 } else { 0.0 };
    let avg_back = if n_diff > 0 { back_sum as f64 / n_diff as f64 } else { 0.0 };
    eprintln!(
        "[parse] diff={n_diff} ({diff_pct:.1}%)  avg_exc/diff={avg_exc:.1}  avg_back={avg_back:.0}"
    );
    let mut sl = String::new();
    for (n, s) in &parts {
        sl.push_str(&format!("{n}={s} "));
    }
    eprintln!("[streams BSC] {sl}");
    eprintln!("[time] enc {enc_s:.0}s dec {dec_s:.0}s  roundtrip={}", if lossless { "OK" } else { "FAIL" });
    assert!(lossless, "MST roundtrip MISMATCH — not lossless");
    println!(
        "MST reads={nreads} bpb={bpb:.3}  total={total}  vs_BSC={:+.1}%  BSC_raw={raw_bsc}({:.3}bpb)  diff={diff_pct:.0}% K={kk} DEPTH={depth}",
        100.0 * (total as f64 - raw_bsc as f64) / raw_bsc as f64,
        raw_bsc as f64 * 8.0 / total_bases as f64,
    );
}
