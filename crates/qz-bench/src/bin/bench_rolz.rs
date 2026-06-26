//! ROLZ-bio prototype: context-anchored REDUCED-OFFSET LZ on the clustered seq blob.
//! Tests whether hash-anchored reach (vs xz's serial big window) can capture the
//! −27% ratio headroom (§8n) at parallel-block speed.
//!
//! Core ROLZ (decodable): context = preceding K bytes → fixed 2^TBITS hash table →
//! DEPTH-ring of recent positions; emitted offset is the RANK in the bucket (few bits),
//! NOT a byte distance, so a long-range match costs the same as a near one.
//!
//! RESIDUAL coding (§8o, design wf_c3d05fd3-ff6 "ROLZ-LZMA-lite"):
//!   QZ_ROLZ_LITCODER=bsc   (default) — literals → one BSC blob (the original baseline)
//!   QZ_ROLZ_LITCODER=range          — literals → an LZMA-style binary range coder with a
//!     DNA literal model (is_escape / is_N / 2-bit ACGT tree), order-2 preceding-base
//!     context. Per-block model reset → block-parallel. Decoder replays the identical
//!     model from already-decoded output (proven by the existing roundtrip assert).
//!   QZ_ROLZ_DIAG=1 — print rep0-potential + literal-entropy diagnostics (decides Stage 2).
//!
//! Env: QZ_ROLZ_K, QZ_ROLZ_DEPTH, QZ_ROLZ_MINMATCH, QZ_ROLZ_MAXMATCH, QZ_ROLZ_TBITS.
//! Usage: bench_rolz <blob> [bytes_limit]

use qz_lib::compression::bsc;
use std::time::Instant;

fn envu(name: &str, def: usize) -> usize {
    std::env::var(name).ok().and_then(|s| s.parse().ok()).unwrap_or(def)
}

const B: u64 = 131; // rolling-hash base

// ACGT → 0..3, everything else (incl. N) → 0xFF. For literal BASE coding.
fn build_base_lut() -> [u8; 256] {
    let mut t = [0xFFu8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t
}
// ACGT → 0..3, everything else → 0. For CONTEXT (sentinel 0).
fn build_ctx_lut() -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t
}
const INV_BASE: [u8; 4] = [b'A', b'C', b'G', b'T'];

// ---------------- LZMA-style binary range coder (carry/cache variant) ----------------
const PROB_INIT: u16 = 1024; // 0.5 in 11-bit
const K_TOP: u32 = 1 << 24;

struct RangeEncoder {
    low: u64,
    range: u32,
    cache: u8,
    cache_size: u64,
    out: Vec<u8>,
}
impl RangeEncoder {
    fn new() -> Self {
        RangeEncoder { low: 0, range: 0xFFFF_FFFF, cache: 0, cache_size: 1, out: Vec::new() }
    }
    #[inline]
    fn shift_low(&mut self) {
        if (self.low >> 32) != 0 || self.low < 0xFF00_0000 {
            let carry = (self.low >> 32) as u8;
            let mut temp = self.cache;
            loop {
                self.out.push(temp.wrapping_add(carry));
                temp = 0xFF;
                self.cache_size -= 1;
                if self.cache_size == 0 {
                    break;
                }
            }
            self.cache = (self.low >> 24) as u8;
        }
        self.cache_size += 1;
        self.low = (self.low << 8) & 0xFFFF_FFFF;
    }
    #[inline]
    fn enc_bit(&mut self, prob: &mut u16, bit: u32) {
        let bound = (self.range >> 11) * (*prob as u32);
        if bit == 0 {
            self.range = bound;
            *prob += (2048 - *prob) >> 5;
        } else {
            self.low += bound as u64;
            self.range -= bound;
            *prob -= *prob >> 5;
        }
        while self.range < K_TOP {
            self.range <<= 8;
            self.shift_low();
        }
    }
    #[inline]
    fn encode_raw(&mut self, raw: &mut [u16; 256], byte: u8) {
        let mut m = 1usize;
        for i in (0..8).rev() {
            let bit = ((byte >> i) & 1) as u32;
            self.enc_bit(&mut raw[m], bit);
            m = (m << 1) | (bit as usize);
        }
    }
    fn finish(&mut self) -> &[u8] {
        for _ in 0..5 {
            self.shift_low();
        }
        &self.out
    }
}

struct RangeDecoder<'a> {
    code: u32,
    range: u32,
    buf: &'a [u8],
    pos: usize,
}
impl<'a> RangeDecoder<'a> {
    fn new(buf: &'a [u8]) -> Self {
        let mut code = 0u32;
        // buf[0] is the ignored leading byte; next 4 prime `code`.
        for i in 1..5 {
            code = (code << 8) | (*buf.get(i).unwrap_or(&0) as u32);
        }
        RangeDecoder { code, range: 0xFFFF_FFFF, buf, pos: 5 }
    }
    #[inline]
    fn dec_bit(&mut self, prob: &mut u16) -> u32 {
        let bound = (self.range >> 11) * (*prob as u32);
        let bit;
        if self.code < bound {
            self.range = bound;
            *prob += (2048 - *prob) >> 5;
            bit = 0;
        } else {
            self.code -= bound;
            self.range -= bound;
            *prob -= *prob >> 5;
            bit = 1;
        }
        while self.range < K_TOP {
            self.range <<= 8;
            self.code = (self.code << 8) | (*self.buf.get(self.pos).unwrap_or(&0) as u32);
            self.pos += 1;
        }
        bit
    }
    #[inline]
    fn decode_raw(&mut self, raw: &mut [u16; 256]) -> u8 {
        let mut m = 1usize;
        for _ in 0..8 {
            let bit = self.dec_bit(&mut raw[m]);
            m = (m << 1) | (bit as usize);
        }
        (m & 0xFF) as u8
    }
}

// ---------------- DNA literal model (Stage 1: FRESH path only) ----------------
// order-k preceding-base context, NCTX = 4^k. (QZ_ROLZ_ORDERK, default 2.)
struct LitModel {
    escape: [u16; 4],          // is_escape, ctx = prev1 base
    isn: Vec<u16>,             // is_N  [NCTX]
    t0: Vec<u16>,              // base bit0 [NCTX]
    t1: Vec<[u16; 2]>,         // base bit1 | bit0 [NCTX]
    raw: [u16; 256],           // order-0 byte tree for escaped bytes
    mask: usize,               // NCTX-1
}
impl LitModel {
    fn reset(orderk: usize) -> Self {
        let nctx = 1usize << (2 * orderk);
        LitModel {
            escape: [PROB_INIT; 4],
            isn: vec![PROB_INIT; nctx],
            t0: vec![PROB_INIT; nctx],
            t1: vec![[PROB_INIT; 2]; nctx],
            raw: [PROB_INIT; 256],
            mask: nctx - 1,
        }
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: bench_rolz <blob> [bytes_limit]");
    let limit = std::env::args().nth(2).and_then(|s| s.parse::<usize>().ok());
    let k = envu("QZ_ROLZ_K", 16);
    let depth = envu("QZ_ROLZ_DEPTH", 8);
    let min_match = envu("QZ_ROLZ_MINMATCH", 32);
    let max_match = envu("QZ_ROLZ_MAXMATCH", 8192);
    let tbits = envu("QZ_ROLZ_TBITS", 22);
    let orderk = envu("QZ_ROLZ_ORDERK", 2).clamp(1, 11);
    let use_range = std::env::var("QZ_ROLZ_LITCODER").map(|s| s == "range").unwrap_or(false);
    let diag = std::env::var("QZ_ROLZ_DIAG").map(|s| s == "1").unwrap_or(false);
    assert!(depth <= 256, "rank must fit a byte");

    let base_lut = build_base_lut();
    let ctx_lut = build_ctx_lut();

    let mut data = std::fs::read(&path).expect("read blob");
    if let Some(l) = limit {
        data.truncate(l);
    }
    let n = data.len();
    eprintln!(
        "[in] {} bytes ({:.2} GB)  K={k} DEPTH={depth} MIN={min_match} MAX={max_match} TBITS={tbits} litcoder={}",
        n,
        n as f64 / 1e9,
        if use_range { "range" } else { "bsc" }
    );

    let roller = Roller::new(k);

    // ---- ENCODE ----
    let te = Instant::now();
    let mut flags: Vec<u8> = Vec::new(); // 1 bit/token, packed
    let mut flag_acc = 0u8;
    let mut flag_bits = 0u8;
    let mut ranks: Vec<u8> = Vec::new();
    let mut lens: Vec<u8> = Vec::new();
    let mut lits: Vec<u8> = Vec::new(); // BSC path
    let mut renc = RangeEncoder::new(); // range path
    let mut model = LitModel::reset(orderk);

    let emit_flag = |flags: &mut Vec<u8>, flag_acc: &mut u8, flag_bits: &mut u8, bit: u8| {
        *flag_acc |= bit << *flag_bits;
        *flag_bits += 1;
        if *flag_bits == 8 {
            flags.push(*flag_acc);
            *flag_acc = 0;
            *flag_bits = 0;
        }
    };
    let put_varint = |v: &mut Vec<u8>, mut x: usize| loop {
        let b = (x & 0x7f) as u8;
        x >>= 7;
        if x == 0 {
            v.push(b);
            break;
        } else {
            v.push(b | 0x80);
        }
    };

    // order-k literal context packed from the preceding bases; b1 = prev1 base (escape ctx)
    let lit_ctx = |buf: &[u8], p: usize, mask: usize| -> (usize, usize) {
        let mut ctx = 0usize;
        for j in 1..=orderk {
            let base = if p >= j { ctx_lut[buf[p - j] as usize] as usize } else { 0 };
            ctx |= base << (2 * (j - 1));
        }
        let b1 = if p >= 1 { ctx_lut[buf[p - 1] as usize] as usize } else { 0 };
        (ctx & mask, b1)
    };
    let enc_lit = |renc: &mut RangeEncoder, m: &mut LitModel, ctx: usize, b1: usize, byte: u8| {
        let esc = base_lut[byte as usize] == 0xFF && byte != b'N';
        renc.enc_bit(&mut m.escape[b1], esc as u32);
        if esc {
            renc.encode_raw(&mut m.raw, byte);
            return;
        }
        let isn = byte == b'N';
        renc.enc_bit(&mut m.isn[ctx], isn as u32);
        if isn {
            return;
        }
        let c = base_lut[byte as usize];
        let bit0 = (c >> 1) as u32;
        let bit1 = (c & 1) as u32;
        renc.enc_bit(&mut m.t0[ctx], bit0);
        renc.enc_bit(&mut m.t1[ctx][bit0 as usize], bit1);
    };

    let mut chains = Chains::new(tbits, depth);
    let mut ins = k;
    let mut rh = if n > k { roller.at(&data, k) } else { 0 };
    let mut matched_bytes = 0usize;
    let mut n_match = 0usize;
    let mut n_lit = 0usize;

    // Stage-0 diagnostics: rep0 potential + literal entropy
    let mut last_src_end: usize = usize::MAX; // source end of the most recent match
    let mut just_matched = false;
    let mut matched_eligible = 0usize; // literals immediately after a match
    let mut rep0_hist = [0usize; 6]; // {0,1-3,4-7,8-15,16-31,32+}
    let mut lit_o0 = [0u64; 256]; // order-0 literal histogram

    let mut p = 0usize;
    while p < n {
        while ins < p {
            let idx = chains.idx(rh);
            chains.insert(idx, ins as u32);
            let out_b = data[ins - k];
            let in_b = data[ins];
            rh = roller.roll(rh, out_b, in_b);
            ins += 1;
        }

        let mut best_len = 0usize;
        let mut best_rank = 0usize;
        let mut best_q = 0usize;
        if p >= k {
            debug_assert_eq!(ins, p);
            let idx = chains.idx(rh);
            for rank in 0..depth {
                match chains.candidate(idx, rank) {
                    None => break,
                    Some(q) => {
                        let q = q as usize;
                        if q >= p {
                            continue;
                        }
                        let l = lcp(&data, p, q, max_match);
                        if l > best_len {
                            best_len = l;
                            best_rank = rank;
                            best_q = q;
                        }
                    }
                }
            }
        }

        if best_len >= min_match {
            emit_flag(&mut flags, &mut flag_acc, &mut flag_bits, 1);
            ranks.push(best_rank as u8);
            put_varint(&mut lens, best_len - min_match);
            matched_bytes += best_len;
            n_match += 1;
            if diag {
                last_src_end = best_q + best_len;
                just_matched = true;
            }
            p += best_len;
        } else {
            emit_flag(&mut flags, &mut flag_acc, &mut flag_bits, 0);
            if use_range {
                let (ctx, b1) = lit_ctx(&data, p, model.mask);
                enc_lit(&mut renc, &mut model, ctx, b1, data[p]);
            } else {
                lits.push(data[p]);
            }
            if diag {
                lit_o0[data[p] as usize] += 1;
                if just_matched && last_src_end != usize::MAX && last_src_end < n {
                    matched_eligible += 1;
                    // rep0 run after skipping the 1 mismatch byte
                    let rl = lcp(&data, p + 1, last_src_end + 1, max_match);
                    let bucket = match rl {
                        0 => 0,
                        1..=3 => 1,
                        4..=7 => 2,
                        8..=15 => 3,
                        16..=31 => 4,
                        _ => 5,
                    };
                    rep0_hist[bucket] += 1;
                }
                just_matched = false;
            }
            n_lit += 1;
            p += 1;
        }
    }
    if flag_bits > 0 {
        flags.push(flag_acc);
    }
    let enc_s = te.elapsed().as_secs_f64();

    // ---- stream sizes ----
    let cz = |v: &[u8]| -> usize {
        if v.is_empty() {
            0
        } else {
            bsc::compress_parallel(v).expect("bsc").len()
        }
    };
    let (cf, cr, cl) = (cz(&flags), cz(&ranks), cz(&lens));
    let cli = if use_range { renc.finish().len() } else { cz(&lits) };
    let total = cf + cr + cl + cli;
    let bpb = total as f64 * 8.0 / n as f64;

    // ---- DECODE (prove lossless) ----
    let range_blob: Vec<u8> = if use_range { renc.out.clone() } else { Vec::new() };
    let td = Instant::now();
    let mut out: Vec<u8> = Vec::with_capacity(n);
    let mut dchains = Chains::new(tbits, depth);
    let mut dins = k;
    let mut drh = 0u64;
    let mut ctx_inited = false;
    let mut dmodel = LitModel::reset(orderk);
    let mut rdec = RangeDecoder::new(&range_blob);
    let mut fi = 0usize;
    let mut ri = 0usize;
    let mut li = 0usize;
    let mut litc = 0usize;
    let read_flag = |flags: &[u8], fi: &mut usize| -> u8 {
        let byte = flags[*fi / 8];
        let bit = (byte >> (*fi % 8)) & 1;
        *fi += 1;
        bit
    };
    let get_varint = |v: &[u8], li: &mut usize| -> usize {
        let mut x = 0usize;
        let mut shift = 0;
        loop {
            let b = v[*li];
            *li += 1;
            x |= ((b & 0x7f) as usize) << shift;
            if b & 0x80 == 0 {
                break;
            }
            shift += 7;
        }
        x
    };
    let dec_lit = |rdec: &mut RangeDecoder, m: &mut LitModel, ctx: usize, b1: usize| -> u8 {
        if rdec.dec_bit(&mut m.escape[b1]) == 1 {
            return rdec.decode_raw(&mut m.raw);
        }
        if rdec.dec_bit(&mut m.isn[ctx]) == 1 {
            return b'N';
        }
        let bit0 = rdec.dec_bit(&mut m.t0[ctx]);
        let bit1 = rdec.dec_bit(&mut m.t1[ctx][bit0 as usize]);
        INV_BASE[((bit0 << 1) | bit1) as usize]
    };

    while out.len() < n {
        let p = out.len();
        if !ctx_inited && p >= k {
            drh = roller.at(&out, k);
            ctx_inited = true;
        }
        while dins < p {
            let idx = dchains.idx(drh);
            dchains.insert(idx, dins as u32);
            drh = roller.roll(drh, out[dins - k], out[dins]);
            dins += 1;
        }
        let fl = read_flag(&flags, &mut fi);
        if fl == 1 {
            let rank = ranks[ri] as usize;
            ri += 1;
            let len = get_varint(&lens, &mut li) + min_match;
            let idx = dchains.idx(drh);
            let q = dchains.candidate(idx, rank).expect("decode: rank out of range") as usize;
            for i in 0..len {
                let b = out[q + i];
                out.push(b);
            }
        } else if use_range {
            let (ctx, b1) = lit_ctx(&out, p, dmodel.mask);
            let b = dec_lit(&mut rdec, &mut dmodel, ctx, b1);
            out.push(b);
        } else {
            out.push(lits[litc]);
            litc += 1;
        }
    }
    let dec_s = td.elapsed().as_secs_f64();
    let lossless = out == data;

    let traw = Instant::now();
    let raw_bsc = bsc::compress_parallel(&data).expect("bsc raw").len();
    let raw_s = traw.elapsed().as_secs_f64();

    let cov = 100.0 * matched_bytes as f64 / n as f64;
    let avg_ml = if n_match > 0 { matched_bytes as f64 / n_match as f64 } else { 0.0 };
    eprintln!("[tokens] matches={n_match} lits={n_lit} coverage={cov:.1}%  avg_match_len={avg_ml:.1}");
    eprintln!(
        "[streams] flags={cf} ranks={cr} lens={cl} lits={cli} ({})",
        if use_range { "RANGE" } else { "BSC" }
    );
    eprintln!("[time] encode {enc_s:.1}s  decode {dec_s:.1}s  roundtrip={}", if lossless { "OK" } else { "FAIL" });
    if diag {
        let lit_total: u64 = lit_o0.iter().sum();
        let mut h0 = 0.0f64;
        for &c in lit_o0.iter() {
            if c > 0 {
                let pr = c as f64 / lit_total as f64;
                h0 -= pr * pr.log2();
            }
        }
        let me_frac = 100.0 * matched_eligible as f64 / n_lit.max(1) as f64;
        eprintln!(
            "[diag] lit_order0_entropy={h0:.3} bits/lit  matched_eligible={me_frac:.1}% of lits  rep0_runlen_hist[0,1-3,4-7,8-15,16-31,32+]={:?}",
            rep0_hist
        );
    }
    assert!(lossless, "ROLZ roundtrip MISMATCH — not lossless");
    println!(
        "ROLZ n={n} bpb={bpb:.3}  total={total}  vs_BSC={:+.1}%  BSC_raw={raw_bsc}({:.3}bpb,{raw_s:.1}s)  cov={cov:.1}%  litcoder={} K={k} DEPTH={depth} MIN={min_match}",
        100.0 * (total as f64 - raw_bsc as f64) / raw_bsc as f64,
        raw_bsc as f64 * 8.0 / n as f64,
        if use_range { "range" } else { "bsc" },
    );
}

// ===================== ROLZ match-finder primitives =====================

struct Roller {
    k: usize,
    pow: u64,
}
impl Roller {
    fn new(k: usize) -> Self {
        let mut pow = 1u64;
        for _ in 0..k - 1 {
            pow = pow.wrapping_mul(B);
        }
        Roller { k, pow }
    }
    fn at(&self, data: &[u8], pos: usize) -> u64 {
        let mut h = 0u64;
        for &b in &data[pos - self.k..pos] {
            h = h.wrapping_mul(B).wrapping_add(b as u64);
        }
        h
    }
    #[inline]
    fn roll(&self, h: u64, out_byte: u8, in_byte: u8) -> u64 {
        h.wrapping_sub((out_byte as u64).wrapping_mul(self.pow))
            .wrapping_mul(B)
            .wrapping_add(in_byte as u64)
    }
}

struct Chains {
    table: Vec<u32>,
    cnt: Vec<u32>,
    depth: usize,
    mask: usize,
}
impl Chains {
    fn new(tbits: usize, depth: usize) -> Self {
        let nb = 1usize << tbits;
        Chains { table: vec![0u32; nb * depth], cnt: vec![0u32; nb], depth, mask: nb - 1 }
    }
    #[inline]
    fn idx(&self, h: u64) -> usize {
        (h as usize) & self.mask
    }
    #[inline]
    fn insert(&mut self, idx: usize, pos: u32) {
        let c = self.cnt[idx] as usize;
        self.table[idx * self.depth + (c % self.depth)] = pos;
        self.cnt[idx] += 1;
    }
    #[inline]
    fn candidate(&self, idx: usize, rank: usize) -> Option<u32> {
        let c = self.cnt[idx] as usize;
        let live = c.min(self.depth);
        if rank >= live {
            return None;
        }
        let slot = (c - 1 - rank) % self.depth;
        Some(self.table[idx * self.depth + slot])
    }
}

#[inline]
fn lcp(data: &[u8], a: usize, b: usize, max: usize) -> usize {
    let n = data.len();
    let mut l = 0;
    while l < max && a + l < n && data[a + l] == data[b + l] {
        l += 1;
    }
    l
}
