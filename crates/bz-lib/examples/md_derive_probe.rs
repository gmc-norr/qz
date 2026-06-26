//! Feasibility spike for MD:Z derivation from the per-chunk consensus.
//!
//! Plan (self-contained, no external reference):
//!   reference base at a covered ref position is recoverable from MD+SEQ+CIGAR
//!   (MD match -> ref = read base; MD sub -> ref = MD letter; MD del -> ref =
//!   deleted letters). Compare to the consensus bz already stores; positions where
//!   ref != consensus form a small EXCEPTION map. On decode the reference =
//!   consensus + exceptions, and MD regenerates byte-exact from ref+CIGAR+SEQ.
//!
//! This probe verifies (1) byte-exact MD regeneration on real reads and (2) the
//! exception-stream size vs the current MD:Z raw+BSC size.
//!
//! Usage: cargo run --release -p bz-lib --example md_derive_probe -- <BAM> [N]

use bz_lib::compression::streams::{ChunkConsensus, write_varint};
use bz_lib::io::bam::{RawBamReader, RawBamRecord};
use qz_lib::compression::bsc;
use std::collections::BTreeMap;
use std::path::Path;

const CIG_M: u8 = 0;
const CIG_I: u8 = 1;
const CIG_D: u8 = 2;
const CIG_N: u8 = 3;
const CIG_S: u8 = 4;
const CIG_EQ: u8 = 7;
const CIG_X: u8 = 8;

fn nib_to_base(n: u8) -> u8 {
    match n {
        1 => b'A',
        2 => b'C',
        4 => b'G',
        8 => b'T',
        _ => b'N',
    }
}

/// Extract the MD:Z value (without trailing NUL) from a record's aux blob.
fn md_of(rec: &RawBamRecord) -> Option<Vec<u8>> {
    let mut a = rec.aux_bytes();
    while a.len() >= 3 {
        let (t0, t1, ty) = (a[0], a[1], a[2]);
        a = &a[3..];
        let vlen = match ty {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            b'Z' | b'H' => a.iter().position(|&b| b == 0)? + 1,
            b'B' => {
                let sub = a[0];
                let cnt = u32::from_le_bytes([a[1], a[2], a[3], a[4]]) as usize;
                let esz = match sub {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    _ => 4,
                };
                5 + cnt * esz
            }
            _ => return None,
        };
        if a.len() < vlen {
            return None;
        }
        if t0 == b'M' && t1 == b'D' && ty == b'Z' {
            return Some(a[..vlen - 1].to_vec()); // strip NUL
        }
        a = &a[vlen..];
    }
    None
}

/// Walk MD+CIGAR+SEQ to recover the reference base at each ref-consuming position.
/// Returns ref bases in ref order starting at rec.pos(), or None if MD/CIGAR
/// disagree (then this read can't be derived and keeps its MD).
fn recover_ref(rec: &RawBamRecord, md: &[u8], seq: &[u8]) -> Option<Vec<u8>> {
    let cigar = rec.cigar_ops();
    let mut out = Vec::new();
    let mut read_idx = 0usize;
    // MD tokenizer state
    let mut mi = 0usize;
    let mut match_left: i64 = -1; // -1 = need to read a number next
    let next_num = |mi: &mut usize| -> i64 {
        let mut v: i64 = 0;
        let mut any = false;
        while *mi < md.len() && md[*mi].is_ascii_digit() {
            v = v * 10 + (md[*mi] - b'0') as i64;
            *mi += 1;
            any = true;
        }
        if any { v } else { -1 }
    };
    macro_rules! ensure_num {
        () => {
            if match_left < 0 {
                match_left = next_num(&mut mi).max(0);
            }
        };
    }

    for (op, len) in cigar {
        let len = len as usize;
        match op {
            CIG_M | CIG_EQ | CIG_X => {
                for _ in 0..len {
                    ensure_num!();
                    if match_left > 0 {
                        // match: ref == read base
                        out.push(seq[read_idx]);
                        match_left -= 1;
                    } else {
                        // mismatch: MD letter is the ref base
                        if mi >= md.len() || !md[mi].is_ascii_alphabetic() {
                            return None;
                        }
                        out.push(md[mi]);
                        mi += 1;
                        match_left = -1; // a number (maybe 0) follows
                    }
                    read_idx += 1;
                }
            }
            CIG_I | CIG_S => read_idx += len,
            CIG_D => {
                ensure_num!();
                // a deletion in MD is ^ followed by `len` ref bases
                if match_left != 0 {
                    // MD should be at 0 before a deletion
                }
                if mi >= md.len() || md[mi] != b'^' {
                    return None;
                }
                mi += 1;
                for _ in 0..len {
                    if mi >= md.len() || !md[mi].is_ascii_alphabetic() {
                        return None;
                    }
                    out.push(md[mi]);
                    mi += 1;
                }
                match_left = -1;
            }
            CIG_N => {
                // ref skip: not represented in MD; emit N placeholders (rare in DNA)
                out.extend(std::iter::repeat_n(b'N', len));
            }
            _ => {}
        }
    }
    Some(out)
}

/// Regenerate MD from reference bases (ref order) + CIGAR + SEQ. Must match the
/// SAM spec regex `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`.
fn regen_md(rec: &RawBamRecord, ref_bases: &[u8], seq: &[u8]) -> Vec<u8> {
    let cigar = rec.cigar_ops();
    let mut md = Vec::new();
    let mut run: u64 = 0;
    let mut read_idx = 0usize;
    let mut ref_idx = 0usize;
    let push_num = |md: &mut Vec<u8>, n: u64| md.extend_from_slice(n.to_string().as_bytes());
    for (op, len) in cigar {
        let len = len as usize;
        match op {
            CIG_M | CIG_EQ | CIG_X => {
                for _ in 0..len {
                    let rb = ref_bases[ref_idx];
                    if seq[read_idx] == rb {
                        run += 1;
                    } else {
                        push_num(&mut md, run);
                        run = 0;
                        md.push(rb);
                    }
                    read_idx += 1;
                    ref_idx += 1;
                }
            }
            CIG_I | CIG_S => read_idx += len,
            CIG_D => {
                push_num(&mut md, run);
                run = 0;
                md.push(b'^');
                for _ in 0..len {
                    md.push(ref_bases[ref_idx]);
                    ref_idx += 1;
                }
            }
            CIG_N => ref_idx += len,
            _ => {}
        }
    }
    push_num(&mut md, run);
    md
}

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: md_derive_probe <BAM> [N]");
    let n: usize = std::env::args()
        .nth(2)
        .map(|s| s.parse().unwrap())
        .unwrap_or(2_500_000);
    let mut reader = RawBamReader::from_path(Path::new(&path))?;
    let recs = reader.read_chunk(n)?;
    let consensus = ChunkConsensus::build(&recs, None);

    let mut had_md = 0usize;
    let mut derive_fail = 0usize;
    let mut conflicts = 0usize;
    let mut md_raw_total = 0usize; // raw MD bytes (the stream we'd drop)
    let mut md_concat = Vec::new(); // for baseline MD:Z BSC sizing
    // exception map: (ref_id, ref_pos) -> ref base (BAM-position keyed)
    let mut exceptions: BTreeMap<(i32, i32), u8> = BTreeMap::new();

    // PASS 1: recover ref bases from each read's MD; build the global exception
    // map (positions where ref != consensus). Detect cross-read conflicts.
    for rec in &recs {
        let md = match md_of(rec) {
            Some(m) => m,
            None => continue,
        };
        had_md += 1;
        md_raw_total += md.len();
        md_concat.extend_from_slice(&md);
        md_concat.push(0);

        let nibbles = rec.unpack_seq_nibbles();
        let seq: Vec<u8> = nibbles.iter().map(|&nb| nib_to_base(nb)).collect();
        let ref_bases = match recover_ref(rec, &md, &seq) {
            Some(r) => r,
            None => {
                derive_fail += 1;
                continue;
            }
        };
        let ref_id = rec.ref_id();
        let start = rec.pos();
        if let Some(seg) = consensus.get_segment_for_pos(ref_id, start) {
            for (k, &rb) in ref_bases.iter().enumerate() {
                let p = start + k as i32;
                let idx = (p - seg.start_pos) as usize;
                let cons = seg
                    .bases
                    .get(idx)
                    .map(|&nb| nib_to_base(nb))
                    .unwrap_or(b'N');
                if cons != rb {
                    if let Some(&prev) = exceptions.get(&(ref_id, p))
                        && prev != rb
                    {
                        conflicts += 1;
                    }
                    exceptions.insert((ref_id, p), rb);
                }
            }
        }
    }

    // PASS 2 (the real decode path): reconstruct reference = consensus +
    // exceptions, regenerate MD from it, and verify byte-exact vs the original.
    let mut regen_ok = 0usize;
    let mut regen_bad = 0usize;
    for rec in &recs {
        let md = match md_of(rec) {
            Some(m) => m,
            None => continue,
        };
        let nibbles = rec.unpack_seq_nibbles();
        let seq: Vec<u8> = nibbles.iter().map(|&nb| nib_to_base(nb)).collect();
        // determine ref span length from CIGAR (ref-consuming ops)
        let mut span = 0usize;
        for (op, len) in rec.cigar_ops() {
            if matches!(op, CIG_M | CIG_D | CIG_N | CIG_EQ | CIG_X) {
                span += len as usize;
            }
        }
        let ref_id = rec.ref_id();
        let start = rec.pos();
        let seg = consensus.get_segment_for_pos(ref_id, start);
        let mut ref_bases = Vec::with_capacity(span);
        for k in 0..span {
            let p = start + k as i32;
            let b = if let Some(&e) = exceptions.get(&(ref_id, p)) {
                e
            } else if let Some(seg) = seg {
                seg.bases
                    .get((p - seg.start_pos) as usize)
                    .map(|&nb| nib_to_base(nb))
                    .unwrap_or(b'N')
            } else {
                b'N'
            };
            ref_bases.push(b);
        }
        let regen = regen_md(rec, &ref_bases, &seq);
        if regen == md {
            regen_ok += 1;
        } else {
            regen_bad += 1;
            if regen_bad <= 3 {
                eprintln!(
                    "  MD regen mismatch: orig={:?} regen={:?}",
                    String::from_utf8_lossy(&md),
                    String::from_utf8_lossy(&regen)
                );
            }
        }
    }

    // Size the exception stream: sorted (delta ref_id, delta pos, base).
    let mut exc = Vec::new();
    let (mut pr, mut pp) = (0i32, 0i32);
    for ((rid, pos), base) in &exceptions {
        write_varint(&mut exc, (*rid - pr).unsigned_abs() as usize);
        let dp = if *rid == pr { *pos - pp } else { *pos };
        write_varint(&mut exc, dp.max(0) as usize);
        exc.push(*base);
        pr = *rid;
        pp = *pos;
    }
    let exc_bsc = bsc::compress_adaptive(&exc)?.len();
    let md_bsc = bsc::compress_adaptive(&md_concat)?.len();
    let mb = |x: usize| x as f64 / 1e6;

    println!("records              = {}", recs.len());
    println!("reads with MD        = {had_md}");
    println!(
        "byte-exact MD regen  = {regen_ok} ok / {regen_bad} bad (END-TO-END via consensus+exceptions)"
    );
    println!("derive-fail / conflicts = {derive_fail} / {conflicts}");
    println!("MD raw               = {:.2} MB", mb(md_raw_total));
    println!(
        "MD:Z raw+BSC (today) = {:.2} MB   [the stream we'd drop]",
        mb(md_bsc)
    );
    println!("exceptions           = {} entries", exceptions.len());
    println!("exception stream+BSC = {:.3} MB", mb(exc_bsc));
    if md_bsc > 0 {
        println!(
            "=> MD saving         = {:.2} MB ({:.1}% of MD)",
            mb(md_bsc.saturating_sub(exc_bsc)),
            100.0 * (md_bsc.saturating_sub(exc_bsc)) as f64 / md_bsc as f64
        );
    }
    Ok(())
}
