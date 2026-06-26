//! Feasibility probe for MD-derived NM:i and UQ:i aux tags.
//!
//! Both are functions of data bz already stores (MD/consensus + CIGAR + quality):
//!   - NM = (# MD substitution letters) + (CIGAR inserted bases) + (deleted bases)
//!   - UQ = sum of base qualities at the MD substitution positions
//! If they regenerate byte-exact we can strip them on compress and rebuild on
//! decode (like MD/MC). The catch is the BAM integer *type byte*: aux ints use
//! the smallest type that fits (c/C/s/S/i/I), so to reproduce the stored bytes we
//! must reproduce both the value AND the source type. This probe measures the
//! derivable fraction, the type reproducibility (is the stored type the minimal
//! fit for the value?), and the current NM/UQ stream sizes.
//!
//! Usage: cargo run --release -p bz-lib --example nm_uq_probe -- <BAM> [N]

// Dev-only probe binary: relax style lints that don't matter for scratch tooling.
#![allow(unused_assignments, clippy::doc_lazy_continuation)]

use bz_lib::io::bam::{RawBamReader, RawBamRecord};
use qz_lib::compression::bsc;
use std::path::Path;

/// Find a 2-letter aux tag; return (type_byte, signed_value) for integer types.
fn int_tag(rec: &RawBamRecord, t0: u8, t1: u8) -> Option<(u8, i64)> {
    let mut a = rec.aux_bytes();
    while a.len() >= 3 {
        let (k0, k1, ty) = (a[0], a[1], a[2]);
        a = &a[3..];
        let vlen = match ty {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            b'Z' | b'H' => a.iter().position(|&b| b == 0)? + 1,
            b'B' => {
                let cnt = u32::from_le_bytes([a[1], a[2], a[3], a[4]]) as usize;
                let esz = match a[0] {
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
        if k0 == t0 && k1 == t1 {
            let v: i64 = match ty {
                b'c' => (a[0] as i8) as i64,
                b'C' => a[0] as i64,
                b's' => i16::from_le_bytes([a[0], a[1]]) as i64,
                b'S' => u16::from_le_bytes([a[0], a[1]]) as i64,
                b'i' => i32::from_le_bytes([a[0], a[1], a[2], a[3]]) as i64,
                b'I' => u32::from_le_bytes([a[0], a[1], a[2], a[3]]) as i64,
                _ => return None,
            };
            return Some((ty, v));
        }
        a = &a[vlen..];
    }
    None
}

const CIG_M: u8 = 0;
const CIG_I: u8 = 1;
const CIG_D: u8 = 2;
const CIG_N: u8 = 3;
const CIG_S: u8 = 4;
const CIG_EQ: u8 = 7;
const CIG_X: u8 = 8;

/// Derive (NM, UQ) from CIGAR + MD + query SEQ (ASCII) + query QUAL (Phred).
/// Returns None if MD/CIGAR disagree (caller treats as non-derivable).
fn derive_nm_uq(cigar: &[(u8, u32)], md: &[u8], qual: &[u8]) -> Option<(i64, i64)> {
    let mut subs: i64 = 0;
    let mut ins: i64 = 0;
    let mut del: i64 = 0;
    let mut uq: i64 = 0;
    let mut read_idx = 0usize;
    let mut mi = 0usize;
    let mut match_left: i64 = -1;

    fn next_num(md: &[u8], mi: &mut usize) -> i64 {
        let mut v: i64 = 0;
        let mut any = false;
        while *mi < md.len() && md[*mi].is_ascii_digit() {
            v = v * 10 + (md[*mi] - b'0') as i64;
            *mi += 1;
            any = true;
        }
        if any { v } else { -1 }
    }

    for &(op, len) in cigar {
        let len = len as usize;
        match op {
            CIG_M | CIG_EQ | CIG_X => {
                for _ in 0..len {
                    if match_left < 0 {
                        match_left = next_num(md, &mut mi).max(0);
                    }
                    if match_left > 0 {
                        match_left -= 1;
                    } else {
                        if mi >= md.len() || !md[mi].is_ascii_alphabetic() {
                            return None;
                        }
                        mi += 1;
                        match_left = -1;
                        subs += 1;
                        if read_idx >= qual.len() {
                            return None;
                        }
                        uq += qual[read_idx] as i64;
                    }
                    read_idx += 1;
                }
            }
            CIG_I => {
                ins += len as i64;
                read_idx += len;
            }
            CIG_S => read_idx += len,
            CIG_D => {
                if match_left < 0 {
                    match_left = next_num(md, &mut mi).max(0);
                }
                if mi >= md.len() || md[mi] != b'^' {
                    return None;
                }
                mi += 1;
                for _ in 0..len {
                    if mi >= md.len() || !md[mi].is_ascii_alphabetic() {
                        return None;
                    }
                    mi += 1;
                    del += 1;
                }
                match_left = -1;
            }
            CIG_N => {}
            _ => {}
        }
    }
    Some((subs + ins + del, uq))
}

/// Minimal BAM integer type byte for a value (htslib bam_aux_append convention:
/// unsigned types for non-negatives, signed for negatives, smallest that fits).
fn minimal_int_type(v: i64) -> u8 {
    if v >= 0 {
        if v <= 0xff {
            b'C'
        } else if v <= 0xffff {
            b'S'
        } else {
            b'I'
        }
    } else if v >= -128 {
        b'c'
    } else if v >= -32768 {
        b's'
    } else {
        b'i'
    }
}

fn md_of(rec: &RawBamRecord) -> Option<&[u8]> {
    let mut a = rec.aux_bytes();
    while a.len() >= 3 {
        let (k0, k1, ty) = (a[0], a[1], a[2]);
        a = &a[3..];
        let vlen = match ty {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            b'Z' | b'H' => a.iter().position(|&b| b == 0)? + 1,
            b'B' => {
                let cnt = u32::from_le_bytes([a[1], a[2], a[3], a[4]]) as usize;
                let esz = match a[0] {
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
        if k0 == b'M' && k1 == b'D' && ty == b'Z' {
            return Some(&a[..vlen - 1]); // strip NUL
        }
        a = &a[vlen..];
    }
    None
}

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: nm_uq_probe <BAM> [N]");
    let n: usize = std::env::args()
        .nth(2)
        .map(|s| s.parse().unwrap())
        .unwrap_or(2_500_000);
    let mut reader = RawBamReader::from_path(Path::new(&path))?;
    let recs = reader.read_chunk(n)?;

    let mut nm_count = 0usize;
    let mut nm_have_md = 0usize;
    let mut nm_val_ok = 0usize;
    let mut nm_type_minimal = 0usize;
    let mut nm_full_ok = 0usize; // value + type both reproducible
    let mut nm_type_hist: std::collections::BTreeMap<u8, usize> = std::collections::BTreeMap::new();
    let mut nm_stream = Vec::new();
    let mut nm_examples = 0usize;

    let mut uq_count = 0usize;
    let mut uq_val_ok = 0usize;
    let mut uq_full_ok = 0usize;
    let mut uq_type_hist: std::collections::BTreeMap<u8, usize> = std::collections::BTreeMap::new();
    let mut uq_stream = Vec::new();
    let mut uq_examples = 0usize;

    for r in &recs {
        let cigar = r.cigar_ops();
        let qual = r.qual_bytes();
        let md = md_of(r);
        let derived = md.and_then(|md| derive_nm_uq(&cigar, md, qual));

        if let Some((ty, v)) = int_tag(r, b'N', b'M') {
            nm_count += 1;
            *nm_type_hist.entry(ty).or_default() += 1;
            nm_stream.extend_from_slice(&v.to_le_bytes());
            if md.is_some() {
                nm_have_md += 1;
            }
            if let Some((dnm, _)) = derived {
                if dnm == v {
                    nm_val_ok += 1;
                    if minimal_int_type(v) == ty {
                        nm_type_minimal += 1;
                        nm_full_ok += 1;
                    } else if nm_examples < 6 {
                        eprintln!(
                            "  NM type mismatch: stored {} val={v} minimal={}",
                            ty as char,
                            minimal_int_type(v) as char
                        );
                        nm_examples += 1;
                    }
                } else if nm_examples < 6 {
                    eprintln!("  NM value mismatch: derived={dnm} stored={v}");
                    nm_examples += 1;
                }
            }
        }

        if let Some((ty, v)) = int_tag(r, b'U', b'Q') {
            uq_count += 1;
            *uq_type_hist.entry(ty).or_default() += 1;
            uq_stream.extend_from_slice(&v.to_le_bytes());
            if let Some((_, duq)) = derived {
                if duq == v {
                    uq_val_ok += 1;
                    if minimal_int_type(v) == ty {
                        uq_full_ok += 1;
                    }
                } else if uq_examples < 6 {
                    eprintln!("  UQ value mismatch: derived={duq} stored={v}");
                    uq_examples += 1;
                }
            }
        }
    }

    let mb = |x: usize| x as f64 / 1e6;
    let pct = |a: usize, b: usize| {
        if b > 0 {
            100.0 * a as f64 / b as f64
        } else {
            0.0
        }
    };
    let hist = |h: &std::collections::BTreeMap<u8, usize>| {
        h.iter()
            .map(|(t, c)| format!("{}:{}", *t as char, c))
            .collect::<Vec<_>>()
            .join(" ")
    };

    println!("records              = {}", recs.len());
    println!("--- NM:i ---");
    println!("  reads with NM      = {nm_count}  (have MD: {nm_have_md})");
    println!("  type histogram     = [{}]", hist(&nm_type_hist));
    println!(
        "  derived == NM      = {nm_val_ok} ({:.2}% of NM)",
        pct(nm_val_ok, nm_count)
    );
    println!(
        "  + type=minimal-fit = {nm_type_minimal} ({:.2}% of NM)",
        pct(nm_type_minimal, nm_count)
    );
    println!(
        "  FULL byte-exact    = {nm_full_ok} ({:.2}% of NM)",
        pct(nm_full_ok, nm_count)
    );
    if nm_count > 0 {
        println!(
            "  NM raw+BSC         = {:.3} MB [droppable for derivable reads]",
            mb(bsc::compress_adaptive(&nm_stream)?.len())
        );
    }
    println!("--- UQ:i ---");
    println!("  reads with UQ      = {uq_count}");
    println!("  type histogram     = [{}]", hist(&uq_type_hist));
    println!(
        "  derived == UQ      = {uq_val_ok} ({:.2}% of UQ)",
        pct(uq_val_ok, uq_count)
    );
    println!(
        "  FULL byte-exact    = {uq_full_ok} ({:.2}% of UQ)",
        pct(uq_full_ok, uq_count)
    );
    if uq_count > 0 {
        println!(
            "  UQ raw+BSC         = {:.3} MB [droppable for derivable reads]",
            mb(bsc::compress_adaptive(&uq_stream)?.len())
        );
    }
    Ok(())
}
