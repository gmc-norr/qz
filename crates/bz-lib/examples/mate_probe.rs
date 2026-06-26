//! Feasibility probe for mate-derived aux/fields: MC:Z (mate CIGAR) and TLEN.
//!
//! Both are functions of the *mate* record, which (in a coordinate-sorted BAM)
//! is usually but not always in the same chunk. This probe pairs primary reads
//! by QNAME within a chunk and checks:
//!   - MC:Z == textual CIGAR of the mate (byte-exact)?
//!   - TLEN == derived signed template length (value-exact, incl. aligner sign
//!     convention on ties)?
//! and reports the in-chunk-mate fraction and the current MC/TLEN stream sizes.
//!
//! Usage: cargo run --release -p bz-lib --example mate_probe -- <BAM> [N]

// Dev-only probe binary: the `tieA`/`tieB` names and the two deliberately-identical
// tie-convention branches document the comparison, so relax those style lints.
#![allow(
    non_snake_case,
    clippy::if_same_then_else,
    clippy::doc_lazy_continuation
)]

use bz_lib::compression::streams::{write_varint, zigzag_encode};
use bz_lib::io::bam::{RawBamReader, RawBamRecord};
use qz_lib::compression::bsc;
use std::collections::HashMap;
use std::path::Path;

const CONSUMES_REF: [bool; 9] = [
    true,  // M
    false, // I
    true,  // D
    true,  // N
    false, // S
    false, // H
    false, // P
    true,  // =
    true,  // X
];

fn ref_span(cigar: &[(u8, u32)]) -> i32 {
    cigar
        .iter()
        .filter(|(op, _)| CONSUMES_REF[*op as usize])
        .map(|(_, l)| *l as i32)
        .sum()
}

fn cigar_text(cigar: &[(u8, u32)]) -> Vec<u8> {
    const OPS: &[u8; 9] = b"MIDNSHP=X";
    let mut out = Vec::new();
    if cigar.is_empty() {
        out.push(b'*');
        return out;
    }
    for &(op, len) in cigar {
        out.extend_from_slice(len.to_string().as_bytes());
        out.push(OPS[op as usize]);
    }
    out
}

fn tag_value(rec: &RawBamRecord, t0: u8, t1: u8) -> Option<(&[u8], u8)> {
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
            return Some((&a[..vlen], ty));
        }
        a = &a[vlen..];
    }
    None
}

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: mate_probe <BAM> [N]");
    let n: usize = std::env::args()
        .nth(2)
        .map(|s| s.parse().unwrap())
        .unwrap_or(2_500_000);
    let mut reader = RawBamReader::from_path(Path::new(&path))?;
    let recs = reader.read_chunk(n)?;

    // Pair primary, paired, mapped reads by QNAME within the chunk.
    let mut by_name: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for (i, r) in recs.iter().enumerate() {
        let f = r.flag();
        let primary = f & 0x900 == 0; // not secondary/supplementary
        let paired = f & 0x1 != 0;
        if primary && paired {
            let nm = r.read_name();
            let nm = if nm.last() == Some(&0) {
                &nm[..nm.len() - 1]
            } else {
                nm
            };
            by_name.entry(nm).or_default().push(i);
        }
    }
    // mate[i] = index of i's mate (read1<->read2) when both primary & in-chunk.
    let mut mate = vec![usize::MAX; recs.len()];
    for idxs in by_name.values() {
        if idxs.len() == 2 {
            let (a, b) = (idxs[0], idxs[1]);
            if (recs[a].flag() & 0x40) != (recs[b].flag() & 0x40) {
                mate[a] = b;
                mate[b] = a;
            }
        }
    }

    // ---- MC:Z derivability ----
    let mut mc_count = 0usize;
    let mut mc_mate_in_chunk = 0usize;
    let mut mc_derivable = 0usize;
    let mut mc_concat = Vec::new();
    // ---- TLEN derivability ----
    let mut tlen_nonzero = 0usize;
    let mut tlen_mate_in_chunk = 0usize;
    let mut tlen_derivable = 0usize;
    let mut tlen_sign_only = 0usize;
    let mut tlen_mag_diff = 0usize;
    let mut tieA_ok = 0usize;
    let mut tieB_ok = 0usize;
    let mut tlen_examples = 0usize;
    let mut tlen_stream = Vec::new();

    for (i, r) in recs.iter().enumerate() {
        // TLEN
        let tlen = r.tlen();
        tlen_stream.extend_from_slice(&zigzag_encode(tlen).to_le_bytes());
        if tlen != 0 {
            tlen_nonzero += 1;
            if mate[i] != usize::MAX {
                tlen_mate_in_chunk += 1;
                let m = &recs[mate[i]];
                if r.ref_id() == m.ref_id() {
                    let sa = r.pos();
                    let ea = sa + ref_span(&r.cigar_ops());
                    let sb = m.pos();
                    let eb = sb + ref_span(&m.cigar_ops());
                    let left = sa.min(sb);
                    let right = ea.max(eb);
                    let span = right - left;
                    let read1 = r.flag() & 0x40 != 0;
                    // two tie conventions for equal start positions
                    let d_tieA = if sa < sb {
                        span
                    } else if sa > sb {
                        -span
                    } else if read1 {
                        span
                    } else {
                        -span
                    };
                    let d_tieB = if sa < sb {
                        span
                    } else if sa > sb {
                        -span
                    } else if read1 {
                        -span
                    } else {
                        span
                    };
                    if d_tieA == tlen {
                        tieA_ok += 1;
                    }
                    if d_tieB == tlen {
                        tieB_ok += 1;
                    }
                    let derived = d_tieA;
                    if derived == tlen {
                        tlen_derivable += 1;
                    } else if derived == -tlen {
                        tlen_sign_only += 1;
                        if tlen_examples < 4 {
                            eprintln!(
                                "  sign-flip: derived={derived} tlen={tlen} sa={sa} sb={sb} ea={ea} eb={eb} read1={}",
                                r.flag() & 0x40 != 0
                            );
                            tlen_examples += 1;
                        }
                    } else {
                        tlen_mag_diff += 1;
                        if tlen_examples < 8 {
                            eprintln!(
                                "  mag-diff:  derived={derived} tlen={tlen} sa={sa} sb={sb} ea={ea} eb={eb}"
                            );
                            tlen_examples += 1;
                        }
                    }
                }
            }
        }
        // MC:Z
        if let Some((val, b'Z')) = tag_value(r, b'M', b'C') {
            mc_count += 1;
            let mc = &val[..val.len() - 1]; // strip NUL
            mc_concat.extend_from_slice(mc);
            mc_concat.push(0);
            if mate[i] != usize::MAX {
                mc_mate_in_chunk += 1;
                let derived = cigar_text(&recs[mate[i]].cigar_ops());
                if derived == mc {
                    mc_derivable += 1;
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
    println!("records              = {}", recs.len());
    println!("--- MC:Z ---");
    println!("  reads with MC      = {mc_count}");
    println!(
        "  mate in-chunk      = {mc_mate_in_chunk} ({:.1}%)",
        pct(mc_mate_in_chunk, mc_count)
    );
    println!(
        "  MC == mate CIGAR   = {mc_derivable} ({:.1}% of MC, {:.1}% of in-chunk)",
        pct(mc_derivable, mc_count),
        pct(mc_derivable, mc_mate_in_chunk)
    );
    if mc_count > 0 {
        println!(
            "  MC:Z raw+BSC       = {:.2} MB [droppable for derivable reads]",
            mb(bsc::compress_adaptive(&mc_concat)?.len())
        );
    }
    println!("--- TLEN ---");
    println!("  tlen != 0          = {tlen_nonzero}");
    println!(
        "  mate in-chunk      = {tlen_mate_in_chunk} ({:.1}%)",
        pct(tlen_mate_in_chunk, tlen_nonzero)
    );
    println!(
        "  derived == tlen    = {tlen_derivable} ({:.1}% of nonzero, {:.1}% of in-chunk)",
        pct(tlen_derivable, tlen_nonzero),
        pct(tlen_derivable, tlen_mate_in_chunk)
    );
    println!(
        "  sign-flip only     = {tlen_sign_only} ({:.1}% of in-chunk)",
        pct(tlen_sign_only, tlen_mate_in_chunk)
    );
    println!(
        "  magnitude diff     = {tlen_mag_diff} ({:.1}% of in-chunk)",
        pct(tlen_mag_diff, tlen_mate_in_chunk)
    );
    let best = tieA_ok.max(tieB_ok);
    println!(
        "  BEST tie-break     = {best} ({:.2}% of in-chunk)  [tieA={tieA_ok} tieB={tieB_ok}]",
        pct(best, tlen_mate_in_chunk)
    );
    println!(
        "  => exceptions      = {} ({:.2}% of nonzero) under best tie-break + cross-chunk",
        tlen_nonzero - best,
        pct(tlen_nonzero - best, tlen_nonzero)
    );
    println!(
        "  tlen raw+BSC       = {:.2} MB",
        mb(bsc::compress_adaptive(&tlen_stream)?.len())
    );
    let _ = write_varint; // (kept for parity with other probes)
    Ok(())
}
