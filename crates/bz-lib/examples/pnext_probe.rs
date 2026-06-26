//! Feasibility probe for mate-derived PNEXT (next_pos) + RNEXT (next_ref_id).
//!
//! For a primary, paired read whose mate is in the same chunk, the SAM spec says
//! PNEXT == mate.POS and RNEXT == mate.ref_id EXACTLY (no tie/sign convention,
//! unlike TLEN). The in-chunk mate map is already rebuilt deterministically on
//! decode (md_codec::build_mate_map), so we can drop the derivable values and
//! reconstruct them, storing only the non-derivable minority.
//!
//! This probe measures, per chunk:
//!   - in-chunk-mate fraction
//!   - byte-exact derivability of next_pos / next_ref_id for mate-in-chunk reads
//!   - CURRENT next_pos + next_ref_id stream size (production delta-encoding + BSC)
//!   - PROJECTED size after derivation (derivable bit + exceptions only, + BSC)
//!
//! Usage: cargo run --release -p bz-lib --example pnext_probe -- <BAM> [N]

use bz_lib::compression::md_codec;
use bz_lib::io::bam::{RawBamReader, RawBamRecord};
use qz_lib::compression::bsc;
use std::path::Path;

fn strip_nul(nm: &[u8]) -> &[u8] {
    if nm.last() == Some(&0) {
        &nm[..nm.len() - 1]
    } else {
        nm
    }
}

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: pnext_probe <BAM> [N]");
    let n: usize = std::env::args()
        .nth(2)
        .map(|s| s.parse().unwrap())
        .unwrap_or(2_500_000);
    let mut reader = RawBamReader::from_path(Path::new(&path))?;
    let recs: Vec<RawBamRecord> = reader.read_chunk(n)?;

    // Rebuild the SAME mate map production uses (deterministic from names+flags).
    let names: Vec<&[u8]> = recs.iter().map(|r| strip_nul(r.read_name())).collect();
    let flags: Vec<u16> = recs.iter().map(|r| r.flag()).collect();
    let mate = md_codec::build_mate_map(&names, &flags);

    let mut mate_in_chunk = 0usize;
    let mut np_derivable = 0usize; // next_pos == mate.pos
    let mut nr_derivable = 0usize; // next_ref_id == mate.ref_id
    let mut both_derivable = 0usize;

    // CURRENT production encoding (delta vs own pos/ref_id, 4 bytes each).
    let mut cur_next_pos = Vec::with_capacity(recs.len() * 4);
    let mut cur_next_ref = Vec::with_capacity(recs.len() * 4);
    // PROJECTED encoding: per-record derivable bit + exception values only.
    let mut derivable_bits = vec![0u8; recs.len().div_ceil(8)];
    let mut exc_next_pos = Vec::new();
    let mut exc_next_ref = Vec::new();

    for (i, r) in recs.iter().enumerate() {
        let ref_id = r.ref_id();
        let pos = r.pos();
        let nref = r.next_ref_id();
        let npos = r.next_pos();

        // Current production streams.
        cur_next_ref.extend_from_slice(&nref.wrapping_sub(ref_id).to_le_bytes());
        cur_next_pos.extend_from_slice(&npos.wrapping_sub(pos).to_le_bytes());

        let mi = mate[i];
        let derivable = if mi != u32::MAX {
            mate_in_chunk += 1;
            let m = &recs[mi as usize];
            let np_ok = npos == m.pos();
            let nr_ok = nref == m.ref_id();
            if np_ok {
                np_derivable += 1;
            }
            if nr_ok {
                nr_derivable += 1;
            }
            if np_ok && nr_ok {
                both_derivable += 1;
            }
            np_ok && nr_ok
        } else {
            false
        };

        if derivable {
            derivable_bits[i / 8] |= 1 << (i % 8);
        } else {
            // Store non-derivable values (delta vs own, same encoding as today).
            exc_next_ref.extend_from_slice(&nref.wrapping_sub(ref_id).to_le_bytes());
            exc_next_pos.extend_from_slice(&npos.wrapping_sub(pos).to_le_bytes());
        }
    }

    let cur_np = bsc::compress_adaptive(&cur_next_pos)?.len();
    let cur_nr = bsc::compress_adaptive(&cur_next_ref)?.len();
    let new_bits = bsc::compress_adaptive(&derivable_bits)?.len();
    let new_np = bsc::compress_adaptive(&exc_next_pos)?.len();
    let new_nr = bsc::compress_adaptive(&exc_next_ref)?.len();

    let mb = |x: usize| x as f64 / 1e6;
    let pct = |a: usize, b: usize| {
        if b > 0 {
            100.0 * a as f64 / b as f64
        } else {
            0.0
        }
    };

    println!("records              = {}", recs.len());
    println!(
        "mate in-chunk        = {mate_in_chunk} ({:.1}%)",
        pct(mate_in_chunk, recs.len())
    );
    println!("--- derivability (of mate-in-chunk reads) ---");
    println!(
        "  next_pos == mate.pos    = {np_derivable} ({:.2}%)",
        pct(np_derivable, mate_in_chunk)
    );
    println!(
        "  next_ref == mate.ref_id = {nr_derivable} ({:.2}%)",
        pct(nr_derivable, mate_in_chunk)
    );
    println!(
        "  BOTH derivable          = {both_derivable} ({:.2}%)",
        pct(both_derivable, mate_in_chunk)
    );
    println!(
        "  => exceptions (not derivable) = {} ({:.2}% of all reads)",
        recs.len() - both_derivable,
        pct(recs.len() - both_derivable, recs.len())
    );
    println!("--- size (BSC, per chunk) ---");
    println!(
        "  CURRENT   next_pos = {:.3} MB   next_ref_id = {:.3} MB   (sum {:.3} MB)",
        mb(cur_np),
        mb(cur_nr),
        mb(cur_np + cur_nr)
    );
    println!(
        "  PROJECTED derivable-bits = {:.3} MB   exc_next_pos = {:.3} MB   exc_next_ref = {:.3} MB   (sum {:.3} MB)",
        mb(new_bits),
        mb(new_np),
        mb(new_nr),
        mb(new_bits + new_np + new_nr)
    );
    let saved = (cur_np + cur_nr) as i64 - (new_bits + new_np + new_nr) as i64;
    println!("  SAVED   = {:.3} MB", saved as f64 / 1e6);
    Ok(())
}
