//! Probe: does qz's columnar header codec beat bz's current raw+BSC on real BAM
//! read names (QNAMEs)? Reads the first N records, strips the trailing NUL from
//! each QNAME, and compares:
//!   - raw stream (varint len + name) + BSC  [what bz does today]
//!   - header_col::compress_headers_columnar_bytes (qz tokenizer)
//! Also verifies header_col round-trips the names losslessly.
//!
//! Usage: cargo run --release -p bz-lib --example readname_probe -- <BAM> [N]

// Dev-only probe binary: relax doc-formatting lints on the free-form notes above.
#![allow(clippy::doc_lazy_continuation)]

use bz_lib::compression::streams::write_varint;
use bz_lib::io::bam::RawBamReader;
use qz_lib::compression::{bsc, header_col};
use std::path::Path;

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: readname_probe <BAM> [N]");
    let n: usize = std::env::args()
        .nth(2)
        .map(|s| s.parse().unwrap())
        .unwrap_or(2_500_000);

    let mut reader = RawBamReader::from_path(Path::new(&path))?;
    let recs = reader.read_chunk(n)?;

    let names: Vec<Vec<u8>> = recs
        .iter()
        .map(|r| {
            let nm = r.read_name();
            let nm = if nm.last() == Some(&0) {
                &nm[..nm.len() - 1]
            } else {
                nm
            };
            nm.to_vec()
        })
        .collect();
    let refs: Vec<&[u8]> = names.iter().map(|v| v.as_slice()).collect();

    // Baseline: raw varint-prefixed stream + BSC (mirrors bz's read_name path).
    let mut raw = Vec::new();
    for nm in &names {
        write_varint(&mut raw, nm.len());
        raw.extend_from_slice(nm);
    }
    let raw_bsc = bsc::compress(&raw)?;
    let raw_bsc_adaptive = bsc::compress_adaptive(&raw)?;

    // Candidate: qz columnar header codec.
    let hc = header_col::compress_headers_columnar_bytes(&refs)?;

    // Lossless check for header_col.
    let back = header_col::decompress_headers_columnar(&hc, names.len())?;
    let lossless = back.len() == names.len() && back.iter().zip(&names).all(|(a, b)| a == b);

    let mb = |x: usize| x as f64 / 1e6;
    println!("reads            = {}", names.len());
    println!("raw stream       = {:.2} MB", mb(raw.len()));
    println!("raw+BSC          = {:.2} MB", mb(raw_bsc.len()));
    println!("raw+BSC adaptive = {:.2} MB", mb(raw_bsc_adaptive.len()));
    println!(
        "header_col       = {:.2} MB   ({:.2}x vs raw+BSC, lossless={})",
        mb(hc.len()),
        raw_bsc.len() as f64 / hc.len() as f64,
        lossless
    );
    Ok(())
}
