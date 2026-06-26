//! Probe: can a tag-major (columnar) transpose of the aux stream beat bz's
//! current record-major raw+BSC? Parses BAM aux encoding, buckets each tag's
//! value bytes into its own stream, BSC-compresses each, and compares to the
//! current concatenated approach. Also reports per-tag sizes so we can see which
//! tags dominate and what dropping a derivable tag (e.g. NM) would save.
//!
//! Usage: cargo run --release -p bz-lib --example aux_probe -- <BAM> [N]

// Dev-only probe binary: relax style lints that don't matter for scratch tooling.
#![allow(unused_parens, clippy::unnecessary_sort_by)]

use bz_lib::compression::streams::write_varint;
use bz_lib::io::bam::RawBamReader;
use qz_lib::compression::bsc;
use std::collections::BTreeMap;
use std::path::Path;

/// Parse one BAM aux blob into (tag_key 3 bytes [t0,t1,type], value_bytes) items.
/// Returns None on malformed input.
fn parse_aux(mut a: &[u8], out: &mut Vec<([u8; 3], Vec<u8>)>) -> Option<()> {
    while !a.is_empty() {
        if a.len() < 3 {
            return None;
        }
        let key = [a[0], a[1], a[2]];
        let ty = a[2];
        a = &a[3..];
        let vlen = match ty {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            b'Z' | b'H' => {
                // include NUL
                (a.iter().position(|&b| b == 0)? + 1)
            }
            b'B' => {
                if a.len() < 5 {
                    return None;
                }
                let sub = a[0];
                let cnt = u32::from_le_bytes([a[1], a[2], a[3], a[4]]) as usize;
                let esz = match sub {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return None,
                };
                5 + cnt * esz
            }
            _ => return None,
        };
        if a.len() < vlen {
            return None;
        }
        out.push((key, a[..vlen].to_vec()));
        a = &a[vlen..];
    }
    Some(())
}

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: aux_probe <BAM> [N]");
    let n: usize = std::env::args()
        .nth(2)
        .map(|s| s.parse().unwrap())
        .unwrap_or(2_500_000);

    let mut reader = RawBamReader::from_path(Path::new(&path))?;
    let recs = reader.read_chunk(n)?;
    let nrec = recs.len();

    // Baseline: record-major varint+blob + BSC (mirrors bz today).
    let mut record_major = Vec::new();
    for r in &recs {
        let a = r.aux_bytes();
        write_varint(&mut record_major, a.len());
        record_major.extend_from_slice(a);
    }
    let base = bsc::compress_adaptive(&record_major)?;

    // Tag-major: one value-byte stream per tag key + a per-record layout stream
    // (sequence of tag keys present, so decode can re-interleave).
    let mut tag_streams: BTreeMap<[u8; 3], Vec<u8>> = BTreeMap::new();
    let mut layout = Vec::new();
    let mut items = Vec::new();
    let mut malformed = 0usize;
    for r in &recs {
        items.clear();
        if parse_aux(r.aux_bytes(), &mut items).is_none() {
            malformed += 1;
            // Fall back: stuff whole blob under a sentinel key.
            tag_streams
                .entry(*b"\0\0\0")
                .or_default()
                .extend_from_slice(r.aux_bytes());
            layout.push(0xFF);
            continue;
        }
        write_varint(&mut layout, items.len());
        for (key, val) in &items {
            layout.extend_from_slice(key);
            tag_streams.entry(*key).or_default().extend_from_slice(val);
        }
    }

    let layout_bsc = bsc::compress_adaptive(&layout)?;
    let mut col_total = layout_bsc.len();
    let mut per_tag: Vec<(String, usize, usize)> = Vec::new();
    let mut concat = Vec::new(); // tag-major order, single BSC stream
    for (key, stream) in &tag_streams {
        let c = bsc::compress_adaptive(stream)?.len();
        col_total += c;
        concat.extend_from_slice(stream);
        let name = format!("{}{}:{}", key[0] as char, key[1] as char, key[2] as char);
        per_tag.push((name, stream.len(), c));
    }
    let single = bsc::compress_adaptive(&concat)?.len() + layout_bsc.len();

    let mb = |x: usize| x as f64 / 1e6;
    println!("records           = {nrec}  (malformed aux = {malformed})");
    println!("record-major +BSC = {:.2} MB   [bz today]", mb(base.len()));
    println!(
        "tag-major perTag  = {:.2} MB   ({:.2}x vs baseline)",
        mb(col_total),
        base.len() as f64 / col_total as f64
    );
    println!(
        "tag-major 1-stream= {:.2} MB   ({:.2}x vs baseline)",
        mb(single),
        base.len() as f64 / single as f64
    );
    println!("  layout stream   = {:.2} MB", mb(layout_bsc.len()));
    per_tag.sort_by(|a, b| b.2.cmp(&a.2));
    println!("  per-tag (raw -> bsc):");
    for (name, raw, c) in &per_tag {
        println!("    {name:<6} {:>8.2} MB -> {:>7.3} MB", mb(*raw), mb(*c));
    }
    Ok(())
}
