//! Spike: would a BOUNDED-STREAMING fqz quality decode erase the O(reads)
//! decompress-RAM penalty measured for the (sorted) deployed-fqz path?
//!
//! The deployed-fqz path can't stream because it GLOBALLY sorts reads by mean
//! quality — read i's quality lives in sub-chunk rank(i), so decode must hold
//! every block + the permutation to scatter back to input order. Drop the sort
//! and use INPUT-ORDER 500K sub-blocks (what reference mode's
//! `fqz_blocks_capped_quals` already builds) and each block decodes to a
//! contiguous run of original-order reads → streamable.
//!
//! This isolates the quality-stream decode (no sequences/headers), so the RSS
//! here is the *delta* the real decompressor would see on top of its seq+output
//! baseline. Run encode once, then decode in two SEPARATE processes for clean
//! peak-RSS attribution:
//!
//!   bench_fqz_streaming_decode encode <out.qbin> --fastq R1 [R2]
//!   /usr/bin/time -v bench_fqz_streaming_decode decode <out.qbin> allatonce <sink>
//!   /usr/bin/time -v bench_fqz_streaming_decode decode <out.qbin> streaming <sink>
//!
//! allatonce = decode every block, hold ALL decoded qualities live at once
//!             (mimics today's materialize-everything fqz decode).
//! streaming = decode block, write it out, drop it (bounded working set).
//! Both write identical bytes to <sink> (use /dev/null to drop, or a path to
//! verify equality with `cmp`). Internal decode wall-time is printed.

use anyhow::{Context, Result, bail};
use qz_lib::compression::fqzcomp;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::time::Instant;

const SUB_BLOCK: usize = 500_000;

fn trim(v: &mut Vec<u8>) {
    while matches!(v.last(), Some(b'\n' | b'\r')) {
        v.pop();
    }
}

fn read_quals(path: &str, out: &mut Vec<Vec<u8>>) -> Result<()> {
    let f = std::fs::File::open(path).with_context(|| format!("open {path}"))?;
    let mut r = BufReader::with_capacity(1 << 22, f);
    let mut line = Vec::new();
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
            break; // @id EOF
        }
        let mut seq = Vec::new();
        if r.read_until(b'\n', &mut seq)? == 0 {
            break;
        }
        line.clear();
        r.read_until(b'\n', &mut line)?; // +
        let mut qual = Vec::new();
        if r.read_until(b'\n', &mut qual)? == 0 {
            break;
        }
        trim(&mut qual);
        out.push(qual);
    }
    Ok(())
}

// Blob file: [magic 4]["FQS1"][num_blocks u32][total_reads u64]
//   then per block: [n_reads u32][blob_len u32][blob bytes]
fn cmd_encode(out_path: &str, fastqs: &[String]) -> Result<()> {
    let mut quals: Vec<Vec<u8>> = Vec::new();
    let t = Instant::now();
    for f in fastqs {
        read_quals(f, &mut quals)?;
    }
    let total = quals.len();
    let total_bytes: usize = quals.iter().map(|q| q.len()).sum();
    eprintln!(
        "encode: loaded {total} reads ({total_bytes} qual bytes) in {:.1}s",
        t.elapsed().as_secs_f64()
    );

    // INPUT-ORDER 500K sub-blocks, no sort.
    use rayon::prelude::*;
    let num_blocks = total.div_ceil(SUB_BLOCK);
    let t = Instant::now();
    let blocks: Vec<(usize, Vec<u8>)> = (0..num_blocks)
        .into_par_iter()
        .map(|i| {
            let s = i * SUB_BLOCK;
            let e = (s + SUB_BLOCK).min(total);
            let refs: Vec<&[u8]> = quals[s..e].iter().map(|q| q.as_slice()).collect();
            let blob = fqzcomp::compress(&refs, 0)?;
            Ok::<_, anyhow::Error>((e - s, blob))
        })
        .collect::<Result<Vec<_>>>()?;
    let comp: usize = blocks.iter().map(|(_, b)| b.len()).sum();
    eprintln!(
        "encode: {num_blocks} input-order blocks, {comp} compressed bytes in {:.1}s ({:.4} bits/q)",
        t.elapsed().as_secs_f64(),
        comp as f64 * 8.0 / total_bytes as f64
    );

    let f = std::fs::File::create(out_path)?;
    let mut w = BufWriter::with_capacity(1 << 22, f);
    w.write_all(b"FQS1")?;
    w.write_all(&(num_blocks as u32).to_le_bytes())?;
    w.write_all(&(total as u64).to_le_bytes())?;
    for (n, blob) in &blocks {
        w.write_all(&(*n as u32).to_le_bytes())?;
        w.write_all(&(blob.len() as u32).to_le_bytes())?;
        w.write_all(blob)?;
    }
    w.flush()?;
    eprintln!("encode: wrote {out_path}");
    Ok(())
}

fn read_exact_vec(r: &mut impl Read, n: usize) -> Result<Vec<u8>> {
    let mut v = vec![0u8; n];
    r.read_exact(&mut v)?;
    Ok(v)
}

fn cmd_decode(qbin: &str, mode: &str, sink: &str) -> Result<()> {
    let f = std::fs::File::open(qbin).with_context(|| format!("open {qbin}"))?;
    let mut r = BufReader::with_capacity(1 << 22, f);
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != b"FQS1" {
        bail!("bad magic");
    }
    let num_blocks = u32::from_le_bytes(read_exact_vec(&mut r, 4)?.try_into().unwrap()) as usize;
    let total = u64::from_le_bytes(read_exact_vec(&mut r, 8)?.try_into().unwrap()) as usize;

    let sink_f = std::fs::File::create(sink)?;
    let mut sw = BufWriter::with_capacity(1 << 22, sink_f);

    let t = Instant::now();
    match mode {
        "streaming" => {
            // Decode one block, write it, drop it. Bounded working set.
            for _ in 0..num_blocks {
                let n = u32::from_le_bytes(read_exact_vec(&mut r, 4)?.try_into().unwrap()) as usize;
                let blen =
                    u32::from_le_bytes(read_exact_vec(&mut r, 4)?.try_into().unwrap()) as usize;
                let blob = read_exact_vec(&mut r, blen)?;
                let (concat, _lengths) = fqzcomp::decompress(&blob, n)?;
                sw.write_all(&concat)?;
                // concat + blob dropped here.
            }
        }
        "allatonce" => {
            // Mimic today's fqz decode: decode ALL blocks, hold every decoded
            // quality buffer live simultaneously, THEN write. Peak = whole stream.
            let mut all: Vec<Vec<u8>> = Vec::with_capacity(num_blocks);
            for _ in 0..num_blocks {
                let n = u32::from_le_bytes(read_exact_vec(&mut r, 4)?.try_into().unwrap()) as usize;
                let blen =
                    u32::from_le_bytes(read_exact_vec(&mut r, 4)?.try_into().unwrap()) as usize;
                let blob = read_exact_vec(&mut r, blen)?;
                let (concat, _lengths) = fqzcomp::decompress(&blob, n)?;
                all.push(concat);
            }
            for concat in &all {
                sw.write_all(concat)?;
            }
        }
        other => bail!("mode must be streaming|allatonce, got {other}"),
    }
    sw.flush()?;
    eprintln!(
        "decode[{mode}]: {num_blocks} blocks / {total} reads in {:.2}s",
        t.elapsed().as_secs_f64()
    );
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage:\n  {0} encode <out.qbin> --fastq R1 [R2]\n  {0} decode <in.qbin> <allatonce|streaming> <sink>", args[0]);
        std::process::exit(2);
    }
    match args[1].as_str() {
        "encode" => {
            let out = args.get(2).context("need out.qbin")?;
            let mut fastqs = Vec::new();
            let mut i = 3;
            while i < args.len() {
                if args[i] == "--fastq" {
                    i += 1;
                    while i < args.len() && !args[i].starts_with("--") {
                        fastqs.push(args[i].clone());
                        i += 1;
                    }
                } else {
                    i += 1;
                }
            }
            if fastqs.is_empty() {
                bail!("need --fastq");
            }
            cmd_encode(out, &fastqs)
        }
        "decode" => {
            let qbin = args.get(2).context("need in.qbin")?;
            let mode = args.get(3).context("need mode")?;
            let sink = args.get(4).context("need sink path")?;
            cmd_decode(qbin, mode, sink)
        }
        other => bail!("unknown subcommand {other}"),
    }
}
