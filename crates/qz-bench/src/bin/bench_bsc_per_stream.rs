//! Per-stream BSC stage analysis.
//!
//! For each of the three FASTQ streams (headers, sequences, qualities), report
//! how much compression each stage of the BSC pipeline contributes:
//!
//!   1. Raw (no compression)
//!   2. LZP only          (bsc_lzp_compress)
//!   3. QLFC only         (bsc_coder_compress, static QLFC)
//!   4. QLFC only adaptive
//!   5. BWT + QLFC static (full pipeline, no LZP)
//!   6. LZP + BWT + QLFC static  (default)
//!   7. LZP + BWT + QLFC adaptive
//!   8. BWT + QLFC adaptive (no LZP)
//!   9. Block-parallel default (compress_parallel)         -- 25 MB blocks
//!  10. Block-parallel adaptive (compress_parallel_adaptive) -- 25 MB blocks
//!  11. Single big block 64 MB
//!  12. Single big block 128 MB
//!  13. Single block: entire stream (whole-file BWT context)
//!
//! Stages 2-8 and 11-13 use a *single* libbsc call on the whole stream (so
//! block-size effects are isolated from per-stage effects). Stages 9 and 10
//! reflect what qz actually does in production today.
//!
//! Usage: cargo run --release --bin bench_bsc_per_stream -- <fastq_path>

use anyhow::Result;
use libc::{c_int, c_uchar};
use qz_lib::compression::bsc as qz_bsc;
use qz_lib::io::fastq::FastqReader;
use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

// ── Direct FFI to libbsc internal modules ────────────────────────────────
//
// libbsc.h exposes only bsc_compress / bsc_decompress, but lzp/lzp.h and
// coder/coder.h declare bsc_lzp_compress / bsc_coder_compress as extern "C".
// Both translation units are linked into the static libbsc.a built by
// qz-lib/build.rs (see build.rs: lzp/lzp.cpp and coder/coder.cpp are in the
// source list), so we can call them directly.

const LIBBSC_FEATURE_FASTMODE: c_int = 1;
const LIBBSC_CODER_QLFC_STATIC: c_int = 1;
const LIBBSC_CODER_QLFC_ADAPTIVE: c_int = 2;
const LIBBSC_BLOCKSORTER_BWT: c_int = 1;

unsafe extern "C" {
    fn bsc_lzp_compress(
        input: *const c_uchar,
        output: *mut c_uchar,
        n: c_int,
        hash_size: c_int,
        min_len: c_int,
        features: c_int,
    ) -> c_int;

    fn bsc_coder_compress(
        input: *const c_uchar,
        output: *mut c_uchar,
        n: c_int,
        coder: c_int,
        features: c_int,
    ) -> c_int;

    fn bsc_compress(
        input: *const c_uchar,
        output: *mut c_uchar,
        n: c_int,
        lzp_hash_size: c_int,
        lzp_min_len: c_int,
        block_sorter: c_int,
        coder: c_int,
        features: c_int,
    ) -> c_int;
}

fn lzp_only(data: &[u8]) -> Result<usize> {
    // Default LZP params (matching libbsc.h: hash=16, minLen=128).
    let mut out = vec![0u8; data.len() + 4096];
    let n = unsafe {
        bsc_lzp_compress(
            data.as_ptr(),
            out.as_mut_ptr(),
            data.len() as c_int,
            16,
            128,
            LIBBSC_FEATURE_FASTMODE,
        )
    };
    if n < 0 {
        // Negative return values from bsc_lzp_compress are "not compressible"
        // (LIBBSC_NOT_COMPRESSIBLE = -3) or genuine errors. For our purposes
        // treat "not compressible" as a no-op and report the original size.
        return Ok(data.len());
    }
    Ok(n as usize)
}

fn coder_only(data: &[u8], coder: c_int) -> Result<usize> {
    let mut out = vec![0u8; data.len() + 4096];
    let n = unsafe {
        bsc_coder_compress(
            data.as_ptr(),
            out.as_mut_ptr(),
            data.len() as c_int,
            coder,
            LIBBSC_FEATURE_FASTMODE,
        )
    };
    if n < 0 {
        anyhow::bail!("bsc_coder_compress failed: {n}");
    }
    Ok(n as usize)
}

/// Single-block full pipeline call (parameterised lzp / coder).
fn full_pipeline(data: &[u8], lzp_hash: i32, lzp_min: i32, coder: c_int) -> Result<usize> {
    let mut out = vec![0u8; data.len() + 4096];
    let n = unsafe {
        bsc_compress(
            data.as_ptr(),
            out.as_mut_ptr(),
            data.len() as c_int,
            lzp_hash as c_int,
            lzp_min as c_int,
            LIBBSC_BLOCKSORTER_BWT,
            coder,
            LIBBSC_FEATURE_FASTMODE,
        )
    };
    if n < 0 {
        anyhow::bail!("bsc_compress failed: {n}");
    }
    Ok(n as usize)
}

/// Chunked single-block call: split into chunks of `chunk_bytes`, call
/// `bsc_compress` on each sequentially, sum the compressed sizes. This is
/// how qz currently splits big streams (25 MB) — we vary the size here to
/// measure BWT context length effect.
fn full_pipeline_chunked(
    data: &[u8],
    chunk_bytes: usize,
    lzp_hash: i32,
    lzp_min: i32,
    coder: c_int,
) -> Result<usize> {
    use rayon::prelude::*;
    let chunks: Vec<&[u8]> = data.chunks(chunk_bytes).collect();
    let sizes: Vec<usize> = chunks
        .par_iter()
        .map(|c| full_pipeline(c, lzp_hash, lzp_min, coder))
        .collect::<Result<Vec<_>>>()?;
    Ok(sizes.into_iter().sum())
}

struct Streams {
    headers: Vec<u8>,
    sequences: Vec<u8>,
    qualities: Vec<u8>,
    num_reads: u64,
    seq_chars: u64,
    qual_chars: u64,
}

fn read_streams(path: &str) -> Result<Streams> {
    let mut reader = FastqReader::from_path(path, false)?;
    let mut headers = Vec::new();
    let mut sequences = Vec::new();
    let mut qualities = Vec::new();
    let mut n = 0u64;
    let mut seq_chars = 0u64;
    let mut qual_chars = 0u64;
    while let Some(rec) = reader.next()? {
        // We concatenate raw stream contents with newline separators to mimic
        // the order/structure libbsc would see if the FASTQ subfields were
        // written linearly. The newline adds one byte per record but matches
        // a per-stream "newline-delimited records" view, which is what qz
        // serialises before BSC in some paths.
        headers.extend_from_slice(&rec.id);
        headers.push(b'\n');
        sequences.extend_from_slice(&rec.sequence);
        sequences.push(b'\n');
        if let Some(q) = rec.quality {
            qualities.extend_from_slice(&q);
            qualities.push(b'\n');
            qual_chars += q.len() as u64;
        }
        seq_chars += rec.sequence.len() as u64;
        n += 1;
    }
    Ok(Streams {
        headers,
        sequences,
        qualities,
        num_reads: n,
        seq_chars,
        qual_chars,
    })
}

fn fmt_mb(b: usize) -> String {
    format!("{:.2} MB", b as f64 / (1024.0 * 1024.0))
}

fn fmt_bits_per_symbol(out: usize, symbols: u64) -> String {
    if symbols == 0 {
        return "—".to_string();
    }
    format!("{:.3}", (out as f64 * 8.0) / symbols as f64)
}

struct Row {
    name: &'static str,
    bytes: usize,
    ms: f64,
}

fn run_row<F: FnOnce() -> Result<usize>>(name: &'static str, f: F) -> Row {
    let t = Instant::now();
    let bytes = f().unwrap_or_else(|e| {
        eprintln!("  [{name}] failed: {e}");
        0
    });
    let ms = t.elapsed().as_secs_f64() * 1000.0;
    // Per-row progress to stderr so long runs (sequences, qualities) show life.
    eprintln!("  . {name}: {bytes} B in {ms:.1}ms");
    Row { name, bytes, ms }
}

fn print_table(label: &str, raw_bytes: usize, symbol_count: u64, rows: &[Row]) {
    println!();
    println!("================ {label} ================");
    println!(
        "Raw stream: {} ({} bytes), symbols (non-newline chars where applicable): {}",
        fmt_mb(raw_bytes),
        raw_bytes,
        symbol_count,
    );
    println!();
    println!(
        "{:<48} {:>14} {:>10} {:>10} {:>14}",
        "Stage", "Compressed B", "Ratio", "Time (ms)", "Bits/symbol",
    );
    println!("{}", "-".repeat(100));
    for r in rows {
        let ratio = if r.bytes > 0 {
            raw_bytes as f64 / r.bytes as f64
        } else {
            0.0
        };
        println!(
            "{:<48} {:>14} {:>9.3}x {:>10.1} {:>14}",
            r.name,
            r.bytes,
            ratio,
            r.ms,
            fmt_bits_per_symbol(r.bytes, symbol_count),
        );
    }
}

fn bench_stream(label: &str, data: &[u8], symbol_count: u64) {
    if data.is_empty() {
        println!("\n[{label}] empty stream, skipping");
        return;
    }
    let raw_len = data.len();

    let mut rows = Vec::new();

    rows.push(Row {
        name: "01. Raw (uncompressed)",
        bytes: raw_len,
        ms: 0.0,
    });

    rows.push(run_row("02. LZP only (hash=16, minLen=128)", || {
        lzp_only(data)
    }));
    rows.push(run_row("03. QLFC static only (no BWT, no LZP)", || {
        coder_only(data, LIBBSC_CODER_QLFC_STATIC)
    }));
    rows.push(run_row("04. QLFC adaptive only (no BWT, no LZP)", || {
        coder_only(data, LIBBSC_CODER_QLFC_ADAPTIVE)
    }));
    rows.push(run_row("05. BWT + QLFC static  (no LZP)", || {
        full_pipeline(data, 0, 0, LIBBSC_CODER_QLFC_STATIC)
    }));
    rows.push(run_row("06. LZP + BWT + QLFC static  (default)", || {
        full_pipeline(data, 16, 128, LIBBSC_CODER_QLFC_STATIC)
    }));
    rows.push(run_row("07. LZP + BWT + QLFC adaptive", || {
        full_pipeline(data, 16, 128, LIBBSC_CODER_QLFC_ADAPTIVE)
    }));
    rows.push(run_row("08. BWT + QLFC adaptive (no LZP)", || {
        full_pipeline(data, 0, 0, LIBBSC_CODER_QLFC_ADAPTIVE)
    }));

    rows.push(run_row("09. qz compress_parallel (25 MB blocks)", || {
        Ok(qz_bsc::compress_parallel(data)?.len())
    }));
    rows.push(run_row(
        "10. qz compress_parallel_adaptive (25 MB blocks)",
        || Ok(qz_bsc::compress_parallel_adaptive(data)?.len()),
    ));

    rows.push(run_row(
        "11. Single-block BWT+QLFC adapt 64 MB chunks",
        || full_pipeline_chunked(data, 64 * 1024 * 1024, 16, 128, LIBBSC_CODER_QLFC_ADAPTIVE),
    ));
    rows.push(run_row(
        "12. Single-block BWT+QLFC adapt 128 MB chunks",
        || full_pipeline_chunked(data, 128 * 1024 * 1024, 16, 128, LIBBSC_CODER_QLFC_ADAPTIVE),
    ));

    // libbsc max single-block input is 1 GiB. Skip whole-file if larger.
    if data.len() <= 1_000_000_000 {
        rows.push(run_row(
            "13. Single-block BWT+QLFC adapt whole stream",
            || full_pipeline(data, 16, 128, LIBBSC_CODER_QLFC_ADAPTIVE),
        ));
    }

    print_table(label, raw_len, symbol_count, &rows);
}

fn main() -> Result<()> {
    // Force libbsc initialisation BEFORE we make any direct FFI calls; qz-lib
    // gates init behind ensure_initialized() inside its own helpers, so the
    // first time we touch them later would be too late.
    let _ = qz_bsc::compress(b"prime").unwrap();

    let fastq = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/ERR3239334_1.1m.fastq".to_string());
    let path = std::path::Path::new(&fastq);
    if !path.exists() {
        anyhow::bail!("input file not found: {fastq}");
    }
    let file_bytes = std::fs::metadata(path)?.len();
    eprintln!(
        "Reading {fastq} ({:.1} MB on disk) ...",
        file_bytes as f64 / (1024.0 * 1024.0)
    );

    let t = Instant::now();
    let streams = read_streams(&fastq)?;
    eprintln!(
        "Parsed {} reads in {:.2}s — headers={} bytes, sequences={} bytes ({} bp), qualities={} bytes ({} chars)",
        streams.num_reads,
        t.elapsed().as_secs_f64(),
        streams.headers.len(),
        streams.sequences.len(),
        streams.seq_chars,
        streams.qualities.len(),
        streams.qual_chars,
    );

    // Touch BufReader type to silence unused-import warning if file ever opened directly.
    let _ = std::any::type_name::<BufReader<File>>();

    bench_stream("HEADERS", &streams.headers, streams.headers.len() as u64);
    bench_stream("SEQUENCES", &streams.sequences, streams.seq_chars);
    bench_stream("QUALITIES", &streams.qualities, streams.qual_chars);

    println!("\nDone.");
    Ok(())
}
