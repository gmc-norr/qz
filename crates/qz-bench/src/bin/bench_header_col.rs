/// Benchmark columnar header compression vs raw BSC.
///
/// Parses Illumina-style FASTQ headers into structured fields,
/// stores each field as a separate binary column, BSC-compresses
/// each column independently, and compares total size vs raw BSC.
///
/// Usage: bench_header_col <fastq_file> [max_reads]
use std::env;
use std::time::Instant;

use qz_lib::compression::bsc;
use qz_lib::io::fastq::FastqReader;

use qz_lib::compression::header_col;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: bench_header_col <fastq_file> [max_reads]");
        std::process::exit(1);
    }
    let path = &args[1];
    let max_reads: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(usize::MAX);

    // ================================================================
    // Load headers
    // ================================================================
    println!("Loading headers from {}...", path);
    let t0 = Instant::now();
    let mut reader = FastqReader::from_path(path, false).expect("Failed to open FASTQ");
    let mut headers: Vec<String> = Vec::new();

    while let Some(rec) = reader.next().expect("read error") {
        headers.push(String::from_utf8(rec.id).expect("non-UTF8 header"));
        if headers.len() >= max_reads {
            break;
        }
    }
    let n_reads = headers.len();
    let raw_size: usize = headers.iter().map(|h| h.len()).sum();
    println!(
        "Loaded {} headers ({} bytes raw, avg {:.0} bytes/header) in {:.1}s\n",
        n_reads,
        raw_size,
        raw_size as f64 / n_reads as f64,
        t0.elapsed().as_secs_f64()
    );

    // ================================================================
    // BSC baseline (raw ASCII, same as current QZ approach)
    // ================================================================
    println!("=== BSC Baseline (raw ASCII) ===");
    let t0 = Instant::now();
    // Reproduce what compress_headers_bsc_with does: varint-prefixed concat
    let mut header_stream = Vec::new();
    for h in &headers {
        let len = h.len();
        let mut v = len;
        while v >= 0x80 {
            header_stream.push(((v & 0x7F) | 0x80) as u8);
            v >>= 7;
        }
        header_stream.push(v as u8);
        header_stream.extend_from_slice(h.as_bytes());
    }
    let bsc_compressed = bsc::compress_parallel_adaptive(&header_stream).expect("BSC failed");
    let bsc_time = t0.elapsed().as_secs_f64();
    let bsc_size = bsc_compressed.len();
    println!(
        "  BSC: {} bytes ({:.2}x from {} raw), {:.2}s\n",
        bsc_size,
        raw_size as f64 / bsc_size as f64,
        raw_size,
        bsc_time
    );

    // ================================================================
    // Columnar compression
    // ================================================================
    println!("=== Columnar Compression ===");
    let header_refs: Vec<&str> = headers.iter().map(|h| h.as_str()).collect();

    let t0 = Instant::now();
    let compressed = header_col::compress_headers_columnar(&header_refs)
        .expect("columnar compress failed");
    let col_time = t0.elapsed().as_secs_f64();
    let col_size = compressed.len();
    println!(
        "  Columnar: {} bytes ({:.2}x from {} raw), {:.2}s",
        col_size,
        raw_size as f64 / col_size as f64,
        raw_size,
        col_time
    );

    // ================================================================
    // Roundtrip verification
    // ================================================================
    println!("\n=== Roundtrip Verification ===");
    let t0 = Instant::now();
    let decompressed = header_col::decompress_headers_columnar(&compressed, n_reads)
        .expect("columnar decompress failed");
    let dec_time = t0.elapsed().as_secs_f64();

    let mut mismatches = 0;
    for i in 0..n_reads {
        if decompressed[i] != headers[i].as_bytes() {
            mismatches += 1;
            if mismatches <= 5 {
                eprintln!(
                    "MISMATCH read {}: expected {:?}, got {:?}",
                    i, &headers[i], &decompressed[i]
                );
            }
        }
    }
    if mismatches == 0 {
        println!("  Roundtrip OK: all {} headers match ({:.2}s)", n_reads, dec_time);
    } else {
        println!("  ROUNDTRIP FAILED: {} / {} mismatches", mismatches, n_reads);
    }

    // ================================================================
    // Summary
    // ================================================================
    let diff = col_size as i64 - bsc_size as i64;
    let pct = 100.0 * diff as f64 / bsc_size as f64;
    println!("\n=== Summary ===");
    println!("  Raw ASCII:   {:>10} bytes", raw_size);
    println!("  BSC:         {:>10} bytes ({:.2}x)", bsc_size, raw_size as f64 / bsc_size as f64);
    println!("  Columnar:    {:>10} bytes ({:.2}x)", col_size, raw_size as f64 / col_size as f64);
    println!(
        "  Delta:       {:>+10} bytes ({:+.1}%)",
        diff, pct
    );
    println!(
        "  Speed: compress {:.2}s, decompress {:.2}s (BSC: {:.2}s)",
        col_time, dec_time, bsc_time
    );
}
