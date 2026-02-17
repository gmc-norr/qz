/// Benchmark context-aware quality compression vs BSC.
///
/// Usage: bench_quality_ctx <fastq_file> [max_reads]
use std::env;
use std::time::Instant;

use qz_lib::compression::quality_ctx;
use qz_lib::compression::bsc;
use qz_lib::io::fastq::FastqReader;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: bench_quality_ctx <fastq_file> [max_reads]");
        std::process::exit(1);
    }
    let path = &args[1];
    let max_reads: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(usize::MAX);

    // Load reads
    println!("Loading reads from {}...", path);
    let t0 = Instant::now();
    let mut reader = FastqReader::from_path(path, false).expect("Failed to open FASTQ");
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut qualities: Vec<Vec<u8>> = Vec::new();

    while let Some(rec) = reader.next().expect("read error") {
        sequences.push(rec.sequence);
        qualities.push(rec.quality.unwrap_or_default());
        if sequences.len() >= max_reads {
            break;
        }
    }
    let n_reads = sequences.len();
    let read_len = sequences[0].len();
    let total_bases = n_reads * read_len;
    println!(
        "Loaded {} reads x {} = {} bases in {:.1}s\n",
        n_reads,
        read_len,
        total_bases,
        t0.elapsed().as_secs_f64()
    );

    let seq_refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_slice()).collect();
    let qual_refs: Vec<&[u8]> = qualities.iter().map(|s| s.as_slice()).collect();

    // ================================================================
    // BSC baseline
    // ================================================================
    println!("=== BSC Compression ===");
    let qual_concat: Vec<u8> = qualities.iter().flat_map(|q| q.iter().copied()).collect();
    let t0 = Instant::now();
    let bsc_compressed = bsc::compress_parallel_adaptive(&qual_concat).expect("BSC compress failed");
    let bsc_time = t0.elapsed().as_secs_f64();
    let bsc_size = bsc_compressed.len();
    println!(
        "  BSC: {} bytes ({:.4} bits/base), {:.2}s",
        bsc_size,
        8.0 * bsc_size as f64 / total_bases as f64,
        bsc_time
    );

    // ================================================================
    // Context-aware compression (pos_bin, prev_q, delta_stable, base_pair)
    // ================================================================
    println!("\n=== Context-Aware Quality (pos_bin, prev_q, delta_stable, base_pair) ===");
    let t0 = Instant::now();
    let ctx_compressed =
        quality_ctx::compress_qualities_ctx(&qual_refs, &seq_refs).expect("ctx compress failed");
    let ctx_time = t0.elapsed().as_secs_f64();
    let ctx_size = ctx_compressed.len();
    println!(
        "  CTX: {} bytes ({:.4} bits/base), {:.2}s",
        ctx_size,
        8.0 * ctx_size as f64 / total_bases as f64,
        ctx_time
    );

    // Roundtrip verification
    let t0 = Instant::now();
    let ctx_dec =
        quality_ctx::decompress_qualities_ctx(&ctx_compressed, &seq_refs, n_reads)
            .expect("ctx decompress failed");
    let ctx_dec_time = t0.elapsed().as_secs_f64();
    let ctx_ok = (0..n_reads).all(|i| ctx_dec[i] == qualities[i]);
    println!(
        "  Roundtrip: {} ({:.2}s)",
        if ctx_ok { "OK" } else { "FAILED" },
        ctx_dec_time
    );
    if !ctx_ok {
        let mut mismatches = 0;
        for i in 0..n_reads {
            if ctx_dec[i] != qualities[i] {
                mismatches += 1;
                if mismatches <= 5 {
                    eprintln!(
                        "  MISMATCH read {}: expected {:?}, got {:?}",
                        i,
                        &qualities[i][..20.min(qualities[i].len())],
                        &ctx_dec[i][..20.min(ctx_dec[i].len())]
                    );
                }
            }
        }
        eprintln!("  Total mismatches: {}/{}", mismatches, n_reads);
    }

    // ================================================================
    // Summary
    // ================================================================
    let diff = ctx_size as i64 - bsc_size as i64;
    let pct = 100.0 * diff as f64 / bsc_size as f64;

    println!("\n=== Summary ===");
    println!("  BSC: {:>12} bytes ({:.4} bits/base)", bsc_size, 8.0 * bsc_size as f64 / total_bases as f64);
    println!("  CTX: {:>12} bytes ({:.4} bits/base)  {:+.1}% vs BSC", ctx_size, 8.0 * ctx_size as f64 / total_bases as f64, pct);
    println!(
        "  Speed CTX: compress {:.2}s, decompress {:.2}s",
        ctx_time, ctx_dec_time
    );
    println!(
        "  Speed BSC: {:.2}s",
        bsc_time
    );
}
