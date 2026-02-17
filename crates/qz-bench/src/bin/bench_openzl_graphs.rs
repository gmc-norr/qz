/// Benchmark: compare OpenZL graph-based compression strategies on DNA sequences
///
/// Extracts raw DNA sequences from a FASTQ file and compresses with various
/// OpenZL graph strategies and BSC (baseline).

use std::io::BufRead;
use std::time::Instant;

fn main() {
    let fastq_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/NA12878_exome_50k.fastq".to_string());

    eprintln!("Reading FASTQ: {}", fastq_path);
    let file = std::fs::File::open(&fastq_path).expect("Cannot open FASTQ file");
    let reader = std::io::BufReader::new(file);

    let mut sequences = Vec::new();
    let mut qualities = Vec::new();
    let mut line_num = 0u64;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        match line_num % 4 {
            1 => sequences.extend_from_slice(line.trim_end().as_bytes()),
            3 => qualities.extend_from_slice(line.trim_end().as_bytes()),
            _ => {}
        }
        line_num += 1;
    }

    let num_reads = line_num / 4;
    eprintln!(
        "Loaded {} reads, {} bytes sequences, {} bytes qualities",
        num_reads,
        sequences.len(),
        qualities.len()
    );

    println!("=== OpenZL Graph Compression Benchmark ===");
    println!(
        "Input: {} reads, {} bytes ({:.1} MB)",
        num_reads,
        sequences.len(),
        sequences.len() as f64 / (1024.0 * 1024.0)
    );
    println!();

    benchmark_stream("SEQUENCES", &sequences);
    println!();
    benchmark_stream("QUALITIES", &qualities);
}

fn try_compress(
    name: &str,
    data: &[u8],
    compress_fn: impl FnOnce(&[u8]) -> anyhow::Result<Vec<u8>>,
) {
    let t = Instant::now();
    let result = compress_fn(data);
    let comp_time = t.elapsed();

    match result {
        Ok(compressed) => {
            let t2 = Instant::now();
            match qz_lib::compression::openzl::decompress(&compressed) {
                Ok(decompressed) => {
                    let roundtrip = t2.elapsed();
                    assert_eq!(data, decompressed.as_slice(), "{} roundtrip mismatch", name);
                    print_row(name, data.len(), compressed.len(), comp_time, roundtrip);
                }
                Err(e) => {
                    println!("{:<30} compressed OK ({} B) but DECOMPRESS FAILED: {}", name, compressed.len(), e);
                }
            }
        }
        Err(e) => {
            println!("{:<30} FAILED: {}", name, e);
        }
    }
}

fn benchmark_stream(label: &str, data: &[u8]) {
    println!(
        "--- {} ({} bytes, {:.1} MB) ---",
        label,
        data.len(),
        data.len() as f64 / (1024.0 * 1024.0)
    );
    println!(
        "{:<30} {:>12} {:>10} {:>12} {:>10}",
        "Method", "Compressed", "Ratio", "Comp Time", "Decomp"
    );
    println!("{}", "-".repeat(76));

    // 1. BSC adaptive (baseline)
    {
        let t = Instant::now();
        let compressed = qz_lib::compression::bsc::compress_parallel_adaptive(data)
            .expect("BSC compress failed");
        let comp_time = t.elapsed();

        let t2 = Instant::now();
        let decompressed = qz_lib::compression::bsc::decompress_parallel(&compressed)
            .expect("BSC decompress failed");
        let roundtrip = t2.elapsed();

        assert_eq!(data, decompressed.as_slice(), "BSC roundtrip mismatch");
        print_row("BSC (adaptive)", data.len(), compressed.len(), comp_time, roundtrip);
    }

    // 2. OpenZL generic (default = compress_generic internally)
    try_compress("OpenZL generic", data, |d| {
        qz_lib::compression::openzl::compress(d)
    });

    // 3. OpenZL ACE (Adaptive Codec Engine)
    try_compress("OpenZL ACE", data, |d| {
        qz_lib::compression::openzl::compress_ace(d)
    });

    // 4. OpenZL DNA numeric (serial -> u8 numeric + entropy)
    try_compress("OpenZL DNA numeric", data, |d| {
        qz_lib::compression::openzl::compress_dna_numeric(d)
    });

    // 5. OpenZL delta + entropy (serial -> u8 -> delta -> entropy)
    try_compress("OpenZL delta+entropy", data, |d| {
        qz_lib::compression::openzl::compress_delta_entropy(d)
    });

    // 6. OpenZL FSE (Finite State Entropy)
    try_compress("OpenZL FSE", data, |d| {
        qz_lib::compression::openzl::compress_fse(d)
    });

    // 7. OpenZL Huffman
    try_compress("OpenZL Huffman", data, |d| {
        qz_lib::compression::openzl::compress_huffman(d)
    });

    // 8. OpenZL transpose + zstd
    try_compress("OpenZL transpose+zstd", data, |d| {
        qz_lib::compression::openzl::compress_transpose_zstd(d)
    });

    // 9. OpenZL field LZ
    try_compress("OpenZL field_lz", data, |d| {
        qz_lib::compression::openzl::compress_field_lz(d)
    });

    // 10. OpenZL clustering
    try_compress("OpenZL clustering", data, |d| {
        qz_lib::compression::openzl::compress_clustering(d)
    });
}

fn print_row(
    name: &str,
    original: usize,
    compressed: usize,
    comp_time: std::time::Duration,
    decomp_time: std::time::Duration,
) {
    let ratio = original as f64 / compressed as f64;
    println!(
        "{:<30} {:>10} B {:>8.2}x {:>10.1}ms {:>8.1}ms",
        name,
        compressed,
        ratio,
        comp_time.as_secs_f64() * 1000.0,
        decomp_time.as_secs_f64() * 1000.0,
    );
}
