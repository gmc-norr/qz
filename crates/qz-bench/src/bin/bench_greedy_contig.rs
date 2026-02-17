/// Focused benchmark: Raw+BSC vs Greedy Contig+BSC
/// Designed to run on larger datasets (5M+ reads) where coverage matters.

use std::io::BufRead;
use std::time::Instant;

fn main() {
    let fastq_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/ERR3239334_1.5m.fastq".to_string());

    eprintln!("Reading FASTQ: {}", fastq_path);
    let file = std::fs::File::open(&fastq_path).expect("Cannot open FASTQ file");
    let reader = std::io::BufReader::new(file);

    let mut sequences: Vec<String> = Vec::new();
    let mut line_num = 0u64;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line_num % 4 == 1 {
            sequences.push(line.trim_end().to_string());
        }
        line_num += 1;
    }

    let num_reads = sequences.len();
    let total_bases: usize = sequences.iter().map(|s| s.len()).sum();
    let first_len = sequences[0].len();

    eprintln!(
        "Loaded {} reads, {} bases ({:.1} MB), read_len={}\n",
        num_reads,
        total_bases,
        total_bases as f64 / (1024.0 * 1024.0),
        first_len,
    );

    println!(
        "{:<40} {:>12} {:>8} {:>12} {:>10}",
        "Method", "Compressed", "Ratio", "Comp Time", "Decomp"
    );
    println!("{}", "-".repeat(86));

    // === 1. Raw ASCII + BSC (baseline) ===
    {
        let raw_bytes: Vec<u8> = sequences.iter().flat_map(|s| s.as_bytes()).copied().collect();
        let t = Instant::now();
        let compressed = qz_lib::compression::bsc::compress_parallel_adaptive(&raw_bytes).unwrap();
        let comp_time = t.elapsed();

        let t2 = Instant::now();
        let _decompressed = qz_lib::compression::bsc::decompress_parallel(&compressed).unwrap();
        let decomp_time = t2.elapsed();

        let ratio = total_bases as f64 / compressed.len() as f64;
        println!(
            "{:<40} {:>10} B {:>7.2}x {:>10.1}ms {:>8.1}ms",
            "1. Raw ASCII + BSC",
            compressed.len(),
            ratio,
            comp_time.as_secs_f64() * 1000.0,
            decomp_time.as_secs_f64() * 1000.0,
        );
    }

    // === 2. Greedy contig + BSC ===
    {
        let t = Instant::now();
        let compressed =
            qz_lib::compression::greedy_contig::compress_sequences_greedy(&sequences).unwrap();
        let comp_time = t.elapsed();

        // Decompress uses same DBG1 format as de Bruijn
        let t2 = Instant::now();
        let decompressed =
            qz_lib::compression::debruijn::decompress_sequences_debruijn(&compressed, num_reads)
                .unwrap();
        let decomp_time = t2.elapsed();

        // Verify roundtrip (spot check first 1000)
        let check_count = num_reads.min(1000);
        for i in 0..check_count {
            assert_eq!(
                sequences[i], decompressed[i],
                "Greedy contig roundtrip mismatch at read {}",
                i
            );
        }

        let ratio = total_bases as f64 / compressed.len() as f64;
        println!(
            "{:<40} {:>10} B {:>7.2}x {:>10.1}ms {:>8.1}ms",
            "2. Greedy contig + BSC",
            compressed.len(),
            ratio,
            comp_time.as_secs_f64() * 1000.0,
            decomp_time.as_secs_f64() * 1000.0,
        );
    }

    println!("{}", "-".repeat(86));
    println!(
        "Input: {} reads x {} bp = {} bytes ({:.1} MB)",
        num_reads,
        first_len,
        total_bases,
        total_bases as f64 / (1024.0 * 1024.0)
    );
}
