/// Benchmark: Template-hybrid at different scales
/// Tests baseline and template-hybrid on the input FASTQ

use std::io::BufRead;
use std::time::Instant;

use qz_lib::compression::template::{
    compress_sequences_template_hybrid,
    TemplateParams,
};

fn main() {
    let fastq_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "real_data/ERR3239334_1.500k.fastq".to_string());

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
    let first_len = sequences.first().map(|s| s.len()).unwrap_or(0);

    eprintln!(
        "Loaded {} reads, {} bases ({:.1} MB), read_len={}\n",
        num_reads,
        total_bases,
        total_bases as f64 / (1024.0 * 1024.0),
        first_len,
    );

    println!(
        "{:<50} {:>12} {:>8} {:>12}",
        "Method", "Compressed", "Ratio", "Time"
    );
    println!("{}", "-".repeat(85));

    // === 1. Baseline: Raw ASCII + BSC ===
    {
        let raw_bytes: Vec<u8> = sequences
            .iter()
            .flat_map(|s| s.as_bytes())
            .copied()
            .collect();
        let t = Instant::now();
        let compressed = qz_lib::compression::bsc::compress_parallel_adaptive(&raw_bytes).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_bases as f64 / compressed.len() as f64;
        println!(
            "{:<50} {:>10} B {:>7.2}x {:>10.1}ms",
            "1. Baseline (raw + BSC)",
            compressed.len(),
            ratio,
            elapsed.as_secs_f64() * 1000.0,
        );
    }

    // === 2. Template-hybrid fast k=31 syncmer ===
    {
        eprintln!("\n--- Template-hybrid fast k=31 syncmer ---");
        let params = TemplateParams {
            k: 31,
            fast: true,
            use_syncmers: true,
            ..TemplateParams::default()
        };
        let t = Instant::now();
        let compressed = compress_sequences_template_hybrid(&sequences, &params).unwrap();
        let elapsed = t.elapsed();
        let ratio = total_bases as f64 / compressed.len() as f64;
        println!(
            "{:<50} {:>10} B {:>7.2}x {:>10.1}ms",
            "2. Template-hybrid fast k=31 syncmer",
            compressed.len(),
            ratio,
            elapsed.as_secs_f64() * 1000.0,
        );
    }

    println!("{}", "-".repeat(85));
    println!(
        "Input: {} reads x {} bp = {} bytes ({:.1} MB)",
        num_reads,
        first_len,
        total_bases,
        total_bases as f64 / (1024.0 * 1024.0)
    );
}
