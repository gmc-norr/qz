/// Quick analysis of syncmer hash frequency distribution in FASTQ data
use std::collections::HashMap;
use std::env;
use qz_lib::io::fastq::FastqReader;
use qz_lib::compression::dna_utils::compute_min_syncmer_hash;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: analyze_hashes <fastq_file> [max_reads]");
        std::process::exit(1);
    }
    let path = &args[1];
    let max_reads: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(usize::MAX);

    let mut reader = FastqReader::from_path(path, false).expect("Failed to open FASTQ");
    let mut hash_counts: HashMap<u64, u32> = HashMap::new();
    let mut total = 0usize;
    let mut no_hash = 0usize;

    while let Some(rec) = reader.next().expect("read error") {
        if total >= max_reads { break; }
        total += 1;
        let h = compute_min_syncmer_hash(&rec.sequence);
        if h == u64::MAX {
            no_hash += 1;
        } else {
            *hash_counts.entry(h).or_insert(0) += 1;
        }
    }

    let num_unique = hash_counts.len();
    let mut freq_dist: HashMap<u32, u32> = HashMap::new(); // count -> how many hashes have that count
    for &c in hash_counts.values() {
        *freq_dist.entry(c).or_insert(0) += 1;
    }

    println!("=== Syncmer Hash Analysis (k=21, s=10 open syncmer) ===");
    println!("Total reads: {}", total);
    println!("Reads with no hash (too short / all N): {}", no_hash);
    println!("Unique hashes: {}", num_unique);
    println!("Hash-to-read ratio: {:.4}", num_unique as f64 / (total - no_hash) as f64);
    println!();

    // Frequency distribution
    let mut freq_vec: Vec<(u32, u32)> = freq_dist.into_iter().collect();
    freq_vec.sort();

    println!("Frequency distribution (count -> num_hashes -> total_reads_covered):");
    let mut cumulative_reads = 0u64;
    for &(count, num_hashes) in &freq_vec {
        let reads_covered = count as u64 * num_hashes as u64;
        cumulative_reads += reads_covered;
        if count <= 10 || count % 100 == 0 || num_hashes > 100 {
            println!("  count={:<6} hashes={:<8} reads={:<10} cum_reads={:<10} ({:.1}%)",
                count, num_hashes, reads_covered, cumulative_reads,
                100.0 * cumulative_reads as f64 / total as f64);
        }
    }
    println!();

    // Top 20 hashes
    let mut top: Vec<(u64, u32)> = hash_counts.into_iter().collect();
    top.sort_by(|a, b| b.1.cmp(&a.1));
    println!("Top 20 most frequent hashes:");
    for (hash, count) in top.iter().take(20) {
        println!("  hash=0x{:016x} count={}", hash, count);
    }

    // Coverage by top N patterns
    println!();
    println!("Cumulative coverage by top N patterns:");
    let mut cum = 0u64;
    for (i, (_hash, count)) in top.iter().enumerate() {
        cum += *count as u64;
        if i < 20 || (i + 1) % 50 == 0 || i + 1 == 100 || i + 1 == 200 || i + 1 == 500 || i + 1 == 1000 {
            println!("  top {:>5}: {:<8} reads ({:.2}%)", i + 1, cum, 100.0 * cum as f64 / total as f64);
        }
    }
}
