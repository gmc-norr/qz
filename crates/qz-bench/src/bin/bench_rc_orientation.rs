/// Benchmark: reverse-complement orientation strategies for sequence compression.
///
/// This keeps record order fixed and asks whether choosing a stored orientation
/// for each read improves the BSC-compressed sequence stream. Orientation
/// strategies that need restoration metadata include the compressed flag cost in
/// their net size.
///
/// Usage:
///   cargo run --release -p qz-bench --bin bench_rc_orientation -- [FASTQ] [max_reads]
use std::io::BufRead;
use std::time::{Duration, Instant};

use anyhow::{Context, Result};
use qz_lib::compression::{bsc, dna_utils};
use rustc_hash::FxHashMap;

const K: usize = 21;
const S: usize = 10;
const DELTA_IDENTITY_THRESHOLD: f64 = 0.60;

fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let fastq_path = args
        .next()
        .unwrap_or_else(|| "real_data/ERR3239334_1.500k.fastq".to_string());
    let max_reads = args
        .next()
        .map(|s| s.parse::<usize>())
        .transpose()
        .context("max_reads must be an integer")?
        .unwrap_or(usize::MAX);

    eprintln!("Reading FASTQ: {fastq_path}");
    let sequences = read_fastq_sequences(&fastq_path, max_reads)?;
    anyhow::ensure!(!sequences.is_empty(), "no sequences read from {fastq_path}");

    let num_reads = sequences.len();
    let total_bases: usize = sequences.iter().map(Vec::len).sum();
    let first_len = sequences[0].len();
    let const_len = sequences
        .iter()
        .all(|s| s.len() == first_len)
        .then_some(first_len);

    eprintln!(
        "Loaded {} reads, {} bases ({:.1} MiB), read_len={}{}\n",
        num_reads,
        total_bases,
        total_bases as f64 / 1_048_576.0,
        first_len,
        if const_len.is_some() {
            " (uniform)"
        } else {
            " (variable)"
        },
    );

    println!(
        "{:<34} {:>12} {:>15} {:>12} {:>8} {:>9} {:>11}  Notes",
        "Strategy", "Seq BSC", "Meta byte/bit", "Net byte", "Ratio", "vs base", "Time"
    );
    println!("{}", "-".repeat(128));

    let baseline = bench_raw_original(&sequences, const_len)?;

    let start = Instant::now();
    let (oriented, flags, stats) = orient_lexicographic(&sequences);
    bench_oriented(
        "lexicographic canonical",
        &oriented,
        &flags,
        const_len,
        total_bases,
        baseline,
        start.elapsed(),
        &stats,
    )?;

    let start = Instant::now();
    let (oriented, flags, stats) = orient_minimizer_strand(&sequences);
    bench_oriented(
        "minimizer-strand canonical",
        &oriented,
        &flags,
        const_len,
        total_bases,
        baseline,
        start.elapsed(),
        &stats,
    )?;

    let start = Instant::now();
    let (oriented, flags, stats) = orient_bucket_representative(&sequences);
    bench_oriented(
        "bucket representative",
        &oriented,
        &flags,
        const_len,
        total_bases,
        baseline,
        start.elapsed(),
        &stats,
    )?;

    if num_reads <= 2_000_000 {
        let start = Instant::now();
        let (oriented, flags, stats) = orient_adjacent_dp(&sequences);
        bench_oriented(
            "adjacent DP optimum",
            &oriented,
            &flags,
            const_len,
            total_bases,
            baseline,
            start.elapsed(),
            &stats,
        )?;
    } else {
        println!(
            "{:<34} {:>12} {:>15} {:>12} {:>8} {:>9} {:>11}  skipped >2M reads",
            "adjacent DP optimum", "-", "-", "-", "-", "-", "-"
        );
    }

    let start = Instant::now();
    let (stream, stats) = build_delta_stream(&sequences, false);
    bench_delta_stream(
        "inline delta (current)",
        &stream,
        total_bases,
        num_reads,
        baseline,
        start.elapsed(),
        &stats,
    )?;

    let start = Instant::now();
    let (stream, stats) = build_delta_stream(&sequences, true);
    bench_delta_stream(
        "inline delta + RC-aware",
        &stream,
        total_bases,
        num_reads,
        baseline,
        start.elapsed(),
        &stats,
    )?;

    println!("\n=== CHUNK-LOCAL CLUSTERING / REORDER STRATEGIES ===\n");

    let start = Instant::now();
    let plan = cluster_by_key(&sequences, |seq| {
        let (hash, _) = min_syncmer_hash_and_strand(seq).unwrap_or((u64::MAX, false));
        (hash.to_be_bytes().to_vec(), false)
    });
    bench_clustered(
        "cluster: sym-minimizer",
        &sequences,
        &plan,
        const_len,
        total_bases,
        baseline,
        start.elapsed(),
    )?;

    let start = Instant::now();
    let plan = cluster_by_key(&sequences, |seq| {
        let (hash, use_rc) = min_syncmer_hash_and_strand(seq).unwrap_or((u64::MAX, false));
        (hash.to_be_bytes().to_vec(), use_rc)
    });
    bench_clustered(
        "cluster: sym-minimizer + RC",
        &sequences,
        &plan,
        const_len,
        total_bases,
        baseline,
        start.elapsed(),
    )?;

    let start = Instant::now();
    let plan = cluster_by_key(&sequences, |seq| {
        let (hashes, use_rc) = top_syncmer_hashes_and_strand(seq, 2);
        let mut key = Vec::with_capacity(16);
        for hash in hashes {
            key.extend_from_slice(&hash.to_be_bytes());
        }
        (key, use_rc)
    });
    bench_clustered(
        "cluster: top2 syncmers + RC",
        &sequences,
        &plan,
        const_len,
        total_bases,
        baseline,
        start.elapsed(),
    )?;

    let start = Instant::now();
    let plan = cluster_by_key(&sequences, |seq| {
        let (hashes, use_rc) = top_syncmer_hashes_and_strand(seq, 2);
        let mut key = Vec::with_capacity(32);
        for hash in hashes {
            key.extend_from_slice(&hash.to_be_bytes());
        }
        key.extend_from_slice(&position_sketch(seq, use_rc, 16));
        (key, use_rc)
    });
    bench_clustered(
        "cluster: top2 + position sketch",
        &sequences,
        &plan,
        const_len,
        total_bases,
        baseline,
        start.elapsed(),
    )?;

    let start = Instant::now();
    let plan = cluster_top2_bucket_representative(&sequences);
    bench_clustered(
        "cluster: top2 bucket-rep RC",
        &sequences,
        &plan,
        const_len,
        total_bases,
        baseline,
        start.elapsed(),
    )?;

    println!(
        "
=== BASE SORT / SEMI-SORT ORDER-CODING STRATEGIES ===
"
    );
    for window in [2usize, 4, 8, 16, 32, 75, 150] {
        let label = format!("semi-sort bases w={window}");
        bench_semisorted_window(&label, &sequences, const_len, window, total_bases, baseline)?;
    }

    Ok(())
}

#[derive(Clone, Copy)]
struct Baseline {
    compressed_len: usize,
}

#[derive(Default)]
struct OrientationStats {
    rc_reads: usize,
    comparable_reads: usize,
    rc_wins: usize,
    delta_reads: usize,
}

struct ClusterPlan {
    perm: Vec<u32>,
    flags: Vec<u8>,
    stats: OrientationStats,
}

fn bench_clustered(
    label: &str,
    sequences: &[Vec<u8>],
    plan: &ClusterPlan,
    const_len: Option<usize>,
    total_bases: usize,
    baseline: Baseline,
    transform_time: Duration,
) -> Result<()> {
    let t = Instant::now();
    let stream = build_clustered_sequence_stream(sequences, &plan.perm, &plan.flags, const_len);
    let packed_flags = pack_bits(&plan.flags);
    let build_time = t.elapsed();

    let t = Instant::now();
    let seq_compressed = compress_bsc(&stream)?;
    let perm_cost = compress_permutation_best(&plan.perm)?;
    let byte_flags = if plan.stats.rc_reads == 0 {
        0
    } else {
        compress_bsc(&plan.flags)?.len()
    };
    let bit_flags = if plan.stats.rc_reads == 0 {
        0
    } else {
        compress_bsc(&packed_flags)?.len()
    };
    let bsc_time = t.elapsed();

    let meta_byte = perm_cost + byte_flags;
    let meta_bit = perm_cost + bit_flags;
    let net_byte = seq_compressed.len() + meta_byte;
    let net_bit = seq_compressed.len() + meta_bit;
    let notes = format!(
        "perm={} rc={:.1}% cmp={} rc_wins={} bit_net={}",
        humanize(perm_cost),
        pct(plan.stats.rc_reads, sequences.len()),
        plan.stats.comparable_reads,
        plan.stats.rc_wins,
        humanize(net_bit),
    );
    print_row(
        label,
        seq_compressed.len(),
        Some((meta_byte, meta_bit)),
        net_byte,
        total_bases,
        baseline,
        transform_time + build_time + bsc_time,
        &notes,
    );
    Ok(())
}

fn cluster_by_key<F>(sequences: &[Vec<u8>], mut make_key: F) -> ClusterPlan
where
    F: FnMut(&[u8]) -> (Vec<u8>, bool),
{
    let mut items: Vec<(Vec<u8>, u32, bool)> = sequences
        .iter()
        .enumerate()
        .map(|(idx, seq)| {
            let (key, use_rc) = make_key(seq);
            (key, idx as u32, use_rc)
        })
        .collect();
    items.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

    let mut stats = OrientationStats::default();
    let mut prev_key: Option<&[u8]> = None;
    let mut repeated_bucket_reads = 0usize;
    for (key, _, use_rc) in &items {
        stats.rc_reads += *use_rc as usize;
        match prev_key {
            Some(prev) if prev == key.as_slice() => {
                repeated_bucket_reads += 1;
            }
            _ => {
                prev_key = Some(key.as_slice());
            }
        }
    }
    stats.comparable_reads = repeated_bucket_reads;

    ClusterPlan {
        perm: items.iter().map(|(_, idx, _)| *idx).collect(),
        flags: items.iter().map(|(_, _, use_rc)| *use_rc as u8).collect(),
        stats,
    }
}

fn cluster_top2_bucket_representative(sequences: &[Vec<u8>]) -> ClusterPlan {
    let mut items: Vec<(Vec<u8>, u32, bool)> = sequences
        .iter()
        .enumerate()
        .map(|(idx, seq)| {
            let (hashes, fallback_rc) = top_syncmer_hashes_and_strand(seq, 2);
            let mut key = Vec::with_capacity(16);
            for hash in hashes {
                key.extend_from_slice(&hash.to_be_bytes());
            }
            (key, idx as u32, fallback_rc)
        })
        .collect();
    items.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

    let mut stats = OrientationStats::default();
    let mut flags = Vec::with_capacity(items.len());
    let mut current_key: Option<Vec<u8>> = None;
    let mut representative: Option<Vec<u8>> = None;

    for (key, idx, fallback_rc) in &items {
        let seq = &sequences[*idx as usize];
        let mut use_rc = *fallback_rc;
        if current_key.as_deref() != Some(key.as_slice()) {
            current_key = Some(key.clone());
            representative = None;
        }

        if let Some(rep) = representative.as_ref()
            && rep.len() == seq.len()
        {
            stats.comparable_reads += 1;
            let fwd = matching_bases(seq, false, rep);
            let rc = matching_bases(seq, true, rep);
            if rc > fwd {
                use_rc = true;
                stats.rc_wins += 1;
            } else if fwd > rc {
                use_rc = false;
            }
        }

        let encoded = if use_rc {
            dna_utils::reverse_complement(seq)
        } else {
            seq.clone()
        };
        if representative.is_none() {
            representative = Some(encoded);
        }
        stats.rc_reads += use_rc as usize;
        flags.push(use_rc as u8);
    }

    ClusterPlan {
        perm: items.iter().map(|(_, idx, _)| *idx).collect(),
        flags,
        stats,
    }
}

fn build_clustered_sequence_stream(
    sequences: &[Vec<u8>],
    perm: &[u32],
    flags: &[u8],
    const_len: Option<usize>,
) -> Vec<u8> {
    let include_len = const_len.is_none();
    let mut out = Vec::with_capacity(
        sequences.iter().map(Vec::len).sum::<usize>()
            + if include_len { sequences.len() * 2 } else { 0 },
    );
    for (&idx, &flag) in perm.iter().zip(flags.iter()) {
        let seq = &sequences[idx as usize];
        if include_len {
            write_varint(&mut out, seq.len());
        }
        if flag != 0 {
            for i in 0..seq.len() {
                out.push(oriented_base(seq, i, true));
            }
        } else {
            out.extend_from_slice(seq);
        }
    }
    out
}

fn compress_permutation_best(perm: &[u32]) -> Result<usize> {
    let forward = compress_permutation_delta(perm)?;
    let inverse = invert_permutation(perm);
    let inverse = compress_permutation_delta(&inverse)?;
    Ok(forward.min(inverse))
}

fn compress_permutation_delta(perm: &[u32]) -> Result<usize> {
    let mut buf = Vec::with_capacity(perm.len() * 3);
    let mut prev = 0i64;
    for &idx in perm {
        let delta = idx as i64 - prev;
        prev = idx as i64;
        let mut value = ((delta << 1) ^ (delta >> 63)) as u64;
        while value >= 0x80 {
            buf.push(((value & 0x7f) as u8) | 0x80);
            value >>= 7;
        }
        buf.push(value as u8);
    }
    Ok(compress_bsc(&buf)?.len())
}

fn invert_permutation(perm: &[u32]) -> Vec<u32> {
    let mut inverse = vec![0u32; perm.len()];
    for (sorted_idx, &original_idx) in perm.iter().enumerate() {
        inverse[original_idx as usize] = sorted_idx as u32;
    }
    inverse
}

fn top_syncmer_hashes_and_strand(seq: &[u8], limit: usize) -> (Vec<u64>, bool) {
    let fwd = top_syncmer_hashes_one_orientation(seq, limit);
    let rc_seq = dna_utils::reverse_complement(seq);
    let rc = top_syncmer_hashes_one_orientation(&rc_seq, limit);
    if rc < fwd { (rc, true) } else { (fwd, false) }
}

fn top_syncmer_hashes_one_orientation(seq: &[u8], limit: usize) -> Vec<u64> {
    let mut hashes = Vec::new();
    if seq.len() >= K {
        for pos in syncmers::find_syncmers_pos(K, S, &[0], seq) {
            if pos + K > seq.len() {
                continue;
            }
            if let Some(hash) = dna_utils::kmer_to_hash(&seq[pos..pos + K]) {
                hashes.push(hash);
            }
        }
    }
    hashes.sort_unstable();
    hashes.dedup();
    hashes.truncate(limit);
    while hashes.len() < limit {
        hashes.push(u64::MAX);
    }
    hashes
}

fn position_sketch(seq: &[u8], rc: bool, samples: usize) -> Vec<u8> {
    if seq.is_empty() || samples == 0 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(samples);
    if samples == 1 {
        out.push(oriented_base(seq, 0, rc));
        return out;
    }
    for i in 0..samples {
        let pos = i * (seq.len() - 1) / (samples - 1);
        out.push(oriented_base(seq, pos, rc));
    }
    out
}

fn read_fastq_sequences(path: &str, max_reads: usize) -> Result<Vec<Vec<u8>>> {
    let file = std::fs::File::open(path).with_context(|| format!("open {path}"))?;
    let reader = std::io::BufReader::new(file);
    let mut sequences = Vec::new();
    for (line_idx, line) in reader.lines().enumerate() {
        let line = line?;
        if line_idx % 4 == 1 {
            sequences.push(line.trim_end().as_bytes().to_vec());
            if sequences.len() >= max_reads {
                break;
            }
        }
    }
    Ok(sequences)
}

fn bench_raw_original(sequences: &[Vec<u8>], const_len: Option<usize>) -> Result<Baseline> {
    let total_bases: usize = sequences.iter().map(Vec::len).sum();
    let t = Instant::now();
    let stream = build_sequence_stream(sequences, const_len);
    let build_time = t.elapsed();

    let t = Instant::now();
    let compressed = compress_bsc(&stream)?;
    let bsc_time = t.elapsed();

    let baseline = Baseline {
        compressed_len: compressed.len(),
    };
    print_row(
        "original orientation",
        compressed.len(),
        None,
        compressed.len(),
        total_bases,
        baseline,
        build_time + bsc_time,
        "",
    );
    Ok(baseline)
}

fn bench_oriented(
    label: &str,
    sequences: &[Vec<u8>],
    flags: &[u8],
    const_len: Option<usize>,
    total_bases: usize,
    baseline: Baseline,
    transform_time: Duration,
    stats: &OrientationStats,
) -> Result<()> {
    let t = Instant::now();
    let stream = build_sequence_stream(sequences, const_len);
    let packed_flags = pack_bits(flags);
    let build_time = t.elapsed();

    let t = Instant::now();
    let seq_compressed = compress_bsc(&stream)?;
    let byte_flags = compress_bsc(flags)?;
    let bit_flags = compress_bsc(&packed_flags)?;
    let bsc_time = t.elapsed();

    let total_byte_flags = seq_compressed.len() + byte_flags.len();
    let total_bit_flags = seq_compressed.len() + bit_flags.len();
    let notes = format!(
        "rc={:.1}% cmp={} rc_wins={} bit_net={}",
        pct(stats.rc_reads, sequences.len()),
        stats.comparable_reads,
        stats.rc_wins,
        humanize(total_bit_flags),
    );
    print_row(
        label,
        seq_compressed.len(),
        Some((byte_flags.len(), bit_flags.len())),
        total_byte_flags,
        total_bases,
        baseline,
        transform_time + build_time + bsc_time,
        &notes,
    );
    Ok(())
}

fn bench_delta_stream(
    label: &str,
    stream: &[u8],
    total_bases: usize,
    num_reads: usize,
    baseline: Baseline,
    transform_time: Duration,
    stats: &OrientationStats,
) -> Result<()> {
    let t = Instant::now();
    let compressed = compress_bsc(stream)?;
    let bsc_time = t.elapsed();
    let notes = format!(
        "delta={:.1}% rc={:.1}% cmp={} rc_wins={}",
        pct(stats.delta_reads, num_reads),
        pct(stats.rc_reads, num_reads),
        stats.comparable_reads,
        stats.rc_wins,
    );
    print_row(
        label,
        compressed.len(),
        None,
        compressed.len(),
        total_bases,
        baseline,
        transform_time + bsc_time,
        &notes,
    );
    Ok(())
}

fn bench_semisorted_window(
    label: &str,
    sequences: &[Vec<u8>],
    const_len: Option<usize>,
    window: usize,
    total_bases: usize,
    baseline: Baseline,
) -> Result<()> {
    let t = Instant::now();
    let (stream, rank_bits) = build_semisorted_sequence_stream(sequences, const_len, window);
    let mut build_time = t.elapsed();

    let order_stream = if window <= 32 {
        let t = Instant::now();
        let stream = build_semisorted_order_rank_stream(sequences, window);
        build_time += t.elapsed();
        Some(stream)
    } else {
        None
    };

    let t = Instant::now();
    let seq_compressed = compress_bsc(&stream)?;
    let mut bsc_time = t.elapsed();

    let oracle_bytes = (rank_bits / 8.0).ceil() as usize;
    if let Some((order_stream, order_bits)) = order_stream {
        let t = Instant::now();
        let order_compressed = compress_bsc(&order_stream)?;
        bsc_time += t.elapsed();
        let order_bytes = order_compressed.len();
        let net = seq_compressed.len() + order_bytes;
        let notes = format!(
            "order_bsc={} raw={} oracle={} raw_bits={} ({:.3} oracle bpb)",
            humanize(order_bytes),
            humanize(order_stream.len()),
            humanize(oracle_bytes),
            order_bits,
            rank_bits / total_bases.max(1) as f64,
        );
        print_row(
            label,
            seq_compressed.len(),
            Some((order_bytes, oracle_bytes)),
            net,
            total_bases,
            baseline,
            build_time + bsc_time,
            &notes,
        );
    } else {
        let net = seq_compressed.len() + oracle_bytes;
        let notes = format!(
            "oracle only={} ({:.3} bpb); actual rank encoder limited to w<=32",
            humanize(oracle_bytes),
            rank_bits / total_bases.max(1) as f64,
        );
        print_row(
            label,
            seq_compressed.len(),
            Some((oracle_bytes, oracle_bytes)),
            net,
            total_bases,
            baseline,
            build_time + bsc_time,
            &notes,
        );
    }
    Ok(())
}

fn print_row(
    label: &str,
    seq_bsc: usize,
    flags: Option<(usize, usize)>,
    net: usize,
    total_bases: usize,
    baseline: Baseline,
    elapsed: Duration,
    notes: &str,
) {
    let ratio = total_bases as f64 / net.max(1) as f64;
    let delta = (net as i64) - (baseline.compressed_len as i64);
    let flag_text = flags
        .map(|(byte, bit)| format!("{}/{}", humanize(byte), humanize(bit)))
        .unwrap_or_else(|| "-".to_string());
    println!(
        "{:<34} {:>12} {:>15} {:>12} {:>8.2} {:>+8.2}% {:>10.2}s  {}",
        label,
        humanize(seq_bsc),
        flag_text,
        humanize(net),
        ratio,
        100.0 * delta as f64 / baseline.compressed_len.max(1) as f64,
        elapsed.as_secs_f64(),
        notes,
    );
}

fn build_sequence_stream(sequences: &[Vec<u8>], const_len: Option<usize>) -> Vec<u8> {
    let include_len = const_len.is_none();
    let mut out = Vec::with_capacity(
        sequences.iter().map(Vec::len).sum::<usize>()
            + if include_len { sequences.len() * 2 } else { 0 },
    );
    for seq in sequences {
        if include_len {
            write_varint(&mut out, seq.len());
        }
        out.extend_from_slice(seq);
    }
    out
}

fn build_semisorted_sequence_stream(
    sequences: &[Vec<u8>],
    const_len: Option<usize>,
    window: usize,
) -> (Vec<u8>, f64) {
    let include_len = const_len.is_none();
    let mut out = Vec::with_capacity(
        sequences.iter().map(Vec::len).sum::<usize>()
            + if include_len { sequences.len() * 2 } else { 0 },
    );
    let mut rank_bits = 0.0;
    let window = window.max(1);

    for seq in sequences {
        if include_len {
            write_varint(&mut out, seq.len());
        }
        for chunk in seq.chunks(window) {
            rank_bits += multinomial_rank_bits_dna(chunk);
            let mut sorted = chunk.to_vec();
            sorted.sort_by_key(|base| dna_sort_key(*base));
            out.extend_from_slice(&sorted);
        }
    }

    (out, rank_bits)
}

fn build_semisorted_order_rank_stream(sequences: &[Vec<u8>], window: usize) -> (Vec<u8>, usize) {
    let mut writer = BitWriter::default();
    let window = window.max(1);
    for seq in sequences {
        for chunk in seq.chunks(window) {
            let (rank, width) = multiset_rank_and_width(chunk);
            writer.write_bits_msb(rank, width);
        }
    }
    writer.into_parts()
}

#[derive(Default)]
struct BitWriter {
    bytes: Vec<u8>,
    bit_len: usize,
}

impl BitWriter {
    fn write_bits_msb(&mut self, value: u128, bits: u32) {
        debug_assert!(bits <= 128);
        for shift in (0..bits).rev() {
            if self.bit_len.is_multiple_of(8) {
                self.bytes.push(0);
            }
            if ((value >> shift) & 1) != 0 {
                let byte_idx = self.bit_len / 8;
                let bit_idx = 7 - (self.bit_len % 8);
                self.bytes[byte_idx] |= 1 << bit_idx;
            }
            self.bit_len += 1;
        }
    }

    fn into_parts(self) -> (Vec<u8>, usize) {
        (self.bytes, self.bit_len)
    }
}

fn multiset_rank_and_width(seq: &[u8]) -> (u128, u32) {
    let mut counts = [0usize; 5];
    for &base in seq {
        counts[dna_sort_key(base) as usize] += 1;
    }
    let total_perms = multinomial_count_u128(&counts);
    let width = ceil_log2_u128(total_perms);
    let rank = multiset_rank_u128(seq, counts);
    debug_assert!(rank < total_perms);
    (rank, width)
}

fn multiset_rank_u128(seq: &[u8], mut counts: [usize; 5]) -> u128 {
    let mut rank = 0u128;
    for &base in seq {
        let key = dna_sort_key(base) as usize;
        for smaller in 0..key {
            if counts[smaller] == 0 {
                continue;
            }
            counts[smaller] -= 1;
            rank += multinomial_count_u128(&counts);
            counts[smaller] += 1;
        }
        debug_assert!(counts[key] > 0);
        counts[key] -= 1;
    }
    rank
}

fn multinomial_count_u128(counts: &[usize; 5]) -> u128 {
    let n: usize = counts.iter().sum();
    debug_assert!(n <= 32);
    let mut value = factorial_u128(n);
    for &count in counts {
        let denom = factorial_u128(count);
        debug_assert_eq!(value % denom, 0);
        value /= denom;
    }
    value
}

fn factorial_u128(n: usize) -> u128 {
    (2..=n).fold(1u128, |acc, value| acc * value as u128)
}

fn ceil_log2_u128(value: u128) -> u32 {
    if value <= 1 {
        0
    } else {
        128 - (value - 1).leading_zeros()
    }
}

fn multinomial_rank_bits_dna(seq: &[u8]) -> f64 {
    let mut counts = [0usize; 5];
    for &base in seq {
        counts[dna_sort_key(base) as usize] += 1;
    }
    let mut bits = log2_factorial(seq.len());
    for count in counts {
        bits -= log2_factorial(count);
    }
    bits
}

fn log2_factorial(n: usize) -> f64 {
    (2..=n).map(|value| (value as f64).log2()).sum()
}

fn dna_sort_key(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

fn compress_bsc(data: &[u8]) -> Result<Vec<u8>> {
    if data.is_empty() {
        Ok(Vec::new())
    } else {
        bsc::compress_parallel_adaptive_no_lzp(data)
    }
}

fn orient_lexicographic(sequences: &[Vec<u8>]) -> (Vec<Vec<u8>>, Vec<u8>, OrientationStats) {
    let mut oriented = Vec::with_capacity(sequences.len());
    let mut flags = Vec::with_capacity(sequences.len());
    let mut stats = OrientationStats::default();
    for seq in sequences {
        let (canon, was_rc) = dna_utils::canonicalize_sequence(seq);
        stats.rc_reads += was_rc as usize;
        flags.push(was_rc as u8);
        oriented.push(canon);
    }
    (oriented, flags, stats)
}

fn orient_minimizer_strand(sequences: &[Vec<u8>]) -> (Vec<Vec<u8>>, Vec<u8>, OrientationStats) {
    let mut oriented = Vec::with_capacity(sequences.len());
    let mut flags = Vec::with_capacity(sequences.len());
    let mut stats = OrientationStats::default();
    for seq in sequences {
        let use_rc = min_syncmer_hash_and_strand(seq)
            .map(|(_, rc)| rc)
            .unwrap_or(false);
        stats.rc_reads += use_rc as usize;
        flags.push(use_rc as u8);
        oriented.push(if use_rc {
            dna_utils::reverse_complement(seq)
        } else {
            seq.clone()
        });
    }
    (oriented, flags, stats)
}

fn orient_bucket_representative(
    sequences: &[Vec<u8>],
) -> (Vec<Vec<u8>>, Vec<u8>, OrientationStats) {
    let mut reps: FxHashMap<u64, Vec<u8>> = FxHashMap::default();
    let mut oriented = Vec::with_capacity(sequences.len());
    let mut flags = Vec::with_capacity(sequences.len());
    let mut stats = OrientationStats::default();

    for seq in sequences {
        let (hash, fallback_rc) = min_syncmer_hash_and_strand(seq).unwrap_or((u64::MAX, false));
        let mut use_rc = fallback_rc;

        if hash != u64::MAX
            && let Some(rep) = reps.get(&hash)
            && rep.len() == seq.len()
        {
            stats.comparable_reads += 1;
            let fwd = matching_bases(seq, false, rep);
            let rc = matching_bases(seq, true, rep);
            if rc > fwd {
                use_rc = true;
                stats.rc_wins += 1;
            } else if fwd > rc {
                use_rc = false;
            }
        }

        stats.rc_reads += use_rc as usize;
        flags.push(use_rc as u8);
        let encoded = if use_rc {
            dna_utils::reverse_complement(seq)
        } else {
            seq.clone()
        };
        if hash != u64::MAX {
            reps.insert(hash, encoded.clone());
        }
        oriented.push(encoded);
    }

    (oriented, flags, stats)
}

fn orient_adjacent_dp(sequences: &[Vec<u8>]) -> (Vec<Vec<u8>>, Vec<u8>, OrientationStats) {
    let n = sequences.len();
    if n == 0 {
        return (Vec::new(), Vec::new(), OrientationStats::default());
    }

    let mut back = vec![[false; 2]; n];
    let mut prev = [0i64, 0i64];

    for i in 1..n {
        let mut next = [0i64, 0i64];
        for curr_state in 0..=1 {
            let score_from_fwd =
                prev[0] + pair_score(&sequences[i - 1], false, &sequences[i], curr_state == 1);
            let score_from_rc =
                prev[1] + pair_score(&sequences[i - 1], true, &sequences[i], curr_state == 1);
            if score_from_rc > score_from_fwd {
                next[curr_state] = score_from_rc;
                back[i][curr_state] = true;
            } else {
                next[curr_state] = score_from_fwd;
                back[i][curr_state] = false;
            }
        }
        prev = next;
    }

    let mut states = vec![false; n];
    states[n - 1] = prev[1] > prev[0];
    for i in (1..n).rev() {
        states[i - 1] = back[i][states[i] as usize];
    }

    let mut oriented = Vec::with_capacity(n);
    let mut flags = Vec::with_capacity(n);
    let mut stats = OrientationStats::default();
    for (seq, use_rc) in sequences.iter().zip(states) {
        stats.rc_reads += use_rc as usize;
        flags.push(use_rc as u8);
        oriented.push(if use_rc {
            dna_utils::reverse_complement(seq)
        } else {
            seq.clone()
        });
    }
    stats.comparable_reads = n.saturating_sub(1);
    (oriented, flags, stats)
}

fn build_delta_stream(sequences: &[Vec<u8>], rc_aware: bool) -> (Vec<u8>, OrientationStats) {
    let mut stream = Vec::with_capacity(sequences.iter().map(Vec::len).sum::<usize>());
    let mut cache: FxHashMap<u64, usize> = FxHashMap::default();
    let mut stored_reads: Vec<Vec<u8>> = Vec::with_capacity(sequences.len());
    let mut stats = OrientationStats::default();

    for seq in sequences {
        let (hash, fallback_rc) = min_syncmer_hash_and_strand(seq).unwrap_or((u64::MAX, false));
        let mut use_rc = rc_aware && fallback_rc;
        let mut ref_idx = None;

        if hash != u64::MAX
            && let Some(&idx) = cache.get(&hash)
        {
            let ref_seq = &stored_reads[idx];
            if ref_seq.len() == seq.len() {
                stats.comparable_reads += 1;
                let fwd_matches = matching_bases(seq, false, ref_seq);
                let rc_matches = if rc_aware {
                    matching_bases(seq, true, ref_seq)
                } else {
                    0
                };
                let (matches, rc_won) = if rc_matches > fwd_matches {
                    (rc_matches, true)
                } else {
                    (fwd_matches, false)
                };
                if rc_won {
                    use_rc = true;
                    stats.rc_wins += 1;
                } else {
                    use_rc = false;
                }
                if (matches as f64 / seq.len().max(1) as f64) > DELTA_IDENTITY_THRESHOLD {
                    ref_idx = Some(idx);
                }
            }
        }

        let encoded = if use_rc {
            dna_utils::reverse_complement(seq)
        } else {
            seq.clone()
        };
        stats.rc_reads += use_rc as usize;

        write_varint(&mut stream, encoded.len());
        match ref_idx {
            Some(idx) => {
                stats.delta_reads += 1;
                // Proposed flag extension: bit0 = delta, bit1 = stored RC.
                stream.push(0x01 | ((use_rc as u8) << 1));
                write_varint(&mut stream, idx);
                let ref_seq = &stored_reads[idx];
                for (base, ref_base) in encoded.iter().zip(ref_seq.iter()) {
                    stream.push(if base == ref_base { 0x00 } else { *base });
                }
            }
            None => {
                stream.push((use_rc as u8) << 1);
                stream.extend_from_slice(&encoded);
            }
        }

        if hash != u64::MAX {
            cache.insert(hash, stored_reads.len());
        }
        stored_reads.push(encoded);
    }

    (stream, stats)
}

fn min_syncmer_hash_and_strand(seq: &[u8]) -> Option<(u64, bool)> {
    if seq.len() < K {
        return None;
    }

    // Open-syncmer selection is orientation-dependent: applying it only to the
    // forward read and canonicalizing each k-mer does not guarantee that a read
    // and its reverse complement land in the same bucket. For orientation
    // experiments we need a strand-symmetric key, so scan both full orientations
    // and choose the smaller minimizer.
    let fwd = min_syncmer_hash_one_orientation(seq);
    let rc_seq = dna_utils::reverse_complement(seq);
    let rc = min_syncmer_hash_one_orientation(&rc_seq);

    match (fwd, rc) {
        (Some(f), Some(r)) if r < f => Some((r, true)),
        (Some(f), Some(_)) => Some((f, false)),
        (None, Some(r)) => Some((r, true)),
        (Some(f), None) => Some((f, false)),
        (None, None) => None,
    }
}

fn min_syncmer_hash_one_orientation(seq: &[u8]) -> Option<u64> {
    let mut min_hash = u64::MAX;
    for pos in syncmers::find_syncmers_pos(K, S, &[0], seq) {
        if pos + K > seq.len() {
            continue;
        }
        let kmer = &seq[pos..pos + K];
        if let Some(hash) = dna_utils::kmer_to_hash(kmer) {
            min_hash = min_hash.min(hash);
        }
    }
    (min_hash != u64::MAX).then_some(min_hash)
}

fn matching_bases(read: &[u8], read_rc: bool, reference: &[u8]) -> usize {
    read.iter()
        .enumerate()
        .filter(|(i, _)| oriented_base(read, *i, read_rc) == reference[*i])
        .count()
}

fn pair_score(a: &[u8], a_rc: bool, b: &[u8], b_rc: bool) -> i64 {
    let len = a.len().min(b.len());
    let matches = (0..len)
        .filter(|&i| oriented_base(a, i, a_rc) == oriented_base(b, i, b_rc))
        .count();
    matches as i64 - a.len().abs_diff(b.len()) as i64
}

fn oriented_base(seq: &[u8], idx: usize, rc: bool) -> u8 {
    if rc {
        complement(seq[seq.len() - 1 - idx])
    } else {
        seq[idx]
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

fn write_varint(out: &mut Vec<u8>, mut value: usize) {
    while value >= 0x80 {
        out.push(((value & 0x7F) | 0x80) as u8);
        value >>= 7;
    }
    out.push(value as u8);
}

fn pack_bits(flags: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8; flags.len().div_ceil(8)];
    for (i, &flag) in flags.iter().enumerate() {
        if flag != 0 {
            out[i / 8] |= 1 << (i % 8);
        }
    }
    out
}

fn humanize(bytes: usize) -> String {
    if bytes >= 1_048_576 {
        format!("{:.1}M", bytes as f64 / 1_048_576.0)
    } else if bytes >= 1024 {
        format!("{:.1}K", bytes as f64 / 1024.0)
    } else {
        bytes.to_string()
    }
}

fn pct(part: usize, total: usize) -> f64 {
    100.0 * part as f64 / total.max(1) as f64
}
