/// Benchmark FASTR-inspired base+quality joint byte streams with BSC.
///
/// Usage:
///   cargo run -p qz-bench --bin bench_fastr_bsc -- <FASTQ> [max_reads]
use std::env;
use std::time::{Duration, Instant};

use anyhow::{Context, Result};
use qz_lib::compression::bsc;
use qz_lib::io::fastq::FastqReader;

fn main() -> Result<()> {
    let mut args = env::args().skip(1);
    let path = args
        .next()
        .unwrap_or_else(|| "real_data/ERR3239334_1.500k.fastq".to_string());
    let max_reads = args
        .next()
        .map(|s| s.parse::<usize>())
        .transpose()
        .context("max_reads must be an integer")?
        .unwrap_or(usize::MAX);

    eprintln!("Loading {path} ...");
    let t = Instant::now();
    let mut reader = FastqReader::from_path(&path, false).context("open FASTQ")?;
    let mut sequences = Vec::new();
    let mut qualities = Vec::new();
    while let Some(record) = reader.next().context("read FASTQ")? {
        sequences.push(record.sequence);
        qualities.push(record.quality.unwrap_or_default());
        if sequences.len() >= max_reads {
            break;
        }
    }

    anyhow::ensure!(!sequences.is_empty(), "no records read");
    anyhow::ensure!(
        sequences.len() == qualities.len(),
        "sequence/quality count mismatch"
    );
    let read_len = sequences[0].len();
    let const_len = sequences.iter().all(|seq| seq.len() == read_len)
        && qualities.iter().all(|qual| qual.len() == read_len);
    let total_bases: usize = sequences.iter().map(Vec::len).sum();
    let max_q = qualities.iter().flatten().copied().max().unwrap_or(0);
    let min_q = qualities.iter().flatten().copied().min().unwrap_or(0);
    let n_bases = sequences
        .iter()
        .flatten()
        .filter(|&&base| base_to_2bit(base).is_none())
        .count();
    let q_over_63 = qualities
        .iter()
        .flatten()
        .filter(|&&q| q.saturating_sub(33) > 63)
        .count();

    eprintln!(
        "Loaded {} reads, {} bases ({:.1} MiB), read_len={} {}, q_ascii={}..{}, N/other={} q>63={} in {:.2}s\n",
        sequences.len(),
        total_bases,
        total_bases as f64 / 1_048_576.0,
        read_len,
        if const_len { "uniform" } else { "variable" },
        min_q,
        max_q,
        n_bases,
        q_over_63,
        t.elapsed().as_secs_f64(),
    );

    println!(
        "{:<34} {:>12} {:>12} {:>9} {:>9}  Notes",
        "Stream", "Raw", "BSC", "bits/base", "Time"
    );
    println!("{}", "-".repeat(96));

    let seq_ascii = concat(&sequences);
    let qual_ascii = concat(&qualities);
    let seq_bsc = bench("seq ASCII", &seq_ascii, total_bases, "")?;
    let qual_bsc = bench("qual ASCII", &qual_ascii, total_bases, "")?;
    println!(
        "{:<34} {:>12} {:>12} {:>9.3} {:>9}  baseline for direct BSC streams",
        "separate seq+qual ASCII",
        humanize(seq_ascii.len() + qual_ascii.len()),
        humanize(seq_bsc.size + qual_bsc.size),
        8.0 * (seq_bsc.size + qual_bsc.size) as f64 / total_bases.max(1) as f64,
        format_duration(seq_bsc.elapsed + qual_bsc.elapsed),
    );

    let (seq_2bit, nmask) = pack_sequences_2bit_with_nmask(&sequences);
    let seq_2bit_bsc = bench("seq 2-bit packed", &seq_2bit, total_bases, "")?;
    let nmask_bsc = bench(
        "seq N-mask",
        &nmask,
        total_bases,
        "1 bit/base-ish; empty if no Ns",
    )?;
    println!(
        "{:<34} {:>12} {:>12} {:>9.3} {:>9}  closer to qz's separated-body idea",
        "separate 2bit-seq+qual ASCII",
        humanize(seq_2bit.len() + nmask.len() + qual_ascii.len()),
        humanize(seq_2bit_bsc.size + nmask_bsc.size + qual_bsc.size),
        8.0 * (seq_2bit_bsc.size + nmask_bsc.size + qual_bsc.size) as f64
            / total_bases.max(1) as f64,
        format_duration(seq_2bit_bsc.elapsed + nmask_bsc.elapsed + qual_bsc.elapsed),
    );

    let joint_q_major = build_joint_stream(&sequences, &qualities, JointMode::QMajor);
    bench(
        "joint q-major row",
        &joint_q_major.main,
        total_bases,
        &joint_q_major.notes(),
    )?;

    let joint_base_major = build_joint_stream(&sequences, &qualities, JointMode::BaseMajor);
    bench(
        "joint base-major row",
        &joint_base_major.main,
        total_bases,
        &joint_base_major.notes(),
    )?;

    let joint_fastr_ranges = build_joint_stream(&sequences, &qualities, JointMode::FastrRanges);
    bench(
        "joint FASTR ranges row",
        &joint_fastr_ranges.main,
        total_bases,
        &joint_fastr_ranges.notes(),
    )?;

    let interleaved = build_interleaved_base_quality(&sequences, &qualities);
    bench(
        "base,qual interleaved row",
        &interleaved,
        total_bases,
        "2 bytes/base",
    )?;

    if const_len {
        let qual_column = transpose_fixed(&qualities, read_len);
        let qual_col_bsc = bench("qual ASCII column-major", &qual_column, total_bases, "")?;

        let seq_column = transpose_fixed(&sequences, read_len);
        let seq_col_bsc = bench("seq ASCII column-major", &seq_column, total_bases, "")?;
        println!(
            "{:<34} {:>12} {:>12} {:>9.3} {:>9}  same data, cycle-major order",
            "separate col seq+qual ASCII",
            humanize(seq_column.len() + qual_column.len()),
            humanize(seq_col_bsc.size + qual_col_bsc.size),
            8.0 * (seq_col_bsc.size + qual_col_bsc.size) as f64 / total_bases.max(1) as f64,
            format_duration(seq_col_bsc.elapsed + qual_col_bsc.elapsed),
        );

        let joint_q_col =
            transpose_joint_fixed(&sequences, &qualities, read_len, JointMode::QMajor);
        bench(
            "joint q-major column",
            &joint_q_col.main,
            total_bases,
            &joint_q_col.notes(),
        )?;

        let joint_base_col =
            transpose_joint_fixed(&sequences, &qualities, read_len, JointMode::BaseMajor);
        bench(
            "joint base-major column",
            &joint_base_col.main,
            total_bases,
            &joint_base_col.notes(),
        )?;
    }

    Ok(())
}

#[derive(Clone, Copy)]
enum JointMode {
    QMajor,
    BaseMajor,
    FastrRanges,
}

struct BenchResult {
    size: usize,
    elapsed: Duration,
}

struct JointStream {
    main: Vec<u8>,
    escapes: usize,
}

impl JointStream {
    fn notes(&self) -> String {
        if self.escapes == 0 {
            "exact for A/C/G/T and q<=63; no escapes".to_string()
        } else {
            format!(
                "{} escaped base/quality pairs not included in net size",
                self.escapes
            )
        }
    }
}

fn bench(label: &str, data: &[u8], total_bases: usize, notes: &str) -> Result<BenchResult> {
    let t = Instant::now();
    let compressed = compress_bsc(data)?;
    let elapsed = t.elapsed();
    println!(
        "{:<34} {:>12} {:>12} {:>9.3} {:>9}  {}",
        label,
        humanize(data.len()),
        humanize(compressed.len()),
        8.0 * compressed.len() as f64 / total_bases.max(1) as f64,
        format_duration(elapsed),
        notes,
    );
    Ok(BenchResult {
        size: compressed.len(),
        elapsed,
    })
}

fn compress_bsc(data: &[u8]) -> Result<Vec<u8>> {
    if data.is_empty() {
        Ok(Vec::new())
    } else {
        bsc::compress_parallel_adaptive(data)
    }
}

fn concat(records: &[Vec<u8>]) -> Vec<u8> {
    records
        .iter()
        .flat_map(|record| record.iter().copied())
        .collect()
}

fn build_joint_stream(
    sequences: &[Vec<u8>],
    qualities: &[Vec<u8>],
    mode: JointMode,
) -> JointStream {
    let mut main = Vec::with_capacity(sequences.iter().map(Vec::len).sum());
    let mut escapes = 0usize;
    for (seq, qual) in sequences.iter().zip(qualities) {
        for (&base, &q_ascii) in seq.iter().zip(qual) {
            if let Some(byte) = encode_joint(base, q_ascii, mode) {
                main.push(byte);
            } else {
                escapes += 1;
                main.push(255);
            }
        }
    }
    JointStream { main, escapes }
}

fn transpose_joint_fixed(
    sequences: &[Vec<u8>],
    qualities: &[Vec<u8>],
    read_len: usize,
    mode: JointMode,
) -> JointStream {
    let mut main = Vec::with_capacity(sequences.len() * read_len);
    let mut escapes = 0usize;
    for pos in 0..read_len {
        for (seq, qual) in sequences.iter().zip(qualities) {
            if let Some(byte) = encode_joint(seq[pos], qual[pos], mode) {
                main.push(byte);
            } else {
                escapes += 1;
                main.push(255);
            }
        }
    }
    JointStream { main, escapes }
}

fn encode_joint(base: u8, q_ascii: u8, mode: JointMode) -> Option<u8> {
    let base_code = base_to_2bit(base)?;
    let q = q_ascii.checked_sub(33)?;
    if q > 63 {
        return None;
    }
    match mode {
        JointMode::QMajor => Some(q * 4 + base_code),
        JointMode::BaseMajor => Some(base_code * 64 + q),
        JointMode::FastrRanges => {
            // FASTR's paper uses base-specific contiguous ranges and reserves 255.
            // Its N behavior is intentionally not modeled here; qz already separates N-masks.
            let lower = match base {
                b'A' | b'a' => 3,
                b'G' | b'g' => 66,
                b'C' | b'c' => 129,
                b'T' | b't' => 192,
                _ => return None,
            };
            Some(lower + q.min(62))
        }
    }
}

fn base_to_2bit(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn build_interleaved_base_quality(sequences: &[Vec<u8>], qualities: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::with_capacity(sequences.iter().map(Vec::len).sum::<usize>() * 2);
    for (seq, qual) in sequences.iter().zip(qualities) {
        for (&base, &q) in seq.iter().zip(qual) {
            out.push(base);
            out.push(q);
        }
    }
    out
}

fn transpose_fixed(records: &[Vec<u8>], read_len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(records.len() * read_len);
    for pos in 0..read_len {
        for record in records {
            out.push(record[pos]);
        }
    }
    out
}

fn pack_sequences_2bit_with_nmask(sequences: &[Vec<u8>]) -> (Vec<u8>, Vec<u8>) {
    let total_bases: usize = sequences.iter().map(Vec::len).sum();
    let mut packed = Vec::with_capacity(total_bases.div_ceil(4));
    let mut nmask = vec![0u8; total_bases.div_ceil(8)];
    let mut current = 0u8;
    let mut in_byte = 0usize;
    let mut absolute = 0usize;

    for seq in sequences {
        for &base in seq {
            let code = base_to_2bit(base).unwrap_or(0);
            current |= code << (6 - 2 * in_byte);
            if base_to_2bit(base).is_none() {
                nmask[absolute / 8] |= 1 << (absolute % 8);
            }
            in_byte += 1;
            absolute += 1;
            if in_byte == 4 {
                packed.push(current);
                current = 0;
                in_byte = 0;
            }
        }
    }
    if in_byte != 0 {
        packed.push(current);
    }
    if nmask.iter().all(|&byte| byte == 0) {
        nmask.clear();
    }
    (packed, nmask)
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

fn format_duration(duration: Duration) -> String {
    if duration.as_secs_f64() >= 1.0 {
        format!("{:.2}s", duration.as_secs_f64())
    } else {
        format!("{:.1}ms", duration.as_secs_f64() * 1000.0)
    }
}
