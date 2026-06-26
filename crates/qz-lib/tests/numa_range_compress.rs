//! `compress_byte_range` single-end: a byte range compresses to a valid part
//! archive whose decoded reads are exactly the records in that range, and a
//! forced const-length violation in a shard is rejected.

use qz_lib::cli::{AdvancedOptions, CompressConfig, DecompressConfig};
use qz_lib::compression::{ByteRange, compress_byte_range, resolve_plan_override};
use std::path::Path;

fn fixed_fastq(n: usize, len: usize, start: usize) -> String {
    let seq = "A".repeat(len);
    let q = "I".repeat(len);
    (0..n).map(|i| format!("@r{}\n{seq}\n+\n{q}\n", start + i)).collect()
}

fn base_cfg(input: &Path, out: &Path, d: &Path) -> CompressConfig {
    CompressConfig {
        input: vec![input.to_path_buf()],
        output: out.to_path_buf(),
        working_dir: d.to_path_buf(),
        threads: 1,
        advanced: AdvancedOptions { chunk_records: 8, ..AdvancedOptions::default() },
        ..CompressConfig::default()
    }
}

fn decompress_text(d: &Path, archive: &Path, name: &str) -> String {
    let out = d.join(format!("{name}.fastq"));
    let cfg = DecompressConfig {
        input: archive.to_path_buf(),
        output: vec![out.clone()],
        working_dir: d.to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&cfg).unwrap();
    std::fs::read_to_string(&out).unwrap()
}

#[test]
fn compress_byte_range_emits_valid_range_archive() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let text = fixed_fastq(40, 50, 0);
    let input = d.join("in.fastq");
    std::fs::write(&input, &text).unwrap();

    // Whole-file range == in-process compress (decoded text identical).
    let filesize = std::fs::metadata(&input).unwrap().len();
    let out = d.join("range.qz");
    let cfg = base_cfg(&input, &out, d);
    let plan = resolve_plan_override(&cfg).unwrap();
    compress_byte_range(&cfg, &[ByteRange { start: 0, end: filesize }], &plan).unwrap();
    assert_eq!(decompress_text(d, &out, "r"), text);
}

#[test]
fn compress_byte_range_forces_first_chunk_const_validation() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    // Parent plan declares const-50 (from a 50bp head). The worker's range starts
    // on a 51bp record — the first-chunk validation MUST reject it (else silent
    // byte-shift). Build an input whose head is 50bp; take a range that begins at a
    // later 51bp record.
    let head = fixed_fastq(8, 50, 0);
    let tail_51 = {
        let seq = "C".repeat(51);
        let q = "I".repeat(51);
        (0..8).map(|i| format!("@t{i}\n{seq}\n+\n{q}\n", )).collect::<String>()
    };
    let text = format!("{head}{tail_51}");
    let input = d.join("mixed.fastq");
    std::fs::write(&input, &text).unwrap();

    let out = d.join("mixed.qz");
    let cfg = base_cfg(&input, &out, d);
    // Pin const-50 (what the parent head resolves).
    let plan = resolve_plan_override(&cfg).unwrap();
    assert_eq!(plan.const_seq_len, 50);

    // A range covering ONLY the 51bp tail must be rejected at first-chunk validation.
    let head_bytes = head.len() as u64;
    let filesize = std::fs::metadata(&input).unwrap().len();
    let err = compress_byte_range(&cfg, &[ByteRange { start: head_bytes, end: filesize }], &plan)
        .unwrap_err()
        .to_string();
    assert!(err.contains("length"), "expected a const-length mismatch, got: {err}");
}
