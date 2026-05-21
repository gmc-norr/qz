//! Smoke tests for the `qz` binary.
//!
//! These exercise the surface the user actually invokes — flag parsing, exit
//! codes, conflict detection, and output overwrite protection — none of which
//! the library-level integration tests cover.

use assert_cmd::Command;
use predicates::prelude::*;
use std::fs;
use tempfile::TempDir;

/// Returns a `Command` that invokes the compiled `qz` binary.
fn qz() -> Command {
    Command::cargo_bin("qz").expect("qz binary not built")
}

/// Standard 2-record FASTQ for quick roundtrips.
const TINY_FASTQ: &str = "@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n";

#[test]
fn help_shows_subcommands() {
    qz()
        .arg("--help")
        .assert()
        .success()
        .stdout(predicate::str::contains("compress").and(predicate::str::contains("decompress")));
}

#[test]
fn version_prints_a_version_string() {
    qz()
        .arg("--version")
        .assert()
        .success()
        .stdout(predicate::str::contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn missing_input_errors_with_helpful_message() {
    qz()
        .args(["compress", "-i", "/nonexistent/path.fastq", "-o", "/tmp/out.qz"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Failed to open file"));
}

#[test]
fn refuses_to_overwrite_without_force() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.fastq");
    let output = dir.path().join("out.qz");
    fs::write(&input, TINY_FASTQ).unwrap();

    // First compress: succeeds.
    qz()
        .args(["compress", "-i"])
        .arg(&input)
        .args(["-o"])
        .arg(&output)
        .args(["-w"])
        .arg(dir.path())
        .args(["-t", "1"])
        .assert()
        .success();

    // Second compress to same output without --force: must fail.
    qz()
        .args(["compress", "-i"])
        .arg(&input)
        .args(["-o"])
        .arg(&output)
        .args(["-w"])
        .arg(dir.path())
        .args(["-t", "1"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("already exists"));

    // With --force: succeeds again.
    qz()
        .args(["compress", "-i"])
        .arg(&input)
        .args(["-o"])
        .arg(&output)
        .args(["-w"])
        .arg(dir.path())
        .args(["-t", "1", "--force"])
        .assert()
        .success();
}

#[test]
fn fasta_and_lossy_quality_conflict() {
    // --fasta means input has no quality data; a lossy --quality-mode is then
    // nonsensical and must be rejected at flag-parse time.
    qz()
        .args(["compress", "--fasta", "--quality-mode", "discard",
               "-i", "/tmp/dummy.fastq", "-o", "/tmp/dummy.qz"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--fasta is incompatible"));
}

#[test]
fn no_quality_and_quality_mode_conflict() {
    qz()
        .args(["compress", "--no-quality", "--quality-mode", "illumina-bin",
               "-i", "/tmp/dummy.fastq", "-o", "/tmp/dummy.qz"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--no-quality is incompatible"));
}

#[test]
fn end_to_end_compress_decompress_roundtrip() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.fastq");
    let archive = dir.path().join("out.qz");
    let restored = dir.path().join("restored.fastq");
    fs::write(&input, TINY_FASTQ).unwrap();

    qz()
        .args(["compress", "-i"])
        .arg(&input)
        .args(["-o"])
        .arg(&archive)
        .args(["-w"])
        .arg(dir.path())
        .args(["-t", "1"])
        .assert()
        .success();

    qz()
        .args(["decompress", "-i"])
        .arg(&archive)
        .args(["-o"])
        .arg(&restored)
        .args(["-w"])
        .arg(dir.path())
        .args(["-t", "1"])
        .assert()
        .success();

    let original_bytes = fs::read(&input).unwrap();
    let restored_bytes = fs::read(&restored).unwrap();
    assert_eq!(original_bytes, restored_bytes, "decompressed bytes must match original");
}

#[test]
fn verify_reports_ok_for_valid_archive() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.fastq");
    let archive = dir.path().join("out.qz");
    fs::write(&input, TINY_FASTQ).unwrap();

    qz()
        .args(["compress", "-i"])
        .arg(&input)
        .args(["-o"])
        .arg(&archive)
        .args(["-w"])
        .arg(dir.path())
        .args(["-t", "1"])
        .assert()
        .success();

    qz()
        .args(["verify", "-i"])
        .arg(&archive)
        .args(["-w"])
        .arg(dir.path())
        .args(["-t", "1"])
        .assert()
        .success()
        .stderr(predicate::str::contains("Status:      OK"));
}
