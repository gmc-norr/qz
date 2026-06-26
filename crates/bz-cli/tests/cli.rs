//! Smoke tests for the `bz` binary.
//!
//! BAM I/O end-to-end isn't covered here (synthesising a valid BAM in-test is
//! heavy); these focus on the CLI surface itself: --help, --version, exit
//! codes on missing input, and overwrite refusal without --force.

use assert_cmd::Command;
use predicates::prelude::*;
use std::fs;
use tempfile::TempDir;

fn bz() -> Command {
    Command::cargo_bin("bz").expect("bz binary not built")
}

#[test]
fn help_lists_all_subcommands() {
    bz().arg("--help").assert().success().stdout(
        predicate::str::contains("compress")
            .and(predicate::str::contains("decompress"))
            .and(predicate::str::contains("extract"))
            .and(predicate::str::contains("verify")),
    );
}

#[test]
fn version_prints_a_version_string() {
    bz().arg("--version")
        .assert()
        .success()
        .stdout(predicate::str::contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn missing_input_errors_out() {
    bz().args(["compress", "-i", "/nonexistent/path.bam", "-o", "/tmp/x.bz"])
        .assert()
        .failure();
}

#[test]
fn refuses_to_overwrite_compress_output_without_force() {
    let dir = TempDir::new().unwrap();
    let dummy_output = dir.path().join("existing.bz");
    fs::write(&dummy_output, b"placeholder").unwrap();

    // Even though compression itself can't proceed (no real BAM input),
    // the existence check fires *before* the library is invoked.
    let dummy_input = dir.path().join("dummy.bam");
    fs::write(&dummy_input, b"placeholder").unwrap();

    bz().args(["compress", "-i"])
        .arg(&dummy_input)
        .args(["-o"])
        .arg(&dummy_output)
        .args(["-w"])
        .arg(dir.path())
        .assert()
        .failure()
        .stderr(predicate::str::contains("already exists"));
}

#[test]
fn refuses_to_overwrite_extract_output_without_force() {
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("sample");
    // Pre-create an _R1.qz to trigger the overwrite check.
    let r1 = dir.path().join("sample_R1.qz");
    fs::write(&r1, b"placeholder").unwrap();

    let dummy_input = dir.path().join("dummy.bam");
    fs::write(&dummy_input, b"placeholder").unwrap();

    bz().args(["extract", "-i"])
        .arg(&dummy_input)
        .args(["-o"])
        .arg(&prefix)
        .args(["-w"])
        .arg(dir.path())
        .assert()
        .failure()
        .stderr(predicate::str::contains("already exists"));
}
