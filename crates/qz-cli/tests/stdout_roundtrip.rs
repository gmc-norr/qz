//! Test that `qz compress -o -` (stdout) produces a valid v5 archive that
//! decompresses byte-identically to the original input.
//!
//! We use `std::process::Command` directly so we can redirect stdout to a file
//! via `Stdio::from(File::create(...))` — `assert_cmd` does not expose that.

use std::fs::{self, File};
use std::process::{Command, Stdio};
use tempfile::TempDir;

/// Path to the compiled `qz` binary (resolved at compile time by Cargo).
const QZ: &str = env!("CARGO_BIN_EXE_qz");

/// Generate a varied FASTQ with `n` records. Each record has a unique sequence
/// and quality string so the input is not byte-identical across records.
fn varied_fastq(n: usize) -> String {
    let bases = b"ACGT";
    let mut out = String::with_capacity(n * 80);
    for i in 0..n {
        // Header
        out.push('@');
        out.push_str(&format!("read{i} len=150 rg=RG1\n"));

        // 150-bp sequence: cycle through all 4 bases, scrambled by index.
        let seq: String = (0..150_usize)
            .map(|j| bases[(i.wrapping_add(j).wrapping_mul(7)) % 4] as char)
            .collect();
        out.push_str(&seq);
        out.push('\n');

        // Separator
        out.push_str("+\n");

        // 150-bp quality: cycle 33..=73 (printable ASCII '!' to 'J'), varied.
        let qual: String = (0..150_usize)
            .map(|j| (33u8 + ((i.wrapping_add(j).wrapping_mul(13)) % 41) as u8) as char)
            .collect();
        out.push_str(&qual);
        out.push('\n');
    }
    out
}

#[test]
fn v5_compress_to_stdout_roundtrips() {
    let dir = TempDir::new().expect("create tempdir");
    let input = dir.path().join("input.fastq");
    let archive = dir.path().join("archive.qz");
    let restored = dir.path().join("restored.fastq");

    // Write a varied 5000-record FASTQ.
    let fastq_data = varied_fastq(5_000);
    fs::write(&input, &fastq_data).expect("write input fastq");

    // --- compress to stdout ---
    // Redirect stdout directly into `archive` so we capture the raw binary stream.
    let archive_file = File::create(&archive).expect("create archive file");
    let status = Command::new(QZ)
        .env("QZ_NO_BANNER", "1") // suppress the banner printed to stderr
        .args([
            "compress",
            "-i",
        ])
        .arg(&input)
        .args(["-o", "-", "-t", "4"])
        .stdout(Stdio::from(archive_file))
        .stderr(Stdio::inherit()) // let progress/errors show in test output
        .status()
        .expect("spawn qz compress");
    assert!(
        status.success(),
        "qz compress -o - failed with exit code {:?}",
        status.code()
    );

    // --- assert archive is a v5 archive ---
    let archive_bytes = fs::read(&archive).expect("read archive");
    assert!(
        archive_bytes.len() >= 10,
        "archive too small ({} bytes)",
        archive_bytes.len()
    );
    // byte[0..2] = "QZ" magic
    assert_eq!(
        &archive_bytes[0..2],
        b"QZ",
        "archive does not start with QZ magic"
    );
    // byte[2] = version 5 (chunk-major)
    assert_eq!(
        archive_bytes[2], 5,
        "expected archive version 5, got {}",
        archive_bytes[2]
    );
    // last 8 bytes = "QZFOOTR1" trailing locator magic
    let tail = archive_bytes
        .len()
        .checked_sub(8)
        .expect("archive shorter than 8 bytes");
    assert_eq!(
        &archive_bytes[tail..],
        b"QZFOOTR1",
        "archive does not end with QZFOOTR1 footer magic"
    );

    // --- decompress and compare ---
    let status = Command::new(QZ)
        .env("QZ_NO_BANNER", "1")
        .args(["decompress", "-i"])
        .arg(&archive)
        .args(["-o"])
        .arg(&restored)
        .args(["-t", "4"])
        .status()
        .expect("spawn qz decompress");
    assert!(
        status.success(),
        "qz decompress failed with exit code {:?}",
        status.code()
    );

    let restored_bytes = fs::read(&restored).expect("read restored fastq");
    assert_eq!(
        fastq_data.as_bytes(),
        restored_bytes.as_slice(),
        "decompressed content does not match original"
    );
}

/// Build a varied paired FASTQ pair of `n` records. R1 and R2 share read ids
/// (mate field differs) so the R2-header delta path is exercised.
fn varied_pair(n: usize) -> (String, String) {
    let bases = b"ACGT";
    let mut r1 = String::with_capacity(n * 80);
    let mut r2 = String::with_capacity(n * 80);
    for i in 0..n {
        let seq1: String = (0..150_usize)
            .map(|j| bases[(i.wrapping_add(j).wrapping_mul(7)) % 4] as char)
            .collect();
        let seq2: String = (0..150_usize)
            .map(|j| bases[(i.wrapping_add(j).wrapping_mul(11)) % 4] as char)
            .collect();
        let qual1: String = (0..150_usize)
            .map(|j| (33u8 + ((i.wrapping_add(j).wrapping_mul(13)) % 41) as u8) as char)
            .collect();
        let qual2: String = (0..150_usize)
            .map(|j| (33u8 + ((i.wrapping_add(j).wrapping_mul(17)) % 41) as u8) as char)
            .collect();
        r1.push_str(&format!("@read{i} 1:N:0:ACGT\n{seq1}\n+\n{qual1}\n"));
        r2.push_str(&format!("@read{i} 2:N:0:ACGT\n{seq2}\n+\n{qual2}\n"));
    }
    (r1, r2)
}

/// `qz compress -i R1 -i R2 -o -` writes the seek-free v5 paired container to
/// stdout; decompressing it to a prefix reconstructs both mates byte-identically.
#[test]
fn paired_v5_compress_to_stdout_roundtrips() {
    let dir = TempDir::new().expect("create tempdir");
    let r1 = dir.path().join("R1.fastq");
    let r2 = dir.path().join("R2.fastq");
    let archive = dir.path().join("paired.qz");
    let prefix = dir.path().join("restored");

    let (s1, s2) = varied_pair(5_000);
    fs::write(&r1, &s1).expect("write R1");
    fs::write(&r2, &s2).expect("write R2");

    // --- compress two inputs to stdout ---
    let archive_file = File::create(&archive).expect("create archive file");
    let status = Command::new(QZ)
        .env("QZ_NO_BANNER", "1")
        .args(["compress", "-i"])
        .arg(&r1)
        .arg("-i")
        .arg(&r2)
        .args(["-o", "-", "-t", "4"])
        .stdout(Stdio::from(archive_file))
        .stderr(Stdio::inherit())
        .status()
        .expect("spawn qz compress (paired stdout)");
    assert!(
        status.success(),
        "qz paired compress -o - failed with exit code {:?}",
        status.code()
    );

    // --- assert it is a v5 paired archive ---
    let b = fs::read(&archive).expect("read archive");
    assert!(b.len() >= 20, "archive too small ({} bytes)", b.len());
    assert_eq!(&b[0..2], b"QZ", "missing QZ magic");
    assert_eq!(b[2], 5, "expected v5, got {}", b[2]);
    assert_eq!(b[3], 1, "expected archive_type 1 (paired), got {}", b[3]);
    assert_eq!(
        &b[b.len() - 8..],
        b"QZFOOTR1",
        "archive does not end with QZFOOTR1 footer magic"
    );

    // --- decompress to a prefix and compare both mates ---
    let status = Command::new(QZ)
        .env("QZ_NO_BANNER", "1")
        .args(["decompress", "-i"])
        .arg(&archive)
        .args(["-o"])
        .arg(&prefix)
        .args(["-t", "4"])
        .status()
        .expect("spawn qz decompress (paired)");
    assert!(
        status.success(),
        "qz paired decompress failed with exit code {:?}",
        status.code()
    );

    let out1 = dir.path().join("restored_R1.fastq");
    let out2 = dir.path().join("restored_R2.fastq");
    assert_eq!(
        fs::read(&out1).expect("read restored R1"),
        s1.as_bytes(),
        "R1 mismatch after paired stdout roundtrip"
    );
    assert_eq!(
        fs::read(&out2).expect("read restored R2"),
        s2.as_bytes(),
        "R2 mismatch after paired stdout roundtrip"
    );
}
