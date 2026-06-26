use std::process::Command;

fn qz() -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_qz"));
    c.env("QZ_NO_BANNER", "1");
    c
}

fn make_ref(dir: &std::path::Path) -> std::path::PathBuf {
    // Deterministic 2 kb ACGT reference.
    let mut x: u64 = 7;
    let mut seq = String::from(">chr0\n");
    for _ in 0..2000 {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        seq.push(b"ACGT"[((x >> 33) & 3) as usize] as char);
    }
    seq.push('\n');
    let p = dir.join("ref.fa");
    std::fs::write(&p, seq).unwrap();
    p
}

fn make_reads(dir: &std::path::Path, refp: &std::path::Path, len: usize) -> std::path::PathBuf {
    // Take forward substrings of the reference as 150 bp reads.
    let refbytes = std::fs::read_to_string(refp).unwrap();
    let seq: String = refbytes.lines().skip(1).collect();
    let seq = seq.as_bytes();
    let mut fq = String::new();
    for (i, st) in [0usize, 200, 400, 600, 800].iter().enumerate() {
        let r = &seq[*st..st + len];
        let q = "I".repeat(len);
        fq.push_str(&format!("@r{i}\n{}\n+\n{q}\n", std::str::from_utf8(r).unwrap()));
    }
    let p = dir.join("reads.fastq");
    std::fs::write(&p, fq).unwrap();
    p
}

#[test]
fn qz_index_writes_canonical_sidecar() {
    let d = tempfile::tempdir().unwrap();
    let rf = make_ref(d.path());
    let out = qz().args(["index"]).arg(&rf).args(["-r", "150", "-t", "1"]).output().unwrap();
    assert!(out.status.success(), "stderr: {}", String::from_utf8_lossy(&out.stderr));
    let sidecar = {
        let mut s = rf.clone().into_os_string();
        s.push(".qz-r150.sti");
        std::path::PathBuf::from(s)
    };
    assert!(sidecar.exists(), "expected {}", sidecar.display());
}

#[test]
fn compress_without_index_errors_with_guidance() {
    let d = tempfile::tempdir().unwrap();
    let rf = make_ref(d.path());
    let reads = make_reads(d.path(), &rf, 150);
    let out = qz()
        .args(["compress", "-i"]).arg(&reads)
        .args(["-o"]).arg(d.path().join("out.qz"))
        .args(["--reference"]).arg(&rf)
        .args(["--numa", "off", "-t", "1", "-f"])
        .output().unwrap();
    assert!(!out.status.success(), "should fail without a prebuilt index");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("qz index"), "actionable message expected: {err}");
}

#[test]
fn compress_with_build_index_succeeds() {
    let d = tempfile::tempdir().unwrap();
    let rf = make_ref(d.path());
    let reads = make_reads(d.path(), &rf, 150);
    let out = qz()
        .args(["compress", "-i"]).arg(&reads)
        .args(["-o"]).arg(d.path().join("out.qz"))
        .args(["--reference"]).arg(&rf)
        .args(["--build-index", "--numa", "off", "-t", "1", "-f"])
        .output().unwrap();
    assert!(out.status.success(), "stderr: {}", String::from_utf8_lossy(&out.stderr));
    assert!(d.path().join("out.qz").exists());
}

#[test]
fn prebuilt_then_compress_succeeds() {
    let d = tempfile::tempdir().unwrap();
    let rf = make_ref(d.path());
    let reads = make_reads(d.path(), &rf, 150);
    // Pre-build via the subcommand.
    let idx = qz().args(["index"]).arg(&rf).args(["-r", "150", "-t", "1"]).output().unwrap();
    assert!(idx.status.success());
    // Compress WITHOUT --build-index now succeeds (index present).
    let out = qz()
        .args(["compress", "-i"]).arg(&reads)
        .args(["-o"]).arg(d.path().join("out.qz"))
        .args(["--reference"]).arg(&rf)
        .args(["--numa", "off", "-t", "1", "-f"])
        .output().unwrap();
    assert!(out.status.success(), "stderr: {}", String::from_utf8_lossy(&out.stderr));
}

/// Truncate the canonical sidecar to its 68-byte header (valid header, no vector
/// payload) — a present, fresh-mtime, but corrupt index.
fn corrupt_sidecar_to_header(rf: &std::path::Path) {
    let sidecar = {
        let mut s = rf.to_path_buf().into_os_string();
        s.push(".qz-r150.sti");
        std::path::PathBuf::from(s)
    };
    let f = std::fs::OpenOptions::new().write(true).open(&sidecar).unwrap();
    f.set_len(68).unwrap(); // 68 bytes = exactly the .sti fixed header
}

#[test]
fn compress_with_corrupt_index_errors_by_default() {
    let d = tempfile::tempdir().unwrap();
    let rf = make_ref(d.path());
    let reads = make_reads(d.path(), &rf, 150);
    let idx = qz().args(["index"]).arg(&rf).args(["-r", "150", "-t", "1"]).output().unwrap();
    assert!(idx.status.success());
    corrupt_sidecar_to_header(&rf);
    // A header-valid but payload-corrupt sidecar must NOT silently rebuild: require
    // by default → error, no archive (regression for the corrupt-payload bypass).
    let out = qz()
        .args(["compress", "-i"]).arg(&reads)
        .args(["-o"]).arg(d.path().join("out.qz"))
        .args(["--reference"]).arg(&rf)
        .args(["--numa", "off", "-t", "1", "-f"])
        .output().unwrap();
    assert!(
        !out.status.success(),
        "corrupt index must error by default; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(!d.path().join("out.qz").exists(), "no archive should be produced for a corrupt index");
}

#[test]
fn compress_with_corrupt_index_and_build_index_rebuilds() {
    let d = tempfile::tempdir().unwrap();
    let rf = make_ref(d.path());
    let reads = make_reads(d.path(), &rf, 150);
    let idx = qz().args(["index"]).arg(&rf).args(["-r", "150", "-t", "1"]).output().unwrap();
    assert!(idx.status.success());
    corrupt_sidecar_to_header(&rf);
    // --build-index self-heals a corrupt sidecar.
    let out = qz()
        .args(["compress", "-i"]).arg(&reads)
        .args(["-o"]).arg(d.path().join("out.qz"))
        .args(["--reference"]).arg(&rf)
        .args(["--build-index", "--numa", "off", "-t", "1", "-f"])
        .output().unwrap();
    assert!(out.status.success(), "stderr: {}", String::from_utf8_lossy(&out.stderr));
    assert!(d.path().join("out.qz").exists());
}
