use std::path::Path;

use qz_lib::cli::{CompressConfig, ReferenceOptions};

fn ref_cfg(input: std::path::PathBuf, out: std::path::PathBuf, tmp: std::path::PathBuf, refp: std::path::PathBuf) -> CompressConfig {
    CompressConfig {
        input: vec![input],
        output: out,
        working_dir: tmp,
        threads: 1,
        force: true,
        reference: Some(ReferenceOptions {
            reference: refp,
            reference_index: None,
            reference_fast: false,
            reference_window: 2,
        }),
        ..CompressConfig::default()
    }
}

#[test]
fn prebuilt_and_autobuild_archives_are_identical() {
    let d = tempfile::tempdir().unwrap();
    let rf = write_ref(d.path());
    // 150 bp forward-substring reads (first read length 150 → canonical 150).
    let seq = make_seq(2000, 7);
    let mut fq = String::new();
    for (i, st) in [0usize, 200, 400, 600, 800].iter().enumerate() {
        let r = &seq[*st..st + 150];
        let q = "I".repeat(150);
        fq.push_str(&format!("@r{i}\n{}\n+\n{q}\n", std::str::from_utf8(r).unwrap()));
    }
    let reads = d.path().join("reads.fastq");
    std::fs::write(&reads, &fq).unwrap();

    // A) prebuild via qz-index path, then compress (loads it).
    qz_lib::compression::build_reference_index(&rf, 150, false, 1).unwrap();
    let out_a = d.path().join("a.qz");
    qz_lib::compression::compress(&ref_cfg(reads.clone(), out_a.clone(), d.path().to_path_buf(), rf.clone())).unwrap();
    let bytes_a = std::fs::read(&out_a).unwrap();

    // B) remove the sidecar, compress again (library auto-builds it).
    std::fs::remove_file(qz_lib::compression::reference_index_sidecar_path(&rf, 150, false)).unwrap();
    let out_b = d.path().join("b.qz");
    qz_lib::compression::compress(&ref_cfg(reads, out_b.clone(), d.path().to_path_buf(), rf)).unwrap();
    let bytes_b = std::fs::read(&out_b).unwrap();

    assert_eq!(bytes_a, bytes_b, "prebuilt-index archive must equal auto-built archive (threads=1)");
}

#[test]
fn sidecar_path_is_profile_canonical() {
    let p = Path::new("/tmp/ref.fa");
    let s = qz_lib::compression::reference_index_sidecar_path(p, 148, false);
    assert!(
        s.to_str().unwrap().ends_with("ref.fa.qz-r150.sti"),
        "got {}",
        s.display()
    );
    let sf = qz_lib::compression::reference_index_sidecar_path(p, 148, true);
    assert!(
        sf.to_str().unwrap().ends_with("ref.fa.qz-r150-fast.sti"),
        "got {}",
        sf.display()
    );
}

fn make_seq(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        v.push(b"ACGT"[((x >> 33) & 3) as usize]);
    }
    v
}

/// Write a small reference FASTA into `dir`, return its path.
fn write_ref(dir: &Path) -> std::path::PathBuf {
    let refseq = make_seq(2000, 7);
    let rf = dir.join("ref.fa");
    let mut s = String::from(">chr0\n");
    s.push_str(std::str::from_utf8(&refseq).unwrap());
    s.push('\n');
    std::fs::write(&rf, s).unwrap();
    rf
}

#[test]
fn build_reference_index_writes_canonical_sidecar() {
    let d = tempfile::tempdir().unwrap();
    let rf = write_ref(d.path());
    let stats = qz_lib::compression::build_reference_index(&rf, 148, false, 1).unwrap();
    assert!(stats.path.to_str().unwrap().ends_with(".qz-r150.sti"), "{}", stats.path.display());
    assert!(stats.path.exists());
    assert!(stats.references >= 1);
    assert!(stats.randstrobes > 0);
}

#[test]
fn build_reference_index_is_deterministic_single_thread() {
    let d = tempfile::tempdir().unwrap();
    let rf = write_ref(d.path());
    let s1 = qz_lib::compression::build_reference_index(&rf, 150, false, 1).unwrap();
    let bytes1 = std::fs::read(&s1.path).unwrap();
    let s2 = qz_lib::compression::build_reference_index(&rf, 150, false, 1).unwrap();
    let bytes2 = std::fs::read(&s2.path).unwrap();
    assert_eq!(bytes1, bytes2, "index build must be deterministic at threads=1");
}

use qz_lib::compression::ReferenceIndexStatus;

#[test]
fn status_missing_then_ready() {
    let d = tempfile::tempdir().unwrap();
    let rf = write_ref(d.path());
    assert_eq!(
        qz_lib::compression::reference_index_status(&rf, 150, false),
        ReferenceIndexStatus::Missing
    );
    qz_lib::compression::build_reference_index(&rf, 150, false, 1).unwrap();
    assert_eq!(
        qz_lib::compression::reference_index_status(&rf, 150, false),
        ReferenceIndexStatus::Ready
    );
}

#[test]
fn ensure_missing_without_build_errors_with_guidance() {
    let d = tempfile::tempdir().unwrap();
    let rf = write_ref(d.path());
    let err = qz_lib::compression::ensure_reference_index(&rf, 150, false, 1, false).unwrap_err();
    let msg = format!("{err:#}");
    assert!(msg.contains("qz index"), "message must guide the user: {msg}");
    assert!(msg.contains("-r 150"), "message must name the profile: {msg}");
}

#[test]
fn ensure_missing_with_build_creates_it() {
    let d = tempfile::tempdir().unwrap();
    let rf = write_ref(d.path());
    let p = qz_lib::compression::ensure_reference_index(&rf, 150, false, 1, true).unwrap();
    assert!(p.exists());
    assert_eq!(
        qz_lib::compression::reference_index_status(&rf, 150, false),
        ReferenceIndexStatus::Ready
    );
}

#[test]
fn ensure_reuses_index_across_a_profile() {
    let d = tempfile::tempdir().unwrap();
    let rf = write_ref(d.path());
    // Build for 150, then a 143 bp file (same profile) must find it without building.
    qz_lib::compression::build_reference_index(&rf, 150, false, 1).unwrap();
    let p = qz_lib::compression::ensure_reference_index(&rf, 143, false, 1, false).unwrap();
    assert!(p.to_str().unwrap().ends_with(".qz-r150.sti"), "{}", p.display());
}

#[test]
fn ensure_stale_without_build_errors() {
    let d = tempfile::tempdir().unwrap();
    let rf = write_ref(d.path());
    qz_lib::compression::build_reference_index(&rf, 150, false, 1).unwrap();
    // Bump the FASTA mtime to the future so the sidecar is older (stale).
    let f = std::fs::File::options().write(true).open(&rf).unwrap();
    f.set_modified(std::time::SystemTime::now() + std::time::Duration::from_secs(60)).unwrap();
    drop(f);
    assert_eq!(
        qz_lib::compression::reference_index_status(&rf, 150, false),
        ReferenceIndexStatus::Stale
    );
    let err = qz_lib::compression::ensure_reference_index(&rf, 150, false, 1, false).unwrap_err();
    assert!(format!("{err:#}").contains("older than the reference"), "{err:#}");
}

#[test]
fn ensure_corrupt_index_without_build_errors_then_rebuilds() {
    let d = tempfile::tempdir().unwrap();
    let rf = write_ref(d.path());
    // Build a valid index, then clobber it with garbage (keeps a fresh mtime, so
    // only header validation — not the staleness check — can catch it).
    let p = qz_lib::compression::build_reference_index(&rf, 150, false, 1).unwrap().path;
    std::fs::write(&p, b"not a real STI file").unwrap();
    let err = qz_lib::compression::ensure_reference_index(&rf, 150, false, 1, false).unwrap_err();
    assert!(format!("{err:#}").contains("corrupt"), "{err:#}");
    // --build-index rebuilds it; afterwards require-mode is satisfied.
    qz_lib::compression::ensure_reference_index(&rf, 150, false, 1, true).unwrap();
    qz_lib::compression::ensure_reference_index(&rf, 150, false, 1, false).unwrap();
}
