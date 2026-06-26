use std::process::Command;

fn qz() -> Command {
    Command::new(env!("CARGO_BIN_EXE_qz"))
}

/// Compress a MULTI-CHUNK archive via the CLI by lowering chunk_records through the
/// `--config` JSON path (chunk_records is not a direct CLI flag).
fn write_config(tmp: &std::path::Path, chunk_records: u64) -> std::path::PathBuf {
    let p = tmp.join("cfg.json");
    std::fs::write(&p, format!(r#"{{"chunk_records": {chunk_records}}}"#)).unwrap();
    p
}

fn synth_fastq(path: &std::path::Path, n: usize) {
    use std::io::Write;
    let mut w = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for i in 0..n { writeln!(w, "@r{i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII").unwrap(); }
}

fn compress_single(tmp: &std::path::Path, reads: usize, chunk_records: u64) -> std::path::PathBuf {
    let fq = tmp.join("in.fastq"); synth_fastq(&fq, reads);
    let cfg = write_config(tmp, chunk_records);
    let arc = tmp.join("a.qz");
    assert!(qz().args(["compress", "-i"]).arg(&fq).arg("-o").arg(&arc)
        .arg("-w").arg(tmp).arg("--config").arg(&cfg)
        .env("QZ_NO_BANNER", "1").status().unwrap().success());
    arc
}

fn compress_paired(tmp: &std::path::Path, reads: usize, chunk_records: u64) -> std::path::PathBuf {
    let r1 = tmp.join("r1.fastq"); let r2 = tmp.join("r2.fastq");
    synth_fastq(&r1, reads); synth_fastq(&r2, reads);
    let cfg = write_config(tmp, chunk_records);
    let arc = tmp.join("p.qz");
    assert!(qz().args(["compress", "-i"]).arg(&r1).arg("-i").arg(&r2).arg("-o").arg(&arc)
        .arg("-w").arg(tmp).arg("--config").arg(&cfg)
        .env("QZ_NO_BANNER", "1").status().unwrap().success());
    arc
}

/// Compress a single-end CLUSTER archive (`--cluster`, archive_type 3). Cluster
/// archives are not NUMA-shardable, so `--numa` must route them to in-process decode.
fn compress_cluster(tmp: &std::path::Path, reads: usize) -> std::path::PathBuf {
    let fq = tmp.join("cl.fastq");
    synth_fastq(&fq, reads);
    let arc = tmp.join("cl.qz");
    assert!(
        qz()
            .args(["compress", "--cluster", "-i"]).arg(&fq).arg("-o").arg(&arc)
            .arg("-w").arg(tmp)
            .env("QZ_NO_BANNER", "1").status().unwrap().success()
    );
    arc
}

/// Decompress `arc` to `out` with the given `--numa` arg; returns the Output.
fn decompress_numa(tmp: &std::path::Path, arc: &std::path::Path, out: &std::path::Path, numa: &str) -> std::process::Output {
    qz()
        .args(["decompress", "-i"]).arg(arc).arg("-o").arg(out)
        .arg("-w").arg(tmp).args(["--numa", numa, "-f"])
        .env("QZ_NO_BANNER", "1").output().unwrap()
}

fn no_part_leak(workdir: &std::path::Path) -> bool {
    std::fs::read_dir(workdir).unwrap().filter_map(|e| e.ok())
        .all(|e| !e.file_name().to_string_lossy().contains(".qz_numa"))
}

/// True iff NO `.part` file remains in `dir` (proves the single-end direct-write path,
/// which writes regions into one shared file — not the part-file assembly path).
fn no_dot_part_files(dir: &std::path::Path) -> bool {
    std::fs::read_dir(dir).unwrap().filter_map(|e| e.ok())
        .all(|e| !e.file_name().to_string_lossy().ends_with(".part"))
}

/// Compress a multi-chunk single-end archive with extra CLI flags (e.g. `--no-quality`
/// or `--quality-mode discard`). Returns the archive path.
fn compress_single_flags(
    tmp: &std::path::Path,
    reads: usize,
    chunk_records: u64,
    extra: &[&str],
) -> std::path::PathBuf {
    let fq = tmp.join("in.fastq"); synth_fastq(&fq, reads);
    let cfg = write_config(tmp, chunk_records);
    let arc = tmp.join("a.qz");
    let mut c = qz();
    c.args(["compress", "-i"]).arg(&fq).arg("-o").arg(&arc)
        .arg("-w").arg(tmp).arg("--config").arg(&cfg).args(extra)
        .env("QZ_NO_BANNER", "1");
    assert!(c.status().unwrap().success());
    arc
}

// Relies on the `QZ_NUMA_FORCE_WORKERS` decode hook to deterministically shard on a
// non-NUMA box; that hook is `#[cfg(debug_assertions)]`-gated in driver.rs (release
// must never honor test env vars), so under `--release` the shard never happens and the
// spawn-count assertions below fail. Skipped (not compiled out) in release so the test +
// its helpers keep compiling and the reason is visible in the test output.
#[cfg_attr(not(debug_assertions), ignore = "forced-shard decode hook is debug-only; release compiles it out")]
#[test]
fn forced_shard_single_end_equals_off_and_spawns_workers() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 30_000, 10_000);
    let off = tmp.path().join("off.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&off)
        .arg("-w").arg(tmp.path()).args(["--numa", "off", "-f"])
        .env("QZ_NO_BANNER", "1").status().unwrap().success());
    let sharded = tmp.path().join("sharded.fastq");
    let spawn_log = tmp.path().join("spawned.txt");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&sharded)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2").env("QZ_NUMA_SPAWN_LOG", &spawn_log)
        .status().unwrap().success());
    assert_eq!(std::fs::read(&off).unwrap(), std::fs::read(&sharded).unwrap());
    assert_eq!(std::fs::read_to_string(&spawn_log).unwrap().trim(), "2");
    assert!(no_part_leak(tmp.path()));
}

#[test]
fn forced_shard_single_direct_write_byte_identical() {
    // Single-end + non-gzip ⇒ direct-write: forced-2-worker `--numa auto` output must
    // be byte-identical to `--numa off`, and NO `.part` files may remain (which would
    // mean the part-file assembly path ran instead of the direct path).
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 30_000, 10_000);
    let off = tmp.path().join("off.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&off)
        .arg("-w").arg(tmp.path()).args(["--numa", "off", "-f"])
        .env("QZ_NO_BANNER", "1").status().unwrap().success());
    let direct = tmp.path().join("direct.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&direct)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2")
        .status().unwrap().success());
    assert_eq!(std::fs::read(&off).unwrap(), std::fs::read(&direct).unwrap());
    assert!(no_dot_part_files(tmp.path()), "direct-write must not create .part files");
    assert!(no_part_leak(tmp.path()));
}

#[test]
fn forced_shard_single_no_quality_and_discard_byte_identical() {
    // Direct-write must be byte-identical to `--numa off` for archives compressed
    // without quality data (both --no-quality and --quality-mode discard).
    for extra in [&["--no-quality"][..], &["--quality-mode", "discard"][..]] {
        let tmp = tempfile::tempdir().unwrap();
        let arc = compress_single_flags(tmp.path(), 30_000, 10_000, extra);
        let off = tmp.path().join("off.fastq");
        assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&off)
            .arg("-w").arg(tmp.path()).args(["--numa", "off", "-f"])
            .env("QZ_NO_BANNER", "1").status().unwrap().success(),
            "off decode failed for {extra:?}");
        let direct = tmp.path().join("direct.fastq");
        assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&direct)
            .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
            .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2")
            .status().unwrap().success(), "direct decode failed for {extra:?}");
        assert_eq!(std::fs::read(&off).unwrap(), std::fs::read(&direct).unwrap(),
            "byte mismatch for {extra:?}");
        assert!(no_dot_part_files(tmp.path()), "direct-write must not create .part files ({extra:?})");
    }
}

#[test]
fn forced_shard_fallback_on_region_mismatch() {
    // QZ_NUMA_FAKE_REGION perturbs a worker's region len → its base/len verification
    // against the table fails → worker exit 4 → auto cleans the temp + decodes
    // in-process → output equals `--numa off`.
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 30_000, 10_000);
    let off = tmp.path().join("off.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&off)
        .arg("-w").arg(tmp.path()).args(["--numa", "off", "-f"])
        .env("QZ_NO_BANNER", "1").status().unwrap().success());
    let out = tmp.path().join("out.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2").env("QZ_NUMA_FAKE_REGION", "1")
        .status().unwrap().success(), "auto must fall back on a region/integrity failure");
    assert_eq!(std::fs::read(&out).unwrap(), std::fs::read(&off).unwrap());
    assert!(no_dot_part_files(tmp.path()));
    assert!(no_part_leak(tmp.path()));
}

#[test]
fn forced_shard_fallback_on_table_underestimate() {
    // Compress with QZ_NUMA_FAKE_UNDERSIZE (the last chunk's decoded size is undersized
    // in the table). Direct-write decode then overruns the bounded region → worker
    // exit 4 → in-process fallback (which ignores the table) → output equals a normal
    // `--numa off` decode of the SAME (undersized-table) archive.
    let tmp = tempfile::tempdir().unwrap();
    let fq = tmp.path().join("in.fastq"); synth_fastq(&fq, 30_000);
    let cfg = write_config(tmp.path(), 10_000);
    let arc = tmp.path().join("a.qz");
    assert!(qz().args(["compress", "-i"]).arg(&fq).arg("-o").arg(&arc)
        .arg("-w").arg(tmp.path()).arg("--config").arg(&cfg)
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FAKE_UNDERSIZE", "1")
        .status().unwrap().success());
    let off = tmp.path().join("off.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&off)
        .arg("-w").arg(tmp.path()).args(["--numa", "off", "-f"])
        .env("QZ_NO_BANNER", "1").status().unwrap().success());
    let out = tmp.path().join("out.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2")
        .status().unwrap().success(), "auto must fall back on a table-underestimate overrun");
    assert_eq!(std::fs::read(&out).unwrap(), std::fs::read(&off).unwrap());
    assert!(no_dot_part_files(tmp.path()));
    assert!(no_part_leak(tmp.path()));
}

#[test]
fn forced_shard_paired_equals_off() {
    // Paired + non-gzip + a present ChunkDecodedSizes(n_mates=2) table ⇒ direct-write:
    // forced-2-worker `--numa auto` writes both mates straight into their pre-sized
    // shared outputs (no part files, no concat), and the R1/R2 outputs must be
    // byte-identical to `--numa off`.
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_paired(tmp.path(), 30_000, 10_000);
    for (m, hook) in [("off", false), ("auto", true)] {
        let prefix = tmp.path().join(m);
        let mut c = qz();
        c.args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&prefix)
            .arg("-w").arg(tmp.path()).args(["--numa", m, "-f"]).env("QZ_NO_BANNER", "1");
        if hook { c.env("QZ_NUMA_FORCE_WORKERS", "2"); }
        assert!(c.status().unwrap().success());
    }
    for sfx in ["_R1.fastq", "_R2.fastq"] {
        let a = path_suffix(&tmp.path().join("off"), sfx);
        let b = path_suffix(&tmp.path().join("auto"), sfx);
        assert_eq!(std::fs::read(&a).unwrap(), std::fs::read(&b).unwrap(), "mate {sfx}");
    }
    assert!(no_dot_part_files(tmp.path()), "paired direct-write must not create .part files");
    assert!(no_part_leak(tmp.path()));
}

#[test]
fn forced_shard_paired_fallback_on_region_mismatch() {
    // Paired direct-write fault path: QZ_NUMA_FAKE_REGION perturbs each worker's region
    // len → its base/len verification against the per-mate table fails → worker exit 4 →
    // auto cleans BOTH pre-sized temps + decodes in-process → R1/R2 equal `--numa off`.
    // Exercises the new 2-region worker args + the paired two-file cleanup-on-failure.
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_paired(tmp.path(), 30_000, 10_000);
    let off = tmp.path().join("off");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&off)
        .arg("-w").arg(tmp.path()).args(["--numa", "off", "-f"])
        .env("QZ_NO_BANNER", "1").status().unwrap().success());
    let out = tmp.path().join("out");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2").env("QZ_NUMA_FAKE_REGION", "1")
        .status().unwrap().success(), "auto must fall back on a paired region/integrity failure");
    for sfx in ["_R1.fastq", "_R2.fastq"] {
        let a = path_suffix(&off, sfx);
        let b = path_suffix(&out, sfx);
        assert_eq!(std::fs::read(&a).unwrap(), std::fs::read(&b).unwrap(), "mate {sfx}");
    }
    assert!(no_dot_part_files(tmp.path()));
    assert!(no_part_leak(tmp.path()));
}

#[test]
fn forced_shard_single_end_gzip_decompresses_identically() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 30_000, 10_000);
    let off = tmp.path().join("off.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&off)
        .arg("-w").arg(tmp.path()).args(["--numa", "off", "-f"])
        .env("QZ_NO_BANNER", "1").status().unwrap().success());
    let gz = tmp.path().join("sharded.fastq.gz");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&gz)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "--gzipped", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2")
        .status().unwrap().success());
    let raw = std::fs::read(&gz).unwrap();
    let mut dec = flate2::read::MultiGzDecoder::new(&raw[..]);
    let mut out = Vec::new();
    std::io::Read::read_to_end(&mut dec, &mut out).unwrap();
    assert_eq!(out, std::fs::read(&off).unwrap());
}

// Debug-only: drives the `QZ_NUMA_FORCE_WORKERS` + `QZ_NUMA_SPAWN_FAIL` decode hooks
// (both `#[cfg(debug_assertions)]`-gated in driver.rs). Inert under `--release`, so the
// forced shard + injected worker failure never occur — skipped in release.
#[cfg_attr(not(debug_assertions), ignore = "forced-shard decode hook is debug-only; release compiles it out")]
#[test]
fn worker_failure_kills_survivor_and_aborts_cleanly() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 30_000, 10_000);
    let out = tmp.path().join("out.fastq");
    let t0 = std::time::Instant::now();
    let st = qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2")
        .env("QZ_NUMA_FAIL", "0").env("QZ_NUMA_SLOW", "1")
        .status().unwrap();
    let elapsed = t0.elapsed();
    assert!(!st.success(), "worker failure must abort");
    assert!(elapsed.as_secs() < 15, "driver must kill the sleeping survivor (took {elapsed:?})");
    assert!(!out.exists(), "no partial output published");
    assert!(no_part_leak(tmp.path()));
}

// Debug-only: drives the `QZ_NUMA_FORCE_WORKERS` + `QZ_NUMA_SPAWN_FAIL` decode hooks
// (both `#[cfg(debug_assertions)]`-gated in driver.rs). Inert under `--release` — skipped.
#[cfg_attr(not(debug_assertions), ignore = "forced-shard decode hook is debug-only; release compiles it out")]
#[test]
fn spawn_failure_aborts_cleanly() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 30_000, 10_000);
    let out = tmp.path().join("out.fastq");
    let st = qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2").env("QZ_NUMA_SPAWN_FAIL", "1")
        .status().unwrap();
    assert!(!st.success(), "spawn failure must abort");
    assert!(!out.exists());
    assert!(no_part_leak(tmp.path()));
}

#[test]
fn pin_failure_auto_falls_back() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 30_000, 10_000);
    let out = tmp.path().join("out.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2").env("QZ_NUMA_PINFAIL", "all")
        .status().unwrap().success(), "auto must fall back on pin failure");
    let off = tmp.path().join("off.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&off)
        .arg("-w").arg(tmp.path()).args(["--numa", "off", "-f"]).env("QZ_NO_BANNER", "1")
        .status().unwrap().success());
    assert_eq!(std::fs::read(&out).unwrap(), std::fs::read(&off).unwrap());
    assert!(no_part_leak(tmp.path()));
}

#[test]
fn numa_one_runs_in_process_anywhere() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 5_000, 2_500);
    let out = tmp.path().join("out.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "1", "-f"])
        .env("QZ_NO_BANNER", "1").status().unwrap().success(), "--numa 1 must run in-process");
}

#[test]
fn paired_gzip_rejected() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_paired(tmp.path(), 5_000, 2_500);
    let st = qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(tmp.path().join("p"))
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "--gzipped", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2")
        .status().unwrap();
    assert!(!st.success(), "paired + gzip must be rejected");
}

#[test]
fn forced_shard_respects_force_flag() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 30_000, 10_000);
    let out = tmp.path().join("out.fastq");
    std::fs::write(&out, b"OLD").unwrap();
    let st = qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2").status().unwrap();
    assert!(!st.success());
    assert_eq!(std::fs::read(&out).unwrap(), b"OLD");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2").status().unwrap().success());
    assert_ne!(std::fs::read(&out).unwrap(), b"OLD");
}

#[test]
fn forced_shard_refuses_existing_output_before_decoding() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 30_000, 10_000);
    let out = tmp.path().join("out.fastq");
    std::fs::write(&out, b"OLD").unwrap();
    let t0 = std::time::Instant::now();
    let st = qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2")
        .status().unwrap();
    let dt = t0.elapsed();
    assert!(!st.success(), "must refuse existing output without --force");
    assert!(dt.as_secs() < 3, "refusal must be pre-decode/fast (took {dt:?})");
    assert_eq!(std::fs::read(&out).unwrap(), b"OLD", "existing output untouched");
    assert!(no_part_leak(tmp.path()));
}

#[test]
fn stdout_output_fixed_errors_auto_falls_back() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_single(tmp.path(), 5_000, 2_500);
    let st = qz().args(["decompress", "-i"]).arg(&arc).args(["-o", "-", "--numa", "2"])
        .arg("-w").arg(tmp.path()).env("QZ_NO_BANNER", "1").status().unwrap();
    assert!(!st.success(), "fixed + stdout must error");
    let outp = qz().args(["decompress", "-i"]).arg(&arc).args(["-o", "-", "--numa", "auto"])
        .arg("-w").arg(tmp.path()).env("QZ_NO_BANNER", "1").output().unwrap();
    assert!(outp.status.success());
    assert!(outp.stdout.starts_with(b"@r0\n"));
}

fn path_suffix(prefix: &std::path::Path, suffix: &str) -> std::path::PathBuf {
    let mut s = prefix.as_os_str().to_owned(); s.push(suffix); std::path::PathBuf::from(s)
}

/// Build a deterministic low-repeat ACGT reference sequence.
fn make_seq_ref(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        v.push(b"ACGT"[((x >> 33) & 3) as usize]);
    }
    v
}

/// Build a reference FASTA + paired R1/R2 reads that mostly map against it.
/// R1 = forward 120bp slice from reference; R2 = reverse-complement of a
/// downstream slice. ~2% are off-reference junk (literal fallback path); a
/// sparse subset carries a single substitution (edit path). Returns (fa, arc).
fn build_reference_archive(
    tmp: &std::path::Path,
    reads: usize,
    chunk_records: u64,
) -> (std::path::PathBuf, std::path::PathBuf) {
    let ref_len = 50_000usize;
    let refseq = make_seq_ref(ref_len, 7);

    let fa = tmp.join("ref.fa");
    {
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        std::fs::write(&fa, &s).unwrap();
    }

    fn revcomp(seq: &[u8]) -> Vec<u8> {
        seq.iter().rev().map(|&b| match b {
            b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', other => other,
        }).collect()
    }

    let rlen = 120usize;
    let max_start = ref_len - rlen - 200;
    let mut rng = 0x1234_5678_9abc_def0u64;
    let mut next = || {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (rng >> 33) as usize
    };

    let q = "I".repeat(rlen);
    let mut s1 = String::with_capacity(reads * (rlen + 50));
    let mut s2 = String::with_capacity(reads * (rlen + 50));

    for i in 0..reads {
        let unmappable = i % 50 == 0;
        let (r1bytes, r2bytes) = if unmappable {
            (make_seq_ref(rlen, 0x9000_0000 + i as u64), make_seq_ref(rlen, 0x7000_0000 + i as u64))
        } else {
            let st1 = next() % (max_start + 1);
            let mut r1 = refseq[st1..st1 + rlen].to_vec();
            let st2 = st1 + 150;
            let r2_fwd = refseq[st2..st2 + rlen].to_vec();
            let mut r2 = revcomp(&r2_fwd);
            if i % 17 == 0 { let p = next() % rlen; r1[p] = if r1[p] == b'A' { b'C' } else { b'A' }; }
            if i % 23 == 0 { let p = next() % rlen; r2[p] = if r2[p] == b'G' { b'T' } else { b'G' }; }
            (r1, r2)
        };
        s1.push_str(&format!("@frag_{i}/1\n{}\n+\n{q}\n", std::str::from_utf8(&r1bytes).unwrap()));
        s2.push_str(&format!("@frag_{i}/2\n{}\n+\n{q}\n", std::str::from_utf8(&r2bytes).unwrap()));
    }

    let r1p = tmp.join("r1.fastq");
    let r2p = tmp.join("r2.fastq");
    std::fs::write(&r1p, &s1).unwrap();
    std::fs::write(&r2p, &s2).unwrap();

    let cfg = write_config(tmp, chunk_records);
    let arc = tmp.join("ref.qz");
    assert!(
        qz().args(["compress", "-i"]).arg(&r1p).arg("-i").arg(&r2p)
            .arg("--reference").arg(&fa).arg("--build-index")
            .arg("-o").arg(&arc).arg("-w").arg(tmp).arg("--config").arg(&cfg)
            .env("QZ_NO_BANNER", "1").status().unwrap().success(),
        "reference compress failed"
    );
    (fa, arc)
}

/// Build a reference FASTA + a SINGLE R1 read file that mostly maps against it, then
/// compress with one `-i` + `--reference` → a single-end reference archive (archive_type
/// 4). Mirrors `build_reference_archive` but emits only R1 and compresses single-input.
/// Returns (r1_fastq, fa, archive).
fn build_se_reference_archive(
    tmp: &std::path::Path,
    reads: usize,
    chunk_records: u64,
) -> (std::path::PathBuf, std::path::PathBuf, std::path::PathBuf) {
    let ref_len = 50_000usize;
    let refseq = make_seq_ref(ref_len, 7);

    let fa = tmp.join("ref.fa");
    {
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        std::fs::write(&fa, &s).unwrap();
    }

    let rlen = 120usize;
    let max_start = ref_len - rlen - 200;
    let mut rng = 0x1234_5678_9abc_def0u64;
    let mut next = || {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (rng >> 33) as usize
    };

    let q = "I".repeat(rlen);
    let mut s1 = String::with_capacity(reads * (rlen + 50));
    for i in 0..reads {
        let unmappable = i % 50 == 0;
        let r1bytes = if unmappable {
            make_seq_ref(rlen, 0x9000_0000 + i as u64)
        } else {
            let st1 = next() % (max_start + 1);
            let mut r1 = refseq[st1..st1 + rlen].to_vec();
            if i % 17 == 0 { let p = next() % rlen; r1[p] = if r1[p] == b'A' { b'C' } else { b'A' }; }
            r1
        };
        s1.push_str(&format!("@frag_{i}/1\n{}\n+\n{q}\n", std::str::from_utf8(&r1bytes).unwrap()));
    }

    let r1p = tmp.join("r1.fastq");
    std::fs::write(&r1p, &s1).unwrap();

    let cfg = write_config(tmp, chunk_records);
    let arc = tmp.join("refse.qz");
    assert!(
        qz().args(["compress", "-i"]).arg(&r1p)
            .arg("--reference").arg(&fa).arg("--build-index")
            .arg("-o").arg(&arc).arg("-w").arg(tmp).arg("--config").arg(&cfg)
            .env("QZ_NO_BANNER", "1").status().unwrap().success(),
        "single-end reference compress failed"
    );
    (r1p, fa, arc)
}

/// `qz decompress --numa 2` (forced-shard hook, unpinned) on a single-end reference
/// (archive_type 4) archive must produce ONE output byte-identical to an in-process
/// decode (single-output assembly path — bare `-o` path, no `_R1`/`_R2` suffix).
///
/// Debug-only: the `QZ_NUMA_FORCE_WORKERS` decode hook is `#[cfg(debug_assertions)]`-gated
/// in driver.rs, so under `--release` the type-4 shard never happens and the spawn-log
/// assertion fails — skipped in release.
#[cfg_attr(not(debug_assertions), ignore = "forced-shard decode hook is debug-only; release compiles it out")]
#[test]
fn numa_forced_shard_single_end_reference_matches_in_process() {
    let tmp = tempfile::tempdir().unwrap();
    let (_r1, _fa, archive) = build_se_reference_archive(tmp.path(), 30_000, 10_000); // 3 chunks

    // In-process reference (single-output: writes ONE file at the verbatim -o path).
    let ref_out = tmp.path().join("ref.fastq");
    assert!(qz().args(["decompress", "-i"]).arg(&archive)
        .arg("-o").arg(&ref_out).arg("-w").arg(tmp.path())
        .args(["--numa", "off", "-f"])
        .env("QZ_NO_BANNER", "1").status().unwrap().success());

    // Forced 2-way shard, children UNPINNED via the debug-only hook.
    let shard_out = tmp.path().join("shard.fastq");
    let spawn_log = tmp.path().join("spawned.txt");
    assert!(qz().args(["decompress", "-i"]).arg(&archive)
        .arg("-o").arg(&shard_out).arg("-w").arg(tmp.path())
        .args(["--numa", "2", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2").env("QZ_NUMA_SPAWN_LOG", &spawn_log)
        .status().unwrap().success());

    assert_eq!(std::fs::read_to_string(&spawn_log).unwrap().trim(), "2", "type-4 must shard");
    assert_eq!(
        std::fs::read(&shard_out).unwrap(),
        std::fs::read(&ref_out).unwrap(),
        "type-4 sharded NUMA decode != in-process decode"
    );
    // type-4 reference now carries a 1-mate ChunkDecodedSizes table ⇒ direct-write
    // (workers seek into the pre-sized output; no part files / concat).
    assert!(no_dot_part_files(tmp.path()), "type-4 reference direct-write must not create .part files");
    assert!(no_part_leak(tmp.path()));
}

// Debug-only: relies on the `QZ_NUMA_FORCE_WORKERS` decode hook
// (`#[cfg(debug_assertions)]`-gated in driver.rs). Inert under `--release`, so the shard
// + spawn-log assertion never fire — skipped in release.
#[cfg_attr(not(debug_assertions), ignore = "forced-shard decode hook is debug-only; release compiles it out")]
#[test]
fn forced_shard_reference_equals_off_and_spawns_workers() {
    let tmp = tempfile::tempdir().unwrap();
    let (_fa, arc) = build_reference_archive(tmp.path(), 30_000, 10_000);

    let off = tmp.path().join("off");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&off)
        .arg("-w").arg(tmp.path()).args(["--numa", "off", "-f"])
        .env("QZ_NO_BANNER", "1").status().unwrap().success());

    let auto = tmp.path().join("auto");
    let spawn_log = tmp.path().join("spawned.txt");
    assert!(qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&auto)
        .arg("-w").arg(tmp.path()).args(["--numa", "auto", "-f"])
        .env("QZ_NO_BANNER", "1").env("QZ_NUMA_FORCE_WORKERS", "2").env("QZ_NUMA_SPAWN_LOG", &spawn_log)
        .status().unwrap().success());

    assert_eq!(std::fs::read_to_string(&spawn_log).unwrap().trim(), "2", "reference must shard");
    for sfx in ["_R1.fastq", "_R2.fastq"] {
        let a = path_suffix(&off, sfx);
        let b = path_suffix(&auto, sfx);
        assert_eq!(std::fs::read(&a).unwrap(), std::fs::read(&b).unwrap(), "mate {sfx}");
    }
    // type-2 reference now carries a 2-mate ChunkDecodedSizes table ⇒ direct-write.
    assert!(no_dot_part_files(tmp.path()), "type-2 reference direct-write must not create .part files");
    assert!(no_part_leak(tmp.path()));
}

#[allow(dead_code)]
pub(crate) fn is_numa_box() -> bool {
    std::fs::read_to_string("/sys/devices/system/node/online")
        .map(|s| s.contains('-') || s.contains(','))
        .unwrap_or(false)
}

#[test]
fn numa_off_decodes_normally() {
    let tmp = tempfile::tempdir().unwrap();
    let fq = tmp.path().join("in.fastq");
    std::fs::write(&fq, "@r0\nACGT\n+\nIIII\n").unwrap();
    let arc = tmp.path().join("a.qz");
    assert!(qz().args(["compress", "-i"]).arg(&fq).arg("-o").arg(&arc).arg("-w").arg(tmp.path())
        .env("QZ_NO_BANNER", "1").status().unwrap().success());
    let out = tmp.path().join("out.fastq");
    let st = qz().args(["decompress", "-i"]).arg(&arc).arg("-o").arg(&out)
        .arg("-w").arg(tmp.path()).args(["--numa", "off", "-f"])
        .env("QZ_NO_BANNER", "1").status().unwrap();
    assert!(st.success());
    assert_eq!(std::fs::read_to_string(&out).unwrap(), "@r0\nACGT\n+\nIIII\n");
}

#[test]
fn numa_rejects_bad_value() {
    let st = qz().args(["decompress", "-i", "x", "-o", "y", "--numa", "banana"])
        .env("QZ_NO_BANNER", "1").status().unwrap();
    assert!(!st.success());
}

/// `--numa auto` on a cluster archive (archive_type 3, not shardable) must decode
/// in-process and produce output byte-identical to `--numa off`. Portable: `auto`
/// always falls back gracefully whether or not the box has NUMA hardware.
#[test]
fn numa_auto_cluster_decodes_in_process() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_cluster(tmp.path(), 3_000);
    let off = tmp.path().join("off.fastq");
    let auto = tmp.path().join("auto.fastq");
    assert!(decompress_numa(tmp.path(), &arc, &off, "off").status.success());
    let o = decompress_numa(tmp.path(), &arc, &auto, "auto");
    assert!(o.status.success(), "--numa auto on cluster must succeed");
    assert_eq!(
        std::fs::read(&off).unwrap(),
        std::fs::read(&auto).unwrap(),
        "cluster decode must be identical under --numa off and auto"
    );
}

/// `--numa N` (Fixed) on a cluster archive must NEVER leak the internal
/// "unknown archive_type 3" error (the regression this fix closes). On NUMA hardware
/// it decodes in-process (output == `--numa off`); on a non-NUMA box it gives a clean
/// NUMA-usage error. Either way: no internal type leak, and any success is correct.
#[test]
fn numa_fixed_cluster_no_internal_crash() {
    let tmp = tempfile::tempdir().unwrap();
    let arc = compress_cluster(tmp.path(), 3_000);
    let off = tmp.path().join("off.fastq");
    let fixed = tmp.path().join("fixed.fastq");
    assert!(decompress_numa(tmp.path(), &arc, &off, "off").status.success());
    let o = decompress_numa(tmp.path(), &arc, &fixed, "2");
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(
        !stderr.contains("unknown archive_type"),
        "--numa N must not leak the internal type error: {stderr}"
    );
    if o.status.success() {
        // Decoded in-process (NUMA hardware path) — must match the off decode.
        assert_eq!(
            std::fs::read(&off).unwrap(),
            std::fs::read(&fixed).unwrap(),
            "in-process fallback decode must match --numa off"
        );
    }
}
