//! NUMA-compress orchestrator E2E via QZ_NUMA_FORCE_WORKERS (debug-only hook).

use std::path::{Path, PathBuf};
use std::process::Command;

fn qz_bin() -> PathBuf { PathBuf::from(env!("CARGO_BIN_EXE_qz")) }

fn write_fixed_fastq(path: &Path, n: usize, len: usize) {
    let seq = "A".repeat(len);
    let q = "I".repeat(len);
    let mut s = String::new();
    for i in 0..n { s.push_str(&format!("@r{i}\n{seq}\n+\n{q}\n")); }
    std::fs::write(path, s).unwrap();
}

fn run(args: &[&str], env: &[(&str, &str)]) -> std::process::Output {
    let mut c = Command::new(qz_bin());
    c.args(args).env("QZ_NO_BANNER", "1");
    for (k, v) in env { c.env(k, v); }
    c.output().unwrap()
}

fn compress(input: &Path, out: &Path, wd: &Path, numa: &str, env: &[(&str, &str)]) -> std::process::Output {
    run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
          "-w", wd.to_str().unwrap(), "-t", "4", "--numa", numa, "-f"], env)
}

fn decompress_text(arc: &Path, wd: &Path) -> String {
    let out = wd.join("dec.fastq");
    let o = run(&["decompress", "-i", arc.to_str().unwrap(), "-o", out.to_str().unwrap(),
                  "-w", wd.to_str().unwrap(), "-f"], &[]);
    assert!(o.status.success(), "decompress failed: {}", String::from_utf8_lossy(&o.stderr));
    std::fs::read_to_string(&out).unwrap()
}

#[test]
fn forced_shard_single_end_roundtrips_byte_identical_to_in_process() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fastq");
    write_fixed_fastq(&input, 400, 60); // many records → splittable into 3
    let expected = std::fs::read_to_string(&input).unwrap();

    let shard_out = d.join("shard.qz");
    let o = compress(&input, &shard_out, d, "auto", &[("QZ_NUMA_FORCE_WORKERS", "3")]);
    assert!(o.status.success(), "forced shard compress failed: {}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(decompress_text(&shard_out, d), expected, "sharded decoded output must equal the input");

    let inproc_out = d.join("inproc.qz");
    let o2 = compress(&input, &inproc_out, d, "off", &[]);
    assert!(o2.status.success());
    let shard_bytes = std::fs::read(&shard_out).unwrap();
    let inproc_bytes = std::fs::read(&inproc_out).unwrap();
    // The v5 front header's on-disk length is a u32 LE at bytes 4..8; it is 20 or 28
    // bytes. Compare EXACTLY that many bytes — a fixed [..32] would bleed into the
    // first block frame, which differs by design (sharded re-chunks).
    let hlen = u32::from_le_bytes(shard_bytes[4..8].try_into().unwrap()) as usize;
    assert_eq!(hlen, u32::from_le_bytes(inproc_bytes[4..8].try_into().unwrap()) as usize,
        "sharded + in-process header sizes must match");
    assert!(hlen == 20 || hlen == 28, "v5 header is 20 or 28 bytes, got {hlen}");
    assert_eq!(&shard_bytes[..hlen], &inproc_bytes[..hlen],
        "sharded + in-process front headers must be byte-identical");
    // Whole-archive byte-identity is FALSE by design (re-chunking) — do NOT assert it.
}

fn write_fixed_paired(r1: &Path, r2: &Path, n: usize, len: usize) {
    let (sq1, sq2) = ("A".repeat(len), "C".repeat(len));
    let (q1, q2) = ("I".repeat(len), "F".repeat(len));
    let mut s1 = String::new();
    let mut s2 = String::new();
    for i in 0..n {
        s1.push_str(&format!("@r{i}/1\n{sq1}\n+\n{q1}\n"));
        s2.push_str(&format!("@r{i}/2\n{sq2}\n+\n{q2}\n"));
    }
    std::fs::write(r1, s1).unwrap();
    std::fs::write(r2, s2).unwrap();
}

#[test]
fn forced_shard_paired_roundtrips_lossless() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let r1 = d.join("r1.fastq");
    let r2 = d.join("r2.fastq");
    write_fixed_paired(&r1, &r2, 400, 60); // splittable into 3
    let exp1 = std::fs::read_to_string(&r1).unwrap();
    let exp2 = std::fs::read_to_string(&r2).unwrap();

    let out = d.join("paired.qz");
    let o = run(
        &["compress", "-i", r1.to_str().unwrap(), "-i", r2.to_str().unwrap(),
          "-o", out.to_str().unwrap(), "-w", d.to_str().unwrap(), "-t", "4", "--numa", "auto", "-f"],
        &[("QZ_NUMA_FORCE_WORKERS", "3")],
    );
    assert!(o.status.success(), "forced paired shard compress failed: {}", String::from_utf8_lossy(&o.stderr));

    // Paired decompress expands `-o <prefix>` → {prefix}_R1.fastq / _R2.fastq.
    let dec = d.join("dec");
    let o2 = run(
        &["decompress", "-i", out.to_str().unwrap(), "-o", dec.to_str().unwrap(), "-w", d.to_str().unwrap(), "-f"],
        &[],
    );
    assert!(o2.status.success(), "paired decompress failed: {}", String::from_utf8_lossy(&o2.stderr));
    assert_eq!(std::fs::read_to_string(d.join("dec_R1.fastq")).unwrap(), exp1, "sharded paired R1 must roundtrip");
    assert_eq!(std::fs::read_to_string(d.join("dec_R2.fastq")).unwrap(), exp2, "sharded paired R2 must roundtrip");

    // Front-header byte-identity with in-process paired; archive_type byte == 1.
    let inproc = d.join("inproc.qz");
    let o3 = run(
        &["compress", "-i", r1.to_str().unwrap(), "-i", r2.to_str().unwrap(),
          "-o", inproc.to_str().unwrap(), "-w", d.to_str().unwrap(), "-t", "4", "--numa", "off", "-f"],
        &[],
    );
    assert!(o3.status.success());
    let sb = std::fs::read(&out).unwrap();
    let ib = std::fs::read(&inproc).unwrap();
    let hlen = u32::from_le_bytes(sb[4..8].try_into().unwrap()) as usize;
    assert_eq!(hlen, u32::from_le_bytes(ib[4..8].try_into().unwrap()) as usize, "paired header sizes match");
    assert_eq!(&sb[..hlen], &ib[..hlen], "paired sharded + in-process front headers byte-identical");
    assert_eq!(sb[3], 1, "archive_type byte must be 1 (paired)");
}

#[test]
fn forced_shard_paired_mate_mismatch_errors_no_output() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let r1 = d.join("r1.fastq");
    let r2 = d.join("r2.fastq");
    // R1 has 400 records, R2 has 399 → lockstep split must reject (cannot pair-align).
    let (seq, q) = ("A".repeat(60), "I".repeat(60));
    let mut s1 = String::new();
    let mut s2 = String::new();
    for i in 0..400 { s1.push_str(&format!("@r{i}/1\n{seq}\n+\n{q}\n")); }
    for i in 0..399 { s2.push_str(&format!("@r{i}/2\n{seq}\n+\n{q}\n")); }
    std::fs::write(&r1, s1).unwrap();
    std::fs::write(&r2, s2).unwrap();

    let out = d.join("p.qz");
    // Fixed --numa N + mismatch → clear error, never a corrupt archive.
    let o = run(
        &["compress", "-i", r1.to_str().unwrap(), "-i", r2.to_str().unwrap(),
          "-o", out.to_str().unwrap(), "-w", d.to_str().unwrap(), "-t", "4", "--numa", "2", "-f"],
        &[("QZ_NUMA_FORCE_WORKERS", "2")],
    );
    assert!(!o.status.success(), "mate-count mismatch under fixed --numa must error");
    assert!(!out.exists(), "no archive must be left on mate mismatch");
}

/// Write a partial AdvancedOptions JSON (relies on `#[serde(default)]`) and return
/// the path, for tests that need a non-default `chunk_records`.
fn write_config(d: &Path, chunk_records: usize) -> PathBuf {
    let p = d.join("cfg.json");
    std::fs::write(&p, format!("{{\"chunk_records\": {chunk_records}}}")).unwrap();
    p
}

#[test]
fn const_length_violation_in_shard_errors_not_corrupts() {
    // chunk_records=8 → BOTH the parent's plan resolution AND the in-process compress
    // resolve const-60 from the uniform 60bp head. A shard that starts on a 61bp
    // record fails the worker's forced first-chunk const-validation (C1 fix) → exit 1.
    // Under auto the in-process fallback hits the SAME genuine const violation and ALSO
    // errors — so qz REJECTS the input rather than emitting a corrupt archive.
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let cfg = write_config(d, 8);
    let input = d.join("var.fastq");
    let mut s = String::new();
    for i in 0..200 { s.push_str(&format!("@h{i}\n{}\n+\n{}\n", "A".repeat(60), "I".repeat(60))); }
    for i in 0..200 { s.push_str(&format!("@t{i}\n{}\n+\n{}\n", "C".repeat(61), "I".repeat(61))); }
    std::fs::write(&input, &s).unwrap();

    for numa in ["auto", "4"] {
        let out = d.join(format!("var_{numa}.qz"));
        let o = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                      "-w", d.to_str().unwrap(), "-t", "4", "--config", cfg.to_str().unwrap(),
                      "--numa", numa, "-f"],
                   &[("QZ_NUMA_FORCE_WORKERS", "4")]);
        assert!(!o.status.success(), "const-length violation under --numa {numa} must error, not corrupt");
        assert!(!out.exists(), "no archive must be left at the final path under --numa {numa}");
    }
}

#[test]
fn variable_length_input_shards_correctly_with_default_chunk() {
    // Default chunk_records (2.5M) → the parent's first chunk spans BOTH lengths →
    // const_seq_len=0 → variable framing → no worker validation failure → shards OK.
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("var.fastq");
    let mut s = String::new();
    for i in 0..200 { s.push_str(&format!("@h{i}\n{}\n+\n{}\n", "A".repeat(60), "I".repeat(60))); }
    for i in 0..200 { s.push_str(&format!("@t{i}\n{}\n+\n{}\n", "C".repeat(61), "I".repeat(61))); }
    std::fs::write(&input, &s).unwrap();
    let expected = s.clone();
    let out = d.join("var.qz");
    let o = compress(&input, &out, d, "auto", &[("QZ_NUMA_FORCE_WORKERS", "4")]);
    assert!(o.status.success(), "variable-length input must shard correctly: {}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(decompress_text(&out, d), expected, "sharded variable-length archive must roundtrip");
}

#[test]
fn fixed_worker_failure_errors_no_output() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fastq");
    write_fixed_fastq(&input, 400, 60);
    let out = d.join("out.qz");
    // QZ_NUMA_FAIL="1" fails worker-id 1 (0-based) — the middle of the 3 forced workers.
    let o = compress(&input, &out, d, "3",
        &[("QZ_NUMA_FORCE_WORKERS", "3"), ("QZ_NUMA_FAIL", "1")]);
    assert!(!o.status.success(), "fixed + worker failure must error");
    assert!(!out.exists(), "no output on failure");
}

#[test]
fn auto_worker_failure_falls_back_to_correct_output() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fastq");
    write_fixed_fastq(&input, 400, 60);
    let expected = std::fs::read_to_string(&input).unwrap();
    let out = d.join("out.qz");
    // QZ_NUMA_FAIL="1" fails worker-id 1 (0-based) — the middle of the 3 forced workers.
    let o = compress(&input, &out, d, "auto",
        &[("QZ_NUMA_FORCE_WORKERS", "3"), ("QZ_NUMA_FAIL", "1")]);
    assert!(o.status.success(), "auto must fall back: {}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(decompress_text(&out, d), expected);
}

#[test]
fn sharded_output_preserves_input_read_order() {
    // part_paths are indexed by splitter-range order (range 0 = earliest records), so
    // the merged archive is in input order regardless of worker completion order.
    // Distinct sequential headers make a reorder observable.
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fastq");
    let mut s = String::new();
    for i in 0..600 { s.push_str(&format!("@read_{i:04}\n{}\n+\n{}\n", "ACGT".repeat(15), "I".repeat(60))); }
    std::fs::write(&input, &s).unwrap();
    let expected = s.clone();
    let out = d.join("out.qz");
    let o = compress(&input, &out, d, "auto", &[("QZ_NUMA_FORCE_WORKERS", "3")]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(decompress_text(&out, d), expected, "merged output must be in input read order");
}

// Mode-routing prechecks: each unsupported mode → auto in-process (correct), fixed error.
#[test]
fn precheck_fasta_auto_in_process_fixed_errors() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fasta");
    let mut s = String::new();
    for i in 0..50 { s.push_str(&format!(">s{i}\nACGTACGTACGT\n")); }
    std::fs::write(&input, &s).unwrap();
    let out = d.join("out.qz");
    let o = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                  "-w", d.to_str().unwrap(), "--fasta", "--numa", "auto", "-f"],
               &[("QZ_NUMA_FORCE_WORKERS", "3")]);
    assert!(o.status.success(), "fasta+auto must run in-process: {}", String::from_utf8_lossy(&o.stderr));
    let o2 = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                   "-w", d.to_str().unwrap(), "--fasta", "--numa", "2", "-f"],
                &[("QZ_NUMA_FORCE_WORKERS", "3")]);
    assert!(!o2.status.success(), "fasta+fixed must error");
}

#[test]
fn precheck_ultra_auto_in_process_fixed_errors() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fastq");
    write_fixed_fastq(&input, 200, 60);
    let out = d.join("out.qz");
    let o = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                  "-w", d.to_str().unwrap(), "--ultra", "1", "--numa", "auto", "-f"],
               &[("QZ_NUMA_FORCE_WORKERS", "2")]);
    assert!(o.status.success(), "ultra+auto must run in-process: {}", String::from_utf8_lossy(&o.stderr));
    let o2 = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                   "-w", d.to_str().unwrap(), "--ultra", "1", "--numa", "2", "-f"],
                &[("QZ_NUMA_FORCE_WORKERS", "2")]);
    assert!(!o2.status.success(), "ultra+fixed must error");
}

#[test]
fn precheck_illumina_bin_quality_auto_in_process_fixed_errors() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fastq");
    write_fixed_fastq(&input, 200, 60);
    let out = d.join("out.qz");
    let o = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                  "-w", d.to_str().unwrap(), "--quality-mode", "illumina-bin", "--numa", "auto", "-f"],
               &[("QZ_NUMA_FORCE_WORKERS", "3")]);
    assert!(o.status.success(), "illumina-bin+auto must run in-process: {}", String::from_utf8_lossy(&o.stderr));
    let o2 = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                   "-w", d.to_str().unwrap(), "--quality-mode", "illumina-bin", "--numa", "2", "-f"],
                &[("QZ_NUMA_FORCE_WORKERS", "3")]);
    assert!(!o2.status.success(), "illumina-bin+fixed must error");
}

#[test]
fn precheck_gzip_input_auto_in_process() {
    use std::io::Write;
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fastq.gz");
    let mut s = String::new();
    for i in 0..50 { s.push_str(&format!("@r{i}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")); }
    let f = std::fs::File::create(&input).unwrap();
    let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
    enc.write_all(s.as_bytes()).unwrap();
    enc.finish().unwrap();
    let out = d.join("out.qz");
    let o = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                  "-w", d.to_str().unwrap(), "--numa", "auto", "-f"],
               &[("QZ_NUMA_FORCE_WORKERS", "3")]);
    assert!(o.status.success(), "gzip input + auto must run in-process: {}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(decompress_text(&out, d), s, "gzip-input in-process roundtrip must be lossless");
}

#[test]
fn precheck_gzip_input_fixed_errors() {
    use std::io::Write;
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fastq.gz");
    let mut s = String::new();
    for i in 0..50 { s.push_str(&format!("@r{i}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")); }
    let f = std::fs::File::create(&input).unwrap();
    let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
    enc.write_all(s.as_bytes()).unwrap();
    enc.finish().unwrap();
    let out = d.join("out.qz");
    let o = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                  "-w", d.to_str().unwrap(), "--numa", "2", "-f"],
               &[("QZ_NUMA_FORCE_WORKERS", "2")]);
    assert!(!o.status.success(), "gzip input + fixed must error (not silently in-process)");
}

#[test]
fn precheck_quality_mode_discard_auto_in_process_fixed_errors() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fastq");
    write_fixed_fastq(&input, 200, 60);
    let out = d.join("out.qz");
    let o = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                  "-w", d.to_str().unwrap(), "--quality-mode", "discard", "--numa", "auto", "-f"],
               &[("QZ_NUMA_FORCE_WORKERS", "3")]);
    assert!(o.status.success(), "discard+auto must run in-process: {}", String::from_utf8_lossy(&o.stderr));
    let o2 = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                   "-w", d.to_str().unwrap(), "--quality-mode", "discard", "--numa", "2", "-f"],
                &[("QZ_NUMA_FORCE_WORKERS", "3")]);
    assert!(!o2.status.success(), "discard+fixed must error");
}

/// Reference compression IS NUMA-shardable now (it was an "unsupported mode" before
/// the reference-aware merge landed): `--reference` + a fixed `--numa 2` must SHARD
/// and produce a lossless paired-reference (`archive_type 2`) archive, not error.
/// (Reads here are all-fallback against a tiny reference — exercises the
/// empty-union-imap merge path.)
#[test]
fn reference_fixed_shards_losslessly() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let r1 = d.join("r1.fastq");
    let r2 = d.join("r2.fastq");
    write_fixed_fastq(&r1, 50, 60);
    write_fixed_fastq(&r2, 50, 60);
    let exp1 = std::fs::read_to_string(&r1).unwrap();
    let exp2 = std::fs::read_to_string(&r2).unwrap();
    let reff = d.join("ref.fa");
    std::fs::write(&reff, ">c\nACGTACGTACGTACGT\n").unwrap();
    let out = d.join("out.qz");
    let o = run(&["compress", "-i", r1.to_str().unwrap(), "-i", r2.to_str().unwrap(),
                  "-o", out.to_str().unwrap(), "-w", d.to_str().unwrap(),
                  "--reference", reff.to_str().unwrap(), "--build-index", "--numa", "2", "-f"],
               &[("QZ_NUMA_FORCE_WORKERS", "2")]);
    assert!(o.status.success(), "reference + fixed must shard (not error): {}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(std::fs::read(&out).unwrap()[3], 2, "archive_type byte must be 2 (paired reference)");

    let dec = d.join("dec");
    let o2 = run(&["decompress", "-i", out.to_str().unwrap(), "-o", dec.to_str().unwrap(), "-w", d.to_str().unwrap(), "-f"], &[]);
    assert!(o2.status.success(), "decompress failed: {}", String::from_utf8_lossy(&o2.stderr));
    assert_eq!(std::fs::read_to_string(d.join("dec_R1.fastq")).unwrap(), exp1, "sharded reference R1 lossless");
    assert_eq!(std::fs::read_to_string(d.join("dec_R2.fastq")).unwrap(), exp2, "sharded reference R2 lossless");
}

#[test]
fn invalid_input_count_routes_in_process_and_errors() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let a = d.join("a.fastq"); let b = d.join("b.fastq"); let c = d.join("c.fastq");
    for p in [&a, &b, &c] { write_fixed_fastq(p, 50, 60); }
    let out = d.join("out.qz");
    for numa in ["auto", "2"] {
        let o = run(&["compress", "-i", a.to_str().unwrap(), "-i", b.to_str().unwrap(),
                      "-i", c.to_str().unwrap(), "-o", out.to_str().unwrap(),
                      "-w", d.to_str().unwrap(), "--numa", numa, "-f"],
                   &[("QZ_NUMA_FORCE_WORKERS", "3")]);
        assert!(!o.status.success(), "3 inputs under --numa {numa} must error");
        assert!(!out.exists(), "no archive on a 3-input config under --numa {numa}");
    }
}

#[test]
fn plain_paired_shards_under_auto_and_fixed() {
    // Paired NUMA sharding IS built: 2-input plain paired shards under BOTH auto and
    // fixed, producing a valid, verifiable, lossless paired archive (was: rejected).
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let r1 = d.join("r1.fastq");
    let r2 = d.join("r2.fastq");
    write_fixed_paired(&r1, &r2, 200, 60);
    let exp1 = std::fs::read_to_string(&r1).unwrap();
    let exp2 = std::fs::read_to_string(&r2).unwrap();
    for numa in ["auto", "2"] {
        let out = d.join("out.qz");
        let o = run(&["compress", "-i", r1.to_str().unwrap(), "-i", r2.to_str().unwrap(),
                      "-o", out.to_str().unwrap(), "-w", d.to_str().unwrap(),
                      "-t", "4", "--numa", numa, "-f"],
                   &[("QZ_NUMA_FORCE_WORKERS", "2")]);
        assert!(o.status.success(), "plain paired + {numa} must shard: {}", String::from_utf8_lossy(&o.stderr));
        let v = run(&["verify", "-i", out.to_str().unwrap()], &[]);
        assert!(v.status.success(), "paired archive ({numa}) must verify: {}", String::from_utf8_lossy(&v.stderr));
        let dec = d.join("dec");
        let od = run(&["decompress", "-i", out.to_str().unwrap(), "-o", dec.to_str().unwrap(),
                       "-w", d.to_str().unwrap(), "-f"], &[]);
        assert!(od.status.success(), "{}", String::from_utf8_lossy(&od.stderr));
        assert_eq!(std::fs::read_to_string(d.join("dec_R1.fastq")).unwrap(), exp1, "{numa} R1 roundtrip");
        assert_eq!(std::fs::read_to_string(d.join("dec_R2.fastq")).unwrap(), exp2, "{numa} R2 roundtrip");
    }
}

#[test]
fn preexisting_output_without_force_refused() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let input = d.join("in.fastq");
    write_fixed_fastq(&input, 400, 60);
    let out = d.join("out.qz");
    std::fs::write(&out, b"existing").unwrap();
    let o = run(&["compress", "-i", input.to_str().unwrap(), "-o", out.to_str().unwrap(),
                  "-w", d.to_str().unwrap(), "-t", "4", "--numa", "auto"],
               &[("QZ_NUMA_FORCE_WORKERS", "3")]);
    assert!(!o.status.success(), "must refuse pre-existing output without --force");
}

// ---------------------------------------------------------------------------
// Reference compress sharding (--reference): byte-range workers each produce a
// reference part archive (type 4 single / type 2 paired), stitched by the
// reference-aware merge (re-derives the coverage globals over the union of covered
// intervals). Oracle = lossless roundtrip of the merged archive vs the original reads.
// ---------------------------------------------------------------------------

/// Deterministic low-repeat ACGT sequence so reads seed reliably against the mapper.
fn ref_seq(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        v.push(b"ACGT"[((x >> 33) & 3) as usize]);
    }
    v
}

/// Write a single-end reference dataset: ~3 kb FASTA + `n` forward substrings
/// (100–120 bp), a few with planted substitutions, the last two random (fallback).
/// Returns `(fasta, reads, reads_text)`.
fn write_ref_single(dir: &Path, n: usize) -> (PathBuf, PathBuf, String) {
    let r = ref_seq(3000, 7);
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">c\n{}\n", std::str::from_utf8(&r).unwrap())).unwrap();
    let mut s = String::new();
    let mut rng = 0x1234_5678_9abc_def0u64;
    let mut next = || { rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (rng >> 33) as usize };
    for i in 0..n {
        let rlen = 100 + (next() % 21);
        let bytes = if i >= n - 2 {
            ref_seq(rlen, 9000 + i as u64)
        } else {
            let st = next() % (r.len() - rlen + 1);
            let mut b = r[st..st + rlen].to_vec();
            if i % 17 == 0 { let p = next() % rlen; b[p] = if b[p] == b'A' { b'C' } else { b'A' }; }
            b
        };
        let q = "I".repeat(bytes.len());
        s.push_str(&format!("@read_{i}\n{}\n+\n{q}\n", std::str::from_utf8(&bytes).unwrap()));
    }
    let reads = dir.join("reads.fastq");
    std::fs::write(&reads, &s).unwrap();
    (fa, reads, s)
}

#[test]
fn forced_shard_single_end_reference_roundtrips_lossless() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let (fa, reads, expected) = write_ref_single(d, 400);

    let out = d.join("ref_shard.qz");
    let o = run(
        &["compress", "-i", reads.to_str().unwrap(), "--reference", fa.to_str().unwrap(),
          "--build-index", "-o", out.to_str().unwrap(), "-w", d.to_str().unwrap(), "-t", "4", "--numa", "auto", "-f"],
        &[("QZ_NUMA_FORCE_WORKERS", "3")],
    );
    assert!(o.status.success(), "forced-shard single-end reference compress failed: {}", String::from_utf8_lossy(&o.stderr));

    // archive_type byte (header byte 3) == 4 (single-end reference).
    let bytes = std::fs::read(&out).unwrap();
    assert_eq!(bytes[3], 4, "archive_type byte must be 4 (single-end reference)");

    // Lossless: the merged sharded archive decodes to the original reads.
    assert_eq!(decompress_text(&out, d), expected, "sharded reference archive must decode to the input reads");

    // No leftover .part / temp files in the working dir.
    let leftover: Vec<_> = std::fs::read_dir(d).unwrap()
        .filter_map(|e| e.ok()).map(|e| e.file_name().into_string().unwrap())
        .filter(|n| n.contains(".part") || n.starts_with(".qz_numa")).collect();
    assert!(leftover.is_empty(), "no part/temp files should remain: {leftover:?}");
}

#[test]
fn forced_shard_paired_reference_roundtrips_lossless() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let r = ref_seq(4000, 11);
    let fa = d.join("refp.fa");
    std::fs::write(&fa, format!(">c\n{}\n", std::str::from_utf8(&r).unwrap())).unwrap();
    let comp = |b: u8| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', _ => b'A' };
    let (mut s1, mut s2) = (String::new(), String::new());
    let mut rng = 0xfeed_face_dead_beefu64;
    let mut next = || { rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (rng >> 33) as usize };
    let n = 400;
    for i in 0..n {
        let rlen = 100 + (next() % 21);
        let (a, b) = if i >= n - 2 {
            (ref_seq(rlen, 7000 + i as u64), ref_seq(rlen, 8000 + i as u64))
        } else {
            let insert = 250usize;
            let st = next() % (r.len() - (insert + rlen) + 1);
            let a = r[st..st + rlen].to_vec();
            let b: Vec<u8> = r[st + insert..st + insert + rlen].iter().rev().map(|&x| comp(x)).collect();
            (a, b)
        };
        let q = "I".repeat(rlen);
        s1.push_str(&format!("@p_{i}/1\n{}\n+\n{q}\n", std::str::from_utf8(&a).unwrap()));
        s2.push_str(&format!("@p_{i}/2\n{}\n+\n{q}\n", std::str::from_utf8(&b).unwrap()));
    }
    let (r1, r2) = (d.join("r1.fastq"), d.join("r2.fastq"));
    std::fs::write(&r1, &s1).unwrap();
    std::fs::write(&r2, &s2).unwrap();

    let out = d.join("ref_paired_shard.qz");
    let o = run(
        &["compress", "-i", r1.to_str().unwrap(), "-i", r2.to_str().unwrap(),
          "--reference", fa.to_str().unwrap(), "--build-index", "-o", out.to_str().unwrap(),
          "-w", d.to_str().unwrap(), "-t", "4", "--numa", "auto", "-f"],
        &[("QZ_NUMA_FORCE_WORKERS", "3")],
    );
    assert!(o.status.success(), "forced-shard paired reference compress failed: {}", String::from_utf8_lossy(&o.stderr));

    let bytes = std::fs::read(&out).unwrap();
    assert_eq!(bytes[3], 2, "archive_type byte must be 2 (paired reference)");

    let dec = d.join("dec");
    let o2 = run(&["decompress", "-i", out.to_str().unwrap(), "-o", dec.to_str().unwrap(), "-w", d.to_str().unwrap(), "-f"], &[]);
    assert!(o2.status.success(), "paired reference decompress failed: {}", String::from_utf8_lossy(&o2.stderr));
    assert_eq!(std::fs::read_to_string(d.join("dec_R1.fastq")).unwrap(), s1, "sharded reference R1 must roundtrip");
    assert_eq!(std::fs::read_to_string(d.join("dec_R2.fastq")).unwrap(), s2, "sharded reference R2 must roundtrip");
}
