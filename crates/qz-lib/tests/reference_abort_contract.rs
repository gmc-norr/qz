//! Runtime regression for the reference (`--reference`) streaming compress
//! abort/shutdown contract (worker error + producer error each abort cleanly).
//!
//! Both tests drive a debug-only failpoint (`QZ_REF_FAIL`, read once in
//! `compression/reference/stream.rs`) that bails from a worker (`run_ref_job`) or
//! the producer (`run_reference_producer`) with a recognizable message. The hook is
//! `#[cfg(debug_assertions)]` so it compiles OUT of release builds entirely — and
//! THESE TESTS ARE `#[cfg(debug_assertions)]` TOO: under `cargo test --release` the
//! hook is inert, the compress would SUCCEED, and the `expect_err` below would
//! panic. (This exact gating bug was already fixed once this session for
//! `decode_chunk_range_bounded_overrun_is_integrity_error`; we do not repeat it.)
//!
//! ENV ISOLATION: `QZ_REF_FAIL` is process-global. Cargo runs each integration test
//! file as its OWN binary/process, so putting these two tests in their own file
//! guarantees the var can never leak into `reference_integration.rs`'s 27 other
//! reference compresses (a different process). Within THIS file, a poison-tolerant
//! lock serializes the set→compress→remove window so the two tests can't overwrite
//! each other's value. `set_var`/`remove_var` are `unsafe` in edition 2024
//! (codebase convention, same as the NUMA tests).
//!
//! The contract asserted (per Increment-3 spec): the injected error is SURFACED
//! (worker error authoritative; producer error surfaced via `w.and(prod)` after the
//! writer's clean abort return), the call RETURNS (no hang/deadlock — the test
//! simply completing is the no-hang assertion), and NO partial archive is left at
//! the output path (outer `compress()`'s O_EXCL temp + drop-guard cleanup).

#![cfg(debug_assertions)]

use qz_lib::cli::{CompressConfig, ReferenceOptions};
use std::path::PathBuf;
use std::sync::Mutex;

/// Serializes the env-var set/compress/remove window so the two tests in this file
/// (run in parallel by the test harness) can't clobber each other's `QZ_REF_FAIL`.
/// Poison-tolerant: a panicking test still leaves the lock usable.
static REF_FAIL_LOCK: Mutex<()> = Mutex::new(());

fn w(p: &std::path::Path, s: &str) {
    std::fs::write(p, s).unwrap();
}

fn cfg_ref(r1: PathBuf, r2: PathBuf, out: PathBuf, tmp: PathBuf, refp: PathBuf) -> CompressConfig {
    CompressConfig {
        input: vec![r1, r2],
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

/// Deterministic low-repeat ACGT sequence (matches mapping.rs's generator so reads
/// seed reliably). Same as `reference_integration.rs`'s `make_seq`.
fn make_seq(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        v.push(b"ACGT"[((x >> 33) & 3) as usize]);
    }
    v
}

/// Build a small multi-chunk reference dataset: a 2000bp reference + `n_pairs`
/// exact-substring read pairs (R1 forward, R2 a downstream reverse-complement).
/// With `chunk_records < n_pairs` the producer dispatches several chunks, so a
/// non-last worker / a mid-stream producer failpoint is reachable.
fn make_dataset(dir: &std::path::Path, n_pairs: usize) -> (PathBuf, PathBuf, PathBuf) {
    let refseq = make_seq(2000, 7);
    let rf = dir.join("ref.fa");
    {
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        w(&rf, &s);
    }

    fn revcomp(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => other,
            })
            .collect()
    }

    let mut s1 = String::new();
    let mut s2 = String::new();
    let mut rng = 0x1234_5678_9abc_def0u64;
    let mut next = || {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (rng >> 33) as usize
    };
    let rlen = 120usize;
    let max_start = refseq.len() - rlen;
    for i in 0..n_pairs {
        let st1 = next() % (max_start + 1);
        let r1 = refseq[st1..st1 + rlen].to_vec();
        let st2 = (st1 + 150).min(max_start);
        let r2 = revcomp(&refseq[st2..st2 + rlen]);
        let q: String = "I".repeat(rlen);
        s1.push_str(&format!(
            "@read_{i}/1\n{}\n+\n{q}\n",
            std::str::from_utf8(&r1).unwrap()
        ));
        s2.push_str(&format!(
            "@read_{i}/2\n{}\n+\n{q}\n",
            std::str::from_utf8(&r2).unwrap()
        ));
    }
    let r1p = dir.join("r1.fastq");
    let r2p = dir.join("r2.fastq");
    w(&r1p, &s1);
    w(&r2p, &s2);
    (rf, r1p, r2p)
}

/// Set `QZ_REF_FAIL=value`, run `compress(cfg)`, then unconditionally clear it — all
/// under `REF_FAIL_LOCK` so the env mutation never races the sibling test.
fn compress_with_ref_fail(cfg: &CompressConfig, value: &str) -> anyhow::Result<()> {
    let _g = REF_FAIL_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    unsafe { std::env::set_var("QZ_REF_FAIL", value) };
    let res = qz_lib::compression::compress(cfg);
    unsafe { std::env::remove_var("QZ_REF_FAIL") };
    res
}

/// A WORKER error (`run_ref_job` bails on a target chunk) is AUTHORITATIVE: its
/// message surfaces (even though the producer may secondarily see "workers exited"),
/// the call returns (no hang), and no partial archive is left.
#[test]
fn reference_worker_error_aborts_clean_no_partial_archive() {
    let d = tempfile::tempdir().unwrap();
    let (refp, r1p, r2p) = make_dataset(d.path(), 200);
    let out = d.path().join("worker_fail.qz");
    let mut c = cfg_ref(r1p, r2p, out.clone(), d.path().to_path_buf(), refp);
    c.advanced.chunk_records = 50; // 200 pairs / 50 = 4 chunks → fail a non-last one

    // Fail chunk index 1 (mid-stream worker, not the last chunk).
    let err = compress_with_ref_fail(&c, "worker:1")
        .expect_err("injected worker failure must abort the compress");
    assert!(
        err.to_string().contains("injected worker failure at chunk 1"),
        "the WORKER error must be authoritative, got: {err}"
    );
    assert!(
        !out.exists(),
        "no partial archive may be left at the output path after a worker abort"
    );
}

/// A PRODUCER error (`run_reference_producer` bails after dispatching N chunks) sets
/// `abort` + drops its senders → the writer drains, sees `abort`, returns Ok WITHOUT
/// finalizing → `w.and(prod)` surfaces the PRODUCER error. The call returns (no
/// hang), and no partial archive is left.
#[test]
fn reference_producer_error_aborts_clean_no_partial_archive() {
    let d = tempfile::tempdir().unwrap();
    let (refp, r1p, r2p) = make_dataset(d.path(), 200);
    let out = d.path().join("producer_fail.qz");
    let mut c = cfg_ref(r1p, r2p, out.clone(), d.path().to_path_buf(), refp);
    c.advanced.chunk_records = 50; // 200 pairs / 50 = 4 chunks

    // Fail after dispatching 2 chunks (mid-stream, before the final globals send).
    let err = compress_with_ref_fail(&c, "producer:2")
        .expect_err("injected producer failure must abort the compress");
    assert!(
        err.to_string().contains("injected producer failure after 2 chunks"),
        "the PRODUCER error must be surfaced via w.and(prod), got: {err}"
    );
    assert!(
        !out.exists(),
        "no partial archive may be left at the output path after a producer abort"
    );
}
