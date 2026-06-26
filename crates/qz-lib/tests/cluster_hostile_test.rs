//! Hostile / structural decode tests for cluster mode (archive_type 3).
//!
//! Mirrors the reference suites' `expect_reject` harness: mutate a REAL on-disk
//! cluster archive (footer flips, dropped globals, inflated counts, oversized
//! spans, CRC/payload tamper, truncation, record_count mismatch) and assert the
//! PUBLIC `decompress()` and `verify()` (fast + deep) return `Err` — never panic —
//! and leak no output file. This exercises the real
//! parse → validate → allocate → atomic-cleanup pipeline that the hand-built
//! `ChunkDirectory` unit tests in `cluster/decode.rs` skip. Cluster drops order,
//! so there is no byte-identical golden net here — structural rejection IS the net.

use qz_lib::cli::{ClusterOptions, CompressConfig};
use std::path::{Path, PathBuf};

// ── v5 footer wire-format view (mode-agnostic; mirrors chunk_directory.rs) ──────
const ENTRY_BYTES: usize = 31; // chunk_index(4)+role(1)+mate(1)+codec(1)+off(8)+len(8)+rc(8)
const V5_LOCATOR_LEN: usize = 20;
const V5_FOOTER_PREFIX: usize = 8 + 4 + 4; // num_reads + num_chunks + n_entries

// Cluster StreamRole discriminants (chunk_directory.rs).
const R_HEADERS: u8 = 0;
const R_CLUSTERED_SEQ: u8 = 19;
const R_CLUSTER_META: u8 = 21;

#[derive(Clone, Copy)]
struct FEntry {
    abs: usize,
    role: u8,
    offset: u64,
}
struct FView {
    entries: Vec<FEntry>,
}
fn read_u32(b: &[u8], o: usize) -> u32 {
    u32::from_le_bytes(b[o..o + 4].try_into().unwrap())
}
fn read_u64(b: &[u8], o: usize) -> u64 {
    u64::from_le_bytes(b[o..o + 8].try_into().unwrap())
}
fn v5_footer_crc(footer: &[u8]) -> u32 {
    let mut c = flate2::Crc::new();
    c.update(footer);
    c.sum()
}
fn v5_footer_span(bytes: &[u8]) -> (usize, usize) {
    let loc = bytes.len() - V5_LOCATOR_LEN;
    assert_eq!(&bytes[loc + 12..], b"QZFOOTR1", "not a v5 archive");
    let footer_len = u64::from_le_bytes(bytes[loc..loc + 8].try_into().unwrap()) as usize;
    (loc - footer_len, footer_len)
}
fn v5_refresh_footer_crc(bytes: &mut [u8]) {
    let (fs, fl) = v5_footer_span(bytes);
    let crc = v5_footer_crc(&bytes[fs..fs + fl]);
    let loc = bytes.len() - V5_LOCATOR_LEN;
    bytes[loc + 8..loc + 12].copy_from_slice(&crc.to_le_bytes());
}
fn fview(bytes: &[u8]) -> FView {
    let (footer_offset, _) = v5_footer_span(bytes);
    let n = read_u32(bytes, footer_offset + 12) as usize;
    let mut entries = Vec::with_capacity(n);
    let mut o = footer_offset + V5_FOOTER_PREFIX;
    for _ in 0..n {
        let role = bytes[o + 4];
        let offset = read_u64(bytes, o + 7);
        entries.push(FEntry { abs: o, role, offset });
        o += ENTRY_BYTES;
    }
    FView { entries }
}
fn find_role(v: &FView, role: u8) -> FEntry {
    *v.entries.iter().find(|e| e.role == role).expect("role present")
}

/// Build one valid single-end cluster archive on disk; return (bytes, archive_path).
fn valid_cluster_archive(d: &Path) -> (Vec<u8>, PathBuf) {
    let alphabet = b"ACGT";
    let mut buf = String::new();
    for i in 0..300usize {
        let seq: String = (0..100)
            .map(|j| alphabet[(i.wrapping_mul(37).wrapping_add(j).wrapping_mul(13)) % 4] as char)
            .collect();
        let qual: String = (0..100u8)
            .map(|j| (b'!' + ((i as u8).wrapping_add(j) % 40)) as char)
            .collect();
        buf.push_str(&format!("@read_{i} lane:1:1:1:{i}\n{seq}\n+\n{qual}\n"));
    }
    let input = d.join("in.fastq");
    std::fs::write(&input, buf.as_bytes()).unwrap();

    let wdir = d.join("wdir");
    std::fs::create_dir_all(&wdir).unwrap();
    let arc = d.join("valid.qz");
    let cfg = CompressConfig {
        input: vec![input],
        output: arc.clone(),
        working_dir: wdir,
        threads: 2,
        cluster: Some(ClusterOptions::default()),
        force: true,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();
    let bytes = std::fs::read(&arc).unwrap();
    assert_eq!(bytes[3], 3, "must be a cluster archive");
    (bytes, arc)
}

/// Apply `f` to a COPY of `bytes`, write it, and assert PUBLIC decompress + verify
/// (fast & deep) all return Err with no leaked output and no panic.
fn expect_reject(d: &Path, case: &str, bytes: &[u8], f: impl FnOnce(&mut Vec<u8>)) {
    let mut b = bytes.to_vec();
    f(&mut b);
    let archive = d.join(format!("corrupt_{case}.qz"));
    std::fs::write(&archive, &b).unwrap();

    let outp = d.join(format!("dec_{case}.fastq"));
    let _ = std::fs::remove_file(&outp);

    let dc = qz_lib::cli::DecompressConfig {
        input: archive.clone(),
        output: vec![outp.clone()],
        working_dir: d.to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        qz_lib::compression::decompress(&dc)
    }));
    match res {
        Err(_) => panic!("[{case}] decompress PANICKED (missing guard)"),
        Ok(r) => assert!(r.is_err(), "[{case}] decompress must return Err"),
    }
    assert!(!outp.exists(), "[{case}] output leaked on failed decode");

    for fast in [true, false] {
        let vc = qz_lib::cli::VerifyConfig {
            input: archive.clone(),
            working_dir: d.to_path_buf(),
            num_threads: 1,
            fast,
        };
        let vr = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            qz_lib::compression::verify(&vc)
        }));
        match vr {
            Err(_) => panic!("[{case}] verify(fast={fast}) PANICKED (missing guard)"),
            Ok(r) => assert!(r.is_err(), "[{case}] verify(fast={fast}) must return Err"),
        }
    }
}

/// A FASTQ record whose sequence and quality lengths differ is rejected by cluster
/// compress (mirroring the standard reader, `io/fastq.rs`). Cluster decode reconstructs
/// each record's sequence boundary from its quality length, so accepting a mismatch would
/// silently corrupt the output or render the archive undecodable.
#[test]
fn cluster_rejects_seq_qual_length_mismatch() {
    let d = tempfile::tempdir().unwrap();
    let mut buf = String::new();
    // Enough well-formed records that cluster mode runs its normal path...
    for i in 0..200usize {
        let seq = "ACGT".repeat(25); // 100 bp
        let qual = "I".repeat(100);
        buf.push_str(&format!("@read_{i}\n{seq}\n+\n{qual}\n"));
    }
    // ...then one malformed record: 100-base sequence, 99-char quality.
    buf.push_str(&format!("@bad_read\n{}\n+\n{}\n", "A".repeat(100), "I".repeat(99)));

    let input = d.path().join("in.fastq");
    std::fs::write(&input, buf.as_bytes()).unwrap();
    let wdir = d.path().join("wdir");
    std::fs::create_dir_all(&wdir).unwrap();
    let arc = d.path().join("out.qz");

    let cfg = CompressConfig {
        input: vec![input],
        output: arc.clone(),
        working_dir: wdir,
        threads: 2,
        cluster: Some(ClusterOptions::default()),
        force: true,
        ..CompressConfig::default()
    };
    let res = qz_lib::compression::compress(&cfg);
    assert!(res.is_err(), "cluster compress must reject a seq/qual length mismatch");
    assert!(!arc.exists(), "no archive must be left behind on rejected malformed input");
}

#[test]
fn cluster_hostile_unknown_archive_type() {
    let d = tempfile::tempdir().unwrap();
    let (bytes, _) = valid_cluster_archive(d.path());
    expect_reject(d.path(), "unknown_type", &bytes, |b| b[3] = 9);
}

#[test]
fn cluster_hostile_footer_and_truncation() {
    let d = tempfile::tempdir().unwrap();
    let (bytes, _) = valid_cluster_archive(d.path());

    // Footer-body byte flipped WITHOUT recomputing the locator CRC → CRC mismatch.
    expect_reject(d.path(), "footer_crc_tamper", &bytes, |b| {
        let (fs, _) = v5_footer_span(b);
        b[fs] ^= 0xFF;
    });

    // Truncate the tail into the locator → magic gone.
    expect_reject(d.path(), "trunc_tail", &bytes, |b| {
        let n = b.len();
        b.truncate(n - 32);
    });

    // Truncate at the footer start (one byte past the last payload).
    expect_reject(d.path(), "trunc_at_footer_start", &bytes, |b| {
        let (fs, _) = v5_footer_span(b);
        b.truncate(fs);
    });
}

#[test]
fn cluster_hostile_directory_structure() {
    let d = tempfile::tempdir().unwrap();
    let (bytes, _) = valid_cluster_archive(d.path());

    // Inflated entry count (→ footer alloc-cap / per-entry span rejection).
    expect_reject(d.path(), "inflated_n_entries", &bytes, |b| {
        let (fs, _) = v5_footer_span(b);
        b[fs + 12..fs + 16].copy_from_slice(&u32::MAX.to_le_bytes());
        v5_refresh_footer_crc(b);
    });

    // Drop the ClusterMeta global: overwrite its role → "missing ClusterMeta".
    expect_reject(d.path(), "dropped_cluster_meta", &bytes, |b| {
        let v = fview(b);
        let e = find_role(&v, R_CLUSTER_META);
        b[e.abs + 4] = R_HEADERS;
        v5_refresh_footer_crc(b);
    });

    // Oversized ClusteredSeq length (→ per-entry span/overflow rejection).
    expect_reject(d.path(), "oversized_clustered_seq", &bytes, |b| {
        let v = fview(b);
        let e = find_role(&v, R_CLUSTERED_SEQ);
        b[e.abs + 15..e.abs + 23].copy_from_slice(&(u64::MAX / 2).to_le_bytes());
        v5_refresh_footer_crc(b);
    });
}

#[test]
fn cluster_hostile_payload_and_record_count() {
    let d = tempfile::tempdir().unwrap();
    let (bytes, _) = valid_cluster_archive(d.path());

    // Flip the first payload byte of a ClusteredSeq frame → per-block CRC mismatch
    // (footer CRC unchanged — the corruption is in the block, not the footer).
    expect_reject(d.path(), "payload_flip", &bytes, |b| {
        let v = fview(b);
        let e = find_role(&v, R_CLUSTERED_SEQ);
        let payload0 = e.offset as usize + 12; // after the 12-byte frame header
        b[payload0] ^= 0xFF;
    });

    // Forge the footer entry record_count so it disagrees with the inline frame
    // record_count → the read_frame / fast-verify cross-check fires.
    expect_reject(d.path(), "record_count_mismatch", &bytes, |b| {
        let v = fview(b);
        let e = find_role(&v, R_HEADERS);
        let rc = read_u64(b, e.abs + 23);
        b[e.abs + 23..e.abs + 31].copy_from_slice(&(rc + 1).to_le_bytes());
        v5_refresh_footer_crc(b);
    });
}
