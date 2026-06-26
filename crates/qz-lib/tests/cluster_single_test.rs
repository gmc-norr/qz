/// End-to-end integration tests for the cluster single-end compression mode.
///
/// These tests drive the mode entirely through the PUBLIC library API:
///   qz_lib::compression::compress / decompress / verify
///   qz_lib::cli::{CompressConfig, DecompressConfig, VerifyConfig, ClusterOptions}
///
/// No pub(crate) cluster internals are called here — this compiles as a separate
/// external test crate and can only see the public surface.
use qz_lib::cli::{ClusterOptions, CompressConfig, DecompressConfig, VerifyConfig};
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

// ─── inline helpers ───────────────────────────────────────────────────────────

/// Write a varied ~2000-read FASTQ to `path`.
///
/// Includes:
/// (a) reads with N bases scattered in,
/// (b) a run of all-identical reads (exercises deterministic tie-break),
/// (c) varied ACGT reads of the standard length,
/// (d) reads of a different length (variable-length support).
/// All quality strings match the length of their sequence.
fn write_fastq_mixed(path: &Path) {
    let mut buf = String::new();

    // (a) 200 reads with N bases interspersed
    let bases_n = b"ACGTNACGTNACGTACGT"; // 18 bp, some N
    for i in 0..200usize {
        let seq: String = bases_n
            .iter()
            .enumerate()
            .map(|(j, &b)| if (i + j) % 7 == 0 { b'N' } else { b } as char)
            .collect();
        let qual = "I".repeat(seq.len());
        buf.push_str(&format!("@read_n_{i} lane:1:2:3:{i}\n{seq}\n+\n{qual}\n"));
    }

    // (b) 200 all-identical reads (deterministic tie-break exercise)
    let identical_seq = "ACGTACGTACGTACGTACGT"; // 20 bp
    let identical_qual = "H".repeat(identical_seq.len());
    for i in 0..200usize {
        buf.push_str(&format!(
            "@read_id_{i} tile:1:1:1:{i}\n{identical_seq}\n+\n{identical_qual}\n"
        ));
    }

    // (c) 1500 varied ACGT reads (standard 100 bp)
    let alphabet = b"ACGT";
    for i in 0..1500usize {
        let seq: String = (0..100)
            .map(|j| alphabet[(i.wrapping_mul(37).wrapping_add(j).wrapping_mul(13)) % 4] as char)
            .collect();
        let qual: String = (0..100u8)
            .map(|j| (b'!' + ((i as u8).wrapping_add(j) % 40)) as char)
            .collect();
        buf.push_str(&format!("@read_varied_{i} rg:1\n{seq}\n+\n{qual}\n"));
    }

    // (d) 100 reads of a different length (50 bp) for variable-length coverage
    for i in 0..100usize {
        let seq: String = (0..50)
            .map(|j| alphabet[(i.wrapping_add(j * 3)) % 4] as char)
            .collect();
        let qual = "B".repeat(50);
        buf.push_str(&format!("@read_short_{i} short\n{seq}\n+\n{qual}\n"));
    }

    fs::write(path, buf.as_bytes()).unwrap();
}

/// Parse `path` as a FASTQ and return a SORTED multiset of (header, seq, qual) byte
/// tuples. Sorting makes comparison order-independent (set-lossless). Because cluster
/// has no byte-identical golden, this helper IS the safety net — so it parses
/// STRICTLY: the file must be a whole number of well-formed 4-line records (`@`
/// header, `+` separator, seq.len() == qual.len()), not a lenient line scan that
/// would silently absorb a line-structure corruption.
fn multiset(path: &Path) -> Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> {
    let data = fs::read(path).unwrap();
    // A well-formed FASTQ ends with a trailing '\n', so the final split element is
    // empty — drop exactly that one, then require a whole number of 4-line records.
    let mut lines: Vec<&[u8]> = data.split(|&b| b == b'\n').collect();
    if lines.last() == Some(&&b""[..]) {
        lines.pop();
    }
    assert!(
        lines.len().is_multiple_of(4),
        "FASTQ {} is not a whole number of 4-line records ({} lines)",
        path.display(),
        lines.len()
    );
    let mut records = Vec::with_capacity(lines.len() / 4);
    for rec in lines.chunks_exact(4) {
        assert_eq!(rec[0].first(), Some(&b'@'), "record header must start with '@'");
        assert_eq!(rec[2].first(), Some(&b'+'), "record separator must be '+'");
        assert_eq!(rec[1].len(), rec[3].len(), "seq/qual length mismatch");
        records.push((rec[0].to_vec(), rec[1].to_vec(), rec[3].to_vec()));
    }
    records.sort_unstable();
    records
}

/// Build a CompressConfig with cluster mode enabled, using a dedicated working_dir subdir.
/// Returns the output archive path after calling compress.
#[allow(dead_code)]
fn compress_cluster_to(workdir: &Path, input: &Path, out_name: &str) -> PathBuf {
    let sub = workdir.join(out_name.replace('.', "_") + "_wdir");
    fs::create_dir_all(&sub).unwrap();
    let output = workdir.join(out_name);
    let cfg = CompressConfig {
        input: vec![input.to_path_buf()],
        output: output.clone(),
        working_dir: sub,
        threads: 2,
        cluster: Some(ClusterOptions::default()),
        force: true,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();
    output
}

/// Decompress a cluster archive using DecompressConfig, return the output path.
fn decompress_cluster_to(workdir: &Path, archive: &Path, out_name: &str) -> PathBuf {
    let output = workdir.join(out_name);
    let dc = DecompressConfig {
        input: archive.to_path_buf(),
        output: vec![output.clone()],
        working_dir: workdir.to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&dc).unwrap();
    output
}

/// Run verify on an archive; returns true iff it succeeds.
fn verify_ok(archive: &Path, workdir: &Path) -> bool {
    let vc = VerifyConfig {
        input: archive.to_path_buf(),
        working_dir: workdir.to_path_buf(),
        num_threads: 1,
        fast: false,
    };
    qz_lib::compression::verify(&vc).is_ok()
}

// ─── tests ────────────────────────────────────────────────────────────────────

/// Test 1: roundtrip set-lossless; byte-identical determinism (two compressions
/// of the same input must produce byte-identical archives).
#[test]
fn roundtrip_multiset_and_determinism() {
    let temp = TempDir::new().unwrap();
    let tp = temp.path();

    // Write mixed input
    let input = tp.join("input.fastq");
    write_fastq_mixed(&input);

    // Compress TWICE with separate working dirs
    let wdir_a = tp.join("wdir_a");
    fs::create_dir_all(&wdir_a).unwrap();
    let a = tp.join("a.qz");
    {
        let cfg = CompressConfig {
            input: vec![input.clone()],
            output: a.clone(),
            working_dir: wdir_a.clone(),
            threads: 2,
            cluster: Some(ClusterOptions::default()),
            force: true,
            ..CompressConfig::default()
        };
        qz_lib::compression::compress(&cfg).unwrap();
    }

    let wdir_b = tp.join("wdir_b");
    fs::create_dir_all(&wdir_b).unwrap();
    let b = tp.join("b.qz");
    {
        let cfg = CompressConfig {
            input: vec![input.clone()],
            output: b.clone(),
            working_dir: wdir_b.clone(),
            threads: 2,
            cluster: Some(ClusterOptions::default()),
            force: true,
            ..CompressConfig::default()
        };
        qz_lib::compression::compress(&cfg).unwrap();
    }

    let bytes_a = fs::read(&a).unwrap();
    let bytes_b = fs::read(&b).unwrap();
    assert_eq!(
        bytes_a, bytes_b,
        "cluster archive must be byte-identical for identical input (determinism)"
    );

    // Byte 3 must be archive_type 3 (cluster)
    assert_eq!(
        bytes_a[3], 3,
        "archive_type byte must be 3 for cluster archive"
    );

    // Decompress and verify set-lossless
    let wdir_dec = tp.join("wdir_dec");
    fs::create_dir_all(&wdir_dec).unwrap();
    let back = decompress_cluster_to(&wdir_dec, &a, "back.fastq");

    let ms_in = multiset(&input);
    let ms_back = multiset(&back);
    assert_eq!(
        ms_in.len(),
        ms_back.len(),
        "read count must be preserved (set-lossless)"
    );
    assert_eq!(ms_in, ms_back, "cluster round-trip is set-lossless");

    // Verify integrity
    assert!(verify_ok(&a, &wdir_dec), "verify must pass on cluster archive");
}

/// Test 2: variable-length and edge cases.
#[test]
fn variable_length_and_edge_cases() {
    let temp = TempDir::new().unwrap();
    let tp = temp.path();

    // --- sub-case: empty input ---
    {
        let sub = tp.join("empty");
        fs::create_dir_all(&sub).unwrap();
        let fq = sub.join("in.fastq");
        fs::write(&fq, b"").unwrap();

        let wdir = sub.join("wdir");
        fs::create_dir_all(&wdir).unwrap();
        let arc = sub.join("out.qz");
        let cfg = CompressConfig {
            input: vec![fq.clone()],
            output: arc.clone(),
            working_dir: wdir.clone(),
            threads: 1,
            cluster: Some(ClusterOptions::default()),
            force: true,
            ..CompressConfig::default()
        };
        qz_lib::compression::compress(&cfg).unwrap();

        let out_fq = sub.join("back.fastq");
        let dc = DecompressConfig {
            input: arc.clone(),
            output: vec![out_fq.clone()],
            working_dir: wdir.clone(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 6,
            force: true,
        };
        qz_lib::compression::decompress(&dc).unwrap();

        assert_eq!(multiset(&fq), multiset(&out_fq), "empty: set-lossless");
        assert!(verify_ok(&arc, &wdir), "empty: verify ok");
    }

    // --- sub-case: single record ---
    {
        let sub = tp.join("single");
        fs::create_dir_all(&sub).unwrap();
        let fq = sub.join("in.fastq");
        fs::write(&fq, b"@r1 desc\nACGTACGT\n+\nIIIIIIII\n").unwrap();

        let wdir = sub.join("wdir");
        fs::create_dir_all(&wdir).unwrap();
        let arc = sub.join("out.qz");
        let cfg = CompressConfig {
            input: vec![fq.clone()],
            output: arc.clone(),
            working_dir: wdir.clone(),
            threads: 1,
            cluster: Some(ClusterOptions::default()),
            force: true,
            ..CompressConfig::default()
        };
        qz_lib::compression::compress(&cfg).unwrap();

        let out_fq = sub.join("back.fastq");
        let dc = DecompressConfig {
            input: arc.clone(),
            output: vec![out_fq.clone()],
            working_dir: wdir.clone(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 6,
            force: true,
        };
        qz_lib::compression::decompress(&dc).unwrap();

        assert_eq!(multiset(&fq), multiset(&out_fq), "single: set-lossless");
        assert!(verify_ok(&arc, &wdir), "single: verify ok");
    }

    // --- sub-case: variable read lengths (50/100/150 bp mixed) ---
    {
        let sub = tp.join("varlen");
        fs::create_dir_all(&sub).unwrap();
        let fq = sub.join("in.fastq");
        let mut data = String::new();
        for i in 0..10usize {
            let len = [50, 100, 150][i % 3];
            let seq: String = std::iter::repeat_n(b"ACGT"[i % 4] as char, len)
                .collect();
            let qual = "I".repeat(len);
            data.push_str(&format!("@varlen_{i} len={len}\n{seq}\n+\n{qual}\n"));
        }
        fs::write(&fq, data.as_bytes()).unwrap();

        let wdir = sub.join("wdir");
        fs::create_dir_all(&wdir).unwrap();
        let arc = sub.join("out.qz");
        let cfg = CompressConfig {
            input: vec![fq.clone()],
            output: arc.clone(),
            working_dir: wdir.clone(),
            threads: 1,
            cluster: Some(ClusterOptions::default()),
            force: true,
            ..CompressConfig::default()
        };
        qz_lib::compression::compress(&cfg).unwrap();

        let out_fq = sub.join("back.fastq");
        let dc = DecompressConfig {
            input: arc.clone(),
            output: vec![out_fq.clone()],
            working_dir: wdir.clone(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 6,
            force: true,
        };
        qz_lib::compression::decompress(&dc).unwrap();

        assert_eq!(multiset(&fq), multiset(&out_fq), "varlen: set-lossless");
        assert!(verify_ok(&arc, &wdir), "varlen: verify ok");
    }
}

/// Test 3: RC-flip path round-trips losslessly.
///
/// `canonical_minimizer_orientation` flips reads whose dominant minimizer k-mers
/// have a smaller RC hash than the forward hash. G=2, T=3 have large 2-bit values
/// while their complements C=1, A=0 are small → a GT-dominant k-mer has a high
/// forward hash and a low RC hash → the RC is "canonical" → flip=true.
/// We use GT-alternating reads to reliably trigger the RC path.
/// The multiset round-trip assertion verifies that flipped reads come back in their
/// ORIGINAL orientation.
#[test]
fn reverse_complement_path_roundtrips() {
    let temp = TempDir::new().unwrap();
    let tp = temp.path();
    let sub = tp.join("rc_flip");
    fs::create_dir_all(&sub).unwrap();

    let fq = sub.join("in.fastq");
    let mut data = String::new();

    // GT-alternating reads of length 100 → all k-mers will be GT-heavy →
    // forward hashes will be large → RC (CA-alternating) has small hashes →
    // canonical = RC → flip = true.
    let gt_base = "GT".repeat(50); // 100 bp, alternating G/T
    for i in 0..100usize {
        // Vary slightly by substituting one position to avoid all-identical
        let mut seq: Vec<u8> = gt_base.as_bytes().to_vec();
        let idx = i % seq.len();
        seq[idx] = b"GT"[i % 2]; // still G or T, keeps the property
        let seq_str = String::from_utf8(seq).unwrap();
        let qual = "I".repeat(seq_str.len());
        data.push_str(&format!("@rcread_{i} rc_test\n{seq_str}\n+\n{qual}\n"));
    }

    // Also include some explicitly all-G reads (guaranteed large forward hashes)
    for i in 0..50usize {
        let seq = "G".repeat(100);
        let qual = "H".repeat(100);
        data.push_str(&format!("@gcread_{i} all_g\n{seq}\n+\n{qual}\n"));
    }

    // And some AC-rich reads that should NOT flip (their forward hash is already small)
    let ac_base = "AC".repeat(50); // 100 bp, alternating A/C
    for i in 0..50usize {
        let mut seq: Vec<u8> = ac_base.as_bytes().to_vec();
        let idx = i % seq.len();
        seq[idx] = b"AC"[i % 2];
        let seq_str = String::from_utf8(seq).unwrap();
        let qual = "J".repeat(seq_str.len());
        data.push_str(&format!("@acread_{i} ac_test\n{seq_str}\n+\n{qual}\n"));
    }

    fs::write(&fq, data.as_bytes()).unwrap();

    let wdir = sub.join("wdir");
    fs::create_dir_all(&wdir).unwrap();
    let arc = sub.join("out.qz");
    let cfg = CompressConfig {
        input: vec![fq.clone()],
        output: arc.clone(),
        working_dir: wdir.clone(),
        threads: 2,
        cluster: Some(ClusterOptions::default()),
        force: true,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();

    let out_fq = sub.join("back.fastq");
    let dc = DecompressConfig {
        input: arc.clone(),
        output: vec![out_fq.clone()],
        working_dir: wdir.clone(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&dc).unwrap();

    let ms_in = multiset(&fq);
    let ms_back = multiset(&out_fq);
    assert_eq!(
        ms_in.len(),
        ms_back.len(),
        "RC-flip round-trip: read count preserved"
    );
    assert_eq!(
        ms_in, ms_back,
        "RC-flip round-trip must be set-lossless (flipped reads restored to original orientation)"
    );
    assert!(verify_ok(&arc, &wdir), "RC-flip: verify ok");
}

/// Test 4: real-data smoke test (gated — skips gracefully if file absent).
///
/// If NA12878_exome_50k.fastq is present, compress it with cluster mode,
/// decompress, and assert set-lossless round-trip, printing the compression ratio.
#[test]
fn real_data_smoke() {
    let manifest = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let real_data = manifest.join("../../real_data/NA12878_exome_50k.fastq");

    if !real_data.exists() {
        eprintln!("skipping real_data smoke — file absent: {}", real_data.display());
        return;
    }

    let temp = TempDir::new().unwrap();
    let tp = temp.path();
    let wdir = tp.join("wdir");
    fs::create_dir_all(&wdir).unwrap();

    let arc = tp.join("smoke.qz");
    let cfg = CompressConfig {
        input: vec![real_data.clone()],
        output: arc.clone(),
        working_dir: wdir.clone(),
        threads: 4,
        cluster: Some(ClusterOptions::default()),
        force: true,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();

    let input_len = fs::metadata(&real_data).unwrap().len();
    let archive_len = fs::metadata(&arc).unwrap().len();
    let ratio = input_len as f64 / archive_len as f64;
    eprintln!("real_data_smoke: input={input_len} archive={archive_len} ratio={ratio:.2}x");

    let out_fq = tp.join("smoke.fastq");
    let dc = DecompressConfig {
        input: arc.clone(),
        output: vec![out_fq.clone()],
        working_dir: wdir.clone(),
        num_threads: 2,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&dc).unwrap();

    let ms_in = multiset(&real_data);
    let ms_back = multiset(&out_fq);
    assert_eq!(
        ms_in.len(),
        ms_back.len(),
        "real_data: read count preserved"
    );
    assert_eq!(ms_in, ms_back, "real_data: set-lossless round-trip");

    let verified = verify_ok(&arc, &wdir);
    eprintln!("real_data_smoke: verify={verified}");
    assert!(verified, "real_data: verify must pass");
}
