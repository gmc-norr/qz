/// End-to-end integration tests for the cluster paired-end compression mode.
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

/// Write paired FASTQ files `r1` and `r2` with `n` pairs, each of `len` bp.
///
/// Pair `i` gets headers `@pair{i}/1` (R1) and `@pair{i}/2` (R2).
/// Sequences are deterministic but DISTINCT between mates:
///   R1 seq: alphabet[(i*37 + j*13) % 4] per position j
///   R2 seq: alphabet[(i*31 + j*17 + 2) % 4] per position j  (offset by 2 = different)
/// Qualities are similarly mate-distinct.
fn write_pair_fastq(r1: &Path, r2: &Path, n: usize, len: usize) {
    let alphabet = b"ACGT";
    let mut buf1 = String::new();
    let mut buf2 = String::new();

    for i in 0..n {
        // R1
        let seq1: String = (0..len)
            .map(|j| alphabet[(i.wrapping_mul(37).wrapping_add(j.wrapping_mul(13))) % 4] as char)
            .collect();
        let qual1: String = (0..len as u8)
            .map(|j| (b'!' + ((i as u8).wrapping_add(j) % 40)) as char)
            .collect();
        buf1.push_str(&format!("@pair{i}/1\n{seq1}\n+\n{qual1}\n"));

        // R2 (distinct: different multipliers + offset)
        let seq2: String = (0..len)
            .map(|j| alphabet[(i.wrapping_mul(31).wrapping_add(j.wrapping_mul(17)).wrapping_add(2)) % 4] as char)
            .collect();
        let qual2: String = (0..len as u8)
            .map(|j| (b'!' + ((i as u8).wrapping_add(j).wrapping_add(5) % 40)) as char)
            .collect();
        buf2.push_str(&format!("@pair{i}/2\n{seq2}\n+\n{qual2}\n"));
    }

    fs::write(r1, buf1.as_bytes()).unwrap();
    fs::write(r2, buf2.as_bytes()).unwrap();
}

/// Parse `path` as a FASTQ and return a SORTED multiset of (header, seq, qual) byte
/// tuples. Sorting makes comparison order-independent (set-lossless). Because cluster
/// has no byte-identical golden, this helper IS the safety net — so it parses
/// STRICTLY: the file must be a whole number of well-formed 4-line records (`@`
/// header, `+` separator, seq.len() == qual.len()), not a lenient line scan that
/// would silently absorb a line-structure corruption.
fn multiset(path: &Path) -> Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> {
    let data = fs::read(path).unwrap();
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

/// Check that read i of `o1` and read i of `o2` share their pair{N} token.
///
/// For each line index (read i), parse the header to find the `pair{N}` token
/// in the prefix (before the `/` mate suffix) and verify both mates agree.
fn pairs_aligned(o1: &Path, o2: &Path) -> bool {
    /// Extract the pair index N from a header like `@pair{N}/1` or `@pair{N}/2`.
    fn pair_token(hdr: &[u8]) -> Option<Vec<u8>> {
        // Strip leading '@'
        let hdr = if hdr.first() == Some(&b'@') { &hdr[1..] } else { hdr };
        // Find the last '/' and strip the mate suffix (e.g. "/1" or "/2")
        let base = if let Some(pos) = hdr.iter().rposition(|&b| b == b'/') {
            &hdr[..pos]
        } else {
            hdr
        };
        // Strip any description (space-separated)
        let name = if let Some(pos) = base.iter().position(|&b| b == b' ') {
            &base[..pos]
        } else {
            base
        };
        Some(name.to_vec())
    }

    let data1 = fs::read(o1).unwrap();
    let data2 = fs::read(o2).unwrap();

    let headers1: Vec<&[u8]> = data1
        .split(|&b| b == b'\n')
        .enumerate()
        .filter(|(i, line)| i % 4 == 0 && !line.is_empty())
        .map(|(_, line)| line)
        .collect();
    let headers2: Vec<&[u8]> = data2
        .split(|&b| b == b'\n')
        .enumerate()
        .filter(|(i, line)| i % 4 == 0 && !line.is_empty())
        .map(|(_, line)| line)
        .collect();

    if headers1.len() != headers2.len() {
        return false;
    }

    for (h1, h2) in headers1.iter().zip(headers2.iter()) {
        let t1 = pair_token(h1);
        let t2 = pair_token(h2);
        if t1 != t2 {
            return false;
        }
    }
    true
}

/// Build a CompressConfig with paired cluster mode, compress, and return the archive path.
///
/// Each call gets its own fresh working_dir subdir named `work_{name}` to avoid
/// scratch-file collisions between sequential calls on the same base directory.
fn compress_paired_to(base: &Path, r1: &Path, r2: &Path, name: &str) -> PathBuf {
    let work_name = format!("work_{}", name.replace('.', "_"));
    let sub = base.join(&work_name);
    fs::create_dir_all(&sub).unwrap();
    let output = base.join(name);
    let cfg = CompressConfig {
        input: vec![r1.to_path_buf(), r2.to_path_buf()],
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

/// Build a DecompressConfig for paired output (two output paths) and decompress.
///
/// Returns `(o1_path, o2_path)`.
fn decompress_paired_to(
    base: &Path,
    archive: &Path,
    o1name: &str,
    o2name: &str,
) -> (PathBuf, PathBuf) {
    let o1 = base.join(o1name);
    let o2 = base.join(o2name);
    let dc = DecompressConfig {
        input: archive.to_path_buf(),
        output: vec![o1.clone(), o2.clone()],
        working_dir: base.to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&dc).unwrap();
    (o1, o2)
}

// ─── tests ────────────────────────────────────────────────────────────────────

/// Test 1: paired roundtrip — set-lossless, pairing intact, byte-identical determinism.
///
/// Compresses the same R1/R2 pair twice (separate working dirs) and asserts:
/// - byte-identical archives (determinism)
/// - archive_type byte == 3 (cluster)
/// - per-mate multiset preserved in each decompressed output
/// - pairing intact: read i of R1-out and read i of R2-out share the same pair token
/// - verify passes
#[test]
fn paired_roundtrip_pairing_and_determinism() {
    let temp = TempDir::new().unwrap();
    let tp = temp.path();

    let r1 = tp.join("r1.fastq");
    let r2 = tp.join("r2.fastq");
    write_pair_fastq(&r1, &r2, 2000, 100);

    // Compress twice with separate working dirs — archives must be byte-identical.
    let a = compress_paired_to(tp, &r1, &r2, "a.qz");
    let b = compress_paired_to(tp, &r1, &r2, "b.qz");

    let bytes_a = fs::read(&a).unwrap();
    let bytes_b = fs::read(&b).unwrap();
    assert_eq!(
        bytes_a, bytes_b,
        "cluster paired archive must be byte-identical for identical input (determinism)"
    );

    // Byte 3 must be archive_type 3 (cluster)
    assert_eq!(
        bytes_a[3], 3,
        "archive_type byte must be 3 for cluster archive"
    );

    // Decompress -> two outputs
    let wdir_dec = tp.join("dec_wdir");
    fs::create_dir_all(&wdir_dec).unwrap();
    let (o1, o2) = decompress_paired_to(&wdir_dec, &a, "o1.fastq", "o2.fastq");

    // Per-mate multiset round-trip
    assert_eq!(
        multiset(&r1),
        multiset(&o1),
        "R1 multiset must be preserved after paired cluster round-trip"
    );
    assert_eq!(
        multiset(&r2),
        multiset(&o2),
        "R2 multiset must be preserved after paired cluster round-trip"
    );

    // Pairing must be intact: read i of o1 and read i of o2 share the same pair token
    assert!(
        pairs_aligned(&o1, &o2),
        "output pairs must be aligned: read i of R1 and read i of R2 share the same pair token"
    );

    // Verify integrity
    let vc = VerifyConfig {
        input: a,
        working_dir: wdir_dec,
        num_threads: 1,
        fast: false,
    };
    assert!(
        qz_lib::compression::verify(&vc).is_ok(),
        "verify must pass on paired cluster archive"
    );
}

/// Test 2: real-data smoke test (deep-coverage ratio validation).
///
/// If HG002 chr20 NovaSeq R1/R2 files are present, compress with cluster mode,
/// decompress, and assert set-lossless round-trip per mate.
///
/// `#[ignore]` by default: the chr20 dataset is ~10.8M pairs (3.7 GB/mate), so a
/// debug-build round-trip takes tens of minutes — too heavy for routine `cargo test`.
/// The synthetic `paired_roundtrip_pairing_and_determinism` is the fast per-run gate;
/// run this deep-coverage validation on demand with
/// `cargo test -p qz-lib --release --test cluster_paired_test -- --ignored`.
#[test]
#[ignore = "heavy: 10.8M-pair chr20 round-trip; run on demand with --ignored"]
fn paired_real_data_smoke() {
    let base = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../real_data");
    let r1 = base.join("HG002.chr20.novaseq.2x151.R1.fastq");
    let r2 = base.join("HG002.chr20.novaseq.2x151.R2.fastq");

    if !r1.exists() || !r2.exists() {
        eprintln!(
            "skipping paired real_data smoke — files absent: {} / {}",
            r1.display(),
            r2.display()
        );
        return;
    }

    let temp = TempDir::new().unwrap();
    let tp = temp.path();

    let a = compress_paired_to(tp, &r1, &r2, "chr20.qz");

    let insz =
        fs::metadata(&r1).unwrap().len() + fs::metadata(&r2).unwrap().len();
    let outsz = fs::metadata(&a).unwrap().len();
    eprintln!(
        "paired chr20 smoke: in={insz} archive={outsz} ratio={:.2}x",
        insz as f64 / outsz as f64
    );

    let dec_wdir = tp.join("dec_wdir");
    fs::create_dir_all(&dec_wdir).unwrap();
    let (o1, o2) = decompress_paired_to(&dec_wdir, &a, "o1.fastq", "o2.fastq");

    assert_eq!(multiset(&r1), multiset(&o1), "chr20 R1: set-lossless");
    assert_eq!(multiset(&r2), multiset(&o2), "chr20 R2: set-lossless");
    assert!(pairs_aligned(&o1, &o2), "chr20: pairing intact");

    let vc = VerifyConfig {
        input: a,
        working_dir: dec_wdir,
        num_threads: 1,
        fast: false,
    };
    assert!(
        qz_lib::compression::verify(&vc).is_ok(),
        "chr20: verify must pass"
    );
}
