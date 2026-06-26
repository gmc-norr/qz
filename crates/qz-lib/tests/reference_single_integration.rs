use qz_lib::cli::{CompressConfig, DecompressConfig, QualityMode, ReferenceOptions, VerifyConfig};
use qz_lib::compression::{compress, decompress, verify, VerifyMode};
use std::path::PathBuf;

// --- CANONICAL shared helpers (defined ONCE here, Increment 3) ---
// Increments 4 and 5 EXTEND this file and REUSE these helpers; they do NOT
// redefine `w` / `make_seq` / `make_synthetic_single` / `cfg_ref_single` /
// `assert_roundtrip_single`, and they do NOT recreate this file.

fn w(p: &std::path::Path, s: &str) {
    std::fs::write(p, s).unwrap();
}

/// Deterministic low-repeat ACGT sequence (matches mapping.rs's generator so
/// reads seed reliably). Identical to the paired suite's `make_seq`.
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

/// Build a synthetic SINGLE-END dataset against a deterministic reference.
/// Returns `(ref_fasta, reads_path, reads_fastq_string)` where the string is the
/// EXACT bytes a byte-identical decompress must produce. ~2 kb reference; 200
/// reads as 100-120 bp forward substrings, a handful with 1-2 planted
/// substitutions, and the last two fully random => literal fallback.
fn make_synthetic_single(dir: &std::path::Path) -> (PathBuf, PathBuf, String) {
    let refseq = make_seq(2000, 7);
    let rf = dir.join("ref.fa");
    {
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        w(&rf, &s);
    }

    let n_reads = 200usize;
    let mut s1 = String::new();
    let mut rng = 0x1234_5678_9abc_def0u64;
    let mut next = || {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (rng >> 33) as usize
    };

    for i in 0..n_reads {
        let unmappable = i >= n_reads - 2; // last two: random => fallback
        let rlen = 100 + (next() % 21); // 100..=120
        let rbytes = if unmappable {
            make_seq(rlen, 9000 + i as u64)
        } else {
            let max_start = refseq.len() - rlen;
            let st1 = next() % (max_start + 1);
            let mut r1 = refseq[st1..st1 + rlen].to_vec();
            if i % 17 == 0 {
                let p = next() % rlen;
                r1[p] = match r1[p] {
                    b'A' => b'C',
                    _ => b'A',
                };
            }
            r1
        };
        let q: String = "I".repeat(rbytes.len());
        s1.push_str(&format!(
            "@read_{i}\n{}\n+\n{q}\n",
            std::str::from_utf8(&rbytes).unwrap()
        ));
    }

    let r1p = dir.join("reads.fastq");
    w(&r1p, &s1);
    (rf, r1p, s1)
}

/// Single-end reference compress config: ONE input file + the reference block.
/// Only the fields differing from `CompressConfig::default()` are set; the rest
/// flow in via the spread, so a new CompressConfig field never touches this helper.
/// The `ReferenceOptions` field set is copied verbatim from `cfg_ref` in
/// `tests/reference_integration.rs`.
fn cfg_ref_single(input: PathBuf, out: PathBuf, tmp: PathBuf, refp: PathBuf) -> CompressConfig {
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

/// Decode the type-4 `archive` to ONE output file at `tmp/name` and assert that
/// verbatim `output[0]` equals `expected` byte-for-byte. Single-end reference
/// writes one bare file (no `_R1` suffix), exactly like single-end default.
/// `name` lets each test pick a distinct output filename. Arg order is
/// `(tmp, archive, name, expected)` (the order Increment 5's call sites use).
fn assert_roundtrip_single(
    tmp: &std::path::Path,
    archive: PathBuf,
    name: &str,
    expected: &str,
) {
    let decp = tmp.join(name);
    let dc = DecompressConfig {
        input: archive,
        output: vec![decp.clone()], // the verbatim output[0] — ONE file, no suffix
        working_dir: tmp.to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    decompress(&dc).unwrap();
    assert_eq!(
        std::fs::read_to_string(&decp).unwrap(),
        expected,
        "single-end reference roundtrip not byte-identical at output[0]"
    );
}

#[test]
fn single_ref_verify_deep_passes() {
    let dir = tempfile::tempdir().unwrap();
    let (ref_fa, reads_fq, _s1) = make_synthetic_single(dir.path());
    let archive = dir.path().join("out.qz");
    compress(&cfg_ref_single(
        reads_fq,
        archive.clone(),
        dir.path().to_path_buf(),
        ref_fa,
    ))
    .unwrap();

    let vc = VerifyConfig { input: archive.clone(), fast: false, ..VerifyConfig::default() };
    let res = verify(&vc).unwrap();
    assert_eq!(res.mode, VerifyMode::Deep);
    assert_eq!(res.num_reads, 200, "single-end: num_reads is NOT doubled");
    assert!(res.crc32 != 0, "deep verify computes a CRC32 over the one FASTQ");
    assert_eq!(res.r2_crc32, None, "single-end reference has no R2 mate");
    assert!(res.total_bytes > 0);
}

#[test]
fn single_ref_verify_fast_passes() {
    let dir = tempfile::tempdir().unwrap();
    let (ref_fa, reads_fq, _s1) = make_synthetic_single(dir.path());
    let archive = dir.path().join("out.qz");
    compress(&cfg_ref_single(
        reads_fq,
        archive.clone(),
        dir.path().to_path_buf(),
        ref_fa,
    ))
    .unwrap();

    let vc = VerifyConfig { input: archive.clone(), fast: true, ..VerifyConfig::default() };
    let res = verify(&vc).unwrap();
    assert_eq!(res.mode, VerifyMode::Fast);
    assert_eq!(res.num_reads, 200, "single-end: num_reads is NOT doubled");
    assert!(res.blocks_verified > 0, "fast verify CRC-walks every block frame");
    assert_eq!(res.crc32, 0, "fast verify does not reconstruct, so no FASTQ CRC");
    assert_eq!(res.r2_crc32, None);
    assert!(res.total_bytes > 0, "total_bytes = sum of compressed entry bytes walked");
}

#[test]
fn single_end_reference_dispatch_roundtrip() {
    let d = tempfile::tempdir().unwrap();
    let (ref_fa, reads_fq, s1) = make_synthetic_single(d.path());
    let out = d.path().join("o.qz");

    compress(&cfg_ref_single(
        reads_fq,
        out.clone(),
        d.path().to_path_buf(),
        ref_fa,
    ))
    .unwrap();

    // archive_type 4 + classifiers: single-end reference IS its own predicate,
    // and is NOT the (two-output) `is_reference_archive` — so it stays decodable
    // through the qz-python single-output API.
    let b = std::fs::read(&out).unwrap();
    assert_eq!(b[3], 4, "single-end reference archive_type");
    assert!(qz_lib::compression::is_reference_single_archive(&out).unwrap());
    assert!(!qz_lib::compression::is_reference_archive(&out).unwrap());

    // Top-level decompress dispatch routes Some(4) to one FASTQ (verbatim output[0]).
    assert_roundtrip_single(d.path(), out, "dec.fastq", &s1);
}

// ============================================================================
// Increment 5: Integration + hostile/structural suite for single-end reference
// (archive_type 4). Mirrors reference_integration.rs (paired suite) but with
// ONE input and ONE output. The paired type-2 path is untouched.
// ============================================================================

// ---- Step 1: Smoke test — verifies archive_type byte and classifiers ----------

// Single-end reference-based mode (archive_type 4). MIRRORS reference_integration.rs
// (the paired suite) but with ONE input and ONE output. The paired type-2 path is
// untouched; this file never asserts against it.
#[test]
fn single_reference_compress_writes_type4_archive() {
    let d = tempfile::tempdir().unwrap();
    let refseq = make_seq(600, 42);
    let rf = d.path().join("ref.fa");
    {
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        w(&rf, &s);
    }
    // A few exact-substring 120 bp reads.
    let mut s1 = String::new();
    for (i, &st) in [0usize, 100, 240, 360, 400].iter().enumerate() {
        let r1 = &refseq[st..st + 120];
        let q: String = "I".repeat(120);
        s1.push_str(&format!(
            "@read_{i}\n{}\n+\n{q}\n",
            std::str::from_utf8(r1).unwrap()
        ));
    }
    let r1p = d.path().join("r1.fastq");
    w(&r1p, &s1);

    let out = d.path().join("o.qz");
    let c = cfg_ref_single(r1p, out.clone(), d.path().to_path_buf(), rf);
    qz_lib::compression::compress(&c).expect("single-end reference compress");

    let b = std::fs::read(&out).expect("archive exists");
    assert!(b.len() > 18, "archive too small: {}", b.len());
    assert_eq!(&b[0..2], b"QZ", "magic");
    assert_eq!(b[2], 5, "v5 chunk-major version");
    assert_eq!(b[3], 4, "single-end reference archive_type");
    // FIX 3: `is_reference_archive` stays type-2-only — a single-end (type-4) archive
    // is its OWN classifier and is NOT a (two-output) reference archive, so it remains
    // decodable through the qz-python single-output API.
    assert!(
        qz_lib::compression::is_reference_single_archive(&out)
            .expect("is_reference_single_archive"),
        "type-4 archive must classify as single-end reference"
    );
    assert!(
        !qz_lib::compression::is_reference_archive(&out).expect("is_reference_archive"),
        "type-4 (single-output) archive must NOT classify as (two-output) reference"
    );
}

// ---- Step 3: Byte-identical roundtrip test (multichunk) ----------------------

#[test]
fn single_reference_roundtrip_byte_identical() {
    let d = tempfile::tempdir().unwrap();
    let (refp, r1p, s1) = make_synthetic_single(d.path());
    let out = d.path().join("s.qz");
    let mut c = cfg_ref_single(r1p, out.clone(), d.path().to_path_buf(), refp);
    c.advanced.chunk_records = 64; // several chunks
    qz_lib::compression::compress(&c).unwrap();
    assert_roundtrip_single(d.path(), out, "dec.fastq", &s1);
}

// ---- Step 4: Fast-seeding roundtrip ------------------------------------------

#[test]
fn single_reference_fast_roundtrip_byte_identical() {
    // --reference-fast (sparser seeds) must still be lossless: the reconstruct-
    // verify falls any unmappable read back to a literal regardless of seeding.
    let d = tempfile::tempdir().unwrap();
    let (refp, r1p, s1) = make_synthetic_single(d.path());
    let out = d.path().join("sf.qz");
    let mut c = cfg_ref_single(r1p, out.clone(), d.path().to_path_buf(), refp);
    c.advanced.chunk_records = 64;
    c.reference.as_mut().unwrap().reference_fast = true;
    qz_lib::compression::compress(&c).unwrap();
    assert_roundtrip_single(d.path(), out, "decf.fastq", &s1);
}

// ---- Step 6: Byte-level v5 footer view + role/codec constants ----------------
// (placed before Step 5's tests so the helpers exist when those tests compile)

const ENTRY_BYTES: usize = 31; // chunk_index(4)+role(1)+mate(1)+codec(1)+off(8)+len(8)+rc(8)
const V5_LOCATOR_LEN: usize = 20;
const V5_FOOTER_PREFIX: usize = 8 + 4 + 4; // num_reads + num_chunks + n_entries

#[derive(Clone, Copy)]
struct FEntry {
    abs: usize,
    role: u8,
    offset: u64,
    length: u64,
}
struct FView {
    footer_offset: usize,
    entry_count_abs: usize,
    entries: Vec<FEntry>,
}
fn read_u64(b: &[u8], o: usize) -> u64 {
    u64::from_le_bytes(b[o..o + 8].try_into().unwrap())
}
fn read_u32(b: &[u8], o: usize) -> u32 {
    u32::from_le_bytes(b[o..o + 4].try_into().unwrap())
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
    let entry_count_abs = footer_offset + 8 + 4;
    let n = read_u32(bytes, entry_count_abs) as usize;
    let mut entries = Vec::with_capacity(n);
    let mut o = footer_offset + V5_FOOTER_PREFIX;
    for _ in 0..n {
        let role = bytes[o + 4];
        let offset = read_u64(bytes, o + 7);
        let length = read_u64(bytes, o + 15);
        entries.push(FEntry { abs: o, role, offset, length });
        o += ENTRY_BYTES;
    }
    FView { footer_offset, entry_count_abs, entries }
}
fn find_role(v: &FView, role: u8) -> FEntry {
    *v.entries.iter().find(|e| e.role == role).expect("role present")
}
fn footer_num_chunks(bytes: &[u8]) -> u32 {
    let v = fview(bytes);
    read_u32(bytes, v.footer_offset + 8)
}
fn sum_perchunk_record_count(bytes: &[u8], role: u8) -> u64 {
    let v = fview(bytes);
    v.entries
        .iter()
        .filter(|e| e.role == role && bytes[e.abs..e.abs + 4] != [0xFF, 0xFF, 0xFF, 0xFF])
        .map(|e| read_u64(bytes, e.abs + 23))
        .sum()
}

// Role byte codes — unified StreamRole discriminants (chunk_directory.rs).
// Some constants are included for completeness / future tests; suppress the lint.
#[allow(dead_code)]
const R_HEADERS: u8 = 0; // StreamRole::Headers
const R_QUAL: u8 = 3; // StreamRole::Qual
const R_CONSENSUS: u8 = 6; // StreamRole::PackedBacking
const R_INTERVALMAP: u8 = 7; // StreamRole::IntervalMap
#[allow(dead_code)]
const R_NBITMAP: u8 = 8; // StreamRole::NBitmap
#[allow(dead_code)]
const R_REFMETA: u8 = 9; // StreamRole::ReferenceMeta
const R_MAPPEDFLAGS: u8 = 10; // StreamRole::MappedFlags
const R_POSITIONS: u8 = 11; // StreamRole::Positions
const R_EDITCOUNT: u8 = 14; // StreamRole::EditCount
#[allow(dead_code)]
const R_EDITPOS: u8 = 15; // StreamRole::EditPos
const R_EDITBASE: u8 = 16; // StreamRole::EditBase
const R_FALLBACKPOOL: u8 = 17; // StreamRole::FallbackPool
const CODEC_FQZCOMP: u8 = 6; // codec_ids.rs (the only quality codec reference emits)

/// Build a valid single-end reference archive once. `chunk_records = 64` forces
/// several chunks; the synthetic dataset's last two reads are unmappable => fallback.
fn valid_archive(d: &std::path::Path) -> (Vec<u8>, PathBuf, String) {
    let (refp, r1p, s1) = make_synthetic_single(d);
    let out = d.join("valid.qz");
    let mut c = cfg_ref_single(r1p, out.clone(), d.to_path_buf(), refp);
    c.advanced.chunk_records = 64;
    qz_lib::compression::compress(&c).unwrap();
    let bytes = std::fs::read(&out).unwrap();
    (bytes, out, s1)
}

#[test]
fn single_reference_fallback_path_is_exercised() {
    let d = tempfile::tempdir().unwrap();
    let (bytes, path, s1) = valid_archive(d.path());
    let v = fview(&bytes);

    // Single-end invariant: NO entry carries mate == 2 (mate at abs+5).
    assert!(
        v.entries.iter().all(|e| bytes[e.abs + 5] != 2),
        "single-end reference must not emit any mate-2 directory entry"
    );

    let mut fallback: u64 = 0;
    for e in &v.entries {
        if e.role == R_FALLBACKPOOL {
            fallback += read_u64(&bytes, e.abs + 23);
        }
    }
    assert!(fallback > 0, "dataset must exercise the fallback path (got {fallback})");

    assert_roundtrip_single(d.path(), path, "ok.fastq", &s1);
}

// ---- Step 5: All-mapped and all-fallback roundtrip tests ---------------------

#[test]
fn single_reference_all_mapped_roundtrip() {
    // Every read is an exact forward substring => all mapped, zero fallback.
    let d = tempfile::tempdir().unwrap();
    let refseq = make_seq(2000, 7);
    let rf = d.path().join("ref.fa");
    {
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        w(&rf, &s);
    }
    let mut s1 = String::new();
    for i in 0..120usize {
        let st = (i * 13) % (refseq.len() - 130);
        let r1 = &refseq[st..st + 120];
        let q: String = "I".repeat(120);
        s1.push_str(&format!("@m_{i}\n{}\n+\n{q}\n", std::str::from_utf8(r1).unwrap()));
    }
    let r1p = d.path().join("reads.fastq");
    w(&r1p, &s1);
    let out = d.path().join("am.qz");
    let mut c = cfg_ref_single(r1p, out.clone(), d.path().to_path_buf(), rf);
    c.advanced.chunk_records = 40; // ≥3 chunks of mapped reads
    qz_lib::compression::compress(&c).unwrap();

    // All-mapped: Σ per-chunk Positions record_count == total reads (every read mapped).
    let bytes = std::fs::read(&out).unwrap();
    assert_eq!(
        sum_perchunk_record_count(&bytes, R_POSITIONS),
        120,
        "all-mapped dataset must have M == N"
    );
    assert_roundtrip_single(d.path(), out, "am.fastq", &s1);
}

#[test]
fn single_reference_all_fallback_zero_mapped_roundtrip() {
    // Every read is random ACGT unrelated to the reference => literal fallback (M==0).
    let d = tempfile::tempdir().unwrap();
    let refseq = make_seq(2000, 7);
    let rf = d.path().join("ref.fa");
    {
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        w(&rf, &s);
    }
    let n_reads = 100usize;
    let mut s1 = String::new();
    for i in 0..n_reads {
        let r1 = make_seq(110, 0xA000_0000 + i as u64);
        let q: String = "I".repeat(110);
        s1.push_str(&format!("@frag_{i}\n{}\n+\n{q}\n", std::str::from_utf8(&r1).unwrap()));
    }
    let r1p = d.path().join("reads.fastq");
    w(&r1p, &s1);
    let out = d.path().join("allfb.qz");
    let mut c = cfg_ref_single(r1p, out.clone(), d.path().to_path_buf(), rf);
    c.advanced.chunk_records = 16; // ≥6 chunks
    qz_lib::compression::compress(&c).unwrap();

    let bytes = std::fs::read(&out).unwrap();
    assert!(footer_num_chunks(&bytes) >= 2, "must span multiple chunks");
    assert_eq!(
        sum_perchunk_record_count(&bytes, R_POSITIONS),
        0,
        "all-fallback dataset must have M==0 (no mapped reads)"
    );
    // MappedFlags record_count is N_chunk (all reads), here all unmapped.
    assert_eq!(
        sum_perchunk_record_count(&bytes, R_MAPPEDFLAGS),
        n_reads as u64,
        "MappedFlags record_count is N_chunk (all reads)"
    );
    assert_eq!(
        sum_perchunk_record_count(&bytes, R_FALLBACKPOOL),
        n_reads as u64,
        "per-(chunk,mate) fallback pools must hold every read"
    );
    assert_roundtrip_single(d.path(), out, "allfb.fastq", &s1);
}

// ---- Step 7: Multichunk (≥3) + verify deep + fast ----------------------------

#[test]
fn single_reference_multichunk_three_plus_verify() {
    let d = tempfile::tempdir().unwrap();
    let (refp, r1p, s1) = make_synthetic_single(d.path());
    let out = d.path().join("mc3.qz");
    let mut c = cfg_ref_single(r1p, out.clone(), d.path().to_path_buf(), refp);
    c.advanced.chunk_records = 50; // 200 reads / 50 => 4 chunks (≥3)
    qz_lib::compression::compress(&c).unwrap();

    let bytes = std::fs::read(&out).unwrap();
    assert!(
        footer_num_chunks(&bytes) >= 3,
        "must span ≥3 chunks (got {})",
        footer_num_chunks(&bytes)
    );

    // Roundtrip + capture the emitted bytes for the CRC contract.
    let outp = d.path().join("mc3.fastq");
    let dc = qz_lib::cli::DecompressConfig {
        input: out.clone(),
        output: vec![outp.clone()],
        working_dir: d.path().to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&dc).unwrap();
    assert_eq!(std::fs::read_to_string(&outp).unwrap(), s1, "not byte-identical");

    // Deep verify: num_reads == 200 (single-end, no ×2), r2_crc32 == None, and the
    // reported crc32 equals CRC32 of the single emitted file (deep verify CRCs the
    // reconstructed stream, which is byte-for-byte what was written).
    let deep = qz_lib::compression::verify(&qz_lib::cli::VerifyConfig {
        input: out.clone(),
        working_dir: d.path().to_path_buf(),
        num_threads: 1,
        fast: false,
    })
    .unwrap();
    assert_eq!(deep.num_reads, 200, "single-end read count (not ×2)");
    assert!(deep.r2_crc32.is_none(), "single-end deep verify must report no R2 CRC");
    let file_crc = v5_footer_crc(&std::fs::read(&outp).unwrap());
    assert_eq!(deep.crc32, file_crc, "deep-verify crc32 must equal CRC32 of the emitted FASTQ");

    // Fast verify walks per-block CRCs.
    let fast = qz_lib::compression::verify(&qz_lib::cli::VerifyConfig {
        input: out,
        working_dir: d.path().to_path_buf(),
        num_threads: 1,
        fast: true,
    })
    .unwrap();
    assert!(fast.blocks_verified > 0);
    assert_eq!(fast.crc32, 0, "fast verify does not reconstruct, so no FASTQ CRC");
    assert!(fast.r2_crc32.is_none(), "single-end has no R2");
}

// ---- Step 8: Quality-is-fqzcomp + decompress_to_writer-rejects ---------------

#[test]
fn single_reference_quality_is_fqzcomp_roundtrip() {
    let d = tempfile::tempdir().unwrap();
    let (refp, r1p, s1) = make_synthetic_single(d.path());
    let out = d.path().join("fqz.qz");
    let mut c = cfg_ref_single(r1p, out.clone(), d.path().to_path_buf(), refp);
    c.advanced.chunk_records = 64;
    qz_lib::compression::compress(&c).unwrap();

    let bytes = std::fs::read(&out).unwrap();
    let v = fview(&bytes);
    let qual_codecs: Vec<u8> = v
        .entries
        .iter()
        .filter(|e| e.role == R_QUAL)
        .map(|e| bytes[e.abs + 6])
        .collect();
    assert!(!qual_codecs.is_empty(), "no quality entries");
    assert!(
        qual_codecs.iter().all(|&x| x == CODEC_FQZCOMP),
        "reference quality entries must all be fqzcomp(6): {qual_codecs:?}"
    );
    assert_roundtrip_single(d.path(), out, "fqz.fastq", &s1);
}

#[test]
fn decompress_to_writer_rejects_single_reference() {
    let d = tempfile::tempdir().unwrap();
    let (refp, r1p, _s1) = make_synthetic_single(d.path());
    let out = d.path().join("s.qz");
    let c = cfg_ref_single(r1p, out.clone(), d.path().to_path_buf(), refp);
    qz_lib::compression::compress(&c).unwrap();
    let mut sink = Vec::new();
    assert!(qz_lib::compression::decompress_to_writer(&out, &mut sink).is_err());
}

// ---- Step 9: Rejected-options matrix -----------------------------------------

#[test]
fn single_reference_rejected_options_matrix() {
    use qz_lib::cli::{HeaderCompressor, QualityCompressor};

    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("reads.fastq");
    w(&r1, "@r\nACGTACGT\n+\nIIIIIIII\n");
    let rf = d.path().join("ref.fa");
    w(&rf, ">c\nACGTACGTACGTACGT\n");

    let mutations: Vec<(&str, Box<dyn Fn(&mut CompressConfig)>)> = vec![
        ("fasta", Box::new(|c: &mut CompressConfig| c.fasta = true)),
        ("no_quality", Box::new(|c: &mut CompressConfig| c.no_quality = true)),
        (
            "lossy_quality",
            Box::new(|c: &mut CompressConfig| c.quality_mode = QualityMode::IlluminaBin),
        ),
        ("ultra", Box::new(|c: &mut CompressConfig| c.ultra = Some(1))),
        (
            "bsc_block_size",
            Box::new(|c: &mut CompressConfig| c.advanced.bsc_block_size_mb = 10),
        ),
        (
            "header_compressor",
            Box::new(|c: &mut CompressConfig| c.advanced.header_compressor = HeaderCompressor::Bsc),
        ),
        (
            "quality_compressor",
            Box::new(|c: &mut CompressConfig| c.advanced.quality_compressor = QualityCompressor::Bsc),
        ),
        ("bsc_static", Box::new(|c: &mut CompressConfig| c.advanced.bsc_static = true)),
        (
            "sequence_hints",
            Box::new(|c: &mut CompressConfig| c.advanced.sequence_hints = true),
        ),
        ("rc_canon", Box::new(|c: &mut CompressConfig| c.advanced.rc_canon = true)),
        (
            "reference_index",
            Box::new(|c: &mut CompressConfig| {
                c.reference.as_mut().unwrap().reference_index = Some(PathBuf::from("idx"))
            }),
        ),
    ];

    for (name, mutate) in &mutations {
        let out = d.path().join(format!("out_{name}.qz"));
        let mut c = cfg_ref_single(r1.clone(), out.clone(), d.path().to_path_buf(), rf.clone());
        mutate(&mut c);
        let res = qz_lib::compression::compress(&c);
        assert!(res.is_err(), "{name} should be rejected but compress succeeded");
        assert!(!out.exists(), "{name}: gate must fire before any output is created");
    }
}

// ---- Step 10: Determinism — parallel matches serial --------------------------

#[test]
fn single_reference_parallel_matches_serial_bytes() {
    let d = tempfile::tempdir().unwrap();
    let (refp, r1p, _s1) = make_synthetic_single(d.path());
    let mk = |out: &std::path::Path, threads: usize| {
        let mut c = cfg_ref_single(
            r1p.clone(),
            out.to_path_buf(),
            d.path().to_path_buf(),
            refp.clone(),
        );
        c.threads = threads;
        c.advanced.chunk_records = 64;
        qz_lib::compression::compress(&c).unwrap();
        std::fs::read(out).unwrap()
    };
    let a = mk(&d.path().join("t1.qz"), 1);
    let b = mk(&d.path().join("t4.qz"), 4);
    let e = mk(&d.path().join("t8.qz"), 8);
    assert_eq!(a, b, "type-4 archive must be byte-identical 1 vs 4 threads");
    assert_eq!(a, e, "type-4 archive must be byte-identical 1 vs 8 threads");
}

// ---- Step 11: Hostile / structural suite -------------------------------------

/// Apply `f` to a COPY of `bytes`, write it, and assert decompress returns Err AND
/// leaves NO output file (single-end => one output, not _R1/_R2). Also asserts
/// verify(fast/deep) Err (never panic).
fn expect_reject(d: &std::path::Path, case: &str, bytes: &[u8], f: impl FnOnce(&mut Vec<u8>)) {
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

#[test]
fn single_hostile_unknown_archive_type() {
    let d = tempfile::tempdir().unwrap();
    let (bytes, _path, _s1) = valid_archive(d.path());
    // byte[3] is archive_type; single-end reference is 4. Use 9 (unknown).
    expect_reject(d.path(), "unknown_archive_type", &bytes, |b| {
        b[3] = 9;
    });
}

#[test]
fn single_hostile_structure_globals_and_perchunk() {
    let d = tempfile::tempdir().unwrap();
    let (bytes, _path, _s1) = valid_archive(d.path());

    // Drop a global (Consensus): overwrite its role with IntervalMap => dup
    // IntervalMap + missing Consensus => validator rejects.
    expect_reject(d.path(), "drop_global", &bytes, |b| {
        let v = fview(b);
        let e = find_role(&v, R_CONSENSUS);
        b[e.abs + 4] = R_INTERVALMAP;
        b[e.abs + 6] = 1; // CODEC_BSC (legal for IntervalMap)
        v5_refresh_footer_crc(b);
    });

    // Duplicate a global: turn IntervalMap into a second Consensus.
    expect_reject(d.path(), "dup_global", &bytes, |b| {
        let v = fview(b);
        let e = find_role(&v, R_INTERVALMAP);
        b[e.abs + 4] = R_CONSENSUS;
        b[e.abs + 6] = 5; // CODEC_PACKED_CONSENSUS (legal for Consensus)
        v5_refresh_footer_crc(b);
    });

    // Drop a per-chunk role: collide a Positions entry with EditCount in its group.
    expect_reject(d.path(), "drop_perchunk_role", &bytes, |b| {
        let v = fview(b);
        let e = find_role(&v, R_POSITIONS);
        b[e.abs + 4] = R_EDITCOUNT;
        v5_refresh_footer_crc(b);
    });
}

#[test]
fn single_hostile_mate2_entry_rejected() {
    // validate_reference_directory_single's "no mate-2 entries" contract: flip a
    // per-read (mate 1) entry's mate byte to 2 => single validator must reject.
    let d = tempfile::tempdir().unwrap();
    let (bytes, _path, _s1) = valid_archive(d.path());
    let v = fview(&bytes);
    // A per-read mate-1 entry (MappedFlags is per-chunk, mate 1). mate @ abs+5.
    let e = v
        .entries
        .iter()
        .find(|e| e.role == R_MAPPEDFLAGS && bytes[e.abs + 5] == 1)
        .expect("a mate-1 MappedFlags entry");
    expect_reject(d.path(), "injected_mate2", &bytes, |b| {
        b[e.abs + 5] = 2; // claim mate 2 => single validator rejects
        v5_refresh_footer_crc(b);
    });
}

#[test]
fn single_hostile_positions_m_greater_than_n() {
    let d = tempfile::tempdir().unwrap();
    let (bytes, _path, _s1) = valid_archive(d.path());
    let v = fview(&bytes);
    let pos = v
        .entries
        .iter()
        .find(|e| e.role == R_POSITIONS && bytes[e.abs..e.abs + 4] != [0xFF, 0xFF, 0xFF, 0xFF])
        .expect("per-chunk Positions entry");
    let m = read_u64(&bytes, pos.abs + 23);
    expect_reject(d.path(), "positions_m_gt_n", &bytes, |b| {
        b[pos.abs + 23..pos.abs + 31].copy_from_slice(&(m + 1).to_le_bytes());
        v5_refresh_footer_crc(b);
    });
}

#[test]
fn single_hostile_fallback_pool_count_mismatch() {
    let d = tempfile::tempdir().unwrap();
    let (bytes, _path, _s1) = valid_archive(d.path());
    let v = fview(&bytes);
    let pool = v
        .entries
        .iter()
        .find(|e| e.role == R_FALLBACKPOOL)
        .expect("a per-(chunk,mate) FallbackPool entry");
    expect_reject(d.path(), "pool_count_mismatch", &bytes, |b| {
        let rc = read_u64(b, pool.abs + 23);
        b[pool.abs + 23..pool.abs + 31].copy_from_slice(&(rc + 1).to_le_bytes());
        v5_refresh_footer_crc(b);
    });
}

#[test]
fn single_hostile_payload_crc_and_truncation() {
    let d = tempfile::tempdir().unwrap();
    let (bytes, _path, _s1) = valid_archive(d.path());
    let v = fview(&bytes);

    // Flip a byte in the EditBase payload => BSC CRC mismatch => Err.
    let eb = find_role(&v, R_EDITBASE);
    expect_reject(d.path(), "editbase_flip", &bytes, |b| {
        if eb.length > 0 {
            let last = (eb.offset + eb.length - 1) as usize;
            b[last] ^= 0xFF;
        }
    });

    // Flip a byte in the IntervalMap global payload => BSC CRC => Err.
    let im = find_role(&v, R_INTERVALMAP);
    expect_reject(d.path(), "intervalmap_flip", &bytes, |b| {
        if im.length > 0 {
            let last = (im.offset + im.length - 1) as usize;
            b[last] ^= 0xFF;
        }
    });

    // Truncate the tail into the footer => Err.
    expect_reject(d.path(), "trunc_tail", &bytes, |b| {
        let n = b.len();
        b.truncate(n - 32);
    });

    // Truncate at the footer start (one byte past the last payload) => Err.
    let foff = v.footer_offset;
    expect_reject(d.path(), "trunc_at_footer_start", &bytes, |b| {
        b.truncate(foff);
    });
}

#[test]
fn single_hostile_inflated_entry_count() {
    // Mirror of `hostile_front_and_footer`'s `inflated_entry_count` in the paired
    // suite: write u32::MAX into the footer's n_entries field (then refresh the
    // footer CRC so the generic CRC check passes), exercising the allocation-cap /
    // entry-count-overflow guard in the single-end reference directory validator.
    let d = tempfile::tempdir().unwrap();
    let (bytes, _path, _s1) = valid_archive(d.path());
    expect_reject(d.path(), "inflated_entry_count", &bytes, |b| {
        let v = fview(b);
        b[v.entry_count_abs..v.entry_count_abs + 4].copy_from_slice(&u32::MAX.to_le_bytes());
        v5_refresh_footer_crc(b);
    });
}
