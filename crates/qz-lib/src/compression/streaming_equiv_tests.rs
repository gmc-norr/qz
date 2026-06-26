//! In-crate streaming-compress golden + roundtrip matrix. The single-end compress
//! engine is now block-granular streaming only (the legacy engine and the Legacy-vs-
//! Streaming A/B seam were removed). These tests pin the public `compress()` output
//! to FNV-1a/64 goldens — the byte-identical legacy fingerprints captured before the
//! deletion — and exercise the full decode roundtrip + clean-abort contract.
use crate::cli::{CompressConfig, QualityCompressor, QualityMode, ReferenceOptions};
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

/// FNV-1a/64 over the whole archive — stable, dependency-free; any byte drift flips it.
fn archive_fingerprint(bytes: &[u8]) -> (usize, u64) {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    (bytes.len(), h)
}

/// Compress `data` through the PUBLIC `compress()` entry with `cfg_fn` overrides;
/// return the archive bytes. For ultra cases the cfg_fn sets `c.ultra = Some(level)`.
fn compress_archive(data: &[u8], cfg_fn: impl FnOnce(&mut CompressConfig)) -> Vec<u8> {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.fastq");
    fs::write(&input, data).unwrap();
    let out = dir.path().join("out.qz");
    let mut cfg = CompressConfig {
        input: vec![input],
        output: out.clone(),
        working_dir: dir.path().to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    cfg_fn(&mut cfg);
    crate::compression::compress(&cfg).unwrap();
    fs::read(&out).unwrap()
}

/// Assert the public `compress()` archive matches the captured FNV golden. A drift
/// means the on-disk format changed: if that was deliberate, regenerate the constant.
fn assert_golden(data: &[u8], want: (usize, u64), cfg_fn: impl FnOnce(&mut CompressConfig)) {
    let got = archive_fingerprint(&compress_archive(data, cfg_fn));
    assert_eq!(
        got, want,
        "archive fingerprint drift — got (len {}, 0x{:016x}), want (len {}, 0x{:016x}); \
         regenerate intentionally if a deliberate format change landed",
        got.0, got.1, want.0, want.1
    );
}

/// Paired analog of `compress_archive`: writes R1/R2 to two files and drives the
/// public `compress()` with TWO inputs (`input.len()==2` ⇒ paired routing), returning
/// the archive bytes. Mirrors `compress_archive` exactly otherwise (threads:1 for
/// determinism).
fn compress_archive_paired(r1: &[u8], r2: &[u8], cfg_fn: impl FnOnce(&mut CompressConfig)) -> Vec<u8> {
    let dir = TempDir::new().unwrap();
    let in1 = dir.path().join("r1.fastq");
    let in2 = dir.path().join("r2.fastq");
    fs::write(&in1, r1).unwrap();
    fs::write(&in2, r2).unwrap();
    let out = dir.path().join("out.qz");
    let mut cfg = CompressConfig {
        input: vec![in1, in2],
        output: out.clone(),
        working_dir: dir.path().to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    cfg_fn(&mut cfg);
    crate::compression::compress(&cfg).unwrap();
    fs::read(&out).unwrap()
}

fn assert_golden_paired(r1: &[u8], r2: &[u8], want: (usize, u64), cfg_fn: impl FnOnce(&mut CompressConfig)) {
    let got = archive_fingerprint(&compress_archive_paired(r1, r2, cfg_fn));
    assert_eq!(
        got, want,
        "paired archive fingerprint drift — got (len {}, 0x{:016x}), want (len {}, 0x{:016x}); \
         regenerate intentionally if a deliberate format change landed",
        got.0, got.1, want.0, want.1
    );
}

/// Two mates with shared read names and independent sequence/quality content. `q_kind`
/// selects the quality regime so a golden can exercise fqz (≥100K reads) vs BSC quality.
/// R2 ids equal R1 ids except the mate flag (so the R2 header-delta path is exercised).
fn gen_paired(n: usize, seq_len: usize, vary: bool) -> (Vec<u8>, Vec<u8>) {
    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    for i in 0..n {
        let l1 = if vary { seq_len + i % 7 } else { seq_len };
        let l2 = if vary { seq_len + (i + 3) % 7 } else { seq_len };
        // R1/R2 ids share the stem; differ only in the trailing mate flag (delta-friendly).
        r1.extend_from_slice(format!("@SRR{i:09}.{i} 1:N:0:ACGT\n").as_bytes());
        r1.extend((0..l1).map(|j| b"ACGT"[(i + j) % 4]));
        r1.extend_from_slice(b"\n+\n");
        r1.extend((0..l1).map(|j| b'!' + ((i + j) % 40) as u8));
        r1.push(b'\n');
        r2.extend_from_slice(format!("@SRR{i:09}.{i} 2:N:0:ACGT\n").as_bytes());
        r2.extend((0..l2).map(|j| b"ACGT"[(i + j + 1) % 4]));
        r2.extend_from_slice(b"\n+\n");
        r2.extend((0..l2).map(|j| b'!' + ((i + j + 7) % 40) as u8));
        r2.push(b'\n');
    }
    (r1, r2)
}

// Paired byte-identity goldens. Originally captured on master c305888 (compress_paired_v5)
// as the streaming-rewrite migration oracle; REGENERATED when paired archives gained the
// `ChunkDecodedSizes(n_mates=2)` global (feeding paired NUMA direct-write decode). The size
// bump over the old goldens is exactly that global (a directory entry + a small per-chunk
// table frame), so these now pin the current paired format and catch any further drift.
#[test]
fn m_paired_default_bsc_quality() {
    // 5K pairs < 100K ⇒ Auto resolves to BSC quality for both mates.
    let (r1, r2) = gen_paired(5_000, 80, false);
    assert_golden_paired(&r1, &r2, (9085, 0xe6d5d2f09da2c061), |c| { c.advanced.chunk_records = 2_000; });
}

#[test]
fn m_paired_fqz_quality() {
    // Explicit fqzcomp quality (codec 6), multi-chunk.
    let (r1, r2) = gen_paired(8_000, 60, false);
    assert_golden_paired(&r1, &r2, (59630, 0xb996fc030a0c0860), |c| {
        c.advanced.chunk_records = 2_000;
        c.advanced.quality_compressor = QualityCompressor::Fqzcomp;
    });
}

#[test]
fn m_paired_variable_length() {
    let (r1, r2) = gen_paired(6_000, 50, true);
    assert_golden_paired(&r1, &r2, (20092, 0x6299e84a3c35be4a), |c| { c.advanced.chunk_records = 1_500; });
}

// NOTE: paired REJECTS `header_compressor=Bsc` (see `paired_integration.rs::
// rejected_options_matrix`), so there is no BSC-headers paired golden — `compress()`
// would return Err and `compress_archive_paired`'s `.unwrap()` would panic.

#[test]
fn m_paired_single_chunk() {
    // One chunk (chunk_records ≥ n) — exercises the first-chunk-only path.
    let (r1, r2) = gen_paired(3_000, 100, false);
    assert_golden_paired(&r1, &r2, (3871, 0x0a8b121266738d10), |c| { c.advanced.chunk_records = 5_000; });
}

/// `gen_paired` makes R2 ids differ from R1 by ONE byte (the mate flag), so the R2-header
/// delta is tiny and the `PairedHeaders` smaller-of decision ALWAYS picks `HeaderDelta`.
/// To also pin the INDEPENDENT-COLUMNAR branch, this generator gives R2 a completely
/// different (but structured) read-name scheme, so the R1→R2 delta is literal-heavy and
/// LARGER than R2's own columnar blob → mate-2 header = `Headers`/`CODEC_COLUMNAR`.
fn gen_paired_independent_r2(n: usize, seq_len: usize) -> (Vec<u8>, Vec<u8>) {
    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    for i in 0..n {
        r1.extend_from_slice(format!("@SRR{i:09}.{i} 1:N:0:ACGT\n").as_bytes());
        r1.extend((0..seq_len).map(|j| b"ACGT"[(i + j) % 4]));
        r1.extend_from_slice(b"\n+\n");
        r1.extend((0..seq_len).map(|j| b'!' + ((i + j) % 40) as u8));
        r1.push(b'\n');
        // Unrelated, still-structured R2 names → big delta, small columnar.
        r2.extend_from_slice(format!("@LANE7:TILE{:05}:X{}:Y{}\n", n - i, (i * 31) % 9973, (i * 17) % 8101).as_bytes());
        r2.extend((0..seq_len).map(|j| b"ACGT"[(i + j + 1) % 4]));
        r2.extend_from_slice(b"\n+\n");
        r2.extend((0..seq_len).map(|j| b'!' + ((i + j + 7) % 40) as u8));
        r2.push(b'\n');
    }
    (r1, r2)
}

#[test]
fn m_paired_independent_columnar_r2_header() {
    // Forces the R2-header INDEPENDENT-COLUMNAR branch (the other goldens force delta).
    let (r1, r2) = gen_paired_independent_r2(5_000, 80);
    assert_golden_paired(&r1, &r2, (11085, 0x7720c12bdb40975c), |c| { c.advanced.chunk_records = 2_000; });
}

/// Synthetic FASTQ generator. `vary` makes per-record lengths differ (variable-length
/// framing exercises the per-record varint path); otherwise lengths are constant.
fn gen_fastq(n: usize, seq_len: usize, vary: bool) -> Vec<u8> {
    let mut s = Vec::new();
    for i in 0..n {
        let l = if vary { seq_len + (i % 7) } else { seq_len };
        s.extend_from_slice(format!("@SRR{:09}.{} read/1\n", i, i).as_bytes());
        s.extend(std::iter::repeat_n(b"ACGT"[i % 4], l));
        s.push(b'\n');
        s.extend_from_slice(b"+\n");
        s.extend(std::iter::repeat_n(b'!' + (i % 40) as u8, l));
        s.push(b'\n');
    }
    s
}

#[test]
fn m_default_fqz() {
    // ≥100K → Auto resolves quality to fqzcomp.
    let d = gen_fastq(120_000, 100, false);
    assert_golden(&d, (20108, 0x9ad924a214f6fc7e), |_| {});
}

#[test]
fn m_default_fqz_multichunk() {
    let d = gen_fastq(120_000, 100, false);
    assert_golden(&d, (21603, 0xe3bec47cc9625a89), |c| {
        c.advanced.chunk_records = 40_000;
    });
}

#[test]
fn m_bsc_quality() {
    let d = gen_fastq(5_000, 80, false);
    assert_golden(&d, (2037, 0x1dad278642ae71f3), |c| {
        c.advanced.quality_compressor = QualityCompressor::Bsc;
    });
}

#[test]
fn m_variable_length() {
    let d = gen_fastq(8_000, 60, true);
    assert_golden(&d, (9926, 0x7bd08f703498a47b), |c| {
        c.advanced.chunk_records = 2_000;
    });
}

#[test]
fn m_bsc_headers() {
    let d = gen_fastq(8_000, 80, false);
    assert_golden(&d, (8768, 0x6843b72832df61ea), |c| {
        c.advanced.header_compressor = crate::cli::HeaderCompressor::Bsc;
        c.advanced.chunk_records = 2_000;
    });
}

#[test]
fn m_no_quality() {
    let d = gen_fastq(8_000, 80, false);
    assert_golden(&d, (3698, 0x32b88cd5702c2e7c), |c| {
        c.no_quality = true;
        c.advanced.chunk_records = 2_000;
    });
}

#[test]
fn m_fasta() {
    let mut s = Vec::new();
    for i in 0..6_000 {
        s.extend_from_slice(format!(">seq{i}\n").as_bytes());
        s.extend(std::iter::repeat_n(b'A', 50));
        s.push(b'\n');
    }
    assert_golden(&s, (1182, 0x60d979b056d42c9a), |c| {
        c.fasta = true;
        c.advanced.chunk_records = 1_500;
    });
}

#[test]
fn m_const_length() {
    let d = gen_fastq(8_000, 100, false); // constant lengths
    assert_golden(&d, (6486, 0x4dbd03dbc4588ef9), |c| {
        c.advanced.chunk_records = 2_000;
    });
}

#[test]
fn m_rc_canon() {
    let d = gen_fastq(8_000, 80, false);
    assert_golden(&d, (6634, 0xf86115962a860761), |c| {
        c.advanced.rc_canon = true;
        c.advanced.sequence_hints = false; // rc_canon + hints is rejected by config guard
        c.advanced.chunk_records = 2_000;
    });
}

#[test]
fn m_sequence_hints() {
    let d = gen_fastq(8_000, 80, false);
    assert_golden(&d, (6366, 0xd2ca964d360f23af), |c| {
        c.advanced.sequence_hints = true;
        c.advanced.chunk_records = 2_000;
    });
}

#[test]
fn m_quality_binning() {
    let d = gen_fastq(8_000, 80, false);
    assert_golden(&d, (1599, 0x6dec5b237b0c6765), |c| {
        c.quality_mode = QualityMode::IlluminaBin;
    });
}

/// Original small multi-chunk case (BSC quality, columnar headers) — kept as a
/// fast smoke test of the engine.
#[test]
fn m_default_multichunk_small() {
    // Distinct ids, constant 16bp seq AND 16-char qual.
    const SMALL: &str = "\
@SRR1.1 a\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n\
@SRR1.2 b\nTGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHH\n\
@SRR1.3 c\nAAAACCCCGGGGTTTT\n+\nBBBBBBBBBBBBBBBB\n\
@SRR1.4 d\nGGGGTTTTAAAACCCC\n+\n????????????????\n";
    assert_golden(SMALL.as_bytes(), (738, 0x43d94296e148281c), |c| {
        c.advanced.chunk_records = 2;
    });
}

/// A Streaming archive must decode (roundtrip) to the original FASTQ content.
/// Compresses through the public `compress()`, then `decompress`es and asserts equality.
#[test]
fn streaming_roundtrip_decodes() {
    // A multi-chunk, fqz-eligible input plus a small BSC-quality input, to cover both
    // quality backends through the full decode path.
    for (label, data, threshold_small) in [
        ("fqz-multichunk", gen_fastq(120_000, 100, false), false),
        ("bsc-small", gen_fastq(3_000, 60, true), true),
    ] {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("in.fastq");
        fs::write(&input, &data).unwrap();
        let arc = dir.path().join("out.qz");
        let mut cfg = CompressConfig {
            input: vec![input],
            output: arc.clone(),
            working_dir: dir.path().to_path_buf(),
            threads: 1,
            ..CompressConfig::default()
        };
        if threshold_small {
            cfg.advanced.chunk_records = 1_000;
        } else {
            cfg.advanced.chunk_records = 40_000;
        }
        crate::compression::compress(&cfg).unwrap();

        let out_fq = dir.path().join("decoded.fastq");
        let dcfg = crate::cli::DecompressConfig {
            input: arc,
            output: vec![out_fq.clone()],
            working_dir: dir.path().to_path_buf(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 0,
            force: true,
        };
        crate::compression::decompress(&dcfg).unwrap();
        let decoded = fs::read(&out_fq).unwrap();
        assert_eq!(
            decoded, data,
            "streaming archive must decode byte-identically for case '{label}'"
        );
    }
}

/// Cleanup contract: a producer error mid-stream (here a const-length mismatch in
/// chunk 1, AFTER the front header is written to the temp) must return Err, leave NO
/// archive at the final path, and leave NO sibling `.tmp` behind (the drop guard
/// fired). Cheap deterministic trigger — no 128 MB oversize record needed.
#[test]
fn streaming_aborts_clean_on_producer_error() {
    // chunk 0 (2 records, len 10) sets const_seq_len=10; r3 (len 11) trips
    // validate_const_one in the producer.
    let data = "\
@r1\nACGTACGTAC\n+\nIIIIIIIIII\n\
@r2\nACGTACGTAC\n+\nIIIIIIIIII\n\
@r3\nACGTACGTACG\n+\nIIIIIIIIIII\n";
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.fastq");
    fs::write(&input, data).unwrap();
    let out = dir.path().join("out.qz");
    let mut cfg = CompressConfig {
        input: vec![input],
        output: out.clone(),
        working_dir: dir.path().to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 2;

    let r = crate::compression::compress(&cfg);
    assert!(r.is_err(), "const-length mismatch must error");
    assert!(!out.exists(), "no partial archive may remain at the final path");
    let leftover: Vec<_> = fs::read_dir(dir.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().contains(".tmp"))
        .collect();
    assert!(
        leftover.is_empty(),
        "the sibling temp must be cleaned up, found {leftover:?}"
    );
}

/// Ultra L1 (big blocks) golden. 120K reads ⇒ fqz quality, so this also exercises
/// fqz-through-ultra (not just BSC big blocks).
#[test]
fn m_ultra_identical() {
    let d = gen_fastq(120_000, 100, false);
    assert_golden(&d, (20108, 0x263a96ed1b5a6428), |c| {
        c.ultra = Some(1);
    });
}

/// Degenerate tiny input through ultra (the MT-BSC degenerate-corruption guard,
/// finding_ultra_mtbsc_degenerate_corruption): compresses via the public ultra path
/// AND decodes losslessly.
#[test]
fn m_ultra_degenerate_tiny_identical() {
    let d = gen_fastq(50, 40, true); // tiny, variable-length
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.fastq");
    fs::write(&input, &d).unwrap();
    let arc = dir.path().join("ud.qz");
    let mut cfg = CompressConfig {
        input: vec![input],
        output: arc.clone(),
        working_dir: dir.path().to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    cfg.ultra = Some(1);
    crate::compression::compress(&cfg).unwrap();

    // The ultra archive must decode losslessly.
    let out_fq = dir.path().join("ud.fastq");
    let dcfg = crate::cli::DecompressConfig {
        input: arc,
        output: vec![out_fq.clone()],
        working_dir: dir.path().to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 0,
        force: true,
    };
    crate::compression::decompress(&dcfg).unwrap();
    assert_eq!(
        fs::read(&out_fq).unwrap(),
        d,
        "tiny ultra archive must decode losslessly"
    );
}

// ===========================================================================
// Reference byte-identity goldens (Increment 3 gate).
//
// These pin `compress_reference`'s archive bytes on the current (unmodified)
// master code, so the upcoming reference-compress convergence onto shared
// concurrency machinery is gated on byte-identity. Mirrors the single-end `m_*`
// and paired `m_paired_*` goldens above. Deterministic — index arithmetic only,
// no `rand` — so the goldens are stable across runs. The generators below produce
// reads that strobealign actually maps (sliced from the reference, R2 reverse-
// complemented), with optional unmapped (literal-fallback) pairs and edits.
// ===========================================================================

/// Deterministic LCG base at `(contig, pos)`. Mirrors the proven `make_seq`
/// generator the reference integration tests use (low-repeat ACGT that
/// strobealign seeds reliably), keyed so each contig is a distinct sequence.
fn ref_base(contig: usize, pos: usize) -> u8 {
    let mut x = (contig as u64)
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(pos as u64)
        .wrapping_add(0x1234_5678_9ABC_DEF0);
    // One LCG step decorrelates adjacent positions (a single multiply alone
    // leaves low-bit structure that maps poorly).
    x = x
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    b"ACGT"[((x >> 33) & 3) as usize]
}

/// Reverse-complement an ACGT slice (non-ACGT passes through).
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            o => o,
        })
        .collect()
}

/// Write a FASTA of `n_contigs` contigs × `len` bp (bases from `ref_base`).
/// Returns its path. The whole reference is also returned via `read_reference`.
fn gen_reference(dir: &std::path::Path, n_contigs: usize, len: usize) -> PathBuf {
    let mut s = String::new();
    for c in 0..n_contigs {
        s.push_str(&format!(">contig{c}\n"));
        let line: Vec<u8> = (0..len).map(|p| ref_base(c, p)).collect();
        s.push_str(std::str::from_utf8(&line).unwrap());
        s.push('\n');
    }
    let p = dir.join("ref.fa");
    fs::write(&p, &s).unwrap();
    p
}

/// Re-read a `gen_reference` FASTA back into per-contig byte vectors (so reads are
/// sliced from EXACTLY the stored bytes).
fn read_reference(reference: &std::path::Path) -> Vec<Vec<u8>> {
    let txt = fs::read_to_string(reference).unwrap();
    let mut contigs: Vec<Vec<u8>> = Vec::new();
    for line in txt.lines() {
        if line.starts_with('>') {
            contigs.push(Vec::new());
        } else if let Some(last) = contigs.last_mut() {
            last.extend_from_slice(line.as_bytes());
        }
    }
    contigs
}

/// Generate a paired FASTQ pair against `reference`. For pair `i`: pick a contig +
/// start by `i`-arithmetic, slice `read_len` bp for R1 (forward); slice a downstream
/// window for R2 and store it reverse-complemented (so the placer's RC path fires).
/// `unmapped_every>0` makes every k-th pair fixed non-mapping literals (NOT sampled
/// from the reference → strobealign won't map → FallbackPool). `edit_every>0` injects
/// a substitution every k bases of the mapped reads (→ EditCount/EditPos/EditBase).
/// A fixed quality string is used. Returns `(r1.fastq, r2.fastq)`.
fn gen_ref_reads(
    dir: &std::path::Path,
    reference: &std::path::Path,
    n_pairs: usize,
    read_len: usize,
    unmapped_every: usize,
    edit_every: usize,
) -> (PathBuf, PathBuf) {
    let contigs = read_reference(reference);
    let n_contigs = contigs.len();
    let q: String = "I".repeat(read_len);
    let mut s1 = String::new();
    let mut s2 = String::new();
    // R2 sits `insert` bp downstream of R1; both must fit inside the contig.
    let insert = read_len / 2;
    for i in 0..n_pairs {
        let unmapped = unmapped_every > 0 && i % unmapped_every == 0;
        let (r1bytes, r2bytes): (Vec<u8>, Vec<u8>) = if unmapped {
            // Deterministic non-reference literals: a different keyspace
            // (huge contig index) than any real contig, so strobealign can't
            // seed these against the reference → literal fallback.
            let r1: Vec<u8> = (0..read_len)
                .map(|p| ref_base(1_000_000 + i, p))
                .collect();
            let r2: Vec<u8> = (0..read_len)
                .map(|p| ref_base(2_000_000 + i, p))
                .collect();
            (r1, r2)
        } else {
            let c = i % n_contigs;
            let contig = &contigs[c];
            let span = read_len + insert; // R1 + downstream offset + R2 footprint
            // Invariant for any future case: the contig must hold a full R1+insert+R2
            // window. Turns a cryptic subtraction underflow into an actionable message.
            assert!(
                contig.len() >= span,
                "contig {c} too short ({} bp) for read_len {read_len} + insert {insert}",
                contig.len()
            );
            let max_start = contig.len() - span;
            // Spread starts deterministically across the mappable range.
            let st1 = (i.wrapping_mul(37)) % (max_start + 1);
            let mut r1 = contig[st1..st1 + read_len].to_vec();
            let st2 = st1 + insert;
            let r2_fwd = contig[st2..st2 + read_len].to_vec();
            let mut r2 = revcomp(&r2_fwd);
            if edit_every > 0 {
                // Inject a substitution every `edit_every` bases of each mate
                // (deterministic positions; flip to a guaranteed-different base).
                let mut p = edit_every.min(read_len.saturating_sub(1));
                while p < read_len {
                    r1[p] = flip_base(r1[p]);
                    p += edit_every;
                }
                let mut p = edit_every.min(read_len.saturating_sub(1));
                while p < read_len {
                    r2[p] = flip_base(r2[p]);
                    p += edit_every;
                }
            }
            (r1, r2)
        };
        s1.push_str(&format!(
            "@read_{i}/1\n{}\n+\n{q}\n",
            std::str::from_utf8(&r1bytes).unwrap()
        ));
        s2.push_str(&format!(
            "@read_{i}/2\n{}\n+\n{q}\n",
            std::str::from_utf8(&r2bytes).unwrap()
        ));
    }
    let r1p = dir.join("r1.fastq");
    let r2p = dir.join("r2.fastq");
    fs::write(&r1p, &s1).unwrap();
    fs::write(&r2p, &s2).unwrap();
    (r1p, r2p)
}

/// Flip an ACGT base to a guaranteed-different ACGT base (deterministic).
fn flip_base(b: u8) -> u8 {
    match b {
        b'A' => b'C',
        b'C' => b'A',
        b'G' => b'T',
        _ => b'G',
    }
}

/// Reference analog of `compress_archive`: drives the public `compress()` with the
/// reference block set (two inputs R1/R2 + a `ReferenceOptions`), returning the
/// archive bytes. `threads:1`, `reference_window:2` for determinism (output is
/// byte-identical across both anyway, but pin them so the golden is reproducible).
fn compress_archive_reference(
    reference: &std::path::Path,
    r1: &std::path::Path,
    r2: &std::path::Path,
    chunk_records: usize,
) -> Vec<u8> {
    let dir = TempDir::new().unwrap();
    let out = dir.path().join("out.qz");
    let mut cfg = CompressConfig {
        input: vec![r1.to_path_buf(), r2.to_path_buf()],
        output: out.clone(),
        working_dir: dir.path().to_path_buf(),
        threads: 1,
        force: true,
        reference: Some(ReferenceOptions {
            reference: reference.to_path_buf(),
            reference_index: None,
            reference_fast: false,
            reference_window: 2,
        }),
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = chunk_records;
    crate::compression::compress(&cfg).unwrap();
    fs::read(&out).unwrap()
}

fn assert_golden_reference(name: &str, bytes: &[u8], expected_fnv: u64) {
    let got = archive_fingerprint(bytes);
    assert_eq!(
        got.1, expected_fnv,
        "reference golden '{name}' drift — got (len {}, 0x{:016x}), want 0x{:016x}; \
         regenerate intentionally if a deliberate format change landed",
        got.0, got.1, expected_fnv
    );
}

#[test]
fn m_ref_default_single_chunk() {
    // ~3000 pairs, 100bp, 2 contigs × 20_000bp, one chunk (chunk_records ≥ n).
    // No unmapped, no edits → all reads map, no FallbackPool, no edits.
    let dir = TempDir::new().unwrap();
    let reference = gen_reference(dir.path(), 2, 20_000);
    let (r1, r2) = gen_ref_reads(dir.path(), &reference, 3_000, 100, 0, 0);
    let bytes = compress_archive_reference(&reference, &r1, &r2, 1_000_000);
    assert_golden_reference("m_ref_default_single_chunk", &bytes, 0xe9b635df275776cd);
}

#[test]
fn m_ref_default_multichunk() {
    // Same shape but chunk_records=1_000 (≈3 chunks).
    let dir = TempDir::new().unwrap();
    let reference = gen_reference(dir.path(), 2, 20_000);
    let (r1, r2) = gen_ref_reads(dir.path(), &reference, 3_000, 100, 0, 0);
    let bytes = compress_archive_reference(&reference, &r1, &r2, 1_000);
    assert_golden_reference("m_ref_default_multichunk", &bytes, 0xa5f20ef635e38db0);
}

#[test]
fn m_ref_with_unmapped_fallback_pool() {
    // unmapped_every=7 ⇒ FallbackPool present (F>0), multi-chunk.
    let dir = TempDir::new().unwrap();
    let reference = gen_reference(dir.path(), 2, 20_000);
    let (r1, r2) = gen_ref_reads(dir.path(), &reference, 3_000, 100, 7, 0);
    let bytes = compress_archive_reference(&reference, &r1, &r2, 1_000);
    assert_golden_reference("m_ref_with_unmapped_fallback_pool", &bytes, 0x8c16f530f4074d89);
}

#[test]
fn m_ref_all_mapped_no_fallback() {
    // unmapped_every=0 + a longer reference (40k vs 20k) so mappability is high enough that
    // EVERY read seeds and maps → FallbackPool absent (F==0). Shrinking it can re-introduce
    // fallbacks and silently flip this branch.
    let dir = TempDir::new().unwrap();
    let reference = gen_reference(dir.path(), 2, 40_000);
    let (r1, r2) = gen_ref_reads(dir.path(), &reference, 3_000, 100, 0, 0);
    let bytes = compress_archive_reference(&reference, &r1, &r2, 1_000);
    assert_golden_reference("m_ref_all_mapped_no_fallback", &bytes, 0x0ceb8f7e378f73dd);
}

#[test]
fn m_ref_with_edits() {
    // edit_every=25 (substitutions) ⇒ non-empty EditPos. Keep edits-per-read below the
    // placer's exact-anchor k so edited reads still MAP (→ EditPos), not fall to FallbackPool.
    let dir = TempDir::new().unwrap();
    let reference = gen_reference(dir.path(), 2, 20_000);
    let (r1, r2) = gen_ref_reads(dir.path(), &reference, 3_000, 100, 0, 25);
    let bytes = compress_archive_reference(&reference, &r1, &r2, 1_000);
    assert_golden_reference("m_ref_with_edits", &bytes, 0xad091724280ec60b);
}
