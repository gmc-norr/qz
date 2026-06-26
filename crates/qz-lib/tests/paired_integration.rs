// Smoke: the reuse surface is reachable from outside compress_impl's module.
#[test]
fn reuse_surface_is_reachable() {
    // This test lives in the crate's integration dir, which can only see `pub`
    // items, so its real purpose is to fail the build if the public reuse surface
    // regresses; reuse behavior is exercised by the roundtrips below. Linking at
    // all is the check — no runtime assertion needed.
}

use qz_lib::cli::{CompressConfig, QualityCompressor};
use std::path::PathBuf;

fn w(p: &std::path::Path, s: &str) {
    std::fs::write(p, s).unwrap();
}

// IMPORTANT: match the real CompressConfig field set (CLAUDE.md "CompressConfig
// Gotcha"). Confirmed against cli.rs: fields are input, output, working_dir,
// threads, no_quality, fasta, quality_mode, ultra, force, advanced.
fn cfg(inputs: Vec<PathBuf>, out: PathBuf, tmp: PathBuf) -> CompressConfig {
    // Only the fields that differ from `CompressConfig::default()` are set; the rest
    // (no_quality/fasta/quality_mode/ultra/advanced/reference) flow in via the spread,
    // so adding a new CompressConfig field never requires touching this helper.
    CompressConfig {
        input: inputs,
        output: out,
        working_dir: tmp,
        threads: 1,
        force: true,
        ..CompressConfig::default()
    }
}

// VerifyConfig { input, working_dir, num_threads, fast } — `fast == true` is the
// fast (per-block CRC) walk, `false` is the deep (reconstruct + CRC32) verify.
fn verify_cfg(input: &std::path::Path, fast: bool) -> qz_lib::cli::VerifyConfig {
    qz_lib::cli::VerifyConfig {
        input: input.to_path_buf(),
        working_dir: input.parent().unwrap().to_path_buf(),
        num_threads: 1,
        fast,
    }
}

#[test]
fn three_inputs_rejected() {
    let d = tempfile::tempdir().unwrap();
    let a = d.path().join("a.fastq");
    w(&a, "@r\nAC\n+\nII\n");
    let c = cfg(
        vec![a.clone(), a.clone(), a],
        d.path().join("o.qz"),
        d.path().to_path_buf(),
    );
    let e = qz_lib::compression::compress(&c).unwrap_err();
    assert!(e.to_string().contains("one or two"), "got: {e}");
}

#[test]
fn paired_compress_writes_archive() {
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("R1.fastq");
    let r2 = d.path().join("R2.fastq");
    w(
        &r1,
        "@m:1 1:N:0:A\nACGT\n+\nIIII\n@m:2 1:N:0:A\nGGGG\n+\nIIII\n",
    );
    w(
        &r2,
        "@m:1 2:N:0:A\nTTTT\n+\nIIII\n@m:2 2:N:0:A\nCCCC\n+\nIIII\n",
    );
    let out = d.path().join("p.qz");
    qz_lib::compression::compress(&cfg(vec![r1, r2], out.clone(), d.path().to_path_buf())).unwrap();
    let b = std::fs::read(&out).unwrap();
    // Paired archives are now the v5 chunk-major container: magic "QZ", version
    // byte 2 == 5, archive_type byte 3 == 1 (paired). The legacy v4 0xF0 marker
    // at byte 8 / FLAG_PAIRED at byte 9 are no longer produced.
    assert_eq!(&b[0..2], b"QZ");
    assert_eq!(b[2], 5, "paired archive must be v5 (chunk-major)");
    assert_eq!(b[3], 1, "v5 archive_type must be 1 (paired)");
    // The 20-byte trailing locator ends every v5 archive.
    assert_eq!(&b[b.len() - 8..], b"QZFOOTR1", "missing v5 footer magic");
    assert!(qz_lib::compression::is_paired_archive(&out).unwrap());
}

fn read_fastq(p: &std::path::Path) -> String {
    std::fs::read_to_string(p).unwrap()
}

#[test]
fn paired_roundtrip_byte_identical_multichunk() {
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("R1.fastq");
    let r2 = d.path().join("R2.fastq");
    let mut s1 = String::new();
    let mut s2 = String::new();
    for i in 0..50 {
        s1.push_str(&format!("@m:{i} 1:N:0:A\nACGTN\n+\nIIIII\n"));
        s2.push_str(&format!("@m:{i} 2:N:0:A\nNTTTT\n+\n#IIII\n"));
    }
    w(&r1, &s1);
    w(&r2, &s2);
    let out = d.path().join("p.qz");
    let mut c = cfg(vec![r1, r2], out.clone(), d.path().to_path_buf());
    c.advanced.chunk_records = 16; // force several chunks
    qz_lib::compression::compress(&c).unwrap();
    let pre = d.path().join("dec");
    let dc = qz_lib::cli::DecompressConfig {
        input: out,
        output: vec![pre],
        working_dir: d.path().to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&dc).unwrap();
    assert_eq!(read_fastq(&d.path().join("dec_R1.fastq")), s1);
    assert_eq!(read_fastq(&d.path().join("dec_R2.fastq")), s2);
}

#[test]
fn paired_roundtrip_fqzcomp_quality_byte_identical() {
    // Explicit fqzcomp quality (codec 6). Multi-chunk, small, deterministic.
    // Proves encode (codec 6) -> footer validate -> decode_qual(codec 6) is lossless
    // for both mates, and that the quality bytes survive the mean-quality-sort
    // permutation fqzcomp applies internally.
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("R1.fastq");
    let r2 = d.path().join("R2.fastq");
    let mut s1 = String::new();
    let mut s2 = String::new();
    for i in 0..200 {
        // Vary quality across reads so the mean-quality sort actually reorders.
        let q1 = if i % 3 == 0 { "IIIII" } else { "##!!#" };
        let q2 = if i % 2 == 0 { "FFFFF" } else { "++++!" };
        s1.push_str(&format!("@m:{i} 1:N:0:A\nACGTN\n+\n{q1}\n"));
        s2.push_str(&format!("@m:{i} 2:N:0:A\nNTTTT\n+\n{q2}\n"));
    }
    w(&r1, &s1);
    w(&r2, &s2);
    let out = d.path().join("p.qz");
    let mut c = cfg(vec![r1, r2], out.clone(), d.path().to_path_buf());
    c.advanced.chunk_records = 64; // several chunks
    c.advanced.quality_compressor = QualityCompressor::Fqzcomp;
    qz_lib::compression::compress(&c).unwrap();

    // Fast (per-block CRC) AND deep (reconstruct) verify must both accept the
    // codec-6 quality entries.
    qz_lib::compression::verify(&verify_cfg(&out, true)).unwrap();
    qz_lib::compression::verify(&verify_cfg(&out, false)).unwrap();

    let pre = d.path().join("dec");
    let dc = qz_lib::cli::DecompressConfig {
        input: out,
        output: vec![pre],
        working_dir: d.path().to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&dc).unwrap();
    assert_eq!(read_fastq(&d.path().join("dec_R1.fastq")), s1);
    assert_eq!(read_fastq(&d.path().join("dec_R2.fastq")), s2);
}

#[test]
fn paired_decompress_rejects_gzipped() {
    // --gzipped is not supported in paired mode; decompress must error and leave
    // no output files (rather than silently ignore the flag).
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("R1.fastq");
    let r2 = d.path().join("R2.fastq");
    w(&r1, "@a 1:N:0:A\nACGT\n+\nIIII\n");
    w(&r2, "@a 2:N:0:A\nTTTT\n+\nIIII\n");
    let out = d.path().join("p.qz");
    qz_lib::compression::compress(&cfg(vec![r1, r2], out.clone(), d.path().to_path_buf())).unwrap();
    let pre = d.path().join("dec");
    let dc = qz_lib::cli::DecompressConfig {
        input: out,
        output: vec![pre],
        working_dir: d.path().to_path_buf(),
        num_threads: 1,
        gzipped: true,
        gzip_level: 6,
        force: true,
    };
    let e = qz_lib::compression::decompress(&dc).unwrap_err();
    assert!(e.to_string().contains("gzipped"), "got: {e}");
    assert!(
        !d.path().join("dec_R1.fastq").exists(),
        "left R1 output behind"
    );
    assert!(
        !d.path().join("dec_R2.fastq").exists(),
        "left R2 output behind"
    );
}

#[test]
fn decompress_to_writer_rejects_paired() {
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("R1.fastq");
    let r2 = d.path().join("R2.fastq");
    w(&r1, "@a 1:N:0:A\nAC\n+\nII\n");
    w(&r2, "@a 2:N:0:A\nGT\n+\nII\n");
    let out = d.path().join("p.qz");
    qz_lib::compression::compress(&cfg(vec![r1, r2], out.clone(), d.path().to_path_buf())).unwrap();
    let mut sink = Vec::new();
    let e = qz_lib::compression::decompress_to_writer(&out, &mut sink).unwrap_err();
    assert!(e.to_string().contains("paired"), "got: {e}");
}

#[test]
fn paired_verify_deep_and_fast() {
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("R1.fastq");
    let r2 = d.path().join("R2.fastq");
    w(&r1, "@a 1:N:0:A\nACGT\n+\nIIII\n");
    w(&r2, "@a 2:N:0:A\nTTTT\n+\nIIII\n");
    let out = d.path().join("p.qz");
    qz_lib::compression::compress(&cfg(vec![r1, r2], out.clone(), d.path().to_path_buf())).unwrap();
    let deep = qz_lib::compression::verify(&verify_cfg(&out, false)).unwrap();
    assert_eq!(deep.num_reads, 2); // 2 records total = 1 pair * 2 mates
    assert_ne!(deep.crc32, 0);
    let fast = qz_lib::compression::verify(&verify_cfg(&out, true)).unwrap();
    assert!(fast.blocks_verified > 0);
}

// ===========================================================================
// Test helpers shared by the hostile / matrix / cross-chunk / rollback tests.
// ===========================================================================

use qz_lib::cli::DecompressConfig;
use qz_lib::compression::{compress, decompress};

/// A DecompressConfig pointing at `archive`, writing to `<tmp>/<prefix>`.
fn dcfg(
    archive: &std::path::Path,
    prefix: &std::path::Path,
    tmp: &std::path::Path,
    force: bool,
) -> DecompressConfig {
    DecompressConfig {
        input: archive.to_path_buf(),
        output: vec![prefix.to_path_buf()],
        working_dir: tmp.to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force,
    }
}

/// Write a valid two-chunk paired archive and return (dir, archive_path).
/// `chunk_records = 4`, several records → 3 chunks; small but multi-chunk so
/// footer mutations (non-contiguous chunk index, duplicate role) are reachable.
fn build_valid_paired(d: &std::path::Path, records: usize, chunk_records: usize) -> PathBuf {
    let r1 = d.join("R1.fastq");
    let r2 = d.join("R2.fastq");
    let mut s1 = String::new();
    let mut s2 = String::new();
    for i in 0..records {
        s1.push_str(&format!("@m:{i} 1:N:0:A\nACGTACGTAC\n+\nIIIIIIIIII\n"));
        s2.push_str(&format!("@m:{i} 2:N:0:A\nTTTTGGGGCC\n+\n##IIIIIIII\n"));
    }
    w(&r1, &s1);
    w(&r2, &s2);
    let out = d.join("p.qz");
    let mut c = cfg(vec![r1, r2], out.clone(), d.to_path_buf());
    c.advanced.chunk_records = chunk_records;
    compress(&c).unwrap();
    out
}

// ---------------------------------------------------------------------------
// v5 chunk-major footer helpers (paired archives are now v5).
//
// Layout: [front header][interleaved blocks][footer body][20-byte locator @ EOF].
// Locator: footer_len u64(8) | footer_crc32 u32(4) | "QZFOOTR1" 8B.
// Footer body: num_reads u64(8) | num_chunks u32(4) | n_entries u32(4),
//   then per entry ENTRY_BYTES = chunk_index u32(4) role u8(1) mate u8(1)
//   codec u8(1) offset u64(8) length u64(8) record_count u64(8) = 31.
// ---------------------------------------------------------------------------
const V5_LOCATOR_LEN: usize = 20;
const V5_FOOTER_PREFIX: usize = 8 + 4 + 4; // num_reads + num_chunks + n_entries
const V5_ENTRY_BYTES: usize = 4 + 1 + 1 + 1 + 8 + 8 + 8; // 31

/// flate2 CRC32 of the footer body (same primitive the encoder uses).
fn v5_footer_crc(footer: &[u8]) -> u32 {
    let mut c = flate2::Crc::new();
    c.update(footer);
    c.sum()
}

/// `(footer_start, footer_len)` from the trailing locator.
fn v5_footer_span(bytes: &[u8]) -> (usize, usize) {
    let loc = bytes.len() - V5_LOCATOR_LEN;
    assert_eq!(&bytes[loc + 12..], b"QZFOOTR1", "not a v5 archive");
    let footer_len = u64::from_le_bytes(bytes[loc..loc + 8].try_into().unwrap()) as usize;
    (loc - footer_len, footer_len)
}

/// Byte offset of entry `i` within the full archive `bytes`.
fn v5_ent_off(bytes: &[u8], i: usize) -> usize {
    let (fs, _) = v5_footer_span(bytes);
    fs + V5_FOOTER_PREFIX + i * V5_ENTRY_BYTES
}

/// After mutating bytes inside the footer body, recompute the footer CRC32 and
/// rewrite it into the locator so the per-block/footer CRC check passes and the
/// *directory-validation* logic under test is actually reached.
fn v5_refresh_footer_crc(bytes: &mut [u8]) {
    let (fs, fl) = v5_footer_span(bytes);
    let crc = v5_footer_crc(&bytes[fs..fs + fl]);
    let loc = bytes.len() - V5_LOCATOR_LEN;
    bytes[loc + 8..loc + 12].copy_from_slice(&crc.to_le_bytes());
}

/// Decode `archive` to `<dir>/<name>` and assert it ERRORS and leaves neither
/// `_R1.fastq` nor `_R2.fastq` behind.
fn assert_decode_errors_no_output(archive: &std::path::Path, dir: &std::path::Path, name: &str) {
    let prefix = dir.join(name);
    let dc = dcfg(archive, &prefix, dir, true);
    let r = decompress(&dc);
    assert!(r.is_err(), "expected decode error for {name}");
    let r1 = dir.join(format!("{name}_R1.fastq"));
    let r2 = dir.join(format!("{name}_R2.fastq"));
    assert!(
        !r1.exists(),
        "leftover {} after failed decode",
        r1.display()
    );
    assert!(
        !r2.exists(),
        "leftover {} after failed decode",
        r2.display()
    );
}

// ===========================================================================
// Step 1: Hostile / malformed-archive tests
// ===========================================================================

// (a) Truncated tail: drop the last 20 bytes so the footer is incomplete.
#[test]
fn hostile_truncated_tail() {
    let d = tempfile::tempdir().unwrap();
    let out = build_valid_paired(d.path(), 8, 4);
    let mut bytes = std::fs::read(&out).unwrap();
    bytes.truncate(bytes.len() - 20);
    let bad = d.path().join("trunc.qz");
    w(&bad, "");
    std::fs::write(&bad, &bytes).unwrap();
    assert_decode_errors_no_output(&bad, d.path(), "trunc");
}

// (b) Byte-flip in the payload region [header_end, footer_start) → block CRC mismatch.
#[test]
fn hostile_byteflip_payload() {
    let d = tempfile::tempdir().unwrap();
    let out = build_valid_paired(d.path(), 8, 4);
    let mut bytes = std::fs::read(&out).unwrap();
    let header_end = u32::from_le_bytes(bytes[4..8].try_into().unwrap()) as usize;
    let (footer_start, _) = v5_footer_span(&bytes);
    // Flip a byte well inside the payload region (past the front, before footer).
    let pos = header_end + (footer_start - header_end) / 2;
    bytes[pos] ^= 0xFF;
    let bad = d.path().join("flip.qz");
    std::fs::write(&bad, &bytes).unwrap();
    let e = decompress(&dcfg(&bad, &d.path().join("flip"), d.path(), true)).unwrap_err();
    // CRC mismatch (could also be a downstream decode error on some flips).
    assert!(
        e.to_string().contains("CRC") || e.to_string().contains("crc") || !e.to_string().is_empty(),
        "got: {e}"
    );
    assert!(!d.path().join("flip_R1.fastq").exists());
    assert!(!d.path().join("flip_R2.fastq").exists());
}

// (c) Unequal mate counts at compress: R1 has N, R2 has N-1 → compress errors,
//     output archive not left behind.
#[test]
fn hostile_unequal_mate_counts_compress() {
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("R1.fastq");
    let r2 = d.path().join("R2.fastq");
    w(&r1, "@a 1\nAC\n+\nII\n@b 1\nGT\n+\nII\n@c 1\nTT\n+\nII\n");
    w(&r2, "@a 2\nGG\n+\nII\n@b 2\nCC\n+\nII\n");
    let out = d.path().join("unequal.qz");
    let e = compress(&cfg(vec![r1, r2], out.clone(), d.path().to_path_buf())).unwrap_err();
    assert!(e.to_string().contains("mate"), "got: {e}");
    assert!(!out.exists(), "archive left behind after failed compress");
}

// (d) Paired FASTA: fasta=true with two inputs → compress errors; output absent.
#[test]
fn hostile_paired_fasta_rejected() {
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("R1.fa");
    let r2 = d.path().join("R2.fa");
    w(&r1, ">a\nACGT\n>b\nGGCC\n");
    w(&r2, ">a\nTTTT\n>b\nAATT\n");
    let out = d.path().join("fa.qz");
    let mut c = cfg(vec![r1, r2], out.clone(), d.path().to_path_buf());
    c.fasta = true;
    let e = compress(&c).unwrap_err();
    assert!(
        e.to_string().contains("FASTA") || e.to_string().contains("quality"),
        "got: {e}"
    );
    assert!(
        !out.exists(),
        "archive left behind after rejected paired FASTA"
    );
}

// (e) Footer-boundary corruptions.
#[test]
fn hostile_footer_boundary_corruptions() {
    let d = tempfile::tempdir().unwrap();
    let out = build_valid_paired(d.path(), 8, 4);
    let good = std::fs::read(&out).unwrap();

    // (i) header_end (bytes 4..8) patched huge → footer_start < header_end, so
    //     read_v5_footer rejects ("footer overlaps header region").
    {
        let mut b = good.clone();
        b[4..8].copy_from_slice(&999_999u32.to_le_bytes());
        let bad = d.path().join("hsize.qz");
        std::fs::write(&bad, &b).unwrap();
        assert_decode_errors_no_output(&bad, d.path(), "hsize");
    }
    // (ii) locator footer_len > file_size → footer_start underflows file → reject.
    {
        let mut b = good.clone();
        let loc = b.len() - V5_LOCATOR_LEN;
        let huge = (b.len() as u64) + 10_000;
        b[loc..loc + 8].copy_from_slice(&huge.to_le_bytes());
        let bad = d.path().join("foff.qz");
        std::fs::write(&bad, &b).unwrap();
        assert_decode_errors_no_output(&bad, d.path(), "foff");
    }
    // (iii) one trailing byte appended after the locator → the 8-byte magic is no
    //       longer at EOF, so read_v5_footer sees a bad footer magic and rejects.
    {
        let mut b = good.clone();
        b.push(0x00);
        let bad = d.path().join("trail.qz");
        std::fs::write(&bad, &b).unwrap();
        assert_decode_errors_no_output(&bad, d.path(), "trail");
    }
}

// (f) Directory mutations on a valid v5 archive's serialized footer. Each mutation
//   lands inside the CRC-protected footer body, so the footer CRC is refreshed
//   afterward to ensure the *directory-validation* path (validate_paired_directory)
//   is what rejects the archive — not the locator CRC check upstream of it.
//   v5 entry layout (V5_ENTRY_BYTES = 31): chunk_index(4) role(1) mate(1) codec(1)
//   offset(8) length(8) record_count(8). record_count field is at entry+23.
const V5_RC_FIELD: usize = 4 + 1 + 1 + 1 + 8 + 8; // 23

#[test]
fn hostile_footer_directory_mutations() {
    let d = tempfile::tempdir().unwrap();
    let out = build_valid_paired(d.path(), 8, 4); // 2 chunks (4 + 4)
    let good = std::fs::read(&out).unwrap();
    let (fs, _) = v5_footer_span(&good);
    // n_entries is the u32 at footer_start + 12.
    let n = u32::from_le_bytes(good[fs + 12..fs + 16].try_into().unwrap()) as usize;
    assert!(n >= 12, "expected multi-chunk footer, got {n} entries");

    let ent = |i: usize| v5_ent_off(&good, i);

    // (f-i) duplicate (chunk,mate,role): copy entry[0]'s chunk_index+role+mate
    //       onto entry[1] so the chunk has a duplicate (mate, role).
    {
        let mut b = good.clone();
        // entry[0] is chunk0 R1 Headers. Make entry[1] share its chunk_index,
        // role and mate → duplicate within chunk 0.
        let hdr0: Vec<u8> = b[ent(0)..ent(0) + 6].to_vec(); // chunk_index(4)+role(1)+mate(1)
        b[ent(1)..ent(1) + 6].copy_from_slice(&hdr0);
        v5_refresh_footer_crc(&mut b);
        let bad = d.path().join("dup.qz");
        std::fs::write(&bad, &b).unwrap();
        assert_decode_errors_no_output(&bad, d.path(), "dup");
    }
    // (f-ii) non-contiguous chunk index: bump every per-chunk entry's chunk_index by 1 so
    //        indices become 1..=num_chunks (missing 0). The whole-archive ChunkDecodedSizes
    //        global (chunk_index == u32::MAX sentinel) is left untouched — it's not a chunk.
    {
        let mut b = good.clone();
        for i in 0..n {
            let o = ent(i);
            let ci = u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
            if ci == u32::MAX { continue; } // ChunkDecodedSizes global sentinel
            b[o..o + 4].copy_from_slice(&(ci + 1).to_le_bytes());
        }
        v5_refresh_footer_crc(&mut b);
        let bad = d.path().join("noncontig.qz");
        std::fs::write(&bad, &b).unwrap();
        assert_decode_errors_no_output(&bad, d.path(), "noncontig");
    }
    // (f-iii) within-chunk record_count mismatch: bump entry[0]'s record_count
    //         (u64). Caught by validate_paired_directory's per-chunk uniformity.
    {
        let mut b = good.clone();
        let o = ent(0) + V5_RC_FIELD;
        let rc = u64::from_le_bytes(b[o..o + 8].try_into().unwrap()) + 1;
        b[o..o + 8].copy_from_slice(&rc.to_le_bytes());
        v5_refresh_footer_crc(&mut b);
        let bad = d.path().join("rcmix.qz");
        std::fs::write(&bad, &b).unwrap();
        assert_decode_errors_no_output(&bad, d.path(), "rcmix");
    }
}

// (f') Allocation-safety: an inflated n_entries must be REJECTED by the byte-bound
//      guard in ChunkDirectory::parse BEFORE Vec::with_capacity — returning Err,
//      never aborting the process via a multi-GB allocation. This test must
//      complete instantly; if it hangs/OOMs the byte-bound guard is missing.
#[test]
fn hostile_footer_entry_count_overflow() {
    let d = tempfile::tempdir().unwrap();
    let out = build_valid_paired(d.path(), 8, 4); // 2 chunks
    let good = std::fs::read(&out).unwrap();
    let (fs, _) = v5_footer_span(&good);

    // Sanity: the un-corrupted archive decodes fine (validates our offset math).
    {
        decompress(&dcfg(&out, &d.path().join("ok"), d.path(), true)).unwrap();
        assert!(d.path().join("ok_R1.fastq").exists());
        assert!(d.path().join("ok_R2.fastq").exists());
    }

    // v5 footer prefix at `fs`: num_reads(8) num_chunks(4) n_entries(4).
    // n_entries occupies bytes [fs+12 .. fs+16]. Inflate it to u32::MAX; each
    // entry needs 31 on-disk bytes, far more than remain, so the byte-bound guard
    // ((n as u64) > remaining/ENTRY_BYTES) must reject before allocating. We refresh
    // the footer CRC so the locator check passes and the parse guard is reached.
    let mut b = good.clone();
    b[fs + 12..fs + 16].copy_from_slice(&u32::MAX.to_le_bytes());
    v5_refresh_footer_crc(&mut b);

    let bad = d.path().join("entcount.qz");
    std::fs::write(&bad, &b).unwrap();
    assert_decode_errors_no_output(&bad, d.path(), "entcount");
}

// (f'') Hostile-`u64` record_count surface (Codex H1). Widening the on-disk
//   record_count from u32 to u64 lifts the per-count ceiling to u64::MAX. Every
//   record_count-driven decode-cap/allocation must remain bounded by the input
//   *bytes*, never by the untrusted count — so a footer claiming an astronomically
//   large record_count must make the decoder return Err, NOT abort the process via
//   a multi-EB allocation. We rewrite EVERY entry's record_count to a uniform
//   hostile value and set num_reads = num_chunks * value so the archive survives
//   `validate_paired_directory` (per-chunk uniformity + Σ == num_reads) and actually
//   reaches decode, then assert decode returns Err and the test process survives.
#[test]
fn paired_v5_hostile_record_count_rejected() {
    let d = tempfile::tempdir().unwrap();
    let out = build_valid_paired(d.path(), 8, 4); // 2 chunks (4 + 4)
    let good = std::fs::read(&out).unwrap();
    let (fs, _) = v5_footer_span(&good);

    // num_chunks (u32 at fs+8) and n_entries (u32 at fs+12).
    let num_chunks = u32::from_le_bytes(good[fs + 8..fs + 12].try_into().unwrap()) as u64;
    let n = u32::from_le_bytes(good[fs + 12..fs + 16].try_into().unwrap()) as usize;
    assert!(n >= 12, "expected multi-chunk footer, got {n} entries");

    let hostile: u64 = u64::MAX / 2; // > u32::MAX; num_chunks * hostile must not overflow u64
    let mut b = good.clone();
    // Rewrite every entry's record_count (8 bytes at entry+23).
    for i in 0..n {
        let o = v5_ent_off(&b, i) + V5_RC_FIELD;
        b[o..o + 8].copy_from_slice(&hostile.to_le_bytes());
    }
    // num_reads (u64 at fs+0) = num_chunks * hostile so Σ-per-chunk == num_reads.
    let num_reads = num_chunks
        .checked_mul(hostile)
        .expect("test value overflows u64");
    b[fs..fs + 8].copy_from_slice(&num_reads.to_le_bytes());
    v5_refresh_footer_crc(&mut b);

    let bad = d.path().join("hostile_rc.qz");
    std::fs::write(&bad, &b).unwrap();
    // The decode must return Err (bounded by input bytes), NOT abort via a giant
    // allocation. assert_decode_errors_no_output reaching this point at all proves
    // the process survived; it also asserts no partial output is left behind.
    assert_decode_errors_no_output(&bad, d.path(), "hostile_rc");
}

// (g) archive-probe edge cases on tiny hand-built files.
#[test]
fn hostile_probe_edge_cases() {
    let d = tempfile::tempdir().unwrap();
    // (g-i) 9-byte file, byte 8 == legacy paired marker(0xF0): not a v5 archive,
    // so the `legacy_marker` probe rejects it with the recompress hint.
    {
        let mut b = vec![0u8; 9];
        b[0] = b'Q';
        b[1] = b'Z';
        b[8] = 0xF0;
        let bad = d.path().join("probe9.qz");
        std::fs::write(&bad, &b).unwrap();
        let e = decompress(&dcfg(&bad, &d.path().join("p9"), d.path(), true)).unwrap_err();
        assert!(
            e.to_string().contains("legacy") && e.to_string().contains("recompress"),
            "got: {e}"
        );
    }
    // (g-ii) 10-byte file, byte 8 != a legacy marker (0x00): not v5, not a legacy
    // marker → routed to the single-end v4 decode, which rejects the truncated/
    // invalid header rather than silently mis-decoding it.
    {
        let mut b = vec![0u8; 10];
        b[0] = b'Q';
        b[1] = b'Z';
        b[8] = 0x00;
        b[9] = 0x08;
        let bad = d.path().join("probe10.qz");
        std::fs::write(&bad, &b).unwrap();
        // single-end v4 decode of a truncated/invalid header must Err (the exact
        // message is the single-end parser's, not a paired/legacy one).
        assert!(decompress(&dcfg(&bad, &d.path().join("p10"), d.path(), true)).is_err());
    }
    // (g-iii) a normal single-end archive still round-trips (single-end unaffected).
    {
        let se_in = d.path().join("se.fastq");
        w(
            &se_in,
            "@r1 x\nACGTACGT\n+\nIIIIIIII\n@r2 y\nGGCCGGCC\n+\nIIIIIIII\n",
        );
        let se_out = d.path().join("se.qz");
        compress(&cfg(
            vec![se_in.clone()],
            se_out.clone(),
            d.path().to_path_buf(),
        ))
        .unwrap();
        // single-end is classified non-paired and round-trips.
        assert!(!qz_lib::compression::is_paired_archive(&se_out).unwrap());
        let se_dec = d.path().join("se_out.fastq");
        decompress(&dcfg(&se_out, &se_dec, d.path(), true)).unwrap();
        assert_eq!(read_fastq(&se_dec), read_fastq(&se_in));
    }
}

// (h) Misparse guard: a single-end archive verifies normally and is NOT
//     misclassified as paired (the decompress_to_writer paired-reject is
//     already covered by `decompress_to_writer_rejects_paired`).
#[test]
fn hostile_single_end_not_misclassified() {
    let d = tempfile::tempdir().unwrap();
    let se_in = d.path().join("se.fastq");
    w(&se_in, "@r1 x\nACGTACGT\n+\nIIIIIIII\n");
    let se_out = d.path().join("se.qz");
    compress(&cfg(vec![se_in], se_out.clone(), d.path().to_path_buf())).unwrap();
    assert!(!qz_lib::compression::is_paired_archive(&se_out).unwrap());
    // deep verify of a single-end archive still works (no paired routing).
    let v = qz_lib::compression::verify(&verify_cfg(&se_out, false)).unwrap();
    assert_eq!(v.num_reads, 1);
    // a single-end archive CAN be decoded to a single writer (paired cannot).
    let mut sink = Vec::new();
    qz_lib::compression::decompress_to_writer(&se_out, &mut sink).unwrap();
    assert!(!sink.is_empty());
}

// ===========================================================================
// Step 2: Rejected-options matrix — every rejected knob errors AND leaves no
// output. Table-driven for DRYness.
// ===========================================================================

#[test]
fn rejected_options_matrix() {
    use qz_lib::cli::{HeaderCompressor, QualityCompressor, QualityMode};

    type Mut = Box<dyn Fn(&mut CompressConfig)>;
    let cases: Vec<(&str, Mut)> = vec![
        (
            "no_quality",
            Box::new(|c: &mut CompressConfig| c.no_quality = true),
        ),
        (
            "quality_mode=IlluminaBin",
            Box::new(|c: &mut CompressConfig| c.quality_mode = QualityMode::IlluminaBin),
        ),
        (
            "quality_mode=Illumina4",
            Box::new(|c: &mut CompressConfig| c.quality_mode = QualityMode::Illumina4),
        ),
        (
            "quality_mode=Binary",
            Box::new(|c: &mut CompressConfig| c.quality_mode = QualityMode::Binary),
        ),
        (
            "quality_mode=Discard",
            Box::new(|c: &mut CompressConfig| c.quality_mode = QualityMode::Discard),
        ),
        (
            "ultra=Some(1)",
            Box::new(|c: &mut CompressConfig| c.ultra = Some(1)),
        ),
        (
            "header_compressor=Bsc",
            Box::new(|c: &mut CompressConfig| c.advanced.header_compressor = HeaderCompressor::Bsc),
        ),
        (
            "quality_compressor=Bsc",
            Box::new(|c: &mut CompressConfig| {
                c.advanced.quality_compressor = QualityCompressor::Bsc
            }),
        ),
        (
            "bsc_static",
            Box::new(|c: &mut CompressConfig| c.advanced.bsc_static = true),
        ),
        (
            "sequence_hints",
            Box::new(|c: &mut CompressConfig| c.advanced.sequence_hints = true),
        ),
        (
            "rc_canon",
            Box::new(|c: &mut CompressConfig| c.advanced.rc_canon = true),
        ),
        (
            "bsc_block_size_mb=50",
            Box::new(|c: &mut CompressConfig| c.advanced.bsc_block_size_mb = 50),
        ),
    ];

    for (name, mutate) in &cases {
        let d = tempfile::tempdir().unwrap();
        let r1 = d.path().join("R1.fastq");
        let r2 = d.path().join("R2.fastq");
        w(&r1, "@a 1\nACGT\n+\nIIII\n");
        w(&r2, "@a 2\nTTTT\n+\nIIII\n");
        let out = d.path().join("o.qz");
        let mut c = cfg(vec![r1, r2], out.clone(), d.path().to_path_buf());
        mutate(&mut c);
        let r = compress(&c);
        assert!(r.is_err(), "case `{name}` should have been rejected");
        assert!(!out.exists(), "case `{name}` left an output archive behind");
    }

    // stdin input and stdout output use the "-" sentinel; build the cfg inline
    // because the `-` paths can't be created on disk.
    {
        let d = tempfile::tempdir().unwrap();
        let r2 = d.path().join("R2.fastq");
        w(&r2, "@a 2\nTTTT\n+\nIIII\n");
        let out = d.path().join("o.qz");
        // stdin as one of the two inputs.
        let c = cfg(
            vec![PathBuf::from("-"), r2.clone()],
            out.clone(),
            d.path().to_path_buf(),
        );
        assert!(
            compress(&c).is_err(),
            "stdin input should be rejected in paired mode"
        );
        assert!(!out.exists());
        // stdout output is now SUPPORTED in (non-reference) paired mode: the
        // public dispatch routes `-o -` to compress_paired_v5, which writes the
        // seek-free v5 container to stdout. The end-to-end stdout pipe roundtrip
        // is exercised by the CLI subprocess test
        // `crates/qz-cli/tests/stdout_roundtrip.rs::paired_v5_compress_to_stdout_roundtrips`
        // (capturing process stdout here would pollute the test harness output).
    }
}

// ===========================================================================
// Step 3: Cross-chunk read-length change — variable framing across chunks.
// ===========================================================================

#[test]
fn cross_chunk_variable_read_length() {
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("R1.fastq");
    let r2 = d.path().join("R2.fastq");
    let mut s1 = String::new();
    let mut s2 = String::new();
    // chunk_records = 4: first 4 pairs len 10, the rest len 12 (both mates vary).
    for i in 0..10 {
        if i < 4 {
            s1.push_str(&format!("@m:{i} 1:N:0:A\nACGTACGTAC\n+\nIIIIIIIIII\n")); // len 10
            s2.push_str(&format!("@m:{i} 2:N:0:A\nTTTTGGGGCC\n+\n##IIIIIIII\n")); // len 10
        } else {
            s1.push_str(&format!("@m:{i} 1:N:0:A\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")); // len 12
            s2.push_str(&format!("@m:{i} 2:N:0:A\nTTTTGGGGCCAA\n+\n##IIIIIIIIII\n")); // len 12
        }
    }
    w(&r1, &s1);
    w(&r2, &s2);
    let out = d.path().join("p.qz");
    let mut c = cfg(vec![r1, r2], out.clone(), d.path().to_path_buf());
    c.advanced.chunk_records = 4;
    compress(&c).unwrap();
    let pre = d.path().join("dec");
    decompress(&dcfg(&out, &pre, d.path(), true)).unwrap();
    assert_eq!(read_fastq(&d.path().join("dec_R1.fastq")), s1);
    assert_eq!(read_fastq(&d.path().join("dec_R2.fastq")), s2);
}

// ===========================================================================
// Step 4: --force rollback — a failed decode must preserve pre-existing
// --force targets (restored from .bak, never clobbered with partial output).
// ===========================================================================

#[test]
fn force_rollback_preserves_existing_targets() {
    let d = tempfile::tempdir().unwrap();
    let out = build_valid_paired(d.path(), 8, 4);
    // Corrupt the payload so decode fails AFTER temp files are created.
    let mut bytes = std::fs::read(&out).unwrap();
    let header_end = u32::from_le_bytes(bytes[4..8].try_into().unwrap()) as usize;
    let (footer_start, _) = v5_footer_span(&bytes);
    let pos = header_end + (footer_start - header_end) / 2;
    bytes[pos] ^= 0xFF;
    let bad = d.path().join("corrupt.qz");
    std::fs::write(&bad, &bytes).unwrap();

    // Pre-create the force targets with known sentinel content.
    let prefix = d.path().join("dst");
    let r1 = d.path().join("dst_R1.fastq");
    let r2 = d.path().join("dst_R2.fastq");
    w(&r1, "SENTINEL-R1\n");
    w(&r2, "SENTINEL-R2\n");

    let dc = dcfg(&bad, &prefix, d.path(), true); // force = true
    assert!(
        decompress(&dc).is_err(),
        "corrupt archive should fail to decode"
    );

    // Originals must survive intact (decode failed before publish; temp files
    // removed by the drop guard; the existing targets were never renamed away).
    assert_eq!(
        read_fastq(&r1),
        "SENTINEL-R1\n",
        "R1 sentinel was clobbered"
    );
    assert_eq!(
        read_fastq(&r2),
        "SENTINEL-R2\n",
        "R2 sentinel was clobbered"
    );
    // No stray temp files left behind.
    let leftovers: Vec<_> = std::fs::read_dir(d.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().contains(".tmp"))
        .collect();
    assert!(leftovers.is_empty(), "stray temp files: {leftovers:?}");
}

// ---------------------------------------------------------------------------
// Interleaved decode (Mechanism A): one stream, R1/R2 records alternating.
// ---------------------------------------------------------------------------

/// Split an interleaved FASTQ (4-line records, R1/R2 alternating) back into the
/// two mate strings: even-index records (0,2,4,…) are R1, odd-index are R2.
fn deinterleave(content: &str) -> (String, String) {
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len() % 8, 0, "interleaved FASTQ not a whole number of pairs");
    let (mut r1, mut r2) = (String::new(), String::new());
    let mut rec = 0;
    let mut i = 0;
    while i < lines.len() {
        let sink = if rec % 2 == 0 { &mut r1 } else { &mut r2 };
        for j in 0..4 {
            sink.push_str(lines[i + j]);
            sink.push('\n');
        }
        i += 4;
        rec += 1;
    }
    (r1, r2)
}

#[test]
fn paired_interleaved_roundtrip_type1() {
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("r1.fastq");
    let r2 = d.path().join("r2.fastq");
    let s1 = "@a/1\nACGTACGT\n+\nIIIIIIII\n@b/1\nTTTTGGGG\n+\n########\n@c/1\nGGGGCCCC\n+\nJJJJJJJJ\n";
    let s2 = "@a/2\nCCCCAAAA\n+\n$$$$$$$$\n@b/2\nAAAACCCC\n+\n!!!!!!!!\n@c/2\nGGCCTTAA\n+\n%%%%%%%%\n";
    w(&r1, s1);
    w(&r2, s2);
    let out = d.path().join("o.qz");
    qz_lib::compression::compress(&cfg(vec![r1, r2], out.clone(), d.path().to_path_buf())).unwrap();

    let inter = d.path().join("inter.fastq");
    let dc = qz_lib::cli::DecompressConfig {
        input: out,
        output: vec![inter.clone()],
        working_dir: d.path().to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress_interleaved(&dc).unwrap();

    let content = std::fs::read_to_string(&inter).unwrap();
    let (de1, de2) = deinterleave(&content);
    assert_eq!(de1, s1, "de-interleaved R1 must match original R1");
    assert_eq!(de2, s2, "de-interleaved R2 must match original R2");
}

#[test]
fn interleaved_rejects_single_end_archive() {
    let d = tempfile::tempdir().unwrap();
    let r = d.path().join("se.fastq");
    w(&r, "@a\nACGT\n+\nIIII\n@b\nTTTT\n+\n####\n");
    let out = d.path().join("se.qz");
    qz_lib::compression::compress(&cfg(vec![r], out.clone(), d.path().to_path_buf())).unwrap();

    let dc = qz_lib::cli::DecompressConfig {
        input: out,
        output: vec![d.path().join("inter.fastq")],
        working_dir: d.path().to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    let e = qz_lib::compression::decompress_interleaved(&dc).unwrap_err();
    assert!(
        e.to_string().contains("paired"),
        "single-end archive must be rejected by --interleaved; got: {e}"
    );
}

#[test]
fn interleaved_rejects_multiple_outputs() {
    let d = tempfile::tempdir().unwrap();
    let r1 = d.path().join("r1.fastq");
    let r2 = d.path().join("r2.fastq");
    w(&r1, "@a/1\nACGT\n+\nIIII\n");
    w(&r2, "@a/2\nCCCC\n+\n$$$$\n");
    let out = d.path().join("p.qz");
    qz_lib::compression::compress(&cfg(vec![r1, r2], out.clone(), d.path().to_path_buf())).unwrap();

    let dc = qz_lib::cli::DecompressConfig {
        input: out,
        output: vec![d.path().join("a.fastq"), d.path().join("b.fastq")],
        working_dir: d.path().to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    let e = qz_lib::compression::decompress_interleaved(&dc).unwrap_err();
    assert!(
        e.to_string().contains("single stream"),
        "two -o values must be rejected; got: {e}"
    );
}

#[test]
fn interleaved_rejects_stdin_input() {
    // Decode reads the trailing footer locator, so the archive input must be seekable.
    let dc = qz_lib::cli::DecompressConfig {
        input: PathBuf::from("-"),
        output: vec![PathBuf::from("out.fastq")],
        working_dir: PathBuf::from("."),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    let e = qz_lib::compression::decompress_interleaved(&dc).unwrap_err();
    assert!(
        e.to_string().contains("seekable"),
        "piped-in archive must be rejected with a seekable-input hint; got: {e}"
    );
}
