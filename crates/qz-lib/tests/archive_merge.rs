use qz_lib::cli::{
    AdvancedOptions, CompressConfig, DecompressConfig, QualityCompressor, ReferenceOptions,
};
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

/// Build a deterministic FASTQ string of `n` reads, each `len` bp.
fn make_fastq(n: usize, len: usize) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut s = String::new();
    for i in 0..n {
        s.push_str(&format!("@read{i}\n"));
        let seq: String = (0..len).map(|j| bases[(i + j) % 4] as char).collect();
        s.push_str(&seq);
        s.push('\n');
        s.push_str("+\n");
        // Vary quality a little so fqz/bsc have something to model.
        let qual: String = (0..len).map(|j| (b'I' - ((i + j) % 5) as u8) as char).collect();
        s.push_str(&qual);
        s.push('\n');
    }
    s
}

/// Compress `fastq_text` to `archive` with a small chunk size so multi-chunk
/// archives are reachable in tests. Returns the archive path.
fn compress_to(dir: &Path, fastq_text: &str, name: &str, chunk_records: usize) -> PathBuf {
    let input = dir.join(format!("{name}.fastq"));
    fs::write(&input, fastq_text).unwrap();
    let archive = dir.join(format!("{name}.qz"));
    let cfg = CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: dir.to_path_buf(),
        threads: 1,
        advanced: AdvancedOptions { chunk_records, ..AdvancedOptions::default() },
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();
    archive
}

/// Compress with a caller-customized config (base: input/output/working_dir/threads=1).
fn compress_cfg(dir: &Path, fastq_text: &str, name: &str, f: impl FnOnce(&mut CompressConfig)) -> PathBuf {
    let input = dir.join(format!("{name}.fastq"));
    fs::write(&input, fastq_text).unwrap();
    let archive = dir.join(format!("{name}.qz"));
    let mut cfg = CompressConfig {
        input: vec![input], output: archive.clone(), working_dir: dir.to_path_buf(),
        threads: 1, ..CompressConfig::default()
    };
    f(&mut cfg);
    qz_lib::compression::compress(&cfg).unwrap();
    archive
}

/// Decompress `archive` and return its FASTQ text.
fn decompress_text(dir: &Path, archive: &Path, name: &str) -> String {
    let out = dir.join(format!("{name}.out.fastq"));
    let cfg = DecompressConfig {
        input: archive.to_path_buf(),
        output: vec![out.clone()],
        working_dir: dir.to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&cfg).unwrap();
    fs::read_to_string(&out).unwrap()
}

/// Split FASTQ text into contiguous record-slices at the given record-count
/// `boundaries` (cumulative; uneven sizes allowed). E.g. `&[6, 10]` → reads
/// [0..6) and [6..10).
fn split_records(fastq_text: &str, boundaries: &[usize]) -> Vec<String> {
    let lines: Vec<&str> = fastq_text.lines().collect();
    let recs: Vec<String> = lines.chunks(4).map(|c| c.join("\n")).collect();
    let mut out = Vec::new();
    let mut start = 0;
    for &end in boundaries {
        let mut blob = recs[start..end].join("\n");
        blob.push('\n');
        out.push(blob);
        start = end;
    }
    out
}

/// Merge `inputs` into `out_path` via the library function.
fn merge_to(inputs: &[PathBuf], out_path: &Path) {
    let mut out = std::io::BufWriter::new(fs::File::create(out_path).unwrap());
    qz_lib::compression::merge_single_end_archives(inputs, &mut out).unwrap();
    use std::io::Write;
    out.flush().unwrap();
}

#[test]
fn merge_two_shards_roundtrips_lossless() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let whole = make_fastq(10, 12);

    // Reference: compress the whole file.
    let whole_arc = compress_to(d, &whole, "whole", 4);
    let whole_text = decompress_text(d, &whole_arc, "whole");

    // Shards: split records [0..6) and [6..10), compress each with chunk_records=4
    // so each shard is multi-chunk, then merge.
    let slices = split_records(&whole, &[6, 10]);
    let s0 = compress_to(d, &slices[0], "s0", 4);
    let s1 = compress_to(d, &slices[1], "s1", 4);
    let merged = d.join("merged.qz");
    merge_to(&[s0, s1], &merged);

    let merged_text = decompress_text(d, &merged, "merged");
    assert_eq!(whole_text, merged_text, "merged archive must decode to the same FASTQ");
}

#[test]
fn merge_single_input_is_identity() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let whole = make_fastq(7, 10);
    let arc = compress_to(d, &whole, "only", 3);
    let whole_text = decompress_text(d, &arc, "only");

    let merged = d.join("merged1.qz");
    merge_to(&[arc], &merged);
    let merged_text = decompress_text(d, &merged, "merged1");
    assert_eq!(whole_text, merged_text);
}

#[test]
fn merge_three_uneven_shards_roundtrips() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let whole = make_fastq(13, 15);
    let whole_arc = compress_to(d, &whole, "w3", 4);
    let whole_text = decompress_text(d, &whole_arc, "w3");

    // Uneven slices: 2 (sub-chunk), 7 (multi-chunk), 4.
    let slices = split_records(&whole, &[2, 9, 13]);
    let s0 = compress_to(d, &slices[0], "t0", 4);
    let s1 = compress_to(d, &slices[1], "t1", 4);
    let s2 = compress_to(d, &slices[2], "t2", 4);
    let merged = d.join("merged3.qz");
    merge_to(&[s0, s1, s2], &merged);
    let merged_text = decompress_text(d, &merged, "merged3");
    assert_eq!(whole_text, merged_text);
}

#[test]
fn merged_properties_sum_reads_and_chunks() {
    use qz_lib::compression::read_chunk_layout;
    let td = TempDir::new().unwrap();
    let d = td.path();
    let whole = make_fastq(10, 12);
    let slices = split_records(&whole, &[6, 10]);
    let s0 = compress_to(d, &slices[0], "p0", 4); // 6 reads → 2 chunks
    let s1 = compress_to(d, &slices[1], "p1", 4); // 4 reads → 1 chunk
    let l0 = read_chunk_layout(&s0).unwrap();
    let l1 = read_chunk_layout(&s1).unwrap();
    let merged = d.join("mergedp.qz");
    merge_to(&[s0, s1], &merged);
    let lm = read_chunk_layout(&merged).unwrap();
    assert_eq!(lm.num_reads, l0.num_reads + l1.num_reads);
    assert_eq!(lm.num_chunks, l0.num_chunks + l1.num_chunks);
}

use qz_lib::cli::VerifyConfig;

fn verify_cfg(archive: &Path, dir: &Path, fast: bool) -> VerifyConfig {
    VerifyConfig {
        input: archive.to_path_buf(),
        working_dir: dir.to_path_buf(),
        num_threads: 1,
        fast,
    }
}

#[test]
fn merged_archive_passes_verify_deep_and_fast() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let whole = make_fastq(11, 14);
    let slices = split_records(&whole, &[5, 11]);
    let s0 = compress_to(d, &slices[0], "v0", 4);
    let s1 = compress_to(d, &slices[1], "v1", 4);
    let merged = d.join("vmerged.qz");
    merge_to(&[s0, s1], &merged);

    // Deep verify (full decompress + FASTQ CRC) and fast verify (per-block CRC walk).
    qz_lib::compression::verify(&verify_cfg(&merged, d, false)).unwrap();
    qz_lib::compression::verify(&verify_cfg(&merged, d, true)).unwrap();
}

#[test]
fn merged_decoded_sizes_table_is_concat_of_shards() {
    use qz_lib::compression::read_chunk_layout;
    let td = TempDir::new().unwrap();
    let d = td.path();
    let whole = make_fastq(12, 13);
    let slices = split_records(&whole, &[5, 12]);
    let s0 = compress_to(d, &slices[0], "z0", 4);
    let s1 = compress_to(d, &slices[1], "z1", 4);
    // read_chunk_layout returns the VALIDATED (CRC-checked) per-chunk decoded sizes.
    let sizes0 = read_chunk_layout(&s0).unwrap().per_chunk_decoded_bytes;
    let sizes1 = read_chunk_layout(&s1).unwrap().per_chunk_decoded_bytes;

    let merged = d.join("zmerged.qz");
    merge_to(&[s0, s1], &merged);
    let sizesm = read_chunk_layout(&merged).unwrap().per_chunk_decoded_bytes;

    let mut expected = sizes0.clone();
    expected.extend_from_slice(&sizes1);
    assert_eq!(sizesm, expected, "merged table must be the shards' tables concatenated");
}

#[test]
fn merge_paired_rebuilds_per_mate_decoded_sizes_table() {
    use qz_lib::compression::{merge_archives, read_chunk_layout};
    let td = TempDir::new().unwrap();
    let d = td.path();
    // Two paired shards; R1 and R2 use DIFFERENT read lengths so the per-mate decoded
    // sizes genuinely differ (a single combined table couldn't reproduce them).
    let r1_0 = make_fastq(6, 12);
    let r2_0 = make_fastq(6, 20);
    let r1_1 = make_fastq(5, 12);
    let r2_1 = make_fastq(5, 20);
    let p0 = compress_paired_to(d, &r1_0, &r2_0, "pp0");
    let p1 = compress_paired_to(d, &r1_1, &r2_1, "pp1");

    let pm0 = read_chunk_layout(&p0).unwrap().decoded_sizes_per_mate;
    let pm1 = read_chunk_layout(&p1).unwrap().decoded_sizes_per_mate;
    assert_eq!(pm0.len(), 2, "paired shard must carry a 2-mate table");
    assert_eq!(pm1.len(), 2);

    let merged = d.join("pp_merged.qz");
    {
        let mut w = std::io::BufWriter::new(fs::File::create(&merged).unwrap());
        merge_archives(&[p0, p1], &mut w).unwrap();
        use std::io::Write;
        w.flush().unwrap();
    }

    // The merged per-mate table must be each shard's per-mate table concatenated — so the
    // merged paired archive stays NUMA direct-write-decodable.
    let pmm = read_chunk_layout(&merged).unwrap().decoded_sizes_per_mate;
    assert_eq!(pmm.len(), 2, "merged paired must keep a 2-mate table");
    for m in 0..2 {
        let mut expected = pm0[m].clone();
        expected.extend_from_slice(&pm1[m]);
        assert_eq!(pmm[m], expected, "merged per-mate[{m}] table != shards concatenated");
    }

    // And it decodes losslessly to the per-mate concatenation of the shards' reads.
    let prefix = d.join("pp_dec");
    qz_lib::compression::decompress(&DecompressConfig {
        input: merged,
        output: vec![prefix.clone()],
        working_dir: d.to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    })
    .unwrap();
    let sfx = |s: &str| {
        let mut p = prefix.as_os_str().to_owned();
        p.push(s);
        PathBuf::from(p)
    };
    assert_eq!(fs::read_to_string(sfx("_R1.fastq")).unwrap(), format!("{r1_0}{r1_1}"));
    assert_eq!(fs::read_to_string(sfx("_R2.fastq")).unwrap(), format!("{r2_0}{r2_1}"));
}

/// Try to merge; return the Result so tests can assert on errors.
fn try_merge(inputs: &[PathBuf], out_path: &Path) -> anyhow::Result<()> {
    let mut out = std::io::BufWriter::new(fs::File::create(out_path).unwrap());
    let r = qz_lib::compression::merge_single_end_archives(inputs, &mut out);
    use std::io::Write;
    let _ = out.flush();
    r
}

#[test]
fn merge_rejects_incompatible_front_headers() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    // One archive with all-equal-length reads (const-lengths flag set), one with
    // variable-length reads (flag clear) → the two front headers differ.
    let const_len = make_fastq(6, 12);
    // Prepend one shorter read so the variable length lands in the FIRST chunk
    // (chunk_records = 4): the const-lengths profile is detected from chunk 0, so a
    // short read there clears the flag and the front header differs. (Appending it
    // instead would leave chunk 0 all-equal-length — flag set — and the later short
    // read would hit the const-length validation bail during compression, not merge.)
    let mut var_len = String::from("@readX\nACGTACGT\n+\nIIIIIIII\n");
    var_len.push_str(&make_fastq(5, 12));
    let a = compress_to(d, &const_len, "cl", 4);
    let b = compress_to(d, &var_len, "vl", 4);
    let out = d.join("bad.qz");
    let err = try_merge(&[a, b], &out).unwrap_err().to_string();
    assert!(err.contains("incompatible archive configs"), "got: {err}");
}

#[test]
fn merge_rejects_empty_input_list() {
    let td = TempDir::new().unwrap();
    let out = td.path().join("none.qz");
    let err = try_merge(&[], &out).unwrap_err().to_string();
    assert!(err.contains("at least one input"), "got: {err}");
}

#[test]
fn merge_rejects_corrupt_input() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let whole = make_fastq(6, 12);
    let a = compress_to(d, &whole, "ok", 4);
    // Truncate the archive: footer/locator gone → read_v5_footer must fail.
    let bytes = fs::read(&a).unwrap();
    let corrupt = d.join("corrupt.qz");
    fs::write(&corrupt, &bytes[..bytes.len() / 2]).unwrap();
    let out = d.join("c.qz");
    // Pin the failure to the footer/locator read (not an unrelated error): a
    // half-truncated archive loses the trailing locator → "bad footer magic".
    let err = try_merge(&[corrupt], &out).unwrap_err().to_string();
    assert!(err.contains("footer"), "expected a footer-read failure, got: {err}");
}

/// Compress two FASTQ files into a paired (archive_type 1) archive.
fn compress_paired_to(dir: &Path, r1_text: &str, r2_text: &str, name: &str) -> PathBuf {
    let r1 = dir.join(format!("{name}_1.fastq"));
    let r2 = dir.join(format!("{name}_2.fastq"));
    fs::write(&r1, r1_text).unwrap();
    fs::write(&r2, r2_text).unwrap();
    let archive = dir.join(format!("{name}.qz"));
    let cfg = CompressConfig {
        input: vec![r1, r2],
        output: archive.clone(),
        working_dir: dir.to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();
    archive
}

#[test]
fn merge_rejects_paired_input() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    // Paired R1/R2 with equal record counts → archive_type 1; merge is single-end only.
    let r1 = make_fastq(50, 12);
    let r2 = make_fastq(50, 12);
    let p = compress_paired_to(d, &r1, &r2, "paired");
    let out = d.join("p_out.qz");
    let err = try_merge(&[p], &out).unwrap_err().to_string();
    assert!(err.contains("is not a single-end archive"), "got: {err}");
}

#[test]
fn merge_to_path_force_and_reject() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let whole = make_fastq(8, 12);
    let a = compress_to(d, &whole, "tp", 4);
    let out = d.join("merged_tp.qz");

    // First write succeeds.
    qz_lib::compression::merge_single_end_archives_to_path(std::slice::from_ref(&a), &out, false)
        .unwrap();
    let text1 = decompress_text(d, &out, "tp1");

    // Second write without force is refused and leaves the file intact.
    let err = qz_lib::compression::merge_single_end_archives_to_path(std::slice::from_ref(&a), &out, false)
        .unwrap_err()
        .to_string();
    assert!(err.contains("already exists"), "got: {err}");
    assert_eq!(decompress_text(d, &out, "tp2"), text1, "refused merge must not change the file");

    // With force it overwrites (still decodes correctly).
    qz_lib::compression::merge_single_end_archives_to_path(&[a], &out, true).unwrap();
    assert_eq!(decompress_text(d, &out, "tp3"), text1);
}

/// Compress the whole input + each shard with the SAME custom config, merge, and
/// assert the merged archive decodes byte-identically to the whole-file archive,
/// then deep- AND fast-verify the merged archive. The lossless-oracle pattern,
/// parameterized over a config customizer so production archive shapes (fqzcomp
/// quality, no-quality, RcCanon) are exercised end-to-end through merge.
fn merge_shape_roundtrips(
    n: usize,
    len: usize,
    boundaries: &[usize],
    tag: &str,
    customize: impl Fn(&mut CompressConfig) + Copy,
) {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let whole = make_fastq(n, len);

    let whole_arc = compress_cfg(d, &whole, &format!("{tag}_whole"), customize);
    let whole_text = decompress_text(d, &whole_arc, &format!("{tag}_whole"));

    let slices = split_records(&whole, boundaries);
    let shards: Vec<PathBuf> = slices
        .iter()
        .enumerate()
        .map(|(i, s)| compress_cfg(d, s, &format!("{tag}_s{i}"), customize))
        .collect();

    let merged = d.join(format!("{tag}_merged.qz"));
    merge_to(&shards, &merged);

    let merged_text = decompress_text(d, &merged, &format!("{tag}_merged"));
    assert_eq!(whole_text, merged_text, "[{tag}] merged archive must decode to the whole-file FASTQ");

    // Both verify modes must accept the merged archive.
    qz_lib::compression::verify(&verify_cfg(&merged, d, false)).unwrap();
    qz_lib::compression::verify(&verify_cfg(&merged, d, true)).unwrap();
}

#[test]
fn merge_fqzcomp_quality_roundtrips() {
    // chunk_records=8, qblock=2 → a chunk spans MULTIPLE Qual sub-blocks: the
    // production NUMA-shard shape the merge must copy verbatim.
    merge_shape_roundtrips(8, 12, &[4, 8], "fqz", |cfg| {
        cfg.advanced.chunk_records = 8;
        cfg.advanced.quality_compressor = QualityCompressor::Fqzcomp;
        cfg.advanced.quality_ctx_block_size = 2;
    });
}

#[test]
fn merge_no_quality_roundtrips() {
    // No Qual role at all in the directory.
    merge_shape_roundtrips(10, 12, &[6, 10], "noqual", |cfg| {
        cfg.no_quality = true;
        cfg.advanced.chunk_records = 4;
    });
}

#[test]
fn merge_rc_canon_roundtrips() {
    // RcCanon → RcFlags role + the RcFlags-iff-has_rc_canon directory self-check.
    merge_shape_roundtrips(10, 12, &[6, 10], "rccanon", |cfg| {
        cfg.advanced.rc_canon = true;
        cfg.advanced.chunk_records = 4;
    });
}

#[test]
fn merge_archives_single_end_roundtrips_via_new_entry() {
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let text = make_fastq(40, 60);
    // `split_records` uses CUMULATIVE boundaries, so the final boundary must be the
    // total read count to cover the whole file: [0..16) [16..32) [32..40) → 3 shards.
    let parts = split_records(&text, &[16, 32, 40]); // 3 shards
    let a = compress_to(d, &parts[0], "a", 8);
    let b = compress_to(d, &parts[1], "b", 8);
    let c = compress_to(d, &parts[2], "c", 8);
    let out = d.join("merged.qz");
    let mut w = std::io::BufWriter::new(std::fs::File::create(&out).unwrap());
    qz_lib::compression::merge_archives(&[a, b, c], &mut w).unwrap();
    use std::io::Write;
    w.flush().unwrap();
    drop(w);
    assert_eq!(decompress_text(d, &out, "dec"), text, "merge_archives must roundtrip single-end losslessly");
}

#[test]
fn merge_archives_rejects_reference_archive_type() {
    // Synthesize an archive_type=2 front-header byte and confirm the dispatcher
    // rejects it by message (no reference fixture needed — we forge the type byte).
    let td = tempfile::TempDir::new().unwrap();
    let d = td.path();
    let a = compress_to(d, &make_fastq(8, 40), "a", 8);
    let mut bytes = std::fs::read(&a).unwrap();
    bytes[3] = 2; // archive_type byte (header byte 3)
    let forged = d.join("forged_ref.qz");
    std::fs::write(&forged, &bytes).unwrap();
    let out = d.join("o.qz");
    let mut w = std::io::BufWriter::new(std::fs::File::create(&out).unwrap());
    let err = qz_lib::compression::merge_archives(&[forged], &mut w).unwrap_err().to_string();
    assert!(err.contains("reference"), "reference merge must be rejected, got: {err}");
}

// ---------------------------------------------------------------------------
// Reference-aware merge (archive_type 2 = paired reference / 4 = single-end
// reference). Unlike single-end/paired merge, the coverage globals
// (PackedBacking/NBitmap/IntervalMap/ReferenceMeta) are re-derived over the UNION
// of the shards' covered intervals from the SAME reference FASTA — so the only
// oracle that matters is lossless roundtrip of the merged archive vs. the original
// reads (the shards' chunk frames relocate verbatim; the globals are rebuilt).
// ---------------------------------------------------------------------------

/// Deterministic low-repeat ACGT sequence (matches the reference suites' generator
/// so reads seed reliably against the mapper).
fn merge_make_seq(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        v.push(b"ACGT"[((x >> 33) & 3) as usize]);
    }
    v
}

/// Build a reference FASTA + a single-end reads FASTQ string: ~3 kb reference, N
/// forward substrings (100–120 bp), a few with planted substitutions, the last two
/// random (literal fallback) so both the mapped-reconstruct and FallbackPool copy
/// paths are exercised.
fn merge_ref_single_dataset(dir: &Path, n_reads: usize) -> (PathBuf, String) {
    let refseq = merge_make_seq(3000, 7);
    let rf = dir.join("mref.fa");
    {
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        fs::write(&rf, &s).unwrap();
    }
    let mut s1 = String::new();
    let mut rng = 0x1234_5678_9abc_def0u64;
    let mut next = || {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (rng >> 33) as usize
    };
    for i in 0..n_reads {
        let unmappable = i >= n_reads - 2;
        let rlen = 100 + (next() % 21);
        let rbytes = if unmappable {
            merge_make_seq(rlen, 9000 + i as u64)
        } else {
            let max_start = refseq.len() - rlen;
            let st = next() % (max_start + 1);
            let mut r = refseq[st..st + rlen].to_vec();
            if i % 17 == 0 {
                let p = next() % rlen;
                r[p] = if r[p] == b'A' { b'C' } else { b'A' };
            }
            r
        };
        let q: String = "I".repeat(rbytes.len());
        s1.push_str(&format!("@read_{i}\n{}\n+\n{q}\n", std::str::from_utf8(&rbytes).unwrap()));
    }
    (rf, s1)
}

/// Compress a single-end reads slice as a type-4 reference archive vs `fasta`.
fn compress_ref_single(dir: &Path, fastq_text: &str, name: &str, fasta: &Path, chunk_records: usize) -> PathBuf {
    let input = dir.join(format!("{name}.fastq"));
    fs::write(&input, fastq_text).unwrap();
    let archive = dir.join(format!("{name}.qz"));
    let cfg = CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: dir.to_path_buf(),
        threads: 1,
        force: true,
        reference: Some(ReferenceOptions {
            reference: fasta.to_path_buf(),
            reference_index: None,
            reference_fast: false,
            reference_window: 2,
        }),
        advanced: AdvancedOptions { chunk_records, ..AdvancedOptions::default() },
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();
    archive
}

/// Decode a type-4 reference archive (ONE bare output file) → FASTQ text.
fn decompress_ref_single(dir: &Path, archive: &Path, name: &str) -> String {
    let out = dir.join(format!("{name}.out.fastq"));
    let cfg = DecompressConfig {
        input: archive.to_path_buf(),
        output: vec![out.clone()],
        working_dir: dir.to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&cfg).unwrap();
    fs::read_to_string(&out).unwrap()
}

#[test]
fn merge_reference_single_end_roundtrips_lossless() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let (fasta, whole) = merge_ref_single_dataset(d, 200);

    // Reference: the whole read set in one type-4 archive.
    let whole_arc = compress_ref_single(d, &whole, "whole_ref", &fasta, 60);
    let whole_text = decompress_ref_single(d, &whole_arc, "whole_ref");
    assert_eq!(whole_text, whole, "reference roundtrip must be byte-identical (sanity)");

    // Shards: split records [0..120) and [120..200), compress each (multi-chunk via
    // chunk_records=60) vs the SAME FASTA, then reference-merge.
    let slices = split_records(&whole, &[120, 200]);
    let s0 = compress_ref_single(d, &slices[0], "rs0", &fasta, 60);
    let s1 = compress_ref_single(d, &slices[1], "rs1", &fasta, 60);
    let merged = d.join("merged_ref_single.qz");
    qz_lib::compression::merge_reference_archives_to_path(&[s0, s1], &fasta, &merged, true).unwrap();

    let merged_text = decompress_ref_single(d, &merged, "merged_ref_single");
    assert_eq!(whole, merged_text, "merged reference archive must decode to the original reads");

    // Deep + fast verify the merged archive (re-derived globals + relocated frames
    // are self-consistent).
    qz_lib::compression::verify(&verify_cfg(&merged, d, false)).unwrap();
    qz_lib::compression::verify(&verify_cfg(&merged, d, true)).unwrap();
}

#[test]
fn merge_reference_rejects_wrong_fasta() {
    // Reference merge re-derives the coverage backing from the supplied FASTA and decode
    // reconstructs bases from it (never re-reading a FASTA). Pointing the merge at the
    // WRONG reference must be rejected via the per-shard ReferenceMeta digest — otherwise
    // a corrupt archive is produced that still passes verify.
    let td = TempDir::new().unwrap();
    let d = td.path();
    let (fasta_a, whole) = merge_ref_single_dataset(d, 200);

    let slices = split_records(&whole, &[120, 200]);
    let s0 = compress_ref_single(d, &slices[0], "wrs0", &fasta_a, 60);
    let s1 = compress_ref_single(d, &slices[1], "wrs1", &fasta_a, 60);

    // A different reference with the SAME contig name + length but different bases, so
    // the front-header config still matches and only the digest can catch the mismatch.
    let fasta_b = d.join("mref_wrong.fa");
    {
        let other = merge_make_seq(3000, 999);
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&other).unwrap());
        s.push('\n');
        fs::write(&fasta_b, &s).unwrap();
    }

    let merged = d.join("merged_wrong.qz");
    let err = qz_lib::compression::merge_reference_archives_to_path(&[s0, s1], &fasta_b, &merged, true)
        .expect_err("merging against the wrong reference must be rejected")
        .to_string();
    assert!(
        err.contains("different reference") || err.contains("digest"),
        "expected a reference-mismatch error, got: {err}"
    );
}

/// Build a reference FASTA + paired R1/R2 reads. R1 = forward substring, R2 = a
/// forward substring downstream (insert ~250) so pairs map; the last two pairs are
/// random (fallback). Returns `(fasta, r1_fastq, r2_fastq)`.
fn merge_ref_paired_dataset(dir: &Path, n_pairs: usize) -> (PathBuf, String, String) {
    let refseq = merge_make_seq(4000, 11);
    let rf = dir.join("mrefp.fa");
    {
        let mut s = String::from(">chrp\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        fs::write(&rf, &s).unwrap();
    }
    let comp = |b: u8| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', _ => b'A' };
    let (mut s1, mut s2) = (String::new(), String::new());
    let mut rng = 0xfeed_face_dead_beefu64;
    let mut next = || {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (rng >> 33) as usize
    };
    for i in 0..n_pairs {
        let unmappable = i >= n_pairs - 2;
        let rlen = 100 + (next() % 21);
        let (r1, r2) = if unmappable {
            (merge_make_seq(rlen, 7000 + i as u64), merge_make_seq(rlen, 8000 + i as u64))
        } else {
            let insert = 250;
            let max_start = refseq.len() - (insert + rlen);
            let st = next() % (max_start + 1);
            let r1 = refseq[st..st + rlen].to_vec();
            // R2 = reverse-complement of a downstream window (typical FR orientation).
            let r2s = st + insert;
            let mut r2: Vec<u8> = refseq[r2s..r2s + rlen].iter().rev().map(|&b| comp(b)).collect();
            if i % 13 == 0 {
                let p = next() % rlen;
                r2[p] = if r2[p] == b'A' { b'C' } else { b'A' };
            }
            (r1, r2)
        };
        let q: String = "I".repeat(rlen);
        s1.push_str(&format!("@pair_{i}/1\n{}\n+\n{q}\n", std::str::from_utf8(&r1).unwrap()));
        s2.push_str(&format!("@pair_{i}/2\n{}\n+\n{q}\n", std::str::from_utf8(&r2).unwrap()));
    }
    (rf, s1, s2)
}

/// Compress a paired reads slice as a type-2 reference archive vs `fasta`.
fn compress_ref_paired(dir: &Path, r1: &str, r2: &str, name: &str, fasta: &Path, chunk_records: usize) -> PathBuf {
    let r1p = dir.join(format!("{name}_R1.fastq"));
    let r2p = dir.join(format!("{name}_R2.fastq"));
    fs::write(&r1p, r1).unwrap();
    fs::write(&r2p, r2).unwrap();
    let archive = dir.join(format!("{name}.qz"));
    let cfg = CompressConfig {
        input: vec![r1p, r2p],
        output: archive.clone(),
        working_dir: dir.to_path_buf(),
        threads: 1,
        force: true,
        reference: Some(ReferenceOptions {
            reference: fasta.to_path_buf(),
            reference_index: None,
            reference_fast: false,
            reference_window: 2,
        }),
        advanced: AdvancedOptions { chunk_records, ..AdvancedOptions::default() },
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();
    archive
}

/// Decode a type-2 reference archive → `(r1_text, r2_text)` (writes `<name>_R1/_R2.fastq`).
fn decompress_ref_paired(dir: &Path, archive: &Path, name: &str) -> (String, String) {
    let prefix = dir.join(name);
    let cfg = DecompressConfig {
        input: archive.to_path_buf(),
        output: vec![prefix.clone()],
        working_dir: dir.to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&cfg).unwrap();
    (
        fs::read_to_string(dir.join(format!("{name}_R1.fastq"))).unwrap(),
        fs::read_to_string(dir.join(format!("{name}_R2.fastq"))).unwrap(),
    )
}

#[test]
fn merge_reference_paired_roundtrips_lossless() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let (fasta, r1, r2) = merge_ref_paired_dataset(d, 200);

    // Reference: the whole pair set in one type-2 archive.
    let whole_arc = compress_ref_paired(d, &r1, &r2, "wholep", &fasta, 60);
    let (w1, w2) = decompress_ref_paired(d, &whole_arc, "wholep");
    assert_eq!(w1, r1, "paired reference R1 roundtrip (sanity)");
    assert_eq!(w2, r2, "paired reference R2 roundtrip (sanity)");

    // Shards: split pairs [0..120) and [120..200) in LOCKSTEP across both mates.
    let r1s = split_records(&r1, &[120, 200]);
    let r2s = split_records(&r2, &[120, 200]);
    let s0 = compress_ref_paired(d, &r1s[0], &r2s[0], "rp0", &fasta, 60);
    let s1 = compress_ref_paired(d, &r1s[1], &r2s[1], "rp1", &fasta, 60);
    let merged = d.join("merged_ref_paired.qz");
    qz_lib::compression::merge_reference_archives_to_path(&[s0, s1], &fasta, &merged, true).unwrap();

    let (m1, m2) = decompress_ref_paired(d, &merged, "merged_ref_paired");
    assert_eq!(r1, m1, "merged reference R1 must decode to the original R1 reads");
    assert_eq!(r2, m2, "merged reference R2 must decode to the original R2 reads");

    qz_lib::compression::verify(&verify_cfg(&merged, d, false)).unwrap();
    qz_lib::compression::verify(&verify_cfg(&merged, d, true)).unwrap();
}
