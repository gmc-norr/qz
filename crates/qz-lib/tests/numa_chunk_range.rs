#[allow(unused_imports)]
use qz_lib::cli::{CompressConfig, DecompressConfig};
use std::path::PathBuf;
use std::sync::Mutex;

/// Process-global lock serializing every `compress()` in this file. Required
/// because `compress_single_with_chunk_records_env` mutates process-global env
/// (`QZ_NUMA_FAKE_UNDERSIZE`) around its compress; without serialization a
/// concurrent test's compress would see the leaked var and emit a corrupt size
/// table. The lock is poison-tolerant (a panicking test still holds correctness
/// for the rest).
static COMPRESS_LOCK: Mutex<()> = Mutex::new(());

/// Run `compress(cfg)` under the process-global `COMPRESS_LOCK`.
fn guarded_compress(cfg: &CompressConfig) -> anyhow::Result<()> {
    let _g = COMPRESS_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    qz_lib::compression::compress(cfg)
}

/// Compress N single-end reads into a multi-chunk archive, then assert the layout
/// reports the right archive_type, chunk count, and per-chunk read totals.
#[test]
fn read_chunk_layout_single_end_multichunk() {
    let tmp = tempfile::tempdir().unwrap();
    let fastq = tmp.path().join("in.fastq");
    write_synthetic_fastq(&fastq, 60_000);
    let archive = tmp.path().join("a.qz");
    let mut cfg = CompressConfig {
        input: vec![fastq.clone()],
        output: archive.clone(),
        working_dir: tmp.path().to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 25_000; // force 3 chunks
    guarded_compress(&cfg).unwrap();

    let layout = qz_lib::compression::read_chunk_layout(&archive).unwrap();
    assert_eq!(layout.archive_type, 0);
    assert_eq!(layout.num_chunks, 3);
    assert_eq!(layout.per_chunk_reads.iter().sum::<u64>(), 60_000);
    assert_eq!(layout.num_reads, 60_000);
    assert!(layout.shardable);
    assert_eq!(layout.reference_resident_bytes, 0);
}

/// A cluster archive (archive_type 3) decodes through a bespoke path that is NOT
/// chunk-range shardable. `read_chunk_layout` must report it as a valid,
/// non-shardable layout (NOT bail "unknown archive_type 3"), so the NUMA driver can
/// route it to in-process decode instead of crashing under `--numa N`.
#[test]
fn read_chunk_layout_cluster_is_non_shardable() {
    let tmp = tempfile::tempdir().unwrap();
    let fastq = tmp.path().join("in.fastq");
    write_synthetic_fastq(&fastq, 5_000);
    let archive = tmp.path().join("c.qz");
    let cfg = CompressConfig {
        input: vec![fastq.clone()],
        output: archive.clone(),
        working_dir: tmp.path().to_path_buf(),
        threads: 1,
        cluster: Some(qz_lib::cli::ClusterOptions::default()),
        ..CompressConfig::default()
    };
    guarded_compress(&cfg).unwrap();

    let layout = qz_lib::compression::read_chunk_layout(&archive).unwrap();
    assert_eq!(layout.archive_type, 3, "cluster archive_type");
    assert!(
        !layout.shardable,
        "cluster archives are not chunk-range shardable"
    );
}

// minimal synthetic FASTQ writer (40bp reads, constant qual)
fn write_synthetic_fastq(path: &std::path::Path, n: usize) {
    use std::io::Write;
    let mut w = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for i in 0..n {
        writeln!(w, "@r{i}\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII").unwrap();
    }
}

/// Decoding [0,k) then [k,n) and concatenating the part files must byte-equal a
/// full single-process decode.
#[test]
fn decode_chunk_range_single_end_slice_equivalence() {
    let tmp = tempfile::tempdir().unwrap();
    let fastq = tmp.path().join("in.fastq");
    write_synthetic_fastq(&fastq, 60_000);
    let archive = tmp.path().join("a.qz");
    let mut cfg = CompressConfig {
        input: vec![fastq.clone()],
        output: archive.clone(),
        working_dir: tmp.path().to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 25_000; // 3 chunks
    guarded_compress(&cfg).unwrap();

    let full = tmp.path().join("full.fastq");
    qz_lib::compression::decompress(&DecompressConfig {
        input: archive.clone(),
        output: vec![full.clone()],
        force: true,
        ..decompress_defaults()
    }).unwrap();

    let p0 = tmp.path().join("p0.fastq");
    let p1 = tmp.path().join("p1.fastq");
    qz_lib::compression::decode_chunk_range(&archive, 0, 1, 4, false, 6, std::slice::from_ref(&p0), &[])
        .unwrap();
    qz_lib::compression::decode_chunk_range(&archive, 1, 3, 4, false, 6, std::slice::from_ref(&p1), &[])
        .unwrap();

    let mut joined = std::fs::read(&p0).unwrap();
    joined.extend_from_slice(&std::fs::read(&p1).unwrap());
    assert_eq!(joined, std::fs::read(&full).unwrap(), "sharded != full");
}

fn decompress_defaults() -> DecompressConfig {
    DecompressConfig {
        input: PathBuf::new(),
        output: vec![],
        working_dir: PathBuf::from("."),
        num_threads: 4,
        gzipped: false,
        gzip_level: 6,
        force: false,
    }
}

#[test]
fn decode_chunk_range_paired_slice_equivalence() {
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("r1.fastq");
    let r2 = tmp.path().join("r2.fastq");
    write_synthetic_fastq(&r1, 60_000);
    write_synthetic_fastq(&r2, 60_000);
    let archive = tmp.path().join("p.qz");
    let mut cfg = CompressConfig {
        input: vec![r1.clone(), r2.clone()],
        output: archive.clone(),
        working_dir: tmp.path().to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 25_000; // 3 chunks
    guarded_compress(&cfg).unwrap();

    let full_prefix = tmp.path().join("full");
    qz_lib::compression::decompress(&DecompressConfig {
        input: archive.clone(),
        output: vec![full_prefix.clone()],
        force: true,
        ..decompress_defaults()
    }).unwrap();

    let p0r1 = tmp.path().join("p0_R1"); let p0r2 = tmp.path().join("p0_R2");
    let p1r1 = tmp.path().join("p1_R1"); let p1r2 = tmp.path().join("p1_R2");
    qz_lib::compression::decode_chunk_range(&archive, 0, 1, 4, false, 6, &[p0r1.clone(), p0r2.clone()], &[]).unwrap();
    qz_lib::compression::decode_chunk_range(&archive, 1, 3, 4, false, 6, &[p1r1.clone(), p1r2.clone()], &[]).unwrap();

    for (parts, full) in [
        (vec![&p0r1, &p1r1], with_suffix_test(&full_prefix, "_R1.fastq")),
        (vec![&p0r2, &p1r2], with_suffix_test(&full_prefix, "_R2.fastq")),
    ] {
        let mut joined = Vec::new();
        for p in parts { joined.extend_from_slice(&std::fs::read(p).unwrap()); }
        assert_eq!(joined, std::fs::read(&full).unwrap());
    }
}

/// Paired direct-write equivalence: two workers seek into TWO pre-sized shared outputs
/// (R1/R2) at their per-mate byte offsets and write their chunk-ranges, then the two
/// files byte-equal a full paired decode — the paired analogue of
/// `decode_chunk_range_direct_write_at_offset`. Exercises the new per-mate
/// `decoded_sizes_per_mate` table and `decode_chunk_range_paired`'s direct branch.
#[test]
fn decode_chunk_range_paired_direct_write_equivalence() {
    use qz_lib::compression::DirectWriteRegion;
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("r1.fastq");
    let r2 = tmp.path().join("r2.fastq");
    write_synthetic_fastq(&r1, 60_000);
    write_synthetic_fastq(&r2, 60_000);
    let archive = tmp.path().join("p.qz");
    let mut cfg = CompressConfig {
        input: vec![r1.clone(), r2.clone()],
        output: archive.clone(),
        working_dir: tmp.path().to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 25_000; // 3 chunks
    guarded_compress(&cfg).unwrap();

    // Reference: a full paired decode.
    let full_prefix = tmp.path().join("full");
    qz_lib::compression::decompress(&DecompressConfig {
        input: archive.clone(),
        output: vec![full_prefix.clone()],
        force: true,
        ..decompress_defaults()
    })
    .unwrap();
    let full_r1 = with_suffix_test(&full_prefix, "_R1.fastq");
    let full_r2 = with_suffix_test(&full_prefix, "_R2.fastq");

    // The per-mate decoded-size table must be present with exactly two mates, and each
    // mate's totals must match the full per-mate output size + the combined per-chunk sum.
    let layout = qz_lib::compression::read_chunk_layout(&archive).unwrap();
    assert_eq!(layout.archive_type, 1);
    assert_eq!(layout.decoded_sizes_per_mate.len(), 2, "paired must expose 2 per-mate size vectors");
    let nchunks = layout.num_chunks as usize;
    for m in 0..2 {
        assert_eq!(layout.decoded_sizes_per_mate[m].len(), nchunks);
    }
    // combined per-chunk == Σ mates per chunk.
    for c in 0..nchunks {
        assert_eq!(
            layout.per_chunk_decoded_bytes[c],
            layout.decoded_sizes_per_mate[0][c] + layout.decoded_sizes_per_mate[1][c]
        );
    }
    let r1_total: u64 = layout.decoded_sizes_per_mate[0].iter().sum();
    let r2_total: u64 = layout.decoded_sizes_per_mate[1].iter().sum();
    assert_eq!(r1_total, std::fs::metadata(&full_r1).unwrap().len());
    assert_eq!(r2_total, std::fs::metadata(&full_r2).unwrap().len());

    // Pre-size the two shared outputs.
    let out_r1 = tmp.path().join("direct_R1.fastq");
    let out_r2 = tmp.path().join("direct_R2.fastq");
    for (p, total) in [(&out_r1, r1_total), (&out_r2, r2_total)] {
        std::fs::OpenOptions::new().write(true).create_new(true).open(p).unwrap().set_len(total).unwrap();
    }

    // Two workers: [0,1) and [1,3). Each computes BOTH mates' regions from the per-mate
    // prefix sums and writes straight into the shared outputs.
    let region = |m: usize, a: usize, b: usize| {
        let off: u64 = layout.decoded_sizes_per_mate[m][..a].iter().sum();
        let len: u64 = layout.decoded_sizes_per_mate[m][a..b].iter().sum();
        DirectWriteRegion { offset: off, len }
    };
    let outs = [out_r1.clone(), out_r2.clone()];
    for (a, b) in [(0usize, 1usize), (1, 3)] {
        qz_lib::compression::decode_chunk_range(
            &archive, a as u32, b as u32, 4, false, 0,
            &outs,
            &[region(0, a, b), region(1, a, b)],
        )
        .unwrap();
    }

    assert_eq!(std::fs::read(&out_r1).unwrap(), std::fs::read(&full_r1).unwrap(), "R1 direct != full");
    assert_eq!(std::fs::read(&out_r2).unwrap(), std::fs::read(&full_r2).unwrap(), "R2 direct != full");
}

/// A paired direct-write region whose len disagrees with the per-mate table is a clean
/// `DirectWriteIntegrityError` (so `auto` can fall back), not silent corruption.
#[test]
fn decode_chunk_range_paired_direct_rejects_bad_region() {
    use qz_lib::compression::{DirectWriteRegion, DirectWriteIntegrityError};
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("r1.fastq");
    let r2 = tmp.path().join("r2.fastq");
    write_synthetic_fastq(&r1, 40_000);
    write_synthetic_fastq(&r2, 40_000);
    let archive = tmp.path().join("p.qz");
    let mut cfg = CompressConfig {
        input: vec![r1.clone(), r2.clone()],
        output: archive.clone(),
        working_dir: tmp.path().to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 20_000; // 2 chunks
    guarded_compress(&cfg).unwrap();

    let layout = qz_lib::compression::read_chunk_layout(&archive).unwrap();
    let r1_total: u64 = layout.decoded_sizes_per_mate[0].iter().sum();
    let r2_total: u64 = layout.decoded_sizes_per_mate[1].iter().sum();
    let out_r1 = tmp.path().join("o_R1.fastq");
    let out_r2 = tmp.path().join("o_R2.fastq");
    for (p, total) in [(&out_r1, r1_total), (&out_r2, r2_total)] {
        std::fs::OpenOptions::new().write(true).create_new(true).open(p).unwrap().set_len(total).unwrap();
    }
    let outs = [out_r1, out_r2];
    let good0 = DirectWriteRegion { offset: 0, len: layout.decoded_sizes_per_mate[0][0] };

    // Wrong R2 len → integrity error.
    let bad1 = DirectWriteRegion { offset: 0, len: layout.decoded_sizes_per_mate[1][0] + 1 };
    let err = qz_lib::compression::decode_chunk_range(&archive, 0, 1, 1, false, 0, &outs, &[good0, bad1]).unwrap_err();
    assert!(err.downcast_ref::<DirectWriteIntegrityError>().is_some(), "bad R2 region must be integrity: {err}");

    // Wrong region COUNT (1 region for a 2-mate archive) → integrity error.
    let err2 = qz_lib::compression::decode_chunk_range(&archive, 0, 1, 1, false, 0, &outs, &[good0]).unwrap_err();
    assert!(err2.downcast_ref::<DirectWriteIntegrityError>().is_some(), "wrong region count must be integrity: {err2}");
}

fn with_suffix_test(prefix: &std::path::Path, suffix: &str) -> std::path::PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    std::path::PathBuf::from(s)
}

fn compress_single_with_chunk_records(fq: &std::path::Path, arc: &std::path::Path, n: usize) {
    let mut cfg = CompressConfig {
        input: vec![fq.to_path_buf()],
        output: arc.to_path_buf(),
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = n;
    guarded_compress(&cfg).unwrap();
}

// Only used by `decode_chunk_range_bounded_overrun_is_integrity_error`, which is itself
// `#[cfg(debug_assertions)]`-gated (it drives the debug-only `QZ_NUMA_FAKE_UNDERSIZE`
// hook). Gate the helper to match so it isn't dead code under `--release`.
#[cfg(debug_assertions)]
fn compress_single_with_chunk_records_env(
    fq: &std::path::Path,
    arc: &std::path::Path,
    n: usize,
    envs: &[(&str, &str)],
) {
    // Hold COMPRESS_LOCK across the whole set-var/compress/remove-var window so no
    // concurrent test's compress can observe the leaked QZ_NUMA_FAKE_UNDERSIZE.
    let _g = COMPRESS_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    for (k, v) in envs {
        unsafe { std::env::set_var(k, v) };
    }
    let mut cfg = CompressConfig {
        input: vec![fq.to_path_buf()],
        output: arc.to_path_buf(),
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = n;
    let res = qz_lib::compression::compress(&cfg);
    for (k, _) in envs {
        unsafe { std::env::remove_var(k) };
    }
    res.unwrap();
}

fn decompress_to(arc: &std::path::Path, out: &std::path::Path) {
    qz_lib::compression::decompress(&DecompressConfig {
        input: arc.to_path_buf(),
        output: vec![out.to_path_buf()],
        force: true,
        ..decompress_defaults()
    })
    .unwrap();
}

#[test]
fn read_chunk_layout_exposes_per_chunk_decoded_bytes() {
    let d = tempfile::tempdir().unwrap();
    let fq = d.path().join("in.fastq");
    let mut s = String::new();
    for i in 0..5 {
        s.push_str(&format!("@read{i}\nACGTACGT\n+\nIIIIIIII\n"));
    }
    std::fs::write(&fq, &s).unwrap();
    let arc = d.path().join("out.qz");
    compress_single_with_chunk_records(&fq, &arc, 2);

    let layout = qz_lib::compression::read_chunk_layout(&arc).unwrap();
    assert_eq!(layout.archive_type, 0);
    assert_eq!(
        layout.per_chunk_decoded_bytes.len(),
        layout.num_chunks as usize
    );
    assert_eq!(layout.per_chunk_decoded_bytes.iter().sum::<u64>(), 5 * 27);

    let out = d.path().join("rt.fastq");
    decompress_to(&arc, &out);
    assert_eq!(std::fs::metadata(&out).unwrap().len(), 5 * 27); // total == actual output bytes
}

#[test]
fn decode_chunk_range_paired_rejects_forged_encoding_type() {
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("r1.fastq");
    let r2 = tmp.path().join("r2.fastq");
    write_synthetic_fastq(&r1, 50_000);
    write_synthetic_fastq(&r2, 50_000);
    let archive = tmp.path().join("p.qz");
    let mut cfg = CompressConfig {
        input: vec![r1.clone(), r2.clone()],
        output: archive.clone(),
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 25_000; // multi-chunk
    guarded_compress(&cfg).unwrap();

    // Sanity: the unmodified archive decodes fine via the range primitive.
    let ok1 = tmp.path().join("ok_R1"); let ok2 = tmp.path().join("ok_R2");
    qz_lib::compression::decode_chunk_range(&archive, 0, 1, 4, false, 6, &[ok1, ok2], &[]).unwrap();

    // Forge encoding_type: body byte 0 lives just past the 8-byte v5 prefix → file offset 8.
    // (magic[2] + version[1] + archive_type[1] + header_size[4] = 8; body[0] = encoding_type.)
    let mut bytes = std::fs::read(&archive).unwrap();
    assert_eq!(bytes[8], 0, "paired archive should have encoding_type 0 at byte 8");
    bytes[8] = 4; // RawWithHints
    let forged = tmp.path().join("forged.qz");
    std::fs::write(&forged, &bytes).unwrap();

    // The range primitive MUST reject the forged archive (matching --numa off).
    let f1 = tmp.path().join("f_R1"); let f2 = tmp.path().join("f_R2");
    let res = qz_lib::compression::decode_chunk_range(&forged, 0, 1, 4, false, 6, &[f1, f2], &[]);
    assert!(res.is_err(), "forged encoding_type must be rejected, not silently decoded");
}

/// Build a deterministic low-repeat ACGT reference sequence (matches the reference
/// integration test generator so reads seed reliably).
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

/// Build a paired reference fixture (lifted from
/// `reference_integration.rs::make_sized_reference_dataset`): a short reference FASTA
/// plus `reads` paired reads that mostly map against it. R1 is a forward 120bp slice,
/// R2 the reverse-complement of a downstream slice; ~2% are off-reference junk
/// (literal fallback) and a sparse subset carries a single substitution (edit path).
/// Returns `(r1_path, r2_path, fasta_path)`.
fn make_reference_inputs(
    dir: &std::path::Path,
    reads: usize,
) -> (PathBuf, PathBuf, PathBuf) {
    let ref_len = 50_000usize;
    let refseq = make_seq(ref_len, 7);

    let rf = dir.join("ref.fa");
    {
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        std::fs::write(&rf, &s).unwrap();
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

    let rlen = 120usize; // fixed read length (constant lengths → tidy framing)
    let max_start = ref_len - rlen - 200; // leave room for the R2 downstream slice
    let mut rng = 0x1234_5678_9abc_def0u64;
    let mut next = || {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (rng >> 33) as usize
    };

    let mut s1 = String::with_capacity(reads * 270);
    let mut s2 = String::with_capacity(reads * 270);

    for i in 0..reads {
        // Deterministic ~2% off-reference junk → literal fallback.
        let unmappable = i % 50 == 0;

        let (r1bytes, r2bytes) = if unmappable {
            (
                make_seq(rlen, 0x9000_0000 + i as u64),
                make_seq(rlen, 0x7000_0000 + i as u64),
            )
        } else {
            let st1 = next() % (max_start + 1);
            let mut r1 = refseq[st1..st1 + rlen].to_vec();
            let st2 = st1 + 150;
            let r2_fwd = refseq[st2..st2 + rlen].to_vec();
            let mut r2 = revcomp(&r2_fwd);
            // Plant a sparse substitution on a subset so the edit path is covered.
            if i % 17 == 0 {
                let p = next() % rlen;
                r1[p] = match r1[p] {
                    b'A' => b'C',
                    _ => b'A',
                };
            }
            if i % 23 == 0 {
                let p = next() % rlen;
                r2[p] = match r2[p] {
                    b'G' => b'T',
                    _ => b'G',
                };
            }
            (r1, r2)
        };

        let q: &str = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\
                       IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"; // 120 'I's
        debug_assert_eq!(q.len(), 120);
        s1.push_str("@frag_");
        s1.push_str(&i.to_string());
        s1.push_str("/1\n");
        s1.push_str(std::str::from_utf8(&r1bytes).unwrap());
        s1.push_str("\n+\n");
        s1.push_str(q);
        s1.push('\n');

        s2.push_str("@frag_");
        s2.push_str(&i.to_string());
        s2.push_str("/2\n");
        s2.push_str(std::str::from_utf8(&r2bytes).unwrap());
        s2.push_str("\n+\n");
        s2.push_str(q);
        s2.push('\n');
    }

    let r1p = dir.join("r1.fastq");
    let r2p = dir.join("r2.fastq");
    std::fs::write(&r1p, &s1).unwrap();
    std::fs::write(&r2p, &s2).unwrap();
    (r1p, r2p, rf)
}

#[test]
fn decode_chunk_range_reference_slice_equivalence() {
    let tmp = tempfile::tempdir().unwrap();
    let (r1, r2, fa) = make_reference_inputs(tmp.path(), 30_000);
    let archive = tmp.path().join("ref.qz");
    let mut cfg = CompressConfig {
        input: vec![r1.clone(), r2.clone()],
        output: archive.clone(),
        reference: Some(qz_lib::cli::ReferenceOptions {
            reference: fa.clone(),
            reference_index: None,
            reference_fast: false,
            reference_window: 4,
        }),
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 10_000; // 3 chunks
    guarded_compress(&cfg).unwrap();

    let full_prefix = tmp.path().join("full");
    qz_lib::compression::decompress(&DecompressConfig {
        input: archive.clone(),
        output: vec![full_prefix.clone()],
        force: true,
        ..decompress_defaults()
    }).unwrap();

    let p0 = [tmp.path().join("p0_R1"), tmp.path().join("p0_R2")];
    let p1 = [tmp.path().join("p1_R1"), tmp.path().join("p1_R2")];
    qz_lib::compression::decode_chunk_range(&archive, 0, 1, 4, false, 6, &p0, &[]).unwrap();
    qz_lib::compression::decode_chunk_range(&archive, 1, 3, 4, false, 6, &p1, &[]).unwrap();

    for (parts, full) in [
        ([&p0[0], &p1[0]], with_suffix_test(&full_prefix, "_R1.fastq")),
        ([&p0[1], &p1[1]], with_suffix_test(&full_prefix, "_R2.fastq")),
    ] {
        let mut joined = Vec::new();
        for p in parts { joined.extend_from_slice(&std::fs::read(p).unwrap()); }
        assert_eq!(joined, std::fs::read(&full).unwrap());
    }
}

/// Paired-reference (type 2) direct-write: two workers seek into pre-sized R1/R2 outputs
/// at per-mate offsets and write their chunk-ranges → byte-equal to a full reference
/// decode. Exercises the reference encoder's `ChunkDecodedSizes(n_mates=2)` global +
/// `decode_reference_range`'s direct branch.
#[test]
fn decode_chunk_range_reference_direct_write_equivalence() {
    use qz_lib::compression::DirectWriteRegion;
    let tmp = tempfile::tempdir().unwrap();
    let (r1, r2, fa) = make_reference_inputs(tmp.path(), 30_000);
    let archive = tmp.path().join("ref.qz");
    let mut cfg = CompressConfig {
        input: vec![r1, r2],
        output: archive.clone(),
        reference: Some(qz_lib::cli::ReferenceOptions {
            reference: fa,
            reference_index: None,
            reference_fast: false,
            reference_window: 4,
        }),
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 10_000; // 3 chunks
    guarded_compress(&cfg).unwrap();

    let full_prefix = tmp.path().join("full");
    qz_lib::compression::decompress(&DecompressConfig {
        input: archive.clone(),
        output: vec![full_prefix.clone()],
        force: true,
        ..decompress_defaults()
    })
    .unwrap();
    let full_r1 = with_suffix_test(&full_prefix, "_R1.fastq");
    let full_r2 = with_suffix_test(&full_prefix, "_R2.fastq");

    let layout = qz_lib::compression::read_chunk_layout(&archive).unwrap();
    assert_eq!(layout.archive_type, 2);
    assert_eq!(layout.decoded_sizes_per_mate.len(), 2, "paired reference must expose 2 per-mate size vectors");
    let r1_total: u64 = layout.decoded_sizes_per_mate[0].iter().sum();
    let r2_total: u64 = layout.decoded_sizes_per_mate[1].iter().sum();
    assert_eq!(r1_total, std::fs::metadata(&full_r1).unwrap().len());
    assert_eq!(r2_total, std::fs::metadata(&full_r2).unwrap().len());

    let out_r1 = tmp.path().join("direct_R1.fastq");
    let out_r2 = tmp.path().join("direct_R2.fastq");
    for (p, total) in [(&out_r1, r1_total), (&out_r2, r2_total)] {
        std::fs::OpenOptions::new().write(true).create_new(true).open(p).unwrap().set_len(total).unwrap();
    }
    let region = |m: usize, a: usize, b: usize| {
        let off: u64 = layout.decoded_sizes_per_mate[m][..a].iter().sum();
        let len: u64 = layout.decoded_sizes_per_mate[m][a..b].iter().sum();
        DirectWriteRegion { offset: off, len }
    };
    let outs = [out_r1.clone(), out_r2.clone()];
    for (a, b) in [(0usize, 1usize), (1, 3)] {
        qz_lib::compression::decode_chunk_range(
            &archive, a as u32, b as u32, 4, false, 0,
            &outs,
            &[region(0, a, b), region(1, a, b)],
        )
        .unwrap();
    }
    assert_eq!(std::fs::read(&out_r1).unwrap(), std::fs::read(&full_r1).unwrap(), "R1 ref direct != full");
    assert_eq!(std::fs::read(&out_r2).unwrap(), std::fs::read(&full_r2).unwrap(), "R2 ref direct != full");
}

/// Single-end reference (type 4) direct-write: one worker per chunk-range writes into the
/// pre-sized output at the table offset → byte-equal to a full decode.
#[test]
fn decode_chunk_range_single_end_reference_direct_write_equivalence() {
    use qz_lib::compression::DirectWriteRegion;
    let tmp = tempfile::tempdir().unwrap();
    let (r1, _r2, fa) = make_reference_inputs(tmp.path(), 30_000);
    let archive = tmp.path().join("refse.qz");
    let mut cfg = CompressConfig {
        input: vec![r1], // single input → single-end reference (type 4)
        output: archive.clone(),
        reference: Some(qz_lib::cli::ReferenceOptions {
            reference: fa,
            reference_index: None,
            reference_fast: false,
            reference_window: 4,
        }),
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 10_000; // 3 chunks
    guarded_compress(&cfg).unwrap();

    let full = tmp.path().join("full.fastq");
    decompress_to(&archive, &full);

    let layout = qz_lib::compression::read_chunk_layout(&archive).unwrap();
    assert_eq!(layout.archive_type, 4);
    assert_eq!(layout.decoded_sizes_per_mate.len(), 1, "single-end reference table is 1-mate");
    let total: u64 = layout.decoded_sizes_per_mate[0].iter().sum();
    assert_eq!(total, std::fs::metadata(&full).unwrap().len());

    let out = tmp.path().join("direct.fastq");
    std::fs::OpenOptions::new().write(true).create_new(true).open(&out).unwrap().set_len(total).unwrap();
    let region = |a: usize, b: usize| {
        let off: u64 = layout.decoded_sizes_per_mate[0][..a].iter().sum();
        let len: u64 = layout.decoded_sizes_per_mate[0][a..b].iter().sum();
        DirectWriteRegion { offset: off, len }
    };
    for (a, b) in [(0usize, 1usize), (1, 3)] {
        qz_lib::compression::decode_chunk_range(
            &archive, a as u32, b as u32, 4, false, 0,
            std::slice::from_ref(&out),
            &[region(a, b)],
        )
        .unwrap();
    }
    assert_eq!(std::fs::read(&out).unwrap(), std::fs::read(&full).unwrap(), "se-ref direct != full");
}

#[test]
fn decode_chunk_range_direct_write_at_offset() {
    use qz_lib::compression::DirectWriteRegion;
    let d = tempfile::tempdir().unwrap();
    let fq = d.path().join("in.fastq");
    let mut s = String::new();
    for i in 0..6 { s.push_str(&format!("@read{i}\nACGTACGT\n+\nIIIIIIII\n")); }
    std::fs::write(&fq, &s).unwrap();
    let arc = d.path().join("out.qz");
    compress_single_with_chunk_records(&fq, &arc, 2); // 3 chunks
    let layout = qz_lib::compression::read_chunk_layout(&arc).unwrap();
    let total: u64 = layout.per_chunk_decoded_bytes.iter().sum();
    let out = d.path().join("direct.fastq");
    { std::fs::OpenOptions::new().write(true).create_new(true).open(&out).unwrap().set_len(total).unwrap(); }
    let base_b: u64 = layout.per_chunk_decoded_bytes[..2].iter().sum();
    let len_a: u64 = layout.per_chunk_decoded_bytes[..2].iter().sum();
    let len_b: u64 = layout.per_chunk_decoded_bytes[2..3].iter().sum();
    qz_lib::compression::decode_chunk_range(&arc, 0, 2, 2, false, 0, std::slice::from_ref(&out), &[DirectWriteRegion { offset: 0, len: len_a }]).unwrap();
    qz_lib::compression::decode_chunk_range(&arc, 2, 3, 2, false, 0, std::slice::from_ref(&out), &[DirectWriteRegion { offset: base_b, len: len_b }]).unwrap();
    let r = d.path().join("ref.fastq"); decompress_to(&arc, &r);
    assert_eq!(std::fs::read(&out).unwrap(), std::fs::read(&r).unwrap());
}

#[test]
fn decode_chunk_range_rejects_wrong_base_and_len() {
    use qz_lib::compression::DirectWriteRegion;
    let d = tempfile::tempdir().unwrap();
    let fq = d.path().join("in.fastq");
    let mut s = String::new(); for i in 0..4 { s.push_str(&format!("@r{i}\nACGT\n+\nIIII\n")); }
    std::fs::write(&fq, &s).unwrap();
    let arc = d.path().join("out.qz");
    compress_single_with_chunk_records(&fq, &arc, 2); // 2 chunks
    let layout = qz_lib::compression::read_chunk_layout(&arc).unwrap();
    let total: u64 = layout.per_chunk_decoded_bytes.iter().sum();
    let out = d.path().join("o.fastq");
    { std::fs::OpenOptions::new().write(true).create_new(true).open(&out).unwrap().set_len(total).unwrap(); }
    let correct_len = layout.per_chunk_decoded_bytes[0];
    assert!(qz_lib::compression::decode_chunk_range(&arc, 0, 1, 1, false, 0, std::slice::from_ref(&out), &[DirectWriteRegion { offset: 7, len: correct_len }]).is_err());           // wrong base
    assert!(qz_lib::compression::decode_chunk_range(&arc, 0, 1, 1, false, 0, std::slice::from_ref(&out), &[DirectWriteRegion { offset: 0, len: correct_len + 1 }]).is_err());        // wrong len
    assert!(qz_lib::compression::decode_chunk_range(&arc, 0, 1, 1, true, 6, std::slice::from_ref(&out), &[DirectWriteRegion { offset: 0, len: correct_len }]).is_err());            // gzip + region
}

// Debug-only: relies on the `QZ_NUMA_FAKE_UNDERSIZE` writer hook, which is itself
// `#[cfg(debug_assertions)]`-gated (the fault injector must not ship in release
// binaries). Under `--release` the hook is inert, the archive is correct, and the
// `unwrap_err()` below would panic — so the test is gated to debug to match its hook.
// The production integrity check it exercises (`DirectWriteIntegrityError`) is compiled
// in all profiles; only the fault-injection trigger is debug-only.
#[cfg(debug_assertions)]
#[test]
fn decode_chunk_range_bounded_overrun_is_integrity_error() {
    use qz_lib::compression::{DirectWriteRegion, DirectWriteIntegrityError};
    let d = tempfile::tempdir().unwrap();
    let fq = d.path().join("in.fastq");
    let mut s = String::new(); for i in 0..6 { s.push_str(&format!("@r{i}\nACGT\n+\nIIII\n")); }
    std::fs::write(&fq, &s).unwrap();
    let arc = d.path().join("out.qz");
    compress_single_with_chunk_records_env(&fq, &arc, 2, &[("QZ_NUMA_FAKE_UNDERSIZE", "1")]); // last chunk size short by 1
    let layout = qz_lib::compression::read_chunk_layout(&arc).unwrap();
    let last = layout.num_chunks - 1;
    let base: u64 = layout.per_chunk_decoded_bytes[..last as usize].iter().sum();
    let len = layout.per_chunk_decoded_bytes[last as usize]; // undersized
    let total: u64 = layout.per_chunk_decoded_bytes.iter().sum();
    let out = d.path().join("o.fastq");
    { std::fs::OpenOptions::new().write(true).create_new(true).open(&out).unwrap().set_len(total).unwrap(); }
    let err = qz_lib::compression::decode_chunk_range(&arc, last, last + 1, 1, false, 0, std::slice::from_ref(&out), &[DirectWriteRegion { offset: base, len }]).unwrap_err();
    assert!(err.downcast_ref::<DirectWriteIntegrityError>().is_some(), "overrun must classify as integrity, got: {err}");
}

/// Compress a FASTA input (no quality lines) into a multi-chunk archive.
fn compress_fasta_with_chunk_records(fa: &std::path::Path, arc: &std::path::Path, n: usize) {
    let mut cfg = CompressConfig {
        input: vec![fa.to_path_buf()],
        output: arc.to_path_buf(),
        fasta: true,
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = n;
    guarded_compress(&cfg).unwrap();
}

/// Larger multichunk direct-write: ~200 reads / ~7 chunks decoded by TWO workers
/// covering lower and upper halves, then joined → byte-identical to full decode.
#[test]
fn direct_write_equivalence_multichunk() {
    use qz_lib::compression::DirectWriteRegion;
    let d = tempfile::tempdir().unwrap();
    let fq = d.path().join("in.fastq");

    // 200 reads, 80bp, so we get plenty of chunks at chunk_records=28 → ceil(200/28)=8 chunks
    let mut s = String::new();
    for i in 0..200 {
        let seq: String = "ACGTACGTACGTACGT".repeat(5); // 80bp
        let qual: String = "I".repeat(80);
        s.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
    }
    std::fs::write(&fq, &s).unwrap();
    let arc = d.path().join("out.qz");
    compress_single_with_chunk_records(&fq, &arc, 28); // ~8 chunks

    let layout = qz_lib::compression::read_chunk_layout(&arc).unwrap();
    let num_chunks = layout.num_chunks as usize;
    assert!(num_chunks >= 7, "expected >= 7 chunks, got {num_chunks}");
    assert!(!layout.per_chunk_decoded_bytes.is_empty(), "decoded sizes table missing");

    // Split at midpoint
    let mid = num_chunks / 2;
    let lower_len: u64 = layout.per_chunk_decoded_bytes[..mid].iter().sum();
    let upper_offset: u64 = lower_len;
    let upper_len: u64 = layout.per_chunk_decoded_bytes[mid..].iter().sum();
    let total = lower_len + upper_len;

    let out = d.path().join("direct.fastq");
    {
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&out)
            .unwrap()
            .set_len(total)
            .unwrap();
    }

    // Worker 1: lower half [0, mid)
    qz_lib::compression::decode_chunk_range(
        &arc, 0, mid as u32, 4, false, 0,
        std::slice::from_ref(&out),
        &[DirectWriteRegion { offset: 0, len: lower_len }],
    ).unwrap();

    // Worker 2: upper half [mid, num_chunks)
    qz_lib::compression::decode_chunk_range(
        &arc, mid as u32, num_chunks as u32, 4, false, 0,
        std::slice::from_ref(&out),
        &[DirectWriteRegion { offset: upper_offset, len: upper_len }],
    ).unwrap();

    let ref_out = d.path().join("ref.fastq");
    decompress_to(&arc, &ref_out);
    assert_eq!(
        std::fs::read(&out).unwrap(),
        std::fs::read(&ref_out).unwrap(),
        "two-worker direct-write != full decode"
    );
}

/// FASTA direct-write end-to-end: compress FASTA records into ~6 chunks, read the
/// layout, pre-size the output, decode the full range in one direct-write call, and
/// assert byte-identical to a full decompress. Validates the `+2` FASTA branch of
/// the decoded-size formula through real compress→direct-decode.
#[test]
fn direct_write_equivalence_fasta() {
    use qz_lib::compression::DirectWriteRegion;
    let d = tempfile::tempdir().unwrap();
    let fa = d.path().join("in.fasta");

    // 50 FASTA records, 60bp each → chunk_records=8 → ceil(50/8)=7 chunks
    let mut s = String::new();
    for i in 0..50 {
        let seq: String = "ACGTACGTACGT".chars().cycle().take(60).collect(); // 60bp
        s.push_str(&format!(">seq{i}\n{seq}\n"));
    }
    std::fs::write(&fa, &s).unwrap();
    let arc = d.path().join("out.qz");
    compress_fasta_with_chunk_records(&fa, &arc, 8); // ~7 chunks

    let layout = qz_lib::compression::read_chunk_layout(&arc).unwrap();
    let num_chunks = layout.num_chunks as usize;
    assert!(num_chunks >= 6, "expected >= 6 chunks, got {num_chunks}");
    assert!(!layout.per_chunk_decoded_bytes.is_empty(), "decoded sizes table missing");

    let total: u64 = layout.per_chunk_decoded_bytes.iter().sum();

    let out = d.path().join("direct.fasta");
    {
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&out)
            .unwrap()
            .set_len(total)
            .unwrap();
    }

    // Decode the full archive in a single direct-write call
    qz_lib::compression::decode_chunk_range(
        &arc, 0, num_chunks as u32, 2, false, 0,
        std::slice::from_ref(&out),
        &[DirectWriteRegion { offset: 0, len: total }],
    ).unwrap();

    let ref_out = d.path().join("ref.fasta");
    decompress_to(&arc, &ref_out);
    assert_eq!(
        std::fs::read(&out).unwrap(),
        std::fs::read(&ref_out).unwrap(),
        "FASTA direct-write != full decode"
    );
}

/// A single-end reference archive (archive_type 4) must read back as a shardable
/// reference layout: type 4, >=3 chunks, resident globals > 0, per-chunk READ table
/// populated, and Σ per-chunk reads == num_reads. NOTE: reference archives write NO
/// `ChunkDecodedSizes` global, so `per_chunk_decoded_bytes` is EMPTY for type 4 — the
/// NUMA driver falls back to part-file assembly for it (asserted in Task 6.3).
#[test]
fn read_chunk_layout_single_end_reference_multichunk() {
    let tmp = tempfile::tempdir().unwrap();
    // make_reference_inputs builds (r1, r2, fa); single-end reference uses r1 alone.
    let (r1, _r2, fa) = make_reference_inputs(tmp.path(), 30_000);
    let archive = tmp.path().join("refse.qz");
    let mut cfg = CompressConfig {
        input: vec![r1.clone()], // ONE input → archive_type 4
        output: archive.clone(),
        reference: Some(qz_lib::cli::ReferenceOptions {
            reference: fa.clone(),
            reference_index: None,
            reference_fast: false,
            reference_window: 4,
        }),
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 10_000; // 3 chunks
    guarded_compress(&cfg).unwrap();

    let layout = qz_lib::compression::read_chunk_layout(&archive).unwrap();
    assert_eq!(layout.archive_type, 4, "single-end reference archive_type");
    assert_eq!(layout.num_chunks, 3);
    assert_eq!(layout.per_chunk_reads.iter().sum::<u64>(), layout.num_reads);
    assert_eq!(layout.num_reads, 30_000);
    assert!(layout.shardable, "type-4 reference must be shardable");
    assert!(layout.reference_resident_bytes > 0, "resident globals must be folded in");
    // Single-end reference now carries a 1-mate ChunkDecodedSizes global (so NUMA decode
    // can direct-write into a pre-sized output). The combined per-chunk table and the
    // per-mate table both have one entry per chunk; for n_mates=1 they are equal.
    assert_eq!(
        layout.per_chunk_decoded_bytes.len(),
        layout.num_chunks as usize,
        "type-4 reference must carry a ChunkDecodedSizes table"
    );
    assert_eq!(layout.decoded_sizes_per_mate.len(), 1, "single-end reference table is 1-mate");
    assert_eq!(layout.decoded_sizes_per_mate[0], layout.per_chunk_decoded_bytes);
    assert_eq!(
        layout.per_chunk_decoded_bytes.iter().sum::<u64>(),
        // Each read decodes to `@id\nseq\n+\nqual\n`; the table total must equal the
        // actual decoded FASTQ size.
        {
            let out = tmp.path().join("rt_probe.fastq");
            decompress_to(&archive, &out);
            std::fs::metadata(&out).unwrap().len()
        },
        "type-4 table total must equal decoded output size"
    );

    // is_reference_archive stays type-2-only: type 4 is a SINGLE-output archive and
    // must remain decodable through the qz-python single-output API, so it must NOT be
    // classified as a (two-output) reference archive here.
    assert!(!qz_lib::compression::is_reference_archive(&archive).unwrap());
}

/// Decoding [0,1) then [1,3) of a single-end reference (type 4) archive and
/// concatenating the part files must byte-equal a full single-process decode.
#[test]
fn decode_chunk_range_single_end_reference_slice_equivalence() {
    let tmp = tempfile::tempdir().unwrap();
    let (r1, _r2, fa) = make_reference_inputs(tmp.path(), 30_000);
    let archive = tmp.path().join("refse.qz");
    let mut cfg = CompressConfig {
        input: vec![r1.clone()], // single-end reference → archive_type 4
        output: archive.clone(),
        reference: Some(qz_lib::cli::ReferenceOptions {
            reference: fa.clone(),
            reference_index: None,
            reference_fast: false,
            reference_window: 4,
        }),
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 10_000; // 3 chunks
    guarded_compress(&cfg).unwrap();

    // Full single-output decode (single-end reference writes ONE file at the
    // verbatim output[0] path — no suffix, like single-end default).
    let full = tmp.path().join("full.fastq");
    qz_lib::compression::decompress(&DecompressConfig {
        input: archive.clone(),
        output: vec![full.clone()],
        force: true,
        ..decompress_defaults()
    }).unwrap();

    // Range decode to part files (single out_part each) and concatenate.
    let p0 = tmp.path().join("p0.fastq");
    let p1 = tmp.path().join("p1.fastq");
    qz_lib::compression::decode_chunk_range(&archive, 0, 1, 4, false, 6, std::slice::from_ref(&p0), &[]).unwrap();
    qz_lib::compression::decode_chunk_range(&archive, 1, 3, 4, false, 6, std::slice::from_ref(&p1), &[]).unwrap();

    let mut joined = std::fs::read(&p0).unwrap();
    joined.extend_from_slice(&std::fs::read(&p1).unwrap());
    assert_eq!(joined, std::fs::read(&full).unwrap(), "type-4 sharded != full");
}
