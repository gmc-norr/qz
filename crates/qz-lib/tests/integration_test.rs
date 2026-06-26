use qz_lib::cli::{
    AdvancedOptions, CompressConfig, DecompressConfig, QualityCompressor, QualityMode,
};
use std::fs;
use tempfile::TempDir;

/// Standard multi-read test data for roundtrip tests.
const MULTI_READ_DATA: &str = "\
@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHH\n\
@read3\nAAAACCCCGGGGTTTT\n+\nBBBBBBBBBBBBBBBB\n";

/// Roundtrip helper: compress with given config overrides, decompress, assert exact match.
/// Returns (original, decompressed) content for tests that need custom assertions.
fn roundtrip(test_data: &str, config_fn: impl FnOnce(&mut CompressConfig)) -> (String, String) {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let mut compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    config_fn(&mut compress_args);

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    (original, decompressed)
}

/// Roundtrip helper that asserts exact line-by-line match.
fn roundtrip_exact(test_data: &str, config_fn: impl FnOnce(&mut CompressConfig)) {
    let (original, decompressed) = roundtrip(test_data, config_fn);
    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (o, d)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(o, d, "mismatch at line {}", i);
    }
}

/// Helper to create default DecompressConfig
fn decompress_args(
    input: std::path::PathBuf,
    output: std::path::PathBuf,
    working_dir: std::path::PathBuf,
) -> DecompressConfig {
    DecompressConfig {
        input,
        output: vec![output],
        working_dir,
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: false,
    }
}

#[test]
fn test_decompress_refuses_to_overwrite_without_force() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("in.fastq");
    fs::write(&input_fastq, "@r1\nACGT\n+\nIIII\n").unwrap();
    let archive = temp_path.join("a.qz");
    let cfg = CompressConfig {
        input: vec![input_fastq],
        output: archive.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();

    let out = temp_path.join("out.fastq");
    fs::write(&out, "existing").unwrap();

    // Without force: refuse and leave the existing file untouched.
    let mut d = decompress_args(archive.clone(), out.clone(), temp_path.clone());
    assert!(qz_lib::compression::decompress(&d).is_err());
    assert_eq!(fs::read_to_string(&out).unwrap(), "existing");

    // With force: overwrite.
    d.force = true;
    qz_lib::compression::decompress(&d).unwrap();
    assert!(fs::read_to_string(&out).unwrap().contains("ACGT"));
}

#[test]
fn test_fast_verify_rejects_tampered_header_size() {
    // The buffered/fast-verify parser used to ignore the stored header_size while
    // the deep/streaming path trusts it as the data offset, so a tampered value
    // could pass fast verify yet fail to decompress. Fast verify must now reject it.
    let temp = TempDir::new().unwrap();
    let tp = temp.path().to_path_buf();
    let input = tp.join("in.fastq");
    fs::write(
        &input,
        "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nTTTTGGGG\n+\nHHHHHHHH\n",
    )
    .unwrap();
    let archive = tp.join("a.qz");
    qz_lib::compression::compress(&CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: tp.clone(),
        threads: 1,
        ..CompressConfig::default()
    })
    .unwrap();

    // Sanity: untampered fast verify succeeds.
    let v_ok = qz_lib::cli::VerifyConfig {
        input: archive.clone(),
        working_dir: tp.clone(),
        num_threads: 1,
        fast: true,
    };
    assert!(qz_lib::compression::verify(&v_ok).is_ok());

    // Tamper the stored header_size (u32 LE at byte offset 4) by +1.
    let mut bytes = fs::read(&archive).unwrap();
    let hs = u32::from_le_bytes(bytes[4..8].try_into().unwrap());
    bytes[4..8].copy_from_slice(&(hs + 1).to_le_bytes());
    let tampered = tp.join("tampered.qz");
    fs::write(&tampered, &bytes).unwrap();

    // Fast verify must now reject, consistent with the deep decode path.
    let v_bad = qz_lib::cli::VerifyConfig {
        input: tampered.clone(),
        working_dir: tp.clone(),
        num_threads: 1,
        fast: true,
    };
    assert!(
        qz_lib::compression::verify(&v_bad).is_err(),
        "fast verify must reject a tampered header_size"
    );
    // And decompression must also fail (no silent fast/deep divergence).
    let d = decompress_args(tampered, tp.join("out.fastq"), tp.clone());
    assert!(qz_lib::compression::decompress(&d).is_err());
}

#[test]
fn test_corrupt_archive_does_not_destroy_existing_output() {
    // Decompress writes to a sibling temp and renames on success, so a corrupt
    // archive must neither leave a partial output nor (with --force) destroy a
    // pre-existing valid file at the target path.
    let temp = TempDir::new().unwrap();
    let tp = temp.path().to_path_buf();
    let input = tp.join("in.fastq");
    fs::write(
        &input,
        "@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@r2\nTTTTGGGGAAAA\n+\nHHHHHHHHHHHH\n",
    )
    .unwrap();
    let archive = tp.join("a.qz");
    qz_lib::compression::compress(&CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: tp.clone(),
        threads: 1,
        ..CompressConfig::default()
    })
    .unwrap();

    // Corrupt a byte deep in the archive body (a compressed stream byte) so decode
    // fails partway, after the output writer would have been opened.
    let mut bytes = fs::read(&archive).unwrap();
    let n = bytes.len();
    bytes[n - 4] ^= 0xFF;
    let corrupt = tp.join("corrupt.qz");
    fs::write(&corrupt, &bytes).unwrap();

    // Case 1: no existing output — a failed decode must leave the target absent.
    let out = tp.join("out.fastq");
    let d = decompress_args(corrupt.clone(), out.clone(), tp.clone());
    assert!(qz_lib::compression::decompress(&d).is_err());
    assert!(
        !out.exists(),
        "corrupt-archive decode left a partial output file"
    );

    // Case 2: a valid existing output with --force must survive a failed decode.
    fs::write(&out, "PRESERVE ME").unwrap();
    let mut d = decompress_args(corrupt, out.clone(), tp.clone());
    d.force = true;
    assert!(qz_lib::compression::decompress(&d).is_err());
    assert_eq!(fs::read_to_string(&out).unwrap(), "PRESERVE ME");
}

#[test]
fn test_fasta_roundtrip_byte_exact() {
    // FASTA input must round-trip byte-exact as FASTA (no quality lines).
    // The archive records a FASTA flag so decompress emits FASTA, not FASTQ.
    let data = ">seq1\nACGTACGT\n>seq2 some description\nTTTTGGGG\n>seq3\nAAAACCCC\n";
    roundtrip_exact(data, |c| {
        c.fasta = true;
    });
}

#[test]
fn test_fasta_illumina_shaped_headers_byte_exact() {
    // Regression: a '>'-prefixed FASTA header that is Casava/Illumina-shaped
    // (7 colon fields, numeric lane/tile/x/y) used to be misclassified by the
    // default columnar header path, which strips only '@' and re-adds a literal
    // '@' on decode — corrupting '>id' into '@>id'. Must round-trip byte-exact.
    let data = ">INSTR:RUN:FLOWCELL:1:1000:2000:3000 desc\nACGTACGT\n\
                >INSTR:RUN:FLOWCELL:1:1001:2001:3001 desc\nTTTTGGGG\n";
    roundtrip_exact(data, |c| {
        c.fasta = true;
    });
}

#[test]
fn test_fasta_ultra_roundtrip_byte_exact() {
    // Regression: the ultra/local-reorder path writes its own archive header and
    // used to hardcode the flags byte to 0, dropping FLAG_FASTA. A FASTA archive
    // produced with --ultra then decompressed as FASTQ (with synthesized empty
    // quality / '+' separator lines). FASTA must round-trip as FASTA on every
    // compress path, including ultra.
    let mut data = String::new();
    for i in 0..2000 {
        data.push_str(&format!(">r{i}\nACGTACGTACGTACGT\n"));
    }
    roundtrip_exact(&data, |c| {
        c.fasta = true;
        c.ultra = Some(0);
    });
}

#[test]
fn test_large_fasta_roundtrip() {
    // Regression: `--fasta` on >= MIN_READS_QUALITY_CTX (100k) records used to
    // crash because the Auto quality compressor selected a context-modeled quality
    // codec, which then errored on the empty quality alphabet ("no quality symbols
    // found"). FASTA forces `no_quality` (→ BSC), so it must compress and decompress
    // (as FASTA) without error at any size, with every header and sequence intact.
    const N: u32 = 100_001;
    let mut data = String::with_capacity(N as usize * 12);
    for i in 0..N {
        data.push_str(">r");
        data.push_str(&i.to_string());
        data.push_str("\nACGTACGT\n");
    }
    let (_original, decompressed) = roundtrip(&data, |c| {
        c.fasta = true;
    });
    // Every FASTA header and its sequence must be present after roundtrip.
    let ids: Vec<&str> = decompressed
        .lines()
        .filter(|l| l.starts_with('>'))
        .collect();
    assert_eq!(ids.len(), N as usize, "all FASTA headers must survive");
    assert_eq!(ids[0], ">r0");
    assert_eq!(ids[N as usize - 1], ">r100000");
    assert_eq!(
        decompressed.matches("ACGTACGT").count(),
        N as usize,
        "all sequences must survive"
    );
}

#[test]
fn test_compress_decompress_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nHHHHHHHH\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let original_lines: Vec<&str> = original.lines().collect();
    let decompressed_lines: Vec<&str> = decompressed.lines().collect();

    assert_eq!(original_lines.len(), decompressed_lines.len());

    for (orig, decomp) in original_lines.iter().zip(decompressed_lines.iter()) {
        assert_eq!(orig.trim(), decomp.trim());
    }
}

#[test]
fn test_decompress_gzipped_output() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq.gz");
    let d_args = DecompressConfig {
        input: archive_path,
        output: vec![output_fastq.clone()],
        working_dir: temp_path,
        num_threads: 1,
        gzipped: true,
        gzip_level: 6,
        force: false,
    };
    qz_lib::compression::decompress(&d_args).unwrap();

    assert!(output_fastq.exists());

    let data = fs::read(&output_fastq).unwrap();
    assert!(data.len() >= 2);
    assert_eq!(data[0], 0x1f);
    assert_eq!(data[1], 0x8b);
}

#[test]
fn test_decompress_no_quality() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        no_quality: true,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    // no_quality discards quality: the 4-line FASTQ shape is preserved with the
    // sequence intact and an EMPTY quality line (matches assemble_fastq_record's
    // QualityInput::None canonical output). The original quality must be gone.
    assert_eq!(content, "@read1\nACGTACGT\n+\n\n");
    assert!(!content.contains("IIIIIIII"));
}

// ========================================
// ERROR HANDLING TESTS
// ========================================

#[test]
fn test_empty_file() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("empty.fastq");
    fs::write(&input_fastq, "").unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };

    // Empty file should either compress successfully to an empty archive,
    // or return a clear error — both are acceptable.
    let result = qz_lib::compression::compress(&compress_args);
    match result {
        Ok(()) => {
            assert!(
                archive_path.exists(),
                "Empty archive should exist after successful compress"
            );
        }
        Err(e) => {
            let msg = e.to_string().to_lowercase();
            assert!(
                msg.contains("empty") || msg.contains("no reads") || msg.contains("0 reads"),
                "Empty file error should mention emptiness, got: {}",
                msg
            );
        }
    }
}

#[test]
fn test_single_read() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("single.fastq");
    fs::write(&input_fastq, "@read1\nACGT\n+\nIIII\n").unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains("read1"));
    assert!(content.contains("ACGT"));
}

#[test]
fn test_long_read() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let long_seq = "ACGT".repeat(250);
    let long_qual = "I".repeat(1000);
    let fastq_data = format!("@long_read\n{}\n+\n{}\n", long_seq, long_qual);

    let input_fastq = temp_path.join("long.fastq");
    fs::write(&input_fastq, fastq_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains(&long_seq));
}

#[test]
fn test_reads_with_n_bases() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let test_data = "@read1\nACGTNNNNACGT\n+\nIIII####IIII\n@read2\nNNNNNNNN\n+\n########\n";
    let input_fastq = temp_path.join("with_n.fastq");
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains("ACGTNNNNA")); // N bases should be preserved
}

// ========================================
// QUALITY MODE TESTS
// ========================================

#[test]
fn test_illumina_binning_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let test_data = "@read1\nACGTACGT\n+\nIIHHGGFF\n";
    let input_fastq = temp_path.join("test.fastq");
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        quality_mode: QualityMode::IlluminaBin,
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    assert!(output_fastq.exists());
    let content = fs::read_to_string(&output_fastq).unwrap();
    // Sequence must survive intact
    assert!(
        content.contains("ACGTACGT"),
        "sequence missing from illumina-binned output"
    );
    // Header must survive
    assert!(
        content.contains("@read1"),
        "header missing from illumina-binned output"
    );
    // Quality must be present (binned, not identical, but same length as sequence)
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 4, "expected 4 lines in output FASTQ");
    assert_eq!(
        lines[1].len(),
        lines[3].len(),
        "quality length ({}) must match sequence length ({})",
        lines[3].len(),
        lines[1].len()
    );
}

#[test]
fn test_binary_quality_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let test_data = "@read1\nACGTACGT\n+\nIIHHGGFF\n";
    let input_fastq = temp_path.join("test.fastq");
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        quality_mode: QualityMode::Binary,
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("output.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    assert!(output_fastq.exists());
    let content = fs::read_to_string(&output_fastq).unwrap();
    // Sequence and header must survive intact
    assert!(
        content.contains("ACGTACGT"),
        "sequence missing from binary quality output"
    );
    assert!(
        content.contains("@read1"),
        "header missing from binary quality output"
    );
    // Quality must be present and same length as sequence
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 4, "expected 4 lines in output FASTQ");
    assert_eq!(
        lines[1].len(),
        lines[3].len(),
        "quality length ({}) must match sequence length ({})",
        lines[3].len(),
        lines[1].len()
    );
}

// ========================================
// MULTI-THREADING TEST
// ========================================

#[test]
fn test_multithreading_deterministic() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n".repeat(100);
    let input_fastq = temp_path.join("test.fastq");
    fs::write(&input_fastq, &test_data).unwrap();

    // Compress with 1 thread
    let archive1 = temp_path.join("test1.qz");
    let compress_args1 = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive1.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args1).unwrap();

    // Compress with 4 threads
    let archive2 = temp_path.join("test2.qz");
    let compress_args2 = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive2.clone(),
        working_dir: temp_path.clone(),
        threads: 4,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args2).unwrap();

    // Decompress both
    let output1 = temp_path.join("output1.fastq");
    let output2 = temp_path.join("output2.fastq");

    qz_lib::compression::decompress(&decompress_args(
        archive1,
        output1.clone(),
        temp_path.clone(),
    ))
    .unwrap();
    qz_lib::compression::decompress(&decompress_args(archive2, output2.clone(), temp_path))
        .unwrap();

    // Both should produce identical output
    let content1 = fs::read_to_string(&output1).unwrap();
    let content2 = fs::read_to_string(&output2).unwrap();
    assert_eq!(
        content1, content2,
        "Multi-threaded compression should produce same output as single-threaded"
    );
}

// ========== ERROR HANDLING TESTS ==========

#[test]
fn test_decompress_truncated_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Create a valid archive
    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nHHHHHHHH\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    // Truncate to half size
    let archive_data = fs::read(&archive_path).unwrap();
    let truncated_path = temp_path.join("truncated.qz");
    fs::write(&truncated_path, &archive_data[..archive_data.len() / 2]).unwrap();

    let output_fastq = temp_path.join("output.fastq");
    let result =
        qz_lib::compression::decompress(&decompress_args(truncated_path, output_fastq, temp_path));
    assert!(
        result.is_err(),
        "Decompressing a truncated archive should fail"
    );
}

#[test]
fn test_decompress_corrupted_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Create a valid archive
    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nHHHHHHHH\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    // Corrupt bytes in the data section
    let mut archive_data = fs::read(&archive_path).unwrap();
    let mid = archive_data.len() / 2;
    for i in mid..std::cmp::min(mid + 32, archive_data.len()) {
        archive_data[i] ^= 0xFF;
    }
    let corrupted_path = temp_path.join("corrupted.qz");
    fs::write(&corrupted_path, &archive_data).unwrap();

    let output_fastq = temp_path.join("output.fastq");
    let result =
        qz_lib::compression::decompress(&decompress_args(corrupted_path, output_fastq, temp_path));
    assert!(
        result.is_err(),
        "Decompressing a corrupted archive should fail"
    );
}

#[test]
fn test_decompress_zero_byte_file() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let empty_archive = temp_path.join("empty.qz");
    fs::write(&empty_archive, b"").unwrap();

    let output_fastq = temp_path.join("output.fastq");
    let result =
        qz_lib::compression::decompress(&decompress_args(empty_archive, output_fastq, temp_path));
    assert!(
        result.is_err(),
        "Decompressing a zero-byte file should fail"
    );
}

#[test]
fn test_fqzcomp_quality_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Multiple reads with varying quality profiles to test mean-quality sorting
    let test_data = "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n\
                     @read2\nTGCATGCATGCATGCA\n+\n!!!!!!!!!!!!!!!!\n\
                     @read3\nAAAACCCCGGGGTTTT\n+\nFFFFFFFFFFFFFFFF\n\
                     @read4\nACGTTGCAACGTTGCA\n+\nHHHHAAAAHHHHAAAA\n\
                     @read5\nGGGGCCCCAAAATTTT\n+\nBBBBBBBBBBBBBBBB\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions {
            quality_compressor: QualityCompressor::Fqzcomp,
            ..Default::default()
        },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let original_lines: Vec<&str> = original.lines().collect();
    let decompressed_lines: Vec<&str> = decompressed.lines().collect();

    assert_eq!(
        original_lines.len(),
        decompressed_lines.len(),
        "Line count mismatch: expected {}, got {}",
        original_lines.len(),
        decompressed_lines.len()
    );

    for (i, (orig, decomp)) in original_lines
        .iter()
        .zip(decompressed_lines.iter())
        .enumerate()
    {
        assert_eq!(
            orig.trim(),
            decomp.trim(),
            "Line {} mismatch:\n  original:     {:?}\n  decompressed: {:?}",
            i,
            orig,
            decomp
        );
    }
}

#[test]
fn test_fqzcomp_quality_record_capped_multiblock_roundtrip() {
    // Force the single-end fqz quality encoder to emit MORE THAN ONE
    // record-capped block by lowering `quality_ctx_block_size` below the read
    // count, then assert byte-identical roundtrip through the real chunked
    // compress/decompress pipeline. This exercises the v4 multi-block decode
    // path over multiple independently-(locally-)sorted fqz blocks.
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // 50 reads, cap 10 -> 5 blocks. Vary quality per read so the per-block
    // local mean-quality sort actually permutes (and must be unsorted on decode).
    let mut data = String::new();
    for i in 0..50u32 {
        let q = (b'!' + (i % 40) as u8) as char;
        let qline: String = std::iter::repeat_n(q, 24).collect();
        data.push_str(&format!("@read{i}\nACGTACGTACGTACGTACGTACGT\n+\n{qline}\n"));
    }
    fs::write(&input_fastq, &data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions {
            quality_compressor: QualityCompressor::Fqzcomp,
            // Force >1 fqz block from only 50 reads.
            quality_ctx_block_size: 10,
            ..Default::default()
        },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    assert_eq!(
        original, decompressed,
        "record-capped multi-block fqz roundtrip must be byte-identical"
    );
}

#[test]
fn test_fqzcomp_quality_variable_length_streaming_roundtrip() {
    // Fqzcomp quality now decodes through the BOUNDED-STREAMING path.
    // With VARIABLE read lengths the archive
    // carries no constant-length optimization (const_seq_len == const_qual_len
    // == 0), so both the sequence and the fqz quality streams are decoded via
    // their per-record varint framing. Multi-block (cap < read count) + varying
    // lengths + varying quality exercises the StreamCodec::Fqz producer and the
    // QualityInput::Decoded raw-ASCII emit. Must be byte-identical.
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("var_input.fastq");
    let mut data = String::new();
    for i in 0..60u32 {
        // Length varies 20..=49, so no constant-length framing.
        let len = 20 + (i % 30) as usize;
        let seq: String = (0..len)
            .map(|j| b"ACGT"[(j + i as usize) % 4] as char)
            .collect();
        // Quality varies per read AND within the read.
        let qline: String = (0..len)
            .map(|j| (b'!' + ((i as usize + j) % 40) as u8) as char)
            .collect();
        data.push_str(&format!("@read{i} extra\n{seq}\n+\n{qline}\n"));
    }
    fs::write(&input_fastq, &data).unwrap();

    let archive_path = temp_path.join("var.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions {
            quality_compressor: QualityCompressor::Fqzcomp,
            quality_ctx_block_size: 16, // 60 reads -> 4 blocks
            ..Default::default()
        },
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("var_out.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    assert_eq!(
        fs::read_to_string(&input_fastq).unwrap(),
        fs::read_to_string(&output_fastq).unwrap(),
        "variable-length multi-block fqz streaming roundtrip must be byte-identical"
    );
}

#[test]
fn test_sequence_hints_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Reads >= 21bp for valid syncmers, plus one short read to test 0xFF fallback
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n+\nBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB\n\
@read4\nACGT\n+\nIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions {
            sequence_hints: true,
            ..Default::default()
        },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify encoding_type byte is 4 (at offset 8, after v2 prefix)
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(
        archive_data[8], 4,
        "encoding_type should be 4 for sequence hints"
    );

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len());
    for (orig, decomp) in orig_lines.iter().zip(decomp_lines.iter()) {
        assert_eq!(orig.trim(), decomp.trim());
    }
}

#[test]
fn test_rc_canon_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();
    let input_fastq = temp_path.join("test.fastq");

    // Create test data with reads that benefit from RC canonicalization
    fs::write(
        &input_fastq,
        "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n\
         @read2\nACGTACGTACGTACGT\n+\nFFFFFFFFFFFFFFFF\n\
         @read3\nTTTTAAAACCCCGGGG\n+\nAAAAAAAAAAAAAAAA\n\
         @read4\nCCCCGGGGTTTTAAAA\n+\nBBBBBBBBBBBBBBBB\n",
    )
    .unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions {
            rc_canon: true,
            ..Default::default()
        },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify encoding_type byte is 6 (at offset 8, after v2 prefix)
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[8], 6, "encoding_type should be 6 for rc_canon");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (orig, decomp)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(orig.trim(), decomp.trim(), "mismatch at line {}", i);
    }
}

#[test]
fn test_rc_canon_roundtrip_lowercase_and_iupac() {
    // Regression: rc_canon previously used a complement that folded lowercase to
    // uppercase and non-ACGT (IUPAC) to N, silently corrupting soft-masked /
    // ambiguous reads when canonicalization reverse-complemented them. The fix
    // uses a case-preserving involution. Use revcomp-favoring reads (start with
    // T/G) so rc_flag=1 actually fires for the lowercase/IUPAC content.
    let data = "@r1\ntttggggcccaaaa\n+\nIIIIIIIIIIIIII\n\
                @r2\nttttacgtacgtRY\n+\nFFFFFFFFFFFFFF\n\
                @r3\ngggannnnttttac\n+\nAAAAAAAAAAAAAA\n\
                @r4\nNRYSWKMacgtTAC\n+\nBBBBBBBBBBBBBB\n";
    roundtrip_exact(data, |c| {
        c.advanced.rc_canon = true;
    });
}

#[test]
fn test_rc_canon_rejects_sequence_hints() {
    // Regression: rc_canon + sequence_hints was an accepted config that emitted a
    // hint byte the type-6 decoder never strips, silently corrupting every
    // sequence. Must be rejected at config validation.
    let data = "@r1\nACGTACGT\n+\nIIIIIIII\n";
    let temp = TempDir::new().unwrap();
    let tp = temp.path().to_path_buf();
    let input = tp.join("in.fastq");
    fs::write(&input, data).unwrap();
    let cfg = CompressConfig {
        input: vec![input.clone()],
        output: tp.join("out.qz"),
        working_dir: tp.clone(),
        threads: 1,
        advanced: AdvancedOptions {
            rc_canon: true,
            sequence_hints: true,
            ..Default::default()
        },
        ..CompressConfig::default()
    };
    assert!(
        qz_lib::compression::compress(&cfg).is_err(),
        "rc_canon + sequence_hints must be rejected"
    );
}

#[test]
fn test_ultra_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read4\nGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n\
@read5\nACGTACGTACGTACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read6\nTGCATGCATGCATGCATGCATGCATGCATGCT\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ultra: Some(3),
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(
        archive_data[8], 10,
        "encoding_type should be 10 (UltraBigBlock) for ultra"
    );

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .unwrap();
    assert!(output_fastq.exists());

    // Ultra preserves original order — verify exact match
    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (o, d)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(o, d, "mismatch at line {}", i);
    }
}

/// Ultra's reorder-local path stores each chunk's sequence stream in BSC blocks
/// capped at `max_seq_block_mb` (750 MB), so a realistic chunk produces a single
/// sequence block whose *decompressed* size far exceeds 64 MiB. Decompression
/// must accept such blocks. (Regression: the `BSC_MAX_BLOCK_LEN` = 64 MiB DoS
/// guard rejected them, making every ultra archive on non-trivial input
/// undecompressable — compress reported success but decompress failed with
/// "decompressed block size … exceeds maximum … (corrupt header)".)
///
/// We generate ~80 MB of sequence (> 64 MiB) in a single ultra level-1 chunk
/// (chunk_size = 1M reads) using long reads to keep the record count — and thus
/// the test's wall-clock — modest. Content is repetitive so BSC compresses it
/// fast; the block's recorded *decompressed* size is what trips the guard.
#[test]
fn test_ultra_large_sequence_block_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // 130k reads × 600 bp = 78 MB of sequence in one chunk — comfortably over
    // the 64 MiB (67,112,960 B) BSC block cap.
    const N_READS: usize = 130_000;
    const READ_LEN: usize = 600;
    let bases = [b'A', b'C', b'G', b'T'];
    let mut data = Vec::with_capacity(N_READS * (READ_LEN * 2 + 32));
    for i in 0..N_READS {
        let seq: Vec<u8> = (0..READ_LEN).map(|j| bases[(i + j) % 4]).collect();
        data.extend_from_slice(format!("@read{i}\n").as_bytes());
        data.extend_from_slice(&seq);
        data.extend_from_slice(b"\n+\n");
        data.extend(std::iter::repeat_n(b'I', READ_LEN));
        data.push(b'\n');
    }

    let input_fastq = temp_path.join("big.fastq");
    fs::write(&input_fastq, &data).unwrap();

    let archive_path = temp_path.join("big.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 4,
        ultra: Some(1),
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());
    assert_eq!(
        fs::read(&archive_path).unwrap()[8],
        10,
        "encoding_type should be 10 (UltraBigBlock)"
    );

    let output_fastq = temp_path.join("big.out.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path,
        output_fastq.clone(),
        temp_path,
    ))
    .expect("ultra decompress of a >64 MiB sequence block must succeed");

    let original = fs::read(&input_fastq).unwrap();
    let decompressed = fs::read(&output_fastq).unwrap();
    assert_eq!(
        original, decompressed,
        "ultra roundtrip must be byte-identical"
    );
}

#[test]
fn ultra_archive_is_ultra_big_block_encoding() {
    let dir = TempDir::new().unwrap();
    let p = dir.path().to_path_buf();
    let input = p.join("in.fastq");
    fs::write(&input, MULTI_READ_DATA).unwrap();
    let archive = p.join("a.qz");
    let cfg = CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: p,
        threads: 2,
        ultra: Some(1),
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();
    assert_eq!(
        fs::read(&archive).unwrap()[8],
        10,
        "ultra must emit encoding_type 10"
    );
}

// ========================================
// BSC STATIC MODE TEST
// ========================================

#[test]
fn test_bsc_static_roundtrip() {
    roundtrip_exact(MULTI_READ_DATA, |c| {
        c.advanced.bsc_static = true;
    });
}

// ========================================
// CORRUPT ARCHIVE TESTS
// ========================================

#[test]
fn test_decompress_invalid_compressor_codes() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, "@r1\nACGT\n+\nIIII\n").unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    // Byte 11 is quality_compressor code (offset 3 in body, after 8-byte v2 prefix) -- set to invalid 99
    let mut data = fs::read(&archive_path).unwrap();
    data[11] = 99;
    let bad_path = temp_path.join("bad.qz");
    fs::write(&bad_path, &data).unwrap();

    let output = temp_path.join("out.fastq");
    let result = qz_lib::compression::decompress(&decompress_args(bad_path, output, temp_path));
    assert!(result.is_err(), "Invalid compressor code should fail");
}

#[test]
fn test_compress_refuses_to_overwrite_without_force() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, "@r1\nACGT\n+\nIIII\n").unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    let original_size = fs::metadata(&archive_path).unwrap().len();

    // Second compress with the same output must refuse, not overwrite.
    let err = qz_lib::compression::compress(&compress_args).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("already exists"), "unexpected error: {msg}");
    assert_eq!(
        fs::metadata(&archive_path).unwrap().len(),
        original_size,
        "existing archive must be left untouched"
    );

    // With --force, it should overwrite cleanly.
    let mut forced = compress_args.clone();
    forced.force = true;
    qz_lib::compression::compress(&forced).unwrap();
    assert!(archive_path.exists());

    // And no .tmp sibling should be left behind on success.
    let tmp = archive_path.with_file_name("test.qz.tmp");
    assert!(!tmp.exists(), "tmp file should be cleaned up on success");
}

#[test]
fn test_decompress_detects_bit_flip_via_block_crc32() {
    // Flip a byte somewhere in the archive. The v4 per-block CRC32 (or the
    // block-length sanity check for header-field flips) must catch this with
    // a clear error — no silent wrong output and no cryptic block_info failure.
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("input.fastq");
    // Need enough data that BSC produces a non-trivial compressed payload.
    let mut data = String::new();
    for i in 0..200 {
        data.push_str(&format!(
            "@r{i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    fs::write(&input_fastq, &data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    // Sanity check: untouched archive must decompress cleanly.
    let output = temp_path.join("ok.fastq");
    qz_lib::compression::decompress(&decompress_args(
        archive_path.clone(),
        output,
        temp_path.clone(),
    ))
    .unwrap();

    // Flip a single byte two-thirds of the way into the archive — inside
    // a compressed stream region. Depending on which field is flipped, the
    // decoder reports either a CRC32 mismatch (payload flip) or a block-size
    // sanity error (header field flip). Both indicate the corruption was caught.
    let mut bytes = fs::read(&archive_path).unwrap();
    let flip_pos = bytes.len() * 2 / 3;
    bytes[flip_pos] ^= 0xFF;
    let bad_path = temp_path.join("bad.qz");
    fs::write(&bad_path, &bytes).unwrap();

    let bad_output = temp_path.join("out.fastq");
    let result = qz_lib::compression::decompress(&decompress_args(bad_path, bad_output, temp_path));
    let err = result.expect_err("bit-flipped archive must error, not silently succeed");
    let msg = format!("{err:#}");
    // The bounded path may surface corruption as a CRC32 mismatch (payload flip),
    // a block-size sanity error (header-field flip: "extends past stream"), or a
    // BSC decompressor error ("truncated" / "corrupt"). All indicate the corruption
    // was caught — accept any of these diagnostic forms.
    assert!(
        msg.contains("CRC32 mismatch")
            || msg.contains("corrupt")
            || msg.contains("truncated")
            || msg.contains("extends past stream")
            || msg.contains("block"),
        "expected corruption-detection error, got: {msg}"
    );
}

#[test]
fn test_decompress_rejects_v2_archive_with_helpful_message() {
    // Construct a minimal byte buffer claiming to be a v2 archive. qz only reads
    // v5 chunk-major archives now, so any older/legacy version must be refused
    // with a clear recompress hint rather than fail deep inside the parser.
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let bad_path = temp_path.join("v2.qz");
    // Magic "QZ", version=2, reserved=0, header_size=60 (enough to look plausible),
    // padded with zeros so the file is large enough to pass the initial size check.
    let mut data = vec![b'Q', b'Z', 2u8, 0u8];
    data.extend_from_slice(&60u32.to_le_bytes());
    data.resize(200, 0); // pad
    fs::write(&bad_path, &data).unwrap();

    let output = temp_path.join("out.fastq");
    let result = qz_lib::compression::decompress(&decompress_args(bad_path, output, temp_path));
    let err = result.expect_err("legacy (non-v5) archive must be rejected");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("no longer supported"),
        "expected a legacy-archive recompress hint, got: {msg}"
    );
}

#[test]
fn test_verify_fast_rc_canon_truncation_fails() {
    // Round-5 fix: verify --fast previously returned OK when an rc_canon
    // (encoding_type=6) archive was truncated exactly before the 8-byte
    // rc_flags length field. That stream is REQUIRED to reconstruct
    // sequences, so missing-or-truncated is corruption, not a legitimate
    // skip-condition. Fast verify must now bail.
    let mut data = String::new();
    for i in 0..50 {
        data.push_str(&format!("@r{i}\nACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFF\n"));
    }
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();
    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, &data).unwrap();

    let archive_path = temp_path.join("out.qz");
    qz_lib::compression::compress(&CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions {
            rc_canon: true,
            ..Default::default()
        },
        ..CompressConfig::default()
    })
    .expect("rc_canon compression should succeed");

    // Truncate the archive just enough that the rc_flags length field is
    // missing. Read the v3 header to compute exactly where qualities end,
    // then truncate to that offset.
    let archive_bytes = fs::read(&archive_path).unwrap();
    // Use the documented v3 layout: prefix(8) + body. Strip enough that the
    // 8-byte rc_flags_len at the end is gone. We don't need a perfect cut —
    // any size where `q_end + 8 > data.len()` triggers the new check, so
    // cutting the trailing 9 bytes guarantees we land in the gap.
    let truncated = &archive_bytes[..archive_bytes.len() - 9];
    fs::write(&archive_path, truncated).unwrap();

    let result = qz_lib::compression::verify(&qz_lib::cli::VerifyConfig {
        input: archive_path,
        working_dir: temp_path,
        num_threads: 1,
        fast: true,
    });
    let err = result.expect_err("truncated rc_canon archive must fail fast verify");
    let msg = format!("{err:#}");
    // v4 archives bail with "rc_canon archive truncated"; v5 chunk-major
    // archives lose their trailing directory locator/footer when the tail is
    // cut, so truncation trips the footer/locator integrity gate instead. Both
    // are corruption — fast verify must bail either way.
    assert!(
        msg.contains("rc_canon archive truncated")
            || msg.contains("rc_flags")
            || msg.contains("footer")
            || msg.contains("locator")
            || msg.contains("too small"),
        "expected rc_canon truncation error, got: {msg}"
    );
}

#[test]
fn test_decompress_too_small_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let small_path = temp_path.join("small.qz");
    fs::write(&small_path, [0u8; 20]).unwrap();

    let output = temp_path.join("out.fastq");
    let result = qz_lib::compression::decompress(&decompress_args(small_path, output, temp_path));
    assert!(result.is_err(), "Archive smaller than header should fail");
}

// ========================================
// EDGE CASE TESTS
// ========================================

#[test]
fn test_all_n_bases() {
    roundtrip_exact(
        "@r1\nNNNNNNNN\n+\nIIIIIIII\n@r2\nNNNNNNNN\n+\nHHHHHHHH\n",
        |_| {},
    );
}

#[test]
fn test_single_base_reads() {
    roundtrip_exact(
        "@r1\nA\n+\nI\n@r2\nC\n+\nH\n@r3\nG\n+\nF\n@r4\nT\n+\nB\n",
        |_| {},
    );
}

#[test]
fn test_mixed_length_reads() {
    roundtrip_exact(
        "@r1\nA\n+\nI\n@r2\nACGT\n+\nIIII\n@r3\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n",
        |_| {},
    );
}

#[test]
fn test_three_input_files_rejected() {
    // Two inputs is now the paired-end happy path; three or more is still rejected.
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let f1 = temp_path.join("r1.fastq");
    let f2 = temp_path.join("r2.fastq");
    let f3 = temp_path.join("r3.fastq");
    fs::write(&f1, "@r\nACGT\n+\nIIII\n").unwrap();
    fs::write(&f2, "@r\nACGT\n+\nIIII\n").unwrap();
    fs::write(&f3, "@r\nACGT\n+\nIIII\n").unwrap();

    let compress_args = CompressConfig {
        input: vec![f1, f2, f3],
        output: temp_path.join("test.qz"),
        working_dir: temp_path,
        threads: 1,
        ..CompressConfig::default()
    };

    let result = qz_lib::compression::compress(&compress_args);
    assert!(result.is_err(), "Three input files should be rejected");
}

/// Paired-end compression through the PUBLIC dispatch now produces a v5
/// chunk-major container (archive_type=1), and decompresses both mates
/// byte-identically. Forcing several small chunks also exercises the bounded
/// per-chunk decode path (Task 2.5 Step 7a + 7c: the multi-chunk roundtrip
/// decodes one chunk at a time, which is the bounded-RAM property — a peak-RSS
/// probe would be flaky under parallel test load, so this is the lighter,
/// deterministic assertion).
#[test]
fn test_paired_public_dispatch_is_v5_and_roundtrips() {
    let temp_dir = TempDir::new().unwrap();
    let tp = temp_dir.path().to_path_buf();

    let r1 = tp.join("R1.fastq");
    let r2 = tp.join("R2.fastq");
    let mut s1 = String::new();
    let mut s2 = String::new();
    for i in 0..40 {
        s1.push_str(&format!("@m:{i} 1:N:0:ACGT\nACGTACGTAN\n+\nIIIIIIIIII\n"));
        s2.push_str(&format!("@m:{i} 2:N:0:ACGT\nNTTTTGGGGC\n+\n#IIIIIIIII\n"));
    }
    fs::write(&r1, &s1).unwrap();
    fs::write(&r2, &s2).unwrap();

    let out = tp.join("paired.qz");
    let mut cfg = CompressConfig {
        input: vec![r1, r2],
        output: out.clone(),
        working_dir: tp.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    cfg.advanced.chunk_records = 8; // 40 pairs / 8 → 5 chunks, multi-chunk decode
    qz_lib::compression::compress(&cfg).unwrap();

    // v5 paired container: "QZ" magic, version 5, archive_type 1, trailing locator.
    let b = fs::read(&out).unwrap();
    assert_eq!(&b[0..2], b"QZ");
    assert_eq!(b[2], 5, "paired archive must be v5");
    assert_eq!(b[3], 1, "v5 archive_type must be 1 (paired)");
    assert_eq!(&b[b.len() - 8..], b"QZFOOTR1", "missing v5 footer magic");
    assert!(qz_lib::compression::is_paired_archive(&out).unwrap());
    assert!(!qz_lib::compression::is_reference_archive(&out).unwrap());

    // Decompress to a prefix; both mates must be byte-identical.
    let prefix = tp.join("dec");
    let dc = DecompressConfig {
        input: out,
        output: vec![prefix],
        working_dir: tp.clone(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&dc).unwrap();
    assert_eq!(fs::read_to_string(tp.join("dec_R1.fastq")).unwrap(), s1);
    assert_eq!(fs::read_to_string(tp.join("dec_R2.fastq")).unwrap(), s2);

    // Deep + fast verify both route through verify_paired_v5 via the public API.
    let vc = qz_lib::cli::VerifyConfig {
        input: tp.join("paired.qz"),
        working_dir: tp.clone(),
        num_threads: 1,
        fast: false,
    };
    let deep = qz_lib::compression::verify(&vc).unwrap();
    assert_eq!(deep.num_reads, 80, "40 pairs * 2 mates");
    assert_ne!(deep.crc32, 0);
    let vc_fast = qz_lib::cli::VerifyConfig { fast: true, ..vc };
    let fast = qz_lib::compression::verify(&vc_fast).unwrap();
    assert!(fast.blocks_verified > 0);
}

// ========== VERIFY TESTS ==========

#[test]
fn test_verify_valid_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, MULTI_READ_DATA).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: archive_path,
        working_dir: temp_path,
        num_threads: 1,
        ..Default::default()
    };
    let result = qz_lib::compression::verify(&verify_config).unwrap();
    assert_eq!(result.num_reads, 3);
    assert!(result.total_bytes > 0);
    assert!(result.crc32 != 0);
    assert_eq!(result.mode, qz_lib::compression::VerifyMode::Deep);
    assert_eq!(result.blocks_verified, 0); // deep mode reports 0 (uses codec, not block walker)
}

#[test]
fn test_verify_fast_valid_archive() {
    // Fast verify: walks per-block CRC32s without invoking BSC decompression.
    // Must succeed cleanly on a healthy archive and report block-count metadata.
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();
    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, MULTI_READ_DATA).unwrap();

    let archive_path = temp_path.join("test.qz");
    qz_lib::compression::compress(&CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    })
    .unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: archive_path,
        working_dir: temp_path,
        num_threads: 1,
        fast: true,
    };
    let result = qz_lib::compression::verify(&verify_config).unwrap();
    assert_eq!(result.mode, qz_lib::compression::VerifyMode::Fast);
    assert_eq!(result.crc32, 0); // not computed in fast mode
    assert!(
        result.blocks_verified > 0,
        "fast verify should report block count"
    );
    assert!(result.total_bytes > 0);
}

#[test]
fn verify_fast_supports_ultra_archive() {
    let dir = TempDir::new().unwrap();
    let p = dir.path().to_path_buf();
    let input = p.join("in.fastq");
    fs::write(&input, MULTI_READ_DATA).unwrap();
    let archive = p.join("a.qz");
    let cfg = CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: p.clone(),
        threads: 2,
        ultra: Some(1),
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();
    let result = qz_lib::compression::verify(&qz_lib::cli::VerifyConfig {
        input: archive,
        working_dir: p,
        num_threads: 2,
        fast: true,
    })
    .expect("fast verify must succeed on an UltraBigBlock archive");
    assert!(result.blocks_verified > 0);
}

#[test]
fn test_verify_fast_covers_rc_canon_flags_stream() {
    // RC canon archives append a v3-framed rc_flags stream after qualities.
    // Fast verify must walk it — otherwise corruption in the RC flags region
    // would silently pass (the worst possible verify-tool failure mode).
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();
    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, MULTI_READ_DATA).unwrap();

    let archive_path = temp_path.join("test.qz");
    qz_lib::compression::compress(&CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions {
            rc_canon: true,
            ..Default::default()
        },
        ..CompressConfig::default()
    })
    .unwrap();

    // Flip a byte at the very end of the archive — that's inside the RC flags
    // payload region (or its CRC), past all other streams.
    let mut bytes = fs::read(&archive_path).unwrap();
    let flip_pos = bytes.len() - 20; // close to end of file
    bytes[flip_pos] ^= 0xFF;
    let bad_path = temp_path.join("bad.qz");
    fs::write(&bad_path, &bytes).unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: bad_path,
        working_dir: temp_path,
        num_threads: 1,
        fast: true,
    };
    let err = qz_lib::compression::verify(&verify_config)
        .expect_err("fast verify must catch RC flags corruption");
    let msg = format!("{err:#}");
    // The corruption lands near the end of the file. For v4 archives that's
    // inside the trailing rc_flags payload/CRC; for v5 chunk-major archives the
    // trailing 20 bytes are the directory locator, so flipping there trips the
    // footer integrity gate. Either way fast verify must reject the archive.
    assert!(
        msg.contains("CRC32 mismatch")
            || msg.contains("rc_flags")
            || msg.contains("footer")
            || msg.contains("locator"),
        "expected RC corruption to be caught, got: {msg}"
    );
}

#[test]
fn test_verify_fast_multi_block_stream() {
    // Confirm fast verify handles streams with >1 BSC block (the per-block
    // loop iteration logic). 25 MB of bases ≈ 1 block by default; force
    // smaller blocks via bsc_block_size_mb=1 to guarantee >1 block.
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();
    let input_fastq = temp_path.join("input.fastq");

    // Generate ~3 MB of FASTQ — split across multiple 1 MB BSC blocks.
    let mut data = String::with_capacity(3 * 1024 * 1024);
    for i in 0..30_000 {
        data.push_str(&format!(
            "@r{i}\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\n\
             IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    fs::write(&input_fastq, &data).unwrap();

    let archive_path = temp_path.join("test.qz");
    qz_lib::compression::compress(&CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions {
            bsc_block_size_mb: 1,
            ..Default::default()
        },
        ..CompressConfig::default()
    })
    .unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: archive_path,
        working_dir: temp_path,
        num_threads: 1,
        fast: true,
    };
    let result = qz_lib::compression::verify(&verify_config).unwrap();
    assert!(
        result.blocks_verified > 4,
        "expected multiple blocks per stream; got {}",
        result.blocks_verified
    );
}

#[test]
fn test_verify_fast_catches_bit_flip() {
    // The same bit-flip the deep path catches must also be caught by fast verify,
    // and faster (no BSC decompression on the bad block, just CRC check).
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();
    let input_fastq = temp_path.join("input.fastq");
    let mut data = String::new();
    for i in 0..100 {
        data.push_str(&format!(
            "@r{i}\nACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    fs::write(&input_fastq, &data).unwrap();

    let archive_path = temp_path.join("test.qz");
    qz_lib::compression::compress(&CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    })
    .unwrap();

    let mut bytes = fs::read(&archive_path).unwrap();
    // Flip a byte in the first third of the file. Block payloads live at the front
    // (front header → per-chunk blocks → ChunkDecodedSizes global → directory footer →
    // locator), so this lands inside a block payload and is caught by the per-block
    // CRC32 walk — not the footer-CRC check at the tail. (2/3 used to land in a block
    // payload but now lands in the footer region after the small ChunkDecodedSizes
    // global block was added before the footer.)
    let flip_pos = bytes.len() / 3;
    bytes[flip_pos] ^= 0xFF;
    let bad_path = temp_path.join("bad.qz");
    fs::write(&bad_path, &bytes).unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: bad_path,
        working_dir: temp_path,
        num_threads: 1,
        fast: true,
    };
    let err = qz_lib::compression::verify(&verify_config)
        .expect_err("fast verify must error on bit-flip");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("CRC32 mismatch"),
        "expected CRC32 mismatch, got: {msg}"
    );
}

#[test]
fn test_verify_corrupted_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, MULTI_READ_DATA).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    // Corrupt bytes in the data section
    let mut archive_data = fs::read(&archive_path).unwrap();
    let mid = archive_data.len() / 2;
    for i in mid..std::cmp::min(mid + 32, archive_data.len()) {
        archive_data[i] ^= 0xFF;
    }
    let corrupted_path = temp_path.join("corrupted.qz");
    fs::write(&corrupted_path, &archive_data).unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: corrupted_path,
        working_dir: temp_path,
        num_threads: 1,
        ..Default::default()
    };
    let result = qz_lib::compression::verify(&verify_config);
    assert!(result.is_err(), "Verifying a corrupted archive should fail");
}

#[test]
fn test_verify_truncated_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, MULTI_READ_DATA).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    let archive_data = fs::read(&archive_path).unwrap();
    let truncated_path = temp_path.join("truncated.qz");
    fs::write(&truncated_path, &archive_data[..archive_data.len() / 2]).unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: truncated_path,
        working_dir: temp_path,
        num_threads: 1,
        ..Default::default()
    };
    let result = qz_lib::compression::verify(&verify_config);
    assert!(result.is_err(), "Verifying a truncated archive should fail");
}

#[test]
fn test_bounded_streaming_broken_pipe_returns_ok() {
    use std::io::{Error, ErrorKind, Write};

    /// A writer that returns BrokenPipe after the first 100 bytes written.
    struct EarlyPipeClose {
        written: usize,
    }
    impl Write for EarlyPipeClose {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.written += buf.len();
            if self.written > 100 {
                Err(Error::new(ErrorKind::BrokenPipe, "consumer closed"))
            } else {
                Ok(buf.len())
            }
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    // Compress a small archive first.
    let mut data = String::new();
    for i in 0..200 {
        data.push_str(&format!("@r{i}\nACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFF\n"));
    }
    let temp = tempfile::TempDir::new().unwrap();
    let input_path = temp.path().join("in.fastq");
    let archive_path = temp.path().join("out.qz");
    std::fs::write(&input_path, &data).unwrap();
    qz_lib::compression::compress(&qz_lib::cli::CompressConfig {
        input: vec![input_path],
        output: archive_path.clone(),
        working_dir: temp.path().to_path_buf(),
        threads: 1,
        ..Default::default()
    })
    .unwrap();

    // Decompress with the EarlyPipeClose sink. Expect Ok(()), not Err.
    let mut writer = EarlyPipeClose { written: 0 };
    qz_lib::compression::decompress_to_writer(&archive_path, &mut writer)
        .expect("BrokenPipe must be swallowed by the library");
}

/// Empty FASTQ input is rejected cleanly at compress time with a diagnostic
/// error message. The implementation currently requires at least one record
/// (compress_impl bails with "Empty input file" before writing any archive).
/// This test verifies no panic occurs and the error is informative.
#[test]
fn test_bounded_streaming_empty_input() {
    let temp = tempfile::TempDir::new().unwrap();
    let input = temp.path().join("in.fastq");
    let archive = temp.path().join("out.qz");
    std::fs::write(&input, "").unwrap();
    let err = qz_lib::compression::compress(&qz_lib::cli::CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: temp.path().to_path_buf(),
        threads: 1,
        ..Default::default()
    })
    .unwrap_err();
    let msg = format!("{err:#}");
    assert!(
        msg.to_lowercase().contains("empty"),
        "expected empty-input error, got: {msg}"
    );
}

#[test]
fn test_bounded_streaming_no_quality() {
    let mut data = String::new();
    for i in 0..200 {
        data.push_str(&format!("@r{i}\nACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFF\n"));
    }
    let temp = tempfile::TempDir::new().unwrap();
    let input = temp.path().join("in.fastq");
    let archive = temp.path().join("out.qz");
    std::fs::write(&input, &data).unwrap();
    qz_lib::compression::compress(&qz_lib::cli::CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: temp.path().to_path_buf(),
        threads: 1,
        no_quality: true,
        ..Default::default()
    })
    .unwrap();
    let decoded_path = temp.path().join("decoded.fastq");
    qz_lib::compression::decompress(&qz_lib::cli::DecompressConfig {
        input: archive,
        output: vec![decoded_path.clone()],
        working_dir: temp.path().to_path_buf(),
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
        force: false,
    })
    .unwrap();
    let decoded = std::fs::read_to_string(decoded_path).unwrap();
    for i in 0..200 {
        assert!(
            decoded.contains(&format!("@r{i}\nACGTACGTACGTACGT\n")),
            "missing record {i}"
        );
    }
}

/// A FASTQ record with a >63 MiB sequence must be rejected at compress time:
/// the bounded streaming decoder's per-block buffer is capped at
/// `BSC_MAX_BLOCK_LEN = 64 MiB + 4 KiB`, and the encoder enforces the cap
/// in `records_to_streams` (chunked path) and the migrated `codecs.rs`
/// helpers (in-memory path).
#[test]
fn test_bounded_streaming_oversized_record_rejected() {
    let huge_seq = vec![b'A'; 64 * 1024 * 1024];
    let huge_qual = vec![b'I'; 64 * 1024 * 1024];
    let data = format!(
        "@big\n{}\n+\n{}\n",
        std::str::from_utf8(&huge_seq).unwrap(),
        std::str::from_utf8(&huge_qual).unwrap()
    );
    let temp = tempfile::TempDir::new().unwrap();
    let input = temp.path().join("in.fastq");
    let archive = temp.path().join("out.qz");
    std::fs::write(&input, &data).unwrap();
    let err = qz_lib::compression::compress(&qz_lib::cli::CompressConfig {
        input: vec![input],
        output: archive,
        working_dir: temp.path().to_path_buf(),
        threads: 1,
        ..Default::default()
    })
    .unwrap_err();
    let msg = format!("{err:#}");
    assert!(
        msg.contains("63") || msg.contains("per-record cap"),
        "expected per-record cap error, got: {msg}"
    );
}

// ---------------------------------------------------------------------------
// Task 21: Memory + first-record latency tests (Linux-only)
// ---------------------------------------------------------------------------


#[test]
#[cfg(target_os = "linux")]
fn test_bounded_streaming_bounded_memory_independent_quality() {
    // Decompress a 100K-record archive in a SUBPROCESS and assert its peak RSS stays
    // well below an O(reads) materialize. Measured via `/usr/bin/time -v` Maximum RSS
    // on an isolated child — NOT in-process sticky VmHWM, which is contaminated by
    // other parallel tests' allocations and by the compress baseline (the streaming
    // compress engine deliberately keeps that baseline low, which an in-process
    // before/after delta would mistake for a decompress regression).
    let qz = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/release/qz");
    if !qz.exists() {
        eprintln!(
            "skipping: {} not built (run `cargo build --release -p qz-cli`)",
            qz.display()
        );
        return;
    }
    let temp = tempfile::TempDir::new().unwrap();
    let p = temp.path();
    let mut data = String::new();
    for i in 0..100_000 {
        data.push_str(&format!(
            "@r{i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    let input = p.join("in.fastq");
    std::fs::write(&input, &data).unwrap();
    let archive = p.join("out.qz");
    let st = std::process::Command::new(&qz)
        .args(["compress", "-i"])
        .arg(&input)
        .arg("-o")
        .arg(&archive)
        .arg("-w")
        .arg(p)
        .args(["-t", "1", "-f"])
        .status()
        .unwrap();
    assert!(st.success(), "compress failed");

    let decoded = p.join("decoded.fastq");
    let timed = std::process::Command::new("/usr/bin/time")
        .arg("-v")
        .arg(&qz)
        .args(["decompress", "-i"])
        .arg(&archive)
        .arg("-o")
        .arg(&decoded)
        .arg("-w")
        .arg(p)
        .args(["-t", "1", "-f"])
        .output()
        .unwrap();
    assert!(timed.status.success(), "decompress failed");
    let stderr = String::from_utf8_lossy(&timed.stderr);
    let max_rss_kb: u64 = stderr
        .lines()
        .find(|l| l.contains("Maximum resident set size"))
        .and_then(|l| l.rsplit(':').next())
        .and_then(|v| v.trim().parse().ok())
        .expect("parse Max RSS from /usr/bin/time -v");
    // 100K reads of ~32-byte sequences ~= 6 MB input; the bounded decode path should
    // stay well under 500 MB absolute even with the process baseline + mimalloc pool.
    assert!(
        max_rss_kb < 500_000,
        "decompress peak RSS = {max_rss_kb} kB exceeds 500 MB"
    );
    assert_eq!(
        std::fs::read(&decoded).unwrap(),
        data.as_bytes(),
        "decoded content mismatch"
    );
}

#[test]
#[cfg(target_os = "linux")]
fn test_bounded_streaming_first_record_latency() {
    // Smoke test: decompress to a custom writer that records when the first
    // FASTQ record (4 lines) arrives, asserting it shows up in bounded
    // wall-clock time rather than the decoder hanging or doing unbounded work
    // first. This is a coarse latency ceiling, not a streaming-vs-buffering
    // discriminator — the no-buffering guarantee is covered by
    // `test_bounded_streaming_bounded_memory_*`. The threshold is deliberately
    // very generous: the suite now includes several heavy big-block ultra tests
    // that run concurrently, so first-record wall-clock can balloon under
    // full-parallel CPU contention (5 s was no longer enough). 20 s still
    // catches a real regression — the decoder hanging or doing unbounded work
    // before emitting the first record would be far worse than 20 s.
    let mut data = String::new();
    for i in 0..100_000 {
        data.push_str(&format!("@r{i}\nACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFF\n"));
    }
    let temp = tempfile::TempDir::new().unwrap();
    let input = temp.path().join("in.fastq");
    let archive = temp.path().join("out.qz");
    std::fs::write(&input, &data).unwrap();
    qz_lib::compression::compress(&qz_lib::cli::CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: temp.path().to_path_buf(),
        threads: 1,
        ..Default::default()
    })
    .unwrap();

    use std::time::Instant;
    struct LatencyWriter {
        start: Instant,
        first_record_at: Option<u128>,
        newlines: usize,
    }
    impl std::io::Write for LatencyWriter {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.newlines += buf.iter().filter(|&&b| b == b'\n').count();
            if self.first_record_at.is_none() && self.newlines >= 4 {
                self.first_record_at = Some(self.start.elapsed().as_millis());
            }
            Ok(buf.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut writer = LatencyWriter {
        start: Instant::now(),
        first_record_at: None,
        newlines: 0,
    };
    qz_lib::compression::decompress_to_writer(&archive, &mut writer).unwrap();
    let latency_ms = writer.first_record_at.expect("never saw first record");
    assert!(
        latency_ms < 20000,
        "first-record latency = {latency_ms} ms > 20000 ms"
    );
}

/// Regression: when chunk 0 has a constant read length but a later chunk
/// dispatched inside the pipeline priming window has a *different* length,
/// the const-length framing optimization (which omits per-record length
/// prefixes) must not silently corrupt the archive. The steady-state loop
/// validates this and bails; the priming loop must do the same.
///
/// Contract under test: compression either refuses (clean Err) or produces a
/// losslessly-decodable archive — never a successfully-written corrupt one.
#[test]
fn test_chunked_priming_window_const_length_mismatch_is_not_silently_corrupt() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    // Chunk 0: 4 reads of 10 bp (constant). Chunk 1: 4 reads of 8 bp (constant
    // but different). With chunk_records=4 and compress_window=2, chunk 1 is
    // dispatched inside the priming loop, bypassing steady-state validation.
    let mut data = String::new();
    for i in 0..4 {
        data.push_str(&format!("@c0r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n"));
    }
    for i in 0..4 {
        data.push_str(&format!("@c1r{i}\nACGTACGT\n+\nIIIIIIII\n"));
    }

    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, &data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let mut compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    compress_args.advanced.chunk_records = 4;
    compress_args.advanced.compress_window = 2;

    match qz_lib::compression::compress(&compress_args) {
        // Acceptable: refused to build a const-length archive it can't honor.
        Err(_) => {}
        // If it claims success, the archive MUST decode losslessly.
        Ok(()) => {
            let output_fastq = temp_path.join("output.fastq");
            qz_lib::compression::decompress(&decompress_args(
                archive_path,
                output_fastq.clone(),
                temp_path,
            ))
            .expect("archive written by compress must be decodable");
            let decompressed = fs::read_to_string(&output_fastq).unwrap();
            assert_eq!(
                data, decompressed,
                "chunked priming-window roundtrip must be lossless (silent corruption)"
            );
        }
    }
}

/// Regression: FASTQ headers contain arbitrary bytes, but the header
/// compressors were String-based and mangled non-UTF-8 bytes to U+FFFD. A
/// header with a non-UTF-8 byte must now round-trip byte-exact under the
/// default config (columnar headers) and under ultra mode.
#[test]
fn test_non_utf8_header_roundtrips_byte_exact_default_and_ultra() {
    fn check(config_fn: impl FnOnce(&mut CompressConfig), label: &str) {
        let temp_dir = TempDir::new().unwrap();
        let temp_path = temp_dir.path().to_path_buf();

        let input_fastq = temp_path.join("input.fastq");
        let mut data: Vec<u8> = Vec::new();
        data.extend_from_slice(b"@read");
        data.push(0xFF); // non-UTF-8 byte in the header ID
        data.extend_from_slice(b"1 some:comment\nACGTACGT\n+\nIIIIIIII\n");
        data.extend_from_slice(b"@read2\nTGCATGCA\n+\nHHHHHHHH\n");
        fs::write(&input_fastq, &data).unwrap();

        let archive_path = temp_path.join("test.qz");
        let mut compress_args = CompressConfig {
            input: vec![input_fastq.clone()],
            output: archive_path.clone(),
            working_dir: temp_path.clone(),
            threads: 1,
            ..CompressConfig::default()
        };
        config_fn(&mut compress_args);

        qz_lib::compression::compress(&compress_args).unwrap();

        let output_fastq = temp_path.join("output.fastq");
        qz_lib::compression::decompress(&decompress_args(
            archive_path,
            output_fastq.clone(),
            temp_path,
        ))
        .unwrap();

        let original = fs::read(&input_fastq).unwrap();
        let decompressed = fs::read(&output_fastq).unwrap();
        assert_eq!(
            original, decompressed,
            "non-UTF-8 header must round-trip byte-exact ({label})"
        );
    }

    check(|_| {}, "default");
    check(|c| c.ultra = Some(1), "ultra");
}

/// `--fast` / `bsc_static` uses libbsc's static QLFC coder for the
/// header/sequence streams (faster encode+decode, small ratio cost). Archives
/// produced this way must still decode byte-exact, with and without N bases.
#[test]
fn test_bsc_static_fast_mode_roundtrips() {
    roundtrip_exact(MULTI_READ_DATA, |c| c.advanced.bsc_static = true);

    let with_n = "\
@r1\nACGTNNNNACGTACGT\n+\nIIIIIIIIIIIIIIII\n\
@r2\nNNNNNNNNNNNNNNNN\n+\nHHHHHHHHHHHHHHHH\n\
@r3\nACGTACGTACGTACGT\n+\nBBBBBBBBBBBBBBBB\n";
    roundtrip_exact(with_n, |c| c.advanced.bsc_static = true);
}

// ---------------------------------------------------------------------------
// Task 22: Property test for random encoder/decoder roundtrip
// ---------------------------------------------------------------------------

proptest::proptest! {
    #[test]
    fn prop_bounded_streaming_roundtrip(
        num_reads in 0usize..500,
        read_length in 1usize..200,
    ) {
        // Skip empty-input case (compress bails — covered by a separate test).
        if num_reads == 0 { return Ok(()); }

        let mut data = String::new();
        for i in 0..num_reads {
            let seq: String = (0..read_length).map(|j| match (i + j) % 4 { 0 => 'A', 1 => 'C', 2 => 'G', _ => 'T' }).collect();
            let qual: String = std::iter::repeat_n('I', read_length).collect();
            data.push_str(&format!("@r{i}\n{seq}\n+\n{qual}\n"));
        }
        let temp = tempfile::TempDir::new().unwrap();
        let input = temp.path().join("in.fastq");
        let archive = temp.path().join("out.qz");
        let decoded = temp.path().join("decoded.fastq");
        std::fs::write(&input, &data).unwrap();
        qz_lib::compression::compress(&qz_lib::cli::CompressConfig {
            input: vec![input],
            output: archive.clone(),
            working_dir: temp.path().to_path_buf(),
            threads: 1,
            ..Default::default()
        }).unwrap();
        qz_lib::compression::decompress(&qz_lib::cli::DecompressConfig {
            input: archive,
            output: vec![decoded.clone()],
            working_dir: temp.path().to_path_buf(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 6,
            force: false,
        }).unwrap();
        let got = std::fs::read_to_string(decoded).unwrap();
        proptest::prop_assert_eq!(got, data);
    }
}

#[test]
fn ultra_decompress_peak_rss_is_bounded() {
    // The qz binary lives in the workspace target dir (qz-cli crate).
    let qz = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/release/qz");
    if !qz.exists() {
        eprintln!(
            "skipping: {} not built (run `cargo build --release -p qz-cli`)",
            qz.display()
        );
        return;
    }
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    // 130k reads x 600 bp => one ~78 MB sequence block (> 64 MiB). Decoding it
    // via the bounded path needs ~block-size RAM, NOT all-records RAM.
    let bases = [b'A', b'C', b'G', b'T'];
    let mut data = Vec::new();
    for i in 0..130_000usize {
        let seq: Vec<u8> = (0..600).map(|j| bases[(i + j) % 4]).collect();
        data.extend_from_slice(format!("@r{i}\n").as_bytes());
        data.extend_from_slice(&seq);
        data.extend_from_slice(b"\n+\n");
        data.extend(std::iter::repeat_n(b'I', 600));
        data.push(b'\n');
    }
    let input = p.join("in.fastq");
    std::fs::write(&input, &data).unwrap();
    let archive = p.join("a.qz");
    let st = std::process::Command::new(&qz)
        .args(["compress", "-i"])
        .arg(&input)
        .arg("-o")
        .arg(&archive)
        .arg("-w")
        .arg(p)
        .args(["-t", "4", "-f", "--ultra", "1"])
        .status()
        .unwrap();
    assert!(st.success(), "ultra compress failed");

    let out = p.join("out.fastq");
    let timed = std::process::Command::new("/usr/bin/time")
        .arg("-v")
        .arg(&qz)
        .args(["decompress", "-i"])
        .arg(&archive)
        .arg("-o")
        .arg(&out)
        .arg("-w")
        .arg(p)
        .args(["-t", "4", "-f"])
        .output()
        .unwrap();
    assert!(timed.status.success(), "ultra decompress failed");
    let stderr = String::from_utf8_lossy(&timed.stderr);
    let max_rss_kb: u64 = stderr
        .lines()
        .find(|l| l.contains("Maximum resident set size"))
        .and_then(|l| l.rsplit(':').next())
        .and_then(|v| v.trim().parse().ok())
        .expect("parse Max RSS from /usr/bin/time -v");
    // Bounded path: a single 78 MB block needs ~block-size RAM, not all-records.
    assert!(
        max_rss_kb < 2_000_000,
        "decompress peak RSS = {max_rss_kb} kB exceeds 2 GB"
    );
    assert_eq!(
        std::fs::read(&out).unwrap(),
        data,
        "ultra roundtrip byte-exact"
    );
}

#[test]
fn ultra_level_validation() {
    use qz_lib::compression;
    let dir = TempDir::new().unwrap();
    let p = dir.path().to_path_buf();
    let input = p.join("in.fastq");
    fs::write(&input, MULTI_READ_DATA).unwrap();
    let mk = |lvl: Option<u8>| CompressConfig {
        input: vec![input.clone()],
        output: p.join("a.qz"),
        working_dir: p.clone(),
        threads: 1,
        ultra: lvl,
        force: true,
        ..CompressConfig::default()
    };
    assert!(
        compression::compress(&mk(Some(4))).is_err(),
        "level 4 must error"
    );
    assert!(
        compression::compress(&mk(Some(0))).is_ok(),
        "0 = auto must work"
    );
    assert!(
        compression::compress(&mk(Some(2))).is_ok(),
        "level 2 must work"
    );
}

#[test]
fn legacy_v4_archive_is_rejected() {
    // Any non-v5 (legacy v4 stream-major) archive is rejected outright now —
    // qz only reads v5 chunk-major. Hand-build a minimal v4 (version byte 4)
    // header and assert the recompress-hint rejection.
    let dir = TempDir::new().unwrap();
    let archive = dir.path().join("old.qz");
    // Minimal v2 header: magic(2) + version(1) + reserved(1) + header_size u32 LE
    // + a 52-byte fixed body (the old v4 body size).
    let mut bytes = vec![b'Q', b'Z']; // ARCHIVE_MAGIC
    bytes.push(4); // legacy v4 version byte
    bytes.push(0); // reserved
    bytes.extend_from_slice(&60u32.to_le_bytes()); // header_size
    let body = [0u8; 52];
    bytes.extend_from_slice(&body);
    std::fs::write(&archive, &bytes).unwrap();
    let out = dir.path().join("out.fastq");
    let err =
        qz_lib::compression::decompress(&decompress_args(archive, out, dir.path().to_path_buf()))
            .unwrap_err();
    assert!(
        err.to_string().contains("no longer supported"),
        "got: {err}"
    );
}

#[test]
fn ultra_const_length_no_hint_byte_roundtrips() {
    // 200k reads, fixed 150 bp (const-length, no hint byte). Single chunk at L1.
    let dir = TempDir::new().unwrap();
    let p = dir.path().to_path_buf();
    let bases = [b'A', b'C', b'G', b'T'];
    let mut data = Vec::new();
    for i in 0..200_000usize {
        let seq: Vec<u8> = (0..150).map(|j| bases[(i + j) % 4]).collect();
        data.extend_from_slice(format!("@r{i}\n").as_bytes());
        data.extend_from_slice(&seq);
        data.extend_from_slice(b"\n+\n");
        data.extend(std::iter::repeat_n(b'I', 150));
        data.push(b'\n');
    }
    let input = p.join("in.fastq");
    fs::write(&input, &data).unwrap();
    let archive = p.join("a.qz");
    let cfg = CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: p.clone(),
        threads: 4,
        ultra: Some(1),
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&cfg).unwrap();
    // Sanity: ultra emits encoding 10.
    assert_eq!(fs::read(&archive).unwrap()[8], 10);
    let out = p.join("out.fastq");
    qz_lib::compression::decompress(&decompress_args(archive, out.clone(), p)).unwrap();
    assert_eq!(
        fs::read(&out).unwrap(),
        data,
        "const-length ultra must roundtrip byte-exact"
    );
}

// ---------------------------------------------------------------------------
// Per-platform generic header codec tests (Task 9)
// ---------------------------------------------------------------------------

/// Build a FASTQ string from header lines, with a fixed 8bp seq/qual each.
fn fastq_from_headers(headers: &[String]) -> String {
    let mut s = String::new();
    for h in headers {
        s.push_str(h);
        s.push_str("\nACGTACGT\n+\nIIIIIIII\n");
    }
    s
}

fn platform_headers() -> Vec<(&'static str, Vec<String>)> {
    let n = 3000;
    vec![
        (
            "nanopore",
            (1..=n)
                .map(|i| {
                    format!(
                        "@d3{:04x}9a read={i} ch={} runid=abcd start_time=2019-09-12T12:00:0{}Z",
                        i % 4096,
                        i % 512,
                        i % 10
                    )
                })
                .collect(),
        ),
        (
            "pacbio",
            (1..=n)
                .map(|i| format!("@m54238_180628/{}/0_{}", 4194300 + i, 900 + i))
                .collect(),
        ),
        (
            "bgi",
            (1..=n)
                .map(|i| format!("@V300047L1C001R{:03}_{}/1", i % 1000, i))
                .collect(),
        ),
        (
            "older_illumina",
            (1..=n)
                .map(|i| format!("@HWUSI-EAS100R:6:73:{}:{}#0/1", 900 + i, 1900 + i))
                .collect(),
        ),
        ("numbered", (1..=n).map(|i| format!("@read{i}")).collect()),
    ]
}

#[test]
fn generic_header_codec_roundtrips_all_platforms_and_emits_0x05() {
    for (name, headers) in platform_headers() {
        let refs: Vec<&str> = headers.iter().map(String::as_str).collect();
        let blob = qz_lib::compression::header_col::compress_headers_columnar(&refs).unwrap();
        assert_eq!(
            blob[0], 0x05,
            "{name}: expected generic tier (0x05), got 0x{:02x}",
            blob[0]
        );
        let dec = qz_lib::compression::header_col::decompress_headers_columnar(&blob, refs.len())
            .unwrap();
        for (i, h) in refs.iter().enumerate() {
            assert_eq!(
                dec[i].as_slice(),
                h.as_bytes(),
                "{name}: header {i} mismatch"
            );
        }
    }
}

#[test]
fn generic_header_codec_end_to_end_byte_identical_all_platforms() {
    for (name, headers) in platform_headers() {
        let input = fastq_from_headers(&headers);
        std::panic::catch_unwind(|| roundtrip_exact(&input, |_| {}))
            .unwrap_or_else(|_| panic!("end-to-end roundtrip failed for platform {name}"));
    }
}

// ---------------------------------------------------------------------------
// Task 1.10 + 1.11: chunk-major v5 archive coverage.
//
// The archive layout (see `old_reorder_local_encoding_is_unsupported`) is:
//   byte[0..2] = magic "QZ"
//   byte[2]    = ARCHIVE_VERSION  (4 = stream-major, 5 = chunk-major)
//   byte[3]    = reserved
//   byte[4..8] = header_size (LE u32)
//   byte[8]    = encoding_type
// v5 chunk-major archives also end with the trailing footer magic "QZFOOTR1".
// ---------------------------------------------------------------------------

// Task 1.10: ultra emits a v5 chunk-major archive and round-trips byte-exact.
#[test]
fn ultra_v5_roundtrips() {
    // Ultra level 1 over a sizeable, VARIED input. The reads are deliberately
    // non-degenerate (per-record sequence + quality variation): libbsc's OpenMP
    // multithreaded BWT (used only by ultra) has a separate, pre-existing
    // intermittent-corruption bug on tiny/degenerate blocks (e.g. tens of
    // thousands of byte-identical short reads) that predates the chunk-major
    // work and reproduces on the v4 ultra binary too — it is NOT a v5 issue.
    // Using varied data exercises the v5 ultra framing/roundtrip reliably while
    // steering clear of that orthogonal MT-BSC edge case.
    let dir = TempDir::new().unwrap();
    let temp_path = dir.path().to_path_buf();
    let input = temp_path.join("in.fastq");
    let mut s = String::new();
    let bases = [b'A', b'C', b'G', b'T'];
    for i in 0..50_000u32 {
        let seq: String = (0..60u32)
            .map(|j| {
                bases[(((i.wrapping_mul(7)).wrapping_add(j.wrapping_mul(3))) % 4) as usize] as char
            })
            .collect();
        let qual: String = (0..60u32)
            .map(|j| {
                (b'!' + (((i.wrapping_mul(13)).wrapping_add(j.wrapping_mul(5))) % 40) as u8) as char
            })
            .collect();
        s.push_str(&format!("@read{i} lane:1:{i}\n{seq}\n+\n{qual}\n"));
    }
    fs::write(&input, &s).unwrap();

    let archive = temp_path.join("out.qz");
    let compress_args = CompressConfig {
        input: vec![input.clone()],
        output: archive.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ultra: Some(1),
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    // Assert v5 on disk, with the trailing footer magic.
    let bytes = fs::read(&archive).unwrap();
    assert_eq!(bytes[2], 5, "ultra must emit chunk-major v5");
    assert_eq!(
        &bytes[bytes.len() - 8..],
        b"QZFOOTR1",
        "missing v5 footer magic"
    );

    // Roundtrip byte-exact.
    let out = temp_path.join("rt.fastq");
    qz_lib::compression::decompress(&decompress_args(archive, out.clone(), temp_path)).unwrap();
    assert_eq!(fs::read(&out).unwrap(), s.into_bytes());
}

/// Shared small-but-multi-record FASTQ body for the v5 semantics matrix.
/// Reads are >= 21bp so syncmer-based sequence hints are exercised.
fn v5_matrix_fastq() -> String {
    let mut s = String::new();
    for i in 0..2_000 {
        s.push_str(&format!(
            "@read{i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFFGGAABBCCDDEEFFGG\n"
        ));
    }
    s
}

/// Compress `data` with `config_fn`, assert the on-disk version byte == `want_version`
/// (and, for v5, the trailing footer magic), then assert a byte-exact roundtrip.
/// Returns the raw archive bytes so callers can make extra assertions (e.g. encoding_type).
fn v5_roundtrip_checked(
    data: &str,
    want_version: u8,
    config_fn: impl FnOnce(&mut CompressConfig),
) -> Vec<u8> {
    let dir = TempDir::new().unwrap();
    let temp_path = dir.path().to_path_buf();
    let input = temp_path.join("in.fastq");
    fs::write(&input, data).unwrap();

    let archive = temp_path.join("out.qz");
    let mut compress_args = CompressConfig {
        input: vec![input.clone()],
        output: archive.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    };
    config_fn(&mut compress_args);
    qz_lib::compression::compress(&compress_args).unwrap();

    let bytes = fs::read(&archive).unwrap();
    assert_eq!(
        bytes[2], want_version,
        "archive version byte mismatch (want {want_version}, got {})",
        bytes[2]
    );
    if want_version == 5 {
        assert_eq!(
            &bytes[bytes.len() - 8..],
            b"QZFOOTR1",
            "v5 archive must end with footer magic"
        );
    }

    let out = temp_path.join("rt.fastq");
    qz_lib::compression::decompress(&decompress_args(archive, out.clone(), temp_path)).unwrap();
    assert_eq!(
        fs::read(&out).unwrap(),
        data.as_bytes(),
        "v{want_version} roundtrip must be byte-exact"
    );
    bytes
}

// Task 1.11: semantics matrix + quality-codec roundtrips.

#[test]
fn v5_raw_with_hints_roundtrips() {
    // sequence_hints -> RawWithHints (encoding_type 4) -> chunk-major v5.
    let bytes = v5_roundtrip_checked(&v5_matrix_fastq(), 5, |c| {
        c.advanced.sequence_hints = true;
    });
    assert_eq!(
        bytes[8], 4,
        "sequence_hints must produce encoding_type 4 (RawWithHints)"
    );
}

#[test]
fn v5_rc_canon_roundtrips() {
    // rc_canon -> RcCanon (encoding_type 6) -> chunk-major v5.
    let bytes = v5_roundtrip_checked(&v5_matrix_fastq(), 5, |c| {
        c.advanced.rc_canon = true;
    });
    assert_eq!(
        bytes[8], 6,
        "rc_canon must produce encoding_type 6 (RcCanon)"
    );
}

#[test]
fn v5_no_quality_roundtrips() {
    // no_quality discards qualities; the FASTQ shape is preserved (empty quality
    // lines) and the archive is chunk-major v5.
    //
    // When no_quality is set, the decompressed FASTQ has the '+' line followed by
    // an empty quality line, so it is NOT byte-identical to the input. Assert the
    // version + footer here, then check the structural roundtrip separately.
    let dir = TempDir::new().unwrap();
    let temp_path = dir.path().to_path_buf();
    let input = temp_path.join("in.fastq");
    let data = v5_matrix_fastq();
    fs::write(&input, &data).unwrap();

    let archive = temp_path.join("out.qz");
    qz_lib::compression::compress(&CompressConfig {
        input: vec![input.clone()],
        output: archive.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        no_quality: true,
        ..CompressConfig::default()
    })
    .unwrap();

    let bytes = fs::read(&archive).unwrap();
    assert_eq!(bytes[2], 5, "no_quality archive must be chunk-major v5");
    assert_eq!(&bytes[bytes.len() - 8..], b"QZFOOTR1");

    let out = temp_path.join("rt.fastq");
    qz_lib::compression::decompress(&decompress_args(archive, out.clone(), temp_path)).unwrap();
    let decoded = fs::read_to_string(&out).unwrap();
    // Headers + sequences are preserved; quality lines are dropped (empty).
    for i in 0..2_000 {
        assert!(
            decoded.contains(&format!("@read{i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n")),
            "no_quality roundtrip missing record {i}"
        );
    }
}

#[test]
fn v5_fasta_roundtrips() {
    // FASTA input must round-trip byte-exact as FASTA and emit a chunk-major v5
    // archive (FASTA reuses the bounded-streamable Raw encoding).
    let data = {
        let mut s = String::new();
        for i in 0..2_000 {
            s.push_str(&format!(
                ">seq{i} desc {i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n"
            ));
        }
        s
    };
    v5_roundtrip_checked(&data, 5, |c| {
        c.fasta = true;
    });
}

#[test]
fn v5_fqzcomp_quality_roundtrips() {
    // Fqzcomp is the post-streaming-fqz default quality codec.
    v5_roundtrip_checked(&v5_matrix_fastq(), 5, |c| {
        c.advanced.quality_compressor = QualityCompressor::Fqzcomp;
    });
}

#[test]
fn v5_bsc_quality_roundtrips() {
    v5_roundtrip_checked(&v5_matrix_fastq(), 5, |c| {
        c.advanced.quality_compressor = QualityCompressor::Bsc;
    });
}

#[test]
fn v5_decompress_to_writer_roundtrips() {
    // Exercise the public `decompress_to_writer` entry (the path Python / pipe
    // consumers use) on a v5 archive; the written FASTQ bytes must equal the input.
    let dir = TempDir::new().unwrap();
    let temp_path = dir.path().to_path_buf();
    let input = temp_path.join("in.fastq");
    let data = v5_matrix_fastq();
    fs::write(&input, &data).unwrap();

    let archive = temp_path.join("out.qz");
    qz_lib::compression::compress(&CompressConfig {
        input: vec![input.clone()],
        output: archive.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        ..CompressConfig::default()
    })
    .unwrap();

    let bytes = fs::read(&archive).unwrap();
    assert_eq!(bytes[2], 5, "default archive must be chunk-major v5");

    let mut writer: Vec<u8> = Vec::new();
    qz_lib::compression::decompress_to_writer(&archive, &mut writer).unwrap();
    assert_eq!(
        writer,
        data.into_bytes(),
        "decompress_to_writer must be byte-exact"
    );
}

// ---------------------------------------------------------------------------
// Task 1.12: adversarial footer / integrity tests.
//
// The v5 archive ends with a 20-byte locator:
//   [footer_len: u64 LE][footer_crc32: u32 LE][magic "QZFOOTR1": 8B]
// so for an archive of length `len`:
//   footer_len   = bytes[len-20 .. len-12]
//   footer_crc32 = bytes[len-12 .. len-8]
//   magic        = bytes[len-8  .. len]
//   footer_start = len - 20 - footer_len
// The footer BODY occupies [footer_start .. len-20] and is:
//   num_reads(8) | num_chunks(4) | n_entries(4)   (16-byte prefix)
//   then n_entries × 31-byte ChunkDirEntry:
//     chunk_index u32(4) | role u8(1) | mate u8(1) | codec u8(1) | offset u64(8)
//     | length u64(8) | record_count u64(8)
// Entry i begins at  footer_start + 16 + i*31.  Field offsets within an entry:
//   chunk_index +0 | role +4 | mate +5 | codec +6 | offset +7 | length +15 | record_count +23 (8 bytes)
//
// The decoder (build_indices_from_footer) validates: locator magic, footer_len
// cap / underflow, footer CRC32, per-entry frame span + offset/length overflow,
// inline-v4-header cross-check (record_count / block_len), global non-overlap,
// and the role/chunk pairing contract (incl. rejecting Nmask frames). Tests that
// mutate the footer BODY and need to reach a DEEPER check than the CRC gate must
// recompute the CRC32 over the mutated body and rewrite [len-12 .. len-8].
// `chunk_directory` is pub(crate), so we recompute inline with `flate2::Crc`.
// ---------------------------------------------------------------------------

/// Compress ~300 records with a small chunk size (chunk_records ≈ 100 → ≥3
/// chunks → multiple directory entries per role) into a v5 archive. Asserts the
/// archive is chunk-major v5 and returns its path. The returned `TempDir` keeps
/// the directory alive for the caller.
fn make_v5(dir: &std::path::Path) -> std::path::PathBuf {
    let input = dir.join("mk.fastq");
    let mut s = String::new();
    for i in 0..300 {
        s.push_str(&format!(
            "@read{i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFFGGAABBCCDDEEFFGG\n"
        ));
    }
    fs::write(&input, &s).unwrap();

    let archive = dir.join("mk.qz");
    let mut compress_args = CompressConfig {
        input: vec![input],
        output: archive.clone(),
        working_dir: dir.to_path_buf(),
        threads: 1,
        ..CompressConfig::default()
    };
    // Small chunks → ≥3 chunks → multiple directory entries per role.
    compress_args.advanced.chunk_records = 100;
    qz_lib::compression::compress(&compress_args).unwrap();

    let bytes = fs::read(&archive).unwrap();
    assert_eq!(bytes[2], 5, "make_v5 must produce a chunk-major v5 archive");
    assert_eq!(
        &bytes[bytes.len() - 8..],
        b"QZFOOTR1",
        "make_v5 archive must end with the footer magic"
    );
    archive
}

/// Read the locator's `footer_len` (u64 LE at [len-20 .. len-12]) and derive the
/// footer-body start offset.
fn v5_footer_start(bytes: &[u8]) -> usize {
    let len = bytes.len();
    let footer_len = u64::from_le_bytes(bytes[len - 20..len - 12].try_into().unwrap()) as usize;
    len - 20 - footer_len
}

/// Recompute the footer CRC32 over the mutated body [footer_start .. len-20] and
/// rewrite the 4-byte CRC field at [len-12 .. len-8] so the body passes the CRC
/// gate and the decoder reaches the deeper validation.
fn v5_refresh_footer_crc(bytes: &mut [u8]) {
    let len = bytes.len();
    let footer_start = v5_footer_start(bytes);
    let mut c = flate2::Crc::new();
    c.update(&bytes[footer_start..len - 20]);
    let crc = c.sum();
    bytes[len - 12..len - 8].copy_from_slice(&crc.to_le_bytes());
}

/// Byte offset of entry `i` inside the footer body, given the body start.
fn v5_entry_offset(footer_start: usize, i: usize) -> usize {
    footer_start + 16 + i * 31
}

/// Write `bytes` to a fresh path in `dir` and attempt to `decompress` it,
/// returning the result (the caller asserts it is an Err).
fn v5_decompress_mutated(dir: &std::path::Path, bytes: &[u8]) -> anyhow::Result<()> {
    let bad = dir.join("bad.qz");
    fs::write(&bad, bytes).unwrap();
    let out = dir.join("bad_out.fastq");
    qz_lib::compression::decompress(&decompress_args(bad, out, dir.to_path_buf()))
}

#[test]
fn v5_truncated_locator_rejected() {
    // Drop the last 8 bytes (the magic) → the locator is short/missing its magic.
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let bytes = fs::read(&archive).unwrap();
    let truncated = &bytes[..bytes.len() - 8];
    let res = v5_decompress_mutated(dir.path(), truncated);
    assert!(
        res.is_err(),
        "truncated locator must be rejected, got Ok: {res:?}"
    );
}

#[test]
fn v5_bad_footer_magic_rejected() {
    // Flip a byte in the trailing magic [len-8 .. len] → "bad footer magic".
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    let len = bytes.len();
    bytes[len - 4] ^= 0xFF; // inside the 8-byte magic.
    let res = v5_decompress_mutated(dir.path(), &bytes);
    let err = res.unwrap_err().to_string();
    assert!(
        err.contains("bad footer magic"),
        "expected bad-footer-magic rejection, got: {err}"
    );
}

#[test]
fn v5_corrupt_footer_crc_rejected() {
    // Flip a byte in the footer body and do NOT refresh the CRC → the footer CRC
    // gate must fire ("footer CRC mismatch").
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    let footer_start = v5_footer_start(&bytes);
    // Flip a byte inside the body (num_reads region of the prefix).
    bytes[footer_start] ^= 0xFF;
    let res = v5_decompress_mutated(dir.path(), &bytes);
    let err = res.unwrap_err().to_string();
    assert!(
        err.contains("footer CRC mismatch"),
        "expected footer-CRC-mismatch rejection, got: {err}"
    );
}

#[test]
fn v5_footer_record_count_mismatch_rejected() {
    // Change an entry's record_count in the footer (REFRESH the CRC so it passes
    // the CRC gate) → the inline-v4-header cross-check must reject it as
    // "inline record_count != footer".
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    let footer_start = v5_footer_start(&bytes);
    // record_count of entry 0 lives at +23 within the 31-byte entry (8 bytes).
    let rc_off = v5_entry_offset(footer_start, 0) + 23;
    let orig = u64::from_le_bytes(bytes[rc_off..rc_off + 8].try_into().unwrap());
    bytes[rc_off..rc_off + 8].copy_from_slice(&orig.wrapping_add(1).to_le_bytes());
    v5_refresh_footer_crc(&mut bytes);
    let res = v5_decompress_mutated(dir.path(), &bytes);
    let err = res.unwrap_err().to_string();
    assert!(
        err.contains("record_count") && err.contains("!= footer"),
        "expected inline-vs-footer record_count mismatch, got: {err}"
    );
}

#[test]
fn v5_block_crc_corruption_rejected() {
    // Flip a byte in a block PAYLOAD in the data region (before footer_start).
    // The inline 12-byte v4 header [block_len|record_count|crc32] is left intact,
    // so the entry passes the span / inline cross-check / overlap gates and the
    // corruption surfaces only as a per-block CRC mismatch during decode.
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    let header_end = u32::from_le_bytes(bytes[4..8].try_into().unwrap()) as usize;
    // First block frame starts at header_end: 12-byte header then payload. Flip a
    // byte inside the payload of the first block.
    let payload_byte = header_end + 12;
    assert!(
        payload_byte < v5_footer_start(&bytes),
        "payload offset sanity"
    );
    bytes[payload_byte] ^= 0xFF;
    let res = v5_decompress_mutated(dir.path(), &bytes);
    assert!(
        res.is_err(),
        "block payload corruption must be rejected (per-block CRC), got Ok: {res:?}"
    );
    let err = res.unwrap_err().to_string();
    assert!(
        err.to_lowercase().contains("crc"),
        "expected a CRC-related error for block payload corruption, got: {err}"
    );
}

#[test]
fn v5_overlapping_spans_rejected() {
    // Make two frame spans overlap. The per-entry inline cross-check (which seeks
    // to each entry's `offset` and demands the inline v4 header's block_len /
    // record_count equal the footer's) runs BEFORE the global non-overlap gate,
    // so a span edit that breaks an entry's own inline match would trip the
    // inline check first. To isolate the overlap gate we duplicate entry 0 into
    // entry 1 (same offset, length, record_count): entry 1 then reads entry 0's
    // real inline header — which matches its (now identical) footer fields — so
    // both entries pass the inline cross-check, and the two identical spans
    // collide → "overlapping frame spans".
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    let footer_start = v5_footer_start(&bytes);
    let e0 = v5_entry_offset(footer_start, 0);
    let e1 = v5_entry_offset(footer_start, 1);
    // Copy entry 0's offset(+7,8) + length(+15,8) + record_count(+23,8) onto
    // entry 1, leaving entry 1's chunk_index/role/mate/codec (+0..+7) intact so the
    // role/chunk-pairing contract is not what trips first.
    let e0_off = bytes[e0 + 7..e0 + 15].to_vec();
    let e0_len = bytes[e0 + 15..e0 + 23].to_vec();
    let e0_rc = bytes[e0 + 23..e0 + 31].to_vec();
    bytes[e1 + 7..e1 + 15].copy_from_slice(&e0_off);
    bytes[e1 + 15..e1 + 23].copy_from_slice(&e0_len);
    bytes[e1 + 23..e1 + 31].copy_from_slice(&e0_rc);
    v5_refresh_footer_crc(&mut bytes);
    let res = v5_decompress_mutated(dir.path(), &bytes);
    let err = res.unwrap_err().to_string();
    assert!(
        err.contains("overlapping frame spans"),
        "expected the overlap gate to fire, got: {err}"
    );
}

#[test]
fn v5_offset_overflow_rejected() {
    // Set an entry's length to u64::MAX so offset+length overflows. Refresh CRC →
    // the overflow guard must fire ("offset+length overflow").
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    let footer_start = v5_footer_start(&bytes);
    // length of entry 0 lives at +15 within the 31-byte entry.
    let len_off = v5_entry_offset(footer_start, 0) + 15;
    bytes[len_off..len_off + 8].copy_from_slice(&u64::MAX.to_le_bytes());
    v5_refresh_footer_crc(&mut bytes);
    let res = v5_decompress_mutated(dir.path(), &bytes);
    let err = res.unwrap_err().to_string();
    assert!(
        err.contains("offset+length overflow"),
        "expected offset+length overflow rejection, got: {err}"
    );
}

#[test]
fn v5_role_corruption_rejected() {
    // Change an entry's role byte to Nmask (2), which single-end v5 never emits.
    // Refresh CRC → the decoder must reject the unexpected Nmask frame.
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    let footer_start = v5_footer_start(&bytes);
    // role byte of entry 0 lives at +4 within the entry.
    let role_off = v5_entry_offset(footer_start, 0) + 4;
    bytes[role_off] = 2; // StreamRole::Nmask
    v5_refresh_footer_crc(&mut bytes);
    let res = v5_decompress_mutated(dir.path(), &bytes);
    let err = res.unwrap_err().to_string();
    assert!(
        err.contains("Nmask")
            || err.contains("unknown role byte")
            || err.contains("record counts disagree")
            || err.contains("!= footer"),
        "expected a role-corruption rejection (Nmask / unknown role / count disagreement), got: {err}"
    );

    // Also exercise an entirely out-of-range role byte (not a defined variant) →
    // ChunkDirectory::parse must reject "unknown role byte".
    let mut bytes2 = fs::read(&archive).unwrap();
    let fs2 = v5_footer_start(&bytes2);
    let role_off2 = v5_entry_offset(fs2, 0) + 4;
    bytes2[role_off2] = 99; // not a valid StreamRole
    v5_refresh_footer_crc(&mut bytes2);
    let res2 = v5_decompress_mutated(dir.path(), &bytes2);
    let err2 = res2.unwrap_err().to_string();
    assert!(
        err2.contains("unknown role byte"),
        "expected unknown-role-byte rejection, got: {err2}"
    );
}

#[test]
fn v5_single_end_rejects_nonzero_mate() {
    // Single-end v5 never emits mate != 0. Set entry 0's mate byte (+5) to 1 and
    // refresh the footer CRC → the single-end validator must reject it.
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    let footer_start = v5_footer_start(&bytes);
    let mate_off = v5_entry_offset(footer_start, 0) + 5;
    bytes[mate_off] = 1;
    v5_refresh_footer_crc(&mut bytes);
    let res = v5_decompress_mutated(dir.path(), &bytes);
    let err = res.unwrap_err().to_string();
    assert!(
        err.contains("mate"),
        "expected a mate-rejection error, got: {err}"
    );
}

#[test]
fn v5_single_end_rejects_headerdelta_role() {
    // Single-end v5 never emits HeaderDelta (role 5). Set entry 0's role byte (+4)
    // to 5 and refresh the footer CRC → the single-end validator must reject it.
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    let footer_start = v5_footer_start(&bytes);
    let role_off = v5_entry_offset(footer_start, 0) + 4;
    bytes[role_off] = 5; // StreamRole::HeaderDelta
    v5_refresh_footer_crc(&mut bytes);
    let res = v5_decompress_mutated(dir.path(), &bytes);
    let err = res.unwrap_err().to_string();
    assert!(
        err.contains("HeaderDelta"),
        "expected a HeaderDelta-rejection error, got: {err}"
    );
}

#[test]
fn v5_single_end_rejects_paired_archive_type() {
    // archive_type byte is byte 3 of the front header. Flipping a single-end v5
    // archive's type to 1 (paired) now ROUTES it to the paired decoder (production
    // dispatch decodes archive_type=1 via decompress_paired_v5). The paired
    // directory validator then rejects it: a single-end footer carries mate-0
    // entries, which are illegal in the paired (archive_type=1) role/codec matrix.
    // So a mistyped archive still cannot be mis-decoded — it errors at validation.
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    bytes[3] = 1;
    let res = v5_decompress_mutated(dir.path(), &bytes);
    // Pin the specific gate: validate_paired_directory rejects the mate-0 entry.
    let err = res.unwrap_err().to_string();
    assert!(
        err.contains("illegal") && err.contains("mate 0"),
        "expected the paired directory validator to reject the mate-0 entry, got: {err}"
    );
}

#[test]
fn v5_overlarge_footer_len_rejected() {
    // Rewrite the locator's footer_len ([len-20 .. len-12]) to a huge value → the
    // MAX_FOOTER_BYTES cap (or the underflow guard) must reject it before the
    // footer is allocated/read. We do NOT refresh the CRC; the size gate fires
    // first.
    let dir = TempDir::new().unwrap();
    let archive = make_v5(dir.path());
    let mut bytes = fs::read(&archive).unwrap();
    let len = bytes.len();
    // A value far above MAX_FOOTER_BYTES (256 MiB) but below u64::MAX so it is the
    // cap check, not the underflow guard, that fires.
    bytes[len - 20..len - 12].copy_from_slice(&(1u64 << 40).to_le_bytes());
    let res = v5_decompress_mutated(dir.path(), &bytes);
    let err = res.unwrap_err().to_string();
    assert!(
        err.contains("footer_len") && err.contains("exceeds cap"),
        "expected footer_len cap rejection, got: {err}"
    );
}

// Regression: libbsc's OpenMP multithreaded BWT (used by ultra's big-block path)
// intermittently produced a CORRUPT, non-decodable block on SMALL, highly
// repetitive inputs (~2/3 of runs on 20000 byte-identical short reads → a ~320KB
// degenerate block). Root cause: the MT gate keyed on the configured block-size
// *target* (188 MiB for ultra), not the *actual* block size, so a tiny block was
// still routed through MT. The fix routes small actual blocks through
// single-thread BSC. Loop several times — every run must roundtrip losslessly.
#[test]
fn ultra_degenerate_small_input_roundtrips() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("degen.fastq");
    let mut s = String::new();
    for i in 0..20_000 {
        s.push_str(&format!("@r{i}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n"));
    }
    fs::write(&input, &s).unwrap();
    for iter in 0..8 {
        let archive = dir.path().join(format!("degen{iter}.qz"));
        let out = dir.path().join(format!("degen{iter}.out"));
        let cfg = CompressConfig {
            input: vec![input.clone()],
            output: archive.clone(),
            working_dir: dir.path().to_path_buf(),
            threads: 1,
            ultra: Some(1),
            ..CompressConfig::default()
        };
        qz_lib::compression::compress(&cfg).unwrap();
        qz_lib::compression::decompress(&decompress_args(
            archive,
            out.clone(),
            dir.path().to_path_buf(),
        ))
        .unwrap();
        assert_eq!(
            fs::read(&out).unwrap(),
            s.as_bytes(),
            "ultra degenerate roundtrip failed on iteration {iter}"
        );
    }
}
