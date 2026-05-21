use qz_lib::cli::{AdvancedOptions, CompressConfig, DecompressConfig, HeaderCompressor, QualityCompressor, QualityMode, ReorderMode, SequenceCompressor};
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
        archive_path, output_fastq.clone(), temp_path,
    )).unwrap();

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
fn decompress_args(input: std::path::PathBuf, output: std::path::PathBuf, working_dir: std::path::PathBuf) -> DecompressConfig {
    DecompressConfig {
        input,
        output: vec![output],
        working_dir,
        num_threads: 1,
        gzipped: false,
        gzip_level: 6,
    }
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
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

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
        advanced: AdvancedOptions {
            quality_compressor: QualityCompressor::Zstd,
            sequence_compressor: SequenceCompressor::Zstd,
            header_compressor: HeaderCompressor::Zstd,
            ..Default::default()
        },
        ..CompressConfig::default()
    };
    qz_lib::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    let content = fs::read_to_string(&output_fastq).unwrap();
    assert!(content.contains("IIIIIIII")); // Dummy quality
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
            assert!(archive_path.exists(), "Empty archive should exist after successful compress");
        }
        Err(e) => {
            let msg = e.to_string().to_lowercase();
            assert!(msg.contains("empty") || msg.contains("no reads") || msg.contains("0 reads"),
                "Empty file error should mention emptiness, got: {}", msg);
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
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

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
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

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
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

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
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());
    let content = fs::read_to_string(&output_fastq).unwrap();
    // Sequence must survive intact
    assert!(content.contains("ACGTACGT"), "sequence missing from illumina-binned output");
    // Header must survive
    assert!(content.contains("@read1"), "header missing from illumina-binned output");
    // Quality must be present (binned, not identical, but same length as sequence)
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 4, "expected 4 lines in output FASTQ");
    assert_eq!(lines[1].len(), lines[3].len(),
        "quality length ({}) must match sequence length ({})", lines[3].len(), lines[1].len());
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
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());
    let content = fs::read_to_string(&output_fastq).unwrap();
    // Sequence and header must survive intact
    assert!(content.contains("ACGTACGT"), "sequence missing from binary quality output");
    assert!(content.contains("@read1"), "header missing from binary quality output");
    // Quality must be present and same length as sequence
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 4, "expected 4 lines in output FASTQ");
    assert_eq!(lines[1].len(), lines[3].len(),
        "quality length ({}) must match sequence length ({})", lines[3].len(), lines[1].len());
}

// ========================================
// COMPRESSION LEVEL TESTS
// ========================================

#[test]
fn test_different_compression_levels() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n".repeat(100);
    let input_fastq = temp_path.join("test.fastq");
    fs::write(&input_fastq, &test_data).unwrap();

    let levels = vec![1, 3, 9, 19];
    let mut sizes = Vec::new();

    for level in levels {
        let archive_path = temp_path.join(format!("test_level_{}.qz", level));
        let compress_args = CompressConfig {
            input: vec![input_fastq.clone()],
            output: archive_path.clone(),
            working_dir: temp_path.clone(),
            threads: 1,
            advanced: AdvancedOptions {
                compression_level: level,
                ..Default::default()
            },
            ..CompressConfig::default()
        };

        qz_lib::compression::compress(&compress_args).unwrap();
        let size = fs::metadata(&archive_path).unwrap().len();
        sizes.push((level, size));
    }

    assert_eq!(sizes.len(), 4);

    for (level, size) in &sizes {
        println!("Level {}: {} bytes", level, size);
    }
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

    qz_lib::compression::decompress(&decompress_args(archive1, output1.clone(), temp_path.clone())).unwrap();
    qz_lib::compression::decompress(&decompress_args(archive2, output2.clone(), temp_path)).unwrap();

    // Both should produce identical output
    let content1 = fs::read_to_string(&output1).unwrap();
    let content2 = fs::read_to_string(&output2).unwrap();
    assert_eq!(content1, content2, "Multi-threaded compression should produce same output as single-threaded");
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
    let result = qz_lib::compression::decompress(&decompress_args(
        truncated_path, output_fastq, temp_path,
    ));
    assert!(result.is_err(), "Decompressing a truncated archive should fail");
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
    let result = qz_lib::compression::decompress(&decompress_args(
        corrupted_path, output_fastq, temp_path,
    ));
    assert!(result.is_err(), "Decompressing a corrupted archive should fail");
}

#[test]
fn test_decompress_zero_byte_file() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let empty_archive = temp_path.join("empty.qz");
    fs::write(&empty_archive, b"").unwrap();

    let output_fastq = temp_path.join("output.fastq");
    let result = qz_lib::compression::decompress(&decompress_args(
        empty_archive, output_fastq, temp_path,
    ));
    assert!(result.is_err(), "Decompressing a zero-byte file should fail");
}

#[test]
fn test_twobit_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGTNNACGT\n+\nIIIIIIIIIIIIII\n@read2\nTGCATGCAATGCAT\n+\nHHHHHHHHHHHHHH\n@read3\nNNNACGTACGTNNN\n+\nAAAAAAAAAAAABB\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions { twobit: true, ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    assert_eq!(original.trim(), decompressed.trim());
}

#[test]
fn test_header_template_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Illumina-style headers with common prefix and varying coordinates
    let test_data = "@INSTRUMENT:123:FLOWCELL:1:2101:1000:2000 1:N:0:ATCG\nACGTACGT\n+\nIIIIIIII\n\
                     @INSTRUMENT:123:FLOWCELL:1:2101:1001:2001 1:N:0:ATCG\nTGCATGCA\n+\nHHHHHHHH\n\
                     @INSTRUMENT:123:FLOWCELL:1:2101:1002:2002 1:N:0:ATCG\nAAAACCCC\n+\nBBBBBBBB\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions { header_template: true, ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    assert_eq!(original.trim(), decompressed.trim());
}

#[test]
fn test_quality_modeling_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nHHHHHHHH\n@read3\nAAAACCCC\n+\nBBBBAAAA\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions { quality_modeling: true, ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    // Quality modeling is lossy (i8 delta clamping), so check sequences match and qualities are close
    let original_content = fs::read_to_string(&input_fastq).unwrap();
    let original_lines: Vec<&str> = original_content.lines().collect();
    let decompressed_content = fs::read_to_string(&output_fastq).unwrap();
    let decompressed_lines: Vec<&str> = decompressed_content.lines().collect();

    // Check same number of lines
    assert_eq!(original_lines.len(), decompressed_lines.len());

    // Check headers and sequences match exactly
    for (i, (orig, decomp)) in original_lines.iter().zip(decompressed_lines.iter()).enumerate() {
        let line_type = i % 4;
        if line_type == 0 || line_type == 1 || line_type == 2 {
            // Header, sequence, separator should match exactly
            assert_eq!(orig.trim(), decomp.trim(), "Line {} mismatch", i);
        }
        // Quality lines may differ slightly due to modeling
    }
}

#[test]
fn test_openzl_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@read2\nTGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHH\n@read3\nAAAACCCCGGGGTTTT\n+\nFFFFFFFFFFFFFFFF\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions {
            quality_compressor: QualityCompressor::OpenZl,
            sequence_compressor: SequenceCompressor::OpenZl,
            header_compressor: HeaderCompressor::OpenZl,
            ..Default::default()
        },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

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
        advanced: AdvancedOptions { quality_compressor: QualityCompressor::Fqzcomp, ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let original_lines: Vec<&str> = original.lines().collect();
    let decompressed_lines: Vec<&str> = decompressed.lines().collect();

    assert_eq!(original_lines.len(), decompressed_lines.len(),
        "Line count mismatch: expected {}, got {}", original_lines.len(), decompressed_lines.len());

    for (i, (orig, decomp)) in original_lines.iter().zip(decompressed_lines.iter()).enumerate() {
        assert_eq!(orig.trim(), decomp.trim(),
            "Line {} mismatch:\n  original:     {:?}\n  decompressed: {:?}", i, orig, decomp);
    }
}

#[test]
fn test_openzl_mixed_compressors() {
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
        advanced: AdvancedOptions {
            quality_compressor: QualityCompressor::OpenZl,
            sequence_compressor: SequenceCompressor::Bsc,
            header_compressor: HeaderCompressor::OpenZl,
            ..Default::default()
        },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

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
fn test_reorder_local_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\
@read4\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions { reorder: Some(ReorderMode::Local), ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    // With reorder, reads may be in a different order but all reads must be present
    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len());

    let mut orig_records: Vec<String> = orig_lines.chunks(4).map(|c| c.join("\n")).collect();
    let mut decomp_records: Vec<String> = decomp_lines.chunks(4).map(|c| c.join("\n")).collect();
    orig_records.sort();
    decomp_records.sort();
    assert_eq!(orig_records, decomp_records);
}

#[test]
fn test_reorder_global_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\
@read4\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions { reorder: Some(ReorderMode::Global), ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path.clone())).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len());

    let mut orig_records: Vec<String> = orig_lines.chunks(4).map(|c| c.join("\n")).collect();
    let mut decomp_records: Vec<String> = decomp_lines.chunks(4).map(|c| c.join("\n")).collect();
    orig_records.sort();
    decomp_records.sort();
    assert_eq!(orig_records, decomp_records);
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
        advanced: AdvancedOptions { sequence_hints: true, ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify encoding_type byte is 4 (at offset 8, after v2 prefix)
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[8], 4, "encoding_type should be 4 for sequence hints");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

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
fn test_sequence_delta_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Include near-duplicate reads to exercise delta path, plus a unique read
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n\
@read3\nACGTACGTACGTACGTACGTACGTATGTACGT\n+\nBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB\n\
@read4\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read5\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n+\nDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\n\
@read6\nACGT\n+\nIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions { sequence_delta: true, ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify encoding_type byte is 5 (at offset 8, after v2 prefix)
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[8], 5, "encoding_type should be 5 for inline delta");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

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
    fs::write(&input_fastq,
        "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n\
         @read2\nACGTACGTACGTACGT\n+\nFFFFFFFFFFFFFFFF\n\
         @read3\nTTTTAAAACCCCGGGG\n+\nAAAAAAAAAAAAAAAA\n\
         @read4\nCCCCGGGGTTTTAAAA\n+\nBBBBBBBBBBBBBBBB\n"
    ).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions { rc_canon: true, ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify encoding_type byte is 6 (at offset 8, after v2 prefix)
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[8], 6, "encoding_type should be 6 for rc_canon");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

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
fn test_local_reorder_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    let test_data = "\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read3\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read4\nGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read5\nACGTACGTACGTACGTACGTACGTACGTACGA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@read6\nTGCATGCATGCATGCATGCATGCATGCATGCT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions { local_reorder: true, ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    let archive_data = fs::read(&archive_path).unwrap();
    // local-reorder now uses ReorderLocal (encoding_type=9) instead of deprecated Delta (8)
    assert_eq!(archive_data[8], 9, "encoding_type should be 9 for local-reorder (ReorderLocal)");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();
    assert!(output_fastq.exists());

    // Local reorder restores original order — verify exact match
    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (o, d)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(o, d, "mismatch at line {}", i);
    }
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
    assert_eq!(archive_data[8], 9, "encoding_type should be 9 for ultra");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();
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

#[test]
fn test_quality_ctx_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Generate 20 fixed-length reads (quality_ctx requires fixed length)
    let mut test_data = String::new();
    let seq_base = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32bp
    let qual_base = "IIHHJJBBCCDDEEFFIIHHJJBBCCDDEEFF"; // 32 chars
    assert_eq!(seq_base.len(), 32);
    assert_eq!(qual_base.len(), 32);
    for i in 0..20 {
        // Rotate sequence to get variety
        let rot = i % seq_base.len();
        let seq: String = seq_base.bytes().cycle().skip(rot).take(32).map(|b| b as char).collect();
        test_data.push_str(&format!("@read{}\n{}\n+\n{}\n", i, seq, qual_base));
    }
    fs::write(&input_fastq, &test_data).unwrap();

    // Force quality_ctx (bypasses the 100K-read auto-selection threshold)
    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions { quality_compressor: QualityCompressor::QualityCtx, ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify archive header has quality_compressor = 4 (QualityCtx), at offset 11 after v2 prefix
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[11], 4, "quality_compressor code should be 4 (QualityCtx)");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();

    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();

    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (o, d)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(o, d, "mismatch at line {}", i);
    }
}

#[test]
fn test_quality_ctx_variable_length_roundtrip() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let input_fastq = temp_path.join("test_input.fastq");
    // Variable-length reads: 20bp, 32bp, 15bp, 40bp, 32bp
    let test_data = "\
@read0\nACGTACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFFIIHH\n\
@read1\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFFIIHHJJBBCCDDEEFF\n\
@read2\nACGTACGTACGTACG\n+\nIIHHJJBBCCDDEEF\n\
@read3\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFFIIHHJJBBCCDDEEFFIIHHJJBB\n\
@read4\nTTTGGGCCCAAAGGGTTTGGGCCCAAAGGGTT\n+\nBBBBAAAABBBBAAAABBBBAAAABBBBAAAA\n";
    fs::write(&input_fastq, test_data).unwrap();

    let archive_path = temp_path.join("test.qz");
    let compress_args = CompressConfig {
        input: vec![input_fastq.clone()],
        output: archive_path.clone(),
        working_dir: temp_path.clone(),
        threads: 1,
        advanced: AdvancedOptions { quality_compressor: QualityCompressor::QualityCtx, ..Default::default() },
        ..CompressConfig::default()
    };

    qz_lib::compression::compress(&compress_args).unwrap();
    assert!(archive_path.exists());

    // Verify archive uses QualityCtx (not fallen back to BSC), at offset 11 after v2 prefix
    let archive_data = fs::read(&archive_path).unwrap();
    assert_eq!(archive_data[11], 4, "quality_compressor should be 4 (QualityCtx), not BSC fallback");

    let output_fastq = temp_path.join("decompressed.fastq");
    qz_lib::compression::decompress(&decompress_args(archive_path, output_fastq.clone(), temp_path)).unwrap();
    assert!(output_fastq.exists());

    let original = fs::read_to_string(&input_fastq).unwrap();
    let decompressed = fs::read_to_string(&output_fastq).unwrap();
    let orig_lines: Vec<&str> = original.lines().collect();
    let decomp_lines: Vec<&str> = decompressed.lines().collect();
    assert_eq!(orig_lines.len(), decomp_lines.len(), "line count mismatch");
    for (i, (o, d)) in orig_lines.iter().zip(decomp_lines.iter()).enumerate() {
        assert_eq!(o, d, "mismatch at line {}", i);
    }
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
// QUALITY ENCODING COMBO TESTS
// ========================================

#[test]
fn test_quality_delta_roundtrip() {
    roundtrip_exact(MULTI_READ_DATA, |c| {
        c.advanced.quality_delta = true;
        c.advanced.sequence_compressor = SequenceCompressor::Zstd;
        c.advanced.quality_compressor = QualityCompressor::Zstd;
        c.advanced.header_compressor = HeaderCompressor::Zstd;
    });
}

#[test]
fn test_dict_training_roundtrip() {
    // Dictionary training needs enough data to build a meaningful dictionary
    let mut data = String::new();
    for i in 0..100 {
        data.push_str(&format!(
            "@read{}\nACGTACGTACGTACGT\n+\nIIHHJJBBCCDDEEFF\n", i
        ));
    }
    roundtrip_exact(&data, |c| {
        c.advanced.dict_training = true;
        c.advanced.dict_size = 32;
        c.advanced.sequence_compressor = SequenceCompressor::Zstd;
        c.advanced.quality_compressor = QualityCompressor::Zstd;
        c.advanced.header_compressor = HeaderCompressor::Zstd;
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
    assert_eq!(fs::metadata(&archive_path).unwrap().len(), original_size,
        "existing archive must be left untouched");

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
    // Flip a byte deep inside a compressed BSC block (not in the header area).
    // The new v3 per-block CRC32 must catch this with a clear "CRC32 mismatch"
    // error instead of bubbling through libbsc as a cryptic block_info failure.
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
    qz_lib::compression::decompress(&decompress_args(archive_path.clone(), output, temp_path.clone())).unwrap();

    // Flip a single byte two-thirds of the way into the archive — well inside
    // a compressed stream payload, past all length/CRC prefixes.
    let mut bytes = fs::read(&archive_path).unwrap();
    let flip_pos = bytes.len() * 2 / 3;
    bytes[flip_pos] ^= 0xFF;
    let bad_path = temp_path.join("bad.qz");
    fs::write(&bad_path, &bytes).unwrap();

    let bad_output = temp_path.join("out.fastq");
    let result = qz_lib::compression::decompress(&decompress_args(bad_path, bad_output, temp_path));
    let err = result.expect_err("bit-flipped archive must error, not silently succeed");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("CRC32 mismatch"),
        "expected CRC32 mismatch error, got: {msg}"
    );
}

#[test]
fn test_decompress_rejects_v2_archive_with_helpful_message() {
    // Construct a minimal byte buffer claiming to be a v2 archive. The decoder
    // must refuse with a clear "please re-encode" message rather than fail
    // somewhere deep inside the block parser.
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
    let err = result.expect_err("v2 archive must be rejected");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("v2") && msg.contains("re-encode"),
        "expected helpful v2-rejection message, got: {msg}"
    );
}

#[test]
fn test_dict_training_rejects_non_zstd_compressor() {
    // dict_training only works with zstd-dict. With other compressors the dict
    // was silently ignored on encode but the archive header still claimed
    // dict_present=true, producing archives that fail to decompress.
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();
    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, "@r1\nACGT\n+\nIIII\n").unwrap();

    let archive_path = temp_path.join("out.qz");
    let result = qz_lib::compression::compress(&CompressConfig {
        input: vec![input_fastq],
        output: archive_path,
        working_dir: temp_path,
        threads: 1,
        advanced: AdvancedOptions {
            dict_training: true,
            dict_size: 4,
            quality_compressor: QualityCompressor::Bsc,
            ..Default::default()
        },
        ..CompressConfig::default()
    });
    let err = result.expect_err("dict_training+Bsc must be rejected, not produce a corrupt archive");
    let msg = format!("{err:#}");
    assert!(msg.contains("--dict-training requires --quality-compressor zstd"),
        "unexpected error: {msg}");
}

#[test]
fn test_quality_modeling_rejects_fqzcomp() {
    // quality_modeling silently substituted BSC for the model deltas when
    // Fqzcomp was selected, but stored Fqzcomp in the header. Result: archive
    // unreadable. Must reject at compress time instead.
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();
    let input_fastq = temp_path.join("input.fastq");
    fs::write(&input_fastq, "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n").unwrap();

    let archive_path = temp_path.join("out.qz");
    let result = qz_lib::compression::compress(&CompressConfig {
        input: vec![input_fastq],
        output: archive_path,
        working_dir: temp_path,
        threads: 1,
        advanced: AdvancedOptions {
            quality_modeling: true,
            quality_compressor: QualityCompressor::Fqzcomp,
            ..Default::default()
        },
        ..CompressConfig::default()
    });
    let err = result.expect_err("quality_modeling+Fqzcomp must be rejected");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("--quality-modeling is incompatible with --quality-compressor fqzcomp"),
        "unexpected error: {msg}"
    );
}

#[test]
fn test_decompress_too_small_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let small_path = temp_path.join("small.qz");
    fs::write(&small_path, &[0u8; 20]).unwrap();

    let output = temp_path.join("out.fastq");
    let result = qz_lib::compression::decompress(&decompress_args(small_path, output, temp_path));
    assert!(result.is_err(), "Archive smaller than header should fail");
}

// ========================================
// EDGE CASE TESTS
// ========================================

#[test]
fn test_pack_dna_2bit_rejects_iupac() {
    // pack_dna_2bit must refuse IUPAC ambiguity codes instead of silently
    // mapping them to A or N (data corruption). Default raw+BSC path stays
    // tolerant; this guard only fires for the 2-bit codec invoked by
    // debruijn/greedy/template hybrids.
    use qz_lib::compression::dna_utils::pack_dna_2bit;
    let pure = vec![b"ACGT".to_vec()];
    assert!(pack_dna_2bit(&pure).is_ok(), "pure ACGT must pack");
    let with_n = vec![b"ACGTN".to_vec()];
    assert!(pack_dna_2bit(&with_n).is_ok(), "ACGTN must pack");
    for &iupac in b"RYSWKMBDHV".iter() {
        let seq = vec![iupac, b'A', b'C', b'G'];
        let err = pack_dna_2bit(&[seq])
            .expect_err(&format!("IUPAC {} must be rejected", iupac as char));
        let msg = format!("{err}");
        assert!(msg.contains("unsupported base"), "unexpected error: {msg}");
    }
}

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
fn test_multiple_input_files_rejected() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path().to_path_buf();

    let f1 = temp_path.join("r1.fastq");
    let f2 = temp_path.join("r2.fastq");
    fs::write(&f1, "@r\nACGT\n+\nIIII\n").unwrap();
    fs::write(&f2, "@r\nACGT\n+\nIIII\n").unwrap();

    let compress_args = CompressConfig {
        input: vec![f1, f2],
        output: temp_path.join("test.qz"),
        working_dir: temp_path,
        threads: 1,
        ..CompressConfig::default()
    };

    let result = qz_lib::compression::compress(&compress_args);
    assert!(result.is_err(), "Multiple input files should be rejected");
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
    }).unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: archive_path,
        working_dir: temp_path,
        num_threads: 1,
        fast: true,
    };
    let result = qz_lib::compression::verify(&verify_config).unwrap();
    assert_eq!(result.mode, qz_lib::compression::VerifyMode::Fast);
    assert_eq!(result.crc32, 0); // not computed in fast mode
    assert!(result.blocks_verified > 0, "fast verify should report block count");
    assert!(result.total_bytes > 0);
}

#[test]
fn test_verify_fast_refuses_ultra_archive() {
    // Ultra-mode archives wrap streams in an outer per-chunk frame with no
    // per-block CRC; fast verify must refuse cleanly rather than mis-walking
    // the frame and reporting spurious "CRC mismatch" on healthy data.
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
        ultra: Some(1),
        ..CompressConfig::default()
    }).unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: archive_path,
        working_dir: temp_path,
        num_threads: 1,
        fast: true,
    };
    let err = qz_lib::compression::verify(&verify_config)
        .expect_err("fast verify must refuse ultra archives");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("does not support encoding_type") && msg.contains("ultra"),
        "unexpected error: {msg}"
    );
}

#[test]
fn test_verify_fast_skips_zstd_quality_stream() {
    // When qualities are stored as raw zstd (no v3 framing), fast verify must
    // skip that stream — not bail with a confusing "num_blocks=4.2 billion"
    // error from misreading the zstd magic as a num_blocks header.
    //
    // The chunked path always coerces qualities to Bsc/Fqzcomp/QualityCtx; to
    // actually produce a Zstd-quality archive we force the in-memory path via
    // quality_delta=true (the chunked-streaming check disables it).
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
            quality_compressor: QualityCompressor::Zstd,
            quality_delta: true, // forces compress_in_memory path
            ..Default::default()
        },
        ..CompressConfig::default()
    }).unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: archive_path,
        working_dir: temp_path,
        num_threads: 1,
        fast: true,
    };
    let result = qz_lib::compression::verify(&verify_config).unwrap();
    assert_eq!(result.mode, qz_lib::compression::VerifyMode::Fast);
    assert!(result.streams_skipped >= 1, "qualities stream must be reported skipped");
    assert!(result.blocks_verified > 0, "v3 streams (headers/sequences/nmasks) must still verify");
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
        advanced: AdvancedOptions { rc_canon: true, ..Default::default() },
        ..CompressConfig::default()
    }).unwrap();

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
    assert!(
        msg.contains("CRC32 mismatch") || msg.contains("rc_flags"),
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
        advanced: AdvancedOptions { bsc_block_size_mb: 1, ..Default::default() },
        ..CompressConfig::default()
    }).unwrap();

    let verify_config = qz_lib::cli::VerifyConfig {
        input: archive_path,
        working_dir: temp_path,
        num_threads: 1,
        fast: true,
    };
    let result = qz_lib::compression::verify(&verify_config).unwrap();
    assert!(result.blocks_verified > 4,
        "expected multiple blocks per stream; got {}", result.blocks_verified);
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
    }).unwrap();

    let mut bytes = fs::read(&archive_path).unwrap();
    let flip_pos = bytes.len() * 2 / 3;
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
    assert!(msg.contains("CRC32 mismatch"), "expected CRC32 mismatch, got: {msg}");
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
