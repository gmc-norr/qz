use bz_lib::{AdvancedOptions, CompressConfig, DecompressConfig, VerifyConfig};
use std::io::{Read, Write};
use tempfile::TempDir;

/// Create a minimal valid BAM file with synthetic records.
fn create_test_bam(path: &std::path::Path, num_records: usize) {
    let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(path).unwrap());

    // BAM magic
    writer.write_all(b"BAM\x01").unwrap();

    // SAM header
    let header = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n";
    writer
        .write_all(&(header.len() as i32).to_le_bytes())
        .unwrap();
    writer.write_all(header).unwrap();

    // Reference dictionary: 1 reference
    writer.write_all(&1i32.to_le_bytes()).unwrap();
    let ref_name = b"chr1\0";
    writer
        .write_all(&(ref_name.len() as i32).to_le_bytes())
        .unwrap();
    writer.write_all(ref_name).unwrap();
    writer.write_all(&248956422i32.to_le_bytes()).unwrap();

    // Write records
    for i in 0..num_records {
        let read_name = format!("read{:06}\0", i);
        let l_read_name = read_name.len() as u8;
        let n_cigar_op: u16 = 1;
        let l_seq: i32 = 10;

        // CIGAR: 10M
        let cigar_op: u32 = (10 << 4) | 0; // 10M

        // Sequence: ACGTACGTAC in BAM 4-bit encoding
        // A=1, C=2, G=4, T=8 → packed pairs: 0x12, 0x48, 0x12, 0x48, 0x12
        let seq_bytes: Vec<u8> = vec![0x12, 0x48, 0x12, 0x48, 0x12];

        // Quality: Phred 30
        let qual_bytes: Vec<u8> = vec![30u8; 10];

        let block_size: i32 =
            32 + l_read_name as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap(); // refID = 0
        writer
            .write_all(&(i as i32 * 100).to_le_bytes())
            .unwrap(); // pos
        writer.write_all(&[l_read_name]).unwrap();
        writer.write_all(&[60u8]).unwrap(); // mapq
        writer.write_all(&0u16.to_le_bytes()).unwrap(); // bin
        writer.write_all(&n_cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap(); // flag
        writer.write_all(&l_seq.to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap(); // next_refID
        writer.write_all(&(-1i32).to_le_bytes()).unwrap(); // next_pos
        writer.write_all(&0i32.to_le_bytes()).unwrap(); // tlen
        writer.write_all(read_name.as_bytes()).unwrap();
        writer.write_all(&cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&seq_bytes).unwrap();
        writer.write_all(&qual_bytes).unwrap();
        // no aux tags
    }

    writer.finish().unwrap();
}

/// Create a BAM file with more realistic overlapping reads (tests consensus-delta).
fn create_overlapping_bam(path: &std::path::Path, num_records: usize) {
    let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(path).unwrap());

    writer.write_all(b"BAM\x01").unwrap();
    let header = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n";
    writer
        .write_all(&(header.len() as i32).to_le_bytes())
        .unwrap();
    writer.write_all(header).unwrap();
    writer.write_all(&1i32.to_le_bytes()).unwrap();
    let ref_name = b"chr1\0";
    writer
        .write_all(&(ref_name.len() as i32).to_le_bytes())
        .unwrap();
    writer.write_all(ref_name).unwrap();
    writer.write_all(&248956422i32.to_le_bytes()).unwrap();

    // A repeating 20-base "reference" pattern: ACGTACGTACGTACGTACGT
    // BAM nibbles: A=1, C=2, G=4, T=8
    let ref_pattern: Vec<u8> = vec![1, 2, 4, 8, 1, 2, 4, 8, 1, 2, 4, 8, 1, 2, 4, 8, 1, 2, 4, 8];

    for i in 0..num_records {
        let read_name = format!("r{:06}\0", i);
        let l_read_name = read_name.len() as u8;
        let l_seq: i32 = 20;
        let n_cigar_op: u16 = 1;
        let cigar_op: u32 = (20 << 4) | 0; // 20M

        // Position: reads overlap significantly (every 5 bases)
        let pos = (i as i32) * 5;

        // Sequence: mostly matches the reference pattern, with occasional SNP
        let mut nibbles: Vec<u8> = Vec::with_capacity(20);
        for j in 0..20usize {
            let base = ref_pattern[j % ref_pattern.len()];
            // Introduce a mismatch at position 10 for every 3rd read
            if j == 10 && i % 3 == 0 {
                nibbles.push(if base == 1 { 2 } else { 1 }); // A→C or X→A
            } else {
                nibbles.push(base);
            }
        }

        // Pack nibbles
        let mut seq_bytes = Vec::with_capacity(10);
        for chunk in nibbles.chunks(2) {
            seq_bytes.push((chunk[0] << 4) | chunk[1]);
        }

        let qual_bytes: Vec<u8> = vec![35u8; 20];

        let block_size: i32 =
            32 + l_read_name as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(&pos.to_le_bytes()).unwrap();
        writer.write_all(&[l_read_name]).unwrap();
        writer.write_all(&[60u8]).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&n_cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&l_seq.to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(read_name.as_bytes()).unwrap();
        writer.write_all(&cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&seq_bytes).unwrap();
        writer.write_all(&qual_bytes).unwrap();
    }

    writer.finish().unwrap();
}

/// Create a BAM with records that have insertions and soft clips.
fn create_indel_bam(path: &std::path::Path) {
    let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(path).unwrap());

    writer.write_all(b"BAM\x01").unwrap();
    let header = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n";
    writer
        .write_all(&(header.len() as i32).to_le_bytes())
        .unwrap();
    writer.write_all(header).unwrap();
    writer.write_all(&1i32.to_le_bytes()).unwrap();
    let ref_name = b"chr1\0";
    writer
        .write_all(&(ref_name.len() as i32).to_le_bytes())
        .unwrap();
    writer.write_all(ref_name).unwrap();
    writer.write_all(&248956422i32.to_le_bytes()).unwrap();

    // Record 1: 5S 10M (soft clip + match)
    {
        let read_name = b"indel1\0";
        let l_seq: i32 = 15;
        let n_cigar_op: u16 = 2;
        let cigar: Vec<u32> = vec![
            (5 << 4) | 4, // 5S
            (10 << 4) | 0, // 10M
        ];
        // Sequence: 15 bases, A=1 for all
        let nibbles: Vec<u8> = vec![1; 15];
        let mut seq_bytes = Vec::new();
        for chunk in nibbles.chunks(2) {
            if chunk.len() == 2 {
                seq_bytes.push((chunk[0] << 4) | chunk[1]);
            } else {
                seq_bytes.push(chunk[0] << 4);
            }
        }
        let qual_bytes: Vec<u8> = vec![30u8; 15];

        let block_size: i32 =
            32 + read_name.len() as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap(); // refID
        writer.write_all(&100i32.to_le_bytes()).unwrap(); // pos
        writer.write_all(&[read_name.len() as u8]).unwrap();
        writer.write_all(&[60u8]).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&n_cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&l_seq.to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(read_name).unwrap();
        for op in &cigar {
            writer.write_all(&op.to_le_bytes()).unwrap();
        }
        writer.write_all(&seq_bytes).unwrap();
        writer.write_all(&qual_bytes).unwrap();
    }

    // Record 2: 5M 2I 8M (match + insertion + match)
    {
        let read_name = b"indel2\0";
        let l_seq: i32 = 15;
        let n_cigar_op: u16 = 3;
        let cigar: Vec<u32> = vec![
            (5 << 4) | 0,  // 5M
            (2 << 4) | 1,  // 2I
            (8 << 4) | 0,  // 8M
        ];
        let nibbles: Vec<u8> = vec![2; 15]; // all C
        let mut seq_bytes = Vec::new();
        for chunk in nibbles.chunks(2) {
            if chunk.len() == 2 {
                seq_bytes.push((chunk[0] << 4) | chunk[1]);
            } else {
                seq_bytes.push(chunk[0] << 4);
            }
        }
        let qual_bytes: Vec<u8> = vec![25u8; 15];

        let block_size: i32 =
            32 + read_name.len() as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(&200i32.to_le_bytes()).unwrap();
        writer.write_all(&[read_name.len() as u8]).unwrap();
        writer.write_all(&[60u8]).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&n_cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&l_seq.to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(read_name).unwrap();
        for op in &cigar {
            writer.write_all(&op.to_le_bytes()).unwrap();
        }
        writer.write_all(&seq_bytes).unwrap();
        writer.write_all(&qual_bytes).unwrap();
    }

    // Record 3: unmapped (no CIGAR, flag 0x4)
    {
        let read_name = b"unmapped\0";
        let l_seq: i32 = 10;
        let n_cigar_op: u16 = 0;
        let nibbles: Vec<u8> = vec![8; 10]; // all T
        let mut seq_bytes = Vec::new();
        for chunk in nibbles.chunks(2) {
            seq_bytes.push((chunk[0] << 4) | chunk[1]);
        }
        let qual_bytes: Vec<u8> = vec![20u8; 10];

        let block_size: i32 =
            32 + read_name.len() as i32 + seq_bytes.len() as i32 + l_seq;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap(); // unmapped refID
        writer.write_all(&(-1i32).to_le_bytes()).unwrap(); // unmapped pos
        writer.write_all(&[read_name.len() as u8]).unwrap();
        writer.write_all(&[0u8]).unwrap(); // mapq 0
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&n_cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&4u16.to_le_bytes()).unwrap(); // flag 0x4 = unmapped
        writer.write_all(&l_seq.to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(read_name).unwrap();
        // no cigar ops
        writer.write_all(&seq_bytes).unwrap();
        writer.write_all(&qual_bytes).unwrap();
    }

    writer.finish().unwrap();
}

/// BGZF-decompress an entire BAM file to raw bytes for comparison.
fn bgzf_decompress_all(path: &std::path::Path) -> Vec<u8> {
    let mut reader =
        noodles::bgzf::io::Reader::new(std::fs::File::open(path).unwrap());
    let mut data = Vec::new();
    reader.read_to_end(&mut data).unwrap();
    data
}

fn roundtrip(
    create_fn: impl Fn(&std::path::Path),
) {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_bam = temp_path.join("input.bam");
    create_fn(&input_bam);

    let archive_path = temp_path.join("test.bz");
    let output_bam = temp_path.join("output.bam");

    // Compress
    bz_lib::compress(&CompressConfig {
        input: input_bam.clone(),
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: bz_lib::AdvancedOptions::default(),
    })
    .unwrap();

    assert!(archive_path.exists(), "Archive file was not created");

    // Decompress
    bz_lib::decompress(&DecompressConfig {
        input: archive_path,
        output: output_bam.clone(),
        working_dir: temp_path.to_path_buf(),
    })
    .unwrap();

    assert!(output_bam.exists(), "Output BAM was not created");

    // Compare BGZF-decompressed content
    let original = bgzf_decompress_all(&input_bam);
    let roundtripped = bgzf_decompress_all(&output_bam);
    assert_eq!(
        original, roundtripped,
        "BAM roundtrip produced different content (original {} bytes, roundtripped {} bytes)",
        original.len(),
        roundtripped.len()
    );
}

#[test]
fn test_roundtrip_small() {
    roundtrip(|p| create_test_bam(p, 10));
}

#[test]
fn test_roundtrip_medium() {
    roundtrip(|p| create_test_bam(p, 1000));
}

#[test]
fn test_roundtrip_empty() {
    roundtrip(|p| create_test_bam(p, 0));
}

#[test]
fn test_roundtrip_single_record() {
    roundtrip(|p| create_test_bam(p, 1));
}

#[test]
fn test_roundtrip_overlapping() {
    roundtrip(|p| create_overlapping_bam(p, 100));
}

#[test]
fn test_roundtrip_indels_and_unmapped() {
    roundtrip(|p| create_indel_bam(p));
}

// --- Edge case tests ---

/// BAM with reads that have all-0xFF quality (unavailable quality).
fn create_bam_with_unavailable_quality(path: &std::path::Path) {
    let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(path).unwrap());

    writer.write_all(b"BAM\x01").unwrap();
    let header = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n";
    writer.write_all(&(header.len() as i32).to_le_bytes()).unwrap();
    writer.write_all(header).unwrap();
    writer.write_all(&1i32.to_le_bytes()).unwrap();
    let ref_name = b"chr1\0";
    writer.write_all(&(ref_name.len() as i32).to_le_bytes()).unwrap();
    writer.write_all(ref_name).unwrap();
    writer.write_all(&248956422i32.to_le_bytes()).unwrap();

    for i in 0..50 {
        let read_name = format!("q{:04}\0", i);
        let l_read_name = read_name.len() as u8;
        let l_seq: i32 = 10;
        let n_cigar_op: u16 = 1;
        let cigar_op: u32 = (10 << 4) | 0;
        let seq_bytes: Vec<u8> = vec![0x12, 0x48, 0x12, 0x48, 0x12];

        // All quality bytes are 0xFF (unavailable)
        let qual_bytes: Vec<u8> = vec![0xFF; 10];

        let block_size: i32 =
            32 + l_read_name as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(&(i as i32 * 100).to_le_bytes()).unwrap();
        writer.write_all(&[l_read_name]).unwrap();
        writer.write_all(&[60u8]).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&n_cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&l_seq.to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(read_name.as_bytes()).unwrap();
        writer.write_all(&cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&seq_bytes).unwrap();
        writer.write_all(&qual_bytes).unwrap();
    }

    writer.finish().unwrap();
}

/// BAM with reads that have aux tags.
fn create_bam_with_aux_tags(path: &std::path::Path) {
    let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(path).unwrap());

    writer.write_all(b"BAM\x01").unwrap();
    let header = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n";
    writer.write_all(&(header.len() as i32).to_le_bytes()).unwrap();
    writer.write_all(header).unwrap();
    writer.write_all(&1i32.to_le_bytes()).unwrap();
    let ref_name = b"chr1\0";
    writer.write_all(&(ref_name.len() as i32).to_le_bytes()).unwrap();
    writer.write_all(ref_name).unwrap();
    writer.write_all(&248956422i32.to_le_bytes()).unwrap();

    for i in 0..100 {
        let read_name = format!("aux{:04}\0", i);
        let l_read_name = read_name.len() as u8;
        let l_seq: i32 = 10;
        let n_cigar_op: u16 = 1;
        let cigar_op: u32 = (10 << 4) | 0;
        let seq_bytes: Vec<u8> = vec![0x12, 0x48, 0x12, 0x48, 0x12];
        let qual_bytes: Vec<u8> = vec![30u8; 10];

        // Aux tags: NM:i:3 + MD:Z:5A4
        let mut aux_bytes = Vec::new();
        // NM tag (type 'C' = u8)
        aux_bytes.extend_from_slice(b"NM");
        aux_bytes.push(b'C'); // type = uint8
        aux_bytes.push(3u8);  // value
        // MD tag (type 'Z' = string)
        aux_bytes.extend_from_slice(b"MD");
        aux_bytes.push(b'Z');
        aux_bytes.extend_from_slice(b"5A4\0");

        let block_size: i32 =
            32 + l_read_name as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq
            + aux_bytes.len() as i32;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(&(i as i32 * 100).to_le_bytes()).unwrap();
        writer.write_all(&[l_read_name]).unwrap();
        writer.write_all(&[60u8]).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&n_cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&0u16.to_le_bytes()).unwrap();
        writer.write_all(&l_seq.to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&(-1i32).to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(read_name.as_bytes()).unwrap();
        writer.write_all(&cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&seq_bytes).unwrap();
        writer.write_all(&qual_bytes).unwrap();
        writer.write_all(&aux_bytes).unwrap();
    }

    writer.finish().unwrap();
}

#[test]
fn test_roundtrip_unavailable_quality() {
    // Tests the BSC fallback path (quality_ctx is skipped when 0xFF is present)
    roundtrip(|p| create_bam_with_unavailable_quality(p));
}

#[test]
fn test_roundtrip_with_aux_tags() {
    roundtrip(|p| create_bam_with_aux_tags(p));
}

#[test]
fn test_roundtrip_bsc_quality_compressor() {
    // Test with quality_compressor=1 (BSC instead of quality_ctx)
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 100);

    let archive_path = temp_path.join("test.bz");
    let output_bam = temp_path.join("output.bam");

    let mut opts = AdvancedOptions::default();
    opts.quality_compressor = 1; // BSC for quality

    bz_lib::compress(&CompressConfig {
        input: input_bam.clone(),
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: opts,
    })
    .unwrap();

    bz_lib::decompress(&DecompressConfig {
        input: archive_path,
        output: output_bam.clone(),
        working_dir: temp_path.to_path_buf(),
    })
    .unwrap();

    let original = bgzf_decompress_all(&input_bam);
    let roundtripped = bgzf_decompress_all(&output_bam);
    assert_eq!(original, roundtripped);
}

#[test]
fn test_roundtrip_zstd_compressors() {
    // Test with zstd for alignment and aux streams
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 100);

    let archive_path = temp_path.join("test.bz");
    let output_bam = temp_path.join("output.bam");

    let mut opts = AdvancedOptions::default();
    opts.alignment_compressor = 1; // zstd
    opts.aux_compressor = 1;       // zstd

    bz_lib::compress(&CompressConfig {
        input: input_bam.clone(),
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: opts,
    })
    .unwrap();

    bz_lib::decompress(&DecompressConfig {
        input: archive_path,
        output: output_bam.clone(),
        working_dir: temp_path.to_path_buf(),
    })
    .unwrap();

    let original = bgzf_decompress_all(&input_bam);
    let roundtripped = bgzf_decompress_all(&output_bam);
    assert_eq!(original, roundtripped);
}

#[test]
fn test_corrupt_archive_magic() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 10);

    let archive_path = temp_path.join("test.bz");
    let output_bam = temp_path.join("output.bam");

    bz_lib::compress(&CompressConfig {
        input: input_bam,
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: AdvancedOptions::default(),
    })
    .unwrap();

    // Corrupt the magic bytes
    let mut data = std::fs::read(&archive_path).unwrap();
    data[0] = b'X';
    std::fs::write(&archive_path, &data).unwrap();

    let result = bz_lib::decompress(&DecompressConfig {
        input: archive_path,
        output: output_bam,
        working_dir: temp_path.to_path_buf(),
    });

    assert!(result.is_err(), "Should fail on corrupt magic");
    let err_msg = format!("{}", result.unwrap_err());
    assert!(err_msg.contains("invalid magic"), "Error should mention invalid magic: {err_msg}");
}

#[test]
fn test_truncated_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 100);

    let archive_path = temp_path.join("test.bz");
    let output_bam = temp_path.join("output.bam");

    bz_lib::compress(&CompressConfig {
        input: input_bam,
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: AdvancedOptions::default(),
    })
    .unwrap();

    // Truncate to half the size
    let data = std::fs::read(&archive_path).unwrap();
    std::fs::write(&archive_path, &data[..data.len() / 2]).unwrap();

    let result = bz_lib::decompress(&DecompressConfig {
        input: archive_path,
        output: output_bam,
        working_dir: temp_path.to_path_buf(),
    });

    assert!(result.is_err(), "Should fail on truncated archive");
}

#[test]
fn test_config_validation() {
    let mut opts = AdvancedOptions::default();

    // Valid config should pass
    assert!(opts.validate().is_ok());

    // Invalid chunk_size
    opts.chunk_size = 0;
    assert!(opts.validate().is_err());
    opts.chunk_size = 2_500_000;

    // Invalid quality_compressor
    opts.quality_compressor = 5;
    assert!(opts.validate().is_err());
    opts.quality_compressor = 0;

    // Invalid alignment_compressor
    opts.alignment_compressor = 2;
    assert!(opts.validate().is_err());
    opts.alignment_compressor = 0;

    // Invalid aux_compressor
    opts.aux_compressor = 255;
    assert!(opts.validate().is_err());
}

// --- Verify tests ---

#[test]
fn test_verify_valid_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 100);

    let archive_path = temp_path.join("test.bz");

    bz_lib::compress(&CompressConfig {
        input: input_bam,
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: AdvancedOptions::default(),
    })
    .unwrap();

    let result = bz_lib::verify(&VerifyConfig {
        input: archive_path,
    })
    .unwrap();

    assert_eq!(result.num_records, 100);
    assert!(result.crc32 != 0, "CRC32 should be non-zero for non-empty archive");
    assert!(result.total_bytes > 0);
    assert!(result.elapsed_secs >= 0.0);
}

#[test]
fn test_verify_consistent_crc() {
    // Verify twice → same CRC32
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_overlapping_bam(&input_bam, 200);

    let archive_path = temp_path.join("test.bz");

    bz_lib::compress(&CompressConfig {
        input: input_bam,
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: AdvancedOptions::default(),
    })
    .unwrap();

    let result1 = bz_lib::verify(&VerifyConfig {
        input: archive_path.clone(),
    })
    .unwrap();

    let result2 = bz_lib::verify(&VerifyConfig {
        input: archive_path,
    })
    .unwrap();

    assert_eq!(result1.crc32, result2.crc32, "CRC32 should be deterministic");
    assert_eq!(result1.num_records, result2.num_records);
    assert_eq!(result1.total_bytes, result2.total_bytes);
}

#[test]
fn test_verify_corrupted_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 100);

    let archive_path = temp_path.join("test.bz");

    bz_lib::compress(&CompressConfig {
        input: input_bam,
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: AdvancedOptions::default(),
    })
    .unwrap();

    // Corrupt bytes in the middle of the compressed data
    let mut data = std::fs::read(&archive_path).unwrap();
    let mid = data.len() / 2;
    for i in mid..std::cmp::min(mid + 16, data.len()) {
        data[i] ^= 0xFF;
    }
    std::fs::write(&archive_path, &data).unwrap();

    let result = bz_lib::verify(&VerifyConfig {
        input: archive_path,
    });

    assert!(result.is_err(), "Verify should fail on corrupted archive");
}

#[test]
fn test_verify_truncated_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 100);

    let archive_path = temp_path.join("test.bz");

    bz_lib::compress(&CompressConfig {
        input: input_bam,
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: AdvancedOptions::default(),
    })
    .unwrap();

    // Truncate to half
    let data = std::fs::read(&archive_path).unwrap();
    std::fs::write(&archive_path, &data[..data.len() / 2]).unwrap();

    let result = bz_lib::verify(&VerifyConfig {
        input: archive_path,
    });

    assert!(result.is_err(), "Verify should fail on truncated archive");
}

#[test]
fn test_verify_empty_archive() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 0);

    let archive_path = temp_path.join("test.bz");

    bz_lib::compress(&CompressConfig {
        input: input_bam,
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: AdvancedOptions::default(),
    })
    .unwrap();

    let result = bz_lib::verify(&VerifyConfig {
        input: archive_path,
    })
    .unwrap();

    assert_eq!(result.num_records, 0);
}

#[test]
fn test_verify_bsc_quality_compressor() {
    // Verify an archive compressed with BSC quality (not quality_ctx)
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 100);

    let archive_path = temp_path.join("test.bz");

    let mut opts = AdvancedOptions::default();
    opts.quality_compressor = 1; // BSC quality

    bz_lib::compress(&CompressConfig {
        input: input_bam,
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: opts,
    })
    .unwrap();

    let result = bz_lib::verify(&VerifyConfig {
        input: archive_path,
    })
    .unwrap();

    assert_eq!(result.num_records, 100);
    assert!(result.crc32 != 0);
}
