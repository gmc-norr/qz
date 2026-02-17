use bz_lib::{CompressConfig, DecompressConfig};
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
        threads: 1,
    })
    .unwrap();

    assert!(archive_path.exists(), "Archive file was not created");

    // Decompress
    bz_lib::decompress(&DecompressConfig {
        input: archive_path,
        output: output_bam.clone(),
        working_dir: temp_path.to_path_buf(),
        threads: 1,
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
