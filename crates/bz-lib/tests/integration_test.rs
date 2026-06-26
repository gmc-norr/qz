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
        let cigar_op: u32 = 10 << 4; // 10M

        // Sequence: ACGTACGTAC in BAM 4-bit encoding
        // A=1, C=2, G=4, T=8 → packed pairs: 0x12, 0x48, 0x12, 0x48, 0x12
        let seq_bytes: Vec<u8> = vec![0x12, 0x48, 0x12, 0x48, 0x12];

        // Quality: Phred 30
        let qual_bytes: Vec<u8> = vec![30u8; 10];

        let block_size: i32 =
            32 + l_read_name as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap(); // refID = 0
        writer.write_all(&(i as i32 * 100).to_le_bytes()).unwrap(); // pos
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

/// Create a BAM whose SEQ fields contain IUPAC ambiguity codes and '=' (every
/// BAM 4-bit nibble value, not just A/C/G/T/N). These must round-trip byte-exact.
fn create_iupac_bam(path: &std::path::Path) {
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

    // SEQ nibbles covering all 16 BAM codes: =,A,C,M,G,R,S,V,T,W,Y,H,K,D,B,N
    // (values 0..=15). Length 16 → 8 packed bytes.
    let seq_nibbles: [u8; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    let mut seq_bytes = Vec::new();
    for pair in seq_nibbles.chunks(2) {
        seq_bytes.push((pair[0] << 4) | pair[1]);
    }
    let l_seq: i32 = 16;

    for i in 0..3i32 {
        let read_name = format!("read{:06}\0", i);
        let l_read_name = read_name.len() as u8;
        let n_cigar_op: u16 = 1;
        let cigar_op: u32 = (l_seq as u32) << 4; // 16M
        let qual_bytes: Vec<u8> = vec![30u8; l_seq as usize];

        let block_size: i32 =
            32 + l_read_name as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(&(i * 100).to_le_bytes()).unwrap();
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
        let cigar_op: u32 = 20 << 4; // 20M

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

/// Create a coordinate-sorted BAM of paired reads to exercise PNEXT/RNEXT
/// mate-derivation (archive v10). Each pair shares a QNAME with opposite
/// read1/read2 bits; for derivable pairs the stored next_pos/next_ref_id point
/// exactly at the mate (so the derivable bit fires). The final pair stores a
/// deliberately *wrong* next_pos on read2 (mate pointer that doesn't match the
/// mate's POS) to force the non-derivable exception path — the roundtrip must
/// still reproduce that stored value byte-exactly.
fn create_paired_bam(path: &std::path::Path, num_pairs: usize) {
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

    // Emit, in coordinate order, read1@(i*200) then read2@(i*200+100) per pair.
    let mut emit = |pos: i32, flag: u16, name: &str, next_pos: i32| {
        let read_name = format!("{name}\0");
        let l_read_name = read_name.len() as u8;
        let n_cigar_op: u16 = 1;
        let l_seq: i32 = 10;
        let cigar_op: u32 = 10 << 4; // 10M
        let seq_bytes: Vec<u8> = vec![0x12, 0x48, 0x12, 0x48, 0x12]; // ACGTACGTAC
        let qual_bytes: Vec<u8> = vec![30u8; 10];
        let block_size: i32 =
            32 + l_read_name as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq;
        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap(); // refID = 0
        writer.write_all(&pos.to_le_bytes()).unwrap();
        writer.write_all(&[l_read_name]).unwrap();
        writer.write_all(&[60u8]).unwrap(); // mapq
        writer.write_all(&0u16.to_le_bytes()).unwrap(); // bin
        writer.write_all(&n_cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&flag.to_le_bytes()).unwrap();
        writer.write_all(&l_seq.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap(); // next_refID = 0 (same ref)
        writer.write_all(&next_pos.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap(); // tlen
        writer.write_all(read_name.as_bytes()).unwrap();
        writer.write_all(&cigar_op.to_le_bytes()).unwrap();
        writer.write_all(&seq_bytes).unwrap();
        writer.write_all(&qual_bytes).unwrap();
    };

    for i in 0..num_pairs {
        let p1 = i as i32 * 200;
        let p2 = p1 + 100;
        let name = format!("pair{:06}", i);
        // read1 (0x1 paired | 0x40 first): mate at p2.
        emit(p1, 0x41, &name, p2);
        // read2 (0x1 paired | 0x80 second): mate at p1, except the last pair
        // stores a wrong pointer to force the non-derivable exception path.
        let r2_next = if i + 1 == num_pairs { 123456 } else { p1 };
        emit(p2, 0x81, &name, r2_next);
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
            (10 << 4),    // 10M
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
            (5 << 4),     // 5M
            (2 << 4) | 1, // 2I
            (8 << 4),     // 8M
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

        let block_size: i32 = 32 + read_name.len() as i32 + seq_bytes.len() as i32 + l_seq;

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
    let mut reader = noodles::bgzf::io::Reader::new(std::fs::File::open(path).unwrap());
    let mut data = Vec::new();
    reader.read_to_end(&mut data).unwrap();
    data
}

fn roundtrip(create_fn: impl Fn(&std::path::Path)) {
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
        original,
        roundtripped,
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
    roundtrip(create_indel_bam);
}

#[test]
fn test_roundtrip_iupac_ambiguity_bases() {
    // SEQ with IUPAC ambiguity codes (M/R/S/V/W/Y/H/K/D/B) and '=' must
    // round-trip byte-exact, not collapse to N.
    roundtrip(create_iupac_bam);
}

#[test]
fn test_lossy_quality_reduction_roundtrips() {
    // Variant-aware quality reduction is lossy on quality but must still
    // round-trip self-consistently: decompress recomputes the per-chunk content
    // CRC over the reconstructed (reduced) records and errors on any mismatch, so
    // a successful decompress proves the reduced quality round-trips byte-exact
    // and that non-quality fields are untouched. (Size is data-dependent — the
    // win is on varied real-world quality, verified on chr20, not on the constant
    // synthetic quality here, where flattening can add entropy.)
    let temp = TempDir::new().unwrap();
    let p = temp.path();
    let input = p.join("in.bam");
    create_overlapping_bam(&input, 400);

    for scheme in [
        bz_lib::FlattenScheme::TwoBin {
            thresh: 20,
            low: 6,
            high: 40,
        },
        bz_lib::FlattenScheme::Coarse8,
        bz_lib::FlattenScheme::Single(40),
    ] {
        let lossy = p.join("lossy.bz");
        let advanced = bz_lib::AdvancedOptions {
            quality_reduction: Some(bz_lib::QualityReduction::level(2, scheme)),
            ..Default::default()
        };
        bz_lib::compress(&CompressConfig {
            input: input.clone(),
            output: lossy.clone(),
            working_dir: p.to_path_buf(),
            advanced,
        })
        .unwrap();

        // Decompress must succeed: the content CRC validates the reduced
        // round-trip. This is the core correctness gate.
        let out = p.join("out.bam");
        bz_lib::decompress(&DecompressConfig {
            input: lossy.clone(),
            output: out.clone(),
            working_dir: p.to_path_buf(),
        })
        .unwrap_or_else(|e| panic!("lossy decompress failed for {scheme:?}: {e}"));
        assert!(
            out.exists(),
            "lossy decompress produced no output for {scheme:?}"
        );
    }
}

#[test]
fn test_compress_rejects_invalid_options_without_panicking() {
    // The public `bz_lib::compress` entry must validate options at the library
    // boundary: a zero `fqz_block_size` previously reached `.step_by(0)`
    // and panicked instead of returning a clean error.
    let temp = TempDir::new().unwrap();
    let p = temp.path();
    let input = p.join("in.bam");
    create_test_bam(&input, 10);

    let advanced = AdvancedOptions {
        fqz_block_size: 0,
        ..Default::default()
    };
    let err = bz_lib::compress(&CompressConfig {
        input,
        output: p.join("out.bz"),
        working_dir: p.to_path_buf(),
        advanced,
    })
    .expect_err("zero fqz_block_size must be rejected");
    assert!(
        err.to_string().contains("fqz_block_size"),
        "expected option-validation error, got: {err}"
    );
}

#[test]
fn test_roundtrip_paired_pnext_rnext() {
    // Exercises PNEXT/RNEXT mate-derivation (v10): mate-in-chunk pairs take the
    // derivable path, the final pair's wrong mate pointer takes the exception
    // path. Full byte-level BAM comparison proves both reconstruct exactly.
    roundtrip(|p| create_paired_bam(p, 50));
}

// --- Edge case tests ---

/// BAM with reads that have all-0xFF quality (unavailable quality).
fn create_bam_with_unavailable_quality(path: &std::path::Path) {
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

    for i in 0i32..50 {
        let read_name = format!("q{:04}\0", i);
        let l_read_name = read_name.len() as u8;
        let l_seq: i32 = 10;
        let n_cigar_op: u16 = 1;
        let cigar_op: u32 = 10 << 4;
        let seq_bytes: Vec<u8> = vec![0x12, 0x48, 0x12, 0x48, 0x12];

        // All quality bytes are 0xFF (unavailable)
        let qual_bytes: Vec<u8> = vec![0xFF; 10];

        let block_size: i32 =
            32 + l_read_name as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(&(i * 100).to_le_bytes()).unwrap();
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

    for i in 0i32..100 {
        let read_name = format!("aux{:04}\0", i);
        let l_read_name = read_name.len() as u8;
        let l_seq: i32 = 10;
        let n_cigar_op: u16 = 1;
        let cigar_op: u32 = 10 << 4;
        let seq_bytes: Vec<u8> = vec![0x12, 0x48, 0x12, 0x48, 0x12];
        let qual_bytes: Vec<u8> = vec![30u8; 10];

        // Aux tags: NM:i:3 + MD:Z:5A4
        let mut aux_bytes = Vec::new();
        // NM tag (type 'C' = u8)
        aux_bytes.extend_from_slice(b"NM");
        aux_bytes.push(b'C'); // type = uint8
        aux_bytes.push(3u8); // value
        // MD tag (type 'Z' = string)
        aux_bytes.extend_from_slice(b"MD");
        aux_bytes.push(b'Z');
        aux_bytes.extend_from_slice(b"5A4\0");

        let block_size: i32 = 32
            + l_read_name as i32
            + 4 * n_cigar_op as i32
            + seq_bytes.len() as i32
            + l_seq
            + aux_bytes.len() as i32;

        writer.write_all(&block_size.to_le_bytes()).unwrap();
        writer.write_all(&0i32.to_le_bytes()).unwrap();
        writer.write_all(&(i * 100).to_le_bytes()).unwrap();
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
    // Tests the BSC fallback path (fqz is skipped when 0xFF is present)
    roundtrip(create_bam_with_unavailable_quality);
}

#[test]
fn test_roundtrip_with_aux_tags() {
    roundtrip(create_bam_with_aux_tags);
}

#[test]
fn test_roundtrip_bsc_quality_compressor() {
    // Test with quality_compressor=1 (BSC instead of fqz)
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 100);

    let archive_path = temp_path.join("test.bz");
    let output_bam = temp_path.join("output.bam");

    let opts = AdvancedOptions {
        quality_compressor: 1, // BSC for quality
        ..Default::default()
    };

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
    assert!(
        err_msg.contains("invalid magic"),
        "Error should mention invalid magic: {err_msg}"
    );
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
    assert!(
        result.crc32 != 0,
        "CRC32 should be non-zero for non-empty archive"
    );
    assert!(result.total_bytes > 0);
    assert!(result.elapsed_secs >= 0.0);
}

#[test]
fn test_verify_paired_bam_with_derived_mate_fields() {
    // v10 PNEXT/RNEXT derivation omits next_ref_id (stream 6) and next_pos
    // (stream 7) for in-chunk mates, so those streams are NOT num_records × 4.
    // The structural validator must account for the next_derivable bitmap rather
    // than rigidly expecting fixed-width streams (the roundtrip already proves the
    // data is intact; verify must not false-positive on it).
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_paired_bam(&input_bam, 50);

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
    });
    assert!(
        result.is_ok(),
        "verify must accept a paired archive with derived mate fields: {:?}",
        result.err()
    );
    assert_eq!(result.unwrap().num_records, 100);
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

    assert_eq!(
        result1.crc32, result2.crc32,
        "CRC32 should be deterministic"
    );
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
    // Verify an archive compressed with BSC quality (not fqz)
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 100);

    let archive_path = temp_path.join("test.bz");

    let opts = AdvancedOptions {
        quality_compressor: 1, // BSC quality
        ..Default::default()
    };

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

/// Decompression must detect a content-fidelity mismatch: the per-chunk
/// content CRC is computed over the original records at compress time and
/// recomputed over the reconstructed records at decompress time. Corrupting the
/// stored content CRC (without touching the compressed payloads or their own
/// CRC) must cause decompression to fail rather than silently emit a BAM that
/// doesn't match the source.
#[test]
fn test_decompress_detects_content_hash_mismatch() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_bam = temp_path.join("input.bam");
    create_test_bam(&input_bam, 5);

    let archive_path = temp_path.join("test.bz");
    bz_lib::compress(&CompressConfig {
        input: input_bam,
        output: archive_path.clone(),
        working_dir: temp_path.to_path_buf(),
        advanced: AdvancedOptions::default(),
    })
    .unwrap();

    // Archive header is 21 bytes + compressed SAM header: magic(2) + version(1) +
    // reserved(1) + flags(1) + num_records(8) + num_chunks(4) + sam_len(4), so the
    // SAM-length field is at bytes 17..21. The first chunk header follows; its
    // content_crc32 sits at offset num_records(4) + chunk_flags(1) + crc32(4) = +9.
    let mut bytes = std::fs::read(&archive_path).unwrap();
    let sam_len = u32::from_le_bytes(bytes[17..21].try_into().unwrap()) as usize;
    let content_crc_off = 21 + sam_len + 9;
    bytes[content_crc_off] ^= 0xFF;
    std::fs::write(&archive_path, &bytes).unwrap();

    let output_bam = temp_path.join("output.bam");
    let result = bz_lib::decompress(&DecompressConfig {
        input: archive_path,
        output: output_bam,
        working_dir: temp_path.to_path_buf(),
    });
    assert!(
        result.is_err(),
        "decompress must reject a content-hash mismatch"
    );
}

// --- NUMA decode-sharding surface (read_bz_chunk_layout / decode_chunk_range) ---

/// Inflate a whole BGZF/BAM file to its uncompressed bytes (magic + header + all
/// records). Robust to BGZF block reframing: two BAMs with identical records but
/// different block boundaries inflate to identical bytes.
fn inflate_bam(path: &std::path::Path) -> Vec<u8> {
    let mut r = noodles::bgzf::io::Reader::new(std::fs::File::open(path).unwrap());
    let mut buf = Vec::new();
    r.read_to_end(&mut buf).unwrap();
    buf
}

#[test]
fn test_numa_decode_chunk_range_layout_and_assembly() {
    // The fixed 28-byte BGZF EOF marker noodles' `finish()` appends. The NUMA driver
    // strips it from every part except the last so the parts concatenate into one
    // valid BAM; the test mirrors that here.
    const BGZF_EOF_LEN: usize = 28;

    let temp = TempDir::new().unwrap();
    let p = temp.path();
    let input = p.join("in.bam");
    // bz's level presets override any explicit chunk_size (smallest = level 1 @ 500K
    // records/chunk), so force multiple chunks with > 500K records at level 1.
    create_test_bam(&input, 600_000);

    let archive = p.join("a.bz");
    bz_lib::compress(&CompressConfig {
        input: input.clone(),
        output: archive.clone(),
        working_dir: p.to_path_buf(),
        advanced: AdvancedOptions {
            level: 1,
            ..Default::default()
        },
    })
    .unwrap();

    // Layout pre-scan: counts sum to the total, offsets strictly increase.
    let layout = bz_lib::compression::read_bz_chunk_layout(&archive).unwrap();
    assert!(layout.num_chunks >= 2, "need multiple chunks to test ranges");
    assert_eq!(layout.per_chunk_reads.len(), layout.num_chunks as usize);
    assert_eq!(layout.chunk_offsets.len(), layout.num_chunks as usize);
    assert_eq!(layout.per_chunk_reads.iter().sum::<u64>(), layout.num_records);
    assert_eq!(layout.num_records, 600_000);
    assert!(
        layout.chunk_offsets.windows(2).all(|w| w[1] > w[0]),
        "chunk offsets must strictly increase"
    );

    // Reference: the full in-process decode.
    let full = p.join("full.bam");
    bz_lib::decompress(&DecompressConfig {
        input: archive.clone(),
        output: full.clone(),
        working_dir: p.to_path_buf(),
    })
    .unwrap();
    let full_bytes = inflate_bam(&full);

    // Two-way split: part A (header + chunks [0,k)), part B (chunks [k,N), no header).
    // Strip A's trailing BGZF EOF marker, concatenate → one valid BAM that must
    // inflate to the same bytes as the full decode (exactly what the NUMA driver does).
    let k = layout.num_chunks / 2;
    let part_a = p.join("a.part");
    let part_b = p.join("b.part");
    let na = bz_lib::compression::decode_chunk_range(&archive, 0, k, 4, &part_a).unwrap();
    let nb = bz_lib::compression::decode_chunk_range(&archive, k, layout.num_chunks, 4, &part_b).unwrap();
    assert_eq!(na + nb, 600_000, "split ranges must cover every record");

    let mut a_bytes = std::fs::read(&part_a).unwrap();
    assert!(a_bytes.len() >= BGZF_EOF_LEN);
    a_bytes.truncate(a_bytes.len() - BGZF_EOF_LEN); // strip A's EOF marker
    let b_bytes = std::fs::read(&part_b).unwrap();
    let assembled = p.join("assembled.bam");
    {
        let mut out = std::fs::File::create(&assembled).unwrap();
        out.write_all(&a_bytes).unwrap();
        out.write_all(&b_bytes).unwrap();
    }
    assert_eq!(
        inflate_bam(&assembled),
        full_bytes,
        "EOF-stripped concatenation of range parts != full decode"
    );
}

#[test]
fn test_decode_chunk_range_rejects_bad_range() {
    let temp = TempDir::new().unwrap();
    let p = temp.path();
    let input = p.join("in.bam");
    create_test_bam(&input, 50);
    let archive = p.join("a.bz");
    bz_lib::compress(&CompressConfig {
        input,
        output: archive.clone(),
        working_dir: p.to_path_buf(),
        advanced: AdvancedOptions {
            chunk_size: 10,
            ..Default::default()
        },
    })
    .unwrap();
    let layout = bz_lib::compression::read_bz_chunk_layout(&archive).unwrap();
    let out = p.join("o.bam");
    assert!(
        bz_lib::compression::decode_chunk_range(&archive, 2, 2, 2, &out).is_err(),
        "start >= end must error"
    );
    assert!(
        bz_lib::compression::decode_chunk_range(&archive, 0, layout.num_chunks + 1, 2, &out).is_err(),
        "end past num_chunks must error"
    );
}

/// NUMA compress-sharding oracle: splitting a BAM into contiguous chunk ranges,
/// compressing each range independently (header in part 0 only), and concatenating
/// the parts must reproduce the single-process archive BYTE-FOR-BYTE. This proves
/// (a) `compress_one_chunk`'s output is independent of where a worker starts
/// (chunk bytes depend only on the chunk's records), and (b) the part layout +
/// pure-concat assembly is correct. A 3-way split also exercises a middle part
/// (chunks-only, not the last range).
#[test]
fn test_numa_compress_chunk_range_byte_identical() {
    let temp = TempDir::new().unwrap();
    let p = temp.path();
    let input = p.join("in.bam");
    // bz level presets override explicit chunk_size (level 1 = 500K/chunk), so
    // > 500K records guarantees multiple chunks.
    create_test_bam(&input, 600_000);

    let make_config = |out: &std::path::Path| CompressConfig {
        input: input.clone(),
        output: out.to_path_buf(),
        working_dir: p.to_path_buf(),
        advanced: AdvancedOptions {
            level: 1,
            ..Default::default()
        },
    };

    // Reference: single-process compress.
    let full = p.join("full.bz");
    bz_lib::compress(&make_config(&full)).unwrap();
    let full_bytes = std::fs::read(&full).unwrap();

    // Prescan: one virtual-offset per chunk, counts sum to the total.
    let cfg = make_config(&full);
    let layout = bz_lib::compression::read_bam_compress_layout(&input, &cfg).unwrap();
    assert_eq!(layout.num_records, 600_000);
    assert!(layout.num_chunks >= 2, "need multiple chunks to test ranges");
    assert_eq!(
        layout.chunk_start_vpos.len(),
        layout.num_chunks as usize,
        "one start offset per chunk"
    );

    // 2-way split: part 0 = header + chunks [0, k); part 1 = chunks [k, N) (no
    // header). This covers both code paths — the header-carrying first part and a
    // headerless, seeked, budget-bounded later part (identical path for any
    // interior or final range).
    let n = layout.num_chunks;
    let k = n / 2;
    assert!(0 < k && k < n, "split point must be interior");

    let bounds = [(0u32, k), (k, n)];
    let mut total = 0u64;
    let mut part_paths = Vec::new();
    for (i, &(a, b)) in bounds.iter().enumerate() {
        let part = p.join(format!("part{i}.bz"));
        let expected_records: u64 =
            layout.per_chunk_reads[a as usize..b as usize].iter().sum();
        let written = bz_lib::compression::compress_chunk_range(
            &input,
            a,
            b,
            layout.chunk_start_vpos[a as usize],
            layout.num_records,
            layout.num_chunks,
            expected_records,
            &cfg,
            &part,
        )
        .unwrap();
        total += written;
        part_paths.push(part);
    }
    assert_eq!(total, 600_000, "split ranges must cover every record exactly once");

    // Pure concat: part 0 carries the global header (with the WHOLE-archive
    // totals), later parts are headerless chunk streams — exactly what the driver
    // does. Result must equal the single-process archive byte-for-byte.
    let assembled = p.join("assembled.bz");
    {
        let mut out = std::fs::File::create(&assembled).unwrap();
        for part in &part_paths {
            out.write_all(&std::fs::read(part).unwrap()).unwrap();
        }
    }
    let assembled_bytes = std::fs::read(&assembled).unwrap();
    assert_eq!(
        assembled_bytes.len(),
        full_bytes.len(),
        "assembled archive size != single-process archive size"
    );
    assert_eq!(
        assembled_bytes, full_bytes,
        "concatenated range parts != single-process archive (chunk bytes are not \
         position-independent, or the part layout is wrong)"
    );

    // Sanity: the assembled archive decodes without error (byte-identity already
    // implies this, but exercise the round-trip explicitly).
    let back = p.join("back.bam");
    bz_lib::decompress(&DecompressConfig {
        input: assembled.clone(),
        output: back.clone(),
        working_dir: p.to_path_buf(),
    })
    .unwrap();
    assert!(
        std::fs::metadata(&back).unwrap().len() > 0,
        "decoded BAM is empty"
    );
}

/// `compress_chunk_range` rejects an empty/inverted chunk range, and bails (no
/// output) when the records it compresses disagree with the prescan's
/// `expected_records` — the defensive guard against a chunk-size divergence
/// silently dropping or duplicating records.
#[test]
fn test_compress_chunk_range_rejects_empty_range() {
    let temp = TempDir::new().unwrap();
    let p = temp.path();
    let input = p.join("in.bam");
    create_test_bam(&input, 1_000);
    let cfg = CompressConfig {
        input: input.clone(),
        output: p.join("unused.bz"),
        working_dir: p.to_path_buf(),
        advanced: AdvancedOptions {
            level: 1,
            ..Default::default()
        },
    };
    let part = p.join("part.bz");
    // Empty/inverted ranges error before any work.
    assert!(
        bz_lib::compression::compress_chunk_range(&input, 2, 2, 0, 1_000, 4, 0, &cfg, &part)
            .is_err(),
        "empty range [2,2) must error"
    );
    assert!(
        bz_lib::compression::compress_chunk_range(&input, 3, 1, 0, 1_000, 4, 0, &cfg, &part)
            .is_err(),
        "inverted range [3,1) must error"
    );
    // A valid single-chunk range (1000 records at level 1) compresses fine when
    // expected_records matches, and BAILS when it does not (the drift guard).
    let layout = bz_lib::compression::read_bam_compress_layout(&input, &cfg).unwrap();
    assert_eq!(layout.num_chunks, 1);
    assert_eq!(layout.num_records, 1_000);
    let vpos0 = layout.chunk_start_vpos[0];
    assert!(
        bz_lib::compression::compress_chunk_range(
            &input, 0, 1, vpos0, 1_000, 1, 1_000, &cfg, &part
        )
        .is_ok(),
        "correct expected_records must succeed"
    );
    let bad_part = p.join("bad.bz");
    assert!(
        bz_lib::compression::compress_chunk_range(
            &input, 0, 1, vpos0, 1_000, 1, 999, &cfg, &bad_part
        )
        .is_err(),
        "wrong expected_records (999 vs 1000) must bail"
    );
    assert!(
        !bad_part.exists(),
        "a bailed compress_chunk_range must leave NO output (atomic temp discarded)"
    );
}
