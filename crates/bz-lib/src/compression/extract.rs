use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use tracing::{info, warn};

use crate::cli::ExtractConfig;
use crate::io::bam::RawBamRecord;

// BAM FLAG bits
const FLAG_PAIRED: u16 = 0x1;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_FIRST_IN_PAIR: u16 = 0x40;
const FLAG_SECOND_IN_PAIR: u16 = 0x80;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// BAM 4-bit encoding → ASCII base lookup.
/// Index 0 = '=', 1 = 'A', 2 = 'C', ..., 15 = 'N'
const NIBBLE_TO_BASE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

/// Complement table for BAM 4-bit nibbles.
/// A(1)↔T(8), C(2)↔G(4), M(3)↔K(11), R(5)↔Y(10), S(6)↔S(6), V(7)↔B(14),
/// W(9)↔W(9), H(12)↔D(13), N(15)↔N(15), =(0)↔=(0)
const COMPLEMENT_NIBBLE: [u8; 16] = [0, 8, 4, 11, 2, 10, 6, 14, 1, 9, 5, 3, 13, 12, 7, 15];

/// Convert a BAM record to FASTQ fields: (header, sequence, quality).
/// Handles reverse-complement for reads mapped to the reverse strand.
fn bam_record_to_fastq(record: &RawBamRecord) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let flag = record.flag();

    // Read name: strip NUL terminator, prepend '@'
    let raw_name = record.read_name();
    let name = if raw_name.last() == Some(&0) {
        &raw_name[..raw_name.len() - 1]
    } else {
        raw_name
    };
    let mut header = Vec::with_capacity(1 + name.len());
    header.push(b'@');
    header.extend_from_slice(name);

    let l_seq = record.l_seq() as usize;
    let nibbles = record.unpack_seq_nibbles();
    let qual_bytes = record.qual_bytes();

    if flag & FLAG_REVERSE != 0 {
        // Reverse complement sequence, reverse quality
        let mut seq = Vec::with_capacity(l_seq);
        let mut qual = Vec::with_capacity(l_seq);
        for i in (0..l_seq).rev() {
            seq.push(NIBBLE_TO_BASE[COMPLEMENT_NIBBLE[nibbles[i] as usize] as usize]);
            let q = qual_bytes[i];
            // 0xFF means quality unavailable in BAM → Phred 0
            qual.push(if q == 0xFF { b'!' } else { q + 33 });
        }
        (header, seq, qual)
    } else {
        // Forward strand
        let mut seq = Vec::with_capacity(l_seq);
        let mut qual = Vec::with_capacity(l_seq);
        for i in 0..l_seq {
            seq.push(NIBBLE_TO_BASE[nibbles[i] as usize]);
            let q = qual_bytes[i];
            qual.push(if q == 0xFF { b'!' } else { q + 33 });
        }
        (header, seq, qual)
    }
}

/// Write a FASTQ record (header, sequence, quality) to a writer.
fn write_fastq_record(w: &mut impl Write, header: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    w.write_all(header)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")?;
    Ok(())
}

/// Extract FASTQ from BAM and compress to QZ format.
///
/// For paired-end data, produces `{prefix}_R1.qz` and `{prefix}_R2.qz` with
/// matching read order. For single-end data, produces `{prefix}.qz`.
///
/// Reads are extracted from coordinate-sorted BAM, preserving genomic locality
/// for better BSC compression. Secondary and supplementary alignments are skipped.
pub fn extract(config: &ExtractConfig) -> Result<()> {
    let threads = if config.threads == 0 {
        qz_lib::cli::num_cpus()
    } else {
        config.threads
    };

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;

    pool.install(|| extract_inner(config, threads))
}

fn extract_inner(config: &ExtractConfig, threads: usize) -> Result<()> {
    let start = std::time::Instant::now();

    // Open BAM
    let worker_count = NonZeroUsize::new(4).unwrap();
    let mut reader = crate::io::bam::RawBamReader::from_path_mt(&config.input, worker_count)
        .with_context(|| format!("Failed to open BAM file: {:?}", config.input))?;

    info!("Reading BAM file: {:?}", config.input);

    // First pass: determine if paired-end and extract reads
    let tmp_r1_path = config.working_dir.join(".qz_extract_R1.fastq.tmp");
    let tmp_r2_path = config.working_dir.join(".qz_extract_R2.fastq.tmp");
    let tmp_se_path = config.working_dir.join(".qz_extract_SE.fastq.tmp");

    // Guard for temp file cleanup
    struct TmpCleanup<'a>(&'a [PathBuf]);
    impl Drop for TmpCleanup<'_> {
        fn drop(&mut self) {
            for p in self.0 {
                let _ = std::fs::remove_file(p);
            }
        }
    }
    let tmp_paths = vec![tmp_r1_path.clone(), tmp_r2_path.clone(), tmp_se_path.clone()];
    let _cleanup = TmpCleanup(&tmp_paths);

    let mut r1_writer = BufWriter::with_capacity(
        4 * 1024 * 1024,
        std::fs::File::create(&tmp_r1_path)?,
    );
    let mut r2_writer = BufWriter::with_capacity(
        4 * 1024 * 1024,
        std::fs::File::create(&tmp_r2_path)?,
    );

    // HashMap for pairing: read_name → (Option<(header, seq, qual)>, Option<(header, seq, qual)>)
    // Using (R1, R2) slots
    type FastqFields = (Vec<u8>, Vec<u8>, Vec<u8>);
    let mut pair_buffer: HashMap<Vec<u8>, (Option<FastqFields>, Option<FastqFields>)> =
        HashMap::new();

    let mut total_records: u64 = 0;
    let mut skipped_secondary: u64 = 0;
    let mut paired_count: u64 = 0;
    let mut single_count: u64 = 0;
    let mut is_paired_end = false;
    let mut pairs_written: u64 = 0;

    loop {
        let record = match reader.next_record()? {
            Some(r) => r,
            None => break,
        };
        total_records += 1;

        let flag = record.flag();

        // Skip secondary and supplementary alignments
        if flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
            skipped_secondary += 1;
            continue;
        }

        let (header, seq, qual) = bam_record_to_fastq(&record);

        if flag & FLAG_PAIRED != 0 {
            is_paired_end = true;
            paired_count += 1;

            // Extract read name for pairing (strip '@' prefix)
            let name = header[1..].to_vec();

            let is_r1 = flag & FLAG_FIRST_IN_PAIR != 0;
            let is_r2 = flag & FLAG_SECOND_IN_PAIR != 0;

            if !is_r1 && !is_r2 {
                // Paired but no R1/R2 flag — treat as R1
                warn!("Paired read without R1/R2 flag: {:?}", String::from_utf8_lossy(&name));
            }

            let entry = pair_buffer.entry(name).or_insert((None, None));
            if is_r2 {
                entry.1 = Some((header, seq, qual));
            } else {
                entry.0 = Some((header, seq, qual));
            }

            // Check if pair is complete — write immediately (coordinate order)
            if entry.0.is_some() && entry.1.is_some() {
                let name_key = entry.0.as_ref().unwrap().0[1..].to_vec();
                let (r1, r2) = pair_buffer.remove(&name_key).unwrap();
                let (h1, s1, q1) = r1.unwrap();
                let (h2, s2, q2) = r2.unwrap();
                write_fastq_record(&mut r1_writer, &h1, &s1, &q1)?;
                write_fastq_record(&mut r2_writer, &h2, &s2, &q2)?;
                pairs_written += 1;
            }
        } else {
            single_count += 1;
        }

        if total_records % 10_000_000 == 0 {
            info!(
                "Processed {} records ({} pairs written, {} buffered)",
                total_records, pairs_written, pair_buffer.len()
            );
        }

        if pair_buffer.len() > 10_000_000 {
            warn!(
                "Pair buffer has {} unpaired reads — BAM may not be properly paired",
                pair_buffer.len()
            );
        }
    }

    // Handle remaining unpaired reads
    let unpaired = pair_buffer.len();
    if unpaired > 0 {
        warn!(
            "{} reads remain unpaired after processing all {} records",
            unpaired, total_records
        );
        for (_name, (r1_opt, r2_opt)) in pair_buffer.drain() {
            match (r1_opt, r2_opt) {
                (Some((h1, s1, q1)), None) => {
                    write_fastq_record(&mut r1_writer, &h1, &s1, &q1)?;
                    let dummy_seq = vec![b'N'; s1.len()];
                    let dummy_qual = vec![b'!'; s1.len()];
                    write_fastq_record(&mut r2_writer, &h1, &dummy_seq, &dummy_qual)?;
                    pairs_written += 1;
                }
                (None, Some((h2, s2, q2))) => {
                    let dummy_seq = vec![b'N'; s2.len()];
                    let dummy_qual = vec![b'!'; s2.len()];
                    write_fastq_record(&mut r1_writer, &h2, &dummy_seq, &dummy_qual)?;
                    write_fastq_record(&mut r2_writer, &h2, &s2, &q2)?;
                    pairs_written += 1;
                }
                _ => {}
            }
        }
    }

    r1_writer.flush()?;
    r2_writer.flush()?;
    drop(r1_writer);
    drop(r2_writer);

    let read_time = start.elapsed();
    info!(
        "BAM reading complete in {:.1}s: {} total records, {} paired, {} single-end, {} pairs written, {} secondary/supplementary skipped",
        read_time.as_secs_f64(), total_records, paired_count, single_count, pairs_written, skipped_secondary
    );

    if !is_paired_end && single_count > 0 {
        // Single-end: rename R1 temp to SE, compress as single file
        std::fs::rename(&tmp_r1_path, &tmp_se_path)?;
        let _ = std::fs::remove_file(&tmp_r2_path);

        let output_path = format!("{}.qz", config.output_prefix);
        info!("Compressing single-end reads to {}", output_path);
        compress_fastq_to_qz(&tmp_se_path, &output_path, &config.working_dir, threads)?;
    } else {
        // Paired-end: compress R1 and R2
        let r1_output = format!("{}_R1.qz", config.output_prefix);
        let r2_output = format!("{}_R2.qz", config.output_prefix);

        info!(
            "Compressing {} paired reads to {} and {}",
            pairs_written, r1_output, r2_output
        );

        // Compress R1 and R2 sequentially (each uses all threads internally)
        info!("Compressing R1...");
        compress_fastq_to_qz(&tmp_r1_path, &r1_output, &config.working_dir, threads)?;
        info!("Compressing R2...");
        compress_fastq_to_qz(&tmp_r2_path, &r2_output, &config.working_dir, threads)?;
    }

    let total_time = start.elapsed();
    info!("Extraction complete in {:.1}s", total_time.as_secs_f64());

    Ok(())
}

/// Compress a FASTQ file to QZ format using default settings.
fn compress_fastq_to_qz(
    input: &std::path::Path,
    output: &str,
    working_dir: &std::path::Path,
    threads: usize,
) -> Result<()> {
    let config = qz_lib::cli::CompressConfig {
        input: vec![input.to_path_buf()],
        output: PathBuf::from(output),
        working_dir: working_dir.to_path_buf(),
        threads,
        ..Default::default()
    };
    qz_lib::compression::compress(&config)?;
    Ok(())
}
