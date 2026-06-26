/// Columnar compression for FASTQ files (patent-free alternative)
///
/// Strategy:
/// 1. Separate FASTQ into 3 streams (headers, sequences, quality)
/// 2. Encode each optimally (2-bit DNA, quality binning)
/// 3. Compress each with BSC
/// 4. No reordering required!
use anyhow::Result;
use std::io::Write;

/// Quality binning schemes
#[derive(Debug, Clone, Copy)]
pub enum QualityBinning {
    /// Illumina 8-level binning (3 bits per quality)
    Illumina8,
    /// Binary threshold (1 bit per quality)
    Binary { threshold: u8 },
    /// 4-level binning (2 bits per quality)
    FourLevel,
    /// No binning (keep original, 7 bits per quality)
    None,
}

impl QualityBinning {
    /// Get number of bits per quality score
    pub fn bits_per_quality(&self) -> usize {
        match self {
            QualityBinning::Illumina8 => 3,
            QualityBinning::Binary { .. } => 1,
            QualityBinning::FourLevel => 2,
            QualityBinning::None => 7,
        }
    }

    /// Encode a quality score
    pub fn encode(&self, phred: u8) -> u8 {
        match self {
            QualityBinning::Illumina8 => {
                // Illumina 8-level binning
                match phred {
                    0..=1 => 0,
                    2..=9 => 1,
                    10..=19 => 2,
                    20..=24 => 3,
                    25..=29 => 4,
                    30..=34 => 5,
                    35..=39 => 6,
                    _ => 7,
                }
            }
            QualityBinning::Binary { threshold } => {
                if phred >= *threshold {
                    1
                } else {
                    0
                }
            }
            QualityBinning::FourLevel => match phred {
                0..=9 => 0,
                10..=19 => 1,
                20..=29 => 2,
                _ => 3,
            },
            QualityBinning::None => phred.min(127),
        }
    }

    /// Decode a quality score
    pub fn decode(&self, encoded: u8) -> u8 {
        match self {
            QualityBinning::Illumina8 => {
                // Use representative values
                match encoded {
                    0 => 6,
                    1 => 6,
                    2 => 15,
                    3 => 22,
                    4 => 27,
                    5 => 33,
                    6 => 37,
                    7 => 40,
                    _ => 40,
                }
            }
            QualityBinning::Binary { threshold } => {
                if encoded == 1 {
                    40
                } else {
                    *threshold / 2
                }
            }
            QualityBinning::FourLevel => match encoded {
                0 => 6,
                1 => 15,
                2 => 25,
                3 => 37,
                _ => 37,
            },
            QualityBinning::None => encoded,
        }
    }
}

/// Pack quality scores with variable bit width.
///
/// For the lossless `None` binning (7-bit), a quality byte whose Phred value
/// exceeds 127 (ASCII > 160 — outside valid Phred+33 FASTQ) cannot be
/// represented and would otherwise be silently clamped; this is rejected so the
/// "lossless" contract is never silently violated. Lossy binnings clamp by
/// design and never error.
pub fn pack_qualities(qualities: &[u8], binning: QualityBinning) -> Result<Vec<u8>> {
    let bits_per_qual = binning.bits_per_quality();

    let mut packed = Vec::with_capacity((qualities.len() * bits_per_qual).div_ceil(8));
    let mut buffer = 0u64;
    let mut bits_in_buffer = 0;

    for &qual_ascii in qualities {
        // Reject sub-33 bytes BEFORE the saturating_sub below silently collapses
        // them to phred 0 (which would round-trip back to '!' = 33 on decode —
        // silent lossless violation). Mirror the high-end guard so the "lossless"
        // contract is never violated silently for out-of-spec / corrupt quality lines.
        if matches!(binning, QualityBinning::None) && qual_ascii < 33 {
            anyhow::bail!(
                "lossless quality packing cannot represent quality byte {} (< 33); \
                 not valid Phred+33 FASTQ quality",
                qual_ascii
            );
        }
        // Convert ASCII to Phred (assuming Phred+33 encoding)
        let phred = qual_ascii.saturating_sub(33);
        if matches!(binning, QualityBinning::None) && phred > 127 {
            anyhow::bail!(
                "lossless quality packing cannot represent quality byte {} (phred {} > 127); \
                 not valid Phred+33 FASTQ quality",
                qual_ascii,
                phred
            );
        }
        let encoded = binning.encode(phred);

        // Add to buffer
        buffer |= (encoded as u64) << bits_in_buffer;
        bits_in_buffer += bits_per_qual;

        // Flush complete bytes
        while bits_in_buffer >= 8 {
            packed.push((buffer & 0xFF) as u8);
            buffer >>= 8;
            bits_in_buffer -= 8;
        }
    }

    // Flush remaining bits
    if bits_in_buffer > 0 {
        packed.push((buffer & 0xFF) as u8);
    }

    Ok(packed)
}

/// Unpack quality scores
pub fn unpack_qualities(packed: &[u8], length: usize, binning: QualityBinning) -> Vec<u8> {
    let bits_per_qual = binning.bits_per_quality();
    let mask = (1u64 << bits_per_qual) - 1;

    let mut qualities = Vec::with_capacity(length);
    let mut buffer = 0u64;
    let mut bits_in_buffer = 0;
    let mut packed_idx = 0;

    for _ in 0..length {
        // Ensure we have enough bits
        while bits_in_buffer < bits_per_qual && packed_idx < packed.len() {
            buffer |= (packed[packed_idx] as u64) << bits_in_buffer;
            bits_in_buffer += 8;
            packed_idx += 1;
        }

        // Extract quality
        let encoded = (buffer & mask) as u8;
        let phred = binning.decode(encoded);
        qualities.push(phred + 33); // Convert to ASCII

        buffer >>= bits_per_qual;
        bits_in_buffer = bits_in_buffer.saturating_sub(bits_per_qual);
    }

    qualities
}

/// Unpack quality scores directly to a writer, avoiding intermediate String allocation
pub fn unpack_qualities_to_writer<W: Write + ?Sized>(
    packed: &[u8],
    length: usize,
    binning: QualityBinning,
    writer: &mut W,
) -> std::io::Result<()> {
    let bits_per_qual = binning.bits_per_quality();
    let mask = (1u64 << bits_per_qual) - 1;

    let mut buffer = 0u64;
    let mut bits_in_buffer = 0;
    let mut packed_idx = 0;

    // Process in chunks to reduce write syscalls
    let mut out_buf = [0u8; 256];
    let mut out_pos = 0;

    for _ in 0..length {
        while bits_in_buffer < bits_per_qual && packed_idx < packed.len() {
            buffer |= (packed[packed_idx] as u64) << bits_in_buffer;
            bits_in_buffer += 8;
            packed_idx += 1;
        }

        let encoded = (buffer & mask) as u8;
        let phred = binning.decode(encoded);
        out_buf[out_pos] = phred + 33;
        out_pos += 1;

        if out_pos == 256 {
            writer.write_all(&out_buf)?;
            out_pos = 0;
        }

        buffer >>= bits_per_qual;
        bits_in_buffer = bits_in_buffer.saturating_sub(bits_per_qual);
    }

    if out_pos > 0 {
        writer.write_all(&out_buf[..out_pos])?;
    }

    Ok(())
}

/// Columnar compression statistics
#[derive(Debug)]
#[allow(dead_code)] // some fields are written for logging but not read back
pub struct ColumnarStats {
    pub num_reads: usize,
    pub total_bases: usize,
    pub original_size: usize,
    pub headers_size: usize,
    pub sequences_size: usize,
    pub qualities_size: usize,
    pub total_compressed: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_qualities_none_is_lossless_or_errors() {
        // Lossless (None) packing of a valid Phred+33 quality string round-trips.
        let quals = b"IIIIB!~"; // ASCII 33..126 range
        let packed = pack_qualities(quals, QualityBinning::None).unwrap();
        let unpacked = unpack_qualities(&packed, quals.len(), QualityBinning::None);
        assert_eq!(unpacked, quals);

        // A quality byte that can't be represented losslessly (phred > 127,
        // i.e. ASCII > 160 — not valid FASTQ) must error, not silently clamp.
        let bad = [200u8]; // phred 167
        assert!(pack_qualities(&bad, QualityBinning::None).is_err());

        // A sub-33 byte (e.g. ASCII space 32, or any control/corrupt byte) used
        // to saturate to phred 0 and round-trip back to '!' (33) — silent data
        // loss. It must now error instead of corrupting.
        for bad_byte in [0u8, 10, 31, 32] {
            assert!(
                pack_qualities(&[bad_byte], QualityBinning::None).is_err(),
                "sub-33 byte {bad_byte} should be rejected, not silently clamped"
            );
        }
    }
}
