//! N-mask encoding for lossless DNA sequences with sparse N bases
//!
//! Strategy: 2-bit encoding + bitmap of N positions
//! - Sequences stored in 2-bit (N→A temporarily): 4 bases/byte
//! - N-mask bitmap (1 bit per position): 8 positions/byte
//! - Total: ~2.125 bits/base for sequences with sparse N
//!
//! This is optimal for real Illumina data where:
//! - 62% of reads contain N bases
//! - But typically only 1-5 N bases per read (sparse)

/// N-mask encoded sequence: 2-bit data + bitmap of N positions
#[derive(Debug, Clone)]
pub struct NMaskEncoding {
    pub sequence_2bit: Vec<u8>,  // 2-bit encoded (N→A)
    pub n_mask: Vec<u8>,          // Bitmap: 1 bit per base (1 = is N)
    pub length: usize,            // Original sequence length
}

/// Lookup table: base ASCII → 2-bit encoding (N/unknown → 0 = A).
static NMASK_BASE_2BIT: [u8; 256] = {
    let mut t = [0u8; 256]; // default 0 covers A, a, N, n, and unknowns
    t[b'C' as usize] = 0b01; t[b'c' as usize] = 0b01;
    t[b'G' as usize] = 0b10; t[b'g' as usize] = 0b10;
    t[b'T' as usize] = 0b11; t[b't' as usize] = 0b11;
    t
};

/// Lookup table: base ASCII → true if N/n.
static NMASK_IS_N: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'N' as usize] = 1;
    t[b'n' as usize] = 1;
    t
};

/// Encode sequence with N-mask: 2-bit + bitmap
pub fn encode_with_n_mask(seq: &[u8]) -> NMaskEncoding {
    let len = seq.len();

    // Allocate 2-bit sequence storage (4 bases per byte)
    let mut sequence_2bit = vec![0u8; (len + 3) / 4];

    // Allocate N-mask bitmap storage (8 bases per byte)
    let mut n_mask = vec![0u8; (len + 7) / 8];

    // Process 4 bases at a time (one packed byte) to avoid per-base division
    let full_chunks = len / 4;
    for chunk_idx in 0..full_chunks {
        let base_off = chunk_idx * 4;
        let b0 = seq[base_off];
        let b1 = seq[base_off + 1];
        let b2 = seq[base_off + 2];
        let b3 = seq[base_off + 3];

        // Pack 4 bases into one byte using LUT
        sequence_2bit[chunk_idx] = NMASK_BASE_2BIT[b0 as usize]
            | (NMASK_BASE_2BIT[b1 as usize] << 2)
            | (NMASK_BASE_2BIT[b2 as usize] << 4)
            | (NMASK_BASE_2BIT[b3 as usize] << 6);

        // Mark N positions in bitmap using LUT
        if NMASK_IS_N[b0 as usize] != 0 { n_mask[base_off / 8] |= 1 << (base_off % 8); }
        if NMASK_IS_N[b1 as usize] != 0 { n_mask[(base_off + 1) / 8] |= 1 << ((base_off + 1) % 8); }
        if NMASK_IS_N[b2 as usize] != 0 { n_mask[(base_off + 2) / 8] |= 1 << ((base_off + 2) % 8); }
        if NMASK_IS_N[b3 as usize] != 0 { n_mask[(base_off + 3) / 8] |= 1 << ((base_off + 3) % 8); }
    }

    // Handle remaining 0-3 bases
    for i in (full_chunks * 4)..len {
        let base = seq[i];
        let byte_idx = i / 4;
        let bit_offset = (i % 4) * 2;
        sequence_2bit[byte_idx] |= NMASK_BASE_2BIT[base as usize] << bit_offset;

        if NMASK_IS_N[base as usize] != 0 {
            n_mask[i / 8] |= 1 << (i % 8);
        }
    }

    NMaskEncoding {
        sequence_2bit,
        n_mask,
        length: len,
    }
}

/// Decode sequence with N-mask: restore original N positions
pub fn decode_with_n_mask(encoding: &NMaskEncoding) -> Vec<u8> {
    let mut result = Vec::with_capacity(encoding.length);

    for i in 0..encoding.length {
        // Get 2-bit base
        let byte_idx = i / 4;
        let bit_offset = (i % 4) * 2;
        let bits = (encoding.sequence_2bit[byte_idx] >> bit_offset) & 0b11;

        // Check if this position is N (from bitmap)
        let mask_byte_idx = i / 8;
        let mask_bit_offset = i % 8;
        let is_n = (encoding.n_mask[mask_byte_idx] >> mask_bit_offset) & 1 == 1;

        let base = if is_n {
            b'N'
        } else {
            match bits {
                0b00 => b'A',
                0b01 => b'C',
                0b10 => b'G',
                0b11 => b'T',
                _ => unreachable!(),
            }
        };

        result.push(base);
    }

    result
}

/// Analyze N-base statistics for a sequence
#[cfg(test)]
fn analyze_n_bases(seq: &str) -> NBaseStats {
    let mut n_count = 0;
    let mut n_positions = Vec::new();

    for (i, c) in seq.chars().enumerate() {
        if c.to_ascii_uppercase() == 'N' {
            n_count += 1;
            n_positions.push(i);
        }
    }

    NBaseStats {
        n_count,
        n_positions,
        n_percentage: if seq.is_empty() {
            0.0
        } else {
            (n_count as f64 / seq.len() as f64) * 100.0
        },
    }
}

#[cfg(test)]
#[derive(Debug, Clone)]
struct NBaseStats {
    n_count: usize,
    n_positions: Vec<usize>,
    n_percentage: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_no_n() {
        let seq = b"ACGTACGT";
        let encoded = encode_with_n_mask(seq);
        let decoded = decode_with_n_mask(&encoded);
        assert_eq!(decoded, b"ACGTACGT");
    }

    #[test]
    fn test_encode_decode_with_n() {
        let seq = b"ACNGTANNCGT";
        let encoded = encode_with_n_mask(seq);
        let decoded = decode_with_n_mask(&encoded);
        assert_eq!(decoded, b"ACNGTANNCGT");
    }

    #[test]
    fn test_encode_decode_all_n() {
        let seq = b"NNNNNNNN";
        let encoded = encode_with_n_mask(seq);
        let decoded = decode_with_n_mask(&encoded);
        assert_eq!(decoded, b"NNNNNNNN");
    }

    #[test]
    fn test_sparse_n() {
        // Typical Illumina: 1-5 N bases in 100bp read
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let encoded = encode_with_n_mask(seq);
        let decoded = decode_with_n_mask(&encoded);
        assert_eq!(decoded, seq.as_slice());

        // Check storage efficiency (sequence is actually 104 bases)
        // 2-bit: 104 bases / 4 = 26 bytes
        // N-mask: 104 bases / 8 = 13 bytes
        // Total: 39 bytes vs 104 bytes (byte-aligned)
        assert_eq!(encoded.sequence_2bit.len(), 26);
        assert_eq!(encoded.n_mask.len(), 13);
    }

    #[test]
    fn test_analyze_n_bases() {
        let seq = "ACNGTNAN";
        let stats = analyze_n_bases(seq);
        assert_eq!(stats.n_count, 3);
        assert_eq!(stats.n_positions, vec![2, 5, 7]);
        assert!((stats.n_percentage - 37.5).abs() < 0.01);
    }
}
