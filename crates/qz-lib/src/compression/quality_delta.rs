//! Quality delta encoding between adjacent reads
//!
//! After reordering reads by sequence similarity, adjacent reads often have
//! similar quality patterns. Instead of storing absolute quality values,
//! we store differences (deltas) from the previous read.
//!
//! Expected improvement: ~20-30% when combined with read reordering

use anyhow::Result;

/// Encode quality strings as deltas from previous read
/// First read is stored as-is, subsequent reads as signed differences
pub fn encode_quality_deltas(quality_strings: &[Vec<u8>]) -> Vec<Vec<i8>> {
    if quality_strings.is_empty() {
        return Vec::new();
    }

    let mut encoded = Vec::with_capacity(quality_strings.len());

    // First quality string is stored as absolute values
    let first_qual = &quality_strings[0];
    let first_deltas: Vec<i8> = first_qual
        .iter()
        .map(|&q| (q.saturating_sub(33)) as i8)
        .collect();
    encoded.push(first_deltas);

    // Subsequent quality strings are stored as deltas
    for i in 1..quality_strings.len() {
        let prev_qual = &quality_strings[i - 1];
        let curr_qual = &quality_strings[i];

        let mut deltas = Vec::with_capacity(curr_qual.len());

        // For each position, compute delta
        for (pos, &curr_byte) in curr_qual.iter().enumerate() {
            let curr_phred = curr_byte.saturating_sub(33) as i16;

            // Get previous quality at this position (or default if lengths differ)
            let prev_phred = if pos < prev_qual.len() {
                prev_qual[pos].saturating_sub(33) as i16
            } else {
                30 // Default quality if previous read is shorter
            };

            // Compute and clamp delta
            let delta = (curr_phred - prev_phred).clamp(-128, 127) as i8;
            deltas.push(delta);
        }

        encoded.push(deltas);
    }

    encoded
}

/// Decode quality strings from delta encoding
pub fn decode_quality_deltas(encoded: &[Vec<i8>]) -> Result<Vec<Vec<u8>>> {
    if encoded.is_empty() {
        return Ok(Vec::new());
    }

    let mut decoded: Vec<Vec<u8>> = Vec::with_capacity(encoded.len());

    // First quality string is stored as absolute values
    let first_qual: Vec<u8> = encoded[0]
        .iter()
        .map(|&phred| {
            let phred_u8 = phred.max(0).min(93) as u8;
            phred_u8 + 33
        })
        .collect();
    decoded.push(first_qual);

    // Decode subsequent quality strings using deltas
    for i in 1..encoded.len() {
        let prev_qual = &decoded[i - 1];
        let deltas = &encoded[i];

        let mut curr_qual = Vec::with_capacity(deltas.len());

        for (pos, &delta) in deltas.iter().enumerate() {
            // Get previous quality at this position
            let prev_phred = if pos < prev_qual.len() {
                prev_qual[pos].saturating_sub(33) as i16
            } else {
                30
            };

            // Reconstruct quality
            let curr_phred = (prev_phred + delta as i16).clamp(0, 93) as u8;
            curr_qual.push(curr_phred + 33);
        }

        decoded.push(curr_qual);
    }

    Ok(decoded)
}

/// Pack delta values efficiently
pub fn pack_deltas(deltas: &[i8]) -> Vec<u8> {
    deltas.iter().map(|&d| d as u8).collect()
}

/// Unpack delta values
pub fn unpack_deltas(packed: &[u8]) -> Vec<i8> {
    packed.iter().map(|&b| b as i8).collect()
}

/// Compute statistics on delta magnitudes
pub fn analyze_deltas(encoded: &[Vec<i8>]) -> DeltaStats {
    let mut total_deltas = 0usize;
    let mut zero_count = 0usize;
    let mut small_count = 0usize;
    let mut max_delta = 0i8;

    for deltas in encoded.iter().skip(1) {
        for &delta in deltas {
            total_deltas += 1;
            let abs_delta = delta.abs();

            if delta == 0 {
                zero_count += 1;
            } else if abs_delta <= 3 {
                small_count += 1;
            }

            if abs_delta > max_delta.abs() {
                max_delta = delta;
            }
        }
    }

    DeltaStats {
        max_delta,
        zero_percent: if total_deltas > 0 {
            (zero_count as f64 / total_deltas as f64) * 100.0
        } else {
            0.0
        },
        small_percent: if total_deltas > 0 {
            (small_count as f64 / total_deltas as f64) * 100.0
        } else {
            0.0
        },
    }
}

pub struct DeltaStats {
    pub max_delta: i8,
    pub zero_percent: f64,
    pub small_percent: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_qualities() {
        let qualities = vec![
            b"IIIIIIIIII".to_vec(),
            b"IIIIIIIIII".to_vec(),
            b"IIIIIIIIII".to_vec(),
        ];

        let encoded = encode_quality_deltas(&qualities);
        let decoded = decode_quality_deltas(&encoded).unwrap();

        assert_eq!(qualities, decoded);

        // All deltas should be 0 for identical quality strings
        for deltas in encoded.iter().skip(1) {
            assert!(deltas.iter().all(|&d| d == 0));
        }
    }

    #[test]
    fn test_similar_qualities() {
        let qualities = vec![
            b"IIIIIHHHHH".to_vec(),
            b"IIIIIHHHHG".to_vec(), // Last base differs by 1
            b"IIIIIHHHHH".to_vec(), // Back to original
        ];

        let encoded = encode_quality_deltas(&qualities);
        let decoded = decode_quality_deltas(&encoded).unwrap();

        assert_eq!(qualities, decoded);

        // Deltas should be small
        let stats = analyze_deltas(&encoded);
        assert!(stats.max_delta.abs() <= 2);
    }

    #[test]
    fn test_different_lengths() {
        let qualities = vec![
            b"IIIII".to_vec(),
            b"IIIIIHHHH".to_vec(), // Longer
            b"III".to_vec(),       // Shorter
        ];

        let encoded = encode_quality_deltas(&qualities);
        let decoded = decode_quality_deltas(&encoded).unwrap();

        assert_eq!(qualities, decoded);
    }

    #[test]
    fn test_empty() {
        let qualities: Vec<Vec<u8>> = vec![];
        let encoded = encode_quality_deltas(&qualities);
        let decoded = decode_quality_deltas(&encoded).unwrap();
        assert_eq!(qualities, decoded);
    }

    #[test]
    fn test_single_quality() {
        let qualities = vec![b"IIIII".to_vec()];
        let encoded = encode_quality_deltas(&qualities);
        let decoded = decode_quality_deltas(&encoded).unwrap();
        assert_eq!(qualities, decoded);
    }
}
