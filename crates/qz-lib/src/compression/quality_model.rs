//! Positional quality score modeling
//!
//! Quality scores degrade predictably along read length in Illumina data.
//! Instead of storing absolute quality values, we:
//! 1. Build a position-specific quality model (median quality at each position)
//! 2. Store only deviations from the model
//! 3. Deviations are smaller and compress better
//!
//! Expected improvement: ~20-30% reduction in quality stream size

use anyhow::Result;
use std::collections::HashMap;

/// Quality model: expected quality score at each position
#[derive(Debug, Clone)]
pub struct QualityModel {
    /// Expected quality at each position (0-indexed)
    pub positional_medians: Vec<u8>,
    /// Read length this model applies to
    pub read_length: usize,
}

/// Build a quality model from a set of reads
pub fn build_quality_model(quality_strings: &[Vec<u8>]) -> Result<QualityModel> {
    if quality_strings.is_empty() {
        anyhow::bail!("Cannot build quality model from empty dataset");
    }

    // Find the most common read length
    let mut length_counts: HashMap<usize, usize> = HashMap::new();
    for qual in quality_strings {
        *length_counts.entry(qual.len()).or_insert(0) += 1;
    }

    let read_length = length_counts
        .into_iter()
        .max_by_key(|(_, count)| *count)
        .map(|(len, _)| len)
        .unwrap_or(150);

    // Collect quality scores at each position
    let mut position_qualities: Vec<Vec<u8>> = vec![Vec::new(); read_length];

    for qual in quality_strings {
        if qual.len() != read_length {
            continue; // Skip reads of different lengths
        }

        for (pos, &qual_ascii) in qual.iter().enumerate() {
            let phred = qual_ascii.saturating_sub(33);
            position_qualities[pos].push(phred);
        }
    }

    // Calculate median quality at each position
    let mut positional_medians = Vec::with_capacity(read_length);
    for mut qualities_at_pos in position_qualities {
        if qualities_at_pos.is_empty() {
            positional_medians.push(30); // Default quality
        } else {
            qualities_at_pos.sort_unstable();
            let median = qualities_at_pos[qualities_at_pos.len() / 2];
            positional_medians.push(median);
        }
    }

    Ok(QualityModel {
        positional_medians,
        read_length,
    })
}

/// Encode quality string using positional model
/// Returns delta-encoded values (signed differences from expected)
pub fn encode_with_model(quality: &[u8], model: &QualityModel) -> Vec<i8> {
    let mut deltas = Vec::with_capacity(quality.len());

    for (pos, &qual_ascii) in quality.iter().enumerate() {
        let phred = qual_ascii.saturating_sub(33) as i16;

        // Get expected quality at this position
        let expected = if pos < model.positional_medians.len() {
            model.positional_medians[pos] as i16
        } else {
            // For reads longer than model, use last position's median
            *model.positional_medians.last().unwrap_or(&30) as i16
        };

        // Store delta (clamped to i8 range)
        let delta = (phred - expected).clamp(-128, 127) as i8;
        deltas.push(delta);
    }

    deltas
}

/// Decode quality string from delta-encoded values
pub fn decode_with_model(deltas: &[i8], model: &QualityModel) -> Vec<u8> {
    let mut quality_bytes = Vec::with_capacity(deltas.len());

    for (pos, &delta) in deltas.iter().enumerate() {
        // Get expected quality at this position
        let expected = if pos < model.positional_medians.len() {
            model.positional_medians[pos] as i16
        } else {
            *model.positional_medians.last().unwrap_or(&30) as i16
        };

        // Reconstruct quality
        let phred = (expected + delta as i16).clamp(0, 93) as u8;
        let qual_ascii = phred + 33;
        quality_bytes.push(qual_ascii);
    }

    quality_bytes
}

/// Serialize quality model to bytes
pub fn serialize_model(model: &QualityModel) -> Vec<u8> {
    let mut bytes = Vec::new();

    // Write read length (2 bytes)
    bytes.extend_from_slice(&(model.read_length as u16).to_le_bytes());

    // Write model length (2 bytes, should match read_length but stored separately)
    bytes.extend_from_slice(&(model.positional_medians.len() as u16).to_le_bytes());

    // Write median quality values (1 byte each)
    bytes.extend_from_slice(&model.positional_medians);

    bytes
}

/// Deserialize quality model from bytes
pub fn deserialize_model(bytes: &[u8]) -> Result<QualityModel> {
    if bytes.len() < 4 {
        anyhow::bail!("Invalid quality model: too short");
    }

    let read_length = u16::from_le_bytes([bytes[0], bytes[1]]) as usize;
    let model_length = u16::from_le_bytes([bytes[2], bytes[3]]) as usize;

    if bytes.len() < 4 + model_length {
        anyhow::bail!("Invalid quality model: truncated");
    }

    let positional_medians = bytes[4..4 + model_length].to_vec();

    Ok(QualityModel {
        positional_medians,
        read_length,
    })
}

/// Pack delta values efficiently
/// Since deltas are typically small, we can use variable-length encoding
pub fn pack_deltas(deltas: &[i8]) -> Vec<u8> {
    // Simple approach: store as signed bytes directly
    // TODO: Could use more sophisticated encoding for better compression
    deltas.iter().map(|&d| d as u8).collect()
}

/// Unpack delta values
pub fn unpack_deltas(packed: &[u8]) -> Vec<i8> {
    packed.iter().map(|&b| b as i8).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_model() {
        let qualities = vec![
            b"IIIIIIIIHHHHGGGGFFFF".to_vec(),
            b"IIIIIIIIHHHHGGGGEEEE".to_vec(),
            b"IIIIIIIIHHHHGGGGFFFF".to_vec(),
        ];

        let model = build_quality_model(&qualities).unwrap();
        assert_eq!(model.read_length, 20);
        assert_eq!(model.positional_medians.len(), 20);

        // Quality degrades: I(40) -> H(39) -> G(38) -> F(37)/E(36)
        assert!(model.positional_medians[0] >= 39); // High quality at start
        assert!(model.positional_medians[19] >= 36); // Lower at end
    }

    #[test]
    fn test_encode_decode() {
        let qualities = vec![
            b"IIIIIHHHHGGGGFFFF".to_vec(),
            b"IIIIIHHHHGGGGEEEE".to_vec(),
        ];

        let model = build_quality_model(&qualities).unwrap();
        let test_qual = b"IIIIIHHHHGGGGFFFF";

        let deltas = encode_with_model(test_qual, &model);
        let decoded = decode_with_model(&deltas, &model);

        assert_eq!(test_qual.to_vec(), decoded);
    }

    #[test]
    fn test_serialize_deserialize() {
        let model = QualityModel {
            positional_medians: vec![40, 40, 39, 38, 37, 36, 35],
            read_length: 7,
        };

        let bytes = serialize_model(&model);
        let restored = deserialize_model(&bytes).unwrap();

        assert_eq!(model.read_length, restored.read_length);
        assert_eq!(model.positional_medians, restored.positional_medians);
    }

    #[test]
    fn test_deltas_are_small() {
        let qualities = vec![
            b"IIIIIHHHHHGGGGGFFFFF".to_vec(),
            b"IIIIIHHHHHGGGGGFFFFF".to_vec(),
            b"IIIIIHHHHHGGGGGFFFFF".to_vec(),
        ];

        let model = build_quality_model(&qualities).unwrap();

        // Most deltas should be 0 or very small
        for qual in &qualities {
            let deltas = encode_with_model(qual, &model);
            let max_delta = deltas.iter().map(|&d| d.abs()).max().unwrap_or(0);
            assert!(max_delta <= 5, "Deltas should be small, got max: {}", max_delta);
        }
    }
}
