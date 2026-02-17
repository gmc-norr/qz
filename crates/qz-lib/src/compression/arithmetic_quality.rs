/// Adaptive arithmetic coding for quality scores using positional context models
///
/// This achieves better compression than simple binning by using context:
/// - Position in read (bin into 20 positions)
/// - Previous quality score
/// - Base call (A/C/G/T/N)
///
/// Expected improvement: 2-3x better than current quality compression

use anyhow::Result;
use constriction::stream::model::DefaultContiguousCategoricalEntropyModel;
use constriction::stream::stack::DefaultAnsCoder;
use constriction::stream::{Decode, Encode};
use constriction::UnwrapInfallible;
use std::collections::HashMap;

const NUM_QUALITY_VALUES: usize = 94; // Phred 0-93
const NUM_POSITION_BINS: usize = 20; // Divide read into 20 position bins
const NUM_BASES: usize = 5; // A, C, G, T, N
const MAX_QUALITY_CONTEXT: u8 = 40; // Bin high qualities together

/// Context for quality score prediction
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct QualityContext {
    position_bin: u8,    // 0-19 (position in read divided into bins)
    prev_quality: u8,    // Previous quality (binned to 0-40)
    base: u8,            // Base call (0=A, 1=C, 2=G, 3=T, 4=N)
}

impl QualityContext {
    fn new(position: usize, read_length: usize, prev_quality: u8, base: u8) -> Self {
        let position_bin = ((position as f64 / read_length as f64) * NUM_POSITION_BINS as f64) as u8;
        Self {
            position_bin: position_bin.min(NUM_POSITION_BINS as u8 - 1),
            prev_quality: prev_quality.min(MAX_QUALITY_CONTEXT),
            base: base.min(NUM_BASES as u8 - 1),
        }
    }

    fn hash(&self) -> u64 {
        ((self.position_bin as u64) << 16)
            | ((self.prev_quality as u64) << 8)
            | (self.base as u64)
    }
}

/// Adaptive quality model with context-based prediction
pub struct AdaptiveQualityModel {
    /// Store frequency counts for each context: context_hash -> [count for each quality value]
    contexts: HashMap<u64, Vec<u32>>,
    max_total: u32,
}

impl AdaptiveQualityModel {
    pub fn new() -> Self {
        Self {
            contexts: HashMap::new(),
            max_total: 1 << 15, // 32K max before rescaling
        }
    }

    /// Get probability distribution for a context
    pub(crate) fn get_probabilities(&mut self, context: &QualityContext) -> Vec<f64> {
        let ctx_hash = context.hash();
        let counts = self
            .contexts
            .entry(ctx_hash)
            .or_insert_with(|| vec![1u32; NUM_QUALITY_VALUES]); // Laplace smoothing

        let total: u32 = counts.iter().sum();
        counts.iter().map(|&c| c as f64 / total as f64).collect()
    }

    /// Update model after observing a symbol
    pub(crate) fn update(&mut self, context: &QualityContext, symbol: u8) {
        let ctx_hash = context.hash();
        let counts = self
            .contexts
            .entry(ctx_hash)
            .or_insert_with(|| vec![1u32; NUM_QUALITY_VALUES]);

        counts[(symbol as usize).min(NUM_QUALITY_VALUES - 1)] += 1;

        // Rescale if needed to prevent overflow
        let total: u32 = counts.iter().sum();
        if total >= self.max_total {
            for c in counts.iter_mut() {
                *c = (*c + 1) / 2; // Halve all counts
            }
        }
    }
}

/// Convert Vec<u32> to Vec<u8> for storage
fn u32_to_u8(data: &[u32]) -> Vec<u8> {
    data.iter()
        .flat_map(|&x| x.to_le_bytes())
        .collect()
}

/// Convert Vec<u8> back to Vec<u32>
fn u8_to_u32(data: &[u8]) -> Vec<u32> {
    data.chunks_exact(4)
        .map(|chunk| u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
        .collect()
}

/// Encode quality scores using adaptive arithmetic coding
///
/// Returns compressed bytes
pub fn encode_qualities_arithmetic(
    qualities: &[String],
    sequences: &[String],
    read_length: usize,
) -> Result<Vec<u8>> {
    if qualities.len() != sequences.len() {
        anyhow::bail!(
            "Qualities and sequences length mismatch: {} != {}",
            qualities.len(),
            sequences.len()
        );
    }

    let mut model = AdaptiveQualityModel::new();
    let mut models_and_symbols = Vec::new();

    // Build all models in forward order (matching decode order)
    for (read_idx, qual_str) in qualities.iter().enumerate() {
        let qual_bytes = qual_str.as_bytes();
        let seq_bytes = sequences[read_idx].as_bytes();

        if qual_bytes.len() != seq_bytes.len() {
            anyhow::bail!(
                "Read {} quality length {} != sequence length {}",
                read_idx,
                qual_bytes.len(),
                seq_bytes.len()
            );
        }

        let mut prev_quality = 0u8; // Start with Phred 0

        for (pos, &qual_char) in qual_bytes.iter().enumerate() {
            let quality = qual_char.saturating_sub(33); // Phred score (clamp to 0 if < 33)
            let base = match seq_bytes[pos] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4, // N or other
            };

            let context = QualityContext::new(pos, read_length, prev_quality, base);

            // Get probabilities from model
            let probs = model.get_probabilities(&context);

            // Create entropy model
            let entropy_model =
                DefaultContiguousCategoricalEntropyModel::from_floating_point_probabilities_fast(
                    &probs, None,
                )
                .map_err(|_| anyhow::anyhow!("Failed to create entropy model"))?;

            models_and_symbols.push((quality, entropy_model));

            // Update model
            model.update(&context, quality);
            prev_quality = quality;
        }
    }

    // Encode all symbols in reverse (ANS requirement)
    let mut encoder = DefaultAnsCoder::new();
    for (symbol, model) in models_and_symbols.iter().rev() {
        encoder.encode_symbol(*symbol as usize, model)?;
    }

    let compressed_u32 = encoder.into_compressed().unwrap_infallible();
    Ok(u32_to_u8(&compressed_u32))
}

/// Decode quality scores from compressed data
pub fn decode_qualities_arithmetic(
    compressed: &[u8],
    sequences: &[String],
    read_length: usize,
    num_reads: usize,
) -> Result<Vec<String>> {
    let compressed_u32 = u8_to_u32(compressed);
    let mut decoder = DefaultAnsCoder::from_compressed(compressed_u32)
        .map_err(|e| anyhow::anyhow!("Failed to create decoder: {:?}", e))?;
    let mut model = AdaptiveQualityModel::new();
    let mut qualities = Vec::with_capacity(num_reads);

    for seq_str in sequences.iter().take(num_reads) {
        let seq_bytes = seq_str.as_bytes();
        let mut qual_bytes = Vec::with_capacity(read_length);
        let mut prev_quality = 0u8;

        for pos in 0..seq_bytes.len().min(read_length) {
            let base = match seq_bytes[pos] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            };

            let context = QualityContext::new(pos, read_length, prev_quality, base);
            let probs = model.get_probabilities(&context);

            let entropy_model =
                DefaultContiguousCategoricalEntropyModel::from_floating_point_probabilities_fast(
                    &probs, None,
                )
                .map_err(|_| anyhow::anyhow!("Failed to create entropy model"))?;

            // Decode symbol
            let quality = decoder.decode_symbol(&entropy_model)? as u8;
            qual_bytes.push(quality + 33); // Convert back to ASCII

            // Update model
            model.update(&context, quality);
            prev_quality = quality;
        }

        qualities.push(String::from_utf8(qual_bytes)?);
    }

    Ok(qualities)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_context() {
        let ctx1 = QualityContext::new(0, 100, 30, 0); // Position 0, Phred 30, base A
        let ctx2 = QualityContext::new(50, 100, 30, 0); // Position 50, Phred 30, base A

        assert_ne!(ctx1.hash(), ctx2.hash(), "Different positions should have different hashes");

        let ctx3 = QualityContext::new(0, 100, 35, 0); // Position 0, Phred 35, base A
        assert_ne!(ctx1.hash(), ctx3.hash(), "Different prev qualities should have different hashes");

        let ctx4 = QualityContext::new(0, 100, 30, 1); // Position 0, Phred 30, base C
        assert_ne!(ctx1.hash(), ctx4.hash(), "Different bases should have different hashes");
    }

    #[test]
    fn test_adaptive_quality_model() {
        let mut model = AdaptiveQualityModel::new();
        let context = QualityContext::new(0, 100, 30, 0);

        // Initially should have uniform-ish distribution
        let probs1 = model.get_probabilities(&context);
        assert_eq!(probs1.len(), NUM_QUALITY_VALUES);

        // After updating with quality 40 multiple times
        for _ in 0..10 {
            model.update(&context, 40);
        }

        let probs2 = model.get_probabilities(&context);
        // Quality 40 should now have higher probability
        assert!(
            probs2[40] > probs1[40],
            "Probability of quality 40 should increase after updates"
        );
    }

    #[test]
    fn test_quality_arithmetic_roundtrip() -> Result<()> {
        // Simple test with 3 reads
        let qualities = vec![
            "IIIIIIIIII".to_string(), // Phred 40 (I = 73, 73-33=40)
            "HHHHHIIIII".to_string(), // Mix of 39 and 40
            "GGGGGHHHHI".to_string(), // Mix of 38, 39, 40
        ];
        let sequences = vec![
            "ACGTACGTAC".to_string(),
            "ACGTACGTAC".to_string(),
            "ACGTACGTAC".to_string(),
        ];

        let compressed = encode_qualities_arithmetic(&qualities, &sequences, 10)?;
        let decoded = decode_qualities_arithmetic(&compressed, &sequences, 10, 3)?;

        assert_eq!(qualities, decoded, "Roundtrip failed");

        // Check compression
        let original_size: usize = qualities.iter().map(|s| s.len()).sum();
        let compression_ratio = original_size as f64 / compressed.len() as f64;

        println!(
            "Quality arithmetic coding: {} bytes -> {} bytes ({:.2}x)",
            original_size,
            compressed.len(),
            compression_ratio
        );

        // Should achieve some compression
        assert!(
            compressed.len() < original_size,
            "Should compress better than original"
        );

        Ok(())
    }

    #[test]
    fn test_varied_qualities() -> Result<()> {
        // Test with more varied quality scores
        let qualities = vec![
            "!\"#$%&'()*".to_string(), // Low qualities (Phred 0-9)
            "0123456789".to_string(),  // Medium qualities (Phred 15-24)
            "ABCDEFGHIJ".to_string(),  // High qualities (Phred 32-41)
        ];
        let sequences = vec![
            "AAAAAAAAAA".to_string(),
            "CCCCCCCCCC".to_string(),
            "GGGGGGGGGG".to_string(),
        ];

        let compressed = encode_qualities_arithmetic(&qualities, &sequences, 10)?;
        let decoded = decode_qualities_arithmetic(&compressed, &sequences, 10, 3)?;

        assert_eq!(qualities, decoded, "Roundtrip failed for varied qualities");

        Ok(())
    }

    #[test]
    fn test_position_dependent_quality() -> Result<()> {
        // Test that position-dependent quality patterns compress well
        // Quality drops at the end (typical Illumina pattern)
        let qualities = vec![
            "IIIIIHHHGG".to_string(), // 40,40,40,40,40,39,39,39,38,38
            "IIIIIHHHGG".to_string(), // Same pattern
            "IIIIIHHHGG".to_string(), // Same pattern
        ];
        let sequences = vec![
            "ACGTACGTAC".to_string(),
            "ACGTACGTAC".to_string(),
            "ACGTACGTAC".to_string(),
        ];

        let compressed = encode_qualities_arithmetic(&qualities, &sequences, 10)?;
        let decoded = decode_qualities_arithmetic(&compressed, &sequences, 10, 3)?;

        assert_eq!(qualities, decoded, "Roundtrip failed");

        let original_size: usize = qualities.iter().map(|s| s.len()).sum();
        let compression_ratio = original_size as f64 / compressed.len() as f64;

        println!(
            "Position-dependent quality: {} bytes -> {} bytes ({:.2}x)",
            original_size,
            compressed.len(),
            compression_ratio
        );

        // Should achieve some compression due to repeated pattern
        // Note: With only 30 bytes, we can't expect huge compression ratios
        assert!(
            compression_ratio > 1.0,
            "Should achieve some compression for repeated patterns (got {:.2}x)",
            compression_ratio
        );

        Ok(())
    }
}
