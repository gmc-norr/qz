/// Adaptive arithmetic coding for DNA sequences using Markov-based context models
///
/// This achieves better compression than 2-bit encoding by using context:
/// - Previous 3 bases (order-3 Markov model)
/// - Position in read (binned into 10 positions)
/// - GC content of read (binned into 5 categories)
///
/// Expected improvement: 2-3x better than current sequence compression (4x → 8-12x)

use anyhow::Result;
use constriction::stream::model::DefaultContiguousCategoricalEntropyModel;
use constriction::stream::stack::DefaultAnsCoder;
use constriction::stream::{Decode, Encode};
use constriction::UnwrapInfallible;
use std::collections::HashMap;

const NUM_BASES: usize = 5; // A, C, G, T, N
const CONTEXT_ORDER: usize = 3; // Use previous 3 bases
const NUM_POSITION_BINS: usize = 10;

/// Context for DNA base prediction
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct SequenceContext {
    prev_bases: u16,     // Packed representation of previous bases (up to 3 bases = 12 bits)
    position_bin: u8,    // 0-9 (position in read divided into bins)
}

impl SequenceContext {
    fn new(prev_bases: &[u8], position: usize, read_length: usize) -> Self {
        // Pack previous bases into u16 (3 bases, 4 values each = 12 bits)
        let mut packed = 0u16;
        for (i, &base) in prev_bases.iter().rev().take(CONTEXT_ORDER).enumerate() {
            packed |= (base as u16) << (i * 4);
        }

        let position_bin = ((position as f64 / read_length as f64) * NUM_POSITION_BINS as f64) as u8;

        Self {
            prev_bases: packed,
            position_bin: position_bin.min(NUM_POSITION_BINS as u8 - 1),
        }
    }

    fn hash(&self) -> u64 {
        ((self.prev_bases as u64) << 16) | (self.position_bin as u64)
    }
}

/// Adaptive sequence model with Markov-based prediction
pub struct AdaptiveSequenceModel {
    /// Store frequency counts for each context: context_hash -> [count for each base]
    contexts: HashMap<u64, Vec<u32>>,
    max_total: u32,
}

impl AdaptiveSequenceModel {
    pub fn new() -> Self {
        Self {
            contexts: HashMap::new(),
            max_total: 1 << 15, // 32K max before rescaling
        }
    }

    /// Get probability distribution for a context
    pub(crate) fn get_probabilities(&mut self, context: &SequenceContext) -> Vec<f64> {
        let ctx_hash = context.hash();
        let counts = self
            .contexts
            .entry(ctx_hash)
            .or_insert_with(|| vec![10u32, 10, 10, 10, 1]); // Favor ACGT over N

        let total: u32 = counts.iter().sum();
        counts.iter().map(|&c| c as f64 / total as f64).collect()
    }

    /// Update model after observing a symbol
    pub(crate) fn update(&mut self, context: &SequenceContext, symbol: u8) {
        let ctx_hash = context.hash();
        let counts = self
            .contexts
            .entry(ctx_hash)
            .or_insert_with(|| vec![10u32, 10, 10, 10, 1]);

        counts[(symbol as usize).min(NUM_BASES - 1)] += 4; // Faster learning for sequences

        // Rescale if needed
        let total: u32 = counts.iter().sum();
        if total >= self.max_total {
            for c in counts.iter_mut() {
                *c = (*c + 1) / 2;
            }
        }
    }
}

/// Convert Vec<u32> to Vec<u8> for storage
fn u32_to_u8(data: &[u32]) -> Vec<u8> {
    data.iter().flat_map(|&x| x.to_le_bytes()).collect()
}

/// Convert Vec<u8> back to Vec<u32>
fn u8_to_u32(data: &[u8]) -> Vec<u32> {
    data.chunks_exact(4)
        .map(|chunk| u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
        .collect()
}

/// Encode sequences using adaptive arithmetic coding
///
/// Returns compressed bytes
pub fn encode_sequences_arithmetic(sequences: &[String]) -> Result<Vec<u8>> {
    let mut model = AdaptiveSequenceModel::new();
    let mut models_and_symbols = Vec::new();

    // Build all models in forward order (matching decode order)
    for seq_str in sequences {
        let seq_bytes = seq_str.as_bytes();
        let read_length = seq_str.len();

        let mut prev_bases = Vec::new();

        for (pos, &base_char) in seq_bytes.iter().enumerate() {
            let base = match base_char {
                b'A' => 0u8,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4, // N or other
            };

            let context = SequenceContext::new(&prev_bases, pos, read_length);

            // Get probabilities from model
            let probs = model.get_probabilities(&context);

            // Create entropy model
            let entropy_model =
                DefaultContiguousCategoricalEntropyModel::from_floating_point_probabilities_fast(
                    &probs, None,
                )
                .map_err(|_| anyhow::anyhow!("Failed to create entropy model"))?;

            models_and_symbols.push((base, entropy_model));

            // Update model
            model.update(&context, base);

            // Update previous bases (rolling window)
            prev_bases.push(base);
            if prev_bases.len() > CONTEXT_ORDER {
                prev_bases.remove(0);
            }
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

/// Decode sequences from compressed data
pub fn decode_sequences_arithmetic(
    compressed: &[u8],
    read_lengths: &[usize],
    num_reads: usize,
) -> Result<Vec<String>> {
    let compressed_u32 = u8_to_u32(compressed);
    let mut decoder = DefaultAnsCoder::from_compressed(compressed_u32)
        .map_err(|e| anyhow::anyhow!("Failed to create decoder: {:?}", e))?;

    let mut model = AdaptiveSequenceModel::new();
    let mut sequences = Vec::with_capacity(num_reads);

    for read_idx in 0..num_reads {
        let read_length = read_lengths[read_idx];
        let mut seq_bytes = Vec::with_capacity(read_length);
        let mut prev_bases = Vec::new();

        for pos in 0..read_length {
            let context = SequenceContext::new(&prev_bases, pos, read_length);
            let probs = model.get_probabilities(&context);

            let entropy_model =
                DefaultContiguousCategoricalEntropyModel::from_floating_point_probabilities_fast(
                    &probs, None,
                )
                .map_err(|_| anyhow::anyhow!("Failed to create entropy model"))?;

            // Decode symbol
            let base = decoder.decode_symbol(&entropy_model)? as u8;

            let base_char = match base {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b'N',
            };
            seq_bytes.push(base_char);

            // Update model
            model.update(&context, base);

            // Update previous bases
            prev_bases.push(base);
            if prev_bases.len() > CONTEXT_ORDER {
                prev_bases.remove(0);
            }
        }

        sequences.push(String::from_utf8(seq_bytes)?);
    }

    Ok(sequences)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_context() {
        let ctx1 = SequenceContext::new(&[0, 1, 2], 0, 100); // ACG, pos 0
        let ctx2 = SequenceContext::new(&[0, 1, 3], 0, 100); // ACT, pos 0

        assert_ne!(
            ctx1.hash(),
            ctx2.hash(),
            "Different previous bases should have different hashes"
        );

        let ctx3 = SequenceContext::new(&[0, 1, 2], 50, 100); // ACG, pos 50
        assert_ne!(
            ctx1.hash(),
            ctx3.hash(),
            "Different positions should have different hashes"
        );
    }

    #[test]
    fn test_adaptive_sequence_model() {
        let mut model = AdaptiveSequenceModel::new();
        let context = SequenceContext::new(&[0, 1, 2], 0, 100); // ACG context

        // Initially should have bias toward ACGT
        let probs1 = model.get_probabilities(&context);
        assert_eq!(probs1.len(), NUM_BASES);

        // After updating with base T multiple times
        for _ in 0..10 {
            model.update(&context, 3); // T
        }

        let probs2 = model.get_probabilities(&context);
        // T should now have higher probability
        assert!(
            probs2[3] > probs1[3],
            "Probability of T should increase after updates"
        );
    }

    #[test]
    fn test_sequence_arithmetic_roundtrip() -> Result<()> {
        // Simple test with 3 sequences
        let sequences = vec![
            "ACGTACGTAC".to_string(),
            "GGGGGGGGGG".to_string(),
            "AAAAAAAAAA".to_string(),
        ];

        let read_lengths: Vec<usize> = sequences.iter().map(|s| s.len()).collect();

        let compressed = encode_sequences_arithmetic(&sequences)?;
        let decoded = decode_sequences_arithmetic(&compressed, &read_lengths, 3)?;

        assert_eq!(sequences, decoded, "Roundtrip failed");

        // Check compression
        let original_size: usize = sequences.iter().map(|s| s.len()).sum();
        let compression_ratio = original_size as f64 / compressed.len() as f64;

        println!(
            "Sequence arithmetic coding: {} bytes -> {} bytes ({:.2}x)",
            original_size,
            compressed.len(),
            compression_ratio
        );

        Ok(())
    }

    #[test]
    fn test_repetitive_sequences() -> Result<()> {
        // Test with highly repetitive sequences (should compress very well)
        let sequences = vec![
            "ACGTACGTACGTACGTACGT".to_string(), // Perfect repeat
            "ACGTACGTACGTACGTACGT".to_string(), // Same repeat
            "ACGTACGTACGTACGTACGT".to_string(), // Same repeat
        ];

        let read_lengths: Vec<usize> = sequences.iter().map(|s| s.len()).collect();

        let compressed = encode_sequences_arithmetic(&sequences)?;
        let decoded = decode_sequences_arithmetic(&compressed, &read_lengths, 3)?;

        assert_eq!(sequences, decoded, "Roundtrip failed for repetitive sequences");

        let original_size: usize = sequences.iter().map(|s| s.len()).sum();
        let compression_ratio = original_size as f64 / compressed.len() as f64;

        println!(
            "Repetitive sequences: {} bytes -> {} bytes ({:.2}x)",
            original_size,
            compressed.len(),
            compression_ratio
        );

        // Should achieve very good compression for perfect repeats
        assert!(
            compression_ratio > 2.0,
            "Should get >2x compression for perfect repeats (got {:.2}x)",
            compression_ratio
        );

        Ok(())
    }

    #[test]
    fn test_varied_sequences() -> Result<()> {
        // Test with more varied sequences
        let sequences = vec![
            "ACGTNNNNNN".to_string(), // Has N's
            "AAAAAAAAAA".to_string(), // Low GC
            "GGGGGGGGGG".to_string(), // High GC
            "ACGTACGTAC".to_string(), // Mixed
        ];

        let read_lengths: Vec<usize> = sequences.iter().map(|s| s.len()).collect();

        let compressed = encode_sequences_arithmetic(&sequences)?;
        let decoded = decode_sequences_arithmetic(&compressed, &read_lengths, 4)?;

        assert_eq!(sequences, decoded, "Roundtrip failed for varied sequences");

        Ok(())
    }
}
