//! Zstandard dictionary training for improved compression
//!
//! Quality scores have repeating patterns that can be exploited with
//! dictionary-based compression. We train a dictionary on sample data,
//! then use it for compression.
//!
//! ## Performance Considerations
//!
//! Dictionary training introduces overhead that may outweigh compression gains
//! for smaller datasets:
//!
//! - **Small datasets (<500K reads)**: Dictionary overhead (typically 64KB)
//!   exceeds compression improvement. NOT RECOMMENDED.
//!
//! - **Large datasets (>1M reads)**: Dictionary overhead is amortized over
//!   more data, yielding 5-10% improvement over standard zstd.
//!
//! - **Multi-file scenarios**: Training once and reusing the dictionary across
//!   multiple files with similar patterns can provide significant benefits.
//!
//! - **Streaming compression**: Shared dictionaries enable better compression
//!   when compressing many independent chunks.
//!
//! ## Benchmark Results (50K reads)
//!
//! - Baseline (no dict): 4.23x compression ratio
//! - With dictionary: 4.02x (WORSE due to 64KB overhead)
//! - Quality modeling: 4.54x
//! - Modeling + dict: 4.37x (WORSE due to overhead)
//!
//! **Recommendation**: Use `--quality-modeling` without dictionary training
//! for typical single-file compression tasks.

use anyhow::Result;
use tracing::info;

/// Train a Zstd dictionary from sample data
///
/// # Arguments
/// * `samples` - Vec of byte slices to train on (quality strings)
/// * `dict_size` - Target dictionary size in bytes (typically 16KB-100KB)
///
/// # Returns
/// Trained dictionary bytes
pub fn train_dictionary(samples: &[&[u8]], dict_size: usize) -> Result<Vec<u8>> {
    if samples.is_empty() {
        anyhow::bail!("Cannot train dictionary on empty samples");
    }

    // Collect training samples with size prefixes
    let mut sample_sizes = Vec::new();
    let mut concatenated = Vec::new();

    for sample in samples {
        if !sample.is_empty() {
            sample_sizes.push(sample.len());
            concatenated.extend_from_slice(sample);
        }
    }

    if concatenated.is_empty() {
        anyhow::bail!("No training data available");
    }

    info!("Training zstd dictionary: {} samples, {} bytes total, target size {} bytes",
        sample_sizes.len(), concatenated.len(), dict_size);

    // Train dictionary using zstd's built-in training
    // The zstd crate uses the official zstd library's training algorithm
    let dict = zstd::dict::from_continuous(&concatenated, &sample_sizes, dict_size)
        .map_err(|e| anyhow::anyhow!("Dictionary training failed: {}", e))?;

    info!("Dictionary trained successfully: {} bytes", dict.len());

    Ok(dict)
}

/// Train dictionary from quality scores for optimal compression
///
/// # Arguments
/// * `quality_strings` - Vector of quality score strings
/// * `dict_size` - Target dictionary size (recommended: 16KB-64KB)
/// * `sample_fraction` - Fraction of data to use for training (0.0-1.0)
///
/// # Returns
/// Trained dictionary
pub fn train_from_quality_scores(
    quality_strings: &[Vec<u8>],
    dict_size: usize,
    sample_fraction: f64,
) -> Result<Vec<u8>> {
    if quality_strings.is_empty() {
        anyhow::bail!("No quality strings provided");
    }

    // Determine how many samples to use
    let max_samples = (quality_strings.len() as f64 * sample_fraction).ceil() as usize;
    let max_samples = max_samples.max(100).min(quality_strings.len());

    // Collect samples evenly distributed
    let samples: Vec<&[u8]> = if quality_strings.len() <= max_samples {
        quality_strings.iter().map(|s| s.as_slice()).collect()
    } else {
        let step = quality_strings.len() / max_samples;
        quality_strings
            .iter()
            .step_by(step)
            .take(max_samples)
            .map(|s| s.as_slice())
            .collect()
    };

    train_dictionary(&samples, dict_size)
}

/// Compress data using Zstd with a trained dictionary
pub fn compress_with_dict(data: &[u8], dictionary: &[u8], level: i32) -> Result<Vec<u8>> {
    use std::io::Write;

    let mut encoder = zstd::stream::write::Encoder::with_dictionary(Vec::new(), level, dictionary)
        .map_err(|e| anyhow::anyhow!("Failed to create encoder with dictionary: {}", e))?;

    encoder.write_all(data)
        .map_err(|e| anyhow::anyhow!("Compression failed: {}", e))?;

    encoder.finish()
        .map_err(|e| anyhow::anyhow!("Failed to finish compression: {}", e))
}

/// Decompress data using Zstd with a dictionary
pub fn decompress_with_dict(compressed: &[u8], dictionary: &[u8]) -> Result<Vec<u8>> {
    use std::io::Read;

    let mut decoder = zstd::stream::read::Decoder::with_dictionary(compressed, dictionary)
        .map_err(|e| anyhow::anyhow!("Failed to create decoder with dictionary: {}", e))?;

    let mut decompressed = Vec::new();
    decoder.read_to_end(&mut decompressed)
        .map_err(|e| anyhow::anyhow!("Decompression failed: {}", e))?;

    Ok(decompressed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_train_simple_dictionary() {
        // Create realistic training data (quality strings with repeating patterns)
        // Need at least ~10KB of training data for a 256-byte dictionary
        let mut samples = Vec::new();
        for _ in 0..200 {
            samples.push(b"IIIIIHHHHGGGGFFFFEEEEEDDDDCCCCBBBBAAAAAA".as_ref());
            samples.push(b"HHHHHGGGGFFFFEEEEEDDDDCCCCBBBBAAAAAIIII".as_ref());
            samples.push(b"GGGGFFFFEEEEEDDDDCCCCBBBBAAAAIIIIHHHHHH".as_ref());
        }
        // Total: 600 samples × 40 bytes = 24KB training data

        let dict = train_dictionary(&samples, 256).unwrap();
        assert!(!dict.is_empty());
        assert!(dict.len() <= 256);
    }

    #[test]
    fn test_compress_decompress() {
        let data = b"IIIIIHHHHGGGGFFFFEEEE";
        let dict = vec![b'I', b'H', b'G', b'F', b'E'];

        let compressed = compress_with_dict(data, &dict, 3).unwrap();
        let decompressed = decompress_with_dict(&compressed, &dict).unwrap();

        assert_eq!(data.as_ref(), decompressed.as_slice());
    }
}
