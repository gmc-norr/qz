use crate::cli::QualityMode;
use std::borrow::Cow;
use tracing::info;

/// Illumina 8-level binning table
const ILLUMINA_BINS: [u8; 94] = [
    33, 33, 33, 33, 33, 33, 33, 33, 33, 33, // 0-9
    39, 39, 39, 39, 39, 39, 39, 39, 39, 39, // 10-19
    48, 48, 48, 48, 48, 48, 48, 48, 48, 48, // 20-29
    55, 55, 55, 55, 55, 55, 55, 55, 55, 55, // 30-39
    60, 60, 60, 60, 60, 60, 60, 60, 60, 60, // 40-49
    66, 66, 66, 66, 66, 66, 66, 66, 66, 66, // 50-59
    70, 70, 70, 70, 70, 70, 70, 70, 70, 70, // 60-69
    73, 73, 73, 73, 73, 73, 73, 73, 73, 73, // 70-79
    73, 73, 73, 73, 73, 73, 73, 73, 73, 73, // 80-89
    73, 73, 73, 73, // 90-93
];

/// Quantize quality bytes based on mode
pub fn quantize_quality<'a>(quality: &'a [u8], mode: QualityMode) -> Cow<'a, [u8]> {
    match mode {
        QualityMode::Lossless => Cow::Borrowed(quality),
        QualityMode::IlluminaBin => Cow::Owned(quantize_illumina(quality)),
        QualityMode::Illumina4 => Cow::Owned(quantize_four_level(quality)),
        QualityMode::Binary => Cow::Owned(quantize_binary(quality, 20, 40, 6)),
        QualityMode::Qvz => {
            // TODO: Implement QVZ compression
            info!("QVZ mode not yet implemented, using lossless");
            Cow::Borrowed(quality)
        }
        QualityMode::Discard => {
            Cow::Owned(vec![b'!'; quality.len()])
        }
    }
}

/// Apply 4-level binning (more aggressive than Illumina 8-level)
fn quantize_four_level(quality: &[u8]) -> Vec<u8> {
    quality
        .iter()
        .map(|&qv| {
            let phred = (qv as i32).saturating_sub(33).max(0).min(93);
            let binned = if phred < 10 {
                6  // Low quality
            } else if phred < 20 {
                15 // Medium quality
            } else if phred < 30 {
                25 // Good quality
            } else {
                37 // High quality
            };
            (binned + 33) as u8
        })
        .collect()
}

/// Apply Illumina 8-level binning
fn quantize_illumina(quality: &[u8]) -> Vec<u8> {
    quality
        .iter()
        .map(|&qv| {
            let idx = (qv as usize).saturating_sub(33).min(93);
            ILLUMINA_BINS[idx]
        })
        .collect()
}

/// Binary thresholding: high if >= threshold, low otherwise
fn quantize_binary(quality: &[u8], threshold: u8, high: u8, low: u8) -> Vec<u8> {
    quality
        .iter()
        .map(|&qv| {
            let phred = qv.saturating_sub(33);
            if phred >= threshold {
                33 + high
            } else {
                33 + low
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lossless() {
        let quality = b"IIIIIIII";
        assert_eq!(quantize_quality(quality, QualityMode::Lossless).as_ref(), quality);
    }

    #[test]
    fn test_illumina_binning() {
        let quality = b"!!!!!!!!"; // Phred 0
        let binned = quantize_illumina(quality);
        assert_eq!(binned[0], b'!'); // Should map to lowest bin
    }

    #[test]
    fn test_binary_threshold() {
        let quality = b"!!!!IIII"; // Phred 0 and 40
        let binned = quantize_binary(quality, 20, 40, 6);
        assert_eq!(&binned[..4], b"''''");
        assert_eq!(&binned[4..], b"IIII");
    }
}
