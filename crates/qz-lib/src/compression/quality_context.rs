/// Context-modeled quality score compression using adaptive arithmetic coding.
///
/// Context = (position_in_read quantized, previous_quality quantized)
/// Each context has its own adaptive frequency model.
/// This exploits the strong positional and local-correlation structure
/// in Illumina quality scores.

use anyhow::Result;

// ── Adaptive frequency model ────────────────────────────────────────────────

const NSYM: usize = 96; // Quality symbols: ASCII 33-128 mapped to 0-95

/// Adaptive frequency model for a single context.
struct AdaptiveModel {
    freq: [u16; NSYM],
    total: u32,
}

impl AdaptiveModel {
    fn new() -> Self {
        let mut m = Self {
            freq: [0; NSYM],
            total: 0,
        };
        // Initialize with count=1 for common quality range only (ASCII 33-73 → 0-40)
        for i in 0..NSYM {
            m.freq[i] = 1;
        }
        m.total = NSYM as u32;
        m
    }

    #[inline]
    fn encode_symbol(&self, symbol: u8) -> (u32, u32, u32) {
        let s = symbol as usize;
        let mut cum = 0u32;
        for i in 0..s {
            cum += self.freq[i] as u32;
        }
        (cum, cum + self.freq[s] as u32, self.total)
    }

    #[inline]
    fn update(&mut self, symbol: u8) {
        self.freq[symbol as usize] += 1;
        self.total += 1;
        if self.total >= 0x7FFF {
            self.total = 0;
            for f in self.freq.iter_mut() {
                *f = (*f >> 1).max(1);
                self.total += *f as u32;
            }
        }
    }

    #[inline]
    fn decode_symbol(&self, target: u32) -> (u8, u32, u32) {
        let mut cum = 0u32;
        for i in 0..NSYM {
            let next = cum + self.freq[i] as u32;
            if target < next {
                return (i as u8, cum, next);
            }
            cum = next;
        }
        ((NSYM - 1) as u8, cum - self.freq[NSYM - 1] as u32, cum)
    }
}

// ── Carryless range coder (Subbotin-style) ──────────────────────────────────

const RANGE_TOP: u64 = 1u64 << 32;
const RANGE_BOT: u64 = 1u64 << 24;

struct RangeEncoder {
    low: u64,
    range: u64,
    cache: u8,
    cache_size: u32,
    output: Vec<u8>,
    first_byte: bool,
}

impl RangeEncoder {
    fn new() -> Self {
        Self {
            low: 0,
            range: RANGE_TOP,
            cache: 0,
            cache_size: 0,
            output: Vec::new(),
            first_byte: true,
        }
    }

    fn shift_low(&mut self) {
        let low_hi = (self.low >> 32) as u8; // carry bit
        if self.low < 0xFF000000u64 || low_hi != 0 {
            if !self.first_byte {
                self.output.push(self.cache.wrapping_add(low_hi));
            }
            self.first_byte = false;
            for _ in 0..self.cache_size {
                self.output.push(0xFFu8.wrapping_add(low_hi));
            }
            self.cache_size = 0;
            self.cache = ((self.low >> 24) & 0xFF) as u8;
        } else {
            self.cache_size += 1;
        }
        self.low = (self.low << 8) & 0xFFFFFFFF;
    }

    #[inline]
    fn encode(&mut self, cum_low: u32, cum_high: u32, total: u32) {
        let r = self.range / total as u64;
        if cum_high < total {
            self.range = r * (cum_high - cum_low) as u64;
        } else {
            self.range -= r * cum_low as u64;
        }
        self.low += r * cum_low as u64;

        while self.range < RANGE_BOT {
            self.range <<= 8;
            self.shift_low();
        }
    }

    fn finish(mut self) -> Vec<u8> {
        for _ in 0..5 {
            self.shift_low();
        }
        self.output
    }
}

struct RangeDecoder<'a> {
    low: u64,
    range: u64,
    code: u64,
    input: &'a [u8],
    pos: usize,
}

impl<'a> RangeDecoder<'a> {
    fn new(input: &'a [u8]) -> Self {
        let mut dec = Self {
            low: 0,
            range: RANGE_TOP,
            code: 0,
            input,
            pos: 0,
        };
        for _ in 0..4 {
            dec.code = (dec.code << 8) | dec.read_byte() as u64;
        }
        dec
    }

    #[inline]
    fn read_byte(&mut self) -> u8 {
        if self.pos < self.input.len() {
            let b = self.input[self.pos];
            self.pos += 1;
            b
        } else {
            0
        }
    }

    #[inline]
    fn get_freq(&self, total: u32) -> u32 {
        ((self.code - self.low) / (self.range / total as u64)) as u32
    }

    #[inline]
    fn decode(&mut self, cum_low: u32, cum_high: u32, total: u32) {
        let r = self.range / total as u64;
        self.low += r * cum_low as u64;
        if cum_high < total {
            self.range = r * (cum_high - cum_low) as u64;
        } else {
            self.range -= r * cum_low as u64;
        }

        while self.range < RANGE_BOT {
            self.code = (self.code << 8) | self.read_byte() as u64;
            self.range <<= 8;
            self.low <<= 8;
        }
    }
}

// ── Context computation ─────────────────────────────────────────────────────

#[inline]
fn context_index(pos: usize, prev_qual: u8, pos_bits: u32, qual_bits: u32) -> usize {
    let max_pos = (1usize << pos_bits) - 1;
    let pos_quant = if pos < 8 {
        pos.min(max_pos)
    } else {
        (8 + (pos >> 3)).min(max_pos)
    };

    let qual_quant = (prev_qual >> (8 - qual_bits)) as usize;
    let qual_mask = (1usize << qual_bits) - 1;

    (pos_quant << qual_bits) | (qual_quant & qual_mask)
}

// ── Public API ──────────────────────────────────────────────────────────────

pub struct QualityContextConfig {
    pub pos_bits: u32,
    pub qual_bits: u32,
}

impl Default for QualityContextConfig {
    fn default() -> Self {
        Self {
            pos_bits: 6,
            qual_bits: 4,
        }
    }
}

impl QualityContextConfig {
    fn num_contexts(&self) -> usize {
        (1 << self.pos_bits) * (1 << self.qual_bits)
    }
}

/// Compress quality scores using context-modeled adaptive arithmetic coding.
pub fn compress(qualities: &[u8], read_lengths: &[u32], config: &QualityContextConfig) -> Result<Vec<u8>> {
    if qualities.is_empty() {
        return Ok(Vec::new());
    }

    let num_contexts = config.num_contexts();
    let mut models: Vec<AdaptiveModel> = (0..num_contexts).map(|_| AdaptiveModel::new()).collect();
    let mut encoder = RangeEncoder::new();

    let mut offset = 0usize;
    for &rlen in read_lengths {
        let len = rlen as usize;
        let end = offset + len;
        if end > qualities.len() { break; }

        let mut prev_qual = 0u8;
        for pos in 0..len {
            let q = qualities[offset + pos];
            // Map ASCII to 0-based symbol
            let sym = q.saturating_sub(33);
            let ctx = context_index(pos, prev_qual, config.pos_bits, config.qual_bits);
            let (cum_low, cum_high, total) = models[ctx].encode_symbol(sym);
            encoder.encode(cum_low, cum_high, total);
            models[ctx].update(sym);
            prev_qual = q;
        }
        offset = end;
    }

    Ok(encoder.finish())
}

/// Decompress quality scores.
pub fn decompress(
    compressed: &[u8],
    read_lengths: &[u32],
    config: &QualityContextConfig,
) -> Result<Vec<u8>> {
    if compressed.is_empty() {
        return Ok(Vec::new());
    }

    let total_len: usize = read_lengths.iter().map(|&l| l as usize).sum();
    let num_contexts = config.num_contexts();
    let mut models: Vec<AdaptiveModel> = (0..num_contexts).map(|_| AdaptiveModel::new()).collect();
    let mut decoder = RangeDecoder::new(compressed);
    let mut output = Vec::with_capacity(total_len);

    for &rlen in read_lengths {
        let len = rlen as usize;
        let mut prev_qual = 0u8;
        for pos in 0..len {
            let ctx = context_index(pos, prev_qual, config.pos_bits, config.qual_bits);
            let target = decoder.get_freq(models[ctx].total);
            let (sym, cum_low, cum_high) = models[ctx].decode_symbol(target);
            decoder.decode(cum_low, cum_high, models[ctx].total);
            models[ctx].update(sym);
            let q = sym + 33; // Map back to ASCII
            output.push(q);
            prev_qual = q;
        }
    }

    Ok(output)
}

/// Compress quality scores in parallel blocks.
pub fn compress_parallel(
    qualities: &[u8],
    read_lengths: &[u32],
    config: &QualityContextConfig,
) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    if qualities.is_empty() {
        return Ok(Vec::new());
    }

    let block_size = 100_000;
    let num_blocks = (read_lengths.len() + block_size - 1) / block_size;

    let mut block_starts: Vec<(usize, usize)> = Vec::with_capacity(num_blocks);
    let mut byte_off = 0usize;
    for block_idx in 0..num_blocks {
        let read_start = block_idx * block_size;
        block_starts.push((read_start, byte_off));
        let read_end = ((block_idx + 1) * block_size).min(read_lengths.len());
        for i in read_start..read_end {
            byte_off += read_lengths[i] as usize;
        }
    }

    let blocks: Vec<Vec<u8>> = block_starts
        .par_iter()
        .enumerate()
        .map(|(block_idx, &(_read_start, byte_start))| {
            let read_start = block_idx * block_size;
            let read_end = ((block_idx + 1) * block_size).min(read_lengths.len());
            let byte_end = if block_idx + 1 < block_starts.len() {
                block_starts[block_idx + 1].1
            } else {
                qualities.len()
            };
            compress(&qualities[byte_start..byte_end], &read_lengths[read_start..read_end], config)
                .unwrap_or_default()
        })
        .collect();

    let mut output = Vec::new();
    output.extend_from_slice(&(num_blocks as u32).to_le_bytes());
    for block in &blocks {
        output.extend_from_slice(&(block.len() as u32).to_le_bytes());
    }
    for block in &blocks {
        output.extend_from_slice(block);
    }

    Ok(output)
}

/// Decompress parallel-compressed quality scores.
pub fn decompress_parallel(
    compressed: &[u8],
    read_lengths: &[u32],
    config: &QualityContextConfig,
) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    if compressed.is_empty() {
        return Ok(Vec::new());
    }

    let block_size = 100_000;

    let num_blocks = super::read_le_u32(compressed, 0)? as usize;
    let mut offset = 4;
    let mut block_lens = Vec::with_capacity(num_blocks);
    for _ in 0..num_blocks {
        let blen = super::read_le_u32(compressed, offset)? as usize;
        block_lens.push(blen);
        offset += 4;
    }

    let mut block_data: Vec<&[u8]> = Vec::with_capacity(num_blocks);
    for &blen in &block_lens {
        block_data.push(&compressed[offset..offset + blen]);
        offset += blen;
    }

    let decoded_blocks: Vec<Vec<u8>> = (0..num_blocks)
        .into_par_iter()
        .map(|block_idx| {
            let read_start = block_idx * block_size;
            let read_end = ((block_idx + 1) * block_size).min(read_lengths.len());
            decompress(block_data[block_idx], &read_lengths[read_start..read_end], config)
                .unwrap_or_default()
        })
        .collect();

    let total_len: usize = decoded_blocks.iter().map(|b| b.len()).sum();
    let mut output = Vec::with_capacity(total_len);
    for block in decoded_blocks {
        output.extend_from_slice(&block);
    }

    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip_basic() {
        let qualities = b"IIIIIIIIIIFFFFFFFFFFF!!!!ABCDEFGHIJ";
        let read_lengths = vec![10u32, 11, 4, 10];
        let config = QualityContextConfig::default();

        let compressed = compress(qualities, &read_lengths, &config).unwrap();
        let decompressed = decompress(&compressed, &read_lengths, &config).unwrap();
        assert_eq!(decompressed, qualities.to_vec());
    }

    #[test]
    fn test_roundtrip_parallel() {
        let qual_str = b"IIIIIIIIIIFFFFFFFFFFF!!!!ABCDEFGHIJ";
        let read_lengths = vec![10u32, 11, 4, 10];
        let config = QualityContextConfig::default();

        let compressed = compress_parallel(qual_str, &read_lengths, &config).unwrap();
        let decompressed = decompress_parallel(&compressed, &read_lengths, &config).unwrap();
        assert_eq!(decompressed, qual_str.to_vec());
    }

    #[test]
    fn test_roundtrip_realistic() {
        // Simulate typical Illumina quality pattern
        let mut qualities = Vec::new();
        let mut read_lengths = Vec::new();
        for _ in 0..100 {
            let mut read_qual = Vec::new();
            for pos in 0..150 {
                // Typical pattern: high quality in middle, lower at edges
                let q = if pos < 5 { 35 } else if pos > 140 { 30 } else { 37 + (pos % 5) as u8 };
                read_qual.push(q + 33); // ASCII offset
            }
            qualities.extend_from_slice(&read_qual);
            read_lengths.push(150);
        }

        let config = QualityContextConfig::default();
        let compressed = compress(&qualities, &read_lengths, &config).unwrap();
        let decompressed = decompress(&compressed, &read_lengths, &config).unwrap();
        assert_eq!(decompressed, qualities);
    }
}
