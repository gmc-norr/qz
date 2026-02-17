/// Context-aware quality compression using adaptive range coding.
///
/// Context: (position_bin5, prev_quality, quality_stable, prev_base, cur_base)
/// ~4K contexts, converges within ~5K reads.
///
/// Custom forward range coder with integer cumulative frequencies — eliminates
/// the per-symbol entropy model creation overhead from the constriction approach.
/// Quality delta bit (stable vs changing) inspired by ENANO.
/// Achieves ~0.151 bits/base vs BSC's ~0.162 bits/base (-7.2% improvement).
use anyhow::Result;

// ============================================================================
// Range coder constants
// ============================================================================
const RC_TOP: u32 = 1 << 24;

// ============================================================================
// Context parameters
// ============================================================================
const POS_BIN_SIZE: usize = 5;
const MAX_POS_BINS: usize = 64; // supports reads up to 320bp
const MAX_QUAL: usize = 50; // max Phred score we handle
const N_BASE_CTX: usize = 25; // 5 prev_base * 5 cur_base (including sentinel)
const N_DELTA: usize = 2; // 0 = stable (|q_prev - q_prev2| <= 2), 1 = changing
const MAX_CONTEXTS: usize = MAX_POS_BINS * MAX_QUAL * N_DELTA * N_BASE_CTX; // 160000
const RESCALE_THRESHOLD: u32 = 1 << 20;

// ============================================================================
// Range Encoder (LZMA-style with carry propagation via 64-bit low)
// ============================================================================
struct RangeEncoder {
    low: u64,
    range: u32,
    cache: u8,
    cache_size: u32,
    output: Vec<u8>,
}

impl RangeEncoder {
    fn new() -> Self {
        Self {
            low: 0,
            range: 0xFFFF_FFFF,
            cache: 0,
            cache_size: 1,
            output: Vec::new(),
        }
    }

    #[inline(always)]
    fn shift_low(&mut self) {
        let low_hi = (self.low >> 32) as u8;
        if low_hi != 0 || (self.low as u32) < 0xFF00_0000 {
            let mut byte = self.cache;
            loop {
                self.output.push(byte.wrapping_add(low_hi));
                byte = 0xFF;
                self.cache_size -= 1;
                if self.cache_size == 0 {
                    break;
                }
            }
            self.cache = ((self.low >> 24) & 0xFF) as u8;
        }
        self.cache_size += 1;
        self.low = ((self.low as u32) << 8) as u64;
    }

    #[inline(always)]
    fn encode(&mut self, cum: u32, freq: u32, total: u32) {
        let r = self.range / total;
        self.low += cum as u64 * r as u64;
        if cum + freq < total {
            self.range = r * freq;
        } else {
            self.range -= r * cum;
        }
        while self.range < RC_TOP {
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

// ============================================================================
// Range Decoder
// ============================================================================
struct RangeDecoder<'a> {
    range: u32,
    code: u32,
    input: &'a [u8],
    pos: usize,
}

impl<'a> RangeDecoder<'a> {
    fn new(input: &'a [u8]) -> Self {
        let mut dec = Self {
            range: 0xFFFF_FFFF,
            code: 0,
            input,
            pos: 0,
        };
        // Skip first byte (always 0x00 from encoder's initial cache)
        if !input.is_empty() {
            dec.pos = 1;
        }
        // Read 4 bytes into code
        for _ in 0..4 {
            dec.code = (dec.code << 8) | dec.next_byte() as u32;
        }
        dec
    }

    #[inline(always)]
    fn next_byte(&mut self) -> u8 {
        if self.pos < self.input.len() {
            let b = self.input[self.pos];
            self.pos += 1;
            b
        } else {
            0
        }
    }

    #[inline(always)]
    fn normalize(&mut self) {
        while self.range < RC_TOP {
            self.code = (self.code << 8) | self.next_byte() as u32;
            self.range <<= 8;
        }
    }

    #[inline(always)]
    fn decode(&mut self, cum_freqs: &[u32], n_symbols: usize, total: u32) -> usize {
        let r = self.range / total;
        let offset = (self.code / r).min(total - 1);

        // Linear scan for symbol (fast for small alphabets ~20)
        let mut sym = 0;
        while sym + 1 < n_symbols && cum_freqs[sym + 1] <= offset {
            sym += 1;
        }

        let cum = cum_freqs[sym];
        let freq = cum_freqs[sym + 1] - cum;

        self.code -= cum * r;
        if cum + freq < total {
            self.range = r * freq;
        } else {
            self.range -= r * cum;
        }

        self.normalize();
        sym
    }
}

// ============================================================================
// Adaptive Model with integer cumulative frequencies
// ============================================================================
struct AdaptiveModel {
    cum_freqs: Vec<u32>, // length = n_symbols + 1, cum_freqs[0] = 0
    n_symbols: usize,
    total: u32,
}

impl AdaptiveModel {
    fn new(n_symbols: usize) -> Self {
        // Uniform init (Laplace smoothing): each symbol has count 1
        let cum_freqs: Vec<u32> = (0..=n_symbols).map(|i| i as u32).collect();
        Self {
            cum_freqs,
            n_symbols,
            total: n_symbols as u32,
        }
    }

    #[inline(always)]
    fn encode_params(&self, sym: usize) -> (u32, u32, u32) {
        let cum = self.cum_freqs[sym];
        let freq = self.cum_freqs[sym + 1] - cum;
        (cum, freq, self.total)
    }

    #[inline(always)]
    fn update(&mut self, sym: usize) {
        for i in (sym + 1)..=self.n_symbols {
            self.cum_freqs[i] += 1;
        }
        self.total += 1;
        if self.total >= RESCALE_THRESHOLD {
            self.rescale();
        }
    }

    fn rescale(&mut self) {
        // Halve all frequencies, maintaining minimum of 1
        let mut cum = 0u32;
        self.cum_freqs[0] = 0;
        for i in 0..self.n_symbols {
            let freq = self.cum_freqs[i + 1] - self.cum_freqs[i];
            let new_freq = (freq >> 1).max(1);
            cum += new_freq;
            self.cum_freqs[i + 1] = cum;
        }
        self.total = cum;
    }
}

// ============================================================================
// Context computation
// ============================================================================

#[inline(always)]
fn context_id(pos: usize, prev_q: u8, prev_q2: u8, prev_base: u8, cur_base: u8) -> usize {
    let pos_bin = (pos / POS_BIN_SIZE).min(MAX_POS_BINS - 1);
    let pq = (prev_q as usize).min(MAX_QUAL - 1);
    let db = if (prev_q as i16 - prev_q2 as i16).unsigned_abs() <= 2 { 0 } else { 1 };
    let base_ctx = (prev_base as usize).min(4) * 5 + (cur_base as usize).min(4);
    let bc = base_ctx.min(N_BASE_CTX - 1);
    pos_bin * (MAX_QUAL * N_DELTA * N_BASE_CTX)
        + pq * (N_DELTA * N_BASE_CTX)
        + db * N_BASE_CTX
        + bc
}

/// Lookup table for base → index (branchless, called ~300M times).
static BASE_TO_IDX_LUT: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0; t[b'a' as usize] = 0;
    t[b'C' as usize] = 1; t[b'c' as usize] = 1;
    t[b'G' as usize] = 2; t[b'g' as usize] = 2;
    t[b'T' as usize] = 3; t[b't' as usize] = 3;
    t
};

#[inline(always)]
fn base_to_idx(b: u8) -> u8 {
    BASE_TO_IDX_LUT[b as usize]
}

// ============================================================================
// Compress
// ============================================================================

/// Compress quality scores using context-adaptive range coding.
///
/// Format: [n_symbols:1B][symbols:nB][read_len:2B][num_reads:4B][data_len:4B][data]
pub fn compress_qualities_ctx(qualities: &[&[u8]], sequences: &[&[u8]]) -> Result<Vec<u8>> {
    if qualities.is_empty() {
        return Ok(Vec::new());
    }

    // Build symbol alphabet from input qualities
    let mut seen = [false; 256];
    for q in qualities {
        for &b in *q {
            seen[(b - 33) as usize] = true;
        }
    }
    let symbols: Vec<u8> = (0u8..=255).filter(|&i| seen[i as usize]).collect();
    let n_symbols = symbols.len();
    let mut sym_map = [0u8; 256];
    for (i, &s) in symbols.iter().enumerate() {
        sym_map[s as usize] = i as u8;
    }

    // Detect variable-length reads; store read_len=0 as sentinel
    let first_len = qualities[0].len();
    let variable = qualities.iter().any(|q| q.len() != first_len);
    let read_len: u16 = if variable { 0 } else { first_len as u16 };
    let num_reads = qualities.len();

    // Context model lookup: flat Vec, u32::MAX = uninitialized
    let mut ctx_slots: Vec<u32> = vec![u32::MAX; MAX_CONTEXTS];
    let mut models: Vec<AdaptiveModel> = Vec::with_capacity(4096);

    let mut encoder = RangeEncoder::new();

    for (qual, seq) in qualities.iter().zip(sequences.iter()) {
        let qb = *qual;
        let sb = *seq;
        let mut prev_q: u8 = 0;
        let mut prev_q2: u8 = 0;

        for j in 0..qb.len() {
            let phred = qb[j] - 33;
            let sym = sym_map[phred as usize] as usize;
            let cur_base = base_to_idx(sb[j]);
            let prev_base = if j > 0 { base_to_idx(sb[j - 1]) } else { 4 };
            let ctx = context_id(j, prev_q, prev_q2, prev_base, cur_base);

            let slot = if ctx_slots[ctx] == u32::MAX {
                let id = models.len() as u32;
                ctx_slots[ctx] = id;
                models.push(AdaptiveModel::new(n_symbols));
                id as usize
            } else {
                ctx_slots[ctx] as usize
            };

            let (cum, freq, total) = models[slot].encode_params(sym);
            encoder.encode(cum, freq, total);
            models[slot].update(sym);
            prev_q2 = prev_q;
            prev_q = phred;
        }
    }

    let compressed = encoder.finish();

    // Build output blob
    let mut result = Vec::new();
    result.push(n_symbols as u8);
    result.extend_from_slice(&symbols);
    result.extend_from_slice(&read_len.to_le_bytes());
    result.extend_from_slice(&(num_reads as u32).to_le_bytes());
    result.extend_from_slice(&(compressed.len() as u32).to_le_bytes());
    result.extend_from_slice(&compressed);
    Ok(result)
}

// ============================================================================
// Decompress
// ============================================================================

/// Decompress quality scores.
///
/// Format: [n_symbols:1B][symbols:nB][read_len:2B][num_reads:4B][data_len:4B][data]
pub fn decompress_qualities_ctx(
    data: &[u8],
    sequences: &[&[u8]],
    num_reads: usize,
) -> Result<Vec<Vec<u8>>> {
    if data.is_empty() || num_reads == 0 {
        return Ok(Vec::new());
    }

    let mut pos = 0;

    let n_symbols = data[pos] as usize;
    pos += 1;
    let symbols: Vec<u8> = data[pos..pos + n_symbols].to_vec();
    pos += n_symbols;

    let read_len = data[pos] as usize | ((data[pos + 1] as usize) << 8);
    pos += 2;
    let _stored_reads =
        u32::from_le_bytes([data[pos], data[pos + 1], data[pos + 2], data[pos + 3]]);
    pos += 4;
    let compressed_len =
        u32::from_le_bytes([data[pos], data[pos + 1], data[pos + 2], data[pos + 3]]) as usize;
    pos += 4;

    let compressed_data = &data[pos..pos + compressed_len];
    let mut decoder = RangeDecoder::new(compressed_data);

    let mut ctx_slots: Vec<u32> = vec![u32::MAX; MAX_CONTEXTS];
    let mut models: Vec<AdaptiveModel> = Vec::with_capacity(4096);
    let mut result = Vec::with_capacity(num_reads);

    for i in 0..num_reads {
        let sb = sequences[i];
        let this_len = if read_len > 0 { read_len } else { sb.len() };
        let mut qual_bytes = Vec::with_capacity(this_len);
        let mut prev_q: u8 = 0;
        let mut prev_q2: u8 = 0;

        for j in 0..this_len {
            let cur_base = base_to_idx(sb[j]);
            let prev_base = if j > 0 { base_to_idx(sb[j - 1]) } else { 4 };
            let ctx = context_id(j, prev_q, prev_q2, prev_base, cur_base);

            let slot = if ctx_slots[ctx] == u32::MAX {
                let id = models.len() as u32;
                ctx_slots[ctx] = id;
                models.push(AdaptiveModel::new(n_symbols));
                id as usize
            } else {
                ctx_slots[ctx] as usize
            };

            let sym = decoder.decode(&models[slot].cum_freqs, n_symbols, models[slot].total);
            models[slot].update(sym);

            let phred = symbols[sym];
            qual_bytes.push(phred + 33);
            prev_q2 = prev_q;
            prev_q = phred;
        }

        result.push(qual_bytes);
    }

    Ok(result)
}

/// Decompress multi-block quality_ctx data.
///
/// Format: `[num_blocks: u32] [block_len: u32] [quality_ctx_blob]...`
/// Each blob is a self-describing quality_ctx archive (contains num_reads in its header).
/// Sequences are consumed in order: first blob's num_reads from the front, etc.
pub fn decompress_quality_ctx_multiblock(
    data: &[u8],
    sequences: &[Vec<u8>],
) -> Result<Vec<Vec<u8>>> {
    if data.is_empty() || sequences.is_empty() {
        return Ok(Vec::new());
    }

    let num_blocks = super::read_le_u32(data, 0)? as usize;
    let mut offset = 4;

    // Phase 1: parse block boundaries and compute sequence ranges (fast, sequential)
    let mut block_info: Vec<(&[u8], usize, usize)> = Vec::with_capacity(num_blocks); // (blob, seq_start, num_reads)
    let mut seq_cursor = 0;

    for _ in 0..num_blocks {
        if offset + 4 > data.len() {
            anyhow::bail!("quality_ctx multiblock: truncated block length");
        }
        let block_len = super::read_le_u32(data, offset)? as usize;
        offset += 4;

        if offset + block_len > data.len() {
            anyhow::bail!("quality_ctx multiblock: truncated block data");
        }
        let blob = &data[offset..offset + block_len];
        offset += block_len;

        // Extract num_reads from blob header:
        // [n_symbols:1B][symbols:nB][read_len:2B][num_reads:4B]
        let n_symbols = blob[0] as usize;
        let nr_offset = 1 + n_symbols + 2;
        let chunk_reads = super::read_le_u32(blob, nr_offset)? as usize;

        block_info.push((blob, seq_cursor, chunk_reads));
        seq_cursor += chunk_reads;
    }

    // Phase 2: decode all blocks in parallel (each block is independent)
    use rayon::prelude::*;
    let decoded_blocks: Vec<Result<Vec<Vec<u8>>>> = block_info
        .into_par_iter()
        .map(|(blob, seq_start, chunk_reads)| {
            let seq_refs: Vec<&[u8]> = sequences[seq_start..seq_start + chunk_reads]
                .iter()
                .map(|s| s.as_slice())
                .collect();
            decompress_qualities_ctx(blob, &seq_refs, chunk_reads)
        })
        .collect();

    // Phase 3: concatenate in order
    let mut all_qualities = Vec::with_capacity(sequences.len());
    for result in decoded_blocks {
        all_qualities.extend(result?);
    }

    Ok(all_qualities)
}

/// Wrap a single quality_ctx blob in multi-block format (num_blocks=1).
pub fn wrap_as_multiblock(blob: Vec<u8>) -> Vec<u8> {
    let mut result = Vec::with_capacity(8 + blob.len());
    result.extend_from_slice(&1u32.to_le_bytes()); // num_blocks = 1
    result.extend_from_slice(&(blob.len() as u32).to_le_bytes()); // block_len
    result.extend(blob);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_range_coder_basic() {
        let mut enc = RangeEncoder::new();
        enc.encode(0, 1, 4); // symbol 0
        enc.encode(1, 1, 4); // symbol 1
        enc.encode(2, 1, 4); // symbol 2
        enc.encode(3, 1, 4); // symbol 3
        enc.encode(0, 1, 4); // symbol 0
        let compressed = enc.finish();

        let mut dec = RangeDecoder::new(&compressed);
        let cum = &[0u32, 1, 2, 3, 4];
        assert_eq!(dec.decode(cum, 4, 4), 0);
        assert_eq!(dec.decode(cum, 4, 4), 1);
        assert_eq!(dec.decode(cum, 4, 4), 2);
        assert_eq!(dec.decode(cum, 4, 4), 3);
        assert_eq!(dec.decode(cum, 4, 4), 0);
    }

    #[test]
    fn test_range_coder_skewed() {
        let mut enc = RangeEncoder::new();
        let total = 100u32;
        let syms: Vec<usize> = (0..1000)
            .map(|i| if i % 33 == 0 { 1 } else { 0 })
            .collect();
        for &s in &syms {
            let (cum, freq) = if s == 0 { (0, 97) } else { (97, 3) };
            enc.encode(cum, freq, total);
        }
        let compressed = enc.finish();

        let mut dec = RangeDecoder::new(&compressed);
        let cum = &[0u32, 97, 100];
        for &expected in &syms {
            let got = dec.decode(cum, 2, total);
            assert_eq!(got, expected);
        }
    }

    #[test]
    fn test_range_coder_many_symbols() {
        let n = 20;
        let mut enc = RangeEncoder::new();
        let total = n as u32;
        let cum: Vec<u32> = (0..=n).map(|i| i as u32).collect();

        let syms: Vec<usize> = (0..5000).map(|i| (i * 7 + 3) % n).collect();
        for &s in &syms {
            enc.encode(cum[s], cum[s + 1] - cum[s], total);
        }
        let compressed = enc.finish();

        let mut dec = RangeDecoder::new(&compressed);
        for &expected in &syms {
            let got = dec.decode(&cum, n, total);
            assert_eq!(got, expected);
        }
    }

    #[test]
    fn test_roundtrip_basic() {
        let sequences: Vec<&[u8]> = vec![b"ACGTACGTAC"; 10];
        let qualities: Vec<&[u8]> = vec![b"IIIIIIIIIB"; 10];

        let compressed = compress_qualities_ctx(&qualities, &sequences).unwrap();

        let decompressed =
            decompress_qualities_ctx(&compressed, &sequences, 10).unwrap();

        for i in 0..10 {
            assert_eq!(decompressed[i], qualities[i], "read {i} mismatch");
        }
    }

    #[test]
    fn test_roundtrip_varied() {
        let sequences: Vec<&[u8]> = vec![
            b"AAACCCGGGTTT",
            b"TTTGGGCCCAAA",
            b"ACGTACGTACGT",
            b"NNNNACGTNNNN",
        ];
        let qualities: Vec<&[u8]> = vec![
            b"IIIIII555555",
            b"555555IIIIII",
            b"I5I5I5I5I5I5",
            b"!!!!IIII!!!!",
        ];

        let compressed = compress_qualities_ctx(&qualities, &sequences).unwrap();

        let decompressed =
            decompress_qualities_ctx(&compressed, &sequences, 4).unwrap();

        for i in 0..4 {
            assert_eq!(decompressed[i], qualities[i], "read {i} mismatch");
        }
    }

    #[test]
    fn test_roundtrip_many_reads() {
        let n = 2000;
        let sequences: Vec<&[u8]> = vec![b"ACGTACGTACGTACGT"; n];
        let qualities: Vec<&[u8]> = vec![b"IIII5555!!!!BBBB"; n];

        let compressed = compress_qualities_ctx(&qualities, &sequences).unwrap();
        let decompressed = decompress_qualities_ctx(&compressed, &sequences, n).unwrap();

        for i in 0..n {
            assert_eq!(decompressed[i], qualities[i], "read {i} mismatch");
        }
    }

    #[test]
    fn test_roundtrip_variable_length() {
        let sequences: Vec<&[u8]> = vec![
            b"ACGTACGTAC",       // 10bp
            b"TTTGGGCCCAAACCC",  // 15bp
            b"ACGTAC",           // 6bp
            b"NNNNACGTNNNNACGT", // 16bp
        ];
        let qualities: Vec<&[u8]> = vec![
            b"IIIIII5555",       // 10
            b"555555IIIIIIIII",  // 15
            b"I5I5I5",           // 6
            b"!!!!IIII!!!!IIII", // 16
        ];

        let compressed = compress_qualities_ctx(&qualities, &sequences).unwrap();

        // Verify read_len == 0 sentinel in header
        // Format: [n_symbols:1B][symbols:nB][read_len:2B]...
        let n_sym = compressed[0] as usize;
        let stored_read_len =
            compressed[1 + n_sym] as u16 | ((compressed[1 + n_sym + 1] as u16) << 8);
        assert_eq!(stored_read_len, 0, "variable-length sentinel should be 0");

        let decompressed = decompress_qualities_ctx(&compressed, &sequences, 4).unwrap();

        for i in 0..4 {
            assert_eq!(decompressed[i], qualities[i], "read {i} mismatch");
        }
    }

    #[test]
    fn test_empty() {
        let empty: Vec<&[u8]> = vec![];
        let compressed = compress_qualities_ctx(&empty, &empty).unwrap();
        assert!(compressed.is_empty());
    }
}
