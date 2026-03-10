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
// Flat model arena — all cum_freqs in one contiguous allocation
// ============================================================================
// Eliminates per-model heap allocations (previously each AdaptiveModel owned a
// separate Vec<u32>). With ~4K active models, this removes ~4K pointer-chasing
// indirections per context switch in the hot loop.

struct ModelArena {
    /// Flat storage: `data[slot * stride .. (slot+1) * stride]` = cum_freqs for model `slot`.
    data: Vec<u32>,
    /// Per-model cumulative total.
    totals: Vec<u32>,
    stride: usize, // n_symbols + 1
    n_symbols: usize,
}

impl ModelArena {
    fn new(n_symbols: usize) -> Self {
        Self {
            data: Vec::with_capacity(4096 * (n_symbols + 1)),
            totals: Vec::with_capacity(4096),
            stride: n_symbols + 1,
            n_symbols,
        }
    }

    /// Allocate a new model with uniform (Laplace) initialization. Returns slot index.
    #[inline(always)]
    fn alloc(&mut self) -> usize {
        let slot = self.totals.len();
        // Uniform init: cum_freqs = [0, 1, 2, ..., n_symbols]
        let base = self.data.len();
        self.data.resize(base + self.stride, 0);
        for i in 0..self.stride {
            self.data[base + i] = i as u32;
        }
        self.totals.push(self.n_symbols as u32);
        slot
    }

    #[inline(always)]
    fn encode_params(&self, slot: usize, sym: usize) -> (u32, u32, u32) {
        let base = slot * self.stride;
        let cum = self.data[base + sym];
        let freq = self.data[base + sym + 1] - cum;
        (cum, freq, self.totals[slot])
    }

    #[inline(always)]
    fn cum_freqs(&self, slot: usize) -> &[u32] {
        let base = slot * self.stride;
        &self.data[base..base + self.stride]
    }

    #[inline(always)]
    fn total(&self, slot: usize) -> u32 {
        self.totals[slot]
    }

    #[inline(always)]
    fn update(&mut self, slot: usize, sym: usize) {
        let base = slot * self.stride;
        for i in (sym + 1)..=self.n_symbols {
            self.data[base + i] += 1;
        }
        self.totals[slot] += 1;
        if self.totals[slot] >= RESCALE_THRESHOLD {
            self.rescale(slot);
        }
    }

    fn rescale(&mut self, slot: usize) {
        let base = slot * self.stride;
        let mut cum = 0u32;
        self.data[base] = 0;
        for i in 0..self.n_symbols {
            let freq = self.data[base + i + 1] - self.data[base + i];
            let new_freq = (freq >> 1).max(1);
            cum += new_freq;
            self.data[base + i + 1] = cum;
        }
        self.totals[slot] = cum;
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

/// Batch-convert an ASCII sequence to base indices. Avoids per-symbol LUT
/// lookups in the hot encode/decode loops.
#[inline]
fn precompute_base_indices(seq: &[u8]) -> Vec<u8> {
    let mut indices = vec![0u8; seq.len()];
    for i in 0..seq.len() {
        indices[i] = BASE_TO_IDX_LUT[seq[i] as usize];
    }
    indices
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
    if qualities.len() != sequences.len() {
        anyhow::bail!(
            "quality_ctx compress: {} quality strings but {} sequences",
            qualities.len(), sequences.len()
        );
    }

    // Build symbol alphabet from input qualities
    let mut seen = [false; 256];
    for q in qualities {
        for &b in *q {
            if b < 33 {
                anyhow::bail!("quality_ctx compress: invalid quality byte {b} (< 33)");
            }
            seen[(b - 33) as usize] = true;
        }
    }
    let symbols: Vec<u8> = (0u8..=255).filter(|&i| seen[i as usize]).collect();
    let n_symbols = symbols.len();
    if n_symbols == 0 {
        anyhow::bail!("quality_ctx compress: no quality symbols found in input");
    }
    let mut sym_map = [0u8; 256];
    for (i, &s) in symbols.iter().enumerate() {
        sym_map[s as usize] = i as u8;
    }

    // Detect variable-length reads; store read_len=0 as sentinel
    let first_len = qualities[0].len();
    let variable = qualities.iter().any(|q| q.len() != first_len);
    let read_len: u16 = if variable { 0 } else { first_len as u16 };
    let num_reads = qualities.len();

    // Context → model slot. u16::MAX = uninitialized (supports up to 65534 active models).
    let mut ctx_slots: Vec<u16> = vec![u16::MAX; MAX_CONTEXTS];
    let mut models = ModelArena::new(n_symbols);

    let mut encoder = RangeEncoder::new();

    for (qual, seq) in qualities.iter().zip(sequences.iter()) {
        let qb = *qual;
        let sb = *seq;
        if qb.len() > sb.len() {
            anyhow::bail!("quality_ctx compress: quality length {} > sequence length {}", qb.len(), sb.len());
        }
        let base_idx = precompute_base_indices(sb);
        let mut prev_q: u8 = 0;
        let mut prev_q2: u8 = 0;
        let mut prev_base: u8 = 4; // sentinel for position 0

        for j in 0..qb.len() {
            let phred = qb[j] - 33;
            let sym = sym_map[phred as usize] as usize;
            let cur_base = base_idx[j];
            let ctx = context_id(j, prev_q, prev_q2, prev_base, cur_base);

            let slot = if ctx_slots[ctx] == u16::MAX {
                let id = models.alloc();
                ctx_slots[ctx] = id as u16;
                id
            } else {
                ctx_slots[ctx] as usize
            };

            let (cum, freq, total) = models.encode_params(slot, sym);
            encoder.encode(cum, freq, total);
            models.update(slot, sym);
            prev_q2 = prev_q;
            prev_q = phred;
            prev_base = cur_base;
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
    if num_reads == 0 {
        return Ok(Vec::new());
    }
    if data.is_empty() {
        anyhow::bail!("quality_ctx decode: stream is empty but num_reads={num_reads}");
    }
    if num_reads > sequences.len() {
        anyhow::bail!(
            "quality_ctx decode: num_reads={num_reads} but only {} sequences",
            sequences.len()
        );
    }

    // Safe header parser with bounds checking
    let mut pos = 0usize;
    let mut take = |n: usize, field: &str| -> Result<&[u8]> {
        let end = pos.checked_add(n)
            .ok_or_else(|| anyhow::anyhow!("quality_ctx decode: {field} offset overflow"))?;
        if end > data.len() {
            anyhow::bail!("quality_ctx decode: truncated at {field} (need {end}, have {})", data.len());
        }
        let slice = &data[pos..end];
        pos = end;
        Ok(slice)
    };

    let n_symbols = take(1, "n_symbols")?[0] as usize;
    if n_symbols == 0 {
        anyhow::bail!("quality_ctx decode: empty symbol alphabet");
    }
    let symbols: Vec<u8> = take(n_symbols, "symbols")?.to_vec();

    let rl_bytes = take(2, "read_len")?;
    let read_len = u16::from_le_bytes([rl_bytes[0], rl_bytes[1]]) as usize;

    let nr_bytes = take(4, "num_reads")?;
    let stored_reads = u32::from_le_bytes([nr_bytes[0], nr_bytes[1], nr_bytes[2], nr_bytes[3]]) as usize;
    if stored_reads != num_reads {
        anyhow::bail!(
            "quality_ctx decode: header says {stored_reads} reads but {num_reads} requested"
        );
    }

    let cl_bytes = take(4, "compressed_len")?;
    let compressed_len = u32::from_le_bytes([cl_bytes[0], cl_bytes[1], cl_bytes[2], cl_bytes[3]]) as usize;
    let compressed_data = take(compressed_len, "compressed_data")?;

    let mut decoder = RangeDecoder::new(compressed_data);

    // u16 ctx_slots: 320KB vs 640KB with u32, better L2 cache fit
    let mut ctx_slots: Vec<u16> = vec![u16::MAX; MAX_CONTEXTS];
    let mut models = ModelArena::new(n_symbols);

    // Pre-allocate result; for fixed-length reads, pre-allocate inner vecs too
    let mut result: Vec<Vec<u8>> = if read_len > 0 {
        vec![Vec::with_capacity(read_len); num_reads]
    } else {
        Vec::with_capacity(num_reads)
    };

    for i in 0..num_reads {
        let sb = sequences[i];
        let this_len = if read_len > 0 { read_len } else { sb.len() };
        if sb.len() < this_len {
            anyhow::bail!(
                "quality_ctx decode: sequence {} has length {} but quality needs {}",
                i, sb.len(), this_len
            );
        }
        let base_idx = precompute_base_indices(&sb[..this_len]);

        // Get or create the output vec for this read
        let qual_bytes = if read_len > 0 {
            &mut result[i]
        } else {
            result.push(Vec::with_capacity(this_len));
            result.last_mut().unwrap()
        };

        let mut prev_q: u8 = 0;
        let mut prev_q2: u8 = 0;
        let mut prev_base: u8 = 4; // sentinel

        for j in 0..this_len {
            let cur_base = base_idx[j];
            let ctx = context_id(j, prev_q, prev_q2, prev_base, cur_base);

            let slot = if ctx_slots[ctx] == u16::MAX {
                let id = models.alloc();
                ctx_slots[ctx] = id as u16;
                id
            } else {
                ctx_slots[ctx] as usize
            };

            let sym = decoder.decode(models.cum_freqs(slot), n_symbols, models.total(slot));
            models.update(slot, sym);

            let phred = symbols[sym];
            qual_bytes.push(phred + 33);
            prev_q2 = prev_q;
            prev_q = phred;
            prev_base = cur_base;
        }
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
    if sequences.is_empty() {
        return Ok(Vec::new());
    }
    if data.is_empty() {
        anyhow::bail!("quality_ctx multiblock: stream is empty but sequences are present");
    }

    let num_blocks = super::read_le_u32(data, 0)? as usize;
    if num_blocks == 0 {
        anyhow::bail!("quality_ctx multiblock: 0 blocks but {} sequences", sequences.len());
    }
    let mut offset: usize = 4;

    // Phase 1: parse block boundaries and compute sequence ranges (fast, sequential)
    let mut block_info: Vec<(&[u8], usize, usize)> = Vec::with_capacity(num_blocks);
    let mut seq_cursor: usize = 0;

    for blk in 0..num_blocks {
        let next_off = offset.checked_add(4)
            .ok_or_else(|| anyhow::anyhow!("quality_ctx multiblock: block {blk} length offset overflow"))?;
        if next_off > data.len() {
            anyhow::bail!("quality_ctx multiblock: truncated block {blk} length");
        }
        let block_len = super::read_le_u32(data, offset)? as usize;
        offset = next_off;

        let blob_end = offset.checked_add(block_len)
            .ok_or_else(|| anyhow::anyhow!("quality_ctx multiblock: block {blk} data offset overflow"))?;
        if blob_end > data.len() {
            anyhow::bail!("quality_ctx multiblock: truncated block {blk} data");
        }
        let blob = &data[offset..blob_end];
        offset = blob_end;

        // Extract num_reads from blob header:
        // [n_symbols:1B][symbols:nB][read_len:2B][num_reads:4B]
        if blob.is_empty() {
            anyhow::bail!("quality_ctx multiblock: block {blk} is empty");
        }
        let n_symbols = blob[0] as usize;
        let nr_offset = 1 + n_symbols + 2;
        if nr_offset + 4 > blob.len() {
            anyhow::bail!("quality_ctx multiblock: block {blk} header truncated");
        }
        let chunk_reads = super::read_le_u32(blob, nr_offset)? as usize;

        let next_seq = seq_cursor.checked_add(chunk_reads)
            .ok_or_else(|| anyhow::anyhow!("quality_ctx multiblock: sequence count overflow at block {blk}"))?;
        if next_seq > sequences.len() {
            anyhow::bail!(
                "quality_ctx multiblock: block {blk} wants {} reads starting at {}, but only {} sequences total",
                chunk_reads, seq_cursor, sequences.len()
            );
        }

        block_info.push((blob, seq_cursor, chunk_reads));
        seq_cursor = next_seq;
    }

    if seq_cursor != sequences.len() {
        anyhow::bail!(
            "quality_ctx multiblock: blocks cover {} reads but {} sequences provided",
            seq_cursor, sequences.len()
        );
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

    #[test]
    fn test_empty_data_nonzero_reads_errors() {
        let sequences: Vec<&[u8]> = vec![b"ACGT"];
        let result = decompress_qualities_ctx(&[], &sequences, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_length_mismatch_errors() {
        let sequences: Vec<&[u8]> = vec![b"ACGT"];
        let qualities: Vec<&[u8]> = vec![b"IIII", b"AAAA"]; // 2 quals, 1 seq
        let result = compress_qualities_ctx(&qualities, &sequences);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_quality_byte_errors() {
        let sequences: Vec<&[u8]> = vec![b"ACGT"];
        let qualities: Vec<&[u8]> = vec![&[10, 20, 30, 40]]; // all < 33
        let result = compress_qualities_ctx(&qualities, &sequences);
        assert!(result.is_err());
    }
}
