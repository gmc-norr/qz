use crate::cli::QualityReduction;
use crate::compression::md_codec;
use crate::io::bam::{RawBamRecord, pack_nibbles};
use anyhow::Result;
use std::collections::HashMap;

/// Per-position "keep full quality" mask for lossy variant-aware reduction.
/// A position is a variant site if a non-consensus, non-N allele has
/// `≥ min_alt_count` reads and `≥ min_alt_frac · depth`; the mask is then
/// dilated by `± window`. `bases` holds the consensus nibble per position.
fn compute_keep_mask(counts: &[BaseCounts], bases: &[u8], r: &QualityReduction) -> Vec<bool> {
    let n = counts.len();
    let mut variant = vec![false; n];
    for i in 0..n {
        let c = &counts[i];
        let cons = bases[i];
        let depth: u32 = c.counts.iter().map(|&x| x as u32).sum();
        if depth == 0 {
            continue;
        }
        let mut alt = 0u16;
        for nib in 0u8..16 {
            if nib == cons || nib == 15 {
                continue; // skip the consensus allele and N
            }
            alt = alt.max(c.counts[nib as usize]);
        }
        if alt >= r.min_alt_count && (alt as f32) >= r.min_alt_frac * depth as f32 {
            variant[i] = true;
        }
    }
    if r.window <= 0 {
        return variant;
    }
    let w = r.window as usize;
    let mut keep = vec![false; n];
    for i in 0..n {
        if variant[i] {
            let lo = i.saturating_sub(w);
            let hi = (i + w + 1).min(n);
            keep[lo..hi].iter_mut().for_each(|k| *k = true);
        }
    }
    keep
}

/// Flatten the quality of confident reference-matching bases in place (Phred+33
/// ascii). A base keeps full resolution if its column is a keep-site (variant
/// ± window), it mismatches the consensus, or it's an inserted base; soft-clips
/// and consensus-matching confident bases are flattened per `r.scheme`.
fn apply_quality_reduction(
    qual_ascii: &mut [u8],
    nibbles: &[u8],
    cigar_ops: &[(u8, u32)],
    ref_id: i32,
    pos: i32,
    consensus: &ChunkConsensus,
    r: &QualityReduction,
) {
    let seg = consensus.get_segment_for_pos(ref_id, pos);
    let flatten = |q_ascii: &mut u8| {
        *q_ascii = r.scheme.flatten(q_ascii.wrapping_sub(33)).wrapping_add(33);
    };
    let mut ref_pos = pos;
    let mut qpos = 0usize;
    for &(op, len) in cigar_ops {
        let len = len as usize;
        match op {
            CIGAR_M | CIGAR_EQ | CIGAR_X => {
                for _ in 0..len {
                    if qpos < qual_ascii.len() {
                        let cons = seg.map_or(0, |s| s.get(ref_pos));
                        let mismatch = nibbles.get(qpos).copied().unwrap_or(0) != cons;
                        let keep = mismatch || seg.is_some_and(|s| s.is_keep_quality(ref_pos));
                        if !keep {
                            flatten(&mut qual_ascii[qpos]);
                        }
                    }
                    ref_pos += 1;
                    qpos += 1;
                }
            }
            CIGAR_I => qpos += len, // insertions kept (variant evidence)
            CIGAR_S => {
                for _ in 0..len {
                    if qpos < qual_ascii.len() {
                        flatten(&mut qual_ascii[qpos]);
                    }
                    qpos += 1;
                }
            }
            CIGAR_D | CIGAR_N => ref_pos += len as i32,
            CIGAR_H | 6 => {} // hard clip / padding: nothing consumed
            _ => {
                qpos += if consumes_query(op) { len } else { 0 };
                ref_pos += if consumes_reference(op) {
                    len as i32
                } else {
                    0
                };
            }
        }
    }
}

// CIGAR op codes that consume query (read) bases
const CIGAR_M: u8 = 0;
const CIGAR_I: u8 = 1;
const CIGAR_D: u8 = 2;
const CIGAR_N: u8 = 3;
const CIGAR_S: u8 = 4;
const CIGAR_H: u8 = 5;
const CIGAR_EQ: u8 = 7;
const CIGAR_X: u8 = 8;

/// Returns true if the CIGAR op consumes query (read) bases.
fn consumes_query(op: u8) -> bool {
    matches!(op, CIGAR_M | CIGAR_I | CIGAR_S | CIGAR_EQ | CIGAR_X)
}

/// Returns true if the CIGAR op consumes reference positions.
fn consumes_reference(op: u8) -> bool {
    matches!(op, CIGAR_M | CIGAR_D | CIGAR_N | CIGAR_EQ | CIGAR_X)
}

/// Returns true if the CIGAR op aligns a query base to a reference position
/// (both consumed = M/=/X).
fn is_alignment_op(op: u8) -> bool {
    matches!(op, CIGAR_M | CIGAR_EQ | CIGAR_X)
}

// --- Varint encoding (unsigned LEB128) ---

pub fn write_varint(buf: &mut Vec<u8>, mut value: usize) {
    loop {
        let byte = (value & 0x7F) as u8;
        value >>= 7;
        if value == 0 {
            buf.push(byte);
            break;
        } else {
            buf.push(byte | 0x80);
        }
    }
}

pub fn read_varint(data: &[u8], offset: &mut usize) -> Result<usize> {
    let mut value: usize = 0;
    let mut shift: u32 = 0;
    loop {
        if *offset >= data.len() {
            anyhow::bail!("Truncated varint at offset {}", *offset);
        }
        let byte = data[*offset];
        *offset += 1;
        // Guard the shift *before* applying it. Both checks must precede the
        // `value |=` so an over-long varint is rejected rather than silently
        // truncating the high bits (which `(chunk << 63)` would drop).
        if shift >= usize::BITS {
            anyhow::bail!("Varint too large at offset {}", *offset);
        }
        let chunk = (byte & 0x7F) as usize;
        if (chunk << shift) >> shift != chunk {
            anyhow::bail!("Varint overflow at offset {}", *offset);
        }
        value |= chunk << shift;
        if byte & 0x80 == 0 {
            break;
        }
        shift += 7;
    }
    Ok(value)
}

// --- Zigzag encoding for signed integers ---

pub fn zigzag_encode(v: i32) -> u32 {
    ((v << 1) ^ (v >> 31)) as u32
}

pub fn zigzag_decode(v: u32) -> i32 {
    ((v >> 1) as i32) ^ -((v & 1) as i32)
}

// --- Nibble pair packing ---

/// Pack nibble values (0-15) into bytes, 2 per byte (high nibble first).
pub fn pack_nibble_pairs(nibbles: &[u8]) -> Vec<u8> {
    let out_len = nibbles.len().div_ceil(2);
    let mut packed = vec![0u8; out_len];
    let pairs = nibbles.len() / 2;
    for i in 0..pairs {
        packed[i] = (nibbles[i * 2] << 4) | (nibbles[i * 2 + 1] & 0x0F);
    }
    if !nibbles.len().is_multiple_of(2) {
        packed[pairs] = nibbles[pairs * 2] << 4;
    }
    packed
}

/// Unpack nibble pairs back to individual values.
pub fn unpack_nibble_pairs(packed: &[u8], count: usize) -> Vec<u8> {
    let mut nibbles = vec![0u8; count];
    let pairs = count / 2;
    for i in 0..pairs {
        nibbles[i * 2] = packed[i] >> 4;
        nibbles[i * 2 + 1] = packed[i] & 0x0F;
    }
    if !count.is_multiple_of(2) {
        nibbles[count - 1] = packed[pairs] >> 4;
    }
    nibbles
}

// --- BAM nibble to ASCII ---

const NIBBLE_TO_ASCII_TABLE: [u8; 16] = [
    b'N', b'A', b'C', b'N', b'G', b'N', b'N', b'N', b'T', b'N', b'N', b'N', b'N', b'N', b'N', b'N',
];

/// Convert BAM 4-bit nibble to ASCII base character.
pub fn nibble_to_ascii(n: u8) -> u8 {
    NIBBLE_TO_ASCII_TABLE[(n & 0x0F) as usize]
}

/// Scalar fallback for nibbles_to_ascii.
#[inline]
fn nibbles_to_ascii_scalar(nibbles: &[u8], ascii: &mut [u8]) {
    for i in 0..nibbles.len() {
        ascii[i] = NIBBLE_TO_ASCII_TABLE[(nibbles[i] & 0x0F) as usize];
    }
}

/// SSSE3 pshufb: 16-byte LUT lookup, 16 bytes at a time.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
unsafe fn nibbles_to_ascii_ssse3(nibbles: &[u8], ascii: &mut [u8]) {
    use std::arch::x86_64::*;
    let table = unsafe { _mm_loadu_si128(NIBBLE_TO_ASCII_TABLE.as_ptr() as *const __m128i) };
    let mask = _mm_set1_epi8(0x0F);

    let len = nibbles.len();
    let chunks = len / 16;
    for i in 0..chunks {
        unsafe {
            let input = _mm_loadu_si128(nibbles.as_ptr().add(i * 16) as *const __m128i);
            let masked = _mm_and_si128(input, mask);
            let result = _mm_shuffle_epi8(table, masked);
            _mm_storeu_si128(ascii.as_mut_ptr().add(i * 16) as *mut __m128i, result);
        }
    }
    // Scalar tail
    for i in (chunks * 16)..len {
        ascii[i] = NIBBLE_TO_ASCII_TABLE[(nibbles[i] & 0x0F) as usize];
    }
}

/// Batch-convert nibble array to ASCII. Uses SSSE3 pshufb when available.
pub fn nibbles_to_ascii(nibbles: &[u8]) -> Vec<u8> {
    let mut ascii = vec![0u8; nibbles.len()];
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("ssse3") {
            unsafe { nibbles_to_ascii_ssse3(nibbles, &mut ascii) };
            return ascii;
        }
    }
    nibbles_to_ascii_scalar(nibbles, &mut ascii);
    ascii
}

// --- BAM quality <-> ASCII Phred+33 ---

/// Convert BAM quality bytes to ASCII Phred+33 in bulk.
/// 0xFF (unavailable) maps to '!' (ASCII 33, Phred 0).
/// Returns true if any 0xFF byte was encountered.
#[inline]
pub fn bam_qual_to_ascii(src: &[u8], dst: &mut [u8]) -> bool {
    debug_assert_eq!(src.len(), dst.len());
    // Branchless: mask q to 0 when 0xFF, then add 33.
    // LLVM vectorizes this to vpcmpeqb + vpand + vpaddb.
    let mut any_ff = 0u8;
    for i in 0..src.len() {
        let q = src[i];
        let is_ff = (q == 0xFF) as u8;
        any_ff |= is_ff;
        let mask = !is_ff.wrapping_neg(); // 0xFF if q != 0xFF, 0x00 if q == 0xFF
        dst[i] = (q & mask).wrapping_add(33);
    }
    any_ff != 0
}

/// Convert ASCII Phred+33 quality back to BAM phred bytes in bulk.
///
/// Uses **wrapping** subtract so this is the EXACT inverse of `bam_qual_to_ascii`'s
/// `wrapping_add(33)` over the full 0..=254 raw range. A saturating subtract here
/// was not a true inverse: raw quality bytes in 223..=254 encode (wrapping_add) to
/// an ASCII value < 33, which saturating_sub then clamped to 0 instead of restoring
/// the original — making a structurally-valid (if out-of-spec) BAM fail to
/// round-trip (caught as a content-CRC mismatch on decode). 0xFF is never seen on
/// this path (chunks with unavailable quality use the byte-exact BSC fallback), so
/// real Phred values (< 223) are unaffected — wrapping_sub == saturating_sub there.
#[inline]
pub fn ascii_qual_to_bam(src: &[u8], dst: &mut [u8]) {
    debug_assert_eq!(src.len(), dst.len());
    for i in 0..src.len() {
        dst[i] = src[i].wrapping_sub(33);
    }
}

// --- Consensus ---

/// Per-position base counts for consensus building.
/// Index 0=padding/unused, 1=A, 2=C, 4=G, 8=T, 15=N in BAM nibble encoding.
/// We use a flat [u16; 16] so any nibble value can be counted directly.
#[derive(Clone)]
struct BaseCounts {
    counts: [u16; 16],
}

impl BaseCounts {
    fn new() -> Self {
        Self { counts: [0; 16] }
    }

    fn add(&mut self, nibble: u8) {
        self.counts[nibble as usize] = self.counts[nibble as usize].saturating_add(1);
    }

    /// Majority-vote consensus nibble (the nibble with the highest count).
    /// Ties broken by preferring lower nibble value.
    fn consensus(&self) -> u8 {
        let mut best = 0u8;
        let mut best_count = 0u16;
        for i in 1..16u8 {
            if self.counts[i as usize] > best_count {
                best_count = self.counts[i as usize];
                best = i;
            }
        }
        best
    }
}

/// A consensus segment covering a contiguous range of reference positions on one ref_id.
pub struct ConsensusSegment {
    pub ref_id: i32,
    pub start_pos: i32,
    /// One BAM nibble value per reference position.
    pub bases: Vec<u8>,
    /// Per-position "keep full quality" flag for lossy variant-aware reduction
    /// (variant column ± window). Empty when reduction is disabled.
    pub keep_quality: Vec<bool>,
}

impl ConsensusSegment {
    pub fn get(&self, ref_pos: i32) -> u8 {
        let idx = (ref_pos - self.start_pos) as usize;
        if idx < self.bases.len() {
            self.bases[idx]
        } else {
            0 // outside range -> XOR with 0 = identity
        }
    }

    /// Whether a reference position is a quality-keep site (variant column or
    /// within its window). False when reduction is disabled (`keep_quality` empty)
    /// or the position is outside this segment.
    pub fn is_keep_quality(&self, ref_pos: i32) -> bool {
        let idx = (ref_pos - self.start_pos) as usize;
        self.keep_quality.get(idx).copied().unwrap_or(false)
    }
}

/// Consensus data for a chunk, potentially spanning multiple chromosomes.
pub struct ChunkConsensus {
    /// Segments ordered by (ref_id, start_pos). Multiple segments may share a ref_id
    /// when reads on that chromosome are separated by a large gap.
    pub segments: Vec<ConsensusSegment>,
    /// Fast lookup: ref_id -> sorted list of segment indices (by start_pos).
    index: HashMap<i32, Vec<usize>>,
}

impl ChunkConsensus {
    /// Build consensus from a chunk of records.
    ///
    /// Reads on the same chromosome that are separated by more than GAP_THRESHOLD
    /// reference bases get separate consensus segments, avoiding huge sparse arrays
    /// for chromosomes with only a few scattered reads (e.g. supplementary alignments).
    pub fn build(records: &[RawBamRecord], reduction: Option<&QualityReduction>) -> Self {
        // For each ref_id: an ordered list of (start_pos, Vec<BaseCounts>) sub-segments.
        const GAP_THRESHOLD: i32 = 50_000;
        let mut accum: HashMap<i32, Vec<(i32, Vec<BaseCounts>)>> = HashMap::new();

        for rec in records {
            let ref_id = rec.ref_id();
            if ref_id < 0 {
                continue; // unmapped — no reference position
            }
            let pos = rec.pos();
            let n_cigar = rec.n_cigar_op();
            if n_cigar == 0 {
                continue; // no alignment
            }

            let nibbles = rec.unpack_seq_nibbles();
            let cigar_ops = rec.cigar_ops();

            // Bound the per-read reference span before growing the consensus accumulator.
            // The `resize_with` below extends `entry.1` to `ref_pos - start`, so a
            // crafted/corrupt CIGAR with huge D/N ops can drive a process-aborting
            // multi-GB allocation (this runs on a worker thread whose join recovers only
            // unwinding panics, not abort — the same hazard decode guards with
            // `MAX_REF_SPAN`). No real linear alignment spans this far; skip such a read's
            // consensus contribution. Its bases are still encoded losslessly downstream
            // (an absent consensus segment XORs to the raw base — see the seq encoder).
            const MAX_CONSENSUS_REF_SPAN: i64 = 64 * 1024 * 1024;
            let ref_span: i64 = cigar_ops
                .iter()
                .filter(|(op, _)| consumes_reference(*op))
                .map(|(_, len)| *len as i64)
                .sum();
            if ref_span > MAX_CONSENSUS_REF_SPAN {
                continue;
            }

            let sub_segs = accum.entry(ref_id).or_default();

            // Decide whether to extend the last sub-segment or start a new one.
            // BAM is coordinate-sorted so pos >= all previous positions on this ref_id.
            if sub_segs.is_empty() {
                sub_segs.push((pos, Vec::new()));
            } else {
                let last = sub_segs.last().unwrap();
                let last_end = last.0 + last.1.len() as i32;
                if pos > last_end + GAP_THRESHOLD {
                    sub_segs.push((pos, Vec::new()));
                }
            }

            let entry = sub_segs.last_mut().unwrap();

            // Handle rare case where this read starts before segment start.
            if pos < entry.0 {
                let prepend = (entry.0 - pos) as usize;
                let mut new_counts = vec![BaseCounts::new(); prepend];
                new_counts.append(&mut entry.1);
                entry.1 = new_counts;
                entry.0 = pos;
            }

            let start = entry.0;
            let mut ref_pos = pos;
            let mut query_pos = 0usize;

            for &(op, len) in &cigar_ops {
                let len = len as usize;
                match op {
                    CIGAR_M | CIGAR_EQ | CIGAR_X => {
                        for _ in 0..len {
                            let idx = (ref_pos - start) as usize;
                            if idx >= entry.1.len() {
                                entry.1.resize_with(idx + 1, BaseCounts::new);
                            }
                            if query_pos < nibbles.len() {
                                entry.1[idx].add(nibbles[query_pos]);
                            }
                            ref_pos += 1;
                            query_pos += 1;
                        }
                    }
                    CIGAR_I | CIGAR_S => {
                        query_pos += len; // insertion/softclip: consume query only
                    }
                    CIGAR_D | CIGAR_N => {
                        ref_pos += len as i32; // deletion/skip: consume reference only
                    }
                    CIGAR_H | 6 => {
                        // hard clip / padding: nothing consumed
                    }
                    _ => {
                        query_pos += if consumes_query(op) { len } else { 0 };
                        ref_pos += if consumes_reference(op) {
                            len as i32
                        } else {
                            0
                        };
                    }
                }
            }
        }

        // Convert accumulated counts to ConsensusSegments.
        let mut segments: Vec<ConsensusSegment> = Vec::new();
        let mut index: HashMap<i32, Vec<usize>> = HashMap::new();
        let mut sorted_refs: Vec<i32> = accum.keys().copied().collect();
        sorted_refs.sort();

        for ref_id in sorted_refs {
            let sub_segs = accum.remove(&ref_id).unwrap();
            let mut seg_indices = Vec::with_capacity(sub_segs.len());
            for (start_pos, counts) in sub_segs {
                let bases: Vec<u8> = counts.iter().map(|c| c.consensus()).collect();
                let keep_quality = match reduction {
                    Some(r) => compute_keep_mask(&counts, &bases, r),
                    None => Vec::new(),
                };
                seg_indices.push(segments.len());
                segments.push(ConsensusSegment {
                    ref_id,
                    start_pos,
                    bases,
                    keep_quality,
                });
            }
            index.insert(ref_id, seg_indices);
        }

        Self { segments, index }
    }

    /// Get consensus nibble at (ref_id, ref_pos). Returns 0 if not covered.
    pub fn get(&self, ref_id: i32, ref_pos: i32) -> u8 {
        self.get_segment_for_pos(ref_id, ref_pos)
            .map_or(0, |seg| seg.get(ref_pos))
    }

    /// Look up the consensus segment that covers (ref_id, pos).
    /// Performs one HashMap lookup + binary search among segments for that ref_id.
    /// Returns None if no segment covers pos (reads there XOR with 0 = no consensus help).
    pub fn get_segment_for_pos(&self, ref_id: i32, pos: i32) -> Option<&ConsensusSegment> {
        let indices = self.index.get(&ref_id)?;
        // Segments for a ref_id are in ascending start_pos order (BAM is sorted).
        // Find the last segment with start_pos <= pos.
        let idx = indices.partition_point(|&i| self.segments[i].start_pos <= pos);
        if idx == 0 {
            return None;
        }
        let seg = &self.segments[indices[idx - 1]];
        // Only return the segment if pos falls within its covered range.
        if pos < seg.start_pos + seg.bases.len() as i32 {
            Some(seg)
        } else {
            None
        }
    }

    /// Serialize consensus data to bytes for archival (nibbles packed 2 per byte).
    pub fn serialize(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&(self.segments.len() as u32).to_le_bytes());
        for seg in &self.segments {
            buf.extend_from_slice(&seg.ref_id.to_le_bytes());
            buf.extend_from_slice(&seg.start_pos.to_le_bytes());
            buf.extend_from_slice(&(seg.bases.len() as u32).to_le_bytes());
            buf.extend_from_slice(&pack_nibble_pairs(&seg.bases));
        }
        buf
    }

    /// Deserialize consensus data from bytes (nibbles packed 2 per byte).
    pub fn deserialize(data: &[u8]) -> Result<Self> {
        if data.len() < 4 {
            anyhow::bail!("Consensus data too short ({} bytes)", data.len());
        }
        let mut offset = 0;
        let num_segments =
            u32::from_le_bytes(data[offset..offset + 4].try_into().unwrap()) as usize;
        offset += 4;

        let mut segments = Vec::with_capacity(num_segments.min(1024));
        let mut index: HashMap<i32, Vec<usize>> = HashMap::new();

        for i in 0..num_segments {
            if offset + 12 > data.len() {
                anyhow::bail!("Truncated consensus segment {i} at offset {offset}");
            }
            let ref_id = i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
            offset += 4;
            let start_pos = i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
            offset += 4;
            let num_bases =
                u32::from_le_bytes(data[offset..offset + 4].try_into().unwrap()) as usize;
            offset += 4;
            let packed_len = num_bases.div_ceil(2);
            if offset + packed_len > data.len() {
                anyhow::bail!(
                    "Truncated consensus bases: need {packed_len} bytes at offset {offset}, have {}",
                    data.len() - offset
                );
            }
            let bases = unpack_nibble_pairs(&data[offset..offset + packed_len], num_bases);
            offset += packed_len;
            index.entry(ref_id).or_default().push(i);
            segments.push(ConsensusSegment {
                ref_id,
                start_pos,
                bases,
                keep_quality: Vec::new(), // decode-side: keep-mask is compress-only
            });
        }

        Ok(Self { segments, index })
    }
}

// --- Columnar streams ---

/// All columnar streams for a chunk of BAM records.
pub struct BamStreams {
    // Per-chunk consensus data (serialized)
    pub consensus: Vec<u8>,

    // Fixed-width per record
    pub ref_id: Vec<u8>,
    pub pos: Vec<u8>,
    pub mapq: Vec<u8>,
    pub bin: Vec<u8>,
    pub flag: Vec<u8>,
    pub next_ref_id: Vec<u8>,
    pub next_pos: Vec<u8>,
    pub tlen: Vec<u8>,

    // Variable-length per record (varint-prefixed)
    pub read_name: Vec<u8>,
    pub cigar: Vec<u8>,
    pub seq_diff: Vec<u8>,
    pub seq_extra: Vec<u8>,
    pub qual: Vec<u8>,
    pub aux: Vec<u8>,

    // MD derivation (v8): aux above stores MD-stripped blobs for reads whose MD
    // was derived; these three streams carry the metadata to regenerate it.
    /// Per-record bit (bit-packed, LSB-first): 1 = MD was stripped & is derived.
    pub md_present: Vec<u8>,
    /// Per stripped-MD record: varint tag-position of MD within its aux blob.
    pub md_index: Vec<u8>,
    /// Chunk-level: serialized (ref≠consensus) exception map.
    pub md_exceptions: Vec<u8>,

    // Mate-derivation (v9). MC:Z is regenerated from the mate's CIGAR; TLEN is
    // derived from the two mates' mapped extents. tlen (stream 8) holds only the
    // non-derivable values; the rest are reconstructed on decode.
    /// Per-record bit (bit-packed): 1 = MC:Z was stripped & is mate-derived.
    pub mc_present: Vec<u8>,
    /// Per derived-MC record: varint tag-position of MC within its (MD-stripped) aux.
    pub mc_index: Vec<u8>,
    /// Per-record bit (bit-packed): 1 = TLEN was derived (omitted from stream 8).
    pub tlen_derivable: Vec<u8>,

    // PNEXT/RNEXT derivation (v10). When the primary mate is in this chunk,
    // next_pos == mate.pos and next_ref_id == mate.ref_id exactly (SAM spec, no
    // tie convention). Such records set the bit and omit both values from streams
    // 6/7; the rest store the delta as before. Decode rebuilds them from the mate.
    /// Per-record bit (bit-packed): 1 = next_pos/next_ref_id are mate-derived.
    pub next_derivable: Vec<u8>,

    // NM:i derivation (v11). NM = MD-substitutions + CIGAR(I+D), derived from the
    // recovered reference (only when MD is also derivable). Stripped reads set the
    // bit and record the aux tag-position + the source integer type byte; the
    // value is fully regenerated on decode.
    /// Per-record bit (bit-packed): 1 = NM:i was stripped & is derived.
    pub nm_present: Vec<u8>,
    /// Per stripped-NM record: varint tag-position of NM within its (MD+MC-stripped) aux.
    pub nm_index: Vec<u8>,
    /// Per stripped-NM record: the source integer type byte (c/C/s/S/i/I) to reproduce.
    pub nm_type: Vec<u8>,

    pub num_records: usize,
}

/// Number of streams in a BamStreams.
pub const NUM_STREAMS: usize = 25;

impl BamStreams {
    /// Get all streams as an array of references (in archive order).
    pub fn as_slices(&self) -> [&[u8]; NUM_STREAMS] {
        [
            &self.consensus,
            &self.ref_id,
            &self.pos,
            &self.mapq,
            &self.bin,
            &self.flag,
            &self.next_ref_id,
            &self.next_pos,
            &self.tlen,
            &self.read_name,
            &self.cigar,
            &self.seq_diff,
            &self.seq_extra,
            &self.qual,
            &self.aux,
            &self.md_present,
            &self.md_index,
            &self.md_exceptions,
            &self.mc_present,
            &self.mc_index,
            &self.tlen_derivable,
            &self.next_derivable,
            &self.nm_present,
            &self.nm_index,
            &self.nm_type,
        ]
    }
}

/// Per-record quality and sequence data for fqz compression.
pub struct ChunkQualityData {
    /// Per-record quality bytes in ASCII Phred+33 format.
    pub qualities_ascii: Vec<Vec<u8>>,
    /// Per-record reverse-strand bit (FLAG & 0x10), used to reorient quality for
    /// the v12 fqzcomp path.
    pub strands: Vec<bool>,
    /// True if any quality byte was 0xFF (unavailable).
    pub has_unavailable_quality: bool,
}

/// Result from records_to_streams: columnar streams + quality data for fqz.
pub struct StreamsResult {
    pub streams: BamStreams,
    pub quality_data: ChunkQualityData,
}

/// Build columnar streams from a chunk of raw BAM records.
/// Two-pass: first builds consensus, then delta-encodes sequences.
/// Returns streams plus per-record quality/sequence data for fqz.
///
/// Pass 2 (per-record extraction) is split into contiguous record ranges and run
/// in parallel; each range appends to range-local streams, then the ranges are
/// concatenated in record order, producing **byte-identical output to a serial
/// scan**. This is safe because: (a) the delta-encoded fields (ref_id/pos) depend
/// only on the previous record's *absolute* value, which each range seeds from
/// `records[start-1]`; (b) range starts are multiples of 8, so the bit-packed
/// streams' bytes are disjoint and concatenate cleanly; and (c) the md_exceptions
/// map is merged in record order (later records override earlier, matching the
/// serial last-write-wins per key).
pub fn records_to_streams(
    records: &[RawBamRecord],
    reduction: Option<&QualityReduction>,
) -> StreamsResult {
    use rayon::prelude::*;

    let n = records.len();

    // Pass 1: build consensus (whole chunk). When `reduction` is set, the
    // consensus also carries a per-position quality-keep mask.
    let consensus = ChunkConsensus::build(records, reduction);
    let consensus_bytes = consensus.serialize();

    // Mate map (deterministic; rebuilt identically on decode from names+flags).
    let names: Vec<&[u8]> = records.iter().map(|r| r.read_name()).collect();
    let flags: Vec<u16> = records.iter().map(|r| r.flag()).collect();
    let mate = md_codec::build_mate_map(&names, &flags);

    // Pass 2: extract fields + delta-encode sequences, in parallel over ranges.
    let ranges = split_record_ranges(n);
    let parts: Vec<RangeStreams> = ranges
        .par_iter()
        .map(|&(start, end)| build_streams_range(records, start, end, &mate, &consensus, reduction))
        .collect();

    assemble_streams(consensus_bytes, n, parts)
}

/// Populate the raw quality stream (per record: varint l_seq + Phred bytes) for
/// the BSC fallback path, in record order. Built lazily — only chunks with
/// unavailable (0xFF) quality use the fallback, and those bytes can't be
/// recovered from `qualities_ascii` (0xFF is masked there), so this reads the
/// originals from the records. Normal chunks use the fqzcomp path and leave the
/// stream empty (slot 13 is replaced by the fqzcomp blob).
pub fn populate_qual_stream(streams: &mut BamStreams, records: &[RawBamRecord]) {
    let total: usize = records.iter().map(|r| r.qual_bytes().len() + 5).sum();
    streams.qual.clear();
    streams.qual.reserve(total);
    for rec in records {
        write_varint(&mut streams.qual, rec.l_seq() as usize);
        streams.qual.extend_from_slice(rec.qual_bytes());
    }
}

/// Split `n` records into contiguous ranges with 8-aligned starts (required so
/// each range's bit-packed bytes are disjoint and concatenate cleanly). Targets
/// ~64 parts for load balancing; small chunks collapse to a single range.
fn split_record_ranges(n: usize) -> Vec<(usize, usize)> {
    if n == 0 {
        return Vec::new();
    }
    const TARGET_PARTS: usize = 64;
    let step = (n / TARGET_PARTS).max(1).next_multiple_of(8);
    let mut ranges = Vec::new();
    let mut start = 0;
    while start < n {
        let end = (start + step).min(n);
        ranges.push((start, end));
        start = end;
    }
    ranges
}

/// Range-local output of [`build_streams_range`], concatenated by [`assemble_streams`].
struct RangeStreams {
    /// Per-record byte streams for this range (consensus empty; num_records = range len;
    /// bit-packed streams sized for the range and concatenated at assemble time).
    streams: BamStreams,
    qualities_ascii: Vec<Vec<u8>>,
    strands: Vec<bool>,
    has_unavailable_quality: bool,
    md_exceptions: std::collections::BTreeMap<(i32, i32), u8>,
}

/// Build range-local streams for `records[start..end]`. `start` must be a
/// multiple of 8. The per-record logic is identical to a serial scan; the only
/// range-awareness is the delta seed (`records[start-1]`) and local bit-packed
/// indexing (`li = rec_i - start`, valid because `start % 8 == 0`).
fn build_streams_range(
    records: &[RawBamRecord],
    start: usize,
    end: usize,
    mate: &[u32],
    consensus: &ChunkConsensus,
    reduction: Option<&QualityReduction>,
) -> RangeStreams {
    let range_len = end - start;

    let mut streams = BamStreams {
        consensus: Vec::new(),
        ref_id: Vec::with_capacity(range_len * 4),
        pos: Vec::with_capacity(range_len * 4),
        mapq: Vec::with_capacity(range_len),
        bin: Vec::with_capacity(range_len * 2),
        flag: Vec::with_capacity(range_len * 2),
        next_ref_id: Vec::with_capacity(range_len * 4),
        next_pos: Vec::with_capacity(range_len * 4),
        tlen: Vec::with_capacity(range_len * 4),
        read_name: Vec::with_capacity(range_len * 32),
        cigar: Vec::with_capacity(range_len * 24),
        seq_diff: Vec::with_capacity(range_len * 80),
        seq_extra: Vec::with_capacity(range_len * 16),
        qual: Vec::new(), // built lazily only for the 0xFF BSC-fallback path
        aux: Vec::with_capacity(range_len * 64),
        md_present: vec![0u8; range_len.div_ceil(8)],
        md_index: Vec::with_capacity(range_len),
        md_exceptions: Vec::new(),
        mc_present: vec![0u8; range_len.div_ceil(8)],
        mc_index: Vec::with_capacity(range_len),
        tlen_derivable: vec![0u8; range_len.div_ceil(8)],
        next_derivable: vec![0u8; range_len.div_ceil(8)],
        nm_present: vec![0u8; range_len.div_ceil(8)],
        nm_index: Vec::with_capacity(range_len),
        nm_type: Vec::with_capacity(range_len),
        num_records: range_len,
    };

    let mut qualities_ascii: Vec<Vec<u8>> = Vec::with_capacity(range_len);
    let mut strands: Vec<bool> = Vec::with_capacity(range_len);
    let mut has_unavailable_quality = false;
    // MD derivation: accumulate (ref_id,pos)->ref base where ref != consensus.
    let mut md_exceptions: std::collections::BTreeMap<(i32, i32), u8> =
        std::collections::BTreeMap::new();

    // Delta seed: previous record's absolute ref_id/pos (0/0 for the first record).
    let mut prev_ref_id: i32 = if start == 0 {
        0
    } else {
        records[start - 1].ref_id()
    };
    let mut prev_pos: i32 = if start == 0 {
        0
    } else {
        records[start - 1].pos()
    };

    for rec_i in start..end {
        let rec = &records[rec_i];
        let li = rec_i - start; // local index for bit-packed streams (start % 8 == 0)

        // Delta-encode ref_id
        let ref_id = rec.ref_id();
        let delta_ref = ref_id.wrapping_sub(prev_ref_id);
        streams.ref_id.extend_from_slice(&delta_ref.to_le_bytes());

        // Delta-encode pos (reset on chromosome boundary)
        let pos = rec.pos();
        let delta_pos = if ref_id == prev_ref_id {
            pos.wrapping_sub(prev_pos)
        } else {
            pos
        };
        streams.pos.extend_from_slice(&delta_pos.to_le_bytes());

        prev_ref_id = ref_id;
        prev_pos = pos;

        // Fixed-width fields
        streams.mapq.push(rec.mapq());
        streams.bin.extend_from_slice(&rec.bin().to_le_bytes());
        streams.flag.extend_from_slice(&rec.flag().to_le_bytes());

        // PNEXT/RNEXT (v10): when the primary mate is in this chunk, both
        // next_pos == mate.pos and next_ref_id == mate.ref_id hold exactly (SAM
        // spec). Set the derivable bit and omit both from streams 6/7; the rest
        // store the delta as before. (The two fields are always derivable
        // together for an in-chunk mate, so one bit governs both.)
        let next_derivable = {
            let mi = mate[rec_i];
            mi != u32::MAX && {
                let m = &records[mi as usize];
                rec.next_pos() == m.pos() && rec.next_ref_id() == m.ref_id()
            }
        };
        if next_derivable {
            streams.next_derivable[li / 8] |= 1 << (li % 8);
        } else {
            // Delta-encode next_ref_id relative to ref_id
            let delta_next_ref = rec.next_ref_id().wrapping_sub(ref_id);
            streams
                .next_ref_id
                .extend_from_slice(&delta_next_ref.to_le_bytes());

            // Delta-encode next_pos relative to pos
            let delta_next_pos = rec.next_pos().wrapping_sub(pos);
            streams
                .next_pos
                .extend_from_slice(&delta_next_pos.to_le_bytes());
        }

        // tlen handled below (after cigar_ops), conditional on mate-derivability.

        // Read name (varint + raw bytes including NUL)
        let rn = rec.read_name();
        write_varint(&mut streams.read_name, rn.len());
        streams.read_name.extend_from_slice(rn);

        // CIGAR (varint n_ops + raw 4-byte ops)
        let n_cigar = rec.n_cigar_op() as usize;
        write_varint(&mut streams.cigar, n_cigar);
        streams.cigar.extend_from_slice(rec.cigar_bytes());

        // Quality stream (raw varint+bytes) is built lazily: only the rare
        // unavailable-quality (0xFF) chunks fall back to BSC and consume it, and
        // those bytes can't be recovered from qualities_ascii (0xFF is masked),
        // so it's populated from the records in compress_one_chunk only when the
        // fallback fires. Normal chunks use the fqzcomp path and never touch it.

        // Sequence nibbles + CIGAR ops (also used by MD stripping below).
        let nibbles = rec.unpack_seq_nibbles();
        let cigar_ops = rec.cigar_ops();

        // TLEN: derive from the two mates' mapped extents when possible; store
        // only non-derivable values (the rest are regenerated on decode).
        let tlen = rec.tlen();
        let mut tlen_derivable = false;
        if tlen != 0 {
            let mi = mate[rec_i];
            if mi != u32::MAX {
                let m = &records[mi as usize];
                if m.ref_id() == ref_id {
                    let derived = md_codec::derive_tlen(
                        pos,
                        &cigar_ops,
                        m.pos(),
                        &m.cigar_ops(),
                        rec.flag() & 0x40 != 0,
                    );
                    tlen_derivable = derived == tlen;
                }
            }
        }
        if tlen_derivable {
            streams.tlen_derivable[li / 8] |= 1 << (li % 8);
        } else {
            streams
                .tlen
                .extend_from_slice(&zigzag_encode(tlen).to_le_bytes());
        }

        // Aux tags: strip MD:Z (consensus-derivable), then MC:Z (mate-CIGAR), then
        // NM:i (MD-substitutions + indels), all regenerated byte-exact on decode.
        // Strip MD, MC, NM in that order; decode inserts NM, MC, MD (reverse) so
        // each tag index stays valid. Non-derivable tags are kept verbatim.
        let aux = rec.aux_bytes();
        let consensus_seg = if n_cigar > 0 {
            consensus.get_segment_for_pos(ref_id, pos)
        } else {
            None
        };
        let mut cur_aux: std::borrow::Cow<[u8]> = std::borrow::Cow::Borrowed(aux);
        // Recovered reference bases (ref order), shared by MD and NM derivation;
        // Some(..) only when MD was derivable, which NM derivation requires.
        let mut md_ref_bases: Option<Vec<u8>> = None;
        if n_cigar > 0 {
            // ASCII read bases feed both MD and NM derivation (NM only fires when MD was
            // derivable). Cache lazily so it is computed at most once per record — never
            // twice, and never at all for a record with no derivable MD/NM tag.
            let mut ascii_cache: Option<Vec<u8>> = None;
            // MD
            if let Some(st) = md_codec::strip_tag(&cur_aux, *b"MD") {
                let ascii: &[u8] = ascii_cache.get_or_insert_with(|| nibbles_to_ascii(&nibbles));
                // Self-check: only strip MD if the decoder's `regen_md` (which is
                // what reconstructs it on decompress) reproduces the original
                // value byte-exactly. `recover_ref` accepts some non-canonical MD
                // strings — e.g. a matched base written as a mismatch (`4A3` where
                // the read base is also `A`), redundant `0` runs, lowercase —
                // that `regen_md` normalizes differently. Stripping those would
                // make the archive fail its content-CRC at decompress
                // (undecodable), so keep MD verbatim in that case instead.
                let derivable = md_codec::recover_ref(&cigar_ops, &st.value, ascii)
                    .filter(|rb| md_codec::regen_md(&cigar_ops, rb, ascii) == st.value);
                if let Some(ref_bases) = derivable {
                    for (k, &rb) in ref_bases.iter().enumerate() {
                        let p = pos + k as i32;
                        let cons = nibble_to_ascii(consensus_seg.map_or(0, |s| s.get(p)));
                        if cons != rb {
                            md_exceptions.insert((ref_id, p), rb);
                        }
                    }
                    streams.md_present[li / 8] |= 1 << (li % 8);
                    write_varint(&mut streams.md_index, st.index as usize);
                    cur_aux = std::borrow::Cow::Owned(st.aux);
                    md_ref_bases = Some(ref_bases);
                }
            }
            // MC = mate's CIGAR text
            let mi = mate[rec_i];
            if mi != u32::MAX
                && let Some(st) = md_codec::strip_tag(&cur_aux, *b"MC")
                && st.value == md_codec::cigar_text(&records[mi as usize].cigar_ops())
            {
                streams.mc_present[li / 8] |= 1 << (li % 8);
                write_varint(&mut streams.mc_index, st.index as usize);
                cur_aux = std::borrow::Cow::Owned(st.aux);
            }
            // NM:i — only when MD was derivable (so decode has the ref bases) and
            // the derived value matches the stored one byte-exact. Preserve the
            // source integer type byte (aligner-specific: C vs c vs ...).
            if let Some(ref ref_bases) = md_ref_bases
                && let Some((ty, val, idx, stripped)) = md_codec::strip_int_tag(&cur_aux, *b"NM")
            {
                let ascii: &[u8] = ascii_cache.get_or_insert_with(|| nibbles_to_ascii(&nibbles));
                if md_codec::derive_nm(&cigar_ops, ref_bases, ascii) == Some(val)
                    && md_codec::encode_int_item(*b"NM", ty, val).is_some()
                {
                    streams.nm_present[li / 8] |= 1 << (li % 8);
                    write_varint(&mut streams.nm_index, idx as usize);
                    streams.nm_type.push(ty);
                    cur_aux = std::borrow::Cow::Owned(stripped);
                }
            }
        }
        write_varint(&mut streams.aux, cur_aux.len());
        streams.aux.extend_from_slice(&cur_aux);

        // Sequence: consensus-delta encoding

        let mut diff_nibbles: Vec<u8> = Vec::new();
        let mut extra_nibbles: Vec<u8> = Vec::new();

        if n_cigar == 0 {
            // Unmapped / no CIGAR: all bases go to extra
            extra_nibbles.extend_from_slice(&nibbles);
        } else {
            // One segment lookup per record instead of one HashMap lookup per base.
            let consensus_seg = consensus.get_segment_for_pos(ref_id, pos);
            let mut ref_pos = pos;
            let mut query_pos = 0usize;

            for &(op, len) in &cigar_ops {
                let len = len as usize;
                if is_alignment_op(op) {
                    // M/=/X: XOR with consensus
                    for _ in 0..len {
                        if query_pos < nibbles.len() {
                            let cons_base = consensus_seg.map_or(0, |seg| seg.get(ref_pos));
                            diff_nibbles.push(nibbles[query_pos] ^ cons_base);
                        }
                        ref_pos += 1;
                        query_pos += 1;
                    }
                } else if op == CIGAR_I || op == CIGAR_S {
                    // Insertion / soft-clip: raw to extra
                    for _ in 0..len {
                        if query_pos < nibbles.len() {
                            extra_nibbles.push(nibbles[query_pos]);
                        }
                        query_pos += 1;
                    }
                } else if op == CIGAR_D || op == CIGAR_N {
                    ref_pos += len as i32;
                } else if consumes_query(op) {
                    query_pos += len;
                    if consumes_reference(op) {
                        ref_pos += len as i32;
                    }
                } else if consumes_reference(op) {
                    ref_pos += len as i32;
                }
                // H, P: nothing consumed
            }
        }

        // Write diff (varint count + packed nibble bytes)
        write_varint(&mut streams.seq_diff, diff_nibbles.len());
        streams
            .seq_diff
            .extend_from_slice(&pack_nibble_pairs(&diff_nibbles));

        // Write extra (varint count + packed nibble bytes)
        write_varint(&mut streams.seq_extra, extra_nibbles.len());
        streams
            .seq_extra
            .extend_from_slice(&pack_nibble_pairs(&extra_nibbles));

        // Collect per-record quality (Phred+33 ASCII) and ASCII sequences for fqz
        let qual_raw = rec.qual_bytes();
        let mut qual_ascii = vec![0u8; qual_raw.len()];
        let has_ff = bam_qual_to_ascii(qual_raw, &mut qual_ascii);
        if has_ff {
            has_unavailable_quality = true;
        }
        // Lossy variant-aware reduction: flatten confident reference-matching
        // bases. Skipped for unavailable-quality (0xFF) reads — those chunks fall
        // back to BSC and discard qualities_ascii. Mutates the Phred+33 ascii in
        // place; kept bases (variant column ± window, mismatches, insertions)
        // retain full resolution.
        if let Some(r) = reduction
            && !has_ff
        {
            apply_quality_reduction(
                &mut qual_ascii,
                &nibbles,
                &cigar_ops,
                ref_id,
                pos,
                consensus,
                r,
            );
        }
        qualities_ascii.push(qual_ascii);

        // Reverse-strand bit (FLAG & 0x10), used to reorient quality for fqzcomp.
        strands.push(rec.flag() & 0x10 != 0);
    }

    RangeStreams {
        streams,
        qualities_ascii,
        strands,
        has_unavailable_quality,
        md_exceptions,
    }
}

/// Concatenate range-local streams (in record order) into the final chunk
/// streams. Byte streams are appended end-to-end; the bit-packed streams
/// concatenate cleanly because every range start is 8-aligned; the md_exceptions
/// maps are merged in range order so later records win, matching the serial scan.
fn assemble_streams(consensus_bytes: Vec<u8>, n: usize, parts: Vec<RangeStreams>) -> StreamsResult {
    // Pre-reserve each final stream by summing the range lengths to avoid
    // repeated reallocation during concatenation.
    macro_rules! total {
        ($field:ident) => {
            parts.iter().map(|p| p.streams.$field.len()).sum::<usize>()
        };
    }
    let mut streams = BamStreams {
        consensus: consensus_bytes,
        ref_id: Vec::with_capacity(total!(ref_id)),
        pos: Vec::with_capacity(total!(pos)),
        mapq: Vec::with_capacity(total!(mapq)),
        bin: Vec::with_capacity(total!(bin)),
        flag: Vec::with_capacity(total!(flag)),
        next_ref_id: Vec::with_capacity(total!(next_ref_id)),
        next_pos: Vec::with_capacity(total!(next_pos)),
        tlen: Vec::with_capacity(total!(tlen)),
        read_name: Vec::with_capacity(total!(read_name)),
        cigar: Vec::with_capacity(total!(cigar)),
        seq_diff: Vec::with_capacity(total!(seq_diff)),
        seq_extra: Vec::with_capacity(total!(seq_extra)),
        qual: Vec::with_capacity(total!(qual)),
        aux: Vec::with_capacity(total!(aux)),
        md_present: Vec::with_capacity(n.div_ceil(8)),
        md_index: Vec::with_capacity(total!(md_index)),
        md_exceptions: Vec::new(),
        mc_present: Vec::with_capacity(n.div_ceil(8)),
        mc_index: Vec::with_capacity(total!(mc_index)),
        tlen_derivable: Vec::with_capacity(n.div_ceil(8)),
        next_derivable: Vec::with_capacity(n.div_ceil(8)),
        nm_present: Vec::with_capacity(n.div_ceil(8)),
        nm_index: Vec::with_capacity(total!(nm_index)),
        nm_type: Vec::with_capacity(total!(nm_type)),
        num_records: n,
    };

    let mut qualities_ascii: Vec<Vec<u8>> = Vec::with_capacity(n);
    let mut strands: Vec<bool> = Vec::with_capacity(n);
    let mut has_unavailable_quality = false;
    let mut md_exceptions: std::collections::BTreeMap<(i32, i32), u8> =
        std::collections::BTreeMap::new();

    for part in parts {
        let RangeStreams {
            streams: ps,
            qualities_ascii: pq,
            strands: pst,
            has_unavailable_quality: phu,
            md_exceptions: pex,
        } = part;
        streams.ref_id.extend_from_slice(&ps.ref_id);
        streams.pos.extend_from_slice(&ps.pos);
        streams.mapq.extend_from_slice(&ps.mapq);
        streams.bin.extend_from_slice(&ps.bin);
        streams.flag.extend_from_slice(&ps.flag);
        streams.next_ref_id.extend_from_slice(&ps.next_ref_id);
        streams.next_pos.extend_from_slice(&ps.next_pos);
        streams.tlen.extend_from_slice(&ps.tlen);
        streams.read_name.extend_from_slice(&ps.read_name);
        streams.cigar.extend_from_slice(&ps.cigar);
        streams.seq_diff.extend_from_slice(&ps.seq_diff);
        streams.seq_extra.extend_from_slice(&ps.seq_extra);
        streams.qual.extend_from_slice(&ps.qual);
        streams.aux.extend_from_slice(&ps.aux);
        streams.md_present.extend_from_slice(&ps.md_present);
        streams.md_index.extend_from_slice(&ps.md_index);
        streams.mc_present.extend_from_slice(&ps.mc_present);
        streams.mc_index.extend_from_slice(&ps.mc_index);
        streams.tlen_derivable.extend_from_slice(&ps.tlen_derivable);
        streams.next_derivable.extend_from_slice(&ps.next_derivable);
        streams.nm_present.extend_from_slice(&ps.nm_present);
        streams.nm_index.extend_from_slice(&ps.nm_index);
        streams.nm_type.extend_from_slice(&ps.nm_type);
        qualities_ascii.extend(pq);
        strands.extend(pst);
        has_unavailable_quality |= phu;
        // Later ranges override earlier on key collision (matches serial scan).
        for (k, v) in pex {
            md_exceptions.insert(k, v);
        }
    }

    // Serialize the chunk's MD exception map (ref != consensus positions).
    streams.md_exceptions = md_codec::serialize_exceptions(&md_exceptions);

    StreamsResult {
        streams,
        quality_data: ChunkQualityData {
            qualities_ascii,
            strands,
            has_unavailable_quality,
        },
    }
}

/// Reconstruct sequence nibbles from consensus + CIGAR + diff/extra nibbles.
/// Used during decompression to get sequences before fqz decompression.
pub fn reconstruct_sequence_nibbles(
    consensus: &ChunkConsensus,
    ref_id: i32,
    pos: i32,
    cigar_bytes: &[u8],
    n_cigar_op: usize,
    diff_nibbles: &[u8],
    extra_nibbles: &[u8],
    l_seq: usize,
) -> Result<Vec<u8>> {
    // SEQ='*' (l_seq=0) carries no bases even when a CIGAR is present; the
    // encoder emits no diff/extra nibbles for it, so there is nothing to
    // reconstruct. Short-circuit before walking the CIGAR to avoid over-reading.
    if l_seq == 0 {
        return Ok(Vec::new());
    }

    let mut seq_nibbles: Vec<u8> = Vec::with_capacity(l_seq);

    if n_cigar_op == 0 {
        if extra_nibbles.len() < l_seq {
            anyhow::bail!(
                "reconstruct_sequence_nibbles: extra_nibbles too short ({} < l_seq {})",
                extra_nibbles.len(),
                l_seq
            );
        }
        seq_nibbles.extend_from_slice(&extra_nibbles[..l_seq]);
    } else {
        let cigar_ops = parse_cigar_ops(cigar_bytes, n_cigar_op);
        let mut diff_pos = 0usize;
        let mut extra_pos = 0usize;
        let mut ref_pos = pos;

        // One segment lookup per record instead of per base.
        let segment = consensus.get_segment_for_pos(ref_id, pos);

        for (op, len) in cigar_ops {
            let len = len as usize;
            if is_alignment_op(op) {
                for _ in 0..len {
                    let cons_base = segment.map_or(0, |seg| seg.get(ref_pos));
                    let diff = *diff_nibbles.get(diff_pos).ok_or_else(|| {
                        anyhow::anyhow!(
                            "reconstruct_sequence_nibbles: diff stream exhausted at {} (CIGAR demands more alignment bases)",
                            diff_pos
                        )
                    })?;
                    seq_nibbles.push(diff ^ cons_base);
                    diff_pos += 1;
                    ref_pos += 1;
                }
            } else if op == CIGAR_I || op == CIGAR_S {
                for _ in 0..len {
                    let extra = *extra_nibbles.get(extra_pos).ok_or_else(|| {
                        anyhow::anyhow!(
                            "reconstruct_sequence_nibbles: extra stream exhausted at {} (CIGAR demands more inserted/clipped bases)",
                            extra_pos
                        )
                    })?;
                    seq_nibbles.push(extra);
                    extra_pos += 1;
                }
            } else if op == CIGAR_D || op == CIGAR_N {
                ref_pos += len as i32;
            }
        }
    }

    // The reconstructed base count must match the record's l_seq, or the packed
    // SEQ would disagree with the BAM header (silently corrupt record).
    if seq_nibbles.len() != l_seq {
        anyhow::bail!(
            "reconstruct_sequence_nibbles: reconstructed {} bases but l_seq is {}",
            seq_nibbles.len(),
            l_seq
        );
    }

    Ok(seq_nibbles)
}

/// Reconstruct a single BAM record's raw bytes from decompressed stream data.
/// Returns the record data (to be passed to BamWriter::write_record).
pub fn reconstruct_record(
    consensus: &ChunkConsensus,
    ref_id: i32,
    pos: i32,
    mapq: u8,
    bin_val: u16,
    flag: u16,
    next_ref_id: i32,
    next_pos: i32,
    tlen: i32,
    read_name: &[u8],
    cigar_bytes: &[u8],
    n_cigar_op: u16,
    diff_nibbles: &[u8],
    extra_nibbles: &[u8],
    qual_bytes: &[u8],
    aux_bytes: &[u8],
    l_seq: usize,
) -> Result<Vec<u8>> {
    // Reconstruct sequence from diff + extra + consensus + CIGAR
    let seq_nibbles = reconstruct_sequence_nibbles(
        consensus,
        ref_id,
        pos,
        cigar_bytes,
        n_cigar_op as usize,
        diff_nibbles,
        extra_nibbles,
        l_seq,
    )?;

    let seq_packed = pack_nibbles(&seq_nibbles);

    // Assemble record bytes
    let record_len = 32
        + read_name.len()
        + cigar_bytes.len()
        + seq_packed.len()
        + qual_bytes.len()
        + aux_bytes.len();
    let mut data = Vec::with_capacity(record_len);

    data.extend_from_slice(&ref_id.to_le_bytes()); // 0..4
    data.extend_from_slice(&pos.to_le_bytes()); // 4..8
    data.push(read_name.len() as u8); // 8
    data.push(mapq); // 9
    data.extend_from_slice(&bin_val.to_le_bytes()); // 10..12
    data.extend_from_slice(&n_cigar_op.to_le_bytes()); // 12..14
    data.extend_from_slice(&flag.to_le_bytes()); // 14..16
    data.extend_from_slice(&(l_seq as i32).to_le_bytes()); // 16..20
    data.extend_from_slice(&next_ref_id.to_le_bytes()); // 20..24
    data.extend_from_slice(&next_pos.to_le_bytes()); // 24..28
    data.extend_from_slice(&tlen.to_le_bytes()); // 28..32
    data.extend_from_slice(read_name);
    data.extend_from_slice(cigar_bytes);
    data.extend_from_slice(&seq_packed);
    data.extend_from_slice(qual_bytes);
    data.extend_from_slice(aux_bytes);

    Ok(data)
}

/// Parse CIGAR bytes into (op, len) pairs.
/// Callers must ensure `bytes.len() >= n_ops * 4`.
pub fn parse_cigar_ops(bytes: &[u8], n_ops: usize) -> Vec<(u8, u32)> {
    assert!(
        bytes.len() >= n_ops * 4,
        "parse_cigar_ops: bytes.len()={} < n_ops*4={}",
        bytes.len(),
        n_ops * 4
    );
    let mut ops = Vec::with_capacity(n_ops);
    for i in 0..n_ops {
        let start = i * 4;
        let val = u32::from_le_bytes(bytes[start..start + 4].try_into().unwrap());
        ops.push(((val & 0xF) as u8, val >> 4));
    }
    ops
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_varint_roundtrip_full_range() {
        for v in [0usize, 1, 127, 128, 16384, 1 << 32, usize::MAX] {
            let mut buf = Vec::new();
            write_varint(&mut buf, v);
            let mut off = 0;
            assert_eq!(read_varint(&buf, &mut off).unwrap(), v, "roundtrip {v}");
            assert_eq!(off, buf.len());
        }
    }

    #[test]
    fn test_read_varint_rejects_overflow() {
        // 9 continuation bytes (shifts 0..56) then a final byte at shift 63 whose
        // value bits don't fit in usize — must error, not silently truncate.
        let mut data = vec![0x80u8; 9];
        data.push(0x7F);
        let mut off = 0;
        assert!(read_varint(&data, &mut off).is_err());

        // An unterminated varint (all continuation bits) must also error rather
        // than shifting past the integer width.
        let data = vec![0x80u8; 16];
        let mut off = 0;
        assert!(read_varint(&data, &mut off).is_err());
    }

    #[test]
    fn test_reconstruct_sequence_nibbles_rejects_short_diff_without_panic() {
        // A corrupt/inconsistent archive can have a CIGAR demanding more
        // alignment bases than the diff stream holds. Reconstruction must
        // return an error rather than index out of bounds and panic (these run
        // inside par_iter; under panic="abort" a panic kills the process).
        let consensus = ChunkConsensus::build(&[], None);
        let cigar = (4u32 << 4).to_le_bytes(); // one op: 4M
        let result = reconstruct_sequence_nibbles(&consensus, 0, 0, &cigar, 1, &[], &[], 4);
        assert!(result.is_err());
    }

    #[test]
    fn test_reconstruct_sequence_nibbles_seq_star_with_cigar_is_empty() {
        // SEQ='*' is encoded as l_seq=0 even when a CIGAR is present. The
        // encoder emits no diff/extra nibbles for such a record, so decode must
        // reconstruct an empty sequence rather than over-reading the (empty)
        // diff stream and panicking.
        let consensus = ChunkConsensus::build(&[], None);
        let cigar = (4u32 << 4).to_le_bytes(); // 4M, but l_seq=0
        let result =
            reconstruct_sequence_nibbles(&consensus, 0, 0, &cigar, 1, &[], &[], 0).unwrap();
        assert!(result.is_empty());
    }
}
