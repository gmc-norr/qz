use crate::io::bam::{RawBamRecord, pack_nibbles};
use std::collections::HashMap;

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

pub fn read_varint(data: &[u8], offset: &mut usize) -> usize {
    let mut value: usize = 0;
    let mut shift = 0;
    loop {
        let byte = data[*offset];
        *offset += 1;
        value |= ((byte & 0x7F) as usize) << shift;
        if byte & 0x80 == 0 {
            break;
        }
        shift += 7;
    }
    value
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
    let mut packed = Vec::with_capacity((nibbles.len() + 1) / 2);
    for chunk in nibbles.chunks(2) {
        let byte = if chunk.len() == 2 {
            (chunk[0] << 4) | (chunk[1] & 0x0F)
        } else {
            chunk[0] << 4
        };
        packed.push(byte);
    }
    packed
}

/// Unpack nibble pairs back to individual values.
pub fn unpack_nibble_pairs(packed: &[u8], count: usize) -> Vec<u8> {
    let mut nibbles = Vec::with_capacity(count);
    for &byte in packed {
        nibbles.push(byte >> 4);
        if nibbles.len() < count {
            nibbles.push(byte & 0x0F);
        }
    }
    nibbles
}

// --- BAM nibble to ASCII ---

const NIBBLE_TO_ASCII_TABLE: [u8; 16] = [
    b'N', b'A', b'C', b'N', b'G', b'N', b'N', b'N',
    b'T', b'N', b'N', b'N', b'N', b'N', b'N', b'N',
];

/// Convert BAM 4-bit nibble to ASCII base character.
pub fn nibble_to_ascii(n: u8) -> u8 {
    NIBBLE_TO_ASCII_TABLE[(n & 0x0F) as usize]
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
}

impl ConsensusSegment {
    fn get(&self, ref_pos: i32) -> u8 {
        let idx = (ref_pos - self.start_pos) as usize;
        if idx < self.bases.len() {
            self.bases[idx]
        } else {
            0 // outside range -> XOR with 0 = identity
        }
    }
}

/// Consensus data for a chunk, potentially spanning multiple chromosomes.
pub struct ChunkConsensus {
    /// Segments keyed by ref_id, in order.
    pub segments: Vec<ConsensusSegment>,
    /// Fast lookup: ref_id -> index in segments.
    index: HashMap<i32, usize>,
}

impl ChunkConsensus {
    /// Build consensus from a chunk of records.
    pub fn build(records: &[RawBamRecord]) -> Self {
        // Accumulate per-position base counts. Key: ref_id, Value: (start_pos, Vec<BaseCounts>).
        let mut accum: HashMap<i32, (i32, Vec<BaseCounts>)> = HashMap::new();

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

            let entry = accum.entry(ref_id).or_insert_with(|| (pos, Vec::new()));

            // Ensure the range covers this read's start
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
                            // Extend if needed
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
                        ref_pos += if consumes_reference(op) { len as i32 } else { 0 };
                    }
                }
            }
        }

        // Convert counts to consensus bases
        let mut segments: Vec<ConsensusSegment> = Vec::new();
        let mut sorted_refs: Vec<i32> = accum.keys().copied().collect();
        sorted_refs.sort();

        let mut index = HashMap::new();
        for ref_id in sorted_refs {
            let (start_pos, counts) = accum.remove(&ref_id).unwrap();
            let bases: Vec<u8> = counts.iter().map(|c| c.consensus()).collect();
            index.insert(ref_id, segments.len());
            segments.push(ConsensusSegment {
                ref_id,
                start_pos,
                bases,
            });
        }

        Self { segments, index }
    }

    /// Get consensus nibble at (ref_id, ref_pos). Returns 0 if not covered.
    pub fn get(&self, ref_id: i32, ref_pos: i32) -> u8 {
        if let Some(&idx) = self.index.get(&ref_id) {
            self.segments[idx].get(ref_pos)
        } else {
            0
        }
    }

    /// Look up the consensus segment for a ref_id (one HashMap lookup per record).
    pub fn get_segment(&self, ref_id: i32) -> Option<&ConsensusSegment> {
        self.index.get(&ref_id).map(|&idx| &self.segments[idx])
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
    pub fn deserialize(data: &[u8]) -> Self {
        let mut offset = 0;
        let num_segments = u32::from_le_bytes(data[offset..offset + 4].try_into().unwrap()) as usize;
        offset += 4;

        let mut segments = Vec::with_capacity(num_segments);
        let mut index = HashMap::new();

        for i in 0..num_segments {
            let ref_id = i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
            offset += 4;
            let start_pos = i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
            offset += 4;
            let num_bases = u32::from_le_bytes(data[offset..offset + 4].try_into().unwrap()) as usize;
            offset += 4;
            let packed_len = (num_bases + 1) / 2;
            let bases = unpack_nibble_pairs(&data[offset..offset + packed_len], num_bases);
            offset += packed_len;
            index.insert(ref_id, i);
            segments.push(ConsensusSegment {
                ref_id,
                start_pos,
                bases,
            });
        }

        Self { segments, index }
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

    pub num_records: usize,
}

/// Number of streams in a BamStreams.
pub const NUM_STREAMS: usize = 15;

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
        ]
    }
}

/// Per-record quality and sequence data for quality_ctx compression.
pub struct ChunkQualityData {
    /// Per-record quality bytes in ASCII Phred+33 format.
    pub qualities_ascii: Vec<Vec<u8>>,
    /// Per-record sequence in ASCII (ACGTN) format.
    pub sequences_ascii: Vec<Vec<u8>>,
    /// True if any quality byte was 0xFF (unavailable).
    pub has_unavailable_quality: bool,
}

/// Result from records_to_streams: columnar streams + quality data for quality_ctx.
pub struct StreamsResult {
    pub streams: BamStreams,
    pub quality_data: ChunkQualityData,
}

/// Build columnar streams from a chunk of raw BAM records.
/// Two-pass: first builds consensus, then delta-encodes sequences.
/// Returns streams plus per-record quality/sequence data for quality_ctx.
pub fn records_to_streams(records: &[RawBamRecord]) -> StreamsResult {
    let n = records.len();

    // Pass 1: build consensus
    let consensus = ChunkConsensus::build(records);
    let consensus_bytes = consensus.serialize();

    // Pass 2: extract fields and delta-encode sequences
    let mut streams = BamStreams {
        consensus: consensus_bytes,
        ref_id: Vec::with_capacity(n * 4),
        pos: Vec::with_capacity(n * 4),
        mapq: Vec::with_capacity(n),
        bin: Vec::with_capacity(n * 2),
        flag: Vec::with_capacity(n * 2),
        next_ref_id: Vec::with_capacity(n * 4),
        next_pos: Vec::with_capacity(n * 4),
        tlen: Vec::with_capacity(n * 4),
        read_name: Vec::with_capacity(n * 32),
        cigar: Vec::with_capacity(n * 24),
        seq_diff: Vec::with_capacity(n * 80),
        seq_extra: Vec::with_capacity(n * 16),
        qual: Vec::with_capacity(n * 150),
        aux: Vec::with_capacity(n * 64),
        num_records: n,
    };

    // Quality data for quality_ctx
    let mut qualities_ascii: Vec<Vec<u8>> = Vec::with_capacity(n);
    let mut sequences_ascii: Vec<Vec<u8>> = Vec::with_capacity(n);
    let mut has_unavailable_quality = false;

    let mut prev_ref_id: i32 = 0;
    let mut prev_pos: i32 = 0;

    for rec in records {
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

        // Delta-encode next_ref_id relative to ref_id
        let delta_next_ref = rec.next_ref_id().wrapping_sub(ref_id);
        streams.next_ref_id.extend_from_slice(&delta_next_ref.to_le_bytes());

        // Delta-encode next_pos relative to pos
        let delta_next_pos = rec.next_pos().wrapping_sub(pos);
        streams.next_pos.extend_from_slice(&delta_next_pos.to_le_bytes());

        // Zigzag-encode tlen
        streams.tlen.extend_from_slice(&zigzag_encode(rec.tlen()).to_le_bytes());

        // Read name (varint + raw bytes including NUL)
        let rn = rec.read_name();
        write_varint(&mut streams.read_name, rn.len());
        streams.read_name.extend_from_slice(rn);

        // CIGAR (varint n_ops + raw 4-byte ops)
        let n_cigar = rec.n_cigar_op() as usize;
        write_varint(&mut streams.cigar, n_cigar);
        streams.cigar.extend_from_slice(rec.cigar_bytes());

        // Quality (varint l_seq + raw bytes) — kept in streams for BSC fallback
        let l_seq = rec.l_seq() as usize;
        write_varint(&mut streams.qual, l_seq);
        streams.qual.extend_from_slice(rec.qual_bytes());

        // Aux tags (varint len + raw bytes)
        let aux = rec.aux_bytes();
        write_varint(&mut streams.aux, aux.len());
        streams.aux.extend_from_slice(aux);

        // Sequence: consensus-delta encoding
        let nibbles = rec.unpack_seq_nibbles();
        let cigar_ops = rec.cigar_ops();

        let mut diff_nibbles: Vec<u8> = Vec::new();
        let mut extra_nibbles: Vec<u8> = Vec::new();

        if n_cigar == 0 {
            // Unmapped / no CIGAR: all bases go to extra
            extra_nibbles.extend_from_slice(&nibbles);
        } else {
            let mut ref_pos = pos;
            let mut query_pos = 0usize;

            for &(op, len) in &cigar_ops {
                let len = len as usize;
                if is_alignment_op(op) {
                    // M/=/X: XOR with consensus
                    for _ in 0..len {
                        if query_pos < nibbles.len() {
                            let cons_base = consensus.get(ref_id, ref_pos);
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
        streams.seq_diff.extend_from_slice(&pack_nibble_pairs(&diff_nibbles));

        // Write extra (varint count + packed nibble bytes)
        write_varint(&mut streams.seq_extra, extra_nibbles.len());
        streams.seq_extra.extend_from_slice(&pack_nibble_pairs(&extra_nibbles));

        // Collect per-record quality (Phred+33 ASCII) and ASCII sequences for quality_ctx
        let qual_raw = rec.qual_bytes();
        let mut qual_ascii = Vec::with_capacity(qual_raw.len());
        for &q in qual_raw {
            if q == 0xFF {
                has_unavailable_quality = true;
                qual_ascii.push(33u8); // map unavailable to Phred 0
            } else {
                qual_ascii.push(q + 33);
            }
        }
        qualities_ascii.push(qual_ascii);

        let ascii_seq: Vec<u8> = nibbles.iter().map(|&n| nibble_to_ascii(n)).collect();
        sequences_ascii.push(ascii_seq);
    }

    StreamsResult {
        streams,
        quality_data: ChunkQualityData {
            qualities_ascii,
            sequences_ascii,
            has_unavailable_quality,
        },
    }
}

/// Reconstruct sequence nibbles from consensus + CIGAR + diff/extra nibbles.
/// Used during decompression to get sequences before quality_ctx decompression.
pub fn reconstruct_sequence_nibbles(
    consensus: &ChunkConsensus,
    ref_id: i32,
    pos: i32,
    cigar_bytes: &[u8],
    n_cigar_op: usize,
    diff_nibbles: &[u8],
    extra_nibbles: &[u8],
    l_seq: usize,
) -> Vec<u8> {
    let mut seq_nibbles: Vec<u8> = Vec::with_capacity(l_seq);

    if n_cigar_op == 0 {
        seq_nibbles.extend_from_slice(&extra_nibbles[..l_seq]);
    } else {
        let cigar_ops = parse_cigar_ops(cigar_bytes, n_cigar_op);
        let mut diff_pos = 0usize;
        let mut extra_pos = 0usize;
        let mut ref_pos = pos;

        // Cache segment lookup: 1 HashMap lookup per record instead of per base
        let segment = consensus.get_segment(ref_id);

        for (op, len) in cigar_ops {
            let len = len as usize;
            if is_alignment_op(op) {
                for _ in 0..len {
                    let cons_base = segment.map_or(0, |seg| seg.get(ref_pos));
                    let base = diff_nibbles[diff_pos] ^ cons_base;
                    seq_nibbles.push(base);
                    diff_pos += 1;
                    ref_pos += 1;
                }
            } else if op == CIGAR_I || op == CIGAR_S {
                for _ in 0..len {
                    seq_nibbles.push(extra_nibbles[extra_pos]);
                    extra_pos += 1;
                }
            } else if op == CIGAR_D || op == CIGAR_N {
                ref_pos += len as i32;
            }
        }
    }

    seq_nibbles
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
) -> Vec<u8> {
    // Reconstruct sequence from diff + extra + consensus + CIGAR
    let seq_nibbles = reconstruct_sequence_nibbles(
        consensus, ref_id, pos, cigar_bytes, n_cigar_op as usize,
        diff_nibbles, extra_nibbles, l_seq,
    );

    let seq_packed = pack_nibbles(&seq_nibbles);

    // Assemble record bytes
    let record_len = 32 + read_name.len() + cigar_bytes.len() + seq_packed.len() + qual_bytes.len() + aux_bytes.len();
    let mut data = Vec::with_capacity(record_len);

    data.extend_from_slice(&ref_id.to_le_bytes());       // 0..4
    data.extend_from_slice(&pos.to_le_bytes());           // 4..8
    data.push(read_name.len() as u8);                     // 8
    data.push(mapq);                                       // 9
    data.extend_from_slice(&bin_val.to_le_bytes());       // 10..12
    data.extend_from_slice(&n_cigar_op.to_le_bytes());    // 12..14
    data.extend_from_slice(&flag.to_le_bytes());          // 14..16
    data.extend_from_slice(&(l_seq as i32).to_le_bytes()); // 16..20
    data.extend_from_slice(&next_ref_id.to_le_bytes());   // 20..24
    data.extend_from_slice(&next_pos.to_le_bytes());      // 24..28
    data.extend_from_slice(&tlen.to_le_bytes());          // 28..32
    data.extend_from_slice(read_name);
    data.extend_from_slice(cigar_bytes);
    data.extend_from_slice(&seq_packed);
    data.extend_from_slice(qual_bytes);
    data.extend_from_slice(aux_bytes);

    data
}

/// Parse CIGAR bytes into (op, len) pairs.
pub fn parse_cigar_ops(bytes: &[u8], n_ops: usize) -> Vec<(u8, u32)> {
    let mut ops = Vec::with_capacity(n_ops);
    for i in 0..n_ops {
        let val = u32::from_le_bytes(bytes[i * 4..(i + 1) * 4].try_into().unwrap());
        ops.push(((val & 0xF) as u8, val >> 4));
    }
    ops
}
