//! Covered-interval coordinate map for the reference codec.
//!
//! The reference codec's coordinate system maps `(ref_id, ref_pos)` to a flat
//! "base index" into the concatenated covered-interval consensus. `IntervalMap`
//! is that translation plus its validation. R2 read positions are stored as a
//! signed zigzag delta from the mate's R1 base index (see `write_svarint`).
//!
//! `deserialize` is hostile-input-facing (fuzzed in a later task): it bounds the
//! interval count before allocating, requires exact byte consumption, and
//! re-validates ordering/overlap via `from_intervals`.

use super::super::dna_utils::{read_varint, write_varint};
use anyhow::{Result, bail};

/// A contiguous covered region on one reference sequence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Interval {
    pub ref_id: u32,
    pub ref_start: u64,
    pub length: u64,
}

/// Translation from `(ref_id, ref_pos)` to a flat base index over the
/// concatenated covered intervals.
///
/// Intervals are stored strictly ordered by `(ref_id, ref_start)` with no
/// same-reference overlap. `prefix[i]` is the base offset where interval `i`
/// begins in the flat coordinate space; `prefix[len]` (the final guard entry)
/// is the total base count.
pub struct IntervalMap {
    intervals: Vec<Interval>,
    /// `prefix[i]` = base start of `intervals[i]`; final entry = total bases.
    prefix: Vec<u64>,
}

impl IntervalMap {
    /// Build from intervals that must be STRICTLY ordered by `(ref_id,
    /// ref_start)` with no same-reference overlap. Errors on disorder, overlap,
    /// or coordinate overflow.
    pub fn from_intervals(intervals: Vec<Interval>) -> Result<Self> {
        let mut prefix = Vec::with_capacity(intervals.len() + 1);
        let mut acc: u64 = 0;
        for (i, cur) in intervals.iter().enumerate() {
            // Reject coordinate overflow on ref_start + length.
            cur.ref_start
                .checked_add(cur.length)
                .ok_or_else(|| anyhow::anyhow!("interval {i} ref_start+length overflows u64"))?;
            if i > 0 {
                let prev = &intervals[i - 1];
                // Strict ordering by (ref_id, ref_start).
                if (prev.ref_id, prev.ref_start) >= (cur.ref_id, cur.ref_start) {
                    bail!(
                        "intervals not strictly ordered by (ref_id, ref_start): \
                         interval {} ({},{}) >= interval {} ({},{})",
                        i - 1,
                        prev.ref_id,
                        prev.ref_start,
                        i,
                        cur.ref_id,
                        cur.ref_start
                    );
                }
                // No same-reference overlap.
                if prev.ref_id == cur.ref_id && prev.ref_start + prev.length > cur.ref_start {
                    bail!(
                        "overlapping intervals on ref {}: [{},{}) overlaps start {}",
                        prev.ref_id,
                        prev.ref_start,
                        prev.ref_start + prev.length,
                        cur.ref_start
                    );
                }
            }
            prefix.push(acc);
            acc = acc
                .checked_add(cur.length)
                .ok_or_else(|| anyhow::anyhow!("total base count overflows u64"))?;
        }
        prefix.push(acc); // final guard entry = total bases
        Ok(IntervalMap { intervals, prefix })
    }

    /// The covered intervals, ordered by `(ref_id, ref_start)`. Used by Task 12
    /// to build the informational reference_meta blob from a finalized map.
    pub fn intervals(&self) -> &[Interval] {
        &self.intervals
    }

    /// Total number of bases across all covered intervals.
    pub fn total_bases(&self) -> u64 {
        *self.prefix.last().unwrap_or(&0)
    }

    /// Translate `(ref_id, ref_pos)` to its flat base index, or `None` if the
    /// position is not inside any covered interval.
    pub fn ref_to_base(&self, ref_id: u32, ref_pos: u64) -> Option<u64> {
        // Binary search for the last interval whose (ref_id, ref_start) <=
        // (ref_id, ref_pos), then check it actually covers ref_pos.
        let idx = self
            .intervals
            .partition_point(|iv| (iv.ref_id, iv.ref_start) <= (ref_id, ref_pos));
        if idx == 0 {
            return None;
        }
        let iv = &self.intervals[idx - 1];
        if iv.ref_id == ref_id && ref_pos >= iv.ref_start && ref_pos < iv.ref_start + iv.length {
            Some(self.prefix[idx - 1] + (ref_pos - iv.ref_start))
        } else {
            None
        }
    }

    /// True iff `[base, base + len)` lies within a single interval.
    ///
    /// An empty span (`len == 0`) is true iff `base` falls inside some interval.
    pub fn span_in_single_interval(&self, base: u64, len: u64) -> bool {
        // Find interval i with prefix[i] <= base < prefix[i] + length.
        // prefix is strictly increasing per interval (zero-length intervals
        // cannot occur because length==0 would make ordering checks pass but
        // partition_point still resolves; guard with the explicit range test).
        let idx = self.prefix.partition_point(|&p| p <= base);
        if idx == 0 || idx > self.intervals.len() {
            return false;
        }
        let i = idx - 1;
        let iv_start = self.prefix[i];
        let iv_end = iv_start + self.intervals[i].length;
        if base >= iv_end {
            return false;
        }
        match base.checked_add(len) {
            Some(end) => end <= iv_end,
            None => false,
        }
    }

    /// Serialize to bytes: varint count, then per interval `ref_id`,
    /// `ref_start`, `length` (all varint).
    #[allow(dead_code)] // used in tests + planned for future archive verify
    pub fn serialize(&self) -> Vec<u8> {
        let mut out = Vec::new();
        write_varint(&mut out, self.intervals.len() as u64);
        for iv in &self.intervals {
            write_varint(&mut out, iv.ref_id as u64);
            write_varint(&mut out, iv.ref_start);
            write_varint(&mut out, iv.length);
        }
        out
    }

    /// Parse bytes written by `serialize`. Bounds the interval count before
    /// allocating, requires exact byte consumption, and re-validates
    /// ordering/overlap via `from_intervals` (so a serialized overlap or
    /// disorder fails to deserialize).
    #[allow(dead_code)] // used in tests + planned for future archive verify
    pub fn deserialize(data: &[u8]) -> Result<Self> {
        let mut offset = 0usize;
        let count = read_varint(data, &mut offset)
            .ok_or_else(|| anyhow::anyhow!("IntervalMap: truncated count"))?;
        // Each interval is at least 3 bytes (three single-byte varints), so a
        // count larger than the remaining bytes / 3 cannot be honest.
        let max_count = (data.len() - offset) as u64 / 3;
        if count > max_count {
            bail!("IntervalMap: claimed {count} intervals but only room for {max_count}");
        }
        let mut intervals = Vec::with_capacity(count as usize);
        for _ in 0..count {
            let ref_id = read_varint(data, &mut offset)
                .ok_or_else(|| anyhow::anyhow!("IntervalMap: truncated ref_id"))?;
            let ref_id: u32 = ref_id
                .try_into()
                .map_err(|_| anyhow::anyhow!("IntervalMap: ref_id {ref_id} exceeds u32"))?;
            let ref_start = read_varint(data, &mut offset)
                .ok_or_else(|| anyhow::anyhow!("IntervalMap: truncated ref_start"))?;
            let length = read_varint(data, &mut offset)
                .ok_or_else(|| anyhow::anyhow!("IntervalMap: truncated length"))?;
            intervals.push(Interval {
                ref_id,
                ref_start,
                length,
            });
        }
        if offset != data.len() {
            bail!(
                "IntervalMap: {} trailing bytes after {} intervals",
                data.len() - offset,
                count
            );
        }
        Self::from_intervals(intervals)
    }
}

/// One window's (or the whole archive's) finalized consensus, stored as a 2-bit
/// `ACGT` buffer plus a separate 1-bit-per-base `N` bitmap (spec §5.2).
///
/// ## Bit conventions (locked — Task 6 framing & Task 13 decode must match)
/// - **2-bit packed**: base codes `A=0, C=1, G=2, T=3`. Base `i` occupies bits
///   `2*(i%4) .. 2*(i%4)+2` of byte `i/4`, **LSB-first** (base 0 in the two
///   lowest bits). `N` positions hold the placeholder code `0` (`A`).
/// - **N-bitmap**: bit `i` is byte `i/8`, bit `i%8`, **LSB-first** (set ⇒
///   position `i` is `N`).
pub struct PackedBacking {
    /// 2-bit ACGT, 4 bases/byte, LSB-first; N slots hold placeholder 0.
    packed: Vec<u8>,
    /// 1 bit/base, LSB-first; bit set ⇒ that base is N.
    n_bitmap: Vec<u8>,
    len: usize,
}

impl PackedBacking {
    /// Number of consensus bases.
    pub fn len(&self) -> usize {
        self.len
    }

    /// True iff this consensus has zero bases.
    #[allow(dead_code)] // used in tests
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// True iff position `i` is an `N` (its N-bitmap bit is set).
    pub fn is_n(&self, i: usize) -> bool {
        debug_assert!(
            i < self.len,
            "is_n index {i} out of bounds (len {})",
            self.len
        );
        (self.n_bitmap[i / 8] >> (i % 8)) & 1 == 1
    }

    /// The 2-bit base code at position `i` (without N overlay).
    #[inline]
    fn code_at(&self, i: usize) -> u8 {
        (self.packed[i / 4] >> (2 * (i % 4))) & 0b11
    }

    /// Reconstruct the full consensus as ASCII `ACGTN` bytes (N overlaid where
    /// the bitmap bit is set).
    #[allow(dead_code)] // used in tests
    pub fn unpack(&self) -> Vec<u8> {
        self.slice(0, self.len)
    }

    /// Reconstruct the sub-range `[base, base+len)` as ASCII `ACGTN` bytes with
    /// `N` overlaid from the bitmap. Caller guarantees `base + len <= self.len`.
    pub fn slice(&self, base: usize, len: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(len);
        self.fill_slice(base, len, &mut out);
        out
    }

    /// Append the ASCII `ACGTN` bytes of `[base, base+len)` to `out` (N overlaid
    /// from the bitmap). Lets a caller reconstruct a read directly into an owned
    /// buffer with no intermediate `slice()` allocation. Caller guarantees
    /// `base + len <= self.len`.
    pub fn fill_slice(&self, base: usize, len: usize, out: &mut Vec<u8>) {
        debug_assert!(
            base + len <= self.len,
            "slice [{base},{}) out of bounds (len {})",
            base + len,
            self.len
        );
        const LUT: [u8; 4] = [b'A', b'C', b'G', b'T'];
        out.reserve(len);
        for i in base..base + len {
            if self.is_n(i) {
                out.push(b'N');
            } else {
                out.push(LUT[self.code_at(i) as usize]);
            }
        }
    }

    /// Raw 2-bit packed bytes (for Task 6 block framing).
    #[allow(dead_code)] // used in tests
    pub fn packed_bytes(&self) -> &[u8] {
        &self.packed
    }

    /// Raw N-bitmap bytes (for Task 6 block framing).
    #[allow(dead_code)] // used in tests
    pub fn n_bitmap_bytes(&self) -> &[u8] {
        &self.n_bitmap
    }

    /// Rebuild from decoded block bytes (Task 6 decode path). `packed` must be
    /// `len.div_ceil(4)` bytes and `n_bitmap` `len.div_ceil(8)` bytes.
    pub fn from_parts(packed: Vec<u8>, n_bitmap: Vec<u8>, len: usize) -> Self {
        debug_assert_eq!(packed.len(), len.div_ceil(4), "from_parts packed length");
        debug_assert_eq!(
            n_bitmap.len(),
            len.div_ceil(8),
            "from_parts n_bitmap length"
        );
        PackedBacking {
            packed,
            n_bitmap,
            len,
        }
    }
}

// ---------------------------------------------------------------------------
// Block-streamed packed-consensus + N-bitmap framing (Task 6).
//
// The consensus can be gigabytes, so it is stored as a sequence of
// independently-BSC-compressed blocks: decode can stream/seek block-by-block
// (Task 13 reads each block from the file at its offset). The 2-bit packed
// stream and the 1-bit N-bitmap use the SAME v4 block framing, differing only
// in their packing density (4 bases/byte vs 8 bases/byte).
//
// ## On-wire block-stream layout
// ```
// [num_blocks: u32 LE]
// repeated num_blocks times:
//     [v4 block header: 12 B]   (block_len LE u32 | record_count LE u32 | crc LE u32)
//     [payload: block_len bytes] (BSC-compressed packed/bitmap bytes for the block)
// ```
// `record_count` holds the **base count** (packed stream) or **bit/base count**
// (bitmap stream) of that block — NOT the compressed byte length — so the sum of
// all `record_count`s equals the total base count (spec §7.1; this is what the
// footer's Consensus/NBitmap `record_count` carries, and Task 3's validator
// asserts `NBitmap.record_count == Consensus.record_count`).
//
// ## Packed-only vs bitmap-overlay decode contract
// `read_packed_backing_blocks` reconstructs ONLY the 2-bit packed bytes and
// returns a `PackedBacking` with an **all-zero N-bitmap** (no N positions).
// The N-bitmap is framed and read SEPARATELY via `write_bitmap_blocks` /
// `read_bitmap_blocks`. Task 13 must overlay the two: read the packed stream and
// the bitmap stream, then build the final consensus via
// `read_packed_backing_with_bitmap` (or `PackedBacking::from_parts` directly).
// ---------------------------------------------------------------------------

use super::super::bsc;

/// On-disk size of one v4 block header (`block_len u32 | record_count u32 | crc
/// u32`). Used as the per-block lower bound for the `num_blocks` DoS guard.
const V4_BLOCK_HEADER_LEN: usize = 12;

/// Write the 2-bit packed consensus as a stream of BSC-compressed v4 blocks.
///
/// `block_bases` is the number of bases per block (last block = remainder); it
/// **must** be a positive multiple of 4 so each block starts/ends on a packed
/// byte boundary (4 bases/byte) — asserted, never silently misaligned.
///
/// Layout: a `[num_blocks: u32 LE]` prefix, then per block a 12-byte v4 header
/// (`block_len`, `record_count` = base count, crc) followed by the BSC payload.
#[allow(dead_code)] // used in tests
pub fn write_packed_backing_blocks(
    out: &mut Vec<u8>,
    c: &PackedBacking,
    block_bases: usize,
) -> Result<()> {
    assert!(block_bases > 0, "block_bases must be > 0");
    assert!(
        block_bases.is_multiple_of(4),
        "block_bases must be a multiple of 4 (got {block_bases})"
    );
    let len = c.len();
    let packed = c.packed_bytes();
    let num_blocks = len.div_ceil(block_bases) as u32;
    out.extend_from_slice(&num_blocks.to_le_bytes());
    let mut base = 0usize;
    while base < len {
        let block_base_count = block_bases.min(len - base);
        let byte_start = base / 4;
        let byte_end = (base + block_base_count).div_ceil(4);
        let block_packed = &packed[byte_start..byte_end];
        let comp = bsc::compress(block_packed)?;
        bsc::write_block_frame_header(out, comp.len() as u32, block_base_count as u32, &comp);
        out.extend_from_slice(&comp);
        base += block_base_count;
    }
    // Edge case: an empty consensus (len == 0) writes num_blocks = 0 and no
    // blocks; decode of total_bases == 0 reconstructs an empty packed buffer.
    Ok(())
}

/// Read a packed-consensus block stream written by `write_packed_backing_blocks`.
///
/// Returns a `PackedBacking` carrying ONLY the 2-bit packed bytes with an
/// all-zero N-bitmap (no N positions); pair it with a separately-read bitmap
/// (see module docs) to recover N positions.
///
/// Hostile-input-facing: bounds every read against the remaining bytes, verifies
/// each block's CRC, checks per-block decoded length == `ceil(record_count/4)`,
/// requires Σ record_counts == `total_bases`, and requires exact consumption.
#[allow(dead_code)] // used in tests
pub fn read_packed_backing_blocks(data: &[u8], total_bases: usize) -> Result<PackedBacking> {
    let packed = read_blocks_into(data, total_bases, /*bits_per_record*/ false)?;
    Ok(PackedBacking::from_parts(
        packed,
        vec![0u8; total_bases.div_ceil(8)],
        total_bases,
    ))
}

/// Write a 1-bit-per-base bitmap (e.g. the N-bitmap) as BSC-compressed v4 blocks.
///
/// `block_bits` is the number of bits/bases per block (last block = remainder);
/// it **must** be a positive multiple of 8 so each block starts/ends on a byte
/// boundary (8 bits/byte) — asserted. `record_count` per block = bit (base)
/// count, so Σ record_counts == `total_bits`. `bitmap_bytes` must be exactly
/// `total_bits.div_ceil(8)` bytes.
#[allow(dead_code)] // used in tests
pub fn write_bitmap_blocks(
    out: &mut Vec<u8>,
    bitmap_bytes: &[u8],
    total_bits: usize,
    block_bits: usize,
) -> Result<()> {
    assert!(block_bits > 0, "block_bits must be > 0");
    assert!(
        block_bits.is_multiple_of(8),
        "block_bits must be a multiple of 8 (got {block_bits})"
    );
    assert_eq!(
        bitmap_bytes.len(),
        total_bits.div_ceil(8),
        "bitmap_bytes length mismatch"
    );
    let num_blocks = total_bits.div_ceil(block_bits) as u32;
    out.extend_from_slice(&num_blocks.to_le_bytes());
    let mut bit = 0usize;
    while bit < total_bits {
        let block_bit_count = block_bits.min(total_bits - bit);
        let byte_start = bit / 8;
        let byte_end = (bit + block_bit_count).div_ceil(8);
        let block_bytes = &bitmap_bytes[byte_start..byte_end];
        let comp = bsc::compress(block_bytes)?;
        bsc::write_block_frame_header(out, comp.len() as u32, block_bit_count as u32, &comp);
        out.extend_from_slice(&comp);
        bit += block_bit_count;
    }
    Ok(())
}

/// Read a bitmap block stream written by `write_bitmap_blocks`, returning the
/// reassembled `total_bits.div_ceil(8)`-byte bitmap.
///
/// Same defensive guarantees as `read_packed_backing_blocks`, but 1 bit/base:
/// per-block decoded length == `ceil(record_count/8)`, Σ record_counts ==
/// `total_bits`, bounded reads, exact consumption, CRC verified per block.
pub fn read_bitmap_blocks(data: &[u8], total_bits: usize) -> Result<Vec<u8>> {
    read_blocks_into(data, total_bits, /*bits_per_record*/ true)
}

/// Build a `PackedBacking` by overlaying a separately-framed N-bitmap onto a
/// separately-framed packed stream (the Task 13 combiner). Decodes both streams
/// and pairs them; the bitmap's total bit count must equal the packed stream's
/// base count (`total_bases`).
pub fn read_packed_backing_with_bitmap(
    packed_data: &[u8],
    bitmap_data: &[u8],
    total_bases: usize,
) -> Result<PackedBacking> {
    let packed = read_blocks_into(packed_data, total_bases, /*bits_per_record*/ false)?;
    let bitmap = read_bitmap_blocks(bitmap_data, total_bases)?;
    Ok(PackedBacking::from_parts(packed, bitmap, total_bases))
}

/// Shared block-stream reader for both the 2-bit packed stream
/// (`bits_per_record == false`, 4 records/byte) and the 1-bit bitmap stream
/// (`bits_per_record == true`, 8 records/byte). Returns the reassembled raw
/// byte buffer (`total_records.div_ceil(records_per_byte)` bytes).
fn read_blocks_into(data: &[u8], total_records: usize, bits_per_record: bool) -> Result<Vec<u8>> {
    let records_per_byte = if bits_per_record { 8usize } else { 4usize };
    let mut cursor: &[u8] = data;
    let hdr_blocks = read_u32_le(&mut cursor)?;
    let num_blocks = hdr_blocks as usize;
    // Upfront DoS bound: every block carries at least a 12-byte v4 header, so a
    // declared `num_blocks` larger than `remaining / 12` cannot be honest. Reject
    // before the loop (which would otherwise iterate up to ~4B times on a hostile
    // num_blocks before the cursor runs dry).
    if num_blocks > cursor.len() / V4_BLOCK_HEADER_LEN {
        bail!(
            "consensus block stream: declared {num_blocks} blocks exceeds {} possible in {} remaining bytes",
            cursor.len() / V4_BLOCK_HEADER_LEN,
            cursor.len()
        );
    }

    // `total_records` derives from an untrusted footer count (the v5 footer CRC is
    // recomputable, not a MAC). BEFORE committing the `out_len`-sized reassembly
    // buffer, walk the block HEADERS ONLY (no decompression) to prove the declared
    // `total_records` is actually achievable by the blocks present:
    //   - each block's payload fits in the remaining bytes (framing is sound),
    //   - each block's decoded size (`block_records / records_per_byte`) is within
    //     the single-block BSC decode cap, so one 12-byte header can't claim
    //     billions of records, and
    //   - Σ block record_counts == total_records exactly.
    // This bounds `out_len` by `num_blocks * BSC_MAX_BLOCK_LEN` and rejects a footer
    // that declares a large `total_records` with too few backing blocks *before* any
    // large allocation. Without it, `try_reserve_exact` could succeed under memory
    // overcommit and the subsequent `resize` would fault every page → OOM-kill.
    // (The remaining O(decoded-size) materialize-all cost is the documented
    // bounded-RAM decode follow-up.)
    let block_decode_cap = crate::compression::decompress_impl::BSC_MAX_BLOCK_LEN;
    {
        let mut scan: &[u8] = cursor;
        let mut acc: usize = 0;
        for b in 0..num_blocks {
            let hdr = bsc::read_block_frame_header(&mut scan)?;
            let block_len = hdr.block_len as usize;
            if block_len > scan.len() {
                bail!(
                    "consensus block {b}: payload len {block_len} exceeds {} remaining bytes",
                    scan.len()
                );
            }
            let block_records = hdr.record_count as usize;
            if block_records.div_ceil(records_per_byte) > block_decode_cap {
                bail!(
                    "consensus block {b}: declared {block_records} records exceed the single-block decode cap"
                );
            }
            scan = &scan[block_len..];
            acc = acc
                .checked_add(block_records)
                .ok_or_else(|| anyhow::anyhow!("consensus block {b}: record count overflow"))?;
            if acc > total_records {
                bail!(
                    "consensus block {b}: cumulative record count {acc} exceeds total {total_records}"
                );
            }
        }
        if !scan.is_empty() {
            bail!(
                "consensus block stream: {} trailing bytes after {num_blocks} blocks",
                scan.len()
            );
        }
        if acc != total_records {
            bail!("consensus block stream: Σ record counts {acc} != total {total_records}");
        }
    }

    // Now `total_records` is validated against the present blocks; allocate fallibly
    // anyway (belt-and-suspenders against a still-large but block-backed count).
    let out_len = total_records.div_ceil(records_per_byte);
    let mut out: Vec<u8> = Vec::new();
    out.try_reserve_exact(out_len).map_err(|_| {
        anyhow::anyhow!(
            "consensus stream: {out_len}-byte buffer alloc failed (hostile record_count?)"
        )
    })?;
    out.resize(out_len, 0u8);
    let mut record_acc: usize = 0;
    for b in 0..num_blocks {
        let hdr = bsc::read_block_frame_header(&mut cursor)?;
        let block_len = hdr.block_len as usize;
        if block_len > cursor.len() {
            bail!(
                "consensus block {b}: payload len {block_len} exceeds {} remaining bytes",
                cursor.len()
            );
        }
        let payload = &cursor[..block_len];
        cursor = &cursor[block_len..];
        bsc::verify_block_frame_crc(hdr.expected_crc, hdr.record_count, payload)?;
        let block_records = hdr.record_count as usize;
        // Σ record_counts must not exceed total.
        let new_acc = record_acc
            .checked_add(block_records)
            .ok_or_else(|| anyhow::anyhow!("consensus block {b}: record count overflow"))?;
        if new_acc > total_records {
            bail!(
                "consensus block {b}: cumulative record count {new_acc} exceeds total {total_records}"
            );
        }
        let decoded = bsc::decompress(payload)?;
        let expected_bytes = block_records.div_ceil(records_per_byte);
        if decoded.len() != expected_bytes {
            bail!(
                "consensus block {b}: decoded {} bytes, expected {expected_bytes} for {block_records} records",
                decoded.len()
            );
        }
        // Because records_per_byte divides block boundaries (block_bases%4==0 /
        // block_bits%8==0), record_acc is a byte-aligned offset for all but the
        // final (remainder) block, where the final partial byte lands exactly at
        // the end of `out`.
        let byte_off = record_acc / records_per_byte;
        out[byte_off..byte_off + decoded.len()].copy_from_slice(&decoded);
        record_acc = new_acc;
    }
    if record_acc != total_records {
        bail!("consensus block stream: Σ record counts {record_acc} != total {total_records}");
    }
    if !cursor.is_empty() {
        bail!(
            "consensus block stream: {} trailing bytes after {num_blocks} blocks",
            cursor.len()
        );
    }
    Ok(out)
}

/// Fast-verify the `[num_blocks u32 LE][v4 header 12B][payload]...` framing used
/// by `write_packed_backing_blocks` (Consensus) and `write_bitmap_blocks`
/// (NBitmap) WITHOUT bsc-decoding: read each block's v4 header, verify the
/// per-block CRC32 over `record_count || payload`, and require exact consumption.
/// Returns `(num_blocks_verified, Σ record_count)`. This is the consensus/bitmap
/// analogue of `streams::decode_block_payloads` — same CRC scope, no codec decode.
pub fn verify_block_stream_crc(data: &[u8]) -> Result<(u32, u64)> {
    let mut cursor: &[u8] = data;
    let num_blocks = read_u32_le(&mut cursor)?;
    // Upfront DoS bound: each block is at least a 12-byte v4 header, so a declared
    // num_blocks larger than `remaining / 12` is impossible — reject before the
    // loop (mirrors `read_blocks_into`).
    if num_blocks as usize > cursor.len() / V4_BLOCK_HEADER_LEN {
        bail!(
            "block stream: declared {num_blocks} blocks exceeds {} possible in {} remaining bytes",
            cursor.len() / V4_BLOCK_HEADER_LEN,
            cursor.len()
        );
    }
    let mut blocks: u32 = 0;
    let mut record_sum: u64 = 0;
    for b in 0..num_blocks {
        let hdr = bsc::read_block_frame_header(&mut cursor)?;
        let block_len = hdr.block_len as usize;
        if block_len > cursor.len() {
            bail!(
                "block stream block {b}: payload len {block_len} exceeds {} remaining bytes",
                cursor.len()
            );
        }
        let payload = &cursor[..block_len];
        cursor = &cursor[block_len..];
        bsc::verify_block_frame_crc(hdr.expected_crc, hdr.record_count, payload)?;
        blocks += 1;
        record_sum = record_sum
            .checked_add(hdr.record_count as u64)
            .ok_or_else(|| anyhow::anyhow!("block stream block {b}: record count overflow"))?;
    }
    if !cursor.is_empty() {
        bail!(
            "block stream: {} trailing bytes after {num_blocks} blocks",
            cursor.len()
        );
    }
    Ok((blocks, record_sum))
}

/// Read a little-endian u32 from the front of `cursor`, advancing it.
fn read_u32_le(cursor: &mut &[u8]) -> Result<u32> {
    if cursor.len() < 4 {
        bail!("consensus block stream: truncated num_blocks prefix");
    }
    let v = u32::from_le_bytes([cursor[0], cursor[1], cursor[2], cursor[3]]);
    *cursor = &cursor[4..];
    Ok(v)
}

// ---------------------------------------------------------------------------
// Streaming backing writers (Task 9).
//
// These emit the SAME `[num_blocks u32 LE][v4 header 12B][BSC payload]...`
// wire format as `write_packed_backing_blocks` / `write_bitmap_blocks`, but
// without materializing a full `PackedBacking` — they walk the `IntervalMap`
// intervals one block at a time, gathering normalized reference bytes directly
// from the `RefSequence` slices.
//
// "Normalized" here means: the byte is taken verbatim from the reference
// sequence (ASCII ACGTN). Non-ACGTN bytes are treated as N.
// ---------------------------------------------------------------------------

/// Block size (in bases) for the packed consensus + N-bitmap streams written
/// into an archive. Must be a multiple of 4 (4 bases/byte packed) and 8 (8
/// bits/byte bitmap); `100_000_000` satisfies both.
///
/// Decode ceiling: a packed block decodes to `block_bases / 4` bytes, which must
/// stay under `BSC_MAX_BLOCK_LEN` (64 MiB + 4096) on the decode side. So
/// `BACKING_BLOCK_BASES` must also NOT exceed ~268M bases (`4 * BSC_MAX_BLOCK_LEN`).
/// At the current 100M a packed block is ~25 MiB — well under the cap.
pub const BACKING_BLOCK_BASES: usize = 100_000_000;

/// Per-block base count for the reference-aware MERGE's parallel backing/N-bitmap
/// re-derivation. Much smaller than `BACKING_BLOCK_BASES` ON PURPOSE: the merge runs
/// serially after the shard workers finish (all cores idle), and the single
/// 100M-base block would BSC-compress single-threaded — the dominant merge tax. At
/// 8M bases (~2 MiB packed) a chr20-sized union (~64 Mbp) splits into ~8 independent
/// blocks that rayon compresses concurrently, cutting that tax ~8×. The smaller
/// block trades a negligible ratio loss (2-bit DNA is already near-incompressible by
/// BSC) for the parallelism, and the per-block framing is internal to the global so
/// the archive still decodes identically. Multiple of 8 (N-bitmap byte alignment).
pub const MERGE_BACKING_BLOCK_BASES: usize = 8_000_000;

/// Normalize a raw reference byte to an ACGTN upper-case byte. Bytes outside
/// the ACGTN alphabet (e.g. lower-case, IUPAC ambiguity codes) map to `N`.
#[inline]
fn normalize_ref_byte(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        v @ (b'A' | b'C' | b'G' | b'T') => v,
        _ => b'N',
    }
}

/// An iterator that yields all normalized reference bases across the
/// `IntervalMap`'s intervals in flat order, one byte at a time.
struct FlatRefIter<'a> {
    refs: &'a [strobealign::io::fasta::RefSequence],
    intervals: &'a [Interval],
    iv_idx: usize, // current interval index
    iv_off: u64,   // offset within the current interval
}

impl<'a> FlatRefIter<'a> {
    fn new(refs: &'a [strobealign::io::fasta::RefSequence], imap: &'a IntervalMap) -> Self {
        FlatRefIter {
            refs,
            intervals: imap.intervals(),
            iv_idx: 0,
            iv_off: 0,
        }
    }

    /// Position the iterator at flat base `start` (`0 ≤ start ≤ total_bases`): the
    /// next `next()` yields the base at flat index `start`. Lets the parallel backing
    /// builders walk a single block's base range without re-scanning from 0. Same
    /// module as `IntervalMap`, so it reads `prefix` directly.
    fn new_at(
        refs: &'a [strobealign::io::fasta::RefSequence],
        imap: &'a IntervalMap,
        start: u64,
    ) -> Self {
        let intervals = imap.intervals();
        // prefix[i] = flat start of interval i (final guard entry = total_bases).
        // The interval containing `start` is the last one whose prefix ≤ start.
        let idx = imap.prefix.partition_point(|&p| p <= start);
        if idx == 0 {
            return FlatRefIter { refs, intervals, iv_idx: 0, iv_off: 0 };
        }
        let i = idx - 1;
        if i >= intervals.len() {
            // start == total_bases → positioned at end (next() yields None).
            return FlatRefIter { refs, intervals, iv_idx: intervals.len(), iv_off: 0 };
        }
        FlatRefIter { refs, intervals, iv_idx: i, iv_off: start - imap.prefix[i] }
    }
}

impl Iterator for FlatRefIter<'_> {
    type Item = u8;
    #[inline]
    fn next(&mut self) -> Option<u8> {
        // Advance past exhausted intervals.
        while self.iv_idx < self.intervals.len()
            && self.iv_off >= self.intervals[self.iv_idx].length
        {
            self.iv_idx += 1;
            self.iv_off = 0;
        }
        if self.iv_idx >= self.intervals.len() {
            return None;
        }
        let iv = &self.intervals[self.iv_idx];
        let pos = (iv.ref_start + self.iv_off) as usize;
        let seq = &self.refs[iv.ref_id as usize].sequence;
        let raw = if pos < seq.len() { seq[pos] } else { b'N' };
        self.iv_off += 1;
        Some(normalize_ref_byte(raw))
    }
}

/// Stream the 2-bit packed backing block-by-block from normalized reference
/// intervals into `w`. Holds ~one block at a time.
///
/// The emitted bytes are byte-identical to what `write_packed_backing_blocks`
/// would produce for the same bases, using the same `[num_blocks u32 LE]` then
/// per-block `[v4 header 12B][BSC payload]` framing.
///
/// Returns the total number of bases written (= `imap.total_bases()`).
pub fn stream_packed_backing(
    w: &mut impl std::io::Write,
    refs: &[strobealign::io::fasta::RefSequence],
    imap: &IntervalMap,
    block_bases: usize,
) -> anyhow::Result<u64> {
    assert!(block_bases > 0, "block_bases must be > 0");
    assert!(
        block_bases.is_multiple_of(4),
        "block_bases must be a multiple of 4 (got {block_bases})"
    );
    // FlatRefIter indexes refs[iv.ref_id]; reject an out-of-range ref_id up front
    // (clean error instead of an index panic on a misused interval set).
    if !imap
        .intervals()
        .iter()
        .all(|iv| (iv.ref_id as usize) < refs.len())
    {
        bail!("interval ref_id out of range for reference set");
    }
    let total_bases = imap.total_bases();
    let num_blocks = total_bases.div_ceil(block_bases as u64) as u32;
    w.write_all(&num_blocks.to_le_bytes())?;

    let mut iter = FlatRefIter::new(refs, imap);
    let mut bases_remaining = total_bases as usize;
    while bases_remaining > 0 {
        let block_base_count = block_bases.min(bases_remaining);
        // Pack block_base_count bases into 2-bit bytes.
        let packed_bytes_count = block_base_count.div_ceil(4);
        let mut packed = vec![0u8; packed_bytes_count];
        for i in 0..block_base_count {
            let b = iter.next().unwrap_or(b'N');
            let code: u8 = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 0, // N → placeholder 0
            };
            packed[i / 4] |= code << (2 * (i % 4));
        }
        let comp = bsc::compress(&packed)?;
        // Write the v4 block header into a small local buffer, then flush.
        let mut hdr_buf = Vec::with_capacity(12);
        bsc::write_block_frame_header(
            &mut hdr_buf,
            comp.len() as u32,
            block_base_count as u32,
            &comp,
        );
        w.write_all(&hdr_buf)?;
        w.write_all(&comp)?;
        bases_remaining -= block_base_count;
    }
    Ok(total_bases)
}

/// Stream the N-bitmap block-stream over the same flat interval walk into `w`.
///
/// The emitted bytes are byte-identical to what `write_bitmap_blocks` would
/// produce for the same bases' N flags, using the same framing.
///
/// `block_bases` is rounded up to the nearest multiple of 8 internally (bitmap
/// requires byte-aligned blocks: 8 bits/byte). The caller may pass any positive
/// multiple of 4; the rounding ensures blocks start/end on byte boundaries.
pub fn stream_nbitmap(
    w: &mut impl std::io::Write,
    refs: &[strobealign::io::fasta::RefSequence],
    imap: &IntervalMap,
    block_bases: usize,
) -> anyhow::Result<()> {
    assert!(block_bases > 0, "block_bases must be > 0");
    // FlatRefIter indexes refs[iv.ref_id]; reject an out-of-range ref_id up front
    // (clean error instead of an index panic on a misused interval set).
    if !imap
        .intervals()
        .iter()
        .all(|iv| (iv.ref_id as usize) < refs.len())
    {
        bail!("interval ref_id out of range for reference set");
    }
    // Round up to the nearest multiple of 8 (bitmap must be byte-aligned).
    let block_bases = block_bases.next_multiple_of(8);
    let total_bases = imap.total_bases();
    let num_blocks = total_bases.div_ceil(block_bases as u64) as u32;
    w.write_all(&num_blocks.to_le_bytes())?;

    let mut iter = FlatRefIter::new(refs, imap);
    let mut bases_remaining = total_bases as usize;
    while bases_remaining > 0 {
        let block_bit_count = block_bases.min(bases_remaining);
        let bitmap_bytes_count = block_bit_count.div_ceil(8);
        let mut bitmap = vec![0u8; bitmap_bytes_count];
        for i in 0..block_bit_count {
            let b = iter.next().unwrap_or(b'N');
            if !matches!(b, b'A' | b'C' | b'G' | b'T') {
                bitmap[i / 8] |= 1u8 << (i % 8);
            }
        }
        let comp = bsc::compress(&bitmap)?;
        let mut hdr_buf = Vec::with_capacity(12);
        bsc::write_block_frame_header(
            &mut hdr_buf,
            comp.len() as u32,
            block_bit_count as u32,
            &comp,
        );
        w.write_all(&hdr_buf)?;
        w.write_all(&comp)?;
        bases_remaining -= block_bit_count;
    }
    Ok(())
}

/// Pack one backing block's `count` bases (from a positioned `FlatRefIter`) into a
/// framed `[12B v4 header][BSC payload]` 2-bit block — the exact per-block bytes
/// `stream_packed_backing` emits for the same base run. Shared by the parallel builder.
fn pack_backing_block(iter: &mut FlatRefIter, count: usize) -> anyhow::Result<Vec<u8>> {
    let mut packed = vec![0u8; count.div_ceil(4)];
    for i in 0..count {
        let code: u8 = match iter.next().unwrap_or(b'N') {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 0, // N → placeholder 0
        };
        packed[i / 4] |= code << (2 * (i % 4));
    }
    let comp = bsc::compress(&packed)?;
    let mut frame = Vec::with_capacity(12 + comp.len());
    bsc::write_block_frame_header(&mut frame, comp.len() as u32, count as u32, &comp);
    frame.extend_from_slice(&comp);
    Ok(frame)
}

/// Pack one N-bitmap block's `count` bases (from a positioned `FlatRefIter`) into a
/// framed `[12B v4 header][BSC payload]` 1-bit block — the exact per-block bytes
/// `stream_nbitmap` emits for the same base run.
fn pack_nbitmap_block(iter: &mut FlatRefIter, count: usize) -> anyhow::Result<Vec<u8>> {
    let mut bitmap = vec![0u8; count.div_ceil(8)];
    for i in 0..count {
        if !matches!(iter.next().unwrap_or(b'N'), b'A' | b'C' | b'G' | b'T') {
            bitmap[i / 8] |= 1u8 << (i % 8);
        }
    }
    let comp = bsc::compress(&bitmap)?;
    let mut frame = Vec::with_capacity(12 + comp.len());
    bsc::write_block_frame_header(&mut frame, comp.len() as u32, count as u32, &comp);
    frame.extend_from_slice(&comp);
    Ok(frame)
}

/// Per-block `(start_base, count)` partition of `[0, total_bases)` at `block_bases`.
fn block_specs(total_bases: u64, block_bases: usize) -> Vec<(u64, usize)> {
    let num_blocks = total_bases.div_ceil(block_bases as u64);
    (0..num_blocks)
        .map(|b| {
            let start = b * block_bases as u64;
            let count = (block_bases as u64).min(total_bases - start) as usize;
            (start, count)
        })
        .collect()
}

/// Build the full PackedBacking block-stream blob (`[num_blocks u32 LE]` then framed
/// blocks) for `imap` over `refs`, compressing the independent blocks in PARALLEL.
///
/// Byte-equivalent to `stream_packed_backing` for the SAME `block_bases` (each block
/// is the same positioned base run), but materializes the whole blob (the merge can
/// afford it) and uses rayon — so the chr20-sized backing re-derivation in the
/// reference merge runs across all idle cores instead of single-threaded. Returns
/// `(blob, total_bases)`.
pub(crate) fn build_packed_backing_parallel(
    refs: &[strobealign::io::fasta::RefSequence],
    imap: &IntervalMap,
    block_bases: usize,
) -> anyhow::Result<(Vec<u8>, u64)> {
    use rayon::prelude::*;
    assert!(block_bases > 0 && block_bases.is_multiple_of(4));
    if !imap.intervals().iter().all(|iv| (iv.ref_id as usize) < refs.len()) {
        bail!("interval ref_id out of range for reference set");
    }
    let total_bases = imap.total_bases();
    let specs = block_specs(total_bases, block_bases);
    let frames: Vec<Vec<u8>> = specs
        .par_iter()
        .map(|&(start, count)| {
            let mut iter = FlatRefIter::new_at(refs, imap, start);
            pack_backing_block(&mut iter, count)
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    Ok((assemble_blockstream(specs.len() as u32, frames), total_bases))
}

/// Build the full N-bitmap block-stream blob for `imap` over `refs`, compressing the
/// independent blocks in PARALLEL. Byte-equivalent to `stream_nbitmap` for the same
/// (multiple-of-8) `block_bases`.
pub(crate) fn build_nbitmap_parallel(
    refs: &[strobealign::io::fasta::RefSequence],
    imap: &IntervalMap,
    block_bases: usize,
) -> anyhow::Result<Vec<u8>> {
    use rayon::prelude::*;
    assert!(block_bases > 0 && block_bases.is_multiple_of(8));
    if !imap.intervals().iter().all(|iv| (iv.ref_id as usize) < refs.len()) {
        bail!("interval ref_id out of range for reference set");
    }
    let total_bases = imap.total_bases();
    let specs = block_specs(total_bases, block_bases);
    let frames: Vec<Vec<u8>> = specs
        .par_iter()
        .map(|&(start, count)| {
            let mut iter = FlatRefIter::new_at(refs, imap, start);
            pack_nbitmap_block(&mut iter, count)
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    Ok(assemble_blockstream(specs.len() as u32, frames))
}

/// Concatenate `[num_blocks u32 LE]` + the framed blocks in order.
fn assemble_blockstream(num_blocks: u32, frames: Vec<Vec<u8>>) -> Vec<u8> {
    let total: usize = 4 + frames.iter().map(|f| f.len()).sum::<usize>();
    let mut out = Vec::with_capacity(total);
    out.extend_from_slice(&num_blocks.to_le_bytes());
    for f in frames {
        out.extend_from_slice(&f);
    }
    out
}

// ---------------------------------------------------------------------------
// Core-plus-halo windowed consensus builder (Task 7).
//
// Compress can't hold all reads in RAM, so each placed read's reference-forward
// bases are bucketed into fixed-size genome "cores" (windows of `window` bases)
// keyed by `(ref_id, core_idx)`. A read straddling a core boundary contributes
// its overlapping slice to EACH core it touches ("core-plus-halo"), but every
// genome POSITION belongs to exactly ONE core, so the final consensus + interval
// map carry NO duplicate positions.
//
// Coordinate rule (the crux): a read at `[ref_pos, ref_pos+len)` contributing to
// core `[core_start, core_end)` stores the contribution's ACTUAL reference start
// `max(ref_pos, core_start)` and only the bases sub-slice lying in
// `[core_start, core_end)` — NOT the read's `ref_pos`. A boundary read's
// second-core slice begins at the core boundary.
//
// Two builders produce BYTE-IDENTICAL `(PackedBacking, IntervalMap)`:
//   * `ConsensusBuilder` — in-memory (unit tests + small inputs).
//   * `CoreBucketSet` + `finalize_buckets` — disk-backed streaming (Task 12).
// They share `finalize_cores` so the boundary/run-splitting logic can't diverge.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Generic record-boundary block-stream helper (Task 5).
//
// Serializes an arbitrary iterator of byte records as a sequence of
// BSC-compressed v4 blocks (via `paired::streams::write_block_stream`), where
// each block ends on a record boundary and carries a per-block record count in
// the v4 header.  The reader reverses the process: v4-frame parse →
// BSC-decompress → varint-length-prefixed record parse.
//
// Record framing inside each raw (pre-BSC) block:
//   [len: varint][bytes: len bytes]  repeated record_count times.
// ---------------------------------------------------------------------------

use super::super::paired::streams as block_streams;

/// Write `records` as a BSC-compressed, record-boundary v4 block-stream into `out`.
///
/// Records accumulate into a raw buffer; when adding a record would push the
/// buffer past `target_bytes`, the current buffer is closed and BSC-compressed
/// into a new block **before** adding the record.  A single record larger than
/// `target_bytes` still goes in its own block (never split across blocks).
///
/// Returns the total number of records written.
pub fn write_record_block_stream<'a>(
    out: &mut Vec<u8>,
    records: impl Iterator<Item = &'a [u8]>,
    target_bytes: usize,
) -> anyhow::Result<u64> {
    use rayon::prelude::*;

    // Phase 1 (serial, cheap): partition records into raw block buffers at
    // `target_bytes` boundaries (never splitting a record across blocks).
    let mut raw_blocks: Vec<(u32, Vec<u8>)> = Vec::new();
    let mut raw: Vec<u8> = Vec::new();
    let mut block_record_count: u32 = 0;
    let mut total: u64 = 0;

    for rec in records {
        // Encoded size of this record: varint(len) + payload.
        let mut len_prefix = Vec::new();
        write_varint(&mut len_prefix, rec.len() as u64);
        let encoded_len = len_prefix.len() + rec.len();

        // Close the current block if it's non-empty AND adding this record
        // would exceed the target.  (An empty block always accepts the record,
        // regardless of size — never split a record across blocks.)
        if !raw.is_empty() && raw.len() + encoded_len > target_bytes {
            raw_blocks.push((block_record_count, std::mem::take(&mut raw)));
            block_record_count = 0;
        }

        // Append [len varint][bytes] to the current raw buffer.
        raw.extend_from_slice(&len_prefix);
        raw.extend_from_slice(rec);
        block_record_count += 1;
        total += 1;
    }

    // Flush the final non-empty block.
    if !raw.is_empty() {
        raw_blocks.push((block_record_count, raw));
    }

    // Phase 2 (parallel): BSC-compress each block concurrently. This runs at
    // finalize, where the per-chunk pipeline is done and the cores are otherwise
    // idle — the fallback pool's serial BSC was the long pole there (~10s on
    // chr20). Byte-identical: same block boundaries, same deterministic
    // `bsc::compress` per block, written in the same order.
    let blocks: Vec<(u32, Vec<u8>)> = raw_blocks
        .into_par_iter()
        .map(|(rc, raw)| bsc::compress(&raw).map(|c| (rc, c)))
        .collect::<anyhow::Result<Vec<_>>>()?;

    block_streams::write_block_stream(out, &blocks);
    Ok(total)
}

/// Write an `IntervalMap` as a record-boundary block-stream.
///
/// One record per interval: `[ref_id:varint][ref_start:varint][length:varint]`.
/// Returns the interval count (= record count written).
pub fn write_intervalmap_blockstream(
    out: &mut Vec<u8>,
    imap: &IntervalMap,
    target_bytes: usize,
) -> anyhow::Result<u64> {
    let records: Vec<Vec<u8>> = imap
        .intervals()
        .iter()
        .map(|iv| {
            let mut rec = Vec::new();
            write_varint(&mut rec, iv.ref_id as u64);
            write_varint(&mut rec, iv.ref_start);
            write_varint(&mut rec, iv.length);
            rec
        })
        .collect();
    write_record_block_stream(out, records.iter().map(|v| v.as_slice()), target_bytes)
}

/// Read an `IntervalMap` block-stream written by `write_intervalmap_blockstream`.
///
/// Parses each record as `[ref_id:varint][ref_start:varint][length:varint]`,
/// rejecting records with trailing bytes, then validates via
/// `IntervalMap::from_intervals` (ordering + no-overlap).
pub fn read_intervalmap_blockstream(data: &[u8]) -> anyhow::Result<IntervalMap> {
    let raw_records = read_record_block_stream(data)?;
    let mut intervals = Vec::with_capacity(raw_records.len());
    for (i, rec) in raw_records.iter().enumerate() {
        let mut off = 0usize;
        let ref_id_u64 = read_varint(rec, &mut off).ok_or_else(|| {
            anyhow::anyhow!("intervalmap blockstream record {i}: truncated ref_id")
        })?;
        let ref_id: u32 = ref_id_u64.try_into().map_err(|_| {
            anyhow::anyhow!("intervalmap blockstream record {i}: ref_id {ref_id_u64} exceeds u32")
        })?;
        let ref_start = read_varint(rec, &mut off).ok_or_else(|| {
            anyhow::anyhow!("intervalmap blockstream record {i}: truncated ref_start")
        })?;
        let length = read_varint(rec, &mut off).ok_or_else(|| {
            anyhow::anyhow!("intervalmap blockstream record {i}: truncated length")
        })?;
        if off != rec.len() {
            anyhow::bail!(
                "intervalmap blockstream record {i}: {} trailing bytes",
                rec.len() - off
            );
        }
        intervals.push(Interval {
            ref_id,
            ref_start,
            length,
        });
    }
    IntervalMap::from_intervals(intervals)
}

/// Read a record-boundary v4 block-stream written by `write_record_block_stream`.
///
/// Parses the v4 frame via `decode_block_payloads`, BSC-decompresses each
/// block, then extracts varint-length-prefixed records.  Validates:
/// - each record's `[len]` stays in-bounds,
/// - after parsing a block you consumed exactly its bytes,
/// - the number of records parsed equals the block's declared `record_count`.
pub fn read_record_block_stream(data: &[u8]) -> anyhow::Result<Vec<Vec<u8>>> {
    let block_payloads = block_streams::decode_block_payloads(data)?;
    let mut all_records: Vec<Vec<u8>> = Vec::new();

    for (block_idx, (record_count, compressed)) in block_payloads.into_iter().enumerate() {
        let raw = bsc::decompress(&compressed)?;
        let mut off = 0usize;
        let mut parsed: u32 = 0;

        while off < raw.len() {
            let len = read_varint(&raw, &mut off).ok_or_else(|| {
                anyhow::anyhow!("record block {block_idx}: truncated length varint")
            })?;
            let len = len as usize;
            // checked_add: `len` is an untrusted varint from attacker-controlled
            // BSC-decompressed bytes; a plain `off + len` would wrap in release (no
            // overflow-checks) and bypass the `> len()` guard, panicking on the slice.
            let end = off
                .checked_add(len)
                .filter(|&e| e <= raw.len())
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "record block {block_idx}: record length {len} out of bounds (remaining {})",
                        raw.len() - off
                    )
                })?;
            all_records.push(raw[off..end].to_vec());
            off = end;
            parsed += 1;
        }

        if off != raw.len() {
            anyhow::bail!(
                "record block {block_idx}: {} trailing bytes after {} records",
                raw.len() - off,
                parsed
            );
        }
        if parsed != record_count {
            anyhow::bail!(
                "record block {block_idx}: parsed {parsed} records but header declared {record_count}"
            );
        }
    }

    Ok(all_records)
}

#[cfg(test)]
mod packstore_tests {
    use super::*;

    /// Build a `PackedBacking` directly from ASCII ACGTN bytes.
    ///
    /// Bit layout matches the spec (§5.2): 2-bit packed A=0,C=1,G=2,T=3 LSB-first,
    /// 4 bases/byte; N-bitmap bit `i` set iff base `i` is `N`.
    fn consensus_from(bases: &[u8]) -> PackedBacking {
        let len = bases.len();
        let mut packed = vec![0u8; len.div_ceil(4)];
        let mut n_bitmap = vec![0u8; len.div_ceil(8)];
        for (i, &b) in bases.iter().enumerate() {
            let code: u8 = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => {
                    // N (or any unrecognized byte) → set bitmap bit, leave 2-bit as 0.
                    n_bitmap[i / 8] |= 1u8 << (i % 8);
                    0
                }
            };
            if b != b'N' {
                packed[i / 4] |= code << (2 * (i % 4));
            }
        }
        PackedBacking::from_parts(packed, n_bitmap, len)
    }

    #[test]
    fn packed_consensus_block_roundtrip() {
        let bases: Vec<u8> = (0..600).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        let c = consensus_from(&bases);
        let mut buf = Vec::new();
        write_packed_backing_blocks(&mut buf, &c, /*block_bases*/ 256).unwrap();
        let decoded = read_packed_backing_blocks(&buf, c.len()).unwrap();
        // ACGT-only data has no N positions, so packed-only decode equals full.
        assert_eq!(decoded.unpack(), c.unpack());
    }

    #[test]
    fn n_bitmap_framing_roundtrip() {
        // Consensus with several N positions.
        let bases: Vec<u8> = (0..600)
            .map(|i| {
                if i % 37 == 0 {
                    b'N'
                } else {
                    b"ACGT"[(i * 7 % 4) as usize]
                }
            })
            .collect();
        let c = consensus_from(&bases);
        // Sanity: there ARE Ns.
        assert!((0..c.len()).any(|i| c.is_n(i)));

        let mut buf = Vec::new();
        write_bitmap_blocks(
            &mut buf,
            c.n_bitmap_bytes(),
            c.len(),
            /*block_bits*/ 256,
        )
        .unwrap();
        let decoded_bitmap = read_bitmap_blocks(&buf, c.len()).unwrap();
        assert_eq!(decoded_bitmap.len(), c.len().div_ceil(8));
        // Decoded bitmap must reproduce is_n for every position.
        let overlaid =
            PackedBacking::from_parts(c.packed_bytes().to_vec(), decoded_bitmap, c.len());
        for i in 0..c.len() {
            assert_eq!(overlaid.is_n(i), c.is_n(i), "is_n mismatch at {i}");
        }
        // And the full unpack matches with the overlay.
        assert_eq!(overlaid.unpack(), c.unpack());
    }

    #[test]
    fn multi_block_boundary_counts() {
        let bases: Vec<u8> = (0..600).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        let c = consensus_from(&bases);
        let mut buf = Vec::new();
        write_packed_backing_blocks(&mut buf, &c, 256).unwrap();
        // num_blocks prefix == 3 (256, 256, 88).
        let num_blocks = u32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]);
        assert_eq!(num_blocks, 3);
        // Decode still produces all 600 bases.
        let decoded = read_packed_backing_blocks(&buf, 600).unwrap();
        assert_eq!(decoded.len(), 600);
        assert_eq!(decoded.unpack(), c.unpack());
    }

    #[test]
    fn block_bases_must_be_multiple_of_four() {
        let c = consensus_from(b"ACGTACGT");
        let mut buf = Vec::new();
        // A valid multiple-of-4 block_bases works.
        write_packed_backing_blocks(&mut buf, &c, 4).unwrap();
        let decoded = read_packed_backing_blocks(&buf, c.len()).unwrap();
        assert_eq!(decoded.unpack(), c.unpack());

        // A non-multiple-of-4 block_bases must panic (assert).
        let mut buf2 = Vec::new();
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = write_packed_backing_blocks(&mut buf2, &c, 6);
        }));
        assert!(r.is_err(), "non-multiple-of-4 block_bases should panic");
    }

    #[test]
    fn inflated_num_blocks_rejected_up_front() {
        // A valid stream whose num_blocks prefix is overwritten with u32::MAX must
        // be rejected BEFORE the per-block loop (the Task 14 upfront DoS bound),
        // not after ~4B iterations.
        let bases: Vec<u8> = (0..600).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        let c = consensus_from(&bases);
        let mut buf = Vec::new();
        write_packed_backing_blocks(&mut buf, &c, 256).unwrap();
        buf[0..4].copy_from_slice(&u32::MAX.to_le_bytes());
        let e = match read_packed_backing_blocks(&buf, 600) {
            Ok(_) => panic!("inflated num_blocks must be rejected"),
            Err(e) => e,
        };
        assert!(e.to_string().contains("exceeds"), "got: {e}");
        // Same guard in the fast-verify CRC walker.
        let e2 = verify_block_stream_crc(&buf).unwrap_err();
        assert!(e2.to_string().contains("exceeds"), "got: {e2}");

        // And the bitmap stream.
        let mut bbuf = Vec::new();
        write_bitmap_blocks(&mut bbuf, c.n_bitmap_bytes(), c.len(), 256).unwrap();
        bbuf[0..4].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(read_bitmap_blocks(&bbuf, 600).is_err());
    }

    #[test]
    fn hostile_block_record_count_rejected_before_allocation() {
        // A single 12-byte block header (empty payload) claiming u32::MAX records,
        // with `total_records` set to match. The header-only pre-pass must reject it
        // via the single-block decode cap BEFORE the `out_len`-sized allocation —
        // u32::MAX/4 ≈ 1.07e9 bytes > the 64 MiB per-block decode cap. Reaching this
        // assert proves the process survived (no resize-driven OOM abort).
        let mut buf = Vec::new();
        buf.extend_from_slice(&1u32.to_le_bytes()); // num_blocks = 1
        bsc::write_block_frame_header(
            &mut buf,
            /*block_len*/ 0,
            /*record_count*/ u32::MAX,
            &[],
        );
        let err = match read_packed_backing_blocks(&buf, u32::MAX as usize) {
            Ok(_) => panic!("hostile block record_count must be rejected"),
            Err(e) => e,
        };
        assert!(
            err.to_string().contains("decode cap"),
            "expected single-block decode-cap rejection, got: {err}"
        );

        // A footer that declares more records than the (valid) blocks provide is
        // also rejected by the pre-pass Σ check, before allocation: build a real
        // 600-base stream but ask for 600 + a huge surplus.
        let bases: Vec<u8> = (0..600).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        let c = consensus_from(&bases);
        let mut good = Vec::new();
        write_packed_backing_blocks(&mut good, &c, 256).unwrap();
        let err2 = match read_packed_backing_blocks(&good, 600 + 1_000_000_000) {
            Ok(_) => panic!("inflated total_records must be rejected"),
            Err(e) => e,
        };
        assert!(
            err2.to_string().contains("Σ record counts")
                || err2.to_string().contains("exceeds total"),
            "expected Σ-mismatch rejection, got: {err2}"
        );
    }

    #[test]
    fn corrupt_block_crc_is_rejected() {
        let bases: Vec<u8> = (0..600).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        let c = consensus_from(&bases);
        let mut buf = Vec::new();
        write_packed_backing_blocks(&mut buf, &c, 256).unwrap();
        // Flip a payload byte (past the 4-byte num_blocks prefix + first 12-byte
        // v4 header → into the first block payload).
        let idx = 4 + 12 + 1;
        buf[idx] ^= 0xFF;
        assert!(read_packed_backing_blocks(&buf, 600).is_err());
    }
}

#[cfg(test)]
mod record_block_stream_tests {
    use super::*;

    #[test]
    fn record_block_stream_roundtrips_and_counts_sum() {
        let recs: Vec<Vec<u8>> = (0..5u8).map(|i| vec![i; 10]).collect();
        let mut out = Vec::new();
        let total =
            write_record_block_stream(&mut out, recs.iter().map(|v| v.as_slice()), 16).unwrap();
        assert_eq!(total, 5);
        let got = read_record_block_stream(&out).unwrap();
        assert_eq!(got, recs);
    }

    #[test]
    fn record_block_stream_single_oversized_record_ok() {
        let recs: Vec<Vec<u8>> = vec![vec![7u8; 100]]; // bigger than target -> own block
        let mut out = Vec::new();
        let total =
            write_record_block_stream(&mut out, recs.iter().map(|v| v.as_slice()), 16).unwrap();
        assert_eq!(total, 1);
        assert_eq!(read_record_block_stream(&out).unwrap(), recs);
    }

    /// Hostile input: a block whose first record's length varint decodes to
    /// `u64::MAX`. A plain `off + len` add would wrap in release (no overflow-checks)
    /// and bypass the bounds check, panicking on the slice. The `checked_add` guard
    /// must reject it with a clean error instead of aborting.
    #[test]
    fn record_block_stream_rejects_overflow_length_varint() {
        use crate::compression::paired::streams::write_block_stream;
        // 10-byte LEB128 0xFF×9,0x01 decodes to u64::MAX, then 2 payload bytes.
        let mut raw = vec![0xFFu8; 9];
        raw.push(0x01);
        raw.extend_from_slice(&[0xAA, 0xBB]);
        let comp = bsc::compress(&raw).unwrap();
        let mut data = Vec::new();
        write_block_stream(&mut data, &[(1u32, comp)]);
        let err = read_record_block_stream(&data).unwrap_err();
        assert!(
            err.to_string().contains("out of bounds"),
            "expected clean rejection, got: {err}"
        );
    }
}

#[cfg(test)]
mod interval_tests {
    use super::*;

    #[test]
    fn interval_map_translates_and_validates() {
        let im = IntervalMap::from_intervals(vec![
            Interval {
                ref_id: 0,
                ref_start: 100,
                length: 10,
            },
            Interval {
                ref_id: 1,
                ref_start: 50,
                length: 10,
            },
        ])
        .unwrap();
        assert_eq!(im.total_bases(), 20);
        assert_eq!(im.ref_to_base(0, 105), Some(5));
        assert_eq!(im.ref_to_base(1, 55), Some(15));
        assert_eq!(im.ref_to_base(0, 120), None); // outside interval
        assert!(im.span_in_single_interval(5, 4)); // base 5..9 in interval 0
        assert!(!im.span_in_single_interval(8, 5)); // 8..13 crosses boundary
    }

    #[test]
    fn interval_map_rejects_overlap_and_disorder() {
        assert!(
            IntervalMap::from_intervals(vec![
                Interval {
                    ref_id: 0,
                    ref_start: 100,
                    length: 20,
                },
                Interval {
                    ref_id: 0,
                    ref_start: 110,
                    length: 5,
                },
            ])
            .is_err()
        );
        assert!(
            IntervalMap::from_intervals(vec![
                Interval {
                    ref_id: 1,
                    ref_start: 0,
                    length: 5,
                },
                Interval {
                    ref_id: 0,
                    ref_start: 0,
                    length: 5,
                },
            ])
            .is_err()
        );
    }

    #[test]
    fn interval_map_serialize_roundtrips() {
        let im = IntervalMap::from_intervals(vec![
            Interval {
                ref_id: 0,
                ref_start: 100,
                length: 10,
            },
            Interval {
                ref_id: 1,
                ref_start: 50,
                length: 10,
            },
            Interval {
                ref_id: 1,
                ref_start: 1_000_000,
                length: 250,
            },
        ])
        .unwrap();
        let bytes = im.serialize();
        let back = IntervalMap::deserialize(&bytes).unwrap();
        assert_eq!(back.intervals, im.intervals);
        assert_eq!(back.total_bases(), im.total_bases());
        assert_eq!(back.serialize(), bytes);
    }

    #[test]
    fn deserialize_revalidates_overlap() {
        // Hand-craft a serialized form with two overlapping same-ref intervals;
        // deserialize must re-run from_intervals validation and reject it.
        let mut bytes = Vec::new();
        write_varint(&mut bytes, 2); // count
        write_varint(&mut bytes, 0); // ref_id
        write_varint(&mut bytes, 100); // ref_start
        write_varint(&mut bytes, 20); // length -> covers [100,120)
        write_varint(&mut bytes, 0); // ref_id
        write_varint(&mut bytes, 110); // ref_start (overlaps previous)
        write_varint(&mut bytes, 5); // length
        assert!(IntervalMap::deserialize(&bytes).is_err());
    }

    #[test]
    fn deserialize_rejects_trailing_and_oversized_count() {
        // Trailing byte after a valid single interval.
        let im = IntervalMap::from_intervals(vec![Interval {
            ref_id: 0,
            ref_start: 0,
            length: 5,
        }])
        .unwrap();
        let mut bytes = im.serialize();
        bytes.push(0xAB);
        assert!(IntervalMap::deserialize(&bytes).is_err());

        // Oversized count vs available bytes must not over-allocate / must error.
        let mut hostile = Vec::new();
        write_varint(&mut hostile, u64::MAX); // absurd count
        assert!(IntervalMap::deserialize(&hostile).is_err());
    }

    #[test]
    fn span_len_zero_inside_interval() {
        let im = IntervalMap::from_intervals(vec![Interval {
            ref_id: 0,
            ref_start: 0,
            length: 10,
        }])
        .unwrap();
        assert!(im.span_in_single_interval(5, 0)); // empty span inside
        assert!(!im.span_in_single_interval(10, 0)); // at total == not inside
    }

    #[test]
    fn intervalmap_blockstream_roundtrip() {
        let ivs = vec![
            Interval {
                ref_id: 0,
                ref_start: 0,
                length: 100,
            },
            Interval {
                ref_id: 0,
                ref_start: 200,
                length: 50,
            },
            Interval {
                ref_id: 1,
                ref_start: 0,
                length: 10,
            },
        ];
        let imap = IntervalMap::from_intervals(ivs.clone()).unwrap();
        let mut out = Vec::new();
        let n = write_intervalmap_blockstream(&mut out, &imap, 8).unwrap(); // tiny target -> multi-block
        assert_eq!(n, 3);
        let back = read_intervalmap_blockstream(&out).unwrap();
        assert_eq!(back.intervals(), imap.intervals());
    }
}

#[cfg(test)]
mod streaming_backing_tests {
    use super::*;

    /// Extract the normalized (ACGTN, non-ACGTN→N) flat reference bytes over
    /// `imap`'s intervals — the same bytes the streaming walk emits and the same
    /// bytes the materialize path would pack.
    fn flat_over_intervals(
        refs: &[strobealign::io::fasta::RefSequence],
        imap: &IntervalMap,
    ) -> Vec<u8> {
        let mut out = Vec::new();
        for iv in imap.intervals() {
            let seq = &refs[iv.ref_id as usize].sequence;
            for off in 0..iv.length {
                let pos = (iv.ref_start + off) as usize;
                let raw = if pos < seq.len() { seq[pos] } else { b'N' };
                out.push(normalize_ref_byte(raw));
            }
        }
        out
    }

    /// Build the materialized packed/bitmap block bytes for the same flat bases
    /// via `write_packed_backing_blocks` / `write_bitmap_blocks` using the same
    /// bit layout (A=0,C=1,G=2,T=3, 4 bases/byte LSB-first; N-bitmap 1 bit/base).
    fn materialize_blocks(flat: &[u8], block_bases: usize) -> (Vec<u8>, Vec<u8>) {
        let len = flat.len();
        let mut packed_bytes = vec![0u8; len.div_ceil(4)];
        let mut bitmap_bytes = vec![0u8; len.div_ceil(8)];
        for (i, &b) in flat.iter().enumerate() {
            let code: u8 = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => {
                    bitmap_bytes[i / 8] |= 1u8 << (i % 8);
                    0
                }
            };
            if b != b'N' {
                packed_bytes[i / 4] |= code << (2 * (i % 4));
            }
        }
        let pb = PackedBacking::from_parts(packed_bytes, bitmap_bytes, len);
        let mut packed = Vec::new();
        write_packed_backing_blocks(&mut packed, &pb, block_bases).unwrap();
        let mut nbm = Vec::new();
        // The bitmap writer rounds nothing; block_bases must be a multiple of 8,
        // matching the streaming writer's internal `next_multiple_of(8)`.
        write_bitmap_blocks(&mut nbm, pb.n_bitmap_bytes(), pb.len(), block_bases).unwrap();
        (packed, nbm)
    }

    /// (a) Byte-identity vs the materializing writers over a 3-interval input
    /// whose total bases are NOT a multiple of block_bases and with a block
    /// boundary falling mid-interval.
    #[test]
    fn stream_backing_byte_identical_to_materialize() {
        // contig 0, len 40.
        let refbytes: Vec<u8> = b"ACGTNNACGTACGTACGTACATCGGATCNTAGCATGCATGC".to_vec();
        let refs = vec![strobealign::io::fasta::RefSequence {
            name: "c".into(),
            sequence: refbytes.clone(),
        }];
        // Three intervals: [0,7), [9,22), [25,38). Total = 7+13+13 = 33 bases.
        // With block_bases = 8: 33 is not a multiple of 8 (final block = 1 base),
        // and the block boundary at flat index 8 lands mid-interval-2 (the second
        // interval starts at flat index 7).
        let imap = IntervalMap::from_intervals(vec![
            Interval {
                ref_id: 0,
                ref_start: 0,
                length: 7,
            },
            Interval {
                ref_id: 0,
                ref_start: 9,
                length: 13,
            },
            Interval {
                ref_id: 0,
                ref_start: 25,
                length: 13,
            },
        ])
        .unwrap();
        let block_bases = 8usize;

        let flat = flat_over_intervals(&refs, &imap);
        assert_eq!(flat.len(), 33);
        assert_ne!(
            flat.len() % block_bases,
            0,
            "total must not be block multiple"
        );

        let mut packed = Vec::new();
        let bases = stream_packed_backing(&mut packed, &refs, &imap, block_bases).unwrap();
        let mut nbm = Vec::new();
        stream_nbitmap(&mut nbm, &refs, &imap, block_bases).unwrap();
        assert_eq!(bases, 33);

        let (exp_packed, exp_nbm) = materialize_blocks(&flat, block_bases);
        assert_eq!(
            packed, exp_packed,
            "streamed packed bytes must match materialize"
        );
        assert_eq!(
            nbm, exp_nbm,
            "streamed n-bitmap bytes must match materialize"
        );

        // And the decoded consensus equals the expected flat ACGTN.
        let pb = read_packed_backing_with_bitmap(&packed, &nbm, bases as usize).unwrap();
        assert_eq!(pb.unpack(), flat);
    }

    /// The PARALLEL merge builders must be byte-identical to the SERIAL streamers for
    /// the same `block_bases` — i.e. `FlatRefIter::new_at` seeking + per-block packing
    /// reproduce the same blocks the serial flat walk emits. Multi-interval with a
    /// block boundary landing mid-interval (the case where seeking is non-trivial).
    #[test]
    fn parallel_backing_byte_identical_to_serial() {
        let refbytes: Vec<u8> = b"ACGTNNACGTACGTACGTACATCGGATCNTAGCATGCATGC".to_vec();
        let refs = vec![strobealign::io::fasta::RefSequence {
            name: "c".into(),
            sequence: refbytes,
        }];
        let imap = IntervalMap::from_intervals(vec![
            Interval { ref_id: 0, ref_start: 0, length: 7 },
            Interval { ref_id: 0, ref_start: 9, length: 13 },
            Interval { ref_id: 0, ref_start: 25, length: 13 },
        ])
        .unwrap();
        // block_bases = 8 (multiple of 8) → 33 bases split into 5 blocks, boundaries
        // landing mid-interval; exercises new_at across interval edges.
        let block_bases = 8usize;

        let mut serial_packed = Vec::new();
        let serial_bases =
            stream_packed_backing(&mut serial_packed, &refs, &imap, block_bases).unwrap();
        let mut serial_nbm = Vec::new();
        stream_nbitmap(&mut serial_nbm, &refs, &imap, block_bases).unwrap();

        let (par_packed, par_bases) =
            build_packed_backing_parallel(&refs, &imap, block_bases).unwrap();
        let par_nbm = build_nbitmap_parallel(&refs, &imap, block_bases).unwrap();

        assert_eq!(par_bases, serial_bases, "parallel base count must match serial");
        assert_eq!(par_packed, serial_packed, "parallel packed blob must match serial");
        assert_eq!(par_nbm, serial_nbm, "parallel n-bitmap blob must match serial");
    }

    /// (b) An N at flat indices adjacent to a block boundary (block_bases-1 and
    /// block_bases) survives the roundtrip.
    #[test]
    fn stream_backing_n_adjacent_to_block_boundary() {
        let block_bases = 8usize;
        // contig 0, len 32: put N at ref positions 7 and 8 (which become flat
        // indices 7 = block_bases-1 and 8 = block_bases for a single interval
        // starting at ref 0).
        let mut refbytes: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        refbytes[7] = b'N';
        refbytes[8] = b'N';
        let refs = vec![strobealign::io::fasta::RefSequence {
            name: "c".into(),
            sequence: refbytes.clone(),
        }];
        // Single interval [0,20) → flat indices 0..20, so flat 7 and 8 are the Ns.
        let imap = IntervalMap::from_intervals(vec![Interval {
            ref_id: 0,
            ref_start: 0,
            length: 20,
        }])
        .unwrap();

        let mut packed = Vec::new();
        let bases = stream_packed_backing(&mut packed, &refs, &imap, block_bases).unwrap();
        let mut nbm = Vec::new();
        stream_nbitmap(&mut nbm, &refs, &imap, block_bases).unwrap();

        let pb = read_packed_backing_with_bitmap(&packed, &nbm, bases as usize).unwrap();
        let flat = flat_over_intervals(&refs, &imap);
        assert_eq!(pb.unpack(), flat);
        assert!(
            pb.is_n(block_bases - 1),
            "flat base block_bases-1 must be N"
        );
        assert!(pb.is_n(block_bases), "flat base block_bases must be N");
        assert!(
            !pb.is_n(block_bases - 2),
            "flat base block_bases-2 must not be N"
        );
        assert!(
            !pb.is_n(block_bases + 1),
            "flat base block_bases+1 must not be N"
        );
    }

    /// (c) The existing 2-interval roundtrip via `read_packed_backing_with_bitmap`.
    #[test]
    fn stream_backing_matches_reader_over_intervals() {
        // contig 0, len 20: "ACGTNNACGTACGTACGTAC"
        let refbytes: Vec<u8> = b"ACGTNNACGTACGTACGTAC".to_vec();
        let refs = vec![strobealign::io::fasta::RefSequence {
            name: "c".into(),
            sequence: refbytes.clone(),
        }];
        // Two intervals: [0,8) and [12,20) — gaps at [8,12) excluded.
        let imap = IntervalMap::from_intervals(vec![
            Interval {
                ref_id: 0,
                ref_start: 0,
                length: 8,
            },
            Interval {
                ref_id: 0,
                ref_start: 12,
                length: 8,
            },
        ])
        .unwrap();

        // Tiny block_bases=4 forces multi-block (16 total bases → 4 blocks).
        let mut packed = Vec::new();
        let bases = stream_packed_backing(&mut packed, &refs, &imap, 4).unwrap();
        let mut nbm = Vec::new();
        stream_nbitmap(&mut nbm, &refs, &imap, 4).unwrap();

        let pb = read_packed_backing_with_bitmap(&packed, &nbm, bases as usize).unwrap();

        // refbytes = "ACGTNNACGTACGTACGTAC" (indices 0-19)
        //   [0..8]  = "ACGTNNAC"
        //   [12..20] = "GTACGTAC"
        let expect: Vec<u8> = [&refbytes[0..8], &refbytes[12..20]].concat();
        assert_eq!(pb.unpack(), expect);
        assert_eq!(bases, 16);

        // Verify N positions at flat indices 4 and 5 (refbytes[4]='N', refbytes[5]='N').
        assert!(pb.is_n(4), "flat base 4 must be N");
        assert!(pb.is_n(5), "flat base 5 must be N");
        assert!(!pb.is_n(0), "flat base 0 must not be N");
        assert!(!pb.is_n(6), "flat base 6 must not be N");
    }

    /// Fix 1: an interval with an out-of-range ref_id yields a clean error, not a
    /// panic, from both streaming writers.
    #[test]
    fn stream_backing_rejects_out_of_range_ref_id() {
        let refs = vec![strobealign::io::fasta::RefSequence {
            name: "c".into(),
            sequence: b"ACGTACGT".to_vec(),
        }];
        // ref_id 1 is out of range (only refs[0] exists).
        let imap = IntervalMap::from_intervals(vec![Interval {
            ref_id: 1,
            ref_start: 0,
            length: 4,
        }])
        .unwrap();
        let mut packed = Vec::new();
        assert!(stream_packed_backing(&mut packed, &refs, &imap, 8).is_err());
        let mut nbm = Vec::new();
        assert!(stream_nbitmap(&mut nbm, &refs, &imap, 8).is_err());
    }
}
