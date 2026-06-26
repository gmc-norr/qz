//! Per-(chunk,mate) FallbackPool block-stream (Task 11; bounded-RAM refactor).
//!
//! When a read cannot be encoded as edits against the consensus it is emitted
//! verbatim as a *fallback literal*. Each (chunk,mate)'s fallback literals are
//! written inline into that chunk's payload as their own FallbackPool block-stream
//! during encode (no global spool). This module provides:
//!   1. [`write_fallback_pool`] — writes one (chunk,mate)'s FallbackPool as a v4
//!      record-boundary block-stream (delegates to `backing`), rejecting any
//!      literal larger than [`PER_RECORD_MAX_BYTES`].
//!   2. [`FallbackPoolReader`] — a lazy reader that decodes exactly ONE block
//!      at a time and yields its literals in order (bounded by one block).

use super::super::bsc;
use super::super::dna_utils::read_varint;
use super::super::paired::streams as block_streams;
use super::backing;
use anyhow::{Result, anyhow, bail};

/// Hard cap on a single fallback literal. Anything larger is rejected before it
/// is written — a literal this size signals corrupt/hostile input, not a read.
pub const PER_RECORD_MAX_BYTES: usize = 8 * 1024 * 1024;

/// Write a per-(chunk,mate) FallbackPool as a v4 record-boundary block-stream
/// from the mate's fallback literals (in read order). The literals ARE the
/// records, so this delegates to `backing::write_record_block_stream`. Returns
/// the literal count written.
///
/// Enforces the [`PER_RECORD_MAX_BYTES`] per-literal cap up front (this used to
/// be the spool's guard; the spool is gone, so the cap lives here): a literal
/// larger than the cap signals corrupt/hostile input, not a read.
pub fn write_fallback_pool(
    out: &mut Vec<u8>,
    literals: &[&[u8]],
    target_bytes: usize,
) -> Result<u64> {
    for (i, lit) in literals.iter().enumerate() {
        if lit.len() > PER_RECORD_MAX_BYTES {
            bail!(
                "fallback literal {i} of {} bytes exceeds PER_RECORD_MAX_BYTES ({PER_RECORD_MAX_BYTES})",
                lit.len()
            );
        }
    }
    backing::write_record_block_stream(out, literals.iter().copied(), target_bytes)
}

/// Lazy reader over a FallbackPool block-stream: decodes ONE block at a time
/// and yields its literals in order. Never materializes more than a single
/// decoded block at once.
pub struct FallbackPoolReader<'a> {
    /// Per-block headers `(record_count, compressed_bytes)` in stream order.
    blocks: Vec<(u32, Vec<u8>)>,
    /// Index of the next block to decode.
    next_block: usize,
    /// Decoded raw bytes of the current block.
    cur: Vec<u8>,
    /// Cursor into `cur`.
    cur_off: usize,
    /// Literals already parsed out of the current block.
    cur_parsed: u32,
    /// Declared record count of the current block.
    cur_count: u32,
    /// Ties the reader to the borrowed input slice (bounded contract: payloads
    /// are decoded one block at a time, never all up front).
    _marker: std::marker::PhantomData<&'a [u8]>,
}

impl<'a> FallbackPoolReader<'a> {
    /// Parse the v4 frame headers (cheap — no BSC decode yet). Holding the
    /// borrowed `&[u8]` keeps the bounded contract: block payloads are only
    /// BSC-decompressed one at a time in [`Self::next_literal`].
    pub fn new(data: &'a [u8]) -> Result<Self> {
        let blocks = block_streams::decode_block_payloads(data)?;
        Ok(Self {
            blocks,
            next_block: 0,
            cur: Vec::new(),
            cur_off: 0,
            cur_parsed: 0,
            cur_count: 0,
            _marker: std::marker::PhantomData,
        })
    }

    /// Decode and load the next block into `cur`. Returns false when no blocks
    /// remain. Validates the previous block was fully consumed.
    fn load_next_block(&mut self) -> Result<bool> {
        // Before advancing, the previous block must have been fully consumed.
        if self.cur_off != self.cur.len() {
            bail!(
                "fallback pool block {}: {} trailing bytes after {} literals",
                self.next_block.saturating_sub(1),
                self.cur.len() - self.cur_off,
                self.cur_parsed
            );
        }
        if self.cur_parsed != self.cur_count {
            bail!(
                "fallback pool block {}: parsed {} literals but header declared {}",
                self.next_block.saturating_sub(1),
                self.cur_parsed,
                self.cur_count
            );
        }
        if self.next_block >= self.blocks.len() {
            return Ok(false);
        }
        let (record_count, compressed) = &self.blocks[self.next_block];
        self.cur = bsc::decompress(compressed)?;
        self.cur_off = 0;
        self.cur_parsed = 0;
        self.cur_count = *record_count;
        self.next_block += 1;
        Ok(true)
    }

    /// Yield the next literal in stream order, or `None` at end of stream.
    /// Decodes block-by-block — at most one decompressed block is resident.
    pub fn next_literal(&mut self) -> Result<Option<Vec<u8>>> {
        loop {
            if self.cur_off < self.cur.len() {
                let block_idx = self.next_block - 1;
                let len = read_varint(&self.cur, &mut self.cur_off).ok_or_else(|| {
                    anyhow!("fallback pool block {block_idx}: truncated length varint")
                })? as usize;
                // checked_add: `len` is an untrusted varint from BSC-decompressed
                // (attacker-controlled) bytes; a plain `cur_off + len` would wrap in
                // release (no overflow-checks) and bypass the `> len()` guard, panicking
                // on the slice below. Compute the end once, bounded to the block.
                let end = self
                    .cur_off
                    .checked_add(len)
                    .filter(|&e| e <= self.cur.len())
                    .ok_or_else(|| {
                        anyhow!(
                            "fallback pool block {block_idx}: literal length {len} out of bounds (remaining {})",
                            self.cur.len() - self.cur_off
                        )
                    })?;
                let lit = self.cur[self.cur_off..end].to_vec();
                self.cur_off = end;
                self.cur_parsed += 1;
                return Ok(Some(lit));
            }
            // Current block exhausted — try to load the next one.
            if !self.load_next_block()? {
                return Ok(None);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Borrow a `&[Vec<u8>]` as the `&[&[u8]]` write_fallback_pool now takes.
    fn as_slices(v: &[Vec<u8>]) -> Vec<&[u8]> {
        v.iter().map(Vec::as_slice).collect()
    }

    #[test]
    fn pool_blockstream_roundtrips_in_group_order() {
        let lits: Vec<Vec<u8>> = vec![b"AAAA".to_vec(), b"CCCCCC".to_vec(), b"G".to_vec()];
        let mut out = Vec::new();
        let total = write_fallback_pool(&mut out, &as_slices(&lits), 8).unwrap(); // tiny target -> multi-block
        assert_eq!(total, 3);
        let mut rdr = FallbackPoolReader::new(&out).unwrap();
        assert_eq!(rdr.next_literal().unwrap().unwrap(), b"AAAA");
        assert_eq!(rdr.next_literal().unwrap().unwrap(), b"CCCCCC");
        assert_eq!(rdr.next_literal().unwrap().unwrap(), b"G");
        assert!(rdr.next_literal().unwrap().is_none());
    }

    #[test]
    fn oversized_literal_rejected_by_pool_write() {
        // Ported from the deleted FallbackSpool::push guard: an over-cap literal
        // is now rejected by write_fallback_pool before any block is written.
        let big = vec![b'A'; PER_RECORD_MAX_BYTES + 1];
        let mut out = Vec::new();
        assert!(write_fallback_pool(&mut out, &[big.as_slice()], 16).is_err());
        // A pool of in-bounds literals still writes fine.
        out.clear();
        let n = write_fallback_pool(&mut out, &[b"ok".as_slice()], 16).unwrap();
        assert_eq!(n, 1);
    }

    #[test]
    fn pool_roundtrips_multiple_literals() {
        let lits = vec![b"TT".to_vec(), b"GGGG".to_vec(), b"A".to_vec()];
        let mut out = Vec::new();
        let n = write_fallback_pool(&mut out, &as_slices(&lits), 16).unwrap();
        assert_eq!(n, 3);
        let mut rdr = FallbackPoolReader::new(&out).unwrap();
        let mut got = Vec::new();
        while let Some(l) = rdr.next_literal().unwrap() {
            got.push(l);
        }
        assert_eq!(got, vec![b"TT".to_vec(), b"GGGG".to_vec(), b"A".to_vec()]);
    }

    /// Hostile input: a pool block whose first literal's length varint decodes to
    /// `u64::MAX`. A plain `cur_off + len` add would wrap in release and bypass the
    /// bounds check, panicking on the slice. The `checked_add` guard must reject it
    /// with a clean error instead of aborting.
    #[test]
    fn next_literal_rejects_overflow_length_varint() {
        use crate::compression::bsc;
        use crate::compression::paired::streams::write_block_stream;
        // 10-byte LEB128 0xFF×9,0x01 decodes to u64::MAX, then 2 payload bytes.
        let mut raw = vec![0xFFu8; 9];
        raw.push(0x01);
        raw.extend_from_slice(&[0xAA, 0xBB]);
        let comp = bsc::compress(&raw).unwrap();
        let mut data = Vec::new();
        write_block_stream(&mut data, &[(1u32, comp)]);
        let mut rdr = FallbackPoolReader::new(&data).unwrap();
        let err = rdr.next_literal().unwrap_err();
        assert!(
            err.to_string().contains("out of bounds"),
            "expected clean rejection, got: {err}"
        );
    }

    /// The 25 MiB `META_TARGET` used by the reference encoder. We pass it directly
    /// as `target_bytes` so this test measures the SAME partition target the real
    /// per-(chunk,mate) encode uses (`REF_POOL_TARGET = META_TARGET`).
    const META_TARGET: usize = 25 * 1024 * 1024;

    /// Deterministically build `n` pseudo-random ACGT reads of `len` bytes each,
    /// emulating unmapped reads (the only reads that take the literal fallback).
    /// Uses a constant-seeded LCG — `rand` is not a dependency of this crate.
    fn synth_literals(n: usize, len: usize, seed: u64) -> Vec<Vec<u8>> {
        let mut x = seed;
        let mut next = || {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (x >> 33) as usize
        };
        (0..n)
            .map(|_| (0..len).map(|_| b"ACGT"[next() & 3]).collect())
            .collect()
    }

    /// Sum the per-group `write_fallback_pool` sizes when `literals` is split into
    /// `k` contiguous read-order groups (exactly how the real per-chunk encode
    /// partitions fallback literals by chunk). Empty groups are skipped (the real
    /// encoder omits an empty (chunk,mate) FallbackPool entirely).
    fn per_chunk_total(literals: &[Vec<u8>], k: usize) -> usize {
        let mut total = 0usize;
        for g in literals.chunks(literals.len().div_ceil(k)) {
            if g.is_empty() {
                continue;
            }
            let mut out = Vec::new();
            write_fallback_pool(&mut out, &as_slices(g), META_TARGET).unwrap();
            total += out.len();
        }
        total
    }

    /// Task 1.4 ratio gate (Test A): characterize the worst-case RELATIVE
    /// fragmentation cost of splitting the fallback literals from ONE global pool
    /// into per-(chunk,mate) pools. The OLD global-pool size is exactly
    /// `write_fallback_pool(all_literals, META_TARGET)`, so both "before" (global)
    /// and "after" (per-chunk) are computed from the SAME literals with the SAME
    /// current code — no old binary needed.
    ///
    /// MEASURED (256 reads × 110 bp, this code, this box):
    ///   K=8   : global = 7670 B, per_chunk_total =  8398 B, overhead =  +9.49%
    ///   K=64  : (all-fallback worst case) per_chunk_total = 12530 B, overhead = +63.36%
    ///   K=256 : (every read its own chunk)               = 25028 B, overhead = +226.31%
    /// The K=8 figure is the realistic chunking; K=64/256 are the unrealistic
    /// all-fallback extremes (reported informationally, NOT gated).
    ///
    /// NOTE: this overhead is POOL-RELATIVE (denominator = fallback-pool bytes).
    /// The ARCHIVE-relative regression is this number scaled DOWN by the pool's
    /// fraction of the whole archive — bounded by Test B
    /// (`reference_fallback_pool_is_tiny_fraction_of_realistic_archive`).
    #[test]
    fn perchunk_pool_fragmentation_overhead_bounded() {
        // 256 pseudo-random unmapped reads ~110 bp — a representative fallback set.
        let literals = synth_literals(256, 110, 0xDEADBEEF);

        // "Before": the OLD single global pool over all literals.
        let mut global_buf = Vec::new();
        write_fallback_pool(&mut global_buf, &as_slices(&literals), META_TARGET).unwrap();
        let global = global_buf.len();

        // "After": per-chunk pools at a realistic chunking (K=8 groups of 32).
        let k = 8;
        let per_chunk = per_chunk_total(&literals, k);
        let overhead = (per_chunk as f64 - global as f64) / global as f64;

        // All-fallback EXTREMES (informational only — never realistic): every read
        // in its own tiny chunk. These are the absolute worst-case fragmentation.
        let extreme_64 = per_chunk_total(&literals, 64);
        let extreme_256 = per_chunk_total(&literals, 256);
        let ov_64 = (extreme_64 as f64 - global as f64) / global as f64;
        let ov_256 = (extreme_256 as f64 - global as f64) / global as f64;

        println!(
            "[fragmentation] global={global} B  per_chunk(K={k})={per_chunk} B  overhead={:.2}%",
            overhead * 100.0
        );
        println!(
            "[fragmentation] all-fallback worst case: K=64 -> {extreme_64} B (+{:.2}%), K=256 -> {extreme_256} B (+{:.2}%)",
            ov_64 * 100.0,
            ov_256 * 100.0
        );

        // BUDGET: the realistic (K=8) pool-relative fragmentation overhead measured
        // at +9.49%. We gate at 0.30 (30%) — ~3.2× headroom — so an unrelated future
        // change has to massively worsen pool fragmentation before this trips.
        const BUDGET: f64 = 0.30;
        assert!(
            overhead <= BUDGET,
            "per-chunk pool fragmentation overhead {:.4} exceeds BUDGET {BUDGET} \
             (global={global} B, per_chunk K={k}={per_chunk} B)",
            overhead
        );
    }
}
