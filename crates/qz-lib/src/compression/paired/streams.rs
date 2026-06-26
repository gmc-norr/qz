use crate::compression::bsc;
use anyhow::{Result, bail};
use flate2::Crc;

const PREFIX: usize = 12; // block_len(4)+record_count(4)+crc32(4)
// Single definition of the attacker-facing per-block size cap, shared with the
// single-end path (and reference, via the same const).
pub(crate) use crate::compression::decompress_impl::BSC_MAX_BLOCK_LEN as MAX_BLOCK;

/// Serialize blocks as `[num_blocks u32]` then per block `write_block_frame_header`
/// (block_len | record_count | crc32) + block bytes. Mirrors compress_impl's
/// `write_blocks_v4`. Empty input writes 4 zero bytes.
pub fn write_block_stream(out: &mut Vec<u8>, blocks: &[(u32, Vec<u8>)]) {
    out.extend_from_slice(&(blocks.len() as u32).to_le_bytes());
    for (rc, block) in blocks {
        bsc::write_block_frame_header(out, block.len() as u32, *rc, block);
        out.extend_from_slice(block);
    }
}

/// Parse the v4 frame and return CRC-verified `(record_count, raw_payload)` per
/// block. Does **no** codec decode — the caller knows the role's codec.
pub fn decode_block_payloads(data: &[u8]) -> Result<Vec<(u32, Vec<u8>)>> {
    if data.len() < 4 {
        bail!("block stream truncated");
    }
    let nb = u32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
    // `nb` is attacker-controlled (up to 4.29e9); each block needs >= PREFIX bytes
    // on disk, so bound it by the remaining bytes BEFORE with_capacity to avoid a
    // multi-GB over-allocation (abort) on hostile input.
    if nb > (data.len() - 4) / PREFIX {
        bail!(
            "block stream claims {nb} blocks but only {} bytes remain",
            data.len() - 4
        );
    }
    let mut o = 4usize;
    let mut out = Vec::with_capacity(nb);
    for _ in 0..nb {
        let hdr = data
            .get(o..o + PREFIX)
            .ok_or_else(|| anyhow::anyhow!("block header truncated"))?;
        let blen = u32::from_le_bytes(hdr[0..4].try_into().unwrap()) as usize;
        let rc = u32::from_le_bytes(hdr[4..8].try_into().unwrap());
        let crc = u32::from_le_bytes(hdr[8..12].try_into().unwrap());
        if blen > MAX_BLOCK {
            bail!("block_len {blen} exceeds cap");
        }
        o += PREFIX;
        let block = data
            .get(o..o + blen)
            .ok_or_else(|| anyhow::anyhow!("block body truncated"))?;
        o += blen;
        // CRC domain: record_count || block — matches compute_block_frame_crc in bsc.rs (flate2::Crc)
        let mut h = Crc::new();
        h.update(&rc.to_le_bytes());
        h.update(block);
        if h.sum() != crc {
            bail!("block CRC32 mismatch");
        }
        out.push((rc, block.to_vec()));
    }
    if o != data.len() {
        bail!("trailing bytes after block stream");
    }
    Ok(out)
}

/// Convenience for **BSC-payload** streams (sequences, BSC quality, the DeltaV1 op
/// stream): verify frames, then BSC-decompress and concatenate. Do NOT use for
/// columnar-header streams.
///
/// `max_total` caps the cumulative decompressed output: each block is already
/// bounded at 64 MiB by `bsc::decompress`, but without a running total a hostile
/// paired archive could chain many valid blocks into a multi-GB allocation. Pass
/// a content-derived bound (e.g. `codecs::stream_decode_cap(record_count)`).
pub fn decode_bsc_stream(data: &[u8], max_total: usize) -> Result<Vec<u8>> {
    use crate::compression::decompress_impl::MAX_PARALLEL_DECOMPRESS;
    use rayon::prelude::*;

    let payloads = decode_block_payloads(data)?;

    // Single-block streams (small columns, the R2 header-delta op stream, etc.):
    // decode inline. A `par_iter` over one block buys nothing and would strand a
    // freed buffer in a rayon worker's mimalloc heap for the pool's lifetime.
    if payloads.len() <= 1 {
        let mut out = Vec::new();
        for (_rc, block) in &payloads {
            let chunk = bsc::decompress(block)?;
            if out.len().saturating_add(chunk.len()) > max_total {
                bail!("paired BSC stream exceeds cap of {max_total} bytes (decompression bomb?)");
            }
            out.extend_from_slice(&chunk);
        }
        return Ok(out);
    }

    // Multi-block streams (sequences = the decode long pole, columnar header
    // columns): BSC-decode blocks in parallel. This is the paired analogue of the
    // single-end bounded producer — without it paired decode ran the whole block
    // loop on one core (~1.2 CPUs). Concurrency is bounded to
    // `MAX_PARALLEL_DECOMPRESS` blocks in flight (each already ≤64 MiB via
    // `bsc::decompress`), and the running `max_total` cap is re-checked between
    // batches so a decompression-bomb archive still bails before unbounded
    // allocation. `collect()` preserves block order ⇒ the concatenation stays in
    // original record order (byte-identical output).
    let mut out = Vec::new();
    for batch in payloads.chunks(MAX_PARALLEL_DECOMPRESS) {
        let decoded: Vec<Vec<u8>> = batch
            .par_iter()
            .map(|(_rc, block)| bsc::decompress(block))
            .collect::<Result<Vec<_>>>()?;
        for chunk in &decoded {
            if out.len().saturating_add(chunk.len()) > max_total {
                bail!("paired BSC stream exceeds cap of {max_total} bytes (decompression bomb?)");
            }
            out.extend_from_slice(chunk);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn block_payloads_roundtrip_and_bsc_concat() {
        let raw_a = b"HELLO-WORLD-AAAA".repeat(40);
        let raw_b = b"second-block-bbbb".repeat(40);
        // Each payload is a SINGLE raw BSC block (bsc::compress, NOT compress_parallel,
        // which self-frames). This matches what compress_parallel_with_breakpoints /
        // compress_stream_to_bsc_blocks_v4 emit per block; decoded by bsc::decompress.
        let blk_a = crate::compression::bsc::compress(&raw_a).unwrap();
        let blk_b = crate::compression::bsc::compress(&raw_b).unwrap();
        let mut buf = Vec::new();
        write_block_stream(&mut buf, &[(40, blk_a), (40, blk_b)]);
        // CRC-verified raw payloads (NOT BSC-decoded):
        let payloads = decode_block_payloads(&buf).unwrap();
        assert_eq!(payloads.len(), 2);
        assert_eq!(payloads[0].0, 40); // record_count
        // BSC-stream convenience decode:
        let decoded = decode_bsc_stream(&buf, usize::MAX).unwrap();
        let mut expect = raw_a.clone();
        expect.extend_from_slice(&raw_b);
        assert_eq!(decoded, expect);
    }

    #[test]
    fn decode_block_payloads_rejects_inflated_count() {
        // num_blocks = 0xFFFFFFFF (4.29e9) with only a couple of trailing bytes:
        // must Err quickly (bounded by remaining bytes) instead of trying to
        // Vec::with_capacity(~4e9) and aborting.
        let buf = [0xFFu8, 0xFF, 0xFF, 0xFF, 0x00, 0x00];
        let err = decode_block_payloads(&buf).unwrap_err();
        assert!(err.to_string().contains("blocks"), "got: {err}");
    }
}
