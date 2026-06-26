//! Stream codec helpers: the shared per-record framing primitives and the
//! fqzcomp quality codec used by the v5 chunk-major compress/decompress paths.
//! Called by compress() and decompress() in the parent module.

use super::*;
use anyhow::Result;

// ── Shared stream-building helpers ────────────────────────────────────────

/// Per-record byte cap enforced by the bounded streaming encoder.
///
/// A single FASTQ record larger than this is rejected at compress time —
/// the bounded decoder's per-block buffer is capped at
/// `BSC_MAX_BLOCK_LEN = 64 MiB + 4 KiB`, and the encoder must guarantee
/// no record exceeds that.
pub(super) const PER_RECORD_MAX_BYTES: usize = 63 * 1024 * 1024;

/// Cumulative decompressed-output cap for a per-record BSC stream decoded from
/// an untrusted archive. A stream of `num_reads` records cannot legitimately
/// exceed `num_reads × PER_RECORD_MAX_BYTES` (every record is capped at 63 MiB
/// on encode), floored at one BSC block so tiny/empty streams still decode.
///
/// This bounds a stacked-block decompression bomb on the **in-memory fallback**
/// decode path, which materializes the whole stream at once. The
/// bounded-streaming path is already capped via its in-flight memory budget;
/// this closes the same hole on the fallback path, where `decompress_parallel`
/// otherwise passes `usize::MAX`.
pub(super) fn stream_decode_cap(num_reads: usize) -> usize {
    num_reads
        .saturating_mul(PER_RECORD_MAX_BYTES)
        .max(super::decompress_impl::BSC_MAX_BLOCK_LEN)
}

// ── fqzcomp quality codec ───────────────────────────────────────────────────

/// Reads per fqzcomp sub-chunk within a single fqz quality blob. The default
/// `quality_ctx_block_size` (see `cli.rs`) is intentionally kept EQUAL to this so
/// every default single-end fqz block is exactly one sub-chunk → the inner
/// `par_iter` in `decompress_qualities_fqzcomp` never fans out on the default path
/// (default single-end parallelism comes from the bounded-streaming Fqz producer
/// decoding multiple blocks concurrently). Compress and decompress MUST agree on
/// this value to reconstruct sub-chunk boundaries, so both sides reference this one
/// const.
const FQZCOMP_SUB_CHUNK: usize = 500_000;

/// Split quality strings into consecutive slices of at most `record_cap`
/// records, fqzcomp-compress each slice independently (each gets its OWN local
/// mean-quality sort permutation, scoped to that ≤`record_cap`-record slice),
/// and return framed `(record_count, blob)` tuples in input order.
///
/// This is the **record-capped** splitter for single-end fqz encode. It exists
/// so each emitted block can be decoded in isolation (the sort permutation is
/// local to the block, not global to the whole chunk) — the prerequisite for
/// bounded-streaming fqz decode. It is intentionally distinct from the
/// **byte-capped** `fqz_blocks_capped_quals` used by reference mode (which
/// bisects by output size, not record count).
///
/// `record_cap` is clamped to ≥1. Empty input yields zero blocks.
///
/// The framed per-block `record_count` is the true number of records in that
/// block (callers must reject `None`-quality records before slicing, so the
/// borrowed `qual_strs` length equals the record count). The bounded-streaming
/// fqz decoder relies on this: it sums each block's `record_count` into the
/// total record count to drive emission.
///
/// Test-only: the production single-end/paired fqz path is the streaming
/// `QualStager::Fqz`, which reimplements this slicing inline; this reference
/// implementation is retained solely as the byte-identity oracle the stager is
/// tested against (`stream_compress::fqz_blocks_match_record_capped`).
#[cfg(test)]
pub(super) fn fqz_record_capped_blocks(
    qual_strs: &[&[u8]],
    record_cap: usize,
) -> Result<Vec<(u32, Vec<u8>)>> {
    use rayon::prelude::*;

    if qual_strs.is_empty() {
        return Ok(Vec::new());
    }
    let cap = record_cap.max(1);
    qual_strs
        .par_chunks(cap)
        .map(|slice| {
            let blob = compress_qualities_fqzcomp_quals(slice)?;
            Ok((slice.len() as u32, blob))
        })
        .collect()
}

/// Quality-slice entry point for the fqzcomp codec: mean-quality stable sort then
/// fqzcomp strat=0. Takes borrowed quality strings directly so callers need not
/// materialize throwaway `FastqRecord`s with reconstructed sequences just to reach
/// the quality field — fqzcomp never reads the sequence.
pub(super) fn compress_qualities_fqzcomp_quals(qual_strs: &[&[u8]]) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    if qual_strs.is_empty() {
        return Ok(Vec::new());
    }

    // htscodecs' fqz encoder rejects any record with len <= 0 (it returns NULL, which
    // would otherwise surface as the opaque "fqz_compress returned NULL"). A valid FASTQ
    // can contain an empty (zero-length) read, so detect it here and fail with a clear,
    // actionable diagnostic. fqzcomp has no zero-length representation; such inputs must
    // use the BSC quality backend instead.
    if let Some(idx) = qual_strs.iter().position(|q| q.is_empty()) {
        anyhow::bail!(
            "fqzcomp quality codec cannot encode an empty (zero-length) quality string \
             (first at record {idx} of this block): the input contains an empty read. The \
             context-modeled fqzcomp backend (auto-selected for lossless inputs with many \
             reads) has no zero-length representation — recompress with the BSC quality \
             backend (set `CompressConfig.quality_compressor = QualityCompressor::Bsc`)."
        );
    }

    // Compute mean quality per read (1 byte each)
    let sort_keys: Vec<u8> = qual_strs
        .iter()
        .map(|q| {
            let sum: u64 = q.iter().map(|&b| b as u64).sum();
            (sum / q.len().max(1) as u64) as u8
        })
        .collect();

    // Sort by mean quality (stable sort)
    let mut indices: Vec<usize> = (0..qual_strs.len()).collect();
    indices.sort_by_key(|&i| sort_keys[i]);

    // Reorder quality strings
    let sorted_quals: Vec<&[u8]> = indices.iter().map(|&i| qual_strs[i]).collect();

    let num_sub_chunks = sorted_quals.len().div_ceil(FQZCOMP_SUB_CHUNK);
    info!(
        "fqzcomp: {} reads, sorting by mean quality, {} sub-chunks, compressing with strat=0",
        qual_strs.len(),
        num_sub_chunks
    );

    // Compress sub-chunks in parallel with rayon, and BSC-compress sort keys concurrently
    let (sub_chunk_blobs, sort_keys_bsc) = rayon::join(
        || -> Result<Vec<Vec<u8>>> {
            sorted_quals
                .par_chunks(FQZCOMP_SUB_CHUNK)
                .map(|chunk| fqzcomp::compress(chunk, 0))
                .collect::<Result<Vec<_>>>()
        },
        || bsc::compress_parallel_adaptive(&sort_keys),
    );
    let sub_chunk_blobs = sub_chunk_blobs?;
    let sort_keys_bsc = sort_keys_bsc?;

    let fqzcomp_total: usize = sub_chunk_blobs.iter().map(|b| b.len()).sum();
    info!(
        "fqzcomp: sort keys {} B (BSC'd), fqzcomp {} B in {} sub-chunks",
        sort_keys_bsc.len(),
        fqzcomp_total,
        sub_chunk_blobs.len()
    );

    // Pack: [num_reads][sort_keys_bsc_len][sort_keys_bsc][num_sub_chunks][sub_chunk_len, data]...
    let mut output =
        Vec::with_capacity(12 + sort_keys_bsc.len() + fqzcomp_total + sub_chunk_blobs.len() * 4);
    output.extend_from_slice(&(qual_strs.len() as u32).to_le_bytes());
    output.extend_from_slice(&(sort_keys_bsc.len() as u32).to_le_bytes());
    output.extend_from_slice(&sort_keys_bsc);
    output.extend_from_slice(&(sub_chunk_blobs.len() as u32).to_le_bytes());
    for blob in &sub_chunk_blobs {
        output.extend_from_slice(&(blob.len() as u32).to_le_bytes());
        output.extend_from_slice(blob);
    }

    Ok(output)
}

/// Decompress fqzcomp quality data with mean-quality unsort.
///
/// Format: [num_reads: u32][sort_keys_bsc_len: u32][sort_keys_bsc]
///         [num_sub_chunks: u32][sub_chunk_len: u32, sub_chunk_data]...
pub(crate) fn decompress_qualities_fqzcomp(compressed: &[u8]) -> Result<Vec<u8>> {
    use rayon::prelude::*;

    // Untrusted-input DoS guard: `num_reads` comes straight from the blob header and
    // sizes several O(num_reads) allocations below (`lengths` 4 B/read, `indices` 8 B,
    // `unsorted_quals` 16 B). A hostile header could otherwise claim up to u32::MAX
    // (~4.3 B) reads and force a ~68 GB allocation → unrecoverable process abort. No
    // legitimate fqzcomp blob carries more than this many reads (reference sub-blocks
    // are ≤ quality_ctx_block_size; single-end blobs are ≤ one chunk), so anything
    // larger is malformed. Bounds the post-check allocations to a clean `Err`.
    // (The separate concern of an over-expanding BSC `sort_keys` decompression is a
    // BSC-layer issue shared by several decode call sites, tracked separately.)
    const FQZ_MAX_READS_PER_BLOB: usize = 1 << 28; // 268_435_456

    if compressed.is_empty() {
        return Ok(Vec::new());
    }

    if compressed.len() < 12 {
        anyhow::bail!("fqzcomp quality data too short");
    }

    let mut offset = 0;

    // Parse header
    let num_reads = read_le_u32(compressed, offset)? as usize;
    if num_reads > FQZ_MAX_READS_PER_BLOB {
        anyhow::bail!(
            "fqzcomp: num_reads {num_reads} exceeds sanity cap {FQZ_MAX_READS_PER_BLOB} (malformed)"
        );
    }
    offset += 4;
    let sort_keys_len = read_le_u32(compressed, offset)? as usize;
    offset += 4;

    info!(
        "fqzcomp decompress: {} reads, sort_keys_bsc {} bytes, total blob {} bytes",
        num_reads,
        sort_keys_len,
        compressed.len()
    );

    if offset + sort_keys_len > compressed.len() {
        anyhow::bail!("fqzcomp: truncated sort keys");
    }

    // BSC-decompress sort keys
    let sort_keys_data = &compressed[offset..offset + sort_keys_len];
    offset += sort_keys_len;

    // Parse sub-chunk table
    if offset + 4 > compressed.len() {
        anyhow::bail!("fqzcomp: truncated sub-chunk header");
    }
    let num_sub_chunks = read_le_u32(compressed, offset)? as usize;
    offset += 4;

    // The sub-chunk count must match what the compressor would have produced
    // for `num_reads`. A mismatched (corrupt) header would otherwise underflow
    // `num_reads - start` below and panic.
    let expected_sub_chunks = num_reads.div_ceil(FQZCOMP_SUB_CHUNK);
    if num_sub_chunks != expected_sub_chunks {
        anyhow::bail!(
            "fqzcomp: sub-chunk count {} inconsistent with {} reads (expected {})",
            num_sub_chunks,
            num_reads,
            expected_sub_chunks
        );
    }

    // Collect sub-chunk slices
    let mut sub_chunk_slices: Vec<&[u8]> = Vec::with_capacity(num_sub_chunks);
    for i in 0..num_sub_chunks {
        if offset + 4 > compressed.len() {
            anyhow::bail!("fqzcomp: truncated sub-chunk {} header", i);
        }
        let chunk_len = read_le_u32(compressed, offset)? as usize;
        offset += 4;
        if offset + chunk_len > compressed.len() {
            anyhow::bail!("fqzcomp: truncated sub-chunk {} data", i);
        }
        sub_chunk_slices.push(&compressed[offset..offset + chunk_len]);
        offset += chunk_len;
    }

    // Check 1: exact outer consumption — no trailing bytes after all sub-chunks.
    if offset != compressed.len() {
        anyhow::bail!(
            "fqzcomp: {} trailing bytes after sub-chunks",
            compressed.len() - offset
        );
    }

    // Compute reads per sub-chunk (must match compression side)
    let reads_per_sub: Vec<usize> = (0..num_sub_chunks)
        .map(|i| {
            let start = i * FQZCOMP_SUB_CHUNK;
            FQZCOMP_SUB_CHUNK.min(num_reads - start)
        })
        .collect();

    info!(
        "fqzcomp decompress: {} sub-chunks, reads per sub: {:?}",
        num_sub_chunks, &reads_per_sub
    );

    // Decompress sort keys and fqzcomp sub-chunks in parallel. The sort_keys stream is
    // exactly one byte per read, so a valid decode is `num_reads` bytes (verified by the
    // `!= num_reads` check below). Cap the BSC decode at `num_reads` (already bounded by
    // FQZ_MAX_READS_PER_BLOB) so a crafted/over-expanding sort_keys block can't force an
    // unbounded allocation before that check runs.
    //
    // Only fan the sub-chunk decode across the rayon pool when there is more than one
    // sub-chunk. A blob with a single sub-chunk (every reference per-chunk-mate quality
    // entry, and any quality block ≤ FQZCOMP_SUB_CHUNK = 500K reads) gains nothing from
    // `par_iter` but DOES pull in an extra worker thread, and each freshly-touched
    // worker pays a one-time per-thread allocator-arena + stack RSS cost that is never
    // reclaimed for the pool's lifetime. Reference decode runs one such blob per chunk,
    // so spreading them across the whole pool made decode peak RSS ramp with chunk count
    // (more chunks → more distinct workers warmed). Decoding the lone sub-chunk inline
    // keeps that work on the current thread, so the warmed-worker set — and thus peak
    // RSS — stays bounded and constant in read/chunk count.
    //
    // The inner par_iter only fans out when a single fqz block holds >1 sub-chunk
    // (i.e. a non-default `quality_ctx_block_size` > FQZCOMP_SUB_CHUNK). On the
    // DEFAULT single-end path every block is exactly one sub-chunk
    // (quality_ctx_block_size == FQZCOMP_SUB_CHUNK == 500K), so this loop never fans
    // out — default single-end parallelism comes from the bounded-streaming Fqz
    // producer decoding multiple blocks concurrently, not from this rayon::join.
    let (sort_keys_result, sub_results): (Result<Vec<u8>>, Vec<Result<(Vec<u8>, Vec<i32>)>>) =
        if num_sub_chunks <= 1 {
            // Inline: no extra worker is recruited for a single sub-chunk.
            let sub: Vec<Result<(Vec<u8>, Vec<i32>)>> = sub_chunk_slices
                .iter()
                .zip(reads_per_sub.iter())
                .map(|(data, &n_reads)| fqzcomp::decompress(data, n_reads))
                .collect();
            (bsc::decompress_parallel_max(sort_keys_data, num_reads), sub)
        } else {
            rayon::join(
                || bsc::decompress_parallel_max(sort_keys_data, num_reads),
                || {
                    sub_chunk_slices
                        .par_iter()
                        .zip(reads_per_sub.par_iter())
                        .map(|(data, &n_reads)| fqzcomp::decompress(data, n_reads))
                        .collect()
                },
            )
        };
    let sort_keys = sort_keys_result?;
    // Check 2: exact sort-key count — reject both too-few and too-many.
    if sort_keys.len() != num_reads {
        anyhow::bail!(
            "fqzcomp: decoded {} sort keys for {} reads",
            sort_keys.len(),
            num_reads
        );
    }

    // Concatenate decoded qualities and lengths from all sub-chunks
    let mut decoded_quals = Vec::new();
    let mut lengths: Vec<i32> = Vec::with_capacity(num_reads);
    for res in sub_results {
        let (quals, lens) = res?;
        decoded_quals.extend_from_slice(&quals);
        lengths.extend_from_slice(&lens);
    }
    // Check 3: exact length count — reject both too-few and too-many.
    if lengths.len() != num_reads {
        anyhow::bail!(
            "fqzcomp: decoded {} read lengths for {} reads",
            lengths.len(),
            num_reads
        );
    }

    // Reconstruct permutation (same stable sort as compression)
    let mut indices: Vec<usize> = (0..num_reads).collect();
    indices.sort_by_key(|&i| sort_keys[i]);

    // Split decoded qualities into individual reads
    let mut sorted_quals: Vec<&[u8]> = Vec::with_capacity(num_reads);
    let mut q_offset = 0;
    for i in 0..num_reads {
        // Check 4: checked length conversion — reject negative lengths without panic.
        let len = usize::try_from(lengths[i])
            .map_err(|_| anyhow::anyhow!("fqzcomp: negative read length {}", lengths[i]))?;
        if q_offset + len > decoded_quals.len() {
            anyhow::bail!("fqzcomp: truncated quality data at read {}", i);
        }
        sorted_quals.push(&decoded_quals[q_offset..q_offset + len]);
        q_offset += len;
    }

    // Check 5: exact quality consumption — no trailing bytes in decoded quality data.
    if q_offset != decoded_quals.len() {
        anyhow::bail!(
            "fqzcomp: {} trailing quality bytes",
            decoded_quals.len() - q_offset
        );
    }

    // Inverse permute: unsorted[indices[i]] = sorted[i]
    let mut unsorted_quals: Vec<&[u8]> = vec![&[]; num_reads];
    for (sorted_idx, &original_idx) in indices.iter().enumerate() {
        unsorted_quals[original_idx] = sorted_quals[sorted_idx];
    }

    // Build varint-prefixed stream (raw ASCII, compatible with QualityBinning::None)
    let total_est: usize = unsorted_quals.iter().map(|q| q.len() + 2).sum();
    let mut output = Vec::with_capacity(total_est);
    for qual in &unsorted_quals {
        write_varint(&mut output, qual.len())?;
        output.extend_from_slice(qual);
    }

    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fqzcomp_decode_rejects_inconsistent_sub_chunk_count_without_panic() {
        // Adversarial blob: header claims 1 read but 5 sub-chunks. The
        // per-sub-chunk read count math (`num_reads - start`) would underflow
        // and panic; a corrupt archive must yield an Err instead.
        let mut blob = Vec::new();
        blob.extend_from_slice(&1u32.to_le_bytes()); // num_reads = 1
        blob.extend_from_slice(&0u32.to_le_bytes()); // sort_keys_len = 0
        blob.extend_from_slice(&5u32.to_le_bytes()); // num_sub_chunks = 5
        for _ in 0..5 {
            blob.extend_from_slice(&0u32.to_le_bytes()); // each sub-chunk len = 0
        }
        assert!(decompress_qualities_fqzcomp(&blob).is_err());
    }

    #[test]
    fn empty_quality_record_reports_clear_actionable_error() {
        // htscodecs' fqz encoder rejects any record with len <= 0 (returns NULL).
        // A valid FASTQ may contain an empty (zero-length) read, which the standard
        // reader accepts. The fqz quality path must fail with a CLEAR, actionable
        // message naming the empty read — not the opaque "fqz_compress returned NULL".
        let quals: Vec<&[u8]> = vec![b"IIIIIIII", b"", b"FFFFFFFF"];
        let err = compress_qualities_fqzcomp_quals(&quals)
            .expect_err("empty quality must be rejected")
            .to_string();
        assert!(
            err.contains("empty") || err.contains("zero-length"),
            "expected a clear empty-quality diagnostic, got: {err}"
        );
        assert!(
            !err.contains("returned NULL"),
            "must not surface the opaque htscodecs NULL error, got: {err}"
        );
    }

    #[test]
    fn test_fqz_record_capped_blocks_splits_and_roundtrips() {
        // 25 reads, cap 10 -> 3 blocks (10, 10, 5), each independently sorted.
        let quals: Vec<Vec<u8>> = (0..25u8)
            .map(|i| {
                // Distinct mean qualities so the per-block local sort actually
                // permutes; verifies decode unsorts per-block correctly.
                let q = b'!' + (i % 40);
                vec![q; 16]
            })
            .collect();
        let refs: Vec<&[u8]> = quals.iter().map(|q| q.as_slice()).collect();

        let blocks = fqz_record_capped_blocks(&refs, 10).expect("encode ok");
        assert_eq!(blocks.len(), 3, "25 reads / cap 10 => 3 blocks");
        assert_eq!(
            blocks.iter().map(|(c, _)| *c).collect::<Vec<_>>(),
            vec![10u32, 10, 5],
            "real per-block record counts"
        );

        // Each block is a self-contained single-blob fqz stream (the unit the v5
        // bounded-streaming Fqz producer decodes). Decoding each block and
        // concatenating in block order reconstructs every quality string in
        // input order (the per-block local sort is reversed inside the decoder).
        let mut decoded: Vec<Vec<u8>> = Vec::new();
        for (_, blob) in &blocks {
            let raw = decompress_qualities_fqzcomp(blob).expect("decode ok");
            // The blob decodes to a varint-length-prefixed stream of its records.
            let mut pos = 0usize;
            while pos < raw.len() {
                let len = read_varint(&raw, &mut pos).expect("varint");
                decoded.push(raw[pos..pos + len].to_vec());
                pos += len;
            }
        }
        assert_eq!(decoded, quals, "lossless input-order roundtrip");
    }

    #[test]
    fn test_fqz_record_capped_blocks_empty_and_single() {
        // Empty -> zero blocks.
        assert!(fqz_record_capped_blocks(&[], 10).unwrap().is_empty());
        // Fewer reads than cap -> exactly one block with the real count.
        let quals: Vec<Vec<u8>> = (0..3).map(|_| vec![b'I'; 8]).collect();
        let refs: Vec<&[u8]> = quals.iter().map(|q| q.as_slice()).collect();
        let blocks = fqz_record_capped_blocks(&refs, 10).unwrap();
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].0, 3);
    }

    #[test]
    fn test_fqzcomp_decode_rejects_implausible_num_reads_without_alloc() {
        // Adversarial blob: header claims ~4.3 billion reads. Without the sanity cap
        // this would drive a multi-GB allocation (and likely process abort) before any
        // real validation. It must return an Err cheaply, never allocating O(num_reads).
        let mut blob = Vec::new();
        blob.extend_from_slice(&u32::MAX.to_le_bytes()); // num_reads = 4_294_967_295
        blob.extend_from_slice(&0u32.to_le_bytes()); // sort_keys_len = 0
        blob.extend_from_slice(&0u32.to_le_bytes()); // num_sub_chunks = 0
        let err = decompress_qualities_fqzcomp(&blob).unwrap_err();
        assert!(
            err.to_string().contains("sanity cap"),
            "expected sanity-cap rejection, got: {err}"
        );
    }
}
