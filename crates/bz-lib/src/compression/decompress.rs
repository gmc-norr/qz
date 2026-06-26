use crate::cli::DecompressConfig;
use crate::compression::archive::{ArchiveHeader, CHUNK_FLAG_FQZ_QUALITY, ChunkHeader};
use crate::compression::md_codec;
use crate::compression::streams::{self, ChunkConsensus, NUM_STREAMS};
use crate::io::bam::BamWriter;
use anyhow::{Context, Result};
use flate2::Crc;
use qz_lib::compression::bsc;
use qz_lib::compression::fqzcomp;
use std::collections::{BTreeMap, VecDeque};
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::time::Instant;
use tracing::info;

// --- BZ_DECODE_TIMING: env-gated cumulative per-phase decode timers ---
// Profiling aid, off by default. Set BZ_DECODE_TIMING=1 to print cumulative
// wall-time per decode phase, summed across all chunk-decode worker threads.
// Because phases run inside per-chunk threads (a window of chunks decode
// concurrently), these are CPU-seconds-per-phase, not critical-path; dividing a
// phase's total by the overall wall gives its effective concurrency, which is
// what pins the serial/narrow bottleneck. Zero overhead when unset.
use std::sync::LazyLock;
use std::sync::atomic::{AtomicU64, Ordering};
static DECODE_TIMING: LazyLock<bool> =
    LazyLock::new(|| std::env::var("BZ_DECODE_TIMING").is_ok_and(|v| v != "0" && !v.is_empty()));
static T_BSC: AtomicU64 = AtomicU64::new(0);
static T_P1: AtomicU64 = AtomicU64::new(0);
static T_P2: AtomicU64 = AtomicU64::new(0);
static T_P3: AtomicU64 = AtomicU64::new(0);
static T_P4: AtomicU64 = AtomicU64::new(0);

/// Time `f`, adding its elapsed micros to `acc` only when `BZ_DECODE_TIMING` is
/// set. The `Instant::now()` is taken unconditionally (≈20 ns) so call sites stay
/// uncluttered; the atomic add (Relaxed, ~once per phase per chunk) is gated.
#[inline]
fn timed<T>(acc: &AtomicU64, f: impl FnOnce() -> T) -> T {
    let t = Instant::now();
    let r = f();
    if *DECODE_TIMING {
        acc.fetch_add(t.elapsed().as_micros() as u64, Ordering::Relaxed);
    }
    r
}

/// Validate decoded per-record alignment fields against the reference
/// dictionary. `ref_id`/`next_ref_id` must be -1 (unmapped) or a valid index
/// `[0, n_ref)`, and `pos` must be >= -1. A reconstructed record that violates
/// these would produce a structurally-invalid BAM that crashes downstream tools
/// (samtools); reject it here with a clear error instead. This is defense in
/// depth — on-disk corruption is already caught by the per-chunk compressed
/// CRC, and reconstruction divergence by the per-chunk content CRC.
#[inline]
fn validate_alignment_fields(
    rec_i: usize,
    ref_id: i32,
    next_ref_id: i32,
    pos: i32,
    n_ref: i32,
) -> Result<()> {
    if ref_id < -1 || ref_id >= n_ref {
        anyhow::bail!("record {rec_i}: ref_id {ref_id} out of range [-1, {n_ref})");
    }
    if next_ref_id < -1 || next_ref_id >= n_ref {
        anyhow::bail!("record {rec_i}: next_ref_id {next_ref_id} out of range [-1, {n_ref})");
    }
    if pos < -1 {
        anyhow::bail!("record {rec_i}: pos {pos} < -1");
    }
    Ok(())
}

/// Read a little-endian i32 from a stream at the given offset, advancing it by 4.
/// Returns an error if the stream is too short.
#[inline]
fn read_i32(stream: &[u8], off: &mut usize, name: &str) -> Result<i32> {
    let end = off
        .checked_add(4)
        .filter(|&e| e <= stream.len())
        .ok_or_else(|| anyhow::anyhow!("Truncated {name} stream at offset {}", *off))?;
    let val = i32::from_le_bytes(stream[*off..end].try_into().unwrap());
    *off = end;
    Ok(val)
}

/// Read a little-endian u16 from a stream at the given offset, advancing it by 2.
#[inline]
fn read_u16(stream: &[u8], off: &mut usize, name: &str) -> Result<u16> {
    let end = off
        .checked_add(2)
        .filter(|&e| e <= stream.len())
        .ok_or_else(|| anyhow::anyhow!("Truncated {name} stream at offset {}", *off))?;
    let val = u16::from_le_bytes(stream[*off..end].try_into().unwrap());
    *off = end;
    Ok(val)
}

/// Read a little-endian u32 from a stream at the given offset, advancing it by 4.
#[inline]
fn read_u32(stream: &[u8], off: &mut usize, name: &str) -> Result<u32> {
    let end = off
        .checked_add(4)
        .filter(|&e| e <= stream.len())
        .ok_or_else(|| anyhow::anyhow!("Truncated {name} stream at offset {}", *off))?;
    let val = u32::from_le_bytes(stream[*off..end].try_into().unwrap());
    *off = end;
    Ok(val)
}

/// Read a single byte from a stream at the given offset, advancing it by 1.
#[inline]
fn read_byte(stream: &[u8], off: &mut usize, name: &str) -> Result<u8> {
    if *off >= stream.len() {
        anyhow::bail!("Truncated {name} stream at offset {}", *off);
    }
    let val = stream[*off];
    *off += 1;
    Ok(val)
}

/// Read all compressed stream data for a chunk from the reader and verify CRC32.
///
/// The CRC32 is computed over the concatenated compressed stream payloads and
/// compared to `chunk_header.crc32`.  A mismatch means the archive has been
/// corrupted on disk; we bail before attempting (potentially slow) BSC decompression.
pub(super) fn read_chunk_streams<R: Read>(
    reader: &mut R,
    chunk_header: &ChunkHeader,
) -> Result<Vec<Vec<u8>>> {
    let mut compressed_streams: Vec<Vec<u8>> = Vec::with_capacity(NUM_STREAMS);
    let mut crc = Crc::new();
    for i in 0..NUM_STREAMS {
        let size = chunk_header.stream_sizes[i] as usize;
        // Read incrementally via take() rather than pre-allocating `size` bytes.
        // `stream_sizes` comes straight from the (untrusted) archive header — a
        // corrupt value claiming a multi-GB stream must not drive a huge
        // speculative allocation before we discover the bytes aren't there.
        let mut data = Vec::new();
        let read = reader.by_ref().take(size as u64).read_to_end(&mut data)?;
        if read != size {
            anyhow::bail!(
                "Truncated chunk stream {i}: expected {size} bytes, got {read} — archive is corrupted or truncated"
            );
        }
        if size > 0 {
            crc.update(&data);
        }
        compressed_streams.push(data);
    }
    let computed = crc.sum();
    if computed != chunk_header.crc32 {
        anyhow::bail!(
            "Chunk CRC32 mismatch: stored {:08x}, computed {:08x} — archive is corrupted",
            chunk_header.crc32,
            computed
        );
    }
    Ok(compressed_streams)
}

/// Decompress a chunk using fqz — parallel sequence reconstruction,
/// then fqz decompression, then record assembly.
/// Parse raw 4-byte BAM CIGAR ops into (op_code, len) pairs.
fn parse_cigar_ops(bytes: &[u8]) -> Vec<(u8, u32)> {
    bytes
        .chunks_exact(4)
        .map(|c| {
            let v = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
            ((v & 0xf) as u8, v >> 4)
        })
        .collect()
}

/// Regenerate any derived aux tags (NM:i, MC:Z, MD:Z) and re-insert them into
/// `stored_aux`, byte-exact. Returns `None` when nothing was derived (caller uses
/// `stored_aux` as-is). Compress strips MD, MC, NM (in that order), so decode
/// inserts the reverse: NM, MC, MD. Advances each index/type cursor only for
/// present tags. `nm_present` implies `md_present` (NM reuses MD's recovered
/// reference), so `ref_bases` is available whenever NM needs it.
#[allow(clippy::too_many_arguments)]
fn regen_aux(
    stored_aux: &[u8],
    md_present: bool,
    md_index_stream: &[u8],
    md_index_off: &mut usize,
    mc_present: bool,
    mc_index_stream: &[u8],
    mc_index_off: &mut usize,
    nm_present: bool,
    nm_index_stream: &[u8],
    nm_index_off: &mut usize,
    nm_type_stream: &[u8],
    nm_type_off: &mut usize,
    mate_cigar_bytes: Option<&[u8]>,
    ref_id: i32,
    pos: i32,
    cigar_bytes: &[u8],
    ascii_seq: &[u8],
    consensus: &ChunkConsensus,
    exceptions: &BTreeMap<(i32, i32), u8>,
) -> Result<Option<Vec<u8>>> {
    if !md_present && !mc_present && !nm_present {
        return Ok(None);
    }
    let cigar = parse_cigar_ops(cigar_bytes);
    // Recovered reference bases (ref order), shared by NM and MD. Computed only
    // when MD is derivable; NM derivation requires it.
    let ref_bases: Option<Vec<u8>> = if md_present {
        let span = md_codec::ref_span(&cigar);
        // Bound the reference span derived from the (untrusted) CIGAR before
        // allocating: ref_span sums ref-consuming op lengths with no cap, so a
        // crafted CIGAR of huge D/N ops can drive Vec::with_capacity(span) to a
        // process-aborting allocation (uncatchable — the chunk runs on a thread
        // whose join only recovers unwinding panics, not abort). No single linear
        // alignment legitimately spans this far (far above the largest gene / any
        // real long read), so reject it as a corrupt archive.
        const MAX_REF_SPAN: usize = 64 * 1024 * 1024;
        if span > MAX_REF_SPAN {
            anyhow::bail!(
                "regen_aux: implausible reference span {span} from CIGAR (corrupt archive)"
            );
        }
        let seg = consensus.get_segment_for_pos(ref_id, pos);
        let end = pos.saturating_add(span as i32);
        let exc: Vec<(i32, u8)> = exceptions
            .range((ref_id, pos)..(ref_id, end))
            .map(|(&(_, p), &b)| (p, b))
            .collect();
        let mut ei = 0usize;
        let mut rb = Vec::with_capacity(span);
        for k in 0..span {
            let p = pos + k as i32;
            let b = if ei < exc.len() && exc[ei].0 == p {
                let v = exc[ei].1;
                ei += 1;
                v
            } else {
                streams::nibble_to_ascii(seg.map_or(0, |s| s.get(p)))
            };
            rb.push(b);
        }
        Some(rb)
    } else {
        None
    };

    let mut cur: std::borrow::Cow<[u8]> = std::borrow::Cow::Borrowed(stored_aux);
    // NM first (stripped last on compress).
    if nm_present {
        let nm_index = streams::read_varint(nm_index_stream, nm_index_off)? as u32;
        let ty = *nm_type_stream
            .get(*nm_type_off)
            .ok_or_else(|| anyhow::anyhow!("NM type stream exhausted during decode"))?;
        *nm_type_off += 1;
        let rb = ref_bases
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("NM derived but MD absent (corrupt archive)"))?;
        let nm = md_codec::derive_nm(&cigar, rb, ascii_seq)
            .ok_or_else(|| anyhow::anyhow!("NM derivation failed during decode"))?;
        let item = md_codec::encode_int_item(*b"NM", ty, nm).ok_or_else(|| {
            anyhow::anyhow!("NM value {nm} doesn't fit stored type {}", ty as char)
        })?;
        let out = md_codec::insert_item(&cur, &item, nm_index)
            .ok_or_else(|| anyhow::anyhow!("NM re-insertion failed during decode"))?;
        cur = std::borrow::Cow::Owned(out);
    }
    // MC second (mate's CIGAR text).
    if mc_present {
        let mc_index = streams::read_varint(mc_index_stream, mc_index_off)? as u32;
        let mate_cig = mate_cigar_bytes
            .ok_or_else(|| anyhow::anyhow!("MC derived but mate not in chunk during decode"))?;
        let mc_text = md_codec::cigar_text(&parse_cigar_ops(mate_cig));
        let out = md_codec::insert_tag(&cur, *b"MC", &mc_text, mc_index)
            .ok_or_else(|| anyhow::anyhow!("MC re-insertion failed during decode"))?;
        cur = std::borrow::Cow::Owned(out);
    }
    // MD third (from consensus + exceptions).
    if md_present {
        let md_index = streams::read_varint(md_index_stream, md_index_off)? as u32;
        let rb = ref_bases
            .as_ref()
            .expect("ref_bases present when md_present");
        let md = md_codec::regen_md(&cigar, rb, ascii_seq);
        let out = md_codec::insert_tag(&cur, *b"MD", &md, md_index)
            .ok_or_else(|| anyhow::anyhow!("MD re-insertion failed during decode"))?;
        cur = std::borrow::Cow::Owned(out);
    }
    Ok(Some(cur.into_owned()))
}

/// Decode a v12 fqzcomp quality blob (multiblock framing from the compress side)
/// into per-read ASCII Phred+33, de-reorienting reverse-strand reads using
/// `strands` (record order). Returns the same Vec<Vec<u8>> contract the old
/// fqz decoder did, so the record-assembly write path is unchanged.
fn decompress_quality_fqz_multiblock(blob: &[u8], strands: &[bool]) -> Result<Vec<Vec<u8>>> {
    use rayon::prelude::*;

    if blob.is_empty() {
        return Ok(Vec::new());
    }
    if blob.len() < 4 {
        anyhow::bail!("fqz quality blob too short ({} bytes)", blob.len());
    }
    let num_blocks = u32::from_le_bytes(blob[0..4].try_into().unwrap()) as usize;

    // Pass 1 (serial, cheap): walk the block framing to locate each block's byte
    // range and its starting global read index. Block byte offsets are sequential
    // (each depends on the prior block_len), so this walk can't be parallel — but
    // it only reads 12-byte headers and does bounds arithmetic, no entropy decode.
    // All DoS guards (record_count vs remaining records, block_len overflow) match
    // the old serial path and stay BEFORE any allocation/decode. The capacity hint
    // is bounded by the smallest possible block (12-byte header) so an untrusted
    // `num_blocks` can't drive an unbounded reservation.
    struct BlockSpec {
        start: usize,
        end: usize,
        record_count: usize,
        crc: u32,
        rd_start: usize,
    }
    let cap_hint = num_blocks.min(blob.len().saturating_sub(4) / 12 + 1);
    let mut specs: Vec<BlockSpec> = Vec::with_capacity(cap_hint);
    let mut off = 4usize;
    let mut rd = 0usize; // global read index, for strand lookup
    for b in 0..num_blocks {
        if off + 12 > blob.len() {
            anyhow::bail!("fqz quality blob truncated at block {b} header");
        }
        let block_len = u32::from_le_bytes(blob[off..off + 4].try_into().unwrap()) as usize;
        let record_count = u32::from_le_bytes(blob[off + 4..off + 8].try_into().unwrap()) as usize;
        let crc = u32::from_le_bytes(blob[off + 8..off + 12].try_into().unwrap());
        off += 12;
        // DoS guard: `record_count` is untrusted and drives `vec![0; record_count]`
        // inside fqzcomp::decompress. The running total can't exceed the chunk's
        // real record count (`strands.len()`, already bounded upstream), so reject
        // an over-large block before any allocation. CRC is recomputable by an
        // attacker, so this bound must be independent of it.
        let remaining = strands.len() - rd;
        if record_count > remaining {
            anyhow::bail!(
                "fqz quality block {b}: record_count {record_count} exceeds remaining \
                 records {remaining} (corrupt archive)"
            );
        }
        let end = off
            .checked_add(block_len)
            .filter(|&e| e <= blob.len())
            .ok_or_else(|| anyhow::anyhow!("fqz quality block {b} overflows blob"))?;
        specs.push(BlockSpec {
            start: off,
            end,
            record_count,
            crc,
            rd_start: rd,
        });
        off = end;
        rd += record_count;
    }
    if rd != strands.len() {
        anyhow::bail!(
            "fqz quality: blocks cover {rd} reads, expected {}",
            strands.len()
        );
    }
    // Reject trailing bytes after the last block (verify enforces the same — a corrupt
    // archive must not slip extra bytes past the multiblock decoder).
    if off != blob.len() {
        anyhow::bail!(
            "fqz quality: {} trailing bytes after {num_blocks} blocks (corrupt archive)",
            blob.len() - off
        );
    }

    // Pass 2 (parallel): each block is fully independent — its own CRC, fqz blob,
    // and a disjoint output slice of reads (`rd_start..rd_start+record_count`) — so
    // decode them concurrently. This is the decode analogue of
    // `compress_quality_fqz_parallel`, which already encodes blocks in parallel; the
    // old serial loop here was the dominant decode CPU cost (~59% of decode work).
    // `collect` preserves block order, so output order is unchanged.
    let decode_block = |(b, spec): (usize, &BlockSpec)| -> Result<Vec<Vec<u8>>> {
        let fqz_blob = &blob[spec.start..spec.end];
        {
            use flate2::Crc;
            let mut h = Crc::new();
            h.update(&(spec.record_count as u32).to_le_bytes());
            h.update(fqz_blob);
            if h.sum() != spec.crc {
                anyhow::bail!("fqz quality block {b} CRC mismatch — archive corrupt");
            }
        }
        let (concat, lengths) = fqzcomp::decompress(fqz_blob, spec.record_count)?;
        if lengths.len() != spec.record_count {
            anyhow::bail!(
                "fqz quality block {b}: {} lengths, expected {}",
                lengths.len(),
                spec.record_count
            );
        }
        // Strand flags for this block's records (exactly `record_count` long).
        let strands_b = &strands[spec.rd_start..spec.rd_start + spec.record_count];
        let mut reads: Vec<Vec<u8>> = Vec::with_capacity(spec.record_count);
        let mut p = 0usize;
        for (i, &len) in lengths.iter().enumerate() {
            if len < 0 {
                anyhow::bail!("fqz quality block {b}: negative read length {len}");
            }
            let len = len as usize;
            let e = p
                .checked_add(len)
                .filter(|&e| e <= concat.len())
                .ok_or_else(|| anyhow::anyhow!("fqz quality block {b}: read overruns buffer"))?;
            let mut q = concat[p..e].to_vec();
            p = e;
            if strands_b[i] {
                q.reverse(); // de-reorient reverse-strand reads back to BAM order
            }
            for v in q.iter_mut() {
                *v = v.wrapping_add(33); // raw Phred -> ASCII Phred+33
            }
            reads.push(q);
        }
        Ok(reads)
    };
    let per_block: Vec<Vec<Vec<u8>>> = specs
        .par_iter()
        .enumerate()
        .map(decode_block)
        .collect::<Result<Vec<_>>>()?;

    let mut out: Vec<Vec<u8>> = Vec::with_capacity(strands.len());
    for block_reads in per_block {
        out.extend(block_reads);
    }
    Ok(out)
}

/// Validate the framing and per-block CRC of a v12 fqzcomp quality blob *without*
/// running the (expensive) fqz entropy decode. `verify` uses this to catch
/// corruption in fqz chunks: the per-block CRC covers `record_count ||
/// fqz_blob`, so a flipped payload byte is detected here even though verify never
/// reconstructs records. Mirrors the framing walk in
/// `decompress_quality_fqz_multiblock`; `expected_records` is the chunk's record
/// count, which the per-block `record_count`s must sum to exactly.
pub(super) fn validate_fqz_multiblock_framing(blob: &[u8], expected_records: usize) -> Result<()> {
    if expected_records == 0 {
        // No quality records: an empty (or absent) blob is the only valid state.
        if !blob.is_empty() && blob != [0u8; 4] {
            anyhow::bail!("fqz quality blob non-empty for a 0-record chunk (corrupt archive)");
        }
        return Ok(());
    }
    if blob.len() < 4 {
        anyhow::bail!("fqz quality blob too short ({} bytes)", blob.len());
    }
    let num_blocks = u32::from_le_bytes(blob[0..4].try_into().unwrap()) as usize;
    let mut off = 4usize;
    let mut total_records = 0usize;
    for b in 0..num_blocks {
        if off + 12 > blob.len() {
            anyhow::bail!("fqz quality blob truncated at block {b} header");
        }
        let block_len = u32::from_le_bytes(blob[off..off + 4].try_into().unwrap()) as usize;
        let record_count = u32::from_le_bytes(blob[off + 4..off + 8].try_into().unwrap()) as usize;
        let crc = u32::from_le_bytes(blob[off + 8..off + 12].try_into().unwrap());
        off += 12;
        // Running total can't exceed the chunk's real record count — reject an
        // over-large block before trusting block_len (CRC is attacker-recomputable,
        // so this bound must be independent of it).
        let remaining = expected_records - total_records;
        if record_count > remaining {
            anyhow::bail!(
                "fqz quality block {b}: record_count {record_count} exceeds remaining \
                 records {remaining} (corrupt archive)"
            );
        }
        let end = off
            .checked_add(block_len)
            .filter(|&e| e <= blob.len())
            .ok_or_else(|| anyhow::anyhow!("fqz quality block {b} overflows blob"))?;
        let fqz_blob = &blob[off..end];
        off = end;
        let mut h = Crc::new();
        h.update(&(record_count as u32).to_le_bytes());
        h.update(fqz_blob);
        if h.sum() != crc {
            anyhow::bail!("fqz quality block {b} CRC mismatch — archive corrupt");
        }
        total_records += record_count;
    }
    if off != blob.len() {
        anyhow::bail!(
            "fqz quality blob has {} trailing bytes after {num_blocks} blocks",
            blob.len() - off
        );
    }
    if total_records != expected_records {
        anyhow::bail!(
            "fqz quality blocks total {total_records} records, expected {expected_records} \
             (corrupt archive)"
        );
    }
    Ok(())
}

/// Reject a chunk header whose `num_records` exceeds what the streams can hold.
/// The mapq stream stores exactly one byte per record, so its length is an exact
/// upper bound; this runs before the speculative per-record allocation so an
/// untrusted `num_records` cannot drive an unbounded `Vec::with_capacity`.
fn check_record_count(num_records: usize, mapq_stream_len: usize) -> Result<()> {
    if num_records > mapq_stream_len {
        anyhow::bail!(
            "chunk claims {num_records} records but the mapq stream holds only \
             {mapq_stream_len} bytes (corrupt archive)"
        );
    }
    Ok(())
}

/// Decode one chunk and write its BAM records. `quality_compressed` selects the
/// quality source: `Some(blob)` = the fqzcomp path (quality in a separate blob,
/// stream 13 empty); `None` = the BSC fallback (quality is a per-record stream in
/// `decompressed[13]`). Everything else — field parse, sequence reconstruction,
/// mate-derived TLEN/PNEXT/RNEXT/aux, and record assembly — is shared.
fn decompress_chunk<W: Write>(
    decompressed: &[Vec<u8>],
    quality_compressed: Option<&[u8]>,
    consensus: &ChunkConsensus,
    num_records: usize,
    content_crc32: u32,
    n_ref: i32,
    bam_writer: &mut BamWriter<W>,
) -> Result<()> {
    use rayon::prelude::*;

    // `num_records` is the chunk header's count (u32, untrusted) and drives a
    // speculative `Vec::with_capacity(num_records)` of ~120-byte RecordInfo below.
    // The mapq stream (index 3) holds exactly one byte per record, so a valid
    // chunk satisfies `num_records == mapq.len()`. Reject any header claiming more
    // records than the streams can hold, before allocating (stops a ~30-byte
    // crafted header from forcing a multi-hundred-GB allocation).
    check_record_count(num_records, decompressed.get(3).map_or(0, |s| s.len()))?;

    // Phase 1 (sequential): parse stream offsets and field values.
    let _t_p1 = Instant::now();
    let mut records: Vec<RecordInfo> = Vec::with_capacity(num_records);

    let mut ref_id_off = 0usize;
    let mut pos_off = 0usize;
    let mut mapq_off = 0usize;
    let mut bin_off = 0usize;
    let mut flag_off = 0usize;
    let mut next_ref_id_off = 0usize;
    let mut next_pos_off = 0usize;
    let mut tlen_off = 0usize;
    let mut read_name_off = 0usize;
    let mut cigar_off = 0usize;
    let mut seq_diff_off = 0usize;
    let mut seq_extra_off = 0usize;
    let mut aux_off = 0usize;
    // BSC fallback only: cursor into the per-record raw quality stream (13).
    let mut qual_off = 0usize;
    let parse_qual_stream = quality_compressed.is_none();

    let mut prev_ref_id: i32 = 0;
    let mut prev_pos: i32 = 0;

    for rec_i in 0..num_records {
        let delta_ref = read_i32(&decompressed[1], &mut ref_id_off, "ref_id")
            .with_context(|| format!("record {rec_i}"))?;
        let ref_id = prev_ref_id.wrapping_add(delta_ref);

        let delta_pos = read_i32(&decompressed[2], &mut pos_off, "pos")
            .with_context(|| format!("record {rec_i}"))?;
        let pos = if delta_ref == 0 {
            prev_pos.wrapping_add(delta_pos)
        } else {
            delta_pos
        };
        prev_ref_id = ref_id;
        prev_pos = pos;

        let mapq = read_byte(&decompressed[3], &mut mapq_off, "mapq")
            .with_context(|| format!("record {rec_i}"))?;

        let bin_val = read_u16(&decompressed[4], &mut bin_off, "bin")
            .with_context(|| format!("record {rec_i}"))?;

        let flag = read_u16(&decompressed[5], &mut flag_off, "flag")
            .with_context(|| format!("record {rec_i}"))?;

        // PNEXT/RNEXT (v10): derivable reads carry next_pos/next_ref_id in the
        // bitmap (stream 21), regenerated from the mate in Phase 4; only
        // non-derivable values sit in streams 6/7. For derivable reads we still
        // validate this record's own ref_id/pos here (next is checked in Phase 4
        // once the mate is known).
        let next_derivable = decompressed[21]
            .get(rec_i / 8)
            .is_some_and(|b| (b >> (rec_i % 8)) & 1 != 0);
        let (next_ref_id, next_pos) = if next_derivable {
            validate_alignment_fields(rec_i, ref_id, -1, pos, n_ref)?;
            (0, 0) // placeholder; resolved from the mate in Phase 4
        } else {
            let delta_next_ref = read_i32(&decompressed[6], &mut next_ref_id_off, "next_ref_id")
                .with_context(|| format!("record {rec_i}"))?;
            let next_ref_id = ref_id.wrapping_add(delta_next_ref);
            let delta_next_pos = read_i32(&decompressed[7], &mut next_pos_off, "next_pos")
                .with_context(|| format!("record {rec_i}"))?;
            let next_pos = pos.wrapping_add(delta_next_pos);
            validate_alignment_fields(rec_i, ref_id, next_ref_id, pos, n_ref)?;
            (next_ref_id, next_pos)
        };

        // TLEN: derivable reads carry it in the bitmap (regenerated in Phase 4);
        // only non-derivable values sit in stream 8.
        let tlen_derivable = decompressed[20]
            .get(rec_i / 8)
            .is_some_and(|b| (b >> (rec_i % 8)) & 1 != 0);
        let tlen = if tlen_derivable {
            0
        } else {
            streams::zigzag_decode(
                read_u32(&decompressed[8], &mut tlen_off, "tlen")
                    .with_context(|| format!("record {rec_i}"))?,
            )
        };

        let rn_len = streams::read_varint(&decompressed[9], &mut read_name_off)?;
        let rn_start = read_name_off;
        read_name_off = read_name_off
            .checked_add(rn_len)
            .filter(|&e| e <= decompressed[9].len())
            .ok_or_else(|| anyhow::anyhow!("read_name overflow at record {rec_i}"))?;

        let n_cigar_op = streams::read_varint(&decompressed[10], &mut cigar_off)?;
        let cigar_start = cigar_off;
        let cigar_byte_len = n_cigar_op
            .checked_mul(4)
            .ok_or_else(|| anyhow::anyhow!("cigar overflow at record {rec_i}"))?;
        cigar_off = cigar_off
            .checked_add(cigar_byte_len)
            .filter(|&e| e <= decompressed[10].len())
            .ok_or_else(|| anyhow::anyhow!("cigar stream overflow at record {rec_i}"))?;

        let diff_count = streams::read_varint(&decompressed[11], &mut seq_diff_off)?;
        let diff_start = seq_diff_off;
        let diff_packed_len = diff_count.checked_add(1).map(|n| n / 2).ok_or_else(|| {
            anyhow::anyhow!("diff_count overflow at record {rec_i}: {diff_count}")
        })?;
        seq_diff_off = seq_diff_off
            .checked_add(diff_packed_len)
            .filter(|&e| e <= decompressed[11].len())
            .ok_or_else(|| anyhow::anyhow!("seq_diff overflow at record {rec_i}"))?;

        let extra_count = streams::read_varint(&decompressed[12], &mut seq_extra_off)?;
        let extra_start = seq_extra_off;
        let extra_packed_len = extra_count.checked_add(1).map(|n| n / 2).ok_or_else(|| {
            anyhow::anyhow!("extra_count overflow at record {rec_i}: {extra_count}")
        })?;
        seq_extra_off = seq_extra_off
            .checked_add(extra_packed_len)
            .filter(|&e| e <= decompressed[12].len())
            .ok_or_else(|| anyhow::anyhow!("seq_extra overflow at record {rec_i}"))?;

        let aux_len = streams::read_varint(&decompressed[14], &mut aux_off)?;
        let aux_start = aux_off;
        aux_off = aux_off
            .checked_add(aux_len)
            .filter(|&e| e <= decompressed[14].len())
            .ok_or_else(|| anyhow::anyhow!("aux stream overflow at record {rec_i}"))?;

        let l_seq = diff_count.checked_add(extra_count).ok_or_else(|| {
            anyhow::anyhow!(
                "l_seq overflow at record {rec_i}: diff={diff_count} extra={extra_count}"
            )
        })?;

        // BSC fallback: quality is a per-record [varint l_seq][l_seq bytes] stream
        // in stream 13. Skip the varint (== l_seq by construction) and record the
        // byte offset; the fqzcomp path leaves stream 13 empty (qual_start unused).
        let qual_start = if parse_qual_stream {
            let _stored_l_seq = streams::read_varint(&decompressed[13], &mut qual_off)?;
            let start = qual_off;
            qual_off = qual_off
                .checked_add(l_seq)
                .filter(|&e| e <= decompressed[13].len())
                .ok_or_else(|| anyhow::anyhow!("quality overflow at record {rec_i}"))?;
            start
        } else {
            0
        };

        records.push(RecordInfo {
            ref_id,
            pos,
            mapq,
            bin_val,
            flag,
            next_ref_id,
            next_pos,
            tlen,
            tlen_derivable,
            next_derivable,
            rn_start,
            rn_len,
            cigar_start,
            n_cigar_op,
            diff_start,
            diff_packed_len,
            diff_count,
            extra_start,
            extra_packed_len,
            extra_count,
            aux_start,
            aux_len,
            l_seq,
            qual_start,
        });
    }

    // Phase 2 (parallel): reconstruct sequences. We produce two forms per record:
    //   - `packed`: the BAM 4-bit packed bytes built directly from the recovered
    //     nibbles, which preserves ALL 16 BAM codes (incl. IUPAC ambiguity and
    //     '=') byte-exact. This is what gets written to the output BAM.
    //   - `ascii`: A/C/G/T/N text (ambiguity codes fold to N) used ONLY for MD/NM
    //     derivation, which is consistent with the compress side (it derives MD/NM
    //     from the same folded ASCII), so non-ACGTN reads simply keep their tags
    //     verbatim. Never use `ascii` for the output SEQ — that loses ambiguity.
    if *DECODE_TIMING {
        T_P1.fetch_add(_t_p1.elapsed().as_micros() as u64, Ordering::Relaxed);
    }

    let _t_p2 = Instant::now();
    let seq_forms: Vec<(Vec<u8>, Vec<u8>)> = records
        .par_iter()
        .map(|rec| {
            let diff_nibbles = streams::unpack_nibble_pairs(
                &decompressed[11][rec.diff_start..rec.diff_start + rec.diff_packed_len],
                rec.diff_count,
            );
            let extra_nibbles = streams::unpack_nibble_pairs(
                &decompressed[12][rec.extra_start..rec.extra_start + rec.extra_packed_len],
                rec.extra_count,
            );
            let seq_nibbles = streams::reconstruct_sequence_nibbles(
                consensus,
                rec.ref_id,
                rec.pos,
                &decompressed[10][rec.cigar_start..rec.cigar_start + rec.n_cigar_op * 4],
                rec.n_cigar_op,
                &diff_nibbles,
                &extra_nibbles,
                rec.l_seq,
            )?;
            let packed = crate::io::bam::pack_nibbles(&seq_nibbles);
            let ascii = streams::nibbles_to_ascii(&seq_nibbles);
            Ok((packed, ascii))
        })
        .collect::<Result<Vec<(Vec<u8>, Vec<u8>)>>>()?;
    if *DECODE_TIMING {
        T_P2.fetch_add(_t_p2.elapsed().as_micros() as u64, Ordering::Relaxed);
    }

    // Phase 3: quality. fqzcomp path (v12) decodes the separate blob, de-reorienting
    // reverse-strand reads via the FLAG stream → ASCII Phred+33 per record. BSC
    // fallback path leaves this None and reads raw quality from stream 13 in Phase 4.
    let _t_p3 = Instant::now();
    let fqz_quality: Option<Vec<Vec<u8>>> = match quality_compressed {
        Some(qc) => {
            let strands: Vec<bool> = records.iter().map(|r| r.flag & 0x10 != 0).collect();
            Some(decompress_quality_fqz_multiblock(qc, &strands)?)
        }
        None => None,
    };
    if *DECODE_TIMING {
        T_P3.fetch_add(_t_p3.elapsed().as_micros() as u64, Ordering::Relaxed);
    }

    // MD (15-17) + MC (18-19) + NM (22-24) derivation streams.
    let _t_p4 = Instant::now();
    let md_present_bits = &decompressed[15];
    let md_index_stream = &decompressed[16];
    let md_exceptions = md_codec::deserialize_exceptions(&decompressed[17])?;
    let mc_present_bits = &decompressed[18];
    let mc_index_stream = &decompressed[19];
    let nm_present_bits = &decompressed[22];
    // The MD/MC/NM present bitmaps are indexed `bits[rec_idx / 8]` in the Phase 4
    // loop; validate they hold one bit per record up front so a truncated bitmap
    // in a corrupt archive is a clean error rather than an index panic in a
    // decode worker. (Stream 15/18/22 lengths are otherwise unconstrained.)
    {
        let need = records.len().div_ceil(8);
        for (name, bits) in [
            ("md", md_present_bits),
            ("mc", mc_present_bits),
            ("nm", nm_present_bits),
        ] {
            if bits.len() < need {
                anyhow::bail!(
                    "{name}_present bitmap too short: {} bytes, need {need} for {} records (corrupt archive)",
                    bits.len(),
                    records.len()
                );
            }
        }
    }
    let nm_index_stream = &decompressed[23];
    let nm_type_stream = &decompressed[24];
    let mut md_index_off = 0usize;
    let mut mc_index_off = 0usize;
    let mut nm_index_off = 0usize;
    let mut nm_type_off = 0usize;

    // Rebuild the mate map (deterministic from names + flags) for MC/TLEN derivation.
    let names: Vec<&[u8]> = records
        .iter()
        .map(|r| &decompressed[9][r.rn_start..r.rn_start + r.rn_len])
        .collect();
    let flags: Vec<u16> = records.iter().map(|r| r.flag).collect();
    let mate = md_codec::build_mate_map(&names, &flags);
    let cigar_of = |r: &RecordInfo| -> &[u8] {
        &decompressed[10][r.cigar_start..r.cigar_start + r.n_cigar_op * 4]
    };

    // Phase 4 (sequential): assemble records and write.
    let mut content_crc = Crc::new();
    let mut data_buf = Vec::with_capacity(512);
    for (rec_idx, rec) in records.iter().enumerate() {
        let seq_packed = &seq_forms[rec_idx].0;

        let read_name = &decompressed[9][rec.rn_start..rec.rn_start + rec.rn_len];
        let cigar_bytes = cigar_of(rec);
        let stored_aux = &decompressed[14][rec.aux_start..rec.aux_start + rec.aux_len];
        let mate_idx = mate[rec_idx];
        let mate_cigar = (mate_idx != u32::MAX).then(|| cigar_of(&records[mate_idx as usize]));

        // TLEN: derive from the mate when the bit says so.
        let tlen = if rec.tlen_derivable {
            if mate_idx == u32::MAX {
                anyhow::bail!(
                    "record {rec_idx}: tlen_derivable set but mate not in chunk (corrupt archive)"
                );
            }
            let m = &records[mate_idx as usize];
            md_codec::derive_tlen(
                rec.pos,
                &parse_cigar_ops(cigar_bytes),
                m.pos,
                &parse_cigar_ops(cigar_of(m)),
                rec.flag & 0x40 != 0,
            )
        } else {
            rec.tlen
        };

        // PNEXT/RNEXT: derive from the mate when the bit says so (next_pos ==
        // mate.pos, next_ref_id == mate.ref_id, exact per SAM spec).
        let (next_ref_id, next_pos) = if rec.next_derivable {
            if mate_idx == u32::MAX {
                anyhow::bail!(
                    "record {rec_idx}: next_derivable set but mate not in chunk (corrupt archive)"
                );
            }
            let m = &records[mate_idx as usize];
            validate_alignment_fields(rec_idx, rec.ref_id, m.ref_id, rec.pos, n_ref)?;
            (m.ref_id, m.pos)
        } else {
            (rec.next_ref_id, rec.next_pos)
        };

        let md_present = (md_present_bits[rec_idx / 8] >> (rec_idx % 8)) & 1 != 0;
        let mc_present = (mc_present_bits[rec_idx / 8] >> (rec_idx % 8)) & 1 != 0;
        let nm_present = (nm_present_bits[rec_idx / 8] >> (rec_idx % 8)) & 1 != 0;
        let regen = regen_aux(
            stored_aux,
            md_present,
            md_index_stream,
            &mut md_index_off,
            mc_present,
            mc_index_stream,
            &mut mc_index_off,
            nm_present,
            nm_index_stream,
            &mut nm_index_off,
            nm_type_stream,
            &mut nm_type_off,
            mate_cigar,
            rec.ref_id,
            rec.pos,
            cigar_bytes,
            &seq_forms[rec_idx].1,
            consensus,
            &md_exceptions,
        )?;
        let aux_bytes: &[u8] = regen.as_deref().unwrap_or(stored_aux);

        let record_len = 32
            + read_name.len()
            + cigar_bytes.len()
            + seq_packed.len()
            + rec.l_seq
            + aux_bytes.len();
        data_buf.clear();
        data_buf.reserve(record_len);

        // BAM format limits: QNAME ≤ 254 bytes, n_cigar_op fits in u16, l_seq fits in i32
        if read_name.len() > 254 {
            anyhow::bail!(
                "record {rec_idx}: read name too long ({} bytes, BAM max 254)",
                read_name.len()
            );
        }
        if rec.n_cigar_op > u16::MAX as usize {
            anyhow::bail!(
                "record {rec_idx}: n_cigar_op {} exceeds BAM u16 limit",
                rec.n_cigar_op
            );
        }
        if rec.l_seq > i32::MAX as usize {
            anyhow::bail!(
                "record {rec_idx}: l_seq {} exceeds BAM i32 limit",
                rec.l_seq
            );
        }
        data_buf.extend_from_slice(&rec.ref_id.to_le_bytes());
        data_buf.extend_from_slice(&rec.pos.to_le_bytes());
        data_buf.push(read_name.len() as u8);
        data_buf.push(rec.mapq);
        data_buf.extend_from_slice(&rec.bin_val.to_le_bytes());
        data_buf.extend_from_slice(&(rec.n_cigar_op as u16).to_le_bytes());
        data_buf.extend_from_slice(&rec.flag.to_le_bytes());
        data_buf.extend_from_slice(&(rec.l_seq as i32).to_le_bytes());
        data_buf.extend_from_slice(&next_ref_id.to_le_bytes());
        data_buf.extend_from_slice(&next_pos.to_le_bytes());
        data_buf.extend_from_slice(&tlen.to_le_bytes());
        data_buf.extend_from_slice(read_name);
        data_buf.extend_from_slice(cigar_bytes);
        data_buf.extend_from_slice(seq_packed);

        // Quality: fqzcomp path holds ASCII Phred+33 (convert back to BAM Phred);
        // BSC fallback holds raw BAM quality bytes in stream 13 (copy directly).
        match &fqz_quality {
            Some(qa) => {
                let qual = &qa[rec_idx];
                // BAM stores exactly l_seq quality bytes. The record header already
                // wrote l_seq and the sequence is sized from it, so a decoded
                // quality of a different length would emit a structurally invalid
                // BAM record. Reject the mismatch explicitly rather than relying on
                // the downstream content-CRC to notice.
                if qual.len() != rec.l_seq {
                    anyhow::bail!(
                        "record {rec_idx}: decoded quality length {} != l_seq {} (corrupt archive)",
                        qual.len(),
                        rec.l_seq
                    );
                }
                let dst = data_buf.len();
                data_buf.resize(dst + qual.len(), 0);
                streams::ascii_qual_to_bam(qual, &mut data_buf[dst..dst + qual.len()]);
            }
            None => {
                data_buf.extend_from_slice(
                    &decompressed[13][rec.qual_start..rec.qual_start + rec.l_seq],
                );
            }
        }

        data_buf.extend_from_slice(aux_bytes);

        content_crc.update(&data_buf);
        bam_writer.write_record(&data_buf)?;
    }
    if *DECODE_TIMING {
        T_P4.fetch_add(_t_p4.elapsed().as_micros() as u64, Ordering::Relaxed);
    }

    verify_content_crc(content_crc.sum(), content_crc32)?;
    Ok(())
}

/// Verify that the records reconstructed for a chunk match the fidelity hash
/// captured over the original records at compress time. A mismatch means the
/// chunk did not round-trip (an encode/decode asymmetry or undetected
/// corruption), so the output BAM cannot be trusted.
fn verify_content_crc(computed: u32, expected: u32) -> Result<()> {
    if computed != expected {
        anyhow::bail!(
            "Chunk content CRC32 mismatch: stored {:08x}, reconstructed {:08x} \
             — decompressed records do not match the original (archive corrupt \
             or a codec bug)",
            expected,
            computed
        );
    }
    Ok(())
}

/// Per-record decode working memory used to bound the auto-selected decompress
/// window (decompressed columnar streams + reconstructed BAM record bytes +
/// fqz working set). Conservative; only used for window sizing.
const PER_RECORD_DECODE_BYTES: u64 = 1400;
/// Fixed decompress overhead outside the per-chunk pipeline (archive reader
/// buffer, BGZF writer workers, allocator retention).
const DECODE_BASE_OVERHEAD_BYTES: u64 = 2 * 1024 * 1024 * 1024;
/// Safety ceiling on the auto-selected decompress window.
const DECODE_WINDOW_CAP: usize = 16;
/// Fallback decompress window when available RAM cannot be read.
const DECODE_WINDOW_FALLBACK: usize = 4;

/// Resolve the decompress pipeline window: how many chunks to decode
/// concurrently. Like the compress side, this is bounded by core count and by
/// keeping estimated peak RSS under ~70% of available RAM; the chunk size is
/// estimated from the archive's record/chunk totals. The window only caps
/// concurrency, so over-provisioning is harmless and it never affects output.
fn resolve_decompress_window(num_records: u64, num_chunks: u32) -> usize {
    if num_chunks <= 1 {
        return 1;
    }
    let est_chunk_records = (num_records / num_chunks as u64).max(1);
    let per_chunk = est_chunk_records
        .saturating_mul(PER_RECORD_DECODE_BYTES)
        .max(1);
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8);
    let mem_window = match crate::compression::compress::available_memory_bytes() {
        Some(avail) => {
            let budget = avail / 10 * 7; // ~70% of available RAM
            let usable = budget.saturating_sub(DECODE_BASE_OVERHEAD_BYTES);
            (usable / per_chunk).max(1) as usize
        }
        None => DECODE_WINDOW_FALLBACK,
    };
    mem_window
        .min(cores)
        .min(num_chunks as usize)
        .clamp(1, DECODE_WINDOW_CAP)
}

/// Decode one chunk into its reconstructed, length-prefixed BAM record bytes.
///
/// Runs off the critical path (one of `window` worker threads): decompresses the
/// chunk's streams, reconstructs every record into an in-memory
/// `BamWriter<Vec<u8>>`, and recomputes the per-chunk content CRC (the lossless
/// fidelity check, which errors on mismatch). The returned bytes are written to
/// the real BGZF output in chunk order by the caller. Each chunk is fully
/// self-contained (within-chunk mate derivation), so chunks decode independently.
fn decode_chunk_to_bytes(
    compressed_streams: &[Vec<u8>],
    chunk_flags: u8,
    num_records: usize,
    content_crc32: u32,
    n_ref: i32,
) -> Result<Vec<u8>> {
    let uses_fqz_quality = (chunk_flags & CHUNK_FLAG_FQZ_QUALITY) != 0;
    let decompressed = timed(&T_BSC, || {
        decompress_streams(compressed_streams, uses_fqz_quality, num_records)
    })?;
    let consensus = ChunkConsensus::deserialize(&decompressed[0])?;

    let mut buf_writer = BamWriter::new(Vec::<u8>::new());
    // fqzcomp path when the chunk used fqz (quality in stream 13's blob);
    // otherwise the BSC fallback reads quality from the decompressed stream 13.
    let quality_compressed = if uses_fqz_quality && num_records > 0 {
        Some(compressed_streams[13].as_slice())
    } else {
        None
    };
    decompress_chunk(
        &decompressed,
        quality_compressed,
        &consensus,
        num_records,
        content_crc32,
        n_ref,
        &mut buf_writer,
    )?;
    Ok(buf_writer.into_inner())
}

/// On-disk size of a `ChunkHeader`: num_records(4) + chunk_flags(1) + crc32(4) +
/// content_crc32(4) + stream_sizes[NUM_STREAMS]·4. A positive lower bound on a
/// chunk's disk footprint, used to cap the layout-scan reservation.
const CHUNK_HEADER_BYTES: u64 = 13 + (NUM_STREAMS as u64) * 4;

/// Parse the archive header + SAM/BAM header block from a reader positioned at the
/// start of the archive. Returns the header, the raw BAM header text, the reference
/// dictionary bytes, and the reference count, leaving the reader at the first chunk
/// header. Shared by `decompress` and `decode_chunk_range` (every NUMA worker
/// reconstructs the BAM header from the archive header).
fn read_header_block<R: Read>(reader: &mut R) -> Result<(ArchiveHeader, Vec<u8>, Vec<u8>, i32)> {
    let header = ArchiveHeader::read_from(reader)?;
    let header_payload = decompress_sam_header(&header.sam_header_compressed)?;
    if header_payload.len() < 4 {
        anyhow::bail!("SAM header payload too short ({} bytes)", header_payload.len());
    }
    let header_raw_len = u32::from_le_bytes(header_payload[0..4].try_into().unwrap()) as usize;
    if 4 + header_raw_len > header_payload.len() {
        anyhow::bail!(
            "SAM header length {} exceeds payload size {}",
            header_raw_len,
            header_payload.len() - 4
        );
    }
    let header_raw = header_payload[4..4 + header_raw_len].to_vec();
    let ref_dict_bytes = header_payload[4 + header_raw_len..].to_vec();
    if ref_dict_bytes.len() < 4 {
        anyhow::bail!(
            "reference dictionary too short ({} bytes)",
            ref_dict_bytes.len()
        );
    }
    let n_ref = i32::from_le_bytes(ref_dict_bytes[0..4].try_into().unwrap());
    if n_ref < 0 {
        anyhow::bail!("invalid reference count {n_ref} in archive");
    }
    Ok((header, header_raw, ref_dict_bytes, n_ref))
}

/// Windowed decode pipeline shared by full `decompress` and NUMA `decode_chunk_range`.
/// The reader must be positioned at the first chunk header to decode; this reads
/// `chunk_count` consecutive chunks, dispatching each chunk's reconstruction to a
/// worker thread (up to `decode_window` concurrent) and draining finished chunks
/// into `bam_writer` in archive order. Returns the number of records written. The
/// caller is responsible for `finish()`ing the writer.
fn run_decode_pipeline<W: Write>(
    reader: &mut BufReader<std::fs::File>,
    bam_writer: &mut BamWriter<W>,
    chunk_count: u32,
    decode_window: usize,
    n_ref: i32,
) -> Result<u64> {
    type DecodeHandle = std::thread::JoinHandle<Result<Vec<u8>>>;
    let mut pending: VecDeque<(usize, DecodeHandle)> = VecDeque::new();
    let mut dispatch_idx: u32 = 0;
    let mut total_written: u64 = 0;

    let mut read_and_dispatch = |reader: &mut BufReader<std::fs::File>,
                                 pending: &mut VecDeque<(usize, DecodeHandle)>|
     -> Result<bool> {
        if dispatch_idx >= chunk_count {
            return Ok(false);
        }
        let chunk_header = ChunkHeader::read_from(reader)?;
        let num_records = chunk_header.num_records as usize;
        let compressed_streams = read_chunk_streams(reader, &chunk_header)?;
        let flags = chunk_header.chunk_flags;
        let crc = chunk_header.content_crc32;
        dispatch_idx += 1;
        let handle = std::thread::spawn(move || {
            decode_chunk_to_bytes(&compressed_streams, flags, num_records, crc, n_ref)
        });
        pending.push_back((num_records, handle));
        Ok(true)
    };

    let decode_result: Result<()> = (|| {
        for _ in 0..decode_window {
            if !read_and_dispatch(reader, &mut pending)? {
                break;
            }
        }
        while let Some((num_records, handle)) = pending.pop_front() {
            let bytes = handle
                .join()
                .unwrap_or_else(|_| Err(anyhow::anyhow!("chunk decode thread panicked")))?;
            bam_writer.write_block_bytes(&bytes)?;
            total_written += num_records as u64;
            read_and_dispatch(reader, &mut pending)?;
        }
        Ok(())
    })();

    // Always join in-flight chunk-decode threads before propagating an error, so no
    // detached decode threads keep running after a failure.
    for (_, h) in pending.drain(..) {
        let _ = h.join();
    }
    decode_result?;
    Ok(total_written)
}

/// Chunk layout of a BZ archive — the byte offset and read count of every chunk.
/// Built by a cheap forward pre-scan that reads each `ChunkHeader` and seeks past
/// its compressed data (sized by `sum(stream_sizes)`) without decompressing, so it
/// is O(num_chunks) header reads + seeks. bz's analogue of qz's `read_chunk_layout`,
/// the basis for NUMA range sharding.
#[derive(Clone, Debug)]
pub struct BzChunkLayout {
    pub num_chunks: u32,
    pub num_records: u64,
    /// Records in each chunk, in archive order (`len == num_chunks`).
    pub per_chunk_reads: Vec<u64>,
    /// Byte offset of each chunk's `ChunkHeader`, in archive order (`len == num_chunks`).
    pub chunk_offsets: Vec<u64>,
}

/// Pre-scan a BZ archive's chunk layout (see `BzChunkLayout`). Reads the archive
/// header, then walks every chunk header recording its offset and skipping its data.
/// Decompresses nothing.
pub fn read_bz_chunk_layout(archive: &std::path::Path) -> Result<BzChunkLayout> {
    let file = std::fs::File::open(archive)?;
    let file_len = file.metadata()?.len();
    let mut reader = BufReader::with_capacity(64 * 1024, file);
    let header = ArchiveHeader::read_from(&mut reader)?;
    let num_chunks = header.num_chunks;

    // `num_chunks` is untrusted; each chunk occupies at least a fixed-size header on
    // disk, so cap the reservation by what the file could actually hold.
    let cap = (num_chunks as u64).min(file_len / CHUNK_HEADER_BYTES + 1) as usize;
    let mut per_chunk_reads = Vec::with_capacity(cap);
    let mut chunk_offsets = Vec::with_capacity(cap);

    let mut offset = reader.stream_position()?;
    for c in 0..num_chunks {
        chunk_offsets.push(offset);
        let ch = ChunkHeader::read_from(&mut reader)
            .with_context(|| format!("reading chunk {c} header during layout scan"))?;
        per_chunk_reads.push(ch.num_records as u64);
        let data_len: u64 = ch.stream_sizes.iter().map(|&s| s as u64).sum();
        let data_start = reader.stream_position()?;
        let next = data_start
            .checked_add(data_len)
            .filter(|&n| n <= file_len)
            .ok_or_else(|| anyhow::anyhow!("chunk {c} data ({data_len} bytes) overruns archive"))?;
        reader.seek(SeekFrom::Start(next))?;
        offset = next;
    }
    Ok(BzChunkLayout {
        num_chunks,
        num_records: header.num_records,
        per_chunk_reads,
        chunk_offsets,
    })
}

/// Bound on one NUMA decode **worker**'s peak RSS, for the gate's per-node memory
/// check. Peak decode RSS is the number of chunks held resident at once (the decode
/// window) × per-chunk resident bytes (measured ~1.4 KB/record on HG002 chr20,
/// rounded up). The key subtlety: a *worker* decodes only its shard — at most ~half
/// the chunks for the smallest useful 2-way split — so the per-worker window is
/// `min(full_window, ceil(num_chunks/2))`, NOT the whole-archive window (which would
/// over-budget by ~2× and wrongly disqualify nodes). More workers hold fewer chunks
/// each, so the 2-way figure is the conservative per-worker peak. The worker also
/// self-limits its window by available RAM, so erring slightly low only risks a
/// slower (correct) self-limited decode, never an OOM.
pub fn decode_peak_rss_bound(num_records: u64, num_chunks: u32) -> u64 {
    const RESIDENT_BYTES_PER_RECORD: u64 = 1_600;
    const BASE_BYTES: u64 = 1 << 30; // 1 GiB
    let full_window = resolve_decompress_window(num_records, num_chunks) as u64;
    let per_worker_chunks = (num_chunks as u64).div_ceil(2).max(1);
    let window = full_window.min(per_worker_chunks);
    let records_per_chunk = num_records / (num_chunks.max(1) as u64);
    BASE_BYTES + window * records_per_chunk.max(1) * RESIDENT_BYTES_PER_RECORD
}

/// Decode chunks `[chunk_start, chunk_end)` of a BZ archive into `out_part` as a
/// self-contained BGZF/BAM stream, for NUMA decode sharding. The first range
/// (`chunk_start == 0`) writes the BAM header; later ranges write records only, so
/// the parts concatenate — after stripping each non-final part's trailing 28-byte
/// BGZF EOF marker — into one valid BAM. Always finishes the BGZF stream (writes its
/// EOF marker). Returns the records written, after verifying the count matches the
/// layout's per-chunk totals (a per-range integrity check on top of the per-chunk
/// content CRC the decode itself enforces).
pub fn decode_chunk_range(
    archive: &std::path::Path,
    chunk_start: u32,
    chunk_end: u32,
    threads: usize,
    out_part: &std::path::Path,
) -> Result<u64> {
    let layout = read_bz_chunk_layout(archive)?;
    if chunk_start >= chunk_end || chunk_end > layout.num_chunks {
        anyhow::bail!(
            "invalid chunk range [{chunk_start}, {chunk_end}) for {} chunks",
            layout.num_chunks
        );
    }

    let file = std::fs::File::open(archive)?;
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);
    let (header, header_raw, ref_dict_bytes, n_ref) = read_header_block(&mut reader)?;

    let bgzf_workers = NonZeroUsize::new(threads.clamp(1, 16)).unwrap();
    let mut bam_writer = BamWriter::from_path_mt(out_part, bgzf_workers)?;
    if chunk_start == 0 {
        bam_writer.write_header(&header_raw, &ref_dict_bytes)?;
    }

    reader.seek(SeekFrom::Start(layout.chunk_offsets[chunk_start as usize]))?;
    let decode_window = resolve_decompress_window(header.num_records, header.num_chunks);
    let total = run_decode_pipeline(
        &mut reader,
        &mut bam_writer,
        chunk_end - chunk_start,
        decode_window,
        n_ref,
    )?;

    let expected: u64 = layout.per_chunk_reads[chunk_start as usize..chunk_end as usize]
        .iter()
        .sum();
    if total != expected {
        anyhow::bail!(
            "chunk range [{chunk_start}, {chunk_end}) decoded {total} records, expected {expected} (corrupt archive)"
        );
    }
    bam_writer.finish()?;
    Ok(total)
}

/// Decompress a BZ archive back to a BAM file.
pub fn decompress(config: &DecompressConfig) -> Result<()> {
    let start_time = Instant::now();

    let file = std::fs::File::open(&config.input)?;
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);

    let (header, header_raw, ref_dict_bytes, n_ref) = read_header_block(&mut reader)?;
    info!(
        "Archive: {} records in {} chunks",
        header.num_records, header.num_chunks,
    );

    // Open BAM output with multi-threaded BGZF writer. The output BGZF deflate
    // is ~24% of decompress CPU (profiled) and the drain feeds it synchronously,
    // so too few writer threads throttle the whole decode pipeline. Cap at 16
    // (4 → 16 measured 76s → 45s on HG002 chr20, byte-identical); beyond ~16 the
    // fqz/BSC decode dominates and more deflate threads stop helping.
    let bgzf_workers = std::thread::available_parallelism()
        .unwrap_or(NonZeroUsize::new(1).unwrap())
        .min(NonZeroUsize::new(16).unwrap());

    // Write the BAM to a randomized sibling temp file created with O_EXCL, then
    // rename atomically on success. The decode does its integrity checks (per-chunk
    // content CRC, record-count total) as it streams, so a corrupt archive would
    // otherwise leave a truncated BAM — or, with --force, destroy a valid existing
    // output. The random name + O_EXCL also avoid collisions between concurrent
    // jobs and refuse to follow a pre-planted symlink. We open the temp file
    // ourselves and hand the descriptor to the writer (no path re-open).
    let (temp_file, temp_path) = tempfile::Builder::new()
        .prefix(".bz-")
        .suffix(".tmp")
        .tempfile_in(crate::compression::archive::output_parent(&config.output))?
        .into_parts();
    let mut bam_writer = BamWriter::from_file_mt(temp_file, bgzf_workers)?;
    bam_writer.write_header(&header_raw, &ref_dict_bytes)?;

    // Windowed decode pipeline (see `run_decode_pipeline`). Reader is positioned at
    // the first chunk header; peak RSS is bounded by `decode_window × chunk size`,
    // constant in archive size — true bounded-memory streaming.
    let decode_window = resolve_decompress_window(header.num_records, header.num_chunks);
    if header.num_chunks > 1 {
        info!(
            "Parallel decode pipeline: window={} chunks decoding simultaneously",
            decode_window
        );
    }
    let total_written = run_decode_pipeline(
        &mut reader,
        &mut bam_writer,
        header.num_chunks,
        decode_window,
        n_ref,
    )?;

    // The per-chunk record counts must add up to the archive header's total.
    // A mismatch means the header (or a chunk header) is corrupt.
    if total_written != header.num_records {
        anyhow::bail!(
            "record count mismatch: archive header declares {} records but chunks total {} — archive is corrupted",
            header.num_records,
            total_written
        );
    }

    bam_writer.finish()?;

    // Commit: atomically rename the finished BAM over the target. `persist`
    // consumes the `TempPath` (its Drop no longer removes the file); on failure
    // the error carries the path and cleanup still happens when it drops.
    temp_path
        .persist(&config.output)
        .map_err(|e| anyhow::anyhow!("Failed to finalize output {:?}: {}", config.output, e))?;

    let elapsed = start_time.elapsed().as_secs_f64();
    info!("Decompressed {} records in {:.2}s", total_written, elapsed);

    if *DECODE_TIMING {
        let us = |a: &AtomicU64| a.load(Ordering::Relaxed) as f64 / 1e6;
        let (bsc, p1, p2, p3, p4) = (
            us(&T_BSC),
            us(&T_P1),
            us(&T_P2),
            us(&T_P3),
            us(&T_P4),
        );
        let sum = bsc + p1 + p2 + p3 + p4;
        eprintln!(
            "\n=== BZ_DECODE_TIMING (cumulative CPU-s across chunk-decode threads) ===\n\
             wall {:.2}s   summed-phase {:.2}s\n\
             {:<26} {:>8} {:>8} {:>9}\n\
             {}",
            elapsed,
            sum,
            "phase",
            "CPU-s",
            "%phase",
            "concurrency",
            [
                ("BSC stream decode", bsc),
                ("P1 parse (serial)", p1),
                ("P2 reconstruct (par_iter)", p2),
                ("P3 fqz quality (serial)", p3),
                ("P4 assemble+write (serial)", p4),
            ]
            .iter()
            .map(|(name, t)| format!(
                "{:<26} {:>8.2} {:>7.1}% {:>9.1}x",
                name,
                t,
                if sum > 0.0 { 100.0 * t / sum } else { 0.0 },
                if elapsed > 0.0 { t / elapsed } else { 0.0 }
            ))
            .collect::<Vec<_>>()
            .join("\n")
        );
    }

    Ok(())
}

/// Per-record intermediate data extracted in single pass.
struct RecordInfo {
    ref_id: i32,
    pos: i32,
    mapq: u8,
    bin_val: u16,
    flag: u16,
    next_ref_id: i32,
    next_pos: i32,
    tlen: i32,
    tlen_derivable: bool,
    next_derivable: bool,
    rn_start: usize,
    rn_len: usize,
    cigar_start: usize,
    n_cigar_op: usize,
    diff_start: usize,
    diff_packed_len: usize,
    diff_count: usize,
    extra_start: usize,
    extra_packed_len: usize,
    extra_count: usize,
    aux_start: usize,
    aux_len: usize,
    l_seq: usize,
    /// Offset of this record's raw quality in stream 13 — only set on the BSC
    /// fallback path (where quality lives in stream 13); 0 on the fqzcomp path
    /// (quality comes from the separate fqz blob).
    qual_start: usize,
}

/// Upper bound on a single record's contribution to one decompressed stream.
/// Generous (covers seq+quality for reads up to ~64 kb and large aux blobs).
const MAX_DECOMPRESSED_BYTES_PER_RECORD: usize = 64 * 1024;
/// Absolute per-stream ceiling, so a header claiming billions of records can't
/// turn the per-record budget into an unbounded allocation request.
const MAX_DECOMPRESSED_BYTES_PER_STREAM: usize = 8 * 1024 * 1024 * 1024; // 8 GiB
/// Floor so tiny-record chunks with a sizeable header/aux stream still decode.
const MIN_STREAM_DECOMPRESS_CAP: usize = 64 * 1024 * 1024; // 64 MiB
/// Absolute cap on a decompressed SAM header (decoded with no record context).
const MAX_SAM_HEADER_BYTES: usize = 512 * 1024 * 1024; // 512 MiB

/// Per-stream cumulative-output cap for untrusted input. `bsc::decompress_parallel`
/// is unbounded in cumulative output, so an archive can stack many high-expansion
/// BSC blocks and exhaust memory before record-count validation runs. Genomic
/// streams legitimately expand by large ratios (constant quality → ~1000×), so a
/// ratio bound would reject valid data; bound by records × a generous per-record
/// budget, clamped by an absolute ceiling. This is a memory guard, not a
/// correctness check — the per-chunk content CRC proves fidelity.
/// `decompress_parallel_max` allocates the actual output and bails on exceed, so a
/// loose cap costs nothing on legitimate data.
fn stream_decompress_cap(num_records: usize) -> usize {
    // clamp(floor, ceil): MIN_STREAM_DECOMPRESS_CAP <= MAX_DECOMPRESSED_BYTES_PER_STREAM
    // by construction, so this never panics.
    num_records
        .saturating_mul(MAX_DECOMPRESSED_BYTES_PER_RECORD)
        .clamp(MIN_STREAM_DECOMPRESS_CAP, MAX_DECOMPRESSED_BYTES_PER_STREAM)
}

/// Decompress a single BSC-compressed stream, bailing if the decompressed output
/// would exceed `max_total` bytes (DoS guard for untrusted input).
fn decompress_one(data: &[u8], max_total: usize) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }
    bsc::decompress_parallel_max(data, max_total)
}

/// Decompress the archive's SAM header with an absolute output cap (it carries no
/// record context). Shared by `decompress` and `verify` so both bound this
/// untrusted-input decode identically.
pub(super) fn decompress_sam_header(data: &[u8]) -> Result<Vec<u8>> {
    bsc::decompress_parallel_max(data, MAX_SAM_HEADER_BYTES)
        .context("Failed to decompress SAM header")
}

/// Decompress all 15 streams in parallel groups. Every stream is BSC.
/// `num_records` (the chunk header's count) bounds each stream's decompressed
/// output against memory-exhausting untrusted input.
pub(super) fn decompress_streams(
    compressed: &[Vec<u8>],
    use_fqz_quality: bool,
    num_records: usize,
) -> Result<Vec<Vec<u8>>> {
    let cap = stream_decompress_cap(num_records);
    let ((group_a, group_b), (group_c, group_d)) = rayon::join(
        || {
            rayon::join(
                || -> Result<Vec<Vec<u8>>> {
                    // consensus (index 0) + alignment (1-5)
                    let mut results = Vec::with_capacity(6);
                    for data in &compressed[0..6] {
                        results.push(decompress_one(data, cap)?);
                    }
                    Ok(results)
                },
                || -> Result<Vec<Vec<u8>>> {
                    // alignment streams 6-8
                    let mut results = Vec::with_capacity(3);
                    for data in &compressed[6..9] {
                        results.push(decompress_one(data, cap)?);
                    }
                    Ok(results)
                },
            )
        },
        || {
            rayon::join(
                || -> Result<Vec<Vec<u8>>> {
                    // read data streams 9-12
                    let mut results = Vec::with_capacity(4);
                    for data in &compressed[9..13] {
                        results.push(decompress_one(data, cap)?);
                    }
                    Ok(results)
                },
                || -> Result<Vec<Vec<u8>>> {
                    if use_fqz_quality {
                        // fqz: pass through raw, decompress aux
                        let qual_passthrough = Vec::new();
                        let aux = decompress_one(&compressed[14], cap)?;
                        Ok(vec![qual_passthrough, aux])
                    } else {
                        let qual = decompress_one(&compressed[13], cap)?;
                        let aux = decompress_one(&compressed[14], cap)?;
                        Ok(vec![qual, aux])
                    }
                },
            )
        },
    );

    // Derivation streams 15-24 (always BSC): MD (15-17), MC (18-19),
    // tlen_derivable (20), next_derivable (21), NM (22-24).
    let md_streams = || -> Result<Vec<Vec<u8>>> {
        let mut results = Vec::with_capacity(10);
        for data in &compressed[15..25] {
            results.push(decompress_one(data, cap)?);
        }
        Ok(results)
    }();

    let mut result = Vec::with_capacity(NUM_STREAMS);
    for v in group_a? {
        result.push(v);
    }
    for v in group_b? {
        result.push(v);
    }
    for v in group_c? {
        result.push(v);
    }
    for v in group_d? {
        result.push(v);
    }
    for v in md_streams? {
        result.push(v);
    }

    Ok(result)
}

#[cfg(test)]
mod alignment_field_tests {
    use super::*;

    #[test]
    fn validate_alignment_fields_accepts_valid_and_rejects_out_of_range() {
        // n_ref = 3: valid ref ids are -1 (unmapped) .. 2.
        assert!(validate_alignment_fields(0, 0, -1, 100, 3).is_ok());
        assert!(validate_alignment_fields(0, 2, 2, 0, 3).is_ok());
        assert!(validate_alignment_fields(0, -1, -1, -1, 3).is_ok());
        // ref_id past the dictionary
        assert!(validate_alignment_fields(0, 3, 0, 0, 3).is_err());
        // next_ref_id past the dictionary
        assert!(validate_alignment_fields(0, 0, 5, 0, 3).is_err());
        // pos below -1
        assert!(validate_alignment_fields(0, 0, 0, -2, 3).is_err());
        // ref_id below -1
        assert!(validate_alignment_fields(0, -2, 0, 0, 3).is_err());
    }

    #[test]
    fn fqz_multiblock_rejects_overlong_record_count() {
        // A block header claiming far more records than the chunk has must be
        // rejected before fqzcomp::decompress does vec![0; record_count].
        let record_count: u32 = 1_000_000; // bogus: far more than the chunk holds
        let fqz_blob: &[u8] = &[]; // empty payload, block_len = 0
        // Compute the genuine block CRC (over record_count || fqz_blob) so the
        // CRC check passes and ONLY the explicit record-count bound can reject it.
        let crc = {
            let mut h = flate2::Crc::new();
            h.update(&record_count.to_le_bytes());
            h.update(fqz_blob);
            h.sum()
        };
        let mut blob = Vec::new();
        blob.extend_from_slice(&1u32.to_le_bytes()); // num_blocks = 1
        blob.extend_from_slice(&(fqz_blob.len() as u32).to_le_bytes()); // block_len
        blob.extend_from_slice(&record_count.to_le_bytes());
        blob.extend_from_slice(&crc.to_le_bytes());
        blob.extend_from_slice(fqz_blob);
        let strands = vec![false; 2]; // only 2 real records in the chunk
        let err = decompress_quality_fqz_multiblock(&blob, &strands)
            .expect_err("overlong record_count must be rejected");
        // The explicit bound must fire BEFORE fqzcomp::decompress allocates
        // vec![0; record_count]; assert its message so we don't false-green on a
        // later check (e.g. the lengths-count mismatch).
        assert!(
            err.to_string().contains("exceeds remaining"),
            "expected record-count bound error, got: {err}"
        );
    }

    /// Build a single-block fqz multiblock blob with a genuine per-block CRC.
    fn fqz_blob_one_block(record_count: u32, payload: &[u8]) -> Vec<u8> {
        let crc = {
            let mut h = flate2::Crc::new();
            h.update(&record_count.to_le_bytes());
            h.update(payload);
            h.sum()
        };
        let mut blob = Vec::new();
        blob.extend_from_slice(&1u32.to_le_bytes()); // num_blocks = 1
        blob.extend_from_slice(&(payload.len() as u32).to_le_bytes()); // block_len
        blob.extend_from_slice(&record_count.to_le_bytes());
        blob.extend_from_slice(&crc.to_le_bytes());
        blob.extend_from_slice(payload);
        blob
    }

    #[test]
    fn validate_fqz_framing_accepts_intact_and_rejects_payload_corruption() {
        // The verify-time framing check must accept a well-formed blob and reject a
        // payload byte flip even though it never runs the fqz entropy decode — this
        // is the gap that previously let `verify` approve a corrupt fqz archive.
        let payload = b"some-fqz-bytes-stand-in";
        let blob = fqz_blob_one_block(3, payload);
        assert!(validate_fqz_multiblock_framing(&blob, 3).is_ok());

        // Flip a byte inside the payload (last byte) without updating the CRC.
        let mut corrupt = blob.clone();
        *corrupt.last_mut().unwrap() ^= 0xFF;
        let err = validate_fqz_multiblock_framing(&corrupt, 3)
            .expect_err("payload corruption must fail the per-block CRC");
        assert!(
            err.to_string().contains("CRC mismatch"),
            "expected CRC error, got: {err}"
        );
    }

    #[test]
    fn validate_fqz_framing_rejects_record_count_mismatch() {
        // Per-block record_counts must sum to the chunk's record count exactly.
        let blob = fqz_blob_one_block(3, b"abc");
        let err = validate_fqz_multiblock_framing(&blob, 4)
            .expect_err("record-count total mismatch must fail");
        assert!(
            err.to_string().contains("expected 4"),
            "expected total-count error, got: {err}"
        );
    }

    #[test]
    fn stream_decompress_cap_is_bounded_and_record_scaled() {
        // Floor for tiny chunks.
        assert_eq!(stream_decompress_cap(0), MIN_STREAM_DECOMPRESS_CAP);
        // Scales with the record count above the floor.
        assert_eq!(
            stream_decompress_cap(100_000),
            100_000 * MAX_DECOMPRESSED_BYTES_PER_RECORD
        );
        // Clamped by the absolute ceiling for an absurd (untrusted) record count,
        // so a bomb can't translate a huge declared count into an unbounded alloc.
        assert_eq!(
            stream_decompress_cap(usize::MAX),
            MAX_DECOMPRESSED_BYTES_PER_STREAM
        );
    }

    #[test]
    fn check_record_count_rejects_overlong_header() {
        // mapq stream holds exactly one byte per record, so it bounds num_records.
        assert!(check_record_count(0, 0).is_ok());
        assert!(check_record_count(100, 100).is_ok());
        assert!(check_record_count(50, 100).is_ok()); // fewer is fine (loop stops early)
        // A header claiming more records than the mapq stream can hold is corrupt —
        // must be rejected before the speculative with_capacity(num_records) alloc.
        assert!(check_record_count(101, 100).is_err());
        assert!(check_record_count(u32::MAX as usize, 10).is_err());
    }
}
