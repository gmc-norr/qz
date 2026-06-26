//! Paired-end FASTQ compression (chunked, footer directory).
use anyhow::{Result, bail};
use std::io::{Read, Seek, SeekFrom, Write};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::cli::{CompressConfig, DecompressConfig, QualityCompressor};
use crate::compression::codec_ids::{CODEC_BSC, CODEC_FQZCOMP, CODEC_QUALITY_CTX};

/// Process-wide monotonic nonce so two concurrent decodes IN ONE PROCESS to the
/// same prefix get distinct temp/backup filenames (PID alone would collide).
pub(crate) static TMP_NONCE: AtomicU64 = AtomicU64::new(0);

pub(crate) mod format;
pub(crate) mod header_delta;
pub(crate) mod stream;
pub(crate) mod streams;

// ---------------------------------------------------------------------------
// Config rejection
// ---------------------------------------------------------------------------

/// Reject any configuration the paired path does not yet support. One `bail!`
/// per concern so Task 10 can assert specific rejection messages.
pub(super) fn reject_unsupported(a: &CompressConfig) -> Result<()> {
    // Shared multi-stream gating; paired allows an explicit Fqzcomp quality codec
    // and stdout output. (stdout is INTENTIONALLY allowed — the v5 paired encoder is
    // seek-free; the public compress dispatch rejects only the reference+paired→stdout
    // combination before calling compress_paired_v5.)
    super::reject_unsupported_common(
        a,
        &super::ModeCaps {
            mode: "paired",
            allows_explicit_fqzcomp: true,
            allows_stdout: true,
        },
    )
}

// ---------------------------------------------------------------------------
// Per-mate quality resolution
// ---------------------------------------------------------------------------

/// Resolve the quality codec for one mate from its first chunk's read count,
/// mirroring single-end. Returns `(resolved_codec, directory_codec_byte)` where the
/// byte is paired's role-local quality codec id: `Bsc = 1`, `Fqzcomp = 6`.
pub(super) fn resolve_mate_quality(a: &CompressConfig, first_chunk_reads: usize) -> (QualityCompressor, u8) {
    let resolved = crate::compression::resolve_quality_compressor(
        a.advanced.quality_compressor,
        first_chunk_reads,
        a.quality_mode,
        a.no_quality,
        false, // paired does not support quality modeling/delta
    );
    let codec = match resolved {
        QualityCompressor::Fqzcomp => CODEC_FQZCOMP,
        _ => CODEC_BSC, // Bsc (and any other backend the paired pipeline stores as BSC)
    };
    (resolved, codec)
}

// ---------------------------------------------------------------------------
// v5 chunk-major paired encoder (seek-free)
// ---------------------------------------------------------------------------

pub(crate) fn compress_paired_v5(a: &CompressConfig) -> Result<()> {
    stream::compress_paired_streaming(a)
}

// ---------------------------------------------------------------------------
// Decode + output helpers (shared by v5 decompress and verify)
// ---------------------------------------------------------------------------

/// Append a suffix to a path's full OsString, producing `<prefix><suffix>`.
pub(crate) fn with_suffix(prefix: &std::path::Path, suffix: &str) -> std::path::PathBuf {
    let mut s = prefix.to_path_buf().into_os_string();
    s.push(suffix);
    std::path::PathBuf::from(s)
}

/// Drop guard: removes both temp files on drop unless disarmed.
pub(crate) struct TmpPair {
    pub(crate) a: std::path::PathBuf,
    pub(crate) b: std::path::PathBuf,
    pub(crate) armed: bool,
}
impl TmpPair {
    pub(crate) fn disarm(&mut self) {
        self.armed = false;
    }
}
impl Drop for TmpPair {
    fn drop(&mut self) {
        if self.armed {
            let _ = std::fs::remove_file(&self.a);
            let _ = std::fs::remove_file(&self.b);
        }
    }
}

/// Split a BSC-decompressed, varint-framed stream into exactly `n` records.
/// Mirrors `StreamCursor::read_one_record_varint` (decompress_impl.rs): each
/// record is `[varint(L)][L bytes]`. Asserts full consumption.
fn split_varint_records(raw: &[u8], n: u64) -> Result<Vec<Vec<u8>>> {
    use crate::compression::dna_utils::read_varint;
    // belt-and-suspenders: a record is >= 1 byte, so the real count can't
    // exceed raw.len(); cap the capacity hint so an untrusted n can't over-alloc.
    let mut out = Vec::with_capacity((n as usize).min(raw.len() + 1));
    let mut o = 0usize;
    for _ in 0..n {
        let len = read_varint(raw, &mut o)
            .ok_or_else(|| anyhow::anyhow!("truncated varint in record stream at {o}"))?
            as usize;
        let end = o
            .checked_add(len)
            .filter(|&e| e <= raw.len())
            .ok_or_else(|| anyhow::anyhow!("record claims {len} bytes past end of stream"))?;
        out.push(raw[o..end].to_vec());
        o = end;
    }
    if o != raw.len() {
        bail!("record stream has {} trailing bytes", raw.len() - o);
    }
    Ok(out)
}

/// Split a BSC-decompressed quality stream into `n` records. Mirrors the
/// single-end Bsc-quality path (`read_one_record_qual_varint` +
/// `unpack_qualities`): each record is `[varint(l_seq)][packed bytes]` where
/// `packed_len = (l_seq * bits_per_qual + 7) / 8`. For `QualityBinning::None`
/// the codec uses **7 bits/qual** (not 8), so `packed_len = ceil(l_seq*7/8)`,
/// which only coincidentally equals `l_seq` at certain lengths — it must be
/// computed, not assumed.
pub(crate) fn split_qual_records(raw: &[u8], n: u64) -> Result<Vec<Vec<u8>>> {
    use crate::compression::columnar::{QualityBinning, unpack_qualities};
    use crate::compression::dna_utils::read_varint;
    let bits_per_qual = QualityBinning::None.bits_per_quality();
    // belt-and-suspenders: a record is >= 1 byte, so the real count can't
    // exceed raw.len(); cap the capacity hint so an untrusted n can't over-alloc.
    let mut out = Vec::with_capacity((n as usize).min(raw.len() + 1));
    let mut o = 0usize;
    for _ in 0..n {
        let l_seq = read_varint(raw, &mut o)
            .ok_or_else(|| anyhow::anyhow!("truncated varint in qual stream at {o}"))?
            as usize;
        let packed_len = l_seq
            .checked_mul(bits_per_qual)
            .and_then(|n| n.checked_add(7))
            .map(|n| n / 8)
            .ok_or_else(|| anyhow::anyhow!("qual length overflow l_seq={l_seq}"))?;
        let end = o
            .checked_add(packed_len)
            .filter(|&e| e <= raw.len())
            .ok_or_else(|| {
                anyhow::anyhow!("qual record claims {packed_len} bytes past end of stream")
            })?;
        let packed = &raw[o..end];
        out.push(unpack_qualities(packed, l_seq, QualityBinning::None));
        o = end;
    }
    if o != raw.len() {
        bail!("qual stream has {} trailing bytes", raw.len() - o);
    }
    Ok(out)
}

/// Decode a chunk's columnar header entry into the record ids.
///
/// The entry is a block-stream of one or more blocks; each block's payload is
/// `[num_reads u32][columnar_blob]` for that block's own id slice. The common case
/// is a single block, but when a chunk's columnar blob would exceed the 64 MiB
/// block cap the encoder (`compress_impl::columnar_blocks_capped`) splits it into
/// self-contained sub-blocks. Decode each independently, concatenate in block
/// order, and verify the per-block counts sum to the directory's `record_count`.
pub(crate) fn decode_columnar_headers(bytes: &[u8], record_count: u64) -> Result<Vec<Vec<u8>>> {
    let pl = streams::decode_block_payloads(bytes)?;
    if pl.is_empty() {
        bail!("columnar header stream has no blocks");
    }
    // No `with_capacity(record_count)`: record_count is attacker-controlled footer
    // data, so growth is left to `extend`, bounded by the CRC-framed, MAX_BLOCK-
    // capped block bytes.
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut total: u64 = 0;
    for (rc, blob) in &pl {
        if blob.len() < 4 {
            bail!("columnar header payload too short");
        }
        let num_reads = u32::from_le_bytes(blob[0..4].try_into().unwrap());
        if num_reads != *rc {
            bail!("columnar header block num_reads {num_reads} != block record_count {rc}");
        }
        let block_ids = crate::compression::header_col::decompress_headers_columnar(
            &blob[4..],
            num_reads as usize,
        )?;
        if block_ids.len() as u32 != num_reads {
            bail!(
                "columnar header block decoded {} ids != {num_reads}",
                block_ids.len()
            );
        }
        total = total
            .checked_add(num_reads as u64)
            .ok_or_else(|| anyhow::anyhow!("columnar header block count overflow"))?;
        ids.extend(block_ids);
    }
    if total != record_count {
        bail!("columnar header total {total} != record_count {record_count}");
    }
    Ok(ids)
}

/// Decode one mate's sequence stream into per-record byte vectors.
fn decode_seq(bytes: &[u8], record_count: u64) -> Result<Vec<Vec<u8>>> {
    let payloads = streams::decode_block_payloads(bytes)?;
    let sum: u64 = payloads.iter().map(|(rc, _)| *rc as u64).sum();
    if sum != record_count {
        bail!("seq block record_count sum {sum} != {record_count}");
    }
    let raw = streams::decode_bsc_stream(
        bytes,
        super::codecs::stream_decode_cap(record_count as usize),
    )?;
    split_varint_records(&raw, record_count)
}

/// Decode one mate's quality stream (codec 1 = Bsc, codec 6 = Fqzcomp — a
/// block-stream of fqzcomp sub-blocks decoded in parallel). Used by the legacy
/// `decode_chunk_v5` path (mixed per-chunk R2-header codec). `_sequences` is
/// retained for signature symmetry but unused now that the sequence-context
/// `quality_ctx` codec (code 4) is removed.
pub(crate) fn decode_qual(
    bytes: &[u8],
    codec: u8,
    record_count: u64,
    _sequences: &[Vec<u8>],
) -> Result<Vec<Vec<u8>>> {
    if codec == CODEC_QUALITY_CTX {
        bail!("quality_ctx quality backend has been removed (code 4) — recompress");
    }
    let q = if codec == CODEC_FQZCOMP {
        // fqzcomp: a block-stream of per-sub-block packed blobs. decode_block_payloads runs
        // the strict per-block framing + CRC + no-trailing-bytes check; each blob decodes to a
        // varint-prefixed `[varint(len), ASCII]*` stream already inverse-permuted to original
        // order, so we split with split_varint_records (NOT split_qual_records). The
        // sub-blocks are independent and decoded IN PARALLEL; `par_iter().collect()`
        // preserves block order, so the final concat is still in original read order.
        use rayon::prelude::*;
        let blocks = streams::decode_block_payloads(bytes)?;
        // Fan across the pool only for genuinely multi-block streams. A single block
        // (every reference per-chunk-mate quality entry) gains nothing from `par_iter`
        // but recruits an extra worker whose one-time per-thread allocator-arena/stack
        // RSS is never reclaimed for the pool's lifetime — so spreading the one block
        // per chunk across the pool made reference decode peak RSS ramp with chunk
        // count. Decode it inline to keep the warmed-worker set (and peak RSS) bounded
        // and constant in read/chunk count.
        let per_block: Vec<Vec<Vec<u8>>> = if blocks.len() <= 1 {
            blocks
                .iter()
                .map(|(rc, blk)| -> Result<Vec<Vec<u8>>> {
                    let raw = crate::compression::codecs::decompress_qualities_fqzcomp(blk)?;
                    split_varint_records(&raw, *rc as u64)
                })
                .collect::<Result<Vec<_>>>()?
        } else {
            blocks
                .par_iter()
                .map(|(rc, blk)| -> Result<Vec<Vec<u8>>> {
                    let raw = crate::compression::codecs::decompress_qualities_fqzcomp(blk)?;
                    split_varint_records(&raw, *rc as u64) // exact rc + full-consumption checks
                })
                .collect::<Result<Vec<_>>>()?
        };
        // `record_count` is an untrusted u64 from the on-disk footer; bound the capacity
        // hint by the input size (each record needs >= 1 input byte) so a hostile archive
        // can't force an OOM-abort here. Mirrors split_varint_records / split_qual_records.
        let mut out: Vec<Vec<u8>> =
            Vec::with_capacity((record_count as usize).min(bytes.len() + 1));
        for blk_recs in per_block {
            out.extend(blk_recs);
        }
        out
    } else {
        let raw = streams::decode_bsc_stream(
            bytes,
            super::codecs::stream_decode_cap(record_count as usize),
        )?;
        split_qual_records(&raw, record_count)?
    };
    if q.len() as u64 != record_count {
        bail!("quality decoded {} records != {}", q.len(), record_count);
    }
    Ok(q)
}

/// Read `length` bytes at absolute `offset` from `file`. v5 entries call this
/// with their `ChunkDirEntry` offset/length.
fn read_entry_at(file: &mut std::fs::File, offset: u64, length: u64) -> Result<Vec<u8>> {
    file.seek(SeekFrom::Start(offset))?;
    let mut bytes = vec![0u8; length as usize];
    file.read_exact(&mut bytes)?;
    Ok(bytes)
}

/// Append `n` reconstructed FASTQ records (`@id\nSEQ\n+\nQUAL\n`) to `w`.
fn emit_records<W: Write>(
    w: &mut W,
    ids: &[Vec<u8>],
    seqs: &[Vec<u8>],
    quals: &[Vec<u8>],
) -> Result<()> {
    if seqs.len() != ids.len() || quals.len() != ids.len() {
        bail!(
            "paired decode: column count mismatch (ids {}, seqs {}, quals {})",
            ids.len(),
            seqs.len(),
            quals.len()
        );
    }
    for i in 0..ids.len() {
        // FASTQ requires the quality line to match the sequence length. Sequence and
        // quality come from independent streams, so a corrupt/crafted archive could
        // disagree — reject rather than emit a structurally invalid record. (The
        // streaming paired decoder enforces the same check; this covers the legacy
        // chunk path and deep verify.)
        if seqs[i].len() != quals[i].len() {
            bail!(
                "paired decode: record {i} quality length {} != sequence length {} (corrupt archive)",
                quals[i].len(),
                seqs[i].len()
            );
        }
        // id already includes the leading '@'.
        w.write_all(&ids[i])?;
        w.write_all(b"\n")?;
        w.write_all(&seqs[i])?;
        w.write_all(b"\n+\n")?;
        w.write_all(&quals[i])?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

/// Emit R1[i] then R2[i] alternately to ONE writer (interleaved FASTQ:
/// `r0/1, r0/2, r1/1, …`). Used by the rare mixed-codec legacy interleaved path;
/// the streaming path interleaves via `PairedSink::Interleaved`.
fn emit_records_interleaved<W: Write + ?Sized>(
    w: &mut W,
    r1_ids: &[Vec<u8>],
    r1_seqs: &[Vec<u8>],
    r1_quals: &[Vec<u8>],
    r2_ids: &[Vec<u8>],
    r2_seqs: &[Vec<u8>],
    r2_quals: &[Vec<u8>],
) -> Result<()> {
    if r1_ids.len() != r2_ids.len() {
        bail!(
            "paired interleaved decode: R1/R2 record count mismatch ({} vs {})",
            r1_ids.len(),
            r2_ids.len()
        );
    }
    for i in 0..r1_ids.len() {
        for (ids, seqs, quals) in [
            (r1_ids, r1_seqs, r1_quals),
            (r2_ids, r2_seqs, r2_quals),
        ] {
            if seqs[i].len() != quals[i].len() {
                bail!(
                    "paired interleaved decode: record {i} quality length {} != sequence length {} (corrupt archive)",
                    quals[i].len(),
                    seqs[i].len()
                );
            }
            w.write_all(&ids[i])?;
            w.write_all(b"\n")?;
            w.write_all(&seqs[i])?;
            w.write_all(b"\n+\n")?;
            w.write_all(&quals[i])?;
            w.write_all(b"\n")?;
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// v5 (unified chunk-major) decode.
// ---------------------------------------------------------------------------

/// One chunk's reconstructed columns: (R1 ids/seq/qual, R2 ids/seq/qual).
type DecodedChunk = (
    (Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>),
    (Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<Vec<u8>>),
);

/// Decode one chunk's six unified streams (v5 paired): R1 headers→seq→qual, then
/// R2 header [indep|delta], seq, qual, using the unified directory's (mate, role)
/// entries.
fn decode_chunk_v5(
    file: &mut std::fs::File,
    ents: &[&crate::compression::chunk_directory::ChunkDirEntry],
    chunk_idx: u32,
) -> Result<DecodedChunk> {
    use crate::compression::chunk_directory::StreamRole;
    let find = |mate: u8, role: StreamRole| {
        ents.iter()
            .copied()
            .find(|e| e.mate == mate && e.role == role)
    };

    // R1: headers (delta needs r1_ids) -> seq -> qual.
    let e_r1h = find(1, StreamRole::Headers)
        .ok_or_else(|| anyhow::anyhow!("chunk {chunk_idx}: no R1 headers"))?;
    let r1_ids = decode_columnar_headers(
        &read_entry_at(file, e_r1h.offset, e_r1h.length)?,
        e_r1h.record_count,
    )?;

    let e_r1s = find(1, StreamRole::Sequence)
        .ok_or_else(|| anyhow::anyhow!("chunk {chunk_idx}: no R1 seq"))?;
    let r1_seq = decode_seq(
        &read_entry_at(file, e_r1s.offset, e_r1s.length)?,
        e_r1s.record_count,
    )?;

    let e_r1q = find(1, StreamRole::Qual)
        .ok_or_else(|| anyhow::anyhow!("chunk {chunk_idx}: no R1 qual"))?;
    let r1_qual = decode_qual(
        &read_entry_at(file, e_r1q.offset, e_r1q.length)?,
        e_r1q.codec,
        e_r1q.record_count,
        &r1_seq,
    )?;

    // R2 header: independent columnar (role Headers) OR delta (role HeaderDelta).
    let r2_ids = if let Some(e) = find(2, StreamRole::Headers) {
        decode_columnar_headers(&read_entry_at(file, e.offset, e.length)?, e.record_count)?
    } else if let Some(e) = find(2, StreamRole::HeaderDelta) {
        let ops = streams::decode_bsc_stream(
            &read_entry_at(file, e.offset, e.length)?,
            super::codecs::stream_decode_cap(e.record_count as usize),
        )?;
        let ids = header_delta::decode(&ops, &r1_ids, e.record_count as usize)?;
        if ids.len() as u64 != e.record_count {
            bail!("R2 delta decoded {} ids != {}", ids.len(), e.record_count);
        }
        ids
    } else {
        bail!("chunk {chunk_idx}: no R2 header role");
    };

    let e_r2s = find(2, StreamRole::Sequence)
        .ok_or_else(|| anyhow::anyhow!("chunk {chunk_idx}: no R2 seq"))?;
    let r2_seq = decode_seq(
        &read_entry_at(file, e_r2s.offset, e_r2s.length)?,
        e_r2s.record_count,
    )?;

    let e_r2q = find(2, StreamRole::Qual)
        .ok_or_else(|| anyhow::anyhow!("chunk {chunk_idx}: no R2 qual"))?;
    let r2_qual = decode_qual(
        &read_entry_at(file, e_r2q.offset, e_r2q.length)?,
        e_r2q.codec,
        e_r2q.record_count,
        &r2_seq,
    )?;

    // Joint length check.
    let rc = e_r1h.record_count as usize;
    for (name, len) in [
        ("r1_ids", r1_ids.len()),
        ("r1_seq", r1_seq.len()),
        ("r1_qual", r1_qual.len()),
        ("r2_ids", r2_ids.len()),
        ("r2_seq", r2_seq.len()),
        ("r2_qual", r2_qual.len()),
    ] {
        if len != rc {
            bail!("chunk {chunk_idx}: {name} len {len} != record_count {rc}");
        }
    }
    Ok(((r1_ids, r1_seq, r1_qual), (r2_ids, r2_seq, r2_qual)))
}

/// Group unified directory entries by chunk index for in-order processing.
fn entries_by_chunk_v5(
    dir: &crate::compression::chunk_directory::ChunkDirectory,
) -> std::collections::BTreeMap<u32, Vec<&crate::compression::chunk_directory::ChunkDirEntry>> {
    let mut by_chunk: std::collections::BTreeMap<
        u32,
        Vec<&crate::compression::chunk_directory::ChunkDirEntry>,
    > = Default::default();
    for e in &dir.entries {
        by_chunk.entry(e.chunk_index).or_default().push(e);
    }
    by_chunk
}

/// Parse + validate a paired (`archive_type` 1) v5 archive's front header and footer,
/// returning `(header, directory, use_legacy)`. `use_legacy` is true for the rare
/// MIXED per-chunk R2-header codec (columnar on some chunks, delta on others), which
/// the single-codec-per-stream streaming engine can't express — real Illumina paired
/// is uniform (all-delta), so it streams. Shared by `decompress_paired_v5` (two files)
/// and `decompress_paired_interleaved_v5` (one interleaved stream) so the header/footer
/// contract + mixed-codec gate live in ONE place.
fn parse_paired_v5_prelude(
    input: &std::path::Path,
) -> Result<(
    crate::compression::FixedHeader,
    crate::compression::chunk_directory::ChunkDirectory,
    bool,
)> {
    use crate::compression::chunk_directory::StreamRole;
    let mut hbuf = [0u8; 64];
    let mut hf = std::fs::File::open(input)?;
    let n = std::io::Read::read(&mut hf, &mut hbuf)?;
    if n < crate::compression::V2_PREFIX_SIZE {
        bail!("paired v5 header truncated");
    }
    let hdr = crate::compression::FixedHeader::parse_v5(&hbuf[..n])?;
    if hdr.archive_type != 1 {
        bail!("not a paired v5 archive (archive_type {})", hdr.archive_type);
    }
    // Paired archives are ALWAYS written with encoding_type 0 (Raw): no per-record
    // sequence-hint byte, default 64 MiB block cap. The streaming decode derives both
    // the block-size DoS cap and `has_sequence_hints` from encoding_type, so a forged
    // non-zero value would mis-set the cap (ultra=10 → 768 MiB) or strip a leading byte
    // from every sequence (RawWithHints=4) → silent wrong output. Reject it here
    // (matches the single-end engine's hardening posture).
    if hdr.encoding_type != 0 {
        bail!(
            "paired v5 archive has unexpected encoding_type {} (paired is always 0/Raw)",
            hdr.encoding_type
        );
    }
    let header_end = u32::from_le_bytes(hbuf[4..8].try_into().unwrap()) as u64;
    let dir = crate::compression::chunk_directory::read_v5_footer(input, header_end)?;
    crate::compression::paired::format::validate_paired_directory(&dir)?;
    let r2_has_columnar = dir
        .entries
        .iter()
        .any(|e| e.mate == 2 && e.role == StreamRole::Headers);
    let r2_has_delta = dir
        .entries
        .iter()
        .any(|e| e.mate == 2 && e.role == StreamRole::HeaderDelta);
    Ok((hdr, dir, r2_has_columnar && r2_has_delta))
}

pub fn decompress_paired_v5(a: &DecompressConfig) -> Result<()> {
    use std::io::BufWriter;

    if a.output.is_empty() {
        bail!("paired decompress needs an output prefix");
    }
    let prefix = &a.output[0];
    if crate::cli::is_stdio_path(prefix) {
        bail!("paired decompress needs a named output prefix (no stdout)");
    }
    // Reject rather than silently ignore --gzipped: paired mode writes plain
    // _R1.fastq/_R2.fastq and does not yet support gz output.
    if a.gzipped {
        bail!("gzipped output not supported in paired mode (writes plain _R1.fastq/_R2.fastq)");
    }
    let out1 = with_suffix(prefix, "_R1.fastq");
    let out2 = with_suffix(prefix, "_R2.fastq");

    if !a.force && (out1.exists() || out2.exists()) {
        bail!(
            "Output file already exists: {} or {}\nUse --force to overwrite",
            out1.display(),
            out2.display()
        );
    }

    // ---- v5 front header + footer + mixed-codec gate (shared with interleaved) ----
    let (hdr, dir, use_legacy) = parse_paired_v5_prelude(&a.input)?;

    // ---- Temp output files + drop guard ----
    // Draw one nonce per call; reuse it for tmp1/tmp2 AND the .bak names in
    // publish_pair so concurrent in-process decodes to the same prefix can't
    // collide (PID alone is shared across threads).
    let pid = std::process::id();
    let nonce = TMP_NONCE.fetch_add(1, Ordering::Relaxed);
    let tmp1 = with_suffix(&out1, &format!(".{pid}.{nonce}.tmp"));
    let tmp2 = with_suffix(&out2, &format!(".{pid}.{nonce}.tmp"));
    let mut guard = TmpPair {
        a: tmp1.clone(),
        b: tmp2.clone(),
        armed: true,
    };
    // O_EXCL parity with the single-end path: the names are already pid+nonce
    // unique, so create_new also rejects a pre-planted symlink/file.
    let mut w1 = BufWriter::new(super::create_new_for_write(&tmp1)?);
    let mut w2 = BufWriter::new(super::create_new_for_write(&tmp2)?);

    // ---- Decode ----
    // The unified block-streaming engine (shared with single-end: parallel, bounded
    // RAM, single-end parity) handles the production path; `use_legacy` (from the
    // prelude) takes the byte-identical serial chunk path only for the rare mixed
    // per-chunk R2-header codec the streaming engine can't express.
    if use_legacy {
        let mut file = std::fs::File::open(&a.input)?;
        let by_chunk = entries_by_chunk_v5(&dir);
        for chunk_idx in 0..dir.num_chunks {
            let ents = by_chunk
                .get(&chunk_idx)
                .ok_or_else(|| anyhow::anyhow!("missing chunk {chunk_idx}"))?;
            let ((r1_ids, r1_seq, r1_qual), (r2_ids, r2_seq, r2_qual)) =
                decode_chunk_v5(&mut file, ents, chunk_idx)?;
            emit_records(&mut w1, &r1_ids, &r1_seq, &r1_qual)?;
            emit_records(&mut w2, &r2_ids, &r2_seq, &r2_qual)?;
        }
    } else {
        crate::compression::decompress_impl::decompress_paired_streaming_v5(
            &a.input,
            &hdr,
            &dir,
            &mut crate::compression::decompress_impl::PairedSink::Two(&mut w1, &mut w2),
        )?;
    }

    w1.flush()?;
    w2.flush()?;
    drop(w1);
    drop(w2);

    // ---- Atomic publish (best-effort backup-and-rollback) ----
    publish_pair(&tmp1, &tmp2, &out1, &out2, nonce)?;
    guard.disarm();
    Ok(())
}

/// Decode a paired (`archive_type` 1) v5 archive as a single **interleaved** FASTQ
/// stream — R1/R2 records alternating — to `a.output[0]` (a file, `-`/stdout, or
/// gzipped). Reuses the streaming two-mate decoder via `PairedSink::Interleaved`;
/// the rare mixed per-chunk R2-header-codec archive interleaves per chunk through
/// the legacy `decode_chunk_v5`. The single-sink dispatch owns the atomic-temp /
/// stdout / gzip output machinery (`force`/exists is checked by the caller).
pub fn decompress_paired_interleaved_v5(a: &DecompressConfig) -> Result<()> {
    use crate::compression::decompress_impl::{decompress_paired_streaming_v5, PairedSink};
    use std::time::Instant;

    if a.output.is_empty() {
        bail!("interleaved paired decompress needs one output (-o <file>|-)");
    }

    // Header/footer + mixed-codec gate via the SAME prelude decompress_paired_v5 uses.
    let (hdr, dir, use_legacy) = parse_paired_v5_prelude(&a.input)?;

    // R1 + R2 records (cosmetic, for the progress log).
    let num_records = dir.num_reads.saturating_mul(2) as usize;

    let input = a.input.clone();
    crate::compression::decompress_impl::decompress_streaming_bsc_dispatch(
        a,
        Instant::now(),
        num_records,
        move |w| {
            if use_legacy {
                let mut file = std::fs::File::open(&input)?;
                let by_chunk = entries_by_chunk_v5(&dir);
                for chunk_idx in 0..dir.num_chunks {
                    let ents = by_chunk
                        .get(&chunk_idx)
                        .ok_or_else(|| anyhow::anyhow!("missing chunk {chunk_idx}"))?;
                    let ((r1_ids, r1_seq, r1_qual), (r2_ids, r2_seq, r2_qual)) =
                        decode_chunk_v5(&mut file, ents, chunk_idx)?;
                    emit_records_interleaved(
                        w, &r1_ids, &r1_seq, &r1_qual, &r2_ids, &r2_seq, &r2_qual,
                    )?;
                }
                Ok(())
            } else {
                decompress_paired_streaming_v5(
                    &input,
                    &hdr,
                    &dir,
                    &mut PairedSink::Interleaved(w),
                )
            }
        },
    )
}

// ---------------------------------------------------------------------------
// Verify
// ---------------------------------------------------------------------------

/// A `Write` sink that feeds bytes into a `flate2::Crc` and counts them,
/// discarding the data. Lets `emit_records` drive the deep-verify CRC over the
/// reconstructed FASTQ bytes without materialising them. `pub(crate)` so the
/// reference deep-verify path can build two per-mate sinks from it.
pub(crate) struct CrcWriter {
    pub(crate) crc: flate2::Crc,
    pub(crate) bytes: u64,
}
impl Write for CrcWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.crc.update(buf);
        self.bytes += buf.len() as u64;
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Verify a v5 (chunk-major unified container) paired archive.
///
/// **Deep** (`fast == false`): stream-decode both mates per chunk via
/// `decode_chunk_v5`, feeding every reconstructed record (R1 then R2, in chunk
/// order) into one CRC32 over the FASTQ bytes.
///
/// **Fast** (`fast == true`): for every directory entry, read its bytes and
/// CRC-verify each inline block frame via `decode_block_payloads` WITHOUT
/// BSC-decoding; accumulate block + byte counts and cross-check the per-entry
/// record-count sum.
/// Cosmetic `encoding_type` reported in a paired archive's `VerifyResult`
/// (paired archives are identified on the wire by v5 `archive_type == 1`, not by
/// this byte; this is purely a display tag for `qz verify`).
const PAIRED_REPORTED_ENCODING: u8 = 0xF0;

pub fn verify_paired_v5(
    input: &std::path::Path,
    fast: bool,
) -> Result<crate::compression::VerifyResult> {
    use crate::compression::chunk_directory::{StreamRole, GLOBAL_SENTINEL};
    use crate::compression::{VerifyMode, VerifyResult};
    let started = std::time::Instant::now();

    // v5 front header → header_end (footer byte-range source) lives in the
    // locator, but read_v5_footer needs header_end as its scan anchor.
    let mut prefix = [0u8; 8];
    {
        let mut hf = std::fs::File::open(input)?;
        let mut got = 0usize;
        while got < 8 {
            match Read::read(&mut hf, &mut prefix[got..])? {
                0 => break,
                k => got += k,
            }
        }
        if got < 8 {
            bail!("paired v5 header truncated");
        }
    }
    let header_end = u32::from_le_bytes(prefix[4..8].try_into().unwrap()) as u64;
    let dir = crate::compression::chunk_directory::read_v5_footer(input, header_end)?;
    crate::compression::paired::format::validate_paired_directory(&dir)?;
    let mut file = std::fs::File::open(input)?;

    // 2 mates per pair; num_reads (footer) is the per-lane pair count.
    let num_reads = dir.num_reads.saturating_mul(2) as usize;

    if fast {
        let mut blocks_verified: u32 = 0;
        let mut total_bytes: u64 = 0;
        let mut headers_len: usize = 0;
        let mut sequences_len: usize = 0;
        let mut qualities_len: usize = 0;

        for e in &dir.entries {
            // The whole-archive ChunkDecodedSizes global is a bare single-end frame, not a
            // paired `write_segment_blob` blob, so it can't be parsed as a paired segment.
            // It's a tiny raw table (re-validated by read_chunk_layout on the direct-write
            // path); skip it here.
            if e.chunk_index == GLOBAL_SENTINEL {
                continue;
            }
            let bytes = read_entry_at(&mut file, e.offset, e.length)?;
            let pls = streams::decode_block_payloads(&bytes)?; // CRC-verifies every block
            blocks_verified += pls.len() as u32;
            total_bytes += bytes.len() as u64;
            let rc_sum: u64 = pls.iter().map(|(rc, _)| *rc as u64).sum();
            if rc_sum != e.record_count {
                bail!(
                    "entry (chunk {}, mate {}, role {:?}) block record_count sum {} != {}",
                    e.chunk_index,
                    e.mate,
                    e.role,
                    rc_sum,
                    e.record_count
                );
            }
            match e.role {
                StreamRole::Headers | StreamRole::HeaderDelta => headers_len += bytes.len(),
                StreamRole::Sequence => sequences_len += bytes.len(),
                StreamRole::Qual => qualities_len += bytes.len(),
                other => bail!("unexpected paired v5 stream role {other:?}"),
            }
        }

        return Ok(VerifyResult {
            num_reads,
            encoding_type: PAIRED_REPORTED_ENCODING, // not a single-end EncodingType; report a paired tag
            header_compressor: crate::cli::HeaderCompressor::Columnar,
            quality_compressor: QualityCompressor::Auto,
            headers_compressed_len: headers_len,
            sequences_compressed_len: sequences_len,
            qualities_compressed_len: qualities_len,
            crc32: 0,
            r2_crc32: None, // fast verify does not separate the R1/R2 reconstruction
            total_bytes,
            blocks_verified,
            mode: VerifyMode::Fast,
            elapsed_secs: started.elapsed().as_secs_f64(),
        });
    }

    // ---- Deep: reconstruct + CRC over FASTQ bytes ----
    let by_chunk = entries_by_chunk_v5(&dir);
    let mut sink = CrcWriter {
        crc: flate2::Crc::new(),
        bytes: 0,
    };
    // Compressed-byte sums by stream family, computed from the directory.
    let mut headers_len: usize = 0;
    let mut sequences_len: usize = 0;
    let mut qualities_len: usize = 0;
    for e in &dir.entries {
        if e.chunk_index == GLOBAL_SENTINEL {
            continue; // whole-archive ChunkDecodedSizes global; no stream family
        }
        match e.role {
            StreamRole::Headers | StreamRole::HeaderDelta => headers_len += e.length as usize,
            StreamRole::Sequence => sequences_len += e.length as usize,
            StreamRole::Qual => qualities_len += e.length as usize,
            other => bail!("unexpected paired v5 stream role {other:?}"),
        }
    }

    // Reconstruct through the SAME dispatch the production decoder uses, so deep verify
    // validates what `decompress` actually runs: the unified block-streaming engine for
    // uniform archives (the common case), and the legacy serial chunk path only for a
    // mixed per-chunk R2-header codec (which the single-codec-per-stream engine cannot
    // express). Each mate folds into its own CRC sink (R1 → crc32, R2 → r2_crc32).
    let mut sink2 = CrcWriter {
        crc: flate2::Crc::new(),
        bytes: 0,
    };
    let r2_has_columnar = dir
        .entries
        .iter()
        .any(|e| e.mate == 2 && e.role == StreamRole::Headers);
    let r2_has_delta = dir
        .entries
        .iter()
        .any(|e| e.mate == 2 && e.role == StreamRole::HeaderDelta);
    if r2_has_columnar && r2_has_delta {
        for chunk_idx in 0..dir.num_chunks {
            let ents = by_chunk
                .get(&chunk_idx)
                .ok_or_else(|| anyhow::anyhow!("missing chunk {chunk_idx}"))?;
            let ((r1_ids, r1_seq, r1_qual), (r2_ids, r2_seq, r2_qual)) =
                decode_chunk_v5(&mut file, ents, chunk_idx)?;
            emit_records(&mut sink, &r1_ids, &r1_seq, &r1_qual)?;
            emit_records(&mut sink2, &r2_ids, &r2_seq, &r2_qual)?;
        }
    } else {
        let (hdr, _header_end) = crate::compression::decompress_impl::read_fixed_header(input)?;
        crate::compression::decompress_impl::decompress_paired_streaming_v5(
            input,
            &hdr,
            &dir,
            &mut crate::compression::decompress_impl::PairedSink::Two(&mut sink, &mut sink2),
        )?;
    }

    Ok(VerifyResult {
        num_reads,
        encoding_type: PAIRED_REPORTED_ENCODING,
        header_compressor: crate::cli::HeaderCompressor::Columnar,
        quality_compressor: QualityCompressor::Auto,
        headers_compressed_len: headers_len,
        sequences_compressed_len: sequences_len,
        qualities_compressed_len: qualities_len,
        crc32: sink.crc.sum(),
        // Each mate has its own CRC sink now (R1 → crc32, R2 → r2_crc32), so deep verify
        // can route the uniform case through the production streaming engine.
        r2_crc32: Some(sink2.crc.sum()),
        total_bytes: sink.bytes + sink2.bytes,
        blocks_verified: 0,
        mode: VerifyMode::Deep,
        elapsed_secs: started.elapsed().as_secs_f64(),
    })
}

/// Atomically publish two temp files to their final paths. Under --force the
/// targets may already exist; back them up first so a mid-publish failure can
/// roll back. Best-effort: a failed backup-delete on the success path is a
/// warning, not an error.
pub(crate) fn publish_pair(
    tmp1: &std::path::Path,
    tmp2: &std::path::Path,
    out1: &std::path::Path,
    out2: &std::path::Path,
    nonce: u64,
) -> Result<()> {
    let pid = std::process::id();
    let bak1 = with_suffix(out1, &format!(".{pid}.{nonce}.bak"));
    let bak2 = with_suffix(out2, &format!(".{pid}.{nonce}.bak"));

    let mut backed1 = false;
    let mut backed2 = false;
    if out1.exists() {
        std::fs::rename(out1, &bak1)?;
        backed1 = true;
    }
    if out2.exists() {
        if let Err(e) = std::fs::rename(out2, &bak2) {
            if backed1 {
                let _ = std::fs::rename(&bak1, out1); // roll back
            }
            return Err(e.into());
        }
        backed2 = true;
    }

    // Commit R1.
    if let Err(e) = std::fs::rename(tmp1, out1) {
        if backed1 {
            let _ = std::fs::rename(&bak1, out1);
        }
        if backed2 {
            let _ = std::fs::rename(&bak2, out2);
        }
        return Err(e.into());
    }
    // Commit R2.
    if let Err(e) = std::fs::rename(tmp2, out2) {
        let _ = std::fs::remove_file(out1); // undo R1 commit
        if backed1 {
            let _ = std::fs::rename(&bak1, out1);
        }
        if backed2 {
            let _ = std::fs::rename(&bak2, out2);
        }
        return Err(e.into());
    }

    // Success: drop the backups (best-effort).
    if backed1 {
        let _ = std::fs::remove_file(&bak1);
    }
    if backed2 {
        let _ = std::fs::remove_file(&bak2);
    }
    Ok(())
}

#[cfg(test)]
mod fqz_quality_wiring_tests {
    use super::*;
    use crate::cli::{CompressConfig, QualityCompressor};

    fn cfg(qc: QualityCompressor) -> CompressConfig {
        let mut c = CompressConfig::default();
        c.advanced.quality_compressor = qc;
        c
    }

    #[test]
    fn explicit_fqzcomp_resolves_to_codec_6() {
        let c = cfg(QualityCompressor::Fqzcomp);
        let (resolved, codec) = resolve_mate_quality(&c, 10);
        assert!(matches!(resolved, QualityCompressor::Fqzcomp));
        assert_eq!(codec, 6, "fqzcomp must map to paired quality codec 6");
    }

    #[test]
    fn auto_small_input_stays_bsc() {
        let c = cfg(QualityCompressor::Auto);
        let (resolved, codec) = resolve_mate_quality(&c, 10);
        assert!(matches!(resolved, QualityCompressor::Bsc));
        assert_eq!(codec, 1);
    }

    #[test]
    fn auto_large_input_picks_fqzcomp() {
        // The production path: Auto with >= MIN_READS_QUALITY_CTX reads/mate now
        // resolves to fqzcomp (codec 6) in paired mode.
        let c = cfg(QualityCompressor::Auto);
        let big = 100_001; // MIN_READS_QUALITY_CTX is 100_000
        let (resolved, codec) = resolve_mate_quality(&c, big);
        assert!(matches!(resolved, QualityCompressor::Fqzcomp));
        assert_eq!(codec, 6);
    }
}

#[cfg(test)]
mod fqzcomp_decode_tests {
    use super::*; // REQUIRED: resolves decode_qual, split_varint_records, streams, etc.

    #[test]
    fn decode_qual_fqzcomp_multiblock_roundtrip() {
        use crate::compression::codecs::compress_qualities_fqzcomp_quals;
        use crate::compression::paired::streams::write_block_stream;
        use crate::io::FastqRecord;

        // 7 reads across 3 framed sub-blocks of 3,3,1 (exercises concat across blocks).
        let quals: Vec<Vec<u8>> = (0..7u8)
            .map(|i| vec![b'I' + (i % 5), b'5' + (i % 3), b'#'])
            .collect();
        let recs: Vec<FastqRecord> = quals
            .iter()
            .map(|q| FastqRecord {
                id: b"@r".to_vec(),
                sequence: b"ACG".to_vec(),
                quality: Some(q.clone()),
            })
            .collect();

        let mut blocks: Vec<(u32, Vec<u8>)> = Vec::new();
        for chunk in quals.chunks(3) {
            let slice: Vec<&[u8]> = chunk.iter().map(|q| q.as_slice()).collect();
            let blob = compress_qualities_fqzcomp_quals(&slice).unwrap();
            blocks.push((chunk.len() as u32, blob));
        }
        let mut framed = Vec::new();
        write_block_stream(&mut framed, &blocks);

        let seqs: Vec<Vec<u8>> = recs.iter().map(|r| r.sequence.clone()).collect();
        let out = decode_qual(&framed, 6, recs.len() as u64, &seqs).unwrap();
        assert_eq!(
            out, quals,
            "fqzcomp multi-block decode must roundtrip in original order"
        );
    }

    #[test]
    fn decode_qual_fqzcomp_rejects_corrupt_inner() {
        use crate::compression::codecs::compress_qualities_fqzcomp_quals;
        use crate::compression::paired::streams::write_block_stream;
        use crate::io::FastqRecord;

        let recs: Vec<FastqRecord> = (0..5u8)
            .map(|i| FastqRecord {
                id: b"@r".to_vec(),
                sequence: b"ACG".to_vec(),
                quality: Some(vec![b'I', b'5', b'0' + i]),
            })
            .collect();
        let quals: Vec<&[u8]> = recs.iter().map(|r| r.quality.as_deref().unwrap()).collect();
        let mut blob = compress_qualities_fqzcomp_quals(&quals).unwrap();
        // Append trailing garbage to the inner blob (outer CRC is recomputed over it, so the
        // framing stays valid; only the inner exact-consumption check should fire).
        blob.extend_from_slice(b"GARBAGE");
        let mut framed = Vec::new();
        write_block_stream(&mut framed, &[(recs.len() as u32, blob)]);

        let seqs: Vec<Vec<u8>> = recs.iter().map(|r| r.sequence.clone()).collect();
        let res = decode_qual(&framed, 6, recs.len() as u64, &seqs);
        assert!(
            res.is_err(),
            "trailing inner bytes must be rejected, not silently accepted"
        );
    }
}

#[cfg(test)]
mod v5_roundtrip_tests {
    //! First real exercise of the v5 paired encode→decode path. Each test builds
    //! R1/R2 FASTQ files in a TempDir, runs `compress_paired_v5` then
    //! `decompress_paired_v5`, and asserts the emitted `<prefix>_R1.fastq` /
    //! `<prefix>_R2.fastq` are BYTE-IDENTICAL to the inputs. The decoded footer is
    //! returned so callers can assert which R2-header path (indep vs delta) ran and
    //! the per-stream codec.
    use super::*;
    use crate::cli::{CompressConfig, DecompressConfig, QualityCompressor};
    use crate::compression::chunk_directory::{
        ChunkDirectory, ENTRY_BYTES, StreamRole, footer_crc32,
    };
    use tempfile::TempDir;

    /// Build one FASTQ string from `(id_without_@, seq, qual)` triples. The stored
    /// id round-trips WITH a leading '@' (FastqReader keeps it, `emit_records`
    /// re-emits it), so we write `@{id}\n{seq}\n+\n{qual}\n` — the exact framing the
    /// decoder produces — to make a byte-identical comparison meaningful.
    fn fastq_string(recs: &[(String, String, String)]) -> Vec<u8> {
        let mut s = String::new();
        for (id, seq, qual) in recs {
            s.push('@');
            s.push_str(id);
            s.push('\n');
            s.push_str(seq);
            s.push_str("\n+\n");
            s.push_str(qual);
            s.push('\n');
        }
        s.into_bytes()
    }

    /// Read the v5 front header's `header_size` (u32 LE at bytes [4..8]) — the
    /// `header_end` argument `read_v5_footer` needs.
    fn header_end_of(archive: &std::path::Path) -> u64 {
        let bytes = std::fs::read(archive).unwrap();
        u32::from_le_bytes(bytes[4..8].try_into().unwrap()) as u64
    }

    /// Full encode→decode roundtrip. Returns the decoded footer `ChunkDirectory`.
    fn roundtrip(
        r1: &[(String, String, String)],
        r2: &[(String, String, String)],
        chunk_records: usize,
        quality: Option<QualityCompressor>,
    ) -> ChunkDirectory {
        let dir = TempDir::new().unwrap();
        let r1_path = dir.path().join("r1.fastq");
        let r2_path = dir.path().join("r2.fastq");
        let archive = dir.path().join("paired.qz");
        let prefix = dir.path().join("out");

        let r1_bytes = fastq_string(r1);
        let r2_bytes = fastq_string(r2);
        std::fs::write(&r1_path, &r1_bytes).unwrap();
        std::fs::write(&r2_path, &r2_bytes).unwrap();

        let mut cfg = CompressConfig {
            input: vec![r1_path, r2_path],
            output: archive.clone(),
            working_dir: dir.path().to_path_buf(),
            threads: 1,
            ..CompressConfig::default()
        };
        cfg.advanced.chunk_records = chunk_records;
        if let Some(q) = quality {
            cfg.advanced.quality_compressor = q;
        }
        super::compress_paired_v5(&cfg).expect("compress_paired_v5 failed");

        // Read the footer for the path/codec assertions before decoding.
        let header_end = header_end_of(&archive);
        let footer =
            crate::compression::chunk_directory::read_v5_footer(&archive, header_end).unwrap();

        let dcfg = DecompressConfig {
            input: archive,
            output: vec![prefix.clone()],
            working_dir: dir.path().to_path_buf(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 6,
            force: true,
        };
        super::decompress_paired_v5(&dcfg).expect("decompress_paired_v5 failed");

        let out_r1 = with_suffix(&prefix, "_R1.fastq");
        let out_r2 = with_suffix(&prefix, "_R2.fastq");
        assert_eq!(
            std::fs::read(&out_r1).unwrap(),
            r1_bytes,
            "R1 output not byte-identical to input"
        );
        assert_eq!(
            std::fs::read(&out_r2).unwrap(),
            r2_bytes,
            "R2 output not byte-identical to input"
        );
        footer
    }

    /// 20 paired reads where R2 ids are `readN/2` against R1 `readN/1` — near
    /// identical, so the header-delta encoding beats independent columnar.
    fn delta_pairs(n: usize) -> (Vec<(String, String, String)>, Vec<(String, String, String)>) {
        let mut r1 = Vec::new();
        let mut r2 = Vec::new();
        for i in 0..n {
            let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT".to_string();
            let qual = "IIHHJJBBCCDDEEFFGGAABBCCDDEEFFGG".to_string();
            r1.push((format!("read{i}/1"), seq.clone(), qual.clone()));
            r2.push((format!("read{i}/2"), seq, qual));
        }
        (r1, r2)
    }

    #[test]
    fn v5_roundtrip_r2_header_delta() {
        let (r1, r2) = delta_pairs(20);
        // Single chunk (chunk_records large).
        let footer = roundtrip(&r1, &r2, 1000, None);
        assert!(
            footer
                .entries
                .iter()
                .any(|e| e.mate == 2 && e.role == StreamRole::HeaderDelta),
            "expected R2 header DELTA path, footer roles: {:?}",
            footer
                .entries
                .iter()
                .map(|e| (e.mate, e.role))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn v5_roundtrip_r2_header_independent() {
        // R1 ids are short `readN/1`; R2 ids are a CLEAN Illumina format whose
        // numeric columns (tile/x/y) march monotonically. The columnar header
        // encoder delta-codes those columns → tiny independent encoding. The
        // delta-vs-R1 stream, by contrast, must store every R2 id VERBATIM as a
        // literal op (R1 and R2 share no structure), so it is large. Independent
        // wins. This actually EXERCISES the indep path (asserted below: a
        // (mate 2, Headers) entry and NO (mate 2, HeaderDelta) entry).
        //
        // 1000 reads: at this scale the columnar per-column BSC framing overhead is
        // amortized and the monotonic numeric columns crush, so independent wins
        // decisively (measured ~1.6 KB indep vs ~2.4 KB delta; at only ~200 reads
        // the columnar overhead dominates and delta narrowly wins — hence 1000).
        let mut r1 = Vec::new();
        let mut r2 = Vec::new();
        for i in 0..1000usize {
            let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT".to_string();
            let qual = "IIHHJJBBCCDDEEFFGGAABBCCDDEEFFGG".to_string();
            r1.push((format!("read{i}/1"), seq.clone(), qual.clone()));
            // Illumina-style: instrument:run:flowcell:lane:tile:x:y — monotonic
            // numeric tail columns compress to almost nothing under columnar delta.
            r2.push((
                format!(
                    "M00001:1:000000000-ABCDE:1:{}:{}:{} 2:N:0:1",
                    1101 + i / 50,
                    1000 + (i % 50) * 17,
                    2000 + i * 13
                ),
                seq,
                qual,
            ));
        }
        let footer = roundtrip(&r1, &r2, 1000, None);
        let has_indep = footer
            .entries
            .iter()
            .any(|e| e.mate == 2 && e.role == StreamRole::Headers);
        let has_delta = footer
            .entries
            .iter()
            .any(|e| e.mate == 2 && e.role == StreamRole::HeaderDelta);
        assert!(
            has_indep && !has_delta,
            "expected R2 header INDEPENDENT path (Headers, not HeaderDelta); roles: {:?}",
            footer
                .entries
                .iter()
                .map(|e| (e.mate, e.role))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn v5_roundtrip_multi_chunk() {
        // 10 reads, chunk_records = 3 → ceil(10/3) = 4 chunks. Similar ids for
        // determinism (delta path). Asserts byte-identical roundtrip ACROSS chunks.
        let (r1, r2) = delta_pairs(10);
        let footer = roundtrip(&r1, &r2, 3, None);
        assert!(
            footer.num_chunks >= 2,
            "expected multiple chunks, got {}",
            footer.num_chunks
        );
        assert_eq!(footer.num_chunks, 4, "10 reads / 3 per chunk = 4 chunks");
    }

    #[test]
    fn v5_roundtrip_fqzcomp_quality() {
        // 50 reads with varied quality strings; explicit fqzcomp quality codec.
        let mut r1 = Vec::new();
        let mut r2 = Vec::new();
        for i in 0..50usize {
            let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT".to_string();
            // Vary quality per read so fqzcomp has real signal.
            let qual: String = (0..32)
                .map(|j| (b'!' + ((i * 3 + j * 5) % 40) as u8) as char)
                .collect();
            r1.push((format!("read{i}/1"), seq.clone(), qual.clone()));
            r2.push((format!("read{i}/2"), seq, qual));
        }
        let footer = roundtrip(&r1, &r2, 1000, Some(QualityCompressor::Fqzcomp));
        let r1_qual_fqz = footer.entries.iter().any(|e| {
            e.mate == 1
                && e.role == StreamRole::Qual
                && e.codec == crate::compression::codec_ids::CODEC_FQZCOMP
        });
        assert!(
            r1_qual_fqz,
            "expected (mate 1, Qual) entry with fqzcomp codec; entries: {:?}",
            footer
                .entries
                .iter()
                .map(|e| (e.mate, e.role, e.codec))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn v5_span_tamper_rejected() {
        // Build a valid archive, then corrupt entry 0's `length` to a huge value,
        // recompute the footer CRC so the body passes the CRC gate, and assert the
        // span-validation guard in `read_v5_footer` rejects it with "out of payload
        // region".
        let dir = TempDir::new().unwrap();
        let r1_path = dir.path().join("r1.fastq");
        let r2_path = dir.path().join("r2.fastq");
        let archive = dir.path().join("paired.qz");

        let (r1, r2) = delta_pairs(20);
        std::fs::write(&r1_path, fastq_string(&r1)).unwrap();
        std::fs::write(&r2_path, fastq_string(&r2)).unwrap();

        let mut cfg = CompressConfig {
            input: vec![r1_path, r2_path],
            output: archive.clone(),
            working_dir: dir.path().to_path_buf(),
            threads: 1,
            ..CompressConfig::default()
        };
        cfg.advanced.chunk_records = 1000;
        super::compress_paired_v5(&cfg).unwrap();

        let mut bytes = std::fs::read(&archive).unwrap();
        let len = bytes.len();
        // Locator: last 20 bytes = [footer_len u64][crc u32][magic 8].
        let footer_len = u64::from_le_bytes(bytes[len - 20..len - 12].try_into().unwrap()) as usize;
        let footer_start = len - 20 - footer_len;
        // Footer body: [num_reads u64][num_chunks u32][n_entries u32] = 16-byte
        // prefix, then ENTRY_BYTES-sized entries. Within an entry, `length` follows
        // chunk_index(4)+role(1)+mate(1)+codec(1)+offset(8) = 15 bytes, and is itself
        // 8 bytes followed by record_count(8) — i.e. at ENTRY_BYTES-8-8 = 15.
        let entry0 = footer_start + 16; // FOOTER_PREFIX (num_reads+num_chunks+n_entries)
        let length_off = entry0 + (ENTRY_BYTES - 8 - 8);
        // A large-but-non-overflowing length: enough that offset+length runs past
        // footer_start (so the span check fires) WITHOUT overflowing u64 (which
        // would trip the earlier overflow guard instead of the region guard).
        bytes[length_off..length_off + 8].copy_from_slice(&(1u64 << 40).to_le_bytes());

        // Re-seal the footer (CRC32 over [footer_start .. len-20]) so the body passes
        // the CRC gate and the decoder reaches the span check.
        let crc = footer_crc32(&bytes[footer_start..len - 20]);
        bytes[len - 12..len - 8].copy_from_slice(&crc.to_le_bytes());

        let bad = dir.path().join("bad.qz");
        std::fs::write(&bad, &bytes).unwrap();

        let dcfg = DecompressConfig {
            input: bad,
            output: vec![dir.path().join("bad_out")],
            working_dir: dir.path().to_path_buf(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 6,
            force: true,
        };
        let res = super::decompress_paired_v5(&dcfg);
        let err = res.expect_err("tampered entry length must be rejected");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("out of payload region"),
            "expected span-validation error 'out of payload region', got: {msg}"
        );
    }
}

#[cfg(test)]
mod columnar_header_decode_tests {
    use super::*;

    fn sample_headers(n: usize) -> Vec<Vec<u8>> {
        (0..n)
            .map(|i| format!("@SRR{:08}.{} read_{}/1", i, i * 7 + 3, i).into_bytes())
            .collect()
    }

    /// A single-block columnar stream (the common case) round-trips unchanged.
    #[test]
    fn decode_columnar_headers_single_block_roundtrip() {
        let hdrs = sample_headers(50);
        let refs: Vec<&[u8]> = hdrs.iter().map(|h| h.as_slice()).collect();
        let blocks =
            crate::compression::compress_impl::columnar_blocks_capped(&refs, usize::MAX).unwrap();
        assert_eq!(blocks.len(), 1, "generous cap → single block");
        let mut buf = Vec::new();
        streams::write_block_stream(&mut buf, &blocks);
        let got = decode_columnar_headers(&buf, hdrs.len() as u64).unwrap();
        assert_eq!(got, hdrs);
    }

    /// A columnar stream split into several capped sub-blocks decodes by
    /// concatenating each block's ids in order — covers the >64 MiB-blob fix
    /// without needing a 64 MiB input (the cap is lowered instead).
    #[test]
    fn decode_columnar_headers_multiblock_roundtrip() {
        let hdrs = sample_headers(300);
        let refs: Vec<&[u8]> = hdrs.iter().map(|h| h.as_slice()).collect();
        let full_len = crate::compression::compress_impl::columnar_blocks_capped(&refs, usize::MAX)
            .unwrap()[0]
            .1
            .len();
        let blocks =
            crate::compression::compress_impl::columnar_blocks_capped(&refs, full_len / 4).unwrap();
        assert!(blocks.len() > 1, "low cap must produce multiple blocks");
        let mut buf = Vec::new();
        streams::write_block_stream(&mut buf, &blocks);
        let got = decode_columnar_headers(&buf, hdrs.len() as u64).unwrap();
        assert_eq!(
            got, hdrs,
            "multi-block decode reconstructs all ids in order"
        );
    }

    /// The directory's record_count must match the summed block counts.
    #[test]
    fn decode_columnar_headers_rejects_record_count_mismatch() {
        let hdrs = sample_headers(20);
        let refs: Vec<&[u8]> = hdrs.iter().map(|h| h.as_slice()).collect();
        let blocks =
            crate::compression::compress_impl::columnar_blocks_capped(&refs, usize::MAX).unwrap();
        let mut buf = Vec::new();
        streams::write_block_stream(&mut buf, &blocks);
        assert!(
            decode_columnar_headers(&buf, 19).is_err(),
            "wrong record_count must be rejected"
        );
    }
}
