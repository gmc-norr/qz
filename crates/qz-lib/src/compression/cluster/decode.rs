//! Decoder for the order-drop cluster mode (v5 `archive_type` = 3).
//!
//! Entry point: [`decompress_cluster`].
//!
//! Output is a *permutation* of the input (set-lossless): reads are emitted in
//! cluster order, NOT the original input order.  Flip bits are used to
//! reverse-complement-undo any RC-canonicalised reads back to their original
//! orientation before writing.

use anyhow::{Result, bail, ensure};

use crate::cli::DecompressConfig;
use crate::compression::bsc;
use crate::compression::chunk_directory::{
    ChunkDirEntry, ChunkDirectory, GLOBAL_SENTINEL, StreamRole,
};
use crate::compression::codecs::decompress_qualities_fqzcomp;
use crate::compression::decompress_impl::BSC_MAX_BLOCK_LEN;
use crate::compression::dna_utils::{read_varint, reverse_complement_canonical};
use crate::compression::header_col::decompress_headers_columnar;

/// 2.5 GiB cap on a single `ClusteredSeq` decompressed frame.
pub(crate) const CLUSTER_SEQ_CAP: usize = 2_560 << 20;

// ── ClusterMeta ────────────────────────────────────────────────────────────────

/// Decoded form of the 28-byte global ClusterMeta block. `version`/`kind`/`level`
/// are part of the on-disk wire format and are validated/retained for
/// forward-compatibility and diagnostics; only `nreads`/`checksum` are consumed by
/// the current validator and verifier.
#[allow(dead_code)]
pub(crate) struct ClusterMeta {
    pub version: u8,
    /// Archive kind: 0 = single-end, 1 = paired-end.
    pub kind: u8,
    pub level: u8,
    pub nreads: u64,
    pub checksum: u128,
}

impl ClusterMeta {
    /// Parse from the 28 decoded meta bytes.
    fn from_bytes(b: &[u8]) -> Result<Self> {
        ensure!(b.len() == 28, "ClusterMeta: expected 28 bytes, got {}", b.len());
        ensure!(b[0] == 1, "ClusterMeta: unsupported version {} (expected 1)", b[0]);
        let nreads = u64::from_le_bytes(b[4..12].try_into().unwrap());
        let checksum = u128::from_le_bytes(b[12..28].try_into().unwrap());
        Ok(ClusterMeta {
            version: b[0],
            kind: b[1],
            level: b[2],
            nreads,
            checksum,
        })
    }
}

// ── validate_cluster_directory ─────────────────────────────────────────────────

/// A cluster archive is paired iff any directory entry is tagged mate 2 (R2).
pub(crate) fn is_paired(dir: &ChunkDirectory) -> bool {
    dir.entries.iter().any(|e| e.mate == 2)
}

/// Validate one (chunk, mate)'s role multiset: exactly 1 `ClusteredSeq`, exactly
/// 1 `StrandBits`, ≥1 `Qual`, ≥1 `Headers`, nothing else (frame-size caps are
/// applied to each entry as it is encountered).  Returns that mate's
/// `ClusteredSeq.record_count` so the caller can accumulate and compare to
/// `meta.nreads`.
///
/// Re-scans `dir.entries` per (chunk, mate) — O(num_chunks · entries) overall —
/// which is fine: a validator runs once and cluster directories are small
/// (≈4 roles × mates × chunks + 1).
fn validate_chunk_mate_roles(dir: &ChunkDirectory, chunk: u32, mate: u8) -> Result<u64> {
    let mut seq_count: u32 = 0;
    let mut strand_count: u32 = 0;
    let mut qual_count: u32 = 0;
    let mut header_count: u32 = 0;
    let mut seq_rec: u64 = 0;

    for e in dir.entries.iter().filter(|e| e.chunk_index == chunk && e.mate == mate) {
        match e.role {
            StreamRole::ClusteredSeq => {
                seq_count += 1;
                seq_rec += e.record_count;
                ensure!(
                    e.length <= CLUSTER_SEQ_CAP as u64,
                    "cluster directory: ClusteredSeq entry length {} exceeds cap {}",
                    e.length,
                    CLUSTER_SEQ_CAP
                );
            }
            StreamRole::StrandBits => {
                strand_count += 1;
                ensure!(
                    e.length <= BSC_MAX_BLOCK_LEN as u64,
                    "cluster directory: StrandBits entry length {} exceeds cap {}",
                    e.length,
                    BSC_MAX_BLOCK_LEN
                );
            }
            StreamRole::Qual => {
                qual_count += 1;
                ensure!(
                    e.length <= BSC_MAX_BLOCK_LEN as u64,
                    "cluster directory: Qual entry length {} exceeds cap {}",
                    e.length,
                    BSC_MAX_BLOCK_LEN
                );
            }
            StreamRole::Headers => {
                header_count += 1;
                ensure!(
                    e.length <= BSC_MAX_BLOCK_LEN as u64,
                    "cluster directory: Headers entry length {} exceeds cap {}",
                    e.length,
                    BSC_MAX_BLOCK_LEN
                );
            }
            other => bail!(
                "cluster directory: chunk {} mate {} has unexpected role {:?}",
                chunk,
                mate,
                other
            ),
        }
    }

    ensure!(
        seq_count == 1,
        "cluster directory: chunk {} mate {} has {} ClusteredSeq entries (expected 1)",
        chunk,
        mate,
        seq_count
    );
    ensure!(
        strand_count == 1,
        "cluster directory: chunk {} mate {} has {} StrandBits entries (expected 1)",
        chunk,
        mate,
        strand_count
    );
    ensure!(
        qual_count >= 1,
        "cluster directory: chunk {} mate {} has no Qual entries",
        chunk,
        mate
    );
    ensure!(
        header_count >= 1,
        "cluster directory: chunk {} mate {} has no Headers entries",
        chunk,
        mate
    );

    Ok(seq_rec)
}

/// Structural validator for a cluster archive directory.
///
/// Enforced:
/// (a) Exactly one global (`GLOBAL_SENTINEL`) entry of role `ClusterMeta`; no
///     other global entries.
/// (b) `dir.num_reads == meta.nreads`.
/// (c) `meta.kind` is consistent with mate-2 presence in the directory.
/// (d) **Single-end** (no mate-2 entries): each chunk `0..dir.num_chunks` (all
///     mate 0) has exactly one `ClusteredSeq`, exactly one `StrandBits`, ≥1
///     `Qual`, ≥1 `Headers`, and no other roles; no mate other than 0 allowed.
///     **Paired** (any mate-2 entry): each chunk `0..dir.num_chunks` has both
///     mate 1 and mate 2 sub-sets, each with the same four-role requirement;
///     no mate-0 chunk entries (only the global ClusterMeta carries mate 0).
/// (e) Sum of `ClusteredSeq.record_count` across all chunks equals `meta.nreads`
///     (for each mate independently in paired mode).
/// (f) Frame-size caps: `ClusteredSeq` entries ≤ `CLUSTER_SEQ_CAP`; all
///     others ≤ `BSC_MAX_BLOCK_LEN`.
pub(crate) fn validate_cluster_directory(
    dir: &ChunkDirectory,
    meta: &ClusterMeta,
) -> Result<()> {
    // ── (a) exactly one global entry of role ClusterMeta ─────────────────────
    let global_entries: Vec<&ChunkDirEntry> = dir
        .entries
        .iter()
        .filter(|e| e.chunk_index == GLOBAL_SENTINEL)
        .collect();
    ensure!(
        global_entries.len() == 1,
        "cluster directory: expected exactly 1 global entry, found {}",
        global_entries.len()
    );
    ensure!(
        global_entries[0].role == StreamRole::ClusterMeta,
        "cluster directory: global entry has wrong role {:?}, expected ClusterMeta",
        global_entries[0].role
    );

    // ── (b) num_reads matches meta.nreads ─────────────────────────────────────
    ensure!(
        dir.num_reads == meta.nreads,
        "cluster directory: dir.num_reads {} != meta.nreads {}",
        dir.num_reads,
        meta.nreads
    );

    // ── (c) archive_kind cross-check ──────────────────────────────────────────
    // The on-disk kind byte (0 = single-end, 1 = paired) must agree with the
    // physical layout (paired ⇔ mate-2 entries present), so a forged/garbled kind
    // can't steer decode down the wrong path.
    let paired = is_paired(dir);
    ensure!(
        (meta.kind == 1) == paired,
        "cluster directory: ClusterMeta archive_kind ({}) disagrees with mate-2 presence ({})",
        meta.kind,
        paired
    );

    let nchunks = dir.num_chunks as usize;

    if paired {
        // ── Paired: mates {1, 2}; no mate-0 chunk entries ────────────────────
        // Reject any non-global entry with an unexpected mate.
        for e in dir.entries.iter().filter(|e| e.chunk_index != GLOBAL_SENTINEL) {
            let ci = e.chunk_index as usize;
            ensure!(
                ci < nchunks,
                "cluster directory: chunk_index {} out of range [0, {})",
                ci,
                nchunks
            );
            ensure!(
                e.mate == 1 || e.mate == 2,
                "cluster directory (paired): chunk {} has mate {} (expected 1 or 2)",
                ci,
                e.mate
            );
        }

        // Per-(chunk, mate) role check + record_count accumulation.
        let mut r1_sum: u64 = 0;
        let mut r2_sum: u64 = 0;
        for ci in 0..dir.num_chunks {
            r1_sum += validate_chunk_mate_roles(dir, ci, 1)?;
            r2_sum += validate_chunk_mate_roles(dir, ci, 2)?;
        }
        ensure!(
            r1_sum == meta.nreads,
            "cluster directory: sum of mate-1 ClusteredSeq record_counts {} != meta.nreads {}",
            r1_sum,
            meta.nreads
        );
        ensure!(
            r2_sum == meta.nreads,
            "cluster directory: sum of mate-2 ClusteredSeq record_counts {} != meta.nreads {}",
            r2_sum,
            meta.nreads
        );
    } else {
        // ── Single-end: all chunk entries must be mate 0 ──────────────────────
        for e in dir.entries.iter().filter(|e| e.chunk_index != GLOBAL_SENTINEL) {
            let ci = e.chunk_index as usize;
            ensure!(
                ci < nchunks,
                "cluster directory: chunk_index {} out of range [0, {})",
                ci,
                nchunks
            );
            ensure!(
                e.mate == 0,
                "cluster directory: chunk {} has mate {} (expected 0)",
                ci,
                e.mate
            );
        }

        let mut seq_rec_sum: u64 = 0;
        for ci in 0..dir.num_chunks {
            seq_rec_sum += validate_chunk_mate_roles(dir, ci, 0)?;
        }
        ensure!(
            seq_rec_sum == meta.nreads,
            "cluster directory: sum of ClusteredSeq record_counts {} != meta.nreads {}",
            seq_rec_sum,
            meta.nreads
        );
    }

    Ok(())
}

// ── SeekReader ─────────────────────────────────────────────────────────────────

/// A file open for positioned reads (no shared seek cursor).
struct SeekReader {
    file: std::fs::File,
    len: u64,
}

impl SeekReader {
    fn open(path: &std::path::Path) -> Result<Self> {
        let file = std::fs::File::open(path)
            .map_err(|e| anyhow::anyhow!("open cluster archive {}: {e}", path.display()))?;
        let len = file
            .metadata()
            .map_err(|e| anyhow::anyhow!("stat cluster archive {}: {e}", path.display()))?
            .len();
        Ok(SeekReader { file, len })
    }

    /// Read exactly `length` bytes at absolute `offset`.
    fn read_at(&self, offset: u64, length: u64) -> Result<Vec<u8>> {
        let end = offset
            .checked_add(length)
            .ok_or_else(|| anyhow::anyhow!("entry offset+length overflow"))?;
        if end > self.len {
            bail!(
                "cluster entry [{offset}, {end}) out of bounds (file len {})",
                self.len
            );
        }
        let len = usize::try_from(length)
            .map_err(|_| anyhow::anyhow!("entry length {length} exceeds usize"))?;
        let mut buf = vec![0u8; len];
        #[cfg(unix)]
        {
            use std::os::unix::fs::FileExt;
            self.file.read_exact_at(&mut buf, offset)?;
        }
        #[cfg(not(unix))]
        {
            use std::io::{Read, Seek, SeekFrom};
            // For non-unix: use seek (not concurrent-safe, but cluster is sequential)
            let mut f = &self.file;
            f.seek(SeekFrom::Start(offset))?;
            f.read_exact(&mut buf)?;
        }
        Ok(buf)
    }
}

// ── Frame read + CRC validate ─────────────────────────────────────────────────

/// Read the bare frame at `entry.offset`/`entry.length`, CRC-validate it, and
/// return the payload bytes (everything after the 12-byte frame header).
fn read_frame(reader: &SeekReader, entry: &ChunkDirEntry) -> Result<Vec<u8>> {
    let frame = reader.read_at(entry.offset, entry.length)?;
    ensure!(
        frame.len() >= 12,
        "cluster frame too short ({} bytes)",
        frame.len()
    );
    let block_len = u32::from_le_bytes(frame[0..4].try_into().unwrap()) as usize;
    let record_count = u32::from_le_bytes(frame[4..8].try_into().unwrap());
    let crc = u32::from_le_bytes(frame[8..12].try_into().unwrap());
    let payload = &frame[12..];
    ensure!(
        payload.len() == block_len,
        "cluster frame: block_len {} != actual payload len {}",
        block_len,
        payload.len()
    );
    // The inline frame record_count and the footer entry.record_count are written
    // from the same value at encode (`emit_block`); cross-check so a forged footer
    // can't drive the per-record loop to over/under-read the stream (single-end and
    // paired do the equivalent check). Both are attacker-controlled, so this guards
    // consistency, not magnitude — the magnitude bound lives at the alloc site.
    ensure!(
        record_count as u64 == entry.record_count,
        "cluster frame: inline record_count {} != footer entry record_count {}",
        record_count,
        entry.record_count
    );
    bsc::verify_block_frame_crc(crc, record_count, payload)?;
    Ok(payload.to_vec())
}

// ── Varint stream parser ───────────────────────────────────────────────────────

/// Parse the varint-prefixed per-record quality stream returned by
/// `decompress_qualities_fqzcomp`.  Returns one `Vec<u8>` per record.
fn parse_fqz_varint_stream(data: &[u8], record_count: usize) -> Result<Vec<Vec<u8>>> {
    // `record_count` comes from the attacker-controlled footer; every record
    // consumes at least one varint length byte, so the output can never exceed
    // `data.len() + 1` entries. Clamp the pre-allocation to that (the per-record
    // loop below still errors cleanly via `read_varint` if the count overruns the
    // data) so a forged huge count can't trigger an allocation-failure abort —
    // matching the paired path's `.min(bytes.len() + 1)` clamp.
    let mut out = Vec::with_capacity(record_count.min(data.len() + 1));
    let mut pos = 0usize;
    for i in 0..record_count {
        let len = read_varint(data, &mut pos)
            .ok_or_else(|| anyhow::anyhow!("fqz varint stream truncated at record {i}"))?
            as usize;
        let end = pos
            .checked_add(len)
            .ok_or_else(|| anyhow::anyhow!("fqz varint length overflow at record {i}"))?;
        ensure!(
            end <= data.len(),
            "fqz varint stream: record {i} len {len} overruns data (pos {pos}, total {})",
            data.len()
        );
        out.push(data[pos..end].to_vec());
        pos = end;
    }
    Ok(out)
}

// ── decompress_cluster ─────────────────────────────────────────────────────────

/// Decompress a cluster archive to a single FASTQ output.
///
/// Output is a permutation of the input (set-lossless): records are emitted in
/// cluster order, with RC-canonicalized reads restored to their original
/// orientation via the stored flip bits.
pub(crate) fn decompress_cluster(args: &DecompressConfig) -> Result<()> {
    // ── Open + parse front header ─────────────────────────────────────────────
    let input = &args.input;
    {
        use std::io::Read;
        let mut hbuf = [0u8; 64];
        let mut file = std::fs::File::open(input)
            .map_err(|e| anyhow::anyhow!("open cluster archive {}: {e}", input.display()))?;
        let n = file.read(&mut hbuf)?;
        if n < crate::compression::V2_PREFIX_SIZE {
            bail!("cluster v5 header truncated");
        }
        let hdr = crate::compression::FixedHeader::parse_v5(&hbuf[..n])?;
        if hdr.archive_type != 3 {
            bail!(
                "not a cluster v5 archive (archive_type {})",
                hdr.archive_type
            );
        }
    }

    // ── Read footer ───────────────────────────────────────────────────────────
    let header_end = {
        use std::io::Read;
        let mut hbuf = [0u8; 8];
        let mut file = std::fs::File::open(input)?;
        file.read_exact(&mut hbuf)?;
        u32::from_le_bytes(hbuf[4..8].try_into().unwrap()) as u64
    };
    let dir = crate::compression::chunk_directory::read_v5_footer(input, header_end)?;

    // ── Open seekable reader ──────────────────────────────────────────────────
    let reader = SeekReader::open(input)?;

    // ── Decode ClusterMeta global ─────────────────────────────────────────────
    let meta_entry = dir
        .entries
        .iter()
        .find(|e| e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::ClusterMeta)
        .ok_or_else(|| anyhow::anyhow!("cluster archive: missing ClusterMeta global entry"))?;
    let meta_payload = read_frame(&reader, meta_entry)?;
    // ClusterMeta is small fixed-shape metadata; cap the fallback decode at one block to
    // stop a stacked-BSC bomb in the (untrusted) global from forcing a huge allocation.
    let meta_raw = bsc::decompress_parallel_max(
        &meta_payload,
        crate::compression::decompress_impl::BSC_MAX_BLOCK_LEN,
    )?;
    let meta = ClusterMeta::from_bytes(&meta_raw)?;

    // ── Validate directory ────────────────────────────────────────────────────
    validate_cluster_directory(&dir, &meta)?;

    // ── Paired path: two-output atomic decode ─────────────────────────────────
    if is_paired(&dir) {
        return decompress_cluster_paired(args, &dir, input);
    }

    // Single-end cluster: exactly one output (mirrors the type-0/1/2/4 guard in the
    // dispatcher and the paired-cluster `== 2` check). Reject extra `-o` LOUDLY rather
    // than silently dropping `output[1..]`.
    if args.output.len() != 1 {
        anyhow::bail!(
            "single-end cluster archive decompresses to a single -o argument; got {} -o values",
            args.output.len()
        );
    }

    // ── Decode all chunks via the atomic output sink ──────────────────────────
    crate::compression::decompress_impl::decompress_streaming_bsc_dispatch(
        args,
        std::time::Instant::now(),
        meta.nreads as usize,
        |out| {
            // Re-open a seekable reader inside the closure (decompress_streaming_bsc_dispatch
            // owns the output side; we own the input side).
            let reader = SeekReader::open(input)?;

            for ci in 0..dir.num_chunks {
                decode_chunk(&reader, &dir, ci, out)?;
            }
            Ok(())
        },
    )
}

/// Decode one chunk, reconstructing each record's `(header, original_seq, qual)` and
/// passing them to `on_record`.  This is the shared decode core used by both
/// `decompress_cluster` (writes FASTQ) and `verify_cluster` (folds checksum).
fn decode_chunk_records<F>(
    reader: &SeekReader,
    dir: &ChunkDirectory,
    chunk_index: u32,
    mate: u8,
    mut on_record: F,
) -> Result<()>
where
    F: FnMut(&[u8], &[u8], &[u8]) -> Result<()>,
{
    // Collect this chunk's entries, sorted by offset (preserves block order).
    let mut chunk_entries: Vec<&ChunkDirEntry> = dir
        .entries
        .iter()
        .filter(|e| e.chunk_index == chunk_index && e.mate == mate)
        .collect();
    chunk_entries.sort_by_key(|e| e.offset);

    // ── Qual blocks → per-record qual Vecs ───────────────────────────────────
    let qual_entries: Vec<&ChunkDirEntry> = chunk_entries
        .iter()
        .copied()
        .filter(|e| e.role == StreamRole::Qual)
        .collect();
    let mut quals: Vec<Vec<u8>> = Vec::new();
    for qe in &qual_entries {
        let payload = read_frame(reader, qe)?;
        let varint_stream = decompress_qualities_fqzcomp(&payload)?;
        let mut block_quals = parse_fqz_varint_stream(&varint_stream, qe.record_count as usize)?;
        quals.append(&mut block_quals);
    }

    // ── Headers blocks → per-record header Vecs ───────────────────────────────
    let header_entries: Vec<&ChunkDirEntry> = chunk_entries
        .iter()
        .copied()
        .filter(|e| e.role == StreamRole::Headers)
        .collect();
    let mut headers: Vec<Vec<u8>> = Vec::new();
    for he in &header_entries {
        let payload = read_frame(reader, he)?;
        ensure!(
            payload.len() >= 4,
            "cluster Headers block too short ({} bytes)",
            payload.len()
        );
        let num_reads = u32::from_le_bytes(payload[0..4].try_into().unwrap()) as usize;
        let blob = &payload[4..];
        let mut block_headers = decompress_headers_columnar(blob, num_reads)?;
        headers.append(&mut block_headers);
    }

    // ── StrandBits → flip flags ────────────────────────────────────────────────
    let strand_entry = chunk_entries
        .iter()
        .find(|e| e.role == StreamRole::StrandBits)
        .ok_or_else(|| anyhow::anyhow!("cluster chunk {chunk_index}: missing StrandBits entry"))?;
    let strand_payload = read_frame(reader, strand_entry)?;
    // StrandBits is a 1-bit-per-record bitvec (ceil(nrec/8) bytes, always far under one
    // block); cap the untrusted decode at the block max to stop a stacked-BSC bomb.
    let bits = bsc::decompress_parallel_max(
        &strand_payload,
        crate::compression::decompress_impl::BSC_MAX_BLOCK_LEN,
    )?;

    // ── ClusteredSeq → per-record canon seq Vecs ─────────────────────────────
    let seq_entry = chunk_entries
        .iter()
        .find(|e| e.role == StreamRole::ClusteredSeq)
        .ok_or_else(|| anyhow::anyhow!("cluster chunk {chunk_index}: missing ClusteredSeq entry"))?;

    let chunk_record_count = seq_entry.record_count as usize;

    // Sanity: headers and quals should match
    ensure!(
        headers.len() == chunk_record_count,
        "cluster chunk {chunk_index}: headers.len() {} != chunk_record_count {}",
        headers.len(),
        chunk_record_count
    );
    ensure!(
        quals.len() == chunk_record_count,
        "cluster chunk {chunk_index}: quals.len() {} != chunk_record_count {}",
        quals.len(),
        chunk_record_count
    );

    // expected_len = sum of qual lengths (qual len == seq len for each record)
    let expected_seq_len: usize = quals.iter().map(|q| q.len()).sum();
    let seq_payload = read_frame(reader, seq_entry)?;
    let flat_seqs = crate::compression::cluster::zstd_codec::decompress_seq_zstd_long(
        &seq_payload,
        expected_seq_len,
    )?;

    // Split flat_seqs into per-record slices using qual lengths, restore orientation,
    // and call the per-record callback.
    let mut seq_pos = 0usize;
    for i in 0..chunk_record_count {
        let seq_len = quals[i].len();
        let end = seq_pos + seq_len;
        ensure!(
            end <= flat_seqs.len(),
            "cluster chunk {chunk_index}: seq underrun at record {i} (pos {seq_pos}, end {end}, total {})",
            flat_seqs.len()
        );
        let canon_seq = &flat_seqs[seq_pos..end];
        seq_pos = end;

        // Unflipped reads (the common case) pass the canonical slice straight through —
        // no per-record allocation. Only flipped reads materialize a reverse complement.
        let flip = (bits.get(i >> 3).copied().unwrap_or(0) >> (i & 7)) & 1 == 1;
        if flip {
            let seq = reverse_complement_canonical(canon_seq);
            on_record(&headers[i], &seq, &quals[i])?;
        } else {
            on_record(&headers[i], canon_seq, &quals[i])?;
        }
    }

    Ok(())
}

/// Write one FASTQ record (`header\nseq\n+\nqual\n`) to `out`.
/// `header` already includes the `@` prefix.
fn write_fastq_record(out: &mut dyn std::io::Write, header: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    out.write_all(header)?;
    out.write_all(b"\n")?;
    out.write_all(seq)?;
    out.write_all(b"\n+\n")?;
    out.write_all(qual)?;
    out.write_all(b"\n")?;
    Ok(())
}

/// Decode one chunk and write its records to `out`.
fn decode_chunk(
    reader: &SeekReader,
    dir: &ChunkDirectory,
    chunk_index: u32,
    out: &mut dyn std::io::Write,
) -> Result<()> {
    decode_chunk_records(reader, dir, chunk_index, 0, |header, seq, qual| {
        write_fastq_record(out, header, seq, qual)
    })
}

/// Paired cluster decode: two outputs (R1 = mate 1, R2 = mate 2), atomic two-file
/// publish mirroring `paired::decompress_paired_v5`. Pair order is preserved because
/// the encoder emitted both mates in the same clustered pair-order.
fn decompress_cluster_paired(
    args: &crate::cli::DecompressConfig,
    dir: &ChunkDirectory,
    input: &std::path::Path,
) -> Result<()> {
    use crate::compression::paired::{publish_pair, with_suffix, TmpPair, TMP_NONCE};
    use std::io::{BufWriter, Write as _};
    use std::sync::atomic::Ordering;

    ensure!(
        args.output.len() == 2,
        "paired cluster decode needs two output paths (got {})",
        args.output.len()
    );
    let out1 = &args.output[0];
    let out2 = &args.output[1];
    for o in [out1, out2] {
        if crate::cli::is_stdio_path(o) {
            bail!("paired cluster decode needs named output files (no stdout)");
        }
    }
    // Gzip output is not supported for paired cluster (matches paired mode); plain only.
    if args.gzipped {
        bail!("gzipped output not supported for paired cluster archives (writes plain FASTQ)");
    }

    // Atomic two-output write: O_EXCL temps + drop guard + publish_pair rename.
    // (The shared decompress() force-check already refused pre-existing out1/out2.)
    let pid = std::process::id();
    let nonce = TMP_NONCE.fetch_add(1, Ordering::Relaxed);
    let tmp1 = with_suffix(out1, &format!(".{pid}.{nonce}.tmp"));
    let tmp2 = with_suffix(out2, &format!(".{pid}.{nonce}.tmp"));
    let mut guard = TmpPair { a: tmp1.clone(), b: tmp2.clone(), armed: true };
    let mut w1 = BufWriter::new(crate::compression::create_new_for_write(&tmp1)?);
    let mut w2 = BufWriter::new(crate::compression::create_new_for_write(&tmp2)?);

    let reader = SeekReader::open(input)?;
    for ci in 0..dir.num_chunks {
        decode_chunk_records(&reader, dir, ci, 1, |h, s, q| write_fastq_record(&mut w1, h, s, q))?;
        decode_chunk_records(&reader, dir, ci, 2, |h, s, q| write_fastq_record(&mut w2, h, s, q))?;
    }

    w1.flush()?;
    w2.flush()?;
    drop(w1);
    drop(w2);
    // Publish FIRST, disarm only after success: if publish_pair fails mid-rename the
    // guard's Drop still cleans up both temps (no orphaned .tmp files left behind).
    publish_pair(&tmp1, &tmp2, out1, out2, nonce)?;
    guard.disarm();
    Ok(())
}

// ── verify_cluster ─────────────────────────────────────────────────────────────

/// Open a cluster archive, parse its front header and footer, and return the
/// seekable reader + parsed directory + ClusterMeta.
fn open_cluster_archive(
    input: &std::path::Path,
) -> Result<(SeekReader, crate::compression::chunk_directory::ChunkDirectory, ClusterMeta)> {
    use std::io::Read;

    // Parse front header
    {
        let mut hbuf = [0u8; 64];
        let mut file = std::fs::File::open(input)
            .map_err(|e| anyhow::anyhow!("open cluster archive {}: {e}", input.display()))?;
        let n = file.read(&mut hbuf)?;
        if n < crate::compression::V2_PREFIX_SIZE {
            bail!("cluster v5 header truncated");
        }
        let hdr = crate::compression::FixedHeader::parse_v5(&hbuf[..n])?;
        if hdr.archive_type != 3 {
            bail!(
                "not a cluster v5 archive (archive_type {})",
                hdr.archive_type
            );
        }
    }

    // Read header_end from front header bytes 4..8
    let header_end = {
        let mut hbuf = [0u8; 8];
        let mut file = std::fs::File::open(input)?;
        file.read_exact(&mut hbuf)?;
        u32::from_le_bytes(hbuf[4..8].try_into().unwrap()) as u64
    };
    let dir = crate::compression::chunk_directory::read_v5_footer(input, header_end)?;
    let reader = SeekReader::open(input)?;

    // Decode ClusterMeta global
    let meta_entry = dir
        .entries
        .iter()
        .find(|e| e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::ClusterMeta)
        .ok_or_else(|| anyhow::anyhow!("cluster archive: missing ClusterMeta global entry"))?;
    let meta_payload = read_frame(&reader, meta_entry)?;
    // ClusterMeta is small fixed-shape metadata; cap the fallback decode at one block to
    // stop a stacked-BSC bomb in the (untrusted) global from forcing a huge allocation.
    let meta_raw = bsc::decompress_parallel_max(
        &meta_payload,
        crate::compression::decompress_impl::BSC_MAX_BLOCK_LEN,
    )?;
    let meta = ClusterMeta::from_bytes(&meta_raw)?;

    Ok((reader, dir, meta))
}

/// Verify a cluster archive — either fast (per-block CRC32 only) or deep
/// (full decode + multiset checksum comparison).
pub(crate) fn verify_cluster(
    input: &std::path::Path,
    fast: bool,
) -> Result<crate::compression::VerifyResult> {
    use crate::cli::{HeaderCompressor, QualityCompressor};
    use crate::compression::{VerifyMode, VerifyResult};

    let start = std::time::Instant::now();

    let (reader, dir, meta) = open_cluster_archive(input)?;
    validate_cluster_directory(&dir, &meta)?;

    if fast {
        // Fast path: walk every directory entry, validate per-block CRC32.
        let mut blocks_verified: u32 = 0;
        let mut total_bytes: u64 = 0;

        for entry in &dir.entries {
            let frame = reader.read_at(entry.offset, entry.length)?;
            ensure!(
                frame.len() >= 12,
                "cluster frame too short ({} bytes) at offset {}",
                frame.len(),
                entry.offset
            );
            let block_len = u32::from_le_bytes(frame[0..4].try_into().unwrap()) as usize;
            let record_count = u32::from_le_bytes(frame[4..8].try_into().unwrap());
            let crc = u32::from_le_bytes(frame[8..12].try_into().unwrap());
            let payload = &frame[12..];
            ensure!(
                payload.len() == block_len,
                "cluster frame: block_len {} != actual payload len {}",
                block_len,
                payload.len()
            );
            // Same inline-vs-footer record_count cross-check `read_frame` does, so
            // fast verify rejects a forged footer count that the decode path catches.
            ensure!(
                record_count as u64 == entry.record_count,
                "cluster frame: inline record_count {} != footer entry record_count {}",
                record_count,
                entry.record_count
            );
            bsc::verify_block_frame_crc(crc, record_count, payload)?;
            blocks_verified += 1;
            total_bytes += entry.length;
        }

        let elapsed_secs = start.elapsed().as_secs_f64();
        return Ok(VerifyResult {
            num_reads: dir.num_reads as usize,
            encoding_type: 11,
            header_compressor: HeaderCompressor::Columnar,
            quality_compressor: QualityCompressor::Fqzcomp,
            headers_compressed_len: 0,
            sequences_compressed_len: 0,
            qualities_compressed_len: 0,
            crc32: 0,
            r2_crc32: None,
            total_bytes,
            blocks_verified,
            mode: VerifyMode::Fast,
            elapsed_secs,
        });
    }

    // Deep path: decode all records and fold the multiset checksum.
    let mut acc: u128 = 0;
    let mut count: u64 = 0;

    if !is_paired(&dir) {
        // Single-end: fold hash_record per record (mate 0).
        for ci in 0..dir.num_chunks {
            decode_chunk_records(&reader, &dir, ci, 0, |header, seq, qual| {
                acc = acc.wrapping_add(crate::compression::cluster::checksum::hash_record(
                    header, seq, qual,
                ));
                count += 1;
                Ok(())
            })?;
        }
    } else {
        // Paired: per chunk, collect ONLY mate-1, then STREAM mate-2 and pair positionally
        // as each record arrives — folding hash_pair without materializing the second
        // mate's chunk. This halves the per-chunk transient vs. holding both mate vectors
        // (the checksum is order-identical: same positional pairing). Still per-chunk (not
        // all-chunks), so RAM stays bounded in read count.
        for ci in 0..dir.num_chunks {
            let mut m1: Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> = Vec::new();
            decode_chunk_records(&reader, &dir, ci, 1, |h, s, q| {
                m1.push((h.to_vec(), s.to_vec(), q.to_vec()));
                Ok(())
            })?;
            let mut i = 0usize;
            decode_chunk_records(&reader, &dir, ci, 2, |h2, s2, q2| {
                ensure!(
                    i < m1.len(),
                    "cluster paired verify: chunk {ci} has more mate-2 records than mate-1"
                );
                let (h1, s1, q1) = &m1[i];
                acc = acc.wrapping_add(crate::compression::cluster::checksum::hash_pair(
                    h1, s1, q1, h2, s2, q2,
                ));
                // `count` is in PAIRS here (one per hash_pair), matching meta.nreads which
                // the encoder stores as n_pairs for paired archives. The single-end branch
                // counts records; both feed the shared `count == meta.nreads` check below.
                count += 1;
                i += 1;
                Ok(())
            })?;
            ensure!(
                i == m1.len(),
                "cluster paired verify: chunk {} mate counts differ ({} vs {})",
                ci,
                m1.len(),
                i
            );
        }
    }

    ensure!(
        count == meta.nreads,
        "cluster verify: read count {count} != {} in metadata",
        meta.nreads
    );
    ensure!(
        acc == meta.checksum,
        "cluster verify: multiset checksum mismatch (archive corrupt or tampered)"
    );

    let elapsed_secs = start.elapsed().as_secs_f64();
    Ok(VerifyResult {
        num_reads: count as usize,
        encoding_type: 11,
        header_compressor: HeaderCompressor::Columnar,
        quality_compressor: QualityCompressor::Fqzcomp,
        headers_compressed_len: 0,
        sequences_compressed_len: 0,
        qualities_compressed_len: 0,
        crc32: acc as u32,
        r2_crc32: None,
        total_bytes: 0,
        blocks_verified: 0,
        mode: VerifyMode::Deep,
        elapsed_secs,
    })
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_fastq(path: &std::path::Path, n: usize, len: usize) {
        let mut f = std::fs::File::create(path).unwrap();
        let bases = b"ACGT";
        for i in 0..n {
            let seq: Vec<u8> = (0..len)
                .map(|j| bases[(i * 13 + j * 7 + (j * j)) % 4])
                .collect();
            writeln!(f, "@read{i} extra:{i}").unwrap();
            f.write_all(&seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            let q: Vec<u8> = (0..len).map(|j| b'!' + ((i + j) % 40) as u8).collect();
            f.write_all(&q).unwrap();
            writeln!(f).unwrap();
        }
    }

    fn multiset_records(p: &std::path::Path) -> Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> {
        let s = std::fs::read(p).unwrap();
        let mut v = vec![];
        let mut lines = s.split(|&b| b == b'\n');
        while let (Some(h), Some(sq), Some(_), Some(q)) =
            (lines.next(), lines.next(), lines.next(), lines.next())
        {
            if h.is_empty() {
                break;
            }
            v.push((h.to_vec(), sq.to_vec(), q.to_vec()));
        }
        v.sort();
        v
    }

    #[test]
    fn cluster_roundtrip_empty_input() {
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("in.fastq");
        std::fs::write(&fq, b"").unwrap();
        let qz = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![fq.clone()], output: qz.clone(), working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        crate::compression::cluster::compress_cluster(&cfg).unwrap();
        let back = dir.path().join("back.fastq");
        let dc = crate::cli::DecompressConfig {
            input: qz, output: vec![back.clone()], working_dir: dir.path().to_path_buf(),
            num_threads: 1, gzipped: false, gzip_level: 6, force: true,
        };
        crate::compression::cluster::decompress_cluster(&dc).unwrap();
        assert_eq!(multiset_records(&fq), multiset_records(&back));
    }

    #[test]
    fn cluster_roundtrip_covers_reverse_complement() {
        use crate::compression::dna_utils::canonical_minimizer_orientation;
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("in.fastq");
        // Strategy: sequences composed entirely of G and T have no A/C s-mers, so
        // their open-syncmer (k=21,s=10,ts=[0]) always selects a k-mer starting with
        // G or T.  In the 2-bit encoding A=0,C=1,G=2,T=3 the RC of G is C (1<2), so
        // rc_hash < fwd_hash => canonical_minimizer_orientation returns true for all
        // such reads.  Mix in some regular A/C reads so the fixture exercises both
        // flip=true and flip=false paths through decode_chunk.
        let mut f = std::fs::File::create(&fq).unwrap();
        use std::io::Write;
        let len = 100usize;
        let mut flipped = 0usize;
        // 200 "GT-only" reads → all flip
        for i in 0..200usize {
            let seq: Vec<u8> = (0..len).map(|j| if (i + j) % 2 == 0 { b'G' } else { b'T' }).collect();
            if canonical_minimizer_orientation(&seq) { flipped += 1; }
            writeln!(f, "@flip_read{i}").unwrap();
            f.write_all(&seq).unwrap(); writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(len)).unwrap();
        }
        // 200 regular reads (ACGT mix) → none flip (verified by the syncmer property)
        let bases = b"ACGT";
        for i in 0..200usize {
            let seq: Vec<u8> = (0..len)
                .map(|j| bases[(i * 13 + j * 7 + (j * j)) % 4])
                .collect();
            // These are NOT expected to flip; just pad the fixture with non-flip reads.
            writeln!(f, "@noflip_read{i}").unwrap();
            f.write_all(&seq).unwrap(); writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "J".repeat(len)).unwrap();
        }
        drop(f);
        assert!(flipped > 0, "test fixture must contain at least one reverse-complemented read; got {flipped}");
        let qz = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![fq.clone()], output: qz.clone(), working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        crate::compression::cluster::compress_cluster(&cfg).unwrap();
        let back = dir.path().join("back.fastq");
        let dc = crate::cli::DecompressConfig {
            input: qz, output: vec![back.clone()], working_dir: dir.path().to_path_buf(),
            num_threads: 1, gzipped: false, gzip_level: 6, force: true,
        };
        crate::compression::cluster::decompress_cluster(&dc).unwrap();
        // set-lossless: every original (header,seq,qual) — including the flipped reads, restored to
        // their ORIGINAL orientation — must come back exactly.
        assert_eq!(multiset_records(&fq), multiset_records(&back));
    }

    #[test]
    fn cluster_roundtrip_is_set_lossless() {
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("in.fastq");
        write_fastq(&fq, 3000, 100);
        let qz = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![fq.clone()],
            output: qz.clone(),
            working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        crate::compression::cluster::compress_cluster(&cfg).unwrap();
        let back = dir.path().join("back.fastq");
        let dc = crate::cli::DecompressConfig {
            input: qz,
            output: vec![back.clone()],
            working_dir: dir.path().to_path_buf(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 6,
            force: true,
        };
        crate::compression::cluster::decompress_cluster(&dc).unwrap();
        assert_eq!(multiset_records(&fq), multiset_records(&back));
    }

    #[test]
    fn cluster_roundtrip_multichunk() {
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("in.fastq");
        write_fastq(&fq, 4000, 80);
        let qz = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![fq.clone()],
            output: qz.clone(),
            working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        // SAFETY: test is single-threaded for this var
        unsafe {
            std::env::set_var("QZ_CLUSTER_CHUNK_BYTES", "40000");
        }
        crate::compression::cluster::compress_cluster(&cfg).unwrap();
        unsafe {
            std::env::remove_var("QZ_CLUSTER_CHUNK_BYTES");
        }
        let back = dir.path().join("back.fastq");
        let dc = crate::cli::DecompressConfig {
            input: qz,
            output: vec![back.clone()],
            working_dir: dir.path().to_path_buf(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 6,
            force: true,
        };
        crate::compression::cluster::decompress_cluster(&dc).unwrap();
        assert_eq!(multiset_records(&fq), multiset_records(&back));
    }

    #[test]
    fn cluster_verify_passes_on_good_archive() {
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("in.fastq");
        write_fastq(&fq, 3000, 100);
        let qz = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![fq], output: qz.clone(), working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        crate::compression::cluster::compress_cluster(&cfg).unwrap();
        // deep verify
        let r = crate::compression::cluster::verify_cluster(&qz, false).unwrap();
        assert_eq!(r.num_reads, 3000);
        assert_eq!(r.mode, crate::compression::VerifyMode::Deep);
        // fast verify
        let rf = crate::compression::cluster::verify_cluster(&qz, true).unwrap();
        assert_eq!(rf.mode, crate::compression::VerifyMode::Fast);
        assert!(rf.blocks_verified > 0);
    }

    #[test]
    fn cluster_verify_detects_corruption() {
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("in.fastq");
        write_fastq(&fq, 1000, 80);
        let qz = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![fq], output: qz.clone(), working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        crate::compression::cluster::compress_cluster(&cfg).unwrap();
        // flip a byte in the middle of the archive payload region (after the header,
        // before the footer) to corrupt a block.
        let mut bytes = std::fs::read(&qz).unwrap();
        let mid = bytes.len() / 2;
        bytes[mid] ^= 0xFF;
        std::fs::write(&qz, &bytes).unwrap();
        // BOTH fast (per-block CRC) and deep verify must report failure (Err).
        assert!(crate::compression::cluster::verify_cluster(&qz, true).is_err(), "fast verify must catch corruption");
        assert!(crate::compression::cluster::verify_cluster(&qz, false).is_err(), "deep verify must catch corruption");
    }

    #[test]
    fn validate_paired_directory_requires_both_mates() {
        use crate::compression::chunk_directory::{ChunkDirEntry, ChunkDirectory, StreamRole, GLOBAL_SENTINEL};
        let meta = ClusterMeta { version: 1, kind: 1, level: 16, nreads: 4, checksum: 7 };
        // one chunk, both mates with all four roles + the global.
        let e = |chunk: u32, mate: u8, role: StreamRole, codec: u8, rc: u64| ChunkDirEntry {
            chunk_index: chunk, role, mate, codec, offset: 0, length: 100, record_count: rc,
        };
        let mut entries = vec![e(GLOBAL_SENTINEL, 0, StreamRole::ClusterMeta, 1, 0)];
        for mate in [1u8, 2u8] {
            entries.push(e(0, mate, StreamRole::Headers, 3, 4));
            entries.push(e(0, mate, StreamRole::ClusteredSeq, 7, 4));
            entries.push(e(0, mate, StreamRole::Qual, 6, 4));
            entries.push(e(0, mate, StreamRole::StrandBits, 1, 4));
        }
        let good = ChunkDirectory { num_reads: 4, num_chunks: 1, entries: entries.clone() };
        assert!(is_paired(&good));
        assert!(validate_cluster_directory(&good, &meta).is_ok());
        // drop R2's ClusteredSeq -> reject
        let bad_entries: Vec<_> = entries.into_iter()
            .filter(|x| !(x.mate == 2 && x.role == StreamRole::ClusteredSeq)).collect();
        let bad = ChunkDirectory { num_reads: 4, num_chunks: 1, entries: bad_entries };
        assert!(validate_cluster_directory(&bad, &meta).is_err());
    }

    #[test]
    fn validate_cluster_directory_kind_must_match_layout() {
        use crate::compression::chunk_directory::{ChunkDirEntry, ChunkDirectory, StreamRole, GLOBAL_SENTINEL};
        let e = |chunk: u32, mate: u8, role: StreamRole, codec: u8, rc: u64| ChunkDirEntry {
            chunk_index: chunk, role, mate, codec, offset: 0, length: 100, record_count: rc,
        };
        // A single-end layout (all mate 0, no mate-2 entries).
        let single_entries = vec![
            e(GLOBAL_SENTINEL, 0, StreamRole::ClusterMeta, 1, 0),
            e(0, 0, StreamRole::Headers, 3, 4),
            e(0, 0, StreamRole::ClusteredSeq, 7, 4),
            e(0, 0, StreamRole::Qual, 6, 4),
            e(0, 0, StreamRole::StrandBits, 1, 4),
        ];
        let single = ChunkDirectory { num_reads: 4, num_chunks: 1, entries: single_entries };
        // kind=0 agrees with single-end layout → ok.
        let meta0 = ClusterMeta { version: 1, kind: 0, level: 16, nreads: 4, checksum: 7 };
        assert!(validate_cluster_directory(&single, &meta0).is_ok());
        // kind=1 (paired) on a single-end layout → reject (cross-check).
        let meta1 = ClusterMeta { version: 1, kind: 1, level: 16, nreads: 4, checksum: 7 };
        assert!(validate_cluster_directory(&single, &meta1).is_err());

        // A paired layout (mates 1 and 2).
        let mut paired_entries = vec![e(GLOBAL_SENTINEL, 0, StreamRole::ClusterMeta, 1, 0)];
        for mate in [1u8, 2u8] {
            paired_entries.push(e(0, mate, StreamRole::Headers, 3, 4));
            paired_entries.push(e(0, mate, StreamRole::ClusteredSeq, 7, 4));
            paired_entries.push(e(0, mate, StreamRole::Qual, 6, 4));
            paired_entries.push(e(0, mate, StreamRole::StrandBits, 1, 4));
        }
        let paired = ChunkDirectory { num_reads: 4, num_chunks: 1, entries: paired_entries };
        // kind=0 (single) on a paired layout → reject (cross-check, the other direction).
        assert!(validate_cluster_directory(&paired, &meta0).is_err());
        // kind=1 agrees → ok.
        assert!(validate_cluster_directory(&paired, &meta1).is_ok());
    }

    /// Write `n` paired FASTQ records tagged with `mate_byte` ('1' or '2') so that
    /// post-roundtrip we can verify both per-mate multiset preservation AND that the
    /// i-th R1 and i-th R2 carry the same pair-index token.
    fn write_fastq_tagged(path: &std::path::Path, n: usize, len: usize, mate_byte: u8) {
        let mut f = std::fs::File::create(path).unwrap();
        let bases = b"ACGT";
        for i in 0..n {
            // header encodes both the pair index and the mate so we can check alignment.
            writeln!(f, "@pair{i}_m{}", mate_byte as char).unwrap();
            let seq: Vec<u8> = (0..len)
                .map(|j| bases[(i.wrapping_mul(13).wrapping_add(j.wrapping_mul(7)).wrapping_add(j.wrapping_mul(j)).wrapping_add(mate_byte as usize)) % 4])
                .collect();
            f.write_all(&seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            let q: Vec<u8> = (0..len).map(|j| b'!' + ((i + j + mate_byte as usize) % 40) as u8).collect();
            f.write_all(&q).unwrap();
            writeln!(f).unwrap();
        }
    }

    /// Assert that for every index i, the i-th record in `o1` and the i-th record
    /// in `o2` both contain the same `pair{N}` token in their header.
    fn pairs_aligned(o1: &std::path::Path, o2: &std::path::Path) -> bool {
        let s1 = std::fs::read(o1).unwrap();
        let s2 = std::fs::read(o2).unwrap();
        let headers1: Vec<&[u8]> = s1.split(|&b| b == b'\n')
            .enumerate()
            .filter(|(i, _)| i % 4 == 0)
            .map(|(_, h)| h)
            .filter(|h| !h.is_empty())
            .collect();
        let headers2: Vec<&[u8]> = s2.split(|&b| b == b'\n')
            .enumerate()
            .filter(|(i, _)| i % 4 == 0)
            .map(|(_, h)| h)
            .filter(|h| !h.is_empty())
            .collect();
        if headers1.len() != headers2.len() {
            return false;
        }
        for (h1, h2) in headers1.iter().zip(headers2.iter()) {
            // Extract the pair-index token: everything after '@' up to '_m'.
            let tok1 = h1.splitn(2, |&b| b == b'_').next().unwrap_or(b"");
            let tok2 = h2.splitn(2, |&b| b == b'_').next().unwrap_or(b"");
            if tok1 != tok2 {
                return false;
            }
        }
        true
    }

    #[test]
    fn cluster_paired_roundtrip_is_pair_lossless() {
        let dir = tempfile::tempdir().unwrap();
        let r1 = dir.path().join("r1.fastq");
        let r2 = dir.path().join("r2.fastq");
        // distinct, position-tagged contents so we can check pairing post-roundtrip.
        write_fastq_tagged(&r1, 2000, 100, b'1');
        write_fastq_tagged(&r2, 2000, 100, b'2');
        let qz = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![r1.clone(), r2.clone()],
            output: qz.clone(),
            working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        crate::compression::cluster::compress_cluster_paired(&cfg).unwrap();
        let o1 = dir.path().join("b1.fastq");
        let o2 = dir.path().join("b2.fastq");
        let dc = crate::cli::DecompressConfig {
            input: qz,
            output: vec![o1.clone(), o2.clone()],
            working_dir: dir.path().to_path_buf(),
            num_threads: 1,
            gzipped: false,
            gzip_level: 6,
            force: true,
        };
        crate::compression::cluster::decompress_cluster(&dc).unwrap();
        // each mate's record multiset is preserved
        assert_eq!(multiset_records(&r1), multiset_records(&o1));
        assert_eq!(multiset_records(&r2), multiset_records(&o2));
        // pairing intact: the i-th R1 read and i-th R2 read share their pair index tag
        assert!(pairs_aligned(&o1, &o2));
    }

    #[test]
    fn cluster_paired_verify_passes_and_detects_corruption() {
        let dir = tempfile::tempdir().unwrap();
        let r1 = dir.path().join("r1.fastq");
        let r2 = dir.path().join("r2.fastq");
        write_fastq_tagged(&r1, 1000, 80, b'1');
        write_fastq_tagged(&r2, 1000, 80, b'2');
        let qz = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![r1, r2], output: qz.clone(), working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        crate::compression::cluster::compress_cluster_paired(&cfg).unwrap();
        // fast verify passes on the good archive
        let rf = crate::compression::cluster::verify_cluster(&qz, true).unwrap();
        assert_eq!(rf.mode, crate::compression::VerifyMode::Fast);
        // deep verify passes
        let r = crate::compression::cluster::verify_cluster(&qz, false).unwrap();
        assert_eq!(r.num_reads, 1000); // pairs
        assert_eq!(r.mode, crate::compression::VerifyMode::Deep);
        // corrupt a payload byte -> both fast and deep verify fail
        let mut bytes = std::fs::read(&qz).unwrap();
        let mid = bytes.len() / 2;
        bytes[mid] ^= 0xFF;
        std::fs::write(&qz, &bytes).unwrap();
        assert!(crate::compression::cluster::verify_cluster(&qz, true).is_err());
        assert!(crate::compression::cluster::verify_cluster(&qz, false).is_err());
    }

    #[test]
    fn validate_cluster_directory_cap_rejects_oversized_headers() {
        // A well-formed directory passes; one with an oversized Headers entry fails.
        let meta = ClusterMeta {
            version: 1,
            kind: 0,
            level: 16,
            nreads: 10,
            checksum: 0,
        };

        // Well-formed: one chunk with all required roles, within caps.
        let good_dir = ChunkDirectory {
            num_reads: 10,
            num_chunks: 1,
            entries: vec![
                ChunkDirEntry {
                    chunk_index: GLOBAL_SENTINEL,
                    role: StreamRole::ClusterMeta,
                    mate: 0,
                    codec: 1,
                    offset: 0,
                    length: 100,
                    record_count: 0,
                },
                ChunkDirEntry {
                    chunk_index: 0,
                    role: StreamRole::ClusteredSeq,
                    mate: 0,
                    codec: 7,
                    offset: 100,
                    length: 200,
                    record_count: 10,
                },
                ChunkDirEntry {
                    chunk_index: 0,
                    role: StreamRole::StrandBits,
                    mate: 0,
                    codec: 1,
                    offset: 300,
                    length: 50,
                    record_count: 10,
                },
                ChunkDirEntry {
                    chunk_index: 0,
                    role: StreamRole::Qual,
                    mate: 0,
                    codec: 6,
                    offset: 350,
                    length: 150,
                    record_count: 10,
                },
                ChunkDirEntry {
                    chunk_index: 0,
                    role: StreamRole::Headers,
                    mate: 0,
                    codec: 3,
                    offset: 500,
                    length: 200,
                    record_count: 10,
                },
            ],
        };
        assert!(
            validate_cluster_directory(&good_dir, &meta).is_ok(),
            "well-formed directory should pass"
        );

        // Bad: Headers entry length = 100 MiB > BSC_MAX_BLOCK_LEN (~64 MiB).
        let oversized_len = 100 << 20; // 100 MiB
        let bad_dir = ChunkDirectory {
            num_reads: 10,
            num_chunks: 1,
            entries: {
                let mut entries = good_dir.entries.clone();
                // Replace the Headers entry with an oversized one.
                for e in entries.iter_mut() {
                    if e.role == StreamRole::Headers && e.chunk_index == 0 {
                        e.length = oversized_len as u64;
                        break;
                    }
                }
                entries
            },
        };
        assert!(
            validate_cluster_directory(&bad_dir, &meta).is_err(),
            "oversized Headers entry should be rejected"
        );
    }

    #[test]
    fn cluster_meta_from_bytes_guards() {
        // 28 well-formed bytes, version 1 → ok.
        let mut good = vec![0u8; 28];
        good[0] = 1; // version
        good[1] = 0; // kind = single-end
        good[2] = 16; // level
        good[4..12].copy_from_slice(&7u64.to_le_bytes()); // nreads
        assert!(ClusterMeta::from_bytes(&good).is_ok());

        // Unsupported version → err (the guard the roundtrip tests never reach).
        let mut bad_ver = good.clone();
        bad_ver[0] = 2;
        assert!(ClusterMeta::from_bytes(&bad_ver).is_err(), "version 2 must be rejected");

        // Wrong length (both short and long) → err.
        assert!(ClusterMeta::from_bytes(&good[..27]).is_err(), "27 bytes must be rejected");
        let mut too_long = good.clone();
        too_long.push(0);
        assert!(ClusterMeta::from_bytes(&too_long).is_err(), "29 bytes must be rejected");
    }

    #[test]
    fn parse_fqz_varint_stream_clamps_hostile_record_count() {
        // A forged footer can claim a huge record_count; the per-record loop must
        // error cleanly on the truncated data WITHOUT pre-allocating from the
        // attacker count (the `.min(data.len() + 1)` clamp). A few length-prefixed
        // records, then a count far larger than the data can satisfy.
        let mut data = Vec::new();
        for _ in 0..3 {
            crate::compression::dna_utils::write_varint(&mut data, 4);
            data.extend_from_slice(b"IIII");
        }
        // Honest count parses fine.
        assert_eq!(parse_fqz_varint_stream(&data, 3).unwrap().len(), 3);
        // Hostile huge count: must Err (truncation), not attempt a u64::MAX alloc.
        assert!(parse_fqz_varint_stream(&data, usize::MAX).is_err());
        assert!(parse_fqz_varint_stream(&data, 1_000_000_000).is_err());
    }
}
