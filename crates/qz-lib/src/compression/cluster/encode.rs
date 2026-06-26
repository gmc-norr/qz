//! Encode pipeline for the order-drop cluster mode (v5 `archive_type` = 3,
//! `EncodingType::ClusterReorder` = 11).
//!
//! Entry points:
//! - [`compress_cluster`]   — full encode pipeline
//! - [`reject_unsupported`] — pre-flight config gate (mirrors reference mode)

use anyhow::Result;
use std::io::Write as _;

use crate::cli::CompressConfig;
use crate::compression::bsc;
use crate::compression::chunk_directory::{
    ChunkDirEntry, ChunkDirectory, GLOBAL_SENTINEL,
};
use crate::compression::codec_ids::{CODEC_BSC, CODEC_COLUMNAR, CODEC_FQZCOMP, CODEC_ZSTD_LONG};
use crate::compression::compress_impl::{
    columnar_blocks_capped, fqz_blocks_capped_quals,
};
use crate::compression::cluster::sort::{cluster_sort_paired, cluster_sort_single, ClusteredRecord};
use crate::compression::cluster::zstd_codec::compress_seq_zstd_long;
use crate::compression::chunk_directory::StreamRole;
use crate::compression::decompress_impl::BSC_MAX_BLOCK_LEN;
use crate::compression::{write_v5_archive_header, ModeCaps};
use crate::compression::columnar::QualityBinning;
use crate::cli::{HeaderCompressor, QualityCompressor};

// ── Gate ──────────────────────────────────────────────────────────────────────

/// Pre-flight rejection gate for cluster mode.
///
/// Mirrors `reference::reject_unsupported` — rejects unsupported options and
/// enforces single-end-only for this release.
pub(crate) fn reject_unsupported(a: &CompressConfig) -> Result<()> {
    if a.cluster.is_none() {
        return Ok(());
    }
    // Cluster accepts 1 (single-end) or 2 (paired-end) inputs; 3+ are rejected.
    if a.input.is_empty() || a.input.len() > 2 {
        anyhow::bail!("--cluster accepts one (single-end) or two (paired-end) input files (got {})", a.input.len());
    }
    crate::compression::reject_unsupported_common(
        a,
        &ModeCaps {
            mode: "cluster",
            allows_explicit_fqzcomp: true,
            allows_stdout: false, // cluster needs a named output for the atomic write
        },
    )?;
    Ok(())
}

// ── Emit helper ───────────────────────────────────────────────────────────────

/// Write one block frame (12-byte header + payload) to `out` and push a
/// `ChunkDirEntry`.  Updates `*cur` by exactly `12 + payload.len()`.
fn emit_block(
    out: &mut impl std::io::Write,
    cur: &mut u64,
    entries: &mut Vec<ChunkDirEntry>,
    chunk_index: u32,
    mate: u8,
    role: StreamRole,
    codec: u8,
    record_count: u32,
    payload: &[u8],
) -> Result<()> {
    let mut framing = Vec::new();
    bsc::write_block_frame_header(&mut framing, payload.len() as u32, record_count, payload);
    let frame_offset = *cur;
    let frame_len = (framing.len() + payload.len()) as u64; // 12 + payload.len()
    out.write_all(&framing)?;
    out.write_all(payload)?;
    entries.push(ChunkDirEntry {
        chunk_index,
        role,
        mate,
        codec,
        offset: frame_offset,
        length: frame_len,
        record_count: record_count as u64,
    });
    *cur += frame_len;
    Ok(())
}

// ── ChunkArenas ───────────────────────────────────────────────────────────────

/// Per-stream flat storage for one chunk's records, replacing `Vec<ClusteredRecord>`.
/// Appending avoids the per-record `Vec` allocation + clone the old buffer paid for every
/// record; the streams are stored contiguously so the flush can compress each directly.
#[derive(Default)]
struct ChunkArenas {
    seq: Vec<u8>,
    qual: Vec<u8>,
    qual_lens: Vec<u32>,
    headers: Vec<u8>,
    hdr_lens: Vec<u32>,
    /// Bitvec, 1 bit/record (bit i set iff record i was reverse-complemented). Length is
    /// kept at exactly `ceil(nrec/8)` bytes — identical packing to the old StrandBits.
    flips: Vec<u8>,
    nrec: u32,
}

impl ChunkArenas {
    fn push(&mut self, r: &ClusteredRecord) {
        self.seq.extend_from_slice(&r.canon_seq);
        self.qual.extend_from_slice(&r.qual);
        self.qual_lens.push(r.qual.len() as u32);
        self.headers.extend_from_slice(&r.header);
        self.hdr_lens.push(r.header.len() as u32);
        let i = self.nrec as usize;
        if (i >> 3) >= self.flips.len() {
            self.flips.push(0); // grows by exactly 1 byte every 8 records -> ceil(nrec/8)
        }
        if r.flip {
            self.flips[i >> 3] |= 1 << (i & 7);
        }
        self.nrec += 1;
    }

    fn is_empty(&self) -> bool {
        self.nrec == 0
    }

    fn clear(&mut self) {
        self.seq.clear();
        self.qual.clear();
        self.qual_lens.clear();
        self.headers.clear();
        self.hdr_lens.clear();
        self.flips.clear();
        self.nrec = 0;
    }
}

/// Iterate the per-record slices packed into `arena` according to `lens`.
fn arena_slices<'a>(arena: &'a [u8], lens: &'a [u32]) -> impl Iterator<Item = &'a [u8]> + 'a {
    let mut off = 0usize;
    lens.iter().map(move |&l| {
        let s = &arena[off..off + l as usize];
        off += l as usize;
        s
    })
}

// ── flush_chunk ───────────────────────────────────────────────────────────────

/// Compress one chunk's records and write all four stream blocks
/// (Headers → ClusteredSeq → Qual → StrandBits) to `out`.
///
/// Free function (not a closure) to avoid a double-borrow of `out` + `entries`.
///
/// This is functionally `compress_mate_blocks(mate = 0, ..)` followed by an emit; it
/// is kept SEPARATE (not delegated) so single-end output stays decoupled from any
/// future change to the paired block producer. The two bodies MUST stay in lockstep —
/// if you change the block order, codecs, `window_log`, or StrandBits packing here,
/// mirror it in `compress_mate_blocks` (and vice versa), or single-end and paired
/// sequences will diverge — and the decode path (`decode_chunk_records` /
/// `validate_chunk_mate_roles` in `decode.rs`) must accept the resulting role set/codecs.
fn flush_chunk(
    out: &mut impl std::io::Write,
    cur: &mut u64,
    entries: &mut Vec<ChunkDirEntry>,
    chunk_index: u32,
    chunk: &ChunkArenas,
    level: i32,
    zstd_workers: u32,
) -> Result<()> {
    let nrec = chunk.nrec;

    // ── Headers (N columnar blocks, codec 3) ──────────────────────────────────
    let hrefs: Vec<&[u8]> = arena_slices(&chunk.headers, &chunk.hdr_lens).collect();
    for (rc, blob) in columnar_blocks_capped(&hrefs, BSC_MAX_BLOCK_LEN)? {
        emit_block(out, cur, entries, chunk_index, 0, StreamRole::Headers, CODEC_COLUMNAR, rc, &blob)?;
    }

    // ── ClusteredSeq (ONE zstd-long block, codec 7) ───────────────────────────
    let window_log = (64 - (chunk.seq.len().max(1) as u64).leading_zeros()).clamp(10, 31);
    let payload = compress_seq_zstd_long(&chunk.seq, level, window_log, zstd_workers)?;
    emit_block(out, cur, entries, chunk_index, 0, StreamRole::ClusteredSeq, CODEC_ZSTD_LONG, nrec, &payload)?;

    // ── Qual (N fqzcomp blocks, codec 6) ─────────────────────────────────────
    let qrefs: Vec<&[u8]> = arena_slices(&chunk.qual, &chunk.qual_lens).collect();
    for (rc, blob) in fqz_blocks_capped_quals(&qrefs, BSC_MAX_BLOCK_LEN)? {
        emit_block(out, cur, entries, chunk_index, 0, StreamRole::Qual, CODEC_FQZCOMP, rc, &blob)?;
    }

    // ── StrandBits (ONE BSC block, codec 1) ──────────────────────────────────
    // `chunk.flips` is already the exact ceil(nrec/8)-byte bitvec.
    let strand_payload = bsc::compress_parallel(&chunk.flips)?;
    emit_block(out, cur, entries, chunk_index, 0, StreamRole::StrandBits, CODEC_BSC, nrec, &strand_payload)?;

    Ok(())
}

// ── compress_mate_blocks / flush_chunk_paired ────────────────────────────────

/// One compressed stream block, ready to write. Produced by `compress_mate_blocks`
/// (possibly off-thread) and emitted later in a fixed order, so the codec work can
/// run concurrently while the archive bytes stay deterministic.
struct PreparedBlock {
    mate: u8,
    role: StreamRole,
    codec: u8,
    record_count: u32,
    payload: Vec<u8>,
}

/// Compress one mate's four streams and return their blocks in the fixed emit order
/// (Headers → ClusteredSeq → Qual → StrandBits), WITHOUT writing to `out`.
///
/// Splitting "compress" from "emit" lets `flush_chunk_paired` run both mates' codec
/// work concurrently and then write the blocks sequentially — the payloads and their
/// order are identical to the old serial `flush_mate`, so the archive is byte-for-byte
/// unchanged. The four streams stay sequential WITHIN a mate; the parallelism is across
/// mates (R1 ∥ R2), where the win is: the two zstd-long seq compressions overlap (each
/// mate gets its own `resolve_zstd_workers` pool sized to `-t / 2` lanes, so the two
/// concurrent mates together ≈ saturate `-t` and overlap their MT scaling), as do the
/// two fqzcomp quality compressions.
///
/// RAM tradeoff: returning a mate's compressed blocks (rather than streaming each to
/// `out`) means the caller holds the compressed block sets resident before emitting, and
/// the two zstd-long contexts run concurrently. To bound this, `chunk` is taken **by value**
/// and each stream's arena is freed as soon as its block(s) are produced
/// (headers → seq → qual), so a mate's raw residency drains during its own flush — which is
/// what keeps the pipelined two-chunks-in-flight peak from doubling.
///
/// Byte-equivalent to single-end `flush_chunk` (mate hard-coded to 0). The two MUST
/// stay in lockstep — block order, codecs, `window_log`, StrandBits packing — and the
/// decode path (`decode_chunk_records` / `validate_chunk_mate_roles` in `decode.rs`)
/// must accept the resulting role set and codecs.
fn compress_mate_blocks(
    mate: u8,
    mut chunk: ChunkArenas,
    level: i32,
    zstd_workers: u32,
) -> Result<Vec<PreparedBlock>> {
    let nrec = chunk.nrec;
    // Common case = exactly 4 blocks (1 per stream); Headers/Qual split only when a
    // stream exceeds BSC_MAX_BLOCK_LEN, so 4 is the floor.
    let mut blocks: Vec<PreparedBlock> = Vec::with_capacity(4);

    // ── Headers (N columnar blocks, codec 3) ──────────────────────────────────
    {
        let hrefs: Vec<&[u8]> = arena_slices(&chunk.headers, &chunk.hdr_lens).collect();
        for (rc, blob) in columnar_blocks_capped(&hrefs, BSC_MAX_BLOCK_LEN)? {
            blocks.push(PreparedBlock { mate, role: StreamRole::Headers, codec: CODEC_COLUMNAR, record_count: rc, payload: blob });
        }
    }
    chunk.headers = Vec::new();
    chunk.hdr_lens = Vec::new();

    // ── ClusteredSeq (ONE zstd-long block, codec 7) ───────────────────────────
    let window_log = (64 - (chunk.seq.len().max(1) as u64).leading_zeros()).clamp(10, 31);
    let seq_payload = compress_seq_zstd_long(&chunk.seq, level, window_log, zstd_workers)?;
    blocks.push(PreparedBlock { mate, role: StreamRole::ClusteredSeq, codec: CODEC_ZSTD_LONG, record_count: nrec, payload: seq_payload });
    chunk.seq = Vec::new();

    // ── Qual (N fqzcomp blocks, codec 6) ─────────────────────────────────────
    {
        let qrefs: Vec<&[u8]> = arena_slices(&chunk.qual, &chunk.qual_lens).collect();
        for (rc, blob) in fqz_blocks_capped_quals(&qrefs, BSC_MAX_BLOCK_LEN)? {
            blocks.push(PreparedBlock { mate, role: StreamRole::Qual, codec: CODEC_FQZCOMP, record_count: rc, payload: blob });
        }
    }
    chunk.qual = Vec::new();
    chunk.qual_lens = Vec::new();

    // ── StrandBits (ONE BSC block, codec 1) ──────────────────────────────────
    // `chunk.flips` is already the exact ceil(nrec/8)-byte bitvec.
    let strand_payload = bsc::compress_parallel(&chunk.flips)?;
    blocks.push(PreparedBlock { mate, role: StreamRole::StrandBits, codec: CODEC_BSC, record_count: nrec, payload: strand_payload });

    Ok(blocks)
}

/// Compress one paired chunk and emit R1 roles under mate 1 then R2 roles under mate 2.
///
/// Both mates' codec work runs concurrently via `rayon::join`; the blocks are then
/// written in the fixed R1-then-R2 order, so the output is byte-identical to the old
/// serial two-`flush_mate` path. Emission happens only after BOTH mates succeed, so a
/// failure never writes a torn chunk. On error, R1's error takes precedence when both
/// fail (`r1_res?` before `r2_res?` — matching the old `flush_mate(1)?;flush_mate(2)?`
/// ordering); `rayon::join` can't cancel an in-flight branch, so R2's codec work runs
/// to completion before its result is discarded — bounded wasted work on the rare error
/// path only.
fn flush_chunk_paired(
    out: &mut impl std::io::Write,
    cur: &mut u64,
    entries: &mut Vec<ChunkDirEntry>,
    chunk_index: u32,
    r1_chunk: ChunkArenas,
    r2_chunk: ChunkArenas,
    level: i32,
    zstd_workers: u32,
) -> Result<()> {
    let (r1_res, r2_res) = rayon::join(
        || compress_mate_blocks(1, r1_chunk, level, zstd_workers),
        || compress_mate_blocks(2, r2_chunk, level, zstd_workers),
    );
    let r1_blocks = r1_res?;
    let r2_blocks = r2_res?;

    for b in r1_blocks.iter().chain(r2_blocks.iter()) {
        emit_block(out, cur, entries, chunk_index, b.mate, b.role, b.codec, b.record_count, &b.payload)?;
    }
    Ok(())
}

// ── Shared encode helpers ─────────────────────────────────────────────────────

/// Upper bound on per-mate zstd workers. zstd-MT seq-compress throughput plateaus around
/// here (measured ~48 workers on a 2×18-core Xeon; the seq stage scales ~1.9× to this knee,
/// ~1.36× whole-compress) and each worker preallocates a job buffer, so the cap bounds
/// per-chunk compress RSS at extreme `-t` while capturing essentially all the MT speedup.
/// `QZ_CLUSTER_ZSTD_WORKERS` overrides it for explicit control.
const ZSTD_WORKERS_MAX: usize = 48;

/// Resolve the **per-mate** zstd worker count for the per-chunk seq block.
///
/// Each chunk's seq is ~one window, but with LDM zstd still parallelizes it into jobs that
/// the long-range table spans, so threading buys compress speed without losing the
/// long-range ratio. Worker count **FOLLOWS `-t`** (the user's parallelism/RAM dial,
/// `0` = auto-detect), split across `concurrent_lanes` simultaneous zstd calls, then clamped
/// to `[1, ZSTD_WORKERS_MAX]` (previously a hard 16, which left the seq stage ~1.9× slower
/// on big boxes). `QZ_CLUSTER_ZSTD_WORKERS` overrides the resolved value with an explicit
/// per-mate count (escape hatch; bypasses the `-t` derivation AND the lane split, so a
/// paired run then uses 2× that many workers in total).
///
/// Determinism: zstd-MT job boundaries depend on `(input, level, window_log)`, not thread
/// scheduling, so output is deterministic for a fixed worker count on a given machine.
/// Cluster's roundtrip oracle is a multiset checksum / `qz verify` (it drops read order), so
/// correctness does not depend on byte-identity across worker counts.
///
/// NOTE on the paired path: `flush_chunk_paired` runs the two mates concurrently
/// (`rayon::join`), so it passes `concurrent_lanes = 2` to split the thread budget across the
/// two simultaneous zstd-long calls, keeping the *total* default worker count ≈ `-t`. The two
/// zstd contexts ~double the per-call memory live at peak (the source of paired's higher RSS),
/// bounded by the per-chunk seq target and constant in read count.
fn resolve_zstd_workers(cfg: &CompressConfig, concurrent_lanes: u32) -> u32 {
    // Explicit override: exact per-mate worker count, bypassing the -t derivation and lane split.
    if let Some(n) = std::env::var("QZ_CLUSTER_ZSTD_WORKERS")
        .ok()
        .and_then(|v| v.trim().parse::<u32>().ok())
    {
        return n;
    }
    let resolved = if cfg.threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        cfg.threads
    };
    zstd_workers_from_threads(resolved, concurrent_lanes)
}

/// Pure per-mate worker derivation (no env read, no thread probe — unit-testable): split
/// `threads` across `concurrent_lanes` simultaneous zstd calls, clamped to `[1, ZSTD_WORKERS_MAX]`.
/// For `threads ≤ 16` the result is identical to the historical hard-16 cap.
fn zstd_workers_from_threads(threads: usize, concurrent_lanes: u32) -> u32 {
    (threads / (concurrent_lanes.max(1) as usize)).clamp(1, ZSTD_WORKERS_MAX) as u32
}

/// Assemble the raw ClusterMeta payload (pre-BSC):
/// version(1) | archive_kind(1) | level(1) | reserved(1) | nreads(8, u64 LE) |
/// checksum(16, u128 LE) = 28 bytes. `kind` is 0 = single-end, 1 = paired.
/// Mirrors `ClusterMeta::from_bytes` exactly.
fn cluster_meta_bytes(kind: u8, level: i32, nreads: u64, checksum: u128) -> Vec<u8> {
    let mut meta = Vec::with_capacity(28);
    meta.push(1u8); // version
    meta.push(kind);
    meta.push(level as u8);
    meta.push(0u8); // reserved
    meta.extend_from_slice(&nreads.to_le_bytes()); // 8 bytes, u64 LE
    meta.extend_from_slice(&checksum.to_le_bytes()); // 16 bytes, u128 LE
    meta
}

// ── Public entry point ────────────────────────────────────────────────────────

/// Encode pipeline for cluster mode.
///
/// Writes directly to `cfg.output` (caller already placed a `.tmp` path there
/// and owns the atomic rename on success).  Does NOT create its own temp file.
pub(crate) fn compress_cluster(cfg: &CompressConfig) -> Result<()> {
    let level = cfg
        .cluster
        .as_ref()
        .map(|c| c.level.zstd_level())
        .unwrap_or(16);

    // Single-end flushes one mate at a time → one zstd lane.
    let zstd_workers = resolve_zstd_workers(cfg, 1);

    // Open output (already a .tmp path from the atomic-write wrapper in compress()).
    let file = std::fs::File::create(&cfg.output)?;
    let mut out = std::io::BufWriter::new(file);

    // Write front header: version=5, archive_type=3, encoding_type=11.
    let mut cur: u64 = write_v5_archive_header(
        &mut out,
        11, // EncodingType::ClusterReorder
        QualityBinning::None,
        QualityCompressor::Fqzcomp,
        HeaderCompressor::Columnar,
        false,
        0, // const_seq_len = 0 (variable)
        0, // const_qual_len = 0 (variable)
        3, // archive_type = 3 (cluster)
    )? as u64;

    // ── Streaming sort → chunk accumulation ──────────────────────────────────
    let mut chunk = ChunkArenas::default();
    let mut chunk_seq_bytes: usize = 0;
    let mut num_chunks: u32 = 0; // also the next chunk_index; incremented once per flush
    let mut entries: Vec<ChunkDirEntry> = Vec::new();

    // 2 GB default chunk; QZ_CLUSTER_CHUNK_BYTES overrides for testing.
    let target: usize = std::env::var("QZ_CLUSTER_CHUNK_BYTES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2_000_000_000usize);

    let (checksum, nreads) = cluster_sort_single(&cfg.input[0], &cfg.working_dir, |r| {
        chunk_seq_bytes += r.canon_seq.len();
        chunk.push(r);
        if chunk_seq_bytes >= target {
            flush_chunk(&mut out, &mut cur, &mut entries, num_chunks, &chunk, level, zstd_workers)?;
            num_chunks += 1;
            chunk.clear();
            chunk_seq_bytes = 0;
        }
        Ok(())
    })?;

    // Flush any remaining records (always runs if input is non-empty and didn't
    // hit the target, including the very common single-chunk case).
    if !chunk.is_empty() {
        flush_chunk(&mut out, &mut cur, &mut entries, num_chunks, &chunk, level, zstd_workers)?;
        num_chunks += 1;
    }

    // ── ClusterMeta global (GLOBAL_SENTINEL, codec 1) ─────────────────────────
    let meta = cluster_meta_bytes(0, level, nreads, checksum); // 0 = single-end
    let meta_payload = bsc::compress_parallel(&meta)?;
    emit_block(
        &mut out,
        &mut cur,
        &mut entries,
        GLOBAL_SENTINEL,
        0,
        StreamRole::ClusterMeta,
        CODEC_BSC,
        0,
        &meta_payload,
    )?;

    // ── Footer + locator ──────────────────────────────────────────────────────
    let dir = ChunkDirectory { num_reads: nreads, num_chunks, entries };
    // Encode self-check: validate the directory we are about to write.
    let self_check_meta = crate::compression::cluster::decode::ClusterMeta {
        version: 1,
        kind: 0,
        level: level as u8,
        nreads,
        checksum,
    };
    crate::compression::cluster::decode::validate_cluster_directory(&dir, &self_check_meta)?;
    crate::compression::chunk_directory::write_footer_and_locator(&mut out, &dir)?;
    out.flush()?;
    Ok(())
}

// ── Paired encode entry point ─────────────────────────────────────────────────

/// Encode pipeline for paired-end cluster mode.
///
/// Mirrors `compress_cluster` structurally.  R1 blocks are tagged mate 1,
/// R2 blocks mate 2.  `ClusterMeta` carries `archive_kind = 1` (paired) and
/// `total_reads = n_pairs`.  Writes directly to `cfg.output` (caller already
/// placed a `.tmp` path there and owns the atomic rename on success); does NOT
/// create its own temp file.
pub(crate) fn compress_cluster_paired(cfg: &CompressConfig) -> Result<()> {
    let level = cfg
        .cluster
        .as_ref()
        .map(|c| c.level.zstd_level())
        .unwrap_or(16);

    // Paired flushes both mates concurrently (rayon::join) → split the budget across 2 lanes.
    let zstd_workers = resolve_zstd_workers(cfg, 2);

    // Open output (already a .tmp path from the atomic-write wrapper in compress()).
    let file = std::fs::File::create(&cfg.output)?;
    let mut out = std::io::BufWriter::new(file);

    // Write front header: version=5, archive_type=3, encoding_type=11.
    let mut cur: u64 = write_v5_archive_header(
        &mut out,
        11, // EncodingType::ClusterReorder
        QualityBinning::None,
        QualityCompressor::Fqzcomp,
        HeaderCompressor::Columnar,
        false,
        0, // const_seq_len = 0 (variable)
        0, // const_qual_len = 0 (variable)
        3, // archive_type = 3 (cluster)
    )? as u64;

    // 2 GB default chunk; QZ_CLUSTER_CHUNK_BYTES overrides for testing.
    let target: usize = std::env::var("QZ_CLUSTER_CHUNK_BYTES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2_000_000_000usize);

    // Optional phase timing (QZ_CLUSTER_TIMING=1): isolates codec-flush wall (now in the
    // background flush thread) from the sort. Zero cost when the var is unset.
    let timing = std::env::var("QZ_CLUSTER_TIMING").is_ok();
    let wall_start = std::time::Instant::now();

    // ── Pipeline: producer (sort) ∥ background flush+write ─────────────────────
    // The sort callback fills per-mate arenas and hands each completed chunk to a
    // background flush thread over a depth-1 channel — backpressure caps in-flight at
    // 2 chunks (one flushing while the next fills). The flush thread is the SOLE writer of
    // `out` and processes chunks in FIFO (= chunk_index) order, so the archive is
    // byte-identical to the old synchronous flush at the same chunk size. First-error-wins:
    // a flush error disconnects the channel (stopping the producer) and is returned in
    // preference to the producer's resulting send-failure symptom; the caller's atomic-`.tmp`
    // drop guard removes a partial archive on any error/panic.
    type ChunkMsg = (u32, ChunkArenas, ChunkArenas);
    let header_len = cur;
    let out_ref = &mut out;
    let (checksum, n_pairs, cur_after, dir_entries, flush_total, flush_count) =
        std::thread::scope(
            |s| -> Result<(u128, u64, u64, Vec<ChunkDirEntry>, std::time::Duration, u32)> {
                let (tx, rx) = std::sync::mpsc::sync_channel::<ChunkMsg>(1);
                let writer = s.spawn(
                    move || -> Result<(u64, Vec<ChunkDirEntry>, std::time::Duration, u32)> {
                        let mut cur = header_len;
                        let mut entries: Vec<ChunkDirEntry> = Vec::new();
                        let mut flush_total = std::time::Duration::ZERO;
                        let mut flush_count: u32 = 0;
                        while let Ok((idx, r1a, r2a)) = rx.recv() {
                            let t = std::time::Instant::now();
                            flush_chunk_paired(out_ref, &mut cur, &mut entries, idx, r1a, r2a, level, zstd_workers)?;
                            flush_total += t.elapsed();
                            flush_count += 1;
                        }
                        Ok((cur, entries, flush_total, flush_count))
                    },
                );

                // Producer runs on this thread; the closure borrows the chunk buffers + tx.
                let mut r1_chunk = ChunkArenas::default();
                let mut r2_chunk = ChunkArenas::default();
                let mut chunk_seq_bytes: usize = 0;
                let mut num_chunks: u32 = 0;
                let tx_cb = tx.clone();
                let producer_res: Result<(u128, u64)> = (|| {
                    let (checksum, n_pairs) = cluster_sort_paired(
                        &cfg.input[0],
                        &cfg.input[1],
                        &cfg.working_dir,
                        |r1, r2| {
                            chunk_seq_bytes += r1.canon_seq.len() + r2.canon_seq.len();
                            r1_chunk.push(r1);
                            r2_chunk.push(r2);
                            if chunk_seq_bytes >= target {
                                let r1a = std::mem::take(&mut r1_chunk);
                                let r2a = std::mem::take(&mut r2_chunk);
                                tx_cb
                                    .send((num_chunks, r1a, r2a))
                                    .map_err(|_| anyhow::anyhow!("cluster flush stage terminated"))?;
                                num_chunks += 1;
                                chunk_seq_bytes = 0;
                            }
                            Ok(())
                        },
                    )?;
                    // Final partial chunk (lockstep: non-empty R1 ⇒ equal-length R2).
                    debug_assert_eq!(r1_chunk.nrec, r2_chunk.nrec);
                    if !r1_chunk.is_empty() {
                        let r1a = std::mem::take(&mut r1_chunk);
                        let r2a = std::mem::take(&mut r2_chunk);
                        tx_cb
                            .send((num_chunks, r1a, r2a))
                            .map_err(|_| anyhow::anyhow!("cluster flush stage terminated"))?;
                    }
                    Ok((checksum, n_pairs))
                })();
                // Close the channel so the writer's `rx.recv()` returns Err and the loop
                // terminates. The writer captures only `rx` (never a sender), so `tx` and
                // `tx_cb` are the ONLY live senders — both must drop here. Then prefer the
                // writer's (root-cause) error over the producer's send-failure symptom.
                drop(tx_cb);
                drop(tx);
                let (cur_after, entries, flush_total, flush_count) = match writer.join() {
                    Ok(r) => r?,
                    Err(_) => anyhow::bail!("cluster flush thread panicked"),
                };
                let (checksum, n_pairs) = producer_res?;
                Ok((checksum, n_pairs, cur_after, entries, flush_total, flush_count))
            },
        )?;
    cur = cur_after;
    let mut entries = dir_entries;
    let num_chunks = flush_count;

    // ── ClusterMeta global (GLOBAL_SENTINEL, codec 1, mate 0) ────────────────
    let meta = cluster_meta_bytes(1, level, n_pairs, checksum); // 1 = paired
    let meta_payload = bsc::compress_parallel(&meta)?;
    emit_block(
        &mut out,
        &mut cur,
        &mut entries,
        GLOBAL_SENTINEL,
        0,
        StreamRole::ClusterMeta,
        CODEC_BSC,
        0,
        &meta_payload,
    )?;

    // ── Footer + locator ──────────────────────────────────────────────────────
    let dir = ChunkDirectory { num_reads: n_pairs, num_chunks, entries };
    // Encode self-check: validate the directory we are about to write.
    let self_check_meta = crate::compression::cluster::decode::ClusterMeta {
        version: 1,
        kind: 1,
        level: level as u8,
        nreads: n_pairs,
        checksum,
    };
    crate::compression::cluster::decode::validate_cluster_directory(&dir, &self_check_meta)?;
    crate::compression::chunk_directory::write_footer_and_locator(&mut out, &dir)?;
    out.flush()?;

    if timing {
        eprintln!(
            "[CLUSTER_TIMING paired] codec_flush(zstd+fqz+bsc)={:.1}s over {} chunk-flush(es); total_compress_wall={:.1}s; n_pairs={n_pairs} chunks={num_chunks}",
            flush_total.as_secs_f64(),
            flush_count,
            wall_start.elapsed().as_secs_f64(),
        );
    }
    Ok(())
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
            let seq: Vec<u8> = (0..len).map(|j| bases[(i * 7 + j * 3) % 4]).collect();
            writeln!(f, "@read{i}").unwrap();
            f.write_all(&seq).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(len)).unwrap();
        }
    }

    #[test]
    fn zstd_worker_derivation_caps_and_lane_splits() {
        // -t ≤ 16: identical to the historical hard-16 cap (small-box behaviour preserved).
        assert_eq!(zstd_workers_from_threads(1, 1), 1, "-t 1 single");
        assert_eq!(zstd_workers_from_threads(1, 2), 1, "-t 1 paired: 1/2=0 → clamp to 1");
        assert_eq!(zstd_workers_from_threads(8, 1), 8, "-t 8 single unchanged");
        assert_eq!(zstd_workers_from_threads(8, 2), 4, "-t 8 paired: 8/2");
        assert_eq!(zstd_workers_from_threads(16, 1), 16, "-t 16 single unchanged");
        // -t > old cap: now follows -t up to the plateau, split across paired lanes.
        assert_eq!(zstd_workers_from_threads(72, 1), 48, "-t 72 single: capped at plateau");
        assert_eq!(zstd_workers_from_threads(72, 2), 36, "-t 72 paired: 72/2 < cap");
        assert_eq!(zstd_workers_from_threads(200, 1), 48, "extreme -t single: capped");
        assert_eq!(zstd_workers_from_threads(200, 2), 48, "extreme -t paired: 100 → capped");
        assert_eq!(
            zstd_workers_from_threads(40, 0),
            40,
            "lanes 0 treated as 1 (no div-by-zero)"
        );
    }

    #[test]
    fn compress_cluster_writes_parseable_v5() {
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("in.fastq");
        write_fastq(&fq, 2000, 100);
        let out = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![fq],
            output: out.clone(),
            working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        compress_cluster(&cfg).unwrap();
        let bytes = std::fs::read(&out).unwrap();
        assert_eq!(bytes[3], 3, "archive_type byte must be 3 (cluster)");
        let hdr = crate::compression::FixedHeader::parse_v5(&bytes).unwrap();
        assert_eq!(hdr.archive_type, 3);
    }

    #[test]
    fn compress_cluster_paired_writes_mate_tagged_v5() {
        let dir = tempfile::tempdir().unwrap();
        let r1 = dir.path().join("r1.fastq");
        let r2 = dir.path().join("r2.fastq");
        write_fastq(&r1, 1500, 100);
        write_fastq(&r2, 1500, 100);
        let out = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![r1, r2],
            output: out.clone(),
            working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        compress_cluster_paired(&cfg).unwrap();
        let bytes = std::fs::read(&out).unwrap();
        assert_eq!(bytes[3], 3, "archive_type byte must be 3 (cluster)");
        let hdr = crate::compression::FixedHeader::parse_v5(&bytes).unwrap();
        assert_eq!(hdr.archive_type, 3);
        // footer must carry mate-1 AND mate-2 entries
        let header_end = u32::from_le_bytes(bytes[4..8].try_into().unwrap()) as u64;
        let d = crate::compression::chunk_directory::read_v5_footer(&out, header_end).unwrap();
        assert!(d.entries.iter().any(|e| e.mate == 1), "has R1 (mate 1) entries");
        assert!(d.entries.iter().any(|e| e.mate == 2), "has R2 (mate 2) entries");
    }

    #[test]
    fn compress_cluster_paired_multichunk_increments_chunk_index() {
        let dir = tempfile::tempdir().unwrap();
        let r1 = dir.path().join("r1.fastq");
        let r2 = dir.path().join("r2.fastq");
        write_fastq(&r1, 4000, 80);
        write_fastq(&r2, 4000, 80);
        let out = dir.path().join("out.qz");
        let cfg = crate::cli::CompressConfig {
            input: vec![r1, r2],
            output: out.clone(),
            working_dir: dir.path().to_path_buf(),
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        // Force small chunks so the producer flushes mid-stream (>=2 chunks), exercising
        // the chunk-index increment + per-mate buffer clear. SAFETY: this var is only
        // set/removed here within this single test body.
        unsafe {
            std::env::set_var("QZ_CLUSTER_CHUNK_BYTES", "40000");
        }
        let res = compress_cluster_paired(&cfg);
        unsafe {
            std::env::remove_var("QZ_CLUSTER_CHUNK_BYTES");
        }
        res.unwrap();
        let bytes = std::fs::read(&out).unwrap();
        let header_end = u32::from_le_bytes(bytes[4..8].try_into().unwrap()) as u64;
        let d = crate::compression::chunk_directory::read_v5_footer(&out, header_end).unwrap();
        assert!(d.num_chunks >= 2, "expected multiple chunks, got {}", d.num_chunks);
        assert_eq!(d.num_reads, 4000, "num_reads == n_pairs");
        // Both mates' ClusteredSeq record_counts must sum to n_pairs across all chunks.
        let seq_sum = |mate: u8| -> u64 {
            d.entries
                .iter()
                .filter(|e| e.mate == mate && e.role == StreamRole::ClusteredSeq)
                .map(|e| e.record_count)
                .sum()
        };
        assert_eq!(seq_sum(1), 4000, "mate-1 ClusteredSeq sum == n_pairs");
        assert_eq!(seq_sum(2), 4000, "mate-2 ClusteredSeq sum == n_pairs");
    }

    #[test]
    fn reject_unsupported_accepts_two_inputs_rejects_three() {
        let two = crate::cli::CompressConfig {
            input: vec!["a.fastq".into(), "b.fastq".into()],
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        assert!(reject_unsupported(&two).is_ok(), "paired (2 inputs) now allowed");
        let three = crate::cli::CompressConfig {
            input: vec!["a".into(), "b".into(), "c".into()],
            cluster: Some(crate::cli::ClusterOptions::default()),
            ..crate::cli::CompressConfig::default()
        };
        assert!(reject_unsupported(&three).is_err(), "3+ inputs rejected");
    }
}
