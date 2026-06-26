//! Reference-direct paired-end sequence ENCODE (`--reference`, spec §4).
//!
//! The compress half: a single streaming pass maps each read against a
//! reference FASTA and encodes it as edits over the covered reference slices
//! (the "backing"); off-reference reads spill to a per-(chunk,mate) fallback
//! literal pool written inline into that chunk's payload. The reference is
//! never stored. The decode/verify half lives in [`super::decode`].
use crate::cli::CompressConfig;
use anyhow::{Result, bail};

use super::backing::IntervalMap;
use super::edits::encode_edits;
use super::format::{
    CODEC_BSC, CODEC_COLUMNAR, CODEC_FQZCOMP, CODEC_PACKED_CONSENSUS, GLOBAL_SENTINEL, RefRole,
    ref_role_to_stream,
};
use super::{backing, coverage, edits, fallback, mapping, positions, refmeta};
use crate::compression::chunk_directory::{
    ChunkDirEntry, GLOBAL_SENTINEL as DIR_GLOBAL_SENTINEL,
};
use crate::compression::dna_utils::write_varint;

/// Reject any configuration the reference path does not support. One `bail!`
/// per concern so later tasks can assert specific rejection messages. The gate
/// fires before any compression work begins, so a rejected config never writes
/// an output file.
pub(crate) fn reject_unsupported(a: &CompressConfig) -> Result<()> {
    // Shared multi-stream gating; reference allows only `Auto` quality and requires a
    // named output file (no stdout).
    crate::compression::reject_unsupported_common(
        a,
        &crate::compression::ModeCaps {
            mode: "reference",
            allows_explicit_fqzcomp: false,
            allows_stdout: false,
        },
    )?;
    // Reference-only extra: a prebuilt index is reserved/rejected in this release.
    if let Some(r) = &a.reference
        && r.reference_index.is_some()
    {
        bail!(
            "--reference-index not supported in this release (index is rebuilt from --reference)"
        );
    }
    Ok(())
}

// === Task 12 tunables (single source of truth) ===========================
/// Placement offset search half-window (bases) around the NAM's provisional
/// base index. Small: the NAM projection is already close; a few bases covers
/// the seed-projection slop.
const PLACE_WINDOW: u64 = 8;
/// Pairs sampled in Pass 0 to estimate the insert-size distribution.
const INSERT_SAMPLE: usize = 100_000;
/// Reference-mode chunk size (pairs). Deliberately 5× smaller than the shared
/// `CompressConfig` default (2.5M) used by the paired path: measured on chr20
/// NovaSeq, 500K cuts the in-flight working set ~4× (peak RSS 14.1 → 3.5 GB) and
/// runs ~7% faster (finer producer/consumer pipeline overlap) for only ~0.5%
/// larger archives. NOTE: those speed/RSS figures were measured on the
/// pre-worker-pool bespoke producer (commit 5b426fc); the current shared-pool
/// engine (in-flight ≈ `reference_window` × chunk) has not been re-measured, so
/// treat 500K as a validated ~0.5% ratio knob, not a validated speed setting.
/// The chunk only affects the ~17% BSC-coded side streams
/// (headers/positions/edits/fallback); qualities use fqzcomp, sub-blocked at
/// `quality_ctx_block_size` (500K records), so they're chunk-insensitive.
/// An explicit `CompressConfig::advanced.chunk_records` overrides this default; the
/// `QZ_REF_CHUNK` env var is a lower-priority experimental override (see
/// `resolve_ref_chunk_records`).
const REF_CHUNK_RECORDS: usize = 500_000;

/// Resolve the reference chunk size (records/chunk). Precedence, highest first:
///   1. An explicitly-configured `chunk_records` — anything other than the shared
///      `CompressConfig` default sentinel (2.5M). This makes the config field honest:
///      callers/tests that set `advanced.chunk_records` directly drive reference
///      chunking (it was previously a silent no-op — reference only read the env var).
///   2. The `QZ_REF_CHUNK` env var — experimental override (benchmarks/timing).
///   3. `REF_CHUNK_RECORDS` (500K) — reference's tuned production default (lower peak
///      RSS than the shared 2.5M default; not CLI-exposed).
///
/// Config wins over the env var so a caller that sets `chunk_records` is immune to a
/// `QZ_REF_CHUNK` that leaked into the process environment — closing the setenv/getenv
/// data race the test-side env guard used to work around. Result is clamped to ≥1.
fn resolve_ref_chunk_records(cfg_chunk_records: usize, env_qz_ref_chunk: Option<&str>) -> usize {
    let resolved = if cfg_chunk_records != crate::cli::AdvancedOptions::default().chunk_records {
        cfg_chunk_records
    } else {
        env_qz_ref_chunk
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(REF_CHUNK_RECORDS)
    };
    resolved.max(1)
}

#[cfg(test)]
mod chunk_records_resolution_tests {
    use super::{REF_CHUNK_RECORDS, resolve_ref_chunk_records};

    fn default_sentinel() -> usize {
        crate::cli::AdvancedOptions::default().chunk_records
    }

    #[test]
    fn explicit_config_wins_even_over_env() {
        // A config that differs from the default sentinel drives chunking, and is
        // immune to a QZ_REF_CHUNK that leaked into the environment.
        assert_eq!(resolve_ref_chunk_records(50, Some("9999")), 50);
        assert_eq!(resolve_ref_chunk_records(50, None), 50);
    }

    #[test]
    fn env_used_only_when_config_is_default() {
        assert_eq!(
            resolve_ref_chunk_records(default_sentinel(), Some("128")),
            128
        );
    }

    #[test]
    fn falls_back_to_ref_default() {
        // Config default + no/garbage env → reference's tuned default.
        assert_eq!(
            resolve_ref_chunk_records(default_sentinel(), None),
            REF_CHUNK_RECORDS
        );
        assert_eq!(
            resolve_ref_chunk_records(default_sentinel(), Some("not-a-number")),
            REF_CHUNK_RECORDS
        );
    }

    #[test]
    fn clamps_to_at_least_one() {
        assert_eq!(resolve_ref_chunk_records(0, None), 1);
        assert_eq!(resolve_ref_chunk_records(default_sentinel(), Some("0")), 1);
    }
}
/// Insert sizes above this (bases) are treated as discordant and excluded from
/// the mu/sigma estimate.
const INSERT_MAX: u64 = 2000;

/// Per-read placement k: max accepted Hamming substitutions. Scales with read
/// length (~10%), floored at 4 so very short reads still place.
#[inline]
fn place_k(read_len: usize) -> usize {
    (read_len / 10).max(4)
}

/// Result of diffing one read against the reference.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum Placed {
    Mapped {
        ref_id: u32,
        ref_pos: u64,
        is_rc: bool,
        read_len: u32,
        edits: Vec<edits::Edit>,
    },
    Fallback,
}

/// Diff a read against the normalized reference near its mapped genome position.
/// Small offset search (PLACE_WINDOW) + reconstruct-verify; any failure -> Fallback.
pub(crate) fn place_against_reference(
    refs: &[strobealign::io::fasta::RefSequence],
    ref_id: u32,
    proj_ref_start: u64,
    seq: &[u8],
    is_rc: bool,
) -> Placed {
    let idx = ref_id as usize;
    if idx >= refs.len() {
        return Placed::Fallback;
    }
    let chrom = &refs[idx].sequence;
    let k = place_k(seq.len());
    if let Some(p) = mapping::place_read(chrom, proj_ref_start, seq, is_rc, k, PLACE_WINDOW) {
        let s = p.base_index as usize;
        // place_read guarantees [s, s+len) is in-bounds of chrom
        let slice = &chrom[s..s + seq.len()];
        if edits::reconstruct_read(slice, &p.edits, is_rc) == seq {
            return Placed::Mapped {
                ref_id,
                ref_pos: p.base_index,
                is_rc,
                read_len: seq.len() as u32,
                edits: p.edits,
            };
        }
    }
    Placed::Fallback
}

/// Pass-0 planning prelude (carved out of the former inline `compress_reference`
/// preamble; the legacy path below still calls it, so it stays byte-identical).
///
/// Opens R1/R2 once, buffers the WHOLE `INSERT_SAMPLE`-pair prefix, profiles
/// `read_len` + insert μ/σ, and builds the (normalized) `Mapper`. Returns the
/// open readers + the whole buffered prefix so the caller (legacy consumer OR
/// the streaming producer) can feed prefix-then-readers, chunked at
/// `chunk_records`, decoding the files exactly once. `chunk_records` is taken
/// (not used here) only so the streaming caller and the legacy caller resolve it
/// identically.
/// Open a reference input reader, optionally bounded to a record-aligned byte range
/// (the NUMA byte-range compress worker passes its shard's range; the whole-file
/// entry passes `None`). Reference inputs are never gzipped (rejected upstream), so
/// the bounded path is plain-file `seek + take`.
fn open_ref_reader(
    path: &std::path::Path,
    range: Option<&crate::compression::ByteRange>,
) -> Result<crate::io::FastqReader<crate::io::fastq::FileReader>> {
    match range {
        Some(r) => crate::io::FastqReader::from_range(path, r.start, r.end, false),
        None => crate::io::FastqReader::from_path(path, false),
    }
}

/// `ranges` (the byte-range worker's shard) is `Some([r1, r2])` for a sharded worker
/// or `None` for a whole-file compress; when set, the readers are bounded to the
/// shard. Returns the open readers + the whole buffered prefix so the streaming
/// producer can feed prefix-then-readers.
#[allow(clippy::type_complexity)]
pub(crate) fn reference_prelude(
    a: &CompressConfig,
    ranges: Option<&[crate::compression::ByteRange]>,
) -> Result<(
    super::mapping::Mapper,
    usize, // read_len
    crate::io::FastqReader<crate::io::fastq::FileReader>,
    crate::io::FastqReader<crate::io::fastq::FileReader>,
    Vec<(crate::io::FastqRecord, crate::io::FastqRecord)>, // whole sampled prefix
    usize,                                                 // qctx_block_size
)> {
    use crate::io::FastqRecord;
    use super::mapping::{Mapper, ReferenceSource};
    use rayon::prelude::*;

    reject_unsupported(a)?;

    let reference_path = a
        .reference
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("reference mode requires --reference"))?
        .reference
        .clone();
    let r1_path = a.input[0].clone();
    let r2_path = a.input[1].clone();
    let threads = if a.threads == 0 {
        rayon::current_num_threads()
    } else {
        a.threads
    };

    if let Some(rs) = ranges
        && rs.len() != 2
    {
        bail!("reference prelude: paired worker needs 2 byte ranges, got {}", rs.len());
    }
    // ================= Pass 0 — open once, profile, build mapper =============
    // Open R1/R2 once (bounded to the shard's range when sharded); buffer the first
    // INSERT_SAMPLE pairs (the "prefix") to profile read_len + insert μ/σ. The prefix
    // is then fed to the chunk loop so the files are decoded exactly once.
    let mut r1 = open_ref_reader(&r1_path, ranges.map(|rs| &rs[0]))?;
    let mut r2 = open_ref_reader(&r2_path, ranges.map(|rs| &rs[1]))?;

    let mut prefix: Vec<(FastqRecord, FastqRecord)> = Vec::new();
    while prefix.len() < INSERT_SAMPLE {
        let (a1, a2) = (r1.next()?, r2.next()?);
        match (a1, a2) {
            (Some(x), Some(y)) => prefix.push((x, y)),
            (None, None) => break,
            _ => bail!("reference mode: R1/R2 have unequal read counts"),
        }
    }
    let read_len = prefix
        .first()
        .map(|(rec1, _)| rec1.sequence.len())
        .ok_or_else(|| anyhow::anyhow!("empty R1 input"))?;
    if read_len == 0 {
        bail!("reference mode: first read has zero length");
    }

    let fast = a
        .reference
        .as_ref()
        .map(|r| r.reference_fast)
        .unwrap_or(false);
    let mut mapper = Mapper::build(
        ReferenceSource::FromFasta(reference_path.clone()),
        read_len,
        threads,
        fast,
        a.require_prebuilt_index,
    )?;

    // Estimate insert size from the buffered prefix (same-ref_id pairs only;
    // deterministic — fixed before the chunk loop, no RNG). Mapping each prefix
    // pair is independent, so run the sample in parallel: this is a SERIAL startup
    // section (it gates the chunk loop), so cutting its latency helps the wall at
    // every core count. mean/var are order-independent → mu/sigma byte-identical.
    let (mu, sigma) = {
        let inserts: Vec<u64> = prefix
            .par_iter()
            .map(|(rec1, rec2)| -> Result<Option<u64>> {
                let (m1, m2) = mapper.map_pair(&rec1.sequence, &rec2.sequence)?;
                Ok(match (m1, m2) {
                    (Some(m1), Some(m2)) if m1.ref_id == m2.ref_id => {
                        let ins = (m1.ref_start as i64 - m2.ref_start as i64).unsigned_abs();
                        (ins < INSERT_MAX).then_some(ins)
                    }
                    _ => None,
                })
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect();
        if inserts.len() < 100 {
            (300.0f32, 100.0f32)
        } else {
            let n = inserts.len() as f64;
            let mean = inserts.iter().map(|&x| x as f64).sum::<f64>() / n;
            let var = inserts
                .iter()
                .map(|&x| {
                    let d = x as f64 - mean;
                    d * d
                })
                .sum::<f64>()
                / n;
            (mean as f32, var.sqrt() as f32)
        }
    };
    mapper = mapper.with_insert(mu, sigma);

    let qctx_block_size = a.advanced.quality_ctx_block_size;
    Ok((mapper, read_len, r1, r2, prefix, qctx_block_size))
}

/// Single-end Pass-0 prelude (spec §4 ②). Mirrors [`reference_prelude`] but opens
/// ONE reader, buffers a prefix of single records to profile `read_len`, and
/// builds the `Mapper` WITHOUT `with_insert` — single-end mapping ignores μ/σ. No
/// insert-size sampling pass. Returns the open reader + the whole buffered prefix
/// so the producer feeds prefix-then-reader, decoding the file exactly once.
#[allow(clippy::type_complexity)]
pub(crate) fn reference_prelude_single(
    a: &CompressConfig,
    range: Option<&crate::compression::ByteRange>,
) -> Result<(
    super::mapping::Mapper,
    usize, // read_len
    crate::io::FastqReader<crate::io::fastq::FileReader>,
    Vec<crate::io::FastqRecord>, // whole sampled prefix
    usize,                       // qctx_block_size
)> {
    use crate::io::FastqRecord;
    use super::mapping::{Mapper, ReferenceSource};

    reject_unsupported(a)?;

    let reference_path = a
        .reference
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("reference mode requires --reference"))?
        .reference
        .clone();
    let r1_path = a.input[0].clone();
    let threads = if a.threads == 0 {
        rayon::current_num_threads()
    } else {
        a.threads
    };

    let mut reader = open_ref_reader(&r1_path, range)?;
    let mut prefix: Vec<FastqRecord> = Vec::new();
    while prefix.len() < INSERT_SAMPLE {
        match reader.next()? {
            Some(x) => prefix.push(x),
            None => break,
        }
    }
    let read_len = prefix
        .first()
        .map(|rec| rec.sequence.len())
        .ok_or_else(|| anyhow::anyhow!("empty input"))?;
    if read_len == 0 {
        bail!("reference mode: first read has zero length");
    }

    let fast = a
        .reference
        .as_ref()
        .map(|r| r.reference_fast)
        .unwrap_or(false);
    // No `with_insert`: map_single never reads μ/σ. The mapper keeps its
    // build-time (300, 100) defaults, which the single-end path ignores.
    let mapper = Mapper::build(
        ReferenceSource::FromFasta(reference_path.clone()),
        read_len,
        threads,
        fast,
        a.require_prebuilt_index,
    )?;

    let qctx_block_size = a.advanced.quality_ctx_block_size;
    Ok((mapper, read_len, reader, prefix, qctx_block_size))
}

/// Resolve the reference chunk size for a config (reads `QZ_REF_CHUNK`). Public to
/// the crate so the streaming entry point resolves it identically to the legacy
/// path. See `resolve_ref_chunk_records` for the precedence.
pub(crate) fn resolve_ref_chunk_records_cfg(a: &CompressConfig) -> usize {
    let env_qz_ref_chunk = std::env::var("QZ_REF_CHUNK").ok();
    resolve_ref_chunk_records(a.advanced.chunk_records, env_qz_ref_chunk.as_deref())
}

/// Map + diff one chunk of pairs against the reference (parallel, order-preserving).
/// Byte-identical to the former inline producer map/diff: a parallel `map_pair`
/// pass collected into a Vec, then a parallel `place_against_reference` diff zipped
/// 1:1 with the maps. The result is a `Vec<(Placed, Placed)>` in input order.
pub(crate) fn map_and_diff(
    mapper: &mapping::Mapper,
    chunk: &[(crate::io::FastqRecord, crate::io::FastqRecord)],
) -> Result<Vec<(Placed, Placed)>> {
    use rayon::prelude::*;
    // Parallel map (order-preserving).
    let maps: Vec<(Option<mapping::Mapping>, Option<mapping::Mapping>)> = chunk
        .par_iter()
        .map(|(a, b)| mapper.map_pair(&a.sequence, &b.sequence))
        .collect::<Result<Vec<_>>>()?;
    // Parallel diff (order-preserving).
    let placed: Vec<(Placed, Placed)> = chunk
        .par_iter()
        .zip(maps.par_iter())
        .map(|((a, b), (m1, m2))| {
            (
                diff_one(mapper, &a.sequence, *m1),
                diff_one(mapper, &b.sequence, *m2),
            )
        })
        .collect();
    Ok(placed)
}

/// Map + diff one chunk of SINGLE-END records against the reference (parallel,
/// order-preserving). Single-end analog of [`map_and_diff`]: a parallel
/// `map_single` pass collected into a Vec, then a parallel `diff_one` zipped 1:1.
/// Result is `Vec<Placed>` in input order.
pub(crate) fn map_and_diff_single(
    mapper: &mapping::Mapper,
    chunk: &[crate::io::FastqRecord],
) -> Result<Vec<Placed>> {
    use rayon::prelude::*;
    let maps: Vec<Option<mapping::Mapping>> = chunk
        .par_iter()
        .map(|r| mapper.map_single(&r.sequence))
        .collect::<Result<Vec<_>>>()?;
    let placed: Vec<Placed> = chunk
        .par_iter()
        .zip(maps.par_iter())
        .map(|(r, m)| diff_one(mapper, &r.sequence, *m))
        .collect();
    Ok(placed)
}

/// Single-pass reference-direct compress (plan Tasks 13/14/15, spec §4.1).
///
/// Reference compress is driven by the shared block-granular worker pool
/// (`spawn_worker_pool`); see [`super::stream::compress_reference_streaming`].
/// The alignment front end (mapping → diff → columns → consensus) stays outside
/// the shared block-streaming core: a thin chunk-ordered writer reorders the
/// per-chunk `encode_chunk` outputs and streams the whole-archive globals in its
/// tail.
pub fn compress_reference(a: &CompressConfig) -> Result<()> {
    super::stream::compress_reference_streaming(a, None)
}

/// Single-pass SINGLE-END reference-direct compress (archive_type 4). Mirrors
/// [`compress_reference`] but drives one reader through the single-end streaming
/// pipeline.
pub fn compress_reference_single(a: &CompressConfig) -> Result<()> {
    super::stream_single::compress_reference_single_streaming(a, None)
}

/// Byte-range reference compress for the NUMA shard worker: compress only the records
/// in the per-mate `ranges` (record-aligned by the splitter) into a reference archive
/// (paired `archive_type 2` for 2 ranges, single-end `archive_type 4` for 1). The
/// resulting part archive is byte-range-disjoint from the other shards' and is
/// stitched by `merge_reference_archives_to_path` (front headers match — reference's
/// header is fixed metadata, not data-derived). `a.output` is the worker's part path.
pub(crate) fn compress_reference_byte_range(
    a: &CompressConfig,
    ranges: &[crate::compression::ByteRange],
) -> Result<()> {
    match ranges.len() {
        2 => super::stream::compress_reference_streaming(a, Some(ranges)),
        1 => super::stream_single::compress_reference_single_streaming(a, Some(&ranges[0])),
        n => anyhow::bail!("reference byte-range compress expects 1 or 2 ranges, got {n}"),
    }
}


/// `G` (gap-merge) cap for `CoverageMap::finalize` (spec §9): `min(read_len, CAP)`.
const GAP_MERGE_CAP: u64 = 1000;

/// Diff one read against the mapper's normalized reference. `m == None` →
/// `Placed::Fallback`; else `place_against_reference` near the mapped position.
fn diff_one(mapper: &mapping::Mapper, seq: &[u8], m: Option<mapping::Mapping>) -> Placed {
    match m {
        None => Placed::Fallback,
        Some(mp) => place_against_reference(
            mapper.references(),
            mp.ref_id,
            mp.ref_start,
            seq,
            mp.is_revcomp,
        ),
    }
}

/// Mark coverage for a `Placed::Mapped` over `[ref_pos, ref_pos + read_len)`;
/// no-op for `Fallback`. `read_len` is the stored placed-read length (the
/// original read length), so the marked span is exactly the read's footprint.
pub(crate) fn mark_cov(cov: &mut coverage::CoverageMap, placed: &Placed) {
    if let Placed::Mapped {
        ref_id,
        ref_pos,
        read_len,
        ..
    } = placed
    {
        cov.mark(*ref_id, *ref_pos, *read_len as u64);
    }
}

/// Write a fully-built role payload to the streaming writer, record a
/// `ChunkDirEntry` from the current offset, and advance `*cur`.
#[allow(clippy::too_many_arguments)]
fn push_entry_stream(
    out: &mut dyn std::io::Write,
    cur: &mut u64,
    entries: &mut Vec<ChunkDirEntry>,
    chunk_index: u32,
    role: RefRole,
    mate: u8,
    codec: u8,
    payload: &[u8],
    record_count: u64,
) -> Result<()> {
    out.write_all(payload)?;
    entries.push(ChunkDirEntry {
        chunk_index,
        role: ref_role_to_stream(role),
        mate,
        codec,
        offset: *cur,
        length: payload.len() as u64,
        record_count,
    });
    *cur += payload.len() as u64;
    Ok(())
}

/// Append a fully-built role payload to a per-chunk `buf`, recording a `ChunkDirEntry`
/// whose `offset` is RELATIVE to the start of `buf`. The in-order writer shifts it
/// by the chunk's base offset before pushing to the archive directory. Returns
/// `Result` only to keep the per-role call sites identical to the streaming path.
#[allow(clippy::too_many_arguments)]
fn push_entry_buf(
    buf: &mut Vec<u8>,
    entries: &mut Vec<ChunkDirEntry>,
    chunk_index: u32,
    role: RefRole,
    mate: u8,
    codec: u8,
    payload: &[u8],
    record_count: u64,
) -> Result<()> {
    let offset = buf.len() as u64;
    buf.extend_from_slice(payload);
    entries.push(ChunkDirEntry {
        chunk_index,
        role: ref_role_to_stream(role),
        mate,
        codec,
        offset,
        length: payload.len() as u64,
        record_count,
    });
    Ok(())
}

/// Record a global `ChunkDirEntry` whose bytes were ALREADY streamed straight through
/// `out` (e.g. the packed backing / N-bitmap). `byte_len` is the exact number of
/// bytes written for the role; `cur` is advanced.
#[allow(clippy::too_many_arguments)]
fn push_global_entry(
    cur: &mut u64,
    entries: &mut Vec<ChunkDirEntry>,
    role: RefRole,
    codec: u8,
    byte_len: u64,
    record_count: u64,
) {
    entries.push(ChunkDirEntry {
        chunk_index: DIR_GLOBAL_SENTINEL,
        role: ref_role_to_stream(role),
        mate: 0,
        codec,
        offset: *cur,
        length: byte_len,
        record_count,
    });
    *cur += byte_len;
}

/// A fully-encoded chunk: the concatenated byte payload (mate-1 roles followed by
/// mate-2 roles, including each mate's inline FallbackPool entry), and the per-role
/// `ChunkDirEntry`s with `offset`s RELATIVE to `bytes`. Self-contained — fallback
/// literals are written inline into `bytes`, so nothing escapes the chunk.
pub(crate) struct EncodedChunk {
    pub(crate) bytes: Vec<u8>,
    pub(crate) entries: Vec<ChunkDirEntry>,
    /// Per-mate decoded FASTQ output byte length for this chunk (`[R1, R2]`; mate-2 is
    /// 0 for single-end reference). Reference is always lossless FASTQ (`reject_unsupported`
    /// forbids `--fasta`/`--no-quality`/lossy), so each read contributes
    /// `fastq_output_len(false, id, seq, qual)` with `qual == seq` — byte-exact to what the
    /// decoder writes. Drives the `ChunkDecodedSizes` global the writer emits so reference
    /// NUMA decode can direct-write into pre-sized outputs.
    pub(crate) decoded_output_bytes: [u64; 2],
}

/// Emit the one whole-archive `ChunkDecodedSizes` global block (shared by the paired
/// and single-end reference writers). `per_chunk[c] = [r1_bytes, r2_bytes]` in chunk
/// order; `n_mates` is 2 (paired) or 1 (single-end, mate-2 ignored). Writes a row-major
/// `sizes[chunk*n_mates + mate]` payload as a `StreamRole::ChunkDecodedSizes`
/// `GLOBAL_SENTINEL` block via the shared `write_role_blocks` — byte-framed identically
/// to single/paired so `read_chunk_layout`'s type-agnostic parser picks it up.
pub(crate) fn emit_decoded_sizes_global<W: std::io::Write + ?Sized>(
    out: &mut W,
    n_mates: u8,
    num_chunks: u32,
    per_chunk: &[[u64; 2]],
    cur: &mut u64,
    entries: &mut Vec<ChunkDirEntry>,
) -> Result<()> {
    debug_assert_eq!(per_chunk.len(), num_chunks as usize);
    let mut sizes = Vec::with_capacity(per_chunk.len() * n_mates as usize);
    for d in per_chunk {
        for m in 0..n_mates as usize {
            sizes.push(d[m]);
        }
    }
    let payload = crate::compression::chunk_directory::encode_decoded_sizes(n_mates, num_chunks, &sizes);
    crate::compression::compress_impl::write_role_blocks(
        &[(num_chunks, payload)],
        crate::compression::chunk_directory::StreamRole::ChunkDecodedSizes,
        0, // mate (the global is always mate 0)
        0, // codec (raw metadata)
        crate::compression::chunk_directory::GLOBAL_SENTINEL,
        out,
        cur,
        entries,
    )
}

/// Decoded FASTQ output bytes for one reference record. Reference is always lossless
/// FASTQ (`@id\nseq\n+\nqual\n`), so this is byte-exact to what the decoder writes —
/// the single source feeding the per-chunk `ChunkDecodedSizes` table.
#[inline]
fn ref_record_decoded_len(rec: &crate::io::FastqRecord) -> u64 {
    let qual = rec.quality.as_deref().map_or(0, |q| q.len());
    crate::compression::decompress_impl::fastq_output_len(false, rec.id.len(), rec.sequence.len(), qual)
}

/// Encode one chunk (both mates) into a self-contained [`EncodedChunk`]. A PURE
/// function of its inputs — no shared writer/spool — so chunks encode in parallel
/// and the caller writes them to the archive IN ORDER. The byte layout and the
/// relative-offset entries are identical to the former two sequential
/// `encode_chunk_mate` calls, so the archive is byte-identical.
pub(crate) fn encode_chunk(
    chunk_index: u32,
    chunk: &[(crate::io::FastqRecord, crate::io::FastqRecord)],
    placed: &[(Placed, Placed)],
    qctx_block_size: usize,
) -> Result<EncodedChunk> {
    let n = placed.len();
    // Headers + qualities are BORROWED from `chunk` (which outlives this fn),
    // not cloned. The columnar header codec reads only the id bytes and fqzcomp
    // reads only the quality bytes, both as `&[u8]` — so cloning every id and
    // quality into owned Vecs per chunk was ~176 MB/chunk of dead DRAM copy that
    // was freed immediately. On a memory-bandwidth-bound encoder, removing bus
    // traffic is the relevant lever. Output is byte-identical (the codecs see the
    // same bytes either way).
    let mut r1_hdrs: Vec<&[u8]> = Vec::with_capacity(n);
    let mut r2_hdrs: Vec<&[u8]> = Vec::with_capacity(n);
    let mut r1_quals: Vec<&[u8]> = Vec::with_capacity(n);
    let mut r2_quals: Vec<&[u8]> = Vec::with_capacity(n);
    // R1 anchor per pair (mate-2 position encoder input); None for unmapped R1.
    let mut r1_anchor: Vec<Option<(u32, u64)>> = Vec::with_capacity(n);
    for ((rec1, rec2), (p1, _p2)) in chunk.iter().zip(placed.iter()) {
        let q1 = rec1
            .quality
            .as_deref()
            .ok_or_else(|| anyhow::anyhow!("reference mode: R1 record missing quality"))?;
        let q2 = rec2
            .quality
            .as_deref()
            .ok_or_else(|| anyhow::anyhow!("reference mode: R2 record missing quality"))?;
        r1_hdrs.push(rec1.id.as_slice());
        r2_hdrs.push(rec2.id.as_slice());
        r1_quals.push(q1);
        r2_quals.push(q2);
        r1_anchor.push(match p1 {
            Placed::Mapped {
                ref_id, ref_pos, ..
            } => Some((*ref_id, *ref_pos)),
            Placed::Fallback => None,
        });
    }

    let mut bytes = Vec::new();
    let mut entries: Vec<ChunkDirEntry> = Vec::new();
    // Pass the pair-slice directly + a mate selector instead of cloning each
    // Placed (and its inner Vec<Edit>) into flat per-mate Vecs.
    // Mate 1's serialized columnar header length is the R1-columnar proxy threaded
    // into mate 2 for the predict-and-skip (≈ what R2's independent columnar costs).
    let r1_header_proxy = encode_mate_buf(
        &mut bytes,
        &mut entries,
        chunk_index,
        1,
        placed,
        None,
        chunk,
        RefRole::R1Headers,
        RefRole::R1Qual,
        &r1_hdrs,
        None,
        None,
        &r1_quals,
        qctx_block_size,
    )?;
    encode_mate_buf(
        &mut bytes,
        &mut entries,
        chunk_index,
        2,
        placed,
        Some(r1_anchor.as_slice()),
        chunk,
        RefRole::R2HeadersIndep,
        RefRole::R2Qual,
        &r2_hdrs,
        // R1 ids enable the delta-vs-R1 header candidate (reference's R2 ids are R1
        // ids with the mate suffix flipped → tiny delta, ~5% smaller archive).
        Some(r1_hdrs.as_slice()),
        // R1-columnar proxy → skip the R2 columnar BWT when delta provably wins.
        Some(r1_header_proxy),
        &r2_quals,
        qctx_block_size,
    )?;
    // Per-mate decoded output bytes (byte-exact to the decoder's FASTQ formatter).
    let mut decoded_output_bytes = [0u64; 2];
    for (r1, r2) in chunk.iter() {
        decoded_output_bytes[0] += ref_record_decoded_len(r1);
        decoded_output_bytes[1] += ref_record_decoded_len(r2);
    }
    Ok(EncodedChunk { bytes, entries, decoded_output_bytes })
}

/// Encode one SINGLE-END chunk (mate = 1) into a self-contained [`EncodedChunk`].
/// Single-end analog of [`encode_chunk`]: a PURE function of its inputs (no shared
/// writer), so chunks encode in parallel and the caller writes them IN ORDER.
/// Headers/qualities are BORROWED from `chunk` (not cloned) and handed to
/// [`encode_mate_single`].
pub(crate) fn encode_chunk_single(
    chunk_index: u32,
    chunk: &[crate::io::FastqRecord],
    placed: &[Placed],
    qctx_block_size: usize,
) -> Result<EncodedChunk> {
    let n = placed.len();
    let mut hdrs: Vec<&[u8]> = Vec::with_capacity(n);
    let mut quals: Vec<&[u8]> = Vec::with_capacity(n);
    for rec in chunk.iter() {
        let q = rec
            .quality
            .as_deref()
            .ok_or_else(|| anyhow::anyhow!("reference mode: record missing quality"))?;
        hdrs.push(rec.id.as_slice());
        quals.push(q);
    }
    let mut bytes = Vec::new();
    let mut entries: Vec<ChunkDirEntry> = Vec::new();
    encode_mate_single(
        &mut bytes,
        &mut entries,
        chunk_index,
        placed,
        chunk,
        &hdrs,
        &quals,
        qctx_block_size,
    )?;
    // Single-end: mate-1 decoded output bytes; mate-2 is 0 (no R2).
    let mut decoded0 = 0u64;
    for rec in chunk.iter() {
        decoded0 += ref_record_decoded_len(rec);
    }
    Ok(EncodedChunk { bytes, entries, decoded_output_bytes: [decoded0, 0] })
}

/// Encode one (chunk, mate) into the per-chunk `buf` (chunk-relative `ChunkDirEntry`
/// offsets), collecting fallback-read literals into `fallbacks` (call order). The
/// 7 per-chunk sequence roles + header + quality roles. No shared writer/spool.
///
/// `placed[i]` is the (R1,R2) diff-result pair for read-pair `i` (N pairs total);
/// this mate's `Placed` is `.0` (mate 1) or `.1` (mate 2). `r1_anchor` is
/// the same slice the mate-2 position encoder needs: for mate 2 it must be
/// aligned 1:1 to the MAPPED R2 reads (one anchor per mapped R2). For mate 1 it
/// is ignored (pass all-None). `chunk` carries the original records (for the
/// fallback literal bytes + the header/qual encoder's sequences).
#[allow(clippy::too_many_arguments)]
fn encode_mate_buf(
    buf: &mut Vec<u8>,
    entries: &mut Vec<ChunkDirEntry>,
    chunk_index: u32,
    mate: u8,
    // The full pair-slice; this mate's Placed is selected per-read by `mate`
    // (`.0` for mate 1, `.1` for mate 2). Borrowed, never cloned.
    placed: &[(Placed, Placed)],
    // None for mate 1 (anchors are unused there); Some(per-pair R1 anchors) for
    // mate 2's delta position encoder. Replaces the old all-None throwaway Vec.
    r1_anchor_full: Option<&[Option<(u32, u64)>]>,
    chunk: &[(crate::io::FastqRecord, crate::io::FastqRecord)],
    header_role: RefRole,
    qual_role: RefRole,
    headers: &[&[u8]],
    // Some(R1 ids) for mate 2 → enables the delta-vs-R1 header candidate; None for mate 1.
    r1_headers: Option<&[&[u8]]>,
    // Some(mate-1's serialized columnar header length) for mate 2 → the R1-columnar
    // proxy that drives the predict-and-skip of the R2 columnar candidate; None for mate 1.
    r1_header_proxy: Option<usize>,
    quals: &[&[u8]],
    qctx_block_size: usize,
) -> Result<usize> {
    let n_chunk = placed.len() as u64;
    debug_assert_eq!(headers.len(), placed.len());
    debug_assert_eq!(quals.len(), placed.len());
    debug_assert_eq!(chunk.len(), placed.len());
    // This (chunk,mate)'s fallback literals, in read order. Written inline below
    // as this group's own FallbackPool block-stream entry (omitted when empty).
    // BORROWED from `chunk` (write_fallback_pool only reads the literals), so the
    // unmapped reads' sequences aren't copied into owned Vecs just to be BSC'd.
    let mut fallbacks: Vec<&[u8]> = Vec::new();

    // --- 1. MappedFlags: 1 bit/read, LSB-first, set = mapped (record_count = N).
    let mut mapped_flags = vec![0u8; placed.len().div_ceil(8)];
    for (i, pp) in placed.iter().enumerate() {
        let p = if mate == 1 { &pp.0 } else { &pp.1 };
        if matches!(p, Placed::Mapped { .. }) {
            mapped_flags[i / 8] |= 1u8 << (i % 8);
        }
    }
    let mut s = Vec::new();
    write_bsc_role(&mut s, &mapped_flags, n_chunk)?;
    push_entry_buf(
        buf,
        entries,
        chunk_index,
        RefRole::MappedFlags,
        mate,
        CODEC_BSC,
        &s,
        n_chunk,
    )?;

    // --- 2. Gather mapped-read data (in read order): positions, strands,
    // read_len, edits. Also the per-mapped-R2 r1_anchor for the mate-2 encoder.
    let mut mapped_positions: Vec<(u32, u64)> = Vec::new();
    let mut mapped_r1_anchor: Vec<Option<(u32, u64)>> = Vec::new();
    let mut strand_bits: Vec<bool> = Vec::new();
    let mut read_lens: Vec<u64> = Vec::new();
    // Edits are BORROWED from `placed` (alive for this call), not cloned —
    // encode_edits only reads them.
    let mut edit_recs: Vec<&[edits::Edit]> = Vec::new();
    for (i, pp) in placed.iter().enumerate() {
        let p = if mate == 1 { &pp.0 } else { &pp.1 };
        if let Placed::Mapped {
            ref_id,
            ref_pos,
            is_rc,
            edits,
            read_len,
        } = p
        {
            mapped_positions.push((*ref_id, *ref_pos));
            // Only mate 2 consumes mapped_r1_anchor (encode_positions_mate2);
            // mate 1 passes None and this stays empty (was filled-but-unused).
            if let Some(anchors) = r1_anchor_full {
                mapped_r1_anchor.push(anchors[i]);
            }
            strand_bits.push(*is_rc);
            read_lens.push(*read_len as u64);
            edit_recs.push(edits.as_slice());
        }
    }
    let m = mapped_positions.len() as u64;

    // --- 2a. Positions (record_count = M). Mate 1: absolute (ref_id, ref_pos).
    // Mate 2: form-bits + delta/escape vs the aligned R1 anchor.
    let positions = if mate == 2 {
        positions::encode_positions_mate2(&mapped_r1_anchor, &mapped_positions)
    } else {
        positions::encode_positions_mate1(mapped_positions.iter().copied())
    };
    let mut s = Vec::new();
    write_bsc_role(&mut s, &positions, m)?;
    push_entry_buf(
        buf,
        entries,
        chunk_index,
        RefRole::Positions,
        mate,
        CODEC_BSC,
        &s,
        m,
    )?;

    // --- 2b. Strands: 1 bit/mapped read (record_count = M).
    let mut strand_bytes = vec![0u8; strand_bits.len().div_ceil(8)];
    for (i, &b) in strand_bits.iter().enumerate() {
        if b {
            strand_bytes[i / 8] |= 1u8 << (i % 8);
        }
    }
    let mut s = Vec::new();
    write_bsc_role(&mut s, &strand_bytes, m)?;
    push_entry_buf(
        buf,
        entries,
        chunk_index,
        RefRole::Strands,
        mate,
        CODEC_BSC,
        &s,
        m,
    )?;

    // --- 2c. ReadLen: varint/mapped read (record_count = M).
    let mut read_lens_raw = Vec::new();
    for rl in &read_lens {
        write_varint(&mut read_lens_raw, *rl);
    }
    let mut s = Vec::new();
    write_bsc_role(&mut s, &read_lens_raw, m)?;
    push_entry_buf(
        buf,
        entries,
        chunk_index,
        RefRole::ReadLen,
        mate,
        CODEC_BSC,
        &s,
        m,
    )?;

    // --- 2d. EditCount / EditPos / EditBase.
    let enc = encode_edits(&edit_recs);
    let total_edits: u64 = edit_recs.iter().map(|e| e.len() as u64).sum();
    let mut s = Vec::new();
    write_bsc_role(&mut s, &enc.edit_count, m)?;
    push_entry_buf(
        buf,
        entries,
        chunk_index,
        RefRole::EditCount,
        mate,
        CODEC_BSC,
        &s,
        m,
    )?;
    let mut s = Vec::new();
    write_bsc_role(&mut s, &enc.edit_pos, total_edits)?;
    push_entry_buf(
        buf,
        entries,
        chunk_index,
        RefRole::EditPos,
        mate,
        CODEC_BSC,
        &s,
        total_edits,
    )?;
    let mut s = Vec::new();
    write_bsc_role(&mut s, &enc.edit_base, total_edits)?;
    push_entry_buf(
        buf,
        entries,
        chunk_index,
        RefRole::EditBase,
        mate,
        CODEC_BSC,
        &s,
        total_edits,
    )?;

    // --- 3. Fallback reads → collect their literal sequence in read order. The
    // literal is the ORIGINAL read sequence (R1 for mate 1, R2 for mate 2).
    for (i, pp) in placed.iter().enumerate() {
        let p = if mate == 1 { &pp.0 } else { &pp.1 };
        if matches!(p, Placed::Fallback) {
            let seq: &[u8] = if mate == 2 {
                &chunk[i].1.sequence
            } else {
                &chunk[i].0.sequence
            };
            fallbacks.push(seq);
        }
    }

    // --- 3a. Write this (chunk,mate)'s FallbackPool inline (omitted when empty,
    // matching the validator's "present iff N−M > 0" contract). The literal count
    // F = N − M is implied by the entry's record_count; decode reads it on demand.
    if !fallbacks.is_empty() {
        let mut pool = Vec::new();
        let rc = fallback::write_fallback_pool(&mut pool, &fallbacks, REF_POOL_TARGET)?;
        push_entry_buf(
            buf,
            entries,
            chunk_index,
            RefRole::FallbackPool,
            mate,
            CODEC_BSC,
            &pool,
            rc,
        )?;
    }

    // --- 4. Headers + quality. Reference mode encodes headers columnar and
    // quality via fqzcomp (sub-blocked, cap-enforced). Neither codec reads the
    // read sequence, so we encode the original headers/qualities
    // directly instead of reconstructing every mapped read into a throwaway
    // record. `encode_ref_headers_qual` produces byte-identical h_blocks/q_blocks
    // to the old `encode_headers_qual_chunk` path while skipping the reconstruct,
    // the per-record clones, and the records_to_streams work (all discarded here).
    let (h_role, h_codec, h_blocks, q_blocks) = encode_ref_headers_qual(
        headers,
        quals,
        qctx_block_size,
        r1_headers,
        r1_header_proxy,
        header_role,
    )?;

    let mut hs = Vec::new();
    crate::compression::paired::streams::write_block_stream(&mut hs, &h_blocks);
    // mate 1's columnar header length is the R1-columnar proxy returned to the caller.
    let header_serialized_len = hs.len();
    push_entry_buf(
        buf,
        entries,
        chunk_index,
        h_role,
        mate,
        h_codec,
        &hs,
        n_chunk,
    )?;

    let mut qs = Vec::new();
    crate::compression::paired::streams::write_block_stream(&mut qs, &q_blocks);
    push_entry_buf(
        buf,
        entries,
        chunk_index,
        qual_role,
        mate,
        CODEC_FQZCOMP,
        &qs,
        n_chunk,
    )?;

    Ok(header_serialized_len)
}

/// Encode one SINGLE-END chunk's reads (mate = 1) into the per-chunk `buf`
/// (chunk-relative `ChunkDirEntry` offsets). The single-end analog of
/// [`encode_mate_buf`]: identical role sequence + leaf codecs, but `placed` is a
/// flat `&[Placed]` (no pair selector), positions are ALWAYS absolute
/// (`encode_positions_mate1`), headers are ALWAYS independent columnar
/// (`r1_headers=None`, `r1_header_proxy=None` → no delta candidate), and every
/// entry is tagged `mate = 1`. Fallback literals come from `records[i].sequence`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn encode_mate_single(
    buf: &mut Vec<u8>,
    entries: &mut Vec<ChunkDirEntry>,
    chunk_index: u32,
    placed: &[Placed],
    records: &[crate::io::FastqRecord],
    headers: &[&[u8]],
    quals: &[&[u8]],
    qctx_block_size: usize,
) -> Result<()> {
    let mate: u8 = 1;
    let n_chunk = placed.len() as u64;
    debug_assert_eq!(headers.len(), placed.len());
    debug_assert_eq!(quals.len(), placed.len());
    debug_assert_eq!(records.len(), placed.len());

    // --- 1. MappedFlags: 1 bit/read, LSB-first, set = mapped (record_count = N).
    let mut mapped_flags = vec![0u8; placed.len().div_ceil(8)];
    for (i, p) in placed.iter().enumerate() {
        if matches!(p, Placed::Mapped { .. }) {
            mapped_flags[i / 8] |= 1u8 << (i % 8);
        }
    }
    let mut s = Vec::new();
    write_bsc_role(&mut s, &mapped_flags, n_chunk)?;
    push_entry_buf(buf, entries, chunk_index, RefRole::MappedFlags, mate, CODEC_BSC, &s, n_chunk)?;

    // --- 2. Gather mapped-read data in read order.
    let mut mapped_positions: Vec<(u32, u64)> = Vec::new();
    let mut strand_bits: Vec<bool> = Vec::new();
    let mut read_lens: Vec<u64> = Vec::new();
    let mut edit_recs: Vec<&[edits::Edit]> = Vec::new();
    for p in placed.iter() {
        if let Placed::Mapped { ref_id, ref_pos, is_rc, edits, read_len } = p {
            mapped_positions.push((*ref_id, *ref_pos));
            strand_bits.push(*is_rc);
            read_lens.push(*read_len as u64);
            edit_recs.push(edits.as_slice());
        }
    }
    let m = mapped_positions.len() as u64;

    // --- 2a. Positions (record_count = M): ALWAYS absolute (ref_id, ref_pos).
    let positions = positions::encode_positions_mate1(mapped_positions.iter().copied());
    let mut s = Vec::new();
    write_bsc_role(&mut s, &positions, m)?;
    push_entry_buf(buf, entries, chunk_index, RefRole::Positions, mate, CODEC_BSC, &s, m)?;

    // --- 2b. Strands: 1 bit/mapped read (record_count = M).
    let mut strand_bytes = vec![0u8; strand_bits.len().div_ceil(8)];
    for (i, &b) in strand_bits.iter().enumerate() {
        if b {
            strand_bytes[i / 8] |= 1u8 << (i % 8);
        }
    }
    let mut s = Vec::new();
    write_bsc_role(&mut s, &strand_bytes, m)?;
    push_entry_buf(buf, entries, chunk_index, RefRole::Strands, mate, CODEC_BSC, &s, m)?;

    // --- 2c. ReadLen: varint/mapped read (record_count = M).
    let mut read_lens_raw = Vec::new();
    for rl in &read_lens {
        write_varint(&mut read_lens_raw, *rl);
    }
    let mut s = Vec::new();
    write_bsc_role(&mut s, &read_lens_raw, m)?;
    push_entry_buf(buf, entries, chunk_index, RefRole::ReadLen, mate, CODEC_BSC, &s, m)?;

    // --- 2d. EditCount / EditPos / EditBase.
    let enc = encode_edits(&edit_recs);
    let total_edits: u64 = edit_recs.iter().map(|e| e.len() as u64).sum();
    let mut s = Vec::new();
    write_bsc_role(&mut s, &enc.edit_count, m)?;
    push_entry_buf(buf, entries, chunk_index, RefRole::EditCount, mate, CODEC_BSC, &s, m)?;
    let mut s = Vec::new();
    write_bsc_role(&mut s, &enc.edit_pos, total_edits)?;
    push_entry_buf(buf, entries, chunk_index, RefRole::EditPos, mate, CODEC_BSC, &s, total_edits)?;
    let mut s = Vec::new();
    write_bsc_role(&mut s, &enc.edit_base, total_edits)?;
    push_entry_buf(buf, entries, chunk_index, RefRole::EditBase, mate, CODEC_BSC, &s, total_edits)?;

    // --- 3. Fallback reads → original read sequence in read order (BORROWED).
    let mut fallbacks: Vec<&[u8]> = Vec::new();
    for (i, p) in placed.iter().enumerate() {
        if matches!(p, Placed::Fallback) {
            fallbacks.push(&records[i].sequence);
        }
    }
    // --- 3a. FallbackPool inline, omitted when empty (present iff N−M > 0).
    if !fallbacks.is_empty() {
        let mut pool = Vec::new();
        let rc = fallback::write_fallback_pool(&mut pool, &fallbacks, REF_POOL_TARGET)?;
        push_entry_buf(buf, entries, chunk_index, RefRole::FallbackPool, mate, CODEC_BSC, &pool, rc)?;
    }

    // --- 4. Headers (independent columnar; never delta) + quality (fqzcomp).
    let (h_role, h_codec, h_blocks, q_blocks) = encode_ref_headers_qual(
        headers,
        quals,
        qctx_block_size,
        None, // r1_headers = None → mate-1 independent-columnar branch
        None, // r1_header_proxy = None
        RefRole::R1Headers,
    )?;
    let mut hs = Vec::new();
    crate::compression::paired::streams::write_block_stream(&mut hs, &h_blocks);
    push_entry_buf(buf, entries, chunk_index, h_role, mate, h_codec, &hs, n_chunk)?;

    let mut qs = Vec::new();
    crate::compression::paired::streams::write_block_stream(&mut qs, &q_blocks);
    push_entry_buf(buf, entries, chunk_index, RefRole::R1Qual, mate, CODEC_FQZCOMP, &qs, n_chunk)?;

    Ok(())
}

/// Encode reference-mode headers + quality directly from borrowed strings.
///
/// Byte-identical to the former `paired::encode_headers_qual_chunk(&recs, ..)`
/// call with `use_columnar_headers + use_fqzcomp + fqzcomp_subblock`, but takes
/// `headers`/`quals` slices so the caller need not reconstruct read sequences
/// into throwaway `FastqRecord`s — the columnar header codec reads only the id
/// bytes and fqzcomp reads only the quality bytes; the sequence is never used.
///
/// - Headers: `header_col` columnar blob, length-prefixed with the read count and
///   cap-split into `<=MAX_BLOCK` sub-blocks (common case = one block) — mirrors
///   `compress_one_chunk`'s `use_columnar_headers` branch.
/// - Quality: split into `qctx_block_size`-record sub-blocks, each fqzcomp'd and
///   bisected to honor the 64 MiB block cap, in parallel — mirrors the
///   `use_fqzcomp + fqzcomp_subblock` branch.
fn encode_ref_headers_qual(
    headers: &[&[u8]],
    quals: &[&[u8]],
    qctx_block_size: usize,
    // Some(R1 ids) for mate 2 → also build the delta-vs-R1 candidate and keep the
    // smaller; None for mate 1 (always columnar). Mirrors the FASTQ paired
    // `PairedHeaders` columnar-vs-delta decision — reference's R2 ids are R1 ids with
    // the mate suffix flipped (…/1 → …/2), so the delta is tiny and ~always wins.
    r1_headers: Option<&[&[u8]]>,
    // Some(R1-columnar proxy length) for mate 2 → enables the predict-and-skip of the
    // R2 columnar candidate when the delta provably wins; None disables the skip
    // (compute both and pick smaller — the prior behavior).
    r1_header_proxy: Option<usize>,
    // Base header role (R1Headers / R2HeadersIndep); becomes R2HeaderDelta when the
    // delta candidate wins for mate 2.
    base_header_role: RefRole,
) -> Result<(RefRole, u8, Vec<(u32, Vec<u8>)>, Vec<(u32, Vec<u8>)>)> {
    use rayon::prelude::*;
    let n = headers.len();

    // Headers ‖ quality run CONCURRENTLY — mirrors compress_one_chunk's
    // `rayon::join(headers || (seq || qual))`. Keeping the join matters: with the
    // default shallow encode window, serializing header-compress before quality
    // starves cores and slows the wall (measured ~10% regression without it),
    // even though total work is unchanged.
    let (h_res, q_res) = rayon::join(
        || -> Result<(RefRole, u8, Vec<(u32, Vec<u8>)>)> {
            let cap = crate::compression::paired::streams::MAX_BLOCK;
            let Some(r1) = r1_headers else {
                // mate 1: independent columnar headers (the only candidate).
                let indep =
                    crate::compression::compress_impl::columnar_blocks_capped(headers, cap)?;
                return Ok((base_header_role, CODEC_COLUMNAR, indep));
            };
            // mate 2 candidate B (CHEAP — compute FIRST): delta of R2 ids vs R1 ids,
            // BSC-compressed. On typical Illumina data R2 ids are R1 ids with the mate
            // suffix flipped (…/1 → …/2), so this delta is a few hundred bytes and
            // ~always wins, while candidate A (independent columnar) is a full BWT over
            // the chunk's headers, thrown away every chunk.
            let (ops, offs) = crate::compression::paired::header_delta::encode(r1, headers);
            let delta = crate::compression::bsc::compress_parallel_with_breakpoints(
                &ops,
                &offs,
                25 * 1024 * 1024,
                false,
            )?;
            let mut delta_s = Vec::new();
            crate::compression::paired::streams::write_block_stream(&mut delta_s, &delta);
            // When the delta is ≥2× below the R1-columnar proxy (mate-1's serialized
            // header length ≈ what R2's independent columnar would cost), it provably
            // beats candidate A, so SKIP the expensive R2 columnar BWT entirely.
            // Byte-identical to the compute-both encoder on typical data (R2 columnar
            // ≈ R1 columnar); mirrors the FASTQ-paired `PairedHeaders` skip.
            if let Some(proxy) = r1_header_proxy
                && delta_s.len().saturating_mul(2) <= proxy {
                    return Ok((RefRole::R2HeaderDelta, CODEC_BSC, delta));
                }
            // Otherwise compute candidate A and keep the smaller serialized segment
            // (tie → independent columnar, `<=`), exactly like FASTQ paired.
            let indep = crate::compression::compress_impl::columnar_blocks_capped(headers, cap)?;
            let mut indep_s = Vec::new();
            crate::compression::paired::streams::write_block_stream(&mut indep_s, &indep);
            if indep_s.len() <= delta_s.len() {
                Ok((RefRole::R2HeadersIndep, CODEC_COLUMNAR, indep))
            } else {
                Ok((RefRole::R2HeaderDelta, CODEC_BSC, delta))
            }
        },
        || -> Result<Vec<(u32, Vec<u8>)>> {
            // Quality — fqzcomp sub-blocks, byte-identical to the fqzcomp_subblock branch.
            let sub = qctx_block_size.max(1);
            let cap = crate::compression::paired::streams::MAX_BLOCK;
            let num_sub = n.div_ceil(sub);
            let nested: Vec<Vec<(u32, Vec<u8>)>> = (0..num_sub)
                .into_par_iter()
                .map(|i| {
                    let start = i * sub;
                    let end = (start + sub).min(n);
                    crate::compression::compress_impl::fqz_blocks_capped_quals(
                        &quals[start..end],
                        cap,
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(nested.into_iter().flatten().collect())
        },
    );

    let (role, codec, h_blocks) = h_res?;
    Ok((role, codec, h_blocks, q_res?))
}

/// `G` (gap-merge) cap applied to `CoverageMap::finalize` for a given read length.
/// Single source of truth so the producer and any caller resolve `g` identically.
#[inline]
pub(crate) fn gap_merge_g(read_len: usize) -> u64 {
    (read_len as u64).min(GAP_MERGE_CAP)
}

/// Finalize the coverage map into an [`IntervalMap`] (coalesce covered runs at the
/// read-length gap-merge cap). The streaming producer calls this after the last
/// chunk's coverage fold and hands the result to the thin writer to STREAM.
pub(crate) fn coverage_to_intervalmap(
    coverage: coverage::CoverageMap,
    read_len: usize,
) -> Result<IntervalMap> {
    let intervals = coverage.finalize(gap_merge_g(read_len));
    IntervalMap::from_intervals(intervals)
}

/// A `Write` adapter that counts bytes flowing through it. Replaces the old
/// `BufWriter<File>::stream_position` byte-counting (the streaming writer's `out`
/// is a `Box<dyn Write + Send>` with no `Seek`); the streamed globals are
/// byte-identical regardless of how their length is measured.
struct CountingWriter<'a> {
    inner: &'a mut dyn std::io::Write,
    written: u64,
}
impl std::io::Write for CountingWriter<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.written += n as u64;
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// Stream the 4 whole-archive globals block-by-block straight through `out`, in the
/// order PackedBacking → NBitmap → IntervalMap → ReferenceMeta, recording each as a
/// `GLOBAL_SENTINEL`/`mate=0` directory entry and advancing `*cur`. Carved verbatim
/// out of the former `finalize_reference`'s global-writing block — byte-identical
/// codecs/order/framing. PackedBacking is NEVER materialized into a buffer (it is
/// large); it streams a block at a time, its written length measured by a counting
/// `Write` wrapper. The caller writes the footer + locator after this returns.
pub(crate) fn stream_globals(
    out: &mut dyn std::io::Write,
    imap: &IntervalMap,
    refs: &[strobealign::io::fasta::RefSequence],
    cur: &mut u64,
    entries: &mut Vec<ChunkDirEntry>,
) -> Result<()> {
    // --- PackedBacking: stream block-at-a-time straight through `out`.
    let bases = {
        let mut cw = CountingWriter { inner: out, written: 0 };
        let b = backing::stream_packed_backing(&mut cw, refs, imap, backing::BACKING_BLOCK_BASES)?;
        let written = cw.written;
        push_global_entry(
            cur,
            entries,
            RefRole::PackedBacking,
            CODEC_PACKED_CONSENSUS,
            written,
            b,
        );
        b
    };

    // --- N-bitmap: stream block-at-a-time straight through `out`.
    {
        let mut cw = CountingWriter { inner: out, written: 0 };
        backing::stream_nbitmap(&mut cw, refs, imap, backing::BACKING_BLOCK_BASES)?;
        let written = cw.written;
        push_global_entry(cur, entries, RefRole::NBitmap, CODEC_BSC, written, bases);
    }

    // --- IntervalMap block-stream (small; build in a Vec then write).
    {
        let mut s = Vec::new();
        let interval_count = backing::write_intervalmap_blockstream(&mut s, imap, META_TARGET)?;
        push_entry_stream(
            out,
            cur,
            entries,
            GLOBAL_SENTINEL,
            RefRole::IntervalMap,
            0,
            CODEC_BSC,
            &s,
            interval_count,
        )?;
    }

    // --- ReferenceMeta block-stream (small; build in a Vec then write).
    {
        let mut s = Vec::new();
        let num_refs = refmeta::write_refmeta_blockstream(&mut s, refs, META_TARGET)?;
        push_entry_stream(
            out,
            cur,
            entries,
            GLOBAL_SENTINEL,
            RefRole::ReferenceMeta,
            0,
            CODEC_BSC,
            &s,
            num_refs,
        )?;
    }

    // (FallbackPool is per-(chunk,mate), written inline during encode — nothing to
    // finalize here.)
    Ok(())
}

/// Parallel sibling of [`stream_globals`] for the reference-aware MERGE
/// ([`super::merge`]). PackedBacking + N-bitmap are re-derived with rayon over
/// independent BSC blocks (the merge runs serially after the shard workers exit, so
/// all cores are idle, and the single 100M-base block would otherwise BSC-compress
/// single-threaded — the dominant merge tax). IntervalMap + ReferenceMeta are tiny
/// and stay serial (identical to `stream_globals`). The blobs are materialized (the
/// merge can afford it) then written; the resulting archive decodes identically (the
/// smaller per-block partition is internal to each global). Pushes the SAME 4 global
/// directory entries `stream_globals` would, in the same order.
pub(crate) fn stream_globals_parallel(
    out: &mut dyn std::io::Write,
    imap: &IntervalMap,
    refs: &[strobealign::io::fasta::RefSequence],
    cur: &mut u64,
    entries: &mut Vec<ChunkDirEntry>,
) -> Result<()> {
    // --- PackedBacking (parallel) ---
    let (backing_blob, bases) =
        backing::build_packed_backing_parallel(refs, imap, backing::MERGE_BACKING_BLOCK_BASES)?;
    out.write_all(&backing_blob)?;
    push_global_entry(
        cur,
        entries,
        RefRole::PackedBacking,
        CODEC_PACKED_CONSENSUS,
        backing_blob.len() as u64,
        bases,
    );

    // --- N-bitmap (parallel) ---
    let nbm_blob = backing::build_nbitmap_parallel(refs, imap, backing::MERGE_BACKING_BLOCK_BASES)?;
    out.write_all(&nbm_blob)?;
    push_global_entry(cur, entries, RefRole::NBitmap, CODEC_BSC, nbm_blob.len() as u64, bases);

    // --- IntervalMap block-stream (small; serial, identical to stream_globals). ---
    {
        let mut s = Vec::new();
        let interval_count = backing::write_intervalmap_blockstream(&mut s, imap, META_TARGET)?;
        push_entry_stream(
            out,
            cur,
            entries,
            GLOBAL_SENTINEL,
            RefRole::IntervalMap,
            0,
            CODEC_BSC,
            &s,
            interval_count,
        )?;
    }

    // --- ReferenceMeta block-stream (small; serial, identical to stream_globals). ---
    {
        let mut s = Vec::new();
        let num_refs = refmeta::write_refmeta_blockstream(&mut s, refs, META_TARGET)?;
        push_entry_stream(
            out,
            cur,
            entries,
            GLOBAL_SENTINEL,
            RefRole::ReferenceMeta,
            0,
            CODEC_BSC,
            &s,
            num_refs,
        )?;
    }
    Ok(())
}

/// One v4 block carrying `bsc::compress(raw)` with `record_count`. The decode
/// side reads it via `decode_bsc_stream`.
fn write_bsc_role(out: &mut Vec<u8>, raw: &[u8], record_count: u64) -> Result<()> {
    let comp = crate::compression::bsc::compress(raw)?;
    crate::compression::paired::streams::write_block_stream(out, &[(record_count as u32, comp)]);
    Ok(())
}

/// 25 MiB record block-stream target for the metadata globals and the
/// per-(chunk,mate) FallbackPool. Single source of truth — `stream_globals`
/// uses this module-scope const directly.
const META_TARGET: usize = 25 * 1024 * 1024;

/// Per-block target for a (chunk,mate)'s FallbackPool block-stream. Defaults to
/// the proven-safe `META_TARGET` (25 MiB) — the exact value the old GLOBAL pool
/// used and decoded fine. Do NOT raise this near `BSC_MAX_BLOCK_LEN` (64 MiB):
/// `write_record_block_stream` partitions by FRAMED size (varint len-prefix +
/// bytes), so a near-cap raw block risks tripping the `bsc::decompress` decode
/// cap. Task 1.4's ratio gate (`perchunk_pool_fragmentation_overhead_bounded` +
/// `reference_fallback_pool_is_tiny_fraction_of_realistic_archive`) confirmed
/// META_TARGET keeps the per-(chunk,mate) split's whole-archive regression ≤0.5%
/// on realistic data, so no larger/safer-ceiling target is needed.
const REF_POOL_TARGET: usize = META_TARGET;

#[cfg(test)]
mod place_against_reference_tests {
    use super::*;

    #[test]
    fn place_against_reference_extracts_edits_or_falls_back() {
        let chrom = b"ACGTACGTACGTACGTACGT".to_vec();
        let refs = vec![strobealign::io::fasta::RefSequence {
            name: "c".into(),
            sequence: chrom.clone(),
        }];
        // read = chrom[4..14] with one substitution at offset 3
        let mut read = chrom[4..14].to_vec();
        read[3] = if read[3] == b'A' { b'C' } else { b'A' };
        match place_against_reference(&refs, 0, 4, &read, false) {
            Placed::Mapped { edits, .. } => {
                assert_eq!(edits.len(), 1);
            }
            _ => panic!("expected mapped"),
        }
        // junk read that can't place within k -> fallback
        let junk = vec![b'T'; 10];
        assert!(matches!(
            place_against_reference(&refs, 0, 4, &junk, false),
            Placed::Fallback
        ));
        // out-of-range ref_id -> fallback, no panic
        assert!(matches!(
            place_against_reference(&refs, 9, 4, &read, false),
            Placed::Fallback
        ));
    }

    /// Pins the mate-2 R2-header predict-and-skip DECISION directly (mirrors the
    /// FASTQ-paired `paired_headers_picks_delta_then_columnar`): with the R1-columnar
    /// proxy threaded in, delta-friendly R2 ids skip straight to HeaderDelta, while
    /// unrelated R2 ids do NOT skip and independent columnar wins. The proxy is the
    /// exact value production passes (mate-1's serialized columnar header length).
    #[test]
    fn ref_r2_header_skip_fires_for_delta_friendly_ids_only() {
        fn r1_columnar_proxy(r1: &[&[u8]]) -> usize {
            let cap = crate::compression::paired::streams::MAX_BLOCK;
            let blocks =
                crate::compression::compress_impl::columnar_blocks_capped(r1, cap).unwrap();
            let mut s = Vec::new();
            crate::compression::paired::streams::write_block_stream(&mut s, &blocks);
            s.len()
        }
        fn decide(r1: &[Vec<u8>], r2: &[Vec<u8>], proxy: Option<usize>) -> RefRole {
            let r1r: Vec<&[u8]> = r1.iter().map(|v| v.as_slice()).collect();
            let r2r: Vec<&[u8]> = r2.iter().map(|v| v.as_slice()).collect();
            let quals: Vec<Vec<u8>> = r2.iter().map(|_| vec![b'I'; 24]).collect();
            let qr: Vec<&[u8]> = quals.iter().map(|v| v.as_slice()).collect();
            let (role, _codec, _h, _q) =
                encode_ref_headers_qual(&r2r, &qr, 64, Some(&r1r), proxy, RefRole::R2HeadersIndep)
                    .unwrap();
            role
        }

        // Delta-friendly: R2 = R1 with the mate suffix flipped → tiny delta.
        let r1: Vec<Vec<u8>> = (0..64)
            .map(|i| format!("INSTR:1:FC:1:{}:{}:{i}/1", i / 7, i % 7).into_bytes())
            .collect();
        let r2: Vec<Vec<u8>> = (0..64)
            .map(|i| format!("INSTR:1:FC:1:{}:{}:{i}/2", i / 7, i % 7).into_bytes())
            .collect();
        let r1r: Vec<&[u8]> = r1.iter().map(|v| v.as_slice()).collect();
        let proxy = r1_columnar_proxy(&r1r);

        // Skip fires → HeaderDelta, WITHOUT computing the R2 columnar candidate.
        assert_eq!(
            decide(&r1, &r2, Some(proxy)),
            RefRole::R2HeaderDelta,
            "delta-friendly R2 + R1-columnar proxy must skip straight to delta"
        );
        // No proxy → compute both; delta still wins on these ids.
        assert_eq!(
            decide(&r1, &r2, None),
            RefRole::R2HeaderDelta,
            "no proxy → compute-both still picks delta on delta-friendly ids"
        );

        // Unrelated R2 ids → large delta → the skip must NOT fire; columnar wins.
        let r2u: Vec<Vec<u8>> = (0..64)
            .map(|i| format!("LANE7:TILE{:05}:X{}:Y{}", 64 - i, i * 31, i * 17).into_bytes())
            .collect();
        assert_eq!(
            decide(&r1, &r2u, Some(proxy)),
            RefRole::R2HeadersIndep,
            "unrelated R2 ids → skip must not fire, independent columnar wins"
        );
    }
}

#[cfg(test)]
mod encode_mate_single_tests {
    use super::*;
    use crate::compression::chunk_directory::StreamRole;
    use crate::io::FastqRecord;

    fn rec(seq: &[u8]) -> FastqRecord {
        FastqRecord { id: b"id".to_vec(), sequence: seq.to_vec(), quality: Some(vec![b'I'; seq.len()]) }
    }

    #[test]
    fn encode_mate_single_emits_mate1_roles_with_expected_codecs() {
        // 3 reads: 2 mapped (one with an edit), 1 fallback → N=3, M=2, F=1.
        let placed = vec![
            Placed::Mapped { ref_id: 0, ref_pos: 10, is_rc: false, read_len: 4, edits: vec![] },
            Placed::Mapped { ref_id: 0, ref_pos: 20, is_rc: true,  read_len: 4,
                             edits: vec![edits::Edit { pos: 1, base: b'C' }] },
            Placed::Fallback,
        ];
        let records = vec![rec(b"ACGT"), rec(b"ACGT"), rec(b"TTTT")];
        let headers: Vec<&[u8]> = records.iter().map(|r| r.id.as_slice()).collect();
        let quals: Vec<&[u8]> = records.iter().map(|r| r.quality.as_deref().unwrap()).collect();

        let mut buf = Vec::new();
        let mut entries = Vec::new();
        encode_mate_single(&mut buf, &mut entries, 0, &placed, &records, &headers, &quals, 64).unwrap();

        // Every per-read entry is mate=1.
        assert!(entries.iter().all(|e| e.mate == 1), "all single-end entries are mate 1");
        // MappedFlags present with record_count = N = 3.
        let mf = entries.iter().find(|e| e.role == StreamRole::MappedFlags).unwrap();
        assert_eq!(mf.record_count, 3);
        assert_eq!(mf.codec, CODEC_BSC);
        // Positions present with record_count = M = 2, codec BSC.
        let pos = entries.iter().find(|e| e.role == StreamRole::Positions).unwrap();
        assert_eq!(pos.record_count, 2);
        assert_eq!(pos.codec, CODEC_BSC);
        // FallbackPool present (F=1), record_count = 1, codec BSC.
        let fb = entries.iter().find(|e| e.role == StreamRole::FallbackPool).unwrap();
        assert_eq!(fb.record_count, 1);
        assert_eq!(fb.codec, CODEC_BSC);
        // Headers columnar (single-end always independent, never delta).
        let h = entries.iter().find(|e| e.role == StreamRole::Headers).unwrap();
        assert_eq!(h.codec, CODEC_COLUMNAR);
        assert_eq!(h.record_count, 3);
        assert!(!entries.iter().any(|e| e.role == StreamRole::HeaderDelta), "no delta header single-end");
        // Quality fqzcomp, record_count = N.
        let q = entries.iter().find(|e| e.role == StreamRole::Qual).unwrap();
        assert_eq!(q.codec, CODEC_FQZCOMP);
        assert_eq!(q.record_count, 3);
        // EditPos == EditBase count (= total edits = 1).
        let ep = entries.iter().find(|e| e.role == StreamRole::EditPos).unwrap();
        let eb = entries.iter().find(|e| e.role == StreamRole::EditBase).unwrap();
        assert_eq!(ep.record_count, 1);
        assert_eq!(eb.record_count, 1);
    }

    #[test]
    fn encode_mate_single_omits_fallback_when_all_mapped() {
        let placed = vec![
            Placed::Mapped { ref_id: 0, ref_pos: 0, is_rc: false, read_len: 4, edits: vec![] },
        ];
        let records = vec![rec(b"ACGT")];
        let headers: Vec<&[u8]> = records.iter().map(|r| r.id.as_slice()).collect();
        let quals: Vec<&[u8]> = records.iter().map(|r| r.quality.as_deref().unwrap()).collect();
        let mut buf = Vec::new();
        let mut entries = Vec::new();
        encode_mate_single(&mut buf, &mut entries, 0, &placed, &records, &headers, &quals, 64).unwrap();
        assert!(!entries.iter().any(|e| e.role == StreamRole::FallbackPool),
                "no FallbackPool when N==M");
    }

    #[test]
    fn encode_chunk_single_returns_self_contained_chunk() {
        let placed = vec![
            Placed::Mapped { ref_id: 0, ref_pos: 0, is_rc: false, read_len: 4, edits: vec![] },
            Placed::Fallback,
        ];
        let chunk = vec![rec(b"ACGT"), rec(b"TTTT")];
        let ec = encode_chunk_single(0, &chunk, &placed, 64).unwrap();
        assert!(!ec.bytes.is_empty());
        // Every entry is mate 1 and offsets are within the chunk byte payload.
        assert!(ec.entries.iter().all(|e| e.mate == 1));
        assert!(ec.entries.iter().all(|e| e.offset + e.length <= ec.bytes.len() as u64));
        // The mate-1 Headers entry's record_count is the chunk's read count (= num_reads source).
        use crate::compression::chunk_directory::StreamRole;
        let h = ec.entries.iter().find(|e| e.role == StreamRole::Headers && e.mate == 1).unwrap();
        assert_eq!(h.record_count, 2);
    }
}

