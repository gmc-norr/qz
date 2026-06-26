//! NUMA-compress worker boundary (CLI-only orchestration lives in qz-cli). The
//! mirror of `decode_chunk_range`: the parent resolves the size-dependent plan once
//! (`resolve_plan_override`) and each node-pinned worker compresses one byte range
//! of the input (`compress_byte_range`) into a part archive. The parts merge cleanly
//! because the pinned plan makes every shard codec-homogeneous (spec §5).

use std::io::{BufReader, Read, Seek, SeekFrom};

use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};

use crate::cli::{CompressConfig, QualityCompressor};

/// Half-open byte range `[start, end)` of an input file, aligned to FASTQ record
/// boundaries by the splitter. Shared across the worker boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ByteRange {
    pub start: u64,
    pub end: u64,
}

/// The size/data-dependent compression decisions the PARENT resolves once from the
/// input head, then pins to every worker so all shards are codec-homogeneous and
/// (single-end) front-header byte-identical → the merge validates (spec §5.3).
///
/// Only genuinely size/data-derived fields live here; config-derived fields
/// (encoding_type, header codec, binning) are recomputed identically in each worker.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlanOverride {
    /// Resolved quality codec for mate 1 (single-end: the only mate).
    pub quality_mate1: QualityCompressor,
    /// Resolved quality codec for mate 2 (paired only; == mate1 for single-end).
    pub quality_mate2: QualityCompressor,
    /// Detected constant sequence length (0 = variable).
    pub const_seq_len: usize,
    /// Detected constant quality length (0 = variable / fqzcomp).
    pub const_qual_len: usize,
}

/// Resolve the canonical `PlanOverride` from the input head (single-end: one input).
/// Reads exactly the first chunk (`advanced.chunk_records`), resolving the quality
/// codec and const lengths the same way `resolve_archive_plan` does.
pub fn resolve_plan_override(args: &CompressConfig) -> Result<PlanOverride> {
    // Preserve the library's 1-or-2-inputs invariant (mod.rs ~856) — the sharded path
    // assumes it. Reject 0 or >2 here so a 3-input config can never reach the
    // single-end branch (which would silently use only input[0]).
    match args.input.len() {
        2 => return resolve_plan_override_paired(args), // Task 12
        1 => {}
        n => bail!("NUMA compress: expected 1 or 2 input files, got {n}"),
    }
    // Plan resolution needs only three facts about the input head: the record count
    // (vs the 100K fqzcomp threshold) and the constant seq/qual lengths. The old path
    // read+materialized the whole first chunk (`chunk_records`, default 2.5M records,
    // ~915 MB) into `FastqRecord`s — 3 String allocs/record — just to compute that.
    // In the NUMA driver this is the dominant SERIAL prelude: it runs once before the
    // workers fork, so it wraps the parallel per-socket core and eats the sharding win
    // (phase-timed at ~1.16 s of a ~10 s `--numa 2` compress). `scan_plan_stats` walks
    // the SAME `chunk_records` window with a cheap single-pass line-length scan (no
    // record materialization) → byte-identical plan, ~9× faster. FASTA records are
    // multiline, so keep the full reader there (defensive — the NUMA precheck already
    // excludes FASTA, so this branch is unreachable from the sharded path).
    let (count, const_seq_len, q) = if args.fasta {
        let mut reader = crate::io::FastqReader::from_path_or_stdin(&args.input[0], args.fasta)?;
        let (first_chunk, _b, _o) = super::read_chunk_records(&mut reader, args.advanced.chunk_records)?;
        let (s, q) = super::detect_const_lengths(&first_chunk, args.no_quality);
        (first_chunk.len(), s, q)
    } else {
        let st = scan_plan_stats(&args.input[0], args.advanced.chunk_records, args.no_quality)?;
        (st.count, st.const_seq_len, st.const_qual_len)
    };
    if count == 0 {
        bail!("Empty input file");
    }
    let resolved = super::resolve_quality_compressor(
        args.advanced.quality_compressor,
        count,
        args.quality_mode,
        args.no_quality,
        false,
    );
    let use_fqzcomp = resolved == QualityCompressor::Fqzcomp;
    let const_qual_len = if use_fqzcomp { 0 } else { q };
    Ok(PlanOverride {
        quality_mate1: resolved,
        quality_mate2: resolved,
        const_seq_len,
        const_qual_len,
    })
}

/// Constant seq/qual lengths and record count over the input head, computed without
/// materializing `FastqRecord`s (the cheap core of `resolve_plan_override`).
struct PlanStats {
    count: usize,
    const_seq_len: usize,
    const_qual_len: usize,
}

/// Stripped (CR/LF-trimmed) length of a raw line, matching `FastqRecord` field lengths.
fn stripped_line_len(b: &[u8]) -> usize {
    let mut n = b.len();
    while n > 0 && (b[n - 1] == b'\n' || b[n - 1] == b'\r') {
        n -= 1;
    }
    n
}

/// Single-pass scan of the first `max_records` records of a plaintext (non-FASTA)
/// FASTQ file, returning `(count, const_seq_len, const_qual_len)` byte-for-byte
/// equivalent to `read_chunk_records(max_records)` + `detect_const_lengths` — but
/// length-only (no per-record String allocation). The NUMA precondition guarantees
/// seekable plaintext FASTQ (4 lines/record); a truncated trailing record (a short
/// read mid-record) is not counted, exactly as the reader would drop it.
fn scan_plan_stats(path: &std::path::Path, max_records: usize, no_quality: bool) -> Result<PlanStats> {
    use std::io::BufRead;
    let f = std::fs::File::open(path)?;
    let mut reader = BufReader::with_capacity(1 << 20, f);
    let mut buf = Vec::with_capacity(256);
    let mut count = 0usize;
    let mut seq_len0: Option<usize> = None;
    let mut seq_const = true;
    let mut qual_len0: Option<usize> = None;
    let mut qual_const = true;
    while count < max_records {
        // Line 0: header. EOF here = clean record boundary.
        buf.clear();
        if reader.read_until(b'\n', &mut buf)? == 0 {
            break;
        }
        // Line 1: sequence.
        buf.clear();
        if reader.read_until(b'\n', &mut buf)? == 0 {
            break; // truncated record → drop
        }
        let sl = stripped_line_len(&buf);
        // Line 2: '+' separator.
        buf.clear();
        if reader.read_until(b'\n', &mut buf)? == 0 {
            break;
        }
        // Line 3: quality.
        buf.clear();
        if reader.read_until(b'\n', &mut buf)? == 0 {
            break;
        }
        let ql = stripped_line_len(&buf);
        // Commit a complete 4-line record.
        count += 1;
        match seq_len0 {
            None => seq_len0 = Some(sl),
            Some(x) if x != sl => seq_const = false,
            _ => {}
        }
        if !no_quality {
            match qual_len0 {
                None => qual_len0 = Some(ql),
                Some(x) if x != ql => qual_const = false,
                _ => {}
            }
        }
    }
    let const_seq_len = if seq_const { seq_len0.unwrap_or(0) } else { 0 };
    let const_qual_len = if no_quality || !qual_const {
        0
    } else {
        qual_len0.unwrap_or(0)
    };
    Ok(PlanStats {
        count,
        const_seq_len,
        const_qual_len,
    })
}

/// Paired plan resolver. Paired forces variable framing (const lengths are always 0), so
/// the plan needs ONLY the per-mate record count vs the 100K fqzcomp threshold. Cap the
/// scan — counting past the threshold is wasted work, and with no const-length detection
/// there is NO sampling edge at all (unlike single-end, where the const window must match
/// in-process). FASTA is excluded upstream by the NUMA precheck.
fn resolve_plan_override_paired(args: &CompressConfig) -> Result<PlanOverride> {
    const PAIRED_PLAN_SAMPLE: usize = 200_000; // 2× the 100K fqzcomp threshold — decision-exact
    let count1 = scan_plan_stats(&args.input[0], PAIRED_PLAN_SAMPLE, args.no_quality)?.count;
    let count2 = scan_plan_stats(&args.input[1], PAIRED_PLAN_SAMPLE, args.no_quality)?.count;
    if count1 == 0 || count2 == 0 {
        bail!("Empty paired input");
    }
    let resolve = |count: usize| {
        super::resolve_quality_compressor(
            args.advanced.quality_compressor,
            count,
            args.quality_mode,
            args.no_quality,
            false, // paired does not support quality modeling/delta
        )
    };
    Ok(PlanOverride {
        quality_mate1: resolve(count1),
        quality_mate2: resolve(count2),
        const_seq_len: 0,
        const_qual_len: 0,
    })
}

/// Compress ONE byte range of the input into a complete, valid part archive whose
/// chunks number from 0 (the merge renumbers — spec §8). Single-end: `ranges` has 1
/// entry. The mirror of `decode_chunk_range` (spec §5.1).
pub fn compress_byte_range(
    args: &CompressConfig,
    ranges: &[ByteRange],
    plan: &PlanOverride,
) -> Result<()> {
    // Honor -t inside this worker process (compress() does this; we bypass it).
    let _ = super::ensure_global_thread_pool(args.threads);

    // Reference shard worker: the per-mate ranges drive a bounded reference compress
    // (paired type 2 / single type 4). The reference front header is fixed metadata
    // (not data-derived), so no PlanOverride is needed and every shard's header
    // matches for the merge's front-header gate — `plan` is unused here.
    if args.reference.is_some() {
        let _ = plan;
        return super::reference::compress_reference_byte_range(args, ranges);
    }

    // 1-or-2-inputs invariant — defense in depth alongside the driver +
    // resolve_plan_override guards.
    match args.input.len() {
        2 => return compress_byte_range_paired(args, ranges, plan), // Task 12
        1 => {}
        n => bail!("compress_byte_range: expected 1 or 2 input files, got {n}"),
    }
    if ranges.len() != 1 {
        bail!("single-end compress_byte_range expects exactly 1 range, got {}", ranges.len());
    }
    let r = ranges[0];
    if r.end < r.start {
        bail!("invalid byte range [{}, {})", r.start, r.end);
    }
    let mode = super::compress_impl::ChunkedMode::Default; // ultra excluded by the precheck
    let chunk_size = args.advanced.chunk_records;
    let archive_plan = super::stream_compress::archive_plan_from_override(args, mode, plan);

    // Bounded reader: seek to start, read at most (end - start) bytes. Because
    // start/end are validated record boundaries, this yields whole records only.
    let mut file = std::fs::File::open(&args.input[0])?;
    file.seek(SeekFrom::Start(r.start))?;
    let bounded = BufReader::with_capacity(4 << 20, file.take(r.end - r.start));
    let mut reader = crate::io::FastqReader::new(bounded, args.fasta);
    let (first_chunk, _b, _o) = super::read_chunk_records(&mut reader, chunk_size)?;
    if first_chunk.is_empty() {
        bail!("empty byte range [{}, {})", r.start, r.end);
    }
    super::compress_impl::run_streaming_compress(args, mode, archive_plan, first_chunk, reader, true)
}

/// Paired byte-range compress: open both mates bounded to their record-aligned ranges,
/// read the prelude first chunk from each, build per-mate plans from the driver-PINNED plan
/// (NOT re-resolved — every shard must agree for the merge's front-header identity gate),
/// and run the shared paired streaming core. The 2-input mirror of `compress_byte_range`.
pub(super) fn compress_byte_range_paired(
    args: &CompressConfig,
    ranges: &[ByteRange],
    plan: &PlanOverride,
) -> Result<()> {
    if ranges.len() != 2 {
        bail!("paired compress_byte_range expects exactly 2 ranges, got {}", ranges.len());
    }
    let chunk_records = args.advanced.chunk_records.max(1);
    let r1 = ranges[0];
    let r2 = ranges[1];
    if r1.end < r1.start || r2.end < r2.start {
        bail!("invalid paired byte range");
    }

    // Bounded readers: seek to each range start, read at most (end - start). The splitter
    // guarantees both starts are at the SAME record index and on record boundaries, so each
    // worker sees whole, paired-aligned records.
    let mut f1 = std::fs::File::open(&args.input[0])?;
    f1.seek(SeekFrom::Start(r1.start))?;
    let mut reader1 =
        crate::io::FastqReader::new(BufReader::with_capacity(4 << 20, f1.take(r1.end - r1.start)), false);
    let mut f2 = std::fs::File::open(&args.input[1])?;
    f2.seek(SeekFrom::Start(r2.start))?;
    let mut reader2 =
        crate::io::FastqReader::new(BufReader::with_capacity(4 << 20, f2.take(r2.end - r2.start)), false);

    let first1 = super::paired::stream::read_up_to(&mut reader1, chunk_records)?;
    let first2 = super::paired::stream::read_up_to(&mut reader2, chunk_records)?;
    match (first1.is_empty(), first2.is_empty()) {
        (true, true) => bail!("empty paired byte range [{}, {}) / [{}, {})", r1.start, r1.end, r2.start, r2.end),
        (a, b) if a != b => bail!("unequal mate counts in byte range"),
        _ => {}
    }

    // Plan is PINNED by the driver — build per-mate ArchivePlan + quality codec from it.
    let plan1 = super::stream_compress::paired_mate_plan(args, plan.quality_mate1);
    let plan2 = super::stream_compress::paired_mate_plan(args, plan.quality_mate2);
    let qcodec = |q: QualityCompressor| {
        if q == QualityCompressor::Fqzcomp {
            crate::compression::codec_ids::CODEC_FQZCOMP
        } else {
            crate::compression::codec_ids::CODEC_BSC
        }
    };
    super::paired::stream::compress_paired_streaming_from_readers(
        args,
        reader1,
        reader2,
        first1,
        first2,
        plan1,
        plan2,
        qcodec(plan.quality_mate1),
        qcodec(plan.quality_mate2),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{AdvancedOptions, CompressConfig};

    fn write_fastq(path: &std::path::Path, n: usize, seq_len: usize) {
        let mut s = String::new();
        let seq = "A".repeat(seq_len);
        let q = "I".repeat(seq_len);
        for i in 0..n {
            s.push_str(&format!("@r{i}\n{seq}\n+\n{q}\n"));
        }
        std::fs::write(path, s).unwrap();
    }

    #[test]
    fn resolves_fqzcomp_when_forced_single_end() {
        use crate::cli::QualityCompressor;
        let td = tempfile::TempDir::new().unwrap();
        let input = td.path().join("in.fastq");
        write_fastq(&input, 20, 50);
        let cfg = CompressConfig {
            input: vec![input],
            output: td.path().join("o.qz"),
            working_dir: td.path().to_path_buf(),
            threads: 1,
            advanced: AdvancedOptions {
                chunk_records: 8,
                quality_compressor: QualityCompressor::Fqzcomp,
                ..AdvancedOptions::default()
            },
            ..CompressConfig::default()
        };
        let p = resolve_plan_override(&cfg).unwrap();
        assert_eq!(p.quality_mate1, QualityCompressor::Fqzcomp);
        assert_eq!(p.quality_mate2, QualityCompressor::Fqzcomp);
        // fqzcomp ⇒ const_qual_len zeroed (const-qual framing is BSC-only).
        assert_eq!(p.const_qual_len, 0);
        // const_seq_len is still detected (sequence framing is independent of quality codec).
        assert_eq!(p.const_seq_len, 50);
    }

    #[test]
    fn resolves_const_lengths_for_fixed_length_single_end() {
        let td = tempfile::TempDir::new().unwrap();
        let input = td.path().join("in.fastq");
        write_fastq(&input, 20, 50);
        let cfg = CompressConfig {
            input: vec![input],
            output: td.path().join("o.qz"),
            working_dir: td.path().to_path_buf(),
            threads: 1,
            advanced: AdvancedOptions { chunk_records: 8, ..AdvancedOptions::default() },
            ..CompressConfig::default()
        };
        let p = resolve_plan_override(&cfg).unwrap();
        assert_eq!(p.const_seq_len, 50);
        // < fqzcomp threshold (20 reads) ⇒ BSC quality ⇒ const_qual_len populated.
        assert_eq!(p.quality_mate1, QualityCompressor::Bsc);
        assert_eq!(p.const_qual_len, 50);
        assert_eq!(p.quality_mate1, p.quality_mate2);
    }

    /// The cheap length-only scanner MUST agree byte-for-byte with the old
    /// materializing path (`read_chunk_records` + `detect_const_lengths`) on every
    /// shape it can see: fixed-length, variable-length, and CRLF — across window
    /// sizes that do and don't truncate the file, with quality on and off.
    #[test]
    fn scan_plan_stats_matches_materialized_path() {
        let fixed: String = (0..30).map(|i| format!("@r{i}\nACGTACGT\n+\nIIIIIIII\n")).collect();
        let variable: String = (0..30)
            .map(|i| {
                let len = if i == 17 { 9 } else { 8 };
                format!("@r{i}\n{}\n+\n{}\n", "A".repeat(len), "I".repeat(len))
            })
            .collect();
        let crlf: String = (0..16).map(|i| format!("@r{i}\r\nACGTACGT\r\n+\r\nIIIIIIII\r\n")).collect();
        for (name, text) in [("fixed", &fixed), ("variable", &variable), ("crlf", &crlf)] {
            let td = tempfile::TempDir::new().unwrap();
            let p = td.path().join("in.fastq");
            std::fs::write(&p, text).unwrap();
            for max in [4usize, 8, 17, 18, 1000] {
                for no_q in [false, true] {
                    // Materialized reference.
                    let mut reader = crate::io::FastqReader::from_path_or_stdin(&p, false).unwrap();
                    let (chunk, _b, _o) = super::super::read_chunk_records(&mut reader, max).unwrap();
                    let (rs, rq) = super::super::detect_const_lengths(&chunk, no_q);
                    // Cheap scanner.
                    let st = scan_plan_stats(&p, max, no_q).unwrap();
                    assert_eq!(st.count, chunk.len(), "{name} max={max} count");
                    assert_eq!(st.const_seq_len, rs, "{name} max={max} const_seq_len");
                    assert_eq!(st.const_qual_len, rq, "{name} max={max} no_q={no_q} const_qual_len");
                }
            }
        }
    }
}
