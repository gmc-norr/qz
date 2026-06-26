//! Block-granular streaming SINGLE-END reference compress (spec §4 ②). Mirrors
//! `stream.rs` but drives ONE reader through a mate-agnostic per-read encode
//! (`encode_chunk_single`); archive_type 4. The paired type-2 path is untouched.

use anyhow::{Result, anyhow};
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{SyncSender, sync_channel};
use std::sync::{Arc, Mutex};

use crate::cli::CompressConfig;
use crate::compression::compress_impl::spawn_worker_pool;
use crate::io::FastqRecord;

use super::mapping::Mapper;
use super::stream::{RefPhase, RefProdTiming, RefResult, write_reference};
use super::{
    coverage_to_intervalmap, map_and_diff_single, mark_cov, reference_prelude_single,
    resolve_ref_chunk_records_cfg,
};
use super::Placed;
use super::encode::encode_chunk_single;

/// Single-end reference pool job — one chunk's per-read column/header/quality encode.
enum RefJobSingle {
    Chunk {
        chunk_index: u32,
        records: Vec<FastqRecord>,
        placed: Vec<Placed>,
        qctx_block_size: usize,
    },
}

/// The pool work fn: encode ONE single-end chunk via `encode_chunk_single`. Returns
/// the SHARED `RefResult` (single-end reuses paired's result enum — identical shape).
fn run_ref_job_single(job: RefJobSingle) -> Result<Vec<RefResult>> {
    let RefJobSingle::Chunk { chunk_index, records, placed, qctx_block_size } = job;
    let ec = encode_chunk_single(chunk_index, &records, &placed, qctx_block_size)?;
    Ok(vec![RefResult::Chunk { chunk_index, ec }])
}

/// Run the single-end reference producer on THIS (main) thread: pull records
/// prefix-then-reader, chunk at `chunk_records`, fold coverage SERIALLY in chunk
/// order (one `mark_cov` per mapped read — IntervalMap byte-identity depends on
/// it), dispatch one `RefJobSingle::Chunk` per chunk. After the last chunk,
/// finalize the `IntervalMap` and send it as the whole-archive globals.
#[allow(clippy::too_many_arguments)]
fn run_reference_producer_single(
    mapper: &Mapper,
    read_len: usize,
    qctx_block_size: usize,
    prefix: Vec<FastqRecord>,
    reader: &mut crate::io::FastqReader<crate::io::fastq::FileReader>,
    chunk_records: usize,
    job_tx: &SyncSender<RefJobSingle>,
    globals_tx: &SyncSender<RefResult>,
    abort: &AtomicBool,
) -> Result<()> {
    let contig_lens: Vec<u64> = mapper
        .references()
        .iter()
        .map(|r| r.sequence.len() as u64)
        .collect();
    let mut coverage = super::coverage::CoverageMap::new(&contig_lens);

    let mut prefix = prefix.into_iter();
    let mut chunk_index: u32 = 0u32;
    let mut buf: Vec<FastqRecord> = Vec::with_capacity(chunk_records);

    let mut next_rec = |reader: &mut crate::io::FastqReader<crate::io::fastq::FileReader>|
     -> Result<Option<FastqRecord>> {
        if let Some(r) = prefix.next() {
            return Ok(Some(r));
        }
        reader.next()
    };

    let mut timing = RefProdTiming::new();

    let emit = |buf: &mut Vec<FastqRecord>,
                ci: &mut u32,
                cov: &mut super::coverage::CoverageMap,
                tm: &mut RefProdTiming|
     -> Result<()> {
        let records = std::mem::take(buf);
        let t = tm.tick();
        let placed = map_and_diff_single(mapper, &records)?;
        tm.record(RefPhase::Map, t);
        // Serial in-order coverage fold (IntervalMap byte-identity depends on it).
        let t = tm.tick();
        for p in &placed {
            mark_cov(cov, p);
        }
        tm.record(RefPhase::Cov, t);
        let t = tm.tick();
        job_tx
            .send(RefJobSingle::Chunk { chunk_index: *ci, records, placed, qctx_block_size })
            .map_err(|_| anyhow!("reference workers exited"))?;
        tm.record(RefPhase::Send, t);
        *ci += 1;
        Ok(())
    };

    while let Some(r) = next_rec(reader)? {
        if abort.load(Ordering::Relaxed) {
            return Ok(());
        }
        buf.push(r);
        if buf.len() == chunk_records {
            emit(&mut buf, &mut chunk_index, &mut coverage, &mut timing)?;
        }
    }
    if !buf.is_empty() {
        emit(&mut buf, &mut chunk_index, &mut coverage, &mut timing)?;
    }

    let t = timing.tick();
    let imap = coverage_to_intervalmap(coverage, read_len)?;
    timing.record(RefPhase::Finalize, t);
    timing.report("single", chunk_index);
    globals_tx
        .send(RefResult::Globals { imap })
        .map_err(|_| anyhow!("reference writer exited"))?;
    Ok(())
}

/// Streaming single-end reference compress entry point. Mirrors
/// `stream::compress_reference_streaming`'s scope/shutdown shape exactly, but uses
/// the single-end prelude/producer/job fn and writes archive_type 4.
pub(crate) fn compress_reference_single_streaming(
    a: &CompressConfig,
    range: Option<&crate::compression::ByteRange>,
) -> Result<()> {
    let (mapper, read_len, mut reader, prefix, qctx_block_size) =
        reference_prelude_single(a, range)?;
    let chunk_records = resolve_ref_chunk_records_cfg(a);

    let file = std::fs::File::create(&a.output)
        .map_err(|e| anyhow!("create reference archive {}: {e}", a.output.display()))?;
    let mut out: Box<dyn Write + Send> = Box::new(std::io::BufWriter::new(file));
    let header_size = crate::compression::write_v5_archive_header(
        &mut out,
        0, // encoding_type (Raw)
        crate::compression::columnar::QualityBinning::None,
        crate::cli::QualityCompressor::Bsc,
        crate::cli::HeaderCompressor::Columnar,
        false, // fasta
        0,     // const_seq_len
        0,     // const_qual_len
        4,     // archive_type = single-end reference
    )? as u64;

    let n_workers = std::env::var("QZ_COMPRESS_WORKERS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n >= 1)
        .unwrap_or_else(|| {
            a.reference
                .as_ref()
                .map(|r| r.reference_window)
                .unwrap_or(4)
                .max(1)
        });

    let (job_tx, job_rx) = sync_channel::<RefJobSingle>(n_workers);
    let (res_tx, res_rx) = sync_channel::<RefResult>(n_workers.max(1) * 2);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let abort = Arc::new(AtomicBool::new(false));

    std::thread::scope(|scope| -> Result<()> {
        spawn_worker_pool(
            scope,
            n_workers,
            &job_rx,
            &res_tx,
            &abort,
            run_ref_job_single,
            RefResult::Err,
        );
        drop(job_rx);

        let globals_tx = res_tx.clone();
        let writer = {
            let abort = abort.clone();
            let mapper = &mapper;
            scope.spawn(move || {
                // Single-end reuses the SHARED writer (`stream::write_reference`) with the
                // mate-1-only directory validator. The paired type-2 writer is unchanged.
                write_reference(
                    out,
                    header_size,
                    mapper,
                    res_rx,
                    &abort,
                    super::validate_reference_directory_single,
                    1, // single-end reference: mate 1 only
                )
            })
        };
        drop(res_tx);

        let prod = run_reference_producer_single(
            &mapper,
            read_len,
            qctx_block_size,
            prefix,
            &mut reader,
            chunk_records,
            &job_tx,
            &globals_tx,
            &abort,
        );
        if prod.is_err() {
            abort.store(true, Ordering::Relaxed);
        }
        drop(job_tx);
        drop(globals_tx);

        let w = writer
            .join()
            .map_err(|_| anyhow!("reference writer panicked"))?;
        w.and(prod)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::ReferenceOptions;
    use std::path::PathBuf;

    fn make_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_add(0x9E3779B97F4A7C15);
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            v.push(b"ACGT"[((x >> 33) & 3) as usize]);
        }
        v
    }

    fn write_single_fixture(dir: &std::path::Path) -> (PathBuf, PathBuf) {
        let refseq = make_seq(600, 42);
        let rf = dir.join("ref.fa");
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        std::fs::write(&rf, &s).unwrap();
        let mut fq = String::new();
        for (i, st) in [0usize, 100, 240, 360, 400].iter().enumerate() {
            let r = &refseq[*st..*st + 120];
            let q = "I".repeat(120);
            fq.push_str(&format!("@read_{i}\n{}\n+\n{q}\n", std::str::from_utf8(r).unwrap()));
        }
        let r1p = dir.join("r1.fastq");
        std::fs::write(&r1p, &fq).unwrap();
        (rf, r1p)
    }

    fn cfg(input: PathBuf, refp: PathBuf, out: PathBuf, tmp: PathBuf) -> CompressConfig {
        CompressConfig {
            input: vec![input],
            output: out,
            working_dir: tmp,
            threads: 1,
            force: true,
            reference: Some(ReferenceOptions {
                reference: refp,
                reference_index: None,
                reference_fast: false,
                reference_window: 2,
            }),
            ..CompressConfig::default()
        }
    }

    #[test]
    fn prelude_single_opens_one_reader_and_profiles_read_len() {
        let d = tempfile::tempdir().unwrap();
        let (rf, r1p) = write_single_fixture(d.path());
        let c = cfg(r1p, rf, d.path().join("o.qz"), d.path().to_path_buf());
        let (mapper, read_len, mut reader, prefix, qbs) =
            reference_prelude_single(&c, None).unwrap();
        assert_eq!(read_len, 120);
        assert!(mapper.num_references() >= 1);
        assert_eq!(prefix.len(), 5, "all 5 reads buffered as the prefix");
        assert!(qbs >= 1);
        // The reader is drained: the prefix consumed every record.
        assert!(reader.next().unwrap().is_none());
    }

    #[test]
    fn map_and_diff_single_is_order_preserving_and_deterministic() {
        use crate::io::FastqRecord;
        use super::super::map_and_diff_single;
        let d = tempfile::tempdir().unwrap();
        let (rf, _r1p) = write_single_fixture(d.path());
        let refseq = make_seq(600, 42);
        let c = CompressConfig {
            input: vec![d.path().join("r1.fastq")],
            output: d.path().join("o.qz"),
            working_dir: d.path().to_path_buf(),
            threads: 1,
            force: true,
            reference: Some(ReferenceOptions {
                reference: rf,
                reference_index: None,
                reference_fast: false,
                reference_window: 2,
            }),
            ..CompressConfig::default()
        };
        let (mapper, _rl, _rd, _pf, _qbs) = reference_prelude_single(&c, None).unwrap();
        // 3 exact-substring mapped reads + 1 junk fallback, in a fixed order.
        let recs: Vec<FastqRecord> = [
            refseq[0..120].to_vec(),
            refseq[240..360].to_vec(),
            vec![b'T'; 120], // junk → Fallback
            refseq[100..220].to_vec(),
        ]
        .into_iter()
        .map(|seq| FastqRecord {
            id: b"r".to_vec(),
            sequence: seq,
            quality: Some(vec![b'I'; 120]),
        })
        .collect();
        let a = map_and_diff_single(&mapper, &recs).unwrap();
        let b = map_and_diff_single(&mapper, &recs).unwrap();
        assert_eq!(a.len(), 4);
        assert_eq!(a, b, "deterministic");
        assert!(matches!(a[2], super::super::Placed::Fallback), "junk read 3 falls back");
        assert!(matches!(a[0], super::super::Placed::Mapped { .. }), "read 0 maps");
    }

    #[test]
    fn producer_single_chunks_in_order_and_sends_globals() {
        
        let d = tempfile::tempdir().unwrap();
        let (rf, _r1p) = write_single_fixture(d.path());
        let c = cfg(d.path().join("r1.fastq"), rf, d.path().join("o.qz"), d.path().to_path_buf());
        let (mapper, read_len, mut reader, prefix, qbs) = reference_prelude_single(&c, None).unwrap();

        // Channels (unbounded-enough for the 5-read fixture at chunk_records = 2).
        let (job_tx, job_rx) = std::sync::mpsc::sync_channel::<RefJobSingle>(8);
        let (res_tx, res_rx) = std::sync::mpsc::sync_channel::<RefResult>(8);
        let abort = AtomicBool::new(false);
        // Drain jobs on a thread so the bounded producer doesn't block.
        let jh = std::thread::spawn(move || {
            let mut chunks: Vec<u32> = Vec::new();
            for j in job_rx {
                let RefJobSingle::Chunk { chunk_index, records, placed, .. } = j;
                assert_eq!(records.len(), placed.len());
                chunks.push(chunk_index);
            }
            chunks
        });
        run_reference_producer_single(
            &mapper, read_len, qbs, prefix, &mut reader, 2, &job_tx, &res_tx, &abort,
        )
        .unwrap();
        drop(job_tx);
        let chunks = jh.join().unwrap();
        assert_eq!(chunks, vec![0, 1, 2], "5 reads / chunk_records 2 → chunks 0,1,2 in order");
        // Drop the producer's result-channel sender so the `res_rx` drain below
        // sees the disconnect after the single buffered Globals message (the
        // producer only BORROWED `&res_tx`, so this scope still owns it).
        drop(res_tx);
        // Globals were sent exactly once.
        let mut got_globals = false;
        for r in res_rx {
            if let RefResult::Globals { .. } = r {
                assert!(!got_globals, "globals sent once");
                got_globals = true;
            }
        }
        assert!(got_globals, "producer sent Globals after last chunk");
        let _ = (read_len, &mapper); // silence unused on some toolchains
    }

    #[test]
    fn compress_single_writes_type4_archive() {
        let d = tempfile::tempdir().unwrap();
        let (rf, r1p) = write_single_fixture(d.path());
        let out = d.path().join("o.qz");
        let c = cfg(r1p, rf, out.clone(), d.path().to_path_buf());
        // Direct call (the mod.rs dispatch is wired in Task 2.8).
        super::compress_reference_single_streaming(&c, None).expect("single-end reference compress");
        let b = std::fs::read(&out).unwrap();
        assert!(b.len() > 18, "archive too small: {}", b.len());
        assert_eq!(&b[0..2], b"QZ", "magic");
        assert_eq!(b[2], 5, "v5 chunk-major version");
        assert_eq!(b[3], 4, "single-end reference archive_type = 4");
    }
}
