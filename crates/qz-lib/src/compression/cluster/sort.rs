//! External bucket sort for the order-drop cluster mode.
//!
//! `cluster_sort_single` (single-end) and `cluster_sort_paired` (paired-end) are the
//! production entry points: they read FASTQ (plain or gzip), partition reads into 4096
//! on-disk buckets keyed by their strand-symmetric open-syncmer hash, then
//! stream the sorted records one at a time via a callback.  `cluster_sort_paired`
//! clusters reads **as pairs** — both mates ride one bucket entry keyed by R1's
//! syncmer (each mate canonicalized independently), so pairing is preserved with no
//! order stream.  **Bounded RAM, constant in read count**: only one bucket is
//! resident at a time and the callback drives downstream work without accumulating
//! all records.  The true high-water mark is the splitter-sampling prelude
//! ([`compute_splitters`]), which transiently holds up to `SAMPLE_READS` sampled
//! sequences (≈1.2 GB at 150 bp) — still a *constant* cap independent of the total
//! read count, not an O(reads) peak.

use anyhow::{anyhow, bail, ensure, Context, Result};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::time::Instant;

use crate::compression::cluster::checksum::{hash_pair, hash_record};
use crate::compression::dna_utils::{
    canonicalize_open_syncmer_strand, reverse_complement_canonical,
    strand_symmetric_open_syncmer_hash_pos,
};

const NBUCKETS: usize = 4096;
const SAMPLE_READS: usize = 8_000_000;
const BATCH: usize = 1_000_000;
/// Paired Pass 1 batch size, in *pairs*. Each pair stages ~2× a single record's bytes
/// (two mates) into the batch arena, so half of [`BATCH`] keeps the per-batch arena +
/// computed footprint on par with single-end's 1M-record batch.
const PAIR_BATCH: usize = BATCH / 2;

// ── Record struct ─────────────────────────────────────────────────────────────

/// One record in cluster-emission order.
#[derive(Clone)]
pub(crate) struct ClusteredRecord {
    /// Raw header line (no trailing newline).
    pub header: Vec<u8>,
    /// Canonicalized sequence (may be RC of original if `flip == true`).
    pub canon_seq: Vec<u8>,
    /// Raw quality string (no trailing newline).
    pub qual: Vec<u8>,
    /// `true` when the original sequence was reverse-complemented to form `canon_seq`.
    pub flip: bool,
    /// Start position of the anchor syncmer within `canon_seq`.
    pub apos: u16,
    /// Strand-symmetric open-syncmer hash (the bucket key).
    pub key: u64,
    /// 0-based read index in the original input (for deterministic tie-breaking).
    pub original_index: u64,
}

impl ClusteredRecord {
    /// Reconstruct the **original** (input) sequence: un-flip when RC-canonicalized.
    /// Used by the sort tests to assert set-losslessness; the decode path does its
    /// own un-flipping from the StrandBits stream, so this is test-only today.
    #[allow(dead_code)]
    pub fn original_seq(&self) -> Vec<u8> {
        if self.flip {
            reverse_complement_canonical(&self.canon_seq)
        } else {
            self.canon_seq.clone()
        }
    }
}

// ── I/O helpers ───────────────────────────────────────────────────────────────

fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut f = File::open(path).with_context(|| format!("open {:?}", path))?;
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic).unwrap_or(0);
    use std::io::Seek;
    f.seek(std::io::SeekFrom::Start(0))?;
    if n == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        Ok(Box::new(BufReader::with_capacity(
            1 << 22,
            MultiGzDecoder::new(BufReader::with_capacity(1 << 22, f)),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 22, f)))
    }
}

/// One parsed 4-line FASTQ record.  The four lines are concatenated into `buf`, and
/// the fields expose precise byte ranges so callers can slice directly:
/// header = `buf[0..h_end]`, seq = `buf[ss..se]`, qual = `buf[qs..qe]`.
struct Record {
    /// Raw bytes of the four lines concatenated into one buffer.
    buf: Vec<u8>,
    /// Exclusive end of the header field (before the '\n' if present).
    h_end: usize,
    /// Start of seq.
    ss: usize,
    /// End of seq (exclusive, stripping '\n').
    se: usize,
    /// Start of qual.
    qs: usize,
    /// End of qual (exclusive, stripping '\n').
    qe: usize,
}

impl Record {
    #[inline]
    fn header(&self) -> &[u8] {
        &self.buf[..self.h_end]
    }
    #[inline]
    fn seq(&self) -> &[u8] {
        &self.buf[self.ss..self.se]
    }
    #[inline]
    fn qual(&self) -> &[u8] {
        &self.buf[self.qs..self.qe]
    }
}

/// Read the next 4-line record into a fresh `Record`, or return `None` at EOF.
fn read_record(r: &mut dyn BufRead) -> Result<Option<Record>> {
    let mut buf = Vec::with_capacity(512);

    // Line 1: header
    if r.read_until(b'\n', &mut buf)? == 0 {
        return Ok(None);
    }
    let h_end = if buf.last() == Some(&b'\n') { buf.len() - 1 } else { buf.len() };

    // Line 2: seq
    let ss = buf.len();
    if r.read_until(b'\n', &mut buf)? == 0 {
        return Err(anyhow!("truncated FASTQ: missing seq line"));
    }
    let se = if buf.last() == Some(&b'\n') { buf.len() - 1 } else { buf.len() };

    // Line 3: '+' line (discard content)
    if r.read_until(b'\n', &mut buf)? == 0 {
        return Err(anyhow!("truncated FASTQ: missing '+' line"));
    }

    // Line 4: qual
    let qs = buf.len();
    if r.read_until(b'\n', &mut buf)? == 0 {
        return Err(anyhow!("truncated FASTQ: missing qual line"));
    }
    let qe = if buf.last() == Some(&b'\n') { buf.len() - 1 } else { buf.len() };

    // Mirror the standard reader (`io/fastq.rs`): sequence and quality lengths MUST match.
    // Cluster decode reconstructs each record's sequence boundary purely from its quality
    // length (`decode.rs`), so a mismatch would silently corrupt the output (wrong, mis-
    // boundaried sequences) or render the archive undecodable. Reject it at ingestion.
    let seq_len = se - ss;
    let qual_len = qe - qs;
    if seq_len != qual_len {
        return Err(anyhow!(
            "Invalid FASTQ: sequence length ({seq_len}) != quality length ({qual_len}) for read {}",
            String::from_utf8_lossy(&buf[..h_end])
        ));
    }

    Ok(Some(Record { buf, h_end, ss, se, qs, qe }))
}

// ── Pass 0: splitter sampling ─────────────────────────────────────────────────

/// Sample up to `SAMPLE_READS` reads, compute strand-symmetric syncmer keys, and pick
/// 4095 balanced quantile splitters that partition the key space into `NBUCKETS`
/// roughly equal buckets.
///
/// RAM note: samples up to `SAMPLE_READS` reads, but hashes them in bounded batches and
/// retains only the `u64` keys, so the high-water mark is one batch of sequences
/// (≈`SPLITTER_BATCH`×len ≈ 77 MB) plus the key vector (≈64 MB) — not the ≈1.2 GB it
/// would take to hold all sampled sequences at once. Both caps are constant in the total
/// read count. Keys are identical to a single-shot hash (same sample reads → same keys →
/// same splitters), so the archive is byte-for-byte unchanged.
fn compute_splitters(path: &Path) -> Result<Vec<u64>> {
    /// Reads hashed per batch — bounds the transient sequence buffer in `compute_splitters`.
    const SPLITTER_BATCH: usize = 1 << 19; // 512K reads ≈ 77 MB at 150 bp
    let mut reader = open_maybe_gz(path)?;
    let mut keys: Vec<u64> = Vec::with_capacity(SAMPLE_READS.min(1 << 20));
    let mut batch: Vec<Vec<u8>> = Vec::with_capacity(SPLITTER_BATCH);
    while keys.len() < SAMPLE_READS {
        let want = (SAMPLE_READS - keys.len()).min(SPLITTER_BATCH);
        batch.clear();
        while batch.len() < want {
            match read_record(&mut *reader)? {
                Some(rec) => batch.push(rec.seq().to_vec()),
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }
        keys.par_extend(
            batch
                .par_iter()
                .map(|s| strand_symmetric_open_syncmer_hash_pos(s).0),
        );
    }
    keys.par_sort_unstable();
    if keys.is_empty() {
        return Ok(Vec::new()); // no reads -> no splitters; every (absent) read trivially buckets to 0
    }
    let n = keys.len().max(1);
    Ok((1..NBUCKETS).map(|i| keys[(i * n / NBUCKETS).min(n - 1)]).collect())
}

// ── Bucket entry layout ───────────────────────────────────────────────────────
// [key:u64 LE][apos:u16 LE][flip:u8][original_index:u64 LE]
// [hlen:u32 LE][header bytes]
// [slen:u32 LE][canon seq bytes]
// [qlen:u32 LE][qual bytes]

/// Serialize one bucket entry directly from borrowed field slices — byte-identical
/// to the old `write_entry(&ClusteredRecord)` but without materializing an owned
/// `ClusteredRecord` (3 heap allocs/record) just to write and drop it. `header`/`qual`
/// borrow the batch arena; `canon_seq` borrows the parallel-stage output.
#[allow(clippy::too_many_arguments)]
fn write_entry_parts(
    w: &mut BufWriter<File>,
    key: u64,
    apos: u16,
    flip: bool,
    original_index: u64,
    header: &[u8],
    canon_seq: &[u8],
    qual: &[u8],
) -> Result<()> {
    w.write_all(&key.to_le_bytes())?;
    w.write_all(&apos.to_le_bytes())?;
    w.write_all(&[flip as u8])?;
    w.write_all(&original_index.to_le_bytes())?;
    w.write_all(&(header.len() as u32).to_le_bytes())?;
    w.write_all(header)?;
    w.write_all(&(canon_seq.len() as u32).to_le_bytes())?;
    w.write_all(canon_seq)?;
    w.write_all(&(qual.len() as u32).to_le_bytes())?;
    w.write_all(qual)?;
    Ok(())
}

/// Parse all entries from a raw bucket blob.
fn parse_bucket(data: &[u8]) -> Result<Vec<ClusteredRecord>> {
    let mut ents = Vec::new();
    let mut p = 0usize;

    // Minimum fixed header: key(8)+apos(2)+flip(1)+original_index(8)+hlen(4) = 23
    while p + 23 <= data.len() {
        let key = u64::from_le_bytes(
            data[p..p + 8].try_into().map_err(|_| anyhow!("bucket parse: key slice"))?,
        );
        let apos = u16::from_le_bytes(
            data[p + 8..p + 10].try_into().map_err(|_| anyhow!("bucket parse: apos slice"))?,
        );
        let flip = data[p + 10] != 0;
        let original_index = u64::from_le_bytes(
            data[p + 11..p + 19]
                .try_into()
                .map_err(|_| anyhow!("bucket parse: original_index slice"))?,
        );
        p += 19;

        let hlen = u32::from_le_bytes(
            data[p..p + 4].try_into().map_err(|_| anyhow!("bucket parse: hlen slice"))?,
        ) as usize;
        p += 4;
        ensure!(p + hlen <= data.len(), "bucket parse: header truncated");
        let header = data[p..p + hlen].to_vec();
        p += hlen;

        let slen = u32::from_le_bytes(
            data[p..p + 4].try_into().map_err(|_| anyhow!("bucket parse: slen slice"))?,
        ) as usize;
        p += 4;
        ensure!(p + slen <= data.len(), "bucket parse: seq truncated");
        let canon_seq = data[p..p + slen].to_vec();
        p += slen;

        let qlen = u32::from_le_bytes(
            data[p..p + 4].try_into().map_err(|_| anyhow!("bucket parse: qlen slice"))?,
        ) as usize;
        p += 4;
        ensure!(p + qlen <= data.len(), "bucket parse: qual truncated");
        let qual = data[p..p + qlen].to_vec();
        p += qlen;

        ents.push(ClusteredRecord { header, canon_seq, qual, flip, apos, key, original_index });
    }
    Ok(ents)
}

// ── Public entry point ────────────────────────────────────────────────────────

/// Cluster-sort all records from `input` (plain or gzip FASTQ) using on-disk
/// buckets under `work/`, then stream each record in cluster order by calling
/// `on_record`.
///
/// Returns `(verify_checksum, nreads)` where `verify_checksum` is the
/// order-independent XOR-wrapping sum of `hash_record(header, original_seq, qual)`
/// over all reads.
///
/// **Bounded RAM**: at most one bucket's raw bytes plus one bucket's parsed
/// `Vec<ClusteredRecord>` are resident at a time.  The callback controls what
/// the caller retains.
pub(crate) fn cluster_sort_single(
    input: &Path,
    work: &Path,
    mut on_record: impl FnMut(&ClusteredRecord) -> Result<()>,
) -> Result<(u128, u64)> {
    std::fs::create_dir_all(work)
        .with_context(|| format!("create work dir {:?}", work))?;

    // All 4096 bucket scratch files live in a dedicated subdir whose `Drop`
    // recursively removes them — on success, on any `?` early-return, AND on a
    // panic (the workspace builds with panic=unwind). This keeps scratch usage
    // bounded and self-cleaning (a WGS run stages ≈ input-size of buckets here),
    // and gives every invocation an isolated subdir so re-running into the same
    // `work` can never mix stale data.
    let bucket_dir = tempfile::TempDir::new_in(work)
        .with_context(|| format!("create cluster bucket scratch dir in {:?}", work))?;
    let bdir = bucket_dir.path();

    // QZ_CLUSTER_TIMING=1 emits a Pass-0/1/2 wall breakdown to stderr at the end.
    // Measured unconditionally (only a handful of `Instant`s per multi-million-record
    // batch — nil overhead); the env var only gates the final eprintln.
    let timing = std::env::var_os("QZ_CLUSTER_TIMING").is_some();
    let mut ns_fill: u128 = 0; // Pass 1: read_record + arena fill
    let mut ns_par: u128 = 0; //  Pass 1: par_iter canon + key
    let mut ns_fold: u128 = 0; // Pass 1: serial checksum-fold + bucket write
    let mut ns_p2_io: u128 = 0; // Pass 2: bucket file read
    let mut ns_p2_parse: u128 = 0; // Pass 2: parse_bucket
    let mut ns_p2_sort: u128 = 0; // Pass 2: within-bucket sort
    let mut ns_p2_emit: u128 = 0; // Pass 2: on_record callback (downstream encode)

    // ── Pass 0: splitters ─────────────────────────────────────────────────────
    let t_p0 = Instant::now();
    let splitters = compute_splitters(input)?;
    let ns_pass0 = t_p0.elapsed().as_nanos();

    // ── Open all bucket writers ───────────────────────────────────────────────
    let mut writers: Vec<BufWriter<File>> = (0..NBUCKETS)
        .map(|i| {
            let path = bdir.join(format!("b{i:04}"));
            File::create(&path)
                .with_context(|| format!("create bucket {:?}", path))
                .map(|f| BufWriter::with_capacity(1 << 18, f))
        })
        .collect::<Result<Vec<_>>>()?;

    // ── Pass 1: parse → canon → key → bucket ─────────────────────────────────
    // Batched: store (header, seq, qual) byte-ranges in a flat arena per batch,
    // compute (key, apos, flip, canon) in parallel over the batch, then fold the
    // checksum and write entries sequentially.
    let mut reader = open_maybe_gz(input)?;
    let mut nreads: u64 = 0;
    let mut checksum: u128 = 0;

    // Arena stores consecutive records: [header_bytes][seq_bytes][qual_bytes].
    // meta: (original_index, h_start, h_end, s_start, s_end, q_start, q_end)
    let mut arena: Vec<u8> = Vec::with_capacity(BATCH * 512);
    let mut meta: Vec<(u64, usize, usize, usize, usize, usize, usize)> =
        Vec::with_capacity(BATCH);

    loop {
        arena.clear();
        meta.clear();

        // Fill a batch
        let t_fill = Instant::now();
        while meta.len() < BATCH {
            match read_record(&mut *reader)? {
                Some(rec) => {
                    let idx = nreads;
                    nreads += 1;

                    let h_start = arena.len();
                    arena.extend_from_slice(rec.header());
                    let h_end = arena.len();

                    let s_start = arena.len();
                    arena.extend_from_slice(rec.seq());
                    let s_end = arena.len();

                    let q_start = arena.len();
                    arena.extend_from_slice(rec.qual());
                    let q_end = arena.len();

                    meta.push((idx, h_start, h_end, s_start, s_end, q_start, q_end));
                }
                None => break,
            }
        }
        ns_fill += t_fill.elapsed().as_nanos();

        if meta.is_empty() {
            break;
        }

        // Parallel: compute (key, apos, flip, canon, record-hash) for each record.
        // The verify checksum (`hash_record` over the ORIGINAL header/seq/qual) used to
        // be folded serially below; it's an associative+commutative `wrapping_add` of
        // per-record hashes, so computing each record's hash HERE (across all cores)
        // and summing the precomputed values serially leaves the checksum bit-identical
        // while moving the dominant per-record cost off the serial path.
        let t_par = Instant::now();
        let computed: Vec<(u64, u16, bool, Vec<u8>, u128)> = meta
            .par_iter()
            .map(|&(_, h0, h1, s0, s1, q0, q1)| {
                let header = &arena[h0..h1];
                let seq = &arena[s0..s1]; // original orientation
                let qual = &arena[q0..q1];
                let (canon, flip, key, apos) = canonicalize_open_syncmer_strand(seq);
                let rec_hash = hash_record(header, seq, qual);
                (key, apos, flip, canon, rec_hash)
            })
            .collect();
        ns_par += t_par.elapsed().as_nanos();

        // Sequential: fold the precomputed checksum (same add-order as before → identical
        // value) + write bucket entries straight from the arena/canon slices (no per-record
        // ClusteredRecord allocation). The bucket writers are `&mut` and shared, so the
        // write itself stays serial; only this cheap route+memcpy remains on the serial path.
        let t_fold = Instant::now();
        for ((key, apos, flip, canon, rec_hash), &(original_index, h0, h1, _s0, _s1, q0, q1)) in
            computed.iter().zip(meta.iter())
        {
            let header = &arena[h0..h1];
            let qual = &arena[q0..q1];

            checksum = checksum.wrapping_add(*rec_hash);

            let bucket = splitters.partition_point(|&s| s <= *key);
            write_entry_parts(
                &mut writers[bucket],
                *key,
                *apos,
                *flip,
                original_index,
                header,
                canon,
                qual,
            )?;
        }
        ns_fold += t_fold.elapsed().as_nanos();
    }

    // Flush and drop all writers before Pass 2 reads the files.
    for w in &mut writers {
        w.flush()?;
    }
    drop(writers);

    // ── Pass 2: within-bucket sort → stream ──────────────────────────────────
    let mut emitted: u64 = 0;
    for i in 0..NBUCKETS {
        let path = bdir.join(format!("b{i:04}"));
        let t_io = Instant::now();
        let data = std::fs::read(&path)
            .with_context(|| format!("read bucket {:?}", path))?;
        ns_p2_io += t_io.elapsed().as_nanos();

        if data.is_empty() {
            continue;
        }

        let t_parse = Instant::now();
        let mut ents = parse_bucket(&data)?;
        ns_p2_parse += t_parse.elapsed().as_nanos();

        // Total-order comparator: key → anchor-suffix → full canon → original_index.
        let t_sort = Instant::now();
        ents.sort_unstable_by(|a, b| {
            a.key.cmp(&b.key).then_with(|| {
                let sa_start = (a.apos as usize).min(a.canon_seq.len());
                let sb_start = (b.apos as usize).min(b.canon_seq.len());
                let sa = &a.canon_seq[sa_start..];
                let sb = &b.canon_seq[sb_start..];
                sa.cmp(sb)
                    .then_with(|| a.canon_seq.cmp(&b.canon_seq))
                    .then_with(|| a.original_index.cmp(&b.original_index))
            })
        });
        ns_p2_sort += t_sort.elapsed().as_nanos();

        let t_emit = Instant::now();
        for rec in &ents {
            on_record(rec)?;
            emitted += 1;
        }
        ns_p2_emit += t_emit.elapsed().as_nanos();

        // Drop ents + data before loading the next bucket (bounded RAM).
        drop(ents);
        drop(data);
    }

    ensure!(emitted == nreads, "emitted {emitted} != nreads {nreads}: bucket accounting error");

    if timing {
        let pass1 = ns_fill + ns_par + ns_fold;
        let pass2 = ns_p2_io + ns_p2_parse + ns_p2_sort + ns_p2_emit;
        let ms = |ns: u128| ns as f64 / 1e6;
        let pct = |ns: u128, den: u128| if den == 0 { 0.0 } else { 100.0 * ns as f64 / den as f64 };
        eprintln!("[QZ_CLUSTER_TIMING] cluster_sort_single  reads={nreads}");
        eprintln!(
            "  Pass0 splitters        {:>8.0} ms",
            ms(ns_pass0)
        );
        eprintln!(
            "  Pass1 TOTAL            {:>8.0} ms",
            ms(pass1)
        );
        eprintln!(
            "    fill (read+arena)    {:>8.0} ms  {:>5.1}% of Pass1   <-- mmap-parse target",
            ms(ns_fill),
            pct(ns_fill, pass1)
        );
        eprintln!(
            "    par_iter canon+key   {:>8.0} ms  {:>5.1}% of Pass1",
            ms(ns_par),
            pct(ns_par, pass1)
        );
        eprintln!(
            "    fold checksum+write  {:>8.0} ms  {:>5.1}% of Pass1",
            ms(ns_fold),
            pct(ns_fold, pass1)
        );
        eprintln!(
            "  Pass2 TOTAL            {:>8.0} ms  (incl. on_record encode)",
            ms(pass2)
        );
        eprintln!(
            "    bucket read          {:>8.0} ms  {:>5.1}% of Pass2",
            ms(ns_p2_io),
            pct(ns_p2_io, pass2)
        );
        eprintln!(
            "    parse_bucket         {:>8.0} ms  {:>5.1}% of Pass2",
            ms(ns_p2_parse),
            pct(ns_p2_parse, pass2)
        );
        eprintln!(
            "    within-bucket sort   {:>8.0} ms  {:>5.1}% of Pass2",
            ms(ns_p2_sort),
            pct(ns_p2_sort, pass2)
        );
        eprintln!(
            "    on_record (encode)   {:>8.0} ms  {:>5.1}% of Pass2",
            ms(ns_p2_emit),
            pct(ns_p2_emit, pass2)
        );
    }

    Ok((checksum, nreads))
}

// ── Paired bucket entry layout ────────────────────────────────────────────────
// R1: [key u64][apos u16][flip u8][pair_index u64]
//     [hlen u32][header][slen u32][canon_seq][qlen u32][qual]
// R2: [flip u8][hlen u32][header][slen u32][canon_seq][qlen u32][qual]

fn write_pair_entry(
    w: &mut BufWriter<File>,
    key: u64,
    apos: u16,
    pair_index: u64,
    r1_flip: bool,
    r1_hdr: &[u8],
    r1_canon: &[u8],
    r1_qual: &[u8],
    r2_flip: bool,
    r2_hdr: &[u8],
    r2_canon: &[u8],
    r2_qual: &[u8],
) -> Result<()> {
    w.write_all(&key.to_le_bytes())?;
    w.write_all(&apos.to_le_bytes())?;
    w.write_all(&[r1_flip as u8])?;
    w.write_all(&pair_index.to_le_bytes())?;
    for field in [r1_hdr, r1_canon, r1_qual] {
        w.write_all(&(field.len() as u32).to_le_bytes())?;
        w.write_all(field)?;
    }
    w.write_all(&[r2_flip as u8])?;
    for field in [r2_hdr, r2_canon, r2_qual] {
        w.write_all(&(field.len() as u32).to_le_bytes())?;
        w.write_all(field)?;
    }
    Ok(())
}

struct PairEntry {
    key: u64,
    apos: u16,
    pair_index: u64,
    r1_flip: bool,
    r1_header: Vec<u8>,
    r1_canon: Vec<u8>,
    r1_qual: Vec<u8>,
    r2_flip: bool,
    r2_header: Vec<u8>,
    r2_canon: Vec<u8>,
    r2_qual: Vec<u8>,
}

fn parse_pair_bucket(data: &[u8]) -> Result<Vec<PairEntry>> {
    let mut ents = Vec::new();
    let mut p = 0usize;

    // Fixed overhead per entry (variable fields' bytes come on top, guarded by the
    // ensure!s below): R1 prefix key(8)+apos(2)+flip(1)+pair_index(8) = 19, R1 length
    // fields hlen(4)+slen(4)+qlen(4) = 12, R2 flip(1) + length fields 12 = 13.
    // 19 + 12 + 13 = 44 bytes of fixed framing before any field payload.
    while p + 44 <= data.len() {
        let key = u64::from_le_bytes(
            data[p..p + 8].try_into().map_err(|_| anyhow!("pair bucket parse: key slice"))?,
        );
        let apos = u16::from_le_bytes(
            data[p + 8..p + 10]
                .try_into()
                .map_err(|_| anyhow!("pair bucket parse: apos slice"))?,
        );
        let r1_flip = data[p + 10] != 0;
        let pair_index = u64::from_le_bytes(
            data[p + 11..p + 19]
                .try_into()
                .map_err(|_| anyhow!("pair bucket parse: pair_index slice"))?,
        );
        p += 19;

        // R1 fields: header, canon_seq, qual
        let hlen = u32::from_le_bytes(
            data[p..p + 4].try_into().map_err(|_| anyhow!("pair bucket parse: r1 hlen"))?,
        ) as usize;
        p += 4;
        ensure!(p + hlen <= data.len(), "pair bucket parse: r1 header truncated");
        let r1_header = data[p..p + hlen].to_vec();
        p += hlen;

        let slen = u32::from_le_bytes(
            data[p..p + 4].try_into().map_err(|_| anyhow!("pair bucket parse: r1 slen"))?,
        ) as usize;
        p += 4;
        ensure!(p + slen <= data.len(), "pair bucket parse: r1 seq truncated");
        let r1_canon = data[p..p + slen].to_vec();
        p += slen;

        let qlen = u32::from_le_bytes(
            data[p..p + 4].try_into().map_err(|_| anyhow!("pair bucket parse: r1 qlen"))?,
        ) as usize;
        p += 4;
        ensure!(p + qlen <= data.len(), "pair bucket parse: r1 qual truncated");
        let r1_qual = data[p..p + qlen].to_vec();
        p += qlen;

        // R2 fields: flip, header, canon_seq, qual
        ensure!(p < data.len(), "pair bucket parse: r2 flip truncated");
        let r2_flip = data[p] != 0;
        p += 1;

        let hlen2 = u32::from_le_bytes(
            data[p..p + 4].try_into().map_err(|_| anyhow!("pair bucket parse: r2 hlen"))?,
        ) as usize;
        p += 4;
        ensure!(p + hlen2 <= data.len(), "pair bucket parse: r2 header truncated");
        let r2_header = data[p..p + hlen2].to_vec();
        p += hlen2;

        let slen2 = u32::from_le_bytes(
            data[p..p + 4].try_into().map_err(|_| anyhow!("pair bucket parse: r2 slen"))?,
        ) as usize;
        p += 4;
        ensure!(p + slen2 <= data.len(), "pair bucket parse: r2 seq truncated");
        let r2_canon = data[p..p + slen2].to_vec();
        p += slen2;

        let qlen2 = u32::from_le_bytes(
            data[p..p + 4].try_into().map_err(|_| anyhow!("pair bucket parse: r2 qlen"))?,
        ) as usize;
        p += 4;
        ensure!(p + qlen2 <= data.len(), "pair bucket parse: r2 qual truncated");
        let r2_qual = data[p..p + qlen2].to_vec();
        p += qlen2;

        ents.push(PairEntry {
            key,
            apos,
            pair_index,
            r1_flip,
            r1_header,
            r1_canon,
            r1_qual,
            r2_flip,
            r2_header,
            r2_canon,
            r2_qual,
        });
    }
    Ok(ents)
}

/// Cluster-sort paired-end records from `r1_path` + `r2_path` (plain or gzip FASTQ)
/// using on-disk buckets under `work/`, keyed by R1's strand-symmetric syncmer hash.
/// Both mates travel together as a unit — pairing is preserved with NO order stream.
/// R1 and R2 are each independently canonicalized (each may flip or not, separately).
///
/// Calls `on_pair(&r1_record, &r2_record)` for each pair in cluster order.
///
/// Returns `(verify_checksum, n_pairs)` where `verify_checksum` is the order-independent
/// wrapping sum of `hash_pair(r1_hdr, r1_orig_seq, r1_qual, r2_hdr, r2_orig_seq, r2_qual)`
/// over all pairs.
///
/// **Bounded RAM**: at most one bucket's raw bytes plus one bucket's parsed entries are
/// resident at a time. Pass 1 streams record-by-record.
pub(crate) fn cluster_sort_paired(
    r1_path: &Path,
    r2_path: &Path,
    work: &Path,
    mut on_pair: impl FnMut(&ClusteredRecord, &ClusteredRecord) -> Result<()>,
) -> Result<(u128, u64)> {
    std::fs::create_dir_all(work)
        .with_context(|| format!("create work dir {:?}", work))?;

    let bucket_dir = tempfile::TempDir::new_in(work)
        .with_context(|| format!("create cluster bucket scratch dir in {:?}", work))?;
    let bdir = bucket_dir.path();

    // Optional phase timing (QZ_CLUSTER_TIMING=1). Zero cost when unset.
    let timing = std::env::var("QZ_CLUSTER_TIMING").is_ok();

    // ── Pass 0: splitters from R1 only ───────────────────────────────────────
    let t_pass0 = std::time::Instant::now();
    let splitters = compute_splitters(r1_path)?;
    let pass0 = t_pass0.elapsed();

    // ── Open all bucket writers ───────────────────────────────────────────────
    let mut writers: Vec<BufWriter<File>> = (0..NBUCKETS)
        .map(|i| {
            let path = bdir.join(format!("b{i:04}"));
            File::create(&path)
                .with_context(|| format!("create bucket {:?}", path))
                .map(|f| BufWriter::with_capacity(1 << 18, f))
        })
        .collect::<Result<Vec<_>>>()?;

    // ── Pass 1: batched lockstep read R1+R2 → parallel canon/key/hash → bucket ──
    // Mirrors single-end's batched-parallel Pass 1 (`cluster_sort_single`): a serial
    // lockstep reader stages one batch of pairs into a flat arena, the per-pair
    // canonicalization + R1 syncmer key + pair checksum are computed in parallel across
    // the batch, then the checksum is folded and bucket entries written sequentially.
    // `pair_index` is assigned in input order, and Pass 2's total-order comparator
    // (key → R1 anchor-suffix → R1 canon → pair_index) makes the clustered emission order
    // independent of intra-bucket write order — so this is byte-identical to the old
    // serial Pass 1 (the parallel canon is pure; the `wrapping_add` checksum fold is
    // commutative).
    let t_pass1 = std::time::Instant::now();
    let mut r1_reader = open_maybe_gz(r1_path)?;
    let mut r2_reader = open_maybe_gz(r2_path)?;
    let mut n_pairs: u64 = 0;
    let mut checksum: u128 = 0;

    // One pair's byte ranges within the batch arena: [r1_h][r1_s][r1_q][r2_h][r2_s][r2_q].
    struct PairRanges {
        pair_index: u64,
        r1_h: (usize, usize),
        r1_s: (usize, usize),
        r1_q: (usize, usize),
        r2_h: (usize, usize),
        r2_s: (usize, usize),
        r2_q: (usize, usize),
    }
    // Per-pair parallel output: both mates' canonicalization + R1 key + pair checksum.
    struct PairComputed {
        key: u64,
        apos: u16,
        r1_flip: bool,
        r1_canon: Vec<u8>,
        r2_flip: bool,
        r2_canon: Vec<u8>,
        hash: u128,
    }

    // Arena stores consecutive pairs' raw bytes; `meta` records each pair's ranges.
    let mut arena: Vec<u8> = Vec::with_capacity(PAIR_BATCH * 1024);
    let mut meta: Vec<PairRanges> = Vec::with_capacity(PAIR_BATCH);

    loop {
        arena.clear();
        meta.clear();

        // Fill a batch by reading both mates in lockstep (serial — sequential file I/O).
        while meta.len() < PAIR_BATCH {
            let r1_rec = read_record(&mut *r1_reader)?;
            let r2_rec = read_record(&mut *r2_reader)?;
            match (r1_rec, r2_rec) {
                (Some(r1), Some(r2)) => {
                    let pair_index = n_pairs;
                    n_pairs += 1;

                    // `push` mut-borrows `arena`; the slices passed in borrow the owned
                    // `Record`s (`r1`/`r2`), not `arena`, so there is no aliasing. The
                    // borrow ends (NLL) at the last `push` call below, before `meta.push`.
                    let mut push = |bytes: &[u8]| {
                        let start = arena.len();
                        arena.extend_from_slice(bytes);
                        (start, arena.len())
                    };
                    let r1_h = push(r1.header());
                    let r1_s = push(r1.seq());
                    let r1_q = push(r1.qual());
                    let r2_h = push(r2.header());
                    let r2_s = push(r2.seq());
                    let r2_q = push(r2.qual());

                    meta.push(PairRanges { pair_index, r1_h, r1_s, r1_q, r2_h, r2_s, r2_q });
                }
                (None, None) => break,
                (Some(_), None) | (None, Some(_)) => {
                    bail!("paired inputs have unequal read counts (diverged at pair {n_pairs})");
                }
            }
        }

        if meta.is_empty() {
            break;
        }

        // Parallel: canonicalize both mates, compute R1's syncmer key, and hash the pair
        // over its ORIGINAL (pre-canon) bytes (order-independent — folded sequentially).
        let computed: Vec<PairComputed> = meta
            .par_iter()
            .map(|m| {
                let r1_seq = &arena[m.r1_s.0..m.r1_s.1];
                let r2_seq = &arena[m.r2_s.0..m.r2_s.1];

                let (r1_canon, r1_flip, key, apos) =
                    canonicalize_open_syncmer_strand(r1_seq);

                let (r2_canon, r2_flip, _r2_key, _r2_apos) =
                    canonicalize_open_syncmer_strand(r2_seq);

                let hash = hash_pair(
                    &arena[m.r1_h.0..m.r1_h.1],
                    r1_seq,
                    &arena[m.r1_q.0..m.r1_q.1],
                    &arena[m.r2_h.0..m.r2_h.1],
                    r2_seq,
                    &arena[m.r2_q.0..m.r2_q.1],
                );

                PairComputed { key, apos, r1_flip, r1_canon, r2_flip, r2_canon, hash }
            })
            .collect();

        // Sequential: fold checksum (commutative) + write bucket entries to shared writers.
        for (c, m) in computed.iter().zip(meta.iter()) {
            checksum = checksum.wrapping_add(c.hash);
            let bucket = splitters.partition_point(|&s| s <= c.key);
            write_pair_entry(
                &mut writers[bucket],
                c.key,
                c.apos,
                m.pair_index,
                c.r1_flip,
                &arena[m.r1_h.0..m.r1_h.1],
                &c.r1_canon,
                &arena[m.r1_q.0..m.r1_q.1],
                c.r2_flip,
                &arena[m.r2_h.0..m.r2_h.1],
                &c.r2_canon,
                &arena[m.r2_q.0..m.r2_q.1],
            )?;
        }
    }

    // Free the batch working set (arena + ranges) before Pass 2 + the codec flush so
    // it isn't dead weight under the chunk-accumulation/flush RSS peak (bounded RAM).
    drop(arena);
    drop(meta);

    // Flush and drop all writers before Pass 2 reads the files.
    for w in &mut writers {
        w.flush()?;
    }
    drop(writers);
    let pass1 = t_pass1.elapsed();

    // ── Pass 2: within-bucket sort → stream pairs ────────────────────────────
    let t_pass2 = std::time::Instant::now();
    let mut emitted: u64 = 0;
    for i in 0..NBUCKETS {
        let path = bdir.join(format!("b{i:04}"));
        let data =
            std::fs::read(&path).with_context(|| format!("read bucket {:?}", path))?;

        if data.is_empty() {
            continue;
        }

        let mut ents = parse_pair_bucket(&data)?;

        // Total-order comparator: key → R1 anchor-suffix → full R1 canon → pair_index.
        // Mirrors single-end's comparator shape exactly.
        ents.sort_unstable_by(|a, b| {
            a.key.cmp(&b.key).then_with(|| {
                let sa_start = (a.apos as usize).min(a.r1_canon.len());
                let sb_start = (b.apos as usize).min(b.r1_canon.len());
                let sa = &a.r1_canon[sa_start..];
                let sb = &b.r1_canon[sb_start..];
                sa.cmp(sb)
                    .then_with(|| a.r1_canon.cmp(&b.r1_canon))
                    .then_with(|| a.pair_index.cmp(&b.pair_index))
            })
        });

        // Consume the parsed entries (move fields into the records — no per-pair clone).
        for e in ents {
            let r1_rec = ClusteredRecord {
                header: e.r1_header,
                canon_seq: e.r1_canon,
                qual: e.r1_qual,
                flip: e.r1_flip,
                apos: e.apos,
                key: e.key,
                original_index: e.pair_index,
            };
            let r2_rec = ClusteredRecord {
                header: e.r2_header,
                canon_seq: e.r2_canon,
                qual: e.r2_qual,
                flip: e.r2_flip,
                apos: 0,
                key: 0,
                original_index: e.pair_index,
            };
            on_pair(&r1_rec, &r2_rec)?;
            emitted += 1;
        }

        // Drop the raw bucket bytes before the next bucket (bounded RAM — one resident).
        drop(data);
    }

    ensure!(emitted == n_pairs, "emitted {emitted} != n_pairs {n_pairs}: bucket accounting error");

    if timing {
        let pass2 = t_pass2.elapsed();
        eprintln!(
            "[CLUSTER_TIMING paired sort] pass0(splitter sample)={:.1}s  pass1(batched-parallel canon+hash+bucket)={:.1}s  pass2(parse+sort+emit_callback)={:.1}s  n_pairs={n_pairs}",
            pass0.as_secs_f64(),
            pass1.as_secs_f64(),
            pass2.as_secs_f64(),
        );
    }

    Ok((checksum, n_pairs))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    #[test]
    fn cluster_sort_preserves_record_multiset_streaming() {
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("in.fastq");
        let mut f = std::fs::File::create(&fq).unwrap();
        for (h, s) in [("@a","ACGTACGTAC"),("@b","GTACGTACGT"),("@c","ACGTACGTAC"),
                       ("@d","TTTTAAAACC"),("@e","GGTTTTAAAA"),("@f","ACGTACGTAC")] {
            writeln!(f, "{h}\n{s}\n+\n{}", "I".repeat(s.len())).unwrap();
        }
        drop(f);
        let work = dir.path().join("work");
        let mut got: Vec<(Vec<u8>,Vec<u8>,Vec<u8>)> = Vec::new();
        let (checksum, n) = cluster_sort_single(&fq, &work, |r| {
            got.push((r.header.clone(), r.original_seq(), r.qual.clone())); Ok(())
        }).unwrap();
        assert_eq!(n, 6);
        got.sort();
        let mut want = vec![
            (b"@a".to_vec(), b"ACGTACGTAC".to_vec(), b"IIIIIIIIII".to_vec()),
            (b"@b".to_vec(), b"GTACGTACGT".to_vec(), b"IIIIIIIIII".to_vec()),
            (b"@c".to_vec(), b"ACGTACGTAC".to_vec(), b"IIIIIIIIII".to_vec()),
            (b"@d".to_vec(), b"TTTTAAAACC".to_vec(), b"IIIIIIIIII".to_vec()),
            (b"@e".to_vec(), b"GGTTTTAAAA".to_vec(), b"IIIIIIIIII".to_vec()),
            (b"@f".to_vec(), b"ACGTACGTAC".to_vec(), b"IIIIIIIIII".to_vec()),
        ];
        want.sort();
        assert_eq!(got, want);
        assert_ne!(checksum, 0);
    }

    #[test]
    fn cluster_sort_handles_empty_input() {
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("empty.fastq");
        std::fs::write(&fq, b"").unwrap();
        let mut count = 0u64;
        let (_checksum, n) = cluster_sort_single(&fq, &dir.path().join("w"), |_r| { count += 1; Ok(()) }).unwrap();
        assert_eq!(n, 0);
        assert_eq!(count, 0);
    }

    #[test]
    fn cluster_sort_paired_preserves_pair_multiset_and_pairing() {
        let dir = tempfile::tempdir().unwrap();
        let r1p = dir.path().join("r1.fastq");
        let r2p = dir.path().join("r2.fastq");
        let mut f1 = std::fs::File::create(&r1p).unwrap();
        let mut f2 = std::fs::File::create(&r2p).unwrap();
        // 6 pairs; R2 seq deterministically tied to its R1 so we can assert pairing.
        let r1s = ["ACGTACGTAC","GTACGTACGT","ACGTACGTAC","TTTTAAAACC","GGTTTTAAAA","ACGTACGTAC"];
        for (i, s) in r1s.iter().enumerate() {
            writeln!(f1, "@p{i}/1\n{s}\n+\n{}", "I".repeat(s.len())).unwrap();
            let r2: String = s.chars().rev().collect(); // R2 = reverse of R1 (just a marker)
            writeln!(f2, "@p{i}/2\n{r2}\n+\n{}", "J".repeat(r2.len())).unwrap();
        }
        drop(f1); drop(f2);
        let mut got: Vec<(Vec<u8>,Vec<u8>,Vec<u8>,Vec<u8>)> = Vec::new();
        let (checksum, n) = cluster_sort_paired(&r1p, &r2p, &dir.path().join("w"), |r1, r2| {
            // record (R1.header, R1.original_seq, R2.header, R2.original_seq) — pairing must hold
            got.push((r1.header.clone(), r1.original_seq(), r2.header.clone(), r2.original_seq()));
            Ok(())
        }).unwrap();
        assert_eq!(n, 6);
        assert_ne!(checksum, 0);
        // every emitted pair has matching indices in its two headers (e.g. @pK/1 with @pK/2)
        for (h1, _s1, h2, _s2) in &got {
            let k1 = std::str::from_utf8(h1).unwrap().trim_start_matches("@p").split('/').next().unwrap().to_string();
            let k2 = std::str::from_utf8(h2).unwrap().trim_start_matches("@p").split('/').next().unwrap().to_string();
            assert_eq!(k1, k2, "R1/R2 must stay paired");
        }
        // pair multiset preserved
        let mut keys: Vec<String> = got.iter()
            .map(|(h1,_,_,_)| std::str::from_utf8(h1).unwrap().to_string()).collect();
        keys.sort();
        assert_eq!(keys, vec!["@p0/1","@p1/1","@p2/1","@p3/1","@p4/1","@p5/1"]);
    }

    #[test]
    fn all_identical_reads_are_deterministically_ordered() {
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("in.fastq");
        let mut f = std::fs::File::create(&fq).unwrap();
        for i in 0..50 { writeln!(f, "@r{i}\nACGTACGTACGTACGTACGT\n+\n{}", "I".repeat(20)).unwrap(); }
        drop(f);
        let run = || {
            let mut hdrs = Vec::new();
            cluster_sort_single(&fq, &dir.path().join("w"), |r| { hdrs.push(r.header.clone()); Ok(()) }).unwrap();
            hdrs
        };
        assert_eq!(run(), run(), "identical-read order must be deterministic");
    }
}
