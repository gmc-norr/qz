//! Record-aligned input splitting for NUMA-compress (spec §4). Single-end:
//! seek-and-validate (near-zero I/O). Paired: record-index lockstep (Task 11).
//! Precondition (enforced by §3.1 prechecks before we run): seekable plaintext
//! FASTQ — never stdin, gzip, or FASTA.

use anyhow::{Result, bail};
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::Path;

pub use qz_lib::compression::ByteRange;

/// A paired split: matching record-index ranges in R1 and R2.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PairedRange {
    pub r1: ByteRange,
    pub r2: ByteRange,
}

/// Sentinel error: the input has too few records to make `n` non-empty ranges.
/// The driver handles it mode-sensitively (auto clamps/falls back; fixed errors).
#[derive(Debug)]
pub struct SplitTooFewRecords;
impl std::fmt::Display for SplitTooFewRecords {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "input has too few records to split into the requested number of ranges")
    }
}
impl std::error::Error for SplitTooFewRecords {}

/// Sentinel error: paired mates have different record counts (cannot lockstep-shard).
/// The driver handles it mode-sensitively (auto → in-process; fixed → clear error).
#[derive(Debug)]
pub struct MateCountMismatch;
impl std::fmt::Display for MateCountMismatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "paired inputs have different record counts (cannot shard)")
    }
}
impl std::error::Error for MateCountMismatch {}

/// Largest plausible single FASTQ record, in bytes (mirrors qz-lib
/// `PER_RECORD_MAX_BYTES` = 63 MiB). Used only as a sanity bound on scans.
const PER_RECORD_MAX_BYTES: u64 = 63 << 20;

/// Split a seekable plaintext FASTQ into EXACTLY `n` contiguous record-aligned
/// ranges (spec §4.2). `n >= 2`.
pub fn split_single_end(path: &Path, n: u32) -> Result<Vec<ByteRange>> {
    debug_assert!(n >= 2, "split_single_end requires n >= 2 (the gate guarantees this)");
    let filesize = std::fs::metadata(path)?.len();
    if filesize == 0 {
        return Err(SplitTooFewRecords.into());
    }
    let mut file = File::open(path)?;

    // Find a strictly-increasing boundary at/after each target t_k = k*filesize/n.
    let mut boundaries: Vec<u64> = Vec::with_capacity(n as usize + 1);
    boundaries.push(0);
    for k in 1..n {
        let target = (k as u64) * filesize / (n as u64);
        let want = target.max(boundaries.last().copied().unwrap() + 1);
        if want >= filesize {
            return Err(SplitTooFewRecords.into());
        }
        let b = find_record_boundary_at_or_after(&mut file, want, filesize)?;
        // Strictly increasing? If a boundary collapses onto the previous one (or to
        // EOF), there aren't enough records for n ranges.
        if b <= boundaries.last().copied().unwrap() || b >= filesize {
            return Err(SplitTooFewRecords.into());
        }
        boundaries.push(b);
    }
    boundaries.push(filesize);

    let ranges: Vec<ByteRange> = boundaries
        .windows(2)
        .map(|w| ByteRange { start: w[0], end: w[1] })
        .collect();
    debug_assert_eq!(ranges.len(), n as usize);
    Ok(ranges)
}

/// Split both mates into EXACTLY `n` record-index-aligned paired ranges. Both mates must
/// have the same record count (else `MateCountMismatch`); too few records →
/// `SplitTooFewRecords`. Record i of R1 and record i of R2 always land in the same range, so
/// a pair is never split across workers. `n >= 2`.
pub fn split_paired(r1_path: &Path, r2_path: &Path, n: u32) -> Result<Vec<PairedRange>> {
    debug_assert!(n >= 2, "split_paired requires n >= 2 (the gate guarantees this)");
    // Scan both mates CONCURRENTLY — each is a full memchr pass (the dominant serial cost in
    // the paired driver). Independent, so a scoped thread halves it (bandwidth-bound → not a
    // full 2×, but worthwhile). The offset Vecs are freed at the end of this fn, before the
    // driver forks workers.
    let (res1, res2) = std::thread::scope(|s| {
        let h1 = s.spawn(|| record_start_offsets(r1_path));
        let r2 = record_start_offsets(r2_path);
        (h1.join(), r2)
    });
    let (off1, fs1) = res1.map_err(|_| anyhow::anyhow!("paired R1 scan thread panicked"))??;
    let (off2, fs2) = res2?;
    let total = off1.len() as u64;
    if off2.len() as u64 != total {
        return Err(MateCountMismatch.into());
    }
    if total < n as u64 {
        return Err(SplitTooFewRecords.into());
    }
    // Interior boundary record indices: round(k/n * total), strictly increasing in (0,total).
    let mut interior: Vec<u64> = Vec::with_capacity(n as usize - 1);
    let mut prev = 0u64;
    for k in 1..n {
        let want = (((k as u128) * (total as u128) + (n as u128) / 2) / (n as u128)) as u64;
        let want = want.max(prev + 1);
        if want >= total {
            return Err(SplitTooFewRecords.into());
        }
        interior.push(want);
        prev = want;
    }
    // Stitch full boundary arrays: [0, off[idx]…, filesize] for each mate (same indices).
    let starts1: Vec<u64> = std::iter::once(0)
        .chain(interior.iter().map(|&i| off1[i as usize]))
        .chain(std::iter::once(fs1))
        .collect();
    let starts2: Vec<u64> = std::iter::once(0)
        .chain(interior.iter().map(|&i| off2[i as usize]))
        .chain(std::iter::once(fs2))
        .collect();
    let mut ranges = Vec::with_capacity(n as usize);
    for i in 0..n as usize {
        ranges.push(PairedRange {
            r1: ByteRange { start: starts1[i], end: starts1[i + 1] },
            r2: ByteRange { start: starts2[i], end: starts2[i + 1] },
        });
    }
    debug_assert_eq!(ranges.len(), n as usize);
    Ok(ranges)
}

/// Byte offsets of every record start in a plaintext FASTQ (record i = line 4i), via a
/// SINGLE `memchr` newline scan — NO per-line `read_until`, NO record materialization (the
/// old paired splitter read both whole files twice with `read_until` per line: the dominant
/// serial cost that tanked paired NUMA compress). Returns `(offsets, filesize)` with
/// `offsets.len()` == record count and `offsets[i]` the byte offset of record i's header
/// line. Errors if the line count is not a multiple of 4 (truncated). The offset Vec is
/// O(records) (~8 bytes each) and transient — the driver frees it before forking workers.
fn record_start_offsets(path: &Path) -> Result<(Vec<u64>, u64)> {
    let file = File::open(path)?;
    let filesize = file.metadata()?.len();
    if filesize == 0 {
        return Ok((Vec::new(), 0));
    }
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let mut offsets: Vec<u64> = vec![0]; // record 0 starts at byte 0
    let mut base: u64 = 0; // absolute byte offset of the current fill_buf window
    let mut nl_count: u64 = 0; // newlines seen so far
    let mut last_line_end: u64 = 0; // absolute byte offset just after the last newline
    loop {
        let buf = reader.fill_buf()?;
        if buf.is_empty() {
            break;
        }
        let n = buf.len();
        for nl in memchr::memchr_iter(b'\n', buf) {
            nl_count += 1;
            let after = base + nl as u64 + 1; // byte just after this newline
            last_line_end = after;
            // Every 4th line start begins a new record; skip a target landing at EOF.
            if nl_count.is_multiple_of(4) && after < filesize {
                offsets.push(after);
            }
        }
        base += n as u64;
        reader.consume(n);
    }
    // Total lines = completed (newline-terminated) lines + a trailing unterminated line, if any.
    let total_lines = nl_count + if last_line_end < filesize { 1 } else { 0 };
    if !total_lines.is_multiple_of(4) {
        bail!("FASTQ {} does not end on a record boundary (truncated?)", path.display());
    }
    let records = (total_lines / 4) as usize;
    if offsets.len() != records {
        // Defensive — unreachable given the %4 invariant above.
        bail!("FASTQ {} record-start scan inconsistency ({} vs {})", path.display(), offsets.len(), records);
    }
    Ok((offsets, filesize))
}

/// Return the byte offset of the next FASTQ record boundary at or after `pos`.
/// A boundary is a line start `L` where: line(L) begins `@`, line(L+2) begins `+`,
/// stripped len(line(L+1)) == stripped len(line(L+3)), AND the same 4-line shape
/// repeats at L+4 (two consecutive records — a single `@`/`+` is not reliable
/// because quality bytes can be `@`/`+`). EOF/last-record contract: accept a
/// 1-record validation when L+4 is EOF. Returns `filesize` if nothing validates.
fn find_record_boundary_at_or_after(file: &mut File, pos: u64, filesize: u64) -> Result<u64> {
    // Read a window starting at the line boundary at/after `pos`. Bounded by twice
    // the max record (validating two records) plus slack.
    let mut l = next_line_start(file, pos, filesize)?;
    while l < filesize {
        match validate_two_records(file, l, filesize)? {
            BoundaryCheck::Ok => return Ok(l),
            BoundaryCheck::Advance(next_l) => {
                if next_l <= l { return Ok(filesize); }
                l = next_l;
            }
        }
    }
    Ok(filesize)
}

enum BoundaryCheck { Ok, Advance(u64) }

/// Read the 4 lines of a record at `l` (and peek the next record's header) and
/// decide whether `l` is a valid boundary.
fn validate_two_records(file: &mut File, l: u64, filesize: u64) -> Result<BoundaryCheck> {
    let lines = read_lines_from(file, l, 8)?; // up to 2 records of 4 lines
    if lines.is_empty() {
        return Ok(BoundaryCheck::Advance(filesize));
    }
    // Helper closures over the (line_start_offset, bytes) pairs in `lines`.
    let line_bytes = |i: usize| -> Option<&Vec<u8>> { lines.get(i).map(|(_, b)| b) };
    let stripped_len = |i: usize| -> Option<usize> {
        line_bytes(i).map(|b| {
            let mut n = b.len();
            while n > 0 && (b[n - 1] == b'\n' || b[n - 1] == b'\r') { n -= 1; }
            n
        })
    };
    let begins = |i: usize, c: u8| line_bytes(i).and_then(|b| b.first().copied()) == Some(c);

    // Record 1 shape.
    if !begins(0, b'@') || !begins(2, b'+') {
        return Ok(BoundaryCheck::Advance(line_offset_after(&lines, 0)));
    }
    if stripped_len(1) != stripped_len(3) {
        return Ok(BoundaryCheck::Advance(line_offset_after(&lines, 0)));
    }
    // Two-record confirmation OR EOF acceptance.
    let after_rec1 = line_offset_after(&lines, 3);
    if after_rec1 >= filesize {
        // L+4 is EOF — accept a 1-record validation (last record not lost).
        return Ok(BoundaryCheck::Ok);
    }
    if begins(4, b'@') && begins(6, b'+') && stripped_len(5) == stripped_len(7) {
        return Ok(BoundaryCheck::Ok);
    }
    Ok(BoundaryCheck::Advance(line_offset_after(&lines, 0)))
}

/// Byte offset of the next line start strictly after the start of line index `i`
/// within `lines` (used to advance the scan by one line). Falls back to filesize.
fn line_offset_after(lines: &[(u64, Vec<u8>)], i: usize) -> u64 {
    match lines.get(i + 1) {
        Some((off, _)) => *off,
        None => match lines.get(i) {
            Some((off, b)) => *off + b.len() as u64,
            None => u64::MAX,
        },
    }
}

/// Offset of the first line start at/after `pos`. If `pos` is ALREADY a line start
/// (file start, or the byte after a `\n`), return `pos` unchanged — otherwise an
/// exact-boundary target (e.g. `n == records`, or a uniform-length file where
/// `k*filesize/n` lands on a record start) would be skipped into the record and
/// spuriously collapse to `SplitTooFewRecords`. Else scan forward to the byte after
/// the next `\n`.
fn next_line_start(file: &mut File, pos: u64, filesize: u64) -> Result<u64> {
    if pos == 0 { return Ok(0); }
    // Peek the byte before `pos`: if it's '\n', `pos` already starts a line.
    file.seek(SeekFrom::Start(pos - 1))?;
    let mut prev = [0u8; 1];
    let n = std::io::Read::read(file, &mut prev)?;
    if n == 1 && prev[0] == b'\n' {
        return Ok(pos);
    }
    // Not a line start: scan forward to the next '\n' and return the byte after it.
    file.seek(SeekFrom::Start(pos))?;
    let mut reader = BufReader::with_capacity(64 << 10, file);
    let mut buf = Vec::new();
    let read = reader.read_until(b'\n', &mut buf)?;
    let start = pos + read as u64;
    Ok(start.min(filesize))
}

/// Read up to `max_lines` lines starting at byte `from`, returning each as
/// `(line_start_offset, bytes_including_newline)`. Stops at EOF.
fn read_lines_from(file: &mut File, from: u64, max_lines: usize) -> Result<Vec<(u64, Vec<u8>)>> {
    file.seek(SeekFrom::Start(from))?;
    let mut reader = BufReader::with_capacity(256 << 10, file);
    let mut out = Vec::with_capacity(max_lines);
    let mut off = from;
    for _ in 0..max_lines {
        let mut buf = Vec::new();
        let n = reader.read_until(b'\n', &mut buf)?;
        if n == 0 { break; }
        if buf.len() as u64 > PER_RECORD_MAX_BYTES {
            bail!("FASTQ line at offset {off} exceeds the {PER_RECORD_MAX_BYTES}-byte record cap");
        }
        let start = off;
        off += n as u64;
        out.push((start, buf));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Brute-force oracle: byte offsets of every record start in a plaintext FASTQ.
    fn record_starts(text: &str) -> Vec<u64> {
        let bytes = text.as_bytes();
        let mut starts = Vec::new();
        let mut i = 0usize;
        let mut off = 0u64;
        let lines: Vec<&[u8]> = {
            // split keeping offsets
            let mut v = Vec::new();
            let mut s = 0;
            for (j, &b) in bytes.iter().enumerate() {
                if b == b'\n' { v.push((s, j + 1)); s = j + 1; }
            }
            if s < bytes.len() { v.push((s, bytes.len())); }
            v.iter().map(|&(a, b)| &bytes[a..b]).collect()
        };
        while i < lines.len() {
            starts.push(off);
            for k in 0..4 { off += lines.get(i + k).map(|l| l.len() as u64).unwrap_or(0); }
            i += 4;
        }
        starts
    }

    fn fastq(records: &[(&str, &str, &str)]) -> String {
        records.iter().map(|(h, s, q)| format!("@{h}\n{s}\n+\n{q}\n")).collect()
    }

    fn assert_partition(ranges: &[ByteRange], n: u32, filesize: u64) {
        assert_eq!(ranges.len() as u32, n);
        assert_eq!(ranges[0].start, 0);
        assert_eq!(ranges.last().unwrap().end, filesize);
        for w in ranges.windows(2) {
            assert_eq!(w[0].end, w[1].start, "contiguous");
        }
        assert!(ranges.iter().all(|r| r.end > r.start), "non-empty");
    }

    fn boundaries_are_record_starts(ranges: &[ByteRange], starts: &[u64]) {
        for r in &ranges[1..] {
            assert!(starts.contains(&r.start), "boundary {} must be a record start", r.start);
        }
    }

    #[test]
    fn splits_fixed_length_fastq_on_record_boundaries() {
        let td = tempfile::TempDir::new().unwrap();
        let p = td.path().join("f.fastq");
        let recs: Vec<(String, String, String)> = (0..40)
            .map(|i| (format!("r{i}"), "ACGTACGTACGT".to_string(), "IIIIIIIIIIII".to_string()))
            .collect();
        let text: String = recs.iter().map(|(h, s, q)| format!("@{h}\n{s}\n+\n{q}\n")).collect();
        std::fs::write(&p, &text).unwrap();
        let fs = text.len() as u64;
        let starts = record_starts(&text);
        for n in [2u32, 3, 4, 5] {
            let r = split_single_end(&p, n).unwrap();
            assert_partition(&r, n, fs);
            boundaries_are_record_starts(&r, &starts);
        }
    }

    #[test]
    fn splits_when_quality_lines_start_with_at_and_plus() {
        let td = tempfile::TempDir::new().unwrap();
        let p = td.path().join("trap.fastq");
        // Quality strings deliberately start with '@' and '+' — a single-char check
        // would mis-fire; the 2-record length-equality validator must not.
        let text = fastq(&[
            ("r0", "ACGTAC", "@IIII@"),
            ("r1", "TTGGCC", "+IabI+"),
            ("r2", "GGCCAA", "@@++@@"),
            ("r3", "AACCTT", "IIIIII"),
            ("r4", "CCGGTT", "++++++"),
            ("r5", "TTAAGG", "@a@a@a"),
        ]);
        std::fs::write(&p, &text).unwrap();
        let fs = text.len() as u64;
        let starts = record_starts(&text);
        let r = split_single_end(&p, 2).unwrap();
        assert_partition(&r, 2, fs);
        boundaries_are_record_starts(&r, &starts);
    }

    #[test]
    fn validator_rejects_quality_line_that_looks_like_a_header() {
        // rec0's QUALITY line begins with '@' (a naive single-'@' check would treat it as a
        // record start). Point the scan AT that quality line and confirm
        // find_record_boundary_at_or_after skips it to the TRUE next record boundary.
        let td = tempfile::TempDir::new().unwrap();
        let p = td.path().join("trap2.fastq");
        let rec0 = format!("@r0\n{}\n+\n@{}\n", "A".repeat(30), "I".repeat(29)); // qual starts '@'
        let rec1 = "@r1\nACGTACGT\n+\nIIIIIIII\n".to_string();
        let text = format!("{rec0}{rec1}");
        std::fs::write(&p, &text).unwrap();
        let filesize = text.len() as u64;
        let starts = record_starts(&text);
        // Offset of rec0's quality line: "@r0\n"=4, "A"*30+"\n"=31, "+\n"=2 → 37.
        let qual_off = 37u64;
        assert_eq!(text.as_bytes()[qual_off as usize], b'@', "quality line must start with '@'");
        let mut f = File::open(&p).unwrap();
        let b = find_record_boundary_at_or_after(&mut f, qual_off, filesize).unwrap();
        assert_eq!(b, starts[1], "scan must skip the quality-'@' line to the real record boundary");
        assert_ne!(b, qual_off, "must NOT mistake the quality-'@' line for a record start");
    }

    #[test]
    fn splits_crlf_fastq() {
        let td = tempfile::TempDir::new().unwrap();
        let p = td.path().join("crlf.fastq");
        let text: String = (0..16)
            .map(|i| format!("@r{i}\r\nACGTACGT\r\n+\r\nIIIIIIII\r\n"))
            .collect();
        std::fs::write(&p, &text).unwrap();
        let fs = text.len() as u64;
        let r = split_single_end(&p, 3).unwrap();
        assert_partition(&r, 3, fs);
        // Each boundary header line must begin with '@'.
        for rg in &r[1..] {
            assert_eq!(text.as_bytes()[rg.start as usize], b'@');
        }
    }

    #[test]
    fn too_few_records_errors() {
        let td = tempfile::TempDir::new().unwrap();
        let p = td.path().join("tiny.fastq");
        let text = fastq(&[("r0", "ACGT", "IIII"), ("r1", "TTGG", "IIII")]);
        std::fs::write(&p, &text).unwrap();
        let err = split_single_end(&p, 4).unwrap_err();
        assert!(err.downcast_ref::<SplitTooFewRecords>().is_some(), "want SplitTooFewRecords, got {err}");
    }

    #[test]
    fn split_target_in_last_record_does_not_spuriously_clamp() {
        // n=2 on a 3-record file where filesize/2 lands INSIDE the long first record;
        // the scan must advance forward to a real boundary, not clamp/lose records.
        let td = tempfile::TempDir::new().unwrap();
        let p = td.path().join("three.fastq");
        let text = fastq(&[
            ("r0", "ACGTACGTACGTACGT", "IIIIIIIIIIIIIIII"),
            ("r1", "TT", "II"),
            ("r2", "GG", "II"),
        ]);
        std::fs::write(&p, &text).unwrap();
        let fs = text.len() as u64;
        let r = split_single_end(&p, 2).unwrap();
        assert_partition(&r, 2, fs);
    }

    #[test]
    fn n_equals_records_yields_singletons_on_exact_boundaries() {
        // Uniform-length records → every target k*filesize/n lands EXACTLY on a
        // record start. The next_line_start exact-boundary fix must not skip
        // forward (which would collapse to SplitTooFewRecords).
        let td = tempfile::TempDir::new().unwrap();
        let p = td.path().join("uniform.fastq");
        let text: String = (0..6).map(|i| format!("@r{i}\nACGTACGT\n+\nIIIIIIII\n")).collect();
        std::fs::write(&p, &text).unwrap();
        let fs = text.len() as u64;
        let starts = record_starts(&text);
        let r = split_single_end(&p, 6).unwrap();
        assert_partition(&r, 6, fs);
        boundaries_are_record_starts(&r, &starts);
        assert!(r.iter().all(|rg| {
            // each range is exactly one record (start..next start)
            starts.contains(&rg.start)
        }));
        // A target landing exactly on a boundary must still split there, not later.
        let r2 = split_single_end(&p, 2).unwrap();
        assert_eq!(r2[0].end, starts[3], "n=2 on 6 uniform records splits at record 3");
    }

    #[test]
    fn record_start_offsets_matches_oracle() {
        let td = tempfile::TempDir::new().unwrap();
        // Plain LF and CRLF (both newline-terminated) must match the brute-force oracle.
        for (name, text) in [
            ("plain", (0..25).map(|i| format!("@r{i}\nACGTAC\n+\nIIIIII\n")).collect::<String>()),
            ("crlf", (0..16).map(|i| format!("@r{i}\r\nACGT\r\n+\r\nIIII\r\n")).collect::<String>()),
        ] {
            let p = td.path().join(format!("{name}.fastq"));
            std::fs::write(&p, &text).unwrap();
            let (offs, fs) = record_start_offsets(&p).unwrap();
            assert_eq!(fs, text.len() as u64, "{name} filesize");
            assert_eq!(offs, record_starts(&text), "{name} offsets");
        }
        // No trailing newline on the last qual line must still count R records.
        let p = td.path().join("notrail.fastq");
        let mut body: String = (0..7).map(|i| format!("@r{i}\nACGT\n+\nIIII\n")).collect();
        body.pop(); // drop the final '\n'
        std::fs::write(&p, &body).unwrap();
        let (offs, _) = record_start_offsets(&p).unwrap();
        assert_eq!(offs.len(), 7, "notrail record count");
        assert_eq!(offs, record_starts(&body), "notrail offsets");
        // A non-multiple-of-4 line count is a truncated input → error.
        let p = td.path().join("trunc.fastq");
        std::fs::write(&p, "@r0\nACGT\n+\nIIII\n@r1\nACGT\n").unwrap(); // 6 lines
        assert!(record_start_offsets(&p).is_err(), "truncated record must error");
    }

    #[test]
    fn split_paired_aligns_on_record_index() {
        // R1 and R2 have the SAME record count but DIFFERENT per-record byte sizes (10bp vs
        // 4bp), so byte offsets diverge while record indices must still align.
        let td = tempfile::TempDir::new().unwrap();
        let p1 = td.path().join("r1.fastq");
        let p2 = td.path().join("r2.fastq");
        let r1: String = (0..40).map(|i| format!("@r{i}/1\nACGTACGTAC\n+\nIIIIIIIIII\n")).collect();
        let r2: String = (0..40).map(|i| format!("@r{i}/2\nACGT\n+\nFFFF\n")).collect();
        std::fs::write(&p1, &r1).unwrap();
        std::fs::write(&p2, &r2).unwrap();
        let s1 = record_starts(&r1);
        let s2 = record_starts(&r2);
        for n in [2u32, 3, 4, 7] {
            let pr = split_paired(&p1, &p2, n).unwrap();
            assert_eq!(pr.len() as u32, n);
            assert_eq!(pr[0].r1.start, 0);
            assert_eq!(pr[0].r2.start, 0);
            assert_eq!(pr.last().unwrap().r1.end, r1.len() as u64);
            assert_eq!(pr.last().unwrap().r2.end, r2.len() as u64);
            for w in pr.windows(2) {
                assert_eq!(w[0].r1.end, w[1].r1.start, "contiguous r1");
                assert_eq!(w[0].r2.end, w[1].r2.start, "contiguous r2");
            }
            // Every interior boundary is a record start in BOTH mates, at the SAME index.
            for b in &pr[1..] {
                let i1 = s1.iter().position(|&o| o == b.r1.start).expect("r1 boundary at a record start");
                let i2 = s2.iter().position(|&o| o == b.r2.start).expect("r2 boundary at a record start");
                assert_eq!(i1, i2, "paired boundaries must be at the same record index");
            }
        }
    }

    #[test]
    fn split_paired_detects_mate_count_mismatch() {
        let td = tempfile::TempDir::new().unwrap();
        let p1 = td.path().join("r1.fastq");
        let p2 = td.path().join("r2.fastq");
        let r1: String = (0..10).map(|i| format!("@r{i}\nACGT\n+\nIIII\n")).collect();
        let r2: String = (0..9).map(|i| format!("@r{i}\nACGT\n+\nIIII\n")).collect(); // one fewer
        std::fs::write(&p1, &r1).unwrap();
        std::fs::write(&p2, &r2).unwrap();
        let err = split_paired(&p1, &p2, 2).unwrap_err();
        assert!(err.downcast_ref::<MateCountMismatch>().is_some(), "want MateCountMismatch, got {err}");
    }

    #[test]
    fn split_paired_too_few_records() {
        let td = tempfile::TempDir::new().unwrap();
        let p1 = td.path().join("r1.fastq");
        let p2 = td.path().join("r2.fastq");
        let r1: String = (0..3).map(|i| format!("@r{i}\nACGT\n+\nIIII\n")).collect();
        std::fs::write(&p1, &r1).unwrap();
        std::fs::write(&p2, &r1).unwrap();
        let err = split_paired(&p1, &p2, 4).unwrap_err();
        assert!(err.downcast_ref::<SplitTooFewRecords>().is_some(), "want SplitTooFewRecords, got {err}");
    }
}
