use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// A single FASTQ record with byte-oriented fields.
///
/// Fields are stored as `Vec<u8>` rather than `String` because FASTQ data is
/// ASCII and most operations work on raw bytes (compression, hashing, output).
/// This avoids unnecessary UTF-8 validation when constructing records from
/// decompressed byte streams.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    pub id: Vec<u8>,
    pub sequence: Vec<u8>,
    pub quality: Option<Vec<u8>>,
}

impl FastqRecord {
    /// Create a new FASTQ record
    pub fn new(id: Vec<u8>, sequence: Vec<u8>, quality: Option<Vec<u8>>) -> Self {
        Self {
            id,
            sequence,
            quality,
        }
    }
}

/// Fast FASTQ reader with buffering and optional gzip support
pub struct FastqReader<R: BufRead> {
    reader: R,
    is_fasta: bool,
    buffer: Vec<u8>,
    bom_checked: bool,
}

// Enum to hold either a plain file reader, gzipped reader, stdin reader, or a plain
// file bounded to a byte range (the NUMA byte-range compress worker).
pub enum FileReader {
    Plain(BufReader<std::fs::File>),
    Gzipped(BufReader<MultiGzDecoder<BufReader<std::fs::File>>>),
    Stdin(BufReader<std::io::Stdin>),
    StdinGzipped(BufReader<MultiGzDecoder<BufReader<std::io::Stdin>>>),
    /// A plain file seeked to a start offset and limited to a byte length — yields the
    /// records in `[start, end)` (the range is record-aligned by the splitter).
    Bounded(BufReader<std::io::Take<std::fs::File>>),
}

impl Read for FileReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            FileReader::Plain(r) => r.read(buf),
            FileReader::Gzipped(r) => r.read(buf),
            FileReader::Stdin(r) => r.read(buf),
            FileReader::StdinGzipped(r) => r.read(buf),
            FileReader::Bounded(r) => r.read(buf),
        }
    }
}

impl BufRead for FileReader {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self {
            FileReader::Plain(r) => r.fill_buf(),
            FileReader::Gzipped(r) => r.fill_buf(),
            FileReader::Stdin(r) => r.fill_buf(),
            FileReader::StdinGzipped(r) => r.fill_buf(),
            FileReader::Bounded(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            FileReader::Plain(r) => r.consume(amt),
            FileReader::Gzipped(r) => r.consume(amt),
            FileReader::Stdin(r) => r.consume(amt),
            FileReader::StdinGzipped(r) => r.consume(amt),
            FileReader::Bounded(r) => r.consume(amt),
        }
    }
}

impl FastqReader<FileReader> {
    /// Open a FASTQ file (auto-detects gzip), or read from stdin if path is `-`.
    pub fn from_path_or_stdin(path: impl AsRef<Path>, is_fasta: bool) -> Result<Self> {
        if crate::cli::is_stdio_path(path.as_ref()) {
            return Self::from_stdin(is_fasta);
        }
        Self::from_path(path, is_fasta)
    }

    /// Open a FASTQ file (auto-detects gzip)
    pub fn from_path(path: impl AsRef<Path>, is_fasta: bool) -> Result<Self> {
        let file = std::fs::File::open(path.as_ref())
            .with_context(|| format!("Failed to open file: {}", path.as_ref().display()))?;

        // Check if file is gzipped by reading magic bytes
        let mut buffered = BufReader::with_capacity(4 * 1024 * 1024, file);
        let is_gzipped = {
            let peek = buffered.fill_buf()?;
            peek.len() >= 2 && peek[0] == 0x1f && peek[1] == 0x8b
        };

        let reader = if is_gzipped {
            let decoder = MultiGzDecoder::new(buffered);
            FileReader::Gzipped(BufReader::new(decoder))
        } else {
            FileReader::Plain(buffered)
        };

        Ok(Self::new(reader, is_fasta))
    }

    /// Open a PLAIN (non-gzipped) FASTQ file bounded to the byte range `[start, end)`.
    ///
    /// The range MUST be record-aligned (the NUMA byte-range splitter guarantees this),
    /// so the reader yields whole records only. Used by the byte-range reference compress
    /// worker — the parallel to single-end/paired's `file.take()` bounded readers. Gzip is
    /// not supported here (gzipped reference sharding is rejected upstream).
    pub fn from_range(
        path: impl AsRef<Path>,
        start: u64,
        end: u64,
        is_fasta: bool,
    ) -> Result<Self> {
        use std::io::Seek;
        let mut file = std::fs::File::open(path.as_ref())
            .with_context(|| format!("Failed to open file: {}", path.as_ref().display()))?;
        file.seek(std::io::SeekFrom::Start(start))
            .with_context(|| format!("seek to {start} in {}", path.as_ref().display()))?;
        let bounded = file.take(end.saturating_sub(start));
        let reader = FileReader::Bounded(BufReader::with_capacity(4 * 1024 * 1024, bounded));
        Ok(Self::new(reader, is_fasta))
    }

    /// Read FASTQ from stdin (auto-detects gzip)
    pub fn from_stdin(is_fasta: bool) -> Result<Self> {
        let mut buffered = BufReader::with_capacity(4 * 1024 * 1024, std::io::stdin());
        let is_gzipped = {
            let peek = buffered.fill_buf()?;
            peek.len() >= 2 && peek[0] == 0x1f && peek[1] == 0x8b
        };

        let reader = if is_gzipped {
            let decoder = MultiGzDecoder::new(buffered);
            FileReader::StdinGzipped(BufReader::new(decoder))
        } else {
            FileReader::Stdin(buffered)
        };

        Ok(Self::new(reader, is_fasta))
    }
}

impl<R: BufRead> FastqReader<R> {
    /// Create a new FASTQ reader from any BufRead type
    pub fn new(reader: R, is_fasta: bool) -> Self {
        Self {
            reader,
            is_fasta,
            buffer: Vec::with_capacity(512), // Pre-allocate for typical read lengths
            bom_checked: false,
        }
    }

    /// Trim trailing \n and \r\n from the buffer in-place, return the trimmed length.
    #[inline]
    fn trim_newline(buf: &mut Vec<u8>) -> usize {
        while buf.last().is_some_and(|&b| b == b'\n' || b == b'\r') {
            buf.pop();
        }
        buf.len()
    }

    /// Read the next FASTQ record.
    ///
    /// Not an `Iterator::next` (it returns `Result<Option<_>>` for fallible,
    /// lazy reads), so the trait-confusion lint doesn't apply.
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Result<Option<FastqRecord>> {
        // Read ID line (read_until avoids UTF-8 validation overhead of read_line).
        // Tolerate blank lines at record boundaries — a trailing newline at EOF
        // (`...\n\n`) or stray blank lines between records — by skipping empty
        // lines until a real record header or true EOF. A valid FASTQ/FASTA header
        // is never empty (it starts with '@'/'>'), so an empty line is
        // unambiguously not a record. This is lenient like seqkit/BioPython and
        // matches how the FASTA continuation path already absorbs a trailing blank
        // line; without it a single trailing newline aborted the whole compress.
        let id = loop {
            self.buffer.clear();
            let bytes_read = self.reader.read_until(b'\n', &mut self.buffer)?;
            if bytes_read == 0 {
                return Ok(None); // EOF
            }
            // Strip a UTF-8 BOM (EF BB BF) from the start of the very first line.
            // Some editors and Windows tools add it; without stripping it lands in
            // the first record's ID and the next record's '+' check then fails
            // with a misleading "expected '+' separator" message.
            if !self.bom_checked {
                self.bom_checked = true;
                if self.buffer.starts_with(&[0xEF, 0xBB, 0xBF]) {
                    self.buffer.drain(..3);
                }
            }
            Self::trim_newline(&mut self.buffer);
            if !self.buffer.is_empty() {
                break self.buffer.clone();
            }
            // Blank line at a record boundary — skip and read the next line.
        };

        // Read sequence line
        self.buffer.clear();
        self.reader
            .read_until(b'\n', &mut self.buffer)
            .context("Invalid FASTQ: missing sequence line")?;
        Self::trim_newline(&mut self.buffer);
        let sequence = self.buffer.clone();

        if self.is_fasta {
            // Multi-line FASTA: a sequence may be wrapped across several lines.
            // Keep appending continuation lines until the next record header
            // ('>') or EOF. We peek the first byte of the next line via the
            // buffered reader rather than consuming it.
            let mut sequence = sequence;
            loop {
                let stop = {
                    let peek = self
                        .reader
                        .fill_buf()
                        .context("Invalid FASTA: error reading sequence continuation")?;
                    peek.is_empty() || peek[0] == b'>'
                };
                if stop {
                    break;
                }
                self.buffer.clear();
                self.reader.read_until(b'\n', &mut self.buffer)?;
                Self::trim_newline(&mut self.buffer);
                sequence.extend_from_slice(&self.buffer);
            }
            return Ok(Some(FastqRecord::new(id, sequence, None)));
        }

        // Read and validate separator line ('+')
        self.buffer.clear();
        let sep_bytes = self
            .reader
            .read_until(b'\n', &mut self.buffer)
            .context("Invalid FASTQ: missing '+' line")?;
        if sep_bytes == 0 {
            anyhow::bail!("Invalid FASTQ: unexpected EOF at '+' separator line");
        }
        if self.buffer.first() != Some(&b'+') {
            Self::trim_newline(&mut self.buffer);
            anyhow::bail!(
                "Invalid FASTQ: expected '+' separator line, got {:?}",
                String::from_utf8_lossy(&self.buffer)
            );
        }

        // Read quality line
        self.buffer.clear();
        let qual_bytes = self
            .reader
            .read_until(b'\n', &mut self.buffer)
            .context("Invalid FASTQ: missing quality line")?;
        if qual_bytes == 0 {
            anyhow::bail!("Invalid FASTQ: unexpected EOF at quality line");
        }
        Self::trim_newline(&mut self.buffer);
        let quality = self.buffer.clone();

        // Validate sequence and quality lengths match
        if quality.len() != sequence.len() {
            anyhow::bail!(
                "Invalid FASTQ: sequence length ({}) != quality length ({}) for read {}",
                sequence.len(),
                quality.len(),
                String::from_utf8_lossy(&id)
            );
        }

        Ok(Some(FastqRecord::new(id, sequence, Some(quality))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_fastq_parsing() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);

        let record1 = reader.next().unwrap().unwrap();
        assert_eq!(record1.id, b"@read1");
        assert_eq!(record1.sequence, b"ACGT");
        assert_eq!(record1.quality, Some(b"IIII".to_vec()));

        let record2 = reader.next().unwrap().unwrap();
        assert_eq!(record2.id, b"@read2");
        assert!(reader.next().unwrap().is_none());
    }

    #[test]
    fn test_trailing_blank_lines_tolerated() {
        // A trailing blank line (double newline at EOF) and stray blank lines
        // between records must NOT abort parsing — every real record is returned.
        let data = b"@read1\nACGT\n+\nIIII\n\n@read2\nTGCA\n+\nJJJJ\n\n\n";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);

        let r1 = reader.next().unwrap().unwrap();
        assert_eq!(r1.id, b"@read1");
        assert_eq!(r1.sequence, b"ACGT");
        let r2 = reader.next().unwrap().unwrap();
        assert_eq!(r2.id, b"@read2");
        assert_eq!(r2.sequence, b"TGCA");
        assert!(reader.next().unwrap().is_none());
    }

    #[test]
    fn test_invalid_separator() {
        let data = b"@read1\nACGT\nBAD_LINE\nIIII\n";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);
        let result = reader.next();
        assert!(
            result.is_err(),
            "Should reject FASTQ with invalid '+' separator"
        );
    }

    #[test]
    fn test_length_mismatch() {
        let data = b"@read1\nACGT\n+\nIII\n"; // quality shorter than sequence
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);
        let result = reader.next();
        assert!(
            result.is_err(),
            "Should reject FASTQ with mismatched lengths"
        );
    }

    #[test]
    fn test_empty_file() {
        let data = b"";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);
        assert!(reader.next().unwrap().is_none());
    }

    #[test]
    fn test_fasta_parsing() {
        let data = b">seq1\nACGT\n>seq2\nTGCA\n";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, true);

        let record1 = reader.next().unwrap().unwrap();
        assert_eq!(record1.id, b">seq1");
        assert_eq!(record1.sequence, b"ACGT");
        assert!(record1.quality.is_none());
    }

    #[test]
    fn test_fasta_multiline_sequence() {
        // Standard wrapped FASTA: a sequence split across multiple lines must be
        // concatenated into one record, and a `>` line starts the next record.
        let data = b">seq1\nACGT\nGGCC\nTT\n>seq2\nAAAA\nCCCC\n";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, true);

        let r1 = reader.next().unwrap().unwrap();
        assert_eq!(r1.id, b">seq1");
        assert_eq!(r1.sequence, b"ACGTGGCCTT");
        assert!(r1.quality.is_none());

        let r2 = reader.next().unwrap().unwrap();
        assert_eq!(r2.id, b">seq2");
        assert_eq!(r2.sequence, b"AAAACCCC");

        assert!(reader.next().unwrap().is_none());
    }

    #[test]
    fn test_fasta_multiline_last_record_no_trailing_newline() {
        let data = b">s\nAC\nGT"; // no trailing newline
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, true);
        let r = reader.next().unwrap().unwrap();
        assert_eq!(r.id, b">s");
        assert_eq!(r.sequence, b"ACGT");
        assert!(reader.next().unwrap().is_none());
    }

    #[test]
    fn test_crlf_line_endings() {
        // Windows-style CRLF must not appear in the stored bytes.
        let data = b"@r1\r\nACGT\r\n+\r\nIIII\r\n@r2\r\nTGCA\r\n+\r\nJJJJ\r\n";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);
        let r1 = reader.next().unwrap().unwrap();
        assert_eq!(r1.id, b"@r1");
        assert_eq!(r1.sequence, b"ACGT");
        assert_eq!(r1.quality, Some(b"IIII".to_vec()));
        let r2 = reader.next().unwrap().unwrap();
        assert_eq!(r2.id, b"@r2");
        assert_eq!(r2.sequence, b"TGCA");
        assert_eq!(r2.quality, Some(b"JJJJ".to_vec()));
        assert!(reader.next().unwrap().is_none());
    }

    #[test]
    fn test_utf8_bom_stripped() {
        // EF BB BF prepended; without stripping, the first ID would start with
        // garbage bytes and the next record's '+' check would fail confusingly.
        let mut data = vec![0xEF, 0xBB, 0xBF];
        data.extend_from_slice(b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n");
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);
        let r1 = reader.next().unwrap().unwrap();
        assert_eq!(r1.id, b"@r1");
        let r2 = reader.next().unwrap().unwrap();
        assert_eq!(r2.id, b"@r2");
    }

    #[test]
    fn test_truncated_at_quality_line() {
        // EOF in the middle of a record must error, not silently emit a partial record.
        let data = b"@r1\nACGT\n+\n";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);
        assert!(reader.next().is_err(), "truncated quality line must error");
    }

    #[test]
    fn test_truncated_at_separator() {
        let data = b"@r1\nACGT\n";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);
        assert!(
            reader.next().is_err(),
            "truncated separator line must error"
        );
    }

    #[test]
    fn reads_all_records_from_multi_member_gzip() {
        use flate2::{Compression, write::GzEncoder};
        use std::io::Write;
        let recs: [&[u8]; 2] = [b"@r1\nACGT\n+\nIIII\n", b"@r2\nTTTT\n+\nIIII\n"];
        let mut blob = Vec::new();
        for rec in recs {
            let mut e = GzEncoder::new(Vec::new(), Compression::default());
            e.write_all(rec).unwrap();
            blob.extend_from_slice(&e.finish().unwrap());
        }
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("multi.fastq.gz");
        std::fs::write(&path, &blob).unwrap();
        let mut reader = FastqReader::from_path(&path, false).unwrap();
        let mut ids = Vec::new();
        while let Some(rec) = reader.next().unwrap() {
            ids.push(rec.id.clone());
        }
        assert_eq!(ids, vec![b"@r1".to_vec(), b"@r2".to_vec()]);
    }
}
