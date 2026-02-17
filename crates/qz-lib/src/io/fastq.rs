use anyhow::{Context, Result};
use flate2::read::GzDecoder;
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
        Self { id, sequence, quality }
    }
}

/// Fast FASTQ reader with buffering and optional gzip support
pub struct FastqReader<R: BufRead> {
    reader: R,
    is_fasta: bool,
    buffer: Vec<u8>,
}

// Enum to hold either a plain file reader, gzipped reader, or stdin reader
pub enum FileReader {
    Plain(BufReader<std::fs::File>),
    Gzipped(BufReader<GzDecoder<BufReader<std::fs::File>>>),
    Stdin(BufReader<std::io::Stdin>),
    StdinGzipped(BufReader<GzDecoder<BufReader<std::io::Stdin>>>),
}

impl Read for FileReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            FileReader::Plain(r) => r.read(buf),
            FileReader::Gzipped(r) => r.read(buf),
            FileReader::Stdin(r) => r.read(buf),
            FileReader::StdinGzipped(r) => r.read(buf),
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
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            FileReader::Plain(r) => r.consume(amt),
            FileReader::Gzipped(r) => r.consume(amt),
            FileReader::Stdin(r) => r.consume(amt),
            FileReader::StdinGzipped(r) => r.consume(amt),
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
            .with_context(|| format!("Failed to open file: {:?}", path.as_ref()))?;

        // Check if file is gzipped by reading magic bytes
        let mut buffered = BufReader::with_capacity(4 * 1024 * 1024, file);
        let is_gzipped = {
            let peek = buffered.fill_buf()?;
            peek.len() >= 2 && peek[0] == 0x1f && peek[1] == 0x8b
        };

        let reader = if is_gzipped {
            let decoder = GzDecoder::new(buffered);
            FileReader::Gzipped(BufReader::new(decoder))
        } else {
            FileReader::Plain(buffered)
        };

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
            let decoder = GzDecoder::new(buffered);
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

    /// Read the next FASTQ record
    pub fn next(&mut self) -> Result<Option<FastqRecord>> {
        // Read ID line (read_until avoids UTF-8 validation overhead of read_line)
        self.buffer.clear();
        let bytes_read = self.reader.read_until(b'\n', &mut self.buffer)?;
        if bytes_read == 0 {
            return Ok(None); // EOF
        }
        Self::trim_newline(&mut self.buffer);
        let id = self.buffer.clone();

        // Read sequence line
        self.buffer.clear();
        self.reader.read_until(b'\n', &mut self.buffer)
            .context("Invalid FASTQ: missing sequence line")?;
        Self::trim_newline(&mut self.buffer);
        let sequence = self.buffer.clone();

        if self.is_fasta {
            return Ok(Some(FastqRecord::new(id, sequence, None)));
        }

        // Read and validate separator line ('+')
        self.buffer.clear();
        let sep_bytes = self.reader.read_until(b'\n', &mut self.buffer)
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
        let qual_bytes = self.reader.read_until(b'\n', &mut self.buffer)
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
    fn test_invalid_separator() {
        let data = b"@read1\nACGT\nBAD_LINE\nIIII\n";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);
        let result = reader.next();
        assert!(result.is_err(), "Should reject FASTQ with invalid '+' separator");
    }

    #[test]
    fn test_length_mismatch() {
        let data = b"@read1\nACGT\n+\nIII\n"; // quality shorter than sequence
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);
        let result = reader.next();
        assert!(result.is_err(), "Should reject FASTQ with mismatched lengths");
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
}
