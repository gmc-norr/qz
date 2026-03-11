use anyhow::Result;
use noodles::bgzf;
use std::io::{Read, Write};
use std::num::NonZeroUsize;

/// A raw BAM record stored as its uncompressed byte representation.
/// Fields are extracted by slicing into `data` at BAM-specified offsets.
///
/// **Invariant**: `data.len() >= 32` and all variable-length fields fit within
/// `data`. This is enforced at construction via [`RawBamRecord::new`]; field
/// accessors may therefore use direct indexing without re-checking bounds.
///
/// BAM record layout (all little-endian, excluding the 4-byte block_size prefix):
///   bytes 0..4:    refID (i32)
///   bytes 4..8:    pos (i32)
///   byte  8:       l_read_name (u8)
///   byte  9:       mapq (u8)
///   bytes 10..12:  bin (u16)
///   bytes 12..14:  n_cigar_op (u16)
///   bytes 14..16:  flag (u16)
///   bytes 16..20:  l_seq (i32)
///   bytes 20..24:  next_refID (i32)
///   bytes 24..28:  next_pos (i32)
///   bytes 28..32:  tlen (i32)
///   bytes 32..:    read_name, cigar, seq, qual, aux (variable)
pub struct RawBamRecord {
    data: Vec<u8>,
}

impl RawBamRecord {
    /// Construct a validated record. Returns an error if `data` is structurally
    /// invalid (too short, negative l_seq, or any field extends past the end).
    pub fn new(data: Vec<u8>) -> Result<Self> {
        let record = Self { data };
        record.validate()?;
        Ok(record)
    }

    /// Return the raw record bytes.
    pub fn as_bytes(&self) -> &[u8] {
        &self.data
    }

    fn validate(&self) -> Result<()> {
        let len = self.data.len();
        if len < 32 {
            anyhow::bail!("BAM record too short: {} bytes (minimum 32)", len);
        }
        let l_seq = self.l_seq();
        if l_seq < 0 {
            anyhow::bail!("BAM record has negative l_seq: {}", l_seq);
        }
        let l_seq = l_seq as usize;
        let read_name_end = 32usize
            .checked_add(self.l_read_name() as usize)
            .filter(|&v| v <= len)
            .ok_or_else(|| anyhow::anyhow!(
                "BAM record read_name extends past end (l_read_name={}, record_len={})",
                self.l_read_name(), len
            ))?;
        let cigar_end = read_name_end
            .checked_add(4 * self.n_cigar_op() as usize)
            .filter(|&v| v <= len)
            .ok_or_else(|| anyhow::anyhow!(
                "BAM record CIGAR extends past end (n_cigar_op={}, record_len={})",
                self.n_cigar_op(), len
            ))?;
        let seq_end = cigar_end
            .checked_add((l_seq + 1) / 2)
            .filter(|&v| v <= len)
            .ok_or_else(|| anyhow::anyhow!(
                "BAM record seq extends past end (l_seq={}, record_len={})",
                l_seq, len
            ))?;
        let _qual_end = seq_end
            .checked_add(l_seq)
            .filter(|&v| v <= len)
            .ok_or_else(|| anyhow::anyhow!(
                "BAM record qual extends past end (l_seq={}, record_len={})",
                l_seq, len
            ))?;
        Ok(())
    }

    pub fn ref_id(&self) -> i32 {
        i32::from_le_bytes(self.data[0..4].try_into().unwrap())
    }

    pub fn pos(&self) -> i32 {
        i32::from_le_bytes(self.data[4..8].try_into().unwrap())
    }

    pub fn l_read_name(&self) -> u8 {
        self.data[8]
    }

    pub fn mapq(&self) -> u8 {
        self.data[9]
    }

    pub fn bin(&self) -> u16 {
        u16::from_le_bytes(self.data[10..12].try_into().unwrap())
    }

    pub fn n_cigar_op(&self) -> u16 {
        u16::from_le_bytes(self.data[12..14].try_into().unwrap())
    }

    pub fn flag(&self) -> u16 {
        u16::from_le_bytes(self.data[14..16].try_into().unwrap())
    }

    pub fn l_seq(&self) -> i32 {
        i32::from_le_bytes(self.data[16..20].try_into().unwrap())
    }

    pub fn next_ref_id(&self) -> i32 {
        i32::from_le_bytes(self.data[20..24].try_into().unwrap())
    }

    pub fn next_pos(&self) -> i32 {
        i32::from_le_bytes(self.data[24..28].try_into().unwrap())
    }

    pub fn tlen(&self) -> i32 {
        i32::from_le_bytes(self.data[28..32].try_into().unwrap())
    }

    fn read_name_end(&self) -> usize {
        32 + self.l_read_name() as usize
    }

    fn cigar_end(&self) -> usize {
        self.read_name_end() + 4 * self.n_cigar_op() as usize
    }

    fn seq_end(&self) -> usize {
        let l_seq = self.l_seq() as usize;
        self.cigar_end() + ((l_seq + 1) / 2)
    }

    fn qual_end(&self) -> usize {
        self.seq_end() + self.l_seq() as usize
    }

    pub fn read_name(&self) -> &[u8] {
        &self.data[32..self.read_name_end()]
    }

    pub fn cigar_bytes(&self) -> &[u8] {
        &self.data[self.read_name_end()..self.cigar_end()]
    }

    /// Raw BAM-packed sequence bytes (2 nibbles per byte, high nibble first).
    pub fn seq_bytes(&self) -> &[u8] {
        &self.data[self.cigar_end()..self.seq_end()]
    }

    pub fn qual_bytes(&self) -> &[u8] {
        &self.data[self.seq_end()..self.qual_end()]
    }

    pub fn aux_bytes(&self) -> &[u8] {
        &self.data[self.qual_end()..]
    }

    /// Unpack BAM 4-bit encoded sequence into individual nibble values.
    pub fn unpack_seq_nibbles(&self) -> Vec<u8> {
        unpack_nibbles(self.seq_bytes(), self.l_seq() as usize)
    }

    /// Iterate CIGAR operations as (op_code, op_len) pairs.
    /// op_code: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
    pub fn cigar_ops(&self) -> Vec<(u8, u32)> {
        let bytes = self.cigar_bytes();
        let n = self.n_cigar_op() as usize;
        let mut ops = Vec::with_capacity(n);
        for i in 0..n {
            let val = u32::from_le_bytes(bytes[i * 4..(i + 1) * 4].try_into().unwrap());
            let op = (val & 0xF) as u8;
            let len = val >> 4;
            ops.push((op, len));
        }
        ops
    }
}

/// Unpack BAM 4-bit encoded sequence nibbles.
pub fn unpack_nibbles(packed: &[u8], l_seq: usize) -> Vec<u8> {
    let mut nibbles = vec![0u8; l_seq];
    let pairs = l_seq / 2;
    for i in 0..pairs {
        nibbles[i * 2] = packed[i] >> 4;
        nibbles[i * 2 + 1] = packed[i] & 0x0F;
    }
    if l_seq % 2 != 0 {
        nibbles[l_seq - 1] = packed[pairs] >> 4;
    }
    nibbles
}

/// Pack individual nibble values back into BAM 4-bit encoding (2 per byte, high first).
pub fn pack_nibbles(nibbles: &[u8]) -> Vec<u8> {
    let out_len = (nibbles.len() + 1) / 2;
    let mut packed = vec![0u8; out_len];
    let pairs = nibbles.len() / 2;
    for i in 0..pairs {
        packed[i] = (nibbles[i * 2] << 4) | nibbles[i * 2 + 1];
    }
    if nibbles.len() % 2 != 0 {
        packed[pairs] = nibbles[pairs * 2] << 4;
    }
    packed
}

/// Reads raw BAM records from a BGZF-compressed BAM file.
/// Generic over the BGZF reader type (single-threaded or multi-threaded).
pub struct RawBamReader<R: Read> {
    reader: R,
    /// Raw SAM header text bytes (for lossless roundtrip).
    pub header_raw: Vec<u8>,
    /// Raw reference dictionary bytes: [n_ref:4B] then per-ref [name_len:4B, name, seq_len:4B].
    pub ref_dict_bytes: Vec<u8>,
}

impl RawBamReader<bgzf::io::Reader<std::io::BufReader<std::fs::File>>> {
    pub fn from_path(path: &std::path::Path) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let buf = std::io::BufReader::with_capacity(4 * 1024 * 1024, file);
        Self::new(bgzf::io::Reader::new(buf))
            .map_err(|e| wrap_bam_open_error(e, path))
    }
}

impl RawBamReader<bgzf::io::MultithreadedReader<std::io::BufReader<std::fs::File>>> {
    pub fn from_path_mt(path: &std::path::Path, worker_count: NonZeroUsize) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let buf = std::io::BufReader::with_capacity(4 * 1024 * 1024, file);
        Self::new(bgzf::io::MultithreadedReader::with_worker_count(worker_count, buf))
            .map_err(|e| wrap_bam_open_error(e, path))
    }
}

/// Wrap low-level BGZF/IO errors with a user-friendly message.
fn wrap_bam_open_error(e: anyhow::Error, path: &std::path::Path) -> anyhow::Error {
    // If it's already our own descriptive error, pass through
    let msg = e.to_string();
    if msg.starts_with("Not a BAM file") || msg.starts_with("Invalid BAM") {
        return e;
    }
    anyhow::anyhow!(
        "Failed to read {:?} as BAM: {}. Is this a valid BAM file?",
        path, msg
    )
}

impl<R: Read> RawBamReader<R> {
    pub fn new(mut reader: R) -> Result<Self> {
        // Read BAM magic "BAM\1"
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != b"BAM\x01" {
            anyhow::bail!("Not a BAM file (invalid magic)");
        }

        // Read SAM header text
        let mut header_len_buf = [0u8; 4];
        reader.read_exact(&mut header_len_buf)?;
        let header_len_i32 = i32::from_le_bytes(header_len_buf);
        if header_len_i32 < 0 {
            anyhow::bail!("Invalid BAM header length: {}", header_len_i32);
        }
        let header_len = header_len_i32 as usize;
        let mut header_raw = vec![0u8; header_len];
        if header_len > 0 {
            reader.read_exact(&mut header_raw)?;
        }

        // Read reference dictionary
        let mut ref_dict_bytes = Vec::new();
        let mut n_ref_buf = [0u8; 4];
        reader.read_exact(&mut n_ref_buf)?;
        let n_ref_i32 = i32::from_le_bytes(n_ref_buf);
        if n_ref_i32 < 0 {
            anyhow::bail!("Invalid BAM reference count: {}", n_ref_i32);
        }
        let n_ref = n_ref_i32 as usize;
        ref_dict_bytes.extend_from_slice(&n_ref_buf);

        for _ in 0..n_ref {
            let mut name_len_buf = [0u8; 4];
            reader.read_exact(&mut name_len_buf)?;
            let name_len_i32 = i32::from_le_bytes(name_len_buf);
            if name_len_i32 < 0 || name_len_i32 > 256 {
                anyhow::bail!(
                    "Invalid BAM reference name length: {}",
                    name_len_i32
                );
            }
            let name_len = name_len_i32 as usize;
            let mut name = vec![0u8; name_len];
            reader.read_exact(&mut name)?;
            let mut seq_len_buf = [0u8; 4];
            reader.read_exact(&mut seq_len_buf)?;
            ref_dict_bytes.extend_from_slice(&name_len_buf);
            ref_dict_bytes.extend_from_slice(&name);
            ref_dict_bytes.extend_from_slice(&seq_len_buf);
        }

        Ok(Self {
            reader,
            header_raw,
            ref_dict_bytes,
        })
    }

    /// Read the next BAM record as raw bytes. Returns None at EOF.
    pub fn next_record(&mut self) -> Result<Option<RawBamRecord>> {
        let mut block_size_buf = [0u8; 4];
        match self.reader.read_exact(&mut block_size_buf) {
            Ok(()) => {}
            Err(ref e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(anyhow::Error::from(e)),
        }
        let block_size_i32 = i32::from_le_bytes(block_size_buf);
        if block_size_i32 < 32 {
            anyhow::bail!(
                "Invalid BAM record: block_size {} (must be >= 32)",
                block_size_i32
            );
        }
        let block_size = block_size_i32 as usize;
        let mut data = vec![0u8; block_size];
        self.reader.read_exact(&mut data)?;
        Ok(Some(RawBamRecord::new(data)?))
    }

    /// Read up to `limit` records into a Vec.
    pub fn read_chunk(&mut self, limit: usize) -> Result<Vec<RawBamRecord>> {
        let mut records = Vec::with_capacity(limit.min(65536));
        for _ in 0..limit {
            match self.next_record()? {
                Some(rec) => records.push(rec),
                None => break,
            }
        }
        Ok(records)
    }
}

/// Writes BAM records to a BGZF-compressed file.
/// Generic over the BGZF writer type (single-threaded or multi-threaded).
pub struct BamWriter<W: Write> {
    writer: W,
}

impl BamWriter<bgzf::io::Writer<std::io::BufWriter<std::fs::File>>> {
    pub fn from_path(path: &std::path::Path) -> Result<Self> {
        let file = std::fs::File::create(path)?;
        let buf = std::io::BufWriter::with_capacity(4 * 1024 * 1024, file);
        Ok(Self { writer: bgzf::io::Writer::new(buf) })
    }

    pub fn finish(self) -> Result<()> {
        self.writer.finish()?;
        Ok(())
    }
}

impl BamWriter<bgzf::io::MultithreadedWriter<std::io::BufWriter<std::fs::File>>> {
    pub fn from_path_mt(path: &std::path::Path, worker_count: NonZeroUsize) -> Result<Self> {
        let file = std::fs::File::create(path)?;
        let buf = std::io::BufWriter::with_capacity(4 * 1024 * 1024, file);
        Ok(Self {
            writer: bgzf::io::MultithreadedWriter::with_worker_count(worker_count, buf),
        })
    }

    pub fn finish(mut self) -> Result<()> {
        self.writer.finish()?;
        Ok(())
    }
}

impl<W: Write> BamWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Write the BAM header (magic + SAM header text + reference dictionary).
    pub fn write_header(&mut self, header_raw: &[u8], ref_dict_bytes: &[u8]) -> Result<()> {
        self.writer.write_all(b"BAM\x01")?;
        self.writer
            .write_all(&(header_raw.len() as i32).to_le_bytes())?;
        self.writer.write_all(header_raw)?;
        self.writer.write_all(ref_dict_bytes)?;
        Ok(())
    }

    /// Write a single BAM record from its raw bytes (as stored in RawBamRecord.data).
    pub fn write_record(&mut self, data: &[u8]) -> Result<()> {
        self.writer
            .write_all(&(data.len() as i32).to_le_bytes())?;
        self.writer.write_all(data)?;
        Ok(())
    }
}
