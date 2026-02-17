use crate::compression::streams::NUM_STREAMS;
use anyhow::Result;
use std::io::{Read, Write};

pub const ARCHIVE_MAGIC: [u8; 2] = *b"BZ";
pub const ARCHIVE_VERSION: u8 = 2;

/// Flags for the archive header.
pub const FLAG_CONSENSUS_DELTA: u8 = 0x01;
pub const FLAG_QUALITY_CTX: u8 = 0x02;

/// Archive header (global, written once).
pub struct ArchiveHeader {
    pub flags: u8,
    pub num_records: u64,
    pub num_chunks: u32,
    pub sam_header_compressed: Vec<u8>,
}

impl ArchiveHeader {
    pub fn write_to<W: Write>(&self, w: &mut W) -> Result<()> {
        // Magic + version + reserved
        w.write_all(&ARCHIVE_MAGIC)?;
        w.write_all(&[ARCHIVE_VERSION])?;
        w.write_all(&[0u8])?; // reserved

        // Flags
        w.write_all(&[self.flags])?;

        // Total record count
        w.write_all(&self.num_records.to_le_bytes())?;

        // Number of chunks
        w.write_all(&self.num_chunks.to_le_bytes())?;

        // SAM header (compressed)
        w.write_all(&(self.sam_header_compressed.len() as u32).to_le_bytes())?;
        w.write_all(&self.sam_header_compressed)?;

        Ok(())
    }

    pub fn read_from<R: Read>(r: &mut R) -> Result<Self> {
        // Magic
        let mut magic = [0u8; 2];
        r.read_exact(&mut magic)?;
        if magic != ARCHIVE_MAGIC {
            anyhow::bail!("Not a BZ archive (invalid magic: {:?})", magic);
        }

        // Version
        let mut version = [0u8; 1];
        r.read_exact(&mut version)?;
        if version[0] != ARCHIVE_VERSION {
            anyhow::bail!(
                "Unsupported BZ archive version {} (expected {})",
                version[0],
                ARCHIVE_VERSION
            );
        }

        // Reserved
        let mut reserved = [0u8; 1];
        r.read_exact(&mut reserved)?;

        // Flags
        let mut flags = [0u8; 1];
        r.read_exact(&mut flags)?;

        // Num records
        let mut buf8 = [0u8; 8];
        r.read_exact(&mut buf8)?;
        let num_records = u64::from_le_bytes(buf8);

        // Num chunks
        let mut buf4 = [0u8; 4];
        r.read_exact(&mut buf4)?;
        let num_chunks = u32::from_le_bytes(buf4);

        // SAM header compressed
        r.read_exact(&mut buf4)?;
        let sam_len = u32::from_le_bytes(buf4) as usize;
        let mut sam_header_compressed = vec![0u8; sam_len];
        if sam_len > 0 {
            r.read_exact(&mut sam_header_compressed)?;
        }

        Ok(Self {
            flags: flags[0],
            num_records,
            num_chunks,
            sam_header_compressed,
        })
    }
}

/// Per-chunk header (written before each chunk's stream data).
pub struct ChunkHeader {
    pub num_records: u32,
    pub stream_sizes: [u32; NUM_STREAMS],
}

impl ChunkHeader {
    pub fn write_to<W: Write>(&self, w: &mut W) -> Result<()> {
        w.write_all(&self.num_records.to_le_bytes())?;
        for &size in &self.stream_sizes {
            w.write_all(&size.to_le_bytes())?;
        }
        Ok(())
    }

    pub fn read_from<R: Read>(r: &mut R) -> Result<Self> {
        let mut buf4 = [0u8; 4];
        r.read_exact(&mut buf4)?;
        let num_records = u32::from_le_bytes(buf4);

        let mut stream_sizes = [0u32; NUM_STREAMS];
        for size in &mut stream_sizes {
            r.read_exact(&mut buf4)?;
            *size = u32::from_le_bytes(buf4);
        }

        Ok(Self {
            num_records,
            stream_sizes,
        })
    }
}
