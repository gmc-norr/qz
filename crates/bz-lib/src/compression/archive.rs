use crate::compression::streams::NUM_STREAMS;
use anyhow::Result;
use std::io::{Read, Write};

/// Directory to place an output's sibling temp file in: the output's parent, or
/// the current directory when the path has no parent component. Used so atomic
/// temp files land on the same filesystem as the final output (rename must be
/// intra-filesystem).
pub(crate) fn output_parent(path: &std::path::Path) -> &std::path::Path {
    match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => std::path::Path::new("."),
    }
}

pub const ARCHIVE_MAGIC: [u8; 2] = *b"BZ";
// v7: fqz v2 context (drop base, full prev_q2).
// v8: MD:Z derivation — MD stripped from aux, regenerated on decode from the
// consensus + an exception map (3 new streams).
// v9: MC:Z + TLEN mate-derivation (3 more streams; tlen stream holds only
// non-derivable values). Older archives aren't decodable — re-compress.
// v10: PNEXT/RNEXT mate-derivation (1 more stream, next_derivable bitmap;
// next_pos/next_ref_id streams hold only non-derivable values).
// v11: NM:i derivation (3 more streams: nm_present bitmap, nm_index, nm_type);
// NM = MD-substitutions + CIGAR(I+D), regenerated from the recovered reference.
// v12: quality codec switched to htscodecs fqzcomp on raw Phred, with
// reverse-strand reads reoriented to sequencing-cycle order (reversed back on
// decode via the FLAG stream). ~-3% smaller archives on full-resolution WGS.
// v13: zstd backend removed — the two per-archive compressor-selector bytes
// (alignment_compressor/aux_compressor) were dropped from the header, so the
// v13 body layout is incompatible with v12. The reader rejects any version !=
// ARCHIVE_VERSION, so old v12 archives fail with a clear version error rather
// than silently misparsing the shorter header (re-compress them with current bz).
pub const ARCHIVE_VERSION: u8 = 13;

/// Flags for the archive header.
pub const FLAG_CONSENSUS_DELTA: u8 = 0x01;
/// Archive-level flag: quality was lossily reduced (variant-aware reduction).
/// Informational — decode emits the stored (reduced) quality regardless.
pub const FLAG_QUALITY_LOSSY: u8 = 0x02;

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

        // SAM header compressed (with size limit to prevent OOM on malicious input)
        const MAX_COMPRESSED_HEADER: usize = 100 * 1024 * 1024; // 100 MB
        r.read_exact(&mut buf4)?;
        let sam_len = u32::from_le_bytes(buf4) as usize;
        if sam_len > MAX_COMPRESSED_HEADER {
            anyhow::bail!(
                "Compressed SAM header size {} exceeds limit ({} MB)",
                sam_len,
                MAX_COMPRESSED_HEADER / (1024 * 1024)
            );
        }
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

/// Per-chunk flag: quality stream uses fqz encoding.
pub const CHUNK_FLAG_FQZ_QUALITY: u8 = 0x01;

/// Per-chunk header (written before each chunk's stream data).
pub struct ChunkHeader {
    pub num_records: u32,
    /// Per-chunk flags (bit 0 = fqz used for this chunk's quality stream).
    pub chunk_flags: u8,
    /// CRC32 (IEEE) over the concatenated compressed stream payloads for this chunk.
    /// Computed at write time; verified before BSC decompression at read time.
    /// Detects on-disk corruption of the compressed bytes.
    pub crc32: u32,
    /// CRC32 (IEEE) over the original BAM record bytes (each record's raw bytes,
    /// concatenated in order) for this chunk. Computed at compress time from the
    /// source records; recomputed at decompress/verify time over the
    /// *reconstructed* records and compared. Detects encode/decode asymmetries
    /// — i.e. proves the chunk round-trips, which the compressed-payload CRC
    /// above does not.
    pub content_crc32: u32,
    pub stream_sizes: [u32; NUM_STREAMS],
}

impl ChunkHeader {
    pub fn write_to<W: Write>(&self, w: &mut W) -> Result<()> {
        w.write_all(&self.num_records.to_le_bytes())?;
        w.write_all(&[self.chunk_flags])?;
        w.write_all(&self.crc32.to_le_bytes())?;
        w.write_all(&self.content_crc32.to_le_bytes())?;
        for &size in &self.stream_sizes {
            w.write_all(&size.to_le_bytes())?;
        }
        Ok(())
    }

    pub fn read_from<R: Read>(r: &mut R) -> Result<Self> {
        let mut buf4 = [0u8; 4];
        r.read_exact(&mut buf4)?;
        let num_records = u32::from_le_bytes(buf4);

        let mut flags_buf = [0u8; 1];
        r.read_exact(&mut flags_buf)?;
        let chunk_flags = flags_buf[0];

        r.read_exact(&mut buf4)?;
        let crc32 = u32::from_le_bytes(buf4);

        r.read_exact(&mut buf4)?;
        let content_crc32 = u32::from_le_bytes(buf4);

        let mut stream_sizes = [0u32; NUM_STREAMS];
        for size in &mut stream_sizes {
            r.read_exact(&mut buf4)?;
            *size = u32::from_le_bytes(buf4);
        }

        Ok(Self {
            num_records,
            chunk_flags,
            crc32,
            content_crc32,
            stream_sizes,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn current_archive_version_is_13() {
        // Bumped from 12 when the zstd backend (and its two header bytes) was
        // removed; the v13 header layout is incompatible with v12.
        assert_eq!(ARCHIVE_VERSION, 13);
    }

    #[test]
    fn read_from_rejects_old_archive_version() {
        // A v12 (pre-zstd-removal) archive must be rejected with a clear version
        // error rather than silently misparsing the shorter v13 header.
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&ARCHIVE_MAGIC);
        bytes.push(12); // old version
        let mut cur = std::io::Cursor::new(bytes);
        let err = ArchiveHeader::read_from(&mut cur).map(|_| ()).unwrap_err();
        assert!(
            err.to_string()
                .contains("Unsupported BZ archive version 12"),
            "expected version-rejection error, got: {err}"
        );
    }
}
