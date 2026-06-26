//! Chunk-major (archive v5) directory footer.
//!
//! A v5 archive is: `[v5 front header][interleaved per-chunk blocks][directory
//! footer][20-byte locator]`. Each block uses the per-block frame
//! (`[block_len u32][record_count u32][crc32 u32][payload]`). The directory maps
//! every block to its stream (role) and absolute frame offset, so the decoder
//! rebuilds the same per-stream `StreamIndex` it would from a contiguous walk.
//!
//! The locator lives at the very end (`footer_len u64 | footer_crc32 u32 |
//! magic "QZFOOTR1"`), so the encoder never seeks the output — it works to a
//! pipe/stdout. Single-end, paired-end, and reference archives all use this
//! container (`archive_type` in the front header distinguishes them: 0=single,
//! 1=paired, 2=reference); none back-patch a footer offset, so all support
//! seek-free output.

use anyhow::{Result, bail};

/// Which stream a directory entry belongs to. Stable on-disk byte values.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub(crate) enum StreamRole {
    Headers = 0,
    Sequence = 1,
    /// Reserved, never emitted: the `n_mask` module was removed and no encoder
    /// writes a Nmask frame. Byte value 2 stays reserved; decode rejects any
    /// Nmask entry. (Reference mode uses a separate `NBitmap` role for N runs.)
    Nmask = 2,
    Qual = 3,
    RcFlags = 4,
    HeaderDelta = 5,
    PackedBacking = 6,
    IntervalMap = 7,
    NBitmap = 8,
    ReferenceMeta = 9,
    MappedFlags = 10,
    Positions = 11,
    Strands = 12,
    ReadLen = 13,
    EditCount = 14,
    EditPos = 15,
    EditBase = 16,
    FallbackPool = 17,
    /// Whole-archive global block (`chunk_index` = `GLOBAL_SENTINEL`) carrying a
    /// row-major table of per-chunk decoded output byte lengths. Used by the NUMA
    /// decode path to pre-partition work without decoding blocks.
    ChunkDecodedSizes = 18,
    /// Cluster mode: the reordered, RC-canonicalized sequence stream for a chunk,
    /// compressed with `CODEC_ZSTD_LONG`.
    ClusteredSeq = 19,
    /// Cluster mode: per-record reverse-complement strand bits for a chunk (1 bit
    /// per read, packed) so decode can un-canonicalize.
    StrandBits = 20,
    /// Cluster mode: whole-archive global block (`chunk_index` = `GLOBAL_SENTINEL`)
    /// carrying cluster metadata (total read count, etc.).
    ClusterMeta = 21,
}

impl StreamRole {
    pub fn from_u8(v: u8) -> Option<StreamRole> {
        Some(match v {
            0 => StreamRole::Headers,
            1 => StreamRole::Sequence,
            2 => StreamRole::Nmask,
            3 => StreamRole::Qual,
            4 => StreamRole::RcFlags,
            5 => StreamRole::HeaderDelta,
            6 => StreamRole::PackedBacking,
            7 => StreamRole::IntervalMap,
            8 => StreamRole::NBitmap,
            9 => StreamRole::ReferenceMeta,
            10 => StreamRole::MappedFlags,
            11 => StreamRole::Positions,
            12 => StreamRole::Strands,
            13 => StreamRole::ReadLen,
            14 => StreamRole::EditCount,
            15 => StreamRole::EditPos,
            16 => StreamRole::EditBase,
            17 => StreamRole::FallbackPool,
            18 => StreamRole::ChunkDecodedSizes,
            19 => StreamRole::ClusteredSeq,
            20 => StreamRole::StrandBits,
            21 => StreamRole::ClusterMeta,
            _ => return None,
        })
    }
}

/// One block's directory entry. `offset`/`length` describe the **full block frame**
/// (the 12-byte header included), so a span check is `[offset, offset+length)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ChunkDirEntry {
    pub chunk_index: u32,
    pub role: StreamRole,
    /// `0` = single-end / not-applicable, `1` = R1, `2` = R2.
    pub mate: u8,
    /// Per-entry codec id, in the **`codec_ids`** namespace (`CODEC_BSC` = 1,
    /// `CODEC_OPENZL` = 2, `CODEC_COLUMNAR` = 3, `CODEC_QUALITY_CTX` = 4,
    /// `CODEC_PACKED_CONSENSUS` = 5, `CODEC_FQZCOMP` = 6). This is **one namespace for
    /// every `archive_type`** — single-end, paired, and reference all stamp it, so a
    /// given byte means the same codec regardless of mode. The single-end
    /// *front-header* codec bytes now use this SAME namespace (identical values, e.g.
    /// Fqzcomp = 6 in both); single-end decode reads its codecs from the front
    /// header, and this directory byte is advisory there — uniformity-checked per
    /// role, never dispatched on.
    pub codec: u8,
    pub offset: u64,
    pub length: u64,
    pub record_count: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ChunkDirectory {
    pub num_reads: u64,
    pub num_chunks: u32,
    pub entries: Vec<ChunkDirEntry>,
}

/// 8-byte trailing magic that ends every v5 archive.
pub(crate) const FOOTER_MAGIC: [u8; 8] = *b"QZFOOTR1";
/// Locator size: footer_len(8) + footer_crc32(4) + magic(8).
pub(crate) const LOCATOR_LEN: usize = 20;
/// Sentinel `chunk_index` marking a whole-archive (global) entry (reference mode).
pub(crate) const GLOBAL_SENTINEL: u32 = u32::MAX;
/// On-disk size of one serialized `ChunkDirEntry`:
/// chunk_index u32(4) + role u8(1) + mate u8(1) + codec u8(1) + offset u64(8)
/// + length u64(8) + record_count u64(8) = 31 bytes.
pub(crate) const ENTRY_BYTES: usize = 4 + 1 + 1 + 1 + 8 + 8 + 8;
/// Footer body fixed prefix: num_reads u64(8) + num_chunks u32(4) + n_entries u32(4).
const FOOTER_PREFIX: usize = 8 + 4 + 4;

/// Serialize a per-chunk decoded-size table into a `ChunkDecodedSizes` block payload.
///
/// Layout: `[n_mates: u8][num_chunks: u32 LE][sizes: u64 LE; num_chunks*n_mates]`,
/// row-major: `sizes[chunk * n_mates + mate_idx]`. Raw/uncompressed (the table is tiny).
pub(crate) fn encode_decoded_sizes(n_mates: u8, num_chunks: u32, sizes: &[u64]) -> Vec<u8> {
    debug_assert_eq!(sizes.len(), num_chunks as usize * n_mates as usize);
    let mut out = Vec::with_capacity(5 + sizes.len() * 8);
    out.push(n_mates);
    out.extend_from_slice(&num_chunks.to_le_bytes());
    for s in sizes {
        out.extend_from_slice(&s.to_le_bytes());
    }
    out
}

/// Parse a `ChunkDecodedSizes` block payload. Rejects bad `n_mates`, length
/// mismatch, and overflow. Returns `(n_mates, num_chunks, sizes)`.
pub(crate) fn parse_decoded_sizes(payload: &[u8]) -> Result<(u8, u32, Vec<u64>)> {
    if payload.len() < 5 {
        bail!("ChunkDecodedSizes payload too short ({})", payload.len());
    }
    let n_mates = payload[0];
    if n_mates == 0 || n_mates > 2 {
        bail!("ChunkDecodedSizes: bad n_mates {n_mates}");
    }
    let num_chunks = u32::from_le_bytes(payload[1..5].try_into().unwrap());
    let count = (num_chunks as usize)
        .checked_mul(n_mates as usize)
        .ok_or_else(|| anyhow::anyhow!("ChunkDecodedSizes count overflow"))?;
    let expected = count
        .checked_mul(8)
        .and_then(|x| x.checked_add(5))
        .ok_or_else(|| anyhow::anyhow!("ChunkDecodedSizes size overflow"))?;
    if payload.len() != expected {
        bail!("ChunkDecodedSizes len {} != {}", payload.len(), expected);
    }
    let mut sizes = Vec::with_capacity(count);
    let mut o = 5;
    for _ in 0..count {
        sizes.push(u64::from_le_bytes(payload[o..o + 8].try_into().unwrap()));
        o += 8;
    }
    Ok((n_mates, num_chunks, sizes))
}
/// Generous cap; real footers are KB–MB. Bounds allocation before parse.
pub(crate) const MAX_FOOTER_BYTES: u64 = 256 * 1024 * 1024;

impl ChunkDirectory {
    /// Serialize the directory body (no locator). Returns the bytes.
    pub fn write(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(FOOTER_PREFIX + self.entries.len() * ENTRY_BYTES);
        out.extend_from_slice(&self.num_reads.to_le_bytes());
        out.extend_from_slice(&self.num_chunks.to_le_bytes());
        out.extend_from_slice(&(self.entries.len() as u32).to_le_bytes());
        for e in &self.entries {
            out.extend_from_slice(&e.chunk_index.to_le_bytes());
            out.push(e.role as u8);
            out.push(e.mate);
            out.push(e.codec);
            out.extend_from_slice(&e.offset.to_le_bytes());
            out.extend_from_slice(&e.length.to_le_bytes());
            out.extend_from_slice(&e.record_count.to_le_bytes());
        }
        out
    }

    /// Parse the directory body. Must consume EXACTLY `data` (no trailing
    /// bytes). Bounds the entry count by the bytes available before allocating.
    pub fn parse(data: &[u8]) -> Result<ChunkDirectory> {
        if data.len() < FOOTER_PREFIX {
            bail!("v5 footer truncated (< {FOOTER_PREFIX} bytes)");
        }
        let num_reads = u64::from_le_bytes(data[0..8].try_into().unwrap());
        let num_chunks = u32::from_le_bytes(data[8..12].try_into().unwrap());
        let n = u32::from_le_bytes(data[12..16].try_into().unwrap()) as usize;

        let remaining = data.len() - FOOTER_PREFIX;
        if (n as u64) > (remaining as u64) / (ENTRY_BYTES as u64) {
            bail!("v5 footer entry_count {n} exceeds remaining {remaining} bytes");
        }
        if remaining != n * ENTRY_BYTES {
            bail!(
                "v5 footer has {} trailing bytes",
                remaining as i64 - (n * ENTRY_BYTES) as i64
            );
        }

        let mut entries = Vec::with_capacity(n);
        let mut o = FOOTER_PREFIX;
        for _ in 0..n {
            let chunk_index = u32::from_le_bytes(data[o..o + 4].try_into().unwrap());
            let role = StreamRole::from_u8(data[o + 4])
                .ok_or_else(|| anyhow::anyhow!("v5 footer: unknown role byte {}", data[o + 4]))?;
            let mate = data[o + 5];
            let codec = data[o + 6];
            let offset = u64::from_le_bytes(data[o + 7..o + 15].try_into().unwrap());
            let length = u64::from_le_bytes(data[o + 15..o + 23].try_into().unwrap());
            let record_count = u64::from_le_bytes(data[o + 23..o + 31].try_into().unwrap());
            entries.push(ChunkDirEntry {
                chunk_index,
                role,
                mate,
                codec,
                offset,
                length,
                record_count,
            });
            o += ENTRY_BYTES;
        }
        Ok(ChunkDirectory {
            num_reads,
            num_chunks,
            entries,
        })
    }
}

/// Read + CRC-validate a v5 archive's trailing footer (locator at EOF) and return
/// the parsed `ChunkDirectory`. This is the single source of truth for v5 footer
/// structural validation — single-end (`build_indices_from_footer`), paired, and
/// reference decode all call it, then layer their own per-mode entry validation on
/// top (it returns only the validated directory). `header_end` is the front-header
/// byte length (footer must start at/after it).
pub(crate) fn read_v5_footer(
    archive_path: &std::path::Path,
    header_end: u64,
) -> Result<ChunkDirectory> {
    use std::io::{Read, Seek, SeekFrom};
    let file_size = std::fs::metadata(archive_path)?.len();
    if file_size < header_end + LOCATOR_LEN as u64 {
        bail!("v5 archive too small for header + locator");
    }
    let mut file = std::fs::File::open(archive_path)?;
    file.seek(SeekFrom::Start(file_size - LOCATOR_LEN as u64))?;
    let mut loc = [0u8; LOCATOR_LEN];
    file.read_exact(&mut loc)?;
    if loc[12..20] != FOOTER_MAGIC {
        bail!("v5 archive truncated/corrupt: bad footer magic");
    }
    let footer_len = u64::from_le_bytes(loc[0..8].try_into().unwrap());
    let stored_crc = u32::from_le_bytes(loc[8..12].try_into().unwrap());
    if footer_len > MAX_FOOTER_BYTES {
        bail!("v5 footer_len {footer_len} exceeds cap {MAX_FOOTER_BYTES}");
    }
    let footer_start = (file_size - LOCATOR_LEN as u64)
        .checked_sub(footer_len)
        .ok_or_else(|| anyhow::anyhow!("v5 footer_len underflows file"))?;
    if footer_start < header_end {
        bail!("v5 footer overlaps header region");
    }
    file.seek(SeekFrom::Start(footer_start))?;
    let mut fbuf = vec![0u8; footer_len as usize];
    file.read_exact(&mut fbuf)?;
    validate_and_parse_footer(&fbuf, stored_crc, footer_start, header_end)
}

/// Per-entry byte-span validation + global non-overlap: every frame must lie
/// wholly within the payload region `[header_end, footer_start)`, and no two
/// frames may overlap. Shared by `validate_and_parse_footer` (decode) and the
/// single-end archive merge (which validates its rebuilt directory before writing
/// the footer).
pub(crate) fn validate_directory_spans(
    entries: &[ChunkDirEntry],
    header_end: u64,
    footer_start: u64,
) -> Result<()> {
    let mut spans: Vec<(u64, u64)> = Vec::with_capacity(entries.len());
    for (i, e) in entries.iter().enumerate() {
        let end = e
            .offset
            .checked_add(e.length)
            .ok_or_else(|| anyhow::anyhow!("v5 entry {i}: offset+length overflow"))?;
        if e.offset < header_end || end > footer_start {
            bail!(
                "v5 entry {i}: frame [{}, {end}) out of payload region [{header_end}, {footer_start})",
                e.offset
            );
        }
        spans.push((e.offset, end));
    }
    spans.sort_by_key(|s| s.0);
    for w in spans.windows(2) {
        if w[0].1 > w[1].0 {
            bail!("v5 directory: overlapping frame spans");
        }
    }
    Ok(())
}

/// CRC-check `footer`, parse the directory, then run the shared span/overlap/bounds
/// validation (`validate_directory_spans`) against `[header_end, footer_start)`.
/// Mode-agnostic; per-mode role/codec legality is the caller's job.
fn validate_and_parse_footer(
    footer: &[u8],
    stored_crc: u32,
    footer_start: u64,
    header_end: u64,
) -> Result<ChunkDirectory> {
    if footer_crc32(footer) != stored_crc {
        bail!("v5 footer CRC mismatch — archive corrupt");
    }
    let dir = ChunkDirectory::parse(footer)?;
    validate_directory_spans(&dir.entries, header_end, footer_start)?;
    Ok(dir)
}

/// Validate + parse a v5 footer from an in-memory archive (whole bytes). Mirrors
/// `read_v5_footer` (locator magic, footer_len cap, underflow, header overlap),
/// then funnels into the shared `validate_and_parse_footer` core. Reference decode
/// used to call this when it held the whole archive in memory; it now decodes
/// seekably (file-based footer via `read_v5_footer`), so this byte-slice core is
/// retained only as a tested utility (`#[cfg(test)]`) for the slice path.
#[cfg(test)]
pub(crate) fn parse_v5_footer_from_slice(
    archive: &[u8],
    header_end: u64,
) -> Result<ChunkDirectory> {
    let file_size = archive.len() as u64;
    if file_size < header_end + LOCATOR_LEN as u64 {
        bail!("v5 archive too small for header + locator");
    }
    let loc = &archive[archive.len() - LOCATOR_LEN..];
    if loc[12..20] != FOOTER_MAGIC {
        bail!("v5 archive truncated/corrupt: bad footer magic");
    }
    let footer_len = u64::from_le_bytes(loc[0..8].try_into().unwrap());
    let stored_crc = u32::from_le_bytes(loc[8..12].try_into().unwrap());
    if footer_len > MAX_FOOTER_BYTES {
        bail!("v5 footer_len {footer_len} exceeds cap {MAX_FOOTER_BYTES}");
    }
    let footer_start = (file_size - LOCATOR_LEN as u64)
        .checked_sub(footer_len)
        .ok_or_else(|| anyhow::anyhow!("v5 footer_len underflows file"))?;
    if footer_start < header_end {
        bail!("v5 footer overlaps header region");
    }
    let footer = &archive[footer_start as usize..(file_size - LOCATOR_LEN as u64) as usize];
    validate_and_parse_footer(footer, stored_crc, footer_start, header_end)
}

/// CRC32 of the footer body, via `flate2::Crc` (same primitive as v4 blocks).
pub(crate) fn footer_crc32(footer: &[u8]) -> u32 {
    use flate2::Crc;
    let mut c = Crc::new();
    c.update(footer);
    c.sum()
}

/// Serialize `dir`, append the footer body + 20-byte trailing locator
/// (`footer_len u64 | footer_crc32 u32 | magic`) to `out`, flush, and return the
/// footer body length in bytes (callers size archive metadata from it). This is the
/// write-side sibling of [`read_v5_footer`]: the single place the v5 footer tail is
/// emitted, shared by the single-end, paired, and reference encoders. Seek-free —
/// the footer's offset lives in the trailing locator, so the encoder never seeks
/// backward and the output works to a pipe/stdout. Mode-agnostic: any per-mode
/// directory self-check (`validate_paired_directory` / `validate_reference_directory`)
/// runs at the call site *before* this.
pub(crate) fn write_footer_and_locator(
    out: &mut dyn std::io::Write,
    dir: &ChunkDirectory,
) -> Result<usize> {
    let footer = dir.write();
    let crc = footer_crc32(&footer);
    out.write_all(&footer)?;
    out.write_all(&(footer.len() as u64).to_le_bytes())?;
    out.write_all(&crc.to_le_bytes())?;
    out.write_all(&FOOTER_MAGIC)?;
    out.flush()?;
    Ok(footer.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> ChunkDirectory {
        ChunkDirectory {
            num_reads: 7,
            num_chunks: 2,
            entries: vec![
                ChunkDirEntry {
                    chunk_index: 0,
                    role: StreamRole::Headers,
                    mate: 0,
                    codec: 1,
                    offset: 60,
                    length: 112,
                    record_count: 4,
                },
                ChunkDirEntry {
                    chunk_index: 0,
                    role: StreamRole::ReferenceMeta,
                    mate: 1,
                    codec: 1,
                    offset: 172,
                    length: 212,
                    record_count: 4,
                },
                ChunkDirEntry {
                    chunk_index: 1,
                    role: StreamRole::HeaderDelta,
                    mate: 2,
                    codec: 1,
                    offset: 384,
                    length: 80,
                    record_count: 3,
                },
            ],
        }
    }

    #[test]
    fn entry_is_31_bytes() {
        assert_eq!(ENTRY_BYTES, 31);
    }

    #[test]
    fn reference_roles_and_global_sentinel_roundtrip() {
        let d = ChunkDirectory {
            num_reads: 9,
            num_chunks: 1,
            entries: vec![
                ChunkDirEntry {
                    chunk_index: GLOBAL_SENTINEL,
                    role: StreamRole::PackedBacking,
                    mate: 0,
                    codec: 5,
                    offset: 60,
                    length: 100,
                    record_count: 5_000_000_000,
                },
                ChunkDirEntry {
                    chunk_index: 0,
                    role: StreamRole::Positions,
                    mate: 1,
                    codec: 1,
                    offset: 160,
                    length: 40,
                    record_count: 9,
                },
            ],
        };
        let bytes = d.write();
        assert_eq!(ChunkDirectory::parse(&bytes).unwrap(), d);
    }

    #[test]
    fn directory_roundtrips() {
        let d = sample();
        let bytes = d.write();
        let parsed = ChunkDirectory::parse(&bytes).unwrap();
        assert_eq!(parsed, d);
    }

    #[test]
    fn rejects_trailing_bytes() {
        let mut bytes = sample().write();
        bytes.push(0xAB);
        assert!(ChunkDirectory::parse(&bytes).is_err());
    }

    #[test]
    fn footer_slice_core_roundtrips() {
        // Build a fake archive: [header_end bytes][...payload...][footer][locator].
        let header_end = 8u64;
        let dir = ChunkDirectory {
            num_reads: 0,
            num_chunks: 0,
            entries: vec![],
        };
        let footer = dir.write();
        let mut archive = vec![0u8; header_end as usize]; // stand-in header
        let footer_start = archive.len() as u64;
        archive.extend_from_slice(&footer);
        archive.extend_from_slice(&(footer.len() as u64).to_le_bytes());
        archive.extend_from_slice(&footer_crc32(&footer).to_le_bytes());
        archive.extend_from_slice(&FOOTER_MAGIC);
        let parsed = parse_v5_footer_from_slice(&archive, header_end).unwrap();
        assert_eq!(parsed, dir);
        let _ = footer_start;
    }

    #[test]
    fn decoded_sizes_payload_roundtrips() {
        let sizes = [1234u64, 0, 999_999_999];
        let p = encode_decoded_sizes(1, 3, &sizes);
        assert_eq!(
            parse_decoded_sizes(&p).unwrap(),
            (1u8, 3u32, sizes.to_vec())
        );
    }

    #[test]
    fn decoded_sizes_payload_rejects_corruption() {
        let p = encode_decoded_sizes(1, 2, &[10, 20]);
        assert!(parse_decoded_sizes(&p[..p.len() - 1]).is_err()); // truncated
        let mut b = p.clone();
        b[0] = 0;
        assert!(parse_decoded_sizes(&b).is_err()); // n_mates 0
        let mut b = p.clone();
        b[0] = 3;
        assert!(parse_decoded_sizes(&b).is_err()); // n_mates 3
        let mut b = p.clone();
        b[1] = 5;
        assert!(parse_decoded_sizes(&b).is_err()); // count mismatch
    }

    #[test]
    fn chunk_decoded_sizes_role_roundtrips() {
        assert_eq!(StreamRole::from_u8(18), Some(StreamRole::ChunkDecodedSizes));
        assert_eq!(StreamRole::ChunkDecodedSizes as u8, 18);
        // 19/20/21 are the cluster-mode roles (see cluster_role_tests); 22 is the
        // first unassigned value.
        assert_eq!(StreamRole::from_u8(22), None);
    }

    #[test]
    fn rejects_overlarge_entry_count() {
        // Hand-craft a header claiming a huge entry count with no entry bytes.
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&0u64.to_le_bytes()); // num_reads
        bytes.extend_from_slice(&1u32.to_le_bytes()); // num_chunks
        bytes.extend_from_slice(&u32::MAX.to_le_bytes()); // n_entries
        assert!(ChunkDirectory::parse(&bytes).is_err());
    }

    #[test]
    fn validate_directory_spans_accepts_and_rejects() {
        // Two adjacent in-region frames: OK.
        let ok = vec![
            ChunkDirEntry { chunk_index: 0, role: StreamRole::Headers, mate: 0, codec: 1, offset: 10, length: 20, record_count: 1 },
            ChunkDirEntry { chunk_index: 0, role: StreamRole::Sequence, mate: 0, codec: 1, offset: 30, length: 20, record_count: 1 },
        ];
        assert!(validate_directory_spans(&ok, 10, 50).is_ok());
        // Out of region (before header_end).
        assert!(validate_directory_spans(&ok, 11, 50).is_err());
        // Out of region (past footer_start).
        assert!(validate_directory_spans(&ok, 10, 49).is_err());
        // Overlapping frames.
        let overlap = vec![
            ChunkDirEntry { chunk_index: 0, role: StreamRole::Headers, mate: 0, codec: 1, offset: 10, length: 25, record_count: 1 },
            ChunkDirEntry { chunk_index: 0, role: StreamRole::Sequence, mate: 0, codec: 1, offset: 30, length: 20, record_count: 1 },
        ];
        assert!(validate_directory_spans(&overlap, 10, 50).is_err());
        // offset+length overflow.
        let overflow = vec![
            ChunkDirEntry { chunk_index: 0, role: StreamRole::Headers, mate: 0, codec: 1, offset: u64::MAX, length: 1, record_count: 1 },
        ];
        assert!(validate_directory_spans(&overflow, 0, u64::MAX).is_err());
    }
}

#[cfg(test)]
mod cluster_role_tests {
    use super::*;
    #[test]
    fn cluster_roles_roundtrip() {
        for (v, r) in [
            (19u8, StreamRole::ClusteredSeq),
            (20, StreamRole::StrandBits),
            (21, StreamRole::ClusterMeta),
        ] {
            assert_eq!(StreamRole::from_u8(v), Some(r));
            assert_eq!(r as u8, v);
        }
    }
}
