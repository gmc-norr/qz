//! Full-contig ReferenceMeta (spec §4.3.5). Stores the complete contig list in
//! ref_id order (incl. uncovered) so `ref_id < num_refs` and `ref_pos+read_len`
//! bounds are exact. Digest is informational (reference-free decode).

use crate::compression::dna_utils::{read_varint, write_varint};
use anyhow::Result;
use strobealign::io::fasta::RefSequence;

const META_VERSION: u8 = 1;
/// ReferenceMeta block-stream version that marks the per-(chunk,mate) FallbackPool
/// archive layout. Bumped from the legacy global-pool layout (version 1) so a
/// decoder can reject old archives that would silently mis-decode their pool.
pub const REFMETA_VERSION_PERCHUNK_POOL: u8 = 2;
const MAX_NUM_REFS: u64 = 1 << 24;
const MAX_NAME_LEN: u64 = 1024;

/// Reject a refmeta record count above the structural cap (`MAX_NUM_REFS`).
/// Split out of `read_refmeta_blockstream` so the boundary stays unit-testable
/// without materialising 2^24 records.
fn ensure_num_refs_within_cap(n: usize) -> Result<()> {
    if n as u64 > MAX_NUM_REFS {
        anyhow::bail!("refmeta blockstream: num_refs {n} > cap {MAX_NUM_REFS}");
    }
    Ok(())
}

fn fnv1a64(seed: u64, bytes: &[u8]) -> u64 {
    let mut h = seed;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(1099511628211);
    }
    h
}

pub fn reference_digest(refs: &[RefSequence]) -> u64 {
    let mut h = 1469598103934665603u64; // FNV offset basis
    for r in refs {
        h = fnv1a64(h, &(r.sequence.len() as u64).to_le_bytes());
        h = fnv1a64(h, &r.sequence);
    }
    h
}

pub struct ReferenceMeta {
    pub num_refs: u64,
    lens: Vec<u64>,
}

impl ReferenceMeta {
    pub fn contig_len(&self, ref_id: u32) -> u64 {
        self.lens[ref_id as usize]
    }
}

/// Write a `ReferenceMeta` as a record-boundary block-stream.
///
/// Wire layout: `[meta_version: 1B = 2][digest: 8B LE]` prefix followed by a
/// record-block-stream whose records are per-contig `[name_len:varint][name][contig_len:varint]`.
/// `num_refs` is the total record count (not stored separately).
/// Bounds: total records ≤ `2^24`.  Returns num_refs (= refs.len()).
///
/// The version byte is [`REFMETA_VERSION_PERCHUNK_POOL`] (2): it doubles as the
/// archive format-revision discriminator for the per-(chunk,mate) FallbackPool
/// layout, so `read_refmeta_blockstream` can reject the legacy global-pool
/// version (1) rather than mis-decode it.
pub fn write_refmeta_blockstream(
    out: &mut Vec<u8>,
    refs: &[RefSequence],
    target_bytes: usize,
) -> Result<u64> {
    // Prefix: meta_version (1 byte) + digest (8 bytes LE).
    out.push(REFMETA_VERSION_PERCHUNK_POOL);
    out.extend_from_slice(&reference_digest(refs).to_le_bytes());

    // Build per-contig records.
    let records: Vec<Vec<u8>> = refs
        .iter()
        .map(|r| {
            let mut rec = Vec::new();
            write_varint(&mut rec, r.name.len() as u64);
            rec.extend_from_slice(r.name.as_bytes());
            write_varint(&mut rec, r.sequence.len() as u64);
            rec
        })
        .collect();

    let n = super::backing::write_record_block_stream(
        out,
        records.iter().map(|v| v.as_slice()),
        target_bytes,
    )?;
    Ok(n)
}

/// Read a `ReferenceMeta` block-stream written by `write_refmeta_blockstream`.
///
/// Reads `meta_version` + 8-byte digest prefix, then `read_record_block_stream`
/// on the remainder. The version byte doubles as the archive format-revision
/// discriminator: it MUST be [`REFMETA_VERSION_PERCHUNK_POOL`] (2). The legacy
/// global-fallback-pool layout (version 1) is rejected with a recompress hint
/// (decoding it under the per-(chunk,mate) reader would mis-attribute literals);
/// any other value is corrupt. Enforces:
/// - total records ≤ `2^24`
/// - `name_len ≤ 1024` per record
/// - no trailing bytes in each record
pub fn read_refmeta_blockstream(data: &[u8]) -> Result<ReferenceMeta> {
    // Parse the fixed prefix: meta_version (1 byte) + digest (8 bytes LE).
    if data.is_empty() {
        anyhow::bail!("refmeta blockstream: truncated version byte");
    }
    let ver = data[0];
    if ver == META_VERSION {
        anyhow::bail!(
            "reference archive uses the legacy global fallback-pool format (meta_version 1); recompress with the current qz"
        );
    }
    if ver != REFMETA_VERSION_PERCHUNK_POOL {
        anyhow::bail!("refmeta blockstream: unsupported version {ver}");
    }
    if data.len() < 9 {
        anyhow::bail!("refmeta blockstream: truncated digest");
    }
    // digest is informational — we skip bytes 1..9.
    let rest = &data[9..];

    let raw_records = super::backing::read_record_block_stream(rest)?;
    ensure_num_refs_within_cap(raw_records.len())?;

    let mut lens = Vec::with_capacity(raw_records.len());

    for (i, rec) in raw_records.iter().enumerate() {
        let mut off = 0usize;
        let nl = read_varint(rec, &mut off)
            .ok_or_else(|| anyhow::anyhow!("refmeta blockstream record {i}: truncated name_len"))?;
        if nl > MAX_NAME_LEN {
            anyhow::bail!("refmeta blockstream record {i}: name_len {nl} > cap {MAX_NAME_LEN}");
        }
        let nl = nl as usize;
        if off + nl > rec.len() {
            anyhow::bail!(
                "refmeta blockstream record {i}: name OOB (need {nl} bytes, have {})",
                rec.len() - off
            );
        }
        // Validate the name bytes are present, but the reference-free FASTQ decode never
        // needs the contig names, so skip over them rather than allocate a String.
        off += nl;
        let cl = read_varint(rec, &mut off).ok_or_else(|| {
            anyhow::anyhow!("refmeta blockstream record {i}: truncated contig_len")
        })?;
        if off != rec.len() {
            anyhow::bail!(
                "refmeta blockstream record {i}: {} trailing bytes",
                rec.len() - off
            );
        }
        lens.push(cl);
    }

    Ok(ReferenceMeta {
        num_refs: raw_records.len() as u64,
        lens,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn refs(v: &[(&str, &[u8])]) -> Vec<strobealign::io::fasta::RefSequence> {
        v.iter()
            .map(|(n, s)| strobealign::io::fasta::RefSequence {
                name: n.to_string(),
                sequence: s.to_vec(),
            })
            .collect()
    }

    #[test]
    fn reference_digest_is_deterministic_and_boundary_safe() {
        // `reference_digest` length-prefixes each contig's sequence, so a split
        // can't collide with its concatenation: [("a","AC"),("b","GT")] must
        // differ from [("a","ACGT")].
        let a = reference_digest(&refs(&[("a", b"AC"), ("b", b"GT")]));
        let c = reference_digest(&refs(&[("a", b"ACGT")]));
        assert_ne!(a, c);
        // Deterministic: same input → same digest.
        let a2 = reference_digest(&refs(&[("a", b"AC"), ("b", b"GT")]));
        assert_eq!(a, a2);
    }

    #[test]
    fn ensure_num_refs_within_cap_enforces_boundary() {
        // The structural cap (2^24) is too large to materialise, so cover the
        // boundary via the extracted helper the live reader calls.
        assert!(ensure_num_refs_within_cap(0).is_ok());
        assert!(ensure_num_refs_within_cap(MAX_NUM_REFS as usize).is_ok());
        assert!(ensure_num_refs_within_cap(MAX_NUM_REFS as usize + 1).is_err());
    }

    #[test]
    fn refmeta_blockstream_rejects_oversized_name() {
        // The live reader caps per-contig name_len at MAX_NAME_LEN (1024). Write a
        // valid block-stream whose single contig has a 2000-char name (the writer
        // imposes no cap), then confirm the reader rejects it rather than reading it.
        let long_name = "x".repeat(2000);
        let r = refs(&[(long_name.as_str(), b"ACGT")]);
        let mut out = Vec::new();
        write_refmeta_blockstream(&mut out, &r, 1 << 20).unwrap();
        assert!(read_refmeta_blockstream(&out).is_err());
    }

    #[test]
    fn refmeta_blockstream_roundtrip() {
        let r = refs(&[("chr1", b"ACGT"), ("chr2", b"NNNN"), ("chr3", b"ACGTACGT")]);
        let mut out = Vec::new();
        let n = write_refmeta_blockstream(&mut out, &r, 8).unwrap();
        assert_eq!(n, 3);
        let meta = read_refmeta_blockstream(&out).unwrap();
        assert_eq!(meta.num_refs, 3);
        assert_eq!(meta.contig_len(2), 8);
    }
}
