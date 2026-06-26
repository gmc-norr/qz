//! Per-(chunk, mate) position-stream encode/decode for genome coordinates (spec §4.3.2).
//!
//! ## Wire format
//!
//! **Mate-1 stream:** for each mapped R1 in order: `[ref_id: varint][ref_pos: varint]`.
//!
//! **Mate-2 stream:** `[form_bits: ceil(M/8) bytes]` (bit j, LSB-first, for the
//! j-th mapped R2), then per mapped R2 in order:
//! - `form=0` (same-contig delta): `[delta: svarint]` = `ref_pos_R2 − ref_pos_R1`.
//!   Used iff the pair's R1 is mapped AND `ref_id_R2 == ref_id_R1`.
//! - `form=1` (absolute escape): `[ref_id: varint][ref_pos: varint]`. Used when
//!   R1 is unmapped or the contigs differ.
//!
//! `M` = `r1_anchor.len()` = number of mapped R2 reads in this group.
//! On decode, unused high bits in the last form byte must be zero.

use crate::compression::dna_utils::{read_svarint, read_varint, write_svarint, write_varint};
use anyhow::{Result, bail};

// ---------------------------------------------------------------------------
// Mate 1: absolute (ref_id, ref_pos) per mapped R1
// ---------------------------------------------------------------------------

/// Encode mapped R1 positions into a byte stream.
///
/// For each `(ref_id, ref_pos)` pair: writes `ref_id` as a varint followed by
/// `ref_pos` as a varint.
pub fn encode_positions_mate1<I: Iterator<Item = (u32, u64)>>(mapped_r1: I) -> Vec<u8> {
    let mut out = Vec::new();
    for (ref_id, ref_pos) in mapped_r1 {
        write_varint(&mut out, ref_id as u64);
        write_varint(&mut out, ref_pos);
    }
    out
}

/// Decode `m` mapped R1 positions from a byte slice.
///
/// Returns `Err` if the slice is too short or contains trailing bytes.
pub fn decode_positions_mate1(bytes: &[u8], m: usize) -> Result<Vec<(u32, u64)>> {
    // Clamp the speculative capacity to the payload size: every record consumes at
    // least one payload byte, so a valid `m` is <= bytes.len(). This bounds the
    // allocation if a caller ever passes an unvalidated record count (defense in
    // depth — the loop below still validates each varint and errors on truncation).
    let mut out = Vec::with_capacity(m.min(bytes.len() + 1));
    let mut o = 0usize;
    for i in 0..m {
        let ref_id_u64 = read_varint(bytes, &mut o)
            .ok_or_else(|| anyhow::anyhow!("mate1 positions: truncated ref_id at record {i}"))?;
        let ref_id = u32::try_from(ref_id_u64).map_err(|_| {
            anyhow::anyhow!("mate1 positions: ref_id {ref_id_u64} exceeds u32 at record {i}")
        })?;
        let ref_pos = read_varint(bytes, &mut o)
            .ok_or_else(|| anyhow::anyhow!("mate1 positions: truncated ref_pos at record {i}"))?;
        out.push((ref_id, ref_pos));
    }
    if o != bytes.len() {
        bail!(
            "mate1 positions: {} trailing bytes after {m} records",
            bytes.len() - o
        );
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Mate 2: leading form-bits bitplane + delta-or-absolute per mapped R2
// ---------------------------------------------------------------------------

/// Encode mapped R2 positions into a byte stream (spec §4.3.2 mate-2 wire).
///
/// `r1_anchor[i]` is the `(ref_id, ref_pos)` of the paired R1 iff that R1 is
/// mapped; `None` if R1 is unmapped. `mapped_r2[i]` is the `(ref_id, ref_pos)`
/// of the i-th mapped R2. `r1_anchor.len()` must equal `mapped_r2.len()` (= M).
///
/// Form selection:
/// - `form=0` (delta): R1 is mapped AND `ref_id_R2 == ref_id_R1`.
/// - `form=1` (absolute): R1 is unmapped OR contigs differ.
pub fn encode_positions_mate2(
    r1_anchor: &[Option<(u32, u64)>],
    mapped_r2: &[(u32, u64)],
) -> Vec<u8> {
    let m = mapped_r2.len();
    debug_assert_eq!(r1_anchor.len(), m);

    // Build form bits: bit j = 0 → delta, bit j = 1 → absolute.
    let form_bytes_len = m.div_ceil(8);
    let mut form_bits = vec![0u8; form_bytes_len];
    for (i, (&(r2_ref_id, _), anchor)) in mapped_r2.iter().zip(r1_anchor.iter()).enumerate() {
        let is_absolute = match anchor {
            Some((r1_ref_id, _)) => r2_ref_id != *r1_ref_id,
            None => true,
        };
        if is_absolute {
            form_bits[i / 8] |= 1u8 << (i % 8);
        }
    }

    let mut out = Vec::new();
    out.extend_from_slice(&form_bits);

    // Per-record payload.
    for (i, &(r2_ref_id, r2_ref_pos)) in mapped_r2.iter().enumerate() {
        let form = (form_bits[i / 8] >> (i % 8)) & 1;
        if form == 0 {
            // Delta: r1_anchor is Some and same contig.
            let (_, r1_ref_pos) = r1_anchor[i].expect("form=0 requires mapped R1");
            let delta = r2_ref_pos as i64 - r1_ref_pos as i64;
            write_svarint(&mut out, delta);
        } else {
            // Absolute escape.
            write_varint(&mut out, r2_ref_id as u64);
            write_varint(&mut out, r2_ref_pos);
        }
    }

    out
}

/// Decode `m` mapped R2 positions from a byte stream (spec §4.3.2 mate-2 wire).
///
/// `r1_anchor[i]` must be the same slice that was passed to `encode_positions_mate2`
/// (one entry per mapped R2). Returns `Err` on:
/// - `r1_anchor.len() != m`
/// - form-bits region shorter than `ceil(m/8)` bytes
/// - non-zero unused high bits in the last form byte
/// - `form=0` record whose `r1_anchor[i]` is `None`
/// - truncated payload
/// - trailing bytes after all records
pub fn decode_positions_mate2(
    bytes: &[u8],
    r1_anchor: &[Option<(u32, u64)>],
    m: usize,
) -> Result<Vec<(u32, u64)>> {
    if r1_anchor.len() != m {
        bail!(
            "mate2 positions: r1_anchor.len()={} != m={m}",
            r1_anchor.len()
        );
    }

    let form_bytes_len = m.div_ceil(8);
    if bytes.len() < form_bytes_len {
        bail!(
            "mate2 positions: need {form_bytes_len} form bytes for m={m}, got {}",
            bytes.len()
        );
    }
    let form_bits = &bytes[..form_bytes_len];

    // Validate unused high bits in the last byte are zero.
    if m > 0 {
        let used_bits_in_last = m % 8;
        if used_bits_in_last != 0 {
            let last = form_bits[form_bytes_len - 1];
            let mask = !((1u8 << used_bits_in_last) - 1);
            if last & mask != 0 {
                bail!("mate2 positions: non-zero unused high bits in last form byte");
            }
        }
    }

    let mut o = form_bytes_len;
    // Clamp the speculative capacity to the payload size: every record consumes at
    // least one payload byte, so a valid `m` is <= bytes.len(). This bounds the
    // allocation if a caller ever passes an unvalidated record count (defense in
    // depth — the loop below still validates each varint and errors on truncation).
    let mut out = Vec::with_capacity(m.min(bytes.len() + 1));

    for i in 0..m {
        let form = (form_bits[i / 8] >> (i % 8)) & 1;
        if form == 0 {
            // Delta: requires a mapped R1 on the same contig.
            let (r1_ref_id, r1_ref_pos) = r1_anchor[i].ok_or_else(|| {
                anyhow::anyhow!(
                    "mate2 positions: form=0 (delta) for record {i} but r1_anchor[{i}] is None"
                )
            })?;
            let delta = read_svarint(bytes, &mut o).ok_or_else(|| {
                anyhow::anyhow!("mate2 positions: truncated delta svarint at record {i}")
            })?;
            let ref_pos = (r1_ref_pos as i64 + delta) as u64;
            out.push((r1_ref_id, ref_pos));
        } else {
            // Absolute escape.
            let ref_id_u64 = read_varint(bytes, &mut o).ok_or_else(|| {
                anyhow::anyhow!("mate2 positions: truncated ref_id varint at record {i}")
            })?;
            let ref_id = u32::try_from(ref_id_u64).map_err(|_| {
                anyhow::anyhow!("mate2 positions: ref_id {ref_id_u64} exceeds u32 at record {i}")
            })?;
            let ref_pos = read_varint(bytes, &mut o).ok_or_else(|| {
                anyhow::anyhow!("mate2 positions: truncated ref_pos varint at record {i}")
            })?;
            out.push((ref_id, ref_pos));
        }
    }

    if o != bytes.len() {
        bail!(
            "mate2 positions: {} trailing bytes after {m} records",
            bytes.len() - o
        );
    }

    Ok(out)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r1_absolute_roundtrip() {
        let mapped = vec![(0u32, 1000u64), (3u32, 42u64)];
        let bytes = encode_positions_mate1(mapped.iter().copied());
        let got = decode_positions_mate1(&bytes, 2).unwrap();
        assert_eq!(got, mapped);
    }

    #[test]
    fn r2_delta_escape_and_unmapped_r1_roundtrip() {
        // 3 mapped R2; their pairs' R1 anchors:
        let r1_anchor = vec![Some((0u32, 1000u64)), None, Some((1u32, 500u64))];
        let mapped_r2 = vec![(0u32, 1100u64), (2u32, 7u64), (4u32, 9u64)];
        // R2[0]: same contig as R1 (0==0) -> delta 100; R2[1]: R1 unmapped -> absolute; R2[2]: diff contig (4!=1) -> absolute
        let bytes = encode_positions_mate2(&r1_anchor, &mapped_r2);
        let got = decode_positions_mate2(&bytes, &r1_anchor, 3).unwrap();
        assert_eq!(got, mapped_r2);
    }

    #[test]
    fn r2_form0_with_unmapped_r1_is_rejected_on_decode() {
        // craft a stream that claims form=0 for a read whose r1_anchor is None
        let r1_anchor = vec![None];
        // form_bits = 1 byte, bit0=0 (form0); then a delta svarint
        let mut bytes = vec![0u8]; // form bit 0 -> form0 (delta)
        crate::compression::dna_utils::write_svarint(&mut bytes, 5);
        assert!(decode_positions_mate2(&bytes, &r1_anchor, 1).is_err());
    }

    #[test]
    fn form_bits_wrong_length_rejected() {
        // m=9 needs ceil(9/8)=2 form bytes; give a stream too short
        let r1_anchor = vec![Some((0u32, 0u64)); 9];
        assert!(decode_positions_mate2(&[0u8], &r1_anchor, 9).is_err());
    }
}
