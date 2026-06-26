//! Columnar substitution-edit codec for reference-based read encoding.
//!
//! Mapped reads are stored as substitutions against a consensus slice: a
//! per-read list of `Edit { pos, base }`. The codec is COLUMNAR — three parallel
//! byte streams that compress independently downstream:
//!   - `edit_count`: one `write_varint(count)` per read (`num_reads` entries).
//!   - `edit_pos`:   concatenated within-read positions, each a `write_varint`.
//!     Positions are **ABSOLUTE within the read** (0-based offset into the read),
//!     written in the order the edits appear (which is ascending, since
//!     `extract_substitutions` scans left-to-right). They are **NOT**
//!     offset-from-previous / delta-encoded. This is a deliberate patent
//!     design-around vs ORA's 1-byte offset-from-previous-mismatch record — do
//!     NOT change positions to deltas.
//!   - `edit_base`:  concatenated 3-bit ACGTN codes, LSB-first bitpacked (see
//!     below).
//!
//! `reconstruct_read` overlays the edits onto the consensus slice, then
//! reverse-complements iff the read mapped to the reverse strand. The byte-exact
//! reconstruct is the losslessness backstop.
//!
//! ## 3-bit ACGTN packing (`edit_base`)
//! Codes: A=0, C=1, G=2, T=3, N=4. Codes 5..=7 are invalid and rejected on
//! decode. The 3 bits of edit `k` (across the whole concatenated `edit_base`
//! stream, k = 0..total_edits) occupy bit positions `3*k .. 3*k+3`, LSB-first:
//! bit `b` of the code lands in byte `(3*k + b) / 8` at intra-byte bit
//! `(3*k + b) % 8` (where intra-byte bit 0 is the least-significant bit). A code
//! therefore straddles a byte boundary whenever `3*k % 8 > 5`. The packed stream
//! is exactly `ceil(total_edits * 3 / 8)` bytes; any unused high bits of the last
//! byte are zero.
//!
//! ## OOB contract for `reconstruct_read`
//! `reconstruct_read` is the TRUSTED compress-side path: it `debug_assert!`s that
//! every `edit.pos < consensus_slice.len()` and will index-panic in release if
//! that is violated. `decode_edits` does NOT (and cannot) validate
//! `pos < read_len` — it has no read length. **The untrusted-decode caller
//! (Task 13) MUST validate `edit.pos < read_len` for every decoded edit BEFORE
//! calling `reconstruct_read`.**

use super::backing::PackedBacking;
use crate::compression::dna_utils::{
    read_varint, reverse_complement, reverse_complement_in_place, write_varint,
};
use anyhow::{Result, bail};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Edit {
    pub pos: u32,
    pub base: u8,
}

pub struct EncodedEdits {
    pub edit_count: Vec<u8>,
    pub edit_pos: Vec<u8>,
    pub edit_base: Vec<u8>,
}

/// 3-bit code for an ACGTN base, or `None` for anything else.
fn base_to_code(b: u8) -> Option<u8> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        b'N' => Some(4),
        _ => None,
    }
}

/// Inverse of `base_to_code`. Codes 0..=4 map to ACGTN; >4 is invalid.
fn code_to_base(c: u8) -> Option<u8> {
    match c {
        0 => Some(b'A'),
        1 => Some(b'C'),
        2 => Some(b'G'),
        3 => Some(b'T'),
        4 => Some(b'N'),
        _ => None,
    }
}

/// True iff `b` is an ACGTN base — the only bases the edit codec can represent.
/// `place_read` uses this to validate the (offset-independent) read once instead
/// of re-checking every base at every candidate offset.
pub fn base_is_acgtn(b: u8) -> bool {
    base_to_code(b).is_some()
}

/// Collect substitutions where `consensus_slice != r_fwd`, scanning left-to-right
/// (so positions are ascending). Equal length required.
///
/// Returns `None` if the lengths differ or any read base is not in ACGTN.
pub fn extract_substitutions(consensus_slice: &[u8], r_fwd: &[u8]) -> Option<Vec<Edit>> {
    if consensus_slice.len() != r_fwd.len() {
        return None;
    }
    let mut edits = Vec::new();
    for (i, (&c, &r)) in consensus_slice.iter().zip(r_fwd.iter()).enumerate() {
        // Reject any read base outside ACGTN up front — the codec can only
        // represent those, and a silent loss here would break losslessness.
        base_to_code(r)?;
        if c != r {
            edits.push(Edit {
                pos: i as u32,
                base: r,
            });
        }
    }
    Some(edits)
}

/// Overlay `edits` onto a copy of `consensus_slice`, then reverse-complement iff
/// `is_rc`.
///
/// TRUSTED path: `debug_assert!`s `pos < len`. Untrusted decoders MUST validate
/// `edit.pos < read_len` before calling (see module docs / `decode_edits`).
pub fn reconstruct_read(consensus_slice: &[u8], edits: &[Edit], is_rc: bool) -> Vec<u8> {
    let mut out = consensus_slice.to_vec();
    for e in edits {
        debug_assert!(
            (e.pos as usize) < out.len(),
            "reconstruct_read: edit pos {} out of bounds (len {}); caller must validate pos < read_len",
            e.pos,
            out.len()
        );
        out[e.pos as usize] = e.base;
    }
    if is_rc { reverse_complement(&out) } else { out }
}

/// Single-allocation equivalent of
/// `reconstruct_read(&backing.slice(base, len), edits, is_rc)`: builds the
/// consensus slice, overlays `edits`, and reverse-complements (iff `is_rc`) all in
/// ONE allocation — versus the two-step path's three (`slice()` Vec + `to_vec()`
/// copy + `reverse_complement()` Vec). Byte-for-byte identical to the two-step
/// path. Same TRUSTED OOB contract as [`reconstruct_read`] (`debug_assert!`s
/// `edit.pos < len`; untrusted callers must pre-validate).
pub fn reconstruct_read_from_backing(
    backing: &PackedBacking,
    base: usize,
    len: usize,
    edits: &[Edit],
    is_rc: bool,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    backing.fill_slice(base, len, &mut out);
    for e in edits {
        debug_assert!(
            (e.pos as usize) < out.len(),
            "reconstruct_read_from_backing: edit pos {} out of bounds (len {}); caller must validate pos < read_len",
            e.pos,
            out.len()
        );
        out[e.pos as usize] = e.base;
    }
    if is_rc {
        reverse_complement_in_place(&mut out);
    }
    out
}

/// Encode per-read edit lists into the three columnar streams.
pub fn encode_edits(recs: &[&[Edit]]) -> EncodedEdits {
    let mut edit_count = Vec::new();
    let mut edit_pos = Vec::new();
    // Bitpack 3-bit base codes LSB-first into edit_base.
    let mut edit_base = Vec::new();
    let mut bit_len: usize = 0; // number of bits written so far

    for rec in recs {
        write_varint(&mut edit_count, rec.len() as u64);
        for e in rec.iter() {
            // ABSOLUTE within-read position — NOT delta-encoded (patent design-around).
            write_varint(&mut edit_pos, e.pos as u64);
            let code = base_to_code(e.base)
                .expect("encode_edits: edit base must be ACGTN (validated at extract time)");
            // Write 3 bits LSB-first at bit offset `bit_len`.
            for b in 0..3 {
                let bit = (code >> b) & 1;
                let byte_idx = bit_len / 8;
                let bit_idx = bit_len % 8;
                if byte_idx == edit_base.len() {
                    edit_base.push(0u8);
                }
                edit_base[byte_idx] |= bit << bit_idx;
                bit_len += 1;
            }
        }
    }

    EncodedEdits {
        edit_count,
        edit_pos,
        edit_base,
    }
}

/// Invert `encode_edits`. Defensive against hostile input: validates that the
/// position stream is consumed exactly (no trailing bytes, no mid-varint
/// truncation), that the base bitstream is exactly `ceil(total*3/8)` bytes, and
/// that every 3-bit base code is a valid ACGTN code (0..=4).
///
/// Does NOT validate `edit.pos < read_len` — it has no read length. The untrusted
/// caller (Task 13) MUST check `edit.pos < read_len` for every decoded edit
/// before calling `reconstruct_read`.
pub fn decode_edits(
    edit_count: &[u8],
    edit_pos: &[u8],
    edit_base: &[u8],
    num_reads: usize,
) -> Result<Vec<Vec<Edit>>> {
    // 1. Read num_reads counts, exact consumption of edit_count.
    //    `num_reads` is an untrusted footer-derived count (the v5 footer CRC is
    //    recomputable, not a MAC). Each count is a varint of >= 1 byte, so the
    //    real count can't exceed `edit_count.len()`; cap the capacity hint so a
    //    hostile `num_reads` can't force a multi-GB over-allocation here.
    let mut co = 0usize;
    let mut counts = Vec::with_capacity(num_reads.min(edit_count.len() + 1));
    let mut total: usize = 0;
    for _ in 0..num_reads {
        let c = read_varint(edit_count, &mut co)
            .ok_or_else(|| anyhow::anyhow!("decode_edits: truncated edit_count stream"))?;
        let c = usize::try_from(c)
            .map_err(|_| anyhow::anyhow!("decode_edits: edit_count overflows usize"))?;
        total = total
            .checked_add(c)
            .ok_or_else(|| anyhow::anyhow!("decode_edits: total edit count overflows"))?;
        counts.push(c);
    }
    if co != edit_count.len() {
        bail!("decode_edits: trailing bytes in edit_count stream");
    }

    // 2. Read exactly `total` positions from edit_pos, exact consumption.
    //    `total` is the sum of the untrusted per-read edit counts above — a single
    //    huge varint in `edit_count` can make it enormous. Each position is a
    //    varint of >= 1 byte, so the achievable count can't exceed `edit_pos.len()`;
    //    cap the capacity hint (the loop below still proves the exact count by
    //    consuming `edit_pos` to its end) so the alloc can't overflow/OOM before
    //    that proof runs.
    let mut po = 0usize;
    let mut positions = Vec::with_capacity(total.min(edit_pos.len() + 1));
    for _ in 0..total {
        let p = read_varint(edit_pos, &mut po)
            .ok_or_else(|| anyhow::anyhow!("decode_edits: truncated edit_pos stream"))?;
        let p = u32::try_from(p)
            .map_err(|_| anyhow::anyhow!("decode_edits: edit position overflows u32"))?;
        positions.push(p);
    }
    if po != edit_pos.len() {
        bail!("decode_edits: trailing bytes in edit_pos stream");
    }

    // 3. Read exactly `total` 3-bit codes from edit_base.
    //    Require the packed stream length is exactly ceil(total*3/8); a trailing
    //    full byte beyond that is garbage and rejected. Pad bits in the final
    //    byte are ignored (they are zero on encode).
    // `total * 3` can overflow usize on a hostile count; use checked arithmetic
    // (a real archive's `total` is bounded by the consensus size).
    let expected_base_len = total
        .checked_mul(3)
        .ok_or_else(|| {
            anyhow::anyhow!("decode_edits: edit count {total} too large (bit overflow)")
        })?
        .div_ceil(8);
    if edit_base.len() != expected_base_len {
        bail!(
            "decode_edits: edit_base stream is {} bytes, expected {} for {} edits",
            edit_base.len(),
            expected_base_len,
            total
        );
    }
    // `total` is now proven == edit_base.len()*8/3 (the exact-length check above),
    // so this hint is bounded by input; cap it anyway for symmetry with steps 1-2.
    let mut bases = Vec::with_capacity(total.min(edit_base.len() + 1));
    for k in 0..total {
        let mut code = 0u8;
        for b in 0..3 {
            let bit_pos = 3 * k + b;
            let byte = edit_base[bit_pos / 8];
            let bit = (byte >> (bit_pos % 8)) & 1;
            code |= bit << b;
        }
        let base = code_to_base(code)
            .ok_or_else(|| anyhow::anyhow!("decode_edits: invalid 3-bit base code {}", code))?;
        bases.push(base);
    }

    // 4. Reassemble Vec<Vec<Edit>> using the per-read counts.
    let mut out = Vec::with_capacity(num_reads);
    let mut idx = 0usize;
    for &c in &counts {
        let mut rec = Vec::with_capacity(c);
        for _ in 0..c {
            rec.push(Edit {
                pos: positions[idx],
                base: bases[idx],
            });
            idx += 1;
        }
        out.push(rec);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a `PackedBacking` from an ASCII `ACGTN` string (test-only).
    fn backing_from_ascii(s: &[u8]) -> PackedBacking {
        let len = s.len();
        let mut packed = vec![0u8; len.div_ceil(4)];
        let mut nbm = vec![0u8; len.div_ceil(8)];
        for (i, &b) in s.iter().enumerate() {
            let code: u8 = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                b'N' => {
                    nbm[i / 8] |= 1 << (i % 8);
                    0
                }
                _ => panic!("backing_from_ascii: non-ACGTN byte {b}"),
            };
            packed[i / 4] |= code << (2 * (i % 4));
        }
        PackedBacking::from_parts(packed, nbm, len)
    }

    #[test]
    fn fused_reconstruct_matches_two_step() {
        // Sweep every (base, len, is_rc, edit-set) combination over a backing that
        // includes Ns, and assert the single-allocation fused path is byte-identical
        // to `reconstruct_read(&backing.slice(..), ..)`.
        let backing = backing_from_ascii(b"ACGTNACGTACGTNNACGT");
        let blen = backing.len();
        let edit_sets: &[&[Edit]] = &[
            &[],
            &[Edit { pos: 0, base: b'C' }],
            &[Edit { pos: 1, base: b'N' }, Edit { pos: 3, base: b'T' }],
        ];
        for base in 0..blen {
            for len in 1..=(blen - base) {
                for &is_rc in &[false, true] {
                    for eset in edit_sets {
                        // Keep only edits whose position is within this read length.
                        let edits: Vec<Edit> =
                            eset.iter().copied().filter(|e| (e.pos as usize) < len).collect();
                        let expected =
                            reconstruct_read(&backing.slice(base, len), &edits, is_rc);
                        let got =
                            reconstruct_read_from_backing(&backing, base, len, &edits, is_rc);
                        assert_eq!(
                            got, expected,
                            "fused != two-step at base={base} len={len} is_rc={is_rc} edits={edits:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn reconstruct_forward_and_reverse() {
        let consensus = b"ACGTAC";
        let edits = vec![Edit { pos: 2, base: b'T' }];
        assert_eq!(
            reconstruct_read(consensus, &edits, /*is_rc*/ false),
            b"ACTTAC".to_vec()
        );
        let r = reconstruct_read(consensus, &edits, true);
        assert_eq!(
            r,
            crate::compression::dna_utils::reverse_complement(b"ACTTAC")
        );
    }

    #[test]
    fn columnar_codec_roundtrip() {
        let recs = vec![vec![Edit { pos: 2, base: b'T' }], vec![]];
        let enc = encode_edits(&recs.iter().map(Vec::as_slice).collect::<Vec<_>>());
        let dec = decode_edits(&enc.edit_count, &enc.edit_pos, &enc.edit_base, recs.len()).unwrap();
        assert_eq!(dec, recs);
    }

    #[test]
    fn extract_edits_and_verify() {
        let consensus = b"ACGTACGT";
        let read = b"ACATACGT"; // pos2 G->A
        let e = extract_substitutions(consensus, read).unwrap();
        assert_eq!(e, vec![Edit { pos: 2, base: b'A' }]);
        assert_eq!(reconstruct_read(consensus, &e, false), read.to_vec());
    }

    #[test]
    fn three_bit_packing_crosses_byte_boundaries() {
        // First read has several edits whose bases mix A/C/G/T/N so the packed
        // base bitstream crosses byte boundaries (6 edits * 3 bits = 18 bits = 3
        // bytes), plus a second read.
        let recs = vec![
            vec![
                Edit { pos: 0, base: b'A' }, // 0
                Edit { pos: 1, base: b'C' }, // 1
                Edit { pos: 2, base: b'G' }, // 2
                Edit { pos: 3, base: b'T' }, // 3
                Edit { pos: 4, base: b'N' }, // 4
                Edit { pos: 5, base: b'T' }, // 3
            ],
            vec![Edit { pos: 7, base: b'N' }, Edit { pos: 9, base: b'A' }],
        ];
        let enc = encode_edits(&recs.iter().map(Vec::as_slice).collect::<Vec<_>>());
        // 8 edits total -> ceil(8*3/8) = 3 bytes.
        assert_eq!(enc.edit_base.len(), 3);
        let dec = decode_edits(&enc.edit_count, &enc.edit_pos, &enc.edit_base, recs.len()).unwrap();
        assert_eq!(dec, recs);
    }

    #[test]
    fn n_in_edits_reconstructs_and_roundtrips() {
        let consensus = b"ACGTACGT";
        let read = b"ACNTACGT"; // pos2 G->N
        let e = extract_substitutions(consensus, read).unwrap();
        assert_eq!(e, vec![Edit { pos: 2, base: b'N' }]);
        assert_eq!(reconstruct_read(consensus, &e, false), read.to_vec());

        let recs = vec![e];
        let enc = encode_edits(&recs.iter().map(Vec::as_slice).collect::<Vec<_>>());
        let dec = decode_edits(&enc.edit_count, &enc.edit_pos, &enc.edit_base, recs.len()).unwrap();
        assert_eq!(dec, recs);
    }

    #[test]
    fn extract_rejects_unequal_lengths() {
        assert_eq!(extract_substitutions(b"ACGT", b"ACG"), None);
        assert_eq!(extract_substitutions(b"ACG", b"ACGT"), None);
    }

    #[test]
    fn extract_rejects_non_acgtn_read_base() {
        // 'R' is a valid IUPAC ambiguity code but not ACGTN — must reject.
        assert_eq!(extract_substitutions(b"ACGT", b"ACRT"), None);
    }

    #[test]
    fn decode_rejects_count_position_mismatch() {
        // edit_count says 1 edit, but edit_pos has 0 positions (truncated).
        let mut edit_count = Vec::new();
        write_varint(&mut edit_count, 1);
        let edit_pos: Vec<u8> = Vec::new();
        let edit_base: Vec<u8> = Vec::new();
        assert!(decode_edits(&edit_count, &edit_pos, &edit_base, 1).is_err());
    }

    #[test]
    fn decode_rejects_truncated_edit_pos() {
        // count=1; edit_pos holds a varint with continuation bit but no follow byte.
        let mut edit_count = Vec::new();
        write_varint(&mut edit_count, 1);
        let edit_pos: Vec<u8> = vec![0x80]; // continuation bit set, truncated
        let edit_base: Vec<u8> = vec![0u8];
        assert!(decode_edits(&edit_count, &edit_pos, &edit_base, 1).is_err());
    }

    #[test]
    fn decode_rejects_hostile_edit_count_without_oom() {
        // A single read declaring ~10^18 edits via one huge varint, but tiny
        // edit_pos/edit_base. Pre-fix this hit `Vec::with_capacity(total)` (capacity
        // overflow / OOM) before the position loop could prove the count is a lie.
        // With the capacity hint bounded by `edit_pos.len()`, decode must instead
        // return `Err` quickly and the process must survive (reaching this assert at
        // all proves no allocation abort).
        let mut edit_count = Vec::new();
        write_varint(&mut edit_count, 1_000_000_000_000_000_000u64);
        let edit_pos: Vec<u8> = Vec::new();
        let edit_base: Vec<u8> = Vec::new();
        let err = decode_edits(&edit_count, &edit_pos, &edit_base, 1).unwrap_err();
        assert!(
            err.to_string().contains("edit_pos"),
            "expected truncated-edit_pos error, got: {err}"
        );
    }

    #[test]
    fn decode_rejects_bad_base_code() {
        // One edit; force the 3-bit base code to 7 (>4) in edit_base.
        let mut edit_count = Vec::new();
        write_varint(&mut edit_count, 1);
        let mut edit_pos = Vec::new();
        write_varint(&mut edit_pos, 0);
        // total=1 -> ceil(3/8)=1 byte; low 3 bits = 0b111 = 7 (invalid).
        let edit_base: Vec<u8> = vec![0b0000_0111];
        assert!(decode_edits(&edit_count, &edit_pos, &edit_base, 1).is_err());
    }

    #[test]
    fn decode_rejects_trailing_garbage_in_pos() {
        // count=1; edit_pos has one valid varint plus a trailing garbage byte.
        let mut edit_count = Vec::new();
        write_varint(&mut edit_count, 1);
        let mut edit_pos = Vec::new();
        write_varint(&mut edit_pos, 5);
        edit_pos.push(0x00); // trailing byte
        let edit_base: Vec<u8> = vec![0u8]; // code 0 = 'A'
        assert!(decode_edits(&edit_count, &edit_pos, &edit_base, 1).is_err());
    }

    #[test]
    fn decode_rejects_trailing_garbage_in_count() {
        // num_reads=1; edit_count has an extra trailing byte.
        let mut edit_count = Vec::new();
        write_varint(&mut edit_count, 0);
        edit_count.push(0x00); // trailing byte
        assert!(decode_edits(&edit_count, &[], &[], 1).is_err());
    }

    #[test]
    fn decode_rejects_wrong_base_stream_length() {
        // total=1 -> expects 1 byte; give 2 bytes (trailing full byte).
        let mut edit_count = Vec::new();
        write_varint(&mut edit_count, 1);
        let mut edit_pos = Vec::new();
        write_varint(&mut edit_pos, 0);
        let edit_base: Vec<u8> = vec![0u8, 0u8];
        assert!(decode_edits(&edit_count, &edit_pos, &edit_base, 1).is_err());
    }
}
