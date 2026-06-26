//! MD:Z derivation from the per-chunk consensus.
//!
//! The MD tag is ~57% of a BAM archive's aux stream, and its content is the
//! reference bases each read covers — redundant with the consensus bz already
//! stores. This module provides the byte-exact primitives to drop MD on compress
//! and regenerate it on decompress:
//!
//!   - [`recover_ref`]: reference bases (ref order) from MD + CIGAR + SEQ.
//!   - [`regen_md`]: MD string from reference bases + CIGAR + SEQ (SAM-spec exact).
//!   - [`strip_md`] / [`insert_md`]: remove/re-insert the MD tag in an aux blob,
//!     byte-exact, preserving its tag position.
//!
//! Feasibility verified end-to-end on 5M real reads (HG002 chr20 + NovaSeq):
//! 100% byte-exact regeneration via consensus+exceptions, 0 conflicts; the
//! ref≠consensus exception map shrinks MD ~99%.

use super::streams::{read_varint, write_varint, zigzag_decode, zigzag_encode};
use anyhow::Result;
use std::collections::{BTreeMap, HashMap};

/// Build `mate[i]` = the in-chunk mate index of record `i` (or `u32::MAX`).
/// A pair is two **primary, paired** reads (FLAG & 0x900 == 0, FLAG & 0x1 set)
/// that share a QNAME and have opposite read1/read2 bits (FLAG & 0x40). This is
/// deterministic from `(names, flags)` alone, so compress and decode rebuild the
/// identical map without storing it. `names[i]` is the raw QNAME bytes (the same
/// bytes — trailing NUL included — on both sides).
pub fn build_mate_map(names: &[&[u8]], flags: &[u16]) -> Vec<u32> {
    let n = names.len();
    let mut mate = vec![u32::MAX; n];
    // (count capped at 3, first index, second index)
    let mut map: HashMap<&[u8], (u8, u32, u32)> = HashMap::with_capacity(n / 2 + 1);
    for i in 0..n {
        let f = flags[i];
        if f & 0x900 != 0 || f & 0x1 == 0 {
            continue; // secondary/supplementary, or not paired
        }
        let e = map.entry(names[i]).or_insert((0, u32::MAX, u32::MAX));
        match e.0 {
            0 => {
                e.1 = i as u32;
                e.0 = 1;
            }
            1 => {
                e.2 = i as u32;
                e.0 = 2;
            }
            _ => e.0 = 3, // >2 with this name: ambiguous, don't pair
        }
    }
    for &(cnt, a, b) in map.values() {
        if cnt == 2 && (flags[a as usize] & 0x40) != (flags[b as usize] & 0x40) {
            mate[a as usize] = b;
            mate[b as usize] = a;
        }
    }
    mate
}

/// Derive the signed TLEN of a read from its own and its mate's mapped extent.
/// Uses the SAM convention (leftmost +, rightmost −) with read1 → + on a tie.
/// Matches the source aligner only for some tie cases; callers verify against the
/// stored value and keep an exception otherwise.
pub fn derive_tlen(
    my_pos: i32,
    my_cigar: &[(u8, u32)],
    mate_pos: i32,
    mate_cigar: &[(u8, u32)],
    read1: bool,
) -> i32 {
    // Compute in i64 so an untrusted/crafted CIGAR with a huge ref-span can't overflow
    // `my_pos + span` (a debug/overflow-checks panic). TLEN is a 32-bit signed SAM field;
    // clamp the result into i32 range. For valid inputs no overflow occurs and the clamp
    // is a no-op, so the derived value is identical; the caller verifies it against the
    // stored TLEN and keeps an exception on any mismatch.
    let sa = my_pos as i64;
    let ea = my_pos as i64 + ref_span(my_cigar) as i64;
    let sb = mate_pos as i64;
    let eb = mate_pos as i64 + ref_span(mate_cigar) as i64;
    let span = ea.max(eb) - sa.min(sb);
    let tlen = if sa < sb {
        span
    } else if sa > sb {
        -span
    } else if read1 {
        span
    } else {
        -span
    };
    tlen.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}

/// Serialize the (ref_id, pos) -> ref-base exception map: varint count, then per
/// entry (sorted) zigzag-delta ref_id, zigzag-delta pos, 1 base byte.
pub fn serialize_exceptions(map: &BTreeMap<(i32, i32), u8>) -> Vec<u8> {
    let mut out = Vec::with_capacity(map.len() * 3 + 4);
    write_varint(&mut out, map.len());
    let (mut pr, mut pp) = (0i32, 0i32);
    for (&(rid, pos), &base) in map {
        write_varint(&mut out, zigzag_encode(rid.wrapping_sub(pr)) as usize);
        let dp = if rid == pr { pos.wrapping_sub(pp) } else { pos };
        write_varint(&mut out, zigzag_encode(dp) as usize);
        out.push(base);
        pr = rid;
        pp = pos;
    }
    out
}

/// Inverse of [`serialize_exceptions`].
pub fn deserialize_exceptions(data: &[u8]) -> Result<BTreeMap<(i32, i32), u8>> {
    let mut off = 0usize;
    let count = read_varint(data, &mut off)?;
    let mut map = BTreeMap::new();
    let (mut pr, mut pp) = (0i32, 0i32);
    for _ in 0..count {
        let d_ref = zigzag_decode(read_varint(data, &mut off)? as u32);
        let rid = pr.wrapping_add(d_ref);
        let d_pos = zigzag_decode(read_varint(data, &mut off)? as u32);
        let pos = if rid == pr {
            pp.wrapping_add(d_pos)
        } else {
            d_pos
        };
        if off >= data.len() {
            anyhow::bail!("md exceptions: truncated base at entry");
        }
        let base = data[off];
        off += 1;
        map.insert((rid, pos), base);
        pr = rid;
        pp = pos;
    }
    Ok(map)
}

const CIG_M: u8 = 0;
const CIG_I: u8 = 1;
const CIG_D: u8 = 2;
const CIG_N: u8 = 3;
const CIG_S: u8 = 4;
const CIG_EQ: u8 = 7;
const CIG_X: u8 = 8;

#[inline]
fn consumes_ref(op: u8) -> bool {
    matches!(op, CIG_M | CIG_D | CIG_N | CIG_EQ | CIG_X)
}

/// Number of reference positions a read spans (sum of ref-consuming CIGAR ops).
pub fn ref_span(cigar: &[(u8, u32)]) -> usize {
    cigar
        .iter()
        .filter(|(op, _)| consumes_ref(*op))
        .map(|(_, l)| *l as usize)
        .sum()
}

/// Recover the reference base at each ref-consuming position (ref order, starting
/// at the read's POS) by walking MD + CIGAR + SEQ together. `seq` is the read
/// sequence as ASCII bases (A/C/G/T/N), one per query base. Returns `None` if MD
/// and CIGAR disagree (caller keeps that read's MD verbatim).
#[allow(unused_assignments)] // match_left reset after a trailing deletion op
pub fn recover_ref(cigar: &[(u8, u32)], md: &[u8], seq: &[u8]) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(ref_span(cigar));
    let mut read_idx = 0usize;
    let mut mi = 0usize;
    let mut match_left: i64 = -1; // -1 = read a number next

    fn next_num(md: &[u8], mi: &mut usize) -> i64 {
        let mut v: i64 = 0;
        let mut any = false;
        while *mi < md.len() && md[*mi].is_ascii_digit() {
            v = v * 10 + (md[*mi] - b'0') as i64;
            *mi += 1;
            any = true;
        }
        if any { v } else { -1 }
    }

    for &(op, len) in cigar {
        let len = len as usize;
        match op {
            CIG_M | CIG_EQ | CIG_X => {
                for _ in 0..len {
                    if match_left < 0 {
                        match_left = next_num(md, &mut mi).max(0);
                    }
                    if match_left > 0 {
                        if read_idx >= seq.len() {
                            return None;
                        }
                        out.push(seq[read_idx]);
                        match_left -= 1;
                    } else {
                        if mi >= md.len() || !md[mi].is_ascii_alphabetic() {
                            return None;
                        }
                        out.push(md[mi]);
                        mi += 1;
                        match_left = -1;
                    }
                    read_idx += 1;
                }
            }
            CIG_I | CIG_S => read_idx += len,
            CIG_D => {
                if match_left < 0 {
                    match_left = next_num(md, &mut mi).max(0);
                }
                // a deletion must sit at a match boundary (match_left == 0)
                if mi >= md.len() || md[mi] != b'^' {
                    return None;
                }
                mi += 1;
                for _ in 0..len {
                    if mi >= md.len() || !md[mi].is_ascii_alphabetic() {
                        return None;
                    }
                    out.push(md[mi]);
                    mi += 1;
                }
                match_left = -1;
            }
            CIG_N => {
                // ref skip (intron): not represented in MD. Emit N placeholders;
                // regen skips N ops too, so this round-trips for DNA (no N ops).
                out.extend(std::iter::repeat_n(b'N', len));
            }
            _ => {} // H, P: no ref/query consumption
        }
    }
    Some(out)
}

/// Regenerate the MD string from reference bases (ref order) + CIGAR + SEQ.
/// Matches the SAM-spec form `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`.
pub fn regen_md(cigar: &[(u8, u32)], ref_bases: &[u8], seq: &[u8]) -> Vec<u8> {
    let mut md = Vec::new();
    let mut run: u64 = 0;
    let mut read_idx = 0usize;
    let mut ref_idx = 0usize;
    let push_num = |md: &mut Vec<u8>, n: u64| {
        let mut buf = [0u8; 20];
        let mut i = buf.len();
        let mut x = n;
        loop {
            i -= 1;
            buf[i] = b'0' + (x % 10) as u8;
            x /= 10;
            if x == 0 {
                break;
            }
        }
        md.extend_from_slice(&buf[i..]);
    };
    for &(op, len) in cigar {
        let len = len as usize;
        match op {
            CIG_M | CIG_EQ | CIG_X => {
                for _ in 0..len {
                    let rb = ref_bases[ref_idx];
                    if seq[read_idx] == rb {
                        run += 1;
                    } else {
                        push_num(&mut md, run);
                        run = 0;
                        md.push(rb);
                    }
                    read_idx += 1;
                    ref_idx += 1;
                }
            }
            CIG_I | CIG_S => read_idx += len,
            CIG_D => {
                push_num(&mut md, run);
                run = 0;
                md.push(b'^');
                for _ in 0..len {
                    md.push(ref_bases[ref_idx]);
                    ref_idx += 1;
                }
            }
            CIG_N => ref_idx += len,
            _ => {}
        }
    }
    push_num(&mut md, run);
    md
}

/// Length in bytes of one BAM aux item starting at `a[0]`, or `None` if malformed.
fn aux_item_len(a: &[u8]) -> Option<usize> {
    if a.len() < 3 {
        return None;
    }
    let ty = a[2];
    let vlen = match ty {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'Z' | b'H' => a[3..].iter().position(|&b| b == 0)? + 1,
        b'B' => {
            if a.len() < 8 {
                return None;
            }
            let sub = a[3];
            let cnt = u32::from_le_bytes([a[4], a[5], a[6], a[7]]) as usize;
            let esz = match sub {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            5 + cnt * esz
        }
        _ => return None,
    };
    Some(3 + vlen)
}

/// Split an aux blob into item byte-ranges. Returns `None` if malformed.
fn aux_items(aux: &[u8]) -> Option<Vec<(usize, usize)>> {
    let mut items = Vec::new();
    let mut off = 0;
    while off < aux.len() {
        let len = aux_item_len(&aux[off..])?;
        if off + len > aux.len() {
            return None;
        }
        items.push((off, len));
        off += len;
    }
    Some(items)
}

/// Result of stripping one Z-type tag from an aux blob.
pub struct StrippedTag {
    /// aux blob with the tag item removed (byte-exact remainder).
    pub aux: Vec<u8>,
    /// The tag value bytes (without tag/type/NUL): e.g. `b"10A5"`.
    pub value: Vec<u8>,
    /// 0-based item index where the tag was (for byte-exact re-insertion).
    pub index: u32,
}

/// Remove a two-letter Z-type tag (e.g. `b"MD"` or `b"MC"`) from an aux blob.
/// Returns `None` if the tag is absent or the aux is malformed (caller stores
/// the aux verbatim then). To strip MD and MC from the same read, strip MD
/// first, then MC from the result, and re-insert in the reverse order (MC then
/// MD) so each index stays valid for its array state.
pub fn strip_tag(aux: &[u8], tag: [u8; 2]) -> Option<StrippedTag> {
    let items = aux_items(aux)?;
    for (idx, &(off, len)) in items.iter().enumerate() {
        if aux[off] == tag[0] && aux[off + 1] == tag[1] && aux[off + 2] == b'Z' {
            let value = aux[off + 3..off + len - 1].to_vec(); // exclude NUL
            let mut out = Vec::with_capacity(aux.len() - len);
            out.extend_from_slice(&aux[..off]);
            out.extend_from_slice(&aux[off + len..]);
            return Some(StrippedTag {
                aux: out,
                value,
                index: idx as u32,
            });
        }
    }
    None
}

/// Re-insert a two-letter Z-type tag (`value`) at item `index` into a stripped
/// aux blob, byte-exact. Returns `None` if malformed or the index is out of range.
pub fn insert_tag(stripped: &[u8], tag: [u8; 2], value: &[u8], index: u32) -> Option<Vec<u8>> {
    let items = aux_items(stripped)?;
    let index = index as usize;
    if index > items.len() {
        return None;
    }
    let insert_at = if index < items.len() {
        items[index].0
    } else {
        stripped.len()
    };
    let mut out = Vec::with_capacity(stripped.len() + value.len() + 4);
    out.extend_from_slice(&stripped[..insert_at]);
    out.push(tag[0]);
    out.push(tag[1]);
    out.push(b'Z');
    out.extend_from_slice(value);
    out.push(0);
    out.extend_from_slice(&stripped[insert_at..]);
    Some(out)
}

/// Render a CIGAR op list to its textual form (`b"150M"`, `b"10S140M"`, ...) for
/// MC:Z regeneration. Empty CIGAR renders as `b"*"`.
pub fn cigar_text(cigar: &[(u8, u32)]) -> Vec<u8> {
    const OPS: &[u8; 9] = b"MIDNSHP=X";
    if cigar.is_empty() {
        return vec![b'*'];
    }
    let mut out = Vec::with_capacity(cigar.len() * 4);
    for &(op, len) in cigar {
        out.extend_from_slice(len.to_string().as_bytes());
        out.push(OPS[(op as usize).min(8)]);
    }
    out
}

// ── NM:i derivation (v11) ────────────────────────────────────────────────────
// NM (edit distance) = (# substitutions) + (inserted bases) + (deleted bases).
// It's fully derivable from the recovered reference + CIGAR + SEQ — the same
// `ref_bases` already used for MD regeneration. The aux integer *type byte*
// (c/C/s/S/i/I) is aligner-specific (bwa-mem2 writes unsigned C, some writers
// signed c), so it's preserved per stripped tag rather than re-derived.

/// Derive NM from the recovered reference (ref order, ref-consuming positions),
/// CIGAR, and query SEQ (ASCII). Returns `None` if the walk runs past either
/// buffer (caller keeps NM verbatim). Mirrors recover_ref/regen_md's walk so the
/// substitution count equals the MD mismatch-letter count exactly.
pub fn derive_nm(cigar: &[(u8, u32)], ref_bases: &[u8], seq: &[u8]) -> Option<i64> {
    let mut subs: i64 = 0;
    let mut ins: i64 = 0;
    let mut del: i64 = 0;
    let mut ref_idx = 0usize;
    let mut read_idx = 0usize;
    for &(op, len) in cigar {
        let len = len as usize;
        match op {
            CIG_M | CIG_EQ | CIG_X => {
                for _ in 0..len {
                    if read_idx >= seq.len() || ref_idx >= ref_bases.len() {
                        return None;
                    }
                    if seq[read_idx] != ref_bases[ref_idx] {
                        subs += 1;
                    }
                    read_idx += 1;
                    ref_idx += 1;
                }
            }
            CIG_I => {
                ins += len as i64;
                read_idx += len;
            }
            CIG_S => read_idx += len,
            CIG_D => {
                del += len as i64;
                ref_idx += len;
            }
            CIG_N => ref_idx += len,
            _ => {}
        }
    }
    Some(subs + ins + del)
}

/// Strip a 2-letter integer aux tag (type c/C/s/S/i/I) from an aux blob.
/// Returns `(type_byte, value, item_index, stripped_aux)`, or `None` if the tag
/// is absent / not an integer type / the aux is malformed.
pub fn strip_int_tag(aux: &[u8], tag: [u8; 2]) -> Option<(u8, i64, u32, Vec<u8>)> {
    let items = aux_items(aux)?;
    for (idx, &(off, len)) in items.iter().enumerate() {
        if aux[off] == tag[0] && aux[off + 1] == tag[1] {
            let ty = aux[off + 2];
            let v: i64 = match ty {
                b'c' => (aux[off + 3] as i8) as i64,
                b'C' => aux[off + 3] as i64,
                b's' => i16::from_le_bytes([aux[off + 3], aux[off + 4]]) as i64,
                b'S' => u16::from_le_bytes([aux[off + 3], aux[off + 4]]) as i64,
                b'i' => {
                    i32::from_le_bytes([aux[off + 3], aux[off + 4], aux[off + 5], aux[off + 6]])
                        as i64
                }
                b'I' => {
                    u32::from_le_bytes([aux[off + 3], aux[off + 4], aux[off + 5], aux[off + 6]])
                        as i64
                }
                _ => return None, // not an integer tag
            };
            let mut out = Vec::with_capacity(aux.len() - len);
            out.extend_from_slice(&aux[..off]);
            out.extend_from_slice(&aux[off + len..]);
            return Some((ty, v, idx as u32, out));
        }
    }
    None
}

/// Encode an integer aux item `[tag0, tag1, type, value-bytes]` in the given
/// type's width. Returns `None` if the value doesn't fit (never happens when the
/// derived value equals the stored one, which the caller verifies).
pub fn encode_int_item(tag: [u8; 2], ty: u8, v: i64) -> Option<Vec<u8>> {
    let mut out = vec![tag[0], tag[1], ty];
    match ty {
        b'c' => {
            if !(-128..=127).contains(&v) {
                return None;
            }
            out.push(v as i8 as u8);
        }
        b'C' => {
            if !(0..=255).contains(&v) {
                return None;
            }
            out.push(v as u8);
        }
        b's' => {
            if !(-32768..=32767).contains(&v) {
                return None;
            }
            out.extend_from_slice(&(v as i16).to_le_bytes());
        }
        b'S' => {
            if !(0..=65535).contains(&v) {
                return None;
            }
            out.extend_from_slice(&(v as u16).to_le_bytes());
        }
        b'i' => {
            if v < i32::MIN as i64 || v > i32::MAX as i64 {
                return None;
            }
            out.extend_from_slice(&(v as i32).to_le_bytes());
        }
        b'I' => {
            if !(0..=u32::MAX as i64).contains(&v) {
                return None;
            }
            out.extend_from_slice(&(v as u32).to_le_bytes());
        }
        _ => return None,
    }
    Some(out)
}

/// Insert raw aux-item bytes at item position `index` into a stripped aux blob,
/// byte-exact. Returns `None` if malformed or the index is out of range.
pub fn insert_item(stripped: &[u8], item: &[u8], index: u32) -> Option<Vec<u8>> {
    let items = aux_items(stripped)?;
    let index = index as usize;
    if index > items.len() {
        return None;
    }
    let at = if index < items.len() {
        items[index].0
    } else {
        stripped.len()
    };
    let mut out = Vec::with_capacity(stripped.len() + item.len());
    out.extend_from_slice(&stripped[..at]);
    out.extend_from_slice(item);
    out.extend_from_slice(&stripped[at..]);
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn non_canonical_md_fails_regen_roundtrip() {
        // The compress-side MD self-check only strips MD when
        // `regen_md(cigar, recover_ref(cigar, md, seq), seq) == md`. This guards
        // against non-canonical MD that `recover_ref` accepts but `regen_md`
        // normalizes differently — which would otherwise yield an *undecodable*
        // archive (content-CRC mismatch at decompress).
        let cigar = vec![(CIG_M, 8u32)];
        // MD "4A3" claims a mismatch at pos 4, but the read base there is also 'A',
        // so regen_md collapses the whole span to "8" — NOT reproducible.
        let md = b"4A3";
        let seq_bad = b"AAAAAAAA";
        let rb = recover_ref(&cigar, md, seq_bad).expect("recover_ref accepts it");
        assert_ne!(
            regen_md(&cigar, &rb, seq_bad),
            md,
            "non-canonical MD must NOT round-trip (self-check keeps it verbatim)"
        );
        // A genuine mismatch (read 'C' vs ref 'A') is canonical and round-trips.
        let seq_ok = b"AAAACAAA";
        let rb_ok = recover_ref(&cigar, md, seq_ok).expect("recover_ref");
        assert_eq!(
            regen_md(&cigar, &rb_ok, seq_ok),
            md,
            "canonical MD round-trips"
        );
    }

    // CIGAR helpers
    fn m(n: u32) -> (u8, u32) {
        (CIG_M, n)
    }
    fn d(n: u32) -> (u8, u32) {
        (CIG_D, n)
    }
    fn i(n: u32) -> (u8, u32) {
        (CIG_I, n)
    }
    fn s(n: u32) -> (u8, u32) {
        (CIG_S, n)
    }

    /// Recover ref then regen MD must reproduce the input MD byte-exact.
    fn roundtrip(cigar: &[(u8, u32)], md: &[u8], seq: &[u8]) {
        let refb = recover_ref(cigar, md, seq).expect("recover_ref");
        let regen = regen_md(cigar, &refb, seq);
        assert_eq!(
            regen,
            md,
            "MD regen mismatch: got {:?} want {:?}",
            String::from_utf8_lossy(&regen),
            String::from_utf8_lossy(md)
        );
    }

    #[test]
    fn all_match() {
        roundtrip(&[m(8)], b"8", b"ACGTACGT");
    }

    #[test]
    fn single_mismatch() {
        // ref T at position 3 (read has A there)
        roundtrip(&[m(8)], b"3T4", b"ACGAACGT");
    }

    #[test]
    fn adjacent_mismatches_zero_run() {
        // two adjacent mismatches -> a 0 between them. read bases (C,A) differ
        // from ref letters (A,C); positions 2..6 (GTAC) match.
        roundtrip(&[m(6)], b"0A0C4", b"CAGTAC");
    }

    #[test]
    fn deletion() {
        // 4M 2D 4M: deletion of ref "GT"
        roundtrip(&[m(4), d(2), m(4)], b"4^GT4", b"ACGTACGT");
    }

    #[test]
    fn deletion_then_mismatch() {
        // after the deletion, first ref base is A but read has T (mismatch).
        roundtrip(&[m(4), d(2), m(4)], b"4^GT0A3", b"ACGTTCGT");
    }

    #[test]
    fn insertion_is_skipped_by_md() {
        // 4M 2I 4M: insertion consumes read but not ref; MD spans 8 ref bases
        roundtrip(&[m(4), i(2), m(4)], b"8", b"ACGTNNACGT");
    }

    #[test]
    fn soft_clip_skipped() {
        roundtrip(&[s(2), m(6)], b"6", b"NNACGTAC");
    }

    #[test]
    fn leading_and_trailing_zero() {
        // mismatch at first base (ref A, read C) -> leading 0; mismatch at last
        // (ref C, read A) -> trailing 0; middle GT match.
        roundtrip(&[m(4)], b"0A2C0", b"CGTA");
    }

    #[test]
    fn strip_insert_roundtrip() {
        // aux: NM:C:1  MD:Z:3T4  AS:C:0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC");
        aux.push(1);
        aux.extend_from_slice(b"MDZ");
        aux.extend_from_slice(b"3T4");
        aux.push(0);
        aux.extend_from_slice(b"ASC");
        aux.push(0);
        let st = strip_tag(&aux, *b"MD").expect("has MD");
        assert_eq!(st.value, b"3T4");
        assert_eq!(st.index, 1);
        let rebuilt = insert_tag(&st.aux, *b"MD", &st.value, st.index).unwrap();
        assert_eq!(rebuilt, aux, "aux re-insertion not byte-exact");
    }

    #[test]
    fn strip_tag_first_and_last_position() {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MDZ");
        aux.extend_from_slice(b"10");
        aux.push(0);
        aux.extend_from_slice(b"NMC");
        aux.push(0);
        let st = strip_tag(&aux, *b"MD").unwrap();
        assert_eq!(st.index, 0);
        assert_eq!(
            insert_tag(&st.aux, *b"MD", &st.value, st.index).unwrap(),
            aux
        );
        let mut aux2 = Vec::new();
        aux2.extend_from_slice(b"NMC");
        aux2.push(0);
        aux2.extend_from_slice(b"MDZ");
        aux2.extend_from_slice(b"10");
        aux2.push(0);
        let st2 = strip_tag(&aux2, *b"MD").unwrap();
        assert_eq!(st2.index, 1);
        assert_eq!(
            insert_tag(&st2.aux, *b"MD", &st2.value, st2.index).unwrap(),
            aux2
        );
    }

    #[test]
    fn strip_md_and_mc_compose() {
        // NM:C:1  MD:Z:3T4  MC:Z:150M  AS:C:0 — strip MD then MC; insert MC then MD.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC");
        aux.push(1);
        aux.extend_from_slice(b"MDZ");
        aux.extend_from_slice(b"3T4");
        aux.push(0);
        aux.extend_from_slice(b"MCZ");
        aux.extend_from_slice(b"150M");
        aux.push(0);
        aux.extend_from_slice(b"ASC");
        aux.push(0);
        let md = strip_tag(&aux, *b"MD").unwrap();
        assert_eq!(md.index, 1);
        let mc = strip_tag(&md.aux, *b"MC").unwrap();
        assert_eq!(mc.index, 1);
        let step = insert_tag(&mc.aux, *b"MC", &mc.value, mc.index).unwrap();
        let full = insert_tag(&step, *b"MD", &md.value, md.index).unwrap();
        assert_eq!(full, aux, "MD+MC compose not byte-exact");
    }

    #[test]
    fn cigar_text_rendering() {
        assert_eq!(cigar_text(&[(0, 150)]), b"150M");
        assert_eq!(cigar_text(&[(4, 10), (0, 140)]), b"10S140M");
        assert_eq!(cigar_text(&[(0, 50), (2, 2), (0, 98)]), b"50M2D98M");
        assert_eq!(cigar_text(&[]), b"*");
    }

    #[test]
    fn exceptions_roundtrip() {
        let mut m = BTreeMap::new();
        m.insert((0i32, 100i32), b'A');
        m.insert((0, 105), b'C');
        m.insert((0, 60_000_000), b'G');
        m.insert((3, 7), b'T');
        let ser = serialize_exceptions(&m);
        let de = deserialize_exceptions(&ser).unwrap();
        assert_eq!(m, de);
    }

    #[test]
    fn no_tag_returns_none() {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NMC");
        aux.push(0);
        assert!(strip_tag(&aux, *b"MD").is_none());
    }

    // ── NM:i derivation ──────────────────────────────────────────────────────

    #[test]
    fn derive_nm_substitutions_and_indels() {
        // ref = recover_ref output (ref order). 8M with 1 mismatch at idx 3.
        let refb = recover_ref(&[m(8)], b"3T4", b"ACGAACGT").unwrap();
        assert_eq!(derive_nm(&[m(8)], &refb, b"ACGAACGT"), Some(1));
        // 4M 2I 4M: insertion adds 2 to NM; 1 substitution => NM 3.
        let cig = [m(4), i(2), m(4)];
        let refb = recover_ref(&cig, b"3A4", b"ACGTAAACGT").unwrap();
        assert_eq!(derive_nm(&cig, &refb, b"ACGTAAACGT"), Some(3));
        // 4M 2D 4M: deletion adds 2 to NM; no subs => NM 2.
        let cig = [m(4), d(2), m(4)];
        let refb = recover_ref(&cig, b"4^GG4", b"ACGTACGT").unwrap();
        assert_eq!(derive_nm(&cig, &refb, b"ACGTACGT"), Some(2));
    }

    /// strip an integer tag then re-insert the (type, value) byte-exact.
    fn int_tag_roundtrip(ty: u8, v: i64, val_bytes: &[u8]) {
        // aux = [XX <other Z>] [NM <ty> <val>] [YY <Z>]
        let mut aux = Vec::new();
        aux.extend_from_slice(b"XXZhi\0");
        aux.extend_from_slice(b"NM");
        aux.push(ty);
        aux.extend_from_slice(val_bytes);
        aux.extend_from_slice(b"YYZbye\0");
        let orig = aux.clone();
        let (sty, sval, idx, stripped) = strip_int_tag(&aux, *b"NM").expect("strip NM");
        assert_eq!(sty, ty);
        assert_eq!(sval, v);
        assert_eq!(idx, 1);
        let item = encode_int_item(*b"NM", sty, sval).expect("encode NM");
        let rebuilt = insert_item(&stripped, &item, idx).expect("insert NM");
        assert_eq!(
            rebuilt, orig,
            "NM int-tag roundtrip mismatch for type {}",
            ty as char
        );
    }

    #[test]
    fn nm_int_tag_roundtrip_all_widths() {
        int_tag_roundtrip(b'C', 3, &[3]);
        int_tag_roundtrip(b'c', 3, &[3]); // signed, same byte for small positives
        int_tag_roundtrip(b'C', 200, &[200]);
        int_tag_roundtrip(b'S', 300, &300u16.to_le_bytes());
        int_tag_roundtrip(b's', 300, &300i16.to_le_bytes());
        int_tag_roundtrip(b'i', 100_000, &100_000i32.to_le_bytes());
        int_tag_roundtrip(b'I', 100_000, &100_000u32.to_le_bytes());
    }

    #[test]
    fn encode_int_item_rejects_overflow() {
        assert!(encode_int_item(*b"NM", b'C', 256).is_none()); // > u8::MAX
        assert!(encode_int_item(*b"NM", b'c', 200).is_none()); // > i8::MAX
        assert!(encode_int_item(*b"NM", b'C', 5).is_some());
    }
}
