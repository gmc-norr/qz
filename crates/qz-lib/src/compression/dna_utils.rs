//! Shared DNA compression utilities: reverse-complement / canonicalization,
//! k-mer hashing, open-syncmer cluster keys (cluster mode), and varint I/O used
//! across the sequence, header, reference, and cluster codecs.

/// Get reverse complement of a sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

/// In-place equivalent of [`reverse_complement`]: reverses `seq` and complements
/// each base, allocating nothing. Byte-for-byte identical to
/// `reverse_complement(seq)` for the same bytes (A/a→T, C/c→G, G/g→C, T/t→A,
/// every other byte → N).
pub fn reverse_complement_in_place(seq: &mut [u8]) {
    #[inline]
    fn complement(b: u8) -> u8 {
        match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        }
    }
    let mut lo = 0usize;
    let mut hi = seq.len();
    while lo < hi {
        hi -= 1;
        if lo == hi {
            // odd-length middle element
            seq[lo] = complement(seq[lo]);
            break;
        }
        let new_lo = complement(seq[hi]);
        let new_hi = complement(seq[lo]);
        seq[lo] = new_lo;
        seq[hi] = new_hi;
        lo += 1;
    }
}

/// Case-preserving **involution** complement used by the rc-canonical sequence
/// transform: `A<->T`, `C<->G` (and lowercase `a<->t`, `c<->g`); every other
/// byte (`N`, IUPAC ambiguity codes, soft-mask case it doesn't recognise) maps
/// to itself. Because it is a true involution, [`reverse_complement_canonical`]
/// is its own inverse, so the rc_canon encoding round-trips losslessly for
/// arbitrary sequence bytes — unlike [`reverse_complement`], which folds
/// lowercase to uppercase and non-ACGT to `N` (lossy, intended only for the
/// assembly/k-mer paths).
#[inline]
fn complement_canonical(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

/// Reverse-complement using the case-preserving involution [`complement_canonical`].
/// `reverse_complement_canonical(reverse_complement_canonical(x)) == x` for ALL
/// inputs, which is what makes the rc_canon transform lossless.
pub fn reverse_complement_canonical(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_canonical(b)).collect()
}

/// Canonicalize a sequence: return the lexicographically smaller of {seq, revcomp(seq)}.
/// Returns (canonical_sequence, was_reversed) where was_reversed=true means revcomp was chosen.
/// Uses early-exit comparison — stops at first differing base.
///
/// Uses [`reverse_complement_canonical`] (a lossless involution) so the stored
/// orientation can be reversed exactly on decode regardless of case / IUPAC
/// content.
pub fn canonicalize_sequence(seq: &[u8]) -> (Vec<u8>, bool) {
    let len = seq.len();
    // Compare seq vs revcomp lexicographically without materializing revcomp
    for i in 0..len {
        let fwd = seq[i];
        let rev_base = complement_canonical(seq[len - 1 - i]);
        if fwd < rev_base {
            // Forward is smaller — keep original
            return (seq.to_vec(), false);
        } else if fwd > rev_base {
            // Revcomp is smaller
            return (reverse_complement_canonical(seq), true);
        }
    }
    // Palindrome — keep original
    (seq.to_vec(), false)
}

/// Encode a k-mer as a 2-bit hash (4 bases → 8 bits, k=20 → 40 bits fits u64).
/// Returns None if the k-mer contains N or other ambiguous bases.
pub fn kmer_to_hash(seq: &[u8]) -> Option<u64> {
    let mut hash = 0u64;
    for &base in seq {
        let val = match base {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => return None,
        };
        hash = (hash << 2) | val;
    }
    Some(hash)
}

/// Reverse complement a k-mer hash (A↔T, C↔G in 2-bit: 0↔3, 1↔2)
pub fn reverse_complement_hash(hash: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    let mut h = hash;
    for _ in 0..k {
        rc = (rc << 2) | (3 - (h & 3));
        h >>= 2;
    }
    rc
}

/// Convert a 0-3 index back to a DNA base
pub fn idx_to_base(i: usize) -> u8 {
    match i {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

/// Convert a 2-bit k-mer hash back to a DNA byte sequence.
/// The first base in the original sequence corresponds to the most-significant bits.
pub fn hash_to_kmer(hash: u64, k: usize) -> Vec<u8> {
    let mut kmer = vec![0u8; k];
    let mut h = hash;
    for i in (0..k).rev() {
        kmer[i] = idx_to_base((h & 3) as usize);
        h >>= 2;
    }
    kmer
}

/// Reverse complement a DNA byte sequence.
pub fn reverse_complement_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            other => other,
        })
        .collect()
}

/// Write a variable-length integer (LEB128 encoding)
pub fn write_varint(output: &mut Vec<u8>, mut value: u64) {
    while value >= 128 {
        output.push((value & 0x7F) as u8 | 0x80);
        value >>= 7;
    }
    output.push(value as u8);
}

/// Read a variable-length integer (LEB128 encoding).
///
/// Returns `None` on a malformed varint — a shift past the width of `u64`, high
/// bits that wouldn't fit, or more than 10 continuation bytes — rather than
/// panicking (debug) / silently masking the high bits (release) on crafted
/// input. This is the workhorse reader for the untrusted reference/cluster/paired
/// decode paths, so it mirrors the guards in `mod.rs`'s `read_varint` and
/// `header_col`. Byte-identical on every legitimate varint; every call site
/// already maps `None` to a clean error.
pub fn read_varint(data: &[u8], offset: &mut usize) -> Option<u64> {
    let mut value = 0u64;
    let mut shift = 0u32;
    let mut bytes_read = 0u32;

    loop {
        if *offset >= data.len() {
            return None;
        }

        let byte = data[*offset];
        *offset += 1;
        bytes_read += 1;

        // Reject a shift at/past u64's width or high bits that wouldn't fit,
        // before doing the shift (the `||` short-circuits, so `<< shift` with
        // shift >= 64 never executes).
        let chunk = (byte & 0x7F) as u64;
        if shift >= u64::BITS || (chunk << shift) >> shift != chunk {
            return None;
        }
        value |= chunk << shift;

        if byte & 0x80 == 0 {
            return Some(value);
        }

        shift += 7;
        if bytes_read >= 10 {
            return None; // too many continuation bytes for a u64 LEB128
        }
    }
}

/// Compute a 1-byte sequence hint from the minimum canonical open syncmer hash.
///
/// Uses k=21, s=10 open syncmers. Returns the top 8 bits of the minimum
/// canonical hash. Reads shorter than k=21 return 0xFF.
pub fn compute_sequence_hint(seq: &[u8]) -> u8 {
    let k = 21usize;
    let s = 10usize;
    let (hash, _pos) = min_open_syncmer_hash_with_pos(seq, k, s);
    if hash == u64::MAX {
        0xFF
    } else {
        // k=21 => 42-bit hash in 2-bit encoding. Top 8 bits = shift right by 34.
        (hash >> 34) as u8
    }
}

/// Compute the minimum canonical open syncmer hash for a sequence.
/// Returns u64::MAX if the read is too short for k=21.
///
/// Note: open-syncmer selection (`ts=[0]`) is strand-asymmetric. This helper is
/// kept for older callers that want the historical "forward scan with canonical
/// k-mers" behavior. Cluster mode uses
/// [`canonicalize_open_syncmer_strand`] so a read and its reverse-complement land
/// on the same key.
pub fn compute_min_syncmer_hash(seq: &[u8]) -> u64 {
    let (hash, _pos) = min_open_syncmer_hash_with_pos(seq, 21, 10);
    hash
}

/// Like [`compute_min_syncmer_hash`] but also returns the start position of the
/// minimizing k=21 open syncmer. The position is an anchor shared by reads
/// covering the same locus, so sorting reads by the bases *following* it aligns
/// shifted overlaps. Returns `(u64::MAX, 0)` for reads too short for k=21.
pub fn compute_min_syncmer_hash_pos(seq: &[u8]) -> (u64, u16) {
    min_open_syncmer_hash_with_pos(seq, 21, 10)
}

/// Find the minimum canonical open syncmer hash and its position in a sequence.
/// Open syncmers use ts=[0] (smallest s-mer at position 0 of the k-mer window only).
fn min_open_syncmer_hash_with_pos(seq: &[u8], k: usize, s: usize) -> (u64, u16) {
    if seq.len() <= k {
        return (u64::MAX, 0);
    }
    let positions = syncmers::find_syncmers_pos(k, s, &[0], seq);
    let mut min_hash = u64::MAX;
    let mut min_pos = 0u16;
    for pos in positions {
        if pos + k > seq.len() {
            continue;
        }
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            let canon = fwd.min(rc);
            if canon < min_hash {
                min_hash = canon;
                min_pos = pos as u16;
            }
        }
    }
    (min_hash, min_pos)
}

/// Find the minimum raw open-syncmer hash and its position when scanning a
/// single concrete strand. Unlike [`min_open_syncmer_hash_with_pos`], this does
/// not canonicalize each k-mer internally. It is the building block for the
/// strand-symmetric cluster key below: scan the read and its full reverse
/// complement, then choose the smaller strand/key pair.
fn min_open_syncmer_hash_one_strand(seq: &[u8], k: usize, s: usize) -> (u64, u16) {
    if seq.len() <= k {
        return (u64::MAX, 0);
    }
    let positions = syncmers::find_syncmers_pos(k, s, &[0], seq);
    let mut min_hash = u64::MAX;
    let mut min_pos = 0u16;
    for pos in positions {
        if pos + k > seq.len() {
            continue;
        }
        if let Some(hash) = kmer_to_hash(&seq[pos..pos + k])
            && hash < min_hash {
                min_hash = hash;
                min_pos = pos as u16;
            }
    }
    (min_hash, min_pos)
}

/// Strand-symmetric open-syncmer key used by cluster mode.
///
/// Open syncmers are orientation-dependent: selecting only k-mers whose minimum
/// s-mer is at offset 0 means a k-mer selected on one strand maps to an end
/// syncmer on the reverse-complement strand. To make the cluster key truly
/// RC-invariant, scan both full orientations and choose the smaller raw
/// open-syncmer hash. Returns `(key, anchor_pos, flip)`, where `anchor_pos` is in
/// the chosen orientation and `flip=true` means the caller should store the
/// reverse-complemented read.
pub fn strand_symmetric_open_syncmer_hash_pos(seq: &[u8]) -> (u64, u16, bool) {
    let k = 21usize;
    let s = 10usize;
    if seq.len() <= k {
        return (u64::MAX, 0, false);
    }

    let (fwd_hash, fwd_pos) = min_open_syncmer_hash_one_strand(seq, k, s);
    let rc_seq = reverse_complement_canonical(seq);
    let (rc_hash, rc_pos) = min_open_syncmer_hash_one_strand(&rc_seq, k, s);

    let choose_rc = if rc_hash < fwd_hash {
        true
    } else if rc_hash > fwd_hash || rc_hash == u64::MAX {
        // rc strictly worse, or no syncmer signal at all (both == u64::MAX): keep forward.
        false
    } else {
        // Rare hash/key tie: use lexicographic canonicalization so exact reverse
        // complements still co-orient deterministically.
        rc_seq.as_slice() < seq
    };

    if choose_rc {
        (rc_hash, rc_pos, true)
    } else {
        (fwd_hash, fwd_pos, false)
    }
}

/// Canonicalize a read according to the strand-symmetric open-syncmer key.
/// Returns `(canonical_sequence, flip, key, anchor_pos)`.
pub fn canonicalize_open_syncmer_strand(seq: &[u8]) -> (Vec<u8>, bool, u64, u16) {
    let k = 21usize;
    let s = 10usize;
    if seq.len() <= k {
        return (seq.to_vec(), false, u64::MAX, 0);
    }

    let (fwd_hash, fwd_pos) = min_open_syncmer_hash_one_strand(seq, k, s);
    let rc_seq = reverse_complement_canonical(seq);
    let (rc_hash, rc_pos) = min_open_syncmer_hash_one_strand(&rc_seq, k, s);

    let choose_rc = if rc_hash < fwd_hash {
        true
    } else if rc_hash > fwd_hash || rc_hash == u64::MAX {
        // rc strictly worse, or no syncmer signal at all (both == u64::MAX): keep forward.
        false
    } else {
        rc_seq.as_slice() < seq
    };

    if choose_rc {
        (rc_seq, true, rc_hash, rc_pos)
    } else {
        (seq.to_vec(), false, fwd_hash, fwd_pos)
    }
}

/// Historical cross-strand canonicalization decision. Returns `true` if `seq`
/// should be reverse-complemented so that the minimum canonical k-mer found by
/// a forward open-syncmer scan reads in forward-canonical orientation.
///
/// Open-syncmer selection is itself strand-asymmetric, so this is not a fully
/// RC-invariant cluster key. New cluster-mode code should use
/// [`canonicalize_open_syncmer_strand`], which scans both full orientations.
pub fn canonical_minimizer_orientation(seq: &[u8]) -> bool {
    let k = 21usize;
    let s = 10usize;
    if seq.len() <= k {
        return false;
    }
    let positions = syncmers::find_syncmers_pos(k, s, &[0], seq);
    let mut min_canon = u64::MAX;
    let mut flip = false;
    for pos in positions {
        if pos + k > seq.len() {
            continue;
        }
        let kmer = &seq[pos..pos + k];
        if let Some(fwd) = kmer_to_hash(kmer) {
            let rc = reverse_complement_hash(fwd, k);
            let canon = fwd.min(rc);
            if canon < min_canon {
                min_canon = canon;
                flip = rc < fwd; // the RC strand produced the canonical k-mer
            }
        }
    }
    flip
}

/// Write a little-endian u64
pub fn write_u64(output: &mut Vec<u8>, value: u64) {
    output.extend_from_slice(&value.to_le_bytes());
}

/// Read a little-endian u64
pub fn read_u64(data: &[u8], offset: &mut usize) -> Option<u64> {
    if *offset + 8 > data.len() {
        return None;
    }
    let bytes: [u8; 8] = data[*offset..*offset + 8].try_into().ok()?;
    *offset += 8;
    Some(u64::from_le_bytes(bytes))
}

/// Zigzag-encode a signed value then LEB128 (reuses write_varint).
/// Zigzag maps small-magnitude signed values to small unsigned values so the
/// LEB128 stays compact for both signs. `wrapping_shl` avoids a debug-mode
/// overflow panic on `i64::MIN`.
pub fn write_svarint(out: &mut Vec<u8>, v: i64) {
    write_varint(out, (v.wrapping_shl(1) ^ (v >> 63)) as u64);
}

/// Read a zigzag+LEB128 signed value written by `write_svarint`.
pub fn read_svarint(data: &[u8], offset: &mut usize) -> Option<i64> {
    let u = read_varint(data, offset)?;
    Some(((u >> 1) as i64) ^ -((u & 1) as i64))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn syncmer_entry_points_handle_exactly_k_length_reads() {
        // `syncmers::find_syncmers_pos` asserts `seq.len() > k` (k=21). A read of EXACTLY
        // 21 bp must be guarded out by every open-syncmer entry point (returning the
        // no-signal sentinel), not passed through to trip that assert and panic a worker.
        // Regression for the `< k` vs `<= k` off-by-one (reached unfiltered from cluster
        // mode on adapter-trimmed / fixed-21-mer libraries).
        let r21: Vec<u8> = b"ACGT".iter().cycle().take(21).copied().collect();
        assert_eq!(r21.len(), 21);

        assert_eq!(compute_min_syncmer_hash(&r21), u64::MAX);
        assert_eq!(compute_min_syncmer_hash_pos(&r21), (u64::MAX, 0));
        assert_eq!(
            strand_symmetric_open_syncmer_hash_pos(&r21),
            (u64::MAX, 0, false)
        );
        assert_eq!(
            canonicalize_open_syncmer_strand(&r21),
            (r21.clone(), false, u64::MAX, 0)
        );
        assert!(!canonical_minimizer_orientation(&r21));

        // A 22-bp read (the first length > k) must still be processed without panicking.
        let r22: Vec<u8> = b"ACGTACGTACGTACGTACGTAC".to_vec();
        assert_eq!(r22.len(), 22);
        let _ = compute_min_syncmer_hash(&r22);
        let _ = strand_symmetric_open_syncmer_hash_pos(&r22);
        let _ = canonicalize_open_syncmer_strand(&r22);
        let _ = canonical_minimizer_orientation(&r22);
    }

    #[test]
    fn rc_in_place_matches_reverse_complement() {
        let cases: &[&[u8]] = &[
            b"",
            b"A",
            b"AC",
            b"ACG",
            b"ACGTN",
            b"NNACGTACGT",
            b"GATTACA",
            b"NNNNN",
        ];
        for &seq in cases {
            let mut buf = seq.to_vec();
            reverse_complement_in_place(&mut buf);
            assert_eq!(
                buf,
                reverse_complement(seq),
                "rc-in-place mismatch for {:?}",
                std::str::from_utf8(seq)
            );
        }
    }

    #[test]
    fn revcomp_canonical_is_involution_over_all_bytes() {
        // The rc_canon transform is only lossless if reverse_complement_canonical
        // is its own inverse for EVERY byte (lowercase soft-mask, IUPAC, N, junk).
        let seqs: &[&[u8]] = &[
            b"ACGT",
            b"acgtACGT",
            b"NRYSWKMBDHVnryswkmbdhv",
            b"acgtNNNNacgt",
            b"ACGTacgtRYSW....----",
            &[0u8, 1, 2, 200, 255, b'A', b'a', b'N'],
        ];
        for s in seqs {
            let once = reverse_complement_canonical(s);
            let twice = reverse_complement_canonical(&once);
            assert_eq!(&twice, s, "not an involution for {s:?}");
        }
    }

    #[test]
    fn canonicalize_roundtrips_lowercase_and_iupac() {
        // canonicalize_sequence + (revcomp_canonical if was_reversed) must recover
        // the exact original bytes — the lossless contract the rc_canon archive path
        // relies on. Previously lowercase folded to uppercase and IUPAC folded to N.
        let seqs: &[&[u8]] = &[
            b"acgtACGT",
            b"TTTTaaaa",
            b"NRYSWacgtKMBDHV",
            b"ggggcccc",
            b"ACGTacgtRYSWN",
        ];
        for s in seqs {
            let (canon, was_reversed) = canonicalize_sequence(s);
            let restored = if was_reversed {
                reverse_complement_canonical(&canon)
            } else {
                canon
            };
            assert_eq!(&restored, s, "rc_canon roundtrip lost data for {s:?}");
        }
    }

    #[test]
    fn strand_symmetric_open_syncmer_is_rc_invariant() {
        let mut seed = 0x1234_5678_9abc_def0u64;
        let bases = b"ACGT";

        for _ in 0..512 {
            let mut seq = Vec::with_capacity(151);
            for _ in 0..151 {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                seq.push(bases[((seed >> 32) & 3) as usize]);
            }

            let rc = reverse_complement_canonical(&seq);
            let (key, _pos, _flip) = strand_symmetric_open_syncmer_hash_pos(&seq);
            let (rc_key, _rc_pos, _rc_flip) = strand_symmetric_open_syncmer_hash_pos(&rc);
            assert_eq!(key, rc_key, "read and RC must share the same cluster key");

            let (canon, _was_rc, canon_key, canon_pos) = canonicalize_open_syncmer_strand(&seq);
            let (one_strand_key, one_strand_pos) = min_open_syncmer_hash_one_strand(&canon, 21, 10);
            assert_eq!(canon_key, one_strand_key);
            assert_eq!(canon_pos, one_strand_pos);
        }
    }

    #[test]
    fn zigzag_varint_roundtrips_signed() {
        for v in [0i64, 1, -1, 2, -2, i64::MAX, i64::MIN, 123456, -123456] {
            let mut b = Vec::new();
            write_svarint(&mut b, v);
            let mut o = 0;
            assert_eq!(read_svarint(&b, &mut o).unwrap(), v);
            assert_eq!(o, b.len());
        }
    }

    #[test]
    fn read_varint_roundtrips_every_width_including_u64_max() {
        // Byte-identity guard: the hostile-input hardening must not change the
        // result for any legitimate value, including the full 10-byte u64::MAX.
        for v in [0u64, 1, 127, 128, 300, 16384, u32::MAX as u64, 1 << 62, u64::MAX] {
            let mut b = Vec::new();
            write_varint(&mut b, v);
            let mut o = 0;
            assert_eq!(read_varint(&b, &mut o), Some(v), "value {v}");
            assert_eq!(o, b.len(), "consumed all bytes for {v}");
        }
    }

    #[test]
    fn read_varint_rejects_overlong_continuation_without_panic() {
        // Crafted untrusted bytes: 10+ continuation bytes drive the shift past
        // u64's width. Must return None (clean error), NOT panic in debug
        // ("attempt to shift left with overflow") or mask high bits in release.
        let hostile = [0x80u8; 12]; // all continuation, never terminates in <=10
        let mut o = 0;
        assert_eq!(read_varint(&hostile, &mut o), None);

        // A continuation byte at shift 63 that would overflow is also rejected.
        let mut overflow = vec![0x80u8; 9]; // shifts 0..=56, then shift=63
        overflow.push(0x7F); // 7 high bits at shift 63 don't fit in u64
        let mut o2 = 0;
        assert_eq!(read_varint(&overflow, &mut o2), None);
    }
}
