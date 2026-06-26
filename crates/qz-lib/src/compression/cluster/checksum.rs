// FNV-1a 64-bit constants (two lanes with distinct offset bases -> 128-bit, no dependency).
const FNV_PRIME: u64 = 0x0000_0100_0000_01B3;
const FNV_OFFSET_A: u64 = 0xCBF2_9CE4_8422_2325;
const FNV_OFFSET_B: u64 = 0x84222325CBF29CE4u64.rotate_left(7); // distinct basis for lane B

fn fnv1a_lane(mut h: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

/// Hash one length-prefixed, domain-tagged unit -> 128-bit. Single-end record = (header,seq,qual).
/// Length prefixes make concatenation unambiguous; domain tag 0 = single record.
pub(crate) fn hash_record(header: &[u8], seq: &[u8], qual: &[u8]) -> u128 {
    let mut a = FNV_OFFSET_A;
    let mut b = FNV_OFFSET_B;
    a = fnv1a_lane(a, &[0u8]); b = fnv1a_lane(b, &[0u8]); // domain tag = single
    for field in [header, seq, qual] {
        let len = (field.len() as u64).to_le_bytes();
        a = fnv1a_lane(a, &len); a = fnv1a_lane(a, field);
        b = fnv1a_lane(b, &len); b = fnv1a_lane(b, field);
    }
    ((a as u128) << 64) | (b as u128)
}

/// Hash one PAIR as a unit -> 128-bit. The six fields are length-prefixed in the fixed
/// order R1.header, R1.seq, R1.qual, R2.header, R2.seq, R2.qual, so an R1/R2 swap or an
/// R2-within-pair permutation changes the hash. Domain tag 1 distinguishes a pair from a
/// single record (tag 0), so single/paired checksums never collide.
pub(crate) fn hash_pair(
    r1_header: &[u8], r1_seq: &[u8], r1_qual: &[u8],
    r2_header: &[u8], r2_seq: &[u8], r2_qual: &[u8],
) -> u128 {
    let mut a = FNV_OFFSET_A;
    let mut b = FNV_OFFSET_B;
    a = fnv1a_lane(a, &[1u8]); b = fnv1a_lane(b, &[1u8]); // domain tag = pair
    for field in [r1_header, r1_seq, r1_qual, r2_header, r2_seq, r2_qual] {
        let len = (field.len() as u64).to_le_bytes();
        a = fnv1a_lane(a, &len); a = fnv1a_lane(a, field);
        b = fnv1a_lane(b, &len); b = fnv1a_lane(b, field);
    }
    ((a as u128) << 64) | (b as u128)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn checksum_is_deterministic_and_commutative() {
        let a = hash_record(b"@r1", b"ACGT", b"IIII");
        let b = hash_record(b"@r2", b"TTTT", b"####");
        // commutative accumulation: order-independent
        let s1 = 0u128.wrapping_add(a).wrapping_add(b);
        let s2 = 0u128.wrapping_add(b).wrapping_add(a);
        assert_eq!(s1, s2);
        // deterministic: same input -> same hash, no random seed
        assert_eq!(a, hash_record(b"@r1", b"ACGT", b"IIII"));
        // field framing: swapping seq/qual bytes must change the hash
        assert_ne!(hash_record(b"@r1", b"ACGT", b"IIII"),
                   hash_record(b"@r1", b"IIII", b"ACGT"));
    }

    #[test]
    fn hash_pair_detects_mate_swap_and_is_commutative() {
        let p1 = hash_pair(b"@a/1", b"ACGT", b"IIII", b"@a/2", b"TTTT", b"####");
        let p2 = hash_pair(b"@b/1", b"GGGG", b"JJJJ", b"@b/2", b"CCCC", b"$$$$");
        // commutative accumulation (order-independent set fold)
        assert_eq!(0u128.wrapping_add(p1).wrapping_add(p2),
                   0u128.wrapping_add(p2).wrapping_add(p1));
        // deterministic
        assert_eq!(p1, hash_pair(b"@a/1", b"ACGT", b"IIII", b"@a/2", b"TTTT", b"####"));
        // swapping R1<->R2 mates changes the hash (pairing/orientation matters)
        assert_ne!(p1, hash_pair(b"@a/2", b"TTTT", b"####", b"@a/1", b"ACGT", b"IIII"));
        // a single record and a pair must not collide (distinct domain tags)
        assert_ne!(hash_record(b"@a/1", b"ACGT", b"IIII"), p1);
    }
}
