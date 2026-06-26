/// FFI bindings to htscodecs fqzcomp_qual — context-modeled quality score compression
///
/// fqzcomp uses an adaptive arithmetic coder with a rich context model
/// (position, previous quality, cumulative delta) to compress quality scores.
use anyhow::Result;

// FFI declarations matching fqzcomp_qual.h
#[repr(C)]
struct FqzSlice {
    num_records: libc::c_int,
    len: *const u32,
    flags: *const u32,
}

#[link(name = "htscodecs", kind = "static")]
unsafe extern "C" {
    fn fqz_compress(
        vers: libc::c_int,
        s: *const FqzSlice,
        input: *const libc::c_char,
        in_size: libc::size_t,
        out_size: *mut libc::size_t,
        strat: libc::c_int,
        gp: *mut libc::c_void, // fqz_gparams*, NULL for auto
    ) -> *mut libc::c_char;

    fn fqz_decompress(
        input: *const libc::c_char,
        in_size: libc::size_t,
        out_size: *mut libc::size_t,
        lengths: *mut libc::c_int,
        nlengths: libc::c_int,
    ) -> *mut libc::c_char;
}

/// Compress quality scores using fqzcomp's context-modeled arithmetic coder.
///
/// `qualities` is a slice of quality strings (one per read).
/// `strat` selects one of htscodecs' fixed fqzcomp parameter PRESETS
/// (`strat_opts[][12]` in `fqzcomp_qual.c`), NOT a compression-effort level.
/// `strat=0` is the general/auto-tuned preset and is OPTIMAL for Illumina short
/// reads; presets 1-3 are tuned for other data shapes and measured 1-4% WORSE on
/// both BAM and FASTQ quality with no speed upside (2026-06-23, `bench_fqz_strat`).
/// Keep production callers at 0.
/// Returns compressed bytes. Equivalent to `compress_with_flags` with all-zero
/// flags — kept as the stable entry point so existing callers are byte-identical.
pub fn compress(qualities: &[&[u8]], strat: i32) -> Result<Vec<u8>> {
    if qualities.is_empty() {
        return Ok(Vec::new());
    }
    let flags = vec![0u32; qualities.len()];
    compress_with_flags(qualities, &flags, strat)
}

/// Compress quality scores, passing per-read fqzcomp flags (e.g. `FQZ_FREAD2`
/// for mate-2 reads, `FQZ_FREVERSE` for reverse-strand reads). The flags let
/// fqzcomp's auto-tuner split the context model by read class when it pays.
///
/// `flags.len()` must equal `qualities.len()`.
///
/// SAFETY/aliasing: htscodecs **mutates** `s->flags` during encode (it stashes
/// auto-tuned selector bits in the high 16 bits, then clears them again before
/// returning). We therefore copy `flags` into an owned scratch `Vec` and hand C
/// a pointer to that — never to the caller's slice. The caller's slice is left
/// untouched (asserted by `test_compress_with_flags_does_not_mutate_caller`).
pub fn compress_with_flags(qualities: &[&[u8]], flags: &[u32], strat: i32) -> Result<Vec<u8>> {
    if qualities.is_empty() {
        return Ok(Vec::new());
    }

    if qualities.len() != flags.len() {
        anyhow::bail!(
            "fqz_compress: flags length {} != qualities length {}",
            flags.len(),
            qualities.len()
        );
    }

    // Build concatenated quality buffer and length array
    let total_bytes: usize = qualities.iter().map(|q| q.len()).sum();
    let mut concat = Vec::with_capacity(total_bytes);
    let mut lengths: Vec<u32> = Vec::with_capacity(qualities.len());

    for q in qualities {
        concat.extend_from_slice(q);
        lengths.push(q.len() as u32);
    }

    // Owned, mutable scratch — C writes into this, NOT the caller's `flags`.
    let mut flags_scratch: Vec<u32> = flags.to_vec();

    if qualities.len() > libc::c_int::MAX as usize {
        anyhow::bail!("fqz_compress: too many reads ({})", qualities.len());
    }
    let slice = FqzSlice {
        num_records: qualities.len() as libc::c_int,
        len: lengths.as_ptr(),
        flags: flags_scratch.as_mut_ptr(),
    };

    let mut out_size: libc::size_t = 0;

    // Version: FQZ_VERS (5) << 8 | strat
    let vers = (5 << 8) | (strat as libc::c_int);

    let compressed_ptr = unsafe {
        fqz_compress(
            vers,
            &slice,
            concat.as_ptr() as *const libc::c_char,
            concat.len(),
            &mut out_size,
            strat as libc::c_int,
            std::ptr::null_mut(), // auto-detect parameters
        )
    };

    // Keep the scratch alive until after the FFI call has fully returned.
    drop(flags_scratch);

    if compressed_ptr.is_null() {
        anyhow::bail!("fqz_compress returned NULL");
    }

    // Copy result into a Rust Vec, then free the C allocation
    let result = unsafe {
        let slice = std::slice::from_raw_parts(compressed_ptr as *const u8, out_size);
        let v = slice.to_vec();
        libc::free(compressed_ptr as *mut libc::c_void);
        v
    };

    Ok(result)
}

/// Read the uncompressed-size hint at the start of an fqzcomp blob.
///
/// Mirrors htscodecs' `var_get_u32` with `BIG_END` defined (big-endian
/// base-128, MSB first, high bit = continuation, at most 5 bytes for a u32).
/// Returns `None` if the header is truncated or overlong. Used only to bound
/// the allocation `fqz_decompress` will make — the C decoder re-reads it.
fn fqz_decompressed_size_hint(buf: &[u8]) -> Option<u32> {
    let mut value: u32 = 0;
    for (i, &byte) in buf.iter().take(5).enumerate() {
        value = (value << 7) | (byte & 0x7f) as u32;
        if byte & 0x80 == 0 {
            return Some(value);
        }
        // 5th byte still had the continuation bit set: overlong for a u32.
        if i == 4 {
            return None;
        }
    }
    None
}

/// Decompress quality scores produced by fqzcomp.
///
/// Returns concatenated quality bytes and per-read lengths.
pub fn decompress(compressed: &[u8], num_reads: usize) -> Result<(Vec<u8>, Vec<i32>)> {
    if compressed.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    // Zero-length lengths Vec gives a dangling-but-aligned ptr. Whether htscodecs
    // dereferences it depends on its internal mode, so short-circuit here.
    if num_reads == 0 {
        return Ok((Vec::new(), Vec::new()));
    }

    let mut out_size: libc::size_t = 0;
    if num_reads > libc::c_int::MAX as usize {
        anyhow::bail!("fqz_decompress: too many reads ({})", num_reads);
    }

    // DoS guard. The blob begins with htscodecs' big-endian base-128 varint
    // giving the uncompressed size, which `fqz_decompress` mallocs without any
    // bound (the only cap in htscodecs is behind a fuzzing-only flag). Pre-read
    // that length and reject anything implausibly large for the read count —
    // quality is one byte per base, and 64 KiB per read is far above any real
    // sequencing read length. This stops a tiny crafted blob from forcing a
    // multi-GB allocation (amplified across rayon workers).
    match fqz_decompressed_size_hint(compressed) {
        Some(claimed) => {
            let max_out = (num_reads as u64).saturating_mul(64 * 1024);
            if claimed as u64 > max_out {
                anyhow::bail!(
                    "fqz_decompress: claimed size {claimed} exceeds bound {max_out} for {num_reads} reads"
                );
            }
        }
        None => anyhow::bail!("fqz_decompress: malformed length header"),
    }

    // SAFETY: `lengths` is sized to exactly `num_reads` and `nlengths` is passed as the
    // same `num_reads`. htscodecs only writes `lengths[rec]` under `rec < nlengths`
    // (`fqzcomp_qual.c`, the `if (lengths && rec < nlengths)` guard), so it can never write
    // past this buffer even if the blob's internal record count disagrees. If that upstream
    // guard is ever removed on a pin bump, this allocation must grow to match.
    let mut lengths: Vec<libc::c_int> = vec![0; num_reads];

    let decompressed_ptr = unsafe {
        fqz_decompress(
            compressed.as_ptr() as *const libc::c_char,
            compressed.len(),
            &mut out_size,
            lengths.as_mut_ptr(),
            num_reads as libc::c_int,
        )
    };

    if decompressed_ptr.is_null() {
        anyhow::bail!("fqz_decompress returned NULL");
    }

    let result = unsafe {
        let slice = std::slice::from_raw_parts(decompressed_ptr as *const u8, out_size);
        let v = slice.to_vec();
        libc::free(decompressed_ptr as *mut libc::c_void);
        v
    };

    Ok((result, lengths))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// FQZ_FREAD2 from fqzcomp_qual.h — marks a read as mate-2.
    const FQZ_FREAD2: u32 = 128;

    #[test]
    fn test_compress_with_flags_does_not_mutate_caller() {
        // htscodecs writes selector bits into s->flags during encode and clears
        // them before returning; if we ever passed the caller's slice straight
        // through, this would observe scribbled high bits. The owned scratch
        // copy must keep the caller's input pristine.
        let qualities: Vec<&[u8]> = vec![
            b"IIIIIIIIIIIIIIII",
            b"FFFFFFFFFF!!!!!!",
            b"IIIIIIIIIIIIIIII",
            b"00000000########",
        ];
        let flags = vec![0u32, FQZ_FREAD2, 0u32, FQZ_FREAD2];
        let before = flags.clone();
        let _ = compress_with_flags(&qualities, &flags, 0).unwrap();
        assert_eq!(flags, before, "caller flags slice was mutated by FFI");
    }

    #[test]
    fn test_compress_with_flags_roundtrips() {
        // Setting FQZ_FREAD2 on alternating reads must still roundtrip losslessly
        // (the selector is stored inside the blob; decode needs no flag input).
        let qualities: Vec<&[u8]> = vec![
            b"IIIIIIIIIIIIIIIIIIIIIIIII",
            b"FFFFFFFFFFFFFFFFFFFF!!!!",
            b"ABCDEFGHIJKLMNOPQRSTUVWX",
            b"5555555555555555555555",
        ];
        let flags = vec![0u32, FQZ_FREAD2, 0u32, FQZ_FREAD2];
        let compressed = compress_with_flags(&qualities, &flags, 0).unwrap();
        let (decompressed, lengths) = decompress(&compressed, qualities.len()).unwrap();
        let mut offset = 0;
        for (i, q) in qualities.iter().enumerate() {
            assert_eq!(lengths[i] as usize, q.len());
            assert_eq!(&decompressed[offset..offset + q.len()], *q);
            offset += q.len();
        }
    }

    #[test]
    fn test_compress_with_flags_length_mismatch_errors() {
        let qualities: Vec<&[u8]> = vec![b"IIII", b"FFFF"];
        let flags = vec![0u32]; // too short
        assert!(compress_with_flags(&qualities, &flags, 0).is_err());
    }

    #[test]
    fn test_compress_equals_compress_with_zero_flags() {
        let qualities: Vec<&[u8]> = vec![b"IIIIIIII", b"FFFF!!!!", b"ABCDEFGH"];
        let a = compress(&qualities, 0).unwrap();
        let zero = vec![0u32; qualities.len()];
        let b = compress_with_flags(&qualities, &zero, 0).unwrap();
        assert_eq!(a, b, "compress must equal compress_with_flags(zeros)");
    }

    #[test]
    fn test_decompress_rejects_implausible_size() {
        // The leading big-endian base-128 varint claims ~10M uncompressed bytes
        // for only 2 reads (bound = 2 * 64 KiB = 128 KiB), so the DoS guard must
        // reject it before htscodecs allocates anything.
        // 10_000_000 encoded big-endian base-128: 0x84 0xE2 0xAD 0x00.
        let blob = vec![0x84, 0xE2, 0xAD, 0x00, 0, 0, 0, 0, 0, 0, 0, 0];
        let err = decompress(&blob, 2).unwrap_err();
        assert!(
            err.to_string().contains("exceeds"),
            "expected size-bound rejection, got: {err}"
        );
    }

    #[test]
    fn test_roundtrip_strat1() {
        let qualities: Vec<&[u8]> = vec![
            b"IIIIIIIIIIIIIIIIIIIIIIIII",
            b"FFFFFFFFFFFFFFFFFFFF!!!!",
            b"ABCDEFGHIJKLMNOPQRSTUVWX",
        ];
        let compressed = compress(&qualities, 1).unwrap();
        let (decompressed, lengths) = decompress(&compressed, qualities.len()).unwrap();

        let mut offset = 0;
        for (i, q) in qualities.iter().enumerate() {
            assert_eq!(lengths[i] as usize, q.len());
            assert_eq!(&decompressed[offset..offset + q.len()], *q);
            offset += q.len();
        }
    }

    #[test]
    fn test_roundtrip_strat0() {
        let qualities: Vec<&[u8]> = vec![
            b"BBBBAAAABBBBAAAA",
            b"DDDDDDDDDDDDDDDD",
            b"FFFFFFFFFFFFFFFF",
            b"HHHHHHHHHHHHHHHH",
            b"IIIIIIIIIIIIIIII",
        ];
        let compressed = compress(&qualities, 0).unwrap();
        let (decompressed, lengths) = decompress(&compressed, qualities.len()).unwrap();

        let mut offset = 0;
        for (i, q) in qualities.iter().enumerate() {
            assert_eq!(
                lengths[i] as usize,
                q.len(),
                "Length mismatch for read {}",
                i
            );
            assert_eq!(
                &decompressed[offset..offset + q.len()],
                *q,
                "Data mismatch for read {}",
                i
            );
            offset += q.len();
        }
    }
}
