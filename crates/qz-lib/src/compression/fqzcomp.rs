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
/// `strat` controls compression strategy (0=fast, 1-3=higher compression).
/// Returns compressed bytes.
pub fn compress(qualities: &[&[u8]], strat: i32) -> Result<Vec<u8>> {
    if qualities.is_empty() {
        return Ok(Vec::new());
    }

    // Build concatenated quality buffer and length array
    let total_bytes: usize = qualities.iter().map(|q| q.len()).sum();
    let mut concat = Vec::with_capacity(total_bytes);
    let mut lengths: Vec<u32> = Vec::with_capacity(qualities.len());
    let flags: Vec<u32> = vec![0; qualities.len()]; // No flags needed

    for q in qualities {
        concat.extend_from_slice(q);
        lengths.push(q.len() as u32);
    }

    let slice = FqzSlice {
        num_records: qualities.len() as libc::c_int,
        len: lengths.as_ptr(),
        flags: flags.as_ptr(),
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

/// Decompress quality scores produced by fqzcomp.
///
/// Returns concatenated quality bytes and per-read lengths.
pub fn decompress(compressed: &[u8], num_reads: usize) -> Result<(Vec<u8>, Vec<i32>)> {
    if compressed.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }

    let mut out_size: libc::size_t = 0;
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
            assert_eq!(lengths[i] as usize, q.len(), "Length mismatch for read {}", i);
            assert_eq!(&decompressed[offset..offset + q.len()], *q, "Data mismatch for read {}", i);
            offset += q.len();
        }
    }
}
