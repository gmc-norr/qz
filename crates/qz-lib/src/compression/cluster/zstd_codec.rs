use anyhow::{bail, Result};
use std::io::Write;

const PRELUDE_VERSION: u8 = 1;

/// Compress concatenated canonical seq bytes with zstd long-distance matching.
/// Returns the block PAYLOAD = 4-byte prelude + zstd frame (the v5 frame header is added
/// by the caller via write_block_frame_header).
///
/// `n_workers` sets zstd's internal multithreading (`ZSTD_c_nbWorkers`): with
/// long-distance matching on, zstd splits the input into parallel jobs while the LDM
/// table still spans them, so the whole-chunk long-range matching (the ratio lever) is
/// preserved. `0` = synchronous single-thread. The output is deterministic for a fixed
/// `n_workers` (job boundaries depend on input + params, not thread scheduling), so a
/// given (input, level, window_log, n_workers) always yields byte-identical bytes.
pub(crate) fn compress_seq_zstd_long(
    seq: &[u8],
    level: i32,
    window_log: u32,
    n_workers: u32,
) -> Result<Vec<u8>> {
    let mut out = vec![PRELUDE_VERSION, level as u8, window_log as u8, 0u8];
    let mut enc = zstd::stream::Encoder::new(Vec::new(), level)?;
    enc.long_distance_matching(true)?;
    enc.window_log(window_log)?;
    if n_workers > 0 {
        enc.multithread(n_workers)?;
    }
    enc.write_all(seq)?;
    let frame = enc.finish()?;
    out.extend_from_slice(&frame);
    Ok(out)
}

/// Inverse of compress_seq_zstd_long. Reads the prelude, sets window_log_max, and inflates with a
/// HARD output cap of `expected_len` (the sum of this chunk's record seq-lengths, known from the
/// quality decode / const-length). Caps the decompression bomb: a malicious frame cannot expand
/// past the bytes the directory says the chunk holds.
pub(crate) fn decompress_seq_zstd_long(payload: &[u8], expected_len: usize) -> Result<Vec<u8>> {
    use std::io::Read;
    if payload.len() < 4 || payload[0] != PRELUDE_VERSION {
        bail!("cluster seq payload: bad prelude");
    }
    let window_log = payload[2] as u32;
    let mut dec = zstd::stream::Decoder::new(&payload[4..])?;
    dec.window_log_max(window_log)?;
    let mut out = Vec::with_capacity(expected_len);
    // read at most expected_len + 1: if the frame yields more, it's corrupt/hostile.
    dec.take(expected_len as u64 + 1).read_to_end(&mut out)?;
    if out.len() != expected_len {
        bail!("cluster seq: decoded {} bytes, expected {}", out.len(), expected_len);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn zstd_long_roundtrips_with_prelude() {
        let seq: Vec<u8> = (0..200_000u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let payload = compress_seq_zstd_long(&seq, 16, 27, 0).unwrap();
        assert_eq!(payload[0], 1, "prelude version");
        assert_eq!(payload[1], 16, "level");
        assert_eq!(payload[2], 27, "window_log");
        // bounded decode: exact expected length succeeds
        let back = decompress_seq_zstd_long(&payload, seq.len()).unwrap();
        assert_eq!(back, seq);
        // a too-small expected length is rejected (decompression-bomb guard)
        assert!(decompress_seq_zstd_long(&payload, seq.len() - 1).is_err());
    }

    #[test]
    fn zstd_long_multithread_roundtrips_and_is_deterministic() {
        // larger input so MT actually splits into jobs; LDM spans them.
        let seq: Vec<u8> = (0..4_000_000u32).map(|i| b"ACGT"[((i * 7 + i / 3) % 4) as usize]).collect();
        let a = compress_seq_zstd_long(&seq, 16, 24, 4).unwrap();
        let b = compress_seq_zstd_long(&seq, 16, 24, 4).unwrap();
        assert_eq!(a, b, "MT output must be deterministic for fixed n_workers");
        let back = decompress_seq_zstd_long(&a, seq.len()).unwrap();
        assert_eq!(back, seq, "MT-compressed frame must round-trip");
    }
}
