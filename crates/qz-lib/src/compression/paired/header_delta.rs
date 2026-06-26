use crate::compression::dna_utils::{read_varint, write_varint};
use anyhow::{Result, anyhow, bail};

/// Encode R2 ids as deltas against R1 ids (paired by index). Returns the RAW op
/// stream and per-record start offsets (one per record + a final sentinel == ops.len()).
/// Panics if the two slices differ in length (caller guarantees equal mate counts).
pub fn encode(r1_ids: &[&[u8]], r2_ids: &[&[u8]]) -> (Vec<u8>, Vec<usize>) {
    assert_eq!(r1_ids.len(), r2_ids.len(), "encode: mate count mismatch");
    let mut ops = Vec::new();
    let mut offsets = Vec::with_capacity(r1_ids.len() + 1);
    for (r1, r2) in r1_ids.iter().zip(r2_ids.iter()) {
        offsets.push(ops.len());
        if r1 == r2 {
            ops.push(0u8); // Same
        } else if r1.len() == r2.len() {
            let mut diff_pos = None;
            let mut diff_count = 0u32;
            for (i, (a, b)) in r1.iter().zip(r2.iter()).enumerate() {
                if a != b {
                    diff_count += 1;
                    diff_pos = Some(i);
                    if diff_count > 1 {
                        break;
                    }
                }
            }
            if diff_count == 1 {
                let off = diff_pos.unwrap();
                ops.push(1u8); // Substitute
                write_varint(&mut ops, off as u64);
                ops.push(r2[off]);
            } else {
                ops.push(2u8); // Literal
                write_varint(&mut ops, r2.len() as u64);
                ops.extend_from_slice(r2);
            }
        } else {
            ops.push(2u8); // Literal
            write_varint(&mut ops, r2.len() as u64);
            ops.extend_from_slice(r2);
        }
    }
    offsets.push(ops.len()); // sentinel
    (ops, offsets)
}

/// Adapter: read a LEB128 varint from `data` at position `*o`, advancing `*o`.
/// Returns Err (never panics) if the data is truncated.
fn read_varint_here(data: &[u8], o: &mut usize) -> Result<u64> {
    read_varint(data, o).ok_or_else(|| anyhow!("varint truncated at offset {}", *o))
}

/// Decode `n` R2 ids from the RAW op stream and the R1 ids. Validates offsets,
/// literal lengths, opcode count == n, and full consumption.
pub fn decode(ops: &[u8], r1_ids: &[Vec<u8>], n: usize) -> Result<Vec<Vec<u8>>> {
    if r1_ids.len() < n {
        bail!("decode: have {} r1 ids, need {n}", r1_ids.len());
    }
    // belt-and-suspenders: each record consumes >= 1 op byte, so the real count
    // can't exceed ops.len(); cap the capacity hint against an untrusted n.
    let mut out = Vec::with_capacity(n.min(ops.len() + 1));
    let mut o = 0usize;
    for i in 0..n {
        out.push(apply_one_op(ops, &mut o, &r1_ids[i]).map_err(|e| anyhow!("{e} at record {i}"))?);
    }
    if o != ops.len() {
        bail!("ops has {} trailing bytes", ops.len() - o);
    }
    Ok(out)
}

/// Apply ONE R2 header-delta op against `r1`, advancing the op cursor `*o`, and
/// return the reconstructed R2 id.
///
/// Opcodes: `0` = copy `r1`; `1` = substitute (varint offset, 1 byte); `2` =
/// literal (varint len, `len` bytes). Each op is self-delimiting and depends only
/// on `r1` plus the bytes it consumes — so the op stream can be applied **record by
/// record while streaming**, against the lockstep R1 id, without materialising a
/// whole chunk. This is the streaming entry point used by the unified paired decode
/// writer; `decode` is the batch wrapper (whole-chunk r1_ids in memory).
pub(crate) fn apply_one_op(ops: &[u8], o: &mut usize, r1: &[u8]) -> Result<Vec<u8>> {
    let opcode = *ops.get(*o).ok_or_else(|| anyhow!("ops truncated"))?;
    *o += 1;
    match opcode {
        0 => Ok(r1.to_vec()),
        1 => {
            let off = read_varint_here(ops, o)? as usize;
            let byte = *ops
                .get(*o)
                .ok_or_else(|| anyhow!("substitute byte truncated"))?;
            *o += 1;
            if off >= r1.len() {
                bail!(
                    "substitute offset {off} out of range (r1 len {})",
                    r1.len()
                );
            }
            let mut rec = r1.to_vec();
            rec[off] = byte;
            Ok(rec)
        }
        2 => {
            let len = read_varint_here(ops, o)? as usize;
            // `*o + len` can overflow usize (len is varint-derived, up to
            // u64::MAX as usize) → debug panic before `.get`. Use checked add.
            let end = o
                .checked_add(len)
                .ok_or_else(|| anyhow!("literal len {len} overflow"))?;
            let bytes = ops
                .get(*o..end)
                .ok_or_else(|| anyhow!("literal len {len} exceeds remaining"))?;
            *o = end;
            Ok(bytes.to_vec())
        }
        x => bail!("unknown opcode {x}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_same_sub_literal() {
        let r1: Vec<Vec<u8>> = vec![b"@m:1 1:N:0:A".to_vec(), b"@x".to_vec(), b"@same".to_vec()];
        let r2: Vec<Vec<u8>> = vec![
            b"@m:1 2:N:0:A".to_vec(),
            b"@totally-different".to_vec(),
            b"@same".to_vec(),
        ];
        let r1r: Vec<&[u8]> = r1.iter().map(|v| v.as_slice()).collect();
        let r2r: Vec<&[u8]> = r2.iter().map(|v| v.as_slice()).collect();
        let (ops, offs) = encode(&r1r, &r2r);
        assert_eq!(offs.len(), 4); // 3 records + sentinel
        assert_eq!(*offs.last().unwrap(), ops.len());
        assert_eq!(decode(&ops, &r1, 3).unwrap(), r2);
    }

    #[test]
    fn decode_rejects_bad_offset() {
        let r1 = vec![b"@abc".to_vec()];
        let mut ops = Vec::new();
        ops.push(1u8); // Substitute
        crate::compression::dna_utils::write_varint(&mut ops, 99u64);
        ops.push(b'Z');
        assert!(
            decode(&ops, &r1, 1)
                .unwrap_err()
                .to_string()
                .contains("offset")
        );
    }

    #[test]
    fn decode_rejects_overflow_literal_len() {
        // opcode 2 (Literal) with a varint len of u64::MAX must return Err, not
        // panic on `o + len` overflow (debug builds panic on integer overflow).
        let r1 = vec![b"@abc".to_vec()];
        let mut ops = Vec::new();
        ops.push(2u8); // Literal
        crate::compression::dna_utils::write_varint(&mut ops, u64::MAX);
        let err = decode(&ops, &r1, 1).unwrap_err();
        assert!(
            err.to_string().contains("overflow") || err.to_string().contains("exceeds"),
            "unexpected error: {err}"
        );
    }
}
