use crate::compression::codec_ids::{CODEC_BSC, CODEC_COLUMNAR, CODEC_FQZCOMP};
use anyhow::{Result, bail};

/// Validate a paired-end v5 unified directory: per-entry codec legality (the
/// archive_type=1 legality matrix), contiguous chunks, per-chunk structural
/// completeness (mandatory roles, exactly one R2-header variant, uniform
/// record_count), and that each mate lane totals `num_reads` (pairs). The
/// byte-span/overlap checks live in the decoder's footer read, where file offsets
/// are known.
///
/// The codec matrix is kept **exactly as tight as the decoder's contract** (the
/// validator is also what `verify --fast` trusts): the Headers role is always
/// columnar-decoded (`paired/mod.rs` / `decompress_impl.rs` ignore the codec byte and
/// call `decode_columnar_headers`), so only `CODEC_COLUMNAR` is legal there; Qual is
/// only `CODEC_BSC`/`CODEC_FQZCOMP` (`CODEC_QUALITY_CTX` = 4 was removed and is rejected
/// at decode); and because streaming decode picks ONE mate-wide quality codec from the
/// first Qual entry, every Qual entry of a given mate must carry that same codec. This
/// stops a forged/legacy archive from passing the validator yet misdecoding.
pub(crate) fn validate_paired_directory(
    dir: &crate::compression::chunk_directory::ChunkDirectory,
) -> Result<()> {
    use crate::compression::chunk_directory::{ChunkDirEntry, StreamRole, GLOBAL_SENTINEL};
    use std::collections::{BTreeMap, BTreeSet};

    // Per-entry legality (archive_type=1 matrix) + group by chunk. `mate_qual_codec`
    // also enforces that a mate's Qual entries all share one codec (streaming decode
    // picks the mate-wide codec from the first Qual entry, so a mixed lane misdecodes).
    let mut by_chunk: BTreeMap<u32, Vec<&ChunkDirEntry>> = BTreeMap::new();
    let mut mate_qual_codec: BTreeMap<u8, u8> = BTreeMap::new();
    for (i, e) in dir.entries.iter().enumerate() {
        // Whole-archive globals (chunk_index == GLOBAL_SENTINEL): the ONLY one paired
        // carries is ChunkDecodedSizes (mate 0, raw), feeding direct-write decode
        // pre-sizing. It lives outside the per-chunk per-mate structure validated below,
        // so accept-and-skip it; reject any other global shape.
        if e.chunk_index == GLOBAL_SENTINEL {
            if e.role != StreamRole::ChunkDecodedSizes || e.mate != 0 {
                bail!("paired entry {i}: unexpected global (role {:?}, mate {})", e.role, e.mate);
            }
            continue;
        }
        let codec_ok = match (e.mate, e.role) {
            // Headers is always columnar-decoded — the codec byte is ignored at decode,
            // so only CODEC_COLUMNAR is a legal tag here.
            (1, StreamRole::Headers) | (2, StreamRole::Headers) => e.codec == CODEC_COLUMNAR,
            (2, StreamRole::HeaderDelta) => e.codec == CODEC_BSC,
            (1, StreamRole::Sequence) | (2, StreamRole::Sequence) => e.codec == CODEC_BSC,
            (1, StreamRole::Qual) | (2, StreamRole::Qual) => {
                // CODEC_QUALITY_CTX (4) is removed and rejected at decode — not legal here.
                e.codec == CODEC_BSC || e.codec == CODEC_FQZCOMP
            }
            // mate 0, Nmask, RcFlags, HeaderDelta-on-mate-1, mate>2, etc. are all illegal.
            _ => false,
        };
        if !codec_ok {
            bail!(
                "paired entry {i}: illegal (mate {}, role {:?}, codec {})",
                e.mate,
                e.role,
                e.codec
            );
        }
        if e.role == StreamRole::Qual {
            match mate_qual_codec.get(&e.mate) {
                None => {
                    mate_qual_codec.insert(e.mate, e.codec);
                }
                Some(&c0) if c0 == e.codec => {}
                Some(&c0) => bail!(
                    "paired mate {}: Qual codec varies across chunks ({c0} vs {}) — streaming \
                     decode uses one mate-wide codec",
                    e.mate,
                    e.codec
                ),
            }
        }
        by_chunk.entry(e.chunk_index).or_default().push(e);
    }

    // Chunk indices must be contiguous 0..num_chunks.
    let keys: BTreeSet<u32> = by_chunk.keys().copied().collect();
    if keys.len() as u64 != dir.num_chunks as u64
        || keys.iter().enumerate().any(|(i, k)| *k as usize != i)
    {
        bail!("paired: chunk indices must be contiguous 0..num_chunks");
    }

    // Per-chunk structure + lane totals. Within a chunk, R1 and R2 share the same
    // pair count, so all entries carry one uniform record_count; summing it once
    // per chunk gives each lane's total.
    let mut total: u64 = 0;
    for (c, ents) in &by_chunk {
        let rc = ents[0].record_count;
        if ents.iter().any(|e| e.record_count != rc) {
            bail!("paired chunk {c}: record_count differs within chunk");
        }
        let mut seen: BTreeSet<(u8, u8)> = BTreeSet::new();
        for e in ents {
            if !seen.insert((e.mate, e.role as u8)) {
                bail!(
                    "paired chunk {c}: duplicate (mate {}, role {:?})",
                    e.mate,
                    e.role
                );
            }
        }
        // Mandatory mate-1 roles.
        for need in [StreamRole::Headers, StreamRole::Sequence, StreamRole::Qual] {
            if !seen.contains(&(1, need as u8)) {
                bail!("paired chunk {c}: missing mate-1 role {:?}", need);
            }
        }
        // Mandatory mate-2 roles (Sequence, Qual) + exactly one R2 header variant.
        for need in [StreamRole::Sequence, StreamRole::Qual] {
            if !seen.contains(&(2, need as u8)) {
                bail!("paired chunk {c}: missing mate-2 role {:?}", need);
            }
        }
        let r2_indep = seen.contains(&(2, StreamRole::Headers as u8));
        let r2_delta = seen.contains(&(2, StreamRole::HeaderDelta as u8));
        if r2_indep == r2_delta {
            let found = if r2_indep { "both" } else { "none" };
            bail!(
                "paired chunk {c}: R2 header: found {found}, expected exactly one (Headers or HeaderDelta)"
            );
        }
        total = total
            .checked_add(rc)
            .ok_or_else(|| anyhow::anyhow!("paired: record_count sum overflow"))?;
    }
    if total != dir.num_reads {
        bail!(
            "paired: Σ record_count {total} != num_reads {} (pairs)",
            dir.num_reads
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- v5 unified directory validator tests ---

    use crate::compression::chunk_directory::{ChunkDirEntry, ChunkDirectory, StreamRole};

    fn uent(role: StreamRole, mate: u8, codec: u8, rc: u64) -> ChunkDirEntry {
        ChunkDirEntry {
            chunk_index: 0,
            role,
            mate,
            codec,
            offset: 0,
            length: 1,
            record_count: rc,
        }
    }

    // Legal single-chunk paired directory: R1 H/S/Q (mate 1) + R2 H/S/Q (mate 2),
    // all record_count = 5, R2 header is the independent (Headers/Columnar) variant.
    fn legal_two_mate_chunk() -> ChunkDirectory {
        ChunkDirectory {
            num_reads: 5,
            num_chunks: 1,
            entries: vec![
                uent(StreamRole::Headers, 1, CODEC_COLUMNAR, 5),
                uent(StreamRole::Sequence, 1, CODEC_BSC, 5),
                uent(StreamRole::Qual, 1, CODEC_FQZCOMP, 5),
                uent(StreamRole::Headers, 2, CODEC_COLUMNAR, 5),
                uent(StreamRole::Sequence, 2, CODEC_BSC, 5),
                uent(StreamRole::Qual, 2, CODEC_FQZCOMP, 5),
            ],
        }
    }

    #[test]
    fn paired_dir_accepts_legal_chunk() {
        assert!(validate_paired_directory(&legal_two_mate_chunk()).is_ok());
    }

    #[test]
    fn paired_dir_accepts_r2_header_delta() {
        // The delta variant (HeaderDelta/Bsc on mate 2) is equally legal.
        let mut dir = legal_two_mate_chunk();
        dir.entries[3] = uent(StreamRole::HeaderDelta, 2, CODEC_BSC, 5);
        assert!(validate_paired_directory(&dir).is_ok());
    }

    #[test]
    fn paired_dir_rejects_record_count_mismatch() {
        let mut dir = legal_two_mate_chunk();
        dir.entries.last_mut().unwrap().record_count += 1; // R2 qual count diverges
        assert!(validate_paired_directory(&dir).is_err());
    }

    #[test]
    fn paired_dir_rejects_illegal_codec_for_role() {
        let mut dir = legal_two_mate_chunk();
        dir.entries[1].codec = 99; // Sequence must be CODEC_BSC
        assert!(validate_paired_directory(&dir).is_err());
    }

    #[test]
    fn paired_dir_rejects_headerdelta_on_r1() {
        let mut dir = legal_two_mate_chunk();
        dir.entries[0] = uent(StreamRole::HeaderDelta, 1, CODEC_BSC, 5); // mate 1 may not be HeaderDelta
        assert!(validate_paired_directory(&dir).is_err());
    }

    #[test]
    fn paired_dir_rejects_both_r2_header_roles_unified() {
        let mut dir = legal_two_mate_chunk();
        // Add a second R2 header variant alongside the existing indep one.
        dir.entries
            .push(uent(StreamRole::HeaderDelta, 2, CODEC_BSC, 5));
        assert!(validate_paired_directory(&dir).is_err());
    }

    #[test]
    fn paired_dir_rejects_mate_zero() {
        let mut dir = legal_two_mate_chunk();
        dir.entries[1].mate = 0; // single-end mate is illegal in a paired archive
        assert!(validate_paired_directory(&dir).is_err());
    }

    #[test]
    fn paired_dir_rejects_bsc_header_codec() {
        // Headers is always columnar-decoded (the codec byte is ignored at decode), so a
        // BSC-tagged Headers entry must be rejected — otherwise it passes the validator
        // and verify --fast, then misparses through decode_columnar_headers.
        let mut dir = legal_two_mate_chunk();
        dir.entries[0].codec = CODEC_BSC; // R1 Headers with BSC is illegal
        assert!(validate_paired_directory(&dir).is_err());
    }

    #[test]
    fn paired_dir_rejects_quality_ctx_qual() {
        // CODEC_QUALITY_CTX (4) was removed and is rejected at decode; the validator must
        // not accept it either (else it passes verify --fast but bails during decode).
        use crate::compression::codec_ids::CODEC_QUALITY_CTX;
        let mut dir = legal_two_mate_chunk();
        dir.entries[2].codec = CODEC_QUALITY_CTX; // R1 Qual
        assert!(validate_paired_directory(&dir).is_err());
    }

    #[test]
    fn paired_dir_rejects_mixed_qual_codec_across_chunks() {
        // Streaming decode picks one mate-wide Qual codec from the first entry, so a mate
        // whose Qual codec differs between chunks would misdecode → must be rejected.
        let mut dir = legal_two_mate_chunk(); // chunk 0: R1 Qual = FQZ, R2 Qual = FQZ
        dir.num_chunks = 2;
        dir.num_reads = 10;
        let mut c1 = vec![
            uent(StreamRole::Headers, 1, CODEC_COLUMNAR, 5),
            uent(StreamRole::Sequence, 1, CODEC_BSC, 5),
            uent(StreamRole::Qual, 1, CODEC_BSC, 5), // R1 Qual differs from chunk 0's FQZ
            uent(StreamRole::Headers, 2, CODEC_COLUMNAR, 5),
            uent(StreamRole::Sequence, 2, CODEC_BSC, 5),
            uent(StreamRole::Qual, 2, CODEC_FQZCOMP, 5),
        ];
        for e in &mut c1 {
            e.chunk_index = 1;
        }
        dir.entries.extend(c1);
        let err = validate_paired_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("Qual codec varies"), "got: {err}");
    }

    #[test]
    fn paired_dir_rejects_missing_mandatory_role() {
        // Drop R1's Qual entry → the mate-1 mandatory-role check must fire.
        let mut dir = legal_two_mate_chunk();
        dir.entries.remove(2); // entries[2] = (Qual, mate 1)
        let err = validate_paired_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("missing mate-1 role"), "got: {err}");
    }

    #[test]
    fn paired_dir_rejects_no_r2_header_role() {
        // Replace the R2 header (entries[3]) with a duplicate R2 Qual so neither
        // Headers nor HeaderDelta is present for mate 2.
        let mut dir = legal_two_mate_chunk();
        dir.entries.remove(3); // remove the R2 Headers entry → zero R2 header roles
        let err = validate_paired_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("found none"), "got: {err}");
    }

    #[test]
    fn paired_dir_rejects_chunk_count_mismatch() {
        // num_chunks claims 2 but only chunk 0 is present → contiguity check fires.
        let mut dir = legal_two_mate_chunk();
        dir.num_chunks = 2;
        let err = validate_paired_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("contiguous"), "got: {err}");
    }

    #[test]
    fn paired_dir_rejects_wrong_num_reads() {
        // Internally-consistent chunk (uniform rc) but num_reads disagrees with the
        // record_count sum → the Σ-vs-num_reads check fires (distinct from the
        // within-chunk record_count mismatch path).
        let mut dir = legal_two_mate_chunk();
        dir.num_reads = 99; // entries all have rc=5, sum=5
        let err = validate_paired_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("num_reads"), "got: {err}");
    }
}
