//! Reference-mode role contract + directory validator (spec §9, "locked").
//!
//! Reference archives ride the shared v5 `ChunkDirectory` footer (`archive_type=2`)
//! but carry a richer role set (`RefRole`) and a much stricter validator: four
//! cardinality-1 **global** entries (consensus / interval map / N-bitmap /
//! reference_meta) plus per-chunk, per-mate sequence roles whose declared
//! `record_count`s have role-specific meanings (§9.2) and must satisfy the
//! cross-stream relations of §9.3(b), and an optional per-(chunk,mate) FallbackPool
//! role (present iff `N − M > 0`). Global entries use sentinels
//! (`chunk_index = GLOBAL_SENTINEL`, `mate = 0`).
//!
//! `RefRole`↔`StreamRole` mapping (`ref_role_to_stream`/`stream_role_to_ref`) lets
//! the reference encoder/validator speak its richer role vocabulary over the shared
//! directory's `StreamRole` byte. The directory's own bounded-allocation /
//! exact-consumption / checked-offset parsing lives in `chunk_directory`.

use crate::compression::chunk_directory::{ChunkDirEntry, ChunkDirectory, StreamRole};
use anyhow::{Result, bail};
use std::collections::{BTreeMap, BTreeSet};

/// `chunk_index` value tagging a global (non-per-chunk) entry.
pub const GLOBAL_SENTINEL: u32 = u32::MAX;
/// Per-chunk entry cap (spec §9.4): 7×2 seq roles + 2 header + 2 quality ≈ 18; 24 with margin.
pub const PER_CHUNK_ENTRY_CAP: u64 = 24;
/// Exactly four global roles, each once (spec §9.4): PackedBacking, IntervalMap,
/// NBitmap, ReferenceMeta. FallbackPool is a per-(chunk,mate) role, not global.
pub const GLOBAL_ENTRY_CAP: u64 = 4;

/// Codec codes used on the wire (per-entry `codec` byte). Shared with paired mode
/// via `compression::codec_ids` (single source of truth); re-exported here so the
/// existing `reference::format::CODEC_*` paths keep working.
pub(crate) use crate::compression::codec_ids::{
    CODEC_BSC, CODEC_COLUMNAR, CODEC_FQZCOMP, CODEC_PACKED_CONSENSUS,
};

#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, Debug, PartialOrd, Ord)]
pub enum RefRole {
    PackedBacking = 0,
    IntervalMap = 1,
    NBitmap = 2,
    ReferenceMeta = 3,
    MappedFlags = 4,
    Positions = 5,
    Strands = 6,
    ReadLen = 7,
    EditCount = 8,
    EditPos = 9,
    EditBase = 10,
    FallbackPool = 11,
    R1Headers = 12,
    R2HeadersIndep = 13,
    R2HeaderDelta = 14,
    R1Qual = 15,
    R2Qual = 16,
}

impl RefRole {
    // Test-only since the v5 switch: live code maps roles via `ref_role_to_stream`
    // and never decodes a raw `RefRole` byte. Kept for the role-table unit test.
    #[allow(dead_code)]
    pub fn from_u8(v: u8) -> Option<RefRole> {
        Some(match v {
            0 => RefRole::PackedBacking,
            1 => RefRole::IntervalMap,
            2 => RefRole::NBitmap,
            3 => RefRole::ReferenceMeta,
            4 => RefRole::MappedFlags,
            5 => RefRole::Positions,
            6 => RefRole::Strands,
            7 => RefRole::ReadLen,
            8 => RefRole::EditCount,
            9 => RefRole::EditPos,
            10 => RefRole::EditBase,
            11 => RefRole::FallbackPool,
            12 => RefRole::R1Headers,
            13 => RefRole::R2HeadersIndep,
            14 => RefRole::R2HeaderDelta,
            15 => RefRole::R1Qual,
            16 => RefRole::R2Qual,
            _ => return None,
        })
    }

    /// The four global roles (codes 0–3, cardinality 1, `chunk_index =
    /// GLOBAL_SENTINEL`): `PackedBacking`/`IntervalMap`/`NBitmap`/`ReferenceMeta`.
    /// `FallbackPool` is per-(chunk,mate), not global.
    pub fn is_global(self) -> bool {
        (self as u8) <= 3
    }

    /// The 7 mandatory per-chunk-per-mate sequence roles (spec §9.1/§9.2).
    /// `FallbackPool` is an optional per-(chunk,mate) role (omitted when N==M),
    /// not in `PER_CHUNK_SEQ`.
    pub const PER_CHUNK_SEQ: [RefRole; 7] = [
        RefRole::MappedFlags,
        RefRole::Positions,
        RefRole::Strands,
        RefRole::ReadLen,
        RefRole::EditCount,
        RefRole::EditPos,
        RefRole::EditBase,
    ];

    /// Per-entry codec legality (enforced in `validate_reference_directory`).
    pub fn codec_valid(self, codec: u8) -> bool {
        match self {
            RefRole::PackedBacking => codec == CODEC_PACKED_CONSENSUS,
            RefRole::FallbackPool => codec == CODEC_BSC,
            RefRole::IntervalMap | RefRole::NBitmap | RefRole::ReferenceMeta => codec == CODEC_BSC,
            // The 7 per-chunk sequence roles are all BSC.
            RefRole::MappedFlags
            | RefRole::Positions
            | RefRole::Strands
            | RefRole::ReadLen
            | RefRole::EditCount
            | RefRole::EditPos
            | RefRole::EditBase => codec == CODEC_BSC,
            RefRole::R1Headers | RefRole::R2HeadersIndep => codec == CODEC_COLUMNAR,
            RefRole::R2HeaderDelta => codec == CODEC_BSC,
            RefRole::R1Qual | RefRole::R2Qual => codec == CODEC_BSC || codec == CODEC_FQZCOMP,
        }
    }
}

/// Map a reference-mode `RefRole` onto the unified v5 `StreamRole`.
///
/// Header/quality roles fold by mate: `R1Headers`/`R2HeadersIndep` → `Headers`,
/// `R2HeaderDelta` → `HeaderDelta`, `R1Qual`/`R2Qual` → `Qual`. The 12
/// reference-specific roles map 1:1 to their same-named `StreamRole`. The `mate`
/// of the resulting unified entry is recovered separately (1 for R1*, 2 for R2*,
/// 0 for the four globals); together `(StreamRole, mate)` round-trips back to the
/// `RefRole` via `stream_role_to_ref`.
pub(crate) fn ref_role_to_stream(r: RefRole) -> StreamRole {
    match r {
        RefRole::R1Headers | RefRole::R2HeadersIndep => StreamRole::Headers,
        RefRole::R2HeaderDelta => StreamRole::HeaderDelta,
        RefRole::R1Qual | RefRole::R2Qual => StreamRole::Qual,
        RefRole::PackedBacking => StreamRole::PackedBacking,
        RefRole::IntervalMap => StreamRole::IntervalMap,
        RefRole::NBitmap => StreamRole::NBitmap,
        RefRole::ReferenceMeta => StreamRole::ReferenceMeta,
        RefRole::MappedFlags => StreamRole::MappedFlags,
        RefRole::Positions => StreamRole::Positions,
        RefRole::Strands => StreamRole::Strands,
        RefRole::ReadLen => StreamRole::ReadLen,
        RefRole::EditCount => StreamRole::EditCount,
        RefRole::EditPos => StreamRole::EditPos,
        RefRole::EditBase => StreamRole::EditBase,
        RefRole::FallbackPool => StreamRole::FallbackPool,
    }
}

/// Inverse of `ref_role_to_stream`: recover the reference `RefRole` from a unified
/// entry's `(StreamRole, mate)`. `mate` disambiguates the only many-to-one cases:
/// `Headers`/`Qual` (mate 1 → R1*, mate 2 → R2*) and `HeaderDelta` (mate 2 only).
/// All other roles map 1:1 regardless of mate (their mate is carried for the
/// per-`(chunk,mate)` grouping but is not needed to recover the role). Returns
/// `Err` for `(role, mate)` combinations that have no reference role — e.g. the
/// single-end-only `Sequence`/`Nmask`/`RcFlags`, a mate-1 `HeaderDelta`, or a
/// `Headers`/`Qual` with mate ∉ {1,2}. This keeps the codec matrix single-sourced
/// (`RefRole::codec_valid`) and rejects roles outside the reference contract.
fn stream_role_to_ref(role: StreamRole, mate: u8) -> Result<RefRole> {
    Ok(match (role, mate) {
        (StreamRole::Headers, 1) => RefRole::R1Headers,
        (StreamRole::Headers, 2) => RefRole::R2HeadersIndep,
        (StreamRole::HeaderDelta, 2) => RefRole::R2HeaderDelta,
        (StreamRole::Qual, 1) => RefRole::R1Qual,
        (StreamRole::Qual, 2) => RefRole::R2Qual,
        (StreamRole::PackedBacking, _) => RefRole::PackedBacking,
        (StreamRole::IntervalMap, _) => RefRole::IntervalMap,
        (StreamRole::NBitmap, _) => RefRole::NBitmap,
        (StreamRole::ReferenceMeta, _) => RefRole::ReferenceMeta,
        (StreamRole::MappedFlags, _) => RefRole::MappedFlags,
        (StreamRole::Positions, _) => RefRole::Positions,
        (StreamRole::Strands, _) => RefRole::Strands,
        (StreamRole::ReadLen, _) => RefRole::ReadLen,
        (StreamRole::EditCount, _) => RefRole::EditCount,
        (StreamRole::EditPos, _) => RefRole::EditPos,
        (StreamRole::EditBase, _) => RefRole::EditBase,
        (StreamRole::FallbackPool, _) => RefRole::FallbackPool,
        _ => bail!("no reference role for (role {:?}, mate {})", role, mate),
    })
}

/// Reference-mode validator for the **unified v5 directory** (`archive_type=2`).
/// Carries the SEMANTIC reference footer checks on `ChunkDirectory`
/// entries: per-entry codec legality (via `RefRole::codec_valid`, recovering the
/// `RefRole` from `(StreamRole, mate)`), the entry-count cap, exactly-four globals
/// (each once, mate 0, `chunk_index == GLOBAL_SENTINEL`), `NBitmap == PackedBacking`
/// record_count, contiguous chunks, per-`(chunk,mate)` mandatory-role completeness,
/// the §9.3(b) cross-stream count relations, per-chunk R1==R2 equality, and the
/// per-`(chunk,mate)` FallbackPool ownership contract (present iff `N_chunk − M > 0`,
/// `record_count == N_chunk − M`, never present-but-empty).
///
/// The byte-level span/overlap/offset checks are intentionally NOT done here —
/// they are owned by the shared `chunk_directory::validate_and_parse_footer`,
/// which runs on the raw footer bytes before this function sees the parsed
/// directory.
pub(crate) fn validate_reference_directory(dir: &ChunkDirectory) -> Result<()> {
    // Entry-count cap (§9.4): num_chunks·PER_CHUNK_ENTRY_CAP + GLOBAL_ENTRY_CAP.
    //
    // Enforced UP FRONT here. The directory is already parsed by the
    // shared `ChunkDirectory::parse` (whose own remaining-bytes guard bounds the
    // allocation), so this cap is a semantic backstop. A valid reference directory
    // has 5 + 19·num_chunks entries (always ≤ the 5 + 24·num_chunks bound), so the
    // only way to trip this is a hostile footer whose declared `num_chunks`
    // understates its entry list — exactly the case we reject early.
    // +1 for the optional cross-cutting `ChunkDecodedSizes` global (one block, GLOBAL_SENTINEL,
    // emitted by the encoder so NUMA decode can direct-write; not a reference RefRole).
    let entry_bound = (dir.num_chunks as u64)
        .checked_mul(PER_CHUNK_ENTRY_CAP)
        .and_then(|m| m.checked_add(GLOBAL_ENTRY_CAP))
        .and_then(|m| m.checked_add(1))
        .ok_or_else(|| anyhow::anyhow!("entry bound overflow"))?;
    if dir.entries.len() as u64 > entry_bound {
        bail!(
            "entry count {} exceeds cap {entry_bound}",
            dir.entries.len()
        );
    }

    // Per-entry codec legality + partition global vs per-chunk. The codec matrix
    // stays single-sourced in `RefRole::codec_valid`; we recover the RefRole from
    // `(role, mate)` (which also rejects roles outside the reference contract,
    // e.g. single-end-only Sequence/Nmask/RcFlags).
    let mut globals: BTreeMap<RefRole, &ChunkDirEntry> = BTreeMap::new();
    // (chunk_index, mate) -> map of RefRole -> entry
    let mut by_group: BTreeMap<(u32, u8), BTreeMap<RefRole, &ChunkDirEntry>> = BTreeMap::new();
    for e in &dir.entries {
        // `ChunkDecodedSizes` is a cross-cutting global (the same StreamRole single/paired
        // use), NOT a reference RefRole — accept it (mate 0, GLOBAL_SENTINEL) and skip the
        // RefRole mapping + the "exactly 4 reference globals" set so that contract still
        // counts only PackedBacking/IntervalMap/NBitmap/ReferenceMeta.
        if e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::ChunkDecodedSizes {
            if e.mate != 0 {
                bail!("ChunkDecodedSizes global must have mate==0 (got {})", e.mate);
            }
            continue;
        }
        let is_global = e.chunk_index == GLOBAL_SENTINEL;
        let ref_role = stream_role_to_ref(e.role, e.mate)?;
        if !ref_role.codec_valid(e.codec) {
            bail!("invalid codec {} for {:?}", e.codec, ref_role);
        }
        if is_global {
            if !ref_role.is_global() {
                bail!("per-chunk role {:?} carries GLOBAL_SENTINEL", ref_role);
            }
            if e.mate != 0 {
                bail!("global {:?} must have mate==0", ref_role);
            }
            if globals.insert(ref_role, e).is_some() {
                bail!("duplicate global role {:?}", ref_role);
            }
        } else {
            if ref_role.is_global() {
                bail!("global role {:?} with non-sentinel chunk_index", ref_role);
            }
            if e.mate != 1 && e.mate != 2 {
                bail!("per-chunk {:?} has bad mate {}", ref_role, e.mate);
            }
            let g = by_group.entry((e.chunk_index, e.mate)).or_default();
            if g.insert(ref_role, e).is_some() {
                bail!(
                    "duplicate (chunk={}, mate={}, role={:?})",
                    e.chunk_index,
                    e.mate,
                    ref_role
                );
            }
        }
    }

    // Exactly the 4 globals, each once.
    for need in [
        RefRole::PackedBacking,
        RefRole::IntervalMap,
        RefRole::NBitmap,
        RefRole::ReferenceMeta,
    ] {
        if !globals.contains_key(&need) {
            bail!("missing global role {:?}", need);
        }
    }
    if globals.len() as u64 != GLOBAL_ENTRY_CAP {
        bail!("expected {GLOBAL_ENTRY_CAP} globals, got {}", globals.len());
    }

    // §9.3(d) declared-count relation: NBitmap bases == backing bases.
    let backing_bases = globals[&RefRole::PackedBacking].record_count;
    if globals[&RefRole::NBitmap].record_count != backing_bases {
        bail!(
            "N-bitmap record_count {} != backing record_count {}",
            globals[&RefRole::NBitmap].record_count,
            backing_bases
        );
    }

    // Chunk indices contiguous 0..num_chunks (per-mate groups collapse to the
    // chunk index set).
    let chunk_indices: BTreeSet<u32> = by_group.keys().map(|(c, _)| *c).collect();
    if chunk_indices.len() as u64 != dir.num_chunks as u64
        || chunk_indices
            .iter()
            .enumerate()
            .any(|(i, k)| *k as usize != i)
    {
        bail!("chunk indices must be contiguous 0..num_chunks");
    }

    // Per-(chunk,mate) structural completeness + per-group cross-stream counts
    // (§9.3a/b). Then per-chunk R1/R2 N_chunk equality + Σ N_chunk. Each group's
    // F_group = checked_sub(N, M) is cross-checked against THAT group's own
    // FallbackPool entry (§4.3.3): present iff F_group > 0, record_count == F_group.
    let mut total_pairs: u64 = 0;
    for chunk in 0..dir.num_chunks {
        let mut n_per_mate: [Option<u64>; 3] = [None, None, None]; // index by mate (1,2)
        for mate in [1u8, 2u8] {
            let group = by_group
                .get(&(chunk, mate))
                .ok_or_else(|| anyhow::anyhow!("chunk {chunk} mate {mate} missing entirely"))?;
            // 7 mandatory per-chunk sequence roles present.
            for need in RefRole::PER_CHUNK_SEQ {
                if !group.contains_key(&need) {
                    bail!("chunk {chunk} mate {mate} missing role {:?}", need);
                }
            }
            // Header role: mate 1 -> R1Headers; mate 2 -> exactly one of
            // {R2HeadersIndep, R2HeaderDelta}. Quality: R1Qual / R2Qual.
            if mate == 1 {
                if !group.contains_key(&RefRole::R1Headers) {
                    bail!("chunk {chunk} mate 1 missing R1Headers");
                }
                if !group.contains_key(&RefRole::R1Qual) {
                    bail!("chunk {chunk} mate 1 missing R1Qual");
                }
                // mate-1 group must not carry mate-2 roles.
                if group.contains_key(&RefRole::R2HeadersIndep)
                    || group.contains_key(&RefRole::R2HeaderDelta)
                    || group.contains_key(&RefRole::R2Qual)
                {
                    bail!("chunk {chunk} mate 1 carries mate-2 role");
                }
            } else {
                let indep = group.contains_key(&RefRole::R2HeadersIndep);
                let delta = group.contains_key(&RefRole::R2HeaderDelta);
                if indep == delta {
                    bail!("chunk {chunk} mate 2: exactly one R2 header role required");
                }
                if !group.contains_key(&RefRole::R2Qual) {
                    bail!("chunk {chunk} mate 2 missing R2Qual");
                }
                if group.contains_key(&RefRole::R1Headers) || group.contains_key(&RefRole::R1Qual) {
                    bail!("chunk {chunk} mate 2 carries mate-1 role");
                }
            }

            // §9.3(b) cross-stream counts within this group.
            let rc = |r: RefRole| group[&r].record_count;
            let n_chunk = rc(RefRole::MappedFlags);
            let m = rc(RefRole::Positions);
            if rc(RefRole::Strands) != m || rc(RefRole::ReadLen) != m || rc(RefRole::EditCount) != m
            {
                bail!(
                    "chunk {chunk} mate {mate}: positions/strands/read_len/edit_count must all == M"
                );
            }
            // M ≤ N_chunk (checked_sub guards against hostile underflow);
            // F_group = N − M is this group's fallback-literal count.
            let f_group = n_chunk.checked_sub(m).ok_or_else(|| {
                anyhow::anyhow!("chunk {chunk} mate {mate}: M ({m}) > N_chunk ({n_chunk})")
            })?;
            // Per-(chunk,mate) FallbackPool ownership: present iff F_group > 0,
            // record_count == F_group, never present-but-empty. All footer
            // record_counts are attacker-controllable, so reject every mismatch.
            match group.get(&RefRole::FallbackPool) {
                Some(fb) => {
                    if fb.record_count == 0 {
                        bail!(
                            "chunk {chunk} mate {mate}: FallbackPool present but empty (must be omitted when N==M)"
                        );
                    }
                    if fb.record_count != f_group {
                        bail!(
                            "chunk {chunk} mate {mate}: FallbackPool.record_count {} != F_group {f_group}",
                            fb.record_count
                        );
                    }
                }
                None => {
                    if f_group != 0 {
                        bail!(
                            "chunk {chunk} mate {mate}: FallbackPool absent but F_group {f_group} > 0"
                        );
                    }
                }
            }
            if rc(RefRole::EditPos) != rc(RefRole::EditBase) {
                bail!("chunk {chunk} mate {mate}: len(edit_pos) != len(edit_base)");
            }
            // Headers and qualities are per-read streams: exactly N_chunk records
            // each (the encoder pushes both with record_count = n_chunk, and decode
            // pulls N_chunk header/qual records per chunk to align with the
            // reconstructed sequences). The structural checks above only confirm
            // these roles are PRESENT — a forged record_count would silently
            // misalign headers/qualities against sequences without it. All footer
            // counts are attacker-controllable, so reject any skew. (Single-end and
            // paired get this for free from their uniform per-chunk record_count;
            // reference's per-role counts differ legitimately for M-streams, so the
            // per-read roles must be checked explicitly.)
            let (header_role, qual_role) = if mate == 1 {
                (RefRole::R1Headers, RefRole::R1Qual)
            } else if group.contains_key(&RefRole::R2HeadersIndep) {
                (RefRole::R2HeadersIndep, RefRole::R2Qual)
            } else {
                (RefRole::R2HeaderDelta, RefRole::R2Qual)
            };
            if rc(header_role) != n_chunk {
                bail!(
                    "chunk {chunk} mate {mate}: header record_count {} != N_chunk {n_chunk}",
                    rc(header_role)
                );
            }
            if rc(qual_role) != n_chunk {
                bail!(
                    "chunk {chunk} mate {mate}: quality record_count {} != N_chunk {n_chunk}",
                    rc(qual_role)
                );
            }
            n_per_mate[mate as usize] = Some(n_chunk);
        }
        let n1 = n_per_mate[1].unwrap();
        let n2 = n_per_mate[2].unwrap();
        if n1 != n2 {
            bail!("chunk {chunk}: R1 N_chunk ({n1}) != R2 N_chunk ({n2})");
        }
        total_pairs = total_pairs
            .checked_add(n1)
            .ok_or_else(|| anyhow::anyhow!("num_pairs overflow"))?;
    }
    if total_pairs != dir.num_reads {
        bail!("Σ N_chunk {total_pairs} != num_reads {}", dir.num_reads);
    }

    // (The §4.3.3 FallbackPool ownership contract is now enforced per-(chunk,mate)
    // inside the loop above; the byte-level span/overlap/offset checks remain owned
    // by `chunk_directory::validate_and_parse_footer`.)

    Ok(())
}

/// Single-end reference-mode validator for the **unified v5 directory**
/// (`archive_type=4`). Mirrors [`validate_reference_directory`] restricted to a
/// single mate (mate 1): the same four cardinality-1 globals, but per-chunk
/// **mate-1-only** roles — the 7 mandatory sequence roles + `R1Headers`
/// (columnar) + `R1Qual`, plus an optional mate-1 `FallbackPool` (present iff
/// `N_chunk − M > 0`). Any mate-2 / mate-0 per-chunk entry, any `HeaderDelta`,
/// or any `R2*` role is rejected. `Σ N_chunk == num_reads` (reads, not pairs).
/// Run at both encode self-check and decode.
pub(crate) fn validate_reference_directory_single(dir: &ChunkDirectory) -> Result<()> {
    // Per-chunk allocation backstop. A legitimate single-end chunk has ≤10 mate-1
    // entries: 7 sequence columns (MappedFlags/Positions/Strands/ReadLen/EditCount/
    // EditPos/EditBase) + R1Headers + R1Qual + optional FallbackPool. Cap=12 gives a
    // small margin above that maximum while staying tight enough to be a meaningful
    // DoS bound. This is a coarse count guard, NOT the structural contract — the
    // per-entry mate-1 check below is what rejects any injected mate-2/mate-0 entry.
    const PER_CHUNK_ENTRY_CAP_SINGLE: u64 = 12;
    // +1 for the optional cross-cutting `ChunkDecodedSizes` global (see paired validator).
    let entry_bound = (dir.num_chunks as u64)
        .checked_mul(PER_CHUNK_ENTRY_CAP_SINGLE)
        .and_then(|m| m.checked_add(GLOBAL_ENTRY_CAP))
        .and_then(|m| m.checked_add(1))
        .ok_or_else(|| anyhow::anyhow!("entry bound overflow"))?;
    if dir.entries.len() as u64 > entry_bound {
        bail!(
            "entry count {} exceeds cap {entry_bound}",
            dir.entries.len()
        );
    }

    // Per-entry codec legality + partition global vs per-chunk. Per-chunk roles
    // must all be mate 1 (single-end has no mate 2). Recovering the RefRole from
    // (role, mate) also rejects roles outside the reference contract.
    let mut globals: BTreeMap<RefRole, &ChunkDirEntry> = BTreeMap::new();
    let mut by_chunk: BTreeMap<u32, BTreeMap<RefRole, &ChunkDirEntry>> = BTreeMap::new();
    for e in &dir.entries {
        // `ChunkDecodedSizes` is a cross-cutting global (the same StreamRole single/paired
        // use), NOT a reference RefRole — accept it (mate 0, GLOBAL_SENTINEL) and skip the
        // RefRole mapping + the "exactly 4 reference globals" set so that contract still
        // counts only PackedBacking/IntervalMap/NBitmap/ReferenceMeta.
        if e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::ChunkDecodedSizes {
            if e.mate != 0 {
                bail!("ChunkDecodedSizes global must have mate==0 (got {})", e.mate);
            }
            continue;
        }
        let is_global = e.chunk_index == GLOBAL_SENTINEL;
        let ref_role = stream_role_to_ref(e.role, e.mate)?;
        if !ref_role.codec_valid(e.codec) {
            bail!("invalid codec {} for {:?}", e.codec, ref_role);
        }
        if is_global {
            if !ref_role.is_global() {
                bail!("per-chunk role {:?} carries GLOBAL_SENTINEL", ref_role);
            }
            if e.mate != 0 {
                bail!("global {:?} must have mate==0", ref_role);
            }
            if globals.insert(ref_role, e).is_some() {
                bail!("duplicate global role {:?}", ref_role);
            }
        } else {
            if ref_role.is_global() {
                bail!("global role {:?} with non-sentinel chunk_index", ref_role);
            }
            if e.mate != 1 {
                bail!(
                    "single-end reference per-chunk {:?} has mate {} (must be 1)",
                    ref_role,
                    e.mate
                );
            }
            // mate 2 / delta roles cannot occur in a single-end archive.
            if matches!(
                ref_role,
                RefRole::R2HeadersIndep | RefRole::R2HeaderDelta | RefRole::R2Qual
            ) {
                bail!("single-end reference carries mate-2 role {:?}", ref_role);
            }
            let g = by_chunk.entry(e.chunk_index).or_default();
            if g.insert(ref_role, e).is_some() {
                bail!(
                    "duplicate (chunk={}, mate=1, role={:?})",
                    e.chunk_index,
                    ref_role
                );
            }
        }
    }

    // Exactly the 4 globals, each once.
    for need in [
        RefRole::PackedBacking,
        RefRole::IntervalMap,
        RefRole::NBitmap,
        RefRole::ReferenceMeta,
    ] {
        if !globals.contains_key(&need) {
            bail!("missing global role {:?}", need);
        }
    }
    if globals.len() as u64 != GLOBAL_ENTRY_CAP {
        bail!("expected {GLOBAL_ENTRY_CAP} globals, got {}", globals.len());
    }
    let backing_bases = globals[&RefRole::PackedBacking].record_count;
    if globals[&RefRole::NBitmap].record_count != backing_bases {
        bail!(
            "N-bitmap record_count {} != backing record_count {}",
            globals[&RefRole::NBitmap].record_count,
            backing_bases
        );
    }

    // Chunk indices contiguous 0..num_chunks.
    let chunk_indices: BTreeSet<u32> = by_chunk.keys().copied().collect();
    if chunk_indices.len() as u64 != dir.num_chunks as u64
        || chunk_indices.iter().enumerate().any(|(i, k)| *k as usize != i)
    {
        bail!("chunk indices must be contiguous 0..num_chunks");
    }

    let mut total_reads: u64 = 0;
    for chunk in 0..dir.num_chunks {
        let group = by_chunk
            .get(&chunk)
            .ok_or_else(|| anyhow::anyhow!("chunk {chunk} missing entirely"))?;
        for need in RefRole::PER_CHUNK_SEQ {
            if !group.contains_key(&need) {
                bail!("chunk {chunk}: missing role {:?}", need);
            }
        }
        if !group.contains_key(&RefRole::R1Headers) {
            bail!("chunk {chunk}: missing R1Headers");
        }
        if !group.contains_key(&RefRole::R1Qual) {
            bail!("chunk {chunk}: missing R1Qual");
        }

        let rc = |r: RefRole| group[&r].record_count;
        let n_chunk = rc(RefRole::MappedFlags);
        let m = rc(RefRole::Positions);
        if rc(RefRole::Strands) != m || rc(RefRole::ReadLen) != m || rc(RefRole::EditCount) != m {
            bail!("chunk {chunk}: positions/strands/read_len/edit_count must all == M");
        }
        let f_group = n_chunk.checked_sub(m).ok_or_else(|| {
            anyhow::anyhow!("chunk {chunk}: M ({m}) > N_chunk ({n_chunk})")
        })?;
        match group.get(&RefRole::FallbackPool) {
            Some(fb) => {
                if fb.record_count == 0 {
                    bail!(
                        "chunk {chunk}: FallbackPool present but empty (must be omitted when N==M)"
                    );
                }
                if fb.record_count != f_group {
                    bail!(
                        "chunk {chunk}: FallbackPool.record_count {} != F_group {f_group}",
                        fb.record_count
                    );
                }
            }
            None => {
                if f_group != 0 {
                    bail!("chunk {chunk}: FallbackPool absent but F_group {f_group} > 0");
                }
            }
        }
        if rc(RefRole::EditPos) != rc(RefRole::EditBase) {
            bail!("chunk {chunk}: len(edit_pos) != len(edit_base)");
        }
        if rc(RefRole::R1Headers) != n_chunk {
            bail!(
                "chunk {chunk}: header record_count {} != N_chunk {n_chunk}",
                rc(RefRole::R1Headers)
            );
        }
        if rc(RefRole::R1Qual) != n_chunk {
            bail!(
                "chunk {chunk}: quality record_count {} != N_chunk {n_chunk}",
                rc(RefRole::R1Qual)
            );
        }
        total_reads = total_reads
            .checked_add(n_chunk)
            .ok_or_else(|| anyhow::anyhow!("num_reads overflow"))?;
    }
    if total_reads != dir.num_reads {
        bail!("Σ N_chunk {total_reads} != num_reads {}", dir.num_reads);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn role_table_round4_shape() {
        assert_eq!(RefRole::from_u8(0).unwrap(), RefRole::PackedBacking);
        assert_eq!(RefRole::from_u8(11).unwrap(), RefRole::FallbackPool);
        assert_eq!(GLOBAL_ENTRY_CAP, 4);
        assert_eq!(RefRole::PER_CHUNK_SEQ.len(), 7);
        assert!(!RefRole::PER_CHUNK_SEQ.contains(&RefRole::FallbackPool));
    }

    #[test]
    fn fqzcomp_codec_allowed_for_quality_roles() {
        // fqzcomp (6) and BSC (1) are valid for quality roles…
        assert!(RefRole::R1Qual.codec_valid(CODEC_FQZCOMP));
        assert!(RefRole::R2Qual.codec_valid(CODEC_FQZCOMP));
        assert!(RefRole::R1Qual.codec_valid(CODEC_BSC));
        // …but the removed quality_ctx codec (4) is NOT…
        assert!(!RefRole::R1Qual.codec_valid(crate::compression::codec_ids::CODEC_QUALITY_CTX));
        // …and fqzcomp is NOT valid for a non-quality role.
        assert!(!RefRole::EditPos.codec_valid(CODEC_FQZCOMP));
    }

    // --- v5 unified-directory reference validator (`validate_reference_directory`) ---
    // ChunkDirEntry/ChunkDirectory/StreamRole come in via `use super::*`.

    /// Unified-directory analog of the test `pc()` helper: a per-chunk entry with
    /// a fixed 10-byte payload at `offset`.
    fn ude(
        chunk: u32,
        mate: u8,
        role: StreamRole,
        codec: u8,
        offset: u64,
        record_count: u64,
    ) -> ChunkDirEntry {
        ChunkDirEntry {
            chunk_index: chunk,
            role,
            mate,
            codec,
            offset,
            length: 10,
            record_count,
        }
    }

    /// Unified-directory analog of `global()`: a whole-archive entry (mate 0,
    /// `chunk_index = GLOBAL_SENTINEL`).
    fn ude_global(role: StreamRole, codec: u8, offset: u64, record_count: u64) -> ChunkDirEntry {
        ChunkDirEntry {
            chunk_index: GLOBAL_SENTINEL,
            role,
            mate: 0,
            codec,
            offset,
            length: 10,
            record_count,
        }
    }

    /// The four globals: PackedBacking=100, IntervalMap=1, NBitmap=100,
    /// ReferenceMeta=0. FallbackPool is no longer global (per-(chunk,mate) now).
    fn ude_globals() -> Vec<ChunkDirEntry> {
        vec![
            ude_global(StreamRole::PackedBacking, CODEC_PACKED_CONSENSUS, 1000, 100),
            ude_global(StreamRole::IntervalMap, CODEC_BSC, 1010, 1),
            ude_global(StreamRole::NBitmap, CODEC_BSC, 1020, 100),
            ude_global(StreamRole::ReferenceMeta, CODEC_BSC, 1030, 0),
        ]
    }

    /// A per-(chunk,mate) FallbackPool entry with `record_count = rc` (= N−M for
    /// that group). Codec must be BSC. The payload offset/length are placeholders
    /// (the semantic validator does not inspect them).
    fn ude_chunk_fallback(chunk: u32, mate: u8, rc: u64) -> ChunkDirEntry {
        ChunkDirEntry {
            chunk_index: chunk,
            role: StreamRole::FallbackPool,
            mate,
            codec: CODEC_BSC,
            offset: 0,
            length: 1,
            record_count: rc,
        }
    }

    /// All mandatory per-chunk roles for one mate (translated from
    /// `chunk_entries`): mapped_flags=5 (N_chunk), positions=strands=read_len=
    /// edit_count=4 (M), edit_pos=edit_base=2 (E). Headers/HeaderDelta/Qual map by
    /// mate. Payload spans `[base + k*10, +10)` (non-overlapping).
    fn ude_chunk_entries(chunk: u32, mate: u8, base: u64) -> Vec<ChunkDirEntry> {
        let mut k = 0u64;
        let mut next = || {
            let off = base + k * 10;
            k += 1;
            off
        };
        let mut v = vec![
            ude(chunk, mate, StreamRole::MappedFlags, CODEC_BSC, next(), 5),
            ude(chunk, mate, StreamRole::Positions, CODEC_BSC, next(), 4),
            ude(chunk, mate, StreamRole::Strands, CODEC_BSC, next(), 4),
            ude(chunk, mate, StreamRole::ReadLen, CODEC_BSC, next(), 4),
            ude(chunk, mate, StreamRole::EditCount, CODEC_BSC, next(), 4),
            ude(chunk, mate, StreamRole::EditPos, CODEC_BSC, next(), 2),
            ude(chunk, mate, StreamRole::EditBase, CODEC_BSC, next(), 2),
        ];
        if mate == 1 {
            v.push(ude(
                chunk,
                1,
                StreamRole::Headers,
                CODEC_COLUMNAR,
                next(),
                5,
            ));
            v.push(ude(chunk, 1, StreamRole::Qual, CODEC_BSC, next(), 5));
        } else {
            v.push(ude(chunk, 2, StreamRole::HeaderDelta, CODEC_BSC, next(), 5));
            v.push(ude(chunk, 2, StreamRole::Qual, CODEC_BSC, next(), 5));
        }
        v
    }

    /// A minimal valid reference directory: 4 globals + one chunk (R1+R2).
    /// num_reads = Σ N_chunk = 5 pairs. Each mate's FallbackPool holds N−M = 1
    /// literal (MappedFlags=5, Positions=4 in `ude_chunk_entries`).
    fn legal_reference_dir() -> ChunkDirectory {
        let mut entries = ude_globals();
        entries.extend(ude_chunk_entries(0, 1, 100));
        entries.extend(ude_chunk_entries(0, 2, 300));
        // N−M = 5−4 = 1 per mate, so each (chunk,mate) carries a FallbackPool.
        entries.push(ude_chunk_fallback(0, 1, 1));
        entries.push(ude_chunk_fallback(0, 2, 1));
        ChunkDirectory {
            num_reads: 5,
            num_chunks: 1,
            entries,
        }
    }

    // --- single-end reference validator (`validate_reference_directory_single`) ---

    /// A minimal valid SINGLE-END reference directory: 4 globals + one chunk
    /// (mate 1 only). num_reads = Σ N_chunk = 5 reads. The mate-1 FallbackPool
    /// holds N−M = 1 literal (MappedFlags=5, Positions=4 in `ude_chunk_entries`).
    /// Note `ude_chunk_entries(_,1,_)` already emits the R1 columnar Headers + Qual.
    fn legal_single_reference_dir() -> ChunkDirectory {
        let mut entries = ude_globals();
        entries.extend(ude_chunk_entries(0, 1, 100));
        entries.push(ude_chunk_fallback(0, 1, 1));
        ChunkDirectory {
            num_reads: 5,
            num_chunks: 1,
            entries,
        }
    }

    #[test]
    fn ref_single_dir_accepts_legal() {
        validate_reference_directory_single(&legal_single_reference_dir()).unwrap();
    }

    #[test]
    fn ref_single_dir_rejects_mate2_entry() {
        // A single-end reference directory must carry NO mate-2 per-chunk roles.
        // Inject exactly one mate-2 entry so the total stays below the count cap
        // (14 legal + 1 = 15 ≤ cap bound 16) and the validator rejects on the
        // mate-1 contract, not the coarse count backstop.
        let mut dir = legal_single_reference_dir();
        dir.entries.push(ude(0, 2, StreamRole::MappedFlags, CODEC_BSC, 500, 5));
        let err = validate_reference_directory_single(&dir)
            .unwrap_err()
            .to_string();
        assert!(err.contains("mate"), "unexpected error: {err}");
    }

    #[test]
    fn ref_single_dir_rejects_missing_global() {
        let mut dir = legal_single_reference_dir();
        dir.entries
            .retain(|e| !(e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::NBitmap));
        let err = validate_reference_directory_single(&dir)
            .unwrap_err()
            .to_string();
        assert!(err.contains("missing global"), "unexpected error: {err}");
    }

    #[test]
    fn ref_single_dir_rejects_header_count_mismatch() {
        let mut dir = legal_single_reference_dir();
        for e in &mut dir.entries {
            if e.chunk_index == 0 && e.mate == 1 && e.role == StreamRole::Headers {
                e.record_count = 4; // != N_chunk (MappedFlags = 5)
            }
        }
        let err = validate_reference_directory_single(&dir)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("header record_count"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn ref_single_dir_rejects_fallback_ownership_violation() {
        let mut dir = legal_single_reference_dir();
        for e in &mut dir.entries {
            if e.chunk_index == 0 && e.mate == 1 && e.role == StreamRole::FallbackPool {
                e.record_count = 7; // != F_group (1)
            }
        }
        let err = validate_reference_directory_single(&dir)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("FallbackPool.record_count"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn ref_single_dir_rejects_sigma_n_not_num_reads() {
        let mut dir = legal_single_reference_dir();
        dir.num_reads = 99; // Σ N_chunk = 5
        let err = validate_reference_directory_single(&dir)
            .unwrap_err()
            .to_string();
        assert!(err.contains("num_reads"), "unexpected error: {err}");
    }

    #[test]
    fn ref_role_stream_role_round_trip() {
        // Every RefRole maps to a StreamRole; recovering it via (StreamRole, mate)
        // returns the original. mate disambiguates only Headers/HeaderDelta/Qual.
        let mate_of = |r: RefRole| -> u8 {
            if r.is_global() {
                0
            } else {
                match r {
                    RefRole::R2HeadersIndep | RefRole::R2HeaderDelta | RefRole::R2Qual => 2,
                    _ => 1,
                }
            }
        };
        for code in 0u8..=16 {
            let r = RefRole::from_u8(code).unwrap();
            let sr = ref_role_to_stream(r);
            let back = stream_role_to_ref(sr, mate_of(r)).unwrap();
            assert_eq!(back, r, "round-trip failed for {r:?} via {sr:?}");
        }
        // Single-end-only roles have no reference role.
        assert!(stream_role_to_ref(StreamRole::Sequence, 1).is_err());
        assert!(stream_role_to_ref(StreamRole::Nmask, 0).is_err());
        assert!(stream_role_to_ref(StreamRole::RcFlags, 1).is_err());
        // mate-1 HeaderDelta and out-of-range mate on Headers/Qual are rejected.
        assert!(stream_role_to_ref(StreamRole::HeaderDelta, 1).is_err());
        assert!(stream_role_to_ref(StreamRole::Headers, 3).is_err());
    }

    #[test]
    fn ref_dir_accepts_legal() {
        validate_reference_directory(&legal_reference_dir()).unwrap();
    }

    #[test]
    fn ref_dir_rejects_missing_global() {
        let mut dir = legal_reference_dir();
        // Drop the NBitmap global entirely.
        dir.entries
            .retain(|e| !(e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::NBitmap));
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("missing global"), "unexpected error: {err}");
    }

    #[test]
    fn ref_dir_rejects_nbitmap_count_mismatch() {
        let mut dir = legal_reference_dir();
        // NBitmap.record_count (100) must equal PackedBacking.record_count (100).
        for e in &mut dir.entries {
            if e.chunk_index == GLOBAL_SENTINEL && e.role == StreamRole::NBitmap {
                e.record_count = 99;
            }
        }
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("N-bitmap"), "unexpected error: {err}");
    }

    #[test]
    fn ref_dir_rejects_bad_codec_for_role() {
        let mut dir = legal_reference_dir();
        // Positions must be BSC; fqzcomp is illegal there.
        for e in &mut dir.entries {
            if e.role == StreamRole::Positions {
                e.codec = CODEC_FQZCOMP;
            }
        }
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("invalid codec"), "unexpected error: {err}");
    }

    #[test]
    fn ref_dir_rejects_both_r2_header() {
        let mut dir = legal_reference_dir();
        // The legal dir already has R2HeaderDelta (Headers/HeaderDelta mate 2).
        // Add an independent R2 header (Headers mate 2) so BOTH are present.
        dir.entries
            .push(ude(0, 2, StreamRole::Headers, CODEC_COLUMNAR, 5000, 5));
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(
            err.contains("exactly one R2 header role required"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn ref_dir_rejects_fallback_ownership_violation() {
        let mut dir = legal_reference_dir();
        // F_group = N−M = 1 for (chunk 0, mate 1) (N=5 MappedFlags, M=4
        // Positions); set its per-(chunk,mate) FallbackPool.record_count to a
        // wrong value.
        for e in &mut dir.entries {
            if e.chunk_index == 0 && e.mate == 1 && e.role == StreamRole::FallbackPool {
                e.record_count = 7;
            }
        }
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(
            err.contains("FallbackPool.record_count"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn ref_dir_rejects_fallback_absent_but_nonzero() {
        let mut dir = legal_reference_dir();
        // Drop (chunk 0, mate 1)'s FallbackPool entry though its F_group = 1 > 0.
        dir.entries
            .retain(|e| !(e.chunk_index == 0 && e.mate == 1 && e.role == StreamRole::FallbackPool));
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(
            err.contains("FallbackPool absent"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn ref_dir_rejects_fallback_present_but_empty() {
        // Build a directory where (chunk 0, mate 1) has N==M (no fallbacks) yet
        // still carries an (empty) FallbackPool entry. Set M-carriers to N=5 so
        // F_group = 0, then attach a record_count==0 FallbackPool.
        let mut entries = ude_globals();
        let mut chunk1 = ude_chunk_entries(0, 1, 100);
        for e in &mut chunk1 {
            if matches!(
                e.role,
                StreamRole::Positions
                    | StreamRole::Strands
                    | StreamRole::ReadLen
                    | StreamRole::EditCount
            ) {
                e.record_count = 5; // M = N = 5 ⇒ F_group = 0
            }
        }
        entries.extend(chunk1);
        entries.extend(ude_chunk_entries(0, 2, 300));
        entries.push(ude_chunk_fallback(0, 1, 0)); // present but empty (N==M)
        entries.push(ude_chunk_fallback(0, 2, 1)); // mate 2 F_group = 5−4 = 1
        let dir = ChunkDirectory {
            num_reads: 5,
            num_chunks: 1,
            entries,
        };
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("present but empty"), "unexpected error: {err}");
    }

    #[test]
    fn ref_dir_rejects_r1_header_count_mismatch() {
        let mut dir = legal_reference_dir();
        // R1 headers stream one record per read; its record_count must equal
        // N_chunk (MappedFlags = 5). A forged lower count would misalign headers
        // against the per-chunk reconstructed sequences.
        for e in &mut dir.entries {
            if e.chunk_index == 0 && e.mate == 1 && e.role == StreamRole::Headers {
                e.record_count = 4;
            }
        }
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(
            err.contains("header record_count"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn ref_dir_rejects_r1_qual_count_mismatch() {
        let mut dir = legal_reference_dir();
        for e in &mut dir.entries {
            if e.chunk_index == 0 && e.mate == 1 && e.role == StreamRole::Qual {
                e.record_count = 4;
            }
        }
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(
            err.contains("quality record_count"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn ref_dir_rejects_r2_header_count_mismatch() {
        let mut dir = legal_reference_dir();
        // mate 2's header role in the legal dir is HeaderDelta (one op per read).
        for e in &mut dir.entries {
            if e.chunk_index == 0 && e.mate == 2 && e.role == StreamRole::HeaderDelta {
                e.record_count = 6;
            }
        }
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(
            err.contains("header record_count"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn ref_dir_rejects_r2_qual_count_mismatch() {
        let mut dir = legal_reference_dir();
        for e in &mut dir.entries {
            if e.chunk_index == 0 && e.mate == 2 && e.role == StreamRole::Qual {
                e.record_count = 6;
            }
        }
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(
            err.contains("quality record_count"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn ref_dir_rejects_over_entry_cap() {
        // Entry-count cap = num_chunks·PER_CHUNK_ENTRY_CAP + GLOBAL_ENTRY_CAP and is
        // enforced up front. Declare num_chunks = 1 (bound = 1·24 + 4 = 28) but
        // ship a 2-chunk directory: 4 globals + 2·(9+9) per-chunk seq/header/qual
        // entries = 40 > 28, so the cap rejects it before any structural walk.
        let mut entries = ude_globals();
        entries.extend(ude_chunk_entries(0, 1, 100));
        entries.extend(ude_chunk_entries(0, 2, 300));
        entries.extend(ude_chunk_entries(1, 1, 500));
        entries.extend(ude_chunk_entries(1, 2, 700));
        let dir = ChunkDirectory {
            num_reads: 10,
            num_chunks: 1, // understates the 2 chunks present → 40 entries over the 28 cap
            entries,
        };
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("exceeds cap"), "unexpected error: {err}");
    }

    /// §9.3(b) M > N_chunk: make all four M-carriers (Positions/Strands/ReadLen/EditCount)
    /// exceed MappedFlags (N_chunk=5) for mate 1 so the cross-stream equality check
    /// passes first (all equal 6 == M) and the `checked_sub` underflow fires.
    #[test]
    fn ref_dir_rejects_m_greater_than_n() {
        let mut dir = legal_reference_dir();
        for e in &mut dir.entries {
            if e.chunk_index == 0
                && e.mate == 1
                && matches!(
                    e.role,
                    StreamRole::Positions
                        | StreamRole::Strands
                        | StreamRole::ReadLen
                        | StreamRole::EditCount
                )
            {
                e.record_count = 6; // M=6 > N_chunk=5
            }
        }
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(
            err.contains("M") && err.contains("N_chunk"),
            "unexpected error: {err}"
        );
    }

    /// §9.3(b) Σ N_chunk != num_reads: perturb `num_reads` so the accumulated
    /// total (5 pairs from the single legal chunk) disagrees with the header field.
    #[test]
    fn ref_dir_rejects_sigma_n_not_num_reads() {
        let mut dir = legal_reference_dir();
        dir.num_reads = 99; // Σ N_chunk = 5; num_reads = 99 → mismatch
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(err.contains("num_reads"), "unexpected error: {err}");
    }

    /// §9.3(b) cross-stream M mismatch: make Strands.record_count != Positions
    /// (M) for mate 1 so the equality check fires on that one stream.
    #[test]
    fn ref_dir_rejects_cross_stream_m_mismatch() {
        let mut dir = legal_reference_dir();
        for e in &mut dir.entries {
            if e.chunk_index == 0 && e.mate == 1 && e.role == StreamRole::Strands {
                e.record_count = 99; // Strands != M (Positions = 4)
            }
        }
        let err = validate_reference_directory(&dir).unwrap_err().to_string();
        assert!(
            err.contains("positions/strands/read_len/edit_count must all == M"),
            "unexpected error: {err}"
        );
    }
}
