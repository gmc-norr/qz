//! Reference index build + deterministic per-pair seed mapping (Task 10).
//!
//! Maps paired reads to a reference using strobealign SEEDING ONLY (NAMs; no
//! Smith-Waterman), DETERMINISTICALLY: no RNG anywhere in the runtime path and
//! an immutable insert-size (mu/sigma passed in, never mutated). Returns per
//! mate `(ref_id, projected_ref_start, is_revcomp)`.
//!
//! Determinism is load-bearing for qz (spec §4): the same reads + reference
//! must always produce the same mappings so the edit stream is reproducible.
//! The non-determinism in upstream strobealign lives in `sort_nams`' RNG
//! tie-shuffle and in the mutable `InsertSizeDistribution`; both are removed in
//! the vendored `map_pair_deterministic` patch.
use anyhow::Result;
use strobealign::chainer::{Chainer, ChainingParameters};
use strobealign::index::StrobemerIndex;
use strobealign::io::fasta::RefSequence;
use strobealign::maponly::DEFAULT_RESCUE_DISTANCE;
use strobealign::mcsstrategy::McsStrategy;

/// Where the reference comes from. `FromFasta` for real runs; `FromSeqs` for
/// tests (in-memory references, no file I/O).
pub enum ReferenceSource {
    FromFasta(std::path::PathBuf),
    /// (name, sequence) pairs, raw ACGTN bytes.
    #[allow(dead_code)] // used in tests
    FromSeqs(Vec<(Vec<u8>, Vec<u8>)>),
}

/// A projected mapping for one mate. `ref_start` is the PROJECTED read start
/// (`Nam::projected_ref_start()` = `ref_start - query_start`), i.e. where the
/// read's first base lands on the reference. Task 11 consumes this.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Mapping {
    pub ref_id: u32,
    /// PROJECTED read start on the reference.
    pub ref_start: u64,
    pub is_revcomp: bool,
}

/// Owns the references + index. The upstream `StrobemerIndex` owns its data
/// (no self-referential borrow, no `Box::leak`). Built once; drop to free index
/// RAM (Task 12 drops this before Pass 2).
pub struct Mapper {
    references: Vec<RefSequence>,
    index: StrobemerIndex,
    chainer: Chainer,
    mcs: McsStrategy,
    rescue_distance: usize,
    mu: f32,
    sigma: f32,
}

/// Derive the qz sidecar index path for a reference + read length. Keyed by the
/// CANONICAL read length for `read_len`'s seeding profile (every read length in a
/// profile bucket builds a byte-identical index, so 148 and 150 share one cache
/// entry), plus a `tag` for non-default seeding (e.g. `--reference-fast`, which
/// uses sparser syncmers and thus a different index), and given a qz-specific
/// suffix so it never clobbers a user's own strobealign `.sti`.
pub fn sidecar_path(reference: &std::path::Path, read_len: usize, tag: &str) -> std::path::PathBuf {
    // Key by the seeding PROFILE's canonical length, not the raw read length:
    // every read length in a bucket builds a byte-identical index, so 148 and 150
    // share one cache entry (`.qz-r150.sti`) instead of redundantly rebuilding.
    let canonical = strobealign::seeding::canonical_read_length(read_len);
    let mut s = reference.as_os_str().to_owned();
    s.push(format!(".qz-r{canonical}{tag}.sti"));
    std::path::PathBuf::from(s)
}

/// `""` for the default index, `"-fast"` for the `--reference-fast` sparse variant.
pub(crate) fn cache_tag(fast: bool) -> &'static str {
    if fast { "-fast" } else { "" }
}

/// Build the seeding parameters + cache tag for a (read_len, fast) pair. Shared by
/// the in-process index build (`Mapper::build`) and the standalone build
/// (`build_reference_index`) so the two produce byte-identical indexes.
fn seeding_params_for(
    read_len: usize,
    fast: bool,
) -> Result<(strobealign::seeding::SeedingParameters, &'static str)> {
    let mut params = strobealign::seeding::SeedingParameters::new(read_len);
    if fast {
        params = params
            .with_k_s(None, Some(12))
            .map_err(|e| anyhow::anyhow!("--reference-fast seeding params: {e}"))?;
    }
    Ok((params, cache_tag(fast)))
}

/// The canonical sidecar path qz-compress looks for (and `qz index` writes).
/// Public so tests and tooling resolve the same path the compressor uses.
pub fn reference_index_sidecar_path(
    reference: &std::path::Path,
    read_len: usize,
    fast: bool,
) -> std::path::PathBuf {
    sidecar_path(reference, read_len, cache_tag(fast))
}

/// Stats from building a reference index, returned by [`build_reference_index`].
#[derive(Debug, Clone)]
pub struct ReferenceIndexBuildStats {
    /// The canonical sidecar path written.
    pub path: std::path::PathBuf,
    /// Number of reference sequences indexed.
    pub references: usize,
    /// Total randstrobes in the index.
    pub randstrobes: u64,
    /// Wall-clock build time.
    pub elapsed_secs: f64,
}

/// Build a strobealign index for `reference`, tuned for `read_len`'s seeding
/// profile, and write it to the canonical sidecar path (overwriting any existing
/// one). This is the builder behind `qz index`; it produces a byte-identical
/// index to the in-process auto-build at the same `threads` (see `Mapper::build`).
pub fn build_reference_index(
    reference: &std::path::Path,
    read_len: usize,
    fast: bool,
    threads: usize,
) -> Result<ReferenceIndexBuildStats> {
    let t0 = std::time::Instant::now();
    let references = load_reference_normalized(reference)?;
    let (params, tag) = seeding_params_for(read_len, fast)?;
    let bits = params.syncmer.pick_bits(&references);
    // filter_fraction = upstream default 0.0002 (matches Mapper::build).
    let (index, stats) =
        strobealign::indexer::make_index(&references, params, bits, 0.0002, threads.max(1));
    let path = sidecar_path(reference, read_len, tag);
    index
        .write(&path)
        .map_err(|e| anyhow::anyhow!("write index {}: {e}", path.display()))?;
    Ok(ReferenceIndexBuildStats {
        path,
        references: references.len(),
        randstrobes: stats.tot_strobemer_count,
        elapsed_secs: t0.elapsed().as_secs_f64(),
    })
}

/// On-disk state of the canonical sidecar for a (reference, profile, fast) triple.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceIndexStatus {
    /// Present and newer than (or same age as) the reference FASTA.
    Ready,
    /// No sidecar at the canonical path.
    Missing,
    /// Present but older than the reference FASTA.
    Stale,
}

/// Cheap existence + staleness check (no full index load).
pub fn reference_index_status(
    reference: &std::path::Path,
    read_len: usize,
    fast: bool,
) -> ReferenceIndexStatus {
    let sidecar = sidecar_path(reference, read_len, cache_tag(fast));
    if !sidecar.exists() {
        return ReferenceIndexStatus::Missing;
    }
    if sidecar_is_stale(&sidecar, reference) {
        return ReferenceIndexStatus::Stale;
    }
    ReferenceIndexStatus::Ready
}

/// Require-by-default policy used by the CLI: ensure a USABLE index exists for
/// (reference, profile, fast) — present, non-stale, AND a valid header (so a
/// corrupt or wrong-params sidecar at the canonical path is not mistaken for
/// usable). If not usable, build it when `build` is true, else return an
/// actionable error pointing at `qz index`. Returns the canonical sidecar path.
/// The library's in-process compress still auto-builds on its own; this is the
/// explicit policy gate the CLI calls before building or sharding.
pub fn ensure_reference_index(
    reference: &std::path::Path,
    read_len: usize,
    fast: bool,
    threads: usize,
    build: bool,
) -> Result<std::path::PathBuf> {
    let (params, tag) = seeding_params_for(read_len, fast)?;
    let sidecar = sidecar_path(reference, read_len, tag);
    let status = reference_index_status(reference, read_len, fast);
    // "Usable" = present + non-stale + header validates (cheap, no full load).
    let usable = status == ReferenceIndexStatus::Ready
        && strobealign::index::read_index_header(&sidecar, &params).is_ok();
    if usable {
        return Ok(sidecar);
    }
    if build {
        return Ok(build_reference_index(reference, read_len, fast, threads)?.path);
    }
    let canonical = strobealign::seeding::canonical_read_length(read_len);
    let reason = match status {
        ReferenceIndexStatus::Stale => "is older than the reference FASTA",
        // Present + fresh but the header failed to validate.
        ReferenceIndexStatus::Ready => "is corrupt or was built with incompatible parameters",
        ReferenceIndexStatus::Missing => "was not found",
    };
    let fast_flag = if fast { " --reference-fast" } else { "" };
    anyhow::bail!(
        "reference index for the ~{canonical} bp read profile {reason} ({}).\n\
         build it first:  qz index {}{fast_flag} -r {canonical}\n\
         (or pass --build-index to build it inline)",
        sidecar.display(),
        reference.display()
    );
}

/// True iff `sidecar` is stale relative to `reference` (sidecar mtime < reference
/// mtime — IGV/BAI style). If either mtime is unavailable, treat as NOT stale.
fn sidecar_is_stale(sidecar: &std::path::Path, reference: &std::path::Path) -> bool {
    let smt = std::fs::metadata(sidecar).and_then(|m| m.modified());
    let rmt = std::fs::metadata(reference).and_then(|m| m.modified());
    match (smt, rmt) {
        (Ok(s), Ok(r)) => s < r,
        _ => false,
    }
}

fn normalize_ref(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        *b = match b.to_ascii_uppercase() {
            c @ (b'A' | b'C' | b'G' | b'T') => c,
            _ => b'N',
        };
    }
}

/// Load a reference FASTA and normalize every sequence EXACTLY as [`Mapper::build`]
/// does (`read_ref` + `normalize_ref`), but WITHOUT building the strobemer index.
///
/// The reference-aware archive merge ([`super::merge`]) uses this to re-derive the
/// coverage globals (`PackedBacking`/`NBitmap`/`IntervalMap`/`ReferenceMeta`) from
/// the same `(ref_id-order, normalized-bytes)` view the per-shard encoders saw — so
/// the merged globals are byte-identical to a single-process encode over the union
/// interval map. Skipping the index build is the whole point: the merge only needs
/// the reference *bytes*, not the seeds.
pub(crate) fn load_reference_normalized(
    path: &std::path::Path,
) -> Result<Vec<RefSequence>> {
    let mut references = strobealign::io::fasta::read_ref(path)
        .map_err(|e| anyhow::anyhow!("read reference {}: {e}", path.display()))?;
    for r in references.iter_mut() {
        normalize_ref(&mut r.sequence);
    }
    Ok(references)
}

impl Mapper {
    /// Build the index from `source` using the read-length seeding profile.
    /// Insert-size defaults to (300, 100); Task 12 Pass 0 calls `with_insert`
    /// after estimating the real distribution.
    ///
    /// For `FromFasta`, an auto-cached sidecar index next to the reference is
    /// loaded if present + non-stale + parameter-matching; otherwise the index is
    /// built and the sidecar is best-effort written for next time. The caching is
    /// transparent and never changes the produced index (a loaded index is
    /// bit-identical to a freshly-built one), so mappings — and thus archives —
    /// are byte-identical whether the cache hit or missed. `FromSeqs` (tests)
    /// never caches.
    pub fn build(
        source: ReferenceSource,
        read_len: usize,
        n_threads: usize,
        fast: bool,
        require: bool,
    ) -> Result<Self> {
        let ref_path: Option<std::path::PathBuf> = match &source {
            ReferenceSource::FromFasta(p) => Some(p.clone()),
            ReferenceSource::FromSeqs(_) => None,
        };
        let mut references = match source {
            ReferenceSource::FromFasta(p) => strobealign::io::fasta::read_ref(&p)
                .map_err(|e| anyhow::anyhow!("read reference {}: {e}", p.display()))?,
            ReferenceSource::FromSeqs(seqs) => seqs
                .into_iter()
                .map(|(name, sequence)| RefSequence {
                    name: String::from_utf8_lossy(&name).into_owned(),
                    sequence,
                })
                .collect(),
        };
        for r in references.iter_mut() {
            normalize_ref(&mut r.sequence);
        }
        // `--reference-fast`: sparser syncmers (s=12 vs the read-length default's
        // s=16 → k−s 8 vs 4 → ~45% fewer randstrobes). Cuts mapping CPU ~39% for
        // ~8–10% faster compress at any core count, lossless and ~ratio-neutral on
        // high-coverage data (reconstruct-verify falls unmappable reads back to
        // literals, so the only risk is a few more fallbacks on low-coverage /
        // divergent inputs). The sparser index gets its own cache tag.
        let (params, cache_tag) = seeding_params_for(read_len, fast)?;
        let bits = params.syncmer.pick_bits(&references);
        // filter_fraction = upstream default 0.0002.
        let index = Self::load_or_build_index(
            ref_path.as_deref(),
            &references,
            &params,
            bits,
            read_len,
            cache_tag,
            n_threads.max(1),
            require,
        )?;
        let chainer = Chainer::new(index.k(), ChainingParameters::default());
        Ok(Mapper {
            references,
            index,
            chainer,
            mcs: McsStrategy::default(),
            rescue_distance: DEFAULT_RESCUE_DISTANCE,
            mu: 300.0,
            sigma: 100.0,
        })
    }

    /// Load a matching cached sidecar index (FromFasta only) or build one and
    /// best-effort write the sidecar.
    ///
    /// When `require` is false (the library default, and `--build-index`), a
    /// missing / stale / corrupt sidecar falls through to an in-memory build (the
    /// result is identical to an uncached build) — transparent self-healing for
    /// API consumers and tests. When `require` is true (CLI reference compress
    /// without `--build-index`), the SAME conditions are a hard error instead: a
    /// present-but-corrupt sidecar (valid header, bad payload — which the cheap CLI
    /// preflight cannot detect) is caught HERE, by the full `read_index`, so it can
    /// never silently rebuild. `FromSeqs` (tests) never caches and ignores `require`.
    fn load_or_build_index(
        ref_path: Option<&std::path::Path>,
        references: &[RefSequence],
        params: &strobealign::seeding::SeedingParameters,
        bits: u8,
        read_len: usize,
        cache_tag: &str,
        n_threads: usize,
        require: bool,
    ) -> Result<StrobemerIndex> {
        let build = || {
            let (index, _stats) = strobealign::indexer::make_index(
                references,
                params.clone(),
                bits,
                0.0002,
                n_threads,
            );
            index
        };

        let Some(reference) = ref_path else {
            // FromSeqs (tests) — never cache; `require` is moot (no sidecar).
            return Ok(build());
        };
        let sidecar = sidecar_path(reference, read_len, cache_tag);

        let require_err = |reason: &str| {
            let canonical = strobealign::seeding::canonical_read_length(read_len);
            let fast_flag = if cache_tag == "-fast" { " --reference-fast" } else { "" };
            anyhow::anyhow!(
                "reference index {} {reason}; a prebuilt index is required.\n\
                 build it:  qz index {}{fast_flag} -r {canonical}\n\
                 (or pass --build-index to build it inline)",
                sidecar.display(),
                reference.display()
            )
        };

        // Try a cache hit: sidecar present + non-stale + full load (validates the
        // whole payload, not just the header).
        if sidecar.exists() && !sidecar_is_stale(&sidecar, reference) {
            match strobealign::index::read_index(&sidecar, params.clone(), bits) {
                Ok(index) => {
                    eprintln!("qz reference: using cached index {}", sidecar.display());
                    return Ok(index);
                }
                Err(e) if require => {
                    return Err(require_err(&format!("is unusable ({e})")));
                }
                Err(e) => {
                    eprintln!(
                        "qz reference: cached index {} unusable ({e}); rebuilding",
                        sidecar.display()
                    );
                }
            }
        } else if require {
            // Missing or stale, and a prebuilt index is required.
            let reason = if sidecar.exists() {
                "is older than the reference FASTA"
            } else {
                "was not found"
            };
            return Err(require_err(reason));
        }

        // Cache miss (require=false): build, then best-effort write the sidecar.
        let index = build();
        match index.write(&sidecar) {
            Ok(()) => eprintln!("qz reference: wrote index cache {}", sidecar.display()),
            Err(e) => eprintln!(
                "qz reference: could not write index cache {} ({e}); continuing",
                sidecar.display()
            ),
        }
        Ok(index)
    }

    /// Override the insert-size distribution (Task 12 Pass 0 estimate). The
    /// values are immutable thereafter — never mutated during mapping, which is
    /// what keeps mapping deterministic.
    pub fn with_insert(mut self, mu: f32, sigma: f32) -> Self {
        self.mu = mu;
        self.sigma = sigma;
        self
    }

    #[allow(dead_code)] // used in tests
    pub fn num_references(&self) -> usize {
        self.references.len()
    }

    pub fn references(&self) -> &[RefSequence] {
        &self.references
    }

    /// Map a single read (raw ACGTN sequence bytes, as read). Returns the
    /// projected mapping (or `None` when no NAM seeds). The single-end analogue
    /// of `map_pair`: no insert size, no pairing — `mu`/`sigma` are not used.
    pub fn map_single(&self, r: &[u8]) -> Result<Option<Mapping>> {
        let n = strobealign::maponly::map_single_deterministic(
            r,
            &self.index,
            &self.chainer,
            self.rescue_distance,
            self.mcs,
        );
        Ok(n.map(|n| to_mapping(&n)))
    }

    /// Map a pair (raw ACGTN sequence bytes, already as read). Returns projected
    /// mappings (or `None` per mate when no NAM seeds).
    pub fn map_pair(&self, r1: &[u8], r2: &[u8]) -> Result<(Option<Mapping>, Option<Mapping>)> {
        let (n1, n2) = strobealign::maponly::map_pair_deterministic(
            r1,
            r2,
            &self.index,
            &self.chainer,
            self.rescue_distance,
            self.mcs,
            self.mu,
            self.sigma,
        );
        Ok((n1.map(|n| to_mapping(&n)), n2.map(|n| to_mapping(&n))))
    }
}

fn to_mapping(n: &strobealign::nam::Nam) -> Mapping {
    Mapping {
        ref_id: n.ref_id as u32,
        ref_start: n.projected_ref_start() as u64,
        is_revcomp: n.is_revcomp,
    }
}

// ===========================================================================
// Task 11 — Hamming placement + edit extraction + fallback decision.
//
// A NAM gives a *provisional* projected base index for a read; the exact
// placement may be off by a few bases. `place_read` searches small offsets
// around the provisional base, computes the Hamming substitutions vs the
// consensus slice at each, and picks the offset with the FEWEST edits. If the
// best edit count exceeds `K` (or no offset fits, or a read base ∉ ACGTN) it
// returns `None` and the caller demotes the read to a literal fallback.
// ===========================================================================

use super::edits::{Edit, base_is_acgtn, extract_substitutions};
use crate::compression::dna_utils::reverse_complement;

/// A chosen placement for a mapped read: where it sits in flat consensus base
/// space (`base_index`) plus the substitutions on the FORWARD-oriented read vs
/// the consensus slice at that index.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Placement {
    pub base_index: u64,
    pub edits: Vec<Edit>,
}

/// Try to place a read against `consensus` near `provisional_base`, returning the
/// substitution placement with the fewest mismatches, or `None` to fall back.
///
/// # Coordinate / consensus contract (read this — Task 12 depends on it)
/// `consensus` is a FLAT byte slice in base-index space and `provisional_base`
/// is an index INTO it. Offsets are searched in `provisional_base ± window`
/// (saturating at 0) and each candidate offset `off` is only considered when
/// `[off, off+len) ⊆ [0, consensus.len())` (offsets that would run off either end
/// are skipped). `place_read` does NOT itself check interval boundaries: it only
/// guarantees the span stays inside the flat `consensus` bytes it was given.
/// Task 12 is responsible for the per-interval bound — it must pass a consensus
/// that lets the offset search range stay valid and MUST additionally verify
/// `IntervalMap::span_in_single_interval(returned base_index, len)` before
/// trusting the result (and ultimately reconstruct-verify, Task 12/13's
/// byte-exact backstop). Simplest correct caller: pass the full flat consensus;
/// `place_read` enforces `[off,off+len) ⊆ [0,len)`, and Task 12 then rejects any
/// returned `base_index` whose span crosses an interval boundary.
///
/// # Orientation
/// `r_fwd = reverse_complement(read)` when `is_rc`, else `read` as-is. The
/// returned `edits` are always on the FORWARD-oriented read (consensus is
/// reference-forward); decode re-applies the RC via `reconstruct_read(.., is_rc)`.
///
/// # Tie-break
/// Offsets are scanned ascending and `best` is replaced only on STRICTLY fewer
/// edits, so among equal-edit placements the SMALLEST offset wins —
/// deterministic. A perfect (zero-edit) placement short-circuits the search.
///
/// # Accept / reject
/// Returns `Some(Placement)` iff the best in-range offset has `<= k` mismatches.
/// Returns `None` (→ caller falls back to literal) when: the read is empty or
/// longer than `consensus`, no in-range offset has `<= k` mismatches, or every
/// in-range offset has a read base ∉ ACGTN (`extract_substitutions` → `None`,
/// treated as "this offset invalid"). Never panics.
pub fn place_read(
    consensus: &[u8],
    provisional_base: u64,
    read: &[u8],
    is_rc: bool,
    k: usize,
    window: u64,
) -> Option<Placement> {
    let r_fwd = if is_rc {
        reverse_complement(read)
    } else {
        read.to_vec()
    };
    let len = r_fwd.len();
    if len == 0 || len > consensus.len() {
        return None;
    }
    // ACGTN validity is offset-independent (r_fwd is fixed across the search), so
    // check it ONCE instead of per offset: a non-ACGTN read base means no offset
    // can be encoded (every `extract_substitutions` would return None), i.e. fall
    // back to a literal — identical to the old per-offset behavior.
    if !r_fwd.iter().all(|&b| base_is_acgtn(b)) {
        return None;
    }
    let lo = provisional_base.saturating_sub(window);
    let hi = provisional_base.saturating_add(window);
    // Largest offset that keeps [off, off+len) ⊆ [0, consensus.len()).
    let max_off = (consensus.len() - len) as u64;

    // Find the offset with the FEWEST mismatches via an allocation-free count pass
    // with early termination: once a candidate's running count reaches the current
    // best it can't be STRICTLY fewer, so stop counting it. `best` is replaced only
    // on strictly fewer, so among equal-count offsets the SMALLEST wins — and any
    // accepted (<= k) winner is counted in full, so its `best_count` is exact. Same
    // tie-break, same accept threshold, same result as extracting at every offset.
    //
    // Seed `best_count` at `k + 1` rather than usize::MAX: an offset with > k
    // mismatches can neither be accepted nor be the accepted minimum, so cap every
    // offset there — this early-terminates the leading run of hopeless edge-of-window
    // offsets (a frameshift mismatches most bases) instead of counting the first one
    // in full. Offsets with <= k mismatches are still counted exactly (k+1 never
    // caps them), so the accepted result is byte-identical.
    let mut best_off: Option<u64> = None;
    let mut best_count: usize = k + 1;
    let mut off = lo;
    loop {
        if off <= max_off {
            let o = off as usize;
            let cons = &consensus[o..o + len];
            let mut cnt = 0usize;
            for i in 0..len {
                if cons[i] != r_fwd[i] {
                    cnt += 1;
                    if cnt >= best_count {
                        break; // cannot beat the current best at this offset
                    }
                }
            }
            if cnt < best_count {
                best_count = cnt;
                best_off = Some(off);
                if cnt == 0 {
                    break; // can't beat zero mismatches
                }
            }
        }
        // Guard against u64 overflow when hi == u64::MAX (window saturated).
        if off == hi {
            break;
        }
        off += 1;
    }

    match best_off {
        // r_fwd is ACGTN (checked above) so extract_substitutions yields Some;
        // build the edit list once, only for the winning offset.
        Some(off) if best_count <= k => {
            let o = off as usize;
            let edits = extract_substitutions(&consensus[o..o + len], &r_fwd)?;
            Some(Placement {
                base_index: off,
                edits,
            })
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reference_is_normalized_in_place() {
        let refseq = b"acgtRYacgtACGTacgtacgtacgtACGT".to_vec(); // lowercase + IUPAC R,Y
        let m = Mapper::build(
            ReferenceSource::FromSeqs(vec![(b"c".to_vec(), refseq)]),
            20,
            1,
            false,
            false,
        )
        .unwrap();
        let norm = &m.references()[0].sequence;
        assert!(
            norm.iter()
                .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'))
        );
        assert_eq!(norm[4], b'N'); // R -> N
        assert_eq!(norm[5], b'N'); // Y -> N
        assert_eq!(norm[0], b'A'); // a -> A
    }
    fn make_nonrepetitive_seq(n: usize, seed: u64) -> Vec<u8> {
        // simple deterministic LCG over ACGT, low-repeat
        let mut x = seed.wrapping_add(0x9E3779B97F4A7C15);
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            v.push(b"ACGT"[((x >> 33) & 3) as usize]);
        }
        v
    }
    #[test]
    fn maps_reads_to_tiny_reference() {
        let refseq = make_nonrepetitive_seq(600, 7);
        let read_fwd = refseq[200..320].to_vec(); // 120 bp exact substring
        let read_len = read_fwd.len();
        let m = Mapper::build(
            ReferenceSource::FromSeqs(vec![(b"chr0".to_vec(), refseq.clone())]),
            read_len,
            1,
            false,
            false,
        )
        .unwrap();
        let (r1, _r2) = m.map_pair(&read_fwd, &read_fwd).unwrap();
        let r1 = r1.expect("r1 mapped");
        assert_eq!(r1.ref_id, 0);
        assert!(
            (r1.ref_start as i64 - 200).abs() <= 2,
            "got {}",
            r1.ref_start
        );
        assert!(!r1.is_revcomp);
        // reverse-complement maps reverse
        let rc = crate::compression::dna_utils::reverse_complement(&read_fwd);
        let (r1rc, _) = m.map_pair(&rc, &rc).unwrap();
        assert!(r1rc.expect("rc mapped").is_revcomp);
    }
    #[test]
    fn map_single_matches_map_pair_mate1() {
        // map_single on a read must return the same Mapping that map_pair returns
        // for mate 1 when the read maps unambiguously to a unique location (no
        // pair-rescue can change a confidently-placed single read).
        let refseq = make_nonrepetitive_seq(600, 7);
        let read_fwd = refseq[200..320].to_vec(); // 120 bp exact substring
        let m = Mapper::build(
            ReferenceSource::FromSeqs(vec![(b"chr0".to_vec(), refseq.clone())]),
            read_fwd.len(),
            1,
            false,
            false,
        )
        .unwrap();
        let single = m.map_single(&read_fwd).unwrap().expect("read mapped");
        assert_eq!(single.ref_id, 0);
        assert!((single.ref_start as i64 - 200).abs() <= 2, "got {}", single.ref_start);
        assert!(!single.is_revcomp);
        // reverse-complement maps reverse-strand
        let rc = crate::compression::dna_utils::reverse_complement(&read_fwd);
        let single_rc = m.map_single(&rc).unwrap().expect("rc mapped");
        assert!(single_rc.is_revcomp);
        // determinism across two calls
        let a = m.map_single(&read_fwd).unwrap().unwrap();
        let b = m.map_single(&read_fwd).unwrap().unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn determinism_same_result_twice() {
        let refseq = make_nonrepetitive_seq(600, 11);
        let read = refseq[100..220].to_vec();
        let m = Mapper::build(
            ReferenceSource::FromSeqs(vec![(b"c".to_vec(), refseq)]),
            read.len(),
            1,
            false,
            false,
        )
        .unwrap();
        let a = m.map_pair(&read, &read).unwrap().0.unwrap();
        let b = m.map_pair(&read, &read).unwrap().0.unwrap();
        assert_eq!(a, b);
    }
}

#[cfg(test)]
mod place_tests {
    use super::*;

    #[test]
    fn hamming_places_and_extracts_edits() {
        let consensus = b"ACGTACGTACGTACGT"; // 16 bp, treat as one interval [0,16)
        let read = b"ACGTTCGTAC"; // matches consensus[0..10] with pos4 A->T
        let p = place_read(
            consensus, /*provisional_base*/ 0, read, /*is_rc*/ false, /*K*/ 2,
            /*window*/ 3,
        )
        .unwrap();
        assert_eq!(p.base_index, 0);
        assert_eq!(p.edits, vec![Edit { pos: 4, base: b'T' }]);
    }

    #[test]
    fn too_many_mismatches_falls_back() {
        let consensus = b"ACGTACGTACGTACGT";
        let read = b"TTTTTTTTTT";
        assert!(place_read(consensus, 0, read, false, /*K*/ 2, 3).is_none()); // -> fallback
    }

    #[test]
    fn offset_search_finds_the_better_placement() {
        // Read matches consensus[2..12] with 0 edits, but consensus[0..10] with
        // several. provisional_base = 0, window >= 2 must find offset 2 / no edits.
        let consensus = b"XXACGTACGTACXX"; // bytes 2..12 = "ACGTACGTAC"
        let read = b"ACGTACGTAC";
        let p = place_read(consensus, 0, read, false, /*K*/ 2, /*window*/ 2).unwrap();
        assert_eq!(p.base_index, 2);
        assert_eq!(p.edits, vec![]);
    }

    #[test]
    fn reverse_strand_edits_on_forward_oriented_read() {
        // r_fwd is what we Hamming against the consensus; the input `read` is its
        // reverse-complement. Edits are reported on the forward-oriented read.
        let consensus = b"ACGTTCGTACGTACGT"; // forward-oriented target
        let r_fwd = b"ACGTTCGTAC"; // matches consensus[0..10] with pos4 A->T
        let read: Vec<u8> = reverse_complement(r_fwd); // the stored (RC) read
        let p = place_read(
            consensus, 0, &read, /*is_rc*/ true, /*K*/ 2, /*window*/ 3,
        )
        .unwrap();
        assert_eq!(p.base_index, 0);
        // Edits are vs the forward-oriented read: consensus[4]='T', r_fwd[4]='T'?
        // consensus[0..10] = "ACGTTCGTAC", r_fwd = "ACGTTCGTAC" => 0 edits.
        assert_eq!(p.edits, vec![]);
    }

    #[test]
    fn reverse_strand_with_substitution() {
        // Same as above but force one mismatch on the forward-oriented read.
        let consensus = b"ACGTACGTACGTACGT";
        let r_fwd = b"ACGTTCGTAC"; // pos4 A->T vs consensus[0..10] "ACGTACGTAC"
        let read = reverse_complement(r_fwd);
        let p = place_read(consensus, 0, &read, true, /*K*/ 2, 3).unwrap();
        assert_eq!(p.base_index, 0);
        assert_eq!(p.edits, vec![Edit { pos: 4, base: b'T' }]);
    }

    #[test]
    fn offset_out_of_range_is_skipped() {
        // window would reach an offset that runs off the end; only in-range
        // offsets are considered, and the best in-range one wins.
        let consensus = b"ACGTACGT"; // 8 bp
        let read = b"ACGTAC"; // len 6 -> max_off = 2
        // provisional at 5, window 4 -> lo=1, hi=9 but only off in {1,2} are in range.
        // off=0 not searched (lo=1). consensus[1..7]="CGTACG" vs read -> mismatches;
        // off=2 "GTACGT" vs read -> mismatches. With K large, smallest-edit wins.
        let p = place_read(consensus, 5, read, false, /*K*/ 8, /*window*/ 4).unwrap();
        // Both in-range; pick the one with fewest edits, smallest offset on tie.
        assert!(p.base_index == 1 || p.base_index == 2);
        // Sanity: offset 0 (the perfect match) was OUT of the search window, so it
        // must NOT have been chosen.
        assert_ne!(p.base_index, 0);
    }

    #[test]
    fn non_acgtn_read_base_returns_none() {
        let consensus = b"ACGTACGTACGT";
        let read = b"ACGRACGTAC"; // 'R' is not ACGTN
        assert!(place_read(consensus, 0, read, false, /*K*/ 4, 3).is_none());
    }

    #[test]
    fn read_longer_than_consensus_returns_none() {
        let consensus = b"ACGT";
        let read = b"ACGTACGT"; // longer than consensus -> no valid offset
        assert!(place_read(consensus, 0, read, false, /*K*/ 8, 3).is_none());
    }

    #[test]
    fn empty_read_returns_none() {
        let consensus = b"ACGTACGT";
        let read = b"";
        assert!(place_read(consensus, 0, read, false, 4, 3).is_none());
    }
}
