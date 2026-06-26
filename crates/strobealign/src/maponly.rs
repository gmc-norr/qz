use crate::chainer::Chainer;
use crate::details::NamDetails;
use crate::index::StrobemerIndex;
use crate::math::normal_pdf;
use crate::mcsstrategy::McsStrategy;
use crate::nam::{Nam, get_nams_by_chaining};

/// Default rescue distance for mapping (upstream `MappingParameters` default).
/// Exposed here because the full `MappingParameters` (a SAM-mapper config) was
/// removed with the alignment stack; this is the only field the map-only path uses.
pub const DEFAULT_RESCUE_DISTANCE: usize = 100;

/// Deterministic, RNG-free, immutable-insert-size paired map-only.
/// Replicates map_paired_end_read's NAM path but (1) replaces sort_nams' rng
/// tie-shuffle with a fixed total-order sort, and (2) takes immutable mu/sigma
/// (no InsertSizeDistribution mutation). Returns the best NAM per mate (or None).
/// qz determinism (spec §4) depends on this. [VENDOR PATCH — record in VENDORED.md]
pub fn map_pair_deterministic(
    r1_seq: &[u8],
    r2_seq: &[u8],
    index: &StrobemerIndex,
    chainer: &Chainer,
    rescue_distance: usize,
    mcs_strategy: McsStrategy,
    mu: f32,
    sigma: f32,
) -> (Option<Nam>, Option<Nam>) {
    let (d1, mut nams1) =
        get_nams_by_chaining(r1_seq, index, chainer, rescue_distance, mcs_strategy);
    let (d2, mut nams2) =
        get_nams_by_chaining(r2_seq, index, chainer, rescue_distance, mcs_strategy);

    if nams1.is_empty() && nams2.is_empty() {
        return (None, None);
    }

    // Pairing (deterministic two-pointer scan; immutable mu/sigma). Note this
    // sorts the orientation-partitioned sub-slices by (ref_id, projected_start),
    // not the whole slices.
    let nam_pairs = get_nam_pairs(&mut nams1, &mut nams2, mu, sigma, &d1, &d2);

    // Replaces sort_nams(.., rng): score desc, then a fixed total order on
    // (ref_id, ref_start, is_revcomp) to break score ties WITHOUT an rng shuffle.
    deterministic_sort(&mut nams1);
    deterministic_sort(&mut nams2);

    // Best-location selection: inline port of get_best_paired_mapping_location
    // WITHOUT the insert_size_distribution.update() side effect. Prefer the best
    // proper pair only if its score beats the (penalized) individual mapping;
    // else return the best individual NAM per mate. nam_pairs is already sorted
    // by score desc; nams{1,2} are sorted so .first() is the best individual.
    let best_nam1 = nams1.first();
    let best_nam2 = nams2.first();
    let individual_score = best_nam1.map_or(0.0, |nam| nam.score as f64)
        + best_nam2.map_or(0.0, |nam| nam.score as f64);

    if let Some(NamPair { nam1, nam2, score }) = nam_pairs.first()
        && *score >= individual_score / 2.0
    {
        (Some(nam1.clone()), Some(nam2.clone()))
    } else {
        (best_nam1.cloned(), best_nam2.cloned())
    }
}

/// Fixed total-order sort replacing sort_nams' rng tie-shuffle: score desc,
/// then (ref_id, ref_start, is_revcomp) ascending to break score ties.
fn deterministic_sort(nams: &mut [Nam]) {
    nams.sort_by(|a, b| {
        b.score.total_cmp(&a.score).then_with(|| {
            (a.ref_id, a.ref_start, a.is_revcomp).cmp(&(b.ref_id, b.ref_start, b.is_revcomp))
        })
    });
}

#[derive(Debug)]
pub struct NamPair {
    pub nam1: Nam,
    pub nam2: Nam,
    pub score: f64,
}

/// Build all plausible forward/revcomp mapping pairings
fn get_nam_pairs(
    nams1: &mut [Nam],
    nams2: &mut [Nam],
    mu: f32,
    sigma: f32,
    details1: &NamDetails,
    details2: &NamDetails,
) -> Vec<NamPair> {
    let mut nam_pairs = vec![];
    if nams1.is_empty() || nams2.is_empty() {
        return nam_pairs;
    }

    let (fwd1, rev1): (&mut [Nam], &mut [Nam]) =
        split_nams_by_orientation_checked(nams1, details1.both_orientations);
    let (fwd2, rev2): (&mut [Nam], &mut [Nam]) =
        split_nams_by_orientation_checked(nams2, details2.both_orientations);

    if !fwd1.is_empty() && !rev2.is_empty() {
        fwd1.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        rev2.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        nam_pairs.extend(find_pairs(fwd1, rev2, mu, sigma, false));
    }
    if !fwd2.is_empty() && !rev1.is_empty() {
        fwd2.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        rev1.sort_unstable_by_key(|nam| (nam.ref_id, nam.projected_ref_start()));
        nam_pairs.extend(find_pairs(fwd2, rev1, mu, sigma, true));
    }

    nam_pairs.sort_unstable_by(|a, b| b.score.total_cmp(&a.score));
    nam_pairs
}

/// Split nams into (forward, revcomp),
/// if only 1 orientation exists, returns it and a empty slice
fn split_nams_by_orientation_checked(nams: &mut [Nam], both: bool) -> (&mut [Nam], &mut [Nam]) {
    if both {
        split_nams_by_orientation(nams)
    } else if nams[0].is_revcomp {
        (&mut [], nams)
    } else {
        (nams, &mut [])
    }
}

/// In-place partition of NAMs by orientation:
/// forward on the left, revcomp on the right.
/// Returns two slices separating (forward, revcomp)
fn split_nams_by_orientation(nams: &mut [Nam]) -> (&mut [Nam], &mut [Nam]) {
    let mut left = 0;
    let mut right = nams.len();

    while left < right {
        if nams[left].is_revcomp {
            right -= 1;
            nams.swap(left, right);
        } else {
            left += 1;
        }
    }

    nams.split_at_mut(left)
}

/// Deterministic, RNG-free single-end map-only — the per-read half of
/// `map_pair_deterministic` WITHOUT pairing and WITHOUT mu/sigma. Replicates the
/// single-mate NAM path: `get_nams_by_chaining` → `deterministic_sort` (fixed
/// total order, no rng tie-shuffle) → best NAM (or `None`). Used by qz's
/// single-end reference mode, where reproducible mappings are required.
/// [VENDOR PATCH — recorded in VENDORED.md]
pub fn map_single_deterministic(
    seq: &[u8],
    index: &StrobemerIndex,
    chainer: &Chainer,
    rescue_distance: usize,
    mcs_strategy: McsStrategy,
) -> Option<Nam> {
    let (_d, mut nams) =
        get_nams_by_chaining(seq, index, chainer, rescue_distance, mcs_strategy);
    deterministic_sort(&mut nams);
    nams.first().cloned()
}

/// Find most forward/revcomp pairs using a two-pointer scan.
/// Assumes both slices are sorted by (ref_id, projected_ref_start).
fn find_pairs(fwd: &[Nam], rev: &[Nam], mu: f32, sigma: f32, swap_order: bool) -> Vec<NamPair> {
    let mut out = Vec::new();
    let max_dist = (mu + 10.0 * sigma).ceil() as usize; // distance cutoff from insert size distribution
    let mut rev_ptr = 0;
    let mut last_paired = None;

    for f in fwd {
        // Advance revcomp pointer to the first possible candidate
        while rev_ptr < rev.len()
            && (rev[rev_ptr].ref_id < f.ref_id
                || rev[rev_ptr].ref_id == f.ref_id
                    && rev[rev_ptr].projected_ref_start() < f.projected_ref_start())
        {
            rev_ptr += 1;
        }
        if rev_ptr == rev.len() {
            break;
        }
        if rev[rev_ptr].ref_id > f.ref_id {
            continue;
        }

        // Scan window of revcomp nams within distance limit.
        let mut best = None;
        let mut i = rev_ptr;
        while i < rev.len()
            && rev[i].ref_id == f.ref_id
            && (rev[i].projected_ref_start() - f.projected_ref_start()) <= max_dist
        {
            let r = &rev[i];
            // The pairing score gets a bonus based on the reference distance of the two chosen nam
            // paired and from our current knowledge of the reference distance distribution
            let x = f.ref_start.abs_diff(r.ref_start);
            let score = f.score as f64
                + r.score as f64
                + 0.001f64.max((normal_pdf(x as f32, mu, sigma) + 1.0).ln() as f64);

            if best.is_none_or(|(_, highest_score)| score > highest_score) {
                best = Some((i, score));
            }
            i += 1;
        }

        // Highest scoring candidate
        let Some((best_id, score)) = best else {
            continue;
        };
        let r = &rev[best_id];

        // If the same revcomp nam was paired previously, keep only the better scoring pair.
        if let Some((last_id, prev_score)) = last_paired
            && last_id == best_id
        {
            if score <= prev_score {
                continue; // keep the previous better pair
            }
            out.pop(); // replace it 
        }

        out.push(if swap_order {
            NamPair {
                nam1: r.clone(),
                nam2: f.clone(),
                score,
            }
        } else {
            NamPair {
                nam1: f.clone(),
                nam2: r.clone(),
                score,
            }
        });

        last_paired = Some((best_id, score));
        rev_ptr = best_id;
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::StrobemerIndex;
    use crate::indexer::make_index;
    use crate::io::fasta::RefSequence;
    use crate::seeding::SeedingParameters;

    /// Deterministic low-repeat ACGT generator (mirrors mapping.rs tests).
    fn make_nonrepetitive_seq(n: usize, seed: u64) -> Vec<u8> {
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

    /// Build a real index over one reference sequence at the given read length.
    fn build_index(refseq: &[u8], read_len: usize) -> (StrobemerIndex, Chainer) {
        let references = vec![RefSequence {
            name: "chr0".to_string(),
            sequence: refseq.to_vec(),
        }];
        let params = SeedingParameters::new(read_len);
        let bits = params.syncmer.pick_bits(&references);
        let (index, _stats) = make_index(&references, params, bits, 0.0002, 1);
        let chainer = Chainer::new(index.k(), crate::chainer::ChainingParameters::default());
        (index, chainer)
    }

    #[test]
    fn map_single_maps_a_known_read() {
        let refseq = make_nonrepetitive_seq(600, 7);
        let read = refseq[200..320].to_vec(); // 120 bp exact substring
        let (index, chainer) = build_index(&refseq, read.len());
        let nam = map_single_deterministic(
            &read,
            &index,
            &chainer,
            DEFAULT_RESCUE_DISTANCE,
            McsStrategy::default(),
        )
        .expect("read should map");
        assert_eq!(nam.ref_id, 0);
        assert!(
            (nam.projected_ref_start() as i64 - 200).abs() <= 2,
            "projected start {} not near 200",
            nam.projected_ref_start()
        );
        assert!(!nam.is_revcomp);
    }

    #[test]
    fn map_single_is_deterministic() {
        let refseq = make_nonrepetitive_seq(600, 11);
        let read = refseq[100..220].to_vec();
        let (index, chainer) = build_index(&refseq, read.len());
        let a = map_single_deterministic(
            &read,
            &index,
            &chainer,
            DEFAULT_RESCUE_DISTANCE,
            McsStrategy::default(),
        )
        .unwrap();
        let b = map_single_deterministic(
            &read,
            &index,
            &chainer,
            DEFAULT_RESCUE_DISTANCE,
            McsStrategy::default(),
        )
        .unwrap();
        assert_eq!(a.ref_id, b.ref_id);
        assert_eq!(a.ref_start, b.ref_start);
        assert_eq!(a.is_revcomp, b.is_revcomp);
        assert_eq!(a.score, b.score);
    }
}
