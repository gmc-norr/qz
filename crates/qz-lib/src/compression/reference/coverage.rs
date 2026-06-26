//! Per-reference coverage bitmap + gap-merge coalescing (single-pass compress).
use super::backing::Interval;

/// One flat `Vec<u64>` bit-array with per-contig base offsets: 1 bit/base
/// (~375 MB for the human genome, bounded by genome size, not read count).
pub struct CoverageMap {
    words: Vec<u64>,
    contig_base: Vec<u64>,
    lens: Vec<u64>,
}

impl CoverageMap {
    pub fn new(contig_lens: &[u64]) -> Self {
        let mut contig_base = Vec::with_capacity(contig_lens.len() + 1);
        let mut total = 0u64;
        for &l in contig_lens {
            contig_base.push(total);
            total += l;
        }
        contig_base.push(total);
        Self {
            words: vec![0u64; (total as usize).div_ceil(64)],
            contig_base,
            lens: contig_lens.to_vec(),
        }
    }

    /// Set bits [start, start+len) on contig `ref_id`. Mutated in the SERIAL
    /// in-order fold (no atomics). Clamp to the contig length defensively.
    ///
    /// Word-aligned: fills whole `u64` words in the interior and masks only the
    /// two boundary words, so a length-L span costs ~O(L/64) word ops instead of
    /// L per-bit index/shift/OR ops. Sets EXACTLY the same bits as the old naive
    /// per-bit loop (proven by `mark_is_word_aligned_equivalent_to_naive_bit_set`),
    /// so `finalize` — which reads only the final bit state — is byte-identical.
    pub fn mark(&mut self, ref_id: u32, start: u64, len: u64) {
        let clen = self.lens[ref_id as usize];
        let s = start.min(clen);
        let e = (start + len).min(clen);
        if s >= e {
            return;
        }
        let base = self.contig_base[ref_id as usize];
        let lo = base + s; // first bit to set (inclusive)
        let hi = base + e; // one past the last bit to set (exclusive)
        let first = (lo / 64) as usize;
        let last = ((hi - 1) / 64) as usize; // word holding the last set bit (hi-1)
        let lo_b = lo % 64; // < 64
        // hi_b = one past the last set bit within `last`, in 1..=64.
        let hi_b = ((hi - 1) % 64) + 1;
        if first == last {
            // All set bits live in one word: [lo_b, hi_b), with hi_b > lo_b.
            let high = if hi_b == 64 { u64::MAX } else { (1u64 << hi_b) - 1 };
            self.words[first] |= high & (u64::MAX << lo_b);
            return;
        }
        // First (partial) word: bits [lo_b, 64).
        self.words[first] |= u64::MAX << lo_b;
        // Interior full words.
        for w in (first + 1)..last {
            self.words[w] = u64::MAX;
        }
        // Last (partial) word: bits [0, hi_b).
        let last_mask = if hi_b == 64 { u64::MAX } else { (1u64 << hi_b) - 1 };
        self.words[last] |= last_mask;
    }

    /// Coalesce covered runs separated by gaps < g into intervals, in
    /// (ref_id, ref_start) order. A merged interval only grows, so it still
    /// contains every contributing span.
    pub fn finalize(self, g: u64) -> Vec<Interval> {
        let mut out = Vec::new();
        for ref_id in 0..self.lens.len() {
            let base = self.contig_base[ref_id];
            let clen = self.lens[ref_id];
            // collect maximal covered runs within [0, clen) for this contig
            let mut run_start: Option<u64> = None;
            // (start, end) of the accumulating interval being built
            let mut cur: Option<(u64, u64)> = None;

            for pos in 0..clen {
                let bit = base + pos;
                let set = (self.words[(bit / 64) as usize] >> (bit % 64)) & 1 == 1;
                match (set, run_start) {
                    (true, None) => run_start = Some(pos),
                    (false, Some(rs)) => {
                        // close run [rs, pos)
                        let (rstart, rend) = (rs, pos);
                        run_start = None;
                        // merge into cur if gap < g, otherwise flush cur and start new
                        let merged = match cur {
                            Some((_, e)) if rstart - e < g => {
                                cur = Some((cur.unwrap().0, rend));
                                true
                            }
                            _ => false,
                        };
                        if !merged {
                            if let Some((s, e)) = cur.take() {
                                out.push(Interval {
                                    ref_id: ref_id as u32,
                                    ref_start: s,
                                    length: e - s,
                                });
                            }
                            cur = Some((rstart, rend));
                        }
                    }
                    _ => {}
                }
            }
            // close a run that reaches the contig end
            if let Some(rs) = run_start {
                let (rstart, rend) = (rs, clen);
                let merged = match cur {
                    Some((_, e)) if rstart - e < g => {
                        cur = Some((cur.unwrap().0, rend));
                        true
                    }
                    _ => false,
                };
                if !merged {
                    if let Some((s, e)) = cur.take() {
                        out.push(Interval {
                            ref_id: ref_id as u32,
                            ref_start: s,
                            length: e - s,
                        });
                    }
                    cur = Some((rstart, rend));
                }
            }
            // flush the final interval for this contig
            if let Some((s, e)) = cur.take() {
                out.push(Interval {
                    ref_id: ref_id as u32,
                    ref_start: s,
                    length: e - s,
                });
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coalesces_gaps_below_g_and_splits_above() {
        let mut cov = CoverageMap::new(&[300]); // one contig, len 300
        cov.mark(0, 0, 50);
        cov.mark(0, 60, 40); // gap of 10 (50->60) before this span
        cov.mark(0, 250, 20); // gap of 150 before this span
        let ivs = cov.finalize(16);
        assert_eq!(
            ivs,
            vec![
                Interval {
                    ref_id: 0,
                    ref_start: 0,
                    length: 100
                }, // first two merge (gap 10 < 16)
                Interval {
                    ref_id: 0,
                    ref_start: 250,
                    length: 20
                }, // third separate
            ]
        );
    }

    #[test]
    fn covering_interval_contains_every_marked_span() {
        let mut cov = CoverageMap::new(&[100]);
        cov.mark(0, 10, 20);
        let ivs = cov.finalize(0);
        assert!(
            ivs.iter()
                .any(|iv| iv.ref_start <= 10 && 30 <= iv.ref_start + iv.length)
        );
    }

    /// The word-aligned `mark` must set the EXACT same bits as the old naive
    /// per-bit loop, for every alignment + boundary case: spans straddling u64
    /// words, spans crossing contig bases that are not 64-aligned, single-bit
    /// spans, full-word spans, a 1-base contig, and clamped over-runs (start or
    /// start+len past the contig end). This is the byte-identity gate: `finalize`
    /// reads only the final bit state, so identical `words` ⇒ identical archive.
    #[test]
    fn mark_is_word_aligned_equivalent_to_naive_bit_set() {
        // Lengths chosen so contig bases land off 64-bit boundaries (200, 330,
        // 394, 395, 1395 …) and include a 1-base contig.
        let lens = [200u64, 130, 64, 1, 1000, 63, 65];
        let total: u64 = lens.iter().sum();
        let mut contig_base = Vec::new();
        {
            let mut t = 0u64;
            for &l in &lens {
                contig_base.push(t);
                t += l;
            }
        }

        // Deterministic LCG (no rng dependency, no Math.random in tests).
        let mut state = 0x9E3779B97F4A7C15u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            state >> 11
        };

        let mut cov = CoverageMap::new(&lens);
        let mut naive = vec![0u64; (total as usize).div_ceil(64)];

        for _ in 0..20_000 {
            let ref_id = (next() % lens.len() as u64) as usize;
            let clen = lens[ref_id];
            let start = next() % (clen + 5); // sometimes >= clen (clamp path)
            let len = next() % (clen + 9); // sometimes overruns the contig end
            cov.mark(ref_id as u32, start, len);
            // Replicate the OLD per-bit semantics exactly.
            let s = start.min(clen);
            let e = (start + len).min(clen);
            let base = contig_base[ref_id];
            for bit in (base + s)..(base + e) {
                naive[(bit / 64) as usize] |= 1u64 << (bit % 64);
            }
        }
        assert_eq!(
            cov.words, naive,
            "word-aligned mark must set exactly the same bits as the naive per-bit loop"
        );
    }

    #[test]
    fn multiple_contigs_emit_in_ref_id_order() {
        let mut cov = CoverageMap::new(&[50, 50]);
        cov.mark(1, 5, 10);
        cov.mark(0, 0, 5);
        let ivs = cov.finalize(0);
        assert_eq!(
            ivs,
            vec![
                Interval {
                    ref_id: 0,
                    ref_start: 0,
                    length: 5
                },
                Interval {
                    ref_id: 1,
                    ref_start: 5,
                    length: 10
                },
            ]
        );
    }
}
