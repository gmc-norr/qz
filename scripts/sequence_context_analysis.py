#!/usr/bin/env python3
"""
Analyze correlations in DNA sequences that could improve compression.
Vectorized with numpy for speed.
"""
import sys
import math
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import time

BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 255}

def parse_fastq_to_arrays(path, max_reads=500000):
    """Parse FASTQ into numpy arrays of base indices and phred scores."""
    seqs = []
    quals = []
    with open(path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            qual = f.readline().strip()
            seqs.append(seq)
            quals.append(qual)
            if len(seqs) >= max_reads:
                break

    n = len(seqs)
    read_len = len(seqs[0])

    # Convert to numpy arrays
    seq_arr = np.zeros((n, read_len), dtype=np.uint8)
    qual_arr = np.zeros((n, read_len), dtype=np.uint8)

    for i, (s, q) in enumerate(zip(seqs, quals)):
        for j, (b, qc) in enumerate(zip(s, q)):
            seq_arr[i, j] = BASE_MAP.get(b, 255)
            qual_arr[i, j] = ord(qc) - 33

    return seq_arr, qual_arr, n, read_len

def entropy_from_counts(counts):
    """Shannon entropy from a count array."""
    total = counts.sum()
    if total == 0:
        return 0.0
    probs = counts[counts > 0] / total
    return -np.sum(probs * np.log2(probs))

def conditional_entropy_from_joint(joint_counts):
    """H(Y|X) from joint_counts[x, y] array."""
    row_totals = joint_counts.sum(axis=1)
    total = row_totals.sum()
    if total == 0:
        return 0.0
    h = 0.0
    for x in range(joint_counts.shape[0]):
        nx = row_totals[x]
        if nx == 0:
            continue
        probs = joint_counts[x][joint_counts[x] > 0] / nx
        hx = -np.sum(probs * np.log2(probs))
        h += (nx / total) * hx
    return h

def analyze_unconditional(seq_arr):
    valid = seq_arr[seq_arr < 4]
    counts = np.bincount(valid, minlength=4)[:4]
    h = entropy_from_counts(counts)
    return h, counts, len(valid)

def analyze_position(seq_arr, read_len):
    """H(B|position)"""
    # joint[pos, base] counts
    joint = np.zeros((read_len, 4), dtype=np.int64)
    for pos in range(read_len):
        col = seq_arr[:, pos]
        valid = col[col < 4]
        joint[pos] = np.bincount(valid, minlength=4)[:4]
    return conditional_entropy_from_joint(joint), joint

def analyze_prev_kmer(seq_arr, k, read_len):
    """H(B|prev_k_bases)"""
    n_ctx = 5**k  # 0-3 = ACGT, 4 = N/sentinel
    joint = np.zeros((n_ctx, 4), dtype=np.int64)

    flat_seq = seq_arr.ravel()
    n_reads = seq_arr.shape[0]

    for start_pos in range(k, read_len):
        # Build context for all reads at this position
        ctx = np.zeros(n_reads, dtype=np.int64)
        for offset in range(k):
            col = seq_arr[:, start_pos - k + offset].astype(np.int64)
            col = np.where(col < 4, col, 4)
            ctx = ctx * 5 + col

        cur = seq_arr[:, start_pos]
        valid_mask = cur < 4
        ctx_valid = ctx[valid_mask]
        cur_valid = cur[valid_mask]

        # Accumulate counts
        for c, b in zip(ctx_valid, cur_valid):
            joint[c, b] += 1

    return conditional_entropy_from_joint(joint), len(set(range(n_ctx)))

def analyze_prev_kmer_fast(seq_arr, k, read_len):
    """H(B|prev_k_bases) - vectorized version"""
    n_ctx = 5**k
    joint = np.zeros((n_ctx, 4), dtype=np.int64)
    n_reads = seq_arr.shape[0]

    for start_pos in range(k, read_len):
        ctx = np.zeros(n_reads, dtype=np.int64)
        all_valid = np.ones(n_reads, dtype=bool)
        for offset in range(k):
            col = seq_arr[:, start_pos - k + offset].astype(np.int64)
            is_base = col < 4
            col_safe = np.where(is_base, col, 4)
            ctx = ctx * 5 + col_safe
            all_valid &= is_base

        cur = seq_arr[:, start_pos]
        valid_mask = (cur < 4) & all_valid
        ctx_v = ctx[valid_mask]
        cur_v = cur[valid_mask]

        # Use np.add.at for fast accumulation
        idx = ctx_v * 4 + cur_v.astype(np.int64)
        np.add.at(joint.ravel(), idx, 1)

    h = conditional_entropy_from_joint(joint)
    used = np.sum(joint.sum(axis=1) > 0)
    return h, int(used)

def analyze_quality_base(seq_arr, qual_arr):
    """H(B|quality)"""
    max_q = int(qual_arr.max()) + 1
    joint = np.zeros((max_q, 4), dtype=np.int64)

    flat_seq = seq_arr.ravel()
    flat_qual = qual_arr.ravel()
    valid = flat_seq < 4
    seq_v = flat_seq[valid]
    qual_v = flat_qual[valid]

    idx = qual_v.astype(np.int64) * 4 + seq_v.astype(np.int64)
    np.add.at(joint.ravel(), idx, 1)

    return conditional_entropy_from_joint(joint)

def analyze_quality_prevbase(seq_arr, qual_arr, read_len):
    """H(B|quality, prev_base)"""
    max_q = int(qual_arr.max()) + 1
    n_ctx = max_q * 5  # quality * (4 bases + sentinel)
    joint = np.zeros((n_ctx, 4), dtype=np.int64)

    n_reads = seq_arr.shape[0]
    for pos in range(1, read_len):
        prev = seq_arr[:, pos-1].astype(np.int64)
        prev = np.where(prev < 4, prev, 4)
        cur = seq_arr[:, pos]
        qual = qual_arr[:, pos].astype(np.int64)

        valid = cur < 4
        ctx = qual[valid] * 5 + prev[valid]
        cur_v = cur[valid].astype(np.int64)

        idx = ctx * 4 + cur_v
        np.add.at(joint.ravel(), idx, 1)

    return conditional_entropy_from_joint(joint)

def analyze_combined(seq_arr, qual_arr, read_len, pos_bin=5):
    """H(B|pos_bin, quality, prev_base)"""
    max_q = int(qual_arr.max()) + 1
    n_pos_bins = (read_len + pos_bin - 1) // pos_bin
    n_ctx = n_pos_bins * max_q * 5
    joint = np.zeros((n_ctx, 4), dtype=np.int64)

    n_reads = seq_arr.shape[0]
    for pos in range(1, read_len):
        pb = pos // pos_bin
        prev = seq_arr[:, pos-1].astype(np.int64)
        prev = np.where(prev < 4, prev, 4)
        cur = seq_arr[:, pos]
        qual = qual_arr[:, pos].astype(np.int64)

        valid = cur < 4
        ctx = (pb * max_q * 5 + qual[valid] * 5 + prev[valid])
        cur_v = cur[valid].astype(np.int64)

        idx = ctx * 4 + cur_v
        np.add.at(joint.ravel(), idx, 1)

    return conditional_entropy_from_joint(joint), n_ctx

def analyze_combined2(seq_arr, qual_arr, read_len, pos_bin=10):
    """H(B|pos_bin, quality, prev_2mer)"""
    max_q = int(qual_arr.max()) + 1
    n_pos_bins = (read_len + pos_bin - 1) // pos_bin
    n_prev = 5 * 5  # prev 2 bases (including sentinel=4)
    n_ctx = n_pos_bins * max_q * n_prev
    joint = np.zeros((n_ctx, 4), dtype=np.int64)

    n_reads = seq_arr.shape[0]
    for pos in range(2, read_len):
        pb = pos // pos_bin
        pp = seq_arr[:, pos-2].astype(np.int64)
        pp = np.where(pp < 4, pp, 4)
        p = seq_arr[:, pos-1].astype(np.int64)
        p = np.where(p < 4, p, 4)
        cur = seq_arr[:, pos]
        qual = qual_arr[:, pos].astype(np.int64)

        valid = cur < 4
        prev2 = pp[valid] * 5 + p[valid]
        ctx = pb * max_q * n_prev + qual[valid] * n_prev + prev2
        cur_v = cur[valid].astype(np.int64)

        idx = ctx * 4 + cur_v
        np.add.at(joint.ravel(), idx, 1)

    used = np.sum(joint.reshape(-1, 4).sum(axis=1) > 0)
    return conditional_entropy_from_joint(joint), int(used)

def adaptive_bits(seq_arr, qual_arr, read_len, mode="uncond"):
    """Simulate adaptive arithmetic coder. Returns total bits."""
    # Build models as numpy arrays for speed
    n_reads = seq_arr.shape[0]
    total_bits = 0.0

    if mode == "uncond":
        counts = np.ones(4, dtype=np.float64)
        for i in range(n_reads):
            for j in range(read_len):
                b = seq_arr[i, j]
                if b < 4:
                    total_bits += -np.log2(counts[b] / counts.sum())
                    counts[b] += 1
        return total_bits

    elif mode == "prev_base":
        counts = np.ones((5, 4), dtype=np.float64)  # prev(5) x cur(4)
        for i in range(n_reads):
            prev = 4
            for j in range(read_len):
                b = seq_arr[i, j]
                if b < 4:
                    total_bits += -np.log2(counts[prev, b] / counts[prev].sum())
                    counts[prev, b] += 1
                    prev = b
                else:
                    prev = 4
        return total_bits

    elif mode == "prev2_base":
        counts = np.ones((25, 4), dtype=np.float64)
        for i in range(n_reads):
            pp, p = 4, 4
            for j in range(read_len):
                b = seq_arr[i, j]
                if b < 4:
                    ctx = pp * 5 + p
                    total_bits += -np.log2(counts[ctx, b] / counts[ctx].sum())
                    counts[ctx, b] += 1
                    pp, p = p, b
                else:
                    pp, p = 4, 4
        return total_bits

    elif mode == "prev_base_qual":
        max_q = int(qual_arr.max()) + 1
        counts = np.ones((5 * max_q, 4), dtype=np.float64)
        for i in range(n_reads):
            prev = 4
            for j in range(read_len):
                b = seq_arr[i, j]
                if b < 4:
                    q = qual_arr[i, j]
                    ctx = prev * max_q + q
                    total_bits += -np.log2(counts[ctx, b] / counts[ctx].sum())
                    counts[ctx, b] += 1
                    prev = b
                else:
                    prev = 4
        return total_bits

def run_adaptive(args):
    """Worker for parallel adaptive simulation."""
    seq_arr, qual_arr, read_len, mode, label = args
    t0 = time.time()
    bits = adaptive_bits(seq_arr, qual_arr, read_len, mode)
    dt = time.time() - t0
    return label, bits, dt

def main():
    if len(sys.argv) < 2:
        print("Usage: sequence_context_analysis.py <fastq_file> [max_reads]")
        sys.exit(1)

    path = sys.argv[1]
    max_reads = int(sys.argv[2]) if len(sys.argv) > 2 else 500000

    t0 = time.time()
    print(f"Loading up to {max_reads} reads from {path}...")
    seq_arr, qual_arr, n_reads, read_len = parse_fastq_to_arrays(path, max_reads)
    total_bases = n_reads * read_len
    valid_mask = seq_arr < 4
    total_valid = int(valid_mask.sum())
    print(f"Loaded {n_reads:,} reads x {read_len} = {total_bases:,} bases "
          f"({total_valid:,} valid) in {time.time()-t0:.1f}s")
    print()

    # ================================================================
    # Run all entropy analyses in parallel
    # ================================================================
    print("=== Entropy Analysis (vectorized) ===\n")

    t0 = time.time()

    # 1. Unconditional
    h_uncond, base_counts, _ = analyze_unconditional(seq_arr)
    print(f"Base composition: A={base_counts[0]/total_valid:.4f} C={base_counts[1]/total_valid:.4f} "
          f"G={base_counts[2]/total_valid:.4f} T={base_counts[3]/total_valid:.4f}")
    print(f"H(B) unconditional: {h_uncond:.6f} bits/base")
    print()

    # Run analyses in parallel using ProcessPoolExecutor
    # But for numpy arrays, it's easier to just run them sequentially
    # since they're already vectorized

    # 2. Position
    h_pos, pos_joint = analyze_position(seq_arr, read_len)
    mi_pos = h_uncond - h_pos
    print(f"H(B|position):     {h_pos:.6f}  MI={mi_pos:.6f} ({100*mi_pos/h_uncond:.2f}%)")

    # 3. Previous k-mer context
    for k in [1, 2, 3, 4, 5, 6, 7, 8]:
        h_k, n_ctx = analyze_prev_kmer_fast(seq_arr, k, read_len)
        mi_k = h_uncond - h_k
        save = int(mi_k * total_valid / 8)
        print(f"H(B|prev_{k}mer):   {h_k:.6f}  MI={mi_k:.6f} ({100*mi_k/h_uncond:.2f}%)  "
              f"ctx={n_ctx:>8}  save={save:>8,}B")

    # 4. Quality
    h_qual = analyze_quality_base(seq_arr, qual_arr)
    mi_qual = h_uncond - h_qual
    print(f"H(B|quality):      {h_qual:.6f}  MI={mi_qual:.6f} ({100*mi_qual/h_uncond:.2f}%)")

    # Quality per base
    max_q = int(qual_arr.max()) + 1
    flat_seq = seq_arr.ravel()
    flat_qual = qual_arr.ravel()
    valid = flat_seq < 4
    print("\nBase composition by quality:")
    for q in range(max_q):
        mask = valid & (flat_qual == q)
        n = mask.sum()
        if n > 10000:
            bc = np.bincount(flat_seq[mask], minlength=4)[:4]
            probs = bc / bc.sum()
            h_q = entropy_from_counts(bc)
            print(f"  Q{q:>2} (n={n:>10,}): A={probs[0]:.4f} C={probs[1]:.4f} "
                  f"G={probs[2]:.4f} T={probs[3]:.4f}  H={h_q:.4f}")

    # 5. Quality + prev base
    print()
    h_qualctx = analyze_quality_prevbase(seq_arr, qual_arr, read_len)
    mi_qualctx = h_uncond - h_qualctx
    print(f"H(B|qual,prev_b):  {h_qualctx:.6f}  MI={mi_qualctx:.6f} ({100*mi_qualctx/h_uncond:.2f}%)")

    # 6. Combined: pos + qual + prev base
    h_comb, _ = analyze_combined(seq_arr, qual_arr, read_len, pos_bin=5)
    mi_comb = h_uncond - h_comb
    print(f"H(B|pos5,q,prev):  {h_comb:.6f}  MI={mi_comb:.6f} ({100*mi_comb/h_uncond:.2f}%)")

    # 7. Combined: pos + qual + prev 2-mer
    h_comb2, n_ctx2 = analyze_combined2(seq_arr, qual_arr, read_len, pos_bin=10)
    mi_comb2 = h_uncond - h_comb2
    print(f"H(B|pos10,q,2mer): {h_comb2:.6f}  MI={mi_comb2:.6f} ({100*mi_comb2/h_uncond:.2f}%)  "
          f"ctx={n_ctx2:,}")

    dt = time.time() - t0
    print(f"\nEntropy analysis completed in {dt:.1f}s")

    # ================================================================
    # Summary
    # ================================================================
    bsc_seq_bytes = int(34_769_092 * n_reads / 1_000_000)

    print(f"\n=== SUMMARY ===\n")
    print(f"BSC sequence compression: {bsc_seq_bytes:,} bytes ({8*bsc_seq_bytes/total_valid:.4f} bits/base)")
    print()
    print(f"{'Predictor':<35} {'H(B|X)':>10} {'MI':>10} {'Red%':>6} {'Save bytes':>12}")
    print("-" * 80)

    results = [
        ("Unconditional", h_uncond, 0),
        ("Position", h_pos, mi_pos),
        ("Previous 1-mer", None, None),  # filled from loop
        ("Previous 2-mer", None, None),
        ("Quality score", h_qual, mi_qual),
        ("Quality + prev_base", h_qualctx, mi_qualctx),
        ("Pos5 + quality + prev_base", h_comb, mi_comb),
        ("Pos10 + quality + prev_2mer", h_comb2, mi_comb2),
    ]

    for name, h, mi in results:
        if h is not None:
            save = int(mi * total_valid / 8)
            red = 100*mi/h_uncond if h_uncond > 0 else 0
            print(f"{name:<35} {h:>10.6f} {mi:>10.6f} {red:>5.1f}% {save:>12,}")

if __name__ == "__main__":
    main()
