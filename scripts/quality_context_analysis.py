#!/usr/bin/env python3
"""
Analyze correlations between sequence context, position, and Phred quality scores.

Questions:
1. How much does position alone predict quality?
2. How much does the local sequence context predict quality?
3. How much does (position + context) jointly predict quality?
4. What's the residual entropy after conditioning on these predictors?
"""
import sys
import numpy as np
from collections import defaultdict
import math

BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def parse_fastq(path, max_reads=500000):
    sequences = []
    qualities = []
    with open(path) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            qual = f.readline().strip()
            if len(seq) == len(qual) and 'N' not in seq:
                sequences.append(seq)
                qualities.append(qual)
            if len(sequences) >= max_reads:
                break
    return sequences, qualities

def entropy(counts_dict):
    """Shannon entropy of a distribution given as {value: count}."""
    total = sum(counts_dict.values())
    if total == 0:
        return 0.0
    h = 0.0
    for c in counts_dict.values():
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h

def conditional_entropy(joint_counts):
    """H(Y|X) given joint_counts[x][y] = count."""
    total = sum(sum(yd.values()) for yd in joint_counts.values())
    h = 0.0
    for x, yd in joint_counts.items():
        nx = sum(yd.values())
        if nx == 0:
            continue
        hx = entropy(yd)
        h += (nx / total) * hx
    return h

def main():
    if len(sys.argv) < 2:
        print("Usage: quality_context_analysis.py <fastq_file> [max_reads]")
        sys.exit(1)

    path = sys.argv[1]
    max_reads = int(sys.argv[2]) if len(sys.argv) > 2 else 200000

    print(f"Loading up to {max_reads} reads from {path}...")
    sequences, qualities = parse_fastq(path, max_reads)
    n_reads = len(sequences)
    read_len = len(sequences[0])
    print(f"Loaded {n_reads:,} reads, length {read_len}\n")

    # ================================================================
    # 1. Unconditional quality distribution
    # ================================================================
    print("=== 1. Unconditional Quality Distribution ===\n")
    qual_counts = defaultdict(int)
    total_bases = 0
    for qual in qualities:
        for q in qual:
            phred = ord(q) - 33
            qual_counts[phred] += 1
            total_bases += 1

    h_uncond = entropy(qual_counts)
    print(f"Total bases: {total_bases:,}")
    print(f"Unique quality values: {len(qual_counts)}")
    print(f"Unconditional entropy H(Q): {h_uncond:.4f} bits/base")
    print(f"  -> {total_bases * h_uncond / 8:,.0f} bytes minimum")
    print()

    # Show distribution
    for q in sorted(qual_counts.keys()):
        c = qual_counts[q]
        pct = 100 * c / total_bases
        bar = '#' * int(pct * 2)
        if pct > 0.5:
            print(f"  Q{q:>2}: {pct:>6.2f}%  {bar}")

    # ================================================================
    # 2. Quality conditioned on position
    # ================================================================
    print("\n=== 2. Quality Conditioned on Position ===\n")
    pos_qual = defaultdict(lambda: defaultdict(int))  # pos -> {q: count}
    for qual in qualities:
        for j, q in enumerate(qual):
            phred = ord(q) - 33
            pos_qual[j][phred] += 1

    h_pos = conditional_entropy(pos_qual)
    mi_pos = h_uncond - h_pos
    print(f"H(Q|position): {h_pos:.4f} bits/base")
    print(f"Mutual info I(Q;position): {mi_pos:.4f} bits/base")
    print(f"Reduction: {100*mi_pos/h_uncond:.1f}%")
    print(f"  -> {total_bases * h_pos / 8:,.0f} bytes after position conditioning")
    print()

    # Position-wise mean quality
    print("Position-wise mean quality (every 10th position):")
    for j in range(0, read_len, 10):
        qs = pos_qual[j]
        total_j = sum(qs.values())
        mean_q = sum(q * c for q, c in qs.items()) / total_j
        std_q = math.sqrt(sum((q - mean_q)**2 * c for q, c in qs.items()) / total_j)
        h_j = entropy(qs)
        print(f"  pos {j:>3}: mean={mean_q:.1f}, std={std_q:.1f}, H={h_j:.2f} bits")

    # ================================================================
    # 3. Quality conditioned on current base
    # ================================================================
    print("\n=== 3. Quality Conditioned on Current Base ===\n")
    base_qual = defaultdict(lambda: defaultdict(int))
    for seq, qual in zip(sequences, qualities):
        for j, (b, q) in enumerate(zip(seq, qual)):
            if b in BASE_MAP:
                phred = ord(q) - 33
                base_qual[b][phred] += 1

    h_base = conditional_entropy(base_qual)
    mi_base = h_uncond - h_base
    print(f"H(Q|base): {h_base:.4f} bits/base")
    print(f"Mutual info I(Q;base): {mi_base:.4f} bits/base")
    print(f"Reduction: {100*mi_base/h_uncond:.1f}%")
    print()
    for b in 'ACGT':
        qs = base_qual[b]
        total_b = sum(qs.values())
        mean_q = sum(q * c for q, c in qs.items()) / total_b
        print(f"  Base {b}: mean_Q={mean_q:.2f}, n={total_b:,}")

    # ================================================================
    # 4. Quality conditioned on (position, current base)
    # ================================================================
    print("\n=== 4. Quality Conditioned on (Position, Current Base) ===\n")
    posbase_qual = defaultdict(lambda: defaultdict(int))
    for seq, qual in zip(sequences, qualities):
        for j, (b, q) in enumerate(zip(seq, qual)):
            if b in BASE_MAP:
                phred = ord(q) - 33
                key = (j, b)
                posbase_qual[key][phred] += 1

    h_posbase = conditional_entropy(posbase_qual)
    mi_posbase = h_uncond - h_posbase
    print(f"H(Q|position,base): {h_posbase:.4f} bits/base")
    print(f"Mutual info I(Q;position,base): {mi_posbase:.4f} bits/base")
    print(f"Reduction: {100*mi_posbase/h_uncond:.1f}%")

    # ================================================================
    # 5. Quality conditioned on k-mer context (previous k bases)
    # ================================================================
    print("\n=== 5. Quality Conditioned on Previous k-mer Context ===\n")
    for k in [1, 2, 3, 4, 5]:
        ctx_qual = defaultdict(lambda: defaultdict(int))
        for seq, qual in zip(sequences, qualities):
            for j in range(k, len(seq)):
                context = seq[j-k:j]
                if all(c in BASE_MAP for c in context):
                    phred = ord(qual[j]) - 33
                    ctx_qual[context][phred] += 1

        h_ctx = conditional_entropy(ctx_qual)
        mi_ctx = h_uncond - h_ctx
        n_contexts = len(ctx_qual)
        print(f"  k={k}: H(Q|prev_{k}mer)={h_ctx:.4f}, MI={mi_ctx:.4f} bits "
              f"({100*mi_ctx/h_uncond:.1f}%), {n_contexts} contexts")

    # ================================================================
    # 6. Quality conditioned on (position + previous k-mer)
    # ================================================================
    print("\n=== 6. Quality Conditioned on (Position + Previous k-mer) ===\n")
    for k in [1, 2, 3]:
        posctx_qual = defaultdict(lambda: defaultdict(int))
        for seq, qual in zip(sequences, qualities):
            for j in range(k, len(seq)):
                context = seq[j-k:j]
                if all(c in BASE_MAP for c in context):
                    phred = ord(qual[j]) - 33
                    key = (j, context)
                    posctx_qual[key][phred] += 1

        h_posctx = conditional_entropy(posctx_qual)
        mi_posctx = h_uncond - h_posctx
        n_contexts = len(posctx_qual)
        print(f"  k={k}: H(Q|pos,prev_{k}mer)={h_posctx:.4f}, MI={mi_posctx:.4f} bits "
              f"({100*mi_posctx/h_uncond:.1f}%), {n_contexts} contexts")

    # ================================================================
    # 7. Quality conditioned on (position + surrounding context)
    # ================================================================
    print("\n=== 7. Quality Conditioned on (Position + Surrounding k-mer) ===\n")
    for left, right in [(1,1), (2,1), (1,2), (2,2), (3,2)]:
        surr_qual = defaultdict(lambda: defaultdict(int))
        for seq, qual in zip(sequences, qualities):
            for j in range(left, len(seq) - right):
                context = seq[j-left:j] + seq[j+1:j+1+right]
                if all(c in BASE_MAP for c in context) and seq[j] in BASE_MAP:
                    phred = ord(qual[j]) - 33
                    key = (j, context)
                    surr_qual[key][phred] += 1

        h_surr = conditional_entropy(surr_qual)
        mi_surr = h_uncond - h_surr
        n_contexts = len(surr_qual)
        print(f"  L={left},R={right}: H(Q|pos,surr)={h_surr:.4f}, MI={mi_surr:.4f} bits "
              f"({100*mi_surr/h_uncond:.1f}%), {n_contexts} contexts")

    # ================================================================
    # 8. Quality conditioned on previous quality + position
    # ================================================================
    print("\n=== 8. Quality Conditioned on Previous Quality ===\n")
    for prev_n in [1, 2]:
        prevq_qual = defaultdict(lambda: defaultdict(int))
        for qual in qualities:
            for j in range(prev_n, len(qual)):
                prev_quals = tuple(ord(qual[j-i-1]) - 33 for i in range(prev_n))
                phred = ord(qual[j]) - 33
                prevq_qual[prev_quals][phred] += 1

        h_prevq = conditional_entropy(prevq_qual)
        mi_prevq = h_uncond - h_prevq
        print(f"  prev_{prev_n}_qual: H(Q|prev_Q)={h_prevq:.4f}, MI={mi_prevq:.4f} bits "
              f"({100*mi_prevq/h_uncond:.1f}%)")

    # Prev quality + position
    prevq_pos_qual = defaultdict(lambda: defaultdict(int))
    for qual in qualities:
        for j in range(1, len(qual)):
            prev_q = ord(qual[j-1]) - 33
            phred = ord(qual[j]) - 33
            key = (j, prev_q)
            prevq_pos_qual[key][phred] += 1

    h_prevqpos = conditional_entropy(prevq_pos_qual)
    mi_prevqpos = h_uncond - h_prevqpos
    print(f"\n  H(Q|pos,prev_Q): {h_prevqpos:.4f}, MI={mi_prevqpos:.4f} bits "
          f"({100*mi_prevqpos/h_uncond:.1f}%)")

    # ================================================================
    # 9. Combined: position + prev_quality + sequence context
    # ================================================================
    print("\n=== 9. Combined Predictors ===\n")
    combined_qual = defaultdict(lambda: defaultdict(int))
    for seq, qual in zip(sequences, qualities):
        for j in range(1, len(seq)):
            if seq[j-1] in BASE_MAP and seq[j] in BASE_MAP:
                prev_q = ord(qual[j-1]) - 33
                prev_base = seq[j-1]
                cur_base = seq[j]
                phred = ord(qual[j]) - 33
                # Quantize position to bins of 5
                pos_bin = j // 5
                key = (pos_bin, prev_q, prev_base, cur_base)
                combined_qual[key][phred] += 1

    h_combined = conditional_entropy(combined_qual)
    mi_combined = h_uncond - h_combined
    n_contexts = len(combined_qual)
    print(f"H(Q|pos_bin5, prev_Q, prev_base, cur_base): {h_combined:.4f} bits")
    print(f"Mutual info: {mi_combined:.4f} bits ({100*mi_combined/h_uncond:.1f}%)")
    print(f"Contexts: {n_contexts:,}")
    print(f"  -> {total_bases * h_combined / 8:,.0f} bytes")
    print()

    # Most aggressive: position + prev 2 quals + prev 2 bases
    combined2_qual = defaultdict(lambda: defaultdict(int))
    for seq, qual in zip(sequences, qualities):
        for j in range(2, len(seq)):
            if all(seq[j-i] in BASE_MAP for i in range(3)):
                prev_q1 = ord(qual[j-1]) - 33
                prev_q2 = ord(qual[j-2]) - 33
                ctx = seq[j-2:j+1]  # 2 prev + current base
                phred = ord(qual[j]) - 33
                pos_bin = j // 5
                key = (pos_bin, prev_q1, prev_q2, ctx)
                combined2_qual[key][phred] += 1

    h_combined2 = conditional_entropy(combined2_qual)
    mi_combined2 = h_uncond - h_combined2
    n_contexts2 = len(combined2_qual)
    print(f"H(Q|pos_bin5, prev_2Q, 3mer_ctx): {h_combined2:.4f} bits")
    print(f"Mutual info: {mi_combined2:.4f} bits ({100*mi_combined2/h_uncond:.1f}%)")
    print(f"Contexts: {n_contexts2:,}")
    print(f"  -> {total_bases * h_combined2 / 8:,.0f} bytes")

    # ================================================================
    # Summary
    # ================================================================
    print("\n=== SUMMARY ===\n")
    print(f"{'Predictor':<45} {'H(Q|X)':>8} {'MI':>8} {'Red%':>6} {'Est bytes':>12}")
    print("-" * 85)
    results = [
        ("Unconditional", h_uncond, 0, 0),
        ("Position only", h_pos, mi_pos, mi_pos/h_uncond*100),
        ("Current base only", h_base, mi_base, mi_base/h_uncond*100),
        ("Position + base", h_posbase, mi_posbase, mi_posbase/h_uncond*100),
        ("Previous 2-mer context", None, None, None),  # placeholder
        ("Position + prev 1-mer", None, None, None),
        ("Previous quality (1)", None, None, None),
        ("Position + prev quality", h_prevqpos, mi_prevqpos, mi_prevqpos/h_uncond*100),
        ("Pos_bin5 + prevQ + prev_base + cur_base", h_combined, mi_combined, mi_combined/h_uncond*100),
        ("Pos_bin5 + prev2Q + 3mer_ctx", h_combined2, mi_combined2, mi_combined2/h_uncond*100),
    ]
    for name, h, mi, red in results:
        if h is not None:
            est_bytes = int(total_bases * h / 8)
            print(f"{name:<45} {h:>8.4f} {mi:>8.4f} {red:>5.1f}% {est_bytes:>12,}")

if __name__ == "__main__":
    main()
