#!/usr/bin/env python3
"""
How quickly does the quality context model converge?
Track bits/base as a function of reads processed.
"""
import sys
import math
import numpy as np
import time

BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}

def parse_fastq(path, max_reads=1000000):
    records = []
    with open(path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            qual = f.readline().strip()
            records.append((seq, [ord(q)-33 for q in qual]))
            if len(records) >= max_reads:
                break
    return records

def main():
    if len(sys.argv) < 2:
        print("Usage: quality_model_convergence.py <fastq_file> [max_reads]")
        sys.exit(1)

    path = sys.argv[1]
    max_reads = int(sys.argv[2]) if len(sys.argv) > 2 else 1000000

    print(f"Loading up to {max_reads} reads...")
    records = parse_fastq(path, max_reads)
    n_reads = len(records)
    read_len = len(records[0][1])
    print(f"Loaded {n_reads:,} reads, length {read_len}\n")

    # Find quality alphabet
    all_quals = set()
    for _, phreds in records:
        all_quals.update(phreds)
    qual_list = sorted(all_quals)
    qual_to_idx = {q: i for i, q in enumerate(qual_list)}
    n_sym = len(qual_list)
    print(f"Quality alphabet: {qual_list} ({n_sym} symbols)\n")

    # ================================================================
    # Track convergence for different models
    # ================================================================
    models_config = {
        "Unconditional": lambda j, prev_q, prev_b, cur_b: 0,
        "Prev quality": lambda j, prev_q, prev_b, cur_b: prev_q,
        "Prev_q + cur_base": lambda j, prev_q, prev_b, cur_b: prev_q * 4 + cur_b,
        "Pos_bin5 + prev_q + bases": lambda j, prev_q, prev_b, cur_b: (j//5) * (50 * 4 * 4) + prev_q * (4*4) + prev_b * 4 + cur_b,
    }

    # Checkpoints for reporting
    checkpoints = [10, 50, 100, 500, 1000, 2000, 5000, 10000, 20000,
                   50000, 100000, 200000, 500000, 1000000]
    checkpoints = [c for c in checkpoints if c <= n_reads]

    print(f"{'Reads':>10}", end="")
    for name in models_config:
        print(f"  {name:>25}", end="")
    print(f"  {'BSC equiv':>12}")
    print("-" * (10 + 27 * len(models_config) + 14))

    # For each model: dict of context -> [counts per symbol]
    model_states = {}
    model_bits = {}
    model_bases = {}
    for name in models_config:
        model_states[name] = {}
        model_bits[name] = 0.0
        model_bases[name] = 0

    # Window tracking for recent convergence rate
    window_bits = {name: 0.0 for name in models_config}
    window_bases = {name: 0 for name in models_config}
    last_checkpoint = 0

    bsc_ref_bpb = 0.1770  # BSC bits/base from our 1M benchmark

    for i, (seq, phreds) in enumerate(records):
        prev_q = 0  # sentinel
        for j, q in enumerate(phreds):
            cur_b = BASE_MAP.get(seq[j], 0)
            prev_b = BASE_MAP.get(seq[j-1], 0) if j > 0 else 0
            sym = qual_to_idx[q]

            for name, ctx_fn in models_config.items():
                ctx = ctx_fn(j, prev_q, prev_b, cur_b)
                if ctx not in model_states[name]:
                    model_states[name][ctx] = [1] * n_sym  # Laplace init

                counts = model_states[name][ctx]
                total = sum(counts)
                bits = -math.log2(counts[sym] / total)
                model_bits[name] += bits
                model_bases[name] += 1
                window_bits[name] += bits
                window_bases[name] += 1
                counts[sym] += 1

            prev_q = q

        # Report at checkpoints
        if (i + 1) in checkpoints:
            print(f"{i+1:>10,}", end="")
            for name in models_config:
                bpb = model_bits[name] / model_bases[name]
                # Also compute recent window rate
                if window_bases[name] > 0:
                    wbpb = window_bits[name] / window_bases[name]
                else:
                    wbpb = bpb
                print(f"  {bpb:>10.4f} ({wbpb:>7.4f})", end="")
                # Reset window
                window_bits[name] = 0.0
                window_bases[name] = 0

            # BSC equivalent bytes for comparison
            total_b = model_bases[list(models_config.keys())[0]]
            bsc_equiv = total_b * bsc_ref_bpb / 8
            print(f"  {bsc_equiv:>12,.0f}")

    # ================================================================
    # Final summary
    # ================================================================
    print(f"\n=== Final Results ({n_reads:,} reads) ===\n")
    total_b = model_bases[list(models_config.keys())[0]]

    print(f"{'Model':<35} {'bits/base':>10} {'Bytes':>12} {'Contexts':>10} {'vs BSC':>10}")
    print("-" * 80)

    bsc_bytes = int(total_b * bsc_ref_bpb / 8)
    print(f"{'BSC (reference)':<35} {bsc_ref_bpb:>10.4f} {bsc_bytes:>12,} {'—':>10} {'—':>10}")

    for name in models_config:
        bpb = model_bits[name] / model_bases[name]
        b = int(model_bits[name] / 8)
        n_ctx = len(model_states[name])
        diff = b - bsc_bytes
        sign = "+" if diff >= 0 else ""
        print(f"{name:<35} {bpb:>10.4f} {b:>12,} {n_ctx:>10,} {sign}{diff:>9,}")

    # ================================================================
    # Context count analysis
    # ================================================================
    print(f"\n=== Context Statistics (best model) ===\n")
    best_name = list(models_config.keys())[-1]
    states = model_states[best_name]
    ctx_sizes = [(ctx, sum(counts) - n_sym) for ctx, counts in states.items()]  # subtract Laplace init
    ctx_sizes.sort(key=lambda x: -x[1])

    print(f"Total contexts created: {len(ctx_sizes):,}")

    # How many contexts have >= N observations?
    for threshold in [0, 1, 5, 10, 50, 100, 500, 1000]:
        n = sum(1 for _, s in ctx_sizes if s >= threshold)
        obs = sum(s for _, s in ctx_sizes if s >= threshold)
        print(f"  >= {threshold:>5} observations: {n:>6,} contexts, {obs:>12,} total obs ({100*obs/total_b:.1f}%)")

    # Median observations per context
    sizes = sorted([s for _, s in ctx_sizes], reverse=True)
    if sizes:
        median_obs = sizes[len(sizes)//2]
        print(f"\n  Median observations per context: {median_obs}")
        print(f"  Mean observations per context: {sum(sizes)/len(sizes):.1f}")
        print(f"  Min: {sizes[-1]}, Max: {sizes[0]}")

    # How many bytes needed to store the model?
    # Each context: n_sym counts (could be stored as deltas from uniform)
    model_bytes = len(ctx_sizes) * n_sym * 2  # 2 bytes per count
    print(f"\n  Model storage (raw): {model_bytes:,} bytes ({model_bytes/1024:.1f} KB)")
    print(f"  Model as fraction of compressed data: {100*model_bytes/int(model_bits[best_name]/8):.2f}%")

if __name__ == "__main__":
    main()
