#!/usr/bin/env python3
"""
Prototype context-aware quality compression.

Builds an adaptive arithmetic-like model conditioned on:
- Position bin (position // 5)
- Previous quality value
- Previous base + current base (2-mer context)

Measures actual compressed size using arithmetic coding entropy estimate
and compares with raw BSC compression.
"""
import sys
import math
import struct
from collections import defaultdict

BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}

def parse_fastq_quals_and_seqs(path, max_reads=1000000):
    """Parse FASTQ, return list of (sequence, quality_phreds)."""
    records = []
    with open(path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            qual = f.readline().strip()
            phreds = [ord(q) - 33 for q in qual]
            records.append((seq, phreds))
            if len(records) >= max_reads:
                break
    return records

class AdaptiveModel:
    """Adaptive frequency model for a context."""
    __slots__ = ['counts', 'total']

    def __init__(self, n_symbols):
        self.counts = [1] * n_symbols  # Laplace smoothing
        self.total = n_symbols

    def prob(self, symbol):
        return self.counts[symbol] / self.total

    def update(self, symbol):
        self.counts[symbol] += 1
        self.total += 1

    def entropy_bits(self, symbol):
        """Bits needed to encode this symbol."""
        return -math.log2(self.prob(symbol))

def compress_context_model(records, pos_bin_size=5, use_prev_qual=True,
                           use_seq_context=True, use_position=True):
    """Simulate context-model compression, return total bits."""
    # Determine quality alphabet
    all_quals = set()
    for _, phreds in records:
        all_quals.update(phreds)
    qual_list = sorted(all_quals)
    qual_to_idx = {q: i for i, q in enumerate(qual_list)}
    n_symbols = len(qual_list)

    # Context -> adaptive model
    models = {}
    total_bits = 0.0
    total_symbols = 0

    for seq, phreds in records:
        prev_q = 0  # start-of-read sentinel
        for j, q in enumerate(phreds):
            # Build context key
            parts = []
            if use_position:
                parts.append(j // pos_bin_size)
            if use_prev_qual:
                parts.append(prev_q)
            if use_seq_context and j > 0:
                prev_base = BASE_MAP.get(seq[j-1], 0)
                cur_base = BASE_MAP.get(seq[j], 0)
                parts.append(prev_base * 4 + cur_base)
            else:
                parts.append(0)

            ctx = tuple(parts)
            if ctx not in models:
                models[ctx] = AdaptiveModel(n_symbols)

            model = models[ctx]
            sym = qual_to_idx[q]
            total_bits += model.entropy_bits(sym)
            model.update(sym)
            total_symbols += 1
            prev_q = q

    return total_bits, total_symbols, len(models)

def compress_prev_qual_only(records):
    """Simple model: condition on previous quality only."""
    all_quals = set()
    for _, phreds in records:
        all_quals.update(phreds)
    qual_list = sorted(all_quals)
    qual_to_idx = {q: i for i, q in enumerate(qual_list)}
    n_symbols = len(qual_list)

    models = {}
    total_bits = 0.0

    for _, phreds in records:
        prev_q = 0
        for q in phreds:
            if prev_q not in models:
                models[prev_q] = AdaptiveModel(n_symbols)
            model = models[prev_q]
            sym = qual_to_idx[q]
            total_bits += model.entropy_bits(sym)
            model.update(sym)
            prev_q = q

    return total_bits

def compress_unconditional(records):
    """Baseline: single adaptive model."""
    all_quals = set()
    for _, phreds in records:
        all_quals.update(phreds)
    qual_list = sorted(all_quals)
    qual_to_idx = {q: i for i, q in enumerate(qual_list)}
    n_symbols = len(qual_list)

    model = AdaptiveModel(n_symbols)
    total_bits = 0.0

    for _, phreds in records:
        for q in phreds:
            sym = qual_to_idx[q]
            total_bits += model.entropy_bits(sym)
            model.update(sym)

    return total_bits

def compress_position_only(records, bin_size=5):
    """Condition on position bin only."""
    all_quals = set()
    for _, phreds in records:
        all_quals.update(phreds)
    qual_list = sorted(all_quals)
    qual_to_idx = {q: i for i, q in enumerate(qual_list)}
    n_symbols = len(qual_list)

    models = {}
    total_bits = 0.0

    for _, phreds in records:
        for j, q in enumerate(phreds):
            pos_bin = j // bin_size
            if pos_bin not in models:
                models[pos_bin] = AdaptiveModel(n_symbols)
            model = models[pos_bin]
            sym = qual_to_idx[q]
            total_bits += model.entropy_bits(sym)
            model.update(sym)

    return total_bits

def main():
    if len(sys.argv) < 2:
        print("Usage: quality_context_compress.py <fastq_file> [max_reads]")
        sys.exit(1)

    path = sys.argv[1]
    max_reads = int(sys.argv[2]) if len(sys.argv) > 2 else 1000000

    print(f"Loading up to {max_reads} reads...")
    records = parse_fastq_quals_and_seqs(path, max_reads)
    n_reads = len(records)
    read_len = len(records[0][1])
    total_bases = sum(len(q) for _, q in records)
    print(f"Loaded {n_reads:,} reads, length {read_len}")
    print(f"Total quality values: {total_bases:,}")
    print()

    # Reference: BSC quality compression
    bsc_bytes = 3_319_534  # from our 1M benchmark
    if n_reads != 1_000_000:
        bsc_bytes = int(bsc_bytes * n_reads / 1_000_000)
    print(f"BSC quality compression (reference): {bsc_bytes:,} bytes ({8*bsc_bytes/total_bases:.4f} bits/base)")
    print()

    # 1. Unconditional adaptive
    print("Computing unconditional model...")
    bits_uncond = compress_unconditional(records)
    bytes_uncond = int(math.ceil(bits_uncond / 8))
    print(f"  Unconditional adaptive: {bytes_uncond:,} bytes ({bits_uncond/total_bases:.4f} bits/base)")
    print()

    # 2. Position only
    print("Computing position-only model (bin=5)...")
    bits_pos = compress_position_only(records, bin_size=5)
    bytes_pos = int(math.ceil(bits_pos / 8))
    print(f"  Position-only: {bytes_pos:,} bytes ({bits_pos/total_bases:.4f} bits/base)")
    print()

    # 3. Previous quality only
    print("Computing prev-quality-only model...")
    bits_prevq = compress_prev_qual_only(records)
    bytes_prevq = int(math.ceil(bits_prevq / 8))
    print(f"  Prev-quality-only: {bytes_prevq:,} bytes ({bits_prevq/total_bases:.4f} bits/base)")
    print()

    # 4. Full context model
    print("Computing full context model (pos_bin5 + prev_q + seq_2mer)...")
    bits_full, n_sym, n_ctx = compress_context_model(records, pos_bin_size=5,
                                                      use_prev_qual=True,
                                                      use_seq_context=True,
                                                      use_position=True)
    bytes_full = int(math.ceil(bits_full / 8))
    print(f"  Full context: {bytes_full:,} bytes ({bits_full/total_bases:.4f} bits/base)")
    print(f"  Contexts used: {n_ctx:,}")
    print()

    # 5. No position (just prev_q + seq)
    print("Computing prev_q + seq_2mer (no position)...")
    bits_nop, _, n_ctx2 = compress_context_model(records, use_prev_qual=True,
                                                  use_seq_context=True,
                                                  use_position=False)
    bytes_nop = int(math.ceil(bits_nop / 8))
    print(f"  PrevQ + seq (no pos): {bytes_nop:,} bytes ({bits_nop/total_bases:.4f} bits/base)")
    print(f"  Contexts used: {n_ctx2:,}")
    print()

    # Summary
    print("=== SUMMARY ===\n")
    print(f"{'Model':<40} {'Bytes':>12} {'bits/base':>10} {'vs BSC':>10}")
    print("-" * 75)
    results = [
        ("BSC (current)", bsc_bytes, 8*bsc_bytes/total_bases),
        ("Adaptive unconditional", bytes_uncond, bits_uncond/total_bases),
        ("Adaptive + position (bin=5)", bytes_pos, bits_pos/total_bases),
        ("Adaptive + prev quality", bytes_prevq, bits_prevq/total_bases),
        ("Adaptive + prevQ + seq_2mer", bytes_nop, bits_nop/total_bases),
        ("Adaptive + pos + prevQ + seq_2mer", bytes_full, bits_full/total_bases),
    ]
    for name, b, bpb in results:
        diff = b - bsc_bytes
        sign = "+" if diff >= 0 else ""
        print(f"{name:<40} {b:>12,} {bpb:>10.4f} {sign}{diff:>9,}")

if __name__ == "__main__":
    main()
