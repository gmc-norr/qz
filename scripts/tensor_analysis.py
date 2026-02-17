#!/usr/bin/env python3
"""
Tensor decomposition analysis for FASTQ data.

Treats reads as a 3D tensor: (reads × positions × features)
Features = 4 one-hot base channels + 1 quality channel = 5 features

Explores CP and Tucker decomposition to find latent structure
that could improve compression.
"""
import sys
import numpy as np
import tensorly as tl
from tensorly.decomposition import parafac, tucker, non_negative_parafac
import time

def parse_fastq(path, max_reads=10000):
    """Parse FASTQ file, return sequences and qualities."""
    sequences = []
    qualities = []
    with open(path) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            qual = f.readline().strip()
            sequences.append(seq)
            qualities.append(qual)
            if len(sequences) >= max_reads:
                break
    return sequences, qualities

BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}

def build_tensor(sequences, qualities, read_len=150):
    """Build (reads × positions × features) tensor.

    Features: [A, C, G, T, quality_normalized]
    """
    n = len(sequences)
    tensor = np.zeros((n, read_len, 5), dtype=np.float32)

    for i, (seq, qual) in enumerate(zip(sequences, qualities)):
        for j, (base, q) in enumerate(zip(seq, qual)):
            if j >= read_len:
                break
            # One-hot base encoding
            b = BASE_MAP.get(base, 0)
            tensor[i, j, b] = 1.0
            # Quality: normalized to [0, 1] (Phred 0-41 typical)
            tensor[i, j, 4] = (ord(q) - 33) / 41.0

    return tensor

def analyze_unfoldings(tensor):
    """SVD analysis on each mode unfolding to understand variance structure."""
    print("=== Mode Unfolding SVD Analysis ===\n")
    mode_names = ["reads", "positions", "features"]

    for mode in range(3):
        unfolded = tl.unfold(tensor, mode)
        # Compute SVD (truncated for large matrices)
        k = min(50, *unfolded.shape)
        U, S, Vt = np.linalg.svd(unfolded, full_matrices=False)
        S = S[:k]

        total_var = np.sum(S**2)
        cum_var = np.cumsum(S**2) / total_var

        print(f"Mode {mode} ({mode_names[mode]}): shape={unfolded.shape}")
        print(f"  Top singular values: {S[:10].round(1)}")
        for target in [0.5, 0.8, 0.9, 0.95, 0.99]:
            rank = np.searchsorted(cum_var, target) + 1
            print(f"  Rank for {target*100:.0f}% variance: {rank}")
        print()

def try_cp_decomposition(tensor, ranks=[5, 10, 20, 50]):
    """Try CP decomposition at various ranks, measure reconstruction error."""
    print("=== CP Decomposition ===\n")

    total_norm = np.linalg.norm(tensor)
    print(f"Tensor norm: {total_norm:.1f}")
    print(f"Tensor shape: {tensor.shape}")
    n_elements = tensor.size
    print(f"Total elements: {n_elements:,}")
    print()

    for rank in ranks:
        t0 = time.time()
        try:
            weights, factors = parafac(tensor, rank=rank, init='random',
                                       n_iter_max=50, tol=1e-4, random_state=42)
            dt = time.time() - t0

            # Reconstruction
            recon = tl.cp_to_tensor((weights, factors))
            residual = tensor - recon
            residual_norm = np.linalg.norm(residual)
            rel_error = residual_norm / total_norm
            variance_explained = 1 - rel_error**2

            # Size of factors
            factor_size = sum(f.size for f in factors) + weights.size

            # Analyze residual
            res_flat = residual.flatten()
            res_std = np.std(res_flat)
            res_max = np.max(np.abs(res_flat))

            # For bases (channels 0-3): how often is argmax preserved?
            orig_bases = np.argmax(tensor[:, :, :4], axis=2)
            recon_bases = np.argmax(recon[:, :, :4], axis=2)
            base_accuracy = np.mean(orig_bases == recon_bases)

            # Quality residual stats
            qual_residual = residual[:, :, 4]
            qual_mae = np.mean(np.abs(qual_residual))
            qual_rmse = np.sqrt(np.mean(qual_residual**2))

            print(f"Rank {rank}: error={rel_error:.4f}, var_explained={variance_explained:.4f}, "
                  f"factors={factor_size:,} floats, time={dt:.1f}s")
            print(f"  Base accuracy (argmax): {base_accuracy:.4f}")
            print(f"  Quality residual: MAE={qual_mae:.4f}, RMSE={qual_rmse:.4f} (in [0,1] scale)")
            print(f"  Residual: std={res_std:.4f}, max={res_max:.4f}")
            print(f"  Compression potential: {n_elements:,} elements -> {factor_size:,} factors "
                  f"({n_elements/factor_size:.1f}x) + residuals")
            print()
        except Exception as e:
            print(f"Rank {rank}: FAILED - {e}")
            print()

def try_tucker_decomposition(tensor, ranks_list=[(10, 10, 5), (20, 20, 5), (50, 50, 5)]):
    """Try Tucker decomposition at various ranks."""
    print("=== Tucker Decomposition ===\n")

    total_norm = np.linalg.norm(tensor)
    n_elements = tensor.size

    for ranks in ranks_list:
        t0 = time.time()
        try:
            core, factors = tucker(tensor, rank=ranks, init='random',
                                   n_iter_max=50, tol=1e-4, random_state=42)
            dt = time.time() - t0

            recon = tl.tucker_to_tensor((core, factors))
            residual = tensor - recon
            residual_norm = np.linalg.norm(residual)
            rel_error = residual_norm / total_norm
            variance_explained = 1 - rel_error**2

            core_size = core.size
            factor_size = sum(f.size for f in factors)
            total_params = core_size + factor_size

            # Base accuracy
            orig_bases = np.argmax(tensor[:, :, :4], axis=2)
            recon_bases = np.argmax(recon[:, :, :4], axis=2)
            base_accuracy = np.mean(orig_bases == recon_bases)

            # Quality residual
            qual_residual = residual[:, :, 4]
            qual_mae = np.mean(np.abs(qual_residual))

            print(f"Ranks {ranks}: error={rel_error:.4f}, var_explained={variance_explained:.4f}, "
                  f"time={dt:.1f}s")
            print(f"  Core: {core.shape} ({core_size:,}), Factors: {factor_size:,}, "
                  f"Total params: {total_params:,}")
            print(f"  Base accuracy: {base_accuracy:.4f}, Quality MAE: {qual_mae:.4f}")
            print(f"  Compression: {n_elements:,} -> {total_params:,} params "
                  f"({n_elements/total_params:.1f}x) + residuals")
            print()
        except Exception as e:
            print(f"Ranks {ranks}: FAILED - {e}")
            print()

def analyze_separate_streams(tensor):
    """Analyze base and quality tensors separately (as currently done)."""
    print("=== Separate Stream Analysis ===\n")

    # Base tensor: reads × positions, integer-valued (0-3)
    bases = np.argmax(tensor[:, :, :4], axis=2)  # (reads, positions)
    quals = tensor[:, :, 4]  # (reads, positions), float in [0,1]

    n_reads, read_len = bases.shape

    # Base statistics per position
    print(f"Base tensor: {bases.shape}, dtype={bases.dtype}")
    base_entropy_per_pos = []
    for j in range(read_len):
        counts = np.bincount(bases[:, j], minlength=4)
        probs = counts / counts.sum()
        probs = probs[probs > 0]
        h = -np.sum(probs * np.log2(probs))
        base_entropy_per_pos.append(h)
    mean_entropy = np.mean(base_entropy_per_pos)
    print(f"  Mean per-position base entropy: {mean_entropy:.4f} bits")
    print(f"  Theoretical min for bases: {n_reads * read_len * mean_entropy / 8:,.0f} bytes")
    print(f"  At 2 bits/base (raw): {n_reads * read_len * 2 / 8:,.0f} bytes")
    print()

    # Quality SVD
    print(f"Quality tensor: {quals.shape}")
    k = min(50, *quals.shape)
    U, S, Vt = np.linalg.svd(quals, full_matrices=False)
    S = S[:k]
    total_var = np.sum(S**2)
    cum_var = np.cumsum(S**2) / total_var
    print(f"  Quality SVD - top singular values: {S[:10].round(2)}")
    for target in [0.5, 0.8, 0.9, 0.95, 0.99]:
        rank = np.searchsorted(cum_var, target) + 1
        print(f"  Rank for {target*100:.0f}% variance: {rank}")

    # Position-wise quality mean and std
    qual_means = np.mean(quals, axis=0)
    qual_stds = np.std(quals, axis=0)
    print(f"\n  Quality mean across positions: min={qual_means.min():.3f}, "
          f"max={qual_means.max():.3f}, overall={np.mean(quals):.3f}")
    print(f"  Quality std across positions: min={qual_stds.min():.3f}, "
          f"max={qual_stds.max():.3f}, overall={np.std(quals):.3f}")

    # After subtracting position means, how much variance remains?
    quals_centered = quals - qual_means[np.newaxis, :]
    var_from_position = 1 - np.var(quals_centered) / np.var(quals)
    print(f"  Variance explained by position means alone: {var_from_position:.4f}")
    print()

def analyze_residual_compressibility(tensor, rank=20):
    """Check if tensor residuals are more compressible than raw data."""
    print(f"=== Residual Compressibility (CP rank={rank}) ===\n")

    weights, factors = parafac(tensor, rank=rank, init='random',
                               n_iter_max=50, tol=1e-4, random_state=42)
    recon = tl.cp_to_tensor((weights, factors))
    residual = tensor - recon

    # Convert to bytes and check entropy
    def byte_entropy(data):
        flat = data.flatten()
        # Quantize to uint8
        q = np.clip(np.round(flat * 255), 0, 255).astype(np.uint8)
        counts = np.bincount(q, minlength=256)
        probs = counts / counts.sum()
        probs = probs[probs > 0]
        return -np.sum(probs * np.log2(probs))

    orig_entropy = byte_entropy(tensor)
    res_entropy = byte_entropy(residual + 0.5)  # shift to [0,1] range

    print(f"Original tensor entropy: {orig_entropy:.2f} bits/byte")
    print(f"Residual entropy: {res_entropy:.2f} bits/byte")
    print(f"Entropy reduction: {orig_entropy - res_entropy:.2f} bits/byte")
    print()

    # More practically: look at base residuals and quality residuals separately
    base_res = residual[:, :, :4]
    qual_res = residual[:, :, 4]

    # Quality residual: quantize to int8 range
    qual_res_q = np.clip(np.round(qual_res * 41), -128, 127).astype(np.int8)
    unique_vals, counts = np.unique(qual_res_q, return_counts=True)
    probs = counts / counts.sum()
    qual_res_entropy = -np.sum(probs * np.log2(probs))
    print(f"Quality residual (int8, Phred scale): {len(unique_vals)} unique values")
    print(f"  Entropy: {qual_res_entropy:.2f} bits/value")
    print(f"  vs raw quality entropy: ", end="")

    raw_quals = np.clip(np.round(tensor[:, :, 4] * 41), 0, 41).astype(np.uint8)
    uv, c = np.unique(raw_quals, return_counts=True)
    p = c / c.sum()
    raw_q_ent = -np.sum(p * np.log2(p))
    print(f"{raw_q_ent:.2f} bits/value")
    print()

def main():
    if len(sys.argv) < 2:
        print("Usage: tensor_analysis.py <fastq_file> [max_reads]")
        sys.exit(1)

    path = sys.argv[1]
    max_reads = int(sys.argv[2]) if len(sys.argv) > 2 else 10000

    print(f"Loading {max_reads} reads from {path}...")
    sequences, qualities = parse_fastq(path, max_reads)

    # Filter to fixed-length reads
    read_len = 150
    filtered = [(s, q) for s, q in zip(sequences, qualities)
                if len(s) == read_len and len(q) == read_len]
    print(f"Loaded {len(filtered)} reads of length {read_len} (from {len(sequences)} total)")
    sequences, qualities = zip(*filtered)

    print(f"Building tensor ({len(sequences)} × {read_len} × 5)...")
    tensor = build_tensor(sequences, qualities, read_len)
    print(f"Tensor shape: {tensor.shape}, memory: {tensor.nbytes / 1024**2:.1f} MB\n")

    # 1. Analyze structure per stream
    analyze_separate_streams(tensor)

    # 2. SVD on mode unfoldings
    analyze_unfoldings(tensor)

    # 3. CP decomposition
    try_cp_decomposition(tensor, ranks=[5, 10, 20, 50])

    # 4. Tucker decomposition
    try_tucker_decomposition(tensor, ranks_list=[
        (10, 10, 5),
        (20, 20, 5),
        (50, 50, 5),
        (100, 50, 5),
    ])

    # 5. Residual compressibility
    analyze_residual_compressibility(tensor, rank=20)

if __name__ == "__main__":
    main()
