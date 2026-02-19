"""Genetic algorithm operators: selection, crossover, mutation, repair."""

from __future__ import annotations

import random
from typing import List, Tuple

from .config import get_active_specs
from .genome import Gene, Genome


def tournament_select(population: List[Genome], k: int = 3) -> Genome:
    """Select one individual via tournament selection."""
    contestants = random.sample(population, min(k, len(population)))
    return max(contestants, key=lambda g: g.fitness or 0.0)


def grouped_crossover(
    parent1: Genome, parent2: Genome, n_swap: int = 2
) -> Tuple[Genome, Genome]:
    """Crossover by swapping entire parameter groups between parents."""
    child1 = parent1.copy()
    child2 = parent2.copy()
    child1.parent_uids = [parent1.uid, parent2.uid]
    child2.parent_uids = [parent1.uid, parent2.uid]

    active_groups = sorted(set(s.group for s in get_active_specs()))
    if not active_groups:
        return child1, child2

    swap_groups = random.sample(active_groups, min(n_swap, len(active_groups)))

    for name in child1.genes:
        gene = child1.genes[name]
        if gene.spec.active and gene.spec.group in swap_groups:
            child1.genes[name].value, child2.genes[name].value = (
                child2.genes[name].value,
                child1.genes[name].value,
            )

    child1.fitness = None
    child2.fitness = None
    return child1, child2


def mutate_gaussian(
    genome: Genome, mutation_rate: float = 0.1, sigma: float = 0.1
) -> None:
    """Apply mutation to active genes.

    Bool: flip. Categorical: random choice. Continuous: Gaussian noise.
    """
    for gene in genome.genes.values():
        if not gene.spec.active:
            continue
        if random.random() > mutation_rate:
            continue

        if gene.spec.dtype == "bool":
            gene.value = 1.0 - gene.value
        elif gene.spec.dtype == "categorical":
            gene.value = float(random.randint(int(gene.spec.low), int(gene.spec.high)))
        elif gene.spec.log_scale and gene.value > 0:
            import math
            log_val = math.log(gene.value + 1e-10)
            log_val += random.gauss(0, sigma)
            gene.value = math.exp(log_val)
        else:
            perturbation = random.gauss(0, sigma * gene.spec.range_size)
            gene.value += perturbation

        gene.clamp()

    genome.fitness = None


def _is_on(genes, name):
    return name in genes and genes[name].to_typed_value()


def _cat_val(genes, name):
    return int(round(genes[name].value)) if name in genes else -1


def repair_genome(genome: Genome) -> None:
    """Enforce QZ parameter constraints after mutation/crossover."""
    genes = genome.genes

    # Clamp and round all values first
    for gene in genes.values():
        if gene.spec.dtype in ("int", "categorical"):
            gene.value = float(round(gene.value))
        gene.clamp()

    # Ultra mode overrides everything: force BSC compressors, disable extras
    if "ultra" in genes and _cat_val(genes, "ultra") >= 0:
        genes["sequence_compressor"].value = 1.0  # BSC
        genes["header_compressor"].value = 1.0  # BSC
        genes["sequence_hints"].value = 0.0
        genes["sequence_delta"].value = 0.0
        genes["rc_canon"].value = 0.0
        genes["twobit"].value = 0.0
        genes["header_template"].value = 0.0
        genes["quality_modeling"].value = 0.0
        genes["quality_delta"].value = 0.0
        genes["dict_training"].value = 0.0
        # Ultra uses its own block/chunk sizes
        if "bsc_block_size_mb" in genes:
            genes["bsc_block_size_mb"].value = genes["bsc_block_size_mb"].spec.default
        if "chunk_records" in genes:
            genes["chunk_records"].value = genes["chunk_records"].spec.default
        return

    # sequence_hints requires BSC sequence compressor
    if _is_on(genes, "sequence_hints") and _cat_val(genes, "sequence_compressor") != 1:
        genes["sequence_compressor"].value = 1.0

    # sequence_delta requires BSC and is incompatible with sequence_hints
    if _is_on(genes, "sequence_delta"):
        if _cat_val(genes, "sequence_compressor") != 1:
            genes["sequence_compressor"].value = 1.0
        if _is_on(genes, "sequence_hints"):
            if random.random() < 0.5:
                genes["sequence_hints"].value = 0.0
            else:
                genes["sequence_delta"].value = 0.0

    # sequence_delta and rc_canon are incompatible
    if _is_on(genes, "sequence_delta") and _is_on(genes, "rc_canon"):
        if random.random() < 0.5:
            genes["rc_canon"].value = 0.0
        else:
            genes["sequence_delta"].value = 0.0

    # Quality mode Discard makes quality params irrelevant
    if _cat_val(genes, "quality_mode") == 2:
        genes["quality_compressor"].value = 1.0
        genes["quality_modeling"].value = 0.0
        genes["quality_delta"].value = 0.0
        genes["dict_training"].value = 0.0

    # QualityCtx requires Lossless quality mode
    if _cat_val(genes, "quality_compressor") == 4 and _cat_val(genes, "quality_mode") != 0:
        genes["quality_mode"].value = 0.0

    # quality_modeling and quality_delta are mutually exclusive
    if _is_on(genes, "quality_modeling") and _is_on(genes, "quality_delta"):
        if random.random() < 0.5:
            genes["quality_modeling"].value = 0.0
        else:
            genes["quality_delta"].value = 0.0

    # dict_size only matters if dict_training is on
    if not _is_on(genes, "dict_training"):
        genes["dict_size"].value = genes["dict_size"].spec.default

    # compression_level only matters for Zstd quality compressor
    if _cat_val(genes, "quality_compressor") != 0:
        genes["compression_level"].value = genes["compression_level"].spec.default

    # quality_ctx_block_size only matters for QualityCtx compressor
    if _cat_val(genes, "quality_compressor") != 4:
        genes["quality_ctx_block_size"].value = genes["quality_ctx_block_size"].spec.default
