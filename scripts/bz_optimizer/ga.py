"""Genetic algorithm operators: selection, crossover, mutation, repair.

BZ-specific version with BAM compression constraints.
"""

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
    """Apply mutation to active genes."""
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
    """Enforce BZ parameter constraints after mutation/crossover."""
    genes = genome.genes

    # Clamp and round all values first
    for gene in genes.values():
        if gene.spec.dtype in ("int", "categorical"):
            gene.value = float(round(gene.value))
        gene.clamp()

    # quality_ctx_block_size only matters when quality_compressor == 0 (quality_ctx)
    if _cat_val(genes, "quality_compressor") != 0:
        genes["quality_ctx_block_size"].value = genes["quality_ctx_block_size"].spec.default
