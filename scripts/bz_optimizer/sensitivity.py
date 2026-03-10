"""One-at-a-time sensitivity analysis for parameter importance ranking.

For boolean params: flip the value.
For categorical params: test each non-default value.
For continuous params: perturb ±fraction of range.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import List

from .config import SEARCH_SPACE, GeneSpec, activate_params
from .fitness import FitnessEvaluator
from .genome import Genome

logger = logging.getLogger(__name__)


@dataclass
class SensitivityResult:
    """Impact of a single parameter perturbation."""

    name: str
    default_value: float
    perturbed_value: float
    baseline_fitness: float
    perturbed_fitness: float
    delta: float
    abs_delta: float


def run_sensitivity_analysis(
    evaluator: FitnessEvaluator,
    perturbation_fraction: float = 0.5,
    max_params: int = 0,
) -> List[SensitivityResult]:
    """Run one-at-a-time sensitivity analysis.

    Args:
        evaluator: Configured FitnessEvaluator.
        perturbation_fraction: Fraction of range to perturb continuous params.
        max_params: Max params to test (0 = all).

    Returns:
        Sorted list of SensitivityResult (highest impact first).
    """
    baseline = Genome.from_defaults()
    logger.info("Evaluating baseline (all defaults)...")
    baseline_result = evaluator.evaluate(baseline)
    baseline_fitness = baseline_result.fitness
    logger.info("Baseline: ratio=%.4f (%.1f MB -> %.1f MB in %.1fs)",
                baseline_fitness,
                baseline_result.original_size / 1e6,
                baseline_result.compressed_size / 1e6,
                baseline_result.compress_time_s)

    specs_to_test = list(SEARCH_SPACE)
    if max_params > 0:
        specs_to_test = specs_to_test[:max_params]

    results: List[SensitivityResult] = []
    total = len(specs_to_test)

    for i, spec in enumerate(specs_to_test):
        logger.info("[%d/%d] Testing %s (default=%s, range=[%s, %s])...",
                    i + 1, total, spec.name, spec.default, spec.low, spec.high)

        best_delta = 0.0
        best_perturbed_value = spec.default
        best_perturbed_fitness = baseline_fitness

        if spec.dtype == "bool":
            # Flip the boolean
            perturbed_value = 1.0 - spec.default
            genome = Genome.from_defaults()
            genome.genes[spec.name].value = perturbed_value
            eval_result = evaluator.evaluate(genome)
            delta = eval_result.fitness - baseline_fitness
            if abs(delta) > abs(best_delta):
                best_delta = delta
                best_perturbed_value = perturbed_value
                best_perturbed_fitness = eval_result.fitness

        elif spec.dtype == "categorical":
            # Test each non-default value
            for val in range(int(spec.low), int(spec.high) + 1):
                if val == int(spec.default):
                    continue
                genome = Genome.from_defaults()
                genome.genes[spec.name].value = float(val)
                eval_result = evaluator.evaluate(genome)
                delta = eval_result.fitness - baseline_fitness
                if abs(delta) > abs(best_delta):
                    best_delta = delta
                    best_perturbed_value = float(val)
                    best_perturbed_fitness = eval_result.fitness

        else:
            # Continuous: perturb ±fraction
            for direction in [+1, -1]:
                perturbed_value = spec.default + direction * perturbation_fraction * spec.range_size
                perturbed_value = max(spec.low, min(spec.high, perturbed_value))
                if spec.dtype == "int":
                    perturbed_value = float(round(perturbed_value))
                if abs(perturbed_value - spec.default) < 1e-10:
                    continue

                genome = Genome.from_defaults()
                genome.genes[spec.name].value = perturbed_value
                eval_result = evaluator.evaluate(genome)
                delta = eval_result.fitness - baseline_fitness
                if abs(delta) > abs(best_delta):
                    best_delta = delta
                    best_perturbed_value = perturbed_value
                    best_perturbed_fitness = eval_result.fitness

        result = SensitivityResult(
            name=spec.name,
            default_value=spec.default,
            perturbed_value=best_perturbed_value,
            baseline_fitness=baseline_fitness,
            perturbed_fitness=best_perturbed_fitness,
            delta=best_delta,
            abs_delta=abs(best_delta),
        )
        results.append(result)

        if abs(best_delta) > 0.001:
            logger.info("  -> delta=%+.4f (%.4f -> %.4f)",
                        best_delta, baseline_fitness, best_perturbed_fitness)
        else:
            logger.info("  -> no significant impact")

    results.sort(key=lambda r: r.abs_delta, reverse=True)
    return results


def activate_top_params(results: List[SensitivityResult], top_n: int = 12,
                        min_delta: float = 0.001) -> List[str]:
    """Activate the top-N most impactful params for optimization."""
    to_activate = []
    for r in results[:top_n]:
        if r.abs_delta >= min_delta:
            to_activate.append(r.name)

    count = activate_params(to_activate)
    logger.info("Activated %d params based on sensitivity analysis "
                "(top %d with delta >= %.4f)", count, top_n, min_delta)
    return to_activate


def print_sensitivity_report(results: List[SensitivityResult], top_n: int = 30) -> None:
    """Print a formatted sensitivity report."""
    from .config import COMPRESSOR_NAMES

    print(f"\n{'='*72}")
    print(f"Sensitivity Analysis Results (top {min(top_n, len(results))} of {len(results)})")
    print(f"{'='*72}")
    print(f"{'Param':<30} {'Default':>10} {'Best':>10} {'Delta':>10} {'Impact':>8}")
    print(f"{'-'*30} {'-'*10} {'-'*10} {'-'*10} {'-'*8}")

    for r in results[:top_n]:
        impact = "HIGH" if r.abs_delta > 0.05 else "MED" if r.abs_delta > 0.01 else "low"
        default_str = f"{r.default_value:.0f}" if r.default_value == int(r.default_value) else f"{r.default_value:.3f}"
        perturbed_str = f"{r.perturbed_value:.0f}" if r.perturbed_value == int(r.perturbed_value) else f"{r.perturbed_value:.3f}"
        print(f"{r.name:<30} {default_str:>10} {perturbed_str:>10} {r.delta:>+10.4f} {impact:>8}")

    print(f"{'='*72}\n")
