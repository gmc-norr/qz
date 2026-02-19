#!/usr/bin/env python3
"""Optimizer for QZ compression parameters — GA with sensitivity analysis.

Usage:
    # GA optimizer (default)
    python -m scripts.optimizer.main \
        --qz-bin target/release/qz \
        --input real_data/ERR3239334_1.1m.fastq \
        --population-size 16 --generations 30 --parallel 4

    # Sensitivity analysis
    python -m scripts.optimizer.main \
        --qz-bin target/release/qz \
        --input real_data/ERR3239334_1.1m.fastq \
        --mode sensitivity --top-params 12

    # GA with speed-aware fitness
    python -m scripts.optimizer.main ... --speed-weight 0.1
"""

from __future__ import annotations

import argparse
import json
import logging
import signal
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List

from .checkpoint import load_checkpoint, save_best_config, save_checkpoint, save_csv_log
from .config import get_active_specs
from .fitness import EvalResult, FitnessEvaluator
from .genome import Genome
from .llm_advisor import LLMAdvisor
from .population import Population

logger = logging.getLogger(__name__)

_shutdown = False


def _handle_sigint(sig, frame):
    global _shutdown
    if _shutdown:
        sys.exit(1)
    _shutdown = True
    logger.info("Shutdown requested — finishing current generation...")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Optimizer for QZ compression parameters (GA + sensitivity analysis)"
    )
    # Required paths
    p.add_argument("--qz-bin", required=True, help="Path to qz binary")
    p.add_argument("--input", required=True, help="Input FASTQ file")

    # Mode selection
    p.add_argument("--mode", choices=["ga", "sensitivity"],
                   default="ga", help="Optimization mode (default: ga)")

    # GA parameters
    p.add_argument("--population-size", type=int, default=16, help="Population size")
    p.add_argument("--generations", type=int, default=30, help="Max generations")
    p.add_argument("--elitism", type=int, default=2, help="Number of elites preserved")
    p.add_argument("--parallel", type=int, default=4,
                   help="Concurrent evaluations")
    p.add_argument("--threads-per-eval", type=int, default=4,
                   help="Threads per qz invocation")
    p.add_argument("--timeout", type=int, default=300,
                   help="Per-evaluation timeout in seconds")

    # Fitness
    p.add_argument("--verify-roundtrip", action="store_true",
                   help="Verify decompression matches original")
    p.add_argument("--speed-weight", type=float, default=0.0,
                   help="Blend ratio with speed (0=ratio only, 1=speed only)")

    # Sensitivity analysis
    p.add_argument("--top-params", type=int, default=12,
                   help="Number of top params to activate after sensitivity")
    p.add_argument("--perturbation", type=float, default=0.5,
                   help="Perturbation fraction for sensitivity analysis")

    # LLM advisor
    p.add_argument("--ollama-model", default="qwen2.5-coder:7b",
                   help="Ollama model name")
    p.add_argument("--ollama-url", default="http://localhost:11434",
                   help="Ollama base URL")
    p.add_argument("--llm-interval", type=int, default=5,
                   help="Consult LLM every N generations")
    p.add_argument("--no-llm", action="store_true", help="Disable LLM advisor")

    # Output
    p.add_argument("--checkpoint", default="optimizer_state.json",
                   help="Checkpoint file")
    p.add_argument("--output", default="best_config.json",
                   help="Best config output (usable with qz --config)")
    p.add_argument("--csv-log", default="optimizer_log.csv",
                   help="CSV log file")

    # Early stopping
    p.add_argument("--stagnation-limit", type=int, default=15,
                   help="Stop after N generations without improvement")
    p.add_argument("--ratio-target", type=float, default=10.0,
                   help="Stop when compression ratio exceeds this value")

    p.add_argument("--verbose", "-v", action="store_true", help="Debug logging")

    return p.parse_args()


def make_evaluator(args: argparse.Namespace) -> FitnessEvaluator:
    """Create a FitnessEvaluator from CLI args."""
    return FitnessEvaluator(
        qz_bin=args.qz_bin,
        input_fastq=args.input,
        threads=args.threads_per_eval,
        timeout=args.timeout,
        verify_roundtrip=args.verify_roundtrip,
        speed_weight=args.speed_weight,
    )


def evaluate_individual(args_tuple):
    """Wrapper for ProcessPoolExecutor (must be top-level picklable)."""
    (qz_bin, input_fastq, threads, timeout,
     verify_roundtrip, speed_weight, genome_dict) = args_tuple
    evaluator = FitnessEvaluator(
        qz_bin=qz_bin,
        input_fastq=input_fastq,
        threads=threads,
        timeout=timeout,
        verify_roundtrip=verify_roundtrip,
        speed_weight=speed_weight,
    )
    genome = Genome.from_dict(genome_dict)
    result = evaluator.evaluate(genome)
    return genome_dict["uid"], result


def evaluate_population_parallel(
    population: Population,
    evaluator_args: tuple,
    max_parallel: int,
) -> int:
    """Evaluate all unevaluated individuals in parallel."""
    needs_eval = population.needs_evaluation()
    if not needs_eval:
        return 0

    (qz_bin, input_fastq, threads, timeout,
     verify_roundtrip, speed_weight) = evaluator_args

    tasks = []
    for genome in needs_eval:
        tasks.append((
            qz_bin, input_fastq, threads, timeout,
            verify_roundtrip, speed_weight, genome.to_dict(),
        ))

    uid_to_genome = {g.uid: g for g in needs_eval}
    evaluated = 0

    with ProcessPoolExecutor(max_workers=max_parallel) as pool:
        futures = {pool.submit(evaluate_individual, t): t[-1]["uid"] for t in tasks}
        for future in as_completed(futures):
            uid = futures[future]
            try:
                result_uid, result = future.result()
                genome = uid_to_genome.get(result_uid)
                if genome:
                    genome.fitness = result.fitness
                    genome.original_size = result.original_size
                    genome.compressed_size = result.compressed_size
                    genome.compress_time_s = result.compress_time_s
                    evaluated += 1
                    if result.error:
                        logger.debug("Eval %s: error=%s", result_uid, result.error)
                    else:
                        logger.debug(
                            "Eval %s: ratio=%.4f (%d -> %d bytes, %.1fs)",
                            result_uid, result.ratio,
                            result.original_size, result.compressed_size,
                            result.compress_time_s,
                        )
            except Exception as e:
                logger.error("Future failed for %s: %s", uid, e)
                genome = uid_to_genome.get(uid)
                if genome:
                    genome.fitness = 0.0

    return evaluated


def run_sensitivity(args: argparse.Namespace) -> None:
    """Run sensitivity analysis mode."""
    from .sensitivity import run_sensitivity_analysis, activate_top_params, print_sensitivity_report

    evaluator = make_evaluator(args)
    results = run_sensitivity_analysis(
        evaluator,
        perturbation_fraction=args.perturbation,
    )
    print_sensitivity_report(results, top_n=args.top_params)

    # Save results
    output_path = args.output.replace(".json", "_sensitivity.json")
    with open(output_path, "w") as f:
        json.dump([{
            "name": r.name,
            "default": r.default_value,
            "perturbed": r.perturbed_value,
            "baseline_fitness": r.baseline_fitness,
            "perturbed_fitness": r.perturbed_fitness,
            "delta": r.delta,
            "abs_delta": r.abs_delta,
        } for r in results], f, indent=2)
    logger.info("Sensitivity results saved to %s", output_path)

    # Activate top params
    activated = activate_top_params(results, top_n=args.top_params)
    logger.info("Recommended active params (%d): %s", len(activated), ", ".join(activated))


def run_ga(args: argparse.Namespace) -> None:
    """Run GA optimization mode."""
    signal.signal(signal.SIGINT, _handle_sigint)

    active = get_active_specs()
    logger.info("GA mode: %d active parameters, %d total in search space",
                len(active), len(Genome.from_defaults().genes))

    # Try to resume from checkpoint
    population = load_checkpoint(args.checkpoint)
    if population is None:
        population = Population(size=args.population_size, elitism=args.elitism)
        population.initialize()
        logger.info(
            "Initialized population: %d individuals, %d params each",
            population.size, len(population.individuals[0].genes),
        )
    else:
        logger.info("Resumed from generation %d", population.generation)

    # LLM advisor
    advisor = None
    if not args.no_llm:
        advisor = LLMAdvisor(
            model=args.ollama_model,
            base_url=args.ollama_url,
            consult_interval=args.llm_interval,
        )

    evaluator_args = (
        args.qz_bin, args.input, args.threads_per_eval,
        args.timeout, args.verify_roundtrip, args.speed_weight,
    )

    start_gen = population.generation
    total_start = time.monotonic()

    for gen_idx in range(start_gen, args.generations):
        if _shutdown:
            break

        gen_start = time.monotonic()

        n_eval = evaluate_population_parallel(
            population, evaluator_args, args.parallel,
        )

        population.sort_by_fitness()
        population.update_hall_of_fame()

        # LLM advisor
        llm_mutations = 0
        stagnation = population.detect_stagnation()
        if advisor and advisor.should_consult(gen_idx, stagnation):
            suggestions = advisor.suggest_mutations(population)
            if suggestions:
                llm_mutations = advisor.apply_suggestions(
                    population.individuals, suggestions, skip_elites=args.elitism,
                )
                if llm_mutations > 0:
                    evaluate_population_parallel(
                        population, evaluator_args, args.parallel,
                    )
                    population.sort_by_fitness()

        gen_time = time.monotonic() - gen_start
        stats = population.record_stats(llm_mutations=llm_mutations, gen_time_s=gen_time)

        best = population.best()
        size_mb = best.compressed_size / 1e6 if best and best.compressed_size else 0
        logger.info(
            "Gen %03d | Best=%.4fx (%.1f MB) | Mean=%.4f | Eval=%d | LLM=%d | %.0fs",
            gen_idx, stats.best_ratio, size_mb,
            stats.mean_ratio, n_eval, llm_mutations, gen_time,
        )
        if best:
            diffs = best.diff_from_defaults()
            if diffs:
                diff_str = ", ".join(f"{k}={v[1]}" for k, v in sorted(diffs.items())[:10])
                logger.info("  Changed from defaults: %s", diff_str)

        # Checkpoint
        save_checkpoint(args.checkpoint, population)
        save_best_config(args.output, population)
        save_csv_log(args.csv_log, population)

        # Early stopping
        if stats.best_ratio >= args.ratio_target:
            logger.info("Ratio target %.4f reached — stopping.", args.ratio_target)
            break
        if stagnation >= args.stagnation_limit:
            logger.info("Stagnation limit (%d gen) reached — stopping.", args.stagnation_limit)
            break

        population.evolve()

    # Final summary
    total_time = time.monotonic() - total_start
    best = population.best()
    logger.info("=" * 60)
    logger.info("Optimization complete: %d generations in %.0fs", population.generation, total_time)
    if best:
        logger.info(
            "Best: ratio=%.4fx (%d -> %d bytes, %.1fs compress)",
            best.fitness or 0,
            best.original_size, best.compressed_size,
            best.compress_time_s,
        )
        logger.info("Best config: %s", best.param_summary())
    logger.info("Config saved to: %s", args.output)

    if population.hall_of_fame:
        logger.info("Hall of Fame:")
        for i, hof in enumerate(population.hall_of_fame[:5]):
            logger.info(
                "  #%d: ratio=%.4fx [gen %d]",
                i + 1, hof.fitness or 0, hof.generation,
            )


def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.speed_weight > 0:
        logger.info("Speed weight: %.2f (blending ratio with speed)", args.speed_weight)

    if args.mode == "sensitivity":
        run_sensitivity(args)
    else:
        run_ga(args)


if __name__ == "__main__":
    main()
