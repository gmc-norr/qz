"""Checkpoint save/resume for optimizer state."""

from __future__ import annotations

import json
import logging
import os
from typing import Optional

from .population import Population

logger = logging.getLogger(__name__)


def save_checkpoint(path: str, population: Population) -> None:
    """Atomically save population state to JSON."""
    state = population.to_dict()
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(state, f, indent=2)
    os.replace(tmp, path)
    logger.debug("Checkpoint saved: gen %d -> %s", population.generation, path)


def load_checkpoint(path: str) -> Optional[Population]:
    """Load population state from checkpoint file. Returns None if not found."""
    if not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            state = json.load(f)
        pop = Population.from_dict(state)
        logger.info(
            "Resumed from checkpoint: gen %d, %d individuals, best ratio=%.4f",
            pop.generation,
            len(pop.individuals),
            pop._best_ever_ratio,
        )
        return pop
    except (json.JSONDecodeError, KeyError, TypeError) as e:
        logger.error("Failed to load checkpoint %s: %s", path, e)
        return None


def save_best_config(path: str, population: Population) -> None:
    """Save the best genome as AdvancedOptions JSON (directly usable with --config)."""
    best = population.best()
    if best is None:
        return

    config = best.to_json_config()
    with open(path, "w") as f:
        json.dump(config, f, indent=2)

    # Save metadata alongside
    meta_path = path.replace(".json", "_meta.json")
    meta = {
        "uid": best.uid,
        "generation": best.generation,
        "fitness": best.fitness,
        "ratio": best.fitness,
        "original_size": best.original_size,
        "compressed_size": best.compressed_size,
        "compress_time_s": best.compress_time_s,
        "all_parameters": {},
    }
    for name, gene in sorted(best.genes.items()):
        meta["all_parameters"][name] = gene.to_typed_value()

    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    logger.info("Best config saved: ratio=%.4f -> %s (+ %s)",
                best.fitness or 0, path, meta_path)


def save_csv_log(path: str, population: Population) -> None:
    """Write generation history as CSV."""
    header = (
        "generation,best_ratio,mean_ratio,worst_ratio,"
        "best_compressed_size,best_compress_time_s,"
        "llm_mutations,gen_time_s\n"
    )
    with open(path, "w") as f:
        f.write(header)
        for h in population.history:
            f.write(
                f"{h.generation},{h.best_ratio:.6f},{h.mean_ratio:.6f},{h.worst_ratio:.6f},"
                f"{h.best_compressed_size},{h.best_compress_time_s:.1f},"
                f"{h.llm_mutations},{h.gen_time_s:.1f}\n"
            )
