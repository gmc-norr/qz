"""Population management: initialization, evolution, hall of fame."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List, Optional

from .config import get_active_specs
from .ga import grouped_crossover, mutate_gaussian, repair_genome, tournament_select
from .genome import Genome

logger = logging.getLogger(__name__)


@dataclass
class GenerationStats:
    """Statistics for one generation."""

    generation: int
    best_ratio: float
    mean_ratio: float
    worst_ratio: float
    best_compressed_size: int
    best_compress_time_s: float
    llm_mutations: int = 0
    gen_time_s: float = 0.0

    def to_dict(self) -> dict:
        return {
            "generation": self.generation,
            "best_ratio": self.best_ratio,
            "mean_ratio": self.mean_ratio,
            "worst_ratio": self.worst_ratio,
            "best_compressed_size": self.best_compressed_size,
            "best_compress_time_s": self.best_compress_time_s,
            "llm_mutations": self.llm_mutations,
            "gen_time_s": self.gen_time_s,
        }

    @classmethod
    def from_dict(cls, d: dict) -> GenerationStats:
        return cls(**d)


class Population:
    """Manages a population of genomes through evolution."""

    def __init__(self, size: int = 16, elitism: int = 2):
        self.size = size
        self.elitism = elitism
        self.individuals: List[Genome] = []
        self.hall_of_fame: List[Genome] = []
        self.generation: int = 0
        self.history: List[GenerationStats] = []
        self._stagnation_count: int = 0
        self._best_ever_ratio: float = 0.0

    def initialize(self, seed: Optional[Genome] = None) -> None:
        """Create initial population from seed (defaults) with mutations."""
        if seed is None:
            seed = Genome.from_defaults()

        self.individuals = [seed]  # Slot 0: exact defaults (baseline)

        for _ in range(1, self.size):
            mutant = seed.copy()
            mutate_gaussian(mutant, mutation_rate=0.4, sigma=0.25)
            repair_genome(mutant)
            self.individuals.append(mutant)

    def sort_by_fitness(self) -> None:
        """Sort population by fitness descending."""
        self.individuals.sort(key=lambda g: g.fitness or 0.0, reverse=True)

    def update_hall_of_fame(self, max_size: int = 10) -> None:
        """Update hall of fame with current best individuals."""
        for ind in self.individuals[:5]:
            if ind.fitness is not None and ind.fitness > 0:
                if not any(h.uid == ind.uid for h in self.hall_of_fame):
                    self.hall_of_fame.append(Genome.from_dict(ind.to_dict()))
        self.hall_of_fame.sort(key=lambda g: g.fitness or 0.0, reverse=True)
        self.hall_of_fame = self.hall_of_fame[:max_size]

    def detect_stagnation(self) -> int:
        """How many generations since last improvement."""
        best = max((g.fitness or 0.0) for g in self.individuals)
        if best > self._best_ever_ratio + 1e-6:
            self._best_ever_ratio = best
            self._stagnation_count = 0
        else:
            self._stagnation_count += 1
        return self._stagnation_count

    def record_stats(self, llm_mutations: int = 0, gen_time_s: float = 0.0) -> GenerationStats:
        """Record statistics for the current generation."""
        fitnesses = [g.fitness or 0.0 for g in self.individuals]
        best_ind = self.individuals[0]  # Assumed sorted
        stats = GenerationStats(
            generation=self.generation,
            best_ratio=max(fitnesses),
            mean_ratio=sum(fitnesses) / len(fitnesses),
            worst_ratio=min(fitnesses),
            best_compressed_size=best_ind.compressed_size,
            best_compress_time_s=best_ind.compress_time_s,
            llm_mutations=llm_mutations,
            gen_time_s=gen_time_s,
        )
        self.history.append(stats)
        return stats

    def evolve(self) -> None:
        """Produce next generation using selection, crossover, mutation."""
        self.sort_by_fitness()
        stagnation = self._stagnation_count

        # Adaptive rates: increase exploration on stagnation
        mutation_rate = min(0.5, 0.15 + 0.05 * stagnation)
        sigma = min(0.35, 0.1 + 0.03 * stagnation)

        if stagnation > 0 and stagnation % 5 == 0:
            logger.info(
                "Stagnation at %d generations (mutation_rate=%.2f, sigma=%.2f)",
                stagnation, mutation_rate, sigma,
            )

        # Elitism: keep top N unchanged
        next_gen = []
        for ind in self.individuals[: self.elitism]:
            elite = Genome.from_dict(ind.to_dict())
            elite.fitness = ind.fitness
            elite.original_size = ind.original_size
            elite.compressed_size = ind.compressed_size
            elite.compress_time_s = ind.compress_time_s
            next_gen.append(elite)

        # Fill rest with offspring
        while len(next_gen) < self.size:
            p1 = tournament_select(self.individuals, k=3)
            p2 = tournament_select(self.individuals, k=3)
            c1, c2 = grouped_crossover(p1, p2)
            mutate_gaussian(c1, mutation_rate, sigma)
            mutate_gaussian(c2, mutation_rate, sigma)
            repair_genome(c1)
            repair_genome(c2)
            c1.generation = self.generation + 1
            c2.generation = self.generation + 1
            next_gen.extend([c1, c2])

        self.individuals = next_gen[: self.size]
        self.generation += 1

    def needs_evaluation(self) -> List[Genome]:
        """Return individuals that need fitness evaluation."""
        return [g for g in self.individuals if g.fitness is None]

    def best(self) -> Optional[Genome]:
        """Return the best individual (by fitness)."""
        evaluated = [g for g in self.individuals if g.fitness is not None]
        if not evaluated:
            return None
        return max(evaluated, key=lambda g: g.fitness)

    def to_dict(self) -> dict:
        return {
            "generation": self.generation,
            "size": self.size,
            "elitism": self.elitism,
            "stagnation_count": self._stagnation_count,
            "best_ever_ratio": self._best_ever_ratio,
            "individuals": [g.to_dict() for g in self.individuals],
            "hall_of_fame": [g.to_dict() for g in self.hall_of_fame],
            "history": [h.to_dict() for h in self.history],
        }

    @classmethod
    def from_dict(cls, d: dict) -> Population:
        pop = cls(size=d["size"], elitism=d["elitism"])
        pop.generation = d["generation"]
        pop._stagnation_count = d.get("stagnation_count", 0)
        pop._best_ever_ratio = d.get("best_ever_ratio", 0.0)
        pop.individuals = [Genome.from_dict(g) for g in d["individuals"]]
        pop.hall_of_fame = [Genome.from_dict(g) for g in d.get("hall_of_fame", [])]
        pop.history = [GenerationStats.from_dict(h) for h in d.get("history", [])]
        return pop
