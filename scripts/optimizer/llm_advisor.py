"""LLM-guided mutation advisor using Ollama for compression optimization."""

from __future__ import annotations

import json
import logging
import random
import re
from typing import Dict, List, Optional

from .config import SPEC_BY_NAME, COMPRESSOR_NAMES, get_active_specs
from .genome import Genome
from .population import Population

logger = logging.getLogger(__name__)


class LLMAdvisor:
    """Queries a local Ollama LLM for intelligent mutation suggestions."""

    def __init__(
        self,
        model: str = "qwen2.5-coder:7b",
        base_url: str = "http://localhost:11434",
        consult_interval: int = 5,
        stagnation_trigger: int = 5,
    ):
        self.model = model
        self.base_url = base_url.rstrip("/")
        self.consult_interval = consult_interval
        self.stagnation_trigger = stagnation_trigger
        self._available: Optional[bool] = None

    def is_available(self) -> bool:
        if self._available is not None:
            return self._available
        try:
            import urllib.request
            req = urllib.request.Request(f"{self.base_url}/api/tags", method="GET")
            with urllib.request.urlopen(req, timeout=5) as resp:
                self._available = resp.status == 200
        except Exception:
            self._available = False
        if not self._available:
            logger.warning("Ollama not available at %s — LLM advisor disabled", self.base_url)
        return self._available

    def should_consult(self, generation: int, stagnation: int) -> bool:
        if not self.is_available():
            return False
        if generation > 0 and generation % self.consult_interval == 0:
            return True
        if stagnation >= self.stagnation_trigger:
            return True
        return False

    def suggest_mutations(self, population: Population) -> List[Dict]:
        if not self.is_available():
            return []

        prompt = self._build_prompt(population)
        response = self._query_ollama(prompt)
        if response is None:
            return []

        suggestions = self._parse_suggestions(response)
        if suggestions:
            logger.info(
                "LLM suggested %d mutations: %s",
                len(suggestions),
                ", ".join(s["param"] for s in suggestions),
            )
        return suggestions

    def apply_suggestions(
        self, individuals: List[Genome], suggestions: List[Dict], skip_elites: int = 2
    ) -> int:
        if not suggestions or len(individuals) <= skip_elites:
            return 0

        n_applied = 0
        for suggestion in suggestions:
            param = suggestion["param"]
            value = suggestion["value"]

            candidates = list(range(skip_elites, len(individuals)))
            targets = random.sample(candidates, min(2, len(candidates)))

            for idx in targets:
                genome = individuals[idx]
                if param in genome.genes:
                    genome.genes[param].value = value
                    genome.genes[param].clamp()
                    genome.fitness = None
                    n_applied += 1

        return n_applied

    def _build_prompt(self, pop: Population) -> str:
        best = pop.best()
        if best is None:
            return ""

        history_snippet = ""
        if len(pop.history) >= 2:
            recent = pop.history[-10:]
            history_snippet = ", ".join(f"{h.best_ratio:.4f}" for h in recent)

        stagnation = pop._stagnation_count

        active_specs = get_active_specs()
        param_lines = []
        for spec in active_specs:
            if spec.name in best.genes:
                gene = best.genes[spec.name]
                v = gene.to_typed_value()
                if spec.name in COMPRESSOR_NAMES:
                    label = COMPRESSOR_NAMES[spec.name].get(v, str(v))
                    param_lines.append(
                        f"  {spec.name}: {v} ({label}) (range {spec.low}-{spec.high}, default {spec.default})"
                    )
                elif spec.dtype == "bool":
                    param_lines.append(
                        f"  {spec.name}: {'on' if v else 'off'} (range 0-1, default {spec.default})"
                    )
                else:
                    param_lines.append(
                        f"  {spec.name}: {v} (range {spec.low}-{spec.high}, default {spec.default})"
                    )
        params_str = "\n".join(param_lines)

        size_mb = best.compressed_size / 1e6 if best.compressed_size else 0

        return f"""You are an expert in FASTQ genomic data compression.

I am optimizing a FASTQ compressor's ~{len(active_specs)} parameters using an evolutionary algorithm.
The compressor (QZ) uses columnar stream separation (headers, sequences, quality scores) with
BSC/BWT block-sorting compression. Quality scores can use a context-adaptive range coder.

Current best performance:
- Compression ratio: {best.fitness:.4f}x ({size_mb:.1f} MB compressed)
- Compress time: {best.compress_time_s:.1f}s
- Population mean ratio: {pop.history[-1].mean_ratio:.4f if pop.history else 0:.4f}
- Stagnation: {stagnation} generations without improvement
- Ratio trend (last 10 gen): [{history_snippet}]

Current best parameter values:
{params_str}

Key compressor options:
- quality_compressor: 0=Zstd, 1=Bsc(BWT), 2=OpenZl, 3=Fqzcomp, 4=QualityCtx(range coder)
- sequence_compressor: 0=Zstd, 1=Bsc(BWT), 2=OpenZl
- header_compressor: 0=Zstd, 1=Bsc(BWT), 2=OpenZl

Based on the current performance and parameter values, suggest 3-5 specific parameter
mutations that could improve the compression ratio. For each suggestion:
1. Use the exact parameter name from the list above
2. Suggest a concrete numeric value within the valid range
3. Give a brief reason

Respond ONLY as a JSON array, no other text:
[{{"param": "param_name", "value": 123, "reason": "brief explanation"}}]"""

    def _query_ollama(self, prompt: str) -> Optional[str]:
        try:
            import urllib.request

            payload = json.dumps({
                "model": self.model,
                "prompt": prompt,
                "stream": False,
                "options": {"temperature": 0.7, "num_predict": 512},
            }).encode()

            req = urllib.request.Request(
                f"{self.base_url}/api/generate",
                data=payload,
                headers={"Content-Type": "application/json"},
                method="POST",
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = json.loads(resp.read())
                return data.get("response", "")

        except Exception as e:
            logger.warning("Ollama query failed: %s", e)
            self._available = None
            return None

    def _parse_suggestions(self, response: str) -> List[Dict]:
        try:
            json_match = re.search(r"\[.*\]", response, re.DOTALL)
            if not json_match:
                logger.debug("No JSON array found in LLM response")
                return []

            suggestions = json.loads(json_match.group())
            valid = []
            for s in suggestions:
                param = s.get("param", "")
                if param not in SPEC_BY_NAME:
                    logger.debug("LLM suggested unknown param: %s", param)
                    continue

                spec = SPEC_BY_NAME[param]
                try:
                    value = float(s["value"])
                except (ValueError, TypeError):
                    continue

                if not (spec.low <= value <= spec.high):
                    logger.debug(
                        "LLM value %.2f out of range [%.2f, %.2f] for %s",
                        value, spec.low, spec.high, param,
                    )
                    continue

                valid.append({"param": param, "value": value, "reason": s.get("reason", "")})

            return valid

        except (json.JSONDecodeError, KeyError, TypeError) as e:
            logger.debug("Failed to parse LLM suggestions: %s", e)
            return []
