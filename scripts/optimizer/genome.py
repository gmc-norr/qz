"""Genome encoding: maps to QZ AdvancedOptions JSON + CLI flags."""

from __future__ import annotations

import json
import uuid
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .config import SEARCH_SPACE, GeneSpec


@dataclass
class Gene:
    """A single parameter value with its specification."""

    spec: GeneSpec
    value: float

    def clamp(self) -> None:
        self.value = max(self.spec.low, min(self.spec.high, self.value))

    def to_typed_value(self):
        if self.spec.dtype == "int":
            return int(round(self.value))
        if self.spec.dtype == "bool":
            return self.value >= 0.5
        if self.spec.dtype == "categorical":
            return int(round(self.value))
        return self.value


@dataclass
class Genome:
    """Full parameter set for one QZ compression invocation."""

    genes: Dict[str, Gene] = field(default_factory=dict)
    fitness: Optional[float] = None
    original_size: int = 0
    compressed_size: int = 0
    compress_time_s: float = 0.0
    generation: int = 0
    uid: str = field(default_factory=lambda: uuid.uuid4().hex[:8])
    parent_uids: List[str] = field(default_factory=list)

    @classmethod
    def from_defaults(cls) -> Genome:
        """Create a genome with all parameter values at defaults."""
        g = cls()
        for spec in SEARCH_SPACE:
            g.genes[spec.name] = Gene(spec=spec, value=spec.default)
        return g

    def copy(self) -> Genome:
        new = Genome(
            genes={k: Gene(spec=v.spec, value=v.value) for k, v in self.genes.items()},
            generation=self.generation,
            parent_uids=[self.uid],
        )
        return new

    def to_cli_args(self, input_path: str, output_path: str,
                    config_json_path: str, threads: int = 4,
                    working_dir: Optional[str] = None) -> List[str]:
        """Build command-line argument list for qz compress."""
        args = ["compress", "-i", input_path, "-o", output_path,
                "-t", str(threads), "--config", config_json_path]

        # Set working directory for temp files (critical for parallel runs —
        # ultra mode writes .qz_harc_*.tmp files that collide if multiple
        # processes share the same working directory)
        if working_dir:
            args.extend(["-w", working_dir])

        # Ultra is a CLI flag, not a JSON param
        if "ultra" in self.genes:
            ultra_val = self.genes["ultra"].to_typed_value()
            if ultra_val >= 0:
                if ultra_val == 0:
                    args.append("--ultra")  # auto level
                else:
                    args.extend(["--ultra", str(ultra_val)])

        # Quality mode is a CLI flag
        if "quality_mode" in self.genes:
            qm = self.genes["quality_mode"].to_typed_value()
            if qm == 1:
                args.extend(["--quality-mode", "illumina-bin"])
            elif qm == 2:
                args.extend(["--quality-mode", "discard"])

        return args

    def to_json_config(self) -> dict:
        """Build flat dict matching AdvancedOptions JSON structure."""
        config = {}
        for name, gene in self.genes.items():
            spec = gene.spec
            if not spec.is_json:
                continue
            config[spec.json_key] = gene.to_typed_value()
        return config

    def write_json_config(self, path: str) -> None:
        """Write JSON config file for --config flag."""
        config = self.to_json_config()
        with open(path, "w") as f:
            json.dump(config, f, indent=2)

    def to_dict(self) -> dict:
        """Serialize for checkpointing."""
        return {
            "uid": self.uid,
            "generation": self.generation,
            "parent_uids": self.parent_uids,
            "fitness": self.fitness,
            "original_size": self.original_size,
            "compressed_size": self.compressed_size,
            "compress_time_s": self.compress_time_s,
            "genes": {k: v.value for k, v in self.genes.items()},
        }

    @classmethod
    def from_dict(cls, d: dict) -> Genome:
        """Deserialize from checkpoint."""
        g = cls.from_defaults()
        g.uid = d.get("uid", uuid.uuid4().hex[:8])
        g.generation = d.get("generation", 0)
        g.parent_uids = d.get("parent_uids", [])
        g.fitness = d.get("fitness")
        g.original_size = d.get("original_size", 0)
        g.compressed_size = d.get("compressed_size", 0)
        g.compress_time_s = d.get("compress_time_s", 0.0)
        for name, val in d.get("genes", {}).items():
            if name in g.genes:
                g.genes[name].value = val
        return g

    def param_summary(self, active_only: bool = True) -> str:
        """Human-readable parameter summary."""
        from .config import COMPRESSOR_NAMES
        parts = []
        for name, gene in sorted(self.genes.items()):
            if active_only and not gene.spec.active:
                continue
            v = gene.to_typed_value()
            if name in COMPRESSOR_NAMES:
                label = COMPRESSOR_NAMES[name].get(v, str(v))
                parts.append(f"{name}={label}")
            elif gene.spec.dtype == "bool":
                parts.append(f"{name}={'on' if v else 'off'}")
            elif gene.spec.dtype == "int" or gene.spec.dtype == "categorical":
                parts.append(f"{name}={v}")
            else:
                parts.append(f"{name}={v:.3f}")
        return ", ".join(parts)

    def diff_from_defaults(self) -> Dict[str, tuple]:
        """Return params that differ from defaults."""
        diffs = {}
        for name, gene in self.genes.items():
            default = gene.spec.default
            current = gene.to_typed_value()
            if gene.spec.dtype == "float":
                if abs(current - default) > 1e-6:
                    diffs[name] = (default, current)
            elif gene.spec.dtype in ("int", "categorical"):
                if current != int(round(default)):
                    diffs[name] = (default, current)
            elif gene.spec.dtype == "bool":
                if current != (default >= 0.5):
                    diffs[name] = (default, current)
        return diffs

    def __repr__(self) -> str:
        f = f"{self.fitness:.4f}" if self.fitness is not None else "?"
        return f"Genome({self.uid}, ratio={f})"
