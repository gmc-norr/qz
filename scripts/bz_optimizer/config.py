"""Search space definition for BZ BAM compression parameters.

Parameters are passed via --config config.json (AdvancedOptions fields).
"""

from dataclasses import dataclass
from typing import Dict, List


@dataclass
class GeneSpec:
    """Specification for a single optimizable parameter."""

    name: str
    dtype: str  # "int", "float", "bool", or "categorical"
    low: float
    high: float
    default: float
    group: str
    json_key: str = ""  # Key in AdvancedOptions JSON config
    log_scale: bool = False
    active: bool = True

    @property
    def range_size(self) -> float:
        return self.high - self.low

    @property
    def is_json(self) -> bool:
        return bool(self.json_key)


# ── Chunking parameters ────────────────────────────────────────────

CHUNKING_PARAMS: List[GeneSpec] = [
    GeneSpec(
        name="chunk_size", json_key="chunk_size",
        dtype="int", low=500_000, high=10_000_000, default=2_500_000, group="chunking",
    ),
    GeneSpec(
        name="bsc_block_size_mb", json_key="bsc_block_size_mb",
        dtype="int", low=5, high=100, default=25, group="chunking",
    ),
]

# ── Quality parameters ──────────────────────────────────────────────

QUALITY_PARAMS: List[GeneSpec] = [
    GeneSpec(
        name="quality_compressor", json_key="quality_compressor",
        dtype="categorical", low=0, high=1, default=0, group="quality",
        # 0=quality_ctx (context-adaptive), 1=bsc (BWT)
    ),
    GeneSpec(
        name="quality_ctx_block_size", json_key="quality_ctx_block_size",
        dtype="int", low=100_000, high=2_000_000, default=500_000, group="quality",
    ),
]

# ── BSC encoding parameters ─────────────────────────────────────────

ENCODING_PARAMS: List[GeneSpec] = [
    GeneSpec(
        name="use_lzp", json_key="use_lzp",
        dtype="bool", low=0, high=1, default=0, group="encoding",
    ),
    GeneSpec(
        name="bsc_adaptive", json_key="bsc_adaptive",
        dtype="bool", low=0, high=1, default=1, group="encoding",
    ),
]

# ── Per-stream compressor routing ────────────────────────────────────

STREAM_PARAMS: List[GeneSpec] = [
    GeneSpec(
        name="alignment_compressor", json_key="alignment_compressor",
        dtype="categorical", low=0, high=1, default=0, group="stream_routing",
        # 0=bsc, 1=zstd
    ),
    GeneSpec(
        name="aux_compressor", json_key="aux_compressor",
        dtype="categorical", low=0, high=1, default=0, group="stream_routing",
        # 0=bsc, 1=zstd
    ),
]

# ── Aggregate search space ───────────────────────────────────────────

SEARCH_SPACE: List[GeneSpec] = (
    CHUNKING_PARAMS + QUALITY_PARAMS + ENCODING_PARAMS + STREAM_PARAMS
)

SPEC_BY_NAME: Dict[str, GeneSpec] = {s.name: s for s in SEARCH_SPACE}

GROUPS: List[str] = sorted(set(s.group for s in SEARCH_SPACE))

COMPRESSOR_NAMES = {
    "quality_compressor": {0: "QualityCtx", 1: "Bsc"},
    "alignment_compressor": {0: "Bsc", 1: "Zstd"},
    "aux_compressor": {0: "Bsc", 1: "Zstd"},
}


def get_active_specs() -> List[GeneSpec]:
    """Return only active parameter specs."""
    return [s for s in SEARCH_SPACE if s.active]


def activate_params(names: List[str]) -> int:
    """Activate specific params by name. Returns count activated."""
    count = 0
    for s in SEARCH_SPACE:
        if s.name in names and not s.active:
            s.active = True
            count += 1
    return count


def deactivate_params(names: List[str]) -> int:
    """Deactivate specific params by name. Returns count deactivated."""
    count = 0
    for s in SEARCH_SPACE:
        if s.name in names and s.active:
            s.active = False
            count += 1
    return count
