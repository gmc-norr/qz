"""Search space definition for QZ compression parameters.

Parameters are passed via --config config.json (AdvancedOptions fields)
or as CLI flags (ultra, quality-mode).
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


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


# ── Quality stream parameters ────────────────────────────────────────

QUALITY_PARAMS: List[GeneSpec] = [
    GeneSpec(
        name="quality_compressor", json_key="quality_compressor",
        dtype="categorical", low=0, high=4, default=1, group="quality",
        # 0=Zstd, 1=Bsc, 2=OpenZl, 3=Fqzcomp, 4=QualityCtx
    ),
    GeneSpec(
        name="quality_ctx_block_size", json_key="quality_ctx_block_size",
        dtype="int", low=100_000, high=2_000_000, default=500_000, group="quality",
    ),
    GeneSpec(
        name="quality_modeling", json_key="quality_modeling",
        dtype="bool", low=0, high=1, default=0, group="quality",
    ),
    GeneSpec(
        name="quality_delta", json_key="quality_delta",
        dtype="bool", low=0, high=1, default=0, group="quality",
    ),
]

# ── Sequence stream parameters ───────────────────────────────────────

SEQUENCE_PARAMS: List[GeneSpec] = [
    GeneSpec(
        name="sequence_compressor", json_key="sequence_compressor",
        dtype="categorical", low=0, high=2, default=1, group="sequence",
        # 0=Zstd, 1=Bsc, 2=OpenZl
    ),
    GeneSpec(
        name="twobit", json_key="twobit",
        dtype="bool", low=0, high=1, default=0, group="sequence",
    ),
    GeneSpec(
        name="sequence_hints", json_key="sequence_hints",
        dtype="bool", low=0, high=1, default=0, group="sequence",
    ),
    GeneSpec(
        name="sequence_delta", json_key="sequence_delta",
        dtype="bool", low=0, high=1, default=0, group="sequence",
    ),
    GeneSpec(
        name="rc_canon", json_key="rc_canon",
        dtype="bool", low=0, high=1, default=0, group="sequence",
    ),
]

# ── Header stream parameters ─────────────────────────────────────────

HEADER_PARAMS: List[GeneSpec] = [
    GeneSpec(
        name="header_compressor", json_key="header_compressor",
        dtype="categorical", low=0, high=2, default=1, group="header",
        # 0=Zstd, 1=Bsc, 2=OpenZl
    ),
    GeneSpec(
        name="header_template", json_key="header_template",
        dtype="bool", low=0, high=1, default=0, group="header",
    ),
]

# ── Shared encoding parameters ───────────────────────────────────────

ENCODING_PARAMS: List[GeneSpec] = [
    GeneSpec(
        name="bsc_static", json_key="bsc_static",
        dtype="bool", low=0, high=1, default=0, group="encoding",
    ),
    GeneSpec(
        name="dict_training", json_key="dict_training",
        dtype="bool", low=0, high=1, default=0, group="encoding",
    ),
    GeneSpec(
        name="dict_size", json_key="dict_size",
        dtype="int", low=32, high=256, default=64, group="encoding",
    ),
    GeneSpec(
        name="compression_level", json_key="compression_level",
        dtype="int", low=1, high=22, default=3, group="encoding",
    ),
    GeneSpec(
        name="bsc_block_size_mb", json_key="bsc_block_size_mb",
        dtype="int", low=5, high=100, default=25, group="encoding",
    ),
    GeneSpec(
        name="chunk_records", json_key="chunk_records",
        dtype="int", low=500_000, high=10_000_000, default=2_500_000, group="encoding",
    ),
]

# ── Ultra mode (CLI-level, not JSON) ─────────────────────────────────

ULTRA_PARAMS: List[GeneSpec] = [
    GeneSpec(
        name="ultra", json_key="",
        dtype="categorical", low=-1, high=5, default=-1, group="ultra",
        # -1=disabled, 0=auto, 1-5=level
    ),
]

# ── Quality mode (CLI-level, not JSON) ───────────────────────────────

QUALITY_MODE_PARAMS: List[GeneSpec] = [
    GeneSpec(
        name="quality_mode", json_key="",
        dtype="categorical", low=0, high=2, default=0, group="quality_mode",
        active=False,  # Lossy modes change semantics; exclude from GA by default
        # 0=Lossless, 1=IlluminaBin, 2=Discard
    ),
]

# ── Aggregate search space ───────────────────────────────────────────

SEARCH_SPACE: List[GeneSpec] = (
    QUALITY_PARAMS + SEQUENCE_PARAMS + HEADER_PARAMS +
    ENCODING_PARAMS + ULTRA_PARAMS + QUALITY_MODE_PARAMS
)

SPEC_BY_NAME: Dict[str, GeneSpec] = {s.name: s for s in SEARCH_SPACE}

GROUPS: List[str] = sorted(set(s.group for s in SEARCH_SPACE))

COMPRESSOR_NAMES = {
    "quality_compressor": {0: "Zstd", 1: "Bsc", 2: "OpenZl", 3: "Fqzcomp", 4: "QualityCtx"},
    "sequence_compressor": {0: "Zstd", 1: "Bsc", 2: "OpenZl"},
    "header_compressor": {0: "Zstd", 1: "Bsc", 2: "OpenZl"},
    "quality_mode": {0: "Lossless", 1: "IlluminaBin", 2: "Discard"},
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
