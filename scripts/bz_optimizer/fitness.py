"""Fitness evaluation: run bz compress, measure compressed size."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from typing import Optional

from .genome import Genome

logger = logging.getLogger(__name__)


@dataclass
class EvalResult:
    """Result of a single fitness evaluation."""

    fitness: float  # compression ratio (original / compressed)
    original_size: int
    compressed_size: int
    ratio: float
    compress_time_s: float
    decompress_time_s: float = 0.0
    roundtrip_ok: bool = True
    error: Optional[str] = None


class FitnessEvaluator:
    """Evaluates genomes by running bz compress and measuring output size."""

    def __init__(
        self,
        bz_bin: str,
        input_bam: str,
        threads: int = 4,
        timeout: int = 600,
        verify_roundtrip: bool = False,
        speed_weight: float = 0.0,
        baseline_time: Optional[float] = None,
    ):
        self.bz_bin = bz_bin
        self.input_bam = input_bam
        self.threads = threads
        self.timeout = timeout
        self.verify_roundtrip = verify_roundtrip
        self.speed_weight = speed_weight
        self.baseline_time = baseline_time
        self._original_size: Optional[int] = None

    @property
    def original_size(self) -> int:
        if self._original_size is None:
            self._original_size = os.path.getsize(self.input_bam)
        return self._original_size

    def evaluate(self, genome: Genome, work_dir: Optional[str] = None) -> EvalResult:
        """Evaluate a single genome by compressing and measuring size."""
        cleanup = work_dir is None
        if work_dir is None:
            work_dir = tempfile.mkdtemp(prefix=f"bz_opt_{genome.uid}_")

        output_path = os.path.join(work_dir, "output.bz")
        config_path = os.path.join(work_dir, "config.json")

        try:
            genome.write_json_config(config_path)

            cli_args = genome.to_cli_args(
                input_path=self.input_bam,
                output_path=output_path,
                config_json_path=config_path,
                threads=self.threads,
                working_dir=work_dir,
            )
            cmd = [self.bz_bin] + cli_args
            env = {**os.environ}

            t0 = time.monotonic()
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                env=env,
            )
            compress_time = time.monotonic() - t0

            if result.returncode != 0:
                stderr = result.stderr.strip()
                logger.warning(
                    "bz failed for %s (rc=%d):\n  CMD: %s\n  STDERR (last 1500 chars): %s",
                    genome.uid, result.returncode,
                    " ".join(cmd),
                    stderr[-1500:],
                )
                return self._error_result(compress_time, f"bz exit {result.returncode}")

            if not os.path.exists(output_path):
                return self._error_result(compress_time, "no output file produced")

            compressed_size = os.path.getsize(output_path)
            if compressed_size == 0:
                return self._error_result(compress_time, "empty output file")

            ratio = self.original_size / compressed_size

            # Optional roundtrip verification
            decompress_time = 0.0
            roundtrip_ok = True
            if self.verify_roundtrip:
                roundtrip_ok, decompress_time = self._verify(output_path, work_dir)
                if not roundtrip_ok:
                    return self._error_result(
                        compress_time, "roundtrip verification failed"
                    )

            # Compute fitness
            fitness = ratio
            if self.speed_weight > 0 and self.baseline_time and self.baseline_time > 0:
                speed_score = self.baseline_time / max(compress_time, 0.1)
                fitness = ratio * (1 - self.speed_weight) + speed_score * self.speed_weight

            return EvalResult(
                fitness=fitness,
                original_size=self.original_size,
                compressed_size=compressed_size,
                ratio=ratio,
                compress_time_s=compress_time,
                decompress_time_s=decompress_time,
                roundtrip_ok=roundtrip_ok,
            )

        except subprocess.TimeoutExpired:
            logger.warning("Timeout for genome %s", genome.uid)
            return self._error_result(float(self.timeout), "timeout")
        except Exception as e:
            logger.error("Unexpected error evaluating %s: %s", genome.uid, e)
            return self._error_result(0.0, str(e))
        finally:
            if cleanup:
                shutil.rmtree(work_dir, ignore_errors=True)

    def _verify(self, bz_path: str, work_dir: str) -> tuple:
        """Decompress and verify output is valid BAM."""
        dec_path = os.path.join(work_dir, "roundtrip.bam")

        t0 = time.monotonic()
        result = subprocess.run(
            [self.bz_bin, "decompress", "-i", bz_path, "-o", dec_path,
             "-t", str(self.threads)],
            capture_output=True, text=True, timeout=self.timeout,
        )
        dec_time = time.monotonic() - t0

        if result.returncode != 0:
            return False, dec_time

        # Check that decompressed file exists and is non-empty
        if not os.path.exists(dec_path) or os.path.getsize(dec_path) == 0:
            return False, dec_time

        return True, dec_time

    def _error_result(self, runtime_s: float, error: str) -> EvalResult:
        return EvalResult(
            fitness=0.0, original_size=self.original_size,
            compressed_size=0, ratio=0.0,
            compress_time_s=runtime_s, error=error,
        )

    def apply_result(self, genome: Genome, result: EvalResult) -> None:
        """Store evaluation result back into the genome."""
        genome.fitness = result.fitness
        genome.original_size = result.original_size
        genome.compressed_size = result.compressed_size
        genome.compress_time_s = result.compress_time_s
