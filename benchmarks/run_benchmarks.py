#!/usr/bin/env python3
"""Benchmark QZ against other FASTQ compressors.

Runs compress + decompress for each tool, captures wall time and peak RAM
via /usr/bin/time -v, verifies byte-identical roundtrip, and writes results
to stdout and an output file.

Every command is printed before execution so results are fully auditable.

Usage:
    python benchmarks/run_benchmarks.py <input.fastq> [--threads N] [-o results.txt]
    python benchmarks/run_benchmarks.py <input.fastq> --skip bzip2 pigz
"""

import argparse
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parent
QZ_BIN = REPO_DIR / "target" / "release" / "qz"
SPRING_BIN = shutil.which("spring") or str(REPO_DIR / ".pixi" / "envs" / "default" / "bin" / "spring")


@dataclass
class BenchResult:
    tool: str
    compressed_bytes: int = 0
    ratio: float = 0.0
    compress_wall: float = 0.0
    compress_ram_kb: int = 0
    decompress_wall: float = 0.0
    decompress_ram_kb: int = 0
    roundtrip_ok: bool | None = None  # None = not checked
    raw_stderr_compress: str = ""
    raw_stderr_decompress: str = ""
    errors: list[str] = field(default_factory=list)


def parse_time_v(stderr: str) -> tuple[float, int]:
    """Parse /usr/bin/time -v output for wall time (seconds) and peak RSS (kbytes)."""
    wall = 0.0
    rss = 0
    for line in stderr.splitlines():
        line = line.strip()
        m = re.match(r"Elapsed \(wall clock\) time.*:\s+(?:(\d+):)?(\d+):(\d+\.?\d*)", line)
        if m:
            hours = int(m.group(1) or 0)
            minutes = int(m.group(2))
            seconds = float(m.group(3))
            wall = hours * 3600 + minutes * 60 + seconds
        m = re.match(r"Maximum resident set size \(kbytes\):\s+(\d+)", line)
        if m:
            rss = int(m.group(1))
    return wall, rss


def file_md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_timed(cmd: list[str], stdout_file: str | None = None, label: str = "") -> tuple[subprocess.CompletedProcess, float, int]:
    """Run a command under /usr/bin/time -v. Returns (result, wall_secs, rss_kb).

    If stdout_file is given, stdout is redirected to that file (for bzip2/pigz
    which write compressed data to stdout).
    """
    full_cmd = ["/usr/bin/time", "-v"] + cmd
    print(f"  $ {shlex.join(cmd)}")
    if stdout_file:
        print(f"    > {stdout_file}")
        with open(stdout_file, "wb") as fout:
            r = subprocess.run(full_cmd, stdout=fout, stderr=subprocess.PIPE, text=True)
    else:
        r = subprocess.run(full_cmd, capture_output=True, text=True)

    wall, rss = parse_time_v(r.stderr)
    exit_code = r.returncode

    # Print key metrics immediately
    print(f"    exit={exit_code}  wall={fmt_time(wall)}  ram={fmt_ram(rss)}")
    if exit_code != 0:
        # Print last few lines of stderr for diagnostics
        err_lines = [l for l in r.stderr.splitlines() if not l.strip().startswith(("Command being timed", "User time", "System time", "Percent", "Elapsed", "Average", "Maximum", "Major", "Minor", "Voluntary", "Involuntary", "Swaps", "File system", "Socket", "Signals", "Page size", "Exit status"))]
        if err_lines:
            for line in err_lines[-5:]:
                print(f"    STDERR: {line}")

    return r, wall, rss


def fmt_time(seconds: float) -> str:
    if seconds >= 60:
        m = int(seconds) // 60
        s = seconds - m * 60
        return f"{m}:{s:04.1f}"
    return f"{seconds:.1f}s"


def fmt_ram(kb: int) -> str:
    if kb < 1024:
        return f"{kb} KB"
    mb = kb / 1024
    if mb < 1024:
        return f"{mb:.1f} MB"
    gb = mb / 1024
    return f"{gb:.1f} GB"


def fmt_bytes(b: int) -> str:
    mb = b / (1024 * 1024)
    if mb < 1024:
        return f"{mb:.0f} MB"
    gb = mb / 1024
    return f"{gb:.1f} GB"


# ── Individual benchmarks ─────────────────────────────────────────────────


def bench_qz_default(input_fastq: str, threads: int, tmpdir: str, input_size: int) -> BenchResult:
    res = BenchResult("QZ default")
    out = os.path.join(tmpdir, "default.qz")
    dec = os.path.join(tmpdir, "default_dec.fastq")

    print("\n--- QZ default (compress) ---")
    r, c_wall, c_rss = run_timed(
        [str(QZ_BIN), "compress", "-i", input_fastq, "-o", out, "-t", str(threads)])
    res.compress_wall = c_wall
    res.compress_ram_kb = c_rss
    res.raw_stderr_compress = r.stderr
    if r.returncode != 0:
        res.errors.append(f"compress failed (exit {r.returncode})")
        return res

    res.compressed_bytes = os.path.getsize(out)
    res.ratio = input_size / res.compressed_bytes if res.compressed_bytes else 0
    print(f"    size={fmt_bytes(res.compressed_bytes)}  ratio={res.ratio:.2f}x")

    print("\n--- QZ default (decompress) ---")
    r, d_wall, d_rss = run_timed(
        [str(QZ_BIN), "decompress", "-i", out, "-o", dec, "-t", str(threads)])
    res.decompress_wall = d_wall
    res.decompress_ram_kb = d_rss
    res.raw_stderr_decompress = r.stderr
    if r.returncode != 0:
        res.errors.append(f"decompress failed (exit {r.returncode})")
        return res

    print("\n--- QZ default (roundtrip check) ---")
    md5_in = file_md5(input_fastq)
    md5_out = file_md5(dec)
    res.roundtrip_ok = md5_in == md5_out
    print(f"    input  md5={md5_in}")
    print(f"    output md5={md5_out}")
    print(f"    roundtrip={'PASS' if res.roundtrip_ok else 'FAIL'}")

    return res


def bench_qz_ultra(input_fastq: str, level: int, threads: int, tmpdir: str, input_size: int) -> BenchResult:
    res = BenchResult(f"QZ ultra {level}")
    out = os.path.join(tmpdir, f"ultra{level}.qz")
    dec = os.path.join(tmpdir, f"ultra{level}_dec.fastq")

    print(f"\n--- QZ ultra {level} (compress) ---")
    r, c_wall, c_rss = run_timed(
        [str(QZ_BIN), "compress", "--ultra", str(level),
         "-i", input_fastq, "-o", out, "-t", str(threads)])
    res.compress_wall = c_wall
    res.compress_ram_kb = c_rss
    res.raw_stderr_compress = r.stderr
    if r.returncode != 0:
        res.errors.append(f"compress failed (exit {r.returncode})")
        return res

    res.compressed_bytes = os.path.getsize(out)
    res.ratio = input_size / res.compressed_bytes if res.compressed_bytes else 0
    print(f"    size={fmt_bytes(res.compressed_bytes)}  ratio={res.ratio:.2f}x")

    print(f"\n--- QZ ultra {level} (decompress) ---")
    r, d_wall, d_rss = run_timed(
        [str(QZ_BIN), "decompress", "-i", out, "-o", dec, "-t", str(threads)])
    res.decompress_wall = d_wall
    res.decompress_ram_kb = d_rss
    res.raw_stderr_decompress = r.stderr
    if r.returncode != 0:
        res.errors.append(f"decompress failed (exit {r.returncode})")
        return res

    print(f"\n--- QZ ultra {level} (roundtrip check) ---")
    md5_in = file_md5(input_fastq)
    md5_out = file_md5(dec)
    res.roundtrip_ok = md5_in == md5_out
    print(f"    input  md5={md5_in}")
    print(f"    output md5={md5_out}")
    print(f"    roundtrip={'PASS' if res.roundtrip_ok else 'FAIL'}")

    return res


def bench_spring(input_fastq: str, threads: int, tmpdir: str, input_size: int) -> BenchResult | None:
    spring = SPRING_BIN
    if not os.path.isfile(spring):
        print(f"\n  SPRING not found at {spring}, skipping")
        return None

    res = BenchResult("SPRING")
    out = os.path.join(tmpdir, "spring.spring")
    dec = os.path.join(tmpdir, "spring_dec.fastq")
    workdir = os.path.join(tmpdir, "spring_work")
    os.makedirs(workdir, exist_ok=True)

    # SPRING preserves read order by default (no -r flag = order preserved)
    print("\n--- SPRING (compress, order preserved) ---")
    r, c_wall, c_rss = run_timed(
        [spring, "-c", "-i", input_fastq, "-o", out,
         "-t", str(threads), "-w", workdir])
    res.compress_wall = c_wall
    res.compress_ram_kb = c_rss
    res.raw_stderr_compress = r.stderr
    if r.returncode != 0:
        res.errors.append(f"compress failed (exit {r.returncode})")
        return res

    res.compressed_bytes = os.path.getsize(out)
    res.ratio = input_size / res.compressed_bytes if res.compressed_bytes else 0
    print(f"    size={fmt_bytes(res.compressed_bytes)}  ratio={res.ratio:.2f}x")

    print("\n--- SPRING (decompress) ---")
    r, d_wall, d_rss = run_timed(
        [spring, "-d", "-i", out, "-o", dec,
         "-t", str(threads), "-w", workdir])
    res.decompress_wall = d_wall
    res.decompress_ram_kb = d_rss
    res.raw_stderr_decompress = r.stderr
    if r.returncode != 0:
        res.errors.append(f"decompress failed (exit {r.returncode})")
        return res

    # With --no-reorder, roundtrip should be byte-identical
    print("\n--- SPRING (roundtrip check) ---")
    md5_in = file_md5(input_fastq)
    md5_out = file_md5(dec)
    res.roundtrip_ok = md5_in == md5_out
    print(f"    input  md5={md5_in}")
    print(f"    output md5={md5_out}")
    print(f"    roundtrip={'PASS' if res.roundtrip_ok else 'FAIL'}")

    return res


def bench_bzip2(input_fastq: str, tmpdir: str, input_size: int) -> BenchResult | None:
    if not shutil.which("bzip2"):
        print("\n  bzip2 not found, skipping")
        return None

    res = BenchResult("bzip2 -9")
    out = os.path.join(tmpdir, "compressed.bz2")
    dec = os.path.join(tmpdir, "bzip2_dec.fastq")

    print("\n--- bzip2 -9 (compress) ---")
    r, c_wall, c_rss = run_timed(
        ["bzip2", "-9", "-c", input_fastq], stdout_file=out)
    res.compress_wall = c_wall
    res.compress_ram_kb = c_rss
    res.raw_stderr_compress = r.stderr
    if r.returncode != 0:
        res.errors.append(f"compress failed (exit {r.returncode})")
        return res

    res.compressed_bytes = os.path.getsize(out)
    res.ratio = input_size / res.compressed_bytes if res.compressed_bytes else 0
    print(f"    size={fmt_bytes(res.compressed_bytes)}  ratio={res.ratio:.2f}x")

    print("\n--- bzip2 -9 (decompress) ---")
    r, d_wall, d_rss = run_timed(
        ["bzip2", "-d", "-c", out], stdout_file=dec)
    res.decompress_wall = d_wall
    res.decompress_ram_kb = d_rss
    res.raw_stderr_decompress = r.stderr
    if r.returncode != 0:
        res.errors.append(f"decompress failed (exit {r.returncode})")
        return res

    print("\n--- bzip2 -9 (roundtrip check) ---")
    md5_in = file_md5(input_fastq)
    md5_out = file_md5(dec)
    res.roundtrip_ok = md5_in == md5_out
    print(f"    input  md5={md5_in}")
    print(f"    output md5={md5_out}")
    print(f"    roundtrip={'PASS' if res.roundtrip_ok else 'FAIL'}")

    return res


def bench_pigz(input_fastq: str, threads: int, tmpdir: str, input_size: int) -> BenchResult | None:
    if not shutil.which("pigz"):
        print("\n  pigz not found, skipping")
        return None

    res = BenchResult("pigz -9")
    out = os.path.join(tmpdir, "compressed.gz")
    dec = os.path.join(tmpdir, "pigz_dec.fastq")

    print("\n--- pigz -9 (compress) ---")
    r, c_wall, c_rss = run_timed(
        ["pigz", "-9", "-p", str(threads), "-c", input_fastq], stdout_file=out)
    res.compress_wall = c_wall
    res.compress_ram_kb = c_rss
    res.raw_stderr_compress = r.stderr
    if r.returncode != 0:
        res.errors.append(f"compress failed (exit {r.returncode})")
        return res

    res.compressed_bytes = os.path.getsize(out)
    res.ratio = input_size / res.compressed_bytes if res.compressed_bytes else 0
    print(f"    size={fmt_bytes(res.compressed_bytes)}  ratio={res.ratio:.2f}x")

    print("\n--- pigz -9 (decompress) ---")
    r, d_wall, d_rss = run_timed(
        ["pigz", "-d", "-p", str(threads), "-c", out], stdout_file=dec)
    res.decompress_wall = d_wall
    res.decompress_ram_kb = d_rss
    res.raw_stderr_decompress = r.stderr
    if r.returncode != 0:
        res.errors.append(f"decompress failed (exit {r.returncode})")
        return res

    print("\n--- pigz -9 (roundtrip check) ---")
    md5_in = file_md5(input_fastq)
    md5_out = file_md5(dec)
    res.roundtrip_ok = md5_in == md5_out
    print(f"    input  md5={md5_in}")
    print(f"    output md5={md5_out}")
    print(f"    roundtrip={'PASS' if res.roundtrip_ok else 'FAIL'}")

    return res


# ── Output formatting ─────────────────────────────────────────────────────


def format_table(results: list[BenchResult], input_size: int) -> str:
    """Format results as a markdown table."""
    lines = []
    lines.append("| Tool | Size (MB) | Ratio | Compress | Comp RAM | Decompress | Dec RAM | Roundtrip |")
    lines.append("|------|-----------|-------|----------|----------|------------|---------|-----------|")
    for r in results:
        sz_mb = r.compressed_bytes / (1024 * 1024)
        rt = "PASS" if r.roundtrip_ok else ("N/A" if r.roundtrip_ok is None else "**FAIL**")
        err = " (ERROR)" if r.errors else ""
        lines.append(
            f"| {r.tool}{err} | {sz_mb:.0f} | {r.ratio:.2f}x | {fmt_time(r.compress_wall)} "
            f"| {fmt_ram(r.compress_ram_kb)} | {fmt_time(r.decompress_wall)} "
            f"| {fmt_ram(r.decompress_ram_kb)} | {rt} |"
        )
    return "\n".join(lines)


def format_full_report(results: list[BenchResult], input_fastq: str, input_size: int, threads: int) -> str:
    """Format a complete report with header, table, and raw timing data."""
    now = datetime.now()
    input_mb = input_size / (1024 * 1024)

    sections = []
    sections.append(f"=== QZ Benchmark ===")
    sections.append(f"Date: {now.strftime('%Y-%m-%d %H:%M')}")
    sections.append(f"Input: {os.path.basename(input_fastq)} ({input_size:,} bytes = {input_mb:.0f} MB)")
    sections.append(f"Threads: {threads}")
    sections.append(f"All benchmarks run sequentially (no concurrent load)")
    sections.append("")

    # Summary table
    sections.append("## Results")
    sections.append("")
    sections.append(format_table(results, input_size))
    sections.append("")

    # Errors
    errors = [(r.tool, e) for r in results for e in r.errors]
    if errors:
        sections.append("## Errors")
        for tool, err in errors:
            sections.append(f"  {tool}: {err}")
        sections.append("")

    # Raw /usr/bin/time output for each tool
    sections.append("## Raw timing data")
    sections.append("")
    for r in results:
        sections.append(f"### {r.tool} — compress")
        # Extract just the /usr/bin/time lines
        for line in r.raw_stderr_compress.splitlines():
            stripped = line.strip()
            if any(stripped.startswith(k) for k in (
                "Command being timed", "User time", "System time",
                "Elapsed (wall clock)", "Maximum resident set size",
                "Exit status")):
                sections.append(f"  {stripped}")
        sections.append(f"### {r.tool} — decompress")
        for line in r.raw_stderr_decompress.splitlines():
            stripped = line.strip()
            if any(stripped.startswith(k) for k in (
                "Command being timed", "User time", "System time",
                "Elapsed (wall clock)", "Maximum resident set size",
                "Exit status")):
                sections.append(f"  {stripped}")
        sections.append("")

    return "\n".join(sections)


# ── Main ──────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark FASTQ compressors (QZ, SPRING, bzip2, pigz)")
    parser.add_argument("input", help="Input FASTQ file")
    parser.add_argument("-t", "--threads", type=int, default=os.cpu_count() or 8,
                        help="Number of threads (default: all cores)")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output results file (default: benchmarks/results_<date>.txt)")
    parser.add_argument("--ultra-levels", type=int, nargs="*", default=[1, 3, 5],
                        help="QZ ultra levels to benchmark (default: 1 3 5)")
    parser.add_argument("--skip", type=str, nargs="*", default=[],
                        help="Tools to skip (e.g. --skip spring bzip2)")
    args = parser.parse_args()

    input_fastq = os.path.abspath(args.input)
    if not os.path.isfile(input_fastq):
        print(f"Error: {input_fastq} not found", file=sys.stderr)
        sys.exit(1)

    if not QZ_BIN.is_file():
        print(f"Error: QZ binary not found at {QZ_BIN}", file=sys.stderr)
        print("Run 'cargo build --release' first", file=sys.stderr)
        sys.exit(1)

    input_size = os.path.getsize(input_fastq)
    input_mb = input_size / (1024 * 1024)
    skip = {s.lower() for s in args.skip}

    print(f"Input: {input_fastq} ({input_size:,} bytes = {input_mb:.0f} MB)")
    print(f"Input MD5: {file_md5(input_fastq)}")
    print(f"Threads: {args.threads}")
    print(f"QZ binary: {QZ_BIN}")
    print(f"SPRING binary: {SPRING_BIN}")
    print(f"Running benchmarks sequentially...")

    results: list[BenchResult] = []

    with tempfile.TemporaryDirectory(prefix="qz_bench_") as tmpdir:
        benchmarks = [
            ("pigz",   lambda: bench_pigz(input_fastq, args.threads, tmpdir, input_size)),
            ("bzip2",  lambda: bench_bzip2(input_fastq, tmpdir, input_size)),
            ("spring", lambda: bench_spring(input_fastq, args.threads, tmpdir, input_size)),
            ("qz",     lambda: bench_qz_default(input_fastq, args.threads, tmpdir, input_size)),
        ]
        for level in args.ultra_levels:
            lv = level  # capture for lambda
            benchmarks.append(
                (f"ultra{lv}", lambda lv=lv: bench_qz_ultra(input_fastq, lv, args.threads, tmpdir, input_size)),
            )

        for key, fn in benchmarks:
            if key in skip:
                print(f"\n  Skipping {key}")
                continue
            result = fn()
            if result:
                results.append(result)

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(format_table(results, input_size))

    # Write full report
    report = format_full_report(results, input_fastq, input_size, args.threads)

    if args.output:
        out_path = args.output
    else:
        out_path = str(SCRIPT_DIR / f"results_{datetime.now().strftime('%Y-%m-%d')}.txt")

    with open(out_path, "w") as f:
        f.write(report)
    print(f"\nFull report saved to {out_path}")


if __name__ == "__main__":
    main()
