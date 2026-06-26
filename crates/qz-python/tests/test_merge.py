#!/usr/bin/env python3
"""Roundtrip tests for qz.merge (the Python binding over `qz merge`).

Self-contained: generates synthetic FASTQ, so it needs no real_data. No pytest
dependency — run directly with the interpreter the `qz` extension is built into:

    maturin develop -m crates/qz-python/Cargo.toml
    python crates/qz-python/tests/test_merge.py

Exits 0 on success, non-zero (with a traceback) on the first failure.
"""

import os
import sys
import tempfile

import qz


def _write_fastq(path, n, *, start=0, read_len=150):
    """Write `n` deterministic single-end reads. `start` offsets the read ids so
    two halves are distinct, ordered, and concatenate into a contiguous file."""
    bases = "ACGT"
    with open(path, "w") as fh:
        for i in range(start, start + n):
            seq = "".join(bases[(i * 7 + j) % 4] for j in range(read_len))
            qual = "I" * read_len
            fh.write(f"@read_{i}\n{seq}\n+\n{qual}\n")


def _read_bytes(path):
    with open(path, "rb") as fh:
        return fh.read()


def test_merge_single_end_roundtrip(work):
    """compress two halves separately -> merge -> the merged archive verifies and
    decompresses to the byte-exact concatenation of the two halves, in order."""
    h1 = os.path.join(work, "half1.fastq")
    h2 = os.path.join(work, "half2.fastq")
    _write_fastq(h1, 2000, start=0)
    _write_fastq(h2, 1500, start=2000)

    a = os.path.join(work, "a.qz")
    b = os.path.join(work, "b.qz")
    qz.compress(h1, a, working_dir=work)
    qz.compress(h2, b, working_dir=work)

    merged = os.path.join(work, "merged.qz")
    qz.merge([a, b], merged)

    info = qz.verify(merged, working_dir=work)
    assert info["ok"] is True, info
    assert info["num_reads"] == 3500, info["num_reads"]

    out = os.path.join(work, "merged.fastq")
    qz.decompress(merged, out, working_dir=work)

    expected = _read_bytes(h1) + _read_bytes(h2)
    got = _read_bytes(out)
    assert got == expected, (
        f"merged decode mismatch: {len(got)} bytes vs expected {len(expected)}"
    )


def test_merge_force_overwrite(work):
    """force=True overwrites an existing output; without it the call errors."""
    h1 = os.path.join(work, "f1.fastq")
    _write_fastq(h1, 500, start=0)
    a = os.path.join(work, "fa.qz")
    qz.compress(h1, a, working_dir=work)

    merged = os.path.join(work, "fmerged.qz")
    qz.merge([a], merged)
    assert os.path.exists(merged)

    # Re-merge over an existing file without force -> RuntimeError.
    try:
        qz.merge([a], merged)
        raise AssertionError("expected merge to refuse overwriting an existing file")
    except RuntimeError:
        pass

    # With force=True it succeeds.
    qz.merge([a], merged, force=True)


def test_merge_empty_inputs_errors(work):
    """An empty input list is rejected, not silently producing an empty archive."""
    merged = os.path.join(work, "empty.qz")
    try:
        qz.merge([], merged)
        raise AssertionError("expected merge of zero inputs to error")
    except RuntimeError:
        pass


def main():
    tests = [
        test_merge_single_end_roundtrip,
        test_merge_force_overwrite,
        test_merge_empty_inputs_errors,
    ]
    failed = 0
    for t in tests:
        with tempfile.TemporaryDirectory(prefix="qz_merge_test_") as work:
            try:
                t(work)
                print(f"PASS {t.__name__}")
            except Exception as e:  # noqa: BLE001 - surface the failure + exit nonzero
                failed += 1
                print(f"FAIL {t.__name__}: {e!r}")
                import traceback

                traceback.print_exc()
    if failed:
        print(f"\n{failed} test(s) FAILED")
        sys.exit(1)
    print(f"\nAll {len(tests)} merge tests passed")


if __name__ == "__main__":
    main()
