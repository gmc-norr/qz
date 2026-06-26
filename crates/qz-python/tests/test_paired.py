#!/usr/bin/env python3
"""Roundtrip tests for paired (R1/R2) I/O through the qz Python bindings.

Covers the list-input compress (one path = single-end, two = paired) and the two
paired-decode shapes: split (separate `_R1.fastq`/`_R2.fastq`) and interleaved
(one stream, R1/R2 records alternating).

Self-contained: synthetic FASTQ, no pytest. Run with the interpreter the `qz`
extension is built into:

    maturin develop -m crates/qz-python/Cargo.toml
    python crates/qz-python/tests/test_paired.py
"""

import os
import sys
import tempfile

import qz


def _write_fastq(path, n, *, mate, read_len=150):
    """Write `n` deterministic reads for one mate. R1 and R2 differ in both the
    header mate field and the sequence, so a mate-swap bug fails the roundtrip."""
    bases = "ACGT"
    with open(path, "w") as fh:
        for i in range(n):
            # vary the sequence by mate so R1 != R2
            seq = "".join(bases[(i * 7 + j + mate) % 4] for j in range(read_len))
            qual = "I" * read_len
            fh.write(f"@read_{i} {mate}:N:0:0\n{seq}\n+\n{qual}\n")


def _read(path):
    with open(path, "rb") as fh:
        return fh.read()


def _records(path):
    """Yield each 4-line FASTQ record (as bytes) from `path`."""
    lines = _read(path).splitlines(keepends=True)
    for i in range(0, len(lines), 4):
        yield b"".join(lines[i : i + 4])


def test_compress_paired_split_roundtrip(work):
    """compress [R1, R2] -> paired archive -> decompress split=True yields
    `<prefix>_R1.fastq` / `<prefix>_R2.fastq` byte-identical to the inputs."""
    r1 = os.path.join(work, "in_R1.fastq")
    r2 = os.path.join(work, "in_R2.fastq")
    _write_fastq(r1, 2000, mate=1)
    _write_fastq(r2, 2000, mate=2)

    arc = os.path.join(work, "sample.qz")
    qz.compress([r1, r2], arc, working_dir=work)

    info = qz.verify(arc, working_dir=work)
    assert info["ok"] is True, info

    prefix = os.path.join(work, "out")
    qz.decompress(arc, prefix, split=True, working_dir=work)

    got1 = _read(prefix + "_R1.fastq")
    got2 = _read(prefix + "_R2.fastq")
    assert got1 == _read(r1), "R1 mismatch after paired split decode"
    assert got2 == _read(r2), "R2 mismatch after paired split decode"


def test_compress_paired_interleaved_roundtrip(work):
    """The same paired archive decodes via interleaved=True to ONE stream with
    R1/R2 records alternating (r0/1, r0/2, r1/1, r1/2, ...)."""
    r1 = os.path.join(work, "i_R1.fastq")
    r2 = os.path.join(work, "i_R2.fastq")
    _write_fastq(r1, 1000, mate=1)
    _write_fastq(r2, 1000, mate=2)

    arc = os.path.join(work, "isample.qz")
    qz.compress([r1, r2], arc, working_dir=work)

    out = os.path.join(work, "interleaved.fastq")
    qz.decompress(arc, out, interleaved=True, working_dir=work)

    expected = b"".join(
        a + b for a, b in zip(_records(r1), _records(r2))
    )
    assert _read(out) == expected, "interleaved decode mismatch"


def test_compress_single_via_one_element_list(work):
    """A one-element list is single-end (same as passing the bare string)."""
    reads = os.path.join(work, "se.fastq")
    _write_fastq(reads, 1000, mate=1)

    arc = os.path.join(work, "se.qz")
    qz.compress([reads], arc, working_dir=work)

    out = os.path.join(work, "se_out.fastq")
    qz.decompress(arc, out, working_dir=work)
    assert _read(out) == _read(reads), "single-end (1-list) roundtrip mismatch"


def test_split_on_single_end_errors(work):
    """split=True is only valid for paired archives; a single-end archive errors
    instead of silently writing one file under a confusing name."""
    reads = os.path.join(work, "se2.fastq")
    _write_fastq(reads, 500, mate=1)
    arc = os.path.join(work, "se2.qz")
    qz.compress(reads, arc, working_dir=work)

    try:
        qz.decompress(arc, os.path.join(work, "p"), split=True, working_dir=work)
        raise AssertionError("expected split=True on a single-end archive to error")
    except RuntimeError:
        pass


def test_interleaved_and_split_mutually_exclusive(work):
    """Passing both interleaved=True and split=True errors (they conflict)."""
    r1 = os.path.join(work, "m_R1.fastq")
    r2 = os.path.join(work, "m_R2.fastq")
    _write_fastq(r1, 500, mate=1)
    _write_fastq(r2, 500, mate=2)
    arc = os.path.join(work, "m.qz")
    qz.compress([r1, r2], arc, working_dir=work)

    try:
        qz.decompress(
            arc, os.path.join(work, "x"), interleaved=True, split=True, working_dir=work
        )
        raise AssertionError("expected interleaved+split to error")
    except RuntimeError:
        pass


def test_paired_without_split_or_interleaved_errors(work):
    """A paired archive through the single-output API, with neither split nor
    interleaved, still errors with a hint (unchanged behavior)."""
    r1 = os.path.join(work, "n_R1.fastq")
    r2 = os.path.join(work, "n_R2.fastq")
    _write_fastq(r1, 500, mate=1)
    _write_fastq(r2, 500, mate=2)
    arc = os.path.join(work, "n.qz")
    qz.compress([r1, r2], arc, working_dir=work)

    try:
        qz.decompress(arc, os.path.join(work, "n_out.fastq"), working_dir=work)
        raise AssertionError("expected paired decode without split/interleaved to error")
    except RuntimeError:
        pass


def main():
    tests = [
        test_compress_paired_split_roundtrip,
        test_compress_paired_interleaved_roundtrip,
        test_compress_single_via_one_element_list,
        test_split_on_single_end_errors,
        test_interleaved_and_split_mutually_exclusive,
        test_paired_without_split_or_interleaved_errors,
    ]
    failed = 0
    for t in tests:
        with tempfile.TemporaryDirectory(prefix="qz_paired_test_") as work:
            try:
                t(work)
                print(f"PASS {t.__name__}")
            except Exception as e:  # noqa: BLE001
                failed += 1
                print(f"FAIL {t.__name__}: {e!r}")
                import traceback

                traceback.print_exc()
    if failed:
        print(f"\n{failed} test(s) FAILED")
        sys.exit(1)
    print(f"\nAll {len(tests)} paired tests passed")


if __name__ == "__main__":
    main()
