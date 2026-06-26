# qz — Python bindings

Python bindings for **qz**, a high-performance, lossless, read-order-preserving
FASTQ/FASTA compressor. The extension covers compress, decompress (single-end and
paired), verify, and merge, and releases the GIL during the work.

```python
import qz

# Compress / decompress (single-end or paired)
qz.compress("reads.fastq", "reads.qz", quality_mode="lossless", threads=0)
qz.compress(["R1.fastq", "R2.fastq"], "sample.qz")     # paired-end (R1, R2)
qz.decompress("reads.qz", "reads.fastq")
qz.decompress("sample.qz", "sample", split=True)       # paired: sample_R1.fastq + sample_R2.fastq

# Verify an archive's integrity without writing any output (works for every
# archive type — single-end, paired, reference, and cluster). Returns a dict of
# stats on success; raises RuntimeError if the archive is corrupt or tampered.
info = qz.verify("reads.qz")            # deep: reconstruct through a CRC32
fast = qz.verify("reads.qz", fast=True) # fast: per-block CRC32 walk, no full decode
print(info["num_reads"], info["mode"])

# Stitch several archives into one without re-encoding (reads in input order).
qz.merge(["a.qz", "b.qz"], "merged.qz")

print(qz.version())
```

- `compress` takes a single FASTQ path (single-end) or a list `[R1, R2]`
  (paired-end), plus `quality_mode`, `ultra`, `fasta`, `no_quality`, `threads`,
  `working_dir`, `force`, and `fast`. A one-element list is single-end.
- `decompress` accepts `working_dir`, `threads`, `gzipped`, `gzip_level`, `force`,
  `interleaved`, and `split`. Paired and paired-reference archives reconstruct
  **two** FASTQ files: decode them with `split=True` (separate
  `<prefix>_R1.fastq` / `<prefix>_R2.fastq`, output used as a prefix) or
  `interleaved=True` (one interleaved FASTQ, R1/R2 records alternating); with
  neither they error. Single-end, single-end-reference, and single-end-cluster
  archives decode to one file normally. Paired-cluster archives use the CLI.

```python
# Two mates as separate files, then the same archive as one interleaved stream:
qz.decompress("paired.qz", "sample", split=True)       # sample_R1.fastq + sample_R2.fastq
qz.decompress("paired.qz", "interleaved.fastq", interleaved=True)
```
- `verify` accepts `fast`, `threads`, and `working_dir`.
- `merge` stitches several archives into one without re-encoding (the `qz merge`
  primitive). It accepts `reference` (the FASTA, required to merge reference
  archives, type 2/4) and `force`. All inputs must share one archive type and
  config; cluster archives cannot be merged.

The first qz call in a process initializes rayon's global thread pool, so a
different `threads` on a later call is ignored.

See the project's top-level `README.md` and `CLAUDE.md` for the full CLI,
archive formats, and performance numbers.
