# QZ · BZ

Free and open source (Apache-2.0): high-performance, lossless, read-order-preserving compressors for sequencing data, usable without restriction in academic, commercial, government, and clinical settings. Lossless and order-preserving *by default*, with read reordering (`--cluster`) and lossy quality as explicit opt-ins:

- **`qz`** compresses **FASTQ** (and FASTA) files. Its optional reference-based mode (`--reference`) uses a reference genome at compress time only, no reference is needed to decompress; it is never stored in the archive.
- **`bz`** compresses coordinate-sorted BAM files, reference-free: no reference sequence is needed to compress *or* decompress (unlike CRAM, which requires the matching reference at both ends).

Both share one foundation: columnar stream decomposition, BSC/BWT block-sorting for sequence and header data, and the [fqzcomp](https://github.com/jkbonfield/fqzcomp) context model for quality scores. Both decode through a bounded-memory streaming engine whose peak RAM is constant in read count.

> **Pre-release software.** Validate before using in production workflows. If a bug is encountered, please open an issue and include the failing archive plus the output of re-running the command with `--debug` (verbose logging, full backtraces, system diagnostics).

|                | QZ                                   | BZ                                  |
| -------------- | ------------------------------------ | ----------------------------------- |
| **Input**      | FASTQ / FASTA (file, gzip, or stdin) | Coordinate-sorted BAM               |
| **Output**     | `.qz` archive                        | `.bz` archive                       |
| **Lossless**   | Yes (lossy quality opt-in)           | Yes (lossy quality opt-in)          |
| **Ratio**      | ~22–24× with `--reference`; ~8.2× plain FASTQ | ~2.0× vs BAM (lossless, no reference) |
| **Key idea**   | Columnar streams + BSC/BWT + fqzcomp; reference mapping for max ratio | Consensus-delta + columnar + fqzcomp|
| **License**    | Apache-2.0 (open source, unrestricted use) | Apache-2.0 (open source, unrestricted use) |
| **Binary**     | `qz`                                 | `bz`                                |

## Contents

- [Features](#features)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [qz: FASTQ compression](#qz-fastq-compression)
- [bz: BAM compression](#bz-bam-compression)
- [Python API](#python-api)
- [Performance](#performance)
- [How it works](#how-it-works)
- [Troubleshooting](#troubleshooting)
- [Project layout](#project-layout)
- [Vendored dependencies](#vendored-dependencies)
- [Status & contributing](#status--contributing)
- [Acknowledgments](#acknowledgments)
- [License](#license)

## Features

- **Open source, permissively licensed.** qz and bz are released under Apache-2.0, with all source and algorithms in this repository and the bundled dependencies permissively licensed too. The license adds no academic-only, field-of-use, or commercial restrictions, so the tools can be run, modified, and redistributed anywhere, including commercial and clinical pipelines. (Several specialized sequencing-data compressors are proprietary or restrict commercial use; qz and bz do not.)
- **Lossless and order-preserving by default.** Output round-trips byte-identically and reads keep their input order. Read reordering (`--cluster`) and lossy quality are opt-in.
- **Reference-based mode (`qz --reference`).** Given a reference genome, qz maps each read and stores its mapped position and per-read edits relative to the reference, giving a FASTQ ratio of ~22–24× versus ~8× without a reference, lossless and order-preserving. The reference is used at compress time only: it is not stored in the archive and is not needed to decompress.
- **Bounded memory.** Compression and decompression stream block-by-block. Peak RAM is a function of `-t` and block size and is constant in read count, not proportional to input size.
- **Thread-scalable, NUMA-aware.** Parallelism follows `-t`. On multi-socket machines, `--numa` shards work across sockets with node-local memory; on single-socket hardware it is a no-op.
- **Self-contained archives.** One binary per tool, no runtime services; each archive is a single file. Reference mode writes a small `.sti` index cache beside the FASTA at compress time, which is not part of the archive.
- **Measured against existing tools.** In reference-based mode, qz reaches 22.16× single-end (SPRING: 21.90×) and 23.98× paired (SPRING: 24.16×), lossless and order-preserving. bz compresses BAM to ~2.00× with no reference (lossless CRAM 3.1: 1.83× default, up to 1.94× with non-default codecs, requiring a matching reference). On plain FASTQ without a reference, SPRING's read reordering reaches a higher ratio than qz's order-preserving output.

## Installation

### Requirements

- **Rust nightly**, the toolchain the project is built and tested with. (The workspace uses edition 2024, stable since Rust 1.85; nightly is the project's pinned toolchain, not an edition requirement.)
- A **C/C++ compiler with OpenMP** (e.g. GCC with `libgomp`) for `build.rs`, which compiles the bundled libbsc and htscodecs sources.
- `git`.
- *Optional:* an NVIDIA **CUDA** toolkit (GPU-accelerated BWT), and **Python ≥ 3.8** + [maturin](https://github.com/PyO3/maturin) for the Python bindings.

### Build

```bash
# 1. Clone the repository
git clone <repo-url> qz && cd qz

# 2. Build (nightly toolchain). All C/C++ dependencies (libbsc, htscodecs,
#    strobealign) are vendored in-tree; nothing else to fetch.
rustup toolchain install nightly
cargo +nightly build --release

# Binaries: target/release/qz  and  target/release/bz
```

### Optional: CUDA GPU acceleration (experimental)

The BWT suffix-array construction can run on an NVIDIA GPU via [libcubwt](https://github.com/IlyaGrebnov/libcubwt). It is experimental and per-workload. The main benefit is `qz --ultra` compress (large serial BWT blocks); at the default block size it is roughly a wash, and neutral-to-harmful elsewhere (a ~2.3× regression for `bz`, which therefore ships GPU-off by default; see [Performance](#performance)). On a CUDA build it falls back to CPU per block when no suitable GPU is available:

```bash
cargo +nightly build --release --features cuda
```

### Optional: Python bindings

```bash
pip install maturin
maturin develop --release -m crates/qz-python/Cargo.toml
```

## Quickstart

```bash
# FASTQ with a reference genome — recommended when you have one: best ratio,
# still lossless and order-preserving, and the reference is NOT needed to decompress.
qz index      ref.fa --like reads.fastq                      # one-time: build the sidecar index
qz compress   -i reads.fastq -o reads.qz --reference ref.fa
qz decompress -i reads.qz    -o reads.fastq

# FASTQ without a reference — general-purpose, one command, no setup:
qz compress   -i reads.fastq -o reads.qz
qz decompress -i reads.qz    -o reads.fastq

# Paired-end — pass -i twice (with or without --reference); decode writes _R1/_R2:
qz compress   -i R1.fastq -i R2.fastq -o sample.qz
qz decompress -i sample.qz -o sample

# Coordinate-sorted BAM — reference-free (no reference to compress or decompress):
bz compress   -i aln.bam -o aln.bz
bz decompress -i aln.bz  -o aln.bam

# Check any archive round-trips, without writing the output:
qz verify -i reads.qz
bz verify -i aln.bz
```

`--debug` (verbose logging, full backtraces, crash diagnostics) works on **every** command. `-t N` sets the thread count (`compress`, `decompress`, `verify`, `index`; default: auto). `-f` overwrites an existing output (`compress`, `decompress`, `merge`).

### Which mode should I use?

**Start with what you have:**

- **Reads with a reference genome** — resequencing, clinical panels, anything aligned to a known genome → **`qz --reference`**. This is the highest-ratio lossless option (~22–24×) and, for most resequencing work, the one you'll reach for. The reference is used only to compress; archives decode standalone.
- **FASTQ with no reference** — de novo, metagenomics, or you just want one command → **`qz` default**. Lossless, order-preserving, zero setup.
- **A coordinate-sorted BAM** → **`bz`**. Reference-free at both ends and lossless.

**Then tune for your constraint** (every qz FASTQ mode is lossless and keeps reads in order; the only exception is `--cluster`, which reorders reads while still preserving every one, see below):

| If you want…                                  | Use                              | Cost / caveat                         |
| --------------------------------------------- | -------------------------------- | ------------------------------------- |
| Best ratio on reads that map to a genome      | `qz compress --reference ref.fa` | build the index once (`qz index`)     |
| A simple, order-preserving default            | `qz compress`                    | no reference, no setup                |
| Faster compress **and** decompress            | `qz compress --fast`             | slightly worse ratio                  |
| Best ratio **without** a reference, order kept | `qz compress --ultra 1..3`      | more time and RAM                     |
| Smallest output when read order doesn't matter | `qz compress --cluster`         | output is a *permutation* of the input |
| Paired-end reads                              | a second `-i` (default, `--reference`, or `--cluster`) | decode writes `_R1`/`_R2`; `--ultra` is single-end only |
| A coordinate-sorted BAM                       | `bz compress`                    | reference-free, lossless              |
| A smaller BAM and quality may change          | `bz compress --reduce-quality N` | lossy: quality only                   |
| More throughput on a multi-socket machine     | add `--numa auto`                | no-op on single-socket hardware       |

Rule of thumb: **have a reference → `qz --reference`; don't → `qz` default (or `--ultra` for a bit more, `--cluster` if you can drop read order); BAM → `bz`.**

### What "lossless" means, and how to check it

- **`qz` default / `--fast` / `--ultra` / `--reference`**: byte-identical FASTQ round-trip.
- **`qz --cluster`**: *set-lossless*: every read is preserved, but **read order changes**. Validate with `qz verify` (a multiset checksum), **not** `md5sum` against the original.
- **`bz` (default, lossless)**: the BAM *record content* round-trips exactly; the BGZF container bytes may differ (re-compressed). Compare with `samtools view`, not `cmp` on the `.bam`.
- **`bz --reduce-quality`**: *lossy*: every field except quality is preserved; quality values change by design.
- **`qz verify` / `bz verify`** check archive integrity without writing the output file, not an external `samtools view` equality check. `qz verify` (deep) reconstructs and CRC-checks the FASTQ; `bz verify` decompresses and structurally validates every stream, while the BAM round-trip CRC is enforced at `bz decompress`.

## qz: FASTQ compression

`qz` has five subcommands: `compress`, `decompress`, `verify`, `index`, and `merge`.

### Default compression

```bash
qz compress -i reads.fastq -o reads.qz
```

Headers, DNA sequences, and quality scores are separated into independent streams and compressed in parallel: sequences and headers with libbsc (BWT + LZP + adaptive QLFC), qualities with fqzcomp. Reads stay in input order. Input may be a plain or gzipped FASTQ, or `-` / a pipe for stdin.

### Quality handling

```bash
qz compress -i reads.fastq -o reads.qz -q lossless     # default: exact quality
qz compress -i reads.fastq -o reads.qz -q illumina-bin # lossy: Illumina 8-bin remap
qz compress -i reads.fastq -o reads.qz -q discard      # drop quality entirely
qz compress -i reads.fastq -o reads.qz --no-quality    # same as discard
qz compress -i ref.fasta   -o ref.qz   --fasta         # FASTA input (no quality)
```

### Fast mode

```bash
qz compress -i reads.fastq -o reads.qz --fast
```

Uses libbsc's static QLFC entropy coder for the header/sequence streams instead of the adaptive coder. Measurably faster to compress and decompress, at a small ratio cost.

### Ultra mode

```bash
qz compress -i reads.fastq -o reads.qz --ultra      # auto level (by RAM/cores)
qz compress -i reads.fastq -o reads.qz --ultra 3    # explicit level 1–3
```

Ultra uses **larger BWT blocks** (more context → better ratio). That is its entire advantage, there is no read reordering. Levels 1–3 trade increasing time and memory for a few tenths of a × more ratio; `--ultra` with no value auto-selects the highest level that fits available RAM.

### Paired-end

```bash
qz compress   -i R1.fastq -i R2.fastq -o sample.qz
qz decompress -i sample.qz -o sample           # writes sample_R1.fastq + sample_R2.fastq
```

Specify `-i` twice to compress a mate pair into one archive. On decode, paired archives use the `-o` value as a **prefix** (`sample` → `sample_R1.fastq` / `sample_R2.fastq`).

### Reference-based mode

For reads that map to a known genome — resequencing, clinical panels, any align-to-a-reference workflow — this is qz's highest-ratio lossless option: **~22–24×**, versus ~8× for the reference-free default, and the mode most resequencing pipelines will want. Instead of compressing each read's bases directly, qz aligns the read to the reference and stores only its mapped position, strand, and the handful of edits (substitutions/indels) where it differs. Reads that **don't** map are kept verbatim, so the output is **fully lossless and order-preserving** no matter how well the reads align — poorly-mapping reads just save less.

The reference is used **only at compress time**. It is never written into the archive and is not needed (it isn't even accepted) at decompress time, so an archive is portable to anyone without shipping or re-fetching the genome.

```bash
qz index      ref.fa --like reads.fastq                                  # one-time: build the sidecar index
qz compress   -i reads.fastq -o reads.qz --reference ref.fa              # single-end
qz compress   -i R1.fastq -i R2.fastq -o sample.qz --reference ref.fa    # paired-end
qz decompress -i reads.qz -o reads.fastq                                 # NO reference needed to decode
```

- **Indexing.** `qz index <FASTA>` builds a strobealign sidecar (`<FASTA>.qz-r<LEN>.sti`) once per reference + read length. Set the length with `-r <LEN>`, or `--like reads.fastq` to infer it from the first read; with neither it uses the 150 bp profile. If no usable index exists, `qz compress --reference` errors with the exact `qz index` command to run — or pass `--build-index` to build it inline.
- **Single- and paired-end.** One `-i` for single-end, two for paired; both decode standalone (paired writes `_R1` / `_R2`). Cannot be combined with `--cluster` or `--ultra`.
- **Tuning.** `--reference-window N` (default 4) trades memory for encode speed (byte-identical output); `--reference-fast` uses sparser seeds for ~5–10% faster mapping at a negligible ratio cost.
- **Changed the FASTA?** Delete the stale `ref.fa.qz-r*.sti` and re-run `qz index`.

### Cluster mode (order-dropping)

```bash
qz compress -i reads.fastq -o reads.qz --cluster          # balanced
qz compress -i reads.fastq -o reads.qz --cluster max
```

Reorders reads by minimizer cluster for a better sequence ratio. The output is a permutation of the input, use it only when read order doesn't matter downstream. Set-losslessness is verified by a multiset checksum. Incompatible with `--ultra` and `--reference`.

### NUMA sharding

```bash
qz compress   -i reads.fastq -o reads.qz --numa auto   # default
qz decompress -i reads.qz    -o reads.fastq --numa 2   # force 2 worker processes
```

On real multi-socket NUMA hardware, `auto` launches one pinned worker process per socket (node-local memory) for higher throughput; off NUMA hardware or for small inputs it runs in-process with zero overhead. `off` forces in-process. Accepted by `compress` and `decompress`.

### Merge

```bash
qz merge part1.qz part2.qz part3.qz -o all.qz
```

Concatenates compatible single-end, paired, or reference archives (in output order) into one. Reference archives require `--reference ref.fa` (the per-shard coverage globals are re-derived over the union). This is the primitive behind NUMA-sharded compression.

### Verify

```bash
qz verify -i reads.qz            # deep:  reconstruct through a CRC32, no FASTQ written
qz verify -i reads.qz --fast     # fast:  per-block CRC32 walk, no full decode
```

`verify` checks an archive without writing the output file, and reports the read count, encoding, and per-stream compressed sizes. Two depths:

- **Deep** (default) decodes the whole archive through a CRC32 hasher instead of to disk, then prints the checksum of the reconstructed FASTQ, a per-mate CRC (R1 / R2) for the byte-lossless two-file archives (paired and paired-reference). It is the end-to-end "does this still round-trip" check, and like decode it streams, so its memory stays constant in read count.
- **`--fast`** re-checks the stored per-block CRC32 of every frame in every stream. It covers the whole archive and catches bit-rot or truncation in a fraction of the deep time, without rebuilding any FASTQ.
- **`--cluster`** archives reorder reads, so a byte CRC against the original is meaningless. Deep verify instead folds a 128-bit, order-independent multiset checksum over the reconstructed records and compares it to the value computed at compress time. A match proves *set-losslessness* (every read present exactly once, and for paired data correctly mated), which a positional CRC cannot.

### Piping

`qz` reads stdin and writes stdout via `-` (or a pipe), in **either** direction:

```bash
zcat reads.fastq.gz | qz compress -i - -o reads.qz        # stdin  → file
qz compress   -i reads.fastq -o - > reads.qz              # file   → stdout
zcat reads.fastq.gz | qz compress -i - -o - > reads.qz    # stdin  → stdout (full pipe)

qz decompress -i reads.qz -o - | head                     # file   → stdout
qz decompress -i reads.qz -o out.fastq.gz --gzipped       # gzip the output
```

**Compression streams natively in both directions**. The encoder never seeks backward (the footer locator sits at EOF), so `-o -` to a pipe is byte-identical to writing a file. **Decompression** writes to a pipe freely, but its *input* must be random-access to read that trailing footer: a seekable file is read in place, while a piped-in archive (`qz decompress -i -`) is first spooled to a temp file in the working directory. So when decoding, prefer giving the archive as a real file path rather than piping it in.

#### Interleaved paired output (`--interleaved`)

A **paired** or **paired-reference** archive reconstructs *two* FASTQ files, which a single stdout can't carry. `--interleaved` merges the mates into one stream, R1/R2 records alternating (`r0/1, r0/2, r1/1, …`) so a paired archive can pipe straight into an aligner that accepts interleaved input:

```bash
qz decompress -i sample.qz --interleaved -o - | bwa mem -p ref.fa - | samtools sort -o sample.bam
qz decompress -i sample.qz --interleaved -o - | strobealign ref.fa --interleaved -
qz decompress -i sample.qz --interleaved -o reads.interleaved.fastq.gz --gzipped   # one gzipped file
```

The single output goes to a file, `-`/stdout, or (`--gzipped`) a gzip stream. Mate identity rides on record adjacency plus the original `/1`,`/2` names, which is what `bwa mem -p`, `bowtie2 --interleaved`, `strobealign --interleaved`, `fastp --interleaved_in`, and BBTools (`int=t`) expect. `--interleaved` runs in-process (it doesn't combine with `--numa` sharding) and only applies to two-mate archives: a single-end, cluster, or single-end-reference archive is rejected with a message.

### Decompression output (`-o`) by archive type

| Archive               | `qz decompress -o` writes                            |
| --------------------- | ---------------------------------------------------- |
| single-end            | one exact output file                                |
| single-end reference  | one exact output file                                |
| single-end cluster    | one exact output file                                |
| paired-end            | a **prefix** → `_R1.fastq` / `_R2.fastq`             |
| paired reference      | a **prefix** → `_R1.fastq` / `_R2.fastq`             |
| paired cluster        | **two** exact `-o` paths (`-o R1.fastq -o R2.fastq`) |
| paired / paired reference, `--interleaved` | **one** stream (file, `-`, or `--gzipped`), R1/R2 interleaved |

Passing the wrong number of `-o` values is rejected, not silently ignored. `--interleaved` turns the two paired modes into a single interleaved stream (see [Piping](#piping)).

### qz CLI reference

| Command      | Purpose                                                       |
| ------------ | ------------------------------------------------------------- |
| `compress`   | Compress FASTQ/FASTA → `.qz` (single, paired, reference, ultra, cluster) |
| `decompress` | Decompress `.qz` → FASTQ (optionally gzipped)                 |
| `verify`     | Check archive integrity (`--fast` for CRC-only)              |
| `index`      | Build a strobealign reference index for reference mode        |
| `merge`      | Merge compatible archives into one                            |

## bz: BAM compression

`bz` compresses **coordinate-sorted** BAM. It has four subcommands: `compress`, `decompress`, `extract`, and `verify`.

### Compress / decompress

```bash
bz compress   -i aln.bam -o aln.bz
bz decompress -i aln.bz  -o aln.bam
```

bz builds a local consensus of the aligned sequences and stores each read as its difference from that consensus (mostly zeros), columnarizes the BAM fields (positions, flags, MAPQ, mate info, CIGAR, tags), and compresses qualities with fqzcomp. Round-trips byte-identically on the BGZF-decompressed BAM content.

### Compression level

```bash
bz compress -i aln.bam -o aln.bz --level 2
```

`--level 0-3` is a memory/speed preset (ratio is roughly flat across levels): `0` = auto (downshifts under low RAM), `1` = lowest memory, `2` = balanced/fastest, `3` = largest chunks / most memory.

### Lossy quality reduction (opt-in)

```bash
bz compress -i aln.bam -o aln.bz --reduce-quality 2
bz compress -i aln.bam -o aln.bz --reduce-quality 3 --quality-scheme coarse8
```

`--reduce-quality 1-3` flattens quality at confident reference-matching positions while preserving full resolution at candidate-variant sites, for a large additional size reduction. `--quality-scheme` picks the flatten scheme (`twobin` default, `coarse8`, or `single`). Off by default, output is lossy only when asked for.

### Extract FASTQ from BAM

```bash
bz extract -i aln.bam -o sample
# → sample_R1.qz + sample_R2.qz  (paired)   or   sample_SE.qz  (single-end)
```

`extract` pulls reads out of a BAM and writes them directly as `qz` FASTQ archives.

### Verify and NUMA

```bash
bz verify -i aln.bz
bz compress -i aln.bam -o aln.bz --numa auto
```

`bz verify` is an integrity check without writing a BAM: it re-checks each chunk's compressed-payload CRC32, decompresses every stream (BSC decode fails on corrupt data), validates the columnar structure, and walks the inner fqz per-block CRCs. It does **not** reconstruct records, so it is not by itself a round-trip proof. That guarantee is the per-chunk content CRC, which `bz decompress` recomputes over the rebuilt records and rejects on mismatch. (bz verify has no `--fast`/`--deep` split; it is one mode.)

`--numa auto|off|N` shards across sockets exactly as in qz; the archive is byte-identical to a single-process compress.

### bz CLI reference

| Command      | Purpose                                                |
| ------------ | ------------------------------------------------------ |
| `compress`   | Compress coordinate-sorted BAM → `.bz`                 |
| `decompress` | Decompress `.bz` → BAM                                 |
| `extract`    | Extract FASTQ from a BAM into `qz` archives            |
| `verify`     | Check archive integrity without writing to disk        |

## Python API

The `qz` Python package (built with maturin) exposes `compress` / `decompress` / `verify` / `merge`:

```python
import qz

qz.compress("reads.fastq", "reads.qz", quality_mode="lossless", threads=0)
qz.compress(["R1.fastq", "R2.fastq"], "sample.qz")    # paired-end (R1, R2)
qz.decompress("reads.qz", "reads.fastq")
qz.decompress("sample.qz", "sample", split=True)      # paired → sample_R1.fastq + sample_R2.fastq

info = qz.verify("reads.qz")              # integrity check → stats dict (raises if corrupt)
print(info["num_reads"], info["mode"])    # e.g. 1000000 'deep'
qz.verify("reads.qz", fast=True)          # per-block CRC32 walk, no full decode

qz.merge(["a.qz", "b.qz"], "merged.qz")   # stitch archives into one (reads in input order)

print(qz.version())
```

`compress` takes a single FASTQ path or `[R1, R2]` for paired-end, plus `quality_mode`, `ultra`, `fasta`, `no_quality`, `threads`, `working_dir`, `force`, and `fast`. `decompress` accepts `working_dir`, `threads`, `gzipped`, `gzip_level`, `force`, `interleaved`, and `split`. `verify` accepts `fast`, `threads`, and `working_dir`; on success it returns a stats dict (read count, mode, encoding, per-stream sizes, CRC32s, timing), and raises `RuntimeError` if the archive is corrupt, truncated, tampered with, or fails any check. `merge` stitches several archives into one without re-encoding (the `qz merge` primitive); it accepts `reference` (the FASTA, required for reference archives type 2/4) and `force`. All inputs must share one archive type and config; cluster archives cannot be merged.

Paired and paired-reference archives reconstruct two FASTQ files: decode them with `split=True` (separate `<prefix>_R1.fastq` / `<prefix>_R2.fastq`, using the output as a prefix) or `interleaved=True` (one interleaved stream); with neither they error, since the API is otherwise single-output. Single-end, single-end-reference, and single-end-cluster archives decode to one file normally. Paired-cluster archives use the CLI. `verify` works for every archive type, paired and reference included, because it never writes output.

The GIL is released during the work. The first qz call in a process initializes rayon's global thread pool, so a later call with a different `threads` is ignored.

## Performance

Measured 2026-06-24 on an Intel Xeon Gold 6254 (2× 18 cores / 36 cores / 72 threads, 2 NUMA sockets), 376 GB RAM, Linux, `-t 72`, idle box. Single-run numbers (treat ±5% as noise). Every row round-trips and is verified: md5 of the reconstructed FASTQ (`OK`); `qz verify` / sorted-record multiset for the permutation outputs `--cluster` and `SPRING -r` (`PASS`); `samtools view` fingerprint for bz (`OK`; lossy `--reduce-quality` = `QUAL_ONLY_DIFF`, everything but quality is byte-identical); `n-a` for CRAM (competitor, not re-verified).

Two FASTQ tiers: **shallow** single-end (`ERR3239334`, 30M reads, 150 bp WGS, 10.99 GB raw) and **deep** chr20 (HG002 NovaSeq 2×151) for paired + reference, plus a BAM tier. qz's highest lossless ratios come from reference-based mode: 22.16× single-end / 23.98× paired on deep chr20, compared with SPRING-single at 21.90× and SPRING-paired at 24.16×, while staying lossless *and* order-preserving (SPRING reorders). The shallow single-end tier (8.18×) is qz's order-preserving result on low-coverage WGS. Every table reports size, ratio, and both compress *and* decompress wall-time, peak RSS, and CPU%, so qz↔SPRING and bz↔CRAM compare directly.

**Columns:** `size` = compressed MB · `ratio` = raw ÷ compressed · `c s`/`d s` = compress/decompress wall seconds · `c GB`/`d GB` = peak RSS · `c%`/`d%` = `/usr/bin/time -v` CPU% (**7200% = all 72 threads saturated**) · `RT` = round-trip.

### Single-end FASTQ: 30M reads, 150 bp WGS (10.99 GB raw)

| tool | size | ratio | c s | c GB | c% | d s | d GB | d% | RT |
| ---- | ---- | ----- | --- | ---- | -- | --- | ---- | -- | -- |
| **qz default** | 1281 MB | **8.18×** | 33.5 | 6.3 | 2047% | 18.5 | 2.7 | 2680% | OK |
| qz default `--numa auto` | 1281 MB | 8.18× | **23.5** | 6.3 | 3381% | 17.6 | 3.2 | 2789% | OK |
| qz fast | 1285 MB | 8.15× | 32.9 | 8.0 | 1825% | 14.4 | 3.0 | 2532% | OK |
| qz ultra 1 | 1250 MB | 8.38× | 192.8† | 2.8 | 507% | 17.4 | 5.4 | 4944% | OK |
| qz ultra 2 | 1242 MB | 8.44× | 197.2† | 4.5 | 503% | 22.6 | 7.0 | 3740% | OK |
| qz ultra 3 | 1226 MB | **8.55×** | 199.1† | 6.4 | 502% | 20.8 | 11.8 | 4285% | OK |
| qz cluster *(reorders)* | 1041 MB | **10.07×** | 181.3 | 11.6 | 2291% | 67.6 | 8.3 | 327% | PASS |
| SPRING (order-preserving) | 1129 MB | 9.28× | 186.1 | 21.2 | 3898% | 40.7 | 19.0 | 2486% | OK |
| SPRING `-r` *(reorders)* | 1040 MB | 10.08× | 234.2 | 14.1 | 2995% | 41.9 | 18.9 | 2563% | PASS |
| pigz -9 | 2083 MB | 5.03× | 31.7 | 0.0 | 6466% | 24.2 | 0.0 | 214% | OK |

† Ultra compress time scales with the larger BWT blocks (hardware/load-dependent). **Order-preserving:** SPRING 9.28× vs qz default 8.18×, at ~3.4× the compress RAM (21.2 vs 6.3 GB). **Reorder:** qz `--cluster` 10.07× vs SPRING `-r` 10.08×, at lower decode RAM (8.3 vs 18.9 GB).

### Paired-end FASTQ: deep chr20 (HG002 NovaSeq 2×151, 7.47 GB raw = R1+R2)

| tool | size | ratio | c s | c GB | c% | d s | d GB | d% | RT |
| ---- | ---- | ----- | --- | ---- | -- | --- | ---- | -- | -- |
| qz default (paired) | 863 MB | 8.26× | 29.2 | 8.7 | 1803% | 17.2 | 5.2 | 2198% | OK |
| qz cluster *(reorders)* | 484 MB | 14.72× | 102.4 | 16.4 | 2497% | 63.5 | 4.2 | 169% | PASS |
| SPRING (order-preserving) | 295 MB | **24.17×** | 65.6 | 13.1 | 1658% | 24.9 | 15.0 | 1547% | OK |
| SPRING `-r` *(reorders)* | 295 MB | 24.15× | 89.6 | 8.0 | 1186% | 26.1 | 15.0 | 1412% | PASS |
| pigz -9 *(baseline)* | 1391 MB | 5.12× | 32.6 | 0.0 | 6945% | 16.3 | 0.0 | 209% | OK |

SPRING's consensus-reorder reaches the highest paired-FASTQ ratio here (even order-preserving, via a stored permutation). For aligned reads, qz's lossless option is reference mode ↓.

### Reference-based FASTQ: deep chr20 (reads mapped to GRCh38 chr20; reference used at compress time only)

Single rows: 3.74 GB raw (R1). Paired rows: 7.47 GB raw (R1+R2).

| tool | size | ratio | c s | c GB | c% | d s | d GB | d% | RT |
| ---- | ---- | ----- | --- | ---- | -- | --- | ---- | -- | -- |
| **qz `--reference` (single R1)** | 161 MB | **22.16×** | 29.4 | 4.0 | 1918% | 10.7 | 1.7 | 952% | OK |
| qz `--reference-fast` (single) | 161 MB | 22.15× | 32.9 | 3.0 | 1495% | 10.6 | 1.6 | 950% | OK |
| qz cluster (single, reorders) | 197 MB | 18.07× | 52.8 | 10.0 | 1662% | 31.5 | 5.9 | 168% | PASS |
| SPRING (single, order-preserving) | 163 MB | 21.90× | 32.7 | 12.4 | 1709% | 12.7 | 10.6 | 1608% | OK |
| SPRING `-r` (single, reorders) | 163 MB | 21.85× | 54.4 | 6.0 | 1018% | 12.7 | 10.6 | 1614% | PASS |
| pigz -9 (single, baseline) | 682 MB | 5.23× | 15.1 | 0.0 | 6953% | 8.0 | 0.0 | 211% | OK |
| qz `--reference` (paired) | 297 MB | **23.98×** | 43.4 | 4.9 | 2340% | 22.9 | 2.1 | 831% | OK |
| qz `--reference-fast` (paired) | 297 MB | 23.99× | 42.3 | 5.6 | 2072% | 20.7 | 2.2 | 928% | OK |
| SPRING (paired, order-preserving) | 295 MB | 24.16× | 61.2 | 13.1 | 1787% | 25.4 | 15.0 | 1465% | OK |
| pigz -9 (paired, baseline) | 1391 MB | 5.12× | 32.6 | 0.0 | 6945% | 16.3 | 0.0 | 209% | OK |

qz reference is comparable to SPRING on ratio (22.16× vs SPRING-single's 21.90×; 23.98× vs SPRING-paired's 24.16×), while staying lossless and order-preserving, at ~1/3 the compress RAM and ~5–7× less decode RAM (1.7–2.2 vs 10.6–15.0 GB), with comparable-or-faster decode. Generic `pigz -9` on the same reads reaches ~5.2× (single) / ~5.1× (paired), about a quarter of reference mode's ratio.

### BAM: HG002 GRCh38 chr20 (2.72 GB raw; `bz` vs samtools CRAM)

| tool | size | ratio | c s | c GB | c% | d s | d GB | d% | RT |
| ---- | ---- | ----- | --- | ---- | -- | --- | ---- | -- | -- |
| **bz** (lossless) | 1294 MB | **2.00×** | 28.6 | 13.2 | 1926% | 28.3 | 16.5 | 1862% | OK |
| bz `-l 1` | 1295 MB | 2.00× | 46.3 | 5.6 | 1114% | 23.9 | 14.5 | 2308% | OK |
| bz `-l 2` | 1294 MB | 2.00× | 26.2 | 13.9 | 2190% | 28.2 | 16.6 | 1816% | OK |
| bz `-l 3` | 1293 MB | 2.01× | 39.9 | 19.9 | 1523% | 39.0 | 17.4 | 1314% | OK |
| CRAM 3.0 (`samtools`) | 1467 MB | 1.77× | 21.1 | 1.0 | 587% | 9.5 | 2.1 | 3272% | n-a |
| CRAM 3.1 (`samtools`, default) | 1418 MB | 1.83× | 16.1 | 1.4 | 921% | 9.3 | 2.1 | 3305% | n-a |
| CRAM 3.1 (`use_fqz,use_tok,use_arith`) | 1340 MB | 1.94× | 21.2 | 6.9 | 3147% | 11.2 | 6.5 | 4281% | n-a |
| pigz -9 *(re-gzip baseline)* | 2595 MB | 1.00× | 2.1 | 0.0 | 4485% | 3.7 | 0.0 | 221% | OK |
| bz `--reduce-quality 1` *(lossy)* | 508 MB | 5.10× | 29.8 | 12.9 | 1864% | 20.7 | 15.5 | 2194% | QUAL_ONLY_DIFF |
| bz `--reduce-quality 2` *(lossy)* | 407 MB | 6.37× | 31.2 | 11.9 | 1733% | 19.7 | 16.4 | 2128% | QUAL_ONLY_DIFF |
| bz `--reduce-quality 3` *(lossy)* | 374 MB | 6.93× | 32.1 | 12.3 | 1691% | 18.6 | 18.4 | 2270% | QUAL_ONLY_DIFF |

bz reaches a higher ratio than lossless CRAM and needs **no reference**: 2.00× vs CRAM 3.1's 1.83× out-of-the-box, or 1.94× with CRAM 3.1's strongest codecs enabled (`use_fqz`/`use_tok`/`use_arith`). That 1.83→1.94× lift is **entirely `use_fqz`** (the name tokeniser and range coder add ~nothing here), and `use_fqz` is *fqzcomp*, the same htscodecs quality codec bz uses by default. So on quality the two are comparable; bz's ratio difference comes from its non-quality streams (consensus-delta sequence + columnar fields).

The no-reference property matters here. This BAM carries no per-contig M5 tags and was aligned to GRCh38 + hs38d1 decoys, so CRAM reaches the ratios above only when handed a reference whose contigs match the header *exactly*. Given a mismatched reference, samtools silently embeds the reference (`embed_ref`), which inflates the CRAM to ~1.29×; bz, needing no reference, avoids this. CRAM, in turn, **decodes ~3× faster** (9–11 s vs bz's ~28 s) and uses far less compress RAM.

A generic `pigz -9` re-gzip of the BGZF BAM lands at ~1.00× (no gain: the data is already DEFLATE-compressed), which is why a BAM-aware tool is needed. Levels are a memory/speed dial (ratio ~flat); quality dominates the archive, so opt-in lossy `--reduce-quality` reaches 6.93×.

*The CRAM rows were re-measured 2026-06-25 against a contig-matching reference (GRCh38 + hs38d1); an earlier pass had fed CRAM a decoy-less reference, triggering `embed_ref` and understating CRAM at 1.25–1.29×. The benchmark harness now aborts any CRAM row that falls back to `embed_ref`. The BAM is fully coordinate-sorted at the record level (only its `@HD SO:` label was stale); sort order is not a factor in these ratios.*

### NUMA: `--numa off` vs `auto` (the shipped multi-socket default)

One workload per mode, compressed **and** decoded at each NUMA mode, so the off→auto delta isolates the 2-socket sharding win. `--numa auto` is a zero-overhead in-process fallback on single-socket hardware.

| mode | numa | ratio | c s | c GB | c% | d s | d GB | d% |
| ---- | ---- | ----- | --- | ---- | -- | --- | ---- | -- |
| single 30M | off | 8.18× | 34.3 | 8.3 | 2030% | 26.9 | 3.4 | 1921% |
| single 30M | **auto** | 8.18× | **22.9** | 6.2 | 3460% | **18.3** | 2.8 | 2694% |
| paired chr20 | off | 8.26× | 28.6 | 8.4 | 1826% | 21.6 | 5.3 | 1536% |
| paired chr20 | **auto** | 8.26× | **24.7** | 7.8 | 2677% | **12.3** | 5.2 | 3152% |
| ref single | off | 22.16× | 27.9 | 4.1 | 2011% | 18.7 | 2.0 | 699% |
| ref single | **auto** | 22.11× | **16.7** | 4.2 | 2911% | **10.4** | 1.7 | 975% |
| ref paired | off | 23.98× | 48.2 | 5.0 | 2110% | 40.6 | 2.3 | 623% |
| ref paired | **auto** | 23.96× | **43.7** | 3.5 | 2121% | **20.7** | 2.0 | 925% |
| bz chr20 | off | 2.00× | 27.6 | 14.0 | 1998% | 42.7 | 25.0 | 1460% |
| bz chr20 | **auto** | 2.00× | 29.7 | 12.5 | 2081% | **27.2** | 16.5 | 1903% |

**NUMA helps decode more than compress**: off→auto **1.47–1.96×** on decode across every mode (single 1.47×, paired 1.75×, ref-single 1.80×, ref-paired 1.96×, bz 1.57×), vs 1.10–1.67× on compress (bz compress slightly regresses to 0.93×, because the 2.7 GB BAM is near the sharding threshold). BSC decode is memory-bandwidth-bound, so the second socket's added memory channels account for the gain.

**Memory.** qz's peak RSS is bounded by the working set and **constant in read count** (not file size): single-end compress ~6.3 GB at `-t 72`, scaling down with `-t` (~3.1 GB at `-t 8`); reference decode 1.7–2.2 GB. SPRING used 12–21 GB on the same inputs. Full per-row data: [`benchmarks/results/SUMMARY-2026-06-24.md`](benchmarks/results/SUMMARY-2026-06-24.md).

### GPU (CUDA, opt-in): where it helps, where it doesn't

GPU is compile-time (`--features cuda`): libbsc's BWT (forward for compress, inverse for decode) runs on the GPU via vendored `libcubwt`; everything else stays on CPU. On such a build the GPU is **runtime-controllable**: `--gpu auto|off|require` on `qz`/`bz` `compress`/`decompress` (`auto` = GPU BWT with per-block CPU fallback, `off` = force CPU, `require` = error at startup if no usable CUDA device), backed by the `QZ_GPU` env var (which propagates to `--numa` worker processes). **`qz` defaults to `auto`; `bz` defaults to `off`** (GPU regresses bz; see the takeaways below). The flag is absent from CPU-only builds. Measured **CPU build vs `--features cuda` build** on an **RTX 2080 Ti (11 GB)** (sm_75, CUDA 13.3). Output is **byte-identical** to the CPU build (GPU only changes the BWT backend), so ratios are unchanged. The default shipped build is CPU-only and portable. `c%`/`d%` = `/usr/bin/time -v` CPU% (7200% = 72 threads); `RT` = round-trip.

**Single-end matrix.** 5M (1.83 GB) + 30M (10.99 GB), `default` + `ultra-2`, four configs each:

| config | ratio | c s | c GB | c% | d s | d GB | d% | RT |
| ------ | ----- | --- | ---- | -- | --- | ---- | -- | -- |
| 5M default · CPU | 8.15× | 7.12 | 5.4 | 1672% | 6.60 | 3.2 | 1318% | OK |
| 5M default · NUMA | 8.15× | 5.87 | 4.5 | 2130% | 3.43 | 2.5 | 2411% | OK |
| 5M default · GPU | 8.15× | 5.94 | 4.3 | 1300% | 5.00 | 2.4 | 1520% | OK |
| 5M default · GPU+NUMA | 8.15× | **5.15** | 3.4 | 1530% | 3.53 | 1.4 | 2049% | OK |
| 5M ultra-2 · CPU | 8.40× | 35.61 | 4.5 | 476% | 8.72 | 3.8 | 1544% | OK |
| 5M ultra-2 · NUMA | 8.40× | 33.76 | 4.6 | 486% | 6.38 | 3.0 | 1882% | OK |
| 5M ultra-2 · GPU | 8.40× | 21.49 | 3.5 | 357% | 5.90 | 2.4 | 1326% | OK |
| 5M ultra-2 · GPU+NUMA | 8.40× | **21.42** | 3.4 | 360% | **4.55** | 1.7 | 1664% | OK |
| 30M default · CPU | 8.18× | 34.81 | 7.2 | 1960% | 27.57 | 3.4 | 1853% | OK |
| 30M default · NUMA | 8.18× | 22.90 | 5.6 | 3431% | 16.66 | 2.8 | 2983% | OK |
| 30M default · GPU | 8.18× | 34.20 | 8.1 | 1361% | 23.46 | 2.9 | 1929% | OK |
| 30M default · GPU+NUMA | 8.18× | **19.12** | 5.1 | 2386% | **16.10** | 2.3 | 2706% | OK |
| 30M ultra-2 · CPU | 8.44× | 195.19 | 4.5 | 503% | 32.33 | 7.3 | 2570% | OK |
| 30M ultra-2 · NUMA | 8.44× | 195.43 | 4.5 | 499% | 19.96 | 7.0 | 4221% | OK |
| 30M ultra-2 · GPU | 8.44× | **114.30** | 3.4 | 382% | 23.85 | 4.7 | 2024% | OK |
| 30M ultra-2 · GPU+NUMA | 8.44× | 113.53 | 3.4 | 380% | **16.52** | 4.1 | 3261% | OK |

> On 30M-default decode, NUMA (16.66 s) and GPU+NUMA (16.10 s) are effectively **equal**: once NUMA has captured the memory-bandwidth gain, the GPU adds little.

**Other modes.** bz (chr20 BAM, 2.72 GB), paired `--cluster` (7.47 GB), reference single/paired (chr20):

| config | ratio | c s | c GB | c% | d s | d GB | d% | RT |
| ------ | ----- | --- | ---- | -- | --- | ---- | -- | -- |
| bz · CPU | 2.00× | **26.3** | 13.5 | 2035% | 36.9 | 25.8 | 1623% | OK |
| bz · NUMA | 2.00× | 27.6 | 11.8 | 2270% | **28.2** | 15.9 | 1883% | OK |
| bz · GPU | 2.00× | 59.4 | 11.9 | 715% | 38.7 | 25.2 | 1427% | OK |
| bz · GPU+NUMA | 2.00× | 64.6 | 11.5 | 785% | 31.4 | 14.2 | 1443% | OK |
| cluster · CPU | 14.72× | 124.0 | 13.7 | 1749% | 62.7 | 4.2 | 169% | PASS |
| cluster · NUMA | 14.72× | 120.0 | 16.1 | 1803% | 63.6 | 4.2 | 168% | PASS |
| cluster · GPU | 14.72× | **117.5** | 15.0 | 1719% | 53.7 | 4.3 | 183% | PASS |
| cluster · GPU+NUMA | 14.72× | 119.2 | 15.3 | 1770% | **53.2** | 4.3 | 183% | PASS |
| ref single · CPU | 22.16× | 26.7 | 4.5 | 2096% | 18.9 | 1.9 | 689% | OK |
| ref single · NUMA | 22.11× | **18.2** | 4.2 | 2667% | 11.0 | 1.7 | 897% | OK |
| ref single · GPU | 22.16× | 22.5 | 4.2 | 2422% | 16.3 | 2.2 | 795% | OK |
| ref single · GPU+NUMA | 22.11× | 26.2 | 2.9 | 1795% | **10.6** | 1.8 | 964% | OK |
| ref paired · CPU | 23.98× | 43.0 | 5.0 | 2349% | 39.7 | 2.2 | 633% | OK |
| ref paired · NUMA | 23.96× | 43.7 | 3.4 | 2135% | **21.5** | 2.2 | 880% | OK |
| ref paired · GPU | 23.98× | 48.3 | 5.4 | 2063% | 39.0 | 2.3 | 654% | OK |
| ref paired · GPU+NUMA | 23.96× | **39.3** | 4.4 | 2305% | 21.9 | 2.3 | 878% | OK |

**Takeaways:**
- **GPU's main benefit is `qz --ultra` compress** (big serial BWT blocks): 30M ultra-2 **195 → 114 s (1.71×)**, 5M ultra-2 35.6 → 21.5 (1.66×). At the **default 25 MB block with `-t 72` it's a wash** (30M 34.8 → 34.2): the CPU already parallelizes the many small BWTs across 72 cores while the single GPU serializes them. GPU also *lowers* compress RSS in ultra (BWT workspace moves to VRAM).
- **bz: GPU is a regression.** Compress 26.3 → **59.4 s (2.3× slower)** (CPU% 2035% → 715%): bz's many small BSC blocks parallelize on CPU but serialize on one GPU. NUMA is bz's lever (decode 1.31×).
- **`--cluster`: GPU does nothing.** Its sequence stream is zstd-long and quality is fqz (neither BSC); only a tiny strand-bit stream is BSC. NUMA is flat too.
- **reference: GPU does not help.** Compress is serialization-bound (the in-order coverage fold), so the GPU BWT speedup is masked; GPU+NUMA can even *regress* compress (GPU contention across NUMA workers). NUMA is the effective lever for reference mode (decode 1.7–1.8×).
- **GPU × NUMA:** they coexist (byte-identical, no VRAM-OOM) but are **sub-additive**: one GPU serializes across the NUMA workers. They stack best on **ultra decode** (30M ultra-2 16.5 s, 1.96×).
- **Bottom line:** **NUMA helps across modes; GPU helps mainly `--ultra` compress**, and is neutral-to-harmful elsewhere.

## How it works

**Stream decomposition.** Each FASTQ record is split into three streams (headers, DNA sequences, quality scores), compressed independently and in parallel. Separating them lets each codec see statistically homogeneous data.

**Sequences & headers (libbsc).** DNA from one organism shares extensive k-mer content; the Burrows–Wheeler Transform clusters those shared substrings into long runs, which Quantized Local Frequency Coding then entropy-codes. This reaches ~1.85 bits/base on raw DNA without any read reordering. Headers compress similarly well as raw text through the same BWT pipeline.

**Quality (fqzcomp).** Quality scores are modeled by [fqzcomp](https://github.com/jkbonfield/fqzcomp), conditioning on read position, previous quality value, a stability flag, and local sequence context. This captures positional and sequential correlations in Illumina quality profiles that block-sorting alone cannot. fqzcomp is the quality backend across all modes (single-end, paired, and reference) in both qz and bz.

**Archive format (v5).** A single chunk-major container serves every mode: single-end, paired, reference (single and paired), `--cluster`, and `--ultra` all share it, with an `archive_type` byte in the header selecting the layout. It has four parts: a **front header** (version, archive type, codec flags, sequence/quality metadata); **interleaved per-block frames**, each `[block_len][record_count][crc32][payload]`, so every block carries its own CRC32; a **directory footer** mapping every block to its stream role (headers / sequence / quality / …), mate, and chunk; and a fixed **20-byte trailing locator** (`footer_len | footer_crc32 | "QZFOOTR1"`) pointing straight at the footer so the decoder never scans. The locator at EOF is also why the encoder never seeks backward, so compression can stream to a pipe. On decode the footer is CRC-checked and structurally validated (locator magic, per-entry spans, no block overlaps, role/chunk pairing) before any payload is read.

**Streaming engine.** Compression runs a producer → worker-pool → ordered-writer pipeline: the producer frames records into raw per-role blocks, workers compress them, and a single writer reorders and emits them, streaming to an atomic temp file that is renamed on success (stdout output bypasses the temp and writes directly). Worker count follows `-t` under a memory budget. Decompression mirrors this with a bounded streaming cursor over the footer index. The result is peak RAM that is constant in read count.

**Reference-based mode.** A vendored, map-only [strobealign](https://github.com/ksahlin/strobealign) aligns reads to the reference; qz then stores per read only its mapped position, strand, and edits versus a local consensus, with unmapped reads kept verbatim. Decoding rebuilds each read from the stored edits and the reference is not required.

**NUMA sharding.** On multi-socket hardware, `--numa` re-execs one worker process per socket, each pinned with CPU affinity and node-local memory policy before the thread pool starts, so its buffers fault locally. Each worker handles a disjoint chunk range and the driver assembles the parts; output is byte-identical to a single-process run.

**bz consensus-delta.** For aligned BAM, bz computes a local consensus over each block of reads and stores each read's sequence as its XOR against that consensus (typically >99% zeros, which the BWT pipeline compresses to almost nothing), alongside columnarized alignment fields and fqzcomp-coded quality.

**bz archive format.** `bz` uses its own container (magic `BZ`, version 13), distinct from qz's v5: a global header (flags, record and chunk counts, the compressed SAM header) followed by a forward sequence of chunks. Each chunk has its own header (record count, per-stream compressed sizes, and **two CRC32s**), then the columnarized streams (positions, flags, MAPQ, mate fields, CIGAR, the consensus and sequence diffs, quality, tags, and the mate/MD/NM derivation bitmaps). The first CRC covers the compressed payloads (on-disk corruption, checked before decompression); the second covers the original BAM records (round-trip fidelity, recomputed and compared at decompress time). Stream sizes are inline, so decode walks chunk by chunk; there is no trailing footer or locator.

## Troubleshooting

| Symptom                                                | Cause and fix                                                                                                                                                          |
| ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Build fails: `omp.h` not found, or `-lgomp` link error | No OpenMP toolchain. Install a C/C++ compiler with OpenMP (e.g. GCC + `libgomp`); `build.rs` compiles libbsc and htscodecs with OpenMP.                                 |
| `qz compress --reference` errors asking for an index   | No usable `.sti` index for that reference + read length. Build one with `qz index ref.fa --like reads.fastq` (or `-r <LEN>`), or pass `--build-index` to build inline. |
| Reference results look wrong after editing the FASTA   | Stale index cache. Delete `ref.fa.qz-r*.sti` and rebuild with `qz index`.                                                                                              |
| `bz compress` rejects the BAM                          | bz requires a **coordinate-sorted** BAM (`@HD … SO:coordinate`). Sort first: `samtools sort -o sorted.bam aln.bam`.                                                     |
| "output exists" / refuses to overwrite                 | Pass `-f` to overwrite the existing output.                                                                                                                            |
| "unsupported archive version" / legacy archive rejected | Pre-release format change: the archive predates this build. Recompress from the source FASTQ/BAM (see the format-stability note under [Status](#status--contributing)). |

## Project layout

A Cargo workspace of eight crates:

```
qz/
├── crates/
│   ├── qz-lib/        FASTQ compression library (all algorithms; builds libbsc + htscodecs)
│   ├── qz-cli/        the `qz` binary (Clap CLI)
│   ├── qz-python/     PyO3/maturin bindings (the `qz` Python package)
│   ├── qz-bench/      development benchmark binaries (not shipped)
│   ├── bz-lib/        BAM compression library
│   ├── bz-cli/        the `bz` binary (Clap CLI)
│   ├── strobealign/   vendored, map-only strobealign (reference mode)
│   └── numa-core/     shared NUMA topology / pinning / sharding (qz + bz)
├── third_party/
│   ├── libbsc/        vendored in-tree (BWT + LZP + QLFC; compiled by build.rs)
│   └── htscodecs/     vendored minimal subset (fqzcomp quality codec, pin 69185ce)
└── benchmarks/        benchmark scripts
```

## Vendored dependencies

- **libbsc** (Apache-2.0, Ilya Grebnov). Vendored in-tree under `third_party/libbsc`. It carries a local security patch that bounds BSC decode writes to the output-buffer capacity; provenance and the patch are documented in [`third_party/libbsc/VENDORED.md`](third_party/libbsc/VENDORED.md). Bundles **libsais** (Apache-2.0) for suffix-array construction.
- **htscodecs** (BSD-3-Clause, Genome Research Limited / samtools). A vendored **minimal subset** under `third_party/htscodecs` (pin `69185ce`): only `fqzcomp_qual.c` + `utils.c` and the headers they include, one of which (`c_range_coder.h`) is public domain, compiled with `-DNO_THREADS`. Provenance and rationale in [`third_party/htscodecs/VENDORED.md`](third_party/htscodecs/VENDORED.md).
- **strobealign** (MIT, Kristoffer Sahlin). A vendored, library-only Rust port used by reference mode; provenance in [`crates/strobealign/VENDORED.md`](crates/strobealign/VENDORED.md).

## Status & contributing

Pre-release. The codebase carries an extensive round-trip integration suite (240+ qz and 28 bz integration tests, plus unit tests):

```bash
cargo +nightly test --release
```

**The archive format is not yet stable.** There are no released archives in the field, so the on-disk format may change between versions, and a newer build can reject an older `.qz` / `.bz` archive (with a recompress hint). Keep the binary you used to write any archive you need to read later, or plan to recompress from the source FASTQ/BAM.

Bug reports are very welcome; the most useful ones include a failing archive and the output of re-running the command with `--debug`.

## Acknowledgments

- [libbsc](https://github.com/IlyaGrebnov/libbsc) and [libsais](https://github.com/IlyaGrebnov/libsais) by Ilya Grebnov: block-sorting compression (BWT + LZP + QLFC) and suffix-array construction.
- [fqzcomp](https://github.com/jkbonfield/fqzcomp) by James Bonfield, via [htscodecs](https://github.com/samtools/htscodecs): context-modeled quality compression (the quality codec in CRAM 3.1).
  Bonfield, J.K. & Mahoney, M.V. (2013). Compression of FASTQ and SAM format sequencing data. *PLoS ONE* 8(3):e59190.
- [strobealign](https://github.com/ksahlin/strobealign) by Kristoffer Sahlin: fast short-read mapping (reference mode).
- [SPRING](https://github.com/shubhamchandak94/SPRING) and CRAM/[samtools](https://github.com/samtools/samtools): comparison baselines.

## License

qz and bz are free and open source under the permissive [Apache-2.0](LICENSE) license, with no restrictions on academic, commercial, or clinical use. The bundled dependencies are likewise permissively licensed: libbsc and libsais (Apache-2.0), htscodecs (BSD-3-Clause, plus one public-domain header), and strobealign (MIT). The Rust crate dependencies are all permissive (MIT / Apache-2.0 / BSD and similar) with no copyleft. Full attribution for every third-party component is in [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).
