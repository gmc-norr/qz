# Changelog

## Unreleased

### Docs: README usage rework + third-party license file

Reworked the README for usability and publish-readiness: the Quickstart now leads with
the reference-based workflow (the most common resequencing path); "Which mode should I
use?" became a "start with what you have, then tune" decision guide with a rule of thumb;
and the reference-based section was expanded (when to use it, unmapped reads kept verbatim
so output stays lossless, reference never stored/needed to decode, indexing + tuning).
Added a top-level `THIRD_PARTY_LICENSES.md` (all vendored components + crate deps, all
permissive, no copyleft), linked from the README, and a license-forward intro.

Post-review fixes (a multi-agent review of the adopted publish-prep batch found only
low/info doc nits — the code compiles and all Python tests pass): removed now-dead `docs/`
references (the design docs moved to a local, git-ignored area); scoped the paired-end row
so it no longer implies `--ultra` composes with paired input (it is single-end only and is
silently ignored on a two-`-i` job); corrected the integration-test count; clarified the
CRAM 3.1 ratio (1.83× default vs up to 1.94× with non-default codecs); and reconciled the
`THIRD_PARTY_LICENSES.md` crate table to the stated 205 (four permissive crates were
omitted).

### Feature: paired (R1/R2) I/O in the Python bindings

`qz.compress` and `qz.decompress` now handle paired-end data, so paired compression
and decompression no longer require the CLI:

```python
qz.compress(["R1.fastq", "R2.fastq"], "sample.qz")    # paired-end compress
qz.decompress("sample.qz", "sample", split=True)      # → sample_R1.fastq + sample_R2.fastq
qz.decompress("sample.qz", "out.fastq", interleaved=True)   # one interleaved stream
```

- **`compress`**: the `input` argument now accepts a single FASTQ path (single-end,
  as before) **or** a list `[R1, R2]` for paired-end. A one-element list is
  single-end. The library validates the count (1 or 2) and dispatches paired by
  input length; this is purely a binding change.
- **`decompress`**: new `split` flag. For a paired or paired-reference archive,
  `split=True` writes the two mates to separate files using `output` as a prefix
  (`<output>_R1.fastq` / `<output>_R2.fastq`, the same convention as
  `qz decompress -o <prefix>`). It is mutually exclusive with `interleaved`, and is
  rejected for single-output archives (with a hint). A two-mate archive now needs
  either `split=True` or `interleaved=True`; with neither it still errors.
- Tests: `crates/qz-python/tests/test_paired.py` (paired split + interleaved
  roundtrips, single-end via one-element list, and the split/interleaved error
  paths).

### Feature: `qz.merge` Python binding

The `qz merge` archive-merge primitive is now exposed in the Python package, so the
Python API covers `compress` / `decompress` / `verify` / `merge`. It stitches several
QZ archives into one without re-encoding (block frames relocated verbatim, directory
rebuilt), with the reads in input order:

```python
qz.merge(["a.qz", "b.qz"], "merged.qz")              # single-end or paired
qz.merge(shards, "merged.qz", reference="ref.fa")    # reference archives (type 2/4)
```

- **Signature:** `merge(inputs: list[str], output: str, *, reference: str | None = None, force: bool = False)`.
  `reference` is the FASTA (required to merge reference archives, whose coverage globals
  are re-derived rather than concatenated); `force` overwrites an existing output.
- All inputs must share one archive type and config; cluster archives cannot be merged.
  These constraints are enforced in the library and surface as `RuntimeError` (the
  binding adds no new validation). Releases the GIL during the work.
- Tests: `crates/qz-python/tests/test_merge.py` (roundtrip — compress two halves, merge,
  verify, decode to the byte-exact concatenation — plus force-overwrite and empty-input
  error paths).

### Benchmarks: corrected CRAM comparison (was understated by a reference mismatch)

The published chr20 `bz` vs CRAM numbers reported CRAM at **1.25–1.29×**, which was a
benchmark-harness bug, not a CRAM property. The chr20 BAM was aligned to
`GRCh38_full_plus_hs38d1` (2580 `@SQ` contigs incl. decoys) and has no per-contig M5
tags, but `bench_bam.sh` fed CRAM the decoy-less `no_alt` reference. samtools could not
resolve every `@SQ`, silently enabled `embed_ref=2` (embedding the reference into the
CRAM), and the file bloated. With a contig-matching reference (verified 2580/2580), CRAM
is **1.77× (3.0) / 1.83× (3.1 default) / 1.94× (3.1 with `use_fqz`/`use_tok`/`use_arith`)**.

- `bz` still wins on ratio (2.00×) and needs no reference, but the honest margin over CRAM
  is ~3–9%, not ~55%. CRAM decodes ~3× faster.
- **Harness fix:** `bench_bam.sh` now uses a per-BAM reference (`REF_FULL`/`REF_CHR20`),
  **aborts any CRAM row that falls back to `embed_ref`** so a mismatch can never silently
  publish a bogus ratio again, adds a CRAM 3.1 best-codecs row + a pigz baseline, and gates
  the 51 GB WGS tier behind `DO_FULL=1`.
- The BAM is genuinely coordinate-sorted (0 out-of-order POS over 18.3M records); the prior
  "`SO:unsorted` depresses CRAM" explanation was wrong and has been removed. README,
  `SUMMARY-2026-06-24.md`, and `CLAUDE.md` updated.

### Feature: `qz decompress --interleaved` — paired archives as one interleaved FASTQ

Paired (`archive_type` 1) and paired-reference (`archive_type` 2) archives can now be decoded
as a **single interleaved FASTQ stream** (R1/R2 records alternating: `r0/1, r0/2, r1/1, …`)
instead of separate `_R1`/`_R2` files. This is the streamable paired form aligners consume
directly — `bwa mem -p`, `bowtie2 --interleaved`, `strobealign --interleaved`,
`fastp --interleaved_in`, BBTools `int=t` — so a paired archive can pipe straight into an
aligner:

```
qz decompress -i sample.qz --interleaved -o - | bwa mem -p ref.fa - | samtools sort -o sample.bam
```

- **CLI:** `--interleaved` on `qz decompress`. The single output goes to a file, `-`/stdout, or
  (`--gzipped`) a gzip stream. Runs in-process (does not combine with `--numa` sharding). Only
  valid for two-mate archives; single-end / cluster / single-end-reference archives are rejected
  with a clear message.
- **Python:** `qz.decompress(input, output, interleaved=True)` — the only way to decode a
  two-mate archive through the single-output Python API (the paired/reference rejection now
  points at `interleaved=True`).
- **Internals:** both decoders gained an output-mode enum (`PairedSink` / `RefSink`) threaded
  through the existing streaming writers; the interleaved branch merges each record-aligned
  R1[k]/R2[k] pair into one `format_batch_parallel` batch and writes a single sink — bounded
  memory (one batch/chunk resident), reusing the same atomic-temp / stdout / gzip output
  dispatch as single-end. The two-file decode paths are byte-for-byte unchanged (full paired +
  reference suites pass). Decode-only this increment; an interleaved *compress* input path is
  not yet wired.

### Python: `qz.verify()` binding + qz-python package README

The Python extension now exposes **`qz.verify(input, fast=False, threads=0, working_dir=".")`**,
wrapping the library's `verify` the same way `compress`/`decompress` are wrapped (GIL released
for the work, panics/errors converted to a catchable `RuntimeError`). It returns a dict of stats
— `ok`, `mode` (`"deep"`/`"fast"`), `num_reads`, `encoding`/`encoding_type`, the header/quality
compressor names, per-stream compressed sizes, `crc32` (deep) / `blocks_verified` (fast),
`r2_crc32` (per-mate for two-mate deep verify, else `None`), `total_bytes`, and `elapsed_secs`
— and **raises `RuntimeError`** if the archive is corrupt, truncated, tampered with, an
unsupported/legacy format, or (for `--cluster`) fails its multiset checksum. Unlike the
single-output `decompress` binding, `verify` works for **every** archive type — single-end,
paired, reference, and cluster — because it never reconstructs output to disk. `run_guarded`
was generalized to forward the `Ok` value (the existing `()` compress/decompress callers are
unchanged). Also added the missing `crates/qz-python/README.md` that `pyproject.toml`'s
`readme = "README.md"` points to, so the documented `maturin develop -m crates/qz-python/Cargo.toml`
build no longer fails with "Failed to read readme".

### Docs: README accuracy pass — lead with reference-based mode; fix stale bz ratio

Brought `README.md` in line with the 2026-06-24 benchmark ground truth and the intended
framing. The headline now leads with **reference-based mode** (qz's strongest lossless result,
~22–24×) rather than the shallow single-end ~8.2× WGS figure (promoted/quantified the
reference Highlights bullet and added a reference-first orienting line to the Performance
intro). Corrected the bz capsule ratio from a stale **~1.8×** to the measured **~2.0×** and
dropped the wrong "250 bp WGS" descriptor (the BAM is HG002 chr20 2×151). Scoped the
"beats SPRING on ratio" claim to reference mode (on plain FASTQ, SPRING's read-reordering
wins). Documented the shipped runtime `--gpu auto|off|require` flag (qz default `auto`, bz
default `off`) and reframed the experimental-GPU blurb to match the measured takeaways.
Disclosed the bz BAM `SO:unsorted`→reheadered caveat, and fixed the integration-test count
(236 → 237 qz). Added `pigz -9` general-purpose baselines to the paired (5.12×), reference
single/paired (5.23× / 5.12×), and BAM (1.00× — re-gzip of already-BGZF data) Performance
tables so every dataset has a generic-gzip floor for comparison (script
`benchmarks/bench_pigz_all.sh`; all round-trip `OK`).

### GPU: `--gpu auto|off|require` flag (cuda builds only) + `QZ_GPU` runtime toggle

A `--features cuda` build now exposes `--gpu` on `qz compress`/`decompress` (and `bz`'s
`compress`/`decompress`/`extract`): `auto` uses the GPU BWT with per-block CPU fallback, `off`
forces the CPU BWT, `require` errors at startup if no usable CUDA device is present. The flag
is **absent from the default CPU build** (`#[cfg(feature = "cuda")]`), so `--help` only
advertises GPU control on a binary that can actually do it, and passing `--gpu` to a CPU build
is a clean "unexpected argument" error. Backed by a `QZ_GPU` env var (`off`/`0`/`false` = CPU
BWT) which the flag sets — so the choice propagates to `--numa` worker *processes* through the
inherited environment. `--gpu auto` (the default) *respects* an externally-exported `QZ_GPU`
(it leaves the env untouched); `--gpu off`/`require` override it. Output is **byte-identical**
across `--gpu off`, `--gpu auto`, and a plain CPU build (verified by md5: GPU only selects the
BWT backend, never the archive bytes).

`qz` defaults to `auto` (GPU on, preserving prior cuda-build behavior). **`bz` defaults to
`off`** — GPU regresses bz (its many small BSC blocks serialize through the single GPU, while
the CPU path parallelizes them across cores), and bz is quality-dominated (fqz, which the GPU
cannot touch), so there is no block-size retune that would make GPU pay off; `--gpu auto`/
`require` opt back in. This is orthogonal to the separate "a cuda binary won't load without the
CUDA runtime" property — the CPU build remains the portable default.

### Reference compress: `QZ_REF_TIMING` phase-split diagnostic + two deferred levers resolved DUD

`QZ_REF_TIMING=1` now prints a per-phase wall-time breakdown of the reference-compress producer
to stderr — `map_and_diff` (strobealign mapping + Hamming diff), `cov_fold` (the serial in-order
coverage fold), `send_blocked` (worker-pool starvation), and `finalize` — for both single-end and
paired (`RefProdTiming` in `reference/stream.rs`, shared by `stream_single.rs`). Pure diagnostic:
timing is observational (encoded bytes unchanged) and the clock is read only per-chunk — never when
the env var is unset.

It settles a long-open question (`docs/reference-paired-hotpath.md` items P1/P2, now both marked
resolved-DUD): the reference-compress gate is **`map_and_diff` (~60–72% of wall)**, NOT the serial
coverage fold (measured **<1%**), and **`send_blocked = 0`** at every worker count — the pool is
hidden under mapping, so neither parallelizing the fold (P1) nor raising `--reference-window` (P2;
byte-identical and wall-flat across 2/4/8/16) buys anything. Mapping is bandwidth-bound (IPC 1.06);
the only levers that move reference compress are a ratio tradeoff (`--reference-fast`) or more memory
domains (`--numa auto`) — both already shipped. No code path or default changed.

### `bz extract`: BGZF reader workers now scale with cores (was hard-coded 4)

`bz extract` opened its BAM with a fixed 4 BGZF inflate workers — the lone reader that ignored
available cores (`bz compress` already uses `min(cores, 16)`). It now matches: `min(cores, 16)`.
Reader worker count never affects output; this is a marginal speedup on the niche FASTQ-from-BAM
extraction path and a consistency fix.

### GPU: `bz` can now build with `--features cuda`

`bz-cli`/`bz-lib` gained a `cuda` feature that passes through to `qz-lib/cuda`, so a
GPU build of `bz` routes its BSC BWT through `libcubwt` like `qz` does. Opt-in; the
default build is unchanged (CPU-only, portable). Benchmarking (see
`benchmarks/results/SUMMARY-2026-06-24.md`) found GPU **regresses** bz compress
(~2.3× slower — bz's many small BSC blocks parallelize across cores on CPU but
serialize through the single GPU), so the feature exists for parity/experiments,
not as a recommended bz config. New benchmark drivers `benchmarks/bench_gpu.sh`
(single-end CPU/NUMA/GPU/GPU+NUMA matrix) and `benchmarks/bench_gpu_modes.sh`
(bz / cluster / reference).

### Review fixes: input-validation, DoS-hardening, and doc accuracy

Findings from a multi-agent code review of the workspace, each fixed with a
regression test where the behavior changed:

**Critical / High — silent corruption & crash-on-valid-input:**
- **Cluster mode now rejects a FASTQ record whose sequence and quality lengths
  differ** (`cluster/sort.rs`), mirroring the standard reader. Cluster decode
  reconstructs each sequence boundary from its quality length, so the previous
  acceptance produced silent corruption or an undecodable archive. The unused
  `seq_lens` accumulator field was removed (the invariant it documented is now enforced).
- **Open-syncmer length guards corrected** (`dna_utils.rs`): a read of *exactly*
  k=21 bp passed the `< k` guard and tripped `find_syncmers_pos`'s `assert!(len > k)`,
  panicking a cluster-mode worker on adapter-trimmed / fixed-21-mer libraries.
  Guards are now `<= k`.
- **Empty-quality records now fail with a clear, actionable error** instead of the
  opaque "fqz_compress returned NULL" (`codecs.rs`): fqzcomp cannot encode a
  zero-length quality string; the message names the cause and the BSC-quality
  workaround. (Full lossless fqz support for empty reads is a separate follow-up.)
- **Header columnar encoders guard the u16 `write_string` length** (`header_col.rs`):
  an SRA name-prefix / pair-suffix or Casava common-comment longer than 64 KiB
  previously truncated its stored length while writing the full bytes, producing a
  silently-undecodable archive. Such inputs now fall back to the byte-exact raw path.
- **Reference merge validates the supplied FASTA against each shard's
  `ReferenceMeta` digest** (`reference/merge.rs`): merging shards against the wrong
  reference produced a corrupt archive that still passed `verify`. This also enforces
  cross-shard reference consistency.

**Medium:**
- **bz consensus builder caps the per-read reference span** (`bz-lib/streams.rs`):
  a crafted CIGAR with huge D/N ops could drive a multi-GB allocation during compress
  (decode already had this guard). Over-span reads are skipped from the consensus
  predictor (still encoded losslessly).
- **Paired `verify --deep` now reconstructs through the same dispatch as production
  decode** (`paired/mod.rs`) — the streaming engine for uniform archives (the common
  case), the legacy serial path only for mixed R2-header codecs — so deep verify
  validates the decoder users actually run. Each mate gets its own CRC sink.
- **strobealign `.sti` index reader validates `main_hash_mask`** before deriving the
  orientation bit (`strobealign/index.rs`): a corrupt/foreign mask no longer
  underflows / overflow-shifts; it is rejected as a parameter mismatch.

**Low (correctness / hardening / cleanup):**
- Reference position decode uses checked `u32::try_from` for `ref_id` (`positions.rs`).
- strobealign `is_too_frequent_forward` bounds-checks with the index parameter it
  actually uses (`cutoff`), closing a latent OOB.
- bz `derive_tlen` computes in i64 and clamps to the i32 SAM field, avoiding overflow
  on an untrusted CIGAR (`md_codec.rs`).
- bz fqz multiblock decode rejects trailing bytes after the last block (`decompress.rs`).
- `--ultra` always runs one block at a time even under `QZ_COMPRESS_WORKERS` (the env
  override previously bypassed the clamp, risking concurrent big-block BWT OOM).
- Bounded decode (`decompress_parallel_max`) on the header raw-fallback and cluster
  metadata/strand globals, closing stacked-BSC-bomb vectors on untrusted input.
- Single-end compress clamps `chunk_records` to >= 1 (library-API robustness).
- qz-python surfaces errors via anyhow's clean cause chain (`{e:#}`) instead of `{:?}`.
- Removed dead utilities (`pack_dna_2bit` / `unpack_dna_2bit` / `shannon_entropy` /
  `hamming_distance_within`) and demoted the reference `positions` / `refmeta`
  submodules (no external consumers); slimmed `ReferenceMeta` to its used fields.
- Documentation accuracy: CLAUDE.md crate count (8, incl. `numa-core`), the live
  `ClusterReorder` (encoding 11), and the `integration_test.rs` test count; corrected
  the `codec_ids.rs` Fqzcomp-code comment, the `ultra.rs` peak-memory model comment,
  the `dna_utils.rs` module doc, and the `qz merge` subcommand help text.

One reviewed finding was **declined**: broadening NUMA-decode `auto` to fall back
in-process on *any* worker failure conflicts with the deliberate, tested design
(`worker_failure_kills_survivor_and_aborts_cleanly`) that aborts loudly on an
unexpected failure rather than masking it behind a slow in-process retry.

### Review fixes (round 2): remaining items + performance verification

Performance-relevant cleanups (verified byte-identical output, measured on HG002
chr20 / ERR3239334):
- **bz consensus aux: `nibbles_to_ascii` computed once per record** (lazy cache),
  not twice (MD + NM), and never for a record with no derivable MD/NM tag
  (`bz-lib/streams.rs`). Measured **~8% faster** bz compress on HG002 chr20
  (40.7 s → 37.5 s, `-t 16 --numa off`), byte-identical archive.
- **Cluster paired deep-verify holds only one mate per chunk** and streams the
  other, pairing positionally (`cluster/decode.rs`) — roughly halves the per-chunk
  verify transient, order-identical checksum.
- The `bz compress --numa auto` post-prescan in-process fallback double-inflate is
  now documented as a known, bounded cost (proper fix = non-inflating BGZF prescan,
  deferred); the sub-threshold and no-topology paths already skip the prescan.

Additional hardening:
- Paired legacy/verify decode (`emit_records`) now enforces per-record quality ==
  sequence length and column-count agreement, matching the streaming path.
- Streaming paired decode rejects an archive whose R1 and R2 quality codecs differ
  (the single `bits_per_qual` assumes uniformity; not encoder-reachable).
- Single-end NUMA layout query runs the cheap structural directory validator for
  parity with paired/reference (rejects a forged footer before the planner uses it).
- Single-end cluster decode rejects extra `-o` arguments instead of silently
  dropping them (parity with the paired-cluster and dispatcher checks).

Verification: full workspace `cargo build`/`clippy` clean (0 warnings); the qz and
bz compress/decompress hot paths produce byte-identical output before/after (the qz
5M-read archive and the bz HG002 chr20 archive both hash-identical to the
pre-change build), so no compression-ratio or correctness change.

Intentionally deferred (lower value / larger refactor / superseded): the
bsc/dna_utils `pub`→`pub(crate)` demotion (largely moot — `qz-bench` is a
legitimate cross-crate consumer of those helpers); the `reference/stream_single`
producer dedup (pure maintainability); the bz-verify derivation-stream (15–24)
walk; a `renameat2(RENAME_NOREPLACE)` for the output TOCTOU window; and several
coverage-gap tests for unchanged code.

### Repo cleanup: removed dead PGO hooks + stale files

- **Removed the C/C++ PGO build hooks** from `crates/qz-lib/build.rs` — the
  `QZ_PGO_GENERATE` / `QZ_PGO_USE` env vars are gone. PGO was measured to *regress*
  compress ~6% (entropy-coder branches are unpredictable, so PGO de-optimized the
  hand-tuned hot loops) and the hooks never linked correctly, so they were dead
  weight. Removing them is build-output-neutral.
- **Deleted stale/dead files**: research prototypes for the removed context-modeled
  quality codec (`qual_ctx_*.py`, `quality_context_compress.py`,
  `quality_model_convergence.py`), the `QZ_REORDER` benchmark kit (encoding 8/9, also
  removed), and assorted junk (a stray 2.6 MB blob, a `rustc-ice` crash dump, scan
  outputs, `__pycache__`). Added `.gitignore` rules so they don't recur.

### libbsc vendored in-tree (build no longer clones it)

`third_party/libbsc` is now a **tracked, vendored copy** of upstream
(`IlyaGrebnov/libbsc` @ `baffa62`, v3.3.12) instead of a gitignored external clone.

- **Why:** the QZ-SECURITY-PATCH (`bsc-decode-output-bound`, the heap-overflow
  hardening that bounds BSC decode writes to the output-buffer capacity) previously
  lived only as uncommitted edits in a per-machine clone — a single point of failure.
  It is now in qz's history, captured verbatim in `third_party/libbsc/QZ-SECURITY-PATCH.diff`
  and documented in `third_party/libbsc/VENDORED.md`.
- **Build:** no `git clone`/submodule step for libbsc anymore — the source builds
  straight from the tree. The build.rs sentinel guard stays as defense-in-depth
  against a future re-vendor silently dropping the patch.
- **Unchanged:** `third_party/htscodecs` remains a gitignored upstream clone
  (pinned `69185ce`), and `crates/strobealign` was already vendored.

### `qz index` subcommand + require-reference-index-by-default

Reference-mode compression now **requires a usable strobealign index up front** instead
of silently building one inline on every run. Build it once with the new subcommand:

```
qz index ref.fa -r 150        # build the sidecar for the ~150 bp read profile
qz compress -i reads.fq -o out.qz --reference ref.fa   # reuses the prebuilt index
```

- **`qz index <FASTA> [-r LEN] [--like FASTQ] [--reference-fast] [-t N]`** builds the
  profile-canonical sidecar (`ref.fa.qz-r<LEN>.sti`) next to the FASTA. `-r` sets the
  read-length profile (or `--like` peeks it from a FASTQ's first read).
- **Compress now errors by default when no usable index exists** (missing, stale, corrupt,
  or built with incompatible parameters), printing an actionable hint that names the exact
  `qz index …` command to run. Corruption is caught by the full index load at compress
  time (not just the cheap up-front header check), so a present-but-truncated/damaged
  sidecar errors rather than silently rebuilding. This makes the index-build cost explicit and one-time
  rather than a hidden per-compress tax, and lets a single prebuilt index serve repeated
  compress runs and all NUMA shards.
- **`qz compress --reference … --build-index`** opts back into inline building (the old
  implicit behavior, now explicit) — the parent process builds the index before
  compressing or sharding.
- **NUMA-safe:** the require/build check runs in the parent for the off / in-process path
  AND once per distinct read-length profile across all shards before any worker spawns, so
  workers stay load-only and never race to write the sidecar. `--build-index` is not
  forwarded to workers.
- The **policy is CLI-only**: the library's own `compress()` still auto-builds for library
  and test consumers; only the `qz` binary enforces require-by-default.
- **Reference index sidecars are now keyed by the seeding *profile* (canonical read
  length)** rather than the raw first-read length, so read lengths in the same bucket
  (e.g. 143 and 150 bp) share one cached index. Old `.qz-r{rawlen}.sti` sidecars that
  don't match the canonical name are simply rebuilt once.

### bz: NUMA-aware compress sharding (`bz compress --numa`, ~1.44× at chr20 scale)

`bz compress` gains a `--numa auto|off|N` flag mirroring the decode side. On NUMA
hardware, `auto` (default) re-execs one pinned worker process per socket; each
compresses a **disjoint chunk range** of the input BAM into a self-contained `.bz`
part, and the driver concatenates the parts into the final archive — **byte-identical
to a single-process `bz compress`**. This breaks compress's per-process ~16-core
futex ceiling (two node-pinned pipelines run independently). Off NUMA, sub-threshold
inputs, or `--numa off` compress in-process with zero overhead. Built on the same
shared `numa-core` crate as the decode path.

- **Measured ~1.44× at chr20** (HG002 chr20, 18.4M records, `-l 1`, `-t 72`,
  interleaved A/B on a loaded box: off ~85.1s → numa-2 ~59.1s; range 1.28–1.61×),
  byte-identical archive. The win is the parallel compress minus a serial prescan tax.
- **The input split needs a one-time prescan** (the only structural difference from
  decode, which got seekability for free via `stream_sizes`). A BGZF BAM has no record
  sync marker, so workers can't resync mid-stream: `read_bam_compress_layout` walks the
  BAM once (multithreaded inflate, no compression), capturing the BGZF **virtual
  position** at each chunk boundary. A worker seeks straight to its assigned start while
  keeping multithreaded inflate (the noodles MT reader is seekable). This prescan is the
  serial Amdahl component of the shard.
- **Pure-concat assembly, no header patching.** The driver knows the whole-archive
  `num_records`/`num_chunks` from the prescan, so worker 0 writes the global header with
  the final totals and later workers write headerless chunk streams — assembly is a byte
  concatenation. Chunk bytes are independent of where a worker starts (`compress_one_chunk`'s
  chunk index is log-only), so a range's chunks equal the single-process run's exactly.
- **Engine unification.** `compress()` and the range worker share one extracted
  `run_compress_pipeline` (the windowed chunk pipeline, now budget-bounded), so the
  shard can't drift from the single-process encoder. New bz-lib surface:
  `read_bam_compress_layout`, `compress_chunk_range`, `compress_peak_rss_bound`,
  `BzCompressLayout`. New bz-cli `_numa-compress-worker` subcommand + `--numa` on
  `compress` (the worker receives the driver's `AdvancedOptions` as JSON for exact
  byte-identity). A `bz-lib` byte-identity test (single-process == split-and-concat) and
  a `bz-cli` forced-shard CLI test (2 workers, byte-identical archive) guard it.
- **Hardening (multi-agent review).** A review caught that the **default `--numa auto`
  with the auto level** (`-l 0`) re-resolved the level — and thus `chunk_size` — in the
  prescan and in *each* re-exec'd worker independently, reading live `/proc/meminfo` at
  different times; near the L1/L2 RAM knee two reads could disagree, so a worker's
  chunking diverged from the prescan's virtual-offset partition (silently dropping
  records, or hard-failing). Fixed by **freezing the resolved level once** in the driver
  (`resolve_compress_level`) before the prescan and worker serialization, so every
  `apply_level` is idempotent and one `chunk_size` holds across the prescan and all
  workers. As defence in depth, each worker now also receives its prescan record count
  and **bails (no output) if it compresses a different number**, so any future drift is a
  hard error rather than silent loss. Two further fixes: the per-worker compress RSS gate
  bound budgeted only half the per-process base (now full — it is the gate's sole compress
  OOM guard, since the window is preset and never RAM-self-limited); and `--numa auto` now
  skips the full-inflate prescan for sub-threshold inputs (no double inflate). Regression
  tests added: an auto-level (`-l 0`) forced-shard round-trip, an `expected_records`-drift
  bail, and `compress_peak_rss_bound`/`resolve_compress_level` unit tests.

### bz: libdeflate backend for BGZF (BAM) I/O — ~5% faster decode, slightly smaller output, lossless

bz now compiles `noodles-bgzf` with its `libdeflate` feature, swapping the BGZF
deflate/inflate implementation from the pure-Rust miniz_oxide (flate2) backend to
libdeflate — the same C library htslib/samtools use. Measured (HG002 chr20, `-t 36`,
interleaved A/B on a loaded box):

- **Decode ~5% faster** (~39.6s → ~37.6s). The output BGZF deflate is ~24% of decode
  CPU but runs in the 16-thread writer pool that overlaps BSC/fqz decode, so
  libdeflate recovers only the exposed tail (not the full 24%). The output BAM is
  **byte-identical** under `samtools view --no-PG` and ~0.36% smaller (libdeflate's
  deflate is more efficient at the same level).
- **Compress unchanged** — the input-BAM inflate is fully hidden behind the BWT
  pipeline, so a faster inflater buys nothing; the `.bz` archive is **byte-identical**
  (inflate is a deterministic decode).
- Removes bz's deflate handicap vs samtools, making BAM-tool speed comparisons
  apples-to-apples.

Wiring: a direct `noodles-bgzf` dependency with `features = ["libdeflate"]` enables
the feature via Cargo unification onto the shared v0.45 instance (the umbrella crate
exposes no passthrough); `bam.rs` imports `noodles_bgzf` directly so the dep is used.
No new build burden — the workspace already compiles C deps (libbsc, htscodecs).

### bz: NUMA-aware decode sharding (`bz decompress --numa`, ~1.4× at chr20 scale)

`bz decompress` gains a `--numa auto|off|N` flag mirroring qz's NUMA decode. On
NUMA hardware, `auto` (default) detects ≥2 asymmetric sockets and, for inputs above
a size threshold, re-execs one pinned worker process per socket
(`sched_setaffinity` + `set_mempolicy(MPOL_BIND)` before rayon builds, so decode
buffers fault node-local). Each worker decodes a disjoint chunk range to a BGZF/BAM
part; the driver concatenates the parts into the final BAM. Off NUMA, sub-threshold
inputs, or `--numa off` decode in-process with zero overhead. This breaks bz
decode's per-process ~16-core futex ceiling: two node-pinned processes run their
pipelines independently (~2× the effective cores).

- **Shared `numa-core` crate.** The format-agnostic half of qz's NUMA machinery —
  topology discovery, CPU/memory pinning, the spawn/poll/cleanup runtime, the
  decision gate, and part-file assembly — was extracted from `qz-cli` into a new
  `numa-core` library crate that both CLIs depend on ("one tool, not a mesh"). qz's
  NUMA behavior is unchanged (its full suite — 22 `numa_cli` + 20 `numa_compress_cli`
  — stays green); only the import paths moved.
- **bz seekability without a format change.** A `ChunkHeader` already records its
  stream sizes, so a chunk's on-disk size is recoverable; `read_bz_chunk_layout`
  pre-scans chunk byte-offsets + read counts by reading each header and seeking past
  its data. `decode_chunk_range(start, end, …, out_part)` seeks to a range and
  decodes it into a self-contained BGZF part (the full `decompress` and the range
  decode now share one `run_decode_pipeline`).
- **BGZF assembly.** Each part ends in the 28-byte BGZF EOF marker; the assembler
  stream-copies the parts and strips that marker from every part but the last
  (verifying it first), yielding one continuous valid BAM. Part 0 carries the BAM
  header; later parts are records only. Output is **record-identical** to the
  in-process decode (the per-chunk content CRC validates each worker's reconstruction),
  though not byte-identical (BGZF reframes at chunk boundaries).
- **Measured ~1.40×** on HG002 chr20 (18.4M records, 13 chunks, `--numa 2` vs `off`,
  interleaved A/B on a loaded box: off 41.2s → numa2 29.4s, 1.28–1.62× range), with
  `numa2` also steadier; byte-identical-lossless (inflated payloads compare equal).
  The win grows at full-file scale (like qz's NUMA). Cluster archives and stdio are
  not shardable (decode in-process); BGZF's variable per-chunk size means bz always
  takes the part-file + assembly path (no direct-write).

### bz: parallel fqzcomp quality decode (~13% faster decompress, byte-identical)

bz decompress decoded the per-chunk fqzcomp quality **blocks serially** (`for b in
0..num_blocks { fqzcomp::decompress(...) }`) even though the compress side already
encodes them in parallel (`compress_quality_fqz_parallel`). Re-profiling after the
fqzcomp quality migration (env-gated `BZ_DECODE_TIMING` phase timers, matching the
`QZ_*_TIMING` convention) showed this serial loop was the single dominant decode cost
— **~59% of all decode work** on HG002 chr20 (18.4M records, 13 chunks), running at
only ~4× effective concurrency on 72 cores. (The old `finding_bz_speed_lever_pipeline_depth`
profile, which blamed the long-gone `quality_ctx` codec, was stale.)

- **Fix:** `decompress_quality_fqz_multiblock` now does a cheap serial framing walk
  (locate each block's byte range + starting read index — the byte offsets are
  sequential, so this can't be parallel, but it only reads 12-byte headers) followed
  by a `par_iter` over the independent blocks (each self-contained: own CRC, fqz blob,
  and disjoint output read-slice). All DoS guards (record-count vs remaining, block
  overflow, per-block CRC) are preserved, and the untrusted `num_blocks` no longer
  drives an unbounded capacity reservation.
- **Measured ~1.15× faster decode** (HG002 chr20, interleaved within-process A/B to
  cancel box load, `-t 72`: serial mean 35.96s → parallel 31.31s, −13%), **byte-identical
  output** + deep-verify pass. The fqz phase's CPU work dropped 141s → 52s.
- The remaining decode wall is bound by the ~16-core pipeline ceiling (BGZF output
  writer + in-order drain), not fqz — a separate follow-up. This fix still matters
  beyond the 13%: it removes 88s of fqz CPU work (shared-box throughput) and clears the
  serial bottleneck that would otherwise re-dominate once the per-process wall is broken
  (e.g. NUMA decode sharding).

### Reference NUMA-compress sharding: automatic `--numa` reference compress (1.72× at scale)

`--reference` compression now shards under `--numa` like single-end and paired — completing
the NUMA-compress matrix (every non-cluster, non-gzip FASTQ mode now shards). The input is
byte-range split, each node-pinned worker compresses its shard into a reference part archive,
and the parts are stitched by the reference-aware merge (see the entry below).

- **Byte-range reference compress.** New bounded `FileReader::Bounded` + `FastqReader::from_range`
  (plain-file seek + `take`, the parallel to single/paired's `file.take()`), threaded through the
  reference prelude so a worker compresses only the records in its `[start, end)` shard. The
  reference front header is fixed metadata (not data-derived), so every shard's header is
  byte-identical and the merge's front-header gate passes with no plan-pinning.
- **Driver + worker.** `compress_byte_range` gains a reference arm (1 range → `archive_type 4`,
  2 → `archive_type 2`); the compress worker re-exec forwards
  `--reference`/`--reference-window`/`--reference-fast`; `run_sharded_core` routes reference to
  `merge_reference_archives_to_path`. The gate's per-worker RSS bound adds the reference size
  (each worker loads the strobemer index).
- **Measured 1.72× at production scale** — HG002 chr20, full 10.8M single-end reference reads,
  real 2-node NUMA, `--numa 2` vs `--numa off`, min-of-2: **16.48s vs 28.30s**,
  byte-identical-lossless + deep-verify. At scale the in-process reference compress is
  bandwidth-bound and the two node-pinned workers break the bandwidth wall (the same mechanism
  as NUMA decode), so the win EXCEEDS the ~1.38× raw estimate. Archive +0.22% (the parallel-merge
  backing blocks). At small scale (~2M reads) it's a wash — the fixed serial tax (split +
  per-worker global re-derivation + merge) isn't amortized — so the gate only shards inputs above
  its workload threshold.
- Coverage: forced-shard single + paired reference compress roundtrip (`numa_compress_cli`); the
  stale "reference + fixed errors" precheck flipped to "reference + fixed shards losslessly."

### Reference archive merge (`qz merge --reference`) — unblocks reference NUMA-compress sharding

Reference archives (`--reference`, paired `archive_type 2` and single-end `archive_type 4`)
can now be **merged** into one via `qz merge --reference <fasta>` (and the library
`merge_reference_archives_to_path`). This was the missing piece for reference NUMA-compress
sharding: the per-shard coverage globals are not concatenable like single-end/paired frames.

- **Per-chunk frames relocate verbatim**, exactly like single-end/paired merge — Positions are
  stored absolute `(ref_id, ref_pos)`, and edits/strands/read-len/mapped-flags/FallbackPool are
  per-read functions of `(reference, reads)`, never of the coverage map.
- **Only the 3 coverage globals are re-derived**, over the UNION of the shards' covered
  intervals. They are reference-direct (PackedBacking/NBitmap are the raw normalized FASTA bases
  over the covered intervals; IntervalMap is that interval set), so re-deriving from
  `(reference, union-imap)` is deterministic; ReferenceMeta is FASTA-only. The FASTA is loaded
  once (`read_ref` + normalize, no index) and shared.
- **The re-derivation runs in PARALLEL** (rayon over independent ~8M-base BSC blocks) — and this
  is load-bearing, not a micro-opt. A serial re-derivation of the full chr20 backing measured
  **~3.6s**, which turned 2-process sharding into a **net loss (0.93×)**; the parallel version is
  **~0.75s**, making reference compress sharding **~1.27× at 2M reads** (real HG002 chr20:
  in-process 9.84s vs 7.0s workers + 0.75s merge), byte-identical-lossless. The smaller backing
  blocks cost ~+0.8% archive size at 2M reads — a fixed cost, so ~+0.15% at production 10.8M scale.
- Validated with the SAME `validate_reference_directory{,_single}` contracts decode runs, plus a
  rebuilt per-mate `ChunkDecodedSizes` table so the merged archive stays NUMA direct-write
  decodable. Coverage: single + paired reference merge roundtrip + deep/fast `verify`
  (`archive_merge` tests).
- **Automatic `--numa` reference *compress* sharding now uses this** (see the entry above) — the
  byte-range reference compress worker calls `merge_reference_archives_to_path` to stitch the
  shards; `qz merge --reference` remains the manual primitive.

### Reference NUMA decode: direct-write completes the matrix (1.27× faster, byte-identical)

Reference decode (`--reference`, both paired `archive_type 2` and single-end-reference
`archive_type 4`) was the last mode on the slower **part-file + concat** NUMA path. It now
**direct-writes** like single-end and paired: each worker seeks straight into its pre-sized
output(s) at the per-mate byte offset and writes its chunk-range — **no part files, no
concat pass**. With this, **every non-gzip mode shares one direct-write decode path**
(single = 1 output, paired/paired-ref = 2, single-ref = 1) via the same N-output
`decode_chunk_range(.., regions: &[DirectWriteRegion])` + `verify_direct_region` +
`BoundedCountingWriter` integrity contract.

- **The reference encoder now emits a `ChunkDecodedSizes` global** — `n_mates=2` for paired
  reference, `n_mates=1` for single-end reference. Reference is always lossless FASTQ
  (`reject_unsupported` forbids `--fasta`/`--no-quality`/lossy), so each read's decoded size
  is exactly `fastq_output_len(false, id, seq, qual)` with `qual == seq` — byte-exact to the
  decoder's output. The reference directory validators + deep/fast `verify` treat it as a
  cross-cutting global (the same `StreamRole` single/paired use), keeping the "exactly four
  reference globals" contract on `PackedBacking`/`IntervalMap`/`NBitmap`/`ReferenceMeta`.
- **Measured 1.27× faster than `--numa off`** (HG002 chr20, 1M-read single-end reference,
  real 2-node NUMA, nvme, min-of-2: direct **2.88s** vs in-process 3.65s), byte-identical.
  Reference decode is memory-bound so the win is smaller than paired's 1.40×, but it's a
  clear improvement (direct-write does strictly less I/O than part-file → cannot regress),
  it drops reference decode from 2× to 1× peak output disk, and it unifies the path.
- The driver's `direct_write` detection widens to `archive_type` `0 | 1 | 2 | 4` (the
  table-present check naturally excludes any archive without the global);
  `run_sharded_direct` pre-sizes 1 output for single-ref and 2 for paired-ref.
- **Reference *compress* is still in-process** (the compress-side asymmetry): it shards
  ~1.38× but needs a reference-aware merge over non-concatenable per-shard
  consensus/IntervalMap globals. That merge now exists (see the `qz merge --reference`
  entry above, ~1.27× net once parallelized); the remaining piece is the byte-range
  reference compress worker. This change is decode-only.
- Coverage: paired + single-end reference direct-write equivalence lib tests
  (`numa_chunk_range`), strengthened reference forced-shard CLI tests asserting the
  direct path (`numa_cli`), and the updated `read_chunk_layout` type-4 table assertion;
  all reference roundtrip + deep/fast verify suites green.

### Paired-end NUMA decode: direct-write into pre-sized R1/R2 outputs (1.40× faster, byte-identical)

Paired-end `--numa` decode previously used the **part-file** path: each worker decoded
its chunk-range to part files, then the driver did a serial concat pass into the final
R1/R2 outputs. Paired now uses **direct-write** (already in place for single-end): every
worker seeks straight into the two pre-sized shared outputs at its per-mate byte offset
and writes its region — **no part files, no concat pass**. Each worker re-derives and
verifies its `offset`/`len` against the authoritative per-mate `ChunkDecodedSizes` table
(a wrong table is a clean `DirectWriteIntegrityError` → `auto` falls back to in-process,
never silent corruption) and bounds its write to that length.

- **Measured 1.40× faster than the part-file path** it replaces (HG002 chr20, 3M pairs,
  real 2-node NUMA, nvme, min-of-3: direct **3.88s** vs part-file 5.45s), and **1.66×
  faster than in-process** `--numa off` (6.43s). Byte-identical R1/R2 across all three
  paths. The win is structural — direct-write eliminates a full read+write concat pass
  (the same mechanism that gave single-end its ~1.65×), so it cannot regress.
- **One path, not a fork.** The direct-write decode path is now generalized to **N
  outputs** (single-end = 1, paired = 2): a shared `verify_direct_region` +
  `BoundedCountingWriter` integrity contract, with `decode_chunk_range(..., regions:
  &[DirectWriteRegion])` carrying one region per output (empty ⇒ part-file path).
  `read_chunk_layout` now exposes `decoded_sizes_per_mate` (per-mate per-chunk sizes).
- **The paired streaming encoder now emits a `ChunkDecodedSizes(n_mates=2)` global**
  (row-major `sizes[chunk*2 + mate]`) — the table direct-write needs. Single-end output
  is byte-identical (its `n_mates=1` table is unchanged); paired adds the one global
  entry. (Supersedes the earlier "paired has no `ChunkDecodedSizes`" note.)
- **`qz merge` (paired) rebuilds the merged per-mate table**, so paired archives produced
  by the **default `qz compress --numa auto`** (which shards → byte-range-compresses →
  merges on NUMA hardware) also decode via direct-write — verified end-to-end: a
  `--numa 2`-compressed paired archive decodes in 3.87s (direct), byte-identical to
  `--numa off`. Previously the paired merge dropped the table, forcing the slower
  part-file decode.
- **Reference (type 2/4) is unchanged** — its columns aren't a flat self-delimiting
  stream, so it carries no `ChunkDecodedSizes` table and keeps the part-file path.
- Coverage: new paired direct-write equivalence + bad-region-integrity lib tests
  (`numa_chunk_range`), a paired forced-shard direct + region-mismatch-fallback CLI test
  (`numa_cli`), and a paired merge per-mate-table-rebuild test (`archive_merge`); single
  + paired streaming FNV goldens green (byte-identity oracle).

### Cluster single-end compress: parallel checksum-fold + alloc-free bucket write (~10% faster, byte-identical)

The single-end cluster sort's serial Pass-1 fold was the bottleneck — a sub-timer
(`QZ_CLUSTER_TIMING=1`) measured it at **71% of Pass-1 (~12s on chr20 R1 10.8M)**,
dominated by a serial double-FNV `hash_record` over every record's bytes plus three
per-record heap allocations. The verify checksum is an associative+commutative
`wrapping_add` of per-record hashes, so the hash now computes **inside the existing
parallel canon/key map** (across all cores, ~free given the memory-bandwidth headroom
there) and the fold just sums precomputed values; bucket entries are written straight
from the arena/canon slices via a new `write_entry_parts` (no intermediate
`ClusteredRecord` allocation). This brings single-end to parity with the paired cluster
path, which already computed its hash in parallel and wrote from slices.

- **Byte-identical output** (same bucket layout, same checksum add-order; verified
  old≡new at 206,729,826 bytes on chr20 R1, deterministic across runs, all 29 cluster
  unit + 5 integration tests green, `verify` OK).
- **Pass-1 sort 1.72× faster** (16.4s → 9.5s; the fold dropped 11.7s → 4.4s, par_iter
  grew only 2.41s → 2.52s). **Full single-end cluster compress wall ~+9.6%**
  (64.3s → 58.2s median, interleaved A/B, chr20 R1 `-t 72`).
- **`QZ_CLUSTER_TIMING=1`** now emits a Pass-0/1/2 + per-sub-phase breakdown for the
  single-end sort too (env-gated, zero production overhead), matching the paired path's
  existing timer.
- **Fixed: `--cluster` + `--numa auto` silently ignored clustering.** On NUMA hardware
  the default `--numa auto` split the input and compressed two **non-clustered** default
  halves (producing an ~8.5× archive instead of the cluster ratio). Clustering reorders
  reads globally and can't be byte-range sharded + merged, so it now forces in-process
  compress under `auto` (and errors loudly under fixed `--numa N`), mirroring the
  existing decode-side guard.

### Review fixes: output-contract guard, dead Qvz removal, clippy-clean, docs

Addresses findings from an external (Codex) code review. No archive-format or
ratio change; all modes byte-identical (reference/paired/cluster/single roundtrip
suites green).

- **Decompress output-contract guard (behavior).** Paired and (non-cluster)
  reference archives decompress to a single `-o` PREFIX (→ `_R1.fastq`/`_R2.fastq`);
  single-end, single-end-reference, and cluster-single take one exact path; only
  paired-cluster takes two exact paths. Passing the wrong number of `-o` values is
  now **rejected with a clear error** instead of silently dropping `output[1..]`
  (which previously produced misnamed files like `R1.fastq_R1.fastq`). The
  `qz decompress` `-o` help text now describes each mode's contract.
- **Removed dead `QualityMode::Qvz` (API).** The public `Qvz` variant reached an
  `unimplemented!()` panic on the single-end Auto path and was never constructible
  from the CLI or Python. Removed the variant and its (panicking / TODO) match arms.
  No archive-format impact — `quality_mode` is a compress-time knob, never serialized.
- **`cargo clippy --workspace --all-targets` is now clean** (the workspace's stated
  gate). Mechanical lint fixes only; behavior-preserving (verified byte-identical).
- **API docs:** `compress`/`decompress` (and the Python bindings) now document that
  `threads` initializes rayon's process-global pool on the first call — later calls
  with a different count are ignored (first-call-wins), which matters for long-lived
  services and Python sessions.
- **README** brought up to date with the v5 single-architecture format (was
  documenting the removed v4 stream-major layout), `--ultra` levels 1–3 (was 1–5;
  ultra's gain is larger BSC blocks, not read reordering), the paired/reference/
  cluster modes, and the paired/reference decompress output-prefix contract.

### Hardening: cluster-mode hostile-decode coverage + shared varint overflow guard

Robustness fixes surfaced by an adversarial project review; no format or ratio change,
all modes byte-identical.

- **`dna_utils::read_varint` shift-overflow guard.** The cluster/shared LEB128 reader
  could shift past 64 bits on a crafted over-long continuation sequence (UB-adjacent
  silent wrap). It now rejects `shift >= u64::BITS`, rejects any chunk that would
  overflow the accumulator, and caps at 10 continuation bytes — matching the three
  sibling readers. Returns `None` (the existing clean error path); never panics.
- **Cluster `Qual` allocation clamp.** `parse_fqz_varint_stream` now sizes its result
  vector with `record_count.min(data.len() + 1)`, so a forged footer `record_count`
  can't drive a multi-gigabyte pre-allocation before the frame is even read (matches
  the paired path).
- **Cluster inline-vs-footer `record_count` cross-check.** Both `read_frame` and the
  `verify --fast` cluster path now assert the per-block framed `record_count` equals
  the footer entry's `record_count`, so a tampered directory can't desync the two.
- **Cluster hostile-decode test suite.** New `cluster_hostile_test.rs` mutates a real
  on-disk cluster archive (unknown `archive_type`, footer-CRC tamper, two truncations,
  inflated `n_entries`, dropped `ClusterMeta`, oversized `ClusteredSeq` span,
  payload-flip, `record_count` mismatch) and asserts `decompress` + `verify` (fast +
  deep) all return `Err` with no leaked output and no panic — closing the gap the
  reference/paired suites already covered.
- **Clean rejection of unsupported archive types** at the few entry points that lacked
  it: `decompress_to_writer` and `qz merge` now reject cluster (type 3) and single-end
  reference (type 4) with actionable messages, and the PyO3 binding rejects *paired*
  cluster archives (new `is_cluster_paired_archive` predicate) while still accepting
  single-output types (single-end cluster, type 4) like single-end default. Plus
  doc-drift fixes (`archive_type` enumeration, stale `verify_reference_single` comment)
  and a stricter 4-line-record parser in the cluster roundtrip test helpers.

### Feature: single-end reference-based compression (`--reference` with one input)

`--reference` now accepts a single `-i` input (single-end), not just two
(paired-end). A single-end reference archive is a new v5 `archive_type` 4: reads are
mapped against the reference FASTA and encoded as edits over the covered slices, with
off-reference reads spilling to a per-chunk fallback literal pool, qualities via
fqzcomp, and headers via independent columnar (no delta — there is no R2). The
reference is never stored and is never needed to decode.

This is built as a PARALLEL single-end path (`reference/stream_single.rs` +
single-end helpers in `encode.rs`) that reuses the mate-agnostic per-read codecs, the
shared worker pool, and the shared globals/footer writers — the paired `archive_type`
2 path is BYTE-FOR-BYTE UNCHANGED (verified: the `streaming_equiv` FNV goldens and the
`paired_*` byte-identity tests all pass, and single-end-default type-0 output is
unchanged). Encode emits a real structural self-check
(`validate_reference_directory_single`, a mate-1-only contract) before the footer, so
a malformed single-end directory is rejected at encode time.

The matching single-end **DECODE** now lands too: the v5 decoder accepts
`archive_type` 4, and `decompress` routes it to `decompress_reference_single` — one
mate reconstructed reference-free (absolute within-read edit positions, no anchors)
through the same bounded-RAM, seekable block-streaming engine as paired reference, and
written to ONE FASTQ atomically (no `_R1`/`_R2` suffix, like single-end default).
Single-end reference gets its own `is_reference_single_archive` classifier (type 4) so
single-output API consumers decode it like single-end default; `is_reference_archive`
stays type-2-only. The paired `archive_type` 2 path remains BYTE-FOR-BYTE UNCHANGED
(`reference_roundtrip_byte_identical` + the `streaming_equiv` FNV goldens still pass).

`qz verify` (deep + fast) and `--numa auto|N` sharded decode now cover type 4 too, at
full parity with the other modes. Verify mirrors the paired path
(`verify_reference_single`: deep reconstructs one FASTQ to a single CRC32 with
`r2_crc32 = None` and `num_reads` NOT doubled; fast walks per-block CRCs without
reconstructing). Under `--numa`, type 4 shards exactly like single-end/paired/reference
— each node-pinned worker loads the whole-archive consensus globals once, decodes a
disjoint chunk-range to one part file (reference archives carry no `ChunkDecodedSizes`,
so the driver uses part-file assembly, not direct-write), and the driver concatenates
parts in read order into ONE output.

Measured on HG002 chr20 R1 alone (single-end, 10.81M reads, 3.74 GB) vs `chr20.fa` on
an idle 72-core box, byte-identical lossless roundtrip (`cmp`) + clean deep/fast verify:
single-end reference **22.16×** vs single-end default **8.48×** — the reference archive
is **~2.6× smaller** (168.6 MB vs 440.3 MB), at ~4.3 GB compress / ~1.6 GB decompress
peak RSS. The lift mirrors the paired reference win (8.26× → 22.73×): reference encoding
pays off on any covered data, not just paired. Reference mode's constraints are unchanged
for single-end: **Auto quality only**, **no stdout**, **no FASTA**, **no
`--reference-index`**; edit positions stay **absolute within-read** (the ORA patent
design-around — never delta).

### Performance: reference compress — skip the wasted R2-header columnar; word-aligned coverage mark

Two byte-identical reference-mode (`--reference`) compress optimizations surfaced by
an audit of the encode path:

- **R2 header predict-and-skip.** The reference encoder computed BOTH an independent
  columnar candidate AND a delta-vs-R1 candidate for every mate-2 header chunk and
  kept the smaller — but on paired Illumina data R2 ids are R1 ids with the mate
  suffix flipped (or identical), so the delta is a few hundred bytes and ~always wins
  while the discarded R2 columnar is a full BWT over the chunk's headers. Reference now
  uses the same cheap-first skip the FASTQ-paired path already had: compute the delta
  first, and when it is ≥2× below the R1-columnar proxy (mate-1's serialized header
  length), skip the R2 columnar BWT entirely. Falls back to compute-both when the delta
  is not a clear win. Verified byte-identical on a 3-chunk HG002 chr20 slice (R2
  headers reduced to a 110-byte delta/chunk vs the ~276 KB columnar that was being
  discarded) and across the reference test suite.
- **Word-aligned coverage mark.** `CoverageMap::mark` set the per-read coverage bitmap
  one bit at a time (~150 ops per 150 bp read) in the serial producer fold; it now
  fills whole `u64` words in the interior and masks only the two boundary words
  (~O(len/64) ops), shrinking the serial dispatch gate. Proven to set exactly the same
  bits as the old loop by a 20 000-span equivalence test over all alignment/boundary
  cases, so `finalize` and the archive are byte-identical.

Both verified byte-identical to the prior encoder: master vs branch produce the same
archive (identical sha256) on the 3-chunk chr20 slice.

### Docs: reference-mode comment/doc drift fixes (no behavior change)

A follow-up audit of the reference-based **compress** path corrected stale
comments that no longer matched the code (all comment-only — output unchanged):
the `reference_window` doc-comment said "2 is the balanced default" but the CLI
default is 4; `coverage_to_intervalmap`'s "(cheap)" comment now notes it is an
O(genome-length) serial pass; the `REF_CHUNK_RECORDS=500K` rationale is flagged
as measured on the pre-worker-pool engine (a ~0.5% ratio knob, not a validated
speed setting); the stale `reference_integration.rs` R2-header test comment now
describes the live columnar-vs-delta decision; and `docs/reference-paired-hotpath.md`
gained a corrections block (notably: `QZ_REF_TIMING` was removed in the worker-pool
refactor and no longer exists, and `mu`/`sigma` feed the mapper's pair-rescue
rather than being unused).

### Performance: `--numa` sharding no longer gated off by reclaimable page cache

The `--numa` decode/compress gate filtered out any node whose strict `MemFree` was below a
worker's peak RSS. On I/O-heavy boxes a node can hold tens of GB of **reclaimable** page
cache (file cache + reclaimable slab), leaving little `MemFree` even though that memory is
freely available to a worker — so `--numa auto` silently fell back to in-process and the
cross-socket sharding speedup was lost.

The gate now estimates per-node **available** memory (`MemFree + Active(file) +
Inactive(file) + SReclaimable` — a per-node analogue of `/proc/meminfo` `MemAvailable`,
which node meminfo lacks) and filters on that. Anonymous/unreclaimable memory is excluded,
so the estimate never over-commits; a node under genuine pressure (low free AND low
reclaimable) is still rejected, and sharding still falls back in-process when fewer than two
nodes qualify.

Measured on a 2-socket box whose node 0 held ~37 GB reclaimable cache (1.6 GB `MemFree`):
reference decode of a 10.8M-pair chr20 archive previously fell back to in-process
(~41.8 s); it now shards 2-way and decodes in **~27.1 s (`--numa auto`) — 1.54× faster**,
byte-identical output. This narrows the decode gap vs SPRING (~22 s on the same data) from
~1.9× to ~1.23×. Applies to single-end, paired, and reference decode and to NUMA compress
sharding.

### Refactor: single-allocation reference-decode read reconstruction

Reference decode reconstructed each mapped read in **three** allocations
(`PackedBacking::slice()` → `reconstruct_read`'s `to_vec()` copy →
`reverse_complement()`'s output `Vec`). It now uses a fused
`reconstruct_read_from_backing` that builds the read in **one** allocation:
`PackedBacking::fill_slice` writes the consensus slice into a single buffer, edits are
overlaid in place, and the reverse-complement is applied in place via
`reverse_complement_in_place`. `slice()` is refactored to delegate to `fill_slice`
(single source of truth).

**Byte-identical** decode output — verified by a unit sweep over every
(base, len, is_rc, edit-set) combination plus the full reference suite (485 qz-lib
tests green) and a sha256-identical decode of a real 10.8M-pair chr20 archive.
Performance is **neutral** (measured ~0 wall change, 40.2 vs 40.3 s in-process, within
noise): reference decode is memory-bandwidth-bound and the reconstruction's
allocations/copies are a small fraction of total decode traffic (already cheap under
mimalloc). Kept as a cleaner, lower-allocation path with two reusable helpers
(`fill_slice`, `reverse_complement_in_place`).

### Performance: paired `--cluster` flush↔sort pipeline + per-field arenas

The paired cluster compressor previously compressed each chunk **synchronously** inside the
sort callback (the producer blocked during zstd, and the last chunk flushed after all
sorting finished — a serial tail). It now runs a **pipeline**: the sort streams completed
chunks to a background flush thread over a depth-1 channel (backpressure caps in-flight at
two chunks), which compresses and writes blocks in chunk order. The per-chunk buffer is also
restructured into per-stream flat arenas (`ChunkArenas`) — eliminating the per-record clone
and ~65M small allocations, and letting each mate's sequence/quality/header arena free as
soon as its block is produced, so the in-flight chunk drains while the next fills.

**Byte-identical** output at a given chunk size (verified `cmp`-equal on HG002 chr20 paired
AND single-end, roundtrip CRC unchanged, across chunk sizes from one chunk to 1400+).
Measured chr20 paired: peak RSS ~13.5 → **12.0 GB** (per-field free offsets the second
in-flight chunk) and compress ~**6% faster**. The speed gain is modest because cluster
compress is memory-bandwidth-bound — overlapping the (bandwidth-bound) flush with the
(bandwidth-bound) sort cannot exceed per-socket bandwidth — so the win is weighted toward
RAM. Single-end keeps its synchronous flush (unchanged). No ratio change.

### `--cluster` zstd sequence compress now scales workers with `-t` (per-mate cap 48, was 16) + `QZ_CLUSTER_ZSTD_WORKERS`

The clustered sequence stream is compressed with zstd long-distance matching, whose per-mate
worker count (`ZSTD_c_nbWorkers`) was hard-capped at 16 — leaving the seq stage ~1.9× slower
than a many-core box can sustain. `resolve_zstd_workers` now derives the per-mate count from
`-t`, split across the concurrent flush lanes (paired flushes both mates via `rayon::join`, so
it passes 2 lanes — keeping the *total* worker count ≈ `-t` and never oversubscribing small
boxes), then clamps to `[1, 48]` (the measured MT-throughput plateau; each worker preallocates
a job buffer, so the cap bounds per-chunk RSS at extreme `-t`). `QZ_CLUSTER_ZSTD_WORKERS=N`
overrides with an explicit per-mate count (escape hatch; bypasses the `-t` derivation and lane
split).

For `-t ≤ 16` the resolved count is unchanged, so small-box behaviour is preserved bit-for-bit.
zstd-MT job boundaries depend only on `(input, level, window_log)`, not the worker count, so the
**archive is byte-identical** across worker counts (verified by md5); cluster's order-dropping
roundtrip oracle (multiset checksum / `qz verify`) is unaffected. Measured at `-t 72`,
byte-identical and verify-PASS both: single-end 5M **49.9 s → 37.6 s (1.33×)**, RSS +0.9 GB;
paired chr20 10.8M **121 s → 102 s (1.19×)**, RSS flat (+0.02 GB).

### Internal: cluster compress/decode memory hygiene (byte-identical)

Byte-identical cleanups that remove redundant copies in the cluster path: the sequence
stream is fed record-by-record into the zstd encoder instead of being re-concatenated into
a flat `Vec` per chunk (`compress_seq_zstd_long_records`); splitter sampling hashes in
bounded batches retaining only `u64` keys (was holding ~1.2 GB of sampled sequences); and
decode passes unflipped reads through as borrowed slices instead of `to_vec` (only flipped
reads materialize a reverse complement). Output is byte-for-byte unchanged (verified on
HG002 chr20 paired). Peak RSS at chr20 scale is unchanged — the cluster peak is bound by
whole-chunk residency, not these transients — but the seq-feed change keeps per-chunk RAM
low for the upcoming flush↔sort pipeline.

### Added: `--cluster` order-drop compression mode (archive_type 3) — single-end and paired-end

`qz compress --cluster [fast|balanced|max]` reorders reads into minimizer clusters and
compresses the concatenated, reverse-complement-canonicalized sequence stream with zstd
long-distance matching instead of BSC. The archive is **set-lossless**: every read
(header, sequence, quality) is preserved exactly, but reads are returned in a
deterministic clustered order, not input order — so use it only when downstream tools do
not depend on read order. `qz decompress` reconstructs the FASTQ; `qz verify` checks
integrity via an order-independent multiset checksum (and a fast per-block CRC mode).
Levels map to zstd 12/16/19 (default `balanced` = 16). Conflicts with `--ultra` and
`--reference`. Compression stages reads through external on-disk buckets in the
`-w`/`--working-dir` directory (scratch ≈ input size, auto-removed on completion or
abort), keeping peak RAM bounded and constant in read count. The zstd-long sequence
codec is multi-threaded (per-mate worker count from `-t`, capped at 48; overridable with
`QZ_CLUSTER_ZSTD_WORKERS`).

**Paired-end** (`qz compress --cluster -i R1 -i R2`): reads are clustered as
R1-keyed pairs — R1 drives the minimizer sort, R2 follows its mate. `qz decompress`
emits two FASTQs (`-o R1.out.fastq -o R2.out.fastq`), each set-lossless relative to
its input, with pairing preserved at read-level: read i of R1-out and read i of
R2-out are the same pair. `qz verify` checks integrity via a pair-tuple multiset
checksum (each (R1-record, R2-record) pair is hashed as a unit). The paired variant
uses `zstd --long` on R1 and R2 sequence independently; header/quality streams are
compressed by the standard BSC path.

The benefit is **coverage-driven** — it pays off where reads overlap. Measured on
deep-coverage data (HG002 chr20, ~25×, 10.8M single-end reads): **16.9×** vs default
**8.5×** (balanced/z16, archive ~halved), compressing in ~76 s. On shallow data (10M
WGS reads, ~0.5× coverage) the gain is marginal (~+7%), so reserve `--cluster` for
high-coverage inputs.

### Improved: RC-invariant cluster key — denser `--cluster` archives (single-end −6.4%, paired −4.7%)

The cluster sort key was not fully reverse-complement-invariant: open-syncmer selection
(`ts=[0]`, the smallest s-mer at offset 0 of the k-mer window) is strand-asymmetric, so a
read and its reverse complement could select different syncmer positions and land under
different keys — failing to co-cluster even though they cover the same locus. The key now
scans **both full orientations** and takes the smaller open-syncmer hash
(`strand_symmetric_open_syncmer_hash_pos` / `canonicalize_open_syncmer_strand` in
`dna_utils.rs`), so a read and its RC are guaranteed the same key. Used for splitter
sampling, single-end records, and both mates of a pair (R1 still drives pair order; each
mate is canonicalized independently). Because `reverse_complement_canonical` is a true
involution, the new key is provably RC-invariant for arbitrary read bytes (regression
test `strand_symmetric_open_syncmer_is_rc_invariant`).

Output remains **set-lossless** (verified: byte-identical multiset checksum before and
after — single-end `b13105e4`, paired `e2a9c7c7`) and FTO-clean (the key is a function of
one read's own two orientations — no cross-read comparison). Measured on HG002 chr20
(deep coverage), balanced/z16, byte-verified lossless: **single-end 16.92× → 18.07×
(−6.36% archive)**, also ~8% faster (better clustering gives zstd longer matches → less
codec work); **paired 14.02× → 14.72× (−4.73% archive)**, ~6% slower compress (the
both-mates key scan roughly doubles, and the smaller per-mate seq saving pays back less of
that than single-end's larger one). Peak RSS unchanged. Not a format change — old cluster
archives still decode; re-compress to get the denser layout. The historical forward-only
helpers (`compute_min_syncmer_hash`, `canonical_minimizer_orientation`) are retained for
non-cluster callers.

### Performance: paired `--cluster` compress ~2× faster (parallel bucketing pass)

The paired clustering pass (`cluster_sort_paired`) previously partitioned reads into
on-disk buckets in a single **serial** pass — canonicalizing both mates, computing R1's
minimizer key, and folding the verify checksum one pair at a time. It now batches pairs
into a flat arena and computes the per-pair canonicalization, key, and checksum **in
parallel** (mirroring the single-end clustering pass), then folds the checksum and writes
buckets sequentially. Measured on HG002 chr20 (10.8M pairs, 72 threads): the bucketing
pass dropped 191 s → 24 s (~8×) and total compress 339 s → 170 s (~2×). Archives are
**byte-identical** to the old serial pass (the clustered order is fixed by a deterministic
within-bucket sort, and the checksum fold is commutative), and peak RAM is unchanged
(~12.5 GB — the batch working set is freed before the codec-flush peak).

### Performance: paired `--cluster` compress ~1.3× faster (concurrent-mate codec flush)

The paired per-chunk codec flush was internally serial — it compressed mate R1's four
streams (clustered sequence via zstd-long, qualities via fqzcomp, headers, strand bits)
and *then* mate R2's, while the dominant zstd-long sequence stream was then capped at 16
worker threads. On a many-core box this left most cores idle. The two mates' codec work now runs
**concurrently** (`rayon::join`); the compressed blocks are then written in the same fixed
R1-then-R2 order, so the output is **byte-identical** (verified `cmp` + matching sha256).
Measured on HG002 chr20 (10.8M pairs, 72 threads): the flush dropped 126 s → 87 s (~1.45×)
and total compress 170 s → 130 s (~1.3×). The tradeoff is peak RAM: holding both mates'
compressed blocks plus two live zstd-long contexts at once raises peak RSS ~12.6 → ~16.4 GB
(+3.8 GB) — still bounded by the chunk size and constant in read count, not the read total.

### Fixed: `qz decompress --numa N` no longer fails on cluster archives

Decoding a cluster archive (`archive_type` 3) with an explicit `--numa N` (N ≥ 2) on
NUMA hardware previously failed with an internal `unknown archive_type 3` error: cluster
archives use a bespoke, order-dropping decode path that is not chunk-range shardable, and
the NUMA planner did not recognize the type. `--numa N` now detects cluster archives and
decodes them in-process (the correct path), printing a one-line note that NUMA sharding
does not apply. `--numa auto` already fell back gracefully and is unchanged. Affects both
single-end and paired cluster archives.

### Changed: faster columnar header compression (~1.16× faster compress on Illumina, +0.15% archive)

The columnar header path was the single biggest remaining serial bottleneck in
single-end/paired compress. Two fixes, found by phase-timing the streaming writer:

1. **Parallelized the generic-token columnar encoder.** `leaf_tokenize_block` (per-read
   tokenization) and the generic-token Phase-1 column classification were serial
   `O(reads)` / `O(columns × reads)` loops — ~3.8 s of the ~5 s per-chunk header job on
   10.8M-read NovaSeq. Both are per-read/per-column independent, so they now run under
   `par_iter` (byte-identical output — same algorithm, parallel execution). Header job
   ~2.1× faster on its own. (This path is hit by *comment-less* Illumina headers
   — `@inst:run:fc:lane:tile:x:y` with no trailing ` 1:N:0:…` — which `detect_format`
   routes to the generic-token codec.)

2. **Made the header path block-granular.** Sequence and quality were already split into
   many block-sized parallel jobs, but the columnar header stayed one coarse per-chunk
   unit (one ~5 s job per 2.5M-record chunk) that straggled at the very end of compress
   while most workers sat idle — phase-timed at **~47 % of the compress wall** spent with
   only the writer + a couple of header stragglers running. `columnar_blocks_capped` now
   splits each chunk's headers into independent **250K-record sub-blocks** compressed in
   parallel (the decoders already concatenate N blocks in order, so the on-disk format is
   unchanged). The header job drops **~5.1 s → ~1.2 s per chunk (4.4×)**, collapsing the
   tail. `QZ_HEADER_SUB` overrides the sub-block size for tuning (0 = legacy whole-chunk).

Net measured (10.8M-read NovaSeq chr20 R1, single-end, quiet 2-socket Xeon): **compress
~14.1 s → ~12.1 s (~1.16×)** and far more wall-stable (the tail-straggler variance is
gone). Cost: **+0.15 % archive size** (each 250K sub-range loses a little cross-boundary
columnar context + rebuilds its own combo dict) — a rounding error on a stream that is
~4.5 % of the archive. Lossless (roundtrip + deep `verify`); archives ≤250K records are
byte-identical to before (no split). Decode unchanged.

### Added: NUMA-aware paired-end compress (`qz compress -i R1 -i R2 --numa`)

`qz compress --numa` now shards **paired-end** FASTQ across NUMA sockets, the same way
single-end and reference-based compress already do. The driver splits both mates on a
shared record-index boundary, spawns `N` node-pinned worker processes (each compresses a
record-aligned `(R1, R2)` byte-range pair to a part archive), then merges the parts into
one paired v5 archive. Previously paired input fell back to a single in-process compress.

The paired half reuses the streaming compress core: `compress_paired_streaming` was
refactored into a generic `compress_paired_streaming_from_readers<R: BufRead>` (the
whole-file path and the worker's bounded `BufReader<Take<File>>` range path share it —
output is byte-identical to the in-process encoder). The serial driver wrapper was kept
small to preserve the sharding win: the input split is a single-pass `memchr` newline scan
over both mates **concurrently** (`std::thread::scope`); the plan prelude counts only a
capped 200K-record sample (paired forces variable-length framing, so there is no
constant-length edge to detect); and `merge_paired_core` mirrors the single-end merge
(footer-order verbatim frame copy + chunk renumber, `validate_paired_directory`, no
`ChunkDecodedSizes` rebuild — paired has none). A mate-count mismatch between R1 and R2 is
detected during the split and fails cleanly with no partial output.

Measured on HG002 chr20 NovaSeq 2×151 (10.8M pairs, RAM-backed, 2-socket Xeon):
`--numa 2` vs `--numa off` = **~1.12–1.14× wall** (steady-state), byte-identical-lossless
roundtrip on both mates + deep `verify`. The serial wrapper is ~1.9 s (split 1.1 + prelude
0.06 + merge 0.83) on top of an ~18.6 s parallel core that captures ~1.27× — recovering the
win the earlier parked branch lost (0.95×) to a heavy serial wrapper. No on-disk format
change.

**RAM:** NUMA-compress uses *more* total RAM, not less — the two workers run concurrently,
so the true subtree peak is **16.5 GB vs 12.3 GB in-process (+35%)**. (A `/usr/bin/time -v`
"max RSS" reads ~8.5 GB because it reports only the largest *single* worker, not the
concurrent sum — don't quote it as a saving.) The one genuine RAM property is **per-node
locality**: each worker's footprint is pinned node-local (`MPOL_BIND`), so it stays on its
socket (~8.5 GB/node) — relevant only if you are per-node memory-constrained. Like
single-end, this is an opt-in wall lever (`--numa off` is the default; `auto` never loses
to in-process), not a memory saving.

### Changed: NUMA-compress plan prelude no longer materializes 2.5M records (~3× faster prelude, byte-identical)

`qz compress --numa` resolves a shared compression plan once in the driver before
forking the node-pinned workers. That resolution (`resolve_plan_override`) read and
fully materialized the entire first chunk — `chunk_records` (default **2.5M**) records,
~915 MB, three `String` allocations each — just to compute three facts: the record
count (vs the 100K fqzcomp threshold) and the constant seq/qual lengths. Because it runs
**serially before the parallel worker core**, it was the dominant Amdahl tax on the
sharded speedup: phase-timed at **~1.16 s of a ~10 s** `--numa 2` single-end compress,
which (with the ~0.37 s merge tail) capped production sharding at ~1.03× even though the
parallel worker core matches the manual-`numactl` recipe (~1.2×).

Plan resolution now walks the **same `chunk_records` window** with a cheap single-pass
line-length scan (`scan_plan_stats`) that never builds `FastqRecord`s — byte-for-byte the
same plan (proven by a parity test vs `read_chunk_records` + `detect_const_lengths` across
fixed/variable/CRLF inputs and window sizes, quality on/off), so the archive is unchanged
and lossless. Prelude **1.16 s → 0.39 s**; single-end `--numa` compress ~1.03× → ~1.13×.
FASTA input (multiline records) keeps the full-reader path (defensive; the NUMA precheck
already excludes FASTA). A new gated `QZ_NUMA_TIMING` env var prints the driver's
split/prelude/workers/merge phase breakdown. No on-disk format change.

### Changed: reference mode R2 headers use delta-vs-R1 (~5% smaller reference archives, lossless)

Reference-based paired compress (`--reference`) now makes the same columnar-vs-delta R2
header decision the FASTQ paired path makes, instead of always encoding R2 headers as
independent columnar. On typical Illumina data R2 ids are R1 ids with the mate suffix
flipped (`…/1` → `…/2`), so the delta-vs-R1 representation is a few hundred bytes per
chunk vs megabytes of columnar — previously reference threw that win away. The encoder
(`reference/encode.rs`) builds both candidates per chunk and keeps the smaller (emitting
`R2HeaderDelta`/BSC when it wins, else `R2HeadersIndep`/columnar); the decoder
(`reference/decode.rs`) reconstructs both mates' ids per chunk against R1 when any chunk
uses delta (bounded-RAM — one chunk's ids resident), keeping the streaming columnar path
for columnar-only archives. The `R2HeaderDelta` role and §9 validator already reserved
this; only the encoder/decoder needed wiring (closing an incidental divergence from
FASTQ-paired header handling — same `header_delta` codec is reused).

Measured HG002 chr20 (10.8M pairs): R2 headers **17.5 MB → 3.4 KB**, reference archive
**328.7 MB → 311.3 MB (−5.32%)**, byte-identical roundtrip both mates + deep `verify` OK;
compress wall unchanged (reference is mapping-bound; header columnar isn't the long pole).
Gated by the 5 `m_ref_*` byte-identity goldens (regenerated) + a `reference_integration`
case asserting both the delta (similar headers) and independent-columnar (unrelated
headers) paths roundtrip. Pre-existing reference archives (independent-columnar R2) still
decode unchanged.

### Changed: paired compress skips the redundant R2-header columnar pass (~25% faster, byte-identical)

The paired `PairedHeaders` job used to compress R2 headers **two ways every chunk** —
independent columnar (a full BWT/BSC over millions of headers) *and* delta-vs-R1 — then
keep the smaller. On typical Illumina data R2 ids are R1 ids with the mate suffix flipped
(`…/1` → `…/2`), so the delta is a few hundred bytes and **always wins**, while the
columnar candidate is ~10⁴× larger output and ~80× slower — computed and thrown away
every chunk. The job now computes the cheap delta first and **skips the R2 columnar pass
when the delta clearly wins** (using the already-computed R1 columnar size as a free,
accurate proxy for what R2 columnar would cost); it falls back to the exact
compute-both-and-compare otherwise, so output is **byte-identical** to the previous
encoder. Measured HG002 chr20 (10.8M pairs, /dev/shm, `-t 72`): paired compress **30.8s →
~22.8s (~1.35×)**, archive byte-for-byte identical (904,106,608 bytes), deep `verify` OK.
Gated by the 5 `m_paired_*` streaming-equiv FNV goldens (incl. the independent-columnar-R2
case) + `paired_integration`.

### Added: NUMA-aware multi-process single-end compress (`--numa`)

- `qz compress --numa auto|off|N`: NUMA-aware multi-process single-end compress
  (shards the input by record range across NUMA nodes, compresses each shard on a
  node-pinned process, merges the part archives). `auto` is strictly never worse
  than in-process. Paired-end support follows.

### Changed: reference compress driven by the shared worker pool (Increment 3, byte-identical)

Reference-based compress (`--reference`, `archive_type 2`) now rides the shared
block-granular concurrency machinery instead of a bespoke pipeline. The per-chunk
`encode_chunk` work is dispatched through the generic `spawn_worker_pool<J,R>` (the
same pool single/paired use, parameterized with reference's own `RefJob`/`RefResult`),
and a thin chunk-ordered writer (`reference/stream.rs`) reorders the encoded chunks
and streams the whole-archive globals in its tail. The bespoke producer-thread +
bounded channel + token-semaphore consumer + inline reorder/flush loop is **deleted**.

Reference's alignment front end stays outside the shared block-streaming core (the
justified divergence, mirroring the decode side): strobealign mapping, diff/placement,
the coordinate/edit columns, the serial in-order coverage fold, and the consensus
globals (PackedBacking/NBitmap/IntervalMap/ReferenceMeta) are unchanged. `encode_chunk`
and the global block-streaming writers are byte-for-byte unchanged; the producer pulls
pairs prefix-then-readers chunked at `chunk_records` and folds coverage serially in
chunk order exactly as before, so archives are **byte-identical** to master (gated by 5
reference byte-identity goldens + the 42-test `reference_integration` suite + the
`qz-cli` NUMA/CLI tests). Worker count follows reference's `--reference-window` (heavy
chunk jobs must not run block-wide); `QZ_COMPRESS_WORKERS` overrides. Single-end,
paired-end, and reference now all share the container, codecs, stagers, worker pool,
and footer — only reference's front end and a thin chunk writer sit outside.

Real-data A/B (HG002 chr20 NovaSeq 2×151, 2M pairs / 4 chunks, vs the master `c53aada`
binary): **byte-identical archives** at `-t 16` (both runs), `verify` deep OK, and a
lossless roundtrip (decompressed R1/R2 byte-identical to the inputs). Peak compress RSS
is ~10% **lower** than the legacy encoder (~2.87 GB vs ~3.3 GB) and flat across `-t 8/16`
(bounded by `--reference-window`, constant in read count); wall is ~15–20% slower — the
expected unification/coordination tax (this increment is one-engine consolidation, not a
speed win; reference's pipeline/bandwidth-bound encode is the lever for a later pass).

### Fixed: harden v5 directory validators + the streaming ordered-buffer against malformed/internal-bug inputs

Three robustness fixes (none change valid-archive output — all bytes/roundtrips identical):

- **Paired directory validator now matches the decoder contract.** `validate_paired_directory`
  (also what `verify --fast` trusts) accepted codecs the decoder doesn't support: a
  `Headers`-role entry tagged `BSC` (Headers is always columnar-decoded), a `Qual` entry
  tagged the removed `QUALITY_CTX` (rejected at decode), and a mate whose `Qual` codec
  varied across chunks (streaming decode picks one mate-wide codec from the first entry).
  Such a forged/legacy archive could pass validation yet misdecode. The validator now allows
  only `Headers`=columnar, `Qual`∈{BSC,FQZ}, and enforces per-mate `Qual`-codec uniformity.
- **`OrderedBlockBuffer` fails cleanly on internal bugs in release.** A duplicate
  `(chunk, mate, role, block_idx)` block or a duplicate header segment was a `debug_assert`
  only (a silent overwrite in release); completion checked block *count*, not that indices
  were exactly `0..n`. A producer/worker bug could therefore ship a corrupt archive instead
  of aborting. These are now hard errors surfaced as a first-error-wins writer abort, and
  completion requires contiguous block indices.
- **Malformed single-end footer can no longer panic.** `validate_single_end_directory`
  indexed `counts[0]` before checking that the mandatory Headers+Sequence streams exist, so a
  forged footer using only non-core roles (e.g. Nmask-only) hit an out-of-bounds panic (a DoS
  on hostile input). The mandatory-stream check now runs first, turning that into a clean error.

### Fixed: `cargo test --release` green (debug-only fault-injection test gated)

`decode_chunk_range_bounded_overrun_is_integrity_error` depended on the
`QZ_NUMA_FAKE_UNDERSIZE` writer hook, which is `#[cfg(debug_assertions)]`-gated (the
fault injector must not ship in release binaries). Under `cargo test --release` the hook
was inert, so the test's `unwrap_err()` panicked and the `numa_chunk_range` suite was red
in release. The test is now `#[cfg(debug_assertions)]`-gated to match its hook; the
production integrity check it exercises (`DirectWriteIntegrityError`) is compiled in all
profiles and stays covered in debug. Also corrected a stale "v4 frame" doc comment in
`write_role_blocks` to "v5 block frame".

### Perf: single-end compress rewritten as a block-granular streaming engine

Single-end compress (default and `--ultra`) is now a block-granular streaming pipeline —
producer frames records and stages per-role raw blocks → bounded worker pool compresses
each → single ordered writer emits the archive — replacing the chunk-granular engine,
which is deleted (one engine, no `CompressEngine` selector). Archives are **byte-identical**
to the old engine (pinned by FNV-1a goldens and verified by a release-binary A/B on real
data); decode, format, and all flags are unchanged.

Peak compress RSS is now the worker-pool working set, **constant in read count** instead of
the read-chunk, and the worker pool **follows the `-t` thread setting** (the user's RAM/speed
dial; `QZ_COMPRESS_WORKERS` overrides; ultra runs one block at a time). Measured on 10M-read
WGS, byte-identical at every setting: `-t 8` → 3.1 GB, `-t 16` → 4.9 GB, `-t 32` → 6.4 GB peak
RSS, versus the old engine's ~11 GB; at `-t 32` streaming is also faster than the old engine.

Compress is now atomic: output streams to a sibling `.tmp` and renames on success, so an
aborted compress (record error, worker error, or panic) leaves no partial archive at the
final path (stdout output bypasses the temp, as before).

### Feat: paired-end compress converged onto the block-granular streaming engine (Increment 2)

Paired-end compress now runs on the **same** block-granular streaming engine as single-end,
generalized to be mate-aware: one `OrderedBlockBuffer` (keyed by `(mate, role)`), a
segment-based `ChunkManifest`, and a `write_ordered` that branches on a `CompressMode` for the
per-mode physical layout — single-end emits one directory entry per block; paired emits one
entry per `(mate, role)` segment (a `write_block_stream` blob, fixed 6-segment order, the R2
header chosen as the smaller of independent-columnar vs delta-vs-R1, no `ChunkDecodedSizes`).
A second producer reads the two mate files **record-by-record** (no chunk-Vec buffering), so
only the first chunk is materialized (the planning prelude that resolves per-mate quality).

Archives are **byte-identical** to the previous paired encoder — proven by 5 FNV goldens
(both R2-header branches, BSC + fqzcomp quality, variable + const length, single + multi
chunk), a real-data release-binary A/B at multi-chunk scale (HG002 chr20, 10.8M pairs / 5
chunks: byte-identical at `-t 1/16/32`, roundtrip byte-identical, verify deep + fast), and the
existing paired roundtrip suite. Decode and verify are unchanged.

Peak compress RSS now **follows the `-t` dial** and is **constant in read count** instead of
the chunk-granular legacy's fixed per-chunk hold. Measured on the full 10.8M-pair dataset
(byte-identical at every setting): **6.6 GB (`-t 8`) / 11.2 GB (`-t 16`) / 11.8 GB (`-t 32`)**,
versus the legacy encoder's **14.2 GB (`-t 16`)**. `-t 8` is the memory-lean setting (≈ the
design's ~5–7 GB paired floor); higher `-t` trades RAM for parallelism. The paired floor sits
above single-end's because paired carries two mates' working set plus the per-segment
compressed-chunk hold the paired directory layout requires.

The legacy `compress_paired_v5` path is removed, along with the now-orphaned chunk-granular
`compress_one_chunk` chain (`ChunkParams`/`ChunkCompressResult`/`compress_stream_to_bsc_blocks_v4`/
`records_to_streams`/`paired_params` + obsolete tests) — single-end and paired now share one
compress engine (only reference still uses a per-chunk sequence path). Single-end output is
unchanged.

### Perf: reference encode borrows headers/qualities instead of per-chunk deep-clone

Reference (`--reference`) encode no longer deep-clones every read's id and quality
into owned `Vec<Vec<u8>>` per chunk. `encode_chunk` now builds `Vec<&[u8]>` borrows
into the chunk it already holds (the columnar header codec reads only id bytes and
fqzcomp reads only quality bytes, both as `&[u8]`), and the dead all-`None`
`r1_anchor` throwaway `Vec` for mate 1 is replaced by `Option<&[…]>` (`None` for
mate 1, whose anchors are unused). This removes ~176 MB/chunk of dead DRAM copy that
was allocated and freed immediately.

**Measured** (2M novaseq 2×151 pairs vs chr20, 4 chunks @ 500K, same-load interleaved
A/B, n=6, T=72): **−11% compress wall (14.45→12.86 s median) and −20% peak RSS
(5.98→4.79 GB)**, byte-identical archive (md5 + size unchanged). The win is larger and
more consistent than a typical clone-elimination because it removes *bulk* copies in
the encode prologue (not per-record allocs in the read phase, which previously
*regressed* on this bandwidth-bound box) — fewer bytes on the bus on a core-starved,
memory-bound encoder. All lib/single/paired/reference tests pass.

Follow-on: the remaining per-chunk clones are now also eliminated. `encode_mate_buf`
takes the full pair-slice and selects this mate's `Placed` per-read (`.0`/`.1` by
`mate`) instead of cloning every `Placed` (and its inner `Vec<Edit>`) into flat
`p1s`/`p2s` Vecs, and `edit_recs` borrows the edit slices (`encode_edits` now takes
`&[&[Edit]]`) instead of cloning them. Incremental (same A/B, n=6): **−1.6% wall and
−4% peak RSS** on top of the headers/quals borrow (combined vs pre-change master:
≈ −12.5% wall / −23% peak RSS for reference compress). Byte-identical; adversarially
reviewed byte-identical + lifetime-sound; all tests pass.

Final per-read clone removed: the fallback-literal collector borrows each unmapped
read's sequence (`fallbacks: Vec<&[u8]>` from `chunk`) instead of `seq.to_vec()`, and
`write_fallback_pool` takes `&[&[u8]]` (it only reads the literals). Byte-identical.
The payoff is **data-dependent and small** — bounded because the literal's first copy
already lives in the `chunk` records: on well-mapping data (few fallbacks) it's
**neutral** (RSS flat); at a pathological ~100%-fallback tiny-reference it's **−4% peak
RSS**. So it helps only on low-mapping/divergent inputs. No `unsafe`; all tests pass
(6 `fallback.rs` test callsites updated via an `as_slices` helper).
### Added
- `qz decompress --numa auto|off|N`: NUMA-aware multi-process decode. On real
  multi-socket hardware `auto` shards the archive's chunk-ranges across one
  node-pinned worker process per NUMA node (node-local memory), assembling the
  part files into the final output in read order — ~1.6–1.7× faster decode on a
  2-socket Xeon. Off NUMA hardware (single socket, laptop, VM, non-Linux) `auto`
  falls back to the existing in-process decode with zero overhead. Single-end and
  paired, AND reference archives are all sharded. Plain FASTQ output is byte-identical
  to single-process decode; single-end gzip output is a valid multi-member stream
  with identical decompressed bytes.
- v5 single-end archives now carry a whole-archive `ChunkDecodedSizes` global
  directory block (one raw, uncompressed per-chunk decoded-output byte table),
  letting the NUMA decode path pre-partition work without decoding blocks. The
  table is `no_quality`/FASTA-aware and FIFO-aligned with chunk order; normal
  single-end decode ignores it (the role filter skips global-sentinel entries),
  so roundtrip output is byte-identical.
- Verified direct-write path for single-end NUMA decode: `decode_chunk_range` now
  takes an optional `DirectWriteRegion { offset, len }` so a worker writes its
  chunk-range straight into a pre-sized shared output file at an offset instead of
  a part file. The worker re-derives `offset`/`len` from the `ChunkDecodedSizes`
  table and rejects any mismatch, hard-caps the decode write at `len` (catching a
  table underestimate as an overrun), and verifies the exact byte count at the end
  — every failure surfaces as a `DirectWriteIntegrityError` so a wrong size table
  is a clean failure, never silent corruption. Direct-write is single-end only
  (gzip, paired, and reference reject a region for now).
- Single-end NUMA decode now uses the direct-write path end-to-end: the driver
  pre-sizes ONE shared output file (`ftruncate` to the table's total), spawns
  workers that seek to their chunk-range's byte offset and write disjoint regions
  in place, then atomically renames the completed temp to the final output — no
  per-worker part files and no concatenation pass. (Single-end gzip, paired, and
  reference still use the part-file assembly path.) On a worker pin failure (exit
  3) or a direct-write integrity failure (exit 4), `--numa auto` removes the temp
  (no output published) and transparently falls back to a correct in-process
  decode; `--numa N` (fixed) aborts with an error. The capacity gate is now
  filesystem-aware: direct-write needs only ~1× the decoded size on the output
  filesystem (plus a corrupt-table sanity cap), the part-file path needs ~2× when
  the working dir and output share one filesystem, or ~1× on each when they differ.
- `qz merge` subcommand and the public library entry points
  `qz_lib::compression::merge_single_end_archives` (streams to any `Write`, e.g.
  stdout) and `merge_single_end_archives_to_path` (atomic temp+rename to a file
  path, reject-existing-unless-`--force`): stitch N single-end v5 archives into one
  valid single-end archive (lossless decode, rebuilt `ChunkDecodedSizes` table so
  NUMA direct-write decode still applies). First building block of NUMA-aware
  compression. Single-end only; rejects paired/reference/incompatible-config inputs.

### Removal: delete the QualityCtx quality codec (fqzcomp is the only context-modeled backend)

The `quality_ctx` (context-adaptive range-coded) quality codec is **removed entirely** —
encoder selection, the `quality_ctx.rs` module (~1.5k lines), the coupled decode writer
(`bounded_write_records_quality_ctx`) and its producer (`bounded_quality_ctx_producer`),
and the `QualityCompressor::QualityCtx` enum variant. fqzcomp has been the production
default for single-end, paired, and reference since the fqz-everywhere migration; QualityCtx
was no longer selectable from the CLI (and `Auto` never chose it), so this only removes a
vestigial code path. This also removes the **last coupled (non-unified) writer** from the
decode engine, leaving every mode on the shared bounded-streaming path.

- **Breaking (decode):** codec id `4` is now *reserved-and-rejected* (mirroring the removed
  zstd/OpenZL codes). Any archive carrying a code-4 quality stream is rejected at decode/verify
  with a "quality_ctx quality backend has been removed (code 4) — recompress" hint. qz has no
  released archives, so no field data is affected.
- `resolve_quality_compressor` drops its now-purposeless `fqz_eligible` parameter (its only
  job was choosing fqz vs quality_ctx); it resolves lossless+at-scale quality to fqzcomp.
- The `--quality-ctx-block-size` flag / `quality_ctx_block_size` config field are **retained**
  (they now drive the fqzcomp quality sub-block size); only the name is legacy.
- Removed the dead `bench_quality_ctx` / `bench_quality_codecs` bench binaries and the
  `qual_fqz_probe` bz-lib example, which depended on the deleted module.

### Hardening: reference decode validates per-(chunk,mate) header/quality record counts

`validate_reference_directory` now checks that each per-(chunk,mate) header role
(R1Headers / R2HeadersIndep / R2HeaderDelta) and quality role (R1Qual / R2Qual) carries
`record_count == N_chunk` (the MappedFlags count), closing a hostile-input gap where a forged
directory could silently **misalign** headers/qualities against the per-chunk reconstructed
sequences. Single-end and paired already get this from their uniform per-chunk record_count;
reference's per-role counts differ legitimately for the M-streams, so the per-read roles are
now checked explicitly. The removed code-4 quality codec is also no longer accepted by
`RefRole::codec_valid`. Zero effect on legitimate archives (the encoder writes equal per-chunk
counts by construction); 4 new reject tests added.

### Refactor: unify reference decode onto the single-end block-streaming engine

Reference (`--reference`, archive_type 2) decompress now reuses the **shared
block-streaming decode engine** for headers + qualities (`spawn_stream_producer` +
`StreamCursor` over flat per-mate `build_role_index_mate` indices) and the
batch-parallel FASTQ assembler (`format_batch_parallel`), the same primitives
single-end/paired use. The bespoke chunk-parallel orchestrator (`Arc<RefReader>`
workers + in-order BTreeMap writer + `QZ_REF_DECODE_INFLIGHT`) is removed. Sequences
remain the one justified reference-specific divergence — reconstructed chunk-by-chunk
from the chunk-granular columns (MappedFlags/Positions/Strands/ReadLen/Edits/
FallbackPool, which can't be flat self-delimiting streams without changing the
format/ratio) — but the per-read reconstruction is now `par_iter`'d
(`reconstruct_mate`; the sequential ReadLen varint run + FallbackPool reader are
drained in serial pre-passes first, all untrusted-input validation preserved).

**Measured tradeoff** (2M novaseq 2×151 pairs vs chr20, 4 chunks, same-load
interleaved A/B): ~6% slower wall (8.33 vs 7.89 s) at **~half the peak RSS** (1.9 vs
3.6 GB) and more cores (5.6 vs 3.5). Root-caused via `perf stat`: instructions +2.9%,
cache-misses ~equal, but **IPC 1.50 → 1.33** — the engine's fine-grained
producer/cursor/batch coordination isn't amortised when reconstruction is light (it
*is* for single/paired, where heavy BSC decode dominates), and sequential-per-chunk
processing loses the old orchestrator's 4-way cross-chunk overlap. Decode is
memory-bound (wall is flat from 4→72 rayon threads), so the extra producer cores
can't buy back the wall. **Kept** for one-tool unification + the RAM win; the 6% is
intrinsic to the architecture swap. Byte-identical roundtrip; all reference/paired/
single/lib tests pass. `QZ_REF_DECODE_INFLIGHT` no longer exists.

### Performance: unify paired decode onto the single-end block-streaming engine

Paired-end decompress now reuses the **single-end block-streaming decode engine**
(per-stream block-parallel producers → bounded channels → cursors → batch-parallel
record assembly) instead of its own serial chunk-materialising decoder. The new
`bounded_write_records_paired` is a thin two-sink adapter over the *same* primitives
single-end uses (`spawn_stream_producer`, `StreamCursor`, `assemble_fastq_record`,
`format_batch_parallel`), running R1+R2 stream sets concurrently and emitting two
record-aligned outputs. R2 headers are decoded either as an independent columnar
stream or, via a streaming `R2DeltaCursor`, as a `HeaderDelta` op stream applied
per-record against the lockstep R1 id (`header_delta::apply_one_op`).

Measured on 3M novaseq 2×151 pairs, T=72 (same box/data): paired decompress
**65.2 s → ~6.9 s (~9.4×)**, CPU **1.2 → ~14 cores** — **single-end parity** (single-end
decode is ~12 cores), up from the 4.5 cores of the prior block-parallel-only fix.
Byte-identical roundtrip; peak RSS ~4–5 GB, block-bounded (constant in read count).

Two cases retain the legacy serial chunk path (byte-identical), the same way
single-end keeps a separate coupled writer for QualityCtx:
1. **QualityCtx** quality — coupled (each blob needs its chunk's sequences buffered).
2. **Mixed per-chunk R2-header codec** — R2 headers are chosen per chunk (columnar
   vs delta) for ratio; a stream whose codec varies by chunk can't be expressed by
   the single-codec-per-stream engine. Real Illumina paired is uniform (R1/R2 differ
   only by the mate field → delta always wins → all chunks delta), so production
   streams; only a mixed (e.g. tiny/synthetic) archive takes the legacy path.

`decode_chunk_v5` is retained for those two fallbacks. Full suites pass: paired (21
integration + 5 unit), reference (42), single-end (94), qz-lib lib (282).

Hardening (from an adversarial multi-agent review of the convergence; honest archives
unaffected, byte-identical): the streaming paired path now (1) rejects a forged
non-zero `encoding_type` (paired is always 0/Raw — a forged 10/4 would mis-set the
block-size DoS cap or strip a leading byte from every sequence); (2) consumes
per-(role,mate) directory entries in **chunk order** (sorted by `chunk_index`) so a
reordered footer can't silently mispair records — parity with the single-end engine's
chunk-order guard; (3) caps the concatenated R2 `HeaderDelta` op stream at a
content-derived bound (no unbounded allocation on a forged delta stream); and (4)
re-checks full op-stream consumption (restores the legacy `o != ops.len()` trailing-
bytes guard).

### Performance: parallelize paired/reference BSC-block decode (`decode_bsc_stream`)

Paired-end decompress ran the per-stream BSC block loop in `decode_bsc_stream`
**serially** — one `bsc::decompress` per block on a single core — so paired decode
used only ~1.2 of 72 CPUs and was ~6× slower than it should be (benchmarked: 10.8M
pairs decompressed in ~230 s vs single-end's 14.6 s for 10M reads at ~12 cores). The
single-end path already streams blocks through a bounded parallel producer; paired
(and the reference per-chunk-mate streams, which share this function) did not. The
sequence stream (75% of the data) and columnar header columns (19%) both decode
through `decode_bsc_stream`, so the whole bulk decoded on one core.

`decode_bsc_stream` now BSC-decodes a stream's blocks **in parallel**, batched at
`MAX_PARALLEL_DECOMPRESS` (16) blocks in flight — the same in-flight bound the
single-end decoder uses (each block is ≤64 MiB), with the `max_total`
decompression-bomb cap re-checked between batches and block order preserved
(`collect()`), so output is **byte-identical**. Single-block streams keep the inline
fast path. Measured paired decompress **65.2 s → 17.9 s on 3M pairs (~3.6×)**, CPU
1.2 → 4.5 cores, roundtrip byte-identical; peak RSS ~3.5 GB. All paired + reference
integration tests pass. (Further headroom remains in cross-stream and cross-chunk
parallelism — not addressed here.)

### Performance: skip the dead header byte-stream build on the columnar header path

`records_to_streams` (the per-chunk serial prologue shared by single-end and
paired) unconditionally built a flat `header_stream` — a `~n*64` allocation plus a
full copy of every header byte plus a varint length per record. But the production
default header codec is **columnar**, which consumes header bytes directly from the
`FastqRecord` slices and ignores `header_stream`/`header_offsets` entirely; only the
non-default BSC header path reads them. So on every default/paired compress this was
pure dead work (~100 MB of allocate-write-discard per 2.5 M-record chunk, in the
serial prologue the BSC streams can't hide). `records_to_streams` now takes a
`build_header_stream` flag (`!use_columnar_headers`) and skips the build when the
columnar codec is active. **Byte-identical output**; measured **−4.5% compress**
(10 M reads, 150 bp WGS: 11.83 s → 11.30 s mean over 4 interleaved A/B rounds, all 4
faster) plus lower peak allocation. Full integration/paired/reference roundtrip
suites pass.

### Security: extended the libbsc decode-overflow patch to the LZP stage

A follow-up code review found the original `QZ-SECURITY-PATCH (bsc-decode-output-bound)`
bounded only libbsc's QLFC entropy decoder — **not** the LZP reconstruct stage that
runs after it. `bsc_lzp_decode_block` copied an attacker-controlled match run length
(`len`, accumulated from `0xFE` count bytes in the compressed stream) into the output
with no capacity bound: `while (output < output + len) *output++ = *reference++;`. The
LZP path is live in qz (columnar headers and fqzcomp sort-keys go through the adaptive
LZP path; decode is content-agnostic, so a crafted block can select LZP mode and forge
a large length, recomputing the qz frame CRC32 + libbsc adler32 to stay
self-consistent) — a heap buffer overflow on a hostile/untrusted archive.

The patch now threads the output capacity into `bsc_lzp_decompress` /
`bsc_lzp_decode_block` and bounds **every** output write (match copies, the
speculative fast-path store, and the tail loop) to the buffer, rejecting an
overrunning block with `LIBBSC_DATA_CORRUPT`; the multi-block path also rejects an
out-of-range per-block size table or forged block count, and the match count-byte
walk is bounded to `inputEnd` (closing a paired OOB read). `build.rs` now requires the
`QZ-SECURITY-PATCH` sentinel in **both** `coder/coder.cpp` and `lzp/lzp.cpp`. Valid
data is unaffected: byte-identical roundtrip + deep verify across single-end (default +
`--ultra`), paired, and reference, with multi-block LZP exercised. See
`docs/security-libbsc-decode-overflow.md`.

### Fixed: reject explicit `fqzcomp` quality compressor with lossy / `--no-quality`

`validate_compress_config` guarded explicit `quality-ctx` against lossy quality modes
but had no symmetric guard for explicit `fqzcomp`. The fqz encode path forces
`stream_quality_binning = None` (it has no lossy/binned path), so an explicit
`fqzcomp` + a lossy mode would **silently store full-resolution qualities**, ignoring
the requested reduction (and `fqzcomp` + `--no-quality` silently dropped the codec
choice). Both are now rejected. Unreachable from the CLI/Python (which leave the
compressor `Auto`, and `Auto` resolves lossy / no-quality to BSC); this guards direct
library callers.

### Hardening: panic-safety + allocation bounds on decode paths

- The single-end/paired bounded-decode `rayon::spawn` closures now wrap their work in
  `catch_unwind` and convert a panic into a normal decode `Err`, matching the
  reference decode path. Without it, a panic escaping a spawned closure hits rayon's
  `AbortIfPanic` → `process::abort`, defeating `panic="unwind"` and PyO3's
  panic→exception conversion (an abort the Python binding cannot catch).
- `reference::positions::decode_positions_mate1/2` clamp their `Vec::with_capacity`
  to the payload size (every record consumes ≥1 byte), bounding the speculative
  allocation if a caller passes an unvalidated record count — matching the sibling
  untrusted-count readers.
- Doc-drift/clarity: corrected stale `code 3`→`code 6` and dead `can_stream_decompress_path`
  → `require_streamable_v5` references in the decode comments, and a now-misleading
  rayon-panic recv-error message; the `refmeta` `MAX_NUM_REFS` cap is split into a
  unit-tested helper so its boundary stays covered.

### Changed: unified the codec-byte namespace (front-header now uses `codec_ids`)

The v5 front header and the footer directory previously used two different codec
numberings for the same codecs — the only real divergence being fqzcomp quality
(`3` in the front header's legacy `compressor_to_code` numbering vs `6` in the
directory's `codec_ids`). The front header now writes the `codec_ids` values, so a
given codec byte means one codec everywhere (front header AND directory, across all
archive types). The `front_codec` module and the `compressor_to_code` /
`seq_compressor_to_code` / `header_compressor_to_code` / `code_to_compressor` /
`code_to_header_compressor` functions are deleted; encode/decode now go through the
`codec_ids` forward/reverse mappers (`codec_id_to_quality_compressor` /
`codec_id_to_header_compressor` added).

**Breaking (no field archives → no practical impact):** the only on-disk change is a
single front-header byte — single-end fqzcomp archives now stamp quality codec `6`
instead of `3`. Single-end archives written before this change won't decode with the
new binary (and vice versa); re-compress them. Paired/reference front headers carry a
`Bsc` placeholder (their authoritative per-mate codec lives in the directory), so
their bytes (and all BSC/quality_ctx single-end archives) are byte-identical. Lossless
roundtrip + deep verify confirmed across single/ultra/paired/reference.

### Cleanup: post-purge dead code, dead deps, doc-drift, and naming hazards removed

Follow-up housekeeping after the v4/OpenZL legacy purge. No on-disk format change —
single-end, ultra, paired, and reference v5 archives are byte-identical to before.

- **Dead dependencies:** dropped the unused `openzl-sys` and `constriction`
  workspace dependencies (and `constriction.workspace = true` from qz-lib), trimming
  the build's dependency tree. Updated the stale `panic = "unwind"` profile comment
  (no more "OpenZL graph fns" example) and the `constriction` mention in
  `quality_ctx.rs`.
- **Dead code:** removed the three test-only `can_stream()` codec methods and their
  parity test (the live `require_streamable_v5` gate is the single source of truth),
  `code_to_seq_compressor`, the `compress_qualities_fqzcomp` record wrapper (tests
  migrated to `compress_qualities_fqzcomp_quals`), the orphaned v4 sequential
  `build_stream_index` walker (its live-producer setup tests now build the
  `StreamIndex` directly), the always-zero `VerifyResult.streams_skipped` field and
  the dead CLI `PARTIAL`/`Skipped` branches that read it, the always-zero
  `role_bytes(StreamRole::Nmask)` addend, and the five vestigial v4 `FixedHeader`
  body fields (the 12-byte v5 body is preserved as reserved zero bytes).
- **Naming:** renamed the live per-block frame helpers `write_v4_block_header` →
  `write_block_frame_header`, `read_v4_block_header` → `read_block_frame_header`,
  `verify_v4_block_crc` → `verify_block_frame_crc`, `compute_v4_crc` →
  `compute_block_frame_crc`, and `struct V4BlockHeader` → `BlockFrameHeader`, so "v4"
  no longer names live code (it now means only the removed/rejected stream-major
  format). Pure rename, zero on-disk change.
- **Dedup:** introduced named front-header codec constants (`front_codec::{QUAL_BSC,
  QUAL_FQZ, QUAL_CTX, SEQ_BSC, HDR_BSC, HDR_COLUMNAR}`) replacing raw integer
  literals in the single-end decode dispatch and `require_streamable_v5`, and made
  `paired::streams::MAX_BLOCK` a re-export of `decompress_impl::BSC_MAX_BLOCK_LEN` so
  the per-block size cap has one definition.
- **Doc-drift:** fixed stale comments referencing removed OpenZL/`qz-v3`/v4 machinery
  and deleted functions; corrected `cli.rs` docs (dropped "reordering", fixed the
  `Auto` quality default to Fqzcomp and `ultra` to "1-3, 0=auto"); dropped the
  `encoding_type_display_name` `5 => "sequence-delta"` arm. Indented four doc-list
  continuations to clear clippy.

### Fixed: unchecked length-varint overflow could panic reference decode on a hostile archive

Two reference-mode decode sites parsed varint-length-prefixed records out of
BSC-decompressed (attacker-controlled) block bytes with a plain `offset + length`
add before the bounds check: `FallbackPoolReader::next_literal`
(`reference/fallback.rs`) and `read_record_block_stream` (`reference/backing.rs`, on
the global IntervalMap + ReferenceMeta decode reached at the start of every reference
decompress/deep-verify). Because the release profile has no `overflow-checks`, a
crafted archive whose record length varint decodes to `u64::MAX` made the add wrap,
bypassing the `> len()` guard, so the subsequent slice panicked (a DoS — caught by
Rust's bounds check, not memory unsafety). Both now use the same
`checked_add().filter(|&e| e <= len)` guard already used by every `StreamCursor`
reader and the paired path, rejecting the input with a clean error. Added two
hostile-input regression tests. No change to valid-archive behavior (identical `end`
for in-bounds lengths).

### Removed (breaking): the entire v4 stream-major archive format and the OpenZL codec

qz now has a **single archive architecture: v5 chunk-major**. The v4 stream-major
format — both its encoder (`compress_in_memory`) and its decoder
(`decompress_v4_legacy` plus the v4 fallthrough branches inside the streaming
decode/verify paths) — has been deleted. There were no field archives, so v4
archives never need to be read again.

- **Breaking — decode:** v4/legacy archives (any non-v5 version byte) are now
  rejected at decompress/verify with a "recompress with current qz" hint instead of
  being decoded. Re-compress any old archive.
- **Breaking — config:** the experimental `AdvancedOptions` fields that only ever
  routed to the v4 in-memory encoder were removed: `quality_modeling`,
  `quality_delta`, `sequence_delta`, `twobit`, `header_template`, and
  `compression_level`. (`bsc_static`, `sequence_hints`, `rc_canon`, and the chunk/
  block/window levers stay — they drive the v5 chunked path.) A `--config` JSON that
  names a removed field is now ignored (serde default) rather than selecting v4.
- **Removed — OpenZL:** the `openzl` module, the `openzl-sys` dependency, and the
  `OpenZl` variants of `QualityCompressor`/`SequenceCompressor`/`HeaderCompressor`
  (reachable only via the v4 path) are gone. The wire codec values they used
  (front-header code `2`, directory `CODEC_OPENZL`) are now reserved and rejected
  with a hint.
- **Removed — codecs/modules:** `EncodingType::SequenceDelta` (encoding_type 5),
  the `quality_model`, `quality_delta`, `n_mask`, and `read_id` modules, and the
  v4-only codec helpers (2-bit sequences, OpenZL/template/model/delta quality and
  header encoders). The `experimental` Cargo feature (now vestigial) was removed.

This is a v5-only-output change with **no effect on v5 archives** — v5 encode is
byte-identical and the full v5 single-end/ultra/paired/reference decode + verify
paths are unchanged. The shared per-block framing (`compress_stream_to_bsc_blocks_v4`
/ `write_v4_block_header`) keeps its historical name but is the v5 chunk framing, not
the removed v4 archive layout.

### Changed: unified the v5 directory `codec` byte to one namespace across all modes

The per-entry `codec` byte in the v5 footer directory (`ChunkDirEntry.codec`) now
uses the single `codec_ids` namespace for **every** archive type. Previously the
field was overloaded: single-end archives stamped the front-header
`compressor_to_code` numbering (fqzcomp quality = `3`) while paired and reference
archives stamped `codec_ids` (where `3` = columnar, fqzcomp = `6`), so the same byte
meant different codecs depending on `archive_type` — a footgun guarded only by a
prose comment. Single-end encode now maps its compressor enums through
`codec_ids::{quality,sequence,header}_compressor_to_codec_id`, so a given directory
`codec` value denotes one codec regardless of mode. The single-end **front-header**
codec bytes are unchanged (single-end decode still reads its codecs from there; the
directory byte is advisory — uniformity-checked, never dispatched on). The only
on-disk change is the single-end directory fqzcomp-quality codec byte (`3` → `6`);
decode is unaffected and cross-decodes both ways (old and new archives decode under
either build, verified lossless). No ratio or behavior change.

### Fixed: columnar header blocks over 64 MiB produced a self-undecodable archive

A chunk's columnar-encoded read-header blob was written as a single uncapped block,
but the decoder rejects any block larger than the 64 MiB `MAX_BLOCK` cap. So a chunk
whose columnar header blob exceeded 64 MiB (reachable only at large `chunk_records`
with long, high-entropy read IDs) produced an archive that **failed its own
decode/verify** — full loss of that input. The columnar header blob is now cap-split
into `≤MAX_BLOCK` self-contained sub-blocks on encode (a new
`columnar_blocks_capped`, mirroring the existing `fqz_blocks_capped`), and the
paired/reference decoder (`decode_columnar_headers`) decodes the resulting one-or-more
blocks and concatenates them. The common case — a blob that fits under the cap — still
emits a single block **byte-identical to before**, so existing archives are unchanged
and there is no ratio cost. Covers the shared single-end, paired, and reference
columnar-header path.

### Fixed: reference mode now honors `chunk_records` from the config

Reference (`--reference`) compress read its chunk size only from the process-global
`QZ_REF_CHUNK` env var and ignored `CompressConfig::advanced.chunk_records`, so a
caller (or test) setting that field had no effect. Reference now resolves the chunk
size with the precedence **explicit `chunk_records` → `QZ_REF_CHUNK` → the tuned
`REF_CHUNK_RECORDS` (500K) default**; production runs (which leave `chunk_records` at
the shared default) are unchanged. Because an explicit config value now wins over the
env var, the former test-side `RefChunkEnv` mutex/guard — which existed only to
serialize process-global `QZ_REF_CHUNK` mutation — is gone, closing a latent
`setenv`/`getenv` data race.

### Reference decode is now bounded-RAM streaming (constant in read count)

Reference (`--reference`) decode and deep verify now run as a bounded-RAM
streaming pipeline whose peak RSS is **constant in read count** — it no longer
scales with the number of chunks/reads in the archive. The work landed in five
parts, in implementation order:

**1. Fallback pool is now per-(chunk,mate), not one global stream.** The
fallback-literal pool — the off-reference reads stored verbatim — moved from a
single whole-archive **global** `FallbackPool` stream to a **per-(chunk,mate)**
`FallbackPool` role written inline into each chunk's payload (omitted when that
mate has no fallbacks, i.e. when `N == M`). This is the structural groundwork for
streaming decode: each mate's literals are read on demand from its own entry and
decoded one BSC block at a time, instead of threading one global
`FallbackPoolReader` across the whole archive.

- **Encoder**: no global fallback spool any more — reference compress touches no
  scratch files at all (the `op_dir`/`TmpDirCleanup` temp directory and the spool
  were its only users, both removed). Each `(chunk,mate)` writes its own pool block
  stream inline; the per-literal 8 MiB cap (`PER_RECORD_MAX_BYTES`) now lives in
  `write_fallback_pool`.
- **Validator**: `GLOBAL_ENTRY_CAP` drops 5 → 4 (`FallbackPool` is no longer
  global). The ownership contract is enforced per `(chunk,mate)`: a `FallbackPool`
  entry is present **iff** `N − M > 0`, its `record_count` must equal `N − M`, and
  it is never present-but-empty — every mismatch in the attacker-controllable
  footer is rejected before any decode work. A `FallbackPool` entry carrying the
  global sentinel is now rejected.
- **Format revision**: the `ReferenceMeta` `meta_version` byte is bumped 1 → 2 and
  doubles as the format-revision discriminator. Archives written by an older qz
  (global-pool layout, `meta_version` 1) are **rejected at decode** with a
  recompress hint rather than silently mis-decoded.
- **Ratio impact measured ≤0.5% on realistic data**: splitting the literals into
  per-chunk pools fragments BSC slightly, but the regression is bounded by two
  committed gates — pool-relative fragmentation ~+9.5% at a realistic K=8 chunking
  (Test A) × the pool's ~2.8% share of a realistic ~1%-fallback archive (Test B)
  ⇒ composed whole-archive regression ≈0.27% (≤0.5%). The default pool target stays
  `META_TARGET` (25 MiB).

**2. Output is streamed through two write sinks.** Decompress and deep verify no
longer build full-FASTQ `r1_out`/`r2_out` Vecs. The reconstruction emits each
chunk's reads through `&mut impl Write` sinks, so peak RAM is bounded by one
chunk's working set plus the whole-archive globals — constant in decoded output
size instead of O(decoded output).

- **Decompress**: each `_R1.fastq`/`_R2.fastq` is written straight into a
  `BufWriter` over an `O_EXCL` temp file; the atomic two-file publish + drop-guard
  are unchanged. Only the produced FASTQ is streamed — never materialised.
- **Deep verify** drives the reconstruction through two per-mate `CrcWriter` sinks
  (no full-FASTQ Vec) and reports **per-mate CRC32s**: `qz verify` prints
  `R1 CRC32:` / `R2 CRC32:` for reference archives (informational — not compared
  against a stored value). Single-end keeps its single `CRC32:` line; paired deep
  verify still folds both mates into one combined CRC; fast verify is unchanged
  (per-block CRC walk, no per-mate split).

**3. The archive is read seekably, with resident globals.** `open_ref` returns a
`RefReader { File, len }`. The whole-archive globals
(`PackedBacking`/`IntervalMap`/`ReferenceMeta`) load **once** into a resident
`Globals`; each per-chunk directory entry is then read on demand via a positioned
`read_exact_at` (`RefReader::read_entry`, `&self`, no shared cursor, bounds-checked
against the once-captured file length). The archive is no longer materialised in
memory.

**4. fqzcomp-decode leak fixed (shared single-end + paired fix).** Decode peak RSS
was still climbing ~linearly with chunk count — 100K pairs ≈ 358 MB, 600K pairs ≈
1.52 GB (**4.35× RSS for 6× reads**) — despite the streaming, seekable path. Root
cause (valgrind massif): htscodecs' `fqz_create_models` allocates a ~25 MB qual
model (`SIMPLE_MODEL[CTX_SIZE]`, `CTX_SIZE = 1<<16`) through a **per-thread TLS
arena** (`htscodecs_tls_alloc`) that returns the buffer to a per-thread pool on
destroy and **never reclaims it to the OS** — it lives until the owning thread
exits. Decoding fqz inside a rayon `par_iter` spread the model across
O(workers-ever-touched) threads, retained for the whole process.

- Fixed by compiling htscodecs `utils.c` with `-DNO_THREADS` (build.rs,
  sentinel-guarded), which switches the TLS arena to plain `calloc`/`free` so each
  model is freed promptly. The arena's reuse-within-a-thread saving did not survive
  QZ's per-call rayon thread hops anyway. The C `calloc`/`free` bypasses Rust's
  mimalloc, and there is no decode-speed regression.
- Plus a residual: spreading each tiny per-chunk-mate quality blob across the whole
  pool warmed extra worker threads (each pays a one-time, never-reclaimed
  per-thread arena/stack RSS cost). Fixed by decoding single-block fqzcomp blobs
  **inline** (no `par_iter`/`rayon::join` fan-out) in `decode_qual`
  (`paired/mod.rs`) and `decompress_qualities_fqzcomp` (`codecs.rs`); the genuinely
  multi-block single-end path (> 500K reads/blob) still parallelises.
- **Shared fix** — single-end and paired fqzcomp decode had the same latent leak.

**5. Decode is bounded-parallel with an in-order writer.** Chunks reconstruct on a
rayon worker pool over `Arc`-shared read-only `RefReader`/`Globals`/
`ChunkDirectory` (each worker reads its per-chunk entries on demand via the
positioned `read_entry`).

- **Byte-identical output preserved.** An in-order writer (`completed` BTreeMap +
  `next_write`) flushes chunk *i*'s R1/R2 bytes before chunk *i+1*'s regardless of
  completion order or thread count, so the reconstructed FASTQ — and the
  deep-verify per-mate CRC32s — are byte-for-byte the sequential decode.
  `max_inflight == 1` is exactly sequential. Verified deterministic across repeated
  parallel runs and byte-identical on real multi-chunk data (HG002 chr20 subset).
- Workers `catch_unwind` and convert panics to `Err` (a panic escaping a detached
  task would abort the process); the orchestrator is first-error-wins (bails on the
  first worker `Err`).
- **Concurrency is bounded by a small fixed cap.** `max_inflight` defaults to
  `REF_DECODE_MAX_INFLIGHT = 4` — the MEASURED decode-speed knee (an inflight sweep
  on a real ~40-chunk chr20 archive flattens past 4 because QZ decode is
  memory-bandwidth-bound, so > 4 buys little). It is further bounded by a memory
  budget (`REF_DECODE_BYTES_PER_READ = 1024`, a realistic upper bound that only
  trims inflight below the cap for very large chunks), and is overridable via the
  `QZ_REF_DECODE_INFLIGHT` env knob (clamped `≥ 1` and to the chunk count; pin `1`
  for a sequential oracle). Peak RAM ≈ globals + ~2 × `max_inflight` × one chunk's
  working set — **constant in read count once `num_chunks ≥ max_inflight`** (the
  concurrency plateau).

**Net result.** The `#[ignore]` oracle
`reference_decode_peak_rss_constant_in_read_count` pins `QZ_REF_DECODE_INFLIGHT`
and compares SMALL vs LARGE inputs that both sit past the cap, verifying decode
peak RSS is constant in read count at the concurrency plateau. Final measured
ratio is **≈1.0–1.1× for 6× reads** (vs the pre-leak-fix 4.35×); the residual is
mimalloc high-water / box-load jitter, not growth, and the oracle asserts under
`FACTOR = 1.6`. Lossless roundtrip is preserved throughout (reference R1+R2,
single-end, and paired all byte-identical; real-data chr20 NovaSeq subset cmp'd
byte-identical). A new `qz-bench` bin `bench_ref_decode_mem` (system allocator,
valgrind-ready) reproduces the property. No wire-format change from parts 2–5;
the only format change is part 1's per-(chunk,mate) pool + `meta_version` 2.

**Breaking:** reference archives produced before this change no longer decode —
recompress them.

### Reference-based archives now use the v5 chunk-major container

Reference-based (`--reference`) FASTQ compression now produces the **v5
chunk-major unified container** (archive version 5, `archive_type = 2`) instead of
the bespoke v4 reference front-header + back-patched footer. Encode, decode,
`verify` (deep and `--fast`), and the public decompress/verify dispatch all route
reference archives through the unified directory. The 12 reference-specific stream
roles plus the five global streams map onto the unified `StreamRole`/`ChunkDirEntry`
directory (with `mate`-tagged per-chunk entries and `GLOBAL_SENTINEL` globals), and
a reference-mode footer validator enforces the §9 structural contract. Benefits,
matching single-end/paired v5:

- **Seek-free**: the directory footer + 20-byte locator live at EOF, so the
  encoder never seeks (no `footer_offset` back-patch). Reference output still
  requires a named output prefix (it reconstructs two FASTQ files).
- **Hostile-input hardened**: untrusted footer record counts are clamped with
  `usize::try_from`; the consensus/bitmap block reader validates its block framing
  in a header-only pre-pass (Σ-count + per-block decode cap) **before** the
  reassembly allocation, then reserves fallibly — a crafted count is rejected with
  an error ahead of any large allocation instead of aborting; and the substitution
  edit decoder caps its capacity hints by the input-stream lengths and uses checked
  `total * 3` arithmetic, so a single huge edit-count varint can no longer
  overflow/OOM. (The remaining materialize-all decode is the documented bounded-RAM
  follow-up.)
- `qz verify` (deep and `--fast`) handles v5 reference archives via the unified
  directory.

**Breaking:** a fresh build no longer produces the v4 reference format, and old v4
reference archives (marker `0xF2`, version 4) are no longer decodable — recompress
them (the dispatch rejects legacy `0xF0`/`0xF1`/`0xF2` markers with a recompress
hint). The dead v4 reference footer/front-header machinery (`RefFooter`/`RefEntry`,
the v4 front-header writers, `archive_kind`/`ArchiveKind`, and the legacy
`FLAG_PAIRED`/`PAIRED_MARKER`/`REFERENCE_SUBFORMAT`/`PAIRED_REFERENCE_*` constants)
has been removed; the public `is_paired_archive`/`is_reference_archive` predicates
are now pure v5 `archive_type` checks. Lossless roundtrip (R1+R2 byte-identical)
verified on synthetic fixtures and a 50K-pair real chr20 NovaSeq dataset, plus
deep + `--fast` verify.

### Paired-end archives now use the v5 chunk-major container

Paired-end FASTQ compression now produces the **v5 chunk-major unified container**
(archive version 5, `archive_type = 1`) instead of the legacy v4 paired format
(`0xF0` marker). Both encode and decode production dispatch route paired archives
through the unified path. Benefits, matching single-end v5:

- **Seek-free, pipeable to stdout**: `qz compress -i R1.fastq -i R2.fastq -o -`
  now streams a paired archive to stdout (the directory footer + 20-byte locator
  live at EOF, so the encoder never seeks). Reference-mode paired still requires a
  named output file.
- **Bounded-RAM decode**: the decoder reconstructs one chunk at a time from the
  footer directory; peak decode memory is bounded by chunk size, constant in read
  count.
- `qz verify` (deep and `--fast`) handles v5 paired archives via the unified
  directory (per-block CRC walk for `--fast`, full reconstruct + CRC32 for deep).

**Breaking:** a fresh build no longer produces the v4 paired format, and the dead
v4 paired encoder/decoder/verifier (footer, directory entry, `Role`, writer, and
chunk decoder) has been removed. qz has no field archives, so v4 paired archives
never existed in the wild; the dispatch still detects the legacy paired-archive
marker byte and now returns a clear "recompress with the current qz" error
instead of decoding it.
Lossless roundtrip (R1+R2 byte-identical) verified through the public dispatch on
v5.

### Paired-end quality now uses fqzcomp (ratio parity)

Paired-end FASTQ now selects **fqzcomp** for quality under `Auto` (≥100k reads/mate),
matching single-end and reference mode; below that threshold it stays on BSC,
unchanged. The change is encode-side wiring only — paired already shared
`compress_one_chunk` (which emits fqzcomp blocks) and the paired decoder already
handled the fqzcomp directory codec (`6`); this flips the resolution/params switches
and widens the footer codec validator to accept `6` for quality roles. Also fixes a
latent bug where an explicit `--quality-compressor fqzcomp` in paired mode was
silently downgraded to BSC (paired now accepts `Auto` or explicit `fqzcomp`).

Measured on 200k read pairs of HG002 chr20 2×151 (4-bin quality): total archive
16,951,238 → 16,678,912 bytes, **−1.6% (−0.26 MB), entirely from the quality stream**
(sequence/headers unchanged). The win is modest here because 4-bin quality is already
highly compressible; richer quality distributions gain more. Lossless roundtrip
verified byte-identical (R1+R2). Old paired archives (quality codec 1/4) still decode;
`quality_ctx` remains a selectable decoder and the codec for those archives.

### Fixed
- `--ultra` could intermittently write a **corrupt, non-decodable archive** on
  small, highly repetitive inputs (compress reported success but decompress
  failed with a BSC error — a data-loss hazard if the source was deleted).
  libbsc's OpenMP multithreaded BWT was selected by the configured block-size
  *target* (188 MiB) rather than the *actual* block size, so a tiny degenerate
  block was routed through the MT path, which has a parallel suffix-sort race on
  near-equal suffixes (reproduced ≤640 KiB; ≥1 MiB and large blocks unaffected;
  single-thread BSC is always correct). Small actual blocks (< 16 MiB) now use
  the reliable single-thread adaptive coder; large blocks keep the MT speedup.
  Pre-existing — reproduced on the prior v4 ultra build.

### Changed

- **Archive format: chunk-major layout (version 5) for single-end and `--ultra`
  FASTQ.** Single-end and `--ultra` archives now use archive **version 5**
  (chunk-major) instead of version 4 (stream-major). Layout: a front header
  followed by interleaved per-chunk blocks, a directory footer, and a 20-byte
  trailing locator (`footer_len u64 | footer_crc32 u32 | "QZFOOTR1" 8B`). The
  locator is at the end of the file so the encoder never seeks — output can be
  piped to stdout or any non-seekable sink. The directory footer maps every
  block to its role (Headers/Sequence/Nmask/Qual/RcFlags), chunk index, codec,
  absolute frame offset/length, and record count; on decode a `DecodePlan` is
  built from it and feeds the existing bounded-streaming cursor core. Each
  block retains the existing v4 per-block framing `[block_len][record_count][crc32][payload]`.
  `version` selects the **layout** (4 = stream-major, 5 = chunk-major);
  `encoding_type` continues to select sequence **semantics** (Raw=0,
  RawWithHints=4, RcCanon=6, UltraBigBlock=10). `SequenceDelta` (type 5) is
  not bounded-streamable and stays on the v4 stream-major path. Compress is
  now a **streaming writer** — each chunk's blocks flow straight to the output
  as they complete, with no scratch/spool temp files. This is a **bounded-RAM
  design by construction** (peak ≈ pipeline window × chunk + the small
  directory, constant in input size); ratio is unchanged vs v4 (+~60 bytes
  footer overhead, measured ~0.003% on a 50k-read slice). Large-scale
  flat-RSS verification on an unloaded box is pending. **v4 archives remain
  fully decodable.** Paired and reference archives are unaffected. `verify`
  and `verify --fast` both handle v5.

### Removed
- Reorder compression modes (`AdvancedOptions.reorder` / `local_reorder`, the
  `ReorderMode` enum, `ChunkedMode::Sorted`, the sort path, and the
  global-reorder bucket-sort encoder). They were library-only (never in the
  production CLI) and compressed worse than identity. **Breaking** for library
  callers that set these fields. No archive-format or CLI impact.

### Fix: `-t` / thread count now honored in paired and reference modes

The rayon global thread pool was only built inside the single-end paths
(`compress_impl::compress` / `decompress_impl::decompress`); paired- and
reference-mode invocations dispatch before those and so silently used rayon's
default (all cores), ignoring `-t`. Thread-pool init is now hoisted to a shared
`ensure_global_thread_pool` called at the top of `compress()` / `decompress()`,
so every mode honors the requested count. Verified: paired compress at `-t 4`
caps at ~360 % CPU (was effectively all-cores). (We keep the process-global pool
rather than an owned pool + `install()`: the chunked pipeline spawns
`std::thread::scope` workers that each run rayon work, and `install()`'s
thread-local pool does not propagate to `std::thread` children — only the global
pool does.)

### Single-end default quality codec is now fqzcomp (Task 3 of streaming-fqz)

For **single-end** lossless FASTQ at scale (≥ `MIN_READS_QUALITY_CTX` reads),
`Auto` quality resolution now picks **fqzcomp** instead of quality_ctx. fqzcomp
beats quality_ctx on every measured dataset and is now bounded-streaming (Tasks
1–2), so the decompress-RAM penalty that previously kept it off the default is
gone — measured end-to-end, fqzcomp-streaming decode peak RSS is within ~2 % of
quality_ctx (was +1.4 GB/+3.6 GB for the old non-streaming fqz), at comparable
decode time, with a slightly smaller archive. **Paired-end is unchanged** (keeps
quality_ctx — its pipeline has no fqzcomp path yet; that is a later phase).
`resolve_quality_compressor` gained `fqz_eligible` (single-end only) and
`modeling_or_delta` parameters; when quality modeling or delta is enabled, `Auto`
now resolves to **BSC** (model deltas need a byte backend — resolving to
fqz/quality_ctx would store a code that mismatches the written bytes). Explicit
`--quality-compressor` choices are still honored verbatim.

### fqz quality decode: bounded-streaming (Task 2 of streaming-fqz)

Fqzcomp quality archives now decode through the **bounded-streaming** path
(`QualityCompressor::Fqzcomp.can_stream()` is true) instead of materializing the
whole quality stream in memory. A new `StreamCodec::Fqz` producer decodes each
record-capped v4 block via `decompress_qualities_fqzcomp` (which already returns
the per-record `[varint(len)][raw ASCII]` stream and carries its own DoS guards),
and the independent writer emits those bytes via `QualityInput::Decoded` (raw
ASCII — fqz is not bit-packed). Per-block decode concurrency is clamped (≤4) so
peak resident memory stays bounded. The three stale "Fqzcomp must use the
in-memory path" pre-checks are removed. Decompress peak RAM for fqz archives is
now bounded by block size (constant in read count) rather than O(reads). Lossless
roundtrip verified on real 150 bp WGS and via new variable-length / multi-block
integration tests; the in-memory multi-block decoder is retained for
non-streamable archive combinations.

### fqz quality encode: record-capped multi-block format (Task 1 of streaming-fqz)

Single-end fqzcomp quality encode now emits **record-capped multi-block** output
instead of one global-sorted blob per chunk. The chunk's quality strings are
partitioned into consecutive slices of at most `quality_ctx_block_size` records
(default 500_000); each slice is fqzcomp-compressed independently (its
mean-quality sort permutation is local to that ≤cap-record slice) and framed in
the existing v4 multi-block layout with its REAL per-block record count. This is
the encode prerequisite for bounded-streaming fqz decode (Task 2 wires up the
streaming decoder). No default changes; explicit fqz archives still decode
through the existing in-memory multi-block path (`decompress_fqzcomp_blocks`,
unchanged). New record-capped splitter `fqz_record_capped_blocks` is distinct
from reference mode's byte-capped `fqz_blocks_capped` (left byte-for-byte
unchanged). Byte-identical roundtrip; new unit + integration tests.
### security/correctness: multi-agent code-review hardening pass

Addressed every confirmed finding from an adversarially-verified code review across
qz-lib, bz-lib, and the vendored strobealign. None affected the happy-path lossless
roundtrip on conformant input (verified byte-identical on real data: default, ultra,
exome, and the full reference/paired suites); the fixes harden the FFI boundary,
the untrusted-archive decode path, and error/panic paths.

**Critical — memory safety:**
- **BSC decode heap-overflow bound (libbsc).** The QLFC coder's embedded output
  length is independent of the block-header `dataSize` that sizes the decode
  buffer, so a *crafted* archive could make libbsc write past the allocation. The
  vendored libbsc decode chain (`bsc_coder_decompress` → `*_decode_block` →
  `bsc_qlfc_*_decode`) now threads the real output-buffer capacity and rejects any
  coded length that would overrun it. `build.rs` fails loudly if the patch is
  absent (libbsc is a gitignored upstream clone). Byte-identical for valid data.

**High:**
- **Quality `< 33` lossless guard** (`columnar.rs`): the default lossless quality
  path silently collapsed quality bytes below 33 to `'!'`; now rejects them like
  the high-end / `quality_ctx` guards instead of corrupting silently.
- **`--sequence-delta` multi-chunk fix**: the per-chunk delta encoder is now routed
  to the single-chunk in-memory path (excluded from `can_stream`), so its
  chunk-local back-references can't be mis-resolved against a global decode cache.
- **Reference decode popcount guard**: a crafted `mapped_flags` whose set-bit count
  exceeds the declared mapped count `M` no longer panics with an out-of-bounds
  index (untrusted-input DoS); rejected up front.
- **bz `ref_span` bound**: a crafted CIGAR of huge ref-consuming ops could drive
  `Vec::with_capacity` to an uncatchable process abort on decompress; now capped.

**Medium / robustness:**
- Reference encode: `catch_unwind` around the detached `rayon::spawn` encode so a
  panic returns a clean error instead of `process::abort()` (which skipped temp-dir
  cleanup and PyO3 exception conversion); and `drop(rx)` on the error path so a
  blocked producer can't deadlock the `thread::scope` join.
- `read_id` packed decode now bounds its pre-allocation by the bytes remaining in
  the stream (untrusted-input OOM amplification).
- FASTQ reader tolerates trailing/stray blank lines at record boundaries instead of
  aborting compression with a misleading EOF error.
- Chunked compress joins the in-flight head worker on a mid-stream read error
  (no detached BSC thread left running).
- bz fqz quality decode uses `wrapping_sub(33)` — an exact inverse of the encoder's
  `wrapping_add(33)` over the full byte range (out-of-spec qualities ≥ 223 now
  round-trip instead of failing the content CRC).
- bz `extract` warns and counts reads dropped due to duplicate-QNAME primaries
  instead of silently overwriting them.

**Low / latent hardening:**
- `arithmetic_quality` and `debruijn` decoders reject out-of-range / truncated
  input instead of corrupting or panicking (both off the production dispatch).
- strobealign: validate `.sti` `bucket_starts` values (not just length) so a corrupt
  sidecar triggers a clean rebuild; validate resolved `w_min`/`w_max` in
  `with_window`; tighten the `with_aux_len` bound (≤ 54) and make the
  `SyncmerIterator` k/s masks overflow-safe at the documented `k=32`.
- quality_ctx v2 encoder now requires exact quality/sequence length equality
  (matching v1) to foreclose a shared-range-coder desync class.

### perf: continuous-pipeline encode for reference compression, byte-identical

Replaced the reference-compress consumer's batch BARRIER (fill `window` chunks →
parallel-encode the batch → serial-write → repeat) with a CONTINUOUS pipeline:
each mapped chunk is encoded as it arrives (in-flight bounded to `window` via a
token semaphore) and results are written in chunk order through a reorder buffer.
The wait-for-the-producer now OVERLAPS the encode instead of alternating with it,
so cores don't idle during `recv`.

Profiling (chr20 / HG002 2×151, 4M-pair subset, 72 threads) showed reference
compress is parallelism/latency-bound, not codec-bound: at the default window the
encode stage is core-starved, and the old barrier idled the consumer during recv.
Measured effect — **at `--reference-window 4`: 24.2s→21.0s and 8.9GB→6.8GB peak
RSS**; that 21.0s matches the old *window=8* speed (20.7s) at **~half the memory**
(6.8 vs 12.5GB). So peak throughput is reached at a much shallower (cheaper)
window. Output is **bit-for-bit identical** (verified: same archive md5 on real
data, full reference roundtrip/determinism suite green) — the coverage/`num_pairs`
fold still runs in chunk order and writes stay strictly in chunk order.

The default `--reference-window` is **raised 2→4** to match the new sweet spot:
default reference compress is ~28% faster (29→21s on the 4M-pair bench) for ~+1.4GB
peak RSS (5.4→6.8GB). Lower it to trade speed back for memory.

Also adds an env-gated `QZ_REF_TIMING=1` diagnostic that prints producer
(read/map/diff) and consumer (work-wait/encode-backpressure/write) phase times.

### perf: faster strobealign mapper (reference-based compression), byte-identical

Two byte-identical optimizations to the vendored strobealign map-only path used
by reference-based compression, found by profiling a new standalone-map bench
(`bench_map_only`). On chr20 / HG002 2×151 (mimalloc, 72 threads): **map
throughput +9% multi-thread** (1.00M→1.09M pairs/s), **+6% single-thread**
(23.4K→24.8K pairs/s). Mapping output is provably unchanged — same mapping
digest across the change, and all reference roundtrip/determinism tests stay
green — so archives are bit-for-bit identical and lossless.

- **Anchor sort materializes its key once.** `Chainer::get_chains` sorted anchors
  with `sort_unstable_by_key(packed u128)`, which recomputes the key on every
  comparison (~n·log n times per read; profiled at ~12% of mapper self-time). It
  now packs each anchor into a `u128` once (a lossless, order-preserving bijection
  of `(ref_id, ref_start, query_start)`), sorts the PODs, dedups, and unpacks —
  identical ordering, key computed once.
- **`Nam` no longer carries its anchor `Vec`.** Nothing downstream of chaining
  read `Nam.anchors`; only its length (for scoring) and a trace line used it.
  Replaced the per-NAM `Vec<Anchor>` with a `num_anchors: usize`, removing a heap
  allocation per NAM and a Vec clone on every NAM clone during pairing — the
  larger of the two wins, especially multi-threaded (less allocator contention).

Note: this does not change reference-compress wall-clock — mapping (~11.5s for the
10.8M-pair chr20 set) already overlaps under the encode stage in the
producer/consumer pipeline. The win is for standalone mapping and
low-coverage/few-core configs where mapping is exposed.

### review: hardening round 3 follow-ups (temp-file O_EXCL parity + OpenZL cap)

Closes the two deferred items from round 3:
- **Scratch temp files** in the chunked/sort and global-reorder compress paths now
  use a process-unique token (`scratch_tmp_token`: pid + counter) instead of a
  bare pid, and are reserved with `create_new_for_write` (O_EXCL) — so two
  concurrent in-process jobs sharing a `working_dir` can't collide, and a
  pre-planted symlink can't be followed. The paired-decompress temp writers got
  the same O_EXCL parity (their names were already pid+nonce unique).
- **OpenZL** gained `decompress_parallel_max`, and the OpenZL header/sequence
  fallback decoders now pass a content-derived `stream_decode_cap(num_reads)` —
  matching the BSC fallback decoders' cumulative decompression-bomb cap.

### review: code-review hardening round 3 (qz + bz + strobealign)

A multi-reviewer pass over the whole workspace turned up a cluster of latent
correctness/safety bugs — none reachable from the `qz`/`bz` CLIs, but live traps
for library-API callers and for untrusted-input handling. All fixed with
regression tests; suites green (qz 84 + unit, bz 34+24, strobealign 34), clippy
`-D warnings` clean.

**Silent data corruption / loss (library-API configs):**
- **`rc_canon` + `sequence_hints`/`sequence_delta`** was an accepted combo that
  emitted a per-record hint/delta byte the type-6 (`RcCanon`) decoder never
  strips — corrupting every sequence. Now rejected in `validate_compress_config`
  (rc_canon also requires the BSC sequence compressor, like delta/hints).
- **`rc_canon` was lossy for lowercase / IUPAC bytes**: the canonicalization
  reverse-complement folded `acgt`→uppercase and non-ACGT→`N`. Replaced with a
  case-preserving **involution** (`reverse_complement_canonical`), so rc_canon
  now round-trips arbitrary sequence bytes exactly. (The shared
  `reverse_complement` used by assembly/reference is unchanged.)
- **`twobit` (2-bit + N-mask) silently uppercased** lowercase input (decode
  always emits uppercase). `encode_with_n_mask` now rejects lowercase/IUPAC,
  matching `pack_dna_2bit` — those inputs must use raw+BSC.
- **`quality_ctx` accepted `quality.len() < sequence.len()`**, which the decoder
  can't reproduce (range-coder desync). Now requires exact equality.

**Untrusted-archive robustness (decode-side DoS):**
- **Unchecked `num_reads * const_seq_len`** in the in-memory sequence fast path
  could wrap in release and slip past the length guard into a slice-panic. Now
  `checked_mul`; sibling `offset + len` checks in the in-memory codecs hardened
  with `checked_add`.
- **Cumulative decompression-bomb cap** added to the in-memory / template / mmap
  fallback decoders (headers, sequences, 2-bit, quality, paired BSC streams, the
  reference BSC roles) via a content-derived `stream_decode_cap(num_reads)`. The
  bounded-streaming path was already capped; the fallback path passed
  `usize::MAX`.
- **`header_col` varint** gained the shift-overflow guard the canonical reader
  has; **`quality_ctx`** rejects out-of-range Phred symbols (> 222) that would
  overflow `phred + 33`.

**strobealign index (`.sti`):**
- **`read_vec` no longer reinterprets a `Vec<u8>` as `Vec<T>`** (an alignment-UB
  dealloc, masked only by mimalloc) — it reads directly into a correctly-aligned
  `Vec<T>`. Length fields are bounded by file size (`checked_mul`), so a corrupt
  or foreign-endian sidecar returns `Err`/rebuilds instead of an OOM-abort.
- **Index build now rejects a single contig > `u32::MAX` bp** (and > `u32::MAX`
  references) — the strobe `position`/`ref_index` are `u32`; the previous
  `BucketIndex::MAX` guard was dead. Human chromosomes are far under the limit.

**bz:**
- **MD derivation is self-checked**: MD is stripped only when
  `regen_md(recover_ref(...)) == original`. Non-canonical MD (e.g. a matched base
  written as a mismatch) is now kept verbatim instead of producing an archive
  that fails its content-CRC at decompress (undecodable).
- **`bz extract --force`** is now honored (threaded into the downstream qz
  `CompressConfig`); it was silently dropped, so re-runs bailed on existing
  `_R1/_R2/_SE.qz`.
- **`bz extract` temp files** use pid + a process-unique nonce and `create_new`
  (O_EXCL), so concurrent runs sharing a `--working-dir` can't corrupt each other.

**concurrency:**
- On a compress **worker-thread panic**, the in-flight chunk threads are now
  joined before re-raising, so a panic no longer leaves detached BSC threads
  burning cores (matters under the Python embedding, where the panic is caught).

### qz: atomic compress output uses a randomized O_EXCL temp file

Brings the qz compress path to parity with the decompress and bz hardening: the
single-end/paired/reference encoders wrote to a deterministic `<output>.tmp`
sibling, which could collide between concurrent jobs and follow a pre-planted
symlink. The temp file is now a randomized name reserved with `create_new`
(O_EXCL) before any encoder opens it. (Centralized the `atomic_tmp_suffix` /
`create_new_for_write` helpers shared with the decompress path.)

### tests: vendor the strobealign `phix` fixtures

`crates/strobealign/tests/phix.fasta` and `phix.1.fastq` (the upstream phiX174
reference + a small read set) were referenced by the strobealign unit tests but
not tracked, so `cargo test -p strobealign` failed ~10 index/seeding tests. Added
the exact upstream fixtures; the full strobealign suite now passes.

### bz: rename the misleading `quality_ctx` identifiers to `fqz`

bz's quality codec is htscodecs **fqzcomp** (since v12), but the code still used
`quality_ctx` names inherited from the older arithmetic coder — easily confused
with qz's *actual* `quality_ctx` arithmetic coder (a different thing). The bz
identifiers were renamed for clarity (no behavior or wire-format change; the
chunk-flag bit value and multiblock layout are identical, so existing archives
still decode):

- `AdvancedOptions::quality_ctx_block_size` → `fqz_block_size` (this is also the
  JSON config key for `bz compress --config`; update old configs).
- `CHUNK_FLAG_QUALITY_CTX` → `CHUNK_FLAG_FQZ_QUALITY` (internal; same `0x01` bit).
- Internal `uses_qctx` / `use_quality_ctx` / `chunk_uses_quality_ctx` →
  `*fqz_quality`, and assorted comments/messages.

qz's `quality_ctx` (the arithmetic coder, still qz's large-input default) is
unchanged.

### Hardening: crash, DoS, and crafted-input fixes (round 2, qz + bz)

A second external review surfaced these; all fixes add regression tests.

- **bz decode/verify bound BSC decompression (DoS).** `decompress_streams` and the
  SAM-header decode used the unbounded `bsc::decompress_parallel`, so a small
  crafted archive could stack many high-expansion BSC blocks and exhaust memory
  before record-count validation. Both now use `decompress_parallel_max` with a
  per-stream cap keyed off the chunk's record count (genomic streams legitimately
  expand by large ratios, so a ratio bound would reject valid data) clamped by an
  absolute ceiling, and a fixed cap for the SAM header.
- **qz `--quality-modeling` + quality-ctx no longer panics.** quality_ctx is a
  standalone per-base coder with no model-delta backend, so the combination hit an
  `unreachable!()`. The explicit flag combo is now rejected at config validation,
  and the modeling backend returns a clean error instead of panicking (covers the
  `Auto`-resolves-to-quality-ctx case at ≥100k reads).
- **bz archive version bumped 12 → 13.** The zstd-backend removal dropped two
  header bytes without changing the version, so old v12 archives were
  version-accepted and then misparsed. v13 is rejected-or-accepted cleanly: old
  v12 archives now fail with a clear version error (re-compress them).
- **Missing `--reference` returns an error, not a panic.** `strobealign`'s
  `read_ref` used `File::open(_).unwrap()`; it now propagates the IO error (which
  also carries the underlying cause).
- **qz quality-delta decode validates its record count.** The delta stream was
  decoded until EOF, so a corrupt archive with too few/many entries silently
  produced FASTQ with blank or dropped quality lines. It now requires exactly one
  entry per record.
- **Atomic output temp files are randomized and O_EXCL (qz + bz).** The
  deterministic `with_extension`/`.tmp` sibling paths could collide between
  concurrent jobs and could follow a pre-planted symlink (truncating an unrelated
  file). bz now uses `tempfile` (random name + O_EXCL + atomic persist); qz
  decompress uses a randomized name created with `create_new`.

### Changed

- **Slimmed the vendored `strobealign` crate to the map-only core we use.**
  Reference mode only needs index build + seeding + chaining → NAMs (qz does
  ungapped-Hamming placement itself and feeds raw read bytes), so the upstream
  alignment stack (`mapper`, `aligner`, `piecewisealigner`, `ssw` + its `ext/ssw`
  C library, `cigar`), SAM output (`io::sam`), and CLI read I/O (`io::reads`,
  `io::fastq`, `io::paf`) were removed — along with the now-dead dependencies
  `block-aligner` (SIMD aligner), the `cc` build script, `clap`, `sigpipe`,
  `memchr`, and `mimalloc`. The one function the map path borrowed from the SAM
  mapper (`mapping_quality`) was only used by the removed PAF paths, and the one
  config field qz read (`MappingParameters::rescue_distance`) is now the constant
  `maponly::DEFAULT_RESCUE_DISTANCE`. Crate source **9,475 → 4,162 LOC (−56%)**;
  mapping output is **byte-identical** (subset + full chr20, exact `cmp`); 406
  tests pass. We were already maintaining a fork (the place_read / anchor-sort
  patches), so this removes ~5k LOC of dead alignment code from the build and the
  untrusted-input audit surface.

### Hardening: integrity-check and untrusted-input fixes (qz + bz)

Addresses an external review of corruption-handling and crafted-archive paths.
All fixes are lossless and add regression tests:

- **bz `verify` now validates fqzcomp quality framing.** For quality_ctx chunks
  `verify` previously CRC'd the compressed quality bytes as-is and skipped the
  inner fqz block framing, so a corrupted fqz payload (with the outer chunk CRC
  patched) could pass verification yet fail to decompress. `verify` now walks the
  multiblock framing and checks each block's per-block CRC — without running the
  expensive fqz entropy decode.
- **qz fast `verify` now validates the stored `header_size`.** The buffered/fast
  parser ignored the stored header size while the deep/streaming path trusts it as
  the data offset, so a tampered value could pass `--fast` verify but fail to
  decompress. The parser now rejects a stored `header_size` that disagrees with the
  recomputed header length, so fast and deep paths agree.
- **Truncated BAM input is no longer silently accepted.** `RawBamReader` treated a
  partial (1–3 byte) record `block_size` prefix the same as a clean EOF, silently
  shortening a truncated BAM during compress/extract. A clean EOF is now strictly
  zero bytes at a record boundary; any partial prefix errors.
- **`bz_lib::compress` validates options at the library boundary.** The public
  entry point did not call `AdvancedOptions::validate()`, so a programmatic caller
  passing `quality_ctx_block_size == 0` hit `.step_by(0)` and panicked instead of
  getting a clean error.
- **Single-end qz and bz decompress write atomically.** Both opened the target
  output directly, so a corrupt archive could leave a partial file or (with
  `--force`) destroy an existing valid result. Both now stream to a sibling `.tmp`
  and rename on success, matching the compress and paired/reference paths.
- **Reconstructed sequence/quality lengths are cross-validated.** qz now rejects a
  record whose decoded quality length differs from its sequence length (the two
  come from independent streams) rather than emitting malformed FASTQ; bz rejects a
  decoded fqz quality whose length differs from `l_seq` rather than emitting a
  structurally invalid BAM record.

### bz: fix `verify` false-positive on paired archives (derived mate fields)

`bz verify`'s structural validator rigidly expected the `next_ref_id` (stream 6)
and `next_pos` (stream 7) streams to be `num_records × 4` bytes. But the v10
PNEXT/RNEXT derivation omits both for any record whose primary mate is in the
chunk (governed by the `next_derivable` bitmap, stream 21), so on real paired
data those streams are far smaller and verify falsely reported corruption — even
though the actual roundtrip was byte-perfect. (The bug was latent because the test
BAMs used unmapped mates, which are never derivable; it surfaced on real paired
data, e.g. HG002 chr20.) The check now derives the expected size from the
next_derivable bitmap's popcount, mirroring how stream 8 (tlen) was already
handled. Added a verify test over a derivation-exercising paired BAM.

### Performance

- **Reference mode: cap `place_read`'s mismatch count at `k + 1`.** The offset
  search seeded its best-mismatch count at `usize::MAX`, so the first scanned
  offset (often a hopeless edge-of-window frameshift) was counted over all bases.
  Seeding it at `k + 1` caps every offset there: an offset with `> k` mismatches
  can neither be accepted nor be the accepted minimum, and offsets with `<= k`
  mismatches are still counted exactly, so the result is byte-identical. ~3%
  faster compress at 20 threads (HG002 chr20, won 3/3 order-swapped trials);
  in the noise at 72 threads — it's a work cut, so it shows in the work-bound
  (core-constrained) regime.
- **Reference mode: packed-key anchor sort (strobealign chaining).** The
  per-read anchor sort in `get_chains` (the largest single producer sort, ~9% of
  compress CPU) compared a 3-field `(ref_id, ref_start, query_start)` tuple. It
  now sorts by a single packed `u128` key — an order-preserving bijection of the
  same three fields — so it's one word-compare per step instead of up to three.
  `Anchor` is exactly those three fields and the sort is followed by `dedup()`,
  so the post-dedup anchor set is independent of sort details → byte-identical
  mapping output. ~2% faster compress (HG002 chr20 full, order-swapped A/B,
  consistent across trials).
- **Reference mode: allocation-free read placement (`place_read`).** Placing a
  read against the reference searched ±8 offsets and, for *every* offset,
  allocated an edit `Vec` and ran a full ACGTN-validating scan via
  `extract_substitutions` — then kept only the lowest-edit offset. A post-#2
  perf profile showed `diff_one` had become the single largest symbol (~18% of
  total CPU; mapping+diff ~62% — the encode is no longer the bottleneck). Now:
  validate ACGTN once (it's offset-independent), count mismatches per offset with
  early termination (stop once a candidate can't beat the current best — no
  allocation), and extract the edit list only for the winning offset. Same
  ascending tie-break (smallest offset among equal counts) and same `<= k` accept
  threshold → **byte-identical output**. Measured on HG002 chr20 (full 10.8M
  pairs, order-swapped A/B): compress wall **~74.1s → ~70.1s (~5%)** on top of the
  encoder change (cumulative vs the pre-batch baseline ~87s → ~70s, ~20%); larger
  on core-constrained boxes where compress is work-bound. 406 tests pass.
- **Reference mode: encode quality/headers without reconstructing sequences.**
  The per-chunk encoder reconstructed every mapped read's sequence from the
  reference (and cloned every record) only to hand it to the shared
  header+quality encoder — but reference mode encodes headers columnar (reads
  only the id) and quality via fqzcomp (reads only the quality); the sequence was
  never used, just rebuilt and discarded along with all of `records_to_streams`'
  output. A dedicated `encode_ref_headers_qual` now encodes the borrowed
  headers/qualities directly (headers ‖ quality kept concurrent via `rayon::join`,
  matching the old path's parallelism). Measured on HG002 chr20 (full 10.8M pairs,
  order-swapped A/B with 30s cooldowns): compress wall **~87s → ~74.8s (−14%)** and
  **peak RSS 7.3 GB → 5.0 GB (−31%)**. Output is **byte-identical** (same columnar
  header blob, same fqzcomp sub-blocks); 406 tests pass, roundtrip verified.
- **Reference mode: parallel insert-size sampling.** The Pass-0 insert-size
  estimate mapped the first 100K pairs in a serial loop that gates the chunk
  pipeline (a startup section that costs at every core count). It now maps the
  sample in parallel; mean/variance are order-independent so mu/sigma — and thus
  the archive — are byte-identical.

### Added

- **`--reference-fast`: sparser-seed reference mapping.** Maps with sparser syncmers
  (s=12 vs the read-length default's s=16 → ~45% fewer randstrobes), **cutting mapping
  CPU ~39%** (the deterministic underlying win). Because mapping is now on the critical
  path (after the windowed-encode and parallel-finalize work), that shows on the wall:
  **~5–10% faster compress**, more on quiet machines and on core-constrained boxes where
  compress is work-bound (measured HG002 chr20, quiet-box clean A/B: 72c 1:12→1:05, 20c
  1:30→1:24; the relative win shrinks under heavy load). **Lossless and ~ratio-neutral
  on high-coverage data** (chr20: 325.18 vs 325.22 MB — marginally smaller); the
  reconstruct-verify falls any unmappable read back to a literal, so the only effect on
  other data is possibly a few more literal fallbacks (slightly larger) on
  low-coverage/divergent inputs. The sparser index is cached separately
  (`.qz-r{len}-fast.sti`) so it never collides with the default one. Opt-in;
  byte-identical roundtrip, reference-free decode unchanged.

### Performance

- **Reference mode: parallel finalize block-stream compression.** The global
  metadata block-streams written at finalize (the FallbackPool — unmapped-read
  literals — plus IntervalMap/ReferenceMeta) were BSC-compressed one block at a
  time on a single thread, while the rest of the machine sat idle (the per-chunk
  pipeline is already done by finalize). The fallback pool was the serial long
  pole: ~10.3s of the chr20 wall. `write_record_block_stream` now partitions
  serially (cheap) then BSC-compresses the blocks concurrently via rayon.
  Measured on HG002 chr20 (window=4): fallback-pool compression **10.3s → 1.85s**,
  total compress **1:14 → 1:06 (~11% faster)**. Output is byte-identical (same
  block boundaries, same per-block compression, same order); roundtrip verified.

### Added

- **`--reference-window N` (reference mode): windowed parallel encode.** The
  reference-compress consumer (fqzcomp quality + BSC headers/positions/edits) was
  the measured wall bottleneck — it ran one chunk at a time on ~1–2 cores while
  read-mapping (the producer) sat idle. The encode now runs `N` chunks
  concurrently and writes them in order; output is **byte-identical** regardless of
  `N`. Measured on HG002 chr20 (10.8M pairs): **window=2 (default) 1.53× faster**
  (2:21 → 1:32, peak RSS 4.5 → 6.4 GB); window=3 1.77× (9.2 GB); window=4 1.85×
  (1:16, 10.6 GB). Higher windows trade peak memory (~`N` × chunk working set) for
  speed and approach the read-mapping floor; on full-WGS references the resident
  index dominates RAM, so the window's cost is proportionally smaller there.
  Default `N=2` (balanced); accepts 1–8. Roundtrip byte-identical; deterministic
  across thread counts and window size.

### Changed

- **Reference-mode chunk size lowered to 500K pairs** (from the shared 2.5M
  `CompressConfig` default the paired path still uses). Measured on HG002 chr20
  (10.8M pairs, fqzcomp qualities): peak compress RSS **14.7 → 4.4 GB** (−70%),
  compress **~7% faster** (finer producer/consumer pipeline overlap), for a
  **+0.52%** archive-size cost; byte-identical roundtrip. The chunk only sizes the
  ~17% BSC-coded side streams (headers/positions/edits/fallback) — qualities
  (fqzcomp, sub-blocked at `quality_ctx_block_size`) and the resident reference
  index are unaffected — so the ratio cost is small while the in-flight working
  set shrinks ~3.4×. `QZ_REF_CHUNK` overrides it for experimentation.

- Reference mode now compresses quality scores with **fqzcomp** (sub-blocked at
  `quality_ctx_block_size`, each sub-block cap-bounded under the 64 MiB block limit and
  decoded in parallel) instead of `quality_ctx`. Measured on HG002 chr20 (9.1M pairs):
  archive **−2.15%** smaller, compress **−2.7%** faster, decompress **−13%** faster, and
  byte-identical (lossless). Self-describing quality codec byte (`CODEC_FQZCOMP = 6`);
  existing reference archives (codec 4) still decode unchanged.

- **Reference mode rewritten as single-pass reference-direct (breaking).** `--reference`
  now encodes reads as substitution edits **directly against the supplied reference**
  in a single streaming, chunk-parallel, input-order pass, storing only the **covered
  (normalized) reference slices** in the archive. This replaces the previous
  consensus-based multi-pass design (map → window-bucket → read-derived consensus →
  place). Decompression remains fully **reference-free** and byte-identical, reads
  paired and in original order; unmappable reads and any read that fails a per-read
  reconstruct-verify fall back to literal encoding, so losslessness never depends on
  mapping or reference quality. The mode stays deterministic (byte-identical across
  thread counts).
  - **Why:** the consensus only saves the donor's germline-variant edits, which on
    real data is ~0.2% of the total archive (qualities dominate); the multi-pass
    structure forced a hard map→consensus→place barrier, a disk spool, and a 3× FASTQ
    decode. The single pass removes all three — faster, more parallel, lower-memory
    (no read pileup, no placement spool), and bounded at WGS scale (backing + footer
    streamed block-at-a-time directly to the output; never materialized).
  - **Format break:** new archive marker `0xF2` (`PAIRED_REFERENCE_DIRECT`) with a
    `reference_subformat` byte; per-read positions are genome `(ref_id, ref_pos)` with
    the R2 mate-delta retained; the read-derived consensus global is replaced by a
    covered-reference-slice backing; fallback literals are pooled into one
    archive-wide block-stream; `IntervalMap`/`ReferenceMeta` are block-streamed and
    `ReferenceMeta` carries the full contig list. **Old consensus-format reference
    archives (`0xF1`) are no longer decodable and are rejected with a "recompress"
    error.** (qz is unreleased; no field archives exist.) Non-reference single/paired
    compression is unchanged. This supersedes the consensus-era `--reference` entries
    below in this Unreleased section.
  - **Performance:** the chunk loop is pipelined (a producer thread runs the next
    chunk's map+diff while the current chunk encodes) and the encoder no longer
    BSC-compresses the sequence stream it discards (sequences are still used as
    quality-model context). Measured ~13.5% faster compress on HG002 chr20
    (2:52→2:29), byte-identical output, at ~+3 GB peak RSS (one extra chunk in flight).

### Added

- **Reference-mode compress speedups (auto-cached index + parallel mapping).**
  Two byte-identical compress accelerations for `--reference` mode:
  (1) An **auto-cached sidecar index** (IGV/BAI-style): the strobealign index is
  written next to the reference FASTA as `<reference>.qz-r<read_len>.sti` on first
  use and transparently reloaded on subsequent runs (keyed by read length, since
  the seeding parameters depend on it). The sidecar is treated as stale and
  rebuilt if its mtime predates the reference's; a corrupt/parameter-mismatched
  sidecar falls through to a fresh build. Writing the cache is best-effort — a
  read-only directory logs a warning and compress continues. A loaded index is
  bit-identical to a freshly-built one, so archives are byte-identical whether the
  cache hit or missed. No new flag (the reserved `--reference-index` stays
  rejected; this is automatic via `--reference`).
  (2) **Order-preserving parallel mapping + placement**: the per-pair seed mapping
  (Pass 1) and per-read Hamming placement (Pass 3) now run across rayon threads in
  batches while preserving original read order (I/O, spooling, consensus bucketing
  and chunk encoding stay serial/in-order). Output is byte-identical to the serial
  path and invariant across thread counts.

- **Reference-mode decode hardening + determinism tests.** Added an enumerated
  hostile-decode test suite for reference archives (truncation, bad header_size,
  out-of-bounds/overflowing footer offsets, inflated entry_count/num_chunks,
  reserved `order`, illegal codecs, unknown roles, dropped/duplicated globals and
  per-chunk roles, truncated/flipped payload streams, inflated consensus/N-bitmap
  `num_blocks`, edit-pos≥read_len, R2-delta-out-of-range): every hostile input
  now returns a clean `Err` (no panic/abort) and leaves no partial
  `_R1`/`_R2` output. Plus a cross-thread-count determinism test (reference
  archive is byte-identical for 1 vs 4 threads). Two guards added as a result:
  (1) the front-header `header_size` field is now validated against `FRONT_LEN`
  on decode/open (previously unchecked); (2) the consensus/N-bitmap block-stream
  readers and the fast-verify CRC walker now bound the declared `num_blocks`
  against remaining bytes up front (rejecting a hostile count before the per-block
  loop instead of iterating).


- **Reference archive verification (`qz verify`).** Reference-mode archives now
  verify through a dedicated path. Deep verify (`--fast` off) fully reconstructs
  both mates (enforcing every content invariant) and computes a CRC32 over the
  FASTQ bytes; fast verify CRC-checks every v4 block without bsc-decoding,
  dispatching by block framing (Consensus/NBitmap custom `[num_blocks][v4...]`
  vs the shared block-stream framing) and cross-checking declared record counts.
  `compression::verify` now classifies single/paired/reference and routes
  accordingly. Reference archives are also explicitly rejected by the
  single-stream `decompress_to_writer` and by the Python decode binding (they
  reconstruct two FASTQ files — use the CLI with `-o <prefix>`).

- **Reference-based paired-end compression (`--reference <FASTA>`).** Opt-in mode
  that maps reads to a supplied reference at compress time (vendored strobealign
  seeding — randstrobes/NAMs, no Smith-Waterman, patched for determinism) and
  stores a **read-derived consensus** plus per-read **columnar substitution edits**.
  Decompression is fully **reference-free** and byte-identical, with reads paired
  and in original order; unmappable reads and any read that fails a per-read
  reconstruct-verify fall back to literal encoding, so losslessness never depends
  on mapping or consensus quality. Compress is a bounded-memory multi-pass pipeline
  (profile + insert-size → map/spool/window-bucket → streaming core-plus-halo
  consensus → emit per-chunk streams); the archive carries no reference bases.
  The mode is deterministic (byte-identical across thread counts). Real-data
  validation (HG002 chr20, 2M GIAB pairs, ~248 bp, index built from the chr20
  FASTA): **313 MB vs the order-preserving paired baseline's 475 MB (−34%, 6.9×
  vs 4.5× on 2.1 GB input)**, byte-identical roundtrip; compress 3:51 / 9.9 GB
  peak RSS, decompress 0:33 / 5.6 GB. `--reference-index` and `--reference-fast`
  are reserved (require `--reference`, rejected in this release). The reference
  FASTA is a compress-time-only input — never stored, never needed to decode.
  Non-reference single/paired-end compression is unchanged.

- **Paired-end FASTQ compression.** `qz compress -i R1.fastq[.gz] -i R2.fastq[.gz] -o out.qz`
  stores two position-paired files in one lossless, order-preserving archive
  (chunked, bounded memory, parallel ingestion); `qz decompress -i out.qz -o prefix`
  writes `prefix_R1.fastq` + `prefix_R2.fastq`. Cross-mate header dedup. Single-file
  compression unchanged.
  - `-i` is repeatable in `compress` (once = single-end, twice = paired). The
    public `compress` wrapper validates input count (0/≥3 → error, 1 → single-end,
    2 → paired) and dispatches paired inputs to the `compression::paired` codec.
  - Compress ingests both mates in parallel (one reader thread per mate over
    bounded channels for back-pressured, memory-bounded streaming), zips them into
    ordered chunk pairs, and compresses each chunk's two mates by reusing the
    single-end per-chunk compressor (`compress_one_chunk`) via `rayon::join`.
    Per-mate quality codec is resolved independently from each mate's first-chunk
    read count (mirrors single-end Auto). R2 headers are stored as the smaller of
    an independent columnar encoding vs a per-record delta against R1.
  - Archives carry an 18-byte front header (`QZ` magic, version, paired marker
    `0xF0`, paired flag, footer offset) and a validated footer directory of
    per-(chunk, role) stream entries. Decompress parses + validates the footer,
    then decodes chunks **in order, one at a time** (memory bounded to one chunk):
    per chunk R1Headers → R1Seq → R1Qual → R2 header → R2Seq → R2Qual, so R1 ids
    feed an R2 header delta and sequences feed QualityCtx quality.
  - Paired output is published atomically: decode goes to per-mate temp files
    (removed by a drop guard on any error or panic), then both are renamed into
    place with best-effort backup-and-rollback (under `--force`, pre-existing
    targets are backed up first and restored if a rename fails). No partial output
    is left on failure. Stdout (`-`) prefix/output rejected.
  - `compression::decompress` classifies the archive up front via a validated
    tri-state `paired::archive_kind`: a well-formed paired archive routes to the
    paired decoder; a malformed-paired prefix (paired flag without the marker,
    marker without the flag, truncated front, version mismatch) **errors** rather
    than silently falling into the single-end parser. Single-end and stdin
    archives are unaffected.
  - `compression::verify` routes paired archives to `paired::verify_paired`:
    **deep** verify (`--fast` off) stream-decodes both mates and folds every
    reconstructed FASTQ record into one CRC32; **fast** verify (`--fast`)
    CRC-verifies every directory entry's v4 blocks without BSC-decoding and
    cross-checks per-entry record counts. `VerifyResult.num_reads` is `2 ×
    num_pairs`.
  - `compression::decompress_to_writer` and the qz-python `decompress` binding
    **reject** paired archives (single-stream output has no paired
    representation). New public helper `compression::is_paired_archive(&Path) ->
    Result<bool>` exposes the classifier without reaching into the private
    `paired` module.

### Fixed

- Multi-member gzip inputs (bcl-convert / bgzip) are read in full
  (`GzDecoder` → `MultiGzDecoder`); previously only the first member was decoded.
- Paired Bsc-coded quality streams were split assuming 8 bits/quality, but
  `QualityBinning::None` packs at **7 bits/quality**; `split_qual_records` now
  computes `packed_len = ceil(l_seq * bits_per_qual / 8)`. Previously, paired
  quality lengths where `ceil(l*7/8) != l` (e.g. 10 bp reads) misparsed on
  decode.

### Remove the zstd compression backend (qz + bz)

zstd was a selectable stream backend in both tools but was never the default,
never exposed by either production CLI, and consistently beaten by BSC on every
genomic stream (and no pre/post-BSC zstd stage ever helped). It was dead weight
and an unnecessary decode-time attack surface, so it has been removed entirely
(pre-production, so no archive-compatibility shim).

- **qz:** dropped the `Zstd` variants of `QualityCompressor`/`SequenceCompressor`/
  `HeaderCompressor` (remaining discriminants unchanged; compressor code `0` now
  decodes to a clear "zstd backend removed" error). Deleted the `zstd_dict`
  module and the `dict_training`/`dict_size` options it served. `quality_delta`
  now runs over BSC.
- **bz:** dropped the `alignment_compressor`/`aux_compressor` options and the two
  bytes they occupied in the archive header (header is 2 bytes shorter; offset
  constants and the content-CRC corruption test updated to match). All streams
  are BSC. Removed the now-moot zstd decompression-bomb guard added earlier.
- Removed the `zstd` crate dependency from `qz-lib` and `bz-lib`. (OpenZL's
  internal zstd graph node is unaffected — it's the OpenZL C library, not the
  Rust crate.)
- Tests: removed zstd-/dict-specific tests; repointed the rest at BSC. Full suite
  green; real-data roundtrips byte-identical (HG002 chr20 551K records, ERR3239334
  500K reads).

### Code-review hardening: losslessness fix + untrusted-input robustness

A thorough code review produced one losslessness fix and a set of
malformed-archive hardening fixes. All are lossless/ratio-neutral and covered by
new tests; the full suite plus real-data roundtrips (HG002 chr20: 551K records
byte-identical; ERR3239334 500K reads default + ultra byte-identical) pass.

- **bz: SEQ with IUPAC ambiguity codes (`M/R/S/V/W/Y/H/K/D/B`) and `=` now
  round-trips byte-exact** (previously corrupted to `N`, which then failed the
  chunk content-CRC and made the *whole archive* undecodable). The original
  nibbles were always stored losslessly (XOR-delta + raw extra); the decode path
  needlessly routed the output SEQ through an A/C/G/T/N ASCII intermediate that
  folded ambiguity bases. Decode now packs the recovered nibbles directly;
  ASCII is still derived (consistently with compress) only for MD/NM regeneration.
  Byte-identical for the common A/C/G/T/N case, decode-only, no ratio/format
  change. Removed the now-dead `ascii_to_packed_seq` path.
- **bz: reject a chunk header whose `num_records` exceeds the streams** before the
  speculative `Vec::with_capacity` (the mapq stream is exactly one byte/record, so
  it's an exact upper bound). Stops a ~30-byte crafted header from forcing a
  multi-hundred-GB allocation.
- **bz: bound the per-block `record_count` in the fqzcomp quality decoder** against
  the chunk's real record count before `vec![0; record_count]`.
- **bz: bound zstd stream decode** (alignment/aux streams) at a generous
  ratio-over-input cap to defend against zstd decompression bombs (the BSC path
  was already block-capped).
- **qz: use `checked_add` for `pos + len` in the bounded-streaming record readers**
  (varint / const / seq-with-hint / qual) — a crafted varint length could wrap the
  bounds check and panic on a malformed archive; now a clean error.
- **qz/bz: error-path resource hygiene** — the qz chunked-compress and bz
  decompress pipelines now drain in-flight worker threads on early error instead
  of detaching them (matching the existing bz-compress pattern), so a mid-stream
  failure no longer leaves BSC/decoder threads burning CPU.
- **bz (minor): guard MD/MC/NM present-bitmap lengths and the TLEN-derive mate
  index** against corrupt archives (clean error instead of an index panic in a
  decode worker); make the lossy-quality content-CRC use `saturating_sub` to stay
  in lockstep with decode. **qz (minor):** `read_le_u16/u32/u64` use `checked_add`
  for the offset.

### bz: widen compress pipeline windows — ~30% faster, still bounded memory

- Profiling (perf stat + thread scaling) showed compress was **pipeline-depth
  bound**, not CPU-bound: at `-t 16` only ~9 of 16 cores were used, and even at
  `-t 72` it plateaued at ~13 cores because the level's `window=3` only let 3
  chunks run concurrently (565K context-switches/s = a starved pipeline). The
  earlier memory work freed ~3-4× RSS, so the per-level windows are widened to
  spend that headroom on parallelism (L1 4→6, L2 3→6, L3 2→3). Ratio-neutral,
  lossless (md5-identical roundtrip).
- HG002 chr20, default L2: `-t 16` **50.3→35.2 s** (RSS 7.5→10.5 GB), `-t 72`
  **37.8→30.0 s** (8.5→13.1 GB, ~17.5 cores). The default is now *faster* than
  the original pre-levels default (39.8 s) at **~2.4× less memory** (10.5 vs
  ~25 GB). `--level 1` stays at 5.6 GB for memory-constrained use.
- Decompress was found to already be near its parallelism ceiling (window =
  num_chunks, ~14 cores; barely scales past `-t 16`) — its remaining limit is the
  serial per-chunk parse/assembly phases, not addressed here. (The smaller level
  chunks already sped decompress up as a side effect: ~55→43 s.)

### bz: build the raw quality stream lazily — lower memory + CPU

- `streams.qual` (raw varint+Phred quality) is only consumed by the BSC fallback
  path (chunks with unavailable `0xFF` quality, or `quality_compressor=1`); normal
  chunks use fqzcomp and replace that slot. It was nonetheless built — a full
  quality-sized copy + per-read varint — on **every** chunk and discarded. Now
  it's built lazily in `compress_one_chunk` from the records (which carry the true
  `0xFF` bytes that `qualities_ascii` masks away) **only when the fallback fires**.
  Byte-identical output. HG002 chr20 peak RSS: L2 9.0→7.6 GB, L1 4.6→4.4 GB.
  Together with the level bounding + `sequences_ascii` removal this cuts the
  default (L2) peak from ~25–30 GB to ~7.6 GB at the same ratio.

### bz: drop dead per-chunk `sequences_ascii` copy — lower memory + CPU

- The v12 fqzcomp quality path doesn't use the read sequence, but the build still
  materialized a full ASCII copy of every read's sequence per chunk (`sequences_ascii`,
  a leftover from the old quality_ctx base context) and ran `nibbles_to_ascii` on
  each read. Removed it — **byte-identical output** (the copy had no consumer),
  saving memory and CPU. HG002 chr20 peak RSS: L1 5.2→4.6 GB, L2 9.9→9.0 GB.

### bz: compression levels (`--level`) — bounded, predictable peak memory

- New `--level 0-3` (lossless) preset over `chunk_size` + `compress_window`, the
  single source of truth being `BZ_LEVELS` (mirrors qz's `ULTRA_LEVELS`). Fixes
  the unbounded auto-window that let peak RSS balloon to ~25–30 GB (it grew the
  in-flight window to ~70% of *whatever* RAM the box had).
- **Key finding:** bz's ratio is **~flat across levels** — every stream is split
  into fixed 25 MB BSC blocks regardless of chunk size, and the per-position
  consensus saturates at a few hundred K reads. So `--level` is a **memory/speed**
  knob, not a ratio knob (unlike qz, where bigger BSC blocks help). Measured on
  HG002 chr20 (all **lossless**, md5-identical):
  - L1 (500K chunk): **5.2 GB** peak, 1.998×, 68 s
  - L2 (1.5M, default): **9.9 GB**, 2.000×, 48.6 s (fastest)
  - L3 (5M): 18.1 GB, 2.001×, 64 s
- `--level 0` (default) = auto: the **balanced** level (L2), **downshifted to L1
  under memory pressure** (when L2's estimated peak exceeds ~70% of available
  RAM). L3 is never auto-selected (more memory, no ratio gain here) but stays
  available explicitly. This caps the memory ceiling at a level's footprint
  instead of "70% of all RAM" — e.g. default peak ~10 GB instead of ~30 GB, and
  `--level 1` gives ~5 GB at the same ratio.

### bz: optional lossy variant-aware quality reduction (Crumble/Quartz-style)

- New **opt-in, lossy** mode that flattens quality scores at confident
  reference-matching positions while preserving full resolution where it matters
  for variant calling. Off by default — bz stays lossless unless asked.
- Uses bz's per-position consensus base counts (previously discarded after the
  majority vote) to mark **variant columns** — a non-consensus allele with
  `≥ min_alt_count` reads and `≥ min_alt_frac · depth` — dilated by a `± window`.
  A base keeps full quality if its column is a keep-site, it mismatches the
  consensus, or it's an inserted base; consensus-matching confident bases and
  soft-clips are flattened. The reduced (low-entropy) quality then rides the v12
  fqzcomp path.
- **Flatten schemes** (selectable): `twobin` (preserve low qualities as a
  bad-base signal, flatten the rest — most aggressive), `coarse8` (8 Illumina-style
  bins — most conservative), `single` (one fixed Q). **Levels 1–3** tune the
  variant-detection thresholds + window (higher = more flattening).
- CLI: `bz compress --reduce-quality <1-3> [--quality-scheme twobin|coarse8|single]`
  (also via the `--config` JSON `quality_reduction` field). Sets an archive
  `FLAG_QUALITY_LOSSY` header bit (informational; decode emits the stored quality).
- HG002 chr20 (vs the 1.35 GB lossless v12 archive), all with **non-quality SAM
  columns byte-identical** to the original and the internal content-CRC validating
  the reduced round-trip:
  - `--reduce-quality 2 --quality-scheme twobin`: **407.6 MB (−69.7%)**.
  - level 1 / 3 twobin: 520 MB (−61%) / 373 MB (−72%).
  - level 2 single: 294.5 MB (−78%, most aggressive); level 2 coarse8: 937 MB
    (−30%, most conservative).
- Correctness: the per-chunk fidelity hash is computed over the **reduced** record
  bytes (spliced quality) so decode's content-CRC still validates; the lossless
  path is byte-identical (hash unchanged when reduction is off). Variant-calling
  concordance validation is a planned follow-up.

### bz: fqzcomp quality with strand reorientation — ~3% smaller WGS archives (v12)

- Quality is ~80% of a bz archive, so it's the dominant lever. The quality codec
  is switched from the hand-rolled `quality_ctx` v2 to **htscodecs fqzcomp**, fed
  **raw Phred** with **reverse-strand reads reoriented to sequencing-cycle order**
  (the QUAL of `FLAG&0x10` reads is reversed before coding, and reversed back on
  decode via the FLAG stream — nothing extra is stored). The reorientation makes
  every read's quality run in cycle order, which fqzcomp's position/history model
  exploits; `quality_ctx` v2 could not (measured +0.1%, no benefit).
- A probe (`examples/qual_fqz_probe.rs`) established the mechanism: fqzcomp and
  quality_ctx are equivalent models on non-reoriented data (±0.1%); the entire win
  is the reorientation, and it requires fqzcomp's model to realize it. (Note: this
  is BAM-specific — FASTQ reads are already in cycle order, so qz's FASTQ quality
  is unaffected and its quality_ctx remains equivalent to fqzcomp there.)
- Real A/B, roundtrip **lossless** (SAM bodies md5-identical), 21 integration tests:
  - HG002 chr20 (18.2M reads, full-resolution): 1,391,233,238 → 1,345,744,057
    bytes, **−45.49 MB (−3.270%)**.
  - NovaSeqX panel D26-1809_T (140.8M reads, binned 4-level quality):
    1,875,940,626 → 1,863,482,345 bytes, **−12.46 MB (−0.664%)** (fqzcomp also
    edges quality_ctx on binned data, independent of the reorientation win).
- BSC fallback for unavailable-quality (0xFF) chunks is unchanged.
- **Breaking:** `ARCHIVE_VERSION` 11 → 12; older archives aren't decodable —
  re-compress.

### bz: derive NM:i from the recovered reference — smaller BAM archives (v11)

- `NM:i` (edit distance) = `MD-substitutions + CIGAR(inserted + deleted bases)`,
  fully derivable from the same recovered reference already used for MD:Z
  regeneration. Reads whose MD is derivable and whose stored NM matches the
  derived value set a per-record `nm_present` bit (new stream 22) and record the
  aux tag-position (stream 23) + the source integer type byte (stream 24); the
  value is regenerated on decode. NM is only stripped when MD is also derivable
  (decode reuses MD's recovered reference). The source type byte is preserved
  rather than re-derived because it's aligner-specific (bwa-mem2 writes unsigned
  `C`, some writers signed `c`).
- Feasibility probe (`examples/nm_uq_probe.rs`) measured **100.00% byte-exact**
  NM derivability on HG002 chr20 (bwa-mem2), the NovaSeqX panel, and the older
  chr20 — value matches on all three. Real A/B, roundtrip **lossless** (SAM bodies
  md5-identical), 21 integration tests + NM unit tests:
  - HG002 chr20 (18.2M reads): 1,393,683,404 → 1,391,233,238 bytes,
    **−2.45 MB (−0.176%)**.
  - NovaSeqX panel D26-1809_T (140.8M reads): 1,882,081,670 → 1,875,940,626
    bytes, **−6.14 MB (−0.326%)**.
  Net gains are more modest than the NM stream's isolated size because BSC was
  already exploiting cross-tag correlation in the aux blob, and the new
  index/type/bitmap streams add a little back.
- **UQ:i was investigated and dropped**: absent from both production sorted BAMs
  (written only by `samtools calmd`), and on the one file that has it only 49% of
  values matched `Σ qual at mismatch positions` — the rest reflect calmd's
  BAQ-adjusted qualities, which depend on an HMM over the reference and are not
  derivable from stored data.
- **Breaking:** `ARCHIVE_VERSION` 10 → 11; older archives aren't decodable —
  re-compress.

### bz: derive PNEXT + RNEXT from the mate — ~1.5% smaller BAM archives (v10)

- For a primary, paired read whose mate is in the same chunk, `PNEXT` (next_pos)
  equals the mate's POS and `RNEXT` (next_ref_id) equals the mate's ref_id
  **exactly** per the SAM spec — pure integer equality, with no tie/sign
  convention (unlike TLEN). Such reads now set a per-record `next_derivable` bit
  (new stream 21) and omit both values from streams 6/7; decode rebuilds them
  from the already-reconstructed in-chunk mate map (the same deterministic map
  used for MC/TLEN derivation). Only cross-chunk-mate reads (the ~1% minority)
  store their values.
- Feasibility probe (`examples/pnext_probe.rs`) measured **100.00% byte-exact
  derivability** of mate-in-chunk reads on both HG002 chr20 (99.9% in-chunk) and
  a NovaSeqXPlus panel (98.9% in-chunk). Real A/B, both roundtrips **lossless**
  (SAM bodies md5-identical), 21/21 tests:
  - HG002 chr20 (18.2M reads, 8 chunks): 1,414,992,180 → 1,393,683,404 bytes,
    **−21.31 MB (−1.506%)**, ~2.66 MB/chunk.
  - NovaSeqXPlus panel D26-1809_T (140.8M reads, 57 chunks): 2,000,232,126 →
    1,882,081,670 bytes, **−118.15 MB (−5.907%)**, ~2.1 MB/chunk. The per-chunk
    saving matches chr20; the higher percentage is because the panel's binned
    4-level quality compresses ~34×, so next_pos is a larger share of the archive.
- **Breaking:** `ARCHIVE_VERSION` 9 → 10; older archives aren't decodable —
  re-compress.

### bz: parallelize records_to_streams — ~13% faster compress (byte-identical)

- Pass 2 of `records_to_streams` (per-record columnar extraction + consensus-delta
  encoding) was a serial per-chunk loop — ~10.5% of compress CPU and the main
  source of the pipeline's serial troughs (off-CPU sampling showed the running
  thread count dropping to ~9 during it). It's now split into contiguous record
  ranges (8-aligned starts) built in parallel via rayon and concatenated in record
  order. Output is **byte-identical** to the serial scan: delta-encoded fields seed
  each range from `records[start-1]`, the bit-packed streams' bytes are disjoint
  across 8-aligned ranges, and the md_exceptions maps merge in record order
  (later-record-wins, as before). Measured on HG002 chr20 (clean interleaved A/B):
  **34.8s → 30.1s (~13%)**, utilisation ~19 → ~22 cores, byte-identical archive,
  20/20 tests.

### bz: raise compress BGZF input-reader workers 4 → 16 — ~5% faster compress

- The compress pipeline reads each chunk synchronously on the main thread between
  dispatches, so the hardcoded 4-worker BGZF input reader left serial troughs
  (off-CPU sampling showed the running-thread count dropping to ~9 during reads).
  Raising the reader pool to `min(cores, 16)` keeps input inflate ahead of the
  encoders. Measured on HG002 chr20 (clean interleaved A/B): **34.3s → 32.6s
  (~5%)**, byte-identical archive (reader thread count never affects output),
  20/20 tests. Mirrors the decompress-side BGZF writer bump.

### bz: raise decompress BGZF output writers 4 → 16 — further ~1.7× faster

- The BAM output BGZF deflate is ~24% of decompress CPU (profiled) and the
  pipeline drain feeds it synchronously, so the hardcoded 4-writer cap throttled
  the whole decode pipeline. Raising it to `min(cores, 16)` lets the deflate keep
  up with the windowed decoder. Measured on HG002 chr20 (after the pipelining
  below): **76s → ~45s (−41%)**, utilisation ~8 → ~14 cores, byte-identical BAM
  output, 20/20 tests. Net decompress this cycle: **126s → ~45s (−64%)**.

### bz: pipelined (bounded-streaming) decompress — ~1.6× faster

- **Decompress now decodes chunks through a windowed pipeline** instead of one at
  a time. Previously the chunk loop was fully serial (read → decode →
  write, repeat), leaving most cores idle: ~4.5 of 72 busy, 126s on HG002 chr20
  (18.3M reads). Now up to `decode_window` chunks reconstruct on worker threads
  concurrently while finished chunks drain to the BGZF output **in archive
  order**, so the BAM is byte-for-byte identical to the serial path. Measured:
  **126s → 76s (−39%)**, ~8 cores busy. Larger archives (more chunks) scale
  further toward the core count, since the window is capped by the chunk count.
- **Bounded-memory streaming preserved.** The archive is still read strictly
  sequentially (no seeks) and peak RSS is bounded by `decode_window × chunk
  size`, constant in archive size. The window is auto-selected from core count
  and available RAM (≤70% of `MemAvailable`, capped at 16), so memory-constrained
  machines pick a smaller window automatically. The per-chunk content-CRC
  fidelity check is unchanged and still runs on every chunk.

### bz: faster compress — mimalloc + auto-scaled pipeline window (ratio-neutral)

- **`compress_window` now auto-scales (default `0` = auto).** It is the dominant
  compress-speed lever: a file has `num_records / chunk_size` chunks, and the
  window caps how many compress concurrently — the previous hardcoded default of
  4 starved the rayon pool on any file with more than 4 chunks (i.e. most real
  inputs). The window is now derived at runtime from core count and available
  RAM, keeping estimated peak RSS under ~70% of `MemAvailable` and capped at 16.
  An explicit non-zero `compress_window` in the advanced-options JSON still pins
  it. **The window never changes the archive bytes** — output is identical to the
  old default; this is purely speed. Measured on HG002 chr20 (18.3M reads, 72
  cores): **57.6s → ~36s (−38%)**, utilisation ~10 → ~19 cores.
- **mimalloc global allocator in the `bz` binary** (mirrors `qz`). Recycles the
  per-BSC-block output buffers in userspace instead of munmap-ing them to the
  kernel each block: **~8% faster compress with lower peak RSS** at scale,
  byte-identical output, no ratio change. (Negligible on tiny inputs — the win
  appears once the multi-chunk pipeline is doing real allocation volume.)
- Combined effect vs the old default: **~38% faster compress, byte-identical
  archives, verified lossless roundtrip.** Also corrected the stale
  `compress_window` doc comment (old claims of "84% utilisation / ~2 GB peak"
  were wrong: real figures are ~14% / ~20 GB at window 4).

### bz: MC:Z + TLEN mate-derivation, ~1.5–7.7% smaller BAM archives

- **MC:Z and TLEN are derived from the mate, not stored (`bz` archive version 9).**
  In a coordinate-sorted BAM a read's mate is usually in the same chunk, so its
  MC:Z (mate CIGAR) and TLEN (template length) are recomputable: MC:Z is the
  textual CIGAR of the mate, and TLEN is the signed span of the two mates' mapped
  extents. On compress, derivable MC:Z tags are stripped from aux and derivable
  TLEN values are omitted from the tlen stream; on decompress they are
  regenerated byte-exact. The mate map is rebuilt deterministically on decode
  from the read-names + flags, so **nothing extra is stored to identify mates** —
  only a per-read MC-present bit (+ tag position) and a TLEN-derivable bit, plus
  a small exception list for the minority of TLENs whose sign on fully-overlapping
  mates is aligner-internal and not recoverable.

  Measured: **−1.45%** on HG002 chr20 (TLEN ~100% derivable: 2.91 MB → 694 bytes
  per chunk; no MC tags) and **−7.71%** on a NovaSeq panel BAM (MC:Z 100%
  derivable + TLEN ~83% derivable; short fragments make TLEN a large share, and
  the aux stream goes 25× → 69×). Both verified byte-identical (SAM body md5)
  after round-trip on ~18 M and ~140 M reads.

  **Breaking:** archive version 8 → 9; older bz archives are no longer decodable —
  re-compress.

### bz: MD:Z derivation from consensus, ~2.7–3.5% smaller BAM archives

- **MD:Z tags are derived, not stored (`bz` archive version 8).** The MD tag is
  the reference bases each read covers — redundant with the per-chunk consensus
  bz already builds. On compress, MD is stripped from the aux blob and only the
  rare positions where the true reference differs from the consensus are kept (a
  per-chunk exception map — ~8.6 K entries / 15 KB for a 2.5 M-read chunk). On
  decompress the reference is reconstructed as consensus + exceptions and MD is
  regenerated byte-exact (SAM-spec form) and re-inserted at its original aux tag
  position. Three new columnar streams carry the metadata (per-read MD-present
  bit, MD tag-position, and the exception map).

  Measured: the aux stream goes from ~13× to ~24×; whole-archive **−2.70%** on
  HG002 chr20 (full-resolution Illumina) and **−3.51%** on a NovaSeq panel BAM
  (which also carries MC tags). Both verified byte-identical (SAM body md5) after
  round-trip, on ~18 M and ~140 M reads respectively. Reads whose MD can't be
  derived (no MD tag, or MD/CIGAR disagree) keep their aux verbatim.

  **Breaking:** archive version 7 → 8; older bz archives are no longer decodable —
  re-compress.

### bz: strand-aware quality context (v2), ~2.6% smaller BAM archives

- **New quality_ctx v2 context for `bz` (aligned BAM compression).** The quality
  stream — ~80% of a BAM archive — is now modeled with the context
  `(position, prev_q, prev_q2, reverse-strand)` instead of the v1
  `(position, prev_q, stability-bit, prev_base, cur_base)`. The base dimension is
  dropped (it predicts nothing for quality on real data), the stability bit is
  promoted to the full previous-previous quality, and the read's reverse-strand
  bit (FLAG & 0x10) is added. Strand is per-read side information derived from the
  FLAG stream (re-derived identically on decode) — it is **not** stored in the
  quality blob, so the change costs no extra bytes.

  The win is a **q2full × strand interaction**: neither full-q2 nor strand helps
  much alone (~−0.2% each), but together they cut the quality stream ~3% and the
  whole archive **~2.6%** on full-resolution Illumina WGS (HG002 chr20: 1,526.6 MB
  → 1,487.0 MB). Mechanism: reverse reads store QUAL reverse-complemented, so the
  full (q1,q2) trajectory plus the strand bit let the model predict which way the
  sequencing-cycle quality trend runs. Binned-quality data (e.g. NovaSeq 4-level)
  is unaffected in size (quality there is already ~34×) and round-trips losslessly.
  Both datasets verified byte-identical (SAM body md5) after round-trip.

  qz's FASTQ quality path is **unchanged** (it keeps the v1 context and is
  byte-identical); v2 is bz-only.

  **Breaking:** the bz archive version is bumped 6 → 7. Older bz archives
  (version ≤ 6, v1 quality context) are no longer decodable — re-compress them.

### New: generic token-columnar header codec

- **Generic token-columnar header codec (new `0x05` header encoding).** Headers
  that are neither SRA- nor Casava-shaped (Oxford Nanopore, PacBio, BGI/MGI,
  older Illumina `#index/mate`, plain-numbered) are now compressed by structure
  instead of falling to raw+BSC. The codec tokenizes each header (delimiter split
  + stable digit/non-digit refinement) into uniform per-position columns —
  constant, delta-coded numeric, or string — each BSC-compressed. It runs as a
  verified second tier (specific SRA/Casava → generic → raw), is size-gated
  against raw+BSC (never larger), and every emitted blob passes an in-memory
  byte-exact reconstruction check. SRA/Casava archives are byte-for-byte unchanged.
  Measured on 2.5M-record fixtures, the generic tier compresses these header sets
  ~100–1400× smaller than the previous raw+BSC fallback (e.g. Nanopore 3.16 MB →
  2.3 KB), while a degenerate high-entropy ID-only input is pre-gated straight to
  raw (no wasted work).

  **Forward-compatibility note:** new archives whose headers use the generic tier
  carry inner header version `0x05`. Older qz builds reject unknown header
  versions and **cannot decode such archives** (this only affects headers that
  previously used raw fallback — common Illumina/SRA archives are unchanged).
  Minimum decoder version for `0x05` archives: this release.

### Bug fixes (correctness / lossless round-trip)

- **`--ultra` archives are decompressable again (data-loss regression).** The
  64 MiB `BSC_MAX_BLOCK_LEN` DoS guard added for decompress robustness rejected
  the large BSC blocks the ultra reorder-local path legitimately produces
  (sequence blocks up to 750 MB, header/quality up to 100 MB). The result:
  `qz compress --ultra` on any non-trivial input (a chunk whose sequence stream
  exceeds 64 MiB — roughly ≥0.5M reads at 150 bp) reported success and a good
  ratio, but `qz decompress` failed with "decompressed block size … exceeds
  maximum … (corrupt header)", leaving the data unrecoverable. The block-size
  cap is now path-aware: the default/bounded path keeps the tight 64 MiB bound
  (its blocks are ≤25 MB), while the ultra decode path allows up to 768 MiB
  (`BSC_MAX_BLOCK_LEN_ULTRA`). Verified by a new multi-chunk ultra roundtrip
  regression test and byte-identical roundtrip of all ultra levels on a 10M-read
  file. Existing ultra archives that previously failed to decompress will now
  decompress correctly with this build.
- **Columnar headers are now verified byte-exact at compress time.** The
  structured SRA/Casava encoder output is decoded and compared to the input
  before it is kept; on any mismatch the encoder falls back to the byte-exact
  raw encoding. This makes lossless header round-trip a structural guarantee
  (mirroring the read-ID template path) rather than relying on the input-shape
  guards being exhaustive — any future reconstruction edge case degrades to a
  safe fallback instead of silently corrupting headers.
- **Columnar headers (default path) now round-trip byte-exact in more cases.**
  The structured SRA/Casava encoders previously sampled the name prefix and the
  `/1`÷`/2` pair suffix from the *first* read only, and reconstructed numeric
  fields from parsed integers. This corrupted (a) interleaved files where mates
  carry different pair suffixes, (b) merged/multi-accession files, and (c) any
  header with zero-padded numeric fields (`...:0123:...` → `...:123:...`). The
  encoders now verify the whole batch is representable and fall back to a
  byte-exact raw encoding otherwise. The structured SRA/Casava encoders also
  now require the FASTQ `@` prefix: a `>`-prefixed FASTA header that happened to
  be Illumina-shaped (7 colon fields with numeric lane/tile/x/y) was previously
  misclassified and reconstructed with a literal `@`, corrupting `>id` into
  `@>id`; such headers now take the byte-exact raw fallback.
- **`--header-template` read-ID path is now byte-exact.** It verifies the
  template encoding decodes back to the exact input and falls back to verbatim
  storage otherwise (e.g. zero-padded X/Y coordinates); the decoder no longer
  does a lossy UTF-8 conversion.
- **Multi-line (wrapped) FASTA now parses correctly.** Previously only the first
  sequence line of each record was read, so wrapped FASTA was silently corrupted.
- **FASTA now round-trips byte-exact as FASTA.** A FASTA flag (bit2 of the
  archive flags byte) is recorded at compress time, so `--fasta` input
  decompresses to valid FASTA (`>id\nseq\n`) instead of the previous malformed
  `>id\nseq\n+\n\n`. Backward-compatible: archives without the bit decode as
  FASTQ exactly as before. This now applies on every compress path: the
  `--ultra` / `--local-reorder` path, which writes its own archive header, was
  dropping the FASTA flag and decompressing FASTA-with-`--ultra` back to FASTQ.
- **`--fasta` no longer fails on large inputs.** FASTA implies no quality data,
  so the Auto quality compressor no longer tries (and errors) to encode an empty
  quality alphabet at ≥100k reads.
- **2-bit / columnar sequence encoders reject non-ACGTN bases** (IUPAC codes,
  etc.) instead of silently emitting `A`, matching `pack_dna_2bit`. (The default
  raw+BSC path was already lossless and is unaffected.)
- **quality_ctx**: fixed a model-slot index that overflowed `u16` past 65 535
  contexts, silently degrading compression on rich-alphabet data (round-trip was
  unaffected).
- **Deterministic archives**: the positional quality model now breaks read-length
  ties deterministically, so identical input always produces identical output.
- **BAM (`bz`) fixes**: `extract` no longer overflows converting an out-of-range
  quality byte from corrupt BAM; the BSC decode path bounds `n_cigar_op`/`l_seq`
  to the BAM field limits; `decompress` and `verify` now check the per-chunk
  record counts against the archive header; the stream varint rejects overflow
  instead of silently truncating.

### Robustness / DoS hardening (untrusted archives)

- fqzcomp decompression bounds the self-reported uncompressed size against the
  read count before allocating, so a tiny crafted block can't force a multi-GB
  allocation.
- Columnar header decode clamps the `num_blocks` allocation hint to the input
  size.
- BAM (`bz`) chunk-stream reads are now incremental (`take`), so a corrupt
  stream-size header can't drive a huge speculative allocation.

### CLI / API

- `qz decompress --force` (and `qz.decompress(..., force=True)` in Python) is
  required to overwrite an existing output file; decompress refuses by default,
  matching `compress`.
- `--ultra` / `ultra=` now reject out-of-range levels (valid range 1–3, or 0
  for auto) instead of silently clamping. Previously the CLI and Python bindings
  wrongly rejected `0` (auto) and wrongly accepted `4` and `5`; both are now
  corrected in all callers.
- `--ultra` / `--local-reorder` now reject `--quality-modeling`,
  `--quality-delta`, and `--dict-training` up front. The ultra path's quality
  stage only supports quality_ctx/BSC, so these options were previously accepted
  and then silently ignored; they now error instead of producing an archive that
  doesn't reflect the requested settings.

### Performance

- **Columnar-header compress is much faster: byte-exact verification no longer
  BSC-decodes (SRA and Casava).** After the allocator fix, the columnar-header
  stage became the compress critical-path long pole (e.g. 7.3s of a ~9s chunk on
  5M reads). Profiling showed a large share of it was the compress-time byte-exact
  verification, which decompressed the *entire* BSC-compressed structured stream
  and reconstructed every header just to confirm a lossless round-trip. But BSC is
  deterministically lossless (and CRC32-guarded on disk), and is already trusted
  un-verified for the far larger sequence/quality streams — so decoding it to
  "verify" is redundant. Both the SRA and Casava verifiers now reconstruct headers
  from the **in-memory pre-BSC columns** (byte-identical layout) and compare to the
  input, catching the only thing that can actually be lossy — the field
  parse/reconstruct — without the BSC round-trip. Measured: SRA header stage
  7.3s → 4.4s (5M, total compress ~9.3s → ~6.8s with mimalloc); Casava total
  compress 10.9s → 6.4s (5M, ~42% faster). Output is **byte-identical** and the
  archive is unchanged (verification method only). Lossy/edge-case headers (e.g.
  zero-padded numeric fields, interleaved pair suffixes, multi-accession prefixes,
  FASTA-shaped `@`/`>` confusion) still fall back to the byte-exact raw encoding
  exactly as before. The old decode-based `verified_or_raw` is removed.
- **Compress is ~18% faster via a pooling global allocator (mimalloc).** Profiling
  (`perf`) showed ~44% of compress CPU time was spent in the kernel page-fault
  path, not compute: each BSC block allocates a fresh ~25 MB output buffer plus
  libsais workspace and frees it, and under glibc malloc those large allocations
  are mmap'd/munmap'd back to the kernel every block, so each one re-faults and
  re-zeros its pages (693K page-faults / 113K context-switches on a 1M-read run).
  The `qz` binary **and the `qz` Python extension** now use mimalloc as their
  `#[global_allocator]`, which recycles freed buffers in userspace. Measured on
  5M reads (median of 3, same load):
  compress 11.36s → 9.34s (~18%), with wall-time variance collapsing
  (context-switches 113K → 2.7K). Output is **byte-identical** and the archive
  size is unchanged (allocator change only). Decompress is compute-bound (QLFC
  inverse coding) and unaffected. Trade-off: peak RSS is ~5–8% higher (mimalloc
  retains freed pages in its pool); the pool is bounded by the working set, not
  the input size, so it stays constant across read count. The same win was
  verified from Python (`import qz; qz.compress(...)`): 5M-read compress 9.0s vs
  ~11.4s glibc, byte-identical roundtrip. In the cdylib the allocator governs
  only Rust-side allocations; CPython keeps its own allocator and PyO3 copies at
  the boundary, so the swap is safe.
- **`--ultra` is ~2–3× faster (compress and decompress), byte-identical output.**
  Ultra's few large BSC blocks starved the rayon inter-block parallelism, leaving
  cores idle while a single-threaded BWT churned each big block. Big blocks
  (≥128 MiB) now use libbsc's OpenMP MT BWT: compression runs big blocks with
  bounded concurrency (≈cores/MT-threads at once, so total OpenMP threads and
  concurrent libsais memory stay bounded), and decompression uses the MT
  inverse-BWT under the existing bounded pipeline. Default/small-block archives
  are unchanged (rayon already saturates the cores). Measured on 10M reads:
  ultra-3 compress ~154s→~52s, decompress ~55s→~20s; ultra-2 ~87s→~30s /
  ~33s→~16s; ultra-1 ~62s→~21s / ~28s→~13s. Archives are byte-identical to the
  previous build (MT only parallelizes the deterministic BWT).

### Changed

- **`--ultra` is now a bounded-streaming big-block encoding (`UltraBigBlock`).**
  Ultra archives decompress with **constant memory** bounded by block size rather
  than materializing every record (the bound is ~8–11 GB at level 3's big blocks,
  but constant in read count — the old path was O(reads) and would exhaust memory
  at higher read counts). Compression ratio is preserved. The dead within-chunk
  read-reorder machinery and the stored permutation were removed (ratio unchanged;
  ~2000 lines deleted). `--local-reorder` is now an alias for `--ultra`.
  Fast-verify (`verify --fast`) now supports ultra archives. **Breaking:** archives
  produced by older `--ultra` builds (encoding_type 8/9, reorder-local) are no
  longer decodable — re-compress them.

### Build

- The bundled C/C++ CPU target is configurable via `QZ_TARGET_CPU` (default
  `native`). Set it to a baseline such as `x86-64-v2` for portable distribution
  builds that must not `-march=native`-SIGILL on older CPUs. libgomp/libgcov are
  now located via the C compiler rather than a hardcoded GCC-version path.

### Known limitation

- The FASTQ `+` separator line is normalized to a bare `+` on output; a
  repeated-ID separator (`+READID`) does not round-trip verbatim. This matches
  common FASTQ tooling.

## 0.3.0

### New: `--fast` compression mode

`qz compress --fast` (and `qz.compress(..., fast=True)` in Python) selects
libbsc's static QLFC entropy coder for the header/sequence streams instead of
the adaptive coder, trading a small ratio cost for faster coding. The decoder
auto-detects the coder per block, so `--fast` archives decode with any build.

Measured effect is **scale-dependent**: on small/single-chunk inputs (≤ ~2.5M
reads) it is a clear win (≈ −20% compress and decompress on 1M reads, +2% size),
because the entropy coder sits on the serial critical path there. At production
scale (10M reads, multi-chunk) the win largely disappears — the coder cost is
hidden behind block/chunk parallelism, so wall-clock is bounded elsewhere
(observed: no reliable speedup, +0.34% size). Use it for small/medium files;
for large inputs it mostly just lowers CPU usage at equal wall-clock.

### Bug fixes

- **Non-UTF-8 FASTQ header IDs now round-trip byte-exact.** FASTQ IDs are
  arbitrary bytes, but the header compressors were `String`-based and mangled
  any non-UTF-8 byte to U+FFFD — including on the default (columnar) path.
  Columnar header encoding now accepts raw bytes
  (`compress_headers_columnar_bytes`): valid-UTF-8 headers still use the
  structured SRA/Casava encoders (Illumina headers are ASCII, so this is
  byte-exact), and any header that isn't valid UTF-8 uses a byte-exact raw
  fallback whose decoder no longer UTF-8-validates. The chunked, in-memory, and
  ultra columnar paths all pass raw bytes now. The non-default template paths
  follow suit: `--header-template` (BSC) falls back to raw BSC for non-UTF-8
  headers, and the template+zstd path rejects them with an actionable error
  instead of corrupting.

- **Robustness hardening (corrupt-input / validation pass):**
  - `qz decompress --gzip-level` is now range-validated (0–9) at parse time
    (CLI and Python) instead of passing an out-of-range value to the deflate
    backend.
  - Lossless quality packing (`QualityBinning::None`) now errors on a quality
    byte it can't represent (Phred > 127) instead of silently clamping — the
    "lossless" contract is never silently violated.
  - `read_varint` rejects a 10-byte varint that overflows `usize` instead of
    silently truncating the high bits.
  - `bz` decompression validates each record's `ref_id`/`next_ref_id` against
    the reference count and `pos >= -1`, so a corrupt archive can't emit a
    structurally-invalid BAM that crashes downstream tools.
  - Decompressed-block allocations are bounded (BSC against the 64 MiB per-block
    maximum; OpenZL against its 25 MiB block size) so a corrupt size header
    can't drive a multi-GB zero-fill across worker threads (OOM amplification).
    OpenZL decode also verifies the decoder filled the buffer.
  - The bounded streaming decompress/verify paths now validate stream lengths
    with checked arithmetic against the file size before seeking/indexing,
    matching the mmap path (previously they used unchecked `u64` adds).

- **Release builds now unwind on panic instead of aborting.** The workspace
  release profile changed from `panic = "abort"` to `panic = "unwind"`. This
  lets the Python bindings turn a Rust panic into a catchable `RuntimeError`
  rather than killing the interpreter, and ensures drop guards (e.g. chunked-mode
  temp-file cleanup) actually run when a worker panics. The PyO3 `compress` /
  `decompress` entry points additionally wrap their work in `catch_unwind` so a
  panic surfaces as an ordinary `RuntimeError` (catchable by `except Exception`)
  instead of PyO3's `PanicException` (a `BaseException`). Trade-off: release
  binaries carry unwinding tables (marginally larger). Panics in `extern "C"`
  callbacks still abort at the FFI boundary, which is safe.

- **bz: decompression now proves per-chunk round-trip fidelity (archive format
  v6).** Each chunk header now stores a `content_crc32` computed over the
  original BAM record bytes at compress time. `decompress` recomputes it over
  the reconstructed records and rejects the archive on mismatch, so an
  encode/decode asymmetry can no longer silently emit a BAM that differs from
  the source. (The previous per-chunk CRC covered only the *compressed* stream
  payloads, i.e. on-disk corruption — it could not detect a reconstruction
  bug.) `verify` remains a faster integrity check; its docstring now states that
  the round-trip guarantee is enforced at decompress time. **Format change:** v5
  archives are not readable by v6 binaries — re-compress.

- **quality_ctx: fixed crash on inputs with a single distinct quality value.**
  The context-adaptive range coder requires a model total >= 2, but a
  one-symbol alphabet (every base sharing one quality value) started the total
  at 1, tripping an assertion (debug) / breaking the coder (release). A phantom
  symbol is now added for single-value inputs so the total starts at 2. This
  affected `bz` BAM compression (quality_ctx is the default) and the qz
  local-reorder / independent-quality bounded-streaming paths.

- **bz: BAM sequence reconstruction no longer panics on corrupt archives, and
  validates reconstructed length.** `reconstruct_sequence_nibbles` /
  `reconstruct_record` now return `Result`: they bail (instead of indexing out
  of bounds inside `par_iter`) when a CIGAR demands more diff/extra nibbles than
  the streams hold, short-circuit `SEQ='*'` (l_seq=0) records that carry a
  CIGAR, and verify the reconstructed base count matches the record's `l_seq`
  before writing.

- **Decode no longer panics on truncated/corrupt archives.** Three reachable
  panic vectors now return clean errors instead of aborting the process (which
  matters especially through the Python bindings, where a panic under
  `panic="abort"` would kill the interpreter):
  - `n_mask::decode_with_n_mask` now returns `Result` and validates the 2-bit
    and mask buffers against the read length; the sequence decoders bail on a
    truncated N-mask stream instead of substituting an empty slice and indexing
    out of bounds.
  - The fqzcomp quality decoder validates the sub-chunk count against the read
    count (previously `num_reads - start` could underflow and panic) and checks
    the decoded sort-key/length vectors before indexing them.
  - Columnar header encoding falls back to raw encoding when a single
    instrument:run:flowcell combo exceeds 255 bytes, instead of panicking in
    `write_combo_dict`.

- **Chunked compression: fixed silent archive corruption when read length
  changes across chunks.** Constant read-length framing (which omits per-record
  length prefixes) is decided from chunk 0, and the steady-state pipeline loop
  validated that later chunks matched — but chunks dispatched inside the
  pipeline priming window (`compress_window > 1`, default 4) bypassed that
  check. A file whose first chunk had a constant length and a later
  priming-window chunk had a different length produced an archive that failed
  to decode. The const-length consistency check now runs for every chunk before
  dispatch; mismatches are reported at compress time instead of corrupting the
  archive.

### Breaking format change

The archive format is bumped to v4. v3 archives are no longer readable —
re-encode with `qz compress`. There is no migration path within v4 binaries;
this is acceptable because qz is not yet in production use.

### New: bounded streaming decompression

`qz decompress` and `qz verify --deep` now route eligible archives through a
bounded streaming decoder that walks per-block record_count metadata in the
v4 framing, holds O(block) memory regardless of archive size, and emits the
first record within ~200 ms.

The v4 framing layout is `[block_len: u32][record_count: u32][crc32: u32][payload]`
where CRC32 covers `record_count || payload` so record-count corruption is
caught at the block layer.

Producer threads (per active stream) feed bounded mpsc channels (cap=3 blocks
per stream); a writer thread emits FASTQ records one at a time. For
`quality_compressor=quality_ctx`, the writer uses a chunk-coupled mode that
buffers `chunk_reads` decoded sequences before invoking
`quality_ctx::decompress_qualities_ctx` per quality_ctx block.

Encoder enforces:
- 63 MiB per-record cap, enforced in both the chunked compress path
  (`records_to_streams`) and the in-memory compress path (`codecs::*_bsc_with`
  helpers). A FASTQ record exceeding the cap is rejected at compress time
  rather than producing an archive the bounded decoder cannot read.
- 1,000,000-record cap on quality_ctx chunks (`quality_ctx::MAX_CHUNK_READS`).

The in-memory compress path's BSC helpers (`compress_headers_bsc_with`,
`compress_sequences_raw_bsc_with`, `compress_qualities_with` BSC branch) now
emit record-aligned v4 multi-block streams via
`bsc::compress_parallel_with_breakpoints`, populating real per-block
`record_count` values. Previously these wrote `record_count = 0` placeholders
and the decoder fell back to the legacy buffered streaming path; now all
bounded-eligible archives flow through the bounded path.

Decoder:
- `StreamIndex` walks v4 framing at decode startup without decompressing
  payloads.
- `validate_indices` cross-checks per-active-stream `total_records ==
  archive_num_reads`.
- Per-codec producers: `bounded_bsc_producer`, `bounded_columnar_header_producer`,
  `bounded_openzl_producer`, `bounded_quality_ctx_producer`.
- `decompress_to_writer` public API: bounded streaming directly to a `&mut Write`,
  with `io::ErrorKind::BrokenPipe` swallowed (returns `Ok(())` for pipe-consumer-closed).

### Excluded from bounded streaming (fall back to mmap path)

- `quality_compressor=fqzcomp` — global sort_keys stream prevents in-order
  record emission.
- `sequence_compressor=zstd` archives — raw zstd format has no per-block
  metadata (round-5 already forced all-zstd; that gating stays).
- `ultra` / `local-reorder` (encoding_type 8/9) — outer per-chunk frame can't
  be bounded below chunk size.

### Removed

- Pre-v4 `BLOCK_PREFIX_SIZE = 8` and the v3 multi-block framing
- The `V4_BLOCK_HEADER_BYTES` constant (subsumed by `BLOCK_PREFIX_SIZE = 12`)

### Known follow-ups in this branch

- Wire `PER_RECORD_MAX_BYTES` 63 MiB cap into the in-memory compress path
  (currently only enforced in `build_varint_stream_with_offsets` which the
  chunked encoder calls; the legacy in-memory path uses
  `bsc::compress_parallel_with` directly and bypasses the cap).
- Delete the legacy `stream_records_to_writer` writer once the in-memory
  compress path is migrated to `compress_parallel_with_breakpoints`. Currently
  retained as a fallback when `probe_streams_record_aligned` detects v3-placeholder
  `record_count=0` archives.

## Unreleased

### Fix: `rc_canon` + `quality_ctx` silent quality corruption

When `rc_canon` is enabled, sequences are stored canonicalized (the lex-smaller
of `seq` vs `revcomp(seq)`). The encoder was passing the ORIGINAL sequences to
`compress_qualities_ctx` while the decoder reads CANONICAL sequences from the
archive and passes those to `decompress_qualities_ctx`. For reads where
canonicalization actually reverses the sequence (rc_flag=1), the encoder's
quality-context model saw a different context than the decoder's, producing
silently corrupted quality scores on roundtrip.

The encoder now canonicalizes sequences before feeding them to
`compress_qualities_ctx` so both sides agree on context. The combination is
not currently exposed via the CLI (`rc_canon` is library-only), so no archives
in the wild are affected. Regression test:
`test_rc_canon_plus_quality_ctx_roundtrip`.

Note: qualities are still stored in original-read direction while the
context is canonical, so for `rc_flag=1` reads the model learns a
reverse-complement context against original-direction qualities. This is
correct (encoder/decoder agree) but sub-optimal compression-wise; future
work could also reverse qualities for `rc_flag=1` reads before encoding
to align context and qualities, at the cost of touching the BSC quality
path for symmetry.

### Performance: parallel block decompression in bounded producers

`bounded_bsc_producer`, `bounded_columnar_header_producer`, and
`bounded_openzl_producer` now dispatch per-block decompression to rayon
instead of running it sequentially on the producer thread. The producer
still reads compressed payloads sequentially from disk (cheap I/O — 5–25 MB
per block), but each block's CRC verify + decompress runs as a rayon task
with up to `MAX_PARALLEL_DECOMPRESS = 16` in flight per stream. Results are
collected via an ordered drain so the writer continues to see records in
the original archive order.

This recovers the within-stream rayon parallelism that the old
`bsc::decompress_parallel` path exploited. On a 10M-read benchmark
decompress wall-clock went from 90.8 s (Phase 0) → 46.5 s (Phase 1, MAX=4)
→ 16.0 s (Phase 2, MAX=16 + parallel writer). Peak RSS ≈ 5.9 GB — still
bounded independent of archive size:
`(BOUNDED_CHANNEL_CAP + MAX_PARALLEL_DECOMPRESS) × max_block_size` per
stream + one transient record-format batch.

### Performance: parallel FASTQ assembly in bounded writers

Both `bounded_write_records_independent` and `bounded_write_records_quality_ctx`
now batch records (`FORMAT_BATCH_SIZE = 25_000` per batch) and assemble
FASTQ bytes in parallel via rayon `par_iter` sub-chunks before writing them
sequentially to the output. The per-record FASTQ-assembly logic (RC apply,
quality unpack / Phred+33 emit, header + sequence + separator + quality
formatting) is extracted into a pure free function `assemble_fastq_record`
with `QualityInput::{None, Packed, Decoded}` variants covering all three
quality codec shapes the bounded path supports.

For the quality_ctx coupled writer the per-block `decompress_qualities_ctx`
decodes are also now pipelined: up to `MAX_INFLIGHT_DECODES = 8` blocks
decode concurrently on rayon while the writer drains the front of the FIFO.
The decoder itself remains single-threaded (it's a global range coder), but
multiple blocks make progress in parallel.

Combined with the producer-side `MAX_PARALLEL_DECOMPRESS` bump (4 → 16)
this closes the gap to the pre-bounded-streaming master baseline on a 10M
quality_ctx archive: 46.5 s → 16.0 s (~2.9× speedup, parity with the
pre-bounded path at 15.13 s). Independent (BSC quality) archives improve
from 25.6 s → ~19 s.


### Post-review fixes (round 5)

A second external (Codex) review of the round-4 tree found three more real correctness bugs in the Zstd-sequence + verify-fast + rc_canon paths, plus a smaller block-size mismatch and one footgun in a public helper. All addressed:

- **`--sequence-compressor zstd` now requires `--quality-compressor zstd`** — the Zstd sequence path routes through `compress_columnar`, which zstd-encodes sequences, n-masks AND qualities. With quality_modeling/quality_delta off, the columnar zstd qualities were written verbatim while the archive header recorded the user's (resolved) quality_compressor. Mismatched codecs → archive unreadable AND `qz verify --fast` misreads zstd magic as a v3 `num_blocks` header. Now rejected at compress time. Regression: `test_zstd_sequence_rejects_non_zstd_quality`.
- **`qz verify --fast` no longer false-fails on Zstd-sequence archives** — the fast verifier previously hardcoded `sequences` and `nmasks` as v3/BSC-framed. For `sequence_compressor=Zstd` archives, the streams are raw zstd: the v3 walker mis-read the zstd magic as `num_blocks=4.25 billion` and bailed with a confusing BSC error on healthy files. Now `seq_is_v3` is computed from `hdr.sequence_compressor` (mirroring the existing `header_is_v3` / `qual_is_v3` dispatch), and zstd-encoded sequence streams are skipped and counted in `streams_skipped`. Regression: `test_verify_fast_handles_zstd_streams`.
- **`qz verify --fast` no longer reports OK on truncated rc_canon archives** — when `encoding_type == 6`, the verifier only walked the rc_flags stream `if q_end + 8 <= data.len()`. A file truncated to exactly `q_end` (or anywhere in the next 8 bytes) silently dropped the required stream and reported `Status: OK`. The rc_flags stream is required to reconstruct sequences in encoding_type=6 archives; truncation there is corruption, not an optional-stream skip. Now bails with a clear `"rc_canon archive truncated"` error. Regression: `test_verify_fast_rc_canon_truncation_fails`.
- **BSC block-size cap aligned with libbsc worst-case output** — the compress side caps BSC input blocks at 64 MiB, but the decompressor refused compressed blocks `> 64 * 1024 * 1024`. libbsc's worst-case output is `input + LIBBSC_HEADER_SIZE` (28 bytes), so a 64 MiB block of fully-incompressible data could exceed the decompress cap by ~28 bytes and become undecompressible. Both BSC and columnar-header decoders now share a `BSC_MAX_BLOCK_LEN = 64 MiB + 4 KiB` margin, matching the conservative `compress_with_params` output-buffer sizing.
- **`quality_ctx::wrap_as_multiblock` now writes v3 framing with CRC32** — the helper was emitting the pre-v3 layout `[num_blocks][block_len][blob]`, but `decompress_quality_ctx_multiblock` was upgraded to v3 expecting `[num_blocks][block_len][crc32][blob]`. Using the writer would have produced an unreadable archive. No internal call sites exist, so this was a latent public-API footgun rather than a regression; the helper now emits the matching v3 layout.
- **`build.rs` clippy hygiene** — six lints (`needless_borrows_for_generic_args`, `manual_flatten`, `collapsible_if`) under `-D warnings` are gone. The remaining workspace-wide clippy lints predate this branch and are out of scope.

### Post-review fixes (round 4)

A regression review of the round-3 fixes caught three real gaps in the new compress-time rejections and the new `num_reads` cap. All addressed:

- **`--dict-training` now requires _explicit_ `--quality-compressor zstd`** — the round-3 fix permitted `Auto` as well, but `Auto` resolves to BSC or QualityCtx (both non-Zstd) for typical inputs, leaving the silently-corrupt-archive bug fully reachable from `dict_training: true` + default config. Now any value other than explicit `Zstd` is rejected. Regression test added.
- **`--quality-delta` now requires explicit `--quality-compressor zstd`** — same shape of bug: `compress_qualities_with_delta` unconditionally emits raw zstd, but the archive header stored the user-typed compressor. With a non-Zstd choice, the decompressor read the wrong codec from the header and `qz verify --fast` misread the zstd magic as a v3 `num_blocks` header. Rejected at compress time with a clear error.
- **`num_reads` cap no longer triggers a 176 GB scan-offset allocation** — round-3 loosened the cap to ~11 billion to accommodate highly-compressible inputs, but `scan_header_offsets` / `scan_seq_offsets` / `scan_qual_offsets` still pre-allocated `Vec::with_capacity(count)` of 16-byte tuples (176 GB at the cap). Each scan now caps its `with_capacity` hint at `count.min(data.len())` via a new `scan_capacity` helper — tight against the actual decompressed stream and unaffected for legitimate inputs. The same `scan_capacity` clamp is now applied to the rest of the reachable post-decompression record-build sites: `decode_header_stream` / `decompress_sequences_raw_bsc` (and its delta cache) / `decompress_sequences_raw_openzl` / `decompress_sequences_2bit_bsc` in `codecs.rs`, `decode_read_ids` in `read_id.rs`, `decompress_raw_fallback` in `header_col.rs`, the zstd-sequences in-memory path in `decompress_impl.rs`, and both ANS arithmetic decoders. Other `Vec::with_capacity(num_reads)` sites (ultra/local-reorder paths, in-codec scratch caches) remain bounded by the parse-time `num_reads` cap; that cap's intent ("≤256 GB of per-record `Vec<u8>` metadata") is now documented on `MAX_RECORD_PREALLOC_BYTES`.
- **`qz verify --fast` Status line distinguishes OK from PARTIAL** — was always "OK" even when streams were skipped. Now prints `Status: PARTIAL (N stream(s) skipped)` when `streams_skipped > 0`, so a downstream `grep "OK"` doesn't conflate partial coverage with full integrity.
- **README ASCII diagram corrected** — the round-3 prose fix mentioned "bounded channels" was wrong, but the diagram still showed `bounded channels (capacity = 2)`. Diagram now describes the actual flow: `verify CRC + decode into per-stream Vec (drain_channel)`.

### Post-review fixes (round 3)

Three external reviewers + Codex independently flagged real bugs in the v3 verify path and in two compress-time config paths. All addressed:

- **`qz verify --fast` no longer false-passes on RC-canon archives** — previously slurped only headers/sequences/nmasks/qualities, leaving the RC flags stream (encoding_type 6) entirely unchecked. Bit-rot in RC flags would report "Status: OK". Now walks the RC flags stream via the same v3 block-CRC verifier.
- **`qz verify --fast` no longer false-fails on raw-zstd streams** — previously read the zstd magic as `num_blocks=0xFD2FB528` and bailed with a confusing "BSC parallel: num_blocks=4.25 billion exceeds remaining payload" on healthy archives. Now detects the per-stream compressor and skips non-v3 streams (counted in the new `VerifyResult::streams_skipped`); the CLI prints how many streams were skipped and tells the user to rerun without `--fast` for full coverage.
- **`qz verify --fast` refuses ultra/local-reorder archives** — the outer per-chunk frame in encoding_type 8/9 has no per-block CRC, so a fast walker would mis-interpret chunk_len as block_len. Now bails up front with `"fast verify does not support encoding_type=8/9 (ultra/local-reorder)"`.
- **`--dict-training` now requires `--quality-compressor zstd`** — previously fell back to BSC/OpenZL/fqzcomp silently while the archive header still claimed `dict_present=1`, producing archives that failed to decompress with "Unknown frame descriptor". Now errors at compress time with a clear "use zstd or disable dict-training" message.
- **`--quality-modeling` + `--quality-compressor fqzcomp` now rejected at compress time** — fqzcomp doesn't support quality modeling, so the encoder silently substituted BSC for the model deltas while the header still advertised Fqzcomp. The resulting archive was unreadable. Now errors with a "pick another compressor or disable quality modeling" message.
- **rayon thread-pool conflicts now logged** — `qz_lib::compression::{compress,decompress,verify}` calls `ThreadPoolBuilder::build_global()`, which can only succeed once per process. Subsequent calls with a different thread count were previously silently ignored (a footgun for `qz-python` users calling `compress()` repeatedly). Now logs a `warn!` with both the requested and the actual thread count. The bz-lib design (init in the binary, not the library) is the cleaner long-term fix.
- **Streaming-decompress comments corrected** — module/function docs claimed "bounded channels" and "~300 MB peak memory" but the channels are `mpsc::channel` (unbounded) and `drain_channel` accumulates each stream's decompressed bytes in a contiguous Vec. Comments now match reality and point at `decompress_to_records` for the truly memory-bounded path.
- **README updated to v3 format** — magic-byte/version examples, decompression-pipeline diagram, and "Bounded decompression channels" line all corrected. New "Multi-block stream layout (v3)" section documents the per-block CRC32 framing shared by all five codecs.
- **Wide-spatial header mixed-batch test added** — explicit test that a Casava batch where one row overflows u16 promotes the whole archive to v0x03 (per-archive width, not per-row — guards against a future per-row refactor breaking the format silently).
- **`block_crc32` factored into a `pub(crate)` helper in `bsc.rs`** — 13 inline `flate2::Crc::new()/update/sum()` copies across 8 files collapsed to a single call site. Eliminates drift risk if the algorithm ever changes.
- **`quality_context` module clearly documented as bench-only** — non-v3 framing; module docstring spells out that its output must never be written into a qz archive.
- **`parse_archive_header` `num_reads` cap tightened to drive downstream allocation** — was `num_reads <= data.len()` (broke roundtrips with very compressible data where the read count legitimately exceeds the archive size); now caps at `MAX_RECORD_PREALLOC_BYTES (256 GB) / sizeof(Vec<u8>) ≈ 11 billion`, which bounds the per-record Vec allocations while still allowing extreme compression ratios.

### Workspace builds and `experimental` feature

The `experimental` feature gates `pub mod template`, the bench-only OpenZL graph variants, and the bench-only `quality_context` module. Under `cargo build --workspace` (or `cargo test --workspace`), Cargo's unified feature resolution **does** enable `experimental` for `qz-cli` and `qz-python` because `qz-bench` declares `qz-lib = { features = ["experimental"] }`. The gate fully isolates experimental code only when building production crates individually: `cargo build -p qz-cli --release` or `cargo build -p qz-python --release`. Production release builds and PyPI publishing should use the per-crate form.

### Archive Format

- **Per-block CRC32 (qz archive v3)** — Every block in every codec stream (BSC, OpenZL, columnar headers, fqzcomp, quality_ctx) now carries a CRC32 (IEEE) over its compressed payload. The CRC is computed at write time and verified before invoking the inner codec, so bit-rot and disk corruption are caught with a clear "block N CRC32 mismatch" error instead of bubbling through libbsc as a cryptic `block_info failed` or producing silently wrong output. Wire layout per block: `[block_len: u32 LE][crc32: u32 LE][payload]`. Matches the pattern bz-lib has used since v4 (commit e2a877b); shares `flate2::Crc` so the two crates compute identical checksums.
- **Breaking change**: v2 archives are not readable by v3 builds. The version-mismatch error includes a clear "please re-encode with `qz compress`, or use an older qz binary to decompress" message.
- **`qz verify --fast`** — New fast verification mode walks the per-block CRC32s without invoking any codec decoder. Catches bit-rot in O(IO + CRC) time instead of O(BWT). Measured locally on a 500K-read Illumina archive (22 MB compressed): fast verify completes in ~10 ms vs deep verify's 3.8 s — roughly **two to three orders of magnitude faster** depending on archive size and number of streams. Reproduce: `qz compress -i sample.fastq -o test.qz && time qz verify -i test.qz --fast` vs `time qz verify -i test.qz`. Fast verify catches the same single-byte flip the deep mode catches; deep verify (with full FASTQ-bytes CRC32) remains the default.
- **`VerifyResult` gained `mode: VerifyMode` and `blocks_verified: u32`** to distinguish deep vs fast results. Now `Debug + Clone`.

### Hardening

- **`parse_archive_header` caps declared stream sizes against the archive size**. A hostile archive declaring `headers_len = u64::MAX` (or any of the other four stream lengths) previously wrapped the body-end check and bypassed truncation detection, then triggered huge Vec allocations downstream. All four stream lengths and `num_reads` are now validated `<= data.len()` up front, and the body-end calculation uses `checked_add` so wraparound is caught.
- **`bsc::decompress_parallel` / `openzl::decompress_parallel` bound `num_blocks * 9` against remaining payload** (was `* 5` pre-CRC). Catches malformed `num_blocks` before `Vec::with_capacity` allocation.
- **`quality_ctx` decoder caps `read_len × num_reads` pre-allocation at 256 GB** with a clear error. Previously a corrupt header could ask the decoder to pre-allocate petabytes and OOM the host before any decoding failure could surface.

### Cross-pollination from bz-lib

- The per-block CRC32 design directly mirrors bz-lib's `ChunkHeader.crc32`. The shared decoder helpers (`bsc::verify_parallel_crcs`, `flate2::Crc`) mean qz and bz now produce checksums comparable byte-for-byte.

### Header Compression

- **Casava v0x03 / SRA v0x04 wide-spatial formats** — The columnar header codec now writes 32-bit `tile`/`x`/`y` fields when any value exceeds `u16::MAX` (NovaSeq X x/y coordinates reach ~70k). Typical Illumina datasets continue to use the existing v0x02 (Casava) / v0x01 (SRA) layouts — byte-identical to pre-widening output — so older archives stay compatible and current archives don't grow. Decoder accepts both widths.

### Experimental Code

- **`pub mod template` gated behind the `experimental` feature** — The entire template-based hybrid sequence encoder and its supporting machinery (TemplateGraph, mapping, etc.) was only ever called by `compress_sequences_template_hybrid` (which was already gated) and `qz-bench`. Production qz-cli / qz-python builds no longer ship the code, eliminating 18 "unused symbol" warnings.
- **`openzl` bench-only graph variants gated behind `experimental`** — `compress_ace`, `compress_clustering`, `compress_delta_entropy`, `compress_dna_numeric`, `compress_field_lz`, `compress_fse`, `compress_huffman`, `compress_transpose_zstd` and their `*_graph_fn` helpers were only called by `bench_openzl_graphs.rs`. They no longer reach production binaries.

### Tests

- **`assert_cmd` smoke tests for `qz` and `bz` binaries** — 8 + 5 = 13 new tests covering: `--help` and `--version`, friendly error on missing input, overwrite refusal without `--force`, end-to-end compress→decompress→cmp roundtrip via the binary, `verify` happy path, and CLI flag-conflict rejection (`--fasta`/`--no-quality` × `--quality-mode`).

### Security / Safety

- **Archive parser hardened against malformed input** — `parse_archive_header` previously indexed raw bytes (`data[offset]`) after offsets were advanced by attacker-controlled `num_lengths * 4`; a crafted archive could trigger an index-out-of-bounds panic (DoS). All single-byte reads now go through a bounds-checked `read_u8` helper; the read-lengths block is validated against the archive size before being skipped.
- **`u8_to_u32` no longer panics on truncated payload** — Both `arithmetic_sequence` and `arithmetic_quality` decoders previously asserted that the compressed payload length was a multiple of 4, panicking the process on corrupt input. They now return a clear `anyhow::Error`.
- **`bsc::decompress_parallel` rejects pathological `num_blocks`** — A hostile archive declaring `num_blocks = 0xFFFFFFFF` would previously try to `Vec::with_capacity(~32 GB)` of slice slots. The block count is now validated against the remaining payload before allocation.
- **`pack_dna_2bit` returns an error on IUPAC bases** — Ambiguity codes (R, Y, S, W, K, M, B, D, H, V) were silently mapped to A by the 2-bit codec and to N by the 4-bit codec, corrupting reads round-trip. The codec now returns a clear error; callers needing IUPAC support must route via the raw+BSC sequence path.
- **`debruijn` decoder bails on out-of-bounds mismatch position** — Previously a mismatch position past the read length was silently dropped while subsequent deltas continued from the bad position, irreversibly corrupting later mismatches. Now errors explicitly.
- **`quality_model` and `quality_delta` reject Phred > 93** — The decoders clamp reconstructions to `0..=93`, so encoding higher values silently corrupted PacBio/ONT data. Encoders now refuse such inputs with a clear message pointing at the offending byte.
- **Lossy quality modes emit a `warn!` log at compression start** — `Discard`, `IlluminaBin`, `Illumina4`, `Binary`, and `--no-quality` previously ran silently; the user-facing log now states clearly that the archived qualities are not recoverable.
- **Mismatch positions widened from `u16` to `usize` in memory** — `debruijn` and `greedy_contig` encoders held mismatch positions as `u16`, which silently truncated for reads > 65535 bp. The wire format is unchanged (varint), so existing archives still decompress.
- **Output atomic-write + `--force` flag** — `qz compress` (and now `bz compress` / `bz decompress` / `bz extract`) refuses to overwrite an existing output unless `--force` is set. Internally writes to a sibling `.tmp` file and atomically renames on success; a crash mid-compress no longer leaves a half-written archive at the user's target path.
- **CLI flag combinations now validated** — `--fasta` + `--quality-mode <non-lossless>` and `--no-quality` + `--quality-mode <non-lossless>` were silently accepted (with surprising semantics). They now error at parse time.

### libbsc Threading

- **`LIBBSC_FEATURE_MULTITHREADING` no longer enabled at `bsc_init`** — Initializing libbsc with the multithreading flag pinned process-wide OpenMP state and contradicted CLAUDE.md's "rayon owns parallelism" rule. Per-call enablement on the two MT entry points works without it. `ensure_initialized` now propagates init failure as a `Result` instead of warning and continuing into UB territory.
- **`compress_with_params` / `compress_mt_inner` now `shrink_to_fit()` the output buffer** — Without this, the worst-case `input_size + 1052` capacity per block could pile up to ~35 GB of unused headroom on 1400-block streams (documented in CLAUDE.md but not previously honored at the source).

### Python bindings (qz)

- **GIL released around long-running `compress`/`decompress` calls** — Previously a multi-minute compression call blocked the entire Python interpreter (no other thread, no signal handling, no Ctrl-C). Now wrapped in `py.allow_threads(|| ...)`.
- **`tracing` subscriber installed on first call** — `info!`/`warn!` messages from qz-lib (including the new lossy-quality warning) used to be silently dropped; they now reach stderr. Honors `RUST_LOG`.
- **`pyproject.toml` enriched** — Added classifiers, readme, license, project URLs, keywords, and explicit `module-name`. Version bumped to 0.2.0 to match the workspace.
- **`compress(..., force=False)` parameter** — Mirrors the CLI flag; defaults to refusing to overwrite an existing output.

### FASTQ Reader

- **UTF-8 BOM stripped from the first record** — Some editors and Windows tools prepend EF BB BF; without stripping it landed in the first ID and the next record's `+` check failed with a misleading "expected '+' separator" error. Now silently consumed.
- **Edge-case tests added** — CRLF line endings, BOM, and truncated record handling are now covered by unit tests.
- **Error message no longer debug-prints paths** — `Failed to open file: "/path/x"` (with quotes) is now `Failed to open file: /path/x`.

### Experimental code

- **`compress_sequences_template_hybrid` gated behind `experimental` feature** — This encoder writes a "THB1" archive but has no corresponding decompressor. The function is now `#[cfg(feature = "experimental")]` so production builds cannot accidentally produce unreadable archives. `qz-bench` enables the feature explicitly.
- **`quality_context` module also gated behind `experimental`** — Bench-only alternative quality coder, distinct from the production `quality_ctx`. Now feature-gated so a future contributor cannot pick the wrong file by accident.

### bz-cli

- **`bz verify` shows unknown compressor codes** — Previously rendered any non-zero byte as `zstd`, hiding archive corruption.
- **All three subcommands accept `-f`/`--force`** — Compress, Decompress, and Extract refuse to overwrite an existing output unless the flag is set. Extract additionally checks all three derived outputs (`{prefix}_R1.qz`, `_R2.qz`, `_SE.qz`).

### bz-lib

- **`extract` pair-buffer hard-capped at 50M unpaired reads** — Previously only warned at 10M and grew unboundedly, OOM-ing the host on name-mismatched BAMs. Now errors with a clear pointer at `samtools sort -n`.

### Build

- **`build.rs` checks for `third_party/libbsc` and `third_party/htscodecs` submodules** — Clearly tells the user to run `git submodule update --init --recursive` instead of failing with an opaque `cc` error mid-compile.

### Internal

- **`fast_ultra` field removed** — Was marked deprecated but still wired; now gone. Use `--ultra 2` instead.
- **Unified thread-default rendering in qz-cli** — `--threads` now defaults to `num_cpus()` in all three subcommands (was inconsistent between `"0"` and `num_cpus()`).
- **`compress_window` documentation matches the default** — Comment claimed default 2; actual default is 4.
- **`fqzcomp::decompress` short-circuits on `num_reads == 0`** — Avoids passing a zero-length-Vec ptr into htscodecs.
- **`encode_reads` (debruijn) sorts unitig IDs before assigning compact IDs** — Was iterating `FxHashSet` directly, producing nondeterministic archive bytes for identical inputs.

### Breaking Changes

- **Archive format v2**: Archives now start with an 8-byte prefix (`QZ\x02\x00` + header_size u32 LE). Archives produced by previous versions are no longer readable.
- **FastqRecord uses `Vec<u8>`**: `id`, `sequence`, and `quality` fields are now byte vectors instead of `String`.
- **`CompressConfig.force: bool` added** — Default `false`; existing constructors using `..CompressConfig::default()` are unaffected. Library callers must opt in if they want to overwrite.
- **`CompressConfig.advanced.fast_ultra` removed** — Use `ultra: Some(2)` instead.
- **`quality_model::encode_with_model` and `quality_delta::encode_quality_deltas` return `Result`** — Bail on Phred > 93 rather than silently corrupting.
- **`dna_utils::pack_dna_2bit` returns `Result`** — Bails on IUPAC bases rather than silently corrupting.
- **`debruijn::encode_reads` and `greedy_contig::encode_reads` return `Result`** — Propagate from the changed `pack_dna_2bit`.

### Added

- **`qz verify` — archive integrity verification** — Fully decompresses all streams without writing to disk, verifying that the archive is intact and decompressible. Reports archive metadata (num reads, encoding type, compressor info, stream sizes) and a CRC32 hash of the reconstructed FASTQ output. The CRC32 and byte count are consistent across compression modes (default, ultra) for the same input. Usage: `qz verify -i archive.qz`.
- **Columnar header compression (new default)** — Illumina headers are parsed into typed binary columns (read_num, instrument combo, lane, tile, x, y) and each column is BSC-compressed independently in parallel. Falls back to raw BSC for non-Illumina formats. Saves ~8.5% on headers (e.g., 84.1 MB vs 91.9 MB on 10M reads), improving overall ratio from 7.99x to 8.13x. Selectable via `--config '{"header_compressor": 3}'` (Columnar) or `--config '{"header_compressor": 1}'` (BSC).
- **`bz verify` — archive integrity verification** — Fully decompresses all streams without writing to disk, verifying that the BZ archive is intact and decompressible. Reports archive metadata (num records, num chunks, compressor types) and a CRC32 hash of the reconstructed BAM record data. The CRC32 is consistent across runs for the same archive. Usage: `bz verify -i archive.bz`. Exits with code 1 on failure.
- **`bz extract` — BAM to QZ extraction** — New `bz extract` subcommand reads a coordinate-sorted BAM file, pairs reads by name (R1/R2), and compresses to paired QZ archives (`{prefix}_R1.qz` and `{prefix}_R2.qz`). Reads are in matched pair order. Secondary and supplementary alignments are skipped. Coordinate-sorted read order provides good BSC compression locality.
- **`--config` JSON flag** — Pass a JSON file to override `AdvancedOptions` at runtime (`qz compress --config opts.json ...`). Only fields present in the JSON are overridden; missing fields use defaults. Enables automated parameter tuning without recompilation.
- **Compression parameter optimizer** — Python-based GA optimizer (`scripts/optimizer/`) searches QZ's 17 tunable parameters (compressor selection, encoding modes, block sizes) to maximize compression ratio. Includes sensitivity analysis, checkpoint/resume, LLM-guided mutations via Ollama, and parallel evaluation. Run with `python3 -m scripts.optimizer.main --qz-bin target/release/qz --input data.fastq`.
- **BZ compression parameter optimizer** — Separate GA optimizer for BAM compression (`scripts/bz_optimizer/`) with 8 tunable parameters across 4 groups: chunking (chunk_size, bsc_block_size_mb), quality (quality_compressor, quality_ctx_block_size), encoding (use_lzp, bsc_adaptive), and per-stream compressor routing (alignment_compressor, aux_compressor). BZ now supports `--config` JSON for custom compression parameters. Run with `python3 -m scripts.bz_optimizer.main --bz-bin target/release/bz --input data.bam`.
- **CUDA GPU acceleration (opt-in)** — Build with `--features cuda` to enable GPU-accelerated BWT via libcubwt. Falls back to CPU gracefully when GPU is unavailable. Warns at startup if GPU VRAM is too small for the configured block size.
- **stdin/stdout piping** — Use `-` for `-i` or `-o` to read FASTQ from stdin or write archives/FASTQ to stdout. Supports full pipe chains (`cat reads.fq | qz compress -i - -o - | qz decompress -i - -o -`). Decompression from stdin spools to a temp file (decompressor needs seeking).
- **Parallel gzip output** — Decompression with `--gzipped` now uses multi-threaded gzip via `gzp` for faster output.
- **Tracing to stderr** — All log/tracing output goes to stderr, keeping stdout clean for piped data.
- **Archive format v2** — Magic bytes (`QZ`), version field, and self-describing `header_size` enable file identification and forward-compatible header evolution.
- **Memory-mapped decompression** — The in-memory decompression path uses `memmap2` instead of `read_to_end`, eliminating heap allocation for the full archive.
- **Constant-length read optimization** — When all reads have the same sequence/quality length, per-read varint framing is skipped, saving ~1 byte per read.
- **Codec/Decoder traits** — Unified `StreamCodec` and `StreamDecoder` traits for stream compression and decompression, replacing ad-hoc dispatch.
- **Unified chunked compression** — Single `compress_chunked` orchestrator handles BSC, fqzcomp, and quality_ctx quality paths, replacing duplicated pipelines.

### Removed

- **Dead experimental encodings** — Encoding types 1 (delta), 2 (RLE), 3 (de Bruijn), 5 (paired-end), 7 (factorized), and 8 (local-reorder-delta) removed. Only types 0, 4, 6, 8, 9 remain.
- **Experimental bench-only config fields** — `CompressConfig` no longer carries fields that were only used by benchmark binaries (e.g., `sequence_delta`, `sort_chunks`, `local_reorder`). These moved to benchmark-local config.

### Performance

- **2.5x faster decompression (10M reads: 33s → 13s on 72 cores)** — Three serial bottlenecks in the in-memory decompression path were parallelized: (1) columnar header block decomposition now processes all 4 chunks via rayon::par_iter instead of sequentially; (2) the 6 BSC column arrays within each header chunk now decompress in parallel; (3) constant-length sequence parsing uses rayon::par_iter instead of a serial memcpy loop. The dominant bottleneck was the 10M serial `format!` calls reconstructing SRA/Illumina read IDs, which took 24s and now takes ~4.5s.

### Changed

- **`QualityCompressor::Auto` is the new default; explicit `Bsc` is now honored** — Previously, the lossless quality path silently overrode an explicit `QualityCompressor::Bsc` to `QualityCtx` whenever a chunk had ≥100k records. The override is now explicit: `Auto` (the new default) does the threshold-based selection, and any explicit choice (including `Bsc`) is respected. Default archives still use `QualityCtx` for ≥100k reads; the user-visible change is only for users who explicitly request `Bsc` via `--config`.
- **`CompressConfig` simplified** — Production-only fields remain in `cli.rs`; experimental flags live in `AdvancedConfig` or benchmark code.
- **Ultra level table tuned** — BSC block sizes now scale with ultra level (25→100 MB) for progressively better BWT context. Level 2 chunk size increased from 2M to 2.5M reads. Level 4 reduced from 2 to 1 parallel chunk for better memory efficiency. Each level now specifies `bsc_block_mb` and `max_seq_block_mb` instead of using hardcoded values.
- **`--local-reorder` uses ReorderLocal** — Previously used Delta encoding (encoding_type=8) which was consistently worse than ReorderLocal with identity reorder. Now maps to ReorderLocal with auto-selected ultra level, matching `--ultra 0` behavior. Existing Delta-encoded archives (encoding_type=8) remain decompressible.

### Fixed

- **`qz decompress -t` had no effect on rayon parallelism** — The `--threads` flag was stored on `DecompressConfig` but never installed as the rayon thread pool, so BSC block decompression and record reconstruction always ran on the global pool (typically all available cores). `qz decompress` and `qz verify` now build a rayon pool sized to `-t` before doing any parallel work.
- **`bsc_block_size_mb` could produce undecompressible archives** — Compression accepted any value but the streaming decompressor rejects compressed blocks larger than 64 MiB. Compression now rejects `bsc_block_size_mb` outside `1..=64` so this can no longer happen.
- **Per-read `Vec<u8>` allocation in `quality_ctx`** — The `precompute_base_indices` helper allocated a fresh `Vec` for every read in both compress and decompress hot loops (millions of small allocations on 10M-read inputs). Replaced with direct LUT lookup in the inner loop; the range coder serializes the work anyway, so there was no vectorisation benefit to lose.
- **`DECOMPRESS_BATCH_SIZE` scales with thread count** — The streaming decompressor batched 8 BSC blocks per stream regardless of CPU count, leaving large machines underused. Batch size is now `current_num_threads().div_ceil(3)` clamped to `[8, 32]`, saturating ~24 cores per stream's batch.
- **`--ultra 0` auto-select capped at level 3** — On 10M-read Illumina WGS, ultra levels 4 and 5 produce identical-size archives to level 3 (8.48x) but decompress 2-3x slower (61s vs 23s on 72 cores) because larger BSC blocks reduce decompress parallelism. Auto-select previously chose level 5 on 32+ core machines; it now caps at level 3, which is the empirical sweet spot. Levels 4 and 5 remain available via explicit `--ultra 4` / `--ultra 5`.
- **Parallel compression temp file collision** — Running multiple `qz compress` processes simultaneously in the same working directory would crash with "No such file or directory" because all processes wrote to identically-named temp files (`.qz_harc_*.tmp`, `.qz_chunked_*.tmp`, `.qz_reorder_*.tmp`). Temp file names now include the process ID to prevent collisions.
- **Silent data corruption on malformed archives** — `parse_order` and `read_varint` in ultra decompression silently returned 0 on parse errors, producing garbage output from corrupted archives. Now propagates errors properly.
- **Archive version not validated** — Decompression accepted archives from any future version without error. Now rejects archives with version > current with a clear "please update qz" message.
- **Zstd decompression 100 MB limit** — `decompress_zstd` had a hardcoded 100 MB cap, causing decompression failure for Zstd-compressed sequences > 100 MB. Switched to streaming decoder with no size limit.
- **Unnecessary ~1.5 GB quality clone** — The in-memory compression path cloned all quality data even when quality mode was Lossless (no-op quantization). Now skips the copy.
- **Missing BufWriter in compress_in_memory and write_chunked_archive** — Output files were written with many small unbuffered syscalls. Now uses BufWriter for all file output paths.
- **Missing stdout support in compress_in_memory** — The in-memory compression path created a literal file named "-" instead of writing to stdout. Now checks for stdio path like other compression paths.
- **Missing shrink_to_fit on ultra seq blocks** — Compressed BSC blocks in the ReorderLocal path retained allocated capacity (BSC allocates output = input + header), wasting memory at scale.
- **quality_ctx always enabled in ultra** — The condition `|| true` made quality compressor selection dead code in ultra mode. Now explicitly documents that quality_ctx is always used for lossless ultra.
- **Ultra level 4 never auto-selected** — Core-count thresholds jumped from level 3 (cores < 16) to level 5 (cores >= 16), skipping level 4. Added level 4 for 16-31 cores.
- **filter_map/map mismatch in quality_ctx chunked path** — `filter_map` on quality refs skipped records without quality, causing misalignment with sequence refs (which used `map`). Now uses `map` with empty default for both.
- **FASTA mode ignored in compress_in_memory** — The in-memory compression path always parsed input as FASTQ regardless of the `--fasta` flag.
