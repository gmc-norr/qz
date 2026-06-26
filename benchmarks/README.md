# qz / bz mode benchmarks

End-to-end benchmarks for the production modes, each measuring **compressed size,
ratio, compress & decompress wall-time, peak RSS, CPU%, and a roundtrip
correctness check**. Competitors are genomics-specific: **SPRING** and **pigz**
for FASTQ, **samtools CRAM** for BAM.

## Coverage tiers & datasets

| tier | dataset | used by |
|---|---|---|
| **shallow** single-end | `real_data/ERR3239334_1.30m.fastq` (30M reads, 150 bp WGS) | `bench_single_end.sh` |
| **deep** paired / reference | `real_data/HG002.chr20.novaseq.2x151.R{1,2}.fastq` (chr20 2×151), reference `real_data/chr20.fa` | `bench_paired_end.sh`, `bench_reference.sh` |
| **full WGS** BAM | `HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam` (~51 GB) | `bench_bam.sh` |
| **chr20** BAM | `real_data/HG002.GRCh38.chr20.bam` | `bench_bam.sh` |

The shallow tier is single-end only (fast, broad codec sweep); the deep tier is
where paired + reference modes are measured. CRAM/reference runs use GRCh38
references; `bench_reference.sh` asserts `chr20.fa`'s chr20 M5 equals the GRCh38
chr20 checksum before running, so a mismatched reference fails loudly.

## Scripts

| script | mode | competitor |
|---|---|---|
| `bench_single_end.sh` | single-end FASTQ (shallow 30M) | SPRING lossless, pigz -9 |
| `bench_paired_end.sh` | paired-end FASTQ (deep chr20) | SPRING paired |
| `bench_reference.sh` | reference-based FASTQ (deep chr20, single + paired) | SPRING single + paired |
| `bench_bam.sh` | BAM (full-WGS headline + chr20 level sweep) | samtools CRAM 3.0 / 3.1 |
| `run_all_modes.sh` | drives all four | — |
| `lib/bench_common.sh` | shared helpers (`emit_row`, `cluster_verify`, …) | — |

## Run

```bash
# all four modes, 72 threads (default)
./run_all_modes.sh

# a subset
./run_all_modes.sh single bam

# override threads / inputs / binaries
T=36 ./bench_single_end.sh
INPUT=/path/reads.fastq QZ=/path/qz ./bench_single_end.sh
```

Results land in `results/<mode>.tsv` (machine-readable) + `results/<mode>.log`
(full `/usr/bin/time -v` reports + tool stdout). `run_all_modes.sh` reads the four
`results/{single,paired,reference,bam}.tsv` files and prints a combined `column -t`
table at the end. The `results/` TSV/log artifacts are gitignored.

## qz/bz variants benchmarked

- **single (shallow):** `default`, `--fast`, `--ultra 1/2/3`, `--cluster` (order
  drop), and a `--numa auto` row paired against `default` (which runs `--numa off`)
  to isolate the multi-socket sharding speedup.
- **paired (deep):** `--cluster` only — deep data skips the plain `default` paired
  row this pass; the cluster row is the paired headline.
- **reference (deep):** single-end (R1) `--reference` / `--reference --reference-fast`
  / `--cluster`, and paired `--reference` / `--reference --reference-fast`. An
  untimed `qz index --like <R1>` preflight builds both index variants keyed to the
  read length so compress finds them.
- **bam (bz):** full-WGS `default` + lossy `--reduce-quality 1/2/3`; chr20 `default`
  + `-l 1/2/3` (level presets are a memory/speed dial, ratio is ~flat).

Codec-comparison rows run `--numa off` explicitly (qz and bz both default to
`--numa auto`), so size/ratio reflect the codec, not sharding.

## Roundtrip oracles

The `roundtrip` column uses a different oracle per variant:

- **md5 of the reconstructed FASTQ** → `OK` / `MISMATCH` (qz default/fast/ultra/
  reference, SPRING, pigz).
- **`qz verify`** (`cluster_verify`) → `PASS` / `FAIL` for `--cluster`: cluster
  output is a *permutation* of the input, so md5 doesn't apply; `qz verify` checks
  the embedded multiset checksum (set-losslessness).
- **`samtools view | md5sum` fingerprint** for bz: lossless rows compare the full
  record stream (`OK`/`MISMATCH`); lossy `--reduce-quality` rows compare a
  **QUAL-stripped** fingerprint (`$11="*"`) → `QUAL_ONLY_DIFF` (everything but
  quality is byte-identical) / `NONQUAL_MISMATCH`.
- **`n-a`** for CRAM (ratio/speed competitor; its roundtrip isn't re-verified).

## TSV schema

`variant, lossless, comp_bytes, raw_bytes, ratio, comp_s, comp_rss_gb, comp_cpu,
decomp_s, decomp_rss_gb, decomp_cpu, roundtrip`.

`raw_bytes` is the per-row denominator for `ratio`: R1 bytes for single-end rows,
R1+R2 for paired rows, the input BAM bytes for bz/CRAM rows.

## Notes / gotchas baked into the harness

- qz/bz refuse to overwrite an existing `-o`; the scripts pass `-f`.
- **Paired/reference decompress:** `-o <prefix>` emits `<prefix>_R1.fastq` +
  `<prefix>_R2.fastq`; roundtrip md5-checks both mates. **Paired-cluster decode**
  takes two explicit `-o` paths instead of the prefix contract.
- **`qz --reference` works for single-end (one input) and paired (R1+R2).** The
  reference is used only at compress time — never stored, never needed to decode.
- **bz roundtrip is record-level lossless but re-frames BGZF** → not byte-identical
  to the input BAM; correctness is the `samtools view` fingerprint, not file md5.
- **bz needs `@HD SO:coordinate`.** `prep_fingerprints` reheaders a mislabeled input
  (fast, header-only) into a cached copy; both benchmark BAMs are already
  coordinate-sorted, so this is normally a no-op.
- `bench_bam.sh` runs a free-space preflight (`MIN_FREE_GB`, default 200) and warns
  before the ~51 GB full-WGS rows.
- All temp archives/outputs are written under `/tmp/qzbench_*` and cleaned up; only
  the small TSV/log files remain.
