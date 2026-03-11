# QZ - High-Performance FASTQ Compression

## Build & Environment
- Edition 2024 requires **nightly Rust**: `rustup install nightly`
- C++ dependency: `third_party/libbsc/` (clone from https://github.com/IlyaGrebnov/libbsc.git)
- Cargo workspace with 4 crates: `qz-lib`, `qz-cli`, `qz-bench`, `qz-python`
- Build: `cargo build --release` (build.rs in qz-lib compiles libbsc as static C++ lib, linked via FFI)
- Release profile: LTO "fat", 1 codegen unit, panic="abort", stripped

## Workspace Structure

```
qz/
├── Cargo.toml                   (workspace root)
├── third_party/                 (libbsc, htscodecs)
├── crates/
│   ├── qz-lib/                  (library: all algorithms)
│   │   ├── Cargo.toml
│   │   ├── build.rs             (compiles libbsc + htscodecs)
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── cli.rs           (CompressConfig, DecompressConfig, enums — NO Clap)
│   │   │   ├── compression/     (all modules)
│   │   │   └── io/              (fastq.rs)
│   │   └── tests/
│   │       └── integration_test.rs
│   ├── qz-cli/                  (binary: slim CLI → produces `qz` executable)
│   │   ├── Cargo.toml
│   │   └── src/main.rs
│   ├── qz-python/               (PyO3 cdylib, built with maturin)
│   │   ├── Cargo.toml
│   │   ├── pyproject.toml
│   │   └── src/lib.rs
│   └── qz-bench/                (bench binaries)
│       ├── Cargo.toml
│       └── src/bin/             (14 bench_*.rs files)
├── benchmarks/                  (benchmark scripts and results)
└── real_data/                   (test data)
```

## Architecture Overview

```
FASTQ input
    ├─ Headers   → Raw + BSC (BWT + LZP + adaptive QLFC)
    ├─ Sequences → Raw + BSC (block-parallel, 25 MB blocks)
    └─ Qualities → BSC (with optional delta/modeling/dictionary)
         ↓
    QZ archive (single file with metadata + block offsets)
```

All three streams compress **in parallel** via `rayon::join` (headers parallel with nested seq+qual).

### Key Files
- `crates/qz-lib/src/cli.rs` - CompressConfig, DecompressConfig, enums (no Clap)
- `crates/qz-lib/src/compression/mod.rs` - Main orchestrator (~2400 lines), compress() and decompress()
- `crates/qz-lib/src/compression/bsc.rs` - FFI bindings to libbsc; `compress_parallel` / `decompress_parallel` split into 25 MB blocks
- `crates/qz-lib/src/io/fastq.rs` - `FastqRecord { id, sequence, quality: Option<String> }`
- `crates/qz-lib/tests/integration_test.rs` - 30 roundtrip tests
- `crates/qz-cli/src/main.rs` - Clap CLI (only production features: default, ultra, illumina-bin, discard)

### Module Visibility
- **Private** (internal, access via `super::`): columnar, delta, n_mask, quality, quality_delta, quality_model, read_id, rle, zstd_dict, factorize, ultra
- **Public**: arithmetic_quality, arithmetic_sequence, bsc, debruijn, dna_utils, fqzcomp, greedy_contig, header_col, openzl, paired_end, quality_context, quality_ctx, template

## Key Design Decisions

### BSC is the default compressor for all streams
- BSC adaptive (BWT + LZP + adaptive QLFC) beats Zstd for structured genomic data
- Headers: Raw+BSC is ~2x better than template+Zstd
- Sequences: Raw+BSC gives ~4.6x (order-preserving); ~1.85 bits/base is near theoretical limit without reordering
- BSC threading: ONLY use `LIBBSC_FEATURE_FASTMODE` (never `MULTITHREADING`) -- rayon handles parallelism

### Block-parallel BSC
- Data split into 25 MB blocks, each compressed/decompressed independently via rayon
- `compress_parallel` / `decompress_parallel` in bsc.rs handle block splitting + reassembly
- `compress_parallel_adaptive()` uses adaptive QLFC (slightly better, ~0.5-0.7%)

### Chunked streaming mode (`--chunked`)
- Reads in 5M-record chunks, BSC-compresses each chunk's streams into 25 MB blocks
- Pipelined: reads next chunk on main thread while compressing previous in `std::thread::scope`
- Streams compressed SEQUENTIALLY (headers -> seq -> qual) to limit BSC working memory
- Temp files (`.qz_chunked_{h,s,q}.tmp`) keep memory constant regardless of input size
- TmpCleanup drop guard ensures cleanup on success, error, or panic
- Peak memory: ~6 GB for 100M reads with 72 threads

### Memory lessons
- BSC allocates output buffer = input_size + header. Call `shrink_to_fit()` on compressed blocks
- Avoid `rayon::join` across all 3 streams simultaneously (70+ concurrent BWT = 14 GB). Process streams sequentially in chunked mode
- Don't accumulate `Vec<Vec<u8>>` blocks in memory for large inputs -- use temp files

### Read order is preserved (no reordering by default)
- SPRING reorders reads for ~5% better compression but destroys original order
- QZ keeps reads in input order for downstream tool compatibility

## Archive Format
```
[encoding_type: 1B] [flags: 1B (bit0=arithmetic, bit1=const_lengths)]
[read_lengths (if flags bit0)]
[quality_binning: 1B] [quality_compressor: 1B] [sequence_compressor: 1B] [header_compressor: 1B]
[quality_model (if enabled)] [quality_delta: 1B] [quality_dict (if enabled)]
[template_prefix + has_comment]
[num_reads: 8B] [headers_len: 8B] [sequences_len: 8B] [nmasks_len: 8B] [qualities_len: 8B]
[const_seq_len: 4B] [const_qual_len: 4B]  (only if flags bit1 set)
[headers data] [sequences data] [nmasks data] [qualities data]
```

Compressor codes: zstd=0, bsc=1. Header compressor: zstd=0 (template+zstd), bsc=1 (raw+bsc).
When flags bit1 is set, sequence and quality streams omit per-record varint length prefixes (lengths are constant).

## CompressConfig Gotcha
When adding new fields to `CompressConfig` in `crates/qz-lib/src/cli.rs`, you must update:
1. The `Default` impl in `cli.rs`
2. ALL constructors in `crates/qz-lib/tests/integration_test.rs` (12+ instances). Search for `CompressConfig {` to find them all.
3. The `into_config()` method in `crates/qz-cli/src/main.rs`

## Compression Performance (10M reads, 150bp WGS, Xeon Gold 6254 72 cores)
| Tool | Size (MB) | Ratio | Compress | Decompress |
|------|-----------|-------|----------|------------|
| QZ default | 427 | 8.17x | 15.6s | 35.1s |
| QZ ultra 1 | 418 | 8.35x | 26.6s | 44.0s |
| QZ ultra 3 | 408 | 8.55x | 41.2s | 32.8s |
| QZ ultra 5 | 408 | 8.56x | 46.0s | 1:02.6 |
| SPRING | 431 | 8.10x | 1:01.4 | 15.4s |
| bzip2 -9 | 542 | 6.44x | 2:47.8 | 1:26.6 |
| pigz -9 | 695 | 5.02x | 9.7s | 7.9s |

## Changelog
All user-visible changes (new features, breaking changes, removals, bug fixes) must be recorded in `CHANGELOG.md` at the project root. Update it as part of the same change, not as a separate step.

## Test Data
If available in `real_data/`:
- `NA12878_exome_50k.fastq` - 50K reads, 101bp exome
- `ERR3239334_1.1m.fastq` - 1M reads, 150bp WGS
- `ERR3239334_1.5m.fastq` - ~5M reads, 150bp WGS
