# QZ

Order-preserving FASTQ compression using columnar stream separation and block-sorting transforms.

## Overview

QZ decomposes each FASTQ record into three streams — headers, DNA sequences, and quality scores — and compresses each independently. The primary compression engine is [libbsc](https://github.com/IlyaGrebnov/libbsc) (Grebnov), which applies the Burrows-Wheeler Transform followed by Quantized Local Frequency Coding (QLFC). Separating the streams allows the BWT to operate on large, statistically homogeneous blocks: DNA sequences from the same organism share extensive k-mer content, and the BWT clusters these shared substrings, producing long runs amenable to entropy coding. On raw DNA sequences, this achieves ~1.85 bits per base without read reordering.

For quality scores, QZ implements a context-adaptive range coder (`quality_ctx`) that conditions on read position, previous quality value, a stability flag, and the local DNA sequence context. This exploits the strong positional and sequential correlations in Illumina quality profiles that block-sorting alone cannot capture, reducing quality stream size by ~7% relative to BSC.

All streams are split into blocks (25 MB for BSC, 500K reads for quality_ctx) and compressed in parallel via [rayon](https://github.com/rayon-rs/rayon). Input is read in 2.5M-record chunks with pipelined I/O. Read order is preserved.

## Results

10M reads, 150 bp, whole-genome sequencing (ERR3239334), 3,492 MB uncompressed. 72-core machine, sequential runs, byte-identical roundtrip verified via MD5.

| Tool | Size (MB) | Ratio | Compress | Comp RAM | Decompress | Dec RAM |
|------|-----------|-------|----------|----------|------------|---------|
| **QZ default** | 435 | **8.03x** | 17.4 s | 5.8 GB | 13.8 s | 6.9 GB |
| **QZ ultra 1** | 426 | **8.21x** | 28.4 s | 7.6 GB | 14.2 s | 8.7 GB |
| **QZ ultra 3** | 416 | **8.39x** | 37.1 s | 13.8 GB | 21.6 s | 8.5 GB |
| **QZ ultra 5** | 416 | **8.39x** | 31.4 s | 14.0 GB | 1:02.5 | 8.5 GB |
| SPRING | 431 | 8.10x | 1:01.4 | 11.9 GB | 15.4 s | 10.0 GB |
| bzip2 -9 | 542 | 6.44x | 2:47.8 | 7.3 MB | 1:26.6 | 4.5 MB |
| pigz -9 | 695 | 5.02x | 9.7 s | 20.8 MB | 7.9 s | 1.7 MB |

SPRING was run without `-r` (read order preserved). Raw timing data in [`benchmarks/results_10m_all_ultra.txt`](benchmarks/results_10m_all_ultra.txt).

## Architecture

### System overview

```
                            ┌──────────────────────────────────────────────┐
                            │              QZ Compression                  │
                            │                                              │
  FASTQ input ──►  FastqReader ──►  Chunk loop (2.5M records) ─────────────┤
  (file, .gz,      (buffered,       │                                      │
   or stdin)        auto-detect     │  ┌────────── rayon::join ─────────┐  │
                    gzip)           │  │                                │  │
                                    │  │  Headers ──► BSC blocks ──┐    │  │
                                    │  │                           │    │  │
                                    │  │  ┌── rayon::join ──────┐  │    │  │
                                    │  │  │                     │  │    │  │
                                    │  │  │ Sequences ► BSC ──┐ │  │    │  │
                                    │  │  │                   │ │  │    │  │
                                    │  │  │ Qualities ► BSC ──┤ │  │    │  │
                                    │  │  │  or quality_ctx   │ │  │    │  │
                                    │  │  └───────────────────┘ │  │    │  │
                                    │  └────────────────────────┘  │    │  │
                                    │                              │    │  │
                                    │         Compressed blocks ◄──┘    │  │
                                    │              │                    │  │
                                    │              ▼                    │  │
                                    │     Temp files or memory          │  │
                                    └──────────────┬────────────────────┘  │
                                                   │                       │
                                                   ▼                       │
                                          Archive assembly                 │
                                          (header + stream blocks)         │
                                                   │                       │
                                                   ▼                       │
                                              .qz archive ───────────────► │
                                          (file or stdout)                 │
                            └──────────────────────────────────────────────┘


                            ┌──────────────────────────────────────────────┐
                            │            QZ Decompression                  │
                            │                                              │
  .qz archive ──►  Parse archive header ──► Spawn 3 decompressor threads   │
  (file or           (magic, version,        │                             │
   stdin→tmpfile)     stream offsets)        │  Thread 1: Headers ──┐      │
                                             │    BSC decompress    │      │
                                             │    (batch=8 blocks)  │      │
                                             │                      │      │
                                             │  Thread 2: Seqs ─────┤      │
                                             │    BSC decompress    │      │
                                             │                      │      │
                                             │  Thread 3: Quals ────┤      │
                                             │    BSC decompress    │      │
                                             │                      │      │
                                             │     bounded channels │      │
                                             │     (capacity = 2)   │      │
                                             │          │           │      │
                                             │          ▼           │      │
                                             │  Main thread:        │      │
                                             │    Reconstruct       │      │
                                             │    FASTQ records     │      │
                                             │    from 3 streams    │      │
                                             │          │           │      │
                                             │          ▼           │      │
                                             │  FASTQ output ──────►│      │
                                             │  (file, .gz,         │      │
                                             │   or stdout)         │      │
                            └──────────────────────────────────────────────┘
```

### Compression pipeline

Compression proceeds in five stages: reading, stream building, parallel compression, block accumulation, and archive assembly.

**Stage 1: Chunked reading.** The FASTQ reader (`FastqReader`) reads records in chunks of 2.5M (default mode) or 5M (ultra mode). Input may be a file, gzipped file (auto-detected via magic bytes), or stdin. Each chunk yields a `Vec<FastqRecord>` where records store raw bytes (`Vec<u8>`) to avoid UTF-8 validation overhead. Reading is pipelined: the main thread reads the next chunk while a background `std::thread::scope` compresses the current one.

```
Main thread:     ┃ Read chunk 0 ┃ Read chunk 1 ┃ Read chunk 2 ┃ ...
                 ┃              ┃              ┃              ┃
Compress thread: ┃              ┃ Compress 0   ┃ Compress 1   ┃ ...
                                ╰──overlapped──╯
```

**Stage 2: Stream building.** Each chunk's records are split into three byte streams by `records_to_streams()`:

| Stream | Format | Content |
|--------|--------|---------|
| Headers | `[varint(len), raw bytes]...` | `@SRR...` identifier lines |
| Sequences | `[varint(len), bases]...` | `ACGTN...` DNA bases |
| Qualities | `[varint(len), packed bytes]...` | Phred+33 scores, optionally bit-packed |

When all reads have the same length (common for Illumina), per-read varint framing is omitted and the constant length is stored once in the archive header, saving ~1 byte per read.

**Stage 3: Parallel compression.** The three streams are compressed concurrently using nested `rayon::join`:

```
rayon::join(
    ║ Headers ──► BSC (25 MB blocks, adaptive QLFC, no LZP)
    ║
    ║ rayon::join(
    ║     ║ Sequences ──► BSC (25 MB blocks)
    ║     ║ Qualities ──► BSC or quality_ctx (500K-read sub-blocks)
    ║ )
)
```

Each stream is split into 25 MB blocks that are compressed independently via `rayon::par_iter`. BSC uses the Burrows-Wheeler Transform followed by adaptive Quantized Local Frequency Coding (QLFC). LZP (Lempel-Ziv Prediction) is disabled in the default path because the BWT already captures repeating patterns; LZP adds 5–10% runtime with negligible compression gain on genomic data.

Compressed blocks are shrunk via `shrink_to_fit()` immediately after compression. Without this, BSC's output buffer allocation (input size + header per block) would waste ~35 GB across 1400 blocks.

**Stage 4: Block accumulation.** Compressed blocks are accumulated either in memory (`Vec<Vec<u8>>` for small inputs) or streamed to temp files (for large inputs or reorder mode). Temp files use a RAII cleanup guard (`TmpCleanup`) that deletes files on drop, including on panic or error.

**Stage 5: Archive assembly.** The archive is written sequentially: v2 header, then header blocks, sequence blocks, quality blocks. No seeking is required, so output can go to stdout. The archive format is:

```
┌──────────────────────────────────────────────────────────────────┐
│  V2 prefix (8 bytes)                                             │
│  ┌──────────┬─────────┬──────────┬────────────────┐              │
│  │ "QZ"     │ ver=02  │ rsvd=00  │ header_size u32│              │
│  └──────────┴─────────┴──────────┴────────────────┘              │
│                                                                  │
│  Header body (variable length)                                   │
│  ┌────────────────────────────────────────────────────────┐      │
│  │ encoding_type, flags, quality_binning                  │      │
│  │ quality_compressor, sequence_compressor, header_comp   │      │
│  │ num_reads (u64), stream lengths (u64 x 4)              │      │
│  │ const_seq_len, const_qual_len (if flag set)            │      │
│  └────────────────────────────────────────────────────────┘      │
│                                                                  │
│  Stream data                                                     │
│  ┌────────────────────────────────────────────────────────┐      │
│  │ Headers:   num_blocks u32, [len u32, data]...          │      │
│  │ Sequences: num_blocks u32, [len u32, data]...          │      │
│  │ Qualities: num_blocks u32, [len u32, data]...          │      │
│  └────────────────────────────────────────────────────────┘      │
└──────────────────────────────────────────────────────────────────┘
```

The header is self-describing: all compressor types, encoding modes, and stream lengths are recorded, so decompression requires no external metadata.

### Decompression pipeline

Decompression uses a streaming architecture with three parallel decompressor threads feeding a single reconstruction thread through bounded channels.

**Header parsing.** The archive header is read and validated (magic bytes `QZ`, version `0x02`). Stream offsets are computed from the recorded lengths: headers start at `data_offset`, sequences at `data_offset + headers_len`, qualities at `data_offset + headers_len + sequences_len`.

**Parallel decompression.** Three background threads are spawned via `std::thread::scope`, each responsible for one stream:

1. Seek to stream offset in the archive file
2. Read `num_blocks` from the stream header
3. In batches of 8: read block headers and data sequentially, decompress the batch in parallel via `rayon::par_iter`, send decompressed blocks through a bounded `SyncSender` channel

The channel capacity is 2 blocks per stream. This backpressures the decompressor threads when the reconstruction thread falls behind, bounding peak memory to ~300 MB (3 streams x ~4 blocks x 25 MB + output buffer).

**Record reconstruction.** The main thread reads from all three channels through `ChannelStreamBuffer` wrappers that provide varint/byte-slice reading over the channel. For each of the `num_reads` records, it reads the header length and bytes, sequence length and bytes, and quality length and packed bytes, unpacks qualities to ASCII, and writes the four FASTQ lines (`@header\nseq\n+\nqual\n`).

**Stdin input.** Since the decompressor needs to seek to three different offsets in the archive, stdin input is spooled to a temp file first, then decompressed normally.

**Gzip output.** When `--gzipped` is requested, output is piped through `gzp::ParCompress` for parallel gzip compression (multi-threaded block compression into a multi-member gzip stream).

### BSC block compression

[libbsc](https://github.com/IlyaGrebnov/libbsc) applies the Burrows-Wheeler Transform (via [libsais](https://github.com/IlyaGrebnov/libsais)) to sort all rotations of the input block, clustering repeated substrings. The BWT output is then compressed by adaptive QLFC, which models local symbol frequencies. The block size is 25 MB.

QZ compiles libbsc and libsais from source via `build.rs` with `-O3 -march=native` and OpenMP support. The complete libbsc pipeline is available: BWT, LZP (Lempel-Ziv Prediction), QLFC (static and adaptive), sort transform, and preprocessing filters. QZ drives it through Rust FFI bindings (`bsc.rs`).

```
Stream (75 MB)
    ├─ Block 0 (25 MB) ──► BWT ──► QLFC ──► compressed (rayon worker 1)
    ├─ Block 1 (25 MB) ──► BWT ──► QLFC ──► compressed (rayon worker 2)
    └─ Block 2 (25 MB) ──► BWT ──► QLFC ──► compressed (rayon worker 3)
```

The default configuration uses adaptive QLFC without LZP (`compress_parallel_adaptive_no_lzp`). Adaptive QLFC learns symbol frequencies during encoding, improving compression by 0.5–0.7% over the static variant. LZP is disabled because the BWT already captures the repeating k-mer structure in genomic data; LZP adds 5–10% runtime with negligible compression benefit.

**Two-level threading model.** libbsc has its own OpenMP-based multithreading for parallelizing the BWT within a single block. QZ adds a second level of parallelism above it using rayon for inter-block parallelism. These two levels interact carefully:

```
                    ┌───────────────────────────────────────┐
                    │       rayon thread pool               │
                    │       (inter-block parallelism)       │
                    │                                       │
                    │  Worker 1: Block 0                    │
                    │    └─ bsc_compress(FASTMODE)          │
                    │       └─ BWT (single-threaded)        │
                    │       └─ QLFC adaptive                │
                    │                                       │
                    │  Worker 2: Block 1                    │
                    │    └─ bsc_compress(FASTMODE)          │
                    │       └─ BWT (single-threaded)        │
                    │       └─ QLFC adaptive                │
                    │                                       │
                    │  Worker N: Block N                    │
                    │    └─ ...                             │
                    └───────────────────────────────────────┘

                    vs.

                    ┌───────────────────────────────────────┐
                    │       Single block, MT mode           │
                    │       (intra-block parallelism)       │
                    │                                       │
                    │  bsc_compress(FASTMODE|MULTITHREADING)│
                    │    └─ BWT via libsais_bwt_omp()       │
                    │       └─ up to 16 OpenMP threads      │
                    │    └─ QLFC adaptive                   │
                    └───────────────────────────────────────┘
```

libbsc supports two feature flags passed per call:

| Flag | Effect |
|------|--------|
| `FASTMODE` | Enables fast decompression; skips expensive data-type detection |
| `MULTITHREADING` | Enables OpenMP parallelism within BWT (libsais) |

The default compression path uses **FASTMODE only**: each 25 MB block gets a single-threaded BWT, but rayon runs many blocks concurrently. This avoids thread explosion (N rayon workers x M OpenMP threads) and gives better throughput than intra-block parallelism for typical genomic workloads with many blocks.

The MT path (`FASTMODE | MULTITHREADING`) is available for single large blocks where intra-block parallelism matters more. When used, OpenMP threads are capped at 12 via `omp_set_num_threads()`.

**libbsc modification.** QZ patches one file in the libbsc source: `bwt.cpp`. Upstream libbsc caps the thread count passed to `libsais_bwt_omp()` at 8. QZ raises this cap to 16 to better utilize high-core machines:

```c
// upstream:  numThreads > 8 ? 8 : numThreads
// QZ:        numThreads > 16 ? 16 : numThreads
```

This affects both `libsais_bwt_aux_omp()` and `libsais_bwt_omp()` calls. On a 72-core system, the effective thread count per BWT call is `min(omp_get_max_threads() / omp_get_num_threads(), 16)`, allowing better scaling when the MT path is used for large individual blocks.

### Quality score compression

QZ automatically selects the quality compressor based on input:

| Condition | Compressor | Rationale |
|-----------|-----------|-----------|
| Lossless, >= 100K reads | `quality_ctx` | Context model outperforms BSC by ~7% once converged |
| Lossless, < 100K reads | BSC | Insufficient data for context model convergence |
| Lossy modes | BSC | Reduced alphabet (3–7 bits) compresses well under BWT |

**Context model.** `quality_ctx` is an LZMA-style forward range coder with 160,000 adaptive contexts. Each quality symbol is coded in a context formed by the cross-product of:

| Feature | States | Description |
|---------|--------|-------------|
| Position bin | 64 | `floor(pos / 5)`, captures positional quality curve |
| Previous quality | 50 | Phred score of preceding base (0–49) |
| Stability | 2 | `\|q[i-1] - q[i-2]\| <= 2` (stable) vs. changing |
| Base pair | 25 | `prev_base * 5 + cur_base` (A/C/G/T/N x A/C/G/T/N) |

Each context maintains a Laplace-smoothed frequency table (~20 symbols for Illumina), updated per symbol and rescaled when the total exceeds 2^20. For large inputs, quality scores are partitioned into 500K-read sub-blocks compressed independently in parallel.

On 150 bp Illumina WGS, `quality_ctx` achieves ~0.15 bits per quality score. Quality data accounts for ~6% of total archive size in lossless mode.

**Lossy quantization.** Quality scores may be quantized before compression: Illumina 8-level binning maps ~40 Phred values to 8 representatives (3 bits/symbol); discard mode replaces all scores with a constant. Quantized scores are bit-packed to the minimum width.

**Relation to prior work.** The context model draws on fqzcomp (Bonfield, 2013; [paper](https://doi.org/10.1093/bioinformatics/btac010)), which demonstrated quality-conditioned arithmetic coding and was adopted as the CRAM 3.1 quality codec, and ENANO (Dufort y Alvarez et al., 2020; [paper](https://doi.org/10.1093/bioinformatics/btaa551)), which introduced DNA sequence context and stability tracking for nanopore data. QZ adapts both for Illumina short reads with a 2-base window and binary stability flag, yielding 160K contexts that converge fast on the narrow Illumina quality distribution.

### Ultra mode

Ultra mode (`--ultra [1-5]`) increases compression by reordering reads to group similar sequences together before compression. Levels control chunk size and parallelism:

| Level | Chunk size | Parallel chunks | Quality sub-block | RAM |
|-------|-----------|-----------------|-------------------|-----|
| 1 | 1M | 4 | 250K | ~8 GB |
| 3 | 5M | 2 | 500K | ~14 GB |
| 5 | 10M | 1 | 1M | ~14 GB |
| auto | varies | varies | varies | fits available RAM |

**Reordering strategy.** Reads are grouped by sequence similarity using center-hash grouping: two 32-bit hashes are computed from the central 32 bases of each read, and reads sharing a hash are placed adjacent in the output. Singletons (unique hashes) remain in input order to preserve BWT-friendly locality. The permutation is stored in the archive so the decompressor can restore original order.

### Memory management

Memory is the primary constraint for high-throughput genomic compression. QZ uses several strategies to keep peak usage bounded:

- **Sequential stream compression.** Within each chunk, headers, sequences, and qualities are compressed sequentially (not all three simultaneously). This prevents 70+ concurrent BWT allocations that would consume ~14 GB.
- **Pipelined I/O.** Reading the next chunk overlaps with compressing the current one, hiding I/O latency without doubling memory.
- **Block shrinking.** `shrink_to_fit()` on each BSC output block releases the unused allocation headroom (BSC allocates output = input + header).
- **Temp file accumulation.** For large inputs, compressed blocks are streamed to disk rather than held in memory. A RAII drop guard ensures cleanup on success, error, or panic.
- **Bounded decompression channels.** Decompressor threads send blocks through channels with capacity 2, preventing unbounded memory growth when the writer is slower than decompression.

## QZ File Format

The QZ archive is a self-describing binary format. All integers are little-endian. The archive consists of a fixed-size prefix, a variable-length header body, and concatenated compressed stream data.

### File layout

```
┌─────────────────────────────────────────────────────────────────┐
│  V2 Prefix (8 bytes, fixed)                                     │
│  Header Body (variable, starts at byte 8)                       │
│  Stream Data (starts at byte header_size)                       │
│    ├── Headers stream   (multi-block)                           │
│    ├── Sequences stream (multi-block)                           │
│    ├── N-masks stream   (multi-block, 0 bytes if unused)        │
│    └── Qualities stream (multi-block)                           │
└─────────────────────────────────────────────────────────────────┘
```

### V2 prefix (8 bytes)

```
Offset  Size  Type    Field         Value
------  ----  ------  -----------   -----
0       2     u8[2]   magic         "QZ" (0x51 0x5A)
2       1     u8      version       2
3       1     u8      reserved      0
4       4     u32     header_size   Total bytes of prefix + body.
                                    Stream data starts at this offset.
```

### Header body (variable length, starts at byte 8)

```
Offset  Size  Type    Field                Description
------  ----  ------  -------------------  -----------
+0      1     u8      encoding_type        Sequence encoding mode (see table)
+1      1     u8      flags                Bit 0: arithmetic (legacy)
                                           Bit 1: const-length fields present
+2      1     u8      quality_binning      0=lossless, 1=illumina-8, 2=binary, 3=four-level
+3      1     u8      quality_compressor   0=zstd, 1=BSC, 2=openzl, 3=fqzcomp, 4=quality_ctx
+4      1     u8      seq_compressor       0=zstd, 1=BSC, 2=openzl
+5      1     u8      header_compressor    0=zstd, 1=BSC, 2=openzl
+6      1     u8      quality_model_flag   0=disabled, 1=enabled (followed by model data)
+7      1     u8      quality_delta_flag   0=disabled, 1=enabled
+8      1     u8      quality_dict_flag    0=absent, non-zero=present (followed by dict)
+9      2     u16     template_prefix_len  Read ID template prefix length
+11     N     bytes   template_prefix      Prefix bytes (N = template_prefix_len)
+11+N   1     u8      has_comment          0=no common comment
+12+N   8     u64     num_reads            Total FASTQ records
+20+N   8     u64     headers_len          Compressed headers stream length (bytes)
+28+N   8     u64     sequences_len        Compressed sequences stream length
+36+N   8     u64     nmasks_len           N-mask stream length (0 if unused)
+44+N   8     u64     qualities_len        Compressed qualities stream length
```

When `flags & 0x02` (const-length fields present), two additional fields follow:

```
+52+N   4     u32     const_seq_len        Constant read length (0 = variable)
+56+N   4     u32     const_qual_len       Constant quality length (0 = variable)
```

When constant lengths are set, per-record varint length prefixes are omitted in the corresponding stream, saving ~1 byte per read.

### Encoding types

| Code | Name | Description |
|------|------|-------------|
| 0 | Raw | No special sequence encoding (default) |
| 4 | SequenceHints | 1-byte syncmer hint per record before sequence data |
| 5 | SequenceDelta | Inline delta encoding against cached similar reads |
| 6 | RcCanon | Reverse-complement canonicalized; RC flags stream appended after qualities |
| 8 | HarcDelta | Ultra mode: grouped reads + previous-read delta, stores permutation |
| 9 | ReorderLocal | Ultra mode: grouped reads + raw BSC, stores permutation |

### Multi-block stream format

Each compressed stream (headers, sequences, qualities) uses block framing:

```
+0      4     u32     num_blocks     Number of compressed blocks
Per block:
  +0    4     u32     block_len      Compressed block size in bytes
  +4    N     bytes   block_data     Compressed payload (BSC, quality_ctx, etc.)
```

BSC blocks contain up to 25 MB of uncompressed data each. The `headers_len`, `sequences_len`, and `qualities_len` values in the header include the `num_blocks` prefix and all `(4 + block_len)` entries.

### Per-record stream internals (after decompression)

Within the decompressed payload of each stream, records are framed with unsigned LEB128 varint length prefixes (omitted when constant lengths are set):

| Stream | Per-record format |
|--------|------------------|
| Headers | `[varint(len)] [header bytes]` — ID line without `@` prefix |
| Sequences | `[varint(len)] [ASCII bases]` — A/C/G/T/N |
| Qualities | `[varint(len)] [bit-packed quals]` — packed to minimum bits per binning mode |

For `quality_ctx` compression, the quality stream is compressed as a whole (one or more sub-blocks for 500K reads each) without per-record framing — the context model implicitly tracks record boundaries.

## Installation

### Requirements

- Rust nightly (edition 2024)
- C++ compiler with OpenMP support (for libbsc)

### Build

```bash
git clone <repo-url> qz && cd qz
git clone https://github.com/IlyaGrebnov/libbsc.git third_party/libbsc
rustup install nightly
rustup run nightly cargo build --release
# Binary: target/release/qz
```

### Python bindings

```bash
cd crates/qz-python
pip install maturin
maturin develop --release
```

## Usage

### Compress

```bash
qz compress -i reads.fastq -o reads.qz                          # lossless (default)
qz compress -i reads.fastq.gz -o reads.qz                       # gzipped input (auto-detected)
qz compress -i reads.fastq -o reads.qz --quality-mode illumina-bin  # lossy 8-level binning
qz compress -i reads.fastq -o reads.qz --quality-mode discard   # discard quality scores
qz compress -i reads.fastq -o reads.qz --ultra 3                # ultra compression, level 3
qz compress -i seqs.fasta -o seqs.qz --fasta                    # FASTA input
```

### Decompress

Decompression auto-detects all settings from the archive header.

```bash
qz decompress -i reads.qz -o reads.fastq
qz decompress -i reads.qz -o reads.fastq.gz --gzipped           # gzipped output
```

### Piping (stdin/stdout)

Use `-` for `-i` or `-o` to read from stdin or write to stdout:

```bash
cat reads.fastq | qz compress -i - -o reads.qz                  # compress from stdin
qz compress -i reads.fastq -o - > reads.qz                      # compress to stdout
qz decompress -i reads.qz -o - > reads.fastq                    # decompress to stdout
cat reads.qz | qz decompress -i - -o reads.fastq                # decompress from stdin
cat reads.fastq | qz compress -i - -o - | qz decompress -i - -o - > out.fastq  # full pipe
```

All log output goes to stderr, so stdout remains clean for piped data. Decompression from stdin spools the archive to a temp file first (the decompressor needs to seek to stream offsets).

### CLI reference

**Compress:**

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input FILE` | Input FASTQ file (gzipped auto-detected, `-` for stdin) | required |
| `-o, --output FILE` | Output QZ archive (`-` for stdout) | required |
| `-w, --working-dir PATH` | Working directory for temp files | `.` |
| `-t, --threads N` | Thread count (0 = auto) | auto |
| `--fasta` | Input is FASTA format | off |
| `-q, --quality-mode MODE` | `lossless`, `illumina-bin`, or `discard` | `lossless` |
| `--no-quality` | Equivalent to `--quality-mode discard` | off |
| `--ultra [LEVEL]` | Ultra compression (1–5, or omit for auto) | off |

**Decompress:**

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input FILE` | Input QZ archive (`-` for stdin) | required |
| `-o, --output FILE` | Output FASTQ file (`-` for stdout) | required |
| `-w, --working-dir PATH` | Working directory for temp files | `.` |
| `-t, --threads N` | Thread count | auto |
| `--gzipped` | Output gzipped FASTQ | off |
| `--gzip-level N` | Gzip level (0–9) | `6` |

### Python

```python
import qz

qz.compress("reads.fastq", "reads.qz")
qz.compress("reads.fastq", "reads.qz", ultra=3)
qz.compress("reads.fastq", "reads.qz", quality_mode="illumina-bin")

qz.decompress("reads.qz", "reads.fastq")
qz.decompress("reads.qz", "reads.fastq.gz", gzipped=True)

print(qz.version())
```

### Library (Rust)

```rust
use qz_lib::cli::{CompressConfig, QualityMode};
use qz_lib::compression;

let config = CompressConfig {
    input: vec!["reads.fastq".into()],
    output: "reads.qz".into(),
    quality_mode: QualityMode::Lossless,
    ..CompressConfig::default()
};
compression::compress(&config)?;
```

## Project Structure

```
qz/
├── Cargo.toml                     workspace root
├── third_party/
│   ├── libbsc/                    libbsc (BWT + QLFC), compiled via build.rs
│   │   └── libbsc/bwt/bwt.cpp     ← patched: BWT thread cap 8→16
│   └── htscodecs/                 htscodecs (fqzcomp quality codec), unmodified
├── crates/
│   ├── qz-lib/                    core library (all algorithms, no CLI deps)
│   │   ├── src/
│   │   │   ├── compression/
│   │   │   │   ├── mod.rs             archive format, I/O helpers, public API
│   │   │   │   ├── compress_impl.rs   chunked compression orchestrator
│   │   │   │   ├── decompress_impl.rs streaming decompression + header parsing
│   │   │   │   ├── codecs.rs          per-stream compress/decompress dispatch
│   │   │   │   ├── bsc.rs             libbsc FFI, block-parallel BSC, threading
│   │   │   │   ├── quality_ctx.rs     context-adaptive range coder (160K contexts)
│   │   │   │   ├── ultra.rs           ultra mode (reorder + quality_ctx)
│   │   │   │   ├── columnar.rs        quality binning + bit-packing
│   │   │   │   ├── fqzcomp.rs         htscodecs FFI for fqzcomp quality codec
│   │   │   │   ├── header_col.rs      columnar header compression
│   │   │   │   ├── n_mask.rs          2-bit DNA encoding + N-bitmap
│   │   │   │   ├── dna_utils.rs       k-mer hashing, reverse complement
│   │   │   │   └── ...                additional codec modules
│   │   │   ├── io/fastq.rs        FASTQ/FASTA reader (buffered, gzip, stdin)
│   │   │   └── cli.rs             CompressConfig, DecompressConfig (no Clap)
│   │   ├── build.rs               compiles libbsc + htscodecs as static C/C++ libs
│   │   └── tests/                 38 roundtrip integration tests
│   ├── qz-cli/                    CLI binary (Clap) → produces `qz` executable
│   ├── qz-python/                 Python bindings (PyO3/maturin)
│   └── qz-bench/                  development benchmark binaries
├── benchmarks/                    benchmark scripts and results
└── real_data/                     test data (not tracked)
```

### C/C++ dependencies

`build.rs` compiles two C/C++ libraries as static archives linked into the final binary:

**libbsc** — Block-sorting compressor. All source files are compiled with `-O3 -march=native -std=c++11 -fopenmp`. The `LIBBSC_OPENMP_SUPPORT` and `LIBSAIS_OPENMP` defines enable OpenMP-parallel BWT/suffix-array construction. One source file is patched (see [BSC block compression](#bsc-block-compression) above).

**htscodecs** — Only `fqzcomp_qual.c` and `utils.c` are compiled (not the full library). These provide the fqzcomp quality compression algorithm, available as an alternative quality backend. No modifications to upstream source.

## Testing

```bash
rustup run nightly cargo test --release
```

131 tests total: 88 unit tests (per-module codec roundtrips) and 43 integration tests covering lossless, lossy quality modes, ultra mode, FASTA, arithmetic/de Bruijn/delta/RLE encodings, error handling on corrupt archives, and edge cases.

## Environment Variables

| Variable | Effect |
|----------|--------|
| `QZ_NO_BANNER=1` | Suppress version banner on stderr |
| `RUST_LOG=debug` | Enable debug logging (tracing) |

## Acknowledgments

- [libbsc](https://github.com/IlyaGrebnov/libbsc) by Ilya Grebnov — block-sorting compression (BWT + LZP + QLFC). Apache 2.0.
- [libsais](https://github.com/IlyaGrebnov/libsais) by Ilya Grebnov — suffix array construction used by libbsc's BWT. Apache 2.0.
- [fqzcomp](https://github.com/jkbonfield/fqzcomp) by James Bonfield — context-modeled quality compression; quality codec in CRAM 3.1. Integrated via [htscodecs](https://github.com/samtools/htscodecs). BSD.
  - Bonfield, J.K. and Mahoney, M.V. (2013). Compression of FASTQ and SAM format sequencing data. *PLoS ONE*, 8(3):e59190.
- [ENANO](https://github.com/guilledufort/EnanoFASTQ) by Guillermo Dufort y Álvarez et al. — sequence-aware quality modeling with stability tracking.
  - Dufort y Álvarez, G. et al. (2020). ENANO: encoder for nanopore FASTQ files. *Bioinformatics*, 36(16):4506–4507.

---

# BZ

Lossless BAM compressor using columnar stream separation, alignment-aware consensus-delta sequence encoding, and context-adaptive quality coding.

## Overview

BZ extends QZ's columnar-decomposition + BSC/BWT approach to BAM files. BAM uses generic BGZF (blocked gzip) compression, which does not exploit the structure of aligned sequencing data. CRAM improves on this with reference-based encoding, but requires an external reference genome. BZ achieves high compression without any external reference by exploiting the coordinate-sorted structure of BAM files: overlapping reads at the same genomic position share most of their bases, so a locally-built consensus reduces sequence data to a sparse diff stream.

**Key innovation: consensus-delta encoding.** Since BAM is coordinate-sorted and reads overlap, BZ builds a local consensus from the reads themselves (no external reference), then XOR-encodes each read's bases against the consensus. The resulting diff stream is >99% zeros. Under BWT, this compresses to ~0.1 bits/base — roughly 60x better than raw sequence data through BSC alone.

## Results

NA12878 chromosome 20, Illumina low-coverage WGS (3,245,545 aligned reads, 100 bp). Byte-identical roundtrip verified via MD5 on BGZF-decompressed BAM content.

### Compression ratio

| Format | Size | Ratio vs BAM |
|--------|------|-------------|
| Uncompressed BAM content | 1,133 MB | — |
| **BAM** (BGZF/gzip) | **298 MB** | **1.00x** |
| **BZ** | **135 MB** | **2.20x** |

### Per-stream breakdown (chunk 0, 2.5M records)

| Stream | Raw | Compressed | Ratio | Notes |
|--------|-----|-----------|-------|-------|
| consensus | 24.3 MB | 10.1 MB | 2.4x | Local consensus, packed nibbles |
| ref_id | 10.0 MB | 364 B | 27,473x | Delta-encoded, single chromosome |
| pos | 10.0 MB | 1.8 MB | 5.5x | Delta-encoded positions |
| mapq | 2.5 MB | 242 KB | 10.3x | Mapping quality |
| bin | 5.0 MB | 20 KB | 244x | BAM bin values |
| flag | 5.0 MB | 770 KB | 6.5x | SAM flags |
| next_ref_id | 10.0 MB | 29 KB | 348x | Delta vs ref_id (mate pairs) |
| next_pos | 10.0 MB | 2.7 MB | 3.7x | Delta vs pos (mate pairs) |
| tlen | 10.0 MB | 2.5 MB | 4.0x | Zigzag-encoded template length |
| read_name | 50.0 MB | 2.2 MB | 23.1x | Varint-framed read names |
| cigar | 15.2 MB | 805 KB | 18.9x | Varint-framed CIGAR ops |
| **seq_diff** | **119.4 MB** | **1.9 MB** | **62.5x** | **XOR vs consensus (>99% zeros)** |
| seq_extra | 13.1 MB | 6.0 MB | 2.2x | Insertion/softclip bases |
| **qual** | **255.0 MB** | **73.3 MB** | **3.5x** | **quality_ctx range coder** |
| aux | 389.9 MB | 6.2 MB | 63.2x | Auxiliary tags |

Quality scores dominate the compressed output (54% of archive size), followed by consensus data (7.4%) and seq_extra (4.4%). The consensus-delta encoding reduces the 119 MB aligned sequence stream to just 1.9 MB.

### Speed

| Operation | Time | Throughput |
|-----------|------|-----------|
| Compress | 22.4 s | 13.3 MB/s (from BAM) |
| Decompress | 19.6 s | 6.9 MB/s (to BAM) |

### Benchmark hardware

| Component | Specification |
|-----------|--------------|
| CPU | Intel Core i5-12400T (6 cores / 12 threads, Alder Lake) |
| RAM | 16 GB DDR4 |
| OS | Linux 6.6.87 (WSL2) |
| Rust | nightly, release profile |

## Architecture

```
BAM input (BGZF-compressed, coordinate-sorted)
    |  MultithreadedReader (parallel inflate)
    v
Raw BAM records
    |  Per-chunk (2.5M records):
    |    Pass 1: build local consensus (majority vote per position)
    |    Pass 2: delta-encode sequences against consensus
    v
15 columnar streams per chunk:
    +-- consensus    (local consensus, packed nibbles)
    +-- ref_id       (delta-encoded i32)  --+
    +-- pos          (delta-encoded i32)    | Fixed-width
    +-- mapq         (u8)                   |
    +-- bin          (u16 LE)               |
    +-- flag         (u16 LE)               |
    +-- next_ref_id  (delta vs ref_id)      |
    +-- next_pos     (delta vs pos)         |
    +-- tlen         (zigzag-encoded i32) --+
    +-- read_name    (varint + bytes)    --+
    +-- cigar        (varint + bytes)      | Variable-length
    +-- seq_diff     (XOR vs consensus)    |
    +-- seq_extra    (ins/softclip bases)  |
    +-- qual         (quality_ctx or BSC)  |
    +-- aux          (varint + bytes)    --+
    |
    |  rayon::join:
    |    BSC parallel (14 streams, 25 MB blocks, BWT + adaptive QLFC)
    |    quality_ctx parallel (500K-read blocks, context-adaptive range coder)
    v
BZ archive (header + N chunks, each with 15 compressed streams)
```

### Consensus-delta encoding

For each chunk of coordinate-sorted reads:

1. **Build consensus (pass 1):** For each reference position covered by any read, count base occurrences (A, C, G, T, N). Consensus = majority vote.
2. **Delta-encode (pass 2):** For each read, walk CIGAR ops:
   - **M/=/X** (alignment): XOR read nibble with consensus nibble -> diff stream (0 = match)
   - **I/S** (insertion/softclip): raw nibbles -> extra stream
   - **D/N/H** (deletion/skip/hardclip): no read bases consumed
   - **No CIGAR** (unmapped): all bases -> extra stream
3. **Reconstruction:** Walk CIGAR, XOR diff nibbles with consensus to recover originals, read insertion bases from extra stream.

Since most aligned bases match the consensus, the diff stream is >99% zeros. BWT clusters these zeros into a single run, achieving ~0.1 bits/base on the aligned sequence data.

### Quality compression

BZ reuses QZ's `quality_ctx` context-adaptive range coder for quality scores. The coder conditions on (position bin, previous quality, stability flag, previous base, current base) — 160,000 adaptive contexts that exploit Illumina's positional quality degradation and sequence-dependent error profiles. Quality data is partitioned into 500K-read blocks and compressed in parallel via `rayon::par_iter`.

When quality data contains unavailable scores (0xFF bytes), BZ falls back to BSC compression for the quality stream.

### Parallelism

- **BGZF I/O:** Multi-threaded BGZF reader (compression) and writer (decompression) via noodles, up to 4 worker threads for parallel inflate/deflate.
- **BSC compression:** 14 non-quality streams compressed in 4 parallel groups via nested `rayon::join`. Each stream uses 25 MB parallel blocks.
- **quality_ctx:** Runs in parallel with BSC via `rayon::join`. Quality data split into 500K-read blocks compressed via `rayon::par_iter`.
- **Decompression:** Parallel BSC decompression in groups, parallel sequence reconstruction across records, multi-threaded BGZF output.

### BZ file format

The BZ archive is a chunk-based binary format. All integers are little-endian. Each chunk contains 15 independently compressed columnar streams from up to 2.5M BAM records.

#### File layout

```
┌─────────────────────────────────────────────────────────────────┐
│  Global Header                                                   │
│    ├── Fixed fields (17 bytes)                                   │
│    └── SAM header payload (variable, BSC-compressed)             │
│  Chunk 0                                                         │
│    ├── Chunk header (64 bytes)                                   │
│    └── 15 compressed streams                                     │
│  Chunk 1                                                         │
│    ├── Chunk header (64 bytes)                                   │
│    └── 15 compressed streams                                     │
│  ...                                                             │
└─────────────────────────────────────────────────────────────────┘
```

#### Global header

```
Offset  Size  Type    Field                  Description
------  ----  ------  ---------------------  -----------
0       2     u8[2]   magic                  "BZ" (0x42 0x5A)
2       1     u8      version                2
3       1     u8      reserved               0
4       1     u8      flags                  Bit 0: consensus-delta encoding
                                             Bit 1: quality_ctx compression
5       8     u64     num_records            Total records across all chunks
13      4     u32     num_chunks             Number of chunks
17      4     u32     sam_header_len         Compressed SAM header size (bytes)
21      N     bytes   sam_header_data        BSC-compressed payload:
                                               [header_raw_len: u32]
                                               [header_raw: SAM header text]
                                               [ref_dict_bytes: BAM ref dict]
```

The SAM header payload contains both the raw SAM header text and the BAM reference dictionary (n_ref count + per-reference name/length entries), compressed together as a single BSC block.

#### Chunk header (64 bytes)

Each chunk begins with a fixed-size header:

```
Offset  Size  Type      Field            Description
------  ----  --------  ---------------  -----------
0       4     u32       num_records      Records in this chunk
4       60    u32[15]   stream_sizes     Compressed size of each stream (bytes)
```

The 15 stream size entries correspond to streams 0–14. The compressed stream data follows immediately, concatenated in order.

#### Stream index

| Index | Stream | Encoding | Description |
|-------|--------|----------|-------------|
| 0 | consensus | Packed nibbles | Local consensus built from overlapping reads |
| 1 | ref_id | Delta-encoded i32 | Reference sequence ID |
| 2 | pos | Delta-encoded i32 | Alignment position |
| 3 | mapq | Raw u8 | Mapping quality |
| 4 | bin | Raw u16 LE | BAM bin value |
| 5 | flag | Raw u16 LE | SAM flags |
| 6 | next_ref_id | Delta vs ref_id | Mate reference ID |
| 7 | next_pos | Delta vs pos | Mate position |
| 8 | tlen | Zigzag-encoded i32 | Template length |
| 9 | read_name | Varint + bytes | Read name (NUL-terminated) |
| 10 | cigar | Varint + bytes | Raw CIGAR op bytes |
| 11 | seq_diff | Packed nibbles | XOR of read bases vs consensus |
| 12 | seq_extra | Packed nibbles | Insertion/softclip bases |
| 13 | qual | quality_ctx or BSC | Quality scores |
| 14 | aux | Varint + bytes | Auxiliary tags |

All streams except quality (index 13) are compressed with BSC (BWT + adaptive QLFC, 25 MB blocks). Quality uses `quality_ctx` context-adaptive range coding when `flags & 0x02` is set, otherwise BSC.

#### Consensus stream format (stream 0, after BSC decompression)

```
+0      4     u32     num_segments     Number of reference segments
Per segment:
  +0    4     i32     ref_id           Reference sequence ID
  +4    4     i32     start_pos        Start position on reference
  +8    4     u32     length           Number of consensus bases
  +12   N     bytes   bases            Packed nibbles (2 bases per byte,
                                       high nibble first), ceil(length/2) bytes
```

#### Variable-length stream format (streams 9, 10, 14, after BSC decompression)

```
Per record:
  varint    field_len     Unsigned LEB128 length
  bytes     field_data    Raw field bytes
```

#### Quality stream format (stream 13)

When using `quality_ctx` (`flags & 0x02`):

```
+0      4     u32     num_blocks       Number of quality_ctx sub-blocks
Per block:
  +0    4     u32     block_len        Compressed block size
  +4    N     bytes   block_data       quality_ctx compressed payload
```

Each block contains up to 500K reads. When falling back to BSC (no `flags & 0x02`), the quality stream is a single BSC-compressed blob of varint-framed quality bytes.

#### Fixed-width stream formats (streams 1–8, after BSC decompression)

Streams 1–8 contain one value per record, concatenated with no framing:

| Stream | Width | Encoding |
|--------|-------|----------|
| ref_id (1) | 4 bytes | i32 LE, delta-encoded (value = prev + stored) |
| pos (2) | 4 bytes | i32 LE, delta-encoded |
| mapq (3) | 1 byte | Raw u8 |
| bin (4) | 2 bytes | Raw u16 LE |
| flag (5) | 2 bytes | Raw u16 LE |
| next_ref_id (6) | 4 bytes | i32 LE, delta vs ref_id (value = ref_id + stored) |
| next_pos (7) | 4 bytes | i32 LE, delta vs pos (value = pos + stored) |
| tlen (8) | 4 bytes | u32 LE, zigzag-encoded (decode: (n >> 1) ^ -(n & 1)) |

## Usage

### Compress

```bash
bz compress -i aligned.bam -o aligned.bz
bz compress -i aligned.bam -o aligned.bz -t 8    # limit to 8 threads
```

### Decompress

```bash
bz decompress -i aligned.bz -o restored.bam
```

### Verify roundtrip

```bash
bz compress -i input.bam -o compressed.bz
bz decompress -i compressed.bz -o roundtrip.bam
# Compare BGZF-decompressed content:
diff <(samtools view input.bam) <(samtools view roundtrip.bam)
```

### CLI reference

**Compress:**

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input FILE` | Input BAM file | required |
| `-o, --output FILE` | Output BZ archive | required |
| `-w, --working-dir PATH` | Working directory | `.` |
| `-t, --threads N` | Thread count (0 = auto) | auto |

**Decompress:**

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input FILE` | Input BZ archive | required |
| `-o, --output FILE` | Output BAM file | required |
| `-w, --working-dir PATH` | Working directory | `.` |
| `-t, --threads N` | Thread count (0 = auto) | auto |

## License

See LICENSE file.
