# Changelog

## Unreleased

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

### Known limitations (not fixed in this pass)

- **Casava header columnar codec uses `u16` for tile/x/y** — NovaSeq X x/y coordinates can exceed 65535. Such headers correctly fall back to raw BSC (no data loss) and now emit a one-line `warn!` so operators can see the slower path. A v0x03 format with `u32` fields would be needed to compress these as efficiently as standard Illumina headers.
- **Bench-only OpenZL graph variants ship in qz-lib** — `openzl::compress_ace/dna_numeric/clustering/...` are called only from `bench_openzl_graphs`. Left in qz-lib for now (no correctness impact; small binary-size cost).

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
