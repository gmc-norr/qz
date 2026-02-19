# Changelog

## Unreleased

### Breaking Changes

- **Archive format v2**: Archives now start with an 8-byte prefix (`QZ\x02\x00` + header_size u32 LE). Archives produced by previous versions are no longer readable.
- **FastqRecord uses `Vec<u8>`**: `id`, `sequence`, and `quality` fields are now byte vectors instead of `String`.

### Added

- **`qz verify` — archive integrity verification** — Fully decompresses all streams without writing to disk, verifying that the archive is intact and decompressible. Reports archive metadata (num reads, encoding type, compressor info, stream sizes) and a CRC32 hash of the reconstructed FASTQ output. The CRC32 and byte count are consistent across compression modes (default, ultra) for the same input. Usage: `qz verify -i archive.qz`.
- **Columnar header compression (new default)** — Illumina headers are parsed into typed binary columns (read_num, instrument combo, lane, tile, x, y) and each column is BSC-compressed independently in parallel. Falls back to raw BSC for non-Illumina formats. Saves ~8.5% on headers (e.g., 84.1 MB vs 91.9 MB on 10M reads), improving overall ratio from 7.99x to 8.13x. Selectable via `--config '{"header_compressor": 3}'` (Columnar) or `--config '{"header_compressor": 1}'` (BSC).
- **`bz extract` — BAM to QZ extraction** — New `bz extract` subcommand reads a coordinate-sorted BAM file, pairs reads by name (R1/R2), and compresses to paired QZ archives (`{prefix}_R1.qz` and `{prefix}_R2.qz`). Reads are in matched pair order. Secondary and supplementary alignments are skipped. Coordinate-sorted read order provides good BSC compression locality.
- **`--config` JSON flag** — Pass a JSON file to override `AdvancedOptions` at runtime (`qz compress --config opts.json ...`). Only fields present in the JSON are overridden; missing fields use defaults. Enables automated parameter tuning without recompilation.
- **Compression parameter optimizer** — Python-based GA optimizer (`scripts/optimizer/`) searches QZ's 17 tunable parameters (compressor selection, encoding modes, block sizes) to maximize compression ratio. Includes sensitivity analysis, checkpoint/resume, LLM-guided mutations via Ollama, and parallel evaluation. Run with `python3 -m scripts.optimizer.main --qz-bin target/release/qz --input data.fastq`.
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

### Changed

- **`CompressConfig` simplified** — Production-only fields remain in `cli.rs`; experimental flags live in `AdvancedConfig` or benchmark code.
- **Ultra level table tuned** — BSC block sizes now scale with ultra level (25→100 MB) for progressively better BWT context. Level 2 chunk size increased from 2M to 2.5M reads. Level 4 reduced from 2 to 1 parallel chunk for better memory efficiency. Each level now specifies `bsc_block_mb` and `max_seq_block_mb` instead of using hardcoded values.
- **`--local-reorder` uses ReorderLocal** — Previously used Delta encoding (encoding_type=8) which was consistently worse than ReorderLocal with identity reorder. Now maps to ReorderLocal with auto-selected ultra level, matching `--ultra 0` behavior. Existing Delta-encoded archives (encoding_type=8) remain decompressible.

### Fixed

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
