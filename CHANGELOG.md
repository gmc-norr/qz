# Changelog

## Unreleased

### Breaking Changes

- **Archive format v2**: Archives now start with an 8-byte prefix (`QZ\x02\x00` + header_size u32 LE). Archives produced by previous versions are no longer readable.
- **FastqRecord uses `Vec<u8>`**: `id`, `sequence`, and `quality` fields are now byte vectors instead of `String`.

### Added

- **`bz extract` — BAM to QZ extraction** — New `bz extract` subcommand reads a coordinate-sorted BAM file, pairs reads by name (R1/R2), and compresses to paired QZ archives (`{prefix}_R1.qz` and `{prefix}_R2.qz`). Reads are in matched pair order. Secondary and supplementary alignments are skipped. Coordinate-sorted read order provides good BSC compression locality.
- **CUDA GPU acceleration (opt-in)** — Build with `--features cuda` to enable GPU-accelerated BWT via libcubwt. Falls back to CPU gracefully when GPU is unavailable.
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
