# Vendored: strobealign

This crate is a vendored, library-only copy of upstream strobealign, owned and
built in-tree as a member of the qz workspace.

- Upstream URL: https://github.com/ksahlin/strobealign
- Branch: `main`
- Exact SHA: `213b379e8f7748d3c869e1b60374f1cdea52409d`
- Version: v0.18.0-alpha
- License: MIT (c) Kristoffer Sahlin (see `LICENSE`, copied verbatim)
- Vendored on: 2026-06-07

## What was copied

Only `src/`, `ext/`, `build.rs`, `Cargo.toml`, and `LICENSE` from the upstream
tarball. NOT copied: `cpp/`, `python/`, `pyproject.toml`, `setup.py`,
`CMakeLists.txt`, `tests/`, `Cargo.lock`, `.github/`, and the various top-level
docs.

`ext/ssw/` contains the SSW C sources that `build.rs` compiles
(`cc::Build::new().file("ext/ssw/ssw.c").compile("ssw")`).

## Local patches applied

To make this a library-only workspace member (not its own crate root):

1. Deleted `src/main.rs` (the binary entry point). The `#[global_allocator]
   static GLOBAL: MiMalloc` lived only in `src/main.rs`, so there is NO global
   allocator set by this crate now — no clash with qz's own mimalloc allocator.
   (`src/logger.rs` was main-only and is not declared in `lib.rs`; it is dead but
   harmless. No `src/bin/` existed.)
2. Added `autobins = false` to `[package]` (lib-only; no auto-discovered bins).
3. Removed the `[workspace]` table and its `python` member — this crate is now a
   member of qz's root workspace, not a workspace root itself.
4. Removed the crate-local `[profile.dev]` / `[profile.release]` /
   `[profile.profiling]` tables (non-root profile tables are ignored with a
   warning; profiles come from qz's root workspace).
5. Added `#![allow(warnings, clippy::all)]` as the very first line of
   `src/lib.rs` so upstream lints don't pollute qz's
   `cargo clippy --workspace --all-targets`.

The `mimalloc` dependency is kept in `Cargo.toml` (it is now unused by the lib,
but harmless; left in place to minimize divergence from upstream). `clap`,
`sigpipe`, and the `dev-dependencies` (assert_cmd/predicates/temp-dir/temp-file)
are likewise retained even though their only users (the binary / upstream
`tests/`) were not vendored.

## Determinism patch (Task 10)

Added to `src/maponly.rs` (same module as `get_nam_pairs` +
`get_best_paired_mapping_location`, so no visibility changes were needed):

- `pub fn map_pair_deterministic(r1_seq, r2_seq, index, chainer, rescue_distance,
  mcs_strategy, mu: f32, sigma: f32) -> (Option<Nam>, Option<Nam>)`
- `pub fn map_single_deterministic(seq, index, chainer, rescue_distance,
  mcs_strategy) -> Option<Nam>` — the per-read half of `map_pair_deterministic`
  (single-end reference mode): `get_nams_by_chaining` → `deterministic_sort` →
  best NAM, WITHOUT pairing and WITHOUT mu/sigma.
- `fn deterministic_sort(nams: &mut [Nam])` helper (shared by both)

This is a deterministic, RNG-free, immutable-insert-size variant of
`map_paired_end_read`'s NAM-only path, used by qz's reference mode where
reproducible mappings are required. It replicates the upstream control flow
(get_nams_by_chaining → get_nam_pairs → best-pair-vs-individual selection) with
exactly two behavioural changes, both removing non-determinism:

1. The two `sort_nams(&mut nams, rng)` calls are replaced by `deterministic_sort`,
   which sorts by `score` descending and then breaks score ties with a fixed
   total order on `(ref_id, ref_start, is_revcomp)` instead of the rng-based
   `shuffle_best` tie-shuffle.
2. The best-location selection inlines `get_best_paired_mapping_location`'s rule
   (prefer the best proper pair when its score `>= individual_score / 2.0`, else
   the best individual NAM per mate) but OMITS the
   `insert_size_distribution.update(...)` side effect — mu/sigma are passed in by
   value and never mutated.

No other upstream functions were modified; `get_nam_pairs`,
`get_best_paired_mapping_location`, `NamPair`, etc. are unchanged. The patch is a
single added `pub fn` + one private helper, so upstream re-sync is a clean append.
