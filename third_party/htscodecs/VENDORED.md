# Vendored: htscodecs (minimal subset)

This directory is a vendored, in-tree **minimal subset** of upstream htscodecs —
just the files qz needs to build the `fqzcomp` quality codec. It is compiled from
source by `crates/qz-lib/build.rs` and linked into qz via FFI. fqzcomp is qz's and
bz's default quality-score compressor.

- Upstream URL: https://github.com/samtools/htscodecs
- Exact SHA: `69185ced605aa1ffde9c78b3d17223a114d34ddb` (`69185ce`)
- Commit: "Fix arith-dynamic encoding of files stored using X_CAT", 2026-02-03
- License: BSD-3-Clause (c) Genome Research Limited, except `c_range_coder.h`
  which is public domain (derived from work by Eugene Shelwien). See `LICENSE.md`,
  copied verbatim, for the per-file exceptions.
- Vendored on: 2026-06-25

Previously this was a gitignored upstream clone that every contributor had to fetch
separately (`git clone … && git checkout 69185ce`) before qz would build. It is now
tracked in qz's history, so a fresh checkout builds with no extra steps and the exact
reviewed source is locked (no re-clone drift). This matches the already-vendored
`third_party/libbsc` and `crates/strobealign`.

## What was copied (minimal subset)

`crates/qz-lib/build.rs` compiles only **two** translation units —
`htscodecs/fqzcomp_qual.c` and `htscodecs/utils.c` — with the include path set to
`htscodecs/` alone. The vendored subset is therefore exactly those two `.c` files
plus the transitive closure of their **local** (`#include "…"`) headers; every other
include they use is a system header. That closure is these 9 files under
`htscodecs/htscodecs/`:

- `fqzcomp_qual.c`, `fqzcomp_qual.h` — the fqzcomp quality codec
- `utils.c`, `utils.h` — htscodecs TLS-alloc helpers (see NO_THREADS below)
- `varint.h`, `varint2.h` — varint encode/decode
- `c_range_coder.h`, `c_simple_model.h` — the range coder + adaptive symbol model
- `config.h` — see "Local customizations" below

Plus top-level `LICENSE.md` and `README.md` (provenance/attribution). Everything else
in upstream htscodecs — the rANS SIMD variants, the name tokeniser, the arith codec,
tests, build scaffolding — is **not** compiled by qz and is intentionally excluded.

## Local customizations (no source patch)

Unlike `third_party/libbsc` (which carries `QZ-SECURITY-PATCH.diff`), no upstream
htscodecs `.c`/`.h` source is patched. Two build-level customizations matter:

1. **`config.h` is a hand-authored 10-line minimal stub**, not the autoconf-generated
   header — it defines `HAVE_BUILTIN_PREFETCH` and leaves `HAVE_PTHREAD` undefined, so
   the subset builds with no `configure` step.
2. **`-DNO_THREADS`** is passed by `build.rs`. It makes htscodecs' `htscodecs_tls_alloc`
   / `htscodecs_tls_free` in `utils.c` fall back to plain `calloc(1, size)` / `free`,
   disabling the per-thread TLS arena that never returns the ~25 MB fqzcomp qual model
   to the OS (the bounded-RAM fqz decode leak — see the long comment in `build.rs`).
   `build.rs` has a **fail-closed sentinel guard** that greps `utils.c` for the
   `NO_THREADS` token *and* the exact `calloc(1, size)` / `free(ptr)` fallback bodies,
   so a re-vendor that changes those semantics trips the build instead of silently
   reintroducing the per-thread leak.

## Re-syncing upstream

Re-fetch upstream at the desired SHA, copy the 9 files above (re-run the
`#include "…"` closure in case it changed) plus `LICENSE.md`/`README.md`, re-author
`config.h` if upstream's autoconf inputs changed, drop any nested `.git`, update the
SHA/date above, and confirm the `build.rs` NO_THREADS sentinel guard still passes.
Re-verify reference/single-end decode peak RSS stays constant in read count.
