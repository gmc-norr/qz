# Vendored: libbsc

This directory is a vendored, in-tree copy of upstream libbsc (the BWT + LZP +
adaptive-QLFC C/C++ codec qz uses as its default sequence/header compressor). It
is compiled from source by `crates/qz-lib/build.rs` and linked into qz via FFI.

- Upstream URL: https://github.com/IlyaGrebnov/libbsc
- Branch: `master`
- Exact SHA: `baffa62c70b6ebbecc9af14ce550e965ea247680`
- Version: `v3.3.12-4-gbaffa62` (VERSION file: 3.3.12)
- License: Apache-2.0 (c) Ilya Grebnov (see `LICENSE`, copied verbatim;
  bundled libsais under `libbsc/bwt/libsais/` carries its own Apache-2.0 `LICENSE`)
- Vendored on: 2026-06-24

Previously this was a gitignored upstream clone, symlinked into each worktree,
with the security patch below living only as uncommitted working-tree edits on
one machine. It is now tracked in qz's history so the patch can never be lost on
a re-clone. (`third_party/htscodecs` is now also vendored in-tree, as a minimal
subset — see `third_party/htscodecs/VENDORED.md`.)

## What was copied

The whole upstream working tree, minus its nested `.git` and the `/build` output
dir. That includes files qz does **not** compile (`bsc.cpp` CLI driver,
`CMakeLists.txt`, `libbsc/bwt/libcubwt/libcubwt.cu` GPU path, `st/`, docs),
retained verbatim to minimize divergence and keep `#include`s self-contained.

`crates/qz-lib/build.rs` compiles only: `libbsc/bwt/libsais/libsais.c` plus the
C++ TUs `libbsc/libbsc.cpp`, `adler32/adler32.cpp`, `bwt/bwt.cpp`,
`coder/coder.cpp`, `coder/qlfc/qlfc.cpp`, `coder/qlfc/qlfc_model.cpp`,
`filters/detectors.cpp`, `filters/preprocessing.cpp`, `lzp/lzp.cpp`. Threading
uses `LIBBSC_FEATURE_FASTMODE` only (never `MULTITHREADING` — rayon owns
parallelism).

## Local patch applied: QZ-SECURITY-PATCH (bsc-decode-output-bound)

A heap-overflow hardening of the decode path, captured verbatim in
`QZ-SECURITY-PATCH.diff` (the full diff vs the upstream SHA above) so the delta
stays auditable and re-syncable. It threads the caller's total output-buffer
**capacity** through the decode chain so a malicious per-block size table in a
crafted archive cannot drive a write past the allocation. Files modified:

- `libbsc/coder/coder.cpp`, `libbsc/coder/coder.h` — `bsc_coder_decode_block` /
  `bsc_coder_decompress` take an `outputSize` capacity argument.
- `libbsc/coder/qlfc/qlfc.cpp`, `libbsc/coder/qlfc/qlfc.h` — the static /
  adaptive / fast QLFC decoders bound their writes to `outputSize`.
- `libbsc/lzp/lzp.cpp`, `libbsc/lzp/lzp.h` — LZP decode is output-bounded.
- `libbsc/libbsc/libbsc.cpp` — passes the real output capacity down.

`crates/qz-lib/build.rs` has a **fail-closed guard**: it greps the compiled
sources for the sentinel `QZ-SECURITY-PATCH (bsc-decode-output-bound)` and fails
the build if it is absent, so a re-vendor that silently drops the patch cannot
produce a vulnerable binary.

## Re-syncing upstream

Re-clone upstream at the desired SHA, re-apply `QZ-SECURITY-PATCH.diff` (or
re-derive the equivalent bound on the decode chain), drop the nested `.git`,
update the SHA/Version/date above, and confirm the build.rs sentinel guard still
passes.
