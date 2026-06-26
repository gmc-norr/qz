# Third-party licenses

qz and bz are licensed under [Apache-2.0](LICENSE). They incorporate third-party
components, all under permissive licenses (Apache-2.0, BSD, MIT, ISC, Zlib,
Unicode-3.0, 0BSD, Unlicense, and public domain). **No component is under a
copyleft license (GPL/LGPL/MPL/etc.).** A binary distribution of qz/bz should
include the notices below.

This file covers two groups: (1) vendored native components whose source ships in
this repository and is statically linked into the binaries, and (2) Rust crate
dependencies fetched from crates.io.

## 1. Vendored components (in-tree, statically linked)

Each component's full license text is retained verbatim at the path shown.

| Component | License | Copyright | License text | Builds |
| --------- | ------- | --------- | ------------ | ------ |
| **libbsc** | Apache-2.0 | Ilya Grebnov | `third_party/libbsc/LICENSE` | always |
| **libsais** | Apache-2.0 | Ilya Grebnov | `third_party/libbsc/libbsc/bwt/libsais/LICENSE` | always |
| **htscodecs** (minimal subset) | BSD-3-Clause | Genome Research Limited | `third_party/htscodecs/LICENSE.md` | always |
| ↳ `c_range_coder.h` | Public Domain | Eugene Shelwien | (noted in `third_party/htscodecs/LICENSE.md`) | always |
| **strobealign** | MIT | Kristoffer Sahlin | `crates/strobealign/LICENSE` | always |
| **libcubwt** | Apache-2.0 | Ilya Grebnov | `third_party/libbsc/libbsc/bwt/libcubwt/LICENSE` | `--features cuda` only |

Notes:
- The vendored htscodecs subset is BSD-3-Clause **except** `c_range_coder.h`, which
  its LICENSE.md lists as public domain (derived from work by Eugene Shelwien). The
  other public-domain files that LICENSE.md lists (`rANS_byte.h`, `rANS_word.h`,
  CC0) are not part of the vendored subset. See `third_party/htscodecs/VENDORED.md`.
- BSD-3-Clause clause 3 (no-endorsement): the htscodecs authors' names and "Genome
  Research Limited" are not used to endorse or promote qz/bz.
- libcubwt is compiled only in CUDA builds (`--features cuda`); the default CPU
  build does not link it.
- None of the Apache-2.0 components ship a `NOTICE` file, so there is no additional
  NOTICE text to propagate.

## 2. Rust dependencies (crates.io)

The Rust dependency tree (205 external crates, all features enabled) is entirely
permissive. License distribution by SPDX expression:

| Count | License expression |
| ----- | ------------------ |
| 109 | MIT OR Apache-2.0 |
| 32 | MIT |
| 16 | MIT/Apache-2.0 |
| 14 | Apache-2.0 WITH LLVM-exception OR Apache-2.0 OR MIT |
| 14 | Apache-2.0 OR MIT |
| 5 | Apache-2.0 |
| 1 | Apache-2.0 WITH LLVM-exception |
| 4 | Unlicense OR MIT |
| 3 | Zlib |
| 2 | BSD-2-Clause OR Apache-2.0 OR MIT |
| 1 each | 0BSD OR MIT OR Apache-2.0 · Zlib OR Apache-2.0 OR MIT · MIT OR Zlib OR Apache-2.0 · (MIT OR Apache-2.0) AND Unicode-3.0 · MIT OR Apache-2.0 OR LGPL-2.1-or-later |

Two expressions warrant a note, neither creating a copyleft obligation:
- **`r-efi`** offers `MIT OR Apache-2.0 OR LGPL-2.1-or-later`. Under the `OR`, qz/bz
  use it under MIT/Apache-2.0; the LGPL option is not exercised. (`r-efi` is a
  UEFI-target crate and is not linked into the Linux qz/bz binaries.)
- **`unicode-ident`** is `(MIT OR Apache-2.0) AND Unicode-3.0`; Unicode-3.0 is a
  permissive license covering bundled Unicode data tables.

For a binary release, regenerate the full per-crate license manifest (with license
texts) using a license tool, e.g.:

```bash
cargo install cargo-about && cargo about generate about.hbs   # full texts
# or, for a quick SPDX listing:
cargo install cargo-license && cargo license
```

The summary above was produced from `cargo metadata` on this revision; re-run after
any dependency change.
