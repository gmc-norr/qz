//! NUMA-aware multi-process decode/compress orchestration (CLI-only — never in
//! qz-lib).
//!
//! The format-agnostic half (topology, pinning, runtime, the decision gate, and
//! part-file assembly) lives in the shared `numa-core` crate and is re-exported
//! here so the existing `crate::numa::{topology,gate,pin,runtime,assemble}` paths
//! keep resolving. The qz-format-specific drivers (decode `driver`, `compress_driver`,
//! `compress_gate`, FASTQ `splitter`) stay here.
//!
//! See docs/superpowers/specs/2026-06-15-numa-decode-design.md.

pub use numa_core::{
    assemble, gate, parse_numa_mode, pin, runtime, topology, NumaMode, MAX_DIRECT_OUTPUT_BYTES,
};

pub mod compress_driver;
pub mod compress_gate;
pub mod driver;
pub mod splitter;
