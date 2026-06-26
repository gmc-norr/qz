//! Reference-direct paired-end sequence encoding (reference-free decode).
//!
//! Opt-in `--reference` mode (spec §4): a single streaming pass maps each read
//! against a reference FASTA and encodes it as edits over the covered reference
//! slices (the "backing"); off-reference reads spill to a per-(chunk,mate)
//! fallback literal pool written inline into that chunk's payload. The reference
//! is never stored and is never required to decode — `decode_reference_file_to_sinks`
//! reconstructs every read from the backing + per-(chunk,mate) pools alone.
//!
//! This module is the root: it declares the reference submodules and the two
//! halves of the codec — [`encode`] (compress) and [`decode`] (decompress +
//! verify) — and re-exports their entry points. The split mirrors single-end's
//! `compress_impl` / `decompress_impl`.

mod backing; // Task 5
mod coverage; // Task 8
mod edits; // Task 8
mod fallback; // Task 11
mod format; // Task 3
mod mapping; // Task 10
mod positions; // Task 10 (crate-internal: no consumer outside `reference`)
mod refmeta; // Task 4 (crate-internal: no consumer outside `reference`)

mod decode;
mod encode;
mod merge;
mod stream;
mod stream_single;

// Codec entry points, re-exported at their original visibility (preserving the
// `qz_lib::compression::reference::{compress,decompress}_reference` public API).
pub use decode::{
    decompress_reference, decompress_reference_interleaved, decompress_reference_single,
};
pub use encode::compress_reference;
pub use encode::compress_reference_single;
pub(crate) use encode::compress_reference_byte_range;
pub(crate) use decode::decode_reference_range;
pub(crate) use decode::decode_reference_range_single;
pub(crate) use decode::verify_reference;
pub(crate) use decode::verify_reference_single;
pub(crate) use encode::reject_unsupported;
// Reference encode pieces shared with the streaming entry point (`stream.rs`).
pub(crate) use encode::Placed;
pub(crate) use encode::{
    coverage_to_intervalmap, map_and_diff, map_and_diff_single, mark_cov,
    reference_prelude, reference_prelude_single, resolve_ref_chunk_records_cfg,
};
// `encode_chunk_single` / `encode_mate_single` are reached directly as
// `super::encode::…` by `stream_single.rs` (mirroring how `stream.rs` imports
// `encode_chunk`), so they intentionally are NOT re-exported here.
pub(crate) use format::{validate_reference_directory, validate_reference_directory_single};
// Reference-aware archive merge (NUMA-compress track) + the normalized-FASTA loader
// it needs to re-derive the coverage globals.
pub(crate) use mapping::load_reference_normalized;
pub use mapping::reference_index_sidecar_path;
pub use mapping::{build_reference_index, ReferenceIndexBuildStats};
pub use mapping::{ensure_reference_index, reference_index_status, ReferenceIndexStatus};
pub(crate) use merge::merge_reference_core;

#[cfg(test)]
mod reachability {
    #[test]
    fn paired_helpers_reachable() {
        // pub(crate) reach from the reference module — compile-time proof.
        let mut buf = Vec::new();
        crate::compression::paired::streams::write_block_stream(&mut buf, &[]);
    }
}
