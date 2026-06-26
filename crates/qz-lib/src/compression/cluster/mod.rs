//! Order-drop (cluster-reorder) compression mode — archive_type = 3. See
//! docs/superpowers/specs/2026-06-21-order-drop-cluster-mode-design.md.
pub(crate) mod checksum;
pub(crate) mod decode;
pub(crate) mod encode;
pub(crate) mod sort;
pub(crate) mod zstd_codec;

pub(crate) use decode::decompress_cluster;
pub(crate) use decode::verify_cluster;
pub(crate) use encode::{compress_cluster, compress_cluster_paired, reject_unsupported};
