//! NUMA-compress gate (spec §6): in-process vs sharded. Pure function, unit-tested
//! like the decode gate. Keyed on the EXACT uncompressed input size (no estimate).

use crate::numa::gate::{
    mem_filtered_nodes, require_shardable_topology, validate_fixed_workers, IoKind,
    MIN_WORKLOAD_BYTES,
};
use crate::numa::topology::{NumaNode, Topology};
use crate::numa::NumaMode;

#[derive(Clone, Debug)]
pub enum CompressPlan {
    InProcess,
    /// `worker_count` node-pinned workers; the splitter produces the ranges.
    Sharded { worker_count: u32, nodes: Vec<NumaNode> },
}

#[derive(Clone, Debug)]
pub struct CompressGateInputs {
    pub input_bytes: u64,
    pub total_threads: usize,
    pub io: IoKind,
    pub allowed_cpus: Vec<u32>,
    pub allowed_mem_nodes: Vec<u32>,
    pub worker_peak_rss_bytes: u64,
    pub workdir_free_bytes: u64,
    pub output_fs_free_bytes: u64,
    pub workdir_output_same_fs: bool,
}

/// Measured qz FASTQ ratio (~8×) → part-archive + merge-temp sizing.
const COMPRESS_RATIO: u64 = 8;

pub fn decide_compress(
    mode: NumaMode,
    topo: Option<&Topology>,
    i: &CompressGateInputs,
) -> anyhow::Result<CompressPlan> {
    if mode == NumaMode::Off || mode == NumaMode::Fixed(1) {
        return Ok(CompressPlan::InProcess);
    }
    let fixed_n = match mode { NumaMode::Fixed(n) => Some(n), _ => None };
    let bail_or_fallback = |msg: &str| -> anyhow::Result<CompressPlan> {
        if fixed_n.is_some() { anyhow::bail!("--numa: {msg}") } else { Ok(CompressPlan::InProcess) }
    };

    let topo = match require_shardable_topology(topo, i.io) {
        Ok(t) => t,
        Err(msg) => return bail_or_fallback(&msg),
    };
    let usable = match mem_filtered_nodes(topo, &i.allowed_cpus, &i.allowed_mem_nodes, i.worker_peak_rss_bytes) {
        Ok(u) => u,
        Err(msg) => return bail_or_fallback(&msg),
    };

    // Capacity: parts (~input/8) land in working_dir; the merge temp (~input/8) lands
    // by the output. Same-fs → 2×; cross-fs → 1× on each.
    let part_bytes = i.input_bytes / COMPRESS_RATIO;
    if i.workdir_output_same_fs {
        if i.output_fs_free_bytes < part_bytes.saturating_mul(2) {
            return bail_or_fallback("filesystem lacks ~2x part size (parts + merge temp share one fs)");
        }
    } else {
        if i.workdir_free_bytes < part_bytes {
            return bail_or_fallback("working dir lacks ~1x part size for shard archives");
        }
        if i.output_fs_free_bytes < part_bytes {
            return bail_or_fallback("output filesystem lacks ~1x part size for the merge temp");
        }
    }

    let cap = (usable.len() as u32).min(i.total_threads.max(1) as u32);
    let worker_count = match fixed_n {
        Some(n) => {
            validate_fixed_workers(n, usable.len(), i.total_threads)?;
            n
        }
        None => {
            if i.input_bytes < MIN_WORKLOAD_BYTES { return Ok(CompressPlan::InProcess); }
            cap
        }
    };
    if worker_count < 2 {
        return Ok(CompressPlan::InProcess);
    }
    let nodes: Vec<NumaNode> = (0..worker_count as usize)
        .map(|idx| usable[idx % usable.len()].clone())
        .collect();
    Ok(CompressPlan::Sharded { worker_count, nodes })
}

/// Per-worker thread budget (parent uses `assign_node_threads` per range when
/// building the worker command). Exposed for the driver.
pub fn per_worker_threads(total_threads: usize, worker_count: u32) -> usize {
    (total_threads / worker_count as usize).max(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn topo2() -> Topology {
        Topology {
            nodes: vec![
                NumaNode { id: 0, cpus: (0..18).collect(), mem_available_bytes: 64 << 30 },
                NumaNode { id: 1, cpus: (18..36).collect(), mem_available_bytes: 64 << 30 },
            ],
            has_real_asymmetry: true,
        }
    }
    fn base() -> CompressGateInputs {
        CompressGateInputs {
            input_bytes: 4 << 30,
            total_threads: 36,
            io: IoKind::Files,
            allowed_cpus: (0..36).collect(),
            allowed_mem_nodes: vec![0, 1],
            worker_peak_rss_bytes: 2 << 30,
            workdir_free_bytes: 100 << 30,
            output_fs_free_bytes: 100 << 30,
            workdir_output_same_fs: false,
        }
    }

    #[test]
    fn off_and_fixed1_in_process() {
        assert!(matches!(decide_compress(NumaMode::Off, Some(&topo2()), &base()).unwrap(), CompressPlan::InProcess));
        assert!(matches!(decide_compress(NumaMode::Fixed(1), Some(&topo2()), &base()).unwrap(), CompressPlan::InProcess));
    }
    #[test]
    fn auto_shards_two_nodes() {
        match decide_compress(NumaMode::Auto, Some(&topo2()), &base()).unwrap() {
            CompressPlan::Sharded { worker_count, nodes } => { assert_eq!(worker_count, 2); assert_eq!(nodes.len(), 2); }
            _ => panic!("expected sharded"),
        }
    }
    #[test]
    fn no_asymmetry_fallback_auto_error_fixed() {
        assert!(matches!(decide_compress(NumaMode::Auto, None, &base()).unwrap(), CompressPlan::InProcess));
        assert!(decide_compress(NumaMode::Fixed(2), None, &base()).is_err());
    }
    #[test]
    fn stdio_fallback_auto_error_fixed() {
        let mut i = base(); i.io = IoKind::Stdio;
        assert!(matches!(decide_compress(NumaMode::Auto, Some(&topo2()), &i).unwrap(), CompressPlan::InProcess));
        assert!(decide_compress(NumaMode::Fixed(2), Some(&topo2()), &i).is_err());
    }
    #[test]
    fn subthreshold_workload_auto_only() {
        let mut i = base(); i.input_bytes = 1 << 20;
        assert!(matches!(decide_compress(NumaMode::Auto, Some(&topo2()), &i).unwrap(), CompressPlan::InProcess));
        assert!(matches!(decide_compress(NumaMode::Fixed(2), Some(&topo2()), &i).unwrap(), CompressPlan::Sharded { .. }));
    }
    #[test]
    fn same_fs_needs_2x_part_size() {
        let mut i = base(); i.workdir_output_same_fs = true;
        i.input_bytes = 16 << 30; // part ≈ 2 GiB → need ≈ 4 GiB
        i.output_fs_free_bytes = 3 << 30;
        assert!(matches!(decide_compress(NumaMode::Auto, Some(&topo2()), &i).unwrap(), CompressPlan::InProcess));
        i.output_fs_free_bytes = 5 << 30;
        assert!(matches!(decide_compress(NumaMode::Auto, Some(&topo2()), &i).unwrap(), CompressPlan::Sharded { .. }));
    }
    #[test]
    fn cross_fs_needs_1x_each() {
        let mut i = base(); i.input_bytes = 16 << 30; // part ≈ 2 GiB
        i.workdir_free_bytes = 100 << 30; i.output_fs_free_bytes = 1 << 30;
        assert!(matches!(decide_compress(NumaMode::Auto, Some(&topo2()), &i).unwrap(), CompressPlan::InProcess));
        i.output_fs_free_bytes = 100 << 30; i.workdir_free_bytes = 1 << 30;
        assert!(matches!(decide_compress(NumaMode::Auto, Some(&topo2()), &i).unwrap(), CompressPlan::InProcess));
        // Both filesystems have >1× each (part≈2 GiB, so 5 GiB each is plenty) →
        // cross-fs must SHARD (proves the path isn't stuck returning InProcess).
        i.workdir_free_bytes = 5 << 30;
        i.output_fs_free_bytes = 5 << 30;
        assert!(matches!(decide_compress(NumaMode::Auto, Some(&topo2()), &i).unwrap(), CompressPlan::Sharded { .. }));
    }
    #[test]
    fn fixed_more_than_nodes_or_threads_errors() {
        assert!(decide_compress(NumaMode::Fixed(3), Some(&topo2()), &base()).is_err());
        let mut i = base(); i.total_threads = 1;
        assert!(decide_compress(NumaMode::Fixed(2), Some(&topo2()), &i).is_err());
    }
    #[test]
    fn insufficient_node_memory_falls_back() {
        let mut t = topo2();
        for n in &mut t.nodes { n.mem_available_bytes = 1 << 20; }
        assert!(matches!(decide_compress(NumaMode::Auto, Some(&t), &base()).unwrap(), CompressPlan::InProcess));
    }
}
