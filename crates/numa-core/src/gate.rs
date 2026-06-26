use crate::topology::Topology;
use crate::NumaMode;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IoKind { Files, Stdio }

#[derive(Clone, Debug)]
pub struct WorkerRange {
    pub node: u32,
    pub chunk_start: u32,
    pub chunk_end: u32,
    pub threads: usize,
}

#[derive(Clone, Debug)]
pub enum Plan {
    InProcess,
    Sharded { ranges: Vec<WorkerRange> },
}

#[derive(Clone, Debug)]
pub struct GateInputs {
    pub archive_type: u8,
    pub num_chunks: u32,
    pub per_chunk_reads: Vec<u64>,
    pub est_decoded_bytes: u64,
    pub total_threads: usize,
    pub io: IoKind,
    pub workdir_free_bytes: u64,
    pub worker_peak_rss_bytes: u64,
    pub allowed_cpus: Vec<u32>,
    pub allowed_mem_nodes: Vec<u32>,
    /// false during Increment A (reference not yet shardable).
    pub increment_b_ready: bool,
    /// false iff the archive cannot be range-decoded (mixed-codec legacy paired).
    pub shardable: bool,
    /// true when this decode takes the single-end direct-write path (workers write
    /// disjoint regions into one pre-sized shared output file — no part files, no
    /// assembly temp): single-end + non-gzip + a present ChunkDecodedSizes table.
    pub direct_write: bool,
    /// Free bytes on the filesystem that holds the FINAL output (statvfs of its
    /// parent dir, since the output itself may not exist yet).
    pub output_fs_free_bytes: u64,
    /// true iff the working dir and the output dir live on the SAME filesystem
    /// (device id). Only meaningful on the part-file path: when same-fs the part
    /// files AND the assembly temp coexist on one fs (~2×D); when cross-fs each
    /// needs ~1×D independently.
    pub workdir_output_same_fs: bool,
}

/// Minimum estimated decoded / input size to bother sharding under `auto` (256 MiB).
/// Shared by the decode gate (keyed on `est_decoded_bytes`) and the compress gate (keyed
/// on `input_bytes`).
pub const MIN_WORKLOAD_BYTES: u64 = 256 * 1024 * 1024;

/// Shared gate preamble: require a real-asymmetry topology and non-stdio I/O. Returns the
/// topology ref, or the fallback reason for the caller's `bail_or_fallback`. (The decode
/// gate interposes its archive-type / shardable checks between this and the memory filter;
/// the compress gate goes straight to the filter.)
pub fn require_shardable_topology(
    topo: Option<&Topology>,
    io: IoKind,
) -> Result<&Topology, String> {
    let topo = match topo {
        Some(t) if t.has_real_asymmetry => t,
        _ => return Err("no usable NUMA topology (need >=2 asymmetric nodes on Linux)".into()),
    };
    if io == IoKind::Stdio {
        return Err("stdin/stdout I/O cannot be sharded".into());
    }
    Ok(topo)
}

/// Shared memory filter: usable nodes whose available memory covers the per-worker RSS
/// bound. Returns the filtered nodes (≥2) or the fallback reason.
pub fn mem_filtered_nodes(
    topo: &Topology,
    allowed_cpus: &[u32],
    allowed_mem_nodes: &[u32],
    worker_peak_rss: u64,
) -> Result<Vec<crate::topology::NumaNode>, String> {
    let mut usable = topo.usable_nodes(allowed_cpus, allowed_mem_nodes);
    usable.retain(|n| n.mem_available_bytes >= worker_peak_rss);
    if usable.len() < 2 {
        return Err("fewer than 2 usable nodes after cpuset/mem/RSS filtering".into());
    }
    Ok(usable)
}

/// Shared fixed-`--numa N` validation: N must not exceed the usable node count or the
/// thread budget. Hard error (only `Fixed` mode reaches it).
pub fn validate_fixed_workers(
    n: u32,
    usable_len: usize,
    total_threads: usize,
) -> anyhow::Result<()> {
    if n as usize > usable_len {
        anyhow::bail!("--numa {n}: only {usable_len} usable nodes");
    }
    if n as usize > total_threads {
        anyhow::bail!("--numa {n}: only {total_threads} threads (-t)");
    }
    Ok(())
}

pub fn decide(mode: NumaMode, topo: Option<&Topology>, i: &GateInputs) -> anyhow::Result<Plan> {
    // 0) off / Fixed(1) -> in-process.
    if mode == NumaMode::Off || mode == NumaMode::Fixed(1) {
        return Ok(Plan::InProcess);
    }
    let fixed_n = match mode { NumaMode::Fixed(n) => Some(n), _ => None };
    let bail_or_fallback = |msg: &str| -> anyhow::Result<Plan> {
        if fixed_n.is_some() { anyhow::bail!("--numa: {msg}") } else { Ok(Plan::InProcess) }
    };

    let topo = match require_shardable_topology(topo, i.io) {
        Ok(t) => t,
        Err(msg) => return bail_or_fallback(&msg),
    };
    if matches!(i.archive_type, 2 | 4) && !i.increment_b_ready {
        return bail_or_fallback("reference NUMA decode requires Increment B");
    }
    if !i.shardable {
        return bail_or_fallback("archive is not shardable (mixed-codec legacy paired path)");
    }
    let usable = match mem_filtered_nodes(topo, &i.allowed_cpus, &i.allowed_mem_nodes, i.worker_peak_rss_bytes) {
        Ok(u) => u,
        Err(msg) => return bail_or_fallback(&msg),
    };
    // Filesystem capacity. Auto falls back to in-process (which needs only ~D at the
    // output) when a check fails, so erring conservative here never makes auto worse
    // than off. (F3)
    //
    //  * direct-write: workers write disjoint regions into the pre-sized shared output
    //    file(s) on the OUTPUT fs (one for single-end, two R1/R2 for paired; their sizes
    //    sum to D) — no part files, no assembly temp — so it needs only ~1×D there. A
    //    sanity cap rejects an implausibly large (corrupt-table) size.
    //  * part-file path, same-fs: parts (~D) and the assembly temp (~D) coexist on the
    //    single shared fs during the final concat, so ~2×D is needed.
    //  * part-file path, cross-fs: the parts live in the working dir (~D) and the
    //    assembly temp lives next to the output (~D), each on its own fs → ~1×D each.
    let d = i.est_decoded_bytes;
    if i.direct_write {
        if d > crate::MAX_DIRECT_OUTPUT_BYTES {
            return bail_or_fallback("decoded size implausibly large (corrupt table?)");
        }
        if i.output_fs_free_bytes < d {
            return bail_or_fallback("output filesystem lacks ~1x decoded size (direct-write)");
        }
    } else if i.workdir_output_same_fs {
        if i.output_fs_free_bytes < d.saturating_mul(2) {
            return bail_or_fallback("filesystem lacks ~2x decoded size (parts + assembly temp share one fs)");
        }
    } else {
        if i.workdir_free_bytes < d {
            return bail_or_fallback("working dir lacks ~1x decoded size for part files");
        }
        if i.output_fs_free_bytes < d {
            return bail_or_fallback("output filesystem lacks ~1x decoded size for assembly temp");
        }
    }
    let cap = (usable.len() as u32)
        .min(i.num_chunks)
        .min(i.total_threads.max(1) as u32);
    let workers = match fixed_n {
        Some(n) => {
            validate_fixed_workers(n, usable.len(), i.total_threads)?;
            n.min(i.num_chunks)
        }
        None => {
            if i.est_decoded_bytes < MIN_WORKLOAD_BYTES { return Ok(Plan::InProcess); }
            cap
        }
    };
    if workers < 2 {
        return Ok(Plan::InProcess);
    }
    let ranges = partition_balanced(&i.per_chunk_reads, workers, &usable, i.total_threads);
    Ok(Plan::Sharded { ranges })
}

/// Assign a worker's NUMA node + thread budget: round-robin node by worker index,
/// `total_threads / workers` threads capped by the node's CPU count (min 1).
/// Shared by decode's `partition_balanced` and compress's `decide_compress`.
pub fn assign_node_threads(
    idx: usize,
    nodes: &[crate::topology::NumaNode],
    per_worker_threads: usize,
) -> (u32, usize) {
    let node = &nodes[idx % nodes.len()];
    (node.id, per_worker_threads.min(node.cpus.len()).max(1))
}

/// Split `[0, n)` chunks into EXACTLY `workers` contiguous, non-empty ranges
/// balanced by read count, assigning each to a usable node. Precondition: `1 <= workers <= n`.
pub fn partition_balanced(
    per_chunk_reads: &[u64],
    workers: u32,
    usable: &[crate::topology::NumaNode],
    total_threads: usize,
) -> Vec<WorkerRange> {
    let n = per_chunk_reads.len() as u32;
    debug_assert!(workers >= 1 && workers <= n, "gate must clamp workers to [1, num_chunks]");
    let total: u64 = per_chunk_reads.iter().sum();
    let target = (total / u64::from(workers)).max(1);
    let per_worker_threads = (total_threads / workers as usize).max(1);

    let mut ranges: Vec<WorkerRange> = Vec::with_capacity(workers as usize);
    let mut start = 0u32;
    let mut acc = 0u64;
    let push = |ranges: &mut Vec<WorkerRange>, start: u32, end: u32| {
        let (node, threads) = assign_node_threads(ranges.len(), usable, per_worker_threads);
        ranges.push(WorkerRange { node, chunk_start: start, chunk_end: end, threads });
    };

    for c in 0..n {
        acc += per_chunk_reads[c as usize];
        let made = ranges.len() as u32;
        if made == workers - 1 {
            break;
        }
        let workers_left_after_close = workers - made - 1;
        let chunks_left_after = n - (c + 1);
        let must_close = chunks_left_after == workers_left_after_close;
        let want_close = acc >= target && chunks_left_after > workers_left_after_close;
        if must_close || want_close {
            push(&mut ranges, start, c + 1);
            start = c + 1;
            acc = 0;
        }
    }
    push(&mut ranges, start, n);
    debug_assert_eq!(ranges.len(), workers as usize);
    debug_assert!(ranges.iter().all(|r| r.chunk_end > r.chunk_start));
    debug_assert_eq!(ranges.first().unwrap().chunk_start, 0);
    debug_assert_eq!(ranges.last().unwrap().chunk_end, n);
    ranges
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::{NumaNode, Topology};

    fn topo2() -> Topology {
        Topology {
            nodes: vec![
                NumaNode { id: 0, cpus: (0..18).collect(), mem_available_bytes: 64 << 30 },
                NumaNode { id: 1, cpus: (18..36).collect(), mem_available_bytes: 64 << 30 },
            ],
            has_real_asymmetry: true,
        }
    }

    fn base_inputs() -> GateInputs {
        GateInputs {
            archive_type: 0,
            num_chunks: 8,
            per_chunk_reads: vec![1000; 8],
            est_decoded_bytes: 4 << 30,
            total_threads: 36,
            io: IoKind::Files,
            workdir_free_bytes: 100 << 30,
            worker_peak_rss_bytes: 2 << 30,
            allowed_cpus: (0..36).collect(),
            allowed_mem_nodes: vec![0, 1],
            increment_b_ready: false,
            shardable: true,
            direct_write: false,
            output_fs_free_bytes: 100 << 30,
            workdir_output_same_fs: false,
        }
    }

    #[test]
    fn non_shardable_falls_back_auto_errors_fixed() {
        let mut i = base_inputs(); i.shardable = false;
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
        assert!(decide(NumaMode::Fixed(2), Some(&topo2()), &i).is_err());
    }

    #[test]
    fn off_is_always_in_process() {
        let p = decide(NumaMode::Off, Some(&topo2()), &base_inputs()).unwrap();
        assert!(matches!(p, Plan::InProcess));
    }

    #[test]
    fn auto_shards_two_nodes() {
        let p = decide(NumaMode::Auto, Some(&topo2()), &base_inputs()).unwrap();
        match p {
            Plan::Sharded { ranges } => {
                assert_eq!(ranges.len(), 2);
                assert_eq!(ranges[0].chunk_start, 0);
                assert_eq!(ranges.last().unwrap().chunk_end, 8);
                assert_eq!(ranges[0].chunk_end, ranges[1].chunk_start);
                assert!(ranges.iter().all(|r| r.threads >= 1 && r.threads <= 18));
            }
            _ => panic!("expected sharded"),
        }
    }

    #[test]
    fn single_node_falls_back() {
        let mut t = topo2();
        t.nodes.truncate(1);
        t.has_real_asymmetry = false;
        assert!(matches!(decide(NumaMode::Auto, Some(&t), &base_inputs()).unwrap(), Plan::InProcess));
    }

    #[test]
    fn non_linux_or_no_topo_falls_back_in_auto_errors_in_fixed() {
        assert!(matches!(decide(NumaMode::Auto, None, &base_inputs()).unwrap(), Plan::InProcess));
        assert!(decide(NumaMode::Fixed(2), None, &base_inputs()).is_err());
    }

    #[test]
    fn stdio_falls_back_auto_errors_fixed() {
        let mut i = base_inputs(); i.io = IoKind::Stdio;
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
        assert!(decide(NumaMode::Fixed(2), Some(&topo2()), &i).is_err());
    }

    #[test]
    fn reference_before_increment_b_falls_back_auto_errors_fixed() {
        let mut i = base_inputs(); i.archive_type = 2;
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
        assert!(decide(NumaMode::Fixed(2), Some(&topo2()), &i).is_err());
    }

    #[test]
    fn single_end_reference_type4_before_increment_b_falls_back_auto_errors_fixed() {
        let mut i = base_inputs();
        i.archive_type = 4;
        i.increment_b_ready = false;
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
        assert!(decide(NumaMode::Fixed(2), Some(&topo2()), &i).is_err());
    }

    #[test]
    fn single_end_reference_type4_shards_when_increment_b_ready() {
        let mut i = base_inputs();
        i.archive_type = 4;
        i.increment_b_ready = true;
        match decide(NumaMode::Auto, Some(&topo2()), &i).unwrap() {
            Plan::Sharded { ranges } => assert_eq!(ranges.len(), 2),
            _ => panic!("type-4 reference must shard once Increment B is ready"),
        }
    }

    #[test]
    fn subthreshold_workload_falls_back_in_auto_only() {
        let mut i = base_inputs(); i.est_decoded_bytes = 1 << 20;
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
        assert!(matches!(decide(NumaMode::Fixed(2), Some(&topo2()), &i).unwrap(), Plan::Sharded { .. }));
    }

    #[test]
    fn fixed_more_than_nodes_errors_more_than_threads_errors() {
        assert!(decide(NumaMode::Fixed(3), Some(&topo2()), &base_inputs()).is_err());
        let mut i = base_inputs(); i.total_threads = 1;
        assert!(decide(NumaMode::Fixed(2), Some(&topo2()), &i).is_err());
    }

    #[test]
    fn fixed_one_is_in_process() {
        assert!(matches!(decide(NumaMode::Fixed(1), Some(&topo2()), &base_inputs()).unwrap(), Plan::InProcess));
        assert!(matches!(decide(NumaMode::Fixed(1), None, &base_inputs()).unwrap(), Plan::InProcess));
    }

    #[test]
    fn fixed_clamped_to_one_chunk_runs_in_process_not_error() {
        let mut i = base_inputs(); i.num_chunks = 1;
        assert!(matches!(decide(NumaMode::Fixed(2), Some(&topo2()), &i).unwrap(), Plan::InProcess));
    }

    #[test]
    fn workers_clamped_to_chunks() {
        let mut i = base_inputs(); i.num_chunks = 1;
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
    }

    #[test]
    fn insufficient_node_memory_falls_back() {
        let mut t = topo2();
        for n in &mut t.nodes { n.mem_available_bytes = 1 << 20; }
        assert!(matches!(decide(NumaMode::Auto, Some(&t), &base_inputs()).unwrap(), Plan::InProcess));
    }

    #[test]
    fn insufficient_workdir_capacity_falls_back() {
        // Cross-fs part-file path (default): the working dir must hold ~1×D of part
        // files. Starve the working dir → fall back even though the output fs is fine.
        let mut i = base_inputs(); i.workdir_free_bytes = 1 << 20;
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
    }

    #[test]
    fn cross_fs_part_path_needs_1x_on_each_fs_independently() {
        // Default base_inputs is cross-fs (workdir_output_same_fs=false): the working
        // dir needs ~1×D for parts AND the output fs needs ~1×D for the assembly temp.
        let mut i = base_inputs();
        i.est_decoded_bytes = 4 << 30; // D = 4 GiB
        // Plenty of workdir, but output fs short of 1×D → fall back.
        i.workdir_free_bytes = 100 << 30;
        i.output_fs_free_bytes = 3 << 30; // < 1×D
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
        // Plenty of output fs, but workdir short of 1×D → fall back.
        i.output_fs_free_bytes = 100 << 30;
        i.workdir_free_bytes = 3 << 30; // < 1×D
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
        // Each fs has >1×D but NOT 2×D → still shards (cross-fs does not need 2×).
        i.workdir_free_bytes = 5 << 30;
        i.output_fs_free_bytes = 5 << 30;
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::Sharded { .. }));
    }

    #[test]
    fn same_fs_part_path_needs_2x_decoded_size() {
        // When workdir and output share one fs, parts (~D) + assembly temp (~D)
        // coexist → ~2×D needed on the single (output) fs.
        let mut i = base_inputs();
        i.workdir_output_same_fs = true;
        i.est_decoded_bytes = 4 << 30;          // D = 4 GiB
        i.output_fs_free_bytes = 5 << 30;        // > 1×D, < 2×D → fall back
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
        i.output_fs_free_bytes = 9 << 30;        // > 2×D → shards
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::Sharded { .. }));
    }

    #[test]
    fn direct_write_needs_only_1x_on_output_fs() {
        let mut i = base_inputs();
        i.direct_write = true;
        i.est_decoded_bytes = 4 << 30;          // D = 4 GiB
        // Working dir empty is irrelevant for direct-write; only the output fs matters.
        i.workdir_free_bytes = 0;
        i.output_fs_free_bytes = 3 << 30;        // < 1×D → fall back
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
        i.output_fs_free_bytes = 5 << 30;        // > 1×D (no 2× needed) → shards
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::Sharded { .. }));
    }

    #[test]
    fn direct_write_rejects_implausible_decoded_size() {
        let mut i = base_inputs();
        i.direct_write = true;
        i.est_decoded_bytes = (1u64 << 46) + 1;  // over the sanity cap
        i.output_fs_free_bytes = u64::MAX;        // space is not the issue
        assert!(matches!(decide(NumaMode::Auto, Some(&topo2()), &i).unwrap(), Plan::InProcess));
        assert!(decide(NumaMode::Fixed(2), Some(&topo2()), &i).is_err());
    }

    fn assert_valid_partition(ranges: &[WorkerRange], workers: usize, n: u32) {
        assert_eq!(ranges.len(), workers, "must produce exactly `workers` ranges");
        assert!(ranges.iter().all(|r| r.chunk_end > r.chunk_start), "every range non-empty");
        assert_eq!(ranges[0].chunk_start, 0);
        assert_eq!(ranges.last().unwrap().chunk_end, n);
        for w in ranges.windows(2) {
            assert_eq!(w[0].chunk_end, w[1].chunk_start, "contiguous, no gap/overlap");
        }
    }

    #[test]
    fn partition_four_chunks_two_workers() {
        let usable = topo2().nodes;
        let r = partition_balanced(&[10, 10, 10, 10], 2, &usable, 36);
        assert_valid_partition(&r, 2, 4);
    }

    #[test]
    fn partition_four_chunks_four_workers_yields_four_singletons() {
        let usable = topo2().nodes;
        let r = partition_balanced(&[10, 10, 10, 10], 4, &usable, 36);
        assert_valid_partition(&r, 4, 4);
        assert!(r.iter().all(|x| x.chunk_end - x.chunk_start == 1));
    }

    #[test]
    fn partition_skewed_reads_still_covers_all() {
        let usable = topo2().nodes;
        let r = partition_balanced(&[1, 1, 1, 1000, 1, 1], 3, &usable, 36);
        assert_valid_partition(&r, 3, 6);
    }

    #[test]
    fn auto_partition_is_valid_end_to_end() {
        let mut i = base_inputs(); i.num_chunks = 5; i.per_chunk_reads = vec![100; 5];
        match decide(NumaMode::Auto, Some(&topo2()), &i).unwrap() {
            Plan::Sharded { ranges } => assert_valid_partition(&ranges, 2, 5),
            _ => panic!("expected sharded"),
        }
    }
}
