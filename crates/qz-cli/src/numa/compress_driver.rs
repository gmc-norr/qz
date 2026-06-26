//! NUMA compress driver (spec §7): prechecks → topology → plan → gate → split →
//! spawn pinned `_numa-compress-worker` procs → poll → merge. Inverse of the decode
//! driver. Single-end in Phase 1; paired branch added in Task 13.

use crate::numa::compress_gate::{self, CompressGateInputs, CompressPlan};
use crate::numa::gate::{assign_node_threads, IoKind};
use crate::numa::runtime::{self, ShardGuard};
use crate::numa::splitter;
use crate::numa::{topology, NumaMode};
use anyhow::Result;
use qz_lib::cli::CompressConfig;
use qz_lib::compression::ByteRange;
use std::path::PathBuf;
use std::process::Command;

/// Entry from the CLI Compress arm. Decides in-process vs sharded, runs it.
pub fn compress_with_numa(mode: NumaMode, config: &CompressConfig, build_index: bool) -> Result<()> {
    // Reference require-by-default preflight (spec §3.3): ensure a usable index for
    // input[0]'s profile BEFORE the in-process / off path runs (those return below
    // without sharding). The SHARDED path additionally ensures each shard's profile
    // in run_sharded_core (Step 5b). A missing index is a clean parent error (not an
    // opaque worker failure); `--build-index` builds here in the parent.
    if let Some(r) = &config.reference
        && let Some(first) = config.input.first() {
            let read_len = qz_lib::compression::peek_first_read_length(first)?;
            qz_lib::compression::ensure_reference_index(
                &r.reference,
                read_len,
                r.reference_fast,
                config.threads,
                build_index,
            )?;
        }

    // 0) off / Fixed(1) → straight to qz-lib.
    if mode == NumaMode::Off || mode == NumaMode::Fixed(1) {
        return qz_lib::compression::compress(config);
    }

    // 1) Invalid input counts are a CONFIG error, not a NUMA matter — route to the
    //    in-process compress so the user gets its canonical "one or two input files"
    //    rejection, for BOTH auto AND fixed. Without this, a 3-input config would
    //    reach the sharded path which only ever inspects input[0]/[1].
    if config.input.is_empty() || config.input.len() > 2 {
        return qz_lib::compression::compress(config);
    }

    // Mode-sensitive bail: auto → correct in-process compress, fixed → hard error.
    let fallback_or_err =
        |msg: &str| runtime::fallback_or_err(mode, msg, || qz_lib::compression::compress(config));

    // 2) UNCONDITIONAL mode-routing prechecks (spec §3.1) — run BEFORE the forced
    //    hook so it can't force-shard an unsupported mode.
    if let Some(reason) = unsupported_mode_reason(config)? {
        return fallback_or_err(&reason);
    }

    // 3) Forced-shard test hook (debug only): N-way UNPINNED shard, bypass topology
    //    + gate (children get QZ_NUMA_NO_PIN=1). Still runs prechecks (above) + split.
    if let Some(n) = runtime::test_hook("QZ_NUMA_FORCE_WORKERS").and_then(|s| s.parse::<u32>().ok())
        && n >= 2
    {
        // One fake node; all N workers map to it via `idx % nodes.len()`, and the hook
        // sets QZ_NUMA_NO_PIN=1 so the node id is never used for pinning.
        let nodes = vec![topology::NumaNode {
            id: 0,
            cpus: (0..config.threads as u32).collect(),
            mem_available_bytes: u64::MAX,
        }];
        return run_sharded(config, n, &nodes, mode, build_index);
    }

    // 4) Topology.
    let topo = topology::detect();
    if !topo.as_ref().map(|t| t.has_real_asymmetry).unwrap_or(false) {
        return fallback_or_err("no usable NUMA topology (need >=2 asymmetric nodes on Linux)");
    }

    // 5) Gate.
    let input_bytes = config.input.iter().filter_map(|p| std::fs::metadata(p).ok().map(|m| m.len())).sum();
    let out_dir = config.output.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .map(std::path::Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));
    let inputs = CompressGateInputs {
        input_bytes,
        total_threads: config.threads,
        io: IoKind::Files, // stdio excluded by the precheck
        allowed_cpus: topology::allowed_cpus(),
        allowed_mem_nodes: topology::allowed_mem_nodes(),
        worker_peak_rss_bytes: worker_peak_rss_bound(config),
        workdir_free_bytes: runtime::free_space(&config.working_dir).unwrap_or(0),
        output_fs_free_bytes: runtime::free_space(&out_dir).unwrap_or(0),
        workdir_output_same_fs: runtime::same_device(&config.working_dir, &out_dir),
    };
    match compress_gate::decide_compress(mode, topo.as_ref(), &inputs)? {
        CompressPlan::InProcess => qz_lib::compression::compress(config),
        CompressPlan::Sharded { worker_count, nodes } => run_sharded(config, worker_count, &nodes, mode, build_index),
    }
}

/// Dispatch to the single-end or paired sharded driver by input count (2 = paired; the
/// `unsupported_mode_reason` precheck has already excluded reference, FASTA, gzip, stdio).
fn run_sharded(config: &CompressConfig, worker_count: u32, nodes: &[topology::NumaNode], mode: NumaMode, build_index: bool) -> Result<()> {
    if config.input.len() == 2 {
        run_sharded_paired(config, worker_count, nodes, mode, build_index)
    } else {
        run_sharded_single(config, worker_count, nodes, mode, build_index)
    }
}

/// §3.1 precheck: return Some(reason) for a currently-supported mode the sharding
/// path does not handle (caller routes auto→in-process / fixed→error). Order
/// matters; cheapest first.
fn unsupported_mode_reason(config: &CompressConfig) -> Result<Option<String>> {
    // stdio
    if qz_lib::cli::is_stdio_path(&config.output)
        || config.input.iter().any(|p| qz_lib::cli::is_stdio_path(p))
    {
        return Ok(Some("stdin/stdout cannot be sharded".into()));
    }
    // gzip input — peek EVERY input file (covers R2-only-gzipped paired).
    for p in &config.input {
        if file_is_gzipped(p)? {
            return Ok(Some(format!("gzipped input {} cannot be byte-range sharded", p.display())));
        }
    }
    // Reference compression (`--reference`, single-end type 4 / paired type 2) IS sharded:
    // the byte-range workers each produce a reference part archive, stitched by
    // `merge_reference_archives_to_path` (re-derives the coverage globals over the union of
    // covered intervals). gzip input is already rejected above (the bounded reference reader
    // is plain-file only); the reference precheck below keeps lossy/cluster/ultra out.
    // --cluster reorders reads GLOBALLY (one bucket-sort over the whole input), so a
    // byte-range split + per-shard compress + merge would interleave two independently
    // clustered halves — destroying the global cluster order AND producing a non-cluster
    // archive. Mirror the decode-side guard (numa/driver.rs): auto falls back to the
    // correct in-process cluster compress; fixed --numa N errors loudly.
    if config.cluster.is_some() {
        return Ok(Some("--cluster reorders reads globally and is not NUMA-shardable".into()));
    }
    // Paired-end (2 inputs, non-reference) IS sharded via run_sharded_paired (record-index
    // lockstep split + paired part archives + paired merge). reference (also 2-input) is
    // rejected above.
    if config.fasta {
        return Ok(Some("FASTA input is not NUMA-sharded (multiline records)".into()));
    }
    if config.ultra.is_some() {
        return Ok(Some("--ultra is not NUMA-sharded in this release".into()));
    }
    // Lossy quality modes (IlluminaBin, Discard) are out of scope for NUMA v1: the
    // worker rebuilds CompressConfig with quality_mode = Lossless, so a lossy mode
    // must NOT reach the sharding path. (no_quality WITH Lossless is fine — it's
    // passed through via --no-quality.) Discard sets no_quality but is also != Lossless.
    if config.quality_mode != qz_lib::cli::QualityMode::Lossless {
        return Ok(Some(format!(
            "quality mode {:?} is not NUMA-sharded (lossless only in this release)",
            config.quality_mode
        )));
    }
    Ok(None)
}

/// True iff the file begins with the gzip magic (0x1f 0x8b).
fn file_is_gzipped(p: &std::path::Path) -> Result<bool> {
    use std::io::Read;
    let mut f = match std::fs::File::open(p) { Ok(f) => f, Err(_) => return Ok(false) };
    let mut magic = [0u8; 2];
    match f.read(&mut magic) {
        Ok(2) => Ok(magic == [0x1f, 0x8b]),
        _ => Ok(false),
    }
}

/// Rough per-worker compress peak RSS bound (working set, constant in read count).
/// Conservative fixed bound for v1.
fn worker_peak_rss_bound(config: &CompressConfig) -> u64 {
    // TODO(numa): scale this with per-worker thread count (real RSS follows -t); a
    // fixed bound only gates WHETHER to shard (auto falls back / fixed errors), so an
    // optimistic value risks an unnecessary shard attempt, never corruption.
    let base = 4u64 << 30;
    // Reference shards each load the strobemer index — roughly the FASTA size again on
    // top of the base working set. Add the reference file size so a large reference
    // gates more conservatively (auto falls back in-process instead of OOM-attempting).
    if let Some(r) = &config.reference {
        let ref_bytes = std::fs::metadata(&r.reference).map(|m| m.len()).unwrap_or(0);
        return base.saturating_add(ref_bytes.saturating_mul(2));
    }
    base
}


fn run_sharded_single(
    config: &CompressConfig,
    worker_count: u32,
    nodes: &[topology::NumaNode],
    mode: NumaMode,
    build_index: bool,
) -> Result<()> {
    let fallback_or_err =
        |msg: &str| runtime::fallback_or_err(mode, msg, || qz_lib::compression::compress(config));
    // Split (mode-sensitive on too-few-records). One range per worker (R1 only).
    let t_split0 = std::time::Instant::now();
    let ranges = match splitter::split_single_end(&config.input[0], worker_count) {
        Ok(r) => r,
        Err(e) if e.downcast_ref::<splitter::SplitTooFewRecords>().is_some() => {
            return fallback_or_err("input has too few records to shard");
        }
        Err(e) => return Err(e),
    };
    let split_elapsed = t_split0.elapsed();
    let ranges_per_worker: Vec<Vec<ByteRange>> = ranges.into_iter().map(|r| vec![r]).collect();
    run_sharded_core(config, &ranges_per_worker, split_elapsed, nodes, mode, build_index)
}

fn run_sharded_paired(
    config: &CompressConfig,
    worker_count: u32,
    nodes: &[topology::NumaNode],
    mode: NumaMode,
    build_index: bool,
) -> Result<()> {
    let fallback_or_err =
        |msg: &str| runtime::fallback_or_err(mode, msg, || qz_lib::compression::compress(config));
    // Record-index lockstep split of BOTH mates (single memchr pass per file). Two ranges
    // per worker (R1, R2). Mismatched counts / too-few-records are mode-sensitive.
    let t_split0 = std::time::Instant::now();
    let pairs = match splitter::split_paired(&config.input[0], &config.input[1], worker_count) {
        Ok(r) => r,
        Err(e) if e.downcast_ref::<splitter::SplitTooFewRecords>().is_some() => {
            return fallback_or_err("input has too few records to shard");
        }
        Err(e) if e.downcast_ref::<splitter::MateCountMismatch>().is_some() => {
            return fallback_or_err("paired inputs have different record counts (cannot shard)");
        }
        Err(e) => return Err(e),
    };
    let split_elapsed = t_split0.elapsed();
    let ranges_per_worker: Vec<Vec<ByteRange>> = pairs.into_iter().map(|p| vec![p.r1, p.r2]).collect();
    run_sharded_core(config, &ranges_per_worker, split_elapsed, nodes, mode, build_index)
}

/// Shared spawn → poll → merge core for single (1 range/worker) and paired (2 ranges/worker).
/// Each `ranges_per_worker[i]` aligns positionally with `config.input` (R1 = index 0, R2 =
/// index 1) so the worker rebuilds the matching `CompressConfig.input`. `split_elapsed` is the
/// caller's measured split time (for QZ_NUMA_TIMING — the paired split is the heaviest serial
/// phase, so it must be attributed).
fn run_sharded_core(
    config: &CompressConfig,
    ranges_per_worker: &[Vec<ByteRange>],
    split_elapsed: std::time::Duration,
    nodes: &[topology::NumaNode],
    mode: NumaMode,
    build_index: bool,
) -> Result<()> {
    let worker_count = ranges_per_worker.len() as u32;

    // Reference: ensure an index for EVERY shard's read-length profile before
    // spawning (workers peek their own shard's first read). ensure_reference_index
    // is idempotent across shards in the same profile; it errors (or builds, under
    // --build-index) once per distinct profile. Within a profile, canonical params
    // + naming mean all shards share one index — only reads crossing bucket
    // boundaries across shards need more than one.
    if let Some(r) = &config.reference {
        for ranges in ranges_per_worker {
            if let Some(r1) = ranges.first()
                && let Some(read_len) =
                    qz_lib::compression::peek_read_length_in_range(&config.input[0], r1.start, r1.end)?
                {
                    qz_lib::compression::ensure_reference_index(
                        &r.reference,
                        read_len,
                        r.reference_fast,
                        config.threads,
                        build_index,
                    )?;
                }
        }
    }

    // Phase timing (QZ_NUMA_TIMING): isolate the serial driver fraction (split +
    // plan-resolve prelude + merge tail) that wraps the parallel worker core. The
    // manual `numactl` recipe pays none of it (each process resolves its plan
    // concurrently and never merges), so this is the Amdahl tax that caps the sharded
    // speedup. Kept as a gated diagnostic for future NUMA work.
    let timing = std::env::var_os("QZ_NUMA_TIMING").is_some();
    let t0 = std::time::Instant::now(); // origin AFTER the split (timed by the caller)

    // Resolve the canonical plan ONCE from the input head (per-mate codecs for paired).
    let plan = match qz_lib::compression::resolve_plan_override(config) {
        Ok(p) => p,
        Err(e) => return match mode { NumaMode::Fixed(_) => Err(e), _ => qz_lib::compression::compress(config) },
    };
    let plan_json = serde_json::to_string(&plan)?;
    let t_plan = t0.elapsed();

    // Preflight: refuse pre-existing output BEFORE spawning (mirror in-process).
    if !config.force && config.output.exists() {
        anyhow::bail!("Output file already exists: {}\nUse --force to overwrite", config.output.display());
    }

    let exe = std::env::current_exe()
        .map_err(|e| anyhow::anyhow!("cannot locate self for re-exec: {e}"))?;
    let advanced_json = serde_json::to_string(&config.advanced)?;

    let tmp = tempfile::Builder::new()
        .prefix(".qz_numa_compress_")
        .tempdir_in(&config.working_dir)
        .map_err(|e| anyhow::anyhow!("create NUMA part dir in {:?}: {e}", config.working_dir))?;
    let part_paths: Vec<PathBuf> = (0..ranges_per_worker.len())
        .map(|i| tmp.path().join(format!("w{i}.qz")))
        .collect();

    let fallback_or_err =
        |msg: &str| runtime::fallback_or_err(mode, msg, || qz_lib::compression::compress(config));

    let per_worker = compress_gate::per_worker_threads(config.threads, worker_count);
    let mut guard = ShardGuard::with_tmp(tmp);
    for (i, ranges) in ranges_per_worker.iter().enumerate() {
        let (node_id, threads) = assign_node_threads(i, nodes, per_worker);
        let mut cmd = worker_command(&exe, config, ranges, node_id, threads, i, &part_paths[i], &plan_json, &advanced_json);
        match cmd.spawn() {
            Ok(child) => guard.children.push(child),
            Err(e) => return Err(anyhow::anyhow!("spawn compress worker {i}: {e}")),
        }
    }
    let t_spawn = t0.elapsed();

    // Poll to completion (kills survivors on the first failure). Unlike the decode driver,
    // pin failures (exit 3) need no special-casing here: under auto `fallback_or_err` falls
    // back in-process on ANY worker failure, under fixed it errors on ANY failure — so all
    // non-zero exits are treated alike (the collected codes are ignored).
    let failed = !runtime::wait_for_workers(&mut guard.children)?.is_empty();
    if failed {
        drop(guard);
        return fallback_or_err("a NUMA compress worker failed");
    }
    let t_workers = t0.elapsed();

    // Merge while the guard is STILL ALIVE (parts present), capture, then disarm.
    // Reference shards merge through the reference-aware path (re-derives the coverage
    // globals over the union of covered intervals from the FASTA); all other modes use
    // the verbatim-frame-copy merge.
    let result = match &config.reference {
        Some(r) => qz_lib::compression::merge_reference_archives_to_path(
            &part_paths,
            &r.reference,
            &config.output,
            config.force,
        ),
        None => {
            qz_lib::compression::merge_archives_to_path(&part_paths, &config.output, config.force)
        }
    };
    let t_merge = t0.elapsed();
    if timing {
        eprintln!(
            "QZ_NUMA_TIMING split={:.3}s plan_prelude={:.3}s spawn={:.3}s workers={:.3}s merge={:.3}s | serial_prelude={:.3}s parallel_core={:.3}s serial_merge={:.3}s total={:.3}s",
            split_elapsed.as_secs_f64(),
            t_plan.as_secs_f64(),
            (t_spawn - t_plan).as_secs_f64(),
            (t_workers - t_spawn).as_secs_f64(),
            (t_merge - t_workers).as_secs_f64(),
            (split_elapsed + t_plan).as_secs_f64(),
            (t_workers - t_spawn).as_secs_f64(),
            (t_merge - t_workers).as_secs_f64(),
            (split_elapsed + t_merge).as_secs_f64(),
        );
    }
    guard.armed = false;
    drop(guard);
    match result {
        Ok(()) => Ok(()),
        Err(e) => match mode {
            NumaMode::Fixed(_) => Err(e),
            _ => qz_lib::compression::compress(config), // merge artifact → fall back
        },
    }
}

#[allow(clippy::too_many_arguments)]
fn worker_command(
    exe: &std::path::Path,
    config: &CompressConfig,
    ranges: &[ByteRange],
    node: u32,
    threads: usize,
    worker_id: usize,
    out_part: &std::path::Path,
    plan_json: &str,
    advanced_json: &str,
) -> Command {
    let r1 = &ranges[0];
    let mut cmd = Command::new(exe);
    cmd.arg("_numa-compress-worker")
        .arg("--input").arg(&config.input[0])
        .arg("--byte-start").arg(r1.start.to_string())
        .arg("--byte-end").arg(r1.end.to_string())
        .arg("--node").arg(node.to_string())
        .arg("--threads").arg(threads.to_string())
        .arg("--output").arg(out_part)
        .arg("--worker-id").arg(worker_id.to_string())
        .arg("--plan").arg(plan_json)
        .arg("--working-dir").arg(&config.working_dir)
        .arg("--advanced").arg(advanced_json);
    // Paired: forward R2 (input[1]) + its byte range. The worker rebuilds the 2-input
    // CompressConfig and dispatches compress_byte_range to the paired branch.
    if ranges.len() == 2 {
        let r2 = &ranges[1];
        cmd.arg("--input2").arg(&config.input[1])
            .arg("--byte-start2").arg(r2.start.to_string())
            .arg("--byte-end2").arg(r2.end.to_string());
    }
    if config.no_quality { cmd.arg("--no-quality"); }
    if config.advanced.bsc_static { cmd.arg("--fast"); }
    // Reference shard: forward the FASTA + window/fast so the worker rebuilds
    // config.reference and compress_byte_range takes the reference arm. The reference
    // front header is fixed metadata, so no plan needs forwarding for identity.
    if let Some(r) = &config.reference {
        cmd.arg("--reference").arg(&r.reference)
            .arg("--reference-window").arg(r.reference_window.to_string());
        if r.reference_fast {
            cmd.arg("--reference-fast");
        }
    }
    runtime::apply_worker_env(&mut cmd, &["QZ_NUMA_FAIL", "QZ_NUMA_SLOW", "QZ_NUMA_PINFAIL"]);
    cmd
}

