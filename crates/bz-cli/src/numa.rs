//! NUMA-aware decode sharding for `bz`, built on the shared `numa-core` crate.
//!
//! On NUMA hardware, `bz decompress --numa auto` re-execs one pinned worker process
//! per socket; each decodes a disjoint chunk range of the archive into a BGZF/BAM
//! part (pinned with `sched_setaffinity` + `set_mempolicy(MPOL_BIND)` so its decode
//! buffers fault node-local), and the driver concatenates the parts into the final
//! BAM. Off NUMA / sub-threshold inputs decode in-process with zero overhead.
//!
//! bz's output is a single BGZF stream of variable compressed size per chunk, so —
//! unlike qz single-end's direct-write — workers always take the **part-file +
//! assembly** path. The format-agnostic half (topology, pinning, the gate, the
//! spawn/poll/cleanup runtime) is reused verbatim from `numa-core`; only the chunk
//! layout, range decode (`bz_lib`), and the BGZF-aware assembly are bz-specific.

use anyhow::Result;
use bz_lib::{AdvancedOptions, CompressConfig, DecompressConfig};
use clap::Parser;
use numa_core::gate::{self, GateInputs, IoKind, Plan, WorkerRange};
use numa_core::runtime::{self, ShardGuard};
use numa_core::{pin, topology, NumaMode};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Hidden `_numa-decode-worker` subcommand args: decode one chunk range into one part.
#[derive(Parser, Debug)]
pub struct NumaWorkerArgs {
    #[arg(long)]
    input: PathBuf,
    #[arg(long)]
    chunk_start: u32,
    #[arg(long)]
    chunk_end: u32,
    #[arg(long)]
    node: u32,
    #[arg(long)]
    threads: usize,
    #[arg(long)]
    out_part: PathBuf,
    #[arg(long)]
    worker_id: u32,
}

/// The 28-byte BGZF EOF marker (an empty BGZF block) that noodles appends on
/// `finish()`. Each worker's part is a complete BGZF stream ending in this marker;
/// the assembler strips it from every part except the last so the parts concatenate
/// into one valid BAM. (SAM spec / htslib canonical value.)
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Rough decoded-BAM-size estimate = compressed archive × this. bz's BAM ratio is
/// ~2× (reference-free); 3 adds margin for the gate's workload threshold + the
/// part-file capacity check (over-estimating only makes the gate more conservative).
const BZ_DECODE_RATIO: u64 = 3;

/// Worker entry (hidden `_numa-decode-worker`): pin to the assigned node BEFORE
/// building the rayon pool (so pool threads fault node-local), then decode the chunk
/// range into the part. Exit codes match the cross-process ABI the driver classifies:
/// 0 success, 1 general error (via `main`), 3 pin failure (from `pin_worker_to_node`).
pub fn run_decode_worker(a: &NumaWorkerArgs) -> Result<()> {
    pin::worker_test_hooks(a.worker_id, "decode worker");
    pin::pin_worker_to_node(a.node);
    if a.threads > 0 {
        // Fresh process: this is the pool's first (and only) build.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(a.threads)
            .build_global();
    }
    bz_lib::compression::decode_chunk_range(
        &a.input,
        a.chunk_start,
        a.chunk_end,
        a.threads,
        &a.out_part,
    )?;
    Ok(())
}

/// CLI entry from the `Decompress` arm. Decides in-process vs sharded and runs it.
/// `total_threads` is the resolved `-t` budget (or detected parallelism); `force`
/// gates the final overwrite (the caller already checked it before any work).
pub fn decompress_with_numa(
    mode: NumaMode,
    config: &DecompressConfig,
    total_threads: usize,
    force: bool,
) -> Result<()> {
    // off / Fixed(1) → straight to the in-process decode.
    if mode == NumaMode::Off || mode == NumaMode::Fixed(1) {
        return bz_lib::decompress(config);
    }

    // TEST HOOK (debug-only): force an N-way shard regardless of topology/workload,
    // children UNPINNED (apply_worker_env sets QZ_NUMA_NO_PIN), so the whole
    // spawn → decode → assemble path is exercisable on any box.
    if let Some(n) =
        runtime::test_hook("QZ_NUMA_FORCE_WORKERS").and_then(|s| s.parse::<u32>().ok())
        && n >= 2
    {
        let layout = bz_lib::compression::read_bz_chunk_layout(&config.input)?;
        if layout.num_chunks >= n {
            let usable = vec![topology::NumaNode {
                id: 0,
                cpus: (0..total_threads.max(1) as u32).collect(),
                mem_available_bytes: u64::MAX,
            }];
            let ranges =
                gate::partition_balanced(&layout.per_chunk_reads, n, &usable, total_threads);
            return run_sharded(&ranges, config, mode, force);
        }
    }

    let fallback_or_err =
        |msg: &str| runtime::fallback_or_err(mode, msg, || bz_lib::decompress(config));

    // CHEAP prechecks before touching the archive: a real-asymmetry NUMA topology.
    let topo = topology::detect();
    if !topo.as_ref().map(|t| t.has_real_asymmetry).unwrap_or(false) {
        return fallback_or_err("no usable NUMA topology (need >=2 asymmetric nodes on Linux)");
    }

    // On NUMA hardware: pre-scan the chunk layout.
    let layout = match bz_lib::compression::read_bz_chunk_layout(&config.input) {
        Ok(l) => l,
        Err(e) => {
            return match mode {
                NumaMode::Fixed(_) => Err(e),
                _ => bz_lib::decompress(config),
            };
        }
    };

    let worker_peak_rss =
        bz_lib::compression::decode_peak_rss_bound(layout.num_records, layout.num_chunks);
    let est_decoded_bytes = std::fs::metadata(&config.input)
        .map(|m| m.len())
        .unwrap_or(0)
        .saturating_mul(BZ_DECODE_RATIO);

    let out_dir = output_dir(&config.output);
    let inputs = GateInputs {
        // bz has no archive_type concept; 0 keeps it out of the reference (2/4) gate
        // branch and on the generic part-file path.
        archive_type: 0,
        num_chunks: layout.num_chunks,
        per_chunk_reads: layout.per_chunk_reads.clone(),
        est_decoded_bytes,
        total_threads,
        io: IoKind::Files, // bz decompress is always file → file.
        workdir_free_bytes: runtime::free_space(&config.working_dir).unwrap_or(0),
        worker_peak_rss_bytes: worker_peak_rss,
        allowed_cpus: topology::allowed_cpus(),
        allowed_mem_nodes: topology::allowed_mem_nodes(),
        increment_b_ready: true,
        shardable: true,
        // BGZF output has a variable compressed size per chunk, so workers can't seek
        // into a pre-sized shared file: always the part-file + assembly path.
        direct_write: false,
        output_fs_free_bytes: runtime::free_space(&out_dir).unwrap_or(0),
        workdir_output_same_fs: runtime::same_device(&config.working_dir, &out_dir),
    };

    match gate::decide(mode, topo.as_ref(), &inputs)? {
        Plan::InProcess => bz_lib::decompress(config),
        Plan::Sharded { ranges } => run_sharded(&ranges, config, mode, force),
    }
}

/// Spawn one pinned worker per range, wait, then assemble the BGZF parts. On a pin
/// failure (`exit 3`) under `auto`, fall back to a correct in-process decode.
fn run_sharded(
    ranges: &[WorkerRange],
    config: &DecompressConfig,
    mode: NumaMode,
    force: bool,
) -> Result<()> {
    let exe = std::env::current_exe()
        .map_err(|e| anyhow::anyhow!("cannot locate self for re-exec: {e}"))?;

    let tmp = tempfile::Builder::new()
        .prefix(".bz_numa_")
        .tempdir_in(&config.working_dir)
        .map_err(|e| anyhow::anyhow!("create NUMA part dir in {:?}: {e}", config.working_dir))?;
    let parts: Vec<PathBuf> = (0..ranges.len())
        .map(|wi| tmp.path().join(format!("w{wi}.part")))
        .collect();

    let mut guard = ShardGuard::with_tmp(tmp);
    for (wi, r) in ranges.iter().enumerate() {
        if runtime::forced_spawn_fail(wi) {
            return Err(anyhow::anyhow!("test hook: forced spawn failure for worker {wi}"));
        }
        let mut cmd = worker_command(&exe, config, r, wi, &parts[wi]);
        match cmd.spawn() {
            Ok(child) => guard.children.push(child),
            Err(e) => return Err(anyhow::anyhow!("spawn worker {wi}: {e}")),
        }
    }

    runtime::write_spawn_log(guard.children.len());

    let fail_codes = runtime::wait_for_workers(&mut guard.children)?;
    if !fail_codes.is_empty() {
        let pin_failure = fail_codes.contains(&Some(3));
        drop(guard);
        if pin_failure && mode == NumaMode::Auto {
            return bz_lib::decompress(config);
        }
        anyhow::bail!("a NUMA decode worker failed; aborted (no output written)");
    }

    let result = assemble_bam_parts(&parts, &config.output, force);
    guard.armed = false;
    drop(guard);
    result
}

/// Build the worker re-exec `Command`.
fn worker_command(
    exe: &Path,
    config: &DecompressConfig,
    r: &WorkerRange,
    wi: usize,
    out_part: &Path,
) -> Command {
    let mut cmd = Command::new(exe);
    cmd.arg("_numa-decode-worker")
        .arg("--input")
        .arg(&config.input)
        .arg("--chunk-start")
        .arg(r.chunk_start.to_string())
        .arg("--chunk-end")
        .arg(r.chunk_end.to_string())
        .arg("--node")
        .arg(r.node.to_string())
        .arg("--threads")
        .arg(r.threads.to_string())
        .arg("--out-part")
        .arg(out_part)
        .arg("--worker-id")
        .arg(wi.to_string());
    runtime::apply_worker_env(&mut cmd, &["QZ_NUMA_FAIL", "QZ_NUMA_SLOW", "QZ_NUMA_PINFAIL"]);
    cmd
}

/// Concatenate the BGZF/BAM parts (in worker order) into the final BAM via a sibling
/// temp + atomic rename. Every part is a complete BGZF stream ending in the 28-byte
/// EOF marker; that marker is stripped from all but the LAST part so the result is one
/// continuous, valid BGZF stream (part 0 carries the BAM header; later parts are
/// records only). The streaming copy never materialises a part in memory.
fn assemble_bam_parts(parts: &[PathBuf], out: &Path, _force: bool) -> Result<()> {
    let parent = output_dir(out);
    let (tmp_file, tmp_path) = tempfile::Builder::new()
        .prefix(".bz_numa_asm_")
        .suffix(".tmp")
        .tempfile_in(&parent)
        .map_err(|e| anyhow::anyhow!("create assembly temp in {parent:?}: {e}"))?
        .into_parts();
    let mut dst = tmp_file;

    let last = parts.len().saturating_sub(1);
    for (i, part) in parts.iter().enumerate() {
        append_part(&mut dst, part, i != last)
            .map_err(|e| anyhow::anyhow!("assembling part {i} ({part:?}): {e}"))?;
    }
    dst.sync_all()?;
    drop(dst);

    tmp_path
        .persist(out)
        .map_err(|e| anyhow::anyhow!("finalize NUMA-assembled BAM {out:?}: {e}"))?;
    Ok(())
}

/// Stream-copy one part into `dst`, optionally dropping its trailing BGZF EOF marker.
/// When stripping, the marker is first verified so a truncated/corrupt part is a clean
/// error rather than silent corruption of the assembled BAM.
fn append_part(dst: &mut File, part: &Path, strip_eof: bool) -> Result<()> {
    let mut src = File::open(part)?;
    let len = src.metadata()?.len();
    let to_copy = if strip_eof {
        let eof = BGZF_EOF.len() as u64;
        if len < eof {
            anyhow::bail!("part is shorter than a BGZF EOF marker ({len} bytes)");
        }
        src.seek(SeekFrom::End(-(eof as i64)))?;
        let mut tail = [0u8; BGZF_EOF.len()];
        src.read_exact(&mut tail)?;
        if tail != BGZF_EOF {
            anyhow::bail!("part does not end in the expected BGZF EOF marker");
        }
        src.seek(SeekFrom::Start(0))?;
        len - eof
    } else {
        len
    };
    let mut limited = src.take(to_copy);
    std::io::copy(&mut limited, dst)?;
    Ok(())
}

/// The directory holding `path` (its parent, or `.` if none) — used to statvfs the
/// output filesystem and to place sibling temp files for atomic rename.
fn output_dir(path: &Path) -> PathBuf {
    path.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."))
}

// ─────────────────────────── compress sharding ───────────────────────────
//
// Mirror of the decode path: re-exec one pinned worker per chunk range, each
// compressing its disjoint range of the input BAM into a self-contained `.bz`
// part, then assemble. Two differences vs decode: (1) the input BAM has no record
// sync marker, so a worker seeks to its chunk's captured BGZF virtual position
// (from the prescan) rather than to a self-delimiting archive offset; (2) the
// parts are raw `.bz` fragments (part 0 carries the global header with the
// whole-archive totals, later parts are headerless chunk streams), so assembly is
// a PURE concatenation — byte-identical to a single-process `bz compress`.

/// Hidden `_numa-compress-worker` args: compress one chunk range into one `.bz` part.
#[derive(Parser, Debug)]
pub struct NumaCompressWorkerArgs {
    #[arg(long)]
    input: PathBuf,
    #[arg(long)]
    chunk_start: u32,
    #[arg(long)]
    chunk_end: u32,
    #[arg(long)]
    start_vpos: u64,
    #[arg(long)]
    total_records: u64,
    #[arg(long)]
    total_chunks: u32,
    /// Records this worker's chunk range must compress (from the prescan); the
    /// worker bails if it produces a different count (defence against chunk-size drift).
    #[arg(long)]
    expected_records: u64,
    #[arg(long)]
    node: u32,
    #[arg(long)]
    threads: usize,
    #[arg(long)]
    out_part: PathBuf,
    #[arg(long)]
    worker_id: u32,
    /// JSON-serialized `AdvancedOptions` so the worker compresses byte-identically
    /// to the driver's single-process path.
    #[arg(long)]
    advanced_json: String,
}

/// Worker entry (hidden `_numa-compress-worker`): pin to the assigned node BEFORE
/// building the rayon pool (so pool threads fault node-local), then compress the
/// chunk range into the part. Exit codes match the cross-process ABI the driver
/// classifies: 0 success, 1 general error (via `main`), 3 pin failure.
pub fn run_compress_worker(a: &NumaCompressWorkerArgs) -> Result<()> {
    pin::worker_test_hooks(a.worker_id, "compress worker");
    pin::pin_worker_to_node(a.node);
    if a.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(a.threads)
            .build_global();
    }
    let advanced: AdvancedOptions = serde_json::from_str(&a.advanced_json)
        .map_err(|e| anyhow::anyhow!("invalid --advanced-json: {e}"))?;
    // `compress_chunk_range` reads only `config.advanced`; the other fields are
    // placeholders (it takes `input`/`out_part` as explicit arguments).
    let config = CompressConfig {
        input: a.input.clone(),
        output: a.out_part.clone(),
        working_dir: output_dir(&a.out_part),
        advanced,
    };
    bz_lib::compression::compress_chunk_range(
        &a.input,
        a.chunk_start,
        a.chunk_end,
        a.start_vpos,
        a.total_records,
        a.total_chunks,
        a.expected_records,
        &config,
        &a.out_part,
    )?;
    Ok(())
}

/// CLI entry from the `Compress` arm. Decides in-process vs sharded and runs it.
/// `total_threads` is the resolved `-t` budget; `force` gates the final overwrite
/// (checked by the caller before any work).
pub fn compress_with_numa(
    mode: NumaMode,
    config: &CompressConfig,
    total_threads: usize,
    force: bool,
) -> Result<()> {
    if mode == NumaMode::Off || mode == NumaMode::Fixed(1) {
        return bz_lib::compress(config);
    }

    // FIX (data-loss): FREEZE the auto-resolved compression level ONCE, before any
    // prescan or worker spawn. With the default --level 0 (auto), apply_level re-reads
    // /proc/meminfo on every call, so the prescan and each re-exec'd worker — separate
    // processes, at different times, with siblings allocating GBs — could resolve
    // DIFFERENT levels -> different chunk_size -> the worker's chunking would diverge
    // from the prescan's virtual-offset partition (silent record loss, or a hard fail).
    // Pinning the concrete level makes every downstream apply_level an idempotent no-op,
    // so one chunk_size holds across the prescan and ALL workers.
    let mut advanced = config.advanced.clone();
    advanced.level = bz_lib::compression::resolve_compress_level(config);
    let frozen = CompressConfig {
        input: config.input.clone(),
        output: config.output.clone(),
        working_dir: config.working_dir.clone(),
        advanced,
    };
    let config = &frozen;

    // Workload + part-capacity proxy: the source BAM size (parts hold the smaller
    // archive, so requiring this much free space is safe).
    let est_bytes = std::fs::metadata(&config.input).map(|m| m.len()).unwrap_or(0);

    // Cheap precheck: an AUTO compress of a sub-threshold BAM is sent in-process by the
    // gate anyway, so skip the (full-inflate) prescan that would precede it — otherwise
    // a < MIN_WORKLOAD_BYTES BAM on NUMA hardware gets inflated TWICE (prescan + the
    // in-process fallback). Fixed(N) must still prescan (it must shard or bail). The
    // forced-shard test hook is exempt (it forces sharding regardless of size).
    if matches!(mode, NumaMode::Auto)
        && est_bytes < gate::MIN_WORKLOAD_BYTES
        && !runtime::test_hook_set("QZ_NUMA_FORCE_WORKERS")
    {
        return bz_lib::compress(config);
    }

    // TEST HOOK (debug-only): force an N-way shard regardless of topology/workload,
    // children UNPINNED, so the whole spawn → compress → assemble path is
    // exercisable on any box.
    if let Some(n) =
        runtime::test_hook("QZ_NUMA_FORCE_WORKERS").and_then(|s| s.parse::<u32>().ok())
        && n >= 2
    {
        let layout = bz_lib::compression::read_bam_compress_layout(&config.input, config)?;
        if layout.num_chunks >= n {
            let usable = vec![topology::NumaNode {
                id: 0,
                cpus: (0..total_threads.max(1) as u32).collect(),
                mem_available_bytes: u64::MAX,
            }];
            let ranges =
                gate::partition_balanced(&layout.per_chunk_reads, n, &usable, total_threads);
            return run_sharded_compress(&ranges, config, &layout, mode, force);
        }
    }

    let fallback_or_err =
        |msg: &str| runtime::fallback_or_err(mode, msg, || bz_lib::compress(config));

    // CHEAP precheck: a real-asymmetry NUMA topology, before touching the BAM.
    let topo = topology::detect();
    if !topo.as_ref().map(|t| t.has_real_asymmetry).unwrap_or(false) {
        return fallback_or_err("no usable NUMA topology (need >=2 asymmetric nodes on Linux)");
    }

    // Prescan the per-chunk virtual-offset layout (one inflate pass, no compression).
    let layout = match bz_lib::compression::read_bam_compress_layout(&config.input, config) {
        Ok(l) => l,
        Err(e) => {
            return match mode {
                NumaMode::Fixed(_) => Err(e),
                _ => bz_lib::compress(config),
            };
        }
    };

    let worker_peak_rss =
        bz_lib::compression::compress_peak_rss_bound(layout.num_records, layout.num_chunks);

    let out_dir = output_dir(&config.output);
    let inputs = GateInputs {
        archive_type: 0,
        num_chunks: layout.num_chunks,
        per_chunk_reads: layout.per_chunk_reads.clone(),
        est_decoded_bytes: est_bytes,
        total_threads,
        io: IoKind::Files,
        workdir_free_bytes: runtime::free_space(&config.working_dir).unwrap_or(0),
        worker_peak_rss_bytes: worker_peak_rss,
        allowed_cpus: topology::allowed_cpus(),
        allowed_mem_nodes: topology::allowed_mem_nodes(),
        increment_b_ready: true,
        shardable: true,
        // `.bz` parts have a variable size per chunk: always part-file + assembly.
        direct_write: false,
        output_fs_free_bytes: runtime::free_space(&out_dir).unwrap_or(0),
        workdir_output_same_fs: runtime::same_device(&config.working_dir, &out_dir),
    };

    match gate::decide(mode, topo.as_ref(), &inputs)? {
        // KNOWN COST: when the gate routes an above-threshold, real-topology input to
        // in-process (e.g. too few chunks to shard usefully, or a tight per-node memory
        // budget), the BAM is inflated twice — once by the prescan above to obtain the
        // chunk layout the gate requires, then again here by `bz_lib::compress`. The cheap
        // short-circuits (sub-threshold workload, no NUMA asymmetry) already skip the
        // prescan; the remaining gate triggers genuinely need the real per-chunk layout,
        // so eliminating this would require a non-inflating BGZF-boundary prescan (future
        // work). It is bounded to this fallback path; the sharded path inflates once.
        Plan::InProcess => bz_lib::compress(config),
        Plan::Sharded { ranges } => run_sharded_compress(&ranges, config, &layout, mode, force),
    }
}

/// Spawn one pinned worker per range, wait, then concatenate the `.bz` parts. On a
/// pin failure (`exit 3`) under `auto`, fall back to a correct in-process compress.
fn run_sharded_compress(
    ranges: &[WorkerRange],
    config: &CompressConfig,
    layout: &bz_lib::compression::BzCompressLayout,
    mode: NumaMode,
    force: bool,
) -> Result<()> {
    let exe = std::env::current_exe()
        .map_err(|e| anyhow::anyhow!("cannot locate self for re-exec: {e}"))?;
    let advanced_json = serde_json::to_string(&config.advanced)
        .map_err(|e| anyhow::anyhow!("serialize advanced options: {e}"))?;

    let tmp = tempfile::Builder::new()
        .prefix(".bz_numa_c_")
        .tempdir_in(&config.working_dir)
        .map_err(|e| anyhow::anyhow!("create NUMA part dir in {:?}: {e}", config.working_dir))?;
    let parts: Vec<PathBuf> = (0..ranges.len())
        .map(|wi| tmp.path().join(format!("w{wi}.part")))
        .collect();

    let mut guard = ShardGuard::with_tmp(tmp);
    for (wi, r) in ranges.iter().enumerate() {
        if runtime::forced_spawn_fail(wi) {
            return Err(anyhow::anyhow!("test hook: forced spawn failure for worker {wi}"));
        }
        let start_vpos = layout.chunk_start_vpos[r.chunk_start as usize];
        // Records this worker must compress = the prescan's count for its chunk range.
        // The worker bails if it produces a different number (defence against any
        // chunk-size drift, so a divergence can never silently drop records).
        let expected_records: u64 = layout.per_chunk_reads
            [r.chunk_start as usize..r.chunk_end as usize]
            .iter()
            .sum();
        let mut cmd = compress_worker_command(
            &exe,
            config,
            r,
            wi,
            &parts[wi],
            start_vpos,
            layout.num_records,
            layout.num_chunks,
            expected_records,
            &advanced_json,
        );
        match cmd.spawn() {
            Ok(child) => guard.children.push(child),
            Err(e) => return Err(anyhow::anyhow!("spawn worker {wi}: {e}")),
        }
    }

    runtime::write_spawn_log(guard.children.len());

    let fail_codes = runtime::wait_for_workers(&mut guard.children)?;
    if !fail_codes.is_empty() {
        let pin_failure = fail_codes.contains(&Some(3));
        drop(guard);
        if pin_failure && mode == NumaMode::Auto {
            return bz_lib::compress(config);
        }
        anyhow::bail!("a NUMA compress worker failed; aborted (no output written)");
    }

    let result = assemble_bz_parts(&parts, &config.output, force);
    guard.armed = false;
    drop(guard);
    result
}

/// Build the compress worker re-exec `Command`.
#[allow(clippy::too_many_arguments)]
fn compress_worker_command(
    exe: &Path,
    config: &CompressConfig,
    r: &WorkerRange,
    wi: usize,
    out_part: &Path,
    start_vpos: u64,
    total_records: u64,
    total_chunks: u32,
    expected_records: u64,
    advanced_json: &str,
) -> Command {
    let mut cmd = Command::new(exe);
    cmd.arg("_numa-compress-worker")
        .arg("--input")
        .arg(&config.input)
        .arg("--chunk-start")
        .arg(r.chunk_start.to_string())
        .arg("--chunk-end")
        .arg(r.chunk_end.to_string())
        .arg("--start-vpos")
        .arg(start_vpos.to_string())
        .arg("--total-records")
        .arg(total_records.to_string())
        .arg("--total-chunks")
        .arg(total_chunks.to_string())
        .arg("--expected-records")
        .arg(expected_records.to_string())
        .arg("--node")
        .arg(r.node.to_string())
        .arg("--threads")
        .arg(r.threads.to_string())
        .arg("--out-part")
        .arg(out_part)
        .arg("--worker-id")
        .arg(wi.to_string())
        .arg("--advanced-json")
        .arg(advanced_json);
    runtime::apply_worker_env(&mut cmd, &["QZ_NUMA_FAIL", "QZ_NUMA_SLOW", "QZ_NUMA_PINFAIL"]);
    cmd
}

/// Concatenate the `.bz` parts (in worker/chunk order) into the final archive via a
/// sibling temp + atomic rename. Compress parts are raw archive fragments — part 0
/// carries the global header (with whole-archive totals), later parts are headerless
/// chunk streams — so assembly is a PURE byte concatenation (no EOF stripping). The
/// result is byte-identical to a single-process `bz compress`.
fn assemble_bz_parts(parts: &[PathBuf], out: &Path, _force: bool) -> Result<()> {
    let parent = output_dir(out);
    let (tmp_file, tmp_path) = tempfile::Builder::new()
        .prefix(".bz_numa_casm_")
        .suffix(".tmp")
        .tempfile_in(&parent)
        .map_err(|e| anyhow::anyhow!("create assembly temp in {parent:?}: {e}"))?
        .into_parts();
    let mut dst = tmp_file;

    for (i, part) in parts.iter().enumerate() {
        let mut src =
            File::open(part).map_err(|e| anyhow::anyhow!("open part {i} ({part:?}): {e}"))?;
        std::io::copy(&mut src, &mut dst)
            .map_err(|e| anyhow::anyhow!("copy part {i} ({part:?}): {e}"))?;
    }
    dst.sync_all()?;
    drop(dst);

    tmp_path
        .persist(out)
        .map_err(|e| anyhow::anyhow!("finalize NUMA-assembled archive {out:?}: {e}"))?;
    Ok(())
}
