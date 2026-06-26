//! NUMA decode driver: gate -> spawn pinned re-exec workers -> assemble -> cleanup.

use crate::numa::gate::{self, GateInputs, IoKind, Plan, WorkerRange};
use crate::numa::runtime::{self, ShardGuard};
use crate::numa::{topology, NumaMode};
use anyhow::Result;
use qz_lib::cli::DecompressConfig;
use std::path::PathBuf;
use std::process::Command;

/// Flipped to true once Increment B (reference decode_chunk_range) has landed.
const INCREMENT_B_READY: bool = true;

/// Entry from the CLI Decompress arm. Decides in-process vs sharded, runs it.
pub fn decompress_with_numa(mode: NumaMode, config: &DecompressConfig) -> Result<()> {
    // 0) off and Fixed(1) -> straight to qz-lib.
    if mode == NumaMode::Off || mode == NumaMode::Fixed(1) {
        return qz_lib::compression::decompress(config);
    }

    // TEST HOOK (debug-only): force an N-way shard regardless of topology/workload,
    // children UNPINNED, so spawn/decode/assemble is exercisable on a non-NUMA box.
    if let Some(n) = runtime::test_hook("QZ_NUMA_FORCE_WORKERS").and_then(|s| s.parse::<u32>().ok())
        && n >= 2
    {
        let io = io_kind(config);
        if io == IoKind::Files {
            let layout = qz_lib::compression::read_chunk_layout(&config.input)?;
            if config.gzipped && layout.archive_type != 0 {
                return qz_lib::compression::decompress(config);
            }
            if layout.shardable
                && (!matches!(layout.archive_type, 2 | 4) || INCREMENT_B_READY)
                && layout.num_chunks >= n
            {
                let usable = vec![crate::numa::topology::NumaNode {
                    id: 0, cpus: (0..config.num_threads as u32).collect(), mem_available_bytes: u64::MAX,
                }];
                let ranges = crate::numa::gate::partition_balanced(
                    &layout.per_chunk_reads, n, &usable, config.num_threads,
                );
                return run_sharded(&ranges, config, layout.archive_type, mode);
            }
        }
    }

    // Mode-sensitive bail: auto falls back in-process, Fixed errors.
    let fallback_or_err =
        |msg: &str| runtime::fallback_or_err(mode, msg, || qz_lib::compression::decompress(config));

    // CHEAP prechecks BEFORE reading the archive footer.
    let io = io_kind(config);
    if io == IoKind::Stdio {
        return fallback_or_err("stdin/stdout cannot be sharded");
    }
    let topo = topology::detect();
    if !topo.as_ref().map(|t| t.has_real_asymmetry).unwrap_or(false) {
        return fallback_or_err("no usable NUMA topology (need >=2 asymmetric nodes on Linux)");
    }

    // On NUMA hardware: read the footer.
    let layout = match qz_lib::compression::read_chunk_layout(&config.input) {
        Ok(l) => l,
        Err(e) => return match mode {
            NumaMode::Fixed(_) => Err(e),
            _ => qz_lib::compression::decompress(config),
        },
    };

    // Cluster archives (archive_type 3) are not chunk-range shardable; decode them
    // in-process for EVERY mode, including Fixed. Unlike the mixed-codec legacy paired
    // non-shardable case (a degenerate edge that errors under Fixed via the gate),
    // `--cluster` is a first-class user mode, so an explicit `--numa N` on a cluster
    // archive must decode correctly rather than fail. Emit a one-line note under Fixed
    // so a user who asked for sharding knows why it did not apply.
    if layout.archive_type == 3 {
        if matches!(mode, NumaMode::Fixed(_)) {
            eprintln!(
                "note: cluster (--cluster) archives are not NUMA-shardable; decoding in-process"
            );
        }
        return qz_lib::compression::decompress(config);
    }

    // gzip is single-end-only; for paired/reference route to in-process.
    if config.gzipped && layout.archive_type != 0 {
        return qz_lib::compression::decompress(config);
    }

    let worker_peak_rss = qz_lib::compression::decode_peak_rss_bound(
        layout.encoding_type,
        layout.archive_type,
        layout.reference_resident_bytes,
    );

    // Non-gzip + a present ChunkDecodedSizes table ⇒ direct-write: every worker seeks into
    // its pre-sized shared output file(s) at the per-mate byte offset (no part files /
    // assembly temp / concat pass). Covers single-end (0, 1 output), paired (1, R1/R2),
    // single-end reference (4, 1 output), and paired reference (2, R1/R2) — all four emit
    // the per-mate table. The table-present check naturally excludes any archive without it.
    let direct_write = matches!(layout.archive_type, 0 | 1 | 2 | 4)
        && !config.gzipped
        && !layout.per_chunk_decoded_bytes.is_empty();

    // est_decoded_bytes: for direct-write the table gives the EXACT total, so use a
    // checked sum (a corrupt-table overflow is a clean fallback). Otherwise the rough
    // compressed×8 heuristic feeds the workload threshold + part-file capacity check.
    let est_decoded_bytes = if direct_write {
        match layout.per_chunk_decoded_bytes.iter().try_fold(0u64, |a, &b| a.checked_add(b)) {
            Some(t) => t,
            None => return fallback_or_err("decoded offset overflow (corrupt size table?)"),
        }
    } else {
        estimate_decoded_bytes(&config.input)
    };

    // The final output may not exist yet, so statvfs its PARENT dir.
    let out_dir = config.output[0]
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .map(std::path::Path::to_path_buf)
        .unwrap_or_else(|| std::path::PathBuf::from("."));

    let inputs = GateInputs {
        archive_type: layout.archive_type,
        num_chunks: layout.num_chunks,
        per_chunk_reads: layout.per_chunk_reads.clone(),
        est_decoded_bytes,
        total_threads: config.num_threads,
        io,
        workdir_free_bytes: runtime::free_space(&config.working_dir).unwrap_or(0),
        worker_peak_rss_bytes: worker_peak_rss,
        allowed_cpus: topology::allowed_cpus(),
        allowed_mem_nodes: topology::allowed_mem_nodes(),
        increment_b_ready: INCREMENT_B_READY,
        shardable: layout.shardable,
        direct_write,
        output_fs_free_bytes: runtime::free_space(&out_dir).unwrap_or(0),
        workdir_output_same_fs: runtime::same_device(&config.working_dir, &out_dir),
    };

    match gate::decide(mode, topo.as_ref(), &inputs)? {
        Plan::InProcess => qz_lib::compression::decompress(config),
        Plan::Sharded { ranges } => run_sharded(&ranges, config, layout.archive_type, mode),
    }
}

/// Build the worker re-exec `Command` shared by BOTH the part-file (`run_sharded`) and
/// direct-write (`run_sharded_direct`) paths, so their arg/env sets can't drift.
///
/// Common to both: the `_numa-decode-worker` subcommand + `--input`/`--chunk-start`/
/// `--chunk-end`/`--node`/`--threads`/`--gzip-level`/`--worker-id`, `QZ_NO_BANNER=1`, the
/// `QZ_NUMA_NO_PIN` env when the forced-shard hook is set, and the full debug-only
/// test-hook env pass-through. Path-specific inputs are passed in:
/// * `out_parts` — N distinct part files (part-file path) or N shared pre-sized temps
///   (direct; single-end: 1, paired: 2).
/// * `regions` — one `(base, len)` per output (direct only); EMPTY for the part-file path.
///   Each adds an `--out-base-offset`/`--out-region-len` pair, positionally matched to
///   `--out-part` by the worker.
///
/// `--gzipped` is added iff `config.gzipped` (direct's config.gzipped is always false, so
/// it's naturally never gzipped). The caller keeps the `QZ_NUMA_SPAWN_FAIL` pre-spawn bail
/// and the `cmd.spawn()` + guard bookkeeping outside this builder.
fn worker_command(
    exe: &std::path::Path,
    config: &DecompressConfig,
    r: &WorkerRange,
    wi: usize,
    out_parts: &[std::path::PathBuf],
    regions: &[(u64, u64)],
) -> Command {
    let mut cmd = Command::new(exe);
    cmd.arg("_numa-decode-worker")
        .arg("--input").arg(&config.input)
        .arg("--chunk-start").arg(r.chunk_start.to_string())
        .arg("--chunk-end").arg(r.chunk_end.to_string())
        .arg("--node").arg(r.node.to_string())
        .arg("--threads").arg(r.threads.to_string())
        .arg("--gzip-level").arg(config.gzip_level.to_string())
        .arg("--worker-id").arg(wi.to_string());
    if config.gzipped {
        cmd.arg("--gzipped");
    }
    for p in out_parts {
        cmd.arg("--out-part").arg(p);
    }
    for (base, len) in regions {
        cmd.arg("--out-base-offset").arg(base.to_string())
            .arg("--out-region-len").arg(len.to_string());
    }
    // QZ_NO_BANNER + forced-hook QZ_NUMA_NO_PIN + the debug-only test-hook pass-through.
    // QZ_NUMA_FAKE_REGION is inert on the part-file path (region=None) but passing it
    // everywhere keeps the two spawn loops from diverging.
    runtime::apply_worker_env(
        &mut cmd,
        &["QZ_NUMA_FAIL", "QZ_NUMA_SLOW", "QZ_NUMA_PINFAIL", "QZ_NUMA_FAKE_REGION"],
    );
    cmd
}

fn run_sharded(
    ranges: &[WorkerRange],
    config: &DecompressConfig,
    archive_type: u8,
    mode: NumaMode,
) -> Result<()> {
    let exe = std::env::current_exe()
        .map_err(|e| anyhow::anyhow!("cannot locate self for re-exec: {e}"))?;
    let final_out = &config.output[0];
    // Single-output modes: single-end default (0) AND single-end reference (4) each
    // write ONE bare output[0] file. Paired (1) and (two-output) reference (2) write
    // two suffixed files.
    let n_out = if matches!(archive_type, 0 | 4) { 1 } else { 2 };

    // Preflight: refuse pre-existing output BEFORE spawning workers (match the
    // in-process path, which checks at entry). Otherwise a known refusal wastes a
    // full sharded decode. (F5)
    if !config.force {
        let outs: Vec<std::path::PathBuf> = if matches!(archive_type, 0 | 4) {
            // Single-end (default OR reference) writes ONE bare output[0] file.
            vec![final_out.clone()]
        } else {
            vec![
                crate::numa::assemble::with_suffix(final_out, "_R1.fastq"),
                crate::numa::assemble::with_suffix(final_out, "_R2.fastq"),
            ]
        };
        for o in &outs {
            if o.exists() {
                anyhow::bail!("Output file already exists: {}\nUse --force to overwrite", o.display());
            }
        }
    }

    // Non-gzip + a present ChunkDecodedSizes table ⇒ direct-write: workers write disjoint
    // regions into pre-sized shared output file(s) (no part files, no concat). Covers
    // single-end (0), paired (1), and reference (2, 4). The forced-shard hook reaches
    // run_sharded directly (bypassing the gate), so detecting it here covers that path too.
    let direct = matches!(archive_type, 0 | 1 | 2 | 4)
        && !config.gzipped
        && qz_lib::compression::read_chunk_layout(&config.input)
            .map(|l| !l.per_chunk_decoded_bytes.is_empty())
            .unwrap_or(false);
    if direct {
        return run_sharded_direct(ranges, config, mode, archive_type, &exe);
    }

    let tmp = tempfile::Builder::new()
        .prefix(".qz_numa_")
        .tempdir_in(&config.working_dir)
        .map_err(|e| anyhow::anyhow!("create NUMA part dir in {:?}: {e}", config.working_dir))?;
    let parts: Vec<Vec<PathBuf>> = (0..ranges.len())
        .map(|wi| (0..n_out).map(|oi| tmp.path().join(format!("w{wi}.o{oi}.part"))).collect())
        .collect();

    let mut guard = ShardGuard::with_tmp(tmp);
    for (wi, r) in ranges.iter().enumerate() {
        if runtime::forced_spawn_fail(wi) {
            return Err(anyhow::anyhow!("test hook: forced spawn failure for worker {wi}"));
        }
        let mut cmd = worker_command(&exe, config, r, wi, &parts[wi], &[]);
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
        // Only a pin failure (exit 3) is a benign, recoverable condition that falls back
        // to in-process under `auto`. A general/unexpected worker failure aborts loudly
        // (kills the survivor, writes no partial output) rather than silently masking a
        // bug or corruption behind a much slower in-process retry — see
        // `worker_failure_kills_survivor_and_aborts_cleanly`.
        if pin_failure && mode == NumaMode::Auto {
            return qz_lib::compression::decompress(config);
        }
        anyhow::bail!("a NUMA decode worker failed; aborted (no output written)");
    }

    let result = if matches!(archive_type, 0 | 4) {
        // Single-end (default OR reference): ONE part per worker, concatenated
        // into the verbatim final_out (output[0]) — no suffix.
        let ordered: Vec<PathBuf> = parts.iter().map(|ps| ps[0].clone()).collect();
        crate::numa::assemble::assemble_single(&ordered, final_out, config.force)
    } else {
        let r1: Vec<PathBuf> = parts.iter().map(|ps| ps[0].clone()).collect();
        let r2: Vec<PathBuf> = parts.iter().map(|ps| ps[1].clone()).collect();
        crate::numa::assemble::assemble_paired(&r1, &r2, final_out, config.force)
    };
    guard.armed = false;
    drop(guard);
    result
}

/// Direct-write (single-end OR paired): pre-size N shared output files (single-end: 1;
/// paired: 2, R1/R2), spawn workers that seek to each mate's `prefix_m[chunk_start]` and
/// write `prefix_m[chunk_end]-prefix_m[chunk_start]` bytes (each verifying base/len
/// against that mate's table), then publish (rename for 1; backup-and-rollback for 2). No
/// part files, no concat. On any worker failure: kill survivors, remove the temps (no
/// output published), and under Auto fall back to a correct in-process decode.
fn run_sharded_direct(
    ranges: &[WorkerRange],
    config: &DecompressConfig,
    mode: NumaMode,
    archive_type: u8,
    exe: &std::path::Path,
) -> Result<()> {
    use crate::numa::assemble::with_suffix;
    let final_out = &config.output[0];

    // Mode-sensitive bail mirroring decompress_with_numa's: auto → correct in-process
    // decode (ignores the table), fixed → hard error.
    let fallback_or_err =
        |msg: &str| runtime::fallback_or_err(mode, msg, || qz_lib::compression::decompress(config));

    // The final output path(s): single-output modes — single-end (0) AND single-end
    // reference (4) — write ONE bare output[0]; two-output modes — paired (1) and paired
    // reference (2) — write two suffixed _R1/_R2 files (same naming the part-file assembler
    // uses).
    let finals: Vec<PathBuf> = if matches!(archive_type, 0 | 4) {
        vec![final_out.clone()]
    } else {
        vec![with_suffix(final_out, "_R1.fastq"), with_suffix(final_out, "_R2.fastq")]
    };

    // 1) Authoritative per-MATE decoded sizes → CHECKED per-mate prefix sums. The driver
    //    needs the per-mate breakdown (each mate has its own output): prefix_m[c] is the
    //    byte where chunk c of mate m begins.
    let layout = qz_lib::compression::read_chunk_layout(&config.input)?;
    let per_mate = &layout.decoded_sizes_per_mate;
    if per_mate.len() != finals.len() {
        return fallback_or_err("ChunkDecodedSizes mate count != output count (corrupt table?)");
    }
    let mut prefixes: Vec<Vec<u64>> = Vec::with_capacity(per_mate.len());
    let mut totals: Vec<u64> = Vec::with_capacity(per_mate.len());
    for sizes in per_mate {
        let mut prefix = vec![0u64; sizes.len() + 1];
        for i in 0..sizes.len() {
            prefix[i + 1] = match prefix[i].checked_add(sizes[i]) {
                Some(v) => v,
                None => return fallback_or_err("decoded offset overflow (corrupt size table?)"),
            };
        }
        totals.push(prefix[sizes.len()]);
        prefixes.push(prefix);
    }

    // 2) Sanity cap on the TOTAL output size across all mates (corrupt-table guard).
    let grand_total: u64 = match totals.iter().try_fold(0u64, |a, &b| a.checked_add(b)) {
        Some(t) => t,
        None => return fallback_or_err("decoded size overflow (corrupt table?)"),
    };
    if grand_total > crate::numa::MAX_DIRECT_OUTPUT_BYTES {
        return fallback_or_err("decoded size implausibly large (corrupt table?)");
    }

    // 3) Output-fs capacity: statvfs the EXISTING parent dir (the outputs don't exist yet;
    //    every mate shares the parent of output[0]).
    let out_dir = final_out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| std::path::Path::new("."));
    if runtime::free_space(out_dir).unwrap_or(0) < grand_total {
        return fallback_or_err("output filesystem lacks space for direct-write");
    }

    // 4) Pre-size the N shared temps; their guard removes all of them on any early return.
    let (tmps, tmp_guard) = crate::numa::assemble::make_direct_temps(&finals, &totals)?;

    // 5) Spawn workers — same shared temp paths + this shard's per-mate [offset, offset+len)
    //    regions. The child guard owns NO tempdir here (the pre-sized outputs are owned by
    //    `tmp_guard`); it only reaps children.
    let mut guard = ShardGuard::children_only();
    for (wi, r) in ranges.iter().enumerate() {
        if runtime::forced_spawn_fail(wi) {
            return Err(anyhow::anyhow!("test hook: forced spawn failure for worker {wi}"));
        }
        let regions: Vec<(u64, u64)> = prefixes
            .iter()
            .map(|prefix| {
                let base = prefix[r.chunk_start as usize];
                let len = prefix[r.chunk_end as usize] - prefix[r.chunk_start as usize];
                (base, len)
            })
            .collect();
        // NB: config.gzipped is false here ⇒ worker_command never adds `--gzipped`.
        let mut cmd = worker_command(exe, config, r, wi, &tmps, &regions);
        match cmd.spawn() {
            Ok(child) => guard.children.push(child),
            Err(e) => return Err(anyhow::anyhow!("spawn worker {wi}: {e}")),
        }
    }

    runtime::write_spawn_log(guard.children.len());

    // 6) Poll to completion (kills survivors on the first failure), then classify BOTH
    //    exit 3 (pin) and exit 4 (integrity).
    let fail_codes = runtime::wait_for_workers(&mut guard.children)?;

    // 7) Failure: drop the child guard (kills survivors) and the tmp guard (removes the
    //    pre-sized temps ⇒ no output). Under Auto a pin/integrity failure falls back to a
    //    correct in-process decode (which ignores the table); anything else is hard.
    if !fail_codes.is_empty() {
        let pin_failure = fail_codes.contains(&Some(3));
        let integrity_failure = fail_codes.contains(&Some(4));
        drop(guard);
        drop(tmp_guard);
        // Pin (3) and integrity/region (4) are benign, recoverable conditions that fall
        // back to the in-process decode (which ignores the size table) under `auto`. A
        // general/unexpected failure aborts loudly rather than masking it behind a slow
        // in-process retry — see `worker_failure_kills_survivor_and_aborts_cleanly`.
        if (pin_failure || integrity_failure) && mode == NumaMode::Auto {
            return qz_lib::compression::decompress(config);
        }
        anyhow::bail!("a NUMA direct-write worker failed; aborted (no output written)");
    }

    // 8) Success: disarm the child guard, publish (rename / two-file backup-rollback).
    guard.armed = false;
    drop(guard);
    crate::numa::assemble::publish_direct_multi(&tmps, &finals, tmp_guard)
}

/// stdio classification (shared by the real precheck and the forced-shard hook).
fn io_kind(config: &DecompressConfig) -> IoKind {
    if config.output.iter().any(|o| qz_lib::cli::is_stdio_path(o))
        || qz_lib::cli::is_stdio_path(&config.input)
    { IoKind::Stdio } else { IoKind::Files }
}

/// Rough decoded-size estimate from the compressed archive size and a fixed ratio.
fn estimate_decoded_bytes(archive: &std::path::Path) -> u64 {
    let compressed = std::fs::metadata(archive).map(|m| m.len()).unwrap_or(0);
    // ~8x is the measured qz FASTQ ratio (CLAUDE.md). Used for the gate's workload
    // threshold AND (×2) the working-dir capacity check; the 2× capacity margin
    // absorbs ratio variance on highly-compressible (binned-quality) inputs.
    compressed.saturating_mul(8)
}
