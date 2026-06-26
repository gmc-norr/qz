//! Shared NUMA driver runtime: the spawn / poll / cleanup / test-hook / capacity
//! scaffolding common to the decode driver (`driver.rs`) and the compress driver
//! (`compress_driver.rs`). These were previously copy-pasted in both files; centralizing
//! them keeps the two paths from drifting (a fix to e.g. `free_space` overflow handling
//! now lives in one place) and unifies the old `ShardGuard` + `ChildKiller` guards.

use crate::NumaMode;
use std::path::Path;
use std::process::{Child, Command};

/// Read a test-hook env var, but ONLY in debug builds. A release binary returns `None`
/// unconditionally, so it never honors test env vars like `QZ_NUMA_FORCE_WORKERS`.
pub fn test_hook(name: &str) -> Option<String> {
    if cfg!(debug_assertions) { std::env::var(name).ok() } else { None }
}
pub fn test_hook_set(name: &str) -> bool { test_hook(name).is_some() }

/// The single definition of the fixed-vs-auto fallback policy shared by both drivers and
/// their sub-paths: `auto` runs the supplied in-process fallback, `fixed --numa N` errors.
pub fn fallback_or_err<F>(mode: NumaMode, msg: &str, in_process: F) -> anyhow::Result<()>
where
    F: FnOnce() -> anyhow::Result<()>,
{
    match mode {
        NumaMode::Fixed(_) => anyhow::bail!("--numa: {msg}"),
        _ => in_process(),
    }
}

/// RAII guard for spawned worker children + (optionally) their part-file tempdir. On drop
/// while `armed`, kill+wait every still-running child; the owned `TempDir` (when present)
/// removes all part files. Disarm (`armed = false`) only after the output is published.
///
/// `tmp = None` covers the direct-write path, where the N shared pre-sized outputs are
/// owned by a separate `MultiTmpGuard` (from `assemble::make_direct_temps`) and this guard
/// only reaps children. This unifies the old part-file `ShardGuard` and direct `ChildKiller`.
pub struct ShardGuard {
    pub children: Vec<Child>,
    _tmp: Option<tempfile::TempDir>,
    pub armed: bool,
}
impl ShardGuard {
    /// Part-file path: owns the tempdir holding all part files.
    pub fn with_tmp(tmp: tempfile::TempDir) -> Self {
        Self { children: Vec::new(), _tmp: Some(tmp), armed: true }
    }
    /// Direct-write path: reaps children only (the shared output temp is owned elsewhere).
    pub fn children_only() -> Self {
        Self { children: Vec::new(), _tmp: None, armed: true }
    }
}
impl Drop for ShardGuard {
    fn drop(&mut self) {
        if self.armed {
            for c in &mut self.children {
                let _ = c.kill();
                let _ = c.wait();
            }
        }
    }
}

/// Poll worker children to completion. On the FIRST poll pass that observes any non-success
/// exit, kill+wait all survivors and return the failing exit codes (for pin = 3 / integrity
/// = 4 classification by the caller); an empty `Vec` means every worker succeeded. `try_wait`
/// IO errors propagate. Replaces the three copy-pasted poll loops.
pub fn wait_for_workers(children: &mut [Child]) -> std::io::Result<Vec<Option<i32>>> {
    loop {
        let mut all_done = true;
        let mut fail_codes: Vec<Option<i32>> = Vec::new();
        for c in children.iter_mut() {
            match c.try_wait()? {
                Some(status) if !status.success() => fail_codes.push(status.code()),
                Some(_) => {}
                None => all_done = false,
            }
        }
        if !fail_codes.is_empty() {
            for c in children.iter_mut() {
                let _ = c.kill();
                let _ = c.wait();
            }
            return Ok(fail_codes);
        }
        if all_done {
            return Ok(Vec::new());
        }
        std::thread::sleep(std::time::Duration::from_millis(5));
    }
}

/// Apply the worker-process env shared by both drivers: `QZ_NO_BANNER`, the forced-hook
/// `QZ_NUMA_NO_PIN`, and the debug-only test-hook pass-through (`extra_hooks`). The
/// per-driver args (subcommand, ranges, output) are built by the caller.
pub fn apply_worker_env(cmd: &mut Command, extra_hooks: &[&str]) {
    cmd.env("QZ_NO_BANNER", "1");
    if test_hook_set("QZ_NUMA_FORCE_WORKERS") {
        cmd.env("QZ_NUMA_NO_PIN", "1");
    }
    for k in extra_hooks {
        if let Some(v) = test_hook(k) {
            cmd.env(k, v);
        }
    }
}

/// Debug-only spawn-log hook: write the spawned-child count to the path in
/// `QZ_NUMA_SPAWN_LOG` (assertable by forced-shard tests). No-op in release.
pub fn write_spawn_log(n: usize) {
    if let Some(path) = test_hook("QZ_NUMA_SPAWN_LOG") {
        let _ = std::fs::write(&path, n.to_string());
    }
}

/// Debug-only spawn-fail hook: true iff `QZ_NUMA_SPAWN_FAIL == wi`, i.e. the test asked
/// worker `wi` to fail before it is spawned. Always false in release.
pub fn forced_spawn_fail(wi: usize) -> bool {
    test_hook("QZ_NUMA_SPAWN_FAIL").as_deref() == Some(&wi.to_string())
}

#[cfg(target_os = "linux")]
pub fn free_space(dir: &Path) -> Option<u64> {
    use std::ffi::CString;
    use std::os::unix::ffi::OsStrExt;
    let c = CString::new(dir.as_os_str().as_bytes()).ok()?;
    let mut st: libc::statvfs = unsafe { std::mem::zeroed() };
    if unsafe { libc::statvfs(c.as_ptr(), &mut st) } != 0 {
        return None;
    }
    // saturating_mul → u64::MAX on the (only-at->16-EiB-free) overflow, i.e. report
    // "plenty" rather than a wrapped small value that would falsely deny.
    Some(st.f_bavail.saturating_mul(st.f_frsize))
}
#[cfg(not(target_os = "linux"))]
pub fn free_space(_dir: &Path) -> Option<u64> { None }

/// True iff `a` and `b` reside on the SAME filesystem (same device id). Returns false if
/// either path's metadata can't be read — the cross-fs branch is the conservative default
/// (it never under-provisions capacity vs the same-fs 2×D branch).
#[cfg(unix)]
pub fn same_device(a: &Path, b: &Path) -> bool {
    use std::os::unix::fs::MetadataExt;
    match (std::fs::metadata(a), std::fs::metadata(b)) {
        (Ok(ma), Ok(mb)) => ma.dev() == mb.dev(),
        _ => false,
    }
}
#[cfg(not(unix))]
pub fn same_device(_a: &Path, _b: &Path) -> bool { false }
