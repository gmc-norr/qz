//! Per-thread NUMA binding (Linux). Binds CPU affinity (sched_setaffinity) and
//! memory policy (set_mempolicy MPOL_BIND via raw syscall — libc has no safe wrapper
//! on all targets). MUST be called on the worker's MAIN thread BEFORE rayon builds
//! its pool, so pool threads inherit the affinity and fault buffers node-local.

use crate::topology::NumaNode;

/// Single-word node mask (supports up to 64 nodes — far beyond any real box).
pub fn node_to_mask(node: u32) -> u64 {
    debug_assert!(node < 64);
    1u64 << node
}

#[cfg(target_os = "linux")]
pub fn bind_current_thread_to_node(node: &NumaNode) -> anyhow::Result<()> {
    use std::mem;
    // Defensive: the single-word node mask only addresses nodes 0..64. A node id
    // >= 64 (never seen on real hardware) would silently bind the wrong node via
    // `1 << (id % 64)` in release, so fail closed instead.
    if node.id >= 64 {
        anyhow::bail!("NUMA node id {} exceeds the 64-node mask", node.id);
    }
    // 1) CPU affinity to the node's CPUs.
    unsafe {
        let mut set: libc::cpu_set_t = mem::zeroed();
        libc::CPU_ZERO(&mut set);
        // CPU_SET indexes a fixed-width cpu_set_t (1024 bits on glibc); a cpu id at or
        // beyond that width is an out-of-bounds write (UB). Skip such ids — only
        // reachable on >1024-logical-CPU hardware — rather than risk memory corruption;
        // the node's in-range CPUs still pin affinity correctly.
        let cpuset_bits = 8 * mem::size_of::<libc::cpu_set_t>();
        for &cpu in &node.cpus {
            if (cpu as usize) < cpuset_bits {
                libc::CPU_SET(cpu as usize, &mut set);
            }
        }
        if libc::sched_setaffinity(0, mem::size_of::<libc::cpu_set_t>(), &set) != 0 {
            return Err(std::io::Error::last_os_error()).map_err(|e| {
                anyhow::anyhow!("sched_setaffinity(node {}) failed: {e}", node.id)
            });
        }
    }
    // 2) Memory policy MPOL_BIND to this node (raw syscall).
    const MPOL_BIND: libc::c_int = 2;
    let mask = node_to_mask(node.id);
    let maxnode: libc::c_ulong = 64; // bits in `mask`
    let ret = unsafe {
        libc::syscall(
            libc::SYS_set_mempolicy,
            MPOL_BIND,
            &mask as *const u64,
            maxnode,
        )
    };
    if ret != 0 {
        return Err(std::io::Error::last_os_error())
            .map_err(|e| anyhow::anyhow!("set_mempolicy(MPOL_BIND, node {}) failed: {e}", node.id));
    }
    Ok(())
}

#[cfg(not(target_os = "linux"))]
pub fn bind_current_thread_to_node(_node: &NumaNode) -> anyhow::Result<()> {
    anyhow::bail!("NUMA binding is only supported on Linux")
}

/// Debug-only worker test hooks, shared by the decode and compress worker entry points
/// (no-ops in release): `QZ_NUMA_SLOW` sleeps, `QZ_NUMA_FAIL` exits 1, `QZ_NUMA_PINFAIL`
/// exits 3. `label` distinguishes decode vs compress in the `QZ_NUMA_FAIL` diagnostic.
pub fn worker_test_hooks(worker_id: u32, label: &str) {
    if !cfg!(debug_assertions) {
        return;
    }
    let me = worker_id.to_string();
    if std::env::var("QZ_NUMA_SLOW").ok().as_deref() == Some(&me) {
        std::thread::sleep(std::time::Duration::from_secs(30));
    }
    if std::env::var("QZ_NUMA_FAIL").ok().as_deref() == Some(&me) {
        eprintln!("QZ_NUMA_TEST: forced failure of {label} {me}");
        std::process::exit(1);
    }
    if std::env::var("QZ_NUMA_PINFAIL").ok().as_deref() == Some(&me)
        || std::env::var("QZ_NUMA_PINFAIL").ok().as_deref() == Some("all")
    {
        eprintln!("QZ_NUMA_PIN_FAILED: forced");
        std::process::exit(3);
    }
}

/// Pin the current (worker main) thread to NUMA `node` BEFORE rayon builds its pool, so the
/// pool threads inherit the affinity and fault node-local. The node's CPUs are intersected
/// with the inherited cpuset so pinning never widens the user's granted affinity (F4). On
/// any failure this prints a `QZ_NUMA_PIN_FAILED:` line and exits 3 (the parent classifies
/// exit 3 as a pin failure → `auto` falls back to in-process). `QZ_NUMA_NO_PIN` (debug-only,
/// set by the forced-shard test hook) skips pinning entirely. Shared by both worker entry
/// points in `main.rs`.
pub fn pin_worker_to_node(node: u32) {
    let skip_pin = cfg!(debug_assertions) && std::env::var_os("QZ_NUMA_NO_PIN").is_some();
    if skip_pin {
        return;
    }
    match crate::topology::detect().and_then(|t| t.nodes.into_iter().find(|n| n.id == node)) {
        Some(n) => {
            // allowed_cpus() reads /proc/self/status Cpus_allowed_list (inherited from the
            // parent); intersect so we never widen beyond what the user granted.
            let allowed = crate::topology::allowed_cpus();
            let cpus: Vec<u32> = if allowed.is_empty() {
                n.cpus.clone() // couldn't read cpuset; fall back to the node's cpus
            } else {
                n.cpus.iter().copied().filter(|c| allowed.contains(c)).collect()
            };
            if cpus.is_empty() {
                eprintln!("QZ_NUMA_PIN_FAILED: node {node} has no CPUs in the process cpuset");
                std::process::exit(3);
            }
            let n = NumaNode { id: n.id, cpus, mem_available_bytes: n.mem_available_bytes };
            if let Err(e) = bind_current_thread_to_node(&n) {
                eprintln!("QZ_NUMA_PIN_FAILED: {e}");
                std::process::exit(3);
            }
        }
        None => {
            eprintln!("QZ_NUMA_PIN_FAILED: node {node} not found");
            std::process::exit(3);
        }
    }
}

#[cfg(all(test, target_os = "linux"))]
mod tests {
    use super::*;
    use crate::topology;

    #[test]
    fn bind_to_first_node_succeeds_on_numa_box_else_skips() {
        if let Some(t) = topology::detect()
            && t.has_real_asymmetry && !t.nodes.is_empty()
        {
            let n = &t.nodes[0];
            bind_current_thread_to_node(n).expect("pin to node 0 should succeed");
        }
    }

    #[test]
    fn nodemask_sets_correct_bit() {
        assert_eq!(node_to_mask(0), 0b1);
        assert_eq!(node_to_mask(1), 0b10);
        assert_eq!(node_to_mask(3), 0b1000);
    }
}
