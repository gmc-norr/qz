//! Parse Linux NUMA topology from /sys, intersected with the process's allowed
//! CPU and memory-node masks. Returns None on non-Linux / absent /sys.

#[derive(Clone, Debug)]
pub struct NumaNode {
    pub id: u32,
    pub cpus: Vec<u32>,
    /// Estimated memory a worker can use on this node WITHOUT swapping: free RAM
    /// plus reclaimable page cache + slab (a per-node analogue of /proc/meminfo
    /// `MemAvailable`). Using strict `MemFree` here would gate a node out whenever
    /// it holds a lot of reclaimable cache (common on I/O-heavy boxes), needlessly
    /// disabling NUMA sharding.
    pub mem_available_bytes: u64,
}

#[derive(Clone, Debug)]
pub struct Topology {
    pub nodes: Vec<NumaNode>,
    /// true iff some remote distance strictly exceeds the local distance.
    pub has_real_asymmetry: bool,
}

const SYS_NODE: &str = "/sys/devices/system/node";

/// Detect topology on the live system. None if not Linux / /sys missing / unparsable.
pub fn detect() -> Option<Topology> {
    parse_sys_tree(std::path::Path::new(SYS_NODE)).ok()
}

/// Parse a node tree rooted at `base` (real or a test fixture).
pub fn parse_sys_tree(base: &std::path::Path) -> anyhow::Result<Topology> {
    let online = std::fs::read_to_string(base.join("online"))
        .map_err(|e| anyhow::anyhow!("no NUMA /sys online list: {e}"))?;
    let ids = parse_id_list(online.trim());
    let mut nodes = Vec::new();
    let mut max_local = 10u32;
    let mut max_remote = 10u32;
    for (row, &id) in ids.iter().enumerate() {
        let nd = base.join(format!("node{id}"));
        let cpulist = std::fs::read_to_string(nd.join("cpulist")).unwrap_or_default();
        let cpus = parse_id_list(cpulist.trim());
        let mem_available_bytes = read_node_mem(&nd.join("meminfo"));
        if let Ok(dist) = std::fs::read_to_string(nd.join("distance")) {
            let vals: Vec<u32> = dist.split_whitespace().filter_map(|x| x.parse().ok()).collect();
            for (col, &v) in vals.iter().enumerate() {
                if col == row { max_local = max_local.max(v); } else { max_remote = max_remote.max(v); }
            }
        }
        nodes.push(NumaNode { id, cpus, mem_available_bytes });
    }
    let has_real_asymmetry = nodes.len() >= 2 && max_remote > max_local;
    Ok(Topology { nodes, has_real_asymmetry })
}

impl Topology {
    /// Nodes usable given the process's allowed CPU set and allowed memory nodes:
    /// intersect each node's CPUs with `allowed_cpus`, drop CPUless results, drop
    /// nodes not in `allowed_mem_nodes`.
    pub fn usable_nodes(&self, allowed_cpus: &[u32], allowed_mem_nodes: &[u32]) -> Vec<NumaNode> {
        self.nodes
            .iter()
            .filter(|n| allowed_mem_nodes.contains(&n.id))
            .filter_map(|n| {
                let cpus: Vec<u32> =
                    n.cpus.iter().copied().filter(|c| allowed_cpus.contains(c)).collect();
                if cpus.is_empty() {
                    None
                } else {
                    Some(NumaNode { id: n.id, cpus, mem_available_bytes: n.mem_available_bytes })
                }
            })
            .collect()
    }
}

/// Parse a Linux cpu/node list like "0-3,7,9-11" into a Vec of ids.
fn parse_id_list(s: &str) -> Vec<u32> {
    let mut out = Vec::new();
    for part in s.split(',').filter(|p| !p.is_empty()) {
        if let Some((a, b)) = part.split_once('-') {
            if let (Ok(a), Ok(b)) = (a.parse::<u32>(), b.parse::<u32>()) {
                out.extend(a..=b);
            }
        } else if let Ok(v) = part.parse::<u32>() {
            out.push(v);
        }
    }
    out
}

/// Estimate per-node available memory (bytes): free RAM plus reclaimable page
/// cache + slab. A per-node analogue of `/proc/meminfo` `MemAvailable` (node
/// meminfo has no `MemAvailable` field). Reclaimable file cache and `SReclaimable`
/// slab can be evicted under pressure, so a worker can use them; anonymous memory
/// is excluded (it cannot), keeping the estimate from over-committing. Returns 0
/// if the file is unreadable.
fn read_node_mem(meminfo: &std::path::Path) -> u64 {
    let s = std::fs::read_to_string(meminfo).unwrap_or_default();
    // MemFree + reclaimable file cache + reclaimable slab. (Node meminfo has no
    // MemAvailable field, so we sum its components; absent fields contribute 0,
    // degrading gracefully to MemFree.)
    field_bytes(&s, "MemFree:")
        + field_bytes(&s, "Active(file):")
        + field_bytes(&s, "Inactive(file):")
        + field_bytes(&s, "SReclaimable:")
}

/// Parse one `key`'d value (in kB) from node meminfo (`Node N <key>  <kb> kB`),
/// returning bytes. 0 if the key is absent. Case-sensitive so `Active(file):`
/// does not match `Inactive(file):`.
fn field_bytes(meminfo: &str, key: &str) -> u64 {
    for line in meminfo.lines() {
        if let Some(rest) = line.split(key).nth(1)
            && let Some(tok) = rest.split_whitespace().next()
            && let Ok(kb) = tok.parse::<u64>()
        {
            return kb * 1024;
        }
    }
    0
}

/// Process's allowed CPUs from /proc/self/status "Cpus_allowed_list:" (cpuset/cgroup-aware).
#[cfg(target_os = "linux")]
pub fn allowed_cpus() -> Vec<u32> {
    read_status_list("Cpus_allowed_list:").unwrap_or_default()
}

/// Process's allowed memory nodes from /proc/self/status "Mems_allowed_list:".
#[cfg(target_os = "linux")]
pub fn allowed_mem_nodes() -> Vec<u32> {
    read_status_list("Mems_allowed_list:").unwrap_or_default()
}

#[cfg(target_os = "linux")]
fn read_status_list(key: &str) -> Option<Vec<u32>> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in s.lines() {
        if let Some(rest) = line.strip_prefix(key) {
            return Some(parse_id_list(rest.trim()));
        }
    }
    None
}

#[cfg(not(target_os = "linux"))]
pub fn allowed_cpus() -> Vec<u32> { Vec::new() }
#[cfg(not(target_os = "linux"))]
pub fn allowed_mem_nodes() -> Vec<u32> { Vec::new() }

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn available_counts_free_plus_reclaimable() {
        // A node with little MemFree but lots of reclaimable file cache + slab is
        // still usable: a worker can reclaim that cache. Available must sum them.
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("meminfo");
        // 1 GiB free + 2 GiB Active(file) + 4 GiB Inactive(file) + 1 GiB SReclaimable
        // = 8 GiB available. (Anonymous/used memory is NOT counted.)
        std::fs::write(
            &p,
            "Node 0 MemTotal:      100000000 kB\n\
             Node 0 MemFree:         1048576 kB\n\
             Node 0 MemUsed:        98951424 kB\n\
             Node 0 Active(file):    2097152 kB\n\
             Node 0 Inactive(file):  4194304 kB\n\
             Node 0 SReclaimable:    1048576 kB\n\
             Node 0 SUnreclaim:     50000000 kB\n",
        )
        .unwrap();
        assert_eq!(read_node_mem(&p), 8u64 << 30);
    }

    #[test]
    fn available_falls_back_to_free_when_no_reclaimable_fields() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("meminfo");
        std::fs::write(&p, "Node 0 MemFree:         1048576 kB\n").unwrap();
        assert_eq!(read_node_mem(&p), 1u64 << 30);
    }

    fn fixture(nodes: &[(u32, &str, u64)], distances: &[&str]) -> tempfile::TempDir {
        let d = tempfile::tempdir().unwrap();
        let base = d.path().join("devices/system/node");
        std::fs::create_dir_all(&base).unwrap();
        let online: Vec<String> = nodes.iter().map(|(id, _, _)| id.to_string()).collect();
        std::fs::write(base.join("online"), online.join(",")).unwrap();
        for (i, (id, cpulist, memfree_kb)) in nodes.iter().enumerate() {
            let nd = base.join(format!("node{id}"));
            std::fs::create_dir_all(&nd).unwrap();
            std::fs::write(nd.join("cpulist"), cpulist).unwrap();
            let mut mi = std::fs::File::create(nd.join("meminfo")).unwrap();
            writeln!(mi, "Node {id} MemFree:       {memfree_kb} kB").unwrap();
            std::fs::write(nd.join("distance"), distances[i]).unwrap();
        }
        d
    }

    #[test]
    fn two_real_nodes_parse_with_asymmetry() {
        let d = fixture(
            &[(0, "0-17,36-53", 64_000_000), (1, "18-35,54-71", 64_000_000)],
            &["10 21", "21 10"],
        );
        let topo = parse_sys_tree(&d.path().join("devices/system/node")).unwrap();
        assert_eq!(topo.nodes.len(), 2);
        assert_eq!(topo.nodes[0].cpus.len(), 36);
        assert!(topo.has_real_asymmetry);
    }

    #[test]
    fn single_node_no_asymmetry() {
        let d = fixture(&[(0, "0-71", 128_000_000)], &["10"]);
        let topo = parse_sys_tree(&d.path().join("devices/system/node")).unwrap();
        assert_eq!(topo.nodes.len(), 1);
        assert!(!topo.has_real_asymmetry);
    }

    #[test]
    fn fake_vnuma_flat_distance_is_not_asymmetric() {
        let d = fixture(
            &[(0, "0-3", 8_000_000), (1, "4-7", 8_000_000)],
            &["10 10", "10 10"],
        );
        let topo = parse_sys_tree(&d.path().join("devices/system/node")).unwrap();
        assert!(!topo.has_real_asymmetry, "flat 10/10 must not count as NUMA");
    }

    #[test]
    fn cpuless_node_excluded() {
        let d = fixture(
            &[(0, "0-17", 64_000_000), (1, "", 64_000_000)],
            &["10 21", "21 10"],
        );
        let topo = parse_sys_tree(&d.path().join("devices/system/node")).unwrap();
        let usable = topo.usable_nodes(&(0..18).collect::<Vec<_>>(), &[0, 1]);
        assert_eq!(usable.len(), 1, "memoryful-but-cpuless node 1 is excluded");
    }

    #[test]
    fn cpuset_intersection_drops_disallowed_cpus() {
        let d = fixture(
            &[(0, "0-17,36-53", 64_000_000), (1, "18-35,54-71", 64_000_000)],
            &["10 21", "21 10"],
        );
        let topo = parse_sys_tree(&d.path().join("devices/system/node")).unwrap();
        let usable = topo.usable_nodes(&(0..18).chain(36..54).collect::<Vec<_>>(), &[0, 1]);
        assert_eq!(usable.len(), 1);
        assert_eq!(usable[0].id, 0);
    }
}
