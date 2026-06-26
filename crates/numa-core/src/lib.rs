//! Format-agnostic NUMA orchestration shared by the `qz` and `bz` CLIs.
//!
//! This crate holds the pieces of NUMA-aware multi-process sharding that carry no
//! archive-format knowledge: topology discovery, per-process CPU/memory pinning,
//! the spawn/poll/cleanup runtime, the decision gate (`gate::decide` over a
//! `GateInputs` the caller fills from its own format), and part-file assembly.
//!
//! Each tool keeps its own format-specific driver (build `GateInputs`, locate
//! chunk ranges, decode a range to a part) on top of this. Originally extracted
//! from `qz-cli::numa`; see docs/superpowers/specs/2026-06-15-numa-decode-design.md.

pub mod assemble;
pub mod gate;
pub mod pin;
pub mod runtime;
pub mod topology;

/// Sanity cap on a direct-write archive's total decoded output size (64 TiB). A
/// decoded-size table claiming more than this is treated as corrupt, so the gate
/// refuses to `ftruncate` a file that big and falls back / errors instead. Single
/// source of truth for the gate (pre-spawn) and any direct-write driver.
pub const MAX_DIRECT_OUTPUT_BYTES: u64 = 1 << 46;

/// `--numa` mode. `Fixed(n)` requests exactly `n` worker processes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NumaMode {
    Off,
    Auto,
    Fixed(u32),
}

/// clap value parser for `--numa auto|off|<N>`.
pub fn parse_numa_mode(s: &str) -> Result<NumaMode, String> {
    match s {
        "off" => Ok(NumaMode::Off),
        "auto" => Ok(NumaMode::Auto),
        other => match other.parse::<u32>() {
            Ok(0) => Err("--numa 0 is invalid (use 'off' or 'auto', or N>=1)".to_string()),
            Ok(n) => Ok(NumaMode::Fixed(n)),
            Err(_) => Err(format!(
                "--numa expects 'auto', 'off', or a positive integer (got '{other}')"
            )),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_keywords_and_numbers() {
        assert_eq!(parse_numa_mode("off").unwrap(), NumaMode::Off);
        assert_eq!(parse_numa_mode("auto").unwrap(), NumaMode::Auto);
        assert_eq!(parse_numa_mode("2").unwrap(), NumaMode::Fixed(2));
        assert!(parse_numa_mode("0").is_err());
        assert!(parse_numa_mode("banana").is_err());
        assert!(parse_numa_mode("-1").is_err());
    }
}
