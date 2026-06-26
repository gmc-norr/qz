//! Assemble worker part files into the final output(s) in chunk-range order.

use std::fs::{File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};

/// Append the entire contents of `src` to the open `dst`. Tries copy_file_range
/// (an optimized in-kernel copy where the kernel/FS supports it) then falls back
/// to a buffered read/write loop.
fn append_file(dst: &mut File, src: &Path) -> io::Result<()> {
    let mut s = File::open(src)?;
    #[cfg(target_os = "linux")]
    {
        let len = s.metadata()?.len();
        let mut remaining = len as usize;
        use std::os::unix::io::AsRawFd;
        let (sfd, dfd) = (s.as_raw_fd(), dst.as_raw_fd());
        let mut off_in: i64 = 0;
        while remaining > 0 {
            let n = unsafe {
                libc::copy_file_range(sfd, &mut off_in, dfd, std::ptr::null_mut(), remaining, 0)
            };
            if n <= 0 { break; }
            remaining -= n as usize;
        }
        if remaining == 0 { return Ok(()); }
        use std::io::Seek;
        s.seek(io::SeekFrom::Start((len as usize - remaining) as u64))?;
    }
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = s.read(&mut buf)?;
        if n == 0 { break; }
        dst.write_all(&buf[..n])?;
    }
    Ok(())
}

/// RAII cleanup for an assembly temp file: removed on any early return unless disarmed.
/// Field `.1` is the "armed" flag; disarm via [`TmpGuard::disarm`] after publish.
pub struct TmpGuard(pub PathBuf, pub bool);
impl Drop for TmpGuard {
    fn drop(&mut self) {
        if self.1 { let _ = std::fs::remove_file(&self.0); }
    }
}
impl TmpGuard {
    /// Disarm so the temp is NOT removed on drop (call after a successful publish).
    pub fn disarm(&mut self) { self.1 = false; }
}

/// RAII cleanup for N direct-write temps: all removed on early drop unless disarmed.
/// Generalizes [`TmpGuard`] to the multi-output (paired) case — single-end passes one
/// path, paired passes two (R1/R2).
pub struct MultiTmpGuard {
    paths: Vec<PathBuf>,
    armed: bool,
}
impl Drop for MultiTmpGuard {
    fn drop(&mut self) {
        if self.armed {
            for p in &self.paths {
                let _ = std::fs::remove_file(p);
            }
        }
    }
}
impl MultiTmpGuard {
    /// Disarm so the temps are NOT removed on drop (call after a successful publish).
    fn disarm(&mut self) {
        self.armed = false;
    }
}

/// Create N pre-sized sibling temps (one per final output) for direct-write — workers
/// seek to disjoint regions and write into them. Returns the temp paths (parallel to
/// `finals`) and a guard that removes ALL of them on early drop (failure). Single-end
/// passes 1 output; paired passes 2 (R1/R2). Output existence / `--force` is handled
/// by `run_sharded`'s preflight, so this does not re-check it.
pub fn make_direct_temps(
    finals: &[PathBuf],
    totals: &[u64],
) -> anyhow::Result<(Vec<PathBuf>, MultiTmpGuard)> {
    assert_eq!(finals.len(), totals.len(), "make_direct_temps: finals/totals length mismatch");
    let mut tmps = Vec::with_capacity(finals.len());
    // The guard is armed from the first temp on, so a mid-loop `?` removes whatever
    // was already created.
    let mut guard = MultiTmpGuard { paths: Vec::new(), armed: true };
    for (out, &total) in finals.iter().zip(totals) {
        let tmp = sibling_tmp(out);
        let f = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&tmp)
            .map_err(|e| anyhow::anyhow!("create direct temp {}: {e}", tmp.display()))?;
        f.set_len(total)
            .map_err(|e| anyhow::anyhow!("ftruncate {}: {e}", tmp.display()))?;
        guard.paths.push(tmp.clone());
        tmps.push(tmp);
    }
    Ok((tmps, guard))
}

/// Publish N completed direct temps to their finals, then disarm the guard. For ONE
/// output: a single atomic rename. For TWO: the backup-and-rollback two-file publish
/// (mirrors qz-lib's paired publish — two files cannot be one atomic unit).
pub fn publish_direct_multi(
    tmps: &[PathBuf],
    finals: &[PathBuf],
    mut guard: MultiTmpGuard,
) -> anyhow::Result<()> {
    assert_eq!(tmps.len(), finals.len(), "publish_direct_multi: tmps/finals length mismatch");
    match tmps.len() {
        1 => {
            std::fs::rename(&tmps[0], &finals[0]).map_err(|e| {
                anyhow::anyhow!("rename {} -> {}: {e}", tmps[0].display(), finals[0].display())
            })?;
        }
        2 => publish_pair_local(&tmps[0], &tmps[1], &finals[0], &finals[1])?,
        n => anyhow::bail!("publish_direct_multi: unsupported output count {n}"),
    }
    guard.disarm();
    Ok(())
}

/// Single-end: concatenate parts into a temp file, then atomic-rename to `out`.
pub fn assemble_single(parts: &[PathBuf], out: &Path, force: bool) -> anyhow::Result<()> {
    if out.exists() && !force {
        anyhow::bail!("Output file already exists: {}\nUse --force to overwrite", out.display());
    }
    let tmp = sibling_tmp(out);
    let mut tg = TmpGuard(tmp.clone(), true);
    {
        let mut dst = OpenOptions::new().write(true).create_new(true).open(&tmp)
            .map_err(|e| anyhow::anyhow!("create assembly temp {}: {e}", tmp.display()))?;
        for p in parts {
            append_file(&mut dst, p).map_err(|e| anyhow::anyhow!("append {}: {e}", p.display()))?;
        }
        dst.flush()?;
    }
    std::fs::rename(&tmp, out)
        .map_err(|e| anyhow::anyhow!("rename {} -> {}: {e}", tmp.display(), out.display()))?;
    tg.disarm();
    for p in parts { let _ = std::fs::remove_file(p); }
    Ok(())
}

/// Paired: assemble BOTH _R1 and _R2 fully into temp files first, then publish via
/// backup-and-rollback (two files cannot be one atomic unit).
pub fn assemble_paired(
    r1_parts: &[PathBuf],
    r2_parts: &[PathBuf],
    prefix: &Path,
    force: bool,
) -> anyhow::Result<()> {
    let out1 = with_suffix(prefix, "_R1.fastq");
    let out2 = with_suffix(prefix, "_R2.fastq");
    if !force && (out1.exists() || out2.exists()) {
        anyhow::bail!("Output exists ({} / {}). Use --force.", out1.display(), out2.display());
    }
    let tmp1 = sibling_tmp(&out1);
    let tmp2 = sibling_tmp(&out2);
    let mut g1 = TmpGuard(tmp1.clone(), true);
    let mut g2 = TmpGuard(tmp2.clone(), true);
    for (parts, tmp) in [(r1_parts, &tmp1), (r2_parts, &tmp2)] {
        let mut dst = OpenOptions::new().write(true).create_new(true).open(tmp)?;
        for p in parts.iter() { append_file(&mut dst, p)?; }
        dst.flush()?;
    }
    publish_pair_local(&tmp1, &tmp2, &out1, &out2)?;
    g1.disarm();
    g2.disarm();
    for p in r1_parts.iter().chain(r2_parts) { let _ = std::fs::remove_file(p); }
    Ok(())
}

/// Local backup-and-rollback two-file publish (mirrors qz-lib paired publish_pair).
fn publish_pair_local(tmp1: &Path, tmp2: &Path, out1: &Path, out2: &Path) -> anyhow::Result<()> {
    let bak1 = sibling_tmp(out1);
    let bak2 = sibling_tmp(out2);
    let had1 = out1.exists();
    let had2 = out2.exists();
    if had1 { std::fs::rename(out1, &bak1)?; }
    if had2 && let Err(e) = std::fs::rename(out2, &bak2) {
        if had1 { let _ = std::fs::rename(&bak1, out1); }
        return Err(anyhow::anyhow!("backup out2 failed: {e}"));
    }
    if let Err(e) = std::fs::rename(tmp1, out1) {
        if had1 { let _ = std::fs::rename(&bak1, out1); }
        if had2 { let _ = std::fs::rename(&bak2, out2); }
        return Err(anyhow::anyhow!("publish out1 failed: {e}"));
    }
    if let Err(e) = std::fs::rename(tmp2, out2) {
        let _ = std::fs::remove_file(out1);
        if had1 { let _ = std::fs::rename(&bak1, out1); }
        if had2 { let _ = std::fs::rename(&bak2, out2); }
        return Err(anyhow::anyhow!("publish out2 failed: {e}"));
    }
    let _ = std::fs::remove_file(&bak1);
    let _ = std::fs::remove_file(&bak2);
    Ok(())
}

/// Append a suffix to a path's filename (e.g. "_R1.fastq").
pub fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

fn sibling_tmp(out: &Path) -> PathBuf {
    let pid = std::process::id();
    let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    with_suffix(out, &format!(".{pid}.{n}.numatmp"))
}

static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write(p: &std::path::Path, bytes: &[u8]) {
        std::fs::File::create(p).unwrap().write_all(bytes).unwrap();
    }

    #[test]
    fn assemble_single_concatenates_in_order() {
        let d = tempfile::tempdir().unwrap();
        let p0 = d.path().join("p0"); let p1 = d.path().join("p1");
        write(&p0, b"AAA"); write(&p1, b"BBB");
        let out = d.path().join("out.fastq");
        assemble_single(&[p0.clone(), p1.clone()], &out, false).unwrap();
        assert_eq!(std::fs::read(&out).unwrap(), b"AAABBB");
        assert!(!p0.exists() && !p1.exists());
    }

    #[test]
    fn assemble_single_refuses_existing_without_force() {
        let d = tempfile::tempdir().unwrap();
        let p0 = d.path().join("p0"); write(&p0, b"X");
        let out = d.path().join("out.fastq"); write(&out, b"OLD");
        assert!(assemble_single(&[p0], &out, false).is_err());
        assert_eq!(std::fs::read(&out).unwrap(), b"OLD");
    }

    #[test]
    fn assemble_single_force_replaces() {
        let d = tempfile::tempdir().unwrap();
        let p0 = d.path().join("p0"); write(&p0, b"NEW");
        let out = d.path().join("out.fastq"); write(&out, b"OLD");
        assemble_single(&[p0], &out, true).unwrap();
        assert_eq!(std::fs::read(&out).unwrap(), b"NEW");
    }

    #[test]
    fn assemble_single_missing_part_errors_and_leaves_no_temp() {
        let d = tempfile::tempdir().unwrap();
        let present = d.path().join("p0"); write(&present, b"AAA");
        let missing = d.path().join("nope");
        let out = d.path().join("out.fastq");
        assert!(assemble_single(&[present, missing], &out, false).is_err());
        assert!(!out.exists(), "no output published on assembly failure");
        let leftover = std::fs::read_dir(d.path()).unwrap().filter_map(|e| e.ok())
            .any(|e| e.file_name().to_string_lossy().contains(".numatmp"));
        assert!(!leftover, "assembly temp leaked");
    }

    #[test]
    fn assemble_paired_missing_part_errors_and_leaves_no_temp() {
        let d = tempfile::tempdir().unwrap();
        let a1 = d.path().join("a1"); write(&a1, b"R1");
        let a2 = d.path().join("a2"); // never created -> R2 assembly fails
        let prefix = d.path().join("out");
        assert!(assemble_paired(&[a1], &[a2], &prefix, false).is_err());
        assert!(!with_suffix(&prefix, "_R1.fastq").exists(), "no R1 published on failure");
        assert!(!with_suffix(&prefix, "_R2.fastq").exists(), "no R2 published on failure");
        let leftover = std::fs::read_dir(d.path()).unwrap().filter_map(|e| e.ok())
            .any(|e| e.file_name().to_string_lossy().contains(".numatmp"));
        assert!(!leftover, "assembly temps leaked");
    }

    #[test]
    fn assemble_paired_writes_both_outputs() {
        let d = tempfile::tempdir().unwrap();
        let a1 = d.path().join("a1"); let b1 = d.path().join("b1");
        let a2 = d.path().join("a2"); let b2 = d.path().join("b2");
        write(&a1, b"R1a"); write(&b1, b"R1b"); write(&a2, b"R2a"); write(&b2, b"R2b");
        let prefix = d.path().join("out");
        assemble_paired(&[a1, b1], &[a2, b2], &prefix, false).unwrap();
        let r1 = with_suffix(&prefix, "_R1.fastq");
        let r2 = with_suffix(&prefix, "_R2.fastq");
        assert_eq!(std::fs::read(&r1).unwrap(), b"R1aR1b");
        assert_eq!(std::fs::read(&r2).unwrap(), b"R2aR2b");
    }
}
