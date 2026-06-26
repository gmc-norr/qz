use std::fs;
use std::process::Command;
use tempfile::TempDir;

fn qz() -> Command {
    // Cargo sets this for the qz-cli binary named `qz`.
    Command::new(env!("CARGO_BIN_EXE_qz"))
}

fn write_fastq(path: &std::path::Path, n: usize) {
    let mut s = String::new();
    for i in 0..n {
        s.push_str(&format!("@r{i}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"));
    }
    fs::write(path, s).unwrap();
}

#[test]
fn qz_merge_cli_roundtrips() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let f0 = d.join("a.fastq");
    let f1 = d.join("b.fastq");
    write_fastq(&f0, 5);
    write_fastq(&f1, 4);

    let a = d.join("a.qz");
    let b = d.join("b.qz");
    for (inp, out) in [(&f0, &a), (&f1, &b)] {
        let st = qz().args(["compress", "-i"]).arg(inp).args(["-o"]).arg(out)
            .args(["-w"]).arg(d).args(["-t", "1"]).status().unwrap();
        assert!(st.success());
    }

    let merged = d.join("merged.qz");
    let st = qz().arg("merge").args(["-o"]).arg(&merged).arg(&a).arg(&b)
        .status().unwrap();
    assert!(st.success(), "qz merge should succeed");

    // Decompress the merged archive and check the read count (9 reads).
    let out = d.join("out.fastq");
    let st = qz().args(["decompress", "-i"]).arg(&merged).args(["-o"]).arg(&out)
        .args(["-w"]).arg(d).args(["-t", "1", "-f"]).status().unwrap();
    assert!(st.success());
    let text = fs::read_to_string(&out).unwrap();
    assert_eq!(text.lines().filter(|l| l.starts_with('@')).count(), 9);
}

#[test]
fn qz_merge_cli_force_overwrite() {
    let td = TempDir::new().unwrap();
    let d = td.path();
    let f0 = d.join("a.fastq");
    write_fastq(&f0, 5);
    let a = d.join("a.qz");
    let st = qz().args(["compress", "-i"]).arg(&f0).args(["-o"]).arg(&a)
        .args(["-w"]).arg(d).args(["-t", "1"]).status().unwrap();
    assert!(st.success());

    let merged = d.join("m.qz");
    fs::write(&merged, b"pre-existing").unwrap();
    // Without -f: must refuse to overwrite.
    let st = qz().arg("merge").args(["-o"]).arg(&merged).arg(&a).status().unwrap();
    assert!(!st.success(), "merge must refuse to overwrite an existing output without --force");
    assert_eq!(fs::read(&merged).unwrap(), b"pre-existing", "refused merge must not touch the file");
    // With -f: must succeed and replace.
    let st = qz().arg("merge").args(["-o"]).arg(&merged).arg("-f").arg(&a).status().unwrap();
    assert!(st.success(), "merge --force must overwrite");
    assert_ne!(fs::read(&merged).unwrap(), b"pre-existing", "merge --force must replace the file");
}
