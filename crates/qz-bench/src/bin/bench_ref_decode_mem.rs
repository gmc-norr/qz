//! In-process reference-decode memory harness (SYSTEM allocator — NO mimalloc, so
//! valgrind/massif can attribute every allocation).
//!
//! Two modes:
//!   bench_ref_decode_mem gen <out_dir> <n_pairs> <ref_len> <chunk>
//!       Generate + compress a reference archive at <out_dir>/ref.qz (sets
//!       QZ_REF_CHUNK=<chunk> for the encode).
//!   bench_ref_decode_mem decode <archive> <out_prefix>
//!       Decode the archive into <out_prefix>_R1.fastq / _R2.fastq, printing
//!       VmRSS (from /proc/self/status) before and after. This is the call to run
//!       under valgrind --tool=massif.
//!
//! DELIBERATELY does NOT set `#[global_allocator] = mimalloc` (unlike the other
//! bench bins) so the glibc allocator is in effect.

use std::path::{Path, PathBuf};

use qz_lib::cli::{
    AdvancedOptions, CompressConfig, DecompressConfig, QualityMode, ReferenceOptions,
};

fn vmrss_kb() -> u64 {
    let s = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    for line in s.lines() {
        if let Some(rest) = line.strip_prefix("VmRSS:") {
            // "VmRSS:\t   12345 kB"
            return rest
                .split_whitespace()
                .next()
                .and_then(|n| n.parse::<u64>().ok())
                .unwrap_or(0);
        }
    }
    0
}

fn make_seq(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        v.push(b"ACGT"[((x >> 33) & 3) as usize]);
    }
    v
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// Mirror of the test's `make_sized_reference_dataset` (mapped + edits + rev-comp
/// R2 + ~2% fallback junk), but with per-read *varied* quality so the fqzcomp
/// quality codec carries real signal.
fn gen_dataset(dir: &Path, n_pairs: usize, ref_len: usize) -> (PathBuf, PathBuf, PathBuf) {
    let refseq = make_seq(ref_len, 7);
    let rf = dir.join("ref.fa");
    {
        let mut s = String::from(">chr0\n");
        s.push_str(std::str::from_utf8(&refseq).unwrap());
        s.push('\n');
        std::fs::write(&rf, s).unwrap();
    }

    let rlen = 120usize;
    let max_start = ref_len - rlen - 200;
    let mut rng = 0x1234_5678_9abc_def0u64;
    let mut next = || {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (rng >> 33) as usize
    };

    let mut s1 = String::with_capacity(n_pairs * 270);
    let mut s2 = String::with_capacity(n_pairs * 270);

    for i in 0..n_pairs {
        let unmappable = i % 50 == 0;
        let (r1bytes, r2bytes) = if unmappable {
            (
                make_seq(rlen, 0x9000_0000 + i as u64),
                make_seq(rlen, 0x7000_0000 + i as u64),
            )
        } else {
            let st1 = next() % (max_start + 1);
            let mut r1 = refseq[st1..st1 + rlen].to_vec();
            let st2 = st1 + 150;
            let r2_fwd = refseq[st2..st2 + rlen].to_vec();
            let mut r2 = revcomp(&r2_fwd);
            if i % 17 == 0 {
                let p = next() % rlen;
                r1[p] = if r1[p] == b'A' { b'C' } else { b'A' };
            }
            if i % 23 == 0 {
                let p = next() % rlen;
                r2[p] = if r2[p] == b'G' { b'T' } else { b'G' };
            }
            (r1, r2)
        };

        // Varied quality so fqzcomp has real per-read signal (not a constant blob).
        let mut q = String::with_capacity(rlen);
        for j in 0..rlen {
            let v = (next().wrapping_add(j).wrapping_add(i)) % 40;
            q.push((b'!' + v as u8) as char);
        }

        s1.push_str("@frag_");
        s1.push_str(&i.to_string());
        s1.push_str("/1\n");
        s1.push_str(std::str::from_utf8(&r1bytes).unwrap());
        s1.push_str("\n+\n");
        s1.push_str(&q);
        s1.push('\n');

        s2.push_str("@frag_");
        s2.push_str(&i.to_string());
        s2.push_str("/2\n");
        s2.push_str(std::str::from_utf8(&r2bytes).unwrap());
        s2.push_str("\n+\n");
        s2.push_str(&q);
        s2.push('\n');
    }

    let r1p = dir.join("r1.fastq");
    let r2p = dir.join("r2.fastq");
    std::fs::write(&r1p, &s1).unwrap();
    std::fs::write(&r2p, &s2).unwrap();
    (rf, r1p, r2p)
}

fn cmd_gen(out_dir: &Path, n_pairs: usize, ref_len: usize, chunk: usize) {
    std::fs::create_dir_all(out_dir).unwrap();
    let (rf, r1p, r2p) = gen_dataset(out_dir, n_pairs, ref_len);
    let archive = out_dir.join("ref.qz");

    // SAFETY: single-threaded harness setup, no other threads read env yet.
    unsafe { std::env::set_var("QZ_REF_CHUNK", chunk.to_string()) };

    let cfg = CompressConfig {
        input: vec![r1p, r2p],
        output: archive.clone(),
        working_dir: out_dir.to_path_buf(),
        threads: 0,
        no_quality: false,
        fasta: false,
        quality_mode: QualityMode::Lossless,
        ultra: None,
        force: true,
        advanced: AdvancedOptions::default(),
        reference: Some(ReferenceOptions {
            reference: rf,
            reference_index: None,
            reference_fast: false,
            reference_window: 2,
        }),
        cluster: None,
        require_prebuilt_index: false,
    };
    qz_lib::compression::compress(&cfg).expect("compress reference dataset");
    let sz = std::fs::metadata(&archive).unwrap().len();
    println!(
        "[gen] {n_pairs} pairs, ref_len {ref_len}, QZ_REF_CHUNK={chunk} -> {} ({sz} bytes)",
        archive.display()
    );
}

fn cmd_decode(archive: &Path, out_prefix: &Path) {
    let before = vmrss_kb();
    println!("[decode] VmRSS before = {before} KB");

    // Allow overriding the decode thread count via QZ_DECODE_THREADS for the
    // diagnosis sweep (0 = auto/all cores).
    let nthreads: usize = std::env::var("QZ_DECODE_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let cfg = DecompressConfig {
        input: archive.to_path_buf(),
        output: vec![out_prefix.to_path_buf()],
        working_dir: out_prefix.parent().unwrap_or(Path::new(".")).to_path_buf(),
        num_threads: nthreads,
        gzipped: false,
        gzip_level: 6,
        force: true,
    };
    qz_lib::compression::decompress(&cfg).expect("decompress reference archive");

    let after = vmrss_kb();
    println!(
        "[decode] VmRSS after  = {after} KB  (delta {} KB)",
        after as i64 - before as i64
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "usage:\n  {bin} gen <out_dir> <n_pairs> <ref_len> <chunk>\n  {bin} decode <archive> <out_prefix>",
            bin = args[0]
        );
        std::process::exit(2);
    }
    match args[1].as_str() {
        "gen" => {
            let out_dir = PathBuf::from(&args[2]);
            let n_pairs: usize = args[3].parse().unwrap();
            let ref_len: usize = args[4].parse().unwrap();
            let chunk: usize = args[5].parse().unwrap();
            cmd_gen(&out_dir, n_pairs, ref_len, chunk);
        }
        "decode" => {
            let archive = PathBuf::from(&args[2]);
            let out_prefix = PathBuf::from(&args[3]);
            cmd_decode(&archive, &out_prefix);
        }
        other => {
            eprintln!("unknown subcommand {other}");
            std::process::exit(2);
        }
    }
}
