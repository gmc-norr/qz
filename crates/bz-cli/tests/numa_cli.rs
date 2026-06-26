//! NUMA decode-sharding CLI test: a forced 2-way shard (debug-only `QZ_NUMA_FORCE_WORKERS`
//! hook, children unpinned) must reconstruct a record-identical BAM versus the in-process
//! decode, exercising the whole driver → re-exec workers → BGZF assembly path on any box.
//!
//! bz's smallest chunk preset is 500K records/chunk (level 1), so forcing more than one
//! chunk needs more than 500K records — hence the larger fixture. The forced-shard hook
//! only fires in debug builds (`cargo test` default), which is also where this runs.

use assert_cmd::Command;
use std::io::{Read, Write};
use std::path::Path;
use tempfile::TempDir;

/// Minimal valid BAM with `num_records` synthetic mapped records (1 reference, 10bp reads).
fn create_test_bam(path: &Path, num_records: usize) {
    let mut w = noodles::bgzf::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_all(b"BAM\x01").unwrap();
    let header = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n";
    w.write_all(&(header.len() as i32).to_le_bytes()).unwrap();
    w.write_all(header).unwrap();
    w.write_all(&1i32.to_le_bytes()).unwrap(); // n_ref
    let ref_name = b"chr1\0";
    w.write_all(&(ref_name.len() as i32).to_le_bytes()).unwrap();
    w.write_all(ref_name).unwrap();
    w.write_all(&248956422i32.to_le_bytes()).unwrap();
    for i in 0..num_records {
        let read_name = format!("read{i:06}\0");
        let l_read_name = read_name.len() as u8;
        let n_cigar_op: u16 = 1;
        let l_seq: i32 = 10;
        let cigar_op: u32 = 10 << 4; // 10M
        let seq_bytes: [u8; 5] = [0x12, 0x48, 0x12, 0x48, 0x12];
        let qual_bytes = [30u8; 10];
        let block_size: i32 =
            32 + l_read_name as i32 + 4 * n_cigar_op as i32 + seq_bytes.len() as i32 + l_seq;
        w.write_all(&block_size.to_le_bytes()).unwrap();
        w.write_all(&0i32.to_le_bytes()).unwrap(); // refID
        w.write_all(&(i as i32 * 100).to_le_bytes()).unwrap(); // pos
        w.write_all(&[l_read_name]).unwrap();
        w.write_all(&[60u8]).unwrap(); // mapq
        w.write_all(&0u16.to_le_bytes()).unwrap(); // bin
        w.write_all(&n_cigar_op.to_le_bytes()).unwrap();
        w.write_all(&0u16.to_le_bytes()).unwrap(); // flag
        w.write_all(&l_seq.to_le_bytes()).unwrap();
        w.write_all(&(-1i32).to_le_bytes()).unwrap(); // next_refID
        w.write_all(&(-1i32).to_le_bytes()).unwrap(); // next_pos
        w.write_all(&0i32.to_le_bytes()).unwrap(); // tlen
        w.write_all(read_name.as_bytes()).unwrap();
        w.write_all(&cigar_op.to_le_bytes()).unwrap();
        w.write_all(&seq_bytes).unwrap();
        w.write_all(&qual_bytes).unwrap();
    }
    w.finish().unwrap();
}

/// Inflate a whole BGZF/BAM file to its uncompressed bytes (robust to BGZF reframing).
fn inflate_bam(path: &Path) -> Vec<u8> {
    let mut r = noodles::bgzf::io::Reader::new(std::fs::File::open(path).unwrap());
    let mut buf = Vec::new();
    r.read_to_end(&mut buf).unwrap();
    buf
}

fn bz() -> Command {
    Command::cargo_bin("bz").expect("bz binary not built")
}

#[test]
fn forced_shard_decode_matches_in_process() {
    let temp = TempDir::new().unwrap();
    let p = temp.path();
    let bam = p.join("in.bam");
    // > 500K records at level 1 (500K/chunk) → multiple chunks, so the 2-way shard
    // splits into real disjoint ranges.
    create_test_bam(&bam, 600_000);

    let archive = p.join("a.bz");
    bz().arg("compress")
        .arg("-i").arg(&bam)
        .arg("-o").arg(&archive)
        .arg("-l").arg("1")
        .arg("-w").arg(p)
        .assert()
        .success();

    // Baseline: in-process decode (--numa off).
    let ref_bam = p.join("ref.bam");
    bz().arg("decompress")
        .arg("--numa").arg("off")
        .arg("-i").arg(&archive)
        .arg("-o").arg(&ref_bam)
        .arg("-w").arg(p)
        .assert()
        .success();

    // Forced 2-way shard: spawns two unpinned re-exec workers, each decoding a disjoint
    // chunk range to a part, then the driver strips EOFs and concatenates.
    let shard_bam = p.join("shard.bam");
    let spawn_log = p.join("spawn.log");
    bz().arg("decompress")
        .arg("--numa").arg("auto")
        .arg("-i").arg(&archive)
        .arg("-o").arg(&shard_bam)
        .arg("-w").arg(p)
        .env("QZ_NUMA_FORCE_WORKERS", "2")
        .env("QZ_NUMA_SPAWN_LOG", &spawn_log)
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&spawn_log).unwrap().trim(),
        "2",
        "forced shard must spawn exactly 2 workers"
    );
    assert_eq!(
        inflate_bam(&ref_bam),
        inflate_bam(&shard_bam),
        "forced-shard decode is not record-identical to the in-process decode"
    );
}

#[test]
fn forced_shard_compress_matches_in_process() {
    let temp = TempDir::new().unwrap();
    let p = temp.path();
    let bam = p.join("in.bam");
    // > 500K records at level 1 (500K/chunk) → multiple chunks, so the 2-way shard
    // splits into real disjoint ranges.
    create_test_bam(&bam, 600_000);

    // Baseline: in-process compress (--numa off).
    let ref_archive = p.join("ref.bz");
    bz().arg("compress")
        .arg("--numa")
        .arg("off")
        .arg("-i")
        .arg(&bam)
        .arg("-o")
        .arg(&ref_archive)
        .arg("-l")
        .arg("1")
        .arg("-w")
        .arg(p)
        .assert()
        .success();

    // Forced 2-way shard: two unpinned re-exec workers each compress a disjoint chunk
    // range to a `.bz` part, then the driver concatenates them. Compress is
    // byte-identical, so the assembled archive must equal the in-process archive
    // byte-for-byte (a stronger oracle than the decode side's record-identity).
    let shard_archive = p.join("shard.bz");
    let spawn_log = p.join("spawn.log");
    bz().arg("compress")
        .arg("--numa")
        .arg("auto")
        .arg("-i")
        .arg(&bam)
        .arg("-o")
        .arg(&shard_archive)
        .arg("-l")
        .arg("1")
        .arg("-w")
        .arg(p)
        .env("QZ_NUMA_FORCE_WORKERS", "2")
        .env("QZ_NUMA_SPAWN_LOG", &spawn_log)
        .assert()
        .success();

    assert_eq!(
        std::fs::read_to_string(&spawn_log).unwrap().trim(),
        "2",
        "forced shard must spawn exactly 2 workers"
    );
    assert_eq!(
        std::fs::read(&ref_archive).unwrap(),
        std::fs::read(&shard_archive).unwrap(),
        "forced-shard compress archive is not byte-identical to the in-process archive"
    );
}

/// Regression test for the auto-level (`--level 0`, the DEFAULT) NUMA-compress path —
/// the path the chunk_size-divergence data-loss bug lived in. The driver must FREEZE
/// the auto-resolved level once, so the prescan and every re-exec'd worker chunk at
/// the identical `chunk_size`; without that, a worker re-resolving the auto level from
/// live `/proc/meminfo` could pick a different `chunk_size`, diverge from the prescan's
/// virtual-offset partition, and silently drop records.
///
/// The fixture is sized above the largest auto chunk_size (L2 = 1.5M records), so the
/// 2-way shard is real (≥2 chunks) regardless of which level the box's RAM auto-picks.
/// We assert the sharded archive round-trips LOSSLESSLY to the original input's records
/// (bz roundtrip is byte-lossless on the BGZF-inflated payload). A chunk_size divergence
/// would either drop records (inflate mismatch) or, with the `expected_records` guard,
/// fail the worker outright (compress would not `.success()`).
#[test]
fn forced_shard_compress_auto_level_roundtrips() {
    let temp = TempDir::new().unwrap();
    let p = temp.path();
    let bam = p.join("in.bam");
    // > L2 chunk_size (1.5M) so auto yields >= 2 chunks on any box (L1 -> 4 chunks).
    create_test_bam(&bam, 1_550_000);

    // Forced 2-way shard at the DEFAULT level (no -l flag => level 0 => auto). The
    // driver freezes the auto level into the worker advanced_json, so both workers
    // chunk at the same frozen chunk_size as the prescan's partition.
    let shard_bz = p.join("shard.bz");
    let spawn_log = p.join("spawn.log");
    bz().arg("compress")
        .arg("--numa").arg("auto")
        .arg("-i").arg(&bam)
        .arg("-o").arg(&shard_bz)
        .arg("-w").arg(p)
        .env("QZ_NUMA_FORCE_WORKERS", "2")
        .env("QZ_NUMA_SPAWN_LOG", &spawn_log)
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&spawn_log).unwrap().trim(),
        "2",
        "auto-level forced shard must spawn exactly 2 workers"
    );

    let shard_bam = p.join("shard.bam");
    bz().arg("decompress")
        .arg("--numa").arg("off")
        .arg("-i").arg(&shard_bz)
        .arg("-o").arg(&shard_bam)
        .arg("-w").arg(p)
        .assert()
        .success();

    // Lossless: the sharded auto-level archive decodes to EXACTLY the original records
    // (no records dropped/duplicated by a chunk_size split across the workers).
    assert_eq!(
        inflate_bam(&bam),
        inflate_bam(&shard_bam),
        "auto-level forced-shard compress did not round-trip the original records losslessly"
    );
}
