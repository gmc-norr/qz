/// Benchmark columnar header compression vs raw BSC.
///
/// Tests 4 approaches:
///  1. Raw ASCII + BSC (current QZ default)
///  2. Columnar + per-column BSC (header_col.rs, 6 BSC calls)
///  3. Packed columnar: all binary columns concatenated → 1 BSC call
///  4. Packed columnar with delta-encoding on X/Y/read_num
///
/// Usage: bench_header_col <fastq_file> [max_reads]
use std::collections::HashMap;
use std::env;
use std::time::Instant;

use qz_lib::compression::bsc;
use qz_lib::compression::header_col;
use qz_lib::io::fastq::FastqReader;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: bench_header_col <fastq_file> [max_reads]");
        std::process::exit(1);
    }
    let path = &args[1];
    let max_reads: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(usize::MAX);

    // ================================================================
    // Load headers
    // ================================================================
    println!("Loading headers from {}...", path);
    let t0 = Instant::now();
    let mut reader = FastqReader::from_path(path, false).expect("Failed to open FASTQ");
    let mut headers: Vec<String> = Vec::new();

    while let Some(rec) = reader.next().expect("read error") {
        headers.push(String::from_utf8(rec.id).expect("non-UTF8 header"));
        if headers.len() >= max_reads {
            break;
        }
    }
    let n_reads = headers.len();
    let raw_size: usize = headers.iter().map(|h| h.len()).sum();
    println!(
        "Loaded {} headers ({} bytes raw, avg {:.0} bytes/header) in {:.1}s\n",
        n_reads,
        raw_size,
        raw_size as f64 / n_reads as f64,
        t0.elapsed().as_secs_f64()
    );

    // Show first header for reference
    if !headers.is_empty() {
        println!("  Example: {}\n", &headers[0]);
    }

    // ================================================================
    // 1. BSC baseline (raw ASCII, same as current QZ approach)
    // ================================================================
    println!("=== 1. Raw ASCII + BSC (current default) ===");
    let t0 = Instant::now();
    let mut header_stream = Vec::new();
    for h in &headers {
        write_varint(&mut header_stream, h.len());
        header_stream.extend_from_slice(h.as_bytes());
    }
    let bsc_compressed = bsc::compress_parallel_adaptive(&header_stream).expect("BSC failed");
    let bsc_time = t0.elapsed().as_secs_f64();
    let bsc_size = bsc_compressed.len();
    println!(
        "  {} bytes ({:.2}x), {:.2}s\n",
        bsc_size,
        raw_size as f64 / bsc_size as f64,
        bsc_time
    );

    // ================================================================
    // 2. Columnar + per-column BSC (header_col.rs)
    // ================================================================
    println!("=== 2. Columnar + per-column BSC ===");
    let header_refs: Vec<&str> = headers.iter().map(|h| h.as_str()).collect();
    let t0 = Instant::now();
    let compressed = header_col::compress_headers_columnar(&header_refs)
        .expect("columnar compress failed");
    let col_time = t0.elapsed().as_secs_f64();
    let col_size = compressed.len();
    println!(
        "  {} bytes ({:.2}x), {:.2}s",
        col_size,
        raw_size as f64 / col_size as f64,
        col_time
    );

    // Verify roundtrip
    let decompressed = header_col::decompress_headers_columnar(&compressed, n_reads)
        .expect("columnar decompress failed");
    let ok = (0..n_reads).all(|i| decompressed[i] == headers[i].as_bytes());
    println!("  Roundtrip: {}\n", if ok { "OK" } else { "FAILED" });

    // ================================================================
    // 3. Packed columnar: all binary columns → 1 BSC call
    // ================================================================
    println!("=== 3. Packed columnar (1 BSC call) ===");
    let packed = pack_sra_columns(&headers);
    if let Some((packed_data, meta_size)) = packed {
        let t0 = Instant::now();
        let packed_compressed = bsc::compress_parallel_adaptive(&packed_data).expect("BSC failed");
        let packed_time = t0.elapsed().as_secs_f64();
        let packed_size = packed_compressed.len() + meta_size;
        println!(
            "  {} bytes ({:.2}x), {:.2}s  [BSC: {}, meta: {}]",
            packed_size,
            raw_size as f64 / packed_size as f64,
            packed_time,
            packed_compressed.len(),
            meta_size
        );
    } else {
        println!("  (skipped — not SRA format)");
    }

    // ================================================================
    // 4. Packed columnar with delta encoding
    // ================================================================
    println!("\n=== 4. Packed columnar + delta (1 BSC call) ===");
    let packed_delta = pack_sra_columns_delta(&headers);
    if let Some((packed_data, meta_size)) = packed_delta {
        let t0 = Instant::now();
        let packed_compressed = bsc::compress_parallel_adaptive(&packed_data).expect("BSC failed");
        let packed_time = t0.elapsed().as_secs_f64();
        let packed_size = packed_compressed.len() + meta_size;
        println!(
            "  {} bytes ({:.2}x), {:.2}s  [BSC: {}, meta: {}]",
            packed_size,
            raw_size as f64 / packed_size as f64,
            packed_time,
            packed_compressed.len(),
            meta_size
        );
    } else {
        println!("  (skipped — not SRA format)");
    }

    // ================================================================
    // 5. Packed columnar with delta + zigzag varint (variable-width)
    // ================================================================
    println!("\n=== 5. Packed columnar + delta + varint (1 BSC call) ===");
    let packed_varint = pack_sra_columns_delta_varint(&headers);
    if let Some((packed_data, meta_size)) = packed_varint {
        let t0 = Instant::now();
        let packed_compressed = bsc::compress_parallel_adaptive(&packed_data).expect("BSC failed");
        let packed_time = t0.elapsed().as_secs_f64();
        let packed_size = packed_compressed.len() + meta_size;
        println!(
            "  {} bytes ({:.2}x), {:.2}s  [BSC: {}, meta: {}, pre-BSC: {}]",
            packed_size,
            raw_size as f64 / packed_size as f64,
            packed_time,
            packed_compressed.len(),
            meta_size,
            packed_data.len()
        );
    } else {
        println!("  (skipped — not SRA format)");
    }

    // ================================================================
    // 6. Columnar + parallel BSC (rayon::join across columns)
    // ================================================================
    println!("\n=== 6. Columnar + parallel BSC (rayon) ===");
    let par_col = pack_sra_columns_parallel_bsc(&headers);
    if let Some((total_size, par_time)) = par_col {
        println!(
            "  {} bytes ({:.2}x), {:.2}s",
            total_size,
            raw_size as f64 / total_size as f64,
            par_time
        );
    } else {
        println!("  (skipped — not SRA format)");
    }

    // ================================================================
    // Summary
    // ================================================================
    println!("\n=== Summary ===");
    println!("  Raw ASCII:       {:>10} bytes", raw_size);
    println!("  1. Raw+BSC:      {:>10} bytes ({:.2}x) [{:.2}s]", bsc_size, raw_size as f64 / bsc_size as f64, bsc_time);
    println!("  2. Columnar:     {:>10} bytes ({:.2}x) [{:.2}s]", col_size, raw_size as f64 / col_size as f64, col_time);
    if let Some((packed_data, meta_size)) = pack_sra_columns(&headers) {
        let packed_compressed = bsc::compress_parallel_adaptive(&packed_data).expect("BSC failed");
        let packed_size = packed_compressed.len() + meta_size;
        println!("  3. Packed:       {:>10} bytes ({:.2}x)", packed_size, raw_size as f64 / packed_size as f64);
    }
    if let Some((packed_data, meta_size)) = pack_sra_columns_delta(&headers) {
        let packed_compressed = bsc::compress_parallel_adaptive(&packed_data).expect("BSC failed");
        let packed_size = packed_compressed.len() + meta_size;
        println!("  4. Packed+delta: {:>10} bytes ({:.2}x)", packed_size, raw_size as f64 / packed_size as f64);
    }
    if let Some((packed_data, meta_size)) = pack_sra_columns_delta_varint(&headers) {
        let packed_compressed = bsc::compress_parallel_adaptive(&packed_data).expect("BSC failed");
        let packed_size = packed_compressed.len() + meta_size;
        println!("  5. Packed+dv:    {:>10} bytes ({:.2}x)", packed_size, raw_size as f64 / packed_size as f64);
    }
    if let Some((total_size, _)) = pack_sra_columns_parallel_bsc(&headers) {
        println!("  6. Par-col BSC:  {:>10} bytes ({:.2}x)", total_size, raw_size as f64 / total_size as f64);
    }
}

// ========================================================================
// Helpers
// ========================================================================

fn write_varint(out: &mut Vec<u8>, mut value: usize) {
    while value >= 0x80 {
        out.push(((value & 0x7F) | 0x80) as u8);
        value >>= 7;
    }
    out.push(value as u8);
}

fn zigzag_encode(v: i32) -> u32 {
    ((v << 1) ^ (v >> 31)) as u32
}

fn write_varint_u32(out: &mut Vec<u8>, mut value: u32) {
    while value >= 0x80 {
        out.push(((value & 0x7F) | 0x80) as u8);
        value >>= 7;
    }
    out.push(value as u8);
}

/// Parse SRA-format header: @PREFIX.READ_NUM COMBO:LANE:TILE:X:Y[/PAIR]
fn parse_sra_header(header: &str) -> Option<(u32, &str, u8, u16, u16, u16)> {
    let without_at = header.strip_prefix('@').unwrap_or(header);
    let (name_part, comment) = without_at.split_once(' ')?;

    let read_num: u32 = name_part.rsplit('.').next()?.parse().ok()?;

    let comment_base = comment.split('/').next().unwrap_or(comment);
    let fields: Vec<&str> = comment_base.split(':').collect();
    if fields.len() < 7 {
        return None;
    }

    let lane: u8 = fields[3].parse().ok()?;
    let tile: u16 = fields[4].parse().ok()?;
    let x: u16 = fields[5].parse().ok()?;
    let y: u16 = fields[6].parse().ok()?;

    Some((read_num, comment_base, lane, tile, x, y))
}

/// Approach 3: Pack all columns into one contiguous buffer (column-major order),
/// then compress with a single BSC call.
/// Returns (packed_data, metadata_size) or None if not SRA format.
fn pack_sra_columns(headers: &[String]) -> Option<(Vec<u8>, usize)> {
    let n = headers.len();
    // Pre-allocate columns
    let mut read_nums = Vec::with_capacity(n * 4);
    let mut lanes = Vec::with_capacity(n);
    let mut tiles = Vec::with_capacity(n * 2);
    let mut xs = Vec::with_capacity(n * 2);
    let mut ys = Vec::with_capacity(n * 2);

    let mut combos: Vec<String> = Vec::new();
    let mut combo_map: HashMap<String, u8> = HashMap::new();
    let mut combo_idxs = Vec::with_capacity(n);

    for h in headers {
        let (read_num, combo_str, lane, tile, x, y) = parse_sra_header(h)?;

        // Extract combo (first 3 colon fields)
        let combo_fields: Vec<&str> = combo_str.split(':').collect();
        let combo_key = format!("{}:{}:{}", combo_fields[0], combo_fields[1], combo_fields[2]);
        let idx = if let Some(&idx) = combo_map.get(&combo_key) {
            idx
        } else {
            if combos.len() >= 255 { return None; }
            let idx = combos.len() as u8;
            combos.push(combo_key.clone());
            combo_map.insert(combo_key, idx);
            idx
        };

        read_nums.extend_from_slice(&read_num.to_le_bytes());
        combo_idxs.push(idx);
        lanes.push(lane);
        tiles.extend_from_slice(&tile.to_le_bytes());
        xs.extend_from_slice(&x.to_le_bytes());
        ys.extend_from_slice(&y.to_le_bytes());
    }

    // Concatenate all columns in order (column-major)
    let mut packed = Vec::with_capacity(read_nums.len() + combo_idxs.len() + lanes.len() + tiles.len() + xs.len() + ys.len());
    packed.extend_from_slice(&read_nums);
    packed.extend_from_slice(&combo_idxs);
    packed.extend_from_slice(&lanes);
    packed.extend_from_slice(&tiles);
    packed.extend_from_slice(&xs);
    packed.extend_from_slice(&ys);

    // Metadata: combo dict + template prefix/suffix (estimate ~200 bytes)
    let meta_size = combos.iter().map(|c| c.len() + 1).sum::<usize>() + 50;

    Some((packed, meta_size))
}

/// Approach 4: Same as 3 but delta-encode read_nums, xs, ys before packing.
fn pack_sra_columns_delta(headers: &[String]) -> Option<(Vec<u8>, usize)> {
    let n = headers.len();
    let mut read_nums = Vec::with_capacity(n * 4);
    let mut lanes = Vec::with_capacity(n);
    let mut tiles = Vec::with_capacity(n * 2);
    let mut xs = Vec::with_capacity(n * 2);
    let mut ys = Vec::with_capacity(n * 2);

    let mut combos: Vec<String> = Vec::new();
    let mut combo_map: HashMap<String, u8> = HashMap::new();
    let mut combo_idxs = Vec::with_capacity(n);

    let mut prev_rn: u32 = 0;
    let mut prev_x: u16 = 0;
    let mut prev_y: u16 = 0;

    for h in headers {
        let (read_num, combo_str, lane, tile, x, y) = parse_sra_header(h)?;

        let combo_fields: Vec<&str> = combo_str.split(':').collect();
        let combo_key = format!("{}:{}:{}", combo_fields[0], combo_fields[1], combo_fields[2]);
        let idx = if let Some(&idx) = combo_map.get(&combo_key) {
            idx
        } else {
            if combos.len() >= 255 { return None; }
            let idx = combos.len() as u8;
            combos.push(combo_key.clone());
            combo_map.insert(combo_key, idx);
            idx
        };

        // Delta encode
        let delta_rn = read_num.wrapping_sub(prev_rn);
        let delta_x = x.wrapping_sub(prev_x);
        let delta_y = y.wrapping_sub(prev_y);
        prev_rn = read_num;
        prev_x = x;
        prev_y = y;

        read_nums.extend_from_slice(&delta_rn.to_le_bytes());
        combo_idxs.push(idx);
        lanes.push(lane);
        tiles.extend_from_slice(&tile.to_le_bytes());
        xs.extend_from_slice(&delta_x.to_le_bytes());
        ys.extend_from_slice(&delta_y.to_le_bytes());
    }

    let mut packed = Vec::with_capacity(read_nums.len() + combo_idxs.len() + lanes.len() + tiles.len() + xs.len() + ys.len());
    packed.extend_from_slice(&read_nums);
    packed.extend_from_slice(&combo_idxs);
    packed.extend_from_slice(&lanes);
    packed.extend_from_slice(&tiles);
    packed.extend_from_slice(&xs);
    packed.extend_from_slice(&ys);

    let meta_size = combos.iter().map(|c| c.len() + 1).sum::<usize>() + 50;
    Some((packed, meta_size))
}

/// Approach 6: Per-column BSC but compress all 6 columns in parallel with rayon.
fn pack_sra_columns_parallel_bsc(headers: &[String]) -> Option<(usize, f64)> {
    let n = headers.len();
    let mut read_nums = Vec::with_capacity(n * 4);
    let mut lanes = Vec::with_capacity(n);
    let mut tiles = Vec::with_capacity(n * 2);
    let mut xs = Vec::with_capacity(n * 2);
    let mut ys = Vec::with_capacity(n * 2);

    let mut combos: Vec<String> = Vec::new();
    let mut combo_map: HashMap<String, u8> = HashMap::new();
    let mut combo_idxs = Vec::with_capacity(n);

    for h in headers {
        let (read_num, combo_str, lane, tile, x, y) = parse_sra_header(h)?;
        let combo_fields: Vec<&str> = combo_str.split(':').collect();
        let combo_key = format!("{}:{}:{}", combo_fields[0], combo_fields[1], combo_fields[2]);
        let idx = if let Some(&idx) = combo_map.get(&combo_key) {
            idx
        } else {
            if combos.len() >= 255 { return None; }
            let idx = combos.len() as u8;
            combos.push(combo_key.clone());
            combo_map.insert(combo_key, idx);
            idx
        };

        read_nums.extend_from_slice(&read_num.to_le_bytes());
        combo_idxs.push(idx);
        lanes.push(lane);
        tiles.extend_from_slice(&tile.to_le_bytes());
        xs.extend_from_slice(&x.to_le_bytes());
        ys.extend_from_slice(&y.to_le_bytes());
    }

    let meta_size = combos.iter().map(|c| c.len() + 1).sum::<usize>() + 50;

    let t0 = Instant::now();
    // Compress all 6 columns in parallel using rayon
    let (c_rn, (c_ci, (c_ln, (c_ti, (c_xs, c_ys))))) = rayon::join(
        || bsc::compress_parallel_adaptive(&read_nums),
        || rayon::join(
            || bsc::compress_parallel_adaptive(&combo_idxs),
            || rayon::join(
                || bsc::compress_parallel_adaptive(&lanes),
                || rayon::join(
                    || bsc::compress_parallel_adaptive(&tiles),
                    || rayon::join(
                        || bsc::compress_parallel_adaptive(&xs),
                        || bsc::compress_parallel_adaptive(&ys),
                    ),
                ),
            ),
        ),
    );
    let elapsed = t0.elapsed().as_secs_f64();

    let total = c_rn.ok()?.len() + c_ci.ok()?.len() + c_ln.ok()?.len()
        + c_ti.ok()?.len() + c_xs.ok()?.len() + c_ys.ok()?.len() + meta_size;

    Some((total, elapsed))
}

/// Approach 5: Delta + zigzag varint encoding for variable-width columns.
/// read_num deltas, x deltas, y deltas use zigzag + varint (most deltas are small).
/// tile uses varint directly. lane and combo stay as raw bytes.
fn pack_sra_columns_delta_varint(headers: &[String]) -> Option<(Vec<u8>, usize)> {
    let n = headers.len();

    let mut combos: Vec<String> = Vec::new();
    let mut combo_map: HashMap<String, u8> = HashMap::new();

    // We'll build separate streams then concatenate
    let mut rn_stream = Vec::with_capacity(n * 2);
    let mut combo_idxs = Vec::with_capacity(n);
    let mut lanes = Vec::with_capacity(n);
    let mut tile_stream = Vec::with_capacity(n * 2);
    let mut x_stream = Vec::with_capacity(n * 2);
    let mut y_stream = Vec::with_capacity(n * 2);

    let mut prev_rn: i32 = 0;
    let mut prev_x: i32 = 0;
    let mut prev_y: i32 = 0;

    for h in headers {
        let (read_num, combo_str, lane, tile, x, y) = parse_sra_header(h)?;

        let combo_fields: Vec<&str> = combo_str.split(':').collect();
        let combo_key = format!("{}:{}:{}", combo_fields[0], combo_fields[1], combo_fields[2]);
        let idx = if let Some(&idx) = combo_map.get(&combo_key) {
            idx
        } else {
            if combos.len() >= 255 { return None; }
            let idx = combos.len() as u8;
            combos.push(combo_key.clone());
            combo_map.insert(combo_key, idx);
            idx
        };

        let delta_rn = (read_num as i32) - prev_rn;
        let delta_x = (x as i32) - prev_x;
        let delta_y = (y as i32) - prev_y;
        prev_rn = read_num as i32;
        prev_x = x as i32;
        prev_y = y as i32;

        write_varint_u32(&mut rn_stream, zigzag_encode(delta_rn));
        combo_idxs.push(idx);
        lanes.push(lane);
        write_varint_u32(&mut tile_stream, tile as u32);
        write_varint_u32(&mut x_stream, zigzag_encode(delta_x));
        write_varint_u32(&mut y_stream, zigzag_encode(delta_y));
    }

    // Concatenate all streams with length prefixes so we can split on decompress
    let mut packed = Vec::new();
    write_varint(&mut packed, rn_stream.len());
    packed.extend_from_slice(&rn_stream);
    packed.extend_from_slice(&combo_idxs); // fixed size = n
    packed.extend_from_slice(&lanes);       // fixed size = n
    write_varint(&mut packed, tile_stream.len());
    packed.extend_from_slice(&tile_stream);
    write_varint(&mut packed, x_stream.len());
    packed.extend_from_slice(&x_stream);
    write_varint(&mut packed, y_stream.len());
    packed.extend_from_slice(&y_stream);

    let meta_size = combos.iter().map(|c| c.len() + 1).sum::<usize>() + 50;
    Some((packed, meta_size))
}
