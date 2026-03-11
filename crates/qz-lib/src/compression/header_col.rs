//! Columnar header compression for Illumina FASTQ headers.
//!
//! Supports two common formats:
//!
//! **SRA format** (version 0x01):
//!   `@PREFIX.READ_NUM INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y[/PAIR]`
//!   6 columns: read_num(u32), combo_idx(u8), lane(u8), tile(u16), x(u16), y(u16)
//!
//! **Casava 1.8 / DRAGEN BCL Convert format** (version 0x02):
//!   `@INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y[:UMI] COMMENT`
//!   5 columns: combo_idx(u8), lane(u8), tile(u16), x(u16), y(u16)
//!   + optional UMI column + comment template or column
//!
//! Falls back to raw BSC (version 0x00) for unrecognized formats.

use anyhow::{Context, Result};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;

use crate::compression::bsc;

/// Header format detected from parsing
#[derive(Debug, Clone, Copy, PartialEq)]
enum HeaderFormat {
    /// @PREFIX.READ_NUM COMBO:LANE:TILE:X:Y[/PAIR]
    Sra,
    /// @COMBO:LANE:TILE:X:Y COMMENT
    Casava,
}

/// Parsed SRA header fields (borrows from original header — zero allocation)
struct SraRow<'h> {
    combo: &'h str,
    read_num: u32,
    lane: u8,
    tile: u16,
    x: u16,
    y: u16,
}

/// Parsed Casava header fields (borrows from original header — zero allocation)
struct CasavaRow<'h> {
    combo: &'h str,
    lane: u8,
    tile: u16,
    x: u16,
    y: u16,
    umi: Option<&'h str>,
    comment: &'h str,
}

/// Compress headers using columnar encoding + BSC.
pub fn compress_headers_columnar(headers: &[&str]) -> Result<Vec<u8>> {
    if headers.is_empty() {
        return Ok(vec![0x00]);
    }

    // Detect format from first header
    let format = detect_format(headers[0]);
    match format {
        Some(HeaderFormat::Sra) => compress_sra(headers),
        Some(HeaderFormat::Casava) => compress_casava(headers),
        None => compress_raw_fallback(headers),
    }
}

/// Decompress columnar-encoded headers.
pub fn decompress_headers_columnar(data: &[u8], num_reads: usize) -> Result<Vec<Vec<u8>>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }

    let strings = match data[0] {
        0x00 => decompress_raw_fallback(&data[1..], num_reads),
        0x01 => decompress_sra(&data[1..], num_reads),
        0x02 => decompress_casava(&data[1..], num_reads),
        v => anyhow::bail!("Unknown header_col version: {}", v),
    }?;
    Ok(strings.into_iter().map(String::into_bytes).collect())
}

// ========================================================================
// Format detection
// ========================================================================

fn detect_format(header: &str) -> Option<HeaderFormat> {
    let without_at = header.strip_prefix('@').unwrap_or(header);
    let (name_part, comment) = without_at.split_once(' ')?;

    // SRA format: name has dots, comment has 7+ colon-separated fields
    if name_part.contains('.') {
        let rn_str = name_part.rsplit('.').next()?;
        if rn_str.parse::<u32>().is_ok() {
            let comment_base = comment.split('/').next().unwrap_or(comment);
            let fields: Vec<&str> = comment_base.split(':').collect();
            if fields.len() >= 7
                && fields[3].parse::<u8>().is_ok()
                && fields[5].parse::<u16>().is_ok()
            {
                return Some(HeaderFormat::Sra);
            }
        }
    }

    // Casava format: name has 7+ colon-separated fields with numeric lane/tile/x/y
    let name_fields: Vec<&str> = name_part.split(':').collect();
    if name_fields.len() >= 7
        && name_fields[3].parse::<u8>().is_ok()
        && name_fields[4].parse::<u16>().is_ok()
        && name_fields[5].parse::<u16>().is_ok()
        && name_fields[6].parse::<u16>().is_ok()
    {
        return Some(HeaderFormat::Casava);
    }

    None
}

// ========================================================================
// SRA format: @PREFIX.READ_NUM COMBO:LANE:TILE:X:Y[/PAIR]
// ========================================================================

fn parse_sra_row(header: &str) -> Option<SraRow<'_>> {
    let without_at = header.strip_prefix('@').unwrap_or(header);
    let (name_part, comment_part) = without_at.split_once(' ')?;
    let read_num: u32 = name_part.rsplit('.').next()?.parse().ok()?;
    let comment_base = comment_part.split('/').next().unwrap_or(comment_part);
    let fields: Vec<&str> = comment_base.split(':').collect();
    if fields.len() < 7 {
        return None;
    }
    // Zero-copy combo slice: fields[0]:fields[1]:fields[2] (no allocation)
    let combo_len = fields[0].len() + 1 + fields[1].len() + 1 + fields[2].len();
    let combo = &comment_base[..combo_len];
    Some(SraRow {
        combo,
        read_num,
        lane: fields[3].parse().ok()?,
        tile: fields[4].parse().ok()?,
        x: fields[5].parse().ok()?,
        y: fields[6].parse().ok()?,
    })
}

fn compress_sra(headers: &[&str]) -> Result<Vec<u8>> {
    let (name_prefix, pair_suffix) = detect_sra_template(headers[0]);
    let n = headers.len();

    // Phase 1 (parallel): parse all headers into Row structs.
    // combo is a zero-copy &str slice into the original header — no allocation per read.
    let rows: Vec<Option<SraRow<'_>>> = headers.par_iter().map(|h| parse_sra_row(h)).collect();

    // Phase 2 (sequential): build combo dict via HashMap (no string parsing, just equality checks).
    let mut combos: Vec<&str> = Vec::new();
    let mut combo_map: HashMap<&str, u8> = HashMap::new();
    let mut read_nums = Vec::with_capacity(n * 4);
    let mut combo_idxs = Vec::with_capacity(n);
    let mut lanes = Vec::with_capacity(n);
    let mut tiles = Vec::with_capacity(n * 2);
    let mut xs = Vec::with_capacity(n * 2);
    let mut ys = Vec::with_capacity(n * 2);

    for row_opt in &rows {
        let row = match row_opt {
            Some(r) => r,
            None => return compress_raw_fallback(headers),
        };
        let idx = if let Some(&i) = combo_map.get(row.combo) {
            i
        } else {
            if combos.len() >= 255 {
                return compress_raw_fallback(headers);
            }
            let i = combos.len() as u8;
            combos.push(row.combo);
            combo_map.insert(row.combo, i);
            i
        };
        read_nums.extend_from_slice(&row.read_num.to_le_bytes());
        combo_idxs.push(idx);
        lanes.push(row.lane);
        tiles.extend_from_slice(&row.tile.to_le_bytes());
        xs.extend_from_slice(&row.x.to_le_bytes());
        ys.extend_from_slice(&row.y.to_le_bytes());
    }

    let combos_owned: Vec<String> = combos.iter().map(|s| s.to_string()).collect();
    compress_sra_columns(&name_prefix, &pair_suffix, combos_owned, combo_idxs, read_nums, lanes, tiles, xs, ys, n)
}

fn compress_sra_columns(
    name_prefix: &str,
    pair_suffix: &str,
    combos: Vec<String>,
    combo_idxs: Vec<u8>,
    read_nums: Vec<u8>,
    lanes: Vec<u8>,
    tiles: Vec<u8>,
    xs: Vec<u8>,
    ys: Vec<u8>,
    n: usize,
) -> Result<Vec<u8>> {
    // BSC-compress all 6 columns in parallel using rayon
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
    let c_rn = c_rn?;
    let c_ci = c_ci?;
    let c_ln = c_ln?;
    let c_ti = c_ti?;
    let c_xs = c_xs?;
    let c_ys = c_ys?;

    let mut out = Vec::new();
    out.push(0x01); // version = SRA columnar
    write_string(&mut out, name_prefix);
    write_string(&mut out, pair_suffix);
    write_combo_dict(&mut out, &combos);
    out.extend_from_slice(&(n as u32).to_le_bytes());
    for col in [&c_rn, &c_ci, &c_ln, &c_ti, &c_xs, &c_ys] as [&Vec<u8>; 6] {
        write_blob(&mut out, col);
    }
    Ok(out)
}

fn decompress_sra(data: &[u8], num_reads: usize) -> Result<Vec<String>> {
    let mut pos = 0;

    let name_prefix = read_string(data, &mut pos)?;
    let pair_suffix = read_string(data, &mut pos)?;
    let combos = read_combo_dict(data, &mut pos)?;
    let stored_n = read_u32(data, &mut pos)? as usize;
    anyhow::ensure!(stored_n == num_reads, "SRA header count mismatch: stored {stored_n}, expected {num_reads}");

    // Read all 6 compressed column blobs (sequential I/O), then decompress in parallel
    let s_rn = read_compressed_slice(data, &mut pos)?;
    let s_ci = read_compressed_slice(data, &mut pos)?;
    let s_ln = read_compressed_slice(data, &mut pos)?;
    let s_ti = read_compressed_slice(data, &mut pos)?;
    let s_xs = read_compressed_slice(data, &mut pos)?;
    let s_ys = read_compressed_slice(data, &mut pos)?;

    let col_blobs: [&[u8]; 6] = [s_rn, s_ci, s_ln, s_ti, s_xs, s_ys];
    let col_results: Vec<Result<Vec<u8>>> = col_blobs.into_par_iter()
        .map(|blob| bsc::decompress_parallel(blob))
        .collect();
    let mut col_iter = col_results.into_iter();
    let d_rn = col_iter.next().unwrap()?;
    let d_ci = col_iter.next().unwrap()?;
    let d_ln = col_iter.next().unwrap()?;
    let d_ti = col_iter.next().unwrap()?;
    let d_xs = col_iter.next().unwrap()?;
    let d_ys = col_iter.next().unwrap()?;

    let headers: Result<Vec<String>> = (0..num_reads).into_par_iter().map(|i| {
        let rn = super::read_le_u32(&d_rn, i * 4)?;
        let ci = *d_ci.get(i).ok_or_else(|| anyhow::anyhow!("SRA: combo index out of bounds at read {i}"))?;
        let combo = combos.get(ci as usize).ok_or_else(|| anyhow::anyhow!("SRA: combo dict entry {} out of range at read {i}", ci))?;
        let lane = *d_ln.get(i).ok_or_else(|| anyhow::anyhow!("SRA: lane out of bounds at read {i}"))?;
        let tile = super::read_le_u16(&d_ti, i * 2)?;
        let x = super::read_le_u16(&d_xs, i * 2)?;
        let y = super::read_le_u16(&d_ys, i * 2)?;

        Ok(if pair_suffix.is_empty() {
            format!("@{}{} {}:{}:{}:{}:{}", name_prefix, rn, combo, lane, tile, x, y)
        } else {
            format!("@{}{} {}:{}:{}:{}:{}{}", name_prefix, rn, combo, lane, tile, x, y, pair_suffix)
        })
    }).collect();

    headers
}

fn detect_sra_template(header: &str) -> (String, String) {
    let without_at = header.strip_prefix('@').unwrap_or(header);
    let name_part = without_at.split(' ').next().unwrap_or(without_at);

    let name_prefix = if let Some(dot_pos) = name_part.rfind('.') {
        name_part[..=dot_pos].to_string()
    } else {
        String::new()
    };

    let comment = without_at.split(' ').nth(1).unwrap_or("");
    let pair_suffix = if let Some(slash_pos) = comment.rfind('/') {
        comment[slash_pos..].to_string()
    } else {
        String::new()
    };

    (name_prefix, pair_suffix)
}

// ========================================================================
// Casava 1.8 format: @INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y COMMENT
// ========================================================================

fn parse_casava_row(header: &str) -> Option<CasavaRow<'_>> {
    let without_at = header.strip_prefix('@').unwrap_or(header);
    let (name_part, comment) = without_at.split_once(' ')?;
    let fields: Vec<&str> = name_part.split(':').collect();
    if fields.len() < 7 {
        return None;
    }
    // Zero-copy combo slice: fields[0]:fields[1]:fields[2] (no allocation)
    let combo_len = fields[0].len() + 1 + fields[1].len() + 1 + fields[2].len();
    let combo = &name_part[..combo_len];
    let lane: u8 = fields[3].parse().ok()?;
    let tile: u16 = fields[4].parse().ok()?;
    let x: u16 = fields[5].parse().ok()?;
    let y: u16 = fields[6].parse().ok()?;
    // UMI: everything after field[6] (may be multi-part with ':')
    let umi = if fields.len() >= 8 {
        let no_umi_len = combo_len
            + 1 + fields[3].len()
            + 1 + fields[4].len()
            + 1 + fields[5].len()
            + 1 + fields[6].len();
        Some(&name_part[no_umi_len + 1..])
    } else {
        None
    };
    Some(CasavaRow { combo, lane, tile, x, y, umi, comment })
}

fn compress_casava(headers: &[&str]) -> Result<Vec<u8>> {
    let n = headers.len();

    // Phase 1 (parallel): parse all headers into Row structs.
    // combo, umi, and comment are zero-copy &str slices into the original headers.
    let rows: Vec<Option<CasavaRow<'_>>> = headers.par_iter().map(|h| parse_casava_row(h)).collect();

    // Phase 2 (sequential): build combo dict, detect common comment, fill columns.
    let first_comment = rows[0].as_ref().map_or("", |r| r.comment);
    let mut all_same_comment = true;
    let mut has_umi = false;
    let mut combos: Vec<&str> = Vec::new();
    let mut combo_map: HashMap<&str, u8> = HashMap::new();
    let mut combo_idxs = Vec::with_capacity(n);
    let mut lanes = Vec::with_capacity(n);
    let mut tiles = Vec::with_capacity(n * 2);
    let mut xs = Vec::with_capacity(n * 2);
    let mut ys = Vec::with_capacity(n * 2);

    for row_opt in &rows {
        let row = match row_opt {
            Some(r) => r,
            None => return compress_raw_fallback(headers),
        };
        if row.comment != first_comment {
            all_same_comment = false;
        }
        if row.umi.is_some() {
            has_umi = true;
        }
        let idx = if let Some(&i) = combo_map.get(row.combo) {
            i
        } else {
            if combos.len() >= 255 {
                return compress_raw_fallback(headers);
            }
            let i = combos.len() as u8;
            combos.push(row.combo);
            combo_map.insert(row.combo, i);
            i
        };
        combo_idxs.push(idx);
        lanes.push(row.lane);
        tiles.extend_from_slice(&row.tile.to_le_bytes());
        xs.extend_from_slice(&row.x.to_le_bytes());
        ys.extend_from_slice(&row.y.to_le_bytes());
    }

    // Build per-read comment and UMI byte columns using already-parsed slices.
    let mut comments: Vec<u8> = Vec::new();
    let mut umis: Vec<u8> = Vec::new();
    if !all_same_comment {
        for row_opt in &rows {
            let row = row_opt.as_ref().unwrap();
            write_varint(&mut comments, row.comment.len());
            comments.extend_from_slice(row.comment.as_bytes());
        }
    }
    if has_umi {
        for row_opt in &rows {
            let row = row_opt.as_ref().unwrap();
            let umi = row.umi.unwrap_or("");
            write_varint(&mut umis, umi.len());
            umis.extend_from_slice(umi.as_bytes());
        }
    }

    // BSC-compress all 5 spatial columns in parallel
    let (c_ci, (c_ln, (c_ti, (c_xs, c_ys)))) = rayon::join(
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
    );
    let c_ci = c_ci?;
    let c_ln = c_ln?;
    let c_ti = c_ti?;
    let c_xs = c_xs?;
    let c_ys = c_ys?;

    // Build output
    let mut out = Vec::new();
    out.push(0x02); // version = Casava columnar
    let combos_owned: Vec<String> = combos.iter().map(|s| s.to_string()).collect();
    write_combo_dict(&mut out, &combos_owned);

    // Common comment flag + data
    if all_same_comment {
        out.push(0x01);
        write_string(&mut out, first_comment);
    } else {
        out.push(0x00);
    }

    // UMI flag
    out.push(if has_umi { 0x01 } else { 0x00 });

    out.extend_from_slice(&(n as u32).to_le_bytes());

    // 5 spatial columns
    for col in [&c_ci, &c_ln, &c_ti, &c_xs, &c_ys] as [&Vec<u8>; 5] {
        write_blob(&mut out, col);
    }

    // Comment column (only if not common)
    if !all_same_comment {
        let c_comments = bsc::compress_parallel_adaptive(&comments)?;
        write_blob(&mut out, &c_comments);
    }

    // UMI column (only if present)
    if has_umi {
        let c_umis = bsc::compress_parallel_adaptive(&umis)?;
        write_blob(&mut out, &c_umis);
    }

    Ok(out)
}

fn decompress_casava(data: &[u8], num_reads: usize) -> Result<Vec<String>> {
    let mut pos = 0;

    let combos = read_combo_dict(data, &mut pos)?;

    // Common comment
    if pos >= data.len() {
        anyhow::bail!("Casava header: truncated at common comment flag");
    }
    let has_common = data[pos];
    pos += 1;
    let common_comment = if has_common == 0x01 {
        Some(read_string(data, &mut pos)?)
    } else {
        None
    };

    // UMI flag
    if pos >= data.len() {
        anyhow::bail!("Casava header: truncated at UMI flag");
    }
    let has_umi = data[pos] == 0x01;
    pos += 1;

    let stored_n = read_u32(data, &mut pos)? as usize;
    anyhow::ensure!(stored_n == num_reads, "Casava header count mismatch: stored {stored_n}, expected {num_reads}");

    // Decompress 5 spatial columns
    let d_ci = read_and_decompress(data, &mut pos)?;
    let d_ln = read_and_decompress(data, &mut pos)?;
    let d_ti = read_and_decompress(data, &mut pos)?;
    let d_xs = read_and_decompress(data, &mut pos)?;
    let d_ys = read_and_decompress(data, &mut pos)?;

    // Comment column
    let comment_data = if common_comment.is_none() {
        Some(read_and_decompress(data, &mut pos)?)
    } else {
        None
    };

    // UMI column
    let umi_data = if has_umi {
        Some(read_and_decompress(data, &mut pos)?)
    } else {
        None
    };

    // Pre-scan variable-length comment offsets (sequential, needed before parallel phase)
    let comment_offsets: Option<Vec<(usize, usize)>> = if let Some(cd) = &comment_data {
        let mut offsets = Vec::with_capacity(num_reads);
        let mut off = 0usize;
        for i in 0..num_reads {
            let len = read_varint_from_slice(cd, &mut off)?;
            let end = off.checked_add(len)
                .filter(|&e| e <= cd.len())
                .ok_or_else(|| anyhow::anyhow!("header_col: truncated comment data for header {i}"))?;
            offsets.push((off, end));
            off = end;
        }
        Some(offsets)
    } else {
        None
    };

    // Pre-scan variable-length UMI offsets
    let umi_offsets: Option<Vec<(usize, usize)>> = if let Some(ud) = &umi_data {
        let mut offsets = Vec::with_capacity(num_reads);
        let mut off = 0usize;
        for i in 0..num_reads {
            let len = read_varint_from_slice(ud, &mut off)?;
            let end = off.checked_add(len)
                .filter(|&e| e <= ud.len())
                .ok_or_else(|| anyhow::anyhow!("header_col: truncated UMI data for header {i}"))?;
            offsets.push((off, end));
            off = end;
        }
        Some(offsets)
    } else {
        None
    };

    let headers: Result<Vec<String>> = (0..num_reads).into_par_iter().map(|i| {
        let ci = *d_ci.get(i).ok_or_else(|| anyhow::anyhow!("Casava: combo index out of bounds at read {i}"))?;
        let combo = combos.get(ci as usize).ok_or_else(|| anyhow::anyhow!("Casava: combo dict entry {} out of range at read {i}", ci))?;
        let lane = *d_ln.get(i).ok_or_else(|| anyhow::anyhow!("Casava: lane out of bounds at read {i}"))?;
        let tile = super::read_le_u16(&d_ti, i * 2)?;
        let x = super::read_le_u16(&d_xs, i * 2)?;
        let y = super::read_le_u16(&d_ys, i * 2)?;

        let comment = match &common_comment {
            Some(cc) => cc.as_str(),
            None => {
                let cd = comment_data.as_ref().unwrap();
                let (s, e) = comment_offsets.as_ref().unwrap()[i];
                std::str::from_utf8(&cd[s..e]).context("invalid comment UTF-8")?
            }
        };

        Ok(if has_umi {
            let ud = umi_data.as_ref().unwrap();
            let (s, e) = umi_offsets.as_ref().unwrap()[i];
            let umi = std::str::from_utf8(&ud[s..e]).context("invalid UMI UTF-8")?;
            format!("@{}:{}:{}:{}:{}:{} {}", combo, lane, tile, x, y, umi, comment)
        } else {
            format!("@{}:{}:{}:{}:{} {}", combo, lane, tile, x, y, comment)
        })
    }).collect();

    headers
}

// ========================================================================
// Shared helpers
// ========================================================================

fn write_string(out: &mut Vec<u8>, s: &str) {
    out.extend_from_slice(&(s.len() as u16).to_le_bytes());
    out.extend_from_slice(s.as_bytes());
}

fn read_string(data: &[u8], pos: &mut usize) -> Result<String> {
    if *pos + 2 > data.len() {
        anyhow::bail!("header_col: truncated string length at offset {}", *pos);
    }
    let len = u16::from_le_bytes([data[*pos], data[*pos + 1]]) as usize;
    *pos += 2;
    if *pos + len > data.len() {
        anyhow::bail!("header_col: truncated string data at offset {}", *pos);
    }
    let s = std::str::from_utf8(&data[*pos..*pos + len])
        .context("invalid string")?
        .to_string();
    *pos += len;
    Ok(s)
}

fn write_combo_dict(out: &mut Vec<u8>, combos: &[String]) {
    // combos.len() <= 254 is guaranteed by get_or_insert_combo's >= 255 guard
    out.push(combos.len() as u8);
    for combo in combos {
        // Individual combo strings must fit in one byte length prefix.
        // Combos are constructed from short SAM tag values; > 255 bytes is a bug.
        let combo_len = u8::try_from(combo.len())
            .expect("combo string exceeds 255 bytes; should have been rejected earlier");
        out.push(combo_len);
        out.extend_from_slice(combo.as_bytes());
    }
}

fn read_combo_dict(data: &[u8], pos: &mut usize) -> Result<Vec<String>> {
    if *pos >= data.len() {
        anyhow::bail!("header_col: truncated combo dict at offset {}", *pos);
    }
    let num_combos = data[*pos] as usize;
    *pos += 1;
    let mut combos = Vec::with_capacity(num_combos);
    for i in 0..num_combos {
        if *pos >= data.len() {
            anyhow::bail!("header_col: truncated combo {i} length at offset {}", *pos);
        }
        let combo_len = data[*pos] as usize;
        *pos += 1;
        if *pos + combo_len > data.len() {
            anyhow::bail!("header_col: truncated combo {i} data at offset {}", *pos);
        }
        let combo = std::str::from_utf8(&data[*pos..*pos + combo_len])
            .context("invalid combo string")?
            .to_string();
        combos.push(combo);
        *pos += combo_len;
    }
    Ok(combos)
}

fn write_blob(out: &mut Vec<u8>, blob: &[u8]) {
    out.extend_from_slice(&(blob.len() as u32).to_le_bytes());
    out.extend_from_slice(blob);
}

fn read_and_decompress(data: &[u8], pos: &mut usize) -> Result<Vec<u8>> {
    let compressed = read_compressed_slice(data, pos)?;
    bsc::decompress_parallel(compressed)
}

/// Read a length-prefixed BSC blob and return a slice without decompressing.
fn read_compressed_slice<'a>(data: &'a [u8], pos: &mut usize) -> Result<&'a [u8]> {
    let len = read_u32(data, pos)? as usize;
    if *pos + len > data.len() {
        anyhow::bail!("header_col: truncated compressed block at offset {}", *pos);
    }
    let slice = &data[*pos..*pos + len];
    *pos += len;
    Ok(slice)
}

fn read_u32(data: &[u8], pos: &mut usize) -> Result<u32> {
    if *pos + 4 > data.len() {
        anyhow::bail!("header_col: truncated u32 at offset {}", *pos);
    }
    let v = u32::from_le_bytes([data[*pos], data[*pos+1], data[*pos+2], data[*pos+3]]);
    *pos += 4;
    Ok(v)
}

fn write_varint(out: &mut Vec<u8>, mut value: usize) {
    while value >= 0x80 {
        out.push(((value & 0x7F) | 0x80) as u8);
        value >>= 7;
    }
    out.push(value as u8);
}

fn read_varint_from_slice(data: &[u8], offset: &mut usize) -> Result<usize> {
    let mut value = 0usize;
    let mut shift = 0;
    loop {
        if *offset >= data.len() {
            anyhow::bail!("header_col: truncated varint at offset {}", *offset);
        }
        let byte = data[*offset];
        *offset += 1;
        value |= ((byte & 0x7F) as usize) << shift;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
        shift += 7;
    }
}

// ========================================================================
// Raw BSC fallback
// ========================================================================

fn compress_raw_fallback(headers: &[&str]) -> Result<Vec<u8>> {
    let mut stream = Vec::new();
    for h in headers {
        write_varint(&mut stream, h.len());
        stream.extend_from_slice(h.as_bytes());
    }
    let compressed = bsc::compress_parallel_adaptive(&stream)?;

    let mut out = Vec::new();
    out.push(0x00);
    out.write_all(&compressed)?;
    Ok(out)
}

fn decompress_raw_fallback(data: &[u8], num_reads: usize) -> Result<Vec<String>> {
    let decompressed = bsc::decompress_parallel(data)?;
    let mut headers = Vec::with_capacity(num_reads);
    let mut offset = 0;

    for i in 0..num_reads {
        let len = read_varint_from_slice(&decompressed, &mut offset)?;
        let end = offset.checked_add(len)
            .filter(|&e| e <= decompressed.len())
            .ok_or_else(|| anyhow::anyhow!("header_col: truncated header data for read {i}"))?;
        let header = std::str::from_utf8(&decompressed[offset..end])
            .context("invalid header UTF-8")?
            .to_string();
        offset = end;
        headers.push(header);
    }

    Ok(headers)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip_sra() {
        let headers = vec![
            "@ERR3239334.100 A00296:43:HCLHLDSXX:4:2174:25834:28823/1",
            "@ERR3239334.200 A00296:43:HCLHLDSXX:4:2174:25835:28824/1",
            "@ERR3239334.300 A00217:77:HFJWFDSXX:2:1101:10000:5000/1",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x01); // SRA format
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_roundtrip_sra_no_pair() {
        let headers = vec![
            "@ERR3239334.100 A00296:43:HCLHLDSXX:4:2174:25834:28823",
            "@ERR3239334.200 A00296:43:HCLHLDSXX:4:2174:25835:28824",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x01);
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_roundtrip_casava() {
        let headers = vec![
            "@HWI-D00119:50:H7AP8ADXX:1:1101:1213:2058 1:N:0:TAAGGCGA",
            "@HWI-D00119:50:H7AP8ADXX:1:1101:1214:2059 1:N:0:TAAGGCGA",
            "@HWI-D00119:50:H7AP8ADXX:2:1201:5000:6000 1:N:0:TAAGGCGA",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x02); // Casava format
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_roundtrip_casava_varying_comments() {
        let headers = vec![
            "@HWI-D00119:50:H7AP8ADXX:1:1101:1213:2058 1:N:0:TAAGGCGA",
            "@HWI-D00119:50:H7AP8ADXX:1:1101:1214:2059 1:N:0:CTTGTAAT",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x02);
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_roundtrip_dragen_umi() {
        // DRAGEN BCL Convert format with UMI as 8th field
        let headers = vec![
            "@SIM:1:FCX:1:2106:15337:1063:GATCTGTACGTC 1:N:0:ATCACG",
            "@SIM:1:FCX:1:2106:15338:1064:AACCTGTACGTC 1:N:0:ATCACG",
            "@SIM:1:FCX:2:1205:8000:9000:TTGGATCAACGT 1:N:0:ATCACG",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x02); // Casava format
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_roundtrip_dragen_umi_varying_comments() {
        // DRAGEN with UMI and varying barcodes
        let headers = vec![
            "@SIM:1:FCX:1:2106:15337:1063:GATCTGTACGTC 1:N:0:ATCACG",
            "@SIM:1:FCX:1:2106:15338:1064:AACCTGTACGTC 1:N:0:CTTGTAAT",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x02);
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_roundtrip_dragen_umi_reverse_complement() {
        // DRAGEN with 'r' prefix on UMI (reverse complement flag)
        let headers = vec![
            "@SIM:1:FCX:1:2106:15337:1063:rGATCTGTACGTC 1:N:0:ATCACG",
            "@SIM:1:FCX:1:2106:15338:1064:rAACCTGTACGTC 1:N:0:ATCACG",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x02);
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_fallback_non_illumina() {
        let headers = vec![
            "@READ_1 some comment",
            "@READ_2 another comment",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x00); // raw fallback
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_empty() {
        let headers: Vec<&str> = vec![];
        let compressed = compress_headers_columnar(&headers).unwrap();
        let decompressed = decompress_headers_columnar(&compressed, 0).unwrap();
        assert!(decompressed.is_empty());
    }
}
