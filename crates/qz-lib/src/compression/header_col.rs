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

/// Parsed SRA header fields (borrows from original header — zero allocation).
///
/// Spatial fields are u32 to support newer instruments (NovaSeq X x/y can reach
/// ~70k, which overflows u16). The on-disk layout still uses u16 (v0x01) when
/// all values fit, and switches to u32 (v0x04) only when they don't, keeping
/// archive size stable for typical Illumina runs.
struct SraRow<'h> {
    combo: &'h str,
    /// Name part up to and including the final '.', e.g. "ERR3239334." — stored
    /// once per archive, so every row must share the same prefix.
    prefix: &'h str,
    /// Trailing "/1", "/2", … on the comment (or ""), stored once per archive,
    /// so every row must share the same suffix.
    pair_suffix: &'h str,
    read_num: u32,
    lane: u8,
    tile: u32,
    x: u32,
    y: u32,
}

/// True iff `s` is the canonical decimal rendering of its integer value — i.e.
/// parsing `s` to an integer and re-emitting it with `Display` reproduces `s`
/// byte-for-byte. Rejects empty strings, a leading '+' ("+5"), and leading
/// zeros ("007"), all of which the columnar encoders would otherwise silently
/// normalize away, corrupting the header on decompress.
#[inline]
fn is_canonical_decimal(s: &str) -> bool {
    !s.is_empty()
        && s.bytes().all(|b| b.is_ascii_digit())
        && (s.len() == 1 || s.as_bytes()[0] != b'0')
}

/// Parsed Casava header fields (borrows from original header — zero allocation).
/// See [`SraRow`] for the u32 rationale and the dual-version (0x02 / 0x03) layout.
struct CasavaRow<'h> {
    combo: &'h str,
    lane: u8,
    tile: u32,
    x: u32,
    y: u32,
    umi: Option<&'h str>,
    comment: &'h str,
}

/// Compress headers (raw bytes) using columnar encoding + BSC.
///
/// FASTQ IDs are arbitrary bytes, but the structured SRA/Casava encoders are
/// only meaningful for valid Illumina headers, which are always ASCII. So:
/// if every header is valid UTF-8 we run the structured/raw String pipeline
/// (byte-exact, since valid UTF-8 round-trips losslessly); otherwise we use the
/// byte-exact raw fallback directly. This preserves non-UTF-8 IDs that would
/// previously be mangled by a lossy UTF-8 conversion at the call site.
pub fn compress_headers_columnar_bytes(headers: &[&[u8]]) -> Result<Vec<u8>> {
    if headers.is_empty() {
        return Ok(vec![0x00]);
    }
    match headers
        .iter()
        .map(|h| std::str::from_utf8(h).ok())
        .collect::<Option<Vec<&str>>>()
    {
        Some(strs) => compress_headers_columnar(&strs),
        None => compress_raw_fallback_bytes(headers),
    }
}

/// Compress headers using columnar encoding + BSC.
///
/// Requires valid UTF-8 headers (callers with arbitrary bytes should use
/// [`compress_headers_columnar_bytes`]).
pub fn compress_headers_columnar(headers: &[&str]) -> Result<Vec<u8>> {
    if headers.is_empty() {
        return Ok(vec![0x00]);
    }

    // Tier 1: format-specific (best ratio for the common case).
    let specific = match detect_format(headers[0]) {
        Some(HeaderFormat::Sra) => compress_sra(headers)?,
        Some(HeaderFormat::Casava) => compress_casava(headers)?,
        None => None,
    };
    if let Some(blob) = specific {
        return Ok(blob);
    }

    // Tier 2: generic token-columnar, size-gated against raw (never worse than raw).
    // The raw fallback is needed in BOTH outcomes — as the size-gate baseline when
    // generic succeeds, and as the result when generic declines — so it is computed
    // unconditionally either way. Run it in parallel with the generic encode (both
    // are internally rayon-parallel BSC work) so the size-gate costs no extra
    // wall-clock: the header stage is gated by max(generic, raw), not their sum.
    let (generic, raw) = rayon::join(
        || super::header_token::compress_generic_token(headers),
        || compress_raw_fallback(headers),
    );
    let raw = raw?;
    match generic? {
        // Strict guarantee preserved: emit generic only when it is not larger.
        Some(blob) if blob.len() <= raw.len() => Ok(blob),
        // Tier 3: byte-exact raw fallback (generic declined or lost the size-gate).
        _ => Ok(raw),
    }
}

/// Decompress columnar-encoded headers.
pub fn decompress_headers_columnar(data: &[u8], num_reads: usize) -> Result<Vec<Vec<u8>>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }

    // 0x01/0x02: legacy u16 spatial fields. 0x03/0x04: u32 spatial fields for
    // NovaSeq X-scale x/y. Both widths share the same column-blob layout.
    // 0x00 (raw fallback) decodes directly to bytes so non-UTF-8 headers
    // survive; the structured paths only ever encode ASCII Illumina headers.
    let bytes: Vec<Vec<u8>> = match data[0] {
        0x00 => decompress_raw_fallback(&data[1..], num_reads)?,
        0x01 => string_vec_into_bytes(decompress_sra(&data[1..], num_reads, false)?),
        0x04 => string_vec_into_bytes(decompress_sra(&data[1..], num_reads, true)?),
        0x02 => string_vec_into_bytes(decompress_casava(&data[1..], num_reads, false)?),
        0x03 => string_vec_into_bytes(decompress_casava(&data[1..], num_reads, true)?),
        0x05 => super::header_token::decompress_generic_token(&data[1..], num_reads)?,
        v => anyhow::bail!("Unknown header_col version: {}", v),
    };
    Ok(bytes)
}

#[inline]
fn string_vec_into_bytes(strings: Vec<String>) -> Vec<Vec<u8>> {
    strings.into_iter().map(String::into_bytes).collect()
}

/// Decompress columnar-encoded headers directly into a varint-prefixed byte
/// stream `[varint(len)][bytes]*` — the format the bounded-streaming consumer
/// (`StreamCursor::read_one_record_varint`) expects.
///
/// Fast path: for SRA (v0x01 / v0x04) the per-record `format!()` + `String`
/// allocation chain in [`decompress_sra`] is replaced by chunked-parallel
/// direct byte writes. ~8K records per chunk worker, then sequential concat —
/// 5M allocations become ~640. This avoids hammering glibc's arena mutex
/// (which serialized all 72 rayon workers and burned 43% of decompression
/// cycles on `native_queued_spin_lock_slowpath` before).
///
/// Other formats fall back to [`decompress_headers_columnar`] + a serial
/// varint re-emit, which still re-allocates per record but is rarely on the
/// hot path for our datasets.
pub fn decompress_headers_columnar_to_varint_stream(
    data: &[u8],
    num_reads: usize,
) -> Result<Vec<u8>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }
    match data[0] {
        0x01 => decompress_sra_into_varint(&data[1..], num_reads, false),
        0x04 => decompress_sra_into_varint(&data[1..], num_reads, true),
        0x02 => decompress_casava_into_varint(&data[1..], num_reads, false),
        0x03 => decompress_casava_into_varint(&data[1..], num_reads, true),
        0x05 => super::header_token::decompress_generic_token_to_varint(&data[1..], num_reads),
        _ => {
            // Slow fallback for 0x00 raw-BSC and any future versions:
            // materialize as Vec<Vec<u8>>, then re-emit.
            let headers = decompress_headers_columnar(data, num_reads)?;
            let total: usize = headers.iter().map(|h| h.len() + 5).sum();
            let mut out = Vec::with_capacity(total);
            for h in &headers {
                write_varint(&mut out, h.len());
                out.extend_from_slice(h);
            }
            Ok(out)
        }
    }
}

// ========================================================================
// Format detection
// ========================================================================

fn detect_format(header: &str) -> Option<HeaderFormat> {
    // Structured SRA/Casava encoding only applies to FASTQ headers, which always
    // begin with '@'. A header lacking it (e.g. a '>'-prefixed FASTA header that
    // happens to look Illumina-shaped) must fall through to the byte-exact raw
    // fallback — the structured decoders re-add a literal '@' and would corrupt
    // '>id' into '@>id'.
    let without_at = header.strip_prefix('@')?;
    let (name_part, comment) = without_at.split_once(' ')?;

    // SRA format: name has dots, comment has 7+ colon-separated fields.
    // u32 parses accept NovaSeq X-scale x/y values (which would overflow u16).
    if name_part.contains('.') {
        let rn_str = name_part.rsplit('.').next()?;
        if rn_str.parse::<u32>().is_ok() {
            let comment_base = comment.split('/').next().unwrap_or(comment);
            let fields: Vec<&str> = comment_base.split(':').collect();
            if fields.len() >= 7
                && fields[3].parse::<u8>().is_ok()
                && fields[5].parse::<u32>().is_ok()
            {
                return Some(HeaderFormat::Sra);
            }
        }
    }

    // Casava format: name has 7+ colon-separated fields with numeric lane/tile/x/y
    let name_fields: Vec<&str> = name_part.split(':').collect();
    if name_fields.len() >= 7
        && name_fields[3].parse::<u8>().is_ok()
        && name_fields[4].parse::<u32>().is_ok()
        && name_fields[5].parse::<u32>().is_ok()
        && name_fields[6].parse::<u32>().is_ok()
    {
        return Some(HeaderFormat::Casava);
    }

    None
}

// ========================================================================
// SRA format: @PREFIX.READ_NUM COMBO:LANE:TILE:X:Y[/PAIR]
// ========================================================================

fn parse_sra_row(header: &str) -> Option<SraRow<'_>> {
    // Require the FASTQ '@' prefix; the decoder re-adds a literal '@', so a
    // header without it (e.g. '>'-prefixed FASTA) must not be parsed here.
    let without_at = header.strip_prefix('@')?;
    let (name_part, comment_part) = without_at.split_once(' ')?;

    // Split the name into "<prefix>." + "<read_num>" at the final dot. Both the
    // prefix and the numeric rendering must round-trip exactly.
    let dot = name_part.rfind('.')?;
    let prefix = &name_part[..=dot];
    let rn_str = &name_part[dot + 1..];
    if !is_canonical_decimal(rn_str) {
        return None;
    }
    let read_num: u32 = rn_str.parse().ok()?;

    // The comment is "<combo:lane:tile:x:y>" optionally followed by a pair
    // suffix beginning at the first '/'. Storing the suffix verbatim and
    // requiring `comment_base + pair_suffix == comment_part` keeps it exact.
    let (comment_base, pair_suffix) = match comment_part.find('/') {
        Some(p) => (&comment_part[..p], &comment_part[p..]),
        None => (comment_part, ""),
    };
    let fields: Vec<&str> = comment_base.split(':').collect();
    // Exactly 7 fields (combo's 3 + lane/tile/x/y). `>= 7` would silently drop
    // any trailing colon-separated field on reconstruction.
    if fields.len() != 7 {
        return None;
    }
    // Numeric fields must render canonically or a parse-then-format round trip
    // would normalize zero-padding away.
    if !is_canonical_decimal(fields[3])
        || !is_canonical_decimal(fields[4])
        || !is_canonical_decimal(fields[5])
        || !is_canonical_decimal(fields[6])
    {
        return None;
    }
    // Zero-copy combo slice: fields[0]:fields[1]:fields[2] (no allocation)
    let combo_len = fields[0].len() + 1 + fields[1].len() + 1 + fields[2].len();
    let combo = &comment_base[..combo_len];
    Some(SraRow {
        combo,
        prefix,
        pair_suffix,
        read_num,
        lane: fields[3].parse().ok()?,
        tile: fields[4].parse().ok()?,
        x: fields[5].parse().ok()?,
        y: fields[6].parse().ok()?,
    })
}

/// Pick the on-disk version for the spatial columns.
/// Returns the width byte (`0x01` SRA v0x01 / `0x02` Casava v0x02 → u16 fields;
/// `0x04` SRA v0x04 / `0x03` Casava v0x03 → u32 fields). Encoders pick the
/// narrower version whenever values fit, so existing Illumina datasets stay
/// byte-identical to pre-widening output.
///
/// Empty input returns `true` (vacuously — no value overflows). Callers can
/// safely choose the narrow format for a 0-read archive: no spatial bytes are
/// written either way, and decoders never read spatial data with num_reads=0.
fn fits_in_u16(values: &[u32]) -> bool {
    values.iter().all(|&v| v <= u16::MAX as u32)
}

fn compress_sra(headers: &[&str]) -> Result<Option<Vec<u8>>> {
    let n = headers.len();

    // Phase 1 (parallel): parse all headers into Row structs.
    // combo is a zero-copy &str slice into the original header — no allocation per read.
    let rows: Vec<Option<SraRow<'_>>> = headers.par_iter().map(|h| parse_sra_row(h)).collect();

    // The name prefix and pair suffix are stored once for the whole archive, so
    // every read must share read 0's. Interleaved /1+/2 mates and merged
    // multi-accession files violate this — fall back to raw rather than
    // reconstruct every read with read 0's prefix/suffix (silent corruption).
    let (name_prefix, pair_suffix) = match rows[0].as_ref() {
        Some(r) => (r.prefix, r.pair_suffix),
        None => return Ok(None),
    };
    if rows.iter().any(|r| {
        r.as_ref()
            .is_none_or(|r| r.prefix != name_prefix || r.pair_suffix != pair_suffix)
    }) {
        return Ok(None);
    }

    // `write_string` length-prefixes with a u16 (see its definition). A name prefix or
    // pair suffix longer than `u16::MAX` would truncate the stored length while writing
    // the full bytes, yielding a blob this very encoder cannot decode. Such a giant
    // shared prefix is pathological but valid input — fall back to the byte-exact raw path.
    if name_prefix.len() > u16::MAX as usize || pair_suffix.len() > u16::MAX as usize {
        return Ok(None);
    }

    // Phase 2 (sequential): build combo dict via HashMap (no string parsing, just equality checks).
    let mut combos: Vec<&str> = Vec::new();
    let mut combo_map: HashMap<&str, u8> = HashMap::new();
    let mut read_nums = Vec::with_capacity(n * 4);
    let mut combo_idxs = Vec::with_capacity(n);
    let mut lanes = Vec::with_capacity(n);
    let mut tiles = Vec::with_capacity(n * 2);
    let mut xs = Vec::with_capacity(n * 2);
    let mut ys = Vec::with_capacity(n * 2);

    // Collect spatial fields as u32 first; the on-disk width is decided after
    // the full scan based on whether any value overflows u16.
    let mut tiles_u32 = Vec::with_capacity(n);
    let mut xs_u32 = Vec::with_capacity(n);
    let mut ys_u32 = Vec::with_capacity(n);

    for row_opt in &rows {
        let row = match row_opt {
            Some(r) => r,
            None => return Ok(None),
        };
        let idx = if let Some(&i) = combo_map.get(row.combo) {
            i
        } else {
            // >= 255 distinct combos won't fit the dict's u8 index, and a single
            // combo > 255 bytes won't fit its u8 length prefix; either way fall
            // back to raw rather than panicking in write_combo_dict.
            if combos.len() >= 255 || row.combo.len() > u8::MAX as usize {
                return Ok(None);
            }
            let i = combos.len() as u8;
            combos.push(row.combo);
            combo_map.insert(row.combo, i);
            i
        };
        read_nums.extend_from_slice(&row.read_num.to_le_bytes());
        combo_idxs.push(idx);
        lanes.push(row.lane);
        tiles_u32.push(row.tile);
        xs_u32.push(row.x);
        ys_u32.push(row.y);
    }

    // Pick the narrowest format that fits. Preserves byte-identical output for
    // typical Illumina headers; only NovaSeq X-scale x/y triggers v0x04.
    let narrow = fits_in_u16(&tiles_u32) && fits_in_u16(&xs_u32) && fits_in_u16(&ys_u32);
    if narrow {
        for v in &tiles_u32 {
            tiles.extend_from_slice(&(*v as u16).to_le_bytes());
        }
        for v in &xs_u32 {
            xs.extend_from_slice(&(*v as u16).to_le_bytes());
        }
        for v in &ys_u32 {
            ys.extend_from_slice(&(*v as u16).to_le_bytes());
        }
    } else {
        for v in &tiles_u32 {
            tiles.extend_from_slice(&v.to_le_bytes());
        }
        for v in &xs_u32 {
            xs.extend_from_slice(&v.to_le_bytes());
        }
        for v in &ys_u32 {
            ys.extend_from_slice(&v.to_le_bytes());
        }
    }

    let combos_owned: Vec<String> = combos.iter().map(|s| s.to_string()).collect();
    let version = if narrow { 0x01u8 } else { 0x04u8 };

    // Verify the structured encoding round-trips byte-exact BEFORE paying for
    // BSC. We reconstruct from the in-memory columns (byte-identical layout to
    // the BSC-decoded columns) and compare to the input — this catches any field
    // parse/reconstruct discrepancy, the only way the columnar encoding can be
    // lossy. It deliberately skips the BSC round-trip: BSC is deterministically
    // lossless and CRC32-guarded on disk, and is already trusted un-verified for
    // the much larger sequence/quality streams. On any mismatch (or reconstruct
    // error), fall back to the byte-exact raw encoding. Replaces the old
    // `verified_or_raw`, which decoded the full BSC stream and dominated the
    // header compress stage.
    let verified = reconstruct_sra_headers(
        name_prefix,
        pair_suffix,
        &combos_owned,
        &read_nums,
        &combo_idxs,
        &lanes,
        &tiles,
        &xs,
        &ys,
        !narrow,
        n,
    )
    .map(|recon| {
        recon.len() == n
            && recon
                .iter()
                .zip(headers)
                .all(|(r, h)| r.as_bytes() == h.as_bytes())
    })
    .unwrap_or(false);
    if !verified {
        return Ok(None);
    }

    compress_sra_columns(
        version,
        name_prefix,
        pair_suffix,
        combos_owned,
        combo_idxs,
        read_nums,
        lanes,
        tiles,
        xs,
        ys,
        n,
    )
    .map(Some)
}

fn compress_sra_columns(
    version: u8,
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
        || {
            rayon::join(
                || bsc::compress_parallel_adaptive(&combo_idxs),
                || {
                    rayon::join(
                        || bsc::compress_parallel_adaptive(&lanes),
                        || {
                            rayon::join(
                                || bsc::compress_parallel_adaptive(&tiles),
                                || {
                                    rayon::join(
                                        || bsc::compress_parallel_adaptive(&xs),
                                        || bsc::compress_parallel_adaptive(&ys),
                                    )
                                },
                            )
                        },
                    )
                },
            )
        },
    );
    let c_rn = c_rn?;
    let c_ci = c_ci?;
    let c_ln = c_ln?;
    let c_ti = c_ti?;
    let c_xs = c_xs?;
    let c_ys = c_ys?;

    let mut out = Vec::new();
    out.push(version); // 0x01 (u16 tile/x/y) or 0x04 (u32 tile/x/y)
    write_string(&mut out, name_prefix);
    write_string(&mut out, pair_suffix);
    write_combo_dict(&mut out, &combos);
    out.extend_from_slice(&(n as u32).to_le_bytes());
    for col in [&c_rn, &c_ci, &c_ln, &c_ti, &c_xs, &c_ys] as [&Vec<u8>; 6] {
        write_blob(&mut out, col);
    }
    Ok(out)
}

fn decompress_sra(data: &[u8], num_reads: usize, wide: bool) -> Result<Vec<String>> {
    let mut pos = 0;

    let name_prefix = read_string(data, &mut pos)?;
    let pair_suffix = read_string(data, &mut pos)?;
    let combos = read_combo_dict(data, &mut pos)?;
    let stored_n = read_u32(data, &mut pos)? as usize;
    anyhow::ensure!(
        stored_n == num_reads,
        "SRA header count mismatch: stored {stored_n}, expected {num_reads}"
    );

    // Read all 6 compressed column blobs (sequential I/O), then decompress in parallel
    let s_rn = read_compressed_slice(data, &mut pos)?;
    let s_ci = read_compressed_slice(data, &mut pos)?;
    let s_ln = read_compressed_slice(data, &mut pos)?;
    let s_ti = read_compressed_slice(data, &mut pos)?;
    let s_xs = read_compressed_slice(data, &mut pos)?;
    let s_ys = read_compressed_slice(data, &mut pos)?;

    let col_blobs: [&[u8]; 6] = [s_rn, s_ci, s_ln, s_ti, s_xs, s_ys];
    let col_results: Vec<Result<Vec<u8>>> = col_blobs
        .into_par_iter()
        .map(bsc::decompress_parallel)
        .collect();
    let mut col_iter = col_results.into_iter();
    let d_rn = col_iter.next().unwrap()?;
    let d_ci = col_iter.next().unwrap()?;
    let d_ln = col_iter.next().unwrap()?;
    let d_ti = col_iter.next().unwrap()?;
    let d_xs = col_iter.next().unwrap()?;
    let d_ys = col_iter.next().unwrap()?;

    reconstruct_sra_headers(
        &name_prefix,
        &pair_suffix,
        &combos,
        &d_rn,
        &d_ci,
        &d_ln,
        &d_ti,
        &d_xs,
        &d_ys,
        wide,
        num_reads,
    )
}

/// Reconstruct SRA headers from decoded column buffers. Shared by the decompress
/// path (columns come from BSC decode) and the compress-time verifier (columns
/// are the in-memory pre-BSC bytes — same layout, so verification skips the
/// provably-lossless BSC round-trip). The byte layout of `d_*` is identical in
/// both cases: `d_rn` = u32-LE read numbers, `d_ci` = u8 combo indices, `d_ln` =
/// u8 lanes, `d_ti`/`d_xs`/`d_ys` = u16-LE (narrow) or u32-LE (wide) spatial.
#[allow(clippy::too_many_arguments)]
fn reconstruct_sra_headers(
    name_prefix: &str,
    pair_suffix: &str,
    combos: &[String],
    d_rn: &[u8],
    d_ci: &[u8],
    d_ln: &[u8],
    d_ti: &[u8],
    d_xs: &[u8],
    d_ys: &[u8],
    wide: bool,
    num_reads: usize,
) -> Result<Vec<String>> {
    let read_spatial = |blob: &[u8], idx: usize| -> Result<u32> {
        if wide {
            super::read_le_u32(blob, idx * 4)
        } else {
            super::read_le_u16(blob, idx * 2).map(|v| v as u32)
        }
    };

    (0..num_reads)
        .into_par_iter()
        .map(|i| {
            let rn = super::read_le_u32(d_rn, i * 4)?;
            let ci = *d_ci
                .get(i)
                .ok_or_else(|| anyhow::anyhow!("SRA: combo index out of bounds at read {i}"))?;
            let combo = combos.get(ci as usize).ok_or_else(|| {
                anyhow::anyhow!("SRA: combo dict entry {} out of range at read {i}", ci)
            })?;
            let lane = *d_ln
                .get(i)
                .ok_or_else(|| anyhow::anyhow!("SRA: lane out of bounds at read {i}"))?;
            let tile = read_spatial(d_ti, i)?;
            let x = read_spatial(d_xs, i)?;
            let y = read_spatial(d_ys, i)?;

            Ok(if pair_suffix.is_empty() {
                format!(
                    "@{}{} {}:{}:{}:{}:{}",
                    name_prefix, rn, combo, lane, tile, x, y
                )
            } else {
                format!(
                    "@{}{} {}:{}:{}:{}:{}{}",
                    name_prefix, rn, combo, lane, tile, x, y, pair_suffix
                )
            })
        })
        .collect()
}

/// Direct-to-varint variant of [`decompress_sra`]. Writes records into
/// per-chunk `Vec<u8>` buffers in parallel, then concatenates. No per-record
/// `String` allocation — see the rationale on
/// [`decompress_headers_columnar_to_varint_stream`].
fn decompress_sra_into_varint(data: &[u8], num_reads: usize, wide: bool) -> Result<Vec<u8>> {
    let mut pos = 0;

    let name_prefix = read_string(data, &mut pos)?;
    let pair_suffix = read_string(data, &mut pos)?;
    let combos = read_combo_dict(data, &mut pos)?;
    let stored_n = read_u32(data, &mut pos)? as usize;
    anyhow::ensure!(
        stored_n == num_reads,
        "SRA header count mismatch: stored {stored_n}, expected {num_reads}"
    );

    let s_rn = read_compressed_slice(data, &mut pos)?;
    let s_ci = read_compressed_slice(data, &mut pos)?;
    let s_ln = read_compressed_slice(data, &mut pos)?;
    let s_ti = read_compressed_slice(data, &mut pos)?;
    let s_xs = read_compressed_slice(data, &mut pos)?;
    let s_ys = read_compressed_slice(data, &mut pos)?;

    let col_blobs: [&[u8]; 6] = [s_rn, s_ci, s_ln, s_ti, s_xs, s_ys];
    let col_results: Vec<Result<Vec<u8>>> = col_blobs
        .into_par_iter()
        .map(bsc::decompress_parallel)
        .collect();
    let mut col_iter = col_results.into_iter();
    let d_rn = col_iter.next().unwrap()?;
    let d_ci = col_iter.next().unwrap()?;
    let d_ln = col_iter.next().unwrap()?;
    let d_ti = col_iter.next().unwrap()?;
    let d_xs = col_iter.next().unwrap()?;
    let d_ys = col_iter.next().unwrap()?;

    let read_spatial = |blob: &[u8], idx: usize| -> Result<u32> {
        if wide {
            super::read_le_u32(blob, idx * 4)
        } else {
            super::read_le_u16(blob, idx * 2).map(|v| v as u32)
        }
    };

    // Chunked-parallel: each worker fills its own buffer with varint-prefixed
    // records. CHUNK_RECORDS≈8K keeps per-chunk buffers near ~1 MB and amortizes
    // rayon dispatch.
    const CHUNK_RECORDS: usize = 8192;
    let num_chunks = num_reads.div_ceil(CHUNK_RECORDS);

    let chunks: Result<Vec<Vec<u8>>> = (0..num_chunks)
        .into_par_iter()
        .map(|c| -> Result<Vec<u8>> {
            let lo = c * CHUNK_RECORDS;
            let hi = ((c + 1) * CHUNK_RECORDS).min(num_reads);
            // ~100 bytes per header is a generous estimate that avoids realloc.
            let mut buf: Vec<u8> = Vec::with_capacity((hi - lo) * 100);
            for i in lo..hi {
                let rn = super::read_le_u32(&d_rn, i * 4)?;
                let ci = *d_ci
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("SRA: combo index out of bounds at read {i}"))?;
                let combo = combos.get(ci as usize).ok_or_else(|| {
                    anyhow::anyhow!("SRA: combo dict entry {} out of range at read {i}", ci)
                })?;
                let lane = *d_ln
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("SRA: lane out of bounds at read {i}"))?;
                let tile = read_spatial(&d_ti, i)?;
                let x = read_spatial(&d_xs, i)?;
                let y = read_spatial(&d_ys, i)?;

                let header_len = 1
                    + name_prefix.len()
                    + u32_decimal_len(rn)
                    + 1
                    + combo.len()
                    + 1
                    + u32_decimal_len(lane as u32)
                    + 1
                    + u32_decimal_len(tile)
                    + 1
                    + u32_decimal_len(x)
                    + 1
                    + u32_decimal_len(y)
                    + pair_suffix.len();

                write_varint(&mut buf, header_len);
                buf.push(b'@');
                buf.extend_from_slice(name_prefix.as_bytes());
                write!(buf, "{rn}").unwrap();
                buf.push(b' ');
                buf.extend_from_slice(combo.as_bytes());
                buf.push(b':');
                write!(buf, "{lane}").unwrap();
                buf.push(b':');
                write!(buf, "{tile}").unwrap();
                buf.push(b':');
                write!(buf, "{x}").unwrap();
                buf.push(b':');
                write!(buf, "{y}").unwrap();
                buf.extend_from_slice(pair_suffix.as_bytes());
            }
            Ok(buf)
        })
        .collect();

    let chunks = chunks?;
    let total: usize = chunks.iter().map(|b| b.len()).sum();
    let mut out = Vec::with_capacity(total);
    for buf in chunks {
        out.extend_from_slice(&buf);
    }
    Ok(out)
}

/// Decimal-digit count of a `u32` written via `Display` (matches `format!("{n}")`
/// length). Zero is one digit.
#[inline]
fn u32_decimal_len(n: u32) -> usize {
    if n == 0 { 1 } else { n.ilog10() as usize + 1 }
}

// ========================================================================
// Casava 1.8 format: @INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y COMMENT
// ========================================================================

fn parse_casava_row(header: &str) -> Option<CasavaRow<'_>> {
    // Require the FASTQ '@' prefix; the decoder re-adds a literal '@', so a
    // header without it (e.g. '>'-prefixed FASTA) must not be parsed here.
    let without_at = header.strip_prefix('@')?;
    let (name_part, comment) = without_at.split_once(' ')?;
    let fields: Vec<&str> = name_part.split(':').collect();
    if fields.len() < 7 {
        return None;
    }
    // Numeric fields must render canonically (no leading zeros / '+') or a
    // parse-then-format round trip would corrupt them on decompress.
    if !is_canonical_decimal(fields[3])
        || !is_canonical_decimal(fields[4])
        || !is_canonical_decimal(fields[5])
        || !is_canonical_decimal(fields[6])
    {
        return None;
    }
    // Zero-copy combo slice: fields[0]:fields[1]:fields[2] (no allocation)
    let combo_len = fields[0].len() + 1 + fields[1].len() + 1 + fields[2].len();
    let combo = &name_part[..combo_len];
    let lane: u8 = fields[3].parse().ok()?;
    let tile: u32 = fields[4].parse().ok()?;
    let x: u32 = fields[5].parse().ok()?;
    let y: u32 = fields[6].parse().ok()?;
    // UMI: everything after field[6] (may be multi-part with ':')
    let umi = if fields.len() >= 8 {
        let no_umi_len = combo_len
            + 1
            + fields[3].len()
            + 1
            + fields[4].len()
            + 1
            + fields[5].len()
            + 1
            + fields[6].len();
        Some(&name_part[no_umi_len + 1..])
    } else {
        None
    };
    Some(CasavaRow {
        combo,
        lane,
        tile,
        x,
        y,
        umi,
        comment,
    })
}

fn compress_casava(headers: &[&str]) -> Result<Option<Vec<u8>>> {
    let n = headers.len();

    // Phase 1 (parallel): parse all headers into Row structs.
    // combo, umi, and comment are zero-copy &str slices into the original headers.
    let rows: Vec<Option<CasavaRow<'_>>> =
        headers.par_iter().map(|h| parse_casava_row(h)).collect();

    // Phase 2 (sequential): build combo dict, detect common comment, fill columns.
    let first_comment = rows[0].as_ref().map_or("", |r| r.comment);
    let mut all_same_comment = true;
    let mut has_umi = false;
    let mut combos: Vec<&str> = Vec::new();
    let mut combo_map: HashMap<&str, u8> = HashMap::new();
    let mut combo_idxs = Vec::with_capacity(n);
    let mut lanes = Vec::with_capacity(n);
    // Collect spatial as u32; pick wire width after the full scan.
    let mut tiles_u32 = Vec::with_capacity(n);
    let mut xs_u32 = Vec::with_capacity(n);
    let mut ys_u32 = Vec::with_capacity(n);

    for row_opt in &rows {
        let row = match row_opt {
            Some(r) => r,
            None => return Ok(None),
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
            // >= 255 distinct combos won't fit the dict's u8 index, and a single
            // combo > 255 bytes won't fit its u8 length prefix; either way fall
            // back to raw rather than panicking in write_combo_dict.
            if combos.len() >= 255 || row.combo.len() > u8::MAX as usize {
                return Ok(None);
            }
            let i = combos.len() as u8;
            combos.push(row.combo);
            combo_map.insert(row.combo, i);
            i
        };
        combo_idxs.push(idx);
        lanes.push(row.lane);
        tiles_u32.push(row.tile);
        xs_u32.push(row.x);
        ys_u32.push(row.y);
    }

    // Pick the narrowest format that fits — v0x02 (u16) for typical Illumina,
    // v0x03 (u32) only when NovaSeq X-scale x/y values overflow.
    let narrow = fits_in_u16(&tiles_u32) && fits_in_u16(&xs_u32) && fits_in_u16(&ys_u32);
    let mut tiles = Vec::with_capacity(n * if narrow { 2 } else { 4 });
    let mut xs = Vec::with_capacity(n * if narrow { 2 } else { 4 });
    let mut ys = Vec::with_capacity(n * if narrow { 2 } else { 4 });
    if narrow {
        for v in &tiles_u32 {
            tiles.extend_from_slice(&(*v as u16).to_le_bytes());
        }
        for v in &xs_u32 {
            xs.extend_from_slice(&(*v as u16).to_le_bytes());
        }
        for v in &ys_u32 {
            ys.extend_from_slice(&(*v as u16).to_le_bytes());
        }
    } else {
        for v in &tiles_u32 {
            tiles.extend_from_slice(&v.to_le_bytes());
        }
        for v in &xs_u32 {
            xs.extend_from_slice(&v.to_le_bytes());
        }
        for v in &ys_u32 {
            ys.extend_from_slice(&v.to_le_bytes());
        }
    }
    let version = if narrow { 0x02u8 } else { 0x03u8 };

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

    // Verify the structured encoding round-trips byte-exact BEFORE paying for
    // BSC, reconstructing from the in-memory pre-BSC columns (same byte layout
    // as the BSC-decoded columns). Skips the redundant, provably-lossless BSC
    // round-trip; see the matching rationale in compress_sra. Falls back to the
    // byte-exact raw encoding on any mismatch (or reconstruct error).
    // u16 write_string guard (see compress_sra): a common comment longer than 64 KiB
    // would truncate its stored length, producing an undecodable blob. Fall back to raw.
    if all_same_comment && first_comment.len() > u16::MAX as usize {
        return Ok(None);
    }

    let combos_owned: Vec<String> = combos.iter().map(|s| s.to_string()).collect();
    let common_comment_ref = if all_same_comment {
        Some(first_comment)
    } else {
        None
    };
    let comment_data_ref = if all_same_comment {
        None
    } else {
        Some(comments.as_slice())
    };
    let umi_data_ref = if has_umi { Some(umis.as_slice()) } else { None };
    let verified = reconstruct_casava_headers(
        &combos_owned,
        common_comment_ref,
        has_umi,
        &combo_idxs,
        &lanes,
        &tiles,
        &xs,
        &ys,
        comment_data_ref,
        umi_data_ref,
        !narrow,
        n,
    )
    .map(|recon| {
        recon.len() == n
            && recon
                .iter()
                .zip(headers)
                .all(|(r, h)| r.as_bytes() == h.as_bytes())
    })
    .unwrap_or(false);
    if !verified {
        return Ok(None);
    }

    // BSC-compress all 5 spatial columns in parallel
    let (c_ci, (c_ln, (c_ti, (c_xs, c_ys)))) = rayon::join(
        || bsc::compress_parallel_adaptive(&combo_idxs),
        || {
            rayon::join(
                || bsc::compress_parallel_adaptive(&lanes),
                || {
                    rayon::join(
                        || bsc::compress_parallel_adaptive(&tiles),
                        || {
                            rayon::join(
                                || bsc::compress_parallel_adaptive(&xs),
                                || bsc::compress_parallel_adaptive(&ys),
                            )
                        },
                    )
                },
            )
        },
    );
    let c_ci = c_ci?;
    let c_ln = c_ln?;
    let c_ti = c_ti?;
    let c_xs = c_xs?;
    let c_ys = c_ys?;

    // Build output
    let mut out = Vec::new();
    out.push(version); // 0x02 (u16 tile/x/y) or 0x03 (u32 tile/x/y)
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

    Ok(Some(out))
}

fn decompress_casava(data: &[u8], num_reads: usize, wide: bool) -> Result<Vec<String>> {
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
    anyhow::ensure!(
        stored_n == num_reads,
        "Casava header count mismatch: stored {stored_n}, expected {num_reads}"
    );

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

    reconstruct_casava_headers(
        &combos,
        common_comment.as_deref(),
        has_umi,
        &d_ci,
        &d_ln,
        &d_ti,
        &d_xs,
        &d_ys,
        comment_data.as_deref(),
        umi_data.as_deref(),
        wide,
        num_reads,
    )
}

/// Reconstruct Casava headers from decoded column buffers. Shared by the
/// decompress path (columns from BSC decode) and the compress-time verifier
/// (columns are the in-memory pre-BSC bytes — same layout, so verification skips
/// the provably-lossless BSC round-trip). `comment_data`/`umi_data` are the
/// varint-length-prefixed column blobs (`[varint(len)][bytes]`*); `common_comment`
/// is `Some` iff every read shares one comment (then `comment_data` is `None`).
#[allow(clippy::too_many_arguments)]
fn reconstruct_casava_headers(
    combos: &[String],
    common_comment: Option<&str>,
    has_umi: bool,
    d_ci: &[u8],
    d_ln: &[u8],
    d_ti: &[u8],
    d_xs: &[u8],
    d_ys: &[u8],
    comment_data: Option<&[u8]>,
    umi_data: Option<&[u8]>,
    wide: bool,
    num_reads: usize,
) -> Result<Vec<String>> {
    // Pre-scan variable-length comment offsets (sequential, needed before parallel phase)
    let comment_offsets: Option<Vec<(usize, usize)>> = if let Some(cd) = comment_data {
        let mut offsets = Vec::with_capacity(num_reads);
        let mut off = 0usize;
        for i in 0..num_reads {
            let len = read_varint_from_slice(cd, &mut off)?;
            let end = off
                .checked_add(len)
                .filter(|&e| e <= cd.len())
                .ok_or_else(|| {
                    anyhow::anyhow!("header_col: truncated comment data for header {i}")
                })?;
            offsets.push((off, end));
            off = end;
        }
        Some(offsets)
    } else {
        None
    };

    // Pre-scan variable-length UMI offsets
    let umi_offsets: Option<Vec<(usize, usize)>> = if let Some(ud) = umi_data {
        let mut offsets = Vec::with_capacity(num_reads);
        let mut off = 0usize;
        for i in 0..num_reads {
            let len = read_varint_from_slice(ud, &mut off)?;
            let end = off
                .checked_add(len)
                .filter(|&e| e <= ud.len())
                .ok_or_else(|| anyhow::anyhow!("header_col: truncated UMI data for header {i}"))?;
            offsets.push((off, end));
            off = end;
        }
        Some(offsets)
    } else {
        None
    };

    let read_spatial = |blob: &[u8], idx: usize| -> Result<u32> {
        if wide {
            super::read_le_u32(blob, idx * 4)
        } else {
            super::read_le_u16(blob, idx * 2).map(|v| v as u32)
        }
    };

    (0..num_reads)
        .into_par_iter()
        .map(|i| {
            let ci = *d_ci
                .get(i)
                .ok_or_else(|| anyhow::anyhow!("Casava: combo index out of bounds at read {i}"))?;
            let combo = combos.get(ci as usize).ok_or_else(|| {
                anyhow::anyhow!("Casava: combo dict entry {} out of range at read {i}", ci)
            })?;
            let lane = *d_ln
                .get(i)
                .ok_or_else(|| anyhow::anyhow!("Casava: lane out of bounds at read {i}"))?;
            let tile = read_spatial(d_ti, i)?;
            let x = read_spatial(d_xs, i)?;
            let y = read_spatial(d_ys, i)?;

            let comment = match common_comment {
                Some(cc) => cc,
                None => {
                    let cd = comment_data.unwrap();
                    let (s, e) = comment_offsets.as_ref().unwrap()[i];
                    std::str::from_utf8(&cd[s..e]).context("invalid comment UTF-8")?
                }
            };

            Ok(if has_umi {
                let ud = umi_data.unwrap();
                let (s, e) = umi_offsets.as_ref().unwrap()[i];
                let umi = std::str::from_utf8(&ud[s..e]).context("invalid UMI UTF-8")?;
                format!(
                    "@{}:{}:{}:{}:{}:{} {}",
                    combo, lane, tile, x, y, umi, comment
                )
            } else {
                format!("@{}:{}:{}:{}:{} {}", combo, lane, tile, x, y, comment)
            })
        })
        .collect()
}

/// Direct-to-varint variant of [`decompress_casava`]. Same chunked-parallel
/// design as [`decompress_sra_into_varint`] — see that function's docstring
/// for the rationale.
fn decompress_casava_into_varint(data: &[u8], num_reads: usize, wide: bool) -> Result<Vec<u8>> {
    let mut pos = 0;

    let combos = read_combo_dict(data, &mut pos)?;

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

    if pos >= data.len() {
        anyhow::bail!("Casava header: truncated at UMI flag");
    }
    let has_umi = data[pos] == 0x01;
    pos += 1;

    let stored_n = read_u32(data, &mut pos)? as usize;
    anyhow::ensure!(
        stored_n == num_reads,
        "Casava header count mismatch: stored {stored_n}, expected {num_reads}"
    );

    let d_ci = read_and_decompress(data, &mut pos)?;
    let d_ln = read_and_decompress(data, &mut pos)?;
    let d_ti = read_and_decompress(data, &mut pos)?;
    let d_xs = read_and_decompress(data, &mut pos)?;
    let d_ys = read_and_decompress(data, &mut pos)?;

    let comment_data = if common_comment.is_none() {
        Some(read_and_decompress(data, &mut pos)?)
    } else {
        None
    };

    let umi_data = if has_umi {
        Some(read_and_decompress(data, &mut pos)?)
    } else {
        None
    };

    // Serial pre-scan of variable-length comment/UMI offsets — same as the
    // string-returning path. Cheap (just integer walks).
    let comment_offsets: Option<Vec<(usize, usize)>> = if let Some(cd) = &comment_data {
        let mut offsets = Vec::with_capacity(num_reads);
        let mut off = 0usize;
        for i in 0..num_reads {
            let len = read_varint_from_slice(cd, &mut off)?;
            let end = off
                .checked_add(len)
                .filter(|&e| e <= cd.len())
                .ok_or_else(|| {
                    anyhow::anyhow!("header_col: truncated comment data for header {i}")
                })?;
            offsets.push((off, end));
            off = end;
        }
        Some(offsets)
    } else {
        None
    };

    let umi_offsets: Option<Vec<(usize, usize)>> = if let Some(ud) = &umi_data {
        let mut offsets = Vec::with_capacity(num_reads);
        let mut off = 0usize;
        for i in 0..num_reads {
            let len = read_varint_from_slice(ud, &mut off)?;
            let end = off
                .checked_add(len)
                .filter(|&e| e <= ud.len())
                .ok_or_else(|| anyhow::anyhow!("header_col: truncated UMI data for header {i}"))?;
            offsets.push((off, end));
            off = end;
        }
        Some(offsets)
    } else {
        None
    };

    let read_spatial = |blob: &[u8], idx: usize| -> Result<u32> {
        if wide {
            super::read_le_u32(blob, idx * 4)
        } else {
            super::read_le_u16(blob, idx * 2).map(|v| v as u32)
        }
    };

    const CHUNK_RECORDS: usize = 8192;
    let num_chunks = num_reads.div_ceil(CHUNK_RECORDS);

    let chunks: Result<Vec<Vec<u8>>> = (0..num_chunks)
        .into_par_iter()
        .map(|c| -> Result<Vec<u8>> {
            let lo = c * CHUNK_RECORDS;
            let hi = ((c + 1) * CHUNK_RECORDS).min(num_reads);
            let mut buf: Vec<u8> = Vec::with_capacity((hi - lo) * 100);
            for i in lo..hi {
                let ci = *d_ci.get(i).ok_or_else(|| {
                    anyhow::anyhow!("Casava: combo index out of bounds at read {i}")
                })?;
                let combo = combos.get(ci as usize).ok_or_else(|| {
                    anyhow::anyhow!("Casava: combo dict entry {} out of range at read {i}", ci)
                })?;
                let lane = *d_ln
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("Casava: lane out of bounds at read {i}"))?;
                let tile = read_spatial(&d_ti, i)?;
                let x = read_spatial(&d_xs, i)?;
                let y = read_spatial(&d_ys, i)?;

                let comment_bytes: &[u8] = match &common_comment {
                    Some(cc) => cc.as_bytes(),
                    None => {
                        let cd = comment_data.as_ref().unwrap();
                        let (s, e) = comment_offsets.as_ref().unwrap()[i];
                        // UTF-8 validation done lazily — we wrote these bytes
                        // unchanged, so they are valid UTF-8 by construction.
                        // Re-validate to keep the same invariant as the slow path.
                        let slice = &cd[s..e];
                        std::str::from_utf8(slice).context("invalid comment UTF-8")?;
                        slice
                    }
                };

                let umi_bytes: Option<&[u8]> = if has_umi {
                    let ud = umi_data.as_ref().unwrap();
                    let (s, e) = umi_offsets.as_ref().unwrap()[i];
                    let slice = &ud[s..e];
                    std::str::from_utf8(slice).context("invalid UMI UTF-8")?;
                    Some(slice)
                } else {
                    None
                };

                // '@' + combo + ':' + d(lane) + ':' + d(tile) + ':' + d(x) + ':' + d(y)
                //   + (':' + umi  if has_umi)
                //   + ' ' + comment
                let mut header_len = 1
                    + combo.len()
                    + 1
                    + u32_decimal_len(lane as u32)
                    + 1
                    + u32_decimal_len(tile)
                    + 1
                    + u32_decimal_len(x)
                    + 1
                    + u32_decimal_len(y);
                if let Some(umi) = umi_bytes {
                    header_len += 1 + umi.len();
                }
                header_len += 1 + comment_bytes.len();

                write_varint(&mut buf, header_len);
                buf.push(b'@');
                buf.extend_from_slice(combo.as_bytes());
                buf.push(b':');
                write!(buf, "{lane}").unwrap();
                buf.push(b':');
                write!(buf, "{tile}").unwrap();
                buf.push(b':');
                write!(buf, "{x}").unwrap();
                buf.push(b':');
                write!(buf, "{y}").unwrap();
                if let Some(umi) = umi_bytes {
                    buf.push(b':');
                    buf.extend_from_slice(umi);
                }
                buf.push(b' ');
                buf.extend_from_slice(comment_bytes);
            }
            Ok(buf)
        })
        .collect();

    let chunks = chunks?;
    let total: usize = chunks.iter().map(|b| b.len()).sum();
    let mut out = Vec::with_capacity(total);
    for buf in chunks {
        out.extend_from_slice(&buf);
    }
    Ok(out)
}

// ========================================================================
// Shared helpers
// ========================================================================

fn write_string(out: &mut Vec<u8>, s: &str) {
    // The length is stored as a u16; callers (compress_sra / compress_casava) MUST guard
    // and fall back to raw for strings longer than u16::MAX, or the prefix truncates while
    // the full body is written, yielding an undecodable blob.
    debug_assert!(
        s.len() <= u16::MAX as usize,
        "write_string: {} bytes exceeds the u16 length prefix; caller must guard",
        s.len()
    );
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
        // Individual combo strings fit in one byte length prefix: the combo
        // insertion sites fall back to raw encoding for any combo > 255 bytes,
        // so this conversion cannot fail for combos that reach here.
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
    let v = u32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]);
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
    let mut shift = 0u32;
    loop {
        if *offset >= data.len() {
            anyhow::bail!("header_col: truncated varint at offset {}", *offset);
        }
        // Reject an over-long varint before the shift can reach/exceed the width
        // of `usize` (a crafted run of 0x80 bytes would otherwise be a debug-mode
        // "shift left with overflow" panic / a silently-masked shift in release).
        if shift >= usize::BITS {
            anyhow::bail!("header_col: varint too long at offset {}", *offset);
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

pub(super) fn compress_raw_fallback(headers: &[&str]) -> Result<Vec<u8>> {
    let bytes: Vec<&[u8]> = headers.iter().map(|h| h.as_bytes()).collect();
    compress_raw_fallback_bytes(&bytes)
}

fn compress_raw_fallback_bytes(headers: &[&[u8]]) -> Result<Vec<u8>> {
    let mut stream = Vec::new();
    for h in headers {
        write_varint(&mut stream, h.len());
        stream.extend_from_slice(h);
    }
    let compressed = bsc::compress_parallel_adaptive(&stream)?;

    let mut out = Vec::new();
    out.push(0x00);
    out.write_all(&compressed)?;
    Ok(out)
}

/// Decode the raw fallback to byte-exact headers. Headers are stored as raw
/// bytes (no UTF-8 validation) so non-UTF-8 IDs survive the round-trip.
fn decompress_raw_fallback(data: &[u8], num_reads: usize) -> Result<Vec<Vec<u8>>> {
    // Bound the in-memory fallback decode against a stacked BSC bomb (tiny payload →
    // huge output): a valid header stream holds `num_reads` varint-prefixed headers, so
    // `stream_decode_cap(num_reads)` is a safe upper bound (it floors at one 64 MiB block).
    let decompressed = bsc::decompress_parallel_max(
        data,
        super::codecs::stream_decode_cap(num_reads),
    )?;
    let mut headers = Vec::with_capacity(super::decompress_impl::scan_capacity(
        num_reads,
        decompressed.len(),
    ));
    let mut offset = 0;

    for i in 0..num_reads {
        let len = read_varint_from_slice(&decompressed, &mut offset)?;
        let end = offset
            .checked_add(len)
            .filter(|&e| e <= decompressed.len())
            .ok_or_else(|| anyhow::anyhow!("header_col: truncated header data for read {i}"))?;
        headers.push(decompressed[offset..end].to_vec());
        offset = end;
    }

    Ok(headers)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: reference varint serialization of `Vec<Vec<u8>>`, used to compare
    /// the fast-path output of [`decompress_headers_columnar_to_varint_stream`].
    fn reference_varint_stream(headers: &[Vec<u8>]) -> Vec<u8> {
        let mut out = Vec::new();
        for h in headers {
            write_varint(&mut out, h.len());
            out.extend_from_slice(h);
        }
        out
    }

    /// Compress then decompress; assert every header round-trips byte-exact.
    fn assert_headers_roundtrip(headers: &[&str]) {
        let compressed = compress_headers_columnar(headers).unwrap();
        let out = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        assert_eq!(out.len(), headers.len(), "record count mismatch");
        for (h, o) in headers.iter().zip(out.iter()) {
            assert_eq!(
                h.as_bytes(),
                o.as_slice(),
                "header mismatch: {:?} != {:?}",
                h,
                String::from_utf8_lossy(o)
            );
        }
    }

    #[test]
    fn sra_huge_name_prefix_roundtrips_via_fallback() {
        // `write_string` length-prefixes with a u16. A shared name prefix longer than
        // 64 KiB would truncate that prefix while writing the full bytes, producing a
        // columnar blob the decoder mis-reads — a silently-undecodable archive. The
        // encoder must instead fall back to the raw path and round-trip byte-exact.
        let big = "P".repeat(70_000);
        let headers: Vec<String> = (1..=300u32)
            .map(|i| format!("@{big}.{i} A00296:43:HCLHLDSXX:4:2174:25834:28823"))
            .collect();
        let refs: Vec<&str> = headers.iter().map(|s| s.as_str()).collect();
        assert_headers_roundtrip(&refs);
    }

    #[test]
    fn test_sra_interleaved_pair_suffix_roundtrips() {
        // Interleaved /1 + /2 mates: the encoder must not assume read 0's pair
        // suffix applies to every read.
        let headers = vec![
            "@ERR1.1 A00296:43:HCLHLDSXX:4:2174:25834:28823/1",
            "@ERR1.2 A00296:43:HCLHLDSXX:4:2174:25834:28823/2",
            "@ERR1.3 A00296:43:HCLHLDSXX:4:2174:25835:28824/1",
            "@ERR1.4 A00296:43:HCLHLDSXX:4:2174:25835:28824/2",
        ];
        assert_headers_roundtrip(&headers);
    }

    #[test]
    fn test_sra_multi_accession_prefix_roundtrips() {
        // Merged accessions: the encoder must not assume read 0's name prefix
        // applies to every read.
        let headers = vec![
            "@ERR1.1 A00296:43:HCLHLDSXX:4:2174:25834:28823",
            "@ERR2.1 A00296:43:HCLHLDSXX:4:2174:25834:28823",
        ];
        assert_headers_roundtrip(&headers);
    }

    #[test]
    fn test_sra_leading_zero_fields_roundtrip() {
        // Zero-padded numeric fields must not be normalized away.
        let headers = vec![
            "@ERR1.007 A00296:43:HCLHLDSXX:04:0073:0941:1973",
            "@ERR1.008 A00296:43:HCLHLDSXX:04:0073:0941:1974",
        ];
        assert_headers_roundtrip(&headers);
    }

    #[test]
    fn test_casava_leading_zero_fields_roundtrip() {
        let headers = vec![
            "@A00296:43:HCLHLDSXX:04:0073:0941:1973 1:N:0:ATCG",
            "@A00296:43:HCLHLDSXX:04:0073:0941:1974 1:N:0:ATCG",
        ];
        assert_headers_roundtrip(&headers);
    }

    /// The fast `_to_varint_stream` path must produce the same bytes as the
    /// slow `Vec<Vec<u8>>` path followed by manual varint reframing — across
    /// SRA (with and without pair), Casava (with common comment, per-record
    /// comment, with UMI, narrow and wide spatial), and the raw fallback.
    #[test]
    fn test_varint_stream_matches_slow_path_sra() {
        // All reads share the same /1 pair suffix (a single-end or R1-only set);
        // a mix of /1 and /2 correctly falls back to raw and is covered by
        // test_sra_interleaved_pair_suffix_roundtrips.
        let headers = vec![
            "@ERR3239334.100 A00296:43:HCLHLDSXX:4:2174:25834:28823/1",
            "@ERR3239334.200 A00296:43:HCLHLDSXX:4:2174:25835:28824/1",
            "@ERR3239334.300 A00217:77:HFJWFDSXX:2:1101:10000:5000/1",
        ];
        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x01);

        let slow = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        let fast =
            decompress_headers_columnar_to_varint_stream(&compressed, headers.len()).unwrap();
        assert_eq!(reference_varint_stream(&slow), fast);
    }

    #[test]
    fn test_varint_stream_matches_slow_path_sra_no_pair() {
        let headers = vec![
            "@ERR3239334.100 A00296:43:HCLHLDSXX:4:2174:25834:28823",
            "@ERR3239334.200 A00296:43:HCLHLDSXX:4:2174:25835:28824",
        ];
        let compressed = compress_headers_columnar(&headers).unwrap();
        let slow = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        let fast =
            decompress_headers_columnar_to_varint_stream(&compressed, headers.len()).unwrap();
        assert_eq!(reference_varint_stream(&slow), fast);
    }

    #[test]
    fn test_varint_stream_matches_slow_path_casava() {
        let headers = vec![
            "@A00296:43:HCLHLDSXX:4:2174:25834:28823 1:N:0:ATCGATCG",
            "@A00296:43:HCLHLDSXX:4:2174:25835:28824 1:N:0:ATCGATCG",
            "@A00296:43:HCLHLDSXX:4:2174:25836:28825 1:N:0:ATCGATCG",
        ];
        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x02);
        let slow = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        let fast =
            decompress_headers_columnar_to_varint_stream(&compressed, headers.len()).unwrap();
        assert_eq!(reference_varint_stream(&slow), fast);
    }

    #[test]
    fn test_varint_stream_matches_slow_path_casava_varying_comments() {
        let headers = vec![
            "@A00296:43:HCLHLDSXX:4:2174:25834:28823 1:N:0:ATCG+CGAT",
            "@A00296:43:HCLHLDSXX:4:2174:25835:28824 1:N:0:GCTA+TAGC",
            "@A00296:43:HCLHLDSXX:4:2174:25836:28825 2:Y:0:ATCG+CGAT",
        ];
        let compressed = compress_headers_columnar(&headers).unwrap();
        let slow = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        let fast =
            decompress_headers_columnar_to_varint_stream(&compressed, headers.len()).unwrap();
        assert_eq!(reference_varint_stream(&slow), fast);
    }

    #[test]
    fn test_varint_stream_matches_slow_path_casava_with_umi() {
        let headers = vec![
            "@A00296:43:HCLHLDSXX:4:2174:25834:28823:ACGTACGT 1:N:0:ATCG",
            "@A00296:43:HCLHLDSXX:4:2174:25835:28824:TTTTTTTT 1:N:0:ATCG",
            "@A00296:43:HCLHLDSXX:4:2174:25836:28825:GGGGGGGG 1:N:0:ATCG",
        ];
        let compressed = compress_headers_columnar(&headers).unwrap();
        let slow = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        let fast =
            decompress_headers_columnar_to_varint_stream(&compressed, headers.len()).unwrap();
        assert_eq!(reference_varint_stream(&slow), fast);
    }

    #[test]
    fn test_varint_stream_matches_slow_path_raw_fallback() {
        // Headers that don't fit SRA or Casava format — exercises the 0x00
        // path which falls through to the slow re-serialization.
        let headers = vec!["@unstructured/foo", "@bar baz", "@whatever"];
        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x00);
        let slow = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        let fast =
            decompress_headers_columnar_to_varint_stream(&compressed, headers.len()).unwrap();
        assert_eq!(reference_varint_stream(&slow), fast);
    }

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
    fn test_roundtrip_casava_novaseq_x_wide() {
        // NovaSeq X x/y can reach ~70k, overflowing u16. The encoder must pick
        // v0x03 (u32 spatial fields) and the decoder must read both widths.
        let headers = vec![
            "@LH00001:1:22F2WLLT4:1:1101:70000:65537 1:N:0:ATCACG",
            "@LH00001:1:22F2WLLT4:1:1101:71234:80000 1:N:0:ATCACG",
            "@LH00001:1:22F2WLLT4:2:1102:80001:90000 1:N:0:ATCACG",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x03, "wide-spatial Casava must use v0x03");
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_roundtrip_casava_mixed_widths_picks_wide() {
        // Most rows fit u16, but ONE row overflows y. The encoder must pick
        // v0x03 for the whole batch (the contract: per-archive width, not
        // per-row). A future refactor that did per-row width selection would
        // break the format silently — this test guards against it.
        let headers = vec![
            "@LH00001:1:22F2WLLT4:1:1101:1000:2000 1:N:0:ATCACG",
            "@LH00001:1:22F2WLLT4:1:1101:1100:2100 1:N:0:ATCACG",
            "@LH00001:1:22F2WLLT4:1:1101:1200:80000 1:N:0:ATCACG", // y overflows u16
            "@LH00001:1:22F2WLLT4:1:1101:1300:2300 1:N:0:ATCACG",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(
            compressed[0], 0x03,
            "any row overflowing u16 must promote the whole archive to v0x03"
        );
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_roundtrip_sra_wide_spatial() {
        // Same wide-spatial test for SRA-style headers — must pick v0x04.
        let headers = vec![
            "@ERR3239334.100 A00296:43:HCLHLDSXX:4:2174:70000:65537/1",
            "@ERR3239334.200 A00296:43:HCLHLDSXX:4:2174:71234:80000/1",
        ];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x04, "wide-spatial SRA must use v0x04");
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_fallback_non_illumina() {
        let headers = vec!["@READ_1 some comment", "@READ_2 another comment"];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x00); // raw fallback
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(orig.as_bytes(), dec.as_slice());
        }
    }

    #[test]
    fn test_valid_casava_stays_structured() {
        // Genuine Casava headers must still use the compact structured encoding —
        // the in-memory verify backstop must not cause a needless raw fallback.
        let headers = vec![
            "@INSTR:43:HCLHLDSXX:4:2174:25834:28823 1:N:0:ATCG",
            "@INSTR:43:HCLHLDSXX:4:2174:25835:28824 1:N:0:ATCG",
        ];
        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_ne!(
            compressed[0], 0x00,
            "valid Casava must stay structured, not fall back"
        );
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (o, d) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(o.as_bytes(), d.as_slice());
        }
    }

    #[test]
    fn test_casava_unreconstructable_falls_back() {
        // A Casava header that PARSES (numeric fields, 7 colon parts) but cannot
        // be reconstructed byte-exact — here a zero-padded tile/x/y, stored as
        // integers and re-emitted without the padding — must be caught by the
        // in-memory verify in `compress_casava` and degrade to the byte-exact
        // generic (0x05) or raw (0x00) fallback, still round-tripping losslessly.
        // This exercises the verify path specifically (parse succeeds, reconstruction
        // differs), the replacement for the old `verified_or_raw` round-trip check.
        let headers = vec![
            "@INSTR:43:HCLHLDSXX:4:0073:0941:1973 1:N:0:ATCG",
            "@INSTR:43:HCLHLDSXX:4:0073:0941:1974 1:N:0:ATCG",
        ];
        let compressed = compress_headers_columnar(&headers).unwrap();
        assert!(
            compressed[0] != 0x02 && compressed[0] != 0x03,
            "unreconstructable Casava must not use the Casava encoding (got 0x{:02x})",
            compressed[0]
        );
        let decoded = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (o, d) in headers.iter().zip(decoded.iter()) {
            assert_eq!(o.as_bytes(), d.as_slice());
        }
    }

    #[test]
    fn test_sra_unreconstructable_falls_back() {
        // SRA counterpart: zero-padded spatial fields parse but don't reconstruct
        // byte-exact, so compress_sra's in-memory verify must fall back to raw.
        let headers = vec![
            "@ERR1.1 A00296:43:HCLHLDSXX:4:0073:0941:1973",
            "@ERR1.2 A00296:43:HCLHLDSXX:4:0073:0941:1974",
        ];
        let compressed = compress_headers_columnar(&headers).unwrap();
        // Must NOT use the SRA structured encodings (0x01/0x04); falls back to
        // generic (0x05) or raw (0x00), both lossless.
        assert!(
            compressed[0] != 0x01 && compressed[0] != 0x04,
            "unreconstructable SRA must not use the SRA encoding (got 0x{:02x})",
            compressed[0]
        );
        let decoded = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (o, d) in headers.iter().zip(decoded.iter()) {
            assert_eq!(o.as_bytes(), d.as_slice());
        }
    }

    #[test]
    fn test_fasta_casava_shaped_header_uses_raw_fallback() {
        // A '>'-prefixed FASTA header that happens to be Casava-shaped (7 colon
        // fields, numeric lane/tile/x/y) must NOT be parsed as structured Casava:
        // the structured decoders re-add a literal '@', corrupting '>id' into
        // '@>id'. Structured encoding only applies to FASTQ ('@') headers, so
        // these must take the byte-exact raw fallback.
        let headers = vec![
            ">INSTR:RUN:FLOWCELL:1:1000:2000:3000 desc",
            ">INSTR:RUN:FLOWCELL:1:1001:2001:3001 desc",
        ];
        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(
            compressed[0], 0x00,
            "'>'-prefixed headers must use raw fallback, not structured Casava"
        );
        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        for (orig, dec) in headers.iter().zip(decompressed.iter()) {
            assert_eq!(
                orig.as_bytes(),
                dec.as_slice(),
                "header must round-trip byte-exact"
            );
        }
    }

    #[test]
    fn test_empty() {
        let headers: Vec<&str> = vec![];
        let compressed = compress_headers_columnar(&headers).unwrap();
        let decompressed = decompress_headers_columnar(&compressed, 0).unwrap();
        assert!(decompressed.is_empty());
    }

    #[test]
    fn test_oversized_combo_falls_back_to_raw_without_panic() {
        // A casava-shaped header whose instrument field is > 255 bytes parses
        // fine (lane/tile/x/y are numeric) but yields a combo string too long
        // for the one-byte combo-dict length prefix. This must fall back to
        // raw encoding, not panic in write_combo_dict (which under
        // panic="abort" would kill the process).
        let long_instrument = "A".repeat(300);
        let h1 = format!("@{long_instrument}:43:HCLHLDSXX:4:2174:25834:28823 1:N:0:ATCG");
        let h2 = format!("@{long_instrument}:43:HCLHLDSXX:4:2174:25835:28824 1:N:0:ATCG");
        let headers = vec![h1.as_str(), h2.as_str()];

        let compressed = compress_headers_columnar(&headers).unwrap();
        assert_eq!(compressed[0], 0x00, "oversized combo must use raw fallback");

        let decompressed = decompress_headers_columnar(&compressed, headers.len()).unwrap();
        assert_eq!(decompressed[0].as_slice(), h1.as_bytes());
        assert_eq!(decompressed[1].as_slice(), h2.as_bytes());
    }

    #[test]
    fn test_columnar_bytes_preserves_non_utf8_headers() {
        // A header containing a non-UTF-8 byte is never a valid Illumina header;
        // it must round-trip byte-exact via the raw fallback, not be mangled by
        // a lossy UTF-8 conversion.
        let mut h1 = b"@read".to_vec();
        h1.push(0xFF);
        h1.extend_from_slice(b"1");
        let h2 = b"@plain_header_2".to_vec();
        let headers: Vec<&[u8]> = vec![h1.as_slice(), h2.as_slice()];

        let compressed = compress_headers_columnar_bytes(&headers).unwrap();
        assert_eq!(
            compressed[0], 0x00,
            "non-UTF-8 headers must use raw fallback"
        );

        let decompressed = decompress_headers_columnar(&compressed, 2).unwrap();
        assert_eq!(decompressed[0], h1);
        assert_eq!(decompressed[1], h2);
    }

    #[test]
    fn test_columnar_bytes_still_uses_structured_for_illumina() {
        // Valid Illumina (ASCII) headers still take the structured path.
        let h1 = b"@A00296:43:HCLHLDSXX:4:2174:25834:28823 1:N:0:ATCG".to_vec();
        let h2 = b"@A00296:43:HCLHLDSXX:4:2174:25835:28824 1:N:0:ATCG".to_vec();
        let headers: Vec<&[u8]> = vec![h1.as_slice(), h2.as_slice()];

        let compressed = compress_headers_columnar_bytes(&headers).unwrap();
        assert_eq!(
            compressed[0], 0x02,
            "Illumina headers should use Casava encoding"
        );

        let decompressed = decompress_headers_columnar(&compressed, 2).unwrap();
        assert_eq!(decompressed[0], h1);
        assert_eq!(decompressed[1], h2);
    }

    #[test]
    fn generic_tier_handles_non_illumina_and_roundtrips() {
        let owned: Vec<String> = (1..=2000)
            .map(|i| format!("@m54238_180628/{}/0_{}", 4194300 + i, 900 + i))
            .collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        let compressed = compress_headers_columnar(&refs).unwrap();
        assert_eq!(
            compressed[0], 0x05,
            "non-Illumina structured headers must use generic tier"
        );
        let out = decompress_headers_columnar(&compressed, refs.len()).unwrap();
        for (i, h) in refs.iter().enumerate() {
            assert_eq!(out[i].as_slice(), h.as_bytes());
        }
    }

    #[test]
    fn sra_verify_fail_routes_to_generic_at_scale() {
        let owned: Vec<String> = (1..=4000)
            .map(|i| {
                format!(
                    "@ERR{}.{i} A00296:43:HCLHLDSXX:4:2174:{}:{}",
                    i % 3,
                    100 + i,
                    200 + i
                )
            })
            .collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        let compressed = compress_headers_columnar(&refs).unwrap();
        assert_eq!(
            compressed[0], 0x05,
            "SRA verify-fail should route to generic at scale"
        );
        let out = decompress_headers_columnar(&compressed, refs.len()).unwrap();
        for (i, h) in refs.iter().enumerate() {
            assert_eq!(out[i].as_slice(), h.as_bytes());
        }
    }

    #[test]
    fn casava_verify_fail_routes_to_generic_at_scale() {
        let owned: Vec<String> = (1..=4000)
            .map(|i| {
                format!(
                    "@INSTR:43:HCLHLDSXX:4:0073:{}:{} 1:N:0:ATCG",
                    900 + i,
                    1900 + i
                )
            })
            .collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        let compressed = compress_headers_columnar(&refs).unwrap();
        assert_eq!(
            compressed[0], 0x05,
            "Casava verify-fail should route to generic at scale"
        );
        let out = decompress_headers_columnar(&compressed, refs.len()).unwrap();
        for (i, h) in refs.iter().enumerate() {
            assert_eq!(out[i].as_slice(), h.as_bytes());
        }
    }
}
