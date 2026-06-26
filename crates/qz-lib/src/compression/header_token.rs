//! Generic token-columnar header codec (blob version 0x05).
//!
//! Format-agnostic fallback used by `header_col` when a header is neither SRA nor
//! Casava: tokenize each header into a uniform leaf-token count, columnarize, and
//! BSC each column. Byte-exact by construction (tokens partition the bytes;
//! numeric tokens are guarded by `parse_canonical_u64`); a mandatory in-memory
//! verify makes the caller fall back to raw on any reconstruction mismatch.

use super::bsc;
use anyhow::Result;

// ── primitive helpers ──────────────────────────────────────────────────────
//
// This module keeps its own varint/blob/canonical-number helpers rather than
// reusing `header_col`'s by design: the codec is intentionally self-contained,
// and the varint here is u128-wide (zig-zag deltas of two `u64` values can exceed
// `i64::MAX`) with stricter overflow validation than `header_col`'s usize varint.
// The trivially-identical pieces (e.g. the u32-LE blob framing) are left
// duplicated to avoid coupling this module to `header_col`'s internals; the
// canonical-decimal predicate is NOT shared with `header_col::is_canonical_decimal`
// because that one is on the SRA/Casava byte-identical path and must not change.

/// Parse a token as a canonical `u64`: all ASCII digits, no leading zero (except
/// the single char "0"), and fits in `u64`. Returns `None` otherwise, so the
/// token is treated as a string and round-trips byte-exact (e.g. `0073`, hex,
/// 20-digit numbers).
pub(super) fn parse_canonical_u64(s: &[u8]) -> Option<u64> {
    if s.is_empty() || !s.iter().all(u8::is_ascii_digit) {
        return None;
    }
    if s.len() > 1 && s[0] == b'0' {
        return None;
    }
    // s is ASCII digits with no leading zero; parse with u64 overflow check.
    let mut v: u64 = 0;
    for &b in s {
        v = v.checked_mul(10)?.checked_add((b - b'0') as u64)?;
    }
    Some(v)
}

/// LEB128 unsigned varint over u128 (covers zig-zag u64 deltas and lengths).
fn write_uvarint(out: &mut Vec<u8>, mut v: u128) {
    while v >= 0x80 {
        out.push(((v & 0x7f) as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

/// Read a u128 LEB128 varint with overflow + length validation (rejects malformed
/// blobs rather than wrapping). Bounded to 19 bytes (shifts 0,7,…,126).
fn read_uvarint(data: &[u8], pos: &mut usize) -> Result<u128> {
    let mut v: u128 = 0;
    let mut shift = 0u32;
    for _ in 0..19 {
        if *pos >= data.len() {
            anyhow::bail!("header_token: truncated varint at {}", *pos);
        }
        let b = data[*pos];
        *pos += 1;
        let payload = (b & 0x7f) as u128;
        // On the final byte (shift == 126) only payload bits 126..=127 fit in u128,
        // so payload must be <= 3. `checked_shl(126)` would NOT catch this — it only
        // rejects shift >= 128 and silently drops bit 128 — so guard explicitly.
        if shift == 126 && payload > 3 {
            anyhow::bail!("header_token: varint overflows u128");
        }
        let chunk = payload
            .checked_shl(shift)
            .ok_or_else(|| anyhow::anyhow!("header_token: varint shift overflow"))?;
        v = v
            .checked_add(chunk)
            .ok_or_else(|| anyhow::anyhow!("header_token: varint value overflow"))?;
        if b & 0x80 == 0 {
            return Ok(v);
        }
        shift += 7;
    }
    anyhow::bail!("header_token: varint exceeds 19 bytes")
}

/// Read a varint and convert to `usize`, rejecting values that don't fit (instead
/// of a silent `as usize` truncation).
fn read_uvarint_usize(data: &[u8], pos: &mut usize) -> Result<usize> {
    let v = read_uvarint(data, pos)?;
    usize::try_from(v).map_err(|_| anyhow::anyhow!("header_token: length {v} exceeds usize"))
}

#[inline]
fn zigzag(d: i128) -> u128 {
    ((d << 1) ^ (d >> 127)) as u128
}

#[inline]
fn unzigzag(z: u128) -> i128 {
    ((z >> 1) as i128) ^ -((z & 1) as i128)
}

/// Append the ASCII decimal rendering of `v` to `out` with no leading zeros
/// ("0" for zero), without allocating a per-call `Vec`. Used on the hot decode
/// paths (numeric columns are render-on-write, never materialized per record).
fn push_u64_decimal(out: &mut Vec<u8>, v: u64) {
    let mut tmp = [0u8; 20]; // u64::MAX is 20 digits
    let mut i = tmp.len();
    let mut n = v;
    loop {
        i -= 1;
        tmp[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    out.extend_from_slice(&tmp[i..]);
}

/// Length-prefixed (u32 LE) byte blob.
fn write_blob(out: &mut Vec<u8>, blob: &[u8]) {
    out.extend_from_slice(&(blob.len() as u32).to_le_bytes());
    out.extend_from_slice(blob);
}

fn read_blob<'a>(data: &'a [u8], pos: &mut usize) -> Result<&'a [u8]> {
    if *pos + 4 > data.len() {
        anyhow::bail!("header_token: truncated blob length at {}", *pos);
    }
    let len =
        u32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]) as usize;
    *pos += 4;
    if *pos + len > data.len() {
        anyhow::bail!("header_token: truncated blob body at {}", *pos);
    }
    let s = &data[*pos..*pos + len];
    *pos += len;
    Ok(s)
}

// ── tokenizer ───────────────────────────────────────────────────────────────

/// Fixed delimiter set. Maximal runs of delimiter chars and of non-delimiter
/// (content) chars alternate as Level-1 tokens.
#[inline]
fn is_delim(b: u8) -> bool {
    matches!(
        b,
        b' ' | b'\t' | b':' | b'/' | b'=' | b'.' | b'#' | b'_' | b'-' | b';' | b','
    )
}

/// Level-1: partition `h` into maximal delimiter/content runs (byte ranges).
fn level1_ranges(h: &[u8]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < h.len() {
        let start = i;
        let delim = is_delim(h[i]);
        while i < h.len() && is_delim(h[i]) == delim {
            i += 1;
        }
        out.push((start, i));
    }
    out
}

/// Sub-split a `[start,end)` content run into maximal digit/non-digit runs.
fn digit_split(h: &[u8], start: usize, end: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut i = start;
    while i < end {
        let s = i;
        let digit = h[i].is_ascii_digit();
        while i < end && h[i].is_ascii_digit() == digit {
            i += 1;
        }
        out.push((s, i));
    }
    out
}

#[inline]
fn run_is_delim(h: &[u8], r: (usize, usize)) -> bool {
    r.0 < r.1 && is_delim(h[r.0])
}

#[inline]
fn run_is_digit(h: &[u8], r: (usize, usize)) -> bool {
    r.0 < r.1 && h[r.0].is_ascii_digit()
}

/// Tokenize a whole block of headers into a UNIFORM per-read leaf-token count.
///
/// Returns `None` (→ caller falls back to raw) when the Level-1 field structure
/// (count or delim/content pattern) varies read-to-read. Level-2 refines each
/// content field by digit/non-digit, but adopts the finer split only when it is
/// *stable* — same sub-count and sub-class pattern across every read — so
/// data-dependent hex fields stay whole instead of going ragged.
fn leaf_tokenize_block(headers: &[&[u8]]) -> Option<Vec<Vec<(usize, usize)>>> {
    use rayon::prelude::*;
    let n = headers.len();
    if n == 0 {
        return None;
    }
    // Level-1 partition per read (independent per header → parallel map).
    let l1: Vec<Vec<(usize, usize)>> = headers.par_iter().map(|h| level1_ranges(h)).collect();
    let f = l1[0].len();

    // Uniform Level-1 field count + delim/content pattern (read-only check, any
    // mismatch → None; order-independent so `.all()` over a parallel iterator is
    // byte-equivalent to the serial early-return).
    let uniform = headers.par_iter().zip(l1.par_iter()).all(|(h, r)| {
        r.len() == f
            && (0..f).all(|k| run_is_delim(h, r[k]) == run_is_delim(headers[0], l1[0][k]))
    });
    if !uniform {
        return None;
    }

    // Decide which content fields to expand (stable digit/non-digit refinement).
    // Each field's decision is independent → parallel map over the `f` fields.
    let expand: Vec<bool> = (0..f)
        .into_par_iter()
        .map(|k| {
            if run_is_delim(headers[0], l1[0][k]) {
                return false; // delimiter field — never expand
            }
            let (s0, e0) = l1[0][k];
            let sub0 = digit_split(headers[0], s0, e0);
            if sub0.len() <= 1 {
                return false; // nothing to gain
            }
            // stable ⇔ every read's content field sub-splits to the same count and
            // digit/non-digit class pattern as read 0 (`.all()` short-circuits like
            // the original labeled-break loop).
            headers.iter().zip(&l1).all(|(h, r)| {
                let (s, e) = r[k];
                let sub = digit_split(h, s, e);
                sub.len() == sub0.len()
                    && (0..sub.len())
                        .all(|j| run_is_digit(h, sub[j]) == run_is_digit(headers[0], sub0[j]))
            })
        })
        .collect();

    // Build uniform leaf ranges per read (independent per header → parallel map;
    // `collect` preserves read order).
    let out: Vec<Vec<(usize, usize)>> = headers
        .par_iter()
        .zip(l1.par_iter())
        .map(|(h, r)| {
            let mut leaves = Vec::new();
            for k in 0..f {
                if expand[k] {
                    let (s, e) = r[k];
                    leaves.extend(digit_split(h, s, e));
                } else {
                    leaves.push(r[k]);
                }
            }
            leaves
        })
        .collect();
    Some(out)
}

// ── column model ──────────────────────────────────────────────────────────

const KIND_CONST_STR: u8 = 0;
const KIND_STR_COL: u8 = 1;
const KIND_CONST_NUM: u8 = 2;
const KIND_DELTA_NUM_COL: u8 = 3;

/// Skip generic if at least this fraction of content bytes land in *varying
/// string* columns (i.e. the structure ≈ the whole header → no better than raw).
const PREGATE_VARYING_STRING_FRAC: f64 = 0.85;

/// Per-position classification produced by Phase 1, BEFORE any BSC. It is kept
/// separate from the Phase-2 `Unit` (which carries the built, to-be-compressed
/// streams) so the pre-gate can reject a block *without* having built or
/// compressed any string column — `StrCol` deliberately holds no data and is
/// re-sliced from `leaves` in Phase 2 only once the block is past the pre-gate.
enum ColPlan<'a> {
    ConstStr(&'a [u8]),
    StrCol, // values taken from leaves at serialize time
    ConstNum(u64),
    DeltaNumCol(Vec<u64>),
}

/// Compress headers with the generic token-columnar codec.
///
/// `Ok(Some(blob))` only when the block is structurable AND passes the in-memory
/// verify; `Ok(None)` on ragged input, `token_count` overflow, pre-gate rejection,
/// or verify mismatch (caller falls back to raw). `Err` only on internal failures
/// (e.g. BSC).
pub(super) fn compress_generic_token(headers: &[&str]) -> Result<Option<Vec<u8>>> {
    let n = headers.len();
    if n == 0 {
        return Ok(None);
    }
    let hb: Vec<&[u8]> = headers.iter().map(|h| h.as_bytes()).collect();
    let leaves = match leaf_tokenize_block(&hb) {
        Some(l) => l,
        None => return Ok(None),
    };
    let token_count = leaves[0].len();
    if token_count > u16::MAX as usize {
        return Ok(None);
    }
    // The record count is stored as u32 in the blob header (mirroring the
    // token_count u16 guard above); reject inputs that would silently truncate it.
    if n > u32::MAX as usize {
        return Ok(None);
    }

    // Phase 1: classify each position (no BSC yet) + pre-gate accounting.
    // Columns are independent → classify them in parallel, then reduce the byte
    // counters and collect the plans in column order. `collect` preserves order
    // and integer addition is associative, so this is byte-equivalent to the
    // serial accumulation below it was replacing.
    let col_results: Vec<(ColPlan, u64, u64)> = (0..token_count)
        .into_par_iter()
        .map(|p| {
            // Single parse pass: collect numeric values until the first non-canonical
            // token (after which the column is a string and `vals` is unused), so each
            // numeric token is parsed exactly once. `col_total` sums this column's
            // content bytes across all reads; `col_varying` is the same sum but only
            // for a *varying string* column (else 0).
            let mut vals: Vec<u64> = Vec::with_capacity(n);
            let mut all_num = true;
            let mut col_total: u64 = 0;
            for read in 0..n {
                let (s, e) = leaves[read][p];
                col_total += (e - s) as u64;
                if all_num {
                    match parse_canonical_u64(&hb[read][s..e]) {
                        Some(v) => vals.push(v),
                        None => all_num = false,
                    }
                }
            }
            if all_num {
                if vals.iter().all(|&v| v == vals[0]) {
                    (ColPlan::ConstNum(vals[0]), col_total, 0)
                } else {
                    (ColPlan::DeltaNumCol(vals), col_total, 0)
                }
            } else {
                let (s0, e0) = leaves[0][p];
                let first = &hb[0][s0..e0];
                let all_eq = (0..n).all(|read| {
                    let (s, e) = leaves[read][p];
                    &hb[read][s..e] == first
                });
                if all_eq {
                    (ColPlan::ConstStr(first), col_total, 0)
                } else {
                    // Varying string column: every read's content bytes are "varying",
                    // so its varying-byte contribution equals its total (`col_total`),
                    // matching the original's separate per-read sum.
                    (ColPlan::StrCol, col_total, col_total)
                }
            }
        })
        .collect();

    let mut plans: Vec<ColPlan> = Vec::with_capacity(token_count);
    let mut varying_string_bytes: u64 = 0;
    let mut total_bytes: u64 = 0;
    for (plan, col_total, col_varying) in col_results {
        total_bytes += col_total;
        varying_string_bytes += col_varying;
        plans.push(plan);
    }

    // Pre-gate: if almost everything is a varying string column, bail (≈ raw).
    if total_bytes > 0
        && (varying_string_bytes as f64 / total_bytes as f64) > PREGATE_VARYING_STRING_FRAC
    {
        return Ok(None);
    }

    // Phase 2: build the uncompressed stream for each VARIABLE column (no BSC
    // yet). Const columns serialize inline; variable columns are collected so
    // their BSC runs in parallel (matching the SRA/Casava per-column parallelism).
    enum Unit<'a> {
        ConstStr(&'a [u8]),
        ConstNum(u64),
        Var(u8, usize), // (kind byte, index into `streams`)
    }
    let mut units: Vec<Unit> = Vec::with_capacity(token_count);
    let mut streams: Vec<Vec<u8>> = Vec::new();
    for (p, plan) in plans.iter().enumerate() {
        match plan {
            ColPlan::ConstStr(b) => units.push(Unit::ConstStr(b)),
            ColPlan::ConstNum(v) => units.push(Unit::ConstNum(*v)),
            ColPlan::DeltaNumCol(vals) => {
                let mut stream = Vec::new();
                let mut prev: i128 = 0;
                for (i, &v) in vals.iter().enumerate() {
                    let delta = if i == 0 { v as i128 } else { v as i128 - prev };
                    write_uvarint(&mut stream, zigzag(delta));
                    prev = v as i128;
                }
                units.push(Unit::Var(KIND_DELTA_NUM_COL, streams.len()));
                streams.push(stream);
            }
            ColPlan::StrCol => {
                let mut stream = Vec::new();
                for read in 0..n {
                    let (s, e) = leaves[read][p];
                    let tok = &hb[read][s..e];
                    write_uvarint(&mut stream, tok.len() as u128);
                    stream.extend_from_slice(tok);
                }
                units.push(Unit::Var(KIND_STR_COL, streams.len()));
                streams.push(stream);
            }
        }
    }

    // Phase 3: compress the variable-column streams in parallel.
    use rayon::prelude::*;
    let comps: Vec<Vec<u8>> = streams
        .par_iter()
        .map(|s| bsc::compress_parallel_adaptive(s))
        .collect::<Result<Vec<_>>>()?;

    // Phase 4: serialize in position order.
    let mut blob = Vec::new();
    blob.push(0x05);
    blob.extend_from_slice(&(n as u32).to_le_bytes());
    blob.extend_from_slice(&(token_count as u16).to_le_bytes());
    for unit in &units {
        match unit {
            Unit::ConstStr(b) => {
                blob.push(KIND_CONST_STR);
                write_uvarint(&mut blob, b.len() as u128);
                blob.extend_from_slice(b);
            }
            Unit::ConstNum(v) => {
                blob.push(KIND_CONST_NUM);
                write_uvarint(&mut blob, *v as u128);
            }
            Unit::Var(kind, idx) => {
                blob.push(*kind);
                write_blob(&mut blob, &comps[*idx]);
            }
        }
    }

    // Verify: reconstruct from the in-memory plans (no BSC round-trip) and compare.
    if !verify_reconstruct(&plans, &leaves, &hb, n, token_count) {
        return Ok(None);
    }
    Ok(Some(blob))
}

/// Reconstruct each header from the in-memory column plans and compare to input.
fn verify_reconstruct(
    plans: &[ColPlan],
    leaves: &[Vec<(usize, usize)>],
    hb: &[&[u8]],
    n: usize,
    token_count: usize,
) -> bool {
    let mut buf = Vec::new();
    for read in 0..n {
        buf.clear();
        for p in 0..token_count {
            match &plans[p] {
                ColPlan::ConstStr(bytes) => buf.extend_from_slice(bytes),
                ColPlan::ConstNum(v) => push_u64_decimal(&mut buf, *v),
                ColPlan::DeltaNumCol(vals) => push_u64_decimal(&mut buf, vals[read]),
                ColPlan::StrCol => {
                    let (s, e) = leaves[read][p];
                    buf.extend_from_slice(&hb[read][s..e]);
                }
            }
        }
        if buf != hb[read] {
            return false;
        }
    }
    true
}

// ── decode ──────────────────────────────────────────────────────────────────

/// A decoded column in a form needing NO per-record allocation: constants are
/// stored once; a variable string column keeps its decompressed `stream` plus
/// per-record `spans` (slice on demand); a numeric column is one contiguous
/// `Vec<u64>` rendered to decimal on write. This is what lets the direct-to-varint
/// path avoid the millions of tiny `Vec`s the per-record form would create.
enum DecodedCol {
    ConstStr(Vec<u8>),
    ConstNum(u64),
    StrCol {
        stream: Vec<u8>,
        spans: Vec<(usize, usize)>,
    },
    NumCol(Vec<u64>),
}

/// First-pass column descriptor: constants resolved immediately; variable columns
/// capture their compressed blob slice for parallel decompression in pass 2.
enum RawCol<'a> {
    ConstStr(&'a [u8]),
    ConstNum(u64),
    StrBlob(&'a [u8]),
    NumBlob(&'a [u8]),
}

/// Parse the blob (already past the 0x05 version byte) into decoded columns,
/// validating `num_reads`, exact stream consumption, and numeric ranges.
fn decode_columns(data: &[u8], num_reads: usize) -> Result<(usize, Vec<DecodedCol>)> {
    if data.len() < 6 {
        anyhow::bail!("header_token: blob too short for header");
    }
    let stored_n = u32::from_le_bytes([data[0], data[1], data[2], data[3]]) as usize;
    anyhow::ensure!(
        stored_n == num_reads,
        "header_token: stored num_reads {stored_n} != framing {num_reads}"
    );
    let token_count = u16::from_le_bytes([data[4], data[5]]) as usize;
    let mut pos = 6usize;

    // Pass 1 (serial, cheap): slice out each column's bytes.
    let mut raw: Vec<RawCol> = Vec::with_capacity(token_count);
    for _ in 0..token_count {
        if pos >= data.len() {
            anyhow::bail!("header_token: truncated at column kind");
        }
        let kind = data[pos];
        pos += 1;
        match kind {
            KIND_CONST_STR => {
                let len = read_uvarint_usize(data, &mut pos)?;
                if pos + len > data.len() {
                    anyhow::bail!("header_token: truncated const_str");
                }
                raw.push(RawCol::ConstStr(&data[pos..pos + len]));
                pos += len;
            }
            KIND_CONST_NUM => {
                let v = read_uvarint(data, &mut pos)?;
                let v = u64::try_from(v)
                    .map_err(|_| anyhow::anyhow!("header_token: const_num exceeds u64"))?;
                raw.push(RawCol::ConstNum(v));
            }
            KIND_STR_COL => raw.push(RawCol::StrBlob(read_blob(data, &mut pos)?)),
            KIND_DELTA_NUM_COL => raw.push(RawCol::NumBlob(read_blob(data, &mut pos)?)),
            other => anyhow::bail!("header_token: unknown column kind {other}"),
        }
    }
    anyhow::ensure!(
        pos == data.len(),
        "header_token: trailing bytes after columns"
    );

    // Pass 2 (parallel): decompress + parse the variable columns.
    use rayon::prelude::*;
    let cols: Vec<DecodedCol> = raw
        .par_iter()
        .map(|rc| -> Result<DecodedCol> {
            match rc {
                RawCol::ConstStr(b) => Ok(DecodedCol::ConstStr(b.to_vec())),
                RawCol::ConstNum(v) => Ok(DecodedCol::ConstNum(*v)),
                RawCol::StrBlob(comp) => {
                    let stream = bsc::decompress_parallel(comp)?;
                    // Each entry consumes >=1 stream byte (its length varint), so a
                    // valid stream has >= num_reads bytes; cap the pre-allocation at
                    // the stream size so a crafted huge `num_reads` can't OOM-abort
                    // before per-entry bounds are hit.
                    let mut spans = Vec::with_capacity(num_reads.min(stream.len() + 1));
                    let mut sp = 0usize;
                    for _ in 0..num_reads {
                        let len = read_uvarint_usize(&stream, &mut sp)?;
                        let end = sp
                            .checked_add(len)
                            .filter(|&e| e <= stream.len())
                            .ok_or_else(|| {
                                anyhow::anyhow!("header_token: truncated str_col entry")
                            })?;
                        spans.push((sp, end));
                        sp = end;
                    }
                    anyhow::ensure!(
                        sp == stream.len(),
                        "header_token: str_col stream not fully consumed"
                    );
                    Ok(DecodedCol::StrCol { stream, spans })
                }
                RawCol::NumBlob(comp) => {
                    let stream = bsc::decompress_parallel(comp)?;
                    // Cap pre-allocation at the stream size (each entry is >=1 varint
                    // byte) — see the StrBlob arm for the rationale.
                    let mut vals = Vec::with_capacity(num_reads.min(stream.len() + 1));
                    let mut sp = 0usize;
                    let mut prev: i128 = 0;
                    for i in 0..num_reads {
                        let z = read_uvarint(&stream, &mut sp)?;
                        let delta = unzigzag(z);
                        let val_i = if i == 0 {
                            delta
                        } else {
                            prev.checked_add(delta).ok_or_else(|| {
                                anyhow::anyhow!("header_token: delta accumulation overflow")
                            })?
                        };
                        let v = u64::try_from(val_i)
                            .map_err(|_| anyhow::anyhow!("header_token: value out of u64 range"))?;
                        vals.push(v);
                        prev = val_i;
                    }
                    anyhow::ensure!(
                        sp == stream.len(),
                        "header_token: delta_num stream not fully consumed"
                    );
                    Ok(DecodedCol::NumCol(vals))
                }
            }
        })
        .collect::<Result<Vec<_>>>()?;

    Ok((token_count, cols))
}

/// Append record `i`'s reconstructed bytes (concatenated columns) to `out`.
/// Numbers are rendered on write; strings are sliced from the flat column stream
/// — no per-record allocation.
#[inline]
fn render_record(out: &mut Vec<u8>, cols: &[DecodedCol], token_count: usize, i: usize) {
    for p in 0..token_count {
        match &cols[p] {
            DecodedCol::ConstStr(b) => out.extend_from_slice(b),
            DecodedCol::ConstNum(v) => push_u64_decimal(out, *v),
            DecodedCol::StrCol { stream, spans } => {
                let (s, e) = spans[i];
                out.extend_from_slice(&stream[s..e]);
            }
            DecodedCol::NumCol(vals) => push_u64_decimal(out, vals[i]),
        }
    }
}

/// Decode to per-record header byte vectors (materialize-all path; one output
/// `Vec` per record is inherent to the return type).
///
/// Note: we do NOT short-circuit `num_reads == 0` before parsing — `decode_columns`
/// must still validate the stored `num_reads`, token count, exact stream
/// consumption, and trailing bytes, so a garbage `0x05` body is rejected rather
/// than silently treated as empty. A well-formed empty body decodes to `[]`.
pub(super) fn decompress_generic_token(data: &[u8], num_reads: usize) -> Result<Vec<Vec<u8>>> {
    let (token_count, cols) = decode_columns(data, num_reads)?;
    use rayon::prelude::*;
    Ok((0..num_reads)
        .into_par_iter()
        .map(|i| {
            let mut out = Vec::new();
            render_record(&mut out, &cols, token_count, i);
            out
        })
        .collect())
}

/// Decode directly into the varint-framed stream `([varint len][bytes])*` that the
/// bounded-streaming consumer expects. Each chunk worker reuses ONE scratch `rec`
/// buffer (cleared per record — not a per-record allocation) that `render_record`
/// fills (numbers rendered on write, strings sliced from the flat column streams),
/// then writes `len + rec` into its `buf`; chunk bufs are concatenated. Mirrors
/// `header_col::decompress_sra_into_varint`.
pub(super) fn decompress_generic_token_to_varint(data: &[u8], num_reads: usize) -> Result<Vec<u8>> {
    // No `num_reads == 0` short-circuit — `decode_columns` must validate the body
    // (see `decompress_generic_token`). A well-formed empty body yields an empty stream.
    let (token_count, cols) = decode_columns(data, num_reads)?;
    use rayon::prelude::*;
    const CHUNK: usize = 8192;
    let chunks: Vec<Vec<u8>> = (0..num_reads)
        .into_par_iter()
        .step_by(CHUNK)
        .map(|start| {
            let end = (start + CHUNK).min(num_reads);
            let mut buf = Vec::new();
            let mut rec = Vec::new(); // reused per worker; cleared each record
            for i in start..end {
                rec.clear();
                render_record(&mut rec, &cols, token_count, i);
                write_uvarint(&mut buf, rec.len() as u128);
                buf.extend_from_slice(&rec);
            }
            buf
        })
        .collect();
    Ok(chunks.concat())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_canonical_u64_rules() {
        assert_eq!(parse_canonical_u64(b"0"), Some(0));
        assert_eq!(parse_canonical_u64(b"72185641"), Some(72185641));
        assert_eq!(parse_canonical_u64(b"18446744073709551615"), Some(u64::MAX));
        assert_eq!(parse_canonical_u64(b""), None); // empty
        assert_eq!(parse_canonical_u64(b"0073"), None); // leading zero
        assert_eq!(parse_canonical_u64(b"00"), None); // leading zero
        assert_eq!(parse_canonical_u64(b"12a"), None); // non-digit
        assert_eq!(parse_canonical_u64(b"18446744073709551616"), None); // u64 overflow
    }

    #[test]
    fn uvarint_roundtrip_full_range() {
        for v in [
            0u128,
            1,
            127,
            128,
            300,
            u64::MAX as u128,
            (u64::MAX as u128) * 2,
        ] {
            let mut b = Vec::new();
            write_uvarint(&mut b, v);
            let mut pos = 0;
            assert_eq!(read_uvarint(&b, &mut pos).unwrap(), v);
            assert_eq!(pos, b.len());
        }
    }

    #[test]
    fn zigzag_roundtrip_spans_u64_delta() {
        // A u64-u64 delta can be as large as +/-(u64::MAX), exceeding i64::MAX.
        for d in [
            0i128,
            1,
            -1,
            i64::MAX as i128,
            i64::MIN as i128,
            u64::MAX as i128,
            -(u64::MAX as i128),
        ] {
            assert_eq!(unzigzag(zigzag(d)), d);
        }
    }

    #[test]
    fn read_uvarint_rejects_malformed() {
        // All-continuation bytes (never terminates) → length-bounded reject.
        let bad = vec![0x80u8; 25];
        let mut pos = 0;
        assert!(read_uvarint(&bad, &mut pos).is_err());
        // Truncated (continuation bit set at end of input).
        let trunc = vec![0x80u8, 0x80];
        let mut pos = 0;
        assert!(read_uvarint(&trunc, &mut pos).is_err());
        // 19th-byte overflow: 18 continuation bytes + terminal payload 0x04 at
        // shift 126 would shift past bit 127 → must be rejected, not silently 0.
        let mut overflow = vec![0x80u8; 18];
        overflow.push(0x04);
        let mut pos = 0;
        assert!(read_uvarint(&overflow, &mut pos).is_err());
        // Boundary: terminal payload 0x03 at shift 126 DOES fit (bits 126,127).
        let mut ok = vec![0x80u8; 18];
        ok.push(0x03);
        let mut pos = 0;
        assert_eq!(read_uvarint(&ok, &mut pos).unwrap(), 3u128 << 126);
    }

    #[test]
    fn push_u64_decimal_matches_to_string() {
        for v in [0u64, 7, 10, 99, 100, 12345, u64::MAX] {
            let mut b = Vec::new();
            push_u64_decimal(&mut b, v);
            assert_eq!(b, v.to_string().into_bytes());
        }
    }

    // Helper: render leaf ranges of read `i` back to byte tokens for assertions.
    fn leaves_of<'a>(h: &'a [u8], ranges: &[(usize, usize)]) -> Vec<&'a [u8]> {
        ranges.iter().map(|&(s, e)| &h[s..e]).collect()
    }

    #[test]
    fn tokenize_numbered_expands_stable_alpha_digit() {
        let hs: Vec<&[u8]> = vec![b"@read1", b"@read2", b"@read10"];
        let leaves = leaf_tokenize_block(&hs).expect("uniform");
        // "@read" + number, uniform 2 leaves per read.
        assert_eq!(leaves[0].len(), 2);
        assert_eq!(leaves_of(hs[0], &leaves[0]), vec![&b"@read"[..], &b"1"[..]]);
        assert_eq!(
            leaves_of(hs[2], &leaves[2]),
            vec![&b"@read"[..], &b"10"[..]]
        );
    }

    #[test]
    fn tokenize_keeps_ont_uuid_whole() {
        // Hex UUIDs sub-split raggedly (data-dependent digit positions), so the
        // stable-refinement test must reject the split and keep them whole.
        let hs: Vec<&[u8]> = vec![b"@a1b2 read=1", b"@cd34 read=2", b"@9f0e read=3"];
        let leaves = leaf_tokenize_block(&hs).expect("uniform fields");
        let toks0 = leaves_of(hs[0], &leaves[0]);
        assert!(
            toks0.contains(&&b"@a1b2"[..]),
            "uuid-ish id kept whole: {:?}",
            toks0
        );
        assert!(toks0.contains(&&b"1"[..]), "read number extracted");
        assert_eq!(leaves[0].len(), leaves[1].len());
        assert_eq!(leaves[1].len(), leaves[2].len());
    }

    #[test]
    fn tokenize_ragged_field_count_returns_none() {
        // Different number of space-delimited fields → not representable by A1.
        let hs: Vec<&[u8]> = vec![b"@x a b", b"@x a"];
        assert!(leaf_tokenize_block(&hs).is_none());
    }

    #[test]
    fn tokenize_concatenation_is_identity() {
        let hs: Vec<&[u8]> = vec![b"@SRR1.1 A:1:2:3:4", b"@SRR1.2 A:1:2:9:8"];
        let leaves = leaf_tokenize_block(&hs).expect("uniform");
        for (i, h) in hs.iter().enumerate() {
            let cat: Vec<u8> = leaves[i]
                .iter()
                .flat_map(|&(s, e)| h[s..e].iter().copied())
                .collect();
            assert_eq!(&cat, h, "leaf tokens must concatenate to the original");
        }
    }

    #[test]
    fn compress_emits_blob_for_structured_block() {
        // Numbered headers → structured → Some(blob), version 0x05.
        let owned: Vec<String> = (1..=2000).map(|i| format!("@read{i}")).collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        let blob = compress_generic_token(&refs).unwrap().expect("structured");
        assert_eq!(blob[0], 0x05);
    }

    #[test]
    fn compress_pregate_skips_high_entropy_ids() {
        // Distinct random-ish single-token IDs with no shared structure: almost
        // all bytes are a varying string column ≈ raw → pre-gate returns None.
        let owned: Vec<String> = (0..2000usize)
            .map(|i| {
                format!(
                    "@{:016x}{:016x}",
                    i.wrapping_mul(0x9e3779b97f4a7c15usize),
                    i
                )
            })
            .collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        assert!(
            compress_generic_token(&refs).unwrap().is_none(),
            "high-entropy id-only block should pre-gate to None"
        );
    }

    fn roundtrip_assert(headers: &[&str]) {
        let blob = compress_generic_token(headers).unwrap();
        let blob = match blob {
            Some(b) => b,
            None => return, // pre-gated/ragged: nothing to assert here
        };
        assert_eq!(blob[0], 0x05);
        // decode receives the slice AFTER the version byte (matches header_col dispatch).
        let out = decompress_generic_token(&blob[1..], headers.len()).unwrap();
        assert_eq!(out.len(), headers.len());
        for (i, h) in headers.iter().enumerate() {
            assert_eq!(out[i], h.as_bytes(), "mismatch at read {i}: {h:?}");
        }
    }

    #[test]
    fn roundtrip_numbered() {
        let owned: Vec<String> = (1..=3000).map(|i| format!("@read{i}")).collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        roundtrip_assert(&refs);
    }

    #[test]
    fn roundtrip_nanopore_like() {
        let owned: Vec<String> = (1..=2000)
            .map(|i| format!("@a1b2 runid=abc read={i} ch={}", i % 512))
            .collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        roundtrip_assert(&refs);
    }

    #[test]
    fn roundtrip_delta_exceeds_i64_max() {
        // Adjacent values whose delta exceeds i64::MAX, to exercise i128 delta.
        let owned = [
            "@x 0".to_string(),
            format!("@x {}", u64::MAX),
            "@x 1".to_string(),
        ];
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        roundtrip_assert(&refs);
    }

    #[test]
    fn roundtrip_mixed_numeric_and_string_at_a_position() {
        // Position after "@x " is canonical-numeric in some reads and a leading-zero
        // string in others → must collapse to a single string column, losslessly.
        let owned = [
            "@x 5".to_string(),
            "@x 007".to_string(),
            "@x 42".to_string(),
            "@x 9".to_string(),
        ];
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        roundtrip_assert(&refs);
    }

    #[test]
    fn roundtrip_trailing_delimiter() {
        // Trailing delimiter must be preserved (tokens still partition the bytes).
        let owned: Vec<String> = (1..=2000).map(|i| format!("@read{i}/")).collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        roundtrip_assert(&refs);
    }

    #[test]
    fn decode_rejects_num_reads_mismatch() {
        let owned: Vec<String> = (1..=500).map(|i| format!("@read{i}")).collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        let blob = compress_generic_token(&refs).unwrap().unwrap();
        // Wrong num_reads must error, not silently mis-decode.
        assert!(decompress_generic_token(&blob[1..], 499).is_err());
    }

    #[test]
    fn decode_rejects_malformed_blobs() {
        let owned: Vec<String> = (1..=500).map(|i| format!("@read{i}")).collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        let blob = compress_generic_token(&refs).unwrap().unwrap();
        let body = &blob[1..]; // past the 0x05 version byte

        // Truncated body → error, never a panic.
        assert!(decompress_generic_token(&body[..body.len() / 2], refs.len()).is_err());

        // Trailing garbage after the columns → error (exact-consumption check).
        let mut trailing = body.to_vec();
        trailing.push(0xAB);
        assert!(decompress_generic_token(&trailing, refs.len()).is_err());

        // Unknown column kind → error. (Byte 6 is the first column's kind byte.)
        let mut bad_kind = body.to_vec();
        bad_kind[6] = 0x7F;
        assert!(decompress_generic_token(&bad_kind, refs.len()).is_err());

        // Empty body → error.
        assert!(decompress_generic_token(&[], refs.len()).is_err());

        // num_reads == 0 must still validate the body (no early-return bypass):
        // stored_n=0, token_count=0, but a trailing byte → trailing-bytes error.
        assert!(decompress_generic_token(&[0, 0, 0, 0, 0, 0, 0xFF], 0).is_err());
        // A well-formed empty body (stored_n=0, token_count=0, no columns) → [].
        assert_eq!(
            decompress_generic_token(&[0, 0, 0, 0, 0, 0], 0)
                .unwrap()
                .len(),
            0
        );
    }

    #[test]
    fn varint_stream_matches_reference() {
        let owned: Vec<String> = (1..=1500)
            .map(|i| format!("@read{i} x={}", i * 3))
            .collect();
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        let blob = compress_generic_token(&refs).unwrap().unwrap();
        let recs = decompress_generic_token(&blob[1..], refs.len()).unwrap();
        // Reference varint stream: [varint len][bytes] per record.
        let mut reference = Vec::new();
        for r in &recs {
            write_uvarint(&mut reference, r.len() as u128);
            reference.extend_from_slice(r);
        }
        let got = decompress_generic_token_to_varint(&blob[1..], refs.len()).unwrap();
        assert_eq!(got, reference);
    }

    fn r#gen(name: &str, n: usize) -> Vec<String> {
        match name {
            "nanopore" => (1..=n)
                .map(|i| format!("@d3{:04x}9a read={i} ch={} runid=abcd", i % 4096, i % 512))
                .collect(),
            "pacbio" => (1..=n)
                .map(|i| format!("@m54238_180628/{}/0_{}", 4194300 + i, 900 + i))
                .collect(),
            "bgi" => (1..=n)
                .map(|i| format!("@V300047L1C001R{:03}_{}/1", i % 1000, i))
                .collect(),
            "numbered" => (1..=n).map(|i| format!("@read{i}")).collect(),
            // degenerate: high-entropy id-only → pre-gate should fire.
            "idonly" => (1..=n)
                .map(|i| {
                    format!(
                        "@{:016x}{:016x}",
                        (i as u64).wrapping_mul(0x9e3779b97f4a7c15),
                        i
                    )
                })
                .collect(),
            _ => unreachable!(),
        }
    }

    #[test]
    #[ignore = "benchmark; run with --ignored --nocapture"]
    fn bench_header_stage() {
        use std::time::Instant;
        let n: usize = std::env::var("QZ_BENCH_N")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(2_500_000);
        for name in ["nanopore", "pacbio", "bgi", "numbered", "idonly"] {
            let owned = r#gen(name, n);
            let refs: Vec<&str> = owned.iter().map(String::as_str).collect();

            // Generic-only timing (the added work) — None = pre-gate/ragged/verify-fail.
            let t = Instant::now();
            let generic = compress_generic_token(&refs).unwrap();
            let gen_t = t.elapsed().as_secs_f64();

            // Raw-only baseline timing (today's fallback for non-SRA/Casava data).
            let t = Instant::now();
            let raw = super::super::header_col::compress_raw_fallback(&refs).unwrap();
            let raw_t = t.elapsed().as_secs_f64();

            let decision = match &generic {
                None => "PRE-GATE/none".to_string(),
                Some(b) if b.len() <= raw.len() => format!("generic-wins {}B", b.len()),
                Some(b) => format!("size-gate-loss {}B (>{}B raw)", b.len(), raw.len()),
            };
            println!(
                "{name:10} n={n}  generic={gen_t:.2}s  raw={raw_t:.2}s  raw={}B  -> {decision}",
                raw.len()
            );
        }
    }

    #[test]
    fn pregate_fires_on_idonly_uniform_input() {
        // idonly is uniform (single token), so `None` here can only be the pre-gate
        // (not ragged, and it would verify) — proves the pre-gate, not a size-gate loss.
        let owned = r#gen("idonly", 5000);
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        assert!(
            compress_generic_token(&refs).unwrap().is_none(),
            "pre-gate must skip high-entropy id-only input"
        );
    }
}
