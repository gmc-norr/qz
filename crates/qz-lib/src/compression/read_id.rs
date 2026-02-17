//! Template-based read ID compression
//!
//! Illumina read IDs follow a structured format:
//! @INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y [COMMENT]
//!
//! Most fields are constant across a file, so we can:
//! 1. Extract the common template once
//! 2. Store only the varying parts (X:Y coordinates) as integers
//! 3. Achieve 80-90% reduction in header size

use anyhow::{Context, Result};
use std::io::{Cursor, Read, Write};

/// Write variable-length integer (same as in mod.rs)
fn write_varint<W: Write>(writer: &mut W, mut value: usize) -> std::io::Result<()> {
    while value >= 0x80 {
        writer.write_all(&[((value & 0x7F) | 0x80) as u8])?;
        value >>= 7;
    }
    writer.write_all(&[value as u8])
}

/// Read variable-length integer (same as in mod.rs)
fn read_varint(data: &[u8], offset: &mut usize) -> Option<usize> {
    let mut value = 0usize;
    let mut shift = 0;

    loop {
        if *offset >= data.len() {
            return None;
        }

        let byte = data[*offset];
        *offset += 1;

        value |= ((byte & 0x7F) as usize) << shift;

        if byte & 0x80 == 0 {
            return Some(value);
        }

        shift += 7;
    }
}

#[derive(Debug, Clone)]
pub struct ReadIdTemplate {
    /// Common prefix (e.g., "@INSTRUMENT:RUN:FLOWCELL:LANE:TILE:")
    pub prefix: String,
    /// Whether read IDs have space-separated comments
    pub has_comment: bool,
    /// Common comment suffix (if all reads have the same comment)
    pub common_comment: Option<String>,
}

#[derive(Debug)]
pub struct EncodedReadIds {
    pub template: ReadIdTemplate,
    /// Encoded coordinates and optional comments
    pub encoded_data: Vec<u8>,
}

/// Analyze read IDs to extract common template
pub fn analyze_read_ids(read_ids: &[String]) -> Result<ReadIdTemplate> {
    if read_ids.is_empty() {
        anyhow::bail!("Cannot analyze empty read ID list");
    }

    // Parse first read ID to get structure
    let first_id = &read_ids[0];

    // Remove @ prefix if present
    let id_without_at = first_id.strip_prefix('@').unwrap_or(first_id);

    // Check if there's a comment (space-separated)
    let has_comment = id_without_at.contains(' ');
    let base_id = if has_comment {
        id_without_at.split(' ').next().unwrap_or(id_without_at)
    } else {
        id_without_at
    };

    // Split by colons to find fields
    let fields: Vec<&str> = base_id.split(':').collect();

    if fields.len() < 7 {
        // Not standard Illumina format, fall back to no template
        return Ok(ReadIdTemplate {
            prefix: String::new(),
            has_comment: false,
            common_comment: None,
        });
    }

    // Standard Illumina format: INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y
    // Extract template (everything except X:Y)
    let template_parts = &fields[..5]; // INSTRUMENT:RUN:FLOWCELL:LANE:TILE
    let prefix = format!("@{}:", template_parts.join(":"));

    // Check if all comments are identical
    let common_comment = if has_comment {
        let first_comment = id_without_at.split(' ').nth(1);
        if let Some(first_comm) = first_comment {
            // Check if all reads have the same comment
            let all_same = read_ids.iter().all(|id| {
                let id_without_at = id.strip_prefix('@').unwrap_or(id);
                id_without_at.split(' ').nth(1) == Some(first_comm)
            });
            if all_same {
                Some(first_comm.to_string())
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    };

    Ok(ReadIdTemplate {
        prefix,
        has_comment,
        common_comment,
    })
}

/// Check if a delta fits in 2 bits (values -1, 0, 1)
fn fits_in_2bits(delta: i32) -> bool {
    delta >= -1 && delta <= 1
}

/// Pack small deltas (2 coordinate pairs per byte)
/// Each delta is encoded as: 00 = -1, 01 = 0, 10 = 1
fn pack_small_deltas(deltas: &[(i32, i32)]) -> Vec<u8> {
    let mut packed = Vec::with_capacity((deltas.len() + 1) / 2);

    for chunk in deltas.chunks(2) {
        let mut byte = 0u8;

        // First pair (lower 4 bits)
        if let Some((dx, dy)) = chunk.get(0) {
            let x_bits = (dx + 1) as u8;
            let y_bits = (dy + 1) as u8;
            byte |= (x_bits & 0b11) | ((y_bits & 0b11) << 2);
        }

        // Second pair (upper 4 bits)
        if let Some((dx, dy)) = chunk.get(1) {
            let x_bits = (dx + 1) as u8;
            let y_bits = (dy + 1) as u8;
            byte |= ((x_bits & 0b11) << 4) | ((y_bits & 0b11) << 6);
        }

        packed.push(byte);
    }

    packed
}

/// Unpack small deltas
fn unpack_small_deltas(packed: &[u8], num_pairs: usize) -> Vec<(i32, i32)> {
    let mut deltas = Vec::with_capacity(num_pairs);

    for &byte in packed.iter() {
        if deltas.len() >= num_pairs {
            break;
        }

        // First pair (lower 4 bits)
        let x_bits = byte & 0b11;
        let y_bits = (byte >> 2) & 0b11;
        let dx = (x_bits as i32) - 1;
        let dy = (y_bits as i32) - 1;
        deltas.push((dx, dy));

        if deltas.len() >= num_pairs {
            break;
        }

        // Second pair (upper 4 bits)
        let x_bits = (byte >> 4) & 0b11;
        let y_bits = (byte >> 6) & 0b11;
        let dx = (x_bits as i32) - 1;
        let dy = (y_bits as i32) - 1;
        deltas.push((dx, dy));
    }

    deltas.truncate(num_pairs);
    deltas
}

/// Encode read IDs using template
pub fn encode_read_ids(read_ids: &[String], template: &ReadIdTemplate) -> Result<Vec<u8>> {
    let mut encoded = Vec::new();

    if template.prefix.is_empty() {
        // No template compression, just store raw IDs with length prefixes
        for id in read_ids {
            let bytes = id.as_bytes();
            encoded.write_all(&(bytes.len() as u16).to_le_bytes())?;
            encoded.write_all(bytes)?;
        }
        return Ok(encoded);
    }

    // First pass: extract coordinates and compute deltas
    let mut coordinates = Vec::new();
    let mut comments = Vec::new();

    for id in read_ids {
        let id_without_at = id.strip_prefix('@').unwrap_or(id);

        // Split off comment if present
        let (base_id, comment) = if template.has_comment {
            if let Some(space_pos) = id_without_at.find(' ') {
                let (base, comm) = id_without_at.split_at(space_pos);
                (base, Some(&comm[1..])) // Skip the space
            } else {
                (id_without_at, None)
            }
        } else {
            (id_without_at, None)
        };

        // Remove template prefix
        let template_without_at = template.prefix.strip_prefix('@').unwrap_or(&template.prefix);
        let remaining = base_id.strip_prefix(template_without_at)
            .context("Read ID doesn't match template")?;

        // Parse X:Y coordinates
        let coords: Vec<&str> = remaining.split(':').collect();
        if coords.len() != 2 {
            anyhow::bail!("Expected X:Y coordinates, got: {}", remaining);
        }

        let x: u32 = coords[0].parse()
            .context("Failed to parse X coordinate")?;
        let y: u32 = coords[1].parse()
            .context("Failed to parse Y coordinate")?;

        coordinates.push((x, y));
        comments.push(comment);
    }

    // Compute deltas
    let mut deltas = Vec::new();
    let mut prev_x = 0i32;
    let mut prev_y = 0i32;

    for &(x, y) in &coordinates {
        let delta_x = (x as i32) - prev_x;
        let delta_y = (y as i32) - prev_y;
        deltas.push((delta_x, delta_y));
        prev_x = x as i32;
        prev_y = y as i32;
    }

    // Check if all deltas fit in 2 bits
    let use_packed = deltas.iter().all(|(dx, dy)| fits_in_2bits(*dx) && fits_in_2bits(*dy));

    // Write encoding flag (0 = varint, 1 = packed)
    encoded.write_all(&[if use_packed { 1 } else { 0 }])?;

    if use_packed {
        // Pack deltas into 4 pairs per byte
        let packed = pack_small_deltas(&deltas);
        encoded.write_all(&packed)?;
    } else {
        // Use zigzag varint encoding
        for (delta_x, delta_y) in deltas {
            let zigzag_x = ((delta_x << 1) ^ (delta_x >> 31)) as usize;
            let zigzag_y = ((delta_y << 1) ^ (delta_y >> 31)) as usize;
            write_varint(&mut encoded, zigzag_x)?;
            write_varint(&mut encoded, zigzag_y)?;
        }
    }

    // Encode comments only if not common
    if template.has_comment && template.common_comment.is_none() {
        for comment in comments {
            if let Some(comm) = comment {
                encoded.write_all(&(comm.len() as u16).to_le_bytes())?;
                encoded.write_all(comm.as_bytes())?;
            } else {
                // No comment for this read, write length 0
                encoded.write_all(&0u16.to_le_bytes())?;
            }
        }
    }

    Ok(encoded)
}

/// Decode read IDs from template encoding
pub fn decode_read_ids(
    encoded_data: &[u8],
    template: &ReadIdTemplate,
    num_reads: usize,
) -> Result<Vec<String>> {
    let mut cursor = Cursor::new(encoded_data);
    let mut read_ids = Vec::with_capacity(num_reads);

    if template.prefix.is_empty() {
        // No template compression, read raw IDs
        for _ in 0..num_reads {
            let mut len_buf = [0u8; 2];
            cursor.read_exact(&mut len_buf)?;
            let len = u16::from_le_bytes(len_buf) as usize;
            let mut bytes = vec![0u8; len];
            cursor.read_exact(&mut bytes)?;
            read_ids.push(String::from_utf8_lossy(&bytes).to_string());
        }
        return Ok(read_ids);
    }

    // Read encoding flag (0 = varint, 1 = packed)
    let mut flag_buf = [0u8; 1];
    cursor.read_exact(&mut flag_buf)?;
    let encoding_flag = flag_buf[0];
    let use_packed = encoding_flag == 1;

    // Decode deltas
    let deltas = if use_packed {
        // Calculate number of packed bytes needed
        let num_bytes = (num_reads + 1) / 2; // 2 pairs per byte
        let mut packed_bytes = vec![0u8; num_bytes];
        cursor.read_exact(&mut packed_bytes)?;
        unpack_small_deltas(&packed_bytes, num_reads)
    } else {
        // Use zigzag varint decoding
        let mut deltas = Vec::new();
        let data = encoded_data;
        let mut offset = cursor.position() as usize;

        for _ in 0..num_reads {
            let zigzag_x = read_varint(data, &mut offset)
                .context("Failed to read X coordinate")?;
            let zigzag_y = read_varint(data, &mut offset)
                .context("Failed to read Y coordinate")?;

            // Decode zigzag
            let delta_x = ((zigzag_x >> 1) as i32) ^ (-((zigzag_x & 1) as i32));
            let delta_y = ((zigzag_y >> 1) as i32) ^ (-((zigzag_y & 1) as i32));

            deltas.push((delta_x, delta_y));
        }

        cursor.set_position(offset as u64);
        deltas
    };

    // Reconstruct coordinates from deltas
    let mut prev_x = 0i32;
    let mut prev_y = 0i32;

    for (delta_x, delta_y) in deltas {
        let x = prev_x + delta_x;
        let y = prev_y + delta_y;

        prev_x = x;
        prev_y = y;

        // Reconstruct base ID
        let base_id = format!("{}{}:{}", template.prefix, x, y);

        // Read comment
        let full_id = if template.has_comment {
            if let Some(ref common_comm) = template.common_comment {
                // Use common comment
                format!("{} {}", base_id, common_comm)
            } else {
                // Read per-read comment
                let mut clen_buf = [0u8; 2];
                cursor.read_exact(&mut clen_buf)?;
                let comment_len = u16::from_le_bytes(clen_buf) as usize;

                if comment_len > 0 {
                    let mut comment_bytes = vec![0u8; comment_len];
                    cursor.read_exact(&mut comment_bytes)?;
                    let comment = String::from_utf8_lossy(&comment_bytes);
                    format!("{} {}", base_id, comment)
                } else {
                    base_id
                }
            }
        } else {
            base_id
        };

        read_ids.push(full_id);
    }

    Ok(read_ids)
}

/// Compress read IDs with template encoding
pub fn compress_read_ids(read_ids: &[String]) -> Result<EncodedReadIds> {
    let template = analyze_read_ids(read_ids)?;
    let encoded_data = encode_read_ids(read_ids, &template)?;

    Ok(EncodedReadIds {
        template,
        encoded_data,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_illumina_format() {
        let read_ids = vec![
            "@INSTRUMENT:123:FLOWCELL:1:2101:10000:1000".to_string(),
            "@INSTRUMENT:123:FLOWCELL:1:2101:10001:1001".to_string(),
            "@INSTRUMENT:123:FLOWCELL:1:2101:10002:1002".to_string(),
        ];

        let template = analyze_read_ids(&read_ids).unwrap();
        assert_eq!(template.prefix, "@INSTRUMENT:123:FLOWCELL:1:2101:");
        assert!(!template.has_comment);

        let encoded = encode_read_ids(&read_ids, &template).unwrap();
        let decoded = decode_read_ids(&encoded, &template, read_ids.len()).unwrap();

        assert_eq!(read_ids, decoded);
    }

    #[test]
    fn test_with_comments() {
        let read_ids = vec![
            "@INSTRUMENT:123:FLOWCELL:1:2101:10000:1000 1:N:0:ATCG".to_string(),
            "@INSTRUMENT:123:FLOWCELL:1:2101:10001:1001 1:N:0:ATCG".to_string(),
        ];

        let template = analyze_read_ids(&read_ids).unwrap();
        assert!(template.has_comment);

        let encoded = encode_read_ids(&read_ids, &template).unwrap();
        let decoded = decode_read_ids(&encoded, &template, read_ids.len()).unwrap();

        assert_eq!(read_ids, decoded);
    }

    #[test]
    fn test_non_standard_format() {
        let read_ids = vec![
            "READ_1".to_string(),
            "READ_2".to_string(),
        ];

        let template = analyze_read_ids(&read_ids).unwrap();
        assert!(template.prefix.is_empty()); // Falls back to no template
        assert!(template.common_comment.is_none());

        let encoded = encode_read_ids(&read_ids, &template).unwrap();
        let decoded = decode_read_ids(&encoded, &template, read_ids.len()).unwrap();

        assert_eq!(read_ids, decoded);
    }

    #[test]
    fn test_compression_ratio() {
        let read_ids = vec![
            "@INSTRUMENT:123:FLOWCELL:1:2101:10000:1000".to_string(),
            "@INSTRUMENT:123:FLOWCELL:1:2101:10001:1001".to_string(),
        ];

        let original_size: usize = read_ids.iter().map(|s| s.len()).sum();

        let template = analyze_read_ids(&read_ids).unwrap();
        let encoded = encode_read_ids(&read_ids, &template).unwrap();

        // Each encoded entry: 4 bytes (X) + 4 bytes (Y) = 8 bytes
        // vs original ~48 bytes per ID
        assert!(encoded.len() < original_size / 4); // At least 4x compression
    }

    #[test]
    fn test_small_delta_packing() {
        // Small deltas (all fit in 2 bits: -1, 0, 1)
        // Start with 0:0 so first delta is (0,0)
        let read_ids = vec![
            "@INSTRUMENT:123:FLOWCELL:1:2101:0:0".to_string(),
            "@INSTRUMENT:123:FLOWCELL:1:2101:1:0".to_string(),  // delta: 1, 0
            "@INSTRUMENT:123:FLOWCELL:1:2101:1:1".to_string(),  // delta: 0, 1
            "@INSTRUMENT:123:FLOWCELL:1:2101:0:1".to_string(),  // delta: -1, 0
            "@INSTRUMENT:123:FLOWCELL:1:2101:0:0".to_string(),  // delta: 0, -1
        ];

        let template = analyze_read_ids(&read_ids).unwrap();
        let encoded = encode_read_ids(&read_ids, &template).unwrap();
        let decoded = decode_read_ids(&encoded, &template, read_ids.len()).unwrap();

        assert_eq!(read_ids, decoded);

        // With small deltas, should use packed encoding (flag = 1)
        assert_eq!(encoded[0], 1); // Encoding flag
    }

    #[test]
    fn test_large_delta_varint() {
        // Large deltas (don't fit in 2 bits)
        let read_ids = vec![
            "@INSTRUMENT:123:FLOWCELL:1:2101:10000:1000".to_string(),
            "@INSTRUMENT:123:FLOWCELL:1:2101:10005:1000".to_string(), // delta: 5, 0 (doesn't fit)
            "@INSTRUMENT:123:FLOWCELL:1:2101:10006:1001".to_string(), // delta: 1, 1
        ];

        let template = analyze_read_ids(&read_ids).unwrap();
        let encoded = encode_read_ids(&read_ids, &template).unwrap();
        let decoded = decode_read_ids(&encoded, &template, read_ids.len()).unwrap();

        assert_eq!(read_ids, decoded);

        // With large deltas, should use varint encoding (flag = 0)
        assert_eq!(encoded[0], 0); // Encoding flag
    }

    #[test]
    fn test_pack_unpack_deltas() {
        let deltas = vec![
            (1, 0),
            (0, 1),
            (-1, 0),
            (0, -1),
            (1, 1),
            (-1, -1),
        ];

        let packed = pack_small_deltas(&deltas);
        let unpacked = unpack_small_deltas(&packed, deltas.len());

        assert_eq!(deltas, unpacked);
    }

    #[test]
    fn test_packed_encoding_size() {
        // Create read IDs with small deltas that should trigger packed encoding
        // Pattern: move in a small grid with deltas always in {-1, 0, 1}
        let mut read_ids = vec!["@INSTRUMENT:123:FLOWCELL:1:2101:0:0".to_string()];

        // Add 99 more reads with small deltas (only -1, 0, or 1)
        let mut x = 0i32;
        let mut y = 0i32;
        for i in 1..100 {
            // Alternate between moving right (1,0) and up (0,1)
            // with occasional diagonal (1,1) or stay (0,0)
            match i % 4 {
                0 => { x += 1; }       // Move right
                1 => { y += 1; }       // Move up
                2 => { x += 1; y += 1; } // Move diagonal
                _ => {}                // Stay in place
            }
            read_ids.push(format!("@INSTRUMENT:123:FLOWCELL:1:2101:{}:{}", x, y));
        }

        let template = analyze_read_ids(&read_ids).unwrap();
        let encoded_packed = encode_read_ids(&read_ids, &template).unwrap();

        // Should use packed encoding (flag = 1)
        assert_eq!(encoded_packed[0], 1);

        // Verify roundtrip
        let decoded = decode_read_ids(&encoded_packed, &template, read_ids.len()).unwrap();
        assert_eq!(read_ids, decoded);

        // Packed encoding should be more efficient than varint for clustered data
        // 100 reads = 50 bytes (2 pairs per byte) + 1 flag byte = 51 bytes for coordinates
        // vs varint which would be ~100-200 bytes
        assert!(encoded_packed.len() < 100); // Should be significantly smaller
    }
}
