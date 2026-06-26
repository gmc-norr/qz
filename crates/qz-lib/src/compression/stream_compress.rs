//! Block-granular streaming single-end compress engine (spec
//! 2026-06-15-compress-flat-stream-unify-design.md, Increment 1). Pure/compute
//! pieces: ArchivePlan, per-stream stagers, per-block compress jobs, and the
//! ordered chunk buffer. The orchestrator lives in compress_impl.rs.

use std::collections::{BTreeMap, HashMap};

use anyhow::Result;
use crate::cli::{CompressConfig, HeaderCompressor, QualityCompressor, QualityMode};
use crate::compression::chunk_directory::StreamRole;
use crate::compression::{EncodingType, QualityBinning};
use crate::compression::compress_impl::ChunkedMode;

/// A raw (uncompressed) BSC block ready for a worker. `offsets` are block-relative,
/// sentinel-terminated (offsets[0]==0, *offsets.last()==bytes.len()).
pub(crate) struct RawBscBlock {
    pub chunk_index: u32,
    pub mate: u8,
    pub role: StreamRole,
    pub block_idx: u32,
    pub bytes: Vec<u8>,
    pub offsets: Vec<usize>,
    /// Number of records in the raw block. The compress job recomputes the framed
    /// record_count from the splitter (single source of truth), so production never
    /// reads this — it's a descriptor the stagers/tests assert against.
    #[allow(dead_code)]
    pub record_count: u32,
}

/// Accumulates one role's framed record bytes and flushes a RawBscBlock at the
/// same end-based breakpoint compress_parallel_with_breakpoints uses (a record
/// joins the current block iff the block stays ≤ target after appending; the first
/// overshoot starts a new block; a single oversize record is emitted alone).
pub(crate) struct SeqQualHeaderStager {
    role: StreamRole,
    mate: u8,
    target: usize,
    bytes: Vec<u8>,
    offsets: Vec<usize>, // block-relative, offsets[0]==0
    next_block_idx: u32,
}

impl SeqQualHeaderStager {
    pub(crate) fn new(role: StreamRole, mate: u8, target: usize) -> Self {
        Self { role, mate, target, bytes: Vec::new(), offsets: vec![0], next_block_idx: 0 }
    }

    /// Append one record's already-framed bytes. Returns the completed block if
    /// appending this record would exceed `target` (the splitter's end-based
    /// partition: the record that would overshoot starts a new block instead).
    /// `chunk_index` is only consumed when a block is born here; it tags the
    /// emitted block and is otherwise ignored (the caller passes the same chunk's
    /// index every call, force-flushing at chunk boundaries so a block never spans
    /// chunks).
    pub(crate) fn stage_record_bytes(&mut self, chunk_index: u32, framed: &[u8]) -> Option<RawBscBlock> {
        let out = if !self.bytes.is_empty() && self.bytes.len() + framed.len() > self.target {
            self.take_block(chunk_index)
        } else { None };
        self.bytes.extend_from_slice(framed);
        self.offsets.push(self.bytes.len());
        out
    }

    /// Flush any buffered records as a final (possibly partial) block. Call at chunk
    /// boundary and at EOF.
    pub(crate) fn flush(&mut self, chunk_index: u32) -> Option<RawBscBlock> {
        if self.bytes.is_empty() { return None; }
        self.take_block(chunk_index)
    }

    fn take_block(&mut self, chunk_index: u32) -> Option<RawBscBlock> {
        if self.bytes.is_empty() { return None; }
        let bytes = std::mem::take(&mut self.bytes);
        let offsets = std::mem::replace(&mut self.offsets, vec![0]);
        let record_count = (offsets.len() - 1) as u32;
        let block_idx = self.next_block_idx;
        self.next_block_idx += 1;
        Some(RawBscBlock { chunk_index, mate: self.mate, role: self.role, block_idx, bytes, offsets, record_count })
    }

    /// Reset per-block numbering at a chunk boundary (after flush()).
    pub(crate) fn reset_for_next_chunk(&mut self) { self.next_block_idx = 0; }
}

#[cfg(test)]
mod stager_tests {
    use super::*; // brings in RawBscBlock, SeqQualHeaderStager, and StreamRole

    /// Build a raw byte stream + sentinel offsets the way the legacy path does for a
    /// role with no per-record varint (constant-length seq), then split with the
    /// breakpoint splitter; the stager fed the same per-record bytes must emit the
    /// SAME (record_count, bytes) blocks.
    fn check_matches_splitter(records: &[Vec<u8>], target: usize) {
        let mut data = Vec::new();
        let mut offs = vec![0usize];
        for r in records { data.extend_from_slice(r); offs.push(data.len()); }
        let legacy = bsc_blocks_unframed(&data, &offs, target);

        let mut st = SeqQualHeaderStager::new(StreamRole::Sequence, 0, target);
        let mut got: Vec<(u32, Vec<u8>)> = Vec::new();
        for r in records {
            if let Some(b) = st.stage_record_bytes(0, r) { got.push((b.record_count, b.bytes)); }
        }
        if let Some(b) = st.flush(0) { got.push((b.record_count, b.bytes)); }
        assert_eq!(legacy, got);
    }

    /// Mirror of compress_parallel_with_breakpoints' block partition, returning raw
    /// (record_count, bytes) blocks (no compression) for boundary comparison.
    /// SYNC: must mirror `compress_parallel_with_breakpoints` in
    /// `crates/qz-lib/src/compression/bsc.rs`. If that splitter's partition loop
    /// changes, update this oracle to match or the stager test verifies nothing.
    fn bsc_blocks_unframed(data: &[u8], offs: &[usize], target: usize) -> Vec<(u32, Vec<u8>)> {
        let n = offs.len() - 1;
        let mut ends = vec![0usize]; let mut cur = 0usize;
        while cur < n {
            let tb = offs[cur] + target;
            let pi = offs[cur..=n].partition_point(|&o| o <= tb);
            let mut nxt = cur + pi.saturating_sub(1);
            if nxt <= cur { nxt = cur + 1; }
            nxt = nxt.min(n); ends.push(nxt); cur = nxt;
        }
        ends.windows(2).map(|w| ((w[1]-w[0]) as u32, data[offs[w[0]]..offs[w[1]]].to_vec())).collect()
    }

    #[test]
    fn stager_constant_length_matches_splitter() {
        let recs: Vec<Vec<u8>> = (0..50).map(|_| vec![b'A'; 10]).collect();
        check_matches_splitter(&recs, 32);
    }

    #[test]
    fn stager_variable_length_matches_splitter() {
        let recs: Vec<Vec<u8>> = (0..50).map(|i| vec![b'A'; 5 + (i % 13)]).collect();
        check_matches_splitter(&recs, 41);
    }

    #[test]
    fn stager_oversize_single_record_is_own_block() {
        let recs = vec![vec![b'A'; 5], vec![b'B'; 100], vec![b'C'; 5]];
        check_matches_splitter(&recs, 10);
    }

    /// Stronger than the oracle tests above: this drives the stager AND the REAL splitter
    /// (`bsc::compress_parallel_with_breakpoints`) and asserts identical block boundaries
    /// (record_count sequence) AND identical compressed payloads. Each staged block is
    /// recompressed alone with an effectively-unbounded target so it stays one block;
    /// since BSC is deterministic, equal payloads ⇒ equal block bytes ⇒ equal boundaries.
    /// This fails automatically if bsc's partition loop ever changes — so it does not rely
    /// on the hand-mirrored `bsc_blocks_unframed` oracle staying in sync.
    #[test]
    fn stager_matches_real_splitter_compressed() {
        use crate::compression::bsc::compress_parallel_with_breakpoints as split;
        let recs: Vec<Vec<u8>> = (0..60).map(|i| vec![b'A' + (i % 4) as u8; 5 + (i % 11)]).collect();
        let target = 37usize;
        const UNBOUNDED: usize = 1 << 30; // bigger than the whole input → one block per call

        // Legacy global buffer + sentinel offsets, partitioned by the real splitter.
        let mut data = Vec::new();
        let mut offs = vec![0usize];
        for r in &recs { data.extend_from_slice(r); offs.push(data.len()); }
        let legacy = split(&data, &offs, target, false).unwrap();

        // Stager blocks, each recompressed alone → exactly one block apiece.
        let mut st = SeqQualHeaderStager::new(StreamRole::Sequence, 0, target);
        let mut staged: Vec<RawBscBlock> = Vec::new();
        for r in &recs { if let Some(b) = st.stage_record_bytes(0, r) { staged.push(b); } }
        if let Some(b) = st.flush(0) { staged.push(b); }
        let got: Vec<(u32, Vec<u8>)> = staged.iter().map(|b| {
            let mut c = split(&b.bytes, &b.offsets, UNBOUNDED, false).unwrap();
            assert_eq!(c.len(), 1, "a staged block must recompress to exactly one block");
            (c[0].0, std::mem::take(&mut c[0].1))
        }).collect();

        assert_eq!(legacy, got, "stager boundaries+payloads must equal the real splitter");
    }
}

/// A fqz quality block: the quality strings for ONE fqz sub-block, in record order.
///
/// INVARIANT: the producer's `QualStager` pre-caps this to ≤ `quality_ctx_block_size`
/// (default 500K) records — i.e. each `RawFqzBlock` is exactly one of the slices the
/// legacy `fqz_record_capped_blocks` produces via `par_chunks(cap)`. That is why
/// `compress_job` compresses it with a single `compress_qualities_fqzcomp_quals` call
/// (the same primitive the legacy path applies per slice): one block == one legacy
/// sub-block ⇒ byte-identical framing. `compress_job` does NOT re-partition; the
/// stager owns the cap.
pub(crate) struct RawFqzBlock {
    pub chunk_index: u32,
    pub mate: u8,
    pub block_idx: u32,
    pub quals: Vec<Vec<u8>>,
}

/// Whole-chunk header set for the columnar codec. `compress_job` BSC-compresses the
/// columnar blob and cap-splits it only if the compressed payload exceeds the block
/// cap (`columnar_blocks_capped` — bisect-after-compress), so it can yield N blocks.
pub(crate) struct ColumnarChunkJob {
    pub chunk_index: u32,
    pub mate: u8,
    pub headers: Vec<Vec<u8>>,
}

pub(crate) enum CompressJob {
    Bsc(RawBscBlock),
    Fqz(RawFqzBlock),
    Columnar(ColumnarChunkJob),
    PairedHeaders(PairedHeaderJob),
}

/// Both mates' header sets for one chunk. A worker compresses R1 columnar, R2 columnar,
/// and the R2-vs-R1 delta, then picks the smaller R2 representation — exactly mirroring
/// the legacy paired encoder's per-chunk decision, so the bytes are identical.
pub(crate) struct PairedHeaderJob {
    pub chunk_index: u32,
    pub r1_ids: Vec<Vec<u8>>,
    pub r2_ids: Vec<Vec<u8>>,
    pub bsc_block_mb: usize,
}

/// One compressed BSC/fqz/columnar block result (role + per-role block_idx).
pub(crate) struct BlockResult {
    pub chunk_index: u32,
    pub mate: u8,
    pub role: StreamRole,
    pub block_idx: u32,
    pub record_count: u32,
    pub compressed: Vec<u8>,
}

/// How a manifest segment's blocks reach the buffer.
pub(crate) enum SegSource {
    /// `n` `BlockResult`s keyed `(mate, role)` — sequence, quality, single-end rc &
    /// BSC-headers. `n` is known at chunk close.
    Blocks(u32),
    /// One whole-segment delivery for this mate's header slot — columnar headers
    /// (single-end mate 0, paired R1) or the paired R2 chosen header. The DELIVERED
    /// role/codec override the spec's placeholder (R2 may be Headers or HeaderDelta,
    /// known only post-compress).
    HeaderSegment,
}

/// One segment of a chunk, in physical write order. Single-end and paired both build a
/// `Vec<SegmentSpec>`; the writer emits them per the active `CompressMode` (§6).
pub(crate) struct SegmentSpec {
    pub mate: u8,
    /// For `Blocks`: the role+codec to emit. For `HeaderSegment`: a placeholder
    /// overridden by the delivered segment.
    pub role: StreamRole,
    pub codec: u8,
    pub source: SegSource,
}

/// Per-chunk completion declaration emitted by a producer at chunk close.
pub(crate) struct ChunkManifest {
    pub chunk_index: u32,
    pub segments: Vec<SegmentSpec>,
    /// Per-mate decoded FASTQ output bytes, feeding the `ChunkDecodedSizes` global.
    /// Single-end fills `[bytes, 0]`; paired fills `[r1_bytes, r2_bytes]`.
    pub decoded_output_bytes: [u64; 2],
    pub records: u64,       // → footer num_reads (Σ over chunks); = pairs for paired
    pub total_bases: u64,   // → stats log line
    pub original_size: u64, // → stats line only
}

/// A fully-compressed segment. Serves double duty: (a) a worker→writer message payload
/// for whole header-segment deliveries (so it carries `chunk_index`), and (b) a drained
/// segment inside `CompletedChunk` (where `chunk_index` is redundant but harmless).
/// `blocks` are in block_idx order.
pub(crate) struct CompletedSegment {
    pub chunk_index: u32,
    pub mate: u8,
    pub role: StreamRole,
    pub codec: u8,
    pub blocks: Vec<(u32, Vec<u8>)>, // (record_count, compressed)
}

/// A fully-compressed chunk, drained in chunk order. `segments` are in manifest order.
pub(crate) struct CompletedChunk {
    pub chunk_index: u32,
    pub records: u64,
    pub segments: Vec<CompletedSegment>,
    /// Per-mate decoded FASTQ output bytes (feeds the `ChunkDecodedSizes` global that the
    /// NUMA direct-write decode path pre-sizes from). Single-end fills `[bytes, 0]` (mate
    /// index 0); paired fills `[r1_bytes, r2_bytes]`.
    pub decoded_output_bytes: [u64; 2],
}

/// A delivered whole header segment (role/codec known at delivery, not from the spec).
struct DeliveredSegment {
    role: StreamRole,
    codec: u8,
    blocks: Vec<(u32, Vec<u8>)>,
}

/// A chunk's blocks accumulate keyed by (mate, role) → block_idx → (record_count, bytes).
/// Header segments arrive whole, keyed by mate. The manifest and the blocks/segments
/// arrive on one unordered channel; either may come first. Completion is checked only
/// once the manifest is present.
#[derive(Default)]
struct PartialChunk {
    manifest: Option<ChunkManifest>,
    blocks: HashMap<(u8, StreamRole), BTreeMap<u32, (u32, Vec<u8>)>>,
    header_segs: HashMap<u8, DeliveredSegment>, // keyed by mate (one header slot per mate)
}

pub(crate) struct OrderedBlockBuffer {
    next_to_emit: u32,
    chunks: HashMap<u32, PartialChunk>,
}

impl OrderedBlockBuffer {
    pub(crate) fn new() -> Self { Self { next_to_emit: 0, chunks: HashMap::new() } }

    fn entry(&mut self, chunk_index: u32) -> &mut PartialChunk {
        self.chunks.entry(chunk_index).or_default()
    }

    pub(crate) fn set_manifest(&mut self, m: ChunkManifest) {
        let ci = m.chunk_index;
        self.entry(ci).manifest = Some(m);
    }

    /// A per-(mate, role) compressed block (BSC sequence/quality/rc/headers, fqz quality).
    /// Inserted even pre-manifest. A duplicate (chunk, mate, role, block_idx) is a hard
    /// error (a producer/worker bug must fail the compress cleanly, not silently overwrite
    /// and ship a corrupt archive) — the writer surfaces this `Err` as a first-error-wins
    /// abort, identical to a worker failure.
    pub(crate) fn add_block(&mut self, r: BlockResult) -> Result<()> {
        let (chunk_index, mate, role, block_idx) = (r.chunk_index, r.mate, r.role, r.block_idx);
        let pc = self.entry(chunk_index);
        let slot = pc.blocks.entry((mate, role)).or_default();
        if slot.insert(block_idx, (r.record_count, r.compressed)).is_some() {
            anyhow::bail!(
                "internal: duplicate block (chunk={chunk_index}, mate={mate}, role={role:?}, block_idx={block_idx})"
            );
        }
        Ok(())
    }

    /// A whole header segment for `(chunk, mate)`. `blocks` MUST be in block_idx (0,1,2…)
    /// order. The delivered `role`/`codec` are authoritative (R2 may be Headers or
    /// HeaderDelta). A second delivery for the same (chunk, mate) is a hard error (same
    /// rationale as `add_block`).
    pub(crate) fn add_header_segment(
        &mut self, chunk_index: u32, mate: u8, role: StreamRole, codec: u8, blocks: Vec<(u32, Vec<u8>)>,
    ) -> Result<()> {
        let pc = self.entry(chunk_index);
        if pc.header_segs.contains_key(&mate) {
            anyhow::bail!("internal: header segment delivered twice for chunk {chunk_index} mate {mate}");
        }
        pc.header_segs.insert(mate, DeliveredSegment { role, codec, blocks });
        Ok(())
    }

    fn is_complete(pc: &PartialChunk) -> bool {
        let Some(m) = pc.manifest.as_ref() else { return false; };
        m.segments.iter().all(|s| match s.source {
            // Require the block indices to be EXACTLY 0..n — count alone would let a
            // misindexed delivery (gap/out-of-range idx) look complete and drain a chunk
            // with the wrong blocks. BTreeMap keys iterate ascending, so `eq(0..n)` checks
            // presence + contiguity in one pass (and implies len == n).
            SegSource::Blocks(n) => {
                pc.blocks.get(&(s.mate, s.role)).map_or(n == 0, |b| b.keys().copied().eq(0..n))
            }
            SegSource::HeaderSegment => pc.header_segs.contains_key(&s.mate),
        })
    }

    /// Emit all now-contiguous completed chunks in chunk order. Segments come out in
    /// manifest order; within `Blocks` segments, BTreeMap → block_idx order.
    pub(crate) fn drain_ready(&mut self) -> Vec<CompletedChunk> {
        let mut out = Vec::new();
        while self.chunks.get(&self.next_to_emit).is_some_and(Self::is_complete) {
            let mut pc = self.chunks.remove(&self.next_to_emit)
                .expect("invariant: is_some_and above proved this chunk present");
            let m = pc.manifest.take().expect("invariant: is_complete ⇒ manifest present");
            let mut segments = Vec::with_capacity(m.segments.len());
            for s in &m.segments {
                let cs = match s.source {
                    SegSource::Blocks(_) => {
                        let blocks = pc.blocks.remove(&(s.mate, s.role))
                            .map(|b| b.into_values().collect())
                            .unwrap_or_default();
                        CompletedSegment { chunk_index: m.chunk_index, mate: s.mate, role: s.role, codec: s.codec, blocks }
                    }
                    SegSource::HeaderSegment => {
                        let d = pc.header_segs.remove(&s.mate)
                            .expect("invariant: is_complete ⇒ header segment present");
                        CompletedSegment { chunk_index: m.chunk_index, mate: s.mate, role: d.role, codec: d.codec, blocks: d.blocks }
                    }
                };
                segments.push(cs);
            }
            out.push(CompletedChunk {
                chunk_index: m.chunk_index,
                records: m.records,
                segments,
                decoded_output_bytes: m.decoded_output_bytes,
            });
            self.next_to_emit += 1;
        }
        out
    }

    /// Whether every chunk drained AND no buffered data is orphaned. The `chunks.is_empty()`
    /// half matters because the buffer accepts pre-manifest blocks for ANY chunk_index: a
    /// MISROUTED block (wrong chunk_index ≥ total_chunks) would create a PartialChunk that
    /// never completes and never drains, yet `next_to_emit` would still reach `total_chunks`
    /// — masking the loss. Requiring an empty map turns that into a clean writer-close error.
    pub(crate) fn all_emitted(&self, total_chunks: u32) -> bool {
        self.next_to_emit == total_chunks && self.chunks.is_empty()
    }
}

/// Per-record framing + codec decisions, resolved ONCE from the first chunk (the
/// "planning prelude"). Mirrors what legacy `compress_chunked` computes from chunk 0.
#[derive(Clone)]
pub(crate) struct ArchivePlan {
    pub use_fqzcomp: bool,
    pub use_columnar_headers: bool,
    pub bsc_static: bool,
    pub bsc_block_mb: usize,
    pub quality_mode: QualityMode,
    /// `None` whenever `use_fqzcomp` (fqzcomp takes raw ASCII, no binning) or
    /// `no_quality` (there is no quality stream to bin at all).
    pub stream_quality_binning: QualityBinning,
    /// True for both `no_quality` AND `use_fqzcomp` — quality is omitted from
    /// `records_to_streams`/`frame_record` (fqz compresses it separately).
    pub skip_quality_in_stream: bool,
    pub sequence_hints: bool,
    pub rc_canon: bool,
    pub no_quality: bool,
    pub fasta: bool,
    pub const_seq_len: usize,
    /// Zero whenever `use_fqzcomp` — const-qual framing only applies to the BSC path.
    pub const_qual_len: usize,
    /// Records per fqz sub-block. Legacy flag name (`--quality-ctx-block-size`); now
    /// governs the fqzcomp sub-block record cap (default 500K).
    pub quality_ctx_block_size: usize,
    /// The `EncodingType::to_u8()` value (kept as `u8` to mirror legacy and feed
    /// `write_v5_archive_header` directly).
    pub encoding_type: u8,
    pub quality_compressor_used: QualityCompressor,
    pub header_compressor: HeaderCompressor,
}

/// Build the full `ArchivePlan` from config + the 3 data-derived inputs that differ
/// between the whole-file resolver (`resolve_archive_plan`, from the first chunk) and
/// the NUMA override builder (`archive_plan_from_override`, from a pinned PlanOverride).
/// ALL config/mode-derived fields are computed HERE — the single source of truth, so
/// the two entry points can never diverge. `const_qual_len_detected` is the raw
/// detected quality length; it is zeroed iff fqzcomp (const-qual framing is BSC-only).
fn build_archive_plan(
    args: &CompressConfig,
    mode: ChunkedMode,
    use_fqzcomp: bool,
    const_seq_len: usize,
    const_qual_len_detected: usize,
) -> ArchivePlan {
    let is_ultra = matches!(mode, ChunkedMode::Ultra { .. });
    let bsc_block_mb = match mode {
        ChunkedMode::Ultra { bsc_block_mb, .. } => bsc_block_mb,
        ChunkedMode::Default => args.advanced.bsc_block_size_mb,
    };
    let no_quality = args.no_quality;
    let quality_mode = args.quality_mode;
    let sequence_hints = if is_ultra { false } else { args.advanced.sequence_hints };
    let rc_canon = if is_ultra { false } else { args.advanced.rc_canon };
    let quality_binning = if no_quality { QualityBinning::None }
        else { super::quality_mode_to_binning(quality_mode) };

    let quality_compressor_used = if use_fqzcomp { QualityCompressor::Fqzcomp } else { QualityCompressor::Bsc };
    let skip_quality_in_stream = no_quality || use_fqzcomp;
    let stream_quality_binning = if use_fqzcomp { QualityBinning::None } else { quality_binning };
    let const_qual_len = if use_fqzcomp { 0 } else { const_qual_len_detected };
    let encoding_type = if is_ultra { EncodingType::UltraBigBlock }
        else if rc_canon { EncodingType::RcCanon }
        else if sequence_hints { EncodingType::RawWithHints }
        else { EncodingType::Raw }.to_u8();
    let quality_ctx_block_size = match mode {
        ChunkedMode::Ultra { quality_sub_block, .. } => quality_sub_block,
        ChunkedMode::Default => args.advanced.quality_ctx_block_size,
    };

    ArchivePlan {
        use_fqzcomp,
        use_columnar_headers: args.advanced.header_compressor == HeaderCompressor::Columnar,
        bsc_static: args.advanced.bsc_static,
        bsc_block_mb,
        quality_mode,
        stream_quality_binning,
        skip_quality_in_stream,
        sequence_hints,
        rc_canon,
        no_quality,
        fasta: args.fasta,
        const_seq_len,
        const_qual_len,
        quality_ctx_block_size,
        encoding_type,
        quality_compressor_used,
        header_compressor: args.advanced.header_compressor,
    }
}

pub(super) fn resolve_archive_plan(
    args: &CompressConfig,
    mode: ChunkedMode,
    first_chunk: &[crate::io::FastqRecord],
) -> Result<ArchivePlan> {
    let resolved = super::resolve_quality_compressor(
        args.advanced.quality_compressor, first_chunk.len(), args.quality_mode, args.no_quality, false);
    let use_fqzcomp = resolved == QualityCompressor::Fqzcomp;
    let (s, q) = super::detect_const_lengths(first_chunk, args.no_quality);
    Ok(build_archive_plan(args, mode, use_fqzcomp, s, q))
}

/// Build an `ArchivePlan` from config + a parent-resolved `PlanOverride` instead of
/// detecting the size/data-dependent fields from this worker's own first chunk
/// (spec §5.3). Config-derived fields are computed identically to
/// `resolve_archive_plan`; only the quality codec + const lengths come from the
/// override. Single-end (uses `quality_mate1`).
///
/// `pub(super)` (== `pub(in crate::compression)`) to match `ChunkedMode`'s
/// visibility — the only caller (`range_compress::compress_byte_range`) lives in
/// `crate::compression`, so this is fully reachable while staying warning-free.
pub(super) fn archive_plan_from_override(
    args: &CompressConfig,
    mode: ChunkedMode,
    ov: &crate::compression::PlanOverride,
) -> ArchivePlan {
    let use_fqzcomp = ov.quality_mate1 == QualityCompressor::Fqzcomp;
    build_archive_plan(args, mode, use_fqzcomp, ov.const_seq_len, ov.const_qual_len)
}

/// Build a per-mate `ArchivePlan` for paired-end compress: columnar headers, NEVER
/// const-length (paired always uses variable framing), no sequence_hints / rc_canon,
/// `Raw` encoding, `QualityBinning::None`, 25 MB blocks, adaptive (`bsc_static=false`).
/// `resolved_quality` comes from the mate's first-chunk read count
/// (`resolve_quality_compressor`).
pub(crate) fn paired_mate_plan(args: &CompressConfig, resolved_quality: QualityCompressor) -> ArchivePlan {
    let use_fqzcomp = resolved_quality == QualityCompressor::Fqzcomp;
    ArchivePlan {
        use_fqzcomp,
        use_columnar_headers: true, // paired headers are ALWAYS columnar (BSC headers rejected)
        bsc_static: false,
        bsc_block_mb: args.advanced.bsc_block_size_mb,
        quality_mode: args.quality_mode,
        stream_quality_binning: QualityBinning::None,
        skip_quality_in_stream: args.no_quality || use_fqzcomp,
        sequence_hints: false,
        rc_canon: false,
        no_quality: args.no_quality,
        fasta: false,
        const_seq_len: 0,
        const_qual_len: 0,
        quality_ctx_block_size: args.advanced.quality_ctx_block_size,
        encoding_type: EncodingType::Raw.to_u8(),
        quality_compressor_used: if use_fqzcomp { QualityCompressor::Fqzcomp } else { QualityCompressor::Bsc },
        header_compressor: HeaderCompressor::Columnar, // paired is always columnar
    }
}

/// Compress worker-pool size. The pool **follows the thread setting** (`-t`) — the
/// user's parallelism / RAM dial — because peak RSS scales with the pool's working set
/// (≈ workers × block × ~7 libsais workspace) and is constant in read count. A memory
/// safety budget caps the count so a huge `-t` on a RAM-constrained box (or a large
/// configured block size) can't OOM; peak RSS grows sublinearly as completed blocks free
/// their workspace, so for the default ~25 MiB block `threads` binds up to ~45 (the
/// single-socket bandwidth plateau — more workers buy no speed).
///
/// **Ultra** is special: its big blocks (188–750 MiB) parallelise INTERNALLY via libbsc
/// OpenMP MT-BWT, and running concurrent forward-BWTs on big blocks blows up libsais
/// memory (see `finding_ultra_mt_bwt_big_blocks`). So ultra always runs ONE block at a
/// time and lets the internal MT use the cores — regardless of `-t`.
///
/// The orchestrator's `QZ_COMPRESS_WORKERS` env var overrides the result entirely.
pub(crate) fn compress_worker_count(threads: usize, block_bytes: usize, is_ultra: bool) -> usize {
    if is_ultra {
        return 1;
    }
    const SAFETY_BUDGET: usize = 8 * 1024 * 1024 * 1024; // generous; only binds at huge -t / big blocks
    const WORKSPACE_FACTOR: usize = 7;
    let per_block = block_bytes.saturating_mul(WORKSPACE_FACTOR).max(1);
    let budget_cap = (SAFETY_BUDGET / per_block).max(1);
    threads.clamp(1, budget_cap)
}

/// Shared test fixture: a default-config plan. NOTE the empty `first_chunk` ⇒
/// `detect_const_lengths` returns (0,0) AND the 0 record count is below the fqzcomp
/// threshold, so this ALWAYS resolves to BSC quality with zero const-lengths. Tests
/// that need the fqz path or const-lengths must call `resolve_archive_plan` with a
/// real first chunk, not this helper.
#[cfg(test)]
pub(crate) fn tests_plan_default() -> ArchivePlan {
    resolve_archive_plan(&CompressConfig::default(), ChunkedMode::Default, &[]).unwrap()
}

/// Quality stager. fqz mode buffers ≤cap owned quality strings and flushes a fqz
/// job every `cap` records (chunk-relative; reset at chunk boundary) — matching
/// fqz_record_capped_blocks' par_chunks(cap). BSC mode delegates to a
/// SeqQualHeaderStager over framed quality bytes.
pub(crate) enum QualStager {
    Fqz { cap: usize, mate: u8, buf: Vec<Vec<u8>>, next_block_idx: u32 },
    Bsc(SeqQualHeaderStager),
    None, // no_quality
}

impl QualStager {
    pub(crate) fn fqz(cap: usize, mate: u8) -> Self {
        QualStager::Fqz { cap: cap.max(1), mate, buf: Vec::new(), next_block_idx: 0 }
    }
    pub(crate) fn bsc(mate: u8, target: usize) -> Self {
        QualStager::Bsc(SeqQualHeaderStager::new(StreamRole::Qual, mate, target))
    }

    /// fqz path: stage one record's quality (caller guarantees Some). Returns a job
    /// when the slice fills. Must only be called on the `Fqz` variant — a wrong-variant
    /// call would silently drop the record, so it is debug-asserted.
    pub(crate) fn stage_fqz(&mut self, chunk_index: u32, qual: &[u8]) -> Option<RawFqzBlock> {
        debug_assert!(
            matches!(self, QualStager::Fqz { .. }),
            "stage_fqz called on a non-Fqz QualStager — the producer must route quality \
             by variant (Fqz→stage_fqz, Bsc→bsc_stager, None→skip)"
        );
        if let QualStager::Fqz { cap, mate, buf, next_block_idx } = self {
            buf.push(qual.to_vec());
            if buf.len() == *cap {
                let quals = std::mem::take(buf);
                let block_idx = *next_block_idx; *next_block_idx += 1;
                return Some(RawFqzBlock { chunk_index, mate: *mate, block_idx, quals });
            }
        }
        Option::None
    }

    /// fqz path: flush the final (partial) slice at chunk boundary / EOF. Like
    /// `stage_fqz`, only valid on the `Fqz` variant.
    pub(crate) fn flush_fqz(&mut self, chunk_index: u32) -> Option<RawFqzBlock> {
        debug_assert!(
            matches!(self, QualStager::Fqz { .. }),
            "flush_fqz called on a non-Fqz QualStager"
        );
        if let QualStager::Fqz { mate, buf, next_block_idx, .. } = self && !buf.is_empty() {
            let quals = std::mem::take(buf);
            let block_idx = *next_block_idx; *next_block_idx += 1;
            return Some(RawFqzBlock { chunk_index, mate: *mate, block_idx, quals });
        }
        Option::None
    }

    /// BSC path: the inner stager the producer routes framed quality bytes through.
    /// `None` for the fqz / no-quality variants.
    pub(crate) fn bsc_stager(&mut self) -> Option<&mut SeqQualHeaderStager> {
        match self {
            QualStager::Bsc(s) => Some(s),
            _ => Option::None,
        }
    }

    pub(crate) fn reset_for_next_chunk(&mut self) {
        match self {
            QualStager::Fqz { next_block_idx, .. } => *next_block_idx = 0,
            QualStager::Bsc(s) => s.reset_for_next_chunk(),
            QualStager::None => {}
        }
    }
}

/// Result of one compress job: either per-(mate,role) blocks, or whole header segments.
pub(crate) enum JobOutput {
    Blocks(Vec<BlockResult>),
    Segments(Vec<CompletedSegment>),
}

/// Compress one job. BSC/fqz → one block (`Blocks`); the columnar / paired-headers jobs
/// deliver whole header segments (`Segments`). The paired-headers arm picks the smaller
/// R2 representation (independent columnar vs delta-vs-R1), tie → columnar.
pub(crate) fn compress_job(job: CompressJob, plan: &ArchivePlan) -> Result<JobOutput> {
    use crate::compression::chunk_directory::StreamRole;
    match job {
        CompressJob::Bsc(b) => {
            let target = plan.bsc_block_mb * 1024 * 1024;
            let mut blocks = crate::compression::bsc::compress_parallel_with_breakpoints(
                &b.bytes, &b.offsets, target, plan.bsc_static)?;
            debug_assert!(blocks.len() <= 1, "staged BSC block must yield ≤1 compressed block");
            let mut out = Vec::with_capacity(1);
            for (rc, mut comp) in blocks.drain(..) {
                comp.shrink_to_fit();
                out.push(BlockResult { chunk_index: b.chunk_index, mate: b.mate, role: b.role,
                    block_idx: b.block_idx, record_count: rc, compressed: comp });
            }
            Ok(JobOutput::Blocks(out))
        }
        CompressJob::Fqz(f) => {
            let quals: Vec<&[u8]> = f.quals.iter().map(|q| q.as_slice()).collect();
            let blob = super::codecs::compress_qualities_fqzcomp_quals(&quals)?;
            Ok(JobOutput::Blocks(vec![BlockResult { chunk_index: f.chunk_index, mate: f.mate,
                role: StreamRole::Qual, block_idx: f.block_idx,
                record_count: f.quals.len() as u32, compressed: blob }]))
        }
        CompressJob::Columnar(c) => {
            let refs: Vec<&[u8]> = c.headers.iter().map(|h| h.as_slice()).collect();
            let blocks = crate::compression::compress_impl::columnar_blocks_capped(
                &refs, crate::compression::paired::streams::MAX_BLOCK)?;
            Ok(JobOutput::Segments(vec![CompletedSegment {
                chunk_index: c.chunk_index, mate: c.mate, role: StreamRole::Headers,
                codec: crate::compression::codec_ids::CODEC_COLUMNAR, blocks }]))
        }
        CompressJob::PairedHeaders(j) => {
            use crate::compression::codec_ids::{CODEC_BSC, CODEC_COLUMNAR};
            use crate::compression::paired::{header_delta, streams};
            let cap = streams::MAX_BLOCK;
            let r1_refs: Vec<&[u8]> = j.r1_ids.iter().map(|h| h.as_slice()).collect();
            let r2_refs: Vec<&[u8]> = j.r2_ids.iter().map(|h| h.as_slice()).collect();
            // mate 1: independent columnar headers (always emitted). Serialize once — its
            // length is also a free, accurate proxy for what R2's independent-columnar
            // candidate would cost (mate read-name columns are near-identical in size).
            let r1_blocks = crate::compression::compress_impl::columnar_blocks_capped(&r1_refs, cap)?;
            let mut r1_bytes = Vec::new();
            streams::write_block_stream(&mut r1_bytes, &r1_blocks);

            // mate 2 candidate B (CHEAP — compute FIRST): delta of R2 ids vs R1 ids,
            // BSC-compressed (bsc_static=false, matching legacy compress_paired_v5_inner;
            // target = R1 plan's block size). On typical paired Illumina data R2 ids are R1
            // ids with the mate suffix flipped (…/1 → …/2), so this delta is a few hundred
            // bytes and always wins, while candidate A (independent columnar) is a full BWT
            // over millions of headers (~10⁴× larger output, ~80× slower) thrown away every
            // chunk. Measured HG002 chr20 10.8M pairs: candidate A = ~28 s of wasted CPU.
            let (delta_ops, delta_offs) = header_delta::encode(&r1_refs, &r2_refs);
            let target = j.bsc_block_mb * 1024 * 1024;
            let delta_blocks = crate::compression::bsc::compress_parallel_with_breakpoints(
                &delta_ops, &delta_offs, target, false)?;
            let mut delta_bytes = Vec::new();
            streams::write_block_stream(&mut delta_bytes, &delta_blocks);

            // When the delta is comfortably below the R1-columnar proxy (≥2× headroom) it
            // provably beats candidate A (whose size ≈ R1's), so SKIP the expensive R2
            // columnar entirely. Otherwise fall back to the exact compute-both-and-compare
            // (tie → independent columnar, `<=`) — byte-identical to the unconditional encoder.
            let r2 = if delta_bytes.len().saturating_mul(2) <= r1_bytes.len() {
                CompletedSegment { chunk_index: j.chunk_index, mate: 2, role: StreamRole::HeaderDelta, codec: CODEC_BSC, blocks: delta_blocks }
            } else {
                let indep_blocks = crate::compression::compress_impl::columnar_blocks_capped(&r2_refs, cap)?;
                let mut indep_bytes = Vec::new();
                streams::write_block_stream(&mut indep_bytes, &indep_blocks);
                if indep_bytes.len() <= delta_bytes.len() {
                    CompletedSegment { chunk_index: j.chunk_index, mate: 2, role: StreamRole::Headers, codec: CODEC_COLUMNAR, blocks: indep_blocks }
                } else {
                    CompletedSegment { chunk_index: j.chunk_index, mate: 2, role: StreamRole::HeaderDelta, codec: CODEC_BSC, blocks: delta_blocks }
                }
            };
            Ok(JobOutput::Segments(vec![
                CompletedSegment { chunk_index: j.chunk_index, mate: 1, role: StreamRole::Headers, codec: CODEC_COLUMNAR, blocks: r1_blocks },
                r2,
            ]))
        }
    }
}

#[cfg(test)]
mod job_tests {
    use super::*;
    use crate::compression::{bsc, chunk_directory::StreamRole};

    #[test]
    fn bsc_block_compresses_byte_identical() {
        let data: Vec<u8> = (0..2000u32).flat_map(|i| (i as u8).to_le_bytes()).collect();
        let offs: Vec<usize> = (0..=data.len()).step_by(8).collect(); // 8-byte records
        let target = 1 << 20;
        let want = bsc::compress_parallel_with_breakpoints(&data, &offs, target, false).unwrap();
        assert_eq!(want.len(), 1, "fits one block");
        let raw = RawBscBlock {
            chunk_index: 0, mate: 0, role: StreamRole::Sequence, block_idx: 0,
            bytes: data.clone(), offsets: offs.clone(),
            record_count: (offs.len()-1) as u32,
        };
        let JobOutput::Blocks(results) = compress_job(CompressJob::Bsc(raw), &plan_bsc(target)).unwrap()
            else { panic!("expected blocks") };
        assert_eq!(results.len(), 1);
        let mut w0 = want[0].1.clone(); w0.shrink_to_fit();
        assert_eq!((results[0].record_count, &results[0].compressed), (want[0].0, &w0));
    }

    /// fqz arm wiring: one RawFqzBlock → one BlockResult, role=Qual, record_count==n,
    /// chunk/block indices threaded through. (Byte-identity vs the legacy fqz slice is
    /// proven in the QualStager task; this pins the dispatch.)
    #[test]
    fn fqz_arm_emits_one_qual_block() {
        let quals: Vec<Vec<u8>> = (0..5u8).map(|i| vec![b'I', b'5', b'#', b'0' + i]).collect();
        let job = CompressJob::Fqz(RawFqzBlock { chunk_index: 3, mate: 0, block_idx: 2, quals: quals.clone() });
        let JobOutput::Blocks(results) = compress_job(job, &plan_bsc(1 << 20)).unwrap()
            else { panic!("expected blocks") };
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!((r.chunk_index, r.block_idx, r.role), (3, 2, StreamRole::Qual));
        assert_eq!(r.record_count, 5);
        assert!(!r.compressed.is_empty());
    }

    /// columnar arm wiring: a small header set fits under the cap → one Headers segment
    /// with one block (record_count == header count), and chunk_index threaded through.
    #[test]
    fn columnar_arm_emits_headers_block() {
        let headers: Vec<Vec<u8>> =
            (0..4).map(|i| format!("@SRR{i:08}.{i}/1").into_bytes()).collect();
        let job = CompressJob::Columnar(ColumnarChunkJob { chunk_index: 7, mate: 0, headers });
        let JobOutput::Segments(segs) = compress_job(job, &plan_bsc(1 << 20)).unwrap()
            else { panic!("expected segments") };
        assert_eq!(segs.len(), 1);
        let s = &segs[0];
        assert_eq!((s.chunk_index, s.mate, s.role), (7, 0, StreamRole::Headers));
        assert_eq!(s.blocks.len(), 1, "4 short headers fit under the block cap");
        assert_eq!(s.blocks[0].0, 4);
    }

    /// Pins the `PairedHeaders` smaller-of DECISION directly (independent of golden
    /// capture), so the rarely-hit independent-columnar branch can't drift silently.
    #[test]
    fn paired_headers_picks_delta_then_columnar() {
        use crate::compression::chunk_directory::StreamRole;
        // R2 differs from R1 by one byte (the mate flag) → tiny delta → HeaderDelta wins.
        let r1: Vec<Vec<u8>> = (0..50).map(|i| format!("@SRR{i:09}.{i} 1:N:0:ACGT").into_bytes()).collect();
        let r2: Vec<Vec<u8>> = (0..50).map(|i| format!("@SRR{i:09}.{i} 2:N:0:ACGT").into_bytes()).collect();
        let JobOutput::Segments(segs) = compress_job(
            CompressJob::PairedHeaders(PairedHeaderJob { chunk_index: 0, r1_ids: r1, r2_ids: r2, bsc_block_mb: 25 }),
            &tests_plan_default()).unwrap() else { panic!("expected segments") };
        assert_eq!(segs.len(), 2);
        assert_eq!((segs[0].mate, segs[1].mate), (1, 2));
        assert_eq!(segs[1].role, StreamRole::HeaderDelta, "delta-friendly ids → delta branch");

        // R2 unrelated to R1 → big literal delta → independent columnar wins.
        let r1u: Vec<Vec<u8>> = (0..50).map(|i| format!("@SRR{i:09}.{i} 1:N:0:ACGT").into_bytes()).collect();
        let r2u: Vec<Vec<u8>> = (0..50).map(|i| format!("@LANE7:TILE{:05}:X{}:Y{}", 50 - i, i * 31, i * 17).into_bytes()).collect();
        let JobOutput::Segments(segs) = compress_job(
            CompressJob::PairedHeaders(PairedHeaderJob { chunk_index: 0, r1_ids: r1u, r2_ids: r2u, bsc_block_mb: 25 }),
            &tests_plan_default()).unwrap() else { panic!("expected segments") };
        assert_eq!(segs[1].role, StreamRole::Headers, "unrelated ids → independent columnar branch");
    }

    fn plan_bsc(target: usize) -> ArchivePlan {
        ArchivePlan { bsc_block_mb: target >> 20, ..tests_plan_default() }
    }
}

/// Buffers a chunk's header byte-copies for the columnar codec (compress-then-bisect
/// is post-compression, so headers are chunk-bound — ~75 MB/chunk). Records aren't
/// held by the streaming producer, so it copies each id in.
pub(crate) struct HeaderColStager {
    mate: u8,
    headers: Vec<Vec<u8>>,
}

impl HeaderColStager {
    pub(crate) fn new(mate: u8) -> Self {
        Self { mate, headers: Vec::new() }
    }
    /// Stage one header id (read name) for the current chunk, copied in.
    pub(crate) fn stage(&mut self, id: &[u8]) {
        self.headers.push(id.to_vec());
    }
    /// Take the chunk's buffered headers as a job (leaving the stager empty for the
    /// next chunk). Returns None if no headers were staged.
    pub(crate) fn take_chunk_job(&mut self, chunk_index: u32) -> Option<ColumnarChunkJob> {
        if self.headers.is_empty() {
            return None;
        }
        Some(ColumnarChunkJob { chunk_index, mate: self.mate, headers: std::mem::take(&mut self.headers) })
    }
}

#[cfg(test)]
mod header_tests {
    use super::*;
    use crate::compression::compress_impl::columnar_blocks_capped;
    use crate::compression::paired::streams::MAX_BLOCK;

    #[test]
    fn columnar_chunk_job_matches_legacy() {
        let hdrs: Vec<Vec<u8>> = (0..40).map(|i| format!("@SRR{:08}.{}/1", i, i * 7).into_bytes()).collect();
        let refs: Vec<&[u8]> = hdrs.iter().map(|h| h.as_slice()).collect();
        let want = columnar_blocks_capped(&refs, MAX_BLOCK).unwrap();

        let mut st = HeaderColStager::new(0);
        for h in &hdrs { st.stage(h); }
        let job = st.take_chunk_job(7).expect("non-empty");
        assert_eq!(job.chunk_index, 7);
        assert_eq!(job.headers, hdrs, "stager buffers the chunk's headers in order");

        let plan = tests_plan_default(); // columnar arm ignores BSC plan fields
        let JobOutput::Segments(segs) = compress_job(CompressJob::Columnar(job), &plan).unwrap()
            else { panic!("expected segments") };
        assert_eq!(segs.len(), 1);
        let got: Vec<(u32, Vec<u8>)> = segs.into_iter().next().unwrap().blocks;
        assert_eq!(want, got, "columnar chunk job + compress == legacy columnar_blocks_capped");
    }

    #[test]
    fn empty_stager_yields_no_job() {
        let mut st = HeaderColStager::new(0);
        assert!(st.take_chunk_job(0).is_none());
    }
}

#[cfg(test)]
mod order_tests {
    use super::*;
    use crate::compression::chunk_directory::StreamRole;

    fn br(ci: u32, mate: u8, role: StreamRole, bi: u32, rc: u32) -> BlockResult {
        BlockResult { chunk_index: ci, mate, role, block_idx: bi, record_count: rc,
            compressed: vec![ci as u8, mate, bi as u8] }
    }
    fn seg_blocks(mate: u8, role: StreamRole, codec: u8, n: u32) -> SegmentSpec {
        SegmentSpec { mate, role, codec, source: SegSource::Blocks(n) }
    }
    fn seg_header(mate: u8) -> SegmentSpec {
        SegmentSpec { mate, role: StreamRole::Headers, codec: 0, source: SegSource::HeaderSegment }
    }
    fn manifest(ci: u32, segments: Vec<SegmentSpec>, records: u64, dob: u64) -> ChunkManifest {
        ChunkManifest { chunk_index: ci, segments, decoded_output_bytes: [dob, 0], records, total_bases: 0, original_size: 0 }
    }

    #[test]
    fn drains_in_chunk_order_under_out_of_order_arrival() {
        let mut buf = OrderedBlockBuffer::new();
        buf.set_manifest(manifest(0, vec![
            seg_header(0),
            seg_blocks(0, StreamRole::Sequence, 1, 2),
            seg_blocks(0, StreamRole::Qual, 6, 1),
        ], 9, 11));
        buf.set_manifest(manifest(1, vec![
            seg_header(0),
            seg_blocks(0, StreamRole::Sequence, 1, 1),
            seg_blocks(0, StreamRole::Qual, 6, 1),
        ], 5, 22));
        buf.add_block(br(1, 0, StreamRole::Sequence, 0, 5)).unwrap();
        buf.add_block(br(0, 0, StreamRole::Sequence, 1, 5)).unwrap();
        buf.add_block(br(0, 0, StreamRole::Qual, 0, 9)).unwrap();
        assert!(buf.drain_ready().is_empty());
        buf.add_header_segment(0, 0, StreamRole::Headers, 3, vec![(9, vec![1])]).unwrap();
        buf.add_block(br(0, 0, StreamRole::Sequence, 0, 4)).unwrap();
        let d = buf.drain_ready();
        assert_eq!(d.len(), 1);
        assert_eq!(d[0].chunk_index, 0);
        // segments in manifest order: header, sequence(2 blocks), qual(1)
        assert_eq!(d[0].segments[0].role, StreamRole::Headers);
        assert_eq!(d[0].segments[1].blocks.iter().map(|(_, b)| b[2]).collect::<Vec<_>>(), vec![0, 1]);
        assert!(buf.drain_ready().is_empty());
        buf.add_header_segment(1, 0, StreamRole::Headers, 3, vec![(5, vec![1])]).unwrap();
        buf.add_block(br(1, 0, StreamRole::Qual, 0, 5)).unwrap();
        let d2 = buf.drain_ready();
        assert_eq!(d2.len(), 1);
        assert_eq!((d2[0].chunk_index, d2[0].decoded_output_bytes[0]), (1, 22));
    }

    #[test]
    fn header_segment_completes_chunk() {
        let mut buf = OrderedBlockBuffer::new();
        buf.set_manifest(manifest(0, vec![
            seg_header(0),
            seg_blocks(0, StreamRole::Sequence, 1, 1),
        ], 1, 3));
        buf.add_block(br(0, 0, StreamRole::Sequence, 0, 1)).unwrap();
        assert!(buf.drain_ready().is_empty(), "header still pending");
        buf.add_header_segment(0, 0, StreamRole::Headers, 3, vec![(2, vec![9]), (1, vec![8])]).unwrap();
        let d = buf.drain_ready();
        assert_eq!(d.len(), 1);
        assert_eq!(d[0].segments[0].blocks.len(), 2);
    }

    #[test]
    fn blocks_arrive_before_manifest() {
        let mut buf = OrderedBlockBuffer::new();
        buf.add_block(br(0, 0, StreamRole::Sequence, 0, 4)).unwrap();
        buf.add_header_segment(0, 0, StreamRole::Headers, 3, vec![(4, vec![1])]).unwrap();
        assert!(buf.drain_ready().is_empty(), "no manifest yet");
        buf.set_manifest(manifest(0, vec![
            seg_header(0),
            seg_blocks(0, StreamRole::Sequence, 1, 1),
        ], 4, 7));
        let d = buf.drain_ready();
        assert_eq!(d.len(), 1);
        assert_eq!((d[0].segments[0].blocks.len(), d[0].segments[1].blocks.len()), (1, 1));
    }

    /// Two mates, R2 header delivered as HeaderDelta — the delivered role overrides the
    /// spec placeholder.
    #[test]
    fn paired_two_mate_chunk_completes_and_orders() {
        let mut buf = OrderedBlockBuffer::new();
        buf.set_manifest(manifest(0, vec![
            seg_header(1),
            seg_blocks(1, StreamRole::Sequence, 1, 1),
            seg_blocks(1, StreamRole::Qual, 6, 1),
            seg_header(2),
            seg_blocks(2, StreamRole::Sequence, 1, 1),
            seg_blocks(2, StreamRole::Qual, 6, 1),
        ], 3, 0));
        buf.add_block(br(0, 1, StreamRole::Sequence, 0, 3)).unwrap();
        buf.add_block(br(0, 1, StreamRole::Qual, 0, 3)).unwrap();
        buf.add_block(br(0, 2, StreamRole::Sequence, 0, 3)).unwrap();
        buf.add_block(br(0, 2, StreamRole::Qual, 0, 3)).unwrap();
        buf.add_header_segment(0, 1, StreamRole::Headers, 3, vec![(3, vec![1])]).unwrap();
        assert!(buf.drain_ready().is_empty(), "mate-2 header still pending");
        buf.add_header_segment(0, 2, StreamRole::HeaderDelta, 1, vec![(3, vec![2])]).unwrap();
        let d = buf.drain_ready();
        assert_eq!(d.len(), 1);
        let roles: Vec<_> = d[0].segments.iter().map(|s| (s.mate, s.role)).collect();
        assert_eq!(roles, vec![
            (1, StreamRole::Headers), (1, StreamRole::Sequence), (1, StreamRole::Qual),
            (2, StreamRole::HeaderDelta), (2, StreamRole::Sequence), (2, StreamRole::Qual),
        ]);
    }

    #[test]
    fn fresh_buffer_drains_empty_and_zero_chunks_all_emitted() {
        let mut buf = OrderedBlockBuffer::new();
        assert!(buf.drain_ready().is_empty());
        assert!(buf.all_emitted(0));
        assert!(!buf.all_emitted(1));
    }

    #[test]
    fn duplicate_block_is_hard_error() {
        let mut buf = OrderedBlockBuffer::new();
        buf.add_block(br(0, 0, StreamRole::Sequence, 0, 5)).unwrap();
        // Same (chunk, mate, role, block_idx) again → hard error in release too, not a
        // silent overwrite that would ship a corrupt archive.
        assert!(buf.add_block(br(0, 0, StreamRole::Sequence, 0, 5)).is_err());
    }

    #[test]
    fn duplicate_header_segment_is_hard_error() {
        let mut buf = OrderedBlockBuffer::new();
        buf.add_header_segment(0, 1, StreamRole::Headers, 3, vec![(3, vec![1])]).unwrap();
        assert!(
            buf.add_header_segment(0, 1, StreamRole::Headers, 3, vec![(3, vec![2])]).is_err()
        );
    }

    #[test]
    fn noncontiguous_block_indices_never_complete() {
        // Two blocks present but at idx {0, 2} (gap at 1). Count == n (2) would wrongly
        // pass; the contiguity check keeps the chunk incomplete so all_emitted fails
        // cleanly instead of draining the wrong blocks.
        let mut buf = OrderedBlockBuffer::new();
        buf.set_manifest(manifest(0, vec![seg_blocks(0, StreamRole::Sequence, 1, 2)], 5, 0));
        buf.add_block(br(0, 0, StreamRole::Sequence, 0, 3)).unwrap();
        buf.add_block(br(0, 0, StreamRole::Sequence, 2, 2)).unwrap(); // idx 2, not 1
        assert!(buf.drain_ready().is_empty(), "gap at idx 1 → not complete");
        assert!(!buf.all_emitted(1));
    }

    #[test]
    fn contiguous_block_indices_complete() {
        // The same 2-block segment with the correct idx {0, 1} drains as one chunk.
        let mut buf = OrderedBlockBuffer::new();
        buf.set_manifest(manifest(0, vec![seg_blocks(0, StreamRole::Sequence, 1, 2)], 5, 0));
        buf.add_block(br(0, 0, StreamRole::Sequence, 0, 3)).unwrap();
        buf.add_block(br(0, 0, StreamRole::Sequence, 1, 2)).unwrap();
        assert_eq!(buf.drain_ready().len(), 1);
        assert!(buf.all_emitted(1));
    }
}

#[cfg(test)]
mod qual_tests {
    use super::*;
    use crate::compression::codecs;

    /// Slicing shape + order: 23 records, cap 5 → [5,5,5,5,3], records preserved in order.
    #[test]
    fn fqz_stager_slices_by_cap_in_order() {
        let quals: Vec<Vec<u8>> = (0..23u8).map(|i| vec![b'!' + (i % 30); 4 + (i % 3) as usize]).collect();
        let mut st = QualStager::fqz(5, 0);
        let mut got_jobs: Vec<RawFqzBlock> = Vec::new();
        for q in &quals { if let Some(j) = st.stage_fqz(0, q) { got_jobs.push(j); } }
        if let Some(j) = st.flush_fqz(0) { got_jobs.push(j); }
        assert_eq!(got_jobs.iter().map(|j| j.quals.len()).collect::<Vec<_>>(), vec![5, 5, 5, 5, 3]);
        let flat: Vec<Vec<u8>> = got_jobs.iter().flat_map(|j| j.quals.clone()).collect();
        assert_eq!(flat, quals, "fqz stager preserves record order across slices");
        // block_idx is sequential 0..N
        assert_eq!(got_jobs.iter().map(|j| j.block_idx).collect::<Vec<_>>(), vec![0, 1, 2, 3, 4]);
    }

    /// Byte-identity: staging + compress_job per slice == legacy fqz_record_capped_blocks.
    #[test]
    fn fqz_blocks_match_record_capped() {
        let quals: Vec<Vec<u8>> = (0..23u8).map(|i| vec![b'!' + (i % 30); 4 + (i % 3) as usize]).collect();
        let refs: Vec<&[u8]> = quals.iter().map(|q| q.as_slice()).collect();
        let want = codecs::fqz_record_capped_blocks(&refs, 5).unwrap();

        let mut st = QualStager::fqz(5, 0);
        let mut got_jobs: Vec<RawFqzBlock> = Vec::new();
        for q in &quals { if let Some(j) = st.stage_fqz(0, q) { got_jobs.push(j); } }
        if let Some(j) = st.flush_fqz(0) { got_jobs.push(j); }
        let plan = tests_plan_default(); // Fqz arm ignores BSC plan fields
        let got: Vec<(u32, Vec<u8>)> = got_jobs.into_iter().map(|j| {
            let JobOutput::Blocks(r) = compress_job(CompressJob::Fqz(j), &plan).unwrap()
                else { panic!("expected blocks") };
            (r[0].record_count, r[0].compressed.clone())
        }).collect();
        assert_eq!(want, got, "fqz stager+compress must be byte-identical to fqz_record_capped_blocks");
    }

    #[test]
    fn bsc_stager_accessor() {
        assert!(QualStager::bsc(0, 1 << 20).bsc_stager().is_some());
        assert!(QualStager::fqz(5, 0).bsc_stager().is_none());
        assert!(QualStager::None.bsc_stager().is_none());
    }

    /// reset_for_next_chunk restarts block_idx at 0 — the producer relies on this for
    /// chunk-relative fqz numbering (matching legacy per-chunk framing).
    #[test]
    fn fqz_reset_restarts_block_idx_per_chunk() {
        let mut st = QualStager::fqz(5, 0);
        // chunk 0: 5 records → one block at idx 0, then flush (no tail) is None.
        let mut j0 = None;
        for q in 0..5u8 { if let Some(j) = st.stage_fqz(0, &[b'I', q]) { j0 = Some(j); } }
        let j0 = j0.expect("cap-5 slice emits at the 5th record");
        assert_eq!((j0.chunk_index, j0.block_idx), (0, 0));
        assert!(st.flush_fqz(0).is_none(), "no partial tail for an exact multiple of cap");

        st.reset_for_next_chunk();

        // chunk 1: 3 records → block_idx must restart at 0 (not continue from chunk 0).
        for q in 0..3u8 { assert!(st.stage_fqz(1, &[b'H', q]).is_none()); }
        let j1 = st.flush_fqz(1).expect("partial tail of 3 records");
        assert_eq!((j1.chunk_index, j1.block_idx, j1.quals.len()), (1, 0, 3));
    }
}

#[cfg(test)]
mod plan_tests {
    use super::*;
    use crate::cli::CompressConfig;
    use crate::io::FastqRecord;

    fn recs(n: usize, seq_len: usize) -> Vec<FastqRecord> {
        (0..n).map(|i| FastqRecord {
            id: format!("@r{i}").into_bytes(),
            sequence: vec![b'A'; seq_len],
            quality: Some(vec![b'I'; seq_len]),
        }).collect()
    }

    #[test]
    fn plan_matches_legacy_const_and_codec() {
        let mut cfg = CompressConfig::default();
        cfg.advanced.chunk_records = 4;
        let first = recs(4, 16); // small → Auto + <100K → BSC quality
        let plan = resolve_archive_plan(&cfg, ChunkedMode::Default, &first).unwrap();
        assert!(!plan.use_fqzcomp);
        assert_eq!(plan.const_seq_len, 16);
        assert_eq!(plan.const_qual_len, 16); // BSC quality path carries it
        assert!(plan.use_columnar_headers); // default header compressor = Columnar
        // Defaults (cli.rs): sequence_hints=false, rc_canon=false → encoding_type = Raw.
        assert_eq!(plan.encoding_type, EncodingType::Raw.to_u8());
    }

    #[test]
    fn plan_variable_length_zeroes_const() {
        let cfg = CompressConfig::default();
        let mut first = recs(3, 16);
        first[1].sequence = vec![b'A'; 20];
        first[1].quality = Some(vec![b'I'; 20]);
        let plan = resolve_archive_plan(&cfg, ChunkedMode::Default, &first).unwrap();
        assert_eq!((plan.const_seq_len, plan.const_qual_len), (0, 0));
    }

    /// Ultra mode forces UltraBigBlock encoding and clears hints/rc_canon — an
    /// invariant invisible to downstream roundtrip tests, pinned here in isolation.
    #[test]
    fn plan_ultra_forces_encoding_and_clears_hints() {
        let mut cfg = CompressConfig::default();
        cfg.advanced.sequence_hints = true; // would be RawWithHints in default mode
        let mode = ChunkedMode::Ultra {
            bsc_block_mb: 188,
            chunk_records: 4_000_000,
            quality_sub_block: 500_000,
            compress_window: 2,
        };
        let plan = resolve_archive_plan(&cfg, mode, &[]).unwrap();
        assert_eq!(plan.encoding_type, EncodingType::UltraBigBlock.to_u8());
        assert!(!plan.sequence_hints, "ultra forces sequence_hints off");
        assert!(!plan.rc_canon, "ultra forces rc_canon off");
        assert_eq!(plan.bsc_block_mb, 188); // ultra overrides the block size
    }

    /// no_quality skips the quality stream and silences binning.
    #[test]
    fn plan_no_quality_skips_quality_stream() {
        let cfg = CompressConfig { no_quality: true, ..CompressConfig::default() };
        let plan = resolve_archive_plan(&cfg, ChunkedMode::Default, &[]).unwrap();
        assert!(plan.no_quality);
        assert!(plan.skip_quality_in_stream);
        assert!(matches!(plan.stream_quality_binning, QualityBinning::None));
    }
}

#[cfg(test)]
mod worker_count_tests {
    use super::compress_worker_count;

    /// Default (25 MiB) block: the pool FOLLOWS the thread setting, capped by the 8 GiB
    /// safety budget at ~45 (8 GiB / (7 × 25 MiB) = 46) so a huge -t can't OOM.
    #[test]
    fn default_block_follows_threads_to_budget_cap() {
        assert_eq!(compress_worker_count(8, 25 * 1024 * 1024, false), 8); // follows -t
        assert_eq!(compress_worker_count(32, 25 * 1024 * 1024, false), 32); // follows -t
        // 8 GiB / (7 × 25 MiB) = 46 → the safety cap binds above ~46 threads.
        assert_eq!(compress_worker_count(72, 25 * 1024 * 1024, false), 46);
    }

    /// Ultra runs ONE block at a time regardless of -t (concurrent forward-BWT on big
    /// blocks blows up libsais memory; the block parallelises internally via OpenMP).
    #[test]
    fn ultra_is_single_worker_regardless_of_threads() {
        assert_eq!(compress_worker_count(72, 750 * 1024 * 1024, true), 1);
        assert_eq!(compress_worker_count(72, 188 * 1024 * 1024, true), 1);
        assert_eq!(compress_worker_count(8, 188 * 1024 * 1024, true), 1);
    }

    /// The safety budget binds for a large configured (non-ultra) block size: a 512 MiB
    /// block ⇒ 8 GiB / (7 × 512 MiB) = 2 workers even with 72 threads.
    #[test]
    fn large_block_size_hits_safety_budget() {
        assert_eq!(compress_worker_count(72, 512 * 1024 * 1024, false), 2);
    }
}
