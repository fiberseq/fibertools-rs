//! MA-spec annotation I/O bridge for fibertools-rs.
//!
//! Reading: prefers MA/AQ/AN; falls back to legacy `ns/nl/as/al/aq` (nuc/msp)
//! and legacy `fs/fl/fa` (fibertig) if no MA tag is present. Returns a
//! populated [`MolecularAnnotations`] (read length and aligned blocks set)
//! even when no annotation tags exist.
//!
//! Writing: emits MA-spec tags only. Legacy `ns/nl/as/al/aq` emission is
//! intentionally not supported — the new FIRE-as-subset semantics is
//! incompatible with the legacy dense `aq` layout. Pre-MA BAMs remain
//! readable; downstream consumers must migrate to the MA tag set.
//!
//! Annotation type names produced by fibertools-rs:
//! - `nuc`  (no strand, no quality)
//! - `msp`  (no strand, `Q` zero-valued — qualities live on `fire`)
//! - `fire` (no strand, `Q` linear precision 0–255)
//! - `fiberseq_callable` (no strand, no quality, one annotation per
//!   processed read): the span of the surviving nuc/MSP calls. A read
//!   that fails the minimums gets a zero-length annotation at position 0
//!   (NotCallable); only never-processed reads have none. Non-default
//!   minimums name the annotation in the AN tag (e.g. "m20a10").
//!
//! Hard-clipped records whose tags describe the full-length read (an aligner
//! copied the primary's tags onto a supplementary; see [`record_frame`]) keep
//! their nuc/msp/fire in that frame, lifted with the leading hard clip, and
//! lose their MM/ML. They are always NotCallable: they have no m6A.
//!
//! `m6a` and `cpg` types may appear *in memory* on a [`MolecularAnnotations`]
//! populated by the library's MM/ML parser. Their on-disk source of truth
//! is `MM`/`ML`; they carry `Encoding::MmMl` (set at construction, whether read
//! by the parser or built by a producer via `Encoding::mm_ml()`), so the
//! library serializes them into MM/ML rather than the MA tag set.

use anyhow::{bail, Result};
use molecular_annotation::{
    hard_clips, ma_family_tags, query_span, AlignedBlocks, Encoding, MolecularAnnotations,
    QualitySpec, Strand,
};
use rust_htslib::bam::{self, record::Aux};

/// Annotation type names used by fibertools-rs.
pub const NUC_TYPE: &str = "nuc";
pub const MSP_TYPE: &str = "msp";
pub const FIRE_TYPE: &str = "fire";
pub const FIBERSEQ_CALLABLE_TYPE: &str = "fiberseq_callable";

/// Molecular-orientation nuc/msp arrays returned by [`extract_nuc_msp_arrays`]:
/// `(nuc_starts, nuc_lengths, msp_starts, msp_lengths, msp_qual)`. `msp_qual`
/// is empty when there are no MSP qualities (pre-FIRE state).
type NucMspArrays = (Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>, Vec<u8>);

/// MSP input to [`build_annotations`]: `(starts, lengths, optional qualities)`
/// borrowed from caller-owned slices.
type MspInput<'a> = (&'a [u32], &'a [u32], Option<&'a [u8]>);

/// Reads annotations from a BAM record's tags.
///
/// Delegates to the library's combined MA-spec + MM/ML parser. If the
/// library returns an empty annotation set and the record carries legacy
/// `ns`/`nl`/`as`/`al`/`aq` tags, falls back to the legacy reader.
/// Tolerant of malformed MM/ML — the library handles those internally
/// without panicking.
///
/// The record's frame ([`record_frame`]) decides what is parsed: a stale
/// record comes back empty, a full-read-frame record keeps its MA/legacy
/// annotations under the full read length (lifted with the leading hard
/// clip) and never parses MM/ML, and everything else parses as before.
pub fn read_record(record: &bam::Record) -> Result<MolecularAnnotations> {
    let mut annot = match record_frame(record) {
        // A record whose tags describe another read is untagged, and parsing
        // its MM/ML would only produce truncation noise: decide before parsing.
        Frame::Stale(why) => {
            let mut annot = MolecularAnnotations::new(0);
            drop_stale_frame(&mut annot, record, why);
            return Ok(annot);
        }
        Frame::FullRead {
            h_lead,
            read_length,
        } => {
            let annot = if ma_family_tags(record).is_some() {
                // The library sets the query offset from the MA read length
                // and skips MM/ML itself.
                MolecularAnnotations::from_record(record)
            } else {
                // Legacy arrays record no frame, and every producer wrote
                // them on the full read: build that frame directly so MM/ML
                // are never decoded against the wrong SEQ.
                let mut annot = MolecularAnnotations::new(read_length);
                annot.set_aligned_blocks_raw(
                    AlignedBlocks::from_record(record).with_query_offset(h_lead),
                    record.is_reverse(),
                );
                read_legacy_nuc_msp_into(record, &mut annot)?;
                read_legacy_fibertig_into(record, &mut annot)?;
                annot
            };
            annot
        }
        Frame::Seq => {
            let mut annot = MolecularAnnotations::from_record(record);
            // The MA tag vouches for the frame, but an MN tag that disagrees
            // with SEQ says the base mods were copied from a longer read: keep
            // the calls, drop only the m6A/CpG (the writer strips MM/ML/MN).
            if mn_disagrees(record) {
                annot.annotation_types.retain(|t| !t.is_mm_ml());
            }
            // If MA tag is absent, also ingest legacy nuc/msp tags. The library
            // already populates basemod types (m6a/cpg) from MM/ML, so we add
            // legacy-derived nuc/msp to whatever the library produced rather
            // than gating on `annotation_types.is_empty()` (which would skip
            // the legacy fallback whenever MM/ML is present).
            if ma_family_tags(record).is_none() {
                // Provenance: the gates are type-checked (`has_legacy_nuc_msp`,
                // `has_legacy_fibertig`), so a foreign tool reusing these
                // two-letter names with a different aux type is never parsed
                // as fibertools data. The write path
                // (`strip_consumed_legacy_tags`) removes legacy tags under
                // pair-level gates mirroring `read_legacy_nuc_msp_into`'s
                // consumption, i.e. only what this reader ingested, so keep
                // the two in sync.
                if has_legacy_nuc_msp(record) {
                    read_legacy_nuc_msp_into(record, &mut annot)?;
                }
                // Legacy fibertig (`fs`/`fl`/`fa`): the pre-MA fibertig wire
                // format, dropped from the writer in favour of the MA-spec
                // `AN` tag. Kept readable so older fibertig BAMs stay
                // consumable by `extract`.
                if has_legacy_fibertig(record) {
                    read_legacy_fibertig_into(record, &mut annot)?;
                }
            }
            annot
        }
    };
    // Backstop on the parsed model: a read length that disagrees with the
    // frame, or any annotation ending past it.
    if let Some(why) = model_frame_reason(&annot, record) {
        drop_stale_frame(&mut annot, record, why);
    } else if model_is_full_read_frame(&annot, record) {
        // Warn and count only when this pass really drops m6A. The writer
        // strips MM/ML, so a second run over ft's own output stays quiet.
        if matches!(record.aux(b"MM"), Ok(Aux::String(_))) {
            note_full_read_frame(record);
        } else {
            log::debug!(
                "full-read frame for {} (no m6A to drop)",
                String::from_utf8_lossy(record.qname())
            );
        }
    }
    Ok(annot)
}

/// True when an MN tag is present and names a different length than SEQ:
/// the MM/ML next to it were written for another read.
fn mn_disagrees(record: &bam::Record) -> bool {
    let seq_len = record.seq_len();
    seq_len > 0
        && matches!(record.aux(b"MM"), Ok(Aux::String(_)))
        && record
            .aux(b"MN")
            .ok()
            .and_then(aux_as_usize)
            .is_some_and(|mn| mn != seq_len)
}

/// How many stale-frame records get a WARN before the rest drop to DEBUG.
/// ONT BAMs can hold thousands of hard-clipped supplementary reads. A total
/// is printed at exit by [`report_stale_frames`].
const STALE_FRAME_WARN_LIMIT: usize = 10;

/// Records whose annotations were dropped because their tags did not fit SEQ.
static STALE_FRAMES: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// What to tell a user who hit a stale frame. Aligners that hard-clip
/// supplementary alignments (minimap2 and dorado aligner without -Y) copy the
/// full-length read's tags onto the clipped record; the tags cannot be
/// recovered, only avoided.
pub const HARD_CLIP_REMEDY: &str = "Realign with soft clipping: pbmm2 align (PacBio; it never \
     hard-clips), dorado aligner --mm2-opts \"-Y\", or minimap2 -Y -y. Or drop supplementary \
     alignments with -F 2048.";

/// Treat a record whose tags come from another frame as untagged (#136, #31).
/// Hard-clipped supplementary alignments keep the full-length read's tags,
/// which would index past SEQ: consumers panic, or emit misplaced and
/// u32-wrapped coordinates. Lives here so every path that parses a record
/// (the fiber reader, convert-tags, strip-basemods, ddda-to-m6a, predict-m6a,
/// fibertig) sees the same thing, and [`write_record`] strips the stale tags
/// on the way out so nothing downstream trusts them.
///
/// Not reframed on purpose: nuc/msp could be shifted by the hard clip, but
/// MM/ML cannot be recovered without the clipped bases, and half a record
/// (nucleosomes, no m6A) is worse than an honest untagged one. The spec's
/// answer is soft clipping (`HARD_CLIP_REMEDY`).
fn drop_stale_frame(annot: &mut MolecularAnnotations, record: &bam::Record, why: String) {
    use std::sync::atomic::Ordering;
    let n = STALE_FRAMES.fetch_add(1, Ordering::Relaxed);
    let qname = String::from_utf8_lossy(record.qname());
    if n == 0 {
        // One message, so the cause never prints after the symptom when
        // several threads hit this at once.
        log::warn!(
            "some records carry Fiber-seq tags that describe the full-length read, not their SEQ \
             (hard-clipped supplementary alignments). Their annotations are dropped and they are \
             treated as untagged reads. {HARD_CLIP_REMEDY}\n\
             dropping annotations for {qname}: {why}"
        );
    } else if n < STALE_FRAME_WARN_LIMIT {
        log::warn!("dropping annotations for {qname}: {why}");
        if n + 1 == STALE_FRAME_WARN_LIMIT {
            log::warn!(
                "further such records are logged at debug level; a total is printed at exit"
            );
        }
    } else {
        log::debug!("dropping annotations for {qname}: {why}");
    }
    annot.annotation_types.clear();
    // An untagged read's frame is its SEQ, or its CIGAR query span without SEQ,
    // so the writers persist a clean `Ma:Z:<len>` instead of the stale length.
    // The liftover follows the frame (no query offset).
    annot.read_length = query_span(record);
    annot.set_aligned_blocks_raw(AlignedBlocks::from_record(record), record.is_reverse());
}

/// Records whose Fiber-seq tags describe the full-length read while they
/// hold a hard-clipped part of it (`Frame::FullRead`).
static FULL_FRAMES: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// What happens to the base mods of a full-read-frame record. Part of the
/// once-per-run warning, so a test can look for it.
pub const FULL_FRAME_M6A_DROPPED: &str =
    "their m6A (MM/ML) describes bases this record does not carry and was dropped";

/// Warn once per run that a record is in the full-read frame; later records
/// are logged at debug level and totalled at exit by [`report_full_frames`].
fn note_full_read_frame(record: &bam::Record) {
    use std::sync::atomic::Ordering;
    let n = FULL_FRAMES.fetch_add(1, Ordering::Relaxed);
    let qname = String::from_utf8_lossy(record.qname());
    if n == 0 {
        log::warn!(
            "some hard-clipped records carry Fiber-seq tags computed on the full-length read \
             (supplementary alignments; the aligner copied them from the primary). Their \
             nucleosome/MSP/FIRE calls are kept and lifted with the hard-clip offset; \
             {FULL_FRAME_M6A_DROPPED}, so these reads are not fiberseq-callable. \
             {HARD_CLIP_REMEDY}\n\
             full-read frame for {qname}"
        );
    } else {
        log::debug!("full-read frame for {qname}");
    }
}

/// Log the number of full-read-frame records. Called once at exit by
/// `main`, next to [`report_stale_frames`].
pub fn report_full_frames() {
    let n = FULL_FRAMES.load(std::sync::atomic::Ordering::Relaxed);
    if n > 0 {
        let records = if n == 1 { "record" } else { "records" };
        log::warn!(
            "kept nucleosome/MSP calls on {n} hard-clipped {records} whose Fiber-seq tags \
             describe the full-length read; they have no m6A and are not fiberseq-callable. \
             {HARD_CLIP_REMEDY}"
        );
    }
}

/// Log the number of records whose annotations were dropped for a stale
/// frame. Called once at exit by `main`.
pub fn report_stale_frames() {
    let n = STALE_FRAMES.load(std::sync::atomic::Ordering::Relaxed);
    if n > 0 {
        let records = if n == 1 { "record" } else { "records" };
        log::warn!(
            "dropped annotations on {n} {records} whose Fiber-seq tags did not match their SEQ \
             (hard-clipped supplementary alignments). {HARD_CLIP_REMEDY}"
        );
    }
}

/// True when the model is in the full-read frame of a hard-clipped record:
/// SEQ is present, the record has hard clips, and the model's read length is
/// SEQ plus both clips. Such a model has no m6A (the reader never parses
/// MM/ML in this frame) and its query coordinates run past SEQ; the lift
/// carries the leading hard clip as `annot.query_offset()`.
pub(crate) fn model_is_full_read_frame(annot: &MolecularAnnotations, record: &bam::Record) -> bool {
    let span = query_span(record) as usize;
    if span == 0 {
        return false;
    }
    let (lead, trail) = hard_clips(record);
    lead + trail > 0 && annot.read_length as usize == span + lead as usize + trail as usize
}

/// True when a present SEQ disagrees with the model's frame: the recorded
/// read length must equal SEQ, or SEQ plus the hard clips for a full-read
/// frame (something else rewrote the read after tagging). SEQ-less records
/// are never stale: the MA read length is the frame.
pub(crate) fn frame_is_stale(annot: &MolecularAnnotations, record: &bam::Record) -> bool {
    let seq_len = record.seq_len();
    seq_len > 0 && annot.read_length as usize != seq_len && !model_is_full_read_frame(annot, record)
}

/// Which read a record's Fiber-seq tags describe; see [`record_frame`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Frame {
    /// The tags describe SEQ (unclipped, soft-clipped, SEQ-less, or computed
    /// after clipping).
    Seq,
    /// The tags describe the full-length read of which SEQ is a hard-clipped
    /// part. `h_lead` is the leading hard clip (BAM/CIGAR orientation, the
    /// same end on both strands), `read_length` the full read's length.
    FullRead { h_lead: u32, read_length: u32 },
    /// The tags describe some other read: drop them.
    Stale(String),
}

/// Which read a record's Fiber-seq tags describe, from signals read straight
/// off the record. Each tag family has its own frame signal:
/// - MA: the tag's own read length. Equal to SEQ (a hard-clipped record whose
///   tags were computed after clipping included): `Seq`. Equal to SEQ plus
///   the hard clips: `FullRead`, the aligner copied the full-length read's
///   tags onto this clipped part of it. Anything else: `Stale`.
/// - legacy `ns`/`nl`/`as`/`al` and fibertig `fs`/`fl`: no frame is recorded
///   and every producer wrote them on the full read, so any hard clip means
///   `FullRead` with `read_length = seq_len + H_lead + H_trail`; a missing
///   SEQ is `Stale`, since nothing then anchors them.
/// - MM/ML (with no MA or legacy tags deciding first): the SAM `MN` tag when
///   present (the spec's frame for exactly this case), else hard clips.
///
/// In the full-read frame the MA-family coordinates are still right (they
/// are molecular coordinates of the full read; only the lift needs the
/// leading hard clip), but MM/ML are SEQ-relative and cannot be recovered:
/// the reader drops them and the writer strips them.
///
/// SEQ-less MA records are never stale: the MA read length is the frame.
/// Both the reader ([`read_record`]) and the writer ([`write_record`]) use
/// this, so a stale record is cleared on the way in and cleaned on the way
/// out.
pub(crate) fn record_frame(record: &bam::Record) -> Frame {
    let seq_len = record.seq_len();
    let (lead, trail) = hard_clips(record);
    let hard_clipped = lead + trail > 0;
    // The bases this record holds: SEQ, or the CIGAR query span when SEQ was
    // dropped. A hard-clipped record's frame is judged against that span so
    // a SEQ-less supplementary still gets the offset lift.
    let span = query_span(record) as usize;
    let full = span + lead as usize + trail as usize;
    if let Some((ma, _, _)) = ma_family_tags(record) {
        if let Some(read_length) = ma.split(';').next().and_then(|s| s.parse::<usize>().ok()) {
            if hard_clipped && span > 0 && read_length == full {
                return Frame::FullRead {
                    h_lead: lead,
                    read_length: read_length as u32,
                };
            }
            if seq_len > 0 && read_length != seq_len {
                if hard_clipped && read_length == full {
                    return Frame::FullRead {
                        h_lead: lead,
                        read_length: read_length as u32,
                    };
                }
                if hard_clipped {
                    return Frame::Stale(format!(
                        "MA read length {read_length} matches neither the {seq_len} bp sequence \
                         nor the {full} bp read it was hard-clipped from"
                    ));
                }
                return Frame::Stale(format!(
                    "MA read length {read_length} does not match the {seq_len} bp sequence"
                ));
            }
            // The MA tag was written for this SEQ, so it vouches for the
            // MM/ML next to it too: fibertools' own writers emit MA and
            // MM/ML together, and a hard-clipped record it produced is fine.
            return Frame::Seq;
        }
    } else if has_legacy_nuc_msp(record) || has_legacy_fibertig(record) {
        if seq_len == 0 {
            return Frame::Stale("legacy nuc/msp tags on a record without SEQ".to_string());
        }
        if hard_clipped {
            // Any MM/ML next to them are SEQ-relative copies from the full
            // read, dropped with the frame; they are not a stale signal.
            return Frame::FullRead {
                h_lead: lead,
                read_length: full as u32,
            };
        }
    }
    if matches!(record.aux(b"MM"), Ok(Aux::String(_))) {
        if seq_len == 0 {
            // Positions are implicit in SEQ; without it there is nothing to
            // decode against (minimap2 -y writes SEQ-less secondaries).
            return Frame::Stale("MM/ML on a record without SEQ".to_string());
        }
        match record.aux(b"MN").ok().and_then(aux_as_usize) {
            Some(mn) => {
                if mn != seq_len {
                    return Frame::Stale(format!(
                        "MN {mn} does not match the {seq_len} bp sequence"
                    ));
                }
            }
            None => {
                if hard_clipped {
                    return Frame::Stale(
                        "MM/ML on a hard-clipped alignment without an MN tag".to_string(),
                    );
                }
            }
        }
    }
    Frame::Seq
}

/// Why a record's tags describe a different read than its SEQ, or `None`
/// when they describe SEQ or the full read it was hard-clipped from.
pub(crate) fn record_frame_reason(record: &bam::Record) -> Option<String> {
    match record_frame(record) {
        Frame::Stale(why) => Some(why),
        _ => None,
    }
}

fn aux_as_usize(aux: Aux) -> Option<usize> {
    match aux {
        Aux::I8(v) => usize::try_from(v).ok(),
        Aux::U8(v) => Some(v as usize),
        Aux::I16(v) => usize::try_from(v).ok(),
        Aux::U16(v) => Some(v as usize),
        Aux::I32(v) => usize::try_from(v).ok(),
        Aux::U32(v) => Some(v as usize),
        _ => None,
    }
}

/// The parsed model disagrees with its frame: a read length that matches
/// neither SEQ nor the full read SEQ was hard-clipped from, or an annotation
/// ending past the frame. Backstop behind [`record_frame`] for tags that
/// carry no frame signal of their own.
fn model_frame_reason(annot: &MolecularAnnotations, record: &bam::Record) -> Option<String> {
    let seq_len = record.seq_len();
    if seq_len == 0 {
        return None;
    }
    if frame_is_stale(annot, record) {
        return Some(format!(
            "MA read length {} does not match the {seq_len} bp sequence",
            annot.read_length
        ));
    }
    // SEQ, or the full read in the full-read frame.
    let frame = annot.read_length as usize;
    annot
        .annotation_types
        .iter()
        .find(|t| {
            t.annotations
                .iter()
                .any(|a| a.start as usize + a.length as usize > frame)
        })
        .map(|t| format!("{} coordinates extend past the {frame} bp read", t.name))
}

/// Why a record's annotations do not fit its SEQ, or `None` when they do:
/// the record-level signals first, then the parsed model as a backstop.
#[cfg(test)]
pub(crate) fn stale_frame_reason(
    annot: &MolecularAnnotations,
    record: &bam::Record,
) -> Option<String> {
    record_frame_reason(record).or_else(|| model_frame_reason(annot, record))
}

/// Make the fiberseq_callable annotation agree with the CLI minimums.
/// Runs right after parsing, before any consumer-side pruning. Derives
/// when the tag is absent (backfill) or the minimums are non-default
/// (recalculation); otherwise the on-disk tag stands.
///
/// A full-read-frame record (hard clips, tags from the full read) has no
/// m6A: the reader dropped its MM/ML. FIRE cannot score it, so it is
/// NotCallable whatever its on-disk tag or the minimums say. That is the
/// only "no m6A" this function knows about; m6A is never inspected, so a
/// read whose MM/ML were stripped on purpose (ft strip-basemods) keeps the
/// verdict from calling time.
pub fn sync_fiberseq_callable(
    annot: &mut MolecularAnnotations,
    record: &bam::Record,
    filters: &crate::utils::input_bam::FiberFilters,
) {
    let (min_msp, min_ave) = filters.callable_minimums();
    if model_is_full_read_frame(annot, record) {
        // No m6A means not callable. Without nuc/msp there is nothing to
        // judge either, so a never-called read gets no marker, but a copied
        // Callable span must not survive the frame it no longer describes.
        if has_calls(annot) || annot.get_type(FIBERSEQ_CALLABLE_TYPE).is_some() {
            set_fiberseq_callable(annot, None, callable_minimums_name(min_msp, min_ave));
        }
        return;
    }
    let needed =
        filters.callable_minimums_are_custom() || annot.get_type(FIBERSEQ_CALLABLE_TYPE).is_none();
    if needed && can_derive_callable(annot, record) {
        derive_fiberseq_callable(annot, min_msp, min_ave);
    }
}

/// True when calling ran: nuc or msp present.
fn has_calls(annot: &MolecularAnnotations) -> bool {
    annot.get_type(NUC_TYPE).is_some() || annot.get_type(MSP_TYPE).is_some()
}

/// True when the callable state can be derived: calling ran (nuc or msp
/// present) and the frame is not stale. Derivation is pure MA-tag
/// arithmetic, so SEQ-less records derive fine.
fn can_derive_callable(annot: &MolecularAnnotations, record: &bam::Record) -> bool {
    !frame_is_stale(annot, record) && has_calls(annot)
}

/// AN name recording non-default minimums, e.g. "m20a10". `None` for the
/// defaults: unnamed means default minimums, so a default-run BAM carries
/// no AN tag. Comma-free (the AN tag joins names with commas).
pub fn callable_minimums_name(min_msp: usize, min_ave_msp_size: i64) -> Option<String> {
    use crate::utils::input_bam::{FIRE_CALLABLE_MIN_AVE_MSP_SIZE, FIRE_CALLABLE_MIN_MSP};
    if min_msp == FIRE_CALLABLE_MIN_MSP && min_ave_msp_size == FIRE_CALLABLE_MIN_AVE_MSP_SIZE {
        None
    } else {
        Some(format!("m{min_msp}a{min_ave_msp_size}"))
    }
}

/// Non-default minimums name the annotation (see [`callable_minimums_name`]),
/// so the AN tag records the provenance on disk; default minimums stay
/// unnamed and adds no AN bytes.
pub fn set_fiberseq_callable(
    annot: &mut MolecularAnnotations,
    span: Option<(u32, u32)>,
    minimums_name: Option<String>,
) {
    annot
        .annotation_types
        .retain(|t| t.name != FIBERSEQ_CALLABLE_TYPE);
    let t = annot.add_annotation_type(FIBERSEQ_CALLABLE_TYPE, QualitySpec::none(), Encoding::Ma);
    let (s, l) = match span {
        Some((s, e)) => (s, e - s),
        None => (0, 0),
    };
    t.add(s, l, Strand::Unknown, vec![], minimums_name);
}

/// Derive the callable span and state from the nuc/MSP annotations; the
/// single source of truth for every writer. The span is the extent of the
/// surviving calls; callable requires at least `min_msp` MSPs with a mean
/// length of at least `min_ave_msp_size`. m6A is never inspected: the caller
/// emits no MSP without m6A, and MA-only derivation is what lets SEQ-less
/// records derive. A full-read-frame record (no m6A) is forced NotCallable
/// by [`sync_fiberseq_callable`], not here.
pub fn derive_fiberseq_callable(
    annot: &mut MolecularAnnotations,
    min_msp: usize,
    min_ave_msp_size: i64,
) {
    let nuc = annot.get_forward_coords(NUC_TYPE).unwrap_or_default();
    let msp = annot.get_forward_coords(MSP_TYPE).unwrap_or_default();

    // min of the first starts, max of the last ends; both vectors are
    // ascending and non-overlapping, and either may be empty.
    let start = [nuc.first(), msp.first()]
        .into_iter()
        .flatten()
        .map(|(s, _)| *s)
        .min();
    let end = [nuc.last(), msp.last()]
        .into_iter()
        .flatten()
        .map(|(_, e)| *e)
        .max();
    let span = start.zip(end);

    let callable = msp.len() >= min_msp
        && !msp.is_empty()
        && msp.iter().map(|(s, e)| (e - s) as i64).sum::<i64>() / msp.len() as i64
            >= min_ave_msp_size;

    set_fiberseq_callable(
        annot,
        if callable { span } else { None },
        callable_minimums_name(min_msp, min_ave_msp_size),
    );
}

/// Writes MA-family tags (MA/AQ/AN) to a BAM record, **preserving the
/// record's existing MM/ML bytes**.
///
/// # Which write function do I call?
///
/// The two write paths differ only in how they treat base modifications
/// (MM/ML). Pick by asking: *does this code path create, modify, or remove
/// base mods?*
///
/// | Your code path | Call | MM/ML behavior |
/// |---|---|---|
/// | Edits only nuc/msp/fire (structural annotations) | [`write_record`] | preserved byte-identically |
/// | Creates / modifies / removes base mods | [`write_record_with_basemods`] | canonically re-encoded from the model |
///
/// The distinction exists because the in-memory model **normalizes** MM/ML on
/// read (splits grouped codes, canonicalizes skip flags), so re-encoding from
/// the model is lossy for spec-legal but non-canonical encodings (grouped
/// multi-code `C+mh`, `N+a` wildcards). Leaving MM/ML untouched is therefore
/// the only way to guarantee byte-identical round-trips for structural editors.
/// Calling the wrong one is silently wrong: a structural editor that re-encodes
/// corrupts non-canonical MM/ML; a producer that preserves emits stale mods.
///
/// A future refactor (see `molecular-annotation/docs/mm-ml-per-group-passthrough.md`)
/// would collapse these into one function by tracking dirtiness per MM group.
///
/// ---
///
/// This is the read/edit write path: it does not re-encode base modifications,
/// so any MM/ML present on the input survives byte-identically (including
/// spec-legal encodings the normalized model can't represent, such as grouped
/// multi-code `C+mh`). Use this for subcommands that add or edit non-basemod
/// annotations (nuc/msp/fire) but do not themselves produce or remove base
/// mods.
///
/// Producers that create or modify base mods must instead call
/// [`write_record_with_basemods`], which canonically re-emits MM/ML.
pub fn write_record(record: &mut bam::Record, annot: &MolecularAnnotations) {
    match record_frame(record) {
        // The reader cleared this record's annotations (drop_stale_frame);
        // leave no stale tag behind for another tool to trust.
        Frame::Stale(_) => strip_all_fiber_tags(record),
        Frame::FullRead { .. } => {
            // The MA tag is rewritten from the model below, in the full
            // read's frame. The SEQ-relative base mods describe bases this
            // record does not carry: the reader dropped them, so leave none
            // behind for another tool to trust.
            for tag in [b"MM", b"ML", b"MN"] {
                record.remove_aux(tag).ok();
            }
            strip_consumed_legacy_tags(record);
        }
        Frame::Seq => {
            if mn_disagrees(record) {
                // The reader dropped these base mods (see read_record).
                for tag in [b"MM", b"ML", b"MN"] {
                    record.remove_aux(tag).ok();
                }
            }
            strip_consumed_legacy_tags(record)
        }
    }
    annot.to_record(record);
}

/// Every Fiber-seq tag fibertools knows how to read: legacy nuc/msp/fibertig
/// arrays, MM/ML/MN base mods, and both spellings of the MA family. Used only
/// for records whose frame is stale (`record_frame`), where none of them
/// describe this SEQ.
fn strip_all_fiber_tags(record: &mut bam::Record) {
    for tag in [
        b"ns", b"nl", b"as", b"al", b"aq", b"fs", b"fl", b"fa", b"MM", b"ML", b"MN", b"Ma", b"Aq",
        b"An", b"MA", b"AQ", b"AN",
    ] {
        record.remove_aux(tag).ok();
    }
}

/// Provenance rule shared by every MA write: strip exactly the legacy tags
/// that [`read_record`] consumed as the source of the annotation model being
/// written — they are superseded by the MA-family tags (v0.9 replace
/// semantics; otherwise legacy readers silently see stale calls forever).
///
/// The gates mirror [`read_legacy_nuc_msp_into`]'s consumption at PAIR level
/// (ns+nl, as+al, aq only inside the msp pair, fa only with fs+fl), which is
/// what makes this provenance-safe: a tag is only removed when the reader
/// ingested it. Records that already carry an MA-family main tag were not
/// read from legacy tags, so their legacy-named aux tags are left untouched:
/// a foreign tool reusing those names is never destroyed. Likewise a lone or
/// wrong-typed tag outside its pair gate (e.g. `as` with no `al`, or `aq`
/// with no msp pair) is left alone.
fn strip_consumed_legacy_tags(record: &mut bam::Record) {
    if ma_family_tags(record).is_some() {
        return;
    }
    // The reader consumes legacy tags only when the ENTIRE legacy parse
    // succeeds: any validation bail (a length-mismatched pair) errors
    // read_record, and the writers that catch that error fall back to
    // writing an empty model — nothing was consumed, so nothing may be
    // removed. Gating on the same parse keeps strip and read consumption
    // identical by construction.
    let mut scratch = MolecularAnnotations::new(0);
    if read_legacy_nuc_msp_into(record, &mut scratch).is_err()
        || read_legacy_fibertig_into(record, &mut scratch).is_err()
    {
        return;
    }
    // Pair-level gates mirroring read_legacy_nuc_msp_into: ns+nl only as a pair,
    // as+al only as a pair, aq only inside a valid msp pair (its count is
    // validated by the parse above). A lone or orphan tag was never
    // consumed and so is never removed.
    if u32_array(record, b"ns").is_some() && u32_array(record, b"nl").is_some() {
        record.remove_aux(b"ns").ok();
        record.remove_aux(b"nl").ok();
    }
    if u32_array(record, b"as").is_some() && u32_array(record, b"al").is_some() {
        record.remove_aux(b"as").ok();
        record.remove_aux(b"al").ok();
        if u8_array(record, b"aq").is_some() {
            record.remove_aux(b"aq").ok();
        }
    }
    // fibertig: fs+fl as a pair; fa is only consumed (and so only stripped)
    // when the pair is non-empty, mirroring read_legacy_fibertig_into's early
    // return on empty fs.
    let fs = u32_array(record, b"fs");
    if fs.is_some() && u32_array(record, b"fl").is_some() {
        let non_empty = fs.as_ref().is_some_and(|v| !v.is_empty());
        record.remove_aux(b"fs").ok();
        record.remove_aux(b"fl").ok();
        if non_empty && matches!(record.aux(b"fa"), Ok(Aux::String(_))) {
            record.remove_aux(b"fa").ok();
        }
    }
}

/// Type-checked consumption gate for the legacy nuc/msp tag set. Presence
/// alone is weak provenance — a foreign tool could reuse these two-letter
/// names with a different aux type — so the set only counts as
/// fibertools-legacy when a gate tag parses as the int array legacy
/// fibertools wrote.
fn has_legacy_nuc_msp(record: &bam::Record) -> bool {
    u32_array(record, b"ns").is_some() || u32_array(record, b"as").is_some()
}

/// Type-checked consumption gate for the legacy fibertig tag set.
fn has_legacy_fibertig(record: &bam::Record) -> bool {
    u32_array(record, b"fs").is_some() && u32_array(record, b"fl").is_some()
}

/// Writes MA-family tags **and** canonically re-encodes MM/ML from the
/// annotation model.
///
/// See [`write_record`] for the "which write function do I call?" contract;
/// this is the base-mod-producer path.
///
/// For producers/synthesizers (`predict_m6a`, `ddda_to_m6a`, `strip_basemods`,
/// and record-synthesizing paths) that genuinely create, modify, or remove
/// base modifications. Any pre-existing MM/ML on the record is replaced by the
/// canonical encoding of the model's `Encoding::MmMl` types; if the model has
/// none, MM/ML are removed.
///
/// Basemod types must already be `Encoding::MmMl` — they are when read from a
/// record, and producers create them that way (e.g. `Encoding::mm_ml()`).
pub fn write_record_with_basemods(record: &mut bam::Record, annot: &MolecularAnnotations) {
    write_record(record, annot);
    annot.write_mm_ml(record);
    // MN is the SAM spec's frame for MM/ML: the SEQ length they were written
    // against. With it, any consumer (samtools, modkit, this reader) can tell
    // when a later hard clip has made them stale.
    record.remove_aux(b"MN").ok();
    let seq_len = record.seq_len();
    if seq_len > 0 && matches!(record.aux(b"MM"), Ok(Aux::String(_))) {
        record.push_aux(b"MN", Aux::I32(seq_len as i32)).ok();
    }
}

/// Read annotations from a BAM record.
///
/// Prefers spec MA tags; falls back to legacy `ns`/`nl`/`as`/`al`/`aq`. Always
/// returns a populated [`MolecularAnnotations`] (with read length and aligned
/// blocks set) — even if no annotation tags are present.
pub fn read_annotations(record: &bam::Record) -> Result<MolecularAnnotations> {
    if let Some(annot) = read_ma_tags(record)? {
        return Ok(annot);
    }
    read_legacy_nuc_msp(record)
}

fn read_ma_tags(record: &bam::Record) -> Result<Option<MolecularAnnotations>> {
    let Some((ma, aq, an)) = ma_family_tags(record) else {
        return Ok(None);
    };
    let mut annot = MolecularAnnotations::from_tags(&ma, aq.as_deref(), an.as_deref())
        .map_err(|e| anyhow::anyhow!("MA tag parse error: {e}"))?;
    // Same frame rule as the library's from_record: a full-read MA tag on a
    // hard-clipped record lifts with the leading hard clip.
    let offset = molecular_annotation::full_read_query_offset(annot.read_length, record);
    annot.set_aligned_blocks_raw(
        AlignedBlocks::from_record(record).with_query_offset(offset.unwrap_or(0)),
        record.is_reverse(),
    );
    Ok(Some(annot))
}

fn read_legacy_nuc_msp(record: &bam::Record) -> Result<MolecularAnnotations> {
    let mut annot = MolecularAnnotations::from_record(record);
    read_legacy_nuc_msp_into(record, &mut annot)?;
    Ok(annot)
}

/// Add the legacy `ns`/`nl`/`as`/`al`/`aq` tags to `annot` as nuc/msp/fire
/// types. The model's frame (read length, aligned blocks) is the caller's.
fn read_legacy_nuc_msp_into(record: &bam::Record, annot: &mut MolecularAnnotations) -> Result<()> {
    let ns = u32_array(record, b"ns");
    let nl = u32_array(record, b"nl");
    let a_starts = u32_array(record, b"as");
    let a_lens = u32_array(record, b"al");
    let aq = u8_array(record, b"aq");

    if let (Some(ns), Some(nl)) = (ns.as_ref(), nl.as_ref()) {
        if ns.len() != nl.len() {
            bail!(
                "legacy ns ({}) and nl ({}) length mismatch",
                ns.len(),
                nl.len()
            );
        }
        if !ns.is_empty() {
            let nuc = annot.add_annotation_type(NUC_TYPE, QualitySpec::none(), Encoding::Ma);
            for (s, l) in ns.iter().zip(nl.iter()) {
                nuc.add(*s, *l, Strand::Unknown, vec![], None);
            }
        }
    }

    if let (Some(starts), Some(lens)) = (a_starts.as_ref(), a_lens.as_ref()) {
        if starts.len() != lens.len() {
            bail!(
                "legacy as ({}) and al ({}) length mismatch",
                starts.len(),
                lens.len()
            );
        }
        if let Some(ref q) = aq {
            if q.len() != starts.len() {
                bail!(
                    "legacy aq ({}) and as ({}) length mismatch",
                    q.len(),
                    starts.len()
                );
            }
        }
        // MSPs always ingest without quality: the legacy `aq` byte is the FIRE
        // precision, which now lives on a separate `fire` type
        if !starts.is_empty() {
            let msp = annot.add_annotation_type(MSP_TYPE, QualitySpec::none(), Encoding::Ma);
            for (s, l) in starts.iter().zip(lens.iter()) {
                msp.add(*s, *l, Strand::Unknown, vec![], None);
            }
        }
        // FIRE is the subset of MSPs whose legacy precision crossed the
        // regulatory-element threshold (precision > 0), carrying that precision
        // as the `fire` quality — the same `p > 0` rule the FIRE producer uses.
        if let Some(ref q) = aq {
            let fire: Vec<(u32, u32, u8)> = starts
                .iter()
                .zip(lens.iter())
                .zip(q.iter())
                .filter(|(_, &p)| p > 0)
                .map(|((s, l), &p)| (*s, *l, p))
                .collect();
            if !fire.is_empty() {
                let q_spec = "Q".parse::<QualitySpec>()?;
                let fire_t = annot.add_annotation_type(FIRE_TYPE, q_spec, Encoding::Ma);
                for (s, l, p) in fire {
                    fire_t.add(s, l, Strand::Unknown, vec![p], None);
                }
            }
        }
    }

    Ok(())
}

/// Read the legacy fibertig `fs`/`fl`/`fa` tags into a [`FIBERTIG_TYPE`]
/// annotation type.
///
/// Pre-MA fibertig BAMs stored contig annotations in three tags: `fs`
/// (molecular starts), `fl` (lengths), and `fa` (a single `|`-separated string
/// of names, where an empty segment means "unnamed"). Emission of these was
/// dropped in favour of the MA-spec `AN` tag; this read-only fallback keeps
/// those older BAMs consumable. Mirrors the in-memory shape the MA path
/// produces (`Strand::Forward`, no quality, `Encoding::Ma`) so downstream
/// liftover to reference coordinates is identical either way.
fn read_legacy_fibertig_into(record: &bam::Record, annot: &mut MolecularAnnotations) -> Result<()> {
    use crate::utils::fibertig::FIBERTIG_TYPE;

    let (Some(fs), Some(fl)) = (u32_array(record, b"fs"), u32_array(record, b"fl")) else {
        return Ok(());
    };
    if fs.len() != fl.len() {
        bail!(
            "legacy fibertig fs ({}) and fl ({}) length mismatch",
            fs.len(),
            fl.len()
        );
    }
    if fs.is_empty() {
        return Ok(());
    }

    // `fa` is optional; when present it must have one `|`-separated segment per
    // annotation. An empty segment decodes to `None` (an unnamed annotation).
    let names: Option<Vec<Option<String>>> = match record.aux(b"fa") {
        Ok(Aux::String(s)) => {
            let parts: Vec<Option<String>> = s
                .split('|')
                .map(|p| (!p.is_empty()).then(|| p.to_string()))
                .collect();
            if parts.len() != fs.len() {
                bail!(
                    "legacy fibertig fa ({}) and fs ({}) length mismatch",
                    parts.len(),
                    fs.len()
                );
            }
            Some(parts)
        }
        _ => None,
    };

    let t = annot.add_annotation_type(FIBERTIG_TYPE, QualitySpec::none(), Encoding::Ma);
    for (i, (s, l)) in fs.iter().zip(fl.iter()).enumerate() {
        let name = names.as_ref().and_then(|v| v[i].clone());
        t.add(*s, *l, Strand::Forward, vec![], name);
    }
    Ok(())
}

/// Convenience for callers that still operate on raw `i64` arrays of
/// nucleosome and MSP coordinates. Returns
/// `(nuc_starts, nuc_lengths, msp_starts, msp_lengths, msp_qual)` derived
/// from whichever tag set is present (MA wins over legacy).
///
/// Coordinates are in molecular orientation — the same convention as the
/// legacy `ns`/`nl`/`as`/`al` tags. `msp_qual` is empty when there are no
/// MSP qualities (pre-FIRE state); otherwise one byte per MSP.
pub fn extract_nuc_msp_arrays(record: &bam::Record) -> Result<NucMspArrays> {
    let annot = read_record(record)?;
    let (nuc_starts, nuc_lengths) = annot
        .get_type(NUC_TYPE)
        .map(|t| {
            (
                t.annotations.iter().map(|a| a.start as i64).collect(),
                t.annotations.iter().map(|a| a.length as i64).collect(),
            )
        })
        .unwrap_or_default();
    let (msp_starts, msp_lengths, msp_qual) = annot
        .get_type(MSP_TYPE)
        .map(|t| {
            let starts: Vec<i64> = t.annotations.iter().map(|a| a.start as i64).collect();
            let lens: Vec<i64> = t.annotations.iter().map(|a| a.length as i64).collect();
            let qs: Vec<u8> = if t.quality_spec.has_quality() {
                t.annotations
                    .iter()
                    .map(|a| crate::utils::bamannotations::primary_qual(&a.qualities, MSP_TYPE))
                    .collect()
            } else {
                Vec::new()
            };
            (starts, lens, qs)
        })
        .unwrap_or_default();
    Ok((nuc_starts, nuc_lengths, msp_starts, msp_lengths, msp_qual))
}

/// Add `nuc` annotations (unknown strand, no quality) to `annot`.
///
/// No-op if `starts` is empty. `starts` and `lens` must be the same length
/// and paired positionally.
pub fn add_nuc_annotations(annot: &mut MolecularAnnotations, starts: &[u32], lens: &[u32]) {
    if starts.is_empty() {
        return;
    }
    let t = annot.add_annotation_type(NUC_TYPE, QualitySpec::none(), Encoding::Ma);
    for (s, l) in starts.iter().zip(lens.iter()) {
        t.add(*s, *l, Strand::Unknown, vec![], None);
    }
}

/// Add `msp` annotations (unknown strand) to `annot`.
///
/// With `qualities` present, the type uses `msp+Q`; without, `msp+` (no
/// quality). No-op if `starts` is empty. All input slices must be the same
/// length and paired positionally.
pub fn add_msp_annotations(
    annot: &mut MolecularAnnotations,
    starts: &[u32],
    lens: &[u32],
    qualities: Option<&[u8]>,
) {
    if starts.is_empty() {
        return;
    }
    let qspec = match qualities {
        Some(_) => "Q".parse::<QualitySpec>().expect("Q parses"),
        None => QualitySpec::none(),
    };
    let t = annot.add_annotation_type(MSP_TYPE, qspec, Encoding::Ma);
    for (i, (s, l)) in starts.iter().zip(lens.iter()).enumerate() {
        let qv = qualities.map(|q| vec![q[i]]).unwrap_or_default();
        t.add(*s, *l, Strand::Unknown, qv, None);
    }
}

/// Add `fire` annotations (no strand, linear precision 0–255) to `annot`.
///
/// No-op if `starts` is empty. All input slices must be the same length and
/// paired positionally.
pub fn add_fire_annotations(
    annot: &mut MolecularAnnotations,
    starts: &[u32],
    lens: &[u32],
    precisions: &[u8],
) {
    if starts.is_empty() {
        return;
    }
    let t = annot.add_annotation_type(
        FIRE_TYPE,
        "Q".parse::<QualitySpec>().expect("Q parses"),
        Encoding::Ma,
    );
    for (i, (s, l)) in starts.iter().zip(lens.iter()).enumerate() {
        t.add(*s, *l, Strand::Unknown, vec![precisions[i]], None);
    }
}

pub fn build_annotations(
    record: &bam::Record,
    nuc: Option<(&[u32], &[u32])>,
    msp: Option<MspInput>,
    fire: Option<(&[u32], &[u32], &[u8])>,
) -> MolecularAnnotations {
    let mut annot = MolecularAnnotations::from_record(record);
    if let Some((starts, lens)) = nuc {
        add_nuc_annotations(&mut annot, starts, lens);
    }
    if let Some((starts, lens, q)) = msp {
        add_msp_annotations(&mut annot, starts, lens, q);
    }
    if let Some((starts, lens, precisions)) = fire {
        add_fire_annotations(&mut annot, starts, lens, precisions);
    }
    annot
}

/// Backwards-compatible wrapper around [`build_annotations`]. Existing
/// nuc/msp producers can keep calling this until they migrate to the
/// fire-aware [`build_annotations`].
pub fn build_nuc_msp_annotations(
    record: &bam::Record,
    nuc_starts: &[u32],
    nuc_lengths: &[u32],
    msp_starts: &[u32],
    msp_lengths: &[u32],
    msp_qual: Option<&[u8]>,
) -> MolecularAnnotations {
    let nuc = (!nuc_starts.is_empty()).then_some((nuc_starts, nuc_lengths));
    let msp = (!msp_starts.is_empty()).then_some((msp_starts, msp_lengths, msp_qual));
    build_annotations(record, nuc, msp, None)
}

/// Convenience: read fire annotations as raw arrays in molecular orientation.
/// Returns `(starts, lengths, precisions)`. All three are empty if the record
/// has no `fire` MA type.
pub fn extract_fire_arrays(record: &bam::Record) -> Result<(Vec<u32>, Vec<u32>, Vec<u8>)> {
    let annot = read_record(record)?;
    Ok(annot
        .get_type(FIRE_TYPE)
        .map(|t| {
            let starts = t.annotations.iter().map(|a| a.start).collect();
            let lens = t.annotations.iter().map(|a| a.length).collect();
            let precisions = t
                .annotations
                .iter()
                .map(|a| crate::utils::bamannotations::primary_qual(&a.qualities, FIRE_TYPE))
                .collect();
            (starts, lens, precisions)
        })
        .unwrap_or_default())
}

fn u32_array(record: &bam::Record, tag: &[u8]) -> Option<Vec<u32>> {
    match record.aux(tag) {
        Ok(Aux::ArrayU32(arr)) => Some(arr.iter().collect()),
        Ok(Aux::ArrayI32(arr)) => Some(arr.iter().map(|v| v as u32).collect()),
        _ => None,
    }
}

fn u8_array(record: &bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match record.aux(tag) {
        Ok(Aux::ArrayU8(arr)) => Some(arr.iter().collect()),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use molecular_annotation::{MolecularAnnotations, QualitySpec, Strand};

    /// The strip must remove exactly what the reader consumed: foreign
    /// (wrong-typed) tags, lone pair members, and orphan aq must all
    /// survive; consumed pairs are removed; MA-family presence blocks all
    /// stripping.
    #[test]
    fn strip_consumed_legacy_tags_mirrors_consumption() {
        // foreign: `ns` as a string, `al` as a float — not fibertools-typed
        let mut foreign = synth_record(b"ACGTACGT");
        foreign.push_aux(b"ns", Aux::String("not-ours")).unwrap();
        foreign.push_aux(b"al", Aux::Float(1.5)).unwrap();
        strip_consumed_legacy_tags(&mut foreign);
        assert!(foreign.aux(b"ns").is_ok(), "foreign string ns was stripped");
        assert!(foreign.aux(b"al").is_ok(), "foreign float al was stripped");

        // consumed nuc pair is stripped; a typed-but-orphan aq (no msp pair)
        // and a lone as (no al) were never ingested and must survive
        let mut legacy = synth_record(b"ACGTACGT");
        legacy
            .push_aux(b"ns", Aux::ArrayU32((&vec![1u32, 5]).into()))
            .unwrap();
        legacy
            .push_aux(b"nl", Aux::ArrayU32((&vec![2u32, 2]).into()))
            .unwrap();
        legacy
            .push_aux(b"aq", Aux::ArrayU8((&vec![200u8, 0]).into()))
            .unwrap();
        legacy
            .push_aux(b"as", Aux::ArrayU32((&vec![3u32]).into()))
            .unwrap();
        strip_consumed_legacy_tags(&mut legacy);
        assert!(legacy.aux(b"ns").is_err(), "consumed ns not stripped");
        assert!(legacy.aux(b"nl").is_err(), "consumed nl not stripped");
        assert!(legacy.aux(b"aq").is_ok(), "orphan aq was stripped");
        assert!(legacy.aux(b"as").is_ok(), "lone as was stripped");

        // a consumed msp pair takes its aq with it
        let mut msp = synth_record(b"ACGTACGT");
        msp.push_aux(b"as", Aux::ArrayU32((&vec![1u32]).into()))
            .unwrap();
        msp.push_aux(b"al", Aux::ArrayU32((&vec![4u32]).into()))
            .unwrap();
        msp.push_aux(b"aq", Aux::ArrayU8((&vec![200u8]).into()))
            .unwrap();
        strip_consumed_legacy_tags(&mut msp);
        assert!(msp.aux(b"as").is_err(), "consumed as not stripped");
        assert!(msp.aux(b"al").is_err(), "consumed al not stripped");
        assert!(msp.aux(b"aq").is_err(), "consumed aq not stripped");

        // a length-mismatched pair fails the reader's validation, so the
        // whole legacy parse errors and NOTHING may be stripped — the
        // writer falls back to an empty model and the tags are the only
        // surviving copy of the data
        let mut mismatched = synth_record(b"ACGTACGT");
        mismatched
            .push_aux(b"as", Aux::ArrayU32((&vec![1u32, 5]).into()))
            .unwrap();
        mismatched
            .push_aux(b"al", Aux::ArrayU32((&vec![4u32]).into()))
            .unwrap();
        mismatched
            .push_aux(b"ns", Aux::ArrayU32((&vec![1u32]).into()))
            .unwrap();
        mismatched
            .push_aux(b"nl", Aux::ArrayU32((&vec![2u32]).into()))
            .unwrap();
        strip_consumed_legacy_tags(&mut mismatched);
        for tag in [b"as" as &[u8], b"al", b"ns", b"nl"] {
            assert!(
                mismatched.aux(tag).is_ok(),
                "{} stripped although the legacy parse bailed",
                String::from_utf8_lossy(tag)
            );
        }

        // a record already carrying an MA-family tag keeps everything
        let mut with_ma = synth_record(b"ACGTACGT");
        with_ma.push_aux(b"MA", Aux::String("8;")).unwrap();
        with_ma
            .push_aux(b"ns", Aux::ArrayU32((&vec![1u32]).into()))
            .unwrap();
        with_ma
            .push_aux(b"nl", Aux::ArrayU32((&vec![2u32]).into()))
            .unwrap();
        strip_consumed_legacy_tags(&mut with_ma);
        assert!(with_ma.aux(b"ns").is_ok(), "legacy tag stripped despite MA");
    }

    /// Build a synthetic record with the given forward sequence. Aligned-blocks
    /// metadata is not needed for these I/O tests because we only exercise
    /// tag-level round-trips, not coordinate liftovers.
    fn synth_record(seq: &[u8]) -> bam::Record {
        let mut record = bam::Record::new();
        let qual = vec![60u8; seq.len()];
        // unmapped record, no cigar — sufficient for tag round-trip
        record.set(b"test_read", None, seq, &qual);
        record
    }

    #[test]
    fn read_write_record_roundtrips_basemod_and_msp() {
        use crate::utils::basemods::{canonical_header, M6A_TYPE};

        // 10bp sequence with A's at positions 0, 4, 8
        let mut record = synth_record(b"ATCGATCGAT");

        // Build annotations: one m6a call at position 0 (A+a header),
        // one msp annotation 5..8 with quality 200.
        let mut annot = MolecularAnnotations::from_record(&record);
        let qspec_q = "Q".parse::<QualitySpec>().expect("Q parses");

        let m6a_type = annot.add_annotation_type(M6A_TYPE, qspec_q.clone(), Encoding::mm_ml());
        let header = canonical_header(M6A_TYPE, b'A').unwrap().to_string();
        m6a_type.add(0, 1, Strand::Forward, vec![240], Some(header));

        let msp_type = annot.add_annotation_type(MSP_TYPE, qspec_q, Encoding::Ma);
        msp_type.add(5, 3, Strand::Unknown, vec![100], None);

        write_record_with_basemods(&mut record, &annot);

        // Round-trip: read back and assert. M6A_TYPE is "a", the library's
        // internal name, so no translation is needed at this boundary.
        let back = read_record(&record).expect("read_record");

        let m6a = back
            .annotation_types
            .iter()
            .find(|t| t.name == M6A_TYPE)
            .expect("m6a type present after round-trip");
        assert_eq!(m6a.annotations.len(), 1);
        assert_eq!(m6a.annotations[0].start, 0);
        assert_eq!(m6a.annotations[0].qualities.to_vec(), vec![240]);

        let msp = back
            .annotation_types
            .iter()
            .find(|t| t.name == MSP_TYPE)
            .expect("msp type present after round-trip");
        assert_eq!(msp.annotations.len(), 1);
        assert_eq!(msp.annotations[0].start, 5);
        assert_eq!(msp.annotations[0].length, 3);
    }

    #[test]
    fn read_record_falls_back_to_legacy_ns_as() {
        use rust_htslib::bam::record::Aux;

        let mut record = synth_record(b"ATCGATCGAT");

        // Push legacy `ns`/`nl` (one nuc 2..4) and `as`/`al` (one msp 6..9).
        // No MA tag, no MM/ML.
        let ns: Vec<u32> = vec![2];
        let nl: Vec<u32> = vec![2];
        let as_starts: Vec<u32> = vec![6];
        let al: Vec<u32> = vec![3];
        record.push_aux(b"ns", Aux::ArrayU32((&ns).into())).unwrap();
        record.push_aux(b"nl", Aux::ArrayU32((&nl).into())).unwrap();
        record
            .push_aux(b"as", Aux::ArrayU32((&as_starts).into()))
            .unwrap();
        record.push_aux(b"al", Aux::ArrayU32((&al).into())).unwrap();

        let annot = read_record(&record).expect("read_record");

        let nuc = annot
            .annotation_types
            .iter()
            .find(|t| t.name == NUC_TYPE)
            .expect("nuc populated from legacy ns/nl");
        assert_eq!(nuc.annotations.len(), 1);
        assert_eq!(nuc.annotations[0].start, 2);
        assert_eq!(nuc.annotations[0].length, 2);

        let msp = annot
            .annotation_types
            .iter()
            .find(|t| t.name == MSP_TYPE)
            .expect("msp populated from legacy as/al");
        assert_eq!(msp.annotations.len(), 1);
        assert_eq!(msp.annotations[0].start, 6);
        assert_eq!(msp.annotations[0].length, 3);
    }

    #[test]
    fn read_record_splits_legacy_aq_into_msp_and_fire() {
        use rust_htslib::bam::record::Aux;

        let mut record = synth_record(b"ATCGATCGATCGATCGATCG"); // 20bp

        // Two MSPs via legacy `as`/`al`; `aq` gives the per-MSP precision that
        // legacy FIRE wrote. First MSP has precision 200 (>0 => a fire call);
        // second has precision 0 (=> an MSP that is not a regulatory element).
        let as_starts: Vec<u32> = vec![2, 10];
        let al: Vec<u32> = vec![3, 4];
        let aq: Vec<u8> = vec![200, 0];
        record
            .push_aux(b"as", Aux::ArrayU32((&as_starts).into()))
            .unwrap();
        record.push_aux(b"al", Aux::ArrayU32((&al).into())).unwrap();
        record.push_aux(b"aq", Aux::ArrayU8((&aq).into())).unwrap();

        let annot = read_record(&record).expect("read_record");

        // msp: both entries, and NO quality — the precision moved to `fire`.
        let msp = annot
            .annotation_types
            .iter()
            .find(|t| t.name == MSP_TYPE)
            .expect("msp populated from legacy as/al");
        assert!(
            !msp.quality_spec.has_quality(),
            "msp must carry no quality after legacy ingest"
        );
        assert_eq!(msp.annotations.len(), 2);
        assert_eq!(msp.annotations[0].start, 2);
        assert_eq!(msp.annotations[1].start, 10);

        // fire: only the aq>0 subset, carrying aq as the precision quality.
        let fire = annot
            .annotation_types
            .iter()
            .find(|t| t.name == FIRE_TYPE)
            .expect("fire synthesized from the aq>0 MSP subset");
        assert!(fire.quality_spec.has_quality());
        assert_eq!(fire.annotations.len(), 1);
        assert_eq!(fire.annotations[0].start, 2);
        assert_eq!(fire.annotations[0].length, 3);
        assert_eq!(fire.annotations[0].qualities.to_vec(), vec![200]);
    }

    #[test]
    fn read_record_falls_back_to_legacy_fibertig_fs_fl_fa() {
        use crate::utils::fibertig::FIBERTIG_TYPE;
        use rust_htslib::bam::record::Aux;

        let mut record = synth_record(b"ATCGATCGATCGATCGATCG"); // 20bp

        // Two legacy fibertig annotations via `fs`/`fl`; `fa` names the first
        // and leaves the second unnamed (empty `|` segment). No MA tag.
        let fs: Vec<u32> = vec![2, 10];
        let fl: Vec<u32> = vec![3, 4];
        record.push_aux(b"fs", Aux::ArrayU32((&fs).into())).unwrap();
        record.push_aux(b"fl", Aux::ArrayU32((&fl).into())).unwrap();
        record.push_aux(b"fa", Aux::String("gene_a|")).unwrap();

        let annot = read_record(&record).expect("read_record");

        let tig = annot
            .annotation_types
            .iter()
            .find(|t| t.name == FIBERTIG_TYPE)
            .expect("fibertig populated from legacy fs/fl/fa");
        assert_eq!(tig.annotations.len(), 2);
        assert_eq!(tig.annotations[0].start, 2);
        assert_eq!(tig.annotations[0].length, 3);
        assert_eq!(tig.annotations[0].name.as_deref(), Some("gene_a"));
        assert_eq!(tig.annotations[1].start, 10);
        assert_eq!(tig.annotations[1].length, 4);
        assert_eq!(tig.annotations[1].name, None);
    }

    #[test]
    fn read_record_rejects_mismatched_legacy_fibertig_lengths() {
        use rust_htslib::bam::record::Aux;

        let mut record = synth_record(b"ATCGATCGAT");

        // `fs` has two entries but `fl` only one — a corrupt legacy record.
        let fs: Vec<u32> = vec![2, 6];
        let fl: Vec<u32> = vec![3];
        record.push_aux(b"fs", Aux::ArrayU32((&fs).into())).unwrap();
        record.push_aux(b"fl", Aux::ArrayU32((&fl).into())).unwrap();

        assert!(
            read_record(&record).is_err(),
            "mismatched fs/fl lengths must be a hard error, not a silent drop"
        );
    }

    #[test]
    fn rewrite_replaces_ma_tag_instead_of_appending() {
        // 300 bp so the msps below fit SEQ; read_record drops annotations
        // that run past the sequence (#136).
        let mut record = synth_record(&b"ATCGATCGAT".repeat(30));
        let qspec_q = "Q".parse::<QualitySpec>().unwrap();

        // First write: one msp at 100..150.
        let mut v1 = MolecularAnnotations::from_record(&record);
        v1.add_annotation_type(MSP_TYPE, qspec_q.clone(), Encoding::Ma)
            .add(100, 50, Strand::Unknown, vec![10], None);
        write_record(&mut record, &v1);

        // Second write to the SAME record: a different msp at 200..260.
        let mut v2 = MolecularAnnotations::from_record(&record);
        v2.annotation_types.clear();
        v2.add_annotation_type(MSP_TYPE, qspec_q, Encoding::Ma).add(
            200,
            60,
            Strand::Unknown,
            vec![20],
            None,
        );
        write_record(&mut record, &v2);

        // Structural invariant: exactly one MA tag on the record.
        let ma_count = record
            .aux_iter()
            .filter_map(Result::ok)
            .filter(|(tag, _)| *tag == b"Ma")
            .count();
        assert_eq!(ma_count, 1, "expected exactly one MA tag, found {ma_count}");

        // Behavioral invariant: we read back v2, not the stale v1.
        let back = read_record(&record).expect("read_record");
        let msp = back
            .annotation_types
            .iter()
            .find(|t| t.name == MSP_TYPE)
            .expect("msp present");
        assert_eq!(msp.annotations.len(), 1);
        assert_eq!(msp.annotations[0].start, 200, "read back stale MA tag");
    }

    /// A mapped synthetic record with the given SEQ, CIGAR and flags.
    fn synth_aligned(seq: &[u8], cigar: &str, flags: u16) -> bam::Record {
        use rust_htslib::bam::record::CigarString;
        let mut record = bam::Record::new();
        let qual = vec![60u8; seq.len()];
        let cigar = CigarString::try_from(cigar).expect("cigar parses");
        record.set(b"frame_test", Some(&cigar), seq, &qual);
        record.set_flags(flags);
        record.set_tid(0);
        record.set_pos(0);
        record
    }

    fn legacy(record: &mut bam::Record, starts: &[u32], lens: &[u32]) {
        record
            .push_aux(b"ns", Aux::ArrayU32(starts.into()))
            .unwrap();
        record.push_aux(b"nl", Aux::ArrayU32(lens.into())).unwrap();
    }

    fn nuc_starts(record: &bam::Record) -> Vec<u32> {
        let annot = read_record(record).expect("read_record");
        annot
            .get_type(NUC_TYPE)
            .map(|t| t.annotations.iter().map(|a| a.start).collect())
            .unwrap_or_default()
    }

    // MA carries its own frame: a mismatch is stale on either strand, a
    // match is fine even with hard clips (tags computed after clipping).
    #[test]
    fn stale_frame_ma_read_length() {
        let seq = b"ACGT".repeat(50); // 200 bp
        for flags in [0u16, 16] {
            let mut r = synth_aligned(&seq, "200M", flags);
            r.push_aux(b"Ma", Aux::String("300;nuc.:10-40")).unwrap();
            assert!(record_frame_reason(&r).is_some(), "flags {flags}");
            assert!(nuc_starts(&r).is_empty());
            assert_eq!(
                read_record(&r).unwrap().read_length,
                200,
                "frame reset to SEQ"
            );
        }
        let mut r = synth_aligned(&seq, "50H200M", 2048);
        r.push_aux(b"Ma", Aux::String("200;nuc.:10-40")).unwrap();
        assert!(record_frame_reason(&r).is_none());
        assert_eq!(nuc_starts(&r), vec![9], "MA text is 1-based");
    }

    // Legacy tags carry no frame and were always written on the full read:
    // a hard clip puts them in the full-read frame, a soft clip does not,
    // and no SEQ is stale.
    #[test]
    fn stale_frame_legacy_uses_hard_clips() {
        let seq = b"ACGT".repeat(50);
        let mut fits = synth_aligned(&seq, "30H200M", 2048);
        legacy(&mut fits, &[10, 100], &[20, 20]);
        assert_eq!(
            record_frame(&fits),
            Frame::FullRead {
                h_lead: 30,
                read_length: 230
            }
        );
        assert!(record_frame_reason(&fits).is_none());
        assert_eq!(
            nuc_starts(&fits),
            vec![10, 100],
            "legacy coords are full-read molecular coords"
        );
        assert_eq!(read_record(&fits).unwrap().read_length, 230);

        let mut soft = synth_aligned(&seq, "30S170M", 2048);
        legacy(&mut soft, &[10, 100], &[20, 20]);
        assert!(record_frame_reason(&soft).is_none());
        assert_eq!(nuc_starts(&soft), vec![10, 100]);

        let mut exact = synth_aligned(&seq, "200M", 0);
        legacy(&mut exact, &[180], &[20]);
        assert!(stale_frame_reason(&read_record(&exact).unwrap(), &exact).is_none());
        assert_eq!(nuc_starts(&exact), vec![180]);

        let mut past = synth_aligned(&seq, "200M", 0);
        legacy(&mut past, &[180], &[21]);
        assert!(nuc_starts(&past).is_empty());

        let mut seqless = synth_aligned(b"", "200M", 256);
        legacy(&mut seqless, &[10], &[20]);
        assert!(record_frame_reason(&seqless).is_some());
        let annot = read_record(&seqless).unwrap();
        assert!(annot.annotation_types.is_empty());
        assert_eq!(
            annot.read_length, 200,
            "frame from the CIGAR when SEQ is absent"
        );
    }

    // MM/ML: MN is the frame when present, hard clips otherwise.
    #[test]
    fn stale_frame_mm_ml_uses_mn_then_hard_clips() {
        let seq = b"ACGT".repeat(50);
        let mm = |r: &mut bam::Record| {
            r.push_aux(b"MM", Aux::String("A+a.,0;")).unwrap();
            r.push_aux(b"ML", Aux::ArrayU8((&[200u8][..]).into()))
                .unwrap();
        };
        let mut mn_mismatch = synth_aligned(&seq, "200M", 0);
        mm(&mut mn_mismatch);
        mn_mismatch.push_aux(b"MN", Aux::I32(500)).unwrap();
        assert!(record_frame_reason(&mn_mismatch).is_some());
        assert!(read_record(&mn_mismatch)
            .unwrap()
            .annotation_types
            .is_empty());

        let mut mn_ok_clipped = synth_aligned(&seq, "50H200M", 2048);
        mm(&mut mn_ok_clipped);
        mn_ok_clipped.push_aux(b"MN", Aux::I32(200)).unwrap();
        assert!(record_frame_reason(&mn_ok_clipped).is_none());
        assert!(!read_record(&mn_ok_clipped)
            .unwrap()
            .annotation_types
            .is_empty());

        let mut no_mn_clipped = synth_aligned(&seq, "50H200M", 2048);
        mm(&mut no_mn_clipped);
        assert!(record_frame_reason(&no_mn_clipped).is_some());

        let mut no_mn_soft = synth_aligned(&seq, "50S150M", 2048);
        mm(&mut no_mn_soft);
        assert!(record_frame_reason(&no_mn_soft).is_none());

        // An MA tag written for this SEQ vouches for the MM/ML next to it:
        // that is the shape fibertools' own writers produce.
        let mut ma_vouches = synth_aligned(&seq, "50H200M", 2048);
        mm(&mut ma_vouches);
        ma_vouches
            .push_aux(b"Ma", Aux::String("200;nuc.:10-40"))
            .unwrap();
        assert!(record_frame_reason(&ma_vouches).is_none());
        let annot = read_record(&ma_vouches).unwrap();
        assert!(annot.get_type(NUC_TYPE).is_some());
        assert!(annot.get_type(crate::utils::basemods::M6A_TYPE).is_some());

        let mut seqless = synth_aligned(b"", "200M", 256);
        mm(&mut seqless);
        assert!(record_frame_reason(&seqless).is_some());
    }

    // fibertools' own base-mod writer records the frame (MN) so its output
    // survives a later hard clip check, including on a hard-clipped record.
    #[test]
    fn write_record_with_basemods_writes_mn_and_reads_back() {
        use crate::utils::basemods::{canonical_header, M6A_TYPE};
        let seq = b"ACGT".repeat(50);
        let mut r = synth_aligned(&seq, "50H200M", 2048);
        let mut annot = read_record(&r).unwrap();
        let qspec = "Q".parse::<QualitySpec>().unwrap();
        let header = canonical_header(M6A_TYPE, b'A').unwrap().to_string();
        annot
            .add_annotation_type(M6A_TYPE, qspec, Encoding::mm_ml())
            .add(0, 1, Strand::Forward, vec![200], Some(header));
        write_record_with_basemods(&mut r, &annot);
        assert!(
            matches!(r.aux(b"MN"), Ok(Aux::I32(200))),
            "MN = {:?}",
            r.aux(b"MN")
        );
        assert!(record_frame_reason(&r).is_none());
        let back = read_record(&r).unwrap();
        assert!(back.get_type(M6A_TYPE).is_some(), "m6a lost on re-read");
    }

    // A stale record leaves the writer as an honest untagged read: no legacy
    // arrays, no MM/ML/MN, and an MA tag whose frame is SEQ.
    #[test]
    fn write_record_strips_every_tag_of_a_stale_record() {
        let seq = b"ACGT".repeat(50);
        let mut r = synth_aligned(&seq, "200M", 0);
        legacy(&mut r, &[10], &[20]);
        r.push_aux(b"Ma", Aux::String("300;nuc.:11-20")).unwrap();
        r.push_aux(b"MM", Aux::String("A+a.,0;")).unwrap();
        r.push_aux(b"ML", Aux::ArrayU8((&[200u8][..]).into()))
            .unwrap();
        let annot = read_record(&r).unwrap();
        assert!(annot.annotation_types.is_empty());
        write_record(&mut r, &annot);
        for tag in [b"ns", b"nl", b"MM", b"ML", b"MN"] {
            assert!(
                r.aux(tag).is_err(),
                "{} survived",
                String::from_utf8_lossy(tag)
            );
        }
        let ma = r.aux(b"Ma");
        assert!(
            matches!(ma, Ok(Aux::String(s)) if s.split(';').next() == Some("200")),
            "Ma = {ma:?}"
        );
    }

    // MA read length decides the frame three ways; the offset is the leading
    // hard clip on both strands.
    #[test]
    fn frame_rule_ma_three_way() {
        let seq = b"ACGT".repeat(50); // 200 bp
        let mut full = synth_aligned(&seq, "50H200M50H", 2048);
        // molecular [60,100), [200,240)
        full.push_aux(b"Ma", Aux::String("300;nuc.:61-40,201-40"))
            .unwrap();
        assert_eq!(
            record_frame(&full),
            Frame::FullRead {
                h_lead: 50,
                read_length: 300
            }
        );
        assert!(record_frame_reason(&full).is_none());
        let annot = read_record(&full).unwrap();
        assert_eq!(annot.read_length, 300);
        assert_eq!(annot.query_offset(), 50);
        assert_eq!(nuc_starts(&full), vec![60, 200], "tag kept as is");
        assert_eq!(
            annot.get_ref_coords(NUC_TYPE).unwrap(),
            vec![
                (60, 100, Some(10), Some(50)),
                (200, 240, Some(150), Some(190))
            ]
        );

        let mut clipped = synth_aligned(&seq, "50H200M50H", 2048);
        clipped
            .push_aux(b"Ma", Aux::String("200;nuc.:11-40"))
            .unwrap();
        assert_eq!(record_frame(&clipped), Frame::Seq);
        let annot = read_record(&clipped).unwrap();
        assert_eq!((annot.read_length, annot.query_offset()), (200, 0));
        assert_eq!(
            annot.get_ref_coords(NUC_TYPE).unwrap(),
            vec![(10, 50, Some(10), Some(50))],
            "clipped frame lifts as today"
        );

        let mut stale = synth_aligned(&seq, "50H200M50H", 2048);
        stale
            .push_aux(b"Ma", Aux::String("250;nuc.:11-40"))
            .unwrap();
        assert!(matches!(record_frame(&stale), Frame::Stale(_)));
        assert!(nuc_starts(&stale).is_empty());
        assert_eq!(read_record(&stale).unwrap().read_length, 200);

        // reverse strand, trailing clip: offset 0; molecular [200,240) -> BAM [60,100)
        let mut rev = synth_aligned(&seq, "200M100H", 2064);
        rev.push_aux(b"Ma", Aux::String("300;nuc.:201-40")).unwrap();
        let annot = read_record(&rev).unwrap();
        assert_eq!(annot.query_offset(), 0);
        assert_eq!(
            annot.get_ref_coords(NUC_TYPE).unwrap(),
            vec![(60, 100, Some(60), Some(100))]
        );
        // reverse strand, leading clip: molecular [100,140) -> BAM [160,200) -> SEQ [60,100)
        let mut rev_lead = synth_aligned(&seq, "100H200M", 2064);
        rev_lead
            .push_aux(b"Ma", Aux::String("300;nuc.:101-40"))
            .unwrap();
        let annot = read_record(&rev_lead).unwrap();
        assert_eq!(annot.query_offset(), 100);
        assert_eq!(
            annot.get_ref_coords(NUC_TYPE).unwrap(),
            vec![(160, 200, Some(60), Some(100))]
        );
        // an annotation past the full read length is still stale
        let mut past = synth_aligned(&seq, "50H200M50H", 2048);
        past.push_aux(b"Ma", Aux::String("300;nuc.:281-40"))
            .unwrap();
        assert!(nuc_starts(&past).is_empty());
    }

    // A full-read frame keeps nuc/msp, drops m6A, and the writer strips
    // MM/ML/MN.
    #[test]
    fn full_frame_drops_mm_ml_and_strips_on_write() {
        use crate::utils::basemods::M6A_TYPE;
        let seq = b"ACGT".repeat(50);
        // MN of the full read (dorado copy) or of SEQ: dropped either way
        for mn in [250i32, 200] {
            let mut r = synth_aligned(&seq, "50H200M", 2048);
            r.push_aux(b"Ma", Aux::String("250;nuc.:61-40;msp.:101-20"))
                .unwrap();
            r.push_aux(b"MM", Aux::String("A+a.,0;")).unwrap();
            r.push_aux(b"ML", Aux::ArrayU8((&[200u8][..]).into()))
                .unwrap();
            r.push_aux(b"MN", Aux::I32(mn)).unwrap();
            let annot = read_record(&r).unwrap();
            assert!(annot.get_type(NUC_TYPE).is_some() && annot.get_type(MSP_TYPE).is_some());
            assert!(
                annot.get_type(M6A_TYPE).is_none(),
                "MN {mn}: m6A must be dropped on a full-read frame"
            );
            write_record(&mut r, &annot);
            for tag in [b"MM", b"ML", b"MN"] {
                assert!(
                    r.aux(tag).is_err(),
                    "MN {mn}: {} survived",
                    String::from_utf8_lossy(tag)
                );
            }
            let ma = r.aux(b"Ma");
            assert!(
                matches!(ma, Ok(Aux::String(s)) if s.starts_with("250;") && s.contains("nuc") && s.contains("msp")),
                "Ma = {ma:?}"
            );
            // write_record_with_basemods on the same model writes no MM/MN either
            write_record_with_basemods(&mut r, &annot);
            assert!(r.aux(b"MM").is_err() && r.aux(b"MN").is_err());
            // and the rewritten record reads back in the same frame
            let back = read_record(&r).unwrap();
            assert_eq!((back.read_length, back.query_offset()), (250, 50));
            assert_eq!(back.get_forward_coords(NUC_TYPE), Some(vec![(60, 100)]));
        }
    }

    // Legacy ns/nl on a hard clip: read_length = seq + H_lead + H_trail,
    // full-read frame.
    #[test]
    fn legacy_on_hard_clip_is_a_full_read_frame() {
        let seq = b"ACGT".repeat(50);
        let mut r = synth_aligned(&seq, "30H200M", 2048);
        // [10,30) before SEQ, [20,60) straddles, [100,120) inside
        legacy(&mut r, &[10, 20, 100], &[20, 40, 20]);
        let annot = read_record(&r).unwrap();
        assert_eq!((annot.read_length, annot.query_offset()), (230, 30));
        assert_eq!(
            annot.get_ref_coords(NUC_TYPE).unwrap(),
            vec![
                (10, 30, None, None),
                (20, 60, Some(0), Some(30)),
                (100, 120, Some(70), Some(90))
            ]
        );
        write_record(&mut r, &annot);
        assert!(
            r.aux(b"ns").is_err() && r.aux(b"nl").is_err(),
            "consumed legacy tags stripped"
        );
        assert!(
            matches!(r.aux(b"Ma"), Ok(Aux::String(s)) if s.starts_with("230;") && s.contains("nuc.:11-20,21-40,101-20")),
            "Ma = {:?}",
            r.aux(b"Ma")
        );
        let mut trail = synth_aligned(&seq, "200M30H", 0);
        legacy(&mut trail, &[100], &[20]);
        let annot = read_record(&trail).unwrap();
        assert_eq!((annot.read_length, annot.query_offset()), (230, 0));
        assert_eq!(
            annot.get_ref_coords(NUC_TYPE).unwrap(),
            vec![(100, 120, Some(100), Some(120))]
        );
        // molecular [100,120) -> BAM [110,130)
        let mut rev = synth_aligned(&seq, "200M30H", 2064);
        legacy(&mut rev, &[100], &[20]);
        assert_eq!(
            read_record(&rev).unwrap().get_ref_coords(NUC_TYPE).unwrap(),
            vec![(110, 130, Some(110), Some(130))]
        );
    }

    // No m6A on a full-read frame means NotCallable, whatever the on-disk
    // span says; an unprocessed full-frame record stays Untagged; a clipped
    // frame keeps its span.
    #[test]
    fn full_frame_records_derive_not_callable() {
        use crate::fiber::{CallableState, FiberseqData};
        let seq = b"ACGT".repeat(50);
        let filters = crate::utils::input_bam::FiberFilters::default();
        // 10 MSPs, mean 15, end 150
        let msp = (0..10)
            .map(|i| format!("{}-15", 1 + i * 15))
            .collect::<Vec<_>>()
            .join(",");
        let state = |r: &bam::Record| {
            FiberseqData::new(r.clone(), None, &filters)
                .callable_state()
                .0
        };
        let marker = |r: &bam::Record| {
            let mut annot = read_record(r).unwrap();
            sync_fiberseq_callable(&mut annot, r, &filters);
            annot
                .get_type(FIBERSEQ_CALLABLE_TYPE)
                .map(|t| (t.annotations[0].start, t.annotations[0].length))
        };
        let mut full = synth_aligned(&seq, "50H200M", 2048);
        full.push_aux(
            b"Ma",
            Aux::String(&format!("250;msp.:{msp};fiberseq_callable.:1-150")),
        )
        .unwrap();
        assert_eq!(
            marker(&full),
            Some((0, 0)),
            "no m6A: NotCallable marker replaces the span"
        );
        assert_eq!(state(&full), CallableState::NotCallable);
        let mut clipped = synth_aligned(&seq, "50H200M", 2048);
        clipped
            .push_aux(
                b"Ma",
                Aux::String(&format!("200;msp.:{msp};fiberseq_callable.:1-150")),
            )
            .unwrap();
        assert_eq!(marker(&clipped), Some((0, 150)));
        assert_eq!(state(&clipped), CallableState::Callable);
        let mut bare = synth_aligned(&seq, "50H200M", 2048);
        bare.push_aux(b"Ma", Aux::String("250")).unwrap();
        assert_eq!(marker(&bare), None);
        assert_eq!(state(&bare), CallableState::Untagged);
    }

    // A SEQ-less hard-clipped record with a full-read MA tag is judged by
    // its CIGAR span, so it gets the offset lift like a record with SEQ.
    #[test]
    fn seqless_hard_clipped_ma_is_full_read_frame() {
        let mut r = synth_aligned(b"", "50H200M", 2048);
        r.push_aux(b"Ma", Aux::String("250;nuc.:60-40")).unwrap();
        assert_eq!(
            record_frame(&r),
            Frame::FullRead {
                h_lead: 50,
                read_length: 250
            }
        );
        let annot = read_record(&r).unwrap();
        assert!(annot.get_type(NUC_TYPE).is_some());
        assert!(model_is_full_read_frame(&annot, &r));
    }

    // An MN tag that disagrees with SEQ next to an MA tag that matches it:
    // the calls stay, only the base mods go, on read and on write.
    #[test]
    fn mn_disagreement_drops_only_base_mods() {
        let seq = b"ACGT".repeat(50);
        let mut r = synth_aligned(&seq, "200M", 0);
        r.push_aux(b"Ma", Aux::String("200;nuc.:10-40")).unwrap();
        r.push_aux(b"MM", Aux::String("A+a.,0;")).unwrap();
        r.push_aux(b"ML", Aux::ArrayU8((&[200u8][..]).into()))
            .unwrap();
        r.push_aux(b"MN", Aux::I32(500)).unwrap();
        assert_eq!(record_frame(&r), Frame::Seq);
        let annot = read_record(&r).unwrap();
        assert!(annot.get_type(NUC_TYPE).is_some());
        assert!(annot.annotation_types.iter().all(|t| !t.is_mm_ml()));
        write_record(&mut r, &annot);
        for tag in [b"MM", b"ML", b"MN"] {
            assert!(
                r.aux(tag).is_err(),
                "{} survived",
                String::from_utf8_lossy(tag)
            );
        }
        assert!(r.aux(b"Ma").is_ok());
    }

    // A Callable span copied onto a full-frame record with no calls of its
    // own is replaced by the NotCallable marker rather than kept.
    #[test]
    fn full_frame_copied_callable_span_is_replaced() {
        let seq = b"ACGT".repeat(50);
        let filters = crate::utils::input_bam::FiberFilters::default();
        let mut r = synth_aligned(&seq, "50H200M", 2048);
        r.push_aux(b"Ma", Aux::String("250;fiberseq_callable.:1-200"))
            .unwrap();
        let mut annot = read_record(&r).unwrap();
        sync_fiberseq_callable(&mut annot, &r, &filters);
        let t = annot.get_type(FIBERSEQ_CALLABLE_TYPE).expect("marker kept");
        assert_eq!(t.annotations[0].length, 0, "copied Callable span survived");
    }
}
