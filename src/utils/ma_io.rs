//! MA-spec annotation I/O bridge for fibertools-rs.
//!
//! Reading: prefers MA/AL/AQ/AN; falls back to legacy `ns/nl/as/al/aq` if no
//! MA tag is present. Returns a populated [`MolecularAnnotations`] (read
//! length and aligned blocks set) even when no annotation tags exist.
//!
//! Writing: always emits MA-spec tags. With `legacy=true` it additionally
//! emits legacy `ns/nl/as/al/aq` for `nuc` and `msp` types so pre-MA
//! consumers (older pyft, IGV decorators) keep working during the
//! migration. With `legacy=false` any pre-existing legacy tags are stripped.
//!
//! Annotation type names produced by fibertools-rs:
//! - `nuc`  (forward strand, no quality)
//! - `msp`  (forward strand, no quality pre-FIRE; `Q` post-FIRE)
//! - `fire` (forward strand, `P` phred quality)
//!
//! `m6a` and `cpg` types may appear *in memory* on a [`MolecularAnnotations`]
//! populated by `BaseMods::populate_ma`, but are NEVER read or written here:
//! their on-disk source of truth is `MM`/`ML`.

use anyhow::{bail, Result};
use molecular_annotation::{Annotation, MolecularAnnotations, QualitySpec, Strand};
use rust_htslib::bam::{self, record::Aux};

/// Annotation type names used by fibertools-rs.
pub const NUC_TYPE: &str = "nuc";
pub const MSP_TYPE: &str = "msp";

/// Legacy fibertools-rs tags. Listed here so callers (and the writer below)
/// have one place to look up which tags this module produces/consumes.
pub const LEGACY_NUC_MSP_TAGS: &[&[u8]] = &[b"ns", b"nl", b"as", b"al", b"aq"];

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
    let ma = match record.aux(b"MA") {
        Ok(Aux::String(s)) => s.to_string(),
        _ => return Ok(None),
    };
    let aq: Option<Vec<u8>> = match record.aux(b"AQ") {
        Ok(Aux::ArrayU8(arr)) => Some(arr.iter().collect()),
        _ => None,
    };
    let an: Option<String> = match record.aux(b"AN") {
        Ok(Aux::String(s)) => Some(s.to_string()),
        _ => None,
    };
    let al: Vec<u32> = match record.aux(b"AL") {
        Ok(Aux::ArrayU32(arr)) => arr.iter().collect(),
        Ok(Aux::ArrayI32(arr)) => arr.iter().map(|v| v as u32).collect(),
        _ => Vec::new(),
    };
    let mut annot = MolecularAnnotations::from_tags(&ma, &al, aq.as_deref(), an.as_deref())
        .map_err(|e| anyhow::anyhow!("MA tag parse error: {e}"))?;
    annot.set_aligned_blocks_raw(
        molecular_annotation::AlignedBlocks::from_record(record),
        record.is_reverse(),
    );
    Ok(Some(annot))
}

fn read_legacy_nuc_msp(record: &bam::Record) -> Result<MolecularAnnotations> {
    let mut annot = MolecularAnnotations::from_record(record);

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
            let nuc = annot.add_annotation_type(NUC_TYPE, QualitySpec::none());
            for (s, l) in ns.iter().zip(nl.iter()) {
                nuc.add(*s, *l, Strand::Forward, vec![], None);
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
        let q_spec = if aq.is_some() {
            "Q".parse::<QualitySpec>()?
        } else {
            QualitySpec::none()
        };
        if let Some(ref q) = aq {
            if q.len() != starts.len() {
                bail!(
                    "legacy aq ({}) and as ({}) length mismatch",
                    q.len(),
                    starts.len()
                );
            }
        }
        if !starts.is_empty() {
            let msp = annot.add_annotation_type(MSP_TYPE, q_spec);
            for (i, (s, l)) in starts.iter().zip(lens.iter()).enumerate() {
                let qualities = aq.as_ref().map(|q| vec![q[i]]).unwrap_or_default();
                msp.add(*s, *l, Strand::Forward, qualities, None);
            }
        }
    }

    Ok(annot)
}

/// Write annotations to a BAM record.
///
/// Always emits the MA-spec tags. When `legacy` is set, also emits the
/// fibertools-rs legacy `ns`/`nl`/`as`/`al`/`aq` tags for `nuc` and `msp`
/// annotation types. When not set, any pre-existing legacy tags are stripped
/// so the output is MA-only.
pub fn write_annotations(record: &mut bam::Record, annot: &MolecularAnnotations, legacy: bool) {
    for tag in LEGACY_NUC_MSP_TAGS {
        record.remove_aux(tag).ok();
    }
    annot.to_record(record);
    if legacy {
        write_legacy_nuc_msp(record, annot);
    }
}

fn write_legacy_nuc_msp(record: &mut bam::Record, annot: &MolecularAnnotations) {
    if let Some(nuc) = annot.get_type(NUC_TYPE) {
        push_starts_lens(record, b"ns", b"nl", &nuc.annotations);
    }
    if let Some(msp) = annot.get_type(MSP_TYPE) {
        push_starts_lens(record, b"as", b"al", &msp.annotations);
        if msp.quality_spec.num_qualities() == 1 && msp.quality_spec.has_quality() {
            let aq: Vec<u8> = msp
                .annotations
                .iter()
                .map(|a| a.qualities.first().copied().unwrap_or(0))
                .collect();
            if !aq.is_empty() {
                record.push_aux(b"aq", Aux::ArrayU8((&aq).into())).ok();
            }
        }
    }
}

fn push_starts_lens(
    record: &mut bam::Record,
    starts_tag: &[u8],
    lens_tag: &[u8],
    items: &[Annotation],
) {
    if items.is_empty() {
        return;
    }
    let starts: Vec<u32> = items.iter().map(|a| a.start).collect();
    let lens: Vec<u32> = items.iter().map(|a| a.length).collect();
    record
        .push_aux(starts_tag, Aux::ArrayU32((&starts).into()))
        .ok();
    record
        .push_aux(lens_tag, Aux::ArrayU32((&lens).into()))
        .ok();
}

/// Convenience for callers that still operate on raw `i64` arrays of
/// nucleosome and MSP coordinates. Returns
/// `(nuc_starts, nuc_lengths, msp_starts, msp_lengths, msp_qual)` derived
/// from whichever tag set is present (MA wins over legacy).
///
/// Coordinates are in molecular orientation — the same convention as the
/// legacy `ns`/`nl`/`as`/`al` tags. `msp_qual` is empty when there are no
/// MSP qualities (pre-FIRE state); otherwise one byte per MSP.
pub fn extract_nuc_msp_arrays(
    record: &bam::Record,
) -> Result<(Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>, Vec<u8>)> {
    let annot = read_annotations(record)?;
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
                    .map(|a| a.qualities.first().copied().unwrap_or(0))
                    .collect()
            } else {
                Vec::new()
            };
            (starts, lens, qs)
        })
        .unwrap_or_default();
    Ok((nuc_starts, nuc_lengths, msp_starts, msp_lengths, msp_qual))
}

/// Build a [`MolecularAnnotations`] for `nuc` and `msp` annotation types from
/// raw `u32` start/length arrays in molecular orientation, plus an optional
/// per-MSP linear quality array. Convenience for producers that have already
/// computed coordinates and want to hand them to [`write_annotations`].
pub fn build_nuc_msp_annotations(
    record: &bam::Record,
    nuc_starts: &[u32],
    nuc_lengths: &[u32],
    msp_starts: &[u32],
    msp_lengths: &[u32],
    msp_qual: Option<&[u8]>,
) -> MolecularAnnotations {
    let mut annot = MolecularAnnotations::from_record(record);
    if !nuc_starts.is_empty() {
        let nuc = annot.add_annotation_type(NUC_TYPE, QualitySpec::none());
        for (s, l) in nuc_starts.iter().zip(nuc_lengths.iter()) {
            nuc.add(*s, *l, Strand::Forward, vec![], None);
        }
    }
    if !msp_starts.is_empty() {
        let q_spec = match msp_qual {
            Some(_) => "Q".parse::<QualitySpec>().expect("Q parses"),
            None => QualitySpec::none(),
        };
        let msp = annot.add_annotation_type(MSP_TYPE, q_spec);
        for (i, (s, l)) in msp_starts.iter().zip(msp_lengths.iter()).enumerate() {
            let q = msp_qual.map(|q| vec![q[i]]).unwrap_or_default();
            msp.add(*s, *l, Strand::Forward, q, None);
        }
    }
    annot
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
