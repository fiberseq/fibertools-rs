//! MA-spec annotation I/O bridge for fibertools-rs.
//!
//! Reading: prefers MA/AL/AQ/AN; falls back to legacy `ns/nl/as/al/aq` if no
//! MA tag is present. Returns a populated [`MolecularAnnotations`] (read
//! length and aligned blocks set) even when no annotation tags exist.
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
//!
//! `m6a` and `cpg` types may appear *in memory* on a [`MolecularAnnotations`]
//! populated by the library's MM/ML parser. Their on-disk source of truth
//! is `MM`/`ML`; they carry `Encoding::MmMl` (set at construction, whether read
//! by the parser or built by a producer via `Encoding::mm_ml()`), so the
//! library serializes them into MM/ML rather than the MA tag set.

use anyhow::{bail, Result};
use molecular_annotation::{Encoding, MolecularAnnotations, QualitySpec, Strand};
use rust_htslib::bam::{self, record::Aux};

/// Annotation type names used by fibertools-rs.
pub const NUC_TYPE: &str = "nuc";
pub const MSP_TYPE: &str = "msp";
pub const FIRE_TYPE: &str = "fire";

/// Legacy fibertools-rs tags consumed by the reader as a fallback when no
/// MA tag is present. These are never emitted; we only ingest them.
pub const LEGACY_READ_TAGS: &[&[u8]] = &[b"ns", b"nl", b"as", b"al", b"aq"];

/// Reads annotations from a BAM record's tags.
///
/// Delegates to the library's combined MA-spec + MM/ML parser. If the
/// library returns an empty annotation set and the record carries legacy
/// `ns`/`nl`/`as`/`al`/`aq` tags, falls back to `read_legacy_nuc_msp`.
/// Tolerant of malformed MM/ML — the library handles those internally
/// without panicking.
pub fn read_record(record: &bam::Record) -> Result<MolecularAnnotations> {
    let mut annot = MolecularAnnotations::from_record(record);
    // If MA tag is absent, also ingest legacy nuc/msp tags. The library
    // already populates basemod types (m6a/cpg) from MM/ML, so we merge
    // legacy-derived nuc/msp into whatever the library produced rather
    // than gating on `annotation_types.is_empty()` (which would skip the
    // legacy fallback whenever MM/ML is present).
    let has_ma = matches!(record.aux(b"MA"), Ok(Aux::String(_)));
    if !has_ma {
        // Provenance caveat: presence of `ns`/`as` is treated as evidence of
        // fibertools-legacy data with no further validation. This path is
        // read-only — we never delete these tags — so a foreign tool reusing
        // those names is at worst mis-parsed here, never destroyed. Any future
        // code that *removes* legacy tags must verify provenance first.
        let has_legacy_nuc = matches!(record.aux(b"ns"), Ok(_));
        let has_legacy_msp = matches!(record.aux(b"as"), Ok(_));
        if has_legacy_nuc || has_legacy_msp {
            let legacy = read_legacy_nuc_msp(record)?;
            for t in legacy.annotation_types.into_iter() {
                if annot.get_type(&t.name).is_some() {
                    continue;
                }
                let new_t = annot.add_annotation_type(&t.name, t.quality_spec.clone(), t.encoding);
                for a in t.annotations.into_iter() {
                    new_t.add_shared(a.start, a.length, a.strand, a.qualities, a.name);
                }
            }
        }
    }
    Ok(annot)
}

/// Writes MA-family tags (MA/AL/AQ/AN) to a BAM record, **preserving the
/// record's existing MM/ML bytes**.
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
    annot.to_record(record);
}

/// Writes MA-family tags **and** canonically re-encodes MM/ML from the
/// annotation model.
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
    annot.to_record(record);
    annot.write_mm_ml(record);
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

    Ok(annot)
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
    msp: Option<(&[u32], &[u32], Option<&[u8]>)>,
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
    fn rewrite_replaces_ma_tag_instead_of_appending() {
        let mut record = synth_record(b"ATCGATCGAT");
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
            .filter(|(tag, _)| *tag == b"MA")
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
}
