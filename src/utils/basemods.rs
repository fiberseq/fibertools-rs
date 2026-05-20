//! MM/ML BAM tag I/O for fiberseq base modifications.
//!
//! Scope: m6A (`"a"`) and 5mC (`"m"`) only. MM groups carrying any other
//! mod-type code are logged and skipped on parse — fibertools-rs's typed
//! API only exposes these two canonical fiberseq modifications. Adding a
//! new canonical mod (e.g. 5hmC, `"h"`) is a ~5-line change across
//! [`ma_type_for_mod_code`], [`is_basemod_type`], and
//! [`canonical_header`].
//!
//! Read path: `parse_mm_ml_into_ma` parses a record's MM/ML tags directly
//! into a [`MolecularAnnotations`]. Each call's MM group header (e.g.
//! `"A+a"`, `"T-a"`, `"N+a"`, `"C+m"`) is stored verbatim on
//! [`Annotation::name`], so headers round-trip through `write_mm_ml`
//! byte-for-byte — including non-standard encodings like `N+a` that the
//! old write path silently normalized.
//!
//! Write path: `write_mm_ml` emits MM/ML by reading each annotation's
//! `name` directly. Basemod annotations *must* have `name` set; the
//! parser sets it from the MM tag, and programmatic producers
//! (`predict_m6a`, `ddda_to_m6a`) set it via [`canonical_header`]. A
//! missing `name` on a basemod annotation is a bug and triggers a panic.
//!
//! MA-tag policy: `ma_io::write_annotations` *strips* `m6a` and `cpg`
//! before MA serialization — MM/ML is the on-disk source of truth, and
//! the MA tag set holds nuc/msp/fire only. This is a policy choice
//! (avoid duplication, single source of truth), not a correctness one.
//! Lifting the strip would expose basemods to MA-native consumers (e.g.
//! the `molecular_annotation` library) at the cost of ~2× disk for
//! basemod data and a "MM/ML wins on read" invariant. The matched
//! invariant on the read side is the idempotency wipe in
//! [`parse_mm_ml_into_ma`]: it discards any basemod-typed annotations
//! present in `annot` before parsing MM/ML, so an MA tag set that
//! wrongly contained basemods doesn't get merged on top of MM/ML.

use crate::utils::bio_io::*;
use bio::alphabets::dna::revcomp;
use lazy_static::lazy_static;
use molecular_annotation::{MolecularAnnotations, QualitySpec, Strand};
use regex::Regex;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use std::collections::BTreeMap;

/// Annotation type names emitted by `parse_mm_ml_into_ma`.
pub const M6A_TYPE: &str = "m6a";
pub const CPG_TYPE: &str = "cpg";

/// Map an MM mod-type code (the trailing alpha/digits in a group like
/// `A+a`) to a fibertools MA type name. Unknown codes return `None` and
/// the caller skips the group.
fn ma_type_for_mod_code(code: &str) -> Option<&'static str> {
    match code {
        "a" => Some(M6A_TYPE),
        "m" => Some(CPG_TYPE),
        _ => None,
    }
}

/// True for annotation types whose on-disk source of truth is MM/ML.
/// Used to gate idempotency wipes (parse), MA-tag emission (ma_io), and
/// MM/ML emission (write).
pub fn is_basemod_type(name: &str) -> bool {
    matches!(name, M6A_TYPE | CPG_TYPE)
}

/// Header for a canonical basemod call landing on `base` (read in
/// forward-strand orientation). Returns `&'static str` for the four
/// possible canonical headers; `None` if the (type, base) pair isn't a
/// canonical fiberseq combination. Used by programmatic producers
/// (`predict_m6a`, `ddda_to_m6a`) to label calls so `write_mm_ml` can
/// emit them under the right MM group.
pub fn canonical_header(type_name: &str, base: u8) -> Option<&'static str> {
    match (type_name, base) {
        (M6A_TYPE, b'A') => Some("A+a"),
        (M6A_TYPE, b'T') => Some("T-a"),
        (CPG_TYPE, b'C') => Some("C+m"),
        (CPG_TYPE, b'G') => Some("G+m"),
        _ => None,
    }
}

lazy_static! {
    // MM:Z:([ACGTUN][-+]([A-Za-z]+|[0-9]+)[.?]?(,[0-9]+)*;)*
    static ref MM_RE: Regex =
        Regex::new(r"((([ACGTUN])([-+])([A-Za-z]+|[0-9]+))[.?]?((,[0-9]+)*;)*)").unwrap();
}

/// Parse the MM/ML tags of `record` directly into `annot`. Positions are
/// stored in molecular (forward) orientation as 1bp annotations with a
/// single `Q` quality score per call. Each annotation's `name` carries
/// the verbatim MM group header it was parsed from (e.g. `"A+a"`,
/// `"T-a"`, `"N+a"`, `"C+h"`), so [`write_mm_ml`] round-trips it
/// unchanged.
///
/// Per-call filtering:
/// - calls with ML quality < `min_ml_score` are dropped.
/// - calls within `strip_at_ends` bp of either end of the forward
///   sequence are dropped (no-op if `strip_at_ends <= 0`).
///
/// Dispatch by MM mod_type (exact-match) controls the MA type bucket:
/// - `"a"` → [`M6A_TYPE`] (collapses `A+a`, `T-a`, `N+a`, etc.)
/// - `"m"` → [`CPG_TYPE`] (collapses `C+m`, `G+m`, etc.)
/// - any other → annotation type named verbatim after the MM group
///   header (e.g. `"C+h"` for 5hmC, `"C+ac"` for acetyl).
///
/// Idempotent: any pre-existing basemod-shaped types on `annot` are
/// dropped first (see [`is_basemod_type`]), so MM/ML always overwrites
/// in-memory state rather than appending. This guards against (a)
/// repeated calls on the same `annot`, and (b) malformed inputs that
/// wrongly emit `m6a` into the MA tag set in addition to MM/ML.
pub fn parse_mm_ml_into_ma(
    record: &bam::Record,
    annot: &mut MolecularAnnotations,
    min_ml_score: u8,
    strip_at_ends: i64,
) {
    // MM/ML is the on-disk source of truth for base modifications; any
    // basemod-shaped type that may already be in `annot` (from a prior
    // pass, an MA-tag leak in a third-party producer, or repeated calls)
    // is discarded so we don't silently append on top of stale data.
    annot.annotation_types.retain(|t| !is_basemod_type(&t.name));

    let ml_tag = get_u8_tag(record, b"ML");
    let Ok(Aux::String(mm_text)) = record.aux(b"MM") else {
        log::trace!("No MM tag found");
        return;
    };

    let forward_bases = if record.is_reverse() {
        convert_seq_uppercase(revcomp(record.seq().as_bytes()))
    } else {
        convert_seq_uppercase(record.seq().as_bytes())
    };
    let seq_len = forward_bases.len();
    let strip = strip_at_ends.max(0) as usize;
    let upper = seq_len.saturating_sub(strip);

    // Per-call accumulator: (ma_type, verbatim_header, pos, qual). The
    // header is owned because the regex capture's lifetime ends at the
    // next iteration; the ma_type is a static interned name.
    let mut calls: Vec<(&'static str, String, u32, u8)> = Vec::new();
    let mut num_mods_seen = 0usize;

    for cap in MM_RE.captures_iter(mm_text) {
        // Capture group 2 is the full `base+strand+modtype` prefix —
        // exactly the MM group header (e.g. "A+a", "T-a", "C+h", "N+12").
        let header = cap.get(2).map_or("", |m| m.as_str()).to_string();
        let mod_base = cap.get(3).map(|m| m.as_str().as_bytes()[0]).unwrap();
        let mod_type_str = cap.get(5).map_or("", |m| m.as_str()).to_string();
        let mod_dists_str = cap.get(6).map_or("", |m| m.as_str());
        let mod_dists: Vec<i64> = mod_dists_str
            .trim_end_matches(';')
            .split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().unwrap())
            .collect();

        if mod_dists.is_empty() {
            continue;
        }

        // Resolve positions in forward sequence via the MM distance
        // encoding: each value is "how many of `mod_base` to skip
        // before the next modified base."
        let mut positions: Vec<u32> = vec![0; mod_dists.len()];
        let mut cur_mod_idx = 0;
        let mut dist_from_last_mod_base = 0i64;
        for (cur_seq_idx, cur_base) in forward_bases.iter().enumerate() {
            if cur_mod_idx >= mod_dists.len() {
                break;
            }
            if (*cur_base == mod_base || mod_base == b'N')
                && dist_from_last_mod_base == mod_dists[cur_mod_idx]
            {
                positions[cur_mod_idx] = cur_seq_idx as u32;
                dist_from_last_mod_base = 0;
                cur_mod_idx += 1;
            } else if *cur_base == mod_base {
                dist_from_last_mod_base += 1;
            }
        }
        let resolved = cur_mod_idx;
        if resolved != mod_dists.len() {
            // MM tag claims more `mod_base` positions than the forward
            // sequence actually has. Warn and emit only the calls we
            // could resolve — the alternative was to panic the whole
            // worker. ML offsets must still advance by the on-disk
            // slot count (`mod_dists.len()`) so the next group reads
            // from the correct ML section.
            log::warn!(
                "MM/ML parser: {} (reverse={}) — MM group ended {} positions short of the forward sequence; truncating",
                String::from_utf8_lossy(record.qname()),
                record.is_reverse(),
                mod_dists.len() - resolved,
            );
        }

        // Pair this group's calls with ML qualities. The on-disk ML
        // section for this group has `mod_dists.len()` entries even if
        // we could only resolve `resolved` of them.
        let group_end = num_mods_seen + mod_dists.len();
        let group_quals: Vec<u8> = if group_end > ml_tag.len() {
            let needed = group_end - ml_tag.len();
            let mut existing = ml_tag[num_mods_seen..].to_vec();
            existing.extend(std::iter::repeat_n(0, needed));
            log::warn!(
                "ML tag is too short for the number of modifications found in the MM tag. Assuming an ML value of 0 after the first {group_end} modifications."
            );
            existing
        } else {
            ml_tag[num_mods_seen..group_end].to_vec()
        };
        num_mods_seen = group_end;
        // Drop ML qualities whose position couldn't be resolved, so the
        // dispatch zip below pairs only the surviving calls.
        positions.truncate(resolved);
        let group_quals = &group_quals[..resolved];

        // Canonical "a"/"m" → m6a/cpg; everything else is dropped with a
        // warning. The per-call header travels into `Annotation.name`
        // so write_mm_ml re-emits it verbatim (preserving A+a vs T-a,
        // and non-standard encodings like N+a, on round-trip).
        let Some(ma_type) = ma_type_for_mod_code(&mod_type_str) else {
            log::warn!(
                "MM/ML parser: unsupported mod-type {:?} (group {:?}); skipping {} call(s)",
                mod_type_str,
                header,
                resolved,
            );
            continue;
        };
        for (pos, &qual) in positions.iter().copied().zip(group_quals) {
            if qual < min_ml_score {
                continue;
            }
            let p = pos as usize;
            if strip > 0 && (p < strip || p >= upper) {
                continue;
            }
            calls.push((ma_type, header.clone(), pos, qual));
        }
    }

    if ml_tag.len() != num_mods_seen {
        log::warn!(
            "ML tag ({}) different number than MM tag ({}).",
            ml_tag.len(),
            num_mods_seen,
        );
    }

    // Group flat call list by MA type so each type's annotations land
    // in the correct bucket on the spec object. Within a type, calls
    // are sorted by position — even though different MM group headers
    // (e.g. A+a and T-a) interleave under m6a, each call's `name`
    // keeps its header for the write path.
    let mut by_type: BTreeMap<&'static str, Vec<(String, u32, u8)>> = BTreeMap::new();
    for (ma_type, header, pos, qual) in calls {
        by_type.entry(ma_type).or_default().push((header, pos, qual));
    }
    for (ma_type, mut group) in by_type {
        group.sort_by_key(|(_, p, _)| *p);
        add_calls(annot, ma_type, &group);
    }
}

fn add_calls(annot: &mut MolecularAnnotations, name: &str, calls: &[(String, u32, u8)]) {
    if calls.is_empty() {
        return;
    }
    let qspec = "Q".parse::<QualitySpec>().expect("Q parses");
    let t = annot.add_annotation_type(name, qspec);
    for (header, pos, qual) in calls {
        t.add(*pos, 1, Strand::Forward, vec![*qual], Some(header.clone()));
    }
}

/// Emit MM/ML tags from `annot` onto `record`. Existing `MM` and `ML`
/// aux are removed first.
///
/// Header for each call is read directly from [`Annotation::name`].
/// Every basemod annotation must carry a header — the parser sets it
/// from the MM tag and programmatic producers set it via
/// [`canonical_header`]. A `None` `name` on a basemod annotation is a
/// caller bug and panics.
///
/// Non-basemod annotation types (`msp`, `nuc`, `fire`, etc.) are
/// skipped — see [`is_basemod_type`] for the gate.
///
/// MM groups are emitted in deterministic (lexicographic header) order.
pub fn write_mm_ml(record: &mut bam::Record, annot: &MolecularAnnotations) {
    let forward_seq = if record.is_reverse() {
        revcomp(record.seq().as_bytes())
    } else {
        record.seq().as_bytes()
    };

    // header (e.g. "A+a", "T-a") → Vec<(pos, qual)>. BTreeMap gives
    // deterministic emission order.
    let mut groups: BTreeMap<String, Vec<(u32, u8)>> = BTreeMap::new();

    for t in &annot.annotation_types {
        if !is_basemod_type(&t.name) {
            continue;
        }
        for a in &t.annotations {
            let header = a.name.clone().unwrap_or_else(|| {
                panic!(
                    "basemod annotation in type {:?} at pos {} has no MM group header; \
                     producers must set Annotation.name via canonical_header()",
                    t.name, a.start,
                )
            });
            let qual = a.qualities.first().copied().unwrap_or(0);
            groups.entry(header).or_default().push((a.start, qual));
        }
    }

    let mut mm_tag = String::new();
    let mut ml_tag: Vec<u8> = Vec::new();
    for (header, mut calls) in groups {
        calls.sort_by_key(|(p, _)| *p);
        ml_tag.extend(calls.iter().map(|&(_, q)| q));
        mm_tag.push_str(&header);
        // First byte of the header is the skip base (`N` counts everything).
        let skip_base = header.as_bytes()[0];
        let mut last_pos = 0usize;
        for (pos, _) in &calls {
            let p = *pos as usize;
            let in_between = if last_pos < p {
                forward_seq[last_pos..p]
                    .iter()
                    .filter(|&&b| skip_base == b'N' || b == skip_base)
                    .count()
            } else {
                0
            };
            last_pos = p + 1;
            mm_tag.push_str(&format!(",{in_between}"));
        }
        mm_tag.push(';');
    }

    record.remove_aux(b"MM").ok();
    record.remove_aux(b"ML").ok();
    record
        .push_aux(b"MM", Aux::String(&mm_tag))
        .expect("push MM tag");
    record
        .push_aux(b"ML", Aux::ArrayU8((&ml_tag).into()))
        .expect("push ML tag");
}

