//! MM/ML BAM tag I/O for fiberseq base modifications.
//!
//! Read path: `parse_mm_ml_into_ma` parses a record's MM/ML tags directly
//! into a [`MolecularAnnotations`] as annotation types. No intermediate
//! `BaseMods` struct — the spec object is the sole in-memory
//! representation.
//!
//! Mod-type dispatch (exact match on the MM mod_type code):
//! - `"a"` (any group, e.g. `A+a`, `T-a`, `N+a`) → `m6a` type.
//! - `"m"` (any group, e.g. `C+m`, `G+m`) → `cpg` type.
//! - any other code (`C+h`, `T+f`, `C+ac`, ChEBI numeric IDs, etc.) → its
//!   own annotation type named verbatim after the MM group header (e.g.
//!   `"C+h"`). Calls round-trip through `write_mm_ml` unchanged.
//!
//! Write path: `write_mm_ml` emits MM/ML from a [`MolecularAnnotations`].
//! Canonical types (`m6a`, `cpg`) have their groups reconstructed from
//! the forward sequence (`A`→`A+a`, `T`→`T-a`, `C`→`C+m`, `G`→`G+m`).
//! Non-canonical types stored under a valid MM group name are emitted
//! verbatim under that header.

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

lazy_static! {
    // MM:Z:([ACGTUN][-+]([A-Za-z]+|[0-9]+)[.?]?(,[0-9]+)*;)*
    static ref MM_RE: Regex =
        Regex::new(r"((([ACGTUN])([-+])([A-Za-z]+|[0-9]+))[.?]?((,[0-9]+)*;)*)").unwrap();
}

/// Parse the MM/ML tags of `record` directly into `annot`. Positions are
/// stored in molecular (forward) orientation as 1bp annotations with a
/// single `Q` quality score per call.
///
/// Per-call filtering:
/// - calls with ML quality < `min_ml_score` are dropped.
/// - calls within `strip_at_ends` bp of either end of the forward
///   sequence are dropped (no-op if `strip_at_ends <= 0`).
///
/// Dispatch by MM mod_type (exact-match):
/// - `"a"` → [`M6A_TYPE`] (collapses `A+a`, `T-a`, `N+a`, etc.)
/// - `"m"` → [`CPG_TYPE`] (collapses `C+m`, `G+m`, etc.)
/// - any other → annotation type named verbatim after the MM group
///   header (e.g. `"C+h"` for 5hmC, `"C+ac"` for acetyl). [`write_mm_ml`]
///   emits these unchanged.
///
/// Idempotent: any pre-existing `m6a` / `cpg` / non-canonical basemod
/// types on `annot` are dropped first, so MM/ML always overwrites in-
/// memory state rather than appending. This guards against (a) repeated
/// calls on the same `annot`, and (b) malformed inputs that wrongly
/// emit `m6a` into the MA tag set in addition to MM/ML.
///
/// For canonical `m6a` / `cpg`, MM's per-group `(canonical_base,
/// strand_sign)` is reconstructed from the forward sequence at write
/// time, so input encodings like `N+a` are normalized to `A+a` / `T-a`
/// on output.
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
    annot
        .annotation_types
        .retain(|t| t.name != M6A_TYPE && t.name != CPG_TYPE && !is_valid_mm_header(&t.name));

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

    let mut m6a_calls: Vec<(u32, u8)> = Vec::new();
    let mut cpg_calls: Vec<(u32, u8)> = Vec::new();
    // Non-canonical groups: keyed by their MM header string (e.g. "C+h").
    let mut other_calls: BTreeMap<String, Vec<(u32, u8)>> = BTreeMap::new();
    let mut num_mods_seen = 0usize;

    for cap in MM_RE.captures_iter(mm_text) {
        let mod_base = cap.get(3).map(|m| m.as_str().as_bytes()[0]).unwrap();
        let strand_sign = cap
            .get(4)
            .map_or("+", |m| m.as_str())
            .chars()
            .next()
            .unwrap_or('+');
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

        // Filter + dispatch by mod_type. The canonical fiberseq codes are
        // *exactly* "a" (m6A) and "m" (5mC) — any longer alpha mod type
        // (`ac` acetyl, `pT` etc.) or ChEBI numeric ID is non-canonical
        // and stored under its verbatim MM header so it round-trips
        // exactly through write_mm_ml.
        let bucket: &mut Vec<(u32, u8)> = match mod_type_str.as_str() {
            "a" => &mut m6a_calls,
            "m" => &mut cpg_calls,
            _ => {
                let key = format!("{}{}{}", mod_base as char, strand_sign, mod_type_str);
                other_calls.entry(key).or_default()
            }
        };
        for (pos, &qual) in positions.iter().copied().zip(group_quals) {
            if qual < min_ml_score {
                continue;
            }
            let p = pos as usize;
            if strip > 0 && (p < strip || p >= upper) {
                continue;
            }
            bucket.push((pos, qual));
        }
    }

    if ml_tag.len() != num_mods_seen {
        log::warn!(
            "ML tag ({}) different number than MM tag ({}).",
            ml_tag.len(),
            num_mods_seen,
        );
    }

    add_calls(annot, M6A_TYPE, &mut m6a_calls);
    add_calls(annot, CPG_TYPE, &mut cpg_calls);
    for (name, mut calls) in other_calls {
        add_calls(annot, &name, &mut calls);
    }
}

fn add_calls(annot: &mut MolecularAnnotations, name: &str, calls: &mut [(u32, u8)]) {
    if calls.is_empty() {
        return;
    }
    calls.sort_by_key(|(p, _)| *p);
    let qspec = "Q".parse::<QualitySpec>().expect("Q parses");
    let t = annot.add_annotation_type(name, qspec);
    for &(pos, qual) in calls.iter() {
        t.add(pos, 1, Strand::Forward, vec![qual], None);
    }
}

/// Emit MM/ML tags from `annot` onto `record`. Existing `MM` and `ML`
/// aux are removed first.
///
/// Sources:
/// - [`M6A_TYPE`] / [`CPG_TYPE`] — canonical fiberseq types. MM groups
///   reconstructed from the forward sequence (`A`→`A+a`, `T`→`T-a`,
///   `C`→`C+m`, `G`→`G+m`). Calls landing on an unexpected base are
///   emitted under that base with `+` strand and a warning.
/// - Any other annotation type whose name is a valid MM group header
///   (e.g. `"C+h"`, `"N+a"`, `"T+12"`) is emitted verbatim under that
///   header. Skip counts are computed against the named base (`N`
///   counts every position).
/// - All other annotation types (`msp`, `nuc`, `fire`, etc.) are
///   skipped — they are not base modifications.
///
/// MM groups are emitted in deterministic (lexicographic header) order.
pub fn write_mm_ml(record: &mut bam::Record, annot: &MolecularAnnotations) {
    let forward_seq = if record.is_reverse() {
        revcomp(record.seq().as_bytes())
    } else {
        record.seq().as_bytes()
    };

    // header (e.g. "A+a", "C+h") → Vec<(pos, qual)>. BTreeMap gives
    // deterministic emission order.
    let mut groups: BTreeMap<String, Vec<(u32, u8)>> = BTreeMap::new();

    // Canonical types: reconstruct (base, strand) per call from forward seq.
    for (type_name, mod_type_char, canonical) in [
        (M6A_TYPE, 'a', &[(b'A', '+'), (b'T', '-')] as &[(u8, char)]),
        (CPG_TYPE, 'm', &[(b'C', '+'), (b'G', '+')] as &[(u8, char)]),
    ] {
        let Some(t) = annot.get_type(type_name) else {
            continue;
        };
        for a in &t.annotations {
            let seq_base = forward_seq.get(a.start as usize).copied().unwrap_or(b'N');
            let (group_base, group_strand) = canonical
                .iter()
                .find(|&&(b, _)| b == seq_base)
                .copied()
                .unwrap_or_else(|| {
                    log::warn!(
                        "{} call at pos {} on unexpected base {:?}; emitting as {}{}{}",
                        type_name,
                        a.start,
                        seq_base as char,
                        seq_base as char,
                        '+',
                        mod_type_char,
                    );
                    (seq_base, '+')
                });
            let header = format!("{}{}{}", group_base as char, group_strand, mod_type_char);
            let qual = a.qualities.first().copied().unwrap_or(0);
            groups.entry(header).or_default().push((a.start, qual));
        }
    }

    // Non-canonical types: name is the MM group header (e.g. "C+h").
    for t in annot.annotation_types.iter() {
        if t.name == M6A_TYPE || t.name == CPG_TYPE {
            continue;
        }
        if !is_valid_mm_header(&t.name) {
            continue;
        }
        for a in &t.annotations {
            let qual = a.qualities.first().copied().unwrap_or(0);
            groups
                .entry(t.name.clone())
                .or_default()
                .push((a.start, qual));
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

/// True if `name` looks like an MM group header: one canonical base
/// (`ACGTUN`), a strand sign, then a non-empty alphanumeric mod type.
fn is_valid_mm_header(name: &str) -> bool {
    let bytes = name.as_bytes();
    if bytes.len() < 3 {
        return false;
    }
    let base_ok = matches!(bytes[0], b'A' | b'C' | b'G' | b'T' | b'U' | b'N');
    let strand_ok = bytes[1] == b'+' || bytes[1] == b'-';
    let mod_ok = bytes[2..].iter().all(|b| b.is_ascii_alphanumeric());
    base_ok && strand_ok && mod_ok
}
