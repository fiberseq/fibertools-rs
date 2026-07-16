//! MM/ML parser.
//!
//! `parse_mm_ml_into` operates purely on byte slices — the MM string, the ML
//! bytes, and the **forward** (original molecular orientation) sequence. It has
//! no rust-htslib dependency, so it is reachable under the `mmml` feature
//! alone. `MolecularAnnotations::from_record` (htslib) extracts those slices
//! off a record and delegates here.

use lazy_static::lazy_static;
use regex::Regex;
use std::sync::Arc;

lazy_static! {
    /// Matches one MM group: `[ACGTUN][-+]([a-z]+|[0-9]+)[.?]?(,[0-9]+)*;`.
    /// Capture groups:
    ///   1 = full group (incl. trailing ';')
    ///   2 = base+strand+modcodes (e.g. "A+a", "C+mh")
    ///   3 = base char
    ///   4 = strand char ('+' or '-')
    ///   5 = mod-code string (e.g. "a", "mh", "76792")
    ///   6 = optional skip-flag ('.' or '?'; may be empty)
    ///   7 = delta list (comma-separated digits, leading comma included)
    pub(crate) static ref MM_RE: Regex = Regex::new(
        r"((([ACGTUN])([-+])([A-Za-z]+|[0-9]+))([.?]?)((?:,[0-9]+)*);)"
    ).unwrap();
}

use crate::{AnnotationType, Encoding, MolecularAnnotations, QualitySpec, SkipFlag, Strand};

/// Parse MM/ML into `annot`, populating one `MmMl`-encoded `AnnotationType`
/// per unique mod code encountered.
///
/// * `mm_text` — the MM:Z tag value (delta-skip encoded, forward orientation).
/// * `ml_tag` — the ML:B,C bytes; may be empty (a warning is logged and missing
///   values are padded with 0).
/// * `forward_seq` — the read sequence in **original molecular orientation**
///   (callers holding a reverse-aligned record must revcomp before calling).
///   Case-insensitive; upper-cased internally.
///
/// No-op if `mm_text` is empty.
///
/// **Idempotency:** appends into existing types via `add_annotation_type`, so
/// this must run on a container with no pre-existing MM/ML types.
pub(crate) fn parse_mm_ml_into(
    annot: &mut MolecularAnnotations,
    mm_text: &str,
    ml_tag: &[u8],
    forward_seq: &[u8],
) {
    if mm_text.is_empty() {
        return;
    }

    let forward_bases: Vec<u8> = forward_seq.iter().map(|b| b.to_ascii_uppercase()).collect();

    let mut ml_cursor = 0usize;

    for cap in MM_RE.captures_iter(mm_text) {
        let skip_base = cap.get(3).unwrap().as_str().as_bytes()[0];
        let strand_char = cap.get(4).unwrap().as_str().as_bytes()[0];
        let mod_codes_str = cap.get(5).unwrap().as_str();
        let flag_str = cap.get(6).map(|m| m.as_str()).unwrap_or("");
        let delta_list = cap.get(7).map(|m| m.as_str()).unwrap_or("");

        let strand = match strand_char {
            b'+' => Strand::Forward,
            b'-' => Strand::Reverse,
            _ => continue,
        };
        let skip_flag = match flag_str {
            "." => SkipFlag::LowProbability,
            "?" => SkipFlag::Unknown,
            _ => SkipFlag::Implicit,
        };

        let mut deltas: Vec<i64> =
            Vec::with_capacity(delta_list.bytes().filter(|&b| b == b',').count());
        deltas.extend(
            delta_list
                .trim_start_matches(',')
                .split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse::<i64>().ok()),
        );

        if deltas.is_empty() {
            continue;
        }

        // Resolve positions in forward sequence by walking skip_base.
        let mut positions: Vec<u32> = vec![0; deltas.len()];
        let mut cur_idx = 0usize;
        let mut dist = 0i64;
        for (i, &b) in forward_bases.iter().enumerate() {
            if cur_idx >= deltas.len() {
                break;
            }
            let is_match = b == skip_base || skip_base == b'N';
            if is_match && dist == deltas[cur_idx] {
                positions[cur_idx] = i as u32;
                dist = 0;
                cur_idx += 1;
            } else if is_match {
                dist += 1;
            }
        }
        let resolved = cur_idx;
        if resolved != deltas.len() {
            log::warn!(
                "MM/ML parser: MM group ended {} short; truncating",
                deltas.len() - resolved
            );
        }

        // Determine mod-code list. If the mod-codes are alphabetic (each char
        // is one code), split per char. If numeric (ChEBI), treat as one code.
        let mod_codes: Vec<String> = if mod_codes_str.chars().all(|c| c.is_ascii_alphabetic()) {
            mod_codes_str.chars().map(|c| c.to_string()).collect()
        } else {
            vec![mod_codes_str.to_string()]
        };
        let codes_per_pos = mod_codes.len();

        // ML slice for this group: deltas.len() positions × codes_per_pos values,
        // interleaved per-position.
        let group_slot_count = deltas.len() * codes_per_pos;
        let group_end = ml_cursor + group_slot_count;
        let group_ml: Vec<u8> = if group_end > ml_tag.len() {
            let needed = group_end - ml_tag.len();
            let mut existing = ml_tag[ml_cursor..].to_vec();
            existing.extend(std::iter::repeat_n(0u8, needed));
            log::warn!(
                "ML tag too short; padding {} values with 0 after offset {}",
                needed,
                ml_cursor
            );
            existing
        } else {
            ml_tag[ml_cursor..group_end].to_vec()
        };
        ml_cursor = group_end;

        // Bucket ML per code: ML[i*codes + c] is for code c at position i.
        for (c_idx, code) in mod_codes.iter().enumerate() {
            // Basemod types are MmMl-encoded. If this code already exists with
            // a different skip flag, `add_annotation_type` keeps the first and
            // warns (the SAM spec lets the same mod code appear in multiple
            // groups).
            let t: &mut AnnotationType = annot.add_annotation_type(
                code,
                "Q".parse::<QualitySpec>().unwrap(),
                Encoding::MmMl { skip_flag },
            );

            let skip_base_name: Arc<str> = Arc::from((skip_base as char).to_string().as_str());
            t.annotations.reserve(resolved);
            for pos_idx in 0..resolved {
                let pos = positions[pos_idx];
                let ml_val = group_ml[pos_idx * codes_per_pos + c_idx];
                t.add_shared(
                    pos,
                    1,
                    strand,
                    smallvec::smallvec![ml_val],
                    Some(skip_base_name.clone()),
                );
            }
        }
    }

    if ml_cursor != ml_tag.len() {
        log::warn!(
            "ML tag has {} extra trailing bytes after parse",
            ml_tag.len().saturating_sub(ml_cursor)
        );
    }

    // Sort each MmMl type's annotations by start, for deterministic ordering.
    for t in annot.mm_ml_types_mut() {
        t.annotations.sort_by_key(|a| a.start);
    }
}
