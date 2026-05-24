//! MM/ML writer helpers. Used by `AnnotationType::to_mm_ml_parts`.

use crate::SkipFlag;

/// Render the skip-flag char (or empty if Implicit).
pub(crate) fn skip_flag_str(f: SkipFlag) -> &'static str {
    match f {
        SkipFlag::Implicit => "",
        SkipFlag::LowProbability => ".",
        SkipFlag::Unknown => "?",
    }
}

/// Compute delta-skip encoding for `positions` (already sorted ascending)
/// against `forward_seq`, counting bases of `skip_base` (or all bases when
/// `skip_base == b'N'`).
pub(crate) fn delta_encode(positions: &[u32], forward_seq: &[u8], skip_base: u8) -> Vec<u32> {
    let mut deltas = Vec::with_capacity(positions.len());
    let mut last_pos: usize = 0;
    for &p in positions {
        let p_usize = p as usize;
        let count = if last_pos < p_usize {
            forward_seq[last_pos..p_usize]
                .iter()
                .filter(|&&b| skip_base == b'N' || b.to_ascii_uppercase() == skip_base)
                .count()
        } else {
            0
        };
        deltas.push(count as u32);
        last_pos = p_usize + 1;
    }
    deltas
}
