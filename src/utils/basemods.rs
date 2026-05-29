//! Type-name constants and header helpers for fiberseq base modifications.
//!
//! Storage type names are spec MM mod codes (`"a"` for m6A, `"m"` for 5mC);
//! display labels in TSV output and CLI aliases use the fibertools-friendly
//! names (`"m6a"`, `"5mC"`).
//!
//! BAM tag I/O is owned by `utils::ma_io`; this module is purely a
//! constants + header-formatting helper module.

/// Annotation type name for m6A calls (spec MM mod code).
pub const M6A_TYPE: &str = "a";
/// Annotation type name for 5mC calls (spec MM mod code).
pub const CPG_TYPE: &str = "m";

/// True for annotation types whose on-disk source of truth is MM/ML.
/// Used to gate MA-tag emission (ma_io) and MM/ML encoding flips.
pub fn is_basemod_type(name: &str) -> bool {
    matches!(name, M6A_TYPE | CPG_TYPE)
}

/// Header for a canonical basemod call landing on `base` (read in
/// forward-strand orientation). Returns `&'static str` for the four
/// possible canonical headers; `None` if the (type, base) pair isn't a
/// canonical fiberseq combination. Used by programmatic producers
/// (`predict_m6a`, `ddda_to_m6a`) to label calls with the right MM
/// group header so the library's MM/ML serializer emits them under the
/// correct group.
pub fn canonical_header(type_name: &str, base: u8) -> Option<&'static str> {
    match (type_name, base) {
        (M6A_TYPE, b'A') => Some("A+a"),
        (M6A_TYPE, b'T') => Some("T-a"),
        (CPG_TYPE, b'C') => Some("C+m"),
        (CPG_TYPE, b'G') => Some("G+m"),
        _ => None,
    }
}
