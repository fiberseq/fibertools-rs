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

/// Canonical MM group header for an m6a call landing on `base` (read in
/// forward-strand orientation): `A+a` on A, `T-a` on T. `None` for any other
/// (type, base) pair.
///
/// m6a only — cpg (5mC) is never producer-synthesized; it travels the
/// read→write passthrough, which preserves the input's MM group verbatim, so
/// no canonical cpg header is ever needed here.
pub fn canonical_header(type_name: &str, base: u8) -> Option<&'static str> {
    match (type_name, base) {
        (M6A_TYPE, b'A') => Some("A+a"),
        (M6A_TYPE, b'T') => Some("T-a"),
        _ => None,
    }
}

/// Decomposed form of [`canonical_header`]: the `(skip_base, strand)` an m6a
/// call serializes under (`("A", Forward)` on A, `("T", Reverse)` on T).
///
/// Producers must store these on the annotation — `name` = skip-base,
/// `strand` = strand — NOT the full header string (e.g. `"T-a"`) in `name`
/// with a `Strand::Forward` placeholder. The MM/ML writer derives the group's
/// `+`/`-` sign from `strand`, not from `name`, so a `Strand::Forward`
/// placeholder silently turns a minus-strand call (`T-a`) into `T+a`. This
/// mirrors the parse path, which stores the skip-base in `name` and the strand
/// in `strand`.
///
/// m6a only — cpg is passthrough (see [`canonical_header`]).
pub fn canonical_basemod(
    type_name: &str,
    base: u8,
) -> Option<(&'static str, molecular_annotation::Strand)> {
    use molecular_annotation::Strand;
    match (type_name, base) {
        (M6A_TYPE, b'A') => Some(("A", Strand::Forward)),
        (M6A_TYPE, b'T') => Some(("T", Strand::Reverse)),
        _ => None,
    }
}
