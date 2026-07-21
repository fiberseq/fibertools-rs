//! MM/ML base-modification tag I/O.
//!
//! This module is gated behind the `mmml` feature. The parse/emit algorithm
//! operates on plain byte slices (`parse::parse_mm_ml_into`,
//! `MolecularAnnotations::to_mm_ml`) so it is reachable without rust-htslib —
//! e.g. from the Python bindings, which supply the MM string, ML bytes, and
//! forward sequence via pysam.
//!
//! With the `htslib` feature also enabled, `MolecularAnnotations::from_record`
//! / `to_record` layer the rust_htslib `Record` bridge on top: they extract
//! those same byte slices off a record (revcomp-ing the sequence for reverse
//! reads) and delegate here.
//!
//! Parsing is not idempotent — `parse_mm_ml_into` appends into whatever types
//! already exist, so callers must run it on a container with no pre-existing
//! MM/ML types (`from_record` guarantees this by starting fresh; the Python
//! wrapper clears MM/ML types first).

pub(crate) mod parse;
pub(crate) mod write;
