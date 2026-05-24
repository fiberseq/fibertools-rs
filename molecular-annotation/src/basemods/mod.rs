//! MM/ML base-modification tag I/O.
//!
//! This module is gated behind the `htslib` feature. Public entry points
//! are exposed only through `MolecularAnnotations::from_record` (parse)
//! and `MolecularAnnotations::to_record` (write). There is intentionally
//! no public `parse_mm_ml` or `write_mm_ml` free function — parsing only
//! happens on construction, eliminating idempotency questions.

pub(crate) mod parse;
pub(crate) mod write;
