//! A group of annotations of the same type, plus the per-type tag-emission
//! helpers and their output structs.

use super::Qualities;
use crate::{Annotation, Encoding, MaEncoding, QualitySpec, Strand};

/// A group of annotations of the same type.
///
/// Annotation type identity within a [`MolecularAnnotations`](crate::MolecularAnnotations)
/// is keyed on `name` alone. Strand is a property of each [`Annotation`];
/// one type may contain annotations on different strands. `quality_spec` is
/// a per-type property but is not part of the identity — adding the same
/// `name` twice with conflicting `quality_spec` is an error.
#[derive(Debug, Clone, PartialEq)]
pub struct AnnotationType {
    /// Name of the annotation type (e.g., "msp", "nuc", "fire")
    pub name: String,
    /// Quality specification (number and scaling of quality values per annotation)
    pub quality_spec: QualitySpec,
    /// Individual annotations of this type
    pub annotations: Vec<Annotation>,
    /// Encoding format for this annotation type's on-disk representation
    pub encoding: Encoding,
}

impl AnnotationType {
    /// Create a new annotation type.
    pub fn new(name: impl Into<String>, quality_spec: QualitySpec) -> Self {
        Self {
            name: name.into(),
            quality_spec,
            annotations: Vec::new(),
            encoding: Encoding::default(),
        }
    }

    /// Set encoding on a fresh type. Panics if `self.annotations` is non-empty.
    pub fn set_encoding(&mut self, encoding: Encoding) -> &mut Self {
        assert!(
            self.annotations.is_empty(),
            "cannot change encoding on AnnotationType {:?}: already has {} annotations",
            self.name,
            self.annotations.len(),
        );
        self.encoding = encoding;
        self
    }

    /// True if this type's encoding is `MmMl`.
    pub fn is_mm_ml(&self) -> bool {
        matches!(self.encoding, Encoding::MmMl { .. })
    }

    /// Returns the on-disk section signature string for a given strand
    /// (e.g., `"msp+P"`, `"nuc-"`, `"fire.PQ"`). Used during serialization
    /// when annotations are grouped by `(name, strand)` into MA sections.
    pub fn type_signature(&self, strand: Strand) -> String {
        format!("{}{}{}", self.name, strand, self.quality_spec)
    }

    /// Add an annotation to this type.
    ///
    /// Accepts 0-based half-open [start, start+length) coordinates.
    /// `qualities` should have length equal to `quality_spec.num_qualities()`.
    pub fn add(
        &mut self,
        start: u32,
        length: u32,
        strand: Strand,
        qualities: Vec<u8>,
        name: Option<String>,
    ) -> &mut Self {
        self.annotations
            .push(Annotation::new(start, length, strand, qualities, name));
        self
    }

    /// Add an annotation with already-inline qualities and a shared name
    /// `Arc`. Intended for hot constructors that build many annotations
    /// sharing one name (e.g. the MM/ML parser): callers clone the `Arc`
    /// (a refcount bump) rather than allocating a fresh `String` per call.
    pub fn add_shared(
        &mut self,
        start: u32,
        length: u32,
        strand: Strand,
        qualities: Qualities,
        name: Option<std::sync::Arc<str>>,
    ) -> &mut Self {
        self.annotations.push(Annotation::with_shared(
            start, length, strand, qualities, name,
        ));
        self
    }
    /// Drop annotations for which `predicate` returns false.
    pub fn retain<F: FnMut(&Annotation) -> bool>(&mut self, predicate: F) {
        self.annotations.retain(predicate);
    }

    /// Per-type MA-tag emission.
    ///
    /// Returns `None` if `self.encoding` is not `Ma`. When emitting, all
    /// annotations of this type contribute, regardless of strand — strand
    /// is part of the section header (computed per-strand internally).
    ///
    /// Within the returned `MaParts`, sections are grouped by strand in
    /// `[Forward, Reverse, Unknown]` order. Within a section, annotation
    /// order follows insertion order in `self.annotations`.
    pub fn to_ma_parts(&self, layout: MaEncoding) -> Option<MaParts> {
        if !matches!(self.encoding, Encoding::Ma) {
            return None;
        }

        const STRANDS: [Strand; 3] = [Strand::Forward, Strand::Reverse, Strand::Unknown];

        let mut ma_section = String::new();
        let mut al_values = Vec::new();
        let mut aq_values = Vec::new();
        let mut an_values = Vec::new();

        for strand in STRANDS {
            let indices: Vec<usize> = self
                .annotations
                .iter()
                .enumerate()
                .filter(|(_, a)| a.strand == strand)
                .map(|(i, _)| i)
                .collect();
            if indices.is_empty() {
                continue;
            }

            let positions: Vec<String> = match layout {
                MaEncoding::Inline => indices
                    .iter()
                    .map(|&i| {
                        let a = &self.annotations[i];
                        format!("{}-{}", a.start + 1, a.length)
                    })
                    .collect(),
                MaEncoding::Separate => indices
                    .iter()
                    .map(|&i| (self.annotations[i].start + 1).to_string())
                    .collect(),
            };

            ma_section.push(';');
            ma_section.push_str(&self.type_signature(strand));
            ma_section.push(':');
            ma_section.push_str(&positions.join(","));

            if matches!(layout, MaEncoding::Separate) {
                al_values.extend(indices.iter().map(|&i| self.annotations[i].length));
            }
            if self.quality_spec.has_quality() {
                for &i in &indices {
                    aq_values.extend(self.annotations[i].qualities.iter().copied());
                }
            }
            for &i in &indices {
                an_values.push(
                    self.annotations[i]
                        .name
                        .as_deref()
                        .unwrap_or("")
                        .to_string(),
                );
            }
        }

        Some(MaParts {
            ma_section,
            al_values,
            aq_values,
            an_values,
        })
    }

    /// Per-type MM/ML emission.
    ///
    /// Returns `None` if `self.encoding` is not `MmMl`. Each MmMl-encoded
    /// annotation is expected to carry its skip-base as a one-char `name`
    /// (e.g. `Some("A")`), a single ML byte in `qualities`, and a 1bp
    /// length starting at the position in the forward sequence.
    ///
    /// Annotations are bucketed by `(skip_base, strand)`, sorted by position,
    /// and delta-encoded against `forward_seq`. ML bytes are emitted in
    /// the same order as the assembled groups.
    #[cfg(feature = "htslib")]
    pub fn to_mm_ml_parts(&self, forward_seq: &[u8]) -> Option<MmMlParts> {
        use crate::basemods::write::{delta_encode, skip_flag_str};
        use std::collections::BTreeMap;

        let skip_flag = match self.encoding {
            Encoding::MmMl { skip_flag } => skip_flag,
            Encoding::Ma => return None,
        };

        // Bucket calls by (skip_base, strand).
        // Use a tuple key with Strand replaced by an index (since Strand may not derive Ord).
        let strand_index = |s: Strand| -> u8 {
            match s {
                Strand::Forward => 0,
                Strand::Reverse => 1,
                Strand::Unknown => 2,
            }
        };
        type Bucket = (Strand, Vec<(u32, u8)>);
        let mut buckets: BTreeMap<(u8, u8), Bucket> = BTreeMap::new();
        for a in &self.annotations {
            let skip_base = a
                .name
                .as_deref()
                .and_then(|s| s.bytes().next())
                .map(|b| b.to_ascii_uppercase())
                .expect("MmMl-encoded annotation must carry skip_base in `name`");
            let qual = a.qualities.first().copied().unwrap_or(0);
            let entry = buckets
                .entry((skip_base, strand_index(a.strand)))
                .or_insert_with(|| (a.strand, Vec::new()));
            entry.1.push((a.start, qual));
        }

        let mut mm_groups = Vec::new();
        let mut ml_bytes_in_order = Vec::new();

        for ((skip_base, _), (strand, mut calls)) in buckets {
            calls.sort_by_key(|&(p, _)| p);
            let positions: Vec<u32> = calls.iter().map(|&(p, _)| p).collect();
            let deltas = delta_encode(&positions, forward_seq, skip_base);

            let strand_char = match strand {
                Strand::Forward => '+',
                Strand::Reverse => '-',
                Strand::Unknown => '+', // basemods should always be + or -; fall back.
            };
            let header = format!(
                "{}{}{}{}",
                skip_base as char,
                strand_char,
                self.name,
                skip_flag_str(skip_flag)
            );

            mm_groups.push(MmGroup {
                header,
                deltas,
                skip_base,
            });
            ml_bytes_in_order.extend(calls.into_iter().map(|(_, q)| q));
        }

        Some(MmMlParts {
            mm_groups,
            ml_bytes_in_order,
        })
    }
}

/// One AnnotationType's contribution to MA/AL/AQ/AN tag output.
#[derive(Debug, Clone, PartialEq)]
pub struct MaParts {
    /// MA tag fragment for this type, e.g. ";msp+P:101-50,201-60".
    pub ma_section: String,
    /// AL values (empty for Inline layout).
    pub al_values: Vec<u32>,
    /// AQ values (empty if `quality_spec` is None).
    pub aq_values: Vec<u8>,
    /// AN values (one per annotation; empty `String` if annotation has no name).
    pub an_values: Vec<String>,
}

/// One MM group's contribution to the assembled MM/ML tags.
#[derive(Debug, Clone, PartialEq)]
pub struct MmGroup {
    /// Group header, e.g. "A+a", "C+m.", "N+76792?".
    pub header: String,
    /// Delta-skip encoded position list.
    pub deltas: Vec<u32>,
    /// Skip-base of the group (useful for deterministic ordering across types).
    pub skip_base: u8,
}

/// One AnnotationType's contribution to MM/ML tag output.
#[derive(Debug, Clone, PartialEq)]
pub struct MmMlParts {
    pub mm_groups: Vec<MmGroup>,
    /// ML bytes for this type's calls, in the order matching `mm_groups`.
    pub ml_bytes_in_order: Vec<u8>,
}
