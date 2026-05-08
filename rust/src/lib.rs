//! Molecular Annotation Library
//!
//! This library provides types and functions for working with molecular annotations
//! according to the MA/AL/AQ/AN tag specification for SAM/BAM files.
//!
//! # Coordinate Conventions
//!
//! This library uses two coordinate systems that you should understand:
//!
//! | Context | System | Example |
//! |---------|--------|---------|
//! | Internal/API | **0-based half-open** `[start, end)` | Position 100 with length 50 → `[100, 150)` |
//! | MA tag string | **1-based closed** `[start, end]` | Same annotation → `101-50` in tag |
//!
//! **Key points:**
//! - All Rust API methods (`add_annotations`, `get_coords`, etc.) use 0-based half-open
//! - The `from_tags`/`to_tags` methods automatically convert between conventions
//! - Liftover coordinates (`lift_to_reference`, `ref_coords`) are 0-based half-open
//!
//! ## Coordinate Orientations
//!
//! For reverse-aligned reads, there are two coordinate orientations:
//!
//! | Method | Orientation | Description |
//! |--------|-------------|-------------|
//! | `get_coords()` | BAM | Coordinates as stored in BAM (flipped for reverse reads) |
//! | `get_forward_coords()` | Molecular | Original read orientation (never flipped) |
//!
//! # Example
//! ```
//! use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec, Encoding};
//!
//! // Build annotations programmatically. Annotation type identity is keyed
//! // on `name` alone — strand is a per-annotation property.
//! let mut annotations = MolecularAnnotations::new(1000);
//!
//! // Add MSP annotations with phred quality scores
//! annotations
//!     .add_annotation_type("msp", "P".parse().unwrap())
//!     .add(100, 50, Strand::Forward, vec![40], None)
//!     .add(200, 60, Strand::Forward, vec![35], None);
//!
//! // Add nucleosome annotations without quality scores
//! annotations
//!     .add_annotation_type("nuc", QualitySpec::none())
//!     .add(150, 147, Strand::Forward, vec![], None)
//!     .add(350, 147, Strand::Forward, vec![], None);
//!
//! // Add FIRE annotation with linear quality
//! annotations
//!     .add_annotation_type("fire", "Q".parse().unwrap())
//!     .add(500, 75, Strand::Unknown, vec![200], Some("enhancer1".to_string()));
//!
//! // Add annotation type with multiple quality values per annotation (PQQP = 4 values each)
//! use molecular_annotation::QualityScaling;
//! let pqqp = QualitySpec::new(vec![
//!     QualityScaling::Phred, QualityScaling::Linear,
//!     QualityScaling::Linear, QualityScaling::Phred,
//! ]);
//! annotations
//!     .add_annotation_type("ctcf", pqqp)
//!     .add(600, 20, Strand::Forward, vec![40, 200, 180, 35], None);
//!
//! // Serialize with inline encoding (start-length pairs in MA string).
//! // Sections are emitted per (name, strand): types with annotations on
//! // multiple strands produce multiple sections sharing the same name.
//! annotations.set_encoding(Encoding::Inline);
//! let ma_inline = annotations.to_ma_string();
//! assert_eq!(ma_inline, "1000;msp+P:101-50,201-60;nuc+:151-147,351-147;fire.Q:501-75;ctcf+PQQP:601-20");
//!
//! // Serialize with separate encoding (starts in MA, lengths in AL array)
//! annotations.set_encoding(Encoding::Separate);
//! let ma_separate = annotations.to_ma_string();
//! let al_array = annotations.to_al_array();
//! assert_eq!(ma_separate, "1000;msp+P:101,201;nuc+:151,351;fire.Q:501;ctcf+PQQP:601");
//! assert_eq!(al_array, vec![50, 60, 147, 147, 75, 20]);
//!
//! // AQ tag: quality scores (only for types with quality spec)
//! // Values are grouped per-annotation: msp(40, 35), fire(200), ctcf(40, 200, 180, 35)
//! let aq_array = annotations.to_aq_array();
//! assert_eq!(aq_array, Some(vec![40, 35, 200, 40, 200, 180, 35]));
//!
//! // AN tag: names (empty for annotations without names)
//! let an_string = annotations.to_an_string();
//! assert_eq!(an_string, Some(",,,,enhancer1,".to_string()));
//! ```
//!
//! # Liftover Support
//!
//! The library supports lifting molecular coordinates to reference coordinates
//! using aligned blocks from BAM/CRAM records. All liftover coordinates use
//! **0-based half-open intervals** `[start, end)`.
//!
//! ```
//! use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec};
//! use molecular_annotation::liftover::AlignedBlocks;
//!
//! let mut annotations = MolecularAnnotations::new(1000);
//! annotations
//!     .add_annotation_type("msp", "P".parse().unwrap())
//!     .add(100, 50, Strand::Forward, vec![40], None);  // query [100, 150) in 0-based half-open
//!
//! // Set aligned blocks for liftover (query positions are forward-oriented)
//! // Second argument is is_reverse: false for forward-aligned reads
//! annotations.set_aligned_blocks(
//!     vec![([0, 500], [1000, 1500])],
//!     false,
//! );
//!
//! // Get reference coordinates on demand (0-based half-open)
//! if let Some(blocks) = annotations.aligned_blocks() {
//!     let annot = &annotations.annotation_types[0].annotations[0];
//!     let (ref_start, ref_end) = annot.ref_coords(blocks);
//!     // ref_start = Some(1100), ref_end = Some(1150)
//! }
//! ```

pub mod liftover;
mod types;

pub use liftover::{AlignedBlock, AlignedBlocks};
pub use types::{
    Annotation, AnnotationInfo, AnnotationType, Encoding, LiftedCoords, ParseError,
    QualityScaling, QualitySpec, Strand,
};

use std::str::FromStr;

/// Container for all molecular annotations on a read
#[derive(Debug, Clone, PartialEq)]
pub struct MolecularAnnotations {
    /// Length of the read at the time annotations were made
    pub read_length: u32,
    /// All annotation types on this read
    pub annotation_types: Vec<AnnotationType>,
    /// Encoding format for serialization
    encoding: Encoding,
    /// Optional aligned blocks for liftover calculations (private)
    aligned_blocks: Option<AlignedBlocks>,
    /// Whether the read is reverse-aligned. When true, MA coordinates
    /// (which are in original molecular orientation) need to be flipped
    /// before lifting to reference coordinates.
    is_reverse_aligned: bool,
}

impl MolecularAnnotations {
    /// Create a new empty molecular annotations container
    pub fn new(read_length: u32) -> Self {
        Self {
            read_length,
            annotation_types: Vec::new(),
            encoding: Encoding::default(),
            aligned_blocks: None,
            is_reverse_aligned: false,
        }
    }

    /// Get the current encoding format
    pub fn encoding(&self) -> Encoding {
        self.encoding
    }

    /// Set the encoding format for serialization (returns &mut Self for chaining)
    pub fn set_encoding(&mut self, encoding: Encoding) -> &mut Self {
        self.encoding = encoding;
        self
    }

    /// Get-or-create the annotation type with the given name.
    ///
    /// If a type with this `name` already exists with a matching
    /// `quality_spec`, returns the existing type. If it exists with a
    /// different `quality_spec`, **panics**. Programmatic callers that need
    /// fallible behavior should use [`try_add_annotation_type`](Self::try_add_annotation_type).
    ///
    /// Annotation type identity is keyed on `name` alone — strand is a
    /// per-annotation property, supplied via [`AnnotationType::add`].
    ///
    /// ```ignore
    /// annotations.add_annotation_type("msp", "P".parse().unwrap())
    ///     .add(100, 50, Strand::Forward, vec![40], None)
    ///     .add(200, 60, Strand::Reverse, vec![35], None);
    /// ```
    pub fn add_annotation_type(
        &mut self,
        name: &str,
        quality_spec: QualitySpec,
    ) -> &mut AnnotationType {
        self.try_add_annotation_type(name, quality_spec)
            .expect("conflicting quality_spec for existing annotation type")
    }

    /// Fallible variant of [`add_annotation_type`](Self::add_annotation_type).
    ///
    /// Returns the existing type if `name` is already present with a matching
    /// `quality_spec`. Returns [`ParseError::ConflictingAnnotationType`] if
    /// the existing type has a different `quality_spec`. Otherwise creates
    /// a new type and returns a mutable reference to it.
    pub fn try_add_annotation_type(
        &mut self,
        name: &str,
        quality_spec: QualitySpec,
    ) -> Result<&mut AnnotationType, ParseError> {
        if let Some(idx) = self.annotation_types.iter().position(|t| t.name == name) {
            if self.annotation_types[idx].quality_spec != quality_spec {
                return Err(ParseError::ConflictingAnnotationType {
                    name: name.to_string(),
                    reason: format!(
                        "existing quality_spec {} differs from requested {}",
                        self.annotation_types[idx].quality_spec, quality_spec
                    ),
                });
            }
            return Ok(&mut self.annotation_types[idx]);
        }
        self.annotation_types.push(AnnotationType::new(name, quality_spec));
        Ok(self.annotation_types.last_mut().unwrap())
    }

    /// Get the names of all annotation types.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{MolecularAnnotations, QualitySpec};
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations.add_annotation_type("msp", "P".parse().unwrap());
    /// annotations.add_annotation_type("nuc", QualitySpec::none());
    ///
    /// assert_eq!(annotations.annotation_type_names(), vec!["msp", "nuc"]);
    /// ```
    pub fn annotation_type_names(&self) -> Vec<&str> {
        self.annotation_types.iter().map(|t| t.name.as_str()).collect()
    }

    /// Get an annotation type by name
    pub fn get_type(&self, name: &str) -> Option<&AnnotationType> {
        self.annotation_types.iter().find(|t| t.name == name)
    }

    /// Get a mutable reference to an annotation type by name
    pub fn get_type_mut(&mut self, name: &str) -> Option<&mut AnnotationType> {
        self.annotation_types.iter_mut().find(|t| t.name == name)
    }

    /// Add multiple annotations of the same type at once.
    ///
    /// Accepts 0-based half-open [start, start+length) coordinates. The
    /// `strand` parameter is applied to every annotation in this batch.
    /// Callers needing mixed strands per type call `add_annotations` once
    /// per strand — both calls land in the same in-memory `AnnotationType`.
    ///
    /// # Arguments
    /// * `name` - Name of the annotation type (e.g., "msp", "nuc", "fire")
    /// * `quality_spec` - Quality specification
    /// * `starts` - 0-based start positions
    /// * `lengths` - Lengths of the annotations in base pairs
    /// * `strand` - Strand applied to every annotation in this batch
    /// * `qualities` - Optional flat quality array. Length must be
    ///   `starts.len() * quality_spec.num_qualities()`.
    /// * `names` - Optional names/labels, must match starts length if provided
    ///
    /// # Errors
    /// Returns an error if array lengths don't match, or if `name` already
    /// exists with a different `quality_spec`.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec};
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations.add_annotations(
    ///     "msp",
    ///     "P".parse().unwrap(),
    ///     &[100, 200, 300],
    ///     &[50, 60, 70],
    ///     Strand::Forward,
    ///     Some(&[40, 35, 30]),
    ///     None,
    /// ).unwrap();
    ///
    /// assert_eq!(annotations.total_annotation_count(), 3);
    /// ```
    #[allow(clippy::too_many_arguments)]
    pub fn add_annotations(
        &mut self,
        name: &str,
        quality_spec: QualitySpec,
        starts: &[u32],
        lengths: &[u32],
        strand: Strand,
        qualities: Option<&[u8]>,
        names: Option<&[String]>,
    ) -> Result<&mut Self, ParseError> {
        if starts.len() != lengths.len() {
            return Err(ParseError::MismatchedArrayLengths {
                expected: starts.len(),
                got: lengths.len(),
            });
        }
        let num_q = quality_spec.num_qualities();
        if let Some(q) = qualities {
            let expected_len = starts.len() * num_q;
            if q.len() != expected_len {
                return Err(ParseError::MismatchedArrayLengths {
                    expected: expected_len,
                    got: q.len(),
                });
            }
        }
        if let Some(n) = names {
            if n.len() != starts.len() {
                return Err(ParseError::MismatchedArrayLengths {
                    expected: starts.len(),
                    got: n.len(),
                });
            }
        }

        let at = self.try_add_annotation_type(name, quality_spec)?;

        for i in 0..starts.len() {
            let q = if num_q > 0 {
                qualities
                    .map(|qs| qs[i * num_q..(i + 1) * num_q].to_vec())
                    .unwrap_or_default()
            } else {
                Vec::new()
            };
            let n = names.map(|n| n[i].clone());
            at.add(starts[i], lengths[i], strand, q, n);
        }

        Ok(self)
    }

    /// Returns the total number of annotations across all types
    pub fn total_annotation_count(&self) -> usize {
        self.annotation_types.iter().map(|t| t.annotations.len()).sum()
    }

    /// Iterate over all annotations across all types
    pub fn iter_all_annotations(&self) -> impl Iterator<Item = (&AnnotationType, &Annotation)> {
        self.annotation_types
            .iter()
            .flat_map(|t| t.annotations.iter().map(move |a| (t, a)))
    }

    /// Iterate over annotations for a specific type with full coordinate information.
    ///
    /// This method provides comprehensive information for each annotation including:
    /// - Type metadata (name, strand, quality spec)
    /// - Coordinates in both BAM and molecular orientations
    /// - Reference coordinates (if aligned blocks are set)
    /// - Quality values and name fields
    ///
    /// Returns `None` if the type doesn't exist.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec};
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations
    ///     .add_annotation_type("msp", "P".parse().unwrap())
    ///     .add(100, 50, Strand::Forward, vec![40], None);
    ///
    /// // Collect into a Vec to avoid lifetime issues
    /// let msp_annotations: Vec<_> = annotations.iter_type("msp")
    ///     .map(|iter| iter.collect::<Vec<_>>())
    ///     .unwrap_or_default();
    ///
    /// for info in &msp_annotations {
    ///     println!("[{}, {}) qualities={:?}",
    ///         info.query_start,
    ///         info.query_end,
    ///         info.qualities
    ///     );
    /// }
    /// ```
    pub fn iter_type(&self, type_name: &str) -> Option<impl Iterator<Item = AnnotationInfo<'_>>> {
        let annot_type = self.get_type(type_name)?;
        Some(self.iter_annotation_type(annot_type))
    }

    /// Internal helper to iterate over an annotation type.
    fn iter_annotation_type<'a>(
        &'a self,
        annot_type: &'a AnnotationType,
    ) -> impl Iterator<Item = AnnotationInfo<'a>> {
        annot_type.annotations.iter().map(move |a| {
            let (query_start, query_end): (u32, u32) = if self.is_reverse_aligned {
                self.flip_range(a.start, a.end())
            } else {
                (a.start, a.end())
            };

            let (ref_start, ref_end) = self
                .aligned_blocks
                .as_ref()
                .map(|b| b.lift_to_reference(query_start, query_end))
                .unwrap_or((None, None));

            AnnotationInfo {
                type_name: &annot_type.name,
                strand: a.strand,
                quality_spec: &annot_type.quality_spec,
                query_start,
                query_end,
                forward_start: a.start,
                forward_end: a.end(),
                ref_start,
                ref_end,
                qualities: &a.qualities,
                name: a.name.as_deref(),
            }
        })
    }

    /// Iterate over all annotations across all types with full coordinate information.
    ///
    /// This method provides comprehensive information for each annotation including:
    /// - Type metadata (name, strand, quality spec)
    /// - Coordinates in both BAM and molecular orientations
    /// - Reference coordinates (if aligned blocks are set)
    /// - Quality values and name fields
    ///
    /// This is the recommended method when you need all annotation details at once,
    /// such as for export or analysis. For iterating over a single type, use
    /// [`iter_type`](Self::iter_type).
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec};
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations
    ///     .add_annotation_type("msp", "P".parse().unwrap())
    ///     .add(100, 50, Strand::Forward, vec![40], None);
    ///
    /// for info in annotations.iter_full() {
    ///     println!("{}: [{}, {}) qualities={:?}",
    ///         info.type_name,
    ///         info.query_start,
    ///         info.query_end,
    ///         info.qualities
    ///     );
    /// }
    /// ```
    pub fn iter_full(&self) -> impl Iterator<Item = AnnotationInfo<'_>> {
        self.annotation_types
            .iter()
            .flat_map(move |annot_type| self.iter_annotation_type(annot_type))
    }

    /// Section emission order for MA / AL / AQ / AN serialization.
    ///
    /// For each annotation type (in insertion order), groups its annotations
    /// by [`Strand`] in enum order (`Forward`, `Reverse`, `Unknown`). Each
    /// non-empty group becomes one section. Within a group, annotation order
    /// follows insertion order.
    ///
    /// Yields `(annotation_type, strand, indices)` where `indices` are
    /// positions into `annotation_type.annotations`.
    fn emission_sections(&self) -> impl Iterator<Item = (&AnnotationType, Strand, Vec<usize>)> + '_ {
        const STRANDS: [Strand; 3] = [Strand::Forward, Strand::Reverse, Strand::Unknown];
        self.annotation_types.iter().flat_map(|t| {
            STRANDS.into_iter().filter_map(move |s| {
                let indices: Vec<usize> = t
                    .annotations
                    .iter()
                    .enumerate()
                    .filter(|(_, a)| a.strand == s)
                    .map(|(i, _)| i)
                    .collect();
                if indices.is_empty() {
                    None
                } else {
                    Some((t, s, indices))
                }
            })
        })
    }

    /// Generate the MA:Z tag string.
    ///
    /// Converts internal 0-based half-open to 1-based closed per MA tag spec
    /// (internal `start=99` → MA tag `100`). Annotations within each type
    /// are grouped by strand into sections; types with annotations on
    /// multiple strands emit multiple sections sharing the same `name`.
    /// Empty types emit no section.
    pub fn to_ma_string(&self) -> String {
        let mut parts = vec![self.read_length.to_string()];

        for (annot_type, strand, indices) in self.emission_sections() {
            let positions: Vec<String> = match self.encoding {
                Encoding::Inline => indices
                    .iter()
                    .map(|&i| {
                        let a = &annot_type.annotations[i];
                        format!("{}-{}", a.start + 1, a.length)
                    })
                    .collect(),
                Encoding::Separate => indices
                    .iter()
                    .map(|&i| (annot_type.annotations[i].start + 1).to_string())
                    .collect(),
            };

            parts.push(format!("{}:{}", annot_type.type_signature(strand), positions.join(",")));
        }

        parts.join(";")
    }

    /// Generate the AL:B:I array (lengths). Order matches MA section order.
    pub fn to_al_array(&self) -> Vec<u32> {
        self.emission_sections()
            .flat_map(|(t, _strand, indices)| {
                indices.into_iter().map(move |i| t.annotations[i].length)
            })
            .collect()
    }

    /// Generate the AQ:B:C array (None if no annotations have quality).
    /// Order matches MA section order.
    pub fn to_aq_array(&self) -> Option<Vec<u8>> {
        let qualities: Vec<u8> = self
            .emission_sections()
            .filter(|(t, _, _)| t.quality_spec.has_quality())
            .flat_map(|(t, _strand, indices)| {
                indices
                    .into_iter()
                    .flat_map(move |i| t.annotations[i].qualities.iter().copied())
            })
            .collect();

        if qualities.is_empty() {
            None
        } else {
            Some(qualities)
        }
    }

    /// Generate the AN:Z tag string (None if no annotations have names).
    /// Order matches MA section order.
    pub fn to_an_string(&self) -> Option<String> {
        let has_any_names = self
            .annotation_types
            .iter()
            .any(|t| t.annotations.iter().any(|a| a.name.is_some()));

        if !has_any_names {
            return None;
        }

        let names: Vec<String> = self
            .emission_sections()
            .flat_map(|(t, _strand, indices)| {
                indices
                    .into_iter()
                    .map(move |i| t.annotations[i].name.as_deref().unwrap_or("").to_string())
            })
            .collect();

        Some(names.join(","))
    }

    /// Generate all BAM tag values at once.
    ///
    /// This is the recommended method for serializing annotations to BAM tags.
    /// Returns a tuple of (MA, AL, AQ, AN) where:
    /// - MA: The MA:Z tag string (always present)
    /// - AL: The AL:B:I array (lengths, empty if using inline encoding)
    /// - AQ: The AQ:B:C array (qualities, None if no annotations have quality)
    /// - AN: The AN:Z tag string (names, None if no annotations have names)
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec, Encoding};
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations
    ///     .add_annotation_type("msp", "P".parse().unwrap())
    ///     .add(99, 50, Strand::Forward, vec![40], None);  // 0-based
    ///
    /// let (ma, al, aq, an) = annotations.to_tags();
    /// assert_eq!(ma, "1000;msp+P:100-50");  // 1-based in tag
    /// assert!(al.is_empty());  // Inline encoding doesn't need AL
    /// assert_eq!(aq, Some(vec![40]));
    /// assert_eq!(an, None);
    /// ```
    pub fn to_tags(&self) -> (String, Vec<u32>, Option<Vec<u8>>, Option<String>) {
        let ma = self.to_ma_string();
        let al = if matches!(self.encoding, Encoding::Separate) {
            self.to_al_array()
        } else {
            Vec::new()
        };
        let aq = self.to_aq_array();
        let an = self.to_an_string();
        (ma, al, aq, an)
    }

    /// Write annotations to a BAM record.
    ///
    /// This sets the MA:Z tag, and optionally AL:B:I, AQ:B:C, and AN:Z tags
    /// depending on the encoding format and whether quality/names are present.
    ///
    /// # Example
    /// ```ignore
    /// use rust_htslib::bam::{self, Read};
    /// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec};
    ///
    /// let mut record = /* from BAM file */;
    /// let mut annotations = MolecularAnnotations::from_record(&record);
    /// annotations
    ///     .add_annotation_type("msp", Strand::Forward, "P".parse().unwrap())
    ///     .add(100, 50, vec![40], None);
    ///
    /// annotations.to_record(&mut record);
    /// ```
    #[cfg(feature = "htslib")]
    pub fn to_record(&self, record: &mut rust_htslib::bam::Record) {
        use rust_htslib::bam::record::Aux;

        let (ma, al, aq, an) = self.to_tags();

        // Set MA tag (always)
        record.push_aux(b"MA", Aux::String(&ma)).ok();

        // Set AL tag (only for separate encoding)
        if !al.is_empty() {
            let al_u32: Vec<u32> = al.iter().map(|&v| v as u32).collect();
            record.push_aux(b"AL", Aux::ArrayU32((&al_u32).into())).ok();
        }

        // Set AQ tag (if present)
        if let Some(ref aq_arr) = aq {
            record.push_aux(b"AQ", Aux::ArrayU8(aq_arr.into())).ok();
        }

        // Set AN tag (if present)
        if let Some(ref an_str) = an {
            record.push_aux(b"AN", Aux::String(an_str)).ok();
        }
    }

    // --- Aligned Blocks / Liftover Methods ---

    /// Set aligned blocks for liftover calculations (from AlignedBlocks object).
    ///
    /// These blocks enable reference coordinate computation via getter methods
    /// on individual annotations. This method accepts a pre-constructed `AlignedBlocks`
    /// object. For convenience, use [`set_aligned_blocks`](Self::set_aligned_blocks)
    /// which accepts raw block pairs.
    ///
    /// # Arguments
    /// * `blocks` - The aligned blocks for coordinate transformation.
    /// * `is_reverse` - Whether the read is reverse-aligned. When true, MA coordinates
    ///   (in original molecular orientation) will be flipped before lifting to reference.
    pub fn set_aligned_blocks_raw(&mut self, blocks: AlignedBlocks, is_reverse: bool) {
        self.aligned_blocks = Some(blocks);
        self.is_reverse_aligned = is_reverse;
    }

    /// Set aligned blocks for liftover calculations.
    ///
    /// Accepts 0-based half-open [start, end) intervals.
    ///
    /// # Arguments
    /// * `blocks` - Vector of `([query_start, query_end], [ref_start, ref_end])` pairs.
    /// * `is_reverse` - Whether the read is reverse-aligned.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::MolecularAnnotations;
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations.set_aligned_blocks(
    ///     vec![([0, 500], [1000, 1500])],
    ///     false,  // not reverse-aligned
    /// );
    /// ```
    pub fn set_aligned_blocks(&mut self, blocks: Vec<([u32; 2], [u32; 2])>, is_reverse: bool) {
        self.aligned_blocks = Some(AlignedBlocks::new(blocks, self.read_length));
        self.is_reverse_aligned = is_reverse;
    }

    /// Clear aligned blocks.
    pub fn clear_aligned_blocks(&mut self) {
        self.aligned_blocks = None;
    }

    /// Check if aligned blocks are set.
    pub fn has_aligned_blocks(&self) -> bool {
        self.aligned_blocks.is_some()
    }

    /// Get reference to aligned blocks.
    pub fn aligned_blocks(&self) -> Option<&AlignedBlocks> {
        self.aligned_blocks.as_ref()
    }

    /// Check if the read is reverse-aligned.
    pub fn is_reverse_aligned(&self) -> bool {
        self.is_reverse_aligned
    }

    /// Set whether the read is reverse-aligned.
    ///
    /// When true, annotation coordinates (which are in original molecular orientation)
    /// will be flipped before lifting to reference coordinates.
    pub fn set_reverse_aligned(&mut self, is_reverse: bool) {
        self.is_reverse_aligned = is_reverse;
    }

    /// Flip a coordinate range for reverse-aligned reads.
    ///
    /// For a read of length L, flips `[start, end)` to `[L - end, L - start)`.
    /// This converts coordinates from original molecular orientation to
    /// forward-oriented alignment coordinates.
    #[inline]
    fn flip_range(&self, start: u32, end: u32) -> (u32, u32) {
        (self.read_length - end, self.read_length - start)
    }

    /// Get coordinates for a specific annotation type in BAM orientation.
    ///
    /// Query coordinates are returned in **BAM orientation** (forward-oriented, matching
    /// the sequence as stored in the BAM file). For reverse-aligned reads, this means
    /// the coordinates are flipped from the original molecular orientation.
    ///
    /// This is analogous to pysam's `modified_bases` which returns positions relative
    /// to the BAM sequence. For original molecular orientation, use
    /// [`get_forward_coords`](Self::get_forward_coords) (alias: `get_molecular_coords`).
    ///
    /// # Returns
    /// Vector of `(query_start, query_end)` tuples as 0-based half-open intervals.
    /// Returns `None` if the type doesn't exist.
    ///
    /// # Aliases
    /// - [`get_bam_coords`](Self::get_bam_coords) - Same method, clearer name
    pub fn get_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        let annot_type = self.get_type(type_name)?;

        Some(
            annot_type
                .annotations
                .iter()
                .map(|a| {
                    if self.is_reverse_aligned {
                        self.flip_range(a.start, a.end())
                    } else {
                        (a.start, a.end())
                    }
                })
                .collect(),
        )
    }

    /// Alias for [`get_coords`](Self::get_coords) with a clearer name.
    ///
    /// Returns coordinates in BAM orientation (flipped for reverse-aligned reads).
    #[inline]
    pub fn get_bam_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        self.get_coords(type_name)
    }

    /// Get coordinates for a specific annotation type in molecular (forward) orientation.
    ///
    /// Returns coordinates relative to the original read sequence, before any reverse
    /// complementation. These are the raw stored coordinates.
    ///
    /// This is analogous to pysam's `modified_bases_forward` which returns positions
    /// relative to the original sequence.
    ///
    /// # Returns
    /// Vector of `(query_start, query_end)` tuples as 0-based half-open intervals.
    /// Returns `None` if the type doesn't exist.
    ///
    /// # Aliases
    /// - [`get_molecular_coords`](Self::get_molecular_coords) - Same method, clearer name
    pub fn get_forward_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        let annot_type = self.get_type(type_name)?;
        Some(
            annot_type
                .annotations
                .iter()
                .map(|a| (a.start, a.end()))
                .collect(),
        )
    }

    /// Alias for [`get_forward_coords`](Self::get_forward_coords) with a clearer name.
    ///
    /// Returns coordinates in molecular orientation (original read sequence, never flipped).
    #[inline]
    pub fn get_molecular_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        self.get_forward_coords(type_name)
    }

    /// Get lifted coordinates for a specific annotation type using internal aligned blocks.
    ///
    /// Query coordinates are returned in **BAM orientation** (forward-oriented, matching
    /// the sequence as stored in the BAM file). For reverse-aligned reads, this means
    /// the coordinates are flipped from the original molecular orientation.
    ///
    /// This is analogous to pysam's `modified_bases` which returns positions relative
    /// to the BAM sequence. For original molecular orientation, access the annotations
    /// directly via `annotation_types`.
    ///
    /// # Returns
    /// Vector of [`LiftedCoords`] tuples `(query_start, query_end, ref_start, ref_end)`,
    /// all as 0-based half-open intervals.
    /// Returns `None` if the type doesn't exist or aligned blocks are not set.
    pub fn get_ref_coords(
        &self,
        type_name: &str,
    ) -> Option<Vec<LiftedCoords>> {
        let blocks = self.aligned_blocks.as_ref()?;
        let coords = self.get_coords(type_name)?;

        Some(
            coords
                .into_iter()
                .map(|(query_start, query_end)| {
                    let (ref_start, ref_end) = blocks.lift_to_reference(query_start, query_end);
                    (query_start, query_end, ref_start, ref_end)
                })
                .collect(),
        )
    }

    /// Lift a range from query to reference using internal aligned blocks.
    ///
    /// **Note:** This method does NOT automatically flip coordinates for reverse-aligned
    /// reads. The caller must provide forward-oriented query coordinates.
    /// For lifting annotation coordinates that automatically handle reverse strand,
    /// use [`get_ref_coords`](Self::get_ref_coords).
    ///
    /// # Arguments
    /// * `start` - 0-based query start position (inclusive), forward-oriented
    /// * `end` - 0-based query end position (exclusive), forward-oriented
    ///
    /// # Returns
    /// Tuple of `(ref_start, ref_end)` as 0-based half-open interval,
    /// or `None` if aligned blocks are not set.
    pub fn lift_to_reference(&self, start: u32, end: u32) -> Option<(Option<u32>, Option<u32>)> {
        let blocks = self.aligned_blocks.as_ref()?;
        Some(blocks.lift_to_reference(start, end))
    }

    /// Lift a range from reference to query using internal aligned blocks.
    ///
    /// # Arguments
    /// * `start` - 0-based reference start position (inclusive)
    /// * `end` - 0-based reference end position (exclusive)
    ///
    /// # Returns
    /// Tuple of `(query_start, query_end)` as 0-based half-open interval,
    /// or `None` if aligned blocks are not set.
    pub fn lift_to_query(&self, start: u32, end: u32) -> Option<(Option<u32>, Option<u32>)> {
        let blocks = self.aligned_blocks.as_ref()?;
        Some(blocks.lift_to_query(start, end))
    }

    /// Create a new MolecularAnnotations from a BAM record.
    ///
    /// This extracts the read length, aligned blocks, and reverse-alignment status
    /// from the record. Annotation types can then be added manually or parsed
    /// from tags using methods like `add_annotation_type`.
    ///
    /// # Example
    /// ```ignore
    /// use rust_htslib::bam::{self, Read};
    /// use molecular_annotation::MolecularAnnotations;
    ///
    /// let mut bam = bam::Reader::from_path("test.bam").unwrap();
    /// for record in bam.records() {
    ///     let record = record.unwrap();
    ///     let annotations = MolecularAnnotations::from_record(&record);
    ///     // annotations now has aligned blocks set and is_reverse_aligned configured
    /// }
    /// ```
    #[cfg(feature = "htslib")]
    pub fn from_record(record: &rust_htslib::bam::Record) -> Self {
        Self {
            read_length: record.seq_len() as u32,
            annotation_types: Vec::new(),
            encoding: Encoding::default(),
            aligned_blocks: Some(AlignedBlocks::from_record(record)),
            is_reverse_aligned: record.is_reverse(),
        }
    }

    /// Parse from tag values.
    ///
    /// Converts 1-based closed MA tag positions to 0-based half-open internally
    /// (MA tag `100` → internal `start=99`).
    ///
    /// # Arguments
    /// * `ma` - The MA:Z tag value (e.g., "1000;msp+P:100,200"). Positions are 1-based.
    /// * `al` - The AL:B:I array values (lengths). Empty if using inline format.
    /// * `aq` - Optional AQ:B:C array values (qualities).
    /// * `an` - Optional AN:Z tag value (names).
    pub fn from_tags(
        ma: &str,
        al: &[u32],
        aq: Option<&[u8]>,
        an: Option<&str>,
    ) -> Result<Self, ParseError> {
        // Parse MA tag
        let parts: Vec<&str> = ma.split(';').collect();
        if parts.is_empty() {
            return Err(ParseError::InvalidFormat("Empty MA tag".to_string()));
        }

        // First part is read length
        let read_length: u32 = parts[0]
            .parse()
            .map_err(|_| ParseError::InvalidFormat(format!("Invalid read length: {}", parts[0])))?;

        // Parse annotation type names (optional, from AN tag)
        let names: Vec<&str> = an.map(|s| s.split(',').collect()).unwrap_or_default();

        // Track position in AL and AQ arrays
        let mut al_idx = 0usize;
        let mut aq_idx = 0usize;
        let mut name_idx = 0usize;

        // Build incrementally so try_add_annotation_type is the dedup
        // choke point: same `name` across multiple sections accumulates
        // into one in-memory AnnotationType, with per-annotation strand
        // carried over from each section header.
        let mut out = MolecularAnnotations::new(read_length);

        // Parse each annotation type section
        for part in &parts[1..] {
            if part.is_empty() {
                continue;
            }

            // Split on ':' to get type info and positions
            let type_and_positions: Vec<&str> = part.splitn(2, ':').collect();
            if type_and_positions.len() != 2 {
                return Err(ParseError::InvalidFormat(format!("Invalid annotation type format: {}", part)));
            }

            let type_info = type_and_positions[0];
            let positions_str = type_and_positions[1];

            // Parse type info: name + strand + quality spec
            // Format: name[+-.]([PQ]*)
            let (name, strand, quality_spec) = parse_type_info(type_info)?;
            let num_q = quality_spec.num_qualities();

            // Parse positions (and optionally lengths if inline format)
            // Detect format: "100-50,200-60" (inline) vs "100,200" (separate)
            // MA tag uses 1-based positions, convert to 0-based internally
            let position_length_pairs: Vec<(u32, Option<u32>)> = positions_str
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    if let Some((start_str, len_str)) = s.split_once('-') {
                        let start: u32 = start_str.parse().map_err(|_| ParseError::InvalidCoordinate(format!("Invalid start: {}", start_str)))?;
                        let len: u32 = len_str.parse().map_err(|_| ParseError::InvalidCoordinate(format!("Invalid length: {}", len_str)))?;
                        Ok((start - 1, Some(len)))  // Convert 1-based to 0-based
                    } else {
                        let start: u32 = s.parse().map_err(|_| ParseError::InvalidCoordinate(format!("Invalid position: {}", s)))?;
                        Ok((start - 1, None))  // Convert 1-based to 0-based
                    }
                })
                .collect::<Result<Vec<_>, ParseError>>()?;

            // Get-or-merge into the existing AnnotationType for `name`.
            let at = out.try_add_annotation_type(&name, quality_spec.clone())?;

            for (pos, inline_length) in position_length_pairs {
                // Get length: from inline format or from AL array
                let length = if let Some(len) = inline_length {
                    len
                } else {
                    if al_idx >= al.len() {
                        return Err(ParseError::MismatchedArrayLengths {
                            expected: al_idx + 1,
                            got: al.len(),
                        });
                    }
                    let len = al[al_idx];
                    al_idx += 1;
                    len
                };

                // Get quality values if this type has quality scores
                let qualities = if num_q > 0 {
                    if let Some(aq_arr) = aq {
                        if aq_idx + num_q > aq_arr.len() {
                            return Err(ParseError::MismatchedArrayLengths {
                                expected: aq_idx + num_q,
                                got: aq_arr.len(),
                            });
                        }
                        let q = aq_arr[aq_idx..aq_idx + num_q].to_vec();
                        aq_idx += num_q;
                        q
                    } else {
                        return Err(ParseError::InvalidFormat(
                            "Quality type specified but no AQ array provided".to_string(),
                        ));
                    }
                } else {
                    Vec::new()
                };

                // Get name if available
                let annot_name = if name_idx < names.len() {
                    let n = names[name_idx];
                    name_idx += 1;
                    if n.is_empty() { None } else { Some(n.to_string()) }
                } else {
                    name_idx += 1;
                    None
                };

                at.add(pos, length, strand, qualities, annot_name);
            }
        }

        // Validate that we consumed all AL values (only if using separate format)
        if al_idx > 0 && al_idx != al.len() {
            return Err(ParseError::MismatchedArrayLengths {
                expected: al_idx,
                got: al.len(),
            });
        }

        // Detect encoding from input format
        out.encoding = if al_idx == 0 && !out.annotation_types.is_empty() {
            Encoding::Inline
        } else {
            Encoding::Separate
        };

        Ok(out)
    }
}

/// Parse the type info string (e.g., "msp+P", "nuc-", "fire.PQ")
pub(crate) fn parse_type_info(s: &str) -> Result<(String, Strand, QualitySpec), ParseError> {
    // Find the strand character (+, -, .)
    let strand_pos = s
        .char_indices()
        .find(|(_, c)| *c == '+' || *c == '-' || *c == '.')
        .map(|(i, _)| i);

    let strand_pos = strand_pos.ok_or_else(|| {
        ParseError::InvalidFormat(format!("No strand indicator found in: {}", s))
    })?;

    let name = &s[..strand_pos];
    if name.is_empty() {
        return Err(ParseError::InvalidFormat("Empty annotation type name".to_string()));
    }

    let strand_char = &s[strand_pos..strand_pos + 1];
    let strand = Strand::from_str(strand_char)?;

    let quality_str = &s[strand_pos + 1..];
    let quality_spec = QualitySpec::from_str(quality_str)?;

    Ok((name.to_string(), strand, quality_spec))
}

#[cfg(test)]
mod tests;
