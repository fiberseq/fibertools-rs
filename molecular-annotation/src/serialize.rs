//! MA/AQ/AN and MM/ML tag emission for [`MolecularAnnotations`].

use crate::{Encoding, MaParts, MolecularAnnotations};

impl MolecularAnnotations {
    /// Iterator over per-type [`MaParts`] for every MA-encoded annotation type.
    ///
    /// Types whose encoding is not `Encoding::Ma` are skipped (they don't
    /// participate in MA-tag emission).
    fn ma_parts_iter(&self) -> impl Iterator<Item = MaParts> + '_ {
        self.annotation_types
            .iter()
            .filter_map(|t| t.to_ma_parts())
    }

    /// Generate the MA:Z tag string.
    ///
    /// Converts internal 0-based half-open to 1-based closed per MA tag spec
    /// (internal `start=99` → MA tag `100`). Annotations within each type
    /// are grouped by strand into sections; types with annotations on
    /// multiple strands emit multiple sections sharing the same `name`.
    /// Empty types emit no section.
    pub fn to_ma_string(&self) -> String {
        let mut out = self.read_length.to_string();
        for p in self.ma_parts_iter() {
            // `ma_section` already starts with `;` (or is empty for an
            // empty type — which can't actually happen since
            // `to_ma_parts` would skip strands with no annotations, but
            // a type with zero annotations would still hit this).
            out.push_str(&p.ma_section);
        }
        out
    }

    /// Generate the AQ:B:C array (None if no annotations have quality).
    /// Order matches MA section order.
    pub fn to_aq_array(&self) -> Option<Vec<u8>> {
        let qualities: Vec<u8> = self.ma_parts_iter().flat_map(|p| p.aq_values).collect();

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
            .filter(|t| matches!(t.encoding, Encoding::Ma))
            .any(|t| t.annotations.iter().any(|a| a.name.is_some()));

        if !has_any_names {
            return None;
        }

        let names: Vec<String> = self.ma_parts_iter().flat_map(|p| p.an_values).collect();

        Some(names.join(","))
    }

    /// Generate all BAM tag values at once.
    ///
    /// This is the recommended method for serializing annotations to BAM tags.
    /// Returns a tuple of (MA, AQ, AN) where:
    /// - MA: The MA:Z tag string (always present; lengths inline as `start-length`)
    /// - AQ: The AQ:B:C array (qualities, None if no annotations have quality)
    /// - AN: The AN:Z tag string (names, None if no annotations have names)
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec, Encoding};
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations
    ///     .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
    ///     .add(99, 50, Strand::Forward, vec![40], None);  // 0-based
    ///
    /// let (ma, aq, an) = annotations.to_tags();
    /// assert_eq!(ma, "1000;msp+P:100-50");  // 1-based in tag
    /// assert_eq!(aq, Some(vec![40]));
    /// assert_eq!(an, None);
    /// ```
    pub fn to_tags(&self) -> (String, Option<Vec<u8>>, Option<String>) {
        let ma = self.to_ma_string();
        let aq = self.to_aq_array();
        let an = self.to_an_string();
        (ma, aq, an)
    }

    /// Write MA-family annotations (MA:Z, and optionally AQ:B:C, AN:Z) to a
    /// BAM record. Lengths are always inline in MA:Z; any stale AL tag left by
    /// the retired separate-length encoding is stripped, never written.
    ///
    /// This writes only `Encoding::Ma` annotation types. It deliberately does
    /// NOT touch the record's MM/ML tags: base-modification types
    /// (`Encoding::MmMl`) are left to [`write_mm_ml`](Self::write_mm_ml),
    /// which producers call explicitly. Leaving MM/ML untouched here means any
    /// read/edit path that doesn't change base mods preserves the record's
    /// original MM/ML bytes verbatim.
    ///
    /// # Example
    /// ```ignore
    /// use rust_htslib::bam::{self, Read};
    /// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec, Encoding};
    ///
    /// let mut record = /* from BAM file */;
    /// let mut annotations = MolecularAnnotations::from_record(&record);
    /// annotations
    ///     .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
    ///     .add(100, 50, Strand::Forward, vec![40], None);
    ///
    /// annotations.to_record(&mut record);
    /// ```
    #[cfg(feature = "htslib")]
    pub fn to_record(&self, record: &mut rust_htslib::bam::Record) {
        use rust_htslib::bam::record::Aux;

        let (ma, aq, an) = self.to_tags();

        // push_aux refuses to overwrite an existing tag, so remove first:
        // re-writing a record that already carries MA-family tags must replace
        // them, not silently no-op and leave the stale values in place. AL is
        // stripped but never re-written: lengths are inline now, so any AL on
        // the record is a leftover that must not survive the rewrite.
        record.remove_aux(b"MA").ok();
        record.remove_aux(b"AL").ok();
        record.remove_aux(b"AQ").ok();
        record.remove_aux(b"AN").ok();

        record.push_aux(b"MA", Aux::String(&ma)).ok();
        if let Some(ref aq_arr) = aq {
            record.push_aux(b"AQ", Aux::ArrayU8(aq_arr.into())).ok();
        }
        if let Some(ref an_str) = an {
            record.push_aux(b"AN", Aux::String(an_str)).ok();
        }
    }

    /// Emit `MmMl`-encoded annotation types as MM:Z + ML:B,C tags.
    ///
    /// Any pre-existing MM/ML on the record are removed first; if there are no
    /// MM/ML-encoded types, neither tag is emitted. This is a destructive,
    /// canonical re-encode intended only for callers that actually produce or
    /// modify base modifications. Read/edit paths that don't change base mods
    /// must NOT call this: leaving the record's original MM/ML bytes untouched
    /// preserves them byte-identically, including spec-legal encodings the
    /// normalized model cannot represent (e.g. grouped multi-code `C+mh`).
    #[cfg(feature = "htslib")]
    pub fn write_mm_ml(&self, record: &mut rust_htslib::bam::Record) {
        use rust_htslib::bam::record::Aux;

        // --- MM/ML assembly ---
        let forward_seq: Vec<u8> = if record.is_reverse() {
            bio::alphabets::dna::revcomp(record.seq().as_bytes())
        } else {
            record.seq().as_bytes()
        };

        // Collect per-type MM parts and pair groups with their ML byte slices.
        let mut all_groups: Vec<(crate::MmGroup, Vec<u8>)> = Vec::new();
        for t in self.annotation_types.iter().filter(|t| t.is_mm_ml()) {
            if let Some(parts) = t.to_mm_ml_parts(&forward_seq) {
                let mut cursor = 0usize;
                for g in parts.mm_groups {
                    let n = g.deltas.len();
                    let slice = parts.ml_bytes_in_order[cursor..cursor + n].to_vec();
                    cursor += n;
                    all_groups.push((g, slice));
                }
            }
        }
        // Deterministic order: by (skip_base, header).
        all_groups
            .sort_by(|(a, _), (b, _)| a.skip_base.cmp(&b.skip_base).then(a.header.cmp(&b.header)));

        record.remove_aux(b"MM").ok();
        record.remove_aux(b"ML").ok();
        if !all_groups.is_empty() {
            let mut mm_string = String::new();
            let mut ml_bytes: Vec<u8> = Vec::new();
            for (group, mls) in all_groups {
                mm_string.push_str(&group.header);
                for d in &group.deltas {
                    mm_string.push(',');
                    mm_string.push_str(&d.to_string());
                }
                mm_string.push(';');
                ml_bytes.extend(mls);
            }
            record.push_aux(b"MM", Aux::String(&mm_string)).ok();
            record
                .push_aux(b"ML", Aux::ArrayU8((&ml_bytes).into()))
                .ok();
        }
    }
}
