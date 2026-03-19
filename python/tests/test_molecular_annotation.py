"""Tests for the molecular_annotation Python bindings."""

import pytest
from molecular_annotation import MolecularAnnotations, from_record, write_to_record


class TestMolecularAnnotations:
    """Tests for the MolecularAnnotations class."""

    def test_new(self):
        """Test creating a new empty MolecularAnnotations."""
        annotations = MolecularAnnotations(1000)
        assert annotations.read_length == 1000
        assert annotations.total_annotation_count() == 0

    def test_add_annotations_basic(self):
        """Test adding annotations without quality."""
        annotations = MolecularAnnotations(1000)
        # API uses 0-based coordinates, MA string uses 1-based per spec
        annotations.add_annotations("nuc", "+", "", starts=[100, 250], lengths=[147, 147])

        assert annotations.total_annotation_count() == 2
        assert annotations.to_ma_string() == "1000;nuc+:101-147,251-147"
        assert annotations.to_aq_array() is None

    def test_add_annotations_with_quality(self):
        """Test adding annotations with Phred quality scores."""
        annotations = MolecularAnnotations(1000)
        # API uses 0-based coordinates, MA string uses 1-based per spec
        annotations.add_annotations(
            "msp", "+", "P", starts=[100, 200], lengths=[50, 60], qualities=[40, 35]
        )

        assert annotations.total_annotation_count() == 2
        assert annotations.to_ma_string() == "1000;msp+P:101-50,201-60"
        assert annotations.to_aq_array() == [40, 35]

    def test_add_annotations_with_linear_quality(self):
        """Test adding annotations with linear quality scores."""
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations(
            "fire", ".", "Q", starts=[500], lengths=[75], qualities=[200]
        )

        assert annotations.to_ma_string() == "1000;fire.Q:501-75"
        assert annotations.to_aq_array() == [200]

    def test_add_annotations_with_names(self):
        """Test adding annotations with names."""
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations(
            "fire",
            ".",
            "Q",
            starts=[500, 700],
            lengths=[75, 80],
            qualities=[200, 180],
            names=["enhancer1", "promoter2"],
        )

        assert annotations.to_an_string() == "enhancer1,promoter2"

    def test_add_multiple_annotation_types(self):
        """Test adding multiple annotation types."""
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations(
            "msp", "+", "P", starts=[100, 200], lengths=[50, 60], qualities=[40, 35]
        )
        annotations.add_annotations("nuc", "+", "", starts=[150, 400], lengths=[147, 147])

        assert annotations.total_annotation_count() == 4
        assert annotations.to_ma_string() == "1000;msp+P:101-50,201-60;nuc+:151-147,401-147"
        assert annotations.to_aq_array() == [40, 35]

    def test_add_annotations_to_existing_type(self):
        """Test adding more annotations to an existing type."""
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations("nuc", "+", "", starts=[100], lengths=[147])
        annotations.add_annotations("nuc", "+", "", starts=[300], lengths=[147])

        assert annotations.total_annotation_count() == 2
        assert annotations.to_ma_string() == "1000;nuc+:101-147,301-147"

    def test_strand_variants(self):
        """Test different strand orientations."""
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations("fwd", "+", "", starts=[100], lengths=[50])
        annotations.add_annotations("rev", "-", "", starts=[200], lengths=[50])
        annotations.add_annotations("unk", ".", "", starts=[300], lengths=[50])

        ma = annotations.to_ma_string()
        assert "fwd+:" in ma
        assert "rev-:" in ma
        assert "unk.:" in ma

    def test_repr(self):
        """Test string representation."""
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations("msp", "+", "P", starts=[100], lengths=[50], qualities=[40])

        repr_str = repr(annotations)
        assert "MolecularAnnotations" in repr_str
        assert "read_length=1000" in repr_str
        assert "annotation_types=1" in repr_str


class TestFromTags:
    """Tests for parsing from tag values."""

    def test_from_tags_inline_format(self):
        """Test parsing inline format (start-length pairs)."""
        annotations = MolecularAnnotations.from_tags(
            "1000;msp+P:100-50,200-60", aq=[40, 35]
        )

        assert annotations.read_length == 1000
        assert annotations.total_annotation_count() == 2
        assert annotations.to_aq_array() == [40, 35]

    def test_from_tags_no_quality(self):
        """Test parsing annotations without quality scores."""
        annotations = MolecularAnnotations.from_tags("1000;nuc+:100-147,250-147")

        assert annotations.total_annotation_count() == 2
        assert annotations.to_aq_array() is None

    def test_from_tags_with_names(self):
        """Test parsing annotations with names."""
        annotations = MolecularAnnotations.from_tags(
            "1000;fire.Q:500-75,700-80", aq=[200, 180], an="enhancer1,promoter2"
        )

        assert annotations.to_an_string() == "enhancer1,promoter2"

    def test_from_tags_mixed_quality(self):
        """Test parsing with mixed quality types."""
        annotations = MolecularAnnotations.from_tags(
            "1000;msp+P:100-50,200-60;nuc+:150-147;fire.Q:500-75",
            aq=[40, 35, 200],
        )

        assert annotations.total_annotation_count() == 4
        # msp has 2 quality scores, fire has 1
        assert annotations.to_aq_array() == [40, 35, 200]

    def test_from_tags_roundtrip(self):
        """Test that parsing and regenerating produces same output."""
        original_ma = "1000;msp+P:100-50,200-60;nuc+:150-147,350-147"
        original_aq = [40, 35]

        annotations = MolecularAnnotations.from_tags(original_ma, aq=original_aq)

        assert annotations.to_ma_string() == original_ma
        assert annotations.to_aq_array() == original_aq


class TestErrorHandling:
    """Tests for error handling."""

    def test_invalid_strand(self):
        """Test that invalid strand raises error."""
        annotations = MolecularAnnotations(1000)
        with pytest.raises(ValueError, match="Invalid strand"):
            annotations.add_annotations("test", "X", "", starts=[100], lengths=[50])

    def test_invalid_quality_spec(self):
        """Test that invalid quality spec raises error."""
        annotations = MolecularAnnotations(1000)
        with pytest.raises(ValueError, match="Invalid quality spec character"):
            annotations.add_annotations("test", "+", "X", starts=[100], lengths=[50])

    def test_mismatched_starts_lengths(self):
        """Test that mismatched starts/lengths raises error."""
        annotations = MolecularAnnotations(1000)
        with pytest.raises(ValueError, match="same length"):
            annotations.add_annotations("test", "+", "", starts=[100, 200], lengths=[50])

    def test_mismatched_qualities(self):
        """Test that mismatched qualities raises error."""
        annotations = MolecularAnnotations(1000)
        with pytest.raises(ValueError, match="qualities must have length"):
            annotations.add_annotations(
                "test", "+", "P", starts=[100, 200], lengths=[50, 60], qualities=[40]
            )

    def test_mismatched_names(self):
        """Test that mismatched names raises error."""
        annotations = MolecularAnnotations(1000)
        with pytest.raises(ValueError, match="same length"):
            annotations.add_annotations(
                "test", "+", "", starts=[100, 200], lengths=[50, 60], names=["only_one"]
            )

    def test_conflicting_annotation_type(self):
        """Test that adding same type with different strand/quality raises error."""
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations("msp", "+", "P", starts=[100], lengths=[50], qualities=[40])

        with pytest.raises(ValueError, match="already exists"):
            # Same name but different strand
            annotations.add_annotations(
                "msp", "-", "P", starts=[200], lengths=[60], qualities=[35]
            )

    def test_from_tags_invalid_format(self):
        """Test that invalid MA format raises error."""
        with pytest.raises(ValueError):
            MolecularAnnotations.from_tags("invalid")

    def test_from_tags_missing_quality(self):
        """Test that missing quality for P type raises error."""
        with pytest.raises(ValueError, match="Quality type specified"):
            MolecularAnnotations.from_tags("1000;msp+P:100-50")


class TestMultiQuality:
    """Tests for multi-quality support."""

    def test_add_annotations_multi_quality(self):
        """Test adding annotations with multiple quality values per annotation."""
        annotations = MolecularAnnotations(1000)
        # PQ = 2 quality values per annotation, 2 annotations = 4 total values
        annotations.add_annotations(
            "ctcf", "+", "PQ", starts=[100, 200], lengths=[50, 60],
            qualities=[40, 255, 30, 200],
        )

        assert annotations.total_annotation_count() == 2
        assert annotations.to_ma_string() == "1000;ctcf+PQ:101-50,201-60"
        assert annotations.to_aq_array() == [40, 255, 30, 200]

    def test_add_annotations_quad_quality(self):
        """Test adding annotations with 4 quality values per annotation."""
        annotations = MolecularAnnotations(1000)
        # PQQP = 4 values per annotation, 1 annotation = 4 total
        annotations.add_annotations(
            "test", ".", "PQQP", starts=[100], lengths=[50],
            qualities=[40, 200, 180, 35],
        )

        assert annotations.to_ma_string() == "1000;test.PQQP:101-50"
        assert annotations.to_aq_array() == [40, 200, 180, 35]

    def test_multi_quality_wrong_count(self):
        """Test error when quality count doesn't match spec * annotations."""
        annotations = MolecularAnnotations(1000)
        with pytest.raises(ValueError, match="qualities must have length 4"):
            # PQ = 2 per annotation, 2 annotations needs 4, but only 2 given
            annotations.add_annotations(
                "ctcf", "+", "PQ", starts=[100, 200], lengths=[50, 60],
                qualities=[40, 30],
            )

    def test_from_tags_multi_quality(self):
        """Test parsing multi-quality from tags."""
        annotations = MolecularAnnotations.from_tags(
            "1000;ctcf+PQ:100-50,200-60", aq=[40, 255, 30, 200]
        )

        assert annotations.total_annotation_count() == 2
        assert annotations.to_aq_array() == [40, 255, 30, 200]

    def test_roundtrip_multi_quality(self):
        """Test multi-quality roundtrip through tags."""
        original_ma = "1000;ctcf+PQ:100-50,200-60;nuc+:150-147"
        original_aq = [40, 255, 30, 200]

        annotations = MolecularAnnotations.from_tags(original_ma, aq=original_aq)
        assert annotations.to_ma_string() == original_ma
        assert annotations.to_aq_array() == original_aq

    def test_iter_full_multi_quality(self):
        """Test iter_full returns correct quality lists for multi-quality."""
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations(
            "ctcf", "+", "PQ", starts=[100, 200], lengths=[50, 60],
            qualities=[40, 255, 30, 200],
        )
        annotations.add_annotations("nuc", "+", "", starts=[150], lengths=[147])

        results = annotations.iter_full()
        assert len(results) == 3

        # First ctcf annotation: qualities [40, 255]
        type_name, strand, qs, qs_start, qe, fs, fe, rs, re, quals, name = results[0]
        assert type_name == "ctcf"
        assert qs == "PQ"
        assert quals == [40, 255]

        # Second ctcf annotation: qualities [30, 200]
        _, _, _, _, _, _, _, _, _, quals2, _ = results[1]
        assert quals2 == [30, 200]

        # nuc annotation: no qualities
        _, _, qs3, _, _, _, _, _, _, quals3, _ = results[2]
        assert qs3 == ""
        assert quals3 == []

    def test_iter_type_multi_quality(self):
        """Test iter_type returns correct quality lists for multi-quality."""
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations(
            "ctcf", "+", "PQ", starts=[100, 200], lengths=[50, 60],
            qualities=[40, 255, 30, 200],
        )

        items = annotations.iter_type("ctcf")
        assert items is not None
        assert len(items) == 2

        qs, qe, fs, fe, rs, re, quals, name = items[0]
        assert quals == [40, 255]

        _, _, _, _, _, _, quals2, _ = items[1]
        assert quals2 == [30, 200]

    def test_mixed_single_and_multi_quality(self):
        """Test mixing single-quality and multi-quality annotation types."""
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations(
            "msp", "+", "P", starts=[100], lengths=[50], qualities=[40]
        )
        annotations.add_annotations(
            "ctcf", "+", "PQ", starts=[200], lengths=[60], qualities=[30, 200]
        )
        annotations.add_annotations("nuc", "+", "", starts=[300], lengths=[147])

        assert annotations.total_annotation_count() == 3
        # AQ array: msp's 1 value, then ctcf's 2 values
        assert annotations.to_aq_array() == [40, 30, 200]


class TestPysamIntegration:
    """Tests for pysam integration."""

    def test_write_and_read_record(self):
        """Test writing to and reading from a pysam record."""
        pysam = pytest.importorskip("pysam")

        # Create a minimal BAM record
        header = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 1000000}]}
        )
        record = pysam.AlignedSegment(header)
        record.query_name = "test_read"
        record.query_sequence = "A" * 1000
        record.flag = 0
        record.reference_id = 0
        record.reference_start = 100
        record.cigar = [(0, 1000)]
        record.query_qualities = pysam.qualitystring_to_array("I" * 1000)

        # Create and write annotations
        annotations = MolecularAnnotations(1000)
        annotations.add_annotations(
            "msp", "+", "P", starts=[100, 200], lengths=[50, 60], qualities=[40, 35]
        )
        annotations.add_annotations("nuc", "+", "", starts=[150, 400], lengths=[147, 147])

        write_to_record(annotations, record)

        # Verify tags were written (MA tag uses 1-based coordinates per spec)
        assert record.get_tag("MA") == "1000;msp+P:101-50,201-60;nuc+:151-147,401-147"
        assert list(record.get_tag("AQ")) == [40, 35]

        # Read them back
        annotations2 = from_record(record)

        assert annotations2.read_length == 1000
        assert annotations2.total_annotation_count() == 4
        assert annotations2.to_ma_string() == annotations.to_ma_string()
        assert annotations2.to_aq_array() == annotations.to_aq_array()

    def test_from_record_missing_ma_tag(self):
        """Test that missing MA tag raises KeyError."""
        pysam = pytest.importorskip("pysam")

        header = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 1000000}]}
        )
        record = pysam.AlignedSegment(header)
        record.query_name = "test_read"
        record.query_sequence = "ACGT"
        record.flag = 0
        record.reference_id = 0
        record.reference_start = 0
        record.cigar = [(0, 4)]

        with pytest.raises(KeyError):
            from_record(record)
