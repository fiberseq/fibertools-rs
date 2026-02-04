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
        assert annotations.to_al_array() == [147, 147]
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
        assert annotations.to_al_array() == [50, 60]
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
        assert annotations.to_al_array() == [50, 60, 147, 147]
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
            "1000;msp+P:100-50,200-60", [], aq=[40, 35]
        )

        assert annotations.read_length == 1000
        assert annotations.total_annotation_count() == 2
        assert annotations.to_al_array() == [50, 60]
        assert annotations.to_aq_array() == [40, 35]

    def test_from_tags_separate_format(self):
        """Test parsing separate format (lengths in AL array)."""
        annotations = MolecularAnnotations.from_tags(
            "1000;msp+P:100,200", [50, 60], aq=[40, 35]
        )

        assert annotations.read_length == 1000
        assert annotations.total_annotation_count() == 2
        assert annotations.to_al_array() == [50, 60]

    def test_from_tags_no_quality(self):
        """Test parsing annotations without quality scores."""
        annotations = MolecularAnnotations.from_tags("1000;nuc+:100-147,250-147", [])

        assert annotations.total_annotation_count() == 2
        assert annotations.to_aq_array() is None

    def test_from_tags_with_names(self):
        """Test parsing annotations with names."""
        annotations = MolecularAnnotations.from_tags(
            "1000;fire.Q:500-75,700-80", [], aq=[200, 180], an="enhancer1,promoter2"
        )

        assert annotations.to_an_string() == "enhancer1,promoter2"

    def test_from_tags_mixed_quality(self):
        """Test parsing with mixed quality types."""
        annotations = MolecularAnnotations.from_tags(
            "1000;msp+P:100-50,200-60;nuc+:150-147;fire.Q:500-75",
            [],
            aq=[40, 35, 200],
        )

        assert annotations.total_annotation_count() == 4
        # msp has 2 quality scores, fire has 1
        assert annotations.to_aq_array() == [40, 35, 200]

    def test_from_tags_roundtrip(self):
        """Test that parsing and regenerating produces same output."""
        original_ma = "1000;msp+P:100-50,200-60;nuc+:150-147,350-147"
        original_aq = [40, 35]

        annotations = MolecularAnnotations.from_tags(original_ma, [], aq=original_aq)

        assert annotations.to_ma_string() == original_ma
        assert annotations.to_aq_array() == original_aq


class TestErrorHandling:
    """Tests for error handling."""

    def test_invalid_strand(self):
        """Test that invalid strand raises error."""
        annotations = MolecularAnnotations(1000)
        with pytest.raises(ValueError, match="Invalid strand"):
            annotations.add_annotations("test", "X", "", starts=[100], lengths=[50])

    def test_invalid_quality_type(self):
        """Test that invalid quality type raises error."""
        annotations = MolecularAnnotations(1000)
        with pytest.raises(ValueError, match="Invalid quality_type"):
            annotations.add_annotations("test", "+", "X", starts=[100], lengths=[50])

    def test_mismatched_starts_lengths(self):
        """Test that mismatched starts/lengths raises error."""
        annotations = MolecularAnnotations(1000)
        with pytest.raises(ValueError, match="same length"):
            annotations.add_annotations("test", "+", "", starts=[100, 200], lengths=[50])

    def test_mismatched_qualities(self):
        """Test that mismatched qualities raises error."""
        annotations = MolecularAnnotations(1000)
        with pytest.raises(ValueError, match="same length"):
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
            MolecularAnnotations.from_tags("invalid", [])

    def test_from_tags_missing_quality(self):
        """Test that missing quality for P type raises error."""
        with pytest.raises(ValueError, match="Quality type specified"):
            MolecularAnnotations.from_tags("1000;msp+P:100-50", [])


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
