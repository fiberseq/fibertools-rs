"""Tests for pyft.fiberdata."""

from pathlib import Path

import pytest

from pyft.fiberdata import Feature, Fiberdata, fetch

# Reuse the fibertools-rs repo fixture (all.bam); absent outside the repo.
_BAM = Path(__file__).parents[4] / "tests" / "data" / "all.bam"


class TestRealFixture:
    """Against a real fiber-seq BAM (all 22 reads are mapped)."""

    pytestmark = pytest.mark.skipif(
        not _BAM.exists(),
        reason="fiber-seq BAM fixture only present inside the fibertools-rs repo",
    )

    def _first(self) -> Fiberdata:
        for fd in fetch(_BAM):
            return fd
        raise AssertionError("no records in fixture")

    def test_fetch_yields_fibers(self):
        assert sum(1 for _ in fetch(_BAM)) == 22

    def test_first_fiber_metadata_and_counts(self):
        fd = self._first()
        # Known values for m54329U_210323_190418/5048829/ccs (reverse-aligned).
        assert fd.name == "m54329U_210323_190418/5048829/ccs"
        assert fd.is_mapped
        assert fd.chrom == "ptg000001l"
        assert fd.strand == "-"
        assert fd.length == 15524
        assert len(fd.m6a) == 1541
        assert len(fd.cpg) == 135
        assert len(fd.nuc) == 76
        assert len(fd.msp) == 77

    def test_features_are_typed_and_lift_to_reference(self):
        fd = self._first()
        m6a = fd.m6a
        assert all(isinstance(f, Feature) for f in m6a)
        assert all(0 <= f.start < fd.length for f in m6a)
        assert len([f for f in m6a if f.ref_start is not None]) > 0
        assert all(0 <= f.qual <= 255 for f in m6a)


class TestUnmappedRead:
    """Fiberdata works on unmapped reads (e.g. unaligned fiber-seq BAMs)."""

    def _unmapped_record(self, pysam):
        import array

        header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.0"}})
        rec = pysam.AlignedSegment(header)
        rec.query_name = "unaligned_read"
        rec.query_sequence = "ACAGAA"  # A at 0,2,4,5
        rec.flag = 4  # unmapped
        rec.set_tag("MA", "6")  # read length only, no structural annotations
        rec.set_tag("MM", "A+a,1,0,0;")  # m6A at 2,4,5
        rec.set_tag("ML", array.array("B", [200, 150, 100]))
        return rec

    def test_unmapped_read_parses_without_reference(self):
        pysam = pytest.importorskip("pysam")
        fd = Fiberdata.from_record(self._unmapped_record(pysam))

        assert fd.is_mapped is False
        assert fd.chrom is None and fd.start is None and fd.end is None
        assert fd.strand == "+"
        # base mods still decode from the molecular sequence...
        assert [f.start for f in fd.m6a] == [2, 4, 5]
        # ...but there is nothing to lift to, so reference coords are None.
        assert all(f.ref_start is None for f in fd.m6a)
