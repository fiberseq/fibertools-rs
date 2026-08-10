"""A thin fiber-seq read model over pysam + the molecular_annotation library.

`Fiberdata` wraps a single BAM record's molecular annotations and exposes the
fiber-seq vocabulary (m6a, 5mC, nuc, msp, fire) on top of the generic
`molecular_annotation` types. 

Coordinates are 0-based half-open. `Feature.start/end` are in molecular
(forward) orientation; `Feature.ref_start/ref_end` are lifted to the reference
(None where the annotation falls outside an aligned block).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

import pysam
from molecular_annotation import from_record as _ma_from_record

if TYPE_CHECKING:
    from collections.abc import Iterator

# MSP precision at or above which an MSP is also a FIRE element. Only used when
# the BAM does not already carry a distinct "fire" annotation type.
FIRE_MIN_QUAL = 200


@dataclass(frozen=True)
class Feature:
    """A single annotation. Coordinates are 0-based half-open."""

    start: int
    end: int
    ref_start: Optional[int]
    ref_end: Optional[int]
    qual: Optional[int]


class Fiberdata:
    """A single fiber-seq read: record metadata plus typed annotation accessors."""

    def __init__(self, record: "pysam.AlignedSegment", annotations):
        self._rec = record
        self._annot = annotations

    @classmethod
    def from_record(cls, record: "pysam.AlignedSegment") -> "Fiberdata":
        return cls(record, _ma_from_record(record))

    # --- read metadata, straight off the pysam record ---

    @property
    def name(self) -> str:
        return self._rec.query_name

    @property
    def strand(self) -> str:
        return "-" if self._rec.is_reverse else "+"

    @property
    def is_mapped(self) -> bool:
        return not self._rec.is_unmapped

    @property
    def chrom(self) -> Optional[str]:
        return None if self._rec.is_unmapped else self._rec.reference_name

    @property
    def start(self) -> Optional[int]:
        return None if self._rec.is_unmapped else self._rec.reference_start

    @property
    def end(self) -> Optional[int]:
        return None if self._rec.is_unmapped else self._rec.reference_end

    @property
    def length(self) -> int:
        return self._annot.read_length

    def tag(self, name: str, default=None):
        try:
            return self._rec.get_tag(name)
        except KeyError:
            return default

    # --- typed accessors: fiber-seq names -> molecular_annotation type strings ---

    def _features(self, type_name: str) -> list[Feature]:
        items = self._annot.iter_type(type_name)
        if not items:
            return []
        # iter_type tuple: (q_start, q_end, f_start, f_end, ref_start, ref_end, quals, name)
        return [
            Feature(f_start, f_end, ref_start, ref_end, quals[0] if quals else None)
            for _qs, _qe, f_start, f_end, ref_start, ref_end, quals, _name in items
        ]

    @property
    def m6a(self) -> list[Feature]:
        return self._features("a")

    @property
    def cpg(self) -> list[Feature]:  # 5mC
        return self._features("m")

    @property
    def nuc(self) -> list[Feature]:
        return self._features("nuc")

    @property
    def msp(self) -> list[Feature]:
        return self._features("msp")

    @property
    def fire(self) -> list[Feature]:
        # FIRE is either its own annotation type, or the high-precision MSP subset.
        if "fire" in self._annot.annotation_type_names():
            return self._features("fire")
        return [f for f in self.msp if f.qual is not None and f.qual >= FIRE_MIN_QUAL]

    def __repr__(self) -> str:
        return (
            f"Fiberdata({self.name} {self.chrom}:{self.start}-{self.end} {self.strand} "
            f"m6a={len(self.m6a)} cpg={len(self.cpg)} nuc={len(self.nuc)} msp={len(self.msp)})"
        )


def fetch(bam_path, region=None) -> "Iterator[Fiberdata]":
    """Yield `Fiberdata` over a BAM file.

    Without `region`, iterates every read in file order -- mapped **and**
    unmapped -- so an unaligned fiber-seq BAM (raw PacBio/ONT output, before
    alignment) works too; filter with `Fiberdata.is_mapped` if you only want
    one. With `region`, uses pysam's indexed fetch, which returns the mapped
    reads overlapping that region.

    Args:
        bam_path: Path to a BAM/CRAM file (indexed if `region` is given).
        region: Optional (chrom, start, end) tuple passed to `pysam.fetch`.
    """
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        source = bam.fetch(*region) if region is not None else bam.fetch(until_eof=True)
        for record in source:
            yield Fiberdata.from_record(record)
