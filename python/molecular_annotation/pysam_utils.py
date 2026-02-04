"""Utility functions for working with pysam AlignedSegment records.

All coordinates use 0-based half-open intervals [start, end).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from molecular_annotation._molecular_annotation import MolecularAnnotations

if TYPE_CHECKING:
    import pysam

__all__ = ["from_record", "to_record", "write_to_record"]


def from_record(
    record: "pysam.AlignedSegment", parse_tags: bool = True
) -> MolecularAnnotations:
    """Read molecular annotations from a pysam AlignedSegment record.

    Extracts alignment information (aligned blocks, strand) and optionally
    parses MA/AL/AQ/AN tags.

    All coordinates are 0-based half-open [start, end).

    Args:
        record: A pysam AlignedSegment object
        parse_tags: If True (default), parse MA/AL/AQ/AN tags from the record.
            If False, only extract alignment info (for building annotations manually).

    Returns:
        MolecularAnnotations object with aligned blocks set for liftover support.

    Raises:
        KeyError: If parse_tags=True and MA tag is missing
        ValueError: If tag format is invalid
    """
    is_reverse = record.is_reverse

    if parse_tags:
        # Get MA tag (required when parsing tags)
        ma = record.get_tag("MA")

        # Get AL tag (optional - empty if inline format)
        try:
            al = list(record.get_tag("AL"))
        except KeyError:
            al = []

        # Get AQ tag (optional)
        try:
            aq = list(record.get_tag("AQ"))
        except KeyError:
            aq = None

        # Get AN tag (optional)
        try:
            an = record.get_tag("AN")
        except KeyError:
            an = None

        annot = MolecularAnnotations.from_tags(ma, al, aq=aq, an=an)
        annot.is_reverse_aligned = is_reverse
    else:
        # Create empty annotations with just read length
        annot = MolecularAnnotations(record.query_length)
        annot.is_reverse_aligned = is_reverse

    # Extract aligned blocks from CIGAR for liftover support
    # pysam's get_blocks() returns list of (ref_start, ref_end) tuples
    # We need to compute query positions from the CIGAR
    if not record.is_unmapped and record.cigartuples:
        aligned_blocks = _extract_aligned_blocks(record)
        if aligned_blocks:
            annot.set_aligned_blocks(aligned_blocks, is_reverse=is_reverse)

    return annot


def _extract_aligned_blocks(
    record: "pysam.AlignedSegment",
) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    """Extract aligned block pairs from a pysam record.

    Returns a list of ((query_start, query_end), (ref_start, ref_end)) tuples.
    All coordinates are 0-based half-open.
    """
    # CIGAR operations
    BAM_CMATCH = 0  # M
    BAM_CINS = 1  # I
    BAM_CDEL = 2  # D
    BAM_CREF_SKIP = 3  # N
    BAM_CSOFT_CLIP = 4  # S
    BAM_CHARD_CLIP = 5  # H
    BAM_CPAD = 6  # P
    BAM_CEQUAL = 7  # =
    BAM_CDIFF = 8  # X

    blocks = []
    query_pos = 0
    ref_pos = record.reference_start

    # Handle leading hard/soft clips
    for op, length in record.cigartuples:
        if op == BAM_CHARD_CLIP:
            continue  # Hard clips don't consume query
        elif op == BAM_CSOFT_CLIP:
            query_pos += length
        else:
            break

    # Process CIGAR operations
    for op, length in record.cigartuples:
        if op in (BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF):
            # Aligned block
            blocks.append((
                (query_pos, query_pos + length),
                (ref_pos, ref_pos + length)
            ))
            query_pos += length
            ref_pos += length
        elif op == BAM_CINS:
            # Insertion - consumes query only
            query_pos += length
        elif op in (BAM_CDEL, BAM_CREF_SKIP):
            # Deletion/skip - consumes reference only
            ref_pos += length
        elif op == BAM_CSOFT_CLIP:
            # Soft clip at end - consumes query only
            query_pos += length
        # Hard clips and padding don't consume either

    return blocks


def to_record(annotations: MolecularAnnotations, record: "pysam.AlignedSegment") -> None:
    """Write molecular annotations to a pysam AlignedSegment record.

    Sets MA:Z tag, and optionally AL:B:I, AQ:B:C, and AN:Z tags depending
    on the encoding format and whether quality/names are present.

    This is the inverse of `from_record()`.

    Args:
        annotations: MolecularAnnotations object
        record: A pysam AlignedSegment object

    Example:
        >>> import pysam
        >>> # Read a record
        >>> bam = pysam.AlignmentFile("input.bam", "rb")
        >>> record = next(bam)
        >>>
        >>> # Create/modify annotations
        >>> annot = from_record(record, parse_tags=False)
        >>> annot.add_annotations("msp", "+", "P", [100, 200], lengths=[50, 60], qualities=[40, 35])
        >>>
        >>> # Write back to record
        >>> to_record(annot, record)
    """
    ma, al, aq, an = annotations.to_tags()

    # Set MA tag (always)
    record.set_tag("MA", ma)

    # Set AL tag (only for separate encoding, will be non-empty)
    if al:
        record.set_tag("AL", al)

    # Set AQ tag if present
    if aq is not None:
        record.set_tag("AQ", aq)

    # Set AN tag if present
    if an is not None:
        record.set_tag("AN", an)


# Backwards compatibility alias
write_to_record = to_record
