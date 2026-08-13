"""Utility functions for working with pysam AlignedSegment records.

All coordinates use 0-based half-open intervals [start, end).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from molecular_annotation._molecular_annotation import MolecularAnnotations

if TYPE_CHECKING:
    import pysam

__all__ = ["from_record", "to_record", "write_to_record"]


_COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def _revcomp(seq: bytes) -> bytes:
    """Reverse-complement an ASCII nucleotide byte string."""
    return seq.translate(_COMPLEMENT)[::-1]


def _get_tag_any(record: "pysam.AlignedSegment", *names: str):
    """Return the first present tag among `names`, or None.

    MM/ML are written as "MM"/"ML" per the current SAM spec, but older files
    use the lowercase "Mm"/"Ml"; accept either.
    """
    for name in names:
        try:
            return record.get_tag(name)
        except KeyError:
            continue
    return None


def _ma_family_tags(record: "pysam.AlignedSegment"):
    """Read the MA tag family (Ma, Aq, An), resolving the spelling ONCE.

    Mirrors the Rust `ma_family_tags`: the spelling is chosen from the main
    tag alone, type-checked (string-typed); uppercase wins when both are
    present (a dual-spelled record can only come from an uppercase-only
    0.10-0.12 tool editing a Ma-spelled file, so its family is the fresher
    write). The sibling Aq/An are then read ONLY in the winning spelling —
    per-tag fallback would pair a fresh main tag with stale siblings from
    the other spelling. A wrong-typed main tag (e.g. a foreign MA:i) never
    wins. Returns (ma, aq, an) or None when neither spelling carries a
    string-typed main tag.
    """

    def _string_tag(name: str):
        try:
            v = record.get_tag(name)
        except KeyError:
            return None
        return v if isinstance(v, str) else None

    if _string_tag("MA") is not None:
        ma_t, aq_t, an_t = "MA", "AQ", "AN"
    elif _string_tag("Ma") is not None:
        ma_t, aq_t, an_t = "Ma", "Aq", "An"
    else:
        return None
    ma = record.get_tag(ma_t)
    try:
        aq = list(record.get_tag(aq_t))
    except (KeyError, TypeError):
        aq = None
    an = _string_tag(an_t)
    return ma, aq, an


def _parse_mm_ml_into(
    annot: MolecularAnnotations, record: "pysam.AlignedSegment", is_reverse: bool
) -> None:
    """Decode the record's MM/ML tags into `annot` (no-op if MM is absent).

    MM is delta-skip encoded against the read in its original (forward)
    molecular orientation, so a reverse-aligned record's stored sequence is
    reverse-complemented before decoding. One annotation type is created per
    mod code ("a" = m6A, "m" = 5mC, etc.).
    """
    mm = _get_tag_any(record, "MM", "Mm")
    if not mm:
        return
    ml = _get_tag_any(record, "ML", "Ml")
    ml = list(ml) if ml is not None else []

    seq = record.query_sequence
    if seq is None:
        return  # nothing to delta-decode against
    forward_seq = seq.encode("ascii")
    if is_reverse:
        forward_seq = _revcomp(forward_seq)

    annot.parse_mm_ml(mm, ml, forward_seq)


def from_record(
    record: "pysam.AlignedSegment", parse_tags: bool = True
) -> MolecularAnnotations:
    """Read molecular annotations from a pysam AlignedSegment record.

    Extracts alignment information (aligned blocks, strand) and optionally
    parses MA/AQ/AN tags (structural annotations) and MM/ML tags (base
    modifications such as m6A and 5mC).

    All coordinates are 0-based half-open [start, end).

    Args:
        record: A pysam AlignedSegment object
        parse_tags: If True (default), parse MA/AQ/AN tags from the record.
            If False, only extract alignment info (for building annotations manually).

    Returns:
        MolecularAnnotations object with aligned blocks set for liftover support.

    Raises:
        KeyError: If parse_tags=True and MA tag is missing
        ValueError: If tag format is invalid
    """
    is_reverse = record.is_reverse

    if parse_tags:
        # Resolve the Ma/Aq/An family atomically (required when parsing
        # tags); see _ma_family_tags for the spelling-precedence rules.
        family = _ma_family_tags(record)
        if family is None:
            raise KeyError("record has no Ma/MA tag")
        ma, aq, an = family

        annot = MolecularAnnotations.from_tags(ma, aq=aq, an=an)
        annot.is_reverse_aligned = is_reverse

        # Base modifications (m6A/5mC/…) live in MM/ML, not the MA tag set.
        _parse_mm_ml_into(annot, record, is_reverse)
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

    # Process CIGAR operations. Leading soft clips are consumed by the
    # BAM_CSOFT_CLIP branch below exactly once — a separate pre-scan would
    # double-count them and shift every query coordinate right by the clip
    # length.
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
            # Soft clip (leading or trailing) - consumes query only
            query_pos += length
        # Hard clips and padding don't consume either

    return blocks


def to_record(annotations: MolecularAnnotations, record: "pysam.AlignedSegment") -> None:
    """Write molecular annotations to a pysam AlignedSegment record.

    Sets the Ma:Z tag, and optionally Aq:B:C and An:Z tags depending
    on whether quality/names are present (both spellings are removed
    first; see samtools/hts-specs#862 for the canonical spelling).

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
    ma, aq, an = annotations.to_tags()

    # Write the canonical Ma/Aq/An spellings (samtools/hts-specs#862),
    # removing both spellings first so a rewrite never leaves a stale copy
    # under the other casing.
    for tag in ("MA", "Ma", "AL", "Al", "AQ", "Aq", "AN", "An"):
        if record.has_tag(tag):
            record.set_tag(tag, None)

    # Set the Ma tag (always)
    record.set_tag("Ma", ma)

    # Set the Aq tag if present
    if aq is not None:
        record.set_tag("Aq", aq)

    # Set the An tag if present
    if an is not None:
        record.set_tag("An", an)


# Backwards compatibility alias
write_to_record = to_record
