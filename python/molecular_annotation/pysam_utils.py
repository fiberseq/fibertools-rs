"""Utility functions for working with pysam AlignedSegment records."""

from molecular_annotation._molecular_annotation import MolecularAnnotations


def from_record(record) -> MolecularAnnotations:
    """Read molecular annotations from a pysam AlignedSegment record.

    Extracts MA:Z, AL:B:I (optional), AQ:B:C (optional), and AN:Z (optional) tags.

    Args:
        record: A pysam AlignedSegment object

    Returns:
        MolecularAnnotations object

    Raises:
        KeyError: If MA tag is missing
        ValueError: If tag format is invalid
    """
    # Get MA tag (required)
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

    return MolecularAnnotations.from_tags(ma, al, aq=aq, an=an)


def write_to_record(annotations: MolecularAnnotations, record) -> None:
    """Write molecular annotations to a pysam AlignedSegment record.

    Sets MA:Z tag, and optionally AQ:B:C and AN:Z tags.
    Note: AL tag is not written when using inline format (the default).

    Args:
        annotations: MolecularAnnotations object
        record: A pysam AlignedSegment object
    """
    # Set MA tag
    record.set_tag("MA", annotations.to_ma_string())

    # Set AQ tag if present
    aq = annotations.to_aq_array()
    if aq is not None:
        record.set_tag("AQ", aq)

    # Set AN tag if present
    an = annotations.to_an_string()
    if an is not None:
        record.set_tag("AN", an)
