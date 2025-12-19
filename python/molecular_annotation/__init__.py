"""Python bindings for molecular annotation tags in SAM/BAM/CRAM files."""

# Re-export from the Rust extension
from molecular_annotation._molecular_annotation import MolecularAnnotations

# Re-export pysam utilities
from molecular_annotation.pysam_utils import from_record, write_to_record

__all__ = ["MolecularAnnotations", "from_record", "write_to_record"]
