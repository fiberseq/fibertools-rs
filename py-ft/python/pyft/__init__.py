"""pyft: an inspection and visualization toolkit for fiber-seq BAMs.

Reads fiber-seq annotations (nuc/msp/fire and m6A/5mC base mods) via pysam and
the molecular_annotation library. The plotting helpers live in `pyft.plot` and
`pyft.utils` and are imported on demand (they pull in pandas/altair).
"""

from .fiberdata import Feature, Fiberdata, fetch

__all__ = ["Feature", "Fiberdata", "fetch"]
