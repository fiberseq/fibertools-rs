"""pyft: an inspection and visualization toolkit for fiber-seq BAMs.

Reads fiber-seq annotations (nuc/msp/fire and m6A/5mC base mods) via pysam and
the molecular_annotation library. The plotting helpers live in `pyft.plot` and
`pyft.utils` and are imported on demand (they pull in pandas/altair).
"""

from importlib.metadata import PackageNotFoundError, version

from .fiberdata import Feature, Fiberdata, fetch

try:
    __version__ = version("pyft")
except PackageNotFoundError:  # not installed (e.g. running from a source tree)
    __version__ = "0.0.0"

__all__ = ["Feature", "Fiberdata", "fetch", "__version__"]
