# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import re
import sys
import types

# The docs build must not depend on compiling the pyft Rust extension:
# Read the Docs was failing in `pip install ./py-ft` because the native
# (maturin/PyO3) build no longer compiles there. Instead, make the
# pure-python sources importable directly and stub out the compiled
# `pyft.pyft` submodule so `import pyft` works without the Rust build.
sys.path.insert(0, os.path.abspath("../python"))  # provides the `pyft` package
sys.path.insert(0, os.path.abspath("../python/pyft"))  # provides `utils` for api.rst

_rust_stub = types.ModuleType("pyft.pyft")
_rust_stub.__doc__ = (
    "Compiled Rust extension of pyft (not available during the docs build)."
)
# `pyft/__init__.py` does `from .pyft import *` and then reads `pyft.__doc__`,
# so the star-import must bind the name `pyft` to this stub.
_rust_stub.pyft = _rust_stub
_rust_stub.__all__ = ["pyft"]
sys.modules["pyft.pyft"] = _rust_stub

# Heavy runtime dependencies of pyft.utils/pyft.plot are not installed for the
# docs build; mock them for autodoc.
autodoc_mock_imports = [
    "pandas",
    "altair",
    "polars",
    "numpy",
    "vegafusion",
    "tqdm",
]

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
# The short X.Y version.
project = "pyft"
copyright = "2023, Mitchell R. Vollger"
author = "Mitchell R. Vollger"

# Read the version from py-ft/Cargo.toml instead of importing the built package.
with open(os.path.join(os.path.dirname(__file__), "..", "Cargo.toml")) as _fh:
    version = re.search(r'^version\s*=\s*"([^"]+)"', _fh.read(), re.M).group(1)
release = version


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    # "sphinx_autodoc_typehints",
    "sphinx.ext.viewcode",
    # "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
    "sphinx.ext.intersphinx",
    # "edit_on_github",
    "nbsphinx",
]

source_suffix = [".rst"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The vignette notebooks ship with saved outputs; never re-execute them during
# the docs build (the compiled pyft extension is not available here).
nbsphinx_execute = "never"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "alabaster"
html_theme = "sphinx_rtd_theme"
# html_permalinks_icon = '<span>#</span>'
# html_theme = 'sphinxawesome_theme'
html_static_path = ["_static"]
# html_css_files = [
# "css/rtd_dark.css",
# ]
html_logo = "_static/img/fiber_tools_grey.png"


# other options
autodoc_member_order = "bysource"
