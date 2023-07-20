# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys

sys.path.insert(0, os.path.abspath("../"))

import py_ft

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
# The short X.Y version.
project = "py-ft"
copyright = "2023, Mitchell R. Vollger"
author = "Mitchell R. Vollger"
version = py_ft.__version__
# The full version, including alpha/beta/rc tags.
release = py_ft.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
    "sphinx.ext.intersphinx",
    # "edit_on_github",
    "m2r2",
]

# source_suffix = '.rst'
source_suffix = [".rst", ".md", ".css", ".png"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "alabaster"
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
# html_css_files = [
#    "static/css/rtd_dark.css",
# ]
html_logo = "../assets/img/fiber_tools_grey.png"
