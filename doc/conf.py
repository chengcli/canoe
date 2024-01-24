# Configuration file for the Sphinx documentation builder.

import os, sys

# -- Project information

project = 'pycanoe'
copyright = '2024, Canoe'
author = 'Cheng Li'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'myst_parser',
    'breathe',
]
breathe_projects = {"Canoe": "doxygen/xml"}
breathe_default_project = "Canoe"

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# -- Options for EPUB output
epub_show_urls = 'footnote'
