# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'UniST'
copyright = '2026, Lan Shui'
author = 'Lan Shui'

release = 'latest'
version = 'latest'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    "sphinx_design",
    "nbsphinx",
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = "_static/unist_logo.png"

# -- Options for EPUB output
epub_show_urls = 'footnote'
