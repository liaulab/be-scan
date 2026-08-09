# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'be_scan'
copyright = '2023, XYH, SPS, BBL'
author = 'XYH, SPS, BBL'
release = '3.0.0'

from unittest.mock import Mock
import sys

MOCK_MODULES = ['statsmodels', 'statsmodels.stats.multitest', 
                'statsmodels.nonparametric.smoothers_lowess', 
                'statsmodels.sandbox.regression.predstd', 
                'plotly.graph_objects', 'plotly.express', 'plotly.subplots', 
                'biopandas.pdb', 'plotly', 
                'plotly.express', 'plotly.subplots', 'plotly.graph_objects', 
                ]
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.duration",
    "sphinx.ext.autosectionlabel",
    "nbsphinx", # or MyST-NB
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
] 

# notebooks are committed with their outputs already computed, so the docs
# build does not need the example data or a working kernel
nbsphinx_execute = 'never'

# autodoc: keep the source order rather than alphabetising, and document
# __init__ alongside the class
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# autosectionlabel otherwise collides whenever a notebook and a module page
# use the same heading text
autosectionlabel_prefix_document = True

templates_path = ['_templates']
exclude_patterns = ['notebooks/.ipynb_checkpoints', '**.ipynb_checkpoints']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
