# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'be_scan'
copyright = '2023, XYH, SPS, BBL'
author = 'XYH, SPS, BBL'
release = '1.0.0'

from unittest.mock import Mock
import sys

MOCK_MODULES = ['statsmodels', 'statsmodels.stats.multitest', 
                'statsmodels.nonparametric.smoothers_lowess', 
                'statsmodels.sandbox.regression.predstd', 
                'plotly.graph_objects', 'plotly.express', 'plotly.subplots', 
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

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
