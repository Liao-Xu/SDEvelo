# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))
# Add to the list of paths where Sphinx looks for source files
sys.path.insert(0, os.path.abspath('tutorials'))
# -- Project information -----------------------------------------------------

project = 'SDEvelo'
copyright = '2024, Xu Liao'
author = 'Xu Liao'
release = 'https://pypi.org/project/sdevelo/'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'nbsphinx',
    'sphinx.ext.napoleon'
]

templates_path = ['_templates']

exclude_patterns = [
    '_build',
    'Thumbs.db',
    '.DS_Store',
    '.ipynb_checkpoints',
    '**/.ipynb_checkpoints',
    '**/.ipynb_checkpoints/**',
    '*-checkpoint.rst',
    '**/*-checkpoint.rst',
    '.*',
    '**/.*',
]
nbsphinx_allow_errors = True
# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Add this line to include your custom CSS
html_css_files = [
    'custom.css',
]