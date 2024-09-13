# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'SDEvelo'
copyright = '2024, Xu Liao'
author = 'Xu Liao'
release = 'https://pypi.org/project/sdevelo/'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
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

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']