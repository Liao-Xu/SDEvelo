# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
# import os
# import sys
# sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SDEvelo'
copyright = '2024, Xu Liao'
author = 'Xu Liao'
release = 'https://pypi.org/project/sdevelo/'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
# import os
# exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**/.ipynb_checkpoints', '**/*-checkpoint.rst']

# # Exclude any file or directory starting with a dot
# exclude_patterns += [f for f in os.listdir('.') if f.startswith('.')]

# Configuration file for the Sphinx documentation builder.
# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
