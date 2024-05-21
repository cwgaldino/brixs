# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../'))
import sphinx_rtd_theme

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'brixs'
copyright = '2024, Carlos Galdino'
author = 'Carlos Galdino'
release = '0.9.7'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.mathjax',
              'sphinx.ext.githubpages',
              'sphinx.ext.doctest',
              'sphinx.ext.autosummary',
              'sphinx_copybutton',
              'sphinx_rtd_theme']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- sphinx_copybutton configuration ------------------------------------------
copybutton_prompt_text = ">>> "
copybutton_remove_prompts = True
copybutton_only_copy_prompt_lines = True

# -- autodoc configuration ----------------------------------------------------
autodoc_member_order = 'bysource'

# autodoc_default_flags = ['members']  # auto generate page from classes
# autosummary_generate = True          # auto generate page from classes
# toc_object_entries_show_parents = 'all'

# html_theme_options = {'collapse_navigation': False,
#                       'navigation_depth': 4,}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']
