#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""is_notebook function"""

# %% ------------------------- Standard Imports --------------------------- %% #
import sys

# %% ============================= Experimental =========================== %% #
from IPython import get_ipython
def is_notebook():
    """[EXPERIMENTAL] Return True if running from a ipython interactive terminal or jupyter notebook."""
    # assert IPythonok, 'is_notebook() cannot check for notebook\nError: python package `IPython` not found\nmaybe install it via ``pip install IPython``' 
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return True  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

