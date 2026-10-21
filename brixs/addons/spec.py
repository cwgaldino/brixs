#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quality-of-life functions for saving spectra as spec files"""

# %% ------------------------- Standard Imports -------------------------- %% #
from pathlib import Path
import numpy as np
import datetime


# %% -------------------------- brixs Imports ---------------------------- %% #
import brixs as br
# %%

# %% =========================== save spectra ============================ %% #
def _save_all_as_specfile(self, filepath):
    """[EXPERIMENTAL] Save spectra as specfile that pymca can read
    
    Args:
        filepath (str or path): filepath to save 
    
    Returns:
        None
    """
    ######################################
    # check if filepath points to a file #
    ######################################
    if not isinstance(filepath, (str, Path)):
        raise TypeError(f'filepath must be a string or Path object, not {type(filepath).__name__}')
    filepath = Path(filepath)
    if not filepath.parent.exists():
        raise FileNotFoundError(f'filepath parent folder does not exist: {filepath.parent}')
    if filepath.exists() and not filepath.is_file():
        raise IsADirectoryError(f'filepath must point to a file, not a folder: {filepath}')

    #######################
    # check x is the same #
    #######################
    try:
        _x = self.check_same_x()
    except ValueError:
        raise ValueError('Cannot save spectra in one file. x axis are different.\nMaybe try interpolating the x axis (Spectra.interp()) or use Spectra.save() to save spectra in multiple files.')

    # initialize file
    text = []
    text.append(f'#F {filepath.name}')
    text.append(f"#D {datetime.datetime.now().ctime()}")

    # organize data
    for i, _s in enumerate(self):
        text.append(f'')
        if hasattr(_s, 'label'):
            text.append(f'#S {i+1}  {_s.label}')
        else:
            text.append(f'#S {i+1}  no label')
        text.append(f'#N 2')
        text.append(f'#L A  B')
        for x, y in zip(_s.x, _s.y):
            text.append(f'{x} {y}')
    text.append(f'')

    # save to file
    text = '\n'.join(text)
    br.save_text(text, filepath, newline='\n')

    return 
br.Spectra.save_all_as_specfile = _save_all_as_specfile  
# %%