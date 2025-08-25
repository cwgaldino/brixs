#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quality-of-life functions for handling hdf5 files"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import numpy as np
import h5py

# %% ------------------------------- tree --------------------------------- %% #
def tree(filepath, initial=None):
    """Returns a text with the structure of a hdf5 file

    Args:
        filepath (str or path): file to hdf5 file
        initial (int): initial hdf5 path

    Returns:
        str
    """
    with h5py.File(Path(filepath), 'r') as f:
        if initial is None:
            text = _tree(f)
        else:
            text = _tree(f[initial])
    return text

def _tree(f, pre=''):
    """Returns a text with the structure of a hdf5 file

    Usage:
        >>> import h5py
        >>> 
        >>> with h5py.File(<filepath>, 'r') as f:
        >>>     text = tree(f)
        >>>
        >>> with h5py.File(<filepath>, 'r') as f:
        >>>     text = tree(f[entry])
    
    Args:
        f (h5 object): h5 file object
        pre (str): For internal use only. This function needs to run recursively

    Returns:
        str
    """
    text = ''
    items = len(f)
    for key, f in f.items():
        items -= 1
        if items == 0:
            if type(f) == h5py._hl.group.Group:
                text += pre + '└── ' + key + '\n'
                text += _tree(f, pre + '    ')
            else:
                try:
                    text += pre + '└── ' + key + ' (%d)' % len(f) + '\n'
                except TypeError:
                    text += pre + '├── ' + key + ' (scalar)' + '\n'
        else:
            if type(f) == h5py._hl.group.Group:
                text += pre + '├── ' + key + '\n'
                text += _tree(f, pre+'│   ')
            else:
                try:
                    text += pre + '├── ' + key + ' (%d)' % len(f) + '\n'
                except TypeError:
                    text += pre + '├── ' + key + ' (scalar)' + '\n'
    return text
# %%

# %% ----------------------------- metadata ------------------------------- %% #
def sort_metadata(f, attrs_dict, verbose=True):
    """get and pre-process metadata from hdf5 file

    Args:
        attrs_dict (dict): nested dictionary indicating the type of pre-processing, 
            attr name, and hdf5 path. Available pre-processing options are:
                
                0) `ignore`: metadata isn't saved (entry is ignored). Useful 
                    only so you can set attrs_dict as a full list of attrs, without 
                    having to left attrs out. attrs_dict can than be used elsewhere 
                    (i.e., building scan lists).
                1) `raw`: no processing
                2) `string`: turn into type `str` 
                3) `bool`: turn into type `bool`
                4) `mean, `min`, `max`, `sigma`: Calculate the mean, min, max, and standard deviation of numeric iterables
                5) `round`, `round1`, `round2`., ... : round number to 0, 1, 2, ... decimal points
                6) `first_column_mean, `first_column_min`, `first_column_max`, `first_column_sigma`: same as `mean, `min`, `max`, `sigma`, but for multi-column iterables
                7) `second_column_mean, `second_column_min`, `second_column_max`, `second_column_sigma`: same as `mean, `min`, `max`, `sigma`, but for multi-column iterables
                
        f (str, Path, or hdf5 object): filepath to hdf5 file or hdf5 object
        verbose (bool, optional): if True, if attr value cannot be extracted from 
            hdf5 file it will print a warning. Default is True.

    Usage: 
        >>> attrs = {'raw':   {'A': 'External/b316a-o01/dia/tco-02/Temperature',
        >>>                    'B': 'External/beamline_energy/position'},
        >>>         'string': {'C': 'startmetadata/Sample',
        >>>                    'D': 'endmetadata/Sample}}
        >>> d = sort_metadata(<filepath>, attrs)

    Returns:
        dict[attr] : value
    """
    if isinstance(f, str) or isinstance(f, Path):
        with h5py.File(Path(f), 'r') as f:
            return _sort_metadata(f=f, attrs_dict=attrs_dict, verbose=verbose)
    else:
        return _sort_metadata(f=f, attrs_dict=attrs_dict, verbose=verbose)

def _sort_metadata(f, attrs_dict, verbose):
    # check unique names
    values = {}
    for _type in attrs_dict:
        if _type != 'ignore': 
            for name in attrs_dict[_type]:
                if name in values:
                    raise KeyError(f"name `{name}` is duplicated in VERITAS attrs list. Names must be unique")
                values[name] = None

    # get attr values
    for _type in attrs_dict:
        for name in attrs_dict[_type]:
            address = attrs_dict[_type][name]
            try: 
                if _type == 'ignore':   
                    pass
                elif _type == 'raw':   
                    values[name] = f[address][()]
                elif _type == 'string':   
                    values[name] = f[address][()].decode("utf-8")
                elif _type == 'bool':     
                    values[name] = bool(f[address][()])

                elif _type == 'mean':   
                    values[name] = np.mean(f[address][()])
                elif _type == 'min':   
                    values[name] = np.min(f[address][()])
                elif _type == 'max':   
                    values[name] = np.max(f[address][()])
                elif _type == 'sigma':   
                    values[name] = np.std(f[address][()])

                elif _type.startswith('round'):
                    n = int(_type.split('round')[-1])
                    values[name] = np.round(f[address][()], n)

                elif _type == 'first_column_mean':   
                    values[name] = np.mean(f[address][:, 0])
                elif _type == 'first_column_min':   
                    values[name] = np.min(f[address][:, 0])
                elif _type == 'first_column_max':   
                    values[name] = np.max(f[address][:, 0])
                elif _type == 'first_column_sigma':   
                    values[name] = np.std(f[address][:, 0])
                    
                elif _type == 'second_column_mean':   
                    values[name] = np.mean(f[address][:, 1])
                elif _type == 'second_column_min':   
                    values[name] = np.min(f[address][:, 1])
                elif _type == 'second_column_max':   
                    values[name] = np.max(f[address][:, 1])
                elif _type == 'second_column_sigma':   
                    values[name] = np.std(f[address][:, 1])

            except Exception as e:
                if verbose: 
                    print(f'get attr `{name}` error: {e}')
    return values
# %%  

