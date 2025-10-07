#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quality-of-life functions for handling hdf5 files"""

# %% ------------------------- Standard Imports -------------------------- %% #
from pathlib import Path
import numpy as np

# %% ========================== Special Imports ========================== %% #
try:
    import h5py
except:
    pass
# %%

# %% =============================== tree ================================ %% #
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

# %% ============================= metadata ============================== %% #
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
                    raise KeyError(f"name `{name}` is duplicated in attrs list. Names must be unique")
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

# %% =========================== save spectra ============================ %% #
# 	Check if the x range for each spectrum is correct before creating HDF file. 
# def save_all_as_nexus_map(self, filepath, extra={}, verbose=False):
#     """[EXPERIMENTAL] Saves an HDF5 file in NeXus NXdata format containing spectra.

#     Args:
#         filepath (str or path): filepath to save files
    
#     Returns:
#         None
#     """    
#     ######################################
#     # check if filepath points to a file #
#     ######################################
#     filepath = Path(filepath)
#     assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
#     if filepath.exists():
#         assert filepath.is_file(), 'filepath must point to a file'

#     #######################
#     # check x is the same #
#     #######################
#     try:
#         _x = self.check_same_x()
#     except ValueError:
#         raise ValueError('Cannot save spectra in one file. x axis are different.\nMaybe try interpolating the x axis (Spectra.interp()) or use Spectra.save() to save spectra in multiple files.')
        

#     # --- Checking spectra consistency ---
#     # E_inc_all = np.array([float(_s.Energy) for _s in self], dtype=np.float64)  # X axis
#     # E_loss_ref = np.asarray(self[0].x, dtype=np.float64)                       # Y axis (reference)

#     # Sort spectra by incident energy
#     # order = np.argsort(E_inc_all)
#     # E_inc = E_inc_all[order]

#     # Stack signals according to sorted incident energies
#     _, ys = ss._gather_ys()
#     # Z = np.stack([np.asarray(self[i].y, dtype=np.float32) for i in order], axis=0)

#     # check extra attrs
#     attrs = {}
#     for attr in extra:
#         assert any([hasattr(_s, attr)==False for _s in self])==False, f'Not all spectrum have extra attr `{attr}`'
#         _temp = [_s.__getattribute__(attr) for _s in self]
#         assert sum([br.is_number(x) for x in _temp]) == len(_temp), f'extra attrs must be numbers. Extra attr `{attr}` is not a number for all Spectra'
#         attrs[attr] = _temp

#     # --- Create NeXus data (NXdata) ---
#     with h5py.File(filepath, 'w') as f:
#         f.attrs['default'] = 'entry'

#         nxentry = f.create_group('entry')
#         nxentry.attrs['NX_class'] = 'NXentry'
#         nxentry.attrs['default']  = 'data'

#         nxdata = nxentry.create_group('data')
#         nxdata.attrs['NX_class'] = 'NXdata'

#         # Axes datasets
#         for attr in attrs:
#             print(attrs[attr])
#             a = np.array([float(_) for _ in attrs[attr]], dtype=np.float64)
#             _temp = nxdata.create_dataset(attr, data=a)        # X axis
#             _temp.attrs['units']     = 'units'
#             _temp.attrs['long_name'] = attr
#         x = nxdata.create_dataset('x', data=_x)   # Y axis
#         x.attrs['units']     = 'units'
#         x.attrs['long_name'] = 'x'

#         # Signal dataset
#         ds = nxdata.create_dataset('Spectra', data=ys.transpose(), chunks=True, compression='gzip', compression_opts=4)
#         ds.attrs['Energy_Incident_indices'] = np.array([0])
#         ds.attrs['Energy_Loss_indices']     = np.array([1])
#         ds.attrs['interpretation'] = 'image'       # 2D image/matrix
#         ds.attrs['long_name'] = 'y'           # Label for the signal

#         # --- Essential NXdata attributes ---
#         nxdata.attrs['signal'] = 'Spectra'
#         nxdata.attrs['axes'] = list(attrs.keys()) + ['x', ]
#         print(nxdata.attrs['axes'])

#     if verbose:
#         print(f"HDF5 file successfully created: {filepath.name}")
#     return 
# br.Spectra.save_all_as_nexus_map = save_all_as_nexus_map

# %% 



