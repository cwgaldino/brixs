#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for file handling."""

from pathlib import Path
import numpy as np
import os

import brixs as br

# %% I21 beamline - Diamond - UK ===============================================
import h5py

def h52dict(f):
    data = {}
    if type(f) == h5py._hl.files.File:
        data = h52dict(f['entry'])
    elif type(f) == h5py._hl.group.Group:
        for key in f.keys():
            data[key] = h52dict(f[key])
    elif type(f) == h5py._hl.dataset.Dataset:
        if f.shape == ():
            return None
        else:
            return f[:]
    else:
        print(f'Cannot recognize type {type(f)}. Skipping {f.name}')
    return data

def read(filepath):
    """Read files from I21 beamline at Diamond light source.

    Args:
        filepath (str or Path object): filepath.

    Returns:
        :py:class:`brixs.Image`

    Last updated: 04/04/2022 by Carlos Galdino
    """
    filepath = Path(filepath)

    f = h5py.File(Path(filepath), 'r')
    f = h52dict(f)

    data = np.zeros(f['andor']['data'].shape[1:])
    for temp in f['andor']['data']:
        data += temp
    _ = f['andor'].pop('data')
    _ = f['instrument']['andor'].pop('data')

    im = br.Image(data=data)
    im.nd = f

    return im

def read1(filepath):
    """Read file from I21 beamline at Diamond light source and return a list.

    Args:
        filepath (str or Path object): filepath.

    Returns:
        list with :py:class:`brixs.Image`

    Last updated: 04/04/2022 by Carlos Galdino
    """
    filepath = Path(filepath)

    f = h5py.File(Path(filepath), 'r')
    f = h52dict(f)

    data = []
    for temp in f['andor']['data']:
        temp = np.zeros(f['andor']['data'].shape[1:])
        data.append(br.Image(data=temp))

    _ = f['andor'].pop('data')
    _ = f['instrument']['andor'].pop('data')
    for im in data:
        im.nd = f

    return data








#
