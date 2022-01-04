#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for file handling."""

import os
import sys
from pathlib import Path
import warnings
import numpy as np

import h5py

import brixs as br
# %%


def image2events(data, cutoff=0):
    """Return a PhotonEvents object from image."""
    X, Y = np.meshgrid(np.arange(data.shape[1]) + 0.5, np.arange(data.shape[0]) + 0.5)
    return np.array([event for event in np.vstack((X.ravel(), Y.ravel(), data.ravel())).transpose() if not event[2]<=cutoff])


def read_ADRESS(filepath):
    """Read files from ADRESS beamline at PSI.

    Example:

        >>> import brixs as br
        >>> pe, s, bad, nd = br.read_ADRESS(filepath)

    Args:
        filepath (str or Path object): filepath to h5 file.

    Returns:
        PhotonEvents, Spectrum, PhotonEvents with bad photon counts, and dictionary
        with instrument parameters.

    Last updated: 26/12/2021 by Carlos Galdino
    """

    f = h5py.File(Path(filepath), 'r')

    nd = {key: val[:] for key, val in f['entry']['instrument']['NDAttributes'].items()}
    x_max = nd['ArraySizeX'][0]
    y_max = nd['ArraySizeY'][0]

    # hist = f['entry']['analysis']['histogram'][:]
    s    = br.Spectrum(f['entry']['analysis']['spectrum'][:])

    pe   = br.PhotonEvents(f['entry']['analysis']['events'][:], x_max=x_max, y_max=y_max)
    # bad  = br.PhotonEvents(br.image2events(f['entry/analysis/bad'][:]), x_max=x_max, y_max=y_max)

    return  pe, s, nd


def get_spectrum_ADRESS(filepath):
    """Read files from ADRESS beamline at PSI and returns only the spectrum.

    Example:

        >>> import brixs as br
        >>> s = br.get_spectrum_ADRESS(filepath)

    Args:
        filepath (str or Path object): filepath to h5 file.

    Returns:
        brixs.Spectrum

    Last updated: 26/12/2021 by Carlos Galdino
    """

    f = h5py.File(Path(filepath), 'r')
    return br.Spectrum(f['entry']['analysis']['spectrum'][:])


def get_pe_ADRESS(filepath):
    """Read files from ADRESS beamline at PSI and return photon events.

    Example:

        >>> import brixs as br
        >>> pe = br.get_pe_ADRESS(filepath)

    Args:
        filepath (str or Path object): filepath to h5 file.

    Returns:
        brixs.PhotonEvents

    Last updated: 01/01/2022 by Carlos Galdino
    """

    f = h5py.File(Path(filepath), 'r')
    x_max = f['entry']['instrument']['NDAttributes']['ArraySizeX'][0]
    y_max = f['entry']['instrument']['NDAttributes']['ArraySizeY'][0]
    return br.PhotonEvents(f['entry']['analysis']['events'][:], x_max=x_max, y_max=y_max)


def get_bad_ADRESS(filepath):
    """Read files from ADRESS beamline at PSI and return bad photon events list.

    Example:

        >>> import brixs as br
        >>> bad = br.get_bad_ADRESS(filepath)

    Args:
        filepath (str or Path object): filepath to h5 file.

    Returns:
        brixs.PhotonEvents

    Last updated: 01/01/2021 by Carlos Galdino
    """

    f = h5py.File(Path(filepath), 'r')
    x_max = f['entry']['instrument']['NDAttributes']['ArraySizeX'][0]
    y_max = f['entry']['instrument']['NDAttributes']['ArraySizeY'][0]
    return br.PhotonEvents(br.image2events(f['entry/analysis/bad'][:]), x_max=x_max, y_max=y_max)


def get_scan_ADRESS(folderpath, prefix, n, zfill=4):
    folderpath = Path(folderpath)
    ss  = br.Spectra()
    pes = [0]*3
    for filepath in [folderpath/(prefix+str(n).zfill(4)+f'_d{i}.h5') for i in (1, 2, 3)]:
        ss.append(get_spectrum_ADRESS(filepath))
        pes.append(get_pe_ADRESS(filepath))
    return ss, pes

def get_scan_spectrum_ADRESS(folderpath, prefix, n, zfill=4):
    folderpath = Path(folderpath)
    ss  = br.Spectra()
    for filepath in [folderpath/(prefix+str(n).zfill(4)+f'_d{i}.h5') for i in (1, 2, 3)]:
        ss.append(get_spectrum_ADRESS(filepath))
    return ss

def get_scan_pe_ADRESS(folderpath, prefix, n, zfill=4):
    folderpath = Path(folderpath)
    pes = [0]*3
    for i, filepath in enumerate([folderpath/(prefix+str(n).zfill(4)+f'_d{i}.h5') for i in (1, 2, 3)]):
        pes[i] = get_pe_ADRESS(filepath)
    return pes


class ReadESRF():
    pass

class ReadIPE():
    pass










#
