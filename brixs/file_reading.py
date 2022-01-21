#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for file handling."""

from pathlib import Path
import numpy as np
import h5py
import matplotlib.pyplot as plt

import brixs as br

# %% General ===================================================================

def image2events(data, cutoff=0):
    """Return a PhotonEvents object from image."""
    X, Y = np.meshgrid(np.arange(data.shape[1]) + 0.5, np.arange(data.shape[0]) + 0.5)
    return np.array([event for event in np.vstack((X.ravel(), Y.ravel(), data.ravel())).transpose() if not event[2]<=cutoff])

# %% ADRESS beamline - PSI - Switzerland =======================================
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
    filepaths = [folderpath/(prefix+str(n).zfill(4)+f'_d{i}.h5') for i in (1, 2, 3)]
    ss  = br.Spectra(3)
    pes = [0]*3
    for i, filepath in enumerate(filepaths):
        ss[i] = get_spectrum_ADRESS(filepath)
        pes[i] = get_pe_ADRESS(filepath)
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

def get_dispersion_ADRESS(folderpath, prefix, start_energy=None, stop_energy=None, energies=None, start_scan=None, stop_scan=None, scans=None, **kwargs):

    if energies is None:
        energies = np.arange(start_energy, stop_energy+1)
    if scans is None:
        scans = np.arange(start_scan, stop_scan+1)
    if len(scans) != len(energies):
        raise ValueError(f'number of energies ({len(energies)}) do not match the number of scans ({len(scans)})')

    disp = [0, 0, 0]
    disp_bin = [0, 0, 0]
    # fwhm = [0, 0, 0]
    sss = [0, 0, 0]
    for ccd in range(3):
        ss = br.Spectra(n=len(scans))
        y_max = None
        y_bin = None
        for i, scan in enumerate(scans):
            ss_temp, pes_temp = br.get_scan_ADRESS(folderpath=folderpath, prefix=prefix, n=scan, zfill=4)
            ss[i] = ss_temp[ccd]
            if y_max is None:
                y_max = pes_temp[ccd].y_max
                y_bin = len(ss_temp[ccd].x)
            else:
                if y_max != pes_temp[ccd].y_max:
                    raise ValueError('detector size is not the same for all scans.')
                if y_bin != len(ss_temp[ccd].x):
                    raise ValueError('binning is not the same for all scans.')

        disp_bin[ccd] = ss.calculate_dispersion(energies=energies, **kwargs)
        disp[ccd] = disp_bin[ccd] / (y_max / y_bin)
        sss[ccd] = ss

    return disp, disp_bin, sss


# %% ID32 beamline - ESRF - France =============================================
def read_ID32(filepath):
    """Read files from ID32 beamline at ESRF.

    Example:

        >>> import brixs as br
        >>> pe, s, bad, nd = br.read_ESRF(filepath)

    Args:
        filepath (str or Path object): filepath to h5 file.

    Returns:
        PhotonEvents, Spectrum, PhotonEvents with bad photon counts, and dictionary
        with instrument parameters.

    Last updated: 26/12/2021 by Carlos Galdino
    """
    raise NotImplementedError('not implemented yet.')

# %% IPE beamline - SIRIUS - Brazil ============================================
def read_IPE(filepath):
    """Read files from IPE beamline at SIRIUS.

    Example:

        >>> import brixs as br
        >>> pe, s, bad, nd = br.read_IPE(filepath)

    Args:
        filepath (str or Path object): filepath to h5 file.

    Returns:
        PhotonEvents, Spectrum, PhotonEvents with bad photon counts, and dictionary
        with instrument parameters.

    Last updated: 26/12/2021 by Carlos Galdino
    """
    raise NotImplementedError('not implemented yet.')











#
