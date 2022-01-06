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

def dispersion_ADRESS(folderpath, prefix, start_energy=None, stop_energy=None, energies=None, start_scan=None, stop_scan=None, scans=None, verbose=False, **kwargs):

    if energies is None:
        energies = np.arange(start_energy, stop_energy+1)
    if scans is None:
        scans = np.arange(start_scan, stop_scan+1)
    if len(scans) != len(energies):
        raise ValueError(f'number of energies ({len(energies)}) do not match the number of scans ({len(scans)})')

    if 'asymmetry' not in kwargs:
        kwargs['asymmetry'] = False

    # spectra data
    spectra_data = {i: {'spectrum': {scan: 0 for scan in scans},
                        'pe': {scan: 0 for scan in scans},
                       'position': {scan: 0 for scan in scans},
                       'y_bin': {scan: 0 for scan in scans},
                       'y_max': {scan: 0 for scan in scans},
                       'fwhm':     {scan: 0 for scan in scans}} for i in range(3)}

    if verbose: print('calculating...')
    for i, scan in enumerate(scans):
        if verbose:  print(f'({i}/{len(scans)-1}) scan = {scan}, energy = {energies[i]}.')
        pes = br.get_scan_pe_ADRESS(folderpath=folderpath, prefix=prefix, n=scan, zfill=4)
        ss = br.get_scan_spectrum_ADRESS(folderpath=folderpath, prefix=prefix, n=scan, zfill=4)
        for j in range(3):
            spectra_data[j]['pe'][scan] = pes[j]
            spectra_data[j]['y_max'][scan] = pes[j].y_max

            spectra_data[j]['spectrum'][scan] = ss[j]
            spectra_data[j]['y_bin'][scan] = len(ss[j].x)
            ss[j].guess_elastic_peak(**kwargs)
            spectra_data[j]['position'][scan] = ss[j].elastic_c
            spectra_data[j]['fwhm'][scan]     = ss[j].elastic_w


    # fit positions
    fit_result = {i: {'func': 0,
                      'slope': 0,
                      'positions': 0,
                      'dispersion': 0,         # energy/bin
                      'dispersion_lenght': 0,  # energy/(detector lenght)
                      'offset': 0} for i in range(3)}

    for i in range(3):
        fit_result[i]['positions'] = [spectra_data[i]['position'][scan] for scan in scans]
        y_max = [spectra_data[i]['y_max'][scan] for scan in scans]
        y_bin = [spectra_data[i]['y_bin'][scan] for scan in scans]

        if y_max.count(y_max[0]) == len(y_max): y_max = y_max[0]
        else: raise ValueError('detector size is not the same for all scans.')
        if y_bin.count(y_bin[0]) == len(y_bin): y_bin = y_bin[0]
        else: raise ValueError('binning is not the same for all scans.')


        popt = np.polyfit(energies, fit_result[i]['positions'], deg=1)
        fit_result[i]['slope'] = popt[0]
        fit_result[i]['dispersion'] = 1/popt[0]
        fit_result[i]['dispersion_lenght'] = 1/popt[0] /(y_max / y_bin)
        fit_result[i]['offset'] = popt[1]
        fit_result[i]['func'] = np.poly1d(popt)



    final_dispersion        = [fit_result[i]['dispersion'] for i in range(3)]
    final_dispersion_lenght = [fit_result[i]['dispersion_lenght'] for i in range(3)]
    final_resolution        = [np.mean([spectra_data[i]['fwhm'][scan]*final_dispersion[i] for scan in scans]) for i in range(3)]


    if verbose:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for i in range(3):
            fwhm = round(final_resolution[i]*1000, 1)
            disp = round(fit_result[i]['dispersion']*1000, 3)
            disp_lenght = round(fit_result[i]['dispersion_lenght']*1000, 3)

            ax.plot(energies, fit_result[i]['positions'], marker='o', lw=0, label=f'ccd {i}, {disp_lenght} meV/det.len. = {disp} mev/bin ({fwhm} meV)')
            x = np.linspace(energies[0], energies[-1], len(energies*10))
            y = fit_result[i]['func'](x)
            ax.plot(x, y, color='black')
        plt.legend()
        plt.ylabel('elastic position (bin)')
        plt.xlabel('energy (eV)')
        plt.title(f'{round(np.mean(final_dispersion)*1000, 3)} meV/bin = {round(np.mean(final_dispersion_lenght)*1000, 3)} meV/det.len. ({round(np.mean(final_resolution)*1000, 1)} meV)')

    return np.array(final_dispersion_lenght), np.array(final_dispersion), np.array(final_resolution), fit_result, spectra_data

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
