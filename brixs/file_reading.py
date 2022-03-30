#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for file handling."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

import brixs as br

try:
    import h5py
except:
    pass

# %% General ===================================================================

def image2events(data, cutoff=0):
    """Return a PhotonEvents object from image."""
    X, Y = np.meshgrid(np.arange(data.shape[1]) + 0.5, np.arange(data.shape[0]) + 0.5)
    return np.array([event for event in np.vstack((X.ravel(), Y.ravel(), data.ravel())).transpose() if not event[2]<=cutoff])

# %% ADRESS beamline - PSI - Switzerland =======================================
def read_ADRESS(*args, **kwargs):
    """Read files from ADRESS beamline at PSI.

    Example:

        >>> import brixs as br
        >>>
        >>> # single file
        >>> s = br.read_ADRESS(filepath)
        >>>
        >>> # data for each ccd
        >>> ss = br.read_ADRESS(folderpath, prefix, n)
        >>>
        >>> # One can also use keyword arguments
        >>> s  = br.read_ADRESS(filepath=filepath)
        >>> ss = br.read_ADRESS(folderpath=folderpath, prefix=prefix, n=n, zfill=zfill)
        >>>
        >>> # experiment parameters can be accesed by:
        >>> print(s.nd)

    Args:
        filepath (str or Path object): filepath to h5 file (for loading a single file).
        folderpath (str or Path object): folderpath to files (for loading a data from 3 cdd files).
        prefix (str): prefix (with or without 'underscore'). Example: for
            'Cu_0001_d1.h5' the prefix is 'Cu' or 'Cu_'.
        n (number): scan number. Eample: for 'Cu_0001_d1.h5' the file number is 1.
            The number of extra zeros (0001) is defined by the argument `zfill`.
        zfill (number, optional): number of digits to fill scan number (Default is 4).

    Returns:
        Spectrum

    Last updated: 21/03/2022 by Carlos Galdino
    """
    filepath, folderpath, prefix, n, zfill = _sort_args_ADRESS(*args, **kwargs)

    if filepath is not None:
        return _read_ADRESS1(filepath)
    if folderpath is not None:
        if folderpath.is_dir() == False:
            raise AttributeError(f'It seems like the folderpath is not a folder\nfolderpath = {folderpath}')
        if prefix is None:
            raise AttributeError('Missing prefix.\nPlease, inform the file prefix.\nExample: for Cu_0001_d1.h5 the prefix is `Cu` or `Cu_`.')
        if n is None:
            raise AttributeError('Missing file number (scan number).\nPlease, inform the file number.\nExample: for Cu_0001_d1.h5 the file number is 1.')
        filepaths = [folderpath/(prefix+str(n).zfill(zfill)+f'_d{i}.h5') for i in (1, 2, 3)]
        ss  = br.Spectra(3)
        nd = [0]*3
        for i, filepath in enumerate(filepaths):
            ss[i] = _read_ADRESS1(filepath)
            ss[i].scan = n
            ss[i].ccd = i
        ss.scan = n
        return ss

def read_pe_ADRESS(*args, **kwargs):
    """Return photon event list from ADRESS beamline at PSI.

    Example:

        >>> import brixs as br
        >>>
        >>> # single file
        >>> pe = br.read_pe_ADRESS(filepath)
        >>>
        >>> # data for each ccd
        >>> pes = br.read_pe_ADRESS(folderpath, prefix, n)
        >>>
        >>> # One can also use keyword arguments
        >>> pe  = br.read_pe_ADRESS(filepath=filepath)
        >>> pes = br.read_pe_ADRESS(folderpath=folderpath, prefix=prefix, n=n, zfill=zfill)
        >>>
        >>> # experiment parameters can be accesed by:
        >>> print(pe.nd)

    Args:
        filepath (str or Path object): filepath to h5 file (for loading a single file).
        folderpath (str or Path object): folderpath to files (for loading a data from 3 cdd files).
        prefix (str): prefix (with or without 'underscore'). Example: for
            'Cu_0001_d1.h5' the prefix is 'Cu' or 'Cu_'.
        n (number): scan number. Eample: for 'Cu_0001_d1.h5' the file number is 1.
            The number of extra zeros (0001) is defined by the argument `zfill`.
        zfill (number, optional): number of digits to fill scan number (Default is 4).

    Returns:
        PhotonEvents object

    Last updated: 21/02/2022 by Carlos Galdino
    """
    filepath, folderpath, prefix, n, zfill = _sort_args_ADRESS(*args, **kwargs)

    if filepath is not None:
        return _read_pe_ADRESS1(filepath)
    if folderpath is not None:
        if folderpath.is_dir() == False:
            raise AttributeError(f'It seems like the folderpath is not a folder\nfolderpath = {folderpath}')
        if prefix is None:
            raise AttributeError('Missing prefix.\nPlease, inform the file prefix.\nExample: for Cu_0001_d1.h5 the prefix is `Cu` or `Cu_`.')
        if n is None:
            raise AttributeError('Missing file number (scan number).\nPlease, inform the file number.\nExample: for Cu_0001_d1.h5 the file number is 1.')
        filepaths = [folderpath/(prefix+'_'+str(n).zfill(zfill)+f'_d{i}.h5') for i in (1, 2, 3)]
        pes = [0]*3
        nd = [0]*3
        for i, filepath in enumerate(filepaths):
            pes[i] = _read_pe_ADRESS1(filepath)
            pes[i].scan = n
            pes[i].ccd = i
        pes.scan = n
        return pes

def _sort_args_ADRESS(*args, **kwargs):
    filepath   = None
    folderpath = None
    prefix     = None
    n          = None
    zfill      = 4
    print(args)
    print(kwargs)
    if len(args) > 0 and len(kwargs) > 0:
        raise AttributeError(f'cannot mix positional arguments with keyword arguents.\nKeyword args: {kwargs}\nPosition args: {args}')

    # kwargs
    for key in kwargs.keys():
        if key not in ['filepath', 'folderpath', 'prefix', 'n', 'zfill']:
            raise AttributeError(f'Cannot identify argument {key}')
    if 'filepath' in kwargs:
        filepath = Path(kwargs['filepath'])
    if 'folderpath' in kwargs:
        folderpath = Path(kwargs['folderpath'])
    if 'prefix' in kwargs:
        prefix = kwargs['prefix']
    if 'n' in kwargs:
        n = kwargs['n']
    if 'zfill' in kwargs:
        zfill = kwargs['zfill']

    # args
    if len(args) > 0:
        if len(args) == 1:
            filepath = Path(args[0])
        elif len(args) == 3:
            folderpath = Path(args[0])
            prefix     = args[1]
            n          = args[2]
        elif len(args) == 4:
            folderpath = Path(args[0])
            prefix     = args[1]
            n          = args[2]
            zfill      = args[3]
        else:
            raise AttributeError(f'Cannot read positional attributes.\nAttributes = {args}')

    if prefix is not None:
        if prefix[-1] != '_':
            prefix += '_'

    return filepath, folderpath, prefix, n, zfill

def _read_ADRESS1(filepath):
    f = h5py.File(filepath, 'r')

    nd = {key: val[:] for key, val in f['entry']['instrument']['NDAttributes'].items()}
    s    = br.Spectrum(f['entry']['analysis']['spectrum'][:])
    s.nd = nd
    return s

def _read_pe_ADRESS1(filepath):
    f = h5py.File(Path(filepath), 'r')

    nd = {key: val[:] for key, val in f['entry']['instrument']['NDAttributes'].items()}
    x_max = nd['ArraySizeX'][0]
    y_max = nd['ArraySizeY'][0]

    pe   = br.PhotonEvents(f['entry']['analysis']['events'][:], x_max=x_max, y_max=y_max)
    pe.nd = nd

    return pe

def read_bad_ADRESS1(filepath):
    """Return image with bad events."""
    f = h5py.File(Path(filepath), 'r')
    nd = {key: val[:] for key, val in f['entry']['instrument']['NDAttributes'].items()}
    x_max = nd['ArraySizeX'][0]
    y_max = nd['ArraySizeY'][0]

    bad = np.array(f['entry/analysis/bad'][:])

    return  bad

def calculate_calib_ADRESS(folderpath, prefix, start_scan=None, stop_scan=None, scans=None, start_value=None, stop_value=None, values=None, verbose=False, **kwargs):
    """Returns the calibration factor (dispersion) for ADRESS data.

    If energy values are not given, it will be read from the file.

    >>> import brixs as br
    >>> disp, sss = br.calculate_calib_ADRESS(folderpath, prefix, start_scan, stop_scan)
    >>>
    >>> # disp is a list with the calculated dispersion for each ccd
    >>> disp_ccd0 = disp[0]
    >>> disp_ccd1 = disp[1]
    >>> disp_ccd2 = disp[2]
    >>> disp = np.mean(disp)
    >>>
    >>> # sss is a list of spectra for each ccd
    >>> ss_ccd0 = sss[0]
    >>> ss_ccd1 = sss[1]
    >>> ss_ccd2 = sss[2]

    Returns:
        dispersion and list of brixs.Spectra object.
    """
    pick_values_from_files = False
    if scans is None:
        scans = np.arange(start_scan, stop_scan+1)
    if values is None:
        if start_value is None and stop_value is None:
            pick_values_from_files = True
            values = np.zeros(len(scans))
        else:
            values = np.linspace(start_value, stop_value, len(scans))
    else:
        if len(scans) != len(values):
            raise ValueError(f'number of values ({len(values)}) do not match the number of scans ({len(scans)})')

    disp     = [0, 0, 0]
    disp_bin = [0, 0, 0]
    sss      = [0, 0, 0]
    for ccd in range(3):
        ss = br.Spectra(n=len(scans))
        y_max = None
        y_bin = None
        values_temp = np.zeros(len(scans))
        for i, scan in enumerate(scans):
            ss_temp = br.read_ADRESS(folderpath=folderpath, prefix=prefix, n=scan, zfill=4)
            ss_temp[ccd].photon_energy = np.mean(ss_temp[ccd].nd['PhotonEnergy'])
            values_temp[i] = ss_temp[ccd].photon_energy
            ss[i] = ss_temp[ccd]

            if y_max is None:
                y_max = ss_temp[ccd].nd['ArraySizeY'][0]
                y_bin = len(ss_temp[ccd].x)

            if y_max != ss_temp[ccd].nd['ArraySizeY'][0]:
                raise ValueError('detector size is not the same for all scans.')
            if y_bin != len(ss_temp[ccd].x):
                raise ValueError('binning is not the same for all scans.')

            if verbose:
                print(f'ccd: {ccd}, scan: {scans[i]}, value: {values[i]}')
        if pick_values_from_files:
            values = values_temp

        disp_bin[ccd] = ss.calculate_calib(values=values, **kwargs)
        disp[ccd] = disp_bin[ccd] / (y_max / y_bin)
        ss.ccd = ccd
        sss[ccd] = ss

    # return disp, disp_bin, sss
    return disp_bin, sss



# def read_ADRESS1_full(filepath):
#     """Read files from ADRESS beamline at PSI.
#
#     Example:
#
#         >>> import brixs as br
#         >>> s, pe, nd = br.read_ADRESS_full(filepath)
#
#     Args:
#         filepath (str or Path object): filepath to h5 file.
#
#     Returns:
#         Spectrum, PhotonEvents, and dictionary with instrument parameters.
#
#     Last updated: 21/02/2022 by Carlos Galdino
#     """
#     f = h5py.File(Path(filepath), 'r')
#
#     nd = {key: val[:] for key, val in f['entry']['instrument']['NDAttributes'].items()}
#     x_max = nd['ArraySizeX'][0]
#     y_max = nd['ArraySizeY'][0]
#
#     # hist = f['entry']['analysis']['histogram'][:]
#     s    = br.Spectrum(f['entry']['analysis']['spectrum'][:])
#
#     pe   = br.PhotonEvents(f['entry']['analysis']['events'][:], x_max=x_max, y_max=y_max)
#     # bad  = br.PhotonEvents(br.image2events(f['entry/analysis/bad'][:]), x_max=x_max, y_max=y_max)
#
#     return  s, pe, nd

# def read_ADRESS_bad(filepath):
#     """Read files from ADRESS beamline at PSI and return bad photon events list.
#
#     Example:
#
#         >>> import brixs as br
#         >>> bad = br.get_bad_ADRESS(filepath)
#
#     Args:
#         filepath (str or Path object): filepath to h5 file.
#
#     Returns:
#         brixs.PhotonEvents
#
#     Last updated: 01/01/2021 by Carlos Galdino
#     """
#
#     f = h5py.File(Path(filepath), 'r')
#     x_max = f['entry']['instrument']['NDAttributes']['ArraySizeX'][0]
#     y_max = f['entry']['instrument']['NDAttributes']['ArraySizeY'][0]
#     return br.PhotonEvents(br.image2events(f['entry/analysis/bad'][:]), x_max=x_max, y_max=y_max)

# def read_ADRESS1(filepath):
#     """Read files from ADRESS beamline at PSI and returns only the spectrum.
#
#     Example:
#
#         >>> import brixs as br
#         >>> s = br.get_spectrum_ADRESS(filepath)
#
#     Args:
#         filepath (str or Path object): filepath to h5 file.
#
#     Returns:
#         brixs.Spectrum
#
#     Last updated: 21/02/2022 by Carlos Galdino
#     """
#     f = h5py.File(Path(filepath), 'r')
#     return br.Spectrum(f['entry']['analysis']['spectrum'][:])
#
# def read_ADRESS_full(folderpath, prefix, n, zfill=4):
#     """Read files from the three ccd's from ADRESS beamline at PSI.
#
#     Example:
#
#         >>> import brixs as br
#         >>> ss, pes, nds = br.read_ADRESS_full(folderpath, prefix, n)
#
#     Args:
#         folderpath (str or Path object): folderpath to files.
#         prefix (str): prefix (without 'underscore')
#         n (number): scan number
#         zfill (number, optional): number of digits to fill scan number.
#
#     Returns:
#         Spectra, list of PhotonEvents, and dictionary with instrument parameters.
#
#     Last updated: 21/02/2022 by Carlos Galdino
#     """
#     folderpath = Path(folderpath)
#     filepaths = [folderpath/(prefix+'_'+str(n).zfill(4)+f'_d{i}.h5') for i in (1, 2, 3)]
#     ss  = br.Spectra(3)
#     pes = [0]*3
#     nd = [0]*3
#     for i, filepath in enumerate(filepaths):
#         ss[i], pes[i], nd[i] = read_ADRESS1_full(filepath)
#         ss[i].scan = n
#         ss[i].ccd = i
#         ss[i].nd = nd[i]
#     ss.scan = n
#     return ss, pes, nd
#
# def read_ADRESS(folderpath, prefix, n, zfill=4):
    """Read files from the three ccd's from ADRESS beamline at PSI.

    Example:

        >>> import brixs as br
        >>> ss = br.read_ADRESS(folderpath, prefix, n)

    Args:
        folderpath (str or Path object): folderpath to files.
        prefix (str): prefix (without 'underscore')
        n (number): scan number
        zfill (number, optional): number of digits to fill scan number.

    Returns:
        Spectra object

    Last updated: 21/02/2022 by Carlos Galdino
    """
    folderpath = Path(folderpath)
    ss  = br.Spectra()
    for filepath in [folderpath/(prefix+'_'+str(n).zfill(4)+f'_d{i}.h5') for i in (1, 2, 3)]:
        ss.append(read_ADRESS1(filepath))
    return ss


# %% PEAXIS beamline - HZB - Germany ===========================================
import os
def ReadAndor(fname, dimensions = (2048, 2048), byte_size=4, data_type='c'):
    """Reads the .sif file produced by the Andor iKon-L CCD camera.
    The dimensions of the pixel array can be changed if necessary.
    """
    size = os.path.getsize(fname)
    source = open(fname, 'rb')
    nrows = dimensions[-1]
    header, data = [],[]
    bindata = None
    total_length = 0
    header = np.fromfile(source, np.byte, size - 4*dimensions[0]*dimensions[1] -2*4, '') # 2772 is close
    # print("Headersize",  size - 4*dimensions[0]*dimensions[1] -2*4)
    data = np.fromfile(source, np.float32, dimensions[0]*dimensions[1], '')
    header = bytes(header)
    lines = header.split(b'\n')
    header = []
    for n, line in enumerate(lines):
        try:
            header.append(line.decode('ascii'))
        except UnicodeDecodeError:
            print('header lines skipped:', n+1, 'with length:', len(line))
    source.close()
    return np.array(data).reshape(dimensions), header

def read_PEAXIS(filepath):
    """Read files from PEAXIS beamline at HZB.

    Example:

        >>> import brixs as br
        >>> pe, s, bad, nd = br.read_PEAXIS(filepath)

    Args:
        filepath (str or Path object): filepath to h5 file.

    Returns:
        PhotonEvents, Spectrum, PhotonEvents with bad photon counts, and dictionary
        with instrument parameters.

    Last updated: 26/12/2021 by Carlos Galdino
    """
    data, header = ReadAndor(filepath)
    return br.Spectrum(np.sum(data, axis=1))

def read_PEAXIS2(filepath, cutoff=0):
    """Read files from PEAXIS beamline at HZB.

    Example:

        >>> import brixs as br
        >>> pe, s, bad, nd = br.read_PEAXIS(filepath)

    Args:
        filepath (str or Path object): filepath to h5 file.

    Returns:
        PhotonEvents, Spectrum, PhotonEvents with bad photon counts, and dictionary
        with instrument parameters.

    Last updated: 26/12/2021 by Carlos Galdino
    """
    data, header = ReadAndor(filepath)
    pe = br.PhotonEvents(br.image2events(data[::-1,:], cutoff=cutoff), y_max=2048, x_max=2048)
    return pe, header

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
