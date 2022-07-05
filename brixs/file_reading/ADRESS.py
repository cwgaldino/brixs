#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for file handling."""

from pathlib import Path
import numpy as np
import os
import h5py

import brixs as br

# %% ADRESS beamline - PSI - Switzerland =======================================
def read(*args, **kwargs):
    """Read files from ADRESS beamline at PSI.

    If a specif filepath is given, this file is imported. If a folderpath, prefix,
    and scan number (n) is given, data from the 3 ccds are imported.

    Example:

        >>> import brixs as br
        >>>
        >>> # single file
        >>> s = br.ADRESS.read(filepath)
        >>>
        >>> # data for each ccd
        >>> ss = br.ADRESS.read(folderpath, prefix, n)
        >>>
        >>> # One can also use keyword arguments
        >>> s  = br.ADRESS.read(filepath=filepath)
        >>> ss = br.ADRESS.read(folderpath=folderpath, prefix=prefix, n=n, zfill=zfill)
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
        :py:class:`brixs.Spectrum` or :py:class:`brixs.Spectra`

    Last updated: 04/04/2022 by Carlos Galdino
    """
    filepath, folderpath, prefix, n, zfill = _sort_args(*args, **kwargs)

    if filepath is not None:
        return _read_1(filepath)
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
            ss[i] = _read_1(filepath)
            ss[i].scan = n
            ss[i].ccd = i
        ss.scan = n
        return ss

def read_pe(*args, **kwargs):
    """Return photon event list from ADRESS beamline at PSI.

    If a specif filepath is given, this file is imported. If a folderpath, prefix,
    and scan number (n) is given, data from the 3 ccds are imported.

    Example:

        >>> import brixs as br
        >>>
        >>> # single file
        >>> pe = br.ADRESS.read_pe(filepath)
        >>>
        >>> # data for each ccd
        >>> pes = br.ADRESS.read_pe(folderpath, prefix, n)
        >>>
        >>> # One can also use keyword arguments
        >>> pe  = br.ADRESS.read_pe(filepath=filepath)
        >>> pes = br.ADRESS.read_pe(folderpath=folderpath, prefix=prefix, n=n, zfill=zfill)
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
        :py:class:`brixs.PhotonEvents`

    Last updated: 04/04/2022 by Carlos Galdino
    """
    filepath, folderpath, prefix, n, zfill = _sort_args(*args, **kwargs)

    if filepath is not None:
        return _read_1_pe(filepath)
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
            pes[i] = _read_1_pe(filepath)
            pes[i].scan = n
            pes[i].ccd = i
        pes.scan = n
        return pes

def read_bad(filepath):
    """Return image with bad events read from ADRESS beamline at PSI.

    Example:

        >>> import brixs as br
        >>> im = br.ADRESS.read_bad(filepath)
        >>> print(im.nd)

    Args:
        filepath (str or Path object): filepath.

    Returns:
        :py:class:`brixs.Image`

    Last updated: 04/04/2022 by Carlos Galdino
    """
    f = h5py.File(Path(filepath), 'r')
    nd = {key: val[:] for key, val in f['entry']['instrument']['NDAttributes'].items()}

    im = br.Image(np.array(f['entry/analysis/bad'][:]))
    im.nd = nd

    return  im

def calculate_calib(folderpath, prefix, mode='cc', start=None, stop=None, scans=None, start_value=None, stop_value=None, values=None, verbose=False):
    """Returns the calibration factor (dispersion) for ADRESS data.

    Args:
        folderpath:
        prefix:
        mode:
        ...

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
    # scans
    if scans is None:
        scans = np.arange(start, stop+1)

    # values
    pick_values_from_files = False
    if values is None:
        if start_value is None and stop_value is None:
            pick_values_from_files = True
            values = np.zeros(len(scans))
        else:
            values = np.linspace(start_value, stop_value, len(scans))
    assert len(scans) == len(values), f'number of values ({len(values)}) do not match the number of scans ({len(scans)})'

    # CALCULATION
    disp     = [0, 0, 0]
    disp_bin = [0, 0, 0]
    sss      = [0, 0, 0]
    for ccd in range(3):
        ss = br.Spectra(n=len(scans))
        y_max = None
        y_bin = None
        for i, scan in enumerate(scans):
            ss_temp = read(folderpath=folderpath, prefix=prefix, n=scan, zfill=4)
            ss_temp[ccd].photon_energy = np.mean(ss_temp[ccd].nd['PhotonEnergy'])
            ss[i] = ss_temp[ccd]
            if pick_values_from_files:
                values[i] = ss_temp[ccd].photon_energy

            if y_max is None:
                y_max = ss_temp[ccd].nd['ArraySizeY'][0]
                y_bin = len(ss_temp[ccd].x)

            if y_max != ss_temp[ccd].nd['ArraySizeY'][0]:
                raise ValueError('detector size is not the same for all scans.')
            if y_bin != len(ss_temp[ccd].x):
                raise ValueError('binning is not the same for all scans.')

            if verbose:
                print(f'ccd: {ccd}, scan: {scans[i]}, value: {values[i]}')

        if mode == 'fitted peaks':
            ss.find_peaks(prominence=50)
            p = [0.0]*len(ss)
            for i in range(len(ss)):
                p[i] = len(ss[i].peaks)
            assert br.backpack.all_equal(p) and p[0]==1, f'Some spectra have more the one peak.\nSpectra should have just the elastic line.\nccd number: {ccd}\nNumber of peaks for each spectra: {p}'
            ss.fit_peaks()
        if mode == 'peak':
            ss.find_peaks(prominence=50)
            p = [0.0]*len(ss)
            for i in range(len(ss)):
                p[i] = len(ss[i].peaks)
            assert br.backpack.all_equal(p) and p[0]==1, f'Some spectra have more the one peak.\nSpectra should have just the elastic line.\nccd number: {ccd}\nNumber of peaks for each spectra: {p}'
        disp_bin[ccd] = ss.calculate_calib(values=values, mode=mode, peak=0)
        disp[ccd] = disp_bin[ccd] / (y_max / y_bin)
        ss.ccd = ccd
        sss[ccd] = ss

    return disp, disp_bin, sss
    # return disp_bin, sss

def read_xas(folderpath, prefix, n, zfill=4):
    """Read XAS files from ADRESS beamline at PSI.

    Args:
        folderpath (str or Path object): folderpath to files.
        prefix (str): prefix (with or without 'underscore'). Example: for
            'Cu_0001_d1.h5' the prefix is 'Cu' or 'Cu_'.
        n (number): scan number. Eample: for 'Cu_0001_d1.h5' the file number is 1.
            The number of extra zeros (0001) is defined by the argument `zfill`.
        zfill (number, optional): number of digits to fill scan number (Default is 4).

    Returns:
        TEY (:py:class:`brixs.Spectrum`), TFY (:py:class:`brixs.Spectrum`)

    Last updated: 20/05/2022 by Carlos Galdino
    """

    if folderpath.is_dir() == False:
        raise AttributeError(f'It seems like the folderpath is not a folder\nfolderpath = {folderpath}')

    if prefix is not None:
        if prefix[-1] != '_':
            prefix += '_'

    filepath = folderpath/(prefix+str(n).zfill(zfill)+'.xas')

    # data
    data = br.backpack.load_data(filepath)
    TEY = br.Spectrum(x=data['E_'+str(n).zfill(zfill)], y=data['D1_'+str(n).zfill(zfill)])
    TFY = br.Spectrum(x=data['E_'+str(n).zfill(zfill)], y=data['D2_'+str(n).zfill(zfill)])
    RMU = br.Spectrum(x=data['E_'+str(n).zfill(zfill)], y=data['D3_'+str(n).zfill(zfill)])

    # header
    header = br.backpack.load_Comments(filepath)[1:-7]
    nd = {}
    for line in header:
        split = line[1:].split(':')
        name = split[0].strip()
        value = split[1].strip()
        nd[name] = value
    TEY.nd = nd
    TFY.nd = nd
    RMU.nd = nd

    return TEY, TFY, RMU


# %% SUPPORT ===================================================================

def _sort_args(*args, **kwargs):
    filepath   = None
    folderpath = None
    prefix     = None
    n          = None
    zfill      = 4
    # print(args)
    # print(kwargs)
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

def _read_1_pe(filepath):
    f = h5py.File(Path(filepath), 'r')

    nd = {key: val[:] for key, val in f['entry']['instrument']['NDAttributes'].items()}
    x_max = nd['ArraySizeX'][0]
    y_max = nd['ArraySizeY'][0]

    pe   = br.PhotonEvents(data=f['entry']['analysis']['events'][:], shape=(y_max, x_max))
    pe.nd = nd

    return pe

def _read_1(filepath):
    f = h5py.File(filepath, 'r')

    nd = {key: val[:] for key, val in f['entry']['instrument']['NDAttributes'].items()}
    s    = br.Spectrum(f['entry']['analysis']['spectrum'][:])
    s.nd = nd
    return s
