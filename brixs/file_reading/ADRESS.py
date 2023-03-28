#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from ADRESS beamline - PSI.

Last edited: Carlos Galdino 03-2023
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import numpy as np
from collections.abc import Iterable

# %% ------------------------- Special Imports ---------------------------- %% #
import brixs as br
try:
    import h5py
except ModuleNotFoundError:
    pass

# %% ----------------------------- Support -------------------------------- %% #
def _sort_args(*args, **kwargs):
    """Sort arguments given to a function.

    Allows for more flexibility when using functions from this and other files.

    Args:
        filepath (str or Path object, optional): filepath. 
        folderpath (str or Path object, optional): folderpath.
        prefix (str, optional): file prefix string.
        n (int, optional): scan number (number following of prefix).
        zfill (int, optional): fill "n" zeros on the left side. Default is 4.

    For positional arguments, the number of args define what they are: 
        1 arg: filepath
        2 args: Raises AttributeError
        3 args: folderpath, prefix, n
        4 args: folderpath, prefix, n, zfill

    Raises:
        AttributeError: positional arguments are mixed with keyword arguments.
        AttributeError: if 2 positional arguments are given.
        AttributeError: if 5 or more positional arguments are given.

    Returns: 
        filepath, folderpath, prefix, n, zfill

    Last updated: Carlos Galdino 03-2023
    """
    filepath   = None
    folderpath = None
    prefix     = None
    n          = None
    zfill      = 4

    # check if positional arguments are mixed with keyword arguments
    if len(args) > 0 and len(kwargs) > 0:
        raise AttributeError(f'cannot mix positional arguments with keyword arguments.\nKeyword args: {kwargs}\nPosition args: {args}')

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

    return filepath, folderpath, prefix, n, zfill

def _read_1(filepath, type_='spectrum'):
    """Read one file from ADRESS beamline.

    Args:
        filepath (str or Path object): filepath.
        type_ (str, optional): data to return (case insensitive). Use 
            'spectrum' or 's' for spectrum; 'photon events' or 'pe', for photon
            events list; and 'bad' or 'b' for bad events image.  
    
    Returns:
        brixs.Spectrum() or brixs.PhotonEvents()
    
    Last updated: Carlos Galdino 03-2023
    """
    # read file
    try:
        f = h5py.File(Path(filepath), 'r')
    except NameError:
        raise ModuleNotFoundError('Cannot find h5py. Make sure h5py is installed.') 

    # get metadata
    nd = {key: val[:] for key, val in f['entry']['instrument']['NDAttributes'].items()}

    if type_.lower() in ['spectrum', 's']:
        # spectrum
        s    = br.Spectrum(f['entry']['analysis']['spectrum'][:])

        # save attr to object
        for attr in nd:
            setattr(s, attr, nd[attr])
        # s.nd = nd
        return s
    elif type_.lower() in ['photon events', 'pe']:
        # get array size
        x_max = nd['ArraySizeX'][0]
        y_max = nd['ArraySizeY'][0]
        # get photon events list
        pe   = br.PhotonEvents(data=f['entry']['analysis']['events'][:], shape=(y_max, x_max))
        
        # save attr to object
        for attr in nd:
            setattr(pe, attr, nd[attr])
        # pe.nd = nd
        return pe
    elif type_.lower() in ['bad', 'b']:
        # bad events image
        im = br.Image(np.array(f['entry/analysis/bad'][:]))
        # save attr to object
        for attr in nd:
            setattr(im, attr, nd[attr])
        im.nd = nd
        return  im
    else:
        raise ValueError("type_ must be one of 'spectrum', 'photon events', or 'bad'.")

def _read_xas(filepath):
    """Read xas file from ADRESS beamline.

    Args:
        filepath (str or Path object): filepath.
    
    Returns:
        TEY, TFY, RMU (:py:class:`brixs.Spectrum`)
    
    Last updated: Carlos Galdino 03-2023
    """
    # data
    data = br.load_data(filepath, force_array=True)
    TEY = br.Spectrum(x=data[:, 0], y=data[:, 1])
    TFY = br.Spectrum(x=data[:, 0], y=data[:, 2])
    RMU = br.Spectrum(x=data[:, 0], y=data[:, 3])

    # metadata
    header = br.load_Comments(filepath)[1:-7]
    nd = {}
    for line in header:
        split = line[1:].split(':')
        name  = split[0].strip()
        value = split[1].strip()
        nd[name] = value
    TEY.nd = nd
    TFY.nd = nd
    RMU.nd = nd

    return TEY, TFY, RMU

# %% ---------------------------------------------------------------------- %% #
def read(*args, **kwargs):
    """Read files from ADRESS beamline at PSI.

    Args:
        filepath (str or Path object): filepath to h5 file (for loading a single file).
        folderpath (str or Path object): folderpath to files (for loading a data from 3 cdd files).
        prefix (str): file prefix.
        n (number): scan number. Example: for 'Cu_0001_d1.h5' the file number is 1.
            The number of extra zeros (0001) is defined by the argument `zfill`.
        zfill (number, optional): number of digits to fill scan number (Default is 4).
        type_ (str, optional): data to return (not case sensitive). Use 
            'spectrum' or 's' for spectrum; 'photon events' or 'pe', for photon
            events list; 'bad' or 'b' for bad events image; and 'xas' or 'x' for
            xas data.

    Example:

        >>> s   = br.ADRESS.read(filepath)
        >>> pe  = br.ADRESS.read(filepath, type_='pe')
        >>> bad = br.ADRESS.read(filepath, type_='bad')

    Experiment parameters can be accessed by:
   
        >>> print(s.nd)

    ADRESS beamline has three ccds. Data from all three ccd's can be extracted
    simultaneously:

        >>> ss   = br.ADRESS.read(folderpath, prefix, scan)
        >>> pes  = br.ADRESS.read(folderpath, prefix, scan, type_='pe')
        >>> bads = br.ADRESS.read(folderpath, prefix, scan, type_='bad')

    X-ray absorption data (XAS) from ADRESS data can be read by:

        >>> TEY, TFY, RMU = br.ADRESS.read(filepath, type_='xas')
        >>> TEY, TFY, RMU = br.ADRESS.read(folderpath, prefix, scan, type_='xas')

    Returns:
        :py:class:`brixs.Spectrum` or :py:class:`brixs.Spectra`

        if ``type_ = 'pe'``: :py:class:`brixs.PhotonEvents` or a list of :py:class:`brixs.PhotonEvents`.

        if ``type_ = 'bad'``: :py:class:`brixs.Image` or a list of :py:class:`brixs.Image`

        if ``type_ = 'xas'``: TEY (:py:class:`brixs.Spectrum`), TFY (:py:class:`brixs.Spectrum`), TEY (:py:class:`brixs.Spectrum`), RMU (:py:class:`brixs.Spectrum`)

    Last updated: Carlos Galdino 03-2023
    """
    # select type
    if 'type_' in kwargs: 
        type_ = kwargs.pop('type_')
    else:
        type_ = 'spectrum'

    # sort args
    filepath, folderpath, prefix, n, zfill = _sort_args(*args, **kwargs)

    # get data
    if filepath is not None:
        # xas data
        if type_.lower() in ['xas', 'x']:
            return _read_xas(filepath)
        # spectrum, photon events, bad events
        return _read_1(filepath, type_=type_)
    if folderpath is not None:
        if folderpath.is_dir() == False:
            raise AttributeError(f'It seems like the folderpath is not a folder\nfolderpath = {folderpath}')
        if prefix is None:
            raise AttributeError('Missing prefix.\nPlease, inform the file prefix.\nExample: for Cu_0001_d1.h5 the prefix is `Cu` or `Cu_`.')
        if n is None:
            raise AttributeError('Missing file number (scan number).\nPlease, inform the file number.\nExample: for Cu_0001_d1.h5 the file number is 1.')
        
        # xas data
        if type_.lower() in ['xas', 'x']:
            filepath = folderpath/(prefix+str(n).zfill(zfill)+'.xas')
            return _read_xas(filepath)
        
        # spectrum, photon events, bad events
        filepaths = [folderpath/(prefix+str(n).zfill(zfill)+f'_d{i}.h5') for i in (1, 2, 3)]

        if type_.lower() in ['spectrum', 's']:
            ss  = br.Spectra(3)
            for i, filepath in enumerate(filepaths):
                ss[i] = _read_1(filepath)
                ss[i].scan = n
                ss[i].ccd = i
            ss.scan = n
            # ss.nd   = [ss[0].nd, ss[1].nd, ss[2].nd]
            return ss
        elif type_.lower() in ['photon events', 'pe']:
            pes = [0]*3
            for i, filepath in enumerate(filepaths):
                pes[i] = _read_1(filepath, type_=type_)
                pes[i].scan = n
                pes[i].ccd = i
            pes.scan = n
            return pes
        elif type_.lower() in ['bad', 'b']:
            bads = [0]*3
            for i, filepath in enumerate(filepaths):
                bads[i] = _read_1(filepath, type_=type_)
                bads[i].scan = n
                bads[i].ccd = i
            bads.scan = n
            return bads
        else:
            raise ValueError("type_ must be one of 'spectrum', 'photon events', or 'bad'.")

def calib(folderpath, prefix, mode='cc', start_scan=None, stop_scan=None, scans=None, start_energy=None, stop_energy=None, energies=None):
    """Returns the calibration factor (dispersion) for ADRESS data.

    Args:
        folderpath (str or Path object): folderpath to files (for loading a data from 3 cdd files).
        prefix (str): file prefix.
        mode (string, optional): method used to calculate the shifts.
                The current options are: 'cross-correlation' ('cc'), 'max',
                'fitted peaks', or 'peak'. For mode='peak', mode will be 
                changed to 'fitted peaks'. If mode = 'fitted peaks', data will 
                be fitted with one peak using function br.Spectrum.fit_peak().
                Default is 'cc'.
        start_scan (number, optional): first scan number. Default is None.
        stop_scan (number, optional): last scan number. Default is None.
        scans (list, optional): list of scan numbers. Overwrites start_scan and
            stop_scan. Default is None.
        start_energy (number, optional): first energy value. Default is None.
        stop_energy (number, optional): last energy value. Default is None.
        energies (list, optional): list with energy values. It overwrites
            start_energy and stop_energy. Default is None. If start_energy and 
            stop_energy are also None, the photon energy values will be
            extracted from the file.   


    Example:

        >>> calib, sss = br.ADRESS.calib(folderpath=folderpath, prefix='Cu_', start_scan=19, stop_scan=29)

    Calibration values for each ccd (eV/bin):

        >>> print(calib)
    
    Spectra:

        >>> for ccd in (0, 1, 2):
        >>>     br.figure()
        >>>     _ = sss[ccd].plot()
        >>>     plt.title(f'ccd {ccd}: {np.round(calib[ccd], 4)} eV/bin')
        >>>     plt.legend(np.round(sss[ccd].E))
        >>>     plt.xlabel('Energy loss (eV)')
        >>>     plt.ylabel('Intensity (arb. units)')
        >>>     plt.grid()

    If mode is 'fitted peak', the fitted curves can be plotted as well:

        >>> for ccd in (0, 1, 2):
        >>>     br.figure()
        >>>     _ = sss[ccd].plot()
        >>>     _ = sss[ccd].fit.plot(color='red', lw=1)
        >>>     plt.title(f'ccd {ccd}: {np.round(calib[ccd], 4)} eV/bin')
        >>>     plt.legend(np.round(sss[ccd].E))
        >>>     plt.xlabel('Energy loss (eV)')
        >>>     plt.ylabel('Intensity (arb. units)')
        >>>     plt.grid()
    
    Shift (or peak center) as a function of energy:

        >>> for ccd in (0, 1, 2):
        >>>     br.figure()
        >>>     plt.title(f'ccd {ccd}: {np.round(calib1[ccd], 4)} eV/bin')
        >>>     sss[ccd].calculated_calib.plot(marker='o', lw=0, color='black')
        >>>     sss[ccd].calculated_calib.fit.plot(color='red', label='fit')
        >>>     plt.legend()
        >>>     plt.xlabel('Photon energy (eV)')
        >>>     plt.ylabel('Shift or peak center (bin)')
        >>>     plt.grid()

    Returns:
        2 lists with three elements. First list contains the calibration values 
            for each ccd and the second list contains the spectra.
    """
    # scan numbers
    if scans is None:
        scans = np.arange(start_scan, stop_scan+1)

    # get data
    sss = [br.Spectra(n=len(scans)), br.Spectra(n=len(scans)), br.Spectra(n=len(scans))]
    for i, scan in enumerate(scans):
        temp = read(folderpath=folderpath, prefix=prefix, n=scan, zfill=4)
        for ccd in (0, 1, 2):
            sss[ccd][i] = temp[ccd]

    # get energies
    if energies is None:
        if start_energy is None and stop_energy is None:
            energies = [np.mean(s.nd['PhotonEnergy']) for s in sss[0]]
        else:
            energies = np.linspace(start_energy, stop_energy, len(scans))
    assert len(scans) == len(energies), f'number of energies ({len(energies)}) do not match the number of scans ({len(scans)})'
    
    # save energies to spectra
    for ccd in (0, 1, 2):
        sss[ccd].E = energies
        for j in range(len(sss[ccd])):
            sss[ccd][j].E = energies[j]

    # calculate calib
    disp = [0, 0, 0]   
    if mode == 'peak':
        mode = 'fitted peaks'
    if mode == 'fitted peaks':
        for ccd in (0, 1, 2):
            sss[ccd].fit_peak()
    for ccd in (0, 1, 2):
        disp[ccd] = sss[ccd].calculate_calib(values=energies, mode=mode, bkg_check=False, deg=1)


    # dispersion as a function of detector size instead of bin
    # y_max = None
    # y_bin = None

    # if y_max is None:
    #     y_max = ss_temp[ccd].nd['ArraySizeY'][0]
    #     y_bin = len(ss_temp[ccd].x)

    # if y_max != ss_temp[ccd].nd['ArraySizeY'][0]:
    #     raise ValueError('detector size is not the same for all scans.')
    # if y_bin != len(ss_temp[ccd].x):
    #     raise ValueError('binning is not the same for all scans.')

    # disp_bin[ccd] = ss.calculate_calib(values=energies, mode=mode, peak=0)
    # disp[ccd] = disp_bin[ccd] / (y_max / y_bin)
    # ss.ccd = ccd
    # sss[ccd] = ss

    # return disp, disp_bin, sss
    return disp, sss

def final(folderpath, prefix, scan, calib=None, zero_mode=None, ref=-1, zfill=4):
    """[EXPERIMENTAL] Align and sum data from all ccd's.
    
    Args:
        filepath (str or Path object): filepath to h5 file (for loading a single file).
        folderpath (str or Path object): folderpath to files (for loading a data from 3 cdd files).
        prefix (str): file prefix.
        n (number): scan number. Example: for 'Cu_0001_d1.h5' the file number is 1.
            The number of extra zeros (0001) is defined by the argument `zfill`.
        calib (number, optional): calibration factor. If None, calibration is
            not performed. If list, the number of elements must be the same as 
            the number of ccd's. Each calib. factor will be applied to it's
            respective ccd in order. After calibration, data is interpolated, 
            aligned and summed up. If number, data is aligned first, then 
            calibrated, and finally summed up. Default is None.
        zero_mode (string, optional): 
        zfill (number, optional): number of digits to fill scan number (Default is 4).
    """
    ss = br.ADRESS.read(folderpath, prefix, scan)

    if calib is not None:
        if isinstance(calib, Iterable):
            assert len(ss) == len(calib), f'Number of calibration factor must match the number of spectra.\nnumber of calibration factors: {len(calib)}\nnumber of spectra: {len(ss)}'
            for i in range(len(ss)):
                ss[i].calib = calib[i]
            ss.interp()
            ss.align()
            s = ss.sum
        else:
            ss.align()
            s = ss.sum
            s.calib = calib
    else:
        ss.align()
        s = ss.sum

    # zero
    if zero_mode is not None:
        if zero_mode == 'fitted peaks':
            s.find_peaks()
            s.fit_peaks()
        elif zero_mode == 'peaks':
            s.find_peaks()
        s.zero(mode=zero_mode, ref=ref)

    # attr
    s.scan = scan

    return s
