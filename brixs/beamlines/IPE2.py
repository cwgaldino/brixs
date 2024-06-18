#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from IPE beamline - Sirius.

Last edited: Carlos Galdino 2024-06-12
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from collections.abc import Iterable
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import datetime
import warnings
import copy

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br
import brixs.addons.h5 as h5

# %% ------------------------- Special Imports ---------------------------- %% #
import h5py
# %%

# %% ========================= useful functions =========================== %% #
def scanlist(folderpath):
    """Return list of scans available in folderpath"""
    folderpath = Path(folderpath)
    assert folderpath.exists(), f'fpath does not exist ({folderpath})'
    return br.parsed_filelist(folderpath, string='.dat', ref=3, return_type='dict')
# %%

# %% =================== metadata support functions ======================= %% #
def _str2datetime(string):
    """convert IPE date/time string pattern to date --> '2022/07/20 21:08:36'

    Args:
        string (str): string with IPE date string

    Return
        datetime.datetime object
    """
    #########
    # split #
    #########
    if 'T' in string:
        date, time = string.split('T')
    else:
        date, time = string.split(' ')
    
    ########
    # date #
    ########
    year, month,  day = (int(_) for _ in date.split('/'))

    ########
    # time #
    ########
    hour, minute, seconds = time.split(':')

    hour    = int(hour)
    minute  = int(minute)
    seconds = int(round(float(seconds)))

    if seconds >= 60:
        minute  = minute + 1
        seconds = 0
    if minute >= 60:
        hour   = hour + 1
        minute = 0

    ############
    # datetime #
    ############
    return datetime.datetime(year=year, month=month, day=day, hour=hour, minute=minute, second=seconds)
# %%

# %% ========================== rixs metadata ============================= %% #
rixs_attrs = {'ignore': {}, 'raw':{}}

h = rixs_attrs['ignore']
h['modified_date'] = ''
# h['filename']      = ''

h = rixs_attrs['raw']
h['Energy']             = 'entry/instrument/NDAttributes/Energy'
h['Energy_SP']          = 'entry/instrument/NDAttributes/Energy_SP'
h['NDArrayEpicsTSSec']  = 'entry/instrument/NDAttributes/NDArrayEpicsTSSec'
h['NDArrayEpicsTSnSec'] = 'entry/instrument/NDAttributes/NDArrayEpicsTSnSec'
h['NDArrayTimeStamp']   = 'entry/instrument/NDAttributes/NDArrayTimeStamp'
h['NDArrayUniqueId']    = 'entry/instrument/NDAttributes/NDArrayUniqueId'
h['PGM_Cff']            = 'entry/instrument/NDAttributes/PGM_Cff'
h['PGM_GR']             = 'entry/instrument/NDAttributes/PGM_GR'
h['PGM_GT']             = 'entry/instrument/NDAttributes/PGM_GT'
h['PGM_MR']             = 'entry/instrument/NDAttributes/PGM_MR'
h['RIXS_Ry']            = 'entry/instrument/NDAttributes/RIXS_Ry'
h['RIXS_X']             = 'entry/instrument/NDAttributes/RIXS_X'
h['RIXS_Y']             = 'entry/instrument/NDAttributes/RIXS_Y'
h['RIXS_Z']             = 'entry/instrument/NDAttributes/RIXS_Z'
h['TempA']              = 'entry/instrument/NDAttributes/TempA'
h['TempB']              = 'entry/instrument/NDAttributes/TempB'
h['Undulator']          = 'entry/instrument/NDAttributes/Undulator'                     

# %% =========================== xas metadata ============================= %% #
# xas_attrs = {'ignore': {}, 'raw':{}, 'string':{}, 'round2':{}}                     
# %%

# %% ============================= read =================================== %% #
def _read(filepath, verbose=True):
    """return PhotonEvents list

    Args:
        filepath (str or path): filepath
        verbose (bool, optional): if True, warns when metadata cannot be read

    Returns:
        pe1, pe2 (photon events for each ccd)
    """
    ##################
    # check filepath #
    ##################
    filepath = Path(filepath)
    assert filepath.exists(), f'filepath does not exist ({filepath})'

    #############
    # open file #
    #############
    with h5py.File(Path(filepath), 'r') as f:
        #############
        # read data #
        #############
        data = (f['entry/data/data'][:])
        x    = data[:, 2]
        y    = data[:, 4]

        #############################
        # Create PhotonEvent object #
        #############################
        pe1 = br.PhotonEvents(x=x, y=y, xlim=(0, 1650),    ylim=(0, 1608)).crop(0,    1650, None, None)
        pe2 = br.PhotonEvents(x=x, y=y, xlim=(1651, 3264), ylim=(0, 1608)).crop(1651, None, None, None)

        #########
        # attrs #
        #########
        metadata = h5.sort_metadata(f=f, attrs_dict=rixs_attrs, verbose=verbose)
        for attr in metadata:
            setattr(pe1, attr, metadata[attr][0])
            setattr(pe1, attr, metadata[attr][0])

        # date
        temp = br.get_modified_date(filepath)
        setattr(pe1, 'modified_date', temp)
        setattr(pe2, 'modified_date', temp)

        # filename
        setattr(pe1, 'filename', filepath.name)
        setattr(pe2, 'filename', filepath.name)

        # ccd
        setattr(pe1, 'ccd', 1)
        setattr(pe2, 'ccd', 2)

    return pe1, pe2

def read(fpath, verbose=True):
    """Return data from folderpath

    Args:
        fpath (filepath or folderpath): filepath for xas and folderpath for RIXS
        verbose (bool, optional): Verbose, default is True.

    Returns:
        pe1, pe2, pe1's, pe2's for RIXS
        TEY, TFY, I0 for XAS
    """
    ####################
    # check folderpath #
    ####################
    fpath = Path(fpath)
    assert fpath.exists(), f'fpath does not exist ({fpath})'

    ########
    # RIXS #
    ########
    if fpath.is_dir():
        ##################
        # get each image #
        ##################
        filelist = br.parsed_filelist(dirpath=fpath, string='.h5', ref=3)
        assert len(filelist) > 0, f'no h5 files found in folderpath: {fpath}'
        
        # collect data
        x1 = list()
        y1 = list()
        x2 = list()
        y2 = list()
        dummy1 = br.Dummy()
        dummy2 = br.Dummy()
        for filepath in filelist:
            _pe1, _pe2 = _read(filepath, verbose)
            x1.extend(_pe1.x)
            y1.extend(_pe1.y)
            x2.extend(_pe2.x)
            y2.extend(_pe2.y)
            dummy1.append(_pe1)
            dummy2.append(_pe2)
        pe1 = br.PhotonEvents(x=x1, y=y1, xlim=_pe1.xlim, ylim=_pe1.ylim)
        pe2 = br.PhotonEvents(x=x2, y=y2, xlim=_pe2.xlim, ylim=_pe2.ylim)
        
        # attrs pe1
        for attr in _pe1.get_attrs():
            temp = [getattr(pe, attr) for pe in dummy1]
            setattr(dummy1, attr, temp)
            if attr not in ('modified_date', 'ccd', 'filename'):
                try:
                    setattr(pe1, attr, np.mean(temp))
                    setattr(pe1, attr + '_min', np.min(temp))
                    setattr(pe1, attr + '_max', np.max(temp))
                    setattr(pe1, attr + '_sigma', np.std(temp))
                except:
                    setattr(pe1, attr, None)
        pe1.ccd = 1
        pe1.modified_date = dummy1[0].modified_date

        # attrs pe2
        for attr in _pe2.get_attrs():
            temp = [getattr(pe, attr) for pe in dummy1]
            setattr(dummy2, attr, temp)
            if attr not in ('modified_date', 'ccd', 'filename'):
                try:
                    setattr(pe2, attr, np.mean(temp))
                    setattr(pe2, attr + '_min', np.min(temp))
                    setattr(pe2, attr + '_max', np.max(temp))
                    setattr(pe2, attr + '_sigma', np.std(temp))
                except:
                    setattr(pe2, attr, None)
        pe2.ccd = 2
        pe2.modified_date = dummy2[0].modified_date

        return pe1, pe2, dummy1, dummy2

    #######
    # XAS #
    #######
    else:
        #############
        # load file #
        #############
        try:
            data = br.load_data(filepath=fpath, force_array=True)
        except IndexError:
            raise ValueError(f'Error loading file {filepath}')
        
        #############
        # sort data #
        #############
        TEY = br.Spectrum(x=data[:, 0], y=data[:, 2])
        TFY = br.Spectrum(x=data[:, 0], y=data[:, 3])
        I0  = br.Spectrum(x=data[:, 0], y=data[:, 4])
        ss  = br.Spectra([TEY, TFY, I0])
        for s in [_ for _ in ss] + [ss]:
            s.PHASE        = data[:, 1]
            s.PD           = data[:, 5]
            s.DVF          = data[:, 6]
            s.SP_ENERGY    = data[:, 7]
            s.SP_PHASE     = data[:, 8]
            s.RING_CURRENT = data[:, 9]
            s.TIMESTAMP    = data[:, 10]

        #########
        # attrs #
        #########
        for line in br.load_comments(fpath):
            if line.startswith('#C # '):
                if '=' in line:
                    temp = line.split('#C # ')[1].split('=')
                    name  = temp[0].strip()
                    try:
                        value = float(temp[-1].strip())
                    except ValueError:
                        value = temp[-1].strip()
                    for s in [_ for _ in ss] + [ss]:
                        s.__setattr__(name, value)
            elif line.startswith('#S '):
                scan = line.split('#S ')[1].split()[0].strip()
                for s in [_ for _ in ss] + [ss]:
                    s.__setattr__('scan', scan)

                command = line.split('#S ')[1].split()[1:].strip()
                for s in [_ for _ in ss] + [ss]:
                    s.__setattr__('command', command)
            elif line.startswith('#D '):
                start_time = _str2datetime(line.split('#D ')[1].strip())
                for s in [_ for _ in ss] + [ss]:
                    s.__setattr__('start_time', start_time)
        return ss
    return
# %%

# %% ============================= RIXS =================================== %% #
def _process(fpath, curv='self', curv_nbins=(20, 1000), sbins=1200, calib=None):
    """
    """
    #############
    # read file #
    #############
    pe1, pe2, pes1, pes2 = read(fpath)

    # curvature
    if curv is not None:
        if type(curv) == str:
            if curv == 'self':
                # curvature correction with broken intervals
                fit, curv = pe1.curvature_correction_with_broken_intervals(curv_nbins=curv_nbins)
                s, fit, popt, R2, model = pe1.calculate_vertical_shift_curvature(ncols, nrows, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):

                # curvature correction for solid intervals (works fine)
                # im = pe3.binning(ncols=curv_nbins[1], nrows=curv_nbins[0])
                # limits = [[m[2], m[3]] for m in mask]
                # _s, fit, popt, R2, model = im.calculate_horizontal_shift_curvature(deg=2, mode='cc', limits=limits)
                # curv = popt

        pe4 = pe3.set_horizontal_shift_via_polyval(p=curv)
    else:
        pe4 = None
    
    # spectrum
    if curv is not None and isinstance(curv, br.Iterable):
        s = pe4.integrated_columns_vs_x_centers(ncols=sbins)
        s.scan = scan

        # normalization
        # s = s.set_factor(1/sum([m[3]-m[2] for m in mask]))
        s = s.set_factor(1/s.exposure_time)
        s = s.set_factor(sbins)
        s = s.set_factor(1000)
    else:
        s = None
        
    # calib ==============================================
    if calib is not None and s is not None:
        s.calib = calib
        s.shift = -s.E

    return {'pe1':pe1, 'pe2':pe2, 'pes1':pes1, 'pes2':pes2,
            'curv':curv, 
            'curv_fit':fit,
            's':s}

def verify(filepath, scan, mask, tcutoff=3e7, tnbins=10000, period=1458, offset=None, twidth=320, tcenter='max', curv='self', curv_nbins=(20, 1000), sbins=1200, calib=None, maximize=True, left=0.04, right=0.99, top=0.97, bottom=0.05):
    """open a figure with step-by-step rixs data reduction

    Args:

    
    Returns:
        pe, pe2, pe3, pe4, curv, s
    """
    ##### process ######
    d = _process(filepath=filepath, scan=scan, mask=mask, tcutoff=tcutoff, tnbins=tnbins, period=period, offset=offset, twidth=twidth, tcenter=tcenter, curv=curv, curv_nbins=curv_nbins, sbins=sbins, calib=calib)
    pe                   = d['pe']
    cutoff               = d['cutoff']
    pe2                  = d['pe2']
    time_bunch_histogram = d['time_bunch_histogram']
    folded_time          = d['folded_time']
    tmask                = d['tmask']
    folded_tmask         = d['folded_tmask']
    bad                  = d['bad']
    pe3                  = d['pe3']
    pe4                  = d['pe4']
    curv_fit             = d['curv_fit']
    s                    = d['s']


    ###### figure ######
    fig, axes = br.subplots(2, 5, figsize=(12, 4))
    plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    if maximize:
        br.maximize()

    ###### axes 0 ######
    ax = axes[0]
    ax.set_title('1: Raw (pe1)')
    # plot pe
    pe.plot(ax, color='black')  
    if tcutoff is not None:
        cutoff.plot(ax, color='magenta', s=.6)
    # plot mask
    for m in pe2.mask:
        br.rectangle(m[:2], m[2:], ax=ax, lw=2, edgecolor='red')
    
    ###### axes 1 ######
    ax = axes[1]
    ax.set_title('2: After mask (pe2)')
    # plot pe2
    pe2.plot(ax, color='dodgerblue')   
    
    ###### axes 2 ######
    ax = axes[2]
    ax.set_title('4a: Time-rejected photons (bad)')
    bad.plot(ax, color='lightcoral')

    ###### axes 3 ######
    ax = axes[3]
    ax.set_title('4b: Time-accepted photons (pe3)')
    pe3.plot(ax, color='green')
    # plot fit
    if pe4 is not None:
        temp  = pe3.binning(nrows=20, ncols=1000).rows[0]
        shift = temp.x[np.argmax(temp.y)]
        curv_fit.switch_xy().flip_x().plot(ax, offset=0, shift=shift, color='red')
 

    ###### axes 4 ######
    ax = axes[4]
    if pe4 is not None:
        ax.set_title('5: curvature corrected (pe4)')
        pe4.plot(ax, color='red')

    ###### axes 5 ######
    ax = axes[5]
    ax.set_title('3a: Folded Time (folded_time)')
    folded_time.plot(ax)
    for m in folded_tmask:
        br.rectangle(m, (0, max(folded_time.y)), ax=ax, lw=2, edgecolor='red')
    
    ###### axes 6 ######
    ax = axes[6]
    ax.set_title('3b: Begining of unfolded Time')
    time_bunch_histogram.plot(ax, color='black', lw=.2, marker='o', ms=1)
    for m in tmask:
        i = time_bunch_histogram.index(m[0])
        f = time_bunch_histogram.index(m[1])
        br.rectangle(m, (0, max(time_bunch_histogram.y[i:f])), ax=ax, lw=2, edgecolor='red')
    if offset is None: offset = 0
    br.zoom(min(time_bunch_histogram.x)-min(time_bunch_histogram.x)*-0.1, period*5+offset, ax)        

    # ###### axes 7 ######
    # ax = axes[7]
    # ax.set_title('3b: End of unfolded Time (time_bunch_histogram)')
    # time_bunch_histogram.plot(ax, color='black', lw=.2, marker='o', ms=1)
    # for m in tmask:
    #     i = time_bunch_histogram.index(m[0])
    #     f = time_bunch_histogram.index(m[1])
    #     br.rectangle(m, (0, max(time_bunch_histogram.y[i:f])), ax=ax, lw=2, edgecolor='red')
    # if offset is None: offset = 0
    # br.zoom(max(time_bunch_histogram.x) - period*5+offset, max(time_bunch_histogram.x)+max(time_bunch_histogram.x)*0.1, ax)  

    ###### axes 7 ######
    ax = axes[7]
    ax.set_title('3b: Full unfolded Time (time_bunch_histogram)')
    time_bunch_histogram.plot(ax, color='black', lw=.2, marker='o', ms=1)
    for m in tmask:
        i = time_bunch_histogram.index(m[0])
        f = time_bunch_histogram.index(m[1])
        br.rectangle(m, (0, max(time_bunch_histogram.y[i:f])), ax=ax, lw=2, edgecolor='red')

    ###### axes 8 ######
    ax = axes[8]
    ax.set_title('3b: Photons per time window')
    get_count_per_time_window(pe2.time_bunch, tmask).plot(ax, color='black', marker='o', ms=2)

    ###### axes 9 ######
    ax = axes[9]
    ax.set_title('6: Final spectrum')
    if s is not None:
        s.plot(ax, color='black', marker='o', ms=2)
    
    ###### axis labels ######
    for i in (0, 1, 2, 3, 4):
        axes[i].set_xlabel('x (mm)')
        axes[i].set_ylabel('y (mm)')
    axes[5].set_xlabel('Folded time (ps)')
    axes[5].set_ylabel('Photon count')
    for i in (6, 7):
        axes[i].set_xlabel('Time (ps)')
        axes[i].set_ylabel('Photon count')
    axes[8].set_xlabel('Time (ps)')
    axes[8].set_ylabel('Photon count per window')
    if calib is None:
        axes[9].set_xlabel('x (mm)')
    else:
        axes[9].set_xlabel('Energy loss (eV)')
    axes[9].set_ylabel('Intensity (arb. units)')

    return {'pe':pe, 'cutoff':cutoff, 'pe2':pe2, 'time_bunch_histogram':time_bunch_histogram, 'folded_time':folded_time, 'tmask':tmask, 'folded_tmask':folded_tmask, 'bad':bad, 'pe3':pe3, 'pe4':pe4, 'curv':curv, 's':s}

@br.finder.track
def process(filepath, scan, mask, tcutoff=3e7, tnbins=10000, period=1458, offset=None, twidth=320, tcenter='max', curv='self', curv_nbins=(20, 1000), sbins=1200, calib=None):
    """
    """
    d = _process(filepath=filepath, scan=scan, mask=mask, tcutoff=tcutoff, tnbins=tnbins, period=period, offset=offset, twidth=twidth, tcenter=tcenter, curv=curv, curv_nbins=curv_nbins, sbins=sbins, calib=calib)
    return d['s']
# %%
