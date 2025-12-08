#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced functions for VERITAS beamline at MAX-IV"""

# %% ------------------------- Standard Imports --------------------------- %% #
import matplotlib.pyplot as plt
from itertools import compress
import numpy as np
import copy

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br
import brixs.addons.fitting
from brixs.beamlines.veritas.core import read
# %%

# %% ============================= RIXS =================================== %% #
def _process(scan, filepath, mask=None, tcutoff=3e7, tnbins=10000, period=1458, offset=None, twidth=320, tcenter='max', curv=None, curv_nbins=(20, 1000), sbins=1200, calib=None, verbose=True):
    """internal process function.

    For a complete description of this function, Please refer to 
    brixs/examples/beamlines/veritas/rixs.py

    TODO:
        For the future we could put an argument like, time_start, time_finish 
        to control if we will read all photons or just photos at the beginning
        or end of a scan. Should be easy if we convert from ps to s. We can do 
        just like "time cutoff".

    """
    #############
    # read file #
    #############
    _pe = read(filepath=filepath, scan=scan)

    ###############
    # time cutoff #
    ###############
    if tcutoff is not None:
        pe, cutoff = _pe.apply_time_cutoff(value=tcutoff)
    else:
        pe = _pe
        cutoff = None

    ########
    # mask #
    ########
    if mask is None:
        mask = [[pe.xlim[0], pe.xlim[1], pe.ylim[0], pe.ylim[1]], ]
    else:
        mask = pe._check_mask(mask)
    pe2      = pe.clip2(mask)
    pe2.mask = copy.deepcopy(mask)

    ###############
    # folded time #
    ###############
    time_bunch_histogram = pe2.calculate_time_bunch_histogram(nbins=tnbins)  # returned spectrum has attr: nbins
    folded_time = time_bunch_histogram.calculate_folded_time(period=period, offset=offset)  # returned spectrum has attrs: nbunchs and period

    ##########
    # t mask #
    ##########
    try:
        folded_tmask, tmask = calculate_folded_tmask_and_tmask(folded_time, time_bunch_histogram, w=twidth, c=tcenter)
    except AssertionError:
        print('error finding suitable c for time window')
        folded_tmask, tmask = calculate_folded_tmask_and_tmask(folded_time, time_bunch_histogram, w=twidth, c='max')

    #########################
    # time rejected photons #
    #########################
    pe3, bad = pe2.apply_tmask(tmask=tmask)
    pe3.mask = copy.deepcopy(mask)

    # curvature
    if curv is not None:
        if type(curv) == str:
            if curv == 'self':
                # get separate PhotonEvents'
                pes = br.Dummy()
                for m in mask:
                    pes.append(pe.clip2(m))

                # curvature correction with broken intervals
                fit, curv = pes.curvature_correction_with_broken_intervals(curv_nbins=curv_nbins)
                
                # curvature correction for solid intervals (works fine)
                # im = pe3.binning(ncols=curv_nbins[1], nrows=curv_nbins[0])
                # limits = [[m[2], m[3]] for m in mask]
                # _s, fit, popt, R2, model = im.calculate_horizontal_shift_curvature(deg=2, mode='cc', limits=limits)
                # curv = popt
            else:
                raise ValueError('curv can only be None, a list, or "self"')
        else:
            fit = None

        pe4 = pe3.set_horizontal_shift_via_polyval(p=curv)
    else:
        pe4 = None
        fit = None
    
    # spectrum
    if curv is not None:
        s = pe4.integrated_columns_vs_x_centers(ncols=sbins)
    else:
        s = pe3.integrated_columns_vs_x_centers(ncols=sbins)
    s.scan = scan

    # normalization
    s = s.set_factor(1/sum([m[3]-m[2] for m in mask]))
    try:
        if _pe.exposure_time > 0:
            s = s.set_factor(1/_pe.exposure_time)
        else:
            pass
            # 
    except AttributeError:
        if verbose: print(f'Exposure time not found for scan {scan}. Skipping exposure time normalization')
    s = s.set_factor(sbins)
    # s = s.set_factor(1000)
        
    # calib ==============================================
    if calib is not None:
        s.calib = calib
        # s.shift = -s.E

    # recover attrs
    for attr in _pe.get_attrs():
        s.__setattr__(attr, _pe.__getattribute__(attr))

    return {'pe':pe, 'cutoff':cutoff, 'pe2':pe2, 
            'time_bunch_histogram':time_bunch_histogram, 
            'folded_time':folded_time, 
            'tmask':tmask, 
            'folded_tmask':folded_tmask, 
            'bad':bad, 
            'pe3':pe3, 
            'pe4':pe4, 
            'curv':curv, 
            'curv_fit':fit,
            's':s}

def verify(scan, filepath, mask, tcutoff=3e7, tnbins=10000, period=1458, offset=None, twidth=320, tcenter='max', curv=None, curv_nbins=(20, 1000), sbins=1200, calib=None, figsize=(12, 4), verbose=True):
    """Open a figure with step-by-step rixs data reduction

    For a complete description of this function, Please refer to 
    brixs/examples/beamlines/veritas/rixs.py

    Returns:
        dictionary
    """
    ##### process ######
    d = _process(filepath=filepath, scan=scan, mask=mask, tcutoff=tcutoff, tnbins=tnbins, period=period, offset=offset, twidth=twidth, tcenter=tcenter, curv=curv, curv_nbins=curv_nbins, sbins=sbins, calib=calib, verbose=verbose)
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
    fig, axes = br.subplots(2, 5, figsize=figsize, layout='constrained')

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
        if curv_fit is not None:
            curv_fit.switch_xy().flip_x().plot(ax, offset=0, shift=shift, color='red')
 

    ###### axes 4 ######
    ax = axes[4]
    ax.set_title('5: curvature corrected (pe4)')
    if pe4 is not None:
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

def process(scan,filepath, mask, tcutoff=3e7, tnbins=10000, period=1458, offset=None, twidth=320, tcenter='max', curv=None, curv_nbins=(20, 1000), sbins=1200, calib=None, verbose=True):
    """return processed rixs spectrum

    For a complete description of this function, Please refer to 
        brixs/examples/beamlines/veritas/rixs.py
    
    Returns:
        spectrum
    """
    d = _process(filepath=filepath, scan=scan, mask=mask, tcutoff=tcutoff, tnbins=tnbins, period=period, offset=offset, twidth=twidth, tcenter=tcenter, curv=curv, curv_nbins=curv_nbins, sbins=sbins, calib=calib, verbose=verbose)
    return d['s']
# %%

# %% ========================== Time support ============================== %% #
def _apply_time_cutoff(self, value=3e7):
    """return good and cutoff photon events

    Args:
        Value (number, optional): cutoff time. Photon events with time higher then
            this value will be removed
    
    Returns:
        good, cutoff
    """
    # cutoff time
    indexes = np.array([i for i in range(len(self.time_bunch)) if self.time_bunch[i] > value]) # remove this one last point
    final  = [True]*len(self.time_bunch)
    for i in indexes:
        final[i] = False
    _cutoff = [not elem for elem in final]
    
    # saving
    cutoff = br.PhotonEvents(x=self.x[_cutoff], y=self.y[_cutoff])
    cutoff.copy_attrs_from(self)
    cutoff.xlim = self.xlim
    cutoff.ylim = self.ylim
    cutoff.time_absolute = self.time_absolute[_cutoff]
    cutoff.time_bunch    = self.time_bunch[_cutoff]
    
    good = br.PhotonEvents(x=self.x[final], y=self.y[final])
    good.copy_attrs_from(self)
    good.xlim = self.xlim
    good.ylim = self.ylim
    good.time_absolute = self.time_absolute[final]
    good.time_bunch    = self.time_bunch[final]

    return good, cutoff
br.PhotonEvents.apply_time_cutoff = _apply_time_cutoff    

def _calculate_time_bunch_histogram(self, nbins=10000):
    """returns the histogram of time_bunch

    Args:
        nbins (int, optional): number of bins 

    Return:
        Spectrum
        returned spectrum has attr: nbins
    """
    # calculate time histogram
    hist, bin_edges = np.histogram(self.time_bunch, nbins)
    time            = br.Spectrum(x=br.moving_average(bin_edges, 2), y=hist)
    
    # self.time.shift = -min(self.time.x)
    time.nbins = nbins
    # self.time_bunch_histogram       = time
    # self.time_bunch_histogram_nbins = nbins
    return time
br.PhotonEvents.calculate_time_bunch_histogram = _calculate_time_bunch_histogram        

def _calculate_folded_time(self, period=1459, offset=None):
    """return folded time

    Args:
        period (int, optional): time period to fold. After this value, time is supposed to repeat itself
        offset (number, optional): use `offset` to shift time before applying period, 
            in case time does not start from 0
    
    Return:
        Spectrum
        returned spectrum has attrs: nbunchs and period
    """    
    if offset is None:
        offset = min(self.x)
    x = (self.x - offset)%period
    
    hist, bin_edges = np.histogram(x, bins=int(period/10), weights=self.y)
    folded_time     = br.Spectrum(x=br.moving_average(bin_edges, 2), y=hist)
    # folded_t        = (self.t-min(self.t))%period  # min(self.t) should not be necessary because this is already set to zero during read()
    folded_time.nbunchs = int(round((max(self.x) - offset)/period))
    folded_time.period  = period
    return folded_time
br.Spectrum.calculate_folded_time = _calculate_folded_time

def calculate_folded_tmask_and_tmask(folded_time, time_bunch_histogram, w=300, c='max'):
    """return a suitable time mask from folded time (a peak -bunch- is expected)

    Args:
        folded_time (br.Spectrum): folded time calculated via br.Spectrum.calculate_folded_time() function. 
            Spectrum must have attr `period` used to calculate the folding
        time_bunch_histogram (br.Spectrum): time bunch histogram calculated via
            br.PhotonEvents.calculate_time_bunch_histogram() function
        w (number, optional): size of the mask
        c (number or str, optional): time center of the bunch or method for guessing the center. 
            options for guessing are 'gauss' and 'max'.
        
    Return:
        folded_tmask, tmask
    """
    # w must be less or equal the full time range
    assert w <= max(folded_time.x) - min(folded_time.x), f'w ({w}) must be equal or less than full time range ({max(self.x) - min(self.x)})'

    # get folded_time peak center
    if c == 'gauss':
        # check if fitting was imported
        if hasattr(folded_time, 'fit_peak') == False and callable(folded_time.fit_peak) == False:
            raise ValueError('cannot calculate shifts via `peaks` because fitting functions are not imported\nPlease import fitting function via `import brixs.addons.fitting`')
        _result = folded_time.fit_peak(fixed_m=0)
        popt = _result['popt']
        c = popt[1]
    elif c == 'max':
        c = folded_time.x[np.argmax(folded_time.y)]
    
    # tmask
    right = max(folded_time.x)
    left  = min(folded_time.x)
    assert c >= left and c <= right, f'c={c} outside range'
    
    # find if mask is outside
    if c+w/2 > right: 
        rest = w - (right - (c-w/2))
        folded_tmask = [(left, left+rest), (c-w/2, right)]
    elif c-w/2 < left: 
        rest = w - ((c+w/2) - left)
        folded_tmask = [(left, c+w/2), (right-rest, right)]
    else:
        folded_tmask = [(c-w/2, c+w/2), ]

    # calculate tmask
    tmask = [(c1-w/2+c, c1+w/2+c) for c1 in np.arange(min(time_bunch_histogram.x), max(time_bunch_histogram.x), folded_time.period)][:-1]

    return folded_tmask, tmask

def _apply_tmask(self, tmask):
    """Return photon events with point only within tmask

    Args:
        tmask (list): time mask with the following format: ((t1_start, t1_stop), (t2_start, t2_stop), ...)

    Note:
        tmask is applied on the regular time (not folded time)

    Return:
        good, bad
        PhotonEvent's have attr: percentage
    """
    indexes = []
    for m in tmask:
        indexes += [i for i, _t in enumerate(self.time_bunch) if _t>=m[0] and _t<=m[1]]
    final = [False]*len(self.time_bunch)

    for i in indexes:
        final[i] = True
    opposite = [not elem for elem in final]

    # final
    good = br.PhotonEvents(x=list(compress(self.x, final)), y=list(compress(self.y, final)))
    good.time_bunch    = np.array(self.time_bunch)[final]    #- min(self.t[final])    # set min value to zero
    good.time_absolute = np.array(self.time_absolute)[final]    #- min(self.t[final])    # set min value to zero
    good.copy_attrs_from(self)
    good.xlim = self.xlim
    good.ylim = self.ylim
    good.percentage = round(len(good.x)/len(self.x)*100, 2)

    bad  = br.PhotonEvents(x=list(compress(self.x, opposite)), y=list(compress(self.y, opposite)))
    bad.time_bunch  = np.array(self.time_bunch)[opposite] #- min(self.t[opposite]) # set min value to zero
    good.time_absolute   = np.array(self.time_absolute)[final]    #- min(self.t[final])    # set min value to zero
    good.copy_attrs_from(self)
    good.xlim = self.xlim
    good.ylim = self.ylim
    bad.percentage  = round(len(bad.x)/len(self.x)*100, 2)
    
    return good, bad
br.PhotonEvents.apply_tmask = _apply_tmask

def get_count_per_time_window(time_bunch, tmask):
    centers = []
    counts  = []
    for m in tmask:
        centers += [m[0] + (m[1] - m[0])/2]
        counts += [sum([1 for t in time_bunch if t > m[0] and t < m[1]])]            
    return br.Spectrum(centers, counts)

def _clip(self, mask):
    """Return a masked copy of the object (VERITAS).

    Note:
        This is different from pe.clip() because this also clips the time.

    Args:
        mask (list): list with rectangular coordinates `(x_start, x_stop, y_start, y_stop)`
            or a list with multiple rectangular coordinates, i.e., `[(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])`

    Returns:
        :py:attr:`PhotonEvents`
    """
    return self.clip(mask=mask, attrs2clip=['time_bunch', 'time_absolute'])
br.PhotonEvents.clip2 = _clip
# %%

# %% ====================== curvature correction ========================== %% #
def _curvature_correction_with_broken_intervals(self, curv_nbins):
    # binning PhotonEvents'
    ims = br.Dummy()
    for i in range(len(self)):
        ims.append(self[i].binning(ncols=curv_nbins[1], nrows=int(curv_nbins[0]/len(self))))

    # collect rows
    ss = br.Spectra()
    y_centers = []
    for i in range(len(self)):
        for y in ims[i].y_centers:
            y_centers.append(y)
        for s in ims[i].rows:
            ss.append(s)

    # calculate shift
    try:
        ss.check_same_x()
    except ValueError:
        ss = ss.interp()
    values = ss.calculate_shift(mode='cc')
    _s = br.Spectrum(y_centers, values)
    _result = _s.polyfit(deg=2)
    fit = _result['fit']
    popt = _result['popt']
    curv = popt

    return fit, curv
br.Dummy.curvature_correction_with_broken_intervals = _curvature_correction_with_broken_intervals

def _find_empty(self, mask, axis=0):
    """Find pixel cols or rows that are outside of mask
    
    args:
        mask (list or None): only photons within the limits of this mask are 
            considered. Mask must be the following format 
            ([x1_start, x1_stop, y1_start, y1_stop], [x2_start, x2_stop, y2_start, y2_stop], ...]
            if None, mask is not applied. Default is None.
        axis (0 or 1): If axis=0, it looks for empty rows. If axis=1, it looks
            for empty cols. Default is 0.
    
    Returns:
        list of booleans. True means within the mask. False, outside of the mask
    """
    if axis == 1:
        x = self.x_centers
    if axis == 0:
        x = self.y_centers
    else:
        raise ValueError('axis must be 0 or 1')
        
    empty = [True]*len(x)
    for i, x in enumerate(x):
        for r in mask:
            if (x > r[2] and x < r[3]):
                empty[i] = False
                
    return empty
br.Image.find_empty = _find_empty

def _get_filtered_spectra(self, mask, axis=0):
    """Return rows/cols that are within mask
    
    args:
        mask (list or None): only photons within the limits of this mask are 
            considered. Mask must be the following format 
            ([x1_start, x1_stop, y1_start, y1_stop], [x2_start, x2_stop, y2_start, y2_stop], ...]
            if None, mask is not applied. Default is None.
        axis (0 or 1): If axis=0, it looks for empty rows. If axis=1, it looks
            for empty cols. Default is 0.
    
    Returns:
        Spectra
    """
    empty = self.find_empty(mask=mask, axis=axis)
    
    assert axis == 1 or axis == 0, 'axis must be 0 or 1'
    if axis == 0:
        ss = self.rows
        x  = self.y_centers
    else:
        ss = self.columns
        x  = self.x_centers
    
    filtered = br.Spectra()
    for i, is_empty in enumerate(empty):
        if not is_empty:
            s = ss[i]
            s.center = x[i]
            filtered.append(s)
            
    filtered.centers = [s.center for s in filtered]
    return filtered
br.Image.get_filtered_spectra = _get_filtered_spectra

def verify_curvature_correction(scan, filepath, nbins=(20, 1000), axis=0, mask=None, tcutoff=3e7, tnbins=10000, period=1458, offset=None, twidth=320, tcenter='max', figsize=(12, 4), verbose=True):
    
    d = _process(filepath=filepath, scan=scan, mask=mask, tcutoff=tcutoff, tnbins=tnbins, period=period, offset=offset, twidth=twidth, tcenter=tcenter, curv=None, verbose=verbose)
    pe3 = d['pe3']

    # initiate figure
    fig, axes = br.subplots(2, 3, sharex=True, figsize=figsize, layout='constrained')
    fig.suptitle('Cuvature correction')

    # plot photon events (raw)
    pe3.plot(axes[0], color='black')

    # binning
    im = pe3.binning(ncols=nbins[1], nrows=nbins[0])

    # get rows that are not empty
    ss = im.get_filtered_spectra(mask=pe3.mask, axis=axis)

    # Calculate shifts
    calculated_shift = ss.calculate_shift()

    # fit shifts
    s = br.Spectrum(x=ss.centers, y=calculated_shift)
    _result = s.polyfit(2)
    fit   = _result['fit']
    curv  = _result['popt']
    model = _result['model']
    # print(f'curv = {list(popt)}')
    x   = np.linspace(min(pe3.x), max(pe3.x), 100)
    fit = br.Spectrum(x=x, y=model(x))

    # plot reduced image
    pos = im.plot(axes[1])

    # plot photon events (raw) with vertical bins
    pe3.plot(axes[2], color='black')
    if axis == 0:
        _ = br.axhlines(ax=axes[2], y=pos.y_edges, color='red', lw=.5)
    else:
        _ = br.axvlines(ax=axes[2], x=pos.x_edges, color='red', lw=.5)

    # plot horizontal integration of each vertical bin
    # ss2 = ss.copy()
    # ss2.switch()
    # ss2.flip()
    _ = ss.plot(axes[4])

    # get max and min y (makes plot nicer)
    if axis == 0:
        _result = ss[0].fit_peak()
        popt = _result['popt']
        model = _result['model']
        xmin = popt[1] - popt[2]*5

        _result = ss[-1].fit_peak()
        popt = _result['popt']
        model = _result['model']
        xmax = popt[1] + popt[2]*5
        for i, ax in enumerate(axes):
            ax.set_xlim(xmin, xmax)
            if i != 4:
                ax.set_ylim(min(pe3.y), max(pe3.y))
    else:
        cols2 = im.columns
        _result = cols2[0].fit_peak()
        popt = _result['popt']
        model = _result['model']
        ymin = popt[1] - popt[2]*5

        _result = cols2[-1].fit_peak()
        popt = _result['popt']
        model = _result['model']
        ymax = popt[1] + popt[2]*5
        for ax in axes:
            ax.set_ylim(ymin, ymax)

    # plot fitting
    pe3.plot(axes[3], color='black')
    if axis == 0:
        fit2 = fit.copy()
        fit2 = fit2.set_factor(-1)
        fit2 = fit2.switch_xy()
        fit2 = fit2.set_shift(ss[0].x[np.argmax(ss[0].y)])
    else:
        raise NotImplementedError('not implemented yet')
    fit2.plot(axes[3], color='red')

    # set shifts
    pe4 = pe3.copy()
    if axis == 0:
        pe4 = pe4.set_horizontal_shift_via_polyval(p=curv)
    else:
        pe4 = pe4.set_vertical_shift_via_polyval(p=curv)
    pe4.plot(axes[5], color='red')
    
    # spectrum ==========================================
    # s = pe4.calculate_spectrum(axis=axis, nbins=sbins)
    # s.set_factor(1/sum([m[3]-m[2] for m in pe3.mask]))
    # s.set_factor(1/s.exposure)
    # s.set_factor(sbins)
    # s.t = None

    return curv
# %%


