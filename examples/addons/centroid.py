#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""The centroid algorithm returns photon hit positions with sub pixel resolution

This is an approximate description of the algorithm:
    1) enhance image (only if enhance is True)
    2) get all positions where pixel intensity is above threshold
    3) for positions that are too close (<=n), keep only the one for the brightest pixel
    4) go back to the original image and find brightest pixel around the positions
    suggested by the enhanced image (<=n) 
    5) if two positions are too close (<=_n2), assign them as double events (only 
    works if threshold2 is not None)
    6) if floor is True, floor the image. See floor argument below.
    7) calculate center of mass of position (<=_n2)
    8) if the pixel intensity of a position is higher than threshold2, assign
    it as double event (only works if threshold2 is not None)

The main function of this implementation is the im.centroid(). This methods 
returns a list of photon hits defined by a threshold. See below the arguments of
the centroid method,

    Args:
        n (int): photon hits candidates that are within n pixels of distance 
            from each other will be considered the same candidate. For better 
            results, set n to be roughly the expected pixel distance between a 
            photon hit and the farthest excited pixel.
        threshold (number): threshold value. Any pixel value above threshold 
            will be considered as a photon hit candidate. If enhance=True, 
            threshold must be given in terms of pixel intensity of the 
            enhanced image [use im.enhance(n=_n, bkg=_bkg) to get enhanced image]. 
        threshold2 (number, optional): upper limit for threshold. If 
            threshold2 is not None, pixel values
            above threshold2 will be considered double events. If enhance=True, 
            threshold2 must be given in terms of pixel intensity of the 
            enhanced image [use im.enhance(n=_n, bkg=_bkg) to get enhanced image]. 
            Note that, even if threshold2 is None, two photon hit candidates 
            that are closer than _n2 pixels
            are still going to be considered double events. Default is None.

        floor (bool, optional): if True, an intensity offset is added to the 
            image before calculating the the center of masses such as the 
            average intensity of the whole image is zero. Default is True. This 
            is unnecessary if the image is already originally floored. Center of
            mass calculation can yield less precise results if image is not floored and
            _n2 >> n.

        enhance (bool, optional): If True, the image will be floored (an offset
            will be applied so avg pixel intensity is zero), squared, 
            and a moving averaged of size n will be applied. Default is True.
        _n (int): Use this to overwrite the size of the moving average window 
            used in the process of enhancing the image, i.e., number of points
            to average. If this is None, the averaging window size will be set 
            to n*2+1. If enhance is False, this has no effect.
        _bkg (number, optional): Use this to overwrite the bkg value for flooring
            the image. If _bkg is not None, the image will be subtracted by _bkg.
            If _bkg is None, bkg will be defined so the average of the whole 
            image is zero. Default is None. If enhance is False, this has no effect.

        _n2 (int): Use this to overwrite the number of neighbors to include for
            calculating the center of mass of a photon hit candidate, e.g., if 
            _n2=1, only first neighbors. _n2 also defines how close two 
            candidates need to be to be considered a double event (only valid
            if threshold2 is not None). If None, _n2 will be same as n. Default 
            is None.     
        spot_zeroing_type (str): when calculating the center of mass of a 
            candidate, pixels around the candidate cannot be negative, otherwise
            center of mass calculation can yield to less precise result (the
            result of the center of mass calculation can even be a detector
            position "outside" the range of the detector). Therefore, if a 
            negative pixel intensity is present around a candidate there is 
            two thing we can do 1) set it to zero (which makes sense because
            if it is negative, we assume it is close to zero), or 2) we apply a 
            small intensity offset to the whole spot around the candidate to 
            make all pixels around it positive or zero. Use
            spot_zeroing_type='zero' for 1 and 'offset' for 2. Default is 'zero'
        
        MAX_NUMBER_OF_CANDIDATES (int, optional): raises error if number of 
            photons to be patched out is larger than MAX_NUMBER_OF_CANDIDATES.
            Useful for preventing too low threshold as patching image is slow.
            Default is 10. 
"""

# %% ========================== standard imports ========================= %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# %% =========================== brixs imports =========================== %% #
import brixs as br
import brixs.beamlines.id32 as id32
import brixs.addons.centroid

# %% ============================== settings ============================= %% #
# matplotlib (optional)
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()

# brixs (optional)
# br.settings.FIGURE_POSITION = (283, 567)
# br.get_window_position()
# %%


# %  ===================================================================== %% #
# %  ======================== ID32 beamline example ====================== %% #
# %% ===================================================================== %% #
TOP = Path(r'C:\Users\galdin_c\github\brixsexampledata\beamlines\id32')

# %% Spike finding and patching (step-by-step)
sample  = 'align_rixs'
dataset = 'align_0001'
scan    = 125
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)
im = ims[1]

# defining 
n   = 6
_n  = 4
_n2 = n
threshold = 1e4
_patch_size = _n2
# n: photon hits candidates that are within n pixels of distance 
# from each other will be considered the same candidate. For better 
# results, set n to be roughly the expected pixel distance between a 
# photon hit and the farthest excited pixel.

# _n: moving average window used in the process of enhancing the image, i.e., 
# number of points to average. For optimal results, _n should be roughly the
# same size of the spike, e.g., if a spike lights up 4x4 grid of pixels, then _n=4
# should be optimal. If _n is too high, the image is washed out. If _n is too 
# little, the bkg does not go to zero so fast and spikes do not stand out as much.

# _n2: number of neighbors to include for calculating the center of mass of a 
# photon hit candidate. e.g., if _n2=1, only first neighbors. Calculating the 
# center of mass for spike removal is not necessary though, here we are 
# calculating it just for reference.

# _patch_size: patch size will be a square of side n3+1.

# threshold: intensity threshold for spike candidates (in this case, note that
# threshold must be in terms pixel intensity in the enhanced image). 
# defining the threshold for detecting spikes is an iterative process by 
# inspecting this image we can now go back and tweak threshold. Note also that 
# one must use multiple images to define this threshold.

# enhance spikes
im2 = im.floor()
im3 = im2.multiply(im2).moving_average(_n)

# plot
fig, axes = br.subplots(1, 3, sharex=True, sharey=True, figsize=(40, 12), layout='constrained')
im2.pcolormesh(ax=axes[0])
im3.pcolormesh(ax=axes[1])

# white dots are is the brightest pixel above threshold in the enhanced image around 6 pixels of distance
pos, _ = im3.get_positions_above_threshold(threshold, threshold2=None, n=n, coordinates='centers')

# plot for verification
for i, (y, x) in enumerate(pos):
    for ax in axes:
        if i == 0: 
            ax.scatter(x, y, color='white', label='Brightest pixel around a spike in the enhanced image')
        else:      
            ax.scatter(x, y, color='white')

# using enhanced image to estimate brightest pixels
# red dots are the brightest pixel above threshold in the normal image around 6 pixels of distance
pos2 = []
for y, x in pos:
    pos2.append(im2.get_brightest_pixel_position(y=y, x=x, n=n, coordinates='centers'))

# plot for verification
for i, (y, x) in enumerate(pos2):
    for ax in axes:
        if i == 0: 
            ax.scatter(x, y, color='red', label='Brightest pixel around a spike in the raw image')
        else:      
            ax.scatter(x, y, color='red')

# get center of mass
pos3 = []
for y, x in pos2:
    spot, _, _, _ = im2.get_spot(y=y, x=x, n=_n2, coordinates='centers')
    pos3.append(spot.get_center_of_mass())

# plot for verification
for i, (y, x) in enumerate(pos3):
    for ax in axes:
        if i == 0: 
            ax.scatter(x, y, color='magenta', label='Center of mass around a spike in the raw image')
        else:      
            ax.scatter(x, y, color='magenta')
        ax.scatter(x, y, s=10, edgecolors='magenta', fc='None', linewidths=10*_n2)

# patch
if len(pos) > 0:
    im4 = im.patch(pos=pos2, n=_patch_size, value=None, coordinates='centers')
    im4.pcolormesh(ax=axes[2])
else:
    im2.pcolormesh(ax=axes[2])

# legend, titles, and labels
axes[1].legend()
axes[0].set_title('Raw (floored)')
axes[1].set_title('Enhanced')
if len(pos) > 0:
    axes[2].set_title('Patched')
else:
    axes[2].set_title('Patched (nothing was patched)')
axes[0].set_xlabel('x pixel')
axes[1].set_xlabel('x pixel')
axes[2].set_xlabel('x pixel')
axes[0].set_ylabel('y pixel')
# %%

# %% Spike finding and patching
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scan    = 24
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)
im = ims[0]

# parameters
n  = 6
_n = 4
threshold = 1e4#3e4  # note that threshold must be interms of enhanced image as enhance=True in im.find_and_patch()
_patch_size = n         

# find and patch
im2, pe = im.find_and_patch(n, threshold, threshold2=None, enhance=True, _n=_n, _bkg=None, _patch_size=_patch_size, _patch_value=None)

# plot for verification
fig, axes = br.subplots(1, 3, sharex=True, sharey=True, figsize=(40, 12), layout='constrained')
im.pcolormesh(ax=axes[0])
im.enhance(n=_n).pcolormesh(vmin=0, vmax=2000, ax=axes[1])
im2.pcolormesh(ax=axes[2])
for ax in axes: 
    pe.plot(ax=ax, s=10, edgecolors='magenta', fc='None', linewidths=5*_patch_size)
    br.labels.detector(ax=ax)
axes[0].set_title('raw')
axes[1].set_title('enhanced')
axes[2].set_title(f'patched ({len(pe)})')
# %%

# %% centroid example 1
sample  = 'align_rixs'
dataset = 'align_0001'
scan    = 125
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)
im = ims[1]

# parameters
cosmic   = dict(n=6, _n=4, threshold=1e4, enhance=True)
centroid = dict(n=1, threshold=400, threshold2=None, floor=True, enhance=True, 
                _n=None, _bkg=None, _n2=None, spot_zeroing_type='zero', 
                MAX_NUMBER_OF_CANDIDATES=1000)

# find and patch
im2, pec = im.find_and_patch(**cosmic)

# centroid
pe, pe2 = im2.floor().centroid(**centroid)

# plot for verification
fig, axes = br.subplots(1, 4, sharex=True, sharey=True, figsize=(40, 12), layout='constrained')
im.plot(ax=axes[0])
im.enhance(n=cosmic['n']).plot(vmin=0, vmax=2000, ax=axes[1])
im2.enhance(n=centroid['n']).plot(vmin=0, vmax=2000, ax=axes[2])
im2.floor().plot(ax=axes[3])
for ax in axes: 
    pec.plot(ax=ax, s=10, edgecolors='magenta', fc='None', linewidths=5*cosmic['n'], label='cosmic-rays')
    pe.plot(ax=ax, s=10, edgecolors='green', fc='None', linewidths=5*centroid['n'], label='photon events')
    pe2.plot(ax=ax, s=10, edgecolors='red', fc='None', linewidths=5*centroid['n'], label='double events')
    br.labels.detector(ax=ax)
axes[0].set_title('raw')
axes[1].set_title('enhanced for cosmic')
axes[2].set_title(f'patched and enhanced for photon events ({len(pe)})')
axes[3].set_title(f'patched and floored (patches={len(pec)})')
# %%

# %% centroid example 2
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scan    = 24
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)
im = ims[0]

# parameters
cosmic   = dict(n=6, _n=4, threshold=1e4, enhance=True)
centroid = dict(n=1, threshold=400, threshold2=900, floor=True, enhance=True, 
                _n=None, _bkg=None, _n2=None, spot_zeroing_type='zero', 
                MAX_NUMBER_OF_CANDIDATES=2000)

# find and patch
im2, pec = im.find_and_patch(**cosmic)

# centroid
pe, pe2 = im2.centroid(**centroid)

# plot for verification
fig, axes = br.subplots(1, 4, sharex=True, sharey=True, figsize=(40, 12), layout='constrained')
im.floor().plot(ax=axes[0])
im.enhance(n=cosmic['n']).plot(vmin=0, vmax=2000, ax=axes[1])
im2.enhance(n=centroid['n']).plot(vmin=0, vmax=2000, ax=axes[2])
im2.floor().plot(ax=axes[3])
for ax in axes: 
    pec.plot(ax=ax, s=10, edgecolors='magenta', fc='None', linewidths=2*cosmic['n'], label='cosmic-rays')
    pe.plot(ax=ax, s=10, edgecolors='green', fc='None', linewidths=5*centroid['n'], label='photon events')
    pe2.plot(ax=ax, s=10, edgecolors='red', fc='None', linewidths=5*centroid['n'], label='double events')
    br.labels.detector(ax=ax)
axes[1].legend()
axes[0].set_title('raw (floored)')
axes[1].set_title('enhanced for cosmic')
axes[2].set_title(f'patched and enhanced for photon events ({len(pe)})')
axes[3].set_title(f'patched and floored (patches={len(pec)})')
# %%

# %% Get spectrum (step-by-step)
sample  = 'align_rixs'
dataset = 'align_0001'
scan    = 125

# centroid parameters
cosmic   = dict(n=6, _n=4, threshold=1e4)
centroid = dict(n=1, threshold=400, threshold2=900, floor=True, enhance=True, 
                _n=None, _bkg=None, _n2=None, spot_zeroing_type='zero', 
                MAX_NUMBER_OF_CANDIDATES=2000)
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = 6000
calib    = 8.5878e-3  # eV/px

# get images
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)

# remove cosmic rays, centroid, and adding photon events
pe  = br.PhotonEvents()
pe2 = br.PhotonEvents()
for _im in ims:
    _im2, _pec = _im.find_and_patch(**cosmic)
    _pe, _pe2  = _im2.floor().centroid(**centroid)
    pe  = pe + _pe
    pe2 = pe2 + _pe2

# curvature correction
pe  = pe.set_vertical_shift_via_polyval(curv)
pe2 = pe2.set_vertical_shift_via_polyval(curv)

# plot for verification
br.figure()
pe.plot(label='photon hits')
pe2.plot(label='double events')
br.labels.detector()

# binning
s  = pe.binning(1, nbins).integrated_rows_vs_y_centers()
s2 = pe2.binning(1, nbins).integrated_rows_vs_y_centers()

# calibrate
s  = s.set_calib(calib).set_shift(-9.48)
s2 = s2.set_calib(calib).set_shift(-9.48)

# plot for verification
br.figure()
s.plot(label='photon hits')
s2.plot(label='double events')
br.labels.rixs()
# %%

# %% comparing processed spectrum with online-processed with id32.process()
sample  = 'align_rixs'
dataset = 'align_0001'
scan    = 125

# get online-processed spectrum for comparisson
s1 = id32.read(TOP, sample, dataset, scan)

# integration mode (no cosmic rays removal)
cosmic   = None
centroid = None
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = None
calib    = 8.5878e-3  # eV/px
include_double_events = True
s2 = id32.process(TOP, sample, dataset, scan, nbins=nbins, curv=curv, calib=calib, cosmic=cosmic, centroid=centroid, include_double_events=include_double_events)
s2 = s2.floor(limits=(None, 5))

# integration mode
cosmic   = dict(n=6, _n=4, threshold=1e4)
centroid = None
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = None
calib    = 8.5878e-3  # eV/px
include_double_events = True
s3 = id32.process(TOP, sample, dataset, scan, nbins=nbins, curv=curv, calib=calib, cosmic=cosmic, centroid=centroid, include_double_events=include_double_events)
s3 = s3.floor(limits=(None, 5))

# centroid 1 (no cosmic rays removal)
cosmic   = None
centroid = dict(n=1, threshold=400, threshold2=None, floor=True, MAX_NUMBER_OF_CANDIDATES=2000)
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = 6000
calib    = 8.5878e-3  # eV/px
include_double_events = True
s4 = id32.process(TOP, sample, dataset, scan, nbins=nbins, curv=curv, calib=calib, cosmic=cosmic, centroid=centroid, include_double_events=include_double_events)

# centroid 1
cosmic   = dict(n=6, _n=4, threshold=1e4)
centroid = dict(n=1, threshold=400, threshold2=None, floor=True, MAX_NUMBER_OF_CANDIDATES=2000)
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = 6000
calib    = 8.5878e-3  # eV/px
include_double_events = True
s5 = id32.process(TOP, sample, dataset, scan, nbins=nbins, curv=curv, calib=calib, cosmic=cosmic, centroid=centroid, include_double_events=include_double_events)

# centroid 2
cosmic   = dict(n=6, _n=4, threshold=1e4)
centroid = dict(n=1, threshold=500, threshold2=None, floor=True, MAX_NUMBER_OF_CANDIDATES=2000)
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = 6000
calib    = 8.5878e-3  # eV/px
include_double_events = True
s6 = id32.process(TOP, sample, dataset, scan, nbins=nbins, curv=curv, calib=calib, cosmic=cosmic, centroid=centroid, include_double_events=include_double_events)

# centroid 3
cosmic   = dict(n=6, _n=4, threshold=1e4)
centroid = dict(n=1, threshold=600, threshold2=None, floor=True, MAX_NUMBER_OF_CANDIDATES=2000)
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = 6000
calib    = 8.5878e-3  # eV/px
include_double_events = True
s7 = id32.process(TOP, sample, dataset, scan, nbins=nbins, curv=curv, calib=calib, cosmic=cosmic, centroid=centroid, include_double_events=include_double_events, return_photon_events=True)

# align and normalize
ss = br.Spectra((s1, s2, s3, s4, s5, s6, s7))
for i, s in enumerate(ss):
    s = s.set_factor(1/max(s))
    smooth, popt, err, f = s.set_shift(-9.48).fit_peak(guess_A=1, guess_c=0, guess_w=0.04, fixed_m=0)
    s = s.set_shift(-9.48-popt[1])
    ss[i] = s

# plot for verification
br.figure()
ss[0].plot(marker='o', label='online-processed')
ss[1].plot(marker='o', label='integration mode (no cosmic rays removal)')
ss[2].plot(label='integration mode')
ss[3].plot(marker='o', label='centroid parameters 1 (no cosmic rays removal)')
ss[4].plot(label='centroid parameters 1')
ss[5].plot(label='centroid parameters 2')
ss[6].plot(label='centroid parameters 3')
br.leg()
br.labels.rixs()
# %%

