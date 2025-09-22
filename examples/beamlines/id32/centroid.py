#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Example script for running the centroid algorithm on images collected at the
at the ID32 beamline of ESRF.

Example data for running this example can be downloaded at this Onedrive link
https://1drv.ms/f/c/6666810d266a3cc9/EjNzJssVTWhFqAOx0-nDCxoBSNt2cvvF8zqOmNdx_zzslQ?e=BOAgx4
Note that some less important files have been removed from this example data to
make the folder lighter and easier to download.

the images from the detector are large. Therefore, the example data only
has images for a few selected spectra
for the rixs branch, the following scans have images,
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scans   = np.arange(19, 30+1)

sample  = 'align_rixs'
dataset = 'align_0001'
scans   = 125, 126

and for the xmcd branch, the follwing scan have images,
sample  = 'sample2_xmcd'
dataset = 'align_0002'
scan    = 18, 19
"""


# %% ========================== standard imports ========================= %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# %% =========================== brixs imports =========================== %% #
#### The following two lines are only necessary if you are importing brixs 
#### from a local file
# import sys
# sys.path.insert(0, 'brixs-main/')
import brixs as br
import brixs.beamlines.id32 as id32
import brixs.addons.centroid
import brixs.addons.fitting

# %% ============================== settings ============================= %% #
# matplotlib (optional)
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()

# brixs (optional)
# br.settings.FIGURE_POSITION = (283, 567)
# br.get_window_position()
# %%

# %% ============================= folderpaths =========================== %% #
# the path to the top/main folder is going to be used constantly, so we are better
# off defining a variable for it
TOP = Path(r'C:\Users\galdin_c\github\brixsexampledata\beamlines\id32')
# %%


# %  ===================================================================== %% #
# %  ===================== raw rixs data (CENTROID) ====================== %% #
# %% ===================================================================== %% #

# %% Spike finding and patching (step-by-step)
sample  = 'align_rixs'
dataset = 'align_0001'
scan    = 125
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)
im = ims[1]

# defining 
n   = 6
_n  = 4
_cm_n = n
threshold = 1e4
_patch_size = _cm_n
# n: photon hits candidates that are within n pixels of distance 
# from each other will be considered the same candidate. For better 
# results, set n to be roughly the expected pixel distance between a 
# photon hit and the farthest excited pixel.

# _n (int): Use this to overwrite the size of the moving average window 
# used in the process of enhancing the image, i.e., number of points
# to average. If this is None, the averaging window size will be set 
# to n*2+1. Note that _n=1 implies no moving average.
# For visually inspecting image, sometimes it is 
# best to make _n large to make it easy to see spikes by eye. For finding spikes 
# optimally with im.find_candidates(), one might want to set up _n roughly the
# same size of the spike, e.g., if a spike lights up 4x4 grid of pixels, then _n=4
# should be optimal. If _n is too high, the image is washed out. If _n is too 
# little, the bkg does not go to zero so fast and spikes do not stand out as much.

# _cm_n: Use this to overwrite the number of neighbors when
# calculating the center of mass of a photon hit candidate, e.g., if 
# _cm_n=1, only first neighbors. _cm_n also defines how close two 
# candidates need to be to be considered a double event. If None, 
# _cm_n will be same as n. Default is None. 

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
im2.plot(ax=axes[0])
im3.plot(ax=axes[1])

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
    spot, _, _, _ = im2.get_spot(y=y, x=x, n=_cm_n, coordinates='centers')
    pos3.append(spot.get_center_of_mass())

# plot for verification
for i, (y, x) in enumerate(pos3):
    for ax in axes:
        if i == 0: 
            ax.scatter(x, y, color='magenta', label='Center of mass around a spike in the raw image')
        else:      
            ax.scatter(x, y, color='magenta')
        ax.scatter(x, y, s=10, edgecolors='magenta', fc='None', linewidths=10*_cm_n)

# patch
if len(pos) > 0:
    im4 = im.patch(pos=pos2, n=_patch_size, value=None, coordinates='centers')
    im4.plot(ax=axes[2])
else:
    im2.plot(ax=axes[2])

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
threshold = 1e4#3e4  # note that threshold must be in terms of enhanced image as _square=True in im.find_and_patch()
_patch_size = n         

# find and patch
im2, pe = im.find_and_patch(n, threshold, threshold2=None, _bkg='auto', _square=True, _n=_n, _patch_size=_patch_size, _patch_value=None)

# plot for verification
fig, axes = br.subplots(1, 3, sharex=True, sharey=True, figsize=(40, 12), layout='constrained')
im.plot(ax=axes[0])
im.enhance(n=_n).plot(vmin=0, vmax=2000, ax=axes[1])
im2.plot(ax=axes[2])
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
cosmic   = dict(n=6, threshold=1e4, _bkg='auto', _square=True, _n=4)
centroid = dict(n=1, threshold=400, threshold2=None, _bkg='auto', _square=True, 
                _n=None, _cm_bkg='auto', _cm_n=None, _cm_spot_zero_type='zero', 
                MAX_NUMBER_OF_CANDIDATES=2000)

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
cosmic   = dict(n=6, threshold=1e4, _bkg='auto', _square=True, _n=4)
centroid = dict(n=1, threshold=400, threshold2=900, _bkg='auto', _square=True, 
                _n=None, _cm_bkg='auto', _cm_n=None, _cm_spot_zero_type='zero', 
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
cosmic   = dict(n=6, threshold=1e4, _n=4)
centroid = dict(n=1, threshold=400, threshold2=900, _bkg='auto', _square=True, 
                _n=None, _cm_bkg='auto', _cm_n=None, _cm_spot_zero_type='zero', 
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

# get online-processed spectrum for comparison
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
cosmic   = dict(n=6, threshold=1e4, _n=4)
centroid = None
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = None
calib    = 8.5878e-3  # eV/px
include_double_events = True
s3 = id32.process(TOP, sample, dataset, scan, nbins=nbins, curv=curv, calib=calib, cosmic=cosmic, centroid=centroid, include_double_events=include_double_events)
s3 = s3.floor(limits=(None, 5))

# centroid 1 (no cosmic rays removal)
cosmic   = None
centroid = dict(n=1, threshold=400, threshold2=None, _bkg='auto', _square=True, _n=None, MAX_NUMBER_OF_CANDIDATES=2000)
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = 6000
calib    = 8.5878e-3  # eV/px
include_double_events = True
s4 = id32.process(TOP, sample, dataset, scan, nbins=nbins, curv=curv, calib=calib, cosmic=cosmic, centroid=centroid, include_double_events=include_double_events)

# centroid 1
cosmic   = dict(n=6, threshold=1e4, _n=4)
centroid = dict(n=1, threshold=500, threshold2=None, _bkg='auto',  _square=True, _n=None, MAX_NUMBER_OF_CANDIDATES=2000)
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = 6000
calib    = 8.5878e-3  # eV/px
include_double_events = True
s5 = id32.process(TOP, sample, dataset, scan, nbins=nbins, curv=curv, calib=calib, cosmic=cosmic, centroid=centroid, include_double_events=include_double_events)

# centroid 2
cosmic   = dict(n=6, _n=4, threshold=1e4)
centroid = dict(n=1, threshold=600, threshold2=None, _bkg='auto',  _square=True, _n=None, MAX_NUMBER_OF_CANDIDATES=2000)
curv     = [-1.376e-06, 7.1524e-02, 0]
nbins    = 6000
calib    = 8.5878e-3  # eV/px
include_double_events = True
s6 = id32.process(TOP, sample, dataset, scan, nbins=nbins, curv=curv, calib=calib, cosmic=cosmic, centroid=centroid, include_double_events=include_double_events)

# centroid 3
cosmic   = dict(n=6, _n=4, threshold=1e4)
centroid = dict(n=1, threshold=700, threshold2=None, _bkg='auto',  _square=True, _n=None, MAX_NUMBER_OF_CANDIDATES=2000)
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