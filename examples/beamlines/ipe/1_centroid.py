#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Example script for running the centroid algorithm on images collected at the
at the IPE beamline of SIRIUS.

Example data for running this example can be downloaded at this Onedrive link
https://1drv.ms/f/c/6666810d266a3cc9/EjNzJssVTWhFqAOx0-nDCxoBSNt2cvvF8zqOmNdx_zzslQ?e=BOAgx4
Note that some less important files have been removed from this example data to
make the folder lighter and easier to download.

This example also requires the cv2 package for reading tif images.
pip install opencv-python

# approximate values
# pixel saturation value ~ 16K300
# pixel value for single photon hit ~700 - 2e3 (bkg ~200)
# enhanced (_n=4)       --> 4e4 (bkg ~4e2)
# raw - dark            --> 600 (bkg ~30)
# (raw - dark) enhanced --> 5e4 (bkg ~ 1e3)
"""

# %% ========================== standard imports ========================= %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# %% =========================== brixs imports =========================== %% #
#### The following two lines are only necessary if you are importing brixs 
#### from a local file
# import sys
# sys.path.insert(0, 'brixs/')
import brixs as br
import brixs.beamlines.ipe as ipe
import brixs.addons.centroid

# %% ============================== settings ============================= %% #
# matplotlib (optional)
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()

# brixs (optional)
br.settings.FIGURE_POSITION = (80, 200)
# br.get_window_position()
# %%

# %% ============================= folderpaths =========================== %% #
# the path to the top/main folder is going to be used constantly, so we are better
# off defining a variable for it
TOP = Path(r'C:\Users\galdin_c\github\brixsexampledata\beamlines\ipe')
# %%

# %  ===================================================================== %% #
# %  ==================== defining supporting functions ================== %% #
# %% ===================================================================== %% #
# as for right now, IPE does not have a dedicated function for reading the image
# files, so we have to quickly define a function for that
def list_rixs_scans_with_images(TOP=TOP):
    """Return list of rixs scans available"""
    TOP = Path(TOP)
    assert TOP.exists(), f'Top directory does not exist TOP="{TOP}"'
    assert 'data' in [_.name for _ in br.filelist(TOP)], f'folder `data` cannot be found in TOP="{TOP}". TOP must be the top/main folderpath of the experiment'
    return list(br.parsed_filelist(TOP/'data/RIXS/TIF', string='*', ref=0, return_type='dict').keys())

def number_of_images_inside_scan(scan, TOP=TOP):
    """return the number of images in a scan"""
    assert scan in list_rixs_scans_with_images(TOP), f'scan={scan} does not have images saved'
    return len(br.filelist(TOP/f'data/RIXS/TIF/{str(scan).zfill(4)}', string='.tif'))

import cv2
def read_rixs_image(scan, index=None, TOP=TOP):
    """Return images from a rixs scan

    Args:
        TOP (str or path): TOP directory is expected to point to
            the top (main) experiment folder where one can find the subfolders:
            'data' and 'proc'.
        scan (int): rixs scan number
        index (int or None, optional): index of the image, e.g., if a scan have
            10 images, one can select which image to import via index. If None,
            all images in a scan will be imported. Default is None.

    Returns:
        list of images if index is None
        single image otherwise
    """
    TOP = Path(TOP)
    assert index >= 0, 'index must be a positive integer or zero'
    assert index < number_of_images_inside_scan(scan, TOP=TOP), f'index not valid.'

    imagelist = br.parsed_filelist(TOP/f'data/RIXS/TIF/{str(scan).zfill(4)}', string='.tif', ref=1, return_type='dict')
    if index is None:
        im = br.Dummy()
        for _index in imagelist:
            im.append(br.Image(data=cv2.imread(imagelist[index], -1)))
    else:
        im = br.Image(data=cv2.imread(imagelist[index], -1))
    return im

def get_spectrum(scan, sbins, index=None, calib=None, norm=True, start=0, stop=None, TOP=TOP):
    """shortcut function to get processed spectrum"""
    assert index >= 0, 'index must be a positive integer or zero'
    assert index < number_of_images_inside_scan(scan, TOP=TOP), f'index not valid.'

    if index is None:
        skip = []
    else:
        skip = list(np.arange(0, number_of_images_inside_scan(scan, TOP=TOP)))
        skip.remove(index)
    return ipe.process(folderpath=Path(TOP)/f'data/RIXS/SPE/{str(scan).zfill(4)}', sbins=sbins, calib=calib, norm=norm, start=start, stop=stop, skip=skip)

def get_photon_events(scan, index=None, verbose=False, start=0, stop=None, curv=False, TOP=TOP):
    """shortcut function to get photon events for a rixs scan"""
    assert index >= 0, 'index must be a positive integer or zero'
    assert index < number_of_images_inside_scan(scan, TOP=TOP), f'index not valid.'

    if index is None:
        skip = []
    else:
        skip = list(np.arange(0, number_of_images_inside_scan(scan, TOP=TOP)))
        skip.remove(index)
    return ipe.read(fpath=Path(TOP)/f'data/RIXS/SPE/{str(scan).zfill(4)}', verbose=verbose, start=start, stop=stop, skip=skip, curv=curv)
# %%


# %  ===================================================================== %% #
# %  ============ enhancing image to visually detect photon hits ========= %% #
# %% ===================================================================== %% #
scan = 100
index = 0
_n = 8

# raw image
im   = read_rixs_image(scan, index=index)
_pe1, _pe2, _, _ = get_photon_events(100, index=index)
pe = _pe1 + _pe2

# enhance image
im2 = im.enhance(_n)

# simulated dark image (remove photon events and cosmic rays from random image)
d = read_rixs_image(scan, index=index+1)
_pe1, _pe2, _, _ = get_photon_events(scan, index=index+1)
d = d.patch([(y, x) for x, y in _pe1 + _pe2], n=6)  # remove photon events
d, _ = d.find_and_patch(n=4, threshold=1e4, _bkg=0, _square=False, _n=1, _patch_size=4)  # remove cosmic rays
im3 = im - d

# enhance image
im4 = im3.enhance(_n, 0)

# detect cosmic rays from the raw image
pec, _ = im.find_candidates(n=4, threshold=5e3, _bkg=0, _square=False, _n=1)

# plot
fig, axes = br.subplots(1, 4, sharex=True, sharey=True, figsize=(46, 12), layout='constrained')
im.pcolormesh(ax=axes[0],  vmin=0, vmax=800)
im2.pcolormesh(ax=axes[1], vmin=0, vmax=1e4)
im3.pcolormesh(ax=axes[2], vmin=0, vmax=1000)
im4.pcolormesh(ax=axes[3], vmin=0, vmax=1e4)
for ax in axes: 
    pe.plot(ax=ax, s=50, edgecolors='red', fc='None', linewidths=1, label='photon events (from xcam centroid)')
for ax in axes: 
    pec.plot(ax=ax, s=70, edgecolors='white', fc='None', linewidths=1, label='cosmic rays (from brixs)')

axes[1].legend(labelcolor='black')
axes[0].set_title('raw')
axes[1].set_title('raw enhanced')
axes[2].set_title(f'raw - dark')
axes[3].set_title(f'(raw - dark) enhanced')
# %%

quit()
ipython
# %%
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import brixs as br
import brixs.beamlines.ipe as ipe
import brixs.addons.centroid
from brixs.addons.sound import make_sound
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()
br.settings.FIGURE_POSITION = (80, 100)
TOP = Path(r'C:\Users\galdin_c\github\brixsexampledata\beamlines\ipe')
# %%

# %  ===================================================================== %% #
# %  ============ enhancing image to visually detect photon hits ========= %% #
# %% ===================================================================== %% # 
scan = 100
index = 0
n = 2   # photon hits candidates that are within n pixels of distance from each other will be considered the same candidate.
threshold = 500  # pixel intensity for detecting candidates
_n = 8  # moving average window for enhancing image 

# raw image
im   = read_rixs_image(scan, index=index)
_pe1, _pe2, _, _ = get_photon_events(100, index=index)
pe = _pe1 + _pe2

# enhance image
im2 = im.enhance(_n)

# simulated dark image (remove photon events and cosmic rays from random image from th same scan)
d = read_rixs_image(scan, index=index+1)
_pe1, _pe2, _, _ = get_photon_events(scan, index=index+1)
d = d.patch([(y, x) for x, y in _pe1 + _pe2], n=6)  # remove photon events
d, _ = d.find_and_patch(n=4, threshold=1e4, _bkg=0, _square=False, _n=1, _patch_size=4)  # remove cosmic rays
im3 = im - d

# enhance image
im4 = im3.enhance(_n, 0)

# detect cosmic rays from the raw image
pec, _ = im.find_candidates(n=4, threshold=5e3, _bkg=0, _square=False, _n=1)

# centroid from the raw image (removed cosmic events)
temp = im.patch(pos=[(_[1], _[0]) for _ in pec], n=_n)
peb, _ = temp.centroid(n=n, threshold=threshold, _bkg=0, _square=False, _n=1)

# plot
fig, axes = br.subplots(1, 4, sharex=True, sharey=True, figsize=(46, 12), layout='constrained')
im.pcolormesh(ax=axes[0],  vmin=0, vmax=800)
im2.pcolormesh(ax=axes[1], vmin=0, vmax=1e4)
im3.pcolormesh(ax=axes[2], vmin=0, vmax=1000)
im4.pcolormesh(ax=axes[3], vmin=0, vmax=1e4)
for ax in axes: 
    pe.plot(ax=ax, s=50, edgecolors='red', fc='None', linewidths=1, label=f'photon events from xcam (counts={len(pe)})')
for ax in axes: 
    pec.plot(ax=ax, s=50, edgecolors='white', fc='None', linewidths=1, label=f'cosmic rays from brixs (counts={len(pec)})')
for ax in axes: 
    peb.plot(ax=ax, s=80, edgecolors='orange', fc='None', linewidths=1, label=f'photon events from raw brixs (counts={len(peb)})')

axes[1].legend(labelcolor='black')
axes[0].set_title('raw')
axes[1].set_title('raw enhanced')
axes[2].set_title(f'raw - dark')
axes[3].set_title(f'(raw - dark) enhanced')
# %%


# %  ===================================================================== %% #
# %  ======================== optimizing centroid ======================== %% #
# %% ===================================================================== %% # 
scan = 100
index = 0

# photon hist are assumed to excite at least the first neighbors. Enhancing the
# image does a moving average which makes sure that a pixel will only became a
# candidate if pixels around it were also light up


peb, _ = temp.centroid(n=n, threshold=500, enhance=True)



n = 2   # photon hits candidates that are within n pixels of distance from each other will be considered the same candidate.
_n = 8  # moving average window for enhancing image, 

# raw image
im   = read_rixs_image(scan, index=index)
_pe1, _pe2, _, _ = get_photon_events(100, index=index)
pe = _pe1 + _pe2

# enhance image
im2 = im.enhance(_n)

# simulated dark image (remove photon events and cosmic rays from random image from th same scan)
d = read_rixs_image(scan, index=index+1)
_pe1, _pe2, _, _ = get_photon_events(scan, index=index+1)
d = d.patch([(y, x) for x, y in _pe1 + _pe2], n=6)  # remove photon events
d, _ = d.find_and_patch(n=6, threshold=10000)       # remove cosmic rays
im3 = im - d

# enhance image
im4 = im3.enhance(_n, 0)

# detect cosmic rays from the raw image
pec, _ = im.find_candidates(n=6, threshold=5e3, enhance=False)

# centroid from the raw image (removed cosmic events)
temp = im.patch(pos=[(_[1], _[0]) for _ in pec], n=_n)
peb, _ = temp.centroid(n=n, threshold=500, enhance=False)

# plot
fig, axes = br.subplots(1, 4, sharex=True, sharey=True, figsize=(46, 12), layout='constrained')
im.pcolormesh(ax=axes[0],  vmin=0, vmax=800)
im2.pcolormesh(ax=axes[1], vmin=0, vmax=1e4)
im3.pcolormesh(ax=axes[2], vmin=0, vmax=1000)
im4.pcolormesh(ax=axes[3], vmin=0, vmax=1e4)
for ax in axes: 
    pe.plot(ax=ax, s=50, edgecolors='red', fc='None', linewidths=1, label=f'photon events from xcam (counts={len(pe)})')
for ax in axes: 
    pec.plot(ax=ax, s=50, edgecolors='white', fc='None', linewidths=1, label=f'cosmic rays from brixs (counts={len(pec)})')
for ax in axes: 
    peb.plot(ax=ax, s=80, edgecolors='orange', fc='None', linewidths=1, label=f'photon events from raw brixs (counts={len(peb)})')

axes[1].legend(labelcolor='black')
axes[0].set_title('raw')
axes[1].set_title('raw enhanced')
axes[2].set_title(f'raw - dark')
axes[3].set_title(f'(raw - dark) enhanced')
# %%







