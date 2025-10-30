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
# br.settings.FIGURE_POSITION = (80, 200)
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
    imagelist = br.parsed_filelist(TOP/f'data/RIXS/TIF/{str(scan).zfill(4)}', string='.tif', ref=1, return_type='dict')
    if index is None:
        im = br.Dummy()
        for _index in imagelist:
            im.append(br.Image(data=cv2.imread(imagelist[_index], -1)))
    else:
        assert index >= 0, 'index must be a positive integer or zero'
        assert index < number_of_images_inside_scan(scan, TOP=TOP), f'index not valid.'
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
n = 2   # photon hits candidates that are within n pixels of distance from each other will be considered the same candidate.
threshold = 500  # pixel intensity for detecting candidates
_n = 12  # moving average window for enhancing image 

# raw image
im   = read_rixs_image(scan, index=index)
_pe1, _pe2, _, _ = get_photon_events(scan, index=index)
pe = _pe1 + _pe2

# enhance image
im2 = im.enhance(nx=_n, ny=_n)

# simulated dark image (remove photon events and cosmic rays from random image from th same scan)
d = read_rixs_image(scan, index=index+1)
_pe1, _pe2, _, _ = get_photon_events(scan, index=index+1)
d = d.patch([(y, x) for x, y in _pe1 + _pe2], nx=6, ny=6)  # remove photon events
d, _ = d.find_and_patch(nx=4, ny=4, threshold=1e4, _bkg=None, _square=False, _nx=1, _ny=1, _patch_x_size=4, _patch_y_size=4)  # remove cosmic rays
im3 = im - d

# enhance image
im4 = im3.enhance(nx=_n, ny=_n, bkg=None)

# detect cosmic rays from the raw image
pec, _ = im.find_candidates(nx=4, ny=4, threshold=5e3, _bkg=None, _square=False, _nx=1, _ny=1)

# centroid from the raw image (removed cosmic events)
temp = im.patch(pos=[(_[1], _[0]) for _ in pec], nx=_n, ny=_n)
peb, _ = temp.centroid(nx=n, ny=n, threshold=threshold, _bkg=None, _square=False, _nx=1, _ny=1)

# plot
fig, axes = br.subplots(1, 4, sharex=True, sharey=True, figsize=(46, 12), layout='constrained')
im.plot(ax=axes[0],  vmin=0, vmax=800, origin='lower')
im2.plot(ax=axes[1], vmin=0, vmax=1e4, origin='lower')
im3.plot(ax=axes[2], vmin=0, vmax=1000, origin='lower')
im4.plot(ax=axes[3], vmin=0, vmax=1e4, origin='lower')

ax = axes[0]
pe.plot(ax=ax, s=50, edgecolors='red', fc='None', linewidths=1, label=f'photon events from xcam (counts={len(pe)})')
pec.plot(ax=ax, s=50, edgecolors='white', fc='None', linewidths=1, label=f'cosmic rays from brixs (counts={len(pec)})')
peb.plot(ax=ax, s=80, edgecolors='orange', fc='None', linewidths=1, label=f'photon events from raw brixs (counts={len(peb)})')

axes[0].legend(labelcolor='black')
axes[0].set_title('raw')
axes[1].set_title('raw enhanced')
axes[2].set_title(f'raw - dark')
axes[3].set_title(f'(raw - dark) enhanced')
# %%


# %  ===================================================================== %% #
# %  ============== centroid parameters 1 (no moving window) ============= %% #
# %% ===================================================================== %% # 
scan = 100
# scan = 130
# scan = 124

# parameters
cosmic   = dict(nx=6, ny=6, threshold=5e3)
centroid = dict(nx=2, ny=2, threshold=500)  # photon hit is assumed to excite at up to second neighbors
enhance  = dict(nx=20, ny=20)

# read image, remove cosmic rays, and centroid
ims = read_rixs_image(scan)
enh = br.Dummy()
pe1, pe2, pec = br.Dummy(), br.Dummy(), br.Dummy()
for _im in ims:
    _temp, _pec = _im.find_and_patch(**cosmic)
    _pe1, _pe2  = _temp.centroid(**centroid)
    enh.append(_im.enhance(**enhance))
    pe1.append(_pe1)
    pe2.append(_pe2)
    pec.append(_pec)

# get photon events from xcam
pex = br.Dummy()
for i in range(len(ims)):
    _pe_a, _pe_b, _, _ = get_photon_events(scan, index=i)
    pex.append(_pe_a + _pe_b)

# plot
fig, axes = br.subplots(4, 6, sharex=True, sharey=True, hspace=0.04, wspace=0.1, figsize=(46, 16))
plt.subplots_adjust(left=0.04, right=0.98, top=0.96, bottom=0.1)
for ax in axes[12:]:
    ax.ymove(-0.05)
for i in range(len(ims[:12])):
    if i > 5: j = i + 6
    else: j = i
    ims[i].plot(ax=axes[j],   vmin=0, vmax=800, origin='lower')
    enh[i].plot(ax=axes[j+6], vmin=0, vmax=5e3, origin='lower')
    pex[i].plot(ax=axes[j], s=40, edgecolors='red',    fc='None', linewidths=1, label=f'photon events from xcam (counts={len(pex[i])})')
    pec[i].plot(ax=axes[j], s=60, edgecolors='white',  fc='None', linewidths=1, label=f'cosmic rays from brixs (counts={len(pec[i])})')
    pe1[i].plot(ax=axes[j], s=80, edgecolors='orange', fc='None', linewidths=1, label=f'photon events from raw brixs (counts={len(pe1[i])})')
    pe2[i].plot(ax=axes[j], s=100, edgecolors='green',  fc='None', linewidths=1, label=f'double events from raw brixs (counts={len(pe2[i])})')

axes[0].legend(labelcolor='black', fontsize='xx-small')
for i, ax in enumerate(axes.rows[0] + axes.rows[2]):
    ax.set_title(f'scan {scan}, image {i}', fontsize='x-small')
    ax.note(f'{len(pex[i])}, {len(pec[i])}, {len(pe1[i])}, {len(pe2[i])}', loc='lower right', color='white', fontsize='xx-small')
for ax in axes.rows[1] + axes.rows[3]:
    ax.note('Enhanced', loc='lower right', color='white', fontsize='x-small')
_ = br.label_axes(axes)
# %%


# %  ===================================================================== %% #
# %  ================ centroid parameters 2 (moving window) ============== %% #
# %% ===================================================================== %% # 
scan = 100
# scan = 130
# scan = 124

# parameters
cosmic   = dict(nx=6, ny=6, threshold=4e3)
centroid = dict(nx=2, ny=2, threshold=260, _nx=3, _ny=3)
enhance  = dict(nx=20, ny=20)

# read image, remove cosmic rays, and centroid
ims = read_rixs_image(scan)
enh = br.Dummy()
pe1, pe2, pec = br.Dummy(), br.Dummy(), br.Dummy()
for _im in ims:
    _temp, _pec = _im.find_and_patch(**cosmic)
    _pe1, _pe2  = _temp.centroid(**centroid)
    enh.append(_im.enhance(**enhance))
    pe1.append(_pe1)
    pe2.append(_pe2)
    pec.append(_pec)

# get photon events from xcam
pex = br.Dummy()
for i in range(len(ims)):
    _pe_a, _pe_b, _, _ = get_photon_events(scan, index=i)
    pex.append(_pe_a + _pe_b)

# plot
fig, axes = br.subplots(4, 6, sharex=True, sharey=True, hspace=0.04, wspace=0.1, figsize=(46, 16))
plt.subplots_adjust(left=0.04, right=0.98, top=0.96, bottom=0.1)
for ax in axes[12:]:
    ax.ymove(-0.05)
for i in range(len(ims[:12])):
    if i > 5: j = i + 6
    else: j = i
    ims[i].plot(ax=axes[j],   vmin=0, vmax=800, origin='lower')
    enh[i].plot(ax=axes[j+6], vmin=0, vmax=5e3, origin='lower')
    pex[i].plot(ax=axes[j], s=40, edgecolors='red',    fc='None', linewidths=1, label=f'photon events from xcam (counts={len(pex[i])})')
    pec[i].plot(ax=axes[j], s=60, edgecolors='white',  fc='None', linewidths=1, label=f'cosmic rays from brixs (counts={len(pec[i])})')
    pe1[i].plot(ax=axes[j], s=80, edgecolors='orange', fc='None', linewidths=1, label=f'photon events from raw brixs (counts={len(pe1[i])})')
    pe2[i].plot(ax=axes[j], s=100, edgecolors='green',  fc='None', linewidths=1, label=f'double events from raw brixs (counts={len(pe2[i])})')

axes[0].legend(labelcolor='black', fontsize='xx-small')
for i, ax in enumerate(axes.rows[0] + axes.rows[2]):
    ax.set_title(f'scan {scan}, image {i}', fontsize='x-small')
    ax.note(f'{len(pex[i])}, {len(pec[i])}, {len(pe1[i])}, {len(pe2[i])}', loc='lower right', color='white', fontsize='xx-small')
for ax in axes.rows[1] + axes.rows[3]:
    ax.note('Enhanced', loc='lower right', color='white', fontsize='x-small')
_ = br.label_axes(axes)
# %%




