#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""The brixs.beamlines.id32 module offers a python solution for data processing
and data analysis for data collected at the ID32 beamline of ESRF.

Example data for running this example can be downloaded at this Onedrive link
https://1drv.ms/f/c/6666810d266a3cc9/EjNzJssVTWhFqAOx0-nDCxoBSNt2cvvF8zqOmNdx_zzslQ?e=BOAgx4
Note that some less important files have been removed from this example data to
make the folder lighter and easier to download.


Besides brixs requirements, brixs.beamlines.id32 module also requires h5py and silx.

Onsite, the path of the top (main) directory is 

/data/visitor/<proposalnumber>/id32/<fulldate>/

for example

/data/visitor/hc1234/id32/20250408/

The folder structure of the top directory is:

GALLERY
jupyter
NOBACKUP
NOTES
PROCESSED_DATA
RAW_DATA
SCRIPTS

Data is stored in RAW_DATA and PROCESSED_DATA. In short, RAW_DATA will contain
rixs detector images, xas, linescans, and mesh scans. PROCESSED_DATA will 
contain the processed rixs spectra.

Inside these two folders, there will be a folder for each defined `sample name`
and for each `sample` there will be multiple `datasets` and each dataset will
have multiple scans.

The structure will be then,

RAW_DATA
    sample1
        dataset1
            scan0038    --> These folders exist only for RIXS scans and they store the detector images
            scan0041
            ...
            dataset1.h5  --> This file contains all data that is not RIXS data
        dataset2
            scan0001
            scan0002
            ...
            dataset2.h5
    sample2
        dataset3
            ...
        dataset4
            ...

PROCESSED_DATA
    sample1
        dataset1
            dhyana95v2    --> folder with detector information for each scan
            workflows     --> folder with data processing for each scan
            Online_analysis_dataset1_dhyana95v2.spec       --> File containing processed RIXS data
            Online_analysis_dataset1_dhyana95v2_info.spec  --> File containing information on rixs scans
        dataset2
            andor1
            workflows
            ...
            Online_analysis_dataset2_andor1.spec
            Online_analysis_dataset2_andor1_info.spec
    sample2
        dataset3
            ...
        dataset4
            ...            

For quick analysis, the important files are:
    RAW_DATA/<sample>/<dataset>/dataset.h5  -> xas, linescans, mesh
    PROCESSED_DATA/<sample>/<dataset>/Online_analysis_dataset2_andor1.spec  -> rixs

This module know where to look for these files. One only has to indicate the 
folderpath of the TOP (main) directory.

One can use this code directly on slurm (jupyter at ESRF servers) through the
following link
https://jupyter-slurm.esrf.fr/

As for today, I do not think one can install brixs via pip directly on your 
jupyter slurm workspace. Therefore, you have to upload a zipped version of brixs
to your folder.

For that, download brixs directly from github

https://github.com/cwgaldino/brixs

click on Code/Download ZIP. The downloaded folder will be named brixs-main.zip.
Click and drag this zipped folder to jupyter slurm directory
To unzip it, open a new tap on jupyter and click on Terminal.
In the terminal, type  

unzip brixs-main.zip

The unzipped folder named `brixs-main` will appear. That's all. 

In you script, make sure to add this folder to your path using the two lines 
below
>>> import sys
>>> sys.path.insert(0, 'brixs-main/')

Note that brixs requires only numpy and matplotlib to work. The brixs.beamlines.id32
modules requires h5py and silx. All packages of which should be installed by 
default at ESRF jupyter slurm.

See examples below.
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
# %  ======================== supporting functions ======================= %% #
# %% ===================================================================== %% #

# help can be used for any function
help(id32.list_samples)

# get list of samples
print(id32.list_samples(TOP))

# get list of datasets
sample  = 'sample1_rixs'
print(id32.list_datasets(TOP, sample))

# get list of scans
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
print(id32.list_scans(TOP, sample, dataset))
print(id32.list_scans(filepath=TOP/'RAW_DATA'/sample/dataset/f'{dataset}.h5'))

# available detector for a scan
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scan = 18
print(id32.list_available_detectors(TOP, sample, dataset, scan))
print(id32.list_available_detectors(filepath=TOP/'RAW_DATA'/sample/dataset/f'{dataset}.h5', scan=scan))

# available detector for a metadata
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scan = 18
print(id32.list_available_metadata(TOP, sample, dataset, scan))
print(id32.list_available_metadata(filepath=TOP/'RAW_DATA'/sample/dataset/f'{dataset}.h5', scan=scan))
# %%


# %  ===================================================================== %% #
# %  =============================== basic =============================== %% #
# %% ===================================================================== %% #

# %% xas (rixs branch) - FULL EXPLANATION
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scan    = 18
ss = id32.read(TOP, sample, dataset, scan)

# ss is a list. The list follows the same order as in the arguments 
# detectors_rixs_branch or detectors_xmcd_branch
# if branch=None, the branch will be selected automatically

# plot for verification
fig, axes = br.subplots(1, 3, layout='constrained')
ss[0].plot(ax=axes[0])
ss[1].plot(ax=axes[1])
ss[2].plot(ax=axes[2])
br.labels.xas('TFY', ax=axes[0])
br.labels.xas('TEY', ax=axes[1])
br.labels.xas('I0',  ax=axes[2])

# this sample/dataset was measured at the rixs branch
# therefore, we know the name of the detectors
# by default, id32.read() will data from the following detectors
# detectors_rixs_branch=['dbig_n', 'sam_n', 'mir_rixs']
# detectors_xmcd_branch=['ifluo_n', 'it_n', 'i0_n_xmcd']
# use the following function to get a list of available detectors 
print(id32.list_available_detectors(TOP, sample, dataset, scan))
# detector are pulled out of the hdf5 address: <scan>/measurement

# for some detectors, there will be three variants:
# 'dbig_raw'  : raw signal
# 'dbig_rixs' : sinal subtracted from a "dark current signal". This is the same
# as _raw, but with an offset
# 'dbig_n'    : signal normalized by I0 -> dbig_rixs/mir_rixs

# for the rixs branch, here are the description of important detectors
# dbig: diode for measuring fluorescence (TFY)
# sam_n: TEY
# mir_rixs: recomposition current of the focusing mirror (I0)

# for the xmcd branch, here are the description of important detectors
# ifluo_n: diode for measuring fluorescence (TFY)
# it_n: TEY
# i0_n_xmcd: recomposition current of the focusing mirror (I0)

# metadata is saved inside each scan
s = ss[0]
s.get_attrs()
print(s.detector)
print(s.moving_motor)
print(s.pol)
print(s.metadata)

# use the following function to get a list of available metadata for a scan
print(id32.list_available_metadata(TOP, sample, dataset, scan))
# metadata are pulled out of the hdf5 address: <scan>/instrument/positioners

# list of all metadata being pulled out is stored in the `metadata` variable
print(id32.metadata['rixs'])  # rixs branch
print(id32.metadata['xmcd'])  # xmcd branch
# id32.metadata[<banch>][<name>] = <name of entry in the file (must match exactly)>
# where name can be anything (something easy to remember)

# new metadata can be added
# simply add item to the `metadata` variable
id32.metadata['rixs']['arm_length'] = 'r1'
# loading data again
ss = id32.read(TOP, sample, dataset, scan)
print(ss[0].metadata['arm_length'])  # new metadata is included

# one can pull detector signal as metadata by using the variable `additional`
print(id32.additional['rixs'])  # rixs branch
print(id32.additional['xmcd'])  # xmcd branch
# These will pulled out from the same place as the "detector": <scan>/measurement
# for XAS data (trigscan's), the `additional` metadata also search for entries
# in the <scan>.2 entry of the scan (not only <scan>.1) as trigscan also save
# some extra metadata in the second entry.

# %%

# %% linescan (rixs branch)
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scan    = 284
ss = id32.read(TOP, sample, dataset, scan)

# plot for verification
fig, axes = br.subplots(1, 3, layout='constrained')
ss[0].plot(ax=axes[0])
ss[1].plot(ax=axes[1])
ss[2].plot(ax=axes[2])
for i, ax in enumerate(axes):
    ax.set_xlabel(ss[i].moving_motor)
    ax.set_ylabel('Intensity (arb. units)')
    ax.grid()
axes[0].set_title('TFY')
axes[1].set_title('TEY')
axes[2].set_title('I0')
# %%

# %% mesh scan (rixs branch)
# loading data
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scan    = 231
ims = id32.read(TOP, sample, dataset, scan)

# plot
br.figure()
ims[1].plot()

# depending on which motor was scanned first, the image can be rotated
print(ims[0].moving_motor)
tfy = ims[0].rotate('counterclockwise').set_factor(1e3)
tey = ims[1].rotate('counterclockwise').set_factor(1e1)

# image can also be flipped
# tfy = ims[0].rotate('counterclockwise').flipx()
# tey = ims[1].rotate('counterclockwise').flipx()

# plot for verification
fig, axes = br.subplots(1, 2, wspace=0.4)
tfy.plot(ax=axes[0], colorbar=True)
tey.plot(ax=axes[1], colorbar=True)
for ax in axes:
    ax.set_xlabel('sample_x')
    ax.set_ylabel('sample_z')
axes[0].set_title('TFY')
axes[1].set_title('TEY')

# %%

# %% rixs (rixs branch)
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scan    = 234
s = id32.read(TOP, sample, dataset, scan)

# plot
br.figure()
s.plot()
br.labels.rixs()

# metadata
s.get_attrs()
s.metadata
# %%

# %% xas (xmcd branch)
sample  = 'sample1_xmcd'
dataset = 'sample1_xmcd_0001'
scan    = 8
ss = id32.read(TOP, sample, dataset, scan)

br.figure()
ss[1].plot()
plt.xlabel(ss[1].moving_motor)
plt.ylabel(ss[1].detector)
plt.grid()

ss[1].get_attrs()
ss[1].metadata
# %%

# %% linescan (xmcd branch)
sample  = 'sample1_xmcd'
dataset = 'sample1_xmcd_0001'
scan    = 4
ss = id32.read(TOP, sample, dataset, scan)

br.figure()
ss[1].plot()
plt.xlabel(ss[1].moving_motor)
plt.ylabel(ss[1].detector)
plt.grid()

ss[1].get_attrs()
ss[1].metadata
# %%

# %% rixs (rixs xmcd)
sample  = 'sample1_xmcd'
dataset = 'sample1_xmcd_0001'
scan    = 17
s = id32.read(TOP, sample, dataset, scan)

# plot
br.figure()
s.set_shift(-5.5).plot()
br.labels.rixs()

# metadata
s.get_attrs()
s.metadata
# %%


# %  ===================================================================== %% #
# %  =========================== rixs advanced =========================== %% #
# %% ===================================================================== %% #

# %% get data from multiple datasets
sample  = 'sample1_xmcd'
energies = [574.5, 574.6, 574.7]
datasets = [f'sample1_GI70_4Tcoil_T20K_{_}0_LH' for _ in energies]

ss = br.Spectra()
for i, dataset in enumerate(datasets):
    _s = id32.read(TOP, sample, dataset, scan=1)
    ss.append(_s)
ss = ss.interp().align().set_shift(-10)

# plot
br.figure()
ss.plot(label=energies)
br.leg()
br.labels.rixs()
# %%

# %% adding rixs spectra
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scans = [19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30]

ss = br.Spectra()
for i, scan in enumerate(scans):
    _s = id32.read(TOP, sample, dataset, scan)
    motor_x = round(_s.metadata['motor_x'], 4)
    _s.label = f'{_s.scan}, x={motor_x}'
    ss.append(_s)

# plot
colors = br.get_colors_from_colormap('rainbow', len(ss))
br.figure()
ss.plot(color=colors)
br.leg()
br.labels.rixs()

# align and calculate sum (or average)
ss2 = ss.interp().align().set_shift(-3.34)
s2  = ss2.calculate_average()
s2.label = 'Average'

# plot
br.figure()
ss2.plot(color=colors)
s2.plot(color='black')
br.leg()
br.labels.rixs()
br.zoom(-1, 6)
# %%

# %% energy map (simulated)
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scans = [19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30]

# let's pretend this scan sequence is a energy map
ss = br.Spectra()
# energies = []  # uncomment this for real energy map     
for i, scan in enumerate(scans):
    _s = id32.read(TOP, sample, dataset, scan)
    # energies.append(round(_s.metadata['E'], 2))  # uncomment this for real energy map
    ss.append(_s)
energies = np.linspace(575, 576, len(scans))  # comment this out for real energy map

# align spectra (sometimes, it must be done mannualy)
# ss2 = ss[i].set_shift(value)  # manually 
ss2 = ss.interp().align().set_shift(-3.34)

# create image
im  = ss2.stack_spectra_as_columns()
im.x_centers = energies

# plot
br.figure()
im.plot()
br.labels.energy_map()
plt.ylim(-1, 6)
# %%

# %% Energy calibration
sample  = 'align_rixs'
dataset = 'align_0001'
scans = [120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132]

# get spectra
ss = br.Spectra()
for i, scan in enumerate(scans):
    _s = id32.read(TOP, sample, dataset, scan)
    _s.x = _s.pixel
    ss.append(_s)

# get energies
energies = [_.metadata['E'] for _ in ss]

# plot for verification
colors = br.get_colors_from_colormap('rainbow', len(ss))
br.figure()
ss.plot(color=colors, label=[str(round(_, 1)) + ' eV' for _ in energies])
br.leg()
plt.xlabel('y (pixel)')
plt.ylabel('Photon count')

# calculate calibration via cross-correlation
out = ss.calculate_calib(values=energies, mode='cc')
calib = -out.popt[0]  # eV/pixel

# plot for verification
br.figure()
out.plot(marker='o', lw=0, color='black', label='data')
out.fit.plot(color='red', label='fit')
br.leg()
plt.xlabel('Shift necessary to align spectra (pixel)')
plt.ylabel('Photon energy (eV)')
# %%


# %  ===================================================================== %% #
# %  =========================== raw rixs data =========================== %% #
# %% ===================================================================== %% #

# the images from the detector are large. Therefore, the example data only
# has images for a few selected spectra
# for the rixs branch, the following scans have images,
# sample  = 'sample1_rixs'
# dataset = 'sample1_rixs_0001'
# scans   = np.arange(19, 30+1)

# sample  = 'align_rixs'
# dataset = 'align_0001'
# scans   = 125, 126

# and for the xmcd branch, the follwing scan have images,
# sample  = 'sample2_xmcd'
# dataset = 'align_0002'
# scan    = 18, 19

# %% rixs (rixs branch) - integration mode
sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scan    = 19
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)

# only one image in this scan
print(ims[0].command)
print(len(ims))

# metadata
ims[0].get_attrs()
ims[0].metadata

# plot
br.figure()
ims[0].plot()

# get curvature poly values from file
_s = id32.read(TOP, sample, dataset, scan, processed_rixs=True)
print(_s.header)
poly = [-9e-07, -0.06906, 0]  # [smile, slope]

# plot curvature
br.figure()
ims[0].plot()
x = np.linspace(0, 2000, 2000)
plt.plot(x, np.polyval(poly, x)+600, color='white')

# fix curvature
im2 = ims[0].set_vertical_shift_via_polyval(-np.array(poly))

# plot curvature corrected image
br.figure()
im2.plot()

# integrate image
s = im2.integrated_rows_vs_y_centers()

# get calibration from file
_s = id32.read(TOP, sample, dataset, scan, processed_rixs=True)
print(_s.header) 
calib = 8.5878e-3  # eV/px

# calibrate spectrum and shit zero (applying factor to make it match with processed data)
s = s.set_calib(calib).set_shift(-4.74).floor(limits=(-2, -1)).set_factor(0.5e-2)

# plot spectrum (integration mode)
br.figure()
s.plot(label='integration')
id32.read(TOP, sample, dataset, scan).set_shift(-3.35).plot(label='online processing (centroid)')
br.leg()
br.labels.rixs()
# %%

# %% rixs (rixs branch) - (sum images) integration mode
sample  = 'align_rixs'
dataset = 'align_0001'
scan    = 125
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)

# ten images
print(ims[0].command)
print(len(ims))

# plot
br.figure()
ims[0].plot()

# sum images
im = ims[0]
for _im in ims[1:]:
    im += _im

# sum images (alternative)
im = ims.calculate_sum()

# plot
br.figure()
im.plot()

# get curvature poly values from file
_s = id32.read(TOP, sample, dataset, scan, processed_rixs=True)
print(_s.header)
poly = [-9e-07, -0.06906, 0]  # [smile, slope]

# plot curvature
br.figure()
im.plot()
x = np.linspace(0, 2000, 2000)
plt.plot(x, np.polyval(poly, x)+1100, color='white')

# fix curvature
im2 = im.set_vertical_shift_via_polyval(-np.array(poly))

# plot curvature corrected image
br.figure()
im2.plot()

# quickly recalculating curvature (see other examples)
_im = im.crop(301, 1850, 800, 1500-1).floor(None, None, 1200, 1500).binning(10, 350)
s = _im.calculate_vertical_shift_curvature()
im2 = im.set_vertical_shift_via_polyval(np.array(s.popt))
br.figure()
im.plot()
s.fit.set_factor(-1).set_offset(1100).plot(color='red')

br.figure()
im2.plot()

# integrate image
s = im2.integrated_rows_vs_y_centers()

# get calibration from file
_s = id32.read(TOP, sample, dataset, scan, processed_rixs=True)
print(_s.header) 
calib = 8.5878e-3  # eV/px

# calibrate spectrum and shit zero (applying factor to make it match with processed data)
s = s.set_calib(calib).set_shift(-9.24).floor(limits=(-2, -1)).set_factor(0.38e-2)

# plot spectrum (integration mode)
br.figure()
s.plot(label='integration')
id32.read(TOP, sample, dataset, scan).set_shift(-0.004).plot(label='online processing (centroid)')
br.leg()
br.labels.rixs()
# %%

# %% rixs (xmcd branch) - curvature from file
sample  = 'sample2_xmcd'
dataset = 'align_0002'
scan    = 18
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)

# ten images in this scan
print(ims[0].command)
print(len(ims))

# plot one image
br.figure()
ims[0].plot()

# sum images
im = ims[0]
for _im in ims[1:]:
    im += _im

# sum images (alternative)
im = ims.calculate_sum()

# plot summed image
br.figure()
im.plot()

# get curvature poly values from file
_s = id32.read(TOP, sample, dataset, scan, processed_rixs=True)
print(_s.header)
poly = [-1.4e-06, -0.0262, 0]  # [smile, slope]

# plot curvature (does not look so good)
br.figure()
im.plot()
x = np.linspace(0, 1790, 2000)
plt.plot(x, np.polyval(poly, x)+1082, color='white')
# %%

# %% rixs (xmcd branch) - recalculating the curvature (step-by-step)
sample  = 'sample2_xmcd'
dataset = 'align_0002'
scan    = 18
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)

# get sum of all images
im = ims.calculate_sum()

# bin image
_im2 = im.binning(14, None)
br.figure()
_im2.plot()

# get columns
_ss2 = _im2.get_columns().floor(limits=(1100, 1500)).crop(800, 1500)
br.figure()
_ss2.plot()

# calculate shift necessary to align spectra
shifts = _ss2.calculate_shift(mode='cc', limits=(800, 1500))
poly = np.polyfit(_im2.x_centers, shifts, deg=2)
x = np.linspace(0, 1790, 2000)

br.figure()
plt.plot(_im2.x_centers, shifts, marker='o', color='black')
plt.plot(x, np.polyval(poly, x), color='red')
plt.xlabel('x (pixel)')
plt.ylabel('y shift (pixels)')

# apply shifts for verification (spectra must be aligned)
br.figure()
_ss2.set_shift(shifts).plot()
plt.xlabel('y (pixel)')

# plot curvature for verification
br.figure()
im.plot()
plt.plot(x, -np.polyval(poly, x)+1082, color='red')
# %%

# %% rixs (xmcd branch) - recalculating the curvature (quicker)
sample  = 'sample2_xmcd'
dataset = 'align_0002'
scan    = 18
ims = id32.read(TOP, sample, dataset, scan, processed_rixs=False)
im = ims.calculate_sum().crop(None, None, 800, 1500).floor(None, None, 1100, 1500)

s = im.binning(14, None).calculate_vertical_shift_curvature()
im2 = im.set_vertical_shift_via_polyval(np.array(s.popt))
x = np.linspace(0, 1790, 2000)

# plot for verification
br.figure()
im.plot()
plt.plot(x, -np.polyval(s.popt, x)+1082, color='red')

br.figure()
im2.plot()
# %%

# %% rixs (xmcd branch) - integration mode
sample  = 'sample2_xmcd'
dataset = 'align_0002'
scan    = 18
im = id32.read(TOP, sample, dataset, scan, processed_rixs=False).calculate_sum()

# fix curvature
curv = [1.760e-06, 2.2670e-02, 0]  # from previous example
im = im.set_vertical_shift_via_polyval(curv)

# integrate image
s = im.integrated_rows_vs_y_centers()

# get calibration from file
_s = id32.read(TOP, sample, dataset, scan, processed_rixs=True)
print(_s.header) 
calib = 79e-3  # eV/px

# calibrate spectrum and shit zero
s = s.set_calib(calib).set_shift(-88.75).floor(limits=(-2, -1)).set_factor(0.35e-2)

# plot spectrum (integration mode)
br.figure()
s.plot(label='integration')
id32.read(TOP, sample, dataset, scan).set_shift(-3.35).plot(label='online processing (centroid)')
br.leg()
br.labels.rixs()
# %%


# %  ===================================================================== %% #
# %  ====================== reading files directly ======================= %% #
# %% ===================================================================== %% #
# one can read the files if the filepath is known. For rixs and images, this is 
# less good then reading using the id32.read() function because metada is not 
# loaded. 

# these functions should be preffered over id32.read() only if the expected 
# file structure of the TOP directory is messed up

scan = 18
ss = id32.readxas(TOP/f'RAW_DATA/sample1_rixs/sample1_rixs_0001/sample1_rixs_0001.h5', scan)

scan = 125
s = id32.readrixs(TOP/'PROCESSED_DATA/align_rixs/align_0001/Online_analysis_align_0001_andor1.spec', scan)
s.plot()

im = id32.readraw(TOP/'RAW_DATA/align_rixs/align_0001/scan0125/andor1_0001.h5')
# %%