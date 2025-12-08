#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Example script for processing rixs data collected at the
at the VERITAS beamline of MAX-IV.

Example data for running this example can be downloaded at this Onedrive link
https://1drv.ms/f/c/6666810d266a3cc9/EjNzJssVTWhFqAOx0-nDCxoBSNt2cvvF8zqOmNdx_zzslQ?e=BOAgx4
Note that some less important files have been removed from this example data to
make the folder lighter and easier to download.

Besides brixs requirements (numpy and matplotlib), brixs.beamlines.veritas 
module also requires h5py.

This is an approximate description of the processing operations for getting a 
rixs spectrum from DLD detector data of VERITAS beamline of MAX-IV

1) photon hits detected after tcutoff are ignored
2) photon hits outside of mask are ignored
3) calculate the histogram of the time bunch to know when bunch is arriving
5) fold time. Useful for finding the time center of the electron bunch
4) apply time mask. Accept only photon hits that are close in time electron photon bunch
5) curvature correction
6) binning of the photon events to get a spectrum
7) energy calibration of the spectrum
8) normalization

The main function of this implementation is the veritas.process() function.

In this function, spectrum is normalized by:

    1) vertical size of the mask
    2) scan exposure_time
    3) number of bins to calculate the spectrum (sbins)

See below the arguments of this function

    Args:
        filepath (str or path): filepath of hdf5 file with rixs data.
        scan (int): scan number
        mask (list or None): only photons within the limits of this mask are 
            considered. Mask must be the following format 
            ([x1_start, x1_stop, y1_start, y1_stop], [x2_start, x2_stop, y2_start, y2_stop], ...]
            if None, mask is not applied. Default is None.

        tcutoff (number, optional): photon hits detected after tcutoff are ignored.
            If None, not time cutoff is applied. Default is 3e7
        tnbins (int, optional): number of bins to calculate the histogram of 
            time_bunch. Default is 10000.
        period (int, optional): electron bunch time period. Default is 1458
        offset (number, optional): if necessary, use this to shift photon hits 
            time before applying the time mask.
        twidth (number, optional): width of the time window around the bunch 
            center. Photons within this window will be accepted. Default is 320
        tcenter (number or string, optional): time "center" of the electron 
            bunch. If number, this number is considered to be the time center. 
            If string, time center is found automatically. Options are: 'gauss' 
            and 'max'. Where one fits a gaussian peak (scipy must be installed)
            and the other simply gets the maximum value over the folded time.
            Default is 'max'

        curv (list or None): 1D array of polynomial coefficients (including 
            coefficients equal to zero) from highest degree to the constant 
            term of the shift values necessary to correct the curvature. If
            None, no curvature correction is applied. Default is None.
        curv_nbins (list): horizontal and vertical number of bins for 
            calculating the curvature. Only used if curv='self'. Default is 
            (20, 1000)

        sbins (int, optional): number of bins for converting photon events to 
            spectrum (number of points in the spectrum). Default is 1200
        calib (number or None, optional): energy calibration value. x axis is 
            multipled by calib. If None, no energy calibration is applied. 
            Default is None.    

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
import brixs.beamlines.veritas as veritas

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
TOP = Path(r'C:\Users\galdin_c\github\brixsexampledata\beamlines\veritas')
RAW = TOP/'raw'

DLD = RAW/'dataset1_DLD.h5'
XAS = RAW/'dataset1.h5'
# %%


# %  ===================================================================== %% #
# %  =========================== rixs advanced =========================== %% #
# %% ===================================================================== %% #

# %% quickly verify detector data processing step-by-step (without curvature correction)
scan = 218
mask0 = ([1041, 6000, 1950, 3780], [1041, 6000, 4330, 5930])  # mask
time  = dict(tcutoff=3e7, tnbins=10000, period=1458, twidth=320, tcenter='max')
sbins = 1200

_ = veritas.verify(scan, DLD, mask=mask0, sbins=sbins, **time)
# %%

# %% calculating the curvature
scan  = 304
mask0 = ([1041, 6000, 1950, 3780], [1041, 6000, 4330, 5930])
time = dict(tcutoff=3e7, tnbins=10000, period=1458, twidth=320, tcenter='max')
nbins = (10, 1000)

curv = veritas.verify_curvature_correction(scan, DLD, mask=mask0, nbins=nbins, **time)
print(f'curv = {list(curv[:2]) + [0]}')
# %%

# %% Energy calibration
scan = 292
mask0 = ([1041, 6000, 1950, 3780], [1041, 6000, 4330, 5930])
time = dict(tcutoff=3e7, tnbins=10000, period=1458, offset=None, twidth=320, tcenter='max')
curv  = [0, -0.04, 0]

# plot for verification
_ = veritas.verify(scan, DLD, mask=mask0, curv=curv, sbins=1200, **time)

# calculate energy calibration via scipy (alternatively, one can also calculate it by hand)
from scipy.signal import find_peaks

# include raw energies as a metadata
veritas.rixs_attrs['raw']['E_raw'] = 'External/beamline_energy/position'

# get spectrum
s = veritas.process(scan, DLD, mask=mask0, curv=curv, sbins=1200, **time)


# delete added metadata so we don't load it anymore for future scans
del veritas.rixs_attrs['raw']['E_raw']

# plot energies
br.figure()
plt.plot(s.E_raw[:, 1])
plt.ylabel('beamline energy')
plt.xlabel('photon hit')

# get energies from plot (or from logbook)
energies = [576, 575, 574, 573, 572, 577, 578, 579, 580]
energies = br.sort(energies, energies)

# find peaks (can also be done manually)
temp, _ = find_peaks(s.y, height=8, distance=100)
centers = [s.x[_] for _ in temp]

# fit peaks to get accurate peak centers
final = []
for center in centers:
    _result = s.fit_peak(limits=(center-200, center+200))
    final.append(_result['popt'][1])

# plot for verification
br.figure()
s.plot()
br.axvlines(final, color='black')
br.labels.rixs()

# calibration value
s = br.Spectrum(x=final, y=energies)
polyfit = s.polyfit(deg=1)
print('calib = ', polyfit['popt'][0])

# plot for verification
br.figure()
s.plot(marker='o', color='black')
polyfit['fit'].plot(color='red')
plt.xlabel('centers')
plt.ylabel('photon energy')
# %%


# %% adding rixs spectra
# define rixs processing parameters
mask0 = ([1041, 6000, 1950, 3780], [1041, 6000, 4330, 5930])
time  = dict(tcutoff=3e7, tnbins=10000, period=1458, twidth=320, tcenter='max')
curv  = [0, -0.04, 0]
calib = 0.002119402
sbins = 1200
parameters = dict(filepath=DLD, mask=mask0, curv=curv, calib=calib, sbins=sbins)
parameters.update(**time)

# read scans
scans = [302, 303, 304]

ss = br.Spectra()
for i, scan in enumerate(scans):
    _s = veritas.process(scan, **parameters)
    motor_x = round(_s.sample_x, 4)
    _s.label = f'{_s.scan}, x={motor_x}'
    ss.append(_s)

# plot
colors = br.get_colors_from_colormap('rainbow', len(ss))
br.figure()
ss.plot(color=colors)
br.leg()
br.labels.rixs()

# align and calculate sum (or average)
ss2 = ss.interp().align().set_shift(-7.853)
s2  = ss2.calculate_average()
s2.label = 'Average'

# plot
br.figure()
ss2.plot(color=colors)
s2.plot(color='black')
br.leg()
br.labels.rixs()
br.zoom(-7, 1)
# %%

# %% energy map (simulated)
# define rixs processing parameters
mask0 = ([1041, 6000, 1950, 3780], [1041, 6000, 4330, 5930])
time  = dict(tcutoff=3e7, tnbins=10000, period=1458, twidth=320, tcenter='max')
curv  = [0, -0.04, 0]
calib = 0.002119402
sbins = 1200
parameters = dict(filepath=DLD, mask=mask0, curv=curv, calib=calib, sbins=sbins)
parameters.update(**time)

# read scans
scans = [302, 303, 304, 306]

# let's pretend this scan sequence is a energy map
ss = br.Spectra()
# energies = []  # uncomment this for real energy map     
for i, scan in enumerate(scans):
    _s = veritas.process(scan, **parameters)
    # energies.append(round(_s.E, 2))  # uncomment this for real energy map
    ss.append(_s)
energies = np.linspace(575, 576, len(scans))  # comment this out for real energy map

# align spectra (sometimes, it must be done mannualy)
# ss2 = ss[i].set_shift(value)  # manually 
ss2 = ss.interp().align().set_shift(-7.853)

# create image
im  = ss2.stack_spectra_as_columns()
im.x_centers = energies

# plot
br.figure()
im.plot()
br.labels.energy_map()
plt.ylim(-6, 1)
# %%