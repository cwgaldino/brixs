#! /usr/bin/env python3
# -*- coding: utf-8 -*-

%load_ext autoreload
%autoreload 2

# standard libraries
import copy
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import importlib
# import brixs as br

# import sys
# sys.path.append(r'C:\Users\galdin_c\Documents\github\brixs')
import brixs as br
br = importlib.reload(br)


plt.ion()
%matplotlib qt5
# %matplotlib inline
# %%

# aug_2021 = Path('/media/galdino/Seagate Expansion Drive/CuSb2O6/CuSb2O6_aug_2021/RIXS')
# oct_2021 = Path('/media/galdino/Seagate Expansion Drive/CuSb2O6/CuSb2O6_oct_2021/RIXS')
aug_2021 = Path(r'D:\CuSb2O6\CuSb2O6_aug_2021\RIXS')
fixtures = Path(r'C:/Users/galdin_c/Documents/github/brixs/local/fixtures/')

prefix = 'Cu'
n = 18
filepath = fixtures/'Cu_0018_d1.h5'

# %%


import brixs as br
import numpy as np


# package is based on three objects
pe = PhotonEvents()
s  = Spectrum()
ss = Spectra()



# ==============================================================================
# %% TLDL (too long didn't read) ===============================================
# ==============================================================================

# %% spectra alignment
filelist = br.backpack.filelist('../fixtures/calib_example')
ss = br.Spectra()
for f in filelist:
    ss.append(data=np.loadtxt(f, delimiter=','))

plt.figure()
ss.plot(vertical_increment=0.1)
ss.calculate_shifts(ref=0, mode='cc')
ss.set_shifts()
plt.figure()
ss.plot(vertical_increment=0.1)
ss.crop(stop=600)
plt.figure()
ss.plot(vertical_increment=0.1)

# ==============================================================================
# %% find and fit peaks (easy)
filepath = '../fixtures/peak_fit/easy.dat'
s = br.Spectrum(data=np.loadtxt(filepath, delimiter=','))
s.find_peaks()
plt.figure()
s.plot()
s.plot_detected_peaks()

s.fit_peak(0)
plt.figure()
s.plot()
s.fit.plot()

s.fit_peak()
plt.figure()
s.plot()
s.fit.plot()

# results
s.fit_data[0]
print('elastic peak at', s.fit_data[0]['c'])
print('first excitation at', s.fit_data[1]['c'])


# ==============================================================================
# %% find and fit peaks (medium)
# shoulder is not detected
filepath = '../fixtures/peak_fit/medium.dat'
s = br.Spectrum(data=np.loadtxt(filepath, delimiter=','))
s.find_peaks()
plt.figure()
s.plot()
s.plot_detected_peaks()

# fit is wrong
s.fit_peak(0)
plt.figure()
s.plot()
s.fit.plot()

# fit is good
s.fit_peak(0, multiplicity=2)
plt.figure()
s.plot()
s.fit.plot()

# results
print('elastic peak at', s.fit_data[0]['c'])

# ==============================================================================
# %% find and fit peaks (hard)
filepath = '../fixtures/peak_fit/hard.dat'
s = br.Spectrum(data=np.loadtxt(filepath, delimiter=','))
s.find_peaks()
plt.figure()
s.plot()
s.plot_detected_peaks()

# fit is good
s.fit_peak(multiplicity={0:2})
plt.figure()
s.plot()
s.fit.plot()

# results
print('elastic peak at', s.fit_data[0]['c'])

# ==============================================================================
# %% energy calibration (answer is 10 meV/point)
filelist = br.backpack.filelist('../fixtures/calib_example')
ss = br.Spectra()
for f in filelist:
    ss.append(data=np.loadtxt(f, delimiter=','))

calib = ss.calculate_calib(start=930, stop=939)
print(calib*1000, ' meV/point')

plt.figure()
ss.plot_disp()

# ==============================================================================
# %% energy calibration using peak fit
filelist = br.backpack.filelist('../fixtures/calib_example')
ss = br.Spectra()
for f in filelist:
    ss.append(data=np.loadtxt(f, delimiter=','))

ss.find_peaks(prominence=0.5)
calib = ss.calculate_calib(start=930, stop=939, mode='peak', idx=0)
print(calib*1000, ' meV/point')

plt.figure()
ss.plot(vertical_increment=0.1)
for i, s in enumerate(ss):
    s.fit.plot(offset=-0.1*i, color='black')

plt.figure()
ss.plot_disp()


# ==============================================================================
# %% fake data =================================================================
# ==============================================================================

# creating a sequence of fake spectra
positions = np.linspace(100, 1000, 10)
energies = np.linspace(930, 939, 10)
ss = br.Spectra()
for c in positions:
    fake = br.fake(amp=1, c=c, fwhm=10)
    s = fake.get_spectrum(0, 1500, n_points=6000, noise=3)
    ss.append(s)

plt.figure()
ss.plot(vertical_increment=0.1)

ss.save(r'../fixtures/calib_example')

# ==============================================================================
# create spectrum with excitations (multiple peaks)
fake1 = br.fake(1, 0, 0.1, [[0.5, 3, 0.2], [0.7, 4, 0.2]])
s1 = fake1.get_spectrum(noise=3)
# big elastic peak with low energy excitation
fake2 = br.fake(1, 0, 0.1, [[0.5, 0.15, 0.2], [0.7, 4, 0.2]])
s2 = fake2.get_spectrum(noise=3)
# small elastic peak with intense low energy excitation
fake3 = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
s3 = fake3.get_spectrum(noise=3)

plt.figure()
s1.plot(label='easy')
s2.plot(offset=0.5, label='medium')
s3.plot(offset=1, label='hard')
plt.legend()

s1.save(r'../fixtures/peak_fit/easy.dat')
s2.save(r'../fixtures/peak_fit/medium.dat')
s3.save(r'../fixtures/peak_fit/hard.dat')




# %% functions to read data from ADRESS
pe, s, nd = br.read_ADRESS(filepath)

s   = br.get_spectrum_ADRESS(filepath)
pe  = br.get_pe_ADRESS(filepath)
bad = br.get_bad_ADRESS(filepath)

ss, pes = br.get_scan_ADRESS(filepath, prefix, n, zfill=4)
ss      = br.get_scan_spectrum_ADRESS(filepath, prefix, n, zfill=4)
pes     = br.get_scan_pe_ADRESS(filepath, prefix, n, zfill=4)


# %% instrument/scan parameters
pe, s0, nd = br.read_ADRESS(filepath)
print(nd.keys())
print(nd['AcquireTime'])
print(nd['PhotonEnergy'])
print(nd['RMUCurrent'])
print('......')


# %% photon events with bad events
pe, s0, nd = br.read_ADRESS(filepath)
ax = pe.plot()
bad = br.get_bad_ADRESS(filepath)
bad.plot(ax, mfc='red', pointsize=2)


# %% Curvature correction (offset correction)
pe = br.get_pe_ADRESS(filepath)
pe.plot()
plt.ylim(820, 920)
plt.title('1) Raw photon events')
pe.bins = (10, 1000)
pe.plot(show_bins=(True, False))
plt.ylim(820, 920)
plt.title('2) X bins')
pe.calculate_offsets(ref=5)
pe.fit_offsets(deg=2)
ax1 = pe.plot_offsets()
pe.plot_fit(ax1)
plt.title('3) Offsets and fit')
pe.plot(show_offsets=True, show_fit=True)
plt.ylim(820, 920)
plt.title('4) Offsets and fit over the image')
pe.offsets_correction()
pe.plot()
plt.ylim(820, 920)
plt.title('5) After correction')

# %% plotting columns
pe = br.get_pe_ADRESS(filepath)
pe.bins = (10, 1000)
pe.plot_columns()
plt.title('1) before correction')
plt.xlim(700, 1000)

pe.plot_columns(vertical_increment=5)
plt.title('2) before correction (cascaded)')
plt.xlim(700, 1000)

pe.calculate_offsets(ref=5)
pe.fit_offsets(deg=2)
pe.offsets_correction()
pe.plot_columns(vertical_increment=5)
plt.title('2) after correction')
plt.xlim(700, 1000)

# %% different ways of setting binning
pe = br.get_pe_ADRESS(filepath)
pe.bins = (10, 1000)
pe.bins = 100
pe.bins_size = (150, 2)
pe.bins_size = 10
pe.set_bins(10, 1000)
pe.set_bins((10, 1000))
pe.set_bins(100)
pe.set_bins_size(150, 2)
pe.set_bins_size(10)

# %% read binnig value
print(pe.bins)
print(pe.get_bins())
print(pe.get_bins_size())

# %% guess binning for curvature correction
# pe.bins = 'guess' will find the minimal binning so two adjacent bins have
# different offsets
pe = br.get_pe_ADRESS(filepath)
pe.bins = (9, 300)
pe.calculate_offsets(ref=0)
pe.fit_offsets(deg=2)
ax1 = pe.plot_offsets()
pe.plot_fit(ax1)
plt.title('1) low binning (repeated offsets)')

pe = br.get_pe_ADRESS(filepath)
pe.bins = (9, 500)
pe.calculate_offsets(ref=0)
pe.fit_offsets(deg=2)
ax1 = pe.plot_offsets()
pe.plot_fit(ax1)
plt.title('2) almost good binning (two repeated offsets)')

pe = br.get_pe_ADRESS(filepath)
pe.bins = 'guess'
pe.calculate_offsets(ref=0)
pe.fit_offsets(deg=2)
ax1 = pe.plot_offsets()
pe.plot_fit(ax1)
plt.title('3) guessed binning (no repeated offsets)')

print('guessed binning: ', pe.bins)

# %% Calculate spectrum
pe = br.get_pe_ADRESS(filepath)
pe.bins = (9, 650)
pe.calculate_offsets(ref=0)
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.plot()

plt.title(f'binning = {pe.spectrum_bins[1]}')
plt.ylabel('Intensity')
plt.xlabel('bins')

# %% calculate spectrum in terms of lenght of the detector
pe = br.get_pe_ADRESS(filepath)
pe.bins = (9, 650)
pe.calculate_offsets(ref=0)
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000, mode='lenght')
s.plot()

plt.title(f'detector size = ({pe.x_max}, {pe.y_max})')
plt.ylabel('Intensity')
plt.xlabel('lenght')

# %% Try to find and fit the elastic line
# first it smooths the data (moving average)
# them rouglhy finds peaks on the smoothed data (find_peaks)
# them finds peaks on the original data
aug_2021 = Path(r'D:\CuSb2O6\CuSb2O6_aug_2021\RIXS')
filepath = aug_2021/'O_0100_d2.h5'
# filepath = fixtures/'Cu_0018_d3.h5'  ## SERGIO

s = br.get_spectrum_ADRESS(filepath)
s.guess_elastic_peak()
s.plot(show_fit=True, show_fit_range=True)

help(br.Spectra.append)
ss.calculate_shifts()
ss.shifts
ss.shifts_length
ss.shift_ranges
ss.plot(show_ranges=True)
# %%
s.plot()

pe = br.get_pe_ADRESS(filepath)
pe.bins = (9, 650)
pe.calculate_offsets(ref=0)
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.plot(marker='o')
x, y = s.guess_elastic_peak(moving_average_winddow=8, asymmetry=True)
ax = s.plot()
ax.plot(x, y)
s.plot(show_fit=True, show_fit_range=True)
print('amplitude: ', s.elastic_amp)
print('fwhm: ', s.elastic_w)
print('position: ', s.elastic_c)

# ==============================================================================
# %% Guess elastic peak position ===============================================
# ==============================================================================

# easy
fake = br.fake(1, 0, 0.1, [[0.5, 3, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left')
s.plot(show_fit=True, show_fit_range=True)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# excitation close to the elastic peak
fake = br.fake(1, 0, 0.1, [[0.5, 0.2, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# same as previous, but fitting only the left part
fake = br.fake(1, 0, 0.1, [[0.5, 0.2, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', fit_full=False, verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# excitation very close to the elastic peak
fake = br.fake(1, 0, 0.1, [[0.5, 0.2, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# excitation very close to the elastic peak (fitting only the left part of the peak)
fake = br.fake(1, 0, 0.1, [[0.5, 0.1, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', fit_full=False, verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# excitation crazy close (uses two peaks)
fake = br.fake(1, 0, 0.1, [[0.5, 0.05, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# same as previous, but fitting only the left part (not as good as the previous one)
fake = br.fake(1, 0, 0.1, [[0.5, 0.05, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', fit_full=False, verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# Very small elastic peak and a very intense close excitation
fake = br.fake(0.5, 0, 0.1, [[1, 0.2, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# same as previous, but fitting only the left part
fake = br.fake(0.5, 0, 0.1, [[1, 0.2, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', fit_full=False, verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# Very small elastic peak and a very intense and very close excitation
fake = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# Very small elastic peak and a very intense and very close excitation (fit only left part)
fake = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', fit_full=False, verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# even closer!!!
fake = br.fake(0.5, 0, 0.1, [[1, 0.1, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# of course, if the elastic peak is too overlaped with the first exctitation
# it becames impossible to separate the two peaks
fake = br.fake(0.5, 0, 0.1, [[2, 0.1, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.guess_elastic_peak(side='left', verbose=True)
ax = s.plot(show_fit=True, show_fit_range=True, marker='o', ms=4)
s.plot_fit_full(ax)
title = f'REAL: amp = {round(fake.elastic_amp, 3)}, position = {round(fake.elastic_c, 3)}, fwhm = {round(fake.elastic_fwhm, 3)}\n ' +\
        f'FIT: amp = {round(s.elastic_amp, 3)}, position = {round(s.elastic_c, 3)}, fwhm = {round(s.elastic_fwhm, 3)}\n '
plt.title(title)

# let's try some real data (very few data points)
pe = br.get_pe_ADRESS(filepath)
pe.bins = (9, 650)
pe.calculate_offsets(ref=0)
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(1000)
s.plot(marker='o')
s.guess_elastic_peak(verbose=True)
s.plot(show_fit=True, show_fit_range=True)

# regular data
s = br.get_spectrum_ADRESS(filepath)
s.plot()
s.guess_elastic_peak(side='r', verbose=True, recursive=False)
s.plot(show_fit=True, show_fit_range=True)

s = br.get_spectrum_ADRESS(filepath)
s.calib = 0.01
s.find_peaks2()
s.peaks
s.min_peaks


fake = br.fake(0.5, 0, 0.1, [[1, 0.5, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.find_peaks()
s.plot()
print(s.peaks)
s.fit_peak(idx=0)
# s.fit_peak(idx=[0, 1])
# s.fit_peak(idx=[0, 1, 2])
# s.fit_peak(idx=[0, 2], )
# s.fit_peak(idx=[0, 2], ranges=((-1.5, 0.2), (1, 6)))
plt.plot(s.x, s.elastic_func_full(s.x))
plt.plot(s.x2fit, s.y2fit)


fake = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.find_peaks()
s.plot()
print(s.peaks)
# s.fit_peak(idx=0)
# s.fit_peak(idx=0, asymmetry=True)
# s.fit_popt
# s.fit_popt
# s.

fake = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.find_peaks()
s.fit_peak(idx=0, multiplicity=2)
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
print(s.fit_data[0]['c'])
print(s.peaks)
s.calib = 2
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
print(s.fit_data[0]['c'])
print(s.peaks)

fake = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.find_peaks()
s.fit_peak(idx=0, multiplicity=2)
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
print(s.fit_data[0]['c'])
print(s.peaks)
s.shift = 2
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
print(s.fit_data[0]['c'])
print(s.peaks)


fake = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.find_peaks()
s.fit_peak(idx=0, multiplicity=2)
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
print(s.fit_data[0]['amp'])
print(s.peaks)
s.offset = 2
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
max(s.y)
print(s.fit_data[0]['amp'])
print(s.peaks)



fake = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.find_peaks()
s.fit_peak(idx=0, multiplicity=2)
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
print(s.fit_data[0]['amp'])
print(s.peaks)
s.factor = 1
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
max(s.y)
print(s.fit_data[0]['amp'])
print(s.peaks)



fake = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.find_peaks()
s.fit_peak(idx=0, multiplicity=2)
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
print(s.fit_data[0]['amp'])
print(s.peaks)
s.offset = 2
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
s.floor()
s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
s.offset

# %% ===========================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================



positions = np.linspace(100, 1000, 10)
energies = np.linspace(930, 939, 10)
# print(positions)  # every 100 points there will be a peak
# print(energies)   # every peak is 1 eV shifted
# print(f'dispersion should be {1/100 * 1000} meV/point')
ss = br.Spectra()
for c in positions:
    fake = br.fake(amp=1, c=c, fwhm=10)
    ss.append(fake.get_spectrum(0, 1500, n_points=6000, noise=2))


# ss.plot(vertical_increment=0.1)
# ss.calculate_shifts(ref=0, mode='cc')
# ss.set_shifts()
# ss.plot(vertical_increment=0.1)
# ss.crop(stop=600)
# ss.plot(vertical_increment=0.1)
# ss.check_same_x()   # same x

ss.plot(vertical_increment=0.1)
# ss[4].x[0] = -1000
ss.calculate_shifts(ref=0, mode='max')
ss.set_shifts()
ss.plot(vertical_increment=0.1)
ss.crop()
ss.plot(vertical_increment=0.1)
# ss.check_same_x()   # not necessarily same x
ss.interp()
ss.check_same_x()   # same x
ss.plot(vertical_increment=0.1)

# ss.plot(vertical_increment=0.1)
# ss.find_peaks()
# ss.calculate_shifts(ref=0, mode='peak')
# ss.set_shifts()
# ss.plot(vertical_increment=0.1)
# ss.crop()
# ss.plot(vertical_increment=0.1)
# # ss.check_same_x()   # not necessarily same x
# ss.interp()
# ss.check_same_x()   # same x
# ss.plot(vertical_increment=0.1)











# %% ===========================================================================

fake = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
s = fake.get_spectrum(noise=1)
s.find_peaks()
# s.plot()
print(s.peaks)
[k for k in s.peaks]
s.fit_peak()
s.fit_peak(idx=0, multiplicity=2)
# s.fit_peak(idx=[0, 1], multiplicity={0:2, 1:1})
# s.fit_peak(idx=[0, 1], multiplicity={0:2})
# s.fit_peak(idx=[0, 1], multiplicity={0:2}, fixed_m={1:False})
# plt.plot(s.x, s.fit_model(s.x))
ax=s.plot(show_fit=True, show_fit_ranges=True, show_fit_partial=True)
# s.mark_peaks(ax, show_width=False)
s.mark_peaks(ax)

print(s.peaks)
s.append_peak(10, -2, 1)

a = [3,2,3]
np.diff(a)
any(a < 0)

max(s.peaks.keys())
 s.fit_data[0]['func'] ==  s.fit_data[1]['func']
plt.plot(s.x, s.fit_data[0]['func'](s.x))
s.fit_data[1]['popt']
s.fit_data[0]['popt']


sum(s.fit_data[0]['func'](s.x) - s.fit_data[1]['func'](s.x))

sum(s.fit_data[0]['_func'](s.x, *s.fit_data[0]['popt']) - s.fit_data[1]['_func'](s.x, *s.fit_data[1]['popt']))


# %% Calibrate and shift spectrum
dispersion = 0.012796  # eV/bin
s = br.get_spectrum_ADRESS(filepath)
s.calib(dispersion)
s.guess_elastic_peak(asymmetry=True)
s.apply_shift(-s.elastic_c, mode='x')
s.plot()

# %% "best" binning for curvature correction
# br.best_bins() will try many binnings and return the one that yields the
# sharpest elastic line
pe = br.get_pe_ADRESS(filepath)
best_bins, best, fwhm = pe.best_bins(asymmetry=False)
print('"best" binning: ', best_bins, 'with elastic line with fwhm= ', round(best, 2))


# %% different binnings comparisson
pe = br.get_pe_ADRESS(filepath)
pe.bins = 'guess'
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(asymmetry=False)
ax = s.plot(label=f'{pe.bins} fwhm= {round(s.elastic_w, 2)}')

pe = br.get_pe_ADRESS(filepath)
pe.bins = 'guess'
pe.bins = pe.bins + pe.bins*0.4
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(asymmetry=False)
s.plot(ax, vertical_increment=5, label=f'{pe.bins} fwhm= {round(s.elastic_w, 2)}')

pe = br.get_pe_ADRESS(filepath)
pe.bins = 'guess'
pe.bins = pe.bins + pe.bins*0.7
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(asymmetry=False)
s.plot(ax,  vertical_increment=10, label=f'{pe.bins} fwhm= {round(s.elastic_w, 2)}')

pe = br.get_pe_ADRESS(filepath)
pe.bins = (15, 6000)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(asymmetry=False)
s.plot(ax, vertical_increment=15, label=f'{pe.bins} fwhm= {round(s.elastic_w, 2)}')

pe = br.get_pe_ADRESS(filepath)
pe.bins = (60, 8000)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(asymmetry=False)
s.plot(ax, vertical_increment=20, label=f'{pe.bins} fwhm= {round(s.elastic_w, 2)}')

plt.xlim(3200, 3600)
plt.legend()

# %% Binning multiplicity times doesn't necessarily yield a good result
# if the first correction is good enough, there's no need to correct it again
pe = br.get_pe_ADRESS(filepath)
pe.bins = (18, 2000)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak()
ax = s.plot(label=f'1 correction fwhm= {round(s.elastic_w, 2)}')

pe = br.get_pe_ADRESS(filepath)
for i in (0, 1):
    pe.bins = (18, 2000)
    pe.calculate_offsets(ref=int(pe.bins[0]/2))
    pe.fit_offsets(deg=2)
    pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak()
s.plot(ax, label=f'2 corrections fwhm= {round(s.elastic_w, 2)}')

pe = br.get_pe_ADRESS(filepath)
for i in (0, 1, 2):
    pe.bins = (18, 2000)
    pe.calculate_offsets(ref=int(pe.bins[0]/2))
    pe.fit_offsets(deg=2)
    pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak()
s.plot(ax, label=f'3 corrections fwhm= {round(s.elastic_w, 2)}')

plt.xlim(3200, 3600)
plt.legend()

# %% comparisson between spectrum from beamline and from photon events list
s = br.get_spectrum_ADRESS(filepath)
s.guess_elastic_peak()
s.plot(show_elastic=True)
plt.xlim(3200, 3600)
plt.title('1) spectrum processed by the beamline')

pe  = br.get_pe_ADRESS(filepath)
pe.bins = (60, 6000)  # calculated best binning
pe.calculate_offsets(ref=30)
pe.fit_offsets(deg=2)
pe.offsets_correction()
s2 = pe.calculate_spectrum(bins=6000)
s2.guess_elastic_peak()
s2.plot(show_elastic=True)
plt.xlim(3200, 3600)
plt.title('2) spectrum from photon events list')

s.calib(dispersion)
s.apply_shift(-s.elastic_c, mode='x')
ax = s.plot(label=f'beamline fwhm= {round(s.elastic_w, 2)}')

s2.calib(dispersion)
s2.apply_shift(-s2.elastic_c, mode='x')
s2.plot(ax, label=f'from pe fwhm= {round(s2.elastic_w, 2)}')
plt.legend()
plt.xlim(-3, 1)
plt.title('3) final spectrum')


# %% we can best binning FOR SPECTRUM CALCULATION
# it tries many binning and measure the sharpness of the elastic line
pe = br.get_pe_ADRESS(filepath)
pe.bins = (15, 6000)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
best_bins, best, fwhm = pe.best_bins_for_spectrum(asymmetry=False)
print('"best" binning: ', best_bins, 'with elastic line with fwhm= ', round(best, 2))


# %% no need to use a binning higher that the calculated one
# values higher than the calculated seems to yield the same result
# gotta try that with other data (TODO)
pe = br.get_pe_ADRESS(filepath)
pe.bins = (60, 6000)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()

for y_bin in [1000, 2000, 3000, 4000, 6000, 7000]:
    s = pe.calculate_spectrum(bins=3000, mode='lenght')
    ax = s.plot(ms=4, marker='o', label = best_bins)
    plt.title(f'comparing binnig = 3000 with binning = {y_bin}')

    s = pe.calculate_spectrum(bins=y_bin, mode='lenght')
    s.y = s.y*y_bin/3000
    if y_bin > 7000:
        s.plot(ax, label=(1, y_bin))
    else:
        s.plot(ax, ms=4, marker='o', label=(1, y_bin))

    plt.legend()
    plt.xlim(878, 890)

# %%  Dealing with multiplicity spectra
# ADRESS beamline has three ccds so we have to add data to have the final spectrum
filenames = ['Cu_0018_d1.h5', 'Cu_0018_d2.h5', 'Cu_0018_d3.h5']
ss = br.Spectra()
for filename in filenames:
    s = br.get_spectrum_ADRESS(fixtures/filename)
    ss.append(s)
ss.plot(vertical_increment=5)

# %% same with photon events
filenames = ['Cu_0018_d1.h5', 'Cu_0018_d2.h5', 'Cu_0018_d3.h5']
pes = []
for filename in filenames:
    pes.append(br.get_pe_ADRESS(fixtures/filename))

# %% builtin functon for multiplicity spectra from same scan
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)

# %% individual spectrum can be accesed via indexing
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
print(ss.get_spectra_count())
ss[0]
ss[1]
ss[2]

ss[0].plot()
plt.title('Ploting the first spectrum')

ss.plot(idx=0)
plt.title('Also ploting the first spectrum')

# %% spectra can be removed
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
print(ss.get_spectra_count())
ss.remove(0)  # this will also remove shifts (if it was calculated)
print(ss.get_spectra_count())

# %% aligning spectra and getting final spectrum
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
ss.plot(vertical_increment=5)
plt.title('1) raw')
ss.calculate_shifts()
ss.apply_shifts()
ss.plot(vertical_increment=5)
plt.title('2) after shift correction')
s0 = ss.calculate_sum()
s0.guess_elastic_peak()
s0.plot(show_elastic=True, label = f'final spectrum fwhm= {round(s0.elastic_w, 2)}')
plt.title('3) sum of all 3 spectra')
plt.legend()

# %% shifts can be calculated using three modes
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
ss.calculate_shifts(mode='cc')  # x point must be equally spaced
ss.apply_shifts()
ss.plot(vertical_increment=5)
plt.title('1) via cross-correlation')

ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
ss.calculate_shifts(mode='max')
ss.apply_shifts(mode='x')
ss.plot(vertical_increment=5)
plt.title('2) via max')

ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
ss.calculate_shifts(mode='elastic')
ss.apply_shifts(mode='x')
ss.plot(vertical_increment=5)
plt.title('3) via elastic line')

# %% three shift modes
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
ss.calculate_shifts(mode='cc')
ss.apply_shifts(mode='x')     # shifts the x axis. y data is kept intact
ss.apply_shifts(mode='y')     # interpolates y data. x axis is kept intact
ss.apply_shifts(mode='roll')  # x and y are preserved, but x data must be equally space
# if mode='roll', shift values must be an integer

# %% one can check is x data is equaly spaced for all spectra
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
ss.check_step_x()  # default max error is 1% of the average data step
print('data is equally spaced with step = ', ss.step)
# if ss.step is defined, calculating shifts with 'max' and 'elastic' will also
# return integer shift values so you can use ss.apply_shift() with mode='roll'
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
ss.calculate_shifts(mode='max')
print(ss.shifts)         # not defined
print(ss.shifts_lenght)  # defined
ss.check_step_x()        # checking x step
print(ss.shifts)         # now integer steps are defined

# %% for adding two spectra, x data must be the same
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
ss.check_same_x()  # default max error is 1% of the data step
# this function automatically calls ss.check_step_x() is step is not defined

# %% built in function for photon events
pes = br.get_scan_pe_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
print('pes it is a list: ', type(pes))
ss = br.Spectra()
for pe in pes:
    # pe.bins = (20, 3000)  # 9.28
    # pe.bins = (15, 4000)  # 9.36
    pe.bins = (10, 6000)  # 9.44
    pe.calculate_offsets(ref=5)
    pe.fit_offsets()
    pe.offsets_correction()
    s = pe.calculate_spectrum(bins=6000)
    ss.append(s)
ss.calculate_shifts()
ss.apply_shifts()
s1 = ss.calculate_sum()
s1.guess_elastic_peak()
s1.plot(show_elastic=True, label = f'final spectrum fwhm= {round(s1.elastic_w, 2)}')
plt.title('final spectrum from photon events list')
plt.legend()

# %% using best bins
filenames = ['Cu_0018_d1.h5', 'Cu_0018_d2.h5', 'Cu_0018_d3.h5']
pes = []
for filename in filenames:
    pes.append(br.get_pe_ADRESS(fixtures/filename))
# best_bins, best, fwhm = pes[0].best_bins(asymmetry=False)  # (15, 6000)
# print(best_bins)
# best_bins, best, fwhm = pes[1].best_bins(asymmetry=False)  # (7, 600)
# print(best_bins)
# best_bins, best, fwhm = pes[2].best_bins(asymmetry=False)  # (20, 8000)
# print(best_bins)

bins = [(15, 6000), (7, 600), (20, 8000)]
ss = br.Spectra()
for i, pe in enumerate(pes):
    pe.bins = bins[i]
    pe.calculate_offsets(ref=int(pe.bins[0]/2))
    pe.fit_offsets()
    pe.offsets_correction()
    s = pe.calculate_spectrum(bins=6000)
    ss.append(s)
ss.calculate_shifts()
ss.apply_shifts()
s2 = ss.calculate_sum()
s2.guess_elastic_peak(asymmetry=False)
s2.plot(show_elastic=True, label = f'final spectrum fwhm= {round(s2.elastic_w, 2)}')
plt.legend()
plt.title('final spectrum from pe using best bins')


# %% removing the 3rd ccd
filenames = ['Cu_0018_d1.h5', 'Cu_0018_d2.h5']
pes = []
for filename in filenames:
    pes.append(br.get_pe_ADRESS(fixtures/filename))
bins = [(15, 6000), (7, 600)]
ss = br.Spectra()
for i, pe in enumerate(pes):
    pe.bins = bins[i]
    pe.calculate_offsets(ref=int(pe.bins[0]/2))
    pe.fit_offsets()
    pe.offsets_correction()
    s = pe.calculate_spectrum(bins=6000)
    ss.append(s)
ss.calculate_shifts()
ss.apply_shifts()
s3 = ss.calculate_sum()
s3.guess_elastic_peak(asymmetry=False)
s3.plot(show_elastic=True, label = f'final spectrum fwhm= {round(s3.elastic_w, 2)}')
plt.legend()
plt.title('final spectrum from pe using best bins')

# %% Comparing different methods of getting the final spectrum
# apparently there is no difference
# but removing the 3rd ccd makes the elastic line a little sharper
ss = br.Spectra([s0, s1])
ss.calculate_shifts()
ss.apply_shifts()
ss.plot()

ss = br.Spectra([s0, s2])
ss.calculate_shifts()
ss.apply_shifts()
ss.plot()

ss = br.Spectra([s0, s3])
ss[1].y = ss[1].y*1.4
ss.calculate_shifts()
ss.apply_shifts()
ss.plot()

# %% energy calibration (Spectra)
dispersion = 0.012796
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
ss.calculate_shifts()
ss.apply_shifts()
ss.calculate_sum()
ss.sum.calib(value=dispersion)
ss.sum.guess_elastic_peak()
elastic_position = ss.sum.elastic_c
ss.sum.apply_shift(value=-elastic_position, mode='x')
ss.sum.plot()
plt.title('2) final calibrated spectrum')

# %% defining macros
# pe, s, and ss are dynamic objects. We can define new methods on the fly.
# lets say we want to repeat the previous calculation for other scans
# then we can define the following function and link it to each scan
def get_final_spectrum(self):
    self.calculate_shifts()
    self.apply_shifts()
    s = self.calculate_sum()
    s.calib(value=0.012796)
    s.guess_elastic_peak()
    elastic_position = s.elastic_c
    s.apply_shift(value=-elastic_position, mode='x')
    return s

# linking macro to ss
br.Spectra.get_final_spectrum = get_final_spectrum

# now we can use the macro
ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
s = ss.get_final_spectrum()
s.plot()

# %% save spectra, save spectrum, save photon events
pe = br.get_pe_ADRESS(fixtures/filename)
pe.save('test_pe.dat')

s = br.get_spectrum_ADRESS(fixtures/filename)
s.save('test_s.dat')

ss = br.get_scan_spectrum_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
ss.save(prefix='test_ss_', suffix='.dat', zfill=4)

# %% Use alignment spectra to fix curvature for all subsequant spectra
# scans 1-4 are alignment scans (using carbon tape)
# scan 4 is the one we used to find the curvature (higher intensity)
scan = 4
pes = br.get_scan_pe_ADRESS(folderpath=fixtures, prefix='Cu_', n=scan, zfill=4)
pes[0].plot() # only difuse scatering

pes[0].bins = 'guess'  # guessed is (9, 600)
pes[0].calculate_offsets()
pes[0].fit_offsets()
pes[0].offsets_correction()
s = pes[0].calculate_spectrum(6000)
s.plot()

# check best
best_bins, best, fwhm = pe.best_bins(asymmetry=False)  # best is (7, 600)

# %% in the end, if we use (7, 600) or (20, 3000) as they yield resonable spectrum
ss = br.Spectra()

pes = br.get_scan_pe_ADRESS(folderpath=fixtures, prefix='Cu_', n=4, zfill=4)
pes[0].bins = (7, 600)
pes[0].calculate_offsets()
pes[0].fit_offsets()
pes[0].plot(show_fit=True)
plt.ylim(820, 850)
pes[0].offsets_correction()
ss.append(pes[0].calculate_spectrum(6000))

pes = br.get_scan_pe_ADRESS(folderpath=fixtures, prefix='Cu_', n=4, zfill=4)
pes[0].bins = (20, 3000)
pes[0].calculate_offsets()
pes[0].fit_offsets()
pes[0].plot(show_fit=True)
plt.ylim(820, 850)
pes[0].offsets_correction()
ss.append(pes[0].calculate_spectrum(6000))

ss.calculate_shifts()
ss.apply_shifts()
ss.plot()
plt.title('both binnings are ok')

# %% Get curvature function and apply to another spectrum
ss = br.Spectra()

pes = br.get_scan_pe_ADRESS(folderpath=fixtures, prefix='Cu_', n=4, zfill=4)
pes[0].bins = (20, 3000)
pes[0].calculate_offsets()
pes[0].fit_offsets()
pes[0].offsets_correction()
f0 = pes[0].offsets_func
print(f0)
print(pes[0].offsets_popt)  # 2nd order degree polynomial yields three parameters

# apply to another scan
pes2 = br.get_scan_pe_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
pes2[0].plot()
plt.ylim(820, 920)
pes2[0].offsets_correction(f=f0)
pes2[0].plot()
plt.ylim(820, 920)
ss.append(pes2[0].calculate_spectrum(6000))

# compare with self offset calculation
pes3 = br.get_scan_pe_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
pes3[0].bins = (20, 3000)
pes3[0].calculate_offsets()
pes3[0].fit_offsets()
pes3[0].offsets_correction()
ss.append(pes3[0].calculate_spectrum(6000))

ss.calculate_shifts()
ss.apply_shifts()
ss.plot()
plt.title('both methods are ok')

# %% save function
# sadly, you can't save the function itself, but you can save its parameters
# so, if you know the function definition, you can reproduce it
# as default, it understandas that the correction function is polynomial
pes = br.get_scan_pe_ADRESS(folderpath=fixtures, prefix='Cu_', n=4, zfill=4)
pes[0].bins = (20, 3000)
pes[0].calculate_offsets()
pes[0].fit_offsets()
print(pes[0].offsets_popt)  # deg=2 polynomial. item from lower to higher coeficient
pes[0].save_popt(filepath='curvature_correction.dat', comments='2nd order polynomial')

pes2 = br.get_scan_pe_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
pes2[0].plot()
plt.ylim(820, 930)
print(pes2[1].offsets_func)
pes2[0].load_popt(filepath='curvature_correction.dat')
print(pes2[0].offsets_func)
pes2[0].offsets_correction()
pes2[0].plot()
plt.ylim(820, 900)

# %% RIXS MAP ==================================================================
# Momentum and incident energy maps are created by the exact same method

# fake map
incident_energies = np.linspace(930, 950, 20)
amps = abs(np.sin(np.linspace(0, np.pi, 20)))
excitations = [[amps[i]*0.5, incident_energies[i]*4-3000, 30] for i in range(20)]
ss = br.Spectra()
for i in range(20):
    fake = br.fake(amp=amps[i], c=500, w=10, excitations=excitations[i])
    ss.append(fake.get_spectrum(0, 1500, n_points=6000))

# check spectra
ss.plot(vertical_increment=0.1)

# plot map
ss.plot_map(values=incident_energies)


# %% calculate dispersion ======================================================
# Use fake spectra to simulate a 10 meV/point dispersion
positions = np.linspace(100, 1000, 10)
energies = np.linspace(930, 939, 10)
print(positions)  # every 100 points there will be a peak
print(energies)   # every peak is 1 eV shifted
print(f'dispersion should be {1/100 * 1000} meV/point')
ss = br.Spectra()
for c in positions:
    fake = br.fake(amp=1, c=c, w=10)
    ss.append(fake.get_spectrum(0, 1500, n_points=6000))

# check spectra
ss.plot(vertical_increment=0.1)

# calculate dispersion
disp = ss.calculate_dispersion(start_energy=930, stop_energy=939)
print(f'dispersion = {round(disp*1000, 3)} meV/point')

# we can check dispersion attributes
ss.disp  # dispersion
ss.disp_postions  # elastic peak positions
ss.disp_energies  # energies
ss.disp_popt  # polynomial coeficients of fit
ss.disp_func # polynomial function of fit

# plot dispersion fit
ss.plot_disp()
plt.xlabel('Energy (eV)')
plt.ylabel('Elastic peak position (points)')




# %%
folderpath = Path(r'D:\CuSb2O6\CuSb2O6_aug_2021\RIXS')
filelist = br.backpack.filemanip.filelist(folderpath)

for filepath in filelist:







# %% Calculate dispersion ADRESS ===============================================
disp, disp_bin, sss = br.get_dispersion_ADRESS(folderpath=fixtures, prefix='Cu_', start_energy=925, stop_energy=935, start_scan=5, stop_scan=15)
print('dispersion = ', np.mean(disp)*1000, 'meV/point')    # energy/(detector lenght)  --> binning indepedent
print('dispersion = ', np.mean(disp_bin)*1000, 'meV/bin')  # energy/bin --> use this for calibrating spectrum with binning 6000

# plot
ax = sss[0].plot_disp(label=f'{round(sss[0].disp*1000, 3)} meV/bin')
sss[1].plot_disp(ax, label=f'{round(sss[1].disp*1000, 3)} meV/bin')
sss[2].plot_disp(ax, label=f'{round(sss[2].disp*1000, 3)} meV/bin')
plt.legend()
plt.title(f'dispersion = {round(np.mean(disp_bin)*1000, 3)} meV/bin')
plt.xlabel('Energy (eV)')
plt.ylabel('Elastic peak position (bins)')



# %% image2events
image = np.loadtxt('g3.dat')
plt.imshow(image)
photon_events = br.image2events(image)
pe = br.PhotonEvents(events = photon_events)
pe.plot()
pe.plot(cutoff=1)

# %% fake data
disp_bin = 0.012796  # eV/bin
disp_len = 0.051175  # eV/det.len.

pes = br.get_scan_pe_ADRESS(folderpath=fixtures, prefix='Cu_', n=18, zfill=4)
pes[0].y_max
pes[0].x_max
pes[0].plot()
plt.ylim(825, 900)
plt.title('1) raw image from experiment')
pes[0].bins = (15, 3000)
pes[0].calculate_offsets()
pes[0].fit_offsets()
pes[0].offsets_correction()
poly_coef = pes[0].offsets_popt
print('polynomial coeficients: ', poly_coef)
s0 = pes[0].calculate_spectrum(6000)
s0.calib(disp_bin)
s0.guess_elastic_peak()
s0.apply_shift(-s0.elastic_c, mode='x')
ax = s0.plot(label='experiment')

excitations = ((38, -0.45, 0.35), (90, -1.1, 0.2), (75, -1.33, 0.3), (5, -5.75, 1.5))
fake = br.fake(amp=100, c=0, w=0.14, excitations=excitations)
s1 = fake.get_spectrum()
s1.plot(ax, label='fake spectrum')

# pe = fake.get_photon_events(dispersion=disp_len)
# pe = fake.get_photon_events(dispersion=disp_len, exposure=100e3)
# pe = fake.get_photon_events(dispersion=disp_len, y_max=1500, x_max=1632, exposure=100e3)
# pe = fake.get_photon_events(dispersion=disp_len, y_max=1500, x_max=1632, exposure=650e3, poly_coef=poly_coef)
# pe = fake.get_photon_events(dispersion=disp_len, elastic_position=870, y_max=1500, x_max=1632, exposure=650e3, poly_coef=poly_coef)
pe = fake.get_photon_events(dispersion=disp_len, noise=1, elastic_position=870, y_max=1500, x_max=1632, exposure=650e3, poly_coef=poly_coef)
pe.bins = (15, 3000)
pe.calculate_offsets()
pe.fit_offsets()
print('measured coeficients: ', pe.offsets_popt)
pe.plot(show_offsets=True, show_fit=True)
pe.offsets_correction()
plt.ylim(825, 900)
plt.title('2) fake data')
s2 = pe.calculate_spectrum(6000)
s2.calib(disp_bin)

s2.guess_elastic_peak()
s2.apply_shift(-s2.elastic_c, mode='x')
s2.plot(ax, label='fake photon events')
ax.legend()


# %% energy map



# %% momentum map



# %% esrf data









# %%
