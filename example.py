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
oct_2021 = Path(r'D:\CuSb2O6\CuSb2O6_oct_2021\RIXS')
prefix = 'Cu'
dispersion = 0.012796


# %% check curvature correction

# processed by the beamline
filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
s0.guess_elastic_peak(prominence=25)
s0.set_x(s0.get_x() - s0.elastic_c)
ax2 = s0.plot()
plt.text(50, 0, 'beamline fwhm: ' + str(round(s0.elastic_w, 1)))

# single correction (best binning)
filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
pe.bins = (60, 6000)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(prominence=25)

s.set_x(s.get_x() - s.elastic_c)
s.plot(ax2, vertical_increment=50)
plt.text(50, 50, str(pe.bins) + ' fwhm: ' + str(round(s.elastic_w, 1)))

plt.xlim(-200, 200)

# %% guess binning
filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
print(pe.bins)
pe.bins = 'guess'
print(pe.bins)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.plot(show_offsets=True, show_fit=True)
plt.ylim(840, 890)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(prominence=25)
print(s.elastic_w)

# %% guess binning and increase its size (20 %)
filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
pe.bins = 'guess'
pe.bins = pe.bins + pe.bins*0.2
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(prominence=25)
print(s.elastic_w)

# %% guess binning and increase its size (30 %)
filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
pe.bins = 'guess'
pe.bins = pe.bins + pe.bins*0.3
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(prominence=25)
print(s.elastic_w)

# %% guess + best binning (it should be (7, 6000) with fwhm = 6.902)
filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
pe.bins = 'guess'
print(pe.bins)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
best_bins, best, fwhm = pe.best_bins(prominence=25, deg=2, spectrum_bins=6000)
pe.bins = best_bins
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.plot(show_offsets=True, show_fit=True)
plt.ylim(840, 890)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(prominence=25)
print(s.elastic_w)

# %% best binning (it should be (7, 6000) with fwhm = 6.902)
filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
best_bins, best, fwhm = pe.best_bins(prominence=25, deg=2, spectrum_bins=6000)
pe.bins = best_bins
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.plot(show_offsets=True, show_fit=True)
plt.ylim(840, 890)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
s.guess_elastic_peak(prominence=25)
print(s.elastic_w)


# %% best binning for spectrum calculation
filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
pe.bins = (60, 6000)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
best_bins, best, fwhm = pe.best_bins_for_spectrum(ref_bin=6000, ref_prominence=25, y_bins=None)
print(best_bins)

for y_bin in [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000]:
    s = pe.calculate_spectrum(bins=best_bins, mode='lenght')
    ax = s.plot(ms=4, marker='o', label = best_bins)

    s = pe.calculate_spectrum(bins=y_bin, mode='lenght')
    s.y = s.y*y_bin/3000
    if y_bin > 7000:
        s.plot(ax, label=(1, y_bin))
    else:
        s.plot(ax, ms=4, marker='o', label=(1, y_bin))

    plt.legend()
    plt.xlim(878, 890)
# plt.close('all')

# s = pe.calculate_spectrum(bins=5000, mode='lenght')
# ax = s.plot(ms=4, marker='o')
# s.guess_elastic_peak(points=1, prominence=25)
# plt.plot(s.x, s.elastic_func(s.x))
# print(s.elastic_w)

# %% Adding 3 ccds
filenames = ['Cu_0018_d1.h5', 'Cu_0018_d2.h5', 'Cu_0018_d3.h5']

ss = br.Spectra()
for filename in filenames:
    _, s, _, _ = br.read_ADRESS(oct_2021/filename)
    print(type(s))
    ss.append(s)

    type(s)
isinstance(s, br.Spectrum)

# %%
for i in []:
    print(i)

print(popt[2])

filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
pe.bins = (7, 6000)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=3000)
smooth, popt, err = s.guess_elastic_peak()
s.plot(ax)
print(popt[2])


filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
pe.bins = (7, 6000)
pe.calculate_offsets(ref=int(pe.bins[0]/2))
pe.fit_offsets(deg=2)
pe.offsets_correction()
s = pe.calculate_spectrum(bins=6000)
smooth, popt, err = s.guess_elastic_peak()
s.plot(ax)
print(popt[2])


# %% photon events with bad events
filepath = oct_2021/'Cu_0018_d1.h5'
pe, s0, bad, nd = br.read_ADRESS(filepath)
ax1 = pe.plot()
bad.plot(ax1, mfc='red', pointsize=2)

# binning x too big
pe.bins = (256, 1000)
pe.plot(show_bins=(True, False))
pe.plot_columns(columns=[0, 10, -1])

# binning x for curvature correction
pe.bins = (20, 1000)
pe.plot(show_bins=(True, False))
pe.plot_columns(columns=[0, 10, -1], vertical_increment=5)

# %% autobinning
x_guess = 9
y_guess = 100
max_total_counter = 1000
max_x_counter = 5
x_step = 1
y_step = 50

pe.bins = (x_guess, y_guess)
pe.calculate_offsets()
diff = np.diff(pe.offsets)
diff_sum = sum(diff)

x_guess -= 1
x_counter = 0
total_counter = 0
while 0 in diff and total_counter < max_total_counter:
    print(pe._bins)
    pe.bins = (pe._bins[0]+x_step, pe._bins[1])
    pe.calculate_offsets()
    diff = np.diff(pe.offsets)

    if x_counter < max_x_counter-1:
        if diff_sum == sum(diff):
            x_counter += 1
    else:
        pe._bins = (x_guess, pe._bins[1]+y_step)
        x_counter = 0

    diff_sum = sum(diff)
    total_counter += 1
if total_counter >= max_total_counter:
    pe.bins = (x_guess, y_guess)
    pe.calculate_offsets()
    raise RuntimeError('cannot find solution.')


pe.bins
# %%

pe.bins = (10, 700)
pe.calculate_offsets()
pe.fit_offsets(deg=2)
ax = pe.plot_offsets()
pe.plot_offsets_fit(ax)
pe.plot(show_offsets=True, show_offsets_fit=True)

pe.bins = (12, 1000)
pe.calculate_offsets()
pe.fit_offsets(deg=2)
ax = pe.plot_offsets()
pe.plot_offsets_fit(ax)
pe.plot(show_offsets=True, show_offsets_fit=True)

# %%

# curvature correction
pe.calculate_offsets()
pe.plot_offsets()
pe.plot(show_bins=(True, False), show_offsets=True)

# %%
pe.bins = (12, 1000)
pe.calculate_offsets()
pe.fit_offsets(deg=2)
ax = pe.plot_offsets()
pe.plot_offsets_fit(ax)
pe.plot(show_offsets=True, show_offsets_fit=True)

# pe.bins = (13, 1000)
# pe.calculate_offsets()
# ax = pe.plot_offsets()

pe.bins = (20, 1000)
pe.calculate_offsets()
pe.fit_offsets(deg=2)
pe.plot_offsets(ax)
pe.plot_offsets_fit(ax)
pe.plot(show_offsets=True, show_offsets_fit=True)

pe.bins = (25, 1000)
pe.calculate_offsets()
pe.fit_offsets(deg=2)
pe.plot_offsets(ax)
pe.plot_offsets_fit(ax)
pe.plot(show_offsets=True, show_offsets_fit=True)


pe.bins = (50, 1000)
pe.calculate_offsets()
pe.fit_offsets(deg=2)
pe.plot_offsets(ax)
pe.plot_offsets_fit(ax)
pe.plot(show_offsets=True, show_offsets_fit=True)

# plt.close('all')
# %%
nd['ArraySizeY']

nd['ArraySizeX']
nd.keys()
pe.plot()
print(pe.x_max)
print(pe.y_max)
pe.x_max = None
pe.y_max = None
print(pe.x_max)
print(pe.y_max)

pe.bins =

bad.plot(mfc='red')
print(bad.x_max)
print(bad.y_max)
bad.x_max = None
bad.y_max = None
print(bad.x_max)
print(bad.y_max)
# %%



# %%

# %%

# %%



image = np.loadtxt('g.dat')
image = image.transpose()
# plt.imshow(image)
events = br.image2events(image)
# plt.scatter(events[:, 0], events[:, 1], s=0.01)
br = importlib.reload(br)
e = br.PhotonEvents(events)
e.events
# del e.events
max(e.events[:, 0])
e.y_max = 1280
print(e.x_max)
# e.save2file('pe.dat')

e.plot(cutoff=1)
# %%
print(e.bins)
print(e.bins_size)
e.bins = 10
print(e.bins)
print(e.bins_size)
e.bins = (10, 11)
print(e.bins)
print(e.bins_size)


e.bins_size = 1
print(e.bins)
print(e.bins_size)

e._bins = (111, 111)
e.bins = (111, 111)
e.bins

e.hist
e.hist[0]
# e.hist = 10
# del e.hist

e.x_edges
e.y_edges
e.x_centers
e.y_centers

e.y_max = None
e.y_max
e.bins=(256, 1280)
e.bins_size

e.y_max = 2000
e.y_max
e.bins=(256, 1280)
e.bins_size

e.y_max = 1280
e.y_max
e.bins=(10, 1280)
e.bins_size

e.x_max
e.y_max

e.plot(cutoff=1)
e.calculate_offsets()
e.fit_offsets(deg=2)
e.offsets_correction()
e.calculate_spectrum(1000)
max(e.events[:, 1])
e.offsets
e.x_max(2)
e.y_max

max(e.events[:, 1])
e.y_max
x = np.arange(0, e.bins[1])
np.vstack((x, sum(e.hist))).transpose()

cross_correlation= np.correlate(e.hist[8], e.hist[0], mode='same')
e.y_centers[np.argmax(cross_correlation)]

plt.figure()
plt.plot(e.y_centers, e.hist[0])
plt.plot(e.y_centers, e.hist[-1])

plt.figure()
plt.scatter(e.x_centers, e.offsets)
x = np.linspace(e.x_centers[0], e.x_centers[-1], 100)
plt.plot(x, e.offsets_func(x))
len(sum(e.hist))

e.bins[1]
x = np.arange(0, e.bins[1])
x
data = np.vstack((x, sum(e.hist))).transpose()
plt.plot(data[:, 0], data[:, 1])
br.Spectrum(data=  np.vstack((x, sum(e.hist))).transpose())
e.offsets

256*(0.00036)
256*(0.00036+0.0104)

0.001999378204345703
0.3396110534667969
0.0004571732133626938
0.000865524634718895

# %%
# %%
a = [[11, 22], [33, 44]]
for i, (x, y) in enumerate(a):
    print(x)
a = []*2
a
0.022007465362548828
2.7662899494171143
0.00039771944284439087
0.010400297120213509


import time
x = np.arange(0, 10000)
y = np.zeros((10000, 2))
y[:, 0] += x
y[:, 1] += -x
ranges = [[0, 5000], [9000, 9999]]

# x = np.arange(0, 10000)
# y= -x
# ranges = [[0, 5000]]

start = time.time()
choose_range = choose(x,ranges)
np.compress(choose_range, y, axis=0)
print(time.time() - start)

start = time.time()
choose_range = choose(x, ranges)
np.extract(choose_range, y)
print(time.time() - start)

start = time.time()
extract(x, y, ranges)
print(time.time() - start)

start = time.time()
x2, y2 = extract2(x, y, ranges)
print(time.time() - start)

x2[:, 0]
y2[:,0]

len(y2[0])

# one y
0.0059964656829833984
0.008994340896606445
0.06681180000305176

# two y
0.006978273391723633
0.009994268417358398
0.12664270401000977

# ten y
0.007992982864379883
0.008994102478027344
0.2948317527770996


# one y two ranges
0.005998849868774414
0.005993843078613281
0.08894968032836914

# two y 2 range
0.00599360466003418
0.011993646621704102
0.014786958694458008
0.0005266666412353516


x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
y = np.zeros((10, 2))
y[:, 0] = x**2
y[:, 0] = x**3
ranges = ((0, 3), (7.5, 9))
x_sliced, y_sliced = extract2(x, y, ranges)
print(x_sliced)
print(y_sliced)
a = [[True, False], [False, False], [False, False]]
np.sum(a, axis=1)


def choose(x, ranges):
    try:
        choose_range = [None]*len(ranges)
        for i, (x_init, x_final) in enumerate(ranges):
            choose_range[i] = np.logical_and(x>=x_init, x<=x_final)
        print(choose_range)
        choose_range = [True if x == 1 else False for x in np.sum(choose_range, axis=0)]
    except TypeError:
        x_init, x_final = ranges
        choose_range = np.logical_and(x>=x_init, x<=x_final)
    return choose_range
def extract(x, y, ranges):
    x = np.array(x)
    y = np.array(y)

    choose_range = choose(x, ranges)
    temp = np.compress(choose_range, np.c_[y.transpose(), x], axis=0)
    if len(temp[0]) > 2:
        return temp[:, -1], temp[:, :-1].transpose()
    else:
        return temp[:, -1], temp[:, 0]

# %%
# %%

# br.__dir__()
# %%
# %%
# %%
import brixs3 as br3
br3 = importlib.reload(br3)

image = np.loadtxt('g.dat')
photon_events = br3.image2events(image.transpose())

e = br3.photon_events(photon_events, x_max=image.shape[0], y_max=image.shape[1])
# e.plot(cutoff=1)
# e.plot(cutoff=2)
e.set_binning(bins=(256, 1280))
# ax1 = e.plot_columns(columns=0)
# e.plot_columns(ax1, columns=255)
spectrum = e.calculate_spectrum()

e.calculate_offsets()
# e.plot(cutoff=1, show_offsets=True)
e.fit_offsets(deg=2)
# e.plot(cutoff=1, show_offsets=True, show_offsets_fit=True)
e.offsets_correction()
# e.plot(cutoff=1,)
spectrum2 = e.calculate_spectrum()

# comparisson before and after
ax = spectrum.plot()
spectrum2.plot(ax)

# useful
spectrum2.get_x()
spectrum2.get_y()

# %%
