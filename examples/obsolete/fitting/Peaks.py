#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Peaks example"""

# %% imports ===================================================================
import brixs as br
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
import copy

from brixs.file_reading import ADRESS

# %% autoreload and matplotlib backend (ignore) ================================
if br.is_notebook():
    from IPython import get_ipython
    get_ipython().run_line_magic('matplotlib', 'qt5')
    get_ipython().run_line_magic('load_ext', 'autoreload')
    get_ipython().run_line_magic('autoreload', '2')
else:
    plt.ion()
    
# %% Peak object is just a lmfit Parameters() object
peak = br.Peaks()
peak.append(amp=1, c=0, w=0.2)
peak.pretty_print()

# %%
print(peak['amp'])
print(peak['c'])
print(peak['w'])


# %% testing (Development)
peak._has_i2()
peak.check_i2()
peak._get_indexes_i1()
peak._get_indexes_i2()      # ERROR
peak._get_names_by_index()  # ERROR
peak._get_names_by_index(i1=0)
peak._get_peaks_by_index()  # ERROR
peak._get_peaks_by_index(i1=0)

peak._is_asymmetric()
peak._get_peaks_iterable()
peak._find_suitable_x()

s = peak.calculate_spectrum()

# %% plot
fig = br.figure()
_ = peak.spectrum.plot(marker='o', label='peak curve')
peak.plot(label='peak position')
plt.legend()

# %% asymmetric peak
peak = br.Peaks()
peak.append(amp=1, c=0, w1=0.2, w2=0.4)

# %%
print(peak['amp'])
print(peak['c'])
print(peak['w1'])
print(peak['w2'])
print(peak['w'])

# %% plot
fig = br.figure()
_ = peak.spectrum.plot(marker='o', label='peak curve')
peak.plot(label='peak position')
plt.legend()

# %% multiple peaks
peaks = br.Peaks()
peaks.append(amp=1, c=0, w1=0.2, w2=0.4)
peaks.append(amp=2, c=1, w=0.6)
peaks.pretty_print()

# %%
print(peaks['amp'])  # ERROR
print(peaks['amp_0'])
peaks[0].pretty_print()
peaks[1].pretty_print()
peaks[1]['amp']

# %%
peaks._get_indexes_i1()
peaks.get_area(0)

# %% plot
fig = br.figure()
_ = peaks.spectrum.plot(marker='o', label='peak curve')
peaks.plot(label='peak position')
plt.legend()

# %% save
peaks.save('peaks.p')

# %% load
peaks2 = br.Peaks('peaks.p')
peaks2.pretty_print()

# %% append
peaks.append(amp=3, c=2, w=0.2, i1=7)
peaks.pretty_print()

# %% plot
fig = br.figure()
_ = peaks.spectrum.plot(marker='o', label='peak curve')
peaks.plot(label='peak position')
plt.legend()

# %% remove
peaks.remove(7)
peaks.pretty_print()




# %% Secondary (spectrum) index i2
peaks = br.Peaks()
peaks.append(amp=1, c=0, w1=0.2, w2=0.4, i2=0)
peaks.append(amp=2, c=1, w=0.6, i2=0)

peaks.append(amp=1.1, c=0.2, w1=0.2, w2=0.4, i2=1)
peaks.append(amp=1.8, c=0.9, w=0.6, i2=1)

peaks.pretty_print()

# %%
s0 = peaks[0]
s0.pretty_print()

s1 = peaks[1]
s1.pretty_print()

# %%
br.figure()
s0.spectrum.plot()  
s0.plot()  

s1.spectrum.plot()  
s1.plot() 




# %% BRIXS =====================================================================
filepath = Path(r'./Cu_0005_d1.h5')
s = ADRESS.raw(filepath)

# %%
br.figure()
s.plot()

# %%
s.fit_peak()
s.peaks.pretty_print()

# %%
fig = br.figure()
_ = s.plot(color='black', marker='o', lw=0)
_ = s.peaks.spectrum.plot(color='red')

# %% finding peaks =============================================================
filepath = Path(r'./Cu_0017_d1.h5')
s = ADRESS.raw(filepath)

# %%
fig = br.figure()
_ = s.plot(color='black')

# %%
s.find_peaks()
s.peaks.pretty_print()

# %%
br.figure()
_ = s.plot(color='black', marker='o')
_ = s.peaks.plot()

# %% reduced prominence
s.find_peaks(prominence=2)
s.peaks.pretty_print()

# %%
fig = br.figure()
_ = s.plot(color='black', marker='o')
_ = s.peaks.plot()

# %% add peak by hand
s.peaks.append(amp=200, w=200, c=3886.7)
s.peaks.pretty_print()

# %%
fig = br.figure()
_ = s.plot(color='black')
_ = s.peaks.plot()


# %% fitting multiple peaks ====================================================
filepath = Path(r'./Cu_0017_d1.h5')
s = ADRESS.raw(filepath)

# %% find peaks
s.find_peaks()
s.peaks.append(amp=7.75, w=200, c=3543)
s.peaks.append(amp=196.15, w=30, c=3884.5)
s.peaks.pretty_print()

# %% plot
fig = br.figure()
s.plot(color='black')
s.peaks.plot()

# %% reorder peaks
# s.peaks.reorder()  # not implemented yet

# %% set bounds
s.peaks[4]['c'].max = 3895.4
s.peaks[4]['c'].min = 3873.7

# %%
fig = br.figure()
s.plot(color='black', label='data')
s.peaks.spectrum.plot(color='green', label='guess')
plt.legend()

# %% fit
s.fit_peaks()
s.peaks.pretty_print()

# %%
fig = br.figure()
s.plot(color='black', label='data')
s.peaks.spectrum.plot(color='red', label='fit')
plt.legend()

# %% plot peak contributions
fig = br.figure()
s.plot(color='black')
for peak in s.peaks._get_peaks_iterable():
    peak.spectrum.plot(lw=1)

# %% subtract peak contribution
s2 = s - s.peaks[0].calculate_spectrum(x=s.x)

fig = br.figure()
s.plot(color='black')
s2.plot(color='red')


# %% imports ===================================================================
import brixs as br
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
import copy

from brixs.file_reading import ADRESS

# %% autoreload and matplotlib backend (ignore) ================================
if br.is_notebook():
    from IPython import get_ipython
    get_ipython().run_line_magic('matplotlib', 'qt5')
    get_ipython().run_line_magic('load_ext', 'autoreload')
    get_ipython().run_line_magic('autoreload', '2')
else:
    plt.ion()

# %% fitting peaks through multiple spectra ====================================
ss = ADRESS.raw('./', 'Cu_', 17)

# %%
br.figure()
_ = ss.plot()

# %% automatically list peaks
ss.find_peaks()
ss.peaks.pretty_print()

# %% plot
br.figure()
_ = ss.plot()
for peaks in ss.peaks.iterable():
    peaks.plot()

# %% plot
br.figure()
for i, s in enumerate(ss):
    s.plot(offset=-300*i)
for i, peaks in enumerate(ss.peaks.iterable()):
    peaks.plot(offset=-300*i)

# %%

