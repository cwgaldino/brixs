# ==============================================================================
# %% EXAMPLE: Spectrum ====================================== 22/02/2022 =======
# ==============================================================================

# %% imports ===================================================================
import brixs as br
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

%matplotlib qt5
%load_ext autoreload
%autoreload 2

# %% initial definitions =======================================================
filepath = Path(r'../fixtures/ADRESS/Cu_0005_d1.h5')
filepath = Path(r'D:\galdino\Documents\CuSb2O6\analysis\rixs\exp\2022_02_PSI\data\CuSb2O6_Feb_2022\RIXS\Cu_0033_d1.h5')

# %% import photonn events =====================================================
s = br.read_ADRESS(filepath)
print(s.x)
print(s.y)

# %% plot  =====================================================================
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()

# %% finding peaks =============================================================
s.find_peaks()
print(s.peaks)
_ = s.peaks.plot()

# reduced prominence
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()
s.find_peaks(prominence=2)
_ = s.peaks.plot()

# add peak by hand
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()
s.find_peaks()
s.peaks.append({'amp':7.75, 'fwhm':200, 'c':3543})
_ = s.peaks.plot()

# split peak
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()
s.find_peaks()
s.peaks.split(0)
_ = s.peaks.plot()



# %% fitting peaks =============================================================
s = br.read_ADRESS(filepath)

# find peaks
s.find_peaks()
s.peaks.append({'amp':7.75, 'fwhm':200, 'c':3543})
s.peaks.split(1)

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot(color='black')
_ = s.peaks.plot(color='green')

# fit
s.fit_peaks()
s.fit.plot(color='red')

# fitted peaks parameters
print(s.fit.peaks)

# if main data is modified via a modifier, peaks and fit are modified too
s.shift  = 100
s.offset = 100
s.calib  = 100
s.factor = 100

fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot(color='black')
_ = s.peaks.plot(color='green')
s.fit.plot(color='red')

# %% residue and initial guess =================================================
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot(color='black')
s.fit.plot(color='red', label='fit')
s.guess.plot(color='blue', label='guess')
s.residue.plot(color='green', label='residue')
plt.title(f'R2: {s.R2}')

# %% if data is negative bound must be adjusted ================================
s = br.read_ADRESS(filepath)
s.offset = -300

fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot(color='black')
s.find_peaks()
_ = s.peaks.plot(color='green')
s.fit_peaks(amp_bounds=(3, -3), verbose=True)
s.fit.plot(color='red')
