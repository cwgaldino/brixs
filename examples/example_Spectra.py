# ==============================================================================
# %% EXAMPLE: Spectra ====================================== 26/04/2022 =======
# ==============================================================================

# %% imports ===================================================================
import brixs as br
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import copy

%matplotlib qt5
%load_ext autoreload
%autoreload 2

# %% Create br.Spectra object ==================================================
filepath = Path(r'../fixtures/ADRESS/Cu_0017_d1.h5')
s0 = br.ADRESS.read(filepath)

filepath = Path(r'../fixtures/ADRESS/Cu_0017_d2.h5')
s1 = br.ADRESS.read(filepath)

filepath = Path(r'../fixtures/ADRESS/Cu_0017_d3.h5')
s2 = br.ADRESS.read(filepath)

ss = br.Spectra([s0, s1, s2])
print(ss)

# %% plot ======================================================================
fig = br.backpack.figure()
_ = ss.plot()

# %% plot with legend ==========================================================
fig = br.backpack.figure()
_ = ss.plot()
plt.legend(['s0', 's1', 's2'])

# %% plot with vertical and horizontal increment ===============================
fig = br.backpack.figure()
_ = ss.plot(vi=10, hi=10)

# %% append and remove =========================================================
ss.remove(2)
ss.append(s2)

# %% spectra in ss are the same in s1, s2, and s3 ==============================
print(ss[2])
print(s2)

# %% modifiers =================================================================
print(ss.offset)
print(ss.calib)
print(ss.factor)
print(ss.shift)
print(ss.shift_roll)
print(ss.shift_interp)

# modifer is applied to all spectra (by default, a relative modification is applied)
ss.offset = 50
print(ss.offset)

ss.set_offset(20, type='relative')
print(ss.offset)

# a list can be used
ss.offset = [10, 20, 30]
print(ss.offset)


# %% save and load =============================================================
ss.save(dirpath='test/spectra', prefix='spectrum_', suffix='.dat')

# load from a folder
ss1 = br.Spectra(dirpath='test/spectra')

# load from a folder filtering filenames
ss1 = br.Spectra()
ss.load(dirpath='test/spectra', string='spectrum')



# %% check methods and attributes ==============================================
print(ss.step)
print(ss.length)
print(ss.x)
print(ss.monotonicity)

ss.check_length()
print(ss.length)

ss.check_monotonicity()
print(ss.monotonicity)

ss.check_same_x()
print(ss.length)
print(ss.x)

ss.check_step_x()
print(ss.length)
print(ss.step)
print(ss.x)

# a shift will reset ss.x, because the x is not the same for all spectra
ss1 = copy.deepcopy(ss)
ss1.shift = [10, 20, 30]
print(ss1[0].x)
print(ss1[1].x)
print(ss1[2].x)
print(ss1.step)
print(ss1.length)
print(ss1.x)
print(ss1.monotonicity)

# %% Area ======================================================================
print(ss.area)

# %% Interpolation =============================================================
ss1 = copy.deepcopy(ss)
ss1.shift = [10, 20, 30]
print(ss1[0].x)
print(ss1[1].x)
print(ss1[2].x)
ss1.interp()
print(ss1[0].x)
print(ss1[1].x)
print(ss1[2].x)

ss2 = copy.deepcopy(ss)
ss2.shift = [10, 20, 30]
print(ss2[0].x)
print(ss2[1].x)
print(ss2[2].x)
ss2.interp(start=30, stop=5000, num=3000)
print(ss2[0].x)
print(ss2[1].x)
print(ss2[2].x)


# %% find and defining peaks ====================================================
ss2 = copy.deepcopy(ss)

# automatically list peaks
ss1.find_peaks()

# automatically list peaks
for s in ss1:
    s.find_peaks()

# check
print(ss1.peaks)

# plot
br.backpack.figure()
_, offsets, shifts = ss1.plot(vi=110)
_ = ss1.peaks.plot(offset=offsets)

# adjusting find peaks parameters
ss1.find_peaks(prominence=1)

# split peak 1 into two peaks
for s in ss1:
    s.peaks.split(1)

# manually add peaks
for s in ss1:
    s.peaks.append({'amp': 200, 'c':3000, 'fwhm':10})

# plot
br.backpack.figure()
_, offsets, shifts = ss1.plot(vi=110)
_ = ss1.peaks.plot(offset=offsets)



# %% peak fitting ==============================================================
ss1 = copy.deepcopy(ss)
ss1.find_peaks()

# add satelite peak (split peak 1 into two peaks)
for s in ss1:
    s.peaks.append({'amp': 200, 'c':s.peaks[0]['c']-10, 'fwhm':s.peaks[0]['fwhm']})

# check fitting boundaries
spectrum_index = 0
peak_index     = 0
ss1[spectrum_index].peaks[peak_index].bounds

# manually change fitting boundaries
spectrum_index = 0
peak_index     = 0
parameter      = 'amp'
ss1[spectrum_index].peaks[peak_index].bounds[parameter] = (0, 1000)

# automatically change fitting boundaries
for s in ss1:
    for peak in s.peaks:
        peak.set_bounds(amp=300, fwhm=200, c=1, type='percentage')

# fit
ss1.fit_peaks()

# plot fit
br.backpack.figure()
_, offsets, shifts = ss1.plot(vi=110)
_, offsets, shifts = ss1.fit.plot(offset=offsets, color='red')

# plot guess
br.backpack.figure()
_, offsets, shifts = ss1.plot(vi=110)
_, offsets, shifts = ss1.guess.plot(offset=offsets, color='green')

# plot residue
br.backpack.figure()
_, offsets, shifts = ss1.plot(vi=110)
_, offsets, shifts = ss1.fit.plot(offset=offsets, color='red')
_, offsets, shifts = ss1.residue.plot(offset=offsets, color='blue')

# plot individual peak contributions
br.backpack.figure()
_ = ss1[0].plot(color='black')
_ = ss1[0].fit.plot(color='red')
for peak in ss1[0].fit.peaks:
    _ = peak.spectrum.plot(color='blue')

# fit R2
print(ss1.fit.R2)

# peak paremeters
print(ss1.fit.peaks)
print(ss1.fit.peaks['amp'][0])
print(ss1.fit.peaks['c'][-1])
print(ss1.fit.peaks.get_errors('amp')[0])
print(ss1.fit.peaks.get_errors('c')[-1])


# %% fitting spectra with modifiers ============================================
ss1 = copy.deepcopy(ss)
ss1.find_peaks()
for s in ss1:
    s.peaks.append({'amp': 200, 'c':s.peaks[0]['c']-10, 'fwhm':s.peaks[0]['fwhm']})
    for peak in s.peaks:
        peak.set_bounds(amp=300, fwhm=200, c=1, type='percentage')
ss1.fit_peaks()

# plot fit
br.backpack.figure()
_, offsets, shifts = ss1.plot(vi=110)
_, offsets, shifts = ss1.fit.plot(offset=offsets, color='red')

ss2 = copy.deepcopy(ss)
ss2.shift = [10, 20, 30]
ss2.find_peaks()
for s in ss2:
    s.peaks.append({'amp': 200, 'c':s.peaks[0]['c']-10, 'fwhm':s.peaks[0]['fwhm']})
    for peak in s.peaks:
        peak.set_bounds(amp=300, fwhm=200, c=1, type='percentage')
ss2.fit_peaks()

# plot fit
br.backpack.figure()
_, offsets, shifts = ss2.plot(vi=110)
_, offsets, shifts = ss2.fit.plot(offset=offsets, color='red')

# shift info is passed to the fit
print(ss1.shift)
print(ss1.fit.shift)
print(ss2.shift)
print(ss2.fit.shift)

# modifier info is not passed to the fitted peaks (this is okay for now)
print(ss2.peaks.shift)
print(ss2.guess.peaks.shift)
print(ss2.fit.peaks.shift)

# %% Calculate shifts ==========================================================

# via max
ss1 = copy.deepcopy(ss)
ss1.calculate_shift(mode='max')
ss1.set_shift()

fig = br.backpack.figure()
_, offsets, shifts = ss1.plot()


# via fitted peaks
ss1 = copy.deepcopy(ss)
ss1.find_peaks()
for s in ss1:
    s.peaks.append({'amp': 200, 'c':s.peaks[0]['c']-10, 'fwhm':s.peaks[0]['fwhm']})
    for peak in s.peaks:
        peak.set_bounds(amp=300, fwhm=200, c=1, type='percentage')
ss1.fit_peaks()
ss1.calculate_shift(mode='fitted peaks', peak=-1)
ss1.set_shift()

fig = br.backpack.figure()
_, offsets, shifts = ss1.plot()

# parameters
print(ss.calculated_shift.y)
print(ss.calculated_shift.mode)
print(ss.calculated_shift.ref_value)

# x axis is NOT necessarily preserved!
ss1.check_same_x()
print(ss1.peaks)
ss1.peaks['c'][-1]

# via cross-correlation
ss1 = copy.deepcopy(ss)
ss1.calculate_shift(mode='cc')
ss1.set_shift()

fig = br.backpack.figure()
_, offsets, shifts = ss1.plot()

# x axis is preserved for cc!
ss1.check_same_x()




# %% map =======================================================================
ss1 = copy.deepcopy(ss)
map = ss1.map

fig = br.backpack.figure()
map.plot()
