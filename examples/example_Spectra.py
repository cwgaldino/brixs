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
ss.save(dirpath='test/', prefix='spectrum_', suffix='.dat')

ss1 = br.Spectra()
ss.load(dirpath='test/')

ss1 = br.Spectra(dirpath='test/')


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
ss.shift = [10, 20, 30]
print(ss[0].x)
print(ss[1].x)
print(ss[2].x)
print(ss.step)
print(ss.length)
print(ss.x)
print(ss.monotonicity)

# %% Interpolation =============================================================
print(ss.area)

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


# %% fit =======================================================================
ss1 = copy.deepcopy(ss)
for s in ss1:
    s.find_peaks()

# plot
br.backpack.figure()
_, offsets, shifts = ss1.plot(vi=110)
_ = ss1.plot_peaks(offset=offsets)

# adjust peaks
for s in ss1:
    s.peaks.split(0)
    for peak in s.peaks:
        peak.set_bounds(amp=3, fwhm=2, c=0.2, type='multiplicative')

# fit
ss1.fit_peaks()

# plot
br.backpack.figure()
_, offsets, shifts = ss1.plot(vi=110)
_, offsets, shifts = ss1.fit.plot(offset=offsets, color='red')

# fit parameters
ss1_fitted  = ss1.fit
print(ss1_fitted.R2)
print(ss1_fitted.pcov)
ss1_residue = ss1.residue
ss1_guess   = ss1.guess

fitted_peaks = ss1_fitted.peaks
print(fitted_peaks[0]['amp'])
print(fitted_peaks[-1]['c'])

fitted_errors = ss1_fitted.errors
print(fitted_errors[0]['amp'])
print(fitted_errors[-1]['c'])


# %% Calculate shifts ==========================================================

# via max
ss1 = copy.deepcopy(ss)
ss1.calculate_shift(mode='max')
ss1.set_shift()

fig = br.backpack.figure()
_, offsets, shifts = ss1.plot()


# via fitted peaks
ss1 = copy.deepcopy(ss)
for s in ss1:
    s.find_peaks()
    s.peaks.split(0)
    for peak in s.peaks:
        peak.set_bounds(amp=3, fwhm=2, c=0.2, type='multiplicative')
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
