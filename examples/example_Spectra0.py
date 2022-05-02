# ==============================================================================
# %% EXAMPLE: Spectra ====================================== 26/04/2022 =======
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
filepath1 = Path(r'D:\galdino\Documents\CuSb2O6\analysis\rixs\exp\2022_02_PSI\data\CuSb2O6_Feb_2022\RIXS\Cu_0033_d1.h5')
filepath2 = Path(r'D:\galdino\Documents\CuSb2O6\analysis\rixs\exp\2022_02_PSI\data\CuSb2O6_Feb_2022\RIXS\Cu_0033_d2.h5')
filepath3 = Path(r'D:\galdino\Documents\CuSb2O6\analysis\rixs\exp\2022_02_PSI\data\CuSb2O6_Feb_2022\RIXS\Cu_0033_d3.h5')


# %% import Spectra ============================================================
s0 = br.read_ADRESS(filepath1)
s1 = br.read_ADRESS(filepath2)
s2 = br.read_ADRESS(filepath3)

ss = br.Spectra([s0, s1, s2])
print(ss)

# %% plot ======================================================================
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot()

# %% plot with legend ==========================================================
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot()
plt.legend(['s0', 's1', 's2'])

# %% append and remove =========================================================
# s2 was removed
ss.remove(2)
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot()

# append back s2
ss.append(s2)
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot()

# the spectrum is acctually copied to the dataset
# ss[2] is a different object than s2
s2.offset = 50

# offset to s2 does not change ss[2]
ss.append(s2)
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot()

# apply offsets to ss[2]
ss[2].offset = 50
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot()

# spectrum can be append as a list
ss.append([s0, s1, s2])
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot()

# %% save and load =============================================================
ss.save(dirpath='test/', suffix='test')
ss.load(dirpath='test/')
ss1 = br.Spectra(dirpath='test/')

fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss1.plot()

# %% other methods =============================================================
ss = br.Spectra([s0, s1, s2])

fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot()

ss2 = ss.extract((3000, 5000))
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss2.plot()

ss3 = ss.extract(((100, 1000), (3000, 5000), (5500, None)))
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss3.plot(marker='o', lw=0, ms=3)

ss.interp()
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot()

ss.crop(3000, 5000)
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot()

# %% Calculate shifts ==========================================================
ss = br.Spectra([s0, s1, s2])
for s in ss:
    s.find_peaks()
    s.peaks.append({'amp':7.75, 'fwhm':200, 'c':3543})
    s.peaks.split(1)
    s.fit_peaks()

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
_ = ss.plot(vi=80)
_ = ss.plot_peaks(vi=80, color='black')

# calculate shifts
ss.calculate_shift(mode='fitted peaks', peak=-1)
ss.calculated_shift.y
ss.calculated_shift.mode
ss.calculated_shift.ref_value

# set shift
ss.set_shift()

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
ss.plot()

# x axis is NOT preserved!
ss.check_same_x()

# %% Calculate shifts and change shift mode ====================================
ss = br.Spectra([s0, s1, s2])
for s in ss:
    s.find_peaks()
    s.peaks.append({'amp':7.75, 'fwhm':200, 'c':3543})
    s.peaks.split(1)
    s.fit_peaks()
ss.calculate_shift(mode='fitted peaks', peak=-1)

# set shift
# ss.set_shift(mode='roll')  # will raise an error
ss.check_step_x()
values = np.round(ss.calculated_shift.y/ss.step)
ss.set_shift(values, mode='roll')

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
ss.plot()

# x axis is preserved!
ss.check_same_x()




s0.calculate_area()

# %% modifiers =================================================================
s.offset
s.calib
s.factor
s.shift
s.shift_roll
s.shift_interp

# offset
s.set_offset(20)
s.plot()

# modifiers are absolute quantities
s.set_offset(0)
s.plot()

# modifiers can be set via the attribute
s.offset = 30
s.plot()

# for relative modifications use assignment operators
s.offset -= 30
s.plot()

# %% hard shift modifier =======================================================
s = br.read_ADRESS(filepath)

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()

# hard ('hard', 'x') shift
# note how the whole spectra moves to the right
# y shape is fully preserved
# x axis is not preserved
s.set_shift(500, mode='hard')
s.plot()

# %% soft shift modifier =======================================================
s = br.read_ADRESS(filepath)

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()

# soft ('soft', 'y') shift
# y shape is interpolated in a different position
# x axis is fully preserved
# in this plot, data note that the left edge of the data is meaningless
s.set_shift(500, mode='soft')
s.plot()

# %% roll shift modifier =======================================================
s = br.read_ADRESS(filepath)

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()

# roll ('roll', 'r') shift
# y shape is fully preserved
# x axis is fully preserved
# the data points that fall outside the x range are put on the other side of the
# spectrum. Therefore, the left edge of this plot is meaningless.
s.set_shift(500, mode='roll')
s.plot()
