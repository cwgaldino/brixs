# ==============================================================================
# %% EXAMPLE: Spectrum ====================================== 26/08/2022 =======
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

# %% Many ways to create a br.Spectrum object ==================================
arr = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]]
s = br.Spectrum(data=arr)

# also works if the array is inverted (transposed)
arr = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]]
arr2 = np.transpose(arr)
s = br.Spectrum(data=arr2)

s = br.Spectrum(x=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], y=[0, 1, 2, 3, 4, 4, 3, 2, 1, 0])

# filepath must point to a two column txt or csv file
# s = br.Spectrum('<filepath>')


# %% data reading functions are defined for some beamlines =====================
filepath = Path(r'../fixtures/ADRESS/Cu_0005_d1.h5')
s = br.ADRESS.read(filepath)


# %% some attributes ===========================================================
print('x axis: ', s.x)
print('y axis: ', s.y)
print('Integrated area: ', s.area)

# exclusive of ADRESS dataset
print('experiment parameters: ' + str(s.nd.keys()))
print(s.nd['SampleTheta'])
print(s.nd['ExitSlit'])


# %% plotting functions ========================================================
fig = br.backpack.figure()
s.plot()

# %% parameters can be set =====================================================
s.temperature = 10
print(s.temperature)


# %% save and load =============================================================
# save
s.save(r'test.txt')

# load
s2 = br.Spectrum(r'test.txt')

# parameters are saved
print(s2.temperature)


# %% Spectrum can be operated on (y axis) ======================================
s_B = br.Spectrum(y=[0.1, 1, 2, 3, 4, 5])
s_C = br.Spectrum(y=[5, 4, 3, 2, 1, 0.1])

s_A = s_B + s_C
s_A = s_B - s_C
s_A = s_B * s_C
s_A = s_B / s_C  # s_C cannot contain zeros in the y axis
s_A = s_B + 2.2
s_A = s_B - 2
s_A = s_B * 2
s_A = s_B / 2

# parameters from the first term are transfered
s_B.temperature = 10
s_A = s_B + s_C
print(s_A.temperature)
s_A = s_C + s_B
print(s_A.temperature)  # ERROR

# %% check x axis ==============================================================
s_A = br.Spectrum(x=[0, 1, 2, 3, 4, 5, 6], y=[1, 1, 1, 1, 1, 1, 1])
print(s_A.step)
s_A.check_step_x()
print(s_A.step)
s_A.check_monotonicity()
print(s_A.monotonicity)

s_A = br.Spectrum(x=[0, 2, 4, 6, 8, 10, 12], y=[1, 1, 1, 1, 1, 1, 1])
s_A.check_step_x()
print(s_A.step)
s_A.check_monotonicity()
print(s_A.monotonicity)

s_A = br.Spectrum(x=[0, 1.1, 2, 2.3, 7, 10, 12], y=[1, 1, 1, 1, 1, 1, 1])
s_A.check_step_x()  # ERROR
s_A.check_monotonicity()
print(s_A.monotonicity)

s_A = br.Spectrum(x=[6, 5, 4, 3, 2, 1, 0], y=[1, 1, 1, 1, 1, 1, 1])
s_A.check_step_x()
print(s_A.step)
s_A.check_monotonicity()
print(s_A.monotonicity)

s_A = br.Spectrum(x=[0, 12, 5, 1, 8, 6, 2], y=[0, 1, 2, 3, 4, 5, 6])
s_A.check_monotonicity()  # ERROR
s_A.fix_monotonicity()
s_A.check_monotonicity()
print(s_A.monotonicity)
print(s_A.x)
print(s_A.y)


# %% modifiers =================================================================
print(s.offset)       # additive factor, y axis
print(s.factor)       # multiplicative factor, y axis
print(s.calib)        # multiplicative factor, x axis
print(s.shift)        # additive factor, x axis (see below)
print(s.shift_roll)   # additive factor, x axis (see below)
print(s.shift_interp) # additive factor, x axis (see below)

# offset
fig = br.backpack.figure()
s.plot()
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

# or use argument type
s.set_offset(10, type='relative')
s.plot()


# %% There are many ways of applying a shift in the x axis =====================

# hard shift. Y axis is fully preserved, x is modified
print('x axis before: ', s.x)
fig = br.backpack.figure()
s.plot(label='before')
s.set_shift(500, mode='hard')
s.plot(label='after')
plt.legend()
print('x axis after: ', s.x)

# soft shift. X axis is fully preserved, y axis is interpolated in a different position
fig = br.backpack.figure()
s.plot(label='before')
s.set_shift(500, mode='soft')
s.plot(label='after')
plt.legend()

# roll shift. both axis are preserved. Only integer shifts allowed.
# Data is "rolled" in the array
# The edges of the data became meaningless
fig = br.backpack.figure()
s.plot(label='before')
s.set_shift(500, mode='roll')
s.plot(label='after')
plt.legend()


# %% manipulation ==============================================================
area = s.calculate_area()

s2 = copy.deepcopy(s)

fig = br.backpack.figure()
s.plot(label='initial')
s2.interp(start=2000, stop=4000, num=2000)
s2.plot(label='interp')

s2.crop(start=2500, stop=None)
s2.plot(label='crop')

s2 = s2 + 50
s2.plot(label='offset')
s2.floor()
s2.plot(label='floor')

s2.flip()
s2.plot(label='flip')

s2.normalize(50, [-2920.68, -2912.48])
s2.plot(label='nomalize')

s3 = s.extract([[0, 1000], [2000, 3000], [5000, None]])
s3.plot(label='extracted', marker='o', lw=0, ms=1)

s3.zero()
s3.plot(label='zero')

plt.legend()






# %% fit peak ==================================================================
fig = br.backpack.figure()
_ = s.plot(color='black', marker='o', lw=0)
s.fit_peak()
print(s.peaks)
print(s.fit.peaks)
print('fit R2: ' + str(s.R2))
_ = s.fit.plot(color='red')

s.save('test.txt')

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
