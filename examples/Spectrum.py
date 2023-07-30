#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Spectrum example"""

# %% imports ===================================================================
import brixs as br
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
import copy

# %% autoreload and matplotlib backend (ignore) ================================
if br.is_notebook():
    from IPython import get_ipython
    get_ipython().run_line_magic('matplotlib', 'qt5')
    get_ipython().run_line_magic('load_ext', 'autoreload')
    get_ipython().run_line_magic('autoreload', '2')
else:
    plt.ion()

# %% example ===================================================================
x = np.linspace(-5, 5, 100)
y = np.sin(x)

s = br.Spectrum(x, y)

plt.figure()
s.plot()
plt.show()

# %% Initialization ============================================================
# s = br.Spectrum(x, y)
# s = br.Spectrum(y)
# s = br.Spectrum(filepath)

s = br.Spectrum(x=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], y=[0, 1, 2, 3, 4, 4, 3, 2, 1, 0])

# filepath must point to a two column txt or csv file
# s = br.Spectrum('<filepath>')

# %% data reading functions are defined for some beamlines =====================
filepath = Path(r'./Cu_0005_d1.h5')
s = br.ADRESS.read(filepath)

# some attributes
print('x axis: ', s.x)
print('y axis: ', s.y)
print('Integrated area: ', s.area)

# exclusive of ADRESS dataset
print(s.th)
print(s.slit)

# plotting 
fig = br.figure()  # br.figure is the same as plt.figure, but has additional functionality
s.plot()

# %% save and load =============================================================
s = br.Spectrum(x=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], y=[0, 1, 2, 3, 4, 4, 3, 2, 1, 0])

# parameters can be set
s.temperature = 10
print(s.temperature)

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

# parameters from the first term are transferred
s_B.temperature = 10
s_A = s_B + s_C
print(s_A.temperature)
s_A = s_C + s_B
print(s_A.temperature)  # ERROR

# %% check x axis ==============================================================
s_A = br.Spectrum(x=[0, 1, 2, 3, 4, 5, 6], y=[1, 1, 1, 1, 1, 1, 1])
print(s_A.step)
s_A.check_step()
print(s_A.step)
s_A.check_monotonicity()
print(s_A.monotonicity)

s_A = br.Spectrum(x=[0, 2, 4, 6, 8, 10, 12], y=[1, 1, 1, 1, 1, 1, 1])
s_A.check_step()
print(s_A.step)
s_A.check_monotonicity()
print(s_A.monotonicity)

s_A = br.Spectrum(x=[0, 1.1, 2, 2.3, 7, 10, 12], y=[1, 1, 1, 1, 1, 1, 1])
s_A.check_step()  # ERROR
s_A.check_monotonicity()
print(s_A.monotonicity)

s_A = br.Spectrum(x=[6, 5, 4, 3, 2, 1, 0], y=[1, 1, 1, 1, 1, 1, 1])
s_A.check_step()
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
s = br.Spectrum(x=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], y=[0, 1, 2, 3, 4, 4, 3, 2, 1, 0])

# offset
fig = br.figure()
s.plot()
s.set_offset(20)
s.plot()

# modifiers can be set via the attribute
s.offset = 30
s.plot()

# %% There are two ways of applying a shift in the x axis =====================
s = br.Spectrum(x=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], y=[0, 1, 2, 3, 4, 4, 3, 2, 1, 0])

# hard shift. Y axis is fully preserved, x is modified
print('x axis before: ', s.x)
fig = br.figure()
s.plot(label='before')
s.set_shift(10)
s.plot(label='after')
plt.legend()
print('x axis after: ', s.x)

# roll shift. both axis are preserved. Only integer shifts allowed.
# Data is "rolled" in the array
# The edges of the data became meaningless
fig = br.backpack.figure()
s.plot(label='before')
s.set_roll(3)
s.plot(label='after')
plt.legend()


# %% manipulation ==============================================================
s = br.Spectrum(x=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], y=[0, 1, 2, 3, 4, 4, 3, 2, 1, 0])
s2 = copy.deepcopy(s)

fig = br.backpack.figure()
s.plot(label='initial')
s2.interp(start=0, stop=10, num=100)
s2.plot(label='interp')

s2.crop(start=2, stop=None)
s2.plot(label='crop')

s2 = s2 + 50
s2.plot(label='offset')
s2.floor()
s2.plot(label='floor')

s2.flip()
s2.plot(label='flip')

s2.normalize(50, [-10, -9])
s2.plot(label='nomalize')

s3 = s.copy([[0, 2], [5, 7], [9, None]])
s3.plot(label='extracted', marker='o')

plt.legend()

