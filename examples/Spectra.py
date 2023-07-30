#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Spectra example"""

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



# %% initial example
x0 = np.linspace(-5, 5, 100)
y0 = np.sin(x0)
s0 = br.Spectrum(x0, y0)

x1 = np.linspace(-5, 5, 100)
y1 = 0.6*np.sin(x1)
s1 = br.Spectrum(x1, y1)

x2 = np.linspace(-5, 5, 100)
y2 = 0.8*np.cos(x2)
s2 = br.Spectrum(x2, y2)

ss = br.Spectra(s0, s1, s2)

# %% plot
plt.figure()
_ = ss.plot()
plt.legend(['s0', 's1', 's2'])

# %% two ways of modify a spectrum
s2.shift = 1
ss[2].factor = 0.5

# %% in fact, spectra in ss are the same in s1, s2, and s3
print(ss[1])
print(s1)

# %% plot
plt.figure()
_ = ss.plot()
plt.legend(['s0', 's1', 's2'])

# %% reorder
ss.reorder(0, 2)

# %% plot
plt.figure()
_ = ss.plot()
plt.legend(['s2', 's1', 's0'])

# %% flip order
ss.flip_order()

# %% plot
plt.figure()
_ = ss.plot()
plt.legend(['s0', 's1', 's2'])

# %% plot with vertical and horizontal increment
fig = br.backpack.figure()
_ = ss.plot(vi=100, hi=100)

# %% append and remove
ss.remove(2)
ss.append(s2)

# %% modifier is applied to all spectra
ss.offset = 5
plt.figure()
_ = ss.plot()
plt.legend(['s1', 's2', 's3'])

# %% another away
ss.set_offset(-5)
plt.figure()
_ = ss.plot()
plt.legend(['s1', 's2', 's3'])

# %% a list can be used
ss.offset = [1, 2, 3]
plt.figure()
_ = ss.plot()
plt.legend(['s1', 's2', 's3'])

# %% save and load
ss.save(folderpath='./', prefix='spectrum_', suffix='.dat')

# %% load from a folder filtering filenames
ss1 = br.Spectra()
ss1.load(folderpath='./', string='.dat')



# %% check methods and attributes
print(ss.step)
print(ss.length)
print(ss.x)
print(ss.monotonicity)

# %%
ss.check_length()
print(ss.length)

# %%
ss.check_monotonicity()
print(ss.monotonicity)

# %%
ss.check_same_x()
print(ss.length)
print(ss.x)

# %%
ss.check_step()
print(ss.length)
print(ss.step)
print(ss.x)

# %% a shift will reset ss.x, because the x is not the same for all spectra
ss1 = ss.copy()
ss1.shift = [10, 20, 30]
print(ss1[0].x)
print(ss1[1].x)
print(ss1[2].x)
print(ss1.step)
print(ss1.length)
print(ss1.x)
print(ss1.monotonicity)

# %% Area
print(ss.area)



# %% Interpolation
ss1 = ss.copy()
ss1.shift = [10, 20, 30]
print(ss1[0].x)
print(ss1[1].x)
print(ss1[2].x)
ss1.interp()
print(ss1[0].x)
print(ss1[1].x)
print(ss1[2].x)

ss2 = ss.copy()
ss2.shift = [10, 20, 30]
print(ss2[0].x)
print(ss2[1].x)
print(ss2[2].x)
ss2.interp(start=30, stop=5000, num=3000)
print(ss2[0].x)
print(ss2[1].x)
print(ss2[2].x)

# %% map =======================================================================
x0 = np.linspace(-5, 5, 100)
y0 = np.sin(x0)
s0 = br.Spectrum(x0, y0)

x1 = np.linspace(-5, 5, 100)
y1 = 0.6*np.sin(x1)
s1 = br.Spectrum(x1, y1)

x2 = np.linspace(-5, 5, 100)
y2 = 0.8*np.cos(x2)
s2 = br.Spectrum(x2, y2)

ss = br.Spectra(s0, s1, s2)

fig = br.figure()
ss.map.plot()





