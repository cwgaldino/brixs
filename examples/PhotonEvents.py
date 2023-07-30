#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""PhotonEvents example"""

# %% imports ===================================================================
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import brixs as br
from brixs.file_reading import ADRESS

# %% autoreload and matplotlib backend (ignore) ================================
if br.is_notebook():
    from IPython import get_ipython
    get_ipython().run_line_magic('matplotlib', 'qt5')
    get_ipython().run_line_magic('load_ext', 'autoreload')
    get_ipython().run_line_magic('autoreload', '2')
else:
    plt.ion()

# %% photon events ============================================================
filepath = Path(r'./Cu_0005_d1.h5')
pe = ADRESS.raw(filepath, type_='pe')

# %% PhotonEvents basic attributes =============================================
print(pe.data)
print(pe.x)
print(pe.y)
print(pe.xlim)
print(pe.ylim)
print(len(pe))  # number of photon events

# %% plotting functions ========================================================
fig = br.figure()
_ = pe.plot()

# %% setting new parameters ====================================================
pe.temperature = 10
print(pe.temperature)

# %% save and load =============================================================
# save
pe.save(r'test.pe')

# load
pe2 = br.PhotonEvents(r'test.pe')

# parameters are saved
print(pe2.temperature)


# %% binning ===================================================================
im = pe.binning(100, 10)

# %%
fig = br.figure()
pos = im.plot(colorbar=True)

# %%
fig = br.figure()
_ = im.columns.plot()

# %% plot bins
fig = br.figure()
_ = pe.plot()
br.vlines(pos.x_edges, color='red')
br.hlines(pos.y_edges, color='red')

# %% calculate curvature =======================================================
filepath = Path(r'./Cu_0005_d1.h5')
pe = ADRESS.raw(filepath, type_='pe')

# %%
im = pe.binning(3000, 10)

# %%
fig = br.figure()
_ = im.plot()

# %%
popt, model = pe.calculate_shift(nbins=(3000, 10))
print(pe.calculated_shift)

s = br.Spectrum(x=pe.x, y=pe.calculated_shift)
popt, model, r2 = s.polyfit(deg=2)

# %%
fig = br.figure()
s.plot(marker='o', lw=0, color='black')
x = np.linspace(pe.xlim[0], pe.xlim[1], 100)
plt.plot(x, model(x), color='red')

# %%
fig = br.backpack.figure()
_ = pe.plot()
plt.plot(x, -model(x)+730.5, color='red')

# %% set shifts
pe.set_shift()

# %%
fig = br.figure()
pe.plot()

# %% spectrum ==================================================================
s = pe.calculate_spectrum(6000)

fig = br.figure()
s.plot()


# %% fix curvature (compact) ===================================================
filepath = Path(r'./Cu_0005_d1.h5')
pe = ADRESS.raw(filepath, type_='pe')

# %%
popt, model = pe.fix_curvature((3000, 10))

# %%
fig = br.figure()
pe.plot()

print(popt)  # this can be used to apply the same correction in other files
print(model)
