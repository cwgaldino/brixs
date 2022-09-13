# ==============================================================================
# %% EXAMPLE: photon events ================================= 22/02/2022 =======
# ==============================================================================

# %% imports ===================================================================
import brixs as br
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

%matplotlib qt5
%load_ext autoreload
%autoreload 2

# %% photonn events ============================================================
filepath = Path(r'../fixtures/ADRESS/Cu_0005_d1.h5')
pe = br.ADRESS.read_pe(filepath)

# %% PhotonEvents basic attributes =============================================
print(pe.shape)  # size of the detector
print(pe.data)
print(pe.x)
print(pe.y)
print(pe.I)
print(len(pe))  # number of photon events
print(pe.nd)  # (exclusive of ADRESS data)

# %% plotting functions ========================================================
fig = br.backpack.figure()
_ = pe.plot()

# %% setting new parameters ====================================================
pe.temperature = 10
print(pe.temperature)

# %% save and load =============================================================
# save
pe.save(r'test.txt')

# load
pe2 = br.PhotonEvents(r'test.txt')

# parameters are saved
print(pe2.temperature)


# %% binning ===================================================================
im = pe.binning(100, 10)

fig = br.backpack.figure()
_ = im.plot(colorbar=True)

fig = br.backpack.figure()
_ = im.columns.plot()

# plot bins
fig = br.backpack.figure()
_ = pe.plot()
plt.vlines(pe.reduced.x_centers+pe.bins_size[1]/2, 0, pe.shape[0], color='red')
plt.hlines(pe.reduced.y_centers+pe.bins_size[0]/2, 0, pe.shape[1], color='red')

# %% calculate, fit, and set shifts ==================================================
pe.nbins = (3000, 10)

fig = br.backpack.figure()
_ = pe.reduced.columns.plot()

pe.calculate_shift()
p, f = pe.calculated_shift.polyfit(deg=2)

fig = br.backpack.figure()
pe.calculated_shift.plot(marker='o', lw=0, color='black')
x = np.linspace(0, 1650, 100)
plt.plot(x, f(x), color='red')

fig = br.backpack.figure()
_ = pe.plot()
plt.plot(x, -f(x)+730.5, color='red')

# set shifts
pe.set_shifts(p=p)

fig = br.backpack.figure()
pe.plot()

# %% spectrum ==================================================================
s = pe.calculate_spectrum(6000)

fig = br.backpack.figure()
s.plot()


# %% fix curvature (compact) ===================================================
pe = br.ADRESS.read_pe(filepath)
pe.binning(3000, 10)
pe.fix_curvature()

fig = br.backpack.figure()
pe.plot()

print(pe.p)  # this can be used to apply the same correction in other files
print(pe.f)
