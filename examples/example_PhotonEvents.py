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

# %% initial definitions =======================================================
filepath = Path(r'../fixtures/ADRESS/Cu_0005_d1.h5')

# %% import photonn events =====================================================
pe = br.read_ADRESS_pe(filepath)

# %% PhotonEvents basic attributes =============================================
pe.shape  # size of the detector
pe.data
pe.x
pe.y
pe.I
len(pe)  # number of photon events

# %% detector data (exclusive of ADRESS data) ==================================
pe.nd

# %% setting new parameters ====================================================
pe.temperature = 10
print(pe.temperature)

# %% plot ======================================================================
plt.figure()
pe.plot()

plt.figure()
plt.scatter(pe.x, pe.y, s=1)

# %% binning ===================================================================
pe.binning(15, 10)

plt.figure()
pe.reduced.plot(colorbar=True)

plt.figure()
pe.reduced.columns.plot()

plt.figure()
pe.reduced.rows.plot()

plt.figure()
pe.plot()
plt.vlines(pe.reduced.x+pe.bins_size[1]/2, 0, pe.shape[0], color='red')
plt.hlines(pe.reduced.y+pe.bins_size[0]/2, 0, pe.shape[1], color='red')

# %% calculate and fit shifts ==================================================
pe.binning(3000, 10)
pe.calculate_shifts()
p, f, s_fit = pe.calculated_shifts.polyfit(deg=2)

plt.figure()
pe.calculated_shifts.plot(marker='o', lw=0, color='black')
s_fit.plot(color='red')

plt.figure()
pe.plot()
s_fit.plot(offset=730, factor=-1, color='red')

# %% set shifts ================================================================
pe.set_shifts(p=p)

plt.figure()
pe.plot()

# %% spectrum ==================================================================
s = pe.calculate_spectrum(6000)

plt.figure()
s.plot()

# %% fix curvature =============================================================
pe  = br.read_ADRESS_pe(folderpath/'Cu_0005_d1.h5')
pe.binning(3000, 10)
pe.fix_curvature()

plt.figure()
pe.plot()

pe.p  # this can be used to apply the same correction in other files
pe.f
pe.shifts

# %% save and load data ========================================================
pe.save(r'test.txt')
pe.temperature = 10

pe2 = br.PhotonEvents()
pe2.load(r'test.txt')
pe2.temperature

# %%
plt.close('all')
