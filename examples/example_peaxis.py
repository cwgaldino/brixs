# ==============================================================================
# %% BRIXS quickstart for PEAXIS beamline at HZB ================ 22/03/2022 ===
# ==============================================================================

# %% Initial imports ===========================================================
import brixs as br
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

%matplotlib qt5
%load_ext autoreload
%autoreload 2

# %% initial definitions ==========================
peaxis_example1 = Path('<path-to-brixs>/fixtures/PEAXIS/calib')
peaxis_example2 = Path('<path-to-brixs>/fixtures/PEAXIS/map')

peaxis_example1 = Path(r'D:\galdino\github\brixs/fixtures/PEAXIS/calib')
peaxis_example2 = Path(r'D:\galdino\github\brixs/fixtures/PEAXIS/map')

# %%
image, header = br.ReadAndor(peaxis_example2/'RIXSmshk-004_R0001.sif')
# print(header)
plt.close('all')

# %% initializing image object
im = br.Image(image)
im.vmax
im.vmin
im.shape

# %% plot image ============================
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.plot(colorbar=True)
_ = plt.title('Original')

# %% plot spectrum ========================
s_start = im.calculate_spectrum()
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s_start.plot()

# %% plot histogram ========================
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.histogram.plot()

# %% parameters can be set ===============
im.temperature=10
print(im.temperature)

# %% save and load ========================
im.save(r'D:\galdino\Documents\test.txt')

im2 = br.Image()
im2.load(r'D:\galdino\Documents\test.txt')
print(im2.temperature)

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.plot(colorbar=True)
_ = plt.title('Loaded image')

# %% fake curvature ==========================
im.shifts_v
im.set_shifts([50*(x-1000)**2/(1000**2) for x in range(im.shape[1])])
im._shifts_v = np.zeros(im.shape[1])

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.plot(colorbar=True)
_ = plt.title('Fake curvature')

# %% binning ============================
im.binning(nbins=(None, 32))
im.bins_size

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.reduced.plot()
_ = plt.title(f'number of bins {im.nbins}')

# %% plot columns ======================
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.reduced.columns.plot()

# %% fix curvature ====================
im.fix_curvature()

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.plot(colorbar=True)
_ = plt.title('After correction 1')

# %% fix curvature again ==============
im.binning(nbins=(None, 8))
im.fix_curvature(deg=1)

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.plot(colorbar=True)
_ = plt.title('After correction 2')

# plot columns
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.reduced.columns.plot()

# %% spectrum ===========================
s_final = im.calculate_spectrum()

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s_start.plot()
s_final.plot()
