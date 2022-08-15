# ==============================================================================
# %% BRIXS quickstart for PEAXIS beamline at HZB ================ 24/04/2022 ===
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
filepath = Path(r'../fixtures/PEAXIS/grazingmidsample_R0004.sif')

# %% importing image ============================
im = br.PEAXIS.read(filepath)

# %% plot histogram ========================
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.histogram.plot(color='black')
plt.xlabel('Intensity (arb. units)')
plt.ylabel('Number of pixels')

# %% plot image ============================
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.plot(colorbar=True)
# _ = plt.title('Original')

fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
plt.imshow(im.data)
plt.colorbar()

# %% plot spectrum ========================
s_start = im.calculate_spectrum()
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s_start.plot()

# %% parameters can be set ===============
im.temperature=10
print(im.temperature)

# %% save and load ========================
im.save(r'test.txt')

im2 = br.Image()
im2.load(r'test.txt')
print(im2.temperature)

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.fastplot(colorbar=True)
_ = plt.title('Loaded image')

# %% fake curvature ==========================
im.shifts_v
im.set_shifts([50*(x-1000)**2/(1000**2) for x in range(im.shape[1])])
im._shifts_v = np.zeros(im.shape[1])

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.fastplot(colorbar=True)
_ = plt.title('Fake curvature')

# %% binning ============================
im.binning(nbins=(None, 32))
im.bins_size

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.reduced.fastplot()
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
im.fastplot(colorbar=True)
_ = plt.title('After correction')

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



# %% TODO
floor
spectrum_h
spectrum_v
bins_size
shifts_h
shifts_v
rows
