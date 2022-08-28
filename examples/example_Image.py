# ==============================================================================
# %% EXAMPLE: Image ============================================= 24/08/2022 ===
# ==============================================================================

# %% Initial imports ===========================================================
import brixs as br
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


# %% initial definitions =======================================================
# filepath = Path(r'../fixtures/PEAXIS/elastic_0001.sif')
# filepath = Path(r'../fixtures/PEAXIS/elastic_0002.sif')
filepath = Path(r'../fixtures/PEAXIS/sample_0001.sif')


# %% importing image ===========================================================
im = br.PEAXIS.read(filepath)


# %% some attributes ===========================================================
print('Maximum intensity value: ', im.vmax)
print('Minimum intensity value: ', im.vmin)
print('Image size (y, x): ', im.shape)
print('\n'.join(im.header))  # exclusive of PEAXIS dataset


# %% plotting functions ========================================================

# show data as an image. Equivalent to plt.imshow()
fig = br.backpack.figure()
im.imshow(colorbar=True)

# another way of calling imshow. Also equivalent to plt.imshow()
fig = br.backpack.figure()
im.plot(colorbar=True)

# any argument defined for plt.imshow works here
fig = br.backpack.figure()
im.imshow(colorbar=True, interpolation='antialiased')

# # pcolormesh is the most accurate plotting method, but can be very slow
fig = br.backpack.figure()
im.pcolormesh(colorbar=True)


# %% parameters can be set =====================================================
im.temperature = 10
print(im.temperature)


# %% save and load =============================================================
# save
im.save(r'test.txt')

# load
im2 = br.Image(r'test.txt')

# parameters are saved
print(im2.temperature)


# %% Histogram =================================================================
h = im.histogram
h = im.calculate_histogram()

fig = br.backpack.figure()
h.plot(color='black')
plt.xlabel('Intensity (arb. units)')
plt.ylabel('Number of pixels')


# %% Spectrum (vertical or horizontal integration) =============================

# many ways to integrate the image in the horizontal direction
s1 = im.calculate_spectrum()
s1 = im.calculate_spectrum(axis=1)
s1 = im.calculate_spectrum(axis='row')
s1 = im.calculate_spectrum(axis='horizontal')
s1 = im.calculate_spectrum(axis='x')
s1 = im.spectrum
s1 = im.spectrum_h

# plotting
fig = br.backpack.figure()
s1.plot()
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Pixel')

# also, many ways to integrate the image in the vertical direction
s2 = im.calculate_spectrum(axis=0)
s2 = im.calculate_spectrum(axis='column')
s2 = im.calculate_spectrum(axis='vertical')
s2 = im.calculate_spectrum(axis='y')
s2 = im.spectrum_v

# plotting
fig = br.backpack.figure()
s2.plot()
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Pixel')



# %% Binning ===================================================================

# many ways to define the binning (defining number of bins)
im.nbins = (16, 8)  # 16 rows and 8 columns
im.nbins = (None, 8)  # 2048 rows and 8 columns
im.nbins = (8)      # same as (8, 8)
im.nbins = 8        # same as (8, 8)
im.binning(nbins=(16, 8))
im.binning(nbins=(8))
im.binning(nbins=8)

# many ways to define the binning (defining bin size)
im.bins_size = (128, 256)
im.bins_size = (None, 256)  # row size 1, column size 256
im.bins_size = (256)
im.bins_size = 256
im.binning(bins_size=(128, 256))
im.binning(bins_size=(256))
im.binning(bins_size=256)

# binned data is another image
im2 = im.reduced

# plot
fig = br.backpack.figure()
im2.plot(colorbar=True)

# %% labeling rows and columns =================================================

# change x labels (monotonic regular step size)
im2.x_centers = [0, 2, 4, 6, 8, 10, 12, 14]

# plot (note how the x axis is Labeled)
fig = br.backpack.figure()
im.reduced.imshow(colorbar=True)

# %% x labels (monotonic irregular step size)
im2.x_centers = [0, 2, 4, 6, 8, 10, 12, 16]

# imshow works, but pixels won't be rescaled (yields a warning)
fig = br.backpack.figure()
im.reduced.imshow(colorbar=True)

# pcolormesh rescales the last pixels column
fig = br.backpack.figure()
im2.pcolormesh(colorbar=True)

# one can also change labels by changing the labels of the bin edges
im2.x_edges = [0, 1, 3, 5, 7, 9, 11, 13, 20]
print(im2.y_edges)

fig = br.backpack.figure()
im2.pcolormesh(colorbar=True)

# %% x labels (non-monotonic irregular labeling)
im2.x_centers = [2, 0, 4, 6, 8, 10, 12, 16]

# note how the first and second columns are inverted
fig = br.backpack.figure()
im2.imshow(colorbar=True)

fig = br.backpack.figure()
im2.pcolormesh(colorbar=True)

# x labels (repeated elements) --> not implemented yet!
im2.x_centers = [0, 0, 4, 6, 8, 10, 12, 16]

# # imshow will yield an error
fig = br.backpack.figure()
im2.imshow(colorbar=True)

# # pcolormesh will yield an error
fig = br.backpack.figure()
im2.pcolormesh(colorbar=True)




# %% simulate huge fake curvature to facilitate visualization ==================
print(im.shifts_v)
im.set_shift([50*(x-1000)**2/(1000**2) for x in range(im.shape[1])])
print(im.shifts_v)
im._shifts_v = np.zeros(im.shape[1])  # erase shifts

# plot
fig = br.backpack.figure()
im.imshow(colorbar=True)

# %% binning ===================================================================
im.binning(nbins=(None, 32))
print(im.bins_size)

# plot
fig = br.backpack.figure()
im.reduced.imshow()

# %% plot columns ==============================================================
fig = br.backpack.figure()
_, _ = im.reduced.columns.plot()
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Pixel')


# %% fix curvature ====================
im.fix_curvature(deg=2, axis=0)

# plot
fig = br.backpack.figure()
im.imshow(colorbar=True)

# plot columns
fig = br.backpack.figure()
_, _ = im.reduced.columns.plot()

# %% spectrum ===========================
s = im.calculate_spectrum()

# plot
fig = br.backpack.figure()
s.plot()


# %% ===========================================================================
plt.close('all')
