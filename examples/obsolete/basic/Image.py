#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Image example"""

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

# %% Many ways to create an br.Image object ====================================
arr = [[0, 1, 2, 0, 1, 0], 
       [0, 1, 3, 0, 2, 0], 
       [0, 1, 4, 0, 3, 0], 
       [0, 1, 5, 0, 4, 0]]

im = br.Image(data=arr)
print(im.shape)

# %% show data as an image. Equivalent to plt.imshow()
fig = br.backpack.figure()
im.imshow(colorbar=True)

# %% another way of calling imshow. Also equivalent to plt.imshow()
fig = br.backpack.figure()
im.plot(colorbar=True)

# %% any argument defined for plt.imshow works here
fig = br.backpack.figure()
im.imshow(colorbar=True, interpolation='gaussian')


# %% parameters can be set =====================================================
im.temperature = 10
print(im.temperature)


# %% save and load =============================================================
# save
im.save(r'im.txt')

# load
im2 = br.Image(r'im.txt')

# parameters are saved
print(im2.temperature)

# %% Images can be operated on =================================================
im_B = br.Image(data=[[0, 1], [1, 0]])
im_C = br.Image(data=[[1, 1], [1, 1]])

im_A = im_B + im_C
im_A = im_B - im_C
im_A = im_B * im_C
im_A = im_B / im_C  # im_C cannot contain zero pixels
im_A = im_B + 2.2
im_A = im_B - 2
im_A = im_B * 2
im_A = im_B / 2

# parameters from the first term are transferred
im_B.temperature = 10
im_A = im_B + im_C
print(im_A.temperature)
im_A = im_C + im_B
print(im_A.temperature)  # ERROR

# %% Histogram =================================================================
h = im.histogram
h = im.calculate_histogram()

# %%
fig = br.backpack.figure()
h.plot(color='black')
plt.xlabel('Intensity (arb. units)')
plt.ylabel('Number of pixels')


# %% Spectrum (vertical or horizontal integration) =============================

# %% integrate the image in the horizontal direction
s1 = im.calculate_spectrum()
s1 = im.calculate_spectrum(axis=1)

# %% plotting
fig = br.backpack.figure()
s1.plot(marker='o')
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Pixel')

# %% also, integrate the image in the vertical direction
s2 = im.calculate_spectrum(axis=0)

# %% plotting
fig = br.backpack.figure()
s2.plot(marker='o')
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Pixel')


# %% labeling rows and columns =================================================
arr = [[0, 1, 2, 0, 1, 0], 
       [0, 1, 3, 0, 2, 0], 
       [0, 1, 4, 0, 3, 0], 
       [0, 1, 5, 0, 4, 0]]

im = br.Image(data=arr)

# %% change x labels (monotonic regular step size)
im.x_centers = [0, 2, 4, 6, 8, 10]

# %% plot (note how the x axis is Labeled)
fig = br.backpack.figure()
im.plot(colorbar=True)

# %% x labels (monotonic irregular step size)
im.x_centers = [0, 2, 4, 6, 9, 13]

# %% imshow works, but pixels won't be rescaled (yields a warning)
fig = br.backpack.figure()
im.plot(colorbar=True)

# %% pcolormesh rescales the last pixels column
fig = br.backpack.figure()
im.pcolormesh(colorbar=True)

# %% x labels (non-monotonic irregular labeling)
im.x_centers = [2, 0, 4, 6, 8, 10]

# %% note how the first and second columns are inverted
fig = br.backpack.figure()
im.imshow(colorbar=True)

# %%
fig = br.backpack.figure()
im.pcolormesh(colorbar=True)

# %% x labels (repeated elements) --> not implemented yet!
im.x_centers = [0, 0, 4, 6, 8, 10]

# %% imshow will yield an error
fig = br.backpack.figure()
im.imshow(colorbar=True)  # ERROR

# %% pcolormesh will yield an error
fig = br.backpack.figure()
im.pcolormesh(colorbar=True)  # ERROR

# %% Binning ===================================================================
im2 = im.binning(nbins=(2, None))  # 2 rows, no binning
im2 = im.binning(nbins=(2))        # 2 rows, 2 columns
im2 = im.binning(nbins=2)          # 2 rows, 2 columns

im2 = im.binning((2, None))  # 2 rows, no binning
im2 = im.binning((2))        # 2 rows, 2 columns
im2 = im.binning(2)          # 2 rows, 2 columns

# %% plot
fig = br.backpack.figure()
im2.plot(colorbar=True)


# %% simulate huge fake curvature to facilitate visualization ==================
arr = [[1, 0, 0, 0, 0, 0, 1],
       [2, 1, 0, 0, 0, 1, 2], 
       [1, 2, 1, 0, 1, 2, 1], 
       [0, 1, 2, 1, 2, 1, 0], 
       [0, 0, 1, 2, 1, 0, 0],
       [0, 0, 0, 1, 0, 0, 0]]
im = br.Image(arr)

# %% plot
fig = br.backpack.figure()
im.pcolormesh(colorbar=True)

# %% plot columns ==============================================================
fig = br.figure()
_ = im.columns.plot(marker='o', vi=100)
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Pixel')

# %% roll
im.set_roll([1, 0, -1, -2, -1, 0, 1])

# %% plot
fig = br.backpack.figure()
im.plot(colorbar=True)

# %% cross-corelation
arr = [[1, 0, 0, 0, 0, 0, 1],
       [2, 1, 0, 0, 0, 1, 2], 
       [1, 2, 1, 0, 1, 2, 1], 
       [0, 1, 2, 1, 2, 1, 0], 
       [0, 0, 1, 2, 1, 0, 0],
       [0, 0, 0, 1, 0, 0, 0]]
im = br.Image(arr)

# %%
im.calculate_roll()
print(im.calculated_roll)
im.set_roll()

# %% plot
fig = br.figure()
im.plot(colorbar=True)

# %% fix curvature ====================
arr = [[1, 0, 0, 0, 0, 0, 1],
       [2, 1, 0, 0, 0, 1, 2], 
       [1, 2, 1, 0, 1, 2, 1], 
       [0, 1, 2, 1, 2, 1, 0], 
       [0, 0, 1, 2, 1, 0, 0],
       [0, 0, 0, 1, 0, 0, 0]]
im = br.Image(arr)
im2 = br.Image(arr)

# %%
p, f = im2.fix_curvature(deg=2, axis=0)

# %% plot
fig = br.figure()
im.imshow(colorbar=True)
plt.plot(im.x_centers, -f(im.x_centers), color='white', lw=10)

# %% plot
fig = br.figure()
im2.imshow(colorbar=True)


# %% image manipulations =======================================================
# intensity flooring
arr = [[10, 11, 11, 10], 
       [10, 11, 11, 10], 
       [10, 11, 11, 10], 
       [10, 11, 11, 10]]
im = br.Image(arr)

# %%
im.floor()
print(im.data)

# %% cropping
print(im.shape)
im2 = im.crop(0, None, 0, 2)
print(im2.shape)
print(im2.data)

# %% ===========================================================================
plt.close('all')
