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
filepath = Path(r'../fixtures/PEAXIS/day3encalhk08_R0001.sif')

# %% importing image ======================================
im = br.PEAXIS.read(filepath)
print('Maximum intensity value: ', im.vmax)
print('Minimum intensity value: ', im.vmin)
print('Image size (y, x): ', im.shape)
print('\n'.join(im.header))

# %% plotting functions ====================================
# imshow uses plt.imshow()
# show data as an image
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.imshow(colorbar=True)
# _ = plt.title('plot using plt.imshow (fast most of the time)')

# imshow allows for data interpolation
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.imshow(colorbar=True, interpolation='antialiased')
# _ = plt.title('imshow with interpolation')

# plotting using imshow with matplotlib default arguments
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
plt.imshow(im.data)
plt.colorbar()

im2.x_centers
im2.y_centers

im2.x_centers
im2.y_centers

# # pcolormesh
# # most acurate plotting method
# # can be very slow sometimes
# fig = br.backpack.figure()
# br.backpack.set_window_position(2048, 232)
# im.plot(colorbar=True)
# _ = plt.title('pcolormesh (slow)')

# %% column and row labeling ===========================

# binning image just to deal with a smaller image
im.nbins = (16, 8)
im2 = im.reduced

# note how the x and y labels and tick locations are the same from the original images
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.imshow(colorbar=True)
_ = plt.title('imshow')

# %% x labels --> regular step size
im2.x_centers = [0, 2, 4, 6, 8, 10, 12, 14]

fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.reduced.imshow(colorbar=True)
_ = plt.title('plot relabeled (regular step size)')


# %% x labels --> irregular step size
im2.x_centers = [0, 2, 4, 6, 8, 10, 12, 16]
im2.y_centers

# imshow works, but pixels won't be rescaled (yields a warning)
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.reduced.imshow(colorbar=True)
_ = plt.title('plot relabeled (irregular step size)')

# plot rescales the last pixels column
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.meshplot(colorbar=True)
_ = plt.title('pcolormesh relabeled (irregular step size)')


# %% x labels --> non-monotonic irregular labeling
im2.x_centers = [2, 0, 4, 6, 8, 10, 12, 16]

# imshow
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.imshow(colorbar=True)
_ = plt.title('imshow relabeled')

# plot works even with non-monotonic labeling
# note how the first and second columns are inverted
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.plot(colorbar=True)
_ = plt.title('pcolormesh relabeled')


# x labels --> repeated elements (not implemented yet)
im2.x_centers = [0, 0, 4, 6, 8, 10, 12, 16]

# # imshow will yield an error
# fig = br.backpack.figure()
# br.backpack.set_window_position(2048, 232)
# im2.imshow(colorbar=True)
# _ = plt.title('imshow (error)')


# # plot will yield an error
# fig = br.backpack.figure()
# br.backpack.set_window_position(2048, 232)
# im2.plot(colorbar=True)
# _ = plt.title('pcolormesh (error)')

# %%
plt.close('all')
