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
im = br.read_PEAXIS(filepath)
print('Maximum intensity value: ', im.vmax)
print('Minimum intensity value: ', im.vmin)
print('Image size (y, x): ', im.shape)
print('\n'.join(im.header))

# %% plotting functions ====================================
# fastplot uses plt.pcolorfast() for quick plotting
# should be the prefered method
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.fastplot(colorbar=True)
_ = plt.title('fastplot (fast)')

# imshow uses plt.imshow()
# show data as an image
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.imshow(colorbar=True)
_ = plt.title('imshow (fast most of the time)')

# imshow allows for data interpolation
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.imshow(colorbar=True, interpolation='antialiased')
_ = plt.title('imshow with interpolation')

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

# fastplot and imshow yield the same x and y labels
# ticks are located roughly in the center of the pixels
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.fastplot(colorbar=True)
_ = plt.title('fastplot')

# fastplot and imshow yield the same x and y labels
# ticks are located roughly in the center of the pixels
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.imshow(colorbar=True)
_ = plt.title('imshow')

# plot calculates better the ticks position
# note how the x and y labels and tick locations are the same from the original images
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.plot(colorbar=True)
_ = plt.title('pcolormesh')

# %% x labels --> regular step size
im2.x_centers = [0, 2, 4, 6, 8, 10, 12, 14]

fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.reduced.fastplot(colorbar=True)
_ = plt.title('fastplot relabeled (regular step size)')

fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.imshow(colorbar=True)
_ = plt.title('imshow relabeled (regular step size)')

fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.plot(colorbar=True)
_ = plt.title('pcolormesh relabeled (regular step size)')


# %% x labels --> irregular step size
im2.x_centers = [0, 2, 4, 6, 8, 10, 12, 16]

# fastplot works, but pixels won't be rescaled
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.reduced.fastplot(colorbar=True)
_ = plt.title('fastplot relabeled (irregular step size)')

# imshow works, but pixels won't be rescaled
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.imshow(colorbar=True)
_ = plt.title('imshow relabeled (irregular step size)')

# plot rescales the last pixels column
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im2.plot(colorbar=True)
_ = plt.title('pcolormesh relabeled (irregular step size)')


# %% x labels --> non-monotonic irregular labeling
im2.x_centers = [2, 0, 4, 6, 8, 10, 12, 16]

# fastplot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
im.reduced.fastplot(colorbar=True)
_ = plt.title('fastplot relabeled')

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

# # fastplot will yield an error
# fig = br.backpack.figure()
# br.backpack.set_window_position(2048, 232)
# im.reduced.fastplot(colorbar=True)
# _ = plt.title('fastplot')

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
