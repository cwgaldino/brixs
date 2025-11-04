#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""How to re-order Spectra"""

# %% ========================== standard imports ========================= %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# %% =========================== brixs imports =========================== %% #
import brixs as br
from brixs.model.model_functions import gaussian_fwhm

# %% ============================== settings ============================= %% #
# matplotlib (optional)
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()

# brixs (optional)
# br.settings.FIGURE_POSITION = (283, 567)
# br.get_window_position()
# %%


# %  ===================================================================== %% #
# %  ================ re-order spectra based on attribute ================ %% #
# %% ===================================================================== %% #

# creating dummy data
x = np.linspace(-4, 10, 800)
ss = br.Spectra()
for i, c in enumerate(np.arange(0, 5)):
    ss.append(br.Spectrum(x=x, y=gaussian_fwhm(x, 8-c, c, 0.5)))

# plot for verification
br.figure()
_ = ss.plot()

# let's associate each single spectrum with an Photon energy and temperature value
energies = [530, 532, 528, 526, 534]
temperatures = [20, 21, 22, 23, 24]
for i, _s in enumerate(ss):
    _s.E = energies[i]
    _s.T = temperatures[i]

# we can check if attribute has been stored inside the Spectrum object
print(ss[0].get_attrs())
print(ss[0].E)

# note that these attributes were stored as Spectrum attributes
# but we can also have them as Spectra attributes
ss.get_attrs()  # empty
ss.E = [s.E for s in ss]
ss.T = [s.T for s in ss]

ss.get_attrs()  # now, E and T are a list of energies and temperatures
print(ss.E)
print(ss.T)

# plot for verification
br.figure()
_ = ss.plot(label=ss.E)
br.leg()

# create energy map 
im = ss.stack_spectra_as_columns()
im.x_centers = ss.E

# plot for verification 
# THIS YIELDS AN ERROR because the energies are out of order and it is not 
# possible to create an image 
print(im.x_centers)
# br.figure()
# im.plot()  # this yields an error

# we can, however, re-order the columns of the image. And that makes it possible
# to plot the image 
im2 = im.fix_x_monotonicity()
print(im2.x_centers)

# plot for verification 
br.figure()
im2.plot()

# we can also reorder the Spectra directly
# use ss.reorder_by_attr to order ss based on the any Spectra attribute values
# help(ss.reorder_by_attr)  # use this to print the help for this method
ss2 = ss.reorder_by_attr(attr='E')

# plot for verification
fig, axes = br.subplots(1, 2)
_ = ss2.plot(ax=axes[0], label=ss2.E)
br.leg(ax=axes[0])

im3 = ss2.stack_spectra_as_columns()
im3.x_centers = ss2.E

im3.plot(ax=axes[1])