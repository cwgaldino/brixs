#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""A PhotonEvents object `pe` stores x and y arrays (just like Spectrum object),
but in `pe`, arrays do not correlate. PhotonEvents can be used to store x and y
coordinates of photon detection

>>> x  = np.random.rand(10)
>>> y  = np.random.rand(10)
>>> pe = br.PhotonEvents(x, y)

One can easily get the x and y values from `pe`

>>> pe.x
>>> pe.y

Methods operate on these x and y values


 ##############################
>>> print(pe.y[1])
>>> pe.set_shift(2)
>>> print(pe.y[1])
 ##############################

help for each method is available through python's built-in help function
 ##############################
>>> help(pe.set_offset)
 ##############################

Attrs can be added to `pe` on the go

>>> pe.new_attr = 10

New methods can also be added on the go (see advanced examples)

>>> br.PhotonEvents.new_method = new_method

See below for a more complete description.
"""

# %% ---------------------------- imports --------------------------------- %% #
import brixs as br
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# %% ---------------------- quick Example --------------------------------- %% #
x  = np.random.rand(100)
y  = np.random.rand(100)
pe = br.PhotonEvents(x, y)

plt.figure()
_ = pe.plot()

# %% ---------------------- Example 1: Initialization --------------------- %% #
# empty PhotonEvents list
pe = br.PhotonEvents()

# PhotonEvents can be initialized like this
x  = np.random.rand(100)
y  = np.random.rand(100)
pe = br.PhotonEvents(x, y)
pe = br.PhotonEvents(x=x, y=y)

# x and y data
print(pe.x)
print(pe.y)
print(pe.data)

# get number of PhotonEvents (length of the list)
print(len(pe))

# get (x, y) pairs
print(pe[0])

# slicing also works
pe2 = pe[0:2]

# items (x, y pairs) can be deleted
del pe[0]

# as for right now, items cannot be changed (this might change in the future)
# pe[2] = (9, 10)  # error

# attributes can be associated with `pe` object
pe.temperature = 10
pe.composition = [9, 3, 5, 6]
pe.string      = 'test'
pe.dictionary  = {'a': 10, 'b':5}
print(pe.temperature)
print(pe.composition)
print(pe.string)
print(pe.dictionary)

# get a list of all defined attrs and remove all attrs
print(pe.get_attrs())
# pe.remove_attrs()

# one can save and load a PhotonEvents' list
pe.save('examples/pe_0.dat')
pe2 = br.PhotonEvents('examples/pe_0.dat')
# any xy-type files can be loaded 
# Comments must be marked with `#` and columns must be separated by `,` (comma) 
# The first column must be the x array

# PhotonEvents' list can also be loaded via load function
pe2 = br.PhotonEvents()
pe2.load('examples/pe_0.dat', comments='#')
# load function have args that allow for more loading options
# for instance the comments tag can be set

# attrs are saved and loaded with spectrum
# Note that, dictionaries and more complex attrs are not saved (this might change in the future)
print(pe.temperature)
# print(pe.dictionary)  # ERROR

# attrs can be copied between spectrum objects
pe2.copy_attrs_from(pe)

# %% ----------------------- Example 2: plotting -------------------------- %% #
x  = np.random.rand(100)
y  = np.random.rand(100)
pe = br.PhotonEvents(x, y)

br.figure()   # same as plt.figure(), but with additional features
pe.plot()     # same as plt.plot(s.x, s.y)

# one can define valid x and y limits for plotting
pe.xlim = (0., 0.8)
pe.ylim = (0.2, 2)

br.figure()
pe.plot() 

# note, opening a new figure with br.figure() allows for sending the x and y 
# values where the mouse cursor is to the the clipboard by left and right 
# mouse button click. Pressing the wheel copies the image. Double left click
# prints the figure coordinates (from 0 to 1).

# %% ----------------- Example 3: binning (2D histograming) --------------- %% #
x  = np.random.rand(100)
y  = np.random.rand(100)
pe = br.PhotonEvents(x, y)

# x and y limits are used for defining the start and stop values for histograming
# if they are not defined, limits will be set to the min and max values of the data 
pe.xlim = (0, 1)
pe.ylim = (0, 1)

# binning (or histograming) creates a br.Image type (see Image example file)
im = pe.binning(100, 10)  # 100 bins vertical, 10 bins horizontal

# plot bins
fig = br.figure()
_   = pe.plot()
br.vlines(im.x_edges, color='red', ls='-', lw=0.5)
br.hlines(im.y_edges, color='red', ls='-', lw=0.5)

# plot image (2D histogram)
fig = br.figure()
pos = im.plot(colorbar=True)

# %% ---------------------- Example 4: Spectrum --------------------------- %% #
x  = np.random.rand(4000)
y  = np.random.normal(0.5, 0.15, 4000)
pe = br.PhotonEvents(x, y)
pe.xlim = (0, 1)
pe.ylim = (0, 1)

# create_spectrum() function will bin `pe` in one direction
# axis = 0 (vertical), axis = 1 (horizontal) following numpy convention 
# Default is axis = 0
sy = pe.calculate_spectrum(nbins=100, axis=0)
sx = pe.calculate_spectrum(nbins=100, axis=1)
sx.switch()  # swap x and y coordinates

# plot
fig, axes = plt.subplots(2, 2)
_   = pe.plot(ax=axes[0][0])
_   = sy.plot(ax=axes[1][0])
_   = sx.plot(ax=axes[0][1])

# %% -------------------- Example 5: base modifiers ----------------------- %% #
x  = np.random.rand(100)
y  = np.random.rand(100)
pe = br.PhotonEvents(x, y)

# shifts point (adds value)
# axis = 0, will add value to the y coordinates and axis = 1 to the x coordinates
# default is axis = 0
pe.set_shift(1, axis=0)

# value can be a array of same length of points. Each point will be shifted by 
# the respective amount given in the array
values = np.linspace(0, 1, 100)
pe.set_shift(values, axis=0)

# instead of values, one can pass a 1D array of polynomial coefficients 
# from highest degree to the constant term. For example, the following line of 
# code will shift y coordinates using the polynomial function 
# y = 4*x**2 -5*x + 9
pe.set_shift(p=[4, -5, 9], axis=0)
# similarly, one can shift x coordinates based on their y counterpart by switching
# the argument axis to 1
# x = 4*y**2 -5*y + 9
pe.set_shift(p=[4, -5, 9], axis=1)

# finally, one can pass a function y=fx(x) for axis=0 or x=fy(y) for axis=1
# pe.set_shift(f=fx, axis=0)
# pe.set_shift(f=fy, axis=1)

# %% ----------------------- Example 6: modifiers ------------------------- %% #
# swap x and y coordinates
# pe.transpose()

# remove (x, y) pairs that fall outside of range
# pe.crop(x_start, x_stop, y_start, y_stop)

# remove (x, y) pairs that falls outside the mask
# mask must be a list with rectangles coordinates (x_start, x_stop, y_start, y_stop)
# mask can be something like, mask = (0, 12, 3, 12) or mask = [(0, 12, 3, 12), (1, 4, 15, 18)]
# pe.clip(mask)

# %% ----------------------- Example 6: copying --------------------------- %% #
x  = np.random.rand(100)
y  = np.random.rand(100)
pe = br.PhotonEvents(x, y)

# coping syntax 1
pe2 = br.PhotonEvents()
pe2.copy(pe)

# coping syntax 2
pe2 = pe.copy()

# %% ---------------------- Example 7: calculate -------------------------- %% #

calculate_shift


fix_curvature


# %% calculate curvature =======================================================
filepath = Path(r'./Cu_0005_d1.h5')
pe = ADRESS.raw(filepath, type_='pe')

# 
im = pe.binning(3000, 10)

# 
fig = br.figure()
_ = im.plot()

# 
popt, model = pe.calculate_shift(nbins=(3000, 10))
print(pe.calculated_shift)

s = br.Spectrum(x=pe.x, y=pe.calculated_shift)
popt, model, r2 = s.polyfit(deg=2)

# 
fig = br.figure()
s.plot(marker='o', lw=0, color='black')
x = np.linspace(pe.xlim[0], pe.xlim[1], 100)
plt.plot(x, model(x), color='red')

# 
fig = br.backpack.figure()
_ = pe.plot()
plt.plot(x, -model(x)+730.5, color='red')

# set shifts
pe.set_shift()

# %%
fig = br.figure()
pe.plot()

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

