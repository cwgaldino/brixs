#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Basic example on brixs functionality

The core brixs functionality only requires numpy and matplotlib

Advanced modules may required other python packages.
"""

# standard imports
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# default brixs import
import brixs as br

# BRIXS is based on four major objects:
im = br.Image()
pe = br.PhotonEvents()
s  = br.Spectrum()
ss = br.Spectra()

# br.Spectrum() is used for storing x and y data where y = f(x)
# br.Spectra() stores a list of br.Spectrum()
# br.PhotonEvents() is used for storing uncorelated x and y arrays
# br.Image() is used for storing 2D arrays

# for this example, we are going to use a Spectrum object, but all is also 
# valid for the other objects

# creating a y=f(x) object
x = np.linspace(0, 10, 100)
s = br.Spectrum(x=x, y=x**2)

# the x and y data are stored as attributes
print(s.x)
print(s.y)

# x and y are intrinsec attributes of the object br.Spectrum
# other intrinsec attributes of br.Spectrum are
# x, y, data, step, monotonicity, has_nan, calib, shift, offset, factor
# to get a list of all intrinsec attrs, 
s.get_core_attrs()

# read about these attributes, use help on a object
help(br.Spectrum)

# objects can also have user defined attributes,
s.T = 10  # for example, the temperature of a sample during a measurement
s.P = 1.2 # for example, the preasure of a sample during a measurement

# one can print all user defined attrs,
s.get_attrs()

# objects also have methods that act on the object
plt.figure()
s.plot()  # plots the x and y data
plt.show()

# plot() acts on s
# s.plot() is the same as plt.plot(s.x, s.y)
# use the help function to get a description of a method,
help(s.plot)

# Methods that act on the data always avoid changing the data 
# Instead, they return a new object with the modifed data
# The methods s.set_shit() returns a spectrum object with x array shifted
s2 = s.set_shift(10)      # shift the x-axis

# plot for verification
plt.figure()
s.plot()
s2.plot()
plt.show()

# one can print all methods of an object,
s.get_methods()

# methods can be pilled up to create one-line operations
# s3 is a new object similar to s where x-axis is shifted by 1, y-axis is 
# multiplied by 2.1 and data is interpolated from 2 to 3 with 200 points
s3 = s.set_shift(1).set_factor(2.1).interp(1, 9, 200)

# note that the user-attributes (s.T and s.P) are also saved to the new s3
s3.get_attrs()

# plot for verification
plt.figure()
s.plot()
s3.plot()
plt.show()

# brixs comes with support functions
# for example, br.has_duplicates returns True if list has duplicated items
br.has_duplicates([1, 2, 3, 3])

# print a full list of supporting functions 
br.get_functions()

