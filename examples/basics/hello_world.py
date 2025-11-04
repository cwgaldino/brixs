#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Basic example on brixs functionality

The core of brixs only requires numpy and matplotlib

Advanced modules may required other python packages.
"""

# standard imports
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# default brixs import
import brixs as br

# one can print all functions defined inside br by using the following line
br.get_functions()

# BRIXS is based on four major objects:
im = br.Image()
pe = br.PhotonEvents()
s  = br.Spectrum()
ss = br.Spectra()

# Once a BRIXS object has been created, one can use the following methods to 
# print all methods (functions) and attrs related to that object (in this case s).
s  = br.Spectrum()
s.get_methods()
s.get_attrs()

# The description of methods can be accessed via the python help() function. See example below,
s  = br.Spectrum()
help(s.set_shift)

# creating a y=f(x) object 
x = np.linspace(0, 10, 1000)
s = br.Spectrum(x=x, y=x**2)

# example of initial processing
s = s.set_shift(10)      # shift the x-axis
s = s.set_factor(2.1)    # apply a multiplicative factor
s = s.interp(2, 3, 200)  # interpolate data

# display data
br.figure()                        # open new figure
s.plot(label=f"T={s.T}, P={s.P}")  # plot data
br.leg()                           # legend
