#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""A Spectra object is a list that stores Spectrum objects. 

Just like a Spectrum object, a Spectra type can have Attrs and methods.
Methods operate on the list of Spectrum types. New methods and attrs can be 
defined on the fly.

"""

# %% ---------------------------- imports --------------------------------- %% #
import brixs as br
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# %% ---------------------- quick Example --------------------------------- %% #
x0 = np.linspace(-5, 5, 100)
y0 = np.sin(x0)
s0 = br.Spectrum(x0, y0)

x1 = np.linspace(-5, 5, 100)
y1 = 0.6*np.sin(x1)
s1 = br.Spectrum(x1, y1)

x2 = np.linspace(-5, 5, 100)
y2 = 0.8*np.cos(x2)
s2 = br.Spectrum(x2, y2)

ss = br.Spectra(s0, s1, s2)

plt.figure()
_ = ss.plot()
plt.legend(['s0', 's1', 's2'])

# %% ---------------------- Example 1: Initialization --------------------- %% #
# empty spectra with length n
ss = br.Spectra(n=10)
ss = br.Spectra(10)

# empty spectra
ss = br.Spectra()

# Spectra can be initialized many ways (run s0, s1, and s2 from quick example above)
ss = br.Spectra(s0, s1, s2)
ss = br.Spectra([s0, s1, s2])
ss = br.Spectra(data=[s0, s1, s2])

# get number of Spectrum item (length of the list)
print(len(ss))

# Spectra is a list and Spectrum items can be accessed via indexing and slicing
print(ss[0])
print(ss[1])
print(ss[2])
print(ss[0:2])

# items can be changed, appended, and removed
ss.append(s1)
ss.remove(-1)

ss.append(s1)
del ss[-1]

ss[-1] = s2

# attributes can be associated with spectrum object
ss.temperature = 10
ss.composition = [9, 3, 5, 6]
ss.string      = 'test'
ss.dictionary  = {'a': 10, 'b':5}
print(ss.temperature)
print(ss.composition)
print(ss.string)
print(ss.dictionary)

# get a list of all defined attrs and remove all attrs
print(ss.get_attrs())
# ss.remove_attrs()

# one can save and load an spectrum
ss.save('examples/')
# ss = br.Spectrum('examples/')  # ERROR
# using this syntax, Spectra can only be loaded if they are the only files in the folder

# spectra can also be loaded via load function
ss = br.Spectra()
br.filelist('examples/', string='spectrum')

ss.load('examples/', string='spectrum', comments='#')
# load function have args that allow for more loading options
# for instance, string selects only files where the filename contains a certain string

# multiple spectra can be saved in one file
ss.interp()  # interpolate data so every spectrum has the same x-axis
ss.save_all_single_file('examples/spectra.dat')

# filepath must must point to a xy-type file, where comments must be marked with
# `#` and columns must be separated by `,` (comma). The first column must be the
# x array
# ss = br.Spectra(filepath='<filepath>')
# ss = br.Spectra('<filepath>')

from scipy.optimize import curve_fit

# %%
x = np.array([54, 51, 47, 42, 37, 32, 27, 22, 17, 12, 7])
y = np.array([610,610,590,530,446,345,235,142,70,24,13])
x = x[::-1]
y = y[::-1]

f  = lambda x, a, b, c, d: a*np.sin(b*x + c) + d
p0 = [300, 0.1, 3, 300]
popt, pcov = curve_fit(f, x, y, p0)#, [0.001]*len(x))
# popt = [300, 0.1, 3, 300]

xfit = np.linspace(0, 60, 100)
yfit = f(xfit, *popt)

br.figure()
plt.plot(x, y, marker='o')
plt.plot(xfit, yfit)


# %% corrected
x = np.array([54, 51, 47, 42, 37, 32, 27, 22, 17, 12, 7])
y2 = np.array([0.862973469, 0.862973469, 0.834679257, 0.749796621, 0.63096093, 0.488075159, 0.332456992, 0.200888906, 0.099029742, 0.033953055, 0.018391238])*100
x = x[::-1]
y2 = y2[::-1]

=H$2  *  SIN(H$3*J2 + H$4)   +    H$5
f  = lambda x, a, b, c, d: a*np.sin(b*x + c) + d
p0 = [300, 0.1, 3, 300]
popt, pcov = curve_fit(f, x, y2, p0)#, [0.001]*len(x))

xfit = np.linspace(0, 60, 100)
yfit = f(xfit, *popt)

br.figure()
plt.plot(x, y2, marker='o')
plt.plot(xfit, yfit)

# plt.close('all')
# %% ---------------------- Example 2: operations ------------------------- %% #

# %% two ways of modify a spectrum
s2.shift = 1
ss[2].factor = 0.5

# %% in fact, spectra in ss are the same in s1, s2, and s3
print(ss[1])
print(s1)

# %% plot
plt.figure()
_ = ss.plot()
plt.legend(['s0', 's1', 's2'])

# %% reorder
ss.reorder(0, 2)

# %% plot
plt.figure()
_ = ss.plot()
plt.legend(['s2', 's1', 's0'])

# %% flip order
ss.flip_order()

# %% plot
plt.figure()
_ = ss.plot()
plt.legend(['s0', 's1', 's2'])

# %% plot with vertical and horizontal increment
fig = br.backpack.figure()
_ = ss.plot(vi=100, hi=100)

# %% append and remove
ss.remove(2)
ss.append(s2)

# %% modifier is applied to all spectra
ss.offset = 5
plt.figure()
_ = ss.plot()
plt.legend(['s1', 's2', 's3'])

# %% another away
ss.set_offset(-5)
plt.figure()
_ = ss.plot()
plt.legend(['s1', 's2', 's3'])

# %% a list can be used
ss.offset = [1, 2, 3]
plt.figure()
_ = ss.plot()
plt.legend(['s1', 's2', 's3'])

# %% save and load
ss.save(folderpath='./', prefix='spectrum_', suffix='.dat')

# %% load from a folder filtering filenames
ss1 = br.Spectra()
ss1.load(folderpath='./', string='.dat')



# %% check methods and attributes
print(ss.step)
print(ss.length)
print(ss.x)
print(ss.monotonicity)

# %%
ss.check_length()
print(ss.length)

# %%
ss.check_monotonicity()
print(ss.monotonicity)

# %%
ss.check_same_x()
print(ss.length)
print(ss.x)

# %%
ss.check_step()
print(ss.length)
print(ss.step)
print(ss.x)

# %% a shift will reset ss.x, because the x is not the same for all spectra
ss1 = ss.copy()
ss1.shift = [10, 20, 30]
print(ss1[0].x)
print(ss1[1].x)
print(ss1[2].x)
print(ss1.step)
print(ss1.length)
print(ss1.x)
print(ss1.monotonicity)

# %% Area
print(ss.area)



# %% Interpolation
ss1 = ss.copy()
ss1.shift = [10, 20, 30]
print(ss1[0].x)
print(ss1[1].x)
print(ss1[2].x)
ss1.interp()
print(ss1[0].x)
print(ss1[1].x)
print(ss1[2].x)

ss2 = ss.copy()
ss2.shift = [10, 20, 30]
print(ss2[0].x)
print(ss2[1].x)
print(ss2[2].x)
ss2.interp(start=30, stop=5000, num=3000)
print(ss2[0].x)
print(ss2[1].x)
print(ss2[2].x)

# %% map =======================================================================
x0 = np.linspace(-5, 5, 100)
y0 = np.sin(x0)
s0 = br.Spectrum(x0, y0)

x1 = np.linspace(-5, 5, 100)
y1 = 0.6*np.sin(x1)
s1 = br.Spectrum(x1, y1)

x2 = np.linspace(-5, 5, 100)
y2 = 0.8*np.cos(x2)
s2 = br.Spectrum(x2, y2)

ss = br.Spectra(s0, s1, s2)

fig = br.figure()
ss.map.plot()





