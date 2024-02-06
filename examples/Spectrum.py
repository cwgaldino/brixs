#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Spectrum object is an object that stores a x and y array 

It has pre-defined methods for operating in this x and y arrays

New methods and attrs can be defined on the go
"""

# %% ---------------------------- imports --------------------------------- %% #
import brixs as br
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# %% ---------------------- quick Example --------------------------------- %% #
x = np.linspace(-5, 5, 100)
y = np.sin(x)

s = br.Spectrum(x, y)

plt.figure()
s.plot()

# %% ---------------------- Example 1: Initialization --------------------- %% #
x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y = [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]

# Spectrum object can be initialized in many ways
s1 = br.Spectrum(x=x, y=y)
s1 = br.Spectrum(x, y)
s1 = br.Spectrum(y)
s1 = br.Spectrum(y=y)

# x and y data
print(s1.x)
print(s1.y)
print(s1.data)

# attributes can be associated with spectrum object
s1.temperature = 10
s1.composition = [9, 3, 5, 6]
s1.string      = 'test'
s1.dictionary  = {'a': 10, 'b':5}
print(s1.temperature)
print(s1.composition)
print(s1.string)
print(s1.dictionary)

# get a list of all defined attrs and remove all attrs
print(s1.get_attrs())
# s1.remove_attrs()

# item selection will return the x, y pair
print(s[0])

# items can be deleted
print(s.x)
del s[0]
print(s.x)

# as for right now, items cannot be changed (this might change in the future)
# s[2] = (9, 10)  # error

# one can save and load an spectrum
s1.save('examples/spectrum_0.dat')
s2 = br.Spectrum('examples/spectrum_0.dat')

# attrs are saved and loaded with spectrum
# Note that, dictionaries and more complex attrs are not saved
print(s2.temperature)
# print(s2.dictionary)  # ERROR

# attrs can be copied between spectrum objects
s2.copy_attrs_from(s1)
print(s2.dictionary)


# %% ---------------------- Example 2: operations ------------------------- %% #
s1 = br.Spectrum(y=[1, 2, 3, 4, 5])
s1.temperature = 10
s2 = br.Spectrum(y=[5, 4, 3, 2, 1])
s2.temperature = 20

# Spectrum can be operated on (y axis)
s3 = s1 + s2
s3 = s1 - s2
s3 = s1 * s2
s3 = s1 / s2  # in this case, s_C cannot contain zeros on the y axis
s3 = s1 + 2.2
s3 = s1 - 2
s3 = s1 * 2
s3 = s1 / 2

# Note that this does not work
# s3 = 2 * s1  # ERROR
# but this works
s3 = s1 * 2

# attrs are not transferred
# print(s3.temperature)  # ERROR

# %% ---------------------- Example 3: x axis verification ---------------- %% #
x = [0, 1, 2, 3, 4, 5, 6]
y = [1, 1, 1, 1, 1, 1, 1]

# s.step returns the separation value between two items of the x-axis
s = br.Spectrum(x, y)  # initialization
print(s.step)          # None
s.check_step()         # calculate step
print(s.step)          # 1

# Irregular steps will raise an error
x = [0, 1.1, 2, 2.3, 7, 10, 12]
s = br.Spectrum(x, y)
# s.check_step()  # ERROR

# s.monotonicity returns the x-axis monotonicity (`increasing`, `decreasing`)
s = br.Spectrum(x, y)   # initialization
print(s.monotonicity)   # None
s.check_monotonicity()  # calculate monotonicity
print(s.monotonicity)   # increasing

# irregular monotonicity will raise and error
x=[0, 12, 5, 1, 8, 6, 2]
s = br.Spectrum(x, y)
# s.check_monotonicity()  # ERROR
s.fix_monotonicity()    # fix monotonicity
s.check_monotonicity()  # calculate monotonicity after fixing
print(s.monotonicity)
print(s.x)

# %% ----------------------- Example 4: bese modifiers -------------------- %% #
x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y = [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]
s = br.Spectrum(x, y)

# spectrum has 3 base modifiers acting on the x-axis
s.calib
s.shift
s.roll

# and 2 base modifiers acting on the y-axis
s.offset
s.factor

# modifiers can be set via the attribute or via function
s.offset = 30
s.set_offset(30)

# calib will multiply the x-axis by a number
# shift will add a value to the x-axis
# roll will roll the x-axis by an integer value (see np.roll())
# offset will add a value to the y-axis
# factor will multiply the y-axis by a number

# %% ----------------------- Example 5: modifiers ------------------------- %% #
x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y = [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]
s = br.Spectrum(x, y)

# these are the defined modifiers
s.floor()      # uses s.offset to make the average value of the y-axis zero
s.flip()       # uses s.calib to flip the x-axis (multiply by -1) 
s.normalize()  # uses s.offset to make the average value of the y-axis goes to value
# s.interp()   # interpolate data
s.crop()       # remove the start and end of the dataset
s.switch()     # switch x- and y-axis

# %% ----------------------- Example 6: copying --------------------------- %% #
x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y = [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]
s1 = br.Spectrum(x, y)

s2 = br.Spectrum()
s2.copy(s1)

s2 = s1.copy()

# %% ---------------------- Example 7: calculate -------------------------- %% #
x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y = [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]
s = br.Spectrum(x, y)

# returns a spectrum type
sder   = s.derivative()
smooth = s.smooth()

# returns a number
area = s.calculate_area()
xsum = s.calculate_x_sum()
ysum = s.calculate_y_sum()
area = s.calculate_area()

# index finding
i  = s.index(x=4)
yv = s.x2y(x=4)

# simple polynomial fitting
popt, fit, R2 = s.polyfit(deg=2)

# %% ----------------------- Example 8: plotting -------------------------- %% #
x = np.linspace(0, 10, 100)
y = br.gaussian(x=x, amp=4, c=4, sigma=0.4) + br.gaussian(x=x, amp=2, c=8, sigma=0.2)
s = br.Spectrum(x, y)

br.figure()  # same as plt.figure(), but with additional features
s.plot()     # same as plt.plot(s.x, s.y)
s.plot(shift=1, offset=1)     # quickly adjustments can be made

# note, opening a new figure with br.figure() allows for sending the x and y 
# values where the mouse cursor is to the the clipboard by left and right 
# mouse button click. Pressing the wheel copies the image. Double left click
# prints the figure coordinates (from 0 to 1).

# %% --------------------- Example 9: peak finding ------------------------ %% #
x = np.linspace(0, 10, 100)
y = br.gaussian_fwhm(x=x, amp=4, c=4, w=1) + br.gaussian_fwhm(x=x, amp=2, c=8, w=1)
s = br.Spectrum(x, y)

# uses scipy.signal.find_peaks()
s.find_peaks()
print(s.peaks)  # print list of peaks found

print(s.peaks[0])
print(s.peaks[1])

# peaks can be removed
s.peaks.remove(0)

# peaks can be added
s.peaks.append(amp=4, c=4, w=1)

# peaks can be modifies
s.peaks[1]['amp'].value = 3

# one can calculate a spectrum from peaks
s2 = s.peaks.calculate_spectrum()

# quick plot
br.figure()  
s.plot(marker='o') 
s.peaks.plot()     # place a marker at every peak found  
s2.plot()

# %% --------------------- Example 10: peak fitting ----------------------- %% #
x = np.linspace(0, 10, 100)
y = br.gaussian_fwhm(x=x, amp=4, c=4, w=1) + br.gaussian_fwhm(x=x, amp=2, c=8, w=1)
s = br.Spectrum(x, y)

# fit one peak
s.fit_peak(ranges=(0, 5))
print(s.peaks)
one_peak = s.peaks.calculate_spectrum()

# find all peaks and fit
s.find_peaks()
s.fit_peaks() 
fit = s.peaks.calculate_spectrum()

# quick plot
br.figure()  
s.plot(marker='o') 
fit.plot()
one_peak.plot()   

# %% ----------------- Example 11: defining new methods ------------------- %% #
x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y = [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]
s1 = br.Spectrum(x, y)

# define new function that acts on x and y
def new_function(self, value):
    """Adds a value to the y-axis"""
    self.y = self.y + value

# attaches function to type Spectrum
br.Spectrum.new_function = new_function

# from now one, every spectrum type will have new_function defined.
# even spectrum types created before defying new_function
s2 = br.Spectrum(x, y)

s1.new_function(5)
print(s1.y)

s2.new_function(2)
print(s2.y)