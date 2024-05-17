#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""A Spectra object `ss` is a list that stores Spectrum objects. 

>>> ss = br.Spectra(s0, s1, s2)

`ss` works as a list and indexing available

>>> s = ss[1]

slicing is also possible

>>> ss2 = ss[0:1]

New Spectrum `s` can be added

>>> ss.append(s)

Methods operate on the Spectrum list

>>> print(s[0].y[1])
>>> ss.set_offset(2)
>>> print(s[0].y[1])

help for each method is available through python's built-in help function

>>> help(ss.set_offset)

Attrs can be added to `ss` on the go

>>> ss.new_attr = 10

New methods can also be added on the go (see advanced examples)

>>> br.Spectra.new_method = new_method

See below for a more complete description.
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

# iteration over each spectra
for s in ss:
    print(s)

# attributes can be associated with Spectra object
ss.temperature = 10
ss.composition = [9, 3, 5, 6]
ss.string      = 'test'
ss.dictionary  = {'a': 10, 'b':5}
print(ss.temperature)
print(ss.composition)
print(ss.string)
print(ss.dictionary)

# attributes can be associated with individual spectrum
ss[0].name = 'spectrum0'
ss[1].name = 'spectrum1'
ss[2].name = 'spectrum2'

# get a list of all defined attrs and remove all attrs
print(ss.get_attrs())
# ss.remove_attrs()

# multiple spectra can be saved in one file
ss.interp()  # interpolate data so every spectrum has the same x-axis
ss.save_all_single_file(filepath='examples/spectra.dat')

# multiple spectra can be saved in multiple files
ss.save('examples/')

# Load spectra from single file
ss = br.Spectrum('examples/spectra.dat')

# Load spectra from multiple files
# ss = br.Spectrum('examples/')  # ERROR
# using this syntax, Spectra can only be loaded if they are the only files in the folder

# spectra can also be loaded via load function, which allow for more options
# loading from multiple files:
ss1 = br.Spectra()
ss1.load(folderpath='examples/', string='spectrum', comments='#')
# loading from single file
ss2 = br.Spectra()
ss2.load_from_single_file(filepath='examples/spectra.dat', comments='#')
# for instance, `string` selects only files where the filename contains a certain string
# Note that attrs are saved to the files
# saving all spectra in one file saves the attrs of the Spectra object
# saving each individual spectra in multiple files saves attrs for each spectrum
# if you want to save both the Spectra and Spectrum attrs you have to save files 
# via both methods (save() and save_all_single_file()), then you can use the
# use the function copy_attrs_from() to copy attrs
ss1.copy_attrs_from(ss2)

# %% ---------------------- Example 2: reordering ------------------------- %% #
# creating spectra object
x0 = np.linspace(-5, 5, 100)
y0 = np.sin(x0)
s0 = br.Spectrum(x0, y0)
s0.T = 10

x1 = np.linspace(-5, 5, 100)
y1 = 0.6*np.sin(x1)
s1 = br.Spectrum(x1, y1)
s1.T = 20

x2 = np.linspace(-5, 5, 100)
y2 = 0.8*np.cos(x2)
s2 = br.Spectrum(x2, y2)
s2.T = 30

ss = br.Spectra(s0, s1, s2)
print(len(ss))

# items can be appended
ss.append(s1)
print(len(ss))

ss.append([s1, s2])
print(len(ss))

# item can be removed
ss.remove(-1)
del ss[-1]

# items can be assigned
ss[-1] = s2

# items can be reordered
ss.reorder(0, 1)

# order of the items can be flipped
ss.flip_order()

# one can reorder spectra based on a attr
ss.create_attr_from_spectra(attr='T')
print(ss.T)
ss.reorder_by_attr(attr='T')
print(ss.T)

# get spectrum based on a attr
s = ss.get_by_attr(attr='T', value=10)

# %% ---------------------- Example 3: x axis verification ---------------- %% #
# creating spectra object
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

# check if all spectra have the same x axis
print(ss.x)
ss.check_same_x()
print(ss.x)

# check if all spectra have the same length
for s in ss:
    print(len(s))
print(ss.length)  # None
ss.check_length()
print(ss.length)  # 100

# check if all spectra has the same separation value between two items of the x-axis
print(ss.step)          # None
s.check_step()         # calculate step
print(s.step)          # 0.10101

# Irregular steps will raise an error
# s2.x[-1] = 5.1   # changing one value just to raise an error
# s2.check_step()  # ERROR
# s2.x[-1] = 5

# s.monotonicity returns the x-axis monotonicity (`increasing`, `decreasing`)
print(ss.monotonicity)   # None
ss.check_monotonicity()  # calculate monotonicity
print(ss.monotonicity)   # increasing

# irregular monotonicity will raise and error
# s2.x[-1] = 0                # changing one value just to raise an error
# s.check_monotonicity()  # ERROR
# s.fix_monotonicity()    # fix monotonicity
# s.check_monotonicity()  # calculate monotonicity after fixing
# print(ss.monotonicity)


# %% ----------------------- Example 4: base modifiers -------------------- %% #
# creating spectra object
ss = br.Spectra()

# spectra has 3 base modifiers acting on the x-axis
ss.calib
ss.shift
ss.roll

# and 2 base modifiers acting on the y-axis
ss.offset
ss.factor

# modifiers are applied to ALL spectra inside Spectra object
# modifiers can be set via the attribute or via function
ss.offset = 30
ss.set_offset(30)
# ss.set_offset([10, 12, 22, ...])  # different value for each spectra

# calib will multiply the x-axis by a number
# shift will add a value to the x-axis
# roll will roll the x-axis by an integer value (see np.roll())
# offset will add a value to the y-axis
# factor will multiply the y-axis by a number

# %% ----------------------- Example 5: modifiers ------------------------- %% #
# ss = br.Spectra()

# these are the defined modifiers
# ss.floor()     # uses s.offset to make the average value of the y-axis zero
# ss.flip()      # uses s.calib to flip the x-axis (multiply by -1) 
# ss.normalize() # uses s.offset to make the average value of the y-axis goes to value
# ss.interp()    # interpolate data
# ss.crop()      # remove the start and end of the dataset
# ss.switch()    # switch x- and y-axis
# ss.align()     # uses cross-correlation and ss.shift to align data (make y-axis match)

# %% ----------------------- Example 6: copying --------------------------- %% #
ss1 = br.Spectra()

# coping syntax 1
ss2 = br.Spectra()
ss2.copy(ss1) 

# coping syntax 2
ss2 = ss1.copy() 

# %% ---------------------- Example 7: calculate -------------------------- %% #
s = br.Spectra()

# returns a spectra type (loops through each spectrum inside Spectra object)
sder   = ss.derivative()
smooth = ss.smooth()

# returns a list (loops through each spectrum inside Spectra object)
area = ss.calculate_area()
xsum = ss.calculate_x_sum()
ysum = ss.calculate_y_sum()

# returns a Spectrum type
s = ss.concatenate() 
s = ss.calculate_sum()
s = ss.calculate_average()

# returns a Image type (2D representation of Spectra), see Image example file
im = ss.calculate_map()

# simple polynomial fitting (loops through each spectrum inside Spectra object)
popt, fit, R2 = ss.polyfit(deg=2)

# %% ------------ Example 8: base modifier calculation values ------------- %% #
# Calculate shift values (see modifiers example) to align datasets
ss.calculate_shift()  
print(ss.calculated_shift)  # list

# Calculate roll values (see modifiers example) to align datasets
ss.calculate_roll()
print(ss.calculated_roll)  # list

# Calculate factor values (see modifiers example) to normalize datasets
ss.calculate_factor()
print(ss.calculated_factor)  # list

# Calculate offset values (see modifiers example) to align datasets vertically
ss.calculate_offset()
print(ss.calculated_offset)  # list

# Calculate a calibration factor for the x-axis
# Assumes each spectrum inside Spectra relates to a `value`
# one can find a calibration factor so x-axis is given in terms of `value`
ss.calculate_calib(values=[...])
print(ss.calculated_shift)  # spectrum type

# Note
# after using any of the base-modifier-calculate functions, one can use the base
# modifier function to directly apply the modification, for example, this
# >>> ss.calculate_shift()  
# >>> ss.set_shift(ss.calculated_shit)
# is the same as this
# >>> ss.calculate_shift()  
# >>> ss.set_shift()  
# the set_shift function without any argument will try to find calculated values 

# %% ----------------------- Example 9: plotting -------------------------- %% #
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

# basic plotting
br.figure()  # same as plt.figure(), but with additional features
lines = ss.plot()    # same as doing plt.plot(s.x, s.y) for each s inside ss

# quick adjustments can be made
br.figure() 
lines = ss.plot(shift=1, offset=1)     

# plot with vertical and horizontal increments (percentage wise)
br.figure()
lines = ss.plot(vertical_increment=50, horizontal_increment=50)

# you can get the modifiers values from lines
print([line.shift for line in lines])
print([line.offset for line in lines])
print([line.calib for line in lines])
print([line.factor for line in lines])

# note, opening a new figure with br.figure() allows for sending the x and y 
# values where the mouse cursor is to the the clipboard by left and right 
# mouse button click. Pressing the wheel copies the image. Double left click
# prints the figure coordinates (from 0 to 1).

# sequential_plot [EXPERIMENTAL]
ss.sequential_plot()  # opens a new figure where one can use the keyboard arrows to flip through spectra

# shift_plot and roll_plot [EXPERIMENTAL]
ss.roll_plot()  
# ss.shift_plot()  
# opens a new figure
# flip through spectra using up and down arrows in the keyboard
# shift or roll spectra using left and right arrows in the keyboard
# shift or roll results are saved in ss.calculated_shift or ss.calculated_roll
print(ss.calculated_roll)
# print(ss.calculated_shift)

