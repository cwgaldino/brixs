#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Model object for fitting.

################
# Introduction #
################

Lmfit is already rather straight forward, but two things are still a pain. On 
top of basic functionality, I needed the fitting module to be able to: (1) fit 
multiple spectra in a sequence passing the output of one fit as an input for the
 next fit (like when we are fitting spectra as a function of temperature), and 
(2) fit multiple spectra at once where fit parameters for specific features are 
the same for all spectra (for example, for when we want to fit the elastic line 
in all spectra but you want the constrain that the elastic line have the same 
width in all spectra).

The ideia is that both br.Spectrum() and br.Spectra() have an attribute called 
'model':

s = br.spectrum()
s.model

ss = br.Spectra()
ss.model

This attr 'model' is nothing more than a lmfit.Parameters object with a few 
extra functions attached to it. Note that, s.model is used for fitting one 
spectrum, while ss.model is used for multiple fittings.

One can define 'components' to a model (that's the best name I could think of). 
So one of the components that I defined is 'Peaks'. Therefore one can do:

s.model.peaks.add(amp=10, c=0, w=3)

And all parameters of a peak are added to the lmfit.Parameters object.


#########
# Model #
#########

`Model` is a lmfit.Parameters object (note the `s` in Parameters). 

The parent of `Model` is always a br.Spectrum or br.Spectra object.


brixs     Model   Component  Entries(i1)   Name    Tag   lmfit-parameter

                                        ┌── amp ── amp1 ── amp1_0_(i2)
                            ┌── entry 0 ├──  c  ──  c1  ──  c1_0_(i2)
                            |           └──  w  ──  w1  ──  w1_0_(i2)
                            |
                            |           ┌── amp ── amp1 ── amp1_0_(i2)
                 ┌──  peak  ├── entry 1 ├──  c  ──  c1  ──  c1_0_(i2)
                 |          |           └──  w  ──  w1  ──  w1_0_(i2)
                 |          └── ...  
                 |
Spectra ┐        |          ┌── entry 0
Spectum ├─ Model ├── arctan ├── entry 1
                 |          └── ...
                 |          ┌── entry 0
                 ├──  step  ├── entry 1
                 |          └── ...
                 └──  ...


#####################
# indexes i2 and i1 #
#####################

i2 is the spectrum number. When using br.Spectra object we have multiple spectra.
The index i2 indicates which entries are associated with which spectrum in br.Spectra.

i1 is the entry number. For example, one spectrum may have multiple peaks. Each 
peak will have a number i1.  


See more notes at: brixs/model/model.py
"""


# %% ========================== Standard imports ========================== %% #
import matplotlib.pyplot as plt
import numpy as np
plt.ion()
# %%

# %% ============================ brixs imports =========================== %% #
import brixs as br
import brixs.model
# %%

# %% ============================= example 1 ============================== %% #
# fit a spectrum with two peaks

#############################################
# create (simulate) a spectrum to be fitted #
#############################################
s = br.Spectrum()

# add two entries from one component
s.model.peaks.add(amp=10, c=-4, w=4)
s.model.peaks.add(amp=5, c=6, w=5)

# simulated two peaks
x = s.model.peaks.get_suitable_x()
f = s.model.peaks.get_model()
y = f(x)

# manually assign x and y
s._x = x[::-1]
s._y = y[::-1] + np.random.normal(loc=0, scale=0.5, size=len(x))

# clear peaks (remove model)
s.model.peaks.clear()

# plot
br.figure()
s.plot()
plt.title('spectrum s has two peaks')

################
# fit spectrum #
################

# model information should be empty -> print model information
s.model.pretty_print()

# add `guess` peaks
s.model.peaks.add(amp=8, c=-3, w=5)
s.model.peaks.add(amp=4, c=5, w=4)

# plot before fitting
br.figure()
s.plot(marker='o', color='black', lw=0, label='data')
s.model.plot(color='green', label='guess')
br.leg()
plt.title('guess (initial) fit values')

# fit
report, out = s.model.fit()
print(report)

# plot after
br.figure()
s.plot(marker='o', color='black', lw=0, label='data')
s.model.plot(color='red', label='fit')
ss = s.model.peaks.calculate_spectrum_for_each_i1()
ss[0].label = 'peak 0'
ss[1].label ='peak 1'
ss.plot(color='black', lw=0.5)
br.leg()
plt.title('fitted curve')

# print parameters
print('='*20)
s.model.pretty_print()
print('='*20)
s.model.peaks.pretty_print()
# %%

# %% ============================= example 2 ============================== %% #
# fit one spectrum two peaks and a linear background

#############################################
# create (simulate) a spectrum to be fitted #
#############################################
s = br.Spectrum()

# add two entries from one component
s.model.peaks.add(amp=10, c=-4, w=4)
s.model.peaks.add(amp=5, c=6, w=5)
s.model.linear.add(linear=0.1, const=2)

# simulated two peaks
x = s.model.get_suitable_x()
f = s.model.get_model()
y = f(x)

# manually assign x and y
s._x = x[::-1]
s._y = y[::-1] + np.random.normal(loc=0, scale=0.5, size=len(x))

# clear peaks (remove model)
s.model.clear()

# plot
br.figure()
s.plot()
plt.title('spectrum s has two peaks and a linear bkg')

################
# fit spectrum #
################
# model information should be empty -> print model information
s.model.pretty_print()

# add `guess` peaks
s.model.peaks.add(amp=8, c=-3, w=5)
s.model.peaks.add(amp=4, c=5, w=4)
s.model.linear.add(linear=0.3, const=1)

# plot before fitting
br.figure()
s.plot(marker='o', color='black', lw=0, label='data')
s.model.plot(color='green', label='guess')
br.leg()
plt.title('guess (initial) fit values')

# fit
report, out = s.model.fit()
print(report)

# plot after
br.figure()
s.plot(marker='o', color='black', lw=0, label='data')
s.model.plot(color='red', label='fit')
ss = s.model.peaks.calculate_spectrum_for_each_i1()
ss[0].label = 'peak 0'
ss[1].label ='peak 1'
ss.plot(color='black', lw=0.5)
br.leg()
plt.title('fitted curve')

# print parameters
print('='*20)
s.model.pretty_print()
print('='*20)
s.model.peaks.pretty_print()
print('='*20)
s.model.linear.pretty_print()
# %%

# %% ============================= example 3 ============================== %% #
# test if base-modifiers work

################################
# create (simulate) a spectrum #
################################
s = br.Spectrum()

# add two entries from one component
s.model.peaks.add(amp=10, c=-4, w=4)
s.model.peaks.add(amp=5, c=6, w=5)
s.model.linear.add(linear=0.1, const=2)

# simulated two peaks
x = s.model.get_suitable_x()
f = s.model.get_model()
y = f(x)

# manually assign x and y
s._x = x[::-1]
s._y = y[::-1] + np.random.normal(loc=0, scale=0.5, size=len(x))

br.figure()
s.plot()
s.model.plot()
plt.title('spectrum s with model')

#############
# set shift #
#############
s1 = s.set_shift(10)

print(s.model.peaks['c_0_0'].value)
print(s1.model.peaks['c_0_0'].value)

br.figure()
s1.plot()
s1.model.plot()
plt.title('shifted spectrum (model also changes)')

##############
# set offset #
##############
s1 = s.set_offset(10)

print(s.model.linear['const_0_0'].value)
print(s1.model.linear['const_0_0'].value)

br.figure()
s1.plot()
s1.model.plot()
plt.title('spectrum with an offset (model also changes)')

##############
# set factor #
##############
s1 = s.set_factor(10)

print(s.model.peaks['amp_0_0'].value)
print(s1.model.peaks['amp_0_0'].value)

br.figure()
s1.plot()
s1.model.plot()
plt.title('spectrum with a multiplicative factor (model also changes)')

#############
# set calib #
#############
s1 = s.set_calib(10)

print(s.model.peaks['w_0_0'].value)
print(s1.model.peaks['w_0_0'].value)

br.figure()
s1.plot()
s1.model.plot()
plt.title('spectrum with a different calibration (model also changes)')
# %%

# %% ============================= example 4 ============================== %% #
# scipy find peaks

#############################################
# create (simulate) a spectrum to be fitted #
#############################################
s = br.Spectrum()

# add two entries from one component
s.model.peaks.add(amp=10, c=-4, w=4)
s.model.peaks.add(amp=5, c=6, w=5)
s.model.linear.add(linear=0.1, const=2)

# simulated two peaks
x = s.model.get_suitable_x()
f = s.model.get_model()
y = f(x)

# manually assign x and y
s._x = x[::-1]
s._y = y[::-1] + np.random.normal(loc=0, scale=0.5, size=len(x))

# clear peaks (remove model)
s.model.clear()

# plot
br.figure()
s.plot()
plt.title('spectrum s has two peaks and a linear bkg')

#################################################################
# estimate number of peaks, peak position, amplitude, and width #
#################################################################
print(s.model.peaks)  # -> no peaks
s.model.peaks.find(width=32)
print(s.model.peaks)  # -> peaks

br.figure()
s.plot()
s.model.peaks.plot()
plt.title('peaks found via find_peaks()')
# %%

# %% ============================= example 5 ============================== %% #
# fit three spectra with two peaks and a linear background each

#############################################
# create (simulate) a spectrum to be fitted #
#############################################
def data(amp1, c1, w1, amp2, c2, w2, linear, const):
    s = br.Spectrum()
    s.model.peaks.add(amp=amp1, c=c1, w=w1)
    s.model.peaks.add(amp=amp2, c=c2, w=w2)
    s.model.linear.add(linear=linear, const=const)

    x = s.model.get_suitable_x()
    f = s.model.get_model()
    y = f(x)

    # manually assign x and y
    s._x = x[::-1]
    s._y = y[::-1] + np.random.normal(loc=0, scale=0.2, size=len(x))

    # clear peaks (remove model)
    s.model.clear()

    return s

s1 = data(amp1=10, c1=1, w1=1, amp2=6, c2=4, w2=2, linear=0.1, const=0.5)
s2 = data(amp1=10, c1=1, w1=1, amp2=8, c2=4, w2=2, linear=0.1, const=0.5)
s3 = data(amp1=10, c1=1, w1=1, amp2=10, c2=4, w2=2, linear=0.1, const=0.5)

ss = br.Spectra([s1, s2, s3])

# plot
br.figure()
ss.plot()
plt.title('Spectra')

###############
# build model #
###############
ss.model.peaks.add(amp=8, c=2, w=2, i2='all')
ss.model.peaks.add(amp=8, c=4, w=2, i2='all')
ss.model.linear.add(linear=0.3, const=0, i2='all')
ss.model.pretty_print()

# plot before fitting
fig, axes = br.subplots(1, 3, sharex=True, sharey=True, wspace=0)
ss[0].plot(ax=axes[0], marker='o', color='black', lw=0, label='data 0')
ss[1].plot(ax=axes[1], marker='o', color='black', lw=0, label='data 1')
ss[2].plot(ax=axes[2], marker='o', color='black', lw=0, label='data 2')

ss.model[0].plot(ax=axes[0], color='green', label='guess')
ss.model[1].plot(ax=axes[1], color='green', label='guess')
ss.model[2].plot(ax=axes[2], color='green', label='guess')

br.leg(ax=axes[0])
br.leg(ax=axes[1])
br.leg(ax=axes[2])

################
# fit spectrum #
################
# fit
report, out = ss.model.fit()
print(report)

# plot before fitting
fig, axes = br.subplots(1, 3, sharex=True, sharey=True, wspace=0)
ss[0].plot(ax=axes[0], marker='o', color='black', lw=0, label='data 0')
ss[1].plot(ax=axes[1], marker='o', color='black', lw=0, label='data 1')
ss[2].plot(ax=axes[2], marker='o', color='black', lw=0, label='data 2')

ss.model[0].plot(ax=axes[0], color='red', label='fit')
ss.model[1].plot(ax=axes[1], color='red', label='fit')
ss.model[2].plot(ax=axes[2], color='red', label='fit')

br.leg(ax=axes[0])
br.leg(ax=axes[1])
br.leg(ax=axes[2])

# print parameters
print('='*20)
ss.model.pretty_print()
print('='*20)
ss.model.peaks.pretty_print()
print('='*20)
ss.model[0].peaks[1]['amp']
ss.model[1].peaks[1]['amp']
ss.model[2].peaks[1]['amp']

# TODO
# a better syntax would be:
# ss.model[0].peaks[1]['amp'] # done -> returns `amp` for i2=0, i1=1
# ss.model[0].peaks['amp']    # -> returns all `amp` parameters for spectrum i2=0
# ss.model.peaks['amp']       # -> returns all `amp` parameters for all spectra
# ss.model.peaks[1]['amp']    # -> returns the `amp` parameter of peak i1=1 for all spectra
# %%

# %% ============================= example 6 ============================== %% #
# fit three spectra with two peaks and a linear background each
# with tied parameters

#############################################
# create (simulate) a spectrum to be fitted #
#############################################
def data(amp1, c1, w1, amp2, c2, w2, linear, const):
    s = br.Spectrum()
    s.model.peaks.add(amp=amp1, c=c1, w=w1)
    s.model.peaks.add(amp=amp2, c=c2, w=w2)
    s.model.linear.add(linear=linear, const=const)

    x = s.model.get_suitable_x()
    f = s.model.get_model()
    y = f(x)

    # manually assign x and y
    s._x = x[::-1]
    s._y = y[::-1] + np.random.normal(loc=0, scale=0.2, size=len(x))

    # clear peaks (remove model)
    s.model.clear()

    return s

s1 = data(amp1=10, c1=1, w1=1, amp2=6, c2=4, w2=2, linear=0.1, const=0.5)
s2 = data(amp1=10, c1=1, w1=1, amp2=8, c2=4, w2=2, linear=0.1, const=0.5)
s3 = data(amp1=10, c1=1, w1=1, amp2=10, c2=4, w2=2, linear=0.1, const=0.5)

ss = br.Spectra([s1, s2, s3])

# plot
br.figure()
ss.plot()
plt.title('Spectra')

###############
# build model #
###############
ss.model.peaks.add(amp=8, c=2, w=2, i2='all')
ss.model.peaks.add(amp=8, c=4, w=2, i2='all')
ss.model.linear.add(linear=0.3, const=0, i2='all')
ss.model.pretty_print()

# tie second peak amplitude
ss.model['amp1_1_1'].expr = 'amp1_1_0'
ss.model['amp1_1_2'].expr = 'amp1_1_0'

# TODO -> raise an error for this syntax
# ss.model.peaks['amp1_1_1'].expr = 'amp1_1_0'

# TODO -> after setting expr, I cannot 'break' ss.model anymore
# i.e., if I do this, ss.model[0], then ss.model stop working
# it is a lmfit thing. To fix it, I will have to set a expr2
# then I copy expr2 into expr right before the fitting and I 
# delete expr right after. For now, we gotta remember to remove
# expr by hand after the fitting. 

################
# fit spectrum #
################
# fit
report, out = ss.model.fit()
print(report)

# remove expr by hand
ss.model['amp1_1_1'].expr = ''
ss.model['amp1_1_2'].expr = ''
ss.model.pretty_print()

# plot before fitting
fig, axes = br.subplots(1, 3, sharex=True, sharey=True, wspace=0)
ss[0].plot(ax=axes[0], marker='o', color='black', lw=0, label='data 0')
ss[1].plot(ax=axes[1], marker='o', color='black', lw=0, label='data 1')
ss[2].plot(ax=axes[2], marker='o', color='black', lw=0, label='data 2')

ss.model[0].plot(ax=axes[0], color='red', label='fit')
ss.model[1].plot(ax=axes[1], color='red', label='fit')
ss.model[2].plot(ax=axes[2], color='red', label='fit')

br.leg(ax=axes[0])
br.leg(ax=axes[1])
br.leg(ax=axes[2])

# print parameters
print('='*20)
ss.model.pretty_print()
print('='*20)
ss.model.peaks.pretty_print()
print('='*20)
print(ss.model[0].peaks['amp_1_0'])
print(ss.model[1].peaks['amp_1_1'])
print(ss.model[2].peaks['amp_1_2'])
# %%
