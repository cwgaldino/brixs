# ==============================================================================
# %% EXAMPLE: Spectrum ====================================== 26/08/2022 =======
# ==============================================================================

# %% imports ===================================================================
import brixs as br
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import copy

%matplotlib qt5
%load_ext autoreload
%autoreload 2

# %% Peak object is just a modifed dictionary ==================================
p = br.Peak(amp=1, c=0, fwhm=0.2)
print(p)

# plot
fig = br.backpack.figure()
_ = p.spectrum.plot(marker='o', label='peak curve')
p.plot(label='peak position')
plt.legend()

# by default, peak asymmetry is set to False
print(p.asymmetry)

# setting fwhm1 or fwhm2 change asymmetry to True
p['fwhm1'] = 0.1
print(p.asymmetry)
print(p)  # fwhm2 is initialy set to fwhm-fwhm1

# fwhm is always fwhm1 + fwhm2
print(p['fwhm'])
p['fwhm2'] = 0.15
print(p['fwhm'])

# setting asymmetry back to false
p.asymmetry = False
print(p)


# save and load
p = br.Peak(amp=1, c=0, fwhm=0.2)
p.save()

p2 = br.Peak(amp=0, c=0, fwhm=1)
p2.load('Untitled.txt')
print(p2)

# %% Peaks object is just a modified list ======================================
p1 = br.Peak(amp=1, c=0, fwhm=0.2)
p2 = br.Peak(amp=2, c=1, fwhm=0.4)
ps = br.Peaks(p1, p2)

p1 = br.Peak(amp=1, c=0, fwhm=0.2)
p2 = br.Peak(amp=2, c=1, fwhm=0.4)
ps = br.Peaks([p1, p2])

p1 = br.Peak(amp=1, c=0, fwhm=0.2)
p2 = br.Peak(amp=2, c=1, fwhm=0.4)
ps = br.Peaks()
ps.append(p1)
ps.append(p2)

temp = [{'amp':1, 'c':0, 'fwhm':0.2}, {'amp':2, 'c':1, 'fwhm':0.4}]
ps = br.Peaks(temp)

# print
print(ps)

# enumeration
print(ps[0])
print(ps[-1])

# the peak with smaller center 'c' always comes first
ps = br.Peaks([p1, p2])
print(ps)
ps = br.Peaks([p2, p1])
print(ps)

# plot
fig = br.backpack.figure()
_ = ps.spectrum.plot(marker='o', label='peak curve')
ps.plot('peak position')
plt.legend()

# append and remove
ps.append(p1)
ps.remove(0)

# split and add_near
print(ps)
ps.split(1)
print(ps)

print(ps)
ps.add_near(0)
print(ps)


# save and load
temp = [{'amp':1, 'c':0, 'fwhm':0.2}, {'amp':2, 'c':1, 'fwhm':0.4}]
ps = br.Peaks(temp)
ps.save()

p2 = br.Peak(amp=0, c=0, fwhm=1)
p2.load('Untitled.txt')
print(p2)





# %% Advanced ==================================================================

# %% Peak methods for fitting purposes =========================================
p = br.Peak(amp=1, c=0, fwhm=0.2)

# bounds
print(p.bounds['c'])
p.set_bounds(c=0.1, type='additive')
print(p.bounds['c'])

print(p.bounds['amp'])
p.bounds['amp'] = (0.5, 2)
print(p.bounds['amp'])

print(p.bounds)

# build guess
guess, bounds_min, bounds_max, decode = p.build_guess()
print(guess)
print(bounds_min)
print(bounds_max)

# fix parameters
p.fixed  = ['m']
guess, bounds_min, bounds_max, decode = p.build_guess()
print(guess)

print(p)
print(p.fixed)
print(p.asymmetry)
f_str, args_str = p.build_model_str()
print(f_str)
print(args_str)

model = p.build_model()
model
# %% Peaks methods for fitting purposes ========================================
p = [{'amp':1, 'c':0, 'fwhm':0.2}, {'amp':2, 'c':1, 'fwhm':0.4}]
ps = br.Peaks(p)

p0, bounds_min, bounds_max, decode = ps.build_guess()
print(p0)

ps[0].fixed.append('amp')
f_str, args_str = ps.build_model_str()
print(f_str)
print(args_str)

model = ps.build_model()
model
