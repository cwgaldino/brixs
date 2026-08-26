#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Basic functions and methods

Author: Carlos Galdino
Last updated 08/09/2025
"""

# %% ========================== Standard imports ========================== %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# %% ============================ brixs imports =========================== %% #
import brixs as br
import brixs.multiplet as multiplet

# %% ============================== settings ============================= %% #
# multiplet
# multiplet.settings.QUANTY_FILEPATH = r'C:\Users\galdin_c\github\quanty\quanty_win\QuantyWin64.exe'
multiplet.settings.QUANTY_FILEPATH = r'/Users/oax12540/github/quanty/2024Spring/QuantyMac'

# matplotlib (optional)
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()
# %%

# %  ===================================================================== %% #
# %  ================================ Basics ============================= %% #
# %% ===================================================================== %% #
# initialization
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')

# experimental geometry
q.plot_geometry()

# list of calculation attributes
q.get_attrs()

# list of methods available
q.get_methods() 

# get parameters in a dictionary
par = q.get_parameters()

# [EXPERIMENTAL] save/load parameters to a file
# thise functions are a little finicky because it used json package
# and this package can raise errors if parameters are not formatted right
# for example, it will raise an error if parameters is saved as a type np.int32
# however, it should be fine
q.save_parameters('test.par')
q2 = multiplet.load_calculation('test.par')

# lua template
print(q.template)
q.save_template('template.lua')  # template can be saved to a file

# lua script
# this method uses the "template" to create a lua script
q.update_lua_script()
print(q.lua_script)

# the lua script can be saved
q.save_lua_script('test.lua')

# run lua scripts independently from the calculation object
out = multiplet.quanty('test.lua')

# run calculation
s, out = q.run()
print(out)

# calculation parameters are stored inside the Spectrum object as attrs
s.get_attrs()
s.initial

# plot calculated spectrum
plt.figure()
s.plot()

# if more than one spectrum is created during the calculation, q.run() returns
# a dictionary
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')
q.polarization = 'linear'
data, out = q.run()
print(data.keys())

br.figure()
data['vv'].plot()
data['hh'].plot()
print(data['vv'].get_attrs())