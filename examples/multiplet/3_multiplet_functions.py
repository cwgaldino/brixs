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
import brixs.addons.broaden
import brixs.multiplet as multiplet

# %% ============================== settings ============================= %% #
# multiplets
multiplet.settings.QUANTY_FILEPATH = r'C:\Users\galdin_c\github\quanty\quanty_win\QuantyWin64.exe'

# matplotlib (optional)
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()
# %%

# %  ===================================================================== %% #
# %  ================================ Basics ============================= %% #
# %% ===================================================================== %% #

# Initialization
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='XAS', edge='L2,3 (2p)')

# experimental geometry
q.plot_geometry()

# list of calculation attributes
q.get_attrs()

# list of methods available
q.get_methods() 

# get parameters in a dictionary
par = q.get_parameters()

# save parameters to a file
# this function is a little finicky because it used json package
# and this package can raise errors if parameters are not formatted right
# for example, it will raise an error if parameters is saved as a type np.int32
# however, it should be fine if one uses only "regular" float's and int's
q.save_parameters('test.par')

# parameters can be loaded by
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
s, out = q.run(update=True)
print(out)

# calculation parameters are saved as attrs
s.get_attrs()
s.initial

# plot calculated spectrum
plt.figure()
s.plot()

# Error messages are designed to inform available options
# for instance, if one does not know what edges are available for an experiment
# one can write nonsense and the error message should guide you
# see below the error message for this calculation
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='asdlkfjldakj;sdlkj')
