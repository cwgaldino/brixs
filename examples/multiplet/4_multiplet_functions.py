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
# del par['initial']
# del par['geometry']
# del par['experiment']
# par.keys()
# import json
# filepath = Path('test.par')
# pretty_print = True
# with open(str(filepath), 'w') as file:
#     if pretty_print:
#         file.write(json.dumps(par, indent=4, sort_keys=False))
#     else:
#         file.write(json.dumps(par))


# [EXPERIMENTAL] save/load parameters to a file
# thise functions are a little finicky because it used json package
# and this package can raise errors if parameters are not formatted right
# for example, it will raise an error if parameters is saved as a type np.int32
# however, it should be fine if one uses only "regular" float's and int's
# q.save_parameters('test.par')
# q2 = multiplet.load_calculation('test.par')
# for now, I will leave them commented out until I ran more tests. However, they
# should work fine and can be used.

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
