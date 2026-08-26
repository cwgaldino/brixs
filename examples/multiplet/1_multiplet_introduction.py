#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""The brixs.multiplet module provides tools for multiplet calculations of XPS,
XES XAS, and RIXS spectra

This module is based on Crispy 0.7.3 written by Marius Retegan. Therefore, it 
is instructive to have some familiarity Crispy. The latest version can be 
downloaded at,

https://www.esrf.fr/computing/scientific/crispy/

This implementation, as well as Crispy, requires a Quanty executable 
from Maurits W. Haverkort. It that can be downloaded at,

https://www.quanty.org/index.html

Therefore, if you use this module on publications (brixs.multiplet). Please, 
cite Crispy 0.7.3 and Quanty. You can refer to their respective webpages to 
learn how to do so.

Note that this module was only tested on windows and for a limited number of 
calculation parameters. Please, report bugs!    

Author: Carlos Galdino
Last updated 08/09/2025
"""

# %% ===================================================================== %% #
# %% ================== 1) import brixs related modules ================== %% #
# %% ===================================================================== %% #
import brixs as br                   # core brixs package (optional)
import brixs.addons.broaden          # used for broadening the calculated spectra
import brixs.multiplet as multiplet  # multiplet calculation module
# %%

# %% ===================================================================== %% #
# %%========= 2) tell multiplet where the Quanty executable is %%========= %% #
# %% ===================================================================== %% #
# multiplet.settings.QUANTY_FILEPATH = r'C:\Users\galdin_c\github\quanty\quanty_win\QuantyWin64.exe'
multiplet.settings.QUANTY_FILEPATH = r'/Users/oax12540/github/quanty/2024Spring/QuantyMac'
# %%

# %% ===================================================================== %% #
# %% ============= 3) initialize a calculation object (`q`) ============== %% #
# %% ===================================================================== %% #
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')

# these are the minimum parameters necessary to initialize a calculation object:
# element, charge, symmetry, experiment, edge
# primary parameters cannot be changed once the `q` is created
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')

# you can also initialize it with the default values
q = multiplet.Calculation()
print(q.element, q.charge, q.symmetry, q.experiment, q.edge)

# do not be afraid of making mistakes when initializing a calculation object
# Error messages are designed to inform available options
# for instance, if one does not know what edges are available for an experiment
# one can write nonsense and the error message should guide you
# see below the error message for this calculation
# q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='asdlkfjldakj;sdlkj')  # this line yields an error

# other parameters can also be set on initialization. Below is the full list
# (A description of all parameters can be found in another example)
q = multiplet.Calculation(element='Cu', 
                          charge='2+', 
                          symmetry='D4h', 
                          experiment='RIXS', 
                          edge='L2,3-M4,5 (2p3d)',
                          polarization='Isotropic', 
                          nPsis=None, 
                          nPsisAuto=1, 
                          k=(0, 0, 1), 
                          denseBorder=2000, 
                          magneticField=0.002, 
                          magneticFieldOrientation=[0, 0, 1], 
                          temperature=10, 
                          tth=130, 
                          R=[], 
                          E=None, 
                          xMin=None, 
                          xMax=None, 
                          xNPoints=None, 
                          gamma1=None, 
                          gamma2=None)
# These are not the only calculation parameters available
# The full list of parameters available can be found by 
print(q.get_attrs())
# %%

# once the calculation object is initialized, one can tweak the calculation 
# parameters. A description of all parameters can be found in another example
q.E = 932.7
q.gamma1 = 0.1

# some parameters are read only
# For instance, the pre-defined resonance energy for a given edge
print(q.resonance)
# %%

# %% ===================================================================== %% #
# %% ====================== 4) run calculation ========================= %% #
# %% ===================================================================== %% #
ss, out = q.run()
# `ss` can be a br.Spectrum or a dictionary depending on the type of calculation
# and `out` is the output text
# since `ss` can be a br.Spectrum, one can save it to a file using ss.save(filepath)


# %% ===================================================================== %% #
# %% ==================== 5) Check available elements ==================== %% #
# %% ===================================================================== %% #
# You can check the settings, if you want to check all the elemets available
print(multiplet.settings.ELEMENTS)

# for each element, you can check the available charges, symmetries, experiments, and edges
element = 'V'
print(multiplet.settings.PARAMETERS['elements'][element]['charges'].keys())

charge = '3+'
print(multiplet.settings.PARAMETERS['elements'][element]['charges'][charge]['symmetries'].keys())

symmetry = 'D4h'
print(multiplet.settings.PARAMETERS['elements'][element]['charges'][charge]['symmetries'][symmetry]['experiments'].keys())

experiment = 'RIXS'
print(multiplet.settings.PARAMETERS['elements'][element]['charges'][charge]['symmetries'][symmetry]['experiments'][experiment]['edges'].keys())

edge = 'L2,3-M4,5 (2p3d)'
print(multiplet.settings.PARAMETERS['elements'][element]['charges'][charge]['symmetries'][symmetry]['experiments'][experiment]['edges'][edge].keys())

# one can check the pre-defined resonance energy for a given edge
print(multiplet.settings.PARAMETERS['elements'][element]['charges'][charge]['symmetries'][symmetry]['experiments'][experiment]['edges'][edge]['axes'][0][4])
# however, note that checking these values (like the resonance energy) is not 
# necessary since the calculation object will automatically retrieve them when 
# the calculation is initialized.
# %%