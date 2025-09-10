#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""The assumed experimental geometry is given below:

###############################################################################
######################### Top view (lab coordinates) ##########################
###############################################################################

                       .
                     .   .
LH y (010)         .       .
    |            .           .
────O────────  .   octahedron  . ────> x (100)
 LV z (001)      .           .
                   .       .
                     .   .
                       .

kin = (1, 0, 0)
LHin = (0, 1, 0)
LVin = (0, 0, 1)          

###############################################################################
###################### Top view (outgoing beam for RIXS) ######################
###############################################################################
For RIXS, the outgoing beam (kout) is then defined based on the tth angle. See below

      kout
        .
         .             
          .
           .    tth
            .        
kin ─────── ||-------------  

By default, the outgoing beam is defined by rotating kin, LVin, LHin around
z by the amount defined in tth. 

kout  = kin*Rz(tth)
LHout = LHin*Rz(tth)
LHout = LHin*Rz(tth)

###############################################################################
########################### octahedral environment ############################
###############################################################################
For a octahedral local environment, the default Hamiltonian H (Crystal Field, 
Exchange direction, etc...) is such that the planar ligands point to the x and 
y coordinates and the apical ligands are laid along the z axis. 

Here, I am calling the coordinate system of the Hamiltonian the cf coordinate system. 
If you want H to be rotated in a different orientation in relation to the 
incoming beam, one can you the parameter R.
    
R is a list of rotations one should do with the local environment before 
running the calculations. For instance R = [['x', 90], ] will rotate the 
a octahedra 90 degrees around the 'x' axis (100). This will make the apical 
ligands to point along the 'y' axis (010). 

###############################################################################
############################## Plotting geometry ##############################
###############################################################################
One can double check the experimental geometry by using the method 

q.plot_geometry()

where the left panel must match you experimental setup and the right panel must
match the left one.

WARNING 1: As for right now, the only "shape" of local environment implemented in
the plotting function is octahedral, however, the calculation of tetrahedral 
environments (and others) is working fine. It is just the plotting that will 
always show an octahedron.

WARNING 2: I do not know how the plot_geometry() function will react for a geometry
that is different from the default geometry [kin = (1, 0, 0), LHin = (0, 1, 0),
and LVin = (0, 0, 1)], but the calculation should be fine.


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
# %  ====================== Plot experimental setup ====================== %% #
# %% ===================================================================== %% #

# Initialization
# q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='XAS', edge='L2,3 (2p)')
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')

# default geometry
# The system in lab coordinates is the same as the CF coordinates
q.plot_geometry()
q.cf

# rotated octahedron
# rotation of 90 deg around x
# note how kin and kout are "rotated" in CF coordinates
q.R = [['x', 90], ]
q.plot_geometry()
q.cf

# composed rotations
# rotation of 90 deg around x and a rotation of 45 deg around y
# note how the incoming beam "enters" the octahedron through an edge and not a 
# corner like the previous plot
q.R = [['x', 90], ['y', 45]]
q.plot_geometry()
q.cf

# manipulator motors
# simulating a motor movement causing a rotation around z
th = 10
q.R = [['x', 90], ['y', 45], ['z', th]]
q.plot_geometry()
q.cf
