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

kin  = (1, 0, 0)
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
kin ─────── ||------------- ────> x (100)

By default, the outgoing beam is defined by rotating kin, LVin, LHin around
z by the amount defined in tth. 

kout  = kin*Rz(tth)
LVout = LVin*Rz(tth)
LHout = LHin*Rz(tth)

###############################################################################
########################### octahedral environment ############################
###############################################################################
For a octahedral local environment, the default Hamiltonian H (Crystal Field, 
Exchange direction, etc...) is such that the planar ligands point to the x and 
y coordinates and the apical ligands are laid along the z axis. 

By definition, the coordinate system of the Hamiltonian H is called the cf coordinate system. 
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

# %% ========================== Standard imports ========================= %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# %% ============================ brixs imports ========================== %% #
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

# %% ===================================================================== %% #
# %% ====================== Plot experimental setup ====================== %% #
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

# WARNING!!!!
# Be careful when appending values to q.R.
# Appending to a list in Python does not call the property's setter
# method. Therefore, the internal geometry calculations are not updated
# automatically and must be triggered explicitly.
q.R = [['x', 90], ['z', 10]]
q.R.append(['z', th])
q._update_geometry()  # Explicitly call the geometry update method.
# %%

# %% ===================================================================== %% #
# %% ================== Geometry of the magnetic field =================== %% #
# %% ===================================================================== %% #
# Initialization
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')

# The magnetic field direction is defined by the hamiltonian
print(q.hamiltonianData['Magnetic Field'])
# By default, a small magnetic field is applied along the z axis (001) in relation to the CF environment
# {'Bx': 0.0
# 'By': 0.0
# 'Bz': 1.1576e-07}
# This is to break the degeneracy of the ground state,
# so make the final wavefuctions have nice quantum numbers.

# One can see the direction of the magnetic field in relation to the 
# CF environment by plotting the geometry
q.plot_geometry()

# One can change the direction of the magnetic field by using the following parameters
q.magneticField_coordinate_system = 'cf'    # default is 'cf'
q.magneticField = 1                         # in Tesla
q.magneticFieldOrientation = [0, 1, 0]      # in terms of lab coordinates
print(q.hamiltonianData['Magnetic Field'])  # converted from Tesla to eV
q.plot_geometry()

# One can also define the magnetic field in terms of lab coordinates
# note that, by default, the coordinate sytem of the magnetic field is given 
# in terms of the cf environment (q.magneticField_coordinate_system = 'cf')
q.magneticField_coordinate_system = 'lab'
q.magneticField = 1                     # in Tesla
q.magneticFieldOrientation = [0, 1, 0]  # in terms of lab coordinates
q.plot_geometry()

# if the magnetic field is given in terms of the lab coordinates, than q.R 
# will also affect the direction of the magnetic field as the script needs to  
# convert the magnetic field from lab to cf coordinates
q.magneticField_coordinate_system = 'lab'
q.magneticField = 1                     # in Tesla
q.magneticFieldOrientation = [0, 1, 0]  # in terms of lab coordinates
print(q.hamiltonianData['Magnetic Field']) 
q.R = [['x', 90], ]
print(q.hamiltonianData['Magnetic Field']) 

# Alternatively, one can change the magnetic field by changing the 
# hamiltonianData directly.
# However, this is slightly more incovenient as one has to give the values in 
# eV and x, y, and z components of the magnetic field separetly in terms of the
# cf enviroment. Keep in mind that, any changes to
# q.magneticField_coordinate_system, q.magneticField, and q.magneticFieldOrientation 
# (also q.R if q.magneticField_coordinate_system='lab') will be reflected in the
# hamiltonianData, however the inverse is not true. The calculation will run
# ok though because it ultimately reads the hamiltonianData when setting up the 
# quanty input file.
for _ in ['Initial Hamiltonian', 'Intermediate Hamiltonian', 'Final Hamiltonian']:
    q.hamiltonianData['Magnetic Field'][_]['Bx'] = 0.0
    q.hamiltonianData['Magnetic Field'][_]['By'] = 1.1576e-07
    q.hamiltonianData['Magnetic Field'][_]['Bz'] = 0.0
print(q.magneticFieldOrientation)  # does not match with hamiltonianData