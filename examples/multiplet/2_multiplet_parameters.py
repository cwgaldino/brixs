#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""List of all calculation parameters:

###############################################################################
################### Primary parameters (cannot be modified) ###################
###############################################################################

element (string, optional): transition metals, lanthanides and actinides.
    default is 'Ni'.
charge (string, optional): suitable oxidation state of the element as
    '1+', '2+' and so on. default depends on the element.
symmetry (string, optional): local symmetry. Possible values are 'Oh',
    'D4h', 'Td', 'C3v', and 'D3h'. default depends on the element.
experiment (string, optional): experiment to simulate spectrum. Possible
    values are 'XAS', 'XES', 'XPS', and 'RIXS'. default is 'XAS'.
edge (string, optional): investigated edge.
    * For 'XAS' and 'XPS', possible values are 'K (1s)', 'L1 (2s)', 'L2,3 (2p)', 'M1 (3s)', and 'M2,3 (3p)'.
    * For 'XES', possible values are 'Ka (1s2p)' and 'Kb (1s3p)'.
    * For 'RIXS', possible values are 'K-L2,3 (1s2p)', 'K-M2,3 (1s3p)', and 'L2,3-M4,5 (2p3d)'
    default is depends on the experiment type and element.

yLabel, xLabel (string, read only): Axis labels.
nElectrons (number, read only): number of electrons
configurations (string, read only): Initial, Intermediate, and final configurations. Not always accurate I think.
block (string, read only): d block, f block, ...
resonance (number, read only): Energy value of the resonance.

templatePath (Path, read only): filepath for lua template.
templateName (str, read only): filename for lua template

###############################################################################
################################# Polarization ################################
###############################################################################

polarization (list, optional): type of spectrum to calculate.
    * For 'XAS', possible values are `isotropic`, `linear`, `circular`.
    * For 'XES' and 'XPS', possible value is `isotropic`.
    * For 'RIXS', possible values are `isotropic`, `linear`, `circular`.
    
###############################################################################
################################# Calculation #################################
###############################################################################

verbosity (string, optional): Verbosity for Quanty. By default it is set to minimal verbosity.
denseBorder (int, optional): Number of determinants in the bases where 
    Quanty switch from dense methods to sparse methods. The higher this 
    number, the more precise are the calculations, but also more expensive. 
    Default is 2000, which should be enough. If in doubt, one can try 
    this to see how the final results depends on it. Final result 
    shouldn't improve much for values above 2000.

E (number, optional): Only for RIXS. Energy of the incident photons. Default is
    the resonance value.

xMin, xMax, xNPoints (number, optional): minimum, maximum, and number of points.
    * For 'XAS', incident photon (in eV).
    * For 'XES', emission photon energy (in eV)
    * For 'XPS', biding energy (in eV)
    * For 'RIXS', energy transfer (in eV).
    If None, a suitable value depending on the element will
    be chosen. default is None.    

nPsis (number, optional): number of wavefunctions (eignenvalues of the 
    initial hamiltonian H_i) to br considered for the
    the ground state. The occupation of each wavefunction is defined by 
    the temperature using a boltzmann distribution. The maximum number 
    of wavefunction can be checked by q.nPsisMax. Note that q.nPsiAuto = 1
    overwrites nPsis, e.g., nPsis does not matter if nPsiAuto=1. 
    If None, nPsiAuto will be set to 1.
nPsisMax (int, read only): maximum possible number of wavefunctions.
nPsisAuto (int, optional): If 1, a suitable value will be picked for 
    nPsis, e.g., It will try and get all eigenvalues of H_i with lowest energy.
    Default is 1. Note that. for very high temperatures, the solution might not
    converge and it is recommended to set nPsisAuto = 0 and nPsis to the 
    maximum value.

nConfigurations (int, optional): number of configurations. If None,
    a suitable value will be calculated. default is None. If LMCT or
    MLCT are False in hamiltonianState, then the nConfigurations is 
    always 1 (there is only one final configuration). If LMCT or MLCT 
    are True, then there could be multiple final states where the ligand
    lends the Metal a electron (or more) and vice-versa. The maximum number
    of final states depends on how many holes the Metal or the ligand has.
    To check the maximum number of final states allowed, check q.nConfigurationsMax.
    Example: for XAS Cu2+ (d9) the initial state is p6d9 and the final 
    state is p5d10. If LMCT is `on`, then we have two possible initial 
    states, p6 d9 L10 and p6 d10 L9. nConfigurations is basically, the 
    maximum number of electrons to be shared between M and L.
nConfigurationsMax (int, read only): Maximum number of final configuration 
    allowed.

###############################################################################
############################# Core-hole lifetime ##############################
###############################################################################

gamma1 (number, optional): 
    * For 'RIXS', 'XAS', 'XES', and 'XPS', core-hole lifetime (defines lorentzian broadening - FWHM
    for RIXS it is the lorentzian broadening of the `photon incident energy`
gamma2 (tuple or list, optional): core-hole lifetime. 
    * For 'XAS', 'XES', and 'XPS', n.a.
    * For 'RIXS', lifetime of the intermediate state (defines lorentzian broadening of the `energy loss`)

###############################################################################
################################# Experiment ##################################
###############################################################################

temperature (number, optional): temperature in Kelvin. default is 10 K.
    If temperature is zero, nPsi is set to 1 and nPsiAuto to 0.
magneticField (number, optional): Magnetic field value in Tesla. If zero
    or None, a very small magnetic field (0.002 T, this is gives the smallest 
    float possible when describing muB Bohr magneton (T) in units of electron 
    Volt) will be used to have nice expected values for observables. To force a
    zero magnetic field you have to turn it off via the q.hamiltonianState. 
    Default is 0.002 T. This can be overwritten by changing the hamiltonian data
    directly.
magneticFieldOrientation (vector, optional): Direction of the magnetic field in 
    relation to the CF coordinate system.
    Default is (001). This can be overwritten by changing the hamiltonian data
    directly. Use q.plot_geometry() to double check the orientation of the 
    magnetic field.

###############################################################################
################################## Geometry ###################################
###############################################################################

k (vector, optional): quantization axis for calculating expected values. 
    It does not influence the calculation, just the table printed at the end.
    Default is (001).
R (list, optional): list of rotations one should do with the local environment before 
    running the calculations. For instance R = [['x', 90], ] will rotate the 
    a octahedra 90 degrees around the 'x' axis (100). This will make the apical 
    ligands to point along the 'y' axis (010). Example: R=[['x', 90], ['y', 45]]. 
    One can double check the experimental geometry by using the method 
    q.plot_geometry(). Default is []. 
tth (number, optional): Only for RIXS. Number from 0 to 180 indicating the 
    2theta scattering angle. This defined the outgoing polarization vectors.
    Default is 130. Note that the outgoing polarization vector can be defined 
    directly by kout, LVout, LHout.

kin, LVin, LHin (list, optional): incoming vectors in lab coordinates. 
    By default, the experimental geometry is defined as 
        kin = (1, 0, 0)
        LHin = (0, 1, 0)
        LVin = (0, 0, 1)
kout, LVout, LHout (list, read only): Only for RIXS. Outgoing vectors in lab 
    coordinates. Note that these are defined by rotating kin, LVin, LHin around
    z by the amount defined in tth. 
kin_cf, LVin_cf, LHin_cf (list, read only): incoming vectors in cf coordinates. 
kout_cf, LVout_cf, LHout_cf (list, read only): Only for RIXS. Outgoing vectors 
    in cf coordinates.
cf (list, read only): cf coordinate system in terms of lab coordinate system.
    
###############################################################################
################################# Hamiltonian #################################
###############################################################################
hamiltonianTerms (list, read only): a list of hamiltonian terms.
hamiltonianState (dict, optional): a dictionary that turns on or off the
    different contributions to the hamiltonian. default is such that
    'Atomic', 'Crystal Field', and 'Magnetic Field' terms will be True.
hamiltonianData (dict, optional): dictionary with values for the strength
    of each different contributions to the hamiltonian. The default value
    ise different for each element and charge stat. Values are in eV.

###############################################################################
################################# Lua script ##################################
###############################################################################
template (string): loaded template
lua_script (string): lua script to be run


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
# %  ======================= Calculation parameters ====================== %% #
# %% ===================================================================== %% #

# Initialization
# q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='XAS', edge='L2,3 (2p)')
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')

# list of all calculation attributes
q.get_attrs()

# Primary parameters (cannot be modified)
q.element         # read only
q.charge          # read only
q.symmetry        # read only
q.experiment      # read only
q.edge            # read only

q.configurations  # read only
q.block           # read only
q.nElectrons      # read only
q.resonance       # read only

q.xLabel          # read only
q.yLabel          # read only

q.templateName    # read only
q.templatePath    # read only

# set polarization type
q.polarization

# Calculation parameters
q.verbosity
q.denseBorder

q.E  # only for RIXS

q.xMin
q.xMax
q.xNPoints

q.nPsis
q.nPsisMax  # read only
q.nPsisAuto

q.nConfigurations 
q.nConfigurationsMax  # read only
  
# Core-hole lifetime (lorentzian broadening)
q.gamma1
q.gamma2  # only for RIXS

# Experiment
q.temperature
q.magneticField
q.magneticFieldOrientation

# geometrical parameters
q.k
q.R
q.tth  # only for RIXS

q.kin
q.LHin
q.LVin

q.kout     # read only
q.LHout    # read only
q.LVout    # read only

q.kin_cf   # read only
q.LVin_cf  # read only
q.LHin_cf  # read only

q.kout_cf   # read only
q.LVout_cf  # read only
q.LHout_cf  # read only

q.cf        # read only

# hamiltonian parameters 
q.hamiltonianTerms   # read only
q.hamiltonianState  
print(q.hamiltonianData)  

# Lua script (can be edited manually if necessary)
print(q.template)
q.update_lua_script()  # Use q.update_lua_script() update lua script
q.lua_script  
# %%


# %  ===================================================================== %% #
# %  ================================ Tips =============================== %% #
# %% ===================================================================== %% #
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')

# Before, editing the hamiltonian parameters, print the hamiltonianData to 
# learn which hamiltonian parameters are available
print(q.hamiltonianData)  

# It is easier to loop through Hamiltonian "states". For instance, for a
# transition metal ion in a D4h configuration, the crystal field is described 
# with three parameters Dq, Ds, and Dt. These parameters are present in the 
# initial and final hamiltonian (and Intermediate hamiltonian for RIXS). Loop
# through each state to change them quickly
Dq = 0.1
Ds = 0.08
Dt = -0.01
for h in ['Initial Hamiltonian', 'Final Hamiltonian']:
    q.hamiltonianData['Crystal Field'][h]['Dq(3d)'] = round(Dq, 5)
    q.hamiltonianData['Crystal Field'][h]['Ds(3d)'] = round(Ds, 5)
    q.hamiltonianData['Crystal Field'][h]['Dt(3d)'] = round(Dt, 5)

# For RIXS, the parameter q.resonance indicates the tabulated value to be 
# considered as the energy of the main transition for the selected edge
print(q.resonance)

# Changing kin, LHin, and LVin is not recommended as the default values 
# should be fine for all calculations
print(q.kin, q.LHin, q.LVin)
# Therefore, the current implementation of this module makes it tricky to 
# change them, because the `q` expects kin, LHin, and LVin to always be 
# perpendicular, however, this condition may not be satisfied while changing 
# these parameters. Note how the following command gives an error
# q.kin = [0, 1, 0]  # this line yields an error
# before we change kin, LHin, and LVin, let's print q.kout to see what we have
print(q.kout)  # [-0.6427876096865394, 0.766044443118978, 0.0]
# this parameters is defined from q.kin and q.tth, where kout = kin*Rz(tth)
# Finally, to change kin, LHin, and LVin, one has to change the variable 
# "internally"
q._kin  = [0, 2, 0]  # note the `_` before kin
q._LHin = [0, 0, 1]  # note the `_` before _LHin
q._LVin = [1, 0, 0]  # note the `_` before _LVin
# not how q.kout is still the same as before because although kin, LHin, and 
# LVin changed, all the verification and updating code was overwritten
print(q.kout)
# now we re-assign in, LHin, and LVin but without the `_` allowing `q` to run all 
# internal functions
q.kin   = q.kin
q.LHin  = q.LHin
q.LVin  = q.LVin
print(q.kin, q.LHin, q.LVin)
print(q.kout)  # [-1.532088886237956, -1.2855752193730787, 0.0]
# note how q.kin was normalized and that kout changed
# also note that changing kin, LHin, and LVin makes the left panel of 
# q.plot_geometry() weird, but the left panel is fine and therefore the 
# calculation is fine. This will be fixed someday.
q.plot_geometry()

# The magnetic field parameters is a shortcut so one does not have to change the
# hamiltonianData directly
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')
q.R = [['x', 90], ]
print(q.magneticField)
print(q.magneticFieldOrientation)
print(q.hamiltonianData['Magnetic Field'])
# if q.magneticField, then hamiltonianData changes accordingly
q.magneticField = 1
print(q.hamiltonianData['Magnetic Field'])
# the magnetic field orientation is in relation to the CF environment
q.magneticFieldOrientation = [0, 0, 1]
q.plot_geometry()
q.magneticFieldOrientation = [0, 1, 0]
q.plot_geometry()
