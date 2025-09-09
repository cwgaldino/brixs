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
    for RIXS it is the  lorentzian broadening of the `photon incident energy`
gamma2 (tuple or list, optional): core-hole lifetime. 
    * For 'XAS', 'XES', and 'XPS', n.a.
    * For 'RIXS', lifetime of the intermediate state (defines lorentzian broadening of the `energy loss`)

###############################################################################
################################# Experiment ##################################
###############################################################################

temperature (number, optional): temperature in Kelvin. default is 10 K.
    If temperature is zero, nPsi is set to 1 and nPsiAuto to 0.
magneticField (number, optional): Magnetic field value in Tesla. If zero
    or None, a very small magnetic field (0.002 T) will be used to
    have nice expected values for observables. To force a zero magnetic
    field you have to turn it off via the q.hamiltonianState. default
    is 0.002 T. This can be overwritten by changing the hamiltonian data
    directly.
magneticFieldOrientation (vector, optional): Direction of the magnetic field.
    Default is (001). This can be overwritten by changing the hamiltonian data
    directly.

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
q.element
q.charge
q.symmetry
q.experiment
q.edge

q.configurations
q.block
q.nElectrons
q.resonance

q.xLabel
q.yLabel

q.templateName
q.templatePath

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
q.nPsisMax
q.nPsisAuto

q.nConfigurations 
q.nConfigurationsMax
  
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

q.kout
q.LHout
q.LVout

q.kin_cf
q.LVin_cf
q.LHin_cf

q.kout_cf
q.LVout_cf
q.LHout_cf

q.cf

# hamiltonian parameters 
q.hamiltonianTerms 
q.hamiltonianState  
print(q.hamiltonianData)  

# Lua script (can be edited manually if necessary)
print(q.template)
q.update_lua_script()  # Use q.update_lua_script() update lua script
q.lua_script  

