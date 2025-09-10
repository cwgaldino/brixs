#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Multiplet calculation

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
# %  =========================== XAS isotropic =========================== %% #
# %% ===================================================================== %% #

# Initialization
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='XAS', edge='L2,3 (2p)')
q.polarization = 'isotropic'

# parameters
F2 = 0.8
F4 = 0.8
F2_2p = 0.8
G1_2p = 0.8
G3_2p = 0.8

Dq = 0.1
Ds = 0.05
Dt = 0.01

so    = 0.8
so_2p = 0.9

# core-hole
q.gamma1 = 0.6

#####################
# Hamiltonian terms #
#####################
q.hamiltonianState['Atomic'] = True
q.hamiltonianState['Crystal Field'] = True
q.hamiltonianState['3d-Ligands Hybridization (LMCT)'] = False
q.hamiltonianState['3d-Ligands Hybridization (MLCT)'] = False
q.hamiltonianState['Magnetic Field'] = True
q.hamiltonianState['Exchange Field'] = False

##########
# Atomic #
##########
for h in ['Initial Hamiltonian', 'Final Hamiltonian']:
    q.hamiltonianData['Atomic'][h]['F2(3d,3d)'][1] = F2
    q.hamiltonianData['Atomic'][h]['F4(3d,3d)'][1] = F4
h = 'Final Hamiltonian'
q.hamiltonianData['Atomic'][h]['F2(2p,3d)'][1] = F2_2p
q.hamiltonianData['Atomic'][h]['G1(2p,3d)'][1] = G1_2p
q.hamiltonianData['Atomic'][h]['G3(2p,3d)'][1] = G3_2p

##############
# Spin-orbit #
##############
for h in ['Initial Hamiltonian', 'Final Hamiltonian']:
    q.hamiltonianData['Atomic'][h]['ζ(3d)'][1] = so
h = 'Final Hamiltonian'
q.hamiltonianData['Atomic'][h]['ζ(2p)'][1] = so_2p

#################
# Crystal Field #
#################
for h in ['Initial Hamiltonian', 'Final Hamiltonian']:
    q.hamiltonianData['Crystal Field'][h]['Dq(3d)'] = round(Dq, 5)
    q.hamiltonianData['Crystal Field'][h]['Ds(3d)'] = round(Ds, 5)
    q.hamiltonianData['Crystal Field'][h]['Dt(3d)'] = round(Dt, 5)

# run calculation
s, out = q.run()

# plot calculated spectrum
plt.figure()
s.plot()
br.labels.xas()

# change spin-orbit
h = 'Final Hamiltonian'
q.hamiltonianData['Atomic'][h]['ζ(2p)'][1] = 1.2
s, out = q.run()
s.plot()

# broadening
plt.figure()
s.plot(label='raw', marker='o')
s.broaden(0.5, m=0).plot(label='0.5 gaussian')
s.broaden(0.5, m=1).plot(label='0.5 lorentzian')
s.broaden(0.5, m=0.5).plot(label='0.5 (lorentzian + gaussian)')
s.broaden(0.5, m=1).broaden(0.5, m=0).plot(label='0.5 lorentzian + 0.5 gaussian')
br.leg()
br.labels.xas()
# %%

# %  ===================================================================== %% #
# %  ====================== XAS Linear polarization ====================== %% #
# %% ===================================================================== %% #

# Initialization
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='XAS', edge='L2,3 (2p)')
q.polarization = 'linear'
q.R = [['x', 90], ['z', 45]]
# q.plot_geometry()

Dq = 0.1
Ds = 0.08
Dt = -0.01

#################
# Crystal Field #
#################
for h in ['Initial Hamiltonian', 'Final Hamiltonian']:
    q.hamiltonianData['Crystal Field'][h]['Dq(3d)'] = round(Dq, 5)
    q.hamiltonianData['Crystal Field'][h]['Ds(3d)'] = round(Ds, 5)
    q.hamiltonianData['Crystal Field'][h]['Dt(3d)'] = round(Dt, 5)

# run calculation
ss, out = q.run()

# plot calculated spectrum
plt.figure()
ss['LV'].plot(label='LV')
ss['LH'].plot(label='LH')
br.labels.xas()
# %%

# %  ===================================================================== %% #
# %  ====================== RIXS Linear polarization ===================== %% #
# %% ===================================================================== %% #

# Initialization
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')
q.polarization = 'linear'
q.R = [['x', 90], ['z', 45]]
q.xMin = -1
q.xMax = 3

Dq = 0.1
Ds = 0.08
Dt = 0.06

#################
# Crystal Field #
#################
for h in ['Initial Hamiltonian', 'Final Hamiltonian']:
    q.hamiltonianData['Crystal Field'][h]['Dq(3d)'] = round(Dq, 5)
    q.hamiltonianData['Crystal Field'][h]['Ds(3d)'] = round(Ds, 5)
    q.hamiltonianData['Crystal Field'][h]['Dt(3d)'] = round(Dt, 5)

# run calculation
ss, out = q.run()

# plot calculated spectrum 1 
plt.figure()
ss['vv'].plot(label='LV in, LV out')
ss['vh'].plot(label='LV in, LH out')
ss['hv'].plot(label='LH in, LV out')
ss['hh'].plot(label='LH in, LH out')
br.labels.xas()

# plot calculated spectrum 2
LV = br.Spectra([ss['vv'], ss['vh']]).calculate_average()
LH = br.Spectra([ss['hv'], ss['hh']]).calculate_average()

plt.figure()
LV.plot(label='LV in')
LH.plot(label='LH in')
br.labels.xas()
# %%

# %  ===================================================================== %% #
# %  ========================== RIXS Energy map ========================== %% #
# %% ===================================================================== %% #

# Initialization
q = multiplet.Calculation(element='Cu', charge='2+', symmetry='D4h', experiment='RIXS', edge='L2,3-M4,5 (2p3d)')
q.polarization = 'isotropic'
q.xMin = -1
q.xMax = 2

Dq = 0.1
Ds = 0.08
Dt = 0.06

q.gamma1 = 1

#################
# Crystal Field #
#################
for h in ['Initial Hamiltonian', 'Final Hamiltonian']:
    q.hamiltonianData['Crystal Field'][h]['Dq(3d)'] = round(Dq, 5)
    q.hamiltonianData['Crystal Field'][h]['Ds(3d)'] = round(Ds, 5)
    q.hamiltonianData['Crystal Field'][h]['Dt(3d)'] = round(Dt, 5)

# %% run calculation (takes a while)
energies = np.linspace(q.resonance-10, q.resonance+30, 300)  # 300 points
ss = br.Spectra()
for E in energies:
    q.E = E
    s, out = q.run()
    ss.append(s)

# get map
im = ss.stack_spectra_as_columns()
im.x_centers = energies

# XAS
_y = ss.calculate_y_sum()
xas = br.Spectrum(x=energies, y=_y/max(_y))

# plot map
br.figure()
im.plot(origin='lower')
xas.plot(color='white')
br.labels.energy_map()
# %%

# %% run same calculation but with larger gamma1 (takes a while)
q.gamma1 = 2
ss2 = br.Spectra()
for E in energies:
    q.E = E
    s, out = q.run()
    ss2.append(s)

# get map
im2 = ss2.stack_spectra_as_columns()
im2.x_centers = energies

# XAS
_y = ss2.calculate_y_sum()
xas2 = br.Spectrum(x=energies, y=_y/max(_y))

# plot map
br.figure()
im2.plot(origin='lower')
xas2.plot(color='white')
br.labels.energy_map()
# %%

# %% broadening x (takes a while)
im3 = im.broaden_x(2, m=0).broaden_y(0.1, m=0)
xas3 = im3.integrated_columns_vs_x_centers()
xas3 = xas3.set_factor(1/max(xas3))

br.figure()
im3.plot(origin='lower')
xas3.plot(color='white')
br.labels.energy_map('broadened 2.0 gaussian x')

# get broadened spectra
ss3 = im3.get_columns(max_number_of_columns=im3.shape[1])
# %%

# %% broadening y (takes a while)
im4 = im.broaden_x(2, m=0).broaden_y(0.1, m=0)
xas4 = im4.integrated_columns_vs_x_centers()
xas4 = xas4.set_factor(1/max(xas3))

br.figure()
im4.plot(origin='lower')
xas4.plot(color='white')
br.labels.energy_map('broadened 2 gaussian x and 0.1 gaussian y')

# get broadened spectra
ss4 = im4.get_columns(max_number_of_columns=im3.shape[1])
# %%