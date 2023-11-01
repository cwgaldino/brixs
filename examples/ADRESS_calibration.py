#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""ADRESS beamlime example"""

# %% imports ===================================================================
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import brixs as br
from brixs.file_reading import ADRESS

# %% autoreload and matplotlib backend (ignore) ================================
if br.is_notebook():
    from IPython import get_ipython
    get_ipython().run_line_magic('matplotlib', 'qt5')
    get_ipython().run_line_magic('load_ext', 'autoreload')
    get_ipython().run_line_magic('autoreload', '2')
else:
    plt.ion()

# %% Initial definitions =======================================================
folderpath = Path(r'calib')
prefix = 'O_'

# %% quick calibration calculation =============================================
# calibration uses cross-correlation by default
popt, sss = ADRESS.calib(folderpath, prefix, 17, 26)

# %% print as string
for i, p in enumerate(popt):
    print('calib'+ str(i) + ' = [' + str(popt[i][0]) + ', ' + str(popt[i][1]) + ']')
print('calib  = [calib0, calib1, calib2]')

# mode cc (ev/bin)
# calib0 = [0.004556853344546129, 516.1896152709334]
# calib1 = [0.0045402249384432535, 515.6958557931846]
# calib2 = [0.0045657749610020166, 515.224302111583]
# calib  = [calib0, calib1, calib2]

# %% plot step-by-step =========================================================
fig, axes = br.subplots(2, 3, sharey='row', sharex='col')
for ccd in (0, 1, 2):
    # plot spectra
    sss[ccd].plot(ax=axes[ccd])
    axes[ccd].set_xlabel('bin')

    # fitting results
    bin = np.linspace(sss[ccd].x[0], sss[ccd].x[-1], 100)
    fit = np.polyval(popt[ccd], bin)
    
    # find peak first spectrum
    sss[ccd][0].fit_peak()
    c = sss[ccd][0].peaks['c']

    # plot fitting results
    axes[ccd+3].plot(-sss[ccd].calculated_shift+c, sss[ccd].E, marker='o', color='black')
    axes[ccd+3].plot(bin, fit, color='red')

# labels
axes[3].set_xlabel('bin')
axes[4].set_xlabel('bin')
axes[5].set_xlabel('bin')
axes[0].set_ylabel('Intensity (arb. units)')
axes[3].set_ylabel('Intensity (arb. units)')


# %% other modes ===============================================================

# %% mode peaks (ev/bin)
popt, sss = ADRESS.calib(folderpath, prefix, 17, 26, mode='peak')
calib0 = [0.004556657343031468, 516.1910040053092]
calib1 = [0.004539722713328706, 515.696903029858]
calib2 = [0.004565239596075693, 515.2262740278696]

# %% mode cc (ev/subpixel)
curvature0 = [2.8024508424253414e-06, -0.018400753443323096, 0.7201982228298014]
curvature1 = [3.769963633262668e-06, -0.01463334137482744, 0.3192278992936906]
curvature2 = [2.8209855702191653e-06, -0.004759568633407643, 0.037546422875370744]
curvature  = [curvature0, curvature1, curvature2] # see ADRESS curvature correction example
 
popt, sss = ADRESS.calib(folderpath, prefix, 17, 26, nbins=2000, curvature=curvature)
calib0_subpixel = [0.018228083877778765, 516.1699816340371]
calib1_subpixel = [0.018147251168973262, 515.696753572465]
calib2_subpixel = [0.018262961090482758, 515.2192036580997]


# %% read ADRESS files with calibration factor =================================
# elastic-line error is less than 0.1 eV for oxygen-K edge
br.figure()
for scan in (17, 18, 19, 20, 21, 22, 23, 24, 25, 26):
    s = ADRESS.read(folderpath, prefix, scan, calib=calib0)
    s.plot()
br.vlines(0)
plt.xlim(-0.5, 0.5)

