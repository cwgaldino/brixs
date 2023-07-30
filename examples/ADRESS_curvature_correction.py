#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""ADRESS curvature correction"""

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
scan = 16

# %% step-by-step curvature correction =========================================

# initiate figure
fig, axes = br.subplots(3, 6, sharey='row')
br.maximize()

# load photon events
pes = ADRESS.raw(folderpath, prefix, scan, type_='pe')
for ccd in (1, 2, 3):
    pe = pes[ccd-1]

    # binning
    im = pe.binning(nbins=(6000, 20))

    # Calculate shifts
    im.calculate_roll()

    # fit shifts
    s = br.Spectrum(x=im.x_centers, y=im.calculated_shift)
    popt, model, R2 = s.polyfit(2)
    x = np.linspace(min(pe.x), max(pe.x), 100)
    fit = br.Spectrum(x=x, y=model(x))

    # plot photon events (raw)
    pe.plot(axes[0+(ccd-1)*6], color='black')

    # plot reduced image
    pos = im.plot(axes[2+(ccd-1)*6])

    # plot photon events (raw) with vertical bins
    pe.plot(axes[1+(ccd-1)*6], color='black')
    br.vlines(ax=axes[1+(ccd-1)*6], x=pos.x_edges, color='red', lw=.5)

    # plot horizontal integration of each vertical bin
    cols = im.columns
    cols.switch()
    cols.flip()
    cols.plot(axes[3+(ccd-1)*6])

    # get max and min y (makes plot nicer)
    cols2 = im.columns
    cols2[0].fit_peak()
    ymin = cols2[0].peaks['c'].value - cols2[0].peaks['w'].value*5

    cols2[-1].fit_peak()
    ymax = cols2[-1].peaks['c'].value + cols2[-1].peaks['w'].value*5

    # plot fitting
    offset = cols[0].y[np.argmin(cols[0].x)]
    pe.plot(axes[4+(ccd-1)*6], color='black')
    fit.plot(axes[4+(ccd-1)*6], factor=-1, offset=offset, color='red')

    # set shifts
    pe.plot(axes[5+(ccd-1)*6], color='black')

    pe.set_shift(p=popt)
    pe.plot(axes[5+(ccd-1)*6], color='red')

    # set y lim
    for i in range(6):
        axes[i+(ccd-1)*6].set_ylim(ymin, ymax)
        axes[i+(ccd-1)*6].set_xlabel(f'x pixels (subpixel)')

    axes[0+(ccd-1)*6].set_ylabel(f'ccd {ccd}: y pixels (subpixel)')

fig.subplots_adjust(top=0.99, bottom=0.05, left=0.05, right=0.99)


# %% quick curvature correction ================================================
pes = ADRESS.raw(folderpath, prefix, scan, type_='pe')

# %% plot before
fig, axes = br.subplots(1, 3, sharey=True)
for i, pe in enumerate(pes):
    pe.plot(ax=axes[i])

# %% fix curvature
for i, pe in enumerate(pes):
    popt, model = pe.fix_curvature(nbins=(6000, 20))
    print('curvature'+ str(i) + ' = [' + str(popt[0]) + ', ' + str(popt[1]) + ', ' + str(popt[2]) + ']')
print('curvature  = [curvature0, curvature1, curvature2]')

# curvature0 = [2.8024508424253414e-06, -0.018400753443323096, 0.7201982228298014]
# curvature1 = [3.769963633262668e-06, -0.01463334137482744, 0.3192278992936906]
# curvature2 = [2.8209855702191653e-06, -0.004759568633407643, 0.037546422875370744]
# curvature  = [curvature0, curvature1, curvature2]

# %% plot after
fig, axes = br.subplots(1, 3, sharey=True)
for i, pe in enumerate(pes):
    pe.plot(ax=axes[i])
