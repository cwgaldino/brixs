"""Example file for verifying the rixs processing of scans taken at i21.
"""
# %%
# default imports
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# brixs imports
import sys
sys.path.append('brixs/')
import brixs as br
import brixs.beamlines.i21 as i21

# matplotlib backend (use widget on jupyter notebook, or qt5 elsewhere)
if br.is_jupyter:
    %matplotlib widget
else:
    %matplotlib qt5
    br.settings.FIGURE_POSITION = (1500, 1400)  # home
    br.settings.FIGURE_POSITION = (60, 140) 
    # br.get_window_position()

# Enable interactive mode
plt.ion()

# Constrained layout is on by default to make figures larger
br.settings.FIGURE_CONSTRAINED_LAYOUT = True

# folderpaths
proposal = 'mm43282-1'
TOP   = Path(r'/dls/i21/data/2026/' + proposal + r'/processing')
DATA  = TOP/'..'
TMP   = TOP/'tmp'
OUT   = TOP/'out'
STORE = TOP/'store'

# finder settings
br.finder.folderpath = TMP
br.finder.verbose    = True
br.finder.search_on  = False  # <<< THIS MUST BE OFF FOR _process()
br.finder.save_on    = False  # <<< THIS MUST BE OFF FOR _process()

# default i21 parameters
i21.settings.FOLDERPATH = DATA
i21.settings.SLOPE = 0.03021685
i21.settings.CALIB = 0.00243995

# %% 0. Processing
scan = 472641
dark = 472610
data = i21._process(scan, dark)

# %% 1. Cosmic ray removal from dark image
number_of_frames = len(data['ds0'])
for i in range(number_of_frames):
    fig, axes = br.subplots(1, 2, sharex=True, sharey=True)
    fig.suptitle(f'Cosmic rays removal from dark images: frame {i+1} of {number_of_frames}' + '\nds0, ds1, dpes1\n' + data['text0'], fontsize=8)
    data['ds0'][0].plot(ax=axes[0])
    data['ds1'][0].plot(ax=axes[1])
    for ax in (axes[0], axes[1]):
        data['dpes1'][0].plot(ax=ax, marker='o', s=20, edgecolor='white', facecolor=None)
    axes[0].set_title('Before')
    axes[1].set_title('After')
    
# %% 2. Averaged out and smooth dark image
fig, axes = br.subplots(1, 4, figsize=(28, 12), sharey=True)
fig.suptitle(f'Smooth dark image ' + '\nd0, d2, d3\n' + data['text1'], fontsize=8)
data['d0'].plot(ax=axes[0])
data['d2'].plot(ax=axes[1])
data['d3'].plot(ax=axes[2])
data['d0'].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[3])
data['d2'].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[3])
data['d3'].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[3])
axes[0].set_title('1. Raw average image', fontsize=8)
axes[1].set_title('2. After cosmic rays removal', fontsize=8)
axes[2].set_title('3. After medial and gaussian smoothing', fontsize=8)
axes[3].set_title('4. Integrated rows', fontsize=8)

# %% 3. Cosmic rays removal
number_of_frames = len(data['ims0'])
for i in range(number_of_frames):
    fig, axes = br.subplots(1, 2, sharex=True, sharey=True)
    fig.suptitle(f'Cosmic rays removal: frame {i+1} of {number_of_frames}' + '\nims0, ims1, pes1\n' + data['text2'], fontsize=8)
    data['ims0'][0].plot(ax=axes[0])
    data['ims1'][0].plot(ax=axes[1])
    for ax in (axes[0], axes[1]):
        data['pes1'][0].plot(ax=ax, marker='o', s=20, edgecolor='white', facecolor=None)
    axes[0].set_title('Before')
    axes[1].set_title('After')

# %% 4. Dark image subtraction
number_of_frames = len(data['ims2'])
for i in range(number_of_frames):
    fig, axes = br.subplots(1, 3, figsize=(28, 12), layout='constrained')
    fig.suptitle(f'Dark image subtraction: frame {i+1} of {number_of_frames}' + '\nims1, ims2\n' + data['text3'], fontsize=8)
    data['ims1'][i].integrated_rows_vs_y_centers().plot(ax=axes[0])
    data['d3'].set_factor(data['dark_factors'][i]).set_offset(data['dark_offsets'][i]).integrated_rows_vs_y_centers().plot(ax=axes[0])    
    data['ims2'][i].integrated_rows_vs_y_centers().plot(ax=axes[1])
    data['ims2'][i].plot(ax=axes[2])
    br.axvlines(data['dark_offset_limits'], ax=axes[0], color='black', ls='--', lw=0.5)
    br.axvlines(data['dark_factor_limits'], ax=axes[0], color='red', ls='--', lw=0.5)
    axes[0].set_title('1. Raw spectrum', fontsize=8)
    axes[1].set_title('2. Spectrum after subtraction', fontsize=8)
    axes[2].set_title('3. Frame after subtraction', fontsize=8)

# %% 5. Curvature correction
fig, axes = br.subplots(1, 2, sharex=True, sharey=True)
fig.suptitle(f'Curvature correction ' + '\nims2, ims3\n' + data['text4'], fontsize=8)
data['im2'].plot(ax=axes[0])
data['im3'].plot(ax=axes[1])
axes[0].set_title('Before')
axes[1].set_title(f'After')
y = data['s0'].get_x_where_y_is_max()
for ax in axes:
    br.axhlines(y, ax=ax, color='red', ls='--', lw=0.5)
