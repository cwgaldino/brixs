#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""The brixs.beamlines.ipe module offers a python solution for data processing
and data analysis for data collected at the IPE beamline of SIRIUS.

Example data for running this example can be downloaded at this Onedrive link
https://1drv.ms/f/c/6666810d266a3cc9/EjNzJssVTWhFqAOx0-nDCxoBSNt2cvvF8zqOmNdx_zzslQ?e=BOAgx4
Note that some less important files have been removed from this example data to
make the folder lighter and easier to download.


Besides brixs requirements, brixs.beamlines.ipe module also requires h5py.

Onsite, the path of the top (main) directory is 

/ibira/lnls/beamlines/ipe/proposals/20251234/

The folder structure of the top directory is:

data
    RIXS
        SPE  --> Folders containing processed images (photon events list)
        TIF  --> Folders containing images from the XCAM detector
proc
    DATA
        _XAS
        ASCAN
        MESH
        XAS
    scripts
"""

# %% ========================== standard imports ========================= %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# %% =========================== brixs imports =========================== %% #
#### The following two lines are only necessary if you are importing brixs 
#### from a local file
# import sys
# sys.path.insert(0, 'brixs-main/')
import brixs as br
import brixs.beamlines.ipe as ipe

# %% ============================== settings ============================= %% #
# matplotlib (optional)
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()

# brixs (optional)
# br.settings.FIGURE_POSITION = (283, 567)
# br.get_window_position()
# %%

# %% ============================= folderpaths =========================== %% #
# the path to the top/main folder is going to be used constantly, so we are better
# off defining a variable for it
TOP = Path(r'C:\Users\galdin_c\github\brixsexampledata\beamlines\ipe')
# %%

# TODO

