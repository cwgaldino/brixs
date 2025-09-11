#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exploratory analysis for CrTe2 beamtime at ESRF HC-6098"""

# %%
cd C:\Users\galdin_c\Documents\work\CrSBr\ares\2025_09_10_ESRF_EDA
ipython
# %%

# %% ========================== standard imports ========================= %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np


# %% =========================== brixs imports =========================== %% #
import brixs as br
import brixs.beamlines.id32 as id32

# %% ============================== settings ============================= %% #
# matplotlib
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()

# %% ============================= folderpaths ============================ %% #
TOP = Path(r'C:\Users\galdin_c\Documents\work\CrSBr\ares\2025_09_10_ESRF_EDA')

OUT    = TOP/'out'   # output
TMP    = TOP/'tmp'   # temporary files (can be deleted)
STORE  = TOP/'store' # storage (not output file, cannot be deleted)




RAW = Path(r'D:\esrf_CrTe2_data')
TOP  = Path(r'/data/visitor/hc6098/id32/20250408')

XAS  = TOP/'RAW_DATA'
RIXS = TOP/'PROCESSED_DATA'
SAMPLE_ALIGN = RIXS/'align/align_0001/Online_analysis_align_0001_dhyana95v2.spec'
XAS1 = XAS/'CrTe2/CrTe2_0001/CrTe2_0001.h5'
# %%



list_of_scans

