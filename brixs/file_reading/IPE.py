#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from IPE beamline - Sirius.

Last edited: Felipe Custódio 08-2023
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import numpy as np
from collections.abc import Iterable
from PIL import Image
import glob

# %% ------------------------- Special Imports ---------------------------- %% #
import brixs as br
import h5py

# %% IPE beamline - SIRIUS - Brazil ============================================
def readSPE(folderpath, prefix, ccd=0, x_min=0, x_max=3300):
    x = list()
    y = list()
    if isinstance(prefix, list):
        for p in prefix:
            for name in glob.glob(folderpath+p+'*.h5'):
                file = h5py.File(name, 'r')
                data = file['entry']['data']['data'][:]
                x.extend(data[:,2])
                y.extend(data[:,4])
                file.close()

    elif isinstance(prefix, str):
        for name in glob.glob(folderpath+prefix+'*.h5'):
            file = h5py.File(name, 'r')
            data = file['entry']['data']['data'][:]
            x.extend(data[:,2])
            y.extend(data[:,4])
            file.close()
    
    else:
        print("prefix must be str or list")

    pe = br.PhotonEvents(x, y)
    pe.ylim = [10, 1600]
    
    if ccd==1:
        x_max = 1645
        pe.label='CCD1'
    if ccd==2:
        x_min = 1655
        pe.label='CCD2'
    
    pe.xlim = [x_min, x_max]    
    
    return pe


#
