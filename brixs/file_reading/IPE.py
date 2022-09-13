#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for file handling."""

from pathlib import Path
import numpy as np
import os
from PIL import Image

import brixs as br


# %% IPE beamline - SIRIUS - Brazil ============================================
def read(filepath):
    """Read files from IPE beamline at SIRIUS.

    Last updated: 26/12/2021 by Carlos Galdino
    """
    im = Image.open(filepath)
    return br.Image(data=np.array(im).astype('float64'))












#
