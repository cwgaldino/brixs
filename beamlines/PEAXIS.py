#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for file handling."""

from pathlib import Path
import numpy as np
import os

import brixs as br

# %% PEAXIS beamline - HZB - Germany ===========================================
def _ReadAndor(fname, dimensions = (2048, 2048), byte_size=4, data_type='c'):
    """Reads the .sif file produced by the Andor iKon-L CCD camera.
    The dimensions of the pixel array can be changed if necessary.
    """
    size = os.path.getsize(fname)
    source = open(fname, 'rb')
    nrows = dimensions[-1]
    header, data = [],[]
    bindata = None
    total_length = 0
    header = np.fromfile(source, np.byte, size - 4*dimensions[0]*dimensions[1] -2*4, '') # 2772 is close
    # print("Headersize",  size - 4*dimensions[0]*dimensions[1] -2*4)
    data = np.fromfile(source, np.float32, dimensions[0]*dimensions[1], '')
    header = bytes(header)
    lines = header.split(b'\n')
    header = []
    for n, line in enumerate(lines):
        try:
            header.append(line.decode('ascii'))
        except UnicodeDecodeError:
            print('header lines skipped:', n+1, 'with length:', len(line))
    source.close()
    return np.array(data).reshape(dimensions), header

def read(filepath):
    """Read files from PEAXIS beamline at HZB.

    Example:

        >>> import brixs as br
        >>> im = br.read_PEAXIS(filepath)
        >>> print(im.header)

    Args:
        filepath (str or Path object): filepath.

    Returns:
        :py:class:`brixs.Image`

    Last updated: 04/04/2022 by Carlos Galdino
    """
    filepath = Path(filepath)
    data, header = _ReadAndor(filepath)
    im = br.Image(data)
    im.header = header
    return im









#
