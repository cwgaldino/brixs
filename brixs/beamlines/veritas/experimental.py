#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Experimental functions for VERITAS beamline at MAX-IV"""

# %% ------------------------- Standard Imports --------------------------- %% #
# import matplotlib.pyplot as plt
# from itertools import compress
# import numpy as np
# import copy

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br
# from brixs.beamlines.veritas.core import read
# %%

# %% ======================== spacial support ============================= %% #
def _plot_a2scan(self):
    assert self.command.startswith('a2scan '), 'this is not a a2scan'

    # save initial x
    x0 = self.x

    # figure
    br.figure()

    # plot with x = 1 (ax1)
    self.select_x_for_a2scan(1, verbose=False)
    ax1 = self.plot().axes

    # plot with x = 2 (ax2)
    ax2 = ax1.twiny()
    self.select_x_for_a2scan(x=2, verbose=False)
    self.plot(ax=ax2, color='red')

    # linear functions for converting to motor 1 to motor 2
    # x1 = (TEY.a2scan_x[1][-1] - TEY.a2scan_x[1][0])/(TEY.a2scan_x[0][-1] - TEY.a2scan_x[0][0])*pt + TEY.a2scan_x[1][0]
    # x2 = (TEY.a2scan_x[2][-1] - TEY.a2scan_x[2][0])/(TEY.a2scan_x[0][-1] - TEY.a2scan_x[0][0])*pt + TEY.a2scan_x[2][0]
    def x1_to_pt(x):
        return (x - self.a2scan_x[1][0])*(self.a2scan_x[0][-1] - self.a2scan_x[0][0])/(self.a2scan_x[1][-1] - self.a2scan_x[1][0])
    def pt_to_x2(pt):
        return (self.a2scan_x[2][-1] - self.a2scan_x[2][0])/(self.a2scan_x[0][-1] - self.a2scan_x[0][0])*pt + self.a2scan_x[2][0]
    def x1_to_x2(x):
        return pt_to_x2(x1_to_pt(x))

    # set ax2 ticks
    ax2.set_xlim([x1_to_x2(x) for x in ax1.get_xlim()])
    ax2.spines['top'].set_color('red') 

    # labels
    ax1.set_xlabel(self.a2scan_labels[1])
    ax2.set_xlabel(self.a2scan_labels[2])

    return ax1, ax2
br.Spectrum.plot_a2scan = _plot_a2scan

def _select_x_for_a2scan(self, x, verbose=True):
    """
    x = int or motor name
    """
    if isinstance(x, str):
        assert x in self.a2scan_labels, f'x=`{x}` is not a valid motor name. Available options: {self._a2scan_x_labels}'
        x = self.a2scan_labels.index(x)
    elif isinstance(x, int):
        assert x >= 0 and x <= 2, f'x must be a int between 0 and 2 or a motor name, not `{x}`'
    else:
        raise ValueError(f'x must be a int or a motor name, not type `{type(x)}`')
    self.x = self.a2scan_x[x]
    if verbose: print(f'x selected: {self.a2scan_labels[x]}')
br.Spectrum.select_x_for_a2scan = _select_x_for_a2scan
# %%



