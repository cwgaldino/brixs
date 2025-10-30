#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for everyday use ---> pdf files"""

# %% ------------------------- Standard Imports -------------------------- %% #
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pathlib import Path

# %% =================== save matplotlib figures to pdf ================== %% #
def savefigs(filepath, figs='all'):
    """Save multiple matplotlib figures in a pdf.

    Args:
        filepath (string or pathlib.Path): filepath
        figs (list, optional): list with the figure numbers to save. Use 'all'
            to save all opened figures.
    """
    # check extension
    if Path(filepath).suffix == '.pdf':
        pass
    else:
        filepath = Path(filepath).with_suffix('.pdf')

    if figs == 'all':
        figs = [plt.figure(n) for n in plt.get_fignums()]

    if len(figs) > 1:
        pp = PdfPages(str(filepath))
        for fig in figs:
            fig.savefig(pp, format='pdf')
        pp.close()
    else:
        plt.savefig(str(filepath), format='pdf')