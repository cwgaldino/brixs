#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Matplotlib plotting support functions."""


# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# backpack
from .backpack.filemanip import filelist
from .backpack.arraymanip import all_equal
from .backpack.model_functions import voigt_fwhm, voigt_area_fwhm, dirac_delta
from .backpack.figmanip import n_digits

# BRIXS
import brixs as br


# %%

def _xylabels(title, ax, xlabel, ylabel, grid):
    """Core function for quickly setting x and y labels in plots."""
    if ax is not None:
        if title is not None:
            ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if grid:
            ax.grid(visible=True)
        else:
            ax.grid(visible=False)
    else:
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if grid:
            plt.grid(visible=True)

    plt.tight_layout()


def label_xas(title=None, ax=None, grid=True):
    """Quickly set x and y labels and title in XAS plots."""
    xlabel = 'Photon energy (eV)'
    ylabel = 'Intensity (arb. units)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)


def label_rixs(title=None, ax=None, grid=True):
    """Quickly set x and y labels and title in RIXS plots."""
    xlabel = 'Energy loss (eV)'
    ylabel = 'Intensity (arb. units)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)

    
def label_emap(title=None, ax=None, grid=True):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = 'Photon energy (eV)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)

def label_emission_map(title=None, ax=None, grid=True):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = 'Emission energy (eV)'
    ylabel = 'Photon energy (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)
    

def label_thmap(title=None, ax=None, grid=True):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = r'$\Theta$ ($^\circ$)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)


def label_qmap_rlu(title=None, ax=None, grid=True):
    """Quickly set x and y labels and title in momentum map plots."""
    xlabel = 'Momentum transfer (RLU)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)