#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Matplotlib plotting support functions."""


# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

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
        plt.title(title)
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

def label_emission(title=None, ax=None, grid=True):
    """Quickly set x and y labels and title in emission rixs plots."""
    xlabel = 'Emission energy (eV)'
    ylabel = 'Intensity (arb. units)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)

    
def label_energy_map(title=None, ax=None, grid=False):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = 'Photon energy (eV)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)

def label_emission_map(title=None, ax=None, grid=False):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = 'Emission energy (eV)'
    ylabel = 'Photon energy (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)
    

def label_th_map(title=None, ax=None, grid=False):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = r'$\Theta$ ($^\circ$)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)

# Experimental
def label_qmap_rlu(title=None, ax=None, grid=False):
    """Quickly set x and y labels and title in momentum map plots."""
    xlabel = 'Momentum transfer (RLU)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)