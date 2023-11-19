#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Matplotlib plotting support functions."""


# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from collections.abc import Iterable
from matplotlib.patches import Rectangle

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
def draw_mask(mask, ax=None,  **kwargs):
    """[Experimental] Draw rectangles.

        mask (list): list with rectangles coordinates (x_start, x_stop, y_start, y_stop).
        **kwargs: kwargs are passed to ``matplotlib.patches.Rectangle()`` 
            that plots the data.
        
    Returns:
        None
    """
    if ax is None:
        ax = plt.gca()

    # assert mask is the right format
    assert isinstance(mask, Iterable), 'mask must be iterable'
    if len(mask) == 4:
        if isinstance(mask[0], Iterable) == False:
            mask = [mask, ]

    # kwargs
    # if 'color' not in kwargs:
    #     kwargs['color'] = 'red'
            
    for m in mask:
        w = m[1] - m[0]
        h = m[3] - m[2]
        rect = Rectangle((m[0], m[2]), w, h, fill=False, **kwargs)
        ax.add_patch(rect)

# Experimental
def label_qmap_rlu(title=None, ax=None, grid=False):
    """Quickly set x and y labels and title in momentum map plots."""
    xlabel = 'Momentum transfer (RLU)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)