#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""common x and y labels for xas and rixs plots."""

# %% ------------------------- Standard Imports --------------------------- %% #
import matplotlib.pyplot as plt
import matplotlib as mpl

# %% ================================ core ================================ %% #
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
    return

# %% =============================== labels =============================== %% #
def xas(ax=None, title=None, grid=True):
    """Quickly set x and y labels and title in XAS plots."""
    xlabel = 'Photon energy (eV)'
    ylabel = 'Intensity (arb. units)'
    return _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)
mpl.axes.Axes.labels_xas = lambda self, title=None, grid=True: xas(self, title, grid)

def rixs(self, title=None, grid=True):
    """Quickly set x and y labels and title in RIXS plots."""
    xlabel = 'Energy loss (eV)'
    ylabel = 'Intensity (arb. units)'
    return _xylabels(title=title, ax=self, xlabel=xlabel, ylabel=ylabel, grid=grid)
mpl.axes.Axes.labels_rixs = lambda self, title=None, grid=True: rixs(self, title, grid)

def rixs_emission(title=None, ax=None, grid=True):
    """Quickly set x and y labels and title in emission rixs plots."""
    xlabel = 'Emission energy (eV)'
    ylabel = 'Intensity (arb. units)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)
mpl.axes.Axes.labels_rixs_emission = lambda self, title=None, grid=True: rixs_emission(self, title, grid)

def energy_map(title=None, ax=None, grid=False):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = 'Photon energy (eV)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)
mpl.axes.Axes.labels_energy_map = lambda self, title=None, grid=True: energy_map(self, title, grid)

def emission_map(title=None, ax=None, grid=False):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = 'Emission energy (eV)'
    ylabel = 'Photon energy (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)
mpl.axes.Axes.labels_emission_map = lambda self, title=None, grid=True: emission_map(self, title, grid)
    
def th_map(title=None, ax=None, grid=False):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = r'$\Theta$ ($^\circ$)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)
mpl.axes.Axes.labels_th_map = lambda self, title=None, grid=True: th_map(self, title, grid)

def qmap_rlu(title=None, ax=None, grid=False):
    """Quickly set x and y labels and title in momentum map plots."""
    xlabel = 'Momentum transfer (RLU)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid)
mpl.axes.Axes.labels_qmap_rlu = lambda self, title=None, grid=True: qmap_rlu(self, title, grid)
