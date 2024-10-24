#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""common x and y labels for xas and rixs plots."""

# %% ------------------------- Standard Imports --------------------------- %% #
import matplotlib.pyplot as plt
import matplotlib as mpl

# %% ================================ core ================================ %% #
def _xylabels(title, ax, xlabel, ylabel, grid, **kwargs):
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
        plt.title(title, **kwargs)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if grid:
            plt.grid(visible=True)
        else:
            plt.grid(visible=False)
        # plt.tight_layout()
    return

# %% =============================== labels =============================== %% #
def xas(title=None, grid=True, ax=None, **kwargs):
    """Quickly set x and y labels and title in XAS plots."""
    xlabel = 'Photon energy (eV)'
    ylabel = 'Intensity (arb. units)'
    return _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid, **kwargs)
mpl.axes.Axes.labels_xas = lambda self, title=None, grid=True: xas(ax=self, title=title, grid=grid)

def rixs(title=None, grid=True, ax=None, **kwargs):
    """Quickly set x and y labels and title in RIXS plots."""
    xlabel = 'Energy loss (eV)'
    ylabel = 'Intensity (arb. units)'
    return _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid, **kwargs)
mpl.axes.Axes.labels_rixs = lambda self, title=None, grid=True: rixs(ax=self, title=title, grid=grid)

def time_trace(title=None, grid=True, ax=None, **kwargs):
    """Quickly set x and y labels and title in time-trace plots."""
    xlabel = 'Delay (ps)'
    ylabel = 'Intensity (arb. units)'
    return _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid, **kwargs)
mpl.axes.Axes.labels_time_trace = lambda self, title=None, grid=True: time_trace(ax=self, title=title, grid=grid)

def detector(title=None, grid=False, ax=None, **kwargs):
    """Quickly set x and y labels and title in photon events plots."""
    xlabel = 'x (pixels)'
    ylabel = 'y (pixels)'
    return _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid, **kwargs)
mpl.axes.Axes.detector = lambda self, title=None, grid=True: detector(ax=self, title=title, grid=grid)

def rixs_emission(title=None, grid=True, ax=None, **kwargs):
    """Quickly set x and y labels and title in emission rixs plots."""
    xlabel = 'Emission energy (eV)'
    ylabel = 'Intensity (arb. units)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid, **kwargs)
mpl.axes.Axes.labels_rixs_emission = lambda self, title=None, grid=True: rixs_emission(ax=self, title=title, grid=grid)

def energy_map(title=None, grid=False, ax=None, **kwargs):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = 'Photon energy (eV)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid, **kwargs)
mpl.axes.Axes.labels_energy_map = lambda self, title=None, grid=True: energy_map(ax=self, title=title, grid=grid)

def emission_map(title=None, grid=False, ax=None, **kwargs):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = 'Emission energy (eV)'
    ylabel = 'Photon energy (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid, **kwargs)
mpl.axes.Axes.labels_emission_map = lambda self, title=None, grid=True: emission_map(ax=self, title=title, grid=grid)
    
def th_map(title=None, grid=False, ax=None, **kwargs):
    """Quickly set x and y labels and title in energy map plots."""
    xlabel = r'$\Theta$ ($^\circ$)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid, **kwargs)
mpl.axes.Axes.labels_th_map = lambda self, title=None, grid=True: th_map(ax=self, title=title, grid=grid)

def qmap_rlu(title=None, grid=False, ax=None, **kwargs):
    """Quickly set x and y labels and title in momentum map plots."""
    xlabel = 'Momentum transfer (RLU)'
    ylabel = 'Energy loss (eV)'
    _xylabels(title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, grid=grid, **kwargs)
mpl.axes.Axes.labels_qmap_rlu = lambda self, title=None, grid=True: qmap_rlu(ax=self, title=title, grid=grid)
