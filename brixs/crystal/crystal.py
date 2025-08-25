#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module with function for calculating momentum transfer in single crystals.
 
It is assumed that the photon hits the crystal surface at a angle th and is 
scattered in a 2th angle as the drawing below

     \      /.
      \    /   .
       \  /     .
   th ( \/       .
   ┌──────────┐  . 2th
   ├  sample  ┤  .
   └──────────┘ .
           \   .
            \.
"""

# %% --------------------------- Standard Imports ------------------------- %% #
import numpy as np
from scipy.constants import h, speed_of_light, physical_constants

# %% --------------------------- special Imports -------------------------- %% #
from pbcpy.base import DirectCell
# %%

# %% =============================== crystal ============================== %% #
def ev2angstrom(value):
    """Converts value from energy (eV) to photon wavelength (angstrom)."""
    return h*speed_of_light/(physical_constants['electron volt-joule relationship'][0]*value) * 10**10

def calculate_q_transfer(energy, theta, two_theta):
    """Returns the momentum transfer q=q_f-q_i of a scattering process (in angstrom^-1).

     \      /.
      \    /   .
       \  /     .
   th ( \/      . 2th
    ---------  . 
          \   .
           \.
    
    Args:
        energy (number): energy (in eV) of the incident beam.
        theta (number): angle (in degrees) between the incident beam and the sample surface.
        two_theta (number): angle (in degrees) between the incident and scattered beam.

    Returns:
        q, q_parallel, q_perpendicular
    """
    # calculate wavelength
    wavelength = ev2angstrom(energy)

    # check angles
    # if sum(x is None for x in [theta_i, theta_f, two_theta]) > 1:
    #     raise ValueError('At least two angles must be defined (theta_i, theta_f, two_theta)')
    # if theta_i is None:
    #     theta_i = two_theta - theta_f
    # elif theta_f is None:
    #     theta_f = two_theta - theta_i
    # elif two_theta is None:
    #     two_theta = theta_i + theta_f
    # else:
    #     if two_theta - (theta_i + theta_f) > two_theta*0.001 :
    #         raise ValueError('theta_i + theta_f must be equal to two_theta.')

    # degrees to radians
    theta_i   = np.radians(theta)
    theta_f   = np.radians(two_theta - theta)
    two_theta = np.radians(two_theta)

    # calculate total momentum transfer
    q = 4*np.pi/wavelength * np.sin(two_theta/2)  # 1/A

    # projected momentum
    # q_parallel      = q*(np.cos(theta_f)-np.cos(theta_i))
    # q_perpendicular = q*(np.sin(theta_f)+np.sin(theta_i))

    q_parallel      = q*(2*np.sin(two_theta/2))**-1   *(np.cos(theta_f)-np.cos(theta_i))
    q_perpendicular = q*(2*np.sin(two_theta/2))**-1  *(np.sin(theta_f)+np.sin(theta_i))

    return q, q_parallel, q_perpendicular

def _lattice(a, b, c, alpha=90, beta=90, gamma=90):
    """Returns the unit cell and the reciprocal unit cell matrix in real space.

    Args:
        a, b, c (number): lattice parameters in angstrom.
        alpha, beta, gamma (number): unit cell angles.

    Returns:
        3x3 unit cell matrix, 3x3 reciprocal unit cell matrix
    """
    if alpha == 90 and beta == 90 and gamma == 90:
        a = [a, 0, 0]
        b = [0, b, 0]
        c = [0, 0, c]
    else:
        alpha = np.radians(alpha)
        beta = np.radians(beta)
        gamma = np.radians(gamma)
        cos_delta = (np.cos(alpha) - np.cos(gamma)*np.cos(beta)) / (np.sin(gamma)*np.sin(beta))
        sin_delta = np.sin(np.arccos(cos_delta))
        a = [a, 0, 0]
        b = [b*np.cos(gamma), b*np.sin(gamma), 0]
        c = [c*np.cos(beta), c*np.sin(beta)*cos_delta, c*np.sin(beta)*sin_delta]

    lattice = np.array([a, b, c])

    try:
        cell1 = DirectCell(lattice=lattice, origin=[0,0,0])
    except NameError:
        raise NameError('pbcpy does not seem to be installed.\npbcpy.base.DirectCell() cannot be found.\nPlease, install pbcpy (pip install pbcpy).')
    reciprocal_cell1 = cell1.get_reciprocal()

    return lattice, reciprocal_cell1.lattice

def _brillouin_zone_size(reciprocal_lattice, h=None, k=None, l=None):
    """Returns the size of the brillouin zone.

    Args:
        reciprocal_lattice (matrix): 3x3 reciprocal unit cell matrix. See :py:func:`lattice`.
        h, k, l (number): miller indices defining a direction.

    Returns:
        number
    """
    inplane_vector = (h, k, l)

    t1 = [reciprocal_lattice[i]**2 for i in range(3)]
    t2 = [np.sqrt(sum(x)) for x in t1]
    t3 = [t2[i]*inplane_vector[i] for i in range(3)]
    dist = np.sqrt(sum([x**2 for x in t3]))
    # dist = 2*np.pi*np.sqrt(inplane_vector[0]**2/a**2 + inplane_vector[1]**2/b**2 + inplane_vector[2]**2/c**2)

    return dist

def momentum2rlu(q, hkl, a, b, c, alpha=90, beta=90, gamma=90):
    """Converts from momentum to Relative Lattice Units (RLU).

    Args:
        q (number): momentum (in angstrom^-1)
        hkl (tuple): list or tuple with 3 elements corresponding to miller
            indices defining a direction.
        a, b, c (number): lattice parameters in angstrom.
        alpha, beta, gamma (number): unit cell angles.

    Returns:
        number
    """
    _, rlatt = _lattice(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
    size = _brillouin_zone_size(rlatt, h=hkl[0], k=hkl[1], l=hkl[2])
    return q/size
