#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module for calculating crystal parameters."""

import numpy as np
from scipy.constants import h, speed_of_light, physical_constants

try:
    from pbcpy.base import DirectCell, ReciprocalCell
except:
    pass

def eV2angstrom(value):
    return h*speed_of_light/(physical_constants['electron volt-joule relationship'][0]*value) * 10**10

def momentum_transfer(energy, theta_i=None, theta_f=None, two_theta=None):

    # calculate wavelength
    wavelength = eV2angstrom(energy)

    # check angles
    if sum(x is None for x in [theta_i, theta_f, two_theta]) > 1:
        raise ValueError('At least two angles must be defined (theta_i, theta_f, two_theta)')
    if theta_i is None:
        theta_i = two_theta - theta_f
    elif theta_f is None:
        theta_f = two_theta - theta_i
    elif two_theta is None:
        two_theta = theta_i + theta_f
    else:
        if two_theta - (theta_i + theta_f) > two_theta*0.001 :
            raise ValueError('theta_i + theta_f must be equal to two_theta.')

    # degrees to radians
    theta_i   = np.radians(theta_i)
    theta_f   = np.radians(theta_f)
    two_theta = np.radians(two_theta)

    # calculate total momentum transfer
    q = 4*np.pi/wavelength * np.sin(two_theta/2)  # 1/A

    # projected momentum
    # q_parallel      = q*(np.cos(theta_f)-np.cos(theta_i))
    # q_perperdicular = q*(np.sin(theta_f)+np.sin(theta_i))

    q_parallel      = q*(2*np.sin(two_theta/2))**-1   *(np.cos(theta_f)-np.cos(theta_i))
    q_perperdicular = q*(2*np.sin(two_theta/2))**-1  *(np.sin(theta_f)+np.sin(theta_i))

    return q, q_parallel, q_perperdicular

def lattice(a, b, c, alpha=90, beta=90, gamma=90):
    """Returns the module of the reciprocal vectors.

    Not fully tested.
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


    cell1 = DirectCell(lattice=lattice, origin=[0,0,0])
    reciprocal_cell1 = cell1.get_reciprocal()

    return lattice, reciprocal_cell1.lattice

def brillouin_zone_size(reciprocal_lattice, hkl=None, h=None, k=None, l=None):
    """
    t1 = [rlatt[i]**2 for i in range(3)]
    t2 = [np.sqrt(sum(x)) for x in t1]
    t3 = [t2[i]*inplane_vector[i] for i in range(3)]
    dist = np.sqrt(sum([x**2 for x in t3]))
    print(dist)

    dist = 2*np.pi*np.sqrt(inplane_vector[0]**2/a**2 + inplane_vector[1]**2/b**2 + inplane_vector[2]**2/c**2)
    print(dist)
    """

    if hkl is not None:
        h = hkl[0]
        k = hkl[0]
        l = hkl[0]

    t1 = [rlatt[i]**2 for i in range(3)]
    t2 = [np.sqrt(sum(x)) for x in t1]
    t3 = [t2[i]*inplane_vector[i] for i in range(3)]
    dist = np.sqrt(sum([x**2 for x in t3]))

    return dist

def momentum2rlu(q, hkl, a, b, c, alpha=90, beta=90, gamma=90):
    _, rlatt = lattice(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
    size = brillouin_zone_size(rlatt, hkl=hkl)
    return q/size
