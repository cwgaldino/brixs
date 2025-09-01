#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for everyday use ---> vector manipulation"""

# %% ------------------------- Standard Imports -------------------------- %% #
import numpy as np

# %% ============================== functions ============================ %% #
def Rx(angle):
    """Returns the rotation matrix x given a angle in degrees"""
    return np.array([[1,                0,                     0              ],
                     [0, np.cos(np.radians(angle)), -np.sin(np.radians(angle))],
                     [0, np.sin(np.radians(angle)),  np.cos(np.radians(angle))]]
                    )

def Ry(angle):
    """Returns the rotation matrix y given a angle in degrees"""
    return np.array([[ np.cos(np.radians(angle)), 0, np.sin(np.radians(angle))],
                     [             0,             1,              0           ],
                     [-np.sin(np.radians(angle)), 0, np.cos(np.radians(angle))]]
                    )

def Rz(angle):
    """Returns the rotation matrix z given a angle in degrees"""
    return np.array([[np.cos(np.radians(angle)), -np.sin(np.radians(angle)), 0],
                     [np.sin(np.radians(angle)),  np.cos(np.radians(angle)), 0],
                     [               0,                      0,              1]]
                    )

def change2prime(x_prime=(1, 0, 0), y_prime=(0, 1, 0), z_prime=(0, 0, 1), vector=(1, 0, 0)):
    """Return vector in terms of a different coordinate system (prime system)

    Args:
        x_prime, y_prime, z_prime (vector, optional): new coordinate system 
            x', y', and z' in terms of old coordinates x, y, z.
        vector (vector, optional): Default is (1, 0, 0).

    Returns:
        vector'
    """
    vector_prime = (np.matmul(vector, x_prime), np.matmul(vector, y_prime), np.matmul(vector, z_prime))
    return np.array(chop_vector(vector_prime))

def rotate_vector(axis, angle, vector):
    """Returns a rotated vector

    Args:
        axis (str): reference rotation axis. Options: `x`, `y`, or `z`
        angle (number): angle in degrees
        vector (array): three item array

    Returns:
        rotated vector (array)
    """
    if 'x' in axis:
        return np.array(chop_vector(np.matmul(Rx(angle), vector)))
    elif 'y' in axis:
        return np.array(chop_vector(np.matmul(Ry(angle), vector)))
    elif 'z' in axis:
        return np.array(chop_vector(np.matmul(Rz(angle), vector)))

def rotate_system(axis, angle, x_initial=(1, 0, 0), y_initial=(0, 1, 0), z_initial=(0, 0, 1)):
    """Return a rotated coordinate system (relative to an old coordinate system)
    
    Args:
        axis (str): reference rotation axis. Options: `x`, `y`, or `z`
        angle (number): angle in degrees
        x_initial, y_initial, z_initial (array): three item array defining the old
            coordinate system, e.g., where x, y, and z points toward.

    Returns:
        rotated vector (array)
    """
    if 'x' in axis:
        return np.array(chop_vector(np.matmul(Rx(angle), x_initial))), np.array(chop_vector(np.matmul(Rx(angle), y_initial))), np.array(chop_vector(np.matmul(Rx(angle), z_initial)))
    elif 'y' in axis:
        return np.array(chop_vector(np.matmul(Ry(angle), x_initial))), np.array(chop_vector(np.matmul(Ry(angle), y_initial))), np.array(chop_vector(np.matmul(Ry(angle), z_initial)))
    elif 'z' in axis:
        return np.array(chop_vector(np.matmul(Rz(angle), x_initial))), np.array(chop_vector(np.matmul(Rz(angle), y_initial))), np.array(chop_vector(np.matmul(Rz(angle), z_initial)))

def chop_vector(vector, epsilon=10e-14):
    """change small values to zero"""
    return [x if np.abs(x) > epsilon else 0 for x in vector]

def is_perpendicular(vector1, vector2, epsilon=10e-14):
    """Return True if vectors are perpendicular"""
    if np.dot(vector1, vector2) > epsilon: 
        return False
    return True

def normalize_vector(vector):
    """Returns normalize vector to 1"""
    if np.linalg.norm(vector) != 1:
        return np.array(vector/np.linalg.norm(vector))
    else:
        return np.array(vector)

def angle_between_vectors(vector1, vector2):
    """Returns the angle between two vectors in degrees"""
    vector1 = normalize(vector1)
    vector2 = normalize(vector2)
    return np.degrees(np.arccos(np.dot(vector1, vector2)))


# %% =============================== drawing ============================= %% #
def draw_octahedron(ax, center=(0, 0, 0), x=(1, 0, 0), y=(0, 1, 0), z=(0, 0, 1), color='blue'):
    """Draw octahedron in a axes"""
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    # xy
    _ = ax.quiver(*np.concatenate((x, y-x)), color=color, arrow_length_ratio=0)
    _ = ax.quiver(*np.concatenate((x, -y-x)), color=color, arrow_length_ratio=0)
    _ = ax.quiver(*np.concatenate((-x, -y+x)), color=color, arrow_length_ratio=0)
    _ = ax.quiver(*np.concatenate((-x, +y+x)), color=color, arrow_length_ratio=0)

    # xz
    _ = ax.quiver(*np.concatenate((x, z-x)), color=color, arrow_length_ratio=0)
    _ = ax.quiver(*np.concatenate((x, -z-x)), color=color, arrow_length_ratio=0)
    _ = ax.quiver(*np.concatenate((-x, -z+x)), color=color, arrow_length_ratio=0)
    _ = ax.quiver(*np.concatenate((-x, +z+x)), color=color, arrow_length_ratio=0)

    # yz
    _ = ax.quiver(*np.concatenate((y, z-y)), color=color, arrow_length_ratio=0)
    _ = ax.quiver(*np.concatenate((y, -z-y)), color=color, arrow_length_ratio=0)
    _ = ax.quiver(*np.concatenate((-y, -z+y)), color=color, arrow_length_ratio=0)
    _ = ax.quiver(*np.concatenate((-y, +z+y)), color=color, arrow_length_ratio=0)

    return