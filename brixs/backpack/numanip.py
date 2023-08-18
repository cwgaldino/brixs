#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Everyday use functions for number manipulation."""

import numpy as np
import copy
import decimal
from functools import reduce
import numbers

def is_integer(n):
    """Returns True if number is integer."""
    if isinstance(n, int):
        return True
    elif isinstance(n, np.int32):
        return True
    elif isinstance(n, float):
        return n.is_integer()
    else:
        return False
    
def is_number(n):
    """Returns True if variable is number."""
    if isinstance(n, int):
        return True
    elif isinstance(n, float):
        return True
    elif isinstance(n, numbers.Number):
        return True
    else:
        False

def round_to_1(x):
    """return the most significant digit"""
    return round(x, -int(np.floor(np.log10(abs(x)))))

def n_decimal_places(number):
    """Return the number of decimal places of a number.

    Args:
        number (float or int): number.

    Returns:
        number of decimal places in number.
    """

    return abs(decimal.Decimal(str(number)).as_tuple().exponent)

def n_digits(number):
    """Return the number of digits of a number.

    Args:
        number (float or int): number.

    Returns:
        a tuple with number of digits and number of decimal places.
    """
    if n_decimal_places(number) != 0:
        return (len(str(int(np.around(number))) ), n_decimal_places(number))
    else:
        return (len(str(int(np.around(number))) ), 0)
    
def factors(n):
    """Return a tuple with all the factors of a number."""
    return set(reduce(list.__add__,
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))
