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

def n_decimal_places(number, count_zero=False):
    """Return the number of decimal places of a number.

    Args:
        number (float or int): number.
        count_zero (bool, optional): if an integer is type float, it will came
            out with a zero after the decimal point, eg, `145.0` instead of `145`.
            count_zero == False will not count zeros after the decimal point.
            Default is False.

    Returns:
        number of decimal places in number.
    """

    n = abs(decimal.Decimal(str(number)).as_tuple().exponent)
    if count_zero == False and n == 1:
        if str(number)[-1] == '0':
            return 0
    return n

def n_digits(number, count_zero=False):
    """Return the number of digits of a number.

    Args:
        number (float or int): number.
        count_zero (bool, optional): if number is type float, it will came
            out with a zero after the decimal point, eg, `145.0` instead of `145`.
            count_zero == False will not count zeros after the decimal point.
            Default is False

    Returns:
        number of digits
    """
    a = len(str(int(number)))
    b = n_decimal_places(number, count_zero=count_zero)

    # fix zero (minus sign when number > -1 and < 0)
    if int(number) == 0 and number < 0:
        a += 1

    # if int return a
    if type(number) == int:
        return a
    
    # if float
    if count_zero == False and b == 1:
        if str(number)[-1] == '0':
            b = 0
        else:
            b += 1
    elif b != 0:
        b += 1

    return a + b
    
def factors(n):
    """Return a tuple with all the factors of a number."""
    return set(reduce(list.__add__,
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))
