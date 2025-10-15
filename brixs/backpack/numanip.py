#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for everyday use ---> int and floats"""

# %% ------------------------- Standard Imports --------------------------- %% #
from functools import reduce
import numpy as np
import decimal
import numbers

# %% ============================= functions ============================== %% #
def is_integer(n, allow_float=True, allow_str=True):
    """Returns True if number is integer. 
    
    Args:
        n (number): number
        allow_float (bool, optional): if True, floats can be considered int if 
            int(number) = number, i.e., number with no decimal places.
        allow_str (bool, optional): if True, it tries to convert str to number.
            then it is considered int if int((float(number)) = float(number), 
            i.e., string is number with no decimal places.

    returns
        bool
    """
    if isinstance(n, int):
        return True
    elif isinstance(n, np.int8):
        return True
    elif isinstance(n, np.int16):
        return True
    elif isinstance(n, np.int32):
        return True
    elif isinstance(n, np.int64):
        return True
    elif isinstance(n, np.uint8):
        return True
    elif isinstance(n, np.uint16):
        return True
    elif isinstance(n, np.uint32):
        return True
    elif isinstance(n, np.uint64):
        return True
    
    if allow_float:
        if isinstance(n, float):
            return n.is_integer()
    if allow_str:
        if isinstance(n, str):
            try: 
                return float(n).is_integer()
            except ValueError:
                return False

    return False
    
def is_number(n):
    """Returns True if variable is number."""
    # backpack developers note!!!!
    # if this function changes, it needs to be copied to these files: figmanip
    if isinstance(n, bool):
        return False
    elif isinstance(n, int):
        return True
    elif isinstance(n, float):
        return True
    elif isinstance(n, numbers.Number):
        return True
    elif isinstance(n, str):
        try: 
            is_number(float(n))
            return True
        except ValueError:
            return False
    else:
        False

def round_to_1(x):
    """return the most significant digit"""
    # backpack developers note!!!!
    # if this function changes, it needs to be copied to these files: figmanip
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
    # backpack developers note!!!!
    # if this function changes, it needs to be copied to these files: figmanip

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
