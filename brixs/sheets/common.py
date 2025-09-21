#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for spreadsheets

Usage:
    >>> from brixs.sheets.common import letter2num, str2num
"""

# %% -------------------------- standard imports -------------------------- %% #
import re

# %% -------------------------- support functions ------------------------- %% #
def letter2num(letter):
    """Returns position of letter in the alphabet starting from 0 
    
    Args:
        letter (str): string with a letter (A, B, ..., Z, AA, AB, ..., ZZ, AAA, ...)

    Example:
        >>> letter2num('A')   # -> 0
        >>> letter2num('AA')  # -> 26
    
    Returns:
        number
    """
    letter   = letter.lower()
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    n = 0
    for idx, s in enumerate(letter):
        n += alphabet.index(s) + (idx)*26
    return n

def str2num(string):
    """Returns col, row number based on string like 'A2', 'AJ20', 'A2:AJ20', ...

    Args:
        string (str): cell or range string

    Example:
        >>> str2num('A2')     # -> (0, 1)
        >>> str2num('A2:C4')  # -> (0, 1, 2, 3)

    Returns:
        col, row for `cell` strings ('A2', 'AJ20', ...)
        col_start, row_start, col_stop, row_stop for `range` strings ('A2:AJ20', ...) 
        col and row numbering starts from 0
    """
    # check input
    if isinstance(string, str) == False:
        raise ValueError(f'cell input must be a string, not `{type(string)}`')
        # msgbox(f'_cell2num error: input must be type str\n\ninput = {string}\ninput type={type(string)}')
        # return
    if ':' in string:
        start = str2num(string.split(':')[0])
        stop  = str2num(string.split(':')[1])
        return start[0], start[1], stop[0], stop[1]
    
    # cell2num
    try:
        temp = re.compile("([a-zA-Z]+)([0-9]+)")
        res  = temp.match(string).groups()
        return letter2num(res[0]), int(res[1])-1
    except AttributeError:
        raise ValueError(f'cannot covert input (`{string}`) to number')
        # msgbox(f'str2num error: cannot covert input to number\n\ninput = {string}')
        # return
# %%
