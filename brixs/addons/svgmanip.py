#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for everyday use ---> svg files"""

# %% ------------------------- Standard Imports -------------------------- %% #
from collections import OrderedDict
from pathlib import Path
import copy

# %% ------------------------- Special Imports --------------------------- %% #
try:
    from bs4 import BeautifulSoup
except:
    pass

# %% ============================== svg ================================== %% #
def ungroup_svg(filepath, outfilepath=None):
    """Ungroup all objects in svg files created by matplotlib.

    Text with latex equations are kept grouped.

    Ticks are cloned. The parent ticks will be at the left upper corner of the figure.

    Args:
        filepath (string or pathlib.Path): path to svg file.
        outfilepath (string or pathlib.Path, optional): path for svg output file.
            If None, it will use the input filepath and the file will be overwrited.

    See Also:
        :py:func:`soft_ungroup`

    Notes:
        Use :py:func:`soft_ungroup` to better preserve the structure
        of the svg file in case :py:func:`ungroup` returns strange results.
    """
    # get soup
    filepath = Path(filepath)
    f = filepath.open()
    soup = BeautifulSoup(f, features = 'xml')
    f.close()

    # file prefix
    script_tags = soup.find_all('svg')
    prefix = script_tags[0].prettify().split('<g')[0]

    # end of file
    endOfFile = script_tags[0].prettify().split('</g>')[-1]

    # find objs
    script_tags = soup.find_all(['use', 'path', 'def', 'text'])
    objs = OrderedDict()  # dict.key() is the layer label, dict.value() is the layer itself
    for i in range(len(script_tags)):

        if script_tags[i].name == 'text':  # keep text with latex equations grouped
            if 'transform' in script_tags[i].parent.attrs:
                objs[i] = str(script_tags[i].parent.prettify())
            else:
                objs[i] = str(script_tags[i].prettify())
        else:
            # include obj
            objs[i] = str(script_tags[i].prettify())

    # prepare output
    output = copy.copy(prefix)
    for id in objs:
        output += objs[id]
    output += endOfFile

    # save
    if outfilepath is None:
        outfilepath = filepath
    else:
        outfilepath = Path(outfilepath)
    f = outfilepath.open('w')
    f.write(output)
    f.close()

def soft_ungroup_svg(filepath, outfilepath=None):
    """Ungroup objects in svg file (each object stands in a group by itself).

    Args:
        filepath (string or pathlib.Path): path to svg file.
        outfilepath (string or pathlib.Path, optional): path for svg output file.
            If None, it will use the input filepath and the file will be overwritten.

    See Also:
        :py:func:`ungroup`

    Notes:
        This differs from :py:func:`ungroup` as it better preserves the structure
        of the svg file, but Objects are still kept in groups.
    """
    # get soup
    filepath = Path(filepath)
    f = filepath.open()
    soup = BeautifulSoup(f, features = 'xml')
    f.close()

    # file prefix
    script_tags = soup.find_all('svg')
    prefix = script_tags[0].prettify().split('<g')[0]

    # end of file
    endOfFile = script_tags[0].prettify().split('</g>')[-1]

    # find objs
    script_tags = soup.find_all('g')
    objs = OrderedDict()  # dict.key() is the layer label, dict.value() is the layer itself
    for i in range(len(script_tags)):
        if 'id' in script_tags[i].attrs:

            if script_tags[i].attrs['id'].startswith(('figure', 'axes', 'xtick', 'matplotlib.axis', 'xtick', 'ytick')):
                pass
            else:
                # print(script_tags[i].attrs['id'])
                del_list = []
                for j in range(len(script_tags[i].contents)):
                    try:
                        if 'id' in script_tags[i].contents[j].attrs:
                            del_list.append(j)
                    except AttributeError:
                        pass
                del_list = [x-n for n,x in enumerate(del_list)]
                for j in del_list:
                    del script_tags[i].contents[j]

                # include obj
                objs[script_tags[i].attrs['id']] = str(script_tags[i])#.prettify()

    # prepare output
    output = copy.copy(prefix)
    for id in objs:
        output += objs[id]
    output += endOfFile

    # save
    if outfilepath is None:
        outfilepath = filepath
    else:
        outfilepath = Path(outfilepath)
    f = outfilepath.open('w')
    f.write(output)
    f.close()
# %%