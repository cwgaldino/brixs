#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Everyday use functions for file handling."""

import os
import sys
from pathlib import Path
import numpy as np
import datetime
from copy import deepcopy
import collections
from .intermanip import query_yes_no
import json
from detect_delimiter import detect
# %%
import warnings
import re
# %%

def rename_files(filelist, pattern, newPattern, ask=True):
    """Change the filename pattern of files.

    Args:
        filelist (list): list of filepaths (string or pathlib.Path object).
        pattern (str): string that represents the filename pattern. Use ``{}`` to collect values within the file name.
        newPattern (str): Each chunk of information assigned with ``{}`` in ``pattern``
            is assigned a number based on its position. The first chunk is
            0, then 1, and so on. Use ``{position}`` to define the new pattern,

            Example:
                If filename is 'data_5K-8.dat', ``pattern = {}_{}K-{}.dat`` and
                ``newPattern={2}_{1}K.dat``. The new filename will be: '8_5K.dat'.

            Tip:
                Use ``\`` as an escape key.

            Tip:
                Use ``{n}`` in ``newPattern`` to inlcude the filenane index number (index within filelist).

        ask (bool): If true, it shows all the new filenames and asks for permission to rename.
    """
    permission = True
    n_infos = pattern.count('{}')
    pattern = pattern.replace('{}', '(.+?)')

    newPattern = newPattern.replace('{n}', '[n]')

    a = re.findall('{.+?}', newPattern)
    a = [int(item.replace('{', '').replace('}', '')) for item in a]
    if n_infos < max(a)+1:
        raise AttributeError("newPattern has some {n} where n is bigger than the number of marked infos {} in pattern.")

    if ask:
        permission = False
        print('\n' + '='*20)
        print('The following files will be renamed:\n')

        for n, filePath in enumerate(filelist):
            filePath = Path(filePath)

            # new name
            a = re.search(pattern, filePath.name)
            a = [a.group(i) for i in range(1, n_infos+1)]

            newPattern_mod_1 = newPattern.replace('[n]', str(n))
            newPattern_mod_2 = newPattern_mod_1.replace("{", "{a[").replace("}", "]}")
            nameNew = eval("f'"+ newPattern_mod_2 + "'")

            # print('\n' + '='*20)
            # print('The following files will be renamed:\n')
            print('OLD NAME = ' + filePath.name)
            print('NEW NAME = ' + nameNew)
            print('--')

        permission = query_yes_no('Change names?', default="yes")


    if permission:
        for n, filePath in enumerate(filelist):
            filePath = Path(filePath)

            # new name
            a = re.search(pattern, filePath.name)
            a = [a.group(i) for i in range(1, n_infos+1)]

            newPattern_mod_1 = newPattern.replace('[n]', str(n))
            newPattern_mod_2 = newPattern_mod_1.replace("{", "{a[").replace("}", "]}")
            nameNew = eval("f'"+ newPattern_mod_2 + "'")

            filePath.rename(filePath.parent / nameNew)

        print('Files renamed!')
    else:
        warnings.warn('Files NOT renamed.')


def rmdir(dirpath):
    """Remove a directory and everyting in it.

    Args:
        dirpath (string or pathlib.Path): directory path.
    """
    dirpath = Path(dirpath)
    for item in dirpath.iterdir():
        if item.is_dir():
            rmdir(item)
        else:
            item.unlink()
    dirpath.rmdir()


def filelist(dirpath='.', string='*'):
    """Returns a list with all the files containg `string` in its name.

    Note:
        List is sorted by the filename.

    Args:
        dirpath (str or pathlib.Path, optional): list with full file directory
        paths.
        string (str, optional): string to look for in file names.

    Return:
        list

    See Also:
        :py:func:`parsed_filelist`
    """
    dirpath = Path(dirpath)

    if '*' not in string:
        string = '*' + string + '*'

    temp = list(dirpath.glob(string))

    temp2 = [filepath.name for filepath in temp]

    return [x for _,x in sorted(zip(temp2,temp))]


def parsed_filelist(dirpath='.', string='*', ref=0, type='int'):
    """Returns a filelist organized in a dictionary.

    I searches for numbers (float and int) within the filenames (or foldernames)
    and uses them as keys for the dictionary. Filenames (or foldernames) must
    have the same pattern.

    Args:
        dirpath (str or pathlib.Path, optional): directory path.
        string (str, optional): string to filter filenames.
        ref (int, optional): index of the reference number to be used as key.
        type (string, optional): if 'int', dict keys are transormed in integers.
            If 'float', dict keys are transformed into float.

    Returns:
        Dictionary. Dict keys are some number found in filename.

    See Also:
        :py:func:`filelist`
    """
    dirpath = Path(dirpath)

    file_list = filelist(dirpath=dirpath, string=string)

    temp = dict()

    p = '[\d]+[.,\d]+|[\d]*[.][\d]+|[\d]+'

    for filepath in file_list:
        if re.search(p, filepath.name) is not None:
            n = []
            for catch in re.finditer(p, filepath.name):
                n.append(catch[0])
            if type=='int':
                temp[int(float((n[ref])))] = filepath
            else:
                temp[float(n[ref])] = filepath


    # ordering
    a = list(temp.keys())
    a.sort()

    parsed_folder = collections.OrderedDict()
    for key in a:
        parsed_folder[key] = temp[key]

    return parsed_folder


def save_text(string, filepath='./Untitled.txt', checkOverwrite=False):
    """Save text to txt file.

    Args:
        string (str): string to be saved.
        filepath (str or pathlib.Path, optional): path to save file. If no path
            is given, current working directory is used.
        checkOverwrite (bool, optional): if True, it will check if file exists
            and ask if user want to overwrite file.

    See Also:
        :py:func:`load_text`
    """
    filepath = Path(filepath)

    if checkOverwrite:
        if filepath.exists() == True:
            if filepath.is_file() == True:
                if query_yes_no('File already exists!! Do you wish to ovewrite it?', 'yes') == True:
                    pass
                else:
                    warnings.warn('File not saved.')
                    return
            else:
                warnings.warn('filepath is pointing to a folder. Saving file as Untitled.txt')
                filepath = filepath/'Untitled.txt'

    f = open(str(filepath), 'w')
    f.write(string)
    f.close()


def load_text(filepath):
    """Load text from txt file.

    Args:
        filepath (str or pathlib.Path): filepath to load.

    Returns:
        string.

    See Also:
        :py:func:`save_text`
    """
    f = Path(filepath).open()
    text = f.read()
    f.close()
    return text


def save_obj(obj, filepath='./Untitled.txt', checkOverwrite=False, prettyPrint=True):
    """Save object (array, dictionary, list, etc...) to a txt file.

    Args:
        obj (object): object to be saved.
        filepath (str or pathlib.Path, optional): path to save file.
        checkOverwrite (bool, optional): if True, it will check if file exists
            and ask if user want to overwrite file.

    See Also:
        :py:func:`load_obj`
    """
    filepath = Path(filepath)

    if checkOverwrite:
        if filepath.exists() == True:
            if filepath.is_file() == True:
                if query_yes_no('File already exists!! Do you wish to ovewrite it?', 'yes') == True:
                    pass
                else:
                    warnings.warn('File not saved because user did not allow overwriting.')
                    return
            else:
                warnings.warn('filepath is pointing to a folder. Saving file as Untitled.txt')
                filepath = filepath/'Untitled.txt'

    with open(str(filepath), 'w') as file:
        if prettyPrint:
            file.write(json.dumps(obj, indent=4, sort_keys=False))
        else:
            file.write(json.dumps(obj))


def _to_int(obj):
    """Change keys of a dictionary from string to int when possible."""
    for key in list(obj.keys()):
        try:
            if float(key).is_integer():
                new_key = int(float(key))
        except:
            new_key = key
        if new_key != key:
            obj[new_key] = obj[key]
            del obj[key]
    return obj


def load_obj(filepath, dict_keys_to_int=False):
    """Load object (array, dictionary, list, etc...) from a txt file.

    Args:
        filepath (str or pathlib.Path): file path to load.
        dict_keys_to_int (bool, optional): If True, it will change ALL
            numeric dict keys (even for key in nested dictionarys to int, e.g.,
            dictObject["0.0"] will turn into dictObject[0].

    Returns:
        object.

    See Also:
        :py:func:`save_obj`
    """
    filepath = Path(filepath)

    with open(str(filepath), 'r') as file:
        if dict_keys_to_int:
            obj = json.load(file, object_hook=_to_int)
        else:
            obj = json.load(file)
    return obj


def load_Comments(filepath, comment_flag='#', stop_flag='#'):
    """Return comments from text files.

    Comments must be indicated at the begining of the line by the comment flag.

    Args:
        filepath (str or pathlib.Path): fullpath to file
        comment_flag (str, optional): string that indicate line is a comment.
        stop_flag (str, optional): string that indicates line to stop looking for
            comments, e.g. `#f`, or `#L`. Use `None` to read all lines in file (useful for reading
            file with comments not only at the beginning). If `stop_flag` is
            equal to `comment_flag` it will read from the first line with
            `comment_flag` and keep reading until `comment_flag` does not apper
            anymore (useful to read comments at the beginning of a file).

    Returns:
        list with comments or False if no comment is found.
    """
    comments = []
    filepath = str(Path(filepath))
    l = len(comment_flag)

    if stop_flag is None:
        with open(filepath) as file:
            for line in file:
                if line[0:l] == comment_flag:
                    comments.append(line[:])
    elif stop_flag == comment_flag:
        with open(filepath) as file:
            comment_started = 0
            for line in file:
                if line[0:l] == comment_flag and comment_started == 0:
                    comments.append(line[:])
                    comment_started = 1
                elif line[0:l] == comment_flag and comment_started == 1:
                    comments.append(line[:])
                elif line[0:l] != comment_flag and comment_started == 1:
                    break
    else:
        with open(fullpath) as file:
            for line in file:
                if line[0:len(stop_flag)] != stop_flag:
                    if line[0:len(comment_flag)] == comment_flag:
                        comments.append(line[:])
                else:
                    comments.append(line[:])
                    break

    if comments == []:
        return False
    else:
        return comments[:]


def save_data(obj, filepath='./untitled.txt', add_labels=True, data_format='% .10e', header='', footer='', delimiter=', ', comment_flag='# ', newline='\n', checkOverwrite=False):
    r"""Save an array or a dictionary in a txt file.

    Args:
        obj (dict, list, or numpy.array): data to be saved to a file. If obj is
        a dictonary, use ``*`` in front of a key to do not save it to the file.
        filepath (str or pathlib.Path, optional): path to save file.
        add_labels (bool, optional): When obj is a dictonary, ``add_labels=True``
            makes the dict keys to be added to the header as label for each data column.
        data_format (string, or list, optional): If obj is a list, fmt can also
            be a list where each fmt element is associated with a column. If
            obj is a dict, fmt can also be a dict with same keys of obj. Then,
            each fmt value is associated with the corresponding column.

            See `np.savetxt <https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html?highlight=savetxt#numpy.savetxt>`_ documentation::

                fmt = (%[flag]width[.precision]specifier)

            * flag can be: '-' for left justify, '+', whch forces to precede

            * result with + or -, or '0' to Left pad the number with zeros
              instead of space (see width).

            * width is the minimum number of characters to be printed.

            * precision is tipically the number of significant digits

            * specifier is the type of notation. Tipically, either 'e' for
              scientific notation of 'f' for decimal floating point.

            * a common fmt strings is: '%.3f' for 3 decimal places.

        header (str, oprional): string that will be written at the beggining of
            the file (comment flag is added automatically).
        footer (str, oprional): string that will be written at the end of the
            file (comment flag is added automatically).
        delimiter (str, optional): The string used to separate values.
        comment_flag (str, optional): string that flag comments.
        newline (str, optional): string to indicate new lines.
        checkOverwrite (bool, optional): if True, it will check if file exists
            and ask if user want to overwrite file.

    See Also:
        :py:func:`load_data`
    """
    filepath = Path(filepath)

    if checkOverwrite:
        if filepath.exists() == True:
            if filepath.is_file() == True:
                if query_yes_no('File already exists!! Do you wish to ovewrite it?', 'yes') == True:
                    pass
                else:
                    warnings.warn('File not saved.')
                    return
            else:
                warnings.warn('filepath is pointing to a folder. Saving file as Untitled.txt')
                filepath = filepath/'Untitled.txt'

    if type(obj) == dict:
        # remove keys that start with star (*)
        obj2 = {key: obj[key] for key in obj if str(key).startswith('*') is False}

        # col labels
        if labels:
            if not header == '' and not header.endswith('\n'):
                header += '\n'
            for key in obj2:
                header += str(key) + f'{delimiter}'
            header = header[:-(len(delimiter))]

        # dict to array
        obj = []
        for key in obj2:
            obj.append(obj2[key])
        obj = np.array(obj).transpose()

    np.savetxt(filepath, obj, fmt=data_format, delimiter=delimiter, newline=newline, header=header, footer=footer, comments=comment_flag)


def load_data(filepath, delimiter=None, comment_flag='#', labels=None, force_array=False):
    """Load data from text file. Data is formated in a dictionary or array.

    The dictionary keys are set as the label of the corresponding data columns, where
    the last comment line before data starts is assumed to have the labels of each data column.
    If column labels cannot be found, data is imported as an array. The file expects
    comments at the begining of the file (file must not have comments elsewhere).

    Args:
        filepath (str or pathlib.Path): path to file
        delimiter (str, optional): The string used to separate data values. If whitespaces are used,
            consecutive whitespaces act as delimiter. Use ``\\t`` for tab. If None,
            the script will read the file and try to guess the appropriate delimiter.
            If the file has a header, and if the header has a line with the column labels, it
            tries to guess the delimiter of this line. If it cannot, it tries
            to use the same delimiter for the data and for the leader in the header.
        comment_flag (str, optional): string indicating comments.
        labels (list, optional): It forces data to be loaded as a dictionary where
            each label is associated with a data column. Its lenght must have the same as the number of
            columns. To avoid importing a column, put an asterisk (*) in front of the corresponding label.
        force_array (bool, optional): If ``force_array=True``, data it will be returned in a array.

    Returns:
        Dictionary or array.

    See Also:
        :py:func:`save_data`.
    """
    filepath = Path(filepath)

    # guess delimiter
    if delimiter == ' ':
        delimiter = None
    elif delimiter is None:
        header = load_Comments(filepath, comment_flag=comment_flag, stop_flag=comment_flag)
        try: line_number = len(header)
        except TypeError: line_number = 0
        f = open(filepath)
        for i, line in enumerate(f):
            if i == line_number:
                delimiter = detect(line)
            elif i > line_number:
                break
        f.close()
        if delimiter is None:
            warnings.warn('Could not figure out the delimiter. Trying space.')

    # get data
    data = np.genfromtxt(str(filepath), delimiter=delimiter, comments=comment_flag)

    # get labels
    header = load_Comments(filepath, comment_flag=comment_flag, stop_flag=comment_flag)
    if labels is None:
        if header is False:
            warnings.warn('Cannot find header. Importing data as an array.')
            return data
        elif force_array:
            return data
        else:
            header_line = header[-1].replace(comment_flag, '').strip()
            header_line = header_line.replace('\n', '')
            header_delimiter = detect(header_line)
            labels = header_line.split(header_delimiter)
            if len(labels) != data.shape[1]:
                labels = header_line.split(delimiter)
                if len(labels) != data.shape[1]:
                    warnings.warn('Cannot find column labels. Importing data as an array.')
                    return data
            # remove empty itens and trailing spaces
            labels = [item.strip() for item in labels if item != '']

    # create dict
    datadict = {labels[i]: data[:, i] for i in range(len(labels)) if not labels[i].startswith('*')}

    # if a column returns only nan, this column is read as string
    for key in datadict:
        if False in [x for x in datadict[key] != datadict[key]]:  # test if the whole col is nan
            pass
        else:
            # key conversion to col number
            col_string = labels.index(key)
            data = np.genfromtxt(filepath, delimiter=None, comments='#', usecols=[-1], dtype='S8')
            datadict[key] = [x.decode('utf-8') for x in data]

    return datadict
