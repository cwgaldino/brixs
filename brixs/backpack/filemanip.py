#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for everyday use ---> files and folders"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import collections
import numpy as np
import datetime
import warnings
import fnmatch
import shutil
import json
import os
import re

# %% ----------------- supporting functions from query -------------------- %% #
# backpack developers note --> if these function change, it needs to be copied to query.py
import platform
def _operating_system():
    """Return string with name of operating system (windows, linux, or mac)."""
    system = platform.system().lower()
    is_windows = system == 'windows'
    is_linux = system == 'linux'
    is_mac = system == 'darwin'
    if is_windows:
        return 'windows'
    elif is_linux:
        return 'linux'
    elif is_mac:
        return 'mac'
    else:
        raise ValueError('OS not recognized')

is_windows = _operating_system() == 'windows'
is_linux   = _operating_system() == 'linux'
is_mac     = _operating_system() == 'mac'

import sys
def _query(question, default="yes"):
    """Ask a yes/no question and return answer.

    Note:
        It accepts many variations of yes and no as answer, like, "y", "YES", "N", ...

    Args:
        question (str): string that is presented to the user.
        default ('yes', 'no' or None): default answer if the user just hits
            <Enter>. If None, an answer is required of the user.

    Returns:
        True for "yes" or False for "no".
    """    
    valid = {"yes": True, "y": True, "ye": True, "Y": True, "YES": True, "YE": True,
             "no": False, "n": False, "No":True, "NO":True, "N":True}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt + '\n')
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "('y' or 'n').\n")

# %%

# %% ================================ files =============================== %% #
def rename_files(filelist, pattern, new_pattern, ask=True):
    r"""Change the filename pattern of files.

    Args:
        filelist (list): list of filepaths (string or pathlib.Path object).
        pattern (str): string that represents the filename pattern. Use ``{}`` to collect values within the file name.
        new_pattern (str): Each chunk of information assigned with ``{}`` in ``pattern``
            is assigned a number based on its position. The first chunk is
            0, then 1, and so on. Use ``{position}`` to define the new pattern,

            Example:
                If filename is 'data_5K-8.dat', ``pattern = {}_{}K-{}.dat`` and
                ``new_pattern={2}_{1}K.dat``. The new filename will be: '8_5K.dat'.

            Note:
                Use ``\`` as an escape key.

            Note:
                Use ``{n}`` in ``new_pattern`` to inlcude the filenane index number (index regarding filelist).

        ask (bool): If true, it shows all the new filenames and asks for permission to rename.
    """
    permission = True
    n_infos = pattern.count('{}')
    pattern = pattern.replace('{}', '(.+?)')

    new_pattern = new_pattern.replace('{n}', '[n]')

    a = re.findall('{.+?}', new_pattern)
    a = [int(item.replace('{', '').replace('}', '')) for item in a]
    if n_infos < max(a)+1:
        raise AttributeError("new_pattern has some {n} where n is bigger than the number of marked infos {} in pattern.")

    if ask:
        permission = False
        print('\n' + '='*20)
        print('The following files will be renamed:\n')

        for n, filePath in enumerate(filelist):
            filePath = Path(filePath)

            # new name
            a = re.search(pattern, filePath.name)
            a = [a.group(i) for i in range(1, n_infos+1)]

            newPattern_mod_1 = new_pattern.replace('[n]', str(n))
            newPattern_mod_2 = newPattern_mod_1.replace("{", "{a[").replace("}", "]}")
            nameNew = eval("f'"+ newPattern_mod_2 + "'")

            # print('\n' + '='*20)
            # print('The following files will be renamed:\n')
            print('OLD NAME = ' + filePath.name)
            print('NEW NAME = ' + nameNew)
            print('--')

        permission = _query('Change names?', default="yes")


    if permission:
        for n, filePath in enumerate(filelist):
            filePath = Path(filePath)

            # new name
            a = re.search(pattern, filePath.name)
            a = [a.group(i) for i in range(1, n_infos+1)]

            newPattern_mod_1 = new_pattern.replace('[n]', str(n))
            newPattern_mod_2 = newPattern_mod_1.replace("{", "{a[").replace("}", "]}")
            nameNew = eval("f'"+ newPattern_mod_2 + "'")

            filePath.rename(filePath.parent / nameNew)

        print('Files renamed!')
    else:
        warnings.warn('Files NOT renamed.')

def copy_file(filepath, folderpath):
    """Copy file at filepath to folderpath

    Args:
        filepath (str or path): file to be copied
        folderpath (str or path): destination folderpath

    Returns:
        None
    """
    folderpath = Path(folderpath)
    filepath = Path(filepath)
    
    assert folderpath.is_dir() == True,  'folderpath must point to a folder'
    assert filepath.is_dir()   == False, 'filepath must point to a file'
    assert filepath.exists()   == True,  'file cannot be found'
    shutil.copy(filepath, folderpath)
    return

def copy_files(target, destination, verbose=True): 
    """copy all files from one directory to another.

    Args:
        target (str or path): target folderpath
        destination (str or path): destination folderpath
        
    Returns:
        None
    """   
    assert target.is_dir()      == True,  'target must point to a folder'
    assert destination.is_dir() == True,  'destination must point to a folder'
    assert target.exists()      == True,  'target cannot be found'
    assert destination.exists() == True,  'destination cannot be found'

    for filepath in filelist(dirpath=target):
        copy_file(filepath, destination)
        if verbose: print(f'copied: {filepath}')
    if verbose: print('done')
    return

# %% ================================ directories ========================= %% #
def mkdir(dirpath):
    """create new directory (it also creates parent folders if they do not exist)"""
    dirpath = Path(dirpath)
    dirpath.mkdir(parents=True, exist_ok=True)

def rmdir(dirpath, not_found_ok=False):
    """Remove a directory and everyting in it.

    Args:
        dirpath (string or pathlib.Path): directory path.
    """
    dirpath = Path(dirpath)
    if not_found_ok:
        if dirpath.exists() == False:
            return            
    for item in dirpath.iterdir():
        if item.is_dir():
            rmdir(item)
        else:
            item.unlink()
    dirpath.rmdir()

def rm(filepath):
    """Remove a file.

    Args:
        filepath (string or pathlib.Path): filepath.
    """
    filepath = Path(filepath)

    if filepath.is_dir(): 
        raise IsADirectoryError('Is a directory: `{filepath}`\nPlease, use br.rmdir() for directories')
    filepath = Path(filepath).unlink()

# %% ========================= file/folder information ==================== %% #
def filelist(dirpath='.', string='*', case_sensitive=True):
    """Returns a list with all the files containing `string` in its name.

    Note:
        List is sorted by the filename.

    Args:
        dirpath (str or pathlib.Path, optional): list with full file directory
            paths.
        string (str, optional): pattern. only filenames with this string will be considered.
            Use '*' for matching anything. Default is '*'.
        case_sensitive (bool or None, optional): Default is True. If None, 
            case_sensitive will matches paths using platform-specific casing rules.
            Typically, case-sensitive on POSIX, and case-insensitive on Windows.

    Return:
        list

    See Also:
        :py:func:`parsed_filelist`
    """
    dirpath = Path(dirpath)

    if '*' not in string:
        string = '*' + string + '*'

    # on linux and mac, glob is naturally case sensitive
    if (is_linux or is_mac) and case_sensitive == False:
        rule = re.compile(fnmatch.translate(string), re.IGNORECASE)
        temp = [dirpath/name for name in os.listdir(dirpath) if rule.match(name)]
    # on windows, glob is naturally case INsensitive
    elif is_windows and case_sensitive:
        temp = list(dirpath.glob(pattern=string))
        match = re.compile(fnmatch.translate(str(dirpath/string))).match
        temp  = [path for path in temp if match(str(path))]
    else:
        temp = list(dirpath.glob(pattern=string))

    # sort and return
    temp2 = [filepath.name for filepath in temp]
    return [x for _,x in sorted(zip(temp2, temp))]

def parsed_filelist(dirpath='.', string='*', ref=0, _type='int', return_type='list'):
    """Returns a filelist organized in a dictionary.

    I searches for numbers (float and int) within the filenames (or foldernames)
    and uses them as keys for the dictionary.

    Args:
        dirpath (str or pathlib.Path, optional): directory path.
        string (str, optional): string to filter filenames by.
        ref (int, optional): index of the reference number to be used as key.
        _type (string, optional): if 'int', dict keys are transformed in integers.
            If 'float', dict keys are transformed into float.

    Returns:
        if return type is 'list'
            returns a list
        if return type is 'dict'
            ordered Dictionary. Dict keys are some number found in filename.

    See Also:
        :py:func:`filelist`
    """
    dirpath = Path(dirpath)

    file_list = filelist(dirpath=dirpath, string=string)

    temp = dict()

    p = r'[\d]+[.,\d]+|[\d]*[.][\d]+|[\d]+'

    for filepath in file_list:
        if re.search(p, filepath.name) is not None:
            n = []
            for catch in re.finditer(p, filepath.with_suffix('').name):
                n.append(catch[0])
            if _type=='int':
                temp[int(float((n[ref])))] = filepath
            else:
                temp[float(n[ref])] = filepath


    # ordering
    a = list(temp.keys())
    a.sort()

    parsed_folder = collections.OrderedDict()
    for key in a:
        parsed_folder[key] = temp[key]

    if return_type == 'list':
        return list(parsed_folder.values())
    elif return_type == 'dict':
        return parsed_folder
    else:
        raise ValueError('return_type not valid')

def get_modified_date(filepath):
    """Return modified date of filepath in a datetime object."""
    return datetime.datetime.fromtimestamp((filepath).stat().st_mtime)

# %% ============================= save/load ============================== %% #
def save_text(string, filepath='./Untitled.txt', check_overwrite=False, **kwargs):
    """Save text to txt file.

    Args:
        string (str): string to be saved.
        filepath (str or pathlib.Path, optional): path to save file. If no path
            is given, current working directory is used.
        check_overwrite (bool, optional): if True, it will check if file exists
            and ask if user want to overwrite file.
        **kwargs are passed to open(filepath, 'w', **kwargs) function.

    Returns:
        None:

    See Also:
        :py:func:`load_text`
    """
    filepath = Path(filepath)

    if check_overwrite:
        if filepath.exists() == True:
            if filepath.is_file() == True:
                if _query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                    pass
                else:
                    warnings.warn('File not saved.')
                    return
            else:
                warnings.warn('filepath is pointing to a folder. Saving file as Untitled.txt')
                filepath = filepath/'Untitled.txt'

    f = open(str(filepath), 'w', **kwargs)
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

def save_obj(obj, filepath='./Untitled.txt', check_overwrite=False, pretty_print=True):
    """Save object (array, dictionary, list, etc...) to a txt file.

    Args:
        obj (object): object to be saved.
        filepath (str or pathlib.Path, optional): path to save file.
        check_overwrite (bool, optional): if True, it will check if file exists
            and ask if user want to overwrite file.

    Returns:
        None:

    See Also:
        :py:func:`load_obj`
    """
    filepath = Path(filepath)

    if check_overwrite:
        if filepath.exists() == True:
            if filepath.is_file() == True:
                if _query('File already exists!! Do you wish to ovewrite it?', 'yes') == True:
                    pass
                else:
                    warnings.warn('File not saved because user did not allow overwriting.')
                    return
            else:
                warnings.warn('filepath is pointing to a folder. Saving file as Untitled.txt')
                filepath = filepath/'Untitled.txt'

    with open(str(filepath), 'w') as file:
        if pretty_print:
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

def load_comments(filepath, comment_flag='#', stop_flag='#'):
    """Return comments from text files.

    Comments must be indicated at the beginning of the line by the comment flag.

    Args:
        filepath (str or pathlib.Path): fullpath to file
        comment_flag (str, optional): string that indicate line is a comment.
        stop_flag (str, optional): string that indicates line to stop looking for
            comments, e.g. `#f`, or `#L`. Use `None` to read all lines in file (useful for reading
            file with comments not only at the beginning). If `stop_flag` is
            equal to `comment_flag` it will read from the first line with
            `comment_flag` and keep reading until `comment_flag` does not appear
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
        with open(filepath) as file:
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

def save_data(obj, filepath='./untitled.txt', add_labels=True, fmt='% .10e', header='', footer='', delimiter=',', comment_flag='# ', newline='\n', check_overwrite=False):
    r"""Save an array or a dictionary in a txt file.

    Args:
        obj (dict, list, or numpy.array): data to be saved to a file. If obj is
        a dictionary, use ``*`` in front of a key to do not save it to the file.
        filepath (str or pathlib.Path, optional): path to save file.
        add_labels (bool, optional): When obj is a dictionary, ``add_labels=True``
            makes the dict keys to be added to the header as label for each data column.
        fmt (string, or list, optional): format for saving data.
            If string, the value is used for x- and y-coordinates. If tuple
            of strings, the first string is used for x-coordinates and the
            second for y-coordinates.

            ``fmt = (%[flag]width[.precision]specifier)``

            * `flag` can be:

                1. '-' for left justify

                2. '+', puts + or - in front of the numbers

                3. '0' to Left pad the number with zeros instead of space (see width).

            * `width` is the minimum number of characters to be printed.

            * `precision` is the number of significant digits.

            * `specifier` is the type of notation. Typically, either 'e' for scientific notation of 'f' for decimal floating point.

            * a common `fmt` strings is: '%.3f' for 3 decimal places.

            *  for more information see `np.savetxt <https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html?highlight=savetxt#numpy.savetxt>`_ documentation::

        header (str, optional): string that will be written at the beginning of
            the file (comment flag is added automatically).
        footer (str, optional): string that will be written at the end of the file (comment flag is added automatically).
        delimiter (str, optional): The string used to separate values.
        comment_flag (str, optional): string that flag comments.
        newline (str, optional): string to indicate new lines.
        check_overwrite (bool, optional): if True, it will check if file exists
            and ask if user want to overwrite file.

    Returns:
        None:

    See Also:
        :py:func:`load_data`
    """
    filepath = Path(filepath)

    if check_overwrite:
        if filepath.exists() == True:
            if filepath.is_file() == True:
                if _query('File already exists!! Do you wish to ovewrite it?', 'yes') == True:
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
        if add_labels:
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

    np.savetxt(filepath, obj, fmt=fmt, delimiter=delimiter, newline=newline, header=header, footer=footer, comments=comment_flag)

def load_data(filepath, labels=None, force_array=False, **kwargs):
    """Load data from text file. Wrapper for `np.genfromtxt <https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html>`_.

    Data is formatted in a dictionary or array.

    The dictionary keys are set as the label of the corresponding data columns, where
    the last comment line before data starts is assumed to have the labels of each data column.
    If column labels cannot be found, data is imported as an array. The file expects
    comments at the beginning of the file (file must not have comments elsewhere).

    If a column is full of NaN's, it's converted to a column full of zeros.

    Args:
        filepath (str or pathlib.Path): path to file
        labels (list, optional): It forces data to be loaded as a dictionary where
            each label is associated with a data column. Its length must have the same as the number of
            columns. To avoid importing a column, put an asterisk (*) in front of the corresponding label.
        force_array (bool, optional): If ``force_array=True``, data it will be 
            returned in a array. Overwrites labels.

    If not specified, the following parameters are passed to `np.savetxt <https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html?highlight=savetxt#numpy.savetxt>`_:

    Args:
        delimiter (str, optional): String or character separating columns.
            Use ``\\t`` for tab. If None,
            the script will read the file and try to guess the appropriate delimiter.
            If the file has a header, and if the header has a line with the column labels, it
            tries to guess the delimiter of this line. If it cannot, it tries
            to use the same delimiter for the data and for the leader in the header.
            Default is None.
        comments (str, optional): The character used to indicate the start of a comment.
            Default is '#'.
        skip_header (int, optional): The number of lines to skip at the beginning of the file.
            Default is 0.


    Returns:
        Dictionary or array.

    See Also:
        :py:func:`save_data`.
    """
    if 'delimiter' not in kwargs:
        kwargs['delimiter'] = None
    # if 'newline' not in kwargs:
    #     kwargs['newline'] = '\n'
    if 'comments' not in kwargs:
        kwargs['comments'] = '#'
    if 'skip_header' not in kwargs:
        kwargs['skip_header'] = 0

    filepath = Path(filepath)

    # guess delimiter
    # if kwargs['delimiter'] is None:
    #     # find header length
    #     header = load_Comments(filepath, comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
    #     try: line_number = len(header)
    #     except TypeError: line_number = 0
    #     if kwargs['skip_header'] > line_number:
    #         line_number = kwargs['skip_header']

    #     # find delimiter from the first data row
    #     f = open(filepath)
    #     for i, line in enumerate(f):
    #         if i == line_number:
    #             kwargs['delimiter'] = detect(line)
    #         elif i > line_number:
    #             break
    #     f.close()
    #     if kwargs['delimiter'] is None:
    #         warnings.warn('Could not figure out the delimiter. Trying space.')
    # if kwargs['delimiter'] == ' ':
    #     kwargs['delimiter'] = None

    # get data
    # print(kwargs)
    data = np.genfromtxt(str(filepath), **kwargs)

    # if a column returns only nan, this column yield 0
    for i in range(data.shape[1]):
        if False in [x for x in data[:, i] != data[:, i]]:  # test if the whole col is nan
            pass
        else:
            # key conversion to col number
            temp = np.genfromtxt(filepath, usecols=i, dtype='S8', **kwargs)
            # data[:, i] = [x.decode('utf-8') for x in temp]
            data[:, i] = [0 for x in temp]

    # # get labels
    # if labels is None and force_array is False:
    #     # print('here')
    #     header = load_Comments(filepath, comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
    #     if header:
    #         header_line = header[-1].replace(kwargs['comments'], '').strip()
    #         header_line = header_line.replace('\n', '')
    #         if header_delimiter is None:
    #             header_delimiter = detect(header_line)
    #         labels = header_line.split(header_delimiter)
    # # print(labels)

    # check labels length
    if labels is not None and force_array is False:
        assert len(labels) == data.shape[1], 'labels must be the same length as the number of columns'
        labels = [item.strip() for item in labels if item != '']  # remove empty items and trailing spaces
        datadict = {labels[i]: data[:, i] for i in range(len(labels)) if not labels[i].startswith('*')}
        return datadict
    else:
        return data

