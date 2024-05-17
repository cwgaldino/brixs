#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for manipulating Libreoffice Calc files.

- Install libreoffice. Follow the instructions.
In this case we are in Ubuntu 18.04 and libreoffice version 7.0.4

tar -xvzf file
cd folder
cd DEBS
sudo dpkg -i *.deb


sudo apt-get install python3-uno  >> import uno
pip install unotools >>>
from unotools import Socket, connect
from unotools.component.calc import Calc
from unotools.unohelper import convert_path_to_url

TO DO:
    - method to delete rows and cols
    - method in sheet: get sheet name
    - conditional formatting is not wroking (it seems)

Used methods from unotools.component.calc.Calc:
    # 'set_rows_str', OK---
    # 'set_columns_str', OK---
    # 'set_rows_value', OK---
    # 'set_columns_value', OK---
    # 'set_rows_formula', OK---
    # 'set_columns_formula', OK---
    # 'get_cell_by_position', OK---
    # 'get_cell_range_by_position', OK---

Not used:
    # 'set_rows_cell_data',
    # 'set_columns_cell_data',
    # 'set_rows',
    # 'set_columns',
    # 'get_cell_range_by_name',
    # 'charts',
    # 'get_charts_count',
    # 'add_charts_new_by_name',
    # 'get_chart_by_index',
    # 'get_chart_by_name',

"""

# standard imports
import numpy as np
from pathlib import Path
import os
import inspect
import psutil
import signal
import subprocess
import time
import warnings
import re
import copy
from collections.abc import Iterable

# import uno
import sys
from unotools import Socket, connect
from unotools.component.calc import Calc
from unotools.unohelper import convert_path_to_url

from .intermanip import query_yes_no


# %%

class soffice():

    def __init__(self, port=8100, norestore=False):
        self.pid_previous = self._libreoffice_pid_list()
        self.port = port
        if norestore:
            self.process = subprocess.Popen([f"soffice --nodefault --norestore --nologo --accept='socket,host=localhost,port={port};urp;'"], shell=True, close_fds=True)
        else:
            self.process = subprocess.Popen([f"soffice --nodefault --nologo --accept='socket,host=localhost,port={port};urp;'"], shell=True, close_fds=True)

        time.sleep(1)
        self.pid = self.process.pid
        self.pid_children = self._get_children_pid(self.pid)
        self.apps = []


    def openCalc(self, filepath=None):
        children_old = self._get_children_pid(self.pid)

        calcObject = calc(filepath, self)
        calcObject.pid = [child for child in self._get_children_pid(self.pid) if child not in children_old]

        self.pid_children += calcObject.pid
        # self.pid_children = [item for sublist in self.pid_children for item in sublist]
        self.apps.append(calcObject)

        return calcObject


    def _get_children_pid(self, pid):
        output = subprocess.check_output(["bash", "-c", f"pstree -p -n  {pid}"])
        pattern = re.compile(r"\((\d+)\)")
        pid_children = pattern.findall(output.decode('utf-8'))
        return [int(pid) for pid in pid_children]


    def _get_pid_by_name(self, string):
        """Get a list of all the PIDs of all the running process whose name contains
        string.

        Args:
            string (str): string.
        """

        listOfProcessObjects = []

        # Iterate over the all the running process
        for proc in psutil.process_iter():
           try:
               pinfo = proc.as_dict(attrs=['pid', 'name', 'create_time'])
               # Check if process name contains the given name string.
               if string.lower() in pinfo['name'].lower() :
                   listOfProcessObjects.append(pinfo)
           except (psutil.NoSuchProcess, psutil.AccessDenied , psutil.ZombieProcess):
               pass

        return [item['pid'] for item in listOfProcessObjects];


    def _libreoffice_pid_list(self):
        """Return a list of processes associated with libreoffice.

        Note:
            This function try to match the name of a process with names tipically
            related with libreoffice processes ('soffice.bin' or 'oosplash').
            Therefore, it might return processes that are not related to
            libreoffice if their name mathces with words: 'soffice'
            and 'oosplash'.

        Returns:
            list.
        """
        process_list = []
        for proc in self._get_pid_by_name('soffice'):
            process_list.append(proc)
        for proc in self._get_pid_by_name('oosplash'):
            process_list.append(proc)
        return [i for sublist in [self._get_children_pid(item) for item in process_list] for i in sublist]


    def kill_libreoffice_processes(self):
        """Kill libreoffice processes.

        Note:
            It will close ALL processes that are related to libreoffice (processes
            that have 'soffice.bin' or 'oosplash' in their name)."""

        pid_list = self._libreoffice_pid_list()
        for proc in process_list:
            os.kill(proc['pid'], signal.SIGKILL)


    def terminate(self, ask=True):

        if ask:
            if len(self.apps)>0:
                print('You appear to have opened libreoffice apps.')
                for app in self.apps:
                    if app.filepath is None:
                        msg = f'{app.object} is not associated with a file. Do you wish to discard this work?'
                    else:
                        msg = f'Progress not saved on {app.filepath} will be discarted. Do you wish to continue?'
                    if not query_yes_no(msg, 'no'):
                        return
            msg = 'All instances of libreoffice (even if not initialized via python) will be discarted and progress not saved will be lost. Do you wish to continue?'
            if not query_yes_no(msg, 'no'):
                return

        for pid in self._libreoffice_pid_list():
            try:
                os.kill(int(pid), 9)
                print(f'{int(pid)} killed.')
            except ProcessLookupError:
                print(f'{int(pid)} does not exist.')
        print('Done!')


class calc():


    def __init__(self, filepath=None, parent=None):
        self.pid = []
        self.object = None
        self.object_parent = parent


        context = connect(Socket('localhost', f'{self.object_parent.port}'))
        if filepath is None:
            self.object = Calc(context)
            self.filepath = None
        else:
            self.filepath = Path(filepath).absolute()
            self.object = Calc(context, convert_path_to_url(str(self.filepath)))


    def save(self, filepath=None):
        """Save xlsx file.

        Note:
            If filepath have no suffix, it adds '.ods' at the end of the filename.

        Args:
            filepath (string or pathlib.Path, optional): filepath to save file.
        """
        if filepath is None and self.filepath is None:
            temporary_path = Path.cwd()/'Untitled.ods'
            if not query_yes_no(f'Filepath not defined. Wish to save at {temporary_path}?'):
                filepath = temporary_path
                return
        elif filepath is not None:
            self.filepath = Path(filepath)


        self.filepath = Path(self.filepath)

        # fix extension
        if self.filepath.suffix != '.ods':
            self.filepath = Path(self.filepath).with_suffix('.ods')

        # save
        url = convert_path_to_url(str(self.filepath))
        self.object.store_as_url(url, 'FilterName')
        print(f'Saved at: {url}')


    def close(self):
        """Close window."""

        self.object.close(True)


    def terminate(self, ask=True):
        self.object_parent.terminate(ask)


    def get_sheets_count(self):
        return self.object.get_sheets_count()


    def get_sheets_name(self):
        """Returns the sheets names in a tuple."""
        return self.object.Sheets.ElementNames


    def insert_sheets(self, name, position=None):
        """position starts from 1. If position = 1, the sheet will be the first one.
        """
        if position is None:
            position = self.get_sheets_count()+1

        if type(name) == str:
            name = [name]

        existing_names = [name2 for name2 in name if name2 in self.get_sheets_name()]
        if len(existing_names) != 0:
            raise SheetNameExistError(existing_names)

        self.object.insert_multisheets_new_by_name(name, position-1)


    def remove_sheets(self, name):
        self.remove_sheets_by_name(name)


    def remove_sheets_by_name(self, name):

        if type(name) == str:
            name = [name]

        not_existing_names = [name2 for name2 in name if name2 not in self.get_sheets_name()]
        if len(not_existing_names) != 0:
            raise SheetNameDoNotExistError(not_existing_names)

        for n in name:
            if len(self.get_sheets_name()) == 1:
                raise SheetRemoveError(n)
            self.object.remove_sheets_by_name(n)


    def remove_sheets_by_position(self, position):
        names = self.get_sheets_name()

        if position > len(names) or position < 1:
            raise IndexError('Position outside range.')

        if len(names) == 1:
            raise SheetRemoveError(names[position-1])

        self.object.remove_sheets_by_name(names[position-1])


    def get_sheets(self, name):
        if type(name) == str:
            name = [name]

        if type(name) == int:
            name = [name]

        sheet_objects = []
        for n in name:
            if type(n) == str:
                sheet_objects.append(self.get_sheets_by_name(n))
            if type(n) == int:
                sheet_objects.append(self.get_sheets_by_position(n))
        return [next(obj) for obj in sheet_objects]


    def get_sheets_by_name(self, name):
        if type(name) != list:
            name = [name]

        for n in name:
            yield sheet(self.object.get_sheet_by_name(str(n)), self)


    def get_sheets_by_position(self, position):
        names = self.get_sheets_name()

        if type(position) == int:
            position = [position]

        outside_range = []
        for p in position:
            if p > len(names) or p < 1:
                outside_range.append(p)
        if len(outside_range) > 0:
            raise IndexError(f'Positions {outside_range} outside range.')

        for p in position:
            yield sheet(self.object.get_sheet_by_index(p-1), self)


class SheetNameExistError(Exception):

    # Constructor or Initializer
    def __init__(self, sheet_names):
        if type(sheet_names) == str:
                sheet_names = [sheet_names]
        self.names = sheet_names

    # __str__ is to print() the value
    def __str__(self):
        msg = ''
        for name in self.names:
            msg += f"'{name}' already exists\n"
        return(msg)


class SheetNameDoNotExistError(Exception):

    # Constructor or Initializer
    def __init__(self, sheet_names):
        if type(sheet_names) == str:
                sheet_names = [sheet_names]
        self.names = sheet_names

    # __str__ is to print() the value
    def __str__(self):
        msg = ''
        for name in self.names:
            msg += f"'{name}' does not exists\n"
        return(msg)


class SheetRemoveError(Exception):

    # Constructor or Initializer
    def __init__(self, sheet_names):
        if type(sheet_names) == str:
                sheet_names = [sheet_names]
        self.names = sheet_names

    # __str__ is to print() the value
    def __str__(self):
        msg = ''
        for name in self.names:
            msg += f"'{name}' cannot be removed because it is the only existing sheet.\n"
        return(msg)


class sheet():

    def __init__(self, object, object_parent):
        self.object = object
        self.object_parent = object_parent

    def get_name(self):
        return self.object.getName()

    def set_name(self, name):
        return self.object.setName(name)

    def get_last_row(self):
        # higher_edited_row = len(self.object.getRowDescriptions()) + self.object.queryVisibleCells().Count -1
        # higher_edited_col = len(self.object.getColDescriptions()) + self.object.queryVisibleCells().Count -1
        #
        # row = higher_edited_row
        # if not any(v for v in self.get_row_values(row=row, col_stop=higher_edited_col)):
        #
        #
        # higher_edited_row
        idx = int(self.object.getRowDescriptions()[-1].split(' ')[-1])
        visible = False
        while visible == False:
            visible = self.object.getRows().getByIndex(idx).IsVisible
            idx += 1

        return idx-1


    def get_last_col(self):
        idx = _letter2num(self.object.getColumnDescriptions()[-1].split(' ')[-1])
        visible = False
        while visible == False:
            visible = self.object.getColumns().getByIndex(idx).IsVisible
            idx += 1


        return idx-1


    def set_col_width(self, width, col):
        col = _check_col_value(col)

        colsObject = self.object.getColumns()
        for c in col:
            colsObject[c].setPropertyValue('Width', width)


    def get_col_width(self, col):

        col = _check_col_value(col)

        colsObject = self.object.getColumns()

        width = []
        for c in col:
            width.append(colsObject[c].Width)
        return width


    def set_row_height(self, height, row):

        row = _check_row_value(row)

        rowsObject = self.object.getRows()
        for r in row:
            rowsObject[r].setPropertyValue('Height', height)


    def get_row_height(self, row):
        row = _check_row_value(row)

        rowsObject = self.object.getRows()

        height = []
        for r in row:
            height.append(rowsObject[r].Height)
        return height


    def set_cell_value(self, value, row, col, format='formula'):
        """
        format='data', 'formula'
        """
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        if format == 'formula':
            self.object.get_cell_by_position(col, row).setFormula(value)
        elif format == 'string':
            self.object.get_cell_by_position(col, row).setString(value)
        elif format == 'number':
            self.object.get_cell_by_position(col, row).setValue(value)
        else:
            raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")


    def get_cell_value(self, row, col, format='string'):
        """
        format='data', 'formula'
        """
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        if format == 'formula':
            value = self.object.get_cell_by_position(col, row).getFormula()
        elif format == 'string':
            value = self.object.get_cell_by_position(col, row).getString()
        elif format == 'number':
            value = self.object.get_cell_by_position(col, row).getValue()
        else:
            raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")

        return value


    def get_cells_value(self, row_start=1, col_start=1, row_stop=None, col_stop=None, format='string'):
        """
        format= formula or data.
        """
        row_start = _check_row_value(row_start)[0]
        col_start = _check_col_value(col_start)[0]

        if row_stop is None:
            row_stop = self.get_last_row()
        if col_stop is None:
            col_stop = self.get_last_col()

        row_stop = _check_row_value(row_stop)[0]
        col_stop = _check_col_value(col_stop)[0]

        if col_stop < col_start:
            raise ValueError('col_start cannot be bigger than col_stop')
        if row_stop < row_start:
            raise ValueError('row_start cannot be bigger than row_stop')

        sheet_data = self.object.get_cell_range_by_position(col_start, row_start, col_stop, row_stop)
        if format == 'formula':
            sheet_data = list(sheet_data.getFormulaArray())
        elif format == 'string':
            sheet_data = list(sheet_data.getDataArray())
        elif format == 'number':
            sheet_data = list(sheet_data.getDataArray())
        else:
            raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")

        # transform in list
        for row_number, row_data in enumerate(sheet_data):
            sheet_data[row_number] = list(row_data)

            # if one column or one row data, transform in vector
            if col_start == col_stop:
                sheet_data[row_number] = row_data[0]
        if row_start == row_stop:
            sheet_data = sheet_data[0]

        return sheet_data


    def get_row_values(self, row, col_start=1, col_stop=None, format='string'):
        return self.get_cells_value(row_start=row, col_start=col_start, row_stop=row, col_stop=col_stop, format=format)


    def get_col_values(self, col, row_start=1, row_stop=None, format='string'):
        return self.get_cells_value(row_start=row_start, col_start=col, row_stop=row_stop, col_stop=col, format=format)


    def set_col_values(self, data, col, row_start=1, format='formula'):
        row_start = _check_row_value(row_start)[0]
        col = _check_col_value(col)[0]

        if type(data)==np.ndarray:
            data = data.tolist()

        if format == 'formula':
            self.object.set_rows_formula(col, row_start, data)
        elif format == 'string':
            self.object.set_rows_str(col, row_start, data)
        elif format == 'number':
            self.object.set_rows_value(col, row_start, data)
        else:
            raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")


    def set_row_values(self, data, row, col_start=1, format='formula'):
        row = _check_row_value(row)[0]
        col_start = _check_col_value(col_start)[0]

        if type(data)==np.ndarray:
            data = data.tolist()

        if format == 'formula':
            self.object.set_columns_formula(col_start, row, data)
        elif format == 'string':
            self.object.set_columns_str(col_start, row, data)
        elif format == 'number':
            self.object.set_columns_value(col_start, row, data)
        else:
            raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")


    def set_cells_value(self, data, row_start=1, col_start=1, format='formula'):
        """
        if data is 1d array or list, data is placed in a row.

        formula set formulas, but also set numbers fine. Dates and time not so much because it changes the formating (if setting date and time iwth formula you might wanna format the
        cell like date or time using copy_cells to copy formatting).

        string (data) works fine with date, time and number, but formulas are set as string. Therefore, formulas do not work.

        value (data_number) works fine for numbers ONLY.

        """
        row_start = _check_row_value(row_start)[0]
        col_start = _check_col_value(col_start)[0]

        try:
            row_count, col_count = np.shape(data)
        except ValueError:
            col_count = len(data)
            row_count = 0

        if type(data) == np.ndarray:
            if row_count >= col_count:
                for idx in range(col_count):
                    if format == 'formula':
                        self.set_col_values(data=data[:, idx], col=col_start+idx+1, row_start=row_start+1, format='formula')
                    elif format == 'string':
                        self.set_col_values(data=data[:, idx], col=col_start+idx+1, row_start=row_start+1, format='string')
                    elif format == 'number':
                        self.set_col_values(data=data[:, idx], col=col_start+idx+1, row_start=row_start+1, format='number')
                    else:
                        raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")

            else:
                for idx in range(row_count):
                    if format == 'formula':
                        self.set_row_values(data=data[idx, :], row=row_start+idx+1, col_start=col_start+1, format='formula')
                    elif format == 'string':
                        self.set_row_values(data=data[idx, :], row=row_start+idx+1, col_start=col_start+1, format='string')
                    elif format == 'number':
                        self.set_row_values(data=data[idx, :], row=row_start+idx+1, col_start=col_start+1, format='number')
                    else:
                        raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")
        else:
            for idx in range(row_count):
                if format == 'formula':
                    self.set_row_values(data=data[idx], row=row_start+idx+1, col_start=col_start+1, format='formula')
                elif format == 'string':
                    self.set_row_values(data=data[idx], row=row_start+idx+1, col_start=col_start+1, format='string')
                elif format == 'number':
                    self.set_row_values(data=data[idx], row=row_start+idx+1, col_start=col_start+1, format='number')
                else:
                    raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")


    def list_properties(self):
        return self.object._show_attributes()


    def list_cell_properties(self, filter=None, row=1, col=1):
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        cellObject = self.object.get_cell_by_position(col, row)
        if filter is None:
            # return cellObject._show_attributes()
            return cellObject.as_raw().__dir__()
        else:
            p = cellObject.as_raw().__dir__()
            matching = [s for s in p if filter.lower() in s.lower()]
            return matching


    def get_cell_property(self, property, row, col):
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        if type(property) == str:
            property = [property]

        object_cell = self.object.get_cell_by_position(col, row)
        attr = getattr(object_cell, property[0])

        for idx in range(1, len(property)):
            attr = getattr(attr, property[idx])

        # if type(attr) != int and type(attr) != float and type(attr) != bool:
        #     return attr, attr.value.__dir__()
        # else:
        #     return attr, None
        try:
            return attr, attr.value.__dir__()
        except AttributeError:
            return attr, None


    def set_cell_property(self, property, value, row, col):
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        if type(property) == str:
            property = [property]

        self.object.get_cell_by_position(col, row).setPropertyValue(property[0], value)


    def get_cells_properties(self, property, row_start, col_start, row_stop, col_stop):
        row_start = _check_row_value(row_start)[0]
        col_start = _check_col_value(col_start)[0]
        row_stop = _check_row_value(row_stop)[0]
        col_stop = _check_col_value(col_stop)[0]

        if type(property) == str:
            property = [property]

        object_cell = self.object.get_cell_range_by_position(col_start, row_start, col_stop, row_stop)
        attr = getattr(object_cell, property[0])

        for idx in range(1, len(property)):
            attr = getattr(attr, property[idx])

        try:
            return attr, attr.value.__dir__()
        except AttributeError:
            return attr, None


    def set_cells_properties(self, property, value, row_start, col_start, row_stop, col_stop):
        row_start = _check_row_value(row_start)[0]
        col_start = _check_col_value(col_start)[0]
        row_stop = _check_row_value(row_stop)[0]
        col_stop = _check_col_value(col_stop)[0]

        if type(property) == str:
            property = [property]

        self.object.get_cell_range_by_position(col_start, row_start, col_stop, row_stop).setPropertyValue(property[0], value)


    def get_cell_object(self, row=1, col=1):
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        return self.object.get_cell_by_position(col, row)


    def get_cell_formating(self, row, col, extra=None):
        """font : bold, font name, font size, italic, color
        text: vertical justify,horizontal justify
        cell: border, color
        conditional formating: conditional formating
        """

        p = ['FormatID']
        p += ['CharWeight', 'CharFontName', 'CharHeight', 'CharPosture', 'CharColor']
        p += ['VertJustify', 'HoriJustify']
        p += ['CellBackColor', 'TableBorder', 'TableBorder2']
        p += ['ConditionalFormat']

        if extra is not None:
            p += extra

        p_list = []
        for property in p:
            obj, _ = self.get_cell_property(property, row=row, col=col)
            p_list.append(obj)

        return p_list


    def set_cell_formating(self, obj_list, row, col, extra=None):
        """font : bold, font name, font size, italic, color
        text: vertical justify,horizontal justify
        cell: border, color
        conditional formating: conditional formating
        """
        p = ['FormatID']
        p += ['CharWeight', 'CharFontName', 'CharHeight', 'CharPosture', 'CharColor']
        p += ['VertJustify', 'HoriJustify']
        p += ['CellBackColor', 'TableBorder', 'TableBorder2']
        p += ['ConditionalFormat']

        if extra is not None:
            p += extra

        for obj, property in zip(obj_list, p):
            self.set_cell_property(property, value=obj, row=row, col=col)


    def get_cells_formatting(self, row_start=1, col_start=1, row_stop=None, col_stop=None, extra=None):
        if row_stop is None:
            row_stop = self.get_last_row()
        if col_stop is None:
            col_stop = self.get_last_col()

        row_stop = _check_row_value(row_stop)[0]
        col_stop = _check_col_value(col_stop)[0]

        cell_formatting_list = []
        for idx, row in enumerate(range(row_start, row_stop+2)):
            cell_formatting_list.append([])
            for col in range(col_start, col_stop+2):
                p_list = self.get_cell_formating(row=row, col=col, extra=extra)
                cell_formatting_list[idx].append(p_list)
        return cell_formatting_list


    def set_cells_formatting(self, cell_formatting_list, row_start=1, col_start=1, extra=None):

        for idx, row in enumerate(range(row_start, row_start+len(cell_formatting_list))):
            for idx2, col in enumerate(range(col_start, col_start+len(cell_formatting_list[idx]))):
                self.set_cell_formating(cell_formatting_list[idx][idx2], row=row, col=col, extra=extra)


    def get_merged(self, ):
        merged_ranges = []
        #################
        start=[]
        stop = []
        #################

        ucf = self.object.getUniqueCellFormatRanges()
        for ranges in ucf:
            rgtest = ranges.getByIndex(0)
            if rgtest.getIsMerged():
                for rg in ranges:
                    oCursor = rg.getSpreadsheet().createCursorByRange(rg)
                    oCursor.collapseToMergedArea()
                    ######################
                    row_start = int(oCursor.getRowDescriptions()[0].split(' ')[-1])-1
                    row_end = int(oCursor.getRowDescriptions()[-1].split(' ')[-1])-1
                    for row in range(row_end-row_start+1):
                        col_start = backpack.libremanip._letter2num(oCursor.getColumnDescriptions()[0].split(' ')[-1])-1
                        col_end = backpack.libremanip._letter2num(oCursor.getColumnDescriptions()[-1].split(' ')[-1])-1
                        for col in range(col_end-col_start+1):
                            if oCursor.getCellByPosition(col, row).IsMerged:
                                # print(row+row_start, col+col_start)
                                start.append([row+row_start, col+col_start])
                            else:
                                stop.append([row+row_start, col+col_start])
                                # print(f'not: {row+row_start}, {col+col_start}')
                    ######################
                    addr = oCursor.getRangeAddress()

                    col_start = addr.StartColumn+1
                    row_start = addr.StartRow+1
                    col_stop = addr.EndColumn+1
                    row_stop = addr.EndRow+1
                    merged_ranges.append([row_start, col_start, row_stop, col_stop])
        return merged_ranges


    def merge(self, row_start, col_start, row_stop, col_stop):
        row_start = _check_row_value(row_start)[0]
        col_start = _check_col_value(col_start)[0]
        row_stop = _check_row_value(row_stop)[0]
        col_stop = _check_col_value(col_stop)[0]

        sheet_data = self.object.get_cell_range_by_position(col_start, row_start, col_stop, row_stop)
        sheet_data.merge(True)


def _iterable(obj):
    return isinstance(obj, Iterable)


def _letter2num(string):

    string = string.lower()
    alphabet = 'abcdefghijklmnopqrstuvwxyz'

    col = 0
    for idx, s in enumerate(string):
        col += alphabet.index(s)+(idx)*26+1

    return(col)


def _check_row_value(row):
    if _iterable(row) == False:
        row2 = [row]
    else:
        row2 = copy.deepcopy(row)

    row2 = [r-1 for r in row2]

    outside_range = []
    for r in row2:
        if r < 0:
            outside_range.append(r)
    if len(outside_range) > 0:
        raise IndexError(f'Row number {outside_range} outside range.')

    return row2


def _check_col_value(col):
    if _iterable(col)==False or type(col)==str:
        col2 = [col]
    else:
        col2 = copy.deepcopy(col)

    for idx, c in enumerate(col2):
        if type(c) == str:
            col2[idx] = _letter2num(c)-1
        else:
            col2[idx] = c-1
            # raise TypeError(f'Col value {c} must be an int or str.')

    outside_range = []
    for c in col2:
        if type(c) == int:
            if c < 1:
                outside_range.append(c)
    if len(outside_range) > 1:
        raise IndexError(f'Col number {outside_range} outside range.')

    return col2


def get_cell_value_from_sheets(sheetObject_list, row, col, format='string'):
    """
    """
    values = []
    for sheetObject in sheetObject_list:
        values.append(sheetObject.get_cell_value(row=row, col=col, format=format))
    return values
