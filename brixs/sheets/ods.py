#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quality-of-life functions for libreoffice Calc macros

Usage:
    >>> from brixs.sheets.ods import get_sheet_count, get_sheet_names 
    >>> from brixs.sheets.ods import get_sheet_by_position, get_sheet_by_name, get_active_sheet
    >>> from brixs.sheets.ods import new_sheet, remove_sheet_by_name, remove_sheet_by_position, move_sheet, copy_sheet
    >>> from brixs.sheets.ods import get_style_names, remove_style, new_style, select_active_sheet
    >>> from brixs.sheets.ods import get_cells, get_last_used_row, get_last_used_col
    >>> from brixs.sheets.ods import delete_rows
    >>> from brixs.sheets.ods import merge, unmerge
    >>> from brixs.sheets.ods import get_row, get_row_data, get_col_data, set_row_data, get_cells_data, set_cells_data
    >>> from brixs.sheets.ods import get_property_names, get_property_value, set_property_value
    >>> from brixs.sheets.ods import get_width, get_height, set_width, set_height
    >>> from brixs.sheets.ods import get_border, set_border
    >>> from brixs.sheets.ods import get_protection, lock_cells, unlock_cells
"""
# %% ------------------------- standard imports --------------------------- %% #
pass

# %% ---------------------------- uno imports ----------------------------- %% #
import uno

# %% -------------------------- brixs imports ----------------------------- %% #
from brixs.sheets.common import str2num
# %%

# %% ------------------------ calc document functions --------------------- %% #
def get_sheet_count(doc):
    """Return number of sheets"""
    return doc.Sheets.getCount()

def get_sheet_names(doc):
    """Return sheet names in order"""
    return doc.Sheets.getElementNames()

def get_sheet_by_position(doc, i):
    """Return sheet object by sheet position starting from 0"""
    assert i > 0, 'position must a positive integer'

    name = doc.Sheets.getElementNames()[int(i)]
    return doc.Sheets.getByName(name)

def get_sheet_by_name(doc, name):
    """Return sheet object by sheet name"""
    return doc.Sheets.getByName(name)

def get_active_sheet(doc):
    """Return sheet object"""
    return doc.getCurrentController().getActiveSheet()

def new_sheet(doc, name, i=None):
    """Create new sheet
    
    Args:
        name (str): sheet name
        i (int, optional): position to put new sheet. Starts from 0.
        If None, sheet will be create in the last position. Default is None.
    
    Returns:
        sheet object
    """
    # position
    if i is None:
        i = get_sheet_count(doc) + 1
    assert i >= 0, 'position i must be integer >= 0'

    # unique name
    if name in get_sheet_names(doc):
        msgbox(f'Cannot create sheet with name {name} because there is a sheet with this name already', type_msg='errorbox')
        return
    
    # create sheet
    doc.Sheets.insertNewByName(name, i)

    return get_sheet_by_name(doc, name)

def remove_sheet_by_name(doc, name):
    """delete sheet by name"""
    if name in get_sheet_names(doc):
        if get_sheet_count(doc) > 1:
            doc.Sheets.removeByName(name)
        else:
            msgbox(f"{name} cannot be removed because it is the only existing sheet", type_msg='errorbox')
            return 
    else:
        raise msgbox(f'cannot remove sheet `{name}` because it does not exists', type_msg='errorbox')
        return 

def remove_sheet_by_position(doc, i):
    """delete sheet at position i starting from 0"""
    assert i >= 0, 'position i must be integer >= 0'
    name = doc.Sheets.getElementNames()[int(i)]
    remove_sheet_by_name(doc, name)

def move_sheet(doc, name, i):
    """move sheet to new position
    Args:
        name (str): sheet name
        i (int): new sheet position

    Returns:
        None
    """
    assert i >= 0, 'position i must be integer >= 0'

    if name in get_sheet_names(doc):
        doc.Sheets.moveByName(name, i)
    else:
        raise msgbox(f'cannot move sheet `{name}` because it does not exists', type_msg='errorbox')
        return 

def copy_sheet(doc, name, new_name, i=None):
    """copy sheet"""
    # position
    if i is None:
        i = get_sheet_count(doc) + 1
    assert i >= 0, 'position i must be integer >= 0'

    # check new name
    names = get_sheet_names(doc)
    if new_name in names:
        if name in get_sheet_names(doc):
            msgbox(f'Cannot create sheet with name {new_name} because there is a sheet with this name already', type_msg='errorbox')
            return

    # copy
    if name in names:
        doc.Sheets.copyByName(name, new_name, i)
    else:
        raise msgbox(f'cannot copy sheet `{name}` because it does not exists', type_msg='errorbox')
        return 
# %%

# %% ---------------------- calc styles document functions ----------------- %% #
def get_style_names(doc):
    """return list with available cell/font styles"""
    return [x.Name for x in doc.getStyleFamilies()['CellStyles']]

def remove_style(doc, name):
    """delete cell/font style"""
    if name in get_style_names(doc):
        doc.getStyleFamilies()['CellStyles'].removeByName(name)
    else:
        msgbox(f'Trying to delete style `{name}`, but it does not exist\n\navailable styles: {get_styles(doc)}', type_msg='errorbox')
        return

def new_style(doc, name, properties):
    """Create new style

    Args:
        name (str): style name
        properties (dict): cell/font properties. See below list of useful properties

            'CellBackColor' = -1, 
            'TopBorder'     = {'Color':16776960}
            'CharFontName'  = 'Liberation Sans',
            'CharHeight'    = 20.0,
            'CharWeight'    = 150,
            'CharPosture'   = 2,
            'VertJustify'   = 2,
            'HoriJustify'   = 2,

    Returns:
        None
    """

    # check name
    if name in get_style_names(doc):
        msgbox(f'cannot create style `{name}` because it already exists', type_msg='errorbox')
        return

    # create new style
    new_style = doc.createInstance('com.sun.star.style.CellStyle')
    doc.getStyleFamilies()['CellStyles'].insertByName(name, new_style)

    keys   = tuple(properties.keys())
    values = tuple(properties.values())
    new_style.setPropertyValues(keys, values)
    return

def select_active_sheet(doc, sheet):
    doc.CurrentController.setActiveSheet(sheet)
# %%

# %% -------------------------- sheet functions --------------------------- %% #
def get_cells(sheet, position):
    """return cell or range object from position

    Args:
        sheet (sheet): sheet object
        position (string or tuple): String of type 'A2' or 'A2:J5' or tuple of 
            ints like (col_number, row_number) or (start_col, start_row, 
            stop_col, stop_row). Int starts from 0.

    Returns:
        cell or range object
    """
    if isinstance(position, str):
        position = str2num(position)
    
    if len(position) == 2:
        return sheet.getCellByPosition(*position)
    elif len(position) == 4:
        return sheet.getCellRangeByPosition(position[0], position[1], position[2], position[3])

def get_last_used_row(sheet):
    """Return the bottom most used row in a sheet. Starts from 0"""
    cursor = sheet.createCursor()
    cursor.gotoEndOfUsedArea(False)
    return cursor.getRangeAddress().EndRow

def get_last_used_col(sheet):
    """Return the right most used col in a sheet. Starts from 0"""
    cursor = sheet.createCursor()
    cursor.gotoEndOfUsedArea(False)
    return cursor.getRangeAddress().EndColumn

def delete_rows(sheet, start, stop):
    """delete rows
    
    Args:
        sheet (sheet): sheet object
        start, stop (int): number of the first and last row to delete. Starts from 
            0.
    
    Returns:
        None
    """
    assert isinstance(start, int), f'start must be an integer, not type `{type(start)}`'
    assert isinstance(stop, int), f'stop must be an integer, not type `{type(stop)}`'
    assert start <= stop, f'start must less or the same as stop'

    oRows = sheet.getRows()
    oRows.removeByIndex(start, stop-start+1)

def merge(cells):
    """merge cells"""
    cells.merge(True)

def unmerge(cells):
    """unmerge cells"""
    cells.merge(False)


def get_row(sheet, row):
    return sheet.getCellRangeByPosition(0, row, get_last_used_col(sheet), row)


def get_row_data(sheet, row):
    """return list with row data (row number starting from 0)"""
    return get_row(sheet, row).getDataArray()[0]

def get_col_data(sheet, col):
    """return list with column data (column number starting from 0)"""
    if isinstance(col, str):
        if ':' in col:
            msgbox(f'get_col_data error: input must be a column position, not range\n\ninput = {position}', type_msg='errorbox')
            return
        position = str2num(position)
    return [_[0] for _ in sheet.getCellRangeByPosition(col, 0, col, get_last_used_row(sheet)).getDataArray()]

def set_row_data(sheet, row, value, start_col=0):
    sheet.getCellRangeByPosition(start_col, row, start_col+len(value)-1, row).setDataArray([value, ])
    return 



def get_cells_data(sheet, position):
    """return cell or range values

    Args:
        sheet (sheet): sheet object
        position (string or tuple): String of type 'A2' or 'A2:J5' or tuple of 
            ints like (col_number, row_number) or (start_col, start_row, 
            stop_col, stop_row). Int starts from 0.

    Returns:
        value (for cell) or list of list (for range)
    """
    if isinstance(position, str):
        position = str2num(position)
    if len(position) == 2:
        return sheet.getCellByPosition(*position).getDataArray()[0][0]
    elif len(position) == 4:
        return sheet.getCellRangeByPosition(*position).getDataArray()

def set_cells_data(sheet, position, value):
    """set cell or range values

    Args:
        sheet (sheet): sheet object
        position (string or tuple): String of type 'A2' or 'A2:J5' or tuple of 
            ints like (col_number, row_number) or (start_col, start_row, 
            stop_col, stop_row). Int starts from 0.
        value (value or list of list): if Cell, a unique value. If Ranges, 
            a list of lists (First index is the row, second index is the col)

    Returns:
        None
    """
    if isinstance(position, str):
        position = str2num(position)
    if len(position) == 2:
        sheet.getCellByPosition(*position).setDataArray([[value, ], ])
    elif len(position) == 4:
        sheet.getCellRangeByPosition(*position).setDataArray(value)
# %%

# %% ----------------------- cell/ranges properties ----------------------- %% #
def get_property_names(cell):
    """return property names from cell or ranges
    
    Args:
        cell (cell, range): Cell or range object

    Returns:
        list
    """
    return [x.Name for x in cell.getPropertySetInfo().Properties]

def get_property_value(cell, name):
    """return cell/ranges property value
    
    Args:
        cell (cell, range): Cell or range object
        name (str): property name
        value: new property value
    
    Returns:
        None
    """
    return cell.getPropertyValue(name)

def set_property_value(cell, name, value):
    """set cell/ranges property value
    
    Args:
        cell (cell, range): Cell or range object
        name (str): property name
        value: new property value
    
    Returns:
        None
    """
    return cell.setPropertyValue(name, value)


def get_width(cells):
    """return cells width"""
    temp = get_property_value(cells, 'Size')
    return temp.Width

def get_height(cells):
    """return cells height"""
    temp = get_property_value(cells, 'Size')
    return temp.Height

def set_width(cells, value):
    """use value = 'optimal'"""
    for col in cells.Columns:
        if value == 'optimal':
            col.OptimalWidth = True
        else:
            col.Width = value

def set_height(cells, value):
    """use value = 'optimal'"""
    for col in cells.Rows:
        col.Height = value


def get_border(cell, border):

    if border == 'top':
        border = 'TopBorder'
    elif border == 'bottom':
        border = 'BottomBorder'
    elif border == 'right':
        border = 'RightBorder'
    elif border == 'left':
        border = 'LeftBorder'
    else:
        raise ValueError('not valid border name: top, bottom, right, left')
    
    border = get_property_value(cell, border)

    prop = {'color':          border.Color, 
            'InnerLineWidth': border.InnerLineWidth,
            'OuterLineWidth': border.OuterLineWidth,
            'LineDistance':   border.LineDistance,
            'LineStyle':      border.LineStyle,
            'LineWidth':      border.LineWidth}
    return prop

def set_border(cells, border, color=None, InnerLineWidth=None, OuterLineWidth=None, LineDistance=None, LineStyle=None, LineWidth=None):
    """'color', 'InnerLineWidth', 'OuterLineWidth', 'LineDistance', 
    'LineStyle', 'LineWidth'

    border = top, bottom, right, left, all
    """

    if border == 'top':
        border = 'TopBorder'
    elif border == 'bottom':
        border = 'BottomBorder'
    elif border == 'right':
        border = 'RightBorder'
    elif border == 'left':
        border = 'LeftBorder'
    elif border == 'all':
        pass
    else:
        raise ValueError('not valid border name: top, bottom, right, left')
    
    temp = uno.createUnoStruct("com.sun.star.table.BorderLine2")
    if color is not None:
        temp.color = color
    if InnerLineWidth is not None:
        temp.InnerLineWidth = InnerLineWidth
    if OuterLineWidth is not None:
        temp.OuterLineWidth = OuterLineWidth
    if LineDistance is not None:
        temp.LineDistance = LineDistance
    if LineStyle is not None:
        temp.LineStyle = LineStyle
    if LineWidth is not None:
        temp.LineWidth = LineWidth

    if border == 'all':
        for border in ('TopBorder', 'BottomBorder', 'RightBorder', 'LeftBorder'):
            set_property_value(cells, border, temp)
    else:
        set_property_value(cells, border, temp)


def get_protection(cell):

    prot = get_property_value(cell, 'CellProtection')

    prop = {'IsLocked':        prot.IsLocked, 
            'IsHidden':        prot.IsHidden,
            'IsFormulaHidden': prot.IsFormulaHidden,
            'IsPrintHidden':   prot.IsPrintHidden}
    return prop

def lock_cells(cells):
    """'Lock cells'"""
    temp = uno.createUnoStruct("com.sun.star.util.CellProtection")
    temp.IsLocked = True

    set_property_value(cells, 'CellProtection', temp)

def unlock_cells(cells):
    """'Unlock cells'"""
    temp = uno.createUnoStruct("com.sun.star.util.CellProtection")
    temp.IsLocked = False

    set_property_value(cells, 'CellProtection', temp)
# %%





