#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""This example shows how to use brixs.sheets to get read/write excel files

As for today (2024), excel does not have a dedicated python API, therefore 
we have to use openpyxl, which is a Python library to read/write Excel.

Note that libreoffice HAS a dedicated python API, see brixs/sheets/ods.py
"""

# %% =========================== brixs imports =========================== %% #
from brixs.sheets.xlsx import xlsx, new_xlsx
# %%

# % ====================================================================== %% #
# % =========================== creating a file ========================== %% #
# %% ===================================================================== %% #

# create new xlsx file or open an existing file
sheet = new_xlsx('test.xlsx', sheetname='brixs')  # new file
# or open an existing file
# sheet = xlsx(filepath) 

# setting data
sheet.set_row([1, 2, 3, 4],    row=2, start_col=4)
sheet.set_row([5, 6, 7, 8],    row=3, start_col=4)
sheet.set_row([9, 10, 11, 12], row=4, start_col=4)

# setting fill and font
sheet.set_row_fill(color='ff0000', row=2, start_col=4, stop_col=8)
sheet.set_row_font(row=2, bold=True, color='000000', start_col=4, stop_col=8)

# setting and getting column width
w = sheet.get_col_width(col=5)
sheet.set_col_width(col=5, width=w*1.5)

# saving file
sheet.save()
# %%

# % ====================================================================== %% #
# % ============================ loading a file ========================== %% #
# %% ===================================================================== %% #
# open a file
sheet2 = xlsx('test.xlsx', sheetname='brixs') 

# get data
data = sheet2.get_col_data(4, start_row=2)
print(data)
data = sheet2.get_col_data('E', start_row=2)
print(data)
data = sheet2.get_row_data(3, start_col=4)
print(data)
# %%


# % ====================================================================== %% #
# % ================================= other ============================== %% #
# %% ===================================================================== %% #
# sheet2 is an openpyxl Worksheet object and comes with all its functionality
# refer to openpyxl for documentation

# open a file
sheet2 = xlsx('test.xlsx', sheetname='brixs') 

# get data from cell
sheet2['D4'].value

# get data from row and column
data = [_[0].value for _ in sheet2['D2:D6']]  # column
data = [_.value for _ in sheet2['D2:F2'][0]]  # row

# list of methods and attrs that can be used
print([_ for _ in sheet2.__dir__() if _.startswith('_')==False])
# %%


