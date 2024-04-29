#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import libreoffice_wrapper as lw
import time
import importlib
importlib.reload(lw)

# start LibreOffice and establish communication
pid = lw.start_soffice()
time.sleep(10)
soffice = lw.soffice()

# Open Calc
calc = soffice.Calc()  # tries to connect with any open Calc instance
# If nothing is open, it will start a new spreadsheet
# calc = soffice.Calc('<path-to-spreadsheet-file>')  # connects/opens specific file
# calc = soffice.Calc(force_new=True)  # open a new file

# Calc info
print(calc.get_title())
print(calc.get_filepath())
print(calc.get_sheets_count())
print(calc.get_sheets_name())

# save
calc.save()
# calc.save('<path-to-save>')

# close Calc
# calc.close()

# insert new sheet
calc.insert_sheet('my_new_sheet')
calc.insert_sheet('sheet_to_be_remove')
calc.insert_sheet('another_sheet_to_be_remove')

# remove sheet
calc.remove_sheets_by_position(3)
calc.remove_sheet('sheet_to_be_remove')

# move sheet
calc.move_sheet(name='my_new_sheet', position=0)

# copy_sheet
calc.copy_sheet(name='my_new_sheet', new_name='copied_sheet', position=2)

# sheet name and position
print(calc.get_sheet_position(name='my_new_sheet'))
print(calc.get_sheet_name_by_position(position=0))

# Styles
print(calc.get_styles())
properties = {'CellBackColor':16776960, 'CharWeight':150}
calc.new_style(name='my_new_style', properties=properties)
calc.remove_style(name='my_new_style')

# get sheet
sheet = calc.get_sheet_by_position(0)
sheet = calc.get_sheet('my_new_sheet')

# sheet name
print(sheet.get_name())
sheet.set_name('new_name')

# visibility
print(sheet.isVisible())

# move
sheet.move(position=2)  # in this case moving to 0 or 1 yields the same result

# remove (delete)
# sheet.remove()

# set/get data (data can be set in many ways)
sheet.set_value('A1', 'name')
print(sheet.get_value('A1'))

sheet.set_value('B', '1', 'color')
print(sheet.get_value('B', '1'))

sheet.set_value('C', 0, 'quantity')
print(sheet.get_value('C', 0))

sheet.set_value(3, 0, 'taste')
print(sheet.get_value(3, 0))

sheet.set_value(4, '1', 'weight')
print(sheet.get_value(4, '1'))

sheet.set_value('A2:C3', [['apple', 'red', 3], ['banana', 'yellow', 6]])
print(sheet.get_value('A2:C3'))

sheet.set_value('A4', 'C5', [['orange', 'orange', 4], ['pineapple', 'yellow', 1]])
print(sheet.get_value('A4', 'C5'))

sheet.set_value('A6', [['grapes', 'purple', 12], ['zuchini', 'green', 6]])
print(sheet.get_value('A6:C7'))

sheet.set_value('A', '8', [['avocado', 'green', 1], ['pear', 'yellow', 10]])
print(sheet.get_value('A8:C9'))

sheet.set_value('A', 9, [['lettuce', 'green', 21], ['watermelon', 'green', 2]])
print(sheet.get_value(0, 9, 2, 10))

sheet.set_value('A', '12', 'C', '13', [['potato', 'yellow', 10], ['carrot', 'orange', 3]])
print(sheet.get_value('A', '12', 'C', '13'))

sheet.set_value(0, 13, 2, 14, [['spinach', 'green', 4], ['lemon', 'green', 2]])
print(sheet.get_value(0, 13, 2, 14))

# If necessary, the cell format can be set to 'formula', 'string', or 'number'
# format = number --> forces values to set as number
# format = string --> forces values to be set as string (text)
# format = formula --> works fine for strings, numbers, and formulas
sheet.set_value('E2', 10, format='string')
sheet.set_value('E3', 10, format='formula')
sheet.set_value('E4', value=20, format='number')

sheet.set_value('E5', '20', format='string')
sheet.set_value('E6', '20', format='formula')
sheet.set_value('E7', '20', format='number')

sheet.set_value('E8', '=E4', format='string')
sheet.set_value('E9', '=E4', format='formula')
# sheet.set_value('D10', '=E4', format='number')  # will raise an error

sheet.set_value('E11', '10/05/2021', format='string')
sheet.set_value('E12', '10/05/2021', format='formula')
# sheet.set_value('E13', '10/05/2021', format='number')  # will raise an error

# default for set is 'formula' which should work fine in most cases
# default for get is 'string' which should work fine in most cases


# set values of entire rows/column
# it clears the row/column before seting new values
sheet.set_row(15, value=['mango', 'red', 3, 'sweet', 10.9])
sheet.set_row('17', value=['papaya', 'yellow', 1, 'sweet', 12.0])
sheet.set_row('B2', value=['red', 2, 'sweet', 40.1])
sheet.set_row('4', column_start='C', value=[6, 'acid', 5.12])

print(sheet.get_row(12))
print(sheet.get_row('12'))
print(sheet.get_row('B16'))
print(sheet.get_row('16', column_start='B'))
print(sheet.get_row('16', column_start=1))

sheet.set_column('E', row_start='2', value=[1.27, 2.23, 1.50, 6.5])
sheet.set_column('E', row_start=5, value=[5.27, 1.28, 2.50, 6.12])
sheet.set_column(4, row_start=9, value=[5.00, 0.28, 9.10, 1.02])
sheet.set_column('E14', value=[5.1, 0.53, 9.11, 1.10])

print(sheet.get_column(0))
print(sheet.get_column('A'))
print(sheet.get_column('A1'))
print(sheet.get_column('B', row_start='1'))
print(sheet.get_column('B', row_start=0))

# last row/column (spreadsheet size)
print(sheet.get_last_row())
print(sheet.get_last_column())

# length of row/column
print(sheet.get_row_length(11))
print(sheet.get_row_length('11'))
print(sheet.get_column_length('B'))
print(sheet.get_column_length('B'))
print(sheet.get_column_length(1))

# column width
print(sheet.get_column_width(2))
print(sheet.get_column_width('C'))
sheet.set_column_width('C', 1500)
sheet.set_column_width(['A', 'B'], [2000, 1500])
sheet.set_column_width([3, 4], 2000)

# row height
print(sheet.get_row_height(0))
print(sheet.get_row_height('1'))
sheet.set_row_height('1', 1000)
sheet.set_row_height(['2', '3'], [1000, 1000])
sheet.set_row_height([1, 2], 452)

# clear cells
# sheet.clear()
# sheet.clear('F1')
# sheet.clear('F1:G1')
# sheet.clear('F2', 'G2')
# sheet.clear('F', '3', 'G', '3')
# sheet.clear(5, 3, 6, 3)
# sheet.clear_row(13)
# sheet.clear_row('B13')
# sheet.clear_row(13, column_start='B')
# sheet.clear_column('B')
# sheet.clear_column('B2')
# sheet.clear_column('B', row_start='2')

# merge
sheet.merge('F1:G1')
sheet.merge('F2', 'G2')
sheet.merge('F', '3', 'G', '3')
sheet.merge(5, 3, 6, 3)

sheet.unmerge('F1:G1')
sheet.unmerge('F2', 'G2')
sheet.unmerge('F', '3', 'G', '3')
sheet.unmerge(5, 3, 6, 3)

# cell properties
print(sheet.cell_properties())
print(sheet.cell_properties('B2'))

print(sheet.get_property(0, 0, 'CellBackColor'))
print(sheet.get_property('A1', 'CellBackColor'))
sheet.set_property(0, 0, 'CellBackColor', int('9c9c9c', 16))  # color must be int
sheet.set_property('A1:E1', 'CellBackColor', int('9c9c9c', 16))  # color must be int

print(sheet.get_property(0, 0, 4, 0, 'CharWeight'))
sheet.set_property(0, 0, 4, 0, 'CharWeight', 150)

sheet.set_property('A1:E1', 'VertJustify', 2)
sheet.set_property('A1:E1', 'HoriJustify', 2)

sheet.set_property('E1:E100', 'HoriJustify', 3)

sheet.get_property('B2', 'TopBorder')
sheet.set_property('B2', 'TopBorder.LineWidth', 100)

d = sheet.get_property('A', 1, 'TopBorder')
d['LineWidth'] = 100
sheet.set_property('A', 1, 'TopBorder', d)

sheet.set_property('C2:E2', 'TopBorder.LineWidth', 100)

# conditional format
print(sheet.get_conditional_formats())
sheet.new_conditional_format('E2:E100', Operator='is greater than', Formula1='1', StyleName='Good')
sheet.new_conditional_format('C2:C100', Operator='is equal to', Formula1='1', StyleName='Bad')

print(sheet.get_conditional_formats())
sheet.remove_conditional_format(range_index=1)

# saving
calc.save()

# close
calc.close()

# close communication
soffice.kill()
