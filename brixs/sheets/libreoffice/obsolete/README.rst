===================
libreoffice-wrapper
===================

Python module for controlling `LibreOffice`_ programs (Writer, Calc, Impress, Draw, Math, and Base). Currently, manipulation of Calc instances is (somewhat) fully implemented and the module supports features such as:

[x] core functionality (open, save, close, ...)

[x] save multiple formats (ods, xlsx, ...)

[x] add/remove styles

[x] insert/delete/move sheets

[x] get/set values from a cell/range/rows/columns

[x] get/set cell/range properties (color, border, ...)

[x] merge/unmerge cells

[x] conditional formatting

 Manipulation of Writer, Impress, Draw, and Math instances is in its early development and the module only allows for basic core functionality such as opening/closing/saving files. Base is not implemented at all and trying to open a LibreOffice Base instance will raise an error.

About
==========

This module uses `tmux`_ to intermediate communication between a Python instance and the LibreOffice's internal python interpreter, which has free access to LibreOffice's Python API that allows controlling the LibreOffice components. This way, you are not limited to the functionality of LibreOffice's internal python, i.e., one is able to manipulate LibreOffice components from any Python terminal and inside Python environments. In addition to that, modifications to a file happen "real time" (no need to reload the file).

It was tested on:

- Linux (Ubuntu 20.04) and LibreOffice Version 7.0.4.2

and it should also work fine on MacOS. Currently, I don't think it will work on Windows since `tmux` is not implement there. However, it might work if one uses Windows Subsystem for Linux (WSL). I'm still trying to make it work.

 DISCLAIMER: At first, this module was built to allow manipulating of Calc spreadsheets without the need to reload the document every time a modification was made. Since, this is done, I'm not sure I will keep working on it in order to extend the functionality to Writer, Impress, etc.. In any case, it should be easy enough to implement code for those since the core functionality is the same.


Dependencies
=============

This module is heavily dependent on `tmux`_ which can be installed via apt-get on Debian or Ubuntu (check `tmux`_ page for installation instruction on other OS)::

  apt install tmux

The package `libtmux`_ (tmux workspace manager in python) is also necessary::

  pip install libtmux

Lastly, some features for dealing with Calc instances uses `numpy`_::

  pip install numpy


Usage (Initialization)
=======================

Firstly, one has to start the office in Listening Mode. This can be done by opening the terminal and issuing the command::

  soffice -accept=socket,host=0,port=8100;urp;

Alternatively, libreoffice-wrapper has a built-in function that starts LibreOffice in Listening Mode,

.. code-block:: python

    import libreoffice_wrapper as lw

    pid = lw.start_soffice()


.. The function :python:`lw.start_soffice()` returns the pid of the process. Note that, this function starts a ``tmux`` session called ``libreoffice-wrapper`` with a window named ``soffice``, which can be accessed on a different terminal via ``tmux``. In addition to that, ```lw.start_soffice()``` searches for LibreOffice in the default folder ``/opt/libreoffice7.0``. If LibreOffice is installed in a different folder, it must be passed as an argument of the function ```lw.start_soffice(folder=<path-to-libreoffice>)```.

Once LibreOffice is up and running in listening mode, one can now establish the communication,

.. code-block:: python

  soffice = lw.soffice()

.. where `lw.soffice()` starts a `tmux` session `'libreoffice-wrapper'` with a window named `'python'`, with opens the internal LibreOffice's Python interpreter. After that, the `soffice` object manages to communicate to LibreOffice through this Python instance opened in this `tmux` window.

To close LibreOffice and its tmux session,

.. code-block:: python

  soffice.kill()

.. which just ends the `tmux` session.

Calc
========

Most functionality regarding Calc spreadsheets can be found in the example below,

.. code-block:: python

  import libreoffice_wrapper as lw
  import time

  # start LibreOffice and establish communication
  pid = lw.start_soffice()
  time.sleep(5)
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




Writer, Impress, Draw, Math and Base
======================================

Manipulation of Writer, Impress, Draw, and Math instances is in its early development and the module only allows for basic core functionality such as opening/closing/saving files. Base is not implemented at all and trying to open a LibreOffice Base instance will raise an error.

.. code-block:: python

  import libreoffice_wrapper as lw

  # start LibreOffice and establish communication
  pid = lw.start_soffice()
  soffice = lw.soffice()

  # Writer
  writer = soffice.Writer()
  writer.save()
  writer.close()

  # Impress
  impress = soffice.Impress()
  impress.save()
  impress.close()

  # Draw
  draw = soffice.Draw()
  draw.save()
  draw.close()

  # Math
  math = soffice.Math()
  math.save()
  math.close()

  # close LibreOffice/communication
  soffice.kill()



.. _tmux: https://github.com/tmux/tmux/wiki
.. _LibreOffice: https://www.libreoffice.org/
.. _libtmux: https://github.com/tmux-python/libtmux
.. _numpy: https://numpy.org/
