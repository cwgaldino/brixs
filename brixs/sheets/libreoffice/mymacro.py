"""Macro for creating an online scanlist for experiments

See below instruction on how to use/install it

Items 1 -- 5   will talk about how to install this macro
Items 6 -- 8  will talk about how to use this macro
Item  9       will talk about how to write a new macro or edit this one
 

#######################################
# 1. Download and install libreoffice #
#######################################

Follow the instructions on the libreoffice webpage, download and install it.

############################################
# 2a. Install pip on libreoffice (WINDOWS) #
############################################

Option 1: via zaz-pip libreoffice extension

1. Go to zaz-pip git page here -> https://git.cuates.net/elmau/zaz-pip/src/branch/master/source
2. Follow the instructions on the git page for installation
3. Follow instructions on the git page for installing packages

Option 2: via pypi bootstrap

1. Download the pypi bootstrap here -> https://bootstrap.pypa.io/get-pip.py
2. Start a command prompt at the LibreOfice installation directory (e.g., cd C:\Program Files\LibreOffice\program)
3. Run .\python.exe <path-to-get-pip>\get-pip.py
4. Use command .\python.exe -m pip install <package-name> for installing packages

For me, regardless of the option, packages were being saved here ->
C://Users//galdin_c//Appdata//Roaming//Python//Python38//site-packages

##########################################
# 2b. Install pip on libreoffice (LINUX) #
##########################################

To do. Must be similar as on Windows.

################################
# 3. Install required packages #
################################

The required packages for brixs are: numpy, matplotlib, and h5py (h5py is not
required for brixs, but will br necessary to use the file reading functionality).
These packages must be installed via pip (see step 2 for instructions on how to 
install pip and how to install packages).

####################
# 4. Install brixs #
####################

brixs can be imported via adding brixs path to sys.path in this file like this:

import sys
sys.path.append(r'<path-to-brixs>')
import brixs as br

or it can be imported via pip by pointing to the newest distribution .tar.gz or 
.whl, for example

pip install <path-to-dist-folder>/'brixs-0.9.6-py3-none-any.whl'

However, I think installing via adding brixs to the path is safer for now because
If changes are made to brixs, then you don't need to install it again is the 
libreoffice python environment.

#######################################
# 5. Install the macro on libreoffice #
#######################################

Place this file in the libreoffice python share folder. On my computer this folder
is located here:

C:\Program Files\LibreOffice\share\Scripts\python\

From now on, every function listed on `g_exportedScripts` variable inside this
script will be accessible as a macro

g_exportedScripts = (function1, function2, ...)

###############################
# 6. run macro on libreoffice #
###############################

Tools > Macros > Run Macro > library:`Application Macros` > <name-of-this-file> > <macro-name> > Run


##########################################
# 7. add a Macro shortcut on libreoffice #
##########################################

View > Toolbars > Customize > Macros > Application Macros > select the desired Macro >
> press the right arrow

Down below there is a panel called `Customize` in which you can press the button 
`insert` to insert a separator. You can also press the button `Modify` to change
the name of the macro or put an icon. 

############################
# 8. How to use this macro #
############################

This macro expects that the Calc file has two Sheets called 'settings' and 
'scanlist'.

The `update` macro will read the files in a folder defined in the 'settings' 
sheet and return the metadata. The files inside that folder must be ordered. The
first row of 'scanlist' is user defined names and can be anything. The second 
row must have names that match exactly the target metadata.

The `metadata` macro down below will open a msgbox with all the available 
metadata of a filepath.

##########################################
# 9. How to modify/adapt this macro file #
##########################################

This file is supposed to serve as a reference. You won't need all the functions
defined here. Keep only what you need.

For sure, you cannot remove the `MACROS` section and the functions that are used 
in the `MACROS` section.

For sure, you will have to edit the section `READ` in which you define the 
function to read the files. The only requirement here is that the read() 
function returns a object with attrs defined. Doesn't even have to be a brixs 
type.

"""
# %% ------------------------- standard imports --------------------------- %% #
from pathlib import Path
import sys

# %% ---------------------------- uno imports ----------------------------- %% #
import uno
from com.sun.star.awt import MessageBoxButtons as MSG_BUTTONS
from com.sun.star.beans import PropertyValue

CTX = uno.getComponentContext()
SM  = CTX.getServiceManager()

# %% -------------------------- brixs imports ----------------------------- %% #
path2brixs = Path(r'C:\Users\galdin_c\github\brixs')
sys.path.append(str(path2brixs))
import brixs as br
import brixs.beamlines.VERITAS as VERITAS
import datetime

# %% -------------------------- support functions ------------------------- %% #
def letter2num(letter):
    """Returns position of letter in the alphabet starting from 0 (A, B, ..., Z, AA, AB, ..., ZZ, AAA, ...)"""
    letter   = letter.lower()
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    n = 0
    for idx, s in enumerate(letter):
        n += alphabet.index(s)+(idx)*26
    return n

import re
def str2num(string):
    """Returns col, row number based on string like 'A2', 'AJ20', 'A2:AJ20', ...

    Returns:
        col, row for `cell` strings ('A2', 'AJ20', ...)
        col_start, row_start, col_stop, row_stop for `range` strings ('A2:AJ20', ...) 
        col and row numbering starts from 0
    """
    # check input
    if isinstance(string, str) == False:
        msgbox(f'_cell2num error: input must be type str\n\ninput = {string}\ninput type={type(string)}')
        return
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
        msgbox(f'_cell2num error: cannot covert input to number\n\ninput = {string}')
        return
# %%

# %% --------------------------- base functions --------------------------- %% #
def create_instance(name, with_context=False):
    if with_context:
        instance = SM.createInstanceWithContext(name, CTX)
    else:
        instance = SM.createInstance(name)
    return instance    
# %%

# %% -------------------------- message box ------------------------------- %% #
def msgbox(message, title='LibreOffice', buttons=MSG_BUTTONS.BUTTONS_OK, type_msg='infobox'):
    """Create message box

    Args:
        message (string): text
        title (string, optional): msgbox title
        buttons ():
        type_msg (string): type of boxmsg. Default is `infobox`. Options are:

            infobox
            warningbox
            errorbox
            querybox
            messbox

    Retuns:
        None

    See also:
        https://api.libreoffice.org/docs/idl/ref/interfacecom_1_1sun_1_1star_1_1awt_1_1XMessageBoxFactory.html
    """
    toolkit = create_instance('com.sun.star.awt.Toolkit')
    parent  = toolkit.getDesktopWindow()
    mb      = toolkit.createMessageBox(parent, type_msg, buttons, title, str(message))
    return mb.execute()
# %%

# %% -------------------------- new document functions -------------------- %% #
def get_current_document():
    """Return document object"""
    return XSCRIPTCONTEXT.getDocument()

def new_calc():
    """open calc and return document object"""
    URL = 'private:factory/scalc'
    desktop = create_instance('com.sun.star.frame.Desktop')
    return desktop.loadComponentFromURL(URL, '_default', 0, ())

def new_writer():
    """open writer and return document object"""
    URL = 'private:factory/swriter'
    desktop = create_instance('com.sun.star.frame.Desktop')
    return desktop.loadComponentFromURL(URL, '_default', 0, ())

def new_impress():
    """open impress and return document object"""
    URL = 'private:factory/SdXImpressDocument'
    desktop = create_instance('com.sun.star.frame.Desktop')
    return desktop.loadComponentFromURL(URL, '_default', 0, ())

def new_draw():
    """open draw and return document object"""
    URL = 'private:factory/sdraw'
    desktop = create_instance('com.sun.star.frame.Desktop')
    return desktop.loadComponentFromURL(URL, '_default', 0, ())

def new_math():
    """open math and return document object"""
    URL = 'private:factory/smath'
    desktop = create_instance('com.sun.star.frame.Desktop')
    return desktop.loadComponentFromURL(URL, '_default', 0, ())

def new_base():
    """open base and return document object"""
    raise NotImplementedError('base not implemented yet')

def open_document(filepath):
    """open libreoffice file and return document object"""
    filepath = Path(filepath)
    URL      = uno.systemPathToFileUrl(str(filepath))

    desktop = create_instance('com.sun.star.frame.Desktop')
    return desktop.loadComponentFromURL(URL, "_default", 0, ())
# %%

# %% ------------------------ base document functions --------------------- %% #
def get_filepath(doc):
    """return doc filepath or empty string '' if no filepath has been linked to file yet"""
    if doc.hasLocation():
        return uno.fileUrlToSystemPath(doc.getURL())
    else:
        return ''

def get_filename(doc):
    """return document filename"""
    return doc.getTitle()

def save(doc, filepath=None, ext=None, verbose=True):
    """Save document (only set for calc)

    Args:
        filepath (str or Path, optional): If None, it will try to get current
            filepath associated with the document. It raises an error/msgbox if 
            document has no filepath associated. User have to save it by hand
            at least once.
        ext (str, optional): file extension. If None, it will try to guess ext from
            filename. Default is None. Available options: 

            calc8/.ods:    ['ods', '.ods', 'calc8', 'calc', 'Calc']
            excel/.xlsx:   ['excel', 'xlsx', '.xlsx']
            Excel 97/.xls: ['xls', '.xls']
            csv/.csv:      ['csv', '.csv', 'text', 'txt', '.txt']

    Returns:
        None
    """
    # filepath
    if filepath is None:
        filepath = get_filepath(doc)
        if filepath == '' or filepath == 'Untitled 1':
            msgbox('filepath not defined. Cannot save file via macro. Please, save document by hand at least once first.', type_msg='errorbox')
            return
        
    # if ext is None, try and get ext from filename
    if ext is None:
        ext = get_filename(doc).split('.')[-1]
        
    # check ext
    if ext in ['ods', '.ods', 'calc8', 'calc', 'Calc']:
        ext = 'calc8'
        extension = '.ods'
    elif ext in ['excel', 'xlsx', '.xlsx']:
        ext = 'Calc MS Excel 2007 XML'
        extension = '.xlsx'
    elif ext in ['xls', '.xls']:
        ext = 'MS Excel 97'
        extension = '.xls'
    elif ext in ['csv', '.csv', 'text', 'txt', '.txt']:
        ext = 'Text - txt - csv (StarCalc)'
        extension = '.csv'
    else:
        msgbox('Error: extension `{ext}` for saving file not recognized', type_msg='errorbox')
        return

    # filepath with extension
    filepath = Path(filepath).with_suffix(extension)
    
    # save document
    URL        = uno.systemPathToFileUrl(str(filepath))
    properties = (PropertyValue('FilterName', 0, ext, 0),)
    doc.storeAsURL(URL, properties)

    if verbose:
        msgbox(f'Saved via macro: {filepath}')

def close(doc):
    """Close document"""
    doc.close(True)

def lock_undo(doc):
    """after locking edits cannot be traced back"""
    doc.getUndoManager().lock()

def unlock_undo(doc):
    """after unlocking edits will be traced"""
    doc.getUndoManager().unlock()
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


def get_current_selection():
    """return a cell or range object"""
    desktop = XSCRIPTCONTEXT.getDesktop()
    return desktop.CurrentComponent.CurrentController.getSelection()
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


def get_size(cells):
    temp = get_property_value(cells, 'Size')
    return temp.Width, temp.Height

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

# %% ----------------------------- Colors --------------------------------- %% #
light_grey_2   = 11711154
light_gold_4   = 16774606
light_indigo_4 = 14605542
gold           = 16750383
# %%

# ============================================================================ #
# ============================= READ (EDIT THIS) ============================= #
# %% ====================================================================== %% #
# the only requirement here is that the read() function returns a object with 
# attrs defined. Doesn't even have to be a brixs type
def read_I21(filepath):

    # read file
    try:
        im, ims = I21.read(filepath)
        im.type = 'rixs'
    except KeyError:
        im, b, c, d = I21.read_xas(filepath)
        im.type = 'not rixs'

    # fix attrs
    im.scan = int(filepath.name.split('-')[-1].split('.')[0])
    for attr in ('start_time', 'end_time', 'checkbeam', 'exposure_time', 'finished'):
        try:
            setattr(im, attr, str(getattr(im, attr)))
        except:
            pass

    return im

def read_spectrum(filepath):
    return br.Spectrum(filepath)

# assign the correct read function to the name `read`
read = read_I21
# %%




# ============================================================================ #
# ======================== MACRO - VERITAS - 2024-04-29 ====================== #
# %% ====================================================================== %% #

##################
# new empty file #
##################
header_rixs = {}
header_rixs['basic'] = ['scan', 'start_time', 'end_time', 'exposure_time', 'exposure_time_calc', 'E', 'pol', 'T', 'th', 'sample_x', 'sample_y','sample_z']
header_rixs['extra'] = ['a_mp1_x', 'a_mp1_y', 'a_mp1_z', 'a_mp1_yaw', 'temperature', 'epu_r3_316_gap', 'epu_r3_316_phase', 'mono_energy_calib', 'veritasarm_energy']
header_rixs['startmetadata'] = ['startmetadata_Sample', 'startmetadata_a10_detmc1_pit', 'startmetadata_a10_detmc1_y', 'startmetadata_a10_pol1_y', 'startmetadata_a7_gr_baff_h_gap', 'startmetadata_a7_gr_baff_h_offset', 'startmetadata_a7_gr_baff_hl', 'startmetadata_a7_gr_baff_hr', 'startmetadata_a7_gr_baff_v_gap', 'startmetadata_a7_gr_baff_v_offset', 'startmetadata_a7_gr_baff_vb', 'startmetadata_a7_gr_baff_vt', 'startmetadata_a9_btl1_y', 'startmetadata_a9_btl2_y', 'startmetadata_a9_tab_z', 'startmetadata_a_mp1_x', 'startmetadata_a_mp1_y', 'startmetadata_a_mp1_yaw', 'startmetadata_a_mp1_z', 'startmetadata_a_slit1_v', 'startmetadata_temperature', 'startmetadata_beamline_energy', 'startmetadata_epu_r3_316_gap', 'startmetadata_epu_r3_316_phase', 'startmetadata_mono_energy', 'startmetadata_mono_energy_calib', 'startmetadata_q_angle', 'startmetadata_rixsdet_theta', 'startmetadata_rixsdet_y', 'startmetadata_rixsdet_z']
header_rixs['endmetadata']   = ['endmetadata_Sample', 'endmetadata_a10_detmc1_pit', 'endmetadata_a10_detmc1_y', 'endmetadata_a10_pol1_y', 'endmetadata_a7_gr_baff_h_gap', 'endmetadata_a7_gr_baff_h_offset', 'endmetadata_a7_gr_baff_hl', 'endmetadata_a7_gr_baff_hr', 'endmetadata_a7_gr_baff_v_gap', 'endmetadata_a7_gr_baff_v_offset', 'endmetadata_a7_gr_baff_vb', 'endmetadata_a7_gr_baff_vt', 'endmetadata_a9_btl1_y', 'endmetadata_a9_btl2_y', 'endmetadata_a9_tab_z', 'endmetadata_a_mp1_x', 'endmetadata_a_mp1_y', 'endmetadata_a_mp1_yaw', 'endmetadata_a_mp1_z', 'endmetadata_a_slit1_v', 'endmetadata_temperature', 'endmetadata_beamline_energy', 'endmetadata_epu_r3_316_gap', 'endmetadata_epu_r3_316_phase', 'endmetadata_mono_energy', 'endmetadata_mono_energy_calib', 'endmetadata_q_angle', 'endmetadata_rixsdet_theta', 'endmetadata_rixsdet_y', 'endmetadata_rixsdet_z']

header_xas = {}
header_xas['basic']    = ['scan', 'scan_type', 'start_time', 'end_time', 'command', 'scanned_motors'] + ['E', 'Pol', 'T', 'exit_slit'] + ['sample_x', 'sample_y', 'sample_z', 'th']
header_xas['pre_scan'] = ['pre_scan_a_m4_lateral',  'pre_scan_a_m4_pitch', 'pre_scan_a_m4_roll', 'pre_scan_a_m4_vertical', 'pre_scan_a_m4_yaw', 'pre_scan_a_mp1_x', 'pre_scan_a_mp1_y', 'pre_scan_a_mp1_yaw', 'pre_scan_a_mp1_z', 'pre_scan_a_slit1_hz', 'pre_scan_m1_yaw', 'pre_scan_m3_yaw', 'macro_start_time', 'entry_identifier']

def _input(sheet, text, text_position, input_position, input_default='', editable=True, text_width='optimal', input_width=False):
    """quickly define input cells"""
    # help text
    if isinstance(text_position, str):
        text_position = str2num(text_position)
    cells = get_cells(sheet, text_position)
    if len(text_position) == 4:
        merge(cells)
    set_cells_data(sheet, text_position[:2], text)
    if text_width:
        set_width(cells, text_width)

    # input
    if isinstance(input_position, str):
        input_position = str2num(input_position)
    cells = get_cells(sheet, input_position)
    if len(input_position) == 4: 
        merge(cells)
    set_cells_data(sheet, input_position[:2], input_default)  # default input
    set_property_value(cells, 'CellBackColor', -1)            # set color as null
    if editable:
        unlock_cells(cells)                                   # unlock cells
    if input_width:
        set_width(cells, input_width)

def _output(sheet, text, text_position, out_position, out_default='', editable=False, text_width='optimal', out_width=False):
    """quickly define output cells"""
    # help text
    if isinstance(text_position, str):
        text_position = str2num(text_position)
    cells = get_cells(sheet, text_position)
    if len(text_position) == 4:
        merge(cells)
    set_cells_data(sheet, text_position[:2], text)
    set_property_value(cells, 'CellBackColor', gold)
    if text_width:
        set_width(cells, text_width)

    # input
    if isinstance(out_position, str):
        out_position = str2num(out_position)
    cells = get_cells(sheet, out_position)
    if len(out_position) == 4: 
        merge(cells)
    set_cells_data(sheet, out_position[:2], out_default)      # default input
    set_property_value(cells, 'CellBackColor', light_gold_4)  # set color
    if editable:
        unlock_cells(cells)                                   # unlock cells
    if out_width:
        set_width(cells, out_width)

def _instruct(sheet, text, position):
    if isinstance(position, str):
        position = str2num(position)

    cells = get_cells(sheet, position)

    if len(position) == 4:
        merge(cells)
    set_property_value(cells, 'CellBackColor', light_gold_4)
    set_property_value(cells, 'IsTextWrapped', True)
    set_property_value(cells, 'VertJustify', 1)
    set_cells_data(sheet, position[:2], text)

def new_empty(*args, **kwargs):
    """Create new template spreadsheet"""
    # open new spreadsheet (Calc)
    doc = new_calc()

    # lock undo
    lock_undo(doc)

    # add new sheets
    settings = new_sheet(doc, 'settings', i=0)
    rixs     = new_sheet(doc, 'RIXS',     i=1)
    xas      = new_sheet(doc, 'XAS',      i=2)
    remove_sheet_by_name(doc, 'Sheet1')
    select_active_sheet(doc, settings)

    # unprotect
    settings.unprotect('123456')

    ##############
    # BACKGROUND #
    ##############
    cells = get_cells(settings, 'A1:F6')
    set_property_value(cells, 'CellBackColor', light_grey_2)

    ##########
    # WIDTHS #
    ##########
    cells = get_cells(settings, 'G1:H1')
    set_width(cells, 1000)

    ################
    # INSTRUCTIONS #
    ################
    if True:
        text  = 'This sheet is password protected to avoid mistakes\n\n'
        text += 'The password is: 123456\n\n'
        text += 'To unlock it:\nTools > Protect Sheet > Type in password'
        _instruct(settings, text, 'I8:J18')

        text  = 'INSTRUCTIONS:\n\n'
        text += '1) Edit setting:\n\n'
        text += '    RIXS filepath:\n        path to h5 file with RIXS scans\n\n'
        text += '    XAS filepath:\n        path to h5 file with XAS scans\n\n'
        text += "    Save after update:\n        if `yes`, spreadsheet will saved after running the `update` macro\n\n" 
        text += "    Intercalate bkg color:\n        if `yes`, lines will have alternating bkg color\n\n" 
        text += '\n2) Run update macro:\n'
        text += '    Tools > Macros > Run Macro > `Application Macros`:`mymacro`:`update` > Run'
        _instruct(settings, text, 'A8:F26')

        text  = 'REMINDERS:\n\n'
        text += '    - Add `update` macro to toolbar:\n        Right click on the toolbar > Customize Toolbar... > Search: `update`\n\n'
        text += '    - Turn grid for colored cells:\n        Tools > Options > LibreOffice Calc > View > Visual Aids > Grid Lines > select: `Show on colored cells`\n\n'
        text += "    - Freeze First rows/cols:\n        View > Freeze Cells > Freeze First column/Row\n\n"
        text += "    - Freeze multiple rows/cols:\n        Select the cell which left/above that cell should freeze > View > Freeze Rows and Columns\n\n"
        text += "    - Move columns:\n        Select a column by clicking the column header > hold down Alt and drag the column by one of its cells (it you try to drag from the header it will not work)\n\n"
        _instruct(settings, text, 'A28:F45')

    ##########
    # Inputs #
    ##########
    if True:
        _input(settings, 'RIXS filepath',         'A1', 'B1:F1', '', True, input_width=4000)
        _input(settings, 'XAS filepath',          'A2', 'B2:F2', '', True)
        _input(settings, 'Save after update',     'A4', 'B4', 'yes', True)
        _input(settings, 'Intercalate bkg color', 'A5', 'B5', 'yes', True)

    ############
    # Internal #
    ############
    if True:
        _output(settings, 'RIXS: number of scans', 'I1', 'J1', 0)
        _output(settings, 'RIXS: next row',        'I2', 'J2', 2)

        _output(settings, 'XAS: number of scans', 'I3', 'J3', 0)
        _output(settings, 'XAS: next row',        'I4', 'J4', 2)

    ###########
    # protect #
    ###########
    settings.protect('123456')

    ##########
    # HEADER #
    ##########
    if True:
        # rixs
        final = ['dataset'] + header_rixs['basic'] + header_rixs['extra'] + ['error', 'Comments']
        set_row_data(rixs, 0, final)
        set_row_data(rixs, 1, final[1:-1], start_col=1)

        cells = get_cells(rixs, position=(0, 0, len(final)+2, 0))
        set_property_value(cells, 'CharWeight', 200)
        set_property_value(cells, 'CharHeight', 12)
        set_property_value(cells, 'CellBackColor', light_grey_2)

        cells = get_cells(rixs, position=(0, 1, len(final)+2, 1))
        set_property_value(cells, 'CellBackColor', light_grey_2)
        set_border(cells, 'bottom', OuterLineWidth=40, LineWidth=40)

        row = get_row(rixs, 0)
        set_width(row, 'optimal')
        for col in row.Columns:
            set_width(col, get_size(col)[0]*1.1)

        # XAS
        final = ['dataset'] + header_xas['basic'] + header_xas['pre_scan'] + ['error', 'Comments']
        set_row_data(xas, 0, final)
        set_row_data(xas, 1, final[1:-1], start_col=1)

        cells = get_cells(xas, position=(0, 0, len(final)+2, 0))
        set_property_value(cells, 'CharWeight', 200)
        set_property_value(cells, 'CharHeight', 12)
        set_property_value(cells, 'CellBackColor', light_grey_2)

        cells = get_cells(xas, position=(0, 1, len(final)+2, 1))
        set_property_value(cells, 'CellBackColor', light_grey_2)
        set_border(cells, 'bottom', OuterLineWidth=40, LineWidth=40)

        row = get_row(xas, 0)
        set_width(row, 'optimal')
        for col in row.Columns:
            set_width(col, get_size(col)[0]*1.1)


    ###############
    # unlock undo #
    ###############
    unlock_undo(doc)
    return

##########
# update #
##########
# C:\Users\galdin_c\Documents\work\CrTe2\apollo\2023_11_VERITAS_20230650\data\raw\CrTe2_DLD.h5
# C:\Users\galdin_c\Documents\work\CrTe2\apollo\2023_11_VERITAS_20230650\data\raw\CrTe2.h5
def check_document():
    """Check if doc is Calc and if Sheets named `settings`, `RIXS`, and `XAS` exist
    
    Args:
        None

    Returns:
        doc (document), settings (sheet), rixs (sheet), xas (sheet)  
    """
    # get document object
    doc = get_current_document()
    
    # check if document is calc
    if 'Sheets' not in doc.__dir__():
        msgbox('ERROR: This macro was supposed to be run in a Calc instance (spreadsheet)', type_msg='errorbox')
        return
    
    # check if calc has sheets
    sheet_names = get_sheet_names(doc)
    if 'settings' not in sheet_names or 'RIXS' not in sheet_names or 'XAS' not in sheet_names:
        msgbox('ERROR: Calc must have sheets named: `settings`, `RIXS`, and `XAS`', type_msg='errorbox')
        return

    # get sheets
    settings = get_sheet_by_name(doc, 'settings')
    rixs     = get_sheet_by_name(doc, 'RIXS')
    xas      = get_sheet_by_name(doc, 'XAS')

    return doc, settings, rixs, xas

def check_settings(settings):
    """Check if option in the settings sheet are sound

    Args:
        settings (sheet): sheet object to the settings sheet
    
    Returns:
        dictionary
    """
    address = {'RIXS filepath':         'B1',
               'XAS filepath':          'B2',
               'Save after update':     'B4',
               'Intercalate bkg color': 'B5',
               'RIXS: number of scans': 'J1',
               'RIXS: next row':        'J2',
               'XAS: number of scans':  'J3',
               'XAS: next row':         'J4'}

    ################
    # get settings #
    ################
    values = {}
    for text in address:
        values[text] = get_cells_data(settings, address[text])

    ###########################
    # validate INPUT settings #
    ###########################
    value = values['RIXS filepath']
    if value != '':
        if Path(value).exists() == False or Path(value).is_file() == False:
            msgbox(f'`RIXS filepath` in settings does not exist or do not point to a file', type_msg='errorbox')
            return
        
    value = values['XAS filepath']
    if value != '':
        if Path(value).exists() == False or Path(value).is_file() == False:
            msgbox(f'`XAS filepath` in settings does not exist or do not point to a file', type_msg='errorbox')
            return

    value = values['Save after update']
    if value not in ('yes', 'no'):
        msgbox(f"`Save after update` in settings cannot be {value}\n\navailable options = ('yes', 'no')", type_msg='errorbox')
        return 

    value = values['Intercalate bkg color']
    if value not in ('yes', 'no'):
        msgbox(f"`Intercalate bkg color` in settings cannot be {value}\n\navailable options = ('yes', 'no')", type_msg='errorbox')
        return 

    return values

def update(*args, **kwargs):
    """update macro for datafiles"""
    # get document object
    doc, settings, rixs, xas = check_document()
    
    # get settings from settings sheet
    values = check_settings(settings)

    ############
    # METADATA #
    ############
    for i, filepath in enumerate((values['RIXS filepath'], values['XAS filepath'])):
        if filepath != '':
            # scanlist
            scanlist = VERITAS.scanlist(filepath=filepath)

            # get sheet and next_row
            if i == 0:
                sheet = rixs
                next_row = int(values['RIXS: next row'])
            else:
                sheet = xas
                next_row = int(values['XAS: next row'])

            # heaser and last col
            header   = get_row_data(sheet, 1)
            last_col = get_last_used_col(sheet)

            # loop each scan and load metadata
            for _row, scan in enumerate(scanlist[int(next_row - 2):]):
                row = next_row

                # color row (intercalated)
                if values['Intercalate bkg color']:
                    temp = sheet.getCellRangeByPosition(0, row, last_col, row)
                    lock_undo(doc)
                    if row%2 == 0:
                        temp.setPropertyValue('CellBackColor', light_gold_4)
                    else:
                        temp.setPropertyValue('CellBackColor', light_indigo_4)
                    unlock_undo(doc)

                # paste on spreadsheet
                try:
                    pe = VERITAS.read(filepath, scan)
                    for col, attr in enumerate(header):
                        if attr != '':
                            if hasattr(pe, attr):
                                value = getattr(pe, attr)
                                if isinstance(value, datetime.datetime):
                                    value = str(value) 
                                elif isinstance(value, list) or isinstance(value, tuple):
                                    value = str(value)
                                lock_undo(doc)
                                set_cells_data(sheet, (col, row), value)
                                unlock_undo(doc)
                            else:
                                pass
                except Exception as e:
                    for col, attr in enumerate(header):
                        if attr != '':
                            if attr == 'scan':
                                value = scan
                            elif attr == 'scan_type':
                                value = 'err'
                            elif attr == 'error':
                                value = str(e)
                            else:
                                value = ''
                            lock_undo(doc)
                            set_cells_data(sheet, (col, row), value)
                            unlock_undo(doc)
                        else:
                            pass
                next_row += 1

            # update outputs
            lock_undo(doc)
            settings.unprotect('123456')
            if i == 0:
                set_cells_data(settings, 'J1', len(scanlist))
                set_cells_data(settings, 'J2', next_row)
            elif i == 1:
                set_cells_data(settings, 'J3', len(scanlist))
                set_cells_data(settings, 'J4', next_row)
            settings.protect('123456')
            unlock_undo(doc)

    #############
    # save file #
    #############
    if values['Save after update'] == 'yes':
        save(doc, verbose=False)

    return

#########
# OTHER #
#########
# C:\Users\galdin_c\Documents\current\2024_04_DLS-yuan\mock\mm30806-1
def build_header(*args, **kwargs):
    # get document object
    doc, settings, scanlist = check_document()

    # get settings from settings sheet
    folderpath, next_row, string, ref = check_settings(settings)

    # get filelist
    filelist = br.parsed_filelist(dirpath=folderpath, ref=ref, string=string)

    # get metadata
    file = filelist[-1]
    s    = read(file)

    # header
    header  = ['dataset', 'name']
    header += s.get_attrs()
    header += ['comments1', 'comments2']

    # set header
    set_row_data(scanlist, 0, header)
    set_row_data(scanlist, 1, header[2:-2], start_col=2)
    select_active_sheet(doc, scanlist)

    # format
    cells1 = get_cells(scanlist, position=(0, 0, len(header)+2, 0))
    set_property_value(cells1, 'CharWeight', 200)
    set_property_value(cells1, 'CharHeight', 12)
    set_property_value(cells1, 'CellBackColor', light_grey_2)

    cells2 = get_cells(scanlist, position=(0, 1, len(header)+2, 1))
    set_property_value(cells2, 'CellBackColor', light_grey_2)
    set_border(cells2, 'bottom', OuterLineWidth=40, LineWidth=40)

    set_width(get_row(scanlist, 0), 'optimal')
    return

def intercalate_row_bkg_colors(doc, scanlist):
    g = True

    last_col = get_last_used_col(scanlist) + 3
    for n in range(2, get_last_used_row(scanlist)):
        row = scanlist.getCellRangeByPosition(0, n, last_col, n)

        if g:
            msgbox(get_last_used_row(scanlist)+1)
            g = False
        if row.getPropertyValue('CellBackColor') == -1:
            lock_undo(doc)
            if n%2 == 0:
                row.setPropertyValue('CellBackColor', light_gold_4)
            else:
                row.setPropertyValue('CellBackColor', light_indigo_4)
            unlock_undo(doc)
    return

def metadata(*args, **kwargs):
    """Open msgbox with available metadata for datafiles"""

    # get document object
    doc, settings, scanlist = check_document()

    # get settings from settings sheet
    folderpath, next_row, string, ref = check_settings(settings)

    # get filelist
    filelist = br.parsed_filelist(dirpath=folderpath, ref=ref, string=string)

    # get metadata
    file = filelist[-1]
    s    = read(file)
    
    # msg box
    msgbox('\n'.join(s.get_attrs()))

########################
# cell property values #
########################
def property_values(*args, **kwargs):
    """Open msg box with important properties about the selection"""
    sel = get_current_selection()
    
    # text
    text = ''
    text += f"AbsoluteName: {get_property_value(sel, 'AbsoluteName')}\n"
    text += '\n'
    text +=  f"CellBackColor: {get_property_value(sel, 'CellBackColor')}\n"
    text +=  f"Size.Width: {get_size(sel)[0]}\n"
    text +=  f"Size.Height: {get_size(sel)[1]}\n"
    text += '\n'
    text +=  f"HoriJustify: {get_property_value(sel, 'HoriJustify')}\n"
    text +=  f"VertJustify: {get_property_value(sel, 'VertJustify')}\n"
    text += '\n'
    text +=  f"NumberFormat:      {get_property_value(sel, 'NumberFormat')}\n"
    # text +=  f"ConditionalFormat: {get_property_value(sel, 'ConditionalFormat')}\n"
    temp = get_protection(sel)['IsLocked']
    text +=  f"CellProtection.IsLocked: {temp}\n"
    text +=  f"CellStyle:         {get_property_value(sel, 'CellStyle')}\n"
    text += '-'*10 + '\n\n'
    
    text += 'Top border\n'
    props = get_border(sel, 'top')
    for name in props:
        text += f'{name}: {props[name]}\n'
    text += '-'*10 + '\n\n'

    text += 'Bottom border\n'
    props = get_border(sel, 'bottom')
    for name in props:
        text += f'{name}: {props[name]}\n'
    text += '-'*10 + '\n\n'

    text += 'Right border\n'
    props = get_border(sel, 'right')
    for name in props:
        text += f'{name}: {props[name]}\n'
    text += '-'*10 + '\n\n'

    text += 'Left border\n'
    props = get_border(sel, 'left')
    for name in props:
        text += f'{name}: {props[name]}\n'
    text += '-'*10 + '\n\n'

    props = ['IsTextWrapped', 'CharHeight', 'CharWeight', 'CharFontName', 'CharColor', 
             'CharStrikeout', 'CharUnderlineHasColor', 'CharUnderlineColor', 
             'CharUnderline']
    for name in props:
        try:
            text += f'{name}: {get_property_value(sel, name)}\n'
        except:
            text += f'{name}: I do not know\n'
    text += '-'*10 + '\n\n'

    # msgbox
    msgbox(text)
    
    return

# %%




# %% ---------------------- Import function as Macro ---------------------- %% #
# Only the specified functions will show in the Tools > Macro > Organize Macro dialog
g_exportedScripts = (metadata, update, new_empty, build_header, property_values)
# g_exportedScripts = (metadata, update, new_empty, build_header)

# %%

# ============================================================================ #
# =========================== MACRO - I21 - 2024-04-27 ======================= #
# %% ====================================================================== %% #
# def property_values(*args, **kwargs):
#     """Open msg box with important properties about the selection"""
#     sel = get_current_selection()
    
#     # text
#     text = ''
#     text += f"AbsoluteName: {get_property_value(sel, 'AbsoluteName')}\n"
#     text += '\n'
#     text +=  f"CellBackColor: {get_property_value(sel, 'CellBackColor')}\n"
#     text +=  f"Size.Width: {get_size(sel)[0]}\n"
#     text +=  f"Size.Height: {get_size(sel)[1]}\n"
#     text += '\n'
#     text +=  f"HoriJustify: {get_property_value(sel, 'HoriJustify')}\n"
#     text +=  f"VertJustify: {get_property_value(sel, 'VertJustify')}\n"
#     text += '\n'
#     text +=  f"NumberFormat:      {get_property_value(sel, 'NumberFormat')}\n"
#     # text +=  f"ConditionalFormat: {get_property_value(sel, 'ConditionalFormat')}\n"
#     temp = get_protection(sel)['IsLocked']
#     text +=  f"CellProtection.IsLocked: {temp}\n"
#     text +=  f"CellStyle:         {get_property_value(sel, 'CellStyle')}\n"
#     text += '-'*10 + '\n\n'
    
#     text += 'Top border\n'
#     props = get_border(sel, 'top')
#     for name in props:
#         text += f'{name}: {props[name]}\n'
#     text += '-'*10 + '\n\n'

#     text += 'Bottom border\n'
#     props = get_border(sel, 'bottom')
#     for name in props:
#         text += f'{name}: {props[name]}\n'
#     text += '-'*10 + '\n\n'

#     text += 'Right border\n'
#     props = get_border(sel, 'right')
#     for name in props:
#         text += f'{name}: {props[name]}\n'
#     text += '-'*10 + '\n\n'

#     text += 'Left border\n'
#     props = get_border(sel, 'left')
#     for name in props:
#         text += f'{name}: {props[name]}\n'
#     text += '-'*10 + '\n\n'

#     props = ['IsTextWrapped', 'CharHeight', 'CharWeight', 'CharFontName', 'CharColor', 
#              'CharStrikeout', 'CharUnderlineHasColor', 'CharUnderlineColor', 
#              'CharUnderline']
#     for name in props:
#         try:
#             text += f'{name}: {get_property_value(sel, name)}\n'
#         except:
#             text += f'{name}: I do not know\n'
#     text += '-'*10 + '\n\n'

#     # msgbox
#     msgbox(text)
    
#     return

# def new_empty(*args, **kwargs):
#     """Create new template spreadsheet"""
#     doc = new_calc()
#     lock_undo(doc)
#     settings = new_sheet(doc, 'settings', i=0)
#     scanlist = new_sheet(doc, 'scanlist', i=1)
#     remove_sheet_by_name(doc, 'Sheet1')
#     select_active_sheet(doc, settings)

#     # doc      = get_current_document()
#     # settings = get_sheet_by_name(doc, 'settings')

#     # unprotect
#     settings.unprotect('123456')

#     ##################
#     # settings sheet #
#     ##################

#     # bkg 
#     cells = get_cells(settings, 'A1:F6')
#     set_property_value(cells, 'CellBackColor', light_grey_2)
#     set_border(cells=cells, border='all', OuterLineWidth=26, LineStyle=0, LineWidth=26)
#     set_width(cells, 4000)

#     cells = get_cells(settings, 'I1:I6')
#     set_property_value(cells, 'CellBackColor', gold)
#     set_width(cells, 3000)
#     set_border(cells=cells, border='all', OuterLineWidth=26, LineStyle=0, LineWidth=26)

#     cells = get_cells(settings, 'J1:J6')
#     set_property_value(cells, 'CellBackColor', light_gold_4)
#     set_border(cells=cells, border='all', OuterLineWidth=26, LineStyle=0, LineWidth=26)

#     cells = get_cells(settings, 'G1:H1')
#     set_width(cells, 1000)

#     # INSTRUCTUIONS
#     cells = get_cells(settings, 'I8:J22')
#     merge(cells)
#     set_property_value(cells, 'CellBackColor', light_gold_4)
#     set_property_value(cells, 'IsTextWrapped', True)
#     set_property_value(cells, 'VertJustify', 1)
#     text  = 'This sheet is password protected to avoid mistakes\n\n'
#     text += 'The password is: 123456\n\n'
#     text += 'To unlock it:\nTools > Protect Sheet > Type in password'
#     set_cells_data(settings, 'I8', text)

#     cells = get_cells(settings, 'A8:F22')
#     merge(cells)
#     set_property_value(cells, 'CellBackColor', light_gold_4)
#     set_property_value(cells, 'IsTextWrapped', True)
#     set_property_value(cells, 'VertJustify', 1)
#     text  = 'folderpath/filepath: path to folder with datafiles or path to h5 file with multiple scans\n\n'
#     text += 'filter string: Only usable for folder with multiple datafiles. Filter string filter files with the defined string in its filename\n\n'
#     text += "filename ref: let's say your datafiles are named `I21-scan134058.nxs`. The number 134058 is the scan number. "
#     text += "you want to order the files in order of scan number. The internal filename parser will isolate every number "
#     text += "sequence from the name. For this filename, the parser will get: [21, 134058]. "
#     text += "You know that the second number that shows up in the filename is the scan number, therefore, the ref number "
#     text += "should be 1 (starts from zero). In your filename is something like this `scan134058.dat`, then ref must be 0.\n\n"
#     text += "Save after update: Spreadsheet is saved after running the `update` macro. This settings cannot be changed for now.\n\n" 
#     text += "Intercalate bkg color: Scan lines have alternating bkg color. This settings cannot be changed for now." 
#     set_cells_data(settings, 'A8', text)

#     cells = get_cells(settings, 'A24:F34')
#     merge(cells)
#     set_property_value(cells, 'CellBackColor', light_indigo_4)
#     set_property_value(cells, 'IsTextWrapped', True)
#     set_property_value(cells, 'VertJustify', 1)
#     text  = 'Reminders:\n\n'
#     text += '- Turn grid for colored cells: Tools > Options > LibreOffice Calc > View > Visual Aids > Grid Lines > select: `Show on colored cells`\n\n'
#     text += "- Freeze First rows/cols: View > Freeze Cells > Freeze First column/Row\n\n"
#     text += "- Freeze multiple rows/cols: Select the cell which left/above that cell should freeze > View > Freeze Rows and Columns\n\n"
#     text += "- Move columns: Select a column by clicking the column header > hold down Alt and drag the column by one of its cells (it you try to drag from the header it will not work)\n\n"
#     set_cells_data(settings, 'A24', text)

#     # folderpath
#     set_cells_data(settings, 'A1', 'folderpath/filepath')
#     cells = get_cells(settings, 'B1:F1')
#     merge(cells)

#     # save after update
#     set_cells_data(settings, 'A2', 'Save after update')
#     set_cells_data(settings, 'B2', 'yes')

#     # Intercalate bkg color
#     set_cells_data(settings, 'A3', 'Intercalate bkg color')
#     set_cells_data(settings, 'B3', 'yes')

#     # filter string
#     set_cells_data(settings, 'D2', 'filter string')
#     set_cells_data(settings, 'E2', '.nxs')

#     # filename ref
#     set_cells_data(settings, 'D3', 'filename ref')
#     set_cells_data(settings, 'E3', 0)

#     # editable
#     for position in ('B1:G1', 'B2', 'B3', 'E2', 'E3'):
#         cells = get_cells(settings, position)
#         unlock_cells(cells)
#         set_property_value(cells, 'CellBackColor', -1)
#     for position in ('B2', 'B3'):
#         cells = get_cells(settings, position)
#         lock_cells(cells) 

#     # non-editable
#     set_cells_data(settings, 'I1', 'number of files')
#     set_cells_data(settings, 'J1', 0)

#     set_cells_data(settings, 'I2', 'next row')
#     set_cells_data(settings, 'J2', 2)

#     # protect
#     settings.protect('123456')

#     # unlock undo
#     unlock_undo(doc)
#     return

# # C:\Users\galdin_c\Documents\current\2024_04_DLS-yuan\mock\mm30806-1
# def build_header(*args, **kwargs):
#     # get document object
#     doc, settings, scanlist = check_document()

#     # get settings from settings sheet
#     folderpath, next_row, string, ref = check_settings(settings)

#     # get filelist
#     filelist = br.parsed_filelist(dirpath=folderpath, ref=ref, string=string)

#     # get metadata
#     file = filelist[-1]
#     s    = read(file)

#     # header
#     header  = ['dataset', 'name']
#     header += s.get_attrs()
#     header += ['comments1', 'comments2']

#     # set header
#     set_row_data(scanlist, 0, header)
#     set_row_data(scanlist, 1, header[2:-2], start_col=2)
#     select_active_sheet(doc, scanlist)

#     # format
#     cells1 = get_cells(scanlist, position=(0, 0, len(header)+2, 0))
#     set_property_value(cells1, 'CharWeight', 200)
#     set_property_value(cells1, 'CharHeight', 12)
#     set_property_value(cells1, 'CellBackColor', light_grey_2)

#     cells2 = get_cells(scanlist, position=(0, 1, len(header)+2, 1))
#     set_property_value(cells2, 'CellBackColor', light_grey_2)
#     set_border(cells2, 'bottom', OuterLineWidth=40, LineWidth=40)

#     set_width(get_row(scanlist, 0), 'optimal')
#     return

# def check_document():
#     """Check if doc is Calc and if Sheets named `settings` and `scalist` exist
    
#     Args:
#         None

#     Returns:
#         doc (document), settings (sheet) scanlist (sheet) objects  
#     """
#     # get document object
#     doc = get_current_document()
    
#     # check if document is calc
#     if 'Sheets' not in doc.__dir__():
#         msgbox('ERROR: This macro was supposed to be run in a Calc instance (spreadsheet)', type_msg='errorbox')
#         return
    
#     # check if calc has two sheets
#     # they must be named settings and scanlist
#     sheet_names = get_sheet_names(doc)
#     if 'settings' not in sheet_names or 'scanlist' not in sheet_names:
#         msgbox('ERROR: Calc must have two sheets named: `settings` and `scanlist`', type_msg='errorbox')
#         return

#     # get sheets
#     settings = get_sheet_by_name(doc, 'settings')
#     scanlist = get_sheet_by_name(doc, 'scanlist')

#     return doc, settings, scanlist

# def check_settings(settings):
#     """Check if option in the settings sheet are sound

#     Args:
#         settings (sheet): sheet object to the settings sheet
    
#     Returns:
#         folderpath, next_row
#     """
#     # get settings
#     folderpath = Path(get_cells_data(settings, 'B1'))
#     # time       = get_cells_data(settings, 'B2')
#     next_row   = get_cells_data(settings, 'J2')
#     string     = get_cells_data(settings, 'E2')
#     ref        = get_cells_data(settings, 'E3')


#     # validate settings
#     if folderpath == Path(''):
#         msgbox(f'folderpath does not exist\n\nfolderpath = {folderpath}', type_msg='errorbox')
#         return
#     if folderpath.exists() == False:
#         msgbox(f'folderpath does not exist\n\nfolderpath = {folderpath}', type_msg='errorbox')
#         return 
    
#     # try: 
#     #     time = float(time)
#     # except ValueError: 
#     #     msgbox(f'update time must be a number\n\nupdate time = {time}')
#     # if time < 5:
#     #     msgbox(f'update time is too little (minimum is 5s)\n\nupdate time = {time}')
#     #     return
    
#     if next_row == '':
#         next_row = 2
#     else:
#         try: 
#             next_row = int(next_row)
#         except ValueError: 
#             msgbox(f'next row must be a number\n\nnext row = {next_row}\n\n')
#     if next_row < 2:
#         msgbox(f'next row must be a number >= 2\n\nnext row = {next_row}\n\n')


#     if ref == '':
#         ref = 0
#     else:
#         try: 
#             ref = int(ref)
#         except ValueError: 
#             msgbox(f'filename ref must be a number\n\nref = {ref}\n\n')
#     if next_row < 0:
#         msgbox(f'filename ref must be a number >= 0\n\nref = {next_row}\n\n')

#     return folderpath, next_row, string, ref

# def intercalate_row_bkg_colors(doc, scanlist):
#     g = True

#     last_col = get_last_used_col(scanlist) + 3
#     for n in range(2, get_last_used_row(scanlist)):
#         row = scanlist.getCellRangeByPosition(0, n, last_col, n)

#         if g:
#             msgbox(get_last_used_row(scanlist)+1)
#             g = False
#         if row.getPropertyValue('CellBackColor') == -1:
#             lock_undo(doc)
#             if n%2 == 0:
#                 row.setPropertyValue('CellBackColor', light_gold_4)
#             else:
#                 row.setPropertyValue('CellBackColor', light_indigo_4)
#             unlock_undo(doc)
#     return

# def update(*args, **kwargs):
#     """update macro for datafiles"""
#     # get document object
#     doc, settings, scanlist = check_document()
    
#     # get settings from settings sheet
#     folderpath, next_row, string, ref = check_settings(settings)

#     # get header
#     header = get_row_data(scanlist, 1)

#     # get filelist
#     if Path(folderpath).is_dir():
#         filelist = br.parsed_filelist(dirpath=folderpath, ref=ref, string=string)
#     else:
#         # open h5 file for scans
#         msgbox('Not implemented error. spreadsheet only work for folderpath, not filepath')
#         return
    
#     # get attrs and place on the table
#     last_col = get_last_used_col(scanlist)
#     for _row, file in enumerate(filelist[next_row-2:]):
#         row = next_row

#         # color row (intercalated)
#         temp = scanlist.getCellRangeByPosition(0, row, last_col, row)
#         lock_undo(doc)
#         if row%2 == 0:
#             temp.setPropertyValue('CellBackColor', light_gold_4)
#         else:
#             temp.setPropertyValue('CellBackColor', light_indigo_4)
#         unlock_undo(doc)


#         # get metadata
#         s = read(file)

#         # paste on spreadsheet
#         for col, attr in enumerate(header):
#             if attr != '':
#                 if hasattr(s, attr):
#                     value = getattr(s, attr)
#                     lock_undo(doc)
#                     set_cells_data(scanlist, (col, row), value)
#                     unlock_undo(doc)
#                 else:
#                     pass
#         next_row += 1

#     lock_undo(doc)
#     settings.unprotect('123456')
#     set_cells_data(settings, 'J1', len(filelist))
#     set_cells_data(settings, 'J2', next_row)
#     settings.protect('123456')
#     unlock_undo(doc)

#     save(doc, verbose=False)

#     return

# def metadata(*args, **kwargs):
#     """Open msgbox with available metadata for datafiles"""

#     # get document object
#     doc, settings, scanlist = check_document()

#     # get settings from settings sheet
#     folderpath, next_row, string, ref = check_settings(settings)

#     # get filelist
#     filelist = br.parsed_filelist(dirpath=folderpath, ref=ref, string=string)

#     # get metadata
#     file = filelist[-1]
#     s    = read(file)
    
#     # msg box
#     msgbox('\n'.join(s.get_attrs()))
# # %%