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

This macro expects that the Calc file has two Sheets called 'options' and 
'scanlist'.

The `update` macro will read the files in a folder defined in the 'options' 
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
path2brixs = Path(r'C:\Users\galdin_c\github\brixs_raw')
sys.path.append(str(path2brixs))
import brixs as br
import brixs.beamlines.I21 as I21

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

    # get filepath
    if filepath is None:
        filepath = get_filepath(doc)
        if filepath == '':
            msgbox('filepath not defined. Cannot save file via macro. Please, save document by hand at least once first.', type_msg='errorbox')
            return
        filepath = Path(filepath).with_suffix(extension)
    else:
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
        None
    """
    # position
    if i is None:
        i = doc.get_sheet_count() + 1
    assert i >= 0, 'position i must be integer >= 0'

    # unique name
    if name in doc.get_sheet_names():
        msgbox(f'Cannot create sheet with name {name} because there is a sheet with this name already', type_msg='errorbox')
        return
    
    # create sheet
    doc.Sheets.insertNewByName(name, i)

    return

def remove_sheet_by_name(doc, name):
    """delete sheet by name"""
    if name in doc.get_sheet_names():
        if doc.get_sheet_count() > 1:
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
        i = doc.get_sheet_count() + 1
    assert i >= 0, 'position i must be integer >= 0'

    # check new name
    names = get_sheet_names(doc)
    if new_name in names:
        if name in doc.get_sheet_names():
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
# %%

# %% -------------------------- sheet functions --------------------------- %% #
def get_cell_by_position(sheet, position):
    """return cell object from position = 'A2' or position = (0, 2) - starts from 0"""
    if isinstance(position, str):
        if ':' in position:
            msgbox(f'get_cell_by_position error: input must be a cell position, not range\n\ninput = {position}', type_msg='errorbox')
            return
        position = str2num(position)
    return sheet.getCellByPosition(*position)

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

def get_row_data(sheet, row):
    """return list with row data (row number starting from 0)"""
    return sheet.getCellRangeByPosition(0, row, get_last_used_col(sheet), row).getDataArray()[0]

def get_col_data(sheet, col):
    """return list with column data (column number starting from 0)"""
    if isinstance(col, str):
        if ':' in col:
            msgbox(f'get_col_data error: input must be a column position, not range\n\ninput = {position}', type_msg='errorbox')
            return
        position = str2num(position)
    return [_[0] for _ in sheet.getCellRangeByPosition(col, 0, col, get_last_used_row(sheet)).getDataArray()]

def get_cell_data(sheet, position):
    """return cell data. Position = 'A2' or (0, 2). Starts from 0"""
    if isinstance(position, str):
        if ':' in position:
            msgbox(f'get_cell_data error: input must be a cell position, not range\n\ninput = {position}')
            return
        position = str2num(position)
    return sheet.getCellRangeByPosition(position[0], position[1], position[0], position[1]).getDataArray()[0][0]

def set_cell_data(sheet, position, value):
    """set cell data. Position = 'A2' or (0, 2). Starts from 0"""
    if isinstance(position, str):
        if ':' in position:
            msgbox(f'get_cell_data error: input must be a cell position, not range\n\ninput = {position}', type_msg='errorbox')
            return
        position = str2num(position)
    return sheet.getCellRangeByPosition(position[0], position[1], position[0], position[1]).setDataArray([[value, ], ])

def get_range_data(sheet, position):
    """return list of list. First index is the row, second index is the col
    
    Position = 'A2"C5' or (col_start, row_start, col_stop, row_stop). Starts from 0
    """
    if isinstance(position, str):
        if ':' not in position:
            msgbox(f'get_range_data error: input must be a cell range\n\ninput = {position}', type_msg='errorbox')
            return
        position = str2num(position)
    return sheet.getCellRangeByPosition(position[0], position[1], position[2], position[3]).getDataArray()
# %%

# %% ----------------------- cell/ranges functions ------------------------ %% #
def get_property_names(cell):
    """return property names from cell or ranges"""
    return [x.Name for x in cell.getPropertySetInfo().Properties]

def get_property_value(cell, name):
    """return cell/ranges property value"""
    return cell.getPropertyValue(name)

def get_current_selection():
    """return a cell or range object"""
    desktop = XSCRIPTCONTEXT.getDesktop()
    return desktop.CurrentComponent.CurrentController.getSelection()
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
    im.scan       = int(filepath.name.split('-')[-1].split('.')[0])
    im.start_time = str(im.start_time)
    im.end_time   = str(im.end_time)

    return im

def read_spectrum(filepath):
    return br.Spectrum(filepath)

# assign the correct read function to the name `read`
read = read_I21
# %%

# ============================================================================ #
# ================================== MACRO =================================== #
# %% ====================================================================== %% #
light_yellow = 16774606
light_purple = 14605542

def create_template():
    msgbox('This was not implemented yet.') 
    return

def check_document():
    """Check if doc is Calc and if Sheets named `options` and `scalist` exist
    
    Args:
        None

    Returns:
        doc (document), options (sheet) scanlist (sheet) objects  
    """
    # get document object
    doc   = get_current_document()
    
    # check if document is calc
    if 'Sheets' not in doc.__dir__():
        msgbox('ERROR: This macro was supposed to be run in a Calc instance (spreadsheet)', type_msg='errorbox')
        return
    
    # check if calc has two sheets
    # they must be named options and scanlist
    sheet_names = get_sheet_names(doc)
    if 'options' not in sheet_names or 'scanlist' not in sheet_names:
        msgbox('ERROR: Calc must have two sheets named: `options` and `scanlist`', type_msg='errorbox')
        return

    # get sheets
    options  = get_sheet_by_name(doc, 'options')
    scanlist = get_sheet_by_name(doc, 'scanlist')

    return doc, options, scanlist

def check_options(options):
    """Check if option in the options sheet are sound

    Args:
        options (sheet): sheet object to the options sheet
    
    Returns:
        folderpath, next_row
    """
    # get options
    folderpath = Path(get_cell_data(options, 'B1'))
    # time       = get_cell_data(options, 'B2')
    next_row   = get_cell_data(options, 'J2')

    # validate options
    if folderpath.exists() == False:
        msgbox(f'folderpath does not exist\n\nfolderpath = {folderpath}')
        return 
    
    # try: 
    #     time = float(time)
    # except ValueError: 
    #     msgbox(f'update time must be a number\n\nupdate time = {time}')
    # if time < 5:
    #     msgbox(f'update time is too little (minimum is 5s)\n\nupdate time = {time}')
    #     return
    
    if next_row == '':
        next_row = 2
    else:
        try: 
            next_row = int(next_row)
        except ValueError: 
            msgbox(f'next row must be a number\n\nnext row = {next_row}\n\n')
    if next_row < 2:
        msgbox(f'next row must be a number >= 2\n\nnext row = {next_row}\n\n')

    return folderpath, next_row

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
                row.setPropertyValue('CellBackColor', light_yellow)
            else:
                row.setPropertyValue('CellBackColor', light_purple)
            unlock_undo(doc)
    return

def update(*args, **kwargs):
    """update macro for datafiles"""
    # get document object
    doc, options, scanlist = check_document()
    
    # get options from options sheet
    folderpath, next_row   = check_options(options)

    # get header
    header = get_row_data(scanlist, 1)

    # get filelist
    filelist = br.parsed_filelist(dirpath=folderpath, ref=1, string='.nxs')
    
    # get attrs and place on the table
    last_col = get_last_used_col(scanlist) + 3
    for _row, file in enumerate(filelist[next_row-2:]):
        row = next_row

        # color row (intercalated)
        temp = scanlist.getCellRangeByPosition(0, row, last_col, row)
        lock_undo(doc)
        if row%2 == 0:
            temp.setPropertyValue('CellBackColor', light_yellow)
        else:
            temp.setPropertyValue('CellBackColor', light_purple)
        unlock_undo(doc)


        # get metadata
        s = read(file)

        # paste on spreadsheet
        for col, attr in enumerate(header):
            if attr != '':
                if hasattr(s, attr):
                    value = getattr(s, attr)
                    lock_undo(doc)
                    set_cell_data(scanlist, (col, row), value)
                    unlock_undo(doc)
                else:
                    pass
        next_row += 1

    lock_undo(doc)
    set_cell_data(options, 'J1', len(filelist))
    set_cell_data(options, 'J2', next_row)
    unlock_undo(doc)

    save(doc, verbose=False)

    return

def metadata(*args, **kwargs):
    """Open msgbox with available metadata for datafiles"""

    # get document object
    doc, options, scanlist = check_document()

    # get options
    folderpath = Path(get_cell_data(options, 'B1'))

    # get filelist
    filelist = br.parsed_filelist(dirpath=folderpath, ref=1, string='.nxs')

    # get metadata
    file = filelist[-1]
    s    = read(file)
    
    # add new attrs
    s.scan = int(file.name.split('-')[-1].split('.')[0])

    # msg box
    msgbox('\n'.join(s.get_attrs()))
# %%

# %% ---------------------- Import function as Macro ---------------------- %% #
# Only the specified functions will show in the Tools > Macro > Organize Macro dialog
g_exportedScripts = (metadata, update)
# %%