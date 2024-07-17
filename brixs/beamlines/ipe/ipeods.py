"""Macro for creating an online scanlist for experiments

###############
# Last edited #
###############
galdino 27-06-2024

############
# contents #
############
See below instruction on how to use/install it

Items 1 -- 5  how to install this macro
Items 6 -- 8  how to use this macro
Item  9       how to write a new macro or edit this one
 
#######################################
# 1. Download and install libreoffice #
#######################################

Follow the instructions on the libreoffice webpage, download and install it.

############################################
# 2a. Install pip on libreoffice (WINDOWS) #
############################################

Option 1: via zaz-pip libreoffice extension

1. Go to zaz-pip git page here -> https://git.cuates.net/elmau/zaz-pip
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

Option 1: via zaz-pip libreoffice extension 

1. Go to zaz-pip git page here -> https://git.cuates.net/elmau/zaz-pip
2. Follow the instructions on the git page for installation
3. Follow instructions on the git page for installing packages

Option 2: via pypi bootstrap [I never tested this option!!]

1. Download the pypi bootstrap here -> https://bootstrap.pypa.io/get-pip.py
2. Start a command prompt at the LibreOfice installation directory (e.g., cd /usr/lib/libreoffice/program)
3. Run python <path-to-get-pip>\get-pip.py
4. Use command python -m pip install <package-name> for installing packages

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
If changes are made to brixs, then you don't need to install it again inside the 
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

View > Toolbars > Customize

Select the Toolbars tab

set Category to Macros:Application Macros (left panel)

Move the desired Macros from the left panel o the right panel (assigned commands)
by pressing the right arrow between the left and right panels.

Down below the Customize windows there is a button called `insert`. Use this to 
insert a separator. You can also press the button `Modify` to change the name of
 the macro or put an icon. 

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
import numpy as np
import datetime
import sys

# %% -------------------------- brixs imports ----------------------------- %% #
# path2brixs = Path(r'/ibira/lnls/beamlines/ipe/proposals/20232623/proc/brixs/')
path2brixs = Path(r'/usr/local/scripts/apps/brixs/')
sys.path.append(str(path2brixs))
import brixs as br
import brixs.beamlines.IPE2 as ipe

# common
from brixs.sheets.common import letter2num, str2num

# libreoffice
from brixs.sheets.libreoffice import create_instance, msgbox, colors
from brixs.sheets.libreoffice import new_calc, new_writer, new_impress, new_draw, new_math, new_base, open_document
from brixs.sheets.libreoffice import get_filepath, get_filename, save, close, lock_undo, unlock_undo

# ods
from brixs.sheets.ods import get_sheet_count, get_sheet_names 
from brixs.sheets.ods import get_sheet_by_position, get_sheet_by_name, get_active_sheet
from brixs.sheets.ods import new_sheet, remove_sheet_by_name, remove_sheet_by_position, move_sheet, copy_sheet
from brixs.sheets.ods import get_style_names, remove_style, new_style, select_active_sheet
from brixs.sheets.ods import get_cells, get_last_used_row, get_last_used_col
from brixs.sheets.ods import delete_rows
from brixs.sheets.ods import merge, unmerge
from brixs.sheets.ods import get_row, get_row_data, get_col_data, set_row_data, get_cells_data, set_cells_data
from brixs.sheets.ods import get_property_names, get_property_value, set_property_value
from brixs.sheets.ods import get_width, get_height, set_width, set_height
from brixs.sheets.ods import get_border, set_border
from brixs.sheets.ods import get_protection, lock_cells, unlock_cells
# %%

import brixs.addons.fitting

# %% -------------------------- uno definitions --------------------------- %% #
def get_current_document():
    """Return document object"""
    return XSCRIPTCONTEXT.getDocument()

def get_current_selection():
    """return a cell or range object"""
    desktop = XSCRIPTCONTEXT.getDesktop()
    return desktop.CurrentComponent.CurrentController.getSelection()
# %%

# %% --------------------------- test functions --------------------------- %% #
def test1(*args, **kwargs):
    msgbox('test 1')

def test2(*args, **kwargs):
    doc = get_current_document()
    msgbox(get_sheet_count(doc))
# %%

# %% ------------------------ support functions --------------------------- %% #
def _input(sheet, text, text_position, input_position, input_default='', editable=True, text_width='optimal', input_width=False):
    """quickly set 'input' cells

    Args:
        sheet (sheet): sheet object
        text (str): help text
        text_position (str): cell address to put help text. If range, cells are merged
        input_position (str): cell address for input. If range, cells are merged
        input_default (str): default text on input cells
        editable (bool): if False, input cells will be protected when sheet is locked
        text_width (number, 'optimal', or False, optional): help text width. 
            If merged, width is applied to each column. Default is 'optimal'
        input_width (number, 'optimal', or False, optional): input cells width. 
            If merged, width is applied to each column. Default is 'False

    Usage:
        _input(sheet=settings, text='RIXS folderpath', text_position='A1', 
               input_position='B1:F1', input_default='', 
               editable=True, 
               input_width=4000)

    Returns:
        None
    """
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
    """quickly set 'output' cells

    Args:
        sheet (sheet): sheet object
        text (str): help text
        text_position (str): cell address to put help text. If range, cells are merged
        out_position (str): cell address for output. If range, cells are merged
        out_default (str): default text on output cells
        editable (bool): if False, input cells will be protected when sheet is locked
        text_width (number, 'optimal', or False, optional): help text width. 
            If merged, width is applied to each column. Default is 'optimal'
        out_width (number, 'optimal', or False, optional): ouput cells width. 
            If merged, width is applied to each column. Default is 'False

    Usage:
        _output(sheet=settings, text='RIXS: number of scans', 
                text_position='I1', out_position='J1', out_default0)

    Returns:
        None
    """
    # help text
    if isinstance(text_position, str):
        text_position = str2num(text_position)
    cells = get_cells(sheet, text_position)
    if len(text_position) == 4:
        merge(cells)
    set_cells_data(sheet, text_position[:2], text)
    set_property_value(cells, 'CellBackColor', colors['gold'])
    if text_width:
        set_width(cells, text_width)

    # input
    if isinstance(out_position, str):
        out_position = str2num(out_position)
    cells = get_cells(sheet, out_position)
    if len(out_position) == 4: 
        merge(cells)
    set_cells_data(sheet, out_position[:2], out_default)      # default input
    set_property_value(cells, 'CellBackColor', colors['light_gold_4'])  # set color
    if editable:
        unlock_cells(cells)                                   # unlock cells
    if out_width:
        set_width(cells, out_width)

def _instruct(sheet, text, position):
    """quickly set 'instruction' cells

    Args:
        sheet (sheet): sheet object
        text (str): help text
        position (str): cell address to put help text. If range, cells are merged

    Usage:
        _output(sheet=settings, text='RIXS: number of scans', 
                text_position='I1', out_position='J1', out_default0)

    Returns:
        None
    """
    if isinstance(position, str):
        position = str2num(position)

    cells = get_cells(sheet, position)

    if len(position) == 4:
        merge(cells)
    set_property_value(cells, 'CellBackColor', colors['light_gold_4'])
    set_property_value(cells, 'IsTextWrapped', True)
    set_property_value(cells, 'VertJustify', 1)
    set_cells_data(sheet, position[:2], text)

def _check_document():
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
    # _1D      = get_sheet_by_name(doc, 'ASCAN')
    # _2D      = get_sheet_by_name(doc, 'MESH')

    return doc, settings, rixs, xas#, _1D, _2D

def _check_settings(settings):
    """Check if options in the settings sheet are sound

    Args:
        settings (sheet): sheet object to the settings sheet
    
    Returns:
        dictionary
    """
    address = {'RIXS folderpath':         'B1',
               'XAS folderpath':          'B2',
               'Save after update':     'B4',
               'Intercalate bkg color': 'B5',
               'Reset data':            'E4',
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
    value = values['RIXS folderpath']
    if value != '':
        if Path(value).exists() == False or Path(value).is_dir() == False:
            msgbox(f'`RIXS folderpath` in settings does not exist or do not point to a folder', type_msg='errorbox')
            return
        
    value = values['XAS folderpath']
    if value != '':
        if Path(value).exists() == False or Path(value).is_dir() == False:
            msgbox(f'`XAS folderpath` in settings does not exist or do not point to a folder', type_msg='errorbox')
            return

    value = values['Save after update']
    if value not in ('yes', 'no'):
        msgbox(f"`Save after update` in settings cannot be {value}\n\navailable options = ('yes', 'no')", type_msg='errorbox')
        return 

    value = values['Intercalate bkg color']
    if value not in ('yes', 'no'):
        msgbox(f"`Intercalate bkg color` in settings cannot be {value}\n\navailable options = ('yes', 'no')", type_msg='errorbox')
        return 
    
    value = values['Reset data']
    if value not in ('yes', 'no'):
        msgbox(f"`Reset data` in settings cannot be {value}\n\navailable options = ('yes', 'no')", type_msg='errorbox')
        return 

    return values
# %%

# ============================================================================ #
# ============================= MACRO - IPE ================================== #
# %% ====================================================================== %% #

# %% rixs metadata list
rixs_attrs = []
for key in ipe.rixs_attrs:
    for attr in ipe.rixs_attrs[key]:
        rixs_attrs.append(attr)
rixs_attrs = list(np.sort(rixs_attrs))

# %% xas metadata list
# xas_attrs = ['scan', 'start_time', 'command', 'energy', 'ring_current', 'cff', 'GR', 'MR', 'GR_offset', 'MR_offset', 'GT', 'MT', 'line_density', 'und_phaseMon', 'FOE_right', 'FOE_left', 'FOE_top', 'FOE_bottom', 'FOE_Vgap', 'FOE_Voffset', 'FOE_Hgap', 'FOE_Hoffset', 'WBS_right', 'WBS_left', 'WBS_top', 'WBS_bottom', 'WBS_Vgap', 'WBS_Voffset', 'WBS_Hgap', 'WBS_Hoffset', 'MPS1A_Vgap', 'MPS1A_Hgap', 'MPS2B_Vgap', 'MPS2B_Hgap', 'WBS_MKS', 'M1_MKS', 'PGM_MKS', 'M4_MKS', 'M5_MKS', 'MPS1A_MKS', 'MPS2B_MKS', 'M7_MKS', 'WBS_pump', 'M1_pump', 'PVS_pump', 'PGM_pump', 'MVS1_pump', 'M4_pump', 'MVS2_pump', 'M5_pump', 'PGshutter_pump', 'TPA_1A_pump', 'TPA_1B_pump', 'TPB_1A_pump', 'TPB_1B_pump', 'MPS1A_pump', 'MPS2B_pump', 'TPA_2A_pump', 'TPA_2B_pump', 'TPB_2A_pump', 'TPB_2B_pump', 'MVS3A_pump', 'MVS4B_pump', 'M7_pump', 'M1_X', 'M1_Y', 'M1_Z', 'M1_Rx', 'M1_Ry', 'M1_Rz', 'M4_X', 'M4_Y', 'M4_Z', 'M4_Rx', 'M4_Ry', 'M4_Rz', 'M4_Uy', 'M4_Ry_piezo', 'M5_X', 'M5_Y', 'M5_Z', 'M5_Rx', 'M5_Ry', 'M5_Rz', 'M5_Uy', 'M5_Ry_piezo', 'M6_X', 'M6_Y', 'M6_Z', 'M6_Rx', 'M6_Ry', 'M6_Rz', 'M6_Uy', 'M6_Rx_piezo', 'M6_Ry_piezo', 'M7_X', 'M7_Y', 'M7_Z', 'M7_Rx', 'M7_Ry', 'M7_Rz', 'M7_Uy', 'M7_Rx_piezo', 'M7_Ry_piezo', 'M1_collector', 'M1_collector_rng', 'MVS1', 'MVS1_rng', 'MVS2', 'MVS2_rng', 'MVS3A', 'MVS3A_rng', 'MVS4B', 'MVS4B_rng', 'MPS1A_HDSL', 'MPS1A_HDSL_rng', 'MPS1A_HDSR', 'MPS1A_HDSR_rng', 'MPS2B_HDSL', 'MPS2B_HDSL_rng', 'MPS2B_HDSR', 'MPS2B_HDSR_rng', 'MBS3A_topC', 'MBS3A_bottomC', 'MBS3A_leftC', 'MBS3A_rightC', 'MBS4B_topC', 'MBS4B_bottomC', 'MBS4B_leftC', 'MBS4B_rightC', 'MVS5A', 'MVS5A_rng', 'MVS6B', 'MVS6B_rng', 'XPS_i0', 'XPS_i0_rngN', 'XPS_tey', 'XPS_tey_rngN', 'XPS_fy', 'XPS_fy_rngN', 'RIXS_i0', 'RIXS_i0_rngN', 'RIXS_tey', 'RIXS_tey_rngN', 'RIXS_fy', 'RIXS_fy_rngN', 'RIXS_pd', 'RIXS_pd_rng', 'RIO_P_AvgTime', 'RIO_R_AvgTime', 'RIXS_X', 'RIXS_Y', 'RIXS_Z', 'RIXS_Ry', 'RIXS_SAMPLE_X', 'RIXS_SAMPLE_Z', 'RIXS_GRAS_Rx1', 'RIXS_GRAS_Rz1', 'RIXS_GRAS_X', 'RIXS_GRAS_Y', 'RIXS_GRAS_Z', 'RIXS_DETS_Y', 'RIXS_DETS_Z', 'XPS_X', 'XPS_Y', 'XPS_Z', 'XPS_Ry']
xas_attrs = ['scan', 'start_time', 'scan_type', 'exposure', 'temperatureA', 'temperatureB', 'samplex', 'sampley', 'samplez', 'th']

# %% new empty file 
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
    # _1D      = new_sheet(doc, 'ASCAN',    i=3)
    # _2D      = new_sheet(doc, 'MESH',     i=4)
    remove_sheet_by_name(doc, 'Sheet1')
    select_active_sheet(doc, settings)

    # unprotect
    settings.unprotect('123456')

    ##############
    # BACKGROUND #
    ##############
    cells = get_cells(settings, 'A1:F6')
    set_property_value(cells, 'CellBackColor', colors['light_grey_2'])

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
        text += '    RIXS folderpath:\n     path to folder with RIXS folder-scans\n\n'
        text += '    XAS folderpath:\n        path to folder with XAS scans\n\n'
        text += "    Save after update:\n        if `yes`, spreadsheet will saved after running the `update` macro\n\n" 
        text += "    Intercalate bkg color:\n        if `yes`, lines will have alternating bkg color\n\n" 
        text += "    Reset data:\n        if `yes`, all data (rows) is deleted and re-read again upon update\n\n" 
        text += '2) Run update macro:\n'
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
        _input(settings, 'RIXS folderpath',       'A1', 'B1:F1', '', True, input_width=4000)
        _input(settings, 'XAS folderpath',        'A2', 'B2:F2', '', True)
        _input(settings, 'Save after update',     'A4', 'B4', 'yes', True)
        _input(settings, 'Intercalate bkg color', 'A5', 'B5', 'yes', True)
        _input(settings, 'Reset data',            'D4', 'E4', 'no', True)

    ############
    # Internal #
    ############
    if True:
        _output(settings, 'RIXS: number of scans', 'I1', 'J1', 0)
        _output(settings, 'RIXS: current row',     'I2', 'J2', 2)

        _output(settings, 'XAS: number of scans', 'I3', 'J3', 0)
        _output(settings, 'XAS: current row',     'I4', 'J4', 2)

    ###########
    # protect #
    ###########
    settings.protect('123456')

    ##########
    # HEADER #
    ##########
    if True:
        # rixs
        final = ['comments', 'error'] + rixs_attrs + ['']
        set_row_data(rixs, 0, final)
        set_row_data(rixs, 1, final[2:-1], start_col=2)

        cells = get_cells(rixs, position=(0, 0, len(final)+2, 0))
        set_property_value(cells, 'CharWeight', 200)
        set_property_value(cells, 'CharHeight', 12)
        set_property_value(cells, 'CellBackColor', colors['light_grey_2'])

        cells = get_cells(rixs, position=(0, 1, len(final)+2, 1))
        set_property_value(cells, 'CellBackColor', colors['light_grey_2'])
        set_border(cells, 'bottom', OuterLineWidth=40, LineWidth=40)

        row = get_row(rixs, 0)
        set_width(row, 'optimal')
        for col in row.Columns:
            set_width(col, get_width(col)*1.1)

        # XAS
        final = ['comments', 'error'] + xas_attrs + ['']
        set_row_data(xas, 0, final)
        set_row_data(xas, 1, final[2:-1], start_col=2)

        cells = get_cells(xas, position=(0, 0, len(final)+2, 0))
        set_property_value(cells, 'CharWeight', 200)
        set_property_value(cells, 'CharHeight', 12)
        set_property_value(cells, 'CellBackColor', colors['light_grey_2'])

        cells = get_cells(xas, position=(0, 1, len(final)+2, 1))
        set_property_value(cells, 'CellBackColor', colors['light_grey_2'])
        set_border(cells, 'bottom', OuterLineWidth=40, LineWidth=40)

        row = get_row(xas, 0)
        set_width(row, 'optimal')
        for col in row.Columns:
            set_width(col, get_width(col)*1.1)

    ###############
    # unlock undo #
    ###############
    unlock_undo(doc)
    return

# %% update
def update(*args, **kwargs):
    """update macro for datafiles"""
    # get document object
    doc, settings, rixs, xas = _check_document()
    
    # get settings from settings sheet
    values = _check_settings(settings)

    #########
    # reset #
    #########
    if values['Reset data'] == 'yes':
        lock_undo(doc)
        if get_last_used_row(rixs) >= 2:
            delete_rows(sheet=rixs, start=2, stop=get_last_used_row(rixs))
        if get_last_used_row(xas) >= 2:
            delete_rows(sheet=xas,  start=2, stop=get_last_used_row(xas))
        settings.unprotect('123456')
        set_cells_data(settings, 'J1', 0)
        set_cells_data(settings, 'J2', 2)
        set_cells_data(settings, 'J3', 0)
        set_cells_data(settings, 'J4', 2)
        settings.protect('123456')
        set_cells_data(settings, 'E4', 'no')
        unlock_undo(doc)
        values = _check_settings(settings)

    ############
    # METADATA #
    ############
    for i, folderpath in enumerate((values['RIXS folderpath'], values['XAS folderpath'])):
        if folderpath != '':
            # scanlist
            scanlist = ipe.scanlist(folderpath=folderpath)

            # get sheet and next_row
            if i == 0:
                sheet    = rixs
                next_row = int(values['RIXS: next row'])
            else:
                sheet    = xas
                next_row = int(values['XAS: next row'])

            # heaser and last col
            header   = get_row_data(sheet, 1)
            last_col = get_last_used_col(sheet)

            # loop each scan and load metadata
            for _row, scan in enumerate(list(scanlist)[int(next_row - 2):]):
                row = next_row

                # color row (intercalated)
                if values['Intercalate bkg color']:
                    temp = sheet.getCellRangeByPosition(0, row, last_col, row)
                    lock_undo(doc)
                    if row%2 == 0:
                        temp.setPropertyValue('CellBackColor', colors['light_gold_4'])
                    else:
                        temp.setPropertyValue('CellBackColor', colors['light_indigo_4'])
                    unlock_undo(doc)

                # paste on spreadsheet
                try:
                    if i == 0:
                        pe1, pe2, pes1, pes2 = ipe.read(scanlist[scan], verbose=False)
                    elif i == 1:
                        pe1, TFY, I0, PD     = ipe.read(scanlist[scan], verbose=False)


                    # sort attrs
                    for col, attr in enumerate(header):
                        if attr != '':
                            if hasattr(pe1, attr):
                                value = getattr(pe1, attr)
                                if isinstance(value, datetime.datetime):
                                    value = str(value) 
                                elif isinstance(value, list) or isinstance(value, tuple):
                                    value = str(value)
                                elif br.is_number(value):
                                    value = float(value)
                                elif value is None:
                                    value = 'None'
                                # else:
                                #     value = str(value)
                                lock_undo(doc)
                                set_cells_data(sheet, (col, row), value)
                                unlock_undo(doc)
                            else:
                                if attr == 'error':
                                    value = ''
                                    lock_undo(doc)
                                    set_cells_data(sheet, (col, row), value)
                                    unlock_undo(doc)
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

# %% ---------------------- Import function as Macro ---------------------- %% #
# Only the specified functions will show in the Tools > Macro > Organize Macro dialog
g_exportedScripts = (update, new_empty, test1, test2)
# %%
