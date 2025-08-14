"""Quality-of-life functions for libreoffice

Usage:
    >>> from brixs.sheets.libreoffice import create_instance, msgbox, colors
    >>> from brixs.sheets.libreoffice import new_calc, new_writer, new_impress, new_draw, new_math, new_base, open_document
    >>> from brixs.sheets.libreoffice import get_filepath, get_filename, save, close, lock_undo, unlock_undo
"""

# %% ------------------------- standard imports --------------------------- %% #
from pathlib import Path

# %% ---------------------------- uno imports ----------------------------- %% #
import uno
from com.sun.star.awt import MessageBoxButtons as MSG_BUTTONS
from com.sun.star.beans import PropertyValue

CTX = uno.getComponentContext()
SM  = CTX.getServiceManager()

# %% -------------------------- brixs imports ----------------------------- %% #
import brixs as br
# %%

# %% ----------------------------- Colors --------------------------------- %% #
colors = {}
colors['light_grey_2']   = 11711154
colors['light_gold_4']   = 16774606
colors['light_indigo_4'] = 14605542
colors['gold']           = 16750383

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
    """Save document (for now, it only works for calc)

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





