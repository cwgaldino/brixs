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
