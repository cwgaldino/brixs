#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Everyday use functions that eases matplotlib figure manipulation."""


# %% ------------------------- Standard Imports --------------------------- %% #
import sys
import copy
from pathlib import Path
import numpy as np
import warnings
from subprocess import Popen, PIPE

from collections import OrderedDict
try:
    from bs4 import BeautifulSoup
except ModuleNotFoundError:
    pass

# %% ------------------------- Special Imports ---------------------------- %% #
import brixs as br

# %% ------------------------- Matplotlib Imports ------------------------- %% #
import matplotlib
from matplotlib.pyplot import get_current_fig_manager
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator

# %% ------------------------- backpack Imports --------------------------- %% #
from .arraymanip import index
from .interact import copy2clipboard, png2clipboard, svg2clipboard, operating_system
from .numanip import n_digits, n_decimal_places, round_to_1, is_integer

# %% ------------------------- Initial definitions ------------------------ %% #
is_windows = operating_system() == 'windows'
is_linux   = operating_system() == 'linux'
is_mac     = operating_system() == 'mac'

def available_colors():
    """Returns matplotlib available colors."""
    from matplotlib import colors as mcolors    
    return dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

def set_cycler(n='tab10', ls=1):
    """Set color cycles for matplotlib.
    
    Args:
        n (int, str, list, optional): name of a matplotlib sequential colormap.      
            If n is a number, n colors will be selected in order from the 
            listing: 'k', 'b', 'g', 'r', 'c', 'm', 'y', tab20b, tab20c. 
            The maximum value is then 47 colors. n can also be a list of colors 
            by name. Default is 'tab10'. Use None to use the default.
        ls (int, list, optional): number of different linestyles to iterate from
            (from 1 to 4) or a list of linestyles by name. Default is 1.

    Returns:
        None
    """
    from cycler import cycler

    # colors
    if n is None:
        n = 'tab10'
    elif type(n) == str:
        if n not in plt.colormaps():
            raise ValueError(f'{n} is not a recognized colormap. To check all available colormaps use plt.colormaps().')
        else:
            colors = plt.get_cmap(n).colors
    elif type(n) == int:
        if n > 47 or n < 1:
            raise ValueError('Permitted values for the number of colors is between 1 and 67.')
        init = ['k', 'b', 'g', 'r', 'c', 'm', 'y']
        init += plt.get_cmap('tab20b').colors + plt.get_cmap('tab20c').colors
        colors = init[0:n]
    elif type(n) == list:
        available_colors = available_colors()
        for color in n:
            if color not in available_colors:
                raise ValueError(f'{color} is not a recognized color. To check all available colors use brixs.available_colors().')
        colors = n
    else:
        raise ValueError(f'"n" must be False, int, str, or a list.\nType(n) = {type(n)}')

    # linestyle
    linestyles_available = ('-', '--', ':', '-.')
    if type(ls) == int:
        if ls > 4 or ls < 1:
            raise ValueError('Permitted values for the number of linestyles is between 1 and 4.')
        linestyles = linestyles_available[0:ls]
    elif type(ls) == list:
        linestyles = ls
    else:
        raise ValueError(f'"ls" must be False, int, or a list.\nType(ls) = {type(ls)}')
    
    # adjust lists
    number_of_colors     = len(colors)
    number_of_linestyles = len(linestyles)
    colors     = colors*number_of_linestyles
    linestyles = list(np.repeat(linestyles, number_of_colors))
    
    # set colormap and linestyles
    plt.rc('axes', prop_cycle=(cycler('color', colors) + cycler('linestyle', linestyles)))

def bring2top():
    """Brings current window to the top.

    This function was not tested for all available matplotlib backends.
    """
    # fig = plt.gcf()
    # fig.canvas.manager.window.raise_()

    backend = matplotlib.get_backend()
    if backend.startswith(('Qt5', 'qt5', 'QT5')):

        from PyQt5 import QtCore
        window = plt.get_current_fig_manager().window
        window.setWindowFlags(window.windowFlags() | QtCore.Qt.WindowStaysOnTopHint)
        plt.show()
        # window.setWindowFlags(window.windowFlags() & ~QtCore.Qt.WindowStaysOnTopHint)
        # plt.show()
    elif backend.startswith(('tk', 'TK', 'Tk')):
        plt.gcf().canvas.get_tk_widget().focus_force() 
    else:
        figManager = get_current_fig_manager()
        figManager.window.raise_()

def set_window_position(*args):
    """Change position of a maptplotlib figure on screen.

    Tipically, (0, 0) is the top left corner.

    Args:
        *args: A tuple like (x, y) or two separate x, y values (in px). If None,
            the default value defined by :py:func:`set_default_window_position`
            will be used. If there is no defined default value, the window's
            position will be preserved.

    See Also:
        :py:func:`get_window_position`
    """
    if len(args) > 1:
        x = int(args[1])
        y = int(args[0])
    elif len(args) == 1 and len(args[0]) == 2:
        x = int(args[0][1])
        y = int(args[0][0])
    elif len(args) == 0:
        if br.settings.FIGURE_POSITION is not None:
            position = (int(br.settings.FIGURE_POSITION[0]), int(br.settings.FIGURE_POSITION[1]))
            set_window_position(position)
        return
    else:
        warnings.warn('Wrong input')
        return

    figManager    = get_current_fig_manager()
    height, width = get_window_size()

    try:  # tested on tKinter backend
        figureGeometry = str(width) + 'x' + str(height) + '+' + str(x) + '+' + str(y)
        figManager.window.wm_geometry(figureGeometry)
    except AttributeError:
        try:  # tested on qt4 and qt5 backends
            figManager.window.setGeometry(int(x), int(y), width, height)
        except AttributeError:
            warnings.warn('Backend not suported.')
    
def get_window_position():
    """Return the position of a matplotlib position on the screen.

    Tipically, (0, 0) is the top left corner of your monitor.

    Returns:
        Tuple with the x and y position.

    See Also:
        :py:func:`set_window_position`
    """
    figManager = get_current_fig_manager()

    try:  # tested under tKinter backend
        return (figManager.window.winfo_y(), figManager.window.winfo_x())
    except AttributeError:  # tested under qt4 and qt5 backends
        try:
            return (figManager.window.geometry().y(), figManager.window.geometry().x())
        except AttributeError:
            warnings.warn('Backend not suported.')
            return (0, 0)

def set_window_size(*args):
    """Change the size of the window of a matplotlib figure.

    Args:
        *args: A tuple like (width, height) or two separate width, height values
            in px.

    See Also:
        :py:func:`get_window_size`
    """
    if len(args) > 1:
        width  = int(args[0])
        height = int(args[1])
    elif len(args) == 1 and len(args[0]) == 2:
        width  = int(args[0][0])
        height = int(args[0][1])
    else:
        return 
    
    figManager = get_current_fig_manager()
    x,y = get_window_position()

    try:  # tested on tKinter backend
        figureGeometry = str(height) + 'x' + str(width) + '+' + str(x) + '+' + str(y)
        figManager.window.wm_geometry(figureGeometry)
    except AttributeError:
        try:  # tested on qt4 and qt5 backends
            figManager.window.setGeometry(x, y, height, width)
        except AttributeError:
            warnings.warn('Backend not suported.')

    # This also works:
    # plt.gcf().set_size_inches(height, width)

def get_window_size():
    """Returns the size of the window of a matplotlib figure.

    Returns:
        Tuple with the width and height values in px.

    See Also:
        :py:func:`set_window_size`
    """
    figManager = get_current_fig_manager()

    try:  # tested on tKinter backend
        return (figManager.window.winfo_height(), figManager.window.winfo_width())

    except AttributeError:  # tested on qt4 and qt5 backends
        try:
            return (figManager.window.geometry().height(), figManager.window.geometry().width())
        except AttributeError:
            warnings.warn('Backend not suported.')
            return (0, 0)
    
    # this also works, but in inches
    # return list(plt.gcf().get_size_inches())

def maximize():
    """Maximize current fig."""
    figManager = plt.get_current_fig_manager()

    try:  # tested on tKinter backend
        figManager.frame.Maximize(True)

    except AttributeError:  # tested on qt4 and qt5 backends
        try:
            figManager.window.showMaximized()
        except AttributeError:
            warnings.warn('Backend not suported.')
            return (0, 0)

def _fix_figure(fig, **kwargs):
    cid1 = fig.canvas.mpl_connect('button_press_event', onclick)
    # cid2 = fig.canvas.mpl_connect('resize_event', onmove)

    # position
    if br.settings.FIGURE_POSITION is not None:
        set_window_position()
    if br.settings.FIGURE_FORCE_ON_TOP:
        bring2top()

    # size and dpi
    if 'dpi' not in kwargs:
        if br.settings.FIGURE_DPI is not None:
            # kwargs['dpi'] = br.settings.FIGURE_DPI
            fig.set_dpi(br.settings.FIGURE_DPI)
    if 'figsize' not in kwargs:
        if br.settings.FIGURE_SIZE is not None:
            # kwargs['figsize'] = br.settings.FIGURE_SIZE
            set_window_size(br.settings.FIGURE_SIZE)

        # grid
        if br.settings.FIGURE_GRID != (1, 1) and br.settings.FIGURE_GRID:
            rows    = br.settings.FIGURE_GRID[0]
            columns = br.settings.FIGURE_GRID[1]

            if (rows > 1 and columns > 0) or (columns > 1 and rows > 0):
                count  = br.settings._figure_count - 1
                row    = int((count/columns)%rows)
                column = count%columns

                if br.settings.FIGURE_SIZE is None:
                    height, width = get_window_size()
                else:
                    height = br.settings.FIGURE_SIZE[0]
                    width  = br.settings.FIGURE_SIZE[1]

                position = (br.settings.FIGURE_POSITION[0]+row*(height+br.settings.FIGURE_GRID_OFFSET[0]), br.settings.FIGURE_POSITION[1]+column*(width+br.settings.FIGURE_GRID_OFFSET[1]))
                set_window_position(position)

                br.settings._figure_count += 1
        else:
            set_window_position()
            br.settings._figure_count = 0
    return fig

def figure(**kwargs):
    """Create figure object.

    This command is the same as ``plt.figure()``, but is attaches the function
    :py:func:`onclick` to the figure so everytime you click on the figure, it
    calls :py:func:`onclick`.

    If figsize is passed as an argument, grid is disabled.

    Args:
        **kwargs: kwargs are passed to ``plt.figure()``.

    See Also:
        :py:func:`onclick`
    """   
    # open figure and initialize event callbacks
    fig = plt.figure(**kwargs)
    fig = _fix_figure(fig, **kwargs)
    return fig

def onmove(event):
    if br.br.settings.FIGURE_FIX_RESOLUTION:
        fig = plt.gcf()
        fig.set_dpi(fig.old_dpi)
        # time.sleep(1)

        # dpi
        

        # adjust size
        # set_window_size(fig.old_size)
        
        # # before
        # old_size = fig.old_size
        # old_dpi  = fig.old_dpi

        # # current
        # size = (event.width, event.height)
        # dpi  = old_dpi
        
        # if old_size == size:
        #     print('same')
        # else:
        #     print('changed')
        #     print(size, old_size)

        #     # new (corrected)
        #     new_size = size
        #     new_dpi  = dpi*old_size[0]/(size[0])
        #     print(new_dpi)

        #     # set
        #     # fig.set_dpi(new_dpi)
        #     fig.old_size = new_size
        #     # set_window_size(new_size)

        # # print('moved')
        # # print(new_dpi)
        # # # print('old size: ' + str(old_size))
        # # # print('new size: ' + str(new_size))
    pass

def set_onclick(format='svg', resolution=300, round_x=2, round_y=2, folder=None):
    """Set the default format for saving figures using :py:func:`onclick()`.

    Args:
        format (string, optional): 'svg' or 'png'.
        resolution (int, optional): resolution for png images in dpi.
        folder (string or pathlib.Path): folderpath to save figures.
    """

    global onclick_fig_format, onclick_resolution, onclick_folder, onclick_round_x, onclick_round_y
    if format == 'png':
        onclick_fig_format = 'png'
        onclick_resolution = resolution
    else:
        onclick_fig_format = format

    if folder is None:
        onclick_folder = Path.cwd()
    else:
        onclick_folder = Path(folder)

    if round_x is None:
        onclick_round_x = None
    else:
        onclick_round_x = round_x

    if round_y is None:
        onclick_round_y = None
    else:
        onclick_round_y = round_y

def onclick(event):
    """This function is called everytime a mouse key is pressed over a figure.

    Right click:
        x value is copied to the clipboard.
    Left click OR (y + Right click):
        y value is copied to the clipboard.
    Middle click:
        Figure is saved as svg or png in the default folder and copied to the clipboard
        (see :py:func:`set_onclick()`).

    Note:
        The matplotlib figure must be started by :py:func:`backpack.figmanip.figure`, and not
        by the default `matplotlib.pyplot.figure() <https://matplotlib.org/3.3.1/api/_as_gen/matplotlib.pyplot.figure.html>`_ fuction.

    Note:
        On linux this function uses xsel and xclip (use ``sudo apt-get install -y xsel`` and
        ``sudo apt-get install -y xclip`` to install them). The middle click
        function is not implemented on windows and mac yet.
    """

    # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
    #       ('double' if event.dblclick else 'single', event.button,
    #        event.x, event.y, event.xdata, event.ydata))
    try:
        onclick_round_y
    except NameError:
        onclick_round_y = 2

    try:
        onclick_round_x
    except NameError:
        onclick_round_x = 2

    if event.key == 'y' or event.button == 3:
        try:
            ax = plt.gca()
            ticks_showing = [y for y in ax.get_yticks() if y >= ax.get_ylim()[0] and y <= ax.get_ylim()[1]]
            delta = ticks_showing[-1] - ticks_showing[0]
            n = n_decimal_places(round_to_1(delta))*2
            y = round(event.ydata, n)
            # if onclick_round_y is not None:
            #     y = round(event.ydata, onclick_round_y)
            # else:
            #     y = event.ydata
            copy2clipboard(str(y))

        except TypeError:
            pass
    else: # left
        try:
            # if onclick_round_x is not None:
            #     x = round(event.xdata, onclick_round_x)
            # else:
            ax = plt.gca()
            ticks_showing = [x for x in ax.get_xticks() if x >= ax.get_xlim()[0] and x <= ax.get_xlim()[1]]
            delta = ticks_showing[-1] - ticks_showing[0]
            n = n_decimal_places(round_to_1(delta))*2
            x = round(event.xdata, n)
            # x = event.xdata
            # print(x)
            copy2clipboard(str(x))
        except TypeError:
            pass
    # double click (put image on clipboard)
    if event.dblclick or event.button == 2:
        global onclick_fig_format, onclick_resolution, onclick_folder

        try:
            onclick_fig_format
        except NameError:
            onclick_fig_format = 'svg'
        try:
            onclick_folder
        except NameError:
            onclick_folder = Path.cwd()
        try:
            onclick_resolution
        except NameError:
            onclick_resolution = 300

        onclick_folder = Path(onclick_folder)
        if onclick_fig_format == 'svg':
            if is_linux:
                plt.savefig(f'{onclick_folder/".temporary_fig.svg"}')
                svg2clipboard(onclick_folder/".temporary_fig.svg")
            elif is_windows:
                print('hh')

        elif onclick_fig_format == 'png':
            if is_linux:
                plt.savefig(f'{onclick_folder/".temporary_fig.png"}', dpi=onclick_resolution)
                png2clipboard(onclick_folder/".temporary_fig.png")

def zoom(start, stop, ymargin=2):
    """Zoom up portion of current figure from start to stop.

    Args:
        start (float or int): initial x value.
        ymargin (number, optional): margin value between data and the edges of
            plot in percentage of the y data range.
    """
    fig = plt.gcf()

    ymax = None
    ymin = None

    for axis in fig.axes:
        for line in axis.get_lines():

            x = line.get_data()[0]
            y = line.get_data()[1]

            _, y = br.extract(x=x, y=y, ranges=(start, stop))
            ymin_temp = min(y)
            ymax_temp = max(y)

            if ymin is None:
                    ymin = copy.copy(ymin_temp)
                    ymax = copy.copy(ymax_temp)
            try:
                if ymax_temp > ymax:
                    ymax = copy.copy(ymax_temp)
                if ymin_temp < ymin:
                    ymin = copy.copy(ymin_temp)
            except UnboundLocalError:
                warnings.warn("All data are outside of the required range. Cannot zoom.")
    m = (ymax-ymin)*ymargin/100
    plt.ylim(ymin-m, ymax+m)
    plt.xlim(start, stop)

def savefigs(filepath, figs='all'):
    """Save multiple matplotlib figures in a pdf.

    Args:
        filepath (string or pathlib.Path): filepath
        figs (list, optional): list with the figure numbers to save. Use 'all'
            to save all opened figures.
    """
    # check extension
    if Path(filepath).suffix == '.pdf':
        pass
    else:
        filepath = Path(filepath).with_suffix('.pdf')

    if figs == 'all':
        figs = [plt.figure(n) for n in plt.get_fignums()]

    if len(figs) > 1:
        pp = PdfPages(str(filepath))
        for fig in figs:
            fig.savefig(pp, format='pdf')
        pp.close()
    else:
        plt.savefig(str(filepath), format='pdf')

def cm2pt(*tupl):
    """Convert values from cm to pt.

    Args:
        *Args: single value or a tuple with values to convert

    Returns:
        A tuple with values converted
    """
    if isinstance(tupl[0], tuple):
        return tuple(i*28.346 for i in tupl[0])
    else:
        return tuple(i*28.346 for i in tupl)

def cm2px(*tupl, dpi=None):
    """Convert values from cm to px.

    Args:
        *Args: A tuple with values to convert
        dpi (int): pixel density or dots per inch.


    Returns:
        A tuple with values converted
    """
    if dpi is None:
        dpi = 100
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch*dpi for i in tupl[0])
    else:
        return tuple(i/inch*dpi for i in tupl)

def cm2inch(*tupl):
    """Convert values from cm to inches.

    Args:
        *Args: single value or a tuple with values to convert

    Returns:
        A tuple with values converted
    """
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def mm2inch(*tupl):
    """Convert values from mm to inches.

    Args:
        *Args: single value or a tuple with values to convert

    Returns:
        A tuple with values converted
    """
    inch = 2.54/10
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def set_ticks(ax=None, axis='x', autoscale=True, **kwargs):
    """Set ticks of a plot.

    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.
        axis (string, optional): possible values are 'x' or 'y'.
        **kwargs:

            #. min_value (float or int)
                start value for ticks (not the plot edge --- see ``pad``)
            #. max_value (float or int)
                stop value for ticks (not the plot edge --- see ``pad``)
            #. n_ticks (int)
                Number of ticks. Ticks separation
                is calculated accordingly and
                this parameter overwrites
                ticks_sep.
            #. ticks_sep (float or int)
                Ticks separation.
            #. n_minor_ticks (int)
                Number of minor ticks between two major ticks.
            #. fontproperties
                Label ticks font. Use ``matplotlib.font_manager.FontProperties``.
            #. pad (float or int)
                Distance between plot edge and the first tick in terms of tick separation.
                Typically, must be something between 0 and 1.
            #. n_decimal_places (int)
                Number of decimal places to use for tick labels.
            #. direction
                default is 'out', possible values are 'in' and 'out'

    Note:
        To set minor and major ticks 'manually' use `xaxis.set_ticks() <https://matplotlib.org/3.2.2/api/_as_gen/matplotlib.axis.XAxis.set_ticks.html>`_, for example::

            _ = ax.xaxis.set_ticks([1, 2, 3, 4, 5, 6], minor=True)

        Set labels manually by using:

            _ = ax.set_yticklabels([1, 2, 3, 4, 5, 6])

    """
    if ax is None:
        ax = plt.gca()
    use_sep = False

    if axis == 'y':
        ticks_showing = [y for y in ax.get_yticks() if y >= ax.get_ylim()[0] and y <= ax.get_ylim()[1]]
    elif axis == 'x':
        ticks_showing = [x for x in ax.get_xticks() if x >= ax.get_xlim()[0] and x <= ax.get_xlim()[1]]
    else:
        raise ValueError("Axis must be either 'x' or 'y'.")

    ## collecting kwargs
    if 'min_value' in kwargs:
        min_value = kwargs.pop('min_value')
        if min_value is None:
            min_value = ticks_showing[-1]
            for line in ax.get_lines():
                if axis == 'x':
                    if min_value > min(line.get_xdata()):
                        min_value = min(line.get_xdata())
                if axis == 'y':
                    if min_value > min(line.get_ydata()):
                        min_value = min(line.get_ydata())
            # min_value = ticks_showing[0]
    else:
        min_value = min(ticks_showing)
        for line in ax.get_lines():
            if axis == 'x':
                if min_value > min(line.get_xdata()):
                    min_value = min(line.get_xdata())
            if axis == 'y':
                if min_value > min(line.get_ydata()):
                    min_value = min(line.get_ydata())
        # min_value = ticks_showing[0]

    if 'max_value' in kwargs:
        max_value = kwargs.pop('max_value')
        if max_value is None:
            max_value = ticks_showing[0]
            for line in ax.get_lines():
                if axis == 'x':
                    if max_value < max(line.get_xdata()):
                        max_value = max(line.get_xdata())
                elif axis == 'y':
                    if max_value < max(line.get_ydata()):
                        max_value = max(line.get_ydata())
            # max_value = ticks_showing[-1]
    else:
        max_value = max(ticks_showing)
        for line in ax.get_lines():
            if axis == 'x':
                if max_value < max(line.get_xdata()):
                    max_value = max(line.get_xdata())
            elif axis == 'y':
                if max_value < max(line.get_ydata()):
                    max_value = max(line.get_ydata())
        # max_value = ticks_showing[-1]

    if 'n_ticks' in kwargs:
        n_ticks = kwargs.pop('n_ticks')
        if n_ticks is None:
            if 'ticks_sep' in kwargs:
                ticks_sep = kwargs.pop('ticks_sep')
                if ticks_sep is None:
                    n_ticks = len(ticks_showing)
                else:
                    use_sep = True
            else:
                n_ticks = len(ticks_showing)
    elif 'ticks_sep' in kwargs:
        ticks_sep = kwargs.pop('ticks_sep')
        if ticks_sep is not None:
            use_sep = True
        else:
            n_ticks = len(ticks_showing)
    else:
        use_sep = True
        ticks_sep = np.mean(np.diff(ticks_showing))

    if 'n_minor_ticks' in kwargs:
        n_minor_ticks = kwargs.pop('n_minor_ticks')
        if n_minor_ticks is None:
            n_minor_ticks = 2
    else:
        n_minor_ticks = 2

    if 'fontproperties' in kwargs:
        fontproperties = kwargs.pop('fontproperties')
    else:
        fontproperties = None

    if 'pad' in kwargs:
        pad = kwargs.pop('pad')
    else:
        pad = None
    
    

    # ticks
    if use_sep:
        # print('hh')
        ticks   = np.arange(min_value, max_value+ticks_sep*0.1, ticks_sep)
        # print(ticks_sep)
        # print(ticks)
    else:
        if n_ticks < 2:
            raise ValuError('n_ticks needs to be bigger or equal than 2.')
        ticks   = np.linspace(min_value, max_value, n_ticks)
    # ticks shift to get better values (include zero)
    if any(x<0 for x in ticks) and any(x>0 for x in ticks) and 0 not in ticks:
        ticks = ticks-ticks[index(ticks, 0)]
    elif max_value-min_value > 5:
        ticks = ticks-(ticks[0]-int(ticks[0]))
    # ticks edges
    # print(ticks[-1])
    # print(ticks[-1]+np.mean(np.diff(ticks))*0.5)
    # print(max_value)
    if ticks[-1]+np.mean(np.diff(ticks))*0.5 < max_value:
        # print('ff')
        ticks = np.append(ticks, ticks[-1]+np.mean(np.diff(ticks)))
    if ticks[0]-np.mean(np.diff(ticks))*0.8 > min_value:
        ticks = np.append(ticks[0]-np.mean(np.diff(ticks)), ticks)

    # limits
    if len(ticks) < 2:
        raise ValueError(f'ticks = {ticks} has only one tick, please, reduce ticks_sep.')
    try:
        if len(pad) == 2:
            min_lim = ticks[0] - (ticks[1]-ticks[0])*pad[0]
            max_lim = ticks[-1] + (ticks[1]-ticks[0])*pad[1]
        else:
            min_lim = ticks[0] - (ticks[1]-ticks[0])*pad[0]
            max_lim = ticks[-1] + (ticks[1]-ticks[0])*pad[0]
    except TypeError:
        if pad is not None:
            min_lim = ticks[0] - (ticks[1]-ticks[0])*pad
            max_lim = ticks[-1] + (ticks[1]-ticks[0])*pad


    # decimal places to show
    if 'n_decimal_places' in kwargs:
        n_decimal_places2 = kwargs.pop('n_decimal_places')
        if n_decimal_places2 is None:
            n_decimal_places2 = 0
            for n in ticks:
                if n_decimal_places(n)>n_decimal_places2:
                    n_decimal_places2 = n_decimal_places(n)
    else:
        n_decimal_places2 = 0
        for n in ticks:
            if n_decimal_places(n)>n_decimal_places2:
                n_decimal_places2 = n_decimal_places(n)

    if len(kwargs.keys()) > 0:
        raise AttributeError(f'Args not recognized: {list(kwargs.keys())}')



    # applying changes ======================
    s = '{' + f'0:.{n_decimal_places2}f' + '}'
    if axis == 'y':
        if fontproperties is None:
            dummy = ax.set_yticks(ticks)
        else:
            dummy = ax.set_yticks(ticks)
            dummy = ax.set_yticklabels([s.format(i) for i in ticks], fontproperties=fontproperties, visible=True)

        # minor ticks
        ax.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))


        # limits
        if pad is not None:
            ax.set_ylim((min_lim, max_lim), auto=False)

    elif axis == 'x':
        if fontproperties is None:
            dummy = ax.set_xticks(ticks)
        else:
            dummy = ax.set_xticks(ticks)
            dummy = ax.set_xticklabels([s.format(i) for i in ticks], fontproperties=fontproperties, visible=True)

        # minor ticks
        ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))

        # limits
        if pad is not None:
            ax.set_xlim((min_lim, max_lim), auto=False)

    # direction
    if 'direction' in kwargs:
        if kwargs['direction'] == 'in':
            direction = 'in'
        elif kwargs['direction'] == 'out':
            direction = 'out'
        else:
            raise ValueError('invalid value for direction.\nMust be either `in` or `out`.')
    else:
        direction = 'out'

    if axis == 'x':
        ax.tick_params(which='major', direction=direction, bottom=True, top=True, labeltop=False)
        ax.tick_params(which='minor', direction=direction, bottom=True, top=True, labeltop=False)
    elif axis == 'y':
        ax.tick_params(which='major', direction=direction, left=True, right=True, labelright=False, labeltop=False)
        ax.tick_params(which='minor', direction=direction, left=True, right=True)



    #
    # # test ==============================
    # if autoscale:
    #     if axis == 'x':
    #         string =''
    #         string +='def on_xlim_changed(ax):\n'
    #             # xlim = ax.get_xlim()
    #         string +='    min_value, max_value = ax.get_xlim()\n'
    #         string +=f'    n_ticks = {len(ticks)}\n'
    #
    #             # # ticks
    #             # ticks   = np.linspace(min_value, max_value, n_ticks)
    #
    #         string +='    print(n_ticks)\n'
    #
    #             # # ticks shift to get better values
    #             # if any(x<0 for x in ticks) and any(x>0 for x in ticks) and 0 not in ticks:
    #             #     ticks = ticks-ticks[index(ticks, 0)]
    #             # elif max_value-min_value > 5:
    #             #     ticks = ticks-(ticks[0]-int(ticks[0]))
    #             # # ticks edges
    #             # # print(ticks[-1])
    #             # # print(ticks[-1]+np.mean(np.diff(ticks))*0.5)
    #             # # print(max_value)
    #             # if ticks[-1]+np.mean(np.diff(ticks))*0.5 < max_value:
    #             #     # print('ff')
    #             #     ticks = np.append(ticks, ticks[-1]+np.mean(np.diff(ticks)))
    #             # if ticks[0]-np.mean(np.diff(ticks))*0.8 > min_value:
    #             #     ticks = np.append(ticks[0]-np.mean(np.diff(ticks)), ticks)
    #             #
    #             # if fontproperties is None:
    #             #     dummy = ax.set_xticks(ticks)
    #             # else:
    #             #     dummy = ax.set_xticks(ticks)
    #             #     dummy = ax.set_xticklabels([s.format(i) for i in ticks], fontproperties=fontproperties, visible=True)
    #             #
    #             # # minor ticks
    #             # ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
    #
    #             # if 'n_ticks' in kwargs:
    #             #     n_ticks = kwargs['n_ticks']
    #             #     if n_ticks is None:
    #             #         if 'ticks_sep' in kwargs:
    #             #             ticks_sep = kwargs['ticks_sep']
    #             #             if ticks_sep is None:
    #             #                 n_ticks = len(ticks_showing)
    #             #             else:
    #             #                 use_sep = True
    #             #         else:
    #             #             n_ticks = len(ticks_showing)
    #             # elif 'ticks_sep' in kwargs:
    #             #     ticks_sep = kwargs['ticks_sep']
    #             #     if ticks_sep is not None:
    #             #         use_sep = True
    #             #     else:
    #             #         n_ticks = len(ticks_showing)
    #             # else:
    #             #     use_sep = True
    #             #     ticks_sep = np.mean(np.diff(ticks_showing))
    #
    #             # for a in ax.figure.axes:
    #             #     # shortcuts: last avoids n**2 behavior when each axis fires event
    #             #     if a is ax or len(a.lines) == 0 or getattr(a, 'xlim', None) == xlim:
    #             #         continue
    #
    #                 # ylim = np.inf, -np.inf
    #                 # for l in a.lines:
    #                 #     x, y = l.get_data()
    #                 #     # faster, but assumes that x is sorted
    #                 #     start, stop = np.searchsorted(x, xlim)
    #                 #     yc = y[max(start-1,0):(stop+1)]
    #                 #     ylim = min(ylim[0], np.nanmin(yc)), max(ylim[1], np.nanmax(yc))
    #             #
    #             #     # TODO: update limits from Patches, Texts, Collections, ...
    #             #
    #             #     # x axis: emit=False avoids infinite loop
    #             #     a.set_xlim(xlim, emit=False)
    #             #
    #             #     # y axis: set dataLim, make sure that autoscale in 'y' is on
    #             #     corners = (xlim[0], ylim[0]), (xlim[1], ylim[1])
    #             #     a.dataLim.update_from_data_xy(corners, ignore=True, updatex=False)
    #             #     a.autoscale(enable=True, axis='y')
    #             #     # cache xlim to mark 'a' as treated
    #             #     a.xlim = xlim
    #         # print(string)
    #         exec(string)
    #         # print(on_xlim_changed)
    #         class Test:
    #             def __init__(self):
    #                 exec(string, globals())
    #                 # on_xlim_changed()
    #                 ax.callbacks.connect('xlim_changed', on_xlim_changed)

def remove_ticks_edge(ax):
    """Remove ticks that fall over the edges of the plot.

    This is useful when ticks are thicker than the plot edge lines.

    Args:
        ax (matplotlib.axes.Axes): axes instance.
    """
    ticks = ax.xaxis.get_major_ticks()

    if ax.get_xticks()[0] == ax.get_xlim()[0]:
        ticks[0].tick1line.set_visible(False)
        ticks[0].tick2line.set_visible(False)

    if ax.get_xticks()[-1] == ax.get_xlim()[-1]:
        ticks[-1].tick1line.set_visible(False)
        ticks[-1].tick2line.set_visible(False)

    ticks = ax.yaxis.get_major_ticks()

    if ax.get_yticks()[0] == ax.get_ylim()[0]:
        ticks[0].tick1line.set_visible(False)
        ticks[0].tick2line.set_visible(False)

    if ax.get_yticks()[-1] == ax.get_ylim()[-1]:
        ticks[-1].tick1line.set_visible(False)
        ticks[-1].tick2line.set_visible(False)

def vlines(x, ymin=None, ymax=None, colors=None, linestyles='solid', label='', ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()

    if ymin is None:
        ymin = ax.get_ylim()[0]
    if ymax is None:
        ymax = ax.get_ylim()[1]

    ax.vlines(x=x, ymin=ymin, ymax=ymax, colors=colors, linestyles=linestyles, label=label, **kwargs)

def hlines(y, xmin=None, xmax=None, colors=None, linestyles='solid', label='', ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()

    if xmin is None:
        xmin = ax.get_xlim()[0]
    if xmax is None:
        xmax = ax.get_xlim()[1]

    ax.hlines(y=y, xmin=xmin, xmax=xmax, colors=colors, linestyles=linestyles, label=label, **kwargs)

import types
def gridspec(nrows, ncols, instructions=None):
    """instructions = row_init, row_final, col_init, col_final, hspace, wspace"""
    fig = plt.gcf()

    if instructions is None:
        instructions = []
        for row in range(nrows):
            for col in range(ncols):
                instructions.append((row, row, col, col, None, None))

    axes = []
    for instruc in instructions:
        row_init  = instruc[0]
        row_final = instruc[1]
        col_init  = instruc[2]
        col_final = instruc[3]
        hspace    = instruc[4]
        wspace    = instruc[5]

        if row_final is None:
            row_final = nrows-1
        if col_final is None:
            col_final = ncols-1

        if hspace is None:
            hspace = .3
        if wspace is None:
            wspace = .2

        gs = fig.add_gridspec(nrows, ncols, hspace=hspace, wspace=wspace)
        axes.append(fig.add_subplot(gs[row_init:row_final+1, col_init:col_final+1]))


        def _remove_xticklabels(self):
            self.tick_params(axis='x', labelbottom=False)
        def _remove_yticklabels(self):
            self.tick_params(axis='y', labelbottom=False)
        def _remove_xticks(self, top=True, bottom=True):
            if bottom: bottom = False
            if top:    top    = False
            self.tick_params(axis='x', which='major', top=top, bottom=bottom)
            self.tick_params(axis='x', which='minor', top=top, bottom=bottom)
        def _remove_yticks(self, left=True, right=True):
            if left:  left  = False
            if right: right = False
            self.tick_params(axis='y', which='major', left=left, right=right)
            self.tick_params(axis='y', which='minor', left=left, right=right)
        def _remove_spline(self, left=False, right=False, top=False, bottom=False):
            if left:
                self.spines['left'].set_visible(False)
            if right:
                self.spines['right'].set_visible(False)
            if top:
                self.spines['top'].set_visible(False)
            if bottom:
                self.spines['bottom'].set_visible(False)
        for ax in axes:
            ax.remove_xticklabels = types.MethodType(_remove_xticklabels, ax)
            ax.remove_yticklabels = types.MethodType(_remove_yticklabels, ax)
            ax.remove_xticks = types.MethodType(_remove_xticks, ax)
            ax.remove_yticks = types.MethodType(_remove_yticks, ax)
            ax.remove_spline = types.MethodType(_remove_spline, ax)

    return axes

def subplots(nrows, ncols, sharex=False, sharey=False, hspace=None, wspace=None, width_ratios=None, height_ratios=None, gridspec_kw=None, **fig_kw):

    class myarray(list):
        def remove_label_from_shared_axis(self, sharex=None, sharey=None):
            nrows = axes.nrows
            ncols = axes.ncols
            if sharex is None:
                sharex = axes.sharex
            if sharey is None:
                sharey = axes.sharey
            if sharey:
                keep = np.arange(0, nrows*ncols, ncols)
                print(keep)
                for i in range(nrows*ncols):
                    if i not in keep:
                        self[i].yaxis.label.set_visible(False)
            if sharex:
                for i in range(0, (nrows-1)*ncols):
                    self[i].xaxis.label.set_visible(False)
    
    axes = myarray([0]*(ncols*nrows))
    axes.nrows  = nrows
    axes.ncols  = ncols
    axes.sharex = sharex
    axes.sharey = sharey

    # ratios
    if width_ratios is None:
        width_ratios=[1]*ncols
    if height_ratios is None:
        height_ratios=[1]*nrows
    
    if hspace is None:
        hspace = .3
    if wspace is None:
        wspace = .3

    if sharex:
        hspace = 0
    if sharey:
        wspace = 0

    # gridspec_kw
    if gridspec_kw is None:
        gridspec_kw = {}
    if 'wspace' not in gridspec_kw:
        gridspec_kw['wspace'] = wspace 
    if 'hspace' not in gridspec_kw:
        gridspec_kw['hspace'] = hspace 
    if 'width_ratios' not in gridspec_kw:
        gridspec_kw['width_ratios'] = width_ratios 
    if 'height_ratios' not in gridspec_kw:
        gridspec_kw['height_ratios'] = height_ratios 
    if 'wspace' not in gridspec_kw:
        gridspec_kw['wspace'] = hspace 

    fig, _axes = plt.subplots(nrows, ncols, 
                              sharex=sharex,
                              sharey=sharey,
                              gridspec_kw=gridspec_kw, **fig_kw)
    fig = _fix_figure(fig, **fig_kw)

    del gridspec_kw

    try:
        for i, ax in enumerate(br.flatten(_axes)):
            axes[i] = ax
    except TypeError:
        axes[0] = _axes

    return fig, axes

def ax_box2fig_box(ax, points):
    """Transform 'bbox like' axis position values to percentage fig position.

    Useful for positioning insets.

    Args:
        ax (matplotlib.axes.Axes): axes instance.
        points (list): list like ``[x_init, y_init, x_final, y_final]``

    Returns:
        'bbox like' figure positions ``[x_init, y_init, delta_x, delta_y]``.
    """
    [x_init, y_init, x_final, y_final] = points

    x_init  = ax_pos2fig_pos(ax, x_init, direction='x')
    y_init  = ax_pos2fig_pos(ax, y_init, direction='y')
    delta_x = ax_pos2fig_pos(ax, x_final, direction='x') - x_init
    delta_y = ax_pos2fig_pos(ax, y_final, direction='y') - y_init

    return [x_init, y_init, delta_x, delta_y]

def ax_pos2fig_pos(ax, value, direction='x'):
    """Transform axis position values to figure percentage position.

    Args:
        ax (matplotlib.axes.Axes): axes instance.
        value (float): value
        direction (string, optional): 'x' or 'y'.

    Returns:
        Figure positions from 0 to 1.
    """
    if direction == 'x':
        point1 = (ax.get_xlim()[0], ax.get_position().xmin)
        point2 = (ax.get_xlim()[1], ax.get_position().xmax)
    else:
        point1 = (ax.get_ylim()[0], ax.get_position().ymin)
        point2 = (ax.get_ylim()[1], ax.get_position().ymax)
    delta = (point2[1]-point1[1])/(point2[0]-point1[0])
    x0 = point2[1] - (delta*point2[0])

    return x0 + delta*value

def ungroup_svg(filepath, outfilepath=None):
    """Ungroup all objects in svg files created by matplotlib.

    Text with latex equations are kept grouped.

    Ticks are cloned. The parent ticks will be at the left upper corner of the figure.

    Args:
        filepath (string or pathlib.Path): path to svg file.
        outfilepath (string or pathlib.Path, optional): path for svg output file.
            If None, it will use the input filepath and the file will be overwrited.

    See Also:
        :py:func:`soft_ungroup`

    Notes:
        Use :py:func:`soft_ungroup` to better preserve the structure
        of the svg file in case :py:func:`ungroup` returns strange results.
    """
    # get soup
    filepath = Path(filepath)
    f = filepath.open()
    soup = BeautifulSoup(f, features = 'xml')
    f.close()

    # file prefix
    script_tags = soup.find_all('svg')
    prefix = script_tags[0].prettify().split('<g')[0]

    # end of file
    endOfFile = script_tags[0].prettify().split('</g>')[-1]

    # find objs
    script_tags = soup.find_all(['use', 'path', 'def', 'text'])
    objs = OrderedDict()  # dict.key() is the layer label, dict.value() is the layer itself
    for i in range(len(script_tags)):

        if script_tags[i].name == 'text':  # keep text with latex equations grouped
            if 'transform' in script_tags[i].parent.attrs:
                objs[i] = str(script_tags[i].parent.prettify())
            else:
                objs[i] = str(script_tags[i].prettify())
        else:
            # include obj
            objs[i] = str(script_tags[i].prettify())

    # prepare output
    output = copy.copy(prefix)
    for id in objs:
        output += objs[id]
    output += endOfFile

    # save
    if outfilepath is None:
        outfilepath = filepath
    else:
        outfilepath = Path(outfilepath)
    f = outfilepath.open('w')
    f.write(output)
    f.close()

def soft_ungroup_svg(filepath, outfilepath=None):
    """Ungroup objects in svg file (each object stands in a group by itself).

    Args:
        filepath (string or pathlib.Path): path to svg file.
        outfilepath (string or pathlib.Path, optional): path for svg output file.
            If None, it will use the input filepath and the file will be overwrited.

    See Also:
        :py:func:`ungroup`

    Notes:
        This differs from :py:func:`ungroup` as it better preserves the structure
        of the svg file, but Objects are still kept in groups.
    """
    # get soup
    filepath = Path(filepath)
    f = filepath.open()
    soup = BeautifulSoup(f, features = 'xml')
    f.close()

    # file prefix
    script_tags = soup.find_all('svg')
    prefix = script_tags[0].prettify().split('<g')[0]

    # end of file
    endOfFile = script_tags[0].prettify().split('</g>')[-1]

    # find objs
    script_tags = soup.find_all('g')
    objs = OrderedDict()  # dict.key() is the layer label, dict.value() is the layer itself
    for i in range(len(script_tags)):
        if 'id' in script_tags[i].attrs:

            if script_tags[i].attrs['id'].startswith(('figure', 'axes', 'xtick', 'matplotlib.axis', 'xtick', 'ytick')):
                pass
            else:
                # print(script_tags[i].attrs['id'])
                del_list = []
                for j in range(len(script_tags[i].contents)):
                    try:
                        if 'id' in script_tags[i].contents[j].attrs:
                            del_list.append(j)
                    except AttributeError:
                        pass
                del_list = [x-n for n,x in enumerate(del_list)]
                for j in del_list:
                    del script_tags[i].contents[j]

                # include obj
                objs[script_tags[i].attrs['id']] = str(script_tags[i])#.prettify()

    # prepare output
    output = copy.copy(prefix)
    for id in objs:
        output += objs[id]
    output += endOfFile

    # save
    if outfilepath is None:
        outfilepath = filepath
    else:
        outfilepath = Path(outfilepath)
    f = outfilepath.open('w')
    f.write(output)
    f.close()


# gradient functions ========================================================
def rgb2hex(rgb, max_rgb_value=1):
    """Return hex value.

    Args:
        rgb (tuple): rgb values ([255 255 255])
        max_rgb_value (int, optional): max rgb possible value. Tipically, 1 or 255.

    Returns:
        hex value."""
    if max(rgb) > max_rgb_value:
        raise ValueError(f'rgb={rgb} maximum value is bigger than max_rgb_value.')
    return mpl.colors.to_hex([x/max_rgb_value for x in rgb])

def hex2rgb(string, max_rgb_value=1):
    """Return rgb value.

    Args:
        string (string): hex value as string.
        max_rgb_value (int, optional): max rgb possible value. Tipically, 1 or 255.

    Returns:
        rgb value (tuple)."""
    return [x*max_rgb_value for x in mpl.colors.to_rgb(string)]

def _linear_gradient(start, end, n=10, max_rgb_value=1):
    """Returns a gradient list of n colors between two hex colors.

    Args:
        start (tuple): rgb start value.
        end (tuple): rgb final value.
        n (int, optional): number of intermediate colors.
        max_rgb_value (int, optional): max rgb possible value. Tipically, 1 or 255.

    Note:
        see https://bsouthga.dev/posts/color-gradients-with-python.

    Returns:
        list with rgb colors.
    """
    # Starting and ending colors in RGB form
    start = hex2rgb(start, 255)
    end = hex2rgb(end, 255)
    # Initilize a list of the output colors with the starting color
    RGB_list = [start]
    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
        temp = [int(start[j] + (float(t)/(n-1))*(end[j]-start[j]))/255*max_rgb_value for j in range(3)]
        # Add it to our list of output colors
        RGB_list.append(temp)
    return RGB_list

def linear_gradient(colors, n, max_rgb_value=1):
    """Returns a list of colors forming linear gradients between colors.

    Note:
        see https://bsouthga.dev/posts/color-gradients-with-python

    Args:
      colors (tuple): list with rgb values.
      n (int, optional): number of intermediate colors.
      max_rgb_value (int, optional): max rgb possible value. Tipically, 1 or 255.

    Returns:
      list with rgb colors.
    """
    # The number of colors per individual linear gradient
    n_out = int(float(n) / (len(colors) - 1))
    # returns dictionary defined by color_dict()
    RGB_list = _linear_gradient(colors[0], colors[1], n_out, max_rgb_value=255)

    if len(colors) > 2:
        for col in range(1, len(colors) - 1):
            next = _linear_gradient(colors[col], colors[col+1], n_out, max_rgb_value=255)
        for color in next[1:]:
            RGB_list.append(color)
    return [[x/255*max_rgb_value for x in color] for color in RGB_list]

# Value cache
fact_cache = {}
def _factorial(n):
    """Memoized factorial function used by the  _bernstein() function."""
    try:
        return fact_cache[n]
    except(KeyError):
        if n == 1 or n == 0:
            result = 1
        else:
            result = n*_factorial(n-1)
        fact_cache[n] = result
        return result

def _bernstein(t, n, i):
    """Bernstein coefficient used by the bezier_gradient() function."""
    binom = _factorial(n)/float(_factorial(i)*_factorial(n - i))
    return binom*((1-t)**(n-i))*(t**i)

def bezier_gradient(colors, n=100, max_rgb_value=1):
    """Returns a list of rgb colors forming a bezier gradient between colors.

    Note:
        see https://bsouthga.dev/posts/color-gradients-with-python.

    Args:
      colors (tuple): list with hex values.
      n (int, optional): number of intermediate colors.
      max_rgb_value (int, optional): max rgb possible value. Tipically, 1 or 255.

    Returns:
      list with rgb colors.
    """
    # RGB vectors for each color, use as control points
    RGB_list = [hex2rgb(color, 255) for color in colors]
    degree = len(RGB_list) - 1

    def bezier_interp(t):
        """Define an interpolation function for this specific curve."""
        # List of all summands
        summands = [list(map(lambda x: int(_bernstein(t, degree, i)*x), c)) for i, c in enumerate(RGB_list)]
        # Output color
        out = [0, 0, 0]
        # Add components of each summand together
        for vector in summands:
            for c in range(3):
                out[c] += vector[c]
        return out

    return [[x/255*max_rgb_value for x in bezier_interp(float(t)/(n-1))] for t in range(n)]
