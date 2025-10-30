#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for everyday use ---> Matplotlib figures"""

# %% ------------------------- Standard Imports --------------------------- %% #
from string import ascii_lowercase
from types import MethodType
from pathlib import Path
import numpy as np
import warnings
import decimal
import copy

# %% ------------------------- Matplotlib Imports ------------------------- %% #
import matplotlib
from matplotlib.pyplot import get_current_fig_manager as _get_current_fig_manager
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
from matplotlib import colormaps
from collections.abc import Iterable
from cycler import cycler

# %% -------------------------- operating system -------------------------- %% #
import platform
def _operating_system():
    """Return string with name of operating system (windows, linux, or mac)."""
    system = platform.system().lower()
    is_windows = system == 'windows'
    is_linux = system == 'linux'
    is_mac = system == 'darwin'
    if is_windows:
        return 'windows'
    elif is_linux:
        return 'linux'
    elif is_mac:
        return 'mac'
    else:
        raise ValueError('OS not recognized')

is_windows = _operating_system() == 'windows'
is_linux   = _operating_system() == 'linux'
is_mac     = _operating_system() == 'mac'

# %% --------------- supporting functions from numanip -------------------- %% #
# backpack developers note --> if these function change, it needs to be copied to numanip.py
import numbers
def is_number(n):
    """Returns True if variable is number."""
    # backpack developers note!!!!
    # if this function changes, it needs to be copied to these files: numanip
    if isinstance(n, bool):
        return False
    elif isinstance(n, int):
        return True
    elif isinstance(n, float):
        return True
    elif isinstance(n, numbers.Number):
        return True
    elif isinstance(n, str):
        try: 
            is_number(float(n))
            return True
        except ValueError:
            return False
    else:
        False
        
def _round_to_1(x):
    """return the most significant digit"""
    return round(x, -int(np.floor(np.log10(abs(x)))))

def _n_decimal_places(number, count_zero=False):
    """Return the number of decimal places of a number.

    Args:
        number (float or int): number.
        count_zero (bool, optional): if an integer is type float, it will came
            out with a zero after the decimal point, eg, `145.0` instead of `145`.
            count_zero == False will not count zeros after the decimal point.
            Default is False.

    Returns:
        number of decimal places in number.
    """
    n = abs(decimal.Decimal(str(number)).as_tuple().exponent)
    if count_zero == False and n == 1:
        if str(number)[-1] == '0':
            return 0
    return n

# %% ------------ supporting functions from arraymanip -------------------- %% #
# backpack developers note --> if these function change, it needs to be copied to arraymanip.py
def _index(x, value, closest=True, roundup=False):
    """Returns the first index of the element in array.

    Args:
        x (list or array): 1D array.
        value (float or int): value.
        closest (book, optional): if True, returns the index of the element in 
            array which is closest to value.
        roundup (bool, optional): if closest=True, and value is exactly midway
            between 2 items in array x, rounup=True will return the index of 
            item in x with highest value. Default is False.

    Returns:
        index (int)
    """
    # backpack developers note!!!!
    # if this function changes, it needs to be copied to these files: figmanip
    if closest:
        _inner1 = np.array(x) - value
        _inner2 = np.ma.masked_array(_inner1, np.isnan(_inner1))
        absv    = np.abs(_inner2)
        vmin    = np.min(absv)
        if np.sum(np.where(absv==vmin, 1, 0)) > 1:
            indexes = [_[0] for _ in np.argwhere(absv==vmin)]
            if roundup:
                if x[indexes[0]] > x[indexes[1]]:
                    return indexes[0]
                else:
                    return indexes[1]
            else:
                if x[indexes[0]] < x[indexes[1]]:
                    return indexes[0]
                else:
                    return indexes[1]
        else:
            return int(np.argmin(np.abs(_inner2)))
    else:
        return np.where(x == value)[0]


def _choose(x, ranges):
    """Return a mask of x values inside range pairs.

    Args:
        x (list or array): 1d array.
        ranges (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x. [[xi, xf], [xi_2, xf_2], ...]

    Returns:
        1d list.
    """
    assert isinstance(x, Iterable), 'input must be a iterable'
    x = np.array(x)

    try:  # ((start, end), )
        choose_range = [None]*len(ranges)
        for i, (x_init, x_final) in enumerate(ranges):
            choose_range[i] = np.logical_and(x>=x_init, x<=x_final)
        choose_range = [True if x == 1 else False for x in np.sum(choose_range, axis=0)]
    except TypeError:  # (start, end)
        x_init, x_final = ranges
        choose_range = np.logical_and(x>=x_init, x<=x_final)
    return choose_range

def _extract(x, y, ranges, invert=False):
    """Returns specific data ranges from x and y.

    Args:
        x (list or array): 1D reference vector.
        y (list or array): 1D y-coordinates or list of several data sets.
        ranges (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x.
        invert (bool, optional): if inverted is True, data outside of the data 
            will be returned. Default is False.

    Returns:
        x and y arrays. If `y` is 1d, the returned `y` is 1d. If `y` is
        a multicolumn array then the returned `y` is also multicolumn


    Examples:

        if `y` is 1d, the returned `y` is 1d:

        >>> x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        >>> y = np.array(x)**2
        >>> ranges = ((0, 3), (7.5, 9))
        >>> x_sliced, y_sliced = am.extract(x, y, ranges)
        >>> print(x_sliced)
        [0 1 2 3 8 9]
        >>> print(y_sliced)
        [0 1 4 9 64 81]

        if `y` is multicolumn, the returned `y` is also multicolumn:

        >>> x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        >>> y = np.zeros((10, 2))
        >>> y[:, 0] = x**2
        >>> y[:, 0] = x**3
        >>> ranges = ((0, 3), (7.5, 9))
        >>> x_sliced, y_sliced = am.extract(x, y, ranges)
        >>> print(x_sliced)
        [0. 1. 2. 3. 8. 9.]
        >>> print(y_sliced)
        [[  0.   0.]
         [  1.   0.]
         [  8.   0.]
         [ 27.   0.]
         [512.   0.]
         [729.   0.]]


    """
    x = np.array(x)
    y = np.array(y)

    # if data is all inside ranges, then nothing is done
    if len(ranges) == 1:
        if ranges[0][0] <= min(x) and ranges[0][1] >= max(x):
            if invert:
                return np.array([]), np.array([])
            else:
                return x, y

    choose_range = _choose(x, ranges)
    # print(choose_range[0:10])
    if invert:
        choose_range = np.invert(choose_range)
    # temp = np.compress(choose_range, np.c_[y.transpose(), x], axis=0)
    # print(choose_range[0:10])
    temp = np.compress(choose_range, np.c_[y, x], axis=0)
    # print(temp)
    if temp.any():
        if len(temp[0]) > 2:
            return temp[:, -1], temp[:, :-1]#.transpose()
        else:
            # print('here')
            return temp[:, -1], temp[:, 0]
    else:
        raise RuntimeError('No data points within the selected range.')

# %% ---------------- supporting functions from query --------------------- %% #
import subprocess
def _copy2clipboard(txt):
    """Copy text to clipboard.

    on linux it uses ``xclip`` package (``sudo apt install xclip``).
    """
    if is_windows:
        try:
            # cmd='echo ' + txt.strip() + ' | clip'
            cmd=f'echo|set /p={txt.strip()}| clip'
            # cmd='echo ' + txt.strip() + '| Set-Clipboard -Value {$_.Trim()}'
            subprocess.check_call(cmd, shell=True)
        except:
            pass
    elif is_linux:
        try:
            # p = subprocess.Popen(['xsel','-bi'], stdin=subprocess.PIPE)
            # p.communicate(input=bytes(txt.strip()).encode())
            # p.communicate(input=bytes(txt.strip(), encoding='utf-8'))
            # p.communicate(input=txt.strip())
            with subprocess.Popen(['xclip','-selection', 'clipboard'], stdin=subprocess.PIPE) as pipe:
                pipe.communicate(input=txt.strip().encode('utf-8'))
        except:
            pass
    elif is_mac:
        try:
            cmd='echo '+ txt.strip() + ' | pbcopy'
            subprocess.check_call(cmd, shell=True)
        except:
            pass
    return
# %%

# %% ================================ colors ============================== %% #
def get_available_colors():
    """Returns matplotlib available colors."""
    from matplotlib import colors as mcolors    
    return dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

def reset_cycler():
    """Resets the linestyle cycle"""
    plt.gca().set_prop_cycle(None)

def set_cycler(n='tab10', ls=1):
    """Set color cycles for matplotlib.
    
    Args:
        n (int, str, list, optional): name of a matplotlib sequential colormap.      
            If n is a number, n colors will be selected in order from the 
            listing: 'k', 'r', 'g', 'b', 'c', 'm', 'y', tab20b, tab20c. 
            The maximum value is then 47 colors. n can also be a list of colors 
            by name or by (rgb alpha). Default is 'tab10'. Use None to use the default.
        ls (int, list, optional): number of different linestyles to iterate from
            (from 1 to 4) or a list of linestyles by name. Default is 1.

    Returns:
        None
    """
    # reset cycler
    reset_cycler()

    # set default
    if n is None:
        n = 'tab10'
    
    # colors
    if type(n) == str:
        if n not in plt.colormaps():
            raise ValueError(f'{n} is not a recognized colormap. To check all available colormaps use plt.colormaps().')
        else:
            colors = plt.get_cmap(n).colors
    elif type(n) == int:
        if n > 47 or n < 1:
            raise ValueError('Permitted values for the number of colors is between 1 and 67.')
        init = ['k', 'r', 'g', 'b', 'c', 'm', 'y']
        init += plt.get_cmap('tab20b').colors + plt.get_cmap('tab20c').colors
        colors = init[0:n]
    elif type(n) == list:
        available_colors = get_available_colors()
        for color in n:
            if isinstance(color, str):
                assert color in available_colors, f'{color} is not a recognized color\nAvailable colors: {available_colors}'
            else:
                assert len(color) == 3 or len(color) == 4, f'{color} is not recognized\ncolors must be rgb-alpha values or color names'        
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
    colors               = colors*number_of_linestyles
    linestyles           = list(np.repeat(linestyles, number_of_colors))
    
    # set colormap and linestyles
    plt.rc('axes', prop_cycle=(cycler('color', colors) + cycler('linestyle', linestyles)))

    # # reset cycler
    # reset_cycler()

def get_available_colormaps():
    """returns a list of available colormaps"""
    return list(colormaps)

def set_cycler_from_colormap(name, n):
    """Sets a color cycler based on a colormap. 

    Args:
        name (str): name of the colormap
        n (int): number of colors to pick from the colormap.

    Returns:
        None
    """
    # available_colormaps = get_available_colormaps()
    # assert name in available_colormaps, f'{name} is not a recognized colormap\nAvailable colormaps: {available_colormaps}'

    # cmap = matplotlib.cm.get_cmap(name)
    # colors = [cmap(j/n) for j in range(n)]
    colors = get_colors_from_colormap(name, n)
    set_cycler(n=colors, ls=1)

def get_colors_from_colormap(name, n):
    """Returns a list of colors from colormap.

    Args:
        name (str): name of the colormap
        n (int): number of colors to pick from the colormap.

    Returns:
        None
    """
    available_colormaps = get_available_colormaps()
    assert name in available_colormaps, f'{name} is not a recognized colormap\nAvailable colormaps: {available_colormaps}'

    cmap   = matplotlib.cm.get_cmap(name)
    colors = [cmap(j/n) for j in range(n)]
    return colors
# %%

# %% ================================ window ============================== %% #
def set_window_position(*args):
    """Change position of a maptplotlib figure on screen.

    Typically, (0, 0) is the top left corner of the display.

    Usage:
        set_window_position((x, y))
        set_window_position(x, y)


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
    else:
        warnings.warn('Wrong input')
        return

    figManager    = _get_current_fig_manager()
    height, width = get_window_size()

    try:  # tested on tKinter backend
        figureGeometry = str(width) + 'x' + str(height) + '+' + str(x) + '+' + str(y)
        figManager.window.wm_geometry(figureGeometry)
    except AttributeError:
        try:  # tested on qt4 and qt5 backends
            figManager.window.setGeometry(int(x), int(y), width, height)
        except AttributeError:
            warnings.warn('Backend not supported.')
    return

def get_window_position():
    """Return the position of a matplotlib position on the screen.

    Typically, (0, 0) is the top left corner of your monitor.

    Returns:
        Tuple with the x and y position.

    See Also:
        :py:func:`set_window_position`
    """
    figManager = _get_current_fig_manager()

    try:  # tested under tKinter backend
        return (figManager.window.winfo_y(), figManager.window.winfo_x())
    except AttributeError:  # tested under qt4 and qt5 backends
        try:
            return (figManager.window.geometry().y(), figManager.window.geometry().x())
        except AttributeError:
            warnings.warn('Backend not suported.')
            return (0, 0)
    return

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
    
    figManager = _get_current_fig_manager()
    x, y       = get_window_position()

    try:  # tested on tKinter backend
        figureGeometry = str(height) + 'x' + str(width) + '+' + str(x) + '+' + str(y)
        figManager.window.wm_geometry(figureGeometry)
    except AttributeError:
        try:  # tested on qt4 and qt5 backends
            # figManager.window.setGeometry(x, y, height, width)
            figManager.window.setGeometry(y, x, height, width)
        except AttributeError:
            warnings.warn('Backend not supported.')

    # This also works:
    # plt.gcf().set_size_inches(height, width)
    return

def get_window_size():
    """Returns the size of the window of a matplotlib figure.

    Returns:
        Tuple with the width and height values in px.

    See Also:
        :py:func:`set_window_size`
    """
    figManager = _get_current_fig_manager()

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
    return

def get_window_dpi():
    """returns figure dpi"""
    return plt.gcf().dpi

def maximize():
    """Maximize current fig."""
    figManager = _get_current_fig_manager()

    try:  # tested on tKinter backend
        figManager.frame.Maximize(True)

    except AttributeError:  # tested on qt4 and qt5 backends
        try:
            figManager.window.showMaximized()
        except AttributeError:
            warnings.warn('Backend not suported.')
            return (0, 0)
    return
# %%
        
# %% ================================ figure ============================== %% #
def figure(*args, **kwargs):
    """Create figure object. Wrapper for `plt.figure()`_.

    `figsize` option must be given in cm.

    br.figure comes equipped with mouse clicking functionality

        Right click:
            x value is copied to the clipboard.
        Left click OR (y + Right click):
            y value is copied to the clipboard.
        Middle click:
            copies cursor position in terms of figure coordinates.
    
    Also, figure created by br.figure() comes with extra methods:

        >>> fig = br.figure()
        >>> fig.grid() # creates a figure grid with grid lines to help axes alignment

    Args:
        *args, **kwargs: args and kwargs are passed to `plt.figure()`.

    Note:
        This function overwrites the behavior of `figsize` parameters. In
        plt.figure(figsize=(w, h)), w and h must be given in inches. However,
        this function gets `w` and `h` in cm. 
    
    Returns:
        figure object
    
    .. _plt.figure(): https://matplotlib.org/stable/api/figure_api.html
    """   
    # figsize in cm
    if 'figsize' in kwargs:
        kwargs['figsize'] = (cm2inch(kwargs['figsize'][0])[0], cm2inch(kwargs['figsize'][1])[0])

    # initialize figure
    fig = plt.figure(*args, **kwargs)

    # event callbacks
    fig = _set_figure_methods(fig)
    return fig

def _set_figure_methods(fig):
    """adds onclick functionality to figures"""

    #  mouse events
    cid1 = fig.canvas.mpl_connect('button_press_event', _onclick)
    # cid2 = fig.canvas.mpl_connect('resize_event', _onmove)

    # grid function
    fig._grid = False
    fig.grid  = MethodType(_grid, fig)

    return fig

def set_onclick_save_defaults(format='svg', resolution=300, folder=None):
    """Set the default format for saving figures using Middle click press on a figure.

    Args:
        format (string, optional): 'svg' or 'png'.
        resolution (int, optional): resolution for png images in dpi.
        folder (string or pathlib.Path): folderpath to save figures.
    
    Returns:
        None
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

def _onclick(event):
    """This function is called every time a mouse key is pressed over a figure.

    Middle click:
        prints cursor position in terms of figure coordinates.
    Right click:
        x value is copied to the clipboard.
    Left click OR (y + Right click):
        y value is copied to the clipboard.

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
    # try:
    #     onclick_round_y
    # except NameError:
    #     onclick_round_y = 2

    # try:
    #     onclick_round_x
    # except NameError:
    #     onclick_round_x = 2

    ################
    # double click #
    ################
    if event.dblclick:
        # w, h = plt.gcf().get_size_inches()

        # x = event.x/(w*100)
        # y = event.y/(h*100)

        # _copy2clipboard([x, y])
        pass
    ###############
    # middle click #
    ###############
    elif event.button == 2:
        w, h = plt.gcf().get_size_inches()

        x = event.x/(w*100)
        y = event.y/(h*100)

        _copy2clipboard(str([x, y]))

        # # get variables        
        # global onclick_fig_format, onclick_resolution, onclick_folder
        # try:
        #     onclick_fig_format
        # except NameError:
        #     onclick_fig_format = 'png'
        # try:
        #     onclick_folder
        # except NameError:
        #     onclick_folder = Path.cwd()
        # try:
        #     onclick_resolution
        # except NameError:
        #     onclick_resolution = 300

        # with tempfile.NamedTemporaryFile("r+b", delete=True) as fd:
        #     if onclick_fig_format == 'svg':
        #         plt.savefig(fd)
        #         _svg2clipboard(fd)
        #     elif onclick_fig_format == 'png':
        #         plt.savefig(fd, dpi=onclick_resolution)
        #         _png2clipboard(fd)
        #     # print(f'figure copied to clipboard as {onclick_fig_format}')
    ###############
    # right click #
    ###############
    # elif event.key == 'y' or event.button == 3:
    elif event.button == 3:
        ax = event.inaxes
        if ax is not None: # check if mouse is over an axes
            try:
                lim   = ax.get_ylim()
                delta = lim[1] - lim[0]
                r21   = _round_to_1(delta)

                n = False
                if r21 < 1: 
                    n = _n_decimal_places(r21) + 2
                # if r21 < 1:
                #     i = 0
                #     while r21 < 1:
                #         r21 = r21*10
                #         i += 1
                #         if i == 11: break
                #     if i < 11:
                #         n = i + 2
                #     else:
                #         n = False
                elif r21 >= 1 and r21 < 2:
                    n = 4
                elif r21 >= 2 and r21 < 100:
                    n = 2
                elif r21 >= 100 and r21 < 800:
                    n = 1
                elif r21 >= 800:
                    n = 0
                
                final = round(event.ydata, n)
                if n == 0:
                    final = int(final)
                _copy2clipboard(str(final))
                # print('y coordinate copied to clipboard')
            except TypeError:
                pass
    ###############
    # left click #
    ###############
    elif event.button == 1:
        ax = event.inaxes
        if ax is not None: # check if mouse is over an axes
            try:
                lim = ax.get_xlim()
                delta = lim[1] - lim[0]
                r21 = _round_to_1(delta)

                if r21 < 1: 
                    n = _n_decimal_places(r21) + 2
                # if r21 < 1:
                #     i = 0
                #     while r21 < 1:
                #         r21 = r21*10
                #         i += 1
                #         if i == 11: break
                #     if i < 11:
                #         n = i + 2
                #     else:
                #         n = False
                elif r21 >= 1 and r21 < 2:
                    n = 4
                elif r21 >= 2 and r21 < 100:
                    n = 2
                elif r21 >= 100 and r21 < 800:
                    n = 1
                elif r21 >= 800:
                    n = 0

                final = round(event.xdata, n)
                if n == 0:
                    final = int(final)
                _copy2clipboard(str(final))
                # print('x coordinate copied to clipboard')
            except TypeError:
                pass
   
def _grid(self, visible=None):
    """show figure grid (not be confused with axes grid)"""
    if visible is None:
        if self._grid == False:
            visible = True
        else:
            visible = False
    elif visible == True:
        if self._grid != False:
            return
    elif visible == False:
        if self._grid == False:
            return

    if visible:
        self._grid = self.add_axes([0, 0, 1, 1], alpha=0.5)
        self._grid.patch.set_alpha(0.5)

        set_xticks(ax=self._grid, start=0, stop=1, step=0.1, pad=(0, 0), n_minor_ticks=1)
        set_yticks(ax=self._grid, start=0, stop=1, step=0.1, pad=(0, 0), n_minor_ticks=1)
        
        self._grid.grid(which='major', color='red', lw=0.5)
        self._grid.grid(which='minor', color='black', ls=':', lw=0.5)
    
        for v in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
            note(ax=self._grid, s=v, x=v, loc='lower', horizontalalignment='center')
            note(ax=self._grid, s=v, y=v, loc='left')
    else:
        self._grid.set_visible(False)
        del self._grid
        self._grid = False
    return

# %% ============================== subplots ============================== %% #
class Axes(list):
    def __init__(self, axes, nrows, ncols):
        """This class facilitates selecting axes generated by br.subplots().

        br.subplots() arrange axes in a grid with nrows and ncols.

        `axes` is expected to to be a list o axes. Indexes should go from left 
        to right, top to bottom, e. g., if nrows = 3 
        and ncols = 2, then, the top left axes will be axes[0], then the top right 
        one will be axes[1]. Second row starts from axes[2] and so on.

        The class adds extra methods to help selecting axes: 
        
        >>> row_i = axes.rows[i] to get a list of all axes in row i 
        >>> col_j = axes.cols[j] to get a list of all axes in column j
        
        >>> first_row = axes.first_row to get a list of all axes in the first row 
        >>> last_row  = axes.last_row to get a list of all axes in the last row 
        >>> first_col = axes.first_row to get a list of all axes in the first column
        >>> last_col  = axes.first_row to get a list of all axes in the last column 

        Args:
            axes (list): list of axes. Axes should 
            nrows, ncols (int): number of cols and rows

        """
        super().__init__(axes)
        self.nrows = nrows
        self.ncols = ncols

    @property
    def cols(self):
        final = []
        for i in range(self.ncols):
            final.append([self[_] for _ in np.arange(0, self.ncols*self.nrows, self.ncols) + i])
        return final
    @cols.setter
    def cols(self, value):
        raise AttributeError('Cannot set object.')
    @cols.deleter
    def cols(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def rows(self):
        final = []
        for i in range(self.nrows):
            final.append([self[_] for _ in np.arange(0, self.ncols) + i*self.ncols])
        return final
    @rows.setter
    def rows(self, value):
        raise AttributeError('Cannot set object.')
    @rows.deleter
    def rows(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def first_row(self):
        return [self[_] for _ in np.arange(0, self.ncols)]
    @first_row.setter
    def first_row(self, value):
        raise AttributeError('Cannot set object.')
    @first_row.deleter
    def first_row(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def last_row(self):
        return [self[_] for _ in np.arange(0, self.ncols) + (self.nrows-1)*self.ncols]
    @last_row.setter
    def last_row(self, value):
        raise AttributeError('Cannot set object.')
    @last_row.deleter
    def last_row(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def first_col(self):
        return [self[_] for _ in np.arange(0, self.ncols*self.nrows, self.ncols)]
    @first_col.setter
    def first_col(self, value):
        raise AttributeError('Cannot set object.')
    @first_col.deleter
    def first_col(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def last_col(self):
        return [self[_] for _ in np.arange(0, self.ncols*self.nrows, self.ncols) + self.ncols-1]
    @last_col.setter
    def last_col(self, value):
        raise AttributeError('Cannot set object.')
    @last_col.deleter
    def last_col(self):
        raise AttributeError('Cannot delete object.')

def subplots(nrows, ncols, sharex=False, sharey=False, hspace=0.3, wspace=0.3, width_ratios=None, height_ratios=None, layout=None, **fig_kw):
    """Create a figure and a set of subplots in a grid. Wrapper for `plt.subplots()`_.

    The difference between this function and plt.subplots is that this function 
    returns the axes in a list all the time, while plt.subplots returns axes as list
    of lists if nrows >= 2 and ncols >= 2).

    The axes indexes go from left to right, top to bottom, e. g., if nrows = 3 
    and ncols = 2, then, the top left axes will be axes[0], then the top right 
    one will be axes[1]. Second row starts from axes[2] and so on.

    The returned `axes` is a list with extra methods to help selecting axes: 
    
        >>> row_i = axes.rows[i] to get a list of all axes in row i 
        >>> col_j = axes.cols[j] to get a list of all axes in column j
        
        >>> first_row = axes.first_row to get a list of all axes in the first row 
        >>> last_row  = axes.last_row to get a list of all axes in the last row 
        >>> first_col = axes.first_row to get a list of all axes in the first column
        >>> last_col  = axes.first_row to get a list of all axes in the last column 

    Args:
        nrows, ncols (int): Number of rows/columns of the subplot grid.
        sharex, sharey (bool or str, optional): Share the x or y axis between 
            axes. Options are: True, False, 'all', 'row', 'col'. Default is False
        hspace, wspace (float, optional): The amount of height/width reserved for space 
            between subplots, expressed as a fraction of the average axis 
            height/width. Default is 0.3.
        width_ratios, height_ratios (list, optional): Defines the relative heights/widths of the 
            ows/columns. Each row/column gets a relative height/width of 
            ratios[i] / sum(ratios). If not given, all rows/columns will have 
            the same width. 
        layout (str or None, optional): use 'constrained' to automatically adjusts 
            subplots so that decorations like tick labels, legends, and 
            colorbars do not overlap, while still preserving the logical layout.
            Default is None. Layout='constrained' overwrites hspace, wspace, 
            width_ratios, height_ratios
        **fig_kw: All additional keyword arguments are passed to the plt.figure call
    
    Returns:
        fig, axes
    
    .. _plt.subplots(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html
    """
    #########
    # check #
    #########
    assert nrows > 0, 'nrows must be > 0'
    assert ncols > 0, 'ncols must be > 0'

    ###############
    # gridspec_kw #
    ###############
    # ratios
    if width_ratios is None:
        width_ratios=[1]*ncols
    if height_ratios is None:
        height_ratios=[1]*nrows

    # final
    if layout is None:
        gridspec_kw = {}
        gridspec_kw['wspace']        = wspace 
        gridspec_kw['hspace']        = hspace 
        gridspec_kw['width_ratios']  = width_ratios 
        gridspec_kw['height_ratios'] = height_ratios 
    else:
        gridspec_kw = None

    ##########
    # fig_kw #
    ##########
    if 'figsize' in fig_kw:
        fig_kw['figsize'] = cm2inch(fig_kw['figsize'][0], fig_kw['figsize'][1])

    ############
    # subplots #
    ############
    fig, _axes = plt.subplots(nrows, ncols, 
                              sharex=sharex,
                              sharey=sharey,
                              layout=layout,
                              gridspec_kw=gridspec_kw, **fig_kw)
    fig = _set_figure_methods(fig)

    #################
    # fatten output #
    #################
    axes = [0]*(ncols*nrows)
    if ncols == 1 and nrows == 1:
        _axes.row = 0
        _axes.col = 0
        axes[0] = _axes
    elif nrows == 1:
        for i, ax in enumerate(_axes):
            ax.row = 0
            ax.col = i
            axes[i] = ax
    elif ncols == 1:
        for i, ax in enumerate(_axes):
            ax.row = i
            ax.col = 0
            axes[i] = ax
    else:
        for i, row in enumerate(_axes):
            for j, ax in enumerate(row):
                ax.row = i
                ax.col = j
                axes[i*ncols+j] = ax
    # try:
    #     for i, ax in enumerate(flatten(_axes)):
    #         axes[i] = ax
    # except TypeError:
    #     axes[0] = _axes

    ####################
    # new axes methods #
    ####################
    for ax in axes:
        ax.remove_xlabel = MethodType(remove_xlabel, ax)
        ax.remove_ylabel = MethodType(remove_ylabel, ax)

        ax.remove_xticklabels = MethodType(remove_xticklabels, ax)
        ax.remove_yticklabels = MethodType(remove_yticklabels, ax)

    return fig, Axes(axes, ncols=ncols, nrows=nrows)
# %%

# %% ==================================== font ============================ %% #
def publication_font(size=9):
    """Set font options for publication quality font (cmr10).

    Change this function accordingly to your preferences.

    Returns:
        None
    """
    # set default font
    matplotlib.rcParams['font.family']          = 'cmr10'
    matplotlib.rcParams['font.size']            = size
    plt.rcParams["axes.formatter.use_mathtext"] = True

    # math text
    matplotlib.rcParams['mathtext.fontset'] = 'cm'  # Change mathtext to CMR ("latex")
    matplotlib.rcParams['svg.fonttype']     = 'none'  # Change text to text, and not paths.

def default_font(size=10):
    """Set font options to Matplotlib's default.

    Returns:
        None
    """
    # set default font
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.size']   = size
    plt.rcParams["axes.formatter.use_mathtext"] = False

    # math text
    matplotlib.rcParams['mathtext.fontset'] = 'dejavusans'
    matplotlib.rcParams['svg.fonttype'] = 'path'
# %%
    
# %% ================================== text ============================== %% #
def rtext(x, s, yoffset=0, xoffset=0, ax=None, copy_color=False, **kwargs):
    """Create text at x coordinate on the right side of the last curve plotted

    Args:
        x (number): x coordinate for the text
        s (str): text string
        xoffset, yoffset (number, optional): offset between the curve and the text.
        ax (axes, optional): axes to put the text on. If None, the current axes 
        will be used. Default is None
        copy_color (bool, optional): if True, the color of the text will match
        the color of the curve. Default is False
        **kwargs: kwargs are passed to `plt.text()`_

    Returns:
        matplotlib.text.Text
    
    .. _plt.text(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    """
    if ax is None:
        ax = plt.gca()

    # get lines
    temp = ax.get_lines()[-1]
    y = temp.get_ydata()[_index(temp.get_xdata(), x)] + yoffset

    # get color
    if 'color' not in kwargs and copy_color:
        kwargs['color'] = temp.get_color()

    # fix align
    if 'ha' or 'horizontalalignment' not in kwargs:
        kwargs['ha'] = 'left'
    if 'va' or 'verticalalignment' not in kwargs:
        kwargs['va'] = 'center'

    # print(f'x={x + xoffset}, y={y}')
    return ax.text(x + xoffset, y, s, kwargs)

def ltext(x, s, yoffset=0, xoffset=0, ax=None, copy_color=False, **kwargs):
    """Create text at x coordinate on the left side of the last curve plotted

    Args:
        x (number): x coordinate for the text
        s (str): text string
        xoffset, yoffset (number, optional): offset between the curve and the text.
        ax (axes, optional): axes to put the text on. If None, the current axes 
        will be used. Default is None
        copy_color (bool, optional): if True, the color of the text will match
        the color of the curve. Default is False
        **kwargs: kwargs are passed to `plt.text()`_

    Returns:
        matplotlib.text.Text
    
    .. _plt.text(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    """
    if ax is None:
        ax = plt.gca()

    # get lines
    temp = ax.get_lines()[-1]
    y = temp.get_ydata()[_index(temp.get_xdata(), x)] + yoffset

    # get color
    if 'color' not in kwargs and copy_color:
        kwargs['color'] = temp.get_color()

    # fix align
    if 'ha' or 'horizontalalignment' not in kwargs:
        kwargs['ha'] = 'right'
    if 'va' or 'verticalalignment' not in kwargs:
        kwargs['va'] = 'center'
    
    # write
    # print(f'x={x + xoffset}, y={y}')
    return ax.text(x + xoffset, y, s, **kwargs)

def note(s, loc='upper left', x=None, y=None, ax=None, coord='axes', **kwargs):
    """Write text to a pre-defined locations (auto-position text). Wrapper for `plt.text()`_.

    Text is written in axes coord (i.e. text does not move with the data). The 
    coordinate system of the Axes; (0, 0) is bottom left of the axes, and (1, 1)
    is top right of the axes.

    Args:
        s (str): text string
        loc (str, optional): location of the text. Supported locations are the
        same ones for plt.legend(): 'upper right', 'upper left', 'lower left', 
        'lower right', 'right', 'center left', 'center right', 'lower center', 
        'upper center', 'center'. Note that 'best' is not a supported value.
        default is 'upper left'
        x, y (number): x and y in axes coordinates (0, 1). Overwrites loc.
        ax (axes, optional): axes to put the text on. If None, the current axes 
        will be used. Default is None
        **kwargs: kwargs are passed to `plt.text()`_

    Returns:
        matplotlib.text.Text
    
    .. _plt.text(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    """
    if ax is None:
        ax = plt.gca()

    if x is not None and y is not None:
        # return ax.text(x, y, s, transform=fig.transAxes, **kwargs)
        return ax.text(x, y, s, transform=ax.transAxes, **kwargs)

    # get figure size
    # w, h = ax.get_figure().get_size_inches()

    # copy legend if any
    leg = ax.get_legend()

    # create ghost legend
    # ax.autoscale(False) # do we really need this?
    # line, = ax.plot([])
    temp  = ax.legend(labels=[''], loc='upper left')
    
    # _ = plt.gca().add_artist(leg)
    # leg  = ax.legend(line, [''], loc=loc)

    # get legend coordinates in display coord (I think)
    # ax = axes[0]
    # leg = axes[0].legend()
    p = temp.get_window_extent(ax.get_figure().canvas.get_renderer())
    _x, _y = (p.p0[0], p.p1[1])

    # remove ghost legend and add old legend back
    temp.remove()
    if leg is not None:
        ax.add_artist(leg)

    # if x, y are defined, transform it from data coord to display coord
    if x is not None:
        # _x, _ = ax.transData.transform((x, 0))
        # _x, _ = ax.transLimits.transform((x, 0))
        if 'ha' not in kwargs and 'horizontalalignment' not in kwargs:
            kwargs['ha'] = 'left'
    if y is not None:
        # _, _y = ax.transData.transform((0, y))
        # _, y = ax.transLimits.transform((0, y))
        if 'va' not in kwargs and 'verticalalignment' not in kwargs:
            kwargs['va'] = 'center'
    # print('_x, + _y')
    # print((_x, _y))
    # print('end')

    # transform from display coord to axes coord
    inv  = ax.transAxes.inverted()
    if y is not None:
        _x, _ = inv.transform((_x, _y))
        _y    = y
    elif x is not None:
        _, _y = inv.transform((_x, _y))
        _x    = x
    else:
        _x, _y = inv.transform((_x, _y))
    _x = round(_x, 4)
    _y = round(_y, 4)
    # (0.021576382861200977, 0.9672013006569351)

    # convert from upper right to loc
    if y is None:
        if 'upper' in loc:
            pass
        elif 'lower' in loc:
            _y = 1-_y
        elif loc.startswith('center'):
            _y = 0.5

    if x is None:
        if 'left' in loc:
            pass
        elif 'right' in loc:
            _x = 1 - _x
        elif loc.endswith('center'):
            _x = 0.5

    # adjust alignment
    if loc == 'best':
        raise ValueError('`best` is not a supported value')
    
    if 'ha' not in kwargs and 'horizontalalignment' not in kwargs:
        if 'right' in loc:
            kwargs['ha'] = 'right'
        elif 'left' in loc:
            kwargs['ha'] = 'left'
        elif loc.endswith('center'):
            kwargs['ha'] = 'center'
        else:
            kwargs['ha'] = 'center'

    if 'va' not in kwargs and 'verticalalignment' not in kwargs:
        if 'upper' in loc:
            kwargs['va'] = 'top'
        elif 'lower' in loc:
            kwargs['va'] = 'bottom'
        elif loc.startswith('center'):
            kwargs['va'] = 'center'
        else:
            kwargs['va'] = 'center'

    # # transform from axes coord to display coord
    # _x, _y = ax.transAxes.transform((x, y))
    # # transform from display coord to data coord (just for printing)
    # inv    = ax.transData.inverted()
    # _x, _y = inv.transform((_x, _y))
    # # print(_x)
    # print(f'x={_x}, y={_y}')

    # ax.autoscale(True)  # do we really need this

    # final
    # print(kwargs)
    return ax.text(_x, _y, s, transform=ax.transAxes, **kwargs)
    # return ax.text(x, y, s, transform=fig.transFigure, **kwargs)
    # return ax.text(x, y, s, **kwargs)

def label_axes(axes, loc='upper left', color=None, start_letter='a', start_n=0):
    """put letters from `(a)` to `(z)` in axes.

    Args:
        axes (axes object or list): axes list
        locs (str or list, optional): The location of the text. If list, list
            must be same length as number of axes. Options are
            'center',
            'upper left', 'upper right', 'lower left', 'lower right', 
            'upper center', 'lower center', 'center left', 'center right'. 
            Default is 'best'.
        color (str, list, or None, optional): color to be used. If None, the
            default color is used. If list, list must be the same lenght as 
            number of axes. Default is None.
        start_letter (string, optional): letter of the first panel. Default is `a`
        start_n (int, optional): number of the first panel. This number is used 
            only if n>0. Default is 0. This is used for figures where the number
            number of axes exceed the number of letters in the alphabet. 

    Returns:
        list with text object
    """
    # fix axes
    if isinstance(axes, Iterable) == False:
        axes = [axes, ]

    # fix locs
    if isinstance(loc, str):
        loc = [loc]*len(axes)

    # fix colors
    if isinstance(color, str) or color is None:
        color = [color]*len(axes)
    
    assert len(axes) == len(loc),   'number of axes must be the same as number of locs'
    assert len(axes) == len(color), 'number of axes must be the same as number of colors'

    n = start_n
    j = ascii_lowercase.index(start_letter)
    jstart = j
    final = []
    for i, ax in enumerate(axes):
        if (i+jstart) % 26 == 0 and i != 0:
            j = 0
            jstart = 0
            n += 1

        if n > 0:
            letter = ascii_lowercase[j] + str(n)
        else:
            letter = ascii_lowercase[j]
        final.append(note(s='(' + letter + ')', ax=ax, loc=loc[i], color=color[i]))
        j += 1
    return final

# %% add to axes
def _rtext(self, x, s, yoffset=0, xoffset=0, copy_color=False, **kwargs):
    """Create text at x coordinate on the right side of the last curve plotted

    Args:
        x (number): x coordinate for the text
        s (str): text string
        xoffset, yoffset (number, optional): offset between the curve and the text.
        copy_color (bool, optional): if True, the color of the text will match
        the color of the curve. Default is False
        **kwargs: kwargs are passed to `plt.text()`_

    Returns:
        matplotlib.text.Text
    
    .. _plt.text(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    """
    return rtext(x=x, s=s, yoffset=yoffset, xoffset=xoffset, ax=self, copy_color=copy_color, **kwargs)
mpl.axes.Axes.rtext = _rtext

def _ltext(self, x, s, yoffset=0, xoffset=0, copy_color=False, **kwargs):
    """Create text at x coordinate on the left side of the last curve plotted

    Args:
        x (number): x coordinate for the text
        s (str): text string
        xoffset, yoffset (number, optional): offset between the curve and the text.
        copy_color (bool, optional): if True, the color of the text will match
        the color of the curve. Default is False
        **kwargs: kwargs are passed to `plt.text()`_

    Returns:
        matplotlib.text.Text
    
    .. _plt.text(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    """
    return ltext(x=x, s=s, yoffset=yoffset, xoffset=xoffset, ax=self, copy_color=copy_color, **kwargs)
mpl.axes.Axes.ltext = _ltext

def _note(self, s, loc='upper left', x=None, y=None, **kwargs):
    """Write text to a pre-defined locations. Wrapper for `plt.text()`_.

    Text is written in axes coord (i.e. text does not move with the data). The 
    coordinate system of the Axes; (0, 0) is bottom left of the axes, and (1, 1)
    is top right of the axes.

    Args:
        s (str): text string
        loc (str, optional): location of the text. Supported locations are the
        same ones for plt.legend(): 'upper right', 'upper left', 'lower left', 
        'lower right', 'right', 'center left', 'center right', 'lower center', 
        'upper center', 'center'. Note that 'best' is not a supported value.
        default is 'upper left'
        x, y (number): x and y in axes coordinates (0, 1). Overwrites loc.
        **kwargs: kwargs are passed to `plt.text()`_

    Returns:
        matplotlib.text.Text
    
    .. _plt.text(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
    """
    return note(s=s, loc=loc, x=x, y=y, ax=self, **kwargs)
mpl.axes.Axes.note = _note
# %%

# %% ================================= label ============================== %% #
def remove_xlabel(ax):
    """remove xlabel from axes"""
    ax.set_xlabel('')
    return
mpl.axes.Axes.remove_xlabel = remove_xlabel

def remove_ylabel(ax):
    """remove ylabel from axes"""
    ax.set_ylabel('')
    return
mpl.axes.Axes.remove_ylabel = remove_ylabel
# %%

# %% ================================= axes =============================== %% #
def merge_axes(ax1, ax2):
    """merge axes created via plt.subplots

    Args:
        ax1, ax2 (axes object): axes to be merged.

    Returns
        merged axes object
    """
    # check if axes came from the same figure
    assert ax1.figure == ax2.figure, 'axes must come from the same figure'

    # check if axes are neighbors and the same size
    xmin1 = ax1._subplotspec.rowspan.start
    xmax1 = ax1._subplotspec.rowspan.stop
    ymin1 = ax1._subplotspec.colspan.start
    ymax1 = ax1._subplotspec.colspan.stop

    xmin2 = ax2._subplotspec.rowspan.start
    xmax2 = ax2._subplotspec.rowspan.stop
    ymin2 = ax2._subplotspec.colspan.start
    ymax2 = ax2._subplotspec.colspan.stop

    assert xmin1 == xmin2 or ymin1 == ymin2, 'axes must be neighbors'
    if xmin1 == xmin2:
        assert (xmax1 - xmin1) == (xmax2 - xmin2), 'axes must have the same height'
    elif ymin1 == ymin2:
        assert (ymax1 - ymin1) == (ymax2 - ymin2), 'axes must have the same width'

    # get new grid coordinates
    pos = ax1.get_position().get_points()
    xmin = min(xmin1, xmin2)
    xmax = max(xmax1, xmax2)
    ymin = min(ymin1, ymin2)
    ymax = max(ymax1, ymax2)

    # add subplot
    gs  = ax1.get_gridspec()
    ax3 = ax1.figure.add_subplot(gs[xmin:xmax, ymin:ymax])

    # remove the underlying axes
    ax1.remove()
    ax2.remove()

    return ax3

def _mergewith(self, ax):
    """merge axes created via plt.subplots

    Args:
        ax (axes object): axes to be merged with.

    Returns
        merged axes object
    """
    return merge_axes(self, ax)
mpl.axes.Axes.mergewith = _mergewith

def add_lspace(value, ax=None):
    """add space to the left of an axes

    Args:
        value (float): value in figure units (0 - 1).
        ax (axes): axes to apply change. If None, current axes is selected.
            Default is None.

    Returns:
        None
    """
    if ax is None:
        ax = plt.gca()

    pos = ax.get_position().get_points()

    final = [pos[0][0] + value, 
             pos[0][1], 
             0, 
             0]
    
    final[2] = pos[1][0] - final[0]
    final[3] = pos[1][1] - final[1]

    ax.set_position(final)
    return

def add_rspace(value, ax=None):
    """add space to the right of an axes

    Args:
        value (float): value in figure units (0 - 1).
        ax (axes): axes to apply change. If None, current axes is selected.
            Default is None.

    Returns:
        None
    """
    if ax is None:
        ax = plt.gca()
    pos = ax.get_position().get_points()

    final = [pos[0][0], 
             pos[0][1], 
             0, 
             0]
    
    final[2] = pos[1][0] - final[0]  + value
    final[3] = pos[1][1] - final[1]

    ax.set_position(final)
    return

def add_tspace(value, ax=None):
    """add space to the top of an axes

    Args:
        value (float): value in figure units (0 - 1).
        ax (axes): axes to apply change. If None, current axes is selected.
            Default is None.

    Returns:
        None
    """
    if ax is None:
        ax = plt.gca()
    pos = ax.get_position().get_points()

    final = [pos[0][0], 
             pos[0][1], 
             0, 
             0]
    
    final[2] = pos[1][0] - final[0]
    final[3] = pos[1][1] - final[1] + value

    ax.set_position(final)
    return

def add_bspace(value, ax=None):
    """add space to the left of an axes

    Args:
        value (float): value in figure units (0 - 1).
        ax (axes): axes to apply change. If None, current axes is selected.
            Default is None.

    Returns:
        None
    """
    if ax is None:
        ax = plt.gca()
    pos = ax.get_position().get_points()

    final = [pos[0][0], 
             pos[0][1] + value, 
             0, 
             0]
    
    final[2] = pos[1][0] - final[0]
    final[3] = pos[1][1] - final[1]

    ax.set_position(final)
    return

def xmove(value, ax=None):
    """move axes in the x direction

    Args:
        value (float): value in figure units (0 - 1).
        ax (axes): axes to apply change. If None, current axes is selected.
            Default is None.

    Returns:
        None
    """
    if ax is None:
        ax = plt.gca()
    add_lspace(ax=ax, value=value)
    add_rspace(ax=ax, value=value)
    return

def ymove(value, ax=None):
    """move axes in the y direction

    Args:
        value (float): value in figure units (0 - 1).
        ax (axes): axes to apply change. If None, current axes is selected.
            Default is None.

    Returns:
        None
    """
    if ax is None:
        ax = plt.gca()
    add_tspace(ax=ax, value=value)
    add_bspace(ax=ax, value=value)
    return

def sharex(axes):
    """Connect x axis. Cannot be used if the x-axis is already being shared with another Axes"""
    assert isinstance(axes, Iterable), 'input must be a list of axes'
    for i, ax in enumerate(axes):
        ax.sharex(axes[0])
    return

def sharey(axes):
    """Connect y axis. Cannot be used if the y-axis is already being shared with another Axes"""
    assert isinstance(axes, Iterable), 'input must be a list of axes'
    for i, ax in enumerate(axes):
        ax.sharey(axes[0])
    return

# %% add to axes
def _add_lspace(self, value):
    """add space to the left of an axes

    Args:
        value (float): value in figure units (0 - 1).

    Returns:
        None
    """
    return add_lspace(value=value, ax=self)
mpl.axes.Axes.add_lspace = _add_lspace

def _add_rspace(self, value):
    """add space to the right of an axes

    Args:
        value (float): value in figure units (0 - 1).

    Returns:
        None
    """
    return add_rspace(value=value, ax=self)
mpl.axes.Axes.add_rspace = _add_rspace

def _add_tspace(self, value):
    """add space to the top of an axes

    Args:
        value (float): value in figure units (0 - 1).

    Returns:
        None
    """
    return add_tspace(value=value, ax=self)
mpl.axes.Axes.add_tspace = _add_tspace

def _add_bspace(self, value):
    """add space to the bottom of an axes

    Args:
        value (float): value in figure units (0 - 1).

    Returns:
        None
    """
    return add_bspace(value=value, ax=self)
mpl.axes.Axes.add_bspace = _add_bspace

def _xmove(self, value):
    """move axes in the x direction

    Args:
        value (float): value in figure units (0 - 1).
        ax (axes): axes to apply change. If None, current axes is selected.
            Default is None.

    Returns:
        None
    """
    return xmove(value=value, ax=self)
mpl.axes.Axes.xmove = _xmove

def _ymove(self, value):
    """move axes in the y direction

    Args:
        value (float): value in figure units (0 - 1).
        ax (axes): axes to apply change. If None, current axes is selected.
            Default is None.

    Returns:
        None
    """
    return ymove(value=value, ax=self)
mpl.axes.Axes.ymove = _ymove
# %%

# %% ================================= ticks ============================== %% #
# get ticks showing
def get_xticks_showing(ax):
    """get x ticks that are sowing in a plot
    
    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.

    Returns:
        list
    """
    return [x for x in ax.get_xticks() if x >= ax.get_xlim()[0] and x <= ax.get_xlim()[1]]

def get_yticks_showing(ax):
    """get y ticks that are sowing in a plot
    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.

    Returns:
        list
    """
    return [y for y in ax.get_yticks() if y >= ax.get_ylim()[0] and y <= ax.get_ylim()[1]]

# set ticks
def set_xticks(ax=None, start=None, stop=None, nticks=None, step=None, ticks=None, labels=None, pad=None, n_minor_ticks=None, minor_step=None, minor_ticks=None, fontproperties=None, **kwargs):
    """Set x ticks of a plot.

    Major and minor ticks position, location, and direction changes are applied to all
    shared axis. However, ticklabel changes are only applyied to the selected
    axes. Ticklabel parameters: `labelsize`, `labelcolor`, `labelfontfamily`, 
    `labelbottom`, `labeltop`, `labelleft`, `labelright`, `labelrotation`.

    Usage:
        For setting up major ticks:

            set_ticks(start, stop, nticks)
            set_ticks(start, stop, step)
            set_ticks(start, nticks, step)  # `stop` is calculated
            set_ticks(nticks)                    # `start` and `stop` are kept from the plot
            set_ticks(ticks)                      # manually define each tick position
            set_ticks(ticks, labels)              # manually define each tick position and label

        For setting up minor ticks:

            set_ticks(n_minor_ticks)             # equally spaced minor ticks
            set_ticks(minor_step)           # equally spaced minor ticks
            set_ticks(minor_ticks)               # manually define each minor tick position

    Args:
        ax (matplotlib.axes, optional): axes object to set ticks.

        *Major ticks*
        start, stop (number, optional): major ticks start/stop value. Use None
            to not keep current start/stop value. Default is None.
        nticks, step (int/number, optional): Number of major ticks or separation between ticks. 
            nticks overwrites step. Default is None.
        ticks, labels (Iterable, optional): Iterable with tick positions and labels.
        Overwrites start, stop, nticks, step.
        
        *Padding (plotting limits)*
        pad (number or tuple, optional): padding between plot edge and 
            closest tick in terms of ticks separation. If tuple, first value 
            applies to the left edge and second value to the right edge. 
            Typically, padding should be between a value between 0 and 1. If pad
            is None, plot limits remain unchanged. However, if pad is None and 
            start/stop are defined, pad is set to (0, 0). Default is None.

        *Minor ticks*
        n_minor_ticks, minor_step (int/number, optional): Number of minor 
        ticks between two major ticks or minor tick separation. n_minor_ticks
        overwrites minor_step. Default is None.
        minor_ticks (Iterable, optional): Iterable with minor tick positions.
        Overwrites n_minor_ticks and minor_step.

        *font*
        fontproperties: Label ticks font. Use ``matplotlib.font_manager.FontProperties``.

        *ticks parameters*
        **kwargs: kwargs are passed to `ax.tick_params()`_ that sets tick parameters.

        If not specified, the following parameters are passed to `ax.tick_params()`_:

        Args:
            which ('major', 'minor', 'both', optional): The group of ticks to which the parameters 
                are applied. Default is 'both'.
            direction ('in', 'out', 'inout', optional): Puts ticks inside the 
                Axes, outside the Axes, or both. Default is 'in'
            bottom, top (bool, optional): Whether to draw the bottom/top ticks.
                Default is True.

    Returns:
        None
    
    .. _ax.tick_params(): https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.tick_params.html
    """
    ############
    # get axis #
    ############
    if ax is None:
        ax = plt.gca()

    ######################
    # default parameters #
    ######################
    if 'axis' not in kwargs:
        kwargs['axis'] = 'x'
    if 'which' not in kwargs:
        kwargs['which'] = 'both'
    # if 'direction' not in kwargs:
    #     kwargs['direction'] = 'in'
    # if 'bottom' not in kwargs:
    #     kwargs['bottom'] = True
    # if 'top' not in kwargs:
    #     kwargs['top'] = True

    #####################
    # major ticks setup #
    #####################
    if ticks is None:
        if nticks is not None:
            assert nticks > 1, 'nticks needs to be > 1.'
        if start is not None or stop is not None or nticks is not None or step is not None:
            if start is not None and stop is not None and nticks is not None:
                ticks = np.linspace(start, stop, nticks)
            elif start is not None and stop is not None and step is not None:
                ticks   = np.arange(start, stop + step*0.1, step) 
            elif start is not None and nticks is not None and step is not None:
                stop = start + nticks*step
                ticks = np.linspace(start, stop, nticks)
            elif nticks is not None:
                ticks_showing = get_xticks_showing(ax=ax)
                start = ticks_showing[0]
                stop  = ticks_showing[-1]
                ticks = np.linspace(start, stop, nticks)
            else:
                raise ValueError('Invalid input. See help(set_xticks) for help.')
            
            # set ticks
            _ = ax.xaxis.set_ticks(ticks)

            # set label font
            if fontproperties is not None:
                _ = ax.set_xticklabels([str(i) for i in ticks], fontproperties=fontproperties, visible=True)
            
            # if start/stop is defined, select pad to (0, 0)
            # if pad is None:
            #     pad = (0, 0)
    else:
        assert isinstance(ticks, Iterable), 'ticks must be iterable'
        assert len(ticks) > 1, 'ticks must have length > 1'
        _ = ax.xaxis.set_ticks(ticks)

        # set labels
        if labels is not None:
            assert isinstance(labels, Iterable), 'ticks must be iterable'
            assert len(ticks) == len(labels), 'len(ticks) must be the same as len(labels)'
            if fontproperties is not None:
                _ = ax.set_xticklabels(labels, fontproperties=fontproperties, visible=True)
            else:
                _ = ax.set_xticklabels(labels)
        else:
            if fontproperties is not None:
                _ = ax.set_xticklabels([str(i) for i in ticks], fontproperties=fontproperties, visible=True)

    ##########
    # limits #
    ##########
    if pad is not None:
        if isinstance(pad, Iterable) == False:
            assert is_number(pad), 'pad must be a number or a tuple/list with length 2'
            pad = (pad, pad)

        # check
        assert len(pad) == 2, 'pad must be a number or a tuple/list with length 2'
        assert is_number(pad[0]), 'pad must be a number or a tuple/list with length 2'
        assert is_number(pad[1]), 'pad must be a number or a tuple/list with length 2'

        # set limits
        step = ticks[1]  - ticks[0]
        min_lim   = ticks[0]  - step*pad[0]
        max_lim   = ticks[-1] + step*pad[1]
        ax.set_xlim((min_lim, max_lim))#, auto=False)

    ###############
    # minor ticks #
    ###############
    if minor_ticks is None:
        if n_minor_ticks is not None:
            ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
        elif minor_step is not None:
            step     = ticks[1] - ticks[0]
            n_minor_ticks = int(round(step/minor_step))
            ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
    else:
        assert isinstance(ticks, Iterable), 'minor_ticks must be iterable'
        assert len(ticks) > 0, 'minor_ticks must have length > 0'
        _ = ax.xaxis.set_ticks(minor_ticks, minor=True)

    ###################
    # tick parameters #
    ###################
    labelparams = ['labelsize', 'labelcolor', 'labelfontfamily', 'labelbottom', 'labeltop', 'labelleft', 'labelright', 'labelrotation']
    labelparams = {param:kwargs[param] for param in labelparams if param in kwargs}
    _ = [kwargs.pop(param) for param in labelparams]
    ax.tick_params(**labelparams)

    for _ax in ax.get_shared_x_axes().get_siblings(ax):
        _ax.tick_params(**kwargs)

    return

def set_yticks(ax=None, start=None, stop=None, nticks=None, step=None, ticks=None, labels=None, pad=None, n_minor_ticks=None, minor_step=None, minor_ticks=None, fontproperties=None, **kwargs):
    """Set y ticks of a plot.

    Major and minor ticks position, location, and direction changes are applied to all
    shared axis. However, ticklabel changes are only applyied to the selected
    axes. Ticklabel parameters: `labelsize`, `labelcolor`, `labelfontfamily`, 
    `labelbottom`, `labeltop`, `labelleft`, `labelright`, `labelrotation`.


    Usage:
        For setting up major ticks:

            set_ticks(start, stop, nticks)
            set_ticks(start, stop, step)
            set_ticks(start, nticks, step)  # `stop` is calculated
            set_ticks(nticks)                    # `start` and `stop` are kept from the plot
            set_ticks(ticks)                      # manually define each tick position
            set_ticks(ticks, labels)              # manually define each tick position and label

        For setting up minor ticks:

            set_ticks(n_minor_ticks)             # equally spaced minor ticks
            set_ticks(minor_step)           # equally spaced minor ticks
            set_ticks(minor_ticks)               # manually define each minor tick position

    Args:
        ax (matplotlib.axes, optional): axes object to set ticks.

        *Major ticks*
        start, stop (number, optional): major ticks start/stop value. Use None
            to not keep current start/stop value. Default is None.
        nticks, step (int/number, optional): Number of major ticks or separation between ticks. 
            nticks overwrites step. Default is None.
        ticks, labels (Iterable, optional): Iterable with tick positions and labels.
        Overwrites start, stop, nticks, step.
        
        *Padding (plotting limits)*
        pad (number or tuple, optional): padding between plot edge and 
            closest tick in terms of ticks separation. If tuple, first value 
            applies to the left edge and second value to the right edge. 
            Typically, padding should be between a value between 0 and 1. If pad
            is None, plot limits remain unchanged. However, if pad is None and 
            start/stop are defined, pad is set to (0, 0). Default is None.

        *Minor ticks*
        n_minor_ticks, minor_step (int/number, optional): Number of minor 
        ticks between two major ticks or minor tick separation. n_minor_ticks
        overwrites minor_step. Default is None.
        minor_ticks (Iterable, optional): Iterable with minor tick positions.
        Overwrites n_minor_ticks and minor_step.

        *font*
        fontproperties: Label ticks font. Use ``matplotlib.font_manager.FontProperties``.

        *ticks parameters*
        **kwargs: kwargs are passed to `ax.tick_params()`_ that sets tick parameters.

        If not specified, the following parameters are passed to `ax.tick_params()`_:

        Args:
            which ('major', 'minor', 'both', optional): The group of ticks to which the parameters 
                are applied. Default is 'both'.
            direction ('in', 'out', 'inout', optional): Puts ticks inside the 
                Axes, outside the Axes, or both. Default is 'in'
            left, right (bool, optional): Whether to draw the left/right ticks.
                Default is True.

    Returns:
        None
    
    .. _ax.tick_params(): https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.tick_params.html
    """
    ############
    # get axis #
    ############
    if ax is None:
        ax = plt.gca()

    ######################
    # default parameters #
    ######################
    if 'axis' not in kwargs:
        kwargs['axis'] = 'y'
    if 'which' not in kwargs:
        kwargs['which'] = 'both'
    # if 'direction' not in kwargs:
    #     kwargs['direction'] = 'in'
    # if 'left' not in kwargs:
    #     kwargs['left'] = True
    # if 'right' not in kwargs:
    #     kwargs['right'] = True

    #####################
    # major ticks setup #
    #####################
    if ticks is None:
        if nticks is not None:
            assert nticks > 1, 'nticks needs to be > 1.'
        if start is not None or stop is not None or nticks is not None or step is not None:
            if start is not None and stop is not None and nticks is not None:
                ticks = np.linspace(start, stop, nticks)
            elif start is not None and stop is not None and step is not None:
                ticks   = np.arange(start, stop + step*0.1, step) 
            elif start is not None and nticks is not None and step is not None:
                stop = start + nticks*step
                ticks = np.linspace(start, stop, nticks)
            elif nticks is not None:
                ticks_showing = get_yticks_showing(ax=ax)
                start = ticks_showing[0]
                stop  = ticks_showing[-1]
                ticks = np.linspace(start, stop, nticks)
            else:
                raise ValueError('Invalid input. See help(set_yticks) for help.')
            
            # set ticks
            _ = ax.yaxis.set_ticks(ticks)

            # set label font
            if fontproperties is not None:
                _ = ax.set_yticklabels([str(i) for i in ticks], fontproperties=fontproperties, visible=True)
            
            # if start/stop is defined, select pad to (0, 0)
            # if pad is None:
            #     pad = (0, 0)
    else:
        assert isinstance(ticks, Iterable), 'ticks must be iterable'
        assert len(ticks) > 1, 'ticks must have length > 1'
        _ = ax.yaxis.set_ticks(ticks)

        # set labels
        if labels is not None:
            assert isinstance(labels, Iterable), 'ticks must be iterable'
            assert len(ticks) == len(labels), 'len(ticks) must be the same as len(labels)'
            if fontproperties is not None:
                _ = ax.set_yticklabels(labels, fontproperties=fontproperties, visible=True)
            else:
                _ = ax.set_yticklabels(labels)
        else:
            if fontproperties is not None:
                _ = ax.set_yticklabels([str(i) for i in ticks], fontproperties=fontproperties, visible=True)

    ##########
    # limits #
    ##########
    if pad is not None:
        if isinstance(pad, Iterable) == False:
            assert is_number(pad), 'pad must be a number or a tuple/list with length 2'
            pad = (pad, pad)

        # check
        assert len(pad) == 2, 'pad must be a number or a tuple/list with length 2'
        assert is_number(pad[0]), 'pad must be a number or a tuple/list with length 2'
        assert is_number(pad[1]), 'pad must be a number or a tuple/list with length 2'

        # set limits
        step = ticks[1]  - ticks[0]
        min_lim   = ticks[0]  - step*pad[0]
        max_lim   = ticks[-1] + step*pad[1]
        ax.set_ylim((min_lim, max_lim), auto=False)

    ###############
    # minor ticks #
    ###############
    if minor_ticks is None:
        if n_minor_ticks is not None:
            ax.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
        elif minor_step is not None:
            step     = ticks[1] - ticks[0]
            n_minor_ticks = int(round(step/minor_step))
            ax.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
    else:
        assert isinstance(ticks, Iterable), 'minor_ticks must be iterable'
        assert len(ticks) > 0, 'minor_ticks must have length > 0'
        _ = ax.yaxis.set_ticks(minor_ticks, minor=True)

    ###################
    # tick parameters #
    ###################
    labelparams = ['labelsize', 'labelcolor', 'labelfontfamily', 'labelbottom', 'labeltop', 'labelleft', 'labelright', 'labelrotation']
    labelparams = {param:kwargs[param] for param in labelparams if param in kwargs}
    _ = [kwargs.pop(param) for param in labelparams]
    ax.tick_params(**labelparams)

    for _ax in ax.get_shared_y_axes().get_siblings(ax):
        _ax.tick_params(**kwargs)
   
    return

# remove tick labels
def remove_xticklabels(ax=None):
    """remove/hide x tick labels from axes
    
    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.

    Returns:
        None

    See Also:
        :py:func:`remove_yticklabels`
        :py:func:`hide_xticklabels`
        :py:func:`hide_yticklabels`
        :py:func:`show_ticklabels`
    """
    if ax is None:
        ax = plt

    return ax.tick_params(labelbottom=False, labeltop=False)  

def remove_yticklabels(ax=None):
    """remove/hide y tick labels from axes
    
    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.

    Returns:
        None
    
    See Also:
        :py:func:`remove_xticklabels`
        :py:func:`hide_xticklabels`
        :py:func:`hide_yticklabels`
        :py:func:`show_ticklabels`
    """
    if ax is None:
        ax = plt

    return ax.tick_params(labelleft=False, labelright=False)  

def hide_xticklabels(ax=None):
    """remove/hide x tick labels from axes
    
    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.

    Returns:
        None
    
    See Also:
        :py:func:`remove_xticklabels`
        :py:func:`remove_yticklabels`
        :py:func:`hide_yticklabels`
        :py:func:`show_ticklabels`
    """
    return remove_xticklabels(ax=ax)

def hide_yticklabels(ax=None):
    """remove/hide y tick labels from axes
    
    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.

    Returns:
        None

    See Also:
        :py:func:`remove_xticklabels`
        :py:func:`remove_yticklabels`
        :py:func:`hide_xticklabels`
        :py:func:`show_ticklabels`
    """
    return remove_yticklabels(ax=ax)

def show_ticklabels(ax=None, left=None, right=None, bottom=None, top=None):
    """show tick labels in axes
    
    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.
        left, right, bottom, top (bool or None, optional): if True (False), tick
          labels will be shown (hidden). If None, status won't change. Default is None.

    Returns:
        None
    
    See Also:
        :py:func:`remove_xticklabels`
        :py:func:`remove_yticklabels`
        :py:func:`hide_xticklabels`
        :py:func:`hide_yticklabels`
    """
    if ax is None:
        ax = plt

    _kwargs = {'left':left, 'right':right, 'bottom':bottom, 'top':top}
    kwargs  = {}
    for name in _kwargs:
        if _kwargs[name] is not None:
            kwargs[name] = _kwargs[name]

    return ax.tick_params(**kwargs)  
 


# %% add to axe
def _get_xticks_showing(self):
    """get x ticks that are sowing in a plot

    Returns:
        list
    """
    return get_xticks_showing(self)
mpl.axes.Axes.get_xticks_showing = _get_xticks_showing

def _get_yticks_showing(self):
    """get y ticks that are sowing in a plot

    Returns:
        list
    """
    return get_yticks_showing(self)
mpl.axes.Axes.get_yticks_showing = _get_yticks_showing



def _remove_xticklabels(self):
    """remove/hide x tick labels from axes
    
    Returns:
        None
    """
    return remove_xticklabels(self)
mpl.axes.Axes.remove_xticklabels = _remove_xticklabels

def _remove_yticklabels(self):
    """remove/hide y tick labels from axes
    
    Returns:
        None
    """
    return remove_yticklabels(self)
mpl.axes.Axes.remove_yticklabels = _remove_yticklabels

def _hide_xticklabels(self):
    """remove/hide x tick labels from axes
    
    Returns:
        None
    """
    return hide_xticklabels(self)
mpl.axes.Axes.hide_xticklabels = _hide_xticklabels

def _hide_yticklabels(self):
    """remove/hide y tick labels from axes
    
    Returns:
        None
    """
    return hide_yticklabels(self)
mpl.axes.Axes.hide_yticklabels = _hide_yticklabels

def _show_ticklabels(self, left=None, right=None, bottom=None, top=None):
    """show tick labels in axes
    
    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.
        left, right, bottom, top (bool or None, optional): if True (False), tick
          labels will be shown (hidden). If None, status won't change. Default is None.    

    Returns:
        None
    """
    return show_ticklabels(ax=self, left=left, right=right, bottom=bottom, top=top)
mpl.axes.Axes.show_ticklabels = _show_ticklabels



def _set_xticks(self, start=None, stop=None, nticks=None, step=None, ticks=None, labels=None, pad=None, n_minor_ticks=None, minor_step=None, minor_ticks=None, fontproperties=None, **kwargs):
    """Set x ticks of a plot.

    Usage:
        For setting up major ticks:

            set_ticks(start, stop, nticks)
            set_ticks(start, stop, step)
            set_ticks(start, nticks, step)  # `stop` is calculated
            set_ticks(nticks)                    # `start` and `stop` are kept from the plot
            set_ticks(ticks)                      # manually define each tick position
            set_ticks(ticks, labels)              # manually define each tick position and label

        For setting up minor ticks:

            set_ticks(n_minor_ticks)             # equally spaced minor ticks
            set_ticks(minor_step)           # equally spaced minor ticks
            set_ticks(minor_ticks)               # manually define each minor tick position

    Args:
        *Major ticks*
        start, stop (number, optional): major ticks start/stop value. Use None
            to not keep current start/stop value. Default is None.
        nticks, step (int/number, optional): Number of major ticks or separation between ticks. 
            nticks overwrites step. Default is None.
        ticks, labels (Iterable, optional): Iterable with tick positions and labels.
        Overwrites start, stop, nticks, step.
        
        *Padding (plotting limits)*
        pad (number or tuple, optional): padding between plot edge and 
            closest tick in terms of ticks separation. If tuple, first value 
            applies to the left edge and second value to the right edge. 
            Typically, padding should be between a value between 0 and 1. If pad
            is None, plot limits remain unchanged. However, if pad is None and 
            start/stop are defined, pad is set to (0, 0). Default is None.

        *Minor ticks*
        n_minor_ticks, minor_step (int/number, optional): Number of minor 
        ticks between two major ticks or minor tick separation. n_minor_ticks
        overwrites minor_step. Default is None.
        minor_ticks (Iterable, optional): Iterable with minor tick positions.
        Overwrites n_minor_ticks and minor_step.

        *font*
        fontproperties: Label ticks font. Use ``matplotlib.font_manager.FontProperties``.

        *ticks parameters*
        **kwargs: kwargs are passed to `ax.tick_params()`_ that sets tick parameters.

        If not specified, the following parameters are passed to `ax.tick_params()`_:

        Args:
            which ('major', 'minor', 'both', optional): The group of ticks to which the parameters 
                are applied. Default is 'both'.
            direction ('in', 'out', 'inout', optional): Puts ticks inside the 
                Axes, outside the Axes, or both. Default is 'in'
            bottom, top (bool, optional): Whether to draw the bottom/top ticks.
                Default is True.

    Returns:
        None
    
    .. _ax.tick_params(): https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.tick_params.html
    """
    return set_xticks(ax=self, start=start, stop=stop, nticks=nticks, step=step, ticks=ticks, labels=labels, pad=pad, n_minor_ticks=n_minor_ticks, minor_step=minor_step, minor_ticks=minor_ticks, fontproperties=fontproperties, **kwargs)
mpl.axes.Axes.xticks = _set_xticks

def _set_yticks(self, start=None, stop=None, nticks=None, step=None, ticks=None, labels=None, pad=None, n_minor_ticks=None, minor_step=None, minor_ticks=None, fontproperties=None, **kwargs):
    """Set y ticks of a plot.

    Usage:
        For setting up major ticks:

            set_ticks(start, stop, nticks)
            set_ticks(start, stop, step)
            set_ticks(start, nticks, step)  # `stop` is calculated
            set_ticks(nticks)                    # `start` and `stop` are kept from the plot
            set_ticks(ticks)                      # manually define each tick position
            set_ticks(ticks, labels)              # manually define each tick position and label

        For setting up minor ticks:

            set_ticks(n_minor_ticks)             # equally spaced minor ticks
            set_ticks(minor_step)           # equally spaced minor ticks
            set_ticks(minor_ticks)               # manually define each minor tick position

    Args:
        *Major ticks*
        start, stop (number, optional): major ticks start/stop value. Use None
            to not keep current start/stop value. Default is None.
        nticks, step (int/number, optional): Number of major ticks or separation between ticks. 
            nticks overwrites step. Default is None.
        ticks, labels (Iterable, optional): Iterable with tick positions and labels.
        Overwrites start, stop, nticks, step.
        
        *Padding (plotting limits)*
        pad (number or tuple, optional): padding between plot edge and 
            closest tick in terms of ticks separation. If tuple, first value 
            applies to the left edge and second value to the right edge. 
            Typically, padding should be between a value between 0 and 1. If pad
            is None, plot limits remain unchanged. However, if pad is None and 
            start/stop are defined, pad is set to (0, 0). Default is None.

        *Minor ticks*
        n_minor_ticks, minor_step (int/number, optional): Number of minor 
        ticks between two major ticks or minor tick separation. n_minor_ticks
        overwrites minor_step. Default is None.
        minor_ticks (Iterable, optional): Iterable with minor tick positions.
        Overwrites n_minor_ticks and minor_step.

        *font*
        fontproperties: Label ticks font. Use ``matplotlib.font_manager.FontProperties``.

        *ticks parameters*
        **kwargs: kwargs are passed to `ax.tick_params()`_ that sets tick parameters.

        If not specified, the following parameters are passed to `ax.tick_params()`_:

        Args:
            which ('major', 'minor', 'both', optional): The group of ticks to which the parameters 
                are applied. Default is 'both'.
            direction ('in', 'out', 'inout', optional): Puts ticks inside the 
                Axes, outside the Axes, or both. Default is 'in'
            left, right (bool, optional): Whether to draw the left/right ticks.
                Default is True.

    Returns:
        None
    
    .. _ax.tick_params(): https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.tick_params.html
    """
    return set_yticks(ax=self, start=start, stop=stop, nticks=nticks, step=step, ticks=ticks, labels=labels, pad=pad, n_minor_ticks=n_minor_ticks, minor_step=minor_step, minor_ticks=minor_ticks, fontproperties=fontproperties, **kwargs)
mpl.axes.Axes.yticks = _set_yticks
# %%

# %% ============================= legend ================================= %% #
def leg(*args, **kwargs):
    """Place a legend on the Axes. Wrapper for `plt.legend()`_

    Args:
        *args, **kwargs: args and kwargs passed to `plt.legend()`_
        ax (axes, optional): axes. If None, the current axes will be used. Default is None
        coord (tuple): (x, y) coordinates in data coordinates to place legend. 
            (x, y) places the corner of the legend specified by `loc` at x, y.
            If None, this is ignored. Default is None. Overwrites `bbox_to_anchor`.
        
    If not specified, the following parameters are passed to `ax.legend()`_:

    Args:
        labelspacing (float, optional): The vertical space between the legend 
            entries, in font-size units. Default is 0.5
        frameon (bool, optional): Whether the legend should be drawn on a patch 
            (frame). Default is False

    Returns:
        legend object

    .. _plt.legend(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
    """
    # get axes
    ax = None
    if 'ax' in kwargs:
        ax = kwargs.pop('ax')
    if ax is None:
        ax = plt.gca()

    # default attrs
    if 'labelspacing' not in kwargs:
        kwargs['labelspacing'] = 0.1
    if 'frameon' not in kwargs:
        kwargs['frameon'] = False

    # get coord
    coord = None
    if 'coord' in kwargs:
        coord = kwargs.pop('coord')

    # legend
    if coord is not None:   
        # check
        assert isinstance(coord, Iterable), 'coord must be a tuple of type (x, y)'
        assert len(coord) == 2, 'coord must be a tuple of type (x, y)'

        # covert to Data coordinates to axes coordinates
        _x, _y = ax.transData.transform(coord)
        x, y = ax.transAxes.inverted().transform((_x, _y)) 

        if 'loc' not in kwargs:
            kwargs['loc'] = "upper left"
        if 'bbox_to_anchor' not in kwargs:
            kwargs['bbox_to_anchor'] = (x, y)

    return ax.legend(*args, **kwargs)

def _leg(self, *args, **kwargs):
    """Place a legend on the Axes. Wrapper for `plt.legend()`_

    Args:
        *args, **kwargs: args and kwargs passed to `plt.legend()`_
        coord (tuple): (x, y) coordinates in data coordinates to place legend. 
            (x, y) places the corner of the legend specified by `loc` at x, y.
            If None, this is ignored. Default is None. Overwrites `bbox_to_anchor`.
        
    If not specified, the following parameters are passed to `ax.legend()`_:

    Args:
        labelspacing (float, optional): The vertical space between the legend 
            entries, in font-size units. Default is 0.5
        frameon (bool, optional): Whether the legend should be drawn on a patch 
            (frame). Default is False

    Returns:
        legend object

    .. _plt.legend(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
    """
    kwargs['ax'] = self
    leg(*args, **kwargs)
mpl.axes.Axes.leg = _leg
# %%

# %% ================================ rectangle =========================== %% #
def rectangle(xlim, ylim, ax=None, **kwargs):
    """plot a rectangle. Wrapper for `plt.Rectangle()`_.

    Args:
        xlim, ylim (tuple): xmin, xmax, ymin, ymax coordinates of the rectangle
        in data coordinates.
        ax (axes, optional): axes. If None, the current axes will be used. Default is None
        **kwargs: kwargs are passed to `plt.Rectangle()`_

    If not specified, the following parameters are passed to `ax.tick_params()`_:

    Args:
        linewidth (number, optional): Set the patch linewidth in points. 
            Default is 0.6.
        facecolor (color or None, optional): Set the patch face color. 
            Default is 'none'
        edgecolor (color or None, optional): Set the patch edge color. 
            Default is gray.
        alpha (number, optional): Set the alpha value used for blending 
            - not supported on all backends. Default is 0.5
        zorder (number, optional): Set the zorder for the artist. 
            Artists with lower zorder values are drawn first. Default is such 
            that rectangle is ploted on top of all other curves.
        
    Returns:
        rectangle object 
    
    .. _plt.Rectangle(): https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Rectangle.html#matplotlib.patches.Rectangle
    """
    if ax is None:
        ax = plt.gca()

    if 'linewidth' not in kwargs and 'lw' not in kwargs:
            kwargs['linewidth'] = 0.6
    if 'facecolor' not in kwargs and 'fc' not in kwargs:
        kwargs['fc'] = 'none'
    if 'edgecolor' not in kwargs:
        kwargs['edgecolor'] = 'gray'
    if 'alpha' not in kwargs:
        kwargs['alpha'] = 0.5
    if 'zorder' not in kwargs:
        kwargs['zorder'] = max([_.zorder for _ in ax.get_children()])

    rect = plt.Rectangle((xlim[0], ylim[0]), xlim[1]-xlim[0], ylim[1]-ylim[0], **kwargs)
    ax.add_patch(rect)

    return rect

def _rectangle(self, xlim, ylim, **kwargs):
    """plot a rectangle. Wrapper for `plt.Rectangle()`_.

    Args:
        xlim, ylim (tuple): xmin, xmax, ymin, ymax coordinates of the rectangle
        in data coordinates.
        **kwargs: kwargs are passed to `plt.Rectangle()`_

    If not specified, the following parameters are passed to `ax.tick_params()`_:

    Args:
        linewidth (number, optional): Set the patch linewidth in points. 
            Default is 0.6.
        facecolor (color or None, optional): Set the patch face color. 
            Default is 'none'
        edgecolor (color or None, optional): Set the patch edge color. 
            Default is gray.
        alpha (number, optional): Set the alpha value used for blending 
            - not supported on all backends. Default is 0.5
        zorder (number, optional): Set the zorder for the artist. 
            Artists with lower zorder values are drawn first. Default is such 
            that rectangle is ploted on top of all other curves.
        
    Returns:
        rectangle object 
    
    .. _plt.Rectangle(): https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Rectangle.html#matplotlib.patches.Rectangle
    """
    return rectangle(xlim=xlim, ylim=ylim, ax=self, **kwargs)
mpl.axes.Axes.rectangle = _rectangle
# %%

# %% ================================== zoom ============================== %% #
def zoom(start, stop, ax=None, ymargin=5):
    """Zoom up portion of current figure from start to stop.

    Args:
        start, stop (float): initial and final x value.
        ax (axes, optional): axes to apply zoom. If None, the current axes 
        will be used. Default is None
        ymargin (number, optional): margin value between data and the edges of
            plot in percentage of the y data range.
    """
    if ax is None:
        ax = plt.gca()

    ymax = None
    ymin = None

    for line in ax.get_lines():

        x = line.get_data()[0]
        y = line.get_data()[1]

        try:
            _, y = _extract(x=x, y=y, ranges=(start, stop))
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
        except RuntimeError:
            pass

    if ymax is None:
        ax.set_xlim(start, stop)
    else:
        m = (ymax-ymin)*ymargin/100
        ax.set_ylim(ymin-m, ymax+m)
        ax.set_xlim(start, stop)
    return

def _zoom(self, start, stop, ymargin=5):
    """Zoom up portion of current figure from start to stop.

    Args:
        start, stop (float): initial and final x value.
        will be used. Default is None
        ymargin (number, optional): margin value between data and the edges of
            plot in percentage of the y data range.
    """
    zoom(start=start, stop=stop, ax=self, ymargin=ymargin)
mpl.axes.Axes.zoom = _zoom
# %%

# %% ================================= inset ============================== %% #
def inset(xlim, ylim, ax=None, xticks_kwargs=None, yticks_kwargs=None, rect=True, rect_xlim=None, rect_ylim=None, ax2put_inset=None, **kwargs):
    """add inset to axes. Inset will have all curves from axes

    Args:
        xlim, ylim (tuple): (xi, xf), (yi, yf) coordinates of inset in data 
            coordinates.
        ax (axes, optional): axes to get inset data from. If None, the current axes 
            will be used. Default is None
        xticks_kwargs, yticks_kwargs (dictionary, optional): kwargs for 
            set_xticks() and set_yticks() functions to set the inset ticks
            and plot limits. Default is None.
        rect (bool, optional): if True, rectangle will be draw in the main axes
            at the plotting limits defined by xticks_kwargs and yticks_kwargs. 
        rect_xlim, rect_ylim (tuple, optional): if defines the rectangle edges.
            overwrites limits defined by xticks_kwargs, yticks_kwargs. 
        ax2put_inset (axex, optional): axes to put inset in. If None, the current axes 
            will be used. Default is None

    Returns:    
        inset axes object
    """
    if ax is None:
        ax = plt.gca()
    if ax2put_inset is None:
        ax2put_inset = ax

    x_init      = xlim[0]
    x_final     = xlim[1]
    y_init      = ylim[0]
    y_final     = ylim[1]

    # create inset
    inset = ax.get_figure().add_axes(_ax_box2fig_box(ax2put_inset, [x_init, y_init, x_final, y_final]))

    for line in ax2put_inset.get_lines():
        line2, = inset.plot(line.get_xdata(), line.get_ydata(),
                            ls=line.get_linestyle(),
                            linewidth=line.get_linewidth(),
                            marker = line.get_marker(),
                            markevery=line.get_markevery(),
                            ms=line.get_markersize(),
                            color=line.get_color(),
                            alpha=line.get_alpha(),
                            markerfacecolor=line.get_markerfacecolor(),
                            markeredgewidth=line.get_markeredgewidth(),
                            label=line.get_label())
    

    if xticks_kwargs is not None:
        set_xticks(ax=inset, **xticks_kwargs)
    if yticks_kwargs is not None:
        set_yticks(ax=inset, **yticks_kwargs)

    # Rectangle
    if rect:
        if 'linewidth' not in kwargs or 'lw' not in kwargs:
            kwargs['linewidth'] = 0.6
        if 'facecolor' not in kwargs or 'fc' not in kwargs:
            kwargs['fc'] = 'none'
        if 'edgecolor' not in kwargs:
            kwargs['edgecolor'] = 'gray'
        if 'alpha' not in kwargs:
            kwargs['alpha'] = 0.5
        if 'zorder' not in kwargs:
            kwargs['zorder'] = max([_.zorder for _ in ax.get_children()])

        if rect_xlim is None:
            rect_xlim = (inset.get_xlim()[0], inset.get_xlim()[1]-inset.get_xlim()[0])
            print(f'rect_xlim={rect_xlim}')
        if rect_ylim is None:
            rect_ylim = (inset.get_ylim()[0],  inset.get_ylim()[1]-inset.get_ylim()[0])
            print(f'rect_ylim={rect_ylim}')
        
        _ = ax2put_inset.spines['left'].get_linewidth()
        rect = plt.Rectangle((rect_xlim[0], rect_ylim[0]), rect_xlim[1], rect_ylim[1], **kwargs)
        ax2put_inset.add_patch(rect)

    return inset

def _ax_box2fig_box(ax, points):
    """Transform 'bbox like' axis position values to percentage fig position.

    Useful for positioning insets.

    Args:
        ax (matplotlib.axes.Axes): axes instance.
        points (list): list like ``[x_init, y_init, x_final, y_final]``

    Returns:
        'bbox like' figure positions ``[x_init, y_init, delta_x, delta_y]``.
    """
    [x_init, y_init, x_final, y_final] = points

    x_init  = _ax_pos2fig_pos(ax, x_init, direction='x')
    y_init  = _ax_pos2fig_pos(ax, y_init, direction='y')
    delta_x = _ax_pos2fig_pos(ax, x_final, direction='x') - x_init
    delta_y = _ax_pos2fig_pos(ax, y_final, direction='y') - y_init

    return [x_init, y_init, delta_x, delta_y]

def _ax_pos2fig_pos(ax, value, direction='x'):
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

def _inset(self, xlim, ylim, xticks_kwargs=None, yticks_kwargs=None, rect=True, rect_xlim=None, rect_ylim=None, ax2putinset=None, **kwargs):
    """add inset to axes. Inset will have all curves from axes

    Args:
        xlim, ylim (tuple): (xi, xf), (yi, yf) coordinates of inset in data 
            coordinates.
        xticks_kwargs, yticks_kwargs (dictionary, optional): kwargs for 
            set_xticks() and set_yticks() functions to set the inset ticks
            and plot limits. Default is None.
        rect (bool, optional): if True, rectangle will be draw in the main axes
            at the plotting limits defined by xticks_kwargs and yticks_kwargs. 
        rect_xlim, rect_ylim (tuple, optional): if defines the rectangle edges.
            overwrites limits defined by xticks_kwargs, yticks_kwargs. 
        ax2putinset (axex, optional): puts inset in a different axes.

    Returns:    
        inset axes object
    """
    return inset(xlim=xlim, ylim=ylim, ax=self, xticks_kwargs=xticks_kwargs, yticks_kwargs=yticks_kwargs, rect=rect, rect_xlim=rect_xlim, rect_ylim=rect_ylim, ax2putinset=ax2putinset, **kwargs)
mpl.axes.Axes.inset = _inset
# %%

# %% ======================== units and conversion ======================== %% #
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

def rgb2hex(rgb, max_rgb_value=1):
    """Return hex value.

    Args:
        rgb (tuple): rgb values ([255 255 255])
        max_rgb_value (int, optional): max rgb possible value. Tipically, 1 or 255.

    Returns:
        hex value
    """
    if max(rgb) > max_rgb_value:
        raise ValueError(f'rgb={rgb} maximum value is bigger than max_rgb_value.')
    return mpl.colors.to_hex([x/max_rgb_value for x in rgb])

def hex2rgb(string, max_rgb_value=1):
    """Return rgb value.

    Args:
        string (string): hex value as string.
        max_rgb_value (int, optional): max rgb possible value. Tipically, 1 or 255.

    Returns:
        rgb value (tuple)
    """
    return [x*max_rgb_value for x in mpl.colors.to_rgb(string)]
# %%

# %% ================================= lines ============================== %% #
def axvlines(x, ymin=None, ymax=None, ax=None, **kwargs):
    """draw vertical lines. Wrapper for `plt.vlines()`_ and `plt.axvline()`_

    Args:
        x (float or array): x coordinates where to plot the lines.
        ymin, ymax (float, optional): Respective beginning and end of 
            each line. If None, infinite lines are plotted. Default is None.
        color, colors, or c (str, number, or list, optional): if str or number, this color will be 
                applied to every spectra. If list, it must have the same length 
                as the number of spectra and each element must be a color (a 
                color can be a str or a 3 element (RGB) list.
                Default is None. 
        linestyles, linestyle, or ls (string or list, optional): linestyles. Default is '--'.
        label or labels (string or list, optional): labels. Default is ''.
            If labels is str, only first line gets labeled. If labels is a list
            each line is labeled independtly.
        ax (axes, optional): axes. If None, the current axes will be used. Default is None
        **kwargs: kwargs are passed to `plt.vlines()`_ or `plt.axvline()`_

    Returns:
        list with Line2D objects.
        
    .. _plt.axvline(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axvline.html  
    .. _plt.vlines(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.vlines.html
    """
    # get axis
    if ax is None:
        ax = plt.gca()

    # fix attrs
    if isinstance(x, Iterable) == False:
        x = [x, ]

    #########
    # label #
    #########
    if 'label' in kwargs or 'labels' in kwargs:
        if 'label' in kwargs:
            _label = kwargs.pop('label')
        else:
            _label = kwargs.pop('labels')

        if isinstance(_label, Iterable) == True and isinstance(_label, str) == False:
            assert len(_label) == len(x), f'label must be a number of a list with length compatible with the number of lines.\nnumber of labels: {len(label)}\nnumber of spectra: {len(x)}'
            label = _label
        else:
            label = [None]*len(x)
            label[0] = _label
    else:
        label = [None]*len(x)

    #########
    # color #
    #########
    if 'color' in kwargs or 'colors' in kwargs or 'c' in kwargs:
        if 'color' in kwargs:
            _color = kwargs.pop('color')
        elif 'colors' in kwargs:
            _color = kwargs.pop('colors')
        else:
            _color = kwargs.pop('c')

        if isinstance(_color, Iterable) == True and isinstance(_color, str) == False:
            assert len(_color) == len(x), f'`color` must be a number of a list with length compatible with the number of lines.\nnumber of colors: {len(_color)}\nnumber of lines: {len(x)}'
            colors = _color
        else:
            colors = [_color]*len(x)
    else:
        colors = [None]*len(x)

    ##############
    # linestyles #
    ##############
    if 'linestyles' in kwargs or 'linestyle' in kwargs or 'ls' in kwargs:
        if 'linestyles' in kwargs:
            _linestyle = kwargs.pop('linestyles')
        elif 'linestyle' in kwargs:
            _linestyle = kwargs.pop('linestyle')
        else:
            _linestyle = kwargs.pop('ls')

        if isinstance(_linestyle, Iterable) == True and isinstance(_linestyle, str) == False:
            assert len(_linestyle) == len(x), f'`linestyles` must be a number of a list with length compatible with the number of lines.\nnumber of linestyles: {len(_linestyle)}\nnumber of lines: {len(x)}'
            linestyles = _linestyle
        else:
            linestyles = [_linestyle]*len(x)
    else:
        linestyles = ['--']*len(x)
        
    # plot
    final = []
    if ymin is None and ymax is None:
        for i, _x in enumerate(x):
            final.append(ax.axvline(x=_x, color=colors[i], linestyle=linestyles[i], label=label[i], **kwargs))
    else:
        for i, _x in enumerate(x):
            final.append(ax.vlines(x=_x, ymin=ymin, ymax=ymax, colors=colors[i], linestyles=linestyles[i], label=label[i], **kwargs))
    return final

def axhlines(y, xmin=None, xmax=None, ax=None, **kwargs):
    """draw horizontal lines. Wrapper for `plt.hlines()`_ and `plt.axhline()`_

    Args:
        y (float or array): y coordinates where to plot the lines.
        xmin, xmax (float, optional): Respective beginning and end of 
            each line. If None, infinite lines are plotted. Default is None.
        color, colors, or c (str, number, or list, optional): if str or number, this color will be 
                applied to every spectra. If list, it must have the same length 
                as the number of spectra and each element must be a color (a 
                color can be a str or a 3 element (RGB) list.
                Default is None. 
        linestyles, linestyle, or ls (string or list, optional): linestyles. Default is '--'.
        label or labels (string or list, optional): labels. Default is ''.
            If labels is str, only first line gets labeled. If labels is a list
            each line is labeled independtly.
        ax (axes, optional): axes. If None, the current axes will be used. Default is None
        **kwargs: kwargs are passed to `plt.hlines()`_ or `plt.axhline()`_

    Returns:
        list with Line2D objects.
        
    .. _plt.axhline(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axhline.html  
    .. _plt.hlines(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hlines.html
    """
    # get axis
    if ax is None:
        ax = plt.gca()

    # fix attrs
    if isinstance(y, Iterable) == False:
        y = [y, ]

    #########
    # label #
    #########
    if 'label' in kwargs or 'labels' in kwargs:
        if 'label' in kwargs:
            _label = kwargs.pop('label')
        else:
            _label = kwargs.pop('labels')

        if isinstance(_label, Iterable) == True and isinstance(_label, str) == False:
            assert len(_label) == len(y), f'label must be a number of a list with length compatible with the number of lines.\nnumber of labels: {len(label)}\nnumber of spectra: {len(y)}'
            label = _label
        else:
            label = [None]*len(y)
            label[0] = _label
    else:
        label = [None]*len(y)

    #########
    # color #
    #########
    if 'color' in kwargs or 'colors' in kwargs or 'c' in kwargs:
        if 'color' in kwargs:
            _color = kwargs.pop('color')
        elif 'colors' in kwargs:
            _color = kwargs.pop('colors')
        else:
            _color = kwargs.pop('c')

        if isinstance(_color, Iterable) == True and isinstance(_color, str) == False:
            assert len(_color) == len(y), f'`color` must be a number of a list with length compatible with the number of lines.\nnumber of colors: {len(_color)}\nnumber of lines: {len(y)}'
            colors = _color
        else:
            colors = [_color]*len(y)
    else:
        colors = [None]*len(y)

    ##############
    # linestyles #
    ##############
    if 'linestyles' in kwargs or 'linestyle' in kwargs or 'ls' in kwargs:
        if 'linestyles' in kwargs:
            _linestyle = kwargs.pop('linestyles')
        elif 'linestyle' in kwargs:
            _linestyle = kwargs.pop('linestyle')
        else:
            _linestyle = kwargs.pop('ls')

        if isinstance(_linestyle, Iterable) == True and isinstance(_linestyle, str) == False:
            assert len(_linestyle) == len(y), f'`linestyles` must be a number of a list with length compatible with the number of lines.\nnumber of linestyles: {len(_linestyle)}\nnumber of lines: {len(y)}'
            linestyles = _linestyle
        else:
            linestyles = [_linestyle]*len(y)
    else:
        linestyles = ['--']*len(y)
        
    # plot
    final = []
    if xmin is None and xmax is None:
        for i, _y in enumerate(y):
            final.append(ax.axhline(y=_y, color=colors[i], linestyle=linestyles[i], label=label[i], **kwargs))
    else:
        for i, _y in enumerate(y):
            final.append(ax.hlines(y=_y, xmin=xmin, xmax=xmax, colors=colors[i], linestyles=linestyles[i], label=label[i], **kwargs))
    return final

def _axvlines(self, x, ymin=None, ymax=None, **kwargs):
    """draw vertical lines. Wrapper for `plt.vlines()`_ and `plt.axvline()`_

    Args:
        x (float or array): x coordinates where to plot the lines.
        ymin, ymax (float, optional): Respective beginning and end of 
            each line. If None, infinite lines are plotted. Default is None.
        color, colors, or c (str, number, or list, optional): if str or number, this color will be 
                applied to every spectra. If list, it must have the same length 
                as the number of spectra and each element must be a color (a 
                color can be a str or a 3 element (RGB) list.
                Default is None. 
        linestyles, linestyle, or ls (string or list, optional): linestyles. Default is '--'.
        label or labels (string or list, optional): labels. Default is ''.
            If labels is str, only first line gets labeled. If labels is a list
            each line is labeled independtly.
        **kwargs: kwargs are passed to `plt.vlines()`_ or `plt.axvline()`_

    Returns:
        list with Line2D objects.
        
    .. _plt.axvline(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axvline.html  
    .. _plt.vlines(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.vlines.html
    """
    return axvlines(x=x, ymin=ymin, ymax=ymax, ax=self, **kwargs)
mpl.axes.Axes.axvlines = _axvlines

def _axhlines(self, y, xmin=None, xmax=None, **kwargs):
    """draw horizontal lines. Wrapper for `plt.hlines()`_ and `plt.axhline()`_

    Args:
        y (float or array): y coordinates where to plot the lines.
        xmin, xmax (float, optional): Respective beginning and end of 
            each line. If None, infinite lines are plotted. Default is None.
        color, colors, or c (str, number, or list, optional): if str or number, this color will be 
                applied to every spectra. If list, it must have the same length 
                as the number of spectra and each element must be a color (a 
                color can be a str or a 3 element (RGB) list.
                Default is None. 
        linestyles, linestyle, or ls (string or list, optional): linestyles. Default is '--'.
        label or labels (string or list, optional): labels. Default is ''.
            If labels is str, only first line gets labeled. If labels is a list
            each line is labeled independtly.
        ax (axes, optional): axes. If None, the current axes will be used. Default is None
        **kwargs: kwargs are passed to `plt.hlines()`_ or `plt.axhline()`_

    Returns:
        list with Line2D objects.
        
    .. _plt.axhline(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axhline.html  
    .. _plt.hlines(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hlines.html
    """
    return axhlines(y=y, xmin=xmin, xmax=xmax, ax=self, **kwargs)
mpl.axes.Axes.axhlines = _axhlines
# %%

# %% =============================== Gradient ============================= %% #
# linear
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

# benzier
fact_cache = {} # Value cache
def _factorial(n):
    """Memoized factorial function used by the _bernstein() function."""
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
# %%

# %% =============================== Experimental ========================= %% #
def _onmove(event):
    """This function is called every time the figure changes position on the screen."""
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
    return

def _bring2top():
    """Brings current window to the top.

    This function was not tested for all available matplotlib backends.
    """
    # fig = plt.gcf()
    # fig.canvas.manager.window.raise_()

    backend = matplotlib.get_backend()
    if backend.startswith(('Qt5', 'qt5', 'QT5')):
        print('ff')

        from PyQt5 import QtCore
        window = plt._get_current_fig_manager().window
        window.activateWindow()
        window.raise_()
        # window.setWindowFlags(window.windowFlags() | QtCore.Qt.WindowStaysOnTopHint)
        # plt.show()

        # window.setWindowFlags(window.windowFlags() & ~QtCore.Qt.WindowStaysOnTopHint)
        # plt.show()
    elif backend.startswith(('tk', 'TK', 'Tk')):
        plt.gcf().canvas.get_tk_widget().focus_force() 
    else:
        figManager = _get_current_fig_manager()
        figManager.window.raise_()
    return

# %% =============================== Obsolete ============================= %% #
def _remove_ticks_edge(ax):
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

def _set_ticks_old(ax=None, axis='x', autoscale=True, **kwargs):
    """Set ticks of a plot.

    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.
        axis (string, optional): possible values are 'x' or 'y'.
        start (float or int): start value for ticks (not the plot edge --- see ``pad``)
        stop (float or int): stop value for ticks (not the plot edge --- see ``pad``)
        nticks (int): Number of ticks. Ticks separation is calculated accordingly and this parameter overwrites step.
        step (float or int): Ticks separation.
        n_minor_ticks (int): Number of minor ticks between two major ticks.
        fontproperties: Label ticks font. Use ``matplotlib.font_manager.FontProperties``.
        pad (float or int): Distance between plot edge and the first tick in terms of tick separation. Typically, must be something between 0 and 1.
        n_decimal_places (int): Number of decimal places to use for tick labels.
        direction (str, optional): default is 'out', possible values are 'in' and 'out'

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
    if 'start' in kwargs:
        start = kwargs.pop('start')
        if start is None:
            start = ticks_showing[-1]
            for line in ax.get_lines():
                if axis == 'x':
                    if start > min(line.get_xdata()):
                        start = min(line.get_xdata())
                if axis == 'y':
                    if start > min(line.get_ydata()):
                        start = min(line.get_ydata())
            # start = ticks_showing[0]
    else:
        start = min(ticks_showing)
        for line in ax.get_lines():
            if axis == 'x':
                if start > min(line.get_xdata()):
                    start = min(line.get_xdata())
            if axis == 'y':
                if start > min(line.get_ydata()):
                    start = min(line.get_ydata())
        # start = ticks_showing[0]

    if 'stop' in kwargs:
        stop = kwargs.pop('stop')
        if stop is None:
            stop = ticks_showing[0]
            for line in ax.get_lines():
                if axis == 'x':
                    if stop < max(line.get_xdata()):
                        stop = max(line.get_xdata())
                elif axis == 'y':
                    if stop < max(line.get_ydata()):
                        stop = max(line.get_ydata())
            # stop = ticks_showing[-1]
    else:
        stop = max(ticks_showing)
        for line in ax.get_lines():
            if axis == 'x':
                if stop < max(line.get_xdata()):
                    stop = max(line.get_xdata())
            elif axis == 'y':
                if stop < max(line.get_ydata()):
                    stop = max(line.get_ydata())
        # stop = ticks_showing[-1]

    if 'nticks' in kwargs:
        nticks = kwargs.pop('nticks')
        if nticks is None:
            if 'step' in kwargs:
                step = kwargs.pop('step')
                if step is None:
                    nticks = len(ticks_showing)
                else:
                    use_sep = True
            else:
                nticks = len(ticks_showing)
    elif 'step' in kwargs:
        step = kwargs.pop('step')
        if step is not None:
            use_sep = True
        else:
            nticks = len(ticks_showing)
    else:
        use_sep = True
        step = np.mean(np.diff(ticks_showing))

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
        ticks   = np.arange(start, stop + step*0.1, step)
        # print(step)
        # print(ticks)
    else:
        if nticks < 2:
            raise ValueError('nticks needs to be bigger or equal than 2.')
        ticks   = np.linspace(start, stop, nticks)
    # ticks shift to get better values (include zero)
    if any(x<0 for x in ticks) and any(x>0 for x in ticks) and 0 not in ticks:
        ticks = ticks-ticks[_index(ticks, 0)]
    elif stop-start > 5:
        ticks = ticks-(ticks[0]-int(ticks[0]))
    # ticks edges
    # print(ticks[-1])
    # print(ticks[-1]+np.mean(np.diff(ticks))*0.5)
    # print(stop)
    if ticks[-1]+np.mean(np.diff(ticks))*0.5 < stop:
        # print('ff')
        ticks = np.append(ticks, ticks[-1]+np.mean(np.diff(ticks)))
    if ticks[0]-np.mean(np.diff(ticks))*0.8 > start:
        ticks = np.append(ticks[0]-np.mean(np.diff(ticks)), ticks)

    # limits
    if len(ticks) < 2:
        raise ValueError(f'ticks = {ticks} has only one tick, please, reduce step.')
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
    #         string +=f'    nticks = {len(ticks)}\n'
    #
    #             # # ticks
    #             # ticks   = np.linspace(min_value, max_value, nticks)
    #
    #         string +='    print(nticks)\n'
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
    #             # if 'nticks' in kwargs:
    #             #     nticks = kwargs['nticks']
    #             #     if nticks is None:
    #             #         if 'step' in kwargs:
    #             #             step = kwargs['step']
    #             #             if step is None:
    #             #                 nticks = len(ticks_showing)
    #             #             else:
    #             #                 use_sep = True
    #             #         else:
    #             #             nticks = len(ticks_showing)
    #             # elif 'step' in kwargs:
    #             #     step = kwargs['step']
    #             #     if step is not None:
    #             #         use_sep = True
    #             #     else:
    #             #         nticks = len(ticks_showing)
    #             # else:
    #             #     use_sep = True
    #             #     step = np.mean(np.diff(ticks_showing))
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

def _set_xticks_old(ax=None, start=None, stop=None, pad=None, nticks=None, step=None, n_minor_ticks=None, fontproperties=None, **kwargs):
    """Set y ticks of a plot.

    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.
        start (float or int): start value for ticks (not the plot edge --- see ``pad``)
        stop (float or int): stop value for ticks (not the plot edge --- see ``pad``)
        nticks (int): Number of ticks. Ticks separation is calculated 
            accordingly. Overwrites step.
        step (float or int): Ticks separation.
        n_minor_ticks (int): Number of minor ticks between two major ticks.
        fontproperties: Label ticks font. Use ``matplotlib.font_manager.FontProperties``.
        pad (float or int): Distance between plot edge and the first tick in terms of tick separation. Typically, must be something between 0 and 1.
        # n_decimal_places (int): Number of decimal places to use for tick labels.
        direction (str, optional): default is 'out', possible values are 'in' and 'out'

    Note:
        To set minor and major ticks 'manually' use `xaxis.set_ticks() <https://matplotlib.org/3.2.2/api/_as_gen/matplotlib.axis.XAxis.set_ticks.html>`_, for example::

            _ = ax.xaxis.set_ticks([1, 2, 3, 4, 5, 6], minor=True)

        Set labels manually by using:

            _ = ax.set_yticklabels([1, 2, 3, 4, 5, 6])

    """
    if ax is None:
        ax = plt.gca()

    # get ticks showing
    # ticks_showing = get_xticks_showing(ax=ax)

    ## collecting kwargs
    # nticks   = None
    # step = None
    # n_minor_ticks  = None
    # fontproperties = None
    # pad            = None
    # n_decimal_places = None

    # # % tick params
    # tick_params = {'bottom': True,
    #                 'top': True,
    #                 # 'right': True,
    #                 # 'left': True,
    #                 'labelbottom': True,
    #                 # 'labelleft': True,
    #                 # 'labelright': False,
    #                 'labeltop': False,
    #                 'direction': 'out'}
    # for attr in tick_params:
    #     if attr in kwargs:
    #         tick_params[attr] = kwargs.pop(attr)
    #         if attr != 'direction':
    #             assert type(tick_params[attr]) == bool, f'{attr} must be True or False'
    # assert tick_params['direction'] in ('in', 'out',  'inout'), 'invalid value for direction.\nMust be either `in` or `out`'
    
    # % start and stop
    # if 'start' in kwargs:
    #     start = kwargs.pop('start')
    # # if start is None:
    # #     start = ticks_showing[-1]
    # #     for line in ax.get_lines():
    # #         if start > min(line.get_xdata()):
    # #             start = min(line.get_xdata())
    # if 'stop' in kwargs:
    #     stop = kwargs.pop('stop')
    # # if stop is None:
    # #     stop = ticks_showing[0]
    # #     for line in ax.get_lines():
    # #         if stop < max(line.get_xdata()):
    # #             stop = max(line.get_xdata())
    
    # # % minor ticks
    # if 'n_minor_ticks' in kwargs:
    #     n_minor_ticks = kwargs.pop('n_minor_ticks')

    # # % font  
    # if 'fontproperties' in kwargs:
    #     fontproperties = kwargs.pop('fontproperties')

    # # % padding
    # if 'pad' in kwargs:
    #     pad = kwargs.pop('pad')

    # # % n decimal places
    # if 'n_decimal_places' in kwargs:
    #     n_decimal_places = kwargs.pop('n_decimal_places')

    
    # # % nticks and tick_sep
    # if 'nticks' in kwargs:
    #     nticks = kwargs.pop('nticks')
    # if 'step' in kwargs:
    #     step = kwargs.pop('step')

    if nticks is not None:
        if nticks < 2:
            raise ValueError('nticks needs to be bigger or equal than 2.')
    if nticks is None and step is None:
        # nticks = len(ticks_showing)
        ticks_showing = get_xticks_showing(ax=ax)
        # print(ticks_showing)
        step     = np.mean(np.diff(ticks_showing))

    
    # # raise error for non recognized args
    # if len(kwargs.keys()) > 0:
    #     raise AttributeError(f'Args not recognized: {list(kwargs.keys())}')


    # % Ticks ======================
    if start is None and stop is None and pad is None:
        # ticks = ticks_showing
        # start = ticks[0]
        # stop  = ticks[-1]
        pass
    else:
        # initial
        if nticks is not None:
            ticks   = np.linspace(start, stop, nticks)
        else:
            ticks   = np.arange(start, stop + step*0.1, step) 

        assert len(ticks) > 1, 'number of ticks needs to be equal or larger than 2'

        # ticks shift to get better values (include zero)
        if any(x < 0 for x in ticks) and any(x > 0 for x in ticks) and 0 not in ticks:
            ticks = ticks-ticks[_index(ticks, 0)]
        elif stop-start > 5: # If zero is not included, round ticks if range is larger than 5
            ticks = ticks-(ticks[0]-int(ticks[0]))

        # ticks edges
        if ticks[-1]+np.mean(np.diff(ticks))*0.5 < stop:
            ticks = np.append(ticks, ticks[-1]+np.mean(np.diff(ticks)))
        if ticks[0]-np.mean(np.diff(ticks))*0.8 > start:
            ticks = np.append(ticks[0]-np.mean(np.diff(ticks)), ticks)

        # limits
        if pad is None:
            min_lim = start
            max_lim = stop
        else:
            try:
                if len(pad) == 2:
                    min_lim = ticks[0] - (ticks[1]-ticks[0])*pad[0]
                    max_lim = ticks[-1] + (ticks[1]-ticks[0])*pad[1]
                else:
                    min_lim = ticks[0] - (ticks[1]-ticks[0])*pad[0]
                    max_lim = ticks[-1] + (ticks[1]-ticks[0])*pad[0]
            except TypeError:
                min_lim = ticks[0] - (ticks[1]-ticks[0])*pad
                max_lim = ticks[-1] + (ticks[1]-ticks[0])*pad

        # # decimal places to show
        # if n_decimal_places is not None:
        #     ticks = [np.round(tick, n_decimal_places) for tick in ticks]


        # applying changes ======================
        # print(ticks)
    

        # set ticks
        if fontproperties is not None:
            dummy = ax.set_xticklabels([str(i) for i in ticks], fontproperties=fontproperties, visible=True)
        else:
            dummy = ax.set_xticks(ticks)

        # limits
        ax.set_xlim((min_lim, max_lim), auto=False)

    # minor ticks
    if n_minor_ticks is not None:
        ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))

    # direction
    for which in ('major', 'minor'):
        ax.tick_params(which=which, axis='x', **kwargs)
    
    return

def _set_yticks_old(ax=None, start=None, stop=None, pad=None, nticks=None, step=None, n_minor_ticks=None, fontproperties=None, **kwargs):
    """Set y ticks of a plot.

    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks. If None,
            last ax will be used.
        start (float or int): start value for ticks (not the plot edge --- see ``pad``)
        stop (float or int): stop value for ticks (not the plot edge --- see ``pad``)
        nticks (int): Number of ticks. Ticks separation is calculated 
            accordingly. Overwrites step.
        step (float or int): Ticks separation.
        n_minor_ticks (int): Number of minor ticks between two major ticks.
        fontproperties: Label ticks font. Use ``matplotlib.font_manager.FontProperties``.
        pad (float or int): Distance between plot edge and the first tick in terms of tick separation. Typically, must be something between 0 and 1.
        # n_decimal_places (int): Number of decimal places to use for tick labels.
        direction (str, optional): default is 'out', possible values are 'in' and 'out'

    Note:
        To set minor and major ticks 'manually' use `xaxis.set_ticks() <https://matplotlib.org/3.2.2/api/_as_gen/matplotlib.axis.XAxis.set_ticks.html>`_, for example::

            _ = ax.xaxis.set_ticks([1, 2, 3, 4, 5, 6], minor=True)

        Set labels manually by using:

            _ = ax.set_yticklabels([1, 2, 3, 4, 5, 6])

    """
    if ax is None:
        ax = plt.gca()

    # get ticks showing
    # ticks_showing = get_xticks_showing(ax=ax)

    ## collecting kwargs
    # nticks   = None
    # step = None
    # n_minor_ticks  = None
    # fontproperties = None
    # pad            = None
    # n_decimal_places = None

    # # % tick params
    # tick_params = {'bottom': True,
    #                 'top': True,
    #                 # 'right': True,
    #                 # 'left': True,
    #                 'labelbottom': True,
    #                 # 'labelleft': True,
    #                 # 'labelright': False,
    #                 'labeltop': False,
    #                 'direction': 'out'}
    # for attr in tick_params:
    #     if attr in kwargs:
    #         tick_params[attr] = kwargs.pop(attr)
    #         if attr != 'direction':
    #             assert type(tick_params[attr]) == bool, f'{attr} must be True or False'
    # assert tick_params['direction'] in ('in', 'out',  'inout'), 'invalid value for direction.\nMust be either `in` or `out`'
    
    # % start and stop
    # if 'start' in kwargs:
    #     start = kwargs.pop('start')
    # # if start is None:
    # #     start = ticks_showing[-1]
    # #     for line in ax.get_lines():
    # #         if start > min(line.get_xdata()):
    # #             start = min(line.get_xdata())
    # if 'stop' in kwargs:
    #     stop = kwargs.pop('stop')
    # # if stop is None:
    # #     stop = ticks_showing[0]
    # #     for line in ax.get_lines():
    # #         if stop < max(line.get_xdata()):
    # #             stop = max(line.get_xdata())
    
    # # % minor ticks
    # if 'n_minor_ticks' in kwargs:
    #     n_minor_ticks = kwargs.pop('n_minor_ticks')

    # # % font  
    # if 'fontproperties' in kwargs:
    #     fontproperties = kwargs.pop('fontproperties')

    # # % padding
    # if 'pad' in kwargs:
    #     pad = kwargs.pop('pad')

    # # % n decimal places
    # if 'n_decimal_places' in kwargs:
    #     n_decimal_places = kwargs.pop('n_decimal_places')

    
    # # % nticks and tick_sep
    # if 'nticks' in kwargs:
    #     nticks = kwargs.pop('nticks')
    # if 'step' in kwargs:
    #     step = kwargs.pop('step')

    if nticks is not None:
        if nticks < 2:
            raise ValueError('nticks needs to be bigger or equal than 2.')
    if nticks is None and step is None:
        # nticks = len(ticks_showing)
        ticks_showing = get_yticks_showing(ax=ax)
        # print(ticks_showing)
        step     = np.mean(np.diff(ticks_showing))

    
    # # raise error for non recognized args
    # if len(kwargs.keys()) > 0:
    #     raise AttributeError(f'Args not recognized: {list(kwargs.keys())}')


    # % Ticks ======================
    if start is None and stop is None and pad is None:
        # ticks = ticks_showing
        # start = ticks[0]
        # stop  = ticks[-1]
        pass
    else:
        # initial
        if nticks is not None:
            ticks   = np.linspace(start, stop, nticks)
        else:
            ticks   = np.arange(start, stop + step*0.1, step) 

        assert len(ticks) > 1, 'number of ticks needs to be equal or larger than 2'

        # ticks shift to get better values (include zero)
        if any(x < 0 for x in ticks) and any(x > 0 for x in ticks) and 0 not in ticks:
            ticks = ticks-ticks[_index(ticks, 0)]
        elif stop-start > 5: # If zero is not included, round ticks if range is larger than 5
            ticks = ticks-(ticks[0]-int(ticks[0]))

        # ticks edges
        if ticks[-1]+np.mean(np.diff(ticks))*0.5 < stop:
            ticks = np.append(ticks, ticks[-1]+np.mean(np.diff(ticks)))
        if ticks[0]-np.mean(np.diff(ticks))*0.8 > start:
            ticks = np.append(ticks[0]-np.mean(np.diff(ticks)), ticks)

        # limits
        if pad is None:
            min_lim = start
            max_lim = stop
        else:
            try:
                if len(pad) == 2:
                    min_lim = ticks[0] - (ticks[1]-ticks[0])*pad[0]
                    max_lim = ticks[-1] + (ticks[1]-ticks[0])*pad[1]
                else:
                    min_lim = ticks[0] - (ticks[1]-ticks[0])*pad[0]
                    max_lim = ticks[-1] + (ticks[1]-ticks[0])*pad[0]
            except TypeError:
                min_lim = ticks[0] - (ticks[1]-ticks[0])*pad
                max_lim = ticks[-1] + (ticks[1]-ticks[0])*pad

        # # decimal places to show
        # if n_decimal_places is not None:
        #     ticks = [np.round(tick, n_decimal_places) for tick in ticks]


        # applying changes ======================
        # print(ticks)
    

        # set ticks
        if fontproperties is not None:
            dummy = ax.set_yticklabels([str(i) for i in ticks], fontproperties=fontproperties, visible=True)
        else:
            dummy = ax.set_yticks(ticks)

        # limits
        ax.set_ylim((min_lim, max_lim), auto=False)

    # minor ticks
    if n_minor_ticks is not None:
        ax.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))

    # direction
    for which in ('major', 'minor'):
        ax.tick_params(which=which, axis='y', **kwargs)
    
    return
# %%