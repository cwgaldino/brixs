#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Everyday use functions that eases matplotlib figure manipulation."""


# standard libraries
import sys
import numpy as np
from pathlib import Path
import copy
import warnings
from subprocess import Popen, PIPE
import decimal

# matplotlib libraries
from matplotlib.pyplot import get_current_fig_manager
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# backpack
from .arraymanip import index


def set_default_window_position(*args):
    """Set the default window position for when :py:func:`setWindowPosition` is called.

    Args:
        *args: A tuple like (x, y) or two separate x, y values (in px).
    """
    if len(args) > 1:
        x = int(args[0])
        y = int(args[1])
    elif len(args) == 1 and len(args[0]) == 2:
        x = int(args[0][0])
        y = int(args[0][1])
    else:
        warnings.warn('Wrong input')
        return

    global p
    p = (x, y)


def set_onclick(format='svg', resolution=300, folder=None):
    """Set the default format for saving figures using :py:func:`onclick()`.

    Args:
        format (string, optional): 'svg' or 'png'.
        resolution (int, optional): resolution for png images in dpi.
        folder (string or pathlib.Path): folderpath to save figures.
    """

    global onclick_fig_format, onclick_resolution, onclick_folder
    if format == 'png':
        onclick_fig_format = 'png'
        onclick_resolution = resolution
    else:
        onclick_fig_format = format

    if folder is None:
        onclick_folder = Path.cwd()
    else:
        onclick_folder = Path(folder)


def onclick(event):
    """This function is called everytime a mouse key is pressed whitin a figure.

    Right click:
        x value is copied to the clipboard (using xsel).
    Left click OR (y + Right click):
        y value is copied to the clipboard (using xsel).
    Middle click:
        Figure is saved as svg or png in the default folder and copied to the clipboard
        (see :py:func:`set_onclick()`).

    Note:
        The matplotlib figure must be started by :py:func:`figure`, and not
        by the default `matplotlib.pyplot.figure() <https://matplotlib.org/3.3.1/api/_as_gen/matplotlib.pyplot.figure.html>`_ fuction.

    Note:
        This function is based on xsel and xclip and therefore works only on linux
        machines (use ``sudo apt-get install -y xsel`` and
        ``sudo apt-get install -y xclip`` to install xsel and xclip, respectivelly).
    """

    # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
    #       ('double' if event.dblclick else 'single', event.button,
    #        event.x, event.y, event.xdata, event.ydata))

    # right
    if event.key == 'y' or event.button == 3:
        p = Popen(['xsel','-bi'], stdin=PIPE)  # ctrl+V
        p.communicate(input=bytes(f'{event.ydata}'.encode()))
    else: # left
        p = Popen(['xsel','-bi'], stdin=PIPE)  # ctrl+V
        p.communicate(input=bytes(f'{event.xdata}'.encode()))

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

        if onclick_fig_format == 'svg':
            plt.savefig('.temporary_fig.svg')
            p = Popen([f'xclip -selection clipboard -t image/svg+xml -i {onclick_folder/".temporary_fig.svg"}'], shell=True)  # ctrl+V
            p = Popen([f'echo {onclick_folder/".temporary_fig.svg"}  |xclip -in -selection primary -target text/uri-list'], shell=True)  # ctrl+V

        elif onclick_fig_format == 'png':
            plt.savefig('.temporary_fig.png', dpi=onclick_resolution)
            p = Popen([f'xclip -selection clipboard -t image/png -i {onclick_folder/".temporary_fig.png"}'], shell=True)  # ctrl+V




def figure(**kwargs):
    """Create figure object attached to :py:func:`onclick`.

    Args:
        **kwargs: kwargs are passed to ``plt.figure()``.
    """
    fig = plt.figure(**kwargs)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    return fig


def setWindowPosition(*args):
    """Change position of a maptplotlib figure on screen.

    Tipically, (0, 0) is the top left corner.

    Args:
        *args: A tuple like (x, y) or two separate x, y values (in px).
    """
    if len(args) > 1:
        x = int(args[0])
        y = int(args[1])
    elif len(args) == 1 and len(args[0]) == 2:
        x = int(args[0][0])
        y = int(args[0][1])
    elif len(args) == 0:
        try:
            global p
            setWindowPosition(p)
            return
        except:
            pass
    else:
        warnings.warn('Wrong input')
        return

    figManager = get_current_fig_manager()
    width, height = getWindowSize()

    try:  # tested on tKinter backend
        figureGeometry = str(width) + 'x' + str(height) + '+' + str(x) + '+' + str(y)
        figManager.window.wm_geometry(figureGeometry)

    except AttributeError:
        try:  # tested on qt4 and qt5 backends
            figManager.window.setGeometry(int(x), int(y), width, height)
        except AttributeError:
            warnings.warn('Backend not suported.')


def setWindowSize(*args):
    """Change the size of the window of a matplotlib figure.

    Args:
        *args: A tuple like (width, height) or two separate width, height  values (in px).
    """
    if len(args) > 1:
        width = int(args[0])
        height = int(args[1])
    elif len(args) == 1 and len(args[0]) == 2:
        width = int(args[0][0])
        height = int(args[0][1])
    else:
        warnings.warn('Wrong input')
        return

    figManager = get_current_fig_manager()
    x,y = getWindowPosition()

    try:  # tested on tKinter backend
        figureGeometry = str(width) + 'x' + str(height) + '+' + str(x) + '+' + str(y)
        figManager.window.wm_geometry(figureGeometry)

    except AttributeError:
        try:  # tested on qt4 and qt5 backends
            figManager.window.setGeometry(x, y, width, height)
        except AttributeError:
            warnings.warn('Backend not suported.')


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


def getWindowPosition():
    """Return the position of a matplotlib position on the screen.

    Tipically, (0, 0) is the top left corner of your monitor.

    Returns:
        Tuple with the x and y position.
    """
    figManager = get_current_fig_manager()

    try:  # tested under tKinter backend
        return (figManager.window.winfo_x(), figManager.window.winfo_y())

    except AttributeError:  # tested under qt4 and qt5 backends
        try:
            return (figManager.window.geometry().x(), figManager.window.geometry().y())
        except AttributeError:
            warnings.warn('Backend not suported.')
            return (0, 0)


def getWindowSize():
    """Returns the size of the window of a matplotlib figure.

    Returns:
        Tuple with the width and height values.
    """
    figManager = get_current_fig_manager()

    try:  # tested on tKinter backend
        return (figManager.window.winfo_width(), figManager.window.winfo_height())

    except AttributeError:  # tested on qt4 and qt5 backends
        try:
            return (figManager.window.geometry().width(), figManager.window.geometry().height())
        except AttributeError:
            warnings.warn('Backend not suported.')
            return (0, 0)


def zoom(start, stop, fig=None, marginx=5, marginy=5):
    """Zoom up portion of a plot from start to stop.

    Args:
        start (float or int): initial x value.
        stop (float or int): final x value.
        fig (int, optional): number of the figure. If None, current figure is used.
        marginx (int, optional): margin value between data and the edges of plot in percentage.
        marginy (int, optional): margin value between data and the edges of plot in percentage.
    """
    if fig is None:
        fig = plt.gcf()

    ymax = 0
    ymin = 0

    for axis in fig.axes:
        for line in axis.get_lines():
            try:
                ymax_temp = max(line.get_data()[1][index(line.get_data()[0], start): index(line.get_data()[0], stop)])
                ymin_temp = min(line.get_data()[1][index(line.get_data()[0], start): index(line.get_data()[0], stop)])
            except ValueError:
                warnings.warn("All points of some data are outside of the required range.")
            try:
                if ymax_temp > ymax:
                    ymax = copy.copy(ymax_temp)
                if ymin_temp < ymin:
                    ymin = copy.copy(ymin_temp)

                m =  (ymax-ymin)*marginy/100
                plt.ylim(ymin-m, ymax+m)

                m =  (stop-start)*marginx/100
                plt.xlim(start, stop)
            except UnboundLocalError:
                warnings.warn("All data are outside of the required range. Cannot zoom.")


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

    if figs is 'all':
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
        *Args: A tuple with values to convert

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
        *Args: A tuple with values to convert

    Returns:
        A tuple with values converted
    """
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def n_decimal_places(number):
    """Return the number of decimal places of a number.

    Args:
        number (float or int): number.
    """

    return abs(decimal.Decimal(str(number)).as_tuple().exponent)


def n_digits(number):
    """Return the number of digits of a number.

    Args:
        number (float or int): number.
    """
    if n_decimal_places(number) != 0:
        return (len(str(int(np.around(number))) ), n_decimal_places(number))
    else:
        return (len(str(int(np.around(number))) ), 0)


def set_ticks(ax, axis='x', **kwargs):
    """Set ticks of a plot.

    Args:
        ax (matplotlib.axes.Axes): The axes of the subplot to set ticks.
        axis (string, optional): possible values are 'x' or 'y'.
        **kwargs:

            =======================================================  =================================
            min_value (float or int)                                 start value for ticks and plot
            max_value (float or int)                                 stop value for ticks and plot
            n_ticks (int)                                            Number of ticks. Ticks separation
                                                                     is calculated accordingly and
                                                                     this parameter overwrites
                                                                     ticks_sep.
            ticks_sep (float or int)                                 Ticks separation.
            n_minor_ticks (int)                                      Number of minor ticks between
                                                                     two major ticks.
            fontproperties (matplotlib.font_manager.FontProperties)  Label ticks font.
            pad (float or int)                                       Distance between plot edge and
                                                                     the first tick
                                                                     in terms of tick separation.
                                                                     Tipically, must be something
                                                                     between 0 and 1.
            n_decimal_places (int)                                   Number of decimal places to use
                                                                     at tick labels.
            =======================================================  =================================

    Note:
        To set minor and major ticks 'manually' use `XAxis.set_ticks() <https://matplotlib.org/3.2.2/api/_as_gen/matplotlib.axis.XAxis.set_ticks.html>`_, for example::

            ax.xaxis.set_ticks([1, 2, 3, 4, 5, 6], minor=True)
    """
    use_sep = False

    if axis == 'y':
        ticks_showing = [y for y in ax.get_yticks() if y >= ax.get_ylim()[0] and y <= ax.get_ylim()[1]]
    elif axis == 'x':
        ticks_showing = [x for x in ax.get_xticks() if x >= ax.get_xlim()[0] and x <= ax.get_xlim()[1]]
    else:
        raise ValueError("Axis must be either 'x' or 'y'.")

    ## collecting kwargs
    if 'min_value' in kwargs:
        min_value = kwargs['min_value']
        if min_value is None:
            min_value = ticks_showing[0]
    else:
        min_value = ticks_showing[0]

    if 'max_value' in kwargs:
        max_value = kwargs['max_value']
        if max_value is None:
            max_value = ticks_showing[-1]
    else:
        max_value = ticks_showing[-1]

    if 'n_ticks' in kwargs:
        n_ticks = kwargs['n_ticks']
        if n_ticks is None:
            if 'ticks_sep' in kwargs:
                ticks_sep = kwargs['ticks_sep']
                if ticks_sep is None:
                    n_ticks = len(ticks_showing)
                else:
                    use_sep = True
            else:
                n_ticks = len(ticks_showing)
    else:
        if 'ticks_sep' in kwargs:
            ticks_sep = kwargs['ticks_sep']
            if ticks_sep is not None:
                use_sep = True
            else:
                n_ticks = len(ticks_showing)
        else:
            n_ticks = len(ticks_showing)

    if 'n_minor_ticks' in kwargs:
        n_minor_ticks = kwargs['n_minor_ticks']
        if n_minor_ticks is None:
            n_minor_ticks = 2
    else:
        n_minor_ticks = 2

    if 'fontproperties' in kwargs:
        fontproperties = kwargs['fontproperties']
    else:
        fontproperties = None

    if 'pad' in kwargs:
        pad = kwargs['pad']
    else:
        pad = None

    # ticks
    if use_sep:
        ticks   = np.arange(min_value, max_value+ticks_sep*0.1, ticks_sep)

    else:
        ticks   = np.linspace(min_value, max_value, n_ticks)

    # limits
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
        n_decimal_places2 = kwargs['n_decimal_places']
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

    # applying changes ======================
    s = '{' + f'0:.{n_decimal_places2}f' + '}'
    if axis == 'y':
        if fontproperties is None:
            pass
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
            pass
        else:
            dummy = ax.set_xticks(ticks)
            dummy = ax.set_xticklabels([s.format(i) for i in ticks], fontproperties=fontproperties, visible=True)

        # minor ticks
        ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))

        # limits
        if pad is not None:
            ax.set_xlim((min_lim, max_lim), auto=False)


def remove_ticks_edge(ax):
    """Remove ticks the fall over the edges of the plot.

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


def axBox2figBox(ax, points):
    """Transform 'bbox like' axes position values to percentage fig position.

    Args:
        ax (matplotlib.axes.Axes): axes instance.
        points (list): [x_init, y_init, x_final, y_final]

    Retuns:
        'bbox like' figure positions.
    """
    [x_init, y_init, x_final, y_final] = points

    x_init  = axPos2figPos(ax, x_init, direction='x')
    y_init  = axPos2figPos(ax, y_init, direction='y')
    delta_x = axPos2figPos(ax, x_final, direction='x') - x_init
    delta_y = axPos2figPos(ax, y_final, direction='y') - y_init

    return [x_init, y_init, delta_x, delta_y]


def axPos2figPos(ax, value, direction='x'):
    """Transform axes position values to figure percentage position.

    Args:
        ax (matplotlib.axes.Axes): axes instance.
        value (float): value
        direction (string, optional): 'x' or 'y'.

    Retuns:
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

# OLD ====================================================

# def getFigureSize(fig=None):
#     if fig is None:
#         fig = plt.gcf()
#     return [x for x in fig.bbox_inches.get_points()[1]]

# def plot_pause(data_list, ax=None, pause=0.5, replace=False, **kwargs):
#     """Call plt.plot() function and wait for keyboard input.
#
#     Args:
#         data_list (list of arrays):
#         ax
#         plot_pause
#         replace
#         **kwargs: kwargs are passed to ``plt.plot()``.
#     """
#
#     if ax is None:
#         fig = figure()
#         ax = fig.add_subplot(111)
#
#     for data in data_list:
#         if replace:
#             try:
#                 print(curve)
#                 curve[-1].set_visible(False)
#             except Exception as e: print(e)
#         curve = ax.plot(data[:, 0], data[:, 1], **kwargs)
#         plt.pause(pause)
#         key = input("Press Enter to continue or q to quit...")
#         if key == 'q':
#             break


# import matplotlib.patches as patches
# rect = patches.Rectangle((16, 0), 10, 10, linewidth=4, edgecolor='r', facecolor='white', zorder=2)
# ax.add_patch(rect)


# def plot_ruler(orientation, height=4, width=4):
#     """Open matplotlib figure which can be used as a ruler.
#
#     Ruler goes from 0 to 1.
#
#     Args:
#         orientation (string): 'v' for vertical ruler and 'h' for horizontal ruler.
#         height (float): height in cm. If orientation is 'v', width is ignored.
#         width (float): width in cm. If orientation is 'h', height is ignored.
#     """
#     aux = plt.figure()
#     aux_ax = plt.subplot(111)
#     temp = np.linspace(0, 1, 11)
#     if orientation == 'v':
#         aux.set_size_inches(1, cm2inch(height)[0])
#         aux.subplots_adjust(left=0, bottom=0, right=.001, top=1)
#         aux_ax.set_yticks(temp, minor=False)  # Major
#         aux_ax.yaxis.set_minor_locator(MultipleLocator(temp[1]/5))  # Minor
#     else:
#         aux.set_size_inches(cm2inch(width)[0], .5)
#         aux.subplots_adjust(left=0, bottom=0, right=1, top=.01)
#         aux_ax.set_xticks(temp, minor=False)  # Major
#         aux_ax.xaxis.set_minor_locator(MultipleLocator(temp[1]/5))  # Minor
#     aux_ax.tick_params(which='major',direction='out',top=True,left=True,right=True, labeltop=True, labelright=True)
#     aux_ax.tick_params(which='minor',direction='out',top=True,left=True,right=True)
