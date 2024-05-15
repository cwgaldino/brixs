#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Experimental functions and methods"""

# %% ------------------------- Standard Imports --------------------------- %% #
import matplotlib.pyplot as plt
import numpy as np
# %%

# %% -------------------------- brixs Imports ----------------------------- %% #
import brixs as br
# %%

# %% Spectra
def sequential_plot(self, xlim=None, ylim=None, keep=None, update_function=None, **kwargs):
    """[EXPERIMENTAL] plot where you can use up and down keys to flip through spectra.

    Warning: 
        Only one plot can be open at a time. Also, some default 
        matplotlib keybidings are changed (should not be a problem)

    Args:
        xlim (tuple, optional): min and max ``x`` value. Default is None (full
            data range).
        ylim (tuple, optional): min and max ``y`` value. Default is None.
        keep (int, optional): index of the spectrum to keep on screen.
        update_function (function, optional): function that is called when
            the left or right arrows are pressed. update_function must be a 
            function of type:

            >>> def update_function(ss, __i):
            >>>     plt.title(__i)
            >>>     ss[__i].plot(color='black', marker='o')
        
            where ``__i`` is updated in every iteraction.
        **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
    
    Returns:
        None
    """
    # vars
    self.__i = 0
    self.__xlim = xlim
    self.__ylim = ylim
    self.__ax   = None

    # change keybindings
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass
    
    # lims
    if self.__xlim is None:
        self.__xlim = [min(self[0].x), max(self[0].x)]
        for s in self:
            if min(s.x) < self.__xlim[0]:
                self.__xlim[0] = min(s.x)
            if max(s.x) > self.__xlim[1]:
                self.__xlim[1] = max(s.x)
    if self.__ylim is None:
        self.__ylim = [min(self[0].y), max(self[0].y)]
        for s in self:
            if min(s.y) < self.__ylim[0]:
                self.__ylim[0] = min(s.y)
            if max(s.y) > self.__ylim[1]:
                self.__ylim[1] = max(s.y)

    # keep
    if keep is not None:
        if isinstance(keep, str):
            assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
        else:
            assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

    # kwargs
    if 'color' not in kwargs:
        kwargs['color'] = 'black'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'

    # core update function
    if update_function is None:
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            plt.title(self.__i)
            if keep is not None:
                if keep == 'next':
                    if self.__i+1 < len(self):
                        ss[self.__i+1].plot(color='red', alpha=0.5)
                elif keep == 'previous':
                    if self.__i-1 >= 0:
                        ss[self.__i-1].plot(color='red', alpha=0.5)
                else:
                    ss[keep].plot(color='red', alpha=0.5)
            ss[self.__i].plot(**kwargs)
            
            if self.__xlim is not None:
                plt.xlim(self.__xlim)
            if self.__ylim is not None:
                plt.ylim(self.__ylim)
    else:
        # add counter and xlim/ylim to update function
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            update_function(ss, __i=self.__i)
            
            if self.__xlim is not None:
                plt.xlim(self.__xlim)
            if self.__ylim is not None:
                plt.ylim(self.__ylim)

    # keyboard events
    def keyboard(event, ss):
        # print(event.key)
        # print('keyboard')
        # print(event.key)
        if event.key == 'right' or event.key == 'up':
            self.__i = self.__i + 1

            # self.__ax.cla()
            # for line in self.__ax.lines:
            #     line.remove()
            self.__ax.lines.clear()
            _update(ss)
            plt.draw()

        elif event.key == 'left' or event.key == 'down':
            self.__i = self.__i - 1

            # self.__ax.cla()
            # for line in self.__ax.lines:
            #     line.remove()
            self.__ax.lines.clear()
            _update(ss)
            plt.draw()

    # mouse events
    def mouse(event):
        """Mouse can be used with a keyboard key"""
        # print('mouse')
        # print(event.key)
        # print(event.button)
        pass

    # axis zoom changes
    def on_xlims_change(event_ax):
        # print("updated xlims: ", event_ax.get_xlim())
        self.__xlim = event_ax.get_xlim()

    def on_ylims_change(event_ax):
        # print("updated ylims: ", event_ax.get_ylim())
        self.__ylim = event_ax.get_ylim()

    # plotting
    fig, self.__ax = plt.subplots(1, 1)
    _update(self)

    # register callbacks
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
    fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))

    self.__ax.callbacks.connect('xlim_changed', on_xlims_change)
    self.__ax.callbacks.connect('ylim_changed', on_ylims_change)
    return
br.Spectra.sequential_plot = sequential_plot

def sequential_plot2(self, xlim=None, ylim=None, keep=None, update_function=None, **kwargs):
    """[EXPERIMENTAL] plot where you can use arrows to flip through spectra.

    Warning:
        Old version, zoom is reset every time. Only one plot can be open at a time.
        Also, some default matplotlib keybidings are changed (should not be a problem)

    Args:
        xlim (tuple, optional): min and max ``x`` value. Default is None (full
            data range).
        ylim (tuple, optional): min and max ``y`` value. Default is None.
        keep (int, optional): index of the spectrum to keep on screen.
        update_function (function, optional): function that is called when
            the left or right arrows are pressed. update_function must be a 
            function of type:

            >>> def update_function(ss, __i):
            >>>     plt.title(__i)
            >>>     ss[__i].plot(color='black', marker='o')
        
            where ``__i`` is updated in every iteraction.
        **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
    
    Returns:
        None
    """
    self.__i = 0

    # change keybindings
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    # keep
    if keep is not None:
        if isinstance(keep, str):
            assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
        else:
            assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

    # kwargs
    if 'color' not in kwargs:
        kwargs['color'] = 'black'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'

    # core update function
    if update_function is None:
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            plt.title(self.__i)
            if keep is not None:
                if keep == 'next':
                    if self.__i+1 < len(self):
                        ss[self.__i+1].plot(color='red', alpha=0.5)
                elif keep == 'previous':
                    if self.__i-1 >= 0:
                        ss[self.__i-1].plot(color='red', alpha=0.5)
                else:
                    ss[keep].plot(color='red', alpha=0.5)
            ss[self.__i].plot(**kwargs)
            
            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)
    else:
        # add counter and xlim/ylim to update function
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            update_function(ss, __i=self.__i)
            
            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)

    # keyboard events
    def keyboard(event, ss):
        # print(event.key)
        # print('keyboard')
        # print(event.key)
        if event.key == 'right' or event.key == 'up':
            self.__i = self.__i + 1

            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'left' or event.key == 'down':
            self.__i = self.__i - 1

            plt.cla()
            _update(ss)
            plt.draw()

    # mouse events
    def mouse(event):
        """Mouse can be used with a keyboard key"""
        # print('mouse')
        # print(event.key)
        # print(event.button)
        pass

    # plotting
    fig = figmanip.figure()

    # register callbacks
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
    fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))

    _update(self)
    return
# br.Spectra.sequential_plot2 = sequential_plot2

def shift_plot(self, xlim=None, ylim=None, step=None, vlines=None, keep=None, update_function=None, **kwargs):
    """[EXPERIMENTAL] flip through spectra (up/down keys), shift spectrum (a/d keys)

    Warning: 
        Only one plot can be open at a time. Also, some default 
        matplotlib keybidings are changed (should not be a problem)

    Args:
        xlim (tuple, optional): min and max ``x`` value. Default is None (full
            data range).
        ylim (tuple, optional): min and max ``y`` value. Default is None.
        step (number, optional): step. If None, 1% of the min datapoint 
            separation is used.
        vlines (list or number, optional): vertical dashed lines for 
            reference, default is None.
        keep (int, optional): index of the spectrum to keep on screen.
        update_function (function, optional): function that is called when
            the left or right arrows are pressed. update_function must be a 
            function of type:

            >>> def update_function(ss, __i):
            >>>     plt.title(__i)
            >>>     ss[__i].plot(shift=ss._calculated_shift[__i])
        
            where ``__i`` is updated in every iteraction.
        **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
    
    Returns:
        None
    """
    self.__i = 0
    self.__xlim = xlim
    self.__ylim = ylim
    self.__ax   = None
    self._calculated_shift = np.array([0.0]*len(self))

    # set step
    if step is None: 
        if self.step is None:
            try:
                self.check_step()
                self.__step = self.step*0.5
            except:
                self.__step = min([np.mean(np.diff(s.x)) for s in self])*0.5
        else:
            self.__step = self.step*0.5
    else:
        self.__step = step

    # change keybindings
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    # keep
    if keep is not None:
        if isinstance(keep, str):
            assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
        else:
            assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

    # kwargs
    if 'color' not in kwargs:
        kwargs['color'] = 'black'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'

    # core update function
    if update_function is None:
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            plt.title(f'{self.__i}: {self._calculated_shift[self.__i]}')
            if keep is not None:
                if keep == 'next':
                    if self.__i+1 < len(self):
                        ss[self.__i+1].plot(shift=self._calculated_shift[self.__i+1], color='red', alpha=0.5)
                elif keep == 'previous':
                    if self.__i-1 >= 0:
                        ss[self.__i-1].plot(shift=self._calculated_shift[self.__i-1], color='red', alpha=0.5)
                else:
                    ss[keep].plot(shift=self._calculated_shift[keep], color='red', alpha=0.5)
            ss[self.__i].plot(shift=self._calculated_shift[self.__i], **kwargs)
            
            if vlines is not None:
                figmanip.vlines(vlines, color='red', ls='--', lw=1)

            if self.__xlim is not None:
                plt.xlim(self.__xlim)
            if self.__ylim is not None:
                plt.ylim(self.__ylim)
    else:
        # add counter and xlim/ylim to update function
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            update_function(ss, __i=self.__i)

            if vlines is not None:
                figmanip.vlines(vlines, color='red', ls='--', lw=1)
            
            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)

    # keyboard events
    def keyboard(event, ss):
        # print(event.key)
        # print('keyboard')
        # print(event.key)
        if event.key == 'up':
            self.__i = self.__i + 1

            self.__ax.callbacks.disconnect(2)
            self.__ax.callbacks.disconnect(3)
            self.__newx = False
            self.__newy = False

            plt.cla()
            _update(ss)
            plt.draw()

            self.__ax.callbacks.connect('xlim_changed', on_xlims_change)
            self.__ax.callbacks.connect('ylim_changed', on_ylims_change)
        elif event.key == 'down':
            self.__i = self.__i - 1
            
            self.__ax.callbacks.disconnect(2)
            self.__ax.callbacks.disconnect(3)
            self.__newx = False
            self.__newy = False

            plt.cla()
            _update(ss)
            plt.draw()

            self.__ax.callbacks.connect('xlim_changed', on_xlims_change)
            self.__ax.callbacks.connect('ylim_changed', on_ylims_change)
        elif event.key == 'right':
            self._calculated_shift[self.__i] += self.__step
            
            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'left':
            self._calculated_shift[self.__i] -= self.__step
            
            plt.cla()
            _update(ss)
            plt.draw()

    # mouse events
    def mouse(event):
        """Mouse can be used with a keyboard key"""
        # print('mouse')
        # print(event.key)
        # print(event.button)
        pass

    # axis zoom changes
    def on_xlims_change(event_ax):
        if self.__newx:
            # print("updated xlims: ", event_ax.get_xlim())
            self.__xlim = event_ax.get_xlim()
        self.__newx = True

    def on_ylims_change(event_ax):
        if self.__newy:
            # print("updated ylims: ", event_ax.get_ylim())
            self.__ylim = event_ax.get_ylim()
        self.__newy = True

    # plotting
    fig, self.__ax = plt.subplots(1, 1)

    # register callbacks
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
    fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))
    
    _update(self)
    return
br.Spectra.shift_plot = shift_plot

def shift_plot2(self, xlim=None, ylim=None, step=None, vlines=None, keep=None, update_function=None, **kwargs):
    """[EXPERIMENTAL] flip through spectra (up/down keys), shift spectrum (left/right keys)

    Warning:
        Old version, zoom is reset every time. Only one plot can be open at a time.
        Also, some default matplotlib keybidings are changed (should not be a problem)

    Args:
        xlim (tuple, optional): min and max ``x`` value. Default is None (full
            data range).
        ylim (tuple, optional): min and max ``y`` value. Default is None.
        step (number, optional): step. If None, 1% of the min datapoint 
            separation is used.
        vlines (list or number, optional): vertical dashed lines for 
            reference, default is None.
        keep (int, optional): index of the spectrum to keep on screen.
        update_function (function, optional): function that is called when
            the left or right arrows are pressed. update_function must be a 
            function of type:

            >>> def update_function(ss, __i):
            >>>     plt.title(__i)
            >>>     ss[__i].plot(color='black', marker='o')
        
            where ``__i`` is updated in every iteraction.
        **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
    
    Returns:
        None
    """
    # vars
    self.__i = 0
    self._calculated_shift = np.array([0.0]*len(self))

    # set step
    if step is None: 
        if self.step is None:
            try:
                self.check_step()
                self.__step = self.step*0.5
            except:
                self.__step = min([np.mean(np.diff(s.x)) for s in self])*0.5
        else:
            self.__step = self.step*0.5
    else:
        self.__step = step

    # change keybindings
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    # keep
    if keep is not None:
        if isinstance(keep, str):
            assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
        else:
            assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

    # kwargs
    if 'color' not in kwargs:
        kwargs['color'] = 'black'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'

    # core update function
    if update_function is None:
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            plt.title(f'{self.__i}: {self._calculated_shift[self.__i]}')
            if keep is not None:
                if keep == 'next':
                    if self.__i+1 < len(self):
                        ss[self.__i+1].plot(shift=self._calculated_shift[self.__i+1], color='red', alpha=0.5)
                elif keep == 'previous':
                    if self.__i-1 >= 0:
                        ss[self.__i-1].plot(shift=self._calculated_shift[self.__i-1], color='red', alpha=0.5)
                else:
                    ss[keep].plot(shift=self._calculated_shift[keep], color='red', alpha=0.5)
            ss[self.__i].plot(shift=self._calculated_shift[self.__i], **kwargs)
            
            if vlines is not None:
                figmanip.vlines(vlines, color='red', ls='--', lw=1)

            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)
    else:
        # add counter and xlim/ylim to update function
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            update_function(ss, __i=self.__i)

            if vlines is not None:
                figmanip.vlines(vlines, color='red', ls='--', lw=1)
            
            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)

    # keyboard events
    def keyboard(event, ss):
        # print(event.key)
        # print('keyboard')
        # print(event.key)
        if event.key == 'up':
            self.__i = self.__i + 1

            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'down':
            self.__i = self.__i - 1
            
            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'right':
            self._calculated_shift[self.__i] += self.__step
            
            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'left':
            self._calculated_shift[self.__i] -= self.__step
            
            plt.cla()
            _update(ss)
            plt.draw()

    # mouse events
    def mouse(event):
        """Mouse can be used with a keyboard key"""
        # print('mouse')
        # print(event.key)
        # print(event.button)
        pass

    # plotting
    fig = figmanip.figure()
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
    fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))
    _update(self)
    return
# br.Spectra.shift_plot2 = shift_plot2

def roll_plot(self, xlim=None, ylim=None, vlines=None, keep=None, update_function=None, **kwargs):
    """[EXPERIMENTAL] flip through spectra (up/down keys), roll spectrum (a/d keys).

    Warning: 
        Only one plot can be open at a time. Also, some default 
        matplotlib keybidings are changed (should not be a problem)

    Args:
        xlim (tuple, optional): min and max ``x`` value. Default is None (full
            data range).
        ylim (tuple, optional): min and max ``y`` value. Default is None.
        vlines (list or number, optional): vertical dashed lines for 
            reference, default is None.
        keep (int or str, optional): index of the spectrum to keep on screen. 
            Use 'previous' or 'next' to show previous or next spectrum.
        update_function (function, optional): function that is called when
            the left or right arrows are pressed. update_function must be a 
            function of type:

            >>> def update_function(ss, __i):
            >>>     plt.title(__i)
            >>>     ss[__i].plot(color='black', marker='o')
        
            where ``__i`` is updated in every iteraction. If update_function
            is given, `keep` is disregarded.
        **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
    
    Returns:
        None
    """
    # check step
    if self.check_step is None:
        try:
            self.check_step()
        except ValueError:
            raise ValueError('cannot apply roll because step is not the same. Use shift_plot().')

    # vars
    self.__i = 0
    self.__xlim = xlim
    self.__ylim = ylim
    self.__ax   = None
    self._calculated_roll = np.array([0]*len(self))

    # change keybindings
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    # lims
    if self.__xlim is None:
        self.__xlim = [min(self[0].x), max(self[0].x)]
        for s in self:
            if min(s.x) < self.__xlim[0]:
                self.__xlim[0] = min(s.x)
            if max(s.x) > self.__xlim[1]:
                self.__xlim[1] = max(s.x)
    if self.__ylim is None:
        self.__ylim = [min(self[0].y), max(self[0].y)]
        for s in self:
            if min(s.y) < self.__ylim[0]:
                self.__ylim[0] = min(s.y)
            if max(s.y) > self.__ylim[1]:
                self.__ylim[1] = max(s.y)

    # check if all data has well defined step
    for i, s in enumerate(self):
        try:
            s.check_step()
        except ValueError:
            raise ValueError(f'Spectrum number {i} has non-uniform x-coord.\nRoll can only be applied with uniform x-coord.\nMaybe used ss.interp() to fix that.')

    # keep
    if keep is not None:
        if isinstance(keep, str):
            assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
        else:
            assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

    # kwargs
    if 'color' not in kwargs:
        kwargs['color'] = 'black'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'

    # core update function
    if update_function is None:
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            plt.title(f'{self.__i}: {self._calculated_roll[self.__i]}')
            if keep is not None:
                if keep == 'next':
                    if self.__i+1 < len(self):
                        ss[self.__i+1].plot(shift=self._calculated_roll[self.__i+1]*ss[self.__i+1].step, color='red', alpha=0.5)
                elif keep == 'previous':
                    if self.__i-1 >= 0:
                        ss[self.__i-1].plot(shift=self._calculated_roll[self.__i-1]*ss[self.__i-1].step, color='red', alpha=0.5)
                else:
                    ss[keep].plot(shift=self._calculated_roll[keep]*ss[keep].step, color='red', alpha=0.5)
            ss[self.__i].plot(shift=self._calculated_roll[self.__i]*ss[self.__i].step, **kwargs)
            
            if vlines is not None:
                figmanip.vlines(vlines, color='red', ls='--', lw=1)

            if self.__xlim is not None:
                plt.xlim(self.__xlim)
            if self.__ylim is not None:
                plt.ylim(self.__ylim)
    else:
        # add counter and xlim/ylim to update function
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            update_function(ss, __i=self.__i)

            if vlines is not None:
                figmanip.vlines(vlines, color='red', ls='--', lw=1)
            
            if self.__xlim is not None:
                plt.xlim(self.__xlim)
            if self.__ylim is not None:
                plt.ylim(self.__ylim)

    # keyboard events
    def keyboard(event, ss):
        # print(event.key)
        # print('keyboard')
        # print(event.key)
        if event.key == 'up':
            self.__i = self.__i + 1

            # self.__ax.cla()
            self.__ax.lines.clear()
            _update(ss)
            plt.draw()
        elif event.key == 'down':
            self.__i = self.__i - 1

            # self.__ax.cla()
            self.__ax.lines.clear()
            _update(ss)
            plt.draw()
        elif event.key == 'right':
            self._calculated_roll[self.__i] += 1

            # self.__ax.cla()
            self.__ax.lines.clear()
            _update(ss)
            plt.draw()
        elif event.key == 'left':
            self._calculated_roll[self.__i] -= 1

            # self.__ax.cla()
            self.__ax.lines.clear()
            _update(ss)
            plt.draw()

    # axis zoom changes
    def on_xlims_change(event_ax):
        # print("updated xlims: ", event_ax.get_xlim())
        self.__xlim = event_ax.get_xlim()

    def on_ylims_change(event_ax):
        # print("updated ylims: ", event_ax.get_ylim())
        self.__ylim = event_ax.get_ylim()

    # plotting
    fig, self.__ax = plt.subplots(1, 1)
    _update(self)

    # register callbacks
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
    # fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))
    self.__ax.callbacks.connect('xlim_changed', on_xlims_change)
    self.__ax.callbacks.connect('ylim_changed', on_ylims_change)

    # _update(self)
    return
br.Spectra.roll_plot = roll_plot

def roll_plot2(self, xlim=None, ylim=None, vlines=None, keep=None, update_function=None, **kwargs):
    """[EXPERIMENTAL] flip through spectra (up/down keys), roll spectrum (left/right keys)

    Warning:
        Old version, zoom is reset every time. Only one plot can be open at a time.
        Also, some default matplotlib keybidings are changed (should not be a problem)

    Args:
        xlim (tuple, optional): min and max ``x`` value. Default is None (full
            data range).
        ylim (tuple, optional): min and max ``y`` value. Default is None.
        vlines (list or number, optional): vertical dashed lines for 
            reference, default is None.
        keep (int, optional): index of the spectrum to keep on screen.
        update_function (function, optional): function that is called when
            the left or right arrows are pressed. update_function must be a 
            function of type:

            >>> def update_function(ss, __i):
            >>>     plt.title(__i)
            >>>     ss[__i].plot(color='black', marker='o')
        
            where ``__i`` is updated in every iteraction.
        **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
    
    Returns:
        None
    """
    # check step
    if self.check_step is None:
        try:
            self.check_step()
        except ValueError:
            raise ValueError('cannot apply roll because step is not the same. Use shift_plot().')

    # vars
    self.__i = 0
    self._calculated_roll = np.array([0]*len(self))

    # change keybindings
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    # check if all data has well defined step
    for i, s in enumerate(self):
        try:
            s.check_step()
        except ValueError:
            raise ValueError(f'Spectrum number {i} has non-uniform x-coord.\nRoll can only be applied with uniform x-coord.\nMaybe used ss.interp() to fix that.')

    # keep
    if keep is not None:
        if isinstance(keep, str):
            assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
        else:
            assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

    # kwargs
    if 'color' not in kwargs:
        kwargs['color'] = 'black'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'

    # core update function
    if update_function is None:
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            plt.title(f'{self.__i}: {self._calculated_roll[self.__i]}')
            if keep is not None:
                if keep == 'next':
                    if self.__i+1 < len(self):
                        ss[self.__i+1].plot(shift=self._calculated_roll[self.__i+1]*ss[self.__i+1].step, color='red', alpha=0.5)
                elif keep == 'previous':
                    if self.__i-1 >= 0:
                        ss[self.__i-1].plot(shift=self._calculated_roll[self.__i-1]*ss[self.__i-1].step, color='red', alpha=0.5)
                else:
                    ss[keep].plot(shift=self._calculated_roll[keep]*ss[keep].step, color='red', alpha=0.5)
            ss[self.__i].plot(shift=self._calculated_roll[self.__i]*ss[self.__i].step, **kwargs)
            
            if vlines is not None:
                figmanip.vlines(vlines, color='red', ls='--', lw=1)

            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)
    else:
        # add counter and xlim/ylim to update function
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            update_function(ss, __i=self.__i)

            if vlines is not None:
                figmanip.vlines(vlines, color='red', ls='--', lw=1)
            
            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)

    # keyboard events
    def keyboard(event, ss):
        # print(event.key)
        # print('keyboard')
        # print(event.key)
        if event.key == 'up':
            self.__i = self.__i + 1

            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'down':
            self.__i = self.__i - 1
            
            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'right':
            self._calculated_roll[self.__i] += 1
            
            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'left':
            self._calculated_roll[self.__i] -= 1
            
            plt.cla()
            _update(ss)
            plt.draw()


    # plotting
    fig = figmanip.figure()
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))

    _update(self)
    return
# br.Spectra.roll_plot2 = roll_plot2

def roll_plot3(self, xlim=None, ylim=None, vlines=None, keep=None, update_function=None, **kwargs):
    """[experimental] flip through spectra (up/down keys), roll spectrum (left/right keys).

    Warning:
        [Experimental] Use mouse to click and drag spectra to the left or right.
        Only one plot can be open at a time. Also, some default matplotlib 
        keybidings are changed (should not be a problem)

    Args:
        xlim (tuple, optional): min and max ``x`` value. Default is None (full
            data range).
        ylim (tuple, optional): min and max ``y`` value. Default is None.
        vlines (list or number, optional): vertical dashed lines for 
            reference, default is None.
        keep (int, optional): index of the spectrum to keep on screen.
        update_function (function, optional): function that is called when
            the left or right arrows are pressed. update_function must be a 
            function of type:

            >>> def update_function(ss, __i):
            >>>     plt.title(__i)
            >>>     ss[__i].plot(color='black', marker='o')
        
            where ``__i`` is updated in every iteraction.
        **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
    
    Returns:
        None
    """
    # check step
    if self.step is None:
        try:
            self.check_step()
        except ValueError:
            raise ValueError('cannot apply roll because step is not the same. Use shift_plot().')

    # raise NotImplementedError('This is not implemented yet.')
    self.__i = 0
    self._calculated_roll = np.array([0]*len(self))
    self._mouse_press   = 0

    # change keybindings
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    # check if all data has well defined step
    for i, s in enumerate(self):
        try:
            s.check_step()
        except ValueError:
            raise ValueError(f'Spectrum number {i} has non-uniform x-coord.\nRoll can only be applied with uniform x-coord.\nMaybe used ss.interp() to fix that.')
    
    # keep
    if keep is not None:
        if isinstance(keep, str):
            assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
        else:
            assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

    # kwargs
    if 'color' not in kwargs:
        kwargs['color'] = 'black'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'

    # core update function
    if update_function is None:
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            plt.title(f'{self.__i}: {self._calculated_roll[self.__i]}')
            if keep is not None:
                if keep == 'next':
                    if self.__i+1 < len(self):
                        ss[self.__i+1].plot(shift=self._calculated_roll[self.__i+1]*ss[self.__i+1].step, color='red', alpha=0.5)
                elif keep == 'previous':
                    if self.__i-1 >= 0:
                        ss[self.__i-1].plot(shift=self._calculated_roll[self.__i-1]*ss[self.__i-1].step, color='red', alpha=0.5)
                else:
                    ss[keep].plot(shift=self._calculated_roll[keep]*ss[keep].step, color='red', alpha=0.5)
            ss[self.__i].plot(shift=self._calculated_roll[self.__i]*ss[self.__i].step, **kwargs)
            
            if vlines is not None:
                figmanip.vlines(vlines, color='red', ls='--', lw=1)

            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)
    else:
        # add counter and xlim/ylim to update function
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            update_function(ss, __i=self.__i)

            if vlines is not None:
                figmanip.vlines(vlines, color='red', ls='--', lw=1)
            
            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)

    # keyboard events
    def keyboard(event, ss):
        # print(event.key)
        # print('keyboard')
        # print(event.key)
        if event.key == 'up':
            self.__i = self.__i + 1

            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'down':
            self.__i = self.__i - 1
            
            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'right':
            self._calculated_roll[self.__i] += 1
            
            plt.cla()
            _update(ss)
            plt.draw()
        elif event.key == 'left':
            self._calculated_roll[self.__i] -= 1
            
            plt.cla()
            _update(ss)
            plt.draw()

    # mouse events
    def mouse_press(event):
        """Mouse can be used with a keyboard key"""
        if event.button is not None:
            if event.button == 1:
                self._mouse_press = event.xdata

    def mouse_release(event, ss):
        """Mouse can be used with a keyboard key"""
        if event.button is not None:
            if event.button == 1:
                self._calculated_roll[self.__i] += int(round((event.xdata - self._mouse_press)/self.step))   
                plt.cla()
                _update(ss)
                plt.draw()

    # plotting
    fig = figmanip.figure()
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
    fig.canvas.mpl_connect('button_press_event', lambda event: mouse_press(event))
    fig.canvas.mpl_connect('button_release_event', lambda event: mouse_release(event, ss=self))

    _update(self)
    return
# br.Spectra.roll_plot3 = roll_plot3
