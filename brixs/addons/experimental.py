#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Experimental functions and methods"""

# %% ------------------------- Standard Imports --------------------------- %% #
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
# %%

# %% -------------------------- brixs Imports ----------------------------- %% #
import brixs as br
# %%


# %% --------------------------------------------------------------------- %% #
try:
    import h5py
except:
    pass
def create_hdf5(self, filename, verbose=False):
    """Create an HDF5 file in NeXus NXdata format containing stacked spectra.
	
	Check if the x range for each spectrum is correct before creating HDF file. 
	
	Parameters
    ----------
    filename : str
        Name of the output file.
    spectra : list
        List of spectra objects. Each spectrum must have attributes:
        - Energy : incident energy value
        - x      : energy-loss axis
        - y      : signal values
    """
    spectra  = spectra.fix_monotonicity().interp()
    # --- Checking spectra consistency ---
    E_inc_all = np.array([float(_s.Energy) for _s in spectra], dtype=np.float64)   # X axis
    E_loss_ref = np.asarray(spectra[0].x, dtype=np.float64)                      # Y axis (reference)

    # Verify that all energy-loss axes have the same shape and values
    for k, _s in enumerate(spectra):
        xk = np.asarray(_s.x, dtype=np.float64)
        if xk.shape != E_loss_ref.shape or not np.allclose(xk, E_loss_ref, rtol=0, atol=1e-12):
            raise ValueError(f"Spectrum {k} has a different energy-loss axis. "
                             "Standardize or interpolate before saving as NXdata.")

    # Sort spectra by incident energy
    order = np.argsort(E_inc_all)
    E_inc = E_inc_all[order]

    # Stack signals according to sorted incident energies
    Z = np.stack([np.asarray(spectra[i].y, dtype=np.float32) for i in order], axis=0)

    # Sanity check for expected dimensions
    if Z.shape != (E_inc.size, E_loss_ref.size):
        raise RuntimeError("Unexpected shape for Z. "
                           f"Expected ({E_inc.size}, {E_loss_ref.size}), got {Z.shape}.")
    
    if not filename.lower().endswith('.h5'): # If it does not end with '.h5', the extension is added.
        filename = f"{filename}.h5"

    #print("Shapes -> E_inc:", E_inc.shape, " E_loss:", E_loss_ref.shape, " Z:", Z.shape)

    # --- Create NeXus data (NXdata) ---
    with h5py.File(filename, 'w') as f:
        f.attrs['default'] = 'entry'

        nxentry = f.create_group('entry')
        nxentry.attrs['NX_class'] = 'NXentry'
        nxentry.attrs['default'] = 'data'

        nxdata = nxentry.create_group('data')
        nxdata.attrs['NX_class'] = 'NXdata'

        # Axes datasets
        dEinc  = nxdata.create_dataset('Energy Incident (eV)', data=E_inc)        # X axis
        dEloss = nxdata.create_dataset('Energy Loss (eV)',    data=E_loss_ref)    # Y axis

        # Signal dataset
        ds = nxdata.create_dataset('spectra', data=Z, chunks=True, compression='gzip', compression_opts=4)

        # --- Essential NXdata attributes ---
        nxdata.attrs['signal'] = 'spectra'
        nxdata.attrs['axes'] = ['Energy Incident (eV)', 'Energy Loss (eV)']

        # Axis-to-dimension mapping
        ds.attrs['Energy_Incident_indices'] = np.array([0])
        ds.attrs['Energy_Loss_indices']     = np.array([1])

        # Viewer hints
        ds.attrs['interpretation'] = 'image'       # 2D image/matrix
        ds.attrs['long_name'] = 'Counts (Arb. Units)'           # Label for the signal

        dEinc.attrs['units'] = 'eV'
        dEinc.attrs['long_name'] = 'Incident energy'

        dEloss.attrs['units'] = 'eV'
        dEloss.attrs['long_name'] = 'Energy loss'

    if verbose:
        print(f"HDF5 file successfully created: {filename}")
br.Spectra.create_hdf5 = create_hdf5
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
    fig = br.figure()

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
                br.axvlines(vlines, colors='red', ls='--', lw=1)

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
                br.axvlines(vlines, colors='red', ls='--', lw=1)
            
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
                br.axvlines(vlines, colors='red', ls='--', lw=1)

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
                br.axvlines(vlines, colors='red', ls='--', lw=1)
            
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
    fig = br.figure()
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
                br.axlines(vlines, colors='red', ls='--', lw=1)

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
                br.axvlines(vlines, colors='red', ls='--', lw=1)
            
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
                br.axvlines(vlines, colors='red', ls='--', lw=1)

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
                br.axvlines(vlines, colors='red', ls='--', lw=1)
            
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
    fig = br.figure()
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
                br.axvlines(vlines, colors='red', ls='--', lw=1)

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
                br.axvlines(vlines, colors='red', ls='--', lw=1)
            
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
    fig = br.figure()
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
    fig.canvas.mpl_connect('button_press_event', lambda event: mouse_press(event))
    fig.canvas.mpl_connect('button_release_event', lambda event: mouse_release(event, ss=self))

    _update(self)
    return
# br.Spectra.roll_plot3 = roll_plot3

# %% IMAGE
def linecuts(self, axis=0, xlim=None, ylim=None, ilim=None, keep=None, **kwargs):
    """[EXPERIMENTAL] Plot image and flip trhough linecuts with keyboard arrows.

    Args:
        axis (int or string, optional): Axis along linecuts.
            By default, linecuts are in the vertical (0) direction.
        xlim, ylim, ilim (tuple, optional): Image ploting limits.
        keep (int or str, optional): index of the spectrum to keep on screen. 
            Use 'previous' or 'next' to show previous or next spectrum.
        **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

    Returns:
        None
    """
    # axis
    if axis == 0:
        ss = self.columns
    elif axis == 1:
        raise NotImplementedError('Not implemented yet for axis = 1')
        # ss = self.rows
    else:
        raise ValueError('axis must be 0 or 1')

    # change keybindings
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass
    
    # vars
    self.__i    = 0
    self.__axes = None
    self.__xlim = xlim
    self.__ylim = ylim
    self.__ilim = ilim

    # lims
    if self.__xlim is None:
        self.__xlim = [min(self._x_centers), max(self._x_centers)]
    if self.__ylim is None:
        self.__ylim = [min(self._y_centers), max(self._y_centers)]
    if self.__ilim is None:
        self.__ilim     = [min([min(l) for l in self.data]), max([max(l) for l in self.data])]

    # keep
    if keep is not None:
        if isinstance(keep, str):
            assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
        else:
            assert abs(keep) < len(ss), 'keep must be a valid spectrum index, or "previous/next"'

    # kwargs
    if 'color' not in kwargs:
        kwargs['color'] = 'black'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'

    # update function
    def update(ss):
        # plot spectrum
        ax = self.__axes[0]
        if axis == 0:
            ax.set_title(f'{self.__i}: {self.x_centers[self.__i]} eV')
        if keep is not None:
            if keep == 'next':
                if self.__i+1 < len(self):
                    ss[self.__i+1].plot(color='red', alpha=0.5)
            elif keep == 'previous':
                if self.__i-1 >= 0:
                    ss[self.__i-1].plot(color='red', alpha=0.5)
            else:
                ss[keep].plot(color='red', alpha=0.5)
        ss[self.__i].plot(ax=ax, **kwargs)

        # lim
        if self.__ylim is not None:
            ax.set_xlim(self.__ylim)
        if self.__ilim is not None:
            ax.set_ylim(self.__ilim)



        # plot map (shouldn't we put this outside of the update function?)
        ax = self.__axes[1]
        
        # lim
        if self.__ilim is not None:
            self.pcolormesh(ax=ax, vmin=self.__ilim[0], vmax=self.__ilim[1])
        else:
            self.pcolormesh(ax=ax)
        if self.__xlim is not None:
            if axis == 0:
                ax.set_xlim(self.__xlim)
            else:
                pass
                # ax.set_xlim(self.__ylim)
        if self.__ylim is not None:
            if axis == 0:
                ax.set_ylim(self.__ylim)
            else:
                pass
                # ax.set_ylim(self.__ylim)

        # vline
        if axis == 0:
            if keep is not None:
                if keep == 'next':
                    if self.__i+1 < len(self):
                        ax.plot([self.x_centers[self.__i+1], self.x_centers[self.__i+1]], [self.__ylim[0], self.__ylim[1]], color='red', alpha=0.5)
                elif keep == 'previous':
                    if self.__i-1 >= 0:
                        ax.plot([self.x_centers[self.__i-1], self.x_centers[self.__i-1]], [self.__ylim[0], self.__ylim[1]], color='red', alpha=0.5)
                else:
                    ax.plot([self.x_centers[keep], self.x_centers[keep]], [self.__ylim[0], self.__ylim[1]], color='red', alpha=0.5)
            ax.plot([self.x_centers[self.__i], self.x_centers[self.__i]], [self.__ylim[0], self.__ylim[1]], color='white')
            # E = self.x_centers[self.__i]
            # br.axvlines(E, ax=ax, color='red', lw=1)
        else:
            E = self.y_centers[self.__i]
            br.axhlines(E, ax=ax, colors='red', lw=1)

        

    def _update(ss):
        if self.__i >= len(ss):
            self.__i = len(ss) - 1
        elif self.__i < 0:
            self.__i = 0

        update(ss)

    # keyboard events
    def keyboard(event, ss):
        # print(event.key)
        # print('keyboard')
        # print(event.key)
        if event.key == 'right' or event.key == 'up':
            self.__i = self.__i + 1

            for ax in self.__axes:
                # ax.cla()
                ax.lines.clear()
            _update(ss)
            plt.draw()
        elif event.key == 'left' or event.key == 'down':
            self.__i = self.__i - 1
            
            for ax in self.__axes:
                # ax.cla()
                ax.lines.clear()
            _update(ss)
            plt.draw()

    # axis zoom changes
    def on_1_xlims_change(event_ax):
        # print("updated xlims: ", event_ax.get_xlim())
        self.__xlim = event_ax.get_xlim()

    def on_1_ylims_change(event_ax):
        # print("updated ylims: ", event_ax.get_ylim())
        self.__ylim = event_ax.get_ylim()

    def on_0_xlims_change(event_ax):
        # print("updated xlims: ", event_ax.get_xlim())
        self.__ylim = event_ax.get_xlim()

    def on_0_ylims_change(event_ax):
        # print("updated ylims: ", event_ax.get_ylim())
        self.__ilim = event_ax.get_ylim()

    # figure
    fig, self.__axes = br.subplots(nrows=2, ncols=1)
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=ss))
    _update(ss)

    self.__axes[0].callbacks.connect('xlim_changed', on_0_xlims_change)
    self.__axes[0].callbacks.connect('ylim_changed', on_0_ylims_change)
    self.__axes[1].callbacks.connect('xlim_changed', on_1_xlims_change)
    self.__axes[1].callbacks.connect('ylim_changed', on_1_ylims_change)
# br.Image.linecuts = linecuts

def linecuts2(self, axis=0, xlim=None, ylim=None):
    """[EXPERIMENTAL] Plot image and flip trhough linecuts with keyboard arrows.

    Warning: old version, with auto rescaling.

    Args:
        axis (int or string, optional): Axis along linecuts.
            By default, linecuts are in the vertical (0) direction.
        xlim, ylim (tuple, optional): spectra ploting limits (not image 
            ploting limits).

    Returns:
        None
    """
    # axis
    if axis == 0:
        ss = self.columns
    elif axis == 1:
        ss = self.rows
    else:
        raise ValueError('axis must be 0 or 1')
    
    # update function
    def update(ss):
        # plot
        ax = axes[0]
        ss[ss.__i].plot(ax=ax, color='black', marker='o')

        # lim
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        # br.label_rixs(ax=ax)

        # plot map (shouldn't we put this outside of the update function?)
        ax = axes[1]
        
        if ylim is not None:
            self.pcolormesh(ax=ax, vmin=ylim[0], vmax=ylim[1])
        else:
            self.pcolormesh(ax=ax)
        if xlim is not None:
            if axis == 0:
                ax.set_ylim(xlim)
            else:
                ax.set_xlim(xlim)

        # vline
        if axis == 0:
            E = self.x_centers[ss.__i]
            br.axvlines(E, ax=ax, colors='red', lw=1)
        else:
            E = self.y_centers[ss.__i]
            br.axhlines(E, ax=ax, colors='red', lw=1)

        # title
        axes[0].set_title(f'{ss.__i}: {E} eV')

    def _update(ss):
        if ss.__i >= len(ss):
            ss.__i = len(ss) - 1
        elif ss.__i < 0:
            ss.__i = 0

        update(ss)

    # keyboard events
    def keyboard(event, ss):
        # print(event.key)
        # print('keyboard')
        # print(event.key)
        if event.key == 'right' or event.key == 'up':
            ss.__i = ss.__i + 1

            for ax in axes:
                ax.cla()
            _update(ss)
            # for ax in axes:
            plt.draw()
        elif event.key == 'left' or event.key == 'down':
            ss.__i = ss.__i - 1
            
            for ax in axes:
                ax.cla()
            _update(ss)
            # for ax in axes:
            plt.draw()

    # figure
    fig, axes = br.subplots(nrows=2, ncols=1)
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=ss))
    _update(ss)
# br.Image.linecuts2 = linecuts2

def roll_plot(self, axis=0, vlines=None, hlines=None, **kwargs):
    """[EXPERIMENTAL] Display data as an image. Wrapper for `matplotlib.pyplot.imshow()`_.

    Warning:
        Pixels are always square. For irregular pixel row/columns, see Image.pcolormesh()

    Args:
        axis (int or string, optional): Axis along roll.
            By default, columns are rolled down or up (0).
        vlines, hlines (list or number, optional): vertical or horizontal
            dashed lines for reference, default is None.
        **kwargs: kwargs are passed to `matplotlib.pyplot.imshow()`_.

    If not specified, the following parameters are passed to `matplotlib.pyplot.imshow()`_:

    Args:
        cmap: The Colormap instance. Default is 'jet'.
        aspect: The aspect ratio of the Axes. Default is 'auto'. If 'equal',
            an aspect ratio of 1 will be used (pixels will be square).
        origin: Location of the [0, 0] index. default is 'lower'.
        interpolation: The interpolation method used. Default is 'none'.
            Supported values are 'none', 'antialiased', 'nearest', 'bilinear',
            'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite',
            'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell',
            'sinc', 'lanczos', 'blackman'.
        extent: minimun and maximum x and y values. Default will be given by
            the Image.x and Image.y attributes.
        vmin: Minimum intensity that the colormap covers. The intensity histogram is
            calculated and vmin is set on the position of the maximum.
        vmax: Maximmum intensity that the colormap covers. The intensity histogram is
            calculated and vmax is set to the value where the 
            intensity drops below 0.01 % of the maximum.
    """
    # axis
    assert axis == 0 or axis == 1, 'axis must be 0 or 1'

    # default arguments
    if 'cmap' not in kwargs:
        kwargs['cmap'] = 'jet'
    if 'aspect' not in kwargs:
        kwargs['aspect'] = 'auto'
    if 'origin' not in kwargs:
        kwargs['origin'] = 'lower'
    if 'interpolation' not in kwargs:
        kwargs['interpolation'] = 'none'
    if 'vmin' not in kwargs or 'vmax' not in kwargs:
        vmin, vmax = self._calculated_vmin_vmax()
        if 'vmin' not in kwargs:
            kwargs['vmin'] = vmin
        if 'vmax' not in kwargs:
            kwargs['vmax'] = vmax

    # vars
    self.__i = 0
    self.__xlim = None
    self.__ylim = None
    self.__ax   = None
    self.__temp = self.copy()
    if axis == 0:
        self._calculated_roll = np.array([0]*self.shape[1])
    else:
        raise NotImplementedError('axis = 1 nor implemented yet')

    # change keybindings
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    # lims
    if self.__xlim is None:
        self.__xlim = [min(self._x_centers), max(self._x_centers)]
    if self.__ylim is None:
        self.__ylim = [min(self._y_centers), max(self._y_centers)]

    # vlines hlines warning
    if vlines is not None or hlines is not None:
        raise NotImplementedError('vlines is not implemented yet')
    
    # core update function
    update_function = None
    if update_function is None:
        def _update(ss):
            if self.__i >= len(self):
                self.__i = len(self) - 1
            elif self.__i < 0:
                self.__i = 0

            plt.title(f'{self.__i}: {self._calculated_roll[self.__i]}')
            im = self.__temp.imshow(ax=self.__ax, verbose=False)

            if axis == 0:
                self.__ax.plot([self.x_centers[self.__i], self.x_centers[self.__i]], [self.__ylim[0], self.__ylim[1]], color='white')
            else:
                raise NotImplementedError('not implemented yet')
            
            # if vlines is not None:
            #     if isinstance(vlines, Iterable) == False:
            #         vlines = [vlines, ]
            #     for vline in vlines:
            #         self.__ax.plot([vline, vline], [self.__ylim[0], self.__ylim[1]], color='red', ls='--', lw=1)
            # if hlines is not None:
            #     if isinstance(hlines, Iterable) == False:
            #         hlines = [hlines, ]
            #     for hline in hlines:
            #         self.__ax.plot([self.__xlim[0], self.__xlim[1]], [hline, hline], color='red', ls='--', lw=1)

            if self.__xlim is not None:
                plt.xlim(self.__xlim)
            if self.__ylim is not None:
                plt.ylim(self.__ylim)

    # keyboard events
    def keyboard(event, ss):
        # print(event.key)
        # print('keyboard')
        # print(event.key)
        if axis == 0:
            if event.key == 'up':
                self._calculated_roll[self.__i] += 1

                if axis == 0:
                    rolls = np.zeros(self.shape[1])
                    rolls[self.__i] = 1
                    self.__temp.set_roll(rolls)
                _update(ss)
                plt.draw()
            elif event.key == 'down':
                self._calculated_roll[self.__i] -= 1

                if axis == 0:
                    rolls = np.zeros(self.shape[1])
                    rolls[self.__i] = -1
                    self.__temp.set_roll(rolls)
                _update(ss)
                plt.draw()
            elif event.key == 'right':
                self.__i = self.__i + 1

                # self.__ax.cla()
                self.__ax.lines.clear()
                _update(ss)
                plt.draw()
            elif event.key == 'left':
                self.__i = self.__i - 1

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
# br.Image.roll_plot = roll_plot

# %% Other
from IPython import get_ipython
def is_notebook():
    """[EXPERIMENTAL] Return True if running from a ipython interactive terminal or jupyter notebook."""
    # assert IPythonok, 'is_notebook() cannot check for notebook\nError: python package `IPython` not found\nmaybe install it via ``pip install IPython``' 
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return True  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

