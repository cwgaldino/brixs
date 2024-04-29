# SPECTRUM
def find_peaks(self, prominence=None, width=4, moving_average_window=4, ranges=None):
    """Find peaks. Wrapper for `scipy.signal.find_peaks()`_.

    Sets :py:attr:`peaks` attribute.

    Args:
        prominence (number, optional): minimum prominence of peaks in percentage
            of the maximum prominence [max(y) - min(y)]. Default is 5.
        width (number, optional): minimum number of data points defining a peak.
        moving_average_window (int, optional): window size for smoothing the
            data for finding peaks. Default is 4.
        ranges (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x. Use None to indicate
            the minimum or maximum x value of the data.

    Returns:
        None

    .. _scipy.signal.find_peaks(): https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
    """
    if ranges is None:
        s = self
    else:
        s = self._extract(ranges)

    # check monotonicity
    if s.monotonicity is None:
        s.check_monotonicity()

    # step
    if s.step is None:
        s.check_step()

    self.peaks.find(x=s.x, y=s.y, prominence=prominence, width=width, moving_average_window=moving_average_window)

def fit_peaks(self, method='least_squares', ranges=None):
    """Fit peaks. Wrapper for `lmfit.minimize()`_.

    Args:
        yerr (array, optional): data uncertainty. 
        method (str, optional): Name of the fitting method to use. See methods
            available on `lmfit.minimize()`_ documentation.
        ranges (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x. Use None to indicate
            the minimum or maximum x value of the data.

    Returns:
        None

    .. _lmfit.minimize(): https://lmfit.github.io/lmfit-py/fitting.html
    """
    if ranges is None:
        self.check_monotonicity()
        x = self.x
        y = self.y
    else:
        ranges = self._validate_ranges(ranges)
        s = self._extract(ranges)
        s.check_monotonicity()
        x = s.x
        y = s.y

    # check if peaks is defined
    if len(self.peaks) == 0:
        raise ValueError('No peaks to fit.\nRun Spectrum.find_peaks() or set s.peaks manually.')

    # fit
    return self.peaks.fit(x, y, method=method, update_peaks=True)

# SPECTRA
def copy_peaks_from_spectra(self):
    """Copy peaks from each spectrum.

    Returns:
        None
    """
    self.peaks._copy_from_spectra(self)

def copy_peaks_to_spectra(self):
    """Copy peaks to each spectrum

    Returns:
        None
    """
    self.peaks._copy_to_spectra(self)

def find_peaks(self, prominence=None, width=4, moving_average_window=8, ranges=None):
    """Find peaks recursively. Wrapper for `scipy.signal.find_peaks()`_.

    Args:
        prominence (number, optional): minimum prominence of peaks in percentage
            of the maximum prominence [max(y) - min(y)]. Default is 5.
        width (number, optional): minimum number of data points defining a peak.
        moving_average_window (int, optional): window size for smoothing the
            data for finding peaks. Default is 4.
        ranges (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x. Use None to indicate
            the minimum or maximum x value of the data.

    Returns:
        None

    .. _scipy.signal.find_peaks(): https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
    """
    for s in self:
        s.find_peaks(prominence=prominence, width=width, moving_average_window=moving_average_window, ranges=ranges)
    self.copy_peaks_from_spectra()

def fit_peak(self, asymmetry=False, moving_average_window=1, method='least_squares', ranges=None, verbose=False):
    """Fit one peak for all spectra recursively. Wrapper for `lmfit.minimize()`_.

    Args:
        asymmetry (bool or dict, optional): if True, fits each half of the
            with a different width. Default is False
        moving_average_window (int, optional): window size for smoothing the
            data for finding the peak. Default is 1.
        method (str, optional): Name of the fitting method to use. See methods
            available on `lmfit.minimize()`_ documentation.
        ranges (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x. Use None to indicate
            the minimum or maximum x value of the data.

    Returns:
        None

    .. _lmfit.minimize(): https://lmfit.github.io/lmfit-py/fitting.html
    """
    for i, s in enumerate(self):
        if verbose:
            print(f'{i+1}/{len(self)}')
        s.fit_peak(asymmetry=asymmetry, moving_average_window=moving_average_window, method=method, ranges=ranges)
    self.copy_peaks_from_spectra()

    # THIS IS FOR FITTING ONE PEAK SIMULTANOUSLY
    # self.peaks.clear()

    # if ranges is None:
    #     temp = self
    # else:
    #     temp = self._extract(ranges=ranges)
    # xs = [s.x for s in temp]
    # ys = [s.y for s in temp]

    # for i2 in range(len(self)):
    #     # guess amp and c
    #     amp = max(ys[i2])
    #     c   = xs[i2][np.argmax(ys[i2])]

    #     # guess fwhm
    #     try:
    #         w1 = xs[i2][np.argmax(ys[i2])] - xs[i2][:np.argmax(ys[i2])][::-1][arraymanip.index(ys[i2][:np.argmax(ys[i2])][::-1], max(ys[i2])/2)]
    #     except ValueError:
    #         w1 = xs[i2][np.argmax(ys[i2]):][arraymanip.index(ys[i2][np.argmax(ys[i2]):], max(ys[i2])/2)] - xs[i2][np.argmax(ys[i2])]
    #     try:
    #         w2 = xs[i2][np.argmax(ys[i2]):][arraymanip.index(ys[i2][np.argmax(ys[i2]):], max(ys[i2])/2)] - xs[i2][np.argmax(ys[i2])]
    #     except ValueError:
    #         w2 = xs[i2][np.argmax(ys[i2])] - xs[i2][:np.argmax(ys[i2])][::-1][arraymanip.index(ys[i2][:np.argmax(ys[i2])][::-1], max(ys[i2])/2)]
    #     w = w1 + w2
    #     if w <= 0:
    #         w = 0.1*(max(xs[i2])-min(xs[i2]))
    #     if w <= 0:
    #         w = 1

    #     # peaks
    #     if asymmetry:
    #         self.peaks.append(i2=i2, amp=amp, c=c, w1=w1, w2=w2)
    #     else:
    #         self.peaks.append(i2=i2, amp=amp, c=c, w=w)

    # self.peaks.fit(xs=xs, ys=ys, method=method)

def fit_peaks(self, method='least_squares', ranges=None):
    """Fit peaks for all spectra SIMULTANEOUSLY. Wrapper for `lmfit.minimize()`_.

    Args:
        method (str, optional): Name of the fitting method to use. See methods
            available on `lmfit.minimize()`_ documentation.
        ranges (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x. Use None to indicate
            the minimum or maximum x value of the data.

    Returns:
        None

    .. _lmfit.minimize(): https://lmfit.github.io/lmfit-py/fitting.html
    """
    if ranges is None:
        temp = self
    else:
        temp = self._extract(ranges=ranges)
    xs = [s.x for s in temp]
    ys = [s.y for s in temp]

    self.peaks.fit(xs=xs, ys=ys, method=method)
