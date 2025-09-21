#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for everyday use ---> array and iterables"""

# %% ------------------------- Standard Imports --------------------------- %% #
from collections.abc import Iterable
import numpy as np

# %% ======================== calculation/modification =================== %% #
def sort(ref, *args):
    """Returns sorted arrays based on a reference array.

    Args:
        ref (list or array): 1D array.
        *args: array to sort.

    Returns:
        sorted arrays."""

    s = []
    for x in args:
        s.append( [x1 for (y,x1) in sorted(zip(ref,x), key=lambda pair: pair[0])])
    if len(args) == 1:
        s = s[0]
    return s

def moving_average(x, n):
    """Returns the moving average of an array.

    Args:
        x (list or array): 1D array.
        n (int): number of points to average.

    Returns:
        array of length given by (len(x)-n+1).

    Example:
        >>> x = [0,1,2,3,4,5,6,7,8,9]
        >>> print(am.moving_average(x, 1))
        [0. 1. 2. 3. 4. 5. 6. 7. 8. 9.]
        >>> print(am.moving_average(x, 2))
        [0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]
        >>> print(am.moving_average(x, 3))
        [1. 2. 3. 4. 5. 6. 7. 8.]
        >>> print(am.moving_average(x, 4))
        [1.5 2.5 3.5 4.5 5.5 6.5 7.5]

    """
    if n < 1:
        raise ValueError('n must be a positive integer (> 1).')
    if isinstance(n, int) == False:
        if n.is_integer() == False:
            raise ValueError('n must be a positive integer (> 1).')
    assert n <= len(x), f'number of points to average (n={n}) must be equal or less than the length of the array ({len(x)})'

    x = np.array(x)
    window = np.ones(int(n))/float(n)

    return np.convolve(x, window, 'valid')

def derivative(x, y, order=1):
    """Returns the derivative of y-coordinates as a function of x-coordinates.

        Args:
            x (list or array): 1D array x-coordinates.
            y (list or array): 1D array y-coordinates.
            order (number, optional): derivative order.

        Returns:
            x and y arrays.
    """

    # check order   
    if order<0:
        raise ValueError('order must be a positive integer.')
    
    # check monotonicity
    if check_monotonicity(x) < 0:
        raise ValueError('x array must be monotonically increasing. Please use brixs.fix_monotonocity().')


    x = np.array(x)
    y = np.array(y)

    x_diff = np.diff(x)
    y_diff = np.diff(y)/x_diff
    for i in range(order-1):
        y_diff = np.diff(y_diff)/x_diff[:len(x_diff)-(i+1)]

    for i in range(order):
        x = moving_average(x, n=2)

    return x, y_diff

def shifted(x, y, value, mode='hard'):
    """Shift (x, y) data.

    Args:
        x (list or array): 1D array.
        y (list or array): 1D array.
        value (float or int): shift value.
        mode (string, optional):
            #. ``mode='x'`` or ``mode='hard'``
                y is fully preserved while x is shifted.
            #. ``mode='y'``, ``'interp'``, or ``'soft'``
                x is preserved while y is interpolated with a shift
            #. ``mode='roll'` or ``r``,
                x and y are preserved and y elements are just rolled along the
                array (in this case ``shift`` value must be an integer).

    Returns:
        Shifted x and y.

    Warning:
        It is always better to use ``mode='hard'`` or ``'roll'`` since the form of y is fully
        preserved (no interpolation). After applying a shift using the ``mode='interp'``,
        one can apply a
        'inverse' shift to retrieve the original data. The diference between the retrieved
        y data and the original data will give an ideia of the information loss
        caused by the interpolation.
    """
    x = np.array(x)
    y = np.array(y)

    if value == 0:
        return x, y

    if mode == 'y' or mode == 'interp' or mode=='soft':
        if check_monotonicity(x) != 1:
            raise ValueError('x array must be increasingly monotonic.')
        y = np.interp(x, x + value, y)

    elif mode == 'x' or mode == 'hard':
        x = np.array(x) + value

    elif mode == 'roll' or mode == 'r':
        if is_integer(value):
            y = np.roll(y, int(value))
        else:
            raise ValueError("value must be an interger for mode='roll'.")
        # try:
        #     if value.is_integer():
        #         y = np.roll(y, int(value))
        #     else:
        #         raise ValueError("value must be an interger for mode='roll'.")
        # except AttributeError:
        #     y = np.roll(y, int(value))
        # if value > 0:
        #     y[:int(value)] = y[int(value)]
        # elif value < 0:
        #     # print(y[int(value):])
        #     # print( y[int(value-1)])
        #     y[int(value):] = y[int(value-1)]
    else:
        raise ValueError("mode not recognized (valid: 'y', 'x', 'roll').")

    return x, y

def transpose(array):
    """Returns transposed lists/arrays.

    Differently from numpy.transpose(), this function also transposes 1D arrays/lists.
    """
    try:
        row_count, col_count = np.shape(array)
        return np.transpose(array)
    except ValueError:
        return [[x] for x in array]

def mask(x, mask):
    """Returns a reduced array based on mask.

    Usage:
        mask('ABCDEF', [1, 0, 1, 0, 1, 1]) --> A C E F

    Args:
        x (Iterable): array
        mask (list): list with bools (True, False, ...)
        
    Returns:
        reduced array
    """
    return [d for d, s in zip(x, mask) if s]

def solve_linear_system(x1, y1, x2, y2):
	"""Return m and b where, y = m*x + b

	Args:
    	x1, y1, x2, y2 (numbers): one needs two x, y pairs
        	to solve the linear system.

	returns:
    	m, b
	"""
	m = (y2-y1)/(x2-x1)
	b = y1 - m*x1
	return m, b


# %% ============================= array check =========================== %% #
# def index(x, value, closest=True):
#     """Returns the first index of the element in array.

#     Args:
#         x (list or array): 1D array.
#         value (float or int): value.
#         closest (book, optional): if True, returns the index of the element in 
#             array which is closest to value.

#     Returns:
#         index (int)
#     """
#     # backpack developers note!!!!
#     # if this function changes, it needs to be copied to these files: figmanip

#     if closest:
#         # return int(np.argmin(np.abs(  np.array(x)-value)   ))
#         _inner1 = np.array(x) - value
#         _inner2 = np.ma.masked_array(_inner1, np.isnan(_inner1))
#         return int(np.argmin(np.abs(_inner2)))
#     else:
#         return np.where(x == value)[0]

def index(x, value, closest=True, roundup=False):
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

def all_equal(array):
    """Returns True if all elements of an array are equal."""
    iterator = iter(array)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)

def has_duplicates(arr, max_error=0):
    """Returns True if given array contains any duplicates.
    
    Args:
        arr (iterable): list of numbers
        max_error (number, optional): max error allowed to consider two numbers
            equal in percentage of the average value of the array. If zero,
            numbers are requires to be "equal" to be considered duplicated. 
            Default is zero.

    Returns:
        bool        
    """

    if max_error == 0:
        if len(arr) == len(set(arr)):
            return False
        else:
            return True
    else:
        if sum(abs(np.diff(arr))) > np.mean(arr)*max_error/100:
            return True
        return False
    
def remove_duplicates(array):
    """Returns the sorted unique elements of an array. Wrapper for np.unique()"""
    temp = np.unique(array)
    if isinstance(array, list):
        temp = list(temp)
    return temp

# %% ============================ monotonicity =========================== %% #
def check_monotonicity(array):
    """return 1 (-1) if increas. (decre.) monotonic or 0 if not monotonic."""
    if np.all(np.diff(array) > 0) == True:
        return 1
    elif np.all(np.diff(array) < 0) == True:
        return -1
    else:
        return 0

def fix_monotonicity(x, y, mode='increasing'):
    """return x, y where the x array is monotonically increasing or decreasing."""
    if mode != 'increasing' and mode != 'decreasing':
        raise ValueError('mode should be "decreasing" or "increasing".')

    if mode == 'increasing':
        if check_monotonicity(x) == 1:
            return x, y
        else:
            unqa, ID, counts = np.unique(x, return_inverse=True, return_counts=True)
            return unqa, np.bincount(ID, y)/counts
    if mode == 'decreasing':
        if check_monotonicity(x) == -1:
            return x, y
        else:
            unqa, ID, counts = np.unique(x, return_inverse=True, return_counts=True)
            return np.fliplr(unqa), np.fliplr(np.bincount(ID, y)/counts)



    if self.monotonicity != mode:
        temp = sort(self.x, self.x, self.y)
        if mode == 'increasing':
            self._x = np.array(temp[0])
            self._y = np.array(temp[1])
        if mode == 'decreasing':
            self._x = np.fliplr(temp[0])
            self._y = np.fliplr(temp[1])
        self._calib = -self.calib
        if self.step is not None:
            self._step = -self.step

        self.check_monotonicity()

# %% ============================== extract ============================== %% #
def choose(x, limits):
    """Return a mask of x values inside range pairs.

    Args:
        x (list or array): 1d array.
        limits (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x. [[xi, xf], [xi_2, xf_2], ...]

    Raises:
        RuntimeError: if limits have overlapping intervals

    Returns:
        1d list.
    """
    # backpack developers note!!!!
    # if this function changes, it needs to be copied to these files: figmanip

    assert isinstance(x, Iterable), 'input must be a iterable'
    x = np.array(x)

    try:  # ((start, end), )

        # check overlapping limits
        for i, lim in enumerate(limits):
            # put start as the lower number
            if lim[0] > lim[1]:
                xf1, xi1 = lim
            else:
                xi1, xf1 = lim
            for j, lim2 in enumerate(limits): 
                # put start as the lower number
                if lim2[0] > lim2[1]:
                    xf2, xi2 = lim2
                else:
                    xi2, xf2 = lim2
                if i != j:
                    if (xf1 >= xi2 and xf1 <= xf2):
                        raise RuntimeError(f'overlapping intervals [{(xi1, xf1)}] and [{(xi2, xf2)}]')

        choose_range = [None]*len(limits)
        for i, lim in enumerate(limits):
            # put start as the lower number
            if lim[0] > lim[1]:
                x_final, x_init = lim
            else:
                x_init, x_final = lim
            choose_range[i] = np.logical_and(x>=x_init, x<=x_final)
        choose_range = [True if x == 1 else False for x in np.sum(choose_range, axis=0)]
    except TypeError:  # (start, end) 
        # put start as the lower number
        if limits[0] > limits[1]:
            x_final, x_init = limits
        else:
            x_init, x_final = limits
        choose_range = np.logical_and(x>=x_init, x<=x_final)
    return choose_range

def extract(x, y, limits, invert=False):
    """Returns specific data limits from x and y.

    Args:
        x (list or array): 1D reference vector.
        y (list or array): 1D y-coordinates or list of several data sets.
        limits (list): a pair of values or a list of pairs. Each pair represents
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
    # backpack developers note!!!!
    # if this function changes, it needs to be copied to these files: figmanip
    
    x = np.array(x)
    y = np.array(y)

    # if data is all inside limits, then nothing is done
    if len(limits) == 1:
        if limits[0][0] <= min(x) and limits[0][1] >= max(x):
            if invert:
                return np.array([]), np.array([])
            else:
                return x, y

    choose_range = choose(x, limits)
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
        return [], []
        # raise RuntimeError('No data points within the selected range.')

# %% ========================== Experimental ============================= %% #
def array2list(arr):
    """recursively transform a numpy array to a list"""
    final = []
    if isinstance(arr, Iterable) and isinstance(arr, str) == False:
        for item in arr:
            if isinstance(item, Iterable):
                final.append(array2list(item))
            else:
                final.append(item)
        return final
    else:
        return arr
    
def flatten(x):
    """[EXPERIMENTAL] Return a copy of the array collapsed into one dimension."""
    assert isinstance(x, Iterable), 'input must be Iterable (list, tuple, array)'

    # empty array
    if len(x) == 0:
        return x

    # if len(np.array(x).shape) == 1:
    #     return x
    return np.concatenate(x).ravel()

def digitize(x, bins, xmin=None, xmax=None):
    """Returns dictionary with bin values (keys) and indexes. Wrapper for np.digitize and np.histogram_bin_edges)
    
    bin values is the central value between two bin edges and the indexes 
    indicate the indexes of x items that goes in that bin

    Args:
        x (list): list with values to place in the bins.
        bins (int, str, list): If bins is a list, it is assumed to be the bins edges. Note that,           
            bins is inclusive at the bin_edge with smaller value, 
            and exclusive at the bin_edge with higher value.
        
            If bins is an int (from np.histogram_bin_edges), it defines the number of equal-width 
            bins in the given range (10, by default). If bins is a sequence, it defines the 
            bin edges, including the rightmost edge, allowing for non-uniform bin widths.

            If bins is a string from the list below (from np.histogram_bin_edges), 
            histogram_bin_edges will use the method 
            chosen to calculate the optimal bin width and consequently the number of bins 
            (see the Notes section for more detail on the estimators) from the data that 
            falls within the requested range. While the bin width will be optimal for the 
            actual data in the range, the number of bins will be computed to fill the entire 
            range, including the empty portions. For visualization, using the 'auto' option 
            is suggested. Weighted data is not supported for automated bin size selection.

            'auto'
            Minimum bin width between the 'sturges' and 'fd' estimators. Provides good 
            all-around performance.

            'fd' (Freedman Diaconis Estimator)
            Robust (resilient to outliers) estimator that takes into account data variability 
            and data size.

            'doane'
            An improved version of Sturges' estimator that works better with non-normal datasets.

            'scott'
            Less robust estimator that takes into account data variability and data size.

            'stone'
            Estimator based on leave-one-out cross-validation estimate of the integrated squared error. 
            Can be regarded as a generalization of Scott's rule.

            'rice'
            Estimator does not take variability into account, only data size. Commonly overestimates 
            number of bins required.

            'sturges'
            R's default method, only accounts for data size. Only optimal for gaussian data and 
            underestimates number of bins for large non-gaussian datasets.

            'sqrt'
            Square root (of data size) estimator, used by Excel and other programs for its speed and simplicity.
        xmin, xmax (int, optional): minimum/maximum value for bins (only used if bins are not explicitly 
            given. Default is None. If None, min and max array values are used. 
            
        Returns:
            indexes (dict), bin_edges (list)
    """
    if isinstance(bins, Iterable)==False or isinstance(bins, str):
        if xmin is None: xmin = np.min(x)
        if xmax is None: xmax = np.max(x)
        assert xmax > xmin, 'xmax must be higher than xmin'
        bin_edges = np.histogram_bin_edges(x, bins=bins, range=(xmin, xmax), weights=None)
    else:
        assert xmin is None, 'if bins is a list (assumed to be the bins edges), than xmin must be None'
        assert xmax is None, 'if bins is a list (assumed to be the bins edges), than xmax must be None'
        bin_edges = bins
    temp = np.digitize(x, bins=bin_edges, right=False)
    
    # building indexes dictionary
    indexes = {}
    for i in range(len(bin_edges[:-1])):
        # indexes[(bin_edges[i]+bin_edges[i+1])/2] = [_ for _ in np.argwhere(temp==i+1)]
        indexes[(bin_edges[i]+bin_edges[i+1])/2] = [int(_[0]) for _ in list(np.argwhere(temp==i+1))]
        
    return indexes, bin_edges
