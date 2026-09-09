# %% ========================== Standard Imports ========================= %% #
from pathlib import Path
import numpy as np

# %% =============================== brixs =============================== %% #
import brixs as br
import brixs.addons.centroid

# %% ========================== Special Imports ========================== %% #
try:
    from scipy.ndimage import median_filter, uniform_filter1d
    from scipy.ndimage import maximum_filter, uniform_filter
except:
    pass

# %% ========================= new image methods ========================= %% #
def _detect_outliers_via_row_median(self, contrast_ratio=8, slope=0.030755125, nrows2bin=2, median_filter_window=100,
                                    moving_average_window=20, minimum_intensity=600, minimum_average=400):
    # 1. Curvature correction
    im1 = self.set_vertical_shift_via_polyval(p=[-slope, 0])

    # 2. Bin vertically
    nrows = np.rint(2048 / nrows2bin)
    im2 = im1.binning(ncols=None, nrows=nrows)

    # 3. Extract binned rows into a 2D NumPy array
    # ss = im2.get_rows(max_number_of_rows=2048)
    # data = np.asarray([s.y for s in ss], dtype=np.float32)
    # nrows, ncols = data.shape

    # 4. Median background removal, vectorized over all rows
    _im3 = median_filter(im2.data, size=(1, median_filter_window), mode="reflect")
    _im4 = np.abs(im2.data - _im3)

    # 5. Moving average normalization, vectorized
    _im5 = uniform_filter1d(_im4, size=moving_average_window, axis=1, mode="reflect")
    eps = np.finfo(np.float32).eps  # This prevents dividing by zero
    _im6 = _im4/np.maximum(_im5, eps)

    # 6. Candidate mask
    candidate_mask = _im6 > contrast_ratio

    # Remove rows affected by slope correction
    edge_rows = int((abs(slope) * 2048 + 10) / nrows2bin)
    edge_rows = min(edge_rows, nrows)
    if slope > 0:
        # invalid region at bottom
        candidate_mask[nrows - edge_rows:, :] = False
    elif slope < 0:
        # invalid region at top
        candidate_mask[:edge_rows, :] = False

    # get rows and cols of all candidates
    candidate_rows, candidate_cols = np.where(candidate_mask)

    if len(candidate_rows) == 0:
        pe1 = br.PhotonEvents()
        accepted_mask = None
    else:
        # 8. Convert candidate coordinates from binned image to original image
        y_indices = np.rint(candidate_rows*nrows2bin).astype(int)
        x_indices = candidate_cols.astype(int)

        # precompute local max and local average images.
        local_size = 3
        im_data    = np.asarray(im1.data, dtype=np.float32)
        local_max  = maximum_filter(im_data, size=local_size, mode="reflect")
        local_mean = uniform_filter(im_data, size=local_size, mode="reflect")

        # Validation at all candidate positions at once.
        accepted_mask = ((local_max[y_indices, x_indices] > minimum_intensity) | (local_mean[y_indices, x_indices] > minimum_average))
        accepted_x = x_indices[accepted_mask]
        accepted_y = y_indices[accepted_mask]

        # 9. Build PhotonEvents
        pe1 = br.PhotonEvents()
        for x, y in zip(accepted_x, accepted_y):
            pe1.append(x=float(x), y=float(y))

        # 10. Undo curvature correction
        pe1 = pe1.set_vertical_shift_via_polyval(p=[slope, 0]).crop(y_start=0, y_stop=2046)

    return {'im1':im1, 'im2':im2, '_im3':_im3, '_im4':_im4, '_im5':_im5, '_im6':_im6, 'pe1':pe1}
br.Image._detect_outliers_via_row_median = _detect_outliers_via_row_median

def detect_outliers_via_row_median(*args, **kwargs):
    """Return list of detected outliers using median-based background removal.

    median_filter_window -> "blurr image". If it is larger it can wipe out 
        sharp features. Increase it if you are detecting photon events as cosmic rays.
    contrast_ratio -> increse contrast_ratio to be less sensitive to outliers.

    Processing steps::

        - Input image (im0)
            Curvature correction (im0, slope -> im1)
            └─ vertical binning (im1, nrows2bin -> im2)
                └─ Calculate median of each row (im2, median_filter_window -> _im3)
                    └─ [im2-_im3]: median filter subtraction (im2, _im3 -> _im4)
                        └─ Calculate moving average of each row (_im4, moving_average_window -> _im5)
                            └─ [_im4/_im5]: Divide by moving average (_im4, _im5 -> _im6)
                                └─ Define candidates: pixels that are higher than contrast_ratio (_im6, contrast_ratio -> candidate_mask)
                                    └─ Check candidates (_im6, candidate_mask, minimum_intensity minimum_average -> pe1)

    Description of steps:
        1) The image is curvature-corrected to align isoenergetic lines .
        2) The image is vertically binned (`nrows2bin`) to improve statistics.
        3) Each row is smoothed using a median filter (`median_filter_window`)
            to remove small outliers
        4) The median-filtered signal is subtracted from the original row to
        obtain the residual:
            residual = |y - median_filter(y)|
        5) A moving average is computed along the rows using `moving_average_window` as a 
            measurement of the local average intensity.
        6) The residual is normalized by its moving average to enhance
            points that exceed the local average intensity.
        7) Pixels in the 'enhanced residual image' with values above 
            `contrast_ratio` are marked as candidates.
        8) Each candidate is validated using the original image:
            - If any neighboring pixel exceeds `minimum_intensity`, it is
                accepted as an outlier.
            - If the average of neighboring pixels exceeds `minimum_average`,
                it is also accepted.
            - Otherwise, the candidate is rejected.

    Args:
        contrast_ratio (float, optional): Minimum factor by which a pixel must
            exceed its local background to be considered an outlier candidate.
            Default is 8.
        slope (float, optional): Slope used for curvature correction of the image.
            Default is 0.030755125.
        nrows2bin (int, optional): Number of rows to bin together.
            Default is 2.
        median_filter_window (int, optional): Window size of the median filter
            used to estimate the background.
            Default is 100.
        moving_average_window (int, optional): Window size for the moving average
            used for normalization.
            Default is 20.
        minimum_intensity (float, optional): Minimum neighboring pixel intensity
            required to accept a candidate. Default is 600.
        minimum_average (float, optional): Minimum neighboring average intensity
            required to accept a candidate. Default is 400.

    Returns:
        br.PhotonEvents: Collection of detected outlier events.
    """
    return _detect_outliers_via_row_median(*args, **kwargs)['pe1']
br.Image.detect_outliers_via_row_median = detect_outliers_via_row_median

def _detect_outliers_via_threshold(self, threshold):
    """Return outliers detected using a simple intensity threshold.

    This method identifies pixels in the image whose intensity exceeds a given
    threshold.

    Args:
        threshold (float): Intensity threshold above which a pixel is considered
            an outlier.

    Returns:
        br.PhotonEvents: Collection of detected outlier events with coordinates
        corresponding to the pixel centers in the image.
    """
    _pos, _ = self.get_positions_above_threshold(threshold=threshold, nx=0, ny=0, coordinates='centers')

    # transform it in photon events
    _temp = np.array(_pos)
    return br.PhotonEvents(x=_temp[:, 1], y=_temp[:, 0])
br.Image.detect_outliers_via_threshold = _detect_outliers_via_threshold

# %% ============================= OBSOLETE ============================== %% #
def _median_filter_1d(y, k):
    """Apply a 1D median filter to a signal.
    
    This is an alternative to scipy.ndimage function median_filter(). The scipy
    version is faster and should be prefered and scipy is installed.

    This function smooths a 1D array by replacing each point with the median
    value of its surrounding window. It is robust to outliers and effectively
    suppresses spikes while preserving the overall shape of the signal.

    Processing steps:
        1) The input signal is padded at both ends using edge values to handle
           boundary conditions.
        2) A sliding window of size `k` is constructed over the padded signal.
        3) The median value within each window is computed.
        4) The resulting array of medians forms the filtered signal.

    Args:
        y (array-like): Input 1D signal to be filtered.
        k (int): Size of the sliding window. Must be a positive odd integer for
            symmetric filtering.

    Returns:
        np.ndarray: Filtered signal with the same length as the input.

    Notes:
        - Larger values of `k` result in stronger smoothing but may suppress
          narrow features.
        - Smaller values of `k` preserve fine detail but may not remove all
          outliers.
    """

    k = int(k)
    pad = k // 2
    y_pad = np.pad(y, pad, mode='edge')

    windows = np.lib.stride_tricks.sliding_window_view(y_pad, k)
    return np.median(windows, axis=-1)
    
def _detect_outliers_via_row_fitting(self, contrast_ratio=8, slope=0.030755125, nrows2bin=2,
                                    median_filter_window=20, maximum_contrast_for_polyfit=0.1,
                                    moving_average_window=20, minimum_intensity=600, 
                                    minimum_average=400):
    """Return list of detected outliers.

    Processing steps:
        1) The image is curvature-corrected to align isoenergetic lines.
        2) The image is vertically binned (every two rows) to improve statistics 
            and enhance contrast.
        3) Each row is smoothed using a median filter (`median_filter_window`) 
            to suppress spikes and noise.
        4) A contrast metric is computed as:
               contrast = (max(y) - median(y)) / median(y)
            This is used to decide the fitting model:
            - If contrast > `maximum_contrast_for_polyfit`, the row is considered
                to have a peak-like (Gaussian-like) structure.
            - Otherwise, the row is treated as background/noise-dominated.
        5) The row is fitted accordingly:
            - Gaussian-like model (asymmetric peak + linear background) for 
                high-contrast rows.
            - Polynomial fit for low-contrast rows.
        6) The fitted curve is subtracted from the row to remove the underlying 
            sample signal:
                residua = row - fit
        7) The absolute value of the residual is taken to avoid negative values.
        8) A moving average is computed using `moving_average_window`.
        9) The residual signal is normalized (divided) by its moving average, providing a 
            relative measure of how many times a pixel exceeds its local background.
        10) Pixels with values above `contrast_ratio` are marked as candidates.
        11) Each candidate is validated using the original (non-processed) image:
            - If any neighboring pixel exceeds `minimum_intensity`, it is flagged 
                as an outlier.
            - If the average of the first neighboring pixels exceeds 
                `minimum_average`, it is also flagged as an outlier.
            - Otherwise, the candidate is rejected.

    Args:
        contrast_ratio (float, optional): Minimum factor by which a pixel must
            exceed its local background (moving average) to be considered an
            outlier candidate. Default is 8.
        slope (float, optional): Slope used for curvature correction of the image.
            Determines the vertical shift applied during alignment.
            Default is 0.030755125.
        nrows2bin (int, optional): Number of rows to bin together when increasing
            statistics. Default is 2.
        median_filter_window (int, optional): Window size used for the median
            filter applied to each row before estimating contrast. This helps
            suppress spikes and outliers when deciding the fitting model.
            Default is 20.
        maximum_contrast_for_polyfit (float, optional): Threshold on the contrast
            metric used to decide the fitting model. Rows with contrast above this
            value are fitted with a Gaussian-like model; otherwise, a polynomial
            fit is used. Default is 0.1.
        moving_average_window (int, optional): Number of neighboring pixels used
            to compute the local moving average (background estimate).
            Default is 20.
        minimum_intensity (float, optional): Minimum absolute intensity required
            in any neighboring pixel for a candidate to be accepted as an outlier.
            Default is 600.
        minimum_average (float, optional): Minimum average intensity of neighboring
            pixels required for a candidate to be accepted as an outlier if the
            individual pixel threshold is not met. Default is 400.

    Returns:
        br.PhotonEvents: Collection of detected outlier events.
    """

    # fix curvature (isoenergetic lines must be in the same pixel row)
    im = self.set_vertical_shift_via_polyval(p=[-slope, 0])

    # binning every 2 rows to increase statistics
    im2 = im.binning(ncols=None, nrows=int(2048/nrows2bin))

    # get rows
    ss = im2.get_rows(max_number_of_rows=2048)

    # loop each row and detect outliers that pop above 
    c1 = br.PhotonEvents()
    c2 = br.PhotonEvents()
    fit = br.Spectra()
    divided = br.Spectra()
    for i, s in enumerate(ss):
        if i < (2048-slope*2048-10)/nrows2bin:  # remove bottom edge because of slope correction
            # reduce noise for determining if fit should be gaussian or poly
            y2 = median_filter(s.y, size=median_filter_window)
            baseline = np.median(y2)
            contrast = (np.max(y2) - baseline) / baseline
            if contrast > maximum_contrast_for_polyfit: gaussian = True
            else: gaussian = False

            # if scipy is not available
            # y2 = median_filter_1d(s.y, median_filter_window)
            # if contrast > maximum_contrast_for_polyfit: gaussian = True
            # else: gaussian = False

            # peak fit
            s.model.clear()
            if gaussian:
                s.model.peaks.add(c=1100, w=500, amp=max(s.y), m=0, m_vary=True, w_max=2000)

                s.model.asymmetricpeaks.add(c=1100, w1=500, w2=500, amp=max(s.y), m1=0, m2=0, 
                                            m1_vary=True, m2_vary=True, w1_min=100, w1_max=2000,
                                            w2_min=100, w2_max=2000)
                s.model.linear.add(linear=0, const=min(s.y))
                _ = s.model.fit(limits=(0, 2048))
                _fit = br.Spectrum(x=s.x, y=s.model.calculate_spectrum(x=s.x))
            else:
                _temp = s.polyfit(deg=10)
                _fit = br.Spectrum(x=s.x, y=_temp['model'](x=s.x))

            # subtract the fitted curve from the row to get a flat signal
            # this removes the signal from the sample
            _raw = s - _fit

            # remove negative values
            temp = br.Spectrum(x=_raw.x, y=np.abs(_raw.y))

            # divide the signal by the average of the 20 points surrounding each
            # point to get a relative signal that pops above the noise
            temp2 = (temp/temp.moving_average(moving_average_window).interp(x=temp.x))
            # divided.append(temp2)

            # find points that pop `contrast_ratio`` times above the background
            # if the any first neighboring pixel has a signal above `minimum_intensity`, it's a cosmic ray
            # if the average of the first neighboring pixels is above `minimum_average`, it's a cosmic ray
            # otherwise, it's not a cosmic ray
            itemindex = np.where(temp2.y > contrast_ratio)
            y = im2.y_centers[i] - 0.5  # because of the binning, the y value is at the center of the 2 rows, so we need to subtract 0.5 to get the original y value
            for _x in itemindex[0]:
                spot = im.get_spot(x=_x, y=y, nx=1, ny=1)
                if (spot[0].data > minimum_intensity).any():
                    c1.append(x=_x, y=y)
                elif spot[0].calculate_average() > minimum_average:
                    c1.append(x=_x, y=y)
                else:
                    c2.append(x=_x, y=y)
        
    # fix back the curvature
    c1 = c1.set_vertical_shift_via_polyval(p=[slope, 0]).crop(y_stop=2046)
    # c1.rejected = c2.set_vertical_shift_via_polyval(p=[slope, 0]).crop(y_stop=2046)

    # save spectra for later verification
    # c1.ss = ss
    # c1.fit = fit
    # c1.divided = divided

    return c1
# br.Image.detect_outliers_via_row_fitting = _detect_outliers_via_row_fitting
