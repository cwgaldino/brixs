#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Centroid algorithm"""

# %% ------------------------- Standard Imports -------------------------- %% #
import numpy as np
import math
# %%

# %% -------------------------- brixs Imports ---------------------------- %% #
import brixs as br
# %%

# %% ============================= support =============================== %% #
def get_positions_above_threshold(self, threshold, n=0, coordinates='centers'):
    """returns list of y, x positions where pixel intensity is above threshold

    Args:
        threshold (number): threshold value
        n (int, optional): If n>0, only the brightest pixel around the first 
            n neighbors of the (y, x) positions will trigger the threshold,
            e.g., if n=0 and all first neighbors around (y, x) are above the 
            threshold, then all positions (y-1, x-1), ..., (y+1, x+1), all 
            positions which intensity is above threshold will be considered.
            However, if n=1 and all first neighbors around (y, x) are above the 
            threshold, only the (yb, xb) with the brightest intensity amongst
            all (y-1, x-1), ..., (y+1, x+1) will be considered. This is helpfull
            when bright spot are more than one pixel wide. Default is n=0.
        coordinates (str, optional): If `pixels`, (y, x) is given in pixel 
            coordinates. If `centers`, (y, x) is given in terms of x_centers
            and y_centers. Default is `centers`.

    Returns:
        list of positions [(y1, x1), (y2, x2), ...]
    """
    assert coordinates in ['pixels', 'centers'], f'coordinates must be `pixels` or `centers`, not `{coordinates}`'
    assert isinstance(n, int), f'n must be of type integer, not {type(n)}'
    assert n >= 0, f'n must be an integer number higher or equal to 0'

    ys, xs = np.where(self.data > threshold)

    central = []
    for y, x in zip(ys, xs):
        if n > 0:
            _y, _x = self.get_brightest_pixel_position(y=y, x=x, n=n, coordinates='pixels')
            bright_pixel = (_y, _x)
            if bright_pixel not in central:
                central.append(bright_pixel)
        else:
            central.append((y, x))
    if coordinates == 'centers':
        central = [self.index2center(y=_[0], x=_[1]) for _ in central]
    return central
br.Image.get_positions_above_threshold = get_positions_above_threshold

def get_brightest_pixel_position(self, y, x, n, coordinates='centers'):
    """Returns the positions of the most intense pixel around y, x

    Args:
        y, x (tuple or list): (y, x) position in terms of pixel value (if 
            coordinates='pixels') or in terms of y_centers and x_centers (if
            coordinates='centers').
        n (int): number of neighbors to include, e.g., if n=1, only first neighbors 
            will be considered and the spot will be a 3x3 array.
        coordinates (str, optional): `pixels` if pixel positions are given in 
            pixel coordinates or `centers` if positions are given in terms of 
            y_centers and x_centers.
    
    Returns:
        y, x
    """
    assert coordinates in ['pixels', 'centers'], f'coordinates must be `pixels` or `centers`, not `{coordinates}`'
    assert isinstance(n, int), f'n must be of type integer, not {type(n)}'
    assert n >= 0, f'n must be an integer number higher or equal to 0'

    spot, yslice, xslice, is_near_edge = self.get_spot(y=y, x=x, n=n, coordinates=coordinates)
    _y, _x = spot.argmax(coordinates='pixels')
    
    y = int(yslice.start+_y)
    x = int(xslice.start+_x)

    if coordinates == 'centers':
        y, x = self.index2center(y=y, x=x)
    return y, x
br.Image.get_brightest_pixel_position = get_brightest_pixel_position

def get_center_of_mass(self, coordinates='centers'):
    """Returns the y, x position of the intensity center of mass of pixels

    Args:
        coordinates (str, optional): If `pixels`, (y, x) is given in pixel 
            coordinates. If `centers`, (y, x) is given in terms of x_centers
            and y_centers. Default is `centers`.

    Returns:
        y, x
    """
    if coordinates == 'centers':
        mx = np.average(self.x_centers, weights=self.data.sum(axis=0))
        my = np.average(self.y_centers, weights=self.data.sum(axis=1))
    else:
        mx = np.average(np.arange(0, self.shape[1]), weights=self.data.sum(axis=0))
        my = np.average(np.arange(0, self.shape[0]), weights=self.data.sum(axis=1))
    return my, mx
br.Image.get_center_of_mass = get_center_of_mass

def enhance(self, n):
    """multiplies image by itself and apply moving average

    This makes outliers stadout more.

    Args:
        n (int): Size of the moving average window, i.e., number of points to average.

    Returns:
        enhanced image
    """
    im = self.floor()
    return im.multiply(im).moving_average(n)
br.Image.enhance = enhance

def find_candidates(self, n, n2, threshold):
    """Returns a list of photon hit candidates defined by a threshold in a enhanced image

    Args:
        n (int): Size of the moving average window used in the process of
            enhancing the image, i.e., number of points to average.
        n2 (int): photon hits candidates that are within n2 pixels of distance 
            from each other will be considered the same candidate. n2 must
            be the pixel distance between a photon hit and the farthest excited
            pixel.
        threshold (number): threshold value in terms of pixel value of the 
            ENHANCED image (one might need to plot the enhanced image to get 
            an idea of the pixel intensity).
        
    Returns:
        photon events list
    """
    # enhance image
    # # image is floored, squared, and moving_averaged
    im = self.enhance(n=n)

    # find photon hit candidates from the enhanced image
    pos = im.get_positions_above_threshold(threshold, n=n2, coordinates='centers')
    
    # get brightst pixel around candidates in the original image
    pos2 = []
    for _pos in pos:
        y, x = self.get_brightest_pixel_position(y=_pos[0], x=_pos[1], n=n2, coordinates='centers')
        if (y, x) not in pos2:
            pos2.append((y, x))
    return br.PhotonEvents(x=[_[1] for _ in pos2], y=[_[0] for _ in pos2])
br.Image.find_candidates = find_candidates

def find_and_patch(self, n, n2, threshold, MAX_NUMBER_OF_CANDIDATES=200):
    """Returns image where photon event are patched out (cosmic rays correction)

    Args:
        n (int): Size of the moving average window used in the process of
            enhancing the image, i.e., number of points to average.
        n2 (int): photon hits candidates that are within n2 pixels of distance 
            from each other will be considered the same candidate. n2 must
            be the pixel distance between a photon hit and the farthest excited
            pixel.
        threshold (number): threshold value in terms of pixel value of the 
            ENHANCED image (one might need to plot the enahnced image to get 
            an idea of the pixel intensity).
        MAX_NUMBER_OF_CANDIDATES (int, optional): raises error if number of 
            photons to be patched out is larger than MAX_NUMBER_OF_CANDIDATES.
            Useful for preventing too low threshold as patching image is slow.
            Default is 10. 
        
    Returns:
        im, pe
        patched image and photon events list
    """
    pe = self.find_candidates(n, n2, threshold)
    if len(pe) > MAX_NUMBER_OF_CANDIDATES:
        raise ValueError(f'number of photon hit candidates ({len(pe)}) is higher than MAX_NUMBER_OF_CANDIDATES={MAX_NUMBER_OF_CANDIDATES}. Maybe try increasing MAX_NUMBER_OF_CANDIDATES, n2, or `threshold`')
    if len(pe) > 0:
        im = self.patch(pos=[(_[1], _[0]) for _ in pe], n=int(n2/2), value=None, coordinates='centers')
    else:
        im = self.copy()
    return im, pe
br.Image.find_and_patch = find_and_patch
# %%

# %% ============================ centroid =============================== %% #
def centroid(self, n, n2, n3, threshold, threshold2=None, floor=True, MAX_NUMBER_OF_CANDIDATES=300):
    """Returns a list of photon hits defined by a threshold in a enhanced image

    Args:
        n (int): Size of the moving average window used in the process of
            enhancing the image, i.e., number of points to average.
        n2 (int): photon hits candidates that are within n2 pixels of distance 
            from each other will be considered the same candidate. n2 must
            be the pixel distance between a photon hit and the farthest excited
            pixel.
        n3 (int): number of neighbors to include for calculating the center of
            mass, e.g., if n=1, only first neighbors.
        threshold (number): threshold value in terms of pixel value of the 
            ENHANCED image (one might need to plot the enahnced image to get 
            an idea of the pixel intensity).
        threshold2 (number, optional): threshold value for double events in 
            terms of the intensity of the original image (NOT the enhanced 
            image like `threshold`). Any photon hit candidate where the 
            brightest pixel in the original image is above threshold2 is 
            considered as a double event. Also, two photon hit candidates that
            are closer than n3 pixels is considered a double event.
        floor (bool, optional): if True, an intensity offset is added to the image such as the average
            intensity of the whole image is zero. Default is False.
        MAX_NUMBER_OF_CANDIDATES (int, optional): raises error if number of 
            photons to be patched out is larger than MAX_NUMBER_OF_CANDIDATES.
            Useful for preventing too low threshold as patching image is slow.
            Default is 10. 

    Returns:
        pe, pe2
        photon events list, and photon events list for double events
    """
    _pe = self.find_candidates(n, n2, threshold)
    if len(_pe) > MAX_NUMBER_OF_CANDIDATES:
        raise ValueError(f'number of photon hit candidates ({len(_pe)}) is higher than MAX_NUMBER_OF_CANDIDATES={MAX_NUMBER_OF_CANDIDATES}. Maybe try increasing MAX_NUMBER_OF_CANDIDATES, n2, or `threshold`')
    
    # discriminate events that are too close together
    if threshold2 is not None: 
        double = []
        for i, _a in enumerate(_pe):
            for j, _b in enumerate(_pe[i+1:]):
                if math.dist(_a, _b) < n3:
                    double.append(i)
                    double.append(j+1)
        double = br.remove_duplicates(double)

        # rebuild photon events list
        pe2 = br.PhotonEvents()
        pe  = br.PhotonEvents()
        for i, _a in enumerate(_pe):
            if i in double:
                pe2.append(_a[0], _a[1])
            else:
                pe.append(_a[0], _a[1])     
    else:
        pe2 = br.PhotonEvents()
        pe = _pe

    # floor before calculating the center of mass
    if floor:
        im = self.floor()
    else:
        im = self
    
    # get center of mass
    x  = []
    y  = []
    for _pos in pe:
        spot, _, _, _ = im.get_spot(y=_pos[1], x=_pos[0], n=n3, coordinates='centers')
        
        # offset the spot to force all pixel values to be positive
        # this is a requirement for the weighted sum (center of mass) to work fine
        # if one mixes negative and positive weights, than the sum of the weights may
        # be zero (or close to zero) and the weighted average x and/or y positions 
        # are wrong
        # vmin = np.min(spot)
        # vmax = np.max(spot)
        # if vmin < 0 and vmax > 0:
        #     spot = spot - vmin
        spot.data[spot.data < 0] = 0
        # I still don't like this, if we are putting an offset in the spot, I think the 
        # threshold should also be offset. I don't know.
        # but also, It does not make much sense that we have negative pixel values

        _y, _x = spot.get_center_of_mass()
        if threshold2 is not None:
            if spot.max() > threshold2:
                pe2.append(_x, _y)
            else:
                x.append(_x)
                y.append(_y)
        else:
            x.append(_x)
            y.append(_y)

    return br.PhotonEvents(x=x, y=y), pe2
br.Image.centroid = centroid  
# %%

# %% ============================ process ================================ %% #
# %%


# %% -------------------------- OBSOLETE ------------------------ %% #
def calculate_threshold_via_photon_energy(photon_energy, avg_multiplication_factor=0.2, avg_multiplication_factor_double=1):
    """return estimated avg_threshold from photon energy (from ESRF toolbox)
    
    Args:
        photon_energy (number): average energy of the photons hitting the detector
            in eV.
        avg_multiplication_factor (number, optional): analog-to-digital (ADU) 
            multiplication factor. Default is 0.2
        avg_multiplication_factor_double (number, optional): analog-to-digital (ADU) 
            multiplication factor for discriminating double from single events. 
            Default is 1
            
    Returns:
        avg_threshold, double_threshold
    """
    # Multiplication factor * ADU/photon
    photons = photon_energy/3.6/1.06

    avg_threshold    = avg_multiplication_factor * photons
    double_threshold = avg_multiplication_factor_double * photons

    return avg_threshold, double_threshold

def calculate_threshold_via_std(nstd, std_pixel_value, average_pixel_value=0):
    """return estimated avg_threshold via standard deviation method (from XFEL toolbox)

    Args:
        nstd (number): number of standard deviations above the averaged pixel
            value to be considered a photon hit candidate.
        std_pixel_value (number): standard deviation of the pixel value.
        average_pixel_value (number, optional): average pixel value. Default is
            zero.

    Returns:
        avg_threshold, double_threshold
    """
    avg_threshold = average_pixel_value + nstd*std_pixel_value
    double_threshold = 1.5*avg_threshold
    return avg_threshold, double_threshold

def _centroid(self, n1, n2, avg_threshold, double_threshold, floor=False, avg_threshold_max=None, include_doubles=True, MAX_PHOTONS=100000):
    """Returns a PhotonEvents object.
    
    This function is adapted from the ESRF RIXS toolbox and from XFEL SCS toolbox.
    
    Warning:
        for the algorithm to work, the background average must be zero. One can accomplish that by
        using the function im.floor(x_start, x_stop, y_start, y_stop), where start and stop args
        can be used to select a region of the image which is only background. One can also use the
        argument floor=True in this function. This will `floor` the whole image (the average of the
        whole image is zero). This should be ok for most images with low number of photons.
    
    Note:
        image x_centers and y_centers must be monotonic
        
    Args:
        n1 (int): number of pixels to average to detect a photon hit candidate, e.g., a photon hit candidate 
            is selected if the average of the intensities within a n1-by-n1 square exceeds avg_threshold.
        n2 (int): For a photon hit candidates the 'center-of-mass' of the photon hit is calculated 
            within a n2-by-n2 square.
        avg_threshold (number): any pixel with intensity higher than avg_threshold in the n1-by-n1 
            averaged spot surrounding said pixel will be selected as a photon-hit-candidate position.
        double_threshold (number): any a photon-hit-candidate position where the sum of the surrounding
            n2-by-n2 square is higher than double_threshold is considered a double hit.
        floor (bool, optional): if True, an intensity offset is added to the image such as the average
            intensity of the whole image is zero. Default is False.
        avg_threshold_max (number or None, optional): any pixel with intensity higher than avg_threshold_max
            in the n1-by-n1 averaged spot surrounding said pixel will be removed from the photon-hit-candidate list. If None,
            this functionality is disabled and only the lower limit avg_threshold is used. Default is None.
        include_doubles (bool, optional): if False, double photon events will be removed from the
            final result
        MAX_PHOTONS (number or False): if number, this function will raise an error if the number of 
            detected photons exceed MAX_PHOTONS
    
    Returns:
        PhotonEvents and double PhotonEvents in terms of x_centers and y_centers

    Note:
        The averaging n1 prevents that noise is counted as a candidate. By averaging, we are
        requiring that, not only the `central` pixel is illuminated, but also the surrounding 
        pixels. n1=2 seems like a good averaging for most cases.
    
    Note:
        n2 must be roughly the same number of pixel which a photon can excite. For example,
        if n2 = 4, we are expecting that a photon will excite at most a 4x4 array of pixels.
    
    Note:
        if n2 is too large, than the algorithm might mistakenly assign double photon hits as
        it will think that multiple photons will be falling within the n2 squares.
        
    Note:
        I think a good metric for setting up the avg_threshold is to choose the highest threshold possible
        such as one has no photon at `negative` energy loss. One can do this by slowly increasing 
        avg_threshold, and/or by visually inspecting the image using
        
        >>> avg_threshold = <increase_value_slowly>
        >>> n1 = 2
        >>> n2 = 2
        >>> pe, pe2 = scs.centroid(n1, n2, avg_threshold, double_threshold=<very_high_number>, floor=True)
        >>>
        >>> br.figure()
        >>> im.moving_average(n1).set_factor(n1).plot()
        >>> pe.plot()

        if the photon count per image is so low that the elastic line cannot be identified and therefore
        the `negative` energy loss region is unknown, then one have to sum many images.  
        
    Note:
        For determining a reasonable double_threshold, I think the best option is to get a pixel which
        for sure is a single hit, then draw a n2 square around it, and sum the intensities inside 
        this square and multiply by 1.5. For example, given the following matrix with a single photon
        hit in the center (note that the center is slightly shifted to the left and top because the 
        length of the square is a even number):
        
            [[ 4,  6, -4,  8],
             [-9, 79, 52, 10],
             [-2, 19, 17, -2],
             [-6, -9,  2, -8]]
        
        we see that the photon hit is on 79, then for n2=4 we have that the sum is 157 (see also example
        below). A suitable double hit threshold would be 1.5 * 157 = 235.
        
        
    Example:
    
        This example goes over the step-by-step algorithm for detecting a candidate and assigning a photon hit.
    
        Given the following matrix:

            [[14, -7,  4,  6, -4,  8],
             [14,  4, -9, 79, 52, 10],
             [ 1, 17, -2, 19, 17, -2],
             [-1, -4, -6, -9,  2, -8],
             [-4,  3,  2,  7, -9, -3]]

        Its approx. rounded averaged counterpart for n1=2 is:

            [[ 6, -4, 18, 33, 16,  1],
             [ 9,  1, 20, 42, 19,  0],
             [ 3,  1,  0,  6,  2, -5],
             [-2, -1, -2, -5, -7, -6],
             [ 0,  4,  4, -2, -6,  0]]
                
        note that in this example we are using a "approx" rounded average to facilitate the visualization
        of the matrix (instead of an "exact" average), but the script does an exact calculation.
                
        For avg_threshold = 25, we have two spots which are candidates: (x, y) = (3, 0) and (3, 1), with
        intensity values 33 and 42, respectively.
        
        Between these two spots, (3, 0) will be disregarded, because the intensity of (3, 1) is 
        brightest. Going back to the original matrix, we see that position (3, 1) yields intensity 79
        (note that the pixels pixels that yielded 42 in the averaged (n1 x n1) matrix are:
        
        [[79, 52],
         [19, 17]]
        
        From now on, we don't need to worry about the averaged (n1 x n1) matrix. This matrix is only
        used to find the candidates. 

        In this example, the only candidate is (3, 1). We then gather pixel rows and cols to the 
        left/right/top/bottom of the candidate so to make a square of size n2. For n2 = 4 we have
        
        [[ 4,  6, -4,  8],
         [-9, 79, 52, 10],
         [-2, 19, 17, -2],
         [-6, -9,  2, -8]]
        
        this is what we call a `spot`. Now we only have to determine if this spot is a double or
        single hit. 
        
        the sum of the spot is 157. Let's say that double_threshold = 255, therefore, this spot is
        a single hit.
        
        Here is an example where it is known to be a double hit
        
        Original:
        [[-6, -1,  8, -4,  6,  0],
         [ 4, -9, 14, 23,  0,  7],
         [ 4,  4, 99, 41, -9, -9],
         [37, 81, 50, 12, 13, -9],
         [-3, 46, 16, -6,  5, -9]]
        
        approx. rounded averaged for n1=2
        [[-4,  2, 10,  6,  3,  1],
         [ 0, 27, 45, 14, -4, -1],
         [31, 59, 51, 14, -6, -3],
         [40, 48, 18,  6, -1, -1],
         [12, 17,  3, -1,  0, -2]]
            
        From the averaged matrix we have many candidates (intensity > avg_threshold = 25). Taking the
        brightest one, we have x, y = (1, 2), which has intensity 59
        
        the spot for this pixel is (given n2=4)
        
        [[ 4, -9, 14, 23],
         [ 4,  4, 99, 41],
         [37, 81, 50, 12],
         [-3, 46, 16, -6]]
            
        For this example, we know that a photon hit threshold should be around 200. Therefore, a 
        suitable double hit threshold would be 1.5 * 200 = 300.
            
        This sum is 413, so it will be regarded as a double photon hit.
        
        The position of a photon hit is then calculated by the 2D weighted sum of the intensities 
        (the mean x and y values within a spot).
    """
    ###################
    # check n1 and n2 #
    ###################
    assert n1 > 0,  f'n1 must be an positive integer'
    assert n2 > 0,  f'n2 must be an positive even integer'
    # assert n2%2 == 0, f'n2 must be an even number to ensure that the `center of mass` of a photon hit is calculated in a square where the pixel with the photon hit is the central one'
    
    #################################
    # check x_centers and y_centers #
    #################################
    if self.x_monotonicity is None:
        self.check_x_monotonicity()
    if self.y_monotonicity is None:
        self.check_y_monotonicity()
    
    ###############
    # floor image #
    ###############
    if floor:
        image = self.floor()
    else:
        image = self.copy()
    
    ##################
    # moving average #
    ##################
    if n1 > 1:
        image = image.moving_average(n1).set_factor(n1)

    ##########################
    # check double threshold #
    ##########################
    assert double_threshold > avg_threshold, 'avg_threshold must be smaller than double_threshold'
    
    ###################
    # find candidates #
    ###################
    # remove the edges of image (because one cannot calculated the `center of mass` at the edges)
    # select pixels with intensity between low_threshold and high_threshold
    if avg_threshold_max is not None:
        assert avg_threshold_max > avg_threshold, 'avg_threshold must be smaller than avg_threshold_max'
        if n2 > 0:
            cp = np.argwhere((image.data[n2//2:-n2//2, n2//2:-n2//2] > avg_threshold)*(image.data[n2//2 : -n2//2, n2//2 : -n2//2] < avg_threshold_max))
        else:
            cp = np.argwhere((image.data > avg_threshold)*(image.data < avg_threshold_max))
    else:
        if n2 > 0:
            cp = np.argwhere((image.data[n2//2:-n2//2, n2//2:-n2//2] > avg_threshold))
        else:
            cp = np.argwhere((image.data > avg_threshold))
    

    #############################################################
    # shift photon position because we removed the edges before # 
    #############################################################
    # when we were looking for candidates, we excluded the edges
    # now, to have the candidates position in terms of positions 
    # in the original image (image) we have to add back the edges
    # that we excluded
    cp += np.array((n2//2, n2//2))
    
    #################################
    # max allowed number of photons #
    #################################
    if MAX_PHOTONS:
        if len(cp) > MAX_PHOTONS:
            raise RuntimeError(f'number of detected photons is too high (> {MAX_PHOTONS}). Maybe try and change MAX_PHOTONS or Threshold.')

    ########################
    # centroid algorithm 1 #
    ########################
    # runs slower than algorithm 2 (I guess)
    # has the problem that the brightest pixel must have n1 x n1 average above threshold for this pixel to be counted
#         res  = []
#         dres = []
#         for i, (y, x) in enumerate(cp):

#             # isolate the spot where photon hit is in the center
#             spot = self.data[y-n2//2:y+n1+n2//2, x-n2//2:x+n1+n2//2]

#             # check if the central spot is the brightest one
#             if (spot > self.data[y, x]).sum() == 0:

#                 # calculate x center of mass
#                 mx = np.average(np.arange(x-n2//2, x+n1+n2//2), weights=spot.sum(axis=0))

#                 # calculate y center of mass
#                 my = np.average(np.arange(y-n2//2, y+n1+n2//2), weights=spot.sum(axis=1))

#                 # check if spot is a double event
#                 # if (spot.sum() >= SpotLOW) and (spot.sum() < SpotHIGH):         
#                 if (spot.sum() <= double_threshold):               
#                     res.append((mx, my))
#                 else:
#                     if include_doubles:
#                         res.append((mx, my))
#                         res.append((mx, my))
#                     dres.append((mx, my))
        
    ########################
    # centroid algorithm 2 #
    ########################
    # an improved version of algorithm 1
    flag = []
    res  = []
    dres = []
    for i, (y, x) in enumerate(cp):       
        
        # isolate the spot where photon hit is in the center
        if n2%2 == 0:
            spot = self.data[y-n2//2+1:y+n2//2+1, x-(n2//2+1):x+n2//2+1]
        else:
            spot = self.data[y-(n2//2):y+1+n2//2, x-(n2//2):x+1+n2//2]
            
        # offset the spot to force all pixel values to be positive
        # this is a requirement for the weighted sum (center of mass) to work fine
        # if one mixes negative and positive weights, than the sum of the weights may
        # be zero (or close to zero) and the weighted average x and/or y positions 
        # are wrong
        # vmin = np.min(spot)
        # vmax = np.max(spot)
        # if vmin < 0 and vmax > 0:
        #     spot = spot - vmin
        spot[spot < 0] = 0
        # I still don't like this, if we are putting an offset in the spot, I think the 
        # threshold should also be offset. I don't know.
        # but also, It does not make much sense that we have negative pixel values
        
        # get brightest pixel
        _y, _x = np.unravel_index(np.argmax(spot), np.array(spot).shape)
        bx = x + (_x - n2//2)
        by = y + (_y - n2//2)
                            
        if (bx, by) not in flag:
            # flag that this point has been accounted for
            flag.append((bx, by))
            
            # calculate center of mass
            if n2%2 == 0:
                mx = np.average(np.arange(x-n2//2+1, x+n2//2+1), weights=spot.sum(axis=0))
                my = np.average(np.arange(y-n2//2+1, y+n2//2+1), weights=spot.sum(axis=1))
            else:
                mx = np.average(np.arange(x-n2//2, x+1+n2//2), weights=spot.sum(axis=0))
                my = np.average(np.arange(y-n2//2, y+1+n2//2), weights=spot.sum(axis=1))
            # print(x+500, y)
            # print(list(spot))
            # print(mx+500)

            # OBSOLETE CODE FOR AVOIDING SUM OF THE WEIGHTS CLOSE TO ZERO
            # allow_hit_outside_of_spot=False
            # allow_hit_outside_of_spot (bool, optional): default is False. Sometimes, if threshold is too low or
            # noise is too high, it may happen that when calculating the x and y weighted average of the pixel 
            # intensity, the sum of all pixels within a spot [(n1+n2)-by-(n1+n2) square] can be close to 
            # zero. This leads to a wrong weighted sum because the x and/or y averaged position of the hit may
            # fall outside of the spot, which does not make sense. If allow_hit_outside_of_spot = False, it
            # will raise an error when this happens.
            # if allow_hit_outside_of_spot == False:
            #     if mx < x-n2//2 or mx > x+n1+n2//2 or my < y-n2//2 or my > y+n1+n2//2:
            #         report  = 'x values = [' + ','.join([str(_) for _ in np.arange(x-n2//2, x+n1+n2//2)]) + ']\n'
            #         report += 'y values = [' + ','.join([str(_) for _ in np.arange(y-n2//2, y+n1+n2//2)]) + ']\n'
            #         report += 'spot = ' + '[' + ', '.join(['\n        [' + ', '.join([str(_2) for _2 in _]) + ']' for _ in spot])[9:] + ']' + '\n'
            #         report += f'x_weights = {spot.sum(axis=0)}  # spot.sum(axis=0)' + '\n'
            #         report += f'y_weights = {spot.sum(axis=1)}  # spot.sum(axis=1)' + '\n'
            #         report += f'sum_of_weights = {sum(spot.sum(axis=0))}' + '\n'
            #         report += f'x_weighted_avg = {mx}' + '\n'
            #         report += f'y_weighted_avg = {my}' + '\n'
            #         raise ValueError(f'spot `{i}` have center of mass outside of the spot' + '\n' + 'possible causes: 1) threshold too low or 2) noise level too high for single photon count\n' + 'Please, try increasing the threshold' + report)

            # check if spot is a double event
            if (spot.sum() <= double_threshold):               
                res.append((mx, my))
            else:
                if include_doubles:
                    res.append((mx, my))
                    res.append((mx, my))
                dres.append((mx, my))

    ############################
    # convert pixel to centers #
    ############################
    data = np.array([[self.x_centers[int(x)] + (self.x_centers[int(x)+1]-self.x_centers[int(x)])*(x-int(x)), self.y_centers[int(y)] + (self.y_centers[int(y)+1]-self.y_centers[int(y)])*(y-int(y))] for x, y in res])        
    if len(data) == 0:
        pe = br.PhotonEvents()
    else:
        pe = br.PhotonEvents(x=data[:, 0], y=data[:, 1], xlim=(min(self.x_centers), max(self.x_centers)), ylim=(min(self.y_centers), max(self.y_centers)))
    pe.copy_attrs_from(self)
    
    data = np.array([[self.x_centers[int(x)] + (self.x_centers[int(x)+1]-self.x_centers[int(x)])*(x-int(x)), self.y_centers[int(y)] + (self.y_centers[int(y)+1]-self.y_centers[int(y)])*(y-int(y))] for x, y in dres])        
    if len(data) == 0:
        bad = br.PhotonEvents()
    else:
        bad = br.PhotonEvents(x=data[:, 0], y=data[:, 1])
    bad.copy_attrs_from(self)

    # return res, bad
    return pe, bad
# br.Image.centroid = _centroid  