#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module for generating fake RIXS spectra.
"""

# .. autosummary::
#
#     fake_spectrum
#     fake_photon_events

# standard libraries
import numpy as np
import copy

# specific libraries
from collections.abc import Iterable

# backpack
from .backpack.model_functions import gaussian_fwhm
from .backpack.arraymanip import index

# brixs
import brixs as br


class _Meta(type):
    """Metaclass to facilitate creation of read-only attributes."""
    def __new__(self, class_name, bases, attrs):

        def lazy_read_only(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    return getattr(self, variable)
                def setter(self, value):
                    raise AttributeError('Attribute is "read only". Cannot set attribute.')
                def deleter(self):
                    raise AttributeError('Attribute is "read only". Cannot delete object.')
            return getter, setter, deleter, 'read only variable'

        def lazy_non_removable(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    return getattr(self, variable)
                def setter(self, value):
                    return setattr(self, variable, value)
                def deleter(self):
                    raise AttributeError('Attribute cannot be deleted.')
            return getter, setter, deleter, 'non removable variable'

        new_attrs = {}
        for name, value in attrs.items():
            if name == '_read_only':
                for attr in value:
                    _property = property(*lazy_read_only(attr))
                    new_attrs[attr] = _property
            elif name == '_non_removable':
                for attr in value:
                    _property = property(*lazy_non_removable(attr))
                    new_attrs[attr] = _property
            else:
                new_attrs[name] = value

        return type(class_name, bases, new_attrs)

class fake(metaclass=_Meta):

    _read_only     = ['func', 'pe', 'spectrum']

    def __init__(self, amp=1, c=0, fwhm=1, excitations=None):
        self.elastic_amp = amp
        self.elastic_c = c
        self.elastic_fwhm = fwhm

        self._func = None
        self._pe = None
        self._spectrum = None

        if excitations is not None:
            self.excitations = excitations
        else:
            self.excitations = []

    @property
    def elastic_amp(self):
        return self._elastic_amp
    @elastic_amp.setter
    def elastic_amp(self, value):
        if value < 0:
            raise ValueError('amplitude (amp) cannot be negative.')
        else:
            self._elastic_amp = value
    @elastic_amp.deleter
    def elastic_amp(self):
        raise AttributeError('Cannot delete object.')

    @property
    def elastic_c(self):
        return self._elastic_c
    @elastic_c.setter
    def elastic_c(self, value):
        self._elastic_c = value
    @elastic_c.deleter
    def elastic_c(self):
        raise AttributeError('Cannot delete object.')

    @property
    def elastic_fwhm(self):
        return self._elastic_fwhm
    @elastic_fwhm.setter
    def elastic_fwhm(self, value):
        if value < 0:
            raise ValueError('width (fwhm) cannot be negative.')
        else:
            self._elastic_fwhm = value
    @elastic_fwhm.deleter
    def elastic_fwhm(self):
        raise AttributeError('Cannot delete object.')

    @property
    def excitations(self):
        return self._excitations
    @excitations.setter
    def excitations(self, value):
        if isinstance(value, Iterable):
            for v in value:
                if isinstance(v, Iterable):  # multiple excitations
                    if len(v) == 3:
                        if v[0] < 0:
                            raise ValueError(f'amp cannot be negative ({value[0]}).')
                        if v[2] < 0:
                            raise ValueError(f'width cannot be negative ({value[2]}).')
                    else:
                        raise ValueError('excitation must be a list (or a list of lists) with amp, c, and width')
                else:
                    if len(value) == 3:
                        if value[0] < 0:
                            raise ValueError(f'amp cannot be negative ({value[0]}).')
                        if value[2] < 0:
                            raise ValueError(f'width cannot be negative ({value[2]}).')
                        self._excitations = [value, ]
                        return
                    else:
                        raise ValueError('excitation must be a list (or a list of lists) with amp, c, and width')
            self._excitations = value
        else:
            raise ValueError('excitation must be a list (or a list of lists) with amp, c, and width')
    @excitations.deleter
    def excitations(self):
        raise AttributeError('Cannot delete object.')

    def get_function(self):
        """Returns a function I(E) of a simulated RIXS spectrum.

        I(E) is a function that returns the intensity of the simulated spectrum. The
        spectrum is composed by an elastic peak and many other peaks refered to as
        excitations. All peaks have a gaussian profile.

        Args:
            c (number, optional): position of elastic peak in energy. Default is 0.
            w (number, optional): fwhm of the elastic peak in energy. Default is 1
            excitations (list, optional): list of excitations. Each element must be a list
                with three elements: the relative area compared to the elastic peak
                area, the distance from the elastic peak in energy units, and the
                relative width compared to the elastic peak width.

        Returns:
            function I(E).

        .. seealso:: :py:func:`get_photon_events`
        """
        I = f'lambda energy: gaussian_fwhm(energy, {self.elastic_amp}, {self.elastic_c}, {self.elastic_fwhm})'

        if self.excitations is not None:
            for excitation in self.excitations:
                if excitation[0] < 0:
                    raise ValueError('amplitude (amp) of excitation cannot be negative.')
                if excitation[2] < 0:
                    raise ValueError('width (fwhm) of excitation cannot be negative.')
                I += f'+gaussian_fwhm(energy, amp={excitation[0]}, c={excitation[1]}, w={excitation[2]})'
        self._func = eval(I)
        return self.func

    def append(self, *args):
        if len(args) == 1:
            if isinstance(args[0], Iterable):
                if len(args[0]) == 3:
                    if value[0] < 0:
                        raise ValueError(f'amp cannot be negative ({args[0][0]}).')
                    if value[2] < 0:
                        raise ValueError(f'width cannot be negative ({args[0][2]}).')
                    self._excitations.append([args[0][0], args[0][1], args[0][2]])
            else:
                raise ValueError('excitation must be a list with amp, c, and width')
        elif len(args) == 3:
            if args[0] < 0:
                raise ValueError(f'amp cannot be negative ({args[0]}).')
            if args[2] < 0:
                raise ValueError(f'width cannot be negative ({args[2]}).')
            self._excitations.append([args[0], args[1], args[2]])
        else:
            raise ValueError('excitation must be a list (or a list of lists) with amp, c, and width')

        self.excitations.append([amp, c, w])

    def remove(self, idx):
        del self.excitation[idx]

    def _get_x_min(self):
        if self.excitations is None:
            return self.elastic_c-self.elastic_fwhm*10
        else:
            energy_min1 = self.elastic_c-self.elastic_fwhm*10
            energy_min2 = min([excitation[1] - excitation[2]*10 for excitation in self.excitations])
            if energy_min2 > energy_min1:
                return energy_min1
            else:
                return energy_min2

    def _get_x_max(self):
        if self.excitations is None:
            return self.elastic_c+self.elastic_fwhm*10
        else:
            energy_max1 = self.elastic_c+self.elastic_fwhm*10
            energy_max2 = max([excitation[1] + excitation[2]*10 for excitation in self.excitations])
            if energy_max2 < energy_max1:
                return energy_max1
            else:
                return energy_max2

    def _get_y_max(self):
        if self.excitations is None:
            return self.elastic_amp
        else:
            max_value1 = self.elastic_amp
            if self.excitations != []:
                max_value2 = max([excitation[0] for excitation in self.excitations])
                if max_value1 > max_value2:
                    return max_value1
                else:
                    return max_value2
            else:
                return max_value1


    def get_spectrum(self, xmin=None, xmax=None, n_points=None, step=None, x=None, noise=0):
        """returns x, y.

            Args:
                n_points (int, optional): number of

            Returns:
                :py:class:`brixs.Spectrum` object.

            Raises:
                ValueError: noise is negative.
        """
        if x is None:
            if xmin is None:
                xmin = self._get_x_min()
            if xmax is None:
                xmax = self._get_x_max()
            if step is None and n_points is None:
                n_points = round((energy_max-energy_min)/self.elastic_fwhm*10)
            elif step is not None:
                if xmax-xmin < step:
                    raise ValueError('step is too big for data')
                x = np.arange(xmin, xmax, step)
            elif n_points is not None:
                if n_points < 1:
                    raise ValueError('n_points must be higher than 1.')
                n_points = int(round(n_points))
                x = np.linspace(xmin, xmax, n_points)
        else:
            x = np.array(x)


        if noise > 0:
            noise = self._get_y_max()*noise/100
            random_noise = np.random.default_rng().normal(-noise, noise, size=int(n_points))
        elif noise < 0:
            raise ValueError('noise is a percentage value and cannot be negative.')
        else:
            random_noise = 0

        if self.func is None:
            _ = self.get_function()

        y = abs(self.func(x) + random_noise)
        self._spectrum = br.Spectrum(x=x, y=y)
        return self._spectrum

    def get_photon_events(self, dispersion, # (energy)/(detector units)
                                x_max=None,  # in detector units (detector size, number of pixels, or number of bins)
                                y_max=None,  # in detector units (detector size, number of pixels, or number of bins)
                                noise=0,       # in percentage of the maximum intensity
                                exposure=100e4,
                                elastic_position=None,  # if None, it will be half of y_max
                                elastic_side = 'right',
                                poly_coef=[0],  # curvature coeficients in decreasing power
                                # psf_fwhm=(0, 0)
                                ):

        right = ['right', 'r', 'RIGHT', 'R', 'max', 'top', 't', 'above']
        left = ['left', 'l', 'LEFT', 'L', 'min', 'bottom', 'below']

        if self.excitations is None:
            max_value = self.elastic_amp
        else:
            max_value = max([excitation[0] for excitation in self.excitations])

        if y_max is None:
            if self.excitations is None:
                energy_max = self.elastic_c+self.elastic_fwhm*10
            else:
                if elastic_side in left:
                    energy_max = max([excitation[1] + excitation[2]*10 for excitation in self.excitations])
                    # energy_min = min([excitation[1] - excitation[2]*10 for excitation in self.excitations])
                    y_max = abs(energy_max/dispersion)
                else:
                    # energy_max = max([excitation[1] + excitation[2]*10 for excitation in self.excitations])
                    energy_min = min([excitation[1] - excitation[2]*10 for excitation in self.excitations])
                    y_max = abs(energy_min/dispersion)

        if x_max is None:
            x_max = y_max

        if elastic_position is None:
            elastic_position = y_max/2

        noise      = max_value*noise/100

        # create random (x, y) positions
        random_y     = (y_max*np.random.rand(int(exposure)))
        random_x     = (x_max*np.random.rand(int(exposure)))
        random_noise = np.random.default_rng().normal(-noise, noise, size=int(exposure))

        # create probability distribution
        if self.func is None:
            _ = self.get_function()
        # energy = (random_y + (random_x-x_max/2) * (np.tan(np.radians(angle))) )
        energy = (random_y  - elastic_position)*dispersion
        prob = self.func(energy) + random_noise
        prob = prob/max(prob)

        # test each position against the probability distribution
        choose   = (prob > np.random.rand(len(prob)))
        y_spatial = random_y[choose]  # m
        x_spatial = random_x[choose]  # m

        if elastic_side in left:
            y_spatial = 2*elastic_position - y_spatial

        # apply curvature
        f_poly = np.poly1d(poly_coef)
        f = lambda x, y: (x, y+f_poly(x))
        x_spatial, y_spatial = f(x_spatial, y_spatial)

        # add intensity of 1 for each position
        I2 = np.ones(len(y_spatial))

        # # apply PSF
        # try:
        #     if len(psf_fwhm) == 1:
        #         x_spread, y_spread = psf_fwhm[0]/dispersion/2.355, psf_fwhm[0]/2.355  # 2.355 is to convert from fwhm to sigma
        #     else:
        #         x_spread, y_spread = psf_fwhm[0]/dispersion/2.355, psf_fwhm[1]/2.355
        # except TypeError:
        #     x_spread, y_spread = psf_fwhm/dispersion/2.355, psf_fwhm/dispersion/2.355
        #
        # if x_spread > 0:
        #     x_spatial = [x+np.random.normal(-x_spread/2, x_spread/2) for x in x_spatial]
        #
        #     # fix position that fall outside the detector
        #     while max(x_spatial) > x_max:
        #         x_spatial[index(x_spatial, max(x_spatial))] = x_max
        # if y_spread > 0:
        #     y_spatial = [y+np.random.normal(-y_spread/2, y_spread/2) for y in y_spatial]
        #
        #     # fix position that fall outside the detector
        #     while max(y_spatial) > y_max:
        #         y_spatial[index(y_spatial, max(y_spatial))] = y_max

        photon_events = np.vstack((x_spatial, y_spatial, I2)).transpose()
        self.pe = br.PhotonEvents(photon_events, x_max=x_max, y_max=y_max)
        return self.pe


# def fake_spectrum(amp=1, c=0, w=1, excitations=None):
#
#     if amp < 0:
#         raise ValueError('amplitude (amp) cannot be negative.')
#     if w < 0:
#         raise ValueError('width (fwhm) cannot be negative.')
#
#     I = f'lambda energy: gaussian_fwhm(energy, {amp}, {c}, {w})'
#
#     if excitations is not None:
#         for excitation in excitations:
#             if excitation[0] < 0:
#                 raise ValueError('amplitude (amp) of excitation cannot be negative.')
#             if excitation[2] < 0:
#                 raise ValueError('width (fwhm) of excitation cannot be negative.')
#             I += f'+gaussian_fwhm(energy, amp={excitation[0]}, c={excitation[1]}, w={excitation[2]})'
#     return eval(I)

# def fake_photon_events(I, background=0,
#                               noise=0,
#                               exposure=100e4,
#                               dispersion= 8.45 * (10**-3 / 10**-6),
#                               x_max=52.22e-3,
#                               y_max=25.73e-3,
#                               elastic_position=None,
#                               angle=0,
#                               psf_fwhm=(0, 0)):
#     r"""Returns a list of simulated photon events, based on a simulated spectrum.
#
#     Args:
#         I (function): a function that returns the intensity of a spectrum as a
#             function of the energy. The maximum of the elastic peak is expected
#             to be equal to 1.
#         background (float, optional): this value will be added to
#             the spectrum. Can also be called offset. The final spectrum is
#             normalized so that the maximum of the elastic peak is equal to 1.
#             Recomended value: from 0 to 1, assuming the  maximum of the elastic
#             peak is equal to 1.
#         noise (float, optional): This
#             adds a random gaussian probability (from -``noise`` to ``noise``) to the
#             probability finding a photon at a certain position. Recomended value:
#             from 0 to 1, assuming the  maximum of the elastic peak is equal to 1.
#         exposure (float, optional): number of random (x, y) positions to test test
#             the probability of finding a photon.
#         dispersion (float, optional): must be given in units of [energy/distance\. This is
#             the value used for converting from energy units to positional units.
#         x_max (float, optional): maximum value of x (detector's size in the x
#             direction).
#         y_max (float, optional): maximum value of y (detector's size in the y
#             direction).
#         y_zero_energy (float, optional): energy value of isoenergetic lines at
#             y=0 (minimum energy value). The maximum energy value is calculated
#             by ``y_max`` :math:`\times` ``dispersion``.
#         angle (float, optional): rotation of the detector.
#         psf_fwhm (float or list, optional): fwhm value of the gaussian error
#             included in the (x, y) position to simulated precision errors of the detector and beamline.
#             If one value is given, this is used for both x and y directions. If two values are given,
#             they are used separetely for the x and y directions, respectively.
#
#     .. warning::
#
#         The script does not expect the user to input
#         the parameter in any specific set units, like energy in eV's or lenght in
#         meters. The only requirement is that the user keeps the units consistent.
#         If meter is used as the unit of lenght of a certain parameter, then meter
#         shall be used in all other parameters.
#
#     Returns:
#         array with three columns, x, y, and I. The latter will always be equal to 1, i.e., only one
#         photon is detected at each position.
#
#     .. seealso:: :py:func:`fake_spectrum`
#
#     This functions creates a 2d probability distribution of detecting a photon
#     onto the detector's active area. Then, it creates a list of random (x, y)
#     positions whithin (0, ``x_max``) and (0, ``y_max``). Finally, it tests each
#     positions against the probability distribution to see if a photon is
#     detected on that position.
#
#     The figure below shows a schematic representation
#     of a RIXS spectrometer. In (a), we can see that the diffraction gratting
#     will disperse scattered photon onto the detector along the y direction. The
#     dispersion given in units of energy per length must be known in advance and
#     assigned to ``dispersion`` parameter. By default, the position ``y=0`` is
#     assigned to energy zero, but one can change that by assinging a different value
#     to ``y_zero_energy``.
#
#     In (a) one can see the that the ``angle`` parameter rotates the detector and
#     desalign the isoenergetic lines with the edge of the detector (See (b) and (c)).
#
#     .. image:: _figs/spectrometer_1.svg
#         :target: _static/spectrometer_1.svg
#         :width: 600
#         :align: center
#
#     The algorithm works in the following manner:
#
#     1) we create a list of random (x, y) positions. The variable ``exposure``
#     gives the lenght of this list. The higher the value of ``exposure``, the more
#     random position will be created and tested. Note that, this parameter emulates
#     the effect of a longer or shorter exposure time, however, it cannot be directly
#     associated with the real exposure time of a experiment.
#
#     2) The variable ``I(E)``, which is a function that returns the intensity of a spectrum as a
#     function of the energy (This can be generated by :py:func:`fake_spectrum`), is
#     proportional to the probability of finding a photon at a certain position in the detector.
#     Therefore, one can calculate the probability of detecting a photon at each of the
#     random (x, y) position created in step 1. The equation below depicts how background and
#     noise affect this probability.
#
#     .. math::
#
#         \text{Prob}(x, y) = (I(\text{energy}  + \text{y_zero_energy}) + \text{background})*(1 + \text{random_noise}) / (1 + \text{background})
#
#     where the (x, y) position is converted into energy by taking into account the
#     ``dispersion`` and the parameter ``angle``.
#
#     3) After we have calculated the probability of finding an photon for each (x, y)
#     position we test each probability, e.g., if the probability of finding a photon
#     in a certain position is 70%, we simulate a riggeed coin flip where the
#     outcome has 70% chance of being heads. After tossing the coin, if the result
#     is heads, the photon is counted and the respetive (x, y) position is included
#     in the ``photon_event`` list. If the result is tails, that position is
#     discarted.
#
#     4) At this point the ``photon_event`` list is ready. As a final step,
#     we can simulate precision errors which tipically comes from the
#     limited resolution of the detector or possible instabilities of the beamline
#     itself. This is done by adding a random gaussian error in each of the (x, y)
#     positions in the ``photon_event`` list, e.g., a certain position :math:`(x, y)` where
#     a photon as detected will became :math:`(x+\Delta x, y+\Delta y)`, where
#     :math:`\Delta x` and :math:`\Delta y` are a random gaussian error.
#     The width of the gaussian error is set by `psf_fwhm`.
#
#     """
#     if elastic_position is None:
#         elastic_position = y_max/2
#
#     # create random (x, y) positions
#     random_y     = (y_max*np.random.rand(int(exposure)))
#     random_x     = (x_max*np.random.rand(int(exposure)))
#     random_noise = np.random.default_rng().normal(-noise, noise, size=int(exposure))
#
#     # create probability distribution
#     energy = (random_y + (random_x-x_max/2) * (np.tan(np.radians(angle))) )
#     print(min(energy))
#     print(max(energy))
#     prob = (I((energy  - elastic_position)*dispersion) + background)*(1 + random_noise) / (1+background)
#     print(max(prob))
#     # test each position against the probability distribution
#     choose   = (prob > np.random.rand(len(prob))*max(prob))
#     y_spatial = random_y[choose]  # m
#     x_spatial = random_x[choose]  # m
#
#     # add intensity of 1 for each position
#     I2 = np.ones(len(y_spatial))
#
#     # apply PSF
#     try:
#         if len(psf_fwhm) == 1:
#             x_spread, y_spread = psf_fwhm[0]/2.355, psf_fwhm[0]/2.355
#         else:
#             x_spread, y_spread = psf_fwhm[0]/2.355, psf_fwhm[1]/2.355
#     except TypeError:
#         x_spread, y_spread = psf_fwhm/2.355, psf_fwhm/2.355
#
#     if x_spread > 0:
#         x_spatial = [x+np.random.normal(0, x_spread) for x in x_spatial]
#
#         # fix position that fall outside the detector
#         while max(x_spatial) > x_max:
#             x_spatial[index(x_spatial, max(x_spatial))] = x_max
#     if y_spread > 0:
#         y_spatial = [y+np.random.normal(0, y_spread) for y in y_spatial]
#
#         # fix position that fall outside the detector
#         while max(y_spatial) > y_max:
#             y_spatial[index(y_spatial, max(y_spatial))] = y_max
#
#     return np.vstack((x_spatial, y_spatial, I2)).transpose()
