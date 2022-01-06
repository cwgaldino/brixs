=======
BRIXS
=======
Python package for analysis of RIXS spectra.

Documentation: https://cwgaldino.github.io/brixs/


Installation
==================

There are two recommended methods:

1) Using pip

.. code-block:: bash

  pip install git+https://github.com/cwgaldino/brixs

or

2) Cloning (or downloading) the GitHub repository

.. code-block:: python

  import sys
  sys.path.append('<path-to-brixs>')
  import brixs as br


Requirements
==================

- numpy
- scipy
- h5py
- matplotlib
- beautifulsoup4
- detect_delimiter


Usage
=================

test


.. code-block:: python

  import brixs as br

  # %% photon events with bad events
  pe, s0, nd = br.read_ADRESS(filepath)
  ax = pe.plot()
  bad = br.get_bad_ADRESS(filepath)
  bad.plot(ax, mfc='red', pointsize=2)



.. code-block:: python

  # %% Curvature correction (offset correction)
  pe = br.get_pe_ADRESS(filepath)
  pe.plot()
  plt.ylim(820, 920)
  plt.title('1) Raw photon events')
  pe.bins = (10, 1000)
  pe.plot(show_bins=(True, False))
  plt.ylim(820, 920)
  plt.title('2) X bins')
  pe.calculate_offsets(ref=5)
  pe.fit_offsets(deg=2)
  ax1 = pe.plot_offsets()
  pe.plot_fit(ax1)
  plt.title('3) Offsets and fit')
  pe.plot(show_offsets=True, show_fit=True)
  plt.ylim(820, 920)
  plt.title('4) Offsets and fit over the image')
  pe.offsets_correction()
  pe.plot()
  plt.ylim(820, 920)
  plt.title('5) After correction')

  # %% plotting columns
  pe = br.get_pe_ADRESS(filepath)
  pe.bins = (10, 1000)
  pe.plot_columns()
  plt.title('1) before correction')
  plt.xlim(700, 1000)

  pe.plot_columns(vertical_increment=5)
  plt.title('2) before correction (cascaded)')
  plt.xlim(700, 1000)

  pe.calculate_offsets(ref=5)
  pe.fit_offsets(deg=2)
  pe.offsets_correction()
  pe.plot_columns(vertical_increment=5)
  plt.title('2) after correction')
  plt.xlim(700, 1000)


Read files from ADRESS beamline at PSI

.. code-block:: python

  import brixs as br

  pe, s, nd = br.read_ADRESS(filepath)

  s   = br.get_spectrum_ADRESS(filepath)
  pe  = br.get_pe_ADRESS(filepath)
  bad = br.get_bad_ADRESS(filepath)

  # get data from same scan (three files, one for each ccd d1, d2, d3)
  ss, pes = br.get_scan_ADRESS(filepath, prefix, n, zfill=4)
  ss      = br.get_scan_spectrum_ADRESS(filepath, prefix, n, zfill=4)
  pes     = br.get_scan_pe_ADRESS(filepath, prefix, n, zfill=4)

  # get all scans from energy calibration
  disp, res, _, _ = dispersion_ADRESS(folderpath, prefix, start_energy, stop_energy, start_scan=None, stop_scan=None)


Simulate a spectrum with two excitations with half the area of the
elastic peak and twice its fwhm.

  >>> import brixs
  >>> import matplotlib.pyplot as plt
  >>> import numpy as np
  >>> I = brixs.dummy_spectrum(0, 0.2, excitations=[[0.5, 2, 2], [0.5, 4, 2]])
  >>> x = np.linspace(-2, 6, 1000)
  >>> plt.figure()
  >>> plt.plot(x, I(x))
  >>> plt.xlabel('Energy (eV)')
  >>> plt.ylabel('Intensity')
  >>> plt.show()

  .. image:: _figs/simulate_spectrum.png
      :target: _static/simulate_spectrum.png
      :width: 600
      :align: center


Simulate a ``photon_event`` list of a generic spectrum. The spectrum
will have two excitations with half the area of the
elastic peak and twice its fwhm.

  >>> import brixs
  >>> import matplotlib.pyplot as plt
  >>> import numpy as np
  >>> # simulating a generic spectrum
  >>> I = brixs.dummy_spectrum(0, 0.2, excitations=[[0.5, 2, 2], [0.5, 4, 2]])
  >>> # simulating the photon_event list(where we're using energy in eV's and length in meters)
  >>> photon_events = brixs.dummy_photon_events(I, background=0.02,
  >>>                                                 noise=0.05,
  >>>                                                 exposure=50e4,
  >>>                                                 dispersion= 8.45 * (10**-3 / 10**-6),
  >>>                                                 x_max=52.22e-3,
  >>>                                                 y_max=25.73e-3,
  >>>                                                 y_zero_energy=-20,
  >>>                                                 angle=0,
  >>>                                                 psf_fwhm=(3e-6, 1e-6))
  >>> print(photon_events)
      [[1.36263387e-02 2.29071963e-03 1.00000000e+00]
       [3.19917559e-02 4.48965702e-04 1.00000000e+00]
       [4.96047073e-02 9.76363776e-05 1.00000000e+00]
       ...
       [2.37889174e-04 1.50658922e-02 1.00000000e+00]
       [1.58122734e-02 4.33273762e-03 1.00000000e+00]
       [4.61469585e-02 2.45566831e-03 1.00000000e+00]]
  >>> # ploting photon_events
  >>> plt.figure()
  >>> plt.plot(photon_events[:, 0]*10**3,
  >>>          photon_events[:, 1]*10**3,
  >>>          linewidth=0,
  >>>          marker='o',
  >>>          ms=1)
  >>> plt.xlabel('x position (mm)')
  >>> plt.ylabel('y position (mm)')
  >>> plt.show()

  .. image:: _figs/photon_events1.png
      :target: _static/photon_events1.png
      :width: 600
      :align: center

  Zooming in, we can clearly see the isoenergetic lines formed by the elastic
  peak and the other two excitations.

  .. image:: _figs/photon_events2.png
      :target: _static/photon_events2.png
      :width: 600
      :align: center
