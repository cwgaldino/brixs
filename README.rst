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


Usage (not an extensive description)
====================================

Package is based on three objects

.. code-block:: python

  import brixs as br

  pe = br.PhotonEvents()
  s  = br.Spectrum()
  ss = br.Spectra()

Spectra/Spectrum manipulation.

.. code-block:: python

  import brixs as br
  import numpy as np

  # ==============================================================================
  # %% spectra alignment
  filelist = br.backpack.filelist('../fixtures/calib_example')
  ss = br.Spectra()
  for f in filelist:
      ss.append(data=np.loadtxt(f, delimiter=','))

  plt.figure()
  ss.plot(vertical_increment=0.1)
  ss.calculate_shifts(ref=0, mode='cc')
  ss.set_shifts()
  plt.figure()
  ss.plot(vertical_increment=0.1)
  ss.crop(stop=600)
  plt.figure()
  ss.plot(vertical_increment=0.1)

  # ==============================================================================
  # %% find and fit peaks (easy)
  filepath = '../fixtures/peak_fit/easy.dat'
  s = br.Spectrum(data=np.loadtxt(filepath, delimiter=','))
  s.find_peaks()
  plt.figure()
  s.plot()
  s.plot_detected_peaks()

  s.fit_peak(0)
  plt.figure()
  s.plot()
  s.fit.plot()

  s.fit_peak()
  plt.figure()
  s.plot()
  s.fit.plot()

  # results
  s.fit_data[0]
  print('elastic peak at', s.fit_data[0]['c'])
  print('first excitation at', s.fit_data[1]['c'])


  # ==============================================================================
  # %% find and fit peaks (medium)
  # shoulder is not detected
  filepath = '../fixtures/peak_fit/medium.dat'
  s = br.Spectrum(data=np.loadtxt(filepath, delimiter=','))
  s.find_peaks()
  plt.figure()
  s.plot()
  s.plot_detected_peaks()

  # fit is wrong
  s.fit_peak(0)
  plt.figure()
  s.plot()
  s.fit.plot()

  # fit is good
  s.fit_peak(0, multiplicity=2)
  plt.figure()
  s.plot()
  s.fit.plot()

  # results
  print('elastic peak at', s.fit_data[0]['c'])

  # ==============================================================================
  # %% find and fit peaks (hard)
  filepath = '../fixtures/peak_fit/hard.dat'
  s = br.Spectrum(data=np.loadtxt(filepath, delimiter=','))
  s.find_peaks()
  plt.figure()
  s.plot()
  s.plot_detected_peaks()

  # fit is good
  s.fit_peak(multiplicity={0:2})
  plt.figure()
  s.plot()
  s.fit.plot()

  # results
  print('elastic peak at', s.fit_data[0]['c'])

  # ==============================================================================
  # %% energy calibration (answer is 10 meV/point)
  filelist = br.backpack.filelist('../fixtures/calib_example')
  ss = br.Spectra()
  for f in filelist:
      ss.append(data=np.loadtxt(f, delimiter=','))

  calib = ss.calculate_calib(start=930, stop=939)
  print(calib*1000, ' meV/point')

  plt.figure()
  ss.plot_disp()

  # ==============================================================================
  # %% energy calibration using peak fit
  filelist = br.backpack.filelist('../fixtures/calib_example')
  ss = br.Spectra()
  for f in filelist:
      ss.append(data=np.loadtxt(f, delimiter=','))

  ss.find_peaks(prominence=0.5)
  calib = ss.calculate_calib(start=930, stop=939, mode='peak', idx=0)
  print(calib*1000, ' meV/point')

  plt.figure()
  ss.plot(vertical_increment=0.1)
  for i, s in enumerate(ss):
      s.fit.plot(offset=-0.1*i, color='black')

  plt.figure()
  ss.plot_disp()


Creating fake data for testing purposes

.. code-block:: python

  # creating a sequence of fake spectra
  positions = np.linspace(100, 1000, 10)
  energies = np.linspace(930, 939, 10)
  ss = br.Spectra()
  for c in positions:
      fake = br.fake(amp=1, c=c, fwhm=10)
      s = fake.get_spectrum(0, 1500, n_points=6000, noise=3)
      ss.append(s)

  plt.figure()
  ss.plot(vertical_increment=0.1)

  ss.save(r'../fixtures/calib_example')

  # ==============================================================================
  # create spectrum with excitations (multiple peaks)
  fake1 = br.fake(1, 0, 0.1, [[0.5, 3, 0.2], [0.7, 4, 0.2]])
  s1 = fake1.get_spectrum(noise=3)
  # big elastic peak with low energy excitation
  fake2 = br.fake(1, 0, 0.1, [[0.5, 0.15, 0.2], [0.7, 4, 0.2]])
  s2 = fake2.get_spectrum(noise=3)
  # small elastic peak with intense low energy excitation
  fake3 = br.fake(0.5, 0, 0.1, [[1, 0.15, 0.2], [0.7, 4, 0.2]])
  s3 = fake3.get_spectrum(noise=3)

  plt.figure()
  s1.plot(label='easy')
  s2.plot(offset=0.5, label='medium')
  s3.plot(offset=1, label='hard')
  plt.legend()

  s1.save(r'../fixtures/peak_fit/easy.dat')
  s2.save(r'../fixtures/peak_fit/medium.dat')
  s3.save(r'../fixtures/peak_fit/hard.dat')
