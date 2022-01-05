=======
brixs
=======
Python package for analysis of RIXS spectra.

Documentation: https://cwgaldino.github.io/brixs/


Installation
==================

1) Using pip straight from GitHub

.. code-block:: bash

  pip install git+https://github.com/cwgaldino/brixs

2) Without installing. Clone (or download) the GitHub repository and use

>>> import sys
>>> sys.path.append('<path-to-brixs>')
>>> import brixs as br


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
