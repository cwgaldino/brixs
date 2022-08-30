# BRIXS

Python package for analysis of XAS/RIXS spectra.

Documentation: https://cwgaldino.github.io/brixs/


## Installation

There are two recommended methods:

1. Using pip

```python
    pip install git+https://github.com/cwgaldino/brixs
```

or

2. Cloning (or downloading) the GitHub repository

```python
    import sys
    sys.path.append('<path-to-brixs>')
    import brixs as br
```

## Requirements

Base (required):
- numpy
- scipy
- matplotlib

Reading files:
- detect_delimiter
- h5py

Reciprocal space calculations:
- pbcpy

## Usage (not an extensive description)

Refer to the examples folder on GitHub for code examples.

BRIXS is based on four objects:

```python
    import brixs as br

    im = br.Image()
    pe = br.PhotonEvents()
    s  = br.Spectrum()
    ss = br.Spectra()
```

### Image attributes and methods

```python
    # basic
    im.data              # np.array
    im.vmin              # float (read only)
    im.vmax              # float (read only)
    im.shape             # tuple (read only)
    im.x_centers         # np.array
    im.y_centers         # np.array
    im.x_edges           # np.array
    im.y_edges           # np.array

    # binning
    im.nbins             # np.array [runs Image.binning()]
    im.bins_size         # np.array [runs Image.binning()]
    im.reduced           # br.Image (read only)

    # shifts
    im.shifts_v          # np.array
    im.shifts_h          # np.array
    im.p                 # np.array      (read only)
    im.f                 # function f(x) (read only)
    im.calculated_shift  # br.Spectrum   (read only)

    # spectrum (All Computed Attributes)
    im.histogram         # br.Spectrum
    im.spectrum          # br.Spectrum
    im.spectrum_v        # br.Spectrum
    im.spectrum_h        # br.Spectrum
    im.columns           # br.Spectra
    im.rows              # br.Spectra

    # methods
    im.save()
    im.load()

    im.floor()
    im.crop()
    
    im.pcolormesh()
    im.imshow()
    im.plot()

    im.binning()
    im.calculate_histogram()
    im.calculate_spectrum()
    im.calculate_shift()
    im.set_shift()
    im.fix_curvature()
```
