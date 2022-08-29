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
    im.data              # np.array type
    im.vmin              # float type
    im.vmax              # float type
    im.shape             # tuple type
    im.histogram         # br.Spectrum type
    im.x_centers         # np.array type
    im.y_centers         # np.array type
    im.x_edges           # np.array type
    im.y_edges           # np.array type

    # binning
    im.nbins             # np.array type
    im.bins_size         # np.array type
    im.reduced           # br.Image type

    # shifts
    im.shifts_v          # np.array type
    im.shifts_h          # np.array type
    im.p                 # np.array type
    im.f                 # function f(x)
    im.calculated_shift  # br.Spectrum type

    # spectrum
    im.spectrum          # br.Spectrum type
    im.spectrum_v        # br.Spectrum type
    im.spectrum_h        # br.Spectrum type
    im.columns           # br.Spectra type
    im.rows              # br.Spectra type

    # methods
    im.save()
    im.load()

    im.pcolormesh()
    im.imshow()
    im.plot()

    im.binning()
    im.calculate_histogram()
    im.calculate_spectrum()
    im.floor()
    im.calculate_shift()
    im.set_shift()
    im.fix_curvature()
```
