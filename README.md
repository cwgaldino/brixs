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

### Spectrum attributes and methods

```python
    # basic
    s.data              # np.array
    s.vmin              # float (read only)
    s.vmax              # float (read only)
    s.shape             # tuple (read only)
    s.x_centers         # np.array
    s.y_centers         # np.array
    s.x_edges           # np.array
    s.y_edges           # np.array

    # binning
    s.nbins             # np.array [runs Image.binning()]
    s.bins_size         # np.array [runs Image.binning()]
    s.reduced           # br.Image (read only)

    # shifts
    s.shifts_v          # np.array
    s.shifts_h          # np.array
    s.p                 # np.array      (read only)
    s.f                 # function f(x) (read only)
    s.calculated_shift  # br.Spectrum   (read only)

    # spectrum (All Computed Attributes)
    s.histogram         # br.Spectrum
    s.spectrum          # br.Spectrum
    s.spectrum_v        # br.Spectrum
    s.spectrum_h        # br.Spectrum
    s.columns           # br.Spectra
    s.rows              # br.Spectra

    # methods
    s.save()
    s.load()

    s.floor()
    s.crop()

    s.pcolormesh()
    s.imshow()
    s.plot()

    s.binning()
    s.calculate_histogram()
    s.calculate_spectrum()
    s.calculate_shift()
    s.set_shift()
    s.fix_curvature()
```
