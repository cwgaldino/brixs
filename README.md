# BRIXS

Python package for analysis of XAS/RIXS spectra.

Documentation: https://cwgaldino.github.io/brixs/

WARNING: The code has been in active development this past few days and some documentation might not be up to date.

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

    # spectrum (Computed Attributes)
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
    s.x                 # np.array
    s.y                 # np.array
    s.area              # float (computed attribute)

    # modifiers
    s.offset        # float
    s.factor        # float
    s.calib         # float
    s.shift         # float
    s.shift_roll    # int
    s.shift_interp  # float

    # check
    s.step          # float
    s.monotonicity  # string

    # fit
    s.fit     # br.Spectrum
    s.residue # br.Spectrum
    s.guess   # br.Spectrum
    s.R2      # float

    # peaks
    s.peaks   # br.peaks

    # methods
    s.save()
    s.load()

    s.plot()

    s.check_step_x()
    s.check_monotonicity()
    s.fix_monotonicity()

    s.set_calib()
    s.set_shift()
    s.set_factor()
    s.set_offset()

    s.interp()
    s.extract()
    s.crop()
    s.floor()
    s.flip()
    s.normalize()
    s.zero()
    s.calculate_area()

    s.find_peaks()
    s.fit_peak()
    s.fit_peaks()
    s.polyfit()
    s.apply_correction()
```


## Information for developers

I will list here some practices that I have been trying to follow while developing brixs. This should facilitate the maintenance and expansion of the package.

- the metaclass `_Meta()` helps other classes creating read-only and non-removable attributes. The variables `_read_only = []` and `_non_removable = []` should be defined at the beginning of a class before __init__(). Attributes will be created based on the names listed inside these lists.
- __init__() methods should initialize *ALL* attributes, except for computed attributes.
- __init__() should call _sort_args().
- _sort_args() should sort the arguments and return a data type or a filepath to __init__(), and nothing else
- f
