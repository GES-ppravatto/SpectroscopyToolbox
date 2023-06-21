(changelog)=
# Release notes

## Version `0.0.2-alpha`

* Added new features to the `UVVisSpectrum` class:
  * Added `wavenumber` and `electrovolt` properties to rescale plots according to wavenumber and energies.
  * Added `peak_search` function.
* Added new options to the `plot_spectrum` function:
  * Added option to select different units for the x-axis.
  * Added option to integrate peak search labels.
* Early implementation of the `FittingEngine` class to perform gaussian fitting with baseline of a `UVVisSpectrum` object.
* Update to documentation and tests (to be improved in the beta version).

## Version `0.0.1`

* Early implementation of the `uvvis` submodule:
  * Defined a `UVVisSpectrum` class to store and manipulate UV-Visible specta:
    * Defined a classmethod to load data contained in ASCII file formatted by JASCO spectrophotometers.
    * Defined basic operations between spectra and a scale function.
    * Defined basic derived properties like wavelength, absorbance, transmittance and timesptamp.
    * Defined basic interpolation functions based on B-spline (scipy)
  * Defined a simple `plot_spectrum` function to plot the spectroscopic data.
  * Early version of the documentation
  * Early version of unit tests