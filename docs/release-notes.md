(changelog)=
# Release notes

## Version `0.0.1`

* Early implementation of the `uvvis` submodule:
  * Defined a `UVVisSpectrum` class to store and manipulate UV-Visible specta:
    * Defined a classmethod to load data contained in ASCII file formatted by JASCO spectrophotometers.
    * Defined basic operations between spectra and a scale function.
    * Defined basic derived properties like wavelength, absorbance, transmittance and timesptamp.
  * Defined a simple `plot_spectrum` function to plot the spectroscopic data.
  