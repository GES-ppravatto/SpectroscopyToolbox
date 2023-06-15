---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(Guide-UVVIS)=
# Using the UV-Visible spectroscopy submodule

The `spectroscopytools.uvvis` submodule provides a simple interface to the analysis of UV-Visible spectra. The core of the module is represented by the `UVVisSpectrum` class. This object can be constructed directly from an experimental file according to a specific classmethod. At the moment the following format are supported:

* ASCII format exported by JASCO instruments (classmethod: `from_JASCO_ASCII`)

As an example, an ASCII file produced by a JASCO V-550 instrument can be loaded using the command:

```{code-cell} python
from spectroscopytools.uvvis import UVVisSpectrum

spectrum = UVVisSpectrum.from_JASCO_ASCII("../utils/iodine.txt")

print(spectrum)
```

The recorded spectrum can be plotted using the `plot_spectrum` function according to:

```{code-cell} python
from spectroscopytools.uvvis import plot_spectrum

plot_spectrum(spectrum, figsize=(10, 6), xrange=(300, 700))
```