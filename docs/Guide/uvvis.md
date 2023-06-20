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

## Performing operations between spectra
The `UVVisSpectrum` class also supports mathematical operations such as the scaling of the spectrum according to a scalar `float` value (`scale`) or the application of binary operations, such as sum (`+`), subtraction (`-`), multiplication (`*`) or division (`/`), between two different `UVVisSpectrum` spectra objects.

:::{admonition} Note
:class: warning
To perfrom binary operations (`+`, `-`, `*`, `/`) between two distinct `UVVisSpectrum` objects, the two objects must be recorded on the same wavelengths scale. If this is not the case the spectra must be either cut [using the `subspectrum()` method](UVVIS-subspectrum) or resampled by interpolation [using the `resample()` method](UVVIS-resample).
:::

The scaling operation can be invoked using the `scale` method according to the command:

```{code-cell} python
scaled = spectrum.scale(0.5)

plot_spectrum([spectrum, scaled], figsize=(10, 6), xrange=(300, 700))
```

Operations between different spectra can be performed by simply using the algebraic symbols (`+`, `-`, `*`, `/`). As an example, the following code will compute the ratio between two spectra:

```{code-cell} python
sp1 = UVVisSpectrum.from_JASCO_ASCII("../utils/example_1.txt")
sp2 = UVVisSpectrum.from_JASCO_ASCII("../utils/example_2.txt")

sp_diff = sp1 / sp2

plot_spectrum([sp1, sp2, sp_diff], figsize=(10, 6), xrange=(300, 700), yrange=(0, 1.2))
```
Please notice how, in both case, the title of the spectrum is automatically updated to keep track of the operation performed. The timestamp is kept unchanged in the case of a simple scaling operation while is set to `None` if a binary operation is considered.

(UVVIS-subspectrum)=
## Extracting a sub-spectrum from a spectrum
Starting from an initialized `UVVisSpectrum` object, a smaller portion of the spectrum, covering only a smaller part of the spectral region can be extracted using the `subspectrum` method. The method can be invoked by simply specifying the range of the new spectrum according to the command:

```{code-cell} python
subspectrum = spectrum.subspectrum(350, 600)

print("Original spectrum:")
print(spectrum)

print("New spectrum:")
print(subspectrum)
```

(UVVIS-interpolate)=
## Interpolating a spectrum
An `UVVisSpectrum` object can be interpolated using a `k`-th order B-spline. The operation will return a `scipy.interpolate.BSpline` object that can than be called to predict absorbance values in a continuous range of values going from the lower to the upper extremes of the wavelength scale of the original object.

```{code-cell} python
bspline = spectrum.interpolate()

print(f"The extimated value of absorbance at 513.4nm is {bspline(513.4):.4f}")
```

(UVVIS-resample)=
## Resample a spectrum
Using the interpolation function, a `UVVisSpectrum` object can be resampled and possibly downsized to match a particular user requirement. A possible application is the use of binary operations: A spectrum recorded with a higher sample rate (smaller pitch) can be "downsampled" to allow operations with other lower resolution spectra. Another application is to approximate the shape of a spectrum starting from a lower resolution one. The resampling operation is conveniently implemented in the `resample` function.

As an example let us show how a low resolution spectra can be resampled using a cubic B-spline:

```{code-cell} python
old = UVVisSpectrum.from_JASCO_ASCII("../utils/iodine_lres.txt")
new = old.resample(lower=350, upper=600, pitch=1)

print(new)

# The following lines are only responsible to plot the results
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(10, 6))
plt.scatter(old.wavelength, old.absorbance, c="black", marker="+", s=100, zorder=3)
plt.plot(new.wavelength, new.absorbance, c="red")
plt.rc("font", **{"size": 18})
plt.xlabel("Wavelength [nm]", size=22)
plt.ylabel("Absorbance [a.u.]", size=22)
plt.grid(which="major", c="#DDDDDD")
plt.grid(which="minor", c="#EEEEEE")
plt.show()
```

## Finding peaks in a spectrum
The `UVVisSpectrum` class provides the built-in `peak_search` method to detect peaks in the spectrum. The function is based around the `scipy.signal.find_peaks` function and operates following a user-defined prominence level. To run a peak search the following command can be used:

```{code-cell} python
spectrum = UVVisSpectrum.from_JASCO_ASCII("../utils/fit_example.txt")
peaks = spectrum.peak_search()
print(peaks)
```

The peaks locaded are returned as a dictionary of tuples encoding wavelength and absorbance associated to an integer key labelling the peak order. The peak finding procedure is also directly implemented in the `plot_spectrum` and can be activated by setting the `peak_prominence` keyword to a float value (we suggest 0.01 absorbance units as a good starting value for UV-Visible spectra). As an example consider the following example:

```{code-cell} python
spectrum = UVVisSpectrum.from_JASCO_ASCII("../utils/fit_example.txt")
plot_spectrum(spectrum, peak_prominence=0.01, xrange=(190, 900), yrange=(-0.1, 4.6))
```

:::{admonition} Tip
:class: info
Axes are not automatically rescaled to fit the peak labels, so you will probably have to manually set the y axis range as shown above with the `yrange` parameter.
:::


# Perfrming a Gaussian fitting

Once a `UVVisSpectrum` object has been created, a Gaussian-based fitting can be performed using the `spectroscopytools.uvvis.FittingEngine`. The fitting routine optimizes a linear combination of gaussian functions combined with a polynobial baseline of different order. The general equation of the composite function $f(x)$ is the following:

$$
  f(x) := \sum_{k=0}^{N_b} c^{(b)}_k x^k + \sum_{k=0}^{N_g} c^{(g)}_k e^{-\frac{(x-x_k)^2}{2\sigma^2}} \qquad \text{with} \qquad \sigma = \frac{w_\mathrm{FWHM}}{2\sqrt{2 \ln(2)}}
$$

where $N_b$ and $N_g$ are set freely by the user. The optimized parameters are the expansion coefficients $c^{(b)}_k$ and $c^{(b)}_k$, the position $x_k$ of each gaussian and its amplitude at half the maximum $w_\mathrm{FWHM}$. A constrained minimization is performed using the `scipy.optimize.curve_fit` function setting the following limits:

* The coefficients $c^{(b)}_k$, $c^{(b)}_k$ and the full width at half maximum $w_\mathrm{FWHM}$ are limited to positive values.
* The center $x_k$ of each gaussian is limited to fall within a finite maximum excursion from the wavelength range included in the spectrum. This can be set by the user (by default $100\mathrm{nm}$ from each end of the spectrum is adopted).

To perform a gaussian fitting using 13 Gaussian function and a zero-order baseline (constant shift), the following code can be used:

```{code-cell} python
from spectroscopytools.uvvis import UVVisSpectrum, FittingEngine

spectrum = UVVisSpectrum.from_JASCO_ASCII("../utils/fit_example.txt")
engine = FittingEngine(spectrum)
engine.fit(13, 0)
engine.plot()
```