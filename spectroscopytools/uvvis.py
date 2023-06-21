from __future__ import annotations

from os.path import isfile
from typing import Tuple, List, Optional, Union, Dict
from datetime import datetime
from warnings import warn
from scipy.interpolate import make_interp_spline, BSpline
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from spectroscopytools.constants import *

import numpy as np
import matplotlib.pyplot as plt


def unitary_gaussian(x: float, x0: float, fwhm: float) -> float:
    r"""
    Unitary height Gaussian function :math:`G(x)`.

    .. math:: G(x):= e^{-\frac{(x-x_0)^2}{2\sigma^2}} \qquad \mathrm{where} \qquad \sigma = \frac{w_\mathrm{FWHM}}{2\sqrt{2\ln(2)}}

    Arguments
    ---------
    x: float
        The float value encoding the current value of the coordinate.
    x0: float
        The float value encoding the position of the center of the Gaussian function.
    fwhm: float
        The float value indicating the full width at half maximum setting the amplitude of the Gaussian function.

    Returns
    -------
    float
        The value of the Gaussian function at the point `x`.
    """
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-0.5 * ((x - x0) / sigma) ** 2)


class UVVisSpectrum:
    """
    The UVVisSpectrum class allows the manipulation of UV-Visible data originated from experimental
    measurements. The class provides a standard constructor returning an empty object and a series
    of classmethods designed to parse specific file formats.
    """

    def __init__(self) -> None:
        self.title: str = None
        self.instrument: Optional[str] = None
        self.__timestamp: Optional[datetime] = None
        self.__wavelength: List[float] = []
        self.__absorbance: List[float] = []

    def __str__(self) -> str:
        msg = f"UV-Visible: {self.title}\n"
        msg += "-----------------------------------------------------------------\n"
        msg += f"Date: {self.__timestamp}\n"
        msg += f"Instrument: {self.instrument}\n"
        msg += f"Wavelength: {max(self.__wavelength):.1f} - {min(self.__wavelength):.1f} nm\n"
        msg += f"Max absorbance: {max(self.__absorbance):.4f}\n"
        msg += f"Min absorbance: {min(self.__absorbance):.4f}\n"
        return msg

    def __repr__(self) -> str:
        return str(self)

    @property
    def timestamp(self) -> datetime:
        """
        The time at which the measurement has been started.

        Returns
        -------
        datetime
            The datetime object encoding the time at which the measurement has been started.
        """
        return self.__timestamp

    @property
    def wavelength(self) -> List[float]:
        """
        The wavelength values in nanometers associated with each datapoint.

        Returns
        -------
        List[float]
            The list of float values encoding the wavelength, in nanometers, associated with each datapoint.
        """
        return self.__wavelength

    @property
    def electronvolt(self) -> List[float]:
        """
        The energy values in electronvolts (eV) associated with each datapoint.

        Returns
        -------
        List[float]
            The list of float values encoding the energy, in eV, associated with each datapoint.
        """
        energy = [J_to_eV * PLANCK_CONSTANT * SPEED_OF_LIGHT / (1e-9 * wl) for wl in self.wavelength]
        return energy

    @property
    def wavenumber(self) -> List[float]:
        """
        The wavenumber in cm^-1 associated with each datapoint.

        Returns
        -------
        List[float]
            The list of float values encoding the wavenumber, in cm^-1, associated with each datapoint.
        """
        wavenumber = [1 / (1e-7 * wl) for wl in self.wavelength]
        return wavenumber

    @property
    def absorbance(self) -> List[float]:
        """
        The absorbance values associated with each datapoint.

        Returns
        -------
        List[float]
            The list of float values encoding absorbance associated with each datapoint.
        """
        return self.__absorbance

    @property
    def transmittance(self) -> List[float]:
        """
        The transmittance, expressed as a percentage value, associated with each datapoint.

        Returns
        -------
        List[float]
            The list of float values encoding transmittance associated with each datapoint.
        """
        return [10 ** (2 - A) for A in self.__absorbance]

    @property
    def pitch(self) -> float:
        """
        The wavelength separation existing between two subsequent data-points.

        Returns
        -------
        float
            The wavelength separation in nanometers existing between to subsequent datapoints. If the points are not
            equally spaced a warning will be raised and the average pitch value will be returned
        """
        delta = [self.wavelength[i] - self.wavelength[i - 1] for i in range(1, len(self))]

        if any([value != delta[0] for value in delta]):
            warn("The spectral points are not equally spaced. The pitch property will return the average pitch value.")
            return abs(sum(delta) / len(delta))

        else:
            return abs(delta[0])

    @classmethod
    def from_JASCO_ASCII(cls, path: str) -> UVVisSpectrum:
        """
        The classmethod designed to parse ASCII data files generated by JASCO instruments.

        Arguments
        ---------
        path: str
            The path to the ASCII file encoding the experimental measurements.

        Raises
        ------
        RuntimeError
            Exception raised if the file does not exist or if it cannot be properly parsed.
        """

        if not isfile(path):
            raise RuntimeError(f"The specified file '{path}' does not exist.")

        obj = cls()

        with open(path, "r") as file:
            npt, xunits, yunits = None, None, None

            for line in file:
                if "TITLE" in line:
                    obj.title = line.split("\t")[-1].strip("\n")

                if "SPECTROMETER" in line:
                    obj.instrument = line.split("\t")[-1].strip("\n")

                if "XUNITS" in line:
                    xunits = line.split("\t")[-1].strip("\n")

                if "YUNITS" in line:
                    yunits = line.split("\t")[-1].strip("\n")

                if "DATE" in line:
                    date_str = line.split("\t")[-1].strip("\n")
                    date_str += " "
                    date_str += file.readline().split("\t")[-1].strip("\n")
                    obj.__timestamp = datetime.strptime(date_str, "%y/%m/%d %H:%M:%S")

                if "NPOINTS" in line:
                    npt = int(line.split("\t")[-1])

                if "XYDATA" in line:
                    for _ in range(npt):
                        data = file.readline().split("\t")

                        xvalue, yvalue = float(data[0]), float(data[1])

                        if xunits == "NANOMETERS":
                            obj.__wavelength.append(xvalue)
                        else:
                            raise RuntimeError(f"Cannot parse unit {xunits}.")

                        if yunits == "ABSORBANCE":
                            obj.__absorbance.append(yvalue)
                        else:
                            raise RuntimeError(f"Cannot parse unit {yunits}.")

        obj.__sort()
        return obj

    def __getitem__(self, i: int) -> Tuple[float, float]:
        if i < 0 or i >= len(self):
            raise ValueError("Index out of bounds")

        return self.__wavelength[i], self.__absorbance[i]

    def __iter__(self) -> Tuple[float, float]:
        for w, a in zip(self.__wavelength, self.__absorbance):
            yield w, a

    def __len__(self) -> int:
        return len(self.__wavelength)

    def __sort(self) -> None:
        """
        Function sorting the spectal datapoint from the smaller wavelength to the highest.
        """
        self.__absorbance = [A for _, A in sorted(zip(self.__wavelength, self.__absorbance), key=lambda pair: pair[0])]
        self.__wavelength.sort()

    def __check_binary_operation(self, obj: UVVisSpectrum) -> None:
        """
        The function checks whether a given UVVisSpectrum object can be used in a binary operation involving
        this UVVisSpectum.

        Raises
        ------
        RuntimeError
            Exception raised if the lenght of the spectrum objects, their pitch or their spectral range is different.
        """
        if len(self) != len(obj):
            raise RuntimeError("Cannot perform binary operation between spectra of different lengths.")

        if self.pitch != obj.pitch:
            raise RuntimeError("Cannot perform binary operation between spectra of different pitch.")

        if not all([w1 == w2 for w1, w2 in zip(self.wavelength, obj.wavelength)]):
            raise RuntimeError("Cannot perform binary operation between spectra with different wavelength ranges")

    def __add__(self, obj: UVVisSpectrum) -> UVVisSpectrum:
        self.__check_binary_operation(obj)

        result = UVVisSpectrum()
        result.title = self.title + " + " + obj.title
        result.__timestamp = None
        result.__wavelength = self.__wavelength
        result.__absorbance = [A1 + A2 for A1, A2 in zip(self.__absorbance, obj.__absorbance)]

        return result

    def __sub__(self, obj: UVVisSpectrum) -> UVVisSpectrum:
        self.__check_binary_operation(obj)

        result = UVVisSpectrum()
        result.title = self.title + " - " + obj.title
        result.__timestamp = None
        result.__wavelength = self.__wavelength
        result.__absorbance = [A1 - A2 for A1, A2 in zip(self.__absorbance, obj.__absorbance)]

        return result

    def __mul__(self, obj: UVVisSpectrum) -> UVVisSpectrum:
        self.__check_binary_operation(obj)

        result = UVVisSpectrum()
        result.title = self.title + " * " + obj.title
        result.__timestamp = None
        result.__wavelength = self.__wavelength
        result.__absorbance = [A1 * A2 for A1, A2 in zip(self.__absorbance, obj.__absorbance)]

        return result

    def __truediv__(self, obj: UVVisSpectrum) -> UVVisSpectrum:
        self.__check_binary_operation(obj)

        result = UVVisSpectrum()
        result.title = self.title + " / " + obj.title
        result.__timestamp = None
        result.__wavelength = self.__wavelength
        result.__absorbance = [A1 / A2 for A1, A2 in zip(self.__absorbance, obj.__absorbance)]

        return result

    def scale(self, value: float, inplace: bool = False) -> Optional[UVVisSpectrum]:
        """
        Scale the absorbance of the spectrum according to a float scalar value.

        Arguments
        ---------
        value: float
            The factor to be used in the multiplication
        inplace: bool
            If set to False (default) will return the spectrum scaled by the specified value. Else, it will update the
            absorbance list of the current spectrum with the updated values.

        Returns
        -------
        Optional[UVVisSpectrum]
            If inplace is set to False (default) will return the scaled spectrum.
        """

        if inplace is True:
            self.title = f"{value}*{self.title}"
            self.__absorbance = [value * A for A in self.__absorbance]

        else:
            result = UVVisSpectrum()
            result.title = f"{value}*{self.title}"
            result.__timestamp = self.__timestamp
            result.__wavelength = self.__wavelength
            result.__absorbance = [value * A for A in self.__absorbance]

            return result

    def subspectrum(self, lower: float, upper: float) -> UVVisSpectrum:
        """
        Generate from the spectrum a new UVVisSpectrum object containing only the datapoints between the user specified
        lower and upper limits.

        Arguments
        ---------
        lower: float
            The lower limit of the new spectrum (included if available)
        upper: float
            The upper limit of the new spectrum (included if available)

        Raises
        ------
        ValueError
            Exception raised if the lower limit is not smaller than the upper one.

        Returns
        -------
        UVVisSpectrum
            The spectrum object containing the experimental point withing the user specified wavelength bounaries.
        """
        if lower >= upper:
            raise ValueError("The lower limit must be smaller than the upper one.")

        subspectrum = UVVisSpectrum()
        subspectrum.title = self.title
        subspectrum.instrument = self.instrument
        subspectrum.__timestamp = self.__timestamp

        for wavelength, absorbance in self:
            if wavelength >= lower and wavelength <= upper:
                subspectrum.__wavelength.append(wavelength)
                subspectrum.__absorbance.append(absorbance)

        return subspectrum

    def interpolate(self, k: int = 3) -> BSpline:
        """
        Compute a k-th order interpolating B-spline for the spectrum.

        Arguments
        ---------
        k: int
            The B-spline degree (default: cubic B-spline k=3)

        Returns
        -------
        BSpline
            The interpolating B-spline function as a `scipy.interpolate.BSpline` object.
        """
        bspline = make_interp_spline(self.__wavelength, self.__absorbance, k=k)
        return bspline

    def resample(self, lower: float, upper: float, pitch: float, k: int = 3) -> UVVisSpectrum:
        """
        Using a B-spline interpolation, generate a new UVVisSpectum object responding the the user specified requirements.

        Arguments
        ---------
        lower: float
            The minimum value of the wavelength scale in nanometers.
        upper: float
            The maximum value of the wavelength scale in nanometers.
        pitch: float
            The value of the wavelength sparation, in nanometers, existing between two subsequent datapoints.
        k: int
            The B-spline degree used during the resample (default: cubic B-spline k=3)

        Raises
        ------
        ValueError
            Exception raised if the `lower` is bigger or equal to the `upper` value or if `pitch` is invalid.
        RuntimeError
            Exception raised if the `lower` and `upper` values specified by the user are outside the original spectum range.

        Returns
        -------
        UVVisSpectrum
            The UVVisSpectrum generated according to the user specified set of parameters.
        """
        if lower >= upper:
            raise ValueError("The lower wavelength value must be lower than the upper one.")

        if pitch <= 0:
            raise ValueError("The pitch value must be a positive float.")

        if lower < min(self.wavelength) or upper > max(self.wavelength):
            raise RuntimeError("Cannot perform interpolation. The required range exceed the range of available data.")

        spectrum = UVVisSpectrum()
        spectrum.title = self.title + " (resampled)"
        spectrum.instrument = f"B-spline (k={k}) interpolation <- original: {self.instrument}"
        spectrum.__timestamp = self.__timestamp

        spectrum.__wavelength = np.arange(start=lower, stop=upper + pitch, step=pitch)

        bspline = self.interpolate(k=k)
        spectrum.__absorbance = bspline(spectrum.__wavelength)
        spectrum.__sort()

        return spectrum

    def peak_search(self, prominence: float = 0.01) -> Dict[int, Tuple[float, float]]:
        """
        Function running the `scipy.signal.find_peaks` funciton to detect peaks in the spectrum.

        Arguments
        ---------
        prominence: float
            The prominence threshold to be used in the peak search (default: 0.01). The prominence of a peak measures
            how much a peak stands out from the surrounding baseline of the signal and is defined as the vertical
            distance between the peak and its lowest contour line.

        Returns
        -------
        Dict[int, Tuple[float, float]]
            The dictionary containing the index of the peak as the key and the tuple containing the wavelength and
            absorbance as the value.
        """
        idx_list = find_peaks(self.__absorbance, prominence=prominence)[0]

        peak_dict = {}
        for i, idx in enumerate(idx_list):
            peak_dict[i] = (self.__wavelength[idx], self.__absorbance[idx])

        return peak_dict


def plot_spectrum(
    spectra: Union[List[UVVisSpectrum], UVVisSpectrum],
    xunit: str = "wavelength",
    transmittance: bool = False,
    xrange: Optional[Tuple[float, float]] = None,
    yrange: Optional[Tuple[float, float]] = None,
    figsize: Tuple[float, float] = (12.0, 8.0),
    savepath: Optional[str] = None,
    peak_prominence: Optional[float] = None,
    show: bool = True,
):
    """
    Function capable of plotting one or more spectra.

    Arguments
    ---------
    spectra: Union[List[UVVisSpectrum], UVVisSpectrum]
        The spectrum or the list of spectra to plot.
    xunit: str
        The unit to be adopted for the x-axis. Available options: "wavelength", "wavenumber", "electronvolt"
    transmittance: bool
        If set to True wil switch the y-axis to transmittance mode.
    xrange: Optional[Tuple[float, float]]
        The range of values to be shown on the x-axis. The meaning of the values depends on the display mode selected.
    yrange: Optional[Tuple[float, float]]
        The range of values to be shown on the y-axis. The meaning of the values depends on the display mode selected.
    figsize: Tuple[float, float]
        The size of the matplotlib figure to be used in plotting the spectrum. (default: (12, 8))
    savepath: Optional[str]
        If set to a value different from None, will specify the path of the file to be saved.
    peak_prominence: Optional[float]
        If set to a value different from None, will run a peak_search on the spectrum with the user provided prominence
        and will reder a marker on the spectrum.
    show: bool
        If set to True (default) will show an interactive window where the plot is displayed.

    Raises
    ------
    ValueError
        Exception raised if the selected xunit option is invalid:
    """
    plt.rc("font", **{"size": 18})

    fig = plt.figure(figsize=figsize)

    spectra: List[UVVisSpectrum] = [spectra] if type(spectra) == UVVisSpectrum else spectra

    for spectrum in spectra:
        xseries = None
        if xunit == "wavelength":
            xseries = spectrum.wavelength
        elif xunit == "wavenumber":
            xseries = spectrum.wavenumber
        elif xunit == "electronvolt":
            xseries = spectrum.electronvolt
        else:
            raise ValueError(f"The xunit option `{xunit}` is invalid.")

        line = plt.plot(
            xseries,
            spectrum.transmittance if transmittance else spectrum.absorbance,
            label=f"{spectrum.title}",
        )

        if peak_prominence is not None:
            color = line[-1].get_color()
            peak_report = spectrum.peak_search(prominence=peak_prominence)

            for key, (x, y) in peak_report.items():
                plt.scatter(x, y + 0.08, c=color, marker="v", s=90)
                plt.text(
                    x,
                    y + 0.18,
                    s=f"{x:.1f}",
                    rotation=90,
                    horizontalalignment="center",
                    verticalalignment="bottom",
                    size=16,
                    c=color,
                )

    if xrange is not None:
        plt.xlim(xrange)

    if yrange is not None:
        plt.ylim(yrange)

    if xunit == "wavelength":
        plt.xlabel("Wavelength [nm]", size=22)
    elif xunit == "wavenumber":
        plt.xlabel(r"Wavenumber [$cm^{-1}$]", size=22)
    elif xunit == "electronvolt":
        plt.xlabel(r"Energy [eV]", size=22)

    plt.ylabel(r"Transmittance [$\%$]" if transmittance else "Absorbance [a.u.]", size=22)

    plt.grid(which="major", c="#DDDDDD")
    plt.grid(which="minor", c="#EEEEEE")

    plt.legend()

    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=600)

    if show is True:
        plt.show()


class FittingEngine:
    """
    Simple class capable of performing Gaussian fitting of UVVisSpectrum objects. The class allows the user to fit the
    spectrum using a linear combination of Gaussians superimposed to a user defined polynomial baseline.

    Arguments
    ---------
    spectrum: UVVisSpectrum
        The spectrum to be fitted
    """

    def __init__(self, spectrum: UVVisSpectrum) -> None:
        self.__ngaussians: int = 0
        self.__baseline_degree: int = None
        self.__spectrum: UVVisSpectrum = spectrum
        self.__popt: List[float] = []
        self.__pcov: List[List[float]] = []
        self.__infodict: dict = (None,)
        self.__mesg: str = (None,)
        self.__ier: int = None

    def __str__(self) -> str:
        msg = ""

        if self.__mesg != "":
            msg += "============================================================================\n"
            msg += "                            FITTING REPORT\n"
            msg += "============================================================================\n"
            msg += "Output message: " + self.__mesg + "\n\n"

        if self.__baseline_degree is not None:
            msg += "Polynomial baseline:\n"
            msg += "  -----------------------\n"
            msg += "  | order | coefficient |\n"
            msg += "  -----------------------\n"
            for n in range(self.__baseline_degree + 1):
                p_n = self.optimized_parameters[n]
                msg += "  | {0:<5} | {1:^11} |\n".format(f"{n:d}", f"{p_n:.4e}")
            msg += "  -----------------------\n"
            msg += "\n"

        msg += "Gaussian functions:\n"
        msg += "  ---------------------------------------------------\n"
        msg += "  | index | coefficient | center (nm) |  FWHM (nm)  |\n"
        msg += "  ---------------------------------------------------\n"

        offset = self.__baseline_degree + 1 if self.__baseline_degree is not None else 0
        for n in range(self.__ngaussians):
            c_n = self.optimized_parameters[3 * n + offset]
            x0_n = self.optimized_parameters[3 * n + offset + 1]
            fwhm_n = self.optimized_parameters[3 * n + offset + 2]
            msg += "  | {0:<5} | {1:^11} | {2:^11} | {3:^11} |\n".format(
                f"{n:d}", f"{c_n:.4e}", f"{x0_n:.2f}", f"{fwhm_n:.4e}"
            )
        msg += "  ---------------------------------------------------\n"

        return msg

    def __repr__(self) -> str:
        return str(self)

    @property
    def optimized_parameters(self) -> List[float]:
        """
        The list of optimized parameters minimizing the squared residuals of the composite function difference to the
        provided spectroscopic data.

        Returns
        -------
        List[float]
            The list of optimized parameters.
        """
        if len(self.__popt) == 0:
            raise RuntimeError(
                "Cannot obtain the optimized parameters. The fit function has failed or it has not been called."
            )
        return self.__popt

    @property
    def estimated_covariance(self) -> List[List[float]]:
        """
        The estimated covariance of the optimized_parameters. The diagonals provide the variance of the parameter estimate.

        Returns
        -------
        List[List[float]]
            The covariance matrix of the optimized_parameters.
        """
        if len(self.__pcov) == 0:
            raise RuntimeError(
                "Cannot obtain the etimated covariance. The fit function has failed or it has not been called."
            )
        return self.__pcov

    def combined_function(self, x: float, *p: Tuple[float]) -> float:
        """
        Combined fitting function incorporating the combination of a polynomial baseline with a linear combination of
        Gaussian functions.

        Arguments
        ---------
        x: float
            The float value associated to the coordinata at which the function must be computed.
        *p: Tuple[float]
            The tuple encoding the float parameters setting the baseline and gaussian shapes. The ordering of the
            parameters is the following: the first `n+1` parameters, where `n` represents the degree of the baseline
            polynomials, sets the coefficients of the various powers of the coordinate. The others represents the
            parameters associated to the gaussians. These are ordered in triplets of: coefficient, center and FWHM.

        Returns
        -------
        float
            The value of combined function at the slected coordinate value.
        """
        sum = 0
        if self.__baseline_degree is not None:
            for n in range(self.__baseline_degree + 1):
                sum += p[n] * (x**n)

        offset = self.__baseline_degree + 1 if self.__baseline_degree is not None else 0

        for n in range(self.__ngaussians):
            sum += p[3 * n + offset] * unitary_gaussian(x, p[3 * n + 1 + offset], p[3 * n + 2 + offset])

        return sum

    def fit(
        self,
        ngaussians: int,
        baseline_degree: Optional[int] = None,
        max_feval: int = 10000000,
        x0_bound: float = 100,
        verbose: bool = True,
    ) -> None:
        """
        Function running the fitting process.

        Arguments
        ---------
        ngaussians: int
            The number of gaussian functions to be used in the optimization.
        baseline: Optional[int]
            The degree of the baseline polynomial.
        max_feval: int
            The maximum number of function evaluations allowed during the optimization (default: 1e7).
        x0_bound: float
            The maximum extension below the minimum and above the maximum of the center of a gaussian function used in
            fitting.
        verbose: bool
            If verbose is set to True (default), A report will be printed at the end of the fitting process.

        Raises
        ------
        ValueError
            Exception raised if the number of gaussians is less than one or if an invalid value is given to the baseline
            degree.
        """
        if ngaussians <= 0:
            raise ValueError("The number of gaussian used in the fitting cannot be smaller than 1.")

        if baseline_degree is not None:
            if baseline_degree < 0:
                raise ValueError("The baseline degree must be a non-negative integer.")

        self.__ngaussians = ngaussians
        self.__baseline_degree = baseline_degree

        wlmin = min(self.__spectrum.wavelength)
        wlmax = max(self.__spectrum.wavelength)

        bounds = ([], [])
        guess = []
        if baseline_degree is not None:
            for n in range(baseline_degree + 1):
                bounds[0].append(0)
                bounds[1].append(np.inf)
                guess.append(0)

        wldelta = (wlmax - wlmin) / self.__ngaussians
        for n in range(ngaussians):
            guess.append(max(self.__spectrum.absorbance) / 2)
            bounds[0].append(0)
            bounds[1].append(np.inf)

            guess.append(0.5 * wldelta + n * wldelta + wlmin)
            bounds[0].append(wlmin - x0_bound)
            bounds[1].append(wlmax + x0_bound)

            guess.append(20)
            bounds[0].append(0)
            bounds[1].append(np.inf)

        self.__popt, self.__pcov, self.__infodict, self.__mesg, self.__ier = curve_fit(
            self.combined_function,
            self.__spectrum.wavelength,
            self.__spectrum.absorbance,
            p0=guess,
            maxfev=max_feval,
            bounds=bounds,
            full_output=True,
        )

        if verbose:
            print(self)

    def plot(
        self,
        figsize: Tuple[float, float] = (16.0, 10.0),
        savepath: Optional[str] = None,
        show: bool = True,
    ) -> None:
        """
        Function to plot the results of the fitting operation.

        Arguments
        ---------
        figsize: Tuple[float, float]
            The size of the matplotlib figure to be used in plotting the spectrum. (default: (16, 10))
        savepath: Optional[str]
            If set to a value different from None, will specify the path of the file to be saved.
        show: bool
            If set to True (default) will show an interactive window where the plot is displayed.

        Raises
        ------
        RuntimeError
            Exception raised if the `fit` function has not been called before.
        """
        if len(self.__popt) == 0:
            raise RuntimeError("Cannot plot the results, the optimized parameters list is empty.")

        plt.rc("font", **{"size": 18})

        fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=figsize, gridspec_kw={"height_ratios": [4, 1]})

        p = self.__popt
        wl = self.__spectrum.wavelength

        ax1.plot(wl, self.__spectrum.absorbance, c="red", linewidth=1.5)

        if self.__baseline_degree is not None:
            baseline = []
            for x in wl:
                sum = 0
                for n in range(self.__baseline_degree + 1):
                    sum += p[n] * (x**n)
                baseline.append(sum)

            ax1.plot(wl, baseline, c="#888888", linewidth=0.5)
            ax1.fill_between(wl, [0 for _ in wl], baseline, alpha=0.5, label="baseline")

        offset = self.__baseline_degree + 1 if self.__baseline_degree is not None else 0
        for n in range(self.__ngaussians):
            gaussian = []
            for x in wl:
                y = p[3 * n + offset] * unitary_gaussian(x, p[3 * n + 1 + offset], p[3 * n + 2 + offset])
                gaussian.append(y)

            ax1.plot(wl, gaussian, c="#888888", linewidth=0.5)
            ax1.fill_between(wl, [0 for _ in wl], gaussian, alpha=0.5, label=str(n))

        fit = [self.combined_function(x, *p) for x in wl]
        error = [100 * (x - y) / x for x, y in zip(self.__spectrum.absorbance, fit)]

        ax1.plot(
            wl,
            fit,
            c="black",
            linestyle="--",
            linewidth=1.5,
            label="fit",
        )

        ax2.plot(wl, error, c="red", linewidth=1.5)

        ax1.set_xlim((min(wl), max(wl)))
        ax2.set_xlim((min(wl), max(wl)))

        ax2.set_xlabel("Wavelength [nm]", size=22)

        ax1.set_ylabel("Absorbance [a.u.]", size=22)
        ax2.set_ylabel("Error [%]", size=22)

        ax1.legend()

        plt.tight_layout()

        if savepath is not None:
            plt.savefig(savepath, dpi=600)

        if show is True:
            plt.show()
