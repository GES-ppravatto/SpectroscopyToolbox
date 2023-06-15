from os.path import abspath, dirname
from datetime import datetime
from numpy.testing import assert_array_almost_equal, assert_almost_equal

from spectroscopytools.uvvis import UVVisSpectrum

# Get the path of the tests directory
TEST_DIR = dirname(abspath(__file__))


# Test UVVisSpectrum class constructor
def test_UVVisSpectrum___init__():
    try:
        UVVisSpectrum()
    except:
        assert False, "Exception occurred on UVVisSpectrum class construction"


# Test UVVisSpectrum from_JASCO_ASCII class metod
def test_UVVisSpectrum_from_JASCO_ASCII():
    try:
        obj = UVVisSpectrum.from_JASCO_ASCII(f"{TEST_DIR}/utils/JASCO_ASCII_example.txt")
    except:
        assert False, "Exception raised on loading JASCO_ASCII"

    assert obj.title == "I2_water"
    assert len(obj) == 21
    assert obj.timestamp == datetime(2023, 6, 14, 15, 14, 38)


# Test UVVisSpectrum __iter__ method
def test_UVVisSpectrum___iter__():

    obj = UVVisSpectrum.from_JASCO_ASCII(f"{TEST_DIR}/utils/JASCO_ASCII_example.txt")

    expected_results = [
        (700.0, 0.00953),
        (680.0, 0.01085),
        (660.0, 0.01255),
        (640.0, 0.01447),
        (620.0, 0.01692),
        (600.0, 0.02068),
        (580.0, 0.02688),
        (560.0, 0.03929),
        (540.0, 0.06284),
        (520.0, 0.10024),
        (500.0, 0.14998),
        (480.0, 0.19426),
        (460.0, 0.21183),
        (440.0, 0.18741),
        (420.0, 0.13148),
        (400.0, 0.07338),
        (380.0, 0.04436),
        (360.0, 0.04166),
        (340.0, 0.04171),
        (320.0, 0.03982),
        (300.0, 0.06099),
    ]

    for idx, (wl, ab) in enumerate(obj):
        assert_almost_equal(expected_results[idx], (wl, ab), decimal=4)


# Test UVVisSpectrum __getitem__ method
def test_UVVisSpectrum___getitem__():

    obj = UVVisSpectrum.from_JASCO_ASCII(f"{TEST_DIR}/utils/JASCO_ASCII_example.txt")

    expected_results = [
        (700.0, 0.00953),
        (680.0, 0.01085),
        (660.0, 0.01255),
        (640.0, 0.01447),
        (620.0, 0.01692),
        (600.0, 0.02068),
        (580.0, 0.02688),
        (560.0, 0.03929),
        (540.0, 0.06284),
        (520.0, 0.10024),
        (500.0, 0.14998),
        (480.0, 0.19426),
        (460.0, 0.21183),
        (440.0, 0.18741),
        (420.0, 0.13148),
        (400.0, 0.07338),
        (380.0, 0.04436),
        (360.0, 0.04166),
        (340.0, 0.04171),
        (320.0, 0.03982),
        (300.0, 0.06099),
    ]
    
    for idx in range(len(obj)):
        assert_almost_equal(expected_results[idx], obj[idx], decimal=4)


# Test UVVisSpectrum properties 
def test_UVVisSpectrum_properties():

    obj = UVVisSpectrum.from_JASCO_ASCII(f"{TEST_DIR}/utils/JASCO_ASCII_example.txt")

    expected_results = [
        (700.0, 0.00953),
        (680.0, 0.01085),
        (660.0, 0.01255),
        (640.0, 0.01447),
        (620.0, 0.01692),
        (600.0, 0.02068),
        (580.0, 0.02688),
        (560.0, 0.03929),
        (540.0, 0.06284),
        (520.0, 0.10024),
        (500.0, 0.14998),
        (480.0, 0.19426),
        (460.0, 0.21183),
        (440.0, 0.18741),
        (420.0, 0.13148),
        (400.0, 0.07338),
        (380.0, 0.04436),
        (360.0, 0.04166),
        (340.0, 0.04171),
        (320.0, 0.03982),
        (300.0, 0.06099),
    ]

    assert_array_almost_equal([d[0] for d in expected_results], obj.wavelength, decimal=4)
    assert_array_almost_equal([d[1] for d in expected_results], obj.absorbance, decimal=4)
    assert_array_almost_equal([10**(2-d[1]) for d in expected_results], obj.transmittance, decimal=4)
