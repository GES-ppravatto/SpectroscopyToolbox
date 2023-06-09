from os.path import abspath, dirname
from datetime import datetime

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
    
    assert len(obj) == 711
    assert obj.timestamp == datetime(2023, 6, 7, 12, 22, 15)