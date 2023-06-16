import setuptools

setuptools.setup(
    name="spectroscopy-toolbox",
    version="0.0.1",
    description="A small collection of python-based tools useful in the analysis of experimental spectroscopic data",
    long_description="A small collection of python-based tools useful in the analysis of experimental spectroscopic data",
    packages=["spectroscopytools"],
    package_data={
        "spectroscopytools": ["*"],
    },
    install_requires=[],
)