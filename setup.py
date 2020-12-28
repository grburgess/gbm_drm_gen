import io
import os

import numpy
from numpy.distutils.core import Extension  # as Numpy_Extension
from numpy.distutils.core import setup
from setuptools import Command, find_packages  # Extension
from setuptools.command.build_ext import build_ext as _build_ext

import versioneer

# Package meta-data.
NAME = "gbm_drm_gen"
DESCRIPTION = "A DRM generator for the Fermi GBM detectors"
URL = "https://github.com/mpe-grb/gbm_drm_gen"
EMAIL = "jmichaelburgess@gmail.com"
AUTHOR = "J. Michael Burgess"
REQUIRES_PYTHON = ">=2.7.0"
VERSION = None

REQUIRED = [
    'numpy',
    'astropy',
    'scipy',
    'h5py',
    'numba',
    #    'threeML',
    'gbmgeometry'
]

# What packages are optional?
EXTRAS = {
    # 'fancy feature': ['django'],
}

here = os.path.abspath(os.path.dirname(__file__))

try:
    with io.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


# Create list of data files
def find_data_files(directory):

    paths = []

    for (path, directories, filenames) in os.walk(directory):

        for filename in filenames:

            paths.append(os.path.join('..', path, filename))

    print(paths)

    return paths


extra_files = find_data_files("gbm_drm_gen/data")

# setup(

#     ext_modules=cythonize(extensions_cython), requires=['numpy']
# )

setup(
    name=NAME,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=("tests",)),
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    package_data={
        '': extra_files,
    },
    license='BSD',

    #   ext_modules=extensions,



)
