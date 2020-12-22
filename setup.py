from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
from setuptools import find_packages
from numpy.distutils.core import Extension
from numpy.distutils.core import setup
import os
import io
import numpy
import versioneer


ext1 = Extension(
    name="ftran",
    sources=["src/sig_ftran.pyf", "src/ftran.f90"],
    libraries=["m"],
    extra_f90_compile_args=["-fimplicit-none", "-ffree-line-length-none", "-Ofast"],
    include_dirs=[numpy.get_include()]
)

# ext2 = Extension(
#     name="at_scat", sources=["src/at_scat.pyx"], include_dirs=[numpy.get_include()]

ext2 = Extension(
    name="at_scat", sources=["src/at_scat.pyx"]
)

extensions = cythonize([ext1, ext2])

# extensions_cython = [ext2]

# Package meta-data.
NAME = "gbm_drm_gen"
DESCRIPTION = "A DRM generator for the Fermi GBM detectors"
URL = "https://github.com/mpe-grb/gbm_drm_gen"
EMAIL = "jmichaelburgess@gmail.com"
AUTHOR = "J. Michael Burgess"
REQUIRES_PYTHON = ">=2.7.0"
VERSION = None

REQUIRED = ["numpy", "astropy", "scipy", "h5py", "responsum", "gbmgeometry"]

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
    package_data={"": extra_files,},
    license="BSD",
    ext_modules=extensions,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
