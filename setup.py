from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
from numpy.distutils.core import Extension # as Numpy_Extension
from numpy.distutils.core import setup
import os
import io
import numpy
import versioneer

ext1 = Extension(
    name='ftran',
    sources=['src/sig_ftran.pyf', 'src/ftran.f90'],
    libraries=['m'],
    extra_f90_compile_args=[
        "-fimplicit-none", "-ffree-line-length-none", "-Ofast"
    ])

ext2 = Extension(name='at_scat',
                 sources=['src/at_scat.pyx'],
                 include_dirs=[numpy.get_include()])

extensions = cythonize([ext1, ext2])



# Create list of data files
def find_data_files(directory):

    paths = []

    for (path, directories, filenames) in os.walk(directory):

        for filename in filenames:

            paths.append(os.path.join('..', path, filename))

    return paths


extra_files = find_data_files('gbm_drm_gen/data')


setup(
    ext_modules=extensions,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    include_package_data=True,
    package_data={"": extra_files},

)
