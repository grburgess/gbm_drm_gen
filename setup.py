import os


from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

ext1 = Extension(name='ftran',
                 sources=['src/sig_ftran.pyf', 'src/ftran.f90'],
                 libraries=['m'],
                 # extra_f77_compile_args=["-fimplicit-none", "-ffree-line-length-none","-O3","-Ofast"],
                 extra_f90_compile_args=["-fimplicit-none", "-ffree-line-length-none", "-Ofast"])
ext2 = Extension(name='at_scat',
                 sources=['src/at_scat.pyx'],
                 include_dirs=[numpy.get_include()])

extensions = [ext1]
extensions_cython = [ext2]



# Create list of data files
def find_data_files(directory):

    paths = []

    for (path, directories, filenames) in os.walk(directory):

        for filename in filenames:

            paths.append(os.path.join('..', path, filename))

    return paths

extra_files = find_data_files('gbm_drm_gen/data')


setup(

    ext_modules=cythonize(extensions_cython), requires=['numpy']
)


setup(

    name="gbm_drm_gen",

    packages=['gbm_drm_gen',
              'gbm_drm_gen/utils',
              'gbm_drm_gen/io'],
    version='v0.1',
    description='GBM DRM generator',
    author='J. Michael Burgess',
    author_email='jmichaelburgess@gmail.com',
    package_data={'': extra_files, },
    include_package_data=True,

    install_requires=[
        'numpy',
        'astropy',
        'scipy',
        'h5py',
        'threeML',
    ],
    # dependency_links = [
    #  "git+git://github.com/giacomov/3ML.git/#egg=threeML",
    # ],

    ext_modules=extensions


)

