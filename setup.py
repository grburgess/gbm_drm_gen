import os

# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if not on_rtd:
    from numpy.distutils.core import setup
    from numpy.distutils.core import Extension
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
    import numpy

    ext1 = Extension(name='ftran',
                     sources=['src/sig_ftran.pyf', 'src/ftran.f90'],
                     libraries=['m'],
                     # extra_f77_compile_args=["-fimplicit-none", "-ffree-line-length-none","-O3","-Ofast"],
                     extra_f90_compile_args=["-fimplicit-none", "-ffree-line-length-none", " -O3", "-Ofast"])
    ext2 = Extension(name='at_scat',
                     sources=['src/at_scat.pyx'],
                     include_dirs=[numpy.get_include()])

    extensions = [ext1]
    extensions_cython = [ext2]

    setup(

        ext_modules=cythonize(extensions_cython), requires=['numpy']
    )


else:
    print "On RTD"
    from setuptools import setup
    from distutils.core import Extension

    extensions = []
    extensions_cython = []

setup(

    name="balrog",
    packages=['balrog'],
    version='v0.1',
    description=' BAyesian Location Reconstruction of GRBs ',
    author='J. Michael Burgess',
    author_email='jmichaelburgess@gmail.com',

    requires=[
        'numpy',
        'matplotlib'
        'astropy',
        'pymultinest',
        'scipy',
    ],

    ext_modules=extensions
)

