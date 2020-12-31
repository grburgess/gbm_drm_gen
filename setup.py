from setuptools import setup
import os


import versioneer



extra_files = ["gbm_drm_gen/data/balrog_db.h5"]

setup(
    version=versioneer.get_version(),
   # include_package_data=True,
    package_data={"": extra_files},
    cmdclass=versioneer.get_cmdclass(),
)
