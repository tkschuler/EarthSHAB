from distutils.core import setup
from setuptools import find_packages

setup(
    name='EarthSHAB',
    version='2.2.0',
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[],
)