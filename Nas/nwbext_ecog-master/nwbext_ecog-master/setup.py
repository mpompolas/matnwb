from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path


here = path.abspath(path.dirname(__file__))


with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='nwbext_ecog',
    version='0.5.0',
    description='Convert data to nwb',
    long_description=long_description,
    author='Ben Dichter',
    author_email='ben.dichter@gmail.com',
    keywords=['nwb', 'extension'],
    packages=find_packages(),
    install_requires=['pynwb'],
    package_data={'': ['ecog.namespace.yaml']},
    include_package_data=True,
    #entry_points={'pynwb.extensions': 'ecog = nwbext_ecog'},
)
