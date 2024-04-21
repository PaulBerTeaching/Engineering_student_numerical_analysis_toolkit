# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='ESNAT',
    version='0.0.1',
    description='Nnumerical analysis toolkit for first years of engineering studies',
    long_description=readme,
    author='Paul BERRAUTE',
    author_email='paul.berraute@hotmail.fr',
    url='https://github.com/PaulBerTeaching',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

