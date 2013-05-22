#!/usr/bin/env python

from distutils.core import setup

setup(
    name='Faps',
    version='2.0.0.0',
    description='Fully Automated Adsorption Analysis in Porous Solids',
    author='Tom Daff',
    author_email='tdaff@uottawa.ca',
    url='http://titan.chem.uottawa.ca/',
    packages=['faps',],
    license='BSD',
    long_description=open('README.rst').read(),
)
