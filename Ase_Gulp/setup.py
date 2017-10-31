#!/usr/bin/env python

from distutils.core import setup

setup (name='ase-gulp',
      version='0.1',
      description='A GULP calculator for ase',
      author='Matthew Dyer',
      author_email='msd30@liv.ac.uk',
      packages=['ase/io','ase/calculators']
)

