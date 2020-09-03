#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='chronosynth',
      version='0.0.1',
      description='ChronoSynth',
      author='Luna Luisa Sanchez Reyes, Emily Jane McTavish',
      author_email='ejmctavish@gmail.com',
      packages=['chronosynth'],
      install_requires=['opentree']
     )
