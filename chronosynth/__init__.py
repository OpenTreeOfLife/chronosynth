#!/usr/bin/env python3
__version__ = "0.0.1"  # sync with setup.py

import logging.config
import sys
import os

#https://medium.com/@narengowda/better-pyramid-logging-with-python-and-log-configurationsyoure-up-and-running-c47203dfaff2

configfile = chronosynth_dir = os.path.dirname(__file__) + '/../log.config'
logging.config.fileConfig(configfile)
