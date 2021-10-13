#!/usr/bin/env python3
__version__ = "0.0.1"  # sync with setup.py

import logging.config
import sys
import os

#https://medium.com/@narengowda/better-pyramid-logging-with-python-and-log-configurationsyoure-up-and-running-c47203dfaff2


if 'CHRONOSYNTH_CONFIG_FILE' in os.environ:
    cfn = os.path.abspath(os.environ['CHRONOSYNTH_CONFIG_FILE'])
else:
    cfn = os.path.expanduser("~/.chrnosynth/config")

if not os.path.isfile(cfn):
    if 'CHRONOSYNTH_CONFIG_FILE' in os.environ:
        sys.stderr.write('Filepath "{}" specified via CHRONOSYNTH_CONFIG_FILE={} was not found'.format(cfn, os.environ[
                    'CHRONOSYNTH_CONFIG_FILE']))

cfn = os.path.dirname(__file__) + '/default.config'

if not os.path.isfile(cfn):
    raise RuntimeError('The peyotl configuration file cascade failed looking for "{}"'.format(cfn))

configfile = os.path.abspath(cfn)

logging.config.fileConfig(configfile)
