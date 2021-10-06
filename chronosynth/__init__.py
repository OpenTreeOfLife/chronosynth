#!/usr/bin/env python3
__version__ = "0.0.1"  # sync with setup.py

import logging
import sys

logging.basicConfig(filename='chronosynth.log', 
                    filemode='w', 
                    level=logging.DEBUG,
                    format='%(name)s - %(levelname)s - %(message)s')

