# ChronoSynth


[![Build
Status](https://travis-ci.org/OpenTreeOfLife/ChronoSynth.svg?branch=main)](https://travis-ci.org/OpenTreeOfLife/ChronoSynth)[![Documentation](https://readthedocs.org/projects/chronosynth/badge/?version=latest&style=flat)](https://chronosynth.readthedocs.io/en/latest/)[![codecov](https://codecov.io/gh/OpenTreeOfLife/ChronoSynth/branch/main/graph/badge.svg)](https://codecov.io/gh/OpenTreeOfLife/ChronoSynth)[![NSF-1759846](https://img.shields.io/badge/NSF-1759846-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1759846) 

Get dates for nodes in synthetic tree from chronograms  in phylesystem

See [ReadTheDocs](https://chronosynth.readthedocs.io/en/latest/index.html) for the API, a work in progress. Currently the main dependency is on [python-opentree](https://github.com/OpenTreeOfLife/python-opentree), which is [used here](https://github.com/OpenTreeOfLife/ChronoSynth/blob/main/chronosynth/chronogram.py#L9) to [return the `OpenTree` singleton](https://github.com/OpenTreeOfLife/python-opentree/blob/main/opentree/ot_object.py#L51).


Configuration settings described in default.config. 
Copy to  "~/.chronosynth/config.ini" st set local config or 
'CHRONOSYNTH_CONFIG_FILE' in environment,
or it will use the values forn defaul.config.
