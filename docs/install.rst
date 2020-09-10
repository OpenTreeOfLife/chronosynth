######################
Installing ChronoSynth
######################


Downloading the Source Code Repository
======================================

The ChronoSynth source code is version-controlled using |Git|_, and available at <https://github.com/OpenTreeOfLife/ChronoSynth>`_.
Make sure you have |Git|_ for the following to work:

    - `Source <http://www.kernel.org/pub/software/scm/git/git-1.6.6.tar.gz>`_
    - `OS X binaries <http://code.google.com/p/git-osx-installer/downloads/list?can=3>`_
    - `Microsoft Windows <http://code.google.com/p/msysgit/downloads/list>`_

Download ChronoSynth by running::

    $ git clone https://github.com/OpenTreeOfLife/ChronoSynth.git

To use ChronoSynth's library directly, install it in developer mode by running::

    $ cd ChronoSynth
    $ python3 -m pip install -e .
