######################
Installing ChronoSynth
######################


Downloading the Source Code Repository
======================================

The ChronoSynth source code is version-controlled using |Git|_, and available on
a public GitHub `repository <https://github.com/OpenTreeOfLife/ChronoSynth>`_.
To download ChronoSynth you will need |Git|_. Install it following instructions from the
following links:

    - `Source <http://www.kernel.org/pub/software/scm/git/git-1.6.6.tar.gz>`_
    - `OS X binaries <http://code.google.com/p/git-osx-installer/downloads/list?can=3>`_
    - `Microsoft Windows <http://code.google.com/p/msysgit/downloads/list>`_

Now, you can download ChronoSynth by running::

    $ git clone https://github.com/OpenTreeOfLife/ChronoSynth.git

Then, create a virtual environment and activate it::

    $ cd ChronoSynth
    $ virtualenv -p python3 venv-chronosynth
    $ source venv-chronosynth/bin/activate


Finally, install ChronoSynth in developer mode by running::

    $ pip install -r requirements.txt
    $ pip install -e  .


After you are finished working with ChronoSynth and you don't want to run it anymore, deactivate the virtual environment with:


    $ deactivate
    
