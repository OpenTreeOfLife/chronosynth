# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import time
from sphinx.ext import autodoc
# from chronosynth import __version__ as PROJECT_VERSION
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'chronosynth'
copyright = '2020, Luna Luisa Sanchez Reyes, Emily Jane McTavish'
author = 'Luna Luisa Sanchez Reyes, Emily Jane McTavish'

# The full version, including alpha/beta/rc tags
release = '0.0.1'

master_doc = 'index'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinxcontrib.napoleon' # requires: pip install sphinxcontrib-napoleon
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'nature'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Sphinx Hackery ------------------------------------------------

# Following allows for a docstring of a method to be inserted "nakedly"
# (without the signature etc.) into the current context by, for example::
#
#       .. autodocstringonly:: dendropy.dataio.newickreader.NewickReader.__init__
#
# Based on:
#
#   http://stackoverflow.com/questions/7825263/including-docstring-in-sphinx-documentation
# Borrowed from Dendropy's docs/conf.py
class DocStringOnlyMethodDocumenter(autodoc.MethodDocumenter):
    objtype = "docstringonly"

    # do not indent the content
    content_indent = "    "

    # do not add a header to the docstring
    def add_directive_header(self, sig):
        pass

    # def add_line(self, line, source, *lineno):
    #     """Append one line of generated reST to the output."""
    #     print self.indent + line
    #     self.directive.result.append(self.indent + line, source, *lineno)

class KeywordArgumentsOnlyMethodDocumenter(autodoc.MethodDocumenter):
    objtype = "keywordargumentsonly"
    priority = 0 # do not override normal autodocumenter

    # do not indent the content
    content_indent = "    "

    # do not add a header to the docstring
    def add_directive_header(self, sig):
        pass

    def add_line(self, line, source, *lineno):
        if ":Keyword Arguments:" in line:
            line = line.replace(":Keyword Arguments:", "                   ")
            self._emit_line = True
        if getattr(self, "_emit_line", False):
            self.directive.result.append(self.indent + line, source, *lineno)


# General information about the project.
project = u'ChronoSynth'
copyright = u'2020-{}, Luna L. Sanchez Reyes and Emily Jane McTavish'.format(time.strftime('%Y'))

# Shortcuts
rst_prolog = """
.. |Git| replace:: Git
.. _Git: http://git-scm.com/
"""
