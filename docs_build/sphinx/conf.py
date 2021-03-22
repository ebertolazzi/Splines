# -*- coding: utf-8 -*-

# pip3 install recommonmark
# pip3 install exhale
# pip3 install breathe
# pip3 install pydata-sphinx-theme
# pip3 install sphinx-markdown-parser
# pip3 install pymdown-extensions
# pip3 install m2r2
# pip3 install sphinxcontrib-email
# pip3 install furo
# pip3 install faculty-sphinx-theme
# pip3 install install sphinx_sizzle_theme
# pip3 install karma_sphinx_theme
# pip3 install sphinx-book-theme
# pip3 install myst-parser
# https://pradyunsg.me/furo/

import os

# The master toctree document.
master_doc = 'index'

# -- Project information -----------------------------------------------------

project   = 'Quartic Roots'
copyright = '2021, Enrico Bertolazzi'
author    = ':email:`Enrico Bertolazzi <enrico.bertolazzi@unitn.it>`'
version   = os.popen('git describe --tags').read()

#rst_epilog =
#rst_prolog =

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The `extensions` list should already be in here from `sphinx-quickstart`
extensions = [
  'breathe',
  'exhale',
  #'recommonmark', # non funziona
  'm2r2',          # funziona!
  #'myst_parser',  # non funziona
  'sphinx.ext.autodoc',
  'sphinx.ext.doctest',
  'sphinx.ext.coverage',
  'sphinx.ext.mathjax',
  'sphinx.ext.ifconfig',
  'sphinx.ext.githubpages',
  'sphinx.ext.intersphinx',
  'sphinxcontrib.email'
]

source_suffix = ['.rst', '.md']

# Setup the breathe extension
breathe_projects = {
  "Splines": "../xml"
}
breathe_default_project = "Splines"


# somewhere in `conf.py`, *BEFORE* declaring `exhale_args`
def specificationsForKind(kind):
    '''
    For a given input ``kind``, return the list of reStructuredText specifications
    for the associated Breathe directive.
    '''
    # Change the defaults for .. doxygenclass:: and .. doxygenstruct::
    if kind == "class" or kind == "struct":
        return [
          ":members:",
          ":protected-members:",
          ":private-members:",
          ":undoc-members:"
        ]
    # Change the defaults for .. doxygenenum::
    elif kind == "namespace":
        return [":no-link:"]
    # Change the defaults for .. doxygenenum::
    elif kind == "enum":
        return [":no-link:"]
    # An empty list signals to Exhale to use the defaults
    else:
        return []



# Setup the exhale extension
exhale_args = {
  # These arguments are required
  "containmentFolder":     "./api",
  "rootFileName":          "library_root.rst",
  "rootFileTitle":         "C/C++ API",
  "doxygenStripFromPath":  "..",
  # Suggested optional arguments
  "createTreeView":        True,
  # TIP: if using the sphinx-bootstrap-theme, you need
  # "treeViewIsBootstrap": True,
  "exhaleExecutesDoxygen": True,
  #"exhaleDoxygenStdin":    "INPUT = ../../src"
  "exhaleDoxygenStdin":
'''
        EXTRACT_ALL         = YES
        SOURCE_BROWSER      = NO
        EXTRACT_STATIC      = YES
        HIDE_SCOPE_NAMES    = NO
        CALLER_GRAPH        = YES
        GRAPHICAL_HIERARCHY = YES
        HAVE_DOT            = YES
        QUIET               = NO
        INPUT               = ../../src
        GENERATE_TREEVIEW   = YES

        XML_PROGRAMLISTING   = YES
        RECURSIVE            = YES
        FULL_PATH_NAMES      = YES
        ENABLE_PREPROCESSING = YES
        MACRO_EXPANSION      = YES
        SKIP_FUNCTION_MACROS = NO
        EXPAND_ONLY_PREDEF   = NO
''',
  'kindsWithContentsDirectives': [] # tolgo contents a tutte! (serve per Furo)
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

pygments_style      = "sphinx"
pygments_dark_style = "monokai"

#html_theme = 'pydata_sphinx_theme'
html_theme = 'furo'
#html_theme = "sizzle"
#html_theme = "karma_sphinx_theme"
#html_theme = "sphinx_book_theme"
html_logo  = '../logo.png'

email_automode = True
