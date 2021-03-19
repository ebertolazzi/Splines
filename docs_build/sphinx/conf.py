# -*- coding: utf-8 -*-

# pip3 install recommonmark
# pip3 install exhale
# pip3 install breathe
# pip3 install pydata-sphinx-theme
# pip3 install sphinx-markdown-parser
# pip3 install pymdown-extensions
# pip3 install m2r2
# pip3 install sphinxcontrib-email

#
#import recommonmark
#from recommonmark.transform import AutoStructify

#source_suffix = {
#  ".rst": "restructuredtext",
#  ".md": "markdown"
#}

#source_parsers = {
#  ".md": "recommonmark.parser.CommonMarkParser"
#}

# The master toctree document.
master_doc = 'index'

# -- Project information -----------------------------------------------------

project   = 'Splines'
copyright = '2021, :email:`Enrico Bertolazzi'
author    = ':email:`Enrico Bertolazzi <enrico.bertolazzi@unitn.it>`'

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

# Setup the exhale extension
exhale_args = {
  # These arguments are required
  "containmentFolder":     "./api",
  "rootFileName":          "library_root.rst",
  "rootFileTitle":         "Splines API",
  "doxygenStripFromPath":  "..",
  # Suggested optional arguments
  "createTreeView":        True,
  # TIP: if using the sphinx-bootstrap-theme, you need
  # "treeViewIsBootstrap": True,
  "exhaleExecutesDoxygen": True,
  #"exhaleDoxygenStdin":    "INPUT = ../../src"
  "exhaleDoxygenStdin":
'''
        EXTRACT_ALL       = YES
        SOURCE_BROWSER    = YES
        EXTRACT_STATIC    = YES
        HIDE_SCOPE_NAMES  = YES
        QUIET             = YES
        INPUT             = ../../src
        EXAMPLE_RECURSIVE = YES
        GENERATE_TREEVIEW = YES
'''
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

html_theme = 'pydata_sphinx_theme'
html_logo  = '../logo.png'

email_automode = True

# app setup hook
#def setup(app):
#    app.add_config_value('recommonmark_config', {
#        #'url_resolver': lambda url: github_doc_root + url,
#        'auto_toc_tree_section': 'Contents',
#        'enable_math': False,
#        'enable_inline_math': False,
#        'enable_eval_rst': True,
#        'enable_auto_doc_ref': True,
#    }, True)
#    app.add_transform(AutoStructify)
