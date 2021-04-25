import os

project   = "Splines"
copyright = '2021, Enrico Bertolazzi'
author    = ':email:`Enrico Bertolazzi <enrico.bertolazzi@unitn.it>`'
version   = os.popen('git describe --tags --abbrev=0').read()

# The master toctree document.
master_doc = 'index'

# Add any paths that contain templates here, relative to this directory.
templates_path = [ '../_templates' ]
html_css_files = [ 'cloud.css' ]
html_static_path = ['../_static']


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The `extensions` list should already be in here from `sphinx-quickstart`
extensions = [
    # standard sphinx extensions
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',

    # 3rd party extensions
    'sphinxcontrib.fulltoc',
    'sphinx.ext.viewcode',  # mainly to make sure layout works properly

    # cloud's extensions
    'cloud_sptheme.ext.autodoc_sections',
    'cloud_sptheme.ext.relbar_links',
    'cloud_sptheme.ext.escaped_samp_literals',
    'cloud_sptheme.ext.issue_tracker',
    'cloud_sptheme.ext.table_styling',
    'cloud_sptheme.ext.role_index',  # NOTE: used only to provide example role index

    'sphinx.ext.mathjax',
    'sphinxcontrib.email',
]

source_suffix = ['.rst', '.md']

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

html_theme = 'cloud'
html_logo  = '../Spline_interpolation.svg'

email_automode       = True
autodoc_member_order = 'bysource'

html_theme_options = {
  "lighter_header_decor" : False,
  "borderless_decor"     : False
}

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
##html_sidebars = {'**': ['searchbox.html', 'globaltoc.html']}
html_sidebars = { '**': ['searchbox.html', 'globaltoc.html'] }

# If false, no module index is generated.
# html_domain_indices = False

# If false, no index is generated.
# html_use_index = False

# If true, the index is split into individual pages for each letter.
# html_split_index = True

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "%s v%s" % (project, version)

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "%s" % (project)

html_copy_source     = False
html_show_sourcelink = False
