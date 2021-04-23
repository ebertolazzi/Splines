# -*- coding: utf-8 -*-
import os
from pathlib import Path

# -- Project information -----------------------------------------------------
exec(open("../project_common.py").read())

# Setup the breathe extension
breathe_projects = { project+"_c": "../xml-c" }
breathe_default_project = project+"_c"

extensions.append('exhale');
extensions.append('breathe');

dir_path = os.path.dirname(os.path.realpath(__file__))+"../../../src"
dir_path = Path(dir_path).resolve()

# Setup the exhale extension
exhale_args = {
  # These arguments are required
  "containmentFolder":     "./api-c",
  "rootFileName":          "library_root.rst",
  "rootFileTitle":         "C API",
  "doxygenStripFromPath":  str(dir_path),
  # Suggested optional arguments
  "createTreeView":        False,
  # TIP: if using the sphinx-bootstrap-theme, you need
  "treeViewIsBootstrap":   False,
  "exhaleExecutesDoxygen": True,
  #"exhaleDoxygenStdin":    "INPUT = ../../src"
  "exhaleDoxygenStdin":
'''
        QUIET               = NO
        INPUT               = ../../src
        GENERATE_TREEVIEW   = YES
        XML_OUTPUT          = xml-c
        SHORT_NAMES         = YES

        PREDEFINED           += protected=private

        EXTRACT_ALL            = YES
        EXTRACT_STATIC         = YES
        SHORT_NAMES            = YES
        INHERIT_DOCS           = YES
        ENABLE_PREPROCESSING   = YES
        MACRO_EXPANSION        = NO
        XML_OUTPUT             = xml-c
        XML_PROGRAMLISTING     = NO
        XML_NS_MEMB_FILE_SCOPE = NO
        SOURCE_BROWSER         = NO
        OPTIMIZE_OUTPUT_FOR_C  = NO
        HIDE_SCOPE_NAMES       = NO
        SEARCH_INCLUDES        = NO
        CALLER_GRAPH           = YES
        GRAPHICAL_HIERARCHY    = NO
        HAVE_DOT               = NO
        SHOW_INCLUDE_FILES     = NO
        GENERATE_TREEVIEW      = YES

        FILE_PATTERNS         = SplinesCinterface.* *.h *.c
        GENERATE_HTML         = NO
''',
  "lexerMapping": { r".*\.m": "MATLAB" }
}

html_theme_options['logotarget'] = "../index"
html_theme_options['roottarget'] = "../index"