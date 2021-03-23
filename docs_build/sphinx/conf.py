# -*- coding: utf-8 -*-
import os
from past.builtins import execfile

# -- Project information -----------------------------------------------------
execfile('../project_common.py')

extensions.append('breathe');
extensions.append('sphinx.ext.intersphinx');

#intersphinx_mapping = {
#  'c_interface': ('../../otherbook/build/html/objects.inv', None),
#}

# Setup the breathe extension
breathe_projects = {
  project: "../xml",
  project+"_c": "../xml-c",
  project+"_cpp": "../xml-cpp",
  project+"_matlab": "../xml-matlab"
}
breathe_default_project = project
