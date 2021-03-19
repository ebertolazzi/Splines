.. Splines documentation master file, created by
   sphinx-quickstart on Fri Mar 19 01:43:44 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Splines
=======

Splines is a set of C++ classes (with MATLAB mex interface) which implements various spline interpolation. The classes are the following:

- ConstantSpline, for piecewise constants functions

- LinearSpline, for piecewise linear interpolation

- CubicSpline, for classical cubic spline interpolation

- AkimaSpline, for Akima "non oscillatory" spline interpolation

- BesselSpline, for Bessel "non oscillatory" spline interpolation

- PchipSpline, Simple cubic spline based on PCHIP

- QuinticSpline, Simple quintic spline based on PCHIP

.. toctree::
   :maxdepth: 2

   readme.rst
   api/library_root.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

License
-------

.. literalinclude:: ../../License.txt
