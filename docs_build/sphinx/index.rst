Splines
=======

Introduction
------------

``Splines`` is a set of C++ classes (with MATLAB mex interface) which
implements varios spline interpolation.

- `Github repository <https://github.com/ebertolazzi/Splines>`__
- `Matlab Toolbox <https://github.com/ebertolazzi/Splines/releases>`__

Matlab Toolbox
--------------

To use in MATLAB install the toolbox ``Splines.mltbx`` then compile the
necessary ``mex`` files running ``CompileSplinesLib``.

C++ Usage
---------

The usage is simple:

.. code:: cpp

   #include "Splines.hh"
   using namespace SplinesLoad;

   // ....

   CubicSpline spline;
   double x[] = {1,2,3,4};
   double y[] = {3,1,1,3};
   spline.build(x,y,4); // build a cubic spline with 4 points

   cout << spline(1.1)     << '\n'; // spline at x = 1.1
   cout << spline.D(1.1)   << '\n'; // spline first derivative at x = 1.1
   cout << spline.DD(1.1)  << '\n'; // spline second derivative at x = 1.1
   cout << spline.DDD(1.1) << '\n'; // spline third derivative at x = 1.1

splines can be built incrementally

.. code:: cpp

  #include "Splines.hh"
  using namespace SplinesLoad;

  // ....

  CubicSpline spline;

  spline.pushBack( 1, 3 );
  spline.pushBack( 2, 1 );
  spline.pushBack( 3, 1 );
  spline.pushBack( 4, 3 );
  spline.build();

  cout << spline(1.1)     << '\n'; // spline at x = 1.1
  cout << spline.D(1.1)   << '\n'; // spline first derivative at x = 1.1
  cout << spline.DD(1.1)  << '\n'; // spline second derivative at x = 1.1
  cout << spline.DDD(1.1) << '\n'; // spline third derivative at x = 1.1

or by using standard vector

.. code:: cpp

  #include "Splines.hh"
  #include <vector>
  using namespace SplinesLoad;
  using namespace std;

  // ....

  CubicSpline spline;
  std::vector x, y;
  x.push_back(1); y.push_back(3);
  x.push_back(2); y.push_back(1);
  x.push_back(3); y.push_back(1);
  x.push_back(4); y.push_back(3);
  spline.build(x,y);

  cout << spline(1.1)     << '\n'; // spline at x = 1.1
  cout << spline.D(1.1)   << '\n'; // spline first derivative at x = 1.1
  cout << spline.DD(1.1)  << '\n'; // spline second derivative at x = 1.1
  cout << spline.DDD(1.1) << '\n'; // spline third derivative at x = 1.1

Compile and tests
-----------------

Using makefile
~~~~~~~~~~~~~~

Edit makefile file to match compiler of your OS and do:

.. code::  bash

  make

Using rakefile
~~~~~~~~~~~~~~

.. code:: sh

  rake build_win    # on windows
  rake build_linux  # on linux
  rake build_osx    # on mac

To run the test

.. code:: sh

  make run     # using makefile
  rake run     # using rake on linux and osx
  rake run_win # using rake on windows

Contents
--------

.. toctree::
  :maxdepth: 2

  api-cpp/root.rst
  api-c/root.rst
  api-matlab/root.rst

.. include:: author.rst

.. include:: references.rst

License
~~~~~~~

.. literalinclude:: ../../License.txt
