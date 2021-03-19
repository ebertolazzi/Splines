Splines
=======

[![Build Status](https://travis-ci.org/ebertolazzi/Splines.svg?branch=master)](https://travis-ci.org/ebertolazzi/Splines)
[![View Splines on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/54481-splines)

<br>

[Splines](https://github.com/ebertolazzi/Splines) 
is a set of C++ classes (with MATLAB mex interface) which 
implements varios spline interpolation.
The classes are the following:
 
  - ConstantSpline, for piecewise constants functions
  - LinearSpline, for piecewise linear interpolation
  - CubicSpline, for classical cubic spline interpolation
  - AkimaSpline, for Akima "non oscillatory" spline interpolation 
  - BesselSpline, for Bessel "non oscillatory" spline interpolation 
  - PchipSpline, 
  - QuinticSpline, Simple quintic spline based on PCHIP

### References

- F.N. Fritsch and R.E. Carlson,
  Monotone Piecewise Cubic Interpolation,<br>
  SIAM Journal of Numerical Analysis, Vol. 17, No. 2, pp. 238-246,
  April 1980.
  
### Matlab

To use in MATLAB install the toolbox `Splines.mltbx` then compile the files running `CompileSplinesLib` (available at [releases](https://github.com/ebertolazzi/Splines/releases))

 
### C++ Usage

The usage is simple:

~~~~~~~~~~~~~
#include "Splines.hh"
using namespace SplinesLoad;

....

CubicSpline spline;
double x[] = {1,2,3,4};
double y[] = {3,1,1,3};
spline.build(x,y,4); // build a cubic spline with 4 points
  
cout << spline(1.1) << '\n';     // spline at x = 1.1
cout << spline.D(1.1) << '\n';   // spline first derivative at x = 1.1
cout << spline.DD(1.1) << '\n';  // spline second derivative at x = 1.1
cout << spline.DDD(1.1) << '\n'; // spline third derivative at x = 1.1
~~~~~~~~~~~~~

splines can be built incrementally 

~~~~~~~~~~~~~
#include "Splines.hh"
using namespace SplinesLoad;

....

CubicSpline spline;
  
spline . pushBack( 1, 3 );
spline . pushBack( 2, 1 );
spline . pushBack( 3, 1 );
spline . pushBack( 4, 3 );
spline . build();
  
cout << spline(1.1) << '\n';     // spline at x = 1.1
cout << spline.D(1.1) << '\n';   // spline first derivative at x = 1.1
cout << spline.DD(1.1) << '\n';  // spline second derivative at x = 1.1
cout << spline.DDD(1.1) << '\n'; // spline third derivative at x = 1.1
~~~~~~~~~~~~~

or by using standard vector 

~~~~~~~~~~~~~
#include "Splines.hh"
#include <vector>
using namespace SplinesLoad;
using namespace std;

....

CubicSpline spline;
std::vector x, y;
x.push_back(1); y.push_back(3);
x.push_back(2); y.push_back(1);
x.push_back(3); y.push_back(1);
x.push_back(4); y.push_back(3);
spline . build(x,y);
  
cout << spline(1.1) << '\n';     // spline at x = 1.1
cout << spline.D(1.1) << '\n';   // spline first derivative at x = 1.1
cout << spline.DD(1.1) << '\n';  // spline second derivative at x = 1.1
cout << spline.DDD(1.1) << '\n'; // spline third derivative at x = 1.1
~~~~~~~~~~~~~

### Compile and tests

**Using makefile**

Edit makefile file to match compiler of your OS and do:

```sh
  make
```

**Using rakefile**

```sh
  rake build_win    # on windows
  rake build_linux  # on linux
  rake build_osx    # on mac
```

To run the test

```sh
  make run     # using makefile
  rake run     # using rake on linux and osx
  rake run_win # using rake on windows
```

To generate documentation (using [DOXYGEN](https://www.doxygen.nl/index.html) and [SPHINX](https://www.sphinx-doc.org/en/master/))

```sh
  make doc
```

### Online Documentation

Available at: [http://ebertolazzi.github.io/Splines](http://ebertolazzi.github.io/Splines)

* * *

Enrico Bertolazzi<br>
Dipartimento di Ingegneria Industriale<br>
Universita` degli Studi di Trento<br>
email: enrico.bertolazzi@unitn.it

* * *
