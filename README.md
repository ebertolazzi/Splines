# Splines

[![](https://travis-ci.org/ebertolazzi/Splines.svg?branch=master)](https://travis-ci.org/ebertolazzi/Splines) [![](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/54481-splines)

`Splines` is a set of `C++` classes (with `MATLAB` mex interface) which
implements various spline interpolation.

The library contains the following objects:

- Univariate curve
  - Linear Spline
  - Akima spline
  - Bessel spline
  - Cubic spline
  - Hermite spline
  - Pchic non oscillatory spline
  - Quintic spline
- Bivariate spline (surface)
  - Bilinear spline
  - Cubic
  - Quintic

Library is written in `C++11` with a `MATLAB` mex interface.
Thus can be used in fast compiled application or in `MATLAB` scripts.

To compile the `C++11` library the easy way require `cmake` and `rake`

```
ruby setup.rb
```

then

```
rake
```

for more details see: **online documentation** at http://ebertolazzi.github.io/Splines/