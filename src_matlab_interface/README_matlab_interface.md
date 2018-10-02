MATLAB interface for C++ Splines library
=======

`Splines` is a set of C++ classes which implements various spline interpolation.

**Compilation**

To compile the mex files run the script `Compile.m` in MATLAB.
After the compilation the scripts `spline1d` for univariate
spline and `spline2d` for bivariate spline are available.

**Univariate splines**

To build a univariate spline use to command

~~~
S = Spline1D( 'type', X, Y ); 

S = Spline1D( 'type' ); 
S.build( X, Y ); 
~~~

- `type` is a string with the type of spline you are building.
         Possibile values are:
         linear, cubic, akima, bessel, pchip, quintic.
- X, Y   two vectors with the x and y coordinates of the point 
         to be interpolated, X must be increasing.

after the build spline can be used as follows

~~~
% build the spline
S = Sline1D('pippo','akima',X,Y);

...
x = [1,2,3];
y = S.eval(X); % compute values
~~~

it is possibile to get also the derivative of the spline

~~~
dy  = S.eval_D(X);
ddy = S.eval_DD(X);
~~~

where `dy` is a vector with the first derivative and 
`ddy` is a vector with the second derivative.

**Bivariate splines**

To build a bivariate spline use to command

~~~
S = Spline2D( 'type', X, Y, Z );

S = Spline2D( 'type' );
S.build( X, Y, Z );
~~~

- `type` is a string with the type of spline you are building.
         Possibile values are:
         bilinear, bicubic, akima, biquintic.
- X, Y   two vectors with the x and y coordinates of the mesh 
         to be interpolated, X and Y must be increasing.
- Z      matrix of dimension `length(X)` x `length(Y)`
         where `Z(i,j)` is the z elevation at point
         `(X(i),Z(j))`.

after the build spline can be used as follows

~~~
% build the spline
S = Spline2D('akima',X,Y,Z);

...
x       = [1,2,3];
y       = [2,3,4];
[xx,yy] = ndgrid(x,y);
zz      = S.eva(xx,yy); % compute values
% z is a matrix 3x3 with z(i,j) the spline S(x(i),y(j))
~~~

it is possibile to get also the derivative of the spline

~~~
dx = S.eval_Dx(X,Y);
dy = S.eval_Dy(X,Y);

dxx = S.eval_Dxx(X,Y);
dxy = S.eval_Dxy(X,Y);
dyy = S.eval_Dyy(X,Y);
~~~

where

- `dx` derivative respect to `x`
- `dy` derivative respect to `y`
- `dxx` second derivative respect to `x`
- `dxy` second derivative respect to `x` and `y`
- `dyy` second derivative respect to `y`

* * *

Enrico Bertolazzi<br>
Dipartimento di Ingegneria Industriale<br>
Universita` degli Studi di Trento<br>
email: enrico.bertolazzi@unitn.it
