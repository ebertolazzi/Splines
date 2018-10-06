/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2016                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                | 
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Splines.hh"
#include <cmath>

/**
 * 
 */

namespace Splines {

  using namespace std; // load standard namspace

  /*\
   |    ___        _       _   _      ____        _ _
   |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
   |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
   |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                       |_|
   |
  \*/

  /*
  static
  inline
  int
  signTest( real_type a, real_type b ) {
    int sa = 0;
    int sb = 0;
    if      ( a > 0 ) sa =  1;
    else if ( a < 0 ) sa = -1;
    if      ( b > 0 ) sb =  1;
    else if ( b < 0 ) sb = -1;
    return sa*sb;
  }
  */

  static
  inline
  real_type
  minMod( real_type a, real_type b ) {
    if      ( a*b < 0 )         return 0;
    else if ( abs(a) < abs(b) ) return a;
    else                        return b;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  quintic_pchip( real_type const X[],
                 real_type const Y[],
                 real_type       Yp[],
                 real_type       Ypp[],
                 integer         n ) {
    pchip( X, Y, Yp, n );
    real_type h1 = X[1] - X[0];
    Ypp[0] = 0*(6*(X[1]-X[0])/h1-4*Yp[0]-2*Yp[1])/h1; // left
    for ( size_t i = 1; i < size_t(n); ++i ) {
      real_type h2 = X[i+1] - X[i];
      real_type Ypp_L = (6*(X[i-1]-X[i])/h1+4*Yp[i]+2*Yp[i-1])/h1;
      real_type Ypp_R = (6*(X[i+1]-X[i])/h2-4*Yp[i]-2*Yp[i+1])/h2;
      Ypp[i] = minMod(Ypp_L,Ypp_R);
      //Ypp[i] = 2*(h2*Y[i-1]+h1*Y[i+1]-(h1+h2)*Y[i])/((h1+h2)*h1*h2);
      h1 = h2;
    }
    Ypp[n] = 0*(6*(X[n-1]-X[n])/h1+4*Yp[n]+2*Yp[n-1])/h1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSpline::build() {
    SPLINE_ASSERT( npts > 1,
                   "QuinticSpline::build(): npts = " << npts <<
                   " not enought points" );
    integer ibegin = 0;
    integer iend   = 0;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < npts && X[iend-1] < X[iend] ) {};
      quintic_pchip( X+ibegin, Y+ibegin, Yp+ibegin, Ypp+ibegin, (iend-ibegin)-1 );
      ibegin = iend;
    } while ( iend < npts );
    
    SPLINE_CHECK_NAN( Yp,  "QuinticSpline::build(): Yp",  npts );
    SPLINE_CHECK_NAN( Ypp, "QuinticSpline::build(): Ypp", npts );
  }
}
