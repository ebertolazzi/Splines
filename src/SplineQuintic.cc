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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
     Sistema lineare da risolvere

     D U
     L D U
       L D U
         L D U
           .....
              L D U
                L D U
                  L D

  \*/

  static
  void
  QuinticSpline_build(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    real_type       Ypp[],
    integer         npts
  ) {

    size_t n = size_t(npts > 0 ? npts-1 : 0);

    vector<real_type> buffer(3*(n+1));
    real_type * ptr = &buffer.front();
    real_type * L = ptr; ptr += npts;
    real_type * D = ptr; ptr += npts;
    real_type * U = ptr;
    real_type * Z = Ypp;

    size_t i;
    for ( i = 1; i < n; ++i ) {
      real_type h__L = X[i] - X[i-1];
      real_type h__R = X[i+1] - X[i];
      real_type p0   = Y[i-1];
      real_type p1   = Y[i];
      real_type p2   = Y[i+1];
      real_type dp0  = Yp[i-1];
      real_type dp1  = Yp[i];
      real_type dp2  = Yp[i+1];
      real_type t12  = 1/h__L;
      real_type t13 = 1/h__R;
      L[i] = -t12;
      D[i] = 3*(t13+t12);
      U[i] = -t13;
      real_type t15 = h__L*h__L;
      real_type t16 = 1/t15;
      real_type t19 = h__R*h__R;
      real_type t20 = 1/t19;
      real_type t27 = t16/h__L;
      real_type t31 = t20/h__R;
      Z[i] = 8*(dp0*t16-dp2*t20) - 12*dp1*(t20-t16) + 20*(p0*t27-p1*(t31+t27)+p2*t31);
    }
    L[0] = U[0] = 0; D[0] = 1;
    L[n] = U[n] = 0; D[n] = 1;

    {
      real_type H = X[1] - X[0];
      real_type base_DD[4];
      Hermite3_DD( 0, H, base_DD );
      Z[0] = base_DD[0] * Y[0]  + base_DD[1] * Y[1] +
             base_DD[2] * Yp[0] + base_DD[3] * Yp[1];
    }
    {
      real_type H = X[n] - X[n-1];
      real_type base_DD[4];
      Hermite3_DD( H, H, base_DD );
      Z[n] = base_DD[0] * Y[n-1]  + base_DD[1] * Y[n] +
             base_DD[2] * Yp[n-1] + base_DD[3] * Yp[n];
    }

    i = 0;
    do {
      Z[i]   /= D[i];
      U[i]   /= D[i];
      D[i+1] -= L[i+1] * U[i];
      Z[i+1] -= L[i+1] * Z[i];
    } while ( ++i < n );

    Z[i] /= D[i];

    do {
      --i;
      Z[i] -= U[i] * Z[i+1];
    } while ( i > 0 );
  }

  /*\
   |    ___        _       _   _      ____        _ _
   |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
   |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
   |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                       |_|
   |
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  Quintic_build(
    QUINTIC_SPLINE_TYPE q_sub_type,
    real_type const     X[],
    real_type const     Y[],
    real_type           Yp[],
    real_type           Ypp[],
    integer             npts
  ) {

    switch ( q_sub_type ) {
    case CUBIC_QUINTIC:
      CubicSpline_build( X, Y, Yp, npts, NOT_A_KNOT, NOT_A_KNOT );
      break;
    case PCHIP_QUINTIC:
      Pchip_build( X, Y, Yp, npts );
      break;
    case AKIMA_QUINTIC:
      Akima_build( X, Y, Yp, npts );
      break;
    case BESSEL_QUINTIC:
      Bessel_build( X, Y, Yp, npts );
      break;
    }
    QuinticSpline_build( X, Y, Yp, Ypp, npts );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSpline::build() {
    SPLINE_ASSERT(
      this->npts > 1,
      "QuinticSpline::build(): npts = " << this->npts << " not enought points"
    )
    integer ibegin = 0;
    integer iend   = 0;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < this->npts && this->X[iend-1] < this->X[iend] ) {}
      Quintic_build(
        this->q_sub_type,
        this->X+ibegin,  this->Y+ibegin,
        this->Yp+ibegin, this->Ypp+ibegin,
        iend - ibegin
      );
      ibegin = iend;
    } while ( iend < this->npts );
    
    SPLINE_CHECK_NAN( this->Yp,  "QuinticSpline::build(): Yp",  this->npts );
    SPLINE_CHECK_NAN( this->Ypp, "QuinticSpline::build(): Ypp", this->npts );
  }
}
