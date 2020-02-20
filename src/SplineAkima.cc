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
/*\
 |     #
 |    # #   #    # # #    #   ##
 |   #   #  #   #  # ##  ##  #  #
 |  #     # ####   # # ## # #    #
 |  ####### #  #   # #    # ######
 |  #     # #   #  # #    # #    #
 |  #     # #    # # #    # #    #
\*/

namespace Splines {

  using namespace std; // load standard namspace

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  real_type
  akima_one(
    real_type epsi,
    real_type di_m2,
    real_type di_m1,
    real_type di,
    real_type di_p1
  ) {
    real_type wl  = std::abs(di_p1 - di);
    real_type wr  = std::abs(di_m1 - di_m2);
    real_type den = wl + wr;
    if ( den <= epsi ) { wl = wr = 0.5; den = 1; } // if epsi == 0
    real_type num = wl * di_m1 + wr * di;
    return num / den;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Akima_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts
  ) {

    if ( npts == 2 ) { // solo 2 punti, niente da fare
      Yp[0] = Yp[1] = (Y[1]-Y[0])/(X[1]-X[0]);
    } else {
      std::vector<real_type> m;
      m.resize( size_t(npts+3) );

      // calcolo slopes (npts-1) intervals + 4
      for ( size_t i = 1; i < size_t(npts); ++i )
        m[i+1] = (Y[i]-Y[i-1])/(X[i]-X[i-1]);

      // extra slope at the boundary
      m[1] = 2*m[2]-m[3];
      m[0] = 2*m[1]-m[2];
      m[size_t(npts+1)] = 2*m[size_t(npts)]-m[size_t(npts-1)];
      m[size_t(npts+2)] = 2*m[size_t(npts+1)]-m[size_t(npts)];

      // minimum delta slope
      real_type epsi = 0;
      for ( size_t i = 0; i < size_t(npts+2); ++i ) {
        real_type dm = std::abs(m[i+1]-m[i]);
        if ( dm > epsi ) epsi = dm;
      }
      epsi *= 1E-8;

      // 0  1  2  3  4---- n-1 n n+1 n+2
      //       +  +  +      +  +
      for ( size_t i = 0; i < size_t(npts); ++i )
        Yp[i] = akima_one( epsi, m[i], m[i+1], m[i+2], m[i+3] );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  AkimaSpline::build() {
    SPLINE_ASSERT(
      npts > 1,
      "AkimaSpline::build(): npts = " << npts << " not enought points"
    )
    integer ibegin = 0;
    integer iend   = 0;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < npts && X[iend-1] < X[iend] ) {}
      Akima_build( X+ibegin, Y+ibegin, Yp+ibegin, iend-ibegin );
      ibegin = iend;
    } while ( iend < npts );
    
    SPLINE_CHECK_NAN( Yp, "AkimaSpline::build(): Yp", npts );
  }

}
