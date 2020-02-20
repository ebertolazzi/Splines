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

/*\
 |   ____                     _ ____        _ _
 |  | __ )  ___  ___ ___  ___| / ___| _ __ | (_)_ __   ___
 |  |  _ \ / _ \/ __/ __|/ _ \ \___ \| '_ \| | | '_ \ / _ \
 |  | |_) |  __/\__ \__ \  __/ |___) | |_) | | | | | |  __/
 |  |____/ \___||___/___/\___|_|____/| .__/|_|_|_| |_|\___|
 |                                   |_|
\*/

namespace Splines {

  using namespace std; // load standard namspace

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Bessel_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts
  ) {

    size_t n = size_t(npts > 0 ? npts-1 : 0);
    std::vector<real_type> m(n+1);

    // calcolo slopes
    for ( size_t i = 0; i < n; ++i )
      m[i] = (Y[i+1]-Y[i])/(X[i+1]-X[i]);

    if ( npts == 2 ) { // caso speciale 2 soli punti

      Yp[0] = Yp[1] = m[0];

    } else {

      for ( size_t i = 1; i < n; ++i ) {
        real_type DL = X[i]   - X[i-1];
        real_type DR = X[i+1] - X[i];
        Yp[i] = (DR*m[i-1]+DL*m[i])/((DL+DR));
      }

      Yp[0] = 1.5*m[0]-0.5*m[1];
      Yp[n] = 1.5*m[n-1]-0.5*m[n-2];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BesselSpline::build (void) {
    SPLINE_ASSERT(
      this->npts > 1,
      "BesselSpline::build(): npts = " << npts << " not enought points"
    )
    integer ibegin = 0;
    integer iend   = 0;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < this->npts && this->X[iend-1] < this->X[iend] ) {}
      Bessel_build( this->X+ibegin, this->Y+ibegin, this->Yp+ibegin, iend-ibegin );
      ibegin = iend;
    } while ( iend < this->npts );
    
    SPLINE_CHECK_NAN( this->Yp, "BesselSpline::build(): Yp", this->npts );
  }
}
