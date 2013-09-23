/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 1998                                                      |
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

/**
 * 
 */

namespace Splines {

  using namespace std ; // load standard namspace

  /*
  //   ____                     _ ____        _ _            
  //  | __ )  ___  ___ ___  ___| / ___| _ __ | (_)_ __   ___ 
  //  |  _ \ / _ \/ __/ __|/ _ \ \___ \| '_ \| | | '_ \ / _ \
  //  | |_) |  __/\__ \__ \  __/ |___) | |_) | | | | | |  __/
  //  |____/ \___||___/___/\___|_|____/| .__/|_|_|_| |_|\___|
  //                                   |_|                   
  */
  void
  BesselSpline::build (void) {
    sizeType n = npts-1 ;
    Yp . resize(npts) ;

    VectorOfValues m(n) ;

    // calcolo slopes
    for ( sizeType i = 0 ; i < n ; ++i )
      m[i] = (Y[i+1]-Y[i])/(X[i+1]-X[i]) ;

    if ( npts == 2 ) { // caso speciale 2 soli punti

      Yp[0] = Yp[1] = m[0] ;

    } else {

      for ( sizeType i = 1 ; i < n ; ++i ) {
        valueType DL = X[i]   - X[i-1] ;
        valueType DR = X[i+1] - X[i]   ;
        Yp[i] = (DR*m[i-1]+DL*m[i])/((DL+DR)) ;
      }

      Yp[0] = 1.5*m[0]-0.5*m[1] ;
      Yp[n] = 1.5*m[n-1]-0.5*m[n-2] ;
    }
  }

}
