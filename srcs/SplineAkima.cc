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
#include <cmath>

/**
 * 
 */

namespace Splines {

  using namespace std ; // load standard namspace

  /*
  //     #                           
  //    # #   #    # # #    #   ##   
  //   #   #  #   #  # ##  ##  #  #  
  //  #     # ####   # # ## # #    # 
  //  ####### #  #   # #    # ###### 
  //  #     # #   #  # #    # #    # 
  //  #     # #    # # #    # #    #
  */
  valueType
  AkimaSpline::akima_one( valueType epsi,
                          valueType di_m2,
                          valueType di_m1,
                          valueType di,
                          valueType di_p1 ) const {
    valueType wl  = std::abs(di_p1 - di) ;
    valueType wr  = std::abs(di_m1 - di_m2);
    valueType den = wl + wr ;
    if ( den < epsi ) { wl = wr = 0.5 ; den = 1 ; }
    valueType num = wl * di_m1 + wr * di;
    return num / den ;
  }

  void
  AkimaSpline::build() {
    sizeType n = npts-1 ;
    Yp.resize(npts) ;

    VectorOfValues m(npts+4) ;

    // calcolo slopes
    for ( sizeType i = 0 ; i < n ; ++i )
      m[i+2] = (Y[i+1]-Y[i])/(X[i+1]-X[i]) ;

    if ( npts == 2 ) { // solo 2 punti, niente da fare
      Yp[0] = Yp[1] = m[2] ;
    } else {
      // extra slope at the boundary
      m[1]   = 2*m[2]-m[3] ;
      m[0]   = 2*m[1]-m[2] ;
      m[n+2] = 2*m[n+1]-m[n] ;
      m[n+3] = 2*m[n+2]-m[n+1] ;

      valueType epsi = 0 ;
      for ( sizeType i = 0 ; i < n+3 ; ++i ) {
        valueType dm = std::abs(m[i+1]-m[i]) ;
        if ( dm > epsi ) epsi = dm ; 
      }
      epsi *= 1E-8 ;

      for ( sizeType i = 0 ; i < npts ; ++i )
        Yp[i] = akima_one( epsi, m[i], m[i+1], m[i+2], m[i+3]) ;
    }
    SPLINE_CHECK_NAN(&Yp.front(),"AkimaSpline::build(): Yp",npts);

  }

}
