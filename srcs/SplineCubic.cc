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
  //   ####  #    # #####  #  ####  
  //  #    # #    # #    # # #    # 
  //  #      #    # #####  # #      
  //  #      #    # #    # # #      
  //  #    # #    # #    # # #    # 
  //   ####   ####  #####  #  ####  
  */
  void
  CubicSpline::build() {
    sizeType n = npts-1 ;
    Yp . resize(npts) ;

    VectorOfValues L(npts), D(npts), U(npts), Z(npts) ;

    sizeType i ;
    for ( i = 1 ; i < n ; ++i ) {
      L[i] = (X[i]-X[i-1])/(X[i+1]-X[i-1]) ;
      U[i] = (X[i+1]-X[i])/(X[i+1]-X[i-1]) ;
      D[i] = 2 ;
      Z[i] = 6 * ( (Y[i+1]-Y[i])/(X[i+1]-X[i]) -
                   (Y[i]-Y[i-1])/(X[i]-X[i-1]) ) / ( X[i+1] - X[i-1] ) ;
    }

    L[0] = 0 ; D[0] = 1 ; U[0] = 0 ; Z[0] = ddy0 ;
    L[n] = 0 ; D[n] = 1 ; U[n] = 0 ; Z[n] = ddyn ;

    if ( n > 2 ) {
      i = 0 ;
      do {
        Z[i]   /= D[i] ;
        U[i]   /= D[i] ;
        D[i+1] -= L[i+1] * U[i] ;
        Z[i+1] -= L[i+1] * Z[i] ;
      } while ( ++i < n ) ;

      Z[i] /= D[i] ;

      do {
        --i ;
        Z[i] -= U[i] * Z[i+1] ;
      } while ( i > 0 ) ;
    }

    for ( i = 0 ; i < n ; ++i ) {
      valueType DX = X[i+1] - X[i] ;
      Yp[i] = (Y[i+1]-Y[i])/DX - (Z[i]/3 + Z[i+1]/6) * DX ;
    }
    valueType DX = X[n] - X[n-1] ;
    Yp[n] = Yp[n-1] + DX * 0.5*(Z[n-1] + Z[n])  ;
    SPLINE_CHECK_NAN(&Yp.front(),"CubicSpline::cubic(): Yp",npts);
  }

}
