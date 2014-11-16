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
  //    ___        _       _   _      ____        _ _            
  //   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___ 
  //  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
  //  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
  //   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
  //                                       |_|                   
  //  
  */
  static
  inline
  int
  signTest( valueType const a, valueType const b ) {
    int sa = 0 ;
    int sb = 0 ;
    if      ( a > 0 ) sa =  1 ;
    else if ( a < 0 ) sa = -1 ;
    if      ( b > 0 ) sb =  1 ;
    else if ( b < 0 ) sb = -1 ;    
    return sa*sb ;
  }

  void
  QuinticSpline::build() {
    indexType ierr = 0 ;
    Yp.resize(npts) ;
    Ypp.resize(npts) ;
    sizeType n = sizeType(npts - 1) ;

    // function definition is ok, go on.
    valueType h1    = X[1] - X[0] ;
    valueType del1  = (Y[1]-Y[0])/h1 ;
    valueType dsave = del1 ;

    // special case n=2 -- use linear interpolation.
    if ( npts == 2 ) {
      Yp[0]  = Yp[1]  = del1 ;
      Ypp[0] = Ypp[1] = 0 ;
      return ;
    }

    valueType h2   = X[2] - X[1] ;
    valueType del2 = (Y[2]-Y[1])/h2 ;

    // Set Yp[0] via non-centered three-point formula, adjusted to be shape-preserving.
    valueType hsum = h1 + h2 ;
    valueType w1   = (h1 + hsum)/hsum ;
    valueType w2   = -h1/hsum ;
    Yp[0] = w1*del1 + w2*del2 ;
    valueType dmin, dmax ;
    if ( signTest(Yp[0],del1) <= 0 ) {
      Yp[0] = 0 ;
    } else if ( signTest(del1,del2) < 0 ) {
      // NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
      dmax = 3*del1 ;
      if ( std::abs(Yp[0]) > std::abs(dmax) ) Yp[0] = dmax ;
    }

    // loop through interior points.
    for ( sizeType i = 1 ; i < n ; ++i ) {
      if ( i > 1 ) {
        h1   = h2 ;
        h2   = X[i+1] - X[i] ;
        hsum = h1 + h2 ;
        del1 = del2 ;
        del2 = (Y[i+1] - Y[i])/h2 ;
      }
      // set Yp[i]=0 unless data are strictly monotonic.
      Yp[i] = 0 ;
      // count number of changes in direction of monotonicity.
      switch ( signTest(del1,del2) ) {
      case -1:
        if ( del2 == 0 ) break ;
        if ( signTest(dsave,del2) < 0 ) ++ierr ;
        dsave = del2 ;
        break ;
      case 0:
        ++ierr ;
        dsave = del2 ;
        break ;
      case 1: // use brodlie modification of butland formula.
        w1    = (1+h1/hsum)/3 ;
        w2    = (1+h2/hsum)/3 ;
        dmax  = std::max( std::abs(del1), std::abs(del2) ) ;
        dmin  = std::min( std::abs(del1), std::abs(del2) ) ;
        valueType drat1 = del1/dmax ;
        valueType drat2 = del2/dmax ;
        Yp[i] = dmin/(w1*drat1 + w2*drat2) ;
        break ;
      }
    }
    // set Yp[n] via non-centered three-point formula, adjusted to be shape-preserving.
    w1 = -h2/hsum ;
    w2 = (h2 + hsum)/hsum ;
    Yp[n] = w1*del1 + w2*del2 ;
    if ( signTest(Yp[n],del2) <= 0 ) {
      Yp[n] = 0 ;
    } else if ( signTest(del1,del2) < 0 ) {
      // need do this check only if monotonicity switches.
      dmax = 3*del2 ;
      if ( abs(Yp[n]) > abs(dmax) ) Yp[n] = dmax ;
    }
    // cout << "ierr = " << ierr << '\n' ;
    
    // calcolo Ypp (per il momento tutto a 0!)
    Ypp[0] = 0 ;
    h1 = X[1] - X[0] ;
    for ( sizeType i = 1 ; i < n ; ++i, h1 = h2 ) {
      h2 = X[i+1] - X[i] ;
      Ypp[i] = 0*2*(h2*Y[i-1]+h1*Y[i+1]-(h1+h2)*Y[i])/((h1+h2)*h1*h2) ;
    }
    Ypp[n] = 0 ;
  }
}
