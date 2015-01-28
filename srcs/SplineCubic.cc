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

/*
//   ####  #    # #####  #  ####
//  #    # #    # #    # # #    #
//  #      #    # #####  # #
//  #      #    # #    # # #
//  #    # #    # #    # # #    #
//   ####   ####  #####  #  ####
*/

namespace Splines {

  using namespace std ; // load standard namspace

  static
  void
  CubicSpline_build( valueType const X[],
                     valueType const Y[],
                     valueType       Yp[],
                     sizeType        npts,
                     valueType       ddy0,
                     valueType       ddyn ) {

    sizeType n = npts-1 ;

    #ifdef SPLINE_USE_ALLOCA
    valueType * L = (valueType*)alloca( npts*sizeof(valueType) ) ;
    valueType * D = (valueType*)alloca( npts*sizeof(valueType) ) ;
    valueType * U = (valueType*)alloca( npts*sizeof(valueType) ) ;
    valueType * Z = (valueType*)alloca( npts*sizeof(valueType) ) ;
    #else
    valueType L[npts] ;
    valueType D[npts] ;
    valueType U[npts] ;
    valueType Z[npts] ;
    #endif


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

  }

  // ---------------------------------------------------------------------------

  void
  CubicSpline::build() {
    SPLINE_ASSERT( npts > 1,"CubicSpline::build(): npts = " << npts << " not enought points" );
    sizeType ibegin = 0 ;
    sizeType iend   = 0 ;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < npts && X[iend-1] < X[iend] ) {} ;
      valueType d0 = ibegin == 0    ? ddy0 : 0 ;
      valueType d1 = iend   == npts ? ddyn : 0 ;
      CubicSpline_build( X+ibegin, Y+ibegin, Yp+ibegin, iend-ibegin, d0, d1 ) ;
      ibegin = iend ;
    } while ( iend < npts ) ;
    
    SPLINE_CHECK_NAN(Yp,"CubicSpline::build(): Yp",npts);
  }
  
  /*
  //    ____  ____   ____                               _
  //   / ___|/ ___| / ___| _   _ _ __  _ __   ___  _ __| |_
  //  | |  _| |     \___ \| | | | '_ \| '_ \ / _ \| '__| __|
  //  | |_| | |___   ___) | |_| | |_) | |_) | (_) | |  | |_
  //   \____|\____| |____/ \__,_| .__/| .__/ \___/|_|   \__|
  //                            |_|   |_|
  */
  #ifdef SPLINES_USE_GENERIC_CONTAINER
  void
  CubicSpline::build( GC::GenericContainer const & gc ) {
    /*
    // gc["x"]
    // gc["y"]
    // gc["y''(begin)"]
    // gc["y''(end)"]
    //
    */
    SPLINE_ASSERT( gc.exists("x"), "[" << _name << "] ConstantSpline::build, missing `x` field!") ;
    SPLINE_ASSERT( gc.exists("y"), "[" << _name << "] ConstantSpline::build, missing `y` field!") ;

    valueType y_DD_0 = 0 ;
    if ( gc.exists("y''(begin)") ) y_DD_0 = gc("y''(begin)").get_number() ;

    valueType y_DD_1 = 0 ;
    if ( gc.exists("y''(end)") ) y_DD_1 = gc("y''(end)").get_number() ;
    
    setbc( y_DD_0, y_DD_1 ) ;
  
    GC::GenericContainer const & gc_x = gc("x") ;
    GC::GenericContainer const & gc_y = gc("y") ;

    SPLINE_ASSERT( GC::GC_VEC_REAL == gc_x.get_type(),
                   "Field `x` expected to be of type `vec_real_type` found: `" <<
                   gc_x.get_type_name() << "`" ) ;

    SPLINE_ASSERT( GC::GC_VEC_REAL == gc_y.get_type(),
                   "Field `y` expected to be of type `vec_real_type` found: `" <<
                   gc_y.get_type_name() << "`" ) ;

    build( gc_x.get_vec_real(), gc_y.get_vec_real() ) ;
  }
  #endif

}
