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
//   ____                     _ ____        _ _
//  | __ )  ___  ___ ___  ___| / ___| _ __ | (_)_ __   ___
//  |  _ \ / _ \/ __/ __|/ _ \ \___ \| '_ \| | | '_ \ / _ \
//  | |_) |  __/\__ \__ \  __/ |___) | |_) | | | | | |  __/
//  |____/ \___||___/___/\___|_|____/| .__/|_|_|_| |_|\___|
//                                   |_|
*/

namespace Splines {

  using namespace std ; // load standard namspace

  static
  void
  Bessel_build( valueType const X[],
                valueType const Y[],
                valueType       Yp[],
                sizeType        npts ) {
    sizeType n = npts-1 ;

    #ifdef SPLINE_USE_ALLOCA
    valueType * m = (valueType*)alloca( n*sizeof(valueType) ) ;
    #else
    valueType m[n] ;
    #endif

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


  void
  BesselSpline::build (void) {
    SPLINE_ASSERT( npts > 1,"BesselSpline::build(): npts = " << npts << " not enought points" );
    sizeType ibegin = 0 ;
    sizeType iend   = 0 ;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < npts && X[iend-1] < X[iend] ) {} ;
      Bessel_build( X+ibegin, Y+ibegin, Yp+ibegin, iend-ibegin ) ;
      ibegin = iend ;
    } while ( iend < npts ) ;
    
    SPLINE_CHECK_NAN(Yp,"BesselSpline::build(): Yp",npts);
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
  BesselSpline::build( GC::GenericContainer const & gc ) {
    /*
    // gc["x"]
    // gc["y"]
    //
    */
    SPLINE_ASSERT( gc.exists("x"), "[" << _name << "] BesselSpline::build, missing `x` field!") ;
    SPLINE_ASSERT( gc.exists("y"), "[" << _name << "] BesselSpline::build, missing `y` field!") ;
  
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
