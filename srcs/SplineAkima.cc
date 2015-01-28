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
/*
//     #
//    # #   #    # # #    #   ##
//   #   #  #   #  # ##  ##  #  #
//  #     # ####   # # ## # #    #
//  ####### #  #   # #    # ######
//  #     # #   #  # #    # #    #
//  #     # #    # # #    # #    #
*/

namespace Splines {

  using namespace std ; // load standard namspace

  static
  valueType
  akima_one( valueType epsi,
             valueType di_m2,
             valueType di_m1,
             valueType di,
             valueType di_p1 ) {
    valueType wl  = std::abs(di_p1 - di) ;
    valueType wr  = std::abs(di_m1 - di_m2);
    valueType den = wl + wr ;
    if ( den < epsi ) { wl = wr = 0.5 ; den = 1 ; }
    valueType num = wl * di_m1 + wr * di;
    return num / den ;
  }

  static
  void
  Akima_build( valueType const X[],
               valueType const Y[],
               valueType       Yp[],
               sizeType        npts ) {

    if ( npts == 2 ) { // solo 2 punti, niente da fare
      Yp[0] = Yp[1] = (Y[1]-Y[0])/(X[1]-X[0]) ;
    } else {
      #ifdef SPLINE_USE_ALLOCA
      valueType * m = (valueType*)alloca( (npts+4)*sizeof(valueType) ) ;
      #else
      valueType m[npts+4] ;
      #endif
      // calcolo slopes
      for ( sizeType i = 1 ; i < npts ; ++i )
        m[i+1] = (Y[i]-Y[i-1])/(X[i]-X[i-1]) ;

      // extra slope at the boundary
      m[1]      = 2*m[2]-m[3] ;
      m[0]      = 2*m[1]-m[2] ;
      m[npts+1] = 2*m[npts]-m[npts-1] ;
      m[npts+2] = 2*m[npts+1]-m[npts+2] ;

      valueType epsi = 0 ;
      for ( sizeType i = 0 ; i < npts+2 ; ++i ) {
        valueType dm = std::abs(m[i+1]-m[i]) ;
        if ( dm > epsi ) epsi = dm ; 
      }
      epsi *= 1E-8 ;

      for ( sizeType i = 0 ; i < npts ; ++i )
        Yp[i] = akima_one( epsi, m[i], m[i+1], m[i+2], m[i+3]) ;
    }
  }

  void
  AkimaSpline::build() {
    SPLINE_ASSERT( npts > 1,"AkimaSpline::build(): npts = " << npts << " not enought points" );
    sizeType ibegin = 0 ;
    sizeType iend   = 0 ;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < npts && X[iend-1] < X[iend] ) {} ;
      Akima_build( X+ibegin, Y+ibegin, Yp+ibegin, iend-ibegin ) ;
      ibegin = iend ;
    } while ( iend < npts ) ;
    
    SPLINE_CHECK_NAN(Yp,"AkimaSpline::build(): Yp",npts);
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
  AkimaSpline::build( GC::GenericContainer const & gc ) {
    /*
    // gc["x"]
    // gc["y"]
    //
    */
    SPLINE_ASSERT( gc.exists("x"), "[" << _name << "] AkimaSpline::build, missing `x` field!") ;
    SPLINE_ASSERT( gc.exists("y"), "[" << _name << "] AkimaSpline::build, missing `y` field!") ;
  
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
