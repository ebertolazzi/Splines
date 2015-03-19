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
#include <iomanip>

/**
 * 
 */

#ifndef isnan
#define isnan(A) ((A) != (A))
#endif

#ifndef isfinite
#define isfinite(A) ((A*0.0) == 0.0)
#endif

#ifndef isregular
#define isregular(A) ( !isnan(A) && isfinite(A) )
#endif

namespace Splines {
  //! Associate a number for each type of splines implemented
  
  char const * spline_type[] = {
    "constant",
    "linear",
    "Akima",
    "Bessel",
    "Pchip",
    "Cubic",
    "Quintic"
  } ;

  /*       _               _    _   _       _   _
  //   ___| |__   ___  ___| | _| \ | | __ _| \ | |
  //  / __| '_ \ / _ \/ __| |/ /  \| |/ _` |  \| |
  // | (__| | | |  __/ (__|   <| |\  | (_| | |\  |
  //  \___|_| |_|\___|\___|_|\_\_| \_|\__,_|_| \_|
  */
  //! check if the vector `pv` os size `DIM` contains only regular floats. If not an error is issued
  void
  checkNaN( valueType const pv[],
            char      const v_name[],
            sizeType  const DIM ) {
    for ( sizeType i = 0 ; i < DIM ; ++i ) {
      if ( isnan(pv[i]) ) {
        std::ostringstream ost ;
        ost << "\nfound NaN at " << v_name << "[" << i << "]" ;
        throw std::runtime_error(ost.str()) ;
      }
      if ( !isfinite(pv[i]) ) {
        std::ostringstream ost ;
        ost << "\nfound Infinity at " << v_name << "[" << i << "]" ;
        throw std::runtime_error(ost.str()) ;
      }
    }
  }
  
  void
  Spline::pushBack( valueType x, valueType y ) {
    if ( npts > 0 ) {
      SPLINE_ASSERT( x >= X[npts-1], // ammetto punti doppi
                     "Spline::pushBack, non monotone insert at insert N. " << npts <<
                     "\nX[ " << npts-1 << "] = " << X[npts-1] <<
                     "\nX[ " << npts   << "] = " << x ) ;
    }
    if ( npts_reserved == 0 ) {
      reserve( 2 ) ;
    } else if ( npts >= npts_reserved ) {
      // riallocazione & copia
      sizeType saved_npts = npts ; // salvo npts perche reserve lo azzera
      #ifdef SPLINE_USE_ALLOCA
        valueType * Xsaved = (valueType*)alloca( (npts)*sizeof(valueType) ) ;
        valueType * Ysaved = (valueType*)alloca( (npts)*sizeof(valueType) ) ;
      #else
        valueType Xsaved[npts], Ysaved[npts] ;
      #endif

      
      std::copy( X, X+npts, Xsaved ) ;
      std::copy( Y, Y+npts, Ysaved ) ;
      reserve( (npts+1) * 2 ) ;
      npts = saved_npts ;
      std::copy( Xsaved, Xsaved+npts, X ) ;
      std::copy( Ysaved, Ysaved+npts, Y ) ;
    }
    X[npts] = x ;
    Y[npts] = y ;
    ++npts ;
  }

  sizeType
  Spline::search( valueType x ) const {
    SPLINE_ASSERT( npts > 0, "\nsearch(" << x << ") empty spline");
    if ( X[lastInterval] < x || X[lastInterval+1] > x ) {
      if ( _check_range ) {
        SPLINE_ASSERT( x >= X[0] && x <= X[npts-1],
                       "method search( " << x << " ) out of range: [" <<
                       X[0] << ", " << X[npts-1] << "]" ) ;
      }
      lastInterval = sizeType(lower_bound( X, X+npts, x ) - X) ;
      if ( lastInterval > 0 ) --lastInterval ;
      if ( X[lastInterval] == X[lastInterval+1] ) ++lastInterval ; // degenerate interval for duplicated nodes
    }
    return lastInterval ;
  }

  ///////////////////////////////////////////////////////////////////////////
  void
  Spline::setOrigin( valueType x0 ) {
    valueType Tx = x0 - X[0] ;
    valueType *ix = X ;
    while ( ix < X+npts ) *ix++ += Tx ;
  }

  void
  Spline::setRange( valueType xmin, valueType xmax ) {
    SPLINE_ASSERT( xmax > xmin, "Spline::setRange( " << xmin << " , " << xmax << " ) bad range ") ;
    valueType S  = (xmax - xmin) / ( X[npts-1] - X[0] ) ;
    valueType Tx = xmin - S * X[0] ;
    for( valueType *ix = X ; ix < X+npts ; ++ix ) *ix = *ix * S + Tx ;
  }

  ///////////////////////////////////////////////////////////////////////////
  void
  Spline::dump( ostream & s, sizeType nintervals, char const header[] ) const {
    s << header << '\n' ;
    valueType dx = (xMax()-xMin())/nintervals ;
    for ( sizeType i = 0 ; i <= nintervals ; ++i ) {
      valueType x = xMin() + i*dx ;
      s << x << '\t' << (*this)(x) << '\n' ;
    }
  }

}
