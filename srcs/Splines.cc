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
    for ( indexType i = 0 ; i < DIM ; ++i ) {
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
      //// DA RISCRIVERE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if ( std::abs(x-X.back()) < 1e-9 ) return ; // workaround per punti doppi
      SPLINE_ASSERT( x > X.back(),
                     "Spline::pushBack, non monotone insert at insert N. " << npts <<
                     "\nX[ " << npts-1 << "] = " << X.back() <<
                     "\nX[ " << npts   << "] = " << x ) ;
    }
    X . push_back( x ) ;
    Y . push_back( y ) ;
    ++npts ;
  }

  ///////////////////////////////////////////////////////////////////////////

  void
  Spline::dropBack() {
    if ( npts > 0 ) {
      --npts ;
      X . pop_back() ;
      Y . pop_back() ;
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  void
  Spline::setOrigin( valueType x0 ) {
    valueType Tx = x0 - X.front() ;
    for ( VectorOfValues::iterator ix = X.begin() ; ix != X.end() ; ++ix ) *ix += Tx ;
  }

  void
  Spline::setRange( valueType xmin, valueType xmax ) {
    SPLINE_ASSERT( xmax > xmin, "Spline::setRange( " << xmin << " , " << xmax << " ) bad range ") ;
    valueType S  = (xmax - xmin) / ( X.back() - X.front() ) ;
    valueType Tx = xmin - S * X.front() ;
    for ( VectorOfValues::iterator ix = X.begin() ; ix != X.end() ; ++ix ) *ix = *ix * S + Tx ;
  }

  //! change X-range of the spline
  void
  CubicSplineBase::setRange( valueType xmin, valueType xmax ) {
    Spline::setRange( xmin, xmax ) ;
    valueType recS = ( X.back() - X.front() ) / (xmax - xmin) ;
    for ( VectorOfValues::iterator iy = Y.begin() ; iy != Y.end() ; ++iy ) *iy *= recS ;
  }

  void
  ConstantsSpline::writeToStream( ostream & s ) const {
    indexType nseg = indexType(Y.size()) ;
    for ( indexType i = 0 ; i < nseg ; ++i )
      s << "segment N." << setw(4) << i
        << " X:[ " << X[i] << ", " << X[i+1] << " ] Y:" << Y[i]
        << '\n' ; 
  }

  void
  LinearSpline::writeToStream( ostream & s ) const {
    indexType nseg = indexType(Y.size()-1) ;
    for ( indexType i = 0 ; i < nseg ; ++i )
      s << "segment N." << setw(4) << i
        << " X:[ " << X[i] << ", " << X[i+1] << " ] Y:[ " << Y[i] << ", " << Y[i+1] 
        << " ] slope: " << (Y[i+1]-Y[i])/(X[i+1]-X[i])
        << '\n' ; 
  }

}
