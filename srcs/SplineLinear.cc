/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 1998-2014                                                 |
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

#include <iomanip>

/**
 * 
 */

namespace Splines {

  using namespace std ; // load standard namspace

  //! Use externally allocated memory for `npts` points
  void
  LinearSpline::reserve_external( sizeType n, valueType *& p_x, valueType *& p_y ) {
    if ( !_external_alloc ) baseValue.free() ;
    npts            = 0 ;
    npts_reserved   = n ;
    _external_alloc = true ;
    X = p_x ;
    Y = p_y ;
  }

  void
  LinearSpline::reserve( sizeType n ) {
    if ( _external_alloc && n <= npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      baseValue.allocate( 2*n ) ;
      npts_reserved   = n ;
      _external_alloc = false ;
      X = baseValue(n) ;
      Y = baseValue(n) ;
    }
    npts = lastInterval = 0 ;
  }

  //void
  //LinearSpline::pushBack( valueType x, valueType y )

  void
  LinearSpline::build ( valueType const x[], sizeType incx,
                        valueType const y[], sizeType incy,
                        sizeType n ) {
    reserve( n ) ;
    for ( sizeType i = 0 ; i < n ; ++i ) X[i] = x[i*incx] ;
    for ( sizeType i = 0 ; i < n ; ++i ) Y[i] = y[i*incy] ;
    npts = n ;
    build() ;
  }

  void
  LinearSpline::build ( valueType const x[],
                        valueType const y[],
                        sizeType n ) {
    reserve( n ) ;
    std::copy( x, x+n, X );
    std::copy( y, y+n, Y );
    npts = n ;
    build() ;
  }

  void
  LinearSpline::build ( vector<valueType> const & x, vector<valueType> const & y ) {
    sizeType n = sizeType(min( x.size(), y.size() )) ;
    reserve( n ) ;
    std::copy( x.begin(), x.begin()+n, X );
    std::copy( y.begin(), y.begin()+n, Y );
    npts = n ;
    build() ;
  }
  
  void
  LinearSpline::clear(void) {
    if ( !_external_alloc ) baseValue.free() ;
    npts = npts_reserved = 0 ;
    _external_alloc = false ;
    X = Y = nullptr ;
  }

  void
  LinearSpline::writeToStream( std::basic_ostream<char> & s ) const {
    sizeType nseg = npts - 1 ;
    for ( sizeType i = 0 ; i < nseg ; ++i )
      s << "segment N." << setw(4) << i
        << " X:[ " << X[i] << ", " << X[i+1] << " ] Y:[ " << Y[i] << ", " << Y[i+1] 
        << " ] slope: " << (Y[i+1]-Y[i])/(X[i+1]-X[i])
        << '\n' ; 
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
  LinearSpline::build( GC::GenericContainer const & gc ) {
    SPLINE_ASSERT( false, "Not Yet Implemented!" ) ;
  }
  #endif

}
