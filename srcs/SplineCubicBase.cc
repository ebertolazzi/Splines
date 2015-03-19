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

  void
  CubicSplineBase::build ( valueType const x[], sizeType incx,
                           valueType const y[], sizeType incy,
                           sizeType n ) {
    reserve( n ) ;
    for ( sizeType i = 0 ; i < n ; ++i ) X[i] = x[i*incx] ;
    for ( sizeType i = 0 ; i < n ; ++i ) Y[i] = y[i*incy] ;
    npts = n ;
    build() ;
  }

  void
  CubicSplineBase::build ( valueType const x[],
                           valueType const y[],
                           sizeType n ) {
    reserve( n ) ;
    std::copy( x, x+n, X );
    std::copy( y, y+n, Y );
    npts = n ;
    build() ;
  }

  void
  CubicSplineBase::build ( vector<valueType> const & x, vector<valueType> const & y ) {
    sizeType n = sizeType(min( x.size(), y.size() )) ;
    reserve( n ) ;
    std::copy( x.begin(), x.begin()+n, X );
    std::copy( y.begin(), y.begin()+n, Y );
    npts = n ;
    build() ;
  }
  
  void
  CubicSplineBase::clear(void) {
    if ( !_external_alloc ) baseValue.free() ;
    npts = npts_reserved = 0 ;
    _external_alloc = false ;
    X = Y = Yp = nullptr ;
  }

  void
  CubicSplineBase::reserve( sizeType n ) {
    if ( _external_alloc && n <= npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      npts_reserved = n ;
      baseValue.allocate(3*n) ;
      X  = baseValue(n) ;
      Y  = baseValue(n) ;
      Yp = baseValue(n) ;
      _external_alloc = false ;
    }
    npts = lastInterval = 0 ;
  }

  void
  CubicSplineBase::reserve_external( sizeType n, valueType *& p_x, valueType *& p_y, valueType *& p_dy ) {
    npts_reserved = n ;
    X    = p_x ;
    Y    = p_y ;
    Yp   = p_dy ;
    npts = lastInterval = 0 ;
    _external_alloc = true ;
  }

  valueType
  CubicSplineBase::operator () ( valueType x ) const { 
    sizeType i = Spline::search( x ) ;
    Hermite3( x-X[i], X[i+1]-X[i], base ) ;
    return base[0] * Y[i]   +
           base[1] * Y[i+1] +
           base[2] * Yp[i]  +
           base[3] * Yp[i+1] ;
  }

  valueType
  CubicSplineBase::D( valueType x ) const {
    sizeType i = Spline::search( x ) ;
    Hermite3_D( x-X[i], X[i+1]-X[i], base_D ) ;
    return base_D[0] * Y[i]   +
           base_D[1] * Y[i+1] +
           base_D[2] * Yp[i]  +
           base_D[3] * Yp[i+1] ;
  }

  valueType
  CubicSplineBase::DD( valueType x ) const {
    sizeType i = Spline::search( x ) ;
    Hermite3_DD( x-X[i], X[i+1]-X[i], base_DD ) ;
    return base_DD[0] * Y[i]   +
           base_DD[1] * Y[i+1] +
           base_DD[2] * Yp[i]  +
           base_DD[3] * Yp[i+1] ;
  }

  valueType
  CubicSplineBase::DDD( valueType x ) const { 
    sizeType i = Spline::search( x ) ;
    Hermite3_DDD( x-X[i], X[i+1]-X[i], base_DDD ) ;
    return base_DDD[0] * Y[i]   +
           base_DDD[1] * Y[i+1] +
           base_DDD[2] * Yp[i]  +
           base_DDD[3] * Yp[i+1] ;
  }

  // Implementation
  void
  CubicSplineBase::copySpline( CubicSplineBase const & S ) {
    CubicSplineBase::reserve(S.npts) ;
    npts = S.npts ;
    std::copy( S.X, S.X+npts, X ) ;
    std::copy( S.Y, S.Y+npts, Y ) ;
    std::copy( S.Yp, S.Yp+npts, Yp ) ;
  }

  //! change X-range of the spline
  void
  CubicSplineBase::setRange( valueType xmin, valueType xmax ) {
    Spline::setRange( xmin, xmax ) ;
    valueType recS = ( X[npts-1] - X[0] ) / (xmax - xmin) ;
    valueType * iy = Y ;
    while ( iy < Y + npts ) *iy++ *= recS ;
  }

  void
  CubicSplineBase::writeToStream( std::basic_ostream<char> & s ) const {
    sizeType nseg = npts-1 ;
    for ( sizeType i = 0 ; i < nseg ; ++i )
      s << "segment N." << setw(4) << i
        << " X:[" << X[i] << ", " << X[i+1]
        << "] Y:[" << Y[i] << ", " << Y[i+1] 
        << "] Yp:[" << Yp[i] << ", " << Yp[i+1] 
        << "] slope: " << (Y[i+1]-Y[i])/(X[i+1]-X[i])
        << '\n' ; 
  }

}
