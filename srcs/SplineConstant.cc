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
  ConstantSpline::reserve_external( sizeType n, valueType *& p_x, valueType *& p_y ) {
    if ( !_external_alloc ) baseValue.free() ;
    npts            = 0 ;
    npts_reserved   = n ;
    _external_alloc = true ;
    X = p_x ;
    Y = p_y ;
  }

  void
  ConstantSpline::reserve( sizeType n ) {
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

  //! Evalute spline value at `x`
  valueType
  ConstantSpline::operator () ( valueType x ) const {
    if ( x < X[0]      ) return Y[0] ;
    if ( x > X[npts-1] ) return Y[npts-1] ;
    return Y[search(x)] ;
  }

  void
  ConstantSpline::build ( valueType const x[], sizeType incx,
                          valueType const y[], sizeType incy,
                          sizeType n ) {
    reserve( n ) ;
    for ( sizeType i = 0 ; i < n   ; ++i ) X[i] = x[i*incx] ;
    for ( sizeType i = 0 ; i < n-1 ; ++i ) Y[i] = y[i*incy] ;
    npts = n ;
    build() ;
  }

  void
  ConstantSpline::build ( valueType const x[], valueType const y[], sizeType n ) {
    reserve( n ) ;
    for ( sizeType i = 0 ; i < n   ; ++i ) X[i] = x[i] ;
    for ( sizeType i = 0 ; i < n-1 ; ++i ) Y[i] = y[i] ; // ultimo y ignorato
    npts = n ;
    build() ;
  }

  void
  ConstantSpline::build ( vector<valueType> const & x, vector<valueType> const & y ) {
    sizeType n = sizeType(min( x.size(), y.size()-1 )) ;
    reserve( n ) ;
    std::copy( x.begin(), x.begin()+n,   X );
    std::copy( y.begin(), y.begin()+n-1, Y );
    npts = n ;
    build() ;
  }

  void
  ConstantSpline::clear() {
    if ( !_external_alloc ) baseValue.free() ;
    npts = npts_reserved = 0 ;
    _external_alloc = false ;
    X = Y = nullptr ;
  }

  void
  ConstantSpline::writeToStream( std::basic_ostream<char> & s ) const {
    sizeType nseg = npts - 1 ;
    for ( sizeType i = 0 ; i < nseg ; ++i )
      s << "segment N." << setw(4) << i
        << " X:[ " << X[i] << ", " << X[i+1] << " ] Y:" << Y[i]
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
  ConstantSpline::build( GC::GenericContainer const & gc ) {
    /*
    // gc["x"]
    // gc["y"]
    //
    */
    SPLINE_ASSERT( gc.exists("x"), "[" << _name << "] ConstantSpline::build, missing `x` field!") ;
    SPLINE_ASSERT( gc.exists("y"), "[" << _name << "] ConstantSpline::build, missing `y` field!") ;
  
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
