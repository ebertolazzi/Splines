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
#include <iomanip>

/**
 * 
 */

namespace Splines {

  using namespace std ; // load standard namspace

  valueType
  QuinticSplineBase::operator () ( valueType x ) const { 
    SPLINE_ASSERT( X.size() > 0 , "QuinticSplineBase(" << x << ") empty spline") ;
    if ( x < X.front() ) return Y.front() + ( x - X.front() ) * ( Yp.front() + ( x - X.front() ) * Ypp.front()/2 ) ;
    if ( x > X.back()  ) return Y.back()  + ( x - X.back() )  * ( Yp.back()  + ( x - X.back()  ) * Ypp.back()/2  ) ;
    sizeType i = Spline::search( x ) ;
    Hermite5( x-X[i], X[i+1]-X[i], base ) ;
    return base[0] * Y[i]    +
           base[1] * Y[i+1]  +
           base[2] * Yp[i]   +
           base[3] * Yp[i+1] +
           base[4] * Ypp[i]  +
           base[5] * Ypp[i+1] ;
  }

  valueType
  QuinticSplineBase::D( valueType x ) const { 
    SPLINE_ASSERT( X.size() > 0 , "QuinticSplineBase.D(" << x << ") empty spline");
    if ( x < X.front() ) return Yp.front() + ( x - X.front() ) * Ypp.front() ;
    if ( x > X.back()  ) return Yp.back()  + ( x - X.back()  ) * Ypp.back()  ;
    sizeType i = Spline::search( x ) ;
    Hermite5_D( x-X[i], X[i+1]-X[i], base_D ) ;
    return base_D[0] * Y[i]    +
           base_D[1] * Y[i+1]  +
           base_D[2] * Yp[i]   +
           base_D[3] * Yp[i+1] +
           base_D[4] * Ypp[i]  +
           base_D[5] * Ypp[i+1] ;
  }

  valueType
  QuinticSplineBase::DD( valueType x ) const { 
    SPLINE_ASSERT( X.size() > 0 , "QuinticSplineBase.DD(" << x << ") empty spline");
    if ( x < X.front() ) return Ypp.front() ;
    if ( x > X.back()  ) return Ypp.back()  ;
    sizeType i = Spline::search( x ) ;
    Hermite5_DD( x-X[i], X[i+1]-X[i], base_DD ) ;
    return base_DD[0] * Y[i]    +
           base_DD[1] * Y[i+1]  +
           base_DD[2] * Yp[i]   +
           base_DD[3] * Yp[i+1] +
           base_DD[4] * Ypp[i]  +
           base_DD[5] * Ypp[i+1] ;
  }

  valueType
  QuinticSplineBase::DDD( valueType x ) const { 
    SPLINE_ASSERT( X.size() > 0 , "QuinticSplineBase.DDD(" << x << ") empty spline");
    if ( x < X.front() ) return 0 ;
    if ( x > X.back()  ) return 0 ;
    sizeType i = Spline::search( x ) ;
    Hermite5_DDD( x-X[i], X[i+1]-X[i], base_DDD ) ;
    return base_DDD[0] * Y[i]    +
           base_DDD[1] * Y[i+1]  +
           base_DDD[2] * Yp[i]   +
           base_DDD[3] * Yp[i+1] +
           base_DDD[4] * Ypp[i]  +
           base_DDD[5] * Ypp[i+1] ;
  }

  valueType
  QuinticSplineBase::DDDD( valueType x ) const { 
    SPLINE_ASSERT( X.size() > 0 , "QuinticSplineBase.DDDD(" << x << ") empty spline");
    if ( x < X.front() ) return 0 ;
    if ( x > X.back()  ) return 0 ;
    sizeType i = Spline::search( x ) ;
    Hermite5_DDDD( x-X[i], X[i+1]-X[i], base_DDDD ) ;
    return base_DDDD[0] * Y[i]    +
           base_DDDD[1] * Y[i+1]  +
           base_DDDD[2] * Yp[i]   +
           base_DDDD[3] * Yp[i+1] +
           base_DDDD[4] * Ypp[i]  +
           base_DDDD[5] * Ypp[i+1] ;
  }

  valueType
  QuinticSplineBase::DDDDD( valueType x ) const { 
    SPLINE_ASSERT( X.size() > 0 , "QuinticSplineBase.DDDDD(" << x << ") empty spline");
    if ( x < X.front() ) return 0 ;
    if ( x > X.back()  ) return 0 ;
    sizeType i = Spline::search( x ) ;
    Hermite5_DDDDD( x-X[i], X[i+1]-X[i], base_DDDDD ) ;
    return base_DDDDD[0] * Y[i]    +
           base_DDDDD[1] * Y[i+1]  +
           base_DDDDD[2] * Yp[i]   +
           base_DDDDD[3] * Yp[i+1] +
           base_DDDDD[4] * Ypp[i]  +
           base_DDDDD[5] * Ypp[i+1] ;
  }

  // Implementation
  void
  QuinticSplineBase::copySpline( QuinticSplineBase const & S ) {
    npts = S . npts ;
    X   . resize( S . X   . size() ) ;
    Y   . resize( S . Y   . size() ) ;
    Yp  . resize( S . Yp  . size() ) ;
    Ypp . resize( S . Ypp . size() ) ;
    std::copy( S . X   . begin(), S . X   . end(), X   . begin() ) ;
    std::copy( S . Y   . begin(), S . Y   . end(), Y   . begin() ) ;
    std::copy( S . Yp  . begin(), S . Yp  . end(), Yp  . begin() ) ;
    std::copy( S . Ypp . begin(), S . Ypp . end(), Ypp . begin() ) ;
  }

  void
  QuinticSplineBase::allocate( valueType const x[], valueType const y[], sizeType n ) {
    X . clear() ; X . reserve(n) ;
    Y . clear() ; Y . reserve(n) ;
    npts = lastInterval = 0 ;
    for ( sizeType i = 0 ; i < n ; ++i ) pushBack( x[i], y[i] ) ;
    SPLINE_ASSERT ( npts > 1, "QuinticSplineBase::allocate, not enought point to define a spline\n" ) ;
  }

  void
  QuinticSplineBase::writeToStream( std::basic_ostream<char> & s ) const {
    sizeType nseg = sizeType(Y.size()-1) ;
    for ( sizeType i = 0 ; i < nseg ; ++i )
      s << "segment N." << setw(4) << i
        << " X:[" << X[i] << ", " << X[i+1]
        << "] Y:[" << Y[i] << ", " << Y[i+1] 
        << "] Yp:[" << Yp[i] << ", " << Yp[i+1]
        << "] Ypp:[" << Ypp[i] << ", " << Ypp[i+1]
        << "] slope: " << (Y[i+1]-Y[i])/(X[i+1]-X[i])
        << '\n' ; 
  }

}
