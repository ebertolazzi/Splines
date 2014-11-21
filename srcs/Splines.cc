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

  sizeType
  Spline::search( valueType x ) const {
    if ( X[lastInterval] < x || X[lastInterval+1] > x ) {
      lastInterval = sizeType(lower_bound( X.begin(), X.end(), x ) - X.begin()) ;
      if ( lastInterval > 0 ) --lastInterval ;
    }
    return lastInterval ;
  }

  sizeType
  SplineSurf::search_x( valueType x ) const {
    if ( X[lastInterval_x] < x || X[lastInterval_x+1] > x ) {
      lastInterval_x = sizeType(lower_bound( X.begin(), X.end(), x ) - X.begin()) ;
      if ( lastInterval_x > 0 ) --lastInterval_x ;
    }
    return lastInterval_x ;
  }

  sizeType
  SplineSurf::search_y( valueType y ) const {
    if ( Y[lastInterval_y] < y || Y[lastInterval_y+1] > y ) {
      lastInterval_y = sizeType(lower_bound( Y.begin(), Y.end(), y ) - Y.begin()) ;
      if ( lastInterval_y > 0 ) --lastInterval_y ;
    }
    return lastInterval_y ;
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
    X.push_back( x ) ;
    Y.push_back( y ) ;
    ++npts ;
  }

  ///////////////////////////////////////////////////////////////////////////

  void
  Spline::dropBack() {
    if ( npts > 0 ) {
      --npts ;
      X.pop_back() ;
      Y.pop_back() ;
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
  ConstantsSpline::writeToStream( std::basic_ostream<char> & s ) const {
    sizeType nseg = sizeType(Y.size()) ;
    for ( sizeType i = 0 ; i < nseg ; ++i )
      s << "segment N." << setw(4) << i
        << " X:[ " << X[i] << ", " << X[i+1] << " ] Y:" << Y[i]
        << '\n' ; 
  }

  void
  LinearSpline::writeToStream( std::basic_ostream<char> & s ) const {
    sizeType nseg = sizeType(Y.size()-1) ;
    for ( sizeType i = 0 ; i < nseg ; ++i )
      s << "segment N." << setw(4) << i
        << " X:[ " << X[i] << ", " << X[i+1] << " ] Y:[ " << Y[i] << ", " << Y[i+1] 
        << " ] slope: " << (Y[i+1]-Y[i])/(X[i+1]-X[i])
        << '\n' ; 
  }
  
  // ---------------------------------------------------------------------------
  void
  SplineSurf::build ( valueType const x[], sizeType incx,
                      valueType const y[], sizeType incy,
                      valueType const z[], sizeType incz,
                      sizeType nx, sizeType ny,
                      bool transpose ) {
    X.resize(nx) ;
    Y.resize(ny) ;
    Z.resize(nx*ny) ;
    for ( sizeType i = 0 ; i < nx ; ++i ) X[i] = x[i*incx] ;
    for ( sizeType i = 0 ; i < ny ; ++i ) Y[i] = y[i*incy] ;
    if ( transpose ) {
      for ( sizeType i = 0 ; i < nx ; ++i )
        for ( sizeType j = 0 ; j < ny ; ++j )
          Z[ipos(i,j)] = z[ipos_transpose(i,j)*incz] ;
    } else {
      for ( sizeType i = 0 ; i < nx*ny ; ++i ) Z[i] = z[i*incz] ;
    }
    makeSpline() ;
  }

  void
  SplineSurf::build ( VectorOfValues const & x,
                      VectorOfValues const & y,
                      VectorOfValues const & z,
                      bool transpose ) {
    SPLINE_ASSERT( z.size() >= x.size()*y.size(),
                   "in BilinearSpline::build insufficient dimension of Z, expeced >= " << x.size()*y.size() <<
                   " found " << z.size() ) ;
    X.resize(x.size()) ;
    Y.resize(y.size()) ;
    Z.resize(x.size()*y.size()) ;
    std::copy( x.begin(), x.end(), X.begin() ) ;
    std::copy( y.begin(), y.end(), Y.begin() ) ;
    if ( transpose ) {
      for ( sizeType i = 0 ; i < x.size() ; ++i )
        for ( sizeType j = 0 ; j < y.size() ; ++j )
          Z[ipos(i,j)] = z[ipos_transpose(i,j)] ;
    } else {
      std::copy( z.begin(), z.begin()+Z.size(), Z.begin() ) ;
    }
    makeSpline() ;
  }

  void
  SplineSurf::build ( valueType const z[], sizeType incz,
                      sizeType nx, sizeType ny,
                      bool transpose ) {
    X.resize(nx) ;
    Y.resize(ny) ;
    Z.resize(nx*ny) ;
    for ( sizeType i = 0 ; i < nx    ; ++i ) X[i] = i ;
    for ( sizeType i = 0 ; i < ny    ; ++i ) Y[i] = i ;
    for ( sizeType i = 0 ; i < nx*ny ; ++i ) Z[i] = z[i*incz] ;
    if ( transpose ) {
      for ( sizeType i = 0 ; i < nx ; ++i )
        for ( sizeType j = 0 ; j < ny ; ++j )
          Z[ipos(i,j)] = z[ipos_transpose(i,j)*incz] ;
    } else {
      for ( sizeType i = 0 ; i < nx*ny ; ++i ) Z[i] = z[i*incz] ;
    }
    makeSpline() ;
  }

  void
  SplineSurf::build ( VectorOfValues const & z,
                      sizeType nx, sizeType ny,
                      bool transpose ) {
    SPLINE_ASSERT( z.size() >= nx*ny,
                   "in SplineSurf::build insufficient dimension of Z, expeced >= " << nx*ny <<
                   " found " << z.size() ) ;
    X.resize(nx) ;
    Y.resize(ny) ;
    Z.resize(nx*ny) ;
    for ( sizeType i = 0 ; i < nx ; ++i ) X[i] = i ;
    for ( sizeType i = 0 ; i < ny ; ++i ) Y[i] = i ;
    if ( transpose ) {
      for ( sizeType i = 0 ; i < nx ; ++i )
        for ( sizeType j = 0 ; j < ny ; ++j )
          Z[ipos(i,j)] = z[ipos_transpose(i,j)] ;
    } else {
      std::copy( z.begin(), z.begin()+nx*ny, Z.begin() ) ;
    }
    makeSpline() ;
  }

  void
  BiCubicSplineBase::load( sizeType i, sizeType j ) const {

    bili[0][0] = Z[ipos(i,j)] ;
    bili[0][1] = Z[ipos(i,j+1)] ;
    bili[1][0] = Z[ipos(i+1,j)];
    bili[1][1] = Z[ipos(i+1,j+1)];

    bili[2][0] = DX[ipos(i,j)] ;
    bili[2][1] = DX[ipos(i,j+1)] ;
    bili[3][0] = DX[ipos(i+1,j)];
    bili[3][1] = DX[ipos(i+1,j+1)];

    bili[0][2] = DY[ipos(i,j)] ;
    bili[0][3] = DY[ipos(i,j+1)] ;
    bili[1][2] = DY[ipos(i+1,j)] ;
    bili[1][3] = DY[ipos(i+1,j+1)];

    bili[2][2] = DXY[ipos(i,j)] ;
    bili[2][3] = DXY[ipos(i,j+1)] ;
    bili[3][2] = DXY[ipos(i+1,j)] ;
    bili[3][3] = DXY[ipos(i+1,j+1)] ;

  }

  valueType
  BiCubicSplineBase::operator () ( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3( x - X[i], X[i+1] - X[i], u ) ;
    Hermite3( y - Y[j], Y[j+1] - Y[j], v ) ;
    load(i,j) ;
    return bilinear( u, bili, v ) ;
  }

  valueType
  BiCubicSplineBase::Dx( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3_D( x - X[i], X[i+1] - X[i], u_D ) ;
    Hermite3  ( y - Y[j], Y[j+1] - Y[j], v   ) ;
    load(i,j) ;
    return bilinear( u_D, bili, v ) ;
  }

  valueType
  BiCubicSplineBase::Dy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3  ( x - X[i], X[i+1] - X[i], u   ) ;
    Hermite3_D( y - Y[j], Y[j+1] - Y[j], v_D ) ;
    load(i,j) ;
    return bilinear( u, bili, v_D ) ;
  }

  valueType
  BiCubicSplineBase::Dxy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3_D( x - X[i], X[i+1] - X[i], u_D ) ;
    Hermite3_D( y - Y[j], Y[j+1] - Y[j], v_D ) ;
    load(i,j) ;
    return bilinear( u_D, bili, v_D ) ;
  }

  valueType
  BiCubicSplineBase::Dxx( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3_DD( x - X[i], X[i+1] - X[i], u_DD ) ;
    Hermite3   ( y - Y[j], Y[j+1] - Y[j], v    ) ;
    load(i,j) ;
    return bilinear( u_DD, bili, v ) ;
  }

  valueType
  BiCubicSplineBase::Dyy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3   ( x - X[i], X[i+1] - X[i], u    ) ;
    Hermite3_DD( y - Y[j], Y[j+1] - Y[j], v_DD ) ;
    load(i,j) ;
    return bilinear( u, bili, v_DD ) ;
  }

  void
  BiCubicSplineBase::D( valueType x, valueType y, valueType d[3] ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3   ( x - X[i], X[i+1] - X[i], u    ) ;
    Hermite3_D ( x - X[i], X[i+1] - X[i], u_D  ) ;
    Hermite3   ( y - Y[j], Y[j+1] - Y[j], v    ) ;
    Hermite3_D ( y - Y[j], Y[j+1] - Y[j], v_D  ) ;
    load(i,j) ;
    d[0] = bilinear( u, bili, v ) ;
    d[1] = bilinear( u_D, bili, v ) ;
    d[2] = bilinear( u, bili, v_D ) ;
  }

  void
  BiCubicSplineBase::DD( valueType x, valueType y, valueType d[6] ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3   ( x - X[i], X[i+1] - X[i], u    ) ;
    Hermite3_D ( x - X[i], X[i+1] - X[i], u_D  ) ;
    Hermite3_DD( x - X[i], X[i+1] - X[i], u_DD ) ;
    Hermite3   ( y - Y[j], Y[j+1] - Y[j], v    ) ;
    Hermite3_D ( y - Y[j], Y[j+1] - Y[j], v_D  ) ;
    Hermite3_DD( y - Y[j], Y[j+1] - Y[j], v_DD ) ;
    load(i,j) ;
    d[0] = bilinear( u, bili, v ) ;
    d[1] = bilinear( u_D, bili, v ) ;
    d[2] = bilinear( u, bili, v_D ) ;
    d[3] = bilinear( u_DD, bili, v ) ;
    d[4] = bilinear( u_D, bili, v_D ) ;
    d[5] = bilinear( u, bili, v_DD ) ;
  }

}
