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
  
  // ---------------------------------------------------------------------------
  void
  SplineSurf::build ( valueType const x[], sizeType incx,
                      valueType const y[], sizeType incy,
                      valueType const z[], sizeType ldZ,
                      sizeType nx, sizeType ny,
                      bool fortran_storage,
                      bool transposed ) {
    X.resize(nx) ;
    Y.resize(ny) ;
    Z.resize(nx*ny) ;
    for ( sizeType i = 0 ; i < nx ; ++i ) X[i] = x[i*incx] ;
    for ( sizeType i = 0 ; i < ny ; ++i ) Y[i] = y[i*incy] ;
    if ( fortran_storage ) {
      if ( transposed ) {
        for ( sizeType i = 0 ; i < nx ; ++i )
          for ( sizeType j = 0 ; j < ny ; ++j )
            Z[ipos_C(i,j)] = z[ipos_F(j,i,ldZ)] ;
      } else {
        for ( sizeType i = 0 ; i < nx ; ++i )
          for ( sizeType j = 0 ; j < ny ; ++j )
            Z[ipos_C(i,j)] = z[ipos_F(i,j,ldZ)] ;
      }
    } else {
      if ( transposed ) {
        for ( sizeType i = 0 ; i < nx ; ++i )
          for ( sizeType j = 0 ; j < ny ; ++j )
            Z[ipos_C(i,j)] = z[ipos_C(j,i,ldZ)] ;
      } else {
        for ( sizeType i = 0 ; i < nx ; ++i )
          for ( sizeType j = 0 ; j < ny ; ++j )
            Z[ipos_C(i,j)] = z[ipos_C(i,j,ldZ)] ;
      }
    }
    Z_max = *std::max_element(Z.begin(),Z.end()) ;
    Z_min = *std::min_element(Z.begin(),Z.end()) ;
    makeSpline() ;
  }

  void
  SplineSurf::build( valueType const z[], sizeType ldZ,
                     sizeType nx, sizeType ny,
                     bool fortran_storage,
                     bool transposed ) {
    VectorOfValues XX(nx), YY(ny) ; // temporary vector
    for ( sizeType i = 0 ; i < nx ; ++i ) XX[i] = i ;
    for ( sizeType i = 0 ; i < ny ; ++i ) YY[i] = i ;
    build ( &XX.front(), 1, &YY.front(), 1, z, ldZ, nx, ny, fortran_storage, transposed ) ;
  }

  void
  BiCubicSplineBase::load( sizeType i, sizeType j ) const {

    bili[0][0] = Z[ipos_C(i,j)] ;
    bili[0][1] = Z[ipos_C(i,j+1)] ;
    bili[1][0] = Z[ipos_C(i+1,j)];
    bili[1][1] = Z[ipos_C(i+1,j+1)];

    bili[2][0] = DX[ipos_C(i,j)] ;
    bili[2][1] = DX[ipos_C(i,j+1)] ;
    bili[3][0] = DX[ipos_C(i+1,j)];
    bili[3][1] = DX[ipos_C(i+1,j+1)];

    bili[0][2] = DY[ipos_C(i,j)] ;
    bili[0][3] = DY[ipos_C(i,j+1)] ;
    bili[1][2] = DY[ipos_C(i+1,j)] ;
    bili[1][3] = DY[ipos_C(i+1,j+1)];

    bili[2][2] = DXY[ipos_C(i,j)] ;
    bili[2][3] = DXY[ipos_C(i,j+1)] ;
    bili[3][2] = DXY[ipos_C(i+1,j)] ;
    bili[3][3] = DXY[ipos_C(i+1,j+1)] ;

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
