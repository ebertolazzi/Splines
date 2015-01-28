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
      if ( _check_range ) {
        SPLINE_ASSERT( x >= X.front() && x <= X.back(),
                       "method search_x( " << x << " ) out of range: [" <<
                       X.front() << ", " << X.back() << "]" ) ;
      }
      lastInterval_x = sizeType(lower_bound( X.begin(), X.end(), x ) - X.begin()) ;
      if ( lastInterval_x > 0 ) --lastInterval_x ;
    }
    return lastInterval_x ;
  }

  sizeType
  SplineSurf::search_y( valueType y ) const {
    if ( Y[lastInterval_y] < y || Y[lastInterval_y+1] > y ) {
      if ( _check_range ) {
        SPLINE_ASSERT( y >= Y.front() && y <= Y.back(),
                       "method search_y( " << y << " ) out of range: [" <<
                       Y.front() << ", " << Y.back() << "]" ) ;
      }
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
    vector<valueType> XX(nx), YY(ny) ; // temporary vector
    for ( sizeType i = 0 ; i < nx ; ++i ) XX[i] = i ;
    for ( sizeType i = 0 ; i < ny ; ++i ) YY[i] = i ;
    build ( &XX.front(), 1, &YY.front(), 1, z, ldZ, nx, ny, fortran_storage, transposed ) ;
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void
  BiCubicSplineBase::load( sizeType i, sizeType j ) const {

    bili3[0][0] = Z[ipos_C(i,j)] ;
    bili3[0][1] = Z[ipos_C(i,j+1)] ;
    bili3[0][2] = DY[ipos_C(i,j)] ;
    bili3[0][3] = DY[ipos_C(i,j+1)] ;

    bili3[1][0] = Z[ipos_C(i+1,j)];
    bili3[1][1] = Z[ipos_C(i+1,j+1)];
    bili3[1][2] = DY[ipos_C(i+1,j)] ;
    bili3[1][3] = DY[ipos_C(i+1,j+1)];

    bili3[2][0] = DX[ipos_C(i,j)] ;
    bili3[2][1] = DX[ipos_C(i,j+1)] ;
    bili3[2][2] = DXY[ipos_C(i,j)] ;
    bili3[2][3] = DXY[ipos_C(i,j+1)] ;

    bili3[3][0] = DX[ipos_C(i+1,j)];
    bili3[3][1] = DX[ipos_C(i+1,j+1)];
    bili3[3][2] = DXY[ipos_C(i+1,j)] ;
    bili3[3][3] = DXY[ipos_C(i+1,j+1)] ;

  }

  valueType
  BiCubicSplineBase::operator () ( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3( x - X[i], X[i+1] - X[i], u ) ;
    Hermite3( y - Y[j], Y[j+1] - Y[j], v ) ;
    load(i,j) ;
    return bilinear3( u, bili3, v ) ;
  }

  valueType
  BiCubicSplineBase::Dx( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3_D( x - X[i], X[i+1] - X[i], u_D ) ;
    Hermite3  ( y - Y[j], Y[j+1] - Y[j], v   ) ;
    load(i,j) ;
    return bilinear3( u_D, bili3, v ) ;
  }

  valueType
  BiCubicSplineBase::Dy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3  ( x - X[i], X[i+1] - X[i], u   ) ;
    Hermite3_D( y - Y[j], Y[j+1] - Y[j], v_D ) ;
    load(i,j) ;
    return bilinear3( u, bili3, v_D ) ;
  }

  valueType
  BiCubicSplineBase::Dxy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3_D( x - X[i], X[i+1] - X[i], u_D ) ;
    Hermite3_D( y - Y[j], Y[j+1] - Y[j], v_D ) ;
    load(i,j) ;
    return bilinear3( u_D, bili3, v_D ) ;
  }

  valueType
  BiCubicSplineBase::Dxx( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3_DD( x - X[i], X[i+1] - X[i], u_DD ) ;
    Hermite3   ( y - Y[j], Y[j+1] - Y[j], v    ) ;
    load(i,j) ;
    return bilinear3( u_DD, bili3, v ) ;
  }

  valueType
  BiCubicSplineBase::Dyy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3   ( x - X[i], X[i+1] - X[i], u    ) ;
    Hermite3_DD( y - Y[j], Y[j+1] - Y[j], v_DD ) ;
    load(i,j) ;
    return bilinear3( u, bili3, v_DD ) ;
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
    d[0] = bilinear3( u, bili3, v ) ;
    d[1] = bilinear3( u_D, bili3, v ) ;
    d[2] = bilinear3( u, bili3, v_D ) ;
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
    d[0] = bilinear3( u, bili3, v ) ;
    d[1] = bilinear3( u_D, bili3, v ) ;
    d[2] = bilinear3( u, bili3, v_D ) ;
    d[3] = bilinear3( u_DD, bili3, v ) ;
    d[4] = bilinear3( u_D, bili3, v_D ) ;
    d[5] = bilinear3( u, bili3, v_DD ) ;
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void
  BiQuinticSplineBase::load( sizeType i, sizeType j ) const {
  
    // H0, H1, dH0, dH1, ddH0, ddH1
    
    // + + + + + +
    // + + + + + +
    // + + + + . .
    // + + + + . .
    // + + . . . .
    // + + . . . .

    bili5[0][0] = Z[ipos_C(i,j)] ;
    bili5[0][1] = Z[ipos_C(i,j+1)] ;
    bili5[0][2] = DY[ipos_C(i,j)] ;
    bili5[0][3] = DY[ipos_C(i,j+1)] ;
    bili5[0][4] = DYY[ipos_C(i,j)] ;
    bili5[0][5] = DYY[ipos_C(i,j+1)] ;

    bili5[1][0] = Z[ipos_C(i+1,j)];
    bili5[1][1] = Z[ipos_C(i+1,j+1)];
    bili5[1][2] = DY[ipos_C(i+1,j)] ;
    bili5[1][3] = DY[ipos_C(i+1,j+1)];
    bili5[1][4] = DYY[ipos_C(i+1,j)] ;
    bili5[1][5] = DYY[ipos_C(i+1,j+1)];

    bili5[2][0] = DX[ipos_C(i,j)] ;
    bili5[2][1] = DX[ipos_C(i,j+1)] ;
    bili5[2][2] = DXY[ipos_C(i,j)] ;
    bili5[2][3] = DXY[ipos_C(i,j+1)] ;
    bili5[2][4] = DXXY[ipos_C(i,j)] ;
    bili5[2][5] = DXXY[ipos_C(i,j+1)] ;

    bili5[3][0] = DX[ipos_C(i+1,j)];
    bili5[3][1] = DX[ipos_C(i+1,j+1)];
    bili5[3][2] = DXY[ipos_C(i+1,j)] ;
    bili5[3][3] = DXY[ipos_C(i+1,j+1)] ;
    bili5[3][4] = DXXY[ipos_C(i+1,j)] ;
    bili5[3][5] = DXXY[ipos_C(i+1,j+1)] ;

    bili5[4][0] = DXX[ipos_C(i,j)] ;
    bili5[4][1] = DXX[ipos_C(i,j+1)] ;
    bili5[4][2] = DXXY[ipos_C(i,j)] ;
    bili5[4][3] = DXXY[ipos_C(i,j+1)] ;
    bili5[4][4] = 0 ; // DXXXY[ipos_C(i,j)] ;
    bili5[4][5] = 0 ; // DXXXY[ipos_C(i,j+1)] ;

    bili5[5][0] = DXX[ipos_C(i+1,j)];
    bili5[5][1] = DXX[ipos_C(i+1,j+1)];
    bili5[5][2] = DXXY[ipos_C(i+1,j)] ;
    bili5[5][3] = DXXY[ipos_C(i+1,j+1)] ;
    bili5[5][4] = 0 ; // DXXXY[ipos_C(i+1,j)] ;
    bili5[5][5] = 0 ; // DXXXY[ipos_C(i+1,j+1)] ;

  }

  valueType
  BiQuinticSplineBase::operator () ( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite5( x - X[i], X[i+1] - X[i], u ) ;
    Hermite5( y - Y[j], Y[j+1] - Y[j], v ) ;
    load(i,j) ;
    return bilinear5( u, bili5, v ) ;
  }

  valueType
  BiQuinticSplineBase::Dx( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite5_D( x - X[i], X[i+1] - X[i], u_D ) ;
    Hermite5  ( y - Y[j], Y[j+1] - Y[j], v   ) ;
    load(i,j) ;
    return bilinear5( u_D, bili5, v ) ;
  }

  valueType
  BiQuinticSplineBase::Dy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite5  ( x - X[i], X[i+1] - X[i], u   ) ;
    Hermite5_D( y - Y[j], Y[j+1] - Y[j], v_D ) ;
    load(i,j) ;
    return bilinear5( u, bili5, v_D ) ;
  }

  valueType
  BiQuinticSplineBase::Dxy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite5_D( x - X[i], X[i+1] - X[i], u_D ) ;
    Hermite5_D( y - Y[j], Y[j+1] - Y[j], v_D ) ;
    load(i,j) ;
    return bilinear5( u_D, bili5, v_D ) ;
  }

  valueType
  BiQuinticSplineBase::Dxx( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite5_DD( x - X[i], X[i+1] - X[i], u_DD ) ;
    Hermite5   ( y - Y[j], Y[j+1] - Y[j], v    ) ;
    load(i,j) ;
    return bilinear5( u_DD, bili5, v ) ;
  }

  valueType
  BiQuinticSplineBase::Dyy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite5   ( x - X[i], X[i+1] - X[i], u    ) ;
    Hermite5_DD( y - Y[j], Y[j+1] - Y[j], v_DD ) ;
    load(i,j) ;
    return bilinear5( u, bili5, v_DD ) ;
  }

  void
  BiQuinticSplineBase::D( valueType x, valueType y, valueType d[3] ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite5   ( x - X[i], X[i+1] - X[i], u    ) ;
    Hermite5_D ( x - X[i], X[i+1] - X[i], u_D  ) ;
    Hermite5   ( y - Y[j], Y[j+1] - Y[j], v    ) ;
    Hermite5_D ( y - Y[j], Y[j+1] - Y[j], v_D  ) ;
    load(i,j) ;
    d[0] = bilinear5( u, bili5, v ) ;
    d[1] = bilinear5( u_D, bili5, v ) ;
    d[2] = bilinear5( u, bili5, v_D ) ;
  }

  void
  BiQuinticSplineBase::DD( valueType x, valueType y, valueType d[6] ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite5   ( x - X[i], X[i+1] - X[i], u    ) ;
    Hermite5_D ( x - X[i], X[i+1] - X[i], u_D  ) ;
    Hermite5_DD( x - X[i], X[i+1] - X[i], u_DD ) ;
    Hermite5   ( y - Y[j], Y[j+1] - Y[j], v    ) ;
    Hermite5_D ( y - Y[j], Y[j+1] - Y[j], v_D  ) ;
    Hermite5_DD( y - Y[j], Y[j+1] - Y[j], v_DD ) ;
    load(i,j) ;
    d[0] = bilinear5( u, bili5, v ) ;
    d[1] = bilinear5( u_D, bili5, v ) ;
    d[2] = bilinear5( u, bili5, v_D ) ;
    d[3] = bilinear5( u_DD, bili5, v ) ;
    d[4] = bilinear5( u_D, bili5, v_D ) ;
    d[5] = bilinear5( u, bili5, v_DD ) ;
  }

}
