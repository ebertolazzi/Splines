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

  void
  SplineSurf::info( ostream & s ) const {
    s << "Bivariate spline [" << name() << "] of type = "
      << type_name()
      << '\n' ;
  }

  #ifdef SPLINES_USE_GENERIC_CONTAINER
  
  using GenericContainerNamepace::GC_VEC_REAL ;
  using GenericContainerNamepace::GC_MAT_REAL ;

  void
  SplineSurf::setup( GenericContainer const & gc ) {
    /*
    // gc["x"]
    // gc["y"]
    // gc["z"]
    //
    */
    SPLINE_ASSERT( gc.exists("x"), "[" << _name << "] SplineSurf::build, missing `x` field!") ;
    SPLINE_ASSERT( gc.exists("y"), "[" << _name << "] SplineSurf::build, missing `y` field!") ;
    SPLINE_ASSERT( gc.exists("z"), "[" << _name << "] SplineSurf::build, missing `z` field!") ;
  
    GenericContainer const & gc_x = gc("x") ;
    GenericContainer const & gc_y = gc("y") ;
    GenericContainer const & gc_z = gc("z") ;

    SPLINE_ASSERT( GC_VEC_REAL == gc_x.get_type(),
                   "Field `x` expected to be of type `vec_real_type` found: `" <<
                   gc_x.get_type_name() << "`" ) ;

    SPLINE_ASSERT( GC_VEC_REAL == gc_y.get_type(),
                   "Field `y` expected to be of type `vec_real_type` found: `" <<
                   gc_y.get_type_name() << "`" ) ;

    SPLINE_ASSERT( GC_MAT_REAL == gc_z.get_type(),
                   "Field `z` expected to be of type `mat_real_type` found: `" <<
                   gc_z.get_type_name() << "`" ) ;

    bool fortran_storage = false ;
    if ( gc.exists("fortran_storage") )
      fortran_storage = gc("fortran_storage").get_bool() ;

    bool transposed = false ;
    if ( gc.exists("transposed") )
      transposed = gc("transposed").get_bool() ;

    sizeType nx = sizeType(gc_x.get_vec_real().size()) ;
    sizeType ny = sizeType(gc_y.get_vec_real().size()) ;

    build ( &gc_x.get_vec_real().front(), 1,
            &gc_y.get_vec_real().front(), 1,
            &gc_z.get_mat_real().front(), gc_z.get_mat_real().numRows(),
            nx, ny, fortran_storage, transposed ) ;

  }
  #endif

  sizeType
  SplineSurf::search_x( valueType x ) const {
    if ( lastInterval_x+1 >= X.size() || X[lastInterval_x] < x || X[lastInterval_x+1] > x ) {
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
    if ( lastInterval_y+1 >= Y.size() || Y[lastInterval_y] < y || Y[lastInterval_y+1] > y ) {
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

    //
    //  1    3
    //
    //  0    2
    //

    sizeType i0 = ipos_C(i,j) ;
    sizeType i1 = ipos_C(i,j+1) ;
    sizeType i2 = ipos_C(i+1,j) ;
    sizeType i3 = ipos_C(i+1,j+1) ;

    bili3[0][0] = Z[i0];   bili3[0][1] = Z[i1];
    bili3[0][2] = DY[i0];  bili3[0][3] = DY[i1];

    bili3[1][0] = Z[i2];   bili3[1][1] = Z[i3];
    bili3[1][2] = DY[i2];  bili3[1][3] = DY[i3];

    bili3[2][0] = DX[i0];  bili3[2][1] = DX[i1];
    bili3[2][2] = DXY[i0]; bili3[2][3] = DXY[i1];

    bili3[3][0] = DX[i2];  bili3[3][1] = DX[i3];
    bili3[3][2] = DXY[i2]; bili3[3][3] = DXY[i3];

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

    sizeType i00 = ipos_C(i,j) ;
    sizeType i01 = ipos_C(i,j+1) ;
    sizeType i10 = ipos_C(i+1,j) ;
    sizeType i11 = ipos_C(i+1,j+1) ;

    //
    //  1    3
    //
    //  0    2
    //
    // H0, H1, dH0, dH1, ddH0, ddH1

    // + + + + + +
    // + + + + + +
    // + + + + . .
    // + + + + . .
    // + + . . . .
    // + + . . . .

    // P
    bili5[0][0] = Z[i00]; bili5[0][1] = Z[i01];
    bili5[1][0] = Z[i10]; bili5[1][1] = Z[i11];

    // DX
    bili5[2][0] = DX[i00]; bili5[2][1] = DX[i01];
    bili5[3][0] = DX[i10]; bili5[3][1] = DX[i11];

    // DXX
    bili5[4][0] = DXX[i00]; bili5[4][1] = DXX[i01];
    bili5[5][0] = DXX[i10]; bili5[5][1] = DXX[i11];

    // DY
    bili5[0][2] = DY[i00]; bili5[0][3] = DY[i01];
    bili5[1][2] = DY[i10]; bili5[1][3] = DY[i11];

    // DYY
    bili5[0][4] = DYY[i00]; bili5[0][5] = DYY[i01];
    bili5[1][4] = DYY[i10]; bili5[1][5] = DYY[i11];

    // DXY
    bili5[2][2] = DXY[i00]; bili5[2][3] = DXY[i01];
    bili5[3][2] = DXY[i10]; bili5[3][3] = DXY[i11];

    // DXXY
    bili5[4][2] = DXXY[i00]; bili5[4][3] = DXXY[i01];
    bili5[5][2] = DXXY[i10]; bili5[5][3] = DXXY[i11];

    // DXYY
    bili5[2][4] = DXYY[i00]; bili5[2][5] = DXYY[i01];
    bili5[3][4] = DXYY[i10]; bili5[3][5] = DXYY[i11];

    // DXXYY
    bili5[4][4] = DXXYY[i00]; bili5[4][5] = DXXYY[i01];
    bili5[5][4] = DXXYY[i10]; bili5[5][5] = DXXYY[i11];
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
