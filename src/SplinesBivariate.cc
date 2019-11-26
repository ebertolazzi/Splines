/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2016                                                      |
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

namespace Splines {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SplineSurf::~SplineSurf()
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSurf::info( ostream & s ) const {
    s << "Bivariate spline [" << name() << "] of type = "
      << type_name()
      << '\n';
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSurf::clear(void) {
    X.clear();
    Y.clear();
    Z.clear();
    Z_min = Z_max = 0;
    {
      std::unique_lock<std::mutex> lck(lastInterval_x_mutex);
      lastInterval_x_by_thread[std::this_thread::get_id()] = 0;
    }
    {
      std::unique_lock<std::mutex> lck(lastInterval_y_mutex);
      lastInterval_y_by_thread[std::this_thread::get_id()] = 0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSurf::build(
    real_type const x[], integer incx,
    real_type const y[], integer incy,
    real_type const z[], integer ldZ,
    integer nx, integer ny,
    bool fortran_storage,
    bool transposed
  ) {
    X.resize( size_t(nx) );
    Y.resize( size_t(ny) );
    Z.resize( size_t(nx*ny) );
    for ( size_t i = 0; i < size_t(nx); ++i ) X[i] = x[i*size_t(incx)];
    for ( size_t i = 0; i < size_t(ny); ++i ) Y[i] = y[i*size_t(incy)];
    if ( (fortran_storage && transposed) || (!fortran_storage && !transposed) ) {
      SPLINE_ASSERT(
        ldZ >= ny,
        "SplineSurf::build, ldZ = " << ldZ << " must be >= of nx = " << ny
      )
      for ( integer i = 0; i < nx; ++i )
        for ( integer j = 0; j < ny; ++j )
          Z[size_t(ipos_C(i,j,ny))] = z[size_t(ipos_C(i,j,ldZ))];
    } else {
      SPLINE_ASSERT(
        ldZ >= nx,
        "SplineSurf::build, ldZ = " << ldZ << " must be >= of ny = " << nx
      )
      for ( integer i = 0; i < nx; ++i )
        for ( integer j = 0; j < ny; ++j )
          Z[size_t(ipos_C(i,j,ny))] = z[size_t(ipos_F(i,j,ldZ))];
    }
    Z_max = *std::max_element(Z.begin(),Z.end());
    Z_min = *std::min_element(Z.begin(),Z.end());
    makeSpline();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSurf::build(
    real_type const z[], integer ldZ,
    integer nx, integer ny,
    bool fortran_storage,
    bool transposed
  ) {
    vector<real_type> XX, YY; // temporary vector
    XX.resize( size_t(nx) );
    YY.resize( size_t(ny) ); // temporary vector
    for ( size_t i = 0; i < size_t(nx); ++i ) XX[i] = real_type(i);
    for ( size_t i = 0; i < size_t(ny); ++i ) YY[i] = real_type(i);
    build(
      &XX.front(), 1, &YY.front(), 1, z, ldZ, nx, ny,
      fortran_storage, transposed
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSplineBase::load( integer i, integer j ) const {

    //
    //  1    3
    //
    //  0    2
    //
    integer ny = integer(Y.size());
    size_t i0 = size_t(ipos_C(i,j,ny));
    size_t i1 = size_t(ipos_C(i,j+1,ny));
    size_t i2 = size_t(ipos_C(i+1,j,ny));
    size_t i3 = size_t(ipos_C(i+1,j+1,ny));

    bili3[0][0] = Z[i0];   bili3[0][1] = Z[i1];
    bili3[0][2] = DY[i0];  bili3[0][3] = DY[i1];

    bili3[1][0] = Z[i2];   bili3[1][1] = Z[i3];
    bili3[1][2] = DY[i2];  bili3[1][3] = DY[i3];

    bili3[2][0] = DX[i0];  bili3[2][1] = DX[i1];
    bili3[2][2] = DXY[i0]; bili3[2][3] = DXY[i1];

    bili3[3][0] = DX[i2];  bili3[3][1] = DX[i3];
    bili3[3][2] = DXY[i2]; bili3[3][3] = DXY[i3];

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::operator () ( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u );
    Hermite3( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v );
    load(i,j);
    return bilinear3( u, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dx( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3_D( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_D );
    Hermite3  ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v   );
    load(i,j);
    return bilinear3( u_D, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dy( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3  ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u   );
    Hermite3_D( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_D );
    load(i,j);
    return bilinear3( u, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dxy( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3_D( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_D );
    Hermite3_D( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_D );
    load(i,j);
    return bilinear3( u_D, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dxx( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3_DD( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_DD );
    Hermite3   ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v    );
    load(i,j);
    return bilinear3( u_DD, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dyy( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3   ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u    );
    Hermite3_DD( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_DD );
    load(i,j);
    return bilinear3( u, bili3, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSplineBase::D( real_type x, real_type y, real_type d[3] ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3   ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u    );
    Hermite3_D ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_D  );
    Hermite3   ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v    );
    Hermite3_D ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_D  );
    load(i,j);
    d[0] = bilinear3( u, bili3, v );
    d[1] = bilinear3( u_D, bili3, v );
    d[2] = bilinear3( u, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSplineBase::DD( real_type x, real_type y, real_type d[6] ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3   ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u    );
    Hermite3_D ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_D  );
    Hermite3_DD( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_DD );
    Hermite3   ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v    );
    Hermite3_D ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_D  );
    Hermite3_DD( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_DD );
    load(i,j);
    d[0] = bilinear3( u, bili3, v );
    d[1] = bilinear3( u_D, bili3, v );
    d[2] = bilinear3( u, bili3, v_D );
    d[3] = bilinear3( u_DD, bili3, v );
    d[4] = bilinear3( u_D, bili3, v_D );
    d[5] = bilinear3( u, bili3, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSplineBase::load( integer i, integer j ) const {

    integer ny = integer(Y.size());
    size_t i00 = size_t(ipos_C(i,j,ny));
    size_t i01 = size_t(ipos_C(i,j+1,ny));
    size_t i10 = size_t(ipos_C(i+1,j,ny));
    size_t i11 = size_t(ipos_C(i+1,j+1,ny));

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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::operator () ( real_type x, real_type y ) const {
    size_t i = size_t(search_x( x ));
    size_t j = size_t(search_y( y ));
    Hermite5( x - X[i], X[i+1] - X[i], u );
    Hermite5( y - Y[j], Y[j+1] - Y[j], v );
    load(integer(i),integer(j));
    return bilinear5( u, bili5, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dx( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5_D( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_D );
    Hermite5  ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v   );
    load(i,j);
    return bilinear5( u_D, bili5, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dy( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5  ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u   );
    Hermite5_D( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_D );
    load(i,j);
    return bilinear5( u, bili5, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dxy( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5_D( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_D );
    Hermite5_D( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_D );
    load(i,j);
    return bilinear5( u_D, bili5, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dxx( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5_DD( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_DD );
    Hermite5   ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v    );
    load(i,j);
    return bilinear5( u_DD, bili5, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dyy( real_type x, real_type y ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5   ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u    );
    Hermite5_DD( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_DD );
    load(i,j);
    return bilinear5( u, bili5, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSplineBase::D( real_type x, real_type y, real_type d[3] ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5   ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u    );
    Hermite5_D ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_D  );
    Hermite5   ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v    );
    Hermite5_D ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_D  );
    load(i,j);
    d[0] = bilinear5( u, bili5, v );
    d[1] = bilinear5( u_D, bili5, v );
    d[2] = bilinear5( u, bili5, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSplineBase::DD( real_type x, real_type y, real_type d[6] ) const {
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5   ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u    );
    Hermite5_D ( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_D  );
    Hermite5_DD( x - X[size_t(i)], X[size_t(i+1)] - X[size_t(i)], u_DD );
    Hermite5   ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v    );
    Hermite5_D ( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_D  );
    Hermite5_DD( y - Y[size_t(j)], Y[size_t(j+1)] - Y[size_t(j)], v_DD );
    load(i,j);
    d[0] = bilinear5( u, bili5, v );
    d[1] = bilinear5( u_D, bili5, v );
    d[2] = bilinear5( u, bili5, v_D );
    d[3] = bilinear5( u_DD, bili5, v );
    d[4] = bilinear5( u_D, bili5, v_D );
    d[5] = bilinear5( u, bili5, v_DD );
  }

  #ifndef SPLINES_DO_NOT_USE_GENERIC_CONTAINER
  using GenericContainerNamespace::GC_VEC_REAL;
  using GenericContainerNamespace::GC_VEC_INTEGER;
  using GenericContainerNamespace::GC_MAT_REAL;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSurf::setup( GenericContainer const & gc ) {
    /*
    // gc["x"]
    // gc["y"]
    // gc["z"]
    //
    */
    SPLINE_ASSERT(
      gc.exists("x"),
      "[SplineSurf[" << _name << "]::setup] missing `x` field!"
    )
    SPLINE_ASSERT(
      gc.exists("y"),
      "[SplineSurf[" << _name << "]::setup] missing `y` field!"
    )
    SPLINE_ASSERT(
      gc.exists("z"),
      "[SplineSurf[" << _name << "]::setup] missing `z` field!"
    )

    GenericContainer const & gc_x = gc("x");
    GenericContainer const & gc_y = gc("y");
    GenericContainer const & gc_z = gc("z");

    vec_real_type x, y;
    gc_x.copyto_vec_real( x, "SplineSurf::setup, field `x'" );
    gc_y.copyto_vec_real( y, "SplineSurf::setup, field `y'" );

    bool fortran_storage = false;
    if ( gc.exists("fortran_storage") )
      fortran_storage = gc("fortran_storage").get_bool();

    bool transposed = false;
    if ( gc.exists("transposed") )
      transposed = gc("transposed").get_bool();

    integer nx = integer(x.size());
    integer ny = integer(y.size());

    if ( GC_MAT_REAL == gc_z.get_type() ) {
      build(
        &x.front(), 1,
        &y.front(), 1,
        &gc_z.get_mat_real().front(),
        integer(gc_z.get_mat_real().numRows()),
        nx, ny, fortran_storage, transposed
      );
    } else if ( GC_VEC_REAL == gc_z.get_type() ) {
      GenericContainer const & gc_ldz = gc("ldz");
      integer ldz = integer(gc_ldz.get_as_uint("SplineSurf::setup, field `ldz` expected to be and integer"));
      vec_real_type z;
      gc_z.copyto_vec_real( z, "SplineSurf::setup, field `z'" );
      build(
        &x.front(), 1,
        &y.front(), 1,
        &z.front(), ldz,
        nx, ny, fortran_storage, transposed
      );
    } else {
      SPLINE_ASSERT(
        false,
        "[SplineSurf[" << _name <<
        "]::setup] field `z` expected to be of type `mat_real_type` or  `vec_real_type` found: `" <<
        gc_z.get_type_name() << "`"
      )
    }

  }
  #endif

}
