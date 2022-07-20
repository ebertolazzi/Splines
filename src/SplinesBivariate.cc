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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

namespace Splines {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SplineSurf::~SplineSurf()
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  SplineSurf::info() const {
    return fmt::format(
      "Bivariate spline [{}] of type = {}", name(), type_name()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineSurf::clear() {
    m_mem.free();
    m_nx = m_ny = 0;
    m_X = m_Y = m_Z = nullptr;
    m_Z_min = m_Z_max = 0;
    this->init_last_interval_x();
    this->init_last_interval_y();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void
  SplineSurf::load_Z(
    real_type const * z,
    integer           ldZ,
    bool              fortran_storage,
    bool              transposed
  ) {
    if ( transposed ) {
      if ( fortran_storage ) {
        UTILS_ASSERT(
          ldZ >= m_nx,
          "SplineSurf[{}]::load_Z[transposed+fortran_storage]\n"
          "ldZ = {} must be >= of nx = {}\n",
          m_name, ldZ, m_nx
        );
        for ( integer i = 0; i < m_nx; ++i )
          for ( integer j = 0; j < m_ny; ++j )
            m_Z[size_t(ipos_C(i,j))] = z[size_t(ipos_C(j,i,ldZ))];
      } else {
        UTILS_ASSERT(
          ldZ >= m_ny,
          "SplineSurf[{}]::load_Z[transposed]\n"
          "ldZ = {} must be >= of ny = {}\n",
          m_name, ldZ, m_ny
        );
        for ( integer i = 0; i < m_nx; ++i )
          for ( integer j = 0; j < m_ny; ++j )
            m_Z[size_t(ipos_C(i,j))] = z[size_t(ipos_F(j,i,ldZ))];
      }
    } else {
      if ( fortran_storage ) {
        UTILS_ASSERT(
          ldZ >= m_ny,
          "SplineSurf[{}]::load_Z[fortran_storage]\n"
          "ldZ = {} must be >= of ny = {}\n",
          m_name, ldZ, m_ny
        );
        for ( integer i = 0; i < m_nx; ++i )
          for ( integer j = 0; j < m_ny; ++j )
            m_Z[size_t(ipos_C(i,j))] = z[size_t(ipos_C(i,j,ldZ))];
      } else {
        UTILS_ASSERT(
          ldZ >= m_nx,
          "SplineSurf[{}]::load_Z\n"
          "ldZ = {} must be >= of nx = {}\n",
          m_name, ldZ, m_nx
        );
        for ( integer i = 0; i < m_nx; ++i )
          for ( integer j = 0; j < m_ny; ++j )
            m_Z[size_t(ipos_C(i,j))] = z[size_t(ipos_F(i,j,ldZ))];
      }
    }
    m_Z_max = *std::max_element(m_Z,m_Z+m_nx*m_ny);
    m_Z_min = *std::min_element(m_Z,m_Z+m_nx*m_ny);
    makeSpline();
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Build a spline surface with data
  //!
  //! - `x` the vector of x-nodes (size `nx`)
  //! - `y` the vector of y-nodes (size `ny`)
  //! - `z` the matrix of z-values, \f$ z_{ij} = f(x_i,y_j) \f$
  //!
  //! \code{.unparsed}
  //!            +----------+
  //!            |          |  internally data is stored by column
  //!   index j  ny  zij    |  zij = data[ i*ny + j ]
  //!            |          |
  //!            +--- nx ---+
  //!              index i
  //! \endcode
  //!
  //! - `fortran_storage` is `true` if matrix `z` is stored by column
  //! - `transposed` is true means that data are stored transposed
  //!
  void
  SplineSurf::build(
    real_type const * x, integer incx,
    real_type const * y, integer incy,
    real_type const * z, integer ldZ,
    integer           nx,
    integer           ny,
    bool              fortran_storage,
    bool              transposed
  ) {
    m_nx = nx;
    m_ny = ny;
    m_mem.reallocate( size_t((nx+1)*(ny+1)) );
    m_X = m_mem( size_t(nx) );
    m_Y = m_mem( size_t(ny) );
    m_Z = m_mem( size_t(nx*ny) );
    for ( size_t i = 0; i < size_t(nx); ++i ) m_X[i] = x[i*size_t(incx)];
    for ( size_t i = 0; i < size_t(ny); ++i ) m_Y[i] = y[i*size_t(incy)];
    load_Z( z, ldZ, fortran_storage, transposed );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Build a spline surface with data
  //!
  //! - `z`    the matrix of z-values, \f$ z_{ij} = f(x_i,y_j) \f$
  //! - `ldZ`  leading dimension of matriz `z`
  //! - `nx`   number of nodes in x-direction
  //! - `ny`   number of nodes in y-direction
  //!
  //! nodes are equispaced in x and y directions.
  //!
  //! \code{.unparsed}
  //!            +----------+
  //!            |          |  internally data is stored by column
  //!   index j  ny  zij    |  zij = data[ i*ny + j ]
  //!            |          |
  //!            +--- nx ---+
  //!              index i
  //! \endcode
  //!
  //! - `fortran_storage` is `true` if matrix `z` is stored by column
  //! - `transposed` is true means that data are stored transposed
  //!
  void
  SplineSurf::build(
    real_type const * z,
    integer           ldZ,
    integer           nx,
    integer           ny,
    bool              fortran_storage,
    bool              transposed
  ) {
    m_nx = nx;
    m_ny = ny;
    m_mem.reallocate( size_t((nx+1)*(ny+1)) );
    m_X = m_mem( size_t(nx) );
    m_Y = m_mem( size_t(ny) );
    m_Z = m_mem( size_t(nx*ny) );
    for ( size_t i = 0; i < size_t(nx); ++i ) m_X[i] = real_type(i);
    for ( size_t i = 0; i < size_t(ny); ++i ) m_Y[i] = real_type(i);
    load_Z( z, ldZ, fortran_storage, transposed );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSplineBase::load(
    integer i, integer j, real_type bili3[4][4]
  ) const {
    //
    //  1    3
    //
    //  0    2
    //
    size_t i0 = size_t(ipos_C(i,j));
    size_t i1 = size_t(ipos_C(i,j+1));
    size_t i2 = size_t(ipos_C(i+1,j));
    size_t i3 = size_t(ipos_C(i+1,j+1));

    bili3[0][0] = m_Z[i0];   bili3[0][1] = m_Z[i1];
    bili3[0][2] = m_DY[i0];  bili3[0][3] = m_DY[i1];

    bili3[1][0] = m_Z[i2];   bili3[1][1] = m_Z[i3];
    bili3[1][2] = m_DY[i2];  bili3[1][3] = m_DY[i3];

    bili3[2][0] = m_DX[i0];  bili3[2][1] = m_DX[i1];
    bili3[2][2] = m_DXY[i0]; bili3[2][3] = m_DXY[i1];

    bili3[3][0] = m_DX[i2];  bili3[3][1] = m_DX[i3];
    bili3[3][2] = m_DXY[i2]; bili3[3][3] = m_DXY[i3];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::operator () ( real_type x, real_type y ) const {
    real_type bili3[4][4], u[4], v[4];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u );
    Hermite3( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v );
    load( i, j, bili3 );
    return bilinear3( u, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dx( real_type x, real_type y ) const {
    real_type bili3[4][4], u_D[4], v[4];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3_D( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_D );
    Hermite3  ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v   );
    load( i, j, bili3 );
    return bilinear3( u_D, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dy( real_type x, real_type y ) const {
    real_type bili3[4][4], u[4], v_D[4];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3  ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u   );
    Hermite3_D( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_D );
    load( i, j, bili3 );
    return bilinear3( u, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dxy( real_type x, real_type y ) const {
    real_type bili3[4][4], u_D[4], v_D[4];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3_D( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_D );
    Hermite3_D( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_D );
    load( i, j, bili3 );
    return bilinear3( u_D, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dxx( real_type x, real_type y ) const {
    real_type bili3[4][4], u_DD[4], v[4];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3_DD( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_DD );
    Hermite3   ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v    );
    load( i, j, bili3 );
    return bilinear3( u_DD, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dyy( real_type x, real_type y ) const {
    real_type bili3[4][4], u[4], v_DD[4];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3   ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u    );
    Hermite3_DD( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_DD );
    load( i, j, bili3 );
    return bilinear3( u, bili3, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSplineBase::D( real_type x, real_type y, real_type d[3] ) const {
    real_type bili3[4][4], u[4], u_D[4], v[3], v_D[4];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3   ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u    );
    Hermite3_D ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_D  );
    Hermite3   ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v    );
    Hermite3_D ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_D  );
    load( i, j, bili3 );
    d[0] = bilinear3( u, bili3, v );
    d[1] = bilinear3( u_D, bili3, v );
    d[2] = bilinear3( u, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSplineBase::DD( real_type x, real_type y, real_type d[6] ) const {
    real_type bili3[4][4], u[4], u_D[4], u_DD[4], v[3], v_D[4], v_DD[4];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite3   ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u    );
    Hermite3_D ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_D  );
    Hermite3_DD( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_DD );
    Hermite3   ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v    );
    Hermite3_D ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_D  );
    Hermite3_DD( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_DD );
    load( i, j, bili3 );
    d[0] = bilinear3( u, bili3, v );
    d[1] = bilinear3( u_D, bili3, v );
    d[2] = bilinear3( u, bili3, v_D );
    d[3] = bilinear3( u_DD, bili3, v );
    d[4] = bilinear3( u_D, bili3, v_D );
    d[5] = bilinear3( u, bili3, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSplineBase::load(
    integer i, integer j, real_type bili5[6][6]
  ) const {

    size_t i00 = size_t(ipos_C(i,j));
    size_t i01 = size_t(ipos_C(i,j+1));
    size_t i10 = size_t(ipos_C(i+1,j));
    size_t i11 = size_t(ipos_C(i+1,j+1));

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
    bili5[0][0] = m_Z[i00]; bili5[0][1] = m_Z[i01];
    bili5[1][0] = m_Z[i10]; bili5[1][1] = m_Z[i11];

    // DX
    bili5[2][0] = m_DX[i00]; bili5[2][1] = m_DX[i01];
    bili5[3][0] = m_DX[i10]; bili5[3][1] = m_DX[i11];

    // DXX
    bili5[4][0] = m_DXX[i00]; bili5[4][1] = m_DXX[i01];
    bili5[5][0] = m_DXX[i10]; bili5[5][1] = m_DXX[i11];

    // DY
    bili5[0][2] = m_DY[i00]; bili5[0][3] = m_DY[i01];
    bili5[1][2] = m_DY[i10]; bili5[1][3] = m_DY[i11];

    // DYY
    bili5[0][4] = m_DYY[i00]; bili5[0][5] = m_DYY[i01];
    bili5[1][4] = m_DYY[i10]; bili5[1][5] = m_DYY[i11];

    // DXY
    bili5[2][2] = m_DXY[i00]; bili5[2][3] = m_DXY[i01];
    bili5[3][2] = m_DXY[i10]; bili5[3][3] = m_DXY[i11];

    // DXXY
    bili5[4][2] = m_DXXY[i00]; bili5[4][3] = m_DXXY[i01];
    bili5[5][2] = m_DXXY[i10]; bili5[5][3] = m_DXXY[i11];

    // DXYY
    bili5[2][4] = m_DXYY[i00]; bili5[2][5] = m_DXYY[i01];
    bili5[3][4] = m_DXYY[i10]; bili5[3][5] = m_DXYY[i11];

    // DXXYY
    bili5[4][4] = m_DXXYY[i00]; bili5[4][5] = m_DXXYY[i01];
    bili5[5][4] = m_DXXYY[i10]; bili5[5][5] = m_DXXYY[i11];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::operator () ( real_type x, real_type y ) const {
    real_type bili5[6][6], u[6], v[6];
    size_t i = size_t(search_x( x ));
    size_t j = size_t(search_y( y ));
    Hermite5( x - m_X[i], m_X[i+1] - m_X[i], u );
    Hermite5( y - m_Y[j], m_Y[j+1] - m_Y[j], v );
    load( integer(i), integer(j), bili5 );
    return bilinear5( u, bili5, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dx( real_type x, real_type y ) const {
    real_type bili5[6][6], u_D[6], v[6];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5_D( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_D );
    Hermite5  ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v   );
    load( integer(i), integer(j), bili5 );
    return bilinear5( u_D, bili5, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dy( real_type x, real_type y ) const {
    real_type bili5[6][6], u[6], v_D[6];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5  ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u   );
    Hermite5_D( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_D );
    load( integer(i), integer(j), bili5 );
    return bilinear5( u, bili5, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dxy( real_type x, real_type y ) const {
    real_type bili5[6][6], u_D[6], v_D[6];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5_D( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_D );
    Hermite5_D( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_D );
    load( integer(i), integer(j), bili5 );
    return bilinear5( u_D, bili5, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dxx( real_type x, real_type y ) const {
    real_type bili5[6][6], u_DD[6], v[6];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5_DD( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_DD );
    Hermite5   ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v    );
    load( integer(i), integer(j), bili5 );
    return bilinear5( u_DD, bili5, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dyy( real_type x, real_type y ) const {
    real_type bili5[6][6], u[6], v_DD[6];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5   ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u    );
    Hermite5_DD( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_DD );
    load( integer(i), integer(j), bili5 );
    return bilinear5( u, bili5, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSplineBase::D( real_type x, real_type y, real_type d[3] ) const {
    real_type bili5[6][6], u[6], u_D[6], v[6], v_D[6];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5   ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u    );
    Hermite5_D ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_D  );
    Hermite5   ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v    );
    Hermite5_D ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_D  );
    load( integer(i), integer(j), bili5 );
    d[0] = bilinear5( u, bili5, v );
    d[1] = bilinear5( u_D, bili5, v );
    d[2] = bilinear5( u, bili5, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSplineBase::DD( real_type x, real_type y, real_type d[6] ) const {
    real_type bili5[6][6], u[6], u_D[6], u_DD[6], v[6], v_D[6], v_DD[6];
    integer i = search_x( x );
    integer j = search_y( y );
    Hermite5   ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u    );
    Hermite5_D ( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_D  );
    Hermite5_DD( x - m_X[size_t(i)], m_X[size_t(i+1)] - m_X[size_t(i)], u_DD );
    Hermite5   ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v    );
    Hermite5_D ( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_D  );
    Hermite5_DD( y - m_Y[size_t(j)], m_Y[size_t(j+1)] - m_Y[size_t(j)], v_DD );
    load( integer(i), integer(j), bili5 );
    d[0] = bilinear5( u, bili5, v );
    d[1] = bilinear5( u_D, bili5, v );
    d[2] = bilinear5( u, bili5, v_D );
    d[3] = bilinear5( u_DD, bili5, v );
    d[4] = bilinear5( u_D, bili5, v_D );
    d[5] = bilinear5( u, bili5, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  SplineSurf::dump_data( ostream_type & s ) const {
    s << "X = [ " << m_X[0];
    for ( integer i = 1; i < m_nx; ++i ) s << ", " << m_X[size_t(i)];
    s << " ]\nY = [ " << m_Y[0];
    for ( integer i = 1; i < m_ny; ++i ) s << ", " << m_Y[size_t(i)];
    s << " ]\nZ = [\n";
    for ( integer j = 0; j < m_ny; ++j ) {
      s << "  [ " << m_Z[size_t(ipos_C(0,j))];
      for ( integer i = 1; i < m_nx; ++i )
        s << ", " << m_Z[size_t(ipos_C(i,j))];
      s << " ]\n";
    }
    s << "\n];\n";
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  using GC_namespace::GC_VECTOR;
  using GC_namespace::GC_VEC_INTEGER;
  using GC_namespace::GC_VEC_LONG;
  using GC_namespace::GC_VEC_REAL;
  using GC_namespace::GC_MAT_REAL;
  using GC_namespace::mat_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Setup a spline surface using a `GenericContainer`
  //!
  //! - `gc("fortran_storage")` if true `zdata` is stored by column, otherwise by rows
  //! - `gc("transposed")`      if true `zdata` is stored transposed
  //!
  //! - `gc("xdata")` vector of mesh point in `x` direction
  //! - `gc("ydata")` vector of mesh point in `y` direction
  //! - `gc("zdata")` may be
  //!    - matrix of size `nx` x `ny` (`ny` x `nx` if data is transposed)
  //!    - vector of size `nx` x `ny` stiring the data.
  //!
  //! - `nx` number of nodes in x-direction the size of vector in `x` direction
  //! - `ny` number of nodes in y-direction the size of vector in `y` direction
  //!
  //! nodes are equispaced in x and y directions.
  //!
  //! \code{.unparsed}
  //!            +----------+
  //!            |          |  internally data is stored by column
  //!   index j  ny  zij    |  zij = data[ i*ny + j ]
  //!            |          |
  //!            +--- nx ---+
  //!              index i
  //! \endcode
  //!
  //! - `fortran_storage` is `true` if matrix `z` is stored by column
  //! - `transposed` is true means that data are stored transposed
  //!
  void
  SplineSurf::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    // gc["zdata"]
    //
    */
    string msg = fmt::format("SplineSurf[{}]::setup( gc ):", m_name );
    UTILS_ASSERT( gc.exists("xdata"), "{}, missing `xdata` field!\n", msg );
    UTILS_ASSERT( gc.exists("ydata"), "{}, missing `ydata` field!\n", msg );
    UTILS_ASSERT( gc.exists("zdata"), "{}, missing `zdata` field!\n", msg );

    GenericContainer const & gc_x = gc("xdata");
    GenericContainer const & gc_y = gc("ydata");
    GenericContainer const & gc_z = gc("zdata");

    m_nx = gc_x.get_num_elements();
    m_ny = gc_y.get_num_elements();
    m_mem.reallocate( size_t((m_nx+1)*(m_ny+1)) );
    m_X = m_mem( size_t(m_nx) );
    m_Y = m_mem( size_t(m_ny) );
    m_Z = m_mem( size_t(m_nx*m_ny) );

    for ( integer i = 0; i < m_nx; ++i ) m_X[size_t(i)] = gc_x.get_number_at(i);
    for ( integer i = 0; i < m_ny; ++i ) m_Y[size_t(i)] = gc_y.get_number_at(i);

    bool fortran_storage = false;
    bool transposed      = false;
    gc.get_if_exists("fortran_storage",fortran_storage);
    gc.get_if_exists("transposed",transposed);


    /*
    //     +------+
    //  j ny      | (xi,yj)
    //     +  nx  +
    //         i
    */

    // cosa mi aspetto in lettura
    integer N, M;
    if ( transposed ) { N = m_ny; M = m_nx; }
    else              { N = m_nx; M = m_ny; }
    integer LD = fortran_storage ? N : M;

    if ( GC_MAT_REAL == gc_z.get_type() ) {
      mat_real_type const & z = gc_z.get_mat_real();
      UTILS_ASSERT(
        unsigned(N) == z.numRows() && unsigned(M) == z.numCols(),
        "{}, field `z` expected to be of size {} x {}, found: {} x {}\n",
        msg, N, M, z.numRows(), z.numCols()
      );
      load_Z( m_Z, LD, fortran_storage, transposed );
    } else if ( GC_VEC_INTEGER == gc_z.get_type() ||
                GC_VEC_LONG    == gc_z.get_type() ||
                GC_VEC_REAL    == gc_z.get_type() ) {

      integer nz  = gc_z.get_num_elements();
      integer nxy = m_nx * m_ny;
      UTILS_ASSERT(
        nz == nxy,
        "{}, field `z` expected to be of size {} = {}x{}, found: `{}`\n",
        msg, nxy, m_nx, m_ny, nz
      );
      for ( integer i = 0; i < nz ; ++i ) m_Z[size_t(i)] = gc_z.get_number_at(i);
      load_Z( m_Z, LD, fortran_storage, transposed );
    } else if ( GC_VECTOR == gc_z.get_type() ) {
      vector_type const & data = gc_z.get_vector();
      vec_real_type tmp;
      UTILS_ASSERT(
        size_t(M) == data.size(),
        "{}, field `zdata` (vector of vector) expected of size {} found of size {}\n",
        msg, M, data.size()
      );
      for ( integer j = 0; j < M; ++j ) {
        GenericContainer const & row = data[size_t(j)];
        string msg1 = fmt::format( "{}, reading row {}\n", msg, j );
        row.copyto_vec_real( tmp, msg1.c_str() );
        UTILS_ASSERT(
          size_t(N) == tmp.size(),
          "{}, row {}-th of size {}, expected {}\n",
          msg, j, tmp.size(), N
        );
        if ( transposed ) {
          for ( integer i = 0; i < N; ++i )
            m_Z[size_t(ipos_C(j,i))] = tmp[size_t(i)];
        } else {
          for ( integer i = 0; i < N; ++i )
            m_Z[size_t(ipos_C(i,j))] = tmp[size_t(i)];
        }
      }
      load_Z( m_Z, m_ny, true, false );
    } else {
      UTILS_ERROR(
        "{}, field `z` expected to be of type"
        " `mat_real_type` or  `vec_real_type` or `vector_type` found: `{}`\n",
        msg, gc_z.get_type_name()
      );
    }

  }

}
