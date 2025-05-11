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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"

#include <cmath>
#include <iomanip>

namespace Splines {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SplineSurf::~SplineSurf()
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  SplineSurf::info() const {
    return fmt::format( "Bivariate spline [{}] of type = {}", name(), type_name() );
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
    real_type const z[],
    integer   const ldZ,
    bool      const fortran_storage,
    bool      const transposed
  ) {
    // 
    //  +--------------+
    //  |              | ny = nr
    //  |              |
    //  +--------------+
    //      nx = nc
    //
    integer const nr{ transposed ? m_nx : m_ny };
    integer const nc{ transposed ? m_ny : m_nx };
    if ( fortran_storage ) {
      UTILS_ASSERT(
        ldZ >= nr,
        "SplineSurf[{}]::load_Z [fortran storage]\n"
        "ldZ = {} must be >= of nr, nr x nc = {} x {}\n",
        m_name, ldZ, nr, nc
      );
    } else {
      UTILS_ASSERT(
        ldZ >= nc,
        "SplineSurf[{}]::load_Z [C storage]\n"
        "ldZ = {} must be >= of nc, nr x nc= {} x {}\n",
        m_name, ldZ, nr, nc
      );
    }
    integer const tf{ (transposed ? 1 : 0) + (fortran_storage ? 2 : 0) };
    switch ( tf ) {
      case 0: // NO transpose NO fortran
        for ( integer ix{0}; ix < m_nx; ++ix )
          for ( integer iy{0}; iy < m_ny; ++iy )
            z_node_ref(ix,iy) = z[ ipos_C(iy,ix,ldZ) ];
      break;
      case 1: // YES transpose NO fortran
        for ( integer ix{0}; ix < m_nx; ++ix )
          for ( integer iy{0}; iy < m_ny; ++iy )
            z_node_ref(ix,iy) = z[ ipos_C(ix,iy,ldZ) ];
      break;
      case 2: // NO transpose YES fortran
        for ( integer ix{0}; ix < m_nx; ++ix )
          for ( integer iy{0}; iy < m_ny; ++iy )
            z_node_ref(ix,iy) = z[ ipos_F(iy,ix,ldZ) ];
      break;
      case 3: // YES transpose YES fortran
        for ( integer ix{0}; ix < m_nx; ++ix )
          for ( integer iy{0}; iy < m_ny; ++iy )
            z_node_ref(ix,iy) = z[ ipos_F(ix,iy,ldZ) ];
      break;
    }
    m_Z_max = *std::max_element(m_Z,m_Z+m_nx*m_ny);
    m_Z_min = *std::min_element(m_Z,m_Z+m_nx*m_ny);
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
    real_type const x[], integer const incx,
    real_type const y[], integer const incy,
    real_type const z[], integer const ldZ,
    integer   const nx,
    integer   const ny,
    bool      const fortran_storage,
    bool      const transposed
  ) {
    m_nx = nx;
    m_ny = ny;
    m_mem.reallocate( (nx+1)*(ny+1) );
    m_X = m_mem( nx );
    m_Y = m_mem( ny );
    m_Z = m_mem( nx*ny );
    for ( integer i{0}; i < nx; ++i ) m_X[i] = x[i*incx];
    for ( integer i{0}; i < ny; ++i ) m_Y[i] = y[i*incy];
    load_Z( z, ldZ, fortran_storage, transposed );
    make_spline();
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
    real_type const z[],
    integer   const ldZ,
    integer   const nx,
    integer   const ny,
    bool      const fortran_storage,
    bool      const transposed
  ) {
    m_nx = nx;
    m_ny = ny;
    m_mem.reallocate( (nx+1)*(ny+1) );
    m_X = m_mem( nx );
    m_Y = m_mem( ny );
    m_Z = m_mem( nx*ny );
    for ( integer i{0}; i < nx; ++i ) m_X[i] = static_cast<real_type>(i);
    for ( integer i{0}; i < ny; ++i ) m_Y[i] = static_cast<real_type>(i);
    load_Z( z, ldZ, fortran_storage, transposed );
    make_spline();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  BiCubicSplineBase::load(
    integer const i,
    integer const j,
    real_type     bili3[4][4]
  ) const {
    //
    //  1    3
    //
    //  0    2
    //
    integer const i0 { ipos_C(i,j) };
    integer const i1 { ipos_C(i,j+1) };
    integer const i2 { ipos_C(i+1,j) };
    integer const i3 { ipos_C(i+1,j+1) };

    bili3[0][0] = m_Z[i0];   bili3[0][1] = m_Z[i1];
    bili3[0][2] = m_DY[i0];  bili3[0][3] = m_DY[i1];

    bili3[1][0] = m_Z[i2];   bili3[1][1] = m_Z[i3];
    bili3[1][2] = m_DY[i2];  bili3[1][3] = m_DY[i3];

    bili3[2][0] = m_DX[i0];  bili3[2][1] = m_DX[i1];
    bili3[2][2] = m_DXY[i0]; bili3[2][3] = m_DXY[i1];

    bili3[3][0] = m_DX[i2];  bili3[3][1] = m_DX[i3];
    bili3[3][2] = m_DXY[i2]; bili3[3][3] = m_DXY[i3];
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::eval( real_type x, real_type y ) const {
    real_type bili3[4][4], u[4], v[4];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite3( x - m_X[i], m_X[i+1] - m_X[i], u );
    Hermite3( y - m_Y[j], m_Y[j+1] - m_Y[j], v );
    load( i, j, bili3 );
    return bilinear3( u, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dx( real_type x, real_type y ) const {
    real_type bili3[4][4], u_D[4], v[4];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite3_D( x - m_X[i], m_X[i+1] - m_X[i], u_D );
    Hermite3  ( y - m_Y[j], m_Y[j+1] - m_Y[j], v   );
    load( i, j, bili3 );
    return bilinear3( u_D, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dy( real_type x, real_type y ) const {
    real_type bili3[4][4], u[4], v_D[4];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite3  ( x - m_X[i], m_X[i+1] - m_X[i], u   );
    Hermite3_D( y - m_Y[j], m_Y[j+1] - m_Y[j], v_D );
    load( i, j, bili3 );
    return bilinear3( u, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dxy( real_type x, real_type y ) const {
    real_type bili3[4][4], u_D[4], v_D[4];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite3_D( x - m_X[i], m_X[i+1] - m_X[i], u_D );
    Hermite3_D( y - m_Y[j], m_Y[j+1] - m_Y[j], v_D );
    load( i, j, bili3 );
    return bilinear3( u_D, bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dxx( real_type x, real_type y ) const {
    real_type bili3[4][4], u_DD[4], v[4];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite3_DD( x - m_X[i], m_X[i+1] - m_X[i], u_DD );
    Hermite3   ( y - m_Y[j], m_Y[j+1] - m_Y[j], v    );
    load( i, j, bili3 );
    return bilinear3( u_DD, bili3, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiCubicSplineBase::Dyy( real_type x, real_type y ) const {
    real_type bili3[4][4], u[4], v_DD[4];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite3   ( x - m_X[i], m_X[i+1] - m_X[i], u    );
    Hermite3_DD( y - m_Y[j], m_Y[j+1] - m_Y[j], v_DD );
    load( i, j, bili3 );
    return bilinear3( u, bili3, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSplineBase::D( real_type x, real_type y, real_type d[3] ) const {
    real_type bili3[4][4], u[4], u_D[4], v[4], v_D[4];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite3   ( x - m_X[i], m_X[i+1] - m_X[i], u    );
    Hermite3_D ( x - m_X[i], m_X[i+1] - m_X[i], u_D  );
    Hermite3   ( y - m_Y[j], m_Y[j+1] - m_Y[j], v    );
    Hermite3_D ( y - m_Y[j], m_Y[j+1] - m_Y[j], v_D  );
    load( i, j, bili3 );
    d[0] = bilinear3( u,   bili3, v   );
    d[1] = bilinear3( u_D, bili3, v   );
    d[2] = bilinear3( u,   bili3, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSplineBase::DD( real_type x, real_type y, real_type d[6] ) const {
    real_type bili3[4][4], u[4], u_D[4], u_DD[4], v[4], v_D[4], v_DD[4];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite3   ( x - m_X[i], m_X[i+1] - m_X[i], u    );
    Hermite3_D ( x - m_X[i], m_X[i+1] - m_X[i], u_D  );
    Hermite3_DD( x - m_X[i], m_X[i+1] - m_X[i], u_DD );
    Hermite3   ( y - m_Y[j], m_Y[j+1] - m_Y[j], v    );
    Hermite3_D ( y - m_Y[j], m_Y[j+1] - m_Y[j], v_D  );
    Hermite3_DD( y - m_Y[j], m_Y[j+1] - m_Y[j], v_DD );
    load( i, j, bili3 );
    d[0] = bilinear3( u,    bili3, v   );
    d[1] = bilinear3( u_D,  bili3, v    );
    d[2] = bilinear3( u,    bili3, v_D  );
    d[3] = bilinear3( u_DD, bili3, v    );
    d[4] = bilinear3( u_D,  bili3, v_D  );
    d[5] = bilinear3( u,    bili3, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  BiQuinticSplineBase::load(
    integer const i,
    integer const j,
    real_type     bili5[6][6]
  ) const {

    integer const i00 { ipos_C(i,j) };
    integer const i01 { ipos_C(i,j+1) };
    integer const i10 { ipos_C(i+1,j) };
    integer const i11 { ipos_C(i+1,j+1) };

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

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::eval( real_type x, real_type y ) const {
    real_type bili5[6][6], u[6], v[6];
    integer const i{ search_x( x ) };
    integer const j{ search_y( y ) };
    Hermite5( x - m_X[i], m_X[i+1] - m_X[i], u );
    Hermite5( y - m_Y[j], m_Y[j+1] - m_Y[j], v );
    load( i, j, bili5 );
    return bilinear5( u, bili5, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dx( real_type x, real_type y ) const {
    real_type bili5[6][6], u_D[6], v[6];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite5_D( x - m_X[i], m_X[i+1] - m_X[i], u_D );
    Hermite5  ( y - m_Y[j], m_Y[j+1] - m_Y[j], v   );
    load( i, j, bili5 );
    return bilinear5( u_D, bili5, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dy( real_type x, real_type y ) const {
    real_type bili5[6][6], u[6], v_D[6];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite5  ( x - m_X[i], m_X[i+1] - m_X[i], u   );
    Hermite5_D( y - m_Y[j], m_Y[j+1] - m_Y[j], v_D );
    load( i, j, bili5 );
    return bilinear5( u, bili5, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dxy( real_type x, real_type y ) const {
    real_type bili5[6][6], u_D[6], v_D[6];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite5_D( x - m_X[i], m_X[i+1] - m_X[i], u_D );
    Hermite5_D( y - m_Y[j], m_Y[j+1] - m_Y[j], v_D );
    load( i, j, bili5 );
    return bilinear5( u_D, bili5, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dxx( real_type x, real_type y ) const {
    real_type bili5[6][6], u_DD[6], v[6];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite5_DD( x - m_X[i], m_X[i+1] - m_X[i], u_DD );
    Hermite5   ( y - m_Y[j], m_Y[j+1] - m_Y[j], v    );
    load( i, j, bili5 );
    return bilinear5( u_DD, bili5, v );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiQuinticSplineBase::Dyy( real_type x, real_type y ) const {
    real_type bili5[6][6], u[6], v_DD[6];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite5   ( x - m_X[i], m_X[i+1] - m_X[i], u    );
    Hermite5_DD( y - m_Y[j], m_Y[j+1] - m_Y[j], v_DD );
    load( i, j, bili5 );
    return bilinear5( u, bili5, v_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSplineBase::D( real_type x, real_type y, real_type d[3] ) const {
    real_type bili5[6][6], u[6], u_D[6], v[6], v_D[6];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite5   ( x - m_X[i], m_X[i+1] - m_X[i], u    );
    Hermite5_D ( x - m_X[i], m_X[i+1] - m_X[i], u_D  );
    Hermite5   ( y - m_Y[j], m_Y[j+1] - m_Y[j], v    );
    Hermite5_D ( y - m_Y[j], m_Y[j+1] - m_Y[j], v_D  );
    load( i, j, bili5 );
    d[0] = bilinear5( u, bili5, v );
    d[1] = bilinear5( u_D, bili5, v );
    d[2] = bilinear5( u, bili5, v_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSplineBase::DD( real_type x, real_type y, real_type d[6] ) const {
    real_type bili5[6][6], u[6], u_D[6], u_DD[6], v[6], v_D[6], v_DD[6];
    integer const i{search_x( x )};
    integer const j{search_y( y )};
    Hermite5   ( x - m_X[i], m_X[i+1] - m_X[i], u    );
    Hermite5_D ( x - m_X[i], m_X[i+1] - m_X[i], u_D  );
    Hermite5_DD( x - m_X[i], m_X[i+1] - m_X[i], u_DD );
    Hermite5   ( y - m_Y[j], m_Y[j+1] - m_Y[j], v    );
    Hermite5_D ( y - m_Y[j], m_Y[j+1] - m_Y[j], v_D  );
    Hermite5_DD( y - m_Y[j], m_Y[j+1] - m_Y[j], v_DD );
    load( i, j, bili5 );
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
    for ( integer i{1}; i < m_nx; ++i ) s << ", " << m_X[i];
    s << " ]\nY = [ " << m_Y[0];
    for ( integer i{1}; i < m_ny; ++i ) s << ", " << m_Y[i];
    s << " ]\nZ = [\n";
    for ( integer j{0}; j < m_ny; ++j ) {
      s << "  [ " << z_node(0,j);
      for ( integer i{1}; i < m_nx; ++i ) s << ", " << z_node(i,j);
      s << " ]\n";
    }
    s << "\n];\n";
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  using GC_namespace::GC_type;
  using GC_namespace::mat_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Set up a spline surface using a `GenericContainer`
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
    string const where{ fmt::format("SplineSurf[{}]::setup( gc ):", m_name ) };
    GenericContainer const & gc_x { gc("xdata",where) };
    GenericContainer const & gc_y { gc("ydata",where) };
    GenericContainer const & gc_z { gc("zdata",where) };

    m_nx = static_cast<integer>(gc_x.get_num_elements());
    m_ny = static_cast<integer>(gc_y.get_num_elements());
    m_mem.reallocate( (m_nx+1)*(m_ny+1) );
    m_X = m_mem( m_nx );
    m_Y = m_mem( m_ny );
    m_Z = m_mem( m_nx*m_ny );

    for ( integer i{0}; i < m_nx; ++i ) m_X[i] = gc_x.get_number_at(i);
    for ( integer i{0}; i < m_ny; ++i ) m_Y[i] = gc_y.get_number_at(i);

    bool fortran_storage { gc.get_map_bool( "fortran_storage", where ) };
    bool transposed      { gc.get_map_bool( "transposed",      where ) };

    /*
    //     +------+
    //  j ny      | (xi,yj)
    //     +  nx  +
    //         i
    */

    // cosa mi aspetto in lettura
    integer NR, NC;
    if ( transposed ) { NC = m_ny; NR = m_nx; }
    else              { NC = m_nx; NR = m_ny; }
    integer const LD{ fortran_storage ? NR : NC };

    if ( GC_type::MAT_REAL    == gc_z.get_type() ||
         GC_type::MAT_INTEGER == gc_z.get_type() ||
         GC_type::MAT_LONG    == gc_z.get_type() ) {

      UTILS_ASSERT(
        static_cast<unsigned>(NR) == gc_z.num_rows() && static_cast<unsigned>(NC) == gc_z.num_cols(),
        "{}, field `zdata` is a matrix expected to be of size {} x {}, found: {} x {}\n",
        where, NR, NC, gc_z.num_rows(), gc_z.num_cols()
      );

      if ( GC_type::MAT_REAL == gc_z.get_type() ) {
        load_Z( gc_z.get_mat_real().data(), LD, fortran_storage, transposed );
      } else {
        GenericContainer::mat_real_type z_tmp;
        gc_z.copyto_mat_real( z_tmp );
        load_Z( z_tmp.data(), LD, fortran_storage, transposed );  
      }

    } else if ( GC_type::VEC_REAL    == gc_z.get_type() || 
                GC_type::VEC_INTEGER == gc_z.get_type() || 
                GC_type::VEC_LONG    == gc_z.get_type() ) {

      integer nz  { static_cast<integer>( gc_z.get_num_elements() ) };
      integer nxy { m_nx * m_ny };
      UTILS_ASSERT(
        nz == nxy,
        "{}, field `zdata` expected to be of size {} = {}x{}, found: `{}`\n",
        where, nxy, m_nx, m_ny, nz
      );

      for ( integer k{0}; k < nxy; ++k ) {
        integer i, j;
        real_type const v{ gc_z.get_number_at( k ) };
        if ( fortran_storage ) { i = k % NR; j = k / NR; }
        else                   { j = k % NC; i = k / NC; }
        if ( transposed ) z_node_ref(j,i) = v;
        else              z_node_ref(i,j) = v;
      }

    } else if ( GC_type::VECTOR == gc_z.get_type() ) {

      vector_type const & data{ gc_z.get_vector() };
      vec_real_type tmp;
      UTILS_ASSERT(
        static_cast<size_t>(NR) == data.size(),
        "{}, field `zdata` (vector of vector) expected of size {} found of size {}\n",
        where, NR, data.size()
      );
      for ( integer j{0}; j < NC; ++j ) {
        GenericContainer const & row{ data[j] };
        string const msg1{ fmt::format( "{}, reading row {}\n", where, j ) };
        row.copyto_vec_real( tmp, msg1 );
        UTILS_ASSERT(
          static_cast<size_t>(NC) == tmp.size(),
          "{}, row {}-th of size {}, expected {}\n",
          where, j, tmp.size(), NC
        );
        if ( transposed ) {
          for ( integer i{0}; i < NC; ++i ) z_node_ref(j,i) = tmp[ i ];
        } else {
          for ( integer i{0}; i < NC; ++i ) z_node_ref(i,j) = tmp[ i ];
        }
      }

    } else {

      UTILS_ERROR(
        "{}, field `zdata` expected to be of type"
        " `mat_real_type` or  `vec_real_type` or `vector_type` found: `{}`\n",
        where, gc_z.get_type_name()
      );

    }

    make_spline();
  }

  void
  SplineSurf::setup( string const & file_name ) {
    GenericContainer gc;
    UTILS_ASSERT( gc.from_file( file_name ), "Spline::setup( '{}' ) failed to read\n", file_name );
    setup( gc );
  }

}
