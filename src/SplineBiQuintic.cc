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
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namspace
#endif

namespace Splines {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSpline::make_spline() {

    integer const dim{ m_nx*m_ny };
    m_mem_biquintic.reallocate( 8*dim );
    m_DX    = m_mem_biquintic( dim );
    m_DY    = m_mem_biquintic( dim );
    m_DXY   = m_mem_biquintic( dim );
    m_DXX   = m_mem_biquintic( dim );
    m_DYY   = m_mem_biquintic( dim );
    m_DXYY  = m_mem_biquintic( dim );
    m_DXXY  = m_mem_biquintic( dim );
    m_DXXYY = m_mem_biquintic( dim );

    auto minmod = [] ( real_type a, real_type b ) -> real_type {
      if ( a*b <= 0 ) return 0;
      if ( a > 0    ) return std::min(a,b);
      return std::max(a,b);
    };

    PchipSpline sp;

    // calcolo derivate DX, DY, DXY
    for ( integer j{0}; j < m_ny; ++j ) {
      sp.build( m_X, 1, &z_node_ref(0,j), m_ny, m_nx );
      for ( integer i{0}; i < m_nx; ++i ) Dx_node_ref(i,j) = sp.yp_node(i);
    }
    for ( integer i{0}; i < m_nx; ++i ) {
      sp.build( m_Y, 1, &z_node_ref(i,0), 1, m_ny );
      for ( integer j{0}; j < m_ny; ++j ) Dy_node_ref(i,j) = sp.yp_node(j);
    }

    // mixed
    for ( integer j{0}; j < m_ny; ++j ) {
      sp.build( m_X, 1, &Dx_node_ref(0,j), m_ny, m_nx );
      for ( integer i{0}; i < m_nx; ++i ) Dxy_node_ref(i,j) = sp.yp_node(i);
    }
    for ( integer i{0}; i < m_nx; ++i ) {
      sp.build( m_Y, 1, &Dy_node_ref(i,0), 1, m_ny );
      for ( integer j{0}; j < m_ny; ++j ) Dxy_node_ref(i,j) = minmod( Dxy_node_ref(i,j), sp.yp_node(j) );
    }

    // calcolo derivate DXX, DYY, DXXY, DXYY, DXXYY
    for ( integer j{0}; j < m_ny; ++j ) {
      sp.build( m_X, 1, &Dx_node_ref(0,j), m_ny, m_nx );
      for ( integer i{0}; i < m_nx; ++i ) Dxx_node_ref(i,j) = sp.yp_node(i);
    }
    for ( integer i{0}; i < m_nx; ++i ) {
      sp.build( m_Y, 1, &Dy_node_ref(i,0), 1, m_ny );
      for ( integer j{0}; j < m_ny; ++j ) Dyy_node_ref(i,j) = sp.yp_node(j);
    }
    for ( integer j{0}; j < m_ny; ++j ) {
      sp.build( m_X, 1, &Dxy_node_ref(0,j), m_ny, m_nx );
      for ( integer i{0}; i < m_nx; ++i ) Dxxy_node_ref(i,j) = sp.yp_node(i);
    }
    for ( integer i{0}; i < m_nx; ++i ) {
      sp.build( m_Y, 1, &Dxy_node_ref(i,0), 1, m_ny );
      for ( integer j{0}; j < m_ny; ++j ) Dxyy_node_ref(i,j) = sp.yp_node(j);
    }
    // mixed
    for ( integer j{0}; j < m_ny; ++j ) {
      sp.build( m_X, 1, &Dxyy_node_ref(0,j), m_ny, m_nx );
      for ( integer i{0}; i < m_nx; ++i ) Dxxyy_node_ref(i,j) = sp.yp_node(i);
    }
    for ( integer i{0}; i < m_nx; ++i ) {
      sp.build( m_Y, 1, &Dxxy_node_ref(i,0), 1, m_ny );
      for ( integer j{0}; j < m_ny; ++j ) Dxxy_node_ref(i,j) = minmod( Dxxy_node_ref(i,j), sp.yp_node(j) );
    }

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSpline::write_to_stream( ostream_type & s ) const {
    fmt::print( s, "Nx = {} Ny = {}\n", m_nx, m_ny );
    for ( integer i{1}; i < m_nx; ++i ) {
      real_type dx{ m_X[i]-m_X[i-1] };
      for ( integer j{1}; j < m_ny; ++j ) {
        integer   const i00 { ipos_C(i-1,j-1) };
        integer   const i10 { ipos_C(i,j-1) };
        integer   const i01 { ipos_C(i-1,j) };
        integer   const i11 { ipos_C(i,j) };
        real_type const dy  { m_Y[j]-m_Y[j-1] };
        fmt::print( s,
          "patch ({},{})\n"
          "  DX    = {:<12.4}  DY    = {:<12.4}\n"
          "  Z00   = {:<12.4}  Z10   = {:<12.4}\n"
          "  Z01   = {:<12.4}  Z11   = {:<12.4}\n"
          "  Dx00  = {:<12.4}  Dx10  = {:<12.4}\n"
          "  Dx01  = {:<12.4}  Dx11  = {:<12.4}\n"
          "  Dy00  = {:<12.4}  Dy10  = {:<12.4}\n"
          "  Dy01  = {:<12.4}  Dy11  = {:<12.4}\n"
          "  Dxy00 = {:<12.4}  Dxy10 = {:<12.4}\n"
          "  Dxy01 = {:<12.4}  Dxy11 = {:<12.4}\n",
          i, j, dx, dy,
          m_Z[i00],   m_Z[i10],   m_Z[i01],   m_Z[i11],
          m_DX[i00],  m_DX[i10],  m_DX[i01],  m_DX[i11],
          m_DY[i00],  m_DY[i10],  m_DY[i01],  m_DY[i11],
          m_DXY[i00], m_DXY[i10], m_DXY[i01], m_DXY[i11]
        );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char const *
  BiQuinticSpline::type_name() const
  { return "BiQuintic"; }

}
