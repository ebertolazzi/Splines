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

  //!
  //! This subroutine estimates three partial derivatives, zx, zy, and
  //! zxy, of a bivariate function, z(x,y), on a rectangular grid in
  //! the x-y plane.  It is based on the revised Akima method that has
  //! the accuracy of a bicubic polynomial.
  //!
  void
  Akima2Dspline::make_spline() {
  
    UTILS_ASSERT(
      m_nx >= 2 && m_ny >= 2,
      "Akima2Dspline::make_spline, nx={}, ny={} mus be nx >= 2 and ny >= 2", m_nx, m_ny
    );

    integer const nn{ m_nx*m_ny };
    integer const mm{ std::max( m_nx, m_ny ) };
    m_mem_bicubic.reallocate( 3*nn );
    m_DX  = m_mem_bicubic( nn );
    m_DY  = m_mem_bicubic( nn );
    m_DXY = m_mem_bicubic( nn );
    
    Malloc_real m_mem("Akima2Dspline::work memory");
    m_mem.allocate( 3*mm );
    
    real_type * m  { m_mem( mm ) };
    real_type * Z  { m_mem( mm ) };
    real_type * Zp { m_mem( mm ) };

    #define DEBUG_AKIMA
    #ifdef DEBUG_AKIMA
    string msg;
    #endif
    
    for ( integer j{0}; j < m_ny; ++j ) {
      for ( integer i{0}; i < m_nx; ++i ) Z[i] = z_node(i,j);

      #ifdef DEBUG_AKIMA
      msg = fmt::format("Akima2Dspline::make_spline Z1 {} in [0,{})", j, m_ny );
      Utils::check_NaN( Z, msg, m_nx, __LINE__, __FILE__ );
      #endif

      Akima_build( m_X, Z, Zp, m, m_nx );

      #ifdef DEBUG_AKIMA
      msg = fmt::format("Akima2Dspline::make_spline Zp1 {} in [0,{})", j, m_ny );
      Utils::check_NaN( Zp, msg, m_nx, __LINE__, __FILE__ );
      #endif

      for ( integer i{0}; i < m_nx; ++i ) Dx_node_ref(i,j) = Zp[i];
    }

    for ( integer i{0}; i < m_nx; ++i ) {
      for ( integer j{0}; j < m_ny; ++j ) Z[j] = z_node(i,j);

      #ifdef DEBUG_AKIMA
      msg = fmt::format("Akima2Dspline::make_spline Z2 {} in [0,{})", i, m_nx );
      Utils::check_NaN( Z, msg, m_ny, __LINE__, __FILE__ );
      #endif

      Akima_build( m_Y, Z, Zp, m, m_ny );

      #ifdef DEBUG_AKIMA
      msg = fmt::format("Akima2Dspline::make_spline Zp2 {} in [0,{})", i, m_nx );
      Utils::check_NaN( Zp, msg, m_ny, __LINE__, __FILE__ );
      #endif

      for ( integer j{0}; j < m_ny; ++j ) Dy_node_ref(i,j) = Zp[j];
    }

    for ( integer j{0}; j < m_ny; ++j ) {
      for ( integer i{0}; i < m_nx; ++i ) Z[i] = Dy_node(i,j);

      #ifdef DEBUG_AKIMA
      msg = fmt::format("Akima2Dspline::make_spline Zp3 {} in [0,{})", j, m_ny );
      Utils::check_NaN( Z, msg, m_nx, __LINE__, __FILE__ );
      #endif

      Akima_build( m_X, Z, Zp, m, m_nx );

      #ifdef DEBUG_AKIMA
      msg = fmt::format("Akima2Dspline::make_spline Zp3 {} in [0,{})", j, m_ny );
      Utils::check_NaN( Zp, msg, m_nx, __LINE__, __FILE__ );
      #endif

      for ( integer i{0}; i < m_nx; ++i ) Dxy_node_ref(i,j) = Zp[i];
    }

    auto minmod = [] ( real_type a, real_type b ) -> real_type {
      if ( a*b <= 0 ) return 0;
      if ( a > 0    ) return std::min(a,b);
      return std::max(a,b);
    };

    for ( integer i{0}; i < m_nx; ++i ) {

      for ( integer j{0}; j < m_ny; ++j ) Z[j] = Dx_node(i,j);

      #ifdef DEBUG_AKIMA
      msg = fmt::format("Akima2Dspline::make_spline Z4 {} in [0,{})", i, m_nx );
      Utils::check_NaN( Z, msg, m_ny, __LINE__, __FILE__ );
      #endif

      Akima_build( m_Y, Z, Zp, m, m_ny );

      #ifdef DEBUG_AKIMA
      msg = fmt::format("Akima2Dspline::make_spline Zp4 {} in [0,{})", i, m_nx );
      Utils::check_NaN( Zp, msg, m_ny, __LINE__, __FILE__ );
      #endif

      for ( integer j{0}; j < m_ny; ++j ) Dxy_node_ref(i,j) = minmod( Dxy_node_ref(i,j), Zp[j] );
    }

    Utils::check_NaN( m_DX,  "Akima2Dspline::make_spline DX ",  nn, __LINE__, __FILE__ );
    Utils::check_NaN( m_DY,  "Akima2Dspline::make_spline DY ",  nn, __LINE__, __FILE__ );
    Utils::check_NaN( m_DXY, "Akima2Dspline::make_spline DXY ", nn, __LINE__, __FILE__ );
  }

  void
  Akima2Dspline::write_to_stream( ostream_type & s ) const {
    fmt::print( s, "Nx={} Ny={}\n", m_nx, m_ny );
    for ( integer i{1}; i < m_nx; ++i ) {
      for ( integer j{1}; j < m_ny; ++j ) {
        integer const i00 { ipos_C(i-1,j-1) };
        integer const i10 { ipos_C(i,j-1) };
        integer const i01 { ipos_C(i-1,j) };
        integer const i11 { ipos_C(i,j) };
        fmt::print( s,
          "patch ({},{)\n"
          "  DX   = {:<12.4}  DY   = {:<12.4}\n"
          "  Z00  = {:<12.4}  Z01  = {:<12.4}  Z10  = {:<12.4}  Z11  = {:<12.4}\n"
          "  Dx00 = {:<12.4}  Dx01 = {:<12.4}  Dx10 = {:<12.4}  Dx11 = {:<12.4}\n"
          "  Dy00 = {:<12.4}  Dy01 = {:<12.4}  Dy10 = {:<12.4}  Dy11 = {:<12.4}\n",
          i, j,
          m_X[i] - m_X[i-1],
          m_Y[j] - m_Y[j-1],
          m_Z[i00],  m_Z[i01],  m_Z[i10],  m_Z[i11],
          m_DX[i00], m_DX[i01], m_DX[i10], m_DX[i11],
          m_DY[i00], m_DY[i01], m_DY[i10], m_DY[i11]
        );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char const *
  Akima2Dspline::type_name() const
  { return "Akima2D"; }

}
