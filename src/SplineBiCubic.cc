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

  BiCubicSplineBase::BiCubicSplineBase( string_view name )
  : SplineSurf( name )
  , m_mem_bicubic( fmt::format("BiCubicSplineBase[{}]",name) )
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSpline::make_spline() {
    integer const nn{ m_nx*m_ny };
    m_mem_bicubic.reallocate( 3*nn );
    m_DX  = m_mem_bicubic( nn );
    m_DY  = m_mem_bicubic( nn );
    m_DXY = m_mem_bicubic( nn );

    make_derivative_x( m_Z, m_DX );
    make_derivative_y( m_Z, m_DY );
    make_derivative_xy( m_DX, m_DY, m_DXY );

    m_search_x.reset();
    m_search_y.reset();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSpline::write_to_stream( ostream_type & s ) const {
    fmt::print( "Nx = {} Ny = {}\n", m_nx, m_ny );
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
  BiCubicSpline::type_name() const
  { return "BiCubic"; }

}
