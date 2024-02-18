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
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namspace
#endif

namespace Splines {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSpline::make_spline() {
    size_t nn = size_t(m_nx*m_ny);
    m_mem_bicubic.reallocate( 3*nn );
    m_DX  = m_mem_bicubic( nn );
    m_DY  = m_mem_bicubic( nn );
    m_DXY = m_mem_bicubic( nn );

    // calcolo derivate
    PchipSpline sp;
    for ( integer j = 0; j < m_ny; ++j ) {
      sp.build( m_X, 1, &m_Z[size_t(this->ipos_C(0,j))], m_ny, m_nx );
      for ( integer i = 0; i < m_nx; ++i )
        m_DX[size_t(this->ipos_C(i,j))] = sp.yp_node(i);
    }
    for ( integer i = 0; i < m_nx; ++i ) {
      sp.build( m_Y, 1, &m_Z[size_t(this->ipos_C(i,0))], 1, m_ny );
      for ( integer j = 0; j < m_ny; ++j )
        m_DY[size_t(this->ipos_C(i,j))] = sp.yp_node(j);
    }
    std::fill_n( m_DXY, nn, 0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSpline::write_to_stream( ostream_type & s ) const {
    fmt::print( "Nx = {} Ny = {}\n", m_nx, m_ny );
    for ( integer i = 1; i < m_nx; ++i ) {
      for ( integer j = 1; j < m_ny; ++j ) {
        size_t i00 = size_t(this->ipos_C(i-1,j-1));
        size_t i10 = size_t(this->ipos_C(i,j-1));
        size_t i01 = size_t(this->ipos_C(i-1,j));
        size_t i11 = size_t(this->ipos_C(i,j));
        fmt::print( s,
          "patch ({},{})\n"
          "DX   = {:<12.4}  DY   = {:<12.4}\n"
          "Z00  = {:<12.4}  Z01  = {:<12.4}  Z10  = {:<12.4}  Z11  = {:<12.4}\n"
          "Dx00 = {:<12.4}  Dx01 = {:<12.4}  Dx10 = {:<12.4}  Dx11 = {:<12.4}\n"
          "Dy00 = {:<12.4}  Dy01 = {:<12.4}  Dy10 = {:<12.4}  Dy11 = {:<12.4}\n",
          i, j,
          m_X[size_t(i)]-m_X[size_t(i-1)],
          m_Y[size_t(j)]-m_Y[size_t(j-1)],
          m_Z[i00], m_Z[i01], m_Z[i10], m_Z[i11],
          m_DX[i00], m_DX[i01], m_DX[i10], m_DX[i11],
          m_DY[i00], m_DY[i01], m_DY[i10], m_DY[i11]
        );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char const *
  BiCubicSpline::type_name() const
  { return "BiCubic"; }

}
