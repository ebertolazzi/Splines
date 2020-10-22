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

/**
 *
 */

namespace Splines {

  using namespace std; // load standard namspace

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSpline::makeSpline() {
    size_t nn = size_t(m_nx*m_ny);
    m_mem_bicubic.allocate( 3*nn );
    m_DX  = m_mem_bicubic( nn );
    m_DY  = m_mem_bicubic( nn );
    m_DXY = m_mem_bicubic( nn );

    // calcolo derivate
    PchipSpline sp;
    for ( integer j = 0; j < m_ny; ++j ) {
      sp.build( m_X, 1, &m_Z[size_t(this->ipos_C(0,j))], m_ny, m_nx );
      for ( integer i = 0; i < m_nx; ++i )
        m_DX[size_t(this->ipos_C(i,j))] = sp.ypNode(i);
    }
    for ( integer i = 0; i < m_nx; ++i ) {
      sp.build( m_Y, 1, &m_Z[size_t(this->ipos_C(i,0))], 1, m_ny );
      for ( integer j = 0; j < m_ny; ++j )
        m_DY[size_t(this->ipos_C(i,j))] = sp.ypNode(j);
    }
    std::fill_n( m_DXY, nn, 0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSpline::writeToStream( ostream_type & s ) const {
    s << "Nx = " << m_nx << " Ny = " << m_ny << '\n';
    for ( integer i = 1; i < m_nx; ++i ) {
      for ( integer j = 1; j < m_ny; ++j ) {
        size_t i00 = size_t(this->ipos_C(i-1,j-1,m_ny));
        size_t i10 = size_t(this->ipos_C(i,j-1,m_ny));
        size_t i01 = size_t(this->ipos_C(i-1,j,m_ny));
        size_t i11 = size_t(this->ipos_C(i,j,m_ny));
        s << "patch (" << i << "," << j
          << ")\n DX = "  << setw(10) << left << m_X[size_t(i)]-m_X[size_t(i-1)]
          <<    " DY = "  << setw(10) << left << m_Y[size_t(j)]-m_Y[size_t(j-1)]
          << "\n Z00  = " << setw(10) << left << m_Z[i00]
          <<   " Z01  = " << setw(10) << left << m_Z[i01]
          <<   " Z10  = " << setw(10) << left << m_Z[i10]
          <<   " Z11  = " << setw(10) << left << m_Z[i11]
          << "\n Dx00 = " << setw(10) << left << m_DX[i00]
          <<   " Dx01 = " << setw(10) << left << m_DX[i01]
          <<   " Dx10 = " << setw(10) << left << m_DX[i10]
          <<   " Dx11 = " << setw(10) << left << m_DX[i11]
          << "\n Dy00 = " << setw(10) << left << m_DY[i00]
          <<   " Dy01 = " << setw(10) << left << m_DY[i01]
          <<   " Dy10 = " << setw(10) << left << m_DY[i10]
          <<   " Dy11 = " << setw(10) << left << m_DY[i11]
          << '\n';
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char const *
  BiCubicSpline::type_name() const
  { return "BiCubic"; }

}
