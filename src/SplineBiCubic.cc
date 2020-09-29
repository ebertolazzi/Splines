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

  using namespace std; // load standard namspace

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSpline::makeSpline() {
    m_DX.resize(m_Z.size());
    m_DY.resize(m_Z.size());
    m_DXY.resize(m_Z.size());
    // calcolo derivate
    integer nx = integer(m_X.size());
    integer ny = integer(m_Y.size());
    PchipSpline sp;
    for ( integer j = 0; j < ny; ++j ) {
      sp.build( &m_X.front(), 1, &m_Z[size_t(this->ipos_C(0,j))], ny, nx );
      for ( integer i = 0; i < nx; ++i )
        m_DX[size_t(this->ipos_C(i,j))] = sp.ypNode(i);
    }
    for ( integer i = 0; i < nx; ++i ) {
      sp.build( &m_Y.front(), 1, &m_Z[size_t(this->ipos_C(i,0))], 1, ny );
      for ( integer j = 0; j < ny; ++j )
        m_DY[size_t(this->ipos_C(i,j))] = sp.ypNode(j);
    }
    std::fill( m_DXY.begin(), m_DXY.end(), 0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSpline::writeToStream( ostream_type & s ) const {
    integer ny = integer(m_Y.size());
    s << "Nx = " << m_X.size() << " Ny = " << m_Y.size() << '\n';
    for ( integer i = 1; i < integer(m_X.size()); ++i ) {
      for ( integer j = 1; j < integer(m_Y.size()); ++j ) {
        size_t i00 = size_t(this->ipos_C(i-1,j-1,ny));
        size_t i10 = size_t(this->ipos_C(i,j-1,ny));
        size_t i01 = size_t(this->ipos_C(i-1,j,ny));
        size_t i11 = size_t(this->ipos_C(i,j,ny));
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
