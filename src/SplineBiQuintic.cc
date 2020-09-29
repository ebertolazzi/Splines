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
  BiQuinticSpline::makeSpline() {
    m_DX.resize(m_Z.size());
    m_DY.resize(m_Z.size());
    m_DXY.resize(m_Z.size());
    m_DXX.resize(m_Z.size());
    m_DYY.resize(m_Z.size());
    m_DXYY.resize(m_Z.size());
    m_DXXY.resize(m_Z.size());
    m_DXXYY.resize(m_Z.size());
    // calcolo derivate
    integer nx = integer(m_X.size());
    integer ny = integer(m_Y.size());
    QuinticSpline sp, sp1;
    //
    for ( integer j = 0; j < ny; ++j ) {
      sp.build( &m_X.front(), 1, &m_Z[size_t(this->ipos_C(0,j))], ny, nx );
      for ( integer i = 0; i < nx; ++i ) {
        size_t ij = size_t(this->ipos_C(i,j));
        m_DX[ij]  = sp.ypNode(i);
        m_DXX[ij] = sp.yppNode(i);
      }
    }
    for ( integer i = 0; i < nx; ++i ) {
      sp.build( &m_Y.front(), 1, &m_Z[size_t(this->ipos_C(i,0))], 1, ny );
      for ( integer j = 0; j < ny; ++j ) {
        size_t ij = size_t(this->ipos_C(i,j));
        m_DY[ij]  = sp.ypNode(j);
        m_DYY[ij] = sp.yppNode(j);
      }
    }
    // interpolate derivative
    for ( integer i = 0; i < nx; ++i ) {
      sp.build( &m_Y.front(), 1, &m_DX[size_t(this->ipos_C(i,0))], 1, ny );
      sp1.build( &m_Y.front(), 1, &m_DXX[size_t(this->ipos_C(i,0))], 1, ny );
      for ( integer j = 0; j < ny; ++j ) {
        size_t ij = size_t(this->ipos_C(i,j));
        m_DXY[ij]   = sp.ypNode(j);
        m_DXYY[ij]  = sp.yppNode(j);
        m_DXXY[ij]  = sp1.ypNode(j);
        m_DXXYY[ij] = sp1.yppNode(j);
      }
    }
    // interpolate derivative again
    for ( integer j = 0; j < ny; ++j ) {
      sp.build( &m_X.front(), 1, &m_DY[size_t(this->ipos_C(0,j))], ny, nx );
      sp1.build( &m_X.front(), 1, &m_DYY[size_t(this->ipos_C(0,j))], ny, nx );
      for ( integer i = 0; i < nx; ++i ) {
        size_t ij = size_t(this->ipos_C(i,j));
        m_DXY[ij]   += sp.ypNode(i);   m_DXY[ij]   /= 2;
        m_DXXY[ij]  += sp.yppNode(i);  m_DXXY[ij]  /= 2;
        m_DXYY[ij]  += sp1.ypNode(i);  m_DXYY[ij]  /= 2;
        m_DXXYY[ij] += sp1.yppNode(i); m_DXXYY[ij] /= 2;
      }
    }

    //std::fill( DXY.begin(), DXY.end(), 0 );
    //std::fill( DXX.begin(), DXX.end(), 0 );
    //std::fill( DYY.begin(), DYY.end(), 0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSpline::writeToStream( ostream_type & s ) const {
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
          <<   " Dx10 = " << setw(10) << left << m_DX[i11]
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
  BiQuinticSpline::type_name() const
  { return "BiQuintic"; }

}
