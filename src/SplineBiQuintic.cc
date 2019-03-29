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
    this->DX.resize(this->Z.size());
    this->DY.resize(this->Z.size());
    this->DXY.resize(this->Z.size());
    this->DXX.resize(this->Z.size());
    this->DYY.resize(this->Z.size());
    this->DXYY.resize(this->Z.size());
    this->DXXY.resize(this->Z.size());
    this->DXXYY.resize(this->Z.size());
    // calcolo derivate
    integer nx = integer(X.size());
    integer ny = integer(Y.size());
    QuinticSpline sp, sp1;
    //
    for ( integer j = 0; j < ny; ++j ) {
      sp.build( &this->X.front(), 1, &this->Z[size_t(this->ipos_C(0,j))], ny, nx );
      for ( integer i = 0; i < nx; ++i ) {
        this->DX[size_t(this->ipos_C(i,j))]  = sp.ypNode(i);
        this->DXX[size_t(this->ipos_C(i,j))] = sp.yppNode(i);
      }
    }
    for ( integer i = 0; i < nx; ++i ) {
      sp.build( &this->Y.front(), 1, &this->Z[size_t(this->ipos_C(i,0))], 1, ny );
      for ( integer j = 0; j < ny; ++j ) {
        this->DY[size_t(this->ipos_C(i,j))]  = sp.ypNode(j);
        this->DYY[size_t(this->ipos_C(i,j))] = sp.yppNode(j);
      }
    }
    // interpolate derivative
    for ( integer i = 0; i < nx; ++i ) {
      sp.build( &this->Y.front(), 1, &this->DX[size_t(this->ipos_C(i,0))], 1, ny );
      sp1.build( &this->Y.front(), 1, &this->DXX[size_t(this->ipos_C(i,0))], 1, ny );
      for ( integer j = 0; j < ny; ++j ) {
        this->DXY[size_t(this->ipos_C(i,j))]   = sp.ypNode(j);
        this->DXYY[size_t(this->ipos_C(i,j))]  = sp.yppNode(j);
        this->DXXY[size_t(this->ipos_C(i,j))]  = sp1.ypNode(j);
        this->DXXYY[size_t(this->ipos_C(i,j))] = sp1.yppNode(j);
      }
    }
    // interpolate derivative again
    for ( integer j = 0; j < ny; ++j ) {
      sp.build( &this->X.front(), 1, &this->DY[size_t(this->ipos_C(0,j))], ny, nx );
      sp1.build( &this->X.front(), 1, &this->DYY[size_t(this->ipos_C(0,j))], ny, nx );
      for ( integer i = 0; i < nx; ++i ) {
        this->DXY[size_t(this->ipos_C(i,j))]   += sp.ypNode(i);   this->DXY[size_t(this->ipos_C(i,j))]   /= 2;
        this->DXXY[size_t(this->ipos_C(i,j))]  += sp.yppNode(i);  this->DXXY[size_t(this->ipos_C(i,j))]  /= 2;
        this->DXYY[size_t(this->ipos_C(i,j))]  += sp1.ypNode(i);  this->DXYY[size_t(this->ipos_C(i,j))]  /= 2;
        this->DXXYY[size_t(this->ipos_C(i,j))] += sp1.yppNode(i); this->DXXYY[size_t(this->ipos_C(i,j))] /= 2;
      }
    }

    //std::fill( DXY.begin(), DXY.end(), 0 );
    //std::fill( DXX.begin(), DXX.end(), 0 );
    //std::fill( DYY.begin(), DYY.end(), 0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiQuinticSpline::writeToStream( ostream_type & s ) const {
    integer ny = integer(this->Y.size());
    s << "Nx = " << this->X.size() << " Ny = " << this->Y.size() << '\n';
    for ( integer i = 1; i < integer(this->X.size()); ++i ) {
      for ( integer j = 1; j < integer(this->Y.size()); ++j ) {
        size_t i00 = size_t(this->ipos_C(i-1,j-1,ny));
        size_t i10 = size_t(this->ipos_C(i,j-1,ny));
        size_t i01 = size_t(this->ipos_C(i-1,j,ny));
        size_t i11 = size_t(this->ipos_C(i,j,ny));
        s << "patch (" << i << "," << j
          << ")\n DX = "  << setw(10) << left << this->X[size_t(i)]-this->X[size_t(i-1)]
          <<    " DY = "  << setw(10) << left << this->Y[size_t(j)]-this->Y[size_t(j-1)]
          << "\n Z00  = " << setw(10) << left << this->Z[i00]
          <<   " Z01  = " << setw(10) << left << this->Z[i01]
          <<   " Z10  = " << setw(10) << left << this->Z[i10]
          <<   " Z11  = " << setw(10) << left << this->Z[i11]
          << "\n Dx00 = " << setw(10) << left << this->DX[i00]
          <<   " Dx01 = " << setw(10) << left << this->DX[i01]
          <<   " Dx10 = " << setw(10) << left << this->DX[i10]
          <<   " Dx10 = " << setw(10) << left << this->DX[i11]
          << "\n Dy00 = " << setw(10) << left << this->DY[i00]
          <<   " Dy01 = " << setw(10) << left << this->DY[i01]
          <<   " Dy10 = " << setw(10) << left << this->DY[i10]
          <<   " Dy11 = " << setw(10) << left << this->DY[i11]
          << '\n';
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char const *
  BiQuinticSpline::type_name() const
  { return "BiQuintic"; }

}
