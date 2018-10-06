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
    DX.resize(Z.size());
    DY.resize(Z.size());
    DXY.resize(Z.size());
    // calcolo derivate
    integer nx = integer(X.size());
    integer ny = integer(Y.size());
    PchipSpline sp;
    for ( integer j = 0; j < ny; ++j ) {
      sp.build( &X.front(), 1, &Z[size_t(ipos_C(0,j))], ny, nx );
      for ( integer i = 0; i < nx; ++i )
        DX[size_t(ipos_C(i,j))] = sp.ypNode(i);
    }
    for ( integer i = 0; i < nx; ++i ) {
      sp.build( &Y.front(), 1, &Z[size_t(ipos_C(i,0))], 1, ny );
      for ( integer j = 0; j < ny; ++j )
        DY[size_t(ipos_C(i,j))] = sp.ypNode(j);
    }
    std::fill( DXY.begin(), DXY.end(), 0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiCubicSpline::writeToStream( ostream_type & s ) const {
    integer ny = integer(Y.size());
    s << "Nx = " << X.size() << " Ny = " << Y.size() << '\n';
    for ( integer i = 1; i < integer(X.size()); ++i ) {
      for ( integer j = 1; j < integer(Y.size()); ++j ) {
        size_t i00 = size_t(ipos_C(i-1,j-1,ny));
        size_t i10 = size_t(ipos_C(i,j-1,ny));
        size_t i01 = size_t(ipos_C(i-1,j,ny));
        size_t i11 = size_t(ipos_C(i,j,ny));
        s << "patch (" << i << "," << j
          << ")\n DX = " << setw(10) << left << X[size_t(i)]-X[size_t(i-1)]
          <<    " DY = " << setw(10) << left << Y[size_t(j)]-Y[size_t(j-1)]
          << "\n Z00  = " << setw(10) << left << Z[i00]
          <<   " Z01  = " << setw(10) << left << Z[i01]
          <<   " Z10  = " << setw(10) << left << Z[i10]
          <<   " Z11  = " << setw(10) << left << Z[i11]
          << "\n Dx00 = " << setw(10) << left << DX[i00]
          <<   " Dx01 = " << setw(10) << left << DX[i01]
          <<   " Dx10 = " << setw(10) << left << DX[i10]
          <<   " Dx11 = " << setw(10) << left << DX[i11]
          << "\n Dy00 = " << setw(10) << left << DY[i00]
          <<   " Dy01 = " << setw(10) << left << DY[i01]
          <<   " Dy10 = " << setw(10) << left << DY[i10]
          <<   " Dy11 = " << setw(10) << left << DY[i11]
          << '\n';
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char const *
  BiCubicSpline::type_name() const
  { return "BiCubic"; }

}
