/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 1998                                                      |
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

  using namespace std ; // load standard namspace
  
  void
  BiCubicSpline::makeSpline() {
    DX.resize(Z.size()) ;
    DY.resize(Z.size()) ;
    DXY.resize(Z.size()) ;
    // calcolo derivate
    sizeType nx = sizeType(X.size()) ;
    sizeType ny = sizeType(Y.size()) ;
    PchipSpline sp ;
    for ( sizeType j = 0 ; j < ny ; ++j ) {
      sp.build( &X.front(), 1, &Z[ipos_C(0,j)], ny, nx ) ;
      for ( sizeType i = 0 ; i < nx ; ++i ) DX[ipos_C(i,j)] = sp.ypNode(i) ;
    }
    for ( sizeType i = 0 ; i < nx ; ++i ) {
      sp.build( &Y.front(), 1, &Z[ipos_C(i,0)], 1, ny ) ;
      for ( sizeType j = 0 ; j < ny ; ++j ) DY[ipos_C(i,j)] = sp.ypNode(j) ;
    }
    std::fill( DXY.begin(), DXY.end(), 0 ) ;
  }

  void
  BiCubicSpline::writeToStream( ostream & s ) const {
    s << "Nx = " << X.size() << " Ny = " << Y.size() << '\n' ;
    for ( sizeType i = 1 ; i < sizeType(X.size()) ; ++i ) {
      for ( sizeType j = 1 ; j < sizeType(Y.size()) ; ++j ) {
        s << "patch (" << setw(2) << i << "," << setw(2) << j
          << "): DX = " << X[i]-X[i-1] << " DY = " << Y[j]-Y[j-1]
          << "\n Z00 = " << Z[ipos_C(i-1,j-1)]
          << " Z01 = " << Z[ipos_C(i-1,j)]
          << " Z10 = " << Z[ipos_C(i,j-1)]
          << " Z11 = " << Z[ipos_C(i,j)]
          << "\n Dx00 = " << DX[ipos_C(i-1,j-1)]
          << " Dx01 = " << DX[ipos_C(i-1,j)]
          << " Dx10 = " << DX[ipos_C(i,j-1)]
          << " ZDx1 = " << DX[ipos_C(i,j)]
          << "\n Dy00 = " << DY[ipos_C(i-1,j-1)]
          << " Dy01 = " << DY[ipos_C(i-1,j)]
          << " Dy10 = " << DY[ipos_C(i,j-1)]
          << " Dy11 = " << DY[ipos_C(i,j)]
          << '\n' ;
      }
    }
  }

  char const *
  BiCubicSpline::type_name() const
  { return "BiCubic" ; }

}
