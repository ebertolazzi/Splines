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
  BiQuinticSpline::makeSpline() {
    DX.resize(Z.size()) ;
    DY.resize(Z.size()) ;
    DXY.resize(Z.size()) ;
    DXX.resize(Z.size()) ;
    DYY.resize(Z.size()) ;
    DXYY.resize(Z.size()) ;
    DXXY.resize(Z.size()) ;
    DXXYY.resize(Z.size()) ;
    // calcolo derivate
    sizeType nx = sizeType(X.size()) ;
    sizeType ny = sizeType(Y.size()) ;
    QuinticSpline sp, sp1 ;
    //
    for ( sizeType j = 0 ; j < ny ; ++j ) {
      sp.build( &X.front(), 1, &Z[ipos_C(0,j)], ny, nx ) ;
      for ( sizeType i = 0 ; i < nx ; ++i ) {
        DX[ipos_C(i,j)]  = sp.ypNode(i) ;
        DXX[ipos_C(i,j)] = sp.yppNode(i) ;
      }
    }
    for ( sizeType i = 0 ; i < nx ; ++i ) {
      sp.build( &Y.front(), 1, &Z[ipos_C(i,0)], 1, ny ) ;
      for ( sizeType j = 0 ; j < ny ; ++j ) {
        DY[ipos_C(i,j)]  = sp.ypNode(j) ;
        DYY[ipos_C(i,j)] = sp.yppNode(j) ;
      }
    }
    // interpolo le derivate
    for ( sizeType i = 0 ; i < nx ; ++i ) {
      sp.build( &Y.front(), 1, &DX[ipos_C(i,0)], 1, ny ) ;
      sp1.build( &Y.front(), 1, &DXX[ipos_C(i,0)], 1, ny ) ;
      for ( sizeType j = 0 ; j < ny ; ++j ) {
        DXY[ipos_C(i,j)]   = sp.ypNode(j) ;
        DXYY[ipos_C(i,j)]  = sp.yppNode(j) ;
        DXXY[ipos_C(i,j)]  = sp1.ypNode(j) ;
        DXXYY[ipos_C(i,j)] = sp1.yppNode(j) ;
      }
    }
    for ( sizeType j = 0 ; j < ny ; ++j ) {
      sp.build( &X.front(), 1, &DY[ipos_C(0,j)], ny, nx ) ;
      sp1.build( &X.front(), 1, &DYY[ipos_C(0,j)], ny, nx ) ;
      for ( sizeType i = 0 ; i < nx ; ++i ) {
        DXY[ipos_C(i,j)]   += sp.ypNode(i)   ; DXY[ipos_C(i,j)] /= 2 ;
        DXXY[ipos_C(i,j)]  += sp.yppNode(i)  ; DXXY[ipos_C(i,j)] /= 2 ;
        DXYY[ipos_C(i,j)]  += sp1.ypNode(i)  ; DXYY[ipos_C(i,j)] /= 2 ;
        DXXYY[ipos_C(i,j)] += sp1.yppNode(i) ; DXXYY[ipos_C(i,j)] /= 2 ;
      }
    }

    //std::fill( DXY.begin(), DXY.end(), 0 ) ;
    //std::fill( DXX.begin(), DXX.end(), 0 ) ;
    //std::fill( DYY.begin(), DYY.end(), 0 ) ;
  }

  void
  BiQuinticSpline::writeToStream( ostream & s ) const {
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
  BiQuinticSpline::type_name() const
  { return "BiQuintic" ; }

}
