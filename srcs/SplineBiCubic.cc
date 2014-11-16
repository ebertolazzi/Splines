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
  BiCubicSpline::build ( valueType const x[], sizeType incx,
                         valueType const y[], sizeType incy,
                         valueType const z[], sizeType incz,
                         sizeType nx, sizeType ny ) {
    X.resize(nx) ;
    Y.resize(ny) ;
    Z.resize(nx*ny) ;
    for ( sizeType i = 0 ; i < nx    ; ++i ) X[i] = x[i*incx] ;
    for ( sizeType i = 0 ; i < ny    ; ++i ) Y[i] = y[i*incy] ;
    for ( sizeType i = 0 ; i < nx*ny ; ++i ) Z[i] = z[i*incz] ;
    makeSpline() ;
  }

  void
  BiCubicSpline::build ( VectorOfValues const & x,
                         VectorOfValues const & y,
                         VectorOfValues const & z ) {
    SPLINE_ASSERT( z.size() >= x.size()*y.size(),
                   "in BilinearSpline::build insufficient dimension of Z, expeced >= " << x.size()*y.size() <<
                   " found " << z.size() ) ;
    X.resize(x.size()) ;
    Y.resize(y.size()) ;
    Z.resize(x.size()*y.size()) ;
    std::copy( x.begin(), x.end(), X.begin() ) ;
    std::copy( y.begin(), y.end(), Y.begin() ) ;
    std::copy( z.begin(), z.begin()+Z.size(), Z.begin() ) ;
    makeSpline() ;
  }

  void
  BiCubicSpline::build ( valueType const z[], sizeType incz, sizeType nx, sizeType ny ) {
    X.resize(nx) ;
    Y.resize(ny) ;
    Z.resize(nx*ny) ;
    for ( sizeType i = 0 ; i < nx    ; ++i ) X[i] = i ;
    for ( sizeType i = 0 ; i < ny    ; ++i ) Y[i] = i ;
    for ( sizeType i = 0 ; i < nx*ny ; ++i ) Z[i] = z[i*incz] ;
    makeSpline() ;
  }

  void
  BiCubicSpline::build ( VectorOfValues const & z, sizeType nx, sizeType ny ) {
    SPLINE_ASSERT( z.size() >= nx*ny,
                   "in BilinearSpline::build insufficient dimension of Z, expeced >= " << nx*ny <<
                   " found " << z.size() ) ;
    X.resize(nx) ;
    Y.resize(ny) ;
    Z.resize(nx*ny) ;
    for ( sizeType i = 0 ; i < nx ; ++i ) X[i] = i ;
    for ( sizeType i = 0 ; i < ny ; ++i ) Y[i] = i ;
    std::copy( z.begin(), z.begin()+nx*ny, Z.begin() ) ;
    makeSpline() ;
  }
  
  void
  BiCubicSpline::makeSpline() {
    DX.resize(Z.size()) ;
    DY.resize(Z.size()) ;
    // calcolo derivate
    sizeType nx = sizeType(X.size()) ;
    sizeType ny = sizeType(Y.size()) ;
    PchipSpline sp ;
    for ( sizeType j = 0 ; j < ny ; ++j ) {
      sp.build( &X.front(), 1, &Z[ipos(0,j)], 1, nx ) ;
      for ( sizeType i = 0 ; i < nx ; ++i ) DX[ipos(i,j)] = sp.ypNode(i) ;
    }
    for ( sizeType i = 0 ; i < nx ; ++i ) {
      sp.build( &Y.front(), 1, &Z[ipos(i,0)], nx, ny ) ;
      for ( sizeType j = 0 ; j < ny ; ++j ) DY[ipos(i,j)] = sp.ypNode(j) ;
    }
  }

  void
  BiCubicSpline::load( sizeType i, sizeType j ) const {

    bili[0][0] = Z[ipos(i,j)] ;
    bili[0][1] = Z[ipos(i,j+1)] ;
    bili[1][0] = Z[ipos(i+1,j)];
    bili[1][1] = Z[ipos(i+1,j+1)];

    bili[2][0] = DX[ipos(i,j)] ;
    bili[2][1] = DX[ipos(i,j+1)] ;
    bili[3][0] = DX[ipos(i+1,j)];
    bili[3][1] = DX[ipos(i+1,j+1)];

    bili[0][2] = DY[ipos(i,j)] ;
    bili[0][3] = DY[ipos(i,j+1)] ;
    bili[1][2] = DY[ipos(i+1,j)] ;
    bili[1][3] = DY[ipos(i+1,j+1)];

    bili[2][2] = 0 ;
    bili[2][3] = 0 ;
    bili[3][2] = 0 ;
    bili[3][3] = 0 ;

  }

  valueType
  BiCubicSpline::operator () ( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3( x - X[i], X[i+1] - X[i], u ) ;
    Hermite3( y - Y[j], Y[j+1] - Y[j], v ) ;
    load(i,j) ;
    return bilinear( u, bili, v ) ;
  }

  valueType
  BiCubicSpline::Dx( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3_D( x - X[i], X[i+1] - X[i], u_D ) ;
    Hermite3  ( y - Y[j], Y[j+1] - Y[j], v   ) ;
    load(i,j) ;
    return bilinear( u_D, bili, v ) ;
  }

  valueType
  BiCubicSpline::Dy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3  ( x - X[i], X[i+1] - X[i], u   ) ;
    Hermite3_D( y - Y[j], Y[j+1] - Y[j], v_D ) ;
    load(i,j) ;
    return bilinear( u, bili, v_D ) ;
  }

  valueType
  BiCubicSpline::Dxy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3_D( x - X[i], X[i+1] - X[i], u_D ) ;
    Hermite3_D( y - Y[j], Y[j+1] - Y[j], v_D ) ;
    load(i,j) ;
    return bilinear( u_D, bili, v_D ) ;
  }

  valueType
  BiCubicSpline::Dxx( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3_DD( x - X[i], X[i+1] - X[i], u_DD ) ;
    Hermite3   ( y - Y[j], Y[j+1] - Y[j], v    ) ;
    load(i,j) ;
    return bilinear( u_DD, bili, v ) ;
  }

  valueType
  BiCubicSpline::Dyy( valueType x, valueType y ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3   ( x - X[i], X[i+1] - X[i], u    ) ;
    Hermite3_DD( y - Y[j], Y[j+1] - Y[j], v_DD ) ;
    load(i,j) ;
    return bilinear( u, bili, v_DD ) ;
  }

  void
  BiCubicSpline::D( valueType x, valueType y, valueType d[3] ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3   ( x - X[i], X[i+1] - X[i], u    ) ;
    Hermite3_D ( x - X[i], X[i+1] - X[i], u_D  ) ;
    Hermite3   ( y - Y[j], Y[j+1] - Y[j], v    ) ;
    Hermite3_D ( y - Y[j], Y[j+1] - Y[j], v_D  ) ;
    load(i,j) ;
    d[0] = bilinear( u, bili, v ) ;
    d[1] = bilinear( u_D, bili, v ) ;
    d[2] = bilinear( u, bili, v_D ) ;
  }

  void
  BiCubicSpline::DD( valueType x, valueType y, valueType d[6] ) const {
    sizeType i = search_x( x ) ;
    sizeType j = search_y( y ) ;
    Hermite3   ( x - X[i], X[i+1] - X[i], u    ) ;
    Hermite3_D ( x - X[i], X[i+1] - X[i], u_D  ) ;
    Hermite3_DD( x - X[i], X[i+1] - X[i], u_DD ) ;
    Hermite3   ( y - Y[j], Y[j+1] - Y[j], v    ) ;
    Hermite3_D ( y - Y[j], Y[j+1] - Y[j], v_D  ) ;
    Hermite3_DD( y - Y[j], Y[j+1] - Y[j], v_DD ) ;
    load(i,j) ;
    d[0] = bilinear( u, bili, v ) ;
    d[1] = bilinear( u_D, bili, v ) ;
    d[2] = bilinear( u, bili, v_D ) ;
    d[3] = bilinear( u_DD, bili, v ) ;
    d[4] = bilinear( u_D, bili, v_D ) ;
    d[5] = bilinear( u, bili, v_DD ) ;
  }

  void
  BiCubicSpline::writeToStream( ostream & s ) const {
    s << "Nx = " << X.size() << " Ny = " << Y.size() << '\n' ;
    for ( sizeType i = 1 ; i < sizeType(X.size()) ; ++i ) {
      for ( sizeType j = 1 ; j < sizeType(Y.size()) ; ++j ) {
        s << "patch (" << setw(2) << i << "," << setw(2) << j
          << "): DX = " << X[i]-X[i-1] << " DY = " << Y[j]-Y[j-1]
          << "\n Z00 = " << Z[ipos(i-1,j-1)]
          << " Z01 = " << Z[ipos(i-1,j)]
          << " Z10 = " << Z[ipos(i,j-1)]
          << " Z11 = " << Z[ipos(i,j)]
          << "\n Dx00 = " << DX[ipos(i-1,j-1)]
          << " Dx01 = " << DX[ipos(i-1,j)]
          << " Dx10 = " << DX[ipos(i,j-1)]
          << " ZDx1 = " << DX[ipos(i,j)]
          << "\n Dy00 = " << DY[ipos(i-1,j-1)]
          << " Dy01 = " << DY[ipos(i-1,j)]
          << " Dy10 = " << DY[ipos(i,j-1)]
          << " Dy11 = " << DY[ipos(i,j)]
          << '\n' ;
      }
    }
  }

  char const *
  BiCubicSpline::type_name() const
  { return "BiCubic" ; }

}
