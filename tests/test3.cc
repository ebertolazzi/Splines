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
#include <fstream>

using namespace SplinesLoad ;
using namespace std ;
using Splines::valueType ;
using Splines::indexType ;
using Splines::sizeType ;

// Test problem for Akima interpolation
// Ref. : Hiroshi Akima, Journal of the ACM, Vol. 17, No. 4, October 1970, pages 589-602.

valueType xx0[] = {  0,  1,  2,  3,  4,  5,    6,  7,  8,  9, 10 } ; 
valueType yy0[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 } ;

valueType xx1[] = {  0,  1,  3,  4,  6,  7,    9, 10, 12, 13, 15 } ;
valueType yy1[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 } ;

valueType xx2[] = {  0,  2,  3,  5,  6,  8,    9, 11, 12, 14, 15 } ;
valueType yy2[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 } ;

// RPN 14
valueType xx3[] = { 7.99, 8.09,       8.19,       8.7,      9.2,      10,       12,       15,       20       } ;
valueType yy3[] = { 0,    2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994 } ;

// Titanium
valueType xx4[] = { 595,   635,   695,   795,   855,   875,   895,   915,   935,   985,   1035,  1075  } ;
valueType yy4[] = { 0.644, 0.652, 0.644, 0.694, 0.907, 1.336, 2.169, 1.598, 0.916, 0.607, 0.603, 0.608 } ;

// toolpath
valueType xx5[] = { 0.11, 0.12, 0.15, 0.16 } ;
valueType yy5[] = { 0.0003, 0.0003, 0.0004, 0.0004 } ;

// monotone
valueType xx6[] = { 0, 1, 2, 3, 4 } ;
valueType yy6[] = { 0, 1, 1.1, 2.0, 2.1 } ;

sizeType  n[]   = { 11, 11, 11, 9, 12, 4, 5 } ;

int
main() {

  SplineSet   ss ;
  ofstream    file, file_D ;

  for ( indexType k = 0 ; k < 7 ; ++ k ) {
    valueType * xx, * yy ;
    switch ( k ) {
      case 0: xx = xx0 ; yy = yy0 ; break ;
      case 1: xx = xx1 ; yy = yy1 ; break ;
      case 2: xx = xx2 ; yy = yy2 ; break ;
      case 3: xx = xx3 ; yy = yy3 ; break ;
      case 4: xx = xx4 ; yy = yy4 ; break ;
      case 5: xx = xx5 ; yy = yy5 ; break ;
      case 6: xx = xx6 ; yy = yy6 ; break ;
    }
    char fname[100] ;
    sprintf( fname, "out/SplineSet%d.txt", k) ;
    file.open(fname) ;

    sprintf( fname, "out/SplineSet%d_D.txt", k) ;
    file_D.open(fname) ;
    
    valueType xmin = xx[0] ;
    valueType xmax = xx[n[k]-1] ;

    indexType nspl = 7 ;
    indexType npts = n[k] ;

    char const *headers[] = {
      "SPLINE_CONSTANT",
      "SPLINE_LINEAR",
      "SPLINE_CUBIC_BASE",
      "SPLINE_CUBIC",
      "SPLINE_AKIMA",
      "SPLINE_BESSEL",
      "SPLINE_PCHIP",
      "SPLINE_QUINTIC"
    } ;
    
    SplineType const stype[] = {
       Splines::CONSTANT_TYPE,
       Splines::LINEAR_TYPE,
       Splines::CUBIC_BASE_TYPE,
       Splines::CUBIC_TYPE,
       Splines::AKIMA_TYPE,
       Splines::BESSEL_TYPE,
       Splines::PCHIP_TYPE,
       Splines::QUINTIC_TYPE
    } ;

    //bool const rp_policy[] = { true, true, true, true, true, true, true } ;

    std::vector<valueType> YpZero(npts) ;
    std::fill(YpZero.begin(), YpZero.end(), 0 ) ;

    valueType const *Y[]  = { yy, yy, yy, yy, yy, yy, yy, yy } ;
    valueType const *Yp[] = { nullptr, nullptr, &YpZero.front(), nullptr, nullptr, nullptr, nullptr, nullptr } ;

    ss.build( nspl, npts, headers, stype, xx, Y, Yp ) ;
    ss.info(cout) ;

    file   << "x" ;
    file_D << "x" ;
    for ( indexType i = 0 ; i < nspl ; ++i ) {
      file   << '\t' << ss.header(i) ;
      file_D << '\t' << ss.header(i) ;
    }
    file   << '\n' ;
    file_D << '\n' ;
    for ( valueType x = xmin-(xmax-xmin)*0.01 ; x <= xmax+(xmax-xmin)*0.01 ; x += (xmax-xmin)/1000 ) {
      file   << x ;
      file_D << x ;
      for ( indexType i = 0 ; i < nspl ; ++i ) {
        file   << '\t' << ss(x,i) ;
        file_D << '\t' << ss.D(x,i) ;
      }
      file   << '\n' ;
      file_D << '\n' ;
    }
    file.close() ;
    file_D.close() ;
  }
  
  cout << "ALL DONE!\n" ;
}
