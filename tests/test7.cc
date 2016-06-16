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

#include "GenericContainer.hh"
#include "Splines.hh"
#include <fstream>

using namespace SplinesLoad ;
using namespace std ;
using Splines::valueType ;
using Splines::indexType ;
using Splines::sizeType ;

// monotone
//valueType xx[] = { 595, 635, 695, 795, 855, 875, 895, 915, 935, 985, 1035, 1075 } ;
//valueType yy[] = { 0.644, 0.652, 0.644, 0.694, 0.907, 1.336, 2.169, 1.598, 0.916, 0.607, 0.603, 0.608 } ;
//valueType xx[] = { 0, 1, 2, 3, 5, 6, 7, 9, 10, 11, 12 } ;
//valueType yy[] = { 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4 } ;
valueType xx[] = { -10, -9, -6, -1, 2, 3, 5, 6, 7, 9, 10, 11, 12 } ;
valueType yy[] = { 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4 } ;
//valueType yy[] = {  595, 635, 695, 795, 855, 875, 895, 915, 935, 985, 1035, 1075 } ;
//valueType xx[] = { 0, 1, 2, 3, 4, 8, 10 } ;
//valueType yy[] = { 0, 1, 2, 3, 4 } ;
//valueType yy[] = { 1, 1, 1, 1 } ;
//valueType yy[] = { 0, 1, 4, 9, 16, 64, 100 } ;
sizeType  n = sizeof(xx)/sizeof(xx[0]) ;

#define SAVE(S) \
cout << #S": n = " << n << '\n' ; \
S.clear() ; \
S.reserve(n) ; \
for ( indexType i = 0 ; i < n ; ++i ) S.pushBack(xx[i],yy[i]) ; \
S.build() ; \
cout << #S": xMin    = " << S.xMin() << '\n' ; \
cout << #S": xMax    = " << S.xMax() << '\n' ; \
cout << #S": xx[0]   = " << xx[0]    << '\n' ; \
cout << #S": xx[end] = " << xx[n-1]  << '\n' ; \
file_##S << "x\ty\tDy\tDDy\n" ; \
for ( valueType x = xmin-(xmax-xmin)*0.01 ; x <= xmax+(xmax-xmin)*0.01 ; x += (xmax-xmin)/1000 ) \
  file_##S << x << '\t' << S(x) << '\t' << S.D(x) << '\t' << S.DD(x) << '\n' ; \
file_##S.close()

int
main() {

  valueType xmin = xx[0] ;
  valueType xmax = xx[n-1] ;

  PchipSpline pc ;
  BSpline<3>  bs ;
  ofstream    file_pc("out/PchipSpline.txt") ;
  ofstream    file_bs("out/BSpline.txt") ;
  ofstream    file_bs1("out/BSpline1.txt") ;
  ofstream    file_bs2("out/_DATA.txt") ;

  SAVE(pc) ;
  SAVE(bs) ;

  file_bs1 << "x" ;
  for ( sizeType i = 0 ; i < n ; ++i ) file_bs1 << "\tb" << i ;
  file_bs1 << '\n' ;

  for ( valueType x = xmin-(xmax-xmin)*0.01 ; x <= xmax+(xmax-xmin)*0.01 ; x += (xmax-xmin)/1000 ) {
    file_bs1 << x ;
    valueType vals[n] ;
    bs.bases( x, vals) ;
    for ( sizeType i = 0 ; i < n ; ++i ) file_bs1 << '\t' << vals[i] ;
    file_bs1 << '\n' ;
  }
  file_bs1.close() ;

  file_bs2 << "x\ty\tz\n" ;
  for ( sizeType i = 0 ; i < n ; ++i ) file_bs2 << xx[i] << '\t' << yy[i] << "\t0\n" ;
  file_bs2.close() ;

  cout << "v = " << bs(xmin-0.1) << '\n' ;
  cout << "v = " << bs(xmin+0.001) << '\n' ;
  cout << "v = " << bs(xmax-0.1) << '\n' ;
  cout << "v = " << bs(xmax) << '\n' ;
  cout << "v = " << bs(xmax+0.1) << '\n' ;

  cout << "ALL DONE!\n" ;
}
