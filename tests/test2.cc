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

valueType x[] = {  0,  1,  2,  3 } ;
valueType y[] = {  0,  1,  2 } ;
valueType z[] = {  10, 10, 10,   10,
                   10, 10, 10.5, 11,
                   11, 12, 13,   8 } ;

int
main() {

  BiCubicSpline  bc ;
  BilinearSpline bl ;
  ofstream       file_bl("bilinear.txt") ;
  ofstream       file_bc("bicubic.txt") ;
  
  bc.build( x, 1, y, 1, z, 1, 4, 3 ) ;
  bl.build( x, 1, y, 1, z, 1, 4, 3 ) ;
  
  bl.writeToStream( cout ) ;
  
  for ( int i = 0 ; i <= 100 ; ++i ) {
    valueType x = bc.xMin() + (bc.xMax()-bc.xMin())*i/100.0 ;
    for ( int j = 0 ; j <= 100 ; ++j ) {
      valueType y = bc.yMin() + (bc.yMax()-bc.yMin())*j/100.0 ;
      file_bc << bc(x,y) << '\t' ;
      file_bl << bl(x,y) << '\t' ;
    }
    file_bc << '\n' ;
    file_bl << '\n' ;
  }
  
  file_bc.close() ;
  file_bl.close() ;
  
  cout << "ALL DONE!\n" ;
}
