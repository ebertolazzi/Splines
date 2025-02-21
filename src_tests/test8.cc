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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"

#include <GenericContainer/GenericContainer.hh>
#include <fstream>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::real_type;
using Splines::integer;

using namespace GenericContainerNamespace;

template <typename STYPE>
void
testSpline( STYPE & sp ) {
  sp.clear();
  sp.push_back(595,0.644);
  sp.push_back(635,0.652);
  sp.push_back(695,0.644);
  sp.push_back(795,0.694);
  sp.push_back(855,0.907);
  sp.push_back(875,1.336);
  sp.push_back(895,2.169);
  sp.push_back(915,1.598);
  sp.push_back(935,0.916);
  sp.push_back(985,0.607);
  sp.push_back(1035,0.603);
  sp.push_back(1075,0.608);
  sp.build();
  sp.info(cout);

  integer npts = 12;
  static real_type xx[] = { -10,  -9,  -6,  -1,   2,   3,   5,   6,   7,   9,   10,   11 };
  static real_type yy[] = { 595, 635, 695, 795, 855, 875, 895, 915, 935, 985, 1035, 1075 };

  sp.build(xx,yy,npts);
  sp.info(cout);

  GenericContainer gc;
  vec_real_type & x = gc["xdata"].set_vec_real( npts );
  vec_real_type & y = gc["ydata"].set_vec_real( npts );
  std::copy_n( xx, npts, x.begin() );
  std::copy_n( yy, npts, y.begin() );
  sp.build(gc);
  sp.info(cout);
}

int
main() {
  cout << "\n\nTEST N.8\n\n";

  CubicSpline    cs;
  AkimaSpline    ak;
  BesselSpline   bs;
  PchipSpline    pc;
  LinearSpline   ls;
  ConstantSpline csts;
  QuinticSpline  qs;

  testSpline( cs );
  testSpline( ak );
  testSpline( bs );
  testSpline( pc );
  testSpline( ls );
  testSpline( csts );
  testSpline( qs );

  cout << "\nALL DONE!\n\n";
}
