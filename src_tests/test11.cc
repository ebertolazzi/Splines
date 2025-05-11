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
 |      Università degli Studi di Trento                                    |
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
using namespace GenericContainerNamespace;
using Splines::real_type;
using Splines::integer;

int
main() {
  cout << "\n\nTEST N.11\n\n";

  Splines::SplineSurf *_p_spline1 = new BilinearSpline();
  Splines::SplineSurf *_p_spline2 = new BiCubicSpline();
  Splines::SplineSurf *_p_spline3 = new BiQuinticSpline();
  Splines::SplineSurf *_p_spline4 = new Akima2Dspline();

  _p_spline1->build( string("AeroMap_SM.json") );
  _p_spline1->write_to_stream(cout);

  _p_spline2->build( string("AeroMap_SM.json") );
  _p_spline2->write_to_stream(cout);

  _p_spline3->build( string("AeroMap_SM.json") );
  _p_spline3->write_to_stream(cout);

  _p_spline4->build( string("AeroMap_SM.json") );
  _p_spline4->write_to_stream(cout);

  cout << "\nALL DONE!\n\n";
}
