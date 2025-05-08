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

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif
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

  cout << "\n\nTEST N.9\n\n";

  constexpr unsigned npts = 10;
  constexpr unsigned nspl = 4;

  GenericContainer gc;

  vec_real_type   & X           = gc["xdata"].set_vec_real();
  vec_string_type & spline_type = gc["spline_type"].set_vec_string();
  vec_string_type & headers     = gc["headers"].set_vec_string();
  mat_real_type   & Y           = gc["ydata"].set_mat_real(npts,nspl);
  map_type        & Yp          = gc["ypdata"].set_map();

  headers.emplace_back("sp1");
  headers.emplace_back("sp2");
  headers.emplace_back("sp3");
  headers.emplace_back("sp4");

  spline_type.emplace_back("cubic");
  spline_type.emplace_back("akima");
  spline_type.emplace_back("bessel");
  spline_type.emplace_back("hermite");

  X.emplace_back(0);
  X.emplace_back(1);
  X.emplace_back(2);
  X.emplace_back(3);
  X.emplace_back(4);
  X.emplace_back(4.1);
  X.emplace_back(4.2);
  X.emplace_back(5);
  X.emplace_back(6);
  X.emplace_back(7);

  vec_real_type & yp = Yp["sp4"].set_vec_real();

  yp.emplace_back(0);
  yp.emplace_back(1);
  yp.emplace_back(2);
  yp.emplace_back(3);
  yp.emplace_back(4);
  yp.emplace_back(4.1);
  yp.emplace_back(4.2);
  yp.emplace_back(5);
  yp.emplace_back(6);
  yp.emplace_back(7);

  SplineSet ss;

  // check constructor
  cout << "build SplineSet\n";
  ss.build(gc);

  cout << "\nALL DONE!\n\n";
}
