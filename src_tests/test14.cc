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

constexpr double x[]{ 0,   5,  10,  15,  20, 25,  30,  40,  50, 70 };
constexpr double y[]{ 0, 0.2, 0.4, 0.6, 0.8,  1, 1.2, 1.4, 1.6 };

constexpr int nx = std::size(x);
constexpr int ny = std::size(y);

int
main() {
  cout << "\n\nTEST N.14\n\n";
  
  GC::GenericContainer gc;

  GC::GenericContainer & S0{ gc["spline0"] };
  GC::GenericContainer & S1{ gc["spline1"] };
  
  GC::vec_real_type & V0x = S0["xdata"].set_vec_real(nx); std::copy_n( x, nx, V0x.data() );
  GC::vec_real_type & V0y = S0["ydata"].set_vec_real(nx); std::copy_n( y, nx, V0y.data() );
  GC::vec_real_type & V1x = S1["xdata"].set_vec_real(nx); std::copy_n( x, nx, V1x.data() );
  GC::vec_real_type & V1y = S1["ydata"].set_vec_real(nx); std::copy_n( y, nx, V1y.data() );
  
  S0["spline_type"]     = "cubic";
  S0["bc_begin"]        = "natural";
  S0["bc_end"]          = "natural";
  S1["spline_type"]     = "quintic";
  S1["spline_sub_type"] = "pchip";
  S1["bc_begin"]        = "natural";
  S1["bc_end"]          = "natural";
  
  Splines::Spline1Dblend S("pippo");
  
  S.build( gc );

  cout << "\nALL DONE!\n\n";
}
