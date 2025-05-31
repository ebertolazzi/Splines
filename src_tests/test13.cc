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

constexpr int nx = std::size(x); // 10
constexpr int ny = std::size(y); // 9

// data is stored in C-like ordering: the element at graphical position (i,j) (which is at the
// i * ny + j  1-dimensional position) refers to the function evaluated at j-th x and i-th y

constexpr double z[] = {24.2, 24.0, 20.3, 17.3, 14.5, 12.2, 10.2,  5.7,  3.4, 0.1,
                        28.0, 24.6, 21.1, 18.1, 15.2, 12.8, 10.7,  6.5,  3.9, 0.2,
                        28.3, 25.2, 21.9, 18.7, 15.9, 13.4, 11.2,  7.3,  4.4, 0.4,
                        30.8, 27.2, 23.8, 20.5, 17.3, 14.7, 12.3,  8.1,  4.9, 0.8,
                        34.5, 30.3, 26.6, 23.2, 19.8, 16.8, 14.1,  9.4,  5.6, 1.1,
                        37.9, 34.3, 30.4, 26.8, 23.3, 19.8, 16.8, 11.2,  6.8, 1.4,
                        36.1, 38.0, 34.9, 31.3, 27.3, 23.6, 20.1, 13.4,  8.3, 1.7,
                        36.1, 36.6, 38.5, 36.1, 31.6, 28.1, 24.2, 16.2, 10.0, 2.2,
                        36.1, 35.2, 42.1, 38.7, 35.7, 32.0, 28.1, 19.3, 11.9, 2.9};

template <typename Tsurf>
static
void
do_test( string const & msg ) {
  using namespace autodiff::detail;

  Tsurf S;

  constexpr int ldZ{ nx };
  S.build(x, 1, y, 1, z, ldZ, nx, ny, false,  false);
  //S.write_to_stream(cout);
  
  double x0{ 40 };
  double y0{ 0.7 };
  
  autodiff::dual2nd X0{x0};
  autodiff::dual2nd Y0{y0};
  
  X0.grad = 0.5;
  Y0.grad = 0.5;
  autodiff::dual2nd V1{ S.eval( X0, Y0 ) };
  X0.grad = 0.5;
  Y0.grad = -0.5;
  autodiff::dual2nd V2{ S.eval( X0, Y0 ) };
  
  autodiff::dual2nd V3{ V1-V2 };

  fmt::print( "{}\n", msg );

  fmt::print( "S    ({},{}) = {:12.6}\n", x0, y0, S.eval( x0, y0 ) );
  fmt::print( "S_x  ({},{}) = {:12.6}\n", x0, y0, S.Dx( x0, y0 ) );
  fmt::print( "S_y  ({},{}) = {:12.6}\n", x0, y0, S.Dy( x0, y0 ) );
  fmt::print( "S_xx ({},{}) = {:12.6}\n", x0, y0, S.Dxx( x0, y0 ) );
  fmt::print( "S_xy ({},{}) = {:12.6}\n", x0, y0, S.Dxy( x0, y0 ) );
  fmt::print( "S_yy ({},{}) = {:12.6}\n", x0, y0, S.Dyy( x0, y0 ) );
  fmt::print( "GRAD ({},{}) = {:12.6}\n", x0, y0, val(V3.grad) );
  fmt::print( "GGRAD({},{}) = {:12.6}\n", x0, y0, V3.grad.grad );
  std::cout << "\n";
}

int
main() {
  cout << "\n\nTEST N.12\n\n";

  do_test<BiQuinticSpline>("BiQuinticSpline");
  do_test<BiCubicSpline>("BiCubicSpline");
  do_test<BilinearSpline>("BilinearSpline");
  do_test<Akima2Dspline>("Akima2Dspline");

  cout << "\nALL DONE!\n\n";
}
