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


// monotone
static real_type xx[] = { 0, 0.9, 2.1, 3, 4.5 };
static real_type yy[] = { 0, 1, 1.1, 2.0, 2.1 };

static integer npt = 5;

template <typename Tspline>
static
void
do_test( string const & msg ) {
  using namespace autodiff::detail;

  Tspline S;
  S.build( xx, yy, npt );
  
  autodiff::dual2nd t{ 1.1 };
  t.grad = 1;
  autodiff::dual2nd ttt{ t*t-t/2 };
  autodiff::dual2nd v{ S( ttt ) };
  
  fmt::print( "t = {}  t' = {}\n", ttt, ttt.grad );
  fmt::print( "{}  S({})   = {}\n",    msg, ttt, v );
  fmt::print( "{}  S'({})  = {}:{}\n", msg, ttt, val(ttt.grad)*S.D( val(ttt) ), v.grad );
  fmt::print( "{}  S''({}) = {}:{}\n\n", msg, ttt, val(ttt.grad.grad)*S.D( val(ttt) )
                                                  +val(ttt.grad*ttt.grad)*S.DD( val(ttt) ),  v.grad.grad );
  //fmt::print( "S''({}) = {}:{}\n", t, S.DD( tt ), v.grad.grad );

}


int
main() {
  cout << "\n\nTEST N.12\n\n";

  do_test<LinearSpline>("LinearSpline");
  do_test<ConstantSpline>("ConstantSpline");
  do_test<AkimaSpline>("AkimaSpline");
  do_test<CubicSpline>("CubicSpline");
  do_test<BesselSpline>("BesselSpline");
  do_test<PchipSpline>("PchipSpline");
  do_test<QuinticSpline>("QuinticSpline");

  cout << "\nALL DONE!\n\n";
}
