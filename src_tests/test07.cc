/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2026                                                      |
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
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

#include "Splines/Splines.hh"
#include "Utils_fmt.hh"

#include <GenericContainer/GenericContainer.hh>

using namespace SplinesLoad;
using namespace std;
using namespace GenericContainerNamespace;
using Splines::integer;
using Splines::real_type;

// ============================================================================
// Color definitions for fmt output
// ============================================================================

#define HEADER_COLOR fmt::fg( fmt::color::steel_blue ) | fmt::emphasis::bold
#define SUCCESS_COLOR fmt::fg( fmt::color::lime_green ) | fmt::emphasis::bold
#define WARNING_COLOR fmt::fg( fmt::color::gold ) | fmt::emphasis::bold
#define ERROR_COLOR fmt::fg( fmt::color::red ) | fmt::emphasis::bold
#define DATA_COLOR fmt::fg( fmt::color::cornflower_blue )
#define VALUE_COLOR fmt::fg( fmt::color::light_green )
#define SPLINE_COLOR fmt::fg( fmt::color::light_sky_blue ) | fmt::emphasis::bold

// ============================================================================
// Test 1: Spline evaluation with automatic differentiation
// ============================================================================

static real_type test1_xx[] = { 0, 0.9, 2.1, 3, 4.5 };
static real_type test1_yy[] = { 0, 1, 1.1, 2.0, 2.1 };
static integer   test1_npt  = 5;

void test_autodiff_spline( const string & name, Spline * S )
{
  using namespace autodiff::detail;

  S->build( test1_xx, test1_yy, test1_npt );

  autodiff::dual2nd t;
  t.val                 = 1.1;
  t.grad                = 1;
  autodiff::dual2nd ttt = t * t - t / 2;
  autodiff::dual2nd v   = S->eval( ttt );

  // Create a table for results
  fmt::print(
    HEADER_COLOR,
    "\n┌{0:─^{2}}┐\n"
    "│{1: ^{2}}│\n"
    "└{0:─^{2}}┘\n",
    "",
    name,
    60 );

  fmt::print( "  {} = {}\n", fmt::styled( "t", DATA_COLOR ), fmt::styled( val( ttt ), VALUE_COLOR ) );
  fmt::print( "  {} = {}\n", fmt::styled( "t'", DATA_COLOR ), fmt::styled( val( ttt.grad ), VALUE_COLOR ) );
  fmt::print( "  {} = {}\n", fmt::styled( "S(t)", DATA_COLOR ), fmt::styled( val( v ), VALUE_COLOR ) );
  fmt::print(
    "  {} = {} (analytic: {})\n",
    fmt::styled( "S'(t)", DATA_COLOR ),
    fmt::styled( val( v.grad ), VALUE_COLOR ),
    fmt::styled( val( ttt.grad ) * S->D( val( ttt ) ), VALUE_COLOR ) );
  fmt::print(
    "  {} = {} (analytic: {})\n",
    fmt::styled( "S''(t)", DATA_COLOR ),
    fmt::styled( val( v.grad.grad ), VALUE_COLOR ),
    fmt::styled(
      val( ttt.grad.grad ) * S->D( val( ttt ) ) + val( ttt.grad * ttt.grad ) * S->DD( val( ttt ) ),
      VALUE_COLOR ) );
}

// ============================================================================
// Test 2: Various spline constructions
// ============================================================================

template <typename STYPE> void test_spline_construction( const string & name, STYPE & sp )
{
  fmt::print(
    HEADER_COLOR,
    "\n┌{0:─^{2}}┐\n"
    "│{1: ^{2}}│\n"
    "└{0:─^{2}}┘\n",
    "",
    name + " Construction Test",
    70 );

  // First dataset
  fmt::print( DATA_COLOR, "\n  Dataset 1 (push_back method):\n" );
  sp.clear();

  vector<pair<real_type, real_type>> dataset1 = { { 595, 0.644 }, { 635, 0.652 }, { 695, 0.644 },  { 795, 0.694 },
                                                  { 855, 0.907 }, { 875, 1.336 }, { 895, 2.169 },  { 915, 1.598 },
                                                  { 935, 0.916 }, { 985, 0.607 }, { 1035, 0.603 }, { 1075, 0.608 } };

  // Print dataset table
  fmt::print( "  ┌────────────┬────────────┐\n" );
  fmt::print( "  │     x      │     y      │\n" );
  fmt::print( "  ├────────────┼────────────┤\n" );
  for ( const auto & p : dataset1 ) { fmt::print( "  │ {:10.3f} │ {:10.3f} │\n", p.first, p.second ); }
  fmt::print( "  └────────────┴────────────┘\n" );

  for ( const auto & p : dataset1 ) { sp.push_back( p.first, p.second ); }
  sp.build();

  // Second dataset
  fmt::print( DATA_COLOR, "\n  Dataset 2 (array build method):\n" );
  integer   npts = 12;
  real_type xx[] = { -10, -9, -6, -1, 2, 3, 5, 6, 7, 9, 10, 11 };
  real_type yy[] = { 595, 635, 695, 795, 855, 875, 895, 915, 935, 985, 1035, 1075 };

  fmt::print( "  ┌────────────┬────────────┐\n" );
  fmt::print( "  │     x      │     y      │\n" );
  fmt::print( "  ├────────────┼────────────┤\n" );
  for ( integer i = 0; i < npts; ++i ) { fmt::print( "  │ {:10.3f} │ {:10.3f} │\n", xx[i], yy[i] ); }
  fmt::print( "  └────────────┴────────────┘\n" );

  sp.build( xx, yy, npts );

  // Third dataset using GenericContainer
  fmt::print( DATA_COLOR, "\n  Dataset 3 (GenericContainer method):\n" );
  GenericContainer gc;
  vec_real_type &  x = gc["xdata"].set_vec_real( npts );
  vec_real_type &  y = gc["ydata"].set_vec_real( npts );
  copy_n( xx, npts, x.begin() );
  copy_n( yy, npts, y.begin() );
  sp.build( gc );

  // Display spline information
  fmt::print( SPLINE_COLOR, "\n  Spline Information:\n" );
  stringstream info_stream;
  sp.info( info_stream );
  fmt::print( "{}\n", info_stream.str() );
}

// ============================================================================
// Test 3: SplineSet construction
// ============================================================================

void test_spline_set()
{
  fmt::print(
    HEADER_COLOR,
    "\n┌{0:─^{2}}┐\n"
    "│{1: ^{2}}│\n"
    "└{0:─^{2}}┘\n",
    "",
    "SplineSet Construction Test",
    70 );

  constexpr unsigned npts = 10;
  constexpr unsigned nspl = 4;

  GenericContainer  gc;
  vec_real_type &   X           = gc["xdata"].set_vec_real();
  vec_string_type & spline_type = gc["spline_type"].set_vec_string();
  vec_string_type & headers     = gc["headers"].set_vec_string();
  mat_real_type &   Y           = gc["ydata"].set_mat_real( npts, nspl );
  map_type &        Yp          = gc["ypdata"].set_map();

  // Headers
  headers.emplace_back( "sp1" );
  headers.emplace_back( "sp2" );
  headers.emplace_back( "sp3" );
  headers.emplace_back( "sp4" );

  // Spline types
  spline_type.emplace_back( "cubic" );
  spline_type.emplace_back( "akima" );
  spline_type.emplace_back( "vanleer" );
  spline_type.emplace_back( "hermite" );

  // X values
  real_type x_vals[] = { 0, 1, 2, 3, 4, 4.1, 4.2, 5, 6, 7 };
  for ( const auto & val : x_vals ) { X.emplace_back( val ); }

  // Y values (random data for demonstration)
  fmt::print(
    DATA_COLOR,
    "\n"
    "  Dataset Configuration:\n"
    "  ┌────────────┬────────────┬────────────┬────────────┬────────────┐\n"
    "  │     x      │   sp1 (y)  │   sp2 (y)  │   sp3 (y)  │   sp4 (y)  │\n"
    "  ├────────────┼────────────┼────────────┼────────────┼────────────┤\n" );

  for ( unsigned i = 0; i < npts; ++i )
  {
    Y( i, 0 ) = sin( x_vals[i] );
    Y( i, 1 ) = cos( x_vals[i] );
    Y( i, 2 ) = exp( -0.1 * x_vals[i] );
    Y( i, 3 ) = log( 1 + x_vals[i] );

    fmt::print(
      "  │ {:10.3f} │ {:10.3f} │ {:10.3f} │ {:10.3f} │ {:10.3f} │\n",
      x_vals[i],
      Y( i, 0 ),
      Y( i, 1 ),
      Y( i, 2 ),
      Y( i, 3 ) );
  }
  fmt::print( "  └────────────┴────────────┴────────────┴────────────┴────────────┘\n" );

  // Derivative data for sp4 (hermite)
  vec_real_type & yp = Yp["sp4"].set_vec_real();
  for ( const auto & val : x_vals )
  {
    yp.emplace_back( 1.0 / ( 1.0 + val ) );  // derivative of log(1+x)
  }

  // Build SplineSet
  fmt::print( SPLINE_COLOR, "\n  Building SplineSet...\n" );
  SplineSet ss;
  ss.build( gc );

  fmt::print( SUCCESS_COLOR, "  ✓ SplineSet successfully built with {} splines\n", nspl );
}

// ============================================================================
// Main test runner
// ============================================================================

int main()
{
  // Banner
  fmt::print(
    fmt::emphasis::bold | fmt::fg( fmt::color::royal_blue ),
    "\n"
    "╔══════════════════════════════════════════════════════════════════╗\n"
    "║                 SPLINE LIBRARY TEST SUITE                        ║\n"
    "║                    Combined Version                              ║\n"
    "╚══════════════════════════════════════════════════════════════════╝\n\n" );

  // Test 2: Spline construction
  fmt::print(
    HEADER_COLOR,
    "\n"
    "┌────────────────────────────────────────────────────────────┐\n"
    "│                     TEST 1: Construction                   │\n"
    "└────────────────────────────────────────────────────────────┘\n" );

  CubicSpline    cs;
  AkimaSpline    ak;
  VanLeerSpline  vl;
  PchipSpline    pc;
  LinearSpline   ls;
  ConstantSpline csts;
  QuinticSpline  qs_cubic( Splines::Spline_sub_type::CUBIC );
  QuinticSpline  qs_akima( Splines::Spline_sub_type::AKIMA );
  QuinticSpline  qs_vanleer( Splines::Spline_sub_type::VANLEER );
  QuinticSpline  qs_pchip( Splines::Spline_sub_type::PCHIP );

  test_spline_construction( "Constant", csts );
  test_spline_construction( "Linear", ls );
  test_spline_construction( "Cubic", cs );
  test_spline_construction( "Akima", ak );
  test_spline_construction( "VanLeer", vl );
  test_spline_construction( "Pchip", pc );
  test_spline_construction( "Quintic", qs_cubic );
  test_spline_construction( "Quintic_akima", qs_akima );
  test_spline_construction( "Quintic_vanleer", qs_vanleer );
  test_spline_construction( "Quintic_pchip", qs_pchip );

  // Test 2: AutoDiff spline evaluation
  fmt::print(
    HEADER_COLOR,
    "┌────────────────────────────────────────────────────────────┐\n"
    "│                     TEST 2: AutoDiff                       │\n"
    "└────────────────────────────────────────────────────────────┘\n" );

  test_autodiff_spline( "Constant", &csts );
  test_autodiff_spline( "Linear", &ls );
  test_autodiff_spline( "Cubic", &cs );
  test_autodiff_spline( "Akima", &ak );
  test_autodiff_spline( "VanLeer", &vl );
  test_autodiff_spline( "Pchip", &pc );
  test_autodiff_spline( "Quintic", &qs_cubic );
  test_autodiff_spline( "Quintic_akima", &qs_akima );
  test_autodiff_spline( "Quintic_vanleer", &qs_vanleer );
  test_autodiff_spline( "Quintic_pchip", &qs_pchip );

  // Test 3: SplineSet
  fmt::print(
    HEADER_COLOR,
    "\n"
    "┌────────────────────────────────────────────────────────────┐\n"
    "│                     TEST 3: SplineSet                      │\n"
    "└────────────────────────────────────────────────────────────┘\n" );

  test_spline_set();

  // Completion message
  fmt::print(
    SUCCESS_COLOR,
    "\n┌{0:─^{2}}┐\n"
    "│{1: ^{2}}│\n"
    "└{0:─^{2}}┘\n",
    "",
    "✓ ALL TESTS COMPLETED SUCCESSFULLY ✓",
    60 );

  return 0;
}
