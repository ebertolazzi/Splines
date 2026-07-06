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
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Splines/Splines.hh"
#include "Utils_fmt.hh"

#include <GenericContainer/GenericContainerInterface_nlohmann.hh>

#include <cmath>
#include <string>
#include <string_view>

using namespace SplinesLoad;
using namespace std;
using Splines::integer;
using Splines::real_type;

using json = nlohmann::json;

static int n_fail = 0;

static void check( bool ok, string_view msg )
{
  if ( ok ) { fmt::print( fg( fmt::color::green ), "  ✓ {}\n", msg ); }
  else
  {
    ++n_fail;
    fmt::print( fg( fmt::color::red ), "  ✗ {}\n", msg );
  }
}

static void check_close( real_type a, real_type b, real_type tol, string_view msg )
{
  real_type const err = abs( a - b );
  if ( err <= tol ) { fmt::print( fg( fmt::color::green ), "  ✓ {} (err = {:.3e})\n", msg, err ); }
  else
  {
    ++n_fail;
    fmt::print( fg( fmt::color::red ), "  ✗ {} (err = {:.3e}, tol = {:.3e}, a = {}, b = {})\n", msg, err, tol, a, b );
  }
}

static void section( string_view title ) { fmt::print( fg( fmt::color::yellow ) | fmt::emphasis::bold, "\n{}\n", title ); }

// ---------------------------------------------------------------------------
// 1. json -> GenericContainer -> json round trip
// ---------------------------------------------------------------------------
static void test_json_gc_json_roundtrip()
{
  section( "1. json -> GenericContainer -> json round trip" );

  json const j1 = {
    { "an_int", 42 },
    { "a_real", 3.14159 },
    { "a_bool", true },
    { "a_null", nullptr },
    { "a_string", "hello" },
    { "vec_int", { 1, 2, 3 } },
    { "vec_real", { 0.5, 1.5, 2.5 } },
    { "vec_mixed", { 1, 2.5, 3 } },  // must collapse to vec_real
    { "vec_string", { "a", "b", "c" } },
    { "nested", { { "inner", { 1.0, 2.0 } }, { "name", "nested_map" } } },
  };

  auto const gc = j1.get<GC::GenericContainer>();
  json const j2 = gc;

  // nlohmann compares integer vs float numerics by value, so the collapse
  // of "vec_mixed" to vec_real does not break equality
  check( j1 == j2, "json -> GC -> json preserves the whole document" );
  check( gc( "an_int" ).get_type() == GC::GC_type::INTEGER, "integer arrives as INTEGER" );
  check( gc( "a_real" ).get_type() == GC::GC_type::REAL, "real arrives as REAL" );
  check( gc( "vec_real" ).get_type() == GC::GC_type::VEC_REAL, "[0.5,1.5,2.5] collapses to VEC_REAL" );
  check( gc( "vec_mixed" ).get_type() == GC::GC_type::VEC_REAL, "[1,2.5,3] collapses to VEC_REAL" );
  check( gc( "a_null" ).get_type() == GC::GC_type::NOTYPE, "null arrives as NOTYPE" );
}

// ---------------------------------------------------------------------------
// 2. GenericContainer -> json -> GenericContainer for GC-specific types
// ---------------------------------------------------------------------------
static void test_gc_json_gc_roundtrip()
{
  section( "2. GenericContainer -> json -> GenericContainer (vec_real, mat_real)" );

  GC::GenericContainer gc;
  GC::vec_real_type &  v = gc["v"].set_vec_real();
  v = { 0.0, 0.9, 2.1, 3.0, 4.5 };

  GC::mat_real_type & M = gc["M"].set_mat_real( 2, 3 );
  for ( unsigned i = 0; i < 2; ++i )
    for ( unsigned jc = 0; jc < 3; ++jc ) M( i, jc ) = 10.0 * i + jc;

  json const j   = gc;
  auto const gc2 = j.get<GC::GenericContainer>();

  check( gc2( "v" ).get_type() == GC::GC_type::VEC_REAL, "vec_real survives round trip" );
  bool v_ok = true;
  {
    GC::vec_real_type const & v2 = gc2( "v" ).get_vec_real();
    v_ok                         = v2.size() == v.size();
    for ( size_t i = 0; v_ok && i < v.size(); ++i ) v_ok = v2[i] == v[i];
  }
  check( v_ok, "vec_real values are identical" );

  check( gc2( "M" ).get_type() == GC::GC_type::MAT_REAL, "mat_real survives round trip (via collapse)" );
  if ( gc2( "M" ).get_type() == GC::GC_type::MAT_REAL )
  {
    GC::mat_real_type const & M2   = gc2( "M" ).get_mat_real();
    bool                      m_ok = M2.num_rows() == 2 && M2.num_cols() == 3;
    for ( unsigned i = 0; m_ok && i < 2; ++i )
      for ( unsigned jc = 0; m_ok && jc < 3; ++jc ) m_ok = M2( i, jc ) == M( i, jc );
    check( m_ok, "mat_real dimensions and values are identical" );
  }
}

// ---------------------------------------------------------------------------
// 3. Build a CubicSpline from a nlohmann::json document
// ---------------------------------------------------------------------------
static void test_cubic_spline_from_json()
{
  section( "3. CubicSpline built from nlohmann::json" );

  json const j = {
    { "spline_type", "cubic" },
    { "xdata", { 0, 0.9, 2.1, 3, 4.5 } },  // mixed int/real on purpose
    { "ydata", { 0, 1, 1.99, 2.0, 2.1 } },
    { "bc_begin", "natural" },
    { "bc_end", "natural" },
  };

  Splines::CubicSpline spline;
  spline.build( j.get<GC::GenericContainer>() );

  real_type const xx[] = { 0, 0.9, 2.1, 3, 4.5 };
  real_type const yy[] = { 0, 1, 1.99, 2.0, 2.1 };
  for ( integer i = 0; i < 5; ++i )
  {
    check_close( spline( xx[i] ), yy[i], 1e-12, fmt::format( "interpolates knot #{} ({}, {})", i, xx[i], yy[i] ) );
  }

  // natural bc: second derivative vanishes at both ends
  check_close( spline.DD( xx[0] ), 0.0, 1e-9, "natural bc: y''(x_begin) = 0" );
  check_close( spline.DD( xx[4] ), 0.0, 1e-9, "natural bc: y''(x_end) = 0" );
}

// ---------------------------------------------------------------------------
// 4. SplineSet from json, cross-checked through the nlohmann bridge
// ---------------------------------------------------------------------------
static void test_spline_set_from_json()
{
  section( "4. SplineSet from json through nlohmann bridge" );

  string const json_text = R"({
    "spline_type": [ "cubic", "akima", "pchip" ],
    "headers":     [ "cubic", "akima", "pchip" ],
    "xdata": [ 0, 0.9, 2.1, 3, 4.5 ],
    "ydata": [ [ 0, 1, 1.99, 2.0, 2.1 ],
               [ 0, 1, 1.99, 2.0, 2.1 ],
               [ 0, 1, 1.99, 2.0, 2.1 ] ]
  })";

  auto const gc_a = json::parse( json_text ).get<GC::GenericContainer>();
  auto const gc_b = json::parse( json_text ).get<GC::GenericContainer>();

  SplineSet ss_a, ss_b;
  ss_a.build( gc_a );
  ss_b.build( gc_b );

  check( ss_a.num_splines() == 3, "SplineSet from json has 3 splines" );

  bool same = true;
  for ( real_type x = 0.0; x <= 4.5; x += 0.045 )
  {
    for ( integer i = 0; i < ss_a.num_splines() && same; ++i ) same = abs( ss_a( x, i ) - ss_b( x, i ) ) <= 1e-14;
    if ( !same ) break;
  }
  check( same, "both routes give identical spline evaluations on [0, 4.5]" );
}

// ---------------------------------------------------------------------------
// 5. Output direction: SplineSet::eval into GC, then GC -> json
// ---------------------------------------------------------------------------
static void test_eval_output_to_json()
{
  section( "5. SplineSet::eval output through GC -> json" );

  json const j = {
    { "spline_type", { "cubic", "pchip" } },
    { "headers", { "cubic", "pchip" } },
    { "xdata", { 0, 0.9, 2.1, 3, 4.5 } },
    { "ydata", { { 0, 1, 1.99, 2.0, 2.1 }, { 0, 1, 1.99, 2.0, 2.1 } } },
  };

  SplineSet ss;
  ss.build( j.get<GC::GenericContainer>() );

  real_type const      x = 1.5;
  GC::GenericContainer gc_out;
  ss.eval( x, gc_out );

  json const j_out = gc_out;
  check( j_out.is_object(), "eval output converts to a json object" );
  for ( integer i = 0; i < ss.num_splines(); ++i )
  {
    string const name{ ss.header( i ) };
    check( j_out.contains( name ), fmt::format( "json output contains key '{}'", name ) );
    if ( j_out.contains( name ) )
      check_close( j_out[name].get<real_type>(), ss( x, i ), 0.0, fmt::format( "value of '{}' matches ss({}, {})", name, x, i ) );
  }
}

int main()
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║                        TEST N.18                         ║\n"
    "║       (nlohmann::json <-> GenericContainer bridge)       ║\n"
    "╚══════════════════════════════════════════════════════════╝\n" );

  test_json_gc_json_roundtrip();
  test_gc_json_gc_roundtrip();
  test_cubic_spline_from_json();
  test_spline_set_from_json();
  test_eval_output_to_json();

  if ( n_fail == 0 )
  {
    fmt::print(
      fg( fmt::color::cyan ) | fmt::emphasis::bold,
      "\n"
      "╔══════════════════════════════════════════════════════════╗\n"
      "║                    ALL TESTS COMPLETED                   ║\n"
      "╚══════════════════════════════════════════════════════════╝\n\n" );
  }
  else
  {
    fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold, "\n{} CHECK(S) FAILED\n\n", n_fail );
  }

  return n_fail == 0 ? 0 : 1;
}
