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

// test09.cc

#include "Splines.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

using namespace Splines;
using namespace std;
using Utils::m_pi;

// Funzione per stampare intestazioni
void print_header( const string & title )
{
  fmt::print(
    fg( fmt::color::steel_blue ) | fmt::emphasis::bold,
    "\n┌{0:─^{2}}┐\n"
    "│{1: ^{2}}│\n"
    "└{0:─^{2}}┘\n",
    "",
    title,
    50 );
}

// Funzione per stampare sottotitoli
void print_subheader( const string & title )
{
  fmt::print(
    fg( fmt::color::cornflower_blue ) | fmt::emphasis::bold,
    "\n├{0:─^{2}}┤\n"
    "│{1: ^{2}}│\n"
    "├{0:─^{2}}┤\n",
    "",
    title,
    50 );
}

// Funzione per stampare risultati in tabella
template <typename T> void print_result_table( const vector<string> & headers, const vector<vector<T>> & data )
{
  vector<size_t> col_widths( headers.size(), 15 );

  // Calcola larghezze colonne
  for ( size_t i = 0; i < headers.size(); i++ ) { col_widths[i] = max( col_widths[i], headers[i].size() ); }
  for ( const auto & row : data )
  {
    for ( size_t i = 0; i < row.size(); i++ )
    {
      col_widths[i] = max( col_widths[i], fmt::format( "{}", row[i] ).size() );
    }
  }

  // Stampa intestazione
  fmt::print( "┌" );
  for ( size_t i = 0; i < col_widths.size(); i++ )
  {
    fmt::print( "{0:─^{1}}", "", col_widths[i] + 2 );
    if ( i < col_widths.size() - 1 ) fmt::print( "┬" );
  }
  fmt::print( "┐\n" );

  fmt::print( "│" );
  for ( size_t i = 0; i < headers.size(); i++ )
  {
    fmt::print( fg( fmt::color::gold ) | fmt::emphasis::bold, " {:^{}} ", headers[i], col_widths[i] );
    fmt::print( "│" );
  }
  fmt::print( "\n" );

  // Separatore
  fmt::print( "├" );
  for ( size_t i = 0; i < col_widths.size(); i++ )
  {
    fmt::print( "{0:─^{1}}", "", col_widths[i] + 2 );
    if ( i < col_widths.size() - 1 ) fmt::print( "┼" );
  }
  fmt::print( "┤\n" );

  // Dati
  for ( const auto & row : data )
  {
    fmt::print( "│" );
    for ( size_t i = 0; i < row.size(); i++ )
    {
      if ( i == 0 && headers[0] == "x" )
      {
        fmt::print( fg( fmt::color::light_green ), " {:^{}.4f} ", row[i], col_widths[i] );
      }
      else
      {
        fmt::print( " {:^{}.6f} ", row[i], col_widths[i] );
      }
      fmt::print( "│" );
    }
    fmt::print( "\n" );
  }

  // Chiusura tabella
  fmt::print( "└" );
  for ( size_t i = 0; i < col_widths.size(); i++ )
  {
    fmt::print( "{0:─^{1}}", "", col_widths[i] + 2 );
    if ( i < col_widths.size() - 1 ) fmt::print( "┴" );
  }
  fmt::print( "┘\n" );
}

void test_constant_spline()
{
  print_header( "📊 Testing ConstantSpline" );

  ConstantSpline cs( "Constant" );

  // Build from arrays
  real_type x[] = { 0, 1, 2, 3, 4 };
  real_type y[] = { 1, 2, 1, 3, 2 };
  cs.build( x, y, 5 );

  fmt::print( fg( fmt::color::light_blue ), "\n📈 Spline Information:\n" );
  fmt::print( "  Type: {}\n", cs.type_name() );
  fmt::print( "  Points: {}\n", cs.num_points() );
  fmt::print( "  X range: [{:.2f}, {:.2f}]\n", cs.x_min(), cs.x_max() );
  fmt::print( "  Y range: [{:.2f}, {:.2f}]\n", cs.y_min(), cs.y_max() );

  // Test evaluation
  fmt::print( fg( fmt::color::light_blue ), "\n🔍 Evaluations:\n" );

  vector<vector<real_type>> eval_data;
  for ( real_type xi = -1; xi <= 5; xi += 0.5 )
  {
    try
    {
      eval_data.push_back( { xi, cs.eval( xi ), cs.D( xi ), cs.DD( xi ) } );
    }
    catch ( exception & e )
    {
      fmt::print( fg( fmt::color::red ), "  ✗ f({:.1f}) ERROR: {}\n", xi, e.what() );
    }
  }

  print_result_table( { "x", "f(x)", "f'(x)", "f''(x)" }, eval_data );

  // Test min/max
  integer   i_min, i_max;
  real_type x_min, y_min, x_max, y_max;
  cs.y_min_max( i_min, x_min, y_min, i_max, x_max, y_max );

  fmt::print( fg( fmt::color::light_green ), "\n📊 Min/Max Analysis:\n" );
  fmt::print( "  Minimum: f({:.2f}) = {:.2f} {} at point {}\n", x_min, y_min, "🔽", i_min );
  fmt::print( "  Maximum: f({:.2f}) = {:.2f} {} at point {}\n", x_max, y_max, "🔼", i_max );
}

void test_linear_spline()
{
  print_header( "📈 Testing LinearSpline" );

  LinearSpline ls( "Linear" );

  // Build from vectors
  vector<real_type> x = { 0, 1, 2, 3, 4 };
  vector<real_type> y = { 0, 1, 0, 2, 1 };
  ls.build( x, y );

  fmt::print( fg( fmt::color::light_blue ), "\n📈 Spline Information:\n" );
  fmt::print( "  Type: {}\n", ls.type_name() );
  fmt::print( "  Points: {}\n", ls.num_points() );

  // Test evaluation
  fmt::print( fg( fmt::color::light_blue ), "\n🔍 Evaluations:\n" );

  vector<vector<real_type>> eval_data;
  for ( real_type xi = 0; xi <= 4; xi += 0.5 )
  {
    eval_data.push_back( { xi, ls.eval( xi ), ls.D( xi ), ls.DD( xi ) } );
  }

  print_result_table( { "x", "f(x)", "f'(x)", "f''(x)" }, eval_data );

  // Test different boundary conditions
  fmt::print( fg( fmt::color::light_blue ), "\n⚙️ Boundary Condition Tests:\n" );

  ls.make_closed();
  fmt::print( "  {} Closed spline at x=5: {:.6f}\n", "🔄", ls.eval( 5 ) );
  ls.make_opened();

  ls.make_bounded();
  try
  {
    fmt::print( "  {} Bounded spline at x=5: {:.6f}\n", "⏹️", ls.eval( 5 ) );
  }
  catch ( exception & e )
  {
    fmt::print( fg( fmt::color::orange ), "  ⚠️ Bounded spline at x=5: ERROR: {}\n", e.what() );
  }
  ls.make_unbounded();
}

void test_cubic_spline()
{
  print_header( "🔷 Testing CubicSpline" );

  CubicSpline cs( "Cubic" );

  // Set boundary conditions
  cs.set_initial_BC( CubicSpline_BC::NATURAL );
  cs.set_final_BC( CubicSpline_BC::NATURAL );

  // Build with data
  real_type x[] = { 0, 1, 2, 3, 4 };
  real_type y[] = { 0, 1, 0, 1, 0 };
  cs.build( x, y, 5 );

  fmt::print( fg( fmt::color::light_blue ), "\n📈 Spline Information:\n" );
  fmt::print( "  Type: {}\n", cs.type_name() );
  fmt::print( "  Points: {}\n", cs.num_points() );
  fmt::print( "  Boundary Conditions: Natural (both ends)\n" );

  // Test evaluation
  fmt::print( fg( fmt::color::light_blue ), "\n🔍 Evaluations:\n" );

  vector<vector<real_type>> eval_data;
  for ( real_type xi = 0; xi <= 4; xi += 0.25 )
  {
    eval_data.push_back( { xi, cs.eval( xi ), cs.D( xi ), cs.DD( xi ), cs.DDD( xi ) } );
  }

  print_result_table( { "x", "f(x)", "f'(x)", "f''(x)", "f'''(x)" }, eval_data );

  // Test with known derivatives
  print_subheader( "With Known Derivatives" );

  real_type yp[] = { 1, 0, -1, 0, 1 };
  cs.build( x, y, yp, 5 );

  fmt::print( fg( fmt::color::light_blue ), "\n🔍 Evaluations with prescribed derivatives:\n" );

  vector<vector<real_type>> deriv_data;
  for ( real_type xi = 0; xi <= 4; xi += 0.5 )
  {
    real_type dd[2];
    cs.D( xi, dd );
    deriv_data.push_back( { xi, dd[0], dd[1] } );
  }

  print_result_table( { "x", "f(x)", "f'(x)" }, deriv_data );
}

void test_akima_spline()
{
  print_header( "🌀 Testing AkimaSpline" );

  AkimaSpline as( "Akima" );

  vector<real_type> x = { 0, 1, 2, 3, 4, 5 };
  vector<real_type> y = { 0, 2, 1, 3, 2, 1 };
  as.build( x, y );

  fmt::print( fg( fmt::color::light_blue ), "\n📈 Spline Information:\n" );
  fmt::print( "  Type: {}\n", as.type_name() );
  fmt::print( "  Points: {}\n", x.size() );

  // Test evaluation
  fmt::print( fg( fmt::color::light_blue ), "\n🔍 Evaluations:\n" );

  vector<vector<real_type>> eval_data;
  for ( real_type xi = 0; xi <= 5; xi += 0.5 ) { eval_data.push_back( { xi, as.eval( xi ), as.D( xi ) } ); }

  print_result_table( { "x", "f(x)", "f'(x)" }, eval_data );
}

void test_vanleer_spline()
{
  print_header( "⚙️ Testing VanLeerSpline" );

  VanLeerSpline bs( "VanLeer" );

  real_type x[] = { 0, 1, 2, 3, 4 };
  real_type y[] = { 0, 1, 0.5, 1.5, 1 };
  bs.build( x, y, 5 );

  fmt::print( fg( fmt::color::light_blue ), "\n📈 Spline Information:\n" );
  fmt::print( "  Type: {}\n", bs.type_name() );

  vector<vector<real_type>> eval_data;
  for ( real_type xi = 0; xi <= 4; xi += 0.5 ) { eval_data.push_back( { xi, bs.eval( xi ), bs.D( xi ) } ); }

  print_result_table( { "x", "f(x)", "f'(x)" }, eval_data );
}

void test_pchip_spline()
{
  print_header( "📐 Testing PchipSpline" );

  PchipSpline ps( "Pchip" );

  vector<real_type> x = { 0, 1, 2, 3, 4 };
  vector<real_type> y = { 0, 3, 1, 2, 2 };
  ps.build( x, y );

  fmt::print( fg( fmt::color::light_blue ), "\n📈 Spline Information:\n" );
  fmt::print( "  Type: {}\n", ps.type_name() );

  vector<vector<real_type>> eval_data;
  for ( real_type xi = 0; xi <= 4; xi += 0.5 ) { eval_data.push_back( { xi, ps.eval( xi ), ps.D( xi ) } ); }

  print_result_table( { "x", "f(x)", "f'(x)" }, eval_data );
}

void test_hermite_spline()
{
  print_header( "🔺 Testing HermiteSpline" );

  HermiteSpline hs( "Hermite" );

  real_type x[]  = { 0, 1, 2, 3, 4 };
  real_type y[]  = { 0, 1, 0, 1, 0 };
  real_type yp[] = { 1, 0, -1, 0, 1 };
  hs.build( x, y, yp, 5 );

  fmt::print( fg( fmt::color::light_blue ), "\n📈 Spline Information:\n" );
  fmt::print( "  Type: {}\n", hs.type_name() );

  vector<vector<real_type>> eval_data;
  for ( real_type xi = 0; xi <= 4; xi += 0.5 )
  {
    eval_data.push_back( { xi, hs.eval( xi ), hs.D( xi ), hs.DD( xi ) } );
  }

  print_result_table( { "x", "f(x)", "f'(x)", "f''(x)" }, eval_data );
}

void test_quintic_spline()
{
  print_header( "⚙️ Testing QuinticSpline" );

  // Test different subtypes
  vector<Spline_sub_type> types = { Spline_sub_type::CUBIC,
                                    Spline_sub_type::PCHIP,
                                    Spline_sub_type::AKIMA,
                                    Spline_sub_type::VANLEER };

  vector<string>     type_names  = { "CUBIC", "PCHIP", "AKIMA", "VANLEER" };
  vector<fmt::color> type_colors = { fmt::color::light_blue,
                                     fmt::color::light_green,
                                     fmt::color::light_pink,
                                     fmt::color::light_yellow };

  vector<real_type> x = { 0, 1, 2, 3, 4 };
  vector<real_type> y = { 0, 1, 0.5, 1.5, 1 };

  for ( size_t i = 0; i < types.size(); i++ )
  {
    auto type = types[i];

    QuinticSpline qs( type, "Quintic" );
    // qs.set_quintic_type( type );
    qs.build( x, y );

    fmt::print(
      fg( type_colors[i] ) | fmt::emphasis::bold,
      "\n📊 Subtype: {} (Type {})\n",
      type_names[i],
      to_string( type ) );

    vector<vector<real_type>> eval_data;
    for ( real_type xi = 0; xi <= 4; xi += 0.5 )
    {
      eval_data.push_back( { xi, qs.eval( xi ), qs.D( xi ), qs.DD( xi ), qs.DDD( xi ) } );
    }

    print_result_table( { "x", "f(x)", "f'(x)", "f''(x)", "f'''(x)" }, eval_data );
  }
}

void test_spline1d()
{
  print_header( "📊 Testing Spline1D Base Class" );

  vector<real_type> x = { 0, 1, 2, 3, 4 };
  vector<real_type> y = { 0, 1, 0, 1, 0 };

  vector<SplineType1D> types = { SplineType1D::LINEAR, SplineType1D::CUBIC, SplineType1D::AKIMA, SplineType1D::PCHIP };

  vector<string>     type_names = { "LINEAR", "CUBIC", "AKIMA", "PCHIP" };
  vector<fmt::color> colors     = { fmt::color::light_green,
                                    fmt::color::light_blue,
                                    fmt::color::light_pink,
                                    fmt::color::light_yellow };

  vector<vector<real_type>> comparison_data;

  for ( size_t i = 0; i < types.size(); i++ )
  {
    auto      type  = types[i];
    real_type value = 0;

    if ( type == SplineType1D::LINEAR )
    {
      LinearSpline ls( "Linear" );
      ls.build( x, y );
      value = ls.eval( 2.5 );
    }
    else if ( type == SplineType1D::CUBIC )
    {
      CubicSpline cs( "Cubic" );
      cs.build( x, y );
      value = cs.eval( 2.5 );
    }
    else if ( type == SplineType1D::AKIMA )
    {
      AkimaSpline as( "Akima" );
      as.build( x, y );
      value = as.eval( 2.5 );
    }
    else if ( type == SplineType1D::PCHIP )
    {
      PchipSpline ps( "Pchip" );
      ps.build( x, y );
      value = ps.eval( 2.5 );
    }

    comparison_data.push_back( { static_cast<real_type>( i ), value } );

    fmt::print( fg( colors[i] ), "  {} f(2.5) = {:.6f}\n", type_names[i], value );
  }

  fmt::print(
    fg( fmt::color::light_blue ) | fmt::emphasis::bold,
    "\n📈 Comparison of Different Spline Types at x = 2.5:\n" );

  print_result_table( { "Type Index", "f(2.5)" }, comparison_data );
}

void test_spline_set()
{
  print_header( "🧩 Testing SplineSet" );

  SplineSet ss( "TestSet" );

  integer npts = 5;
  integer nspl = 3;

  real_type X[] = { 0, 1, 2, 3, 4 };

  real_type Y0[] = { 0, 0.2, 0.5, 0.8, 1 };
  real_type Y1[] = { 1, 0, 1, 0, 1 };
  real_type Y2[] = { 0, 0.5, 1, 0.5, 0 };

  real_type * Y[] = { Y0, Y1, Y2 };

  char const * headers[] = { "Curve1 📈", "Curve2 📉", "Curve3 🔄" };
  SplineType1D stype[]   = { SplineType1D::CUBIC, SplineType1D::CUBIC, SplineType1D::CUBIC };

  ss.build( nspl, npts, headers, stype, X, Y );

  fmt::print( fg( fmt::color::light_blue ), "\n📊 Spline Set Information:\n" );
  fmt::print( "  Name: {}\n", ss.name() );
  fmt::print( "  Number of splines: {}\n", ss.num_splines() );
  fmt::print( "  Number of points: {}\n", ss.num_points() );

  // Evaluate all splines at a point
  real_type x = 2.5;
  fmt::print( fg( fmt::color::light_blue ), "\n🔍 Evaluating at x = {:.2f}:\n", x );

  vector<vector<real_type>> eval_data;
  for ( integer i = 0; i < nspl; i++ )
  {
    eval_data.push_back( { static_cast<real_type>( i ), ss.eval( x, i ) } );
    fmt::print( "  {}: {:.6f}\n", ss.header( i ), ss.eval( x, i ) );
  }

  // Evaluate to vector
  vector<real_type> vals;
  ss.eval( x, vals );
  fmt::print( fg( fmt::color::light_blue ), "\n📋 Vector evaluation:\n" );
  for ( size_t i = 0; i < vals.size(); i++ ) { fmt::print( "  Curve {}: {:.6f}\n", i + 1, vals[i] ); }

  // Test with independent spline
  fmt::print( fg( fmt::color::light_blue ), "\n🔄 Using Curve1 as independent variable:\n" );
  real_type zeta = 0.5;
  real_type x_val;
  try
  {
    real_type result = ss.eval2( zeta, 0, x_val, 1 );
    fmt::print(
      fg( fmt::color::light_green ),
      "  ✅ When Curve1 = {:.2f}, x = {:.6f}, Curve2 = {:.6f}\n",
      zeta,
      x_val,
      result );
  }
  catch ( exception & e )
  {
    fmt::print( fg( fmt::color::orange ), "  ⚠️ ERROR: {}\n", e.what() );
  }
}

void test_spline_vec()
{
  print_header( "🎯 Testing SplineVec (3D Curve)" );

  SplineVec sv( "VectorSpline" );

  constexpr integer dim  = 2;
  constexpr integer npts = 20;

  real_type Y0[npts];
  real_type Y1[npts];
  for ( integer k = 0; k < npts; ++k )
  {
    Y0[k] = cos( k * 0.1 * m_pi );
    Y1[k] = sin( k * 0.1 * m_pi );
  }

  real_type * Y[] = { Y0, Y1 };

  sv.setup( dim, npts, Y );
  sv.set_knots_chord_length();
  sv.catmull_rom();

  fmt::print( fg( fmt::color::light_blue ), "\n📊 Vector Spline Information:\n" );
  fmt::print( "  Dimension: {}\n", sv.dimension() );
  fmt::print( "  Points: {}\n", sv.num_points() );
  fmt::print( "  Parameterization: Chord Length\n" );
  fmt::print( "  Spline Type: Catmull-Rom\n" );

  // Evaluate at points
  fmt::print( fg( fmt::color::light_blue ), "\n🔍 Evaluating along curve:\n" );

  vector<vector<real_type>> eval_data;
  for ( real_type t = sv.x_min(); t <= sv.x_max(); t += 0.1 )
  {
    vector<real_type> vals;
    sv.eval( t, vals );
    eval_data.push_back( { t, vals[0], vals[1] } );

    fmt::print( "  t = {:.1f}: ({:.4f}, {:.4f})\n", t, vals[0], vals[1] );
  }

  fmt::print( fg( fmt::color::light_blue ), "\n📊 Curvature Analysis:\n" );
  for ( real_type t = sv.x_min(); t <= sv.x_max(); t += 0.1 )
  {
    fmt::print( "  t = {:.1f}: κ = {:.6f}\n", t, sv.curvature( t ) );
  }
}

void test_spline1dblend()
{
  print_header( "🌈 Testing Spline1Dblend" );

  Spline1Dblend blend( "BlendedSpline" );

  vector<real_type> x0 = { 0, 1, 2, 3, 4 };
  vector<real_type> y0 = { 0, 1, 0, 1, 0 };

  vector<real_type> x1 = { 0, 1, 2, 3, 4 };
  vector<real_type> y1 = { 0, 0, 1, 0, 0 };

  blend.build( SplineType1D::CUBIC, x0.data(), y0.data(), 5, SplineType1D::CUBIC, x1.data(), y1.data(), 5 );

  fmt::print( fg( fmt::color::light_blue ), "\n📊 Blended Spline Information:\n" );
  fmt::print( "  Spline 0 type: {}\n", blend.type_name0() );
  fmt::print( "  Spline 1 type: {}\n", blend.type_name1() );

  // Test different blend factors
  vector<real_type> s_values = { 0, 0.25, 0.5, 0.75, 1 };
  real_type         x        = 2.5;

  fmt::print( fg( fmt::color::light_blue ), "\n🔍 Blending at x = {:.2f}:\n", x );

  vector<vector<real_type>> blend_data;
  for ( real_type s : s_values )
  {
    real_type f  = blend.eval( x, s );
    real_type df = blend.D( x, s );
    blend_data.push_back( { s, f, df } );

    fmt::print( "  s = {:.2f}: f({:.2f}) = {:.6f}, f'({:.2f}) = {:.6f}\n", s, x, f, x, df );
  }

  print_result_table( { "Blend s", "f(x)", "f'(x)" }, blend_data );
}

void test_utility_functions()
{
  print_header( "🔧 Testing Utility Functions" );

  integer   dim    = 2;
  integer   npts   = 5;
  real_type pnts[] = { 0, 0, 1, 1, 2, 0, 3, -1, 4, 0 };

  real_type t[5];

  fmt::print( fg( fmt::color::light_blue ), "\n📐 Node Distribution Methods:\n" );

  // Uniform distribution
  fmt::print( fg( fmt::color::light_green ), "\n📏 Uniform Distribution:\n" );
  uniform( dim, npts, pnts, dim, t );
  vector<vector<real_type>> uniform_data;
  for ( int i = 0; i < npts; i++ ) { uniform_data.push_back( { static_cast<real_type>( i ), t[i] } ); }
  print_result_table( { "i", "t[i]" }, uniform_data );

  // Chordal distribution
  fmt::print( fg( fmt::color::light_blue ), "\n📐 Chordal Distribution:\n" );
  chordal( dim, npts, pnts, dim, t );
  vector<vector<real_type>> chordal_data;
  for ( int i = 0; i < npts; i++ ) { chordal_data.push_back( { static_cast<real_type>( i ), t[i] } ); }
  print_result_table( { "i", "t[i]" }, chordal_data );

  // Centripetal distribution
  fmt::print( fg( fmt::color::light_pink ), "\n🌀 Centripetal Distribution (α=0.5):\n" );
  centripetal( dim, npts, pnts, dim, 0.5, t );
  vector<vector<real_type>> centripetal_data;
  for ( int i = 0; i < npts; i++ ) { centripetal_data.push_back( { static_cast<real_type>( i ), t[i] } ); }
  print_result_table( { "i", "t[i]" }, centripetal_data );

  // Test Hermite functions
  fmt::print( fg( fmt::color::light_yellow ), "\n🔺 Hermite Basis Functions (H=1):\n" );
  real_type base[4];
  Splines::Hermite3( 0.5, 1.0, base );

  vector<vector<real_type>> hermite_data;
  hermite_data.push_back( { 0.5, base[0], base[1], base[2], base[3] } );

  fmt::print(
    "  h₀(0.5) = {:.6f}, h₁(0.5) = {:.6f}, h₂(0.5) = {:.6f}, h₃(0.5) = {:.6f}\n",
    base[0],
    base[1],
    base[2],
    base[3] );

  print_result_table( { "ξ", "h₀(ξ)", "h₁(ξ)", "h₂(ξ)", "h₃(ξ)" }, hermite_data );
}

int main()
{
  auto start_time = chrono::high_resolution_clock::now();

  try
  {
    fmt::print(
      fg( fmt::color::steel_blue ) | fmt::emphasis::bold,
      "┌────────────────────────────────────────────────────────────┐\n"
      "│               🚀 Testing Splines 1D Library                │\n"
      "└────────────────────────────────────────────────────────────┘\n\n" );

    // Test basic spline types
    using DATA_type         = pair<string, function<void()>>;
    vector<DATA_type> tests = { { "Constant Spline", test_constant_spline },
                                { "Linear Spline", test_linear_spline },
                                { "Cubic Spline", test_cubic_spline },
                                { "Akima Spline", test_akima_spline },
                                { "VanLeer Spline", test_vanleer_spline },
                                { "Pchip Spline", test_pchip_spline },
                                { "Hermite Spline", test_hermite_spline },
                                { "Quintic Spline", test_quintic_spline },
                                { "Spline1D Base", test_spline1d },
                                { "Spline Set", test_spline_set },
                                { "Spline Vector", test_spline_vec },
                                { "Spline Blend", test_spline1dblend },
                                { "Utility Functions", test_utility_functions } };

    int test_count   = 0;
    int passed_count = 0;

    for ( const auto & [name, test_func] : tests )
    {
      test_count++;
      fmt::print( fg( fmt::color::dark_gray ), "\n[{}/{}] Running {}...\n", test_count, tests.size(), name );

      try
      {
        test_func();
        passed_count++;
        fmt::print( fg( fmt::color::light_green ), "  ✅ {} completed successfully\n", name );
      }
      catch ( const exception & e )
      {
        fmt::print( fg( fmt::color::red ), "  ❌ {} failed: {}\n", name, e.what() );
      }
    }

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration<double>( end_time - start_time );

    fmt::print(
      fg( fmt::color::steel_blue ) | fmt::emphasis::bold,
      "\n┌{0:─^{2}}┐\n"
      "│{1: ^{2}}│\n"
      "├{0:─^{2}}┤\n"
      "│{3: ^{2}}│\n"
      "│{4: ^{2}}│\n"
      "└{0:─^{2}}┘\n",
      "",
      "📊 Test Summary",
      60,
      fmt::format( "Tests: {}/{} passed", passed_count, test_count ),
      fmt::format( "Time: {:.3f} seconds", duration.count() ) );

    if ( passed_count == test_count )
    {
      fmt::print( fg( fmt::color::light_green ) | fmt::emphasis::bold, "\n🎉 All tests completed successfully!\n" );
    }
    else
    {
      fmt::print(
        fg( fmt::color::orange ) | fmt::emphasis::bold,
        "\n⚠️  {} test(s) failed!\n",
        test_count - passed_count );
      return 1;
    }
  }
  catch ( exception & e )
  {
    fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold, "\n❌ Critical Error: {}\n", e.what() );
    return 1;
  }

  return 0;
}
