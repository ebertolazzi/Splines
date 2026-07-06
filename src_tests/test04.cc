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
#endif

#include "Splines/Splines.hh"
#include "Utils_fmt.hh"

#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::integer;
using Splines::real_type;

// Test problem for Akima interpolation
// Ref. : Hiroshi Akima, Journal of the ACM, Vol. 17, No. 4, October 1970, pages 589-602.

static real_type xx0[] = { 0, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 10 };
static real_type yy0[] = { 10, 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

static real_type xx1[] = { 0, 1, 3, 4, 6, 7, 7, 9, 10, 12, 13, 15 };
static real_type yy1[] = { 10, 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

static real_type xx2[] = { 0, 2, 3, 5, 6, 8, 8, 9, 11, 12, 14, 15 };
static real_type yy2[] = { 10, 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

// RPN 14
static real_type xx3[] = { 7.99, 8.09, 8.19, 8.7, 8.7, 9.2, 10, 12, 15, 20 };
static real_type yy3[] = { 0,        2.76429e-5, 4.37498e-2, 0.169183, 0.169183,
                           0.469428, 0.943740,   0.998636,   0.999919, 0.999994 };

// Titanium
static real_type xx4[] = { 595, 635, 695, 795, 855, 875, 875, 895, 915, 935, 985, 1035, 1075 };
static real_type yy4[] = { 0.644, 0.652, 0.644, 0.694, 0.907, 1.336, 1.336, 2.169, 1.598, 0.916, 0.607, 0.603, 0.608 };

// toolpath
static real_type xx5[] = { 0.11, 0.12, 0.15, 0.16 };
static real_type yy5[] = { 0.0003, 0.0003, 0.0004, 0.0004 };

static integer nn[] = { 12, 12, 12, 10, 13, 4 };

// ===========================================================================
// IMPROVED SUPPORT FUNCTIONS
// ===========================================================================

// Function for central finite differences
real_type finite_diff_central( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{ return ( ss( x + h, i ) - ss( x - h, i ) ) / ( 2 * h ); }

// Function for forward/backward finite differences
real_type finite_diff_forward( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{ return ( ss( x + h, i ) - ss( x, i ) ) / h; }

real_type finite_diff_backward( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{ return ( ss( x, i ) - ss( x - h, i ) ) / h; }

// Adaptive finite difference function
real_type finite_diff_adaptive(
  SplineSet const & ss,
  real_type         x,
  integer           i,
  real_type         xmin,
  real_type         xmax,
  real_type         h = 1e-6 )
{
  if ( x <= xmin + h ) { return finite_diff_forward( ss, x, i, h ); }
  else if ( x >= xmax - h ) { return finite_diff_backward( ss, x, i, h ); }
  else
  {
    return finite_diff_central( ss, x, i, h );
  }
}

// Compute error robustly
pair<real_type, real_type> compute_error( real_type exact, real_type approx, real_type abs_tol = 1e-12 )
{
  real_type abs_err = abs( exact - approx );
  real_type rel_err = 0.0;

  if ( abs( exact ) < abs_tol && abs( approx ) < abs_tol ) { rel_err = 0.0; }
  else
  {
    real_type denom = max( abs( exact ), abs( approx ) );
    if ( denom > abs_tol ) { rel_err = 100.0 * abs_err / denom; }
    else
    {
      rel_err = ( abs_err > abs_tol ) ? 100.0 : 0.0;
    }
  }

  return make_pair( abs_err, min( rel_err, 100.0 ) );
}

// ===========================================================================
// TABLE FUNCTIONS
// ===========================================================================

void print_table_header( vector<string> const & headers )
{
  fmt::print( "┌{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < headers.size() - 1; ++i ) { fmt::print( "┬{}", fmt::format( "{:─^{}}", "", 16 ) ); }
  fmt::print( "┬{:─^{}}┐\n", "", 25 );

  fmt::print( "│{:^12}", headers[0] );
  for ( size_t i = 1; i < headers.size() - 1; ++i ) { fmt::print( "│{:^16}", headers[i] ); }
  fmt::print( "│{:^25}│\n", headers.back() );

  fmt::print( "├{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < headers.size() - 1; ++i ) { fmt::print( "┼{}", fmt::format( "{:─^{}}", "", 16 ) ); }
  fmt::print( "┼{:─^{}}┤\n", "", 25 );
}

template <typename COLOR> void print_table_row( vector<string> const & values, COLOR color )
{
  size_t i = 0;
  fmt::print( "│{:^12}", values[i] );
  for ( ++i; i < values.size() - 2; ++i ) { fmt::print( "│{:^16}", values[i] ); }
  fmt::print( "│" );
  fmt::print( fg( color ), "{:^16}", values[i] );
  fmt::print( "│{:^25}│\n", values.back() );
}

void print_table_footer( size_t ncols )
{
  fmt::print( "└{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < ncols - 1; ++i ) { fmt::print( "┴{}", fmt::format( "{:─^{}}", "", 16 ) ); }
  fmt::print( "┴{:─^{}}┘\n", "", 25 );
}

// ===========================================================================
// DUPLICATE KNOT HANDLING
// ===========================================================================

bool is_duplicate_knot( real_type const * xx, integer npts, real_type x, integer & first_index, real_type eps = 1e-12 )
{
  first_index   = -1;
  integer count = 0;
  for ( integer i = 0; i < npts; ++i )
  {
    if ( abs( xx[i] - x ) < eps )
    {
      if ( first_index == -1 ) first_index = i;
      ++count;
      if ( count > 1 ) return true;
    }
  }
  return false;
}

bool too_close_to_duplicate( real_type const * xx, integer npts, real_type x, real_type safety_margin = 1e-4 )
{
  for ( integer i = 0; i < npts; ++i )
  {
    if ( i > 0 && abs( xx[i] - xx[i - 1] ) < 1e-12 )
    {
      if ( abs( x - xx[i] ) < safety_margin ) { return true; }
    }
  }
  return false;
}

bool is_near_knot( real_type x, real_type const * knots, integer n, real_type eps = 1e-10 )
{
  for ( integer i = 0; i < n; ++i )
  {
    if ( abs( x - knots[i] ) <= eps ) return true;
  }
  return false;
}

bool is_spline_differentiable_at_knots( string const & spline_name )
{ return !( spline_name == "CONSTANT" || spline_name == "LINEAR" ); }

// ===========================================================================
// MAIN FUNCTION
// ===========================================================================

int main()
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║                         TEST N.4                         ║\n"
    "║                  (with duplicate knots)                  ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  SplineSet ss;
  ofstream  file, file_D;

  // Configuration parameters
  const real_type h_fd           = 1e-6;
  const real_type abs_tol        = 1e-8;
  const real_type rel_tol        = 1e-4;
  const real_type knot_tolerance = 1e-10;

  for ( integer k = 0; k < 6; ++k )
  {
    real_type const * xx{ nullptr };
    real_type const * yy{ nullptr };
    string            dataset_name;
    switch ( k )
    {
      case 0:
        xx           = xx0;
        yy           = yy0;
        dataset_name = "Akima Test 0 (with duplicate)";
        break;
      case 1:
        xx           = xx1;
        yy           = yy1;
        dataset_name = "Akima Test 1 (with duplicate)";
        break;
      case 2:
        xx           = xx2;
        yy           = yy2;
        dataset_name = "Akima Test 2 (with duplicate)";
        break;
      case 3:
        xx           = xx3;
        yy           = yy3;
        dataset_name = "RPN 14 (with duplicate)";
        break;
      case 4:
        xx           = xx4;
        yy           = yy4;
        dataset_name = "Titanium (with duplicate)";
        break;
      case 5:
        xx           = xx5;
        yy           = yy5;
        dataset_name = "Toolpath";
        break;
    }

    fmt::print(
      fg( fmt::color::yellow ) | fmt::emphasis::bold,
      "\n"
      "┌──────────────────────────────────────────────────────────┐\n"
      "│ Dataset {}: {:45} │\n"
      "└──────────────────────────────────────────────────────────┘\n",
      k,
      dataset_name );

    string const fname{ fmt::format( "out/SplineSet{}.txt", k ) };
    file.open( fname.data() );

    string const fname_D{ fmt::format( "out/SplineSet{}_D.txt", k ) };
    file_D.open( fname_D.data() );

    real_type const xmin{ xx[0] };
    real_type const xmax{ xx[nn[k] - 1] };
    integer const   npts{ nn[k] };

    // Define 10 splines including all quintic subtypes
    constexpr integer nspl      = 10;
    char const *      headers[] = { "CONSTANT", "LINEAR",        "CUBIC",         "AKIMA",           "VANLEER",
                                    "PCHIP",    "QUINTIC_CUBIC", "QUINTIC_AKIMA", "QUINTIC_VANLEER", "QUINTIC_PCHIP" };

    SplineType1D const stype[] = { SplineType1D::CONSTANT,        SplineType1D::LINEAR,
                                   SplineType1D::CUBIC,           SplineType1D::AKIMA,
                                   SplineType1D::VANLEER,         SplineType1D::PCHIP,
                                   SplineType1D::QUINTIC_CUBIC,   SplineType1D::QUINTIC_AKIMA,
                                   SplineType1D::QUINTIC_VANLEER, SplineType1D::QUINTIC_PCHIP };

    // Prepare data for all splines
    vector<real_type const *> Y( nspl );
    for ( integer i = 0; i < nspl; ++i ) Y[i] = yy;

    // Build spline set
    ss.build( nspl, npts, headers, stype, xx, Y.data() );

    // Write spline values file
    file << "x";
    for ( integer i = 0; i < nspl; ++i ) file << '\t' << ss.header( i );
    file << '\n';
    for ( real_type x{ xmin }; x <= xmax; x += ( xmax - xmin ) / 1000 )
    {
      file << x;
      for ( integer i = 0; i < nspl; ++i ) file << '\t' << ss( x, i );
      file << '\n';
    }
    file.close();

    // Write spline derivatives file
    file_D << "x";
    for ( integer i = 0; i < nspl; ++i ) file_D << '\t' << ss.header( i );
    file_D << '\n';
    for ( real_type x{ xmin }; x <= xmax; x += ( xmax - xmin ) / 1000 )
    {
      file_D << x;
      for ( integer i = 0; i < nspl; ++i ) file_D << '\t' << ss.D( x, i );
      file_D << '\n';
    }
    file_D.close();

    // ========== DERIVATIVE CHECK WITH FINITE DIFFERENCES ==========

    // Identify valid intervals (where x[i] != x[i+1])
    vector<pair<real_type, real_type>> valid_intervals;
    for ( integer i = 0; i < npts - 1; ++i )
    {
      if ( abs( xx[i + 1] - xx[i] ) > 1e-12 ) { valid_intervals.emplace_back( xx[i], xx[i + 1] ); }
    }

    // Generate test points in valid intervals
    vector<real_type> test_points;
    for ( const auto & interval : valid_intervals )
    {
      real_type x1 = interval.first;
      real_type x2 = interval.second;
      real_type dx = ( x2 - x1 ) / 10.0;

      for ( integer j = 1; j < 10; ++j )
      {
        real_type x = x1 + j * dx;
        test_points.push_back( x );
      }
    }

    // Add single knots (not duplicates) if not too close to duplicates
    for ( integer i = 0; i < npts; ++i )
    {
      integer first_index;
      if ( !is_duplicate_knot( xx, npts, xx[i], first_index ) && !too_close_to_duplicate( xx, npts, xx[i] ) )
      {
        test_points.push_back( xx[i] );
      }
    }

    // Sort and remove approximate duplicates
    sort( test_points.begin(), test_points.end() );
    test_points.erase(
      unique( test_points.begin(), test_points.end(), []( real_type a, real_type b ) { return abs( a - b ) < 1e-12; } ),
      test_points.end() );

    // Table for derivative check
    fmt::print(
      fg( fmt::color::green ),
      "\nDerivative check - Dataset {} ({} valid test points)\n",
      k,
      test_points.size() );

    vector<string> table_headers = { "x", "Spline", "D(x)", "FinDiff", "AbsErr", "Note" };
    print_table_header( table_headers );

    integer total_points_tested = 0;
    integer points_with_error   = 0;
    integer points_skipped      = 0;

    for ( integer spline_idx = 0; spline_idx < nspl; ++spline_idx )
    {
      string spline_name    = string( ss.header( spline_idx ) );
      bool   differentiable = is_spline_differentiable_at_knots( spline_name );

      for ( size_t pt_idx = 0; pt_idx < test_points.size(); ++pt_idx )
      {
        real_type x = test_points[pt_idx];

        // Check for duplicate knots
        integer first_index;
        bool    is_dup_knot = is_duplicate_knot( xx, npts, x, first_index );
        bool    is_knot     = is_near_knot( x, xx, npts, knot_tolerance );

        // Skip duplicate knots and points too close to them
        if ( is_dup_knot || too_close_to_duplicate( xx, npts, x ) )
        {
          ++points_skipped;
          continue;
        }

        // Skip non-differentiable splines at knots
        if ( !differentiable && is_knot )
        {
          ++points_skipped;
          continue;
        }

        // Calculate spline derivative
        real_type D_spline = ss.D( x, spline_idx );

        // Calculate finite difference derivative
        real_type D_fd;
        string    note = is_knot ? "knot" : "interior";

        try
        {
          D_fd = finite_diff_adaptive( ss, x, spline_idx, xmin, xmax, h_fd );
        }
        catch ( ... )
        {
          if ( x <= xmin + h_fd ) { D_fd = finite_diff_forward( ss, x, spline_idx, h_fd ); }
          else if ( x >= xmax - h_fd ) { D_fd = finite_diff_backward( ss, x, spline_idx, h_fd ); }
          else
          {
            D_fd = finite_diff_central( ss, x, spline_idx, h_fd );
          }
        }

        // Calculate error
        auto [abs_err, rel_err] = compute_error( D_spline, D_fd, abs_tol );

        ++total_points_tested;
        if ( abs_err > abs_tol && rel_err > rel_tol ) ++points_with_error;

        // Print only: knots, significant errors, or sampled points
        bool should_print = ( is_knot || abs_err > abs_tol || pt_idx % 5 == 0 );

        if ( should_print )
        {
          vector<string> row_values = { fmt::format( "{:.6f}", x ),        spline_name,
                                        fmt::format( "{:.6e}", D_spline ), fmt::format( "{:.6e}", D_fd ),
                                        fmt::format( "{:.2e}", abs_err ),  note };

          // Color coding based on absolute error
          fmt::color color_code;
          if ( abs_err < abs_tol || rel_err < rel_tol ) { color_code = fmt::color::green; }
          else if ( abs_err < 10.0 * abs_tol && rel_err < 10.0 * rel_tol ) { color_code = fmt::color::yellow; }
          else
          {
            color_code = fmt::color::red;
          }

          print_table_row( row_values, color_code );
        }
      }

      // Separator between different splines
      if ( spline_idx < nspl - 1 )
      {
        fmt::print(
          "├{}┼{}┼{}┼{}┼{}┼{}┤\n",
          fmt::format( "{:─^{}}", "", 12 ),
          fmt::format( "{:─^{}}", "", 16 ),
          fmt::format( "{:─^{}}", "", 16 ),
          fmt::format( "{:─^{}}", "", 16 ),
          fmt::format( "{:─^{}}", "", 16 ),
          fmt::format( "{:─^{}}", "", 25 ) );
      }
    }

    print_table_footer( table_headers.size() );

    // Statistics
    fmt::print(
      fg( fmt::color::blue ),
      "\n📊 Statistics Dataset {}:\n"
      "   • Valid test points: {}\n"
      "   • Skipped points: {}\n"
      "   • Points with error > tolerance: {} ({:.1f}%)\n\n",
      k,
      total_points_tested,
      points_skipped,
      points_with_error,
      total_points_tested > 0 ? 100.0 * points_with_error / total_points_tested : 0.0 );
  }

  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║                    ALL TESTS COMPLETED                   ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  return 0;
}
