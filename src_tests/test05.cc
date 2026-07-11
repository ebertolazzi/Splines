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

#include "Splines.hh"
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

// monotone
static real_type xx[] = { 0, 0.9, 2.1, 3, 4.5 };
static real_type yy[] = { 0, 1, 1.1, 2.0, 2.1 };

static integer npt = 5;

// Function for central finite differences
real_type finite_diff_central( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{ return ( ss( x + h, i ) - ss( x - h, i ) ) / ( 2 * h ); }

// Function for forward/backward finite differences
real_type finite_diff_forward( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{
  real_type si   = ss( x, i );
  real_type hh   = h * 0.5;
  real_type D_h  = ( ss( x + h, i ) - si ) / h;
  real_type D_h2 = ( ss( x + hh, i ) - si ) / hh;
  return 2.0 * D_h2 - D_h;
}

real_type finite_diff_backward( SplineSet const & ss, real_type x, integer i, real_type h = 1e-6 )
{
  real_type si   = ss( x, i );
  real_type hh   = h * 0.5;
  real_type D_h  = ( si - ss( x - h, i ) ) / h;
  real_type D_h2 = ( si - ss( x - hh, i ) ) / hh;
  return 2.0 * D_h2 - D_h;
}

// Function for central finite differences for eval2 (inverse evaluation)
real_type finite_diff_central_eval2(
  SplineSet const & ss,
  integer           spline_idx,
  real_type         y,
  integer           i,
  real_type         h = 1e-6 )
{
  real_type val_plus[8], val_minus[8], x;
  ss.eval2( spline_idx, y + h, x, val_plus );
  ss.eval2( spline_idx, y - h, x, val_minus );
  return ( val_plus[i] - val_minus[i] ) / ( 2 * h );
}

// Print table with Unicode borders and colors
void print_table_header( vector<string> const & headers )
{
  // Top line
  fmt::print( "┌{}", fmt::format( "{:─^{}}", "", 12 ) );
  for ( size_t i = 1; i < headers.size() - 1; ++i ) { fmt::print( "┬{}", fmt::format( "{:─^{}}", "", 16 ) ); }
  fmt::print( "┬{:─^{}}┐\n", "", 25 );

  // Headers
  fmt::print( "│{:^12}", headers[0] );
  for ( size_t i = 1; i < headers.size() - 1; ++i ) { fmt::print( "│{:^16}", headers[i] ); }
  fmt::print( "│{:^25}│\n", headers.back() );

  // Separator
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

int main()
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║                         TEST N.5                         ║\n"
    "║                  (Monotone Spline Test)                  ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  SplineSet ss;
  ofstream  file, file_D, fileR, fileR_D;

  file.open( "out/SplineSet.txt" );
  file_D.open( "out/SplineSet_D.txt" );
  fileR.open( "out/SplineSetR.txt" );
  fileR_D.open( "out/SplineSetR_D.txt" );

  real_type xmin = xx[0];
  real_type xmax = xx[npt - 1];

  integer   nspl = 7;
  integer   npts = npt;
  real_type val[8], val_D[8];

  char const * headers[] = { "CONSTANT", "LINEAR",  "CUBIC",         "AKIMA",           "VANLEER",
                             "PCHIP",    "QUINTIC", "QUINTIC_AKIMA", "QUINTIC_VANLEER", "QUINTIC_PCHIP" };

  constexpr SplineType1D stype[]{ SplineType1D::CONSTANT,        SplineType1D::LINEAR,
                                  SplineType1D::CUBIC,           SplineType1D::AKIMA,
                                  SplineType1D::VANLEER,         SplineType1D::PCHIP,
                                  SplineType1D::QUINTIC_CUBIC,   SplineType1D::QUINTIC_AKIMA,
                                  SplineType1D::QUINTIC_VANLEER, SplineType1D::QUINTIC_PCHIP };

  Utils::Malloc<real_type> mem( "test5" );
  mem.allocate( npts );
  real_type * YpZero = mem( npts );
  std::fill_n( YpZero, npts, 0 );

  real_type const * Y[] = { yy, yy, yy, yy, yy, yy, yy, yy };

  ss.build( nspl, npts, headers, stype, xx, Y );
  ss.info( cout );

  // Display spline positions with colors
  fmt::print( fg( fmt::color::yellow ), "\nSpline positions in SplineSet:\n" );
  for ( integer i = 0; i < nspl; ++i )
  {
    fmt::print( "  {}: position = {}\n", ss.header( i ), ss.get_position( ss.header( i ) ) );
  }

  // Write spline values file
  file << "x";
  file_D << "x";
  for ( integer i = 0; i < nspl; ++i )
  {
    file << '\t' << ss.header( i );
    file_D << '\t' << ss.header( i );
  }
  file << '\n';
  file_D << '\n';
  for ( real_type x{ xmin }; x <= xmax; x += ( xmax - xmin ) / 1000 )
  {
    file << x;
    file_D << x;
    ss.eval( x, val );
    ss.eval_D( x, val_D );
    for ( integer i = 0; i < nspl; ++i )
    {
      file << '\t' << val[i];
      file_D << '\t' << val_D[i];
    }
    file << '\n';
    file_D << '\n';
  }
  file.close();
  file_D.close();

  // ========== DERIVATIVE CHECK WITH FINITE DIFFERENCES ==========

  // Generate test points: knots + interior points
  vector<real_type> test_points;

  // Add all knots
  for ( integer i = 0; i < npts; ++i ) { test_points.push_back( xx[i] ); }

  // Add interior points between knots
  for ( integer i = 0; i < npts - 1; ++i )
  {
    real_type x1 = xx[i];
    real_type x2 = xx[i + 1];
    // Add 3 interior points per interval
    test_points.push_back( x1 + 0.25 * ( x2 - x1 ) );
    test_points.push_back( x1 + 0.5 * ( x2 - x1 ) );
    test_points.push_back( x1 + 0.75 * ( x2 - x1 ) );
  }

  // Sort and remove duplicates
  sort( test_points.begin(), test_points.end() );
  test_points.erase(
    unique( test_points.begin(), test_points.end(), []( real_type a, real_type b ) { return abs( a - b ) < 1e-12; } ),
    test_points.end() );

  // Table for derivative check
  fmt::print(
    fg( fmt::color::green ),
    "\nDerivative check - Direct evaluation ({} test points)\n",
    test_points.size() );

  vector<string> table_headers = { "x", "Spline", "D(x)", "FinDiff", "Rel. Err%", "Note" };
  print_table_header( table_headers );

  real_type const h                   = 1e-5;
  real_type const tol                 = 1e-5;
  integer         total_points_tested = 0;
  integer         points_with_error   = 0;

  for ( integer spline_idx = 0; spline_idx < nspl; ++spline_idx )
  {
    string spline_name = string( ss.header( spline_idx ) );

    for ( size_t pt_idx = 0; pt_idx < test_points.size(); ++pt_idx )
    {
      real_type x = test_points[pt_idx];

      // Determine if it's a knot
      bool is_knot = ( find( xx, xx + npts, x ) != xx + npts );

      // Calculate spline derivative
      real_type D_spline = ss.D( x, spline_idx );

      // Calculate appropriate finite difference
      real_type D_fd;
      string    note = is_knot ? "knot" : "interior";

      // For knots, use appropriate one-sided finite differences
      if ( is_knot )
      {
        if ( x <= xmin + h )
        {
          D_fd = finite_diff_forward( ss, x, spline_idx, h );
          note += " (left boundary)";
        }
        else if ( x >= xmax - h )
        {
          D_fd = finite_diff_backward( ss, x, spline_idx, h );
          note += " (right boundary)";
        }
        else
        {
          // For interior knots, check if it's start or end of interval
          bool is_start_of_interval = false;
          for ( integer i = 0; i < npts - 1; ++i )
          {
            if ( abs( xx[i] - x ) < 1e-12 && abs( xx[i + 1] - xx[i] ) > 1e-12 )
            {
              is_start_of_interval = true;
              break;
            }
          }
          if ( is_start_of_interval )
          {
            D_fd = finite_diff_forward( ss, x, spline_idx, h );
            note += " (interval start)";
          }
          else
          {
            D_fd = finite_diff_backward( ss, x, spline_idx, h );
            note += " (interval end)";
          }
        }
      }
      else
      {
        // For interior points, use central finite difference
        D_fd = finite_diff_central( ss, x, spline_idx, h );
      }

      // Calculate relative error
      real_type abs_err = abs( D_spline - D_fd );
      real_type rel_err = 100.0 * abs_err / max( { abs( D_spline ), abs( D_fd ), 1e-4 } );

      ++total_points_tested;
      if ( rel_err > tol ) ++points_with_error;

      // Print only if: knot, significant error
      bool should_print = ( is_knot || rel_err > tol );

      if ( should_print )
      {
        vector<string> row_values = { fmt::format( "{:.6f}", x ),        spline_name,
                                      fmt::format( "{:.6f}", D_spline ), fmt::format( "{:.6f}", D_fd ),
                                      fmt::format( "{:.2e}", rel_err ),  note };
        auto color = rel_err < tol ? fmt::color::green : ( rel_err < 100 * tol ? fmt::color::yellow : fmt::color::red );
        print_table_row( row_values, color );
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
    "\n📊 Statistics for direct evaluation:\n"
    "   • Total test points: {}\n"
    "   • Points with error > {:.0e}: {} ({:.1f}%)\n\n",
    total_points_tested,
    tol,
    points_with_error,
    total_points_tested > 0 ? 100.0 * points_with_error / total_points_tested : 0.0 );

  // ========== INVERSE EVALUATION (eval2) ==========

  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "\n"
    "┌──────────────────────────────────────────────────────────┐\n"
    "│               Inverse Evaluation (eval2)                 │\n"
    "└──────────────────────────────────────────────────────────┘\n" );

  real_type ymin = yy[0];
  real_type ymax = yy[npt - 1];

  // Write inverse spline values file
  fileR << "x";
  fileR_D << "x";
  for ( integer i = 0; i < nspl; ++i )
  {
    fileR << '\t' << ss.header( i );
    fileR_D << '\t' << ss.header( i );
  }
  fileR << '\n';
  fileR_D << '\n';

  // Generate test points for inverse evaluation
  vector<real_type> test_points_inverse;

  // Add all y-values (original y data)
  for ( integer i = 0; i < npts; ++i ) { test_points_inverse.push_back( yy[i] ); }

  // Add interior points between y-values
  for ( integer i = 0; i < npts - 1; ++i )
  {
    real_type y1 = yy[i];
    real_type y2 = yy[i + 1];
    // Add 3 interior points per interval
    test_points_inverse.push_back( y1 + 0.25 * ( y2 - y1 ) );
    test_points_inverse.push_back( y1 + 0.5 * ( y2 - y1 ) );
    test_points_inverse.push_back( y1 + 0.75 * ( y2 - y1 ) );
  }

  // Sort and remove duplicates
  sort( test_points_inverse.begin(), test_points_inverse.end() );
  test_points_inverse.erase(
    unique(
      test_points_inverse.begin(),
      test_points_inverse.end(),
      []( real_type a, real_type b ) { return abs( a - b ) < 1e-12; } ),
    test_points_inverse.end() );

  // Test inverse evaluation for each spline
  for ( integer spline_idx = 0; spline_idx < nspl; ++spline_idx )
  {
    string spline_name = string( ss.header( spline_idx ) );

    if ( ss.is_monotone( spline_idx ) <= 0 )
    {
      fmt::print( fg( fmt::color::orange ), "\nInverse evaluation skip {} because is not monotone:\n", spline_name );
      continue;
    }

    fmt::print( fg( fmt::color::green ), "\nInverse evaluation for {}:\n", spline_name );

    for ( real_type y{ ymin }; y <= ymax; y += ( ymax - ymin ) / 1000 )
    {
      real_type x;
      fileR << y;
      fileR_D << y;
      ss.eval2( spline_idx, y, x, val );
      ss.eval_D( x, val_D );
      for ( integer i = 0; i < nspl; ++i )
      {
        fileR << '\t' << val[i];
        fileR_D << '\t' << val_D[i];
      }
      fileR << '\n';
      fileR_D << '\n';
    }

    // Check inverse evaluation for a few sample points
    fmt::print( "  Sample inverse evaluations:\n" );
    for ( size_t i = 0; i < test_points_inverse.size(); ++i )
    {
      real_type x, y = test_points_inverse[i];
      ss.eval2( spline_idx, y, x, val );
      fmt::print( "    y = {:.3f} -> x = {:.3f} -> y(x)-y = {:3f} \n", y, x, val[spline_idx] - y );
    }
  }
  fileR.close();
  fileR_D.close();

  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║                    ALL TESTS COMPLETED                   ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  return 0;
}
