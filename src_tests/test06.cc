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

#include "Splines.hh"
#include "Utils_fmt.hh"

#include <GenericContainer/GenericContainer.hh>
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
static real_type yy[] = { 0, 1, 1.99, 2.0, 2.1 };

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

// Check derivative for constant spline (should be 0 in interior)
void check_constant_spline_derivative(
  SplineSet const & ss,
  real_type         x,
  integer           spline_idx,
  real_type const * knots,
  integer           npts,
  real_type &       D_spline,
  real_type &       D_fd,
  string &          note )
{
  // For constant spline, derivative should be 0 everywhere except at knots (undefined)
  D_spline = ss.D( x, spline_idx );

  // Check if x is exactly at a knot
  bool is_knot = false;
  for ( integer i = 0; i < npts; ++i )
  {
    if ( abs( knots[i] - x ) < 1e-12 )
    {
      is_knot = true;
      break;
    }
  }

  if ( is_knot )
  {
    // At knots, derivative of constant spline is undefined
    note = "knot (undefined)";
  }
  else
  {
    // In interior, derivative should be 0
    D_fd = 0.0;  // Expected derivative
    note = "interior (should be 0)";
  }
}

int main()
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║                         TEST N.6                         ║\n"
    "║              (GenericContainer Spline Test)              ║\n"
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
  real_type val[7], val_D[7];

  char const * headers[] = { "constant", "linear", "cubic", "akima", "vanleer", "pchip", "quintic_pchip" };

  real_type const * Y[] = { yy, yy, yy, yy, yy, yy, yy, yy };

  GC::GenericContainer gc;

  GC::vec_string_type & t    = gc["spline_type"].set_vec_string();
  GC::vec_string_type & hvec = gc["headers"].set_vec_string();  // Changed name from h to hvec
  t.resize( static_cast<size_t>( nspl ) );
  hvec.resize( static_cast<size_t>( nspl ) );
  std::copy_n( headers, nspl, hvec.begin() );
  std::copy_n( headers, nspl, t.begin() );

  GC::vector_type & data = gc["ydata"].set_vector();
  data.resize( static_cast<size_t>( nspl ) );
  for ( integer i = 0; i < nspl; ++i )
  {
    GC::GenericContainer & di = data[i];
    GC::vec_real_type &    v  = di.set_vec_real();
    if ( i == 0 )
    {
      // constant spline has 1 less point
      v.resize( static_cast<size_t>( npts - 1 ) );
      std::copy_n( Y[i], npts - 1, v.begin() );
    }
    else
    {
      v.resize( static_cast<size_t>( npts ) );
      std::copy_n( Y[i], npts, v.begin() );
    }
  }

  GC::vec_real_type & xdata = gc["xdata"].set_vec_real();
  xdata.resize( static_cast<size_t>( npts ) );
  std::copy_n( xx, npts, xdata.begin() );

  fmt::print( fg( fmt::color::yellow ), "\nGenericContainer data structure:\n" );
  gc.print( cout );

  ss.build( gc );  // nspl, npts, headers, stype, xx, Y, Yp );

  fmt::print( fg( fmt::color::green ), "\nSplineSet information:\n" );
  ss.info( cout );

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
  real_type const h_fd                = 1e-6;  // Changed name from h to h_fd
  real_type const tol                 = 1e-4;
  integer         total_points_tested = 0;
  integer         points_with_error   = 0;

  for ( integer spline_idx = 0; spline_idx < nspl; ++spline_idx )
  {
    fmt::print(
      "├{}┼{}┼{}┼{}┼{}┼{}┤\n",
      fmt::format( "{:─^{}}", "", 12 ),
      fmt::format( "{:─^{}}", "", 16 ),
      fmt::format( "{:─^{}}", "", 16 ),
      fmt::format( "{:─^{}}", "", 16 ),
      fmt::format( "{:─^{}}", "", 16 ),
      fmt::format( "{:─^{}}", "", 25 ) );

    string spline_name        = string( ss.header( spline_idx ) );
    bool   is_constant_spline = ( spline_idx == 0 );  // constant spline is at index 0

    for ( size_t pt_idx = 0; pt_idx < test_points.size(); ++pt_idx )
    {
      real_type x = test_points[pt_idx];

      bool      should_print = false;
      real_type D_spline, D_fd, abs_err, rel_err;
      string    note;

      // Special handling for constant spline
      if ( is_constant_spline )
      {
        check_constant_spline_derivative( ss, x, spline_idx, xx, npts, D_spline, D_fd, note );

        // Calculate relative error for constant spline
        abs_err = abs( D_spline - D_fd );
        rel_err = 100.0 * abs_err / max( { abs( D_fd ), abs( D_spline ), 1e-5 } );

        ++total_points_tested;
        if ( rel_err > tol ) ++points_with_error;

        // Print only if error or for sampling
        should_print = ( rel_err > tol || pt_idx % 4 == 0 );
      }
      else
      {
        // For non-constant splines
        // Determine if it's a knot
        bool is_knot = ( find( xx, xx + npts, x ) != xx + npts );

        // Calculate spline derivative
        D_spline = ss.D( x, spline_idx );

        // Calculate appropriate finite difference
        note = is_knot ? "knot" : "interior";

        // For knots, use appropriate one-sided finite differences
        if ( is_knot )
        {
          if ( x <= xmin + h_fd )
          {
            D_fd = finite_diff_forward( ss, x, spline_idx, h_fd );
            note += " (left boundary)";
          }
          else if ( x >= xmax - h_fd )
          {
            D_fd = finite_diff_backward( ss, x, spline_idx, h_fd );
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
              D_fd = finite_diff_forward( ss, x, spline_idx, h_fd );
              note += " (interval start)";
            }
            else
            {
              D_fd = finite_diff_backward( ss, x, spline_idx, h_fd );
              note += " (interval end)";
            }
          }
        }
        else
        {
          // For interior points, use central finite difference
          D_fd = finite_diff_central( ss, x, spline_idx, h_fd );
        }

        // Calculate relative error
        abs_err = abs( D_spline - D_fd );
        rel_err = 100.0 * abs_err / max( { abs( D_spline ), abs( D_fd ), 1e-5 } );

        ++total_points_tested;
        if ( rel_err > tol ) ++points_with_error;

        // Print only if: knot, significant error, or for reduced sampling
        should_print = ( is_knot || rel_err > tol || pt_idx % 4 == 0 );
      }

      if ( should_print )
      {
        vector<string> row_values = { fmt::format( "{:.6f}", x ),        spline_name,
                                      fmt::format( "{:.6f}", D_spline ), fmt::format( "{:.6f}", D_fd ),
                                      fmt::format( "{:.2e}", rel_err ),  note };
        print_table_row(
          row_values,
          rel_err < tol ? fmt::color::green : ( rel_err < 10.0 * tol ? fmt::color::yellow : fmt::color::red ) );
      }
    }
  }

  print_table_footer( table_headers.size() );

  // Statistics
  fmt::print(
    fg( fmt::color::blue ),
    "\n📊 Statistics for derivative check:\n"
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

  // Test inverse evaluation for quintic spline (index 5)
  integer quintic_spline_idx = 5;
  fmt::print( fg( fmt::color::green ), "\nInverse evaluation for {} spline:\n", ss.header( quintic_spline_idx ) );

  // Generate test points for inverse evaluation
  vector<real_type> test_points_inverse;

  // Use a range that makes sense for inverse evaluation
  real_type inv_ymin = yy[0];
  real_type inv_ymax = yy[npt - 1];

  for ( real_type y = inv_ymin; y <= inv_ymax; y += ( inv_ymax - inv_ymin ) / 20 )
  {
    test_points_inverse.push_back( y );
  }

  fmt::print( "  Sample inverse evaluations (y -> x):\n" );
  for ( real_type y : test_points_inverse )
  {
    real_type x;
    real_type tmp_val[7];
    ss.eval2( quintic_spline_idx, y, x, tmp_val );
    fmt::print( "    y = {:.3f} -> x = {:.3f}\n", y, x );
  }

  for ( real_type y = inv_ymin; y <= inv_ymax; y += ( inv_ymax - inv_ymin ) / 1000 )
  {
    real_type x;
    fileR << y;
    fileR_D << y;
    ss.eval2( quintic_spline_idx, y, x, val );
    ss.eval2_D( quintic_spline_idx, y, x, val_D );
    for ( integer i = 0; i < nspl; ++i )
    {
      fileR << '\t' << val[i];
      fileR_D << '\t' << val_D[i];
    }
    fileR << '\n';
    fileR_D << '\n';
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
