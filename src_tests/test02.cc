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
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <sstream>
#include <random>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::integer;
using Splines::real_type;
using Splines::Spline_sub_type;

// static real_type x_data1[] = { 0, 1, 2, 3, 4, 5 };
// static real_type y_data1[] = { 0, 1, 2, 3, 4, 5 };
// static real_type z_data1[] = { 10, 10, 10, 10, 10, 11, 10, 10, 10.5, 11, 9, 9, 11, 12, 13, 8, 7, 9,
//                               11, 12, 13, 8,  7,  9,  11, 12, 13,   8,  7, 9, 11, 12, 13, 8, 7, 9 };

// Enum per le diverse regioni di test
enum TestRegion
{
  INTERIOR_PATCH,         // Dentro le patch, lontano dai nodi
  X_KNOTS_Y_INTERIOR,     // Su nodi x, tra due nodi y
  Y_KNOTS_X_INTERIOR,     // Su nodi y, tra due nodi x
  X_Y_KNOTS_INTERSECTION  // Intersezioni nodi x,y
};

// Structure to hold 2D spline information
struct Spline2DInfo
{
  string    name;
  real_type x_min;
  real_type x_max;
  real_type y_min;
  real_type y_max;
  real_type z_min;
  real_type z_max;
  real_type z_avg;
};

// Structure to hold derivative errors statistics for different regions
struct DerivativeErrorStats
{
  real_type max_abs_error;
  real_type avg_abs_error;
  size_t    points_checked;

  DerivativeErrorStats() : max_abs_error( 0 ), avg_abs_error( 0 ), points_checked( 0 ) {}

  void update( real_type error )
  {
    max_abs_error = max( max_abs_error, error );
    avg_abs_error += error;
    points_checked++;
  }

  void finalize()
  {
    if ( points_checked > 0 ) { avg_abs_error /= static_cast<real_type>( points_checked ); }
  }
};

// Structure to hold all derivative errors for all regions
struct AllDerivativeErrors
{
  DerivativeErrorStats interior_Dx;
  DerivativeErrorStats interior_Dy;
  DerivativeErrorStats interior_Dxx;
  DerivativeErrorStats interior_Dyy;
  DerivativeErrorStats interior_Dxy;

  DerivativeErrorStats x_knots_Dx;
  DerivativeErrorStats x_knots_Dy;
  DerivativeErrorStats x_knots_Dxx;
  DerivativeErrorStats x_knots_Dyy;
  DerivativeErrorStats x_knots_Dxy;

  DerivativeErrorStats y_knots_Dx;
  DerivativeErrorStats y_knots_Dy;
  DerivativeErrorStats y_knots_Dxx;
  DerivativeErrorStats y_knots_Dyy;
  DerivativeErrorStats y_knots_Dxy;

  DerivativeErrorStats knots_intersection_Dx;
  DerivativeErrorStats knots_intersection_Dy;
  DerivativeErrorStats knots_intersection_Dxx;
  DerivativeErrorStats knots_intersection_Dyy;
  DerivativeErrorStats knots_intersection_Dxy;
};

// Helper function to format array as string
string format_array( const real_type * arr, size_t n )
{
  stringstream ss;
  for ( size_t i = 0; i < n; ++i )
  {
    if ( i > 0 ) ss << ", ";
    ss << fixed << setprecision( 2 ) << arr[i];
  }
  return ss.str();
}

// Function to print colored header
void print_header( const string & title )
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════╗\n"
    "║{:^46}║\n"
    "╚══════════════════════════════════════════════╝\n",
    title );
}

// Function to print a table of results for 2D splines
void print_2D_spline_table( const vector<Spline2DInfo> & results )
{
  // Table header
  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "\n"
    "┌──────────────────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐\n"
    "│ {:^20} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │\n"
    "├──────────────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤\n",
    "Spline Type",
    "x_min",
    "x_max",
    "y_min",
    "y_max",
    "z_min",
    "z_max",
    "z_avg" );

  // Table rows with alternating colors
  for ( size_t i = 0; i < results.size(); ++i )
  {
    const auto & r         = results[i];
    auto         row_color = ( i % 2 == 0 ) ? fg( fmt::color::light_green ) : fg( fmt::color::light_blue );

    fmt::print( row_color, "│ {:<20} ", r.name );
    fmt::print( row_color, "│ {:>8.3f} ", r.x_min );
    fmt::print( row_color, "│ {:>8.3f} ", r.x_max );
    fmt::print( row_color, "│ {:>8.3f} ", r.y_min );
    fmt::print( row_color, "│ {:>8.3f} ", r.y_max );

    auto z_min_color = ( r.z_min < 0 ) ? fg( fmt::color::red ) : fg( fmt::color::green );
    auto z_max_color = ( r.z_max < 0 ) ? fg( fmt::color::red ) : fg( fmt::color::green );
    auto z_avg_color = fg( fmt::color::orange );

    fmt::print( z_min_color, "│ {:>8.3f} ", r.z_min );
    fmt::print( z_max_color, "│ {:>8.3f} ", r.z_max );
    fmt::print( z_avg_color, "│ {:>8.3f} │\n", r.z_avg );
  }

  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "└──────────────────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘\n" );
}

// Function to format error with color based on tolerance
string format_error_cell( real_type max_err, real_type avg_err, size_t points, real_type tol_max, real_type tol_avg )
{
  stringstream ss;
  ss << fmt::format( "max:{:>7.1e} avg:{:>7.1e} ({:>2})", max_err, avg_err, points );

  // Determine color based on errors and tolerances
  if ( max_err < tol_max && avg_err < tol_avg )
  {
    return fmt::format( fg( fmt::color::green ) | fmt::emphasis::bold, "{:<30}", ss.str() );
  }
  else if ( max_err < tol_max * 10 && avg_err < tol_avg * 10 )
  {
    return fmt::format( fg( fmt::color::yellow ), "{:<30}", ss.str() );
  }
  else
  {
    return fmt::format( fg( fmt::color::red ) | fmt::emphasis::bold, "{:<30}", ss.str() );
  }
}

// Function to print detailed derivative table for all regions
void print_detailed_derivative_table(
  const string &                                    derivative_name,
  const vector<pair<string, AllDerivativeErrors>> & errors_by_spline,
  real_type                                         tol_max_interior,
  real_type                                         tol_avg_interior,
  real_type                                         tol_max_knots,
  real_type                                         tol_avg_knots,
  real_type                                         tol_max_intersection,
  real_type                                         tol_avg_intersection )
{
  // Table header - stringhe più compatte
  fmt::print(
    fg( fmt::color::light_blue ) | fmt::emphasis::bold,
    "\n"
    "┌─{:─^152}─┐\n"
    "│ {:^152} │\n"
    "├─{:─^20}─┬─{:─^30}─┬─{:─^30}─┬─{:─^30}─┬─{:─^30}─┤\n"
    "│ {:^20} │ {:^30} │ {:^30} │ {:^30} │ {:^30} │\n"
    "├─{:─^20}─┼─{:─^30}─┼─{:─^30}─┼─{:─^30}─┼─{:─^30}─┤\n",
    "",
    fmt::format( "DERIVATIVE {}", derivative_name ),
    "",
    "",
    "",
    "",
    "",
    "Spline",
    "Interior",
    "X-knots",
    "Y-knots",
    "Intersection",
    "",
    "",
    "",
    "",
    "" );

  // Function to get the appropriate error stats for the derivative and region
  auto get_error_stats =
    [&]( const AllDerivativeErrors & err, const string & deriv, TestRegion region ) -> const DerivativeErrorStats &
  {
    if ( deriv == "Dx" )
    {
      switch ( region )
      {
        case INTERIOR_PATCH: return err.interior_Dx;
        case X_KNOTS_Y_INTERIOR: return err.x_knots_Dx;
        case Y_KNOTS_X_INTERIOR: return err.y_knots_Dx;
        case X_Y_KNOTS_INTERSECTION: return err.knots_intersection_Dx;
      }
    }
    else if ( deriv == "Dy" )
    {
      switch ( region )
      {
        case INTERIOR_PATCH: return err.interior_Dy;
        case X_KNOTS_Y_INTERIOR: return err.x_knots_Dy;
        case Y_KNOTS_X_INTERIOR: return err.y_knots_Dy;
        case X_Y_KNOTS_INTERSECTION: return err.knots_intersection_Dy;
      }
    }
    else if ( deriv == "Dxx" )
    {
      switch ( region )
      {
        case INTERIOR_PATCH: return err.interior_Dxx;
        case X_KNOTS_Y_INTERIOR: return err.x_knots_Dxx;
        case Y_KNOTS_X_INTERIOR: return err.y_knots_Dxx;
        case X_Y_KNOTS_INTERSECTION: return err.knots_intersection_Dxx;
      }
    }
    else if ( deriv == "Dyy" )
    {
      switch ( region )
      {
        case INTERIOR_PATCH: return err.interior_Dyy;
        case X_KNOTS_Y_INTERIOR: return err.x_knots_Dyy;
        case Y_KNOTS_X_INTERIOR: return err.y_knots_Dyy;
        case X_Y_KNOTS_INTERSECTION: return err.knots_intersection_Dyy;
      }
    }
    else if ( deriv == "Dxy" )
    {
      switch ( region )
      {
        case INTERIOR_PATCH: return err.interior_Dxy;
        case X_KNOTS_Y_INTERIOR: return err.x_knots_Dxy;
        case Y_KNOTS_X_INTERIOR: return err.y_knots_Dxy;
        case X_Y_KNOTS_INTERSECTION: return err.knots_intersection_Dxy;
      }
    }
    static DerivativeErrorStats empty{};
    return empty;
  };

  // Table rows
  for ( size_t i = 0; i < errors_by_spline.size(); ++i )
  {
    const auto & [name, err] = errors_by_spline[i];
    auto row_color           = ( i % 2 == 0 ) ? fg( fmt::color::light_green ) : fg( fmt::color::light_blue );

    // Get stats for each region
    const auto & interior_stats     = get_error_stats( err, derivative_name, INTERIOR_PATCH );
    const auto & x_knots_stats      = get_error_stats( err, derivative_name, X_KNOTS_Y_INTERIOR );
    const auto & y_knots_stats      = get_error_stats( err, derivative_name, Y_KNOTS_X_INTERIOR );
    const auto & intersection_stats = get_error_stats( err, derivative_name, X_Y_KNOTS_INTERSECTION );

    // Format cells for each region
    string interior_cell = format_error_cell(
      interior_stats.max_abs_error,
      interior_stats.avg_abs_error,
      interior_stats.points_checked,
      tol_max_interior,
      tol_avg_interior );

    string x_knots_cell = format_error_cell(
      x_knots_stats.max_abs_error,
      x_knots_stats.avg_abs_error,
      x_knots_stats.points_checked,
      tol_max_knots,
      tol_avg_knots );

    string y_knots_cell = format_error_cell(
      y_knots_stats.max_abs_error,
      y_knots_stats.avg_abs_error,
      y_knots_stats.points_checked,
      tol_max_knots,
      tol_avg_knots );

    string intersection_cell = format_error_cell(
      intersection_stats.max_abs_error,
      intersection_stats.avg_abs_error,
      intersection_stats.points_checked,
      tol_max_intersection,
      tol_avg_intersection );

    // Print row
    fmt::print( row_color, "│ {:<20} ", name );
    fmt::print( "│ {} ", interior_cell );
    fmt::print( "│ {} ", x_knots_cell );
    fmt::print( "│ {} ", y_knots_cell );
    fmt::print( "│ {} │\n", intersection_cell );

    if ( i < errors_by_spline.size() - 1 )
    {
      fmt::print(
        fg( fmt::color::light_blue ),
        "├─{:─^20}─┼─{:─^30}─┼─{:─^30}─┼─{:─^30}─┼─{:─^30}─┤\n",
        "",
        "",
        "",
        "",
        "" );
    }
  }

  fmt::print(
    fg( fmt::color::light_blue ) | fmt::emphasis::bold,
    "└─{:─^20}─┴─{:─^30}─┴─{:─^30}─┴─{:─^30}─┴─{:─^30}─┘\n",
    "",
    "",
    "",
    "",
    "" );
}

// Finite difference approximations
namespace FiniteDifferences
{
  using real_type = double;

  namespace detail
  {
    constexpr real_type eps = std::numeric_limits<real_type>::epsilon();

    inline real_type h_first( real_type x, real_type y )
    {
      real_type scale = std::max( { real_type( 1 ), std::abs( x ), std::abs( y ) } );
      return std::cbrt( eps ) * scale;
    }

    // inline real_type h_second( real_type x, real_type y )
    //{
    //   real_type scale = std::max( { real_type( 1 ), std::abs( x ), std::abs( y ) } );
    //   return std::sqrt( std::sqrt( eps ) ) * scale;
    // }
    inline real_type h_second( real_type x, real_type y )
    {
      real_type scale = std::max( { real_type( 1 ), std::abs( x ), std::abs( y ) } );
      // Usa eps^(1/3) invece di sqrt(sqrt(eps)) per bilanciare meglio
      return std::cbrt( eps ) * scale;  // ~6e-6 per double
    }

  }  // namespace detail

  template <typename SplineType> real_type dx( const SplineType & spline, real_type x, real_type y )
  {
    real_type h = detail::h_first( x, y );
    return ( -spline( x + 2 * h, y ) + 8 * spline( x + h, y ) - 8 * spline( x - h, y ) + spline( x - 2 * h, y ) ) /
           ( 12 * h );
  }

  template <typename SplineType> real_type dy( const SplineType & spline, real_type x, real_type y )
  {
    real_type h = detail::h_first( x, y );
    return ( -spline( x, y + 2 * h ) + 8 * spline( x, y + h ) - 8 * spline( x, y - h ) + spline( x, y - 2 * h ) ) /
           ( 12 * h );
  }
  template <typename SplineType> real_type dxx( const SplineType & spline, real_type x, real_type y )
  {
    real_type h = detail::h_second( x, y );
    // Formula a 5 punti, O(h^4)
    return ( -spline( x + 2 * h, y ) + 16 * spline( x + h, y ) - 30 * spline( x, y ) + 16 * spline( x - h, y ) -
             spline( x - 2 * h, y ) ) /
           ( 12 * h * h );
  }

  template <typename SplineType> real_type dyy( const SplineType & spline, real_type x, real_type y )
  {
    real_type h = detail::h_second( x, y );
    // Formula a 5 punti, O(h^4)
    return ( -spline( x, y + 2 * h ) + 16 * spline( x, y + h ) - 30 * spline( x, y ) + 16 * spline( x, y - h ) -
             spline( x, y - 2 * h ) ) /
           ( 12 * h * h );
  }

  template <typename SplineType> real_type dxy( const SplineType & spline, real_type x, real_type y )
  {
    real_type h = detail::h_second( x, y );
    return ( spline( x + h, y + h ) - spline( x + h, y - h ) - spline( x - h, y + h ) + spline( x - h, y - h ) ) /
           ( 4 * h * h );
  }
}  // namespace FiniteDifferences

// Function to check derivatives in all four regions
template <typename SplineType> AllDerivativeErrors check_derivatives_by_region(
  SplineType &      spline,
  const string &    spline_name,
  const real_type * X_knots,
  integer           nx_knots,
  const real_type * Y_knots,
  integer           ny_knots,
  integer           points_per_region = 20 )
{
  AllDerivativeErrors errors;

  // Setup random number generator (C++17)
  std::random_device                        rd;
  std::mt19937                              gen( rd() );
  std::uniform_real_distribution<real_type> dis( 0.2, 0.8 );

  // Generate test points for each region
  vector<pair<real_type, real_type>> test_points_interior;
  vector<pair<real_type, real_type>> test_points_x_knots;
  vector<pair<real_type, real_type>> test_points_y_knots;
  vector<pair<real_type, real_type>> test_points_intersection;

  // 1. Points in interior of patches (away from knots)
  integer patches_x = nx_knots - 1;
  integer patches_y = ny_knots - 1;

  for ( integer patch_i = 0; patch_i < patches_x; ++patch_i )
  {
    for ( integer patch_j = 0; patch_j < patches_y; ++patch_j )
    {
      real_type patch_x_min = X_knots[patch_i];
      real_type patch_x_max = X_knots[patch_i + 1];
      real_type patch_y_min = Y_knots[patch_j];
      real_type patch_y_max = Y_knots[patch_j + 1];

      // Generate points in the middle of the patch
      for ( integer k = 0; k < points_per_region / ( patches_x * patches_y ); ++k )
      {
        real_type x = patch_x_min + ( patch_x_max - patch_x_min ) * dis( gen );
        real_type y = patch_y_min + ( patch_y_max - patch_y_min ) * dis( gen );

        test_points_interior.emplace_back( x, y );
      }
    }
  }

  // 2. Points on x-knots but between y-knots
  for ( integer i = 0; i < nx_knots; ++i )
  {
    real_type x = X_knots[i];
    for ( integer patch_j = 0; patch_j < patches_y; ++patch_j )
    {
      real_type patch_y_min = Y_knots[patch_j];
      real_type patch_y_max = Y_knots[patch_j + 1];

      for ( integer k = 0; k < points_per_region / nx_knots; ++k )
      {
        real_type y = patch_y_min + ( patch_y_max - patch_y_min ) * dis( gen );
        test_points_x_knots.emplace_back( x, y );
      }
    }
  }

  // 3. Points on y-knots but between x-knots
  for ( integer j = 0; j < ny_knots; ++j )
  {
    real_type y = Y_knots[j];
    for ( integer patch_i = 0; patch_i < patches_x; ++patch_i )
    {
      real_type patch_x_min = X_knots[patch_i];
      real_type patch_x_max = X_knots[patch_i + 1];

      for ( integer k = 0; k < points_per_region / ny_knots; ++k )
      {
        real_type x = patch_x_min + ( patch_x_max - patch_x_min ) * dis( gen );
        test_points_y_knots.emplace_back( x, y );
      }
    }
  }

  // 4. Points at intersection of x and y knots (grid points)
  for ( integer i = 0; i < nx_knots; ++i )
  {
    for ( integer j = 0; j < ny_knots; ++j ) { test_points_intersection.emplace_back( X_knots[i], Y_knots[j] ); }
  }

  // Test points in each region
  auto test_region = [&]( const vector<pair<real_type, real_type>> & points, TestRegion region )
  {
    for ( const auto & [x_val, y_val] : points )
    {
      // Get analytical derivatives
      real_type Dx_analytic  = spline.Dx( x_val, y_val );
      real_type Dy_analytic  = spline.Dy( x_val, y_val );
      real_type Dxx_analytic = spline.Dxx( x_val, y_val );
      real_type Dyy_analytic = spline.Dyy( x_val, y_val );
      real_type Dxy_analytic = spline.Dxy( x_val, y_val );

      // Get finite difference approximations
      real_type Dx_fd  = FiniteDifferences::dx( spline, x_val, y_val );
      real_type Dy_fd  = FiniteDifferences::dy( spline, x_val, y_val );
      real_type Dxx_fd = FiniteDifferences::dxx( spline, x_val, y_val );
      real_type Dyy_fd = FiniteDifferences::dyy( spline, x_val, y_val );
      real_type Dxy_fd = FiniteDifferences::dxy( spline, x_val, y_val );

      // Calculate absolute errors
      real_type abs_error_Dx  = fabs( Dx_analytic - Dx_fd );
      real_type abs_error_Dy  = fabs( Dy_analytic - Dy_fd );
      real_type abs_error_Dxx = fabs( Dxx_analytic - Dxx_fd );
      real_type abs_error_Dyy = fabs( Dyy_analytic - Dyy_fd );
      real_type abs_error_Dxy = fabs( Dxy_analytic - Dxy_fd );

      // Update appropriate statistics based on region
      switch ( region )
      {
        case INTERIOR_PATCH:
          errors.interior_Dx.update( abs_error_Dx );
          errors.interior_Dy.update( abs_error_Dy );
          errors.interior_Dxx.update( abs_error_Dxx );
          errors.interior_Dyy.update( abs_error_Dyy );
          errors.interior_Dxy.update( abs_error_Dxy );
          break;

        case X_KNOTS_Y_INTERIOR:
          errors.x_knots_Dx.update( abs_error_Dx );
          errors.x_knots_Dy.update( abs_error_Dy );
          errors.x_knots_Dxx.update( abs_error_Dxx );
          errors.x_knots_Dyy.update( abs_error_Dyy );
          errors.x_knots_Dxy.update( abs_error_Dxy );
          break;

        case Y_KNOTS_X_INTERIOR:
          errors.y_knots_Dx.update( abs_error_Dx );
          errors.y_knots_Dy.update( abs_error_Dy );
          errors.y_knots_Dxx.update( abs_error_Dxx );
          errors.y_knots_Dyy.update( abs_error_Dyy );
          errors.y_knots_Dxy.update( abs_error_Dxy );
          break;

        case X_Y_KNOTS_INTERSECTION:
          errors.knots_intersection_Dx.update( abs_error_Dx );
          errors.knots_intersection_Dy.update( abs_error_Dy );
          errors.knots_intersection_Dxx.update( abs_error_Dxx );
          errors.knots_intersection_Dyy.update( abs_error_Dyy );
          errors.knots_intersection_Dxy.update( abs_error_Dxy );
          break;
      }
    }
  };

  fmt::print( fg( fmt::color::blue ), "   Checking derivatives for {} by region...\n", spline_name );

  test_region( test_points_interior, INTERIOR_PATCH );
  test_region( test_points_x_knots, X_KNOTS_Y_INTERIOR );
  test_region( test_points_y_knots, Y_KNOTS_X_INTERIOR );
  test_region( test_points_intersection, X_Y_KNOTS_INTERSECTION );

  // Finalize all statistics (calculate averages)
  errors.interior_Dx.finalize();
  errors.interior_Dy.finalize();
  errors.interior_Dxx.finalize();
  errors.interior_Dyy.finalize();
  errors.interior_Dxy.finalize();

  errors.x_knots_Dx.finalize();
  errors.x_knots_Dy.finalize();
  errors.x_knots_Dxx.finalize();
  errors.x_knots_Dyy.finalize();
  errors.x_knots_Dxy.finalize();

  errors.y_knots_Dx.finalize();
  errors.y_knots_Dy.finalize();
  errors.y_knots_Dxx.finalize();
  errors.y_knots_Dyy.finalize();
  errors.y_knots_Dxy.finalize();

  errors.knots_intersection_Dx.finalize();
  errors.knots_intersection_Dy.finalize();
  errors.knots_intersection_Dxx.finalize();
  errors.knots_intersection_Dyy.finalize();
  errors.knots_intersection_Dxy.finalize();

  fmt::print(
    fg( fmt::color::green ),
    "   ✓ Tested: {} interior, {} x-knots, {} y-knots, {} intersection points\n",
    test_points_interior.size(),
    test_points_x_knots.size(),
    test_points_y_knots.size(),
    test_points_intersection.size() );

  return errors;
}

// Function to generate grid data and save to files
template <typename SplineType> void generate_grid_data(
  SplineType &   spline,
  const string & base_filename,
  const string & spline_name,
  integer        grid_size = 100 )
{
  // Array of filenames
  array<string, 6> filenames = {
    fmt::format( "out/{}.txt", base_filename ),     fmt::format( "out/{}_Dx.txt", base_filename ),
    fmt::format( "out/{}_Dy.txt", base_filename ),  fmt::format( "out/{}_Dxx.txt", base_filename ),
    fmt::format( "out/{}_Dyy.txt", base_filename ), fmt::format( "out/{}_Dxy.txt", base_filename )
  };

  // Array of output streams
  array<ofstream, 6> files;

  // Open all files
  for ( size_t idx = 0; idx < filenames.size(); ++idx )
  {
    files[idx].open( filenames[idx].data() );
    if ( !files[idx].is_open() )
    {
      fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold, "❌ Error opening file: {}\n", filenames[idx] );
      return;
    }
  }

  real_type x_min = spline.x_min();
  real_type x_max = spline.x_max();
  real_type y_min = spline.y_min();
  real_type y_max = spline.y_max();

  real_type dx = ( x_max - x_min ) / grid_size;
  real_type dy = ( y_max - y_min ) / grid_size;

  fmt::print( fg( fmt::color::blue ), "   Generating {}x{} grid for {}...", grid_size + 1, grid_size + 1, spline_name );

  for ( integer i = 0; i <= grid_size; ++i )
  {
    real_type xxx = x_min + i * dx;
    for ( integer j = 0; j <= grid_size; ++j )
    {
      real_type yyy = y_min + j * dy;

      fmt::print( files[0], "{}\t", spline( xxx, yyy ) );
      fmt::print( files[1], "{}\t", spline.Dx( xxx, yyy ) );
      fmt::print( files[2], "{}\t", spline.Dy( xxx, yyy ) );
      fmt::print( files[3], "{}\t", spline.Dxx( xxx, yyy ) );
      fmt::print( files[4], "{}\t", spline.Dyy( xxx, yyy ) );
      fmt::print( files[5], "{}\t", spline.Dxy( xxx, yyy ) );
    }

    for ( auto & file : files ) { file << '\n'; }
  }

  for ( auto & file : files ) { file.close(); }

  fmt::print( fg( fmt::color::green ), " ✓\n" );
}

// Function to analyze a 2D spline and collect statistics
template <typename SplineType>
Spline2DInfo analyze_2D_spline( SplineType & spline, const string & name, integer sample_points = 50 )
{
  Spline2DInfo info;
  info.name  = name;
  info.x_min = spline.x_min();
  info.x_max = spline.x_max();
  info.y_min = spline.y_min();
  info.y_max = spline.y_max();

  real_type z_min = numeric_limits<real_type>::max();
  real_type z_max = numeric_limits<real_type>::lowest();
  real_type z_sum = 0;

  real_type dx = ( info.x_max - info.x_min ) / sample_points;
  real_type dy = ( info.y_max - info.y_min ) / sample_points;

  for ( integer i = 0; i <= sample_points; ++i )
  {
    real_type xxx = info.x_min + i * dx;
    for ( integer j = 0; j <= sample_points; ++j )
    {
      real_type yyy = info.y_min + j * dy;
      real_type zzz = spline( xxx, yyy );

      z_sum += zzz;
      z_min = min( z_min, zzz );
      z_max = max( z_max, zzz );
    }
  }

  info.z_min = z_min;
  info.z_max = z_max;
  info.z_avg = z_sum / ( ( sample_points + 1 ) * ( sample_points + 1 ) );

  return info;
}

int main()
{
  print_header( "2D SPLINE INTERPOLATION TEST SUITE" );

  static real_type X[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  static real_type Y[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  static real_type Z[] = {
    0,      0,      0,      0,      0,      0, 0, 0,      0,      -0.1000, 0,       0,      0,      0, 0, 0, 0, 0,
    0,      0,      0,      0,      0,      0, 0, 0,      0,      0.5000,  0.5000,  0,      0,      0, 0, 0, 0, 0,
    0,      0.5000, 0.5000, 0,      0,      0, 0, 0,      0,      0.2000,  0.2500,  0,      0,      0, 0, 0, 0, 0,
    0,      0.2500, 0.2000, 0,      0,      0, 0, 0.4286, 0.4286, 0.4286,  0.4286,  0,      0,      0, 0, 0, 0, 0.4286,
    0.4286, 0.4286, 0.4286, 0,      0,      0, 0, 0,      0,      0.4286,  0.4286,  0.4286, 0.4286, 0, 0, 0, 0, 0,
    0,      0.4286, 0.4286, 0.4286, 0.4286, 0, 0, 0,      0,      0,       -0.1000, 0,      0,      0, 0, 0, 0, 0,
    0,      0
  };

  integer NX = 10;
  integer NY = 11;
  integer LD = NX;
  bool    f  = true;  // fortran
  bool    t  = true;  // transpose

  fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n📊 Dataset Information:\n" );
  fmt::print( fg( fmt::color::gray ), "   Grid size: {}x{} points\n", NX, NY );

  string x_knots_str = format_array( X, NX );
  string y_knots_str = format_array( Y, NY );

  fmt::print( fg( fmt::color::gray ), "   x knots: {}\n", x_knots_str );
  fmt::print( fg( fmt::color::gray ), "   y knots: {}\n", y_knots_str );
  fmt::print( fg( fmt::color::gray ), "   Domain: x ∈ [{}, {}], y ∈ [{}, {}]\n", X[0], X[NX - 1], Y[0], Y[NY - 1] );

  BiQuinticSpline bq_cubic( Spline_sub_type::CUBIC );
  BiQuinticSpline bq_akima( Spline_sub_type::AKIMA );
  BiQuinticSpline bq_vanleer( Spline_sub_type::VANLEER );
  BiQuinticSpline bq_pchip( Spline_sub_type::PCHIP );
  BiCubicSpline   bc_cubic( Spline_sub_type::CUBIC );
  BiCubicSpline   bc_akima( Spline_sub_type::AKIMA );
  BiCubicSpline   bc_vanleer( Spline_sub_type::VANLEER );
  BiCubicSpline   bc_pchip( Spline_sub_type::PCHIP );
  BilinearSpline  bl;

  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔨 Building 2D splines...\n" );

  bq_cubic.build( X, 1, Y, 1, Z, LD, NX, NY, f, t );
  fmt::print( fg( fmt::color::green ), "   ✓ BiQuinticSpline built\n" );

  bq_akima.build( X, 1, Y, 1, Z, LD, NX, NY, f, t );
  fmt::print( fg( fmt::color::green ), "   ✓ BiQuinticSpline[akima] built\n" );

  bq_vanleer.build( X, 1, Y, 1, Z, LD, NX, NY, f, t );
  fmt::print( fg( fmt::color::green ), "   ✓ BiQuinticSpline[vanleer] built\n" );

  bq_pchip.build( X, 1, Y, 1, Z, LD, NX, NY, f, t );
  fmt::print( fg( fmt::color::green ), "   ✓ BiQuinticSpline[pchip] built\n" );

  bc_cubic.build( X, 1, Y, 1, Z, LD, NX, NY, f, t );
  fmt::print( fg( fmt::color::green ), "   ✓ BiCubicSpline built\n" );

  bc_akima.build( X, 1, Y, 1, Z, LD, NX, NY, f, t );
  fmt::print( fg( fmt::color::green ), "   ✓ BiCubicSpline[akima] built\n" );

  bc_vanleer.build( X, 1, Y, 1, Z, LD, NX, NY, f, t );
  fmt::print( fg( fmt::color::green ), "   ✓ BiCubicSpline[vanleer] built\n" );

  bc_pchip.build( X, 1, Y, 1, Z, LD, NX, NY, f, t );
  fmt::print( fg( fmt::color::green ), "   ✓ BiCubicSpline[pchip] built\n" );

  bl.build( X, 1, Y, 1, Z, LD, NX, NY, f, t );
  fmt::print( fg( fmt::color::green ), "   ✓ BilinearSpline built\n" );

  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n📈 Analyzing spline properties...\n" );

  vector<Spline2DInfo> spline_results;
  spline_results.push_back( analyze_2D_spline( bc_cubic, "BiCubic" ) );
  spline_results.push_back( analyze_2D_spline( bc_akima, "BiCubic[Akima]" ) );
  spline_results.push_back( analyze_2D_spline( bc_vanleer, "BiCubic[VanLeer]" ) );
  spline_results.push_back( analyze_2D_spline( bc_pchip, "BiCubic[Pchip]" ) );
  spline_results.push_back( analyze_2D_spline( bq_cubic, "BiQuintic" ) );
  spline_results.push_back( analyze_2D_spline( bq_akima, "BiQuintic[Akima]" ) );
  spline_results.push_back( analyze_2D_spline( bq_vanleer, "BiQuintic[VanLeer]" ) );
  spline_results.push_back( analyze_2D_spline( bq_pchip, "BiQuintic[Pchip]" ) );
  spline_results.push_back( analyze_2D_spline( bl, "Bilinear" ) );

  print_2D_spline_table( spline_results );

  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n🔍 Checking derivatives by region...\n" );

  vector<pair<string, AllDerivativeErrors>> errors_by_region;

  integer NPT = 30;

  errors_by_region.emplace_back(
    "BiCubic",
    check_derivatives_by_region( bc_cubic, "BiCubicSpline", X, NX, Y, NY, NPT ) );
  errors_by_region.emplace_back(
    "BiCubic[Akima]",
    check_derivatives_by_region( bc_akima, "BiCubicSpline[Akima]", X, NX, Y, NY, NPT ) );
  errors_by_region.emplace_back(
    "BiCubic[VanLeer]",
    check_derivatives_by_region( bc_vanleer, "BiCubicSpline[VanLeer]", X, NX, Y, NY, NPT ) );
  errors_by_region.emplace_back(
    "BiCubic[Pchip]",
    check_derivatives_by_region( bc_pchip, "BiCubicSpline[Pchip]", X, NX, Y, NY, NPT ) );

  errors_by_region.emplace_back(
    "BiQuintic",
    check_derivatives_by_region( bq_cubic, "BiQuinticSpline", X, NX, Y, NY, NPT ) );
  errors_by_region.emplace_back(
    "BiQuintic[Akima]",
    check_derivatives_by_region( bq_akima, "BiQuinticSpline[Akima]", X, NX, Y, NY, NPT ) );
  errors_by_region.emplace_back(
    "BiQuintic[VanLeer]",
    check_derivatives_by_region( bq_vanleer, "BiQuinticSpline[VanLeer]", X, NX, Y, NY, NPT ) );
  errors_by_region.emplace_back(
    "BiQuintic[Pchip]",
    check_derivatives_by_region( bq_pchip, "BiQuinticSpline[Pchip]", X, NX, Y, NY, NPT ) );

  errors_by_region.emplace_back( "Bilinear", check_derivatives_by_region( bl, "BilinearSpline", X, NX, Y, NY, NPT ) );

  fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n📊 DERIVATIVE ERROR ANALYSIS BY REGION:\n" );

  real_type first_deriv_tol_max_interior     = 1e-6;
  real_type first_deriv_tol_avg_interior     = 1e-7;
  real_type first_deriv_tol_max_knots        = 1e-5;
  real_type first_deriv_tol_avg_knots        = 1e-6;
  real_type first_deriv_tol_max_intersection = 1e-4;
  real_type first_deriv_tol_avg_intersection = 1e-5;

  real_type second_deriv_tol_max_interior     = 1e-3;
  real_type second_deriv_tol_avg_interior     = 1e-4;
  real_type second_deriv_tol_max_knots        = 1e-2;
  real_type second_deriv_tol_avg_knots        = 1e-3;
  real_type second_deriv_tol_max_intersection = 5e-2;
  real_type second_deriv_tol_avg_intersection = 1e-2;

  print_detailed_derivative_table(
    "Dx",
    errors_by_region,
    first_deriv_tol_max_interior,
    first_deriv_tol_avg_interior,
    first_deriv_tol_max_knots,
    first_deriv_tol_avg_knots,
    first_deriv_tol_max_intersection,
    first_deriv_tol_avg_intersection );

  print_detailed_derivative_table(
    "Dy",
    errors_by_region,
    first_deriv_tol_max_interior,
    first_deriv_tol_avg_interior,
    first_deriv_tol_max_knots,
    first_deriv_tol_avg_knots,
    first_deriv_tol_max_intersection,
    first_deriv_tol_avg_intersection );

  print_detailed_derivative_table(
    "Dxx",
    errors_by_region,
    second_deriv_tol_max_interior,
    second_deriv_tol_avg_interior,
    second_deriv_tol_max_knots,
    second_deriv_tol_avg_knots,
    second_deriv_tol_max_intersection,
    second_deriv_tol_avg_intersection );

  print_detailed_derivative_table(
    "Dyy",
    errors_by_region,
    second_deriv_tol_max_interior,
    second_deriv_tol_avg_interior,
    second_deriv_tol_max_knots,
    second_deriv_tol_avg_knots,
    second_deriv_tol_max_intersection,
    second_deriv_tol_avg_intersection );

  print_detailed_derivative_table(
    "Dxy",
    errors_by_region,
    second_deriv_tol_max_interior,
    second_deriv_tol_avg_interior,
    second_deriv_tol_max_knots,
    second_deriv_tol_avg_knots,
    second_deriv_tol_max_intersection,
    second_deriv_tol_avg_intersection );

  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n💾 Writing spline data to files...\n" );

  generate_grid_data( bc_cubic, "bicubic", "BiCubicSpline" );
  generate_grid_data( bc_akima, "bicubic[Akima]", "BiCubicSpline[Akima]" );
  generate_grid_data( bc_vanleer, "bicubic[VanLeer]", "BiCubicSpline[VanLeer]" );
  generate_grid_data( bc_pchip, "bicubic[Pchip]", "BiCubicSpline[Pchip]" );

  generate_grid_data( bq_cubic, "biquintic", "BiQuinticSpline" );
  generate_grid_data( bq_akima, "biquintic[Akima]", "BiQuinticSpline[Akima]" );
  generate_grid_data( bq_vanleer, "biquintic[VanLeer]", "BiQuinticSpline[VanLeer]" );
  generate_grid_data( bq_pchip, "biquintic[Pchip]", "BiQuinticSpline[Pchip]" );

  generate_grid_data( bl, "bilinear", "BilinearSpline" );

  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n📁 Generated Files Summary:\n" );

  vector<string> spline_types = { "bicubic", "biquintic", "akima2d", "bilinear" };
  vector<string> derivatives  = { "", "_Dx", "_Dy", "_Dxx", "_Dyy", "_Dxy" };

  for ( const auto & spline_type : spline_types )
  {
    fmt::print( fg( fmt::color::gray ), "   {}:\n", spline_type );
    for ( const auto & deriv : derivatives )
    {
      fmt::print( fg( fmt::color::gray ), "     - out/{}{}.txt\n", spline_type, deriv );
    }
  }

  print_header( "2D SPLINE TEST COMPLETED" );
  fmt::print( fg( fmt::color::light_green ) | fmt::emphasis::bold, "\n🎉 All 2D spline analyses completed! 🎉\n" );
  fmt::print( fg( fmt::color::gray ), "   Grid data saved to: out/*.txt\n" );
  fmt::print( fg( fmt::color::gray ), "   Derivative analysis completed for 4 regions:\n" );
  fmt::print( fg( fmt::color::gray ), "     1. Interior of patches\n" );
  fmt::print( fg( fmt::color::gray ), "     2. X-knots (between Y knots)\n" );
  fmt::print( fg( fmt::color::gray ), "     3. Y-knots (between X knots)\n" );
  fmt::print( fg( fmt::color::gray ), "     4. X-Y knots intersection\n\n" );

  return 0;
}
