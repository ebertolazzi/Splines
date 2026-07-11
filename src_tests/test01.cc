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
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"
#include "Utils_FD.hh"

using namespace SplinesLoad;
using namespace std;
using Splines::integer;
using Splines::real_type;
using Splines::Spline_sub_type;
using Utils::m_pi;

// ============================================================================
// TEST DATASETS
// ============================================================================

static real_type xx0[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
static real_type yy0[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

static real_type xx1[] = { 0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15 };
static real_type yy1[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

static real_type xx2[] = { 0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15 };
static real_type yy2[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };

// RPN 14
static real_type xx3[] = { 7.99, 8.09, 8.19, 8.7, 9.2, 10, 12, 15, 20 };
static real_type yy3[] = { 0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994 };

// Titanium
static real_type xx4[] = { 595, 635, 695, 795, 855, 875, 895, 915, 935, 985, 1035, 1075 };
static real_type yy4[] = { 0.644, 0.652, 0.644, 0.694, 0.907, 1.336, 2.169, 1.598, 0.916, 0.607, 0.603, 0.608 };

// toolpath
static real_type xx5[] = { 0.11, 0.12, 0.15, 0.16 };
static real_type yy5[] = { 0.0003, 0.0003, 0.0004, 0.0004 };

static integer nn[] = { 11, 11, 11, 9, 12, 4 };

// ============================================================================
// CONVERGENCE TEST FUNCTIONS AND STRUCTURES
// ============================================================================

// Definizione del tipo per le funzioni test
typedef real_type ( *TestFunctionPtr )( real_type );

// Test 1: Polinomio di grado 1: f(x) = 2x + 1
real_type test_poly1( real_type x )
{ return 2.0 * x + 1.0; }
real_type test_poly1_D( real_type )
{ return 2.0; }
real_type test_poly1_DD( real_type )
{ return 0.0; }

// Test 2: Polinomio di grado 2: f(x) = x^2 - 3x + 2
real_type test_poly2( real_type x )
{ return x * x - 3.0 * x + 2.0; }
real_type test_poly2_D( real_type x )
{ return 2.0 * x - 3.0; }
real_type test_poly2_DD( real_type )
{ return 2.0; }

// Test 3: Polinomio di grado 3: f(x) = x^3 - 2x^2 + x - 5
real_type test_poly3( real_type x )
{ return x * x * x - 2.0 * x * x + x - 5.0; }
real_type test_poly3_D( real_type x )
{ return 3.0 * x * x - 4.0 * x + 1.0; }
real_type test_poly3_DD( real_type x )
{ return 6.0 * x - 4.0; }

// Test 4: Polinomio di grado 4: f(x) = x^4 + x^3 - 2x^2 + 3x - 1
real_type test_poly4( real_type x )
{
  real_type x2 = x * x;
  return x2 * x2 + x * x2 - 2.0 * x2 + 3.0 * x - 1.0;
}
real_type test_poly4_D( real_type x )
{ return 4.0 * x * x * x + 3.0 * x * x - 4.0 * x + 3.0; }
real_type test_poly4_DD( real_type x )
{ return 12.0 * x * x + 6.0 * x - 4.0; }

// Test 5: Funzione esponenziale/trigonometrica: f(x) = exp(x) + sin(2x)
real_type test_exp_trig( real_type x )
{ return exp( x ) + sin( 2.0 * x ); }
real_type test_exp_trig_D( real_type x )
{ return exp( x ) + 2.0 * cos( 2.0 * x ); }
real_type test_exp_trig_DD( real_type x )
{ return exp( x ) - 4.0 * sin( 2.0 * x ); }

// Struttura per informazioni sulle funzioni test
struct TestFunctionInfo
{
  string          name;
  TestFunctionPtr func;
  TestFunctionPtr deriv;
  TestFunctionPtr deriv2;
  real_type       a;  // dominio minimo
  real_type       b;  // dominio massimo

  TestFunctionInfo(
    const string &  n,
    TestFunctionPtr f,
    TestFunctionPtr d,
    TestFunctionPtr d2,
    real_type       min_domain,
    real_type       max_domain )
    : name( n ), func( f ), deriv( d ), deriv2( d2 ), a( min_domain ), b( max_domain )
  {
  }
};

// Array di funzioni test
vector<TestFunctionInfo> test_functions = {
  { "Polynomial deg 1", test_poly1, test_poly1_D, test_poly1_DD, -2.0, 2.0 },
  { "Polynomial deg 2", test_poly2, test_poly2_D, test_poly2_DD, -2.0, 2.0 },
  { "Polynomial deg 3", test_poly3, test_poly3_D, test_poly3_DD, -2.0, 2.0 },
  { "Polynomial deg 4", test_poly4, test_poly4_D, test_poly4_DD, -2.0, 2.0 },
  { "Exp + sin(2x)", test_exp_trig, test_exp_trig_D, test_exp_trig_DD, 0.0, 2.0 }
};

// Generate uniform mesh
vector<real_type> generate_uniform_mesh( real_type a, real_type b, integer N )
{
  vector<real_type> mesh( N );
  real_type         s = ( b - a ) / ( N - 1 );
  for ( integer i = 0; i < N; ++i ) mesh[i] = a + i * s;
  return mesh;
}

// Generate non-uniform mesh (random with refinement)
vector<real_type> generate_nonuniform_mesh( real_type a, real_type b, integer N )
{
  vector<real_type> mesh( N );

  // Seed for random number generation
  static std::random_device rd;
  static std::mt19937       gen( rd() );

  if ( N <= 2 )
  {
    mesh[0] = a;
    if ( N > 1 ) mesh[N - 1] = b;
    return mesh;
  }

  // Start with endpoints
  mesh[0]     = a;
  mesh[N - 1] = b;

  // Generate random interior points
  std::uniform_real_distribution<real_type> dis( a, b );
  for ( integer i = 1; i < N - 1; ++i ) { mesh[i] = dis( gen ); }

  // Sort the points
  sort( mesh.begin(), mesh.end() );

  return mesh;
}

// ============================================================================
// FUNCTIONS FOR HIERARCHICAL MESH REFINEMENT (FOR CONVERGENCE TESTS)
// ============================================================================

// Generate a base non-uniform mesh
vector<real_type> generate_base_nonuniform_mesh( real_type a, real_type b, integer base_N )
{
  vector<real_type> mesh( base_N );

  if ( base_N <= 2 )
  {
    mesh[0] = a;
    if ( base_N > 1 ) mesh[base_N - 1] = b;
    return mesh;
  }

  // Start with endpoints
  mesh[0]          = a;
  mesh[base_N - 1] = b;

  // Use Chebyshev-like distribution for interior points
  // This creates a non-uniform mesh with more points near boundaries
  for ( integer i = 1; i < base_N - 1; ++i )
  {
    real_type theta = M_PI * i / ( base_N - 1 );
    mesh[i]         = a + 0.5 * ( b - a ) * ( 1 - cos( theta ) );
  }

  return mesh;
}

// Refine mesh by inserting midpoints of all intervals
vector<real_type> refine_mesh_by_insertion( const vector<real_type> & old_mesh )
{
  integer           old_N = static_cast<integer>( old_mesh.size() );
  integer           new_N = 2 * old_N - 1;  // Inserting midpoints
  vector<real_type> new_mesh( new_N );

  for ( integer i = 0; i < old_N - 1; ++i )
  {
    // Original point
    new_mesh[2 * i] = old_mesh[i];
    // Midpoint
    new_mesh[2 * i + 1] = 0.5 * ( old_mesh[i] + old_mesh[i + 1] );
  }
  // Last point
  new_mesh[new_N - 1] = old_mesh[old_N - 1];

  return new_mesh;
}

// Calculate characteristic mesh size (average interval length)
real_type compute_mesh_size( const vector<real_type> & X )
{
  integer   N            = static_cast<integer>( X.size() );
  real_type total_length = 0.0;

  for ( integer i = 0; i < N - 1; ++i ) { total_length += ( X[i + 1] - X[i] ); }

  return total_length / ( N - 1 );  // Average interval length
}

// Calculate maximum interval length
real_type compute_max_interval_length( const vector<real_type> & X )
{
  integer   N          = static_cast<integer>( X.size() );
  real_type max_length = 0.0;

  for ( integer i = 0; i < N - 1; ++i )
  {
    real_type length = X[i + 1] - X[i];
    if ( length > max_length ) max_length = length;
  }

  return max_length;
}

// Structure to hold convergence results
struct ConvergenceResult
{
  string            spline_name;
  string            test_func_name;
  vector<integer>   N_values;
  vector<real_type> errors;
  vector<real_type> orders;
  vector<real_type> mesh_sizes;  // Added to store actual mesh sizes

  ConvergenceResult( const string & spline, const string & test_func )
    : spline_name( spline ), test_func_name( test_func )
  {
  }
};

// ============================================================================
// STRUCTURES FOR DERIVATIVE ERROR ANALYSIS
// ============================================================================

// Structure to hold derivative error statistics
struct DerivativeErrorStats
{
  real_type max_abs_error_interior;
  real_type avg_abs_error_interior;
  size_t    points_checked_interior;

  real_type max_abs_error_boundary;
  real_type avg_abs_error_boundary;
  size_t    points_checked_boundary;

  DerivativeErrorStats()
    : max_abs_error_interior( 0 )
    , avg_abs_error_interior( 0 )
    , points_checked_interior( 0 )
    , max_abs_error_boundary( 0 )
    , avg_abs_error_boundary( 0 )
    , points_checked_boundary( 0 )
  {
  }
};

// Structure to hold all derivative errors for 1D spline
struct DerivativeErrors1D
{
  DerivativeErrorStats D;   // First derivative
  DerivativeErrorStats DD;  // Second derivative

  DerivativeErrors1D() : D(), DD() {}
};

// Structure to hold spline min/max results
struct SplineResult
{
  string    name;
  real_type x_min;
  real_type x_max;
  real_type y_min;
  real_type y_max;
  real_type x_at_min;
  real_type x_at_max;
};

// ============================================================================
// FINITE DIFFERENCE FUNCTIONS FOR 1D
// ============================================================================

namespace FiniteDifferences1D
{
  constexpr real_type eps = std::numeric_limits<real_type>::epsilon();

  // Adaptive step size for first derivative
  inline real_type h_first( real_type x )
  {
    real_type scale = std::max( real_type( 1 ), std::abs( x ) );
    return std::cbrt( eps ) * scale;  // ~ 6e-6 for double
  }

  // Adaptive step size for second derivative
  inline real_type h_second( real_type x )
  {
    real_type scale = std::max( real_type( 1 ), std::abs( x ) );
    return std::sqrt( std::sqrt( eps ) ) * scale;  // ~ 1e-4 for double
  }

  // Fourth-order central difference for first derivative
  template <typename SplineType> real_type dx( const SplineType & spline, real_type x )
  {
    real_type h  = h_first( x );
    real_type x0 = x;
    real_type x1 = x + h / 2;
    real_type x2 = x + h;
    real_type x3 = x - h / 2;
    real_type x4 = x - h;
    real_type y0 = spline( x0 );
    real_type y1 = spline( x1 );
    real_type y2 = spline( x2 );
    real_type y3 = spline( x3 );
    real_type y4 = spline( x4 );
    return Utils::first_derivative_5p( x0, y0, x1, y1, x2, y2, x3, y3, x4, y4 );
  }

  // Second-order central difference for second derivative
  template <typename SplineType> real_type dxx( const SplineType & spline, real_type x )
  {
    real_type h  = h_second( x );
    real_type x0 = x;
    real_type x1 = x + h / 2;
    real_type x2 = x + h;
    real_type x3 = x - h / 2;
    real_type x4 = x - h;
    real_type y0 = spline( x0 );
    real_type y1 = spline( x1 );
    real_type y2 = spline( x2 );
    real_type y3 = spline( x3 );
    real_type y4 = spline( x4 );
    return Utils::second_derivative_5p( x0, y0, x1, y1, x2, y2, x3, y3, x4, y4 );
  }
}  // namespace FiniteDifferences1D

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

// Function to print colored header
void print_header( const string & title )
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "╔═══════════════════════════════════════════════════════════╗\n"
    "║{:^59}║\n"
    "╚═══════════════════════════════════════════════════════════╝\n",
    title );
}

// Helper function to format error with color and location
string format_cell_with_stats(
  real_type      max_err,
  real_type      avg_err,
  size_t         points,
  const string & location,
  real_type      tol_max,
  real_type      tol_avg )
{
  // Determine color based on max error
  if ( max_err < tol_max && avg_err < tol_avg )
  {
    return fmt::format(
      fg( fmt::color::green ) | fmt::emphasis::bold,
      "{:52}",
      fmt::format( "{:<6} max:{:>8.1e} avg:{:>8.1e} ({:>3} pts)", location, max_err, avg_err, points ) );
  }
  else if ( max_err < 1e-4 && avg_err < 1e-5 )
  {
    return fmt::format(
      fg( fmt::color::light_green ),
      "{:52}",
      fmt::format( "{:<6} max:{:>8.1e} avg:{:>8.1e} ({:>3} pts)", location, max_err, avg_err, points ) );
  }
  else if ( max_err < 1e-2 && avg_err < 1e-3 )
  {
    return fmt::format(
      fg( fmt::color::yellow ),
      "{:52}",
      fmt::format( "{:<6} max:{:>8.1e} avg:{:>8.1e} ({:>3} pts)", location, max_err, avg_err, points ) );
  }
  else if ( max_err < 1.0 && avg_err < 0.1 )
  {
    return fmt::format(
      fg( fmt::color::orange ),
      "{:52}",
      fmt::format( "{:<6} max:{:>8.1e} avg:{:>8.1e} ({:>3} pts)", location, max_err, avg_err, points ) );
  }
  else
  {
    return fmt::format(
      fg( fmt::color::red ) | fmt::emphasis::bold,
      "{:52}",
      fmt::format( "{:<6} max:{:>8.1e} avg:{:>8.1e} ({:>3} pts)", location, max_err, avg_err, points ) );
  }
}

// Function to print derivative error table
void print_derivative_table_1D(
  const string &                                   derivative_name,
  const vector<pair<string, DerivativeErrors1D>> & errors,
  integer                                          dataset )
{
  auto color = fg( fmt::color::light_blue ) | fmt::emphasis::bold;
  fmt::print(
    color,
    "\n"
    "┌─{:─^132}─┐\n"
    "│ {:^132} │\n"
    "├─{:─^22}─┬─{:─^52}─┬─{:─^52}─┤\n"
    "│ {:^22} │ {:^52} │ {:^52} │\n"
    "├─{:─^22}─┼─{:─^52}─┼─{:─^52}─┤\n",
    "",
    fmt::format( "DATASET {} - DERIVATIVE {} - FINITE DIFFERENCE ERRORS", dataset, derivative_name ),
    "",
    "",
    "",
    "Spline Type",
    "Interior Points (within segments)",
    "Boundary Points (at knots)",
    "",
    "",
    "" );

  // Tolerances
  real_type first_deriv_tol_max_interior = 1e-6;
  real_type first_deriv_tol_avg_interior = 1e-7;
  real_type first_deriv_tol_max_boundary = 1e-5;
  real_type first_deriv_tol_avg_boundary = 1e-6;

  real_type second_deriv_tol_max_interior = 1e-3;
  real_type second_deriv_tol_avg_interior = 1e-4;
  real_type second_deriv_tol_max_boundary = 1e-2;
  real_type second_deriv_tol_avg_boundary = 1e-3;

  // Table rows
  for ( size_t i = 0; i < errors.size(); ++i )
  {
    const auto & [name, err]           = errors[i];
    const DerivativeErrorStats & stats = ( derivative_name == "D" ) ? err.D : err.DD;

    auto row_color = ( i % 2 == 0 ) ? fg( fmt::color::light_green ) : fg( fmt::color::light_blue );

    // Choose tolerances based on derivative
    real_type tol_max_int = ( derivative_name == "D" ) ? first_deriv_tol_max_interior : second_deriv_tol_max_interior;
    real_type tol_avg_int = ( derivative_name == "D" ) ? first_deriv_tol_avg_interior : second_deriv_tol_avg_interior;
    real_type tol_max_bnd = ( derivative_name == "D" ) ? first_deriv_tol_max_boundary : second_deriv_tol_max_boundary;
    real_type tol_avg_bnd = ( derivative_name == "D" ) ? first_deriv_tol_avg_boundary : second_deriv_tol_avg_boundary;

    // Interior errors
    string interior_cell = format_cell_with_stats(
      stats.max_abs_error_interior,
      stats.avg_abs_error_interior,
      stats.points_checked_interior,
      "Int",
      tol_max_int,
      tol_avg_int );

    // Boundary errors
    string boundary_cell = format_cell_with_stats(
      stats.max_abs_error_boundary,
      stats.avg_abs_error_boundary,
      stats.points_checked_boundary,
      "Bnd",
      tol_max_bnd,
      tol_avg_bnd );

    fmt::print( row_color, "│ {:<22} ", name );
    fmt::print( "│ {} ", interior_cell );
    fmt::print( "│ {} │\n", boundary_cell );

    if ( i < errors.size() - 1 ) { fmt::print( color, "├─{:─^22}─┼─{:─^52}─┼─{:─^52}─┤\n", "", "", "" ); }
  }

  fmt::print( color, "└─{:─^22}─┴─{:─^52}─┴─{:─^52}─┘\n", "", "", "" );
}

// ============================================================================
// DERIVATIVE CHECK FUNCTION FOR 1D SPLINES
// ============================================================================

template <typename SplineType> DerivativeErrors1D check_derivatives_1D(
  SplineType &      spline,
  const string &    spline_name,
  const real_type * knots,
  integer           n_knots,
  integer           points_per_segment = 10 )
{
  DerivativeErrors1D errors;

  // Initialize statistics
  auto init_stats = []() -> DerivativeErrorStats
  {
    DerivativeErrorStats stats;
    stats.max_abs_error_interior  = 0;
    stats.avg_abs_error_interior  = 0;
    stats.points_checked_interior = 0;
    stats.max_abs_error_boundary  = 0;
    stats.avg_abs_error_boundary  = 0;
    stats.points_checked_boundary = 0;
    return stats;
  };

  errors.D  = init_stats();
  errors.DD = init_stats();

  real_type x_min = spline.x_min();
  real_type x_max = spline.x_max();

  // Tolerance for boundary detection (relative to segment length)
  real_type boundary_tol = 1e-3 * ( x_max - x_min );

  fmt::print(
    fg( fmt::color::blue ),
    "   Checking derivatives for {} ({} knots, {} segments)...",
    spline_name,
    n_knots,
    n_knots - 1 );

  size_t total_points_checked = 0;

  // Test points in each segment
  for ( integer seg = 0; seg < n_knots - 1; ++seg )
  {
    real_type seg_start  = knots[seg];
    real_type seg_end    = knots[seg + 1];
    real_type seg_length = seg_end - seg_start;

    // Generate test points in this segment
    for ( integer p = 0; p <= points_per_segment; ++p )
    {
      real_type x;
      if ( p == 0 )
        x = seg_start + boundary_tol;  // Just after knot
      else if ( p == points_per_segment )
        x = seg_end - boundary_tol;  // Just before next knot
      else
        x = seg_start + p * seg_length / points_per_segment;

      // Check if point is interior or boundary
      bool is_boundary = false;
      for ( integer k = 0; k < n_knots; ++k )
      {
        if ( std::abs( x - knots[k] ) < boundary_tol )
        {
          is_boundary = true;
          break;
        }
      }

      // Ensure we're within spline domain
      if ( x < x_min || x > x_max ) continue;

      // Get analytical derivatives
      real_type D_analytic  = spline.D( x );
      real_type DD_analytic = spline.DD( x );

      // Get finite difference approximations
      real_type D_fd  = FiniteDifferences1D::dx( spline, x );
      real_type DD_fd = FiniteDifferences1D::dxx( spline, x );

      // Calculate absolute errors
      real_type abs_error_D  = std::abs( D_analytic - D_fd );
      real_type abs_error_DD = std::abs( DD_analytic - DD_fd );

      // Update statistics based on boundary status
      if ( is_boundary )
      {
        // First derivative boundary stats
        errors.D.max_abs_error_boundary = std::max( errors.D.max_abs_error_boundary, abs_error_D );
        errors.D.avg_abs_error_boundary += abs_error_D;
        errors.D.points_checked_boundary++;

        // Second derivative boundary stats
        errors.DD.max_abs_error_boundary = std::max( errors.DD.max_abs_error_boundary, abs_error_DD );
        errors.DD.avg_abs_error_boundary += abs_error_DD;
        errors.DD.points_checked_boundary++;
      }
      else
      {
        // First derivative interior stats
        errors.D.max_abs_error_interior = std::max( errors.D.max_abs_error_interior, abs_error_D );
        errors.D.avg_abs_error_interior += abs_error_D;
        errors.D.points_checked_interior++;

        // Second derivative interior stats
        errors.DD.max_abs_error_interior = std::max( errors.DD.max_abs_error_interior, abs_error_DD );
        errors.DD.avg_abs_error_interior += abs_error_DD;
        errors.DD.points_checked_interior++;
      }

      total_points_checked++;
    }
  }

  // Calculate averages
  auto calculate_averages = []( DerivativeErrorStats & stats )
  {
    if ( stats.points_checked_interior > 0 ) stats.avg_abs_error_interior /= static_cast<real_type>( stats.points_checked_interior );
    if ( stats.points_checked_boundary > 0 ) stats.avg_abs_error_boundary /= static_cast<real_type>( stats.points_checked_boundary );
  };

  calculate_averages( errors.D );
  calculate_averages( errors.DD );

  fmt::print(
    fg( fmt::color::green ),
    " ✓ ({} points total, {} interior, {} boundary)\n",
    total_points_checked,
    errors.D.points_checked_interior + errors.DD.points_checked_interior,
    errors.D.points_checked_boundary + errors.DD.points_checked_boundary );

  return errors;
}

// ============================================================================
// CONVERGENCE TEST FUNCTIONS - MODIFICATE
// ============================================================================

// Compute maximum error between spline and exact function
real_type compute_max_error(
  Splines::Spline *         spline,
  const vector<real_type> & eval_points,
  TestFunctionPtr           exact_func )
{
  real_type max_error = 0.0;
  for ( const auto & x : eval_points )
  {
    real_type error = abs( spline->eval( x ) - exact_func( x ) );
    if ( error > max_error ) max_error = error;
  }
  return max_error;
}

// Perform convergence test with hierarchical refinement for non-uniform meshes
ConvergenceResult test_convergence_nonuniform_hierarchical(
  Spline *        spline,
  const string &  spline_name,
  const string &  test_func_name,
  real_type       a,
  real_type       b,
  TestFunctionPtr test_func,
  integer         base_N     = 5,  // Starting mesh size
  integer         num_levels = 7,  // Number of refinement levels (5, 9, 17, 33, 65, 129, 257)
  integer         Neval      = 10000 )
{
  ConvergenceResult result( spline_name, test_func_name );

  // Generate base non-uniform mesh
  vector<real_type> current_mesh = generate_base_nonuniform_mesh( a, b, base_N );

  // Generate evaluation points once
  vector<real_type> eval_points = generate_uniform_mesh( a, b, Neval );

  for ( integer level = 0; level < num_levels; ++level )
  {
    integer N = static_cast<integer>( current_mesh.size() );

    // Evaluate function at mesh points
    vector<real_type> Y( N );
    for ( integer i = 0; i < N; ++i ) { Y[i] = test_func( current_mesh[i] ); }

    // Build spline
    spline->build( current_mesh, Y );

    // Compute maximum error
    real_type error = compute_max_error( spline, eval_points, test_func );

    // Compute mesh size (average interval length)
    real_type h = compute_mesh_size( current_mesh );

    // Store results
    result.N_values.push_back( N );
    result.errors.push_back( error );
    result.mesh_sizes.push_back( h );

    // Refine mesh for next level (except for the last iteration)
    if ( level < num_levels - 1 ) { current_mesh = refine_mesh_by_insertion( current_mesh ); }
  }

  // Compute convergence orders using actual mesh sizes
  auto eps = 1e-12;
  for ( size_t i = 0; i < result.errors.size() - 1; ++i )
  {
    real_type h1 = result.mesh_sizes[i];
    real_type h2 = result.mesh_sizes[i + 1];

    // Check for very small errors
    if ( result.errors[i] < eps && result.errors[i + 1] < eps )
    {
      result.orders.push_back( std::numeric_limits<real_type>::infinity() );
    }
    else if ( result.errors[i + 1] < eps ) { result.orders.push_back( std::numeric_limits<real_type>::infinity() ); }
    else if ( result.errors[i] < eps ) { result.orders.push_back( std::numeric_limits<real_type>::quiet_NaN() ); }
    else
    {
      // Compute order using actual mesh sizes
      real_type order = log( result.errors[i] / result.errors[i + 1] ) / log( h1 / h2 );
      result.orders.push_back( order );
    }
  }
  // Add NaN for the last point
  result.orders.push_back( std::numeric_limits<real_type>::quiet_NaN() );

  return result;
}

// Perform convergence test for a specific spline type and test function
ConvergenceResult test_convergence(
  Spline *        spline,
  const string &  spline_name,
  const string &  test_func_name,
  bool            uniform_mesh,
  real_type       a,
  real_type       b,
  TestFunctionPtr test_func,
  integer         Neval = 10000 )
{
  if ( uniform_mesh )
  {
    // Original uniform mesh test
    ConvergenceResult result( spline_name, test_func_name );

    // Sequence of mesh sizes (increasing)
    vector<integer> N_seq = { 5, 9, 17, 33, 65, 129, 257 };

    // Generate evaluation points once - using more points for accurate error estimation
    vector<real_type> eval_points = generate_uniform_mesh( a, b, Neval );

    for ( integer N : N_seq )
    {
      // Generate mesh
      vector<real_type> X = generate_uniform_mesh( a, b, N );

      // Evaluate function at mesh points
      vector<real_type> Y( N );
      for ( integer i = 0; i < N; ++i ) { Y[i] = test_func( X[i] ); }

      // Build spline
      spline->build( X, Y );

      // Compute maximum error
      real_type error = compute_max_error( spline, eval_points, test_func );

      // Compute mesh size (average interval length)
      real_type h = compute_mesh_size( X );

      result.N_values.push_back( N );
      result.errors.push_back( error );
      result.mesh_sizes.push_back( h );
    }

    // Compute convergence orders with protection against division by zero
    auto eps = 1e-12;
    for ( size_t i = 0; i < result.errors.size() - 1; ++i )
    {
      real_type h1 = result.mesh_sizes[i];
      real_type h2 = result.mesh_sizes[i + 1];

      // Check for very small errors (machine precision)
      if ( result.errors[i] < eps && result.errors[i + 1] < eps )
      {
        // Both errors are essentially zero, convergence is perfect
        result.orders.push_back( std::numeric_limits<real_type>::infinity() );
      }
      else if ( result.errors[i + 1] < eps )
      {
        // Only the finer mesh error is zero, perfect convergence
        result.orders.push_back( std::numeric_limits<real_type>::infinity() );
      }
      else if ( result.errors[i] < eps )
      {
        // Coarse mesh error is zero but fine mesh is not (shouldn't happen)
        result.orders.push_back( std::numeric_limits<real_type>::quiet_NaN() );
      }
      else
      {
        // Normal case
        real_type order = log( result.errors[i] / result.errors[i + 1] ) / log( h1 / h2 );
        result.orders.push_back( order );
      }
    }
    // Add NaN for the last point (no next point to compare with)
    result.orders.push_back( std::numeric_limits<real_type>::quiet_NaN() );

    return result;
  }
  else
  {
    // Use hierarchical refinement for non-uniform meshes
    return test_convergence_nonuniform_hierarchical(
      spline,
      spline_name,
      test_func_name,
      a,
      b,
      test_func,
      5,
      7,
      Neval );
  }
}

// Print convergence table for a specific test function
void print_convergence_table_for_test(
  const vector<ConvergenceResult> & results,
  const string &                    mesh_type,
  const TestFunctionInfo &          test_func,
  [[maybe_unused]] integer          func_index )
{
  print_header( fmt::format( "CONVERGENCE TEST - {} - {} MESH", test_func.name, mesh_type ) );

  // Determine if we should show mesh sizes
  bool show_mesh_size = true;  // Always show mesh size now

  if ( show_mesh_size )
  {
    // Table with mesh size column
    fmt::print(
      fg( fmt::color::yellow ) | fmt::emphasis::bold,
      "\n"
      "┌────────────────────────┬──────────────┬──────────────┬──────────────┬──────────────┐\n"
      "│ {:^22} │ {:^12} │ {:^12} │ {:^12} │ {:^12} │\n"
      "├────────────────────────┼──────────────┼──────────────┼──────────────┼──────────────┤\n",
      "Spline Type",
      "N points",
      "Mesh Size",
      "Max Error",
      "Order" );
  }
  else
  {
    // Table without mesh size column
    fmt::print(
      fg( fmt::color::yellow ) | fmt::emphasis::bold,
      "\n"
      "┌────────────────────────┬──────────────┬──────────────┬──────────────┐\n"
      "│ {:^22} │ {:^12} │ {:^12} │ {:^12} │\n"
      "├────────────────────────┼──────────────┼──────────────┼──────────────┤\n",
      "Spline Type",
      "N points",
      "Max Error",
      "Order" );
  }

  // Table rows
  for ( size_t i = 0; i < results.size(); ++i )
  {
    const auto & res = results[i];

    for ( size_t j = 0; j < res.N_values.size(); ++j )
    {
      auto row_color = ( j % 2 == 0 ) ? fg( fmt::color::light_green ) : fg( fmt::color::light_blue );

      fmt::print( row_color, "│ {:<22} ", ( j == 0 ) ? res.spline_name : "" );
      fmt::print( row_color, "│ {:>12} ", res.N_values[j] );

      if ( show_mesh_size )
      {
        // Show mesh size
        real_type h = res.mesh_sizes[j];
        fmt::print( row_color, "│ {:>12.3e} ", h );
      }

      // Color code error based on magnitude
      real_type       error       = res.errors[j];
      fmt::text_style error_style = fg( fmt::color::white );
      if ( error < 1e-10 )
        error_style = fg( fmt::color::green ) | fmt::emphasis::bold;
      else if ( error < 1e-8 )
        error_style = fg( fmt::color::light_green );
      else if ( error < 1e-6 )
        error_style = fg( fmt::color::yellow );
      else if ( error < 1e-4 )
        error_style = fg( fmt::color::orange );
      else
        error_style = fg( fmt::color::red ) | fmt::emphasis::bold;

      fmt::print( error_style, "│ {:>12.3e} ", error );

      // Print order (if available)
      if ( !isnan( res.orders[j] ) )
      {
        if ( isinf( res.orders[j] ) )
        {
          fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "│ {:>12} ", "INF" );
        }
        else
        {
          real_type       order       = res.orders[j];
          fmt::text_style order_style = fg( fmt::color::white );
          if ( order > 5.0 )
            order_style = fg( fmt::color::green ) | fmt::emphasis::bold;
          else if ( order > 3.0 )
            order_style = fg( fmt::color::light_green );
          else if ( order > 1.5 )
            order_style = fg( fmt::color::yellow );
          else if ( order > 0.5 )
            order_style = fg( fmt::color::orange );
          else
            order_style = fg( fmt::color::red ) | fmt::emphasis::bold;

          fmt::print( order_style, "│ {:>12.2f} ", order );
        }
      }
      else
      {
        fmt::print( row_color, "│ {:>12} ", "---" );
      }

      fmt::print( "│\n" );
    }

    if ( i < results.size() - 1 )
    {
      if ( show_mesh_size )
      {
        fmt::print(
          fg( fmt::color::yellow ),
          "├────────────────────────┼──────────────┼──────────────┼──────────────┼──────────────┤\n" );
      }
      else
      {
        fmt::print(
          fg( fmt::color::yellow ),
          "├────────────────────────┼──────────────┼──────────────┼──────────────┤\n" );
      }
    }
  }

  if ( show_mesh_size )
  {
    fmt::print(
      fg( fmt::color::yellow ) | fmt::emphasis::bold,
      "└────────────────────────┴──────────────┴──────────────┴──────────────┴──────────────┘\n" );
  }
  else
  {
    fmt::print(
      fg( fmt::color::yellow ) | fmt::emphasis::bold,
      "└────────────────────────┴──────────────┴──────────────┴──────────────┘\n" );
  }

  // Summary statistics
  fmt::print( fg( fmt::color::cyan ), "\n📊 SUMMARY STATISTICS FOR {}:\n", test_func.name );

  for ( const auto & res : results )
  {
    // Compute average order (excluding first, last, and infinite/NaN values)
    real_type avg_order = 0.0;
    integer   count     = 0;
    for ( size_t i = 1; i < res.orders.size() - 2; ++i )
    {
      if ( !isnan( res.orders[i] ) && !isinf( res.orders[i] ) )
      {
        avg_order += res.orders[i];
        count++;
      }
    }

    if ( count > 0 )
    {
      avg_order /= count;
      fmt::print( fg( fmt::color::white ), "  {}: Average convergence order = {:.3f}\n", res.spline_name, avg_order );
    }
    else if ( res.errors.back() < 1e-10 )
    {
      // If the last error is very small, the spline interpolates exactly
      fmt::print(
        fg( fmt::color::green ) | fmt::emphasis::bold,
        "  {}: Exact interpolation (error < 1e-10)\n",
        res.spline_name );
    }
  }

  // Theoretical expectations based on polynomial degree
  fmt::print( fg( fmt::color::cyan ), "\n🎯 THEORETICAL EXPECTATIONS:\n" );

  // Determine polynomial degree from function name
  if ( test_func.name.find( "deg 1" ) != string::npos )
  {
    fmt::print(
      fg( fmt::color::white ),
      "  Linear polynomial (deg 1):\n"
      "    • Linear splines should interpolate exactly (error ~0)\n"
      "    • Cubic and higher splines should interpolate exactly (error ~0)\n" );
  }
  else if ( test_func.name.find( "deg 2" ) != string::npos )
  {
    fmt::print(
      fg( fmt::color::white ),
      "  Quadratic polynomial (deg 2):\n"
      "    • Linear splines: O(h²) convergence\n"
      "    • Cubic and higher splines should interpolate exactly (error ~0)\n" );
  }
  else if ( test_func.name.find( "deg 3" ) != string::npos )
  {
    fmt::print(
      fg( fmt::color::white ),
      "  Cubic polynomial (deg 3):\n"
      "    • Linear splines: O(h²) convergence\n"
      "    • Cubic splines should interpolate exactly (error ~0)\n"
      "    • Quintic and higher splines should interpolate exactly (error ~0)\n" );
  }
  else if ( test_func.name.find( "deg 4" ) != string::npos )
  {
    fmt::print(
      fg( fmt::color::white ),
      "  Quartic polynomial (deg 4):\n"
      "    • Linear splines: O(h²) convergence\n"
      "    • Cubic splines: O(h⁴) convergence\n"
      "    • Quintic splines should interpolate exactly (error ~0)\n" );
  }
  else if ( test_func.name.find( "Exp" ) != string::npos )
  {
    fmt::print(
      fg( fmt::color::white ),
      "  Exponential + trigonometric function:\n"
      "    • Linear splines: O(h²) convergence\n"
      "    • Cubic splines: O(h⁴) convergence\n"
      "    • Quintic splines: O(h⁶) convergence\n" );
  }
}

// ============================================================================
// FUNCTION TO PRINT RESULTS TABLE
// ============================================================================

// Function to print a table of results
void print_results_table( const vector<SplineResult> & results, integer dataset )
{
  // Table header
  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "\n"
    "┌────────────────────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐\n"
    "│ {:^22} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │\n"
    "├────────────────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤\n",
    "Spline Type",
    "x_min",
    "x_max",
    "y_min",
    "y_max",
    "x@min",
    "x@max" );

  // Table rows with alternating colors
  for ( size_t i = 0; i < results.size(); ++i )
  {
    const auto & r         = results[i];
    auto         row_color = ( i % 2 == 0 ) ? fg( fmt::color::light_green ) : fg( fmt::color::light_blue );

    fmt::print( row_color, "│ {:<22} ", r.name );
    fmt::print( row_color, "│ {:>8.3f} ", r.x_min );
    fmt::print( row_color, "│ {:>8.3f} ", r.x_max );

    // Highlight y_min and y_max with different colors based on value
    auto y_min_color = ( r.y_min < 0 ) ? fg( fmt::color::red ) : fg( fmt::color::green );
    auto y_max_color = ( r.y_max < 0 ) ? fg( fmt::color::red ) : fg( fmt::color::green );

    fmt::print( y_min_color, "│ {:>8.3f} ", r.y_min );
    fmt::print( y_max_color, "│ {:>8.3f} ", r.y_max );
    fmt::print( row_color, "│ {:>8.3f} ", r.x_at_min );
    fmt::print( row_color, "│ {:>8.3f} │\n", r.x_at_max );
  }

  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "└────────────────────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘\n" );

  fmt::print(
    fg( fmt::color::gray ) | fmt::emphasis::italic,
    "Dataset {} - {} points analyzed\n",
    dataset,
    nn[dataset] );
}

// ============================================================================
// DERIVATIVE CONSISTENCY TEST: SPLINE VS AUTODIFF
// ============================================================================

void test_derivative_consistency_for_function(
  const TestFunctionInfo & test_func,
  const string &           spline_name,
  Splines::Spline *        spline_ptr,
  integer                  N_fine = 64 )
{
  // Genera mesh e valori della funzione
  auto              X_fine = generate_uniform_mesh( test_func.a, test_func.b, N_fine );
  vector<real_type> Y_fine( N_fine );
  for ( integer i = 0; i < N_fine; ++i ) { Y_fine[i] = test_func.func( X_fine[i] ); }

  // Costruisce la spline
  spline_ptr->build( X_fine, Y_fine );

  fmt::print(
    fg( fmt::color::white ),
    "\nDerivative consistency test ({} spline, {}, N={}):\n",
    spline_name,
    test_func.name,
    N_fine );

  fmt::print(
    fg( fmt::color::cyan ),
    "{:>10} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n",
    "x",
    "D(spline)",
    "D(autodiff)",
    "|diff|",
    "DD(spline)",
    "DD(autodiff)",
    "|diff|" );

  // Punti di test nel dominio
  integer           num_test_points = 7;
  vector<real_type> test_points( num_test_points );
  for ( integer i = 0; i < num_test_points; ++i )
  {
    test_points[i] = test_func.a + i * ( test_func.b - test_func.a ) / ( num_test_points - 1 );
  }

  // Statistiche delle differenze
  real_type max_D_diff = 0.0, max_DD_diff = 0.0;
  real_type avg_D_diff = 0.0, avg_DD_diff = 0.0;
  integer   num_valid_points = 0;

  for ( real_type x : test_points )
  {
    // Valori dalla spline usando eval_D ed eval_DD
    real_type D_spline  = spline_ptr->eval_D( x );
    real_type DD_spline = spline_ptr->eval_DD( x );

    // Usa autodiff::dual1st per la derivata prima
    autodiff::dual1st xd1, yd1;
    xd1.val              = x;
    xd1.grad             = 1;
    yd1                  = spline_ptr->eval( xd1 );
    real_type D_autodiff = yd1.grad;

    // Usa autodiff::dual2nd per la derivata seconda
    autodiff::dual2nd xd2, yd2;
    xd2.val.val           = x;
    xd2.val.grad          = 1;
    xd2.grad.val          = 1;
    xd2.grad.grad         = 0;
    yd2                   = spline_ptr->eval( xd2 );
    real_type DD_autodiff = yd2.grad.grad;

    // Calcola differenze
    real_type D_diff  = abs( D_spline - D_autodiff );
    real_type DD_diff = abs( DD_spline - DD_autodiff );

    // Aggiorna statistiche
    max_D_diff  = max( max_D_diff, D_diff );
    max_DD_diff = max( max_DD_diff, DD_diff );
    avg_D_diff += D_diff;
    avg_DD_diff += DD_diff;
    num_valid_points++;

    // Stile per le differenze
    auto D_diff_style = ( D_diff < 1e-12 )  ? fg( fmt::color::green )
                        : ( D_diff < 1e-9 ) ? fg( fmt::color::yellow )
                        : ( D_diff < 1e-6 ) ? fg( fmt::color::orange )
                                            : fg( fmt::color::red );

    auto DD_diff_style = ( DD_diff < 1e-10 )  ? fg( fmt::color::green )
                         : ( DD_diff < 1e-7 ) ? fg( fmt::color::yellow )
                         : ( DD_diff < 1e-4 ) ? fg( fmt::color::orange )
                                              : fg( fmt::color::red );

    // Stampa la riga
    fmt::print( fg( fmt::color::white ), "{:>10.3f} ", x );

    // Derivata prima
    fmt::print( fg( fmt::color::light_blue ), "{:>15.6f} ", D_spline );
    fmt::print( fg( fmt::color::light_green ), "{:>15.6f} ", D_autodiff );
    fmt::print( D_diff_style, "{:>15.2e} ", D_diff );

    // Derivata seconda
    fmt::print( fg( fmt::color::light_blue ), "{:>15.6f} ", DD_spline );
    fmt::print( fg( fmt::color::light_green ), "{:>15.6f} ", DD_autodiff );
    fmt::print( DD_diff_style, "{:>15.2e}\n", DD_diff );
  }

  // Calcola medie
  if ( num_valid_points > 0 )
  {
    avg_D_diff /= num_valid_points;
    avg_DD_diff /= num_valid_points;
  }

  // Riepilogo delle differenze
  fmt::print( fg( fmt::color::white ), "\n  Consistency summary:\n" );

  // Controlla la coerenza delle derivate prime
  fmt::print( fg( fmt::color::cyan ), "    First derivative (D):\n" );
  if ( max_D_diff < 1e-12 )
  {
    fmt::print( fg( fmt::color::green ), "      PERFECT match! Max diff = {:.2e}\n", max_D_diff );
  }
  else if ( max_D_diff < 1e-9 )
  {
    fmt::print( fg( fmt::color::green ), "      EXCELLENT consistency. Max diff = {:.2e}\n", max_D_diff );
  }
  else if ( max_D_diff < 1e-6 )
  {
    fmt::print( fg( fmt::color::yellow ), "      GOOD consistency. Max diff = {:.2e}\n", max_D_diff );
  }
  else if ( max_D_diff < 1e-3 )
  {
    fmt::print( fg( fmt::color::orange ), "      FAIR consistency. Max diff = {:.2e}\n", max_D_diff );
  }
  else
  {
    fmt::print( fg( fmt::color::red ), "      POOR consistency. Max diff = {:.2e}\n", max_D_diff );
  }
  fmt::print( fg( fmt::color::white ), "      Average diff = {:.2e}\n", avg_D_diff );

  // Controlla la coerenza delle derivate seconde
  fmt::print( fg( fmt::color::cyan ), "    Second derivative (DD):\n" );
  if ( max_DD_diff < 1e-10 )
  {
    fmt::print( fg( fmt::color::green ), "      PERFECT match! Max diff = {:.2e}\n", max_DD_diff );
  }
  else if ( max_DD_diff < 1e-7 )
  {
    fmt::print( fg( fmt::color::green ), "      EXCELLENT consistency. Max diff = {:.2e}\n", max_DD_diff );
  }
  else if ( max_DD_diff < 1e-4 )
  {
    fmt::print( fg( fmt::color::yellow ), "      GOOD consistency. Max diff = {:.2e}\n", max_DD_diff );
  }
  else if ( max_DD_diff < 1e-1 )
  {
    fmt::print( fg( fmt::color::orange ), "      FAIR consistency. Max diff = {:.2e}\n", max_DD_diff );
  }
  else
  {
    fmt::print( fg( fmt::color::red ), "      POOR consistency. Max diff = {:.2e}\n", max_DD_diff );
  }
  fmt::print( fg( fmt::color::white ), "      Average diff = {:.2e}\n", avg_DD_diff );

  // Confronto finale
  fmt::print( fg( fmt::color::cyan ), "    Overall assessment:\n" );
  if ( max_D_diff < 1e-9 && max_DD_diff < 1e-7 )
  {
    fmt::print( fg( fmt::color::green ), "      ✓ Spline derivatives are fully consistent with autodiff\n" );
  }
  else if ( max_D_diff < 1e-6 && max_DD_diff < 1e-4 )
  {
    fmt::print( fg( fmt::color::yellow ), "      ~ Spline derivatives are reasonably consistent with autodiff\n" );
  }
  else
  {
    fmt::print( fg( fmt::color::red ), "      ✗ Significant differences between spline and autodiff derivatives\n" );

    // Suggerimenti diagnostici
    if ( max_DD_diff > 10 * max_D_diff )
    {
      fmt::print( fg( fmt::color::orange ), "        Note: Second derivatives show larger discrepancies\n" );
    }
    if ( spline_ptr->type() == Splines::SplineType1D::LINEAR )
    {
      fmt::print( fg( fmt::color::orange ), "        Note: Linear splines have zero second derivative\n" );
    }
  }
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================

int main()
{
  print_header( "SPLINE INTERPOLATION TEST SUITE" );

  // Create spline objects
  LinearSpline   li;
  ConstantSpline co;
  AkimaSpline    ak;
  CubicSpline    cs;
  VanLeerSpline  vl;
  PchipSpline    pc;
  QuinticSpline  qs_cubic( Spline_sub_type::CUBIC );
  QuinticSpline  qs_pchip( Spline_sub_type::PCHIP );
  QuinticSpline  qs_akima( Spline_sub_type::AKIMA );
  QuinticSpline  qs_vanleer( Spline_sub_type::VANLEER );

  // Open output files
  ofstream file_li, file_co, file_ak, file_cs, file_vl, file_pc, file_qs_cubic, file_qs_pchip, file_qs_akima,
    file_qs_vanleer;

  for ( integer k = 0; k < 6; ++k )
  {
    fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n\n📊 DATASET {} ANALYSIS\n", k );

    real_type * xx{ nullptr };
    real_type * yy{ nullptr };

    switch ( k )
    {
      case 0:
        xx = xx0;
        yy = yy0;
        break;
      case 1:
        xx = xx1;
        yy = yy1;
        break;
      case 2:
        xx = xx2;
        yy = yy2;
        break;
      case 3:
        xx = xx3;
        yy = yy3;
        break;
      case 4:
        xx = xx4;
        yy = yy4;
        break;
      case 5:
        xx = xx5;
        yy = yy5;
        break;
    }

    // Open output files
    string fname;
    fname = fmt::format( "out/Linear{}.txt", k );
    file_li.open( fname.data() );
    fname = fmt::format( "out/Constant{}.txt", k );
    file_co.open( fname.data() );
    fname = fmt::format( "out/Akima{}.txt", k );
    file_ak.open( fname.data() );
    fname = fmt::format( "out/Cubic{}.txt", k );
    file_cs.open( fname.data() );
    fname = fmt::format( "out/VanLeer{}.txt", k );
    file_vl.open( fname.data() );
    fname = fmt::format( "out/Pchip{}.txt", k );
    file_pc.open( fname.data() );
    fname = fmt::format( "out/Quintic_CUBIC{}.txt", k );
    file_qs_cubic.open( fname.data() );
    fname = fmt::format( "out/Quintic_PCHIP{}.txt", k );
    file_qs_pchip.open( fname.data() );
    fname = fmt::format( "out/Quintic_AKIMA{}.txt", k );
    file_qs_akima.open( fname.data() );
    fname = fmt::format( "out/Quintic_VANLEER{}.txt", k );
    file_qs_vanleer.open( fname.data() );

    real_type xmin{ xx[0] };
    real_type xmax{ xx[nn[k] - 1] };

    vector<SplineResult>                     results;
    vector<pair<string, DerivativeErrors1D>> derivative_errors;

    // Lambda function to process each spline type
    auto process_spline = [&]( auto & spline, const string & name, ofstream & file ) -> DerivativeErrors1D
    {
      fmt::print( fg( fmt::color::cyan ), "\n🔧 Processing {}...\n", name );

      // Build spline
      spline.clear();
      spline.reserve( nn[k] );
      for ( integer i = 0; i < integer( nn[k] ); ++i ) spline.push_back( xx[i], yy[i] );
      spline.build();

      // Get min/max information
      integer   i_min_pos, i_max_pos;
      real_type x_min_pos, x_max_pos, y_min, y_max;
      spline.y_min_max( i_min_pos, x_min_pos, y_min, i_max_pos, x_max_pos, y_max );

      // Save results
      results.push_back( { name, spline.x_min(), spline.x_max(), y_min, y_max, x_min_pos, x_max_pos } );

      // Write spline data to file
      file << "x\ty\tDy\tDDy\n";
      for ( real_type x = xmin; x <= xmax; x += ( xmax - xmin ) / 1000 )
        fmt::print( file, "{}\t{}\t{}\t{}\n", x, spline( x ), spline.D( x ), spline.DD( x ) );
      file.close();

      // Check derivatives
      auto errors = check_derivatives_1D( spline, name, xx, nn[k], 10 );
      derivative_errors.emplace_back( name, errors );

      fmt::print( fg( fmt::color::green ), "   ✓ {} completed\n", name );
      return errors;
    };

    // Process all spline types
    fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n📈 Building splines and checking derivatives...\n" );

    process_spline( li, "Linear", file_li );
    process_spline( co, "Constant", file_co );
    process_spline( ak, "Akima", file_ak );
    process_spline( cs, "Cubic", file_cs );
    process_spline( vl, "VanLeer", file_vl );
    process_spline( pc, "Pchip", file_pc );
    process_spline( qs_cubic, "Quintic (CUBIC)", file_qs_cubic );
    process_spline( qs_pchip, "Quintic (PCHIP)", file_qs_pchip );
    process_spline( qs_akima, "Quintic (AKIMA)", file_qs_akima );
    process_spline( qs_vanleer, "Quintic (VANLEER)", file_qs_vanleer );

    // Print min/max results table
    print_results_table( results, k );

    // Print derivative error tables
    fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n📊 DERIVATIVE ERROR ANALYSIS - DATASET {}:\n", k );

    print_derivative_table_1D( "D", derivative_errors, k );
    print_derivative_table_1D( "DD", derivative_errors, k );

    // Summary of derivative checks
    fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n📋 DERIVATIVE CHECK SUMMARY - DATASET {}:\n", k );

    for ( const auto & [name, err] : derivative_errors )
    {
      fmt::print( fg( fmt::color::white ) | fmt::emphasis::bold, "   {}:\n", name );

      // First derivative summary
      bool first_deriv_interior_ok = ( err.D.max_abs_error_interior < 1e-6 && err.D.avg_abs_error_interior < 1e-7 );

      bool first_deriv_boundary_ok = ( err.D.max_abs_error_boundary < 1e-5 && err.D.avg_abs_error_boundary < 1e-6 );

      // Second derivative summary
      bool second_deriv_interior_ok = ( err.DD.max_abs_error_interior < 1e-3 && err.DD.avg_abs_error_interior < 1e-4 );

      bool second_deriv_boundary_ok = ( err.DD.max_abs_error_boundary < 1e-2 && err.DD.avg_abs_error_boundary < 1e-3 );

      // Print results with color coding
      auto print_check = []( bool ok, const string & text )
      {
        if ( ok )
          fmt::print( fg( fmt::color::green ), "     ✓ {}\n", text );
        else
          fmt::print( fg( fmt::color::yellow ), "     ⚠ {}\n", text );
      };

      print_check( first_deriv_interior_ok, "First derivative (interior)" );
      print_check( first_deriv_boundary_ok, "First derivative (boundary)" );
      print_check( second_deriv_interior_ok, "Second derivative (interior)" );
      print_check( second_deriv_boundary_ok, "Second derivative (boundary)" );
    }

    fmt::print( fg( fmt::color::gray ), "   Output files saved to: out/*{}.txt\n", k );
  }

  // ==========================================================================
  // CONVERGENCE TESTS - ALL SPLINE TYPES FOR 5 DIFFERENT TEST FUNCTIONS
  // ==========================================================================

  print_header( "CONVERGENCE ANALYSIS - 5 TEST FUNCTIONS" );

  // Test each function
  integer SZ = static_cast<integer>( test_functions.size() );
  for ( integer func_idx = 0; func_idx < SZ; ++func_idx )
  {
    const auto & test_func = test_functions[func_idx];

    fmt::print(
      fg( fmt::color::magenta ) | fmt::emphasis::bold,
      "\n\n🔬 TEST FUNCTION {}: {}\n",
      func_idx + 1,
      test_func.name );

    // Test uniform mesh convergence - ALL TYPES
    fmt::print( fg( fmt::color::cyan ), "\n📈 UNIFORM MESH CONVERGENCE:\n" );

    vector<ConvergenceResult> uniform_results;


    // Test ALL spline types on uniform mesh
    uniform_results.push_back(
      test_convergence( &li, "Linear", test_func.name, true, test_func.a, test_func.b, test_func.func ) );

    uniform_results.push_back(
      test_convergence( &co, "Constant", test_func.name, true, test_func.a, test_func.b, test_func.func ) );

    uniform_results.push_back(
      test_convergence( &cs, "Cubic", test_func.name, true, test_func.a, test_func.b, test_func.func ) );

    uniform_results.push_back(
      test_convergence( &ak, "Akima", test_func.name, true, test_func.a, test_func.b, test_func.func ) );

    uniform_results.push_back(
      test_convergence( &vl, "VanLeer", test_func.name, true, test_func.a, test_func.b, test_func.func ) );

    uniform_results.push_back(
      test_convergence( &pc, "Pchip", test_func.name, true, test_func.a, test_func.b, test_func.func ) );

    uniform_results.push_back( test_convergence(
      &qs_cubic,
      "Quintic (CUBIC)",
      test_func.name,
      true,
      test_func.a,
      test_func.b,
      test_func.func ) );

    uniform_results.push_back( test_convergence(
      &qs_pchip,
      "Quintic (PCHIP)",
      test_func.name,
      true,
      test_func.a,
      test_func.b,
      test_func.func ) );

    uniform_results.push_back( test_convergence(
      &qs_akima,
      "Quintic (AKIMA)",
      test_func.name,
      true,
      test_func.a,
      test_func.b,
      test_func.func ) );

    uniform_results.push_back( test_convergence(
      &qs_vanleer,
      "Quintic (VANLEER)",
      test_func.name,
      true,
      test_func.a,
      test_func.b,
      test_func.func ) );

    print_convergence_table_for_test( uniform_results, "UNIFORM", test_func, func_idx );

    // Test non-uniform mesh convergence - ALL TYPES
    fmt::print( fg( fmt::color::cyan ), "\n📈 NON-UNIFORM MESH CONVERGENCE:\n" );

    vector<ConvergenceResult> nonuniform_results;

    // Test ALL spline types on non-uniform mesh
    nonuniform_results.push_back(
      test_convergence( &li, "Linear", test_func.name, false, test_func.a, test_func.b, test_func.func ) );

    nonuniform_results.push_back(
      test_convergence( &co, "Constant", test_func.name, false, test_func.a, test_func.b, test_func.func ) );

    nonuniform_results.push_back(
      test_convergence( &cs, "Cubic", test_func.name, false, test_func.a, test_func.b, test_func.func ) );

    nonuniform_results.push_back(
      test_convergence( &ak, "Akima", test_func.name, false, test_func.a, test_func.b, test_func.func ) );

    nonuniform_results.push_back(
      test_convergence( &vl, "VanLeer", test_func.name, false, test_func.a, test_func.b, test_func.func ) );

    nonuniform_results.push_back(
      test_convergence( &pc, "Pchip", test_func.name, false, test_func.a, test_func.b, test_func.func ) );

    // Stesso pattern per i test non uniformi:
    nonuniform_results.push_back( test_convergence(
      &qs_cubic,
      "Quintic (CUBIC)",
      test_func.name,
      false,
      test_func.a,
      test_func.b,
      test_func.func ) );

    nonuniform_results.push_back( test_convergence(
      &qs_pchip,
      "Quintic (PCHIP)",
      test_func.name,
      false,
      test_func.a,
      test_func.b,
      test_func.func ) );

    nonuniform_results.push_back( test_convergence(
      &qs_akima,
      "Quintic (AKIMA)",
      test_func.name,
      false,
      test_func.a,
      test_func.b,
      test_func.func ) );

    nonuniform_results.push_back( test_convergence(
      &qs_vanleer,
      "Quintic (VANLEER)",
      test_func.name,
      false,
      test_func.a,
      test_func.b,
      test_func.func ) );

    print_convergence_table_for_test( nonuniform_results, "NON-UNIFORM", test_func, func_idx );

    // Test derivative accuracy for this function (only for selected splines)
    if ( func_idx == SZ - 1 )
    {  // Test only for the last function (exp+trig)
      fmt::print( fg( fmt::color::magenta ), "\n🔍 DERIVATIVE ACCURACY FOR {}:\n", test_func.name );

      // Test derivative accuracy for different spline types
      test_derivative_consistency_for_function( test_func, "Cubic", new Splines::CubicSpline() );
      test_derivative_consistency_for_function( test_func, "Akima", new Splines::AkimaSpline() );
      test_derivative_consistency_for_function( test_func, "VanLeer", new Splines::VanLeerSpline() );
      test_derivative_consistency_for_function( test_func, "Pchip", new Splines::PchipSpline() );
      test_derivative_consistency_for_function(
        test_func,
        "Quintic",
        new Splines::QuinticSpline( Spline_sub_type::CUBIC ) );
      test_derivative_consistency_for_function(
        test_func,
        "Quintic[Akima]",
        new Splines::QuinticSpline( Spline_sub_type::AKIMA ) );
      test_derivative_consistency_for_function(
        test_func,
        "Quintic[VanLeer]",
        new Splines::QuinticSpline( Spline_sub_type::VANLEER ) );
      test_derivative_consistency_for_function(
        test_func,
        "Quintic[Pchi]",
        new Splines::QuinticSpline( Spline_sub_type::PCHIP ) );
    }
  }

  // ========================================================================
  // GLOBAL SUMMARY AND THEORETICAL EXPECTATIONS
  // ========================================================================

  print_header( "GLOBAL SUMMARY AND THEORETICAL EXPECTATIONS" );

  fmt::print(
    fg( fmt::color::cyan ),
    "\n"
    "THEORETICAL CONVERGENCE ORDERS SUMMARY:\n"
    "  • Linear splines:     O(h²) for function, O(h) for derivative\n"
    "  • Constant splines:   O(h)  for function (piecewise constant)\n"
    "  • Cubic splines:      O(h⁴) for function, O(h³) for derivative\n"
    "  • Akima splines:      O(h²) for function, O(h) for derivative\n"
    "  • VanLeer splines:    O(h²) for function, O(h) for derivative\n"
    "  • Pchip splines:      O(h⁴) for function (shape-preserving)\n"
    "  • Quintic splines:    O(h⁶) for function, O(h⁵) for derivative\n"
    "\n"
    "EXACT INTERPOLATION PROPERTIES:\n"
    "  • Linear splines:     Exact for polynomials of degree ≤ 1\n"
    "  • Cubic splines:      Exact for polynomials of degree ≤ 3\n"
    "  • Quintic splines:    Exact for polynomials of degree ≤ 5\n"
    "\n"
    "NON-UNIFORM MESH EFFECTS:\n"
    "  • Non-uniform meshes generally maintain theoretical orders\n"
    "  • May affect constant factors in error bounds\n"
    "  • Shape-preserving splines (Pchip, Akima) may be more sensitive\n" );

  // ========================================================================
  // SPECIAL NOTES ON NON-UNIFORM MESH CONVERGENCE
  // ========================================================================

  print_header( "SPECIAL NOTES ON NON-UNIFORM MESH CONVERGENCE" );

  fmt::print(
    fg( fmt::color::cyan ),
    "\n"
    "IMPROVEMENTS IN NON-UNIFORM MESH CONVERGENCE TESTING:\n"
    "  1. Hierarchical refinement: Start from a base mesh and refine by inserting midpoints\n"
    "  2. Controlled h-ratio: Each refinement step approximately halves the mesh size\n"
    "  3. Accurate mesh size computation: Use average interval length instead of 1/N\n"
    "  4. Chebyshev-like base mesh: Ensures genuine non-uniformity with clustering at boundaries\n"
    "\n"
    "COMPARISON METHODOLOGY:\n"
    "  • Uniform mesh: N = 5, 9, 17, 33, 65, 129, 257 points\n"
    "  • Non-uniform mesh: Same sequence via hierarchical refinement\n"
    "  • Mesh size h computed as average interval length\n"
    "  • Order computed as log(error₁/error₂)/log(h₁/h₂)\n"
    "\n"
    "EXPECTED RESULTS:\n"
    "  • Orders should be similar between uniform and non-uniform meshes\n"
    "  • Non-uniform may show slightly different orders due to point distribution\n"
    "  • Exact interpolation properties should hold for both mesh types\n" );

  // ========================================================================
  // SPECIAL NOTES
  // ========================================================================

  print_header( "SPECIAL NOTES" );

  fmt::print(
    fg( fmt::color::cyan ),
    "\n"
    "NOTES ON TEST FUNCTIONS:\n"
    "  1. Polynomial deg 1: Tests exact interpolation capability\n"
    "  2. Polynomial deg 2: Tests quadratic approximation\n"
    "  3. Polynomial deg 3: Tests cubic exact interpolation\n"
    "  4. Polynomial deg 4: Tests beyond cubic interpolation\n"
    "  5. Exp + sin(2x): Tests convergence for non-polynomial functions\n"
    "\n"
    "INTERPRETING RESULTS:\n"
    "  • 'INF' order means error reached machine precision (exact interpolation)\n"
    "  • Orders near theoretical values indicate proper implementation\n"
    "  • Lower orders may indicate issues or limitations of the method\n"
    "  • Non-uniform mesh results should be comparable to uniform mesh\n"
    "  • Mesh size column shows actual average interval length\n" );

  print_header( "ALL TESTS COMPLETED SUCCESSFULLY" );

  fmt::print( fg( fmt::color::light_green ) | fmt::emphasis::bold, "\n🎉 All spline analyses completed! 🎉\n" );
  fmt::print(
    fg( fmt::color::gray ),
    "Summary:\n"
    "  • Original 6 datasets analyzed and results saved to 'out/' directory\n"
    "  • Convergence analysis completed for 5 different test functions\n"
    "  • 7 spline types tested on both uniform and non-uniform meshes\n"
    "  • Improved non-uniform convergence test with hierarchical refinement\n"
    "  • Tables show error vs mesh refinement and estimated convergence orders\n"
    "  • Derivative accuracy verified for key functions\n"
    "  • Theoretical expectations compared with numerical results\n\n" );

  return 0;
}
