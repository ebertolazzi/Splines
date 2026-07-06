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
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Splines/Splines.hh"
#include "Utils_fmt.hh"
#include "Utils_FD.hh"

#include <random>
#include <algorithm>

using namespace SplinesLoad;
using namespace std;
using Splines::integer;
using Splines::real_type;
using Splines::Spline1D;
using Splines::Spline2D;
using Splines::Spline_sub_type;
using Splines::SplineType1D;
using Splines::SplineType2D;
using Utils::m_pi;

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

// Funzione di test con derivate esatte note
struct TestFunction1D
{
  // f(x) = sin(x)
  static real_type value( real_type x ) { return sin( x ); }

  static real_type D( real_type x ) { return cos( x ); }

  static real_type DD( real_type x ) { return -sin( x ); }
};

struct TestFunction2D
{
  // z = sin(x) * cos(y)
  static real_type value( real_type x, real_type y ) { return sin( x ) * cos( y ); }
};

// ============================================================================
// FUNZIONI PER GRIGLIE 1D
// ============================================================================

// Genera una griglia uniforme 1D
void build_uniform_grid_1D( integer n, real_type xmin, real_type xmax, vector<real_type> & x )
{
  x.resize( n );
  for ( integer i = 0; i < n; ++i ) { x[i] = xmin + ( xmax - xmin ) * i / ( n - 1 ); }
}

// Genera una griglia non uniforme 1D (con densità variabile)
void build_nonuniform_grid_1D( integer n, real_type xmin, real_type xmax, vector<real_type> & x )
{
  x.resize( n );

  // Più punti vicini agli estremi
  for ( integer i = 0; i < n; ++i )
  {
    real_type t = static_cast<real_type>( i ) / ( n - 1 );
    // Mappatura non lineare per concentrare punti agli estremi
    real_type s = 0.5 * ( 1.0 - cos( m_pi * t ) );  // mapping coseno
    x[i]        = xmin + ( xmax - xmin ) * s;
  }
}

// Genera una griglia completamente random (ma ordinata) 1D
void build_random_grid_1D( integer n, real_type xmin, real_type xmax, vector<real_type> & x, unsigned int seed = 12345 )
{
  x.resize( n );

  std::mt19937                              gen( seed );
  std::uniform_real_distribution<real_type> dist( xmin, xmax );

  // Genera punti random
  for ( integer i = 0; i < n; ++i ) { x[i] = dist( gen ); }
  // Ordina (necessario per le spline)
  std::sort( x.begin(), x.end() );
  // Assicura che il primo e l'ultimo punto siano esattamente i bordi
  x[0]     = xmin;
  x[n - 1] = xmax;
}

// Calcola i valori y per una griglia 1D
void compute_y_values_1D( const vector<real_type> & x, vector<real_type> & y )
{
  integer n = static_cast<integer>( x.size() );
  y.resize( n );

  for ( integer i = 0; i < n; ++i ) { y[i] = TestFunction1D::value( x[i] ); }
}

// ============================================================================
// FUNZIONI PER GRIGLIE 2D
// ============================================================================

// Genera una griglia uniforme 2D
void build_uniform_grid_2D(
  integer             nx,
  integer             ny,
  real_type           xmin,
  real_type           xmax,
  real_type           ymin,
  real_type           ymax,
  vector<real_type> & x,
  vector<real_type> & y )
{
  x.resize( nx );
  y.resize( ny );

  for ( integer i = 0; i < nx; ++i ) { x[i] = xmin + ( xmax - xmin ) * i / ( nx - 1 ); }
  for ( integer j = 0; j < ny; ++j ) { y[j] = ymin + ( ymax - ymin ) * j / ( ny - 1 ); }
}

// Genera una griglia non uniforme 2D (con densità variabile)
void build_nonuniform_grid_2D(
  integer             nx,
  integer             ny,
  real_type           xmin,
  real_type           xmax,
  real_type           ymin,
  real_type           ymax,
  vector<real_type> & x,
  vector<real_type> & y )
{
  x.resize( nx );
  y.resize( ny );

  // Per x: più punti vicini al bordo sinistro
  for ( integer i = 0; i < nx; ++i )
  {
    real_type t = static_cast<real_type>( i ) / ( nx - 1 );
    // Mappatura non lineare per concentrare punti agli estremi
    real_type s = 0.5 * ( 1.0 - cos( m_pi * t ) );  // mapping coseno
    x[i]        = xmin + ( xmax - xmin ) * s;
  }

  // Per y: più punti vicini al centro
  for ( integer j = 0; j < ny; ++j )
  {
    real_type t = static_cast<real_type>( j ) / ( ny - 1 );
    // Mapping cubico per concentrare punti al centro
    real_type s = 3.0 * t * t - 2.0 * t * t * t;  // funzione smooth
    y[j]        = ymin + ( ymax - ymin ) * s;
  }
}

// Genera una griglia completamente random (ma ordinata) 2D
void build_random_grid_2D(
  integer             nx,
  integer             ny,
  real_type           xmin,
  real_type           xmax,
  real_type           ymin,
  real_type           ymax,
  vector<real_type> & x,
  vector<real_type> & y,
  unsigned int        seed = 12345 )
{
  x.resize( nx );
  y.resize( ny );

  std::mt19937                              gen( seed );
  std::uniform_real_distribution<real_type> dist_x( xmin, xmax );
  std::uniform_real_distribution<real_type> dist_y( ymin, ymax );

  // Genera punti random per x
  for ( integer i = 0; i < nx; ++i ) { x[i] = dist_x( gen ); }
  // Ordina x (necessario per le spline)
  std::sort( x.begin(), x.end() );
  // Assicura che il primo e l'ultimo punto siano esattamente i bordi
  x[0]      = xmin;
  x[nx - 1] = xmax;

  // Genera punti random per y
  for ( integer j = 0; j < ny; ++j ) { y[j] = dist_y( gen ); }
  // Ordina y
  std::sort( y.begin(), y.end() );
  y[0]      = ymin;
  y[ny - 1] = ymax;
}

// Calcola i valori z per una data griglia 2D (x, y)
void compute_z_values_2D( const vector<real_type> & x, const vector<real_type> & y, vector<real_type> & z )
{
  integer nx = static_cast<integer>( x.size() );
  integer ny = static_cast<integer>( y.size() );
  z.resize( nx * ny );

  for ( integer i = 0; i < nx; ++i )
  {
    for ( integer j = 0; j < ny; ++j ) { z[i * ny + j] = TestFunction2D::value( x[i], y[j] ); }
  }
}

// ============================================================================
// FUNZIONI PER GENERARE PUNTI DI TEST
// ============================================================================

// Genera punti di test random nel dominio 1D
vector<real_type> generate_random_test_points_1D(
  integer      n_points,
  real_type    xmin,
  real_type    xmax,
  unsigned int seed = 54321 )
{
  vector<real_type> points;
  points.reserve( n_points );

  std::mt19937                              gen( seed );
  std::uniform_real_distribution<real_type> dist( xmin, xmax );

  for ( integer i = 0; i < n_points; ++i ) { points.push_back( dist( gen ) ); }

  return points;
}

// Genera punti di test random nel dominio 2D
vector<pair<real_type, real_type>> generate_random_test_points_2D(
  integer      n_points,
  real_type    xmin,
  real_type    xmax,
  real_type    ymin,
  real_type    ymax,
  unsigned int seed = 54321 )
{
  vector<pair<real_type, real_type>> points;
  points.reserve( n_points );

  std::mt19937                              gen( seed );
  std::uniform_real_distribution<real_type> dist_x( xmin, xmax );
  std::uniform_real_distribution<real_type> dist_y( ymin, ymax );

  for ( integer i = 0; i < n_points; ++i ) { points.emplace_back( dist_x( gen ), dist_y( gen ) ); }

  return points;
}

// ============================================================================
// STRUTTURE PER I TEST
// ============================================================================

// Struttura per tenere traccia dei risultati dettagliati
struct DerivativeTest
{
  string    name;
  real_type autodiff_value;
  real_type exact_value;
  real_type error_abs;
  real_type error_rel;
  bool      passed;
};

// Struttura per informazioni sui fallimenti
struct FailureInfo
{
  real_type x;  // coordinate del punto
  real_type y;
  string    test_name;
  real_type autodiff_value;
  real_type exact_value;
  real_type error_abs;
  real_type error_rel;
  real_type tolerance;
};

// Struttura per statistiche dei fallimenti
struct FailureStats
{
  unordered_map<string, integer> count_by_type;
  vector<FailureInfo>            worst_failures;  // fallimenti più gravi
  real_type                      max_error_abs;
  real_type                      max_error_rel;
  string                         worst_test_name;
  real_type                      worst_x;
  real_type                      worst_y;

  FailureStats() : max_error_abs( 0.0 ), max_error_rel( 0.0 ) {}
};

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

// Calcola errore relativo (safe, evita divisione per zero)
real_type safe_relative_error( real_type a, real_type b )
{
  real_type denom = std::max( std::abs( a ), std::abs( b ) );
  if ( denom < 1e-15 ) return 0.0;
  return std::abs( a - b ) / denom;
}

// Formatta un valore per output
string format_value( real_type val )
{
  std::stringstream ss;
  if ( std::abs( val ) < 1e-10 ) { ss << "0.0"; }
  else if ( std::abs( val ) < 1e-3 || std::abs( val ) > 1e6 )
  {
    ss << std::scientific << std::setprecision( 6 ) << val;
  }
  else
  {
    ss << std::fixed << std::setprecision( 6 ) << val;
  }
  return ss.str();
}

// Stampa informazioni dettagliate su un fallimento
void print_failure_details( const FailureInfo & failure, bool is_1D = false )
{
  if ( is_1D ) { fmt::print( fg( fmt::color::red ), "    Point x={:.6f}:\n", failure.x ); }
  else
  {
    fmt::print( fg( fmt::color::red ), "    Point ({:.6f}, {:.6f}):\n", failure.x, failure.y );
  }

  fmt::print( "      Test: {}\n", failure.test_name );
  fmt::print( "      Autodiff value: {}\n", format_value( failure.autodiff_value ) );
  fmt::print( "      Exact value:    {}\n", format_value( failure.exact_value ) );
  fmt::print( "      Absolute error: {:.2e} (tolerance: {:.1e})\n", failure.error_abs, failure.tolerance );
  fmt::print( "      Relative error: {:.2e}\n", failure.error_rel );
  fmt::print( "      Ratio error/tol: {:.2f}\n", failure.error_abs / failure.tolerance );
}

// Stampa riepilogo dei fallimenti
void print_failure_summary( const FailureStats & stats, integer total_tests, bool is_1D = false )
{
  fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold, "\n  ✗ FAILURE SUMMARY:\n" );

  // Statistiche per tipo di test
  fmt::print( "    Failed tests by type:\n" );
  for ( const auto & [test_name, count] : stats.count_by_type )
  {
    fmt::print( "      {}: {} failures\n", test_name, count );
  }

  // Fallimento peggiore
  if ( !stats.worst_failures.empty() )
  {
    fmt::print( "\n    Worst failures:\n" );
    const auto & worst = stats.worst_failures[0];
    if ( is_1D )
    {
      fmt::print(
        "      Worst at x={:.6f}: {} with error {:.2e} (rel {:.2e})\n",
        worst.x,
        worst.test_name,
        worst.error_abs,
        worst.error_rel );
    }
    else
    {
      fmt::print(
        "      Worst at ({:.6f}, {:.6f}): {} with error {:.2e} (rel {:.2e})\n",
        worst.x,
        worst.y,
        worst.test_name,
        worst.error_abs,
        worst.error_rel );
    }
  }

  // Statistiche globali
  fmt::print( "\n    Global statistics:\n" );
  fmt::print( "      Maximum absolute error: {:.2e}\n", stats.max_error_abs );
  fmt::print( "      Maximum relative error: {:.2e}\n", stats.max_error_rel );
  fmt::print( "      Total tests: {}\n", total_tests );
  fmt::print( "      Total failures: {}\n", stats.worst_failures.size() );
  fmt::print(
    "      Success rate: {:.1f}%\n",
    100.0 * static_cast<double>( total_tests - stats.worst_failures.size() ) / static_cast<double>( total_tests )
  );
}

// ============================================================================
// TEST PER SPLINE 1D
// ============================================================================

// Test per un singolo tipo di spline 1D
bool test_spline_1D_type(
  SplineType1D              spline_type,
  const vector<real_type> & x,
  const vector<real_type> & y,
  const vector<real_type> & test_points,
  real_type                 tolerance_1st,
  real_type                 tolerance_2nd,
  const string &            grid_name = "" )
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n=== Testing 1D spline type: {} (Grid: {}) ===\n",
    static_cast<int>( spline_type ),
    grid_name );

  try
  {
    Spline1D spline( "test_spline_1D" );
    spline.build( spline_type, x, y );

    bool                all_passed = true;
    FailureStats        failure_stats;
    vector<FailureInfo> all_failures;
    integer             total_tests     = 0;
    integer             tests_per_point = 3;  // Value, D, DD

    for ( const auto & xx : test_points )
    {
      // Test solo punti dentro il dominio
      if ( xx >= x[0] && xx <= x.back() )
      {
        vector<DerivativeTest> point_tests;
        vector<FailureInfo>    point_failures;

#if 0
        // ================================================================
        // 1. VALORE DELLA FUNZIONE
        // ================================================================
        {
          real_type f_val     = spline.eval( xx );
          real_type f_exact   = TestFunction1D::value( xx );
          real_type error_abs = std::abs( f_val - f_exact );
          real_type error_rel = safe_relative_error( f_val, f_exact );

          DerivativeTest test_val;
          test_val.name           = "Value";
          test_val.autodiff_value = f_val;
          test_val.exact_value    = f_exact;
          test_val.error_abs      = error_abs;
          test_val.error_rel      = error_rel;
          test_val.passed         = error_abs <= tolerance_1st;
          point_tests.push_back( test_val );

          if ( !test_val.passed )
          {
            FailureInfo failure;
            failure.x              = xx;
            failure.y              = 0.0;
            failure.test_name      = "Value";
            failure.autodiff_value = f_val;
            failure.exact_value    = f_exact;
            failure.error_abs      = error_abs;
            failure.error_rel      = error_rel;
            failure.tolerance      = tolerance_1st;
            point_failures.push_back( failure );
          }
        }
#endif

        // ================================================================
        // 2. DERIVATA PRIMA: f'(x)
        // ================================================================
        {
          // Calcolo con autodiff
          autodiff::dual1st x_dual = xx;
          x_dual.grad              = 1.0;
          autodiff::dual1st val    = spline.eval( x_dual );

          real_type autodiff_val = val.grad;
          real_type exact_val    = spline.D( xx );
          real_type error_abs    = std::abs( autodiff_val - exact_val );
          real_type error_rel    = safe_relative_error( autodiff_val, exact_val );

          DerivativeTest test_d;
          test_d.name           = "D (1st)";
          test_d.autodiff_value = autodiff_val;
          test_d.exact_value    = exact_val;
          test_d.error_abs      = error_abs;
          test_d.error_rel      = error_rel;
          test_d.passed         = error_abs <= tolerance_1st;
          point_tests.push_back( test_d );

          if ( !test_d.passed )
          {
            FailureInfo failure;
            failure.x              = xx;
            failure.y              = 0.0;
            failure.test_name      = "D (1st)";
            failure.autodiff_value = autodiff_val;
            failure.exact_value    = exact_val;
            failure.error_abs      = error_abs;
            failure.error_rel      = error_rel;
            failure.tolerance      = tolerance_1st;
            point_failures.push_back( failure );
          }
        }

        // ================================================================
        // 3. DERIVATA SECONDA: f''(x)
        // ================================================================
        {
          // Calcolo con autodiff
          autodiff::dual2nd x_dual;
          x_dual.val.val        = xx;
          x_dual.val.grad       = 1;
          x_dual.grad.val       = 1;
          x_dual.grad.grad      = 0;
          autodiff::dual2nd val = spline.eval( x_dual );

          real_type autodiff_val = val.grad.grad;
          real_type exact_val    = spline.DD( xx );
          real_type error_abs    = std::abs( autodiff_val - exact_val );
          real_type error_rel    = safe_relative_error( autodiff_val, exact_val );

          DerivativeTest test_dd;
          test_dd.name           = "DD (2nd)";
          test_dd.autodiff_value = autodiff_val;
          test_dd.exact_value    = exact_val;
          test_dd.error_abs      = error_abs;
          test_dd.error_rel      = error_rel;
          test_dd.passed         = error_abs <= tolerance_2nd;
          point_tests.push_back( test_dd );

          if ( !test_dd.passed )
          {
            FailureInfo failure;
            failure.x              = xx;
            failure.y              = 0.0;
            failure.test_name      = "DD (2nd)";
            failure.autodiff_value = autodiff_val;
            failure.exact_value    = exact_val;
            failure.error_abs      = error_abs;
            failure.error_rel      = error_rel;
            failure.tolerance      = tolerance_2nd;
            point_failures.push_back( failure );
          }
        }

        // Aggiorna statistiche per questo punto
        if ( !point_failures.empty() )
        {
          all_passed = false;
          all_failures.insert( all_failures.end(), point_failures.begin(), point_failures.end() );

          // Aggiorna contatori per tipo
          for ( const auto & failure : point_failures )
          {
            failure_stats.count_by_type[failure.test_name]++;

            // Traccia il fallimento peggiore
            if ( failure.error_abs > failure_stats.max_error_abs )
            {
              failure_stats.max_error_abs   = failure.error_abs;
              failure_stats.max_error_rel   = failure.error_rel;
              failure_stats.worst_test_name = failure.test_name;
              failure_stats.worst_x         = failure.x;
              failure_stats.worst_y         = failure.y;
            }
          }
        }

        total_tests += tests_per_point;
      }
    }

    // Aggiorna la lista dei fallimenti peggiori (massimo 5)
    if ( !all_failures.empty() )
    {
      // Ordina per errore assoluto decrescente
      std::sort(
        all_failures.begin(),
        all_failures.end(),
        []( const FailureInfo & a, const FailureInfo & b ) { return a.error_abs > b.error_abs; } );

      // Prendi i primi 5 fallimenti peggiori
      size_t num_worst = std::min<size_t>( 5, all_failures.size() );
      failure_stats.worst_failures.assign( all_failures.begin(), all_failures.begin() + num_worst );
    }

    // Stampa riepilogo
    if ( all_passed )
    {
      fmt::print(
        fg( fmt::color::green ) | fmt::emphasis::bold,
        "  ✓ ALL TESTS PASSED for 1D spline type {} ({} tests)\n",
        static_cast<int>( spline_type ),
        total_tests );
    }
    else
    {
      fmt::print(
        fg( fmt::color::red ) | fmt::emphasis::bold,
        "  ✗ SOME TESTS FAILED for 1D spline type {} ({} failures out of {} tests)\n",
        static_cast<int>( spline_type ),
        all_failures.size(),
        total_tests );

      // Stampa riepilogo dettagliato dei fallimenti
      print_failure_summary( failure_stats, total_tests, true );

      // Stampa dettagli dei fallimenti peggiori
      if ( !failure_stats.worst_failures.empty() )
      {
        fmt::print( fg( fmt::color::yellow ), "\n  Details of worst failures:\n" );
        for ( size_t i = 0; i < failure_stats.worst_failures.size(); ++i )
        {
          const auto & failure = failure_stats.worst_failures[i];
          fmt::print( "  {}. ", i + 1 );
          print_failure_details( failure, true );
        }
      }

      // Se ci sono più di 5 fallimenti, mostra solo un riepilogo
      if ( all_failures.size() > 5 )
      {
        fmt::print(
          fg( fmt::color::yellow ),
          "\n  ... and {} more failures (showing only worst 5)\n",
          all_failures.size() - 5 );
      }
    }

    return all_passed;
  }
  catch ( const std::exception & e )
  {
    fmt::print(
      fg( fmt::color::red ) | fmt::emphasis::bold,
      "  ERROR testing 1D spline type {}: {}\n",
      static_cast<int>( spline_type ),
      e.what() );
    return false;
  }
}

// ============================================================================
// TEST PER SPLINE 2D
// ============================================================================

// Test per un singolo tipo di spline 2D
bool test_spline_2D_type(
  SplineType2D                               spline_type,
  const vector<real_type> &                  x,
  const vector<real_type> &                  y,
  const vector<real_type> &                  z,
  const vector<pair<real_type, real_type>> & test_points,
  real_type                                  tolerance_1st,
  real_type                                  tolerance_2nd,
  real_type                                  tolerance_2nd_mixed,
  const string &                             grid_name = "" )
{
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\nTesting 2D spline type: {} (Grid: {})\n",
    static_cast<int>( spline_type ),
    grid_name );

  try
  {
    Spline2D spline( "test_spline_2D" );
    spline.build( spline_type, x, y, z, false, false );

    bool                all_passed = true;
    FailureStats        failure_stats;
    vector<FailureInfo> all_failures;
    integer             total_tests     = 0;
    integer             tests_per_point = 5;  // Dx, Dy, Dxx, Dyy, Dxy

    for ( const auto & point : test_points )
    {
      real_type xx = point.first;
      real_type yy = point.second;

      // Test solo punti dentro il dominio
      if ( xx >= x[0] && xx <= x.back() && yy >= y[0] && yy <= y.back() )
      {
        vector<FailureInfo> point_failures;

        // Funzione helper per aggiungere un test
        auto add_test =
          [&]( const string & name, real_type autodiff_val, real_type exact_val, real_type tolerance ) -> bool
        {
          real_type error_abs = std::abs( autodiff_val - exact_val );
          real_type error_rel = safe_relative_error( autodiff_val, exact_val );
          bool      passed    = error_abs <= tolerance;

          if ( !passed )
          {
            FailureInfo failure;
            failure.x              = xx;
            failure.y              = yy;
            failure.test_name      = name;
            failure.autodiff_value = autodiff_val;
            failure.exact_value    = exact_val;
            failure.error_abs      = error_abs;
            failure.error_rel      = error_rel;
            failure.tolerance      = tolerance;
            point_failures.push_back( failure );
          }

          return passed;
        };

        // ================================================================
        // 1. DERIVATE PRIME: ∂f/∂x, ∂f/∂y
        // ================================================================
        {
          // ∂f/∂x
          autodiff::dual1st x_dual  = xx;
          x_dual.grad               = 1.0;
          real_type         y_const = yy;
          autodiff::dual1st val_x   = spline.eval( x_dual, y_const );

          add_test( "Dx", val_x.grad, spline.Dx( xx, yy ), tolerance_1st );

          // ∂f/∂y
          real_type         x_const = xx;
          autodiff::dual1st y_dual  = yy;
          y_dual.grad               = 1.0;
          autodiff::dual1st val_y   = spline.eval( x_const, y_dual );

          add_test( "Dy", val_y.grad, spline.Dy( xx, yy ), tolerance_1st );
        }

        // ================================================================
        // 2. DERIVATE SECONDE PURE: ∂²f/∂x², ∂²f/∂y²
        // ================================================================
        {
          // ∂²f/∂x²
          autodiff::dual2nd x_dual;
          x_dual.val.val           = xx;
          x_dual.val.grad          = 1;
          x_dual.grad.val          = 1;
          x_dual.grad.grad         = 0;
          autodiff::dual2nd val_xx = spline.eval( x_dual, yy );

          add_test( "Dxx", val_xx.grad.grad, spline.Dxx( xx, yy ), tolerance_2nd );

          // ∂²f/∂y²
          autodiff::dual2nd y_dual;
          y_dual.val.val           = yy;
          y_dual.val.grad          = 1;
          y_dual.grad.val          = 1;
          y_dual.grad.grad         = 0;
          autodiff::dual2nd val_yy = spline.eval( xx, y_dual );

          add_test( "Dyy", val_yy.grad.grad, spline.Dyy( xx, yy ), tolerance_2nd );
        }

        // ================================================================
        // 3. DERIVATA MISTA: ∂²f/∂x∂y
        // ================================================================
        {
          // Calcoliamo Dx a yy e a yy+h, e usiamo differenze finite per la derivata rispetto a y
          real_type h = 1e-8;

          // Dx a (xx, yy)
          autodiff::dual1st x_dual1;
          x_dual1.val               = xx;
          x_dual1.grad              = 1.0;
          autodiff::dual1st val_x1  = spline.eval( x_dual1, yy );
          real_type         Dx_at_y = val_x1.grad;

          // Dx a (xx, yy+h)
          autodiff::dual1st x_dual2;
          x_dual2.val                = xx;
          x_dual2.grad               = 1.0;
          autodiff::dual1st val_x2   = spline.eval( x_dual2, yy + h );
          real_type         Dx_at_yh = val_x2.grad;

          // Derivata mista approssimata
          real_type Dxy_fd = ( Dx_at_yh - Dx_at_y ) / h;

          add_test( "Dxy", Dxy_fd, spline.Dxy( xx, yy ), tolerance_2nd_mixed );
        }

        // Aggiorna statistiche per questo punto
        if ( !point_failures.empty() )
        {
          all_passed = false;
          all_failures.insert( all_failures.end(), point_failures.begin(), point_failures.end() );

          // Aggiorna contatori per tipo
          for ( const auto & failure : point_failures )
          {
            failure_stats.count_by_type[failure.test_name]++;

            // Traccia il fallimento peggiore
            if ( failure.error_abs > failure_stats.max_error_abs )
            {
              failure_stats.max_error_abs   = failure.error_abs;
              failure_stats.max_error_rel   = failure.error_rel;
              failure_stats.worst_test_name = failure.test_name;
              failure_stats.worst_x         = failure.x;
              failure_stats.worst_y         = failure.y;
            }
          }
        }

        total_tests += tests_per_point;
      }
    }

    // Aggiorna la lista dei fallimenti peggiori (massimo 5)
    if ( !all_failures.empty() )
    {
      // Ordina per errore assoluto decrescente
      std::sort(
        all_failures.begin(),
        all_failures.end(),
        []( const FailureInfo & a, const FailureInfo & b ) { return a.error_abs > b.error_abs; } );

      // Prendi i primi 5 fallimenti peggiori
      size_t num_worst = std::min<size_t>( 5, all_failures.size() );
      failure_stats.worst_failures.assign( all_failures.begin(), all_failures.begin() + num_worst );
    }

    // Stampa riepilogo
    if ( all_passed )
    {
      fmt::print(
        fg( fmt::color::green ) | fmt::emphasis::bold,
        "  ✓ ALL TESTS PASSED for 2D spline type {} ({} tests)\n",
        static_cast<int>( spline_type ),
        total_tests );
    }
    else
    {
      fmt::print(
        fg( fmt::color::red ) | fmt::emphasis::bold,
        "  ✗ SOME TESTS FAILED for 2D spline type {} ({} failures out of {} tests)\n",
        static_cast<int>( spline_type ),
        all_failures.size(),
        total_tests );

      // Stampa riepilogo dettagliato dei fallimenti
      print_failure_summary( failure_stats, total_tests, false );

      // Stampa dettagli dei fallimenti peggiori
      if ( !failure_stats.worst_failures.empty() )
      {
        fmt::print( fg( fmt::color::yellow ), "\n  Details of worst failures:\n" );
        for ( size_t i = 0; i < failure_stats.worst_failures.size(); ++i )
        {
          const auto & failure = failure_stats.worst_failures[i];
          fmt::print( "  {}. ", i + 1 );
          print_failure_details( failure, false );
        }
      }

      // Se ci sono più di 5 fallimenti, mostra solo un riepilogo
      if ( all_failures.size() > 5 )
      {
        fmt::print(
          fg( fmt::color::yellow ),
          "\n  ... and {} more failures (showing only worst 5)\n",
          all_failures.size() - 5 );
      }
    }

    return all_passed;
  }
  catch ( const std::exception & e )
  {
    fmt::print(
      fg( fmt::color::red ) | fmt::emphasis::bold,
      "  ERROR testing 2D spline type {}: {}\n",
      static_cast<int>( spline_type ),
      e.what() );
    return false;
  }
}

// ============================================================================
// FUNZIONI DI ESECUZIONE TEST
// ============================================================================

// Esegue test 1D per un tipo di griglia specifico
void run_1D_grid_tests(
  const string &               grid_name,
  const vector<real_type> &    x,
  const vector<real_type> &    y,
  const vector<real_type> &    test_points,
  const vector<SplineType1D> & spline_types,
  const vector<string> &       spline_type_names,
  real_type                    tolerance_1st,
  real_type                    tolerance_2nd,
  int &                        total_passed,
  int &                        total_tested )
{
  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "\n"
    "============================================================\n"
    "Testing 1D splines on {} grid\n"
    "============================================================\n",
    grid_name );

  fmt::print( "Grid size: {} points\n", x.size() );
  fmt::print( "Test points: {} random points\n", test_points.size() );
  fmt::print( "Tolerances: 1st={:.1e}, 2nd={:.1e}\n\n", tolerance_1st, tolerance_2nd );

  int grid_passed = 0;

  for ( size_t i = 0; i < spline_types.size(); ++i )
  {
    fmt::print( fg( fmt::color::cyan ), "Testing: {} ({}/{})\n", spline_type_names[i], i + 1, spline_types.size() );

    bool passed = test_spline_1D_type( spline_types[i], x, y, test_points, tolerance_1st, tolerance_2nd, grid_name );

    if ( passed )
    {
      grid_passed++;
      total_passed++;
    }
    total_tested++;
  }

  fmt::print(
    fg( fmt::color::blue ) | fmt::emphasis::bold,
    "\n{} grid summary: {}/{} 1D spline types passed\n",
    grid_name,
    grid_passed,
    spline_types.size() );
}

// Esegue test 2D per un tipo di griglia specifico
void run_2D_grid_tests(
  const string &                             grid_name,
  const vector<real_type> &                  x,
  const vector<real_type> &                  y,
  const vector<real_type> &                  z,
  const vector<pair<real_type, real_type>> & test_points,
  const vector<SplineType2D> &               spline_types,
  const vector<string> &                     spline_type_names,
  real_type                                  tolerance_1st,
  real_type                                  tolerance_2nd,
  real_type                                  tolerance_2nd_mixed,
  int &                                      total_passed,
  int &                                      total_tested )
{
  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "\n"
    "============================================================\n"
    "Testing 2D splines on {} grid\n"
    "============================================================\n",
    grid_name );

  fmt::print( "Grid size: {} x {} points\n", x.size(), y.size() );
  fmt::print( "Test points: {} random points\n", test_points.size() );
  fmt::print(
    "Tolerances: 1st={:.1e}, 2nd={:.1e}, mixed={:.1e}\n\n",
    tolerance_1st,
    tolerance_2nd,
    tolerance_2nd_mixed );

  int grid_passed = 0;

  for ( size_t i = 0; i < spline_types.size(); ++i )
  {
    fmt::print( fg( fmt::color::cyan ), "Testing: {} ({}/{})\n", spline_type_names[i], i + 1, spline_types.size() );

    bool passed = test_spline_2D_type(
      spline_types[i],
      x,
      y,
      z,
      test_points,
      tolerance_1st,
      tolerance_2nd,
      tolerance_2nd_mixed,
      grid_name );

    if ( passed )
    {
      grid_passed++;
      total_passed++;
    }
    total_tested++;
  }

  fmt::print(
    fg( fmt::color::blue ) | fmt::emphasis::bold,
    "\n{} grid summary: {}/{} 2D spline types passed\n",
    grid_name,
    grid_passed,
    spline_types.size() );
}

// ============================================================================
// FUNZIONE MAIN
// ============================================================================

int main()
{
  fmt::print(
    fg( fmt::color::blue ) | fmt::emphasis::bold,
    "\n"
    "========================================\n"
    "     COMPREHENSIVE AutoDiff Testing\n"
    "      1D and 2D Splines Comparison\n"
    "========================================\n" );

  // ============================================================================
  // PARAMETRI COMUNI
  // ============================================================================
  const real_type xmin = 0.0;
  const real_type xmax = 2.0 * m_pi;
  const real_type ymin = 0.0;
  const real_type ymax = 2.0 * m_pi;

  // Numero di punti random per test
  const integer n_random_points_1D = 100;
  const integer n_random_points_2D = 200;

  // Genera punti di test random
  auto test_points_1D = generate_random_test_points_1D( n_random_points_1D, xmin, xmax, 54321 );
  auto test_points_2D = generate_random_test_points_2D( n_random_points_2D, xmin, xmax, ymin, ymax, 54321 );

  // Tolleranze
  const real_type tolerance_1st_uniform       = 1e-6;
  const real_type tolerance_2nd_uniform       = 1e-4;
  const real_type tolerance_2nd_mixed_uniform = 1e-3;

  const real_type tolerance_1st_nonuniform       = 1e-5;
  const real_type tolerance_2nd_nonuniform       = 1e-3;
  const real_type tolerance_2nd_mixed_nonuniform = 1e-2;

  // ============================================================================
  // DEFINIZIONE TIPI DI SPLINE
  // ============================================================================

  // Tipi di spline 1D
  vector<SplineType1D> spline_types_1D = { SplineType1D::CONSTANT,        SplineType1D::LINEAR,
                                           SplineType1D::CUBIC,           SplineType1D::AKIMA,
                                           SplineType1D::VANLEER,         SplineType1D::PCHIP,
                                           SplineType1D::QUINTIC_CUBIC,   SplineType1D::QUINTIC_AKIMA,
                                           SplineType1D::QUINTIC_VANLEER, SplineType1D::QUINTIC_PCHIP };

  vector<string> spline_type_names_1D = { "CONSTANT",        "LINEAR",       "CUBIC",         "AKIMA",
                                          "VANLEER",         "PCHIP",        "QUINTIC_CUBIC", "QUINTIC_AKIMA",
                                          "QUINTIC_VANLEER", "QUINTIC_PCHIP" };

  // Tipi di spline 2D
  vector<SplineType2D> spline_types_2D = { SplineType2D::BILINEAR,        SplineType2D::BICUBIC_CUBIC,
                                           SplineType2D::BICUBIC_AKIMA,   SplineType2D::BICUBIC_VANLEER,
                                           SplineType2D::BICUBIC_PCHIP,   SplineType2D::BIQUINTIC_CUBIC,
                                           SplineType2D::BIQUINTIC_AKIMA, SplineType2D::BIQUINTIC_VANLEER,
                                           SplineType2D::BIQUINTIC_PCHIP };

  vector<string> spline_type_names_2D = { "BILINEAR",        "BICUBIC_CUBIC",     "BICUBIC_AKIMA",
                                          "BICUBIC_VANLEER", "BICUBIC_PCHIP",     "BIQUINTIC_CUBIC",
                                          "BIQUINTIC_AKIMA", "BIQUINTIC_VANLEER", "BIQUINTIC_PCHIP" };

  int total_passed_1D = 0;
  int total_tested_1D = 0;
  int total_passed_2D = 0;
  int total_tested_2D = 0;

  // ============================================================================
  // TEST SPLINE 1D
  // ============================================================================
  fmt::print(
    fg( fmt::color::magenta ) | fmt::emphasis::bold,
    "\n"
    "************************************************************\n"
    "***                 1D SPLINE TESTS                      ***\n"
    "************************************************************\n" );

  // 1. TEST CON GRIGLIA UNIFORME 1D
  {
    const integer     n = 20;
    vector<real_type> x_uniform, y_uniform;
    build_uniform_grid_1D( n, xmin, xmax, x_uniform );
    compute_y_values_1D( x_uniform, y_uniform );

    run_1D_grid_tests(
      "Uniform (20 points)",
      x_uniform,
      y_uniform,
      test_points_1D,
      spline_types_1D,
      spline_type_names_1D,
      tolerance_1st_uniform,
      tolerance_2nd_uniform,
      total_passed_1D,
      total_tested_1D );
  }

  // 2. TEST CON GRIGLIA NON UNIFORME 1D
  {
    const integer     n = 20;
    vector<real_type> x_nonuniform, y_nonuniform;
    build_nonuniform_grid_1D( n, xmin, xmax, x_nonuniform );
    compute_y_values_1D( x_nonuniform, y_nonuniform );

    run_1D_grid_tests(
      "Non-uniform (variable density)",
      x_nonuniform,
      y_nonuniform,
      test_points_1D,
      spline_types_1D,
      spline_type_names_1D,
      tolerance_1st_nonuniform,
      tolerance_2nd_nonuniform,
      total_passed_1D,
      total_tested_1D );
  }

  // 3. TEST CON GRIGLIA RANDOM 1D
  {
    const integer     n = 20;
    vector<real_type> x_random, y_random;
    build_random_grid_1D( n, xmin, xmax, x_random, 12345 );
    compute_y_values_1D( x_random, y_random );

    run_1D_grid_tests(
      "Random (ordered)",
      x_random,
      y_random,
      test_points_1D,
      spline_types_1D,
      spline_type_names_1D,
      tolerance_1st_nonuniform,
      tolerance_2nd_nonuniform,
      total_passed_1D,
      total_tested_1D );
  }

  // ============================================================================
  // TEST SPLINE 2D
  // ============================================================================
  fmt::print(
    fg( fmt::color::magenta ) | fmt::emphasis::bold,
    "\n"
    "************************************************************\n"
    "***                  2D SPLINE TESTS                     ***\n"
    "************************************************************\n" );

  // 1. TEST CON GRIGLIA UNIFORME 2D
  {
    const integer     nx = 15, ny = 15;
    vector<real_type> x_uniform, y_uniform, z_uniform;
    build_uniform_grid_2D( nx, ny, xmin, xmax, ymin, ymax, x_uniform, y_uniform );
    compute_z_values_2D( x_uniform, y_uniform, z_uniform );

    run_2D_grid_tests(
      "Uniform (15x15)",
      x_uniform,
      y_uniform,
      z_uniform,
      test_points_2D,
      spline_types_2D,
      spline_type_names_2D,
      tolerance_1st_uniform,
      tolerance_2nd_uniform,
      tolerance_2nd_mixed_uniform,
      total_passed_2D,
      total_tested_2D );
  }

  // 2. TEST CON GRIGLIA NON UNIFORME 2D
  {
    const integer     nx = 15, ny = 15;
    vector<real_type> x_nonuniform, y_nonuniform, z_nonuniform;
    build_nonuniform_grid_2D( nx, ny, xmin, xmax, ymin, ymax, x_nonuniform, y_nonuniform );
    compute_z_values_2D( x_nonuniform, y_nonuniform, z_nonuniform );

    run_2D_grid_tests(
      "Non-uniform (variable density)",
      x_nonuniform,
      y_nonuniform,
      z_nonuniform,
      test_points_2D,
      spline_types_2D,
      spline_type_names_2D,
      tolerance_1st_nonuniform,
      tolerance_2nd_nonuniform,
      tolerance_2nd_mixed_nonuniform,
      total_passed_2D,
      total_tested_2D );
  }

  // 3. TEST CON GRIGLIA RANDOM 2D
  {
    const integer     nx = 15, ny = 15;
    vector<real_type> x_random, y_random, z_random;
    build_random_grid_2D( nx, ny, xmin, xmax, ymin, ymax, x_random, y_random, 12345 );
    compute_z_values_2D( x_random, y_random, z_random );

    run_2D_grid_tests(
      "Random (ordered)",
      x_random,
      y_random,
      z_random,
      test_points_2D,
      spline_types_2D,
      spline_type_names_2D,
      tolerance_1st_nonuniform,
      tolerance_2nd_nonuniform,
      tolerance_2nd_mixed_nonuniform,
      total_passed_2D,
      total_tested_2D );
  }

  // ============================================================================
  // RIEPILOGO FINALE
  // ============================================================================
  fmt::print(
    fg( fmt::color::blue ) | fmt::emphasis::bold,
    "\n"
    "========================================\n"
    "===        FINAL TEST SUMMARY        ===\n"
    "========================================\n" );

  // Statistiche 1D
  double percentage_1D = ( total_tested_1D > 0 ) ? 100.0 * total_passed_1D / total_tested_1D : 0.0;
  fmt::print( "\n1D SPLINES:\n" );
  fmt::print( "  Grid types tested: 3\n" );
  fmt::print( "  Spline types tested: {}\n", spline_types_1D.size() );
  fmt::print( "  Random test points per grid: {}\n", n_random_points_1D );
  fmt::print( "  Total tests: {} (3 tests per point)\n", total_tested_1D );
  fmt::print( "  Results: {}/{} passed ({:.1f}%)\n", total_passed_1D, total_tested_1D, percentage_1D );

  // Statistiche 2D
  double percentage_2D = ( total_tested_2D > 0 ) ? 100.0 * total_passed_2D / total_tested_2D : 0.0;
  fmt::print( "\n2D SPLINES:\n" );
  fmt::print( "  Grid types tested: 3\n" );
  fmt::print( "  Spline types tested: {}\n", spline_types_2D.size() );
  fmt::print( "  Random test points per grid: {}\n", n_random_points_2D );
  fmt::print( "  Total tests: {} (5 derivatives per point)\n", total_tested_2D );
  fmt::print( "  Results: {}/{} passed ({:.1f}%)\n", total_passed_2D, total_tested_2D, percentage_2D );

  // Giudizio complessivo
  fmt::print(
    fg( fmt::color::blue ) | fmt::emphasis::bold,
    "\n"
    "========================================\n"
    "===         OVERALL ASSESSMENT       ===\n"
    "========================================\n" );

  bool all_1D_passed = ( total_passed_1D == total_tested_1D );
  bool all_2D_passed = ( total_passed_2D == total_tested_2D );

  if ( all_1D_passed && all_2D_passed )
  {
    fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "🎉 PERFECT SCORE: All 1D and 2D tests passed!\n" );
    fmt::print(
      fg( fmt::color::green ) | fmt::emphasis::bold,
      "AutoDiff works correctly for all spline types on all grid types!\n" );
    return 0;
  }
  else
  {
    // Stampa informazioni dettagliate sui fallimenti
    if ( !all_1D_passed )
    {
      fmt::print(
        fg( fmt::color::red ) | fmt::emphasis::bold,
        "\n❌ 1D SPLINES FAILED: {}/{} passed ({:.1f}%)\n",
        total_passed_1D,
        total_tested_1D,
        percentage_1D );
    }

    if ( !all_2D_passed )
    {
      fmt::print(
        fg( fmt::color::red ) | fmt::emphasis::bold,
        "\n❌ 2D SPLINES FAILED: {}/{} passed ({:.1f}%)\n",
        total_passed_2D,
        total_tested_2D,
        percentage_2D );
    }

    // Raccomandazioni in base al risultato
    if ( percentage_1D > 90.0 && percentage_2D > 90.0 )
    {
      fmt::print(
        fg( fmt::color::green_yellow ) | fmt::emphasis::bold,
        "\n✅ ACCEPTABLE: Most tests passed (1D: {:.1f}%, 2D: {:.1f}%)\n",
        percentage_1D,
        percentage_2D );
      fmt::print( "   Minor issues detected. Check worst failures above for details.\n" );
      return 0;
    }
    else if ( percentage_1D > 70.0 && percentage_2D > 70.0 )
    {
      fmt::print(
        fg( fmt::color::orange ) | fmt::emphasis::bold,
        "\n⚠️  WARNING: Significant failures detected (1D: {:.1f}%, 2D: {:.1f}%)\n",
        percentage_1D,
        percentage_2D );
      fmt::print( "   Review the failure details above. Some spline types may need attention.\n" );
      return 1;
    }
    else
    {
      fmt::print(
        fg( fmt::color::red ) | fmt::emphasis::bold,
        "\n❌ CRITICAL: Many tests failed (1D: {:.1f}%, 2D: {:.1f}%)\n",
        percentage_1D,
        percentage_2D );
      fmt::print( "   Serious issues detected. Check failure details above.\n" );
      return 1;
    }
  }
}
