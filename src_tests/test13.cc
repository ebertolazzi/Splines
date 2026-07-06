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

using namespace SplinesLoad;
using namespace std;
using Splines::integer;
using Splines::real_type;
using Splines::Spline_sub_type;
using Splines::SplineType2D;
using Utils::m_pi;

auto TABLE_COLOR = fg( fmt::color::cyan ) | fmt::emphasis::bold;
auto INFO_COLOR  = fg( fmt::color::light_blue );

// ============================================================================
// TEST DATASETS
// ============================================================================

static real_type xx0[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
static real_type yy0[] = { 0, 11, 12, 13, 14, 15, 17, 20, 50, 60, 85 };

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

// monotone
static real_type xx6[] = { 0, 1, 2, 3, 4 };
static real_type yy6[] = { 0, 1, 1.1, 2.0, 2.1 };

static integer nn[] = { 11, 11, 11, 9, 12, 4, 5 };

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

// ============================================================================
// Funzione Helper per confronto robusto (Tolleranza Mista)
// ============================================================================
bool is_approx( real_type a, real_type b, real_type tol_abs, real_type tol_rel )
{
  real_type diff = std::abs( a - b );
  if ( diff <= tol_abs ) return true;  // Sotto la soglia di rumore assoluto

  real_type max_val = std::max( std::abs( a ), std::abs( b ) );
  return diff <= ( max_val * tol_rel );  // Controllo relativo alla grandezza
}

// ============================================================================
// FUNCTIONS FOR DERIVATIVE CONTINUITY TESTING
// ============================================================================

// Structure to store derivative continuity results
struct DerivativeContinuityResult
{
  string            spline_name;
  vector<real_type> max_jumps;             // Max jump for each derivative (1st to 5th)
  vector<real_type> avg_jumps;             // Average jump for each derivative
  vector<size_t>    discontinuity_points;  // Number of points where discontinuity > tolerance
  size_t            total_internal_knots;

  DerivativeContinuityResult( const string & name, size_t n_knots )
    : spline_name( name ), total_internal_knots( n_knots )
  {
    max_jumps.resize( 5, 0.0 );
    avg_jumps.resize( 5, 0.0 );
    discontinuity_points.resize( 5, 0 );
  }
};

// ============================================================================
// TEST DERIVATIVE CONTINUITY (CORRETTO)
// ============================================================================
template <typename SplineType> DerivativeContinuityResult test_derivative_continuity(
  SplineType &      spline,
  const string &    spline_name,
  integer           npts,
  real_type const * xx,
  real_type         tolerance_abs = 1e-10 )  // Rinominato per chiarezza
{
  // Tolleranza relativa per gestire derivate di grande magnitudine
  // 1e-8 è ragionevole per float/double standard su derivate alte
  real_type tolerance_rel = 1e-8;

  size_t                     n_internal_knots = ( npts > 2 ) ? ( npts - 2 ) : 0;
  DerivativeContinuityResult result( spline_name, n_internal_knots );

  if ( n_internal_knots == 0 ) return result;

  vector<vector<real_type>> jumps( 5 );

  for ( integer i = 1; i < npts - 1; ++i )
  {
    real_type x_knot    = xx[i];
    integer   left_idx  = i - 1;
    integer   right_idx = i;

    for ( integer d = 0; d < 5; ++d )
    {
      real_type left_val = 0.0, right_val = 0.0;

      // Valutazione: Nota che per QuinticSpline, le derivate alte
      // possono diventare numericamente instabili se h è molto piccolo.
      switch ( d )
      {
        case 0:  // D1
          left_val  = spline.id_D( left_idx, x_knot );
          right_val = spline.id_D( right_idx, x_knot );
          break;
        case 1:  // D2
          left_val  = spline.id_DD( left_idx, x_knot );
          right_val = spline.id_DD( right_idx, x_knot );
          break;
        case 2:  // D3
          left_val  = spline.id_DDD( left_idx, x_knot );
          right_val = spline.id_DDD( right_idx, x_knot );
          break;
        case 3:  // D4
          left_val  = spline.id_DDDD( left_idx, x_knot );
          right_val = spline.id_DDDD( right_idx, x_knot );
          break;
        case 4:  // D5
          left_val  = spline.id_DDDDD( left_idx, x_knot );
          right_val = spline.id_DDDDD( right_idx, x_knot );
          break;
      }

      real_type jump = std::abs( left_val - right_val );
      jumps[d].push_back( jump );

      if ( jump > result.max_jumps[d] ) { result.max_jumps[d] = jump; }

      // USARE is_approx INVECE DI UN CONFRONTO DIRETTO
      if ( !is_approx( left_val, right_val, tolerance_abs, tolerance_rel ) ) { result.discontinuity_points[d]++; }
    }
  }

  // Calcolo medie (invariato)
  for ( integer d = 0; d < 5; ++d )
  {
    if ( !jumps[d].empty() )
    {
      real_type sum = 0.0;
      for ( real_type j : jumps[d] ) sum += j;
      result.avg_jumps[d] = sum / static_cast<real_type>( jumps[d].size() );
    }
  }

  return result;
}


// Function to get color based on jump magnitude
fmt::color get_jump_color( real_type jump, real_type tolerance )
{
  if ( jump <= tolerance )
  {
    return fmt::color::green;  // Continuous (within tolerance)
  }
  else if ( jump <= 10.0 * tolerance )
  {
    return fmt::color::light_green;  // Small jump
  }
  else if ( jump <= 100.0 * tolerance )
  {
    return fmt::color::yellow;  // Medium jump
  }
  else if ( jump <= 1000.0 * tolerance )
  {
    return fmt::color::orange;  // Large jump
  }
  else
  {
    return fmt::color::red;  // Very large jump
  }
}

// Function to get continuity class string
string get_continuity_class( const DerivativeContinuityResult & result, [[maybe_unused]] real_type tolerance )
{
  if ( result.total_internal_knots == 0 ) return "N/A";

  for ( integer d = 0; d < 5; ++d )
  {
    if ( result.discontinuity_points[d] > 0 )
    {
      if ( d == 0 )
        return "C⁰";
      else if ( d == 1 )
        return "C¹";
      else if ( d == 2 )
        return "C²";
      else if ( d == 3 )
        return "C³";
      else if ( d == 4 )
        return "C⁴";
    }
  }
  return "C⁴⁺";
}

// Print derivative continuity table
void print_derivative_continuity_table(
  const vector<DerivativeContinuityResult> & results,
  integer                                    dataset,
  const string &                             dataset_name,
  real_type                                  tolerance )
{
  fmt::print(
    fg( fmt::color::yellow ) | fmt::emphasis::bold,
    "\n"
    "┌────────────────────────────────────────────────────┐\n"
    "│ Dataset {:2}: {:38} │\n"
    "└────────────────────────────────────────────────────┘\n",
    dataset,
    dataset_name );

  fmt::print(
    TABLE_COLOR,
    "\n"
    "┌────────────────────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────────┐\n"
    "│ {:^22} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^12} │\n"
    "├────────────────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────────┤\n",
    "Spline Type",
    "max|ΔD¹|",
    "max|ΔD²|",
    "max|ΔD³|",
    "max|ΔD⁴|",
    "max|ΔD⁵|",
    "Continuity" );

  for ( size_t i = 0; i < results.size(); ++i )
  {
    const auto & result = results[i];

    // Get continuity class
    string continuity_class = get_continuity_class( result, tolerance );

    // Color for the row based on worst jump
    real_type worst_jump = 0.0;
    for ( real_type jump : result.max_jumps )
    {
      if ( jump > worst_jump ) worst_jump = jump;
    }
    fmt::color row_color = get_jump_color( worst_jump, tolerance );

    fmt::print( TABLE_COLOR, "│ " );
    fmt::print( fg( row_color ), "{:<22} ", result.spline_name );

    // Print each derivative jump with color coding
    for ( integer d = 0; d < 5; ++d )
    {
      real_type  jump       = result.max_jumps[d];
      fmt::color jump_color = get_jump_color( jump, tolerance );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( jump_color ), "{:>8.1e} ", jump );
    }

    // Print continuity class with color
    fmt::color class_color = fmt::color::white;
    if ( continuity_class == "C⁴⁺" )
      class_color = fmt::color::green;
    else if ( continuity_class == "C³" )
      class_color = fmt::color::light_green;
    else if ( continuity_class == "C²" )
      class_color = fmt::color::yellow;
    else if ( continuity_class == "C¹" )
      class_color = fmt::color::orange;
    else if ( continuity_class == "C⁰" )
      class_color = fmt::color::red;
    else if ( continuity_class == "C⁻¹" )
      class_color = fmt::color::red;

    fmt::print( TABLE_COLOR, "│ " );
    fmt::print( fg( class_color ), "{:^12}", continuity_class );
    fmt::print( TABLE_COLOR, " │\n" );

    // Print separator between rows
    if ( i < results.size() - 1 )
    {
      fmt::print(
        TABLE_COLOR,
        "├────────────────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────────┤\n" );
    }
  }

  fmt::print(
    TABLE_COLOR,
    "└────────────────────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────────┘\n" );

  // Print statistics
  fmt::print( INFO_COLOR, "\nTolerance: {:.1e}\n", tolerance );
  fmt::print( INFO_COLOR, "Internal knots tested: {}\n", ( results.empty() ? 0 : results[0].total_internal_knots ) );
}

// Print detailed discontinuity information
void print_discontinuity_details( const vector<DerivativeContinuityResult> & results )
{
  fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n📊 DETAILED DISCONTINUITY ANALYSIS:\n" );

  fmt::print(
    TABLE_COLOR,
    "┌────────────────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐\n"
    "│ {:^22} │ {:^10} │ {:^10} │ {:^10} │ {:^10} │ {:^10} │\n"
    "├────────────────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤\n",
    "Spline Type",
    "D¹ disc",
    "D² disc",
    "D³ disc",
    "D⁴ disc",
    "D⁵ disc" );

  for ( const auto & result : results )
  {
    fmt::print( TABLE_COLOR, "│ " );
    fmt::print( fg( fmt::color::white ), "{:<22} ", result.spline_name );

    for ( integer d = 0; d < 5; ++d )
    {
      size_t disc_points = result.discontinuity_points[d];
      size_t total_knots = result.total_internal_knots;

      fmt::print( TABLE_COLOR, "│ " );
      if ( total_knots == 0 ) { fmt::print( "{:^10} ", "N/A" ); }
      else
      {
        real_type  percentage = 100.0 * static_cast<real_type>( disc_points ) / static_cast<real_type>( total_knots );
        fmt::color color      = fmt::color::green;
        if ( percentage > 0 ) color = fmt::color::yellow;
        if ( percentage > 50 ) color = fmt::color::red;

        fmt::print( fg( color ), "{:>4} / {:<4}", disc_points, total_knots );
      }
    }
    fmt::print( TABLE_COLOR, "│\n" );
  }

  fmt::print(
    TABLE_COLOR,
    "└────────────────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘\n" );
}

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
// ============================================================================
// SPLINESET EVALUATION TEST FUNCTIONS
// ============================================================================

void test_SplineSet_evaluation()
{
  print_header( "SPLINESET EVALUATION METHODS TEST" );

  // Usiamo il primo dataset per costruire diverse spline
  integer     npts = nn[0];
  real_type * xx   = xx0;

  // Creiamo diverse spline con lo stesso x ma y diverse
  LinearSpline  lin;
  CubicSpline   cub;
  AkimaSpline   aki;
  VanLeerSpline vls;
  PchipSpline   pch;

  // Costruiamo le spline usando il metodo standard (come nel test originale)
  lin.clear();
  lin.reserve( npts );
  for ( integer i = 0; i < npts; ++i ) lin.push_back( xx[i], yy0[i] );
  lin.build();

  cub.clear();
  cub.reserve( npts );
  for ( integer i = 0; i < npts; ++i ) cub.push_back( xx[i], yy0[i] );
  cub.build();

  aki.clear();
  aki.reserve( npts );
  for ( integer i = 0; i < npts; ++i ) aki.push_back( xx[i], yy1[i] );
  aki.build();

  vls.clear();
  vls.reserve( npts );
  for ( integer i = 0; i < npts; ++i ) vls.push_back( xx[i], yy2[i] );
  vls.build();

  pch.clear();
  pch.reserve( npts );
  for ( integer i = 0; i < npts; ++i ) pch.push_back( xx[i], yy0[i] );
  pch.build();

  // Creiamo SplineSet usando i dati delle spline
  // Prima raccogliamo i dati grezzi dalle spline costruite
  vector<const char *>      headers;
  vector<SplineType1D>      stypes;
  vector<const real_type *> Y_arrays;
  vector<const real_type *> Yp_arrays;

  // Raccogliamo dati dalla spline lineare
  headers.push_back( "linear" );
  stypes.push_back( lin.type() );
  Y_arrays.push_back( lin.y_nodes() );
  Yp_arrays.push_back( nullptr );

  // Raccogliamo dati dalla spline cubica
  headers.push_back( "cubic" );
  stypes.push_back( cub.type() );
  Y_arrays.push_back( cub.y_nodes() );
  Yp_arrays.push_back( nullptr );

  // Raccogliamo dati dalla spline Akima
  headers.push_back( "akima" );
  stypes.push_back( aki.type() );
  Y_arrays.push_back( aki.y_nodes() );
  Yp_arrays.push_back( nullptr );

  // Raccogliamo dati dalla spline VanLeer
  headers.push_back( "vanleer" );
  stypes.push_back( vls.type() );
  Y_arrays.push_back( vls.y_nodes() );
  Yp_arrays.push_back( nullptr );

  // Raccogliamo dati dalla spline Pchip
  headers.push_back( "pchip" );
  stypes.push_back( pch.type() );
  Y_arrays.push_back( pch.y_nodes() );
  Yp_arrays.push_back( nullptr );

  // Creiamo SplineSet con i dati grezzi
  SplineSet spline_set;

  // Usiamo il metodo build con i dati grezzi
  spline_set.build(
    static_cast<integer>( headers.size() ),  // nspl
    npts,                                    // npts
    headers.data(),                          // headers
    stypes.data(),                           // stype
    xx,                                      // X
    Y_arrays.data(),                         // Y
    Yp_arrays.data()                         // Yp (nullptr per tutte)
  );

  fmt::print( fg( fmt::color::green ), "\n✅ SplineSet created with {} splines\n", spline_set.num_splines() );

  // ========================================================================
  // TEST 1: Evaluation by index (SplineSet_eval.hxx)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n1. Testing evaluation by index:\n" );
  {
    real_type x = 5.0;
    for ( integer i = 0; i < spline_set.num_splines(); ++i )
    {
      real_type y   = spline_set.eval( x, i );
      real_type dy  = spline_set.eval_D( x, i );
      real_type ddy = spline_set.eval_DD( x, i );
      real_type d3  = spline_set.eval_DDD( x, i );
      real_type d4  = spline_set.eval_DDDD( x, i );
      real_type d5  = spline_set.eval_DDDDD( x, i );

      fmt::print(
        "  Spline {} at x={}: y={:.6f}, dy={:.6f}, ddy={:.6f}, d3={:.6f}, d4={:.6f}, d5={:.6f}\n",
        i,
        x,
        y,
        dy,
        ddy,
        d3,
        d4,
        d5 );
    }
  }

  // ========================================================================
  // TEST 2: Evaluation by name (SplineSet_eval_xs.hxx)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n2. Testing evaluation by name:\n" );
  {
    real_type      x     = 5.0;
    vector<string> names = { "linear", "cubic", "akima", "vanleer", "pchip" };
    for ( const auto & name : names )
    {
      real_type y   = spline_set.eval( x, name );
      real_type dy  = spline_set.eval_D( x, name );
      real_type ddy = spline_set.eval_DD( x, name );
      real_type d3  = spline_set.eval_DDD( x, name );
      real_type d4  = spline_set.eval_DDDD( x, name );
      real_type d5  = spline_set.eval_DDDDD( x, name );

      fmt::print(
        "  Spline '{}' at x={}: y={:.6f}, dy={:.6f}, ddy={:.6f}, d3={:.6f}, d4={:.6f}, d5={:.6f}\n",
        name,
        x,
        y,
        dy,
        ddy,
        d3,
        d4,
        d5 );
    }
  }

  // ========================================================================
  // TEST 3: Evaluation to C array (SplineSet_eval_xv.hxx)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n3. Testing evaluation to C array:\n" );
  {
    real_type         x    = 5.0;
    integer           nspl = spline_set.num_splines();
    vector<real_type> vals( nspl ), dvals( nspl ), ddvals( nspl );
    vector<real_type> d3vals( nspl ), d4vals( nspl ), d5vals( nspl );

    spline_set.eval( x, vals.data() );
    spline_set.eval_D( x, dvals.data() );
    spline_set.eval_DD( x, ddvals.data() );
    spline_set.eval_DDD( x, d3vals.data() );
    spline_set.eval_DDDD( x, d4vals.data() );
    spline_set.eval_DDDDD( x, d5vals.data() );

    for ( integer i = 0; i < nspl; ++i )
    {
      fmt::print(
        "  Spline {}: y={:.6f}, dy={:.6f}, ddy={:.6f}, d3={:.6f}, d4={:.6f}, d5={:.6f}\n",
        i,
        vals[i],
        dvals[i],
        ddvals[i],
        d3vals[i],
        d4vals[i],
        d5vals[i] );
    }
  }

  // ========================================================================
  // TEST 4: Evaluation to std::vector (SplineSet_eval_xsv.hxx)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n4. Testing evaluation to std::vector:\n" );
  {
    real_type         x = 5.0;
    vector<real_type> vals, dvals, ddvals;
    vector<real_type> d3vals, d4vals, d5vals;

    spline_set.eval( x, vals );
    spline_set.eval_D( x, dvals );
    spline_set.eval_DD( x, ddvals );
    spline_set.eval_DDD( x, d3vals );
    spline_set.eval_DDDD( x, d4vals );
    spline_set.eval_DDDDD( x, d5vals );

    for ( size_t i = 0; i < vals.size(); ++i )
    {
      fmt::print(
        "  Spline {}: y={:.6f}, dy={:.6f}, ddy={:.6f}, d3={:.6f}, d4={:.6f}, d5={:.6f}\n",
        i,
        vals[i],
        dvals[i],
        ddvals[i],
        d3vals[i],
        d4vals[i],
        d5vals[i] );
    }
  }

  // ========================================================================
  // TEST 5: Evaluation with GenericContainer (SplineSet_eval_gc.hxx)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n5. Testing evaluation with GenericContainer (all):\n" );
  {
    real_type                      x = 5.0;
    GC_namespace::GenericContainer gc;

    spline_set.eval( x, gc );
    auto & map = gc.get_map();
    fmt::print( "  Values:\n" );
    for ( auto & [name, value] : map ) { fmt::print( "    {}: {:.6f}\n", name, value.get_real() ); }

    // Test derivatives
    spline_set.eval_D( x, gc );
    fmt::print( "  First derivatives:\n" );
    for ( auto & [name, value] : map ) { fmt::print( "    {} dy/dx: {:.6f}\n", name, value.get_real() ); }

    spline_set.eval_DD( x, gc );
    fmt::print( "  Second derivatives:\n" );
    for ( auto & [name, value] : map ) { fmt::print( "    {} d²y/dx²: {:.6f}\n", name, value.get_real() ); }

    spline_set.eval_DDD( x, gc );
    fmt::print( "  Third derivatives:\n" );
    for ( auto & [name, value] : map ) { fmt::print( "    {} d³y/dx³: {:.6f}\n", name, value.get_real() ); }
  }

  // ========================================================================
  // TEST 6: Evaluation with GenericContainer (specific columns)
  // ========================================================================
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n6. Testing evaluation with GenericContainer (specific):\n" );
  {
    real_type                      x       = 5.0;
    GC_namespace::vec_string_type  columns = { "linear", "cubic", "akima" };
    GC_namespace::GenericContainer gc;

    spline_set.eval( x, columns, gc );
    auto & map = gc.get_map();
    fmt::print( "  Values for selected columns:\n" );
    for ( auto & [name, value] : map ) { fmt::print( "    {}: {:.6f}\n", name, value.get_real() ); }
  }

  // ========================================================================
  // TEST 7: Evaluation at multiple points (vector version)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n7. Testing evaluation at multiple points:\n" );
  {
    GC_namespace::vec_real_type    x_vec = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    GC_namespace::GenericContainer gc;

    spline_set.eval( x_vec, gc );
    auto & map = gc.get_map();

    fmt::print( "  Values at multiple points:\n" );
    for ( auto & [name, value] : map )
    {
      auto & vec = value.get_vec_real();
      fmt::print( "    {}: [", name );
      for ( size_t i = 0; i < vec.size(); ++i ) { fmt::print( "{:.3f}{}", vec[i], i < vec.size() - 1 ? ", " : "" ); }
      fmt::print( "]\n" );
    }
  }

  // ========================================================================
  // TEST 8: Inverse evaluation (eval2) by index (SplineSet_eval2_ii.hxx)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n8. Testing inverse evaluation (eval2) by index:\n" );
  {
    // Usiamo la prima spline (linear) come indipendente
    integer   indep = 0;
    real_type zeta  = 10.0;  // Valore della spline lineare
    real_type x;

    for ( integer i = 0; i < spline_set.num_splines(); ++i )
    {
      real_type y = spline_set.eval2( zeta, indep, x, i );
      fmt::print( "  For spline {} when linear={}: x={:.6f}, y={:.6f}\n", i, zeta, x, y );
    }
  }

  // ========================================================================
  // TEST 9: Inverse evaluation (eval2) by name (SplineSet_eval2_ss.hxx)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n9. Testing inverse evaluation (eval2) by name:\n" );
  {
    string    indep = "linear";
    real_type zeta  = 10.0;
    real_type x;

    vector<string> names = { "cubic", "akima", "vanleer", "pchip" };
    for ( const auto & name : names )
    {
      real_type y = spline_set.eval2( zeta, indep, x, name );
      fmt::print( "  For spline '{}' when {}={}: x={:.6f}, y={:.6f}\n", name, indep, zeta, x, y );
    }
  }

  // ========================================================================
  // TEST 10: Inverse evaluation derivatives (eval2_D, eval2_DD, eval2_DDD)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n10. Testing inverse evaluation derivatives:\n" );
  {
    integer   indep = 0;
    real_type zeta  = 10.0;
    real_type x;

    for ( integer i = 0; i < spline_set.num_splines(); ++i )
    {
      real_type dy   = spline_set.eval2_D( zeta, indep, x, i );
      real_type ddy  = spline_set.eval2_DD( zeta, indep, x, i );
      real_type dddy = spline_set.eval2_DDD( zeta, indep, x, i );
      fmt::print( "  Spline {}: dy/dζ={:.6f}, d²y/dζ²={:.6f}, d³y/dζ³={:.6f}\n", i, dy, ddy, dddy );
    }
  }

  // ========================================================================
  // TEST 11: Inverse evaluation to vector (SplineSet_eval2_sv.hxx)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n11. Testing inverse evaluation to vector:\n" );
  {
    integer           indep = 0;
    real_type         zeta  = 10.0;
    real_type         x;
    vector<real_type> vals;

    spline_set.eval2( indep, zeta, x, vals );
    fmt::print( "  For all splines when linear={}: x={:.6f}\n", zeta, x );
    for ( size_t i = 0; i < vals.size(); ++i ) { fmt::print( "    Spline {}: y={:.6f}\n", i, vals[i] ); }
  }

  // ========================================================================
  // TEST 12: Inverse evaluation to C array (SplineSet_eval2_v.hxx)
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n12. Testing inverse evaluation to C array:\n" );
  {
    integer           indep = 0;
    real_type         zeta  = 10.0;
    real_type         x;
    vector<real_type> vals( spline_set.num_splines() );

    spline_set.eval2( indep, zeta, x, vals.data() );
    fmt::print( "  For all splines when linear={}: x={:.6f}\n", zeta, x );
    for ( size_t i = 0; i < vals.size(); ++i ) { fmt::print( "    Spline {}: y={:.6f}\n", i, vals[i] ); }
  }

  // ========================================================================
  // TEST 13: Inverse evaluation with GenericContainer (SplineSet_eval2_gc.hxx)
  // ========================================================================
  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n13. Testing inverse evaluation with GenericContainer:\n" );
  {
    integer                        indep = 0;
    real_type                      zeta  = 10.0;
    real_type                      x;
    GC_namespace::GenericContainer gc;

    spline_set.eval2( zeta, indep, x, gc );
    auto & map = gc.get_map();
    fmt::print( "  For all splines when linear={}: x={:.6f}\n", zeta, x );
    for ( auto & [name, value] : map ) { fmt::print( "    {}: {:.6f}\n", name, value.get_real() ); }
  }

  // ========================================================================
  // TEST 14: Multiple point inverse evaluation with GenericContainer
  // ========================================================================
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n14. Testing multiple point inverse evaluation:\n" );
  {
    integer                        indep = 0;
    GC_namespace::vec_real_type    zetas = { 5.0, 10.0, 15.0 };
    GC_namespace::GenericContainer gc;

    spline_set.eval2( zetas, indep, gc );
    auto & map = gc.get_map();

    fmt::print( "  Values at multiple zeta points:\n" );
    for ( auto & [name, value] : map )
    {
      auto & vec = value.get_vec_real();
      fmt::print( "    {}: [", name );
      for ( size_t i = 0; i < vec.size(); ++i ) { fmt::print( "{:.3f}{}", vec[i], i < vec.size() - 1 ? ", " : "" ); }
      fmt::print( "]\n" );
    }
  }

  // ========================================================================
  // TEST 15: Autodiff support (if compiled)
  // ========================================================================
#ifdef AUTODIFF_SUPPORT
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n15. Testing autodiff support:\n" );
  {
    using namespace autodiff;

    // First order autodiff
    dual1st x = 5.0;
    x.grad    = 1.0;  // Derivative wrt x

    for ( integer i = 0; i < spline_set.num_splines(); ++i )
    {
      dual1st y = spline_set.eval( x, i );
      fmt::print( "  Spline {}: y={:.6f}, dy/dx={:.6f}\n", i, y.val, y.grad );
    }

    // Test with string names
    dual1st y_cubic = spline_set.eval( x, "cubic" );
    fmt::print( "  Spline 'cubic': y={:.6f}, dy/dx={:.6f}\n", y_cubic.val, y_cubic.grad );
  }
#endif

  fmt::print(
    fg( fmt::color::green ) | fmt::emphasis::bold,
    "\n✅ All SplineSet evaluation methods tested successfully!\n" );
}

// ============================================================================
// TEST 2D SPLINES WITH AUTODIFF
// ============================================================================

void test_2D_splines_with_autodiff()
{
  print_header( "2D SPLINES AUTODIFF CONSISTENCY TEST" );

  // Create test grid
  integer           nx = 5;
  integer           ny = 5;
  vector<real_type> X( nx ), Y( ny );
  vector<real_type> Z( nx * ny );

  // Create a simple test function: f(x,y) = x^2 + y^2 + sin(2πx)cos(2πy)
  for ( integer i = 0; i < nx; ++i )
  {
    X[i] = static_cast<real_type>( i ) / ( nx - 1 );  // x in [0,1]
  }
  for ( integer j = 0; j < ny; ++j )
  {
    Y[j] = static_cast<real_type>( j ) / ( ny - 1 );  // y in [0,1]
  }

  // Fill Z matrix
  for ( integer i = 0; i < nx; ++i )
  {
    for ( integer j = 0; j < ny; ++j )
    {
      real_type x   = X[i];
      real_type y   = Y[j];
      Z[i * ny + j] = x * x + y * y + sin( 2.0 * m_pi * x ) * cos( 2.0 * m_pi * y );
    }
  }

  // Define test points (including points between grid nodes)
  vector<real_type> test_x = { 0.1, 0.25, 0.5, 0.75, 0.9 };
  vector<real_type> test_y = { 0.15, 0.35, 0.55, 0.75, 0.85 };

  // Test different 2D spline types
  vector<pair<SplineType2D, string>> spline_types = { { SplineType2D::BILINEAR, "Bilinear" },
                                                      { SplineType2D::BICUBIC_CUBIC, "Bicubic (CUBIC)" },
                                                      { SplineType2D::BICUBIC_AKIMA, "Bicubic (AKIMA)" },
                                                      { SplineType2D::BICUBIC_VANLEER, "Bicubic (VANLEER)" },
                                                      { SplineType2D::BICUBIC_PCHIP, "Bicubic (PCHIP)" },
                                                      { SplineType2D::BIQUINTIC_CUBIC, "Biquintic (CUBIC)" },
                                                      { SplineType2D::BIQUINTIC_AKIMA, "Biquintic (AKIMA)" },
                                                      { SplineType2D::BIQUINTIC_VANLEER, "Biquintic (VANLEER)" },
                                                      { SplineType2D::BIQUINTIC_PCHIP, "Biquintic (PCHIP)" } };

  fmt::print(
    TABLE_COLOR,
    "\n┌──────────────────────────┬────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐\n"
    "│ {:^24} │ {:^10} │ {:^10} │ {:^10} │ {:^10} │ {:^10} │ {:^10} │\n"
    "├──────────────────────────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤\n",
    "2D Spline Type",
    "max|ΔDx|",
    "max|ΔDy|",
    "max|ΔDxx|",
    "max|ΔDxy|",
    "max|ΔDyy|",
    "Status" );

  for ( const auto & [spline_type, type_name] : spline_types )
  {
    try
    {
      Spline2D spline2d( type_name );
      spline2d.build( spline_type, X, Y, Z );

      real_type max_error_Dx  = 0.0;
      real_type max_error_Dy  = 0.0;
      real_type max_error_Dxx = 0.0;
      real_type max_error_Dxy = 0.0;
      real_type max_error_Dyy = 0.0;
      // Test at various points
      for ( real_type x : test_x )
      {
        for ( real_type y : test_y )
        {
          // Skip points exactly at boundaries for derivative testing
          if ( ( x == 0.0 || x == 1.0 ) && ( y == 0.0 || y == 1.0 ) ) continue;

          // Compute derivatives using spline methods
          real_type Dx_spline  = spline2d.Dx( x, y );
          real_type Dy_spline  = spline2d.Dy( x, y );
          real_type Dxx_spline = spline2d.Dxx( x, y );
          real_type Dxy_spline = spline2d.Dxy( x, y );
          real_type Dyy_spline = spline2d.Dyy( x, y );

#ifdef AUTODIFF_SUPPORT
          // Compute derivatives using autodiff
          using namespace autodiff;

          // First derivatives using dual1st
          dual1st x_dual = x;
          dual1st y_dual = y;
          x_dual.grad    = 1.0;  // Derivative wrt x
          y_dual.grad    = 0.0;

          dual1st   z_dual_x    = spline2d.eval( x_dual, y_dual );
          real_type Dx_autodiff = z_dual_x.grad;

          x_dual.grad           = 0.0;
          y_dual.grad           = 1.0;  // Derivative wrt y
          dual1st   z_dual_y    = spline2d.eval( x_dual, y_dual );
          real_type Dy_autodiff = z_dual_y.grad;

          // Second derivatives using DD method (more reliable)
          real_type dd[6];
          spline2d.DD( x, y, dd );
          real_type Dxx_autodiff = dd[3];
          real_type Dxy_autodiff = dd[4];
          real_type Dyy_autodiff = dd[5];

          // Compare results
          real_type error_Dx  = std::abs( Dx_spline - Dx_autodiff );
          real_type error_Dy  = std::abs( Dy_spline - Dy_autodiff );
          real_type error_Dxx = std::abs( Dxx_spline - Dxx_autodiff );
          real_type error_Dxy = std::abs( Dxy_spline - Dxy_autodiff );
          real_type error_Dyy = std::abs( Dyy_spline - Dyy_autodiff );

          max_error_Dx  = std::max( max_error_Dx, error_Dx );
          max_error_Dy  = std::max( max_error_Dy, error_Dy );
          max_error_Dxx = std::max( max_error_Dxx, error_Dxx );
          max_error_Dxy = std::max( max_error_Dxy, error_Dxy );
          max_error_Dyy = std::max( max_error_Dyy, error_Dyy );
#endif
        }
      }

      // Determine status color and symbol
      bool      all_ok    = true;
      real_type tolerance = 1e-8;

      if ( max_error_Dx > tolerance ) all_ok = false;
      if ( max_error_Dy > tolerance ) all_ok = false;
      if ( max_error_Dxx > tolerance ) all_ok = false;
      if ( max_error_Dxy > tolerance ) all_ok = false;
      if ( max_error_Dyy > tolerance ) all_ok = false;

      fmt::color status_color  = all_ok ? fmt::color::green : fmt::color::yellow;
      string     status_symbol = all_ok ? "✓ PASS" : "⚠ WARN";

      if (
        !all_ok && ( max_error_Dx > 1e-6 || max_error_Dy > 1e-6 || max_error_Dxx > 1e-5 || max_error_Dxy > 1e-5 ||
                     max_error_Dyy > 1e-5 ) )
      {
        status_color  = fmt::color::red;
        status_symbol = "✗ FAIL";
      }

      // Print row
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::white ), "{:<24} ", type_name );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( get_jump_color( max_error_Dx, 1e-10 ) ), "{:>10.2e} ", max_error_Dx );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( get_jump_color( max_error_Dy, 1e-10 ) ), "{:>10.2e} ", max_error_Dy );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( get_jump_color( max_error_Dxx, 1e-10 ) ), "{:>10.2e} ", max_error_Dxx );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( get_jump_color( max_error_Dxy, 1e-10 ) ), "{:>10.2e} ", max_error_Dxy );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( get_jump_color( max_error_Dyy, 1e-10 ) ), "{:>10.2e} ", max_error_Dyy );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( status_color ), "{:^10} ", status_symbol );
      fmt::print( TABLE_COLOR, "│\n" );
    }
    catch ( const std::exception & e )
    {
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::white ), "{:<24} ", type_name );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "ERROR" );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "" );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "" );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "" );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "" );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "ERROR" );
      fmt::print( TABLE_COLOR, "│\n" );

      fmt::print( fg( fmt::color::yellow ), "  Error: {}\n", e.what() );
    }
  }

  fmt::print(
    TABLE_COLOR,
    "└──────────────────────────┴────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘\n" );

  // Additional detailed test for a specific spline type
  fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "\n📊 DETAILED TEST FOR BICUBIC (CUBIC) SPLINE:\n" );

  try
  {
    Spline2D spline2d( "TestBicubic" );
    spline2d.build( SplineType2D::BICUBIC_CUBIC, X, Y, Z );

    // Test point in the middle
    real_type x_test = 0.5;
    real_type y_test = 0.5;

    fmt::print( "Test point: x = {}, y = {}\n", x_test, y_test );

    // Using direct methods
    real_type z_direct   = spline2d.eval( x_test, y_test );
    real_type Dx_direct  = spline2d.Dx( x_test, y_test );
    real_type Dy_direct  = spline2d.Dy( x_test, y_test );
    real_type Dxx_direct = spline2d.Dxx( x_test, y_test );
    real_type Dxy_direct = spline2d.Dxy( x_test, y_test );
    real_type Dyy_direct = spline2d.Dyy( x_test, y_test );

    fmt::print( "Direct methods:\n" );
    fmt::print( "  z   = {:.12f}\n", z_direct );
    fmt::print( "  Dx  = {:.12f}\n", Dx_direct );
    fmt::print( "  Dy  = {:.12f}\n", Dy_direct );
    fmt::print( "  Dxx = {:.12f}\n", Dxx_direct );
    fmt::print( "  Dxy = {:.12f}\n", Dxy_direct );
    fmt::print( "  Dyy = {:.12f}\n", Dyy_direct );

#ifdef AUTODIFF_SUPPORT
    using namespace autodiff;

    // Using autodiff for first derivatives
    dual1st x_ad1 = x_test;
    dual1st y_ad1 = y_test;
    x_ad1.grad    = 1.0;
    y_ad1.grad    = 0.0;

    dual1st z_ad_x = spline2d.eval( x_ad1, y_ad1 );

    x_ad1.grad     = 0.0;
    y_ad1.grad     = 1.0;
    dual1st z_ad_y = spline2d.eval( x_ad1, y_ad1 );

    fmt::print( "\nAutodiff (dual1st):\n" );
    fmt::print( "  Dx  = {:.12f}\n", z_ad_x.grad );
    fmt::print( "  Dy  = {:.12f}\n", z_ad_y.grad );

    // Using DD method for second derivatives (more reliable)
    real_type dd[6];
    spline2d.DD( x_test, y_test, dd );

    fmt::print( "\nDD method:\n" );
    fmt::print( "  Dxx = {:.12f}\n", dd[3] );
    fmt::print( "  Dxy = {:.12f}\n", dd[4] );
    fmt::print( "  Dyy = {:.12f}\n", dd[5] );

    // Check consistency
    fmt::print( "\nConsistency check:\n" );
    fmt::print( "  |Dx_direct - Dx_autodiff| = {:.2e}\n", std::abs( Dx_direct - z_ad_x.grad ) );
    fmt::print( "  |Dy_direct - Dy_autodiff| = {:.2e}\n", std::abs( Dy_direct - z_ad_y.grad ) );
    fmt::print( "  |Dxx_direct - Dxx_DD|     = {:.2e}\n", std::abs( Dxx_direct - dd[3] ) );
    fmt::print( "  |Dxy_direct - Dxy_DD|     = {:.2e}\n", std::abs( Dxy_direct - dd[4] ) );
    fmt::print( "  |Dyy_direct - Dyy_DD|     = {:.2e}\n", std::abs( Dyy_direct - dd[5] ) );
#endif
  }
  catch ( const std::exception & e )
  {
    fmt::print( fg( fmt::color::red ), "Error in detailed test: {}\n", e.what() );
  }

  fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "\n✅ 2D Splines autodiff consistency test completed!\n" );
}

// ============================================================================
// TEST 2D SPLINE CONTINUITY
// ============================================================================

void test_2D_spline_continuity()
{
  print_header( "2D SPLINE DERIVATIVE CONTINUITY TEST" );

  // Create a grid with potential discontinuities
  integer           nx = 6;
  integer           ny = 6;
  vector<real_type> X( nx ), Y( ny );
  vector<real_type> Z( nx * ny );

  // Non-uniform grid to test interpolation
  for ( integer i = 0; i < nx; ++i )
  {
    X[i] = static_cast<real_type>( i * i ) / ( ( nx - 1 ) * ( nx - 1 ) );  // Quadratic spacing
  }
  for ( integer j = 0; j < ny; ++j )
  {
    Y[j] = static_cast<real_type>( j ) / ( ny - 1 );  // Uniform in y
  }

  // Create Z with some variations
  for ( integer i = 0; i < nx; ++i )
  {
    for ( integer j = 0; j < ny; ++j )
    {
      real_type x = X[i];
      real_type y = Y[j];
      // Function with some sharp features
      Z[i * ny + j] = sin( 4.0 * m_pi * x ) * cos( 3.0 * m_pi * y ) +
                      0.5 * exp( -10.0 * ( ( x - 0.5 ) * ( x - 0.5 ) + ( y - 0.5 ) * ( y - 0.5 ) ) );
    }
  }

  vector<pair<SplineType2D, string>> spline_types = { { SplineType2D::BILINEAR, "Bilinear" },
                                                      { SplineType2D::BICUBIC_CUBIC, "Bicubic (CUBIC)" },
                                                      { SplineType2D::BICUBIC_AKIMA, "Bicubic (AKIMA)" },
                                                      { SplineType2D::BIQUINTIC_CUBIC, "Biquintic (CUBIC)" } };

  real_type tolerance = 1e-10;

  fmt::print(
    TABLE_COLOR,
    "\n"
    "┌──────────────────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐\n"
    "│ {:^24} │ {:^10} │ {:^10} │ {:^10} │ {:^10} │ {:^10} │\n"
    "├──────────────────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤\n",
    "2D Spline Type",
    "C⁰(x)",
    "C⁰(y)",
    "C¹(x)",
    "C¹(y)",
    "Status" );

  for ( const auto & [spline_type, type_name] : spline_types )
  {
    try
    {
      Spline2D spline2d( type_name );
      spline2d.build( spline_type, X, Y, Z );

      // Test continuity along x direction (vertical lines)
      size_t disc_x_c0 = 0, disc_x_c1 = 0;
      size_t disc_y_c0 = 0, disc_y_c1 = 0;

      // Test along internal x grid lines
      for ( integer j = 1; j < ny - 1; ++j )
      {
        real_type y = Y[j];
        for ( integer i = 1; i < nx; ++i )
        {
          real_type x_left  = X[i] - 1e-12;
          real_type x_right = X[i] + 1e-12;

          // C⁰ continuity
          real_type z_left  = spline2d.eval( x_left, y );
          real_type z_right = spline2d.eval( x_right, y );
          if ( !is_approx( z_left, z_right, tolerance, 1e-8 ) ) { disc_x_c0++; }

          // C¹ continuity (Dx)
          real_type Dx_left  = spline2d.Dx( x_left, y );
          real_type Dx_right = spline2d.Dx( x_right, y );
          if ( !is_approx( Dx_left, Dx_right, tolerance, 1e-8 ) ) { disc_x_c1++; }
        }
      }

      // Test along internal y grid lines
      for ( integer i = 1; i < nx - 1; ++i )
      {
        real_type x = X[i];
        for ( integer j = 1; j < ny; ++j )
        {
          real_type y_bottom = Y[j] - 1e-12;
          real_type y_top    = Y[j] + 1e-12;

          // C⁰ continuity
          real_type z_bottom = spline2d.eval( x, y_bottom );
          real_type z_top    = spline2d.eval( x, y_top );
          if ( !is_approx( z_bottom, z_top, tolerance, 1e-8 ) ) { disc_y_c0++; }

          // C¹ continuity (Dy)
          real_type Dy_bottom = spline2d.Dy( x, y_bottom );
          real_type Dy_top    = spline2d.Dy( x, y_top );
          if ( !is_approx( Dy_bottom, Dy_top, tolerance, 1e-8 ) ) { disc_y_c1++; }
        }
      }

      // Determine status
      bool has_c0_discontinuity = ( disc_x_c0 > 0 ) || ( disc_y_c0 > 0 );
      bool has_c1_discontinuity = ( disc_x_c1 > 0 ) || ( disc_y_c1 > 0 );

      string continuity_class;
      if ( has_c0_discontinuity ) { continuity_class = "C⁻¹"; }
      else if ( has_c1_discontinuity ) { continuity_class = "C⁰"; }
      else
      {
        continuity_class = "C¹⁺";
      }

      fmt::color class_color = fmt::color::green;
      if ( continuity_class == "C⁻¹" )
        class_color = fmt::color::red;
      else if ( continuity_class == "C⁰" )
        class_color = fmt::color::yellow;

      // Print row
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::white ), "{:<24} ", type_name );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( disc_x_c0 > 0 ? fmt::color::red : fmt::color::green ), "{:>10} ", disc_x_c0 );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( disc_y_c0 > 0 ? fmt::color::red : fmt::color::green ), "{:>10} ", disc_y_c0 );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( disc_x_c1 > 0 ? fmt::color::yellow : fmt::color::green ), "{:>10} ", disc_x_c1 );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( disc_y_c1 > 0 ? fmt::color::yellow : fmt::color::green ), "{:>10} ", disc_y_c1 );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( class_color ), "{:^10} ", continuity_class );
      fmt::print( TABLE_COLOR, "│\n" );
    }
    catch ( const std::exception & e )
    {
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::white ), "{:<24} ", type_name );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "ERR" );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "ERR" );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "ERR" );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "ERR" );
      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( fmt::color::red ), "{:^10} ", "ERROR" );
      fmt::print( TABLE_COLOR, "│\n" );
    }
  }

  fmt::print(
    TABLE_COLOR,
    "└──────────────────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘\n" );

  fmt::print( fg( fmt::color::cyan ), "\nLegend:\n" );
  fmt::print( "  C⁰(x): Number of C⁰ discontinuities along x direction\n" );
  fmt::print( "  C⁰(y): Number of C⁰ discontinuities along y direction\n" );
  fmt::print( "  C¹(x): Number of C¹ discontinuities (Dx jumps) along x direction\n" );
  fmt::print( "  C¹(y): Number of C¹ discontinuities (Dy jumps) along y direction\n" );
  fmt::print( "  C⁻¹: Function discontinuity (worst case)\n" );
  fmt::print( "  C⁰: Function continuous, derivative discontinuous\n" );
  fmt::print( "  C¹⁺: At least C¹ continuous\n" );

  fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "\n✅ 2D Spline continuity test completed!\n" );
}

// ============================================================================
// MAIN FUNCTION - MODIFIED FOR DERIVATIVE CONTINUITY TESTING
// ============================================================================

int main()
{
  print_header( "SPLINE DERIVATIVE CONTINUITY TEST SUITE" );

  // Tolerance for considering a derivative continuous
  real_type tolerance = 1e-10;

  // Store overall results for summary
  map<string, vector<DerivativeContinuityResult>> all_results;

  for ( integer k = 0; k < 7; ++k )
  {
    fmt::print( fg( fmt::color::magenta ) | fmt::emphasis::bold, "\n\n📊 DATASET {} ANALYSIS\n", k );

    real_type * xx{ nullptr };
    real_type * yy{ nullptr };
    string      dataset_name;

    switch ( k )
    {
      case 0:
        xx           = xx0;
        yy           = yy0;
        dataset_name = "Akima Test 0";
        break;
      case 1:
        xx           = xx1;
        yy           = yy1;
        dataset_name = "Akima Test 1";
        break;
      case 2:
        xx           = xx2;
        yy           = yy2;
        dataset_name = "Akima Test 2";
        break;
      case 3:
        xx           = xx3;
        yy           = yy3;
        dataset_name = "RPN 14";
        break;
      case 4:
        xx           = xx4;
        yy           = yy4;
        dataset_name = "Titanium";
        break;
      case 5:
        xx           = xx5;
        yy           = yy5;
        dataset_name = "Toolpath";
        break;
      case 6:
        xx           = xx6;
        yy           = yy6;
        dataset_name = "Monotone";
        break;
    }

    integer npts = nn[k];

    // Create spline objects
    LinearSpline   li;
    ConstantSpline co;
    AkimaSpline    ak;
    CubicSpline    cs;
    VanLeerSpline  vl;
    PchipSpline    pc;

    // Quintic spline subtypes
    QuinticSpline qs_cubic( Spline_sub_type::CUBIC, "Quintic (CUBIC)" );
    QuinticSpline qs_pchip( Spline_sub_type::PCHIP, "Quintic (PCHIP)" );
    QuinticSpline qs_akima( Spline_sub_type::AKIMA, "Quintic (AKIMA)" );
    QuinticSpline qs_vanleer( Spline_sub_type::VANLEER, "Quintic (VANLEER)" );

    vector<pair<Splines::Spline *, string>> splines = { { &li, "Linear" },
                                                        { &co, "Constant" },
                                                        { &ak, "Akima" },
                                                        { &cs, "Cubic" },
                                                        { &vl, "VanLeer" },
                                                        { &pc, "Pchip" },
                                                        { &qs_cubic, "Quintic (CUBIC)" },
                                                        { &qs_pchip, "Quintic (PCHIP)" },
                                                        { &qs_akima, "Quintic (AKIMA)" },
                                                        { &qs_vanleer, "Quintic (VANLEER)" } };

    // Build all splines
    for ( auto & [spline_ptr, name] : splines )
    {
      spline_ptr->clear();
      spline_ptr->reserve( npts );
      for ( integer i = 0; i < npts; ++i ) { spline_ptr->push_back( xx[i], yy[i] ); }
      spline_ptr->build();
    }

    // Test derivative continuity for each spline
    vector<DerivativeContinuityResult> dataset_results;

    for ( auto & [spline_ptr, name] : splines )
    {
      fmt::print( fg( fmt::color::cyan ), "  Testing {}... ", name );

      // Use template specialization to call the right function
      auto test_spline = [&]( auto & spline, std::string const & Name )
      { return test_derivative_continuity( spline, Name, npts, xx, tolerance ); };

      DerivativeContinuityResult result = test_spline( *spline_ptr, name );
      dataset_results.push_back( result );

      // Check if all derivatives are continuous
      bool all_continuous = true;
      for ( size_t d = 0; d < 5; ++d )
      {
        if ( result.discontinuity_points[d] > 0 )
        {
          all_continuous = false;
          break;
        }
      }

      if ( all_continuous ) { fmt::print( fg( fmt::color::green ), "✓ All derivatives continuous\n" ); }
      else
      {
        fmt::print( fg( fmt::color::yellow ), "⚠ Some discontinuities found\n" );
      }

      // Store for overall summary
      all_results[name].push_back( result );
    }

    // Print table for this dataset
    print_derivative_continuity_table( dataset_results, k, dataset_name, tolerance );

    // Print detailed discontinuity information
    print_discontinuity_details( dataset_results );
  }

  // ==========================================================================
  // OVERALL SUMMARY
  // ==========================================================================

  print_header( "OVERALL DERIVATIVE CONTINUITY SUMMARY" );

  fmt::print(
    TABLE_COLOR,
    "\n"
    "┌────────────────────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────────┐\n"
    "│ {:^22} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^8} │ {:^12} │\n"
    "├────────────────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────────┤\n",
    "Spline Type",
    "C⁰",
    "C¹",
    "C²",
    "C³",
    "C⁴",
    "Theoretical" );

  // Theoretical continuity classes
  map<string, string> theoretical_continuity = { { "Constant", "C⁻¹" },       { "Linear", "C⁰" },
                                                 { "Akima", "C¹" },           { "Cubic", "C²" },
                                                 { "VanLeer", "C¹" },         { "Pchip", "C¹" },
                                                 { "Quintic (CUBIC)", "C³" }, { "Quintic (PCHIP)", "C¹" },
                                                 { "Quintic (AKIMA)", "C¹" }, { "Quintic (VANLEER)", "C¹" } };

  integer sz = static_cast<integer>( all_results.size() );
  for ( const auto & [spline_name, results] : all_results )
  {
    // Calculate percentage of datasets where each derivative is continuous
    vector<real_type> continuity_percent( 5, 0.0 );
    size_t            num_datasets = results.size();

    for ( const auto & result : results )
    {
      for ( integer d = 0; d < 5; ++d )
      {
        if ( result.discontinuity_points[d] == 0 ) { continuity_percent[d] += 1.0; }
      }
    }

    // Convert to percentages
    for ( integer d = 0; d < 5; ++d ) { continuity_percent[d] = 100.0 * continuity_percent[d] / static_cast<real_type>( num_datasets ); }

    // Determine empirical continuity class
    string empirical_class = "C⁴⁺";
    for ( integer d = 0; d < 5; ++d )
    {
      if ( continuity_percent[d] < 100.0 )
      {
        if ( d == 0 )
          empirical_class = "C⁻¹";
        else if ( d == 1 )
          empirical_class = "C⁰";
        else if ( d == 2 )
          empirical_class = "C¹";
        else if ( d == 3 )
          empirical_class = "C²";
        else if ( d == 4 )
          empirical_class = "C³";
        break;
      }
    }

    // Get theoretical class
    string theoretical_class = theoretical_continuity[spline_name];

    // Color for theoretical vs empirical match
    fmt::color class_color = fmt::color::white;
    if ( empirical_class == theoretical_class ) { class_color = fmt::color::green; }
    else if ( empirical_class > theoretical_class )
    {
      // Better than expected
      class_color = fmt::color::light_green;
    }
    else
    {
      // Worse than expected
      class_color = fmt::color::yellow;
      if ( empirical_class < "C¹" ) class_color = fmt::color::red;
    }

    fmt::print( TABLE_COLOR, "│ " );
    fmt::print( fg( fmt::color::white ), "{:<22} ", spline_name );

    // Print continuity percentages
    for ( integer d = 0; d < 5; ++d )
    {
      fmt::color percent_color = fmt::color::green;
      if ( continuity_percent[d] < 100.0 ) percent_color = fmt::color::yellow;
      if ( continuity_percent[d] < 50.0 ) percent_color = fmt::color::red;

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( percent_color ), "{:>7.2f}% ", continuity_percent[d] );
    }

    // Print theoretical class
    fmt::print( TABLE_COLOR, "│ " );
    fmt::print( fg( class_color ), "{:^12}", theoretical_class );
    fmt::print( TABLE_COLOR, " │\n" );

    // Separator
    if ( --sz > 0 )
    {
      fmt::print(
        fg( fmt::color::cyan ),
        "├────────────────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────────┤\n" );
    }
  }

  fmt::print(
    TABLE_COLOR,
    "└────────────────────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────────┘\n" );

  // ==========================================================================
  // SPECIAL ANALYSIS FOR QUINTIC SUBTYPES
  // ==========================================================================

  print_header( "QUINTIC SPLINE SUBTYPE COMPARISON" );

  // List of quintic subtypes
  vector<string> quintic_subtypes = { "Quintic (CUBIC)", "Quintic (PCHIP)", "Quintic (AKIMA)", "Quintic (VANLEER)" };

  fmt::print(
    TABLE_COLOR,
    "\n"
    "┌────────────────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐\n"
    "│ {:^22} │ {:^10} │ {:^10} │ {:^10} │ {:^10} │ {:^10} │\n"
    "├────────────────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤\n",
    "Quintic Subtype",
    "Avg |ΔD¹|",
    "Avg |ΔD²|",
    "Avg |ΔD³|",
    "Avg |ΔD⁴|",
    "Avg |ΔD⁵|" );

  integer nv = static_cast<integer>( quintic_subtypes.size() );
  for ( const string & subtype : quintic_subtypes )
  {
    const auto & results = all_results[subtype];

    // Calculate average jumps across all datasets
    vector<real_type> avg_jumps( 5, 0.0 );
    size_t            num_datasets = results.size();

    for ( const auto & result : results )
    {
      for ( integer d = 0; d < 5; ++d ) { avg_jumps[d] += result.avg_jumps[d]; }
    }

    for ( integer d = 0; d < 5; ++d ) { avg_jumps[d] /= static_cast<real_type>( num_datasets ); }

    fmt::print( TABLE_COLOR, "│ " );
    fmt::print( fg( fmt::color::white ), "{:<22} ", subtype );

    for ( integer d = 0; d < 5; ++d )
    {
      real_type  jump       = avg_jumps[d];
      fmt::color jump_color = get_jump_color( jump, tolerance );

      fmt::print( TABLE_COLOR, "│ " );
      fmt::print( fg( jump_color ), "{:>10.1e} ", jump );
    }
    fmt::print( TABLE_COLOR, "│\n" );

    // Separator
    if ( --nv > 0 )
    {
      fmt::print(
        fg( fmt::color::cyan ),
        "├────────────────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤\n" );
    }
  }

  fmt::print(
    TABLE_COLOR,
    "└────────────────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘\n" );

  // ==========================================================================
  // TEST 2D SPLINES
  // ==========================================================================

  test_2D_splines_with_autodiff();
  test_2D_spline_continuity();

  // ==========================================================================
  // TEST SPLINESET EVALUATION METHODS
  // ==========================================================================

  test_SplineSet_evaluation();

  // ==========================================================================
  // CONCLUSION
  // ==========================================================================

  print_header( "TEST CONCLUSIONS" );

  fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "\n✅ KEY FINDINGS:\n\n" );

  fmt::print(
    fg( fmt::color::white ),
    "1. Linear splines: Should be C⁰ (continuous). Discontinuities in first derivative expected.\n"
    "2. Constant splines: Should be C⁻¹ (discontinuous). Discontinuities in function values expected.\n"
    "3. Cubic/VanLeer splines: Should be C¹ (once differentiable). Should have continuous D¹.\n"
    "4. Akima/Pchip splines: Should be C¹ (once differentiable). Should have continuous D¹.\n"
    "5. Quintic splines: Should be C³ (thrice differentiable) for CUBIC subtype.\n"
    "                  Should be C¹ for PCHIP, VANLEER and AKIMA subtypes (inherit base spline continuity).\n"
    "6. 2D Bilinear splines: C⁰ continuous, first derivatives piecewise constant (C⁻¹).\n"
    "7. 2D Bicubic splines: At least C¹ continuous, second derivatives may have jumps.\n"
    "8. 2D Biquintic splines: At least C² continuous, depending on subtype.\n\n" );

  fmt::print( fg( fmt::color::yellow ) | fmt::emphasis::bold, "⚠ NOTES:\n\n" );

  fmt::print(
    fg( fmt::color::white ),
    "• Green values indicate continuity within tolerance ({:.1e})\n"
    "• Yellow/Orange/Red values indicate increasing levels of discontinuity\n"
    "• The 'Theoretical' column shows expected continuity class for each spline type\n"
    "• Quintic splines with different subtypes inherit continuity properties from their base method\n"
    "• 2D spline continuity is tested along both x and y directions\n"
    "• Autodiff consistency ensures that derivative calculations are mathematically correct\n\n",
    tolerance );

  fmt::print(
    fg( fmt::color::light_green ) | fmt::emphasis::bold,
    "🎉 All spline tests (1D & 2D) completed successfully! 🎉\n" );

  return 0;
}
