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
using Utils::m_pi;

auto TABLE_COLOR = fg( fmt::color::cyan ) | fmt::emphasis::bold;
auto INFO_COLOR  = fg( fmt::color::light_blue );

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
string get_continuity_class( const DerivativeContinuityResult & result )
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
    string continuity_class = get_continuity_class( result );

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
                                                 { "Quintic (AKIMA)", "C¹" }, { "Quintic (VANLEER)", "C²" } };

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
    "                  Should be C¹ for PCHIP VANLEER and AKIMA subtypes (inherit base spline continuity).\n"
    "\n\n" );

  fmt::print( fg( fmt::color::yellow ) | fmt::emphasis::bold, "⚠ NOTES:\n\n" );

  fmt::print(
    fg( fmt::color::white ),
    "• Green values indicate continuity within tolerance ({:.1e})\n"
    "• Yellow/Orange/Red values indicate increasing levels of discontinuity\n"
    "• The 'Theoretical' column shows expected continuity class for each spline type\n"
    "• Quintic splines with different subtypes inherit continuity properties from their base method\n\n",
    tolerance );

  fmt::print(
    fg( fmt::color::light_green ) | fmt::emphasis::bold,
    "🎉 Derivative continuity analysis completed successfully! 🎉\n" );

  return 0;
}
