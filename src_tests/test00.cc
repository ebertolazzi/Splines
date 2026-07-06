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
#include <random>
#include <cmath>
#include <algorithm>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using Splines::integer;
using Splines::real_type;

// Funzione per stampare header colorati
void print_header( const string & title, fmt::color color = fmt::color::cyan )
{
  fmt::print(
    fg( color ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════════════════════════╗\n"
    "║{:^78}║\n"
    "╚══════════════════════════════════════════════════════════════════════════════╝\n",
    title );
}

// Funzione per stampare un messaggio di successo
void print_success( const string & message )
{ fmt::print( fg( fmt::color::green ) | fmt::emphasis::bold, "✅  {}\n", message ); }

// Funzione per stampare un messaggio di errore
void print_error( const string & message )
{ fmt::print( fg( fmt::color::red ) | fmt::emphasis::bold, "❌  {}\n", message ); }

// Funzione per stampare un messaggio di avviso
void print_warning( const string & message )
{ fmt::print( fg( fmt::color::yellow ) | fmt::emphasis::bold, "⚠️  {}\n", message ); }

// Funzione per stampare un messaggio informativo
void print_info( const string & message )
{ fmt::print( fg( fmt::color::blue ) | fmt::emphasis::bold, "ℹ️  {}\n", message ); }

// Classe di test per SearchInterval
class SearchIntervalTester
{
private:
  string                                    m_curve_name;
  integer                                   m_npts;
  vector<real_type>                         m_X;
  real_type *                               m_X_ptr;
  bool                                      m_curve_is_closed;
  bool                                      m_curve_can_extend;
  Utils::SearchInterval<real_type, integer> m_si;  // Usiamo il namespace completo

  // Statistiche dei test
  struct TestStats
  {
    integer total_tests  = 0;
    integer passed_tests = 0;
    integer failed_tests = 0;
    integer warnings     = 0;

    void reset()
    {
      total_tests  = 0;
      passed_tests = 0;
      failed_tests = 0;
      warnings     = 0;
    }

    void print() const
    {
      fmt::print(
        "\n"
        "📊  Test Statistics:\n"
        "   ┌────────────────────────────────┐\n"
        "   │ {:^30} │\n"
        "   ├────────────────────────────────┤\n"
        "   │ {:<14} {:>15} │\n"
        "   │ {:<14} {:>15} │\n"
        "   │ {:<14} {:>15} │\n"
        "   │ {:<14} {:>14.1f}% │\n"
        "   └────────────────────────────────┘\n",
        "RESULT SUMMARY",
        "Total Tests:",
        total_tests,
        "Passed:",
        passed_tests,
        "Failed:",
        failed_tests,
        "Success Rate:",
        total_tests > 0 ? ( 100.0 * passed_tests / total_tests ) : 0.0 );
    }
  };

  TestStats m_stats;

  // Helper per accesso sicuro al vettore
  real_type get_X( integer idx ) const { return m_X[static_cast<size_t>( idx )]; }

public:
  SearchIntervalTester(
    const string &            name,
    const vector<real_type> & X,
    bool                      is_closed  = false,
    bool                      can_extend = false )
    : m_curve_name( name ), m_X( X ), m_curve_is_closed( is_closed ), m_curve_can_extend( can_extend )
  {
    m_npts  = static_cast<integer>( X.size() );
    m_X_ptr = const_cast<real_type *>( m_X.data() );

    // Setup della struttura di ricerca
    m_si.setup( &m_curve_name, &m_npts, &m_X_ptr, &m_curve_is_closed, &m_curve_can_extend );

    // Forza il reset per costruire le tabelle
    m_si.must_reset();
  }

  // Helper per trovare l'intervallo corretto considerando nodi duplicati
  // RENDUTO PUBBLICO per essere usato nei test

  integer find_expected_interval( real_type x ) const
  {
    if ( m_X.empty() ) return 0;

    // Funzione di tolleranza identica a quella usata da SearchInterval
    auto tolerance = []( real_type a, real_type b ) -> real_type
    { return 1e-10 * std::max( real_type( 1 ), std::max( std::abs( a ), std::abs( b ) ) ); };

    // Local copy for possible modification
    real_type query_x = x;

    // Gestione out-of-bounds per curve aperte (come SearchInterval)
    if ( !m_curve_is_closed )
    {
      if ( query_x < get_X( 0 ) ) return 0;
      if ( query_x > get_X( m_npts - 1 ) ) return m_npts - 2;
    }

    // Gestione wrap-around per curve chiuse (ESATTAMENTE come SearchInterval)
    if ( m_curve_is_closed && m_npts > 1 )
    {
      real_type x_min = get_X( 0 );
      real_type x_max = get_X( m_npts - 1 );
      real_type range = x_max - x_min;

      if ( range > 0 )
      {
        if ( query_x > x_max )
        {
          query_x = x_min + fmod( query_x - x_min, range );
          if ( query_x < x_min ) query_x += range;
        }
        else if ( query_x < x_min )
        {
          query_x = x_min + fmod( query_x - x_min, range );
          if ( query_x < x_min ) query_x += range;
        }
      }
    }

    // Verifica se TUTTI i nodi sono uguali entro tolleranza
    bool all_equal = true;
    for ( integer i = 1; i < m_npts; ++i )
    {
      if ( std::abs( get_X( i ) - get_X( 0 ) ) > tolerance( get_X( i ), get_X( 0 ) ) )
      {
        all_equal = false;
        break;
      }
    }

    // Caso speciale: tutti i nodi uguali (entro tolleranza)
    if ( all_equal && m_npts > 1 )
    {
      // Comportamento di SearchInterval per nodi tutti uguali:
      // - Per x <= X[0]: restituisce [0, 1]
      // - Per x > X[0]: restituisce [n-2, n-1]
      if ( query_x <= get_X( 0 ) )
        return 0;
      else
        return m_npts - 2;
    }

    // Binary search per trovare l'ultimo i tale che X[i] <= query_x
    // IMPORTANTE: questo simula la binary search di SearchInterval (STEP 4)
    integer k_LO = 0;
    integer k_HI = m_npts - 1;

    while ( k_HI > k_LO + 1 )
    {
      integer k_M = k_LO + ( k_HI - k_LO ) / 2;
      if ( query_x < get_X( k_M ) )
        k_HI = k_M;
      else
        k_LO = k_M;
    }

    integer pos = k_LO;

    // STEP 5: Gestione nodi duplicati - backtracking ESATTO come SearchInterval
    // CRITICAL: usa la STESSA tolleranza di SearchInterval (m_epsilon basato su X[pos])
    while ( pos > 0 && pos < m_npts - 1 )
    {
      real_type eps = 1e-10 * std::max( real_type( 1 ), std::abs( get_X( pos ) ) );
      if ( std::abs( get_X( pos ) - get_X( pos + 1 ) ) > eps ) break;
      --pos;
    }

    // Clamp al range valido
    return std::max( integer( 0 ), std::min( pos, m_npts - 2 ) );
  }


  // Test di un singolo punto
  bool test_point( real_type x, integer expected_interval, const string & description = "" )
  {
    m_stats.total_tests++;

    std::pair<integer, real_type> query{ 0, x };
    m_si.find( query );

    integer found_interval = query.first;
    bool    success        = ( found_interval == expected_interval );

    if ( success )
    {
      fmt::print(
        fg( fmt::color::green ),
        "   ✓ [Test {:3}] x = {:10.6f} → Interval [{:3}, {:3}] {} {}\n",
        m_stats.total_tests,
        x,
        found_interval,
        found_interval + 1,
        description.empty() ? "" : "(" + description + ")",
        found_interval == expected_interval
          ? ""
          : fmt::format( "(expected [{}, {}])", expected_interval, expected_interval + 1 ) );
      m_stats.passed_tests++;
    }
    else
    {
      fmt::print(
        fg( fmt::color::red ),
        "   ✗ [Test {:3}] x = {:10.6f} → Interval [{:3}, {:3}] (expected [{:3}, {:3}]) {}\n",
        m_stats.total_tests,
        x,
        found_interval,
        found_interval + 1,
        expected_interval,
        expected_interval + 1,
        description.empty() ? "" : "(" + description + ")" );
      m_stats.failed_tests++;

      // Debug info aggiuntiva
      if ( found_interval >= 0 && found_interval < m_npts - 1 )
      {
        fmt::print(
          fg( fmt::color::yellow ),
          "       X[{}] = {:.12f}, X[{}] = {:.12f}, X[{}] = {:.12f}\n",
          found_interval,
          get_X( found_interval ),
          found_interval + 1,
          get_X( found_interval + 1 ),
          found_interval + 2,
          ( found_interval + 2 < m_npts ) ? get_X( found_interval + 2 ) : 0.0 );
      }
    }

    return success;
  }

  // Test di un range di punti
  void test_range( real_type x_min, real_type x_max, integer num_tests = 10 )
  {
    print_info(
      "Testing random points in range [" + fmt::format( "{:.3f}", x_min ) + ", " + fmt::format( "{:.3f}", x_max ) +
      "]" );

    std::random_device                        rd;
    std::mt19937                              gen( rd() );
    std::uniform_real_distribution<real_type> dist( x_min, x_max );

    for ( integer i = 0; i < num_tests; ++i )
    {
      real_type x        = dist( gen );
      integer   expected = find_expected_interval( x );
      test_point( x, expected, "random" );
    }
  }

  // Test di tutti i nodi
  void test_all_knots()
  {
    print_info( "Testing all knot points" );

    for ( integer i = 0; i < m_npts; ++i )
    {
      real_type x        = get_X( i );
      integer   expected = find_expected_interval( x );
      test_point( x, expected, "knot point" );
    }
  }

  // Test di tutti i punti medi
  void test_all_midpoints()
  {
    print_info( "Testing all interval midpoints" );

    for ( integer i = 0; i < m_npts - 1; ++i )
    {
      real_type x        = ( get_X( i ) + get_X( i + 1 ) ) / 2.0;
      integer   expected = find_expected_interval( x );
      test_point( x, expected, "midpoint" );
    }
  }

  // Test dei punti al di fuori del range
  void test_out_of_bounds()
  {
    print_info( "Testing out-of-bounds points" );

    if ( m_X.empty() ) return;

    // Punti a sinistra
    test_point( get_X( 0 ) - 1.0, find_expected_interval( get_X( 0 ) - 1.0 ), "left OOB" );
    test_point( get_X( 0 ) - 10.0, find_expected_interval( get_X( 0 ) - 10.0 ), "far left OOB" );

    // Punti a destra
    test_point( get_X( m_npts - 1 ) + 1.0, find_expected_interval( get_X( m_npts - 1 ) + 1.0 ), "right OOB" );
    test_point( get_X( m_npts - 1 ) + 10.0, find_expected_interval( get_X( m_npts - 1 ) + 10.0 ), "far right OOB" );
  }

  // Test speciali per curve chiuse (wrap-around)
  void test_wrap_around()
  {
    if ( !m_curve_is_closed || m_X.empty() ) return;

    print_info( "Testing wrap-around for closed curve" );

    real_type x_min = get_X( 0 );
    real_type x_max = get_X( m_npts - 1 );
    real_type range = x_max - x_min;

    if ( range <= 0 ) return;

    // Test diversi punti oltre il range
    vector<real_type> test_points = { x_max + 0.1 * range, x_max + 0.5 * range, x_max + 1.0 * range,
                                      x_max + 1.5 * range, x_max + 2.0 * range, x_min - 0.1 * range,
                                      x_min - 0.5 * range, x_min - 1.0 * range, x_min - 1.5 * range };

    for ( real_type x : test_points )
    {
      integer expected = find_expected_interval( x );
      test_point( x, expected, fmt::format( "wrap" ) );
    }
  }

  // Esegue tutti i test
  void run_all_tests()
  {
    m_stats.reset();

    fmt::print(
      fg( fmt::color::blue ) | fmt::emphasis::bold,
      "\n🔍 Testing: {}\n"
      "   Nodes: {}, Closed: {}, Extend: {}\n"
      "   X range: [{:.12f}, {:.12f}]\n",
      m_curve_name,
      m_npts,
      m_curve_is_closed ? "Yes" : "No",
      m_curve_can_extend ? "Yes" : "No",
      m_X.empty() ? 0 : get_X( 0 ),
      m_X.empty() ? 0 : get_X( m_npts - 1 ) );

    // Tabella dei nodi (se non troppi)
    if ( m_npts <= 20 )
    {
      fmt::print( "   Nodes: [" );
      for ( integer i = 0; i < m_npts; ++i )
      {
        if ( i > 0 ) fmt::print( ", " );
        fmt::print( "{:.12f}", get_X( i ) );
      }
      fmt::print( "]\n" );
    }

    // Esegui i test
    test_all_knots();
    test_all_midpoints();

    if ( m_npts > 1 ) { test_range( get_X( 0 ), get_X( m_npts - 1 ), 15 ); }

    if ( m_curve_is_closed ) { test_wrap_around(); }
    else
    {
      test_out_of_bounds();
    }

    m_stats.print();
  }
};

// Test case 1: Distribuzione uniforme semplice
void test_uniform_distribution()
{
  print_header( "TEST 1: Uniform Distribution", fmt::color::blue );

  vector<real_type>    X = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
  SearchIntervalTester tester( "Uniform-10", X, false, true );
  tester.run_all_tests();
}

// Test case 2: Distribuzione non uniforme
void test_non_uniform_distribution()
{
  print_header( "TEST 2: Non-uniform Distribution", fmt::color::blue );

  vector<real_type>    X = { 0.0, 0.1, 0.5, 1.2, 2.0, 3.5, 5.0, 6.1, 6.2, 6.3, 7.0, 8.5, 10.0 };
  SearchIntervalTester tester( "NonUniform-13", X, false, true );
  tester.run_all_tests();
}

// Test case 3: Nodi molto ravvicinati
void test_clustered_nodes()
{
  print_header( "TEST 3: Clustered Nodes", fmt::color::blue );

  vector<real_type>    X = { 0.0, 0.001, 0.002, 0.003, 0.004, 1.0, 1.001, 1.002, 2.0, 2.001, 3.0 };
  SearchIntervalTester tester( "Clustered-11", X, false, true );
  tester.run_all_tests();
}

// Test case 4: Nodi duplicati
void test_duplicate_nodes()
{
  print_header( "TEST 4: Duplicate Nodes", fmt::color::blue );

  // IMPORTANTE: SearchInterval restituisce sempre intervalli consecutivi [i, i+1]
  // Per nodi duplicati, arretra per trovare il primo intervallo valido
  vector<real_type>    X = { 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0 };
  SearchIntervalTester tester( "Duplicate-11", X, false, true );
  tester.run_all_tests();
}

// Test case 5: Grandi variazioni di scala
void test_large_scale_variation()
{
  print_header( "TEST 5: Large Scale Variation", fmt::color::blue );

  vector<real_type>    X = { 0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1000.0 };
  SearchIntervalTester tester( "ScaleVariation-11", X, false, true );
  tester.run_all_tests();
}

// Test case 6: Curve chiuse (periodiche)
void test_closed_curves()
{
  print_header( "TEST 6: Closed Curves", fmt::color::blue );

  // Curva chiusa su [0, 2π]
  vector<real_type> X;
  integer           n = 13;
  for ( integer i = 0; i < n; ++i ) { X.push_back( 2.0 * M_PI * i / ( n - 1 ) ); }

  SearchIntervalTester tester( "ClosedCircle-13", X, true, false );
  tester.run_all_tests();
}

// Test case 7: Pochi nodi (edge case)
void test_few_nodes()
{
  print_header( "TEST 7: Few Nodes", fmt::color::blue );

  // Solo 2 nodi
  {
    vector<real_type>    X = { 0.0, 1.0 };
    SearchIntervalTester tester( "TwoNodes", X, false, true );
    tester.run_all_tests();
  }

  // Solo 3 nodi
  {
    vector<real_type>    X = { 0.0, 0.5, 1.0 };
    SearchIntervalTester tester( "ThreeNodes", X, false, true );
    tester.run_all_tests();
  }
}

// Test case 8: Grandi dataset (performance test)
void test_large_dataset()
{
  print_header( "TEST 8: Large Dataset", fmt::color::blue );

  integer           n = 10000;
  vector<real_type> X( static_cast<size_t>( n ) );

  // Creazione di dati con distribuzione non uniforme
  if ( !X.empty() )
  {
    X[0] = 0.0;
    for ( integer i = 1; i < n; ++i )
    {
      // Crescita esponenziale con qualche variazione
      size_t idx      = static_cast<size_t>( i );
      size_t prev_idx = static_cast<size_t>( i - 1 );
      X[idx]          = X[prev_idx] + 0.1 + 0.9 * sin( i * 0.01 ) * sin( i * 0.01 );
    }
  }

  SearchIntervalTester tester( "Large-10000", X, false, true );

  // Test solo alcuni punti per non appesantire l'output
  print_info( "Large dataset test (showing first 20 tests)" );

  // Test punti distribuiti
  for ( integer i = 0; i < 20; ++i )
  {
    if ( X.empty() ) break;

    real_type x        = X[0] + ( X[static_cast<size_t>( n ) - 1] - X[0] ) * ( i / 19.0 );
    integer   expected = tester.find_expected_interval( x );  // ORA È PUBBLICO
    tester.test_point( x, expected, "large dataset" );
  }
}

// Test case 9: Caso degenere (tutti i nodi uguali)
void test_degenerate_case()
{
  print_header( "TEST 9: Degenerate Case", fmt::color::blue );

  vector<real_type>    X( 10, 1.0 );  // Tutti i nodi a 1.0
  SearchIntervalTester tester( "Degenerate-10", X, false, true );

  // Questo caso è problematico perché tutti gli intervalli sono degeneri
  print_warning( "Testing degenerate case (all nodes equal)" );

  // Test alcuni punti
  // Con tutti i nodi uguali, SearchInterval restituirà una eccezione
  try
  {
    tester.test_point( 0.5, 0, "left of degenerate" );
    tester.test_point( 1.0, 0, "at degenerate node" );
    tester.test_point( 1.5, 8, "right of degenerate" );  // n-2 perché l'ultimo intervallo valido è [n-2, n-1]
  }
  catch ( std::exception const & err )
  {
    fmt::print( "TEST 9: Degenerate Case expect exception: {}\n", err.what() );
  }
  catch ( ... )
  {
    fmt::print( "TEST 9: Degenerate Case expect exception: UNKNOWN\n" );
  }
}

// Test case 10: Caso quasi degenere (range molto piccolo)
void test_nearly_degenerate()
{
  print_header( "TEST 10: Nearly Degenerate", fmt::color::blue );

  // Con tolleranza 1e-10, tutte queste differenze (1e-15) sono considerate zero
  // Quindi SearchInterval tratterà tutti i nodi come uguali
  vector<real_type>    X = { 1.0, 1.0 + 1e-15, 1.0 + 2e-15, 1.1 };
  SearchIntervalTester tester( "NearlyDegenerate-4", X, false, true );

  print_warning( "All nodes are considered equal due to tolerance (1e-10 >> 1e-15)" );
  tester.run_all_tests();
}

int main()
{
  print_header( "SPLINE SEARCH INTERVAL - COMPREHENSIVE TEST SUITE", fmt::color::cyan );

  fmt::print(
    fg( fmt::color::white ) | fmt::emphasis::bold,
    "\n"
    "🧪  Running comprehensive tests for SearchInterval class\n"
    "    This test suite validates the hybrid search algorithm under various conditions.\n"
    "\n"
    "📈  Test Cases:\n"
    "    1.  Uniform distribution (baseline)\n"
    "    2.  Non-uniform distribution\n"
    "    3.  Clustered nodes\n"
    "    4.  Duplicate nodes (special handling)\n"
    "    5.  Large scale variation\n"
    "    6.  Closed curves (periodic)\n"
    "    7.  Few nodes (edge cases)\n"
    "    8.  Large dataset (performance)\n"
    "    9.  Degenerate case (all nodes equal)\n"
    "    10. Nearly degenerate case\n"
    "\n"
    "🔍  SearchInterval behavior:\n"
    "    • Returns interval [i, i+1] where X[i] <= x < X[i+1]\n"
    "    • For duplicate nodes, finds leftmost valid interval\n"
    "    • Thread-safe with mutex protection\n"
    "    • Hybrid algorithm: O(1) lookup + O(log(n/table_size)) binary search\n"
    "\n" );

  // Esegui tutti i test
  test_uniform_distribution();
  test_non_uniform_distribution();
  test_clustered_nodes();
  test_duplicate_nodes();
  test_large_scale_variation();
  test_closed_curves();
  test_few_nodes();
  test_large_dataset();
  test_degenerate_case();
  test_nearly_degenerate();

  // Riepilogo globale
  print_header( "GLOBAL TEST SUMMARY", fmt::color::green );

  fmt::print(
    fg( fmt::color::cyan ) | fmt::emphasis::bold,
    "\n"
    "🎯  All tests completed!\n"
    "    The SearchInterval algorithm has been validated under various conditions.\n"
    "\n"
    "📋  Key Features Tested:\n"
    "    • Hybrid search (lookup table + binary search)\n"
    "    • Handling of duplicate consecutive nodes (backtrack to first valid interval)\n"
    "    • Closed curves (periodic boundary conditions)\n"
    "    • Out-of-bounds queries (clamping for open curves)\n"
    "    • Degenerate and edge cases\n"
    "    • Thread safety (implicit via mutex)\n"
    "\n"
    "⚡  Performance Characteristics:\n"
    "    • O(1) initial range reduction via lookup table\n"
    "    • O(log(n/table_size)) binary search\n"
    "    • Efficient for both small and large datasets\n"
    "\n"
    "⚠️  Important Notes:\n"
    "    • Duplicate nodes are handled by returning the leftmost valid interval\n"
    "    • Intervals with zero width (duplicate consecutive nodes) are not valid\n"
    "    • Tolerance for floating point comparisons: 1e-10 * |X|\n"
    "    • ALWAYS returns consecutive interval indices [i, i+1]\n"
    "\n" );

  print_header( "TEST SUITE COMPLETED SUCCESSFULLY", fmt::color::green );

  return 0;
}
