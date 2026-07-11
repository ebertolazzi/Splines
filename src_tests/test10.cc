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

#include "Splines.hh"
#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

using namespace Splines;
using namespace std;

// Utility per stampare risultati
#define TEST_SECTION( name ) cout << "\n=== " << name << " ===\n"
#define TEST_RESULT( desc, value ) cout << "  " << desc << ": " << value << "\n"
#define TEST_ARRAY( desc, arr, n )  \
  cout << "  " << desc << ": [";    \
  for ( integer i = 0; i < n; ++i ) \
  {                                 \
    if ( i > 0 ) cout << ", ";      \
    cout << arr[i];                 \
  }                                 \
  cout << "]\n"


void test_splinevec()
{
  TEST_SECTION( "SplineVec - Test Completo" );

  // Crea una curva 2D (parabola)
  SplineVec curve( "TestCurve" );

  integer dim  = 2;
  integer npts = 5;

  // Punti: (t, t^2) per t = 0..4
  vector<real_type> X   = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  vector<real_type> Y_x = { 0.0, 1.0, 2.0, 3.0, 4.0 };   // x = t
  vector<real_type> Y_y = { 0.0, 1.0, 4.0, 9.0, 16.0 };  // y = t^2

  real_type const * Y[2] = { Y_x.data(), Y_y.data() };

  curve.setup( dim, npts, Y );
  curve.set_knots( X.data() );
  curve.catmull_rom();

  TEST_RESULT( "Dimension", curve.dimension() );
  TEST_RESULT( "Num points", curve.num_points() );
  TEST_RESULT( "x_min", curve.x_min() );
  TEST_RESULT( "x_max", curve.x_max() );

  // Valutazione in un punto
  real_type t  = 2.5;
  real_type x  = curve.eval( t, 0 );
  real_type y  = curve.eval( t, 1 );
  real_type dx = curve.eval_D( t, 0 );
  real_type dy = curve.eval_D( t, 1 );

  cout << "  t = " << t << ":\n";
  TEST_RESULT( "    x", x );
  TEST_RESULT( "    y", y );
  TEST_RESULT( "    dx/dt", dx );
  TEST_RESULT( "    dy/dt", dy );

  // Valutazione di tutti i componenti
  vector<real_type> point;
  curve.eval( t, point );
  cout << "  Vector evaluation: [" << point[0] << ", " << point[1] << "]\n";

  // Test curvatura (per curve 2D)
  if ( curve.dimension() >= 2 )
  {
    real_type k  = curve.curvature( t );
    real_type dk = curve.curvature_D( t );
    TEST_RESULT( "Curvature", k );
    TEST_RESULT( "Curvature derivative", dk );
  }

  // Test parametrizzazioni alternative
  SplineVec curve2( "TestCurve2" );
  curve2.setup( dim, npts, Y );
  curve2.set_knots_chord_length();
  TEST_RESULT( "Chord length parametrization - x_min", curve2.x_min() );
  TEST_RESULT( "Chord length parametrization - x_max", curve2.x_max() );

  SplineVec curve3( "TestCurve3" );
  curve3.setup( dim, npts, Y );
  curve3.set_knots_centripetal();
  TEST_RESULT( "Centripetal parametrization - x_min", curve3.x_min() );
  TEST_RESULT( "Centripetal parametrization - x_max", curve3.x_max() );

  // Test curve chiuse
  curve.make_closed();
  TEST_RESULT( "After make_closed", curve.is_closed() );

  curve.make_open();
  TEST_RESULT( "After make_open", curve.is_closed() );
}

// ============================================================================
// TEST SPLINE1D
// ============================================================================
void test_spline1d()
{
  TEST_SECTION( "Spline1D - Test Completo" );

  // Dati di test
  vector<real_type> x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
  vector<real_type> y = { 0.0, 1.0, 4.0, 9.0, 16.0, 25.0 };

  Spline1D spline( "TestSpline" );

  // Test 1: Costruzione con vari tipi di spline
  TEST_SECTION( "1. Build e Setup" );

  for ( auto type : { SplineType1D::LINEAR, SplineType1D::CUBIC, SplineType1D::AKIMA, SplineType1D::PCHIP } )
  {
    spline.build( type, x, y );
    cout << "  Built " << spline.type_name() << " spline\n";
    TEST_RESULT( "  Type", static_cast<int>( spline.type() ) );
    TEST_RESULT( "  Order", spline.order() );
    TEST_RESULT( "  Num points", spline.num_points() );
  }

  // Ricostruiamo con CUBIC per i test successivi
  spline.build( SplineType1D::CUBIC, x, y );

  // Test 2: Proprietà Open/Close/Bounded
  TEST_SECTION( "2. Open/Close/Bounded Properties" );
  TEST_RESULT( "Is closed", spline.is_closed() );
  TEST_RESULT( "Is bounded", spline.is_bounded() );
  TEST_RESULT( "Is extended constant", spline.is_extended_constant() );

  spline.make_closed();
  TEST_RESULT( "After make_closed", spline.is_closed() );
  spline.make_opened();
  TEST_RESULT( "After make_opened", spline.is_closed() );

  spline.make_bounded();
  TEST_RESULT( "After make_bounded", spline.is_bounded() );
  spline.make_unbounded();
  TEST_RESULT( "After make_unbounded", spline.is_bounded() );

  spline.make_extended_constant();
  TEST_RESULT( "After make_extended_constant", spline.is_extended_constant() );
  spline.make_extended_not_constant();
  TEST_RESULT( "After make_extended_not_constant", spline.is_extended_constant() );

  // Test 3: Accesso ai nodi
  TEST_SECTION( "3. Node Access" );
  TEST_RESULT( "x_begin", spline.x_begin() );
  TEST_RESULT( "x_end", spline.x_end() );
  TEST_RESULT( "y_begin", spline.y_begin() );
  TEST_RESULT( "y_end", spline.y_end() );
  TEST_RESULT( "x_min", spline.x_min() );
  TEST_RESULT( "x_max", spline.x_max() );
  TEST_RESULT( "y_min", spline.y_min() );
  TEST_RESULT( "y_max", spline.y_max() );

  cout << "  Nodes:\n";
  for ( integer i = 0; i < spline.num_points(); ++i )
  {
    cout << "    Node " << i << ": (" << spline.x_node( i ) << ", " << spline.y_node( i ) << ")\n";
  }

  // Test 4: Valutazione
  TEST_SECTION( "4. Evaluation" );
  vector<real_type> test_x = { 0.0, 0.5, 1.0, 2.5, 4.0, 5.0 };

  for ( auto xi : test_x )
  {
    real_type val   = spline.eval( xi );
    real_type d1    = spline.D( xi );
    real_type d2    = spline.DD( xi );
    real_type d3    = spline.DDD( xi );
    real_type d4    = spline.DDDD( xi );
    real_type d5    = spline.DDDDD( xi );
    real_type val2  = spline( xi );  // operator()
    real_type vald  = spline.eval_D( xi );
    real_type valdd = spline.eval_DD( xi );

    cout << "  x = " << xi << ":\n";
    TEST_RESULT( "    eval", val );
    TEST_RESULT( "    D", d1 );
    TEST_RESULT( "    DD", d2 );
    TEST_RESULT( "    DDD", d3 );
    TEST_RESULT( "    DDDD", d4 );
    TEST_RESULT( "    DDDDD", d5 );
    TEST_RESULT( "    operator()", val2 );
    TEST_RESULT( "    eval_D", vald );
    TEST_RESULT( "    eval_DD", valdd );
  }

  // Test 5: Valutazione con array
  TEST_SECTION( "5. Array Evaluation" );
  real_type dd[2], ddd[3];
  spline.D( 2.5, dd );
  TEST_ARRAY( "D array", dd, 2 );
  spline.DD( 2.5, ddd );
  TEST_ARRAY( "DD array", ddd, 3 );

  // Test 6: Valutazione con indice noto
  TEST_SECTION( "6. Indexed Evaluation" );
  for ( integer i = 0; i < spline.num_points() - 1; ++i )
  {
    real_type xi  = x[static_cast<size_t>( i )] +
                    0.5 * ( x[static_cast<size_t>( i + 1 )] - x[static_cast<size_t>( i )] );
    real_type val = spline.id_eval( i, xi );
    real_type d1  = spline.id_D( i, xi );
    real_type d2  = spline.id_DD( i, xi );
    real_type d3  = spline.id_DDD( i, xi );
    // real_type d4  = spline.id_DDDD( i, xi );  // Not used
    // real_type d5  = spline.id_DDDDD( i, xi ); // Not used

    cout << "  Segment " << i << " at x = " << xi << ":\n";
    TEST_RESULT( "    id_eval", val );
    TEST_RESULT( "    id_D", d1 );
    TEST_RESULT( "    id_DD", d2 );
    TEST_RESULT( "    id_DDD", d3 );
  }

  // Test 7: Manipolazione range
  TEST_SECTION( "7. Range Manipulation" );
  TEST_RESULT( "Original x_min", spline.x_min() );
  TEST_RESULT( "Original x_max", spline.x_max() );

  spline.set_origin( 10.0 );
  TEST_RESULT( "After set_origin(10.0) - x_min", spline.x_min() );
  TEST_RESULT( "After set_origin(10.0) - x_max", spline.x_max() );

  spline.set_range( 0.0, 10.0 );
  TEST_RESULT( "After set_range(0.0, 10.0) - x_min", spline.x_min() );
  TEST_RESULT( "After set_range(0.0, 10.0) - x_max", spline.x_max() );

  // Test 8: Push/Pop
  TEST_SECTION( "8. Incremental Build" );
  Spline1D spline2( "IncrementalSpline" );

  // CORREZIONE: Prima costruiamo una spline vuota o con dati minimali
  // Opzione 1: Costruire prima con dati minimali, poi fare push_back
  vector<real_type> x_init = { 0.0, 1.0 };
  vector<real_type> y_init = { 0.0, 1.0 };
  spline2.build( SplineType1D::LINEAR, x_init, y_init );

  // Ora possiamo fare push_back
  for ( integer i = 2; i < 5; ++i )
  {
    spline2.push_back( static_cast<real_type>( i ), static_cast<real_type>( i * i ) );
  }
  TEST_RESULT( "After push_back - num_points", spline2.num_points() );

  // Test drop_back
  spline2.drop_back();
  TEST_RESULT( "After drop_back - num_points", spline2.num_points() );

  // Ricostruiamo per il test successivo
  spline2.build( SplineType1D::LINEAR, x, y );

  // Test 9: Coefficienti
  TEST_SECTION( "9. Coefficients" );
  integer           n_seg = spline.num_points() - 1;
  integer           order = spline.order();
  vector<real_type> cfs( static_cast<size_t>( n_seg * order ) );
  vector<real_type> nodes( static_cast<size_t>( n_seg + 1 ) );

  spline.coeffs( cfs.data(), nodes.data(), false );
  cout << "  Coefficients (first segment):\n";
  for ( integer i = 0; i < order; ++i ) { cout << "    c[" << i << "] = " << cfs[static_cast<size_t>( i )] << "\n"; }

  // Test 10: Info e dump
  TEST_SECTION( "10. Info and Dump" );
  cout << "  " << spline.info() << "\n";
  spline.dump( cout, 3, "X\tY" );

  // Test 11: Clear
  TEST_SECTION( "11. Clear" );
  spline2.clear();
  TEST_RESULT( "After clear - num_points", spline2.num_points() );

#ifdef AUTODIFF_SUPPORT
  // Test 12: Autodiff
  TEST_SECTION( "12. Autodiff Support" );
  autodiff::dual1st xd1;
  autodiff::dual2nd xd2;
  xd1.val       = 2.5;
  xd1.grad      = 1;
  xd2.val.val   = 2.5;
  xd2.val.grad  = 1;
  xd2.grad.val  = 1;
  xd2.grad.grad = 0;

  auto vd1 = spline.eval( xd1 );
  auto vd2 = spline.eval( xd2 );

  TEST_RESULT( "  dual1st eval", vd1 );
  TEST_RESULT( "  dual2nd eval", vd2 );
#endif
}

// ============================================================================
// TEST SPLINE1DBLEND
// ============================================================================
void test_spline1dblend()
{
  TEST_SECTION( "Spline1Dblend - Test Completo" );

  // Dati di test per due spline
  vector<real_type> x0 = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  vector<real_type> y0 = { 0.0, 1.0, 4.0, 9.0, 16.0 };

  vector<real_type> x1 = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  vector<real_type> y1 = { 0.0, 2.0, 8.0, 18.0, 32.0 };

  Spline1Dblend blend( "TestBlend" );

  // Test 1: Build
  TEST_SECTION( "1. Build" );
  blend.build( SplineType1D::CUBIC, x0, y0, SplineType1D::CUBIC, x1, y1 );

  TEST_RESULT( "Num points spline0", blend.num_points0() );
  TEST_RESULT( "Num points spline1", blend.num_points1() );
  TEST_RESULT( "Order spline0", blend.order0() );
  TEST_RESULT( "Order spline1", blend.order1() );
  cout << "  Type spline0: " << blend.type_name0() << "\n";
  cout << "  Type spline1: " << blend.type_name1() << "\n";

  // Test 2: Accesso alle spline
  TEST_SECTION( "2. Access to Individual Splines" );
  auto const & sp0 = blend.get_spline0();
  auto const & sp1 = blend.get_spline1();

  TEST_RESULT( "Spline0 x_begin", sp0.x_begin() );
  TEST_RESULT( "Spline1 x_begin", sp1.x_begin() );
  TEST_RESULT( "Spline0 y_begin", sp0.y_begin() );
  TEST_RESULT( "Spline1 y_begin", sp1.y_begin() );

  // Test 3: Blend begin/end
  TEST_SECTION( "3. Blend Begin/End Values" );
  vector<real_type> s_values = { 0.0, 0.25, 0.5, 0.75, 1.0 };

  for ( auto s : s_values )
  {
    cout << "  s = " << s << ":\n";
    TEST_RESULT( "    x_begin", blend.x_begin( s ) );
    TEST_RESULT( "    x_end", blend.x_end( s ) );
    TEST_RESULT( "    y_begin", blend.y_begin( s ) );
    TEST_RESULT( "    y_end", blend.y_end( s ) );
  }

  // Test 4: Valutazione blend
  TEST_SECTION( "4. Blend Evaluation" );
  vector<real_type> x_test = { 0.0, 1.0, 2.0, 3.0, 4.0 };

  for ( auto xi : x_test )
  {
    cout << "  x = " << xi << ":\n";
    for ( auto s : s_values )
    {
      real_type val = blend.eval( xi, s );
      real_type d1  = blend.D( xi, s );
      real_type d2  = blend.DD( xi, s );
      // real_type d3   = blend.DDD( xi, s );   // Not used
      // real_type d4   = blend.DDDD( xi, s );  // Not used
      // real_type d5   = blend.DDDDD( xi, s ); // Not used
      real_type val2 = blend( xi, s );

      cout << "    s = " << s << ": eval = " << val << ", D = " << d1 << ", DD = " << d2;
      cout << ", operator() = " << val2 << "\n";
    }
  }

  // Test 5: Valutazione aliases
  TEST_SECTION( "5. Evaluation Aliases" );
  real_type xi = 2.0, s = 0.5;
  TEST_RESULT( "eval", blend.eval( xi, s ) );
  TEST_RESULT( "eval_D", blend.eval_D( xi, s ) );
  TEST_RESULT( "eval_DD", blend.eval_DD( xi, s ) );
  TEST_RESULT( "eval_DDD", blend.eval_DDD( xi, s ) );
  TEST_RESULT( "eval_DDDD", blend.eval_DDDD( xi, s ) );
  TEST_RESULT( "eval_DDDDD", blend.eval_DDDDD( xi, s ) );

  // Test 6: Array evaluation
  TEST_SECTION( "6. Array Evaluation" );
  real_type dd[2], ddd[3];
  blend.D( xi, s, dd );
  TEST_ARRAY( "D array", dd, 2 );
  blend.DD( xi, s, ddd );
  TEST_ARRAY( "DD array", ddd, 3 );

  // Test 7: Manipolazione range
  TEST_SECTION( "7. Range Manipulation" );
  blend.set_origin( 5.0 );
  TEST_RESULT( "After set_origin - spline0 x_begin", blend.get_spline0().x_begin() );
  TEST_RESULT( "After set_origin - spline1 x_begin", blend.get_spline1().x_begin() );

  blend.set_range( 0.0, 10.0 );
  TEST_RESULT( "After set_range - spline0 x_min", blend.get_spline0().x_min() );
  TEST_RESULT( "After set_range - spline0 x_max", blend.get_spline0().x_max() );

  // Test 8: Verifica linearità del blend
  TEST_SECTION( "8. Linearity Check" );
  xi            = 2.0;
  real_type v0  = sp0.eval( xi );
  real_type v1  = sp1.eval( xi );
  real_type vb0 = blend.eval( xi, 0.0 );
  real_type vb1 = blend.eval( xi, 1.0 );
  real_type vbh = blend.eval( xi, 0.5 );

  TEST_RESULT( "Spline0 at x=2.0", v0 );
  TEST_RESULT( "Spline1 at x=2.0", v1 );
  TEST_RESULT( "Blend(x=2.0, s=0.0)", vb0 );
  TEST_RESULT( "Blend(x=2.0, s=1.0)", vb1 );
  TEST_RESULT( "Blend(x=2.0, s=0.5)", vbh );
  TEST_RESULT( "Expected (v0+v1)/2", ( v0 + v1 ) / 2.0 );
  TEST_RESULT( "Error", abs( vbh - ( v0 + v1 ) / 2.0 ) );

  // Test 10: Build con GenericContainer
  TEST_SECTION( "10. Build with GenericContainer" );
  try
  {
    GenericContainer gc;
    gc["spline0"]["spline_type"] = "cubic";
    gc["spline0"]["xdata"]       = x0;
    gc["spline0"]["ydata"]       = y0;
    gc["spline1"]["spline_type"] = "cubic";
    gc["spline1"]["xdata"]       = x1;
    gc["spline1"]["ydata"]       = y1;

    Spline1Dblend blend2( "TestBlend2" );
    blend2.build( gc );
    TEST_RESULT( "Build from GenericContainer - num_points0", blend2.num_points0() );
    TEST_RESULT( "Build from GenericContainer - num_points1", blend2.num_points1() );
  }
  catch ( exception const & e )
  {
    cout << "  Exception: " << e.what() << "\n";
  }
}

// ============================================================================
// MAIN
// ============================================================================
int main()
{
  cout << fixed << setprecision( 6 );

  try
  {
    test_splinevec();
    test_spline1d();
    test_spline1dblend();

    cout << "\n\n=== ALL TESTS COMPLETED ===\n";
  }
  catch ( exception const & e )
  {
    cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
