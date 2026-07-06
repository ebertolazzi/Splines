#include "Splines/Splines.hh"
#include "Splines/SplinesCinterface.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace Splines;

namespace
{
  int constexpr    EXIT_OK  = 0;
  int constexpr    EXIT_KO  = 1;
  real_type constexpr EPS   = 1e-12;
  real_type constexpr EPS_D = 1e-10;

  int failures = 0;

  void
  check_true( bool const cond, std::string const & msg )
  {
    if ( cond ) return;
    ++failures;
    std::cerr << "FAIL: " << msg << '\n';
  }

  void
  check_near( real_type const got, real_type const expected, real_type const tol, std::string const & msg )
  {
    if ( std::abs( got - expected ) <= tol ) return;
    ++failures;
    std::cerr << "FAIL: " << msg << " got=" << got << " expected=" << expected << " tol=" << tol << '\n';
  }

  void
  check_nonempty_cstr( char const * s, std::string const & msg )
  {
    check_true( s != nullptr && std::strlen( s ) > 0, msg );
  }

  void
  check_map_real(
    map_type const &      map,
    std::string const &   key,
    real_type const       expected,
    real_type const       tol,
    std::string const &   msg )
  {
    auto const it = map.find( key );
    check_true( it != map.end(), msg + " missing key '" + key + "'" );
    if ( it == map.end() ) return;
    check_near( it->second.get_real(), expected, tol, msg + " key '" + key + "'" );
  }

  void
  check_map_vec(
    map_type const &                map,
    std::string const &             key,
    std::vector<real_type> const &  expected,
    real_type const                 tol,
    std::string const &             msg )
  {
    auto const it = map.find( key );
    check_true( it != map.end(), msg + " missing key '" + key + "'" );
    if ( it == map.end() ) return;
    auto const & got = it->second.get_vec_real();
    check_true( got.size() == expected.size(), msg + " wrong size for key '" + key + "'" );
    if ( got.size() != expected.size() ) return;
    for ( size_t i = 0; i < got.size(); ++i )
      check_near( got[i], expected[i], tol, msg + " key '" + key + "' idx " + std::to_string( i ) );
  }

  std::vector<real_type>
  make_plane_data(
    std::vector<real_type> const & x,
    std::vector<real_type> const & y,
    real_type (*f)( real_type, real_type ) )
  {
    std::vector<real_type> z( x.size() * y.size() );
    for ( size_t i = 0; i < x.size(); ++i )
      for ( size_t j = 0; j < y.size(); ++j ) z[i * y.size() + j] = f( x[i], y[j] );
    return z;
  }

  real_type
  plane0( real_type const x, real_type const y )
  {
    return x + 2 * y;
  }

  real_type
  plane1( real_type const x, real_type const y )
  {
    return 10 + 3 * x - y;
  }

  void
  test_c_interface()
  {
    std::cout << "test16: C interface regressions\n";

    check_true( SPLINE_new( "bad_type", "does_not_exist" ) == -1, "SPLINE_new must reject invalid spline types" );
    check_true( SPLINE_mem_ptr( "missing_id" ) == nullptr, "SPLINE_mem_ptr must return nullptr for unknown ids" );

    real_type const x[] = { 0, 1, 2, 3 };
    real_type const y[] = { 0, 1, 2, 3 };

    check_true( SPLINE_new( "first", "linear" ) == 0, "SPLINE_new linear spline" );
    check_true( SPLINE_build2( x, y, 4 ) == 0, "SPLINE_build2 linear spline" );
    check_true( std::string( SPLINE_get_type_name() ) == "SPLINE_LINEAR", "SPLINE_get_type_name after build" );
    check_near( SPLINE_eval( 1.5 ), 1.5, EPS, "SPLINE_eval on linear spline" );
    check_near( SPLINE_eval_D( 1.5 ), 1.0, EPS, "SPLINE_eval_D on linear spline" );
    check_near( SPLINE_eval_DD( 1.5 ), 0.0, EPS, "SPLINE_eval_DD on linear spline" );
    check_near( SPLINE_eval_DDD( 1.5 ), 0.0, EPS, "SPLINE_eval_DDD on linear spline" );
    check_near( SPLINE_eval_DDDD( 1.5 ), 0.0, EPS, "SPLINE_eval_DDDD on linear spline" );
    check_near( SPLINE_eval_DDDDD( 1.5 ), 0.0, EPS, "SPLINE_eval_DDDDD on linear spline" );
    check_true( SPLINE_print() == 0, "SPLINE_print on selected spline" );
    check_true( SPLINE_mem_ptr( "first" ) != nullptr, "SPLINE_mem_ptr must expose the internal spline pointer" );

    real_type const y_cubic[] = { 0, 1, 4, 9 };
    check_true( SPLINE_new( "second", "cubic" ) == 0, "SPLINE_new cubic spline" );
    check_true( SPLINE_build2( x, y_cubic, 4 ) == 0, "SPLINE_build2 cubic spline" );
    check_true( SPLINE_select( "first" ) == 0, "SPLINE_select first spline" );
    check_true( SPLINE_delete( "second" ) == 0, "SPLINE_delete on unselected spline" );
    check_near( SPLINE_eval( 1.5 ), 1.5, EPS, "Deleting another spline must not reset the selected head" );

    check_true( SPLINE_new( "incremental", "linear" ) == 0, "SPLINE_new incremental spline" );
    check_true( SPLINE_init() == 0, "SPLINE_init on incremental spline" );
    check_true( SPLINE_push( 0, 0 ) == 0, "SPLINE_push p0" );
    check_true( SPLINE_push( 1, 2 ) == 0, "SPLINE_push p1" );
    check_true( SPLINE_push( 2, 4 ) == 0, "SPLINE_push p2" );
    check_true( SPLINE_build() == 0, "SPLINE_build incremental spline" );
    check_near( SPLINE_eval( 1.5 ), 3.0, EPS, "SPLINE_build after push_back" );

    check_true( SPLINE_delete( "incremental" ) == 0, "SPLINE_delete incremental spline" );
    check_true( SPLINE_delete( "first" ) == 0, "SPLINE_delete first spline" );
  }

  void
  test_spline1dblend_api()
  {
    std::cout << "test16: Spline1Dblend API coverage\n";

    real_type const x0[] = { 0, 1, 2, 3 };
    real_type const y0[] = { 0, 1, 2, 3 };
    real_type const x1[] = { 0, 1, 2, 3 };
    real_type const y1[] = { 10, 12, 14, 16 };

    real_type const x0_stride[] = { 0, -99, 1, -99, 2, -99, 3, -99 };
    real_type const y0_stride[] = { 0, -99, 1, -99, 2, -99, 3, -99 };
    real_type const x1_stride[] = { 0, -99, 1, -99, 2, -99, 3, -99 };
    real_type const y1_stride[] = { 10, -99, 12, -99, 14, -99, 16, -99 };

    real_type const x = 1.5;
    real_type const s = 0.25;

    real_type const expected_value = ( 1 - s ) * x + s * ( 10 + 2 * x );
    real_type const expected_D     = ( 1 - s ) * 1.0 + s * 2.0;

    Spline1Dblend blend_raw( "blend_raw" );
    blend_raw.build( SplineType1D::LINEAR, x0, y0, 4, SplineType1D::LINEAR, x1, y1, 4 );

    check_true( blend_raw.get_spline0().name() == "blend_raw_0", "Spline1Dblend get_spline0" );
    check_true( blend_raw.get_spline1().name() == "blend_raw_1", "Spline1Dblend get_spline1" );
    check_true( blend_raw.num_points0() == 4, "Spline1Dblend num_points0" );
    check_true( blend_raw.num_points1() == 4, "Spline1Dblend num_points1" );
    check_true( blend_raw.order0() == 2, "Spline1Dblend order0" );
    check_true( blend_raw.order1() == 2, "Spline1Dblend order1" );
    check_nonempty_cstr( blend_raw.type_name0(), "Spline1Dblend type_name0" );
    check_nonempty_cstr( blend_raw.type_name1(), "Spline1Dblend type_name1" );
    check_true( blend_raw.type0() == SplineType1D::LINEAR, "Spline1Dblend type0" );
    check_true( blend_raw.type1() == SplineType1D::LINEAR, "Spline1Dblend type1" );
    check_near( blend_raw.x_begin( s ), 0.0, EPS, "Spline1Dblend x_begin" );
    check_near( blend_raw.x_end( s ), 3.0, EPS, "Spline1Dblend x_end" );
    check_near( blend_raw.y_begin( s ), 2.5, EPS, "Spline1Dblend y_begin" );
    check_near( blend_raw.y_end( s ), 6.25, EPS, "Spline1Dblend y_end" );

    check_near( blend_raw.eval( x, s ), expected_value, EPS, "Spline1Dblend eval" );
    check_near( blend_raw( x, s ), expected_value, EPS, "Spline1Dblend operator()" );
    check_near( blend_raw.D( x, s ), expected_D, EPS, "Spline1Dblend D" );
    check_near( blend_raw.DD( x, s ), 0.0, EPS, "Spline1Dblend DD" );
    check_near( blend_raw.DDD( x, s ), 0.0, EPS, "Spline1Dblend DDD" );
    check_near( blend_raw.DDDD( x, s ), 0.0, EPS, "Spline1Dblend DDDD" );
    check_near( blend_raw.DDDDD( x, s ), 0.0, EPS, "Spline1Dblend DDDDD" );
    check_near( blend_raw.eval_D( x, s ), expected_D, EPS, "Spline1Dblend eval_D" );
    check_near( blend_raw.eval_DD( x, s ), 0.0, EPS, "Spline1Dblend eval_DD" );
    check_near( blend_raw.eval_DDD( x, s ), 0.0, EPS, "Spline1Dblend eval_DDD" );
    check_near( blend_raw.eval_DDDD( x, s ), 0.0, EPS, "Spline1Dblend eval_DDDD" );
    check_near( blend_raw.eval_DDDDD( x, s ), 0.0, EPS, "Spline1Dblend eval_DDDDD" );

    real_type d[2], dd[3];
    blend_raw.D( x, s, d );
    blend_raw.DD( x, s, dd );
    check_near( d[0], expected_value, EPS, "Spline1Dblend D array value" );
    check_near( d[1], expected_D, EPS, "Spline1Dblend D array derivative" );
    check_near( dd[0], expected_value, EPS, "Spline1Dblend DD array value" );
    check_near( dd[1], expected_D, EPS, "Spline1Dblend DD array derivative" );
    check_near( dd[2], 0.0, EPS, "Spline1Dblend DD array second derivative" );

#ifdef AUTODIFF_SUPPORT
    float xf = static_cast<float>( x );
    check_near( blend_raw.eval( xf, s ), expected_value, EPS, "Spline1Dblend generic eval<float>" );
#endif

    blend_raw.set_origin( 5.0 );
    check_near( blend_raw.get_spline0().x_begin(), 5.0, EPS, "Spline1Dblend set_origin spline0" );
    check_near( blend_raw.get_spline1().x_begin(), 5.0, EPS, "Spline1Dblend set_origin spline1" );
    blend_raw.set_range( -2.0, 2.0 );
    check_near( blend_raw.get_spline0().x_begin(), -2.0, EPS, "Spline1Dblend set_range spline0 begin" );
    check_near( blend_raw.get_spline0().x_end(), 2.0, EPS, "Spline1Dblend set_range spline0 end" );
    check_near( blend_raw.get_spline1().x_begin(), -2.0, EPS, "Spline1Dblend set_range spline1 begin" );
    check_near( blend_raw.get_spline1().x_end(), 2.0, EPS, "Spline1Dblend set_range spline1 end" );

    Spline1Dblend blend_stride( "blend_stride" );
    blend_stride.build(
      SplineType1D::LINEAR,
      x0_stride,
      2,
      y0_stride,
      2,
      4,
      SplineType1D::LINEAR,
      x1_stride,
      2,
      y1_stride,
      2,
      4 );
    check_near( blend_stride.eval( x, s ), expected_value, EPS, "Spline1Dblend build stride overload" );

    Spline1Dblend blend_vec( "blend_vec" );
    blend_vec.build(
      SplineType1D::LINEAR,
      std::vector<real_type>{ 0, 1, 2, 3 },
      std::vector<real_type>{ 0, 1, 2, 3 },
      SplineType1D::LINEAR,
      std::vector<real_type>{ 0, 1, 2, 3 },
      std::vector<real_type>{ 10, 12, 14, 16 } );
    check_near( blend_vec.eval( x, s ), expected_value, EPS, "Spline1Dblend build vector overload" );

    GenericContainer gc;
    gc["spline0"]["spline_type"] = std::string( "linear" );
    gc["spline0"]["xdata"]       = std::vector<real_type>{ 0, 1, 2, 3 };
    gc["spline0"]["ydata"]       = std::vector<real_type>{ 0, 1, 2, 3 };
    gc["spline1"]["spline_type"] = std::string( "linear" );
    gc["spline1"]["xdata"]       = std::vector<real_type>{ 0, 1, 2, 3 };
    gc["spline1"]["ydata"]       = std::vector<real_type>{ 10, 12, 14, 16 };

    Spline1Dblend blend_gc( "blend_gc" );
    blend_gc.build( gc );
    check_near( blend_gc.eval( x, s ), expected_value, EPS, "Spline1Dblend build(gc)" );
  }

  void
  test_spline2d_and_surface_setup()
  {
    std::cout << "test16: Spline2D setup defaults and surface guards\n";

    std::vector<real_type> const x = { 0, 1, 2 };
    std::vector<real_type> const y = { 0, 1, 2 };
    std::vector<real_type> const z = make_plane_data( x, y, plane0 );

    GenericContainer gc;
    gc["spline_type"] = std::string( "bilinear" );
    gc["xdata"]       = x;
    gc["ydata"]       = y;
    gc["zdata"]       = z;

    Spline2D surface( "surface" );
    surface.build( gc );

    check_true( surface.name() == "surface", "Spline2D name" );
    check_true( surface.num_point_x() == 3, "Spline2D num_point_x" );
    check_true( surface.num_point_y() == 3, "Spline2D num_point_y" );
    check_near( surface.x_node( 2 ), 2.0, EPS, "Spline2D x_node" );
    check_near( surface.y_node( 1 ), 1.0, EPS, "Spline2D y_node" );
    check_near( surface.z_node( 2, 1 ), plane0( 2, 1 ), EPS, "Spline2D z_node" );
    check_near( surface.x_min(), 0.0, EPS, "Spline2D x_min" );
    check_near( surface.x_max(), 2.0, EPS, "Spline2D x_max" );
    check_near( surface.y_min(), 0.0, EPS, "Spline2D y_min" );
    check_near( surface.y_max(), 2.0, EPS, "Spline2D y_max" );
    check_near( surface.z_min(), plane0( 0, 0 ), EPS, "Spline2D z_min" );
    check_near( surface.z_max(), plane0( 2, 2 ), EPS, "Spline2D z_max" );
    check_near( surface.eval( 0.5, 1.25 ), plane0( 0.5, 1.25 ), EPS, "Spline2D eval" );
    check_near( surface( 0.5, 1.25 ), plane0( 0.5, 1.25 ), EPS, "Spline2D operator()" );
    check_near( surface.Dx( 0.5, 1.25 ), 1.0, EPS, "Spline2D Dx" );
    check_near( surface.Dy( 0.5, 1.25 ), 2.0, EPS, "Spline2D Dy" );
    check_near( surface.eval_D_1( 0.5, 1.25 ), 1.0, EPS, "Spline2D eval_D_1" );
    check_near( surface.eval_D_2( 0.5, 1.25 ), 2.0, EPS, "Spline2D eval_D_2" );
    check_near( surface.Dxx( 0.5, 1.25 ), 0.0, EPS, "Spline2D Dxx" );
    check_near( surface.Dxy( 0.5, 1.25 ), 0.0, EPS, "Spline2D Dxy" );
    check_near( surface.Dyy( 0.5, 1.25 ), 0.0, EPS, "Spline2D Dyy" );
    check_near( surface.eval_D_1_1( 0.5, 1.25 ), 0.0, EPS, "Spline2D eval_D_1_1" );
    check_near( surface.eval_D_1_2( 0.5, 1.25 ), 0.0, EPS, "Spline2D eval_D_1_2" );
    check_near( surface.eval_D_2_2( 0.5, 1.25 ), 0.0, EPS, "Spline2D eval_D_2_2" );

    real_type d[3], dd[6];
    surface.D( 0.5, 1.25, d );
    surface.DD( 0.5, 1.25, dd );
    check_near( d[0], plane0( 0.5, 1.25 ), EPS, "Spline2D D value" );
    check_near( d[1], 1.0, EPS, "Spline2D D dx" );
    check_near( d[2], 2.0, EPS, "Spline2D D dy" );
    check_near( dd[0], plane0( 0.5, 1.25 ), EPS, "Spline2D DD value" );
    check_near( dd[1], 1.0, EPS, "Spline2D DD dx" );
    check_near( dd[2], 2.0, EPS, "Spline2D DD dy" );
    check_near( dd[3], 0.0, EPS, "Spline2D DD dxx" );
    check_near( dd[4], 0.0, EPS, "Spline2D DD dxy" );
    check_near( dd[5], 0.0, EPS, "Spline2D DD dyy" );

    surface.make_x_closed();
    surface.make_y_closed();
    surface.make_x_bounded();
    surface.make_y_bounded();
    check_true( surface.is_x_closed(), "Spline2D make_x_closed" );
    check_true( surface.is_y_closed(), "Spline2D make_y_closed" );
    check_true( surface.is_x_bounded(), "Spline2D make_x_bounded" );
    check_true( surface.is_y_bounded(), "Spline2D make_y_bounded" );
    surface.make_x_opened();
    surface.make_y_opened();
    surface.make_x_unbounded();
    surface.make_y_unbounded();
    check_true( !surface.is_x_closed(), "Spline2D make_x_opened" );
    check_true( !surface.is_y_closed(), "Spline2D make_y_opened" );
    check_true( !surface.is_x_bounded(), "Spline2D make_x_unbounded" );
    check_true( !surface.is_y_bounded(), "Spline2D make_y_unbounded" );

    std::ostringstream os;
    os << surface.info() << '\n';
    surface.info( os );
    surface.write_to_stream( os );
    surface.dump_data( os );
    check_true( !os.str().empty(), "Spline2D info/write/dump output" );
    check_true( std::string( surface.type_name() ) == "bilinear", "Spline2D type_name" );

#ifdef AUTODIFF_SUPPORT
    float xf = 0.5f;
    float yf = 1.25f;
    check_near(
      surface.eval( xf, yf ),
      plane0( static_cast<real_type>( xf ), static_cast<real_type>( yf ) ),
      EPS,
      "Spline2D generic eval<float,float>" );
    check_near(
      surface( xf, yf ),
      plane0( static_cast<real_type>( xf ), static_cast<real_type>( yf ) ),
      EPS,
      "Spline2D generic operator<float,float>" );
#endif

    Spline2D surface_vec( "surface_vec" );
    surface_vec.build( SplineType2D::BILINEAR, x, y, z, false, false );
    check_near( surface_vec.eval( 0.5, 1.25 ), plane0( 0.5, 1.25 ), EPS, "Spline2D vector build overload" );
    surface_vec.clear();
    surface_vec.build( SplineType2D::BILINEAR, x, y, z, false, false );
    check_near( surface_vec.eval( 0.5, 1.25 ), plane0( 0.5, 1.25 ), EPS, "Spline2D clear and rebuild" );

    Spline2D surface_raw( "surface_raw" );
    surface_raw.build( SplineType2D::BILINEAR, x.data(), 1, y.data(), 1, z.data(), 3, 3, 3, false, false );
    check_near( surface_raw.eval( 0.5, 1.25 ), plane0( 0.5, 1.25 ), EPS, "Spline2D raw build overload" );

    Spline2D surface_uniform_raw( "surface_uniform_raw" );
    surface_uniform_raw.build( SplineType2D::BILINEAR, z.data(), 3, 3, 3, false, false );
    check_near(
      surface_uniform_raw.eval( 0.5, 1.25 ),
      plane0( 0.5, 1.25 ),
      EPS,
      "Spline2D raw z-only build overload" );

    Spline2D surface_uniform_vec( "surface_uniform_vec" );
    surface_uniform_vec.build( SplineType2D::BILINEAR, z, 3, 3, false, false );
    check_near(
      surface_uniform_vec.eval( 0.5, 1.25 ),
      plane0( 0.5, 1.25 ),
      EPS,
      "Spline2D vector z-only build overload" );

    bool bad_z_caught = false;
    try
    {
      Spline2D bad_surface( "bad_surface" );
      std::vector<real_type> bad_z = z;
      bad_z.pop_back();
      bad_surface.build( SplineType2D::BILINEAR, x, y, bad_z, false, false );
    }
    catch ( std::exception const & )
    {
      bad_z_caught = true;
    }
    check_true( bad_z_caught, "SplineSurf vector build must reject z.size() != nx*ny" );
  }

  void
  test_spline2dblend_api()
  {
    std::cout << "test16: Spline2Dblend API coverage\n";

    integer const nx = 3;
    integer const ny = 3;
    std::vector<real_type> z0;
    std::vector<real_type> z1;
    z0.reserve( static_cast<size_t>( nx * ny ) );
    z1.reserve( static_cast<size_t>( nx * ny ) );

    for ( integer i = 0; i < nx; ++i )
      for ( integer j = 0; j < ny; ++j )
      {
        z0.push_back( plane0( i, j ) );
        z1.push_back( plane1( i, j ) );
      }

    Spline2Dblend blend( "blend2d" );
    blend.build(
      SplineType2D::BILINEAR,
      z0,
      nx,
      ny,
      false,
      false,
      SplineType2D::BILINEAR,
      z1,
      nx,
      ny,
      false,
      false );

    real_type const x = 1.25;
    real_type const y = 0.5;
    real_type const s = 0.25;
    real_type const expected_value = ( 1 - s ) * plane0( x, y ) + s * plane1( x, y );
    real_type const expected_dx    = ( 1 - s ) * 1.0 + s * 3.0;
    real_type const expected_dy    = ( 1 - s ) * 2.0 + s * ( -1.0 );

    check_true( blend.num_point_x0() == nx, "Spline2Dblend num_point_x0" );
    check_true( blend.num_point_x1() == nx, "Spline2Dblend num_point_x1" );
    check_true( blend.num_point_y0() == ny, "Spline2Dblend num_point_y0" );
    check_true( blend.num_point_y1() == ny, "Spline2Dblend num_point_y1" );
    check_near( blend.x_min(), 0.0, EPS, "Spline2Dblend x_min" );
    check_near( blend.x_max(), 2.0, EPS, "Spline2Dblend x_max" );
    check_near( blend.y_min(), 0.0, EPS, "Spline2Dblend y_min" );
    check_near( blend.y_max(), 2.0, EPS, "Spline2Dblend y_max" );
    check_true( std::string( blend.get_surf0().type_name() ) == "bilinear", "Spline2Dblend surf0 type" );
    check_true( std::string( blend.get_surf1().type_name() ) == "bilinear", "Spline2Dblend surf1 type" );

    check_near( blend.eval( x, y, s ), expected_value, EPS, "Spline2Dblend eval" );
    check_near( blend( x, y, s ), expected_value, EPS, "Spline2Dblend operator()" );
    check_near( blend.Dx( x, y, s ), expected_dx, EPS, "Spline2Dblend Dx" );
    check_near( blend.Dy( x, y, s ), expected_dy, EPS, "Spline2Dblend Dy" );
    check_near( blend.eval_D_1( x, y, s ), expected_dx, EPS, "Spline2Dblend eval_D_1" );
    check_near( blend.eval_D_2( x, y, s ), expected_dy, EPS, "Spline2Dblend eval_D_2" );
    check_near( blend.Dxx( x, y, s ), 0.0, EPS, "Spline2Dblend Dxx" );
    check_near( blend.Dxy( x, y, s ), 0.0, EPS, "Spline2Dblend Dxy" );
    check_near( blend.Dyy( x, y, s ), 0.0, EPS, "Spline2Dblend Dyy" );
    check_near( blend.eval_D_1_1( x, y, s ), 0.0, EPS, "Spline2Dblend eval_D_1_1" );
    check_near( blend.eval_D_1_2( x, y, s ), 0.0, EPS, "Spline2Dblend eval_D_1_2" );
    check_near( blend.eval_D_2_2( x, y, s ), 0.0, EPS, "Spline2Dblend eval_D_2_2" );

    real_type d[3], dd[6];
    blend.D( x, y, s, d );
    blend.DD( x, y, s, dd );
    check_near( d[0], expected_value, EPS, "Spline2Dblend D value" );
    check_near( d[1], expected_dx, EPS, "Spline2Dblend D dx" );
    check_near( d[2], expected_dy, EPS, "Spline2Dblend D dy" );
    check_near( dd[0], expected_value, EPS, "Spline2Dblend DD value" );
    check_near( dd[1], expected_dx, EPS, "Spline2Dblend DD dx" );
    check_near( dd[2], expected_dy, EPS, "Spline2Dblend DD dy" );
    check_near( dd[3], 0.0, EPS, "Spline2Dblend DD dxx" );
    check_near( dd[4], 0.0, EPS, "Spline2Dblend DD dxy" );
    check_near( dd[5], 0.0, EPS, "Spline2Dblend DD dyy" );

    float xf = static_cast<float>( x );
    float yf = static_cast<float>( y );
    auto  vf = blend.eval( xf, yf, s );
    check_near( vf, expected_value, EPS, "Spline2Dblend generic eval<T1,T2>" );

    GenericContainer gc;
    gc["surf0"]["spline_type"] = std::string( "bilinear" );
    gc["surf0"]["zdata"]       = z0;
    gc["surf1"]["spline_type"] = std::string( "bilinear" );
    gc["surf1"]["zdata"]       = z1;
    gc["surf0"]["xdata"]       = std::vector<real_type>{ 0, 1, 2 };
    gc["surf0"]["ydata"]       = std::vector<real_type>{ 0, 1, 2 };
    gc["surf1"]["xdata"]       = std::vector<real_type>{ 0, 1, 2 };
    gc["surf1"]["ydata"]       = std::vector<real_type>{ 0, 1, 2 };

    Spline2Dblend blend_gc( "blend2d_gc" );
    blend_gc.build( gc );
    check_near( blend_gc.eval( x, y, s ), expected_value, EPS, "Spline2Dblend build(gc)" );

    std::vector<real_type> const grid = { 0, 1, 2 };

    Spline2Dblend blend_vec_full( "blend2d_vec_full" );
    blend_vec_full.build(
      SplineType2D::BILINEAR,
      grid,
      grid,
      z0,
      false,
      false,
      SplineType2D::BILINEAR,
      grid,
      grid,
      z1,
      false,
      false );
    check_near( blend_vec_full.eval( x, y, s ), expected_value, EPS, "Spline2Dblend vector xyz build overload" );

    Spline2Dblend blend_raw_full( "blend2d_raw_full" );
    blend_raw_full.build(
      SplineType2D::BILINEAR,
      grid.data(),
      1,
      grid.data(),
      1,
      z0.data(),
      ny,
      nx,
      ny,
      false,
      false,
      SplineType2D::BILINEAR,
      grid.data(),
      1,
      grid.data(),
      1,
      z1.data(),
      ny,
      nx,
      ny,
      false,
      false );
    check_near( blend_raw_full.eval( x, y, s ), expected_value, EPS, "Spline2Dblend raw xyz build overload" );

    Spline2Dblend blend_raw_uniform( "blend2d_raw_uniform" );
    blend_raw_uniform.build(
      SplineType2D::BILINEAR,
      z0.data(),
      ny,
      nx,
      ny,
      false,
      false,
      SplineType2D::BILINEAR,
      z1.data(),
      ny,
      nx,
      ny,
      false,
      false );
    check_near( blend_raw_uniform.eval( x, y, s ), expected_value, EPS, "Spline2Dblend raw z-only build overload" );
  }

  void
  test_splineset_regressions()
  {
    std::cout << "test16: SplineSet regressions\n";

    integer const           nspl = 1;
    integer const           npts = 4;
    char const * const      headers[] = { "step" };
    SplineType1D const      types[]   = { SplineType1D::CONSTANT };
    real_type const         x[]       = { 0, 1, 2, 3 };
    real_type const         y_step[]  = { 10, 20, 30 };
    real_type const * const y_ptr[]   = { y_step };

    SplineSet raw( "raw_constant" );
    raw.build( nspl, npts, headers, types, x, y_ptr );

    check_true( raw.name() == "raw_constant", "SplineSet name" );
    check_true( raw.type() == SplineType1D::SPLINE_SET, "SplineSet type" );
    check_true( raw.num_splines() == 1, "SplineSet raw num_splines" );
    check_true( raw.num_points() == 4, "SplineSet raw num_points" );
    check_true( raw.x_nodes() != nullptr, "SplineSet x_nodes pointer" );
    check_true( raw.y_nodes( 0 ) != nullptr, "SplineSet y_nodes pointer" );
    check_true( raw.header( 0 ) == "step", "SplineSet raw header" );
    check_true( raw.get_position( "step" ) == 0, "SplineSet raw get_position" );
    check_near( raw.x_node( 2 ), 2.0, EPS, "SplineSet x_node" );
    check_near( raw.x_min(), 0.0, EPS, "SplineSet x_min" );
    check_near( raw.x_max(), 3.0, EPS, "SplineSet x_max" );
    check_near( raw.y_node( 3, 0 ), 30.0, EPS, "SplineSet raw constant tail copy" );
    check_near( raw.eval( 2.5, 0 ), 30.0, EPS, "SplineSet raw constant eval" );
    check_near( raw( 1.5, 0 ), 20.0, EPS, "SplineSet raw operator()" );
    check_near( raw.D( 1.5, 0 ), 0.0, EPS, "SplineSet raw D" );
    check_near( raw.DD( 1.5, 0 ), 0.0, EPS, "SplineSet raw DD" );
    check_near( raw.DDD( 1.5, 0 ), 0.0, EPS, "SplineSet raw DDD" );
    check_near( raw.DDDD( 1.5, 0 ), 0.0, EPS, "SplineSet raw DDDD" );
    check_near( raw.DDDDD( 1.5, 0 ), 0.0, EPS, "SplineSet raw DDDDD" );
    check_near( raw.y_min( 0 ), 10.0, EPS, "SplineSet raw y_min(index)" );
    check_near( raw.y_max( 0 ), 30.0, EPS, "SplineSet raw y_max(index)" );
    check_near( raw.y_min( "step" ), 10.0, EPS, "SplineSet raw y_min(name)" );
    check_near( raw.y_max( "step" ), 30.0, EPS, "SplineSet raw y_max(name)" );
    check_true(
      raw.is_monotone( 0 ) >= -2 && raw.is_monotone( 0 ) <= 1,
      "SplineSet raw is_monotone classification range" );
    check_true( raw.get_spline( 0 ) != nullptr, "SplineSet raw get_spline(index)" );
    check_true( raw.get_spline( "step" ) != nullptr, "SplineSet raw get_spline(name)" );
    check_true( raw.get_spline( "step" )->type() == SplineType1D::CONSTANT, "SplineSet raw spline type" );

    std::vector<std::string> got_headers;
    raw.get_headers( got_headers );
    check_true( got_headers.size() == 1 && got_headers[0] == "step", "SplineSet get_headers" );
    check_true( raw.name_list().find( "step" ) != std::string::npos, "SplineSet name_list" );

    SplineSet raw_copy( "raw_copy" );
    raw.deep_copy_to( raw_copy );
    check_near( raw_copy.eval( 2.5, 0 ), 30.0, EPS, "SplineSet deep_copy_to" );

    std::ostringstream os;
    os << raw.info() << '\n';
    raw.info( os );
    raw.dump_table( os, 4 );
    check_true( !os.str().empty(), "SplineSet info/dump_table output" );

    GenericContainer gc;
    {
      auto & spline_type = gc["spline_type"].set_vec_string();
      spline_type.emplace_back( "linear" );
      spline_type.emplace_back( "constant" );
    }
    {
      auto & hdrs = gc["headers"].set_vec_string();
      hdrs.emplace_back( "z_linear" );
      hdrs.emplace_back( "a_constant" );
    }
    gc["xdata"] = std::vector<real_type>{ 0, 1, 2, 3 };

    map_type & ymap = gc["ydata"].set_map();
    ymap["a_constant"] = std::vector<real_type>{ 10, 20, 30 };
    ymap["z_linear"]   = std::vector<real_type>{ 0, 1, 2, 3 };

    SplineSet ordered( "ordered_map" );
    ordered.build( gc );

    check_true( ordered.header( 0 ) == "z_linear", "SplineSet MAP must honor explicit headers order (header 0)" );
    check_true( ordered.header( 1 ) == "a_constant", "SplineSet MAP must honor explicit headers order (header 1)" );
    check_true( ordered.get_spline( "z_linear" )->type() == SplineType1D::LINEAR, "SplineSet MAP linear type association" );
    check_true( ordered.get_spline( "a_constant" )->type() == SplineType1D::CONSTANT, "SplineSet MAP constant type association" );
    check_near( ordered.eval( 1.5, 0 ), 1.5, EPS, "SplineSet MAP linear eval" );
    check_near( ordered.eval( 2.5, 1 ), 30.0, EPS, "SplineSet MAP constant eval" );
    check_near( ordered.y_node( 3, 1 ), 30.0, EPS, "SplineSet MAP constant tail copy" );
    check_true( ordered.x_nodes() != nullptr, "SplineSet MAP x_nodes pointer" );
    check_true( ordered.y_nodes( 0 ) != nullptr, "SplineSet MAP y_nodes pointer" );
    check_near( ordered.x_node( 2 ), 2.0, EPS, "SplineSet MAP x_node" );
    check_near( ordered.x_min(), 0.0, EPS, "SplineSet MAP x_min" );
    check_near( ordered.x_max(), 3.0, EPS, "SplineSet MAP x_max" );
    check_true( ordered.is_monotone( 0 ) > 0, "SplineSet MAP linear monotonicity" );
    check_true(
      ordered.is_monotone( 1 ) >= -2 && ordered.is_monotone( 1 ) <= 1,
      "SplineSet MAP constant monotonicity classification range" );

    auto check_all_values = [&]( map_type const & map, std::string const & label ) {
      check_true( map.size() == 2, label + " map size" );
      check_map_real( map, "z_linear", 1.5, EPS, label + " linear" );
      check_map_real( map, "a_constant", 20.0, EPS, label + " constant" );
    };

    auto check_all_derivatives = [&]( map_type const & map, std::string const & label, real_type const lin_expected ) {
      check_true( map.size() == 2, label + " map size" );
      check_map_real( map, "z_linear", lin_expected, EPS, label + " linear" );
      check_map_real( map, "a_constant", 0.0, EPS, label + " constant" );
    };

    auto check_col_values = [&]( map_type const & map, std::string const & label ) {
      check_true( map.size() == 1, label + " map size" );
      check_map_real( map, "a_constant", 20.0, EPS, label + " constant" );
    };

    auto check_all_vectors = [&]( map_type const & map, std::string const & label ) {
      check_true( map.size() == 2, label + " map size" );
      check_map_vec( map, "z_linear", { 0.5, 1.5, 2.5 }, EPS, label + " linear" );
      check_map_vec( map, "a_constant", { 10.0, 20.0, 30.0 }, EPS, label + " constant" );
    };

    auto check_all_zero_vectors = [&]( map_type const & map, std::string const & label, real_type const lin_expected ) {
      check_true( map.size() == 2, label + " map size" );
      check_map_vec( map, "z_linear", { lin_expected, lin_expected, lin_expected }, EPS, label + " linear" );
      check_map_vec( map, "a_constant", { 0.0, 0.0, 0.0 }, EPS, label + " constant" );
    };

    auto check_col_vectors = [&]( map_type const & map, std::string const & label ) {
      check_true( map.size() == 1, label + " map size" );
      check_map_vec( map, "a_constant", { 10.0, 20.0, 30.0 }, EPS, label + " constant" );
    };

    vec_string_type const cols  = { "a_constant" };
    vec_real_type const   xs    = { 0.5, 1.5, 2.5 };
    real_type const       zeta  = 1.5;
    real_type             xhit  = -1;
    GenericContainer      gc_out;
    auto clear_gc = [&]( GenericContainer & gc ) { gc.set_map().clear(); };

    clear_gc( gc_out );
    ordered.eval( 1.5, gc_out );
    check_all_values( gc_out.get_map(), "SplineSet eval(gc) all" );
    clear_gc( gc_out );
    ordered.eval_D( 1.5, gc_out );
    check_all_derivatives( gc_out.get_map(), "SplineSet eval_D(gc) all", 1.0 );
    clear_gc( gc_out );
    ordered.eval_DD( 1.5, gc_out );
    check_all_derivatives( gc_out.get_map(), "SplineSet eval_DD(gc) all", 0.0 );
    clear_gc( gc_out );
    ordered.eval_DDD( 1.5, gc_out );
    check_all_derivatives( gc_out.get_map(), "SplineSet eval_DDD(gc) all", 0.0 );

    clear_gc( gc_out );
    ordered.eval( 1.5, cols, gc_out );
    check_col_values( gc_out.get_map(), "SplineSet eval(gc) columns" );
    clear_gc( gc_out );
    ordered.eval_D( 1.5, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval_D(gc) columns map size" );
    check_map_real( gc_out.get_map(), "a_constant", 0.0, EPS, "SplineSet eval_D(gc) columns" );
    clear_gc( gc_out );
    ordered.eval_DD( 1.5, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval_DD(gc) columns map size" );
    check_map_real( gc_out.get_map(), "a_constant", 0.0, EPS, "SplineSet eval_DD(gc) columns" );
    clear_gc( gc_out );
    ordered.eval_DDD( 1.5, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval_DDD(gc) columns map size" );
    check_map_real( gc_out.get_map(), "a_constant", 0.0, EPS, "SplineSet eval_DDD(gc) columns" );

    clear_gc( gc_out );
    ordered.eval( xs, gc_out );
    check_all_vectors( gc_out.get_map(), "SplineSet eval(vec,gc) all" );
    clear_gc( gc_out );
    ordered.eval_D( xs, gc_out );
    check_all_zero_vectors( gc_out.get_map(), "SplineSet eval_D(vec,gc) all", 1.0 );
    clear_gc( gc_out );
    ordered.eval_DD( xs, gc_out );
    check_all_zero_vectors( gc_out.get_map(), "SplineSet eval_DD(vec,gc) all", 0.0 );
    clear_gc( gc_out );
    ordered.eval_DDD( xs, gc_out );
    check_all_zero_vectors( gc_out.get_map(), "SplineSet eval_DDD(vec,gc) all", 0.0 );

    clear_gc( gc_out );
    ordered.eval( xs, cols, gc_out );
    check_col_vectors( gc_out.get_map(), "SplineSet eval(vec,cols,gc)" );
    clear_gc( gc_out );
    ordered.eval_D( xs, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval_D(vec,cols,gc) map size" );
    check_map_vec( gc_out.get_map(), "a_constant", { 0.0, 0.0, 0.0 }, EPS, "SplineSet eval_D(vec,cols,gc)" );
    clear_gc( gc_out );
    ordered.eval_DD( xs, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval_DD(vec,cols,gc) map size" );
    check_map_vec( gc_out.get_map(), "a_constant", { 0.0, 0.0, 0.0 }, EPS, "SplineSet eval_DD(vec,cols,gc)" );
    clear_gc( gc_out );
    ordered.eval_DDD( xs, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval_DDD(vec,cols,gc) map size" );
    check_map_vec( gc_out.get_map(), "a_constant", { 0.0, 0.0, 0.0 }, EPS, "SplineSet eval_DDD(vec,cols,gc)" );

    clear_gc( gc_out );
    ordered.eval2( zeta, 0, xhit, gc_out );
    check_near( xhit, 1.5, EPS, "SplineSet eval2 index xhit" );
    check_all_values( gc_out.get_map(), "SplineSet eval2(gc) index" );
    clear_gc( gc_out );
    ordered.eval2_D( zeta, 0, xhit, gc_out );
    check_all_derivatives( gc_out.get_map(), "SplineSet eval2_D(gc) index", 1.0 );
    clear_gc( gc_out );
    ordered.eval2_DD( zeta, 0, xhit, gc_out );
    check_all_derivatives( gc_out.get_map(), "SplineSet eval2_DD(gc) index", 0.0 );
    clear_gc( gc_out );
    ordered.eval2_DDD( zeta, 0, xhit, gc_out );
    check_all_derivatives( gc_out.get_map(), "SplineSet eval2_DDD(gc) index", 0.0 );

    clear_gc( gc_out );
    ordered.eval2( zeta, 0, xhit, cols, gc_out );
    check_near( xhit, 1.5, EPS, "SplineSet eval2 columns index xhit" );
    check_col_values( gc_out.get_map(), "SplineSet eval2(gc) index columns" );
    clear_gc( gc_out );
    ordered.eval2_D( zeta, 0, xhit, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_D(gc) index columns map size" );
    check_map_real( gc_out.get_map(), "a_constant", 0.0, EPS, "SplineSet eval2_D(gc) index columns" );
    clear_gc( gc_out );
    ordered.eval2_DD( zeta, 0, xhit, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_DD(gc) index columns map size" );
    check_map_real( gc_out.get_map(), "a_constant", 0.0, EPS, "SplineSet eval2_DD(gc) index columns" );
    clear_gc( gc_out );
    ordered.eval2_DDD( zeta, 0, xhit, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_DDD(gc) index columns map size" );
    check_map_real( gc_out.get_map(), "a_constant", 0.0, EPS, "SplineSet eval2_DDD(gc) index columns" );

    clear_gc( gc_out );
    ordered.eval2( xs, 0, gc_out );
    check_all_vectors( gc_out.get_map(), "SplineSet eval2(vec,gc) index" );
    clear_gc( gc_out );
    ordered.eval2_D( xs, 0, gc_out );
    check_all_zero_vectors( gc_out.get_map(), "SplineSet eval2_D(vec,gc) index", 1.0 );
    clear_gc( gc_out );
    ordered.eval2_DD( xs, 0, gc_out );
    check_all_zero_vectors( gc_out.get_map(), "SplineSet eval2_DD(vec,gc) index", 0.0 );
    clear_gc( gc_out );
    ordered.eval2_DDD( xs, 0, gc_out );
    check_all_zero_vectors( gc_out.get_map(), "SplineSet eval2_DDD(vec,gc) index", 0.0 );

    clear_gc( gc_out );
    ordered.eval2( xs, 0, cols, gc_out );
    check_col_vectors( gc_out.get_map(), "SplineSet eval2(vec,cols,gc) index" );
    clear_gc( gc_out );
    ordered.eval2_D( xs, 0, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_D(vec,cols,gc) index map size" );
    check_map_vec( gc_out.get_map(), "a_constant", { 0.0, 0.0, 0.0 }, EPS, "SplineSet eval2_D(vec,cols,gc) index" );
    clear_gc( gc_out );
    ordered.eval2_DD( xs, 0, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_DD(vec,cols,gc) index map size" );
    check_map_vec( gc_out.get_map(), "a_constant", { 0.0, 0.0, 0.0 }, EPS, "SplineSet eval2_DD(vec,cols,gc) index" );
    clear_gc( gc_out );
    ordered.eval2_DDD( xs, 0, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_DDD(vec,cols,gc) index map size" );
    check_map_vec( gc_out.get_map(), "a_constant", { 0.0, 0.0, 0.0 }, EPS, "SplineSet eval2_DDD(vec,cols,gc) index" );

    clear_gc( gc_out );
    ordered.eval2( zeta, "z_linear", xhit, gc_out );
    check_near( xhit, 1.5, EPS, "SplineSet eval2 name xhit" );
    check_all_values( gc_out.get_map(), "SplineSet eval2(gc) name" );
    clear_gc( gc_out );
    ordered.eval2_D( zeta, "z_linear", xhit, gc_out );
    check_all_derivatives( gc_out.get_map(), "SplineSet eval2_D(gc) name", 1.0 );
    clear_gc( gc_out );
    ordered.eval2_DD( zeta, "z_linear", xhit, gc_out );
    check_all_derivatives( gc_out.get_map(), "SplineSet eval2_DD(gc) name", 0.0 );
    clear_gc( gc_out );
    ordered.eval2_DDD( zeta, "z_linear", xhit, gc_out );
    check_all_derivatives( gc_out.get_map(), "SplineSet eval2_DDD(gc) name", 0.0 );

    clear_gc( gc_out );
    ordered.eval2( zeta, "z_linear", xhit, cols, gc_out );
    check_near( xhit, 1.5, EPS, "SplineSet eval2 name columns xhit" );
    check_col_values( gc_out.get_map(), "SplineSet eval2(gc) name columns" );
    clear_gc( gc_out );
    ordered.eval2_D( zeta, "z_linear", xhit, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_D(gc) name columns map size" );
    check_map_real( gc_out.get_map(), "a_constant", 0.0, EPS, "SplineSet eval2_D(gc) name columns" );
    clear_gc( gc_out );
    ordered.eval2_DD( zeta, "z_linear", xhit, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_DD(gc) name columns map size" );
    check_map_real( gc_out.get_map(), "a_constant", 0.0, EPS, "SplineSet eval2_DD(gc) name columns" );
    clear_gc( gc_out );
    ordered.eval2_DDD( zeta, "z_linear", xhit, cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_DDD(gc) name columns map size" );
    check_map_real( gc_out.get_map(), "a_constant", 0.0, EPS, "SplineSet eval2_DDD(gc) name columns" );

    clear_gc( gc_out );
    ordered.eval2( xs, "z_linear", gc_out );
    check_all_vectors( gc_out.get_map(), "SplineSet eval2(vec,gc) name" );
    clear_gc( gc_out );
    ordered.eval2_D( xs, "z_linear", gc_out );
    check_all_zero_vectors( gc_out.get_map(), "SplineSet eval2_D(vec,gc) name", 1.0 );
    clear_gc( gc_out );
    ordered.eval2_DD( xs, "z_linear", gc_out );
    check_all_zero_vectors( gc_out.get_map(), "SplineSet eval2_DD(vec,gc) name", 0.0 );
    clear_gc( gc_out );
    ordered.eval2_DDD( xs, "z_linear", gc_out );
    check_all_zero_vectors( gc_out.get_map(), "SplineSet eval2_DDD(vec,gc) name", 0.0 );

    clear_gc( gc_out );
    ordered.eval2( xs, "z_linear", cols, gc_out );
    check_col_vectors( gc_out.get_map(), "SplineSet eval2(vec,cols,gc) name" );
    clear_gc( gc_out );
    ordered.eval2_D( xs, "z_linear", cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_D(vec,cols,gc) name map size" );
    check_map_vec( gc_out.get_map(), "a_constant", { 0.0, 0.0, 0.0 }, EPS, "SplineSet eval2_D(vec,cols,gc) name" );
    clear_gc( gc_out );
    ordered.eval2_DD( xs, "z_linear", cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_DD(vec,cols,gc) name map size" );
    check_map_vec( gc_out.get_map(), "a_constant", { 0.0, 0.0, 0.0 }, EPS, "SplineSet eval2_DD(vec,cols,gc) name" );
    clear_gc( gc_out );
    ordered.eval2_DDD( xs, "z_linear", cols, gc_out );
    check_true( gc_out.get_map().size() == 1, "SplineSet eval2_DDD(vec,cols,gc) name map size" );
    check_map_vec( gc_out.get_map(), "a_constant", { 0.0, 0.0, 0.0 }, EPS, "SplineSet eval2_DDD(vec,cols,gc) name" );
  }

  void
  test_splinevec_api()
  {
    std::cout << "test16: SplineVec API smoke\n";

    integer const dim  = 2;
    integer const npts = 4;

    real_type const col_major[] = {
      0, 0,
      1, 1,
      2, 4,
      3, 9
    };
    real_type const knots[] = { 0, 1, 2, 3 };

    SplineVec curve( "curve_matrix" );
    curve.setup( dim, npts, col_major, dim );
    curve.set_knots( knots );
    curve.catmull_rom();

    check_true( curve.name() == "curve_matrix", "SplineVec name" );
    check_true( curve.dimension() == dim, "SplineVec dimension" );
    check_true( curve.num_points() == npts, "SplineVec num_points" );
    check_true( curve.x_nodes() != nullptr, "SplineVec x_nodes pointer" );
    check_true( curve.y_nodes( 0 ) != nullptr, "SplineVec y_nodes pointer" );
    check_near( curve.x_node( 1 ), 1.0, EPS, "SplineVec x_node" );
    check_near( curve.y_node( 2, 1 ), 4.0, EPS, "SplineVec y_node" );
    check_near( curve.x_min(), 0.0, EPS, "SplineVec x_min" );
    check_near( curve.x_max(), 3.0, EPS, "SplineVec x_max" );

    curve.make_bounded();
    check_true( !curve.can_extend(), "SplineVec make_bounded" );
    curve.make_unbounded();
    check_true( curve.can_extend(), "SplineVec make_unbounded" );
    curve.make_closed();
    check_true( curve.is_closed(), "SplineVec make_closed" );
    curve.make_open();
    check_true( !curve.is_closed(), "SplineVec make_open" );

    real_type const t = 1.5;
    check_near( curve.eval( t, 0 ), curve( t, 0 ), EPS, "SplineVec operator()(x,i)" );
    check_true( std::isfinite( curve.curvature( t ) ), "SplineVec curvature" );
    check_true( std::isfinite( curve.curvature_D( t ) ), "SplineVec curvature_D" );
    {
      real_type const h  = 1e-6;
      real_type const fd = ( curve.curvature( t + h ) - curve.curvature( t - h ) ) / ( 2 * h );
      check_near( curve.curvature_D( t ), fd, 1e-4, "SplineVec curvature_D finite-difference consistency" );
    }

    real_type vals[2], d1[2], d2[2], d3[2], d4[2], d5[2];
    curve.eval( t, vals, 1 );
    curve.eval_D( t, d1, 1 );
    curve.eval_DD( t, d2, 1 );
    curve.eval_DDD( t, d3, 1 );
    curve.eval_DDDD( t, d4, 1 );
    curve.eval_DDDDD( t, d5, 1 );

    check_near( vals[0], curve.eval( t, 0 ), EPS, "SplineVec eval array component 0" );
    check_near( vals[1], curve.eval( t, 1 ), EPS, "SplineVec eval array component 1" );
    check_near( d1[0], curve.eval_D( t, 0 ), EPS_D, "SplineVec eval_D array component 0" );
    check_near( d1[1], curve.eval_D( t, 1 ), EPS_D, "SplineVec eval_D array component 1" );
    check_near( d2[0], curve.eval_DD( t, 0 ), EPS_D, "SplineVec eval_DD array component 0" );
    check_near( d3[1], curve.eval_DDD( t, 1 ), EPS_D, "SplineVec eval_DDD array component 1" );
    check_near( d4[0], 0.0, EPS, "SplineVec eval_DDDD array component 0" );
    check_near( d5[1], 0.0, EPS, "SplineVec eval_DDDDD array component 1" );

    std::vector<real_type> vv;
    curve.eval( t, vv );
    curve.eval_D( t, vv );
    curve.eval_DD( t, vv );
    curve.eval_DDD( t, vv );
    curve.eval_DDDD( t, vv );
    curve.eval_DDDDD( t, vv );
    check_true( vv.size() == static_cast<size_t>( dim ), "SplineVec vector evaluation size" );

    GenericContainer gc_vals;
    curve.eval( t, gc_vals );
    curve.eval_D( t, gc_vals );
    curve.eval_DD( t, gc_vals );
    curve.eval_DDD( t, gc_vals );
    curve.eval_DDDD( t, gc_vals );
    curve.eval_DDDDD( t, gc_vals );
    check_true( gc_vals.get_vec_real().size() == static_cast<size_t>( dim ), "SplineVec GC evaluation size" );

    vec_real_type ts = { 0.0, 1.0, 2.0 };
    GenericContainer batch;
    curve.eval( ts, batch );
    check_true(
      batch.get_mat_real().num_rows() == static_cast<size_t>( dim ) &&
      batch.get_mat_real().num_cols() == ts.size(),
      "SplineVec batch eval matrix shape" );
    check_near( batch.get_mat_real()( 0, 1 ), curve.eval( ts[1], 0 ), EPS, "SplineVec batch eval value" );
    check_near( batch.get_mat_real()( 1, 2 ), curve.eval( ts[2], 1 ), EPS, "SplineVec batch eval value second component" );
    curve.eval_D( ts, batch );
    check_near( batch.get_mat_real()( 0, 1 ), curve.eval_D( ts[1], 0 ), EPS_D, "SplineVec batch eval_D value" );
    curve.eval_DD( ts, batch );
    check_near( batch.get_mat_real()( 1, 2 ), curve.eval_DD( ts[2], 1 ), EPS_D, "SplineVec batch eval_DD value" );
    curve.eval_DDD( ts, batch );
    check_near( batch.get_mat_real()( 0, 0 ), curve.eval_DDD( ts[0], 0 ), EPS_D, "SplineVec batch eval_DDD value" );
    curve.eval_DDDD( ts, batch );
    check_near( batch.get_mat_real()( 1, 1 ), 0.0, EPS, "SplineVec batch eval_DDDD value" );
    curve.eval_DDDDD( ts, batch );
    check_near( batch.get_mat_real()( 0, 2 ), 0.0, EPS, "SplineVec batch eval_DDDDD value" );

    SplineVec copy( "curve_copy" );
    curve.deep_copy_to( copy );
    copy.set_knots( knots );
    copy.catmull_rom();
    check_near( copy.eval( t, 0 ), curve.eval( t, 0 ), EPS_D, "SplineVec deep_copy_to eval component 0" );
    check_near( copy.eval( t, 1 ), curve.eval( t, 1 ), EPS_D, "SplineVec deep_copy_to eval component 1" );

    real_type const comp0[] = { 0, 1, 2, 3 };
    real_type const comp1[] = { 0, 1, 4, 9 };
    real_type const * row_ptrs[] = { comp0, comp1 };
    SplineVec         curve_ptr( "curve_ptr" );
    curve_ptr.setup( dim, npts, row_ptrs );
    curve_ptr.set_knots( knots );
    curve_ptr.catmull_rom();
    check_near( curve_ptr.eval( t, 0 ), curve.eval( t, 0 ), EPS_D, "SplineVec setup(pointer array) component 0" );
    check_near( curve_ptr.eval( t, 1 ), curve.eval( t, 1 ), EPS_D, "SplineVec setup(pointer array) component 1" );

    GenericContainer gc_curve;
    auto &           data = gc_curve["data"].set_mat_real( npts, dim );
    data( 0, 0 )         = 0;
    data( 0, 1 )         = 0;
    data( 1, 0 )         = 1;
    data( 1, 1 )         = 1;
    data( 2, 0 )         = 2;
    data( 2, 1 )         = 4;
    data( 3, 0 )         = 3;
    data( 3, 1 )         = 9;
    gc_curve["transposed"] = false;

    SplineVec curve_gc( "curve_gc" );
    curve_gc.build( gc_curve );
    curve_gc.set_knots_chord_length();
    curve_gc.catmull_rom();
    check_true( curve_gc.dimension() == dim, "SplineVec setup(gc) dimension" );

    SplineVec curve_cent( "curve_cent" );
    curve_cent.setup( dim, npts, col_major, dim );
    curve_cent.set_knots_centripetal();
    curve_cent.catmull_rom();
    check_true( curve_cent.x_max() > curve_cent.x_min(), "SplineVec set_knots_centripetal" );

    SplineVec curve_foley( "curve_foley" );
    curve_foley.setup( dim, npts, col_major, dim );
    curve_foley.set_knots_foley();
    curve_foley.catmull_rom();
    check_true( curve_foley.x_max() > curve_foley.x_min(), "SplineVec set_knots_foley" );
  }
}  // namespace

int
main()
{
  try
  {
    test_c_interface();
    test_spline1dblend_api();
    test_spline2d_and_surface_setup();
    test_spline2dblend_api();
    test_splineset_regressions();
    test_splinevec_api();
  }
  catch ( std::exception const & exc )
  {
    std::cerr << "Unhandled exception in test16: " << exc.what() << '\n';
    return EXIT_KO;
  }

  if ( failures == 0 )
  {
    std::cout << "test16: all regression and API smoke checks passed\n";
    return EXIT_OK;
  }

  std::cerr << "test16: " << failures << " checks failed\n";
  return EXIT_KO;
}
