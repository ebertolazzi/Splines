/*
 * Test export_csv/export_json/export_yaml for cubic and quintic spline bases.
 */

#include "Splines/Splines.hh"
#include "GenericContainer/GenericContainerInterface_nlohmann.hh"

#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace Splines;

namespace
{
  int failures = 0;

  class TestCubicExportSpline : public CubicSplineBase
  {
  public:
    explicit TestCubicExportSpline( string_view name = "TestCubicExportSpline" ) : CubicSplineBase( name ) {}

    SplineType1D type() const override { return SplineType1D::HERMITE; }
    void         build() override { m_search.must_reset(); }
    void         setup( GenericContainer const & ) override {}

    void
    set_data( std::vector<real_type> const & x, std::vector<real_type> const & y, std::vector<real_type> const & yp )
    {
      integer const n = static_cast<integer>( x.size() );
      reserve( n );
      for ( integer i = 0; i < n; ++i )
      {
        m_X[i]  = x[static_cast<size_t>( i )];
        m_Y[i]  = y[static_cast<size_t>( i )];
        m_Yp[i] = yp[static_cast<size_t>( i )];
      }
      m_npts = n;
      m_search.must_reset();
    }
  };

  class TestQuinticExportSpline : public QuinticSplineBase
  {
  public:
    explicit TestQuinticExportSpline( string_view name = "TestQuinticExportSpline" )
    : QuinticSplineBase( Spline_sub_type::CUBIC, name ) {}

    void build() override { m_search.must_reset(); }
    void setup( GenericContainer const & ) override {}

    void
    set_data(
      std::vector<real_type> const & x,
      std::vector<real_type> const & y,
      std::vector<real_type> const & yp,
      std::vector<real_type> const & ypp )
    {
      integer const n = static_cast<integer>( x.size() );
      reserve( n );
      for ( integer i = 0; i < n; ++i )
      {
        m_X[i]   = x[static_cast<size_t>( i )];
        m_Y[i]   = y[static_cast<size_t>( i )];
        m_Yp[i]  = yp[static_cast<size_t>( i )];
        m_Ypp[i] = ypp[static_cast<size_t>( i )];
      }
      m_npts = n;
      m_search.must_reset();
    }
  };

  void
  check_true( bool const cond, std::string const & msg )
  {
    if ( cond ) return;
    ++failures;
    std::cerr << "FAIL: " << msg << '\n';
  }

  bool
  parse_json( std::string const & text, GenericContainer & gc )
  {
    try
    {
      gc = nlohmann::json::parse( text ).get<GenericContainer>();
      return true;
    }
    catch ( std::exception const & e )
    {
      std::cout << "[FAIL] JSON parse: " << e.what() << "\n";
      return false;
    }
  }

  void
  check_near( real_type got, real_type expected, real_type tol, std::string const & msg )
  {
    if ( std::abs( got - expected ) <= tol ) return;
    ++failures;
    std::cerr << "FAIL: " << msg << " got=" << got << " expected=" << expected << " tol=" << tol << '\n';
  }

  void
  check_bool( GenericContainer const & gc, bool expected, std::string const & msg )
  {
    check_true( gc.get_bool( msg ) == expected, msg );
  }

  real_type
  gc_as_real( GenericContainer const & gc, std::string const & where )
  {
    switch ( gc.get_type() )
    {
      case GC_type::INTEGER: return static_cast<real_type>( gc.get_int( where ) );
      case GC_type::LONG: return static_cast<real_type>( gc.get_long( where ) );
      case GC_type::REAL: return gc.get_real( where );
      default:
        throw std::runtime_error(
          "Expected numeric value at `" + where + "`, found `" + std::string( gc.get_type_name() ) + "`" );
    }
  }

  std::string
  read_file( std::string const & fname )
  {
    std::ifstream file( fname );
    std::ostringstream ss;
    ss << file.rdbuf();
    return ss.str();
  }

  std::vector<std::vector<real_type>>
  parse_csv( std::string const & text, std::string const & expected_header )
  {
    std::istringstream                     input( text );
    std::string                            line;
    std::vector<std::vector<real_type>>    rows;

    check_true( static_cast<bool>( std::getline( input, line ) ), "CSV must contain header" );
    if ( line != expected_header )
      check_true( false, "Unexpected CSV header: `" + line + "` expected `" + expected_header + "`" );

    while ( std::getline( input, line ) )
    {
      if ( line.empty() ) continue;
      std::vector<real_type> row;
      std::istringstream     row_stream( line );
      std::string            cell;
      while ( std::getline( row_stream, cell, ',' ) ) row.push_back( std::stod( cell ) );
      rows.emplace_back( std::move( row ) );
    }
    return rows;
  }

  real_type
  eval_poly( std::vector<real_type> const & coeffs, real_type x0, real_type x )
  {
    real_type const t = x - x0;
    real_type       y = 0;
    real_type       p = 1;
    for ( real_type const c : coeffs )
    {
      y += c * p;
      p *= t;
    }
    return y;
  }

  void
  check_segment_map(
    GenericContainer const & seg,
    std::vector<real_type> const & expected,
    real_type const x0,
    real_type const x1,
    integer const idx,
    std::string const & label )
  {
    auto const & m = seg.get_map();
    check_true( m.at( "index" ).get_int() == idx, label + " index" );
    check_near( gc_as_real( m.at( "x" ), label + ".x" ), x0, 1e-14, label + " x" );
    check_near( gc_as_real( m.at( "x_next" ), label + ".x_next" ), x1, 1e-14, label + " x_next" );
    check_near( gc_as_real( m.at( "h" ), label + ".h" ), x1 - x0, 1e-14, label + " h" );

    static char const * const names[] = { "a", "b", "c", "d", "e", "f" };
    for ( size_t i = 0; i < expected.size(); ++i )
      check_near( gc_as_real( m.at( names[i] ), label + "." + names[i] ), expected[i], 1e-12, label + " coeff " + names[i] );
  }

  void
  test_cubic_export()
  {
    std::cout << "test17: cubic export\n";

    std::vector<real_type> const x  = { 0, 1, 3 };
    std::vector<real_type> const y  = { 1, 10, 54 };
    std::vector<real_type> const yp = { 2, 20, 16 };

    std::vector<std::vector<real_type>> const expected = {
      { 1, 2, 3, 4 },
      { 10, 20, 5, -2 }
    };

    TestCubicExportSpline spline( "cubic_export" );
    spline.set_data( x, y, yp );

    std::ostringstream csv, json, yaml;
    spline.export_csv( csv );
    spline.export_json( json );
    spline.export_yaml( yaml );

    auto const rows = parse_csv( csv.str(), "x,a,b,c,d" );
    check_true( rows.size() == expected.size() + 1, "Cubic CSV row count" );
    for ( size_t i = 0; i < expected.size() && i < rows.size(); ++i )
    {
      check_true( rows[i].size() == 5, "Cubic CSV column count row " + std::to_string( i ) );
      if ( rows[i].size() != 5 ) continue;
      check_near( rows[i][0], x[i], 1e-14, "Cubic CSV x row " + std::to_string( i ) );
      for ( size_t j = 0; j < 4; ++j )
        check_near( rows[i][j + 1], expected[i][j], 1e-12, "Cubic CSV coeff row " + std::to_string( i ) );
    }
    if ( rows.size() == expected.size() + 1 )
    {
      auto const & row = rows.back();
      check_true( row.size() == 5, "Cubic CSV final dummy row column count" );
      if ( row.size() == 5 )
      {
        check_near( row[0], x.back(), 1e-14, "Cubic CSV final dummy row x" );
        check_near( row[1], y.back(), 1e-14, "Cubic CSV final dummy row a" );
        check_near( row[2], 0, 1e-14, "Cubic CSV final dummy row b" );
        check_near( row[3], 0, 1e-14, "Cubic CSV final dummy row c" );
        check_near( row[4], 0, 1e-14, "Cubic CSV final dummy row d" );
      }
    }

    GenericContainer gc_json, gc_yaml;
    check_true( parse_json( json.str(), gc_json ), "Cubic JSON parse" );
    check_true( gc_yaml.from_yaml( yaml.str() ), "Cubic YAML parse" );
    check_true( nlohmann::json( gc_json ) == nlohmann::json( gc_yaml ), "Cubic JSON/YAML content equality" );

    auto const & map = gc_json.get_map();
    check_true( map.at( "name" ).get_string() == "cubic_export", "Cubic JSON name" );
    check_true( map.at( "spline_type" ).get_string() == "SPLINE_HERMITE", "Cubic JSON type" );
    check_true( map.at( "polynomial_order" ).get_int() == 3, "Cubic JSON polynomial_order" );
    check_true( map.at( "num_points" ).get_int() == 3, "Cubic JSON num_points" );
    check_true( map.at( "num_segments" ).get_int() == 2, "Cubic JSON num_segments" );
    check_near( gc_as_real( map.at( "xmin" ), "cubic.xmin" ), x.front(), 1e-14, "Cubic JSON xmin" );
    check_near( gc_as_real( map.at( "xmax" ), "cubic.xmax" ), x.back(), 1e-14, "Cubic JSON xmax" );
    check_bool( map.at( "closed" ), false, "Cubic JSON closed" );
    check_bool( map.at( "cyclic" ), false, "Cubic JSON cyclic" );
    check_bool( map.at( "bounded" ), false, "Cubic JSON bounded" );
    check_bool( map.at( "can_extend" ), true, "Cubic JSON can_extend" );
    check_bool( map.at( "extended_constant" ), false, "Cubic JSON extended_constant" );

    auto const & segments = map.at( "segments" ).get_vector();
    check_true( segments.size() == 2, "Cubic JSON segments size" );
    if ( segments.size() == 2 )
    {
      check_segment_map( segments[0], expected[0], x[0], x[1], 0, "Cubic segment 0" );
      check_segment_map( segments[1], expected[1], x[1], x[2], 1, "Cubic segment 1" );
    }

    for ( size_t i = 0; i < expected.size(); ++i )
    {
      real_type const x0 = x[i];
      real_type const x1 = x[i + 1];
      real_type const xm = 0.5 * ( x0 + x1 );
      check_near( eval_poly( expected[i], x0, x0 ), spline.eval( x0 ), 1e-12, "Cubic polynomial at left node" );
      check_near( eval_poly( expected[i], x0, xm ), spline.eval( xm ), 1e-12, "Cubic polynomial at midpoint" );
      check_near( eval_poly( expected[i], x0, x1 ), spline.eval( x1 ), 1e-12, "Cubic polynomial at right node" );
    }

    std::string const csv_file  = "test17_cubic.csv";
    std::string const json_file = "test17_cubic.json";
    std::string const yaml_file = "test17_cubic.yaml";
    spline.export_csv( csv_file );
    spline.export_json( json_file );
    spline.export_yaml( yaml_file );
    check_true( read_file( csv_file ) == csv.str(), "Cubic CSV file overload" );
    check_true( read_file( json_file ) == json.str(), "Cubic JSON file overload" );
    check_true( read_file( yaml_file ) == yaml.str(), "Cubic YAML file overload" );
    //std::remove( csv_file.c_str() );
    //std::remove( json_file.c_str() );
    //std::remove( yaml_file.c_str() );
  }

  void
  test_quintic_export()
  {
    std::cout << "test17: quintic export\n";

    std::vector<real_type> const x   = { 0, 1, 3 };
    std::vector<real_type> const y   = { 1, 21, 557 };
    std::vector<real_type> const yp  = { 2, 70, 438 };
    std::vector<real_type> const ypp = { 6, 210, 110 };

    std::vector<std::vector<real_type>> const expected = {
      { 1, 2, 3, 4, 5, 6 },
      { 21, 70, 105, -3, 2, -1 }
    };

    TestQuinticExportSpline spline( "quintic_export" );
    spline.set_data( x, y, yp, ypp );
    spline.make_closed();
    spline.make_extended_constant();

    std::ostringstream csv, json, yaml;
    spline.export_csv( csv );
    spline.export_json( json );
    spline.export_yaml( yaml );

    auto const rows = parse_csv( csv.str(), "x,a,b,c,d,e,f" );
    check_true( rows.size() == expected.size() + 1, "Quintic CSV row count" );
    for ( size_t i = 0; i < expected.size() && i < rows.size(); ++i )
    {
      check_true( rows[i].size() == 7, "Quintic CSV column count row " + std::to_string( i ) );
      if ( rows[i].size() != 7 ) continue;
      check_near( rows[i][0], x[i], 1e-14, "Quintic CSV x row " + std::to_string( i ) );
      for ( size_t j = 0; j < 6; ++j )
        check_near( rows[i][j + 1], expected[i][j], 1e-11, "Quintic CSV coeff row " + std::to_string( i ) );
    }
    if ( rows.size() == expected.size() + 1 )
    {
      auto const & row = rows.back();
      check_true( row.size() == 7, "Quintic CSV final dummy row column count" );
      if ( row.size() == 7 )
      {
        check_near( row[0], x.back(), 1e-14, "Quintic CSV final dummy row x" );
        check_near( row[1], y.back(), 1e-14, "Quintic CSV final dummy row a" );
        check_near( row[2], 0, 1e-14, "Quintic CSV final dummy row b" );
        check_near( row[3], 0, 1e-14, "Quintic CSV final dummy row c" );
        check_near( row[4], 0, 1e-14, "Quintic CSV final dummy row d" );
        check_near( row[5], 0, 1e-14, "Quintic CSV final dummy row e" );
        check_near( row[6], 0, 1e-14, "Quintic CSV final dummy row f" );
      }
    }

    GenericContainer gc_json, gc_yaml;
    check_true( parse_json( json.str(), gc_json ), "Quintic JSON parse" );
    check_true( gc_yaml.from_yaml( yaml.str() ), "Quintic YAML parse" );
    check_true( nlohmann::json( gc_json ) == nlohmann::json( gc_yaml ), "Quintic JSON/YAML content equality" );

    auto const & map = gc_json.get_map();
    check_true( map.at( "name" ).get_string() == "quintic_export", "Quintic JSON name" );
    check_true( map.at( "spline_type" ).get_string() == "SPLINE_QUINTIC", "Quintic JSON type" );
    check_true( map.at( "polynomial_order" ).get_int() == 5, "Quintic JSON polynomial_order" );
    check_true( map.at( "num_points" ).get_int() == 3, "Quintic JSON num_points" );
    check_true( map.at( "num_segments" ).get_int() == 2, "Quintic JSON num_segments" );
    check_near( gc_as_real( map.at( "xmin" ), "quintic.xmin" ), x.front(), 1e-14, "Quintic JSON xmin" );
    check_near( gc_as_real( map.at( "xmax" ), "quintic.xmax" ), x.back(), 1e-14, "Quintic JSON xmax" );
    check_bool( map.at( "closed" ), true, "Quintic JSON closed" );
    check_bool( map.at( "cyclic" ), true, "Quintic JSON cyclic" );
    check_bool( map.at( "bounded" ), false, "Quintic JSON bounded" );
    check_bool( map.at( "can_extend" ), true, "Quintic JSON can_extend" );
    check_bool( map.at( "extended_constant" ), true, "Quintic JSON extended_constant" );

    auto const & segments = map.at( "segments" ).get_vector();
    check_true( segments.size() == 2, "Quintic JSON segments size" );
    if ( segments.size() == 2 )
    {
      check_segment_map( segments[0], expected[0], x[0], x[1], 0, "Quintic segment 0" );
      check_segment_map( segments[1], expected[1], x[1], x[2], 1, "Quintic segment 1" );
    }

    for ( size_t i = 0; i < expected.size(); ++i )
    {
      real_type const x0  = x[i];
      real_type const x1  = x[i + 1];
      real_type const xq1 = x0 + ( x1 - x0 ) / 3;
      real_type const xq2 = x0 + 2 * ( x1 - x0 ) / 3;
      check_near( eval_poly( expected[i], x0, x0 ), spline.eval( x0 ), 1e-11, "Quintic polynomial at left node" );
      check_near( eval_poly( expected[i], x0, xq1 ), spline.eval( xq1 ), 1e-10, "Quintic polynomial at first third" );
      check_near( eval_poly( expected[i], x0, xq2 ), spline.eval( xq2 ), 1e-10, "Quintic polynomial at second third" );
      check_near( eval_poly( expected[i], x0, x1 ), spline.eval( x1 ), 1e-10, "Quintic polynomial at right node" );
    }

    std::string const csv_file  = "test17_quintic.csv";
    std::string const json_file = "test17_quintic.json";
    std::string const yaml_file = "test17_quintic.yaml";
    spline.export_csv( csv_file );
    spline.export_json( json_file );
    spline.export_yaml( yaml_file );
    check_true( read_file( csv_file ) == csv.str(), "Quintic CSV file overload" );
    check_true( read_file( json_file ) == json.str(), "Quintic JSON file overload" );
    check_true( read_file( yaml_file ) == yaml.str(), "Quintic YAML file overload" );
    //std::remove( csv_file.c_str() );
    //std::remove( json_file.c_str() );
    //std::remove( yaml_file.c_str() );
  }
}  // namespace

int
main()
{
  test_cubic_export();
  test_quintic_export();

  if ( failures == 0 )
  {
    std::cout << "test17: all export checks passed\n";
    return 0;
  }

  std::cerr << "test17: failures = " << failures << '\n';
  return 1;
}
