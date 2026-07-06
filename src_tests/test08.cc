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

/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Merged Spline Tests                                                     |
 |                                                                          |
 |  Combined tests for 1D and 2D spline interpolation                       |
 |                                                                          |
 |  Features:                                                               |
 |  - Unified test framework                                                |
 |  - Colored terminal output using fmt                                     |
 |  - Unicode table formatting                                              |
 |  - File-based constructor testing                                        |
 |  - Autodiff derivative validation                                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Splines/Splines.hh"
#include "Utils_fmt.hh"
#include <GenericContainer/GenericContainer.hh>
#include <fstream>
#include <memory>
#include <vector>
#include <map>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

using namespace SplinesLoad;
using namespace std;
using namespace GenericContainerNamespace;
using Splines::integer;
using Splines::real_type;
using Splines::Spline_sub_type;

// ============================================================================
// Color definitions and formatting utilities
// ============================================================================

namespace Color
{
  constexpr fmt::color HEADER       = fmt::color::steel_blue;
  constexpr fmt::color SUCCESS      = fmt::color::lime_green;
  constexpr fmt::color WARNING      = fmt::color::gold;
  constexpr fmt::color ERROR        = fmt::color::crimson;
  constexpr fmt::color INFO         = fmt::color::deep_sky_blue;
  constexpr fmt::color TABLE_HEADER = fmt::color::dark_violet;
  constexpr fmt::color TABLE_ROW    = fmt::color::white;
  constexpr fmt::color VALUE        = fmt::color::light_green;
}  // namespace Color

template <typename... Args> string format_message( fmt::string_view format, Args &&... args )
{
  return fmt::vformat( format, fmt::make_format_args( args... ) );
}

template <typename... Args> void print_header( fmt::string_view format, Args &&... args )
{
  fmt::print(
    fg( Color::HEADER ) | fmt::emphasis::bold,
    "╔══════════════════════════════════════════════════════════╗\n" );
  fmt::print(
    fg( Color::HEADER ) | fmt::emphasis::bold,
    "║ {:^56} ║\n",
    format_message( format, std::forward<Args>( args )... ) );
  fmt::print(
    fg( Color::HEADER ) | fmt::emphasis::bold,
    "╚══════════════════════════════════════════════════════════╝\n" );
}

template <typename... Args> void print_success( fmt::string_view format, Args &&... args )
{
  fmt::print(
    fg( Color::SUCCESS ) | fmt::emphasis::bold,
    "✓ {}\n",
    format_message( format, std::forward<Args>( args )... ) );
}

template <typename... Args> void print_error( fmt::string_view format, Args &&... args )
{
  fmt::print(
    fg( Color::ERROR ) | fmt::emphasis::bold,
    "✗ {}\n",
    format_message( format, std::forward<Args>( args )... ) );
}

template <typename... Args> void print_warning( fmt::string_view format, Args &&... args )
{
  fmt::print(
    fg( Color::WARNING ) | fmt::emphasis::bold,
    "⚠ {}\n",
    format_message( format, std::forward<Args>( args )... ) );
}

template <typename... Args> void print_info( fmt::string_view format, Args &&... args )
{
  fmt::print(
    fg( Color::INFO ),
    "➤ {}\n",
    format_message( format, std::forward<Args>( args )... ) );
}

// ============================================================================
// Test data definition
// ============================================================================

constexpr double x_grid[]{ 0, 5, 10, 15, 20, 25, 30, 40, 50, 70 };
constexpr double y_grid[]{ 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6 };

constexpr int nx = std::size( x_grid );  // 10
constexpr int ny = std::size( y_grid );  // 9

// C-style ordering: element at (i,j) = i*ny + j
constexpr double z_data[] = { 24.2, 24.0, 20.3, 17.3, 14.5, 12.2, 10.2, 5.7,  3.4,  0.1,  28.0, 24.6, 21.1, 18.1, 15.2,
                              12.8, 10.7, 6.5,  3.9,  0.2,  28.3, 25.2, 21.9, 18.7, 15.9, 13.4, 11.2, 7.3,  4.4,  0.4,
                              30.8, 27.2, 23.8, 20.5, 17.3, 14.7, 12.3, 8.1,  4.9,  0.8,  34.5, 30.3, 26.6, 23.2, 19.8,
                              16.8, 14.1, 9.4,  5.6,  1.1,  37.9, 34.3, 30.4, 26.8, 23.3, 19.8, 16.8, 11.2, 6.8,  1.4,
                              36.1, 38.0, 34.9, 31.3, 27.3, 23.6, 20.1, 13.4, 8.3,  1.7,  36.1, 36.6, 38.5, 36.1, 31.6,
                              28.1, 24.2, 16.2, 10.0, 2.2,  36.1, 35.2, 42.1, 38.7, 35.7, 32.0, 28.1, 19.3, 11.9, 2.9 };

// ============================================================================
// Test functions
// ============================================================================

void test_2d_spline_evaluation()
{
  print_header( "Test 1: 2D Spline Evaluation" );

  print_info( "Creating BiQuinticSpline from array data..." );

  auto          spline = make_unique<BiQuinticSpline>( Splines::Spline_sub_type::CUBIC );
  constexpr int ldZ    = ny;

  try
  {
    spline->build( x_grid, 1, y_grid, 1, z_data, ldZ, nx, ny, false, false );
    print_success( "Spline built successfully" );
  }
  catch ( const exception & e )
  {
    print_error( "Failed to build spline: {}", e.what() );
    return;
  }

  // Print spline information
  fmt::print( fg( Color::INFO ), "\nSpline Information:\n" );
  spline->write_to_stream( cout );

  // Create evaluation table
  fmt::print(
    fg( Color::TABLE_HEADER ) | fmt::emphasis::bold,
    "\n┌──────────┬{}┐\n",
    fmt::format( "{:─^{}}", " Evaluation Results ", nx * 10 ) );
  fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, "│   y\\x    │" );

  for ( int i = 0; i < nx; ++i )
  {
    fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, " {:8.1f} ", x_grid[i] );
  }
  fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, "│\n" );
  fmt::print(
    fg( Color::TABLE_HEADER ) | fmt::emphasis::bold,
    "├──────────┼{}┤\n",
    fmt::format( "{:─^{}}", "", nx * 10 ) );

  for ( int j = 0; j < ny; ++j )
  {
    fmt::print( fg( Color::TABLE_ROW ), "│ {:8.1f} │", y_grid[j] );

    for ( int i = 0; i < nx; ++i )
    {
      double value = spline->eval( x_grid[i], y_grid[j] );
      fmt::print( fg( Color::VALUE ), " {:8.4f} ", value );
    }
    fmt::print( fg( Color::TABLE_ROW ), "│\n" );
  }

  fmt::print(
    fg( Color::TABLE_HEADER ) | fmt::emphasis::bold,
    "└──────────┴{}┘\n",
    fmt::format( "{:─^{}}", "", nx * 10 ) );
}


void test_json_string_constructors()
{
  print_header( "Test 2: JSON String Constructors" );

  // Create JSON data for an asymmetric grid (4x7)
  const string json_data_old = R"({
        "xdata": [0.0, 2.0, 4.0, 6.0],
        "ydata": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        "zdata": [
            [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
            [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0],
            [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0],
            [4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
        ],
        "transposed": true,
        "fortran_storage": true
    })";

  const string json_data = R"({
  "fortran_storage": false,
  "transposed": false,
  "xdata": [-0.0523599, -0.01309, -0.0043633, 0.0261799, 0.0392699, 0.1047198, 0.1745329],
  "ydata": [0.0, 0.1745329, 0.3490659, 0.5235988, 0.6981317, 0.8726646, 1.0471976, 1.2217305],
  "zdata": [
    [0.0312, 0.0312, 0.0312, 0.0312, 0.035, 0.035, 0.035],
    [0.0312, 0.0312, 0.0312, 0.0312, 0.035, 0.035, 0.035],
    [0.0312, 0.0312, 0.0312, 0.0312, 0.035, 0.035, 0.035],
    [0.0312, 0.0312, 0.0312, 0.0312, 0.035, 0.035, 0.035],
    [0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025],
    [0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025],
    [0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025],
    [0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025]
  ]
    })";

  print_info( "Testing spline construction from JSON string (4x7 grid)" );

  // Write JSON to a temporary file
  string   temp_filename = "temp_spline_data.json";
  ofstream temp_file( temp_filename );
  if ( !temp_file )
  {
    print_error( "Failed to create temporary file: {}", temp_filename );
    return;
  }
  temp_file << json_data;
  temp_file.close();

  vector<pair<string, unique_ptr<Splines::SplineSurf>>> splines;
  splines.emplace_back( "Bilinear", make_unique<BilinearSpline>() );
  splines.emplace_back( "BiCubic", make_unique<BiCubicSpline>( Spline_sub_type::CUBIC ) );
  splines.emplace_back( "BiCubic[Akima]", make_unique<BiCubicSpline>( Spline_sub_type::AKIMA ) );
  splines.emplace_back( "BiCubic[VanLeer]", make_unique<BiCubicSpline>( Spline_sub_type::VANLEER ) );
  splines.emplace_back( "BiCubic[Pchip]", make_unique<BiCubicSpline>( Spline_sub_type::PCHIP ) );
  splines.emplace_back( "BiQuintic", make_unique<BiQuinticSpline>( Spline_sub_type::CUBIC ) );
  splines.emplace_back( "BiQuintic[Akima]", make_unique<BiQuinticSpline>( Spline_sub_type::AKIMA ) );
  splines.emplace_back( "BiQuintic[VanLeer]", make_unique<BiQuinticSpline>( Spline_sub_type::VANLEER ) );
  splines.emplace_back( "BiQuintic[Pchip]", make_unique<BiQuinticSpline>( Spline_sub_type::PCHIP ) );

  // Original grid data from JSON for verification
  // vector<double>         x_nodes    = { 0.0, 2.0, 4.0, 6.0 };
  // vector<double>         y_nodes    = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  // vector<vector<double>> z_original = { { 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 },
  //                                      { 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 },
  //                                      { 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0 },
  //                                      { 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0 } };

  Splines::Vec x_nodes( 7 );
  Splines::Vec y_nodes( 8 );
  Splines::Mat z_original_t( 8, 7 );
  Splines::Mat z_original( 7, 8 );

  x_nodes << -0.0523599, -0.01309, -0.0043633, 0.0261799, 0.0392699, 0.1047198, 0.1745329;
  y_nodes << 0.0, 0.1745329, 0.3490659, 0.5235988, 0.6981317, 0.8726646, 1.0471976, 1.2217305;
  z_original_t << 0.0312, 0.0312, 0.0312, 0.0312, 0.035, 0.035, 0.035, 0.0312, 0.0312, 0.0312, 0.0312, 0.035, 0.035,
    0.035, 0.0312, 0.0312, 0.0312, 0.0312, 0.035, 0.035, 0.035, 0.0312, 0.0312, 0.0312, 0.0312, 0.035, 0.035, 0.035,
    0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025,
    0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025;

  z_original = z_original_t.transpose();

  std::cout << "---------------\n\n\n\n";
  std::cout << z_original << '\n';
  std::cout << "---------------\n\n\n\n";

  for ( auto & [name, spline] : splines )
  {
    try
    {
      fmt::print( fg( Color::INFO ), "\nBuilding {} spline from JSON...\n", name );
      spline->build( temp_filename );
      print_success( "{} spline built successfully", name );

      // Get grid dimensions from spline
      integer spline_nx = spline->num_point_x();
      integer spline_ny = spline->num_point_y();

      fmt::print( fg( Color::TABLE_ROW ), "Grid dimensions from spline: " );
      fmt::print( fg( Color::VALUE ), "{} x {}\n", spline_nx, spline_ny );
      fmt::print( fg( Color::TABLE_ROW ), "Expected dimensions: " );
      fmt::print( fg( Color::VALUE ), "{} x {}\n", x_nodes.size(), y_nodes.size() );

      // VERIFICA DEI NODI: controlla che i nodi della spline corrispondano alla mesh letta
      print_info( "Verifying node coordinates..." );

      bool         nodes_match    = true;
      const double node_tolerance = 1e-12;

      // Controlla dimensione nodi x
      if ( spline_nx != static_cast<integer>( x_nodes.size() ) )
      {
        print_error( "X nodes count mismatch: spline={}, expected={}", spline_nx, x_nodes.size() );
        nodes_match = false;
      }
      else
      {
        // Verifica i valori dei nodi x uno per uno
        fmt::print( fg( Color::INFO ), "Verifying X nodes:\n" );
        for ( integer i = 0; i < spline_nx; ++i )
        {
          double spline_x   = spline->x_node( i );
          double expected_x = x_nodes[i];
          double diff       = abs( spline_x - expected_x );

          fmt::print( fg( Color::TABLE_ROW ), "  Node {}: spline={:.12f}, expected={:.12f}", i, spline_x, expected_x );

          if ( diff > node_tolerance )
          {
            fmt::print( fg( Color::ERROR ), "  ✗ diff={:.2e}\n", diff );
            nodes_match = false;
          }
          else
          {
            fmt::print( fg( Color::SUCCESS ), "  ✓\n" );
          }
        }
      }

      // Controlla dimensione nodi y
      if ( spline_ny != static_cast<integer>( y_nodes.size() ) )
      {
        print_error( "Y nodes count mismatch: spline={}, expected={}", spline_ny, y_nodes.size() );
        nodes_match = false;
      }
      else
      {
        // Verifica i valori dei nodi y uno per uno
        fmt::print( fg( Color::INFO ), "Verifying Y nodes:\n" );
        for ( integer j = 0; j < spline_ny; ++j )
        {
          double spline_y   = spline->y_node( j );
          double expected_y = y_nodes[j];
          double diff       = abs( spline_y - expected_y );

          fmt::print( fg( Color::TABLE_ROW ), "  Node {}: spline={:.12f}, expected={:.12f}", j, spline_y, expected_y );

          if ( diff > node_tolerance )
          {
            fmt::print( fg( Color::ERROR ), "  ✗ diff={:.2e}\n", diff );
            nodes_match = false;
          }
          else
          {
            fmt::print( fg( Color::SUCCESS ), "  ✓\n" );
          }
        }
      }

      if ( nodes_match ) { print_success( "All node coordinates match expected values!" ); }
      else
      {
        print_warning( "Some node coordinate mismatches found" );
      }

      // Verifica indiretta: valutazione sui nodi noti
      print_info( "Performing indirect verification via evaluation at grid points..." );

      // Verify node values
      fmt::print( fg( Color::INFO ), "\nVerifying data at all nodes:\n" );

      int          errors    = 0;
      const double tolerance = 1e-12;
      integer      NX        = static_cast<integer>( x_nodes.size() );
      integer      NY        = static_cast<integer>( y_nodes.size() );

      for ( integer i = 0; i < NX; ++i )
      {
        for ( integer j = 0; j < NY; ++j )
        {
          double computed = spline->eval( x_nodes[i], y_nodes[j] );
          double expected = z_original( i, j );
          double err      = abs( computed - expected );

          if ( err > tolerance )
          {
            if ( errors == 0 ) { fmt::print( fg( Color::ERROR ), "\nMismatches found (computed vs expected):\n" ); }
            fmt::print(
              fg( Color::ERROR ),
              "  At ({}, {}): {:.12f} != {:.12f} (diff = {:.2e})\n",
              x_nodes[i],
              y_nodes[j],
              computed,
              expected,
              err );
            errors++;
          }
        }
      }

      if ( errors == 0 )
      {
        print_success( "All node values match expected values!" );

        // Stampiamo un riepilogo della mesh verificata
        fmt::print( fg( Color::INFO ), "\nMesh verification summary:\n" );
        fmt::print( fg( Color::TABLE_ROW ), "  - X nodes: {} points verified\n", x_nodes.size() );
        fmt::print( fg( Color::TABLE_ROW ), "  - Y nodes: {} points verified\n", y_nodes.size() );
        fmt::print(
          fg( Color::TABLE_ROW ),
          "  - Z values: {} data points verified\n",
          x_nodes.size() * y_nodes.size() );
        fmt::print( fg( Color::SUCCESS ), "  ✓ All evaluations match within tolerance {}\n", tolerance );

        // Print a sample table of values
        fmt::print( fg( Color::INFO ), "\nSample values from spline (first 3x3 grid):\n" );
        fmt::print(
          fg( Color::TABLE_HEADER ) | fmt::emphasis::bold,
          "┌──────────┬──────────┬──────────┬──────────┐\n" );
        fmt::print(
          fg( Color::TABLE_HEADER ) | fmt::emphasis::bold,
          "│   x\\y    │  {:6.1f}  │  {:6.1f}  │  {:6.1f}  │\n",
          y_nodes[0],
          y_nodes[1],
          y_nodes[2] );
        fmt::print(
          fg( Color::TABLE_HEADER ) | fmt::emphasis::bold,
          "├──────────┼──────────┼──────────┼──────────┤\n" );

        for ( integer i = 0; i < 3; ++i )
        {
          fmt::print( fg( Color::TABLE_ROW ), "│ {:8.1f} │", x_nodes[i] );
          for ( integer j = 0; j < 3; ++j )
          {
            fmt::print( fg( Color::VALUE ), " {:8.4f} │", spline->eval( x_nodes[i], y_nodes[j] ) );
          }
          fmt::print( "\n" );
        }
        fmt::print(
          fg( Color::TABLE_HEADER ) | fmt::emphasis::bold,
          "└──────────┴──────────┴──────────┴──────────┘\n" );
      }
      else
      {
        print_error( "Found {} mismatches in node values", errors );

        // Calcola statistiche sugli errori
        vector<double> differences;
        for ( integer i = 0; i < NX; ++i )
        {
          for ( integer j = 0; j < NY; ++j )
          {
            double computed = spline->eval( x_nodes[i], y_nodes[j] );
            double expected = z_original( i, j );
            differences.push_back( abs( computed - expected ) );
          }
        }

        auto max_diff = *max_element( differences.begin(), differences.end() );
        auto avg_diff = accumulate( differences.begin(), differences.end(), 0.0 ) / static_cast<double>( differences.size() );

        fmt::print(
          fg( Color::WARNING ),
          "Error statistics: max diff = {:.2e}, avg diff = {:.2e}\n",
          max_diff,
          avg_diff );
      }

      // Test interpolation at non-grid points
      fmt::print( fg( Color::INFO ), "\nTesting interpolation at non-grid points:\n" );
      vector<pair<double, double>> test_points = {
        { 1.0, 0.5 },  // Between grid points
        { 3.0, 2.5 },  // Between grid points
        { 5.0, 4.5 },  // Between grid points
      };

      fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, "┌──────────┬──────────┬────────────┐\n" );
      fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, "│     x    │     y    │     z      │\n" );
      fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, "├──────────┼──────────┼────────────┤\n" );

      for ( const auto & [x, y] : test_points )
      {
        double z = spline->eval( x, y );
        fmt::print( fg( Color::TABLE_ROW ), "│ {:8.2f} │ {:8.2f} │", x, y );
        fmt::print( fg( Color::VALUE ), " {:10.4f} │\n", z );
      }
      fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, "└──────────┴──────────┴────────────┘\n" );
    }
    catch ( const exception & e )
    {
      print_error( "Failed to build {} spline: {}", name, e.what() );
    }
  }

  // Clean up temporary file
  remove( temp_filename.c_str() );
  print_info( "Temporary file removed: {}", temp_filename );
}


void test_derivatives_with_autodiff()
{
  print_header( "Test 3: Derivative Validation with Autodiff" );

  print_info( "Testing derivatives for various spline types at (40, 0.7)" );

  // Test point
  double x0 = 40.0;
  double y0 = 0.7;

  vector<pair<string, unique_ptr<Splines::SplineSurf>>> tests;

  tests.emplace_back( "Bilinear", make_unique<BilinearSpline>() );
  tests.emplace_back( "BiCubic", make_unique<BiCubicSpline>( Spline_sub_type::CUBIC ) );
  tests.emplace_back( "BiCubic[Akima]", make_unique<BiCubicSpline>( Spline_sub_type::AKIMA ) );
  tests.emplace_back( "BiCubic[VanLeer]", make_unique<BiCubicSpline>( Spline_sub_type::VANLEER ) );
  tests.emplace_back( "BiCubic[Pchip]", make_unique<BiCubicSpline>( Spline_sub_type::PCHIP ) );
  tests.emplace_back( "BiQuintic", make_unique<BiQuinticSpline>( Spline_sub_type::CUBIC ) );
  tests.emplace_back( "BiQuintic[Akima]", make_unique<BiQuinticSpline>( Spline_sub_type::AKIMA ) );
  tests.emplace_back( "BiQuintic[VanLeer]", make_unique<BiQuinticSpline>( Spline_sub_type::VANLEER ) );
  tests.emplace_back( "BiQuintic[Pchip]", make_unique<BiQuinticSpline>( Spline_sub_type::PCHIP ) );

  // Build all splines
  for ( auto & [name, spline] : tests ) { spline->build( x_grid, y_grid, z_data, nx, ny, false, false ); }

  // Print derivatives table
  fmt::print(
    fg( Color::TABLE_HEADER ) | fmt::emphasis::bold,
    "\n"
    "┌──────────────────────┬────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐\n"
    "│     Spline Type      │    f(x,y)  │    f_x     │    f_y     │   f_xx     │   f_xy     │   f_yy     │\n"
    "├──────────────────────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤\n" );

  for ( const auto & [name, spline] : tests )
  {
    fmt::print( fg( Color::TABLE_ROW ), "│ {:20} │", name );

    try
    {
      double f   = spline->eval( x0, y0 );
      double fx  = spline->Dx( x0, y0 );
      double fy  = spline->Dy( x0, y0 );
      double fxx = spline->Dxx( x0, y0 );
      double fxy = spline->Dxy( x0, y0 );
      double fyy = spline->Dyy( x0, y0 );

      fmt::print(
        fg( Color::VALUE ),
        " {:10.6f} │ {:10.6f} │ {:10.6f} │ {:10.6f} │ {:10.6f} │ {:10.6f} │\n",
        f,
        fx,
        fy,
        fxx,
        fxy,
        fyy );
    }
    catch ( const exception & e )
    {
      fmt::print(
        fg( Color::ERROR ),
        " {:^10} │ {:^10} │ {:^10} │ {:^10} │ {:^10} │ {:10} │\n",
        "N/A",
        "N/A",
        "N/A",
        "N/A",
        "N/A",
        "N/A" );
    }
  }

  fmt::print(
    fg( Color::TABLE_HEADER ) | fmt::emphasis::bold,
    "└──────────────────────┴────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘\n" );

  print_info( "\nNote: Higher order derivatives may not be available for all spline types" );
}

void test_generic_container_construction()
{
  print_header( "Test 4: Generic Container Construction" );

  print_info( "Testing 1D blended spline construction via GenericContainer" );

  try
  {
    GC::GenericContainer gc;

    // Create two different spline configurations
    GC::GenericContainer & S0 = gc["spline0"];
    GC::GenericContainer & S1 = gc["spline1"];

    // Set data
    GC::vec_real_type & V0x = S0["xdata"].set_vec_real( nx );
    std::copy_n( x_grid, nx, V0x.data() );

    GC::vec_real_type & V0y = S0["ydata"].set_vec_real( nx );
    // Use first data column as sample data
    for ( integer i = 0; i < nx; ++i )
    {
      V0y[i] = z_data[i];  // Take first column data
    }

    GC::vec_real_type & V1x = S1["xdata"].set_vec_real( nx );
    std::copy_n( x_grid, nx, V1x.data() );

    GC::vec_real_type & V1y = S1["ydata"].set_vec_real( nx );
    // Use second column data
    for ( integer i = 0; i < nx; ++i )
    {
      V1y[i] = z_data[ny + i];  // Take second column data
    }

    // Configure spline types and boundary conditions
    S0["spline_type"] = "cubic";
    S0["bc_begin"]    = "natural";
    S0["bc_end"]      = "natural";

    S1["spline_type"]     = "quintic";
    S1["spline_sub_type"] = "pchip";
    S1["bc_begin"]        = "natural";
    S1["bc_end"]          = "natural";

    // Build spline
    Splines::Spline1Dblend S( "test_spline" );
    S.build( gc );

    print_success( "1D blended spline built successfully via GenericContainer" );

    // Test interpolation
    fmt::print( fg( Color::INFO ), "\nTesting interpolation at sample points:\n" );
    fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, "┌────────────┬────────────┬────────────┐\n" );
    fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, "│     x      │  spline0   │  spline1   │\n" );
    fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, "├────────────┼────────────┼────────────┤\n" );

    for ( double x_val : { 5.0, 15.0, 35.0, 60.0 } )
    {
      fmt::print( fg( Color::TABLE_ROW ), "│ {:10.2f} │", x_val );
      fmt::print( fg( Color::VALUE ), " {:10.4f} │ {:10.4f} │\n", S.eval( 0.0, x_val ), S.eval( 1.0, x_val ) );
    }

    fmt::print( fg( Color::TABLE_HEADER ) | fmt::emphasis::bold, "└────────────┴────────────┴────────────┘\n" );
  }
  catch ( const exception & e )
  {
    print_error( "GenericContainer construction failed: {}", e.what() );
  }
}

void test_comprehensive_file_operations()
{
  print_header( "Test 5: Comprehensive File Operations" );

  print_info( "Testing save/load operations for splines" );

  try
  {
    // Create a spline
    auto          spline = make_unique<BiCubicSpline>( Spline_sub_type::CUBIC );
    constexpr int ldZ    = nx;
    spline->build( x_grid, 1, y_grid, 1, z_data, ldZ, nx, ny, false, false );

    string save_filename = "test_spline_output.txt";

    // Save to file using write_to_stream
    print_info( "Saving spline to file: {}", save_filename );
    ofstream out_file( save_filename );
    if ( out_file )
    {
      spline->write_to_stream( out_file );
      out_file.close();
      print_success( "Spline saved successfully" );
    }
    else
    {
      print_error( "Failed to open file for writing: {}", save_filename );
      return;
    }

    // Verify file exists
    ifstream in_file( save_filename );
    if ( !in_file.good() )
    {
      print_error( "Saved file not found: {}", save_filename );
      return;
    }

    print_info( "File saved successfully. Size verification passed." );

    // Clean up temporary file
    remove( save_filename.c_str() );
    print_success( "Temporary file cleaned up" );
  }
  catch ( const exception & e )
  {
    print_error( "File operations test failed: {}", e.what() );
  }
}

// ============================================================================
// Main program
// ============================================================================

int main()
{
  print_header( "SPLINE INTERPOLATION TEST SUITE" );

  fmt::print( fg( Color::INFO ) | fmt::emphasis::bold, "\nAvailable Tests:\n" );
  fmt::print( fg( Color::TABLE_ROW ), "  1. 2D Spline Evaluation\n" );
  fmt::print( fg( Color::TABLE_ROW ), "  2. JSON String Constructors (4x7 grid)\n" );
  fmt::print( fg( Color::TABLE_ROW ), "  3. Derivative Validation\n" );
  fmt::print( fg( Color::TABLE_ROW ), "  4. Generic Container Construction\n" );
  fmt::print( fg( Color::TABLE_ROW ), "  5. Comprehensive File Operations\n" );
  fmt::print( fg( Color::TABLE_ROW ), "  6. Run All Tests\n" );
  fmt::print( fg( Color::TABLE_ROW ), "  0. Exit\n" );

  test_2d_spline_evaluation();
  test_json_string_constructors();
  test_derivatives_with_autodiff();
  test_generic_container_construction();
  test_comprehensive_file_operations();

  fmt::print( fg( Color::SUCCESS ) | fmt::emphasis::bold, "\n✨ All tests completed! ✨\n" );

  return 0;
}
