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
#include "Utils_fmt.hh"
#include "Utils_eigen.hh"
#include "Utils_autodiff.hh"

using Utils::m_pi;

#include <random>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

namespace SplinesTest
{

  using namespace Splines;
  using namespace std;
  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  // Alias per tipi autodiff
  using autodiff::dual1st;
  using autodiff::dual2nd;

  // ===========================================================================
  // TEST FUNCTIONS - Well-behaved smooth functions for spline testing
  // ===========================================================================

  // Function 1: Simple 2D polynomial (for basic testing)
  real_type poly_function( real_type x, real_type y )
  { return 1.0 + 0.5 * x - 0.3 * y + 0.2 * x * x - 0.1 * x * y + 0.05 * y * y; }

  // Function 2: Smooth exponential function (for more challenging tests)
  real_type exp_function( real_type x, real_type y )
  { return exp( -0.1 * ( x * x + y * y ) ) * cos( x ) * sin( y ); }

  // Function 3: Gentle sine-cosine product (for accuracy tests)
  real_type trig_function( real_type x, real_type y )
  { return sin( 0.5 * x ) * cos( 0.5 * y ); }

  // ===========================================================================
  // UTILITY FUNCTIONS - AUTODIFF
  // ===========================================================================

  // Funzioni per calcolare derivate usando autodiff
  struct AutoDiffDerivatives
  {
    real_type Dx;   // derivata prima rispetto a x
    real_type Dy;   // derivata prima rispetto a y
    real_type Dxx;  // derivata seconda rispetto a x
    real_type Dyy;  // derivata seconda rispetto a y
    real_type Dxy;  // derivata mista
  };

  // Calcola tutte le derivate usando autodiff
  template <typename Func> AutoDiffDerivatives compute_autodiff_derivatives( Func const & f, real_type x, real_type y )
  {
    AutoDiffDerivatives result;

    // Derivata prima rispetto a x
    {
      autodiff::dual1st xd;
      xd.val    = x;
      xd.grad   = 1;
      auto Dx   = f.eval( xd, y );
      result.Dx = Dx.grad;
    }

    // Derivata prima rispetto a y
    {
      autodiff::dual1st yd;
      yd.val    = y;
      yd.grad   = 1;
      auto Dy   = f.eval( x, yd );
      result.Dy = Dy.grad;
    }

    // Derivata seconda rispetto a x
    {
      autodiff::dual2nd xdd;
      xdd.val.val   = x;
      xdd.val.grad  = 1;
      xdd.grad.val  = 1;
      xdd.grad.grad = 0;
      auto Dxx      = f.eval( xdd, y );
      result.Dxx    = Dxx.grad.grad;
    }

    // Derivata seconda rispetto a y
    {
      autodiff::dual2nd ydd;
      ydd.val.val   = y;
      ydd.val.grad  = 1;
      ydd.grad.val  = 1;
      ydd.grad.grad = 0;
      auto Dyy      = f.eval( x, ydd );
      result.Dyy    = Dyy.grad.grad;
    }

    // Derivata mista: prima rispetto a x, poi rispetto a y
    // fxy = d/dy (df/dx)
    {
      autodiff::dual2nd xdd, ydd;

      xdd.val.val   = x;
      xdd.val.grad  = 1;
      xdd.grad.val  = 0;
      xdd.grad.grad = 0;

      ydd.val.val   = y;
      ydd.val.grad  = 0;
      ydd.grad.val  = 1;
      ydd.grad.grad = 0;

      auto Dxy   = f.eval( xdd, ydd );
      result.Dxy = Dxy.grad.grad;
    }

    if ( abs( result.Dxy - f.Dxy( x, y ) ) > 1E-10 )
    {
      autodiff::dual2nd xdd, ydd;

      xdd.val.val   = x;
      xdd.val.grad  = 1;
      xdd.grad.val  = 0;
      xdd.grad.grad = 0;

      ydd.val.val   = y;
      ydd.val.grad  = 0;
      ydd.grad.val  = 1;
      ydd.grad.grad = 0;

      auto Dxy1 = f.eval( xdd, ydd );
      auto Dxy2 = f.Dxy( x, y );
      fmt::print( "f_xy({:.3},{:.3}) = {:.3}, {:.3}, diff={:.3}\n", x, y, Dxy1.grad.grad, Dxy2, Dxy1.grad.grad - Dxy2 );
    }

    return result;
  }

  // Statistics structure con punto di errore massimo
  struct DerivativeStats
  {
    string    name;
    real_type max_abs_error = 0;
    real_type rms_error     = 0;
    real_type max_rel_error = 0;
    integer   count         = 0;
    bool      is_derivative = false;

    // Punti dove si verificano gli errori massimi
    real_type max_error_x   = 0;
    real_type max_error_y   = 0;
    real_type exact_at_max  = 0;
    real_type approx_at_max = 0;

    // Nuovi campi per memorizzare i valori a confronto
    real_type spline_value_at_max = 0;
    real_type ad_value_at_max     = 0;

    void update( real_type x, real_type y, real_type exact, real_type approx )
    {
      real_type abs_err = abs( exact - approx );

      // Aggiorna punto di errore massimo
      if ( abs_err > max_abs_error )
      {
        max_abs_error = abs_err;
        max_error_x   = x;
        max_error_y   = y;
        exact_at_max  = exact;
        approx_at_max = approx;

        // Memorizza i valori a confronto
        spline_value_at_max = exact;   // Valore dalla spline
        ad_value_at_max     = approx;  // Valore da autodiff
      }

      real_type rel_err = 0.0;
      if ( abs( exact ) > 1e-12 ) { rel_err = abs_err / abs( exact ); }
      max_rel_error = max( max_rel_error, rel_err );

      rms_error += abs_err * abs_err;
      ++count;
    }

    void finalize()
    {
      if ( count > 0 ) { rms_error = sqrt( rms_error / static_cast<real_type>( count ) ); }
    }

    void print_row() const
    {
      fmt::color color = max_abs_error < 1e-6   ? fmt::color::green
                         : max_abs_error < 1e-3 ? fmt::color::yellow
                                                : fmt::color::red;

      if ( is_derivative )
      {
        fmt::print(
          fg( color ),
          "│ {:<14} │ {:^12.2e} │ {:^12.2e} │ {:^12.2e} │ {:18} │\n",
          name,
          max_abs_error,
          rms_error,
          max_rel_error,
          fmt::format( "({:.3f},{:.3f})", max_error_x, max_error_y ) );
      }
      else
      {
        fmt::print(
          fg( color ),
          "│ {:<14} │ {:^12.2e} │ {:^12.2e} │ {:^12.2e} │ {:18} │\n",
          name,
          max_abs_error,
          rms_error,
          max_rel_error,
          fmt::format( "({:.3f},{:.3f})", max_error_x, max_error_y ) );
      }
    }

    void print_max_error_details() const
    {
      if ( count > 0 )
      {
        if ( is_derivative )
        {
          fmt::print(
            fg( fmt::color::light_gray ),
            "    Max error at ({:.4f}, {:.4f}):\n"
            "      Spline: {:.6e}\n"
            "      Autodiff: {:.6e}\n"
            "      Difference: {:.2e} (rel: {:.2e}%)\n",
            max_error_x,
            max_error_y,
            spline_value_at_max,
            ad_value_at_max,
            max_abs_error,
            max_rel_error * 100 );
        }
        else
        {
          fmt::print(
            fg( fmt::color::light_gray ),
            "    Max error at ({:.4f}, {:.4f}): exact={:.6e}, approx={:.6e}, diff={:.2e}\n",
            max_error_x,
            max_error_y,
            exact_at_max,
            approx_at_max,
            max_abs_error );
        }
      }
    }
  };

  struct SplineTestResults
  {
    string          spline_name;
    DerivativeStats value{ "f(x,y)", 0, 0, 0, 0, false };
    DerivativeStats dx_ad{ "∂f/∂x (AD)", 0, 0, 0, 0, true };  // Consistency with Autodiff
    DerivativeStats dy_ad{ "∂f/∂y (AD)", 0, 0, 0, 0, true };
    DerivativeStats dxx_ad{ "∂²f/∂x² (AD)", 0, 0, 0, 0, true };
    DerivativeStats dyy_ad{ "∂²f/∂y² (AD)", 0, 0, 0, 0, true };
    DerivativeStats dxy_ad{ "∂²f/∂x∂y (AD)", 0, 0, 0, 0, true };

    // Statistiche per interpolazione sui punti di griglia
    DerivativeStats grid_interp{ "Grid Interp", 0, 0, 0, 0, false };
    real_type       grid_max_interp_error{ 0 };
    real_type       grid_max_interp_x{ 0 };
    real_type       grid_max_interp_y{ 0 };

    SplineTestResults( const string & name ) : spline_name( name ) {}

    void finalize()
    {
      value.finalize();
      dx_ad.finalize();
      dy_ad.finalize();
      dxx_ad.finalize();
      dyy_ad.finalize();
      dxy_ad.finalize();
      grid_interp.finalize();
    }

    void print_table() const
    {
      fmt::print(
        fg( fmt::color::cyan ),
        "\n┌{0:─^{2}}┐\n"
        "│{1: ^{2}}│\n"
        "└{0:─^{2}}┘\n",
        "",
        spline_name,
        78 );

      fmt::print(
        "┌────────────────┬──────────────┬──────────────┬──────────────┬────────────────────┐\n"
        "│ Derivative     │ Max Abs Err  │ RMS Error    │ Max Rel Err  │    Max Err At      │\n"
        "├────────────────┼──────────────┼──────────────┼──────────────┼────────────────────┤\n" );

      value.print_row();

      fmt::print(
        "├────────────────┼──────────────┼──────────────┼──────────────┼────────────────────┤\n"
        "│ AD Consistency │              │              │              │                    │\n" );

      dx_ad.print_row();
      dy_ad.print_row();

      if ( dxx_ad.count > 0 )
      {
        dxx_ad.print_row();
        dyy_ad.print_row();
        dxy_ad.print_row();
      }

      // Aggiungi riga per interpolazione sui punti di griglia
      fmt::print( "├────────────────┼──────────────┼──────────────┼──────────────┼────────────────────┤\n" );
      grid_interp.print_row();

      fmt::print( "└────────────────┴──────────────┴──────────────┴──────────────┴────────────────────┘\n" );

      // Print summary con dettagli errori massimi
      fmt::print( "\n📋 Detailed comparison for {}:\n", spline_name );

      // Dettagli errori massimi per ogni tipo
      fmt::print( "  Max error details:\n" );
      value.print_max_error_details();

      // Dettagli per derivata AD
      if ( dx_ad.max_abs_error > 1e-10 )
      {
        fmt::print( "  ∂f/∂x comparison:\n" );
        dx_ad.print_max_error_details();
      }
      if ( dy_ad.max_abs_error > 1e-10 )
      {
        fmt::print( "  ∂f/∂y comparison:\n" );
        dy_ad.print_max_error_details();
      }
      if ( dxx_ad.max_abs_error > 1e-10 )
      {
        fmt::print( "  ∂²f/∂x² comparison:\n" );
        dxx_ad.print_max_error_details();
      }
      if ( dyy_ad.max_abs_error > 1e-10 )
      {
        fmt::print( "  ∂²f/∂y² comparison:\n" );
        dyy_ad.print_max_error_details();
      }
      if ( dxy_ad.max_abs_error > 1e-10 )
      {
        fmt::print( "  ∂²f/∂x∂y comparison:\n" );
        dxy_ad.print_max_error_details();
      }

      // Dettagli per interpolazione griglia
      if ( grid_interp.count > 0 )
      {
        fmt::print( "  Grid interpolation:\n" );
        grid_interp.print_max_error_details();
      }

      bool all_good = true;
      if ( value.max_abs_error > 1e-4 )
      {
        fmt::print( fg( fmt::color::yellow ), "  ⚠️  Function error: {:.2e}\n", value.max_abs_error );
        all_good = false;
      }

      if ( grid_interp.max_abs_error > 1e-10 )
      {
        fmt::print( fg( fmt::color::yellow ), "  ⚠️  Grid interpolation error: {:.2e}\n", grid_interp.max_abs_error );
        all_good = false;
      }

      if ( dx_ad.max_abs_error > 1e-6 )
      {
        fmt::print( fg( fmt::color::red ), "  ❌ ∂f/∂x AD consistency: {:.2e}\n", dx_ad.max_abs_error );
        all_good = false;
      }

      if ( dy_ad.max_abs_error > 1e-6 )
      {
        fmt::print( fg( fmt::color::red ), "  ❌ ∂f/∂y AD consistency: {:.2e}\n", dy_ad.max_abs_error );
        all_good = false;
      }

      if ( dxx_ad.count > 0 )
      {
        if ( dxx_ad.max_abs_error > 1e-4 )
        {
          fmt::print( fg( fmt::color::yellow ), "  ⚠️  ∂²f/∂x² AD consistency: {:.2e}\n", dxx_ad.max_abs_error );
          all_good = false;
        }

        if ( dyy_ad.max_abs_error > 1e-4 )
        {
          fmt::print( fg( fmt::color::yellow ), "  ⚠️  ∂²f/∂y² AD consistency: {:.2e}\n", dyy_ad.max_abs_error );
          all_good = false;
        }

        if ( dxy_ad.max_abs_error > 1e-4 )
        {
          fmt::print( fg( fmt::color::yellow ), "  ⚠️  ∂²f/∂x∂y AD consistency: {:.2e}\n", dxy_ad.max_abs_error );
          all_good = false;
        }
      }

      if ( all_good ) { fmt::print( fg( fmt::color::green ), "  ✅ All tests passed!\n" ); }
    }
  };

  // ===========================================================================
  // TEST FUNCTION - MODIFICATA PER USARE AUTODIFF
  // ===========================================================================

  void test_spline_type(
    SplineType2D     spline_type,
    const VectorXd & x_grid,
    const VectorXd & y_grid,
    const MatrixXd & z_data,
    integer          n_test_points = 500 )
  {
    string type_name;
    switch ( spline_type )
    {
      case SplineType2D::BILINEAR: type_name = "Bilinear"; break;
      case SplineType2D::BICUBIC_CUBIC: type_name = "Bicubic"; break;
      case SplineType2D::BICUBIC_AKIMA: type_name = "Bicubic[Akima]"; break;
      case SplineType2D::BICUBIC_VANLEER: type_name = "Bicubic[VanLeer]"; break;
      case SplineType2D::BICUBIC_PCHIP: type_name = "Bicubic[Pchip]"; break;
      case SplineType2D::BIQUINTIC_CUBIC: type_name = "Biquintic"; break;
      case SplineType2D::BIQUINTIC_AKIMA: type_name = "Biquintic[Akima]"; break;
      case SplineType2D::BIQUINTIC_VANLEER: type_name = "Biquintic[VanLeer]"; break;
      case SplineType2D::BIQUINTIC_PCHIP: type_name = "Biquintic[Pchip]"; break;
      default: type_name = "Unknown"; break;
    }

    // Create and build spline
    Spline2D spline( type_name );
    integer  nx = static_cast<integer>( x_grid.size() );
    integer  ny = static_cast<integer>( y_grid.size() );
    spline.build( spline_type, x_grid.data(), 1, y_grid.data(), 1, z_data.data(), nx, nx, ny, true, false );

    // Get bounds
    real_type x_min = spline.x_min();
    real_type x_max = spline.x_max();
    real_type y_min = spline.y_min();
    real_type y_max = spline.y_max();

    // Test interpolation at grid points
    SplineTestResults results( type_name );

    // Prima testiamo l'interpolazione sui punti di griglia
    fmt::print( fg( fmt::color::light_blue ), "\n  Testing interpolation at grid points for {}...\n", type_name );

    for ( integer i = 0; i < nx; ++i )
    {
      for ( integer j = 0; j < ny; ++j )
      {
        real_type x          = x_grid[i];
        real_type y          = y_grid[j];
        real_type exact_val  = z_data( i, j );
        real_type spline_val = spline.eval( x, y );
        results.grid_interp.update( x, y, exact_val, spline_val );
      }
    }

    fmt::print( fg( fmt::color::light_green ), "    Grid points tested: {} × {} = {}\n", nx, ny, nx * ny );

    // Generate test points using Eigen (avoid boundaries)
    std::mt19937                              rng( 1000 + static_cast<unsigned>( spline_type ) );
    std::uniform_real_distribution<real_type> dist_x( x_min * 0.9, x_max * 0.9 );
    std::uniform_real_distribution<real_type> dist_y( y_min * 0.9, y_max * 0.9 );

    // Use trigonometric function for testing (smooth and well-behaved)
    auto & func = trig_function;

    // Test points - use Eigen for generating test points
    VectorXd x_test = VectorXd::NullaryExpr( n_test_points, [&]() { return dist_x( rng ); } );
    VectorXd y_test = VectorXd::NullaryExpr( n_test_points, [&]() { return dist_y( rng ); } );

    for ( integer i_test = 0; i_test < n_test_points; ++i_test )
    {
      real_type x = x_test( i_test );
      real_type y = y_test( i_test );

      // Avoid points too close to boundaries
      real_type margin = 0.05 * ( x_max - x_min );
      if ( x < x_min + margin || x > x_max - margin || y < y_min + margin || y > y_max - margin ) continue;

      // Exact values
      real_type exact_val = func( x, y );

      // Spline values
      real_type spline_val = spline.eval( x, y );
      real_type spline_dx  = spline.Dx( x, y );
      real_type spline_dy  = spline.Dy( x, y );
      real_type spline_dxx = spline.Dxx( x, y );
      real_type spline_dyy = spline.Dyy( x, y );
      real_type spline_dxy = spline.Dxy( x, y );

      // Update accuracy stats for function value
      results.value.update( x, y, exact_val, spline_val );

      // Autodiff consistency check
      // Calcoliamo le derivate con autodiff
      auto ad_derivs = compute_autodiff_derivatives( spline, x, y );

      // Confronto derivate spline con autodiff
      results.dx_ad.update( x, y, spline_dx, ad_derivs.Dx );
      results.dy_ad.update( x, y, spline_dy, ad_derivs.Dy );

      if ( spline_type != SplineType2D::BILINEAR )
      {
        results.dxx_ad.update( x, y, spline_dxx, ad_derivs.Dxx );
        results.dyy_ad.update( x, y, spline_dyy, ad_derivs.Dyy );
        results.dxy_ad.update( x, y, spline_dxy, ad_derivs.Dxy );
      }
    }

    results.finalize();
    results.print_table();
  }

  // Test interpolation at grid points
  void test_interpolation_accuracy( const VectorXd & x_grid, const VectorXd & y_grid, const MatrixXd & z_data )
  {
    fmt::print(
      fg( fmt::color::light_blue ) | fmt::emphasis::bold,
      "\n"
      "╔══════════════════════════════════════════════════════════╗\n"
      "║          Interpolation Accuracy at Grid Points           ║\n"
      "╚══════════════════════════════════════════════════════════╝\n" );

    // Test with Bicubic spline (should interpolate exactly for all types)
    Spline2D spline( "TestSpline" );
    integer  nx = static_cast<integer>( x_grid.size() );
    integer  ny = static_cast<integer>( y_grid.size() );
    spline.build(
      SplineType2D::BIQUINTIC_CUBIC,
      x_grid.data(),
      1,
      y_grid.data(),
      1,
      z_data.data(),
      nx,
      nx,
      ny,
      true,
      false );

    real_type max_error         = 0.0;
    real_type max_error_x       = 0.0;
    real_type max_error_y       = 0.0;
    real_type spline_val_at_max = 0.0;
    real_type exact_val_at_max  = 0.0;
    real_type avg_error         = 0.0;
    integer   count             = 0;

    for ( integer i = 0; i < nx; ++i )
    {
      for ( integer j = 0; j < ny; ++j )
      {
        real_type exact  = z_data( i, j );
        real_type interp = spline.eval( x_grid[i], y_grid[j] );
        real_type error  = abs( interp - exact );

        if ( error > max_error )
        {
          max_error         = error;
          max_error_x       = x_grid[i];
          max_error_y       = y_grid[j];
          spline_val_at_max = interp;
          exact_val_at_max  = exact;
        }

        avg_error += error;
        ++count;
      }
    }

    avg_error /= static_cast<real_type>( count );

    fmt::print( "\nGrid interpolation test ({}×{} points):\n", nx, ny );
    fmt::print( "  Maximum error: {:12.2e} at ({:.4f}, {:.4f})\n", max_error, max_error_x, max_error_y );
    fmt::print( "  Spline value at max error: {:.6e}\n", spline_val_at_max );
    fmt::print( "  Exact value at max error: {:.6e}\n", exact_val_at_max );
    fmt::print( "  Difference: {:.2e}\n", max_error );
    fmt::print( "  Average error: {:12.2e}\n", avg_error );

    if ( max_error < 1e-10 ) { fmt::print( fg( fmt::color::green ), "  ✅ Interpolation is exact at grid points\n" ); }
    else if ( max_error < 1e-6 )
    {
      fmt::print( fg( fmt::color::yellow ), "  ⚠️  Small interpolation errors (numerical precision)\n" );
    }
    else
    {
      fmt::print( fg( fmt::color::red ), "  ❌ Large interpolation errors - check implementation!\n" );
    }
  }

  // ===========================================================================
  // NUOVA FUNZIONE: TEST CONSISTENZA DERIVATE CON AUTODIFF
  // ===========================================================================

  void test_autodiff_consistency( Spline2D & spline, const VectorXd & x_grid, const VectorXd & y_grid )
  {
    fmt::print(
      fg( fmt::color::light_blue ) | fmt::emphasis::bold,
      "\n"
      "╔══════════════════════════════════════════════════════════╗\n"
      "║          Autodiff Derivative Consistency Test            ║\n"
      "╚══════════════════════════════════════════════════════════╝\n" );

    std::mt19937                              rng( 1234 );
    std::uniform_real_distribution<real_type> dist_x( x_grid( 0 ) * 0.8, x_grid.tail( 1 )( 0 ) * 0.8 );
    std::uniform_real_distribution<real_type> dist_y( y_grid( 0 ) * 0.8, y_grid.tail( 1 )( 0 ) * 0.8 );

    // Statistiche
    struct ADConsistencyStats
    {
      real_type max_Dx_error = 0, max_Dy_error = 0;
      real_type max_Dxx_error = 0, max_Dyy_error = 0, max_Dxy_error = 0;
      real_type avg_Dx_error = 0, avg_Dy_error = 0;
      real_type avg_Dxx_error = 0, avg_Dyy_error = 0, avg_Dxy_error = 0;
      integer   count = 0;

      real_type max_Dx_x = 0, max_Dx_y = 0;
      real_type max_Dy_x = 0, max_Dy_y = 0;
      real_type max_Dxx_x = 0, max_Dxx_y = 0;
      real_type max_Dyy_x = 0, max_Dyy_y = 0;
      real_type max_Dxy_x = 0, max_Dxy_y = 0;

      // Valori a confronto nel punto di errore massimo
      real_type spline_Dx_at_max = 0, ad_Dx_at_max = 0;
      real_type spline_Dy_at_max = 0, ad_Dy_at_max = 0;
      real_type spline_Dxx_at_max = 0, ad_Dxx_at_max = 0;
      real_type spline_Dyy_at_max = 0, ad_Dyy_at_max = 0;
      real_type spline_Dxy_at_max = 0, ad_Dxy_at_max = 0;
    } stats;

    // Generiamo punti di test
    VectorXd x_test = VectorXd::NullaryExpr( 100, [&]() { return dist_x( rng ); } );
    VectorXd y_test = VectorXd::NullaryExpr( 100, [&]() { return dist_y( rng ); } );

    for ( integer i = 0; i < 100; ++i )
    {
      real_type x = x_test( i );
      real_type y = y_test( i );

      // Skip boundary regions
      real_type margin = 0.1 * ( x_grid.tail( 1 )( 0 ) - x_grid( 0 ) );
      if (
        x < x_grid( 0 ) + margin || x > x_grid.tail( 1 )( 0 ) - margin || y < y_grid( 0 ) + margin ||
        y > y_grid.tail( 1 )( 0 ) - margin )
        continue;

      // Valori dalla spline
      real_type spline_Dx  = spline.Dx( x, y );
      real_type spline_Dy  = spline.Dy( x, y );
      real_type spline_Dxx = spline.Dxx( x, y );
      real_type spline_Dyy = spline.Dyy( x, y );
      real_type spline_Dxy = spline.Dxy( x, y );

      // Valori da autodiff
      auto ad_derivs = compute_autodiff_derivatives( spline, x, y );

      // Calcola differenze
      real_type Dx_diff  = abs( spline_Dx - ad_derivs.Dx );
      real_type Dy_diff  = abs( spline_Dy - ad_derivs.Dy );
      real_type Dxx_diff = abs( spline_Dxx - ad_derivs.Dxx );
      real_type Dyy_diff = abs( spline_Dyy - ad_derivs.Dyy );
      real_type Dxy_diff = abs( spline_Dxy - ad_derivs.Dxy );

      // Aggiorna statistiche
      stats.avg_Dx_error += Dx_diff;
      stats.avg_Dy_error += Dy_diff;
      stats.avg_Dxx_error += Dxx_diff;
      stats.avg_Dyy_error += Dyy_diff;
      stats.avg_Dxy_error += Dxy_diff;

      if ( Dx_diff > stats.max_Dx_error )
      {
        stats.max_Dx_error     = Dx_diff;
        stats.max_Dx_x         = x;
        stats.max_Dx_y         = y;
        stats.spline_Dx_at_max = spline_Dx;
        stats.ad_Dx_at_max     = ad_derivs.Dx;
      }
      if ( Dy_diff > stats.max_Dy_error )
      {
        stats.max_Dy_error     = Dy_diff;
        stats.max_Dy_x         = x;
        stats.max_Dy_y         = y;
        stats.spline_Dy_at_max = spline_Dy;
        stats.ad_Dy_at_max     = ad_derivs.Dy;
      }
      if ( Dxx_diff > stats.max_Dxx_error )
      {
        stats.max_Dxx_error     = Dxx_diff;
        stats.max_Dxx_x         = x;
        stats.max_Dxx_y         = y;
        stats.spline_Dxx_at_max = spline_Dxx;
        stats.ad_Dxx_at_max     = ad_derivs.Dxx;
      }
      if ( Dyy_diff > stats.max_Dyy_error )
      {
        stats.max_Dyy_error     = Dyy_diff;
        stats.max_Dyy_x         = x;
        stats.max_Dyy_y         = y;
        stats.spline_Dyy_at_max = spline_Dyy;
        stats.ad_Dyy_at_max     = ad_derivs.Dyy;
      }
      if ( Dxy_diff > stats.max_Dxy_error )
      {
        stats.max_Dxy_error     = Dxy_diff;
        stats.max_Dxy_x         = x;
        stats.max_Dxy_y         = y;
        stats.spline_Dxy_at_max = spline_Dxy;
        stats.ad_Dxy_at_max     = ad_derivs.Dxy;
      }

      ++stats.count;
    }

    // Calcola medie
    if ( stats.count > 0 )
    {
      stats.avg_Dx_error /= stats.count;
      stats.avg_Dy_error /= stats.count;
      stats.avg_Dxx_error /= stats.count;
      stats.avg_Dyy_error /= stats.count;
      stats.avg_Dxy_error /= stats.count;
    }

    fmt::print( "\nAutodiff consistency test ({} points):\n", stats.count );

    // Stampa con confronto dettagliato dei valori
    fmt::print( "  ∂f/∂x comparison:\n" );
    fmt::print( "    Max diff: {:12.2e} at ({:.4f}, {:.4f})\n", stats.max_Dx_error, stats.max_Dx_x, stats.max_Dx_y );
    fmt::print( "    Spline value:  {:.6e}\n", stats.spline_Dx_at_max );
    fmt::print( "    Autodiff value: {:.6e}\n", stats.ad_Dx_at_max );
    fmt::print( "    Difference:    {:.2e} (avg: {:.2e})\n", stats.max_Dx_error, stats.avg_Dx_error );

    fmt::print( "\n  ∂f/∂y comparison:\n" );
    fmt::print( "    Max diff: {:12.2e} at ({:.4f}, {:.4f})\n", stats.max_Dy_error, stats.max_Dy_x, stats.max_Dy_y );
    fmt::print( "    Spline value:  {:.6e}\n", stats.spline_Dy_at_max );
    fmt::print( "    Autodiff value: {:.6e}\n", stats.ad_Dy_at_max );
    fmt::print( "    Difference:    {:.2e} (avg: {:.2e})\n", stats.max_Dy_error, stats.avg_Dy_error );

    fmt::print( "\n  ∂²f/∂x² comparison:\n" );
    fmt::print( "    Max diff: {:12.2e} at ({:.4f}, {:.4f})\n", stats.max_Dxx_error, stats.max_Dxx_x, stats.max_Dxx_y );
    fmt::print( "    Spline value:  {:.6e}\n", stats.spline_Dxx_at_max );
    fmt::print( "    Autodiff value: {:.6e}\n", stats.ad_Dxx_at_max );
    fmt::print( "    Difference:    {:.2e} (avg: {:.2e})\n", stats.max_Dxx_error, stats.avg_Dxx_error );

    fmt::print( "\n  ∂²f/∂y² comparison:\n" );
    fmt::print( "    Max diff: {:12.2e} at ({:.4f}, {:.4f})\n", stats.max_Dyy_error, stats.max_Dyy_x, stats.max_Dyy_y );
    fmt::print( "    Spline value:  {:.6e}\n", stats.spline_Dyy_at_max );
    fmt::print( "    Autodiff value: {:.6e}\n", stats.ad_Dyy_at_max );
    fmt::print( "    Difference:    {:.2e} (avg: {:.2e})\n", stats.max_Dyy_error, stats.avg_Dyy_error );

    fmt::print( "\n  ∂²f/∂x∂y comparison:\n" );
    fmt::print( "    Max diff: {:12.2e} at ({:.4f}, {:.4f})\n", stats.max_Dxy_error, stats.max_Dxy_x, stats.max_Dxy_y );
    fmt::print( "    Spline value:  {:.6e}\n", stats.spline_Dxy_at_max );
    fmt::print( "    Autodiff value: {:.6e}\n", stats.ad_Dxy_at_max );
    fmt::print( "    Difference:    {:.2e} (avg: {:.2e})\n", stats.max_Dxy_error, stats.avg_Dxy_error );

    // Valutazione complessiva
    fmt::print( "\n  Overall assessment:\n" );

    bool all_good = true;

    if ( stats.max_Dx_error < 1e-12 ) { fmt::print( fg( fmt::color::green ), "    ✓ ∂f/∂x: PERFECT match\n" ); }
    else if ( stats.max_Dx_error < 1e-9 )
    {
      fmt::print( fg( fmt::color::green ), "    ✓ ∂f/∂x: EXCELLENT consistency\n" );
    }
    else if ( stats.max_Dx_error < 1e-6 ) { fmt::print( fg( fmt::color::yellow ), "    ~ ∂f/∂x: GOOD consistency\n" ); }
    else if ( stats.max_Dx_error < 1e-3 ) { fmt::print( fg( fmt::color::orange ), "    ~ ∂f/∂x: FAIR consistency\n" ); }
    else
    {
      fmt::print( fg( fmt::color::red ), "    ✗ ∂f/∂x: POOR consistency\n" );
      all_good = false;
    }

    if ( stats.max_Dy_error < 1e-12 ) { fmt::print( fg( fmt::color::green ), "    ✓ ∂f/∂y: PERFECT match\n" ); }
    else if ( stats.max_Dy_error < 1e-9 )
    {
      fmt::print( fg( fmt::color::green ), "    ✓ ∂f/∂y: EXCELLENT consistency\n" );
    }
    else if ( stats.max_Dy_error < 1e-6 ) { fmt::print( fg( fmt::color::yellow ), "    ~ ∂f/∂y: GOOD consistency\n" ); }
    else if ( stats.max_Dy_error < 1e-3 ) { fmt::print( fg( fmt::color::orange ), "    ~ ∂f/∂y: FAIR consistency\n" ); }
    else
    {
      fmt::print( fg( fmt::color::red ), "    ✗ ∂f/∂y: POOR consistency\n" );
      all_good = false;
    }

    if ( stats.max_Dxx_error < 1e-10 ) { fmt::print( fg( fmt::color::green ), "    ✓ ∂²f/∂x²: PERFECT match\n" ); }
    else if ( stats.max_Dxx_error < 1e-7 )
    {
      fmt::print( fg( fmt::color::green ), "    ✓ ∂²f/∂x²: EXCELLENT consistency\n" );
    }
    else if ( stats.max_Dxx_error < 1e-4 )
    {
      fmt::print( fg( fmt::color::yellow ), "    ~ ∂²f/∂x²: GOOD consistency\n" );
    }
    else if ( stats.max_Dxx_error < 1e-1 )
    {
      fmt::print( fg( fmt::color::orange ), "    ~ ∂²f/∂x²: FAIR consistency\n" );
    }
    else
    {
      fmt::print( fg( fmt::color::red ), "    ✗ ∂²f/∂x²: POOR consistency\n" );
      all_good = false;
    }

    if ( stats.max_Dyy_error < 1e-10 ) { fmt::print( fg( fmt::color::green ), "    ✓ ∂²f/∂y²: PERFECT match\n" ); }
    else if ( stats.max_Dyy_error < 1e-7 )
    {
      fmt::print( fg( fmt::color::green ), "    ✓ ∂²f/∂y²: EXCELLENT consistency\n" );
    }
    else if ( stats.max_Dyy_error < 1e-4 )
    {
      fmt::print( fg( fmt::color::yellow ), "    ~ ∂²f/∂y²: GOOD consistency\n" );
    }
    else if ( stats.max_Dyy_error < 1e-1 )
    {
      fmt::print( fg( fmt::color::orange ), "    ~ ∂²f/∂y²: FAIR consistency\n" );
    }
    else
    {
      fmt::print( fg( fmt::color::red ), "    ✗ ∂²f/∂y²: POOR consistency\n" );
      all_good = false;
    }

    if ( stats.max_Dxy_error < 1e-9 )
    {
      fmt::print( fg( fmt::color::green ), "    ✓ ∂²f/∂x∂y: EXCELLENT consistency\n" );
    }
    else if ( stats.max_Dxy_error < 1e-6 )
    {
      fmt::print( fg( fmt::color::yellow ), "    ~ ∂²f/∂x∂y: GOOD consistency\n" );
    }
    else if ( stats.max_Dxy_error < 1e-3 )
    {
      fmt::print( fg( fmt::color::orange ), "    ~ ∂²f/∂x∂y: FAIR consistency\n" );
    }
    else
    {
      fmt::print( fg( fmt::color::red ), "    ✗ ∂²f/∂x∂y: POOR consistency\n" );
      all_good = false;
    }

    if ( all_good ) { fmt::print( fg( fmt::color::green ), "\n  ✅ All autodiff consistency tests passed!\n" ); }
    else
    {
      fmt::print( fg( fmt::color::yellow ), "\n  ⚠️  Some autodiff consistency issues detected\n" );
    }
  }

}  // namespace SplinesTest

int main()
{
  using namespace SplinesTest;

  fmt::print(
    fg( fmt::color::light_blue ) | fmt::emphasis::bold,
    "\n"
    "╔══════════════════════════════════════════════════════════╗\n"
    "║            2D Spline Library - Autodiff Test             ║\n"
    "║          (Comparing spline derivatives with AD)          ║\n"
    "╚══════════════════════════════════════════════════════════╝\n\n" );

  // Create a moderately sized uniform grid using Eigen
  integer nx = 51;  // Reduced for faster testing
  integer ny = 63;

  // Domain for trigonometric function (gentle oscillations)
  real_type XMIN = -m_pi;
  real_type XMAX = m_pi;
  real_type YMIN = -m_pi;
  real_type YMAX = m_pi;

  // Create uniform grid using Eigen
  VectorXd x_grid = VectorXd::LinSpaced( nx, XMIN, XMAX );
  VectorXd y_grid = VectorXd::LinSpaced( ny, YMIN, YMAX );

  // Generate data using trigonometric function (well-behaved)
  MatrixXd z_data = MatrixXd::NullaryExpr(
    nx,
    ny,
    [&]( Eigen::Index i, Eigen::Index j ) { return trig_function( x_grid( i ), y_grid( j ) ); } );

  fmt::print(
    "Test configuration:\n"
    "  Grid: {} × {} points\n"
    "  Domain: x ∈ [{:.2f}, {:.2f}], y ∈ [{:.2f}, {:.2f}]\n"
    "  Test function: sin(0.5x) * cos(0.5y) (smooth, well-behaved)\n"
    "  Test points: 500 per spline type\n\n",
    nx,
    ny,
    x_grid.minCoeff(),
    x_grid.maxCoeff(),
    y_grid.minCoeff(),
    y_grid.maxCoeff() );

  // Test all spline types
  vector<SplineType2D> spline_types = { SplineType2D::BILINEAR,        SplineType2D::BICUBIC_CUBIC,
                                        SplineType2D::BICUBIC_AKIMA,   SplineType2D::BICUBIC_VANLEER,
                                        SplineType2D::BICUBIC_PCHIP,   SplineType2D::BIQUINTIC_CUBIC,
                                        SplineType2D::BIQUINTIC_AKIMA, SplineType2D::BIQUINTIC_VANLEER,
                                        SplineType2D::BIQUINTIC_PCHIP };

  for ( auto type : spline_types ) { test_spline_type( type, x_grid, y_grid, z_data, 500 ); }

  // Test interpolation accuracy
  test_interpolation_accuracy( x_grid, y_grid, z_data );

  // Test autodiff consistency
  Spline2D test_spline( "AutodiffConsistencyTest" );

  test_spline
    .build( SplineType2D::BIQUINTIC_CUBIC, x_grid.data(), 1, y_grid.data(), 1, z_data.data(), nx, nx, ny, true, false );
  test_autodiff_consistency( test_spline, x_grid, y_grid );

  // Final summary
  fmt::print(
    fg( fmt::color::light_green ) | fmt::emphasis::bold,
    "\n\n"
    "════════════════════════════════════════════════════════════\n"
    "                    TEST COMPLETED                          \n"
    "════════════════════════════════════════════════════════════\n" );

  fmt::print(
    "\nInterpretation guide:\n"
    "  ✅ Excellent (error < 1e-6)\n"
    "  ⚠️  Acceptable (error < 1e-3)\n"
    "  ❌ Problematic (error ≥ 1e-3)\n"
    "\n"
    "Note: Autodiff (AD) consistency errors measure how well\n"
    "the spline's analytical derivatives match autodiff derivatives\n"
    "of the spline itself. This is a key consistency check.\n"
    "\n"
    "Grid interpolation test: Checks if the spline exactly matches\n"
    "the input data at the grid points (as expected for interpolating splines).\n" );

  return 0;
}
