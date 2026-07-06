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

#include "Splines/Splines.hh"
#include "Utils_fmt.hh"

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace SplinesLoad;
using namespace std;

using Splines::integer;
using Splines::real_type;
using Splines::Spline_sub_type;

namespace {

  auto const HEADER_COLOR  = fg( fmt::color::steel_blue ) | fmt::emphasis::bold;
  auto const OK_COLOR      = fg( fmt::color::lime_green ) | fmt::emphasis::bold;
  auto const FAIL_COLOR    = fg( fmt::color::crimson ) | fmt::emphasis::bold;
  auto const INFO_COLOR    = fg( fmt::color::deep_sky_blue );
  constexpr real_type VALUE_TOL      = 1e-12;
  constexpr real_type DERIVATIVE_TOL = 5e-10;

  real_type reference_surface( real_type x, real_type y )
  {
    return 0.35 * sin( 1.3 * x - 0.2 * y ) +
           0.25 * cos( 0.8 * x + 1.1 * y ) +
           0.08 * x * x -
           0.05 * y * y +
           0.12 * x * y;
  }

  vector<real_type> make_x_grid()
  { return { -1.25, -0.5, 0.0, 0.8, 1.7, 2.5 }; }

  vector<real_type> make_y_grid()
  { return { -2.0, -1.1, -0.2, 0.6, 1.4 }; }

  vector<real_type> make_z_grid( vector<real_type> const & x, vector<real_type> const & y )
  {
    vector<real_type> z( x.size() * y.size() );
    for ( integer ix = 0; ix < integer( x.size() ); ++ix )
      for ( integer iy = 0; iy < integer( y.size() ); ++iy )
        z[size_t( ix ) * y.size() + size_t( iy )] = reference_surface( x[size_t( ix )], y[size_t( iy )] );
    return z;
  }

  struct Stats {
    real_type max_value_diff = 0;
    real_type max_dx_diff    = 0;
    real_type max_dy_diff    = 0;
    real_type max_dxx_diff   = 0;
    real_type max_dxy_diff   = 0;
    real_type max_dyy_diff   = 0;
    real_type min_value      = std::numeric_limits<real_type>::max();
    real_type max_value      = std::numeric_limits<real_type>::lowest();
  };

  template <typename SplineBuilder>
  bool run_test_case( string const & name, SplineBuilder && make_spline ) {
    auto const x = make_x_grid();
    auto const y = make_y_grid();
    auto const z = make_z_grid( x, y );

    auto spline_from_vectors  = make_spline();
    auto spline_from_callable = make_spline();

    spline_from_vectors->build( x, y, z, false, false );
    spline_from_callable->build(
      [&]( integer i ) -> real_type { return x[size_t( i )]; },
      [&]( integer j ) -> real_type { return y[size_t( j )]; },
      [&]( integer ix, integer iy ) -> real_type { return z[size_t( ix ) * y.size() + size_t( iy )]; },
      integer( x.size() ),
      integer( y.size() )
    );

    Stats st;

    for ( integer ix = 0; ix < integer( x.size() ); ++ix ) {
      for ( integer iy = 0; iy < integer( y.size() ); ++iy ) {
        real_type const xv = x[size_t( ix )];
        real_type const yv = y[size_t( iy )];
        real_type const vv = (*spline_from_vectors)( xv, yv );
        st.max_value_diff  = max( st.max_value_diff, abs( vv - (*spline_from_callable)( xv, yv ) ) );
        st.min_value       = min( st.min_value, vv );
        st.max_value       = max( st.max_value, vv );
      }
    }

    vector<pair<real_type, real_type>> const samples = {
      { -0.9, -1.4 }, { -0.2, -0.3 }, { 0.4, 0.1 }, { 1.1, 0.9 }, { 1.9, 0.2 }
    };

    for ( auto const & [xv, yv] : samples ) {
      real_type const vv = spline_from_vectors->eval( xv, yv );
      st.max_value_diff  = max( st.max_value_diff, abs( vv - spline_from_callable->eval( xv, yv ) ) );
      st.max_dx_diff    = max( st.max_dx_diff, abs( spline_from_vectors->Dx( xv, yv ) - spline_from_callable->Dx( xv, yv ) ) );
      st.max_dy_diff    = max( st.max_dy_diff, abs( spline_from_vectors->Dy( xv, yv ) - spline_from_callable->Dy( xv, yv ) ) );
      st.max_dxx_diff   = max( st.max_dxx_diff, abs( spline_from_vectors->Dxx( xv, yv ) - spline_from_callable->Dxx( xv, yv ) ) );
      st.max_dxy_diff   = max( st.max_dxy_diff, abs( spline_from_vectors->Dxy( xv, yv ) - spline_from_callable->Dxy( xv, yv ) ) );
      st.max_dyy_diff   = max( st.max_dyy_diff, abs( spline_from_vectors->Dyy( xv, yv ) - spline_from_callable->Dyy( xv, yv ) ) );
      st.min_value      = min( st.min_value, vv );
      st.max_value      = max( st.max_value, vv );
    }

    bool const ok =
      st.max_value_diff <= VALUE_TOL &&
      st.max_dx_diff    <= DERIVATIVE_TOL &&
      st.max_dy_diff    <= DERIVATIVE_TOL &&
      st.max_dxx_diff   <= DERIVATIVE_TOL &&
      st.max_dxy_diff   <= DERIVATIVE_TOL &&
      st.max_dyy_diff   <= DERIVATIVE_TOL;

    fmt::print(
      ok ? OK_COLOR : FAIL_COLOR,
      "{:<20} range=[{:.3e},{:.3e}] value={:.3e} dx={:.3e} dy={:.3e} dxx={:.3e} dxy={:.3e} dyy={:.3e}\n",
      name,
      st.min_value,
      st.max_value,
      st.max_value_diff,
      st.max_dx_diff,
      st.max_dy_diff,
      st.max_dxx_diff,
      st.max_dxy_diff,
      st.max_dyy_diff
    );

    return ok;
  }

}  // namespace

int
main() {
  fmt::print( HEADER_COLOR, "\n=== test15: generic callable build for surface splines ===\n\n" );
  fmt::print( INFO_COLOR, "Comparing build-from-vectors against build-from-callables\n\n" );

  bool ok = true;

  ok = run_test_case(
         "Bilinear",
         []() { return std::make_unique<Splines::BilinearSpline>(); }
       ) && ok;

  ok = run_test_case(
         "BiCubic[CUBIC]",
         []() { return std::make_unique<Splines::BiCubicSpline>( Spline_sub_type::CUBIC ); }
       ) && ok;

  ok = run_test_case(
         "BiCubic[AKIMA]",
         []() { return std::make_unique<Splines::BiCubicSpline>( Spline_sub_type::AKIMA ); }
       ) && ok;

  ok = run_test_case(
         "BiQuintic[CUBIC]",
         []() { return std::make_unique<Splines::BiQuinticSpline>( Spline_sub_type::CUBIC ); }
       ) && ok;

  if ( ok ) {
    fmt::print( OK_COLOR, "\nAll callable-build tests passed.\n" );
    return 0;
  }

  fmt::print( FAIL_COLOR, "\nCallable-build regression detected.\n" );
  return 1;
}
