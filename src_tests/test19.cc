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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// test19: the two axes of a surface spline carry their own units.
//
// A tensor product spline is INVARIANT under an affine change of either
// coordinate: built on (x, y, Z) and on (x, sigma*y, Z) it must describe the
// same surface, with the y-derivative simply divided by sigma.  Nothing in the
// construction may look at one axis while working along the other.
//
// `SplineSurf::make_derivative_y` did exactly that: it built the 1D spline of
// each ROW of Z -- a row runs along y -- on the abscissas `mX`.  The
// y-derivatives came out multiplied by the ratio of the two grid steps, so the
// surface still interpolated every node exactly and oscillated by orders of
// magnitude in between.  It stayed hidden because every test grid had the two
// axes similarly spaced, and because `Spline::build(x,y)` silently takes
// min(|x|,|y|) points: with nx != ny the rows were truncated or read past the
// end instead of raising.
//
// The table that exposed it is here as a test case: thrust of the Bryson
// supersonic aircraft against Mach number (step 0.2) and altitude in METRES
// (step 1524, and not uniform at the top), from section 6.4 of the GPOPS-II
// User's Guide.  Ratio of the steps: 7620.  Evaluated along a climb the
// surface swung between -2.7e7 and +5.0e7 N on a table whose values run from
// 445 to 2.0e5 N.
//

#include "Splines.hh"
#include "Utils_fmt.hh"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace std;

using Splines::integer;
using Splines::real_type;
using Splines::Spline_sub_type;
using Splines::SplineSurf;

namespace {

  auto const HEADER_COLOR = fg( fmt::color::steel_blue ) | fmt::emphasis::bold;
  auto const OK_COLOR     = fg( fmt::color::lime_green ) | fmt::emphasis::bold;
  auto const FAIL_COLOR   = fg( fmt::color::crimson ) | fmt::emphasis::bold;
  auto const INFO_COLOR   = fg( fmt::color::deep_sky_blue );

  //! Interpolation at the nodes, relative to the largest |z|.
  constexpr real_type NODE_TOL = 1e-9;
  //! Agreement between a surface and its rescaled twin, relative.
  constexpr real_type SCALE_TOL = 1e-9;
  //! How far outside the data range the interpolant may wander, as a fraction
  //! of that range.  Deliberately generous -- a cubic through a steep table is
  //! allowed its overshoot; what this catches is the 250x of a derivative built
  //! on the wrong axis.
  constexpr real_type RANGE_SLACK = 1.0;

  using Factory = function<unique_ptr<SplineSurf>()>;

  /*\
   |   ____        _
   |  |  _ \  __ _| |_ __ _
   |  | | | |/ _` | __/ _` |
   |  | |_| | (_| | || (_| |
   |  |____/ \__,_|\__\__,_|
  \*/

  //! Mach number: uniform step 0.2.
  vector<real_type> thrust_mach() {
    return {
      0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8
    };
  }

  //! Altitude [m]: step 1524 six times, then 3048, 3048, 6096.
  vector<real_type> thrust_altitude() {
    return {
      0, 1524, 3048, 4572, 6096, 7620, 9144, 12192, 15240, 21336
    };
  }

  //! Thrust [N], z[ix*ny+iy] at Mach = mach[ix], altitude = altitude[iy].
  vector<real_type> thrust_table() {
    return {
      107646.9724, 106757.328, 90298.9066, 76954.2406, 64499.219,
      54268.3084, 45371.8644, 25354.8654, 15123.9548, 444.8222,
      124550.216, 109426.2612, 93857.4842, 80512.8182, 67612.9744,
      56937.2416, 47595.9754, 28913.443, 17348.0658, 889.6444,
      125884.6826, 112095.1944, 97416.0618, 83181.7514, 70726.7298,
      59606.1748, 49820.0864, 32472.0206, 19572.1768, 1779.2888,
      137005.2376, 120991.6384, 105867.6836, 91188.551, 76954.2406,
      65388.8634, 54713.1306, 36030.5982, 21796.2878, 3558.5776,
      153463.659, 134781.1266, 118322.7052, 103198.7504, 88074.7956,
      74730.1296, 62719.9302, 41813.2868, 24910.0432, 4893.0442,
      168587.6138, 152574.0146, 135225.9488, 119212.3496, 103643.5726,
      88074.7956, 74730.1296, 49820.0864, 30247.9096, 6227.5108,
      160580.8142, 169032.436, 155242.9478, 139229.3486, 121436.4606,
      104978.0392, 89409.2622, 59606.1748, 36920.2426, 7561.9774,
      160580.8142, 162804.9252, 171256.547, 160580.8142, 140563.8152,
      124995.0382, 107646.9724, 72061.1964, 44482.22, 9786.0884,
      160580.8142, 156577.4144, 187270.1462, 172146.1914, 158801.5254,
      142343.104, 124995.0382, 85850.6846, 52933.8418, 12899.8438,
      160580.8142, 150349.9036, 203283.7454, 183711.5686, 177039.2356,
      153908.4812, 138339.7042, 96526.4174, 59161.3526, 13789.4882
    };
  }

  //! A smooth surface with no particular symmetry between the two axes.
  real_type reference_surface( real_type const x, real_type const y ) {
    return 0.35 * sin( 1.3 * x - 0.2 * y ) + 0.25 * cos( 0.8 * x + 1.1 * y ) + 0.08 * x * x -
           0.05 * y * y + 0.12 * x * y;
  }

  //! Deliberately NON uniform, so a bug cannot hide behind a constant step.
  vector<real_type> stretched_grid( real_type const a, real_type const b, integer const n ) {
    vector<real_type> v( static_cast<size_t>( n ) );
    for ( integer i = 0; i < n; ++i ) {
      real_type const t = static_cast<real_type>( i ) / static_cast<real_type>( n - 1 );
      v[static_cast<size_t>( i )] = a + ( b - a ) * t * ( 0.65 + 0.35 * t );
    }
    return v;
  }

  vector<real_type> sampled_surface( vector<real_type> const & x, vector<real_type> const & y ) {
    vector<real_type> z( x.size() * y.size() );
    for ( size_t ix = 0; ix < x.size(); ++ix )
      for ( size_t iy = 0; iy < y.size(); ++iy )
        z[ix * y.size() + iy] = reference_surface( x[ix], y[iy] );
    return z;
  }

  //! Interior points of every cell: this is where a wrong derivative shows.
  vector<pair<real_type, real_type>> cell_samples( vector<real_type> const & x, vector<real_type> const & y ) {
    vector<pair<real_type, real_type>> s;
    for ( size_t i = 0; i + 1 < x.size(); ++i ) {
      for ( size_t j = 0; j + 1 < y.size(); ++j ) {
        for ( real_type const tx : { 0.25, 0.5, 0.75 } ) {
          for ( real_type const ty : { 0.25, 0.5, 0.75 } ) {
            s.emplace_back( x[i] + tx * ( x[i + 1] - x[i] ), y[j] + ty * ( y[j + 1] - y[j] ) );
          }
        }
      }
    }
    return s;
  }

  real_type max_abs( vector<real_type> const & v ) {
    real_type m = 0;
    for ( real_type const z : v ) m = max( m, abs( z ) );
    return m > 0 ? m : 1;
  }

  struct Case {
    string            name;
    vector<real_type> x, y, z;
  };

  /*\
   |    ____ _               _
   |   / ___| |__   ___  ___| | _____
   |  | |   | '_ \ / _ \/ __| |/ / __|
   |  | |___| | | |  __/ (__|   <\__ \
   |   \____|_| |_|\___|\___|_|\_\___/
  \*/

  //! Every surface spline interpolates: this passes even with a broken
  //! derivative, and is here to prove that a failure below is NOT bad data.
  real_type node_error( Factory const & make, Case const & c ) {
    auto T = make();
    T->build( c.x, c.y, c.z );
    real_type const scale = max_abs( c.z );
    real_type       worst = 0;
    for ( size_t i = 0; i < c.x.size(); ++i )
      for ( size_t j = 0; j < c.y.size(); ++j )
        worst = max( worst, abs( T->eval( c.x[i], c.y[j] ) - c.z[i * c.y.size() + j] ) / scale );
    return worst;
  }

  //! Largest excursion outside the data range, as a fraction of that range.
  //! Zero when the interpolant stays inside it.
  real_type overshoot( Factory const & make, Case const & c ) {
    auto T = make();
    T->build( c.x, c.y, c.z );
    real_type zmin = numeric_limits<real_type>::max();
    real_type zmax = numeric_limits<real_type>::lowest();
    for ( real_type const z : c.z ) { zmin = min( zmin, z ); zmax = max( zmax, z ); }
    real_type const range = zmax > zmin ? zmax - zmin : 1;
    real_type       worst = 0;
    for ( auto const & [xv, yv] : cell_samples( c.x, c.y ) ) {
      real_type const v = T->eval( xv, yv );
      if ( v > zmax ) worst = max( worst, ( v - zmax ) / range );
      if ( v < zmin ) worst = max( worst, ( zmin - v ) / range );
    }
    return worst;
  }

  //! Rebuild the same surface with one axis multiplied by `sigma` and compare.
  //! The value must not move and the derivative along that axis must scale by
  //! 1/sigma -- a property of the spline, with no reference implementation to
  //! drift against.
  pair<real_type, real_type> axis_scaling_error(
    Factory const & make, Case const & c, real_type const sigma, bool const on_y
  ) {
    auto S1 = make();
    S1->build( c.x, c.y, c.z );

    vector<real_type> x2 = c.x, y2 = c.y;
    for ( real_type & v : ( on_y ? y2 : x2 ) ) v *= sigma;
    auto S2 = make();
    S2->build( x2, y2, c.z );

    real_type const scale = max_abs( c.z );
    real_type       worst_v = 0, worst_d = 0, dscale = 0;
    for ( auto const & [xv, yv] : cell_samples( c.x, c.y ) ) {
      real_type const X2 = on_y ? xv : sigma * xv;
      real_type const Y2 = on_y ? sigma * yv : yv;
      worst_v = max( worst_v, abs( S1->eval( xv, yv ) - S2->eval( X2, Y2 ) ) / scale );
      real_type const d1 = on_y ? S1->Dy( xv, yv ) : S1->Dx( xv, yv );
      real_type const d2 = on_y ? S2->Dy( X2, Y2 ) : S2->Dx( X2, Y2 );
      dscale  = max( dscale, abs( d1 ) );
      worst_d = max( worst_d, abs( d1 - sigma * d2 ) );
    }
    return { worst_v, dscale > 0 ? worst_d / dscale : worst_d };
  }

  bool run( string const & spline_name, Factory const & make, Case const & c ) {
    real_type const nodes = node_error( make, c );
    real_type const over  = overshoot( make, c );
    auto const [vy, dy]   = axis_scaling_error( make, c, 1.0e-3, true );
    auto const [vx, dx]   = axis_scaling_error( make, c, 1.0e+3, false );

    bool const ok = nodes <= NODE_TOL && over <= RANGE_SLACK && vy <= SCALE_TOL && dy <= SCALE_TOL &&
                    vx <= SCALE_TOL && dx <= SCALE_TOL;

    fmt::print(
      ok ? OK_COLOR : FAIL_COLOR,
      "  {:<18} {:<22} nodes={:8.2e}  overshoot={:8.2e}  y-rescale=({:8.2e},{:8.2e})  "
      "x-rescale=({:8.2e},{:8.2e})\n",
      spline_name, c.name, nodes, over, vy, dy, vx, dx
    );
    return ok;
  }

}  // namespace

int
main() {
  fmt::print( HEADER_COLOR, "\n=== test19: surface splines, one unit per axis ===\n\n" );
  fmt::print(
    INFO_COLOR,
    "nodes     : interpolation at the nodes (any bug in the derivatives passes this)\n"
    "overshoot : largest excursion outside the data range, in units of that range\n"
    "rescale   : (value, derivative) change after multiplying one axis by a constant\n\n"
  );

  vector<Case> cases;
  {
    auto const mach = thrust_mach();
    auto const alt  = thrust_altitude();
    cases.push_back( { "GPOPS2 thrust 10x10", mach, alt, thrust_table() } );
  }
  for ( auto const & [nx, ny] : vector<pair<integer, integer>>{ { 7, 7 }, { 6, 11 }, { 11, 6 } } ) {
    auto const x = stretched_grid( -1.25, 2.5, nx );
    auto const y = stretched_grid( -2.0, 1.4, ny );
    cases.push_back( { fmt::format( "analytic {}x{}", nx, ny ), x, y, sampled_surface( x, y ) } );
  }

  vector<pair<string, Factory>> const splines{
    { "Bilinear", []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BilinearSpline>(); } },
    { "BiCubic[CUBIC]",
      []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BiCubicSpline>( Spline_sub_type::CUBIC ); } },
    { "BiCubic[AKIMA]",
      []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BiCubicSpline>( Spline_sub_type::AKIMA ); } },
    { "BiCubic[PCHIP]",
      []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BiCubicSpline>( Spline_sub_type::PCHIP ); } },
    { "BiQuintic[CUBIC]",
      []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BiQuinticSpline>( Spline_sub_type::CUBIC ); } }
  };

  bool ok = true;
  for ( auto const & c : cases ) {
    fmt::print( INFO_COLOR, "{}\n", c.name );
    for ( auto const & [name, make] : splines ) ok = run( name, make, c ) && ok;
    fmt::print( "\n" );
  }

  if ( ok ) {
    fmt::print( OK_COLOR, "All surface spline tests passed.\n" );
    return 0;
  }

  fmt::print(
    FAIL_COLOR,
    "Surface spline regression detected.\n"
    "A failing `y-rescale` with `nodes` and `x-rescale` clean means the y-derivatives\n"
    "are being built on the x abscissas: see SplineSurf::make_derivative_y.\n"
  );
  return 1;
}
