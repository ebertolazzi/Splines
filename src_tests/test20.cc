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
// test20: value and every derivative of a surface spline, checked four ways.
//
// A surface spline exposes six numbers at each point -- value, Dx, Dy, Dxx,
// Dxy, Dyy -- through three different entry points (the scalar accessors, the
// packed `D` and the packed `DD`).  Each is a separate piece of arithmetic, and
// a slip in any of them is invisible to a test that only looks at the value at
// the nodes.  Four independent checks, none of which needs a reference
// implementation to compare against:
//
//   exact    every family reproduces BILINEAR data exactly, so the six numbers
//            are known in closed form: Dx = B + D*y, Dy = C + D*x, Dxy = D and
//            the two pure second derivatives vanish.  This pins units, signs
//            and the x/y roles at once.
//   fd       on generic data, each derivative against a Richardson extrapolated
//            central difference of the spline's OWN value, with every stencil
//            kept inside one cell (across a cell edge the second derivative of
//            a bicubic jumps, and the comparison would be meaningless).  Dxy is
//            taken both as d(Dy)/dx and as d(Dx)/dy.
//   cov      covariance under rescaling of one axis: x -> sigma*x turns Dx into
//            Dx/sigma, Dxx into Dxx/sigma^2, Dxy into Dxy/sigma and leaves Dy
//            and Dyy alone.  Order by order, this is what catches a derivative
//            built on the wrong abscissas.  \see test19
//   packed   `D` and `DD` against the scalar accessors: same point, same
//            numbers, three code paths.
//
// The grids are deliberately non uniform, and one of them is the Mach number by
// altitude in METRES of the GPOPS-II thrust table: two axes four orders of
// magnitude apart, which is the configuration under which a derivative built on
// the wrong axis stops being a rounding detail.
//

#include "Splines.hh"
#include "Utils_fmt.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <span>
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

  //! Reproduction of bilinear data: exact but for rounding.
  constexpr real_type EXACT_TOL = 1e-10;
  //! Relative part of the finite difference criterion.  The criterion itself is
  //! mixed -- see `check_finite_differences` -- and its margin is what gets
  //! compared against one.
  constexpr real_type FD_REL = 1e-8;
  //! How many rounding units of the difference quotient are tolerated on a
  //! derivative whose exact value is zero.
  constexpr real_type FD_NOISE = 100;
  //! Covariance under rescaling of an axis: a property, not an approximation.
  constexpr real_type COV_TOL = 1e-9;
  //! `D`/`DD` against the scalar accessors: the same arithmetic twice.
  constexpr real_type PACKED_TOL = 1e-12;

  using Factory = function<unique_ptr<SplineSurf>()>;

  /*\
   |    ____      _     _
   |   / ___|_ __(_) __| |___
   |  | |  _| '__| |/ _` / __|
   |  | |_| | |  | | (_| \__ \
   |   \____|_|  |_|\__,_|___/
  \*/

  struct Grid {
    string            name;
    vector<real_type> x, y;
  };

  //! Non uniform on purpose: a constant step hides an index slip.
  vector<real_type> stretched( real_type const a, real_type const b, integer const n ) {
    vector<real_type> v( static_cast<size_t>( n ) );
    for ( integer i = 0; i < n; ++i ) {
      real_type const t = static_cast<real_type>( i ) / static_cast<real_type>( n - 1 );
      v[static_cast<size_t>( i )] = a + ( b - a ) * t * ( 0.65 + 0.35 * t );
    }
    return v;
  }

  vector<Grid> grids() {
    return {
      { "unit 7x7", stretched( -1.25, 2.5, 7 ), stretched( -2.0, 1.4, 7 ) },
      { "unit 6x11", stretched( -1.25, 2.5, 6 ), stretched( -2.0, 1.4, 11 ) },
      { "unit 11x6", stretched( -1.25, 2.5, 11 ), stretched( -2.0, 1.4, 6 ) },
      // Mach number by altitude [m], the axes of the GPOPS-II thrust table.
      { "Mach x m 10x10",
        { 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 },
        { 0, 1524, 3048, 4572, 6096, 7620, 9144, 12192, 15240, 21336 } }
    };
  }

  real_type max_abs( vector<real_type> const & v ) {
    real_type m = 0;
    for ( real_type const z : v ) m = max( m, abs( z ) );
    return m > 0 ? m : 1;
  }

  real_type min_step( vector<real_type> const & v ) {
    real_type m = numeric_limits<real_type>::max();
    for ( size_t i = 0; i + 1 < v.size(); ++i ) m = min( m, v[i + 1] - v[i] );
    return m;
  }

  //! A smooth surface of the NORMALIZED coordinates, so that the data stays of
  //! order one whatever the units of the two axes.
  real_type generic_surface( real_type const u, real_type const v ) {
    return 0.35 * sin( 1.3 * u - 0.2 * v ) + 0.25 * cos( 0.8 * u + 1.1 * v ) + 0.08 * u * u -
           0.05 * v * v + 0.12 * u * v;
  }

  vector<real_type> generic_data( Grid const & g ) {
    real_type const sx = max_abs( g.x ), sy = max_abs( g.y );
    vector<real_type> z( g.x.size() * g.y.size() );
    for ( size_t i = 0; i < g.x.size(); ++i )
      for ( size_t j = 0; j < g.y.size(); ++j )
        z[i * g.y.size() + j] = generic_surface( 3 * g.x[i] / sx, 3 * g.y[j] / sy );
    return z;
  }

  //! Interior points of every cell, as a fraction of the cell.  Kept away from
  //! the edges so that a finite difference stencil never crosses a node.
  vector<pair<real_type, real_type>> cell_samples( Grid const & g, real_type const lo, real_type const hi ) {
    vector<pair<real_type, real_type>> s;
    for ( size_t i = 0; i + 1 < g.x.size(); ++i ) {
      for ( size_t j = 0; j + 1 < g.y.size(); ++j ) {
        for ( real_type const tx : { lo, 0.5, hi } ) {
          for ( real_type const ty : { lo, 0.5, hi } ) {
            s.emplace_back( g.x[i] + tx * ( g.x[i + 1] - g.x[i] ), g.y[j] + ty * ( g.y[j + 1] - g.y[j] ) );
          }
        }
      }
    }
    return s;
  }

  //! Worst relative error of a set of (got, expected, scale) comparisons.
  struct Worst {
    real_type value = 0;
    void add( real_type const got, real_type const want, real_type const scale ) {
      value = max( value, abs( got - want ) / ( scale > 0 ? scale : 1 ) );
    }
  };

  /*\
   |    ____ _               _
   |   / ___| |__   ___  ___| | _____
   |  | |   | '_ \ / _ \/ __| |/ / __|
   |  | |___| | | |  __/ (__|   <\__ \
   |   \____|_| |_|\___|\___|_|\_\___/
  \*/

  //! Bilinear data is reproduced exactly by every family, so all six numbers
  //! are known in closed form.
  real_type check_exact( Factory const & make, Grid const & g ) {
    constexpr real_type A = 3.7, B = -1.9, C = 0.6, D = 2.3;

    vector<real_type> z( g.x.size() * g.y.size() );
    for ( size_t i = 0; i < g.x.size(); ++i )
      for ( size_t j = 0; j < g.y.size(); ++j )
        z[i * g.y.size() + j] = A + B * g.x[i] + C * g.y[j] + D * g.x[i] * g.y[j];

    auto S = make();
    S->build( g.x, g.y, z );

    real_type const Xm = max_abs( g.x ), Ym = max_abs( g.y );
    real_type const s_v  = abs( A ) + abs( B ) * Xm + abs( C ) * Ym + abs( D ) * Xm * Ym;
    real_type const s_dx = abs( B ) + abs( D ) * Ym;
    real_type const s_dy = abs( C ) + abs( D ) * Xm;
    // The scale a spurious pure second derivative would have.
    real_type const s_dxx = s_dx / min_step( g.x );
    real_type const s_dyy = s_dy / min_step( g.y );

    Worst w;
    for ( auto const & [xv, yv] : cell_samples( g, 0.25, 0.75 ) ) {
      w.add( S->eval( xv, yv ), A + B * xv + C * yv + D * xv * yv, s_v );
      w.add( S->Dx( xv, yv ), B + D * yv, s_dx );
      w.add( S->Dy( xv, yv ), C + D * xv, s_dy );
      w.add( S->Dxx( xv, yv ), 0, s_dxx );
      w.add( S->Dxy( xv, yv ), D, abs( D ) );
      w.add( S->Dyy( xv, yv ), 0, s_dyy );
    }
    return w.value;
  }

  //! Every derivative against a central difference of the spline's OWN value,
  //! extrapolated twice (h, h/2, h/4): both the h^2 and the h^4 terms of the
  //! error cancel, so on a piecewise polynomial of degree <= 5 -- which every
  //! family here is -- the comparison is exact but for rounding.
  //!
  //! The criterion is MIXED and not relative.  A difference quotient divides by
  //! h, so it carries a rounding floor of its own, and a derivative that is
  //! exactly zero -- Dxx of a bilinear patch -- has nothing to be relative to:
  //! comparing it with that floor would read as a 100% error forever.  So each
  //! derivative is allowed `FD_REL` of its own magnitude PLUS `FD_NOISE` units
  //! of the floor, and what is returned is the margin: one or less passes.
  //!
  //! One number per derivative, Dx Dy Dxx Dyy Dxy: a single worst case would
  //! not say WHICH one drifted.
  array<real_type, 5> check_finite_differences( Factory const & make, Grid const & g ) {
    auto const z = generic_data( g );
    auto       S = make();
    S->build( g.x, g.y, z );

    // Every stencil, down to h/4, stays inside one cell: across a cell edge the
    // second derivative of a bicubic jumps and the comparison is meaningless.
    real_type const hx = 0.15 * min_step( g.x );
    real_type const hy = 0.15 * min_step( g.y );

    auto const extrapolate = []( auto && Dh, real_type const h ) {
      real_type const r1 = ( 4 * Dh( h / 2 ) - Dh( h ) ) / 3;
      real_type const r2 = ( 4 * Dh( h / 4 ) - Dh( h / 2 ) ) / 3;
      return ( 16 * r2 - r1 ) / 15;
    };

    real_type           f_max = 0;
    array<real_type, 5> diff{}, scale{};
    for ( auto const & [xv, yv] : cell_samples( g, 0.35, 0.65 ) ) {
      array<real_type, 5> const got{ S->Dx( xv, yv ),  S->Dy( xv, yv ), S->Dxx( xv, yv ),
                                     S->Dyy( xv, yv ), S->Dxy( xv, yv ) };
      array<real_type, 5> const fd{
        extrapolate( [&]( real_type const h ) { return ( S->eval( xv + h, yv ) - S->eval( xv - h, yv ) ) / ( 2 * h ); }, hx ),
        extrapolate( [&]( real_type const h ) { return ( S->eval( xv, yv + h ) - S->eval( xv, yv - h ) ) / ( 2 * h ); }, hy ),
        extrapolate(
          [&]( real_type const h ) {
            return ( S->eval( xv + h, yv ) - 2 * S->eval( xv, yv ) + S->eval( xv - h, yv ) ) / ( h * h );
          },
          hx ),
        extrapolate(
          [&]( real_type const h ) {
            return ( S->eval( xv, yv + h ) - 2 * S->eval( xv, yv ) + S->eval( xv, yv - h ) ) / ( h * h );
          },
          hy ),
        // The mixed one as d(Dy)/dx: it ties Dxy to Dy, which the two entries
        // above have just tied to the value.
        extrapolate( [&]( real_type const h ) { return ( S->Dy( xv + h, yv ) - S->Dy( xv - h, yv ) ) / ( 2 * h ); }, hx )
      };
      f_max = max( f_max, abs( S->eval( xv, yv ) ) );
      for ( size_t k = 0; k < 5; ++k ) {
        diff[k]  = max( diff[k], abs( got[k] - fd[k] ) );
        scale[k] = max( { scale[k], abs( got[k] ), abs( fd[k] ) } );
      }
    }

    // Rounding floor of the quotients, taken at the shortest step actually used.
    real_type const eps = numeric_limits<real_type>::epsilon();
    real_type const ax = hx / 4, ay = hy / 4;
    array<real_type, 5> const floor_{ FD_NOISE * eps * f_max / ax,
                                      FD_NOISE * eps * f_max / ay,
                                      FD_NOISE * eps * f_max / ( ax * ax ),
                                      FD_NOISE * eps * f_max / ( ay * ay ),
                                      FD_NOISE * eps * scale[1] / ax };

    array<real_type, 5> margin{};
    for ( size_t k = 0; k < 5; ++k ) margin[k] = diff[k] / ( FD_REL * scale[k] + floor_[k] );
    return margin;
  }

  //! x -> sigma*x (or y -> sigma*y) and every order rescaled accordingly.
  real_type check_covariance( Factory const & make, Grid const & g, real_type const sigma, bool const on_y ) {
    auto const z = generic_data( g );

    auto S1 = make();
    S1->build( g.x, g.y, z );

    Grid g2 = g;
    for ( real_type & v : ( on_y ? g2.y : g2.x ) ) v *= sigma;
    auto S2 = make();
    S2->build( g2.x, g2.y, z );

    real_type const s_v = max_abs( z );
    real_type       s_dx = 0, s_dy = 0, s_dxx = 0, s_dxy = 0, s_dyy = 0;
    for ( auto const & [xv, yv] : cell_samples( g, 0.25, 0.75 ) ) {
      s_dx  = max( s_dx, abs( S1->Dx( xv, yv ) ) );
      s_dy  = max( s_dy, abs( S1->Dy( xv, yv ) ) );
      s_dxx = max( s_dxx, abs( S1->Dxx( xv, yv ) ) );
      s_dxy = max( s_dxy, abs( S1->Dxy( xv, yv ) ) );
      s_dyy = max( s_dyy, abs( S1->Dyy( xv, yv ) ) );
    }

    real_type const kx = on_y ? 1 : sigma;  // factor turning S2's x-derivative back into S1's
    real_type const ky = on_y ? sigma : 1;

    Worst w;
    for ( auto const & [xv, yv] : cell_samples( g, 0.25, 0.75 ) ) {
      real_type const X2 = on_y ? xv : sigma * xv;
      real_type const Y2 = on_y ? sigma * yv : yv;
      w.add( kx * S2->Dx( X2, Y2 ), S1->Dx( xv, yv ), s_dx );
      w.add( ky * S2->Dy( X2, Y2 ), S1->Dy( xv, yv ), s_dy );
      w.add( kx * kx * S2->Dxx( X2, Y2 ), S1->Dxx( xv, yv ), s_dxx );
      w.add( kx * ky * S2->Dxy( X2, Y2 ), S1->Dxy( xv, yv ), s_dxy );
      w.add( ky * ky * S2->Dyy( X2, Y2 ), S1->Dyy( xv, yv ), s_dyy );
      w.add( S2->eval( X2, Y2 ), S1->eval( xv, yv ), s_v );
    }
    return w.value;
  }

  //! The packed `D` and `DD` must agree with the scalar accessors.
  real_type check_packed( Factory const & make, Grid const & g ) {
    auto const z = generic_data( g );
    auto       S = make();
    S->build( g.x, g.y, z );

    real_type s_v = 0, s_1 = 0, s_2 = 0;
    for ( auto const & [xv, yv] : cell_samples( g, 0.25, 0.75 ) ) {
      s_v = max( s_v, abs( S->eval( xv, yv ) ) );
      s_1 = max( { s_1, abs( S->Dx( xv, yv ) ), abs( S->Dy( xv, yv ) ) } );
      s_2 = max( { s_2, abs( S->Dxx( xv, yv ) ), abs( S->Dxy( xv, yv ) ), abs( S->Dyy( xv, yv ) ) } );
    }

    Worst w;
    for ( auto const & [xv, yv] : cell_samples( g, 0.25, 0.75 ) ) {
      array<real_type, 3> d{};
      S->D( xv, yv, span<real_type, 3>{ d } );
      w.add( d[0], S->eval( xv, yv ), s_v );
      w.add( d[1], S->Dx( xv, yv ), s_1 );
      w.add( d[2], S->Dy( xv, yv ), s_1 );

      array<real_type, 6> dd{};
      S->DD( xv, yv, span<real_type, 6>{ dd } );
      w.add( dd[0], S->eval( xv, yv ), s_v );
      w.add( dd[1], S->Dx( xv, yv ), s_1 );
      w.add( dd[2], S->Dy( xv, yv ), s_1 );
      w.add( dd[3], S->Dxx( xv, yv ), s_2 );
      w.add( dd[4], S->Dxy( xv, yv ), s_2 );
      w.add( dd[5], S->Dyy( xv, yv ), s_2 );
    }
    return w.value;
  }

  bool run( string const & spline_name, Factory const & make, Grid const & g ) {
    real_type const ex  = check_exact( make, g );
    auto const      fdv = check_finite_differences( make, g );
    real_type const fd  = *max_element( fdv.begin(), fdv.end() );
    real_type const cy  = check_covariance( make, g, 0.2, true );
    real_type const cx  = check_covariance( make, g, 5.0, false );
    real_type const pk  = check_packed( make, g );

    bool const ok = ex <= EXACT_TOL && fd <= 1 && cy <= COV_TOL && cx <= COV_TOL && pk <= PACKED_TOL;

    fmt::print(
      ok ? OK_COLOR : FAIL_COLOR,
      "  {:<18} exact={:8.2e}  fd={:8.2e}  cov-y={:8.2e}  cov-x={:8.2e}  packed={:8.2e}\n",
      spline_name, ex, fd, cy, cx, pk
    );
    if ( fd > 1 ) {
      fmt::print(
        FAIL_COLOR,
        "  {:<18}   fd detail: Dx={:8.2e} Dy={:8.2e} Dxx={:8.2e} Dyy={:8.2e} Dxy={:8.2e}\n",
        "", fdv[0], fdv[1], fdv[2], fdv[3], fdv[4]
      );
    }
    return ok;
  }

}  // namespace

int
main() {
  fmt::print( HEADER_COLOR, "\n=== test20: value and derivatives of surface splines ===\n\n" );
  fmt::print(
    INFO_COLOR,
    "exact  : bilinear data, all six numbers against their closed form\n"
    "fd     : derivatives against central differences of the spline itself, as a\n"
    "         margin on a mixed criterion (relative + rounding floor): <= 1 passes\n"
    "cov    : one axis rescaled, every order rescaled with it\n"
    "packed : D() and DD() against the scalar accessors\n\n"
  );

  vector<pair<string, Factory>> const splines{
    { "Bilinear", []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BilinearSpline>(); } },
    { "BiCubic[CUBIC]",
      []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BiCubicSpline>( Spline_sub_type::CUBIC ); } },
    { "BiCubic[AKIMA]",
      []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BiCubicSpline>( Spline_sub_type::AKIMA ); } },
    { "BiCubic[PCHIP]",
      []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BiCubicSpline>( Spline_sub_type::PCHIP ); } },
    { "BiQuintic[CUBIC]",
      []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BiQuinticSpline>( Spline_sub_type::CUBIC ); } },
    { "BiQuintic[AKIMA]",
      []() -> unique_ptr<SplineSurf> { return make_unique<Splines::BiQuinticSpline>( Spline_sub_type::AKIMA ); } }
  };

  bool ok = true;
  for ( auto const & g : grids() ) {
    fmt::print( INFO_COLOR, "{}\n", g.name );
    for ( auto const & [name, make] : splines ) ok = run( name, make, g ) && ok;
    fmt::print( "\n" );
  }

  if ( ok ) {
    fmt::print( OK_COLOR, "All surface derivative tests passed.\n" );
    return 0;
  }

  fmt::print(
    FAIL_COLOR,
    "Surface derivative regression detected.\n"
    "`exact` alone points at units or at the x/y roles, `fd` at the arithmetic of one\n"
    "derivative, `cov` at a derivative built on the wrong abscissas, `packed` at the\n"
    "order in which D()/DD() fill their span.\n"
  );
  return 1;
}
