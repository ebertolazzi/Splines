/*\
 |  Splines micro-benchmark — scalar vs. batch evaluation hot paths.
 |
 |  Baseline + validation harness for the "faster" work in MODERNIZATION.md
 |  (Tier 3). Times the per-point scalar path (operator(), D, DD) against the
 |  batch overloads eval/D/DD(span,span), for Linear/Cubic/Akima/Pchip/Quintic
 |  splines, over both RANDOM and SORTED query arrays (the batch's interval-
 |  reuse fast path only helps monotone/clustered inputs — the common
 |  resample/plot case). Best-of-N with a checksum sink to defeat DCE.
 |
 |  Also VERIFIES the batch results are bit-identical to the scalar path on
 |  adversarial inputs (exact knots, just-off knots, out-of-domain, sorted,
 |  and a closed spline that exercises periodic wrapping) — the batch API's
 |  hard correctness contract. Exits non-zero on any mismatch.
 |
 |  Usage: bench_eval [N_knots] [M_queries] [reps]
 |         defaults: 1000 knots, 2,000,000 queries, best of 5
\*/

#include "Splines/Splines.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <span>
#include <vector>

using Splines::integer;
using Splines::real_type;

namespace
{
  using clock_type = std::chrono::steady_clock;

  template <typename F>
  double best_ns_per_op( F const & fn, std::size_t n_ops, int reps, double & sink )
  {
    double best = std::numeric_limits<double>::infinity();
    for ( int r = 0; r < reps; ++r )
    {
      auto const   t0 = clock_type::now();
      sink += fn();
      auto const   t1 = clock_type::now();
      double const ns = std::chrono::duration<double, std::nano>( t1 - t0 ).count();
      best            = std::min( best, ns / static_cast<double>( n_ops ) );
    }
    return best;
  }

  // Batch must be bit-identical to the scalar path, for every op, on the
  // given queries. Returns false (and reports) on the first divergence.
  bool verify(
    char const *                   name,
    char const *                   qset,
    Splines::Spline const &        s,
    std::vector<real_type> const & xq )
  {
    std::size_t const          M = xq.size();
    std::span<real_type const> xs{ xq };
    std::vector<real_type>     ys( M ), yb( M );
    bool                       ok = true;

    auto one = [&]( char const * op, auto scalar_fn, auto batch_fn ) {
      for ( std::size_t i = 0; i < M; ++i ) ys[i] = scalar_fn( xq[i] );
      batch_fn( xs, std::span<real_type>{ yb } );
      for ( std::size_t i = 0; i < M; ++i )
        if ( ys[i] != yb[i] )
        {
          std::printf(
            "  MISMATCH %s[%s].%s at i=%zu (x=%.17g): scalar=%.17g batch=%.17g\n",
            name, qset, op, i, xq[i], ys[i], yb[i] );
          ok = false;
          return;
        }
    };

    one( "eval", [&]( real_type x ) { return s( x ); },     [&]( auto a, auto b ) { s.eval( a, b ); } );
    one( "D",    [&]( real_type x ) { return s.D( x ); },   [&]( auto a, auto b ) { s.D( a, b ); } );
    one( "DD",   [&]( real_type x ) { return s.DD( x ); },  [&]( auto a, auto b ) { s.DD( a, b ); } );
    return ok;
  }

  void bench_row(
    char const *                   name,
    Splines::Spline &              s,
    std::vector<real_type> const & xq,
    int                            reps,
    double &                       sink )
  {
    std::size_t const          M = xq.size();
    std::span<real_type const> xs{ xq };
    std::vector<real_type>     yb( M );

    double const s_eval = best_ns_per_op(
      [&] { double a = 0; for ( real_type const x : xq ) a += s( x ); return a; }, M, reps, sink );
    double const b_eval = best_ns_per_op(
      [&] { s.eval( xs, std::span<real_type>{ yb } ); return yb[M / 2]; }, M, reps, sink );
    double const s_DD = best_ns_per_op(
      [&] { double a = 0; for ( real_type const x : xq ) a += s.DD( x ); return a; }, M, reps, sink );
    double const b_DD = best_ns_per_op(
      [&] { s.DD( xs, std::span<real_type>{ yb } ); return yb[M / 2]; }, M, reps, sink );

    std::printf(
      "%-9s %9.2f %9.2f %7.2fx %9.2f %9.2f %7.2fx\n",
      name, s_eval, b_eval, s_eval / b_eval, s_DD, b_DD, s_DD / b_DD );
  }
}  // namespace

int main( int argc, char ** argv )
{
  integer     N    = 1000;
  std::size_t M    = 2'000'000;
  int         reps = 5;
  if ( argc > 1 ) N    = std::atoi( argv[1] );
  if ( argc > 2 ) M    = static_cast<std::size_t>( std::atoll( argv[2] ) );
  if ( argc > 3 ) reps = std::atoi( argv[3] );

  std::vector<real_type>                    X( static_cast<std::size_t>( N ) );
  std::vector<real_type>                    Y( static_cast<std::size_t>( N ) );
  std::mt19937_64                           rng( 12345 );
  std::uniform_real_distribution<real_type> gap( 0.2, 1.0 );
  real_type                                 acc = 0;
  for ( integer i = 0; i < N; ++i )
  {
    X[static_cast<std::size_t>( i )] = acc;
    Y[static_cast<std::size_t>( i )] = std::sin( 0.3 * acc ) + 0.1 * acc;
    acc += gap( rng );
  }

  // Query sets. `random`: uniform in-range. `sorted`: the same values sorted
  // (interval-reuse target). `adversarial`: exact knots, just-off knots,
  // midpoints, and out-of-domain points on both sides -- stresses the reuse
  // guard's boundary behaviour.
  std::vector<real_type>                    xrand( M );
  std::uniform_real_distribution<real_type> q( X.front(), X.back() );
  for ( real_type & v : xrand ) v = q( rng );
  std::vector<real_type> xsort( xrand );
  std::sort( xsort.begin(), xsort.end() );

  std::vector<real_type> xadv;
  real_type const        span = X.back() - X.front();
  for ( std::size_t k = 0; k + 1 < X.size(); ++k )
  {
    real_type const xk = X[k];
    xadv.push_back( xk );                                    // exact knot
    xadv.push_back( std::nextafter( xk, X.back() ) );        // just above
    xadv.push_back( std::nextafter( xk, X.front() ) );       // just below
    xadv.push_back( 0.5 * ( X[k] + X[k + 1] ) );             // midpoint
  }
  xadv.push_back( X.back() );
  xadv.push_back( X.front() - 0.3 * span );                  // out-of-domain left
  xadv.push_back( X.back() + 0.3 * span );                   // out-of-domain right
  xadv.push_back( X.front() - 2.7 * span );                  // far left (wrap test)
  xadv.push_back( X.back() + 2.7 * span );                   // far right (wrap test)

  // ---- build the splines ------------------------------------------------
  Splines::LinearSpline  li;  li.build( X, Y );
  Splines::CubicSpline   cs;  cs.build( X, Y );
  Splines::AkimaSpline   ak;  ak.build( X, Y );
  Splines::PchipSpline   pc;  pc.build( X, Y );
  Splines::QuinticSpline qs( Splines::Spline_sub_type::CUBIC ); qs.build( X, Y );
  Splines::CubicSpline   csc; csc.build( X, Y ); csc.make_closed();  // periodic

  struct Entry { char const * name; Splines::Spline * s; };
  Entry const all[] = {
    { "Linear", &li }, { "Cubic", &cs }, { "Akima", &ak },
    { "Pchip", &pc }, { "Quintic", &qs }, { "Cubic(closed)", &csc } };

  // ---- correctness: batch == scalar, bit-for-bit ------------------------
  bool ok = true;
  for ( Entry const & e : all )
  {
    ok = verify( e.name, "rand", *e.s, xrand ) && ok;
    ok = verify( e.name, "sort", *e.s, xsort ) && ok;
    ok = verify( e.name, "adv",  *e.s, xadv )  && ok;
  }
  std::printf( "batch==scalar (rand+sort+adversarial, incl. closed): %s\n\n", ok ? "OK" : "FAILED" );
  if ( !ok ) return 1;

  // ---- timing ------------------------------------------------------------
  double sink = 0;
  std::printf( "N=%d knots, M=%zu queries, best of %d, ns/op\n", N, M, reps );

  std::printf( "-- RANDOM queries (no interval locality) --\n" );
  std::printf( "%-9s %9s %9s %8s %9s %9s %8s\n", "type", "eval_s", "eval_b", "spd", "DD_s", "DD_b", "spd" );
  for ( Entry const & e : all ) if ( e.s != &csc ) bench_row( e.name, *e.s, xrand, reps, sink );

  std::printf( "-- SORTED queries (interval-reuse fast path) --\n" );
  std::printf( "%-9s %9s %9s %8s %9s %9s %8s\n", "type", "eval_s", "eval_b", "spd", "DD_s", "DD_b", "spd" );
  for ( Entry const & e : all ) if ( e.s != &csc ) bench_row( e.name, *e.s, xsort, reps, sink );

  std::printf( "\nchecksum: %.6g\n", sink );
  return 0;
}
