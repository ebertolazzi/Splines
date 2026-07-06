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

namespace Splines
{

  void Akima_build(
    real_type const X[],  // pointers to data
    real_type const Y[],
    real_type       Yp[],
    real_type       m[],
    int             N )
  {
    assert( N >= 2 && "Akima_build requires at least 2 points" );

    if ( N == 2 )
    {
      Yp[0] = Yp[1] = ( Y[1] - Y[0] ) / ( X[1] - X[0] );
      return;
    }

    const int n = N - 1;

    for ( integer i = 0; i < n; ++i ) m[i] = ( Y[i + 1] - Y[i] ) / ( X[i + 1] - X[i] );

    auto slope = [&]( int i ) -> real_type
    {
      if ( i == -2 ) return 2 * m[0] - m[1];
      if ( i == -1 ) return 2 * ( 2 * m[0] - m[1] ) - m[0];
      if ( i >= 0 && i <= n - 1 ) return m[i];
      if ( i == n ) return 2 * m[n - 1] - m[n - 2];
      if ( i == n + 1 ) return 2 * ( 2 * m[n - 1] - m[n - 2] ) - m[n - 1];
      return 0.0;
    };

    real_type m_i   = slope( 0 );
    real_type m_im1 = slope( 0 - 1 );
    real_type m_im2 = slope( 0 - 2 );
    real_type m_ip1 = slope( 0 + 1 );
    for ( integer i = 0; i <= n; ++i )
    {
      real_type w1    = std::abs( m_ip1 - m_i ) + 0.5 * std::abs( m_ip1 + m_i );
      real_type w2    = std::abs( m_im1 - m_im2 ) + 0.5 * std::abs( m_im1 + m_im2 );
      real_type sum_w = w1 + w2;

      if ( sum_w > 1e-12 )
        Yp[i] = ( w1 * m_im1 + w2 * m_i ) / sum_w;
      else
        Yp[i] = 0.5 * ( m_im1 + m_i );  // stable fallback
      m_im2 = m_im1;
      m_im1 = m_i;
      m_i   = m_ip1;
      m_ip1 = slope( i + 2 );
    }
  }

  void VanLeer_build( real_type const X[], real_type const Y[], real_type Yp[], integer const npts )
  {
    UTILS_ASSERT( npts >= 2, "VanLeer_build: npts={} >= 2 required\n", npts );

    if ( npts == 2 )
    {  // only 2 points, linear interpolation
      Yp[0] = Yp[1] = ( Y[1] - Y[0] ) / ( X[1] - X[0] );
      return;
    }

    auto vanLeer = []( real_type a, real_type b, real_type c ) -> real_type
    {
      if ( ( a > 0 && b > 0 && c > 0 ) || ( a < 0 && b < 0 && c < 0 ) ) return 3 / ( 1 / a + 1 / b + 1 / c );
      return 0;
    };

    // ---- left boundary ----
    {
      real_type L  = ( Y[1] - Y[0] ) / ( X[1] - X[0] );
      real_type LL = Utils::first_derivative_3p( X[0], Y[0], X[1], Y[1], X[2], Y[2] );
      Yp[0]        = vanLeer( L, L, LL );
    }

    // ---- interior points ----
    for ( integer i = 1; i < npts - 1; ++i )
    {
      real_type L = ( Y[i] - Y[i - 1] ) / ( X[i] - X[i - 1] );
      real_type C = Utils::first_derivative_3p( X[i], Y[i], X[i - 1], Y[i - 1], X[i + 1], Y[i + 1] );
      real_type R = ( Y[i + 1] - Y[i] ) / ( X[i + 1] - X[i] );
      Yp[i]       = vanLeer( L, C, R );
    }

    // ---- right boundary ----
    {
      real_type R = ( Y[npts - 1] - Y[npts - 2] ) / ( X[npts - 1] - X[npts - 2] );
      real_type RR =
        Utils::first_derivative_3p( X[npts - 1], Y[npts - 1], X[npts - 2], Y[npts - 2], X[npts - 3], Y[npts - 3] );
      Yp[npts - 1] = vanLeer( RR, R, R );
    }
  }

  void Hyman_filter( real_type const X[], real_type const Y[], real_type Yp[], integer const npts )
  {
    auto minmod = []( real_type a, real_type b ) -> real_type
    {
      if ( a > 0 && b > 0 ) return std::min( a, b );
      if ( a < 0 && b < 0 ) return std::max( a, b );
      return 0;
    };

    integer const n = npts - 1;

    for ( integer i = 1; i < n; ++i )
    {
      real_type const hL = X[i] - X[i - 1];
      real_type const hR = X[i + 1] - X[i];

      real_type const delL = ( Y[i] - Y[i - 1] ) / hL;
      real_type const delR = ( Y[i + 1] - Y[i] ) / hR;

      if ( delL * delR <= 0 )
      {
        Yp[i] = 0;
        continue;
      }

      // Hyman bounds
      real_type const dmin = 3 * minmod( delL, delR );
      real_type const dmax = dmin;

      // clamp derivative
      Yp[i] = minmod( Yp[i], dmax );
    }
  }

  void Pchip_build( real_type const X[], real_type const Y[], real_type Yp[], integer const npts )
  {
    UTILS_ASSERT( npts >= 2, "Pchip_build: npts must be >= 2" );

    integer const n = npts - 1;

    // --- Special case: only two points → linear interpolation
    if ( npts == 2 )
    {
      real_type const h = X[1] - X[0];
      real_type const d = ( Y[1] - Y[0] ) / h;
      Yp[0]             = d;
      Yp[1]             = d;
      return;
    }

    // ------------------------------------------------------------
    // 1. Compute secant slopes
    // ------------------------------------------------------------
    std::vector<real_type> h( n ), delta( n );

    for ( integer i = 0; i < n; ++i )
    {
      h[i]     = X[i + 1] - X[i];
      delta[i] = ( Y[i + 1] - Y[i] ) / h[i];
    }

    // ------------------------------------------------------------
    // 2. Derivatives at interior points (Fritsch–Butland)
    // ------------------------------------------------------------
    for ( integer i = 1; i < n; ++i )
    {
      real_type const delL = delta[i - 1];
      real_type const delR = delta[i];

      if ( delL == 0 || delR == 0 || delL * delR < 0 ) { Yp[i] = 0; }
      else
      {
        real_type const wL = h[i - 1];
        real_type const wR = h[i];

        // weighted harmonic mean
        Yp[i] = ( wL + wR ) / ( ( wL / delL ) + ( wR / delR ) );
      }
    }

    // ------------------------------------------------------------
    // 3. Left boundary
    // ------------------------------------------------------------
    {
      real_type const h0 = h[0];
      real_type const h1 = h[1];
      real_type const d0 = delta[0];
      real_type const d1 = delta[1];

      Yp[0] = ( ( 2 * h0 + h1 ) * d0 - h0 * d1 ) / ( h0 + h1 );

      // monotonicity clamp
      if ( Yp[0] * d0 <= 0 )
        Yp[0] = 0;
      else if ( std::abs( Yp[0] ) > 3 * std::abs( d0 ) )
        Yp[0] = 3 * d0;
    }

    // ------------------------------------------------------------
    // 4. Right boundary
    // ------------------------------------------------------------
    {
      real_type const hn1 = h[n - 1];
      real_type const hn2 = h[n - 2];
      real_type const dn1 = delta[n - 1];
      real_type const dn2 = delta[n - 2];

      Yp[n] = ( ( 2 * hn1 + hn2 ) * dn1 - hn1 * dn2 ) / ( hn1 + hn2 );

      // monotonicity clamp
      if ( Yp[n] * dn1 <= 0 )
        Yp[n] = 0;
      else if ( std::abs( Yp[n] ) > 3 * std::abs( dn1 ) )
        Yp[n] = 3 * dn1;
    }

    Hyman_filter( X, Y, Yp, npts );
  }

  void Quintic_build(
    Spline_sub_type const q_sub_type,
    real_type const       X[],
    real_type const       Y[],
    real_type             Yp[],
    real_type             Ypp[],
    integer const         npts )
  {
    UTILS_ASSERT( npts >= 2, "Quintic_build, npts={} must be >= 2\n", npts );

    // ========================================================================
    // PHASE 1: First derivative estimation
    // ========================================================================
    switch ( q_sub_type )
    {
      case Spline_sub_type::CUBIC:
        CubicSpline_build( X, Y, Yp, Ypp, npts, CubicSpline_BC::EXTRAPOLATE, CubicSpline_BC::EXTRAPOLATE );
        break;
      case Spline_sub_type::PCHIP: Pchip_build( X, Y, Yp, npts ); break;
      case Spline_sub_type::AKIMA:
      {
        Vec work( npts );
        Akima_build( X, Y, Yp, work.data(), npts );
      }
      break;
      case Spline_sub_type::VANLEER: VanLeer_build( X, Y, Yp, npts ); break;
      default: UTILS_ERROR( "Unknown QuinticSpline_sub_type value\n" );
    }

    // ========================================================================
    // TRIVIAL CASE: 2 points
    // ========================================================================
    if ( npts == 2 )
    {
      Ypp[0] = Ypp[1] = 0;
      return;
    }

    integer const n = npts - 1;

    // ========================================================================
    // OPTIMIZATION 2: Single allocation for all interval data
    // ========================================================================
    Malloc_real mem_intervals( "QuinticSpline::intervals" );
    real_type * interval_data = mem_intervals.malloc( 4 * n );  // h, h2, h3, inv_h

    real_type * h     = interval_data;
    real_type * h2    = interval_data + n;
    real_type * h3    = interval_data + 2 * n;
    real_type * inv_h = interval_data + 3 * n;

    // Compute all interval terms in one pass (improves cache locality)
    for ( integer i = 0; i < n; ++i )
    {
      h[i] = X[i + 1] - X[i];
      UTILS_ASSERT( h[i] > 0, "X must be strictly increasing" );

      inv_h[i] = 1.0 / h[i];
      h2[i]    = h[i] * h[i];
      h3[i]    = h2[i] * h[i];
    }

    // ========================================================================
    // OPTIMIZATION 3: Tridiagonal system with optimized allocation
    // ========================================================================
    Utils::TridiagonalSolver<real_type> solver;
    solver.resize( npts );

    // Contiguous allocation for system vectors (improves cache performance)
    Malloc_real mem_system( "QuinticSpline::system" );
    real_type * system_data = mem_system.malloc( 3 * npts - 2 );

    // Direct mapping to allocated arrays (zero overhead)
    Eigen::Map<Vec> a{ system_data, npts - 1 };
    Eigen::Map<Vec> b{ system_data + npts - 1, npts };
    Eigen::Map<Vec> c{ system_data + 2 * npts - 1, npts - 1 };

    // d maps directly to Ypp (zero-copy!)
    Eigen::Map<Vec> d{ Ypp, npts };

    // ========================================================================
    // LEFT BOUNDARY (i=0)
    // ========================================================================
    {
      real_type const h0_inv = inv_h[0];
      real_type const dY     = Y[1] - Y[0];

      d( 0 ) = -10.0 * dY * h0_inv + 6.0 * Yp[0] + 4.0 * Yp[1];
      b( 0 ) = 2.0 * h[0];
      c( 0 ) = h[0];
    }

    // ========================================================================
    // INTERNAL EQUATIONS (i=1..n-1)
    // OPTIMIZATION 4: Reduced operations and variable reuse
    // ========================================================================
    for ( integer i = 1; i < n; ++i )
    {
      integer const iL = i - 1;

      real_type const hL_inv  = inv_h[iL];
      real_type const hR_inv  = inv_h[i];
      real_type const hL2_inv = hL_inv * hL_inv;
      real_type const hR2_inv = hR_inv * hR_inv;
      real_type const hL3_inv = hL2_inv * hL_inv;
      real_type const hR3_inv = hR2_inv * hR_inv;

      // Matrix coefficients
      a( iL ) = -hL_inv;
      b( i )  = 3.0 * ( hL_inv + hR_inv );
      c( i )  = -hR_inv;

      // RHS: factorize common calculations
      real_type const dY_L = Y[i] - Y[iL];
      real_type const dY_R = Y[i + 1] - Y[i];

      real_type const termY = 20.0 * ( dY_R * hR3_inv - dY_L * hL3_inv );

      // Compute derivative terms with fewer operations
      real_type const termYp_R = ( 3.0 * Yp[i] + 2.0 * Yp[i + 1] ) * hR2_inv;
      real_type const termYp_L = ( 2.0 * Yp[iL] + 3.0 * Yp[i] ) * hL2_inv;

      d( i ) = termY - 4.0 * ( termYp_R - termYp_L );
    }

    // ========================================================================
    // RIGHT BOUNDARY (i=n)
    // ========================================================================
    {
      real_type const hn_inv = inv_h[n - 1];
      real_type const dY     = Y[n] - Y[n - 1];

      a( n - 1 ) = h[n - 1];
      b( n )     = 2.0 * h[n - 1];
      d( n )     = 10.0 * dY * hn_inv - 4.0 * Yp[n - 1] - 6.0 * Yp[n];
    }

    // ========================================================================
    // SOLUTION: solve_inplace modifies d which is already mapped to Ypp
    // ========================================================================
    solver.factorize( a, b, c );
    solver.solve_inplace( a, b, d );
  }

  void WENO5_build( real_type const X[], real_type const Y[], real_type Yp[], integer const npts )
  {
    UTILS_ASSERT( npts >= 2, "WENO_build requires at least 2 points" );

    constexpr real_type eps = 1e-12;

    // --- Lambda: WENO5 derivative su 5 nodi non uniformi ---
    auto weno5_derivative = [&]( real_type x[5], real_type y[5] ) -> real_type
    {
      real_type d0 = ( -3 * y[2] + 4 * y[3] - y[4] ) / ( x[4] - x[2] );
      real_type d1 = ( y[3] - y[1] ) / ( x[3] - x[1] );
      real_type d2 = ( 3 * y[2] - 4 * y[1] + y[0] ) / ( x[2] - x[0] );

      auto beta = [&]( real_type a, real_type b, real_type c, real_type dx ) -> real_type
      {
        real_type val = a - 2 * b + c;
        return val * val / ( dx * dx );
      };

      real_type beta0 = beta( y[2], y[3], y[4], x[4] - x[2] );
      real_type beta1 = beta( y[1], y[2], y[3], x[3] - x[1] );
      real_type beta2 = beta( y[0], y[1], y[2], x[2] - x[0] );

      constexpr real_type g0 = 0.1, g1 = 0.6, g2 = 0.3;
      real_type           a0  = g0 / ( ( eps + beta0 ) * ( eps + beta0 ) );
      real_type           a1  = g1 / ( ( eps + beta1 ) * ( eps + beta1 ) );
      real_type           a2  = g2 / ( ( eps + beta2 ) * ( eps + beta2 ) );
      real_type           sum = a0 + a1 + a2;
      a0 /= sum;
      a1 /= sum;
      a2 /= sum;

      return a0 * d0 + a1 * d1 + a2 * d2;
    };

    // --- Lambda: fallback PCHIP/lineare ---
    auto fallback = [&]( integer i ) -> real_type
    {
      if ( npts == 2 ) return ( Y[1] - Y[0] ) / ( X[1] - X[0] );
      if ( npts == 3 )
      {
        if ( i == 0 ) return ( Y[1] - Y[0] ) / ( X[1] - X[0] );
        if ( i == 1 ) return ( Y[2] - Y[0] ) / ( X[2] - X[0] );
        return ( Y[2] - Y[1] ) / ( X[2] - X[1] );
      }
      if ( npts == 4 )
      {
        if ( i == 0 ) return ( Y[1] - Y[0] ) / ( X[1] - X[0] );
        if ( i == 1 ) return ( Y[2] - Y[0] ) / ( X[2] - X[0] );
        if ( i == 2 ) return ( Y[3] - Y[1] ) / ( X[3] - X[1] );
        return ( Y[3] - Y[2] ) / ( X[3] - X[2] );
      }
      // bordi
      integer im1 = std::max<integer>( 0, i - 1 );
      integer ip1 = std::min<integer>( npts - 1, i + 1 );
      return ( Y[ip1] - Y[im1] ) / ( X[ip1] - X[im1] );
    };

    // --- Lambda: Hyman clamp per monotonia ---
    auto hyman_clamp = [&]( integer i, real_type d ) -> real_type
    {
      integer   im1  = std::max<integer>( 0, i - 1 );
      integer   ip1  = std::min<integer>( npts - 1, i + 1 );
      real_type dL   = ( Y[i] - Y[im1] ) / ( X[i] - X[im1] );
      real_type dR   = ( Y[ip1] - Y[i] ) / ( X[ip1] - X[i] );
      real_type dmin = std::min( dL, dR );
      real_type dmax = std::max( dL, dR );
      return std::min( std::max( d, dmin ), dmax );
    };

    // --- main loop ---
    for ( integer i = 0; i < npts; ++i )
    {
      real_type d;
      if ( npts < 5 || i < 2 || i > npts - 3 ) { d = fallback( i ); }
      else
      {
        real_type x5[5] = { X[i - 2], X[i - 1], X[i], X[i + 1], X[i + 2] };
        real_type y5[5] = { Y[i - 2], Y[i - 1], Y[i], Y[i + 1], Y[i + 2] };
        d               = weno5_derivative( x5, y5 );
      }
      Yp[i] = hyman_clamp( i, d );
    }
  }

  void WENO5_build2( real_type const X[], real_type const Y[], real_type Ypp[], integer const npts )
  {
    UTILS_ASSERT( npts >= 2, "Need at least 2 points" );

    constexpr real_type eps = 1e-12;

    auto second_derivative_candidate = []( real_type x[3], real_type y[3] ) -> real_type
    {
      return 2.0 * ( y[0] / ( ( x[0] - x[1] ) * ( x[0] - x[2] ) ) + y[1] / ( ( x[1] - x[0] ) * ( x[1] - x[2] ) ) +
                     y[2] / ( ( x[2] - x[0] ) * ( x[2] - x[1] ) ) );
    };

    auto beta_func = []( real_type d, real_type dx ) -> real_type { return d * d * dx * dx; };

    auto fallback = [&]( integer i ) -> real_type
    {
      if ( npts == 2 ) return 0;  // solo due punti: seconda derivata ignota
      if ( npts == 3 )
      {
        integer im1 = std::max<integer>( 0, i - 1 );
        integer ip1 = std::min<integer>( npts - 1, i + 1 );
        return 2.0 * ( Y[im1] - 2 * Y[i] + Y[ip1] ) / ( ( X[ip1] - X[im1] ) * ( X[ip1] - X[im1] ) );
      }
      if ( npts == 4 )
      {
        integer im1 = std::max<integer>( 0, i - 1 );
        integer ip1 = std::min<integer>( npts - 1, i + 1 );
        return 2.0 * ( Y[im1] - 2 * Y[i] + Y[ip1] ) / ( ( X[ip1] - X[im1] ) * ( X[ip1] - X[im1] ) );
      }
      // bordi
      integer im1 = std::max<integer>( 0, i - 1 );
      integer ip1 = std::min<integer>( npts - 1, i + 1 );
      return 2.0 * ( Y[im1] - 2 * Y[i] + Y[ip1] ) / ( ( X[ip1] - X[im1] ) * ( X[ip1] - X[im1] ) );
    };

    // Hyman clamp per seconda derivata (concavità)
    auto hyman_clamp = [&]( integer i, real_type d ) -> real_type
    {
      integer   im1 = std::max<integer>( 0, i - 1 );
      integer   ip1 = std::min<integer>( npts - 1, i + 1 );
      real_type dL  = ( Y[i] - Y[im1] ) / ( X[i] - X[im1] );
      real_type dR  = ( Y[ip1] - Y[i] ) / ( X[ip1] - X[i] );
      // preserva concavità tra punti
      if ( dL * dR >= 0 )
      {
        real_type sign = ( dL > 0 ? 1.0 : -1.0 );
        real_type dmax = 2.0 * std::max( std::abs( dL ), std::abs( dR ) );
        return std::min( std::abs( d ), dmax ) * sign;
      }
      return d;
    };

    for ( integer i = 0; i < npts; ++i )
    {
      real_type d;
      if ( npts < 5 || i < 2 || i > npts - 3 ) { d = fallback( i ); }
      else
      {
        real_type x0[3] = { X[i - 2], X[i - 1], X[i] };
        real_type x1[3] = { X[i - 1], X[i], X[i + 1] };
        real_type x2[3] = { X[i], X[i + 1], X[i + 2] };
        real_type y0[3] = { Y[i - 2], Y[i - 1], Y[i] };
        real_type y1[3] = { Y[i - 1], Y[i], Y[i + 1] };
        real_type y2[3] = { Y[i], Y[i + 1], Y[i + 2] };

        real_type d0 = second_derivative_candidate( x0, y0 );
        real_type d1 = second_derivative_candidate( x1, y1 );
        real_type d2 = second_derivative_candidate( x2, y2 );

        real_type beta0 = beta_func( d0, X[i] - X[i - 2] );
        real_type beta1 = beta_func( d1, X[i + 1] - X[i - 1] );
        real_type beta2 = beta_func( d2, X[i + 2] - X[i] );

        constexpr real_type g0 = 0.1, g1 = 0.6, g2 = 0.3;
        real_type           a0  = g0 / ( ( eps + beta0 ) * ( eps + beta0 ) );
        real_type           a1  = g1 / ( ( eps + beta1 ) * ( eps + beta1 ) );
        real_type           a2  = g2 / ( ( eps + beta2 ) * ( eps + beta2 ) );
        real_type           sum = a0 + a1 + a2;
        a0 /= sum;
        a1 /= sum;
        a2 /= sum;

        d = a0 * d0 + a1 * d1 + a2 * d2;
      }
      Ypp[i] = hyman_clamp( i, d );
    }
  }

  void Quintic_build2(
    Spline_sub_type const q_sub_type,
    real_type const       X[],
    real_type const       Y[],
    real_type             Yp[],
    real_type             Ypp[],
    integer const         npts )
  {
    UTILS_ASSERT( npts >= 2, "Quintic_build2, npts={} must be >= 2\n", npts );

    // ========================================================================
    // PHASE 1: First derivative estimation
    // ========================================================================
    switch ( q_sub_type )
    {
      case Spline_sub_type::CUBIC:
        CubicSpline_build( X, Y, Yp, Ypp, npts, CubicSpline_BC::EXTRAPOLATE, CubicSpline_BC::EXTRAPOLATE );
        break;
      case Spline_sub_type::PCHIP: Pchip_build( X, Y, Yp, npts ); break;
      case Spline_sub_type::AKIMA:
      {
        Vec work( npts );
        Akima_build( X, Y, Yp, work.data(), npts );
      }
      break;
      case Spline_sub_type::VANLEER: VanLeer_build( X, Y, Yp, npts ); break;
      default: UTILS_ERROR( "Unknown QuinticSpline_sub_type value\n" );
    }
    // WENO5_build2( X, Y, Ypp, npts );
    WENO5_build( X, Yp, Ypp, npts );
  }


}  // namespace Splines

//
// EOF: SplineBuild.cc
//
