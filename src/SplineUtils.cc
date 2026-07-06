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

  integer check_cubic_spline_monotonicity(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    integer const   npts )
  {
    integer flag = 1;  // Assume strictly monotone initially

    // Check data monotonicity (X is assumed monotonically increasing)
    for ( integer i = 1; i < npts; ++i )
    {
      if ( Y[i] < Y[i - 1] ) return -2;  // Non-monotone data detected

      // Check for flat regions (non-strict monotonicity)
      if ( flag > 0 && Utils::is_zero( Y[i] - Y[i - 1] ) && X[i] > X[i - 1] )
      {
        flag = 0;  // Monotone but not strictly
      }
    }

    // Check spline monotonicity using Fritsch-Carlson conditions
    // The conditions ensure that the cubic interpolant preserves monotonicity
    for ( integer i = 1; i < npts; ++i )
    {
      const real_type dx = X[i] - X[i - 1];
      if ( dx <= real_type( 0 ) ) continue;  // Skip duplicate or invalid x-values

      const real_type dy = Y[i] - Y[i - 1];
      const real_type dd = dy / dx;  // Average slope of the interval

      // Normalize derivatives by the average slope
      const real_type m0 = Yp[i - 1] / dd;
      const real_type m1 = Yp[i] / dd;

      // Condition 1: Both normalized derivatives must be non-negative
      if ( m0 < 0 || m1 < 0 ) return -1;

      // Condition 2: For monotonicity, m0 and m1 should satisfy specific bounds
      if ( m0 <= 3 && m1 <= 3 )
      {
        // Simple case: both normalized derivatives are in [0, 3]
        // This is sufficient for monotonicity

        if ( flag > 0 )
        {
          // Check for boundary conditions that indicate non-strict monotonicity
          if (
            ( i > 1 && ( Utils::is_zero( m0 ) || Utils::is_zero( m0 - 3 ) ) ) ||
            ( i < npts - 1 && ( Utils::is_zero( m1 ) || Utils::is_zero( m1 - 3 ) ) ) )
          {
            flag = 0;  // Non-strict monotonicity
          }
        }
      }
      else
      {
        // Complex case: at least one normalized derivative exceeds 3
        // Apply the general Fritsch-Carlson discriminant test

        const real_type t1   = 2 * m0 + m1 - 3;
        const real_type t2   = 2 * ( m0 + m1 - 2 );
        const real_type disc = m0 * t2 - t1 * t1;

        // The spline is monotone if disc and t2 have the same sign or disc = 0
        if ( ( t2 >= 0 && disc < 0 ) || ( t2 < 0 && disc > 0 ) )
        {
          return -1;  // Monotonicity violation
        }

        // Check for non-strict monotonicity
        if ( flag > 0 && Utils::is_zero( disc ) ) { flag = 0; }
      }
    }

    return flag;
  }

  integer check_quintic_spline_monotonicity(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    real_type const Ypp[],
    integer const   npts )
  {
    /**
     * @brief Helper lambda to check monotonicity of a single quintic interval
     *
     * @param[in] m0  Normalized first derivative at interval start
     * @param[in] m1  Normalized first derivative at interval end
     * @param[in] c0  Normalized second derivative at interval start
     * @param[in] c1  Normalized second derivative at interval end
     *
     * @return true if the interval is monotone, false otherwise
     *
     * @details For quintic Hermite interpolation, the derivative is a quartic
     *          polynomial. This function applies sufficient conditions to ensure
     *          the derivative remains non-negative throughout [0,1].
     */
    auto check_quintic_interval_monotonicity =
      []( real_type const m0, real_type const m1, real_type const c0, real_type const c1 ) -> bool
    {
      // Condition 1: First derivatives must be non-negative
      if ( m0 < 0 || m1 < 0 ) return false;

      // Condition 2: Second derivatives must be bounded
      // Excessively large second derivatives can cause the derivative to become negative
      constexpr real_type UPPER_BOUND_C = 10.0;
      if ( std::abs( c0 ) > UPPER_BOUND_C || std::abs( c1 ) > UPPER_BOUND_C ) return false;

      // Condition 3: When first derivatives are large, constrain second derivatives
      // Large positive first derivatives with large negative second derivatives
      // can create a dip that violates monotonicity
      if ( m0 > 3 || m1 > 3 )
      {
        constexpr real_type MAX_NEGATIVE_C = -2.0;
        if ( c0 < MAX_NEGATIVE_C || c1 < MAX_NEGATIVE_C ) return false;
      }

      // Condition 4: Check derivative at midpoint as a heuristic
      // The quintic basis functions for the derivative at t=0.5:
      constexpr real_type t            = 0.5;
      constexpr real_type t2           = t * t;
      constexpr real_type t3           = t * t2;
      constexpr real_type one_minus_t  = 1.0 - t;
      constexpr real_type one_minus_t2 = one_minus_t * one_minus_t;
      constexpr real_type one_minus_t3 = one_minus_t * one_minus_t2;

      // Evaluate normalized derivative at t = 0.5
      // This uses the derivative of the quintic Hermite basis functions
      const real_type deriv_at_midpoint = m0 * one_minus_t3 * ( 1 + 3 * t + 6 * t2 ) +
                                          m1 * t3 * ( 10 - 15 * t + 6 * t2 ) +
                                          c0 * 3 * t * one_minus_t2 * ( 1 + 2 * t ) +
                                          c1 * 3 * t2 * ( t - 1 ) * ( 4 - 3 * t );

      constexpr real_type TOLERANCE = 1e-6;
      if ( deriv_at_midpoint < -TOLERANCE ) return false;

      return true;
    };

    integer flag = 1;  // Assume strictly monotone initially

    // Check data monotonicity (X is assumed monotonically increasing)
    for ( integer i = 1; i < npts; ++i )
    {
      if ( Y[i] < Y[i - 1] ) return -2;  // Non-monotone data

      // Check for flat regions
      if ( flag > 0 && Utils::is_zero( Y[i] - Y[i - 1] ) && X[i] > X[i - 1] )
      {
        flag = 0;  // Non-strict monotonicity
      }
    }

    // Check spline monotonicity interval by interval
    for ( integer i = 1; i < npts; ++i )
    {
      const real_type dx = X[i] - X[i - 1];
      if ( dx <= real_type( 0 ) ) continue;  // Skip invalid intervals

      const real_type dy = Y[i] - Y[i - 1];

      // Handle constant segments specially
      if ( Utils::is_zero( dy ) )
      {
        // For flat segments, all derivatives should be zero for monotonicity
        if (
          !Utils::is_zero( Yp[i - 1] ) || !Utils::is_zero( Yp[i] ) || !Utils::is_zero( Ypp[i - 1] ) ||
          !Utils::is_zero( Ypp[i] ) )
        {
          return -1;  // Non-zero derivatives on flat segment violate monotonicity
        }
        continue;
      }

      const real_type dd = dy / dx;  // Average slope

      // Normalize derivatives by average slope
      const real_type m0 = Yp[i - 1] / dd;
      const real_type m1 = Yp[i] / dd;
      const real_type c0 = Ypp[i - 1] * dx / ( dd * 2 );
      const real_type c1 = Ypp[i] * dx / ( dd * 2 );

      // Basic positivity check
      if ( m0 < 0 || m1 < 0 ) return -1;

      // Apply monotonicity conditions based on derivative magnitudes
      if ( m0 <= 3 && m1 <= 3 )
      {
        // Simple case: normalized first derivatives in [0, 3]

        // Check for strict monotonicity
        if ( flag > 0 )
        {
          if (
            ( i > 1 && ( Utils::is_zero( m0 ) || Utils::is_zero( m0 - 3 ) ) ) ||
            ( i < npts - 1 && ( Utils::is_zero( m1 ) || Utils::is_zero( m1 - 3 ) ) ) )
          {
            flag = 0;
          }
        }

        // For quintic splines, second derivatives must also be bounded
        constexpr real_type MAX_C = 4.0;
        if ( std::abs( c0 ) > MAX_C || std::abs( c1 ) > MAX_C )
        {
          // Need more thorough check when second derivatives are large
          if ( !check_quintic_interval_monotonicity( m0, m1, c0, c1 ) ) { return -1; }
        }
      }
      else
      {
        // Complex case: at least one normalized derivative exceeds 3
        // Use comprehensive interval check
        if ( !check_quintic_interval_monotonicity( m0, m1, c0, c1 ) ) { return -1; }

        // Check for strict monotonicity
        if ( flag > 0 )
        {
          if ( Utils::is_zero( m0 - 3 ) || Utils::is_zero( m1 - 3 ) || Utils::is_zero( m0 ) || Utils::is_zero( m1 ) )
          {
            flag = 0;
          }
        }
      }
    }

    return flag;
  }

}  // namespace Splines

//
// EOF: SplineUtils.cc
//
