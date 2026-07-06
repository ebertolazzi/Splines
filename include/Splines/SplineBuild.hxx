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

#pragma once

#ifndef SPLINE_BUILD_HXX
#define SPLINE_BUILD_HXX

namespace Splines
{

  /**
   * @brief Build Akima spline derivative estimates.
   *
   * This function computes the first derivatives for Akima spline interpolation
   * using local slope estimates with tension control to reduce overshoot.
   *
   * Algorithm:
   * - Compute slopes between consecutive points
   * - For each point, consider 4 adjacent slopes (2 left, 2 right)
   * - Apply weighted harmonic mean of adjacent slopes where weights are
   *   based on the difference between slopes
   * - Provides C1 continuity with local control
   *
   * Reference:
   * - H. Akima, "A new method of interpolation and smooth curve fitting
   *   based on local procedures", Journal of the ACM, 1970.
   *
   * @param[in]  X  Array of strictly increasing x-coordinates
   * @param[in]  Y  Array of y-coordinates
   * @param[out] Yp Array to store computed first derivatives
   * @param[out] m  Temporary workspace (size N-1 for slopes)
   * @param[in]  N  Number of data points (must be >= 2)
   */
  void Akima_build(
    real_type const X[],  // pointers to data
    real_type const Y[],
    real_type       Yp[],
    real_type       m[],
    int             N );

  /**
   * @brief Van Leer cubic derivative reconstruction (true order 3).
   *
   * First derivatives computed as derivatives of the local cubic
   * interpolating polynomial using 4-point stencils.
   *
   * Algorithm:
   * - For interior points: 5-point central difference (order 4)
   * - For boundary points: 4-point forward/backward difference
   * - For very small datasets (npts <= 3): fallback to 3-point formulas
   *
   * Requirements:
   * - X strictly monotonic
   * - npts >= 4 for full accuracy
   *
   * Reference:
   * - F. B. Hildebrand, "Introduction to Numerical Analysis", McGraw-Hill, 1974.
   * - M. K. Horn, "Numerical Methods for Scientists and Engineers", 1990.
   *
   * @param[in]  X     Array of strictly increasing x-coordinates
   * @param[in]  Y     Array of y-coordinates
   * @param[out] Yp    Array to store computed first derivatives
   * @param[in]  npts  Number of data points (must be >= 2)
   */
  void VanLeer_build( real_type const X[], real_type const Y[], real_type Yp[], integer const npts );

  /**
   * @brief Hyman monotonicity filter for derivative estimates.
   *
   * Applies minmod limiter to enforce monotonicity preservation.
   * If left and right slopes have opposite signs, sets derivative to zero.
   * Otherwise, clamps derivative to 3× the minmod of adjacent slopes.
   *
   * Algorithm:
   * - For each interior point, check if adjacent slopes have same sign
   * - If signs differ: set derivative to 0 (local extremum)
   * - Otherwise: apply minmod limiter with factor 3
   *
   * Reference:
   * - J. M. Hyman, "Accurate monotonicity preserving cubic interpolation",
   *   SIAM Journal on Scientific and Statistical Computing, 1983.
   *
   * @param[in]     X     Array of x-coordinates
   * @param[in]     Y     Array of y-coordinates
   * @param[in,out] Yp    Array of first derivatives (filtered in-place)
   * @param[in]     npts  Number of data points
   */
  void Hyman_filter( real_type const X[], real_type const Y[], real_type Yp[], integer const npts );

  /**
   * @brief Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) derivative construction.
   *
   * Computes first derivatives that preserve monotonicity of the data.
   * Uses weighted harmonic mean of adjacent slopes with special boundary treatments.
   *
   * Algorithm:
   * 1. Compute secant slopes between consecutive points
   * 2. For interior points: weighted harmonic mean (Fritsch–Butland formula)
   * 3. For boundaries: quadratic extrapolation with monotonicity clamping
   * 4. Apply Hyman filter for additional monotonicity enforcement
   *
   * Properties:
   * - C1 continuous
   * - Shape preserving (monotonicity)
   * - Local support
   *
   * Reference:
   * - F. N. Fritsch and J. Butland, "A method for constructing local
   *   monotone piecewise cubic interpolants", SIAM Journal on Scientific
   *   and Statistical Computing, 1984.
   * - F. N. Fritsch and R. E. Carlson, "Monotone piecewise cubic interpolation",
   *   SIAM Journal on Numerical Analysis, 1980.
   *
   * @param[in]  X     Array of strictly increasing x-coordinates
   * @param[in]  Y     Array of y-coordinates
   * @param[out] Yp    Array to store computed first derivatives
   * @param[in]  npts  Number of data points (must be >= 2)
   */
  void Pchip_build( real_type const X[], real_type const Y[], real_type Yp[], integer const npts );

  /**
   * @brief Construct a complete quintic spline from cubic Hermite data.
   *
   * Builds a quintic spline (C2 continuous) by estimating second derivatives
   * from first derivative estimates and solving a tridiagonal system.
   *
   * Algorithm:
   * 1. Estimate first derivatives using specified method (Cubic, PCHIP, Akima, VanLeer)
   * 2. Set up tridiagonal system for second derivatives from quintic continuity conditions
   * 3. Solve system to obtain C2-continuous second derivatives
   *
   * Mathematical formulation:
   * - Quintic spline requires continuity of function, first, and second derivatives
   * - System derived from Hermite interpolation conditions and C2 continuity at knots
   * - Boundary conditions: natural (zero second derivative extrapolation)
   *
   * Reference:
   * - C. de Boor, "A Practical Guide to Splines", Springer, 2001.
   * - D. S. Watkins, "Fundamentals of Matrix Computations", Wiley, 2002.
   *
   * @param[in]  q_sub_type  Method for first derivative estimation
   * @param[in]  X           Array of x-coordinates
   * @param[in]  Y           Array of y-coordinates
   * @param[out] Yp          Array to store first derivatives
   * @param[out] Ypp         Array to store second derivatives
   * @param[in]  npts        Number of data points
   */
  void Quintic_build(
    Spline_sub_type const q_sub_type,
    real_type const       X[],
    real_type const       Y[],
    real_type             Yp[],
    real_type             Ypp[],
    integer const         npts );

  /**
   * @brief Weighted Essentially Non-Oscillatory (WENO5) derivative estimator with Hyman clamp.
   *
   * Computes first derivatives using a fifth-order WENO scheme tailored for non-uniform grids,
   * combined with a Hyman-type monotonicity clamping to preserve the shape of monotone data.
   *
   * This function:
   * - Constructs three candidate stencils (left-biased, central, right-biased) of 3 points.
   * - Uses local smoothness indicators based on second differences to generate nonlinear weights.
   * - Combines candidate derivatives with nonlinear weights for high-order accuracy in smooth regions.
   * - Applies a Hyman monotonicity clamp limiting the result to the range implied by local slopes.
   * - Falls back to simple finite differences for small npts or boundary points.
   *
   * The underlying WENO methodology originates from:
   * - Liu, Osher & Chan: Weighted Essentially Non-Oscillatory Schemes (ENO → WENO)
   * :contentReference[oaicite:1]{index=1}
   * - Jiang & Shu: Efficient Implementation of Weighted ENO Schemes, J. Comput. Phys. 126 (1996)
   * :contentReference[oaicite:2]{index=2}
   *
   * Additional modern improvements to smoothness indicators and weights include:
   * - Modified nonlinear weights for optimal accuracy near critical points :contentReference[oaicite:3]{index=3}
   *
   * The Hyman clamp is adapted from classical slope-limiter strategies to ensure
   * shape preservation similar to the PCHIP/Hyman filtering used in monotone spline methods.
   *
   * @param[in]  X     Array of x-coordinates (must be strictly increasing)
   * @param[in]  Y     Array of y-coordinates
   * @param[out] Yp    Array to store computed first derivatives
   * @param[in]  npts  Number of data points (must be >= 2)
   */

  void WENO5_build( real_type const X[], real_type const Y[], real_type Yp[], integer const npts );

  /**
   * @brief Weighted Essentially Non-Oscillatory (WENO5) second derivative estimator with Hyman clamp.
   *
   * Estimates second derivatives on a non-uniform grid using a fifth-order WENO reconstruction
   * with Hyman-type limiting to preserve local concavity behavior in monotone data.
   *
   * The procedure:
   * - For each interior point, constructs three candidate second derivatives via
   *   second derivatives of local 3-point Lagrange polynomials.
   * - Computes smoothness indicators as squared differences scaled by local grid spacing.
   * - Builds nonlinear weights to favor the smoothest stencils (WENO-JS weights).
   * - Combines candidate second derivatives with nonlinear weights.
   * - Applies a Hyman-type clamp to limit concavity based on adjacent first differences.
   * - Uses fallback finite difference approximations for small npts or at boundaries.
   *
   * The WENO strategy is based on classic literature:
   * - Liu, Osher & Chan: Weighted Essentially Non-Oscillatory Schemes (ENO → WENO)
   * :contentReference[oaicite:4]{index=4}
   * - Jiang & Shu: Efficient Implementation of Weighted ENO Schemes, J. Comput. Phys. 126 (1996)
   * :contentReference[oaicite:5]{index=5}
   *
   * Various research has improved WENO smoothness indicators and weights:
   * - Modified nonlinear weights for improved accuracy at critical points :contentReference[oaicite:6]{index=6}
   * - Mapped WENO and alternative indicators for Hamilton–Jacobi problems :contentReference[oaicite:7]{index=7}
   *
   * Incorporating a Hyman filter on second derivatives draws inspiration from shape-
   * preserving techniques in spline interpolation and CFD limiters.
   *
   * @param[in]  X     Array of x-coordinates (strictly increasing)
   * @param[in]  Y     Array of y-coordinates
   * @param[out] Ypp   Array to store computed second derivatives
   * @param[in]  npts  Number of data points (must be >= 2)
   */
  void WENO5_build2( real_type const X[], real_type const Y[], real_type Ypp[], integer const npts );

  /**
   * @brief Alternative quintic spline construction with WENO-based second derivatives.
   *
   * Similar to Quintic_build but uses WENO5 for second derivative estimation
   * instead of solving a tridiagonal system. Provides potentially better
   * handling of non-smooth data.
   *
   * @param[in]  q_sub_type  Method for first derivative estimation
   * @param[in]  X           Array of x-coordinates
   * @param[in]  Y           Array of y-coordinates
   * @param[out] Yp          Array to store first derivatives
   * @param[out] Ypp         Array to store second derivatives
   * @param[in]  npts        Number of data points
   */
  void Quintic_build2(
    Spline_sub_type const q_sub_type,
    real_type const       X[],
    real_type const       Y[],
    real_type             Yp[],
    real_type             Ypp[],
    integer const         npts );

}  // namespace Splines

#endif

//
// EOF: SplineBuild.cc
//
