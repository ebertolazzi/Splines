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

/**
 * @file SplineUtils.hxx
 * @brief Utility functions for checking monotonicity of cubic and quintic splines
 * @author Enrico Bertolazzi
 * @date 2026
 */

#pragma once

#ifndef SPLINE_UTILS_HXX
#define SPLINE_UTILS_HXX

namespace Splines
{
  /**
   * @brief Check if a cubic spline is monotone
   *
   * This function verifies monotonicity of a cubic Hermite spline using the
   * Fritsch-Carlson conditions. It checks both the input data and the resulting
   * spline interpolant.
   *
   * @param[in] X     Array of x-coordinates (assumed monotonically increasing)
   * @param[in] Y     Array of y-coordinates
   * @param[in] Yp    Array of first derivatives at each point
   * @param[in] npts  Number of data points
   *
   * @return Monotonicity status:
   *         - -2: Non-monotone data (Y values are not monotone)
   *         - -1: Non-monotone spline (violations of Fritsch-Carlson conditions)
   *         -  0: Monotone but not strictly (some flat regions exist)
   *         -  1: Strictly monotone (increasing throughout)
   *
   * @note Based on Fritsch & Carlson (1980), "Monotone Piecewise Cubic Interpolation"
   *       and related work in "Methods of Shape-Preserving Spline Approximation"
   *
   * @see Fritsch, F. N., & Carlson, R. E. (1980). SIAM Journal on Numerical Analysis
   */
  integer check_cubic_spline_monotonicity(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    integer const   npts );

  /**
   * @brief Check if a quintic spline is monotone
   *
   * This function verifies monotonicity of a quintic Hermite spline by checking
   * extended Fritsch-Carlson-type conditions that account for both first and
   * second derivatives.
   *
   * @param[in] X     Array of x-coordinates (assumed monotonically increasing)
   * @param[in] Y     Array of y-coordinates
   * @param[in] Yp    Array of first derivatives at each point
   * @param[in] Ypp   Array of second derivatives at each point
   * @param[in] npts  Number of data points
   *
   * @return Monotonicity status:
   *         - -2: Non-monotone data (Y values are not monotone)
   *         - -1: Non-monotone spline (derivative becomes negative in some interval)
   *         -  0: Monotone but not strictly (some flat regions exist)
   *         -  1: Strictly monotone (increasing throughout)
   *
   * @note Based on shape-preserving quintic Hermite interpolation theory.
   *       For quintic splines, the derivative is a quartic polynomial, making
   *       monotonicity analysis more complex than for cubic splines.
   *
   * @see Botsch & Kobbelt, "Shape preserving quintic Hermite interpolation"
   * @see Nielson, G. M., "Monotonicity preserving interpolation of data by quintic splines"
   */
  integer check_quintic_spline_monotonicity(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    real_type const Ypp[],
    integer const   npts );

}  // namespace Splines

#endif

//
// EOF: SplineUtils.hxx
//
