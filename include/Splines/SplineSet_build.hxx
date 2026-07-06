///////////////////////////////////////////////////////////////////////////
/**
 * \brief Build the spline set from raw data arrays
 *
 * Constructs multiple splines sharing the same x-nodes. Each spline can be
 * of a different type (linear, cubic, etc.) and may have optional derivative data.
 *
 * \param[in] nspl    Number of splines
 * \param[in] npts    Number of points per spline
 * \param[in] headers Array of spline names (size = nspl)
 * \param[in] stype   Array of spline types (size = nspl)
 * \param[in] X       Array of x-node values (size = npts)
 * \param[in] Y       Array of pointers to y-value arrays (size = nspl)
 * \param[in] Yp      Optional array of pointers to derivative arrays (size = nspl, default: nullptr)
 *
 * \par Memory Management:
 * Allocates internal storage and copies input data. Input arrays can be freed
 * after calling build().
 *
 * \par Supported Spline Types:
 * - CONSTANT: Piecewise constant
 * - LINEAR: Linear interpolation
 * - CUBIC: Cubic spline (natural)
 * - AKIMA: Akima spline
 * - VANLEER: Val Leer limited spline
 * - PCHIP: Piecewise Cubic Hermite Interpolating Polynomial
 * - HERMITE: Hermite spline (requires Yp)
 * - QUINTIC: Quintic spline
 *
 * \throw UTILS_ASSERT if nspl ≤ 0 or npts ≤ 1
 * \throw UTILS_ERROR for unsupported spline types
 *
 * \note For HERMITE splines, Yp must be provided.
 * \note Automatically checks monotonicity for applicable spline types.
 */
void build(
  integer const           nspl,
  integer const           npts,
  char const * const      headers[],
  SplineType1D const      stype[],
  real_type const         X[],
  real_type const * const Y[],
  real_type const * const Yp[] = nullptr );
