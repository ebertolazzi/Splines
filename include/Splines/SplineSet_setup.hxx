/**
 * \brief Build spline set from GenericContainer
 *
 * Constructs a set of splines from a GenericContainer with validated input data.
 * Supports multiple input formats for flexibility and provides comprehensive error checking.
 *
 * \param[in] gc GenericContainer containing spline set configuration
 *
 * \par Required Fields:
 * - "spline_type": vec_string_type - Array of spline type names (e.g., "linear", "cubic")
 * - "xdata": vec_real_type - X-coordinates shared by all splines (must be strictly increasing)
 * - "ydata": Multiple formats supported:
 *   * vec_real_type: Single spline data
 *   * mat_real_type: Matrix where each column is a spline
 *   * vector_type: Vector of vectors, one per spline
 *   * map_type: Map of {spline_name: y_values}
 *
 * \par Optional Fields:
 * - "headers": vec_string_type - Names for splines (required unless ydata is map_type)
 * - "ypdata": map_type - Derivative values for Hermite splines {spline_name: yp_values}
 * - "boundary": vector_type - Boundary conditions per spline (closed, extend, etc.)
 *
 * \par Boundary Configuration:
 * Each boundary element can contain:
 * - "closed": bool - If true, spline wraps around
 * - "extend"/"can_extend": bool - Allow extrapolation beyond data range
 * - "extend_constant": bool - Use constant extrapolation (vs linear)
 *
 * \throw UTILS_ASSERT if required fields missing or data dimensions inconsistent
 * \throw UTILS_ERROR if unsupported data types provided
 *
 * \note CONSTANT splines use n-1 points (interval-based values)
 * \note Warns about unused configuration keys to catch typos
 *
 * \see build() for the underlying construction mechanism
 */
void setup( GenericContainer const & gc );
