/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2016                                                      |
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

/*\
 |   ____        _ _          __     __
 |  / ___| _ __ | (_)_ __   __\ \   / /__  ___
 |  \___ \| '_ \| | | '_ \ / _ \ \ / / _ \/ __|
 |   ___) | |_) | | | | | |  __/\ V /  __/ (__
 |  |____/| .__/|_|_|_| |_|\___| \_/ \___|\___|
 |        |_|
\*/

#pragma once

#ifndef SPLINEVEC_HXX
#define SPLINEVEC_HXX

/**
 * \file SplineVec.hxx
 * \brief Multi-dimensional vector spline interpolation using cubic Hermite basis functions
 *
 * This file implements the SplineVec class, which provides piecewise cubic Hermite
 * interpolation for parametric curves in n-dimensional space. The implementation is
 * optimized for 2D and 3D curves but supports arbitrary dimensions.
 *
 * \section theory Theoretical Background
 *
 * The SplineVec class represents a parametric curve γ(t) in ℝⁿ using cubic Hermite
 * interpolation. Each component is interpolated independently using the same parameter t.
 *
 * For a curve with points P₀, P₁, ..., Pₙ and derivatives D₀, D₁, ..., Dₙ, each
 * segment [tᵢ, tᵢ₊₁] uses four Hermite basis functions:
 *
 * H₀(s) = 2s³ - 3s² + 1       (value at start)
 * H₁(s) = -2s³ + 3s²          (value at end)
 * H₂(s) = s³ - 2s² + s        (derivative at start)
 * H₃(s) = s³ - s²             (derivative at end)
 *
 * where s = (t - tᵢ)/(tᵢ₊₁ - tᵢ) ∈ [0,1] is the normalized parameter.
 *
 * \section parametrization Knot Parametrization Schemes
 *
 * The quality of the parametric spline heavily depends on the choice of knot values.
 * This implementation supports several parametrization methods:
 *
 * - **Uniform**: Equal spacing between knots
 * - **Chord Length**: Distance between consecutive points
 * - **Centripetal**: Square root of chord length (reduces loops and cusps)
 * - **Foley**: Advanced method considering angles between segments
 *
 * \section references References
 *
 * -# Catmull, E., & Rom, R. (1974). "A class of local interpolating splines."
 *    In Computer Aided Geometric Design (pp. 317-326). Academic Press.
 *    DOI: 10.1016/B978-0-12-079050-0.50020-5
 *
 * -# Lee, E. T. (1989). "Choosing nodes in parametric curve interpolation."
 *    Computer-Aided Design, 21(6), 363-370.
 *    DOI: 10.1016/0010-4485(89)90003-1
 *
 * -# Foley, T. A., & Nielson, G. M. (1989). "Knot selection for parametric
 *    spline interpolation." Mathematical Methods in Computer Aided Geometric Design.
 *
 * -# de Boor, C. (2001). "A Practical Guide to Splines" (Revised Edition).
 *    Springer-Verlag, New York. ISBN: 978-0-387-95366-3
 *
 * -# Yuksel, C., Schaefer, S., & Keyser, J. (2011). "Parameterization and
 *    applications of Catmull-Rom curves." Computer-Aided Design, 43(7), 747-755.
 *    DOI: 10.1016/j.cad.2010.08.008
 *
 * \author Enrico Bertolazzi
 * \date 2016
 * \version 1.0
 */

namespace Splines
{

  /// Type alias for real matrix from GenericContainer
  using GC_namespace::mat_real_type;

  /**
   * \class SplineVec
   * \brief Vector spline interpolation class for parametric curves in n-dimensional space
   *
   * This class manages a collection of cubic Hermite splines that share the same
   * parameter values (knots) but interpolate different components of a vector-valued
   * function. It is particularly useful for representing smooth parametric curves in
   * 2D, 3D, or higher-dimensional spaces.
   *
   * \par Key Features:
   * - Cubic Hermite interpolation for C¹ continuity
   * - Support for arbitrary dimensions (optimized for 2D and 3D)
   * - Multiple knot parametrization schemes (chord length, centripetal, Foley)
   * - Catmull-Rom automatic derivative computation
   * - Efficient interval search for evaluation
   * - Support for closed curves with periodic boundary conditions
   * - Automatic differentiation support (when AUTODIFF_SUPPORT is defined)
   *
   * \par Memory Layout:
   * The class uses a compact memory layout with custom allocators:
   * - m_X[npts]: Knot values (parameter values)
   * - m_Y[dim][npts]: Point coordinates for each dimension
   * - m_Yp[dim][npts]: Derivative values for each dimension
   *
   * \par Thread Safety:
   * The class is not thread-safe for concurrent modifications, but multiple threads
   * can safely evaluate the spline simultaneously after construction.
   *
   * \par Example Usage:
   * \code{.cpp}
   * // Create a 3D parametric curve through 5 points
   * SplineVec curve("my_curve");
   *
   * // Setup data (5 points in 3D)
   * real_type X[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
   * real_type Y_x[5] = {0.0, 1.0, 2.0, 1.5, 0.5};
   * real_type Y_y[5] = {0.0, 0.5, 1.0, 1.5, 2.0};
   * real_type Y_z[5] = {0.0, 0.2, 0.4, 0.3, 0.1};
   * real_type const* Y[3] = {Y_x, Y_y, Y_z};
   *
   * curve.setup(3, 5, Y);
   * curve.set_knots(X);
   * curve.catmull_rom();  // Compute derivatives automatically
   *
   * // Evaluate curve at parameter t = 1.5
   * real_type x = curve.eval(1.5, 0);  // x-component
   * real_type y = curve.eval(1.5, 1);  // y-component
   * real_type z = curve.eval(1.5, 2);  // z-component
   *
   * // Or evaluate all components at once
   * std::vector<real_type> point;
   * curve.eval(1.5, point);
   * \endcode
   *
   * \see Hermite3, Hermite3_D, Hermite3_DD, Hermite3_DDD for basis function details
   */
  class SplineVec final
  {
  protected:
    /// Spline identifier name for error reporting and debugging
    string const m_name;

    /// Custom memory allocator for real-type data (contiguous storage)
    Utils::Malloc<real_type> m_mem;

    /// Custom memory allocator for real-type pointers (array of arrays)
    Utils::Malloc<real_type *> m_mem_p;

    /// Number of spatial dimensions (components) in the vector spline
    integer m_dim = 0;

    /// Number of interpolation points (knots + 1)
    integer m_npts = 0;

    /// Flag indicating if the curve forms a closed loop (periodic boundary conditions)
    bool m_curve_is_closed = false;

    /// Flag indicating if extrapolation is allowed outside the knot range
    bool m_curve_can_extend = true;

    /// Knot vector (parameter values): size = m_npts
    real_type * m_X = nullptr;

    /// Point coordinates: m_Y[i][j] = i-th component of j-th point, size = m_dim × m_npts
    real_type ** m_Y = nullptr;

    /// Derivative values: m_Yp[i][j] = i-th component of j-th derivative, size = m_dim × m_npts
    real_type ** m_Yp = nullptr;

    /// Efficient interval search structure for finding the correct spline segment
    Utils::SearchInterval<real_type, integer> m_search;

    /**
     * \brief Allocate memory for spline data structures
     *
     * This method allocates contiguous memory blocks for storing spline data.
     * It uses custom allocators to ensure efficient memory usage and cache locality.
     *
     * \param[in] dim  Number of spatial dimensions (must be > 0)
     * \param[in] npts Number of interpolation points (must be > 1)
     *
     * \throw UTILS_ASSERT if dim ≤ 0 or npts ≤ 1
     *
     * \par Memory Allocation Strategy:
     * - Single allocation for all knots and values: (2*dim + 1)*npts real_type values
     * - Single allocation for all pointers: 2*dim real_type* pointers
     * - Memory layout: [X | Y[0] | Y[1] | ... | Y[dim-1] | Yp[0] | ... | Yp[dim-1]]
     *
     * \par Complexity:
     * Time: O(dim × npts) for initialization
     * Space: O(dim × npts)
     */
    void allocate( integer const dim, integer const npts );

    /**
     * \brief Compute chord lengths between consecutive points
     *
     * This method calculates the Euclidean distance between consecutive points
     * and stores them in m_X. This is used as an intermediate step in parametrization
     * schemes like chord length and centripetal methods.
     *
     * The implementation is optimized for common cases (2D and 3D) using std::hypot
     * for improved numerical stability.
     *
     * \par Algorithm:
     * For each segment [i, i+1]:
     * \f[
     * d_i = \sqrt{\sum_{k=0}^{dim-1} (Y_k[i+1] - Y_k[i])^2}
     * \f]
     *
     * \par Numerical Stability:
     * - Uses std::hypot for 2D and 3D cases to avoid overflow/underflow
     * - Direct computation for higher dimensions
     *
     * \par Complexity:
     * Time: O(dim × npts)
     * Space: O(1) auxiliary
     *
     * \note After calling this method, m_X contains chord lengths, not final knot values.
     *       Call set_knots_chord_length() or set_knots_centripetal() to normalize.
     */
    void compute_chords();

  public:
    /// Deleted copy constructor (non-copyable)
    SplineVec( SplineVec const & ) = delete;

    /// Deleted copy assignment operator (non-copyable)
    SplineVec const & operator=( SplineVec const & ) = delete;

    //!
    //! \name Constructors and Destructor
    //!
    ///@{

    /**
     * \brief Construct an empty SplineVec object
     *
     * Creates an uninitialized spline vector. You must call setup() or build()
     * before using evaluation methods.
     *
     * \param[in] name Identifier name for debugging and error messages (default: "SplineVec")
     *
     * \par Post-conditions:
     * - m_dim = 0, m_npts = 0
     * - Memory allocators are initialized but empty
     * - Search interval structure is configured
     */
    explicit SplineVec( string_view name = "SplineVec" )
      : m_name( name )
      , m_mem( fmt::format( "SplineVec[{}]::m_mem", name ) )
      , m_mem_p( fmt::format( "SplineVec[{}]::m_mem_p", name ) )
    { m_search.setup( &m_name, &m_npts, &m_X, &m_curve_is_closed, &m_curve_can_extend ); }

    /**
     * \brief Destructor - releases all allocated memory
     *
     * Automatically frees all memory allocated by the custom allocators.
     * No manual cleanup is required.
     */
    ~SplineVec()
    {
      m_mem.free();
      m_mem_p.free();
    }

    ///@}


    /**
     * \brief Get the spline identifier name
     *
     * \return String view of the spline name set in the constructor
     */
    [[nodiscard]] string_view name() const noexcept { return m_name; }

    //!
    //! \name Curve Topology Control
    //!
    ///@{

    /**
     * \brief Check if the curve is closed (periodic)
     *
     * \return true if the curve forms a closed loop with periodic boundary conditions
     *
     * \par Closed Curves:
     * When a curve is closed, evaluation at parameter values outside [t₀, tₙ]
     * wraps around cyclically: t → t mod (tₙ - t₀)
     */
    [[nodiscard]] bool is_closed() const noexcept { return m_curve_is_closed; }

    /**
     * \brief Set the curve as closed (periodic)
     *
     * Enables periodic boundary conditions. Parameters outside the knot range
     * will be wrapped cyclically before evaluation.
     *
     * \par Use Cases:
     * - Closed contours and loops
     * - Periodic motion paths
     * - Torus-like geometries
     */
    void make_closed() { m_curve_is_closed = true; }

    /**
     * \brief Set the curve as open (non-periodic)
     *
     * Disables periodic boundary conditions. Evaluation outside the knot range
     * will either extrapolate (if can_extend is true) or raise an error.
     */
    void make_open() { m_curve_is_closed = false; }

    /**
     * \brief Check if extrapolation is allowed
     *
     * \return true if evaluation can extend beyond the knot range
     */
    [[nodiscard]] bool can_extend() const noexcept { return m_curve_can_extend; }

    /**
     * \brief Enable extrapolation beyond knot range
     *
     * When enabled, parameters outside [t₀, tₙ] use the tangent at the
     * nearest endpoint for linear extrapolation.
     *
     * \warning Extrapolation quality degrades rapidly away from the knot range
     */
    void make_unbounded() { m_curve_can_extend = true; }

    /**
     * \brief Disable extrapolation (bounded curve)
     *
     * When disabled, evaluation outside [t₀, tₙ] raises an error.
     * This is the safe default for most applications.
     */
    void make_bounded() { m_curve_can_extend = false; }
    ///@}


    //!
    //! \name Spline Information
    //!
    ///@{

    /**
     * \brief Get the number of interpolation points
     *
     * \return Number of knots/control points defining the spline
     */
    [[nodiscard]] integer num_points() const noexcept { return m_npts; }

    /**
     * \brief Get the spatial dimension
     *
     * \return Number of components in the vector-valued function
     */
    [[nodiscard]] integer dimension() const noexcept { return m_dim; }

    /**
     * \brief Get pointer to knot vector
     *
     * \return Pointer to array of m_npts knot values (read-only)
     *
     * \note The returned pointer is valid until the spline is modified or destroyed
     */
    [[nodiscard]] real_type const * x_nodes() const noexcept { return m_X; }

    /**
     * \brief Get a specific knot value
     *
     * \param[in] npt Knot index (must be in [0, m_npts-1])
     * \return Parameter value at the specified knot
     */
    [[nodiscard]] real_type x_node( integer const npt ) const noexcept { return m_X[npt]; }

    /**
     * \brief Get pointer to point coordinates for a specific dimension
     *
     * \param[in] j Dimension index (must be in [0, m_dim-1])
     * \return Pointer to array of m_npts coordinate values for dimension j
     */
    [[nodiscard]] real_type const * y_nodes( integer const j ) const noexcept { return m_Y[j]; }

    /**
     * \brief Get a specific point coordinate
     *
     * \param[in] npt Point index (must be in [0, m_npts-1])
     * \param[in] j   Dimension index (must be in [0, m_dim-1])
     * \return Coordinate value: Y_j[npt]
     */
    [[nodiscard]] real_type y_node( integer const npt, integer const j ) const noexcept { return m_Y[j][npt]; }

    /**
     * \brief Get minimum parameter value
     *
     * \return First knot value (lower bound of parameter domain)
     */
    [[nodiscard]] real_type x_min() const noexcept { return m_X[0]; }

    /**
     * \brief Get maximum parameter value
     *
     * \return Last knot value (upper bound of parameter domain)
     */
    [[nodiscard]] real_type x_max() const noexcept { return m_X[m_npts - 1]; }

    ///@}

    //!
    //! \name Single Component Evaluation
    //!
    ///@{

    /**
     * \brief Evaluate a single component of the spline
     *
     * Computes the i-th component of the curve at parameter t using cubic Hermite
     * interpolation. The method automatically finds the correct spline segment
     * using an efficient search algorithm.
     *
     * \param[in] x Parameter value (in the knot range or extrapolated if allowed)
     * \param[in] i Component index (must be in [0, m_dim-1])
     * \return Interpolated value of the i-th component at parameter x
     *
     * \par Algorithm:
     * 1. Find segment index k such that X[k] ≤ x < X[k+1]
     * 2. Normalize parameter: s = (x - X[k])/(X[k+1] - X[k])
     * 3. Compute Hermite basis functions H₀(s), H₁(s), H₂(s), H₃(s)
     * 4. Return: H₀·Y[i][k] + H₁·Y[i][k+1] + H₂·Yp[i][k] + H₃·Yp[i][k+1]
     *
     * \par Complexity:
     * Time: O(log npts) average case for segment search
     * Space: O(1)
     *
     * \see operator(), eval_D(), eval_DD()
     */
    [[nodiscard]] real_type eval( real_type const x, integer const i ) const;

    /**
     * \brief Function call operator for single component evaluation
     *
     * Convenience operator that calls eval(x, i).
     *
     * \param[in] x Parameter value
     * \param[in] i Component index
     * \return curve(x)[i]
     */
    [[nodiscard]] real_type operator()( real_type const x, integer const i ) const { return this->eval( x, i ); }

    /**
     * \brief Evaluate first derivative of a single component at parameter x
     *
     * Computes the first derivative (velocity) of the i-th component of the
     * vector-valued spline at parameter x using the derivative of cubic Hermite
     * basis functions.
     *
     * \param[in] x Parameter value within the knot range (or extrapolated if allowed)
     * \param[in] i Component index (0 ≤ i < m_dim)
     * \return First derivative dY_i/dt at parameter x
     *
     * \par Mathematical Details:
     * For segment [t_k, t_{k+1}] with normalized parameter s = (x - t_k)/Δt:
     * \f[
     * \frac{dY_i}{dt} = H_0'(s)·Y_i[k] + H_1'(s)·Y_i[k+1] + H_2'(s)·Y'_i[k] + H_3'(s)·Y'_i[k+1]
     * \f]
     * where H_j'(s) are derivatives of Hermite basis functions scaled by 1/Δt.
     *
     * \par Complexity:
     * Time: O(log npts) for interval search
     * Space: O(1)
     *
     * \see eval(), DD(), DDD()
     */
    [[nodiscard]] real_type D( real_type const x, integer const i ) const;

    /**
     * \brief Evaluate first derivative of a single component
     *
     * \copydetails D()
     * This is an alias for D() method.
     */
    [[nodiscard]] real_type eval_D( real_type const x, integer const i ) const { return this->D( x, i ); }

    /**
     * \brief Evaluate second derivative of a single component at parameter x
     *
     * Computes the second derivative (acceleration) of the i-th component of the
     * vector-valued spline at parameter x using the second derivative of cubic
     * Hermite basis functions.
     *
     * \param[in] x Parameter value within the knot range (or extrapolated if allowed)
     * \param[in] i Component index (0 ≤ i < m_dim)
     * \return Second derivative d²Y_i/dt² at parameter x
     *
     * \par Mathematical Details:
     * For segment [t_k, t_{k+1}] with normalized parameter s = (x - t_k)/Δt:
     * \f[
     * \frac{d^2Y_i}{dt^2} = H_0''(s)·Y_i[k] + H_1''(s)·Y_i[k+1] + H_2''(s)·Y'_i[k] + H_3''(s)·Y'_i[k+1]
     * \f]
     * where H_j''(s) are second derivatives of Hermite basis functions scaled by 1/Δt².
     *
     * \par Complexity:
     * Time: O(log npts) for interval search
     * Space: O(1)
     *
     * \see eval(), D(), DDD()
     */
    [[nodiscard]] real_type DD( real_type const x, integer const i ) const;

    /**
     * \brief Evaluate second derivative of a single component
     *
     * \copydetails DD()
     * This is an alias for DD() method.
     */
    [[nodiscard]] real_type eval_DD( real_type const x, integer const i ) const { return this->DD( x, i ); }

    /**
     * \brief Evaluate third derivative of a single component at parameter x
     *
     * Computes the third derivative (jerk) of the i-th component of the
     * vector-valued spline at parameter x using the third derivative of cubic
     * Hermite basis functions.
     *
     * \param[in] x Parameter value within the knot range (or extrapolated if allowed)
     * \param[in] i Component index (0 ≤ i < m_dim)
     * \return Third derivative d³Y_i/dt³ at parameter x
     *
     * \par Mathematical Details:
     * For segment [t_k, t_{k+1}] with normalized parameter s = (x - t_k)/Δt:
     * \f[
     * \frac{d^3Y_i}{dt^3} = H_0'''(s)·Y_i[k] + H_1'''(s)·Y_i[k+1] + H_2'''(s)·Y'_i[k] + H_3'''(s)·Y'_i[k+1]
     * \f]
     * where H_j'''(s) are third derivatives of Hermite basis functions scaled by 1/Δt³.
     *
     * \par Complexity:
     * Time: O(log npts) for interval search
     * Space: O(1)
     *
     * \see eval(), D(), DD()
     */
    [[nodiscard]] real_type DDD( real_type const x, integer const i ) const;

    /**
     * \brief Evaluate third derivative of a single component
     *
     * \copydetails DDD()
     * This is an alias for DDD() method.
     */
    [[nodiscard]] real_type eval_DDD( real_type const x, integer const i ) const { return this->DDD( x, i ); }

    /**
     * \brief Evaluate fourth derivative of a single component at parameter x
     *
     * Returns the fourth derivative of the i-th component, which is always zero
     * for cubic Hermite splines (piecewise cubic polynomials).
     *
     * \param[in] x Parameter value (unused, included for interface consistency)
     * \param[in] i Component index (unused, included for interface consistency)
     * \return Always returns 0.0
     *
     * \par Mathematical Justification:
     * Cubic Hermite splines are piecewise cubic polynomials, whose fourth derivative
     * is identically zero everywhere except at knots (where it's undefined).
     *
     * \note This method exists for interface completeness with higher-order spline
     *       classes.
     */
    [[nodiscard]] real_type DDDD( [[maybe_unused]] real_type const x, [[maybe_unused]] integer const i ) const
    {
      // Hermite3 (cubic) spline has zero 4th derivative
      return 0.0;
    }

    /**
     * \brief Evaluate fourth derivative of a single component
     *
     * \copydetails DDDD()
     * This is an alias for DDDD() method.
     */
    [[nodiscard]] real_type eval_DDDD( [[maybe_unused]] real_type const x, [[maybe_unused]] integer const i ) const
    { return this->DDDD( x, i ); }

    /**
     * \brief Evaluate fifth derivative of a single component at parameter x
     *
     * Returns the fifth derivative of the i-th component, which is always zero
     * for cubic Hermite splines (piecewise cubic polynomials).
     *
     * \param[in] x Parameter value (unused, included for interface consistency)
     * \param[in] i Component index (unused, included for interface consistency)
     * \return Always returns 0.0
     *
     * \par Mathematical Justification:
     * Cubic Hermite splines are piecewise cubic polynomials, whose fifth derivative
     * is identically zero everywhere (derivative of the constant zero fourth derivative).
     *
     * \note This method exists for interface completeness with higher-order spline
     *       classes.
     */
    [[nodiscard]] real_type DDDDD( [[maybe_unused]] real_type const x, [[maybe_unused]] integer const i ) const
    {
      // Hermite3 (cubic) spline has zero 5th derivative
      return 0.0;
    }

    /**
     * \brief Evaluate fifth derivative of a single component
     *
     * \copydetails DDDDD()
     * This is an alias for DDDDD() method.
     */
    [[nodiscard]] real_type eval_DDDDD( real_type const x, integer const i ) const { return this->DDDDD( x, i ); }

    ///@}

#ifdef AUTODIFF_SUPPORT
    //!
    //! \name Automatic Differentiation Support
    //!
    ///@{

    /**
     * \brief Evaluate with first-order automatic differentiation
     *
     * Computes both the function value and its derivative with respect to the
     * parameter using autodiff::dual1st types.
     *
     * \param[in] x Dual number parameter (contains value and gradient)
     * \param[in] i Component index
     * \return Dual number result with value and gradient
     *
     * \par Usage:
     * \code{.cpp}
     * autodiff::dual1st t = 1.5;
     * autodiff::dual1st y = spline.eval(t, 0);
     * // y.val = curve value at t=1.5
     * // y.grad = derivative dY/dt at t=1.5
     * \endcode
     */
    autodiff::dual1st eval( autodiff::dual1st const & x, integer const i ) const
    {
      real_type         xv = x.val;
      autodiff::dual1st res;
      res.val  = eval( xv, i );
      res.grad = eval_D( xv, i ) * x.grad;
      return res;
    }

    /**
     * \brief Evaluate with second-order automatic differentiation
     *
     * Computes the function value and both first and second derivatives.
     *
     * \param[in] x Dual number parameter (second order)
     * \param[in] i Component index
     * \return Dual number with value, first derivative, and second derivative
     *
     * \par Implementation:
     * Uses chain rule for second derivatives:
     * \f[
     * \frac{d^2y}{dx^2} = \frac{dy}{dt}\frac{d^2t}{dx^2} + \frac{d^2y}{dt^2}\left(\frac{dt}{dx}\right)^2
     * \f]
     */
    autodiff::dual2nd eval( autodiff::dual2nd const & x, integer const i ) const
    {
      // Estrai il valore reale di x e la sua derivata prima
      real_type xv = x.val.val;   // valore reale di x
      real_type xg = x.grad.val;  // derivata prima di x: dx/dX

      // Calcola la funzione e le sue derivate nel punto xv
      real_type f_val = eval( xv, i );     // valore della funzione f(x)
      real_type dfx   = eval_D( xv, i );   // f'(x) = df/dx
      real_type ddfx  = eval_DD( xv, i );  // f''(x) = d²f/dx²

      // Costruisci l'oggetto dual2nd risultante
      autodiff::dual2nd res;

      // Imposta il valore della funzione
      res.val.val = f_val;

      // Calcola la derivata prima: dF/dX = f'(x) * (dx/dX)
      real_type dF_dX = dfx * xg;

      // Imposta la derivata prima
      res.val.grad = dF_dX;  // per coerenza interna
      res.grad.val = dF_dX;  // derivata prima

      // Calcola la derivata seconda:
      // d²F/dX² = f''(x) * (dx/dX)² + f'(x) * d²x/dX²
      real_type d2F_dX2 = ddfx * ( xg * xg ) + dfx * x.grad.grad;

      // Imposta la derivata seconda
      res.grad.grad = d2F_dX2;

      return res;
    }

    /**
     * \brief Generic template for automatic differentiation
     *
     * Automatically selects the appropriate evaluation method based on the
     * type of the input parameter.
     *
     * \tparam T Input type (arithmetic or autodiff dual type)
     * \param[in] x Parameter value (can be double, dual1st, dual2nd, etc.)
     * \param[in] i Component index
     * \return Appropriate result type matching input
     *
     * \par Type Detection:
     * - If T is arithmetic (int, float, double): promotes to real_type
     * - Otherwise: uses autodiff automatic type deduction
     */
    template <typename T> auto eval( T const & x, integer const i ) const
    {
      if constexpr ( std::is_arithmetic<T>::value )
      {
        // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
        return eval( static_cast<real_type>( x ), i );
      }
      else
      {
        // Altrimenti deduce automaticamente il tipo duale appropriato
        return eval( autodiff::detail::to_dual( x ), i );
      }
    }

    /**
     * \brief Generic function call operator for autodiff
     *
     * \tparam T Parameter type
     * \param[in] x Parameter value
     * \param[in] i Component index
     * \return Result matching input type
     */
    template <typename T> auto operator()( T const & x, integer const i ) const -> decltype( eval( x, i ) )
    { return eval( x, i ); }
///@}
#endif

    //!
    //! \name Vector Evaluation (All Components)
    //!
    ///@{

    /**
     * \brief Evaluate all components at once with custom stride
     *
     * Efficient method to evaluate all dimensions of the vector spline at a
     * single parameter value. Results are stored in a pre-allocated array
     * with configurable stride for flexible memory layouts.
     *
     * \param[in]  x    Parameter value
     * \param[out] vals Output array (must have space for at least dim*inc values)
     * \param[in]  inc  Stride between consecutive components (default: 1)
     *
     * \par Memory Layout:
     * - inc = 1: vals = [Y₀, Y₁, ..., Y_{dim-1}] (contiguous)
     * - inc = n: vals[i*n] = Y_i (every n-th element)
     *
     * \par Use Cases:
     * - Interleaved arrays: inc = 1
     * - Structure-of-arrays: inc = row_stride
     * - Matrix columns: inc = num_rows
     *
     * \par Performance:
     * Single segment search followed by dim evaluations.
     * Time: O(log npts + dim)
     */
    void eval( real_type const x, real_type vals[], integer const inc ) const;

    /**
     * \brief Evaluate first derivatives of all components
     *
     * \param[in]  x    Parameter value
     * \param[out] vals Output array for derivatives
     * \param[in]  inc  Stride between components
     *
     * \par Output:
     * vals[i*inc] = dY_i/dt for i = 0, ..., dim-1
     */
    void eval_D( real_type const x, real_type vals[], integer const inc ) const;

    /**
     * \brief Evaluate second derivatives of all components with stride
     *
     * Computes the second derivatives (accelerations) of all vector components
     * at parameter x and stores them in an output array with specified stride.
     *
     * \param[in]  x    Parameter value within the knot range
     * \param[out] vals Output array for second derivatives (size ≥ dim×inc)
     * \param[in]  inc  Stride between consecutive components in output array
     *
     * \par Output Layout:
     * vals[i*inc] = d²Y_i/dt² for i = 0, ..., dim-1
     *
     * \par Complexity:
     * Time: O(log npts + dim) - single search followed by dim evaluations
     * Space: O(1) auxiliary
     */
    void eval_DD( real_type const x, real_type vals[], integer const inc ) const;

    /**
     * \brief Evaluate third derivative of all components with stride
     *
     * Computes the third derivatives (jerks) of all vector components
     * at parameter x and stores them in an output array with specified stride.
     *
     * \param[in]  x    Parameter value within the knot range
     * \param[out] vals Output array for third derivatives (size ≥ dim×inc)
     * \param[in]  inc  Stride between consecutive components in output array
     *
     * \par Output Layout:
     * vals[i*inc] = d³Y_i/dt³ for i = 0, ..., dim-1
     *
     * \par Complexity:
     * Time: O(log npts + dim) - single search followed by dim evaluations
     * Space: O(1) auxiliary
     */
    void eval_DDD( real_type const x, real_type vals[], integer const inc ) const;

    /**
     * \brief Evaluate fourth derivatives (all zeros)
     *
     * \param[in]  x    Parameter value (unused)
     * \param[out] vals Output array filled with zeros
     * \param[in]  inc  Stride between components
     */
    void eval_DDDD( [[maybe_unused]] real_type const x, [[maybe_unused]] real_type vals[], integer const inc ) const
    {
      // Hermite3 (cubic) spline has zero 4th derivative
      real_type * v = vals;
      for ( integer j = 0; j < m_dim; ++j, v += inc ) *v = 0.0;
    }

    /**
     * \brief Evaluate fifth derivatives (all zeros)
     *
     * \param[in]  x    Parameter value (unused)
     * \param[out] vals Output array filled with zeros
     * \param[in]  inc  Stride between components
     */
    void eval_DDDDD( [[maybe_unused]] real_type const x, [[maybe_unused]] real_type vals[], integer const inc ) const
    {
      // Hermite3 (cubic) spline has zero 5th derivative
      real_type * v = vals;
      for ( integer j = 0; j < m_dim; ++j, v += inc ) *v = 0.0;
    }
    ///@}

    //!
    //! \name STL Vector Interface
    //!
    ///@{

    /**
     * \brief Evaluate all components into std::vector
     *
     * Convenience method that automatically resizes the output vector and
     * fills it with all component values.
     *
     * \param[in]  x    Parameter value
     * \param[out] vals Output vector (automatically resized to dim)
     *
     * \par Example:
     * \code{.cpp}
     * std::vector<double> point;
     * spline.eval(1.5, point);
     * // point.size() == spline.dimension()
     * // point[i] == i-th component value
     * \endcode
     */
    void eval( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_dim );
      eval( x, vals.data(), 1 );
    }

    /**
     * \brief Evaluate first derivatives into std::vector
     *
     * \param[in]  x    Parameter value
     * \param[out] vals Output vector (resized and filled with derivatives)
     */
    void eval_D( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_dim );
      eval_D( x, vals.data(), 1 );
    }

    /**
     * \brief Evaluate the second derivative of all the splines at `x` and
     *        store values in `vals`.
     */
    void eval_DD( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_dim );
      eval_DD( x, vals.data(), 1 );
    }

    /**
     * \brief Evaluate the third derivative of all the splines at `x` and
     *        store values in `vals`.
     */
    void eval_DDD( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_dim );
      eval_DDD( x, vals.data(), 1 );
    }

    /**
     * \brief Evaluate the 4th derivative of all the splines at `x` and
     *        store values in `vals`.
     */
    void eval_DDDD( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_dim );
      eval_DDDD( x, vals.data(), 1 );
    }

    /**
     * \brief Evaluate the 5th derivative of all the splines at `x` and
     *        store values in `vals`.
     */
    void eval_DDDDD( real_type const x, vector<real_type> & vals ) const
    {
      vals.resize( m_dim );
      eval_DDDDD( x, vals.data(), 1 );
    }
    ///@}

    //!
    //! \name GenericContainer Interface
    //!
    ///@{

    /**
     * \brief Evaluate at single point into GenericContainer
     *
     * \param[in]  x    Parameter value
     * \param[out] vals GenericContainer (set as vector of reals with dim elements)
     */
    void eval( real_type const x, GenericContainer & vals ) const { eval( x, vals.set_vec_real( m_dim ) ); }

    /**
     * \brief Evaluate first derivatives at `x` and fill a `GenericContainer`
     */
    void eval_D( real_type const x, GenericContainer & vals ) const { eval_D( x, vals.set_vec_real( m_dim ) ); }

    /**
     * \brief Evaluate second derivatives at `x` and fill a `GenericContainer`
     */
    void eval_DD( real_type const x, GenericContainer & vals ) const { eval_DD( x, vals.set_vec_real( m_dim ) ); }

    /**
     * \brief Evaluate third derivatives at `x` and fill a `GenericContainer`
     */
    void eval_DDD( real_type const x, GenericContainer & vals ) const { eval_DDD( x, vals.set_vec_real( m_dim ) ); }

    /**
     * \brief Evaluate 4th derivatives at `x` and fill a `GenericContainer`
     */
    void eval_DDDD( real_type const x, GenericContainer & vals ) const { eval_DDDD( x, vals.set_vec_real( m_dim ) ); }

    /**
     * \brief Evaluate 5th derivatives at `x` and fill a `GenericContainer`
     */
    void eval_DDDDD( real_type const x, GenericContainer & vals ) const { eval_DDDDD( x, vals.set_vec_real( m_dim ) ); }
    ///@}

    //!
    //! \name Batch Evaluation
    //!
    ///@{

    /**
     * \brief Evaluate at multiple parameters
     *
     * Efficiently evaluates the spline at multiple parameter values and stores
     * results in a matrix format within a GenericContainer.
     *
     * \param[in]  x    Vector of parameter values
     * \param[out] vals GenericContainer (set as dim × x.size() matrix)
     *
     * \par Output Matrix Layout:
     * - Row i: All values for component i across all parameters
     * - Column j: All component values at parameter x[j]
     * - vals(i, j) = Y_i(x[j])
     *
     * \par Performance:
     * Time: O(x.size() × (log npts + dim))
     *
     * \par Use Cases:
     * - Plotting curves with many sample points
     * - Batch processing for optimization
     * - Monte Carlo simulations
     */
    void eval( vec_real_type const & x, GenericContainer & vals ) const;

    /**
     * \brief Batch evaluate first derivatives
     *
     * \param[in]  x    Vector of parameter values
     * \param[out] vals GenericContainer (dim × x.size() matrix of derivatives)
     */
    void eval_D( vec_real_type const & x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate second derivative at `x` (x is a vector with many values)
     *        and fill a GenericContainer
     */
    void eval_DD( vec_real_type const & x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate third derivative at `x` (x is a vector with many values)
     *        and fill a GenericContainer
     */
    void eval_DDD( vec_real_type const & x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate 4th derivative at `x` (x is a vector with many values)
     *        and fill a GenericContainer
     */
    void eval_DDDD( vec_real_type const & x, GenericContainer & vals ) const;

    /**
     * \brief Evaluate 5th derivative at `x` (x is a vector with many values)
     *        and fill a GenericContainer
     */
    void eval_DDDDD( vec_real_type const & x, GenericContainer & vals ) const;
    ///@}

    //!
    //! \name Spline Construction and Setup
    //!
    ///@{

    /**
     * \brief Initialize spline with array of pointers
     *
     * Sets up the spline data structure with point coordinates provided as
     * an array of pointers, one pointer per dimension.
     *
     * \param[in] dim  Number of spatial dimensions
     * \param[in] npts Number of interpolation points
     * \param[in] Y    Array of pointers: Y[i] points to npts values for dimension i
     *
     * \par Memory:
     * Data is copied into internal storage. The input arrays can be freed after setup.
     *
     * \par Example:
     * \code{.cpp}
     * double x_coords[5] = {0, 1, 2, 3, 4};
     * double y_coords[5] = {0, 1, 0, -1, 0};
     * double z_coords[5] = {0, 0, 1, 0, -1};
     * double const* data[3] = {x_coords, y_coords, z_coords};
     *
     * SplineVec curve;
     * curve.setup(3, 5, data);  // 3D curve with 5 points
     * \endcode
     *
     * \note You must call set_knots() and catmull_rom() (or compute derivatives
     *       manually) before evaluation.
     */
    void setup( integer const dim, integer const npts, real_type const * const Y[] );

    /**
     * \brief Initialize spline with matrix data
     *
     * Sets up the spline from a matrix stored in column-major (Fortran) format.
     *
     * \param[in] dim  Number of spatial dimensions
     * \param[in] npts Number of interpolation points
     * \param[in] Y    Matrix data in column-major format
     * \param[in] ldY  Leading dimension (row stride) of matrix
     *
     * \par Matrix Layout:
     * Y is a dim × npts matrix stored column-major:
     * - Y[i + j*ldY] = j-th point, i-th component
     * - ldY ≥ dim (may have padding rows)
     *
     * \par Use Cases:
     * - Interfacing with FORTRAN/LAPACK libraries
     * - Reading from column-major matrices
     */
    void setup( integer const dim, integer const npts, real_type const Y[], integer const ldY );

    /**
     * \brief Deep copy to another SplineVec
     *
     * Creates a complete independent copy of this spline in the target object.
     * Only copies point coordinates; derivatives and knots must be recomputed
     * in the target object.
     *
     * \param[out] S Target SplineVec (will be overwritten)
     *
     * \note This method performs a deep copy of point coordinates only.
     *       The target spline will need to have knots set and derivatives
     *       computed before it can be evaluated.
     */
    void deep_copy_to( SplineVec & S ) const { S.setup( m_dim, m_npts, m_Y ); }

    /**
     * \brief Manually set knot values
     *
     * Explicitly provides the parameter values for each interpolation point.
     *
     * \param[in] X Array of npts knot values (must be strictly increasing)
     *
     * \warning Knots must be in strictly increasing order: X[i] < X[i+1]
     *
     * \par When to Use:
     * - Custom parametrization schemes
     * - Time-based parameters
     * - Arc-length parametrization
     *
     * \see set_knots_chord_length(), set_knots_centripetal()
     */
    void set_knots( real_type const X[] )
    {
      if ( m_npts > 0 ) std::memcpy( m_X, X, m_npts * sizeof( *m_X ) );
      m_search.must_reset();
    }

    /**
     * \brief Compute knots using chord length parametrization
     *
     * Sets knots proportional to the Euclidean distance between consecutive points.
     * This is the most common parametrization for general curves.
     *
     * \par Algorithm:
     * 1. Compute chord lengths: d_i = ||P_{i+1} - P_i||
     * 2. Set cumulative distances: t_i = Σ_{j=0}^{i-1} d_j
     * 3. Normalize to [0, 1]: t_i ← t_i / t_{n-1}
     *
     * \par Properties:
     * - Simple and intuitive
     * - Works well for smooth curves
     * - May create loops for curves with sharp corners
     *
     * \par Reference:
     * Lee, E. T. (1989). "Choosing nodes in parametric curve interpolation."
     * Computer-Aided Design, 21(6), 363-370.
     */
    void set_knots_chord_length();

    /**
     * \brief Compute knots using centripetal parametrization
     *
     * Sets knots proportional to the square root of chord length. This method
     * reduces loops and cusps compared to chord length parametrization.
     *
     * \par Algorithm:
     * 1. Compute chord lengths: d_i = ||P_{i+1} - P_i||
     * 2. Apply square root: l_i = √d_i
     * 3. Set cumulative values: t_i = Σ_{j=0}^{i-1} l_j
     * 4. Normalize to [0, 1]: t_i ← t_i / t_{n-1}
     *
     * \par Properties:
     * - Better for curves with varying point density
     * - Reduces self-intersections and loops
     * - Recommended for data from physical measurements
     * - More robust than chord length for sharp corners
     *
     * \par Reference:
     * Lee, E. T. (1989). "Choosing nodes in parametric curve interpolation."
     * Computer-Aided Design, 21(6), 363-370.
     *
     * \par Mathematical Note:
     * The exponent 0.5 (square root) is optimal for reducing cusps in typical
     * applications. This is known as the "centripetal" scheme.
     */
    void set_knots_centripetal();

    /**
     * \brief Compute knots using Foley's method
     *
     * Advanced parametrization that considers both chord lengths and angles
     * between consecutive segments for optimal knot placement.
     *
     * \par Algorithm:
     * 1. Compute chord lengths: d_i = ||P_{i+1} - P_i||
     * 2. Compute turning angles at interior points
     * 3. Compute weighted chord lengths using angle weights
     * 4. Set cumulative values and normalize to [0, 1]
     *
     * \par Formula:
     * t_i = Σ_{j=0}^{i-1} (d_j * w_j)
     * where w_j = 1 + (3/2) * (θ_j * d_{j-1})/(d_{j-1} + d_j)
     *              + (3/2) * (θ_{j+1} * d_j)/(d_j + d_{j+1})
     *
     * \par Reference:
     * Foley, T. A., & Nielson, G. M. (1989). "Knot selection for parametric
     * spline interpolation." Mathematical Methods in Computer Aided Geometric Design.
     */
    void set_knots_foley();

    /**
     * \brief Compute derivatives using Catmull-Rom algorithm
     *
     * Computes derivatives at each knot using the Catmull-Rom scheme,
     * which produces smooth C¹ continuous curves that pass through all
     * control points.
     *
     * \par Algorithm:
     * For internal points (1 ≤ i ≤ n-2):
     * \f[
     * P'_i = \frac{\Delta t_i}{\Delta t_{i-1} + \Delta t_i} \cdot \frac{P_{i+1} - P_i}{\Delta t_i}
     *      + \frac{\Delta t_{i-1}}{\Delta t_{i-1} + \Delta t_i} \cdot \frac{P_i - P_{i-1}}{\Delta t_{i-1}}
     * \f]
     *
     * For endpoints, special formulas are used to maintain smoothness.
     *
     * \par Properties:
     * - Local control: each derivative depends only on neighboring points
     * - Guaranteed C¹ continuity
     * - No overshoot or oscillations for well-behaved data
     *
     * \throw UTILS_ASSERT if number of points < 2
     *
     * \par Reference:
     * Catmull, E., & Rom, R. (1974). "A class of local interpolating splines."
     * In Computer Aided Geometric Design (pp. 317-326). Academic Press.
     */
    void catmull_rom();

    /**
     * \brief Build a spline using data in `GenericContainer`
     *
     * \param[in] gc GenericContainer containing spline data
     *
     * \par Expected GC Structure:
     * \code
     * {
     *   "data": matrix of points (rows: points, columns: dimensions, or transposed),
     *   "transposed": optional boolean (default: false)
     * }
     * \endcode
     *
     * \see setup(GenericContainer const &)
     */
    void build( GenericContainer const & gc ) { setup( gc ); }

    /**
     * \brief Build a spline using data in `GenericContainer`
     *
     * \param[in] gc GenericContainer containing spline configuration
     *
     * \par Data Format:
     * - The "data" key must contain a matrix of point coordinates
     * - If "transposed" is true: matrix is dim × npts (columns are points)
     * - If "transposed" is false or absent: matrix is npts × dim (rows are points)
     *
     * \throw UTILS_ASSERT if required data is missing or malformed
     */
    void setup( GenericContainer const & gc );

    ///@}

    /**
     * \brief Compute spline curvature at parameter x
     *
     * Computes the signed curvature κ(x) of the parametric curve at parameter x.
     * For 2D curves, this gives the signed curvature (positive for counter-clockwise
     * bending, negative for clockwise).
     *
     * \param[in] x Parameter value within the knot range
     * \return Curvature κ(x) = (x'·y'' - y'·x'') / (x'² + y'²)^{3/2}
     *
     * \par Mathematical Definition:
     * For a 2D parametric curve γ(t) = (x(t), y(t)):
     * \f[
     * \kappa(t) = \frac{x'(t)y''(t) - y'(t)x''(t)}{(x'(t)^2 + y'(t)^2)^{3/2}}
     * \f]
     *
     * \note This method assumes at least 2D spline. For higher dimensions,
     *       only the first two components are considered (projection onto xy-plane).
     *
     * \throw UTILS_ASSERT if spline dimension < 2
     *
     * \see curvature_D() for curvature derivative
     */
    real_type curvature( real_type x ) const;

    /**
     * \brief Compute derivative of curvature with respect to parameter x
     *
     * Computes dκ/dt, the rate of change of curvature along the curve.
     * Useful for analyzing smoothness and motion planning applications.
     *
     * \param[in] x Parameter value within the knot range
     * \return Curvature derivative dκ/dt at parameter x
     *
     * \par Mathematical Formula:
     * For a 2D parametric curve:
     * \f[
     * \frac{d\kappa}{dt} = \frac{(x'y''' - y'x''')(x'^2+y'^2) - 3(x'y'' - y'x'')(x'x'' + y'y'')}{(x'^2+y'^2)^{5/2}}
     * \f]
     *
     * \note This method assumes at least 2D spline. For higher dimensions,
     *       only the first two components are considered.
     *
     * \throw UTILS_ASSERT if spline dimension < 2
     *
     * \see curvature() for curvature value
     */
    real_type curvature_D( real_type x ) const;

    /**
     * \brief Get the spline type identifier
     *
     * \return SplineType1D::SPLINE_VEC enum value indicating this is a vector spline
     *
     * \par Usage:
     * This static method allows runtime type checking and dynamic dispatch
     * when working with spline hierarchies.
     */
    static SplineType1D type() { return SplineType1D::SPLINE_VEC; }

    /**
     * \brief Generate informational string about the spline
     *
     * Creates a human-readable string containing key spline properties:
     * name, number of points, and dimension.
     *
     * \return Formatted string: "SplineVec[name] n.points=N dim=D"
     *
     * \par Example Output:
     * "SplineVec[my_curve] n.points=10 dim=3"
     */
    string info() const { return fmt::format( "SplineVec[{}] n.points={}  dim={}", name(), m_npts, m_dim ); }

    /**
     * \brief Print spline information to output stream
     *
     * Writes the informational string (from info()) followed by newline
     * to the specified output stream.
     *
     * \param[out] stream Output stream (e.g., std::cout, std::ofstream)
     *
     * \par Usage:
     * Useful for debugging and logging spline properties.
     */
    void info( ostream_type & stream ) const { stream << this->info() << '\n'; }

    /**
     * \brief Generate table of spline values for plotting
     *
     * Samples the spline uniformly over its parameter domain and writes
     * a tab-separated table to the output stream. The first column contains
     * parameter values, subsequent columns contain component values.
     *
     * \param[out] stream     Output stream for the table
     * \param[in]  num_points Number of sample points (including endpoints)
     *
     * \par Output Format:
     * ```
     * s    0    1    2    ...  # Header line (dim columns)
     * 0.0  v00  v01  v02  ...  # First sample
     * 0.1  v10  v11  v12  ...  # Second sample
     * ...  ...  ...  ...  ...  # More samples
     * 1.0  vn0  vn1  vn2  ...  # Last sample
     * ```
     *
     * \par Use Cases:
     * - Exporting for visualization tools (Gnuplot, MATLAB, Excel)
     * - Debugging spline shape and properties
     *
     * \note Uses uniform sampling which may not capture all features for
     *       non-uniform knot distributions.
     */
    void dump_table( ostream_type & stream, integer const num_points ) const;

#ifdef SPLINES_BACK_COMPATIBILITY
    /**
     * \brief Backward compatibility: get number of points
     * \return Number of interpolation points
     */
    integer numPoints() const { return m_npts; }

    /**
     * \brief Backward compatibility: get pointer to knot vector
     * \return Pointer to array of knot values
     */
    real_type const * xNodes() const { return m_X; }

    /**
     * \brief Backward compatibility: get specific knot value
     * \param[in] npt Knot index
     * \return Knot value at index npt
     */
    real_type xNode( integer npt ) const { return m_X[npt]; }

    /**
     * \brief Backward compatibility: get pointer to point coordinates
     * \param[in] j Dimension index
     * \return Pointer to coordinate values for dimension j
     */
    real_type const * yNodes( integer j ) const { return m_Y[j]; }

    /**
     * \brief Backward compatibility: get specific point coordinate
     * \param[in] npt Point index
     * \param[in] j   Dimension index
     * \return Coordinate value Y_j[npt]
     */
    real_type yNode( integer npt, integer j ) const { return y_node( npt, j ); }

    /**
     * \brief Backward compatibility: get minimum parameter
     * \return Minimum knot value
     */
    real_type xMin() const { return m_X[0]; }

    /**
     * \brief Backward compatibility: get maximum parameter
     * \return Maximum knot value
     */
    real_type xMax() const { return m_X[m_npts - 1]; }

    /**
     * \brief Backward compatibility: set knot values
     * \param[in] X Array of knot values
     */
    void setKnots( real_type const X[] ) { this->set_knots( X ); }

    /**
     * \brief Backward compatibility: chord length parametrization
     */
    void setKnotsChordLength() { this->set_knots_chord_length(); }

    /**
     * \brief Backward compatibility: centripetal parametrization
     */
    void setKnotsCentripetal() { this->set_knots_centripetal(); }

    /**
     * \brief Backward compatibility: Foley parametrization
     */
    void setKnotsFoley() { this->set_knots_foley(); }

    /**
     * \brief Backward compatibility: Catmull-Rom derivative computation
     */
    void CatmullRom() { this->catmull_rom(); }
#endif
  };

}  // namespace Splines

#endif

// EOF: SplineVec.hxx
