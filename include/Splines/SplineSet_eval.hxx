/**
 * \brief Evaluate a spline by index at given x
 *
 * \param[in] x Parameter value (in the domain of x-nodes)
 * \param[in] spl Index of the spline to evaluate
 * \return Value of spline `spl` at parameter `x`
 *
 * \par Complexity:
 * - Time: O(log npts) for interval search plus spline evaluation
 * - Space: O(1)
 */
real_type operator()( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->eval( x );
}

/**
 * \brief Evaluate a spline by index at given x
 *
 * \copydetails operator()()
 * This is an alias for operator()().
 */
real_type eval( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->eval( x );
}

/**
 * \brief Evaluate first derivative of a spline by index
 *
 * \param[in] x Parameter value
 * \param[in] spl Index of the spline
 * \return First derivative of spline `spl` at `x`
 */
real_type D( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->D( x );
}

/**
 * \brief Evaluate first derivative of a spline by index
 *
 * \copydetails D()
 * This is an alias for D().
 */
real_type eval_D( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->D( x );
}

/**
 * \brief Evaluate second derivative of a spline by index
 *
 * \param[in] x Parameter value
 * \param[in] spl Index of the spline
 * \return Second derivative of spline `spl` at `x`
 */
real_type DD( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->DD( x );
}

/**
 * \brief Evaluate second derivative of a spline by index
 *
 * \copydetails DD()
 * This is an alias for DD().
 */
real_type eval_DD( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->DD( x );
}

/**
 * \brief Evaluate third derivative of a spline by index
 *
 * \param[in] x Parameter value
 * \param[in] spl Index of the spline
 * \return Third derivative of spline `spl` at `x`
 */
real_type DDD( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->DDD( x );
}

/**
 * \brief Evaluate third derivative of a spline by index
 *
 * \copydetails DDD()
 * This is an alias for DDD().
 */
real_type eval_DDD( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->DDD( x );
}

/**
 * \brief Evaluate fourth derivative of a spline by index
 *
 * \param[in] x Parameter value
 * \param[in] spl Index of the spline
 * \return Fourth derivative of spline `spl` at `x`
 *
 * \note For cubic splines, returns zero beyond third derivative.
 */
real_type DDDD( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->DDDD( x );
}

/**
 * \brief Evaluate fourth derivative of a spline by index
 *
 * \copydetails DDDD()
 * This is an alias for DDDD().
 */
real_type eval_DDDD( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->DDDD( x );
}

/**
 * \brief Evaluate fifth derivative of a spline by index
 *
 * \param[in] x Parameter value
 * \param[in] spl Index of the spline
 * \return Fifth derivative of spline `spl` at `x`
 *
 * \note For cubic splines, returns zero beyond third derivative.
 */
real_type DDDDD( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->DDDDD( x );
}

/**
 * \brief Evaluate fifth derivative of a spline by index
 *
 * \copydetails DDDDD()
 * This is an alias for DDDDD().
 */
real_type eval_DDDDD( real_type const x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->DDDDD( x );
}

#ifdef AUTODIFF_SUPPORT
//!
//! \name Autodiff
//!
///@{
/**
 * \brief Evaluate spline with first-order automatic differentiation
 *
 * \param[in] x Dual number parameter (contains value and gradient)
 * \param[in] spl Index of the spline
 * \return Dual number containing function value and derivative
 *
 * \par Usage:
 * \code{.cpp}
 * autodiff::dual1st t = 1.5;
 * autodiff::dual1st y = spline_set.eval(t, 0);
 * // y.val = value at t=1.5
 * // y.grad = derivative dY/dt at t=1.5
 * \endcode
 */
autodiff::dual1st eval( autodiff::dual1st const & x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->eval( x );
}

/**
 * \brief Evaluate spline with second-order automatic differentiation
 *
 * \param[in] x Dual number parameter (second order)
 * \param[in] spl Index of the spline
 * \return Dual number containing function value, first and second derivatives
 */
autodiff::dual2nd eval( autodiff::dual2nd const & x, integer const spl ) const
{
  Spline const * S = this->get_spline( spl );
  return S->eval( x );
}

/**
 * \brief Generic template for automatic differentiation
 *
 * Automatically selects the appropriate evaluation method based on input type.
 *
 * \tparam T Input type (arithmetic or autodiff dual type)
 * \param[in] x Parameter value
 * \param[in] spl Index of the spline
 * \return Result matching input type
 *
 * \par Type Detection:
 * - Arithmetic types (int, float, double): promotes to real_type
 * - Autodiff dual types: uses appropriate autodiff evaluation
 */
template <typename T> auto eval( T const & x, integer const spl ) const
{
  if constexpr ( std::is_arithmetic<T>::value )
  {
    // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
    return eval( static_cast<real_type>( x ), spl );
  }
  else
  {
    // Altrimenti deduce automaticamente il tipo duale appropriato
    return eval( autodiff::detail::to_dual( x ), spl );
  }
}

/**
 * \brief Generic function call operator for autodiff
 *
 * \tparam T Parameter type
 * \param[in] x Parameter value
 * \param[in] spl Index of the spline
 * \return Result matching input type
 */
template <typename T> auto operator()( T const & x, integer const spl ) const -> decltype( eval( x, spl ) )
{ return eval( x, spl ); }

#endif
