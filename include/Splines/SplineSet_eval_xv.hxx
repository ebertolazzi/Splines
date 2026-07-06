/**
 * \brief Evaluate all splines at given x and store in C array with stride
 *
 * \param[in] x Parameter value
 * \param[out] vals Output array (must have space for at least m_nspl*incy elements)
 * \param[in] incy Stride between consecutive values in output array (default: 1)
 *
 * \par Memory Layout:
 * - incy = 1: vals[0..m_nspl-1] = spline values
 * - incy = k: vals[0, k, 2k, ...] = spline values
 */
void eval( real_type const x, real_type vals[], integer const incy = 1 ) const
{
  integer ii = 0;
  for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->eval( x );
}

/**
 * \brief Evaluate first derivatives of all splines to C array with stride
 *
 * \param[in] x Parameter value
 * \param[out] vals Output array for derivatives
 * \param[in] incy Stride between consecutive values (default: 1)
 */
void eval_D( real_type const x, real_type vals[], integer const incy = 1 ) const
{
  size_t ii = 0;
  for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->D( x );
}

/**
 * \brief Evaluate second derivatives of all splines to C array with stride
 *
 * \param[in] x Parameter value
 * \param[out] vals Output array for second derivatives
 * \param[in] incy Stride between consecutive values (default: 1)
 */
void eval_DD( real_type const x, real_type vals[], integer const incy = 1 ) const
{
  size_t ii = 0;
  for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->DD( x );
}

/**
 * \brief Evaluate third derivatives of all splines to C array with stride
 *
 * \param[in] x Parameter value
 * \param[out] vals Output array for third derivatives
 * \param[in] incy Stride between consecutive values (default: 1)
 */
void eval_DDD( real_type const x, real_type vals[], integer const incy = 1 ) const
{
  size_t ii = 0;
  for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->DDD( x );
}

/**
 * \brief Evaluate fourth derivatives of all splines to C array with stride
 *
 * \param[in] x Parameter value
 * \param[out] vals Output array for fourth derivatives
 * \param[in] incy Stride between consecutive values (default: 1)
 */
void eval_DDDD( real_type const x, real_type vals[], integer const incy = 1 ) const
{
  size_t ii = 0;
  for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->DDDD( x );
}

/**
 * \brief Evaluate fifth derivatives of all splines to C array with stride
 *
 * \param[in] x Parameter value
 * \param[out] vals Output array for fifth derivatives
 * \param[in] incy Stride between consecutive values (default: 1)
 */
void eval_DDDDD( real_type const x, real_type vals[], integer const incy = 1 ) const
{
  size_t ii = 0;
  for ( integer i = 0; i < m_nspl; ++i, ii += incy ) vals[ii] = m_splines[i]->DDDDD( x );
}
