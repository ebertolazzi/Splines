/**
 * \brief Evaluate all splines at given x and store in std::vector
 *
 * \param[in] x Parameter value
 * \param[out] vals Vector to store results (automatically resized to m_nspl)
 *
 * \par Post-conditions:
 * - vals.size() == m_nspl
 * - vals[i] == value of i-th spline at x
 */
void eval( real_type const x, vector<real_type> & vals ) const
{
  vals.resize( m_nspl );
  for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->eval( x );
}

/**
 * \brief Evaluate first derivatives of all splines and store in std::vector
 *
 * \param[in] x Parameter value
 * \param[out] vals Vector to store derivatives (automatically resized)
 */
void eval_D( real_type const x, vector<real_type> & vals ) const
{
  vals.resize( m_nspl );
  for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->D( x );
}

/**
 * \brief Evaluate second derivatives of all splines and store in std::vector
 *
 * \param[in] x Parameter value
 * \param[out] vals Vector to store second derivatives (automatically resized)
 */
void eval_DD( real_type const x, vector<real_type> & vals ) const
{
  vals.resize( m_nspl );
  for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->DD( x );
}

/**
 * \brief Evaluate third derivatives of all splines and store in std::vector
 *
 * \param[in] x Parameter value
 * \param[out] vals Vector to store third derivatives (automatically resized)
 */
void eval_DDD( real_type const x, vector<real_type> & vals ) const
{
  vals.resize( m_nspl );
  for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->DDD( x );
}

/**
 * \brief Evaluate fourth derivatives of all splines and store in std::vector
 *
 * \param[in] x Parameter value
 * \param[out] vals Vector to store fourth derivatives (automatically resized)
 */
void eval_DDDD( real_type const x, vector<real_type> & vals ) const
{
  vals.resize( m_nspl );
  for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->DDDD( x );
}

/**
 * \brief Evaluate fifth derivatives of all splines and store in std::vector
 *
 * \param[in] x Parameter value
 * \param[out] vals Vector to store fifth derivatives (automatically resized)
 */
void eval_DDDDD( real_type const x, vector<real_type> & vals ) const
{
  vals.resize( m_nspl );
  for ( integer i = 0; i < m_nspl; ++i ) vals[i] = m_splines[i]->DDDDD( x );
}
