
/**
 * \brief Evaluate all splines using a specific spline as independent variable
 *
 * Instead of evaluating at a given x, evaluates at the x where spline `spl`
 * equals `zeta`. This is useful for inverse evaluation or cross-section analysis.
 *
 * \param[in] spl Index of spline to use as independent variable (must be monotone)
 * \param[in] zeta Target value of the independent spline
 * \param[out] x Computed x-value such that spline[spl](x) = zeta
 * \param[out] vals Vector to store all spline values at computed x
 *
 * \throw UTILS_ASSERT if spline `spl` is not monotone
 * \throw UTILS_ASSERT if zeta is outside range of spline `spl`
 */
void eval2( integer const spl, real_type const zeta, real_type & x, vector<real_type> & vals ) const
{
  vals.resize( m_nspl );
  this->eval2( spl, zeta, x, vals.data(), 1 );
}

/**
 * \brief Evaluate first derivatives using a specific spline as independent variable
 *
 * \param[in] spl Index of independent spline
 * \param[in] zeta Target value of independent spline
 * \param[out] x Computed x-value
 * \param[out] vals Vector to store first derivatives at computed x
 */
void eval2_D( integer const spl, real_type const zeta, real_type & x, vector<real_type> & vals ) const
{
  vals.resize( m_nspl );
  this->eval2_D( spl, zeta, x, vals.data(), 1 );
}

/**
 * \brief Evaluate second derivatives using a specific spline as independent variable
 *
 * \param[in] spl Index of independent spline
 * \param[in] zeta Target value of independent spline
 * \param[out] x Computed x-value
 * \param[out] vals Vector to store second derivatives at computed x
 */
void eval2_DD( integer const spl, real_type const zeta, real_type & x, vector<real_type> & vals ) const
{
  vals.resize( m_nspl );
  this->eval2_DD( spl, zeta, x, vals.data(), 1 );
}

/**
 * \brief Evaluate third derivatives using a specific spline as independent variable
 *
 * \param[in] spl Index of independent spline
 * \param[in] zeta Target value of independent spline
 * \param[out] x Computed x-value
 * \param[out] vals Vector to store third derivatives at computed x
 */
void eval2_DDD( integer const spl, real_type const zeta, real_type & x, vector<real_type> & vals ) const
{
  vals.resize( m_nspl );
  this->eval2_DDD( spl, zeta, x, vals.data(), 1 );
}

/**
 * \brief Evaluate all splines using a specific spline as independent variable
 *
 * Instead of evaluating at a given x, evaluates at the x where spline `spl`
 * equals `zeta`. This is useful for inverse evaluation or cross-section analysis.
 *
 * \param[in] spl Index of spline to use as independent variable (must be monotone)
 * \param[in] zeta Target value of the independent spline
 * \param[out] vals Vector to store all spline values at computed x
 *
 * \throw UTILS_ASSERT if spline `spl` is not monotone
 * \throw UTILS_ASSERT if zeta is outside range of spline `spl`
 */
void eval2( integer const spl, real_type const zeta, vector<real_type> & vals ) const
{
  real_type x;
  return this->eval2( spl, zeta, x, vals );
}

/**
 * \brief Evaluate first derivatives using a specific spline as independent variable
 *
 * \param[in] spl Index of independent spline
 * \param[in] zeta Target value of independent spline
 * \param[out] vals Vector to store first derivatives at computed x
 */
void eval2_D( integer const spl, real_type const zeta, vector<real_type> & vals ) const
{
  real_type x;
  return this->eval2_D( spl, zeta, x, vals );
}

/**
 * \brief Evaluate second derivatives using a specific spline as independent variable
 *
 * \param[in] spl Index of independent spline
 * \param[in] zeta Target value of independent spline
 * \param[out] vals Vector to store second derivatives at computed x
 */
void eval2_DD( integer const spl, real_type const zeta, vector<real_type> & vals ) const
{
  real_type x;
  return this->eval2_DD( spl, zeta, x, vals );
}

/**
 * \brief Evaluate third derivatives using a specific spline as independent variable
 *
 * \param[in] spl Index of independent spline
 * \param[in] zeta Target value of independent spline
 * \param[out] vals Vector to store third derivatives at computed x
 */
void eval2_DDD( integer const spl, real_type const zeta, vector<real_type> & vals ) const
{
  real_type x;
  return this->eval2_DDD( spl, zeta, x, vals );
}
