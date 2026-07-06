/**
 * \brief Evaluate all splines at x and store results in GenericContainer
 *
 * Evaluates all splines at the given parameter value and stores results
 * in a map where keys are spline names and values are the evaluated results.
 *
 * \param[in] x Parameter value
 * \param[out] gc GenericContainer to store results (map: spline_name -> value)
 */
void eval( real_type const x, GenericContainer & gc ) const;

/**
 * \brief Evaluate first derivatives of all splines at x and store in GenericContainer
 *
 * Evaluates the first derivative dy/dx for all splines at the given parameter value.
 *
 * \param[in] x Parameter value
 * \param[out] gc GenericContainer with derivative values (map: spline_name -> dy/dx)
 */
void eval_D( real_type const x, GenericContainer & gc ) const;

/**
 * \brief Evaluate second derivatives of all splines at x and store in GenericContainer
 *
 * Evaluates the second derivative d²y/dx² for all splines at the given parameter value.
 *
 * \param[in] x Parameter value
 * \param[out] gc GenericContainer with second derivative values (map: spline_name -> d²y/dx²)
 */
void eval_DD( real_type const x, GenericContainer & gc ) const;

/**
 * \brief Evaluate third derivatives of all splines at x and store in GenericContainer
 *
 * Evaluates the third derivative d³y/dx³ for all splines at the given parameter value.
 *
 * \param[in] x Parameter value
 * \param[out] gc GenericContainer with third derivative values (map: spline_name -> d³y/dx³)
 */
void eval_DDD( real_type const x, GenericContainer & gc ) const;

/**
 * \brief Evaluate all splines at multiple x-values and store in GenericContainer
 *
 * For each spline, creates a vector of values corresponding to each x in the input vector.
 * Results are stored as a map where keys are spline names and values are vectors.
 *
 * \param[in] vec Vector of x-values at which to evaluate
 * \param[out] gc GenericContainer with results (map: spline_name -> vector of values)
 */
void eval( vec_real_type const & vec, GenericContainer & gc ) const;

/**
 * \brief Evaluate first derivatives of all splines at multiple x-values
 *
 * For each spline, creates a vector of first derivatives dy/dx corresponding
 * to each x in the input vector.
 *
 * \param[in] vec Vector of x-values at which to evaluate
 * \param[out] gc GenericContainer with vectors of derivatives (map: spline_name -> vector of dy/dx)
 */
void eval_D( vec_real_type const & vec, GenericContainer & gc ) const;

/**
 * \brief Evaluate second derivatives of all splines at multiple x-values
 *
 * For each spline, creates a vector of second derivatives d²y/dx² corresponding
 * to each x in the input vector.
 *
 * \param[in] vec Vector of x-values at which to evaluate
 * \param[out] gc GenericContainer with vectors of second derivatives (map: spline_name -> vector of d²y/dx²)
 */
void eval_DD( vec_real_type const & vec, GenericContainer & gc ) const;

/**
 * \brief Evaluate third derivatives of all splines at multiple x-values
 *
 * For each spline, creates a vector of third derivatives d³y/dx³ corresponding
 * to each x in the input vector.
 *
 * \param[in] vec Vector of x-values at which to evaluate
 * \param[out] gc GenericContainer with vectors of third derivatives (map: spline_name -> vector of d³y/dx³)
 */
void eval_DDD( vec_real_type const & vec, GenericContainer & gc ) const;

/**
 * \brief Evaluate specific splines at x and store in GenericContainer
 *
 * Evaluates only the splines specified in the columns parameter.
 *
 * \param[in] x Parameter value
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with requested spline values (map: spline_name -> value)
 */
void eval( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const;

/**
 * \brief Evaluate first derivatives of specific splines at x
 *
 * Evaluates first derivatives dy/dx only for the splines specified in the columns parameter.
 *
 * \param[in] x Parameter value
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with derivative values (map: spline_name -> dy/dx)
 */
void eval_D( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const;

/**
 * \brief Evaluate second derivatives of specific splines at x
 *
 * Evaluates second derivatives d²y/dx² only for the splines specified in the columns parameter.
 *
 * \param[in] x Parameter value
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with second derivative values (map: spline_name -> d²y/dx²)
 */
void eval_DD( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const;

/**
 * \brief Evaluate third derivatives of specific splines at x
 *
 * Evaluates third derivatives d³y/dx³ only for the splines specified in the columns parameter.
 *
 * \param[in] x Parameter value
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with third derivative values (map: spline_name -> d³y/dx³)
 */
void eval_DDD( real_type const x, vec_string_type const & columns, GenericContainer & gc ) const;

/**
 * \brief Evaluate specific splines at multiple x-values and store in GenericContainer (OPTIMIZED)
 *
 * For each specified spline, creates a vector of values corresponding to each x in the input vector.
 *
 * \note This version pre-caches spline pointers to avoid repeated lookups in the inner loop.
 *
 * \param[in] vec Vector of x-values at which to evaluate
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of values (map: spline_name -> vector of values)
 */
void eval( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const;

/**
 * \brief Evaluate first derivatives of specific splines at multiple x-values (OPTIMIZED)
 *
 * For each specified spline, creates a vector of first derivatives dy/dx corresponding
 * to each x in the input vector.
 *
 * \note This version pre-caches spline pointers to avoid repeated lookups in the inner loop.
 *
 * \param[in] vec Vector of x-values at which to evaluate
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of derivatives (map: spline_name -> vector of dy/dx)
 */
void eval_D( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const;

/**
 * \brief Evaluate second derivatives of specific splines at multiple x-values (OPTIMIZED)
 *
 * For each specified spline, creates a vector of second derivatives d²y/dx² corresponding
 * to each x in the input vector.
 *
 * \note This version pre-caches spline pointers to avoid repeated lookups in the inner loop.
 *
 * \param[in] vec Vector of x-values at which to evaluate
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of second derivatives (map: spline_name -> vector of d²y/dx²)
 */
void eval_DD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const;

/**
 * \brief Evaluate third derivatives of specific splines at multiple x-values (OPTIMIZED)
 *
 * For each specified spline, creates a vector of third derivatives d³y/dx³ corresponding
 * to each x in the input vector.
 *
 * \note This version pre-caches spline pointers to avoid repeated lookups in the inner loop.
 *
 * \param[in] vec Vector of x-values at which to evaluate
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of third derivatives (map: spline_name -> vector of d³y/dx³)
 */
void eval_DDD( vec_real_type const & vec, vec_string_type const & columns, GenericContainer & gc ) const;
