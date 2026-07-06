/**
 * \brief Evaluate all splines using independent spline and store in GenericContainer
 *
 * This method finds the x-value where the independent spline equals zeta,
 * then evaluates all splines at that x-value.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] x Computed x-value from intersection
 * \param[out] gc GenericContainer to store results (map: spline_name -> value)
 */
void eval2( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval( x );
}

/**
 * \brief Evaluate first derivatives using independent spline
 *
 * This method finds the x-value where the independent spline equals zeta,
 * then evaluates the first derivative of all splines at that x-value.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] x Computed x-value from intersection
 * \param[out] gc GenericContainer with derivative values (map: spline_name -> derivative)
 */
void eval2_D( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_D( x );
}

/**
 * \brief Evaluate second derivatives using independent spline
 *
 * This method finds the x-value where the independent spline equals zeta,
 * then evaluates the second derivative of all splines at that x-value.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] x Computed x-value from intersection
 * \param[out] gc GenericContainer with second derivative values (map: spline_name -> d²y/dx²)
 */
void eval2_DD( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DD( x );
}

/**
 * \brief Evaluate third derivatives using independent spline
 *
 * This method finds the x-value where the independent spline equals zeta,
 * then evaluates the third derivative of all splines at that x-value.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] x Computed x-value from intersection
 * \param[out] gc GenericContainer with third derivative values (map: spline_name -> d³y/dx³)
 */
void eval2_DDD( real_type const zeta, integer const indep, real_type & x, GenericContainer & gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst] = m_splines[snd]->eval_DDD( x );
}

/**
 * \brief Evaluate all splines at multiple zeta values using independent spline
 *
 * For each zeta value, finds the corresponding x-value via intersection,
 * then evaluates all splines at those x-values.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] gc GenericContainer with results (map: spline_name -> vector of values)
 */
void eval2( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // Pre-allocation
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v = vals[fst].get_vec_real();
      v[i]              = m_splines[snd]->eval( x );
    }
  }
}

/**
 * \brief Evaluate first derivatives at multiple zeta values using independent spline
 *
 * For each zeta value, finds the corresponding x-value via intersection,
 * then evaluates the first derivative of all splines at those x-values.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] gc GenericContainer with vectors of derivatives (map: spline_name -> vector of dy/dx)
 */
void eval2_D( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // Pre-allocation
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v = vals[fst].get_vec_real();
      v[i]              = m_splines[snd]->eval_D( x );
    }
  }
}

/**
 * \brief Evaluate second derivatives at multiple zeta values using independent spline
 *
 * For each zeta value, finds the corresponding x-value via intersection,
 * then evaluates the second derivative of all splines at those x-values.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] gc GenericContainer with vectors of second derivatives (map: spline_name -> vector of d²y/dx²)
 */
void eval2_DD( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // Pre-allocation
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v = vals[fst].get_vec_real();
      v[i]              = m_splines[snd]->eval_DD( x );
    }
  }
}

/**
 * \brief Evaluate third derivatives at multiple zeta values using independent spline
 *
 * For each zeta value, finds the corresponding x-value via intersection,
 * then evaluates the third derivative of all splines at those x-values.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] gc GenericContainer with vectors of third derivatives (map: spline_name -> vector of d³y/dx³)
 */
void eval2_DDD( vec_real_type const & zetas, integer const indep, GenericContainer & gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // Pre-allocation
  for ( auto const & [fst, snd] : m_header_to_position ) vals[fst].set_vec_real( npts );

  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );
    for ( auto const & [fst, snd] : m_header_to_position )
    {
      vec_real_type & v = vals[fst].get_vec_real();
      v[i]              = m_splines[snd]->eval_DDD( x );
    }
  }
}

/**
 * \brief Evaluate specific splines using independent spline and store in GenericContainer
 *
 * Evaluates only the splines specified in the columns parameter.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] x Computed x-value from intersection
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer to store results (map: spline_name -> value)
 */
void eval2(
  real_type const         zeta,
  integer const           indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval( x );
  }
}

/**
 * \brief Evaluate first derivatives of specific splines using independent spline
 *
 * Evaluates first derivatives only for the splines specified in the columns parameter.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] x Computed x-value from intersection
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with derivative values (map: spline_name -> dy/dx)
 */
void eval2_D(
  real_type const         zeta,
  integer const           indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval_D( x );
  }
}

/**
 * \brief Evaluate second derivatives of specific splines using independent spline
 *
 * Evaluates second derivatives only for the splines specified in the columns parameter.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] x Computed x-value from intersection
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with second derivative values (map: spline_name -> d²y/dx²)
 */
void eval2_DD(
  real_type const         zeta,
  integer const           indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval_DD( x );
  }
}

/**
 * \brief Evaluate third derivatives of specific splines using independent spline
 *
 * Evaluates third derivatives only for the splines specified in the columns parameter.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[out] x Computed x-value from intersection
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with third derivative values (map: spline_name -> d³y/dx³)
 */
void eval2_DDD(
  real_type const         zeta,
  integer const           indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  map_type & vals = gc.set_map();
  intersect( indep, zeta, x );
  for ( auto const & S : columns )
  {
    Spline const * p_spl = get_spline( S );
    vals[S]              = p_spl->eval_DDD( x );
  }
}

/**
 * \brief Evaluate specific splines at multiple zeta values using independent spline (OPTIMIZED)
 *
 * For each zeta value, finds the corresponding x-value via intersection,
 * then evaluates the specified splines at those x-values.
 *
 * \note This version pre-caches spline pointers to avoid repeated lookups in the inner loop.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with results (map: spline_name -> vector of values)
 */
void eval2( vec_real_type const & zetas, integer const indep, vec_string_type const & columns, GenericContainer & gc )
  const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // Pre-allocation and cache spline pointers
  std::vector<Spline const *> spline_ptrs;
  spline_ptrs.reserve( columns.size() );

  for ( auto const & S : columns )
  {
    vals[S].set_vec_real( npts );
    spline_ptrs.push_back( get_spline( S ) );
  }

  // Evaluate at all points
  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );

    for ( size_t j = 0; j < columns.size(); ++j )
    {
      vec_real_type & v = vals[columns[j]].get_vec_real();
      v[i]              = spline_ptrs[j]->eval( x );
    }
  }
}

/**
 * \brief Evaluate first derivatives of specific splines at multiple zeta values (OPTIMIZED)
 *
 * For each zeta value, finds the corresponding x-value via intersection,
 * then evaluates the first derivative of the specified splines at those x-values.
 *
 * \note This version pre-caches spline pointers to avoid repeated lookups in the inner loop.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of derivatives (map: spline_name -> vector of dy/dx)
 */
void eval2_D( vec_real_type const & zetas, integer const indep, vec_string_type const & columns, GenericContainer & gc )
  const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // Pre-allocation and cache spline pointers
  std::vector<Spline const *> spline_ptrs;
  spline_ptrs.reserve( columns.size() );

  for ( auto const & S : columns )
  {
    vals[S].set_vec_real( npts );
    spline_ptrs.push_back( get_spline( S ) );
  }

  // Evaluate at all points
  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );

    for ( size_t j = 0; j < columns.size(); ++j )
    {
      vec_real_type & v = vals[columns[j]].get_vec_real();
      v[i]              = spline_ptrs[j]->eval_D( x );
    }
  }
}

/**
 * \brief Evaluate second derivatives of specific splines at multiple zeta values (OPTIMIZED)
 *
 * For each zeta value, finds the corresponding x-value via intersection,
 * then evaluates the second derivative of the specified splines at those x-values.
 *
 * \note This version pre-caches spline pointers to avoid repeated lookups in the inner loop.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of second derivatives (map: spline_name -> vector of d²y/dx²)
 */
void eval2_DD(
  vec_real_type const &   zetas,
  integer const           indep,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // Pre-allocation and cache spline pointers
  std::vector<Spline const *> spline_ptrs;
  spline_ptrs.reserve( columns.size() );

  for ( auto const & S : columns )
  {
    vals[S].set_vec_real( npts );
    spline_ptrs.push_back( get_spline( S ) );
  }

  // Evaluate at all points
  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );

    for ( size_t j = 0; j < columns.size(); ++j )
    {
      vec_real_type & v = vals[columns[j]].get_vec_real();
      v[i]              = spline_ptrs[j]->eval_DD( x );
    }
  }
}

/**
 * \brief Evaluate third derivatives of specific splines at multiple zeta values (OPTIMIZED)
 *
 * For each zeta value, finds the corresponding x-value via intersection,
 * then evaluates the third derivative of the specified splines at those x-values.
 *
 * \note This version pre-caches spline pointers to avoid repeated lookups in the inner loop.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Index of independent spline (0-based)
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of third derivatives (map: spline_name -> vector of d³y/dx³)
 */
void eval2_DDD(
  vec_real_type const &   zetas,
  integer const           indep,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{
  integer const npts = static_cast<integer>( zetas.size() );
  map_type &    vals = gc.set_map();

  // Pre-allocation and cache spline pointers
  std::vector<Spline const *> spline_ptrs;
  spline_ptrs.reserve( columns.size() );

  for ( auto const & S : columns )
  {
    vals[S].set_vec_real( npts );
    spline_ptrs.push_back( get_spline( S ) );
  }

  // Evaluate at all points
  for ( integer i = 0; i < npts; ++i )
  {
    real_type x;
    intersect( indep, zetas[i], x );

    for ( size_t j = 0; j < columns.size(); ++j )
    {
      vec_real_type & v = vals[columns[j]].get_vec_real();
      v[i]              = spline_ptrs[j]->eval_DDD( x );
    }
  }
}

// ============================================================================
// String-based interface (delegates to index-based versions)
// ============================================================================

/**
 * \brief Evaluate all splines using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value from intersection
 * \param[out] gc GenericContainer to store results (map: spline_name -> value)
 */
void eval2( real_type const zeta, string_view indep, real_type & x, GenericContainer & gc ) const
{ this->eval2( zeta, this->get_position( indep ), x, gc ); }

/**
 * \brief Evaluate first derivatives using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value from intersection
 * \param[out] gc GenericContainer with derivative values (map: spline_name -> dy/dx)
 */
void eval2_D( real_type const zeta, string_view indep, real_type & x, GenericContainer & gc ) const
{ this->eval2_D( zeta, this->get_position( indep ), x, gc ); }

/**
 * \brief Evaluate second derivatives using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value from intersection
 * \param[out] gc GenericContainer with second derivative values (map: spline_name -> d²y/dx²)
 */
void eval2_DD( real_type const zeta, string_view indep, real_type & x, GenericContainer & gc ) const
{ this->eval2_DD( zeta, this->get_position( indep ), x, gc ); }

/**
 * \brief Evaluate third derivatives using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value from intersection
 * \param[out] gc GenericContainer with third derivative values (map: spline_name -> d³y/dx³)
 */
void eval2_DDD( real_type const zeta, string_view indep, real_type & x, GenericContainer & gc ) const
{ this->eval2_DDD( zeta, this->get_position( indep ), x, gc ); }

/**
 * \brief Evaluate all splines at multiple zeta values using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[out] gc GenericContainer with results (map: spline_name -> vector of values)
 */
void eval2( vec_real_type const & zetas, string_view indep, GenericContainer & gc ) const
{ this->eval2( zetas, this->get_position( indep ), gc ); }

/**
 * \brief Evaluate first derivatives at multiple zeta values using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[out] gc GenericContainer with vectors of derivatives (map: spline_name -> vector of dy/dx)
 */
void eval2_D( vec_real_type const & zetas, string_view indep, GenericContainer & gc ) const
{ this->eval2_D( zetas, this->get_position( indep ), gc ); }

/**
 * \brief Evaluate second derivatives at multiple zeta values using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[out] gc GenericContainer with vectors of second derivatives (map: spline_name -> vector of d²y/dx²)
 */
void eval2_DD( vec_real_type const & zetas, string_view indep, GenericContainer & gc ) const
{ this->eval2_DD( zetas, this->get_position( indep ), gc ); }

/**
 * \brief Evaluate third derivatives at multiple zeta values using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[out] gc GenericContainer with vectors of third derivatives (map: spline_name -> vector of d³y/dx³)
 */
void eval2_DDD( vec_real_type const & zetas, string_view indep, GenericContainer & gc ) const
{ this->eval2_DDD( zetas, this->get_position( indep ), gc ); }

/**
 * \brief Evaluate specific splines using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value from intersection
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer to store results (map: spline_name -> value)
 */
void eval2(
  real_type const         zeta,
  string_view             indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{ this->eval2( zeta, this->get_position( indep ), x, columns, gc ); }

/**
 * \brief Evaluate first derivatives of specific splines using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value from intersection
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with derivative values (map: spline_name -> dy/dx)
 */
void eval2_D(
  real_type const         zeta,
  string_view             indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{ this->eval2_D( zeta, this->get_position( indep ), x, columns, gc ); }

/**
 * \brief Evaluate second derivatives of specific splines using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value from intersection
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with second derivative values (map: spline_name -> d²y/dx²)
 */
void eval2_DD(
  real_type const         zeta,
  string_view             indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{ this->eval2_DD( zeta, this->get_position( indep ), x, columns, gc ); }

/**
 * \brief Evaluate third derivatives of specific splines using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value from intersection
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with third derivative values (map: spline_name -> d³y/dx³)
 */
void eval2_DDD(
  real_type const         zeta,
  string_view             indep,
  real_type &             x,
  vec_string_type const & columns,
  GenericContainer &      gc ) const
{ this->eval2_DDD( zeta, this->get_position( indep ), x, columns, gc ); }

/**
 * \brief Evaluate specific splines at multiple zeta values using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with results (map: spline_name -> vector of values)
 */
void eval2( vec_real_type const & zetas, string_view indep, vec_string_type const & columns, GenericContainer & gc )
  const
{ this->eval2( zetas, this->get_position( indep ), columns, gc ); }

/**
 * \brief Evaluate first derivatives of specific splines at multiple zeta values using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of derivatives (map: spline_name -> vector of dy/dx)
 */
void eval2_D( vec_real_type const & zetas, string_view indep, vec_string_type const & columns, GenericContainer & gc )
  const
{ this->eval2_D( zetas, this->get_position( indep ), columns, gc ); }

/**
 * \brief Evaluate second derivatives of specific splines at multiple zeta values using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of second derivatives (map: spline_name -> vector of d²y/dx²)
 */
void eval2_DD( vec_real_type const & zetas, string_view indep, vec_string_type const & columns, GenericContainer & gc )
  const
{ this->eval2_DD( zetas, this->get_position( indep ), columns, gc ); }

/**
 * \brief Evaluate third derivatives of specific splines at multiple zeta values using independent spline by name
 *
 * Convenience wrapper that accepts the independent spline name instead of index.
 *
 * \param[in] zetas Vector of target values for independent spline
 * \param[in] indep Name of independent spline
 * \param[in] columns Names of splines to evaluate
 * \param[out] gc GenericContainer with vectors of third derivatives (map: spline_name -> vector of d³y/dx³)
 */
void eval2_DDD( vec_real_type const & zetas, string_view indep, vec_string_type const & columns, GenericContainer & gc )
  const
{ this->eval2_DDD( zetas, this->get_position( indep ), columns, gc ); }
