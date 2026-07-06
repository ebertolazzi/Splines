/**
 * \brief Evaluate a specific spline using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] spl Index of spline to evaluate
 * \return Value of spline `spl` at computed x
 */
real_type eval2( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
{
  intersect( indep, zeta, x );
  return m_splines[spl]->eval( x );
}

/**
 * \brief Evaluate first derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] spl Index of spline to evaluate
 * \return First derivative with respect to independent spline value
 */
real_type eval2_D( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
{
  Spline const * S = intersect( indep, zeta, x );
  return m_splines[spl]->D( x ) / S->D( x );
}

/**
 * \brief Evaluate second derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] spl Index of spline to evaluate
 * \return Second derivative with respect to independent spline value
 */
real_type eval2_DD( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
{
  Spline const *  S   = intersect( indep, zeta, x );
  real_type const dt  = 1 / S->D( x );
  real_type const dt2 = dt * dt;
  real_type const ddt = -S->DD( x ) * ( dt * dt2 );
  auto &          SPL = m_splines[spl];
  return SPL->DD( x ) * dt2 + SPL->D( x ) * ddt;
}

/**
 * \brief Evaluate third derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value
 * \param[in] spl Index of spline to evaluate
 * \return Third derivative with respect to independent spline value
 */
real_type eval2_DDD( real_type const zeta, integer const indep, real_type & x, integer const spl ) const
{
  Spline const *  S    = intersect( indep, zeta, x );
  real_type const dt   = 1 / S->D( x );
  real_type const dt3  = dt * dt * dt;
  real_type const ddt  = -S->DD( x ) * dt3;
  real_type const dddt = 3 * ( ddt * ddt ) / dt - S->DDD( x ) * ( dt * dt3 );

  auto & SPL = m_splines[spl];
  return SPL->DDD( x ) * dt3 + 3 * SPL->DD( x ) * dt * ddt + SPL->D( x ) * dddt;
}

/**
 * \brief Evaluate a specific spline using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Value of spline `spl` at computed x
 */
real_type eval2( real_type const zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2( zeta, indep, x, spl );
}

/**
 * \brief Evaluate first derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return First derivative with respect to independent spline value
 */
real_type eval2_D( real_type const zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2_D( zeta, indep, x, spl );
}

/**
 * \brief Evaluate second derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Second derivative with respect to independent spline value
 */
real_type eval2_DD( real_type const zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2_DD( zeta, indep, x, spl );
}

/**
 * \brief Evaluate third derivative using another spline as independent variable (index version)
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Third derivative with respect to independent spline value
 */
real_type eval2_DDD( real_type const zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2_DDD( zeta, indep, x, spl );
}

#ifdef AUTODIFF_SUPPORT

/**
 * \brief Evaluate with first-order automatic differentiation (index version)
 *
 * \param[in] zeta Dual number target value
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] spl Index of spline to evaluate
 * \return Dual number containing value and derivative
 */
autodiff::dual1st eval2( autodiff::dual1st const & zeta, integer const indep, real_type & x, integer const spl ) const
{
  real_type         zv = zeta.val;
  autodiff::dual1st res;
  res.val  = eval2( zv, indep, x, spl );
  res.grad = eval_D( x, spl ) * zeta.grad;  // x gia calcolato
  return res;
}

/**
 * \brief Evaluate with first-order automatic differentiation (index version)
 *
 * \param[in] zeta Dual number target value
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Dual number containing value and derivative
 */
autodiff::dual1st eval2( autodiff::dual1st const & zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2( zeta, indep, x, spl );
}

/**
 * \brief Evaluate with second-order automatic differentiation (index version)
 *
 * \param[in] zeta Dual number target value
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] spl Index of spline to evaluate
 * \return Dual number containing value, first and second derivatives
 */
autodiff::dual2nd eval2( autodiff::dual2nd const & zeta, integer const indep, real_type & x, integer const spl ) const
{
  // Estrai il valore reale di zeta e la sua derivata prima
  real_type zv = zeta.val.val;   // valore di zeta
  real_type zg = zeta.grad.val;  // dzeta/dX (derivata prima di zeta)

  // Calcola la derivata di f rispetto a zeta (o x?) e aggiorna eventualmente x
  real_type dfx = eval2_D( zv, indep, x, spl );

  // Calcola la derivata seconda di f rispetto a x
  real_type dxx = eval_DD( x, spl );

  // Valore di f nel punto x
  real_type f_val = eval( x, spl );

  // Costruisci l'oggetto dual2nd correttamente
  autodiff::dual2nd res;

  // Imposta il valore della funzione
  res.val.val = f_val;

  // Calcola la derivata prima: df/dX = (df/dzeta) * (dzeta/dX)
  // oppure: df/dX = (df/dx) * (dx/dX) a seconda di cosa rappresenta dfx
  real_type dF_dX = dfx * zg;

  // Imposta la derivata prima
  res.val.grad = dF_dX;  // per coerenza interna
  res.grad.val = dF_dX;  // derivata prima

  // Calcola la derivata seconda
  // d²f/dX² = (df/dzeta) * (d²zeta/dX²) + (d²f/dzeta²) * (dzeta/dX)²
  // oppure se dfx è df/dx: d²f/dX² = (df/dx) * (d²x/dX²) + (d²f/dx²) * (dx/dX)²
  // Qui assumiamo che dxx sia la derivata seconda appropriata
  res.grad.grad = dfx * zeta.grad.grad + dxx * ( zg * zg );

  return res;
}

/**
 * \brief Evaluate with second-order automatic differentiation (index version)
 *
 * \param[in] zeta Dual number target value
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Dual number containing value, first and second derivatives
 */
autodiff::dual2nd eval2( autodiff::dual2nd const & zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2( zeta, indep, x, spl );
}

/**
 * \brief Generic template for automatic differentiation (index version)
 *
 * \tparam T Input type
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] spl Index of spline to evaluate
 * \return Result matching input type
 */
template <typename T> auto eval2( T const & zeta, integer const indep, real_type & x, integer const spl ) const
{
  if constexpr ( std::is_arithmetic<T>::value )
  {
    // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
    return eval2( static_cast<real_type>( zeta ), indep, x, spl );
  }
  else
  {
    // Altrimenti deduce automaticamente il tipo duale appropriato
    return eval2( autodiff::detail::to_dual( zeta ), indep, x, spl );
  }
}

/**
 * \brief Generic template for automatic differentiation (index version)
 *
 * \tparam T Input type
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Index of independent spline
 * \param[in] spl Index of spline to evaluate
 * \return Result matching input type
 */
template <typename T> auto eval2( T const & zeta, integer const indep, integer const spl ) const
{
  real_type x;
  return eval2( zeta, indep, x, spl );
}

#endif
