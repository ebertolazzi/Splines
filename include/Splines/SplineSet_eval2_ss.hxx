
/**
 * \brief Evaluate a specific spline using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] name Name of spline to evaluate
 * \return Value of spline `name` at computed x
 */
real_type eval2( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
{ return this->eval2( zeta, this->get_position( indep ), x, this->get_position( name ) ); }

/**
 * \brief Evaluate first derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] name Name of spline to evaluate
 * \return First derivative of spline `name` with respect to independent spline value
 */
real_type eval2_D( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
{ return this->eval2_D( zeta, this->get_position( indep ), x, this->get_position( name ) ); }

/**
 * \brief Evaluate second derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] name Name of spline to evaluate
 * \return Second derivative of spline `name` with respect to independent spline value
 */
real_type eval2_DD( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
{ return this->eval2_DD( zeta, this->get_position( indep ), x, this->get_position( name ) ); }

/**
 * \brief Evaluate third derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value
 * \param[in] name Name of spline to evaluate
 * \return Third derivative of spline `name` with respect to independent spline value
 */
real_type eval2_DDD( real_type const zeta, string_view const indep, real_type & x, string_view const name ) const
{ return this->eval2_DDD( zeta, this->get_position( indep ), x, this->get_position( name ) ); }

/**
 * \brief Evaluate a specific spline using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Value of spline `name` at computed x
 */
real_type eval2( real_type const zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return this->eval2( zeta, indep, x, name );
}

/**
 * \brief Evaluate first derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return First derivative of spline `name` with respect to independent spline value
 */
real_type eval2_D( real_type const zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return this->eval2_D( zeta, indep, x, name );
}

/**
 * \brief Evaluate second derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Second derivative of spline `name` with respect to independent spline value
 */
real_type eval2_DD( real_type const zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return this->eval2_DD( zeta, indep, x, name );
}

/**
 * \brief Evaluate third derivative using another spline as independent variable
 *
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Third derivative of spline `name` with respect to independent spline value
 */
real_type eval2_DDD( real_type const zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return this->eval2_DDD( zeta, indep, x, name );
}

#ifdef AUTODIFF_SUPPORT

/**
 * \brief Evaluate with first-order automatic differentiation using independent spline
 *
 * \param[in] zeta Dual number target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] name Name of spline to evaluate
 * \return Dual number containing value and derivative
 */
autodiff::dual1st eval2(
  autodiff::dual1st const & zeta,
  string_view const         indep,
  real_type &               x,
  string_view const         name ) const
{
  real_type         zv = zeta.val;
  autodiff::dual1st res;
  res.val  = eval2( zv, indep, x, name );
  res.grad = eval_D( x, name ) * zeta.grad;  // x gia calcolato
  return res;
}

/**
 * \brief Evaluate with first-order automatic differentiation using independent spline
 *
 * \param[in] zeta Dual number target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Dual number containing value and derivative
 */
autodiff::dual1st eval2( autodiff::dual1st const & zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return eval2( zeta, indep, x, name );
}

/**
 * \brief Evaluate with second-order automatic differentiation using independent spline
 *
 * \param[in] zeta Dual number target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] name Name of spline to evaluate
 * \return Dual number containing value, first and second derivatives
 */
autodiff::dual2nd eval2(
  autodiff::dual2nd const & zeta,
  string_view const         indep,
  real_type &               x,
  string_view const         name ) const
{
  // Estrazione dei valori dalla variabile dual2nd
  real_type zv = zeta.val.val;   // valore reale di zeta
  real_type zg = zeta.grad.val;  // derivata prima: dzeta/dX

  // Calcolo delle derivate della funzione
  real_type dfx = eval2_D( zv, indep, x, name );  // derivata prima rispetto a zeta
  real_type dxx = eval_DD( x, name );             // derivata seconda rispetto a x

  // Valore della funzione
  real_type f_val = eval( x, name );

  // Costruzione dell'oggetto dual2nd
  autodiff::dual2nd res;

  // Imposta il valore della funzione
  res.val.val = f_val;

  // Calcolo della derivata prima: df/dX = (df/dzeta) * (dzeta/dX)
  real_type dF_dX = dfx * zg;

  // Imposta la derivata prima
  res.val.grad = dF_dX;  // per coerenza interna
  res.grad.val = dF_dX;  // derivata prima

  // Calcolo della derivata seconda:
  // d²f/dX² = (df/dzeta) * (d²zeta/dX²) + (d²f/dzeta²) * (dzeta/dX)²
  // Nel codice originale, dxx è probabilmente la derivata seconda di f rispetto a zeta
  res.grad.grad = dfx * zeta.grad.grad + dxx * ( zg * zg );

  return res;
}

/**
 * \brief Evaluate with second-order automatic differentiation using independent spline
 *
 * \param[in] zeta Dual number target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Dual number containing value, first and second derivatives
 */
autodiff::dual2nd eval2( autodiff::dual2nd const & zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return eval2( zeta, indep, x, name );
}

/**
 * \brief Generic template for automatic differentiation using independent spline
 *
 * \tparam T Input type
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[out] x Computed x-value (real type)
 * \param[in] name Name of spline to evaluate
 * \return Result matching input type
 */
template <typename T> auto eval2( T const & zeta, string_view const indep, real_type & x, string_view const name ) const
{
  if constexpr ( std::is_arithmetic<T>::value )
  {
    // Se T è un tipo numerico (int, float, double, etc.), promuovi a real_type
    return eval2( static_cast<real_type>( zeta ), indep, x, name );
  }
  else
  {
    // Altrimenti deduce automaticamente il tipo duale appropriato
    return eval2( autodiff::detail::to_dual( zeta ), indep, x, name );
  }
}

/**
 * \brief Generic template for automatic differentiation using independent spline
 *
 * \tparam T Input type
 * \param[in] zeta Target value of independent spline
 * \param[in] indep Name of independent spline
 * \param[in] name Name of spline to evaluate
 * \return Result matching input type
 */
template <typename T> auto eval2( T const & zeta, string_view const indep, string_view const name ) const
{
  real_type x;
  return eval2<T>( zeta, indep, x, name );
}

#endif
