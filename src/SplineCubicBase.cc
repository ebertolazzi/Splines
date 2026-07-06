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
 |    ____      _     _        ____        _ _              ____
 |   / ___|   _| |__ (_) ___  / ___| _ __ | (_)_ __   ___  | __ )  __ _ ___  ___
 |  | |  | | | | '_ \| |/ __| \___ \| '_ \| | | '_ \ / _ \ |  _ \ / _` / __|/ _ \
 |  | |__| |_| | |_) | | (__   ___) | |_) | | | | | |  __/ | |_) | (_| \__ \  __/
 |   \____\__,_|_.__/|_|\___| |____/| .__/|_|_|_| |_|\___| |____/ \__,_|___/\___|
 |                                  |_|
\*/

#define AUTODIFF_SUPPORT
#include "Splines/Splines.hh"
#include "GenericContainer/GenericContainerInterface_nlohmann.hh"

namespace Splines
{
  namespace
  {
    void
    fill_cubic_export_data( CubicSplineBase const & spline, GenericContainer & gc )
    {
      integer const nseg = spline.num_points() > 0 ? spline.num_points() - 1 : 0;

      auto & root = gc.set_map();
      root["name"].set_string( spline.name() );
      root["spline_type"].set_string( spline.type_name() );
      root["polynomial_order"].set_int( 3 );
      root["num_points"].set_int( spline.num_points() );
      root["num_segments"].set_int( nseg );
      root["closed"].set_bool( spline.is_closed() );
      root["cyclic"].set_bool( spline.is_closed() );
      root["bounded"].set_bool( spline.is_bounded() );
      root["can_extend"].set_bool( !spline.is_bounded() );
      root["extended_constant"].set_bool( spline.is_extended_constant() );
      if ( spline.num_points() > 0 )
      {
        root["xmin"].set_real( spline.x_min() );
        root["xmax"].set_real( spline.x_max() );
      }
      else
      {
        root["xmin"].set_real( 0 );
        root["xmax"].set_real( 0 );
      }

      auto & segments = root["segments"].set_vector( static_cast<unsigned>( nseg ) );
      if ( nseg <= 0 ) return;

      std::vector<real_type> coeffs( static_cast<size_t>( 4 * nseg ) );
      std::vector<real_type> nodes( static_cast<size_t>( spline.num_points() ) );
      spline.coeffs( coeffs.data(), nodes.data(), true );

      for ( integer i = 0; i < nseg; ++i )
      {
        size_t ii  = i;
        size_t ii4 = 4*ii;
        
        auto & seg = segments[static_cast<unsigned>( i )].set_map();
        seg["index"].set_int( i );
        seg["x"].set_real( nodes[ii] );
        seg["x_next"].set_real( nodes[ii+1] );
        seg["h"].set_real( nodes[ii+1] - nodes[ii] );
        seg["a"].set_real( coeffs[ii4+0] );
        seg["b"].set_real( coeffs[ii4+1] );
        seg["c"].set_real( coeffs[ii4+2] );
        seg["d"].set_real( coeffs[ii4+3] );
      }
    }
  }  // namespace


  void CubicSplineBase::copy_spline( CubicSplineBase const & S )
  {
    // 1. Protezione contro l'auto-assegnazione (Cruciale!)
    if ( this == &S ) return;

    // 2. Alloca memoria
    CubicSplineBase::reserve( S.m_npts );
    m_npts = S.m_npts;

    // 3. Copia in blocco (Bulk Copy)
    if ( m_npts > 0 )
    {
      // Calcolo size in bytes una sola volta
      size_t const byte_size = m_npts * sizeof( real_type );

      std::memcpy( m_X, S.m_X, byte_size );
      std::memcpy( m_Y, S.m_Y, byte_size );
      std::memcpy( m_Yp, S.m_Yp, byte_size );
    }

    copy_flags( S );
  }

  void CubicSplineBase::set_range( real_type xmin, real_type xmax )
  {
    // 1. Calcola il fattore di scala PRIMA che m_X venga modificato.
    //    Se stiamo scalando le derivate (y' = dy/dx), il fattore è OldRange / NewRange.
    real_type const dx_old = m_X[m_npts - 1] - m_X[0];
    real_type const dx_new = xmax - xmin;

    // Protezione divisione per zero e calcolo rapporto
    // Se dx_new è 0, evitiamo NaN.
    real_type const recS = ( std::abs( dx_new ) > 1e-12 ) ? dx_new / dx_old : 1.0;

    // 2. Chiama la classe base per scalare m_X
    Spline::set_range( xmin, xmax );

    // 3. Check ottimizzazione: Se il fattore è 1 (nessun cambio scala), esci.
    if ( std::abs( recS - 1.0 ) < 1e-12 ) return;

    // 4. Applica il fattore di scala usando Eigen (SIMD)
    Eigen::Map<Vec> map_vals{ m_X, m_npts };
    map_vals -= m_X[0];
    map_vals *= recS;
    map_vals += xmin;
    m_X[m_npts - 1] = xmax;
  }


  void CubicSplineBase::y_min_max(
    integer &   i_min_pos,
    real_type & x_min_pos,
    real_type & y_min,
    integer &   i_max_pos,
    real_type & x_max_pos,
    real_type & y_max ) const
  {
    UTILS_ASSERT( m_npts > 0, "CubicSplineBase[{}]::y_min_max() empty spline!", m_name );
    // find max min alongh the nodes
    i_min_pos = i_max_pos = 0;
    x_min_pos = x_max_pos = m_X[0];
    y_min = y_max = m_Y[0];
    PolynomialRoots::Quadratic q;
    for ( integer i = 1; i < m_npts; ++i )
    {
      real_type const & X0  = m_X[i - 1];
      real_type const & X1  = m_X[i];
      real_type const & P0  = m_Y[i - 1];
      real_type const & P1  = m_Y[i];
      real_type const & DP0 = m_Yp[i - 1];
      real_type const & DP1 = m_Yp[i];
      real_type const   H   = X1 - X0;
      real_type         A, B, C, D;
      Hermite3_to_poly( H, P0, P1, DP0, DP1, A, B, C, D );
      q.setup( 3 * A, 2 * B, C );
      real_type     r[2];
      integer const nr = q.getRootsInOpenRange( 0, H, r );
      for ( integer j = 0; j < nr; ++j )
      {
        real_type const rr = r[j];
        real_type const yy = ( ( ( A * rr ) + B ) * rr + C ) * rr + D;
        if ( yy > y_max )
        {
          y_max     = yy;
          x_max_pos = X0 + rr;
          i_max_pos = i;
        }
        else if ( yy < y_min )
        {
          y_min     = yy;
          x_min_pos = X0 + rr;
          i_min_pos = i;
        }
      }
      if ( P1 > y_max )
      {
        y_max     = P1;
        x_max_pos = X1;
        i_max_pos = i;
      }
      else if ( P1 < y_min )
      {
        y_min     = P1;
        x_min_pos = X1;
        i_min_pos = i;
      }
    }
  }

  void CubicSplineBase::y_min_max(
    vector<integer> &   i_min_pos,
    vector<real_type> & x_min_pos,
    vector<real_type> & y_min,
    vector<integer> &   i_max_pos,
    vector<real_type> & x_max_pos,
    vector<real_type> & y_max ) const
  {
    constexpr real_type epsi = 1e-8;
    i_min_pos.clear();
    i_max_pos.clear();
    x_min_pos.clear();
    x_max_pos.clear();
    y_min.clear();
    y_max.clear();
    UTILS_ASSERT( m_npts > 0, "CubicSplineBase[{}]::y_min_max() empty spline!", m_name );
    // find max min along the nodes
    if ( m_Yp[0] >= 0 )
    {
      y_min.emplace_back( m_Y[0] );
      x_min_pos.emplace_back( m_X[0] );
      i_min_pos.emplace_back( 0 );
    }
    if ( m_Yp[0] <= 0 )
    {
      y_max.emplace_back( m_Y[0] );
      x_max_pos.emplace_back( m_X[0] );
      i_max_pos.emplace_back( 0 );
    }
    PolynomialRoots::Quadratic q;
    for ( integer i = 1; i < m_npts; ++i )
    {
      real_type const & X0  = m_X[i - 1];
      real_type const & X1  = m_X[i];
      real_type const & P0  = m_Y[i - 1];
      real_type const & P1  = m_Y[i];
      real_type const & DP0 = m_Yp[i - 1];
      real_type const & DP1 = m_Yp[i];
      real_type const   H   = X1 - X0;
      real_type         A, B, C, D;
      Hermite3_to_poly( H, P0, P1, DP0, DP1, A, B, C, D );
      q.setup( 3 * A, 2 * B, C );
      real_type     r[2];
      integer const nr = q.getRootsInOpenRange( 0, H, r );
      for ( integer j = 0; j < nr; ++j )
      {
        real_type const rr  = r[j];
        real_type const yy  = ( ( ( A * rr ) + B ) * rr + C ) * rr + D;
        real_type const ddy = 3 * A * rr + B;
        if ( ddy > 0 )
        {
          y_min.emplace_back( yy );
          x_min_pos.emplace_back( X0 + rr );
          i_min_pos.emplace_back( i );
        }
        else if ( ddy < 0 )
        {
          y_max.emplace_back( yy );
          x_max_pos.emplace_back( X0 + rr );
          i_max_pos.emplace_back( i );
        }
      }
      if ( i + 1 >= m_npts ) continue;
      if ( std::abs( DP1 ) > ( m_X[i + 1] - m_X[i - 1] ) * epsi ) continue;
      real_type const & X2  = m_X[i + 1];
      real_type const & P2  = m_Y[i + 1];
      real_type const & DP2 = m_Yp[i + 1];
      real_type         A1, B1, C1, D1;
      Hermite3_to_poly( X2 - X1, P1, P2, DP1, DP2, A1, B1, C1, D1 );
      real_type const DD = 2 * A * H + B;
      if ( DD >= 0 && B1 >= 0 )
      {
        y_min.emplace_back( P1 );
        x_min_pos.emplace_back( X1 );
        i_min_pos.emplace_back( i );
      }
      else if ( DD <= 0 && B1 <= 0 )
      {
        y_max.emplace_back( P1 );
        x_max_pos.emplace_back( X1 );
        i_max_pos.emplace_back( i );
      }
    }
    if ( m_Yp[m_npts - 1] <= 0 )
    {
      y_min.emplace_back( m_Y[m_npts - 1] );
      x_min_pos.emplace_back( m_X[m_npts - 1] );
      i_min_pos.emplace_back( 0 );
    }
    if ( m_Yp[m_npts - 1] >= 0 )
    {
      y_max.emplace_back( m_Y[m_npts - 1] );
      x_max_pos.emplace_back( m_X[m_npts - 1] );
      i_max_pos.emplace_back( m_npts - 1 );
    }
  }

  void CubicSplineBase::D( real_type const x, real_type dd[2] ) const
  {
    std::pair<integer, real_type> res( 0, x );
    m_search.find( res );
    integer const   ni = res.first;
    real_type const X  = res.second;
    real_type       dx = X - m_X[ni];
    real_type       DX = m_X[ni + 1] - m_X[ni];
    real_type       base[4], base_D[4];
    Splines::Hermite3( dx, DX, base );
    Splines::Hermite3_D( dx, DX, base_D );
    dd[0] = base[0] * m_Y[ni] + base[1] * m_Y[ni + 1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni + 1];
    dd[1] = base_D[0] * m_Y[ni] + base_D[1] * m_Y[ni + 1] + base_D[2] * m_Yp[ni] + base_D[3] * m_Yp[ni + 1];
  }

  void CubicSplineBase::DD( real_type const x, real_type dd[3] ) const
  {
    std::pair<integer, real_type> res( 0, x );
    m_search.find( res );
    integer const   ni = res.first;
    real_type const X  = res.second;
    real_type       base[4], base_D[4], base_DD[4];
    real_type       dx = X - m_X[ni];
    real_type       DX = m_X[ni + 1] - m_X[ni];
    Splines::Hermite3( dx, DX, base );
    Splines::Hermite3_D( dx, DX, base_D );
    Splines::Hermite3_DD( dx, DX, base_DD );
    dd[0] = base[0] * m_Y[ni] + base[1] * m_Y[ni + 1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni + 1];
    dd[1] = base_D[0] * m_Y[ni] + base_D[1] * m_Y[ni + 1] + base_D[2] * m_Yp[ni] + base_D[3] * m_Yp[ni + 1];
    dd[2] = base_DD[0] * m_Y[ni] + base_DD[1] * m_Y[ni + 1] + base_DD[2] * m_Yp[ni] + base_DD[3] * m_Yp[ni + 1];
  }

  real_type CubicSplineBase::id_eval( integer const ni, real_type const x ) const
  {
    if ( m_curve_can_extend && m_curve_extended_constant )
    {
      // estendo solo quando esco effettivamente
      if ( x < m_X[0] ) return m_Y[0];
      if ( x > m_X[m_npts - 1] ) return m_Y[m_npts - 1];
    }
    real_type base[4];
    Splines::Hermite3( x - m_X[ni], m_X[ni + 1] - m_X[ni], base );
    return base[0] * m_Y[ni] + base[1] * m_Y[ni + 1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni + 1];
  }

  real_type CubicSplineBase::id_D( integer const ni, real_type const x ) const
  {
    if ( m_curve_can_extend && m_curve_extended_constant )
    {
      // estendo solo quando esco effettivamente
      if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
    }
    real_type base_D[4];
    Splines::Hermite3_D( x - m_X[ni], m_X[ni + 1] - m_X[ni], base_D );
    return base_D[0] * m_Y[ni] + base_D[1] * m_Y[ni + 1] + base_D[2] * m_Yp[ni] + base_D[3] * m_Yp[ni + 1];
  }

  real_type CubicSplineBase::id_DD( integer const ni, real_type const x ) const
  {
    if ( m_curve_can_extend && m_curve_extended_constant )
    {
      // estendo solo quando esco effettivamente
      if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
    }
    real_type base_DD[4];
    Splines::Hermite3_DD( x - m_X[ni], m_X[ni + 1] - m_X[ni], base_DD );
    return base_DD[0] * m_Y[ni] + base_DD[1] * m_Y[ni + 1] + base_DD[2] * m_Yp[ni] + base_DD[3] * m_Yp[ni + 1];
  }

  real_type CubicSplineBase::id_DDD( integer const ni, real_type const x ) const
  {
    if ( m_curve_can_extend && m_curve_extended_constant )
    {
      // estendo solo quando esco effettivamente
      if ( x < m_X[0] || x > m_X[m_npts - 1] ) return 0;
    }
    real_type base_DDD[4];
    Splines::Hermite3_DDD( x - m_X[ni], m_X[ni + 1] - m_X[ni], base_DDD );
    return base_DDD[0] * m_Y[ni] + base_DDD[1] * m_Y[ni + 1] + base_DDD[2] * m_Yp[ni] + base_DDD[3] * m_Yp[ni + 1];
  }

  void CubicSplineBase::write_to_stream( ostream_type & s ) const
  {
    integer const nseg = m_npts > 0 ? m_npts - 1 : 0;
    for ( integer i = 0; i < nseg; ++i )
      fmt::print(
        s,
        "segment N.{:4} X:[{:.5},{:.5}] Y:[{:.5},{:.5}] Yp:[{:.5},{:.5}] slope: {}\n",
        i,
        m_X[i],
        m_X[i + 1],
        m_Y[i],
        m_Y[i + 1],
        m_Yp[i],
        m_Yp[i + 1],
        ( m_Y[i + 1] - m_Y[i] ) / ( m_X[i + 1] - m_X[i] ) );
  }

  void CubicSplineBase::export_csv( ostream_type & s ) const
  {
    fmt::print( s, "x,a,b,c,d\n" );

    if ( m_npts <= 0 ) return;

    integer const nseg = m_npts > 0 ? m_npts - 1 : 0;
    if ( nseg > 0 )
    {
      std::vector<real_type> coeffs( static_cast<size_t>( 4 * nseg ) );
      std::vector<real_type> nodes( static_cast<size_t>( m_npts ) );
      this->coeffs( coeffs.data(), nodes.data(), true );

      for ( integer i = 0; i < nseg; ++i ) {
        size_t ii  = i;
        size_t ii4 = ii*4;
        fmt::print(
          s,
          "{:.17g},{:.17g},{:.17g},{:.17g},{:.17g}\n",
          nodes[ii], coeffs[ii4+0], coeffs[ii4+1], coeffs[ii4+2], coeffs[ii4+3] );
      }
    }

    fmt::print(
      s,
      "{:.17g},{:.17g},{:.17g},{:.17g},{:.17g}\n",
      m_X[m_npts - 1],
      m_Y[m_npts - 1],
      0.0,
      0.0,
      0.0 );
  }

  void CubicSplineBase::export_json( ostream_type & s ) const
  {
    GenericContainer gc;
    fill_cubic_export_data( *this, gc );
    nlohmann::json const json = gc;
    s << json.dump( 2 ) << '\n';
  }

  void CubicSplineBase::export_yaml( ostream_type & s ) const
  {
    GenericContainer gc;
    fill_cubic_export_data( *this, gc );
    gc.to_yaml( s );
  }

  void CubicSplineBase::build(
    real_type const x[],
    integer const   incx,
    real_type const y[],
    integer const   incy,
    real_type const yp[],
    integer const   incyp,
    integer         n )
  {
    this->reserve( n );
    for ( integer i = 0; i < n; ++i )
    {
      m_X[i]  = x[i * incx];
      m_Y[i]  = y[i * incy];
      m_Yp[i] = yp[i * incyp];
    }
    m_npts = n;
    m_search.must_reset();
  }

  void CubicSplineBase::reserve( integer npts )
  {
    if ( m_external_alloc && npts <= m_npts_reserved )
    {
      // nothing to do!, already allocated
    }
    else
    {
      m_npts_reserved = npts;
      m_mem_cubic.reallocate( 3 * npts );
      m_X              = m_mem_cubic( npts );
      m_Y              = m_mem_cubic( npts );
      m_Yp             = m_mem_cubic( npts );
      m_external_alloc = false;
    }
    m_npts = 0;
  }

  void CubicSplineBase::clear()
  {
    if ( !m_external_alloc ) m_mem_cubic.free();
    m_npts = m_npts_reserved = 0;
    m_external_alloc         = false;
    m_X = m_Y = m_Yp = nullptr;
  }

  integer  // order
  CubicSplineBase::coeffs( real_type cfs[], real_type nodes[], bool transpose ) const
  {
    UTILS_ASSERT( m_npts >= 2, "CubicSplineBase::coeffs, npts={} must be >= 2\n", m_npts );

    integer const n = m_npts - 1;
    for ( integer i = 0; i < n; ++i )
    {
      real_type const H = m_X[i + 1] - m_X[i];
      real_type       a, b, c, d;

      Hermite3_to_poly( H, m_Y[i], m_Y[i + 1], m_Yp[i], m_Yp[i + 1], a, b, c, d );

      if ( transpose )
      {
        cfs[4 * i + 3] = a;
        cfs[4 * i + 2] = b;
        cfs[4 * i + 1] = c;
        cfs[4 * i + 0] = d;
      }
      else
      {
        cfs[i + 3 * n] = a;
        cfs[i + 2 * n] = b;
        cfs[i + 1 * n] = c;
        cfs[i + 0 * n] = d;
      }
    }
    if ( m_npts > 0 ) std::memcpy( nodes, m_X, m_npts * sizeof( real_type ) );
    return 4;
  }

  autodiff::dual1st CubicSplineBase::eval( autodiff::dual1st const & x ) const
  {
    if ( m_curve_can_extend && m_curve_extended_constant )
    {
      // Estendo solo quando esco effettivamente
      if ( x.val < m_X[0] )
      {
        autodiff::dual1st res;
        res.val  = m_Y[0];
        res.grad = 0;
        return res;
      }
      if ( x.val > m_X[m_npts - 1] )
      {
        autodiff::dual1st res;
        res.val  = m_Y[m_npts - 1];
        res.grad = 0;
        return res;
      }
    }

    std::pair<integer, real_type> X( 0, x.val );
    m_search.find( X );
    integer ni  = X.first;
    integer ni1 = ni + 1;

    autodiff::dual1st XX;
    XX.val  = X.second - m_X[ni];
    XX.grad = x.grad;

    autodiff::dual1st base[4];
    Splines::Hermite3<autodiff::dual1st>( XX, m_X[ni1] - m_X[ni], base );

    return base[0] * m_Y[ni] + base[1] * m_Y[ni1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  autodiff::dual2nd CubicSplineBase::eval( autodiff::dual2nd const & x ) const
  {
    if ( m_curve_can_extend && m_curve_extended_constant )
    {
      // Estendo solo quando esco effettivamente
      if ( x.val.val < m_X[0] )
      {
        autodiff::dual2nd res;
        res.val.val   = m_Y[0];
        res.val.grad  = 0;
        res.grad.val  = 0;
        res.grad.grad = 0;
        return res;
      }
      if ( x.val.val > m_X[m_npts - 1] )
      {
        autodiff::dual2nd res;
        res.val.val   = m_Y[m_npts - 1];
        res.val.grad  = 0;
        res.grad.val  = 0;
        res.grad.grad = 0;
        return res;
      }
    }

    std::pair<integer, real_type> X( 0, x.val.val );
    m_search.find( X );
    integer ni  = X.first;
    integer ni1 = ni + 1;

    autodiff::dual2nd XX;
    // Costruzione corretta di XX = (x - m_X[ni])
    // Dobbiamo impostare tutti i campi
    XX.val.val  = X.second - m_X[ni];  // valore
    XX.val.grad = x.val.grad;          // derivata prima del valore

    // La derivata prima di XX è la stessa di x
    XX.grad.val = x.grad.val;  // derivata prima
    // La derivata seconda di XX è la stessa di x
    XX.grad.grad = x.grad.grad;  // derivata seconda

    autodiff::dual2nd base[4];
    Splines::Hermite3<autodiff::dual2nd>( XX, m_X[ni1] - m_X[ni], base );

    return base[0] * m_Y[ni] + base[1] * m_Y[ni1] + base[2] * m_Yp[ni] + base[3] * m_Yp[ni1];
  }

}  // namespace Splines
