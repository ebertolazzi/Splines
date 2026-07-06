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

#include "Splines/Splines.hh"

namespace Splines
{

  void SplineVec::allocate( integer const dim, integer const npts )
  {
    UTILS_ASSERT( dim > 0, "SplineVec[{}]::build expected positive dim = {}\n", m_name, dim );
    UTILS_ASSERT( npts > 1, "SplineVec[{}]::build expected npts = {} greather than 1\n", m_name, npts );
    m_dim  = dim;
    m_npts = npts;

    // Allocate contiguous memory: (2*dim + 1)*npts for X, Y, Yp arrays
    m_mem.reallocate( ( 2 * dim + 1 ) * npts );
    m_mem_p.reallocate( 2 * dim );

    // Setup pointer arrays
    m_Y  = m_mem_p( m_dim );  // Pointers to Y components
    m_Yp = m_mem_p( m_dim );  // Pointers to Yp components
    m_X  = m_mem( m_npts );   // Knot vector

    // Distribute memory blocks to each component
    for ( integer spl = 0; spl < m_dim; ++spl )
    {
      m_Y[spl]  = m_mem( m_npts );  // Points for dimension spl
      m_Yp[spl] = m_mem( m_npts );  // Derivatives for dimension spl
    }

    // Verify all allocated memory has been assigned
    m_mem.must_be_empty( "SplineVec::build, baseValue" );
    m_mem_p.must_be_empty( "SplineVec::build, basePointer" );
  }

  void SplineVec::compute_chords()
  {
    integer const nn = m_npts - 1;
    switch ( m_dim )
    {
      case 2:
        // Optimized 2D case using hypot(x, y)
        for ( integer j = 0; j < nn; ++j )
        {
          real_type const dx = m_Y[0][j + 1] - m_Y[0][j];
          real_type const dy = m_Y[1][j + 1] - m_Y[1][j];
          m_X[j]             = std::hypot( dx, dy );
        }
        break;
      case 3:
        // Optimized 3D case using nested hypot for numerical stability
        for ( integer j = 0; j < nn; ++j )
        {
          real_type const dx = m_Y[0][j + 1] - m_Y[0][j];
          real_type const dy = m_Y[1][j + 1] - m_Y[1][j];
          real_type const dz = m_Y[2][j + 1] - m_Y[2][j];
          m_X[j]             = std::hypot( dx, dy, dz );
        }
        break;
      default:
        // General n-dimensional case
        for ( integer j = 0; j < nn; ++j )
        {
          real_type l = 0;
          for ( integer k = 0; k < m_dim; ++k )
          {
            real_type const d = m_Y[k][j + 1] - m_Y[k][j];
            l += d * d;
          }
          m_X[j] = std::sqrt( l );
        }
        break;
    }
    m_search.must_reset();  // Invalidate cached search data
  }

  real_type SplineVec::eval( real_type const x, integer const i ) const
  {
    std::pair<integer, real_type> res( 0, x );
    m_search.find( res );
    real_type     base[4];
    integer const idx = res.first;
    Splines::Hermite3( res.second - m_X[idx], m_X[idx + 1] - m_X[idx], base );
    real_type const * Y  = m_Y[i];
    real_type const * Yp = m_Yp[i];
    return base[0] * Y[idx] + base[1] * Y[idx + 1] + base[2] * Yp[idx] + base[3] * Yp[idx + 1];
  }

  real_type SplineVec::D( real_type const x, integer const i ) const
  {
    std::pair<integer, real_type> res( 0, x );
    m_search.find( res );
    real_type     base_D[4];
    integer const idx = res.first;
    Splines::Hermite3_D( res.second - m_X[idx], m_X[idx + 1] - m_X[idx], base_D );
    real_type const * Y  = m_Y[i];
    real_type const * Yp = m_Yp[i];
    return base_D[0] * Y[idx] + base_D[1] * Y[idx + 1] + base_D[2] * Yp[idx] + base_D[3] * Yp[idx + 1];
  }

  real_type SplineVec::DD( real_type const x, integer const i ) const
  {
    std::pair<integer, real_type> res( 0, x );
    m_search.find( res );
    real_type     base_DD[4];
    integer const idx = res.first;
    Splines::Hermite3_DD( res.second - m_X[idx], m_X[idx + 1] - m_X[idx], base_DD );
    real_type const * Y  = m_Y[i];
    real_type const * Yp = m_Yp[i];
    return base_DD[0] * Y[idx] + base_DD[1] * Y[idx + 1] + base_DD[2] * Yp[idx] + base_DD[3] * Yp[idx + 1];
  }

  real_type SplineVec::DDD( real_type const x, integer const i ) const
  {
    std::pair<integer, real_type> res( 0, x );
    m_search.find( res );
    real_type     base_DDD[4];
    integer const idx = res.first;
    Splines::Hermite3_DDD( res.second - m_X[idx], m_X[idx + 1] - m_X[idx], base_DDD );
    real_type const * Y  = m_Y[i];
    real_type const * Yp = m_Yp[i];
    return base_DDD[0] * Y[idx] + base_DDD[1] * Y[idx + 1] + base_DDD[2] * Yp[idx] + base_DDD[3] * Yp[idx + 1];
  }

  void SplineVec::eval( real_type const x, real_type vals[], integer const inc ) const
  {
    std::pair<integer, real_type> res( 0, x );
    m_search.find( res );
    real_type     base[4];
    integer const idx = res.first;
    Splines::Hermite3( res.second - m_X[idx], m_X[idx + 1] - m_X[idx], base );
    real_type * v = vals;
    for ( integer j = 0; j < m_dim; ++j, v += inc )
    {
      real_type const * Y  = m_Y[j];
      real_type const * Yp = m_Yp[j];
      *v                   = base[0] * Y[idx] + base[1] * Y[idx + 1] + base[2] * Yp[idx] + base[3] * Yp[idx + 1];
    }
  }

  void SplineVec::eval_D( real_type const x, real_type vals[], integer const inc ) const
  {
    std::pair<integer, real_type> res( 0, x );
    m_search.find( res );
    real_type     base_D[4];
    integer const idx = res.first;
    Splines::Hermite3_D( res.second - m_X[idx], m_X[idx + 1] - m_X[idx], base_D );
    real_type * v = vals;
    for ( integer j = 0; j < m_dim; ++j, v += inc )
    {
      real_type const * Y  = m_Y[j];
      real_type const * Yp = m_Yp[j];
      *v = base_D[0] * Y[idx] + base_D[1] * Y[idx + 1] + base_D[2] * Yp[idx] + base_D[3] * Yp[idx + 1];
    }
  }

  void SplineVec::eval_DD( real_type const x, real_type vals[], integer const inc ) const
  {
    std::pair<integer, real_type> res( 0, x );
    m_search.find( res );
    real_type     base_DD[4];
    integer const idx = res.first;
    Splines::Hermite3_DD( res.second - m_X[idx], m_X[idx + 1] - m_X[idx], base_DD );
    real_type * v = vals;
    for ( integer j = 0; j < m_dim; ++j, v += inc )
    {
      real_type const * Y  = m_Y[j];
      real_type const * Yp = m_Yp[j];
      *v = base_DD[0] * Y[idx] + base_DD[1] * Y[idx + 1] + base_DD[2] * Yp[idx] + base_DD[3] * Yp[idx + 1];
    }
  }

  void SplineVec::eval_DDD( real_type const x, real_type vals[], integer const inc ) const
  {
    std::pair<integer, real_type> res( 0, x );
    m_search.find( res );
    real_type     base_DDD[4];
    integer const idx = res.first;
    Splines::Hermite3_DDD( res.second - m_X[idx], m_X[idx + 1] - m_X[idx], base_DDD );
    real_type * v = vals;
    for ( integer j = 0; j < m_dim; ++j, v += inc )
    {
      real_type const * Y  = m_Y[j];
      real_type const * Yp = m_Yp[j];
      *v = base_DDD[0] * Y[idx] + base_DDD[1] * Y[idx + 1] + base_DDD[2] * Yp[idx] + base_DDD[3] * Yp[idx + 1];
    }
  }

  void SplineVec::eval( vec_real_type const & x, GenericContainer & vals ) const
  {
    mat_real_type &   m  = vals.set_mat_real( static_cast<unsigned>( m_dim ), static_cast<unsigned>( x.size() ) );
    real_type *       v  = m.data();
    real_type const * px = x.data();
    integer           j  = static_cast<integer>( x.size() );
    while ( j-- > 0 )
    {
      eval( *px++, v, 1 );
      v += m_dim;
    }
  }

  void SplineVec::eval_D( vec_real_type const & x, GenericContainer & vals ) const
  {
    mat_real_type &   m  = vals.set_mat_real( static_cast<unsigned>( m_dim ), static_cast<unsigned>( x.size() ) );
    real_type *       v  = m.data();
    real_type const * px = x.data();
    integer           j  = static_cast<integer>( x.size() );
    while ( j-- > 0 )
    {
      eval_D( *px++, v, 1 );
      v += m_dim;
    }
  }

  void SplineVec::eval_DD( vec_real_type const & x, GenericContainer & vals ) const
  {
    mat_real_type &   m  = vals.set_mat_real( static_cast<unsigned>( m_dim ), static_cast<unsigned>( x.size() ) );
    real_type *       v  = m.data();
    real_type const * px = x.data();
    integer           j  = static_cast<integer>( x.size() );
    while ( j-- > 0 )
    {
      eval_DD( *px++, v, 1 );
      v += m_dim;
    }
  }

  void SplineVec::eval_DDD( vec_real_type const & x, GenericContainer & vals ) const
  {
    mat_real_type &   m  = vals.set_mat_real( static_cast<unsigned>( m_dim ), static_cast<unsigned>( x.size() ) );
    real_type *       v  = m.data();
    real_type const * px = x.data();
    integer           j  = static_cast<integer>( x.size() );
    while ( j-- > 0 )
    {
      eval_DDD( *px++, v, 1 );
      v += m_dim;
    }
  }

  void SplineVec::eval_DDDD( vec_real_type const & x, GenericContainer & vals ) const
  {
    mat_real_type &   m  = vals.set_mat_real( static_cast<unsigned>( m_dim ), static_cast<unsigned>( x.size() ) );
    real_type *       v  = m.data();
    real_type const * px = x.data();
    integer           j  = static_cast<integer>( x.size() );
    while ( j-- > 0 )
    {
      eval_DDDD( *px++, v, 1 );
      v += m_dim;
    }
  }

  void SplineVec::eval_DDDDD( vec_real_type const & x, GenericContainer & vals ) const
  {
    mat_real_type &   m  = vals.set_mat_real( static_cast<unsigned>( m_dim ), static_cast<unsigned>( x.size() ) );
    real_type *       v  = m.data();
    real_type const * px = x.data();
    integer           j  = static_cast<integer>( x.size() );
    while ( j-- > 0 )
    {
      eval_DDDDD( *px++, v, 1 );
      v += m_dim;
    }
  }

  void SplineVec::setup( integer const dim, integer const npts, real_type const * const Y[] )
  {
    // 1. Allocazione (imposta m_dim e m_npts internamente)
    allocate( dim, npts );

    // 2. Controllo Early Exit
    if ( m_npts <= 0 || m_dim <= 0 ) return;

    // 3. Pre-calcolo della dimensione in byte della riga.
    size_t const row_size_bytes = m_npts * sizeof( real_type );

    // 4. Loop di copia ottimizzato
    for ( integer spl = 0; spl < m_dim; ++spl ) std::memcpy( m_Y[spl], Y[spl], row_size_bytes );
  }

  void SplineVec::setup( integer const dim, integer const npts, real_type const Y[], integer const ldY )
  {
    allocate( dim, npts );
    for ( integer spl = 0; spl < m_dim; ++spl )
      for ( integer j = 0; j < m_npts; ++j ) m_Y[spl][j] = Y[spl + j * ldY];
  }

  void SplineVec::set_knots_chord_length()
  {
    compute_chords();
    integer const nn  = m_npts - 1;
    real_type     acc = 0;
    for ( integer j = 0; j <= nn; ++j )
    {
      real_type const l = m_X[j];
      m_X[j]            = acc;
      acc += l;
    }
    for ( integer j = 1; j < nn; ++j ) m_X[j] /= acc;
    m_X[nn] = 1;
    m_search.must_reset();
  }

  void SplineVec::set_knots_centripetal()
  {
    compute_chords();
    integer const nn  = m_npts - 1;
    real_type     acc = 0;
    for ( integer j = 0; j <= nn; ++j )
    {
      real_type const l = std::sqrt( m_X[j] );
      m_X[j]            = acc;
      acc += l;
    }
    for ( integer j = 1; j < nn; ++j ) m_X[j] /= acc;
    m_X[nn] = 1;
    m_search.must_reset();
  }

  void SplineVec::set_knots_foley()
  {
    compute_chords();  // m_X now contains chord lengths d_i

    integer const n = m_npts - 1;  // number of intervals
    if ( n < 2 )
    {
      // For very few points, fall back to uniform parametrization
      for ( integer i = 0; i < m_npts; ++i ) { m_X[i] = static_cast<real_type>( i ) / static_cast<real_type>( n ); }
      m_search.must_reset();
      return;
    }

    // Step 1: Store chord lengths in a separate vector
    std::vector<real_type> d( n );
    for ( integer i = 0; i < n; ++i ) { d[i] = m_X[i]; }

    // Step 2: Compute turning angles at interior points
    std::vector<real_type> theta( n - 1, real_type( 0 ) );

    for ( integer i = 1; i < n; ++i )
    {
      // Vectors for consecutive segments
      // P_i = m_Y[dim][i], where dim is the dimension (0 for x, 1 for y, etc.)
      // We need to compute the Euclidean distance between points

      // Compute vector v1 = P_i - P_{i-1}
      real_type v1x = m_Y[0][i] - m_Y[0][i - 1];
      real_type v1y = ( m_dim > 1 ) ? ( m_Y[1][i] - m_Y[1][i - 1] ) : real_type( 0 );
      // Add more dimensions if needed

      // Compute vector v2 = P_{i+1} - P_i
      real_type v2x = m_Y[0][i + 1] - m_Y[0][i];
      real_type v2y = ( m_dim > 1 ) ? ( m_Y[1][i + 1] - m_Y[1][i] ) : real_type( 0 );

      // Compute norms
      real_type norm1 = std::sqrt( v1x * v1x + v1y * v1y );
      real_type norm2 = std::sqrt( v2x * v2x + v2y * v2y );

      if ( norm1 > real_type( 0 ) && norm2 > real_type( 0 ) )
      {
        // Compute dot product
        real_type dot = v1x * v2x + v1y * v2y;

        // Compute angle between vectors (in radians)
        real_type cos_theta = dot / ( norm1 * norm2 );
        cos_theta           = std::max( real_type( -1 ), std::min( real_type( 1 ), cos_theta ) );
        theta[i - 1]        = std::acos( cos_theta );
      }
      // else: theta remains 0 for degenerate cases
    }

    // Step 3: Compute weighted chord lengths
    std::vector<real_type> weighted_d( n, real_type( 0 ) );

    // Handle first segment (i = 0)
    if ( n > 0 )
    {
      real_type w0 = real_type( 1 );
      if ( n > 1 )
      {
        // Only right angle contributes
        w0 += real_type( 1.5 ) * ( theta[0] * d[0] ) / ( d[0] + d[1] );
      }
      weighted_d[0] = d[0] * w0;
    }

    // Handle interior segments (i = 1 to n-2)
    for ( integer i = 1; i < n - 1; ++i )
    {
      real_type w = real_type( 1 );

      // Left angle contribution
      w += real_type( 1.5 ) * ( theta[i - 1] * d[i - 1] ) / ( d[i - 1] + d[i] );

      // Right angle contribution
      w += real_type( 1.5 ) * ( theta[i] * d[i] ) / ( d[i] + d[i + 1] );

      weighted_d[i] = d[i] * w;
    }

    // Handle last segment (i = n-1)
    if ( n > 1 )
    {
      real_type w_last = real_type( 1 );
      // Only left angle contributes
      w_last += real_type( 1.5 ) * ( theta[n - 2] * d[n - 2] ) / ( d[n - 2] + d[n - 1] );
      weighted_d[n - 1] = d[n - 1] * w_last;
    }

    // Step 4: Compute cumulative parameter values
    real_type acc = 0;
    for ( integer i = 0; i <= n; ++i )
    {
      if ( i == 0 ) { m_X[0] = acc; }
      else
      {
        m_X[i] = acc;
        acc += weighted_d[i - 1];
      }
    }

    // Step 5: Normalize to [0, 1]
    if ( acc > real_type( 0 ) )
    {
      for ( integer i = 1; i < n; ++i ) { m_X[i] /= acc; }
    }
    else
    {
      // Fallback to uniform if sum is zero
      for ( integer i = 0; i < m_npts; ++i ) { m_X[i] = static_cast<real_type>( i ) / static_cast<real_type>( n ); }
    }

    m_X[n] = real_type( 1 );
    m_search.must_reset();
  }

  void SplineVec::catmull_rom()
  {
    UTILS_ASSERT( m_npts >= 2, "catmull_rom, npts={} must be >= 2\n", m_npts );

    integer const n = m_npts - 1;
    integer const d = m_dim;

    real_type l1, l2, ll, a, b;
    for ( integer j = 1; j < n; ++j )
    {
      l1 = m_X[j] - m_X[j - 1];
      l2 = m_X[j + 1] - m_X[j];
      ll = l1 + l2;
      a  = ( l2 / l1 ) / ll;
      b  = ( l1 / l2 ) / ll;
      for ( integer k = 0; k < d; ++k )
        m_Yp[k][j] = a * ( m_Y[k][j] - m_Y[k][j - 1] ) + b * ( m_Y[k][j + 1] - m_Y[k][j] );
    }

    l1 = m_X[1] - m_X[0];
    l2 = m_X[2] - m_X[1];
    ll = l1 + l2;
    a  = ll / ( l1 * l2 );
    b  = ( l1 / l2 ) / ll;
    for ( integer k = 0; k < d; ++k ) m_Yp[k][0] = a * ( m_Y[k][1] - m_Y[k][0] ) - b * ( m_Y[k][2] - m_Y[k][0] );

    l1 = m_X[n] - m_X[n - 1];
    l2 = m_X[n - 1] - m_X[n - 2];
    ll = l1 + l2;
    a  = ll / ( l1 * l2 );
    b  = ( l1 / l2 ) / ll;
    for ( integer k = 0; k < d; ++k )
      m_Yp[k][n] = b * ( m_Y[k][n - 2] - m_Y[k][n] ) - a * ( m_Y[k][n - 1] - m_Y[k][n] );
  }

  void SplineVec::setup( GenericContainer const & gc )
  {
    string const where = fmt::format( "SplineVec[{}]::setup( gc ):", m_name );

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map( where ) ) { keywords.insert( pair.first ); }

    GenericContainer const & data = gc( "data", where );
    keywords.erase( "data" );

    mat_real_type Y;
    data.copyto_mat_real( Y, where );

    bool transposed = false;
    gc.get_if_exists( "transposed", transposed );
    keywords.erase( "transposed" );

    if ( transposed )
    {
      integer const dim  = static_cast<integer>( Y.num_rows() );
      integer const npts = static_cast<integer>( Y.num_cols() );
      allocate( dim, npts );
      for ( integer spl = 0; spl < m_dim; ++spl )
        for ( integer j = 0; j < m_npts; ++j ) m_Y[spl][j] = Y( spl, j );
    }
    else
    {
      integer const dim  = static_cast<integer>( Y.num_cols() );
      integer const npts = static_cast<integer>( Y.num_rows() );
      allocate( dim, npts );
      for ( integer spl = 0; spl < m_dim; ++spl )
        for ( integer j = 0; j < m_npts; ++j ) m_Y[spl][j] = Y( j, spl );
    }

    UTILS_WARNING(
      keywords.empty(),
      "{}: unused keys\n{}\n",
      where,
      [&keywords]() -> string
      {
        string res;
        for ( auto const & it : keywords )
        {
          res += it;
          res += ' ';
        };
        return res;
      }() );
  }

  real_type SplineVec::curvature( real_type x ) const
  {
    UTILS_ASSERT( m_dim == 2, "SplineVec::curvature(x={}) defined only for dim=2 found {}", x, m_dim );
    real_type D[2], DD[2];
    this->eval_D( x, D, 1 );
    this->eval_DD( x, DD, 1 );

    real_type const dx  = D[0];
    real_type const dy  = D[1];
    real_type const ddx = DD[0];
    real_type const ddy = DD[1];

    real_type const speed2 = dx * dx + dy * dy;
    real_type const speed  = std::sqrt( speed2 );

    return ( dx * ddy - dy * ddx ) / ( speed2 * speed );
  }

  real_type SplineVec::curvature_D( real_type x ) const
  {
    UTILS_ASSERT( m_dim == 2, "SplineVec::curvature(x={}) defined only for dim=2 found {}", x, m_dim );
    real_type D[2], DD[2], DDD[2];
    this->eval_D( x, D, 1 );
    this->eval_DD( x, DD, 1 );
    this->eval_DDD( x, DDD, 1 );

    const real_type dx   = D[0];
    const real_type dy   = D[1];
    const real_type ddx  = DD[0];
    const real_type ddy  = DD[1];
    const real_type dddx = DDD[0];
    const real_type dddy = DDD[1];

    const real_type dx2 = dx * dx;
    const real_type dy2 = dy * dy;

    const real_type speed2 = dx2 + dy2;  // |r'|^2
    const real_type speed  = std::sqrt( speed2 );

    const real_type a = dddy * dx - dy * dddx;
    const real_type b = 3 * ddy * ddx;

    const real_type num = dx2 * ( a - b ) + dy2 * ( a + b ) + 3 * dx * dy * ( ddx * ddx - ddy * ddy );

    return num / ( speed * speed2 * speed2 );
  }

  void SplineVec::dump_table( ostream_type & stream, integer const num_points ) const
  {
    vector<real_type> vals;
    stream << 's';
    for ( integer i = 0; i < m_dim; ++i ) stream << '\t' << i;
    stream << '\n';
    for ( integer j = 0; j < num_points; ++j )
    {
      real_type const s = x_min() + ( ( x_max() - x_min() ) * j ) / ( num_points - 1 );
      this->eval( s, vals );
      stream << s;
      for ( integer i = 0; i < m_dim; ++i ) stream << '\t' << vals[i];
      stream << '\n';
    }
  }

}  // namespace Splines

//
// EOF: SplineVec.cc
//
