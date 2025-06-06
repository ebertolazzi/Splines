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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"

#include <limits>
#include <cmath>
#include <set>

namespace Splines {

  using std::abs;
  using std::sqrt;
  using std::copy_n;

  using GC_namespace::mat_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! spline constructor
  SplineVec::SplineVec( string_view const name )
  : m_name(name)
  , m_mem( fmt::format( "SplineVec[{}]::m_mem", name ) )
  , m_mem_p( fmt::format( "SplineVec[{}]::m_mem_p", name ) )
  {
    m_search.setup( &m_name, &m_npts, &m_X, &m_curve_is_closed, &m_curve_can_extend );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! spline destructor
  SplineVec::~SplineVec() {
    m_mem.free();
    m_mem_p.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  SplineVec::info() const {
    return fmt::format( "SplineVec[{}] n.points={}  dim={}", name(), m_npts, m_dim );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::dump_table( ostream_type & stream, integer const num_points ) const {
    vector<real_type> vals;
    stream << 's';
    for ( integer i{0}; i < m_dim; ++i ) stream << '\t' << i;
    stream << '\n';
    for ( integer j{0}; j < num_points; ++j ) {
      real_type const s = x_min() + ((x_max()-x_min())*j)/(num_points-1);
      this->eval( s, vals );
      stream << s;
      for ( integer i{0}; i < m_dim; ++i )
        stream << '\t' << vals[i];
      stream << '\n';
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  SplineVec::allocate( integer const dim, integer const npts ) {

    UTILS_ASSERT( dim  > 0, "SplineVec[{}]::build expected positive dim = {}\n", m_name, dim );
    UTILS_ASSERT( npts > 1, "SplineVec[{}]::build expected npts = {} greather than 1\n", m_name, npts );
    m_dim  = dim;
    m_npts = npts;

    m_mem.reallocate( (2*dim+1)*npts );
    m_mem_p.reallocate( 2*dim );

    m_Y  = m_mem_p ( m_dim  );
    m_Yp = m_mem_p ( m_dim  );
    m_X  = m_mem   ( m_npts );

    for ( integer spl{0}; spl < m_dim; ++spl ) {
      m_Y[spl]  = m_mem( m_npts );
      m_Yp[spl] = m_mem( m_npts );
    }

    m_mem.must_be_empty( "SplineVec::build, baseValue" );
    m_mem_p.must_be_empty( "SplineVec::build, basePointer" );
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setup(
    integer   const         dim,
    integer   const         npts,
    real_type const * const Y[]
  ) {
    allocate( dim, npts );
    for ( integer spl{0}; spl < m_dim; ++spl )
      copy_n( Y[spl], m_npts, m_Y[spl] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setup(
    integer   const dim,
    integer   const npts,
    real_type const Y[],
    integer   const ldY
  ) {
    allocate( dim, npts );
    for ( integer spl{0}; spl < m_dim; ++spl )
      for ( integer j{0}; j < m_npts; ++j )
        m_Y[spl][j] = Y[spl+j*ldY];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::deep_copy_to( SplineVec & S ) const {
    S.setup( m_dim, m_npts, m_Y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::set_knots( real_type const X[] ) {
    copy_n( X, m_npts, m_X );
    m_search.must_reset();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  SplineVec::compute_chords() {
    integer const nn{ m_npts-1 };
    switch ( m_dim ) {
    case 2:
      for ( integer j{0}; j < nn; ++j ) {
        real_type const dx{ m_Y[0][j+1] - m_Y[0][j] };
        real_type const dy{ m_Y[1][j+1] - m_Y[1][j] };
        m_X[j] = hypot( dx, dy );
      }
      break;
    case 3:
      for ( integer j{0}; j < nn; ++j ) {
        real_type const dx{ m_Y[0][j+1] - m_Y[0][j] };
        real_type const dy{ m_Y[1][j+1] - m_Y[1][j] };
        real_type const dz{ m_Y[2][j+1] - m_Y[2][j] };
        m_X[j] = hypot( dx, hypot( dy, dz ) );
      }
      break;
    default:
      for ( integer j{0}; j < nn; ++j ) {
        real_type l{0};
        for ( integer k{0}; k < m_dim; ++k ) {
          real_type const d{ m_Y[k][j+1] - m_Y[k][j] };
          l += d*d;
        }
        m_X[j] = sqrt(l);
      }
      break;
    }
    m_search.must_reset();
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void
  SplineVec::set_knots_chord_length() {
    compute_chords();
    integer const nn{ m_npts-1 };
    real_type acc{0};
    for ( integer j{0}; j <= nn; ++j ) {
      real_type const l{ m_X[j] };
      m_X[j] = acc;
      acc += l;
    }
    for ( integer j{1}; j < nn; ++j ) m_X[j] /= acc;
    m_X[nn] = 1;
    m_search.must_reset();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::set_knots_centripetal() {
    compute_chords();
    integer const nn{ m_npts-1 };
    real_type acc{0};
    for ( integer j{0}; j <= nn; ++j ) {
      real_type const l{ sqrt(m_X[j]) };
      m_X[j] = acc;
      acc += l;
    }
    for ( integer j{1}; j < nn; ++j ) m_X[j] /= acc;
    m_X[nn] = 1;
    m_search.must_reset();
  }

  void
  SplineVec::catmull_rom() {

    UTILS_ASSERT( m_npts >= 2, "catmull_rom, npts={} must be >= 2\n", m_npts );

    integer const n{ m_npts-1 };
    integer const d{ m_dim    };

    real_type l1, l2, ll, a, b;
    for ( integer j{1}; j < n; ++j ) {
      l1 = m_X[j]   - m_X[j-1];
      l2 = m_X[j+1] - m_X[j];
      ll = l1+l2;
      a  = (l2/l1)/ll;
      b  = (l1/l2)/ll;
      for ( integer k{0}; k < d; ++k )
        m_Yp[k][j] = a*( m_Y[k][j]   - m_Y[k][j-1] ) +
                     b*( m_Y[k][j+1] - m_Y[k][j] );
    }

    l1 = m_X[1] - m_X[0];
    l2 = m_X[2] - m_X[1];
    ll = l1+l2;
    a  = ll/(l1*l2);
    b  = (l1/l2)/ll;
    for ( integer k{0}; k < d; ++k )
      m_Yp[k][0] = a*( m_Y[k][1] - m_Y[k][0] ) -
                   b*( m_Y[k][2] - m_Y[k][0] );

    l1 = m_X[n]   - m_X[n-1];
    l2 = m_X[n-1] - m_X[n-2];
    ll = l1+l2;
    a  = ll/(l1*l2);
    b  = (l1/l2)/ll;
    for ( integer k{0}; k < d; ++k )
      m_Yp[k][n] = b*( m_Y[k][n-2] - m_Y[k][n] ) -
                   a*( m_Y[k][n-1] - m_Y[k][n] );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::eval( real_type const x, integer const j ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    real_type base[4];
    integer const i{res.first};
    Hermite3( res.second-m_X[i], m_X[i+1]-m_X[i], base );
    return base[0] * m_Y [j][i]   +
           base[1] * m_Y [j][i+1] +
           base[2] * m_Yp[j][i]   +
           base[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::D( real_type const x, integer const j ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    real_type base_D[4];
    integer const i{res.first};
    Hermite3_D( res.second-m_X[i], m_X[i+1]-m_X[i], base_D );
    return base_D[0] * m_Y [j][i]   +
           base_D[1] * m_Y [j][i+1] +
           base_D[2] * m_Yp[j][i]  +
           base_D[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::DD( real_type const x, integer const j ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    real_type base_DD[4];
    integer const i{res.first};
    Hermite3_DD( res.second-m_X[i], m_X[i+1]-m_X[i], base_DD );
    return base_DD[0] * m_Y [j][i]   +
           base_DD[1] * m_Y [j][i+1] +
           base_DD[2] * m_Yp[j][i]  +
           base_DD[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::DDD( real_type const x, integer const j ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    real_type base_DDD[4];
    integer const i{res.first};
    Hermite3_DDD( res.second-m_X[i], m_X[i+1]-m_X[i], base_DDD );
    return base_DDD[0] * m_Y [j][i]   +
           base_DDD[1] * m_Y [j][i+1] +
           base_DDD[2] * m_Yp[j][i]  +
           base_DDD[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval(
    real_type const x,
    real_type       vals[],
    integer   const inc
  ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    real_type base[4];
    integer const i{res.first};
    Hermite3( res.second-m_X[i], m_X[i+1]-m_X[i], base );
    real_type * v{vals};
    for ( integer j{0}; j < m_dim; ++j, v += inc )
      *v = base[0] * m_Y[j][i]   +
           base[1] * m_Y[j][i+1] +
           base[2] * m_Yp[j][i]  +
           base[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D(
    real_type const x,
    real_type       vals[],
    integer   const inc
  ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    real_type base_D[4];
    integer const i{res.first};
    Hermite3_D( res.second-m_X[i], m_X[i+1]-m_X[i], base_D );
    real_type * v{vals};
    for ( integer j = 0; j < m_dim; ++j, v += inc )
      *v = base_D[0] * m_Y[j][i]   +
           base_D[1] * m_Y[j][i+1] +
           base_D[2] * m_Yp[j][i]  +
           base_D[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD(
    real_type const x,
    real_type       vals[],
    integer   const inc
  ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    real_type base_DD[4];
    integer const i{res.first};
    Hermite3_DD( res.second-m_X[i], m_X[i+1]-m_X[i], base_DD );
    real_type * v{vals};
    for ( integer j = 0; j < m_dim; ++j, v += inc )
      *v = base_DD[0] * m_Y[j][i]   +
           base_DD[1] * m_Y[j][i+1] +
           base_DD[2] * m_Yp[j][i]  +
           base_DD[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD(
    real_type const x,
    real_type       vals[],
    integer   const inc
  ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    real_type base_DDD[4];
    integer const i{res.first};
    Hermite3_DDD( res.second-m_X[i], m_X[i+1]-m_X[i], base_DDD );
    real_type * v{vals};
    for ( integer j = 0; j < m_dim; ++j, v += inc )
      *v = base_DDD[0] * m_Y[j][i]   +
           base_DDD[1] * m_Y[j][i+1] +
           base_DDD[2] * m_Yp[j][i]  +
           base_DDD[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval( real_type const x, vector<real_type> & vals ) const {
    vals.resize( m_dim );
    eval( x, vals.data(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D( real_type const x, vector<real_type> & vals ) const {
    vals.resize( m_dim );
    eval_D( x, vals.data(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD( real_type const x, vector<real_type> & vals ) const {
    vals.resize( m_dim );
    eval_DD( x, vals.data(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD( real_type const x, vector<real_type> & vals ) const {
    vals.resize( m_dim );
    eval_DDD( x, vals.data(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::curvature( real_type const x ) const {
    real_type D[10], DD[10];
    this->eval_D( x, D, 1 );
    this->eval_DD( x, DD, 1 );
    real_type const x_1 = D[0];
    real_type const x_2 = DD[0];
    real_type const y_1 = D[1];
    real_type const y_2 = DD[1];
    real_type const t4  = x_1 * x_1;
    real_type const t5  = y_1 * y_1;
    real_type const t6  = t4 + t5;
    real_type const t7  = sqrt(t6);
    return (x_1 * y_2 - y_1 * x_2) / ( t6 * t7 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::curvature_D( real_type const x ) const {
    real_type D[10], DD[10], DDD[10];
    this->eval_D( x, D, 1 );
    this->eval_DD( x, DD, 1 );
    this->eval_DDD( x, DDD, 1 );
    real_type const x_1 = D[0];
    real_type const x_2 = DD[0];
    real_type const x_3 = DDD[0];
    real_type const y_1 = D[1];
    real_type const y_2 = DD[1];
    real_type const y_3 = DDD[1];
    real_type const t1  = x_1 * x_1;
    real_type const t9  = x_2 * x_2;
    real_type const t13 = y_1 * y_1;
    real_type const t17 = y_2 * y_2;
    real_type const t26 = t1 + t13;
    real_type const t27 = t26 * t26;
    real_type const t28 = sqrt(t26);
    real_type const aa  = y_3 * x_1 - y_1 * x_3;
    real_type const bb  = 3 * y_2 * x_2;
    return ( t1 * ( aa - bb ) + t13 * ( aa + bb )
           + 3 * x_1 * y_1 * ( t9 - t17 ) ) / ( t28 * t27 );
  }

  //!
  //! Evaluate at `x` and fill a GenericContainer
  //!
  void
  SplineVec::eval(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type   & m  = vals.set_mat_real( static_cast<unsigned>(m_dim), static_cast<unsigned>(x.size()) );
    real_type       * v  = m.data();
    real_type const * px = x.data();
    integer j{static_cast<integer>(x.size())};
    while ( j-- > 0 ) { eval( *px++, v, 1 ); v += m_dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type   & m  = vals.set_mat_real( static_cast<unsigned>(m_dim), static_cast<unsigned>(x.size()) );
    real_type       * v  = m.data();
    real_type const * px = x.data();
    integer j{static_cast<integer>(x.size())};
    while ( j-- > 0 ) { eval_D( *px++, v, 1 ); v += m_dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type   & m  = vals.set_mat_real( static_cast<unsigned>(m_dim), static_cast<unsigned>(x.size()) );
    real_type       * v  = m.data();
    real_type const * px = x.data();
    integer j{static_cast<integer>(x.size())};
    while ( j-- > 0 ) { eval_DD( *px++, v, 1 ); v += m_dim; }

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type   & m  = vals.set_mat_real( static_cast<unsigned>(m_dim), static_cast<unsigned>(x.size()) );
    real_type       * v  = m.data();
    real_type const * px = x.data();
    integer j{static_cast<integer>(x.size())};
    while ( j-- > 0 ) { this->eval_DDD( *px++, v, 1 ); v += m_dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setup( GenericContainer const & gc ) {

    string const where{ fmt::format("SplineVec[{}]::setup( gc ):", m_name ) };

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map(where) ) { keywords.insert(pair.first); }

    GenericContainer const & data{ gc("data",where) };

    mat_real_type Y;
    data.copyto_mat_real(Y,where);

    bool transposed{false};
    gc.get_if_exists("transposed",transposed);

    if ( transposed ) {
      integer const dim { static_cast<integer>(Y.num_rows()) };
      integer const npts{ static_cast<integer>(Y.num_cols()) };
      allocate( dim, npts );
      for ( integer spl{0}; spl < m_dim; ++spl )
        for ( integer j{0}; j < m_npts; ++j )
          m_Y[spl][j] = Y(spl,j);
    } else {
      integer const dim { static_cast<integer>(Y.num_cols()) };
      integer const npts{ static_cast<integer>(Y.num_rows()) };
      allocate( dim, npts );
      for ( integer spl{0}; spl < m_dim; ++spl )
        for ( integer j{0}; j < m_npts; ++j )
          m_Y[spl][j] = Y(j,spl);
    }

    UTILS_WARNING(
      keywords.empty(), "{}: unused keys\n{}\n", where,
      [&keywords]()->string {
        string res;
        for ( auto const & it : keywords ) { res += it; res += ' '; };
        return res;
      }()
    );

    // manca build Yp
  }
}
