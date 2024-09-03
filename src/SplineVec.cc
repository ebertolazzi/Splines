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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Splines.hh"
#include "Utils_fmt.hh"

#include <limits>
#include <cmath>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

namespace Splines {

  using std::abs;
  using std::sqrt;

  using GC_namespace::mat_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! spline constructor
  SplineVec::SplineVec( string const & name )
  : m_name(name)
  , m_mem( fmt::format( "SplineVec[{}]::m_mem", name ) )
  , m_mem_p( fmt::format( "SplineVec[{}]::m_mem_p", name ) )
  {
    init_last_interval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! spline destructor
  SplineVec::~SplineVec() {
    m_mem.free();
    m_mem_p.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::search( std::pair<integer,real_type> & res ) const {
    UTILS_ASSERT( m_npts > 0, "in SplineVec[{}]::search(...), npts == 0!", m_name );
    #ifdef SPLINES_USE_THREADS
    std::unique_lock<std::mutex> lock(m_last_interval_mutex);
    auto id = std::this_thread::get_id();
    auto it = m_last_interval.find(id);
    if ( it == m_last_interval.end() ) {
      it = m_last_interval.insert( {id,std::make_shared<integer>()} ).first;
      *it->second.get() = 0;
    }
    integer & last_interval{ *it->second.get() };
    lock.unlock();
    #else
    integer & last_interval{ m_last_interval };
    #endif
    Utils::search_interval(
      m_npts,
      m_X,
      res.second,
      last_interval,
      m_curve_is_closed,
      m_curve_can_extend
    );
    res.first = last_interval;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  SplineVec::init_last_interval() {
    #ifdef SPLINES_USE_THREADS
    std::unique_lock<std::mutex> lock(m_last_interval_mutex);
    auto id = std::this_thread::get_id();
    auto it = m_last_interval.find(id);
    if ( it == m_last_interval.end() ) it = m_last_interval.insert( {id,std::make_shared<integer>()} ).first;
    integer & last_interval{ *it->second.get() };
    #else
    integer & last_interval{ m_last_interval };
    #endif
    last_interval = 0;
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  SplineVec::info() const {
    return fmt::format( "SplineVec[{}] n.points = {}  dim = {}", name(), m_npts, m_dim );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::dump_table( ostream_type & stream, integer num_points ) const {
    vector<real_type> vals;
    stream << 's';
    for ( integer i = 0; i < m_dim; ++i ) stream << '\t' << i;
    stream << '\n';
    for ( integer j = 0; j < num_points; ++j ) {
      real_type s = x_min() + ((x_max()-x_min())*j)/(num_points-1);
      this->eval( s, vals );
      stream << s;
      for ( integer i = 0; i < m_dim; ++i )
        stream << '\t' << vals[size_t(i)];
      stream << '\n';
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  SplineVec::allocate( integer dim, integer npts ) {

    UTILS_ASSERT(
      dim > 0, "SplineVec[{}]::build expected positive dim = {}\n", m_name, dim
    );
    UTILS_ASSERT(
      npts > 1, "SplineVec[{}]::build expected npts = {} greather than 1\n", m_name, npts
    );
    m_dim  = dim;
    m_npts = npts;

    m_mem.reallocate( size_t((2*dim+1)*npts) );
    m_mem_p.reallocate( size_t(2*dim) );

    m_Y  = m_mem_p(size_t(m_dim));
    m_Yp = m_mem_p(size_t(m_dim));
    m_X  = m_mem(size_t(m_npts));

    for ( size_t spl{0}; spl < size_t(m_dim); ++spl ) {
      m_Y[spl]  = m_mem(size_t(m_npts));
      m_Yp[spl] = m_mem(size_t(m_npts));
    }

    m_mem.must_be_empty( "SplineVec::build, baseValue" );
    m_mem_p.must_be_empty( "SplineVec::build, basePointer" );

    // reset last interval search
    init_last_interval();

  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setup(
    integer           dim,
    integer           npts,
    real_type const * Y[]
  ) {
    allocate( dim, npts );
    for ( size_t spl{0}; spl < size_t(m_dim); ++spl )
      std::copy_n( Y[spl], m_npts, m_Y[spl] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setup(
    integer         dim,
    integer         npts,
    real_type const Y[],
    integer         ldY
  ) {
    allocate( dim, npts );
    for ( size_t spl{0}; spl < size_t(m_dim); ++spl )
      for ( size_t j{0}; j < size_t(m_npts); ++j )
        m_Y[spl][j] = Y[spl+j*size_t(ldY)];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::set_knots( real_type const X[] ) {
    std::copy_n( X, m_npts, m_X );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  SplineVec::compute_chords() {
    size_t nn{ size_t(m_npts-1) };
    switch ( m_dim ) {
    case 2:
      for ( size_t j{0}; j < nn; ++j ) {
        real_type dx{ m_Y[0][j+1] - m_Y[0][j] };
        real_type dy{ m_Y[1][j+1] - m_Y[1][j] };
        m_X[j] = hypot( dx, dy );
      }
      break;
    case 3:
      for ( size_t j{0}; j < nn; ++j ) {
        real_type dx{ m_Y[0][j+1] - m_Y[0][j] };
        real_type dy{ m_Y[1][j+1] - m_Y[1][j] };
        real_type dz{ m_Y[2][j+1] - m_Y[2][j] };
        m_X[j] = hypot( dx, hypot( dy, dz ) );
      }
      break;
    default:
      for ( size_t j{0}; j < nn; ++j ) {
        real_type l{0};
        for ( size_t k{0}; k < size_t(m_dim); ++k ) {
          real_type d{ m_Y[k][j+1] - m_Y[k][j] };
          l += d*d;
        }
        m_X[j] = sqrt(l);
      }
      break;
    }
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void
  SplineVec::set_knots_chord_length() {
    compute_chords();
    size_t nn{ size_t(m_npts-1) };
    real_type acc{0};
    for ( size_t j{0}; j <= nn; ++j ) {
      real_type l{ m_X[j] };
      m_X[j] = acc;
      acc += l;
    }
    for ( size_t j{1}; j < nn; ++j ) m_X[j] /= acc;
    m_X[nn] = 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::set_knots_centripetal() {
    compute_chords();
    size_t nn{ size_t(m_npts-1) };
    real_type acc{0};
    for ( size_t j{0}; j <= nn; ++j ) {
      real_type l{ sqrt(m_X[j]) };
      m_X[j] = acc;
      acc += l;
    }
    for ( size_t j = 1; j < nn; ++j ) m_X[j] /= acc;
    m_X[nn] = 1;
  }

  void
  SplineVec::catmull_rom() {
    size_t n{ size_t(m_npts-1) };
    size_t d{ size_t(m_dim) };
    real_type l1, l2, ll, a, b;
    for ( size_t j{1}; j < n; ++j ) {
      l1 = m_X[j]   - m_X[j-1];
      l2 = m_X[j+1] - m_X[j];
      ll = l1+l2;
      a  = (l2/l1)/ll;
      b  = (l1/l2)/ll;
      for ( size_t k{0}; k < d; ++k )
        m_Yp[k][j] = a*( m_Y[k][j]   - m_Y[k][j-1] ) +
                     b*( m_Y[k][j+1] - m_Y[k][j] );
    }

    l1 = m_X[1] - m_X[0];
    l2 = m_X[2] - m_X[1];
    ll = l1+l2;
    a  = ll/(l1*l2);
    b  = (l1/l2)/ll;
    for ( size_t k{0}; k < d; ++k )
      m_Yp[k][0] = a*( m_Y[k][1] - m_Y[k][0] ) -
                   b*( m_Y[k][2] - m_Y[k][0] );

    l1 = m_X[n]   - m_X[n-1];
    l2 = m_X[n-1] - m_X[n-2];
    ll = l1+l2;
    a  = ll/(l1*l2);
    b  = (l1/l2)/ll;
    for ( size_t k{0}; k < d; ++k )
      m_Yp[k][n] = b*( m_Y[k][n-2] - m_Y[k][n] ) -
                   a*( m_Y[k][n-1] - m_Y[k][n] );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::eval( real_type x, integer j ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    real_type base[4];
    integer i{res.first};
    Hermite3( res.second-m_X[i], m_X[i+1]-m_X[i], base );
    return base[0] * m_Y[size_t(j)][i]   +
           base[1] * m_Y[size_t(j)][i+1] +
           base[2] * m_Yp[size_t(j)][i]  +
           base[3] * m_Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::D( real_type x, integer j ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    real_type base_D[4];
    integer i{res.first};
    Hermite3_D( res.second-m_X[i], m_X[i+1]-m_X[i], base_D );
    return base_D[0] * m_Y[size_t(j)][i]   +
           base_D[1] * m_Y[size_t(j)][i+1] +
           base_D[2] * m_Yp[size_t(j)][i]  +
           base_D[3] * m_Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::DD( real_type x, integer j ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    real_type base_DD[4];
    integer i{res.first};
    Hermite3_DD( res.second-m_X[i], m_X[i+1]-m_X[i], base_DD );
    return base_DD[0] * m_Y[size_t(j)][i]   +
           base_DD[1] * m_Y[size_t(j)][i+1] +
           base_DD[2] * m_Yp[size_t(j)][i]  +
           base_DD[3] * m_Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::DDD( real_type x, integer j ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    real_type base_DDD[4];
    integer i{res.first};
    Hermite3_DDD( res.second-m_X[i], m_X[i+1]-m_X[i], base_DDD );
    return base_DDD[0] * m_Y[size_t(j)][i]   +
           base_DDD[1] * m_Y[size_t(j)][i+1] +
           base_DDD[2] * m_Yp[size_t(j)][i]  +
           base_DDD[3] * m_Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval(
    real_type x,
    real_type vals[],
    integer   inc
  ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    real_type base[4];
    integer i{res.first};
    Hermite3( res.second-m_X[i], m_X[i+1]-m_X[i], base );
    real_type * v{vals};
    for ( size_t j{0}; j < size_t(m_dim); ++j, v += inc )
      *v = base[0] * m_Y[j][i]   +
           base[1] * m_Y[j][i+1] +
           base[2] * m_Yp[j][i]  +
           base[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D(
    real_type x,
    real_type vals[],
    integer   inc
  ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    real_type base_D[4];
    integer i{res.first};
    Hermite3_D( res.second-m_X[i], m_X[i+1]-m_X[i], base_D );
    real_type * v{vals};
    for ( size_t j = 0; j < size_t(m_dim); ++j, v += inc )
      *v = base_D[0] * m_Y[j][i]   +
           base_D[1] * m_Y[j][i+1] +
           base_D[2] * m_Yp[j][i]  +
           base_D[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD(
    real_type x,
    real_type vals[],
    integer   inc
  ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    real_type base_DD[4];
    integer i{res.first};
    Hermite3_DD( res.second-m_X[i], m_X[i+1]-m_X[i], base_DD );
    real_type * v{vals};
    for ( size_t j = 0; j < size_t(m_dim); ++j, v += inc )
      *v = base_DD[0] * m_Y[j][i]   +
           base_DD[1] * m_Y[j][i+1] +
           base_DD[2] * m_Yp[j][i]  +
           base_DD[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD(
    real_type x,
    real_type vals[],
    integer   inc
  ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    real_type base_DDD[4];
    integer i{res.first};
    Hermite3_DDD( res.second-m_X[i], m_X[i+1]-m_X[i], base_DDD );
    real_type * v{vals};
    for ( size_t j = 0; j < size_t(m_dim); ++j, v += inc )
      *v = base_DDD[0] * m_Y[j][i]   +
           base_DDD[1] * m_Y[j][i+1] +
           base_DDD[2] * m_Yp[j][i]  +
           base_DDD[3] * m_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(m_dim));
    eval( x, &vals.front(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(m_dim));
    eval_D( x, &vals.front(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(m_dim));
    eval_DD( x, &vals.front(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(m_dim));
    eval_DDD( x, &vals.front(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::curvature( real_type x ) const {
    real_type D[10], DD[10];
    this->eval_D( x, D, 1 );
    this->eval_DD( x, DD, 1 );
    real_type x_1 = D[0];
    real_type x_2 = DD[0];
    real_type y_1 = D[1];
    real_type y_2 = DD[1];
    real_type t4  = x_1 * x_1;
    real_type t5  = y_1 * y_1;
    real_type t6  = t4 + t5;
    real_type t7  = sqrt(t6);
    return (x_1 * y_2 - y_1 * x_2) / ( t6 * t7 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::curvature_D( real_type x ) const {
    real_type D[10], DD[10], DDD[10];
    this->eval_D( x, D, 1 );
    this->eval_DD( x, DD, 1 );
    this->eval_DDD( x, DDD, 1 );
    real_type x_1 = D[0];
    real_type x_2 = DD[0];
    real_type x_3 = DDD[0];
    real_type y_1 = D[1];
    real_type y_2 = DD[1];
    real_type y_3 = DDD[1];
    real_type t1  = x_1 * x_1;
    real_type t9  = x_2 * x_2;
    real_type t13 = y_1 * y_1;
    real_type t17 = y_2 * y_2;
    real_type t26 = t1 + t13;
    real_type t27 = t26 * t26;
    real_type t28 = sqrt(t26);
    real_type aa  = y_3 * x_1 - y_1 * x_3;
    real_type bb  = 3 * y_2 * x_2;
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
    mat_real_type   & m  = vals.set_mat_real( unsigned(m_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { eval( *px++, v, 1 ); v += m_dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type   & m  = vals.set_mat_real( unsigned(m_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { eval_D( *px++, v, 1 ); v += m_dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type   & m  = vals.set_mat_real( unsigned(m_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { eval_DD( *px++, v, 1 ); v += m_dim; }

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type & m = vals.set_mat_real( unsigned(m_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { this->eval_DDD( *px++, v, 1 ); v += m_dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setup( GenericContainer const & ) {
    //allocate( integer dim, integer npts );
    // DA COMPLETARE
    UTILS_ERROR("SplineVec::setup not yet implemented");
  }

}
