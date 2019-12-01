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
#include <limits>
#include <cmath>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wimplicit-fallthrough"
#endif

/**
 * 
 */

namespace Splines {

  using std::abs;
  using std::sqrt;

  #ifndef SPLINES_DO_NOT_USE_GENERIC_CONTAINER
  using GenericContainerNamespace::mat_real_type;
  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! spline constructor
  SplineVec::SplineVec( string const & name )
  : _name(name)
  , baseValue(name+"_values")
  , basePointer(name+"_pointers")
  , _dim(0)
  , _npts(0)
  , _check_range(true)
  , _X(nullptr)
  , _Y(nullptr)
  , _Yp(nullptr)
  {
    std::lock_guard<std::mutex> lck(lastInterval_mutex);
    lastInterval_by_thread[std::this_thread::get_id()] = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! spline destructor
  SplineVec::~SplineVec() {
    baseValue.free();
    basePointer.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::info( ostream_type & s ) const {
    s << "SplineVec[" << name() << "] n.points = "
      << _npts << " dim = " << _dim << '\n';
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::dump_table( ostream_type & stream, integer num_points ) const {
    vector<real_type> vals;
    stream << 's';
    for ( integer i = 0; i < _dim; ++i ) stream << '\t' << i;
    stream << '\n';
    for ( integer j = 0; j < num_points; ++j ) {
      real_type s = xMin() + ((xMax()-xMin())*j)/(num_points-1);
      this->eval( s, vals );
      stream << s;
      for ( integer i = 0; i < _dim; ++i )
        stream << '\t' << vals[size_t(i)];
      stream << '\n';
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::allocate( integer dim, integer npts ) {

    SPLINE_ASSERT(
      dim > 0, "SplineVec::build expected positive dim = " << dim
    )
    SPLINE_ASSERT(
      npts > 1, "SplineVec::build expected npts = " << npts << " greather than 1"
    )
    _dim  = dim;
    _npts = npts;

    baseValue   . allocate( size_t((2*dim+1)*npts) );
    basePointer . allocate( size_t(2*dim) );

    _Y  = basePointer(size_t(_dim));
    _Yp = basePointer(size_t(_dim));
    _X  = baseValue(size_t(_npts));

    for ( size_t spl = 0; spl < size_t(_dim); ++spl ) {
      _Y[size_t(spl)]  = baseValue(size_t(_npts));
      _Yp[size_t(spl)] = baseValue(size_t(_npts));
    }

    baseValue   . must_be_empty( "SplineVec::build, baseValue" );
    basePointer . must_be_empty( "SplineVec::build, basePointer" );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setup(
    integer           dim,
    integer           npts,
    real_type const * Y[]
  ) {
    allocate( dim, npts );
    for ( size_t spl = 0; spl < size_t(_dim); ++spl )
      std::copy( Y[spl], Y[spl]+_npts, _Y[spl] );
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
    for ( size_t spl = 0; spl < size_t(_dim); ++spl )
      for ( size_t j = 0; j < size_t(_npts); ++j )
        _Y[spl][j] = Y[spl+j*size_t(ldY)];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setKnots( real_type const X[] ) {
    std::copy( X, X+_npts, _X );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void
  SplineVec::computeChords() {
    size_t nn = size_t(_npts-1);
    switch ( _dim ) {
    case 2:
      for ( size_t j = 0; j < nn; ++j ) {
        real_type dx = _Y[0][j+1] - _Y[0][j];
        real_type dy = _Y[1][j+1] - _Y[1][j];
        _X[j] = hypot( dx, dy );
      }
      break;
    case 3:
      for ( size_t j = 0; j < nn; ++j ) {
        real_type dx = _Y[0][j+1] - _Y[0][j];
        real_type dy = _Y[1][j+1] - _Y[1][j];
        real_type dz = _Y[2][j+1] - _Y[2][j];
        _X[j] = hypot( dx, hypot( dy, dz ) );
      }
      break;
    default:
      for ( size_t j = 0; j < nn; ++j ) {
        real_type l = 0;
        for ( size_t k = 0; k < size_t(_dim); ++k ) {
          real_type d = _Y[k][j+1] - _Y[k][j];
          l += d*d;
        }
        _X[j] = sqrt(l);
      }
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void
  SplineVec::setKnotsChordLength() {
    computeChords();
    size_t    nn  = size_t(_npts-1);
    real_type acc = 0;
    for ( size_t j = 0; j <= nn; ++j ) {
      real_type l = _X[j];
      _X[j] = acc;
      acc += l;
    }
    for ( size_t j = 1; j < nn; ++j ) _X[j] /= acc;
    _X[nn] = 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setKnotsCentripetal() {
    computeChords();
    size_t    nn  = size_t(this->_npts-1);
    real_type acc = 0;
    for ( size_t j = 0; j <= nn; ++j ) {
      real_type l = sqrt(this->_X[j]);
      this->_X[j] = acc;
      acc += l;
    }
    for ( size_t j = 1; j < nn; ++j ) this->_X[j] /= acc;
    this->_X[nn] = 1;
  }

  void
  SplineVec::CatmullRom() {
    size_t n = size_t(_npts-1);
    size_t d = size_t(_dim);
    real_type l1, l2, ll, a, b;
    for ( size_t j = 1; j < n; ++j ) {
      l1 = this->_X[j]   - this->_X[j-1];
      l2 = this->_X[j+1] - this->_X[j];
      ll = l1+l2;
      a  = (l2/l1)/ll;
      b  = (l1/l2)/ll;
      for ( size_t k = 0; k < d; ++k )
        this->_Yp[k][j] = a*( this->_Y[k][j] - this->_Y[k][j-1] ) +
                          b*( this->_Y[k][j+1] - this->_Y[k][j] );
    }

    l1 = this->_X[1] - this->_X[0];
    l2 = this->_X[2] -this-> _X[1];
    ll = l1+l2;
    a  = ll/(l1*l2);
    b  = (l1/l2)/ll;
    for ( size_t k = 0; k < d; ++k )
      this->_Yp[k][0] = a*( this->_Y[k][1] - this->_Y[k][0] ) -
                        b*( this->_Y[k][2] - this->_Y[k][0] );

    l1 = this->_X[n]   - this->_X[n-1];
    l2 = this->_X[n-1] - this->_X[n-2];
    ll = l1+l2;
    a  = ll/(l1*l2);
    b  = (l1/l2)/ll;
    for ( size_t k = 0; k < d; ++k )
      this->_Yp[k][n] = b*( this->_Y[k][n-2] - this->_Y[k][n] ) -
                        a*( this->_Y[k][n-1] - this->_Y[k][n] );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::operator () ( real_type x, integer j ) const {
    size_t i = size_t(search( x ));
    real_type base[4];
    Hermite3( x-this->_X[i], this->_X[i+1]-this->_X[i], base );
    return base[0] * this->_Y[size_t(j)][i]   +
           base[1] * this->_Y[size_t(j)][i+1] +
           base[2] * this->_Yp[size_t(j)][i]  +
           base[3] * this->_Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::D( real_type x, integer j ) const {
    size_t i = size_t(search( x ));
    real_type base_D[4];
    Hermite3_D( x-this->_X[i], this->_X[i+1]-this->_X[i], base_D );
    return base_D[0] * this->_Y[size_t(j)][i]   +
           base_D[1] * this->_Y[size_t(j)][i+1] +
           base_D[2] * this->_Yp[size_t(j)][i]  +
           base_D[3] * this->_Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::DD( real_type x, integer j ) const {
    size_t i = size_t(search( x ));
    real_type base_DD[4];
    Hermite3_DD( x-this->_X[i], this->_X[i+1]-this->_X[i], base_DD );
    return base_DD[0] * this->_Y[size_t(j)][i]   +
           base_DD[1] * this->_Y[size_t(j)][i+1] +
           base_DD[2] * this->_Yp[size_t(j)][i]  +
           base_DD[3] * this->_Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::DDD( real_type x, integer j ) const {
    size_t i = size_t(search( x ));
    real_type base_DDD[4];
    Hermite3_DDD( x-this->_X[i], this->_X[i+1]-this->_X[i], base_DDD );
    return base_DDD[0] * this->_Y[size_t(j)][i]   +
           base_DDD[1] * this->_Y[size_t(j)][i+1] +
           base_DDD[2] * this->_Yp[size_t(j)][i]  +
           base_DDD[3] * this->_Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval(
    real_type x,
    real_type vals[],
    integer   inc
  ) const {
    size_t i = size_t(search( x ));
    real_type base[4];
    Hermite3( x-this->_X[i], this->_X[i+1]-this->_X[i], base );
    real_type * v = vals;
    for ( size_t j = 0; j < size_t(this->_dim); ++j, v += inc )
      *v = base[0] * this->_Y[j][i]   +
           base[1] * this->_Y[j][i+1] +
           base[2] * this->_Yp[j][i]  +
           base[3] * this->_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D(
    real_type x,
    real_type vals[],
    integer   inc
  ) const {
    size_t i = size_t(search( x ));
    real_type base_D[4];
    Hermite3_D( x-this->_X[i], this->_X[i+1]-this->_X[i], base_D );
    real_type * v = vals;
    for ( size_t j = 0; j < size_t(this->_dim); ++j, v += inc )
      *v = base_D[0] * this->_Y[j][i]   +
           base_D[1] * this->_Y[j][i+1] +
           base_D[2] * this->_Yp[j][i]  +
           base_D[3] * this->_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD(
    real_type x,
    real_type vals[],
    integer   inc
  ) const {
    size_t i = size_t(search( x ));
    real_type base_DD[4];
    Hermite3_DD( x-this->_X[i], this->_X[i+1]-this->_X[i], base_DD );
    real_type * v = vals;
    for ( size_t j = 0; j < size_t(this->_dim); ++j, v += inc )
      *v = base_DD[0] * this->_Y[j][i]   +
           base_DD[1] * this->_Y[j][i+1] +
           base_DD[2] * this->_Yp[j][i]  +
           base_DD[3] * this->_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD(
    real_type x,
    real_type vals[],
    integer   inc
  ) const {
    size_t i = size_t(search( x ));
    real_type base_DDD[4];
    Hermite3_DDD( x-this->_X[i], this->_X[i+1]-this->_X[i], base_DDD );
    real_type * v = vals;
    for ( size_t j = 0; j < size_t(this->_dim); ++j, v += inc )
      *v = base_DDD[0] * this->_Y[j][i]   +
           base_DDD[1] * this->_Y[j][i+1] +
           base_DDD[2] * this->_Yp[j][i]  +
           base_DDD[3] * this->_Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(this->_dim));
    eval( x, &vals.front(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(this->_dim));
    eval_D( x, &vals.front(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(this->_dim));
    eval_DD( x, &vals.front(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(this->_dim));
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

  #ifndef SPLINES_DO_NOT_USE_GENERIC_CONTAINER
  /*!
   | Evaluate at `x` and fill a GenericContainer
  \*/
  void
  SplineVec::eval(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type & m = vals.set_mat_real( unsigned(this->_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { eval( *px++, v, 1 ); v += this->_dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type & m = vals.set_mat_real( unsigned(this->_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { eval_D( *px++, v, 1 ); v += this->_dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type & m = vals.set_mat_real( unsigned(this->_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { eval_DD( *px++, v, 1 ); v += this->_dim; }

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD(
    vec_real_type const & x,
    GenericContainer    & vals
  ) const {
    mat_real_type & m = vals.set_mat_real( unsigned(this->_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { this->eval_DDD( *px++, v, 1 ); v += this->_dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setup( GenericContainer const & ) {
    //allocate( integer dim, integer npts );
    // DA COMPLETARE
  }

  #endif

}
