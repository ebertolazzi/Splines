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

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wc++98-compat"
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
  , lastInterval(0)
  {}

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

    SPLINE_ASSERT( dim > 0,
                   "SplineVec::build expected positive dim = " << dim );
    SPLINE_ASSERT( npts > 1,
                   "SplineVec::build expected npts = " << npts <<
                   " greather than 1" );
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
  SplineVec::setup( integer           dim,
                    integer           npts,
                    real_type const * Y[] ) {
    allocate( dim, npts );
    for ( size_t spl = 0; spl < size_t(_dim); ++spl )
      std::copy( Y[spl], Y[spl]+_npts, _Y[spl] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setup( integer         dim,
                    integer         npts,
                    real_type const Y[],
                    integer         ldY ) {
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
    size_t    nn  = size_t(_npts-1);
    real_type acc = 0;
    for ( size_t j = 0; j <= nn; ++j ) {
      real_type l = sqrt(_X[j]);
      _X[j] = acc;
      acc += l;
    }
    for ( size_t j = 1; j < nn; ++j ) _X[j] /= acc;
    _X[nn] = 1;
  }

  void
  SplineVec::CatmullRom() {
    size_t n = size_t(_npts-1);
    size_t d = size_t(_dim);
    real_type l1, l2, ll, a, b;
    for ( size_t j = 1; j < n; ++j ) {
      l1 = _X[j]   - _X[j-1];
      l2 = _X[j+1] - _X[j];
      ll = l1+l2;
      a  = (l2/l1)/ll;
      b  = (l1/l2)/ll;
      for ( size_t k = 0; k < d; ++k )
        _Yp[k][j] = a*( _Y[k][j] - _Y[k][j-1] ) + b*( _Y[k][j+1] - _Y[k][j] );
    }

    l1 = _X[1] - _X[0];
    l2 = _X[2] - _X[1];
    ll = l1+l2;
    a  = ll/(l1*l2);
    b  = (l1/l2)/ll;
    for ( size_t k = 0; k < d; ++k )
      _Yp[k][0] = a*( _Y[k][1] - _Y[k][0] ) - b*( _Y[k][2] - _Y[k][0] );

    l1 = _X[n]   - _X[n-1];
    l2 = _X[n-1] - _X[n-2];
    ll = l1+l2;
    a  = ll/(l1*l2);
    b  = (l1/l2)/ll;
    for ( size_t k = 0; k < d; ++k )
      _Yp[k][n] = b*( _Y[k][n-2] - _Y[k][n] ) - a*( _Y[k][n-1] - _Y[k][n] );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::operator () ( real_type x, integer j ) const {
    size_t i = size_t(search( x ));
    real_type base[4];
    Hermite3( x-_X[i], _X[i+1]-_X[i], base );
    return base[0] * _Y[size_t(j)][i]   +
           base[1] * _Y[size_t(j)][i+1] +
           base[2] * _Yp[size_t(j)][i]  +
           base[3] * _Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::D( real_type x, integer j ) const {
    size_t i = size_t(search( x ));
    real_type base_D[4];
    Hermite3_D( x-_X[i], _X[i+1]-_X[i], base_D );
    return base_D[0] * _Y[size_t(j)][i]   +
           base_D[1] * _Y[size_t(j)][i+1] +
           base_D[2] * _Yp[size_t(j)][i]  +
           base_D[3] * _Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::DD( real_type x, integer j ) const {
    size_t i = size_t(search( x ));
    real_type base_DD[4];
    Hermite3_DD( x-_X[i], _X[i+1]-_X[i], base_DD );
    return base_DD[0] * _Y[size_t(j)][i]   +
           base_DD[1] * _Y[size_t(j)][i+1] +
           base_DD[2] * _Yp[size_t(j)][i]  +
           base_DD[3] * _Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  SplineVec::DDD( real_type x, integer j ) const {
    size_t i = size_t(search( x ));
    real_type base_DDD[4];
    Hermite3_DDD( x-_X[i], _X[i+1]-_X[i], base_DDD );
    return base_DDD[0] * _Y[size_t(j)][i]   +
           base_DDD[1] * _Y[size_t(j)][i+1] +
           base_DDD[2] * _Yp[size_t(j)][i]  +
           base_DDD[3] * _Yp[size_t(j)][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval( real_type x,
                   real_type vals[],
                   integer   inc ) const {
    size_t i = size_t(search( x ));
    real_type base[4];
    Hermite3( x-_X[i], _X[i+1]-_X[i], base );
    real_type * v = vals;
    for ( size_t j = 0; j < size_t(_dim); ++j, v += inc )
      *v = base[0] * _Y[j][i]   +
           base[1] * _Y[j][i+1] +
           base[2] * _Yp[j][i]  +
           base[3] * _Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D( real_type x,
                     real_type vals[],
                     integer   inc ) const {
    size_t i = size_t(search( x ));
    real_type base_D[4];
    Hermite3_D( x-_X[i], _X[i+1]-_X[i], base_D );
    real_type * v = vals;
    for ( size_t j = 0; j < size_t(_dim); ++j, v += inc )
      *v = base_D[0] * _Y[j][i]   +
           base_D[1] * _Y[j][i+1] +
           base_D[2] * _Yp[j][i]  +
           base_D[3] * _Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD( real_type x,
                      real_type vals[],
                      integer   inc ) const {
    size_t i = size_t(search( x ));
    real_type base_DD[4];
    Hermite3_DD( x-_X[i], _X[i+1]-_X[i], base_DD );
    real_type * v = vals;
    for ( size_t j = 0; j < size_t(_dim); ++j, v += inc )
      *v = base_DD[0] * _Y[j][i]   +
           base_DD[1] * _Y[j][i+1] +
           base_DD[2] * _Yp[j][i]  +
           base_DD[3] * _Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD( real_type x,
                       real_type vals[],
                       integer   inc ) const {
    size_t i = size_t(search( x ));
    real_type base_DDD[4];
    Hermite3_DDD( x-_X[i], _X[i+1]-_X[i], base_DDD );
    real_type * v = vals;
    for ( size_t j = 0; j < size_t(_dim); ++j, v += inc )
      *v = base_DDD[0] * _Y[j][i]   +
           base_DDD[1] * _Y[j][i+1] +
           base_DDD[2] * _Yp[j][i]  +
           base_DDD[3] * _Yp[j][i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(_dim));
    eval( x, &vals.front(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(_dim));
    eval_D( x, &vals.front(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(_dim));
    eval_DD( x, &vals.front(), 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD( real_type x, vector<real_type> & vals ) const {
    vals.resize(size_t(_dim));
    eval_DDD( x, &vals.front(), 1 );
  }

  #ifndef SPLINES_DO_NOT_USE_GENERIC_CONTAINER
  /*!
   | Evaluate at `x` and fill a GenericContainer
  \*/
  void
  SplineVec::eval( vec_real_type const & x,
                   GenericContainer    & vals ) const {
    mat_real_type & m = vals.set_mat_real( unsigned(_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { eval( *px++, v, 1 ); v += _dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_D( vec_real_type const & x,
                     GenericContainer    & vals ) const {
    mat_real_type & m = vals.set_mat_real( unsigned(_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { eval_D( *px++, v, 1 ); v += _dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DD( vec_real_type const & x,
                      GenericContainer    & vals ) const {
    mat_real_type & m = vals.set_mat_real( unsigned(_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { eval_DD( *px++, v, 1 ); v += _dim; }

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::eval_DDD( vec_real_type const & x,
                       GenericContainer    & vals ) const {
    mat_real_type & m = vals.set_mat_real( unsigned(_dim), unsigned(x.size()) );
    real_type       * v  = &m.front();
    real_type const * px = &x.front();
    integer j = integer(x.size());
    while ( j-- > 0 ) { eval_DDD( *px++, v, 1 ); v += _dim; }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  SplineVec::setup( GenericContainer const & gc ) {
    //allocate( integer dim, integer npts );
    // DA COMPLETARE
  }

  #endif

}
