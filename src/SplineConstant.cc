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
#include <iomanip>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

/**
 *
 */

namespace Splines {

  using namespace std; // load standard namspace

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::reserve_external(
    integer      n,
    real_type *& p_x,
    real_type *& p_y
  ) {
    if ( !m_external_alloc ) m_baseValue.free();
    m_npts           = 0;
    m_npts_reserved  = n;
    m_external_alloc = true;
    m_X              = p_x;
    m_Y              = p_y;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::reserve( integer n ) {
    if ( m_external_alloc && n <= m_npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      m_baseValue.allocate( size_t(2*n) );
      m_npts_reserved  = n;
      m_external_alloc = false;
      m_X              = m_baseValue( size_t(n) );
      m_Y              = m_baseValue( size_t(n) );
    }
    initLastInterval();
    m_npts = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Evalute spline value at `x`
  real_type
  ConstantSpline::operator () ( real_type x ) const {
    return m_Y[this->search(x)];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Evalute spline value at `x`
  real_type
  ConstantSpline::id_eval( integer ni, real_type ) const {
    return m_Y[ni];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::build(
    real_type const x[], integer incx,
    real_type const y[], integer incy,
    integer n
  ) {
    reserve( n );
    for ( size_t i = 0; i   < size_t(n); ++i ) m_X[i] = x[i*size_t(incx)];
    for ( size_t i = 0; i+1 < size_t(n); ++i ) m_Y[i] = y[i*size_t(incy)];
    m_npts = n;
    build();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::clear() {
    if ( !m_external_alloc ) m_baseValue.free();
    m_npts = m_npts_reserved = 0;
    m_external_alloc = false;
    m_X = m_Y = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::writeToStream( ostream_type & s ) const {
    size_t nseg = size_t(m_npts > 0 ? m_npts - 1 : 0);
    for ( size_t i = 0; i < nseg; ++i )
      fmt::print( s,
        "segment N. {:4} X:[{},{}] Y:{}\n",
        i,  m_X[i], m_X[i+1], m_Y[i]
      );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  ConstantSpline::coeffs( real_type cfs[], real_type nodes[], bool ) const {
    size_t nseg = size_t(m_npts > 0 ? m_npts - 1 : 0);
    for ( size_t i = 0; i < nseg; ++i ) {
      nodes[i] = m_X[i];
      cfs[i]   = m_Y[i];
    }
    return 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  ConstantSpline::order() const
  { return 1; }

  using GenericContainerNamespace::GC_VEC_REAL;
  using GenericContainerNamespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string msg = fmt::format("ConstantSpline[{}]::setup( gc ):", m_name );
    UTILS_ASSERT(
      gc.exists("xdata"),
      "{} missing `xdata` field!\n", msg
    );
    UTILS_ASSERT(
      gc.exists("ydata"),
      "{} missing `ydata` field!\n", msg
    );

    GenericContainer const & gc_x = gc("xdata");
    GenericContainer const & gc_y = gc("ydata");

    vec_real_type x, y;
    {
      std::string ff = fmt::format( "{}, field `xdata'", msg );
      gc_x.copyto_vec_real ( x, ff.c_str() );
    }
    {
      std::string ff = fmt::format( "{}, field `ydata'", msg );
      gc_y.copyto_vec_real ( y, ff.c_str() );
    }
    this->build( x, y );
  }

}
