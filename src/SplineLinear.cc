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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namspace
#endif

namespace Splines {

  real_type
  LinearSpline::id_eval( integer i, real_type x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0]        ) return m_Y[0];
      if ( x >= m_X[m_npts-1] ) return m_Y[m_npts-1];
    }
    real_type s = (x-m_X[i])/(m_X[i+1] - m_X[i]);
    return (1-s)*m_Y[i] + s * m_Y[i+1];
  }

  real_type
  LinearSpline::operator () ( real_type x ) const {
    integer idx = this->search( x ); // eval idx can modify x
    return this->id_eval( idx, x );
  }

  real_type
  LinearSpline::id_D( integer i, real_type x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    return ( m_Y[i+1] - m_Y[i] ) / ( m_X[i+1] - m_X[i] );
  }

  real_type
  LinearSpline::D( real_type x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    integer idx = this->search( x ); // eval idx can modify x
    return this->id_D( idx, x );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Use externally allocated memory for `npts` points
  void
  LinearSpline::reserve_external(
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
  LinearSpline::reserve( integer n ) {
    if ( m_external_alloc && n <= m_npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      m_baseValue.reallocate( size_t(2*n) );
      m_npts_reserved  = n;
      m_external_alloc = false;
      m_X              = m_baseValue( size_t(n) );
      m_Y              = m_baseValue( size_t(n) );
    }
    initLastInterval();
    m_npts = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::clear(void) {
    if ( !m_external_alloc ) m_baseValue.free();
    m_npts = m_npts_reserved = 0;
    m_external_alloc = false;
    m_X = m_Y = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::writeToStream( ostream_type & s ) const {
    integer nseg = m_npts > 0 ? m_npts - 1 : 0;
    for ( integer i = 0; i < nseg; ++i )
      fmt::print( s,
        "segment N.{:4} X:[{},{}] Y:[{},{}] slope: {}\n",
        i, m_X[i], m_X[i+1], m_Y[i], m_Y[i+1],
        (m_Y[i+1]-m_Y[i])/(m_X[i+1]-m_X[i])
      );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  LinearSpline::coeffs(
    real_type * const cfs,
    real_type * const nodes,
    bool              transpose
  ) const {
    integer n = m_npts > 0 ? m_npts-1 : 0;
    for ( integer i = 0; i < n; ++i ) {
      nodes[i] = m_X[i];
      real_type a = m_Y[i];
      real_type b = (m_Y[i+1]-m_Y[i])/(m_X[i+1]-m_X[i]);
      if ( transpose ) {
        cfs[2*i+1] = a;
        cfs[2*i+0] = b;
      } else {
        cfs[n+i] = a;
        cfs[i]   = b;
      }
    }
    return 2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  LinearSpline::order( ) const { return 2; }

  using GC_namespace::GC_VEC_REAL;
  using GC_namespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string msg = fmt::format("LinearSpline[{}]::setup( gc ):", m_name );
    UTILS_ASSERT( gc.exists("xdata"), "{} missing `xdata` field!\n", msg );
    UTILS_ASSERT( gc.exists("ydata"), "{} missing `ydata` field!\n", msg );

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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::y_min_max(
    integer   & i_min_pos,
    real_type & x_min_pos,
    real_type & y_min,
    integer   & i_max_pos,
    real_type & x_max_pos,
    real_type & y_max
  ) const {
    UTILS_ASSERT0(
      m_npts > 0, "LinearSpline::y_min_max() empty spline!"
    );
    // find max min alongh the nodes
    i_min_pos = i_max_pos = 0;
    x_min_pos = x_max_pos = m_X[0];
    y_min = y_max = m_Y[0];
    for ( integer i = 1; i < m_npts; ++i ) {
      real_type const & P1 = m_Y[i];
      if ( P1 > y_max ) {
        y_max     = P1;
        x_max_pos = m_X[i];
        i_max_pos = i;
      } else if ( P1 < y_min ) {
        y_min     = P1;
        x_min_pos = m_X[i];
        i_min_pos = i;
      }
    }
  }
}
