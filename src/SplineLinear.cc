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
#endif

/**
 * 
 */

namespace Splines {

  using namespace std; // load standard namspace

  real_type
  LinearSpline::id_eval( integer i, real_type x ) const {
    real_type s = (x-m_X[i])/(m_X[i+1] - m_X[i]);
    return (1-s)*m_Y[i] + s * m_Y[i+1];
  }

  real_type
  LinearSpline::operator () ( real_type x ) const {
    integer idx = this->search( x ); // eval idx can modify x
    return this->id_eval( idx, x );
  }

  real_type
  LinearSpline::id_D( integer i, real_type ) const {
    return ( m_Y[i+1] - m_Y[i] ) / ( m_X[i+1] - m_X[i] );
  }

  real_type
  LinearSpline::D( real_type x ) const {
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
      s << "segment N." << setw(4) << i
        << " X:[ "      << m_X[i] << ", " << m_X[i+1]
        << " ] Y:[ "    << m_Y[i] << ", " << m_Y[i+1]
        << " ] slope: " << (m_Y[i+1]-m_Y[i])/(m_X[i+1]-m_X[i])
        << '\n';
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  LinearSpline::coeffs( real_type cfs[], real_type nodes[], bool transpose ) const {
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

  using GenericContainerNamespace::GC_VEC_REAL;
  using GenericContainerNamespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    SPLINE_ASSERT(
      gc.exists("xdata"),
      "LinearSpline[" << m_name << "]::setup missing `xdata` field!"
    )
    SPLINE_ASSERT(
      gc.exists("ydata"),
      "LinearSpline[" << m_name << "]::setup missing `ydata` field!"
    )

    GenericContainer const & gc_x = gc("xdata");
    GenericContainer const & gc_y = gc("ydata");

    vec_real_type x, y;
    {
      std::ostringstream ost;
      ost << "LinearSpline[" << m_name << "]::setup, field `xdata'";
      gc_x.copyto_vec_real ( x, ost.str().c_str() );
    }
    {
      std::ostringstream ost;
      ost << "LinearSpline[" << m_name << "]::setup, field `ydata'";
      gc_y.copyto_vec_real ( y, ost.str().c_str() );
    }
    this->build( x, y );
  }

}
