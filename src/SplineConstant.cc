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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

#include "Splines.hh"
#include "Utils_fmt.hh"

#include <iomanip>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namspace
#endif

namespace Splines {

  using std::copy_n;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ConstantSpline::ConstantSpline( string_view name )
  : Spline(name)
  , m_mem_constant( fmt::format("ConstantSpline[{}]",name) )
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::reserve_external(
    integer const n,
    real_type *&  p_x,
    real_type *&  p_y
  ) {
    if ( !m_external_alloc ) m_mem_constant.free();
    m_npts           = 0;
    m_npts_reserved  = n;
    m_external_alloc = true;
    m_X              = p_x;
    m_Y              = p_y;
    init_last_interval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::reserve( integer npts ) {
    if ( m_external_alloc && npts <= m_npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      m_mem_constant.reallocate( 2*npts );
      m_npts_reserved  = npts;
      m_external_alloc = false;
      m_X              = m_mem_constant( npts );
      m_Y              = m_mem_constant( npts );
    }
    init_last_interval();
    m_npts = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Evalute spline value at `x`
  real_type
  ConstantSpline::eval( real_type x ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    return m_Y[res.first];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Evalute spline value at `x`
  real_type
  ConstantSpline::id_eval( integer const ni, real_type ) const {
    return m_Y[ni];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::build(
    real_type const x[], integer const incx,
    real_type const y[], integer const incy,
    integer const n
  ) {
    reserve( n );
    for ( integer i{0}; i   < n; ++i ) m_X[i] = x[i*incx];
    for ( integer i{0}; i+1 < n; ++i ) m_Y[i] = y[i*incy];
    m_npts = n;
    build();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::clear() {
    if ( !m_external_alloc ) m_mem_constant.free();
    m_npts = m_npts_reserved = 0;
    m_external_alloc = false;
    m_X = m_Y = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::write_to_stream( ostream_type & s ) const {
    integer const nseg{ m_npts > 0 ? m_npts - 1 : 0 };
    for ( integer i{0}; i < nseg; ++i )
      fmt::print( s,
        "segment N. {:4} X:[{},{}] Y:{}\n",
        i,  m_X[i], m_X[i+1], m_Y[i]
      );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  ConstantSpline::coeffs(
    real_type cfs[],
    real_type nodes[],
    bool
  ) const {
    integer const nseg{ m_npts > 0 ? m_npts - 1 : 0 };
    copy_n( m_X, m_npts, nodes );
    copy_n( m_Y, nseg,   cfs   );
    return 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  ConstantSpline::order() const
  { return 1; }

  using GC_namespace::GC_type;
  using GC_namespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where{ fmt::format("ConstantSpline[{}]::setup( gc ):", m_name ) };
    GenericContainer const & gc_x{ gc("xdata",where) };
    GenericContainer const & gc_y{ gc("ydata",where) };

    vec_real_type x, y;
    {
      string const ff{ fmt::format( "{}, field `xdata'", where ) };
      gc_x.copyto_vec_real ( x, ff );
    }
    {
      string const ff{ fmt::format( "{}, field `ydata'", where ) };
      gc_y.copyto_vec_real ( y, ff );
    }
    this->build( x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::y_min_max(
    integer   & i_min_pos,
    real_type & x_min_pos,
    real_type & y_min,
    integer   & i_max_pos,
    real_type & x_max_pos,
    real_type & y_max
  ) const {
    UTILS_ASSERT(
      m_npts > 0, "ConstantSpline[{}]::y_min_max() empty spline!", m_name
    );
    // find max min alongh the nodes
    i_min_pos = i_max_pos = 0;
    x_min_pos = x_max_pos = m_X[0];
    y_min = y_max = m_Y[0];
    for ( integer i{1}; i < m_npts-1; ++i ) {
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::y_min_max(
    vector<integer>   & i_min_pos,
    vector<real_type> & x_min_pos,
    vector<real_type> & y_min,
    vector<integer>   & i_max_pos,
    vector<real_type> & x_max_pos,
    vector<real_type> & y_max
  ) const {
    i_min_pos.clear();
    i_max_pos.clear();
    x_min_pos.clear();
    x_max_pos.clear();
    y_min.clear();
    y_max.clear();
    UTILS_ASSERT(
      m_npts > 0, "ConstantSpline[{}]::y_min_max() empty spline!", m_name
    );
    // find max min along the nodes
    for ( integer i{1}; i < m_npts-1; ++i ) {
      real_type const & P0 = m_Y[i-1];
      real_type const & P1 = m_Y[i];
      real_type const & P2 = m_Y[i+1];
      if ( P1 > P0 && P1 > P2 ) {
        y_max.emplace_back(P1);
        x_max_pos.emplace_back(m_X[i]);
        i_max_pos.emplace_back(i);
      } else if ( P1 < P0 && P1 < P2 ) {
        y_min.emplace_back(P1);
        x_min_pos.emplace_back(m_X[i]);
        i_min_pos.emplace_back(i);
      }
    }
  }
}
