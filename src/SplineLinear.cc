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

  LinearSpline::LinearSpline( string_view name )
  : Spline(name)
  , m_mem_linear( fmt::format( "LinearSpline[{}]", name ) )
  {
    m_curve_extended_constant = true; // by default linear spline extend constant
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  LinearSpline::id_eval( integer const ni, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0]        ) return m_Y[0];
      if ( x >= m_X[m_npts-1] ) return m_Y[m_npts-1];
    }
    real_type const s = (x-m_X[ni])/(m_X[ni+1] - m_X[ni]);
    return (1-s)*m_Y[ni] + s * m_Y[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  LinearSpline::eval( real_type x ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    return this->id_eval( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  LinearSpline::id_D( integer const i, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    return ( m_Y[i+1] - m_Y[i] ) / ( m_X[i+1] - m_X[i] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  LinearSpline::D( real_type x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    return this->id_D( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Use externally allocated memory for `npts` points
  void
  LinearSpline::reserve_external(
    integer const n,
    real_type *&  p_x,
    real_type *&  p_y
  ) {
    if ( !m_external_alloc ) m_mem_linear.free();
    m_npts           = 0;
    m_npts_reserved  = n;
    m_external_alloc = true;
    m_X              = p_x;
    m_Y              = p_y;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::reserve( integer const npts ) {
    if ( m_external_alloc && npts <= m_npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      m_mem_linear.reallocate( 2*npts );
      m_npts_reserved  = npts;
      m_external_alloc = false;
      m_X              = m_mem_linear( npts );
      m_Y              = m_mem_linear( npts );
    }
    m_npts = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::clear() {
    if ( !m_external_alloc ) m_mem_linear.free();
    m_npts = m_npts_reserved = 0;
    m_external_alloc = false;
    m_X = m_Y = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::write_to_stream( ostream_type & s ) const {
    integer const nseg{ m_npts > 0 ? m_npts - 1 : 0 };
    for ( integer i{0}; i < nseg; ++i )
      fmt::print( s,
        "segment N.{:4} X:[{:.5},{:.5}] Y:[{:.5},{:.5}] slope: {:.5}\n",
        i, m_X[i], m_X[i+1], m_Y[i], m_Y[i+1],
        (m_Y[i+1]-m_Y[i])/(m_X[i+1]-m_X[i])
      );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  LinearSpline::coeffs(
    real_type  cfs[],
    real_type  nodes[],
    bool const transpose
  ) const {

    UTILS_ASSERT( m_npts >= 2, "LinearSpline::coeffs, npts={} must be >= 2\n", m_npts );

    integer const n{ m_npts-1 };
    for ( integer i{0}; i < n; ++i ) {
      real_type const a { m_Y[i] };
      real_type const b { (m_Y[i+1]-m_Y[i])/(m_X[i+1]-m_X[i]) };
      if ( transpose ) {
        cfs[2*i+1] = a;
        cfs[2*i+0] = b;
      } else {
        cfs[n+i] = a;
        cfs[i]   = b;
      }
    }
    copy_n( m_X, m_npts, nodes );
    return 2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  LinearSpline::order( ) const { return 2; }

  using GC_namespace::GC_type;
  using GC_namespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where{ fmt::format("LinearSpline[{}]::setup( gc ):", m_name ) };
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
  LinearSpline::y_min_max(
    integer   & i_min_pos,
    real_type & x_min_pos,
    real_type & y_min,
    integer   & i_max_pos,
    real_type & x_max_pos,
    real_type & y_max
  ) const {
    UTILS_ASSERT( m_npts > 0, "LinearSpline[{}]::y_min_max() empty spline!", m_name );
    // find max min along the nodes
    i_min_pos = i_max_pos = 0;
    x_min_pos = x_max_pos = m_X[0];
    y_min = y_max = m_Y[0];
    for ( integer i{1}; i < m_npts; ++i ) {
      real_type const & P1{ m_Y[i] };
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
  LinearSpline::y_min_max(
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
    UTILS_ASSERT( m_npts > 0, "LinearSpline[{}]::y_min_max() empty spline!", m_name );
    // find max min along the nodes
    for ( integer i{1}; i < m_npts-1; ++i ) {
      real_type const & P0 { m_Y[i-1] };
      real_type const & P1 { m_Y[i]   };
      real_type const & P2 { m_Y[i+1] };
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
