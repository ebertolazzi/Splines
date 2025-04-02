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
#include "PolynomialRoots.hh"
#include "Utils_fmt.hh"

#include <iomanip>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namspace

namespace Splines {

  using std::copy_n;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  CubicSplineBase::CubicSplineBase( string_view name )
  : Spline(name)
  , m_mem_cubic( fmt::format("CubicSplineBase[{}]::m_mem_cubic", name ) )
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::build(
    real_type const x[],  integer const incx,
    real_type const y[],  integer const incy,
    real_type const yp[], integer const incyp,
    integer   const n
  ) {
    this->reserve( n );
    for ( integer i{0}; i < n; ++i ) {
      m_X[i]  = x[i*incx];
      m_Y[i]  = y[i*incy];
      m_Yp[i] = yp[i*incyp];
    }
    m_npts = n;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::build(
    vector<real_type> const & x,
    vector<real_type> const & y,
    vector<real_type> const & yp
  ) {
    integer N{ static_cast<integer>(x.size()) };
    if ( N > static_cast<integer>(y.size())  ) N = static_cast<integer>(y.size());
    if ( N > static_cast<integer>(yp.size()) ) N = static_cast<integer>(yp.size());
    this->build (
      &x.front(),  1,
      &y.front(),  1,
      &yp.front(), 1,
      N
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::clear() {
    if ( !m_external_alloc ) m_mem_cubic.free();
    m_npts = m_npts_reserved = 0;
    m_external_alloc = false;
    m_X = m_Y = m_Yp = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::reserve( integer const npts ) {
    if ( m_external_alloc && npts <= m_npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      m_npts_reserved = npts;
      m_mem_cubic.reallocate( 3*npts );
      m_X  = m_mem_cubic( npts );
      m_Y  = m_mem_cubic( npts );
      m_Yp = m_mem_cubic( npts );
      m_external_alloc = false;
    }
    init_last_interval();
    m_npts = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::reserve_external(
    integer const n,
    real_type *&  p_x,
    real_type *&  p_y,
    real_type *& p_dy
  ) {
    m_npts_reserved  = n;
    m_X              = p_x;
    m_Y              = p_y;
    m_Yp             = p_dy;
    m_external_alloc = true;
    init_last_interval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::id_eval( integer const ni, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0]        ) return m_Y[0];
      if ( x >= m_X[m_npts-1] ) return m_Y[m_npts-1];
    }
    real_type base[4];
    Hermite3( x-m_X[ni], m_X[ni+1]-m_X[ni], base );
    return base[0] * m_Y[ni]   +
           base[1] * m_Y[ni+1] +
           base[2] * m_Yp[ni]  +
           base[3] * m_Yp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::eval( real_type x ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    return this->id_eval( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::id_D( integer const ni, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_D[4];
    Hermite3_D( x-m_X[ni], m_X[ni+1]-m_X[ni], base_D );
    return base_D[0] * m_Y[ni]   +
           base_D[1] * m_Y[ni+1] +
           base_D[2] * m_Yp[ni]  +
           base_D[3] * m_Yp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::D( real_type x ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    return this->id_D( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::id_DD( integer const ni, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_DD[4];
    Hermite3_DD( x-m_X[ni], m_X[ni+1]-m_X[ni], base_DD );
    return base_DD[0] * m_Y[ni]   +
           base_DD[1] * m_Y[ni+1] +
           base_DD[2] * m_Yp[ni]  +
           base_DD[3] * m_Yp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::DD( real_type x ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    return this->id_DD( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::id_DDD( integer const ni, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_DDD[4];
    Hermite3_DDD( x-m_X[ni], m_X[ni+1]-m_X[ni], base_DDD );
    return base_DDD[0] * m_Y[ni]   +
           base_DDD[1] * m_Y[ni+1] +
           base_DDD[2] * m_Yp[ni]  +
           base_DDD[3] * m_Yp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::DDD( real_type x ) const {
    std::pair<integer,real_type> res(0,x);
    this->search( res );
    return this->id_DDD( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  CubicSplineBase::coeffs(
    real_type  cfs[],
    real_type  nodes[],
    bool const transpose
  ) const {

    UTILS_ASSERT( m_npts >= 2, "CubicSplineBase::coeffs, npts={} must be >= 2\n", m_npts );

    integer const n{ m_npts-1 };
    for ( integer i{0}; i < n; ++i ) {
      real_type const H{ m_X[i+1]-m_X[i] };
      real_type a, b, c, d;

      Hermite3_to_poly( H, m_Y[i], m_Y[i+1], m_Yp[i], m_Yp[i+1], a, b, c, d );

      if ( transpose ) {
        cfs[4*i+3] = a;
        cfs[4*i+2] = b;
        cfs[4*i+1] = c;
        cfs[4*i+0] = d;
      } else {
        cfs[i+3*n] = a;
        cfs[i+2*n] = b;
        cfs[i+1*n] = c;
        cfs[i+0*n] = d;
      }
    }
    copy_n( m_X, m_npts, nodes );
    return 4;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  CubicSplineBase::order() const
  { return 4; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Implementation
  void
  CubicSplineBase::copy_spline( CubicSplineBase const & S ) {
    CubicSplineBase::reserve(S.m_npts);
    m_npts = S.m_npts;
    copy_n( S.m_X,  m_npts, m_X  );
    copy_n( S.m_Y,  m_npts, m_Y  );
    copy_n( S.m_Yp, m_npts, m_Yp );
    copy_flags( S );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! change X-range of the spline
  void
  CubicSplineBase::set_range( real_type const xmin, real_type const xmax ) {
    Spline::set_range( xmin, xmax );
    real_type const recS = ( m_X[m_npts-1] - m_X[0] ) / (xmax - xmin);
    real_type * iy = m_Y;
    while ( iy < m_Y + m_npts ) *iy++ *= recS;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::write_to_stream( ostream_type & s ) const {
    integer const nseg{ m_npts > 0 ? m_npts - 1 : 0 };
    for ( integer i{0}; i < nseg; ++i )
      fmt::print( s,
        "segment N.{:4} X:[{:.5},{:.5}] Y:[{:.5},{:.5}] Yp:[{:.5},{:.5}] slope: {}\n",
        i, m_X[i], m_X[i+1], m_Y[i], m_Y[i+1], m_Yp[i], m_Yp[i+1],
        (m_Y[i+1]-m_Y[i])/(m_X[i+1]-m_X[i])
      );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::y_min_max(
    integer   & i_min_pos,
    real_type & x_min_pos,
    real_type & y_min,
    integer   & i_max_pos,
    real_type & x_max_pos,
    real_type & y_max
  ) const {
    UTILS_ASSERT( m_npts > 0, "CubicSplineBase[{}]::y_min_max() empty spline!", m_name );
    // find max min alongh the nodes
    i_min_pos = i_max_pos = 0;
    x_min_pos = x_max_pos = m_X[0];
    y_min     = y_max     = m_Y[0];
    PolynomialRoots::Quadratic q;
    for ( integer i{1}; i < m_npts; ++i ) {
      real_type const & X0  { m_X[i-1]  };
      real_type const & X1  { m_X[i]    };
      real_type const & P0  { m_Y[i-1]  };
      real_type const & P1  { m_Y[i]    };
      real_type const & DP0 { m_Yp[i-1] };
      real_type const & DP1 { m_Yp[i]   };
      real_type const   H   { X1 - X0   };
      real_type A, B, C, D;
      Hermite3_to_poly( H, P0, P1, DP0, DP1, A, B, C, D );
      q.setup( 3*A, 2*B, C );
      real_type r[2];
      integer const nr{ q.getRootsInOpenRange( 0, H, r ) };
      for ( integer j{0}; j < nr; ++j ) {
        real_type const rr = r[j];
        real_type const yy = (((A*rr)+B)*rr+C)*rr+D;
        if      ( yy > y_max ) { y_max = yy; x_max_pos = X0+rr; i_max_pos = i; }
        else if ( yy < y_min ) { y_min = yy; x_min_pos = X0+rr; i_min_pos = i; }
      }
      if      ( P1 > y_max ) { y_max = P1; x_max_pos = X1; i_max_pos = i; }
      else if ( P1 < y_min ) { y_min = P1; x_min_pos = X1; i_min_pos = i; }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::y_min_max(
    vector<integer>   & i_min_pos,
    vector<real_type> & x_min_pos,
    vector<real_type> & y_min,
    vector<integer>   & i_max_pos,
    vector<real_type> & x_max_pos,
    vector<real_type> & y_max
  ) const {
    constexpr real_type epsi = 1e-8;
    i_min_pos.clear();
    i_max_pos.clear();
    x_min_pos.clear();
    x_max_pos.clear();
    y_min.clear();
    y_max.clear();
    UTILS_ASSERT( m_npts > 0, "CubicSplineBase[{}]::y_min_max() empty spline!", m_name );
    // find max min along the nodes
    if ( m_Yp[0] >= 0 ) {
      y_min.emplace_back(m_Y[0]);
      x_min_pos.emplace_back(m_X[0]);
      i_min_pos.emplace_back(0);
    }
    if ( m_Yp[0] <= 0 ) {
      y_max.emplace_back(m_Y[0]);
      x_max_pos.emplace_back(m_X[0]);
      i_max_pos.emplace_back(0);
    }
    PolynomialRoots::Quadratic q;
    for ( integer i{1}; i < m_npts; ++i ) {
      real_type const & X0  { m_X[i-1]  };
      real_type const & X1  { m_X[i]    };
      real_type const & P0  { m_Y[i-1]  };
      real_type const & P1  { m_Y[i]    };
      real_type const & DP0 { m_Yp[i-1] };
      real_type const & DP1 { m_Yp[i]   };
      real_type const   H   { X1 - X0   };
      real_type A, B, C, D;
      Hermite3_to_poly( H, P0, P1, DP0, DP1, A, B, C, D );
      q.setup( 3*A, 2*B, C );
      real_type r[2];
      integer const nr{ q.getRootsInOpenRange( 0, H, r ) };
      for ( integer j{0}; j < nr; ++j ) {
        real_type const rr  { r[j] };
        real_type const yy  { (((A*rr)+B)*rr+C)*rr+D };
        real_type const ddy { 3*A*rr+B };
        if ( ddy > 0 ) {
          y_min.emplace_back(yy);
          x_min_pos.emplace_back(X0+rr);
          i_min_pos.emplace_back(i);
        } else if ( ddy < 0 ) {
          y_max.emplace_back(yy);
          x_max_pos.emplace_back(X0+rr);
          i_max_pos.emplace_back(i);
        }
      }
      if ( i+1 >= m_npts ) continue;
      if ( abs(DP1) > (m_X[i+1]-m_X[i-1])*epsi ) continue;
      real_type const & X2  { m_X[i+1]  };
      real_type const & P2  { m_Y[i+1]  };
      real_type const & DP2 { m_Yp[i+1] };
      real_type A1, B1, C1, D1;
      Hermite3_to_poly( X2-X1, P1, P2, DP1, DP2, A1, B1, C1, D1 );
      real_type const DD = 2*A*H+B;
      if ( DD >= 0 && B1 >= 0 ) {
        y_min.emplace_back(P1);
        x_min_pos.emplace_back(X1);
        i_min_pos.emplace_back(i);
      } else if ( DD <= 0 && B1 <= 0 ) {
        y_max.emplace_back(P1);
        x_max_pos.emplace_back(X1);
        i_max_pos.emplace_back(i);
      }
    }
    if ( m_Yp[m_npts-1] <= 0 ) {
      y_min.emplace_back(m_Y[m_npts-1]);
      x_min_pos.emplace_back(m_X[m_npts-1]);
      i_min_pos.emplace_back(0);
    }
    if ( m_Yp[m_npts-1] >= 0 ) {
      y_max.emplace_back(m_Y[m_npts-1]);
      x_max_pos.emplace_back(m_X[m_npts-1]);
      i_max_pos.emplace_back(m_npts-1);
    }
  }
}

#endif
