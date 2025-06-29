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

using namespace std; // load standard namspace

namespace Splines {

  using std::copy_n;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::reserve_external(
    integer const n,
    real_type * & p_x,
    real_type * & p_y,
    real_type * & p_Yp,
    real_type * & p_Ypp
  ) {
    if ( !m_external_alloc ) m_base_quintic.free();
    m_npts           = 0;
    m_npts_reserved  = n;
    m_external_alloc = true;
    m_X              = p_x;
    m_Y              = p_y;
    m_Yp             = p_Yp;
    m_Ypp            = p_Ypp;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::reserve( integer const npts ) {
    if ( m_external_alloc && npts <= m_npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      m_base_quintic.reallocate( 4*npts );
      m_npts_reserved  = npts;
      m_external_alloc = false;
      m_X              = m_base_quintic( npts );
      m_Y              = m_base_quintic( npts );
      m_Yp             = m_base_quintic( npts );
      m_Ypp            = m_base_quintic( npts );
    }
    m_npts = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::clear() {
    if ( !m_external_alloc ) m_base_quintic.free();
    m_npts = m_npts_reserved = 0;
    m_external_alloc = false;
    m_X = m_Y = m_Yp = m_Ypp = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_eval( integer const ni, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0]        ) return m_Y[0];
      if ( x >= m_X[m_npts-1] ) return m_Y[m_npts-1];
    }
    real_type base[6];
    real_type const x0 { m_X[ni] };
    real_type const H  { m_X[ni+1] - x0 };
    Hermite5( x-x0, H, base );
    return base[0] * m_Y[ni]   + base[1] * m_Y[ni+1]  +
           base[2] * m_Yp[ni]  + base[3] * m_Yp[ni+1] +
           base[4] * m_Ypp[ni] + base[5] * m_Ypp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::eval( real_type const x ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    return this->id_eval( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_D( integer const ni, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_D[6];
    real_type const x0 { m_X[ni] };
    real_type const H  { m_X[ni+1] - x0 };
    Hermite5_D( x-x0, H, base_D );
    return base_D[0] * m_Y[ni]   + base_D[1] * m_Y[ni+1]  +
           base_D[2] * m_Yp[ni]  + base_D[3] * m_Yp[ni+1] +
           base_D[4] * m_Ypp[ni] + base_D[5] * m_Ypp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::D( real_type const x ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    return this->id_D( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_DD( integer const ni, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_DD[6];
    real_type const x0 = m_X[ni];
    real_type const H  = m_X[ni+1] - x0;
    Hermite5_DD( x-x0, H, base_DD );
    return base_DD[0] * m_Y[ni]   + base_DD[1] * m_Y[ni+1]  +
           base_DD[2] * m_Yp[ni]  + base_DD[3] * m_Yp[ni+1] +
           base_DD[4] * m_Ypp[ni] + base_DD[5] * m_Ypp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DD( real_type const x ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    return this->id_DD( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_DDD( integer const i, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_DDD[6];
    real_type const x0 { m_X[i] };
    real_type const H  { m_X[i+1] - x0 };
    Hermite5_DDD( x-x0, H, base_DDD );
    return base_DDD[0] * m_Y[i]   + base_DDD[1] * m_Y[i+1]  +
           base_DDD[2] * m_Yp[i]  + base_DDD[3] * m_Yp[i+1] +
           base_DDD[4] * m_Ypp[i] + base_DDD[5] * m_Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDD( real_type const x ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    return this->id_DDD( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_DDDD( integer const ni,  real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_DDDD[6];
    real_type const x0 { m_X[ni] };
    real_type const H  { m_X[ni+1] - x0 };
    Hermite5_DDDD( x-x0, H, base_DDDD );
    return base_DDDD[0] * m_Y[ni]   + base_DDDD[1] * m_Y[ni+1]  +
           base_DDDD[2] * m_Yp[ni]  + base_DDDD[3] * m_Yp[ni+1] +
           base_DDDD[4] * m_Ypp[ni] + base_DDDD[5] * m_Ypp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDDD( real_type const x ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    return this->id_DDDD( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_DDDDD( integer const ni, real_type const x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_DDDDD[6];
    real_type const x0 { m_X[ni] };
    real_type const H  { m_X[ni+1] - x0 };
    Hermite5_DDDDD( x-x0, H, base_DDDDD );
    return base_DDDDD[0] * m_Y[ni]   + base_DDDDD[1] * m_Y[ni+1]  +
           base_DDDDD[2] * m_Yp[ni]  + base_DDDDD[3] * m_Yp[ni+1] +
           base_DDDDD[4] * m_Ypp[ni] + base_DDDDD[5] * m_Ypp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDDDD( real_type const x ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    return this->id_DDDDD( res.first, res.second );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::D( real_type const x, real_type dd[2] ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    integer   const ni { res.first  };
    real_type const X  { res.second };
    real_type base[6], base_D[6];
    real_type dx{ X - m_X[ni] };
    real_type DX{ m_X[ni+1]-m_X[ni] };
    Hermite5   ( dx, DX, base    );
    Hermite5_D ( dx, DX, base_D  );

    dd[0] = base[0] * m_Y[ni]   + base[1] * m_Y[ni+1]  +
            base[2] * m_Yp[ni]  + base[3] * m_Yp[ni+1] +
            base[4] * m_Ypp[ni] + base[5] * m_Ypp[ni+1];

    dd[1] = base_D[0] * m_Y[ni]   + base_D[1] * m_Y[ni+1]  +
            base_D[2] * m_Yp[ni]  + base_D[3] * m_Yp[ni+1] +
            base_D[4] * m_Ypp[ni] + base_D[5] * m_Ypp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::DD( real_type const x, real_type dd[3] ) const {
    std::pair<integer,real_type> res(0,x);
    m_search.find( res );
    integer   const ni { res.first  };
    real_type const X  { res.second };
    real_type base[6], base_D[6], base_DD[6];
    real_type dx{ X - m_X[ni] };
    real_type DX{ m_X[ni+1]-m_X[ni] };
    Hermite5    ( dx, DX, base    );
    Hermite5_D  ( dx, DX, base_D  );
    Hermite5_DD ( dx, DX, base_DD );

    dd[0] = base[0] * m_Y[ni]   + base[1] * m_Y[ni+1]  +
            base[2] * m_Yp[ni]  + base[3] * m_Yp[ni+1] +
            base[4] * m_Ypp[ni] + base[5] * m_Ypp[ni+1];

    dd[1] = base_D[0] * m_Y[ni]   + base_D[1] * m_Y[ni+1]  +
            base_D[2] * m_Yp[ni]  + base_D[3] * m_Yp[ni+1] +
            base_D[4] * m_Ypp[ni] + base_D[5] * m_Ypp[ni+1];

    dd[2] = base_DD[0] * m_Y[ni]   + base_DD[1] * m_Y[ni+1]  +
            base_DD[2] * m_Yp[ni]  + base_DD[3] * m_Yp[ni+1] +
            base_DD[4] * m_Ypp[ni] + base_DD[5] * m_Ypp[ni+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifdef AUTODIFF_SUPPORT
  autodiff::dual1st
  QuinticSplineBase::eval( autodiff::dual1st const & x ) const {
    using autodiff::dual1st;
    using autodiff::detail::val;
    real_type dd[2];
    D( val(x), dd );
    dual1st res { dd[0] };
    res.grad = dd[1] * x.grad;
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  autodiff::dual2nd
  QuinticSplineBase::eval( autodiff::dual2nd const & x ) const {
    using autodiff::dual2nd;
    using autodiff::detail::val;
    real_type dd[3], xg{ val(x.grad) };
    DD( val(x), dd );
    dual2nd res { dd[0] };
    res.grad      = dd[1] * xg;
    res.grad.grad = dd[1] * x.grad.grad + dd[2] * (xg*xg);
    return res;
  }
  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  QuinticSplineBase::coeffs(
    real_type  cfs[],
    real_type  nodes[],
    bool const transpose
  ) const {

    UTILS_ASSERT( m_npts >= 2, "QuinticSplineBase, npts={} must be >= 2\n", m_npts );

    integer const n{ m_npts-1 };

    for ( integer i{0}; i < n; ++i ) {
      real_type const H{ m_X[i+1]-m_X[i] };
      real_type a, b, c, d, e, f;
      Hermite5_to_poly(
        H,
        m_Y[i], m_Y[i+1],
        m_Yp[i], m_Yp[i+1],
        m_Ypp[i], m_Ypp[i+1],
        a, b, c, d, e, f
      );
      if ( transpose ) {
        cfs[6*i+5] = a;
        cfs[6*i+4] = b;
        cfs[6*i+3] = c;
        cfs[6*i+2] = d;
        cfs[6*i+1] = e;
        cfs[6*i+0] = f;
      } else {
        cfs[i+5*n] = a;
        cfs[i+4*n] = b;
        cfs[i+3*n] = c;
        cfs[i+2*n] = d;
        cfs[i+1*n] = e;
        cfs[i+0*n] = f;
      }
    }
    copy_n( m_X, m_npts, nodes );
    return 6;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  QuinticSplineBase::order( ) const { return 6; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Implementation
  void
  QuinticSplineBase::copy_spline( QuinticSplineBase const & S ) {
    QuinticSplineBase::reserve(S.m_npts);
    m_npts = S.m_npts;
    copy_n( S.m_X,   m_npts, m_X   );
    copy_n( S.m_Y,   m_npts, m_Y   );
    copy_n( S.m_Yp,  m_npts, m_Yp  );
    copy_n( S.m_Ypp, m_npts, m_Ypp );
    copy_flags( S );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::write_to_stream( ostream_type & s ) const {
    integer const nseg{ m_npts > 0 ? m_npts - 1 : 0 };
    for ( integer i{0}; i < nseg; ++i )
      fmt::print( s,
        "segment N.{:4} X:[{:.5},{:.5}] Y:[{:.5},{:.5}] Yp:[{:.5},{:.5}] Ypp:[{:.5},{:.5}] slope: {:.5}\n",
        i,
        m_X[i],   m_X[i+1],
        m_Y[i],   m_Y[i+1],
        m_Yp[i],  m_Yp[i+1],
        m_Ypp[i], m_Ypp[i+1],
        (m_Y[i+1]-m_Y[i])/(m_X[i+1]-m_X[i])
      );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::y_min_max(
    integer   & i_min_pos,
    real_type & x_min_pos,
    real_type & y_min,
    integer   & i_max_pos,
    real_type & x_max_pos,
    real_type & y_max
  ) const {
    UTILS_ASSERT( m_npts > 0, "QuinticSplineBase[{}]::y_min_max() empty spline!", m_name );
    // find max min alongh the nodes
    i_min_pos = i_max_pos = 0;
    x_min_pos = x_max_pos = m_X[0];
    y_min = y_max = m_Y[0];
    PolynomialRoots::Quartic q;
    for ( integer i{1}; i < m_npts; ++i ) {
      real_type const & X0   { m_X[i-1]   };
      real_type const & X1   { m_X[i]     };
      real_type const & P0   { m_Y[i-1]   };
      real_type const & P1   { m_Y[i]     };
      real_type const & DP0  { m_Yp[i-1]  };
      real_type const & DP1  { m_Yp[i]    };
      real_type const & DDP0 { m_Ypp[i-1] };
      real_type const & DDP1 { m_Ypp[i]   };
      real_type const   H    { X1 - X0    };
      real_type A, B, C, D, E, F;
      Hermite5_to_poly( H, P0, P1, DP0, DP1, DDP0, DDP1, A, B, C, D, E, F );
      q.setup( 5*A, 4*B, 3*C, 2*D, E );
      real_type r[4];
      integer const nr{ q.getRootsInOpenRange( 0, H, r ) };
      for ( integer j{0}; j < nr; ++j ) {
        real_type const rr = r[j];
        real_type yy = (((((A*rr)+B)*rr+C)*rr+D)*rr+E)*rr+F;
        if      ( yy > y_max ) { y_max = yy; x_max_pos = X0+rr; i_max_pos = i; }
        else if ( yy < y_min ) { y_min = yy; x_min_pos = X0+rr; i_min_pos = i; }
      }
      if      ( P1 > y_max ) { y_max = P1; x_max_pos = X1; i_max_pos = i; }
      else if ( P1 < y_min ) { y_min = P1; x_min_pos = X1; i_min_pos = i; }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::y_min_max(
    vector<integer>   & i_min_pos,
    vector<real_type> & x_min_pos,
    vector<real_type> & y_min,
    vector<integer>   & i_max_pos,
    vector<real_type> & x_max_pos,
    vector<real_type> & y_max
  ) const {
    constexpr real_type epsi{1e-8};
    i_min_pos.clear();
    i_max_pos.clear();
    x_min_pos.clear();
    x_max_pos.clear();
    y_min.clear();
    y_max.clear();
    UTILS_ASSERT( m_npts > 0, "QuinticSplineBase[{}]::y_min_max() empty spline!", m_name );
    // find max min alongh the nodes
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
    PolynomialRoots::Quartic q;
    for ( integer i{1}; i < m_npts; ++i ) {
      real_type const & X0   { m_X[i-1]   };
      real_type const & X1   { m_X[i]     };
      real_type const & P0   { m_Y[i-1]   };
      real_type const & P1   { m_Y[i]     };
      real_type const & DP0  { m_Yp[i-1]  };
      real_type const & DP1  { m_Yp[i]    };
      real_type const & DDP0 { m_Ypp[i-1] };
      real_type const & DDP1 { m_Ypp[i]   };
      real_type const   H    { X1 - X0    };
      real_type A, B, C, D, E, F;
      Hermite5_to_poly( H, P0, P1, DP0, DP1, DDP0, DDP1, A, B, C, D, E, F );
      q.setup( 5*A, 4*B, 3*C, 2*D, E );
      real_type r[4];
      integer const nr = q.get_roots_in_open_range( 0, H, r );
      for ( integer j{0}; j < nr; ++j ) {
        real_type const rr  = r[j];
        real_type       yy  = (((((A*rr)+B)*rr+C)*rr+D)*rr+E)*rr+F;
        real_type       ddy = (((20*A*rr)+12*B)*rr+6*C)*rr+2*D;
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
      real_type const & X2   { m_X[i+1]   };
      real_type const & P2   { m_Y[i+1]   };
      real_type const & DP2  { m_Yp[i+1]  };
      real_type const & DDP2 { m_Ypp[i+1] };
      real_type A1, B1, C1, D1, E1, F1;
      Hermite5_to_poly( X2-X1, P1, P2, DP1, DP2, DDP1, DDP2, A1, B1, C1, D1, E1, F1 );
      real_type DD = (((20*A*H)+12*B)*H+6*C)*H+2*D;
      if ( DD >= 0 && D1 >= 0 ) {
        y_min.emplace_back(P1);
        x_min_pos.emplace_back(X1);
        i_min_pos.emplace_back(i);
      } else if ( DD <= 0 && D1 <= 0 ) {
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
