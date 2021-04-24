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
#include "PolynomialRoots.hh"
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::reserve_external(
    integer       n,
    real_type * & p_X,
    real_type * & p_Y,
    real_type * & p_Yp,
    real_type * & p_Ypp
  ) {
    if ( !m_external_alloc ) m_baseValue.free();
    m_npts           = 0;
    m_npts_reserved  = n;
    m_external_alloc = true;
    m_X              = p_X;
    m_Y              = p_Y;
    m_Yp             = p_Yp;
    m_Ypp            = p_Ypp;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::reserve( integer n ) {
    if ( m_external_alloc && n <= m_npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      m_baseValue.reallocate( size_t(4*n) );
      m_npts_reserved  = n;
      m_external_alloc = false;
      m_X              = m_baseValue( size_t(n) );
      m_Y              = m_baseValue( size_t(n) );
      m_Yp             = m_baseValue( size_t(n) );
      m_Ypp            = m_baseValue( size_t(n) );
    }
    initLastInterval();
    m_npts = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::clear(void) {
    if ( !m_external_alloc ) m_baseValue.free();
    m_npts = m_npts_reserved = 0;
    m_external_alloc = false;
    m_X = m_Y = m_Yp = m_Ypp = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_eval( integer i, real_type x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0]        ) return m_Y[0];
      if ( x >= m_X[m_npts-1] ) return m_Y[m_npts-1];
    }
    real_type base[6];
    real_type x0 = m_X[i];
    real_type H  = m_X[i+1] - x0;
    Hermite5( x-x0, H, base );
    return base[0] * m_Y[i]   + base[1] * m_Y[i+1]  +
           base[2] * m_Yp[i]  + base[3] * m_Yp[i+1] +
           base[4] * m_Ypp[i] + base[5] * m_Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::operator () ( real_type x ) const {
    integer idx = this->search( x ); // eval idx can modify x
    return this->id_eval( idx, x );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_D( integer i, real_type x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_D[6];
    real_type x0 = m_X[i];
    real_type H  = m_X[i+1] - x0;
    Hermite5_D( x-x0, H, base_D );
    return base_D[0] * m_Y[i]   + base_D[1] * m_Y[i+1]  +
           base_D[2] * m_Yp[i]  + base_D[3] * m_Yp[i+1] +
           base_D[4] * m_Ypp[i] + base_D[5] * m_Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::D( real_type x ) const {
    integer idx = this->search( x ); // eval idx can modify x
    return this->id_D( idx, x );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_DD( integer i, real_type x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_DD[6];
    real_type x0 = m_X[i];
    real_type H  = m_X[i+1] - x0;
    Hermite5_DD( x-x0, H, base_DD );
    return base_DD[0] * m_Y[i]   + base_DD[1] * m_Y[i+1]  +
           base_DD[2] * m_Yp[i]  + base_DD[3] * m_Yp[i+1] +
           base_DD[4] * m_Ypp[i] + base_DD[5] * m_Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DD( real_type x ) const {
    integer idx = this->search( x ); // eval idx can modify x
    return this->id_DD( idx, x );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_DDD( integer i, real_type x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_DDD[6];
    real_type x0 = m_X[i];
    real_type H  = m_X[i+1] - x0;
    Hermite5_DDD( x-x0, H, base_DDD );
    return base_DDD[0] * m_Y[i]   + base_DDD[1] * m_Y[i+1]  +
           base_DDD[2] * m_Yp[i]  + base_DDD[3] * m_Yp[i+1] +
           base_DDD[4] * m_Ypp[i] + base_DDD[5] * m_Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDD( real_type x ) const {
    integer idx = this->search( x ); // eval idx can modify x
    return this->id_DDD( idx, x );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_DDDD( integer i,  real_type x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_DDDD[6];
    real_type x0 = m_X[i];
    real_type H  = m_X[i+1] - x0;
    Hermite5_DDDD( x-x0, H, base_DDDD );
    return base_DDDD[0] * m_Y[i]   + base_DDDD[1] * m_Y[i+1]  +
           base_DDDD[2] * m_Yp[i]  + base_DDDD[3] * m_Yp[i+1] +
           base_DDDD[4] * m_Ypp[i] + base_DDDD[5] * m_Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDDD( real_type x ) const {
    integer idx = this->search( x ); // eval idx can modify x
    return this->id_DDDD( idx, x );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::id_DDDDD( integer i, real_type x ) const {
    if ( m_curve_can_extend && m_curve_extended_constant ) {
      if ( x <= m_X[0] || x >= m_X[m_npts-1] ) return 0;
    }
    real_type base_DDDDD[6];
    real_type x0 = m_X[i];
    real_type H  = m_X[i+1] - x0;
    Hermite5_DDDDD( x-x0, H, base_DDDDD );
    return base_DDDDD[0] * m_Y[i]   + base_DDDDD[1] * m_Y[i+1]  +
           base_DDDDD[2] * m_Yp[i]  + base_DDDDD[3] * m_Yp[i+1] +
           base_DDDDD[4] * m_Ypp[i] + base_DDDDD[5] * m_Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDDDD( real_type x ) const {
    integer idx = this->search( x ); // eval idx can modify x
    return this->id_DDDDD( idx, x );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  QuinticSplineBase::coeffs(
    real_type * const cfs,
    real_type * const nodes,
    bool              transpose
  ) const {
    size_t n = size_t(m_npts > 0 ? m_npts-1 : 0);
    for ( size_t i = 0; i < n; ++i ) {
      nodes[i] = m_X[i];
      real_type H = m_X[i+1]-m_X[i];
      real_type a = m_Y[i];
      real_type b = m_Yp[i];
      real_type c = m_Ypp[i]/2;
      real_type d = ((10*(m_Y[i+1]-m_Y[i])/H-6*m_Yp[i]-4*m_Yp[i+1])/H-1.5*m_Ypp[i]+0.5*m_Ypp[i+1])/H;
      real_type e = ((15*(m_Y[i]-m_Y[i+1])/H+8*m_Yp[i]+7*m_Yp[i+1])/H+1.5*m_Ypp[i]-m_Ypp[i+1])/(H*H);
      real_type f = ((6*(m_Y[i+1]-m_Y[i])/H-3*(m_Yp[i]+m_Yp[i+1]))/H-0.5*m_Ypp[i]+0.5*m_Ypp[i+1])/(H*H*H);
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
    return 6;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  QuinticSplineBase::order( ) const { return 6; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Implementation
  void
  QuinticSplineBase::copySpline( QuinticSplineBase const & S ) {
    QuinticSplineBase::reserve(S.m_npts);
    m_npts = S.m_npts;
    std::copy_n( S.m_X,   m_npts, m_X   );
    std::copy_n( S.m_Y,   m_npts, m_Y   );
    std::copy_n( S.m_Yp,  m_npts, m_Yp  );
    std::copy_n( S.m_Ypp, m_npts, m_Ypp );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::writeToStream( ostream_type & s ) const {
    size_t nseg = size_t(m_npts > 0 ? m_npts - 1 : 0);
    for ( size_t i = 0; i < nseg; ++i )
      fmt::print( s,
        "segment N.{:4} X:[{},{}] Y:[{},{}] Yp:[{},{}] Ypp:[{},{}] slope: {}\n",
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
    cout << "QuinticSplineBase entra\n" << flush;
    UTILS_ASSERT0(
      m_npts > 0, "QuinticSplineBase::y_min_max() empty spline!"
    );
    cout << "QuinticSplineBase m_X " << m_X << "\n" << flush;
    cout << "QuinticSplineBase m_Y " << m_Y << "\n" << flush;
    cout << "QuinticSplineBase m_Yp " << m_Yp << "\n" << flush;
    cout << "QuinticSplineBase m_Yp " << m_Yp << "\n" << flush;
    // find max min alongh the nodes
    i_min_pos = i_max_pos = 0;
    x_min_pos = x_max_pos = m_X[0];
    y_min = y_max = m_Y[0];
    PolynomialRoots::Quartic q;
    for ( integer i = 1; i < m_npts; ++i ) {
      cout << "QuinticSplineBase i = " << i << '\n' << flush;
      real_type const & X0   = m_X[i-1];
      real_type const & X1   = m_X[i];
      real_type const & P0   = m_Y[i-1];
      real_type const & P1   = m_Y[i];
      real_type const & DP0  = m_Yp[i-1];
      real_type const & DP1  = m_Yp[i];
      real_type const & DDP0 = m_Ypp[i-1];
      real_type const & DDP1 = m_Ypp[i];
      real_type H = X1 - X0;
      real_type A, B, C, D, E, F;
      cout << "QuinticSplineBase i = " << i << '\n' << flush;
      Hermite5toPoly( H, P0, P1, DP0, DP1, DDP0, DDP1, A, B, C, D, E, F );
      cout << "QuinticSplineBase i = " << i << '\n' << flush;
      cout.precision(30);
      cout << "A = " << 5*A << endl;
      cout << "B = " << 4*B << endl;
      cout << "C = " << 3*C << endl;
      cout << "D = " << 2*D << endl;
      cout << "E = " << E << endl;
      q.setup( 5*A, 4*B, 3*C, 2*D, E );
      cout << "QuinticSplineBase i = " << i << '\n' << flush;
      real_type r[4];
      cout << "QuinticSplineBase i = " << i << '\n' << flush;
      integer nr = q.getRootsInOpenRange( 0, H, r );
      cout << "QuinticSplineBase nr = " << nr << '\n' << flush;
      for ( integer j = 0; j < nr; ++j ) {
        real_type rr = r[j];
        cout << "QuinticSplineBase rr = " << rr << '\n' << flush;
        real_type yy = (((((A*rr)+B)*rr+C)*rr+D)*rr+E)*rr+F;
        if      ( yy > y_max ) { y_max = yy; x_max_pos = X0+rr; i_max_pos = i; }
        else if ( yy < y_min ) { y_min = yy; x_min_pos = X0+rr; i_min_pos = i; }
      }
      cout << "QuinticSplineBase i = " << i << '\n' << flush;
      if      ( P1 > y_max ) { y_max = P1; x_max_pos = X1; i_max_pos = i; }
      else if ( P1 < y_min ) { y_min = P1; x_min_pos = X1; i_min_pos = i; }
    }
  }
}
