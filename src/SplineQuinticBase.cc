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
      m_baseValue.allocate( size_t(4*n) );
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
    real_type cfs[],
    real_type nodes[],
    bool      transpose
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
      s << "segment N." << setw(4) << i
        << " X:["      << m_X[i]   << ", " << m_X[i+1]
        << "] Y:["     << m_Y[i]   << ", " << m_Y[i+1]
        << "] Yp:["    << m_Yp[i]  << ", " << m_Yp[i+1]
        << "] Ypp:["   << m_Ypp[i] << ", " << m_Ypp[i+1]
        << "] slope: " << (m_Y[i+1]-m_Y[i])/(m_X[i+1]-m_X[i])
        << '\n';
  }

}
