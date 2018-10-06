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

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wc++98-compat"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

/**
 * 
 */

namespace Splines {

  using namespace std; // load standard namspace

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::reserve_external( integer      n,
                                       real_type *& p_X,
                                       real_type *& p_Y,
                                       real_type *& p_Yp,
                                       real_type *& p_Ypp ) {
    if ( !_external_alloc ) baseValue.free();
    npts            = 0;
    npts_reserved   = n;
    _external_alloc = true;
    X   = p_X;
    Y   = p_Y;
    Yp  = p_Yp;
    Ypp = p_Ypp;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::reserve( integer n ) {
    if ( _external_alloc && n <= npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      baseValue.allocate( size_t(4*n) );
      npts_reserved   = n;
      _external_alloc = false;
      X   = baseValue( size_t(n) );
      Y   = baseValue( size_t(n) );
      Yp  = baseValue( size_t(n) );
      Ypp = baseValue( size_t(n) );
    }
    npts = lastInterval = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::clear(void) {
    if ( !_external_alloc ) baseValue.free();
    npts = npts_reserved = 0;
    _external_alloc = false;
    X = Y = Yp = Ypp = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::operator () ( real_type x ) const {
    size_t i = size_t(Spline::search( x ));
    Hermite5( x-X[i], X[i+1]-X[i], base );
    return base[0] * Y[i]    +
           base[1] * Y[i+1]  +
           base[2] * Yp[i]   +
           base[3] * Yp[i+1] +
           base[4] * Ypp[i]  +
           base[5] * Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::D( real_type x ) const {
    size_t i = size_t(Spline::search( x ));
    Hermite5_D( x-X[i], X[i+1]-X[i], base_D );
    return base_D[0] * Y[i]    +
           base_D[1] * Y[i+1]  +
           base_D[2] * Yp[i]   +
           base_D[3] * Yp[i+1] +
           base_D[4] * Ypp[i]  +
           base_D[5] * Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DD( real_type x ) const {
    size_t i = size_t(Spline::search( x ));
    Hermite5_DD( x-X[i], X[i+1]-X[i], base_DD );
    return base_DD[0] * Y[i]    +
           base_DD[1] * Y[i+1]  +
           base_DD[2] * Yp[i]   +
           base_DD[3] * Yp[i+1] +
           base_DD[4] * Ypp[i]  +
           base_DD[5] * Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDD( real_type x ) const {
    size_t i = size_t(Spline::search( x ));
    Hermite5_DDD( x-X[i], X[i+1]-X[i], base_DDD );
    return base_DDD[0] * Y[i]    +
           base_DDD[1] * Y[i+1]  +
           base_DDD[2] * Yp[i]   +
           base_DDD[3] * Yp[i+1] +
           base_DDD[4] * Ypp[i]  +
           base_DDD[5] * Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDDD( real_type x ) const {
    size_t i = size_t(Spline::search( x ));
    Hermite5_DDDD( x-X[i], X[i+1]-X[i], base_DDDD );
    return base_DDDD[0] * Y[i]    +
           base_DDDD[1] * Y[i+1]  +
           base_DDDD[2] * Yp[i]   +
           base_DDDD[3] * Yp[i+1] +
           base_DDDD[4] * Ypp[i]  +
           base_DDDD[5] * Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDDDD( real_type x ) const {
    size_t i = size_t(Spline::search( x ));
    Hermite5_DDDDD( x-X[i], X[i+1]-X[i], base_DDDDD );
    return base_DDDDD[0] * Y[i]    +
           base_DDDDD[1] * Y[i+1]  +
           base_DDDDD[2] * Yp[i]   +
           base_DDDDD[3] * Yp[i+1] +
           base_DDDDD[4] * Ypp[i]  +
           base_DDDDD[5] * Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  QuinticSplineBase::coeffs(
    real_type cfs[],
    real_type nodes[],
    bool      transpose
  ) const {
    size_t n = size_t(npts > 0 ? npts-1 : 0);
    for ( size_t i = 0; i < n; ++i ) {
      nodes[i] = X[i];
      real_type H = X[i+1]-X[i];
      real_type a = Y[i];
      real_type b = Yp[i];
      real_type c = Ypp[i]/2;
      real_type d = ((10*(Y[i+1]-Y[i])/H-6*Yp[i]-4*Yp[i+1])/H-1.5*Ypp[i]+0.5*Ypp[i+1])/H;
      real_type e = ((15*(Y[i]-Y[i+1])/H+8*Yp[i]+7*Yp[i+1])/H+1.5*Ypp[i]-Ypp[i+1])/(H*H);
      real_type f = ((6*(Y[i+1]-Y[i])/H-3*(Yp[i]+Yp[i+1]))/H-0.5*Ypp[i]+0.5*Ypp[i+1])/(H*H*H);
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
    QuinticSplineBase::reserve(S.npts);
    npts = S.npts;
    std::copy( S.X,   S.X+npts,   X   );
    std::copy( S.Y,   S.Y+npts,   Y   );
    std::copy( S.Yp,  S.Yp+npts,  Yp  );
    std::copy( S.Ypp, S.Ypp+npts, Ypp );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::writeToStream( ostream_type & s ) const {
    size_t nseg = size_t(npts > 0 ? npts - 1 : 0);
    for ( size_t i = 0; i < nseg; ++i )
      s << "segment N." << setw(4) << i
        << " X:["      << X[i]   << ", " << X[i+1]
        << "] Y:["     << Y[i]   << ", " << Y[i+1]
        << "] Yp:["    << Yp[i]  << ", " << Yp[i+1]
        << "] Ypp:["   << Ypp[i] << ", " << Ypp[i+1]
        << "] slope: " << (Y[i+1]-Y[i])/(X[i+1]-X[i])
        << '\n';
  }

}
