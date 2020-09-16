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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::reserve_external(
    integer       n,
    real_type * & p_X,
    real_type * & p_Y,
    real_type * & p_Yp,
    real_type * & p_Ypp
  ) {
    if ( !this->_external_alloc ) this->baseValue.free();
    this->npts            = 0;
    this->npts_reserved   = n;
    this->_external_alloc = true;
    this->X               = p_X;
    this->Y               = p_Y;
    this->Yp              = p_Yp;
    this->Ypp             = p_Ypp;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::reserve( integer n ) {
    if ( this->_external_alloc && n <= this->npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      this->baseValue.allocate( size_t(4*n) );
      this->npts_reserved   = n;
      this->_external_alloc = false;
      this->X               = this->baseValue( size_t(n) );
      this->Y               = this->baseValue( size_t(n) );
      this->Yp              = this->baseValue( size_t(n) );
      this->Ypp             = this->baseValue( size_t(n) );
    }
    initLastInterval();
    this->npts = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSplineBase::clear(void) {
    if ( !this->_external_alloc ) this->baseValue.free();
    this->npts = this->npts_reserved = 0;
    this->_external_alloc = false;
    this->X = this->Y = this->Yp = this->Ypp = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::operator () ( real_type x, integer i ) const {
    real_type base[6];
    real_type x0 = this->X[i];
    real_type H  = this->X[i+1] - x0;
    Hermite5( x-x0, H, base );
    return base[0] * this->Y[i]   + base[1] * this->Y[i+1]  +
           base[2] * this->Yp[i]  + base[3] * this->Yp[i+1] +
           base[4] * this->Ypp[i] + base[5] * this->Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::operator () ( real_type x ) const {
    return this->operator () ( x,this->search( x ) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::D( real_type x, integer i ) const {
    real_type base_D[6];
    real_type x0 = this->X[i];
    real_type H  = this->X[i+1] - x0;
    Hermite5_D( x-x0, H, base_D );
    return base_D[0] * this->Y[i]   + base_D[1] * this->Y[i+1]  +
           base_D[2] * this->Yp[i]  + base_D[3] * this->Yp[i+1] +
           base_D[4] * this->Ypp[i] + base_D[5] * this->Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::D( real_type x ) const {
    return this->D( x, this->search( x ) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DD( real_type x, integer i ) const {
    real_type base_DD[6];
    real_type x0 = this->X[i];
    real_type H  = this->X[i+1] - x0;
    Hermite5_DD( x-x0, H, base_DD );
    return base_DD[0] * this->Y[i]   + base_DD[1] * this->Y[i+1]  +
           base_DD[2] * this->Yp[i]  + base_DD[3] * this->Yp[i+1] +
           base_DD[4] * this->Ypp[i] + base_DD[5] * this->Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DD( real_type x ) const {
    return this->DD( x, this->search( x ) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDD( real_type x, integer i ) const {
    real_type base_DDD[6];
    real_type x0 = this->X[i];
    real_type H  = this->X[i+1] - x0;
    Hermite5_DDD( x-x0, H, base_DDD );
    return base_DDD[0] * this->Y[i]   + base_DDD[1] * this->Y[i+1]  +
           base_DDD[2] * this->Yp[i]  + base_DDD[3] * this->Yp[i+1] +
           base_DDD[4] * this->Ypp[i] + base_DDD[5] * this->Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDD( real_type x ) const {
    return this->DDD( x, this->search( x ) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDDD( real_type x, integer i ) const {
    real_type base_DDDD[6];
    real_type x0 = this->X[i];
    real_type H  = this->X[i+1] - x0;
    Hermite5_DDDD( x-x0, H, base_DDDD );
    return base_DDDD[0] * this->Y[i]   + base_DDDD[1] * this->Y[i+1]  +
           base_DDDD[2] * this->Yp[i]  + base_DDDD[3] * this->Yp[i+1] +
           base_DDDD[4] * this->Ypp[i] + base_DDDD[5] * this->Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDDD( real_type x ) const {
    return this->DDDD( x, this->search( x ) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDDDD( real_type x, integer i ) const {
    real_type base_DDDDD[6];
    real_type x0 = this->X[i];
    real_type H  = this->X[i+1] - x0;
    Hermite5_DDDDD( x-x0, H, base_DDDDD );
    return base_DDDDD[0] * this->Y[i]   + base_DDDDD[1] * this->Y[i+1]  +
           base_DDDDD[2] * this->Yp[i]  + base_DDDDD[3] * this->Yp[i+1] +
           base_DDDDD[4] * this->Ypp[i] + base_DDDDD[5] * this->Ypp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  QuinticSplineBase::DDDDD( real_type x ) const {
    return this->DDDDD( x, this->search( x ) );
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
