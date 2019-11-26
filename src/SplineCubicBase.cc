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

  using namespace std; // load standard namespace

  // build spline without computation

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::build(
    real_type const x[],  integer incx,
    real_type const y[],  integer incy,
    real_type const yp[], integer incyp,
    integer n
  ) {
    this->reserve( n );
    for ( size_t i = 0; i < size_t(n); ++i ) {
      this->X[i]  = x[i*size_t(incx)];
      this->Y[i]  = y[i*size_t(incy)];
      this->Yp[i] = yp[i*size_t(incyp)];
    }
    this->npts = n;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::clear(void) {
    if ( !this->_external_alloc ) this->baseValue.free();
    this->npts = this->npts_reserved = 0;
    this->_external_alloc = false;
    this->X = this->Y = this->Yp = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::reserve( integer n ) {
    if ( _external_alloc && n <= npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      this->npts_reserved = n;
      this->baseValue.allocate( size_t(3*n) );
      this->X  = this->baseValue( size_t(n) );
      this->Y  = this->baseValue( size_t(n) );
      this->Yp = this->baseValue( size_t(n) );
      this->_external_alloc = false;
    }
    integer & lastInterval = lastInterval_by_thread[std::this_thread::get_id()];
    this->npts = lastInterval = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::reserve_external(
    integer      n,
    real_type *& p_x,
    real_type *& p_y,
    real_type *& p_dy
  ) {
    this->npts_reserved = n;
    this->X    = p_x;
    this->Y    = p_y;
    this->Yp   = p_dy;
    integer & lastInterval = lastInterval_by_thread[std::this_thread::get_id()];
    this->npts = lastInterval = 0;
    this->_external_alloc = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::operator () ( real_type x ) const {
    size_t i = size_t(Spline::search( x ));
    Hermite3( x-this->X[i], this->X[i+1]-this->X[i], this->base );
    return this->base[0] * this->Y[i]   +
           this->base[1] * this->Y[i+1] +
           this->base[2] * this->Yp[i]  +
           this->base[3] * this->Yp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::D( real_type x ) const {
    size_t i = size_t(Spline::search( x ));
    Hermite3_D( x-this->X[i], this->X[i+1]-this->X[i], this->base_D );
    return this->base_D[0] * this->Y[i]   +
           this->base_D[1] * this->Y[i+1] +
           this->base_D[2] * this->Yp[i]  +
           this->base_D[3] * this->Yp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::DD( real_type x ) const {
    size_t i = size_t(Spline::search( x ));
    Hermite3_DD( x-this->X[i], this->X[i+1]-this->X[i], this->base_DD );
    return this->base_DD[0] * this->Y[i]   +
           this->base_DD[1] * this->Y[i+1] +
           this->base_DD[2] * this->Yp[i]  +
           this->base_DD[3] * this->Yp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CubicSplineBase::DDD( real_type x ) const {
    size_t i = size_t(Spline::search( x ));
    Hermite3_DDD( x-this->X[i], this->X[i+1]-this->X[i], this->base_DDD );
    return this->base_DDD[0] * this->Y[i]   +
           this->base_DDD[1] * this->Y[i+1] +
           this->base_DDD[2] * this->Yp[i]  +
           this->base_DDD[3] * this->Yp[i+1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  CubicSplineBase::coeffs(
    real_type cfs[],
    real_type nodes[],
    bool      transpose
  ) const {
    size_t n = size_t(this->npts > 0 ? this->npts-1 : 0);
    for ( size_t i = 0; i < n; ++i ) {
      nodes[i] = this->X[i];
      real_type H  = this->X[i+1]-this->X[i];
      real_type DY = (this->Y[i+1]-this->Y[i])/H;
      real_type a = this->Y[i];
      real_type b = this->Yp[i];
      real_type c = (3*DY-2*this->Yp[i]-this->Yp[i+1])/H;
      real_type d = (this->Yp[i+1]+this->Yp[i]-2*DY)/(H*H);
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
    return 4;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  CubicSplineBase::order() const
  { return 4; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Implementation
  void
  CubicSplineBase::copySpline( CubicSplineBase const & S ) {
    CubicSplineBase::reserve(S.npts);
    this->npts = S.npts;
    std::copy( S.X,  S.X+npts,  this->X  );
    std::copy( S.Y,  S.Y+npts,  this->Y  );
    std::copy( S.Yp, S.Yp+npts, this->Yp );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! change X-range of the spline
  void
  CubicSplineBase::setRange( real_type xmin, real_type xmax ) {
    Spline::setRange( xmin, xmax );
    real_type recS = ( this->X[npts-1] - this->X[0] ) / (xmax - xmin);
    real_type * iy = this->Y;
    while ( iy < this->Y + this->npts ) *iy++ *= recS;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSplineBase::writeToStream( ostream_type & s ) const {
    size_t nseg = size_t( this->npts > 0 ? this->npts - 1 : 0 );
    for ( size_t i = 0; i < nseg; ++i )
      s << "segment N." << setw(4) << i
        << " X:[" << this->X[i] << ", " << this->X[i+1]
        << "] Y:[" << this->Y[i] << ", " << this->Y[i+1]
        << "] Yp:[" << this->Yp[i] << ", " << this->Yp[i+1]
        << "] slope: " << (this->Y[i+1]-this->Y[i])/(this->X[i+1]-this->X[i])
        << '\n';
  }

}
