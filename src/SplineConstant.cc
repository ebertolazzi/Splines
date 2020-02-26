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
  ConstantSpline::reserve_external(
    integer      n,
    real_type *& p_x,
    real_type *& p_y
  ) {
    if ( !this->_external_alloc ) baseValue.free();
    this->npts            = 0;
    this->npts_reserved   = n;
    this->_external_alloc = true;
    this->X               = p_x;
    this->Y               = p_y;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::reserve( integer n ) {
    if ( this->_external_alloc && n <= this->npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      baseValue.allocate( size_t(2*n) );
      this->npts_reserved   = n;
      this->_external_alloc = false;
      this->X               = baseValue( size_t(n) );
      this->Y               = baseValue( size_t(n) );
    }
    {
      std::lock_guard<std::mutex> lck(lastInterval_mutex);
      lastInterval_by_thread[std::this_thread::get_id()] = 0;
    }
    this->npts = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! Evalute spline value at `x`
  real_type
  ConstantSpline::operator () ( real_type x ) const {
    if ( x < X[0] ) return Y[0];
    if ( this->npts > 0 && x > this->X[this->npts-1] ) return this->Y[this->npts-1];
    return this->Y[this->search(x)];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::build(
    real_type const x[], integer incx,
    real_type const y[], integer incy,
    integer n
  ) {
    reserve( n );
    for ( size_t i = 0; i   < size_t(n); ++i ) this->X[i] = x[i*size_t(incx)];
    for ( size_t i = 0; i+1 < size_t(n); ++i ) this->Y[i] = y[i*size_t(incy)];
    this->npts = n;
    build();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::clear() {
    if ( !_external_alloc ) baseValue.free();
    this->npts = this->npts_reserved = 0;
    this->_external_alloc = false;
    this->X = this->Y = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::writeToStream( ostream_type & s ) const {
    size_t nseg = size_t(this->npts > 0 ? this->npts - 1 : 0);
    for ( size_t i = 0; i < nseg; ++i )
      s << "segment N." << setw(4) << i
        << " X:[ " << X[i] << ", " << X[i+1] << " ] Y:" << Y[i]
        << '\n';
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  ConstantSpline::coeffs( real_type cfs[], real_type nodes[], bool ) const {
    size_t nseg = size_t(this->npts > 0 ? this->npts - 1 : 0);
    for ( size_t i = 0; i < nseg; ++i ) {
      nodes[i] = this->X[i];
      cfs[i]   = this->Y[i];
    }
    return 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  ConstantSpline::order() const
  { return 1; }

  #ifndef SPLINES_DO_NOT_USE_GENERIC_CONTAINER

  using GenericContainerNamespace::GC_VEC_REAL;
  using GenericContainerNamespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ConstantSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    SPLINE_ASSERT(
      gc.exists("xdata"),
      "ConstantSpline[" << this->_name << "]::setup missing `xdata` field!"
    )
    SPLINE_ASSERT(
      gc.exists("ydata"),
      "ConstantSpline[" << this->_name << "]::setup missing `ydata` field!"
    )

    GenericContainer const & gc_x = gc("xdata");
    GenericContainer const & gc_y = gc("ydata");

    vec_real_type x, y;
    {
      std::ostringstream ost;
      ost << "ConstantSpline[" << this->_name << "]::setup, field `xdata'";
      gc_x.copyto_vec_real ( x, ost.str().c_str() );
    }
    {
      std::ostringstream ost;
      ost << "ConstantSpline[" << this->_name << "]::setup, field `ydata'";
      gc_y.copyto_vec_real ( y, ost.str().c_str() );
    }
    this->build( x, y );
  }
  #endif

}
