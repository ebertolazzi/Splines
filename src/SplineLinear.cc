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

  //! Use externally allocated memory for `npts` points
  void
  LinearSpline::reserve_external(
    integer      n,
    real_type *& p_x,
    real_type *& p_y
  ) {
    if ( !this->_external_alloc ) this->baseValue.free();
    this->npts            = 0;
    this->npts_reserved   = n;
    this->_external_alloc = true;
    this->X               = p_x;
    this->Y               = p_y;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::reserve( integer n ) {
    if ( this->_external_alloc && n <= this->npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      this->baseValue.allocate( size_t(2*n) );
      this->npts_reserved   = n;
      this->_external_alloc = false;
      this->X               = this->baseValue( size_t(n) );
      this->Y               = this->baseValue( size_t(n) );
    }
    {
      std::unique_lock<std::mutex> lck(lastInterval_mutex);
      lastInterval_by_thread[std::this_thread::get_id()] = 0;
    }
    this->npts = 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::clear(void) {
    if ( !this->_external_alloc ) this->baseValue.free();
    this->npts = this->npts_reserved = 0;
    this->_external_alloc = false;
    this->X = this->Y = nullptr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LinearSpline::writeToStream( ostream_type & s ) const {
    integer nseg = this->npts > 0 ? this->npts - 1 : 0;
    for ( integer i = 0; i < nseg; ++i )
      s << "segment N." << setw(4) << i
        << " X:[ " << this->X[i] << ", " << this->X[i+1]
        << " ] Y:[ " << this->Y[i] << ", " << this->Y[i+1]
        << " ] slope: " << (this->Y[i+1]-this->Y[i])/(this->X[i+1]-this->X[i])
        << '\n';
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer // order
  LinearSpline::coeffs( real_type cfs[], real_type nodes[], bool transpose ) const {
    integer n = this->npts > 0 ? this->npts-1 : 0;
    for ( integer i = 0; i < n; ++i ) {
      nodes[i] = this->X[i];
      real_type a = this->Y[i];
      real_type b = (this->Y[i+1]-this->Y[i])/(this->X[i+1]-this->X[i]);
      if ( transpose ) {
        cfs[2*i+1] = a;
        cfs[2*i+0] = b;
      } else {
        cfs[n+i] = a;
        cfs[i]   = b;
      }
    }
    return 2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  LinearSpline::order( ) const { return 2; }

}
