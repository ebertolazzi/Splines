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

  //! Use externally allocated memory for `npts` points
  void
  LinearSpline::reserve_external(
    integer      n,
    real_type *& p_x,
    real_type *& p_y
  ) {
    if ( !_external_alloc ) baseValue.free();
    npts            = 0;
    npts_reserved   = n;
    _external_alloc = true;
    X = p_x;
    Y = p_y;
  }

  void
  LinearSpline::reserve( integer n ) {
    if ( _external_alloc && n <= npts_reserved ) {
      // nothing to do!, already allocated
    } else {
      baseValue.allocate( size_t(2*n) );
      npts_reserved   = n;
      _external_alloc = false;
      X = baseValue( size_t(n) );
      Y = baseValue( size_t(n) );
    }
    npts = lastInterval = 0;
  }
  
  void
  LinearSpline::clear(void) {
    if ( !_external_alloc ) baseValue.free();
    npts = npts_reserved = 0;
    _external_alloc = false;
    X = Y = nullptr;
  }

  void
  LinearSpline::writeToStream( ostream_type & s ) const {
    integer nseg = npts > 0 ? npts - 1 : 0;
    for ( integer i = 0; i < nseg; ++i )
      s << "segment N." << setw(4) << i
        << " X:[ " << X[i] << ", " << X[i+1] << " ] Y:[ " << Y[i] << ", " << Y[i+1] 
        << " ] slope: " << (Y[i+1]-Y[i])/(X[i+1]-X[i])
        << '\n';
  }

  integer // order
  LinearSpline::coeffs( real_type cfs[], real_type nodes[], bool transpose ) const {
    integer n = npts > 0 ? npts-1 : 0;
    for ( integer i = 0; i < n; ++i ) {
      nodes[i] = X[i];
      real_type a = Y[i];
      real_type b = (Y[i+1]-Y[i])/(X[i+1]-X[i]);
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
  
  integer
  LinearSpline::order( ) const { return 2; }

}
