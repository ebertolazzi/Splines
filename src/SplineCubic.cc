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

/*\
 |   ####  #    # #####  #  ####
 |  #    # #    # #    # # #    #
 |  #      #    # #####  # #
 |  #      #    # #    # # #
 |  #    # #    # #    # # #    #
 |   ####   ####  #####  #  ####
\*/

namespace Splines {

  using namespace std; // load standard namspace

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  CubicSpline_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts,
    real_type       ddy0,
    real_type       ddyn
  ) {

    size_t n = size_t(npts > 0 ? npts-1 : 0);

    vector<real_type> L(n+1);
    vector<real_type> D(n+1);
    vector<real_type> U(n+1);
    vector<real_type> Z(n+1);

    size_t i;
    for ( i = 1; i < n; ++i ) {
      L[i] = (X[i]-X[i-1])/(X[i+1]-X[i-1]);
      U[i] = (X[i+1]-X[i])/(X[i+1]-X[i-1]);
      D[i] = 2;
      Z[i] = 6 * ( (Y[i+1]-Y[i])/(X[i+1]-X[i]) -
                   (Y[i]-Y[i-1])/(X[i]-X[i-1]) ) / ( X[i+1] - X[i-1] );
    }

    L[0] = 0; D[0] = 1; U[0] = 0; Z[0] = ddy0;
    L[n] = 0; D[n] = 1; U[n] = 0; Z[n] = ddyn;

    if ( n > 2 ) {
      i = 0;
      do {
        Z[i]   /= D[i];
        U[i]   /= D[i];
        D[i+1] -= L[i+1] * U[i];
        Z[i+1] -= L[i+1] * Z[i];
      } while ( ++i < n );

      Z[i] /= D[i];

      do {
        --i;
        Z[i] -= U[i] * Z[i+1];
      } while ( i > 0 );
    }

    for ( i = 0; i < n; ++i ) {
      real_type DX = X[i+1] - X[i];
      Yp[i] = (Y[i+1]-Y[i])/DX - (Z[i]/3 + Z[i+1]/6) * DX;
    }
    real_type DX = X[n] - X[n-1];
    Yp[n] = Yp[n-1] + DX * 0.5*(Z[n-1] + Z[n]);

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSpline::build() {
    SPLINE_ASSERT(
      this->npts > 1,
      "CubicSpline::build(): npts = " << this->npts <<
      " not enought points"
    )
    integer ibegin = 0;
    integer iend   = 0;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < this->npts && this->X[iend-1] < this->X[iend] ) {}
      real_type d0 = ibegin == 0          ? ddy0 : 0;
      real_type d1 = iend   == this->npts ? ddyn : 0;
      CubicSpline_build(
        this->X+ibegin,
        this->Y+ibegin,
        this->Yp+ibegin,
        iend-ibegin, d0, d1
      );
      ibegin = iend;
    } while ( iend < this->npts );
    
    SPLINE_CHECK_NAN(this->Yp,"CubicSpline::build(): Yp",this->npts);
  }

  #ifndef SPLINES_DO_NOT_USE_GENERIC_CONTAINER

  using GenericContainerNamespace::GC_VEC_REAL;
  using GenericContainerNamespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["x"]
    // gc["y"]
    //
    */
    SPLINE_ASSERT(
      gc.exists("x"),
      "CubicSpline[" << this->_name << "]::setup missing `x` field!"
    )
    SPLINE_ASSERT(
      gc.exists("y"),
      "CubicSpline[" << this->_name << "]::setup missing `y` field!"
    )

    GenericContainer const & gc_x = gc("x");
    GenericContainer const & gc_y = gc("y");

    vec_real_type x, y;
    {
      std::ostringstream ost;
      ost << "CubicSpline[" << this->_name << "]::setup, field `x'";
      gc_x.copyto_vec_real ( x, ost.str().c_str() );
    }
    {
      std::ostringstream ost;
      ost << "CubicSpline[" << this->_name << "]::setup, field `y'";
      gc_y.copyto_vec_real ( y, ost.str().c_str() );
    }
    ddy0 = ddyn = 0;
    gc.get_if_exists( "ddy_begin", ddy0 );
    gc.get_if_exists( "ddy_end",   ddyn );

    this->build( x, y );
  }
  #endif
}
