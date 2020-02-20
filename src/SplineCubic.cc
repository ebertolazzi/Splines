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

  /*\
     Sistema lineare da risolvere

     D U UU
     L D U
       L D U
         L D U
           .....
              L D U
                L D U
                LL L D

  \*/

  void
  CubicSpline_build(
    real_type const      X[],
    real_type const      Y[],
    real_type            Yp[],
    integer              npts,
    CUBIC_SPLINE_TYPE_BC bc0,
    CUBIC_SPLINE_TYPE_BC bcn
  ) {

    size_t n = size_t(npts > 0 ? npts-1 : 0);

    vector<real_type> buffer(4*(n+1));
    real_type * ptr = &buffer.front();
    real_type * L = ptr; ptr += npts;
    real_type * D = ptr; ptr += npts;
    real_type * U = ptr; ptr += npts;
    real_type * Z = ptr;

    size_t i;
    for ( i = 1; i < n; ++i ) {
      real_type HL = X[i] - X[i-1];
      real_type HR = X[i+1] - X[i];
      real_type HH = HL+HR;
      L[i] = HL/HH;
      U[i] = HR/HH;
      D[i] = 2;
      Z[i] = 6 * ( (Y[i+1]-Y[i])/HR - (Y[i]-Y[i-1])/HL ) / HH;
    }

    real_type UU = 0, LL = 0;

    switch ( bc0 ) {
    case NATURAL_BC:
      L[0] = 0; D[0] = 1; U[0] = 0; Z[0] = 0;
      break;
    case PARABOLIC_RUNOUT_BC:
      L[0] = 0; D[0] = 1; U[0] = -1; Z[0] = 0;
      break;
    case NOT_A_KNOT:
      {
        real_type r = (X[1] - X[0])/(X[2] - X[1]);
        // v0 - v1*(1+r) + r*v2 == 0
        L[0] = 0;
        D[0] = 1;
        U[0] = -(1+r);
        UU   = r;
        Z[0] = 0;
      }
      break;
    }

    switch ( bcn ) {
    case NATURAL_BC:
      L[n] = 0;  D[n] = 1; U[n] = 0; Z[n] = 0;
      break;
    case PARABOLIC_RUNOUT_BC:
      L[n] = -1; D[n] = 1; U[n] = 0; Z[n] = 0;
      break;
    case NOT_A_KNOT:
      {
        real_type r = (X[n-1] - X[n-2])/(X[n] - X[n-1]);
        // r*v0 - v1*(1+r) + v2 == 0
        U[n] = 0;
        D[n] = 1;
        L[n] = -(1+r);
        LL   = r;
        Z[n] = 0;
      }
      break;
    }

    if ( n > 2 ) {
      Z[0] /= D[0];
      U[0] /= D[0];
      UU   /= D[0];
      D[1] -= L[1] * U[0];
      U[1] -= L[1] * UU;
      Z[1] -= L[1] * Z[0];
      i = 1;
      do {
        Z[i]   /= D[i];
        U[i]   /= D[i];
        D[i+1] -= L[i+1] * U[i];
        Z[i+1] -= L[i+1] * Z[i];
      } while ( ++i < n );

      D[i] -= LL * U[i-2];
      Z[i] -= LL * Z[i-2];

      Z[i] /= D[i];

      do {
        --i;
        Z[i] -= U[i] * Z[i+1];
      } while ( i > 0 );

      Z[0] -= UU * Z[2];
    }

    for ( i = 0; i < n; ++i ) {
      real_type DX = X[i+1] - X[i];
      Yp[i] = (Y[i+1]-Y[i])/DX - (2*Z[i] + Z[i+1]) * (DX/6);
    }
    real_type DX2 = (X[n] - X[n-1])/2;
    Yp[n] = Yp[n-1] + DX2 * (Z[n-1] + Z[n]);

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSpline::build() {
    SPLINE_ASSERT(
      this->npts > 1,
      "CubicSpline::build(): npts = " << this->npts << " not enought points"
    )
    integer ibegin = 0;
    integer iend   = 0;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < this->npts && this->X[iend-1] < this->X[iend] ) {}
      CUBIC_SPLINE_TYPE_BC seg_bc0 = NOT_A_KNOT;
      CUBIC_SPLINE_TYPE_BC seg_bcn = NOT_A_KNOT;
      if ( ibegin == 0         ) seg_bc0 = this->bc0;
      if ( iend  == this->npts ) seg_bcn = this->bcn;
      CubicSpline_build(
        this->X+ibegin,
        this->Y+ibegin,
        this->Yp+ibegin,
        iend - ibegin,
        seg_bc0, seg_bcn
      );
      ibegin = iend;
    } while ( iend < this->npts );
    
    SPLINE_CHECK_NAN( this->Yp, "CubicSpline::build(): Yp", this->npts );
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
    if ( gc.exists("bc_begin") ) {
      std::string const & bc = gc("bc_begin").get_string();
      if      ( bc == "natural"    ) this->bc0 = NATURAL_BC;
      else if ( bc == "parabolic"  ) this->bc0 = PARABOLIC_RUNOUT_BC;
      else if ( bc == "not_a_knot" ) this->bc0 = NOT_A_KNOT;
      else {
        SPLINE_DO_ERROR(
          "CubicSpline[" << this->_name << "]::setup unknow initial bc:" << bc
        )
      }
    }

    if ( gc.exists("bc_end") ) {
      std::string const & bc = gc("bc_end").get_string();
      if      ( bc == "natural"    ) this->bcn = NATURAL_BC;
      else if ( bc == "parabolic"  ) this->bcn = PARABOLIC_RUNOUT_BC;
      else if ( bc == "not_a_knot" ) this->bcn = NOT_A_KNOT;
      else {
        SPLINE_DO_ERROR(
          "CubicSpline[" << this->_name << "]::setup unknow final bc:" << bc
        )
      }
    }
    this->build( x, y );
  }
  #endif
}
