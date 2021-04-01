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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

/*\
 |   ####  #    # #####  #  ####
 |  #    # #    # #    # # #    #
 |  #      #    # #####  # #
 |  #      #    # #    # # #
 |  #    # #    # #    # # #    #
 |   ####   ####  #####  #  ####
\*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namspace
#endif

namespace Splines {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  inline
  real_type
  deriv2_3p_L(
    real_type SR,  real_type hR,
    real_type SRR, real_type hRR
  ) {
    return 2*(SRR-SR)/(hR+hRR);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  real_type
  deriv2_4p_L(
    real_type SR,   real_type hR,
    real_type SRR,  real_type hRR,
    real_type SRRR, real_type hRRR
  ) {
    real_type t5  = hRR*hRR;
    real_type t12 = SR-SRR;
    real_type t13 = t12*hRRR;
    real_type t19 = hR*hR;
    real_type t23 = hRRR*hRRR;
    real_type bot = (hRR+hRRR)*(hR+hRR+hRRR)*(hR+hRR);
    real_type top = t5*(6*SRR-4*SR-2*SRRR) + 6*hRR*((2*SRR-SR-SRRR)*hR-t13) +
                    4*t19*(SRR-SRRR) - 6*hR*t13 - 2*t12*t23;
    return top/bot;
  }

  static
  real_type
  deriv2_5p_L(
    real_type SR,    real_type hR,
    real_type SRR,   real_type hRR,
    real_type SRRR,  real_type hRRR,
    real_type SRRRR, real_type hRRRR
  ) {
    real_type t1  = SRRRR-SRRR;
    real_type t3  = hRR*hRR;
    real_type t4  = t3*t3;
    real_type t10 = 6*SRRRR;
    real_type t16 = SR-2*SRR+SRRR;
    real_type t17 = t16*hRRRR;
    real_type t25 = hRRR*hRRR;
    real_type t29 = 60*SRRR;
    real_type t30 = 30*SRRRR;
    real_type t35 = SR-(8.0/5.0)*SRR+(3.0/5.0)*SRRR;
    real_type t41 = hR*hR;
    real_type t49 = hRRRR*hRRRR;
    real_type t55 = 24*SR;
    real_type t56 = 30*SRR;
    real_type t60 = t25*hRRR;
    real_type t75 = 12*SR;
    real_type t76 = 60*SRR;
    real_type t91 = t41*hR;
    real_type t102 = t49*hRRRR;
    real_type t105 = SR-(3.0/2.0)*SRR+SRRR/2;
    real_type t110 = SRR-SR;
    real_type t112 = t25*t25;
    real_type t118 = -t110;
    real_type t119 = t118*hRRRR;
    real_type t134 = t118*t49;
    real_type t142 = SR-(10.0/3.0)*SRR+(7.0/3.0)*SRRR;
    real_type t149 = t118*t102;
    real_type t156 = SRR-SRRR;
    real_type t159 = (5.0/2.0)*t41*t156*hRRRR;
    real_type MapleGenVar2 = 2*t4*hRR*t1 + t4*( hRRR*(12*SRR-6*SR-12*SRRR+t10) + 12*hR*t1 - 6*t17 );
    real_type MapleGenVar3 = MapleGenVar2 +
                             t3*hRR*(
                               t25*(32*SRR-20*SR-18*SRRR+t10) +
                               hRRR*(hR*(t30-18*SR+48*SRR-t29)-30*hRRRR*t35) +
                               24*t41*t1 - 10*t49*t35 + hRRRR*(48*SRR-18*SR-30*SRRR)*hR
                             );
    real_type MapleGenVar1 = MapleGenVar3 +
                             t3*(
                               t60*(2*SRRRR-t55+t56-8*SRRR) +
                               t25*( (-48*SR+96*SRR-72*SRRR+24*SRRRR)*hR - hRRRR*(48*SR-60*SRR+12*SRRR) ) +
                               hRRR*(t41*(48*SRRRR-96*SRRR-t75+t76)- 72*hR*t17-t49*(28*SR-36*SRR+8*SRRR)) +
                               20*t91*t1 - 12*t41*(SR-5*SRR+4*SRRR)*hRRRR - 24*hR*t16*t49 - 4*t105*t102
                             ) +
                             hRR*( 12*t112*t110 + t60*(hR*(t76-42*SR-24*SRRR+t10)-30*t119) +
                                   t25*(t41*(80*SRR-84*SRRR+28*SRRRR-t55) - (84*SR-120*SRR+36*SRRR)*hRRRR*hR - 24*t134 ) +
                                   hRRR*( t91*(t56-t29+t30) - 36*t41*hRRRR*t142 - 48*hR*t105*t49 - 6*t149 )
                                   -( 6*(t16*t102 - t91*t1) + 12*(hR*t142*t49 - t159) )*hR
                             );
    real_type top = ((28*SRR-12*SR-16*SRRR)*t49*hRRR + 6*t41*t156*hRRRR)*t41 +
                    (((6*SRR-12*SRRR+6*SRRRR)*t91+12*t159-6*t149)*hRRR
                      -12*t112*t118-24*t25*t134-30*t60*t119+10*t41*t156*t49)*hR +
                    2*t112*hRRR*t110 +
                    ((4*SRRRR-t75-16*SRRR+24*SRR)*t60-24*t25*t17+4*t156*t102)*t41 +
                    (t91*(20*SRR-30*SRRR+10*SRRRR)-2*t149)*t25 -
                    6*t112*t118*hRRRR+MapleGenVar1-6*t60*t134;
    real_type bot = (hR+hRR)*(hRRR+hRRRR)*(hRR+hRRR)*(hRR+hRRR+hRRRR)*(hR+hRR+hRRR)*(hR+hRR+hRRR+hRRRR);
    return top/bot;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  inline
  real_type
  deriv2_3p_R(
    real_type SL,   real_type hL,
    real_type SLL,  real_type hLL
  ) {
    return 2*(SL-SLL)/(hL+hLL);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  real_type
  deriv2_4p_R(
    real_type SL,   real_type hL,
    real_type SLL,  real_type hLL,
    real_type SLLL, real_type hLLL
  ) {
    real_type t5  = hLL*hLL;
    real_type t12 = SL-SLL;
    real_type t13 = t12*hLLL;
    real_type t19 = hL*hL;
    real_type t23 = hLLL*hLLL;
    real_type bot = (hLL+hLLL)*(hL+hLL+hLLL)*(hL+hLL);
    real_type top = t5*(4*SL-6*SLL+2*SLLL) + 6*hLL*((SL+SLLL-2*SLL)*hL+t13) +
                    4*t19*(SLLL-SLL) + 6*hL*t13+2*t12*t23;
    return top/bot;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  real_type
  deriv2_5p_R(
    real_type SL,    real_type hL,
    real_type SLL,   real_type hLL,
    real_type SLLL,  real_type hLLL,
    real_type SLLLL, real_type hLLLL
  ) {
    real_type t1  = SLLL-SLLLL;
    real_type t3  = hLL*hLL;
    real_type t4  = t3*t3;
    real_type t10 = 6*SLLLL;
    real_type t16 = SL-2*SLL+SLLL;
    real_type t17 = t16*hLLLL;
    real_type t25 = hLLL*hLLL;
    real_type t29 = 60*SLLL;
    real_type t30 = 30*SLLLL;
    real_type t35 = SL-(8.0/5.0)*SLL+(3.0/5.0)*SLLL;
    real_type t41 = hL*hL;
    real_type t49 = hLLLL*hLLLL;
    real_type t55 = 24*SL;
    real_type t56 = 30*SLL;
    real_type t60 = t25*hLLL;
    real_type t75 = 12*SL;
    real_type t76 = 60*SLL;
    real_type t91 = t41*hL;
    real_type t102 = t49*hLLLL;
    real_type t105 = SL-(3.0/2.0)*SLL+SLLL/2;
    real_type t110 = SL-SLL;
    real_type t112 = t25*t25;
    real_type t118 = t110*hLLLL;
    real_type t133 = t110*t49;
    real_type t141 = SL-(10.0/3.0)*SLL+(7.0/3.0)*SLLL;
    real_type t148 = t110*t102;
    real_type t156 = t41*(SLL-SLLL)*hLLLL;
    real_type t157 = (5.0/2.0)*t156;
    real_type MapleGenVar2 = 2*t4*hLL*t1 + t4*(hLLL*(6*SL+12*(SLLL-SLL)-t10)+ 12*hL*t1+6*t17);
    real_type MapleGenVar3 = MapleGenVar2 +
                             t3*hLL*(t25*(20*SL-32*SLL+18*SLLL-t10)+
                                     hLLL*(hL*(18*SL-48*SLL+t29-t30)+30*hLLLL*t35)+
                                     24*t41*t1+10*t49*t35+
                                     (18*SL-48*SLL+30*SLLL)*hLLLL*hL);
    real_type MapleGenVar1 = MapleGenVar3 +
                             t3*( t60*(t55-t56+8*SLLL-2*SLLLL) +
                                  t25*((48*SL-96*SLL+72*SLLL-24*SLLLL)*hL +
                                        12*(4*SL-5*SLL+SLLL)*hLLLL) +
                                  hLLL*(t41*(t75-t76+96*SLLL-48*SLLLL) +
                                        72*hL*t17+(28*SL-36*SLL+8*SLLL)*t49) +
                                  20*t91*t1 + 24*hL*t16*t49 + 4*t105*t102 +
                                  12*t41*(SL-5*SLL+4*SLLL)*hLLLL
                                ) +
                             hLL*(
                               12*t112*t110 +
                               t60*(hL*(42*SL-t76+24*SLLL-t10)+30*t118) +
                               t25*(
                                 t41*(t55-80*SLL+84*SLLL-28*SLLLL) +
                                 12*hLLLL*(7*SL-10*SLL+3*SLLL)*hL+24*t133
                               ) +
                               hLLL*(t91*(t29-t56-t30) + 36*t41*t141*hLLLL + 48*hL*t105*t49 + 6*t148) +
                               6*(t91*t1  + 2*(hL*t141*t49-t157) + t16*t102)*hL
                            );
    real_type top = MapleGenVar1 + 2*t112*hLLL*t110 +
                    6*t112*t110*(2*hL+hLLLL) - t156*(hL+hLLLL)*(6*hL+4*hLLLL) +
                    t60*( t41*(t75-24*SLL+16*SLLL-4*SLLLL) + 30*hL*t118+6*t133 ) +
                    t25*( t91*(30*SLLL-20*SLL-10*SLLLL)+24*(t41*t17+hL*t133)+2*t148 ) +
                    hLLL*hL*( 6*t91*(2*SLLL-SLL-SLLLL)-12*t157 + hL*(12*SL-28*SLL+16*SLLL)*t49+6*t148);
    real_type bot = (hL+hLL)*(hLLL+hLLLL)*(hLL+hLLL)*(hLL+hLLL+hLLLL)*(hL+hLL+hLLL)*(hL+hLL+hLLL+hLLLL);
    return top/bot;
  }

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
    real_type const * X,
    real_type const * Y,
    real_type       * Yp,
    real_type       * Ypp,
    real_type       * L,
    real_type       * D,
    real_type       * U,
    integer           npts,
    CUBIC_SPLINE_TYPE_BC bc0,
    CUBIC_SPLINE_TYPE_BC bcn
  ) {

    size_t n = size_t(npts > 0 ? npts-1 : 0);
    real_type * Z = Ypp;

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
    case EXTRAPOLATE_BC:
      L[0] = 0; D[0] = 1; U[0] = 0;
      if ( npts == 2 ) {
        Z[0] = 0;
      } else if ( npts == 3 ) {
        real_type hR  = X[1] - X[0];
        real_type hRR = X[2] - X[1];
        real_type SR  = (Y[1] - Y[0])/hR;
        real_type SRR = (Y[2] - Y[1])/hRR;
        Z[0] = deriv2_3p_L( SR, hR, SRR, hRR );
      } else if ( npts == 4 ) {
        real_type hR   = X[1] - X[0];
        real_type hRR  = X[2] - X[1];
        real_type hRRR = X[3] - X[2];
        real_type SR   = (Y[1] - Y[0])/hR;
        real_type SRR  = (Y[2] - Y[1])/hRR;
        real_type SRRR = (Y[3] - Y[2])/hRRR;
        Z[0] = deriv2_4p_L( SR, hR, SRR, hRR, SRRR, hRRR );
      } else {
        real_type hR    = X[1] - X[0];
        real_type hRR   = X[2] - X[1];
        real_type hRRR  = X[3] - X[2];
        real_type hRRRR = X[4] - X[3];
        real_type SR    = (Y[1] - Y[0])/hR;
        real_type SRR   = (Y[2] - Y[1])/hRR;
        real_type SRRR  = (Y[3] - Y[2])/hRRR;
        real_type SRRRR = (Y[4] - Y[3])/hRRRR;
        Z[0] = deriv2_5p_L( SR, hR, SRR, hRR, SRRR, hRRR, SRRRR, hRRRR );
      }
      break;
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
    case EXTRAPOLATE_BC:
      L[n] = 0;  D[n] = 1; U[n] = 0;
      if ( npts == 2 ) {
        Z[n] = 0;
      } else if ( npts == 3 ) {
        real_type hL  = X[n] - X[n-1];
        real_type hLL = X[n-1] - X[n-2];
        real_type SL  = (Y[n] - Y[n-1])/hL;
        real_type SLL = (Y[n-1] - Y[n-2])/hLL;
        Z[n] = deriv2_3p_R( SL, hL, SLL, hLL );
      } else if ( npts == 4 ) {
        real_type hL   = X[n] - X[n-1];
        real_type hLL  = X[n-1] - X[n-2];
        real_type hLLL = X[n-2] - X[n-3];
        real_type SL   = (Y[n] - Y[n-1])/hL;
        real_type SLL  = (Y[n-1] - Y[n-2])/hLL;
        real_type SLLL = (Y[n-2] - Y[n-3])/hLLL;
        Z[n] = deriv2_4p_R(  SL, hL, SLL, hLL, SLLL, hLLL );
      } else {
        real_type hL    = X[n] - X[n-1];
        real_type hLL   = X[n-1] - X[n-2];
        real_type hLLL  = X[n-2] - X[n-3];
        real_type hLLLL = X[n-3] - X[n-4];
        real_type SL    = (Y[n] - Y[n-1])/hL;
        real_type SLL   = (Y[n-1] - Y[n-2])/hLL;
        real_type SLLL  = (Y[n-2] - Y[n-3])/hLLL;
        real_type SLLLL = (Y[n-3] - Y[n-4])/hLLLL;
        Z[n] = deriv2_5p_R(  SL, hL, SLL, hLL, SLLL, hLLL, SLLLL, hLLLL );
      }
      break;
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
  CubicSpline_build(
    real_type const * X,
    real_type const * Y,
    real_type       * Yp,
    integer           npts,
    CUBIC_SPLINE_TYPE_BC bc0,
    CUBIC_SPLINE_TYPE_BC bcn
  ) {
    Utils::Malloc<real_type> mem("CubicSpline_build");
    mem.allocate( size_t(4*npts) );
    real_type * L = mem( size_t( npts ) );
    real_type * D = mem( size_t( npts ) );
    real_type * U = mem( size_t( npts ) );
    real_type * Z = mem( size_t( npts ) );
    CubicSpline_build( X, Y, Yp, Z, L, D, U, npts, bc0, bcn );
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CubicSpline::build() {
    string msg = fmt::format("CubicSpline[{}]::build():", m_name );
    UTILS_ASSERT(
      m_npts > 1,
      "{} npts = {} not enought points\n",
      msg, m_npts
    );
    Utils::checkNaN( m_X, (msg+" X").c_str(), m_npts, __LINE__, __FILE__ );
    Utils::checkNaN( m_Y, (msg+" Y").c_str(), m_npts, __LINE__, __FILE__ );
    integer ibegin = 0;
    integer iend   = 0;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < m_npts && m_X[iend-1] < m_X[iend] ) {}
      CUBIC_SPLINE_TYPE_BC seg_bc0 = NOT_A_KNOT;
      CUBIC_SPLINE_TYPE_BC seg_bcn = NOT_A_KNOT;
      if ( ibegin == 0      ) seg_bc0 = m_bc0;
      if ( iend   == m_npts ) seg_bcn = m_bcn;
      CubicSpline_build(
        m_X+ibegin,
        m_Y+ibegin,
        m_Yp+ibegin,
        iend - ibegin,
        seg_bc0, seg_bcn
      );
      ibegin = iend;
    } while ( iend < m_npts );

    Utils::checkNaN( m_Yp, (msg+" Yp").c_str(), m_npts, __LINE__, __FILE__ );
  }

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  using GenericContainerNamespace::GC_VEC_REAL;
  using GenericContainerNamespace::vec_real_type;
  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*!
   *
   * Setup a spline using a `GenericContainer`
   *
   * - gc("xdata") vector with the `x` coordinate of the data
   * - gc("ydata") vector with the `y` coordinate of the data
   *
   * may contain
   * - gc("bc_begin") and/or gc("bc_end") 
   *   - "extrapolate" extrapolate the boundary condition
   *   - "natural"     make second derivative 0 at the border
   *   - "parabolic"   make third derivative 0 at the border
   *   - "not_a_knot"  not a knot condition of De Boor 
   */
  void
  CubicSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string msg = fmt::format("CubicSpline[{}]::setup( gc ):", m_name );
    UTILS_ASSERT( gc.exists("xdata"), "{} missing `xdata` field!\n", msg );
    UTILS_ASSERT( gc.exists("ydata"), "{} missing `y`data field!\n", msg );

    GenericContainer const & gc_x = gc("xdata");
    GenericContainer const & gc_y = gc("ydata");

    vec_real_type x, y;
    {
      std::string ff = fmt::format( "{}, field `xdata'", msg );
      gc_x.copyto_vec_real ( x, ff.c_str() );
    }
    {
      std::string ff = fmt::format( "{}, field `ydata'", msg );
      gc_y.copyto_vec_real ( y, ff.c_str() );
    }
    if ( gc.exists("bc_begin") ) {
      std::string const & bc = gc("bc_begin").get_string();
      if      ( bc == "extrapolate" ) m_bc0 = EXTRAPOLATE_BC;
      else if ( bc == "natural"     ) m_bc0 = NATURAL_BC;
      else if ( bc == "parabolic"   ) m_bc0 = PARABOLIC_RUNOUT_BC;
      else if ( bc == "not_a_knot"  ) m_bc0 = NOT_A_KNOT;
      else {
        UTILS_ERROR( "{} unknow initial bc: {}\n", msg, bc );
      }
    } else {
      UTILS_WARNING( false,
        "{}, missing field `bc_begin` using `extrapolate` as default value\n", msg
      );
    }

    if ( gc.exists("bc_end") ) {
      std::string const & bc = gc("bc_end").get_string();
      if      ( bc == "extrapolate" ) m_bcn = EXTRAPOLATE_BC;
      else if ( bc == "natural"     ) m_bcn = NATURAL_BC;
      else if ( bc == "parabolic"   ) m_bcn = PARABOLIC_RUNOUT_BC;
      else if ( bc == "not_a_knot"  ) m_bcn = NOT_A_KNOT;
      else {
        UTILS_ERROR( "{} unknow final bc: {}\n", msg, bc );
      }
    } else {
      UTILS_WARNING( false,
        "{}, missing field `bc_begin` using `extrapolate` as default value\n", msg
      );
    }
    this->build( x, y );
  }

}
