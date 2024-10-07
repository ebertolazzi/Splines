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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

#include "SplinesUtils.hh"
#include "Utils_fmt.hh"

namespace Splines {

  using std::abs;

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  /*
  //    __ _       _ _             _ _  __  __
  //   / _(_)_ __ (_) |_ ___    __| (_)/ _|/ _| ___ _ __ ___ _ __   ___ ___
  //  | |_| | '_ \| | __/ _ \  / _` | | |_| |_ / _ \ '__/ _ \ '_ \ / __/ _ \
  //  |  _| | | | | | ||  __/ | (_| | |  _|  _|  __/ | |  __/ | | | (_|  __/
  //  |_| |_|_| |_|_|\__\___|  \__,_|_|_| |_|  \___|_|  \___|_| |_|\___\___|
  */

  real_type
  first_deriv4p_R(
    real_type DLLL, real_type hLLL,
    real_type DLL,  real_type hLL,
    real_type DL,   real_type hL
  ) {
    real_type t2  = hL*hL;
    real_type t19 = hLL*hLL;
    real_type t22 = DL-DLL/2;
    real_type t26 = hLLL*hLLL;
    real_type t32 = hLL+hLLL;
    real_type t33 = t32*t32;
    real_type top = ((DLLL-DLL)*t2+4*t19*DL-3*t19*DLL+t19*DLLL+6*hLL*t22*hLLL+2*t22*t26)*hL +
                    ((3*DL-4*DLL+2*DLLL)*hLL+3*DL*hLLL-2*DLL*hLLL)*t2+t33*hLL*DL;
    real_type bot = t32*(hL+hLL+hLLL)*(hL+hLL);
    return top/bot;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  first_deriv4p_L(
    real_type DR,   real_type hR,
    real_type DRR,  real_type hRR,
    real_type DRRR, real_type hRRR
  ) {
    real_type t2  = hR*hR;
    real_type t19 = hRR*hRR;
    real_type t22 = DR-DRR/2;
    real_type t26 = hRRR*hRRR;
    real_type t32 = hRR+hRRR;
    real_type t33 = t32*t32;
    real_type top = ( (DRRR-DRR)*t2 + t19*(4*DR-3*DRR+DRRR) + 6*hRR*t22*hRRR+2*t22*t26)*hR +
                    ((3*DR-4*DRR+2*DRRR)*hRR+3*DR*hRRR-2*DRR*hRRR)*t2+t33*hRR*DR;
    real_type bot = t32*(hR+hRR+hRRR)*(hR+hRR);
    return top/bot;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  first_deriv5p_C(
    real_type DLL, real_type hLL,
    real_type DL,  real_type hL,
    real_type DR,  real_type hR,
    real_type DRR, real_type hRR
  ) {
    real_type hhLL = hL+hLL;
    real_type hhRR = hR+hRR;
    real_type DDRR = (DR*hR+DRR*hRR)/hhRR;
    real_type DDLL = (DL*hL+DLL*hLL)/hhLL;
    real_type t3  = DL*hR+DR*hL;
    real_type t4  = hhRR*hhRR;
    real_type t6  = hR*hR;
    real_type t8  = hL*hL;
    real_type t10 = -t6*DL+t8*DR;
    real_type t12 = hR*hL;
    real_type t13 = hR+hL;
    real_type t14 = t13*DDRR;
    real_type t17 = hhLL*hhLL;
    real_type t20 = t4*hhRR;
    real_type t22 = t6*hR;
    real_type t24 = t8*hL;
    real_type t26 = -t22*DL-t24*DR;
    real_type t41 = -hhRR+hR;
    real_type t43 = hhRR+hL;
    real_type top = ( (hhRR*t10-t12*t14+t3*t4)*hhLL + hhRR*t26 + (hR*t24-hL*t22)*DDRR + t20*t3)*t17 +
                    hhLL*(t14*t6*t8-t10*t20+t26*t4) + t13*t43*t41*hhRR*DDLL*t12;
    real_type bot = t41*(hhRR+hhLL)*t43*(hR+hhLL)*t13*(hL-hhLL);
    return top/bot;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  first_deriv5p_L(
    real_type DLL, real_type hLL,
    real_type DL,  real_type hL,
    real_type DR,  real_type hR,
    real_type DRR, real_type hRR
  ) {
    real_type hhLL = hL+hLL;
    real_type hhRR = hR+hRR;
    real_type DDRR = (DR*hR+DRR*hRR)/hhRR;
    real_type DDLL = (DL*hL+DLL*hLL)/hhLL;

    real_type t1   = -DR+DDRR;
    real_type t3   = -DR+DDLL;
    real_type t5   = DDLL-DDRR;
    real_type t8   = hL*hL;
    real_type t9   = t8*t8;
    real_type t12  = hhLL*hhLL;
    real_type t14  = hhRR*hhRR;
    real_type t16  = hR*hR;
    real_type t21  = t12*hhLL;
    real_type t25  = 3*DR-4*DL;
    real_type t28  = DL-(3.0/4.0)*DDRR;
    real_type t38  = t14*hhRR;
    real_type t41  = DL-(3.0/4.0)*DDLL;
    real_type t48  = t16*hR;
    real_type t58  = DL-(2.0/3.0)*DDRR;
    real_type t70  = -hhRR+hR;
    real_type t80  = 2*DL-DR;
    real_type t83  = DL-DDRR/2;
    real_type t103 = hhRR+hhLL;
    real_type t105 = hR+hhLL;
    real_type top, bot;
    {
      real_type t51 = (3*DL-2*DDLL)*t70;
      real_type t52 = hhRR*hhRR;
      real_type t55 = hR*hR;
      top = (((4*t16*t41+t12*t25)*hhRR+4*(t12*t28-t14*t41)*hR+
               hhLL*(t14*t25+4*t16*t28)+t21*t1+t38*t3-t5*t48)*t8+
               hhRR*t3*t9-hR*t5*t9+(2*DL-DDLL)*t16*t14*t70+
               hhLL*t1*t9+t21*t14*t80-2*t21*t16*t83+t12*(t38*t80-2*t48*t83))*hL+
               t8*(t52*hR*t51+hhRR*(t55*t21+t55*t51)
               -3*t21*hR*t58+(-t38*t55+3*t48*t58)*hhLL)
               +2*t9*(t14*t3-t5*t16-t12*t1)-t105*t70*t103*hhRR*DL*hR*hhLL;
    }
    bot = (hhRR+hL)*(hL-hhLL)*(hR+hL)*t70*t103*t105;
    return top/bot;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  first_deriv5p_R(
    real_type DLL, real_type hLL,
    real_type DL,  real_type hL,
    real_type DR,  real_type hR,
    real_type DRR, real_type hRR
  ) {
    real_type hhLL = hL+hLL;
    real_type hhRR = hR+hRR;
    real_type DDRR = (DR*hR+DRR*hRR)/hhRR;
    real_type DDLL = (DL*hL+DLL*hLL)/hhLL;
    real_type t1   = -DL+DDRR;
    real_type t3   = DDLL-DL;
    real_type t5   = DDLL-DDRR;
    real_type t8   = hR*hR;
    real_type t9   = t8*t8;
    real_type t12  = hhLL*hhLL;
    real_type t14  = hhRR*hhRR;
    real_type t16  = hL*hL;
    real_type t21  = t12*hhLL;
    real_type t25  = 3*DL-4*DR;
    real_type t28  = DR-(3.0/4.0)*DDRR;
    real_type t38  = t14*hhRR;
    real_type t41  = DR-(3.0/4.0)*DDLL;
    real_type t48  = t16*hL;
    real_type t55  = 2*DL-3*DR;
    real_type t58  = DR-(2.0/3.0)*DDRR;
    real_type t73  = hhRR+hL;
    real_type t80  = 2*DR-DL;
    real_type t83  = DR-DDRR/2;
    real_type t103 = hhRR+hhLL;
    real_type t105 = -hhLL+hL;
    real_type top, bot;
    {
      real_type t47 = 2*DDLL-3*DR;
      real_type t49 = hhRR*hhRR;
      real_type t53 = hL*hL;
      top = (((4*t16*t41+t12*t25)*hhRR+(-4*t12*t28+4*t14*t41)*hL+
            hhLL*(t14*t25+4*t16*t28)+t21*t1+t38*t3+t5*t48)*t8+hhRR*t3*t9+hL*t5*t9+(
            -2*DR+DDLL)*t16*t14*t73+hhLL*t1*t9+t21*t14*t80
            -2*t21*t16*t83+t12*(t38*t80+2*t48*t83))*hR+
            t8*(t49*hL*t73*t47+hhRR*(-t53*t73*t47+t55*t21)
            -3*t21*t58*hL+(-t38*t55+3*t48*t58)*hhLL)
            -t105*t73*t103*hhRR*DR*hL*hhLL
            +2*t9*(t5*t16+t12*t1-t14*t3);
    }
    bot = (hR-hhRR)*(hR+hhLL)*(hR+hL)*t103*t73*t105;
    return top/bot;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  deriv_left(
    real_type   SR,
    real_type   hR,
    real_type   SRR,
    real_type   hRR,
    real_type   DyRR,
    real_type & dy0,
    real_type & dy1
  ) {
    real_type t1 = hRR*hRR;
    real_type t2 = t1*SR;
    real_type t3 = DyRR/3;
    real_type t8 = hR*hR;
    real_type t14 = hR+hRR;
    real_type t24 = t14*t14;
    real_type hRRR = 3*hRR*hR;
    dy0 = (hRRR*(SR-SRR+t3)-t8*SRR+t8*DyRR+t2)/(hRR*t14);
    dy1 = ((hRRR+2*t8)*SRR-hRRR*t3-t8*DyRR+t2)/t24;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  deriv_right(
    real_type   SL,
    real_type   hL,
    real_type   SLL,
    real_type   hLL,
    real_type   DyLL,
    real_type & dyn,
    real_type & dynm1
  ) {
    real_type t1 = hLL*hLL;
    real_type t2 = t1*SL;
    real_type t3 = DyLL/3;
    real_type t8 = hL*hL;
    real_type t14 = hL+hLL;
    real_type t24 = t14*t14;
    dyn   = (3*hLL*(SL-SLL+t3)*hL-t8*SLL+t8*DyLL+t2)/(hLL*t14);
    dynm1 = ((3*hL*hLL+2*t8)*SLL-3*hL*hLL*t3-t8*DyLL+t2)/t24;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  second_deriv3p_C(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dpL,
    real_type dp0,
    real_type dpR
  ) {
    real_type hLR = hL*hR;
    real_type hL2 = hL*hL;
    real_type hL4 = hL2*hL2;
    real_type hR2 = hR*hR;
    real_type hR4 = hR2*hR2;
    real_type top = hL4*(6*SR-4*dp0-2*dpR) + hR4*(4*dp0+2*dpL-6*SL) +
                    hLR*( hL2*(10*SR-8*dp0-2*dpR) - hR2*(10*SL-8*dp0-2*dpL) );
    real_type tmp = hR+hL;
    real_type bot = tmp*tmp*tmp*hLR;
    return top/bot;
  }

  real_type
  second_deriv3p_C(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dp0
  ) {
    real_type hh    = hL+hR;
    real_type alpha = hL/hh;
    real_type beta  = hR/hh;
    real_type DDR   = (SR-dp0)/hR;
    real_type DDL   = (dp0-SL)/hL;
    return 2*(alpha*DDR+beta*DDL);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  second_deriv3p_L(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dpL,
    real_type dp0,
    real_type dpR
  ) {
    real_type hLR = hL * hR;
    real_type hL2 = hL * hL;
    real_type hL4 = hL2 * hL2;
    real_type hR2 = hR * hR;
    real_type hR4 = hR2 * hR2;
    real_type top = hL4 * (4 * SR - 2 * (dp0+dpR) )
                  + hR4 * (6 * SL - 2 * dp0 - 4 * dpL)
                  + hL2 * hR2 * (20*SL - 12*dp0 - 8*dpL)
                  + hLR * hR2 * (20*SL - 8*dp0 - 12*dpL)
                  + hL2 * hLR * (10*SR - 8*dp0 - 2*dpR);
    real_type hRpL = hR + hL;
    real_type bot = hR2 * hL * (hRpL*hRpL);
    return top/bot;
  }

  real_type
  second_deriv3p_L(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dpL
  ) {
    real_type hh    = hL+hR;
    real_type alpha = hL/hh;
    real_type beta  = hR/hh;
    real_type t1    = (6*SL-2*SR-4*dpL)/hL;
    real_type t3    = 2*(SL-dpL)/hL;
    real_type t2    = 3*t3;
    return t1*alpha*alpha+t2*alpha*beta+t3*beta*beta;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  second_deriv3p_R(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dpL,
    real_type dp0,
    real_type dpR
  ) {
    real_type hLR = hL * hR;
    real_type hL2 = hL * hL;
    real_type hL4 = hL2 * hL2;
    real_type hR2 = hR * hR;
    real_type hR4 = hR2 * hR2;
    real_type top = hL4 * (2 * dp0 + 4 * dpR - 6 * SR)
                  + hR4 * (2 * (dp0+dpL) - 4*SL)
                  - hL2 * hLR * (20*SR - 8*dp0 - 12*dpR)
                  - hL2 * hR2 * (20*SR - 12*dp0 - 8*dpR)
                  - hR2 * hLR * (10*SL - 8*dp0 - 2*dpL);
    real_type hRpL = hR + hL;
    real_type bot  = hR * hL2 * (hRpL*hRpL);

    return top/bot;
  }

  real_type
  second_deriv3p_R(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dpR
  ) {
    real_type hh    = hL+hR;
    real_type alpha = hL/hh;
    real_type beta  = hR/hh;
    real_type t1    = 2*(dpR-SR)/hR;
    real_type t2    = 3*t1;
    real_type t3    = (2*SL-6*SR+4*dpR)/hR;
    return t1*alpha*alpha+t2*alpha*beta+t3*beta*beta;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  first_derivative_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts
  ) {

    size_t n{ npts > 0 ? size_t( npts - 1 ) : 0 };

    // special case n=2 -- use linear interpolation.
    {
      real_type hL = X[1] - X[0];
      real_type SL = (Y[1] - Y[0])/hL;
      switch ( npts ) {
      case 2:
        Yp[0] = Yp[1] = SL;
        return;
      case 3:
        {
          real_type hR = X[2] - X[1];
          real_type SR = (Y[2] - Y[1])/hR;
          Yp[0] = first_deriv3p_L( SL, hL, SR, hR );
          Yp[1] = first_deriv3p_C( SL, hL, SR, hR );
          Yp[2] = first_deriv3p_R( SL, hL, SR, hR );
        }
        return;
      case 4:
        {
          real_type hC = X[2] - X[1];
          real_type SC = (Y[2] - Y[1])/hC;
          real_type hR = X[3] - X[2];
          real_type SR = (Y[3] - Y[2])/hR;
          Yp[0] = first_deriv4p_L( SL, hL, SC, hC, SR, hR );
          Yp[1] = first_deriv3p_C( SL, hL, SC, hC );
          Yp[2] = first_deriv3p_C( SC, hC, SR, hR );
          Yp[3] = first_deriv4p_R( SL, hL, SC, hC, SR, hR );
        }
        return;
      }
    }

    // loop through interior points.
    for ( size_t i{2}; i < n-1 ; ++i ) {

      real_type hLL = X[i-1] - X[i-2];
      real_type hL  = X[i+0] - X[i-1];
      real_type hR  = X[i+1] - X[i+0];
      real_type hRR = X[i+2] - X[i+1];

      real_type SLL = (Y[i-1]-Y[i-2])/hLL;
      real_type SL  = (Y[i+0]-Y[i-1])/hL;
      real_type SR  = (Y[i+1]-Y[i+0])/hR;
      real_type SRR = (Y[i+2]-Y[i+1])/hRR;

      Yp[i] = first_deriv5p_C( SLL, hLL, SL, hL, SR, hR, SRR, hRR );
    }
    {
      real_type hL  = X[n]   - X[n-1];
      real_type hLL = X[n-1] - X[n-2];
      real_type SL  = (Y[n]   - Y[n-1])/hL;
      real_type SLL = (Y[n-1] - Y[n-2])/hLL;
      deriv_right( SL, hL, SLL, hLL, Yp[n-2], Yp[n], Yp[n-1] );
    }
    {
      real_type hR  = X[1] - X[0];
      real_type hRR = X[2] - X[1];
      real_type SR  = (Y[1] - Y[0])/hR;
      real_type SRR = (Y[2] - Y[1])/hRR;
      deriv_left( SR, hR, SRR, hRR, Yp[2], Yp[0], Yp[1] );
    }

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  second_derivative_build(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    real_type       Ypp[],
    integer         npts
  ) {

    size_t n{ npts > 0 ? size_t( npts - 1 ) : 0 };

    // special case n=2 -- use linear interpolation.
    switch ( npts ) {
    case 2:
      Ypp[0] = Ypp[1] = 0;
      return;
    }

    // loop through interior points.
    for ( size_t i{1}; i < n ; ++i ) {

      real_type hL = X[i+0] - X[i-1];
      real_type hR = X[i+1] - X[i+0];
      real_type SL = (Y[i+0]-Y[i-1])/hL;
      real_type SR = (Y[i+1]-Y[i+0])/hR;

      real_type dp0 = Yp[i];
      real_type dpL = Yp[i-1];
      real_type dpR = Yp[i+1];

      Ypp[i] = second_deriv3p_C( SL, hL, SR, hR, dp0, dpL, dpR );
    }
    {
      real_type hL  = X[1] - X[0];
      real_type hR  = X[2] - X[1];
      real_type SL  = (Y[1]-Y[0])/hL;
      real_type SR  = (Y[2]-Y[1])/hR;
      real_type dp0 = Yp[1];
      real_type dpL = Yp[0];
      real_type dpR = Yp[2];
      Ypp[0] = second_deriv3p_L( SL, hL, SR, hR, dp0, dpL, dpR );
    }
    {
      real_type hL  = X[n-1] - X[n-2];
      real_type hR  = X[n] - X[n-1];
      real_type SL  = (Y[n-1]-Y[n-2])/hL;
      real_type SR  = (Y[n]-Y[n-1])/hR;
      real_type dp0 = Yp[n-1];
      real_type dpL = Yp[n-2];
      real_type dpR = Yp[n];
      Ypp[n] = second_deriv3p_R( SL, hL, SR, hR, dp0, dpL, dpR );
    }

  }

  #endif

}
