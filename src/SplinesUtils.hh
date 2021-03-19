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
/****************************************************************************\
Copyright (c) 2016, Enrico Bertolazzi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
\****************************************************************************/
#pragma once

#ifndef SPLINES_UTILS_HH
#define SPLINES_UTILS_HH

#include "Splines.hh"

//! Various kind of splines
namespace Splines {

  using std::fpclassify;

  /*
  //    __ _       _ _             _ _  __  __
  //   / _(_)_ __ (_) |_ ___    __| (_)/ _|/ _| ___ _ __ ___ _ __   ___ ___
  //  | |_| | '_ \| | __/ _ \  / _` | | |_| |_ / _ \ '__/ _ \ '_ \ / __/ _ \
  //  |  _| | | | | | ||  __/ | (_| | |  _|  _|  __/ | |  __/ | | | (_|  __/
  //  |_| |_|_| |_|_|\__\___|  \__,_|_|_| |_|  \___|_|  \___|_| |_|\___\___|
  */

  /*\
   |  Given
   |  DL = (y[k]-y[k-1])/hL; hL = x[k]-x[k-1];
   |  DR = (y[k+1]-y[k])/hR; hR = x[k+1]-x[k];
   |  approximate the first derivative at x[k]
   */
  static
  inline
  real_type
  first_deriv3p_C(
    real_type DL, real_type hL,
    real_type DR, real_type hR
  ) {
    return (DL*hR+DR*hL)/(hR+hL);
  }

  /*\
   |  Given
   |  DL = (y[k]-y[k-1])/hL; hL = x[k]-x[k-1];
   |  DR = (y[k+1]-y[k])/hR; hR = x[k+1]-x[k];
   |  approximate the first derivative at x[k-1]
   */
  static
  inline
  real_type
  first_deriv3p_L(
    real_type DL, real_type hL,
    real_type DR, real_type hR
  ) {
    return ((2*DL-DR)*hL+DL*hR)/(hR+hL);
  }

  /*\
   |  Given
   |  DL = (y[k]-y[k-1])/hL; hL = x[k]-x[k-1];
   |  DR = (y[k+1]-y[k])/hR; hR = x[k+1]-x[k];
   |  approximate the first derivative at x[k+1]
   */
  static
  inline
  real_type
  first_deriv3p_R(
    real_type DL, real_type hL,
    real_type DR, real_type hR
  ) {
    return ((2*DR-DL)*hR+DR*hL)/(hR+hL);
  }

  /*\
   |  Given
   |  DLLL = (y[k-2]-y[k-3])/hLLL; hLLL = x[k-2]-x[k-3];
   |  DLL  = (y[k-1]-y[k-2])/hLL;  hLL  = x[k-1]-x[k-2];
   |  DL   = (y[k]-y[k-1])/hL;     hL   = x[k]-x[k-1];
   |  approximate the first derivative at x[k]
   */
  real_type
  first_deriv4p_R(
    real_type DLLL, real_type hLLL,
    real_type DLL,  real_type hLL,
    real_type DL,   real_type hL
  );

  /*\
   |  Given
   |  DR   = (y[k+1]-y[k])/hR;     hR   = x[k+1]-x[k];
   |  DRR  = (y[k+2]-y[k+1])/hRR;  hRR  = x[k+2]-x[k+1];
   |  DRRR = (y[k+3]-y[k+2])/hRRR; hRRR = x[k+3]-x[k+2];
   |  approximate the first derivative at x[k]
   */
  real_type
  first_deriv4p_L(
    real_type DR,   real_type hR,
    real_type DRR,  real_type hRR,
    real_type DRRR, real_type hRRR
  );

  /*\
   |  Given
   |  DLL = (y[k-1]-y[k-2])/hLL; hLL = x[k-1]-x[k-2];
   |  DL  = (y[k]-y[k-1])/hL;    hL  = x[k]-x[k-1];
   |  DR  = (y[k+1]-y[k])/hR;    hR  = x[k+1]-x[k];
   |  DRR = (y[k+2]-y[k-1])/hRR; hRR = x[k+2]-x[k+1];
   |  approximate the first derivative at x[k]
   */
  real_type
  first_deriv5p_C(
    real_type DLL, real_type hLL,
    real_type DL,  real_type hL,
    real_type DR,  real_type hR,
    real_type DRR, real_type hRR
  );

  /*\
   |  Given
   |  DLL = (y[k-1]-y[k-2])/hLL; hLL = x[k-1]-x[k-2];
   |  DL  = (y[k]-y[k-1])/hL;    hL  = x[k]-x[k-1];
   |  DR  = (y[k+1]-y[k])/hR;    hR  = x[k+1]-x[k];
   |  DRR = (y[k+2]-y[k-1])/hRR; hRR = x[k+2]-x[k+1];
   |  approximate the first derivative at x[k-1]
   */
  real_type
  first_deriv5p_L(
    real_type DLL, real_type hLL,
    real_type DL,  real_type hL,
    real_type DR,  real_type hR,
    real_type DRR, real_type hRR
  );

  /*\
   |  Given
   |  DLL = (y[k-1]-y[k-2])/hLL; hLL = x[k-1]-x[k-2];
   |  DL  = (y[k]-y[k-1])/hL;    hL  = x[k]-x[k-1];
   |  DR  = (y[k+1]-y[k])/hR;    hR  = x[k+1]-x[k];
   |  DRR = (y[k+2]-y[k-1])/hRR; hRR = x[k+2]-x[k+1];
   |  approximate the first derivative at x[k+1]
   */
  real_type
  first_deriv5p_R(
    real_type DLL, real_type hLL,
    real_type DL,  real_type hL,
    real_type DR,  real_type hR,
    real_type DRR, real_type hRR
  );

  real_type
  second_deriv3p_C(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dpL,
    real_type dp0,
    real_type dpR
  );

  real_type
  second_deriv3p_C(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dp0
  );

  real_type
  second_deriv3p_L(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dpL,
    real_type dp0,
    real_type dpR
  );

  real_type
  second_deriv3p_L(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dpL
  );

  real_type
  second_deriv3p_R(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dpL,
    real_type dp0,
    real_type dpR
  );

  real_type
  second_deriv3p_R(
    real_type SL,
    real_type hL,
    real_type SR,
    real_type hR,
    real_type dpR
  );

  void
  first_derivative_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts
  );

  void
  second_derivative_build(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    real_type       Ypp[],
    integer         npts
  );

}

#endif
