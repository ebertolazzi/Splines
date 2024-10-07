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

#include "Splines.hh"
#include "SplinesUtils.hh"
#include "Utils_fmt.hh"

#include <cmath>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namspace
#endif

namespace Splines {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  QuinticSplineBase::QuinticSplineBase( string const & name )
  : Spline(name)
  , m_base_quintic( fmt::format( "QuinticSplineBase[{}]", name ) )
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
     Sistema lineare da risolvere

     D U
     L D U
       L D U
         L D U
           .....
              L D U
                L D U
                  L D

  \*/

  static
  void
  QuinticSpline_Yppp_continuous(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    real_type       Ypp[],
    integer         npts,
    bool            setbc
  ) {

    size_t n{ size_t(npts > 0 ? npts-1 : 0) };

    Malloc_real mem("QuinticSpline_Yppp_continuous");
    mem.allocate( size_t(3*(n+1)) );
    real_type * L = mem( size_t( n+1 ) );
    real_type * D = mem( size_t( n+1 ) );
    real_type * U = mem( size_t( n+1 ) );
    real_type * Z = Ypp;

    size_t i;
    for ( i = 1; i < n; ++i ) {
      real_type hL  = X[i] - X[i-1];
      real_type hL2 = hL*hL;
      real_type hL3 = hL*hL2;
      real_type hR  = X[i+1] - X[i];
      real_type hR2 = hR*hR;
      real_type hR3 = hR*hR2;
      real_type DL  = 60*(Y[i]-Y[i-1])/hL3;
      real_type DR  = 60*(Y[i+1]-Y[i])/hR3;
      real_type DDL = (36*Yp[i]+24*Yp[i-1])/hL2;
      real_type DDR = (36*Yp[i]+24*Yp[i+1])/hR2;
      L[i] = -3/hL;
      D[i] = 9/hL+9/hR;
      U[i] = -3/hR;
      Z[i] = DR-DL+DDL-DDR;
    }
    L[0] = U[0] = 0; D[0] = 1;
    L[n] = U[n] = 0; D[n] = 1;
    if ( setbc ) {
      {
        real_type hL = X[1] - X[0];
        real_type hR = X[2] - X[1];
        real_type SL = (Y[1]-Y[0])/hL;
        real_type SR = (Y[2]-Y[1])/hR;
        real_type dp0 = Yp[1];
        real_type dpL = Yp[0];
        real_type dpR = Yp[2];
        Z[0] = second_deriv3p_L( SL, hL, SR, hR, dpL, dp0, dpR );
      }
      {
        real_type hL = X[n-1] - X[n-2];
        real_type hR = X[n] - X[n-1];
        real_type SL = (Y[n-1]-Y[n-2])/hL;
        real_type SR = (Y[n]-Y[n-1])/hR;
        real_type dp0 = Yp[n-1];
        real_type dpL = Yp[n-2];
        real_type dpR = Yp[n];
        Z[n] = second_deriv3p_R( SL, hL, SR, hR, dpL, dp0, dpR );
      }
    }

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

  static
  void
  QuinticSpline_Ypp_build(
    real_type const X[],
    real_type const Y[],
    real_type const Yp[],
    real_type       Ypp[],
    integer         npts
  ) {

    size_t n{ size_t(npts > 0 ? npts-1 : 0) };

    if ( n == 1 ) { Ypp[0] = Ypp[1] = 0; return; }

    {
      real_type hL = X[1] - X[0];
      real_type hR = X[2] - X[1];
      real_type SL = (Y[1] - Y[0])/hL;
      real_type SR = (Y[2] - Y[1])/hR;
      //Ypp[0] = (2*SL-SR)*al+SL*be;
      //Ypp[0] = second_deriv3p_L( SL, hL, SR, hR, Yp[0] );
      Ypp[0] = second_deriv3p_L( SL, hL, SR, hR, Yp[0], Yp[1], Yp[2] );
    }
    {
      real_type hL = X[n-1] - X[n-2];
      real_type hR = X[n] - X[n-1];
      real_type SL = (Y[n-1] - Y[n-2])/hL;
      real_type SR = (Y[n] - Y[n-1])/hR;
      //Ypp[n] = (2*SR-SL)*be+SR*al;
      //Ypp[n] = second_deriv3p_R( SL, hL, SR, hR, Yp[n] );
      Ypp[n] = second_deriv3p_R( SL, hL, SR, hR, Yp[n-2], Yp[n-1], Yp[n] );
    }

    size_t i;
    for ( i = 1; i < n; ++i ) {
      real_type hL = X[i] - X[i-1];
      real_type hR = X[i+1] - X[i];
      real_type SL = (Y[i] - Y[i-1])/hL;
      real_type SR = (Y[i+1] - Y[i])/hR;
      //Ypp[i] = second_deriv3p_C( SL, hL, SR, hR, Yp[i] );
      real_type ddC = second_deriv3p_C( SL, hL, SR, hR, Yp[i-1], Yp[i], Yp[i+1] );
      if ( i > 1 ) {
        real_type hLL = X[i-1] - X[i-2];
        real_type SLL = (Y[i-1] - Y[i-2])/hLL;
        //real_type dd  = second_deriv3p_R( SLL, hLL, SL, hR, Yp[i] );
        real_type ddL = second_deriv3p_R( SLL, hLL, SL, hR, Yp[i-2], Yp[i-1], Yp[i] );
        if      ( ddL * ddC < 0 ) ddC = 0;
        else if ( std::abs(ddL) < std::abs(ddC) ) ddC = ddL;
      }
      if ( i < n-1 ) {
        real_type hRR = X[i+2] - X[i+1];
        real_type SRR = (Y[i+2] - Y[i+1])/hRR;
        //real_type dd = second_deriv3p_L( SR, hR, SRR, hRR, Yp[i] );
        real_type ddR = second_deriv3p_L( SR, hR, SRR, hRR, Yp[i], Yp[i+1], Yp[i+2] );
        if      ( ddR * ddC < 0 ) ddC = 0;
        else if ( std::abs(ddR) < std::abs(ddC) ) ddC = ddR;
      }
      Ypp[i] = ddC;
    }
  }

  /*\
   |    ___        _       _   _      ____        _ _
   |   / _ \ _   _(_)_ __ | |_(_) ___/ ___| _ __ | (_)_ __   ___
   |  | | | | | | | | '_ \| __| |/ __\___ \| '_ \| | | '_ \ / _ \
   |  | |_| | |_| | | | | | |_| | (__ ___) | |_) | | | | | |  __/
   |   \__\_\\__,_|_|_| |_|\__|_|\___|____/| .__/|_|_|_| |_|\___|
   |                                       |_|
   |
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  Quintic_build(
    QuinticSpline_sub_type q_sub_type,
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    real_type       Ypp[],
    integer         npts
  ) {
    switch ( q_sub_type ) {
    case QuinticSpline_sub_type::CUBIC:
      {
        size_t n{ size_t(npts > 0 ? npts-1 : 0) };

        Malloc_real mem("QuinticSpline_Yppp_continuous");
        mem.allocate( size_t(3*(n+1)) );
        real_type * L = mem( size_t( n+1 ) );
        real_type * D = mem( size_t( n+1 ) );
        real_type * U = mem( size_t( n+1 ) );
        CubicSpline_build(
          X, Y, Yp, Ypp, L, D, U, npts,
          CubicSpline_BC::EXTRAPOLATE,
          CubicSpline_BC::EXTRAPOLATE
        );
        mem.free();
        QuinticSpline_Yppp_continuous( X, Y, Yp, Ypp, npts, false );
      }
      return;
    case QuinticSpline_sub_type::PCHIP:
      Pchip_build( X, Y, Yp, npts );
      break;
    case QuinticSpline_sub_type::AKIMA:
      Akima_build( X, Y, Yp, npts );
      break;
    case QuinticSpline_sub_type::BESSEL:
      Bessel_build( X, Y, Yp, npts );
      break;
    }
    QuinticSpline_Ypp_build( X, Y, Yp, Ypp, npts );
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSpline::build() {
    string msg{ fmt::format("QuinticSpline[{}]::build():", m_name ) };
    UTILS_ASSERT(
      m_npts > 1,
      "{} npts = {} not enought points\n",
      msg, m_npts
    );
    Utils::check_NaN( m_X, (msg+" X").c_str(), m_npts, __LINE__, __FILE__ );
    Utils::check_NaN( m_Y, (msg+" Y").c_str(), m_npts, __LINE__, __FILE__ );
    integer ibegin{0};
    integer iend{0};
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < m_npts && m_X[iend-1] < m_X[iend] ) {}
      Quintic_build(
        m_q_sub_type,
        m_X+ibegin,  m_Y+ibegin,
        m_Yp+ibegin, m_Ypp+ibegin,
        iend - ibegin
      );
      ibegin = iend;
    } while ( iend < m_npts );

    Utils::check_NaN( m_Yp,  (msg+" Yp").c_str(),  m_npts, __LINE__, __FILE__ );
    Utils::check_NaN( m_Ypp, (msg+" Ypp").c_str(), m_npts, __LINE__, __FILE__ );
  }

  using GC_namespace::GC_type;
  using GC_namespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  QuinticSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string where{ fmt::format("QuinticSpline[{}]::setup( gc ):", m_name ) };
    GenericContainer const & gc_x{ gc("xdata",where.c_str()) };
    GenericContainer const & gc_y{ gc("ydata",where.c_str()) };

    vec_real_type x, y;
    {
      std::string ff{ fmt::format( "{}, field `xdata'", where ) };
      gc_x.copyto_vec_real ( x, ff.c_str() );
    }
    {
      std::string ff{ fmt::format( "{}, field `ydata'", where ) };
      gc_y.copyto_vec_real ( y, ff.c_str() );
    }
    if ( gc.exists("spline_sub_type") ) {
      std::string const & st{ gc.get_map_string("spline_sub_type",where.c_str()) };
      if      ( st == "cubic"  ) m_q_sub_type = QuinticSpline_sub_type::CUBIC;
      else if ( st == "pchip"  ) m_q_sub_type = QuinticSpline_sub_type::PCHIP;
      else if ( st == "akima"  ) m_q_sub_type = QuinticSpline_sub_type::AKIMA;
      else if ( st == "bessel" ) m_q_sub_type = QuinticSpline_sub_type::BESSEL;
      else {
        UTILS_ERROR( "{} unknow sub type: {}\n", where, st );
      }
    } else {
      UTILS_WARNING( false,
        "{}, missing field `spline_sub_type` using `cubic` as default value\n", where
      );
    }
    this->build( x, y );
  }

}
