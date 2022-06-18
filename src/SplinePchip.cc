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
#include "SplinesUtils.hh"

#include <cmath>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namspace
#endif

namespace Splines {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  typedef enum {
    region_A = 0,
    region_B,
    region_C,
    region_D,
    region_E,
    region_M
  } REGION_ABCDEM;

  static
  REGION_ABCDEM
  get_region( real_type alpha, real_type beta ) {
    // assuming alpha >= 0, beta >= 0
    real_type apb = alpha+beta;
    if ( apb <= 3 ) return region_M;
    real_type r = alpha * alpha + beta * beta + alpha * beta - 6*apb + 9;
    if ( r <= 0 ) return region_M;
    if ( apb <= 4 ) {
      if ( beta >= alpha ) return region_A;
      else                 return region_E;
    }
    if ( alpha >= 3 && beta >= 3 ) return region_C;
    if ( beta >= alpha ) return region_B;
    else                 return region_D;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  inline
  real_type
  max_abs( real_type a, real_type b ) {
    real_type res = std::abs(a);
    if ( res < std::abs(b) ) res = std::abs(b);
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  inline
  real_type
  min_abs( real_type a, real_type b ) {
    real_type res = std::abs(a);
    if ( res > std::abs(b) ) res = std::abs(b);
    return res;
  }

  /*
  //   ____      _     _      ____        _ _
  //  |  _ \ ___| |__ (_)_ __/ ___| _ __ | (_)_ __   ___
  //  | |_) / __| '_ \| | '_ \___ \| '_ \| | | '_ \ / _ \
  //  |  __/ (__| | | | | |_) |__) | |_) | | | | | |  __/
  //  |_|   \___|_| |_|_| .__/____/| .__/|_|_|_| |_|\___|
  //                    |_|        |_|
  */

  static
  inline
  int
  signTest( real_type const a, real_type const b ) {
    int sa = 0;
    int sb = 0;
    if      ( a > 0 ) sa =  1;
    else if ( a < 0 ) sa = -1;
    if      ( b > 0 ) sb =  1;
    else if ( b < 0 ) sb = -1;
    return sa*sb;
  }

  //!
  //! References:
  //! ==========
  //!
  //! F.N. Fritsch, R.E. Carlson:
  //! Monotone Piecewise Cubic Interpolation,
  //! SIAM J. Numer. Anal. Vol 17, No. 2, April 1980
  //!
  //! F.N. Fritsch and J. Butland:
  //! A method for constructing local monotone piecewise cubic interpolants,
  //! SIAM Journal on Scientific and Statistical Computing 5, 2 (June 1984), pp. 300-304.
  //!
  void
  Pchip_build(
    real_type const * X,
    real_type const * Y,
    real_type       * Yp,
    integer           npts
  ) {

    size_t n = npts > 0 ? size_t( npts - 1 ) : 0;

    integer ierr = 0;

    // function definition is ok, go on.
    real_type h1    = X[1] - X[0];
    real_type del1  = (Y[1]-Y[0])/h1;
    real_type dsave = del1;

    // special case n=2 -- use linear interpolation.
    if ( n == 1 ) { Yp[0] = Yp[1] = del1; return; }

    real_type h2   = X[2] - X[1];
    real_type del2 = (Y[2]-Y[1])/h2;

    // Set Yp[0] via non-centered three-point formula, adjusted to be shape-preserving.
    real_type hsum = h1 + h2;
    real_type w1   = (h1 + hsum)/hsum;
    real_type w2   = -h1/hsum;
    Yp[0] = w1*del1 + w2*del2;
    real_type dmin, dmax;
    if ( signTest(Yp[0],del1) <= 0 ) {
      Yp[0] = 0;
    } else if ( signTest(del1,del2) < 0 ) {
      // NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
      dmax = 3*del1;
      if ( std::abs(Yp[0]) > std::abs(dmax) ) Yp[0] = dmax;
    }

    // loop through interior points.
    for ( size_t i = 1; i < n; ++i ) {
      if ( i > 1 ) {
        h1   = h2;
        h2   = X[i+1] - X[i];
        hsum = h1 + h2;
        del1 = del2;
        del2 = (Y[i+1] - Y[i])/h2;
      }
      // set Yp[i]=0 unless data are strictly monotonic.
      Yp[i] = 0;
      // count number of changes in direction of monotonicity.
      switch ( signTest(del1,del2) ) {
      case -1:
        if ( Utils::is_zero(del2) ) break;
        if ( signTest(dsave,del2) < 0 ) ++ierr;
        dsave = del2;
        break;
      case 0:
        ++ierr;
        dsave = del2;
        break;
      case 1: // use brodlie modification of butland formula.
        w1   = (1+h1/hsum)/3;
        w2   = (1+h2/hsum)/3;
        dmax = max_abs( del1, del2 );
        dmin = min_abs( del1, del2 );
        real_type drat1 = del1/dmax;
        real_type drat2 = del2/dmax;
        Yp[i] = dmin/(w1*drat1 + w2*drat2);
        break;
      }
    }
    // set Yp[n] via non-centered three-point formula, adjusted to be shape-preserving.
    w1 = -h2/hsum;
    w2 = (h2 + hsum)/hsum;
    Yp[n] = w1*del1 + w2*del2;
    if ( signTest(Yp[n],del2) <= 0 ) {
      Yp[n] = 0;
    } else if ( signTest(del1,del2) < 0 ) {
      // need do this check only if monotonicity switches.
      dmax = 3*del2;
      if ( abs(Yp[n]) > abs(dmax) ) Yp[n] = dmax;
    }
    // cout << "ierr = " << ierr << '\n';
  }

  static // non usata per ora
  void
  Pchip_build_new(
    real_type const * X,
    real_type const * Y,
    real_type       * Yp,
    integer           npts
  ) {

    first_derivative_build( X, Y, Yp, npts );

    size_t n = npts > 0 ? size_t( npts - 1 ) : 0;

    // loop through interior points.
    for ( size_t i = 1; i < n ; ++i ) {

      real_type hL = X[i+0] - X[i-1];
      real_type hR = X[i+1] - X[i+0];

      real_type SL = (Y[i+0]-Y[i-1])/hL;
      real_type SR = (Y[i+1]-Y[i+0])/hR;

      real_type fp = Yp[i];

      real_type sigma = 0;
      if ( SL*SR > 0 ) sigma = SR > 0 ? 1 : -1;
      real_type absSL = SL < 0 ? -SL: SL;
      real_type absSR = SR < 0 ? -SR: SR;
      real_type Delta = 3*( absSL < absSR ? absSL : absSR );
      if ( sigma > 0 ) {
        if ( fp < 0     ) fp = 0;
        if ( fp > Delta ) fp = Delta;
      } else if ( sigma < 0 ) {
        if ( fp > 0      ) fp = 0;
        if ( fp < -Delta ) fp = -Delta;
      } else {
        fp = 0;
      }
      Yp[i] = fp;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PchipSpline::build() {
    string msg = fmt::format("PchipSpline[{}]::build():", m_name );
    UTILS_ASSERT(
      m_npts > 1,
      "{} npts = {} not enought points\n",
      msg, m_npts
    );
    Utils::check_NaN( m_X, (msg+" X").c_str(), m_npts, __LINE__, __FILE__ );
    Utils::check_NaN( m_Y, (msg+" Y").c_str(), m_npts, __LINE__, __FILE__ );
    integer ibegin = 0;
    integer iend   = 0;
    do {
      // cerca intervallo monotono strettamente crescente
      while ( ++iend < m_npts && m_X[iend-1] < m_X[iend] ) {}
      Pchip_build(
        m_X+ibegin,
        m_Y+ibegin,
        m_Yp+ibegin,
        iend-ibegin
      );
      ibegin = iend;
    } while ( iend < m_npts );

    Utils::check_NaN( m_Yp, (msg+" Yp").c_str(), m_npts, __LINE__, __FILE__ );
  }

  using GC_namespace::GC_VEC_REAL;
  using GC_namespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PchipSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string msg = fmt::format("PchipSpline[{}]::setup( gc ):", m_name );
    UTILS_ASSERT( gc.exists("xdata"), "{} missing `xdata` field!\n", msg );
    UTILS_ASSERT( gc.exists("ydata"), "{} missing `ydata` field!\n", msg );

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
    this->build( x, y );
  }

}
