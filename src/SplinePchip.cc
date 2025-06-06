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
 |      Università degli Studi di Trento                                    |
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
#include "Utils_fmt.hh"

#include <cmath>
#include <set>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namspace
#endif

namespace Splines {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if 0
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
  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  real_type
  max_abs( real_type const a, real_type const b ) {
    real_type res{ std::abs(a) };
    if ( res < std::abs(b) ) res = std::abs(b);
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  real_type
  min_abs( real_type const a, real_type const b ) {
    real_type res{ std::abs(a) };
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
  int
  signTest( real_type const a, real_type const b ) {
    int sa{0};
    int sb{0};
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
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer   const npts
  ) {

    UTILS_ASSERT( npts >= 2, "Pchip_build, npts={} must be >= 2\n", npts );

    integer const n{ npts - 1 };

    // function definition is ok, go on.
    real_type h1   { X[1] - X[0]    };
    real_type del1 { (Y[1]-Y[0])/h1 };

    // special case n=2 -- use linear interpolation.
    if ( n == 1 ) { Yp[0] = Yp[1] = del1; return; }

    real_type h2   { X[2] - X[1]    };
    real_type del2 { (Y[2]-Y[1])/h2 };

    // Set Yp[0] via non-centered three-point formula, adjusted to be shape-preserving.
    real_type hsum { h1 + h2 };
    real_type w1   { 1+h1/hsum };
    real_type w2   { -h1/hsum };
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
    for ( integer i{1}; i < n; ++i ) {
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
        break;
      case 0:
        break;
      case 1: // use brodlie modification of butland formula.
        w1   = (1+h1/hsum)/3;
        w2   = (1+h2/hsum)/3;
        dmax = max_abs( del1, del2 );
        dmin = min_abs( del1, del2 );
        real_type const drat1{ del1/dmax };
        real_type const drat2{ del2/dmax };
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PchipSpline::build() {
    string msg{ fmt::format("PchipSpline[{}]::build():", m_name ) };
    UTILS_ASSERT( m_npts > 1, "{} npts = {} not enought points\n", msg, m_npts );

    Utils::check_NaN( m_X, msg+" X", m_npts, __LINE__, __FILE__ );
    Utils::check_NaN( m_Y, msg+" Y", m_npts, __LINE__, __FILE__ );

    integer ibegin { 0 };
    integer iend   { 0 };

    do {
      // cerca intervallo monotono strettamente crescente
      for ( ++iend; iend < m_npts && m_X[iend-1] < m_X[iend]; ++iend ) {}
      Pchip_build( m_X+ibegin, m_Y+ibegin, m_Yp+ibegin, iend-ibegin );
      ibegin = iend;
    } while ( iend < m_npts );

    Utils::check_NaN( m_Yp, msg+" Yp", m_npts, __LINE__, __FILE__ );
    m_search.must_reset();
  }

  using GC_namespace::GC_type;
  using GC_namespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PchipSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where{ fmt::format("PchipSpline[{}]::setup( gc ):", m_name ) };

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map(where) ) { keywords.insert(pair.first); }

    GenericContainer const & gc_x{ gc("xdata",where) }; keywords.erase("xdata");
    GenericContainer const & gc_y{ gc("ydata",where) }; keywords.erase("ydata");

    vec_real_type x, y;
    {
      string const ff{ fmt::format( "{}, field `xdata'", where ) };
      gc_x.copyto_vec_real ( x, ff );
    }
    {
      string const ff{ fmt::format( "{}, field `ydata'", where ) };
      gc_y.copyto_vec_real ( y, ff );
    }

    UTILS_WARNING(
      keywords.empty(), "{}: unused keys\n{}\n", where,
      [&keywords]()->string {
        string res;
        for ( auto const & it : keywords ) { res += it; res += ' '; };
        return res;
      }()
    );

    this->build( x, y );
  }

}
