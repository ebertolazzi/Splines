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
 |   ____                     _ ____        _ _
 |  | __ )  ___  ___ ___  ___| / ___| _ __ | (_)_ __   ___
 |  |  _ \ / _ \/ __/ __|/ _ \ \___ \| '_ \| | | '_ \ / _ \
 |  | |_) |  __/\__ \__ \  __/ |___) | |_) | | | | | |  __/
 |  |____/ \___||___/___/\___|_|____/| .__/|_|_|_| |_|\___|
 |                                   |_|
\*/

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

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  void
  Bessel_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer         npts
  ) {

    size_t n{ size_t(npts > 0 ? npts-1 : 0) };

    Malloc_real mem("Bessel_build");
    real_type * m{ mem.malloc( size_t(n+1) ) };

    // calcolo slopes
    for ( size_t i{0}; i < n; ++i )
      m[i] = (Y[i+1]-Y[i])/(X[i+1]-X[i]);

    if ( npts == 2 ) { // caso speciale 2 soli punti

      Yp[0] = Yp[1] = m[0];

    } else {

      for ( size_t i{1}; i < n; ++i ) {
        real_type DL = X[i]   - X[i-1];
        real_type DR = X[i+1] - X[i];
        Yp[i] = (DR*m[i-1]+DL*m[i])/((DL+DR));
      }

      Yp[0] = 1.5*m[0]-0.5*m[1];
      Yp[n] = 1.5*m[n-1]-0.5*m[n-2];
    }
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BesselSpline::build() {
    string msg{ fmt::format("BesselSpline[{}]::build():", m_name ) };
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
      Bessel_build( m_X+ibegin, m_Y+ibegin, m_Yp+ibegin, iend-ibegin );
      ibegin = iend;
    } while ( iend < m_npts );

    Utils::check_NaN( m_Yp, (msg+" Yp").c_str(), m_npts, __LINE__, __FILE__ );
  }

  using GC_namespace::GC_type;
  using GC_namespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BesselSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string where{ fmt::format("BesselSpline[{}]::setup( gc ):", m_name ) };
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
    this->build( x, y );
  }

}
