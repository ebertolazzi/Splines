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
#include <set>

/*\
 |   ____                     _ ____        _ _
 |  | __ )  ___  ___ ___  ___| / ___| _ __ | (_)_ __   ___
 |  |  _ \ / _ \/ __/ __|/ _ \ \___ \| '_ \| | | '_ \ / _ \
 |  | |_) |  __/\__ \__ \  __/ |___) | |_) | | | | | |  __/
 |  |____/ \___||___/___/\___|_|____/| .__/|_|_|_| |_|\___|
 |                                   |_|
\*/

using namespace std; // load standard namspace

namespace Splines {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Bessel_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    integer   const npts
  ) {
  
    UTILS_ASSERT( npts >= 2, "Bessel_build, npts={} must be >= 2\n", npts );

    integer const n{ npts-1 };

    Malloc_real mem("Bessel_build");
    real_type * m{ mem.malloc( npts ) };

    // calcolo slopes
    for ( integer i{0}; i < n; ++i )
      m[i] = (Y[i+1]-Y[i])/(X[i+1]-X[i]);

    if ( npts == 2 ) { // caso speciale 2 o 3 punti

      Yp[0] = Yp[1] = m[0];

    } else if ( npts == 3 ) { // caso speciale 2 o 3 punti

      Yp[0] = m[0];
      Yp[n] = m[1];

    } else {

      for ( integer i{1}; i < n; ++i ) {
        real_type const DL { X[i]   - X[i-1] };
        real_type const DR { X[i+1] - X[i]   };
        Yp[i] = (DR*m[i-1]+DL*m[i])/(DL+DR);
      }

      Yp[0] = 1.5*m[0]-0.5*m[1];
      Yp[n] = 1.5*m[n-1]-0.5*m[n-2];
    }

    mem.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BesselSpline::build() {
    string msg{ fmt::format("BesselSpline[{}]::build():", m_name ) };
    UTILS_ASSERT(
      m_npts > 1,
      "{} npts={} not enought points\n",
      msg, m_npts
    );
    Utils::check_NaN( m_X, msg+" X", m_npts, __LINE__, __FILE__ );
    Utils::check_NaN( m_Y, msg+" Y", m_npts, __LINE__, __FILE__ );
    integer ibegin{0};
    integer iend{0};
    do {
      // cerca intervallo monotono strettamente crescente
      for ( ++iend; iend < m_npts && m_X[iend-1] < m_X[iend]; ++iend ) {}
      Bessel_build( m_X+ibegin, m_Y+ibegin, m_Yp+ibegin, iend-ibegin );
      ibegin = iend;
    } while ( iend < m_npts );

    Utils::check_NaN( m_Yp, msg+" Yp", m_npts, __LINE__, __FILE__ );
    m_search.must_reset();
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
    string const where{ fmt::format("BesselSpline[{}]::setup( gc ):", m_name ) };

    std::set<std::string> keywords;
    for ( auto const & pair : gc.get_map(where) ) { keywords.insert(pair.first); }
    keywords.erase("spline_type");

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
