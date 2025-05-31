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

/**
 *
 */
/*\
 |     #
 |    # #   #    # # #    #   ##
 |   #   #  #   #  # ##  ##  #  #
 |  #     # ####   # # ## # #    #
 |  ####### #  #   # #    # ######
 |  #     # #   #  # #    # #    #
 |  #     # #    # # #    # #    #
\*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
using namespace std; // load standard namspace
#endif

namespace Splines {

  using std::abs;

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Akima_build(
    real_type const X[],
    real_type const Y[],
    real_type       Yp[],
    real_type       m[],   // work vector dimension N
    integer         N
  ) {

    UTILS_ASSERT( N >= 2, "Akima_build require at least 2 points" );

    if ( N == 2 ) { // solo 2 punti, niente da fare
      Yp[0] = Yp[1] = (Y[1]-Y[0])/(X[1]-X[0]);
      return;
    }
    
    // 1. Calcola le pendenze m_i = (Y[i+1] - Y[i]) / (X[i+1] - X[i])
    for ( integer i{0}; i < N-1; ++i ) {
      UTILS_ASSERT(
        X[i+1] > X[i],
        "Akima_build, X must be strictly increasing X[{}] = {}, X[{}] = {}",
        i, X[i], i+1, X[i+1]
      );
      m[i] = (Y[i+1] - Y[i]) / (X[i+1] - X[i]);
    }

    // 2. Calcolo epsi
    real_type epsi{0};
    for ( integer i{0}; i < N-2; ++i ) {
      real_type const dm{ std::abs(m[i+1]-m[i]) };
      if ( dm > epsi ) epsi = dm;
    }
    epsi *= 1E-8;

    // 3. Calcola le derivate nei punti interni (i = 1..N-2)
    for ( integer i{1}; i < N-1; ++i ) {
      // Estrapola le pendenze se necessario (ai bordi)
      real_type const m_im2 { (i >= 2) ? m[i-2] : 2*m[i-1] - m[i]  }; // m_{i-2}
      real_type const m_im1 { m[i-1]                               }; // m_{i-1}
      real_type const m_i   { m[i]                                 }; // m_i
      real_type const m_ip1 { (i < N-2) ? m[i+1] : 2*m[i] - m[i-1] }; // m_{i+1}

      // Pesatura di Akima (come in MATLAB makima)
      // https://blogs.mathworks.com/cleve/2019/04/29/makima-piecewise-cubic-interpolation/#d9a97978-0b09-4a1f-a6a5-504d088631d0
      real_type const w_left  { std::abs( m_ip1 - m_i )   + std::abs( m_ip1 + m_i   ) / 2 };
      real_type const w_right { std::abs( m_im1 - m_im2 ) + std::abs( m_im1 + m_im2 ) / 2 };
      real_type const sum_w   { w_left + w_right };

      if ( sum_w > epsi ) {
        Yp[i] = (w_right * m_im1 + w_left * m_i) / sum_w;
      } else {
        Yp[i] = 0.5*(m_im1 + m_i);  // Caso speciale
      }
    }

    // 4. Derivate ai bordi (i=0 e i=N-1)
    // Estrapolazione quadratica
    Yp[0]   = m[0] + (m[0] - m[1]) * (X[0] - X[1]) / (X[2] - X[0]);
    Yp[N-1] = m[N-2] + (m[N-2] - m[N-3]) * (X[N-1] - X[N-2]) / (X[N-1] - X[N-3]);
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  AkimaSpline::build() {
    string const msg{ fmt::format("AkimaSpline[{}]::build():", m_name ) };
    UTILS_ASSERT( m_npts > 1, "{} npts = {} not enought points\n", msg, m_npts );
    Utils::check_NaN( m_X, msg+" X ", m_npts, __LINE__, __FILE__ );
    Utils::check_NaN( m_Y, msg+" Y ", m_npts, __LINE__, __FILE__ );
    integer ibegin{0};
    integer iend{0};

    Malloc_real m_mem("AkimaSpline::work memory");
    real_type * m_work{ m_mem.malloc( m_npts ) };

    do {
      // cerca intervallo monotono strettamente crescente
      for ( ++iend; iend < m_npts && m_X[iend-1] < m_X[iend]; ++iend ) {}
      Akima_build( m_X+ibegin, m_Y+ibegin, m_Yp+ibegin, m_work, iend-ibegin );
      ibegin = iend;
    } while ( iend < m_npts );

    Utils::check_NaN( m_Yp, msg+" Yp", m_npts, __LINE__, __FILE__ );
    m_search.must_reset();
  }

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  using GC_namespace::GC_type;
  using GC_namespace::vec_real_type;
  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  AkimaSpline::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where{ fmt::format("AkimaSpline[{}]::setup():", m_name ) };

    GenericContainer const & gc_x{ gc("xdata",where) };
    GenericContainer const & gc_y{ gc("ydata",where) };

    vec_real_type x, y;
    {
      string const ff{ fmt::format( "{}, field `xdata'", where ) };
      gc_x.copyto_vec_real( x, ff );
    }
    {
      string const ff{ fmt::format( "{}, field `ydata'", where ) };
      gc_y.copyto_vec_real( y, ff );
    }
    this->build( x, y );
  }

}
