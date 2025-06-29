/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
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

namespace Splines {

  /*
  //    ____  ____   ____                               _
  //   / ___|/ ___| / ___| _   _ _ __  _ __   ___  _ __| |_
  //  | |  _| |     \___ \| | | | '_ \| '_ \ / _ \| '__| __|
  //  | |_| | |___   ___) | |_| | |_) | |_) | (_) | |  | |_
  //   \____|\____| |____/ \__,_| .__/| .__/ \___/|_|   \__|
  //                            |_|   |_|
  */

  using GC_namespace::GC_type;
  using GC_namespace::vec_real_type;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //!
  //! Setup a spline using a `GenericContainer`
  //!
  //! - gc("spline_type")
  //!   - "constant"
  //!   - "linear"
  //!   - "cubic"
  //!   - "akima"
  //!   - "bessel"
  //!   - "pchip"
  //!   - "quintic"
  //!
  void
  Spline2Dblend::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where{ "Spline2Dblend[{}]::setup( gc )" };
    GenericContainer const & gc_surf0{ gc("surf0",where) };
    GenericContainer const & gc_surf1{ gc("surf1",where) };
    m_surf0.setup( gc_surf0 );
    m_surf1.setup( gc_surf1 );
    check_compatibility();
  }

  void
  Spline2Dblend::check_compatibility() const {
    // check compatibility
    UTILS_ASSERT(
      m_surf0.x_min() == m_surf1.x_min(),
      "Spline2Dblend must have the same initial-x, x_min0={}, x_min1={}, difference={}\n",
      m_surf0.x_min(), m_surf1.x_min(), m_surf0.x_min()-m_surf1.x_min()
    );
    UTILS_ASSERT(
      m_surf0.x_max() == m_surf1.x_max(),
      "Spline2Dblend must have the same final-x, x_max0={}, x_max1={}, difference={}\n",
      m_surf0.x_max(), m_surf1.x_max(), m_surf0.x_max()-m_surf1.x_max()
    );
    // check compatibility
    UTILS_ASSERT(
      m_surf0.y_min() == m_surf1.y_min(),
      "Spline2Dblend must have the same initial-y, y_min0={}, y_min1={}, difference={}\n",
      m_surf0.y_min(), m_surf1.y_min(), m_surf0.y_min()-m_surf1.y_min()
    );
    UTILS_ASSERT(
      m_surf0.y_max() == m_surf1.y_max(),
      "Spline2Dblend must have the same final-y, y_max0={}, y_max1={}, difference={}\n",
      m_surf0.y_max(), m_surf1.y_max(), m_surf0.y_max()-m_surf1.y_max()
    );
  }

}
