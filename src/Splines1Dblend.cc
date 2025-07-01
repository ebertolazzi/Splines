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
  Spline1Dblend::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where{ "Spline1Dblend[{}]::setup( gc )" };
    GenericContainer const & gc_spline0{ gc("spline0",where) };
    GenericContainer const & gc_spline1{ gc("spline1",where) };
    m_spline0.setup( gc_spline0 );
    m_spline1.setup( gc_spline1 );
    check_compatibility();
  }

  void
  Spline1Dblend::check_compatibility() const {
    // check compatibility
    UTILS_ASSERT(
      m_spline0.x_min() == m_spline1.x_min(),
      "Spline1Dblend must have the same initial-x, x_min0={}, x_min1={}, difference={}\n",
      m_spline0.x_min(), m_spline1.x_min(), m_spline0.x_min()-m_spline1.x_min()
    );
    UTILS_ASSERT(
      m_spline0.x_max() == m_spline1.x_max(),
      "Spline1Dblend must have the same final-x, x_max0={}, x_max1={}, difference={}\n",
      m_spline0.x_max(), m_spline1.x_max(), m_spline0.x_max()-m_spline1.x_max()
    );
  }

}
