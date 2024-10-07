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
#include "Utils_fmt.hh"

namespace Splines {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  static
  Spline *
  new_Spline1D( string const & _name, SplineType1D tp ) {
    switch ( tp ) {
    case SplineType1D::CONSTANT:   return new ConstantSpline(_name);
    case SplineType1D::LINEAR:     return new LinearSpline(_name);
    case SplineType1D::CUBIC:      return new CubicSpline(_name);
    case SplineType1D::AKIMA:      return new AkimaSpline(_name);
    case SplineType1D::BESSEL:     return new BesselSpline(_name);
    case SplineType1D::PCHIP:      return new PchipSpline(_name);
    case SplineType1D::QUINTIC:    return new QuinticSpline(_name);
    case SplineType1D::HERMITE:    break;
    case SplineType1D::SPLINE_SET: break;
    case SplineType1D::SPLINE_VEC: break;
    }
    return nullptr;
  }

  void
  Spline1D::build(
    SplineType1D tp,
    real_type const x[], integer incx,
    real_type const y[], integer incy,
    integer n
  ) {
    if ( m_pSpline != nullptr ) delete m_pSpline;
    m_pSpline = new_Spline1D(m_name,tp);
    UTILS_ASSERT( m_pSpline != nullptr, "Spline1D[{}]::build, failed\n", m_name );
    m_pSpline->build( x, incx, y, incy, n );
  }

  #endif

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
  Spline1D::setup( GenericContainer const & gc ) {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string where{ fmt::format("Spline1D[{}]::setup( gc ):", m_name ) };
    std::string const & spl_type{ gc.get_map_string("spline_type",where.c_str()) };

    SplineType1D tp;
    if ( spl_type == "constant" ) {
      tp = SplineType1D::CONSTANT;
    } else if ( spl_type == "linear" ) {
      tp = SplineType1D::LINEAR;
    } else if ( spl_type == "cubic" ) {
      tp = SplineType1D::CUBIC;
    } else if ( spl_type == "akima" ) {
      tp = SplineType1D::AKIMA;
    } else if ( spl_type == "bessel" ) {
      tp = SplineType1D::BESSEL;
    } else if ( spl_type == "pchip" ) {
      tp = SplineType1D::PCHIP;
    } else if ( spl_type == "quintic" ) {
      tp = SplineType1D::QUINTIC;
    } else {
      UTILS_ERROR(
       "Spline1D::setup[{}] unknown type {}, not in "
       "[constant,linear,cubic,akima,bessel,pchip,quintic]\n",
       m_name, spl_type
      );
    }
    if ( m_pSpline != nullptr ) delete m_pSpline;
    m_pSpline = new_Spline1D( m_name, tp );
    m_pSpline->build( gc );
  }

}
