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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

namespace Splines {

  static
  Spline *
  new_Spline1D( string const & _name, SplineType1D tp ) {
    switch ( tp ) {
    case CONSTANT_TYPE:   return new ConstantSpline(_name);
    case LINEAR_TYPE:     return new LinearSpline(_name);
    case CUBIC_TYPE:      return new CubicSpline(_name);
    case AKIMA_TYPE:      return new AkimaSpline(_name);
    case BESSEL_TYPE:     return new BesselSpline(_name);
    case PCHIP_TYPE:      return new PchipSpline(_name);
    case QUINTIC_TYPE:    return new QuinticSpline(_name);
    case HERMITE_TYPE:    break;
    case SPLINE_SET_TYPE: break;
    case SPLINE_VEC_TYPE: break;
    }
    return nullptr;
  }

  void
  Spline1D::build(
    SplineType1D tp,
    real_type const * x, integer incx,
    real_type const * y, integer incy,
    integer n
  ) {
    if ( m_pSpline != nullptr ) delete m_pSpline;
    m_pSpline = new_Spline1D(m_name,tp);
    UTILS_ASSERT0( m_pSpline != nullptr, "Spline1D::build, failed\n" );
    m_pSpline->build( x, incx, y, incy, n );
  }

  /*
  //    ____  ____   ____                               _
  //   / ___|/ ___| / ___| _   _ _ __  _ __   ___  _ __| |_
  //  | |  _| |     \___ \| | | | '_ \| '_ \ / _ \| '__| __|
  //  | |_| | |___   ___) | |_| | |_) | |_) | (_) | |  | |_
  //   \____|\____| |____/ \__,_| .__/| .__/ \___/|_|   \__|
  //                            |_|   |_|
  */

  using GC_namespace::GC_VEC_REAL;
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
    std::string spl_type = gc("spline_type").get_string(
      "Spline1D::setup, spline_type expected to be a string\n"
    );
    SplineType1D tp;
    if ( spl_type == "constant" ) {
      tp = CONSTANT_TYPE;
    } else if ( spl_type == "linear" ) {
      tp = LINEAR_TYPE;
    } else if ( spl_type == "cubic" ) {
      tp = CUBIC_TYPE;
    } else if ( spl_type == "akima" ) {
      tp = AKIMA_TYPE;
    } else if ( spl_type == "bessel" ) {
      tp = BESSEL_TYPE;
    } else if ( spl_type == "pchip" ) {
      tp = PCHIP_TYPE;
    } else if ( spl_type == "quintic" ) {
      tp = QUINTIC_TYPE;
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
