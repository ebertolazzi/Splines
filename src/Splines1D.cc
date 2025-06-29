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

namespace Splines {

  static
  std::unique_ptr<Spline>
  new_Spline1D( string_view const _name, SplineType1D const tp ) {
    switch ( tp ) {
    case SplineType1D::CONSTANT:   return std::make_unique<ConstantSpline>(_name);
    case SplineType1D::LINEAR:     return std::make_unique<LinearSpline>(_name);
    case SplineType1D::CUBIC:      return std::make_unique<CubicSpline>(_name);
    case SplineType1D::AKIMA:      return std::make_unique<AkimaSpline>(_name);
    case SplineType1D::BESSEL:     return std::make_unique<BesselSpline>(_name);
    case SplineType1D::PCHIP:      return std::make_unique<PchipSpline>(_name);
    case SplineType1D::QUINTIC:    return std::make_unique<QuinticSpline>(_name);
    case SplineType1D::HERMITE:    break;
    case SplineType1D::SPLINE_SET: break;
    case SplineType1D::SPLINE_VEC: break;
    }
    return nullptr;
  }

  void
  Spline1D::build(
    SplineType1D const tp,
    real_type const x[], integer const incx,
    real_type const y[], integer const incy,
    integer const n
  ) {
    m_spline = new_Spline1D(m_name,tp);
    m_spline->build( x, incx, y, incy, n );
  }

  #ifdef AUTODIFF_SUPPORT
  autodiff::dual1st
  Spline1D::eval( autodiff::dual1st const & x ) const {
    using autodiff::dual1st;
    using autodiff::detail::val;
    real_type dd[2];
    D( val(x), dd );
    dual1st res { dd[0] };
    res.grad = dd[1] * x.grad;
    return res;
  }

  autodiff::dual2nd
  Spline1D::eval( autodiff::dual2nd const & x ) const {
    using autodiff::dual2nd;
    using autodiff::detail::val;
    real_type dd[3], xg{ val(x.grad) };
    DD( val(x), dd );
    dual2nd res { dd[0] };
    res.grad      = dd[1] * xg;
    res.grad.grad = dd[1] * x.grad.grad + dd[2] * (xg*xg);
    return res;
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
    string const where{ fmt::format("Spline1D[{}]::setup( gc ):", m_name ) };

    string_view spl_type{ gc.get_map_string("spline_type",where) };

    SplineType1D tp;
    if      ( spl_type == "constant" ) tp = SplineType1D::CONSTANT;
    else if ( spl_type == "linear"   ) tp = SplineType1D::LINEAR;
    else if ( spl_type == "cubic"    ) tp = SplineType1D::CUBIC;
    else if ( spl_type == "akima"    ) tp = SplineType1D::AKIMA;
    else if ( spl_type == "bessel"   ) tp = SplineType1D::BESSEL;
    else if ( spl_type == "pchip"    ) tp = SplineType1D::PCHIP;
    else if ( spl_type == "quintic"  ) tp = SplineType1D::QUINTIC;
    else {
      UTILS_ERROR(
       "Spline1D::setup[{}] unknown type {}, not in "
       "[constant,linear,cubic,akima,bessel,pchip,quintic]\n",
       m_name, spl_type
      );
    }
    m_spline = new_Spline1D( m_name, tp );
    m_spline->build( gc );
  }
}
