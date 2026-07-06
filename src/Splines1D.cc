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

#define AUTODIFF_SUPPORT
#include "Splines/Splines.hh"

namespace Splines
{

  /*\
   |   ____        _ _            _ ____
   |  / ___| _ __ | (_)_ __   ___/ |  _ \
   |  \___ \| '_ \| | | '_ \ / _ \ | | | |
   |   ___) | |_) | | | | | |  __/ | |_| |
   |  |____/| .__/|_|_|_| |_|\___|_|____/
   |        |_|
  \*/

  std::unique_ptr<Spline> Spline1D::new_Spline1D( string_view const _name, SplineType1D const tp )
  {
    switch ( tp )
    {
      case SplineType1D::CONSTANT: return std::make_unique<ConstantSpline>( _name );
      case SplineType1D::LINEAR: return std::make_unique<LinearSpline>( _name );
      case SplineType1D::CUBIC: return std::make_unique<CubicSpline>( _name );
      case SplineType1D::AKIMA: return std::make_unique<AkimaSpline>( _name );
      case SplineType1D::VANLEER: return std::make_unique<VanLeerSpline>( _name );
      case SplineType1D::PCHIP: return std::make_unique<PchipSpline>( _name );
      case SplineType1D::QUINTIC_CUBIC: return std::make_unique<QuinticSpline>( Spline_sub_type::CUBIC, _name );
      case SplineType1D::QUINTIC_AKIMA: return std::make_unique<QuinticSpline>( Spline_sub_type::AKIMA, _name );
      case SplineType1D::QUINTIC_VANLEER: return std::make_unique<QuinticSpline>( Spline_sub_type::VANLEER, _name );
      case SplineType1D::QUINTIC_PCHIP: return std::make_unique<QuinticSpline>( Spline_sub_type::PCHIP, _name );
      case SplineType1D::HERMITE: break;
      case SplineType1D::SPLINE_SET: break;
      case SplineType1D::SPLINE_VEC: break;
    }
    return nullptr;
  }

  void Spline1D::setup( GenericContainer const & gc )
  {
    /*
    // gc["xdata"]
    // gc["ydata"]
    //
    */
    string const where = fmt::format( "Spline1D[{}]::setup( gc ):", m_name );

    string_view spl_type = gc.get_map_string( "spline_type", where );

    SplineType1D tp;
    if ( spl_type == "constant" )
      tp = SplineType1D::CONSTANT;
    else if ( spl_type == "linear" )
      tp = SplineType1D::LINEAR;
    else if ( spl_type == "cubic" )
      tp = SplineType1D::CUBIC;
    else if ( spl_type == "akima" )
      tp = SplineType1D::AKIMA;
    else if ( spl_type == "vanleer" )
      tp = SplineType1D::VANLEER;
    else if ( spl_type == "pchip" )
      tp = SplineType1D::PCHIP;
    else if ( spl_type == "quintic_cubic" )
      tp = SplineType1D::QUINTIC_CUBIC;
    else if ( spl_type == "quintic_akima" )
      tp = SplineType1D::QUINTIC_AKIMA;
    else if ( spl_type == "quintic_vanleer" )
      tp = SplineType1D::QUINTIC_VANLEER;
    else if ( spl_type == "quintic_pchip" || spl_type == "quintic" )
      tp = SplineType1D::QUINTIC_PCHIP;
    else
    {
      UTILS_ERROR(
        "Spline1D::setup[{}] unknown type {}, not in "
        "[constant,linear,cubic,akima,vanleer,pchip,quintic]\n",
        m_name,
        spl_type );
    }
    m_spline = new_Spline1D( m_name, tp );
    m_spline->build( gc );
  }

  autodiff::dual1st Spline1D::eval( autodiff::dual1st const & x ) const
  {
    real_type dd[2];
    D( x.val, dd );
    autodiff::dual1st res;
    res.val  = dd[0];
    res.grad = dd[1] * x.grad;
    return res;
  }

  autodiff::dual2nd Spline1D::eval( autodiff::dual2nd const & x ) const
  {
    real_type dd[3];
    DD( x.val.val, dd );
    autodiff::dual2nd res;

    // Impostazione del valore
    res.val.val = dd[0];  // f(x)

    // Calcolo della derivata prima
    real_type xg    = x.grad.val;  // x' = derivata prima di x
    real_type dF_dX = dd[1] * xg;  // f'(x) * x'

    // Impostazione della derivata prima
    res.val.grad = dF_dX;
    res.grad.val = dF_dX;

    // Calcolo della derivata seconda
    // d²F/dX² = f''(x) * (x')² + f'(x) * x''
    real_type d2F_dX2 = dd[2] * ( xg * xg ) + dd[1] * x.grad.grad;
    res.grad.grad     = d2F_dX2;

    return res;
  }

}  // namespace Splines

//
// EOF: Spline1D.cc
//
